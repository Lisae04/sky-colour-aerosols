"""
Photo -> linear RGB -> sun-disk photometry.

The retrieval only needs RATIOS between channels within one fixed-white-balance
sequence, so absolute calibration is not required. What IS required:

* linear pixel values (RAW is best; JPEG is linearised with the sRGB curve,
  which is only approximate because cameras apply their own tone curves);
* an unsaturated sun disk (use a fast shutter, low ISO, or an ND filter);
* the same white balance for every frame (RAW, or WB locked to "Daylight").
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from PIL import Image
from scipy import ndimage

RAW_EXTENSIONS = {".cr2", ".cr3", ".nef", ".arw", ".dng", ".raf", ".orf", ".rw2", ".pef", ".srw"}
RASTER_EXTENSIONS = {".jpg", ".jpeg", ".png", ".tif", ".tiff"}
SUPPORTED_EXTENSIONS = RAW_EXTENSIONS | RASTER_EXTENSIONS

SATURATION_LEVEL = 0.98  # fraction of full scale treated as clipped


def srgb_to_linear(v: np.ndarray) -> np.ndarray:
    """Inverse sRGB transfer function (IEC 61966-2-1); input and output in [0, 1]."""
    v = np.clip(np.asarray(v, dtype=np.float64), 0.0, 1.0)
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)


def load_linear_rgb(path) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """Load an image as linear RGB in [0, 1] plus a per-pixel saturation mask.

    Returns (rgb[H, W, 3] float64, saturated[H, W] bool, info dict).
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    ext = path.suffix.lower()
    if ext not in SUPPORTED_EXTENSIONS:
        raise ValueError(f"Unsupported file type '{ext}'. Supported: {sorted(SUPPORTED_EXTENSIONS)}")

    if ext in RAW_EXTENSIONS:
        try:
            import rawpy  # optional dependency
        except ImportError as exc:
            raise ImportError("Reading RAW files needs 'rawpy' (pip install rawpy).") from exc
        with rawpy.imread(str(path)) as raw:
            # Linear, no auto-brightening, fixed unit white balance, sensor-native
            # colour space: every frame is processed identically.
            rgb16 = raw.postprocess(
                gamma=(1, 1), no_auto_bright=True, output_bps=16,
                use_camera_wb=False, use_auto_wb=False, user_wb=[1.0, 1.0, 1.0, 1.0],
                output_color=rawpy.ColorSpace.raw, half_size=True,
            )
        rgb = rgb16.astype(np.float64) / 65535.0
        saturated = np.any(rgb >= SATURATION_LEVEL, axis=-1)
        return rgb, saturated, {"source": "raw", "linearisation": "native"}

    with Image.open(path) as im:
        im.load()
        mode = im.mode
        if mode in ("I;16", "I;16B", "I"):
            raise ValueError(f"{path.name}: single-channel 16-bit image, expected RGB.")
        arr = np.asarray(im.convert("RGB") if mode not in ("RGB",) else im)
    if arr.ndim != 3 or arr.shape[-1] < 3:
        raise ValueError(f"{path.name}: expected an RGB image, got shape {arr.shape}.")
    arr = arr[..., :3]
    full_scale = 65535.0 if arr.dtype == np.uint16 else 255.0
    encoded = arr.astype(np.float64) / full_scale
    saturated = np.any(encoded >= SATURATION_LEVEL, axis=-1)
    return srgb_to_linear(encoded), saturated, {"source": "srgb", "linearisation": "sRGB curve (approx.)"}


def find_sun(
    rgb: np.ndarray,
    roi: Optional[Tuple[int, int, int, int]] = None,
    min_contrast: float = 5.0,
) -> Dict:
    """Locate the sun disk as the brightest connected blob.

    roi : optional (x0, y0, x1, y1) pixel box to search in, to avoid bright clouds
          or reflections being mistaken for the sun.
    Returns a dict with keys found, cx, cy, radius_px, contrast, reason.
    """
    if rgb.ndim != 3 or rgb.shape[-1] != 3:
        raise ValueError("rgb must be an H x W x 3 array.")
    lum = 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]
    search = np.zeros_like(lum, dtype=bool)
    if roi is None:
        search[:] = True
    else:
        x0, y0, x1, y1 = (int(v) for v in roi)
        if not (0 <= x0 < x1 <= lum.shape[1] and 0 <= y0 < y1 <= lum.shape[0]):
            raise ValueError(f"roi {roi} outside image of size {lum.shape[1]}x{lum.shape[0]}.")
        search[y0:y1, x0:x1] = True

    smooth = ndimage.uniform_filter(lum, size=5)
    smooth_search = np.where(search, smooth, -np.inf)
    py, px = np.unravel_index(np.argmax(smooth_search), smooth.shape)
    peak = smooth[py, px]
    background = float(np.median(lum[search]))
    contrast = peak / max(background, 1e-9)
    if not np.isfinite(peak) or contrast < min_contrast:
        return {"found": False, "cx": np.nan, "cy": np.nan, "radius_px": np.nan,
                "contrast": float(contrast), "reason": f"peak/background contrast {contrast:.1f} < {min_contrast}"}

    threshold = background + 0.5 * (peak - background)
    labels, _ = ndimage.label((smooth >= threshold) & search)
    blob = labels == labels[py, px]
    area = int(blob.sum())
    weights = np.where(blob, lum - background, 0.0)
    cy, cx = ndimage.center_of_mass(np.clip(weights, 0, None))
    return {"found": True, "cx": float(cx), "cy": float(cy), "radius_px": float(np.sqrt(area / np.pi)),
            "contrast": float(contrast), "reason": ""}


def sun_photometry(rgb: np.ndarray, saturated: np.ndarray, cx: float, cy: float, radius_px: float) -> Dict:
    """Background-subtracted integrated sun flux per channel, plus sky-ring colour.

    Aperture radius = max(1.5 r, r + 3) px; background = per-channel median of
    the annulus between 2 and 3 aperture radii (this is also reported as the
    near-sun sky colour).
    """
    if not (np.isfinite(cx) and np.isfinite(cy) and np.isfinite(radius_px) and radius_px > 0):
        raise ValueError("cx, cy and radius_px must be finite and radius_px > 0.")
    h, w, _ = rgb.shape
    yy, xx = np.ogrid[:h, :w]
    dist = np.hypot(xx - cx, yy - cy)
    r_ap = max(1.5 * radius_px, radius_px + 3.0)
    aperture = dist <= r_ap
    ring = (dist >= 2 * r_ap) & (dist <= 3 * r_ap)

    out = {"aperture_px": float(r_ap), "n_aperture": int(aperture.sum()), "n_ring": int(ring.sum()),
           "aperture_clipped_by_edge": bool(cx - 3 * r_ap < 0 or cy - 3 * r_ap < 0
                                             or cx + 3 * r_ap >= w or cy + 3 * r_ap >= h),
           "saturated_frac": float(saturated[aperture].mean()) if aperture.any() else np.nan}
    if out["n_ring"] < 20:
        raise ValueError("Background annulus too small; sun is too close to the image edge.")
    n_ap, n_ring = out["n_aperture"], out["n_ring"]
    for i, ch in enumerate("RGB"):
        ring_vals = rgb[..., i][ring]
        bg = float(np.median(ring_vals))
        # Robust per-pixel noise from the ring, floored (~ one dark 8-bit sRGB
        # step) so a perfectly flat or clipped ring cannot fake a huge SNR.
        sigma = max(1.4826 * float(np.median(np.abs(ring_vals - bg))), 1e-4)
        flux = float(rgb[..., i][aperture].sum() - bg * n_ap)
        flux_err = sigma * np.sqrt(n_ap + n_ap ** 2 / n_ring)  # pixel noise + background-level error
        out[f"sky_{ch}"] = bg
        out[f"flux_{ch}"] = flux
        out[f"snr_{ch}"] = flux / flux_err if flux_err > 0 else np.nan
        out[f"black_frac_{ch}"] = float((ring_vals <= 1e-6).mean())
    return out
