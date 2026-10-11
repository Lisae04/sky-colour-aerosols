"""
Synthetic sunset sequences with known aerosol (beta, alpha), written as real
JPEGs with EXIF time/GPS. Used to test the whole pipeline end to end: if the
retrieval cannot recover the beta and alpha it was given here, it will not
work on real photos either.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage

from .physics import DEFAULT_WAVELENGTHS_UM, STANDARD_PRESSURE_HPA, forward_signals, relative_air_mass
from .solar import solar_position

# Relative camera response x exo-atmospheric sun colour, per channel. Arbitrary:
# it lands in the intercept and must NOT affect the retrieved slopes.
CAMERA_GAIN = {"R": 1.00, "G": 1.15, "B": 0.80}


def _linear_to_srgb(v):
    v = np.clip(v, 0, 1)
    return np.where(v <= 0.0031308, 12.92 * v, 1.055 * v ** (1 / 2.4) - 0.055)


def _deg_to_dms_rationals(value: float):
    value = abs(value)
    d = int(value)
    m = int((value - d) * 60)
    s = ((value - d) * 60 - m) * 60
    from PIL.TiffImagePlugin import IFDRational
    return (IFDRational(d, 1), IFDRational(m, 1), IFDRational(round(s * 1000), 1000))


def _exif(t: pd.Timestamp, lat: float, lon: float) -> Image.Exif:
    exif = Image.Exif()
    exif[306] = t.strftime("%Y:%m:%d %H:%M:%S")
    sub = exif.get_ifd(0x8769)
    sub[36867] = t.strftime("%Y:%m:%d %H:%M:%S")
    off = t.strftime("%z")
    sub[36881] = f"{off[:3]}:{off[3:]}"
    gps = exif.get_ifd(0x8825)
    gps[1] = "N" if lat >= 0 else "S"
    gps[2] = _deg_to_dms_rationals(lat)
    gps[3] = "E" if lon >= 0 else "W"
    gps[4] = _deg_to_dms_rationals(lon)
    return exif


def make_sunset_sequence(
    out_dir,
    beta: float = 0.10,
    alpha: float = 1.3,
    lat: float = 45.5017,
    lon: float = -73.5673,
    date: str = "2026-10-09",
    tz: str = "America/Toronto",
    start_elevation_deg: float = 15.0,
    end_elevation_deg: float = 0.5,
    interval_s: int = 60,
    pressure_hpa: float = STANDARD_PRESSURE_HPA,
    wavelengths_um: Optional[Dict[str, float]] = None,
    noise_rel: float = 0.003,
    size=(480, 320),
    sun_radius_px: float = 9.0,
    fmt: str = "jpg",
    seed: int = 0,
) -> pd.DataFrame:
    """Write one image per `interval_s` while the sun goes from start to end elevation.

    Returns the ground-truth table (time, elevation, air mass, true signals).
    """
    if beta < 0:
        raise ValueError("beta must be >= 0.")
    if fmt not in ("jpg", "png"):
        raise ValueError("fmt must be 'jpg' or 'png'.")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    # Find the evening times covering the requested elevation window.
    day = pd.date_range(f"{date} 12:00", f"{date} 23:59", freq=f"{interval_s}s", tz=tz)
    sp = solar_position(day, lat, lon)
    el = sp["elevation_true_deg"].to_numpy()
    sel = (el <= start_elevation_deg) & (el >= end_elevation_deg)
    if sel.sum() < 5:
        raise ValueError("Fewer than 5 frames fall in the elevation window; check date/location.")
    times = day[sel]
    el = el[sel]
    m = relative_air_mass(el)
    sig = forward_signals(m, beta, alpha, wavelengths_um, pressure_hpa, f0=CAMERA_GAIN)

    w, h = size
    yy, xx = np.mgrid[:h, :w]
    truth = []
    for i, t in enumerate(times):
        # Sun drifts down and slightly right across the frame as it sets.
        frac = i / max(len(times) - 1, 1)
        cx, cy = 0.35 * w + 0.3 * w * frac, 0.25 * h + 0.5 * h * frac
        disk = (np.hypot(xx - cx, yy - cy) <= sun_radius_px).astype(float)
        disk = ndimage.gaussian_filter(disk, 1.2)
        s = np.array([sig[c][i] for c in "RGB"])
        exposure = 0.75 / s.max()  # photographer keeps the disk below clipping
        img = np.zeros((h, w, 3))
        for k in range(3):
            sky = 0.02 * s[k] ** 0.5 * (1.2 - 0.6 * yy / h) / s.max() ** 0.5  # smooth glow
            img[..., k] = exposure * s[k] * disk + sky
        img += rng.normal(0, noise_rel, img.shape)
        enc = np.round(_linear_to_srgb(img) * 255).astype(np.uint8)
        name = f"frame_{i:03d}.{fmt}"
        pil = Image.fromarray(enc, "RGB")
        ex = _exif(t, lat, lon)
        if fmt == "jpg":
            pil.save(out / name, quality=95, exif=ex)
        else:
            pil.save(out / name, exif=ex)
        truth.append({"file": name, "time": t, "elevation_true_deg": el[i], "airmass": m[i],
                      **{f"true_{c}": sig[c][i] for c in "RGB"}})
    gt = pd.DataFrame(truth)
    gt.to_csv(out / "_ground_truth.csv", index=False)
    return gt
