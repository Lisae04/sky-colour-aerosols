"""
Folder of sunset/sunrise photos -> per-frame table -> retrieval summary.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from .fit import DEFAULT_AIRMASS_RANGE, retrieve_from_frames
from .image import SUPPORTED_EXTENSIONS, find_sun, load_linear_rgb, sun_photometry
from .metadata import exif_time_to_timestamp, load_manifest, read_exif
from .physics import STANDARD_PRESSURE_HPA, relative_air_mass
from .solar import solar_position

FRAME_COLUMNS = [
    "file", "time", "lat", "lon", "elevation_true_deg", "azimuth_deg", "airmass",
    "found", "cx", "cy", "radius_px", "contrast", "saturated_frac",
    "flux_R", "flux_G", "flux_B", "sky_R", "sky_G", "sky_B",
    "snr_R", "snr_G", "snr_B", "black_frac_max",
    "CI_RB", "CI_GB", "sky_CI_RB", "source", "qc_ok", "qc_reason",
]


def list_photos(folder) -> list:
    folder = Path(folder)
    if not folder.is_dir():
        raise NotADirectoryError(f"Photo folder not found: {folder}")
    files = sorted(p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in SUPPORTED_EXTENSIONS)
    if not files:
        raise FileNotFoundError(f"No supported images in {folder} (extensions: {sorted(SUPPORTED_EXTENSIONS)})")
    return files


def process_folder(
    folder,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    tz: Optional[str] = None,
    manifest: Optional[str] = None,
    clock_offset_s: float = 0.0,
    roi: Optional[Tuple[int, int, int, int]] = None,
    max_saturated_frac: float = 0.002,
    min_elevation_deg: float = 0.5,
    min_snr: float = 30.0,
    max_black_frac: float = 0.10,
    verbose: bool = True,
) -> pd.DataFrame:
    """Measure every photo in `folder`. Never raises on a single bad frame;
    problems are recorded in `qc_reason` so the table is a complete audit trail.
    """
    files = list_photos(folder)
    man = load_manifest(manifest, tz) if manifest else None
    rows = []
    for i, path in enumerate(files, 1):
        row = {c: np.nan for c in FRAME_COLUMNS}
        row.update(file=path.name, found=False, qc_ok=False, qc_reason="", time=pd.NaT)
        reasons = []
        try:
            # --- time and place -------------------------------------------------
            f_lat, f_lon, t = lat, lon, pd.NaT
            if man is not None and path.name in man.index:
                rec = man.loc[path.name]
                t = rec["time"]
                f_lat = rec.get("lat", f_lat) if pd.notna(rec.get("lat", np.nan)) else f_lat
                f_lon = rec.get("lon", f_lon) if pd.notna(rec.get("lon", np.nan)) else f_lon
            else:
                ex = read_exif(path)
                t = exif_time_to_timestamp(ex["datetime"], ex["offset"], tz)
                if f_lat is None and np.isfinite(ex["lat"]):
                    f_lat = ex["lat"]
                if f_lon is None and np.isfinite(ex["lon"]):
                    f_lon = ex["lon"]
            if pd.isna(t):
                reasons.append("no timestamp (add OffsetTimeOriginal, pass --tz, or use a manifest)")
            if f_lat is None or f_lon is None:
                reasons.append("no location (pass --lat/--lon or use GPS-tagged photos)")
            if not reasons:
                t = t + pd.Timedelta(seconds=float(clock_offset_s))
                sp = solar_position(pd.DatetimeIndex([t]), float(f_lat), float(f_lon)).iloc[0]
                row.update(time=t, lat=float(f_lat), lon=float(f_lon),
                           elevation_true_deg=float(sp["elevation_true_deg"]),
                           azimuth_deg=float(sp["azimuth_deg"]))
                row["airmass"] = float(relative_air_mass(row["elevation_true_deg"]))
                if row["elevation_true_deg"] < min_elevation_deg:
                    reasons.append(f"sun elevation {row['elevation_true_deg']:.2f} deg < {min_elevation_deg}")

            # --- photometry -----------------------------------------------------
            rgb, sat, info = load_linear_rgb(path)
            row["source"] = info["source"]
            sun = find_sun(rgb, roi=roi)
            row.update(found=sun["found"], cx=sun["cx"], cy=sun["cy"],
                       radius_px=sun["radius_px"], contrast=sun["contrast"])
            if not sun["found"]:
                reasons.append(f"sun not found: {sun['reason']}")
            else:
                ph = sun_photometry(rgb, sat, sun["cx"], sun["cy"], sun["radius_px"])
                for k in ("flux_R", "flux_G", "flux_B", "sky_R", "sky_G", "sky_B", "saturated_frac",
                          "snr_R", "snr_G", "snr_B"):
                    row[k] = ph[k]
                row["black_frac_max"] = max(ph[f"black_frac_{c}"] for c in "RGB")
                weakest = min("RGB", key=lambda c: ph[f"snr_{c}"])
                if not ph[f"snr_{weakest}"] >= min_snr:
                    reasons.append(f"low signal-to-noise in {weakest} (SNR {ph[f'snr_{weakest}']:.0f} < {min_snr:g})")
                if row["black_frac_max"] > max_black_frac:
                    reasons.append(f"background black-clipped ({100 * row['black_frac_max']:.0f}% of sky ring at 0); "
                                   "channel underexposed")
                if ph["aperture_clipped_by_edge"]:
                    reasons.append("sun too close to image edge")
                if ph["saturated_frac"] > max_saturated_frac:
                    reasons.append(f"sun disk saturated ({100 * ph['saturated_frac']:.1f}% of aperture)")
                fluxes = np.array([ph["flux_R"], ph["flux_G"], ph["flux_B"]])
                if np.all(fluxes > 0):
                    row["CI_RB"] = float(np.log(ph["flux_R"] / ph["flux_B"]))
                    row["CI_GB"] = float(np.log(ph["flux_G"] / ph["flux_B"]))
                else:
                    reasons.append("non-positive background-subtracted flux")
                if ph["sky_R"] > 0 and ph["sky_B"] > 0:
                    row["sky_CI_RB"] = float(np.log(ph["sky_R"] / ph["sky_B"]))
        except Exception as exc:  # keep going; record the failure
            reasons.append(f"error: {type(exc).__name__}: {exc}")

        row["qc_ok"] = not reasons
        row["qc_reason"] = "; ".join(reasons)
        rows.append(row)
        if verbose:
            status = "ok" if row["qc_ok"] else f"skip ({row['qc_reason']})"
            print(f"[{i}/{len(files)}] {path.name}: {status}")

    df = pd.DataFrame(rows, columns=FRAME_COLUMNS)
    return df.sort_values("time", na_position="last").reset_index(drop=True)


def run(
    folder,
    out_dir,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    tz: Optional[str] = None,
    manifest: Optional[str] = None,
    clock_offset_s: float = 0.0,
    pressure_hpa: float = STANDARD_PRESSURE_HPA,
    airmass_range: Tuple[float, float] = DEFAULT_AIRMASS_RANGE,
    alpha_fixed: Optional[float] = None,
    wavelengths_um: Optional[Dict[str, float]] = None,
    roi=None,
    n_boot: int = 1000,
    slope_sys_rel: Optional[float] = None,
    wavelength_sd_um: float = 0.010,
    make_plot: bool = True,
    verbose: bool = True,
) -> Dict:
    """End-to-end: measure frames, fit, retrieve, write frames.csv / summary.json / plot."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frames = process_folder(folder, lat, lon, tz, manifest, clock_offset_s, roi=roi, verbose=verbose)
    frames.to_csv(out / "frames.csv", index=False)
    if slope_sys_rel is None:
        # Default systematic slope error by input type (see fit.retrieve_from_frames).
        jpeg_like = (frames["source"] == "srgb").any()
        slope_sys_rel = 0.02 if jpeg_like else 0.005
    summary = retrieve_from_frames(frames, wavelengths_um, pressure_hpa, airmass_range, alpha_fixed,
                                   n_boot, slope_sys_rel=slope_sys_rel,
                                   wavelength_sd_um=wavelength_sd_um)
    summary["input_folder"] = str(folder)
    summary["clock_offset_s"] = float(clock_offset_s)
    with open(out / "summary.json", "w") as fh:
        json.dump(_jsonable(summary), fh, indent=2)
    if make_plot:
        try:
            from .plots import plot_colour_index
            plot_colour_index(frames, summary, out / "colour_index_vs_airmass.png")
        except Exception as exc:
            warnings.warn(f"Plot failed: {exc}")
    return summary


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj
