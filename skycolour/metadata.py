"""
Photo timestamps and location.

Order of precedence for each frame:
1. a manifest CSV (columns: file, time[, lat, lon]) if given;
2. EXIF DateTimeOriginal (+ OffsetTimeOriginal, + GPS if present);
3. otherwise the frame is reported with time = NaT and skipped downstream.

Camera clocks drift. Near the horizon, 1 minute of clock error is ~0.15 deg of
solar elevation and changes air mass by >10 %, so correct it with
`clock_offset_s` (photograph a phone showing GPS time once per session).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd


def _ratio_to_float(v) -> float:
    if hasattr(v, "num") and hasattr(v, "den"):  # exifread.Ratio
        return float(v.num) / float(v.den) if v.den else np.nan
    if isinstance(v, tuple) and len(v) == 2:
        return float(v[0]) / float(v[1]) if v[1] else np.nan
    return float(v)


def _dms_to_deg(values, ref) -> float:
    vals = list(values)
    if len(vals) != 3:
        return np.nan
    deg = _ratio_to_float(vals[0]) + _ratio_to_float(vals[1]) / 60 + _ratio_to_float(vals[2]) / 3600
    return -deg if str(ref).strip().upper() in ("S", "W") else deg


def read_exif(path) -> Dict:
    """Return {'datetime': str|None, 'offset': str|None, 'lat': float, 'lon': float}."""
    path = Path(path)
    out = {"datetime": None, "offset": None, "lat": np.nan, "lon": np.nan}
    try:
        import exifread
        with open(path, "rb") as fh:
            tags = exifread.process_file(fh, details=False)
        dt = tags.get("EXIF DateTimeOriginal") or tags.get("Image DateTime")
        out["datetime"] = str(dt) if dt else None
        off = tags.get("EXIF OffsetTimeOriginal") or tags.get("EXIF OffsetTime")
        out["offset"] = str(off).strip() if off else None
        if "GPS GPSLatitude" in tags and "GPS GPSLongitude" in tags:
            out["lat"] = _dms_to_deg(tags["GPS GPSLatitude"].values, tags.get("GPS GPSLatitudeRef", "N"))
            out["lon"] = _dms_to_deg(tags["GPS GPSLongitude"].values, tags.get("GPS GPSLongitudeRef", "E"))
        return out
    except ImportError:
        pass
    except Exception:  # unreadable EXIF: fall back to PIL below
        pass
    try:
        from PIL import Image
        with Image.open(path) as im:
            exif = im.getexif()
        sub = exif.get_ifd(0x8769)
        out["datetime"] = sub.get(36867) or exif.get(306)
        out["offset"] = sub.get(36881)
        gps = exif.get_ifd(0x8825)
        if 2 in gps and 4 in gps:
            out["lat"] = _dms_to_deg(gps[2], gps.get(1, "N"))
            out["lon"] = _dms_to_deg(gps[4], gps.get(3, "E"))
    except Exception:
        pass
    return out


def exif_time_to_timestamp(dt_str: Optional[str], offset: Optional[str], tz: Optional[str]) -> pd.Timestamp:
    """Combine an EXIF 'YYYY:MM:DD HH:MM:SS' string with an offset or IANA zone."""
    if not dt_str:
        return pd.NaT
    try:
        naive = pd.to_datetime(str(dt_str).strip(), format="%Y:%m:%d %H:%M:%S")
    except (ValueError, TypeError):
        return pd.NaT
    if offset:
        try:
            return pd.Timestamp(f"{naive.isoformat()}{offset}")
        except ValueError:
            pass
    if tz:
        return naive.tz_localize(tz)
    return pd.NaT


def load_manifest(path, tz: Optional[str] = None) -> pd.DataFrame:
    """Read a manifest CSV with at least columns 'file' and 'time'."""
    df = pd.read_csv(path)
    missing = {"file", "time"} - set(df.columns)
    if missing:
        raise KeyError(f"Manifest {path} is missing columns: {sorted(missing)}")
    times = pd.to_datetime(df["time"], errors="coerce", utc=False)
    if getattr(times.dt, "tz", None) is None:
        if not tz:
            raise ValueError("Manifest times have no UTC offset; pass tz (e.g. 'America/Toronto').")
        times = times.dt.tz_localize(tz)
    df["time"] = times
    for col in ("lat", "lon"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["file"] = df["file"].astype(str)
    return df.set_index("file")
