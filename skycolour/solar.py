"""
Solar position from time and location (NOAA solar calculator equations,
after Meeus, 1991). Accuracy is about 0.01 deg in elevation for 1800-2100,
which is far better than camera-clock accuracy, the real limiting factor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _to_utc_index(times) -> pd.DatetimeIndex:
    idx = pd.DatetimeIndex(pd.to_datetime(times))
    if idx.tz is None:
        raise ValueError(
            "times must be timezone-aware (e.g. pd.Timestamp('2026-10-09 18:20', tz='America/Toronto'))."
        )
    return idx.tz_convert("UTC")


def solar_position(times, lat_deg: float, lon_deg: float) -> pd.DataFrame:
    """Return true elevation, apparent (refracted) elevation and azimuth in degrees.

    Parameters
    ----------
    times    : timezone-aware timestamps (scalar or array-like)
    lat_deg  : latitude, + north
    lon_deg  : longitude, + east
    """
    if not (-90 <= lat_deg <= 90):
        raise ValueError("lat_deg must be in [-90, 90].")
    if not (-180 <= lon_deg <= 180):
        raise ValueError("lon_deg must be in [-180, 180].")
    t = _to_utc_index(np.atleast_1d(times) if not isinstance(times, pd.DatetimeIndex) else times)

    jd = t.to_julian_date().to_numpy()
    jc = (jd - 2451545.0) / 36525.0

    l0 = np.mod(280.46646 + jc * (36000.76983 + jc * 0.0003032), 360.0)
    m = 357.52911 + jc * (35999.05029 - 0.0001537 * jc)
    e = 0.016708634 - jc * (0.000042037 + 0.0000001267 * jc)
    mr = np.radians(m)
    c = (np.sin(mr) * (1.914602 - jc * (0.004817 + 0.000014 * jc))
         + np.sin(2 * mr) * (0.019993 - 0.000101 * jc)
         + np.sin(3 * mr) * 0.000289)
    true_long = l0 + c
    omega = 125.04 - 1934.136 * jc
    app_long = true_long - 0.00569 - 0.00478 * np.sin(np.radians(omega))
    mean_obl = 23 + (26 + (21.448 - jc * (46.815 + jc * (0.00059 - jc * 0.001813))) / 60) / 60
    obl = mean_obl + 0.00256 * np.cos(np.radians(omega))
    decl = np.degrees(np.arcsin(np.sin(np.radians(obl)) * np.sin(np.radians(app_long))))

    y = np.tan(np.radians(obl / 2)) ** 2
    l0r = np.radians(l0)
    eq_time = 4 * np.degrees(
        y * np.sin(2 * l0r) - 2 * e * np.sin(mr) + 4 * e * y * np.sin(mr) * np.cos(2 * l0r)
        - 0.5 * y * y * np.sin(4 * l0r) - 1.25 * e * e * np.sin(2 * mr)
    )  # minutes

    minutes_utc = (t.hour * 60 + t.minute + t.second / 60 + t.microsecond / 6e7).to_numpy()
    true_solar_time = np.mod(minutes_utc + eq_time + 4 * lon_deg, 1440.0)
    hour_angle = true_solar_time / 4 - 180.0

    lat_r, decl_r, ha_r = np.radians(lat_deg), np.radians(decl), np.radians(hour_angle)
    cos_z = np.clip(np.sin(lat_r) * np.sin(decl_r) + np.cos(lat_r) * np.cos(decl_r) * np.cos(ha_r), -1, 1)
    zenith = np.degrees(np.arccos(cos_z))
    elevation = 90.0 - zenith

    # Azimuth, clockwise from north
    sin_z = np.sin(np.radians(zenith))
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_az = (np.sin(lat_r) * cos_z - np.sin(decl_r)) / (np.cos(lat_r) * sin_z)
    az = np.degrees(np.arccos(np.clip(cos_az, -1, 1)))
    azimuth = np.where(hour_angle > 0, np.mod(az + 180, 360), np.mod(540 - az, 360))

    return pd.DataFrame(
        {
            "elevation_true_deg": elevation,
            "elevation_apparent_deg": elevation + atmospheric_refraction_deg(elevation),
            "azimuth_deg": azimuth,
        },
        index=t,
    )


def atmospheric_refraction_deg(elevation_deg):
    """NOAA approximation of atmospheric refraction (deg) for standard conditions."""
    el = np.asarray(elevation_deg, dtype=float)
    te = np.tan(np.radians(el))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(
            el > 85, 0.0,
            np.where(
                el > 5, 58.1 / te - 0.07 / te ** 3 + 0.000086 / te ** 5,
                np.where(
                    el > -0.575,
                    1735 + el * (-518.2 + el * (103.4 + el * (-12.79 + el * 0.711))),
                    -20.772 / te,
                ),
            ),
        )
    return r / 3600.0
