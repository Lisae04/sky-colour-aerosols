import numpy as np
import pandas as pd
import pytest

from skycolour.solar import solar_position


def test_solstice_noon_montreal():
    # Max elevation at 45.5 N on the June solstice = 90 - 45.50 + 23.44 = 67.94 deg.
    t = pd.date_range("2026-06-21 12:30", "2026-06-21 13:45", freq="1min", tz="America/Toronto")
    el = solar_position(t, 45.5017, -73.5673)["elevation_true_deg"].max()
    assert el == pytest.approx(67.94, abs=0.05)


def test_requires_timezone():
    with pytest.raises(ValueError):
        solar_position(pd.DatetimeIndex(["2026-10-09 18:00"]), 45.5, -73.6)


def test_rejects_bad_coordinates():
    t = pd.DatetimeIndex(["2026-10-09 18:00"]).tz_localize("UTC")
    with pytest.raises(ValueError):
        solar_position(t, 95, 0)


def test_west_at_sunset_east_at_sunrise():
    t = pd.DatetimeIndex(["2026-10-09 06:50", "2026-10-09 18:10"]).tz_localize("Asia/Beirut")
    az = solar_position(t, 33.89, 35.50)["azimuth_deg"].to_numpy()
    assert 60 < az[0] < 120 and 240 < az[1] < 300


def test_matches_pvlib_spa():
    pvlib = pytest.importorskip("pvlib")
    t = pd.date_range("2026-01-01", "2026-12-31", freq="7h", tz="Europe/Paris")
    a = solar_position(t, 48.8566, 2.3522)
    b = pvlib.solarposition.get_solarposition(t, 48.8566, 2.3522, method="nrel_numpy")
    assert np.abs(a["elevation_true_deg"].to_numpy() - b["elevation"].to_numpy()).max() < 0.02
