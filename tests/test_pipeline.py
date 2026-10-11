import json

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from skycolour.cli import main
from skycolour.fit import fit_colour_index, retrieve_from_frames
from skycolour.image import find_sun, load_linear_rgb, srgb_to_linear
from skycolour.metadata import exif_time_to_timestamp, load_manifest
from skycolour.physics import colour_index_slope
from skycolour.pipeline import process_folder, run
from skycolour.synthetic import make_sunset_sequence

TRUE_BETA, TRUE_ALPHA = 0.10, 1.3
TRUE_K_RB = colour_index_slope(0.610, 0.465, TRUE_BETA, TRUE_ALPHA)
TRUE_K_GB = colour_index_slope(0.535, 0.465, TRUE_BETA, TRUE_ALPHA)


@pytest.fixture(scope="module")
def png_sequence(tmp_path_factory):
    d = tmp_path_factory.mktemp("png")
    make_sunset_sequence(d, TRUE_BETA, TRUE_ALPHA, fmt="png", noise_rel=0.0, interval_s=90)
    return d


@pytest.fixture(scope="module")
def jpg_sequence(tmp_path_factory):
    d = tmp_path_factory.mktemp("jpg")
    make_sunset_sequence(d, TRUE_BETA, TRUE_ALPHA, fmt="jpg", noise_rel=0.003, interval_s=90, seed=3)
    return d


def test_srgb_linearisation_endpoints():
    assert srgb_to_linear(0.0) == 0.0
    assert srgb_to_linear(1.0) == pytest.approx(1.0)
    assert srgb_to_linear(0.5) == pytest.approx(0.214, abs=1e-3)


def test_exif_time_parsing():
    t = exif_time_to_timestamp("2026:10:09 18:05:00", "-04:00", None)
    assert t.tz is not None and t.tz_convert("UTC").hour == 22
    t2 = exif_time_to_timestamp("2026:10:09 18:05:00", None, "Europe/Paris")
    assert str(t2.tz) == "Europe/Paris"
    assert pd.isna(exif_time_to_timestamp("2026:10:09 18:05:00", None, None))
    assert pd.isna(exif_time_to_timestamp("garbage", None, "UTC"))


def test_png_round_trip_recovers_slopes(png_sequence):
    frames = process_folder(png_sequence, verbose=False)
    assert frames["qc_ok"].all()
    s = retrieve_from_frames(frames, n_boot=100, slope_sys_rel=0.0, wavelength_sd_um=0.0)
    assert s["fit_RB"]["slope"] == pytest.approx(TRUE_K_RB, rel=0.003)
    assert s["fit_GB"]["slope"] == pytest.approx(TRUE_K_GB, rel=0.003)
    assert s["retrieval"]["alpha"] == pytest.approx(TRUE_ALPHA, abs=0.15)
    # AOD for the true alpha brackets: alpha=1 above truth, alpha=1.5 below.
    true_aod = TRUE_BETA * 0.55 ** -TRUE_ALPHA
    g = s["aod_550_given_alpha"]
    assert g["1.5"] < true_aod < g["1"]


def test_jpeg_bias_is_covered_by_default_uncertainty(jpg_sequence, tmp_path):
    s = run(jpg_sequence, tmp_path / "out", n_boot=300, make_plot=False, verbose=False)
    # JPEG biases the slope by ~2 %: the point estimate is off...
    assert abs(s["fit_RB"]["slope"] / TRUE_K_RB - 1) < 0.04
    # ...but the reported 95 % interval for the aerosol colour slope still contains the truth.
    lo, hi = s["bootstrap"]["ci95"]["aerosol_colour_slope_RB"]
    assert lo < TRUE_K_RB - s["clean_sky_slope_RB"] < hi
    # And the free alpha inversion must not claim to be well constrained from a JPEG.
    assert s["retrieval"]["reliability"] != "constrained"
    assert (tmp_path / "out" / "frames.csv").exists()
    json.loads((tmp_path / "out" / "summary.json").read_text())


def test_saturated_sun_is_rejected(tmp_path):
    img = np.full((200, 300, 3), 40, np.uint8)
    yy, xx = np.mgrid[:200, :300]
    img[np.hypot(xx - 150, yy - 100) < 10] = 255
    Image.fromarray(img).save(tmp_path / "sat.png")
    frames = process_folder(tmp_path, lat=45.5, lon=-73.6, verbose=False)
    assert not frames["qc_ok"].iloc[0]
    assert "saturated" in frames["qc_reason"].iloc[0]
    assert "no timestamp" in frames["qc_reason"].iloc[0]


def test_no_sun_in_frame(tmp_path):
    Image.fromarray(np.full((100, 100, 3), 120, np.uint8)).save(tmp_path / "flat.png")
    rgb, _, _ = load_linear_rgb(tmp_path / "flat.png")
    assert not find_sun(rgb)["found"]


def test_manifest_overrides_exif(png_sequence, tmp_path):
    files = sorted(p.name for p in png_sequence.glob("*.png"))[:3]
    man = tmp_path / "manifest.csv"
    pd.DataFrame({"file": files, "time": ["2026-10-09 17:00:00"] * 3,
                  "lat": [48.8566] * 3, "lon": [2.3522] * 3}).to_csv(man, index=False)
    m = load_manifest(man, tz="Europe/Paris")
    assert str(m["time"].dt.tz) == "Europe/Paris"
    with pytest.raises(ValueError):
        load_manifest(man)  # naive times without tz


def test_fit_needs_enough_points():
    r = fit_colour_index([2, 3], [1, 2])
    assert r["status"] != "ok"


def test_retrieve_requires_columns():
    with pytest.raises(KeyError):
        retrieve_from_frames(pd.DataFrame({"airmass": [1, 2]}))


def test_missing_folder_raises():
    with pytest.raises(NotADirectoryError):
        process_folder("/definitely/not/here")


def test_cli_end_to_end(png_sequence, tmp_path, capsys):
    rc = main(["process", str(png_sequence), "--out", str(tmp_path / "cli"), "--n-boot", "100", "--quiet"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Aerosol colour slope" in out
    assert (tmp_path / "cli" / "colour_index_vs_airmass.png").exists()
