import numpy as np
import pytest

from skycolour.physics import (aod_from_aerosol_slope, colour_index_slope, forward_signals,
                               pressure_from_altitude, rayleigh_optical_depth, relative_air_mass,
                               retrieve_beta_alpha)


def test_rayleigh_reference_values():
    # Cross-check against the independent Hansen & Travis (1974) formula,
    # tau = 0.008569 l^-4 (1 + 0.0113 l^-2 + 0.00013 l^-4), within 1 %.
    for lam in (0.44, 0.465, 0.535, 0.55, 0.61, 0.675):
        ht = 0.008569 * lam ** -4 * (1 + 0.0113 * lam ** -2 + 0.00013 * lam ** -4)
        assert rayleigh_optical_depth(lam) == pytest.approx(ht, rel=0.01)
    assert rayleigh_optical_depth(0.55) == pytest.approx(0.0973, abs=5e-4)
    # Scales linearly with pressure.
    assert rayleigh_optical_depth(0.55, 506.625) == pytest.approx(rayleigh_optical_depth(0.55) / 2)


def test_rayleigh_rejects_bad_input():
    with pytest.raises(ValueError):
        rayleigh_optical_depth(-0.5)
    with pytest.raises(ValueError):
        rayleigh_optical_depth(0.5, 0)


def test_air_mass_limits():
    assert relative_air_mass(90) == pytest.approx(1.0, abs=1e-3)
    assert relative_air_mass(30) == pytest.approx(1.994, abs=5e-3)
    assert relative_air_mass(0) == pytest.approx(37.9, abs=0.2)
    assert np.isnan(relative_air_mass(-5))


def test_pressure_from_altitude():
    assert pressure_from_altitude(0) == pytest.approx(1013.25)
    assert pressure_from_altitude(1800) == pytest.approx(815, abs=5)


def test_forward_model_is_linear_in_airmass():
    m = np.linspace(1, 20, 30)
    s = forward_signals(m, beta=0.12, alpha=1.4)
    ci = np.log(s["R"] / s["B"])
    k, _ = np.polyfit(m, ci, 1)
    assert k == pytest.approx(colour_index_slope(0.610, 0.465, 0.12, 1.4), rel=1e-10)


@pytest.mark.parametrize("beta,alpha", [(0.05, 0.8), (0.1, 1.3), (0.3, 1.8), (0.2, 2.5)])
def test_exact_inversion_round_trip(beta, alpha):
    k_rb = colour_index_slope(0.610, 0.465, beta, alpha)
    k_gb = colour_index_slope(0.535, 0.465, beta, alpha)
    r = retrieve_beta_alpha(k_rb, k_gb)
    assert r.status == "ok"
    assert r.alpha == pytest.approx(alpha, abs=1e-6)
    assert r.beta == pytest.approx(beta, rel=1e-6)


def test_clean_sky_is_flagged_not_inverted():
    k_rb = colour_index_slope(0.610, 0.465, 0.0)
    k_gb = colour_index_slope(0.535, 0.465, 0.0)
    r = retrieve_beta_alpha(k_rb, k_gb)
    assert r.status != "ok" and np.isnan(r.alpha)


def test_alpha_fixed_branch():
    k_rb = colour_index_slope(0.610, 0.465, 0.1, 1.3)
    r = retrieve_beta_alpha(k_rb, 0.0 + colour_index_slope(0.535, 0.465), alpha_fixed=1.3)
    assert r.beta == pytest.approx(0.1, rel=1e-6)


def test_grey_aerosol_is_invisible_in_colour():
    # alpha -> 0: large AOD, almost no colour signal.
    slope_grey = colour_index_slope(0.610, 0.465, 0.5, 0.0) - colour_index_slope(0.610, 0.465, 0.0)
    assert abs(slope_grey) < 1e-9
    # And the implied AOD for a given colour slope grows as alpha falls.
    assert aod_from_aerosol_slope(0.05, 0.25) > aod_from_aerosol_slope(0.05, 1.5)
