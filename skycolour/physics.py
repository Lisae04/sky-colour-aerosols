"""
Stage A physics: clean-sky baseline, colour index vs air mass, and aerosol retrieval.

Model (direct sun, Beer-Lambert):

    S_c(m) = K_c * F0_c * exp[-m * (tau_R(lambda_c) + tau_a(lambda_c))]
    tau_a(lambda) = beta * lambda**(-alpha)          (Angstrom, lambda in um)

For two camera channels X (longer wavelength) and Y (shorter wavelength):

    CI_XY(m) = ln(S_X / S_Y) = c_XY + m * k_XY
    k_XY     = [tau_R(Y) - tau_R(X)] + beta * [Y**-alpha - X**-alpha]

c_XY absorbs the unknown camera response, white balance and exo-atmospheric
solar colour, so it never needs to be known (Langley-style retrieval). Only
the SLOPE k_XY carries atmospheric information.

With two independent slopes (R/B and G/B) we can solve for both beta and alpha.

Units
-----
wavelength: micrometres (um)
pressure:   hPa
angles:     degrees
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, Optional, Sequence

import numpy as np
from scipy.optimize import brentq

# Effective wavelengths (um) of typical Bayer RGB channels. These are
# approximations: broadband channels shift their effective wavelength as the
# spectrum reddens. Override with values for your sensor if you have them.
DEFAULT_WAVELENGTHS_UM: Dict[str, float] = {"R": 0.610, "G": 0.535, "B": 0.465}

STANDARD_PRESSURE_HPA = 1013.25
ALPHA_SEARCH_RANGE = (-1.0, 4.0)  # physically plausible Angstrom exponents


# ---------------------------------------------------------------------------
# Rayleigh and air mass
# ---------------------------------------------------------------------------

def rayleigh_optical_depth(wavelength_um, pressure_hpa: float = STANDARD_PRESSURE_HPA):
    """Rayleigh optical depth (Bodhaine et al., 1999, eq. 30), scaled by pressure.

    Valid roughly 0.25-1.0 um. Returns an array with the shape of `wavelength_um`.
    """
    lam = np.asarray(wavelength_um, dtype=float)
    if np.any(lam <= 0):
        raise ValueError("wavelength_um must be positive (micrometres).")
    if not np.isfinite(pressure_hpa) or pressure_hpa <= 0:
        raise ValueError("pressure_hpa must be a positive finite number.")
    l2 = lam ** 2
    tau = 0.0021520 * (1.0455996 - 341.29061 / l2 - 0.90230850 * l2) / (
        1.0 + 0.0027059889 / l2 - 85.968563 * l2
    )
    return tau * (pressure_hpa / STANDARD_PRESSURE_HPA)


def pressure_from_altitude(altitude_m: float) -> float:
    """Standard-atmosphere surface pressure (hPa) for a site altitude (m).

    Useful for mountain sites (e.g. ~1800 m gives ~815 hPa).
    """
    if not np.isfinite(altitude_m):
        raise ValueError("altitude_m must be finite.")
    return STANDARD_PRESSURE_HPA * (1.0 - 2.25577e-5 * altitude_m) ** 5.25588


def relative_air_mass(elevation_deg):
    """Relative optical air mass (Kasten & Young, 1989) from TRUE solar elevation.

    About 1 at zenith and ~38 at the horizon. Returns NaN for elevations below
    -1 deg, where the formula is not meaningful.
    """
    el = np.asarray(elevation_deg, dtype=float)
    z = 90.0 - el
    with np.errstate(invalid="ignore", divide="ignore"):
        m = 1.0 / (np.cos(np.radians(z)) + 0.50572 * (96.07995 - z) ** (-1.6364))
    return np.where(el >= -1.0, m, np.nan)


# ---------------------------------------------------------------------------
# Aerosol (Angstrom) and forward model
# ---------------------------------------------------------------------------

def aerosol_optical_depth(wavelength_um, beta: float, alpha: float):
    """Angstrom law: tau_a = beta * lambda**(-alpha), lambda in um."""
    lam = np.asarray(wavelength_um, dtype=float)
    return beta * lam ** (-alpha)


def _spectral_diff(alpha: float, lam_short: float, lam_long: float) -> float:
    """lam_short**-alpha - lam_long**-alpha, numerically safe near alpha = 0."""
    a, b = np.log(lam_short), np.log(lam_long)
    if abs(alpha) < 1e-8:
        return alpha * (b - a)  # first-order Taylor term
    return float(np.exp(-alpha * a) - np.exp(-alpha * b))


def colour_index_slope(
    lam_long_um: float,
    lam_short_um: float,
    beta: float = 0.0,
    alpha: float = 1.3,
    pressure_hpa: float = STANDARD_PRESSURE_HPA,
) -> float:
    """Theoretical slope k of CI = ln(S_long/S_short) against air mass."""
    if lam_long_um <= lam_short_um:
        raise ValueError("lam_long_um must be longer than lam_short_um.")
    k_ray = float(
        rayleigh_optical_depth(lam_short_um, pressure_hpa)
        - rayleigh_optical_depth(lam_long_um, pressure_hpa)
    )
    return k_ray + beta * _spectral_diff(alpha, lam_short_um, lam_long_um)


def forward_signals(
    airmass,
    beta: float,
    alpha: float,
    wavelengths_um: Optional[Dict[str, float]] = None,
    pressure_hpa: float = STANDARD_PRESSURE_HPA,
    f0: Optional[Dict[str, float]] = None,
) -> Dict[str, np.ndarray]:
    """Direct-sun signal per channel (arbitrary units) for given air masses."""
    wl = dict(DEFAULT_WAVELENGTHS_UM if wavelengths_um is None else wavelengths_um)
    f0 = f0 or {"R": 1.0, "G": 1.0, "B": 0.9}
    m = np.asarray(airmass, dtype=float)
    out = {}
    for ch, lam in wl.items():
        tau = rayleigh_optical_depth(lam, pressure_hpa) + aerosol_optical_depth(lam, beta, alpha)
        out[ch] = f0.get(ch, 1.0) * np.exp(-m * tau)
    return out


# ---------------------------------------------------------------------------
# Inversion: slopes -> (beta, alpha)
# ---------------------------------------------------------------------------

def aod_from_aerosol_slope(
    aerosol_slope: float,
    alpha: float,
    lam_long_um: float = DEFAULT_WAVELENGTHS_UM["R"],
    lam_short_um: float = DEFAULT_WAVELENGTHS_UM["B"],
    lam_out_um: float = 0.55,
) -> float:
    """AOD at `lam_out_um` implied by an aerosol colour slope, for an ASSUMED alpha.

    Colour only sees the *spectral contrast* of aerosol extinction,
    beta * (lam_short**-alpha - lam_long**-alpha). Converting that into a total
    AOD needs alpha. As alpha -> 0 (coarse dust, sea salt) the contrast vanishes
    and the implied AOD diverges: grey aerosol is nearly invisible in colour.
    """
    d = _spectral_diff(alpha, lam_short_um, lam_long_um)
    if abs(d) < 1e-9:
        return np.nan
    return float(aerosol_slope / d * lam_out_um ** (-alpha))

@dataclass
class AerosolRetrieval:
    beta: float                  # Angstrom turbidity (AOD at 1 um)
    alpha: float                 # Angstrom exponent
    aod_550: float               # AOD at 550 nm, the usual reporting wavelength
    aerosol_slope_rb: float      # aerosol part of k_RB
    aerosol_slope_gb: float      # aerosol part of k_GB
    rayleigh_slope_rb: float     # clean-sky (beta = 0) baseline slope k_RB
    rayleigh_slope_gb: float
    status: str                  # "ok" or a reason the retrieval is unreliable

    def to_dict(self):
        return asdict(self)


def retrieve_beta_alpha(
    slope_rb: float,
    slope_gb: float,
    wavelengths_um: Optional[Dict[str, float]] = None,
    pressure_hpa: float = STANDARD_PRESSURE_HPA,
    alpha_fixed: Optional[float] = None,
) -> AerosolRetrieval:
    """Invert measured colour-index slopes into Angstrom beta and alpha.

    If `alpha_fixed` is given, only the R/B slope is used and beta is solved
    for that alpha (useful when the G/B slope is noisy).
    """
    wl = dict(DEFAULT_WAVELENGTHS_UM if wavelengths_um is None else wavelengths_um)
    for key in ("R", "G", "B"):
        if key not in wl:
            raise KeyError(f"wavelengths_um must contain '{key}'.")
    lr, lg, lb = wl["R"], wl["G"], wl["B"]
    if not (lr > lg > lb):
        raise ValueError("Expected wavelengths R > G > B.")
    for name, v in (("slope_rb", slope_rb), ("slope_gb", slope_gb)):
        if v is None or not np.isfinite(v):
            raise ValueError(f"{name} must be a finite number.")

    k_ray_rb = colour_index_slope(lr, lb, 0.0, 0.0, pressure_hpa)
    k_ray_gb = colour_index_slope(lg, lb, 0.0, 0.0, pressure_hpa)
    a_rb = slope_rb - k_ray_rb
    a_gb = slope_gb - k_ray_gb

    def pack(beta, alpha, status):
        aod = float(beta * 0.55 ** (-alpha)) if np.isfinite(beta) and np.isfinite(alpha) else np.nan
        return AerosolRetrieval(float(beta), float(alpha), aod, float(a_rb), float(a_gb),
                                float(k_ray_rb), float(k_ray_gb), status)

    if alpha_fixed is not None:
        d = _spectral_diff(alpha_fixed, lb, lr)
        if abs(d) < 1e-9:
            return pack(np.nan, alpha_fixed, "alpha_fixed gives no spectral contrast")
        return pack(a_rb / d, alpha_fixed, "ok (alpha fixed)")

    # Spectrally flat or negligible aerosol: beta is not separable from the
    # colour slopes, because grey extinction does not change colour.
    if a_rb <= 1e-4:
        return pack(np.nan, np.nan,
                    "aerosol colour signal below detection (clean sky or spectrally flat aerosol)")
    if a_gb <= 0:
        return pack(np.nan, np.nan, "inconsistent slopes (G/B aerosol part <= 0)")

    target = a_rb / a_gb

    def ratio(alpha):
        return _spectral_diff(alpha, lb, lr) / _spectral_diff(alpha, lb, lg)

    lo, hi = ALPHA_SEARCH_RANGE
    f_lo, f_hi = ratio(lo) - target, ratio(hi) - target
    if f_lo * f_hi > 0:
        return pack(np.nan, np.nan,
                    f"slope ratio {target:.3f} outside the range reachable for "
                    f"alpha in [{lo}, {hi}] ({ratio(lo):.3f}-{ratio(hi):.3f})")
    alpha = brentq(lambda a: ratio(a) - target, lo, hi, xtol=1e-10)
    beta = a_rb / _spectral_diff(alpha, lb, lr)
    return pack(beta, alpha, "ok")
