"""
Colour-index vs air-mass regression (Langley-style) and bootstrap retrieval.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .physics import (DEFAULT_WAVELENGTHS_UM, STANDARD_PRESSURE_HPA, aod_from_aerosol_slope,
                      colour_index_slope, retrieve_beta_alpha)

DEFAULT_AIRMASS_RANGE = (1.5, 20.0)  # ~42 deg down to ~1.9 deg true elevation
# Assumed Angstrom exponents for the conditional AOD table:
# ~0.25 coarse dust/sea salt, ~1.0 mixed, ~1.5 urban/continental, ~2.0 fresh smoke.
ALPHA_GRID = (0.25, 0.5, 1.0, 1.5, 2.0)
ALPHA_CI_WIDTH_OK = 0.6  # free alpha retrieval is trusted only below this 95% CI width


def fit_colour_index(
    airmass: Sequence[float],
    ci: Sequence[float],
    airmass_range: Tuple[float, float] = DEFAULT_AIRMASS_RANGE,
    clip_sigma: float = 3.0,
    max_iter: int = 5,
    min_points: int = 5,
) -> Dict:
    """Ordinary least squares CI = c + k*m with iterative sigma clipping.

    Returns slope, intercept, slope standard error, R^2, number of points used
    and a boolean mask of the points kept (aligned with the inputs).
    """
    m = np.asarray(airmass, dtype=float)
    y = np.asarray(ci, dtype=float)
    if m.shape != y.shape:
        raise ValueError("airmass and ci must have the same length.")
    lo, hi = airmass_range
    keep = np.isfinite(m) & np.isfinite(y) & (m >= lo) & (m <= hi)

    result = {"slope": np.nan, "intercept": np.nan, "slope_se": np.nan, "r2": np.nan,
              "n": int(keep.sum()), "n_rejected": 0, "mask": keep.copy(),
              "status": "ok"}
    if keep.sum() < min_points:
        result["status"] = f"only {int(keep.sum())} usable frames in air mass {lo}-{hi} (need {min_points})"
        return result
    if np.ptp(m[keep]) < 1.0:
        result["status"] = "air-mass span < 1; shoot over a longer stretch of the sunset"
        return result

    for _ in range(max_iter):
        A = np.column_stack([np.ones(keep.sum()), m[keep]])
        coef, *_ = np.linalg.lstsq(A, y[keep], rcond=None)
        resid_all = y - (coef[0] + coef[1] * m)
        s = np.std(resid_all[keep], ddof=2) if keep.sum() > 2 else 0.0
        new_keep = keep & (np.abs(resid_all) <= clip_sigma * max(s, 1e-12))
        if new_keep.sum() < min_points or new_keep.sum() == keep.sum():
            break
        keep = new_keep

    A = np.column_stack([np.ones(keep.sum()), m[keep]])
    coef, *_ = np.linalg.lstsq(A, y[keep], rcond=None)
    resid = y[keep] - A @ coef
    n = int(keep.sum())
    dof = max(n - 2, 1)
    sigma2 = float(resid @ resid) / dof
    cov = sigma2 * np.linalg.inv(A.T @ A)
    ss_tot = float(((y[keep] - y[keep].mean()) ** 2).sum())
    result.update(
        slope=float(coef[1]), intercept=float(coef[0]), slope_se=float(np.sqrt(cov[1, 1])),
        r2=float(1 - (resid @ resid) / ss_tot) if ss_tot > 0 else np.nan,
        n=n, n_rejected=int(result["n"] - n), mask=keep,
    )
    return result


def retrieve_from_frames(
    frames: pd.DataFrame,
    wavelengths_um: Optional[Dict[str, float]] = None,
    pressure_hpa: float = STANDARD_PRESSURE_HPA,
    airmass_range: Tuple[float, float] = DEFAULT_AIRMASS_RANGE,
    alpha_fixed: Optional[float] = None,
    n_boot: int = 1000,
    seed: int = 42,
    slope_sys_rel: float = 0.01,
    wavelength_sd_um: float = 0.010,
) -> Dict:
    """Fit both colour indices, invert to (beta, alpha), and propagate uncertainty.

    `frames` needs columns: airmass, CI_RB, CI_GB, qc_ok.

    Uncertainty = frame bootstrap (random error) combined with Monte Carlo
    systematic error: each draw scales both slopes by independent N(1,
    slope_sys_rel) factors (tone curve, JPEG compression, residual
    non-linearity) and perturbs the effective channel wavelengths by
    N(0, wavelength_sd_um). On synthetic data, JPEG output biased slopes by
    about 2 %, RAW/PNG by < 0.5 %; the bootstrap alone cannot see that.
    """
    if slope_sys_rel < 0 or wavelength_sd_um < 0:
        raise ValueError("slope_sys_rel and wavelength_sd_um must be >= 0.")
    required = {"airmass", "CI_RB", "CI_GB", "qc_ok"}
    missing = required - set(frames.columns)
    if missing:
        raise KeyError(f"frames is missing columns: {sorted(missing)}")
    wl = dict(DEFAULT_WAVELENGTHS_UM if wavelengths_um is None else wavelengths_um)

    good = frames[frames["qc_ok"].astype(bool)]
    m = good["airmass"].to_numpy(float)
    ci_rb = good["CI_RB"].to_numpy(float)
    ci_gb = good["CI_GB"].to_numpy(float)

    fit_rb = fit_colour_index(m, ci_rb, airmass_range)
    fit_gb = fit_colour_index(m, ci_gb, airmass_range)

    clean_rb = colour_index_slope(wl["R"], wl["B"], 0.0, 0.0, pressure_hpa)
    clean_gb = colour_index_slope(wl["G"], wl["B"], 0.0, 0.0, pressure_hpa)

    summary = {
        "n_frames_total": int(len(frames)),
        "n_frames_qc_ok": int(len(good)),
        "airmass_range": list(airmass_range),
        "pressure_hpa": float(pressure_hpa),
        "wavelengths_um": wl,
        "fit_RB": {k: v for k, v in fit_rb.items() if k != "mask"},
        "fit_GB": {k: v for k, v in fit_gb.items() if k != "mask"},
        "clean_sky_slope_RB": clean_rb,
        "clean_sky_slope_GB": clean_gb,
        "systematics": {"slope_sys_rel": float(slope_sys_rel), "wavelength_sd_um": float(wavelength_sd_um)},
    }
    if fit_rb["status"] != "ok" or fit_gb["status"] != "ok":
        summary["retrieval"] = {"status": f"fit failed: RB={fit_rb['status']}; GB={fit_gb['status']}"}
        return summary

    # 1) Robust, model-light product: aerosol colour slope = measured - clean sky.
    a_rb = fit_rb["slope"] - clean_rb
    summary["aerosol_colour_slope_RB"] = a_rb
    summary["excess_reddening_RB_pct"] = 100 * a_rb / clean_rb

    # 2) AOD(550) conditional on an assumed alpha (colour cannot see grey aerosol).
    summary["aod_550_given_alpha"] = {
        f"{a:g}": aod_from_aerosol_slope(a_rb, a, wl["R"], wl["B"]) for a in ALPHA_GRID
    }

    # 3) Free three-channel (beta, alpha) inversion: ill-conditioned, so judged
    #    by its own bootstrap spread before it is trusted.
    point = retrieve_beta_alpha(fit_rb["slope"], fit_gb["slope"], wl, pressure_hpa, alpha_fixed)
    summary["retrieval"] = point.to_dict()

    # Bootstrap over frames (keeps RB/GB pairing so their correlation is kept).
    rng = np.random.default_rng(seed)
    idx_all = np.flatnonzero(fit_rb["mask"] & fit_gb["mask"])
    draws = {"slope_RB": [], "slope_GB": [], "aerosol_colour_slope_RB": [],
             "beta": [], "alpha": [], "aod_550": []}
    grid_draws = {f"{a:g}": [] for a in ALPHA_GRID}
    n_ok = 0
    for _ in range(int(n_boot)):
        idx = rng.choice(idx_all, size=idx_all.size, replace=True)
        if np.ptp(m[idx]) < 1.0:
            continue
        k_rb = np.polyfit(m[idx], ci_rb[idx], 1)[0] * rng.normal(1.0, slope_sys_rel)
        k_gb = np.polyfit(m[idx], ci_gb[idx], 1)[0] * rng.normal(1.0, slope_sys_rel)
        wl_d = {c: wl[c] + rng.normal(0.0, wavelength_sd_um) for c in "RGB"}
        if not (wl_d["R"] > wl_d["G"] > wl_d["B"] > 0.3):
            continue
        clean_d = colour_index_slope(wl_d["R"], wl_d["B"], 0.0, 0.0, pressure_hpa)
        draws["slope_RB"].append(k_rb)
        draws["slope_GB"].append(k_gb)
        draws["aerosol_colour_slope_RB"].append(k_rb - clean_d)
        for a in ALPHA_GRID:
            grid_draws[f"{a:g}"].append(aod_from_aerosol_slope(k_rb - clean_d, a, wl_d["R"], wl_d["B"]))
        r = retrieve_beta_alpha(k_rb, k_gb, wl_d, pressure_hpa, alpha_fixed)
        if r.status.startswith("ok"):
            n_ok += 1
            draws["beta"].append(r.beta)
            draws["alpha"].append(r.alpha)
            draws["aod_550"].append(r.aod_550)

    def pct(v):
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v.size >= 20 else [np.nan, np.nan]

    ci95 = {k: pct(v) for k, v in draws.items()}
    frac_ok = n_ok / max(int(n_boot), 1)
    summary["bootstrap"] = {
        "n_boot": int(n_boot),
        "fraction_successful": frac_ok,
        "ci95": ci95,
        "ci95_aod_550_given_alpha": {k: pct(v) for k, v in grid_draws.items()},
    }

    # Reliability flag for the free inversion.
    width = ci95["alpha"][1] - ci95["alpha"][0] if np.all(np.isfinite(ci95["alpha"])) else np.inf
    if alpha_fixed is not None:
        summary["retrieval"]["reliability"] = "alpha fixed by user"
    elif frac_ok >= 0.9 and width <= ALPHA_CI_WIDTH_OK:
        summary["retrieval"]["reliability"] = "constrained"
    else:
        summary["retrieval"]["reliability"] = (
            f"poorly constrained (alpha 95% CI width {width:.2f}, {100 * frac_ok:.0f}% of bootstrap "
            f"draws invertible); report AOD for an assumed alpha instead"
        )
    return summary
