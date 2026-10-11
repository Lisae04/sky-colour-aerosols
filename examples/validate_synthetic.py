"""
Synthetic validation of Stage A: known (beta, alpha) -> images -> retrieval.

Self-contained: run from the repo root with
    python examples/validate_synthetic.py
Writes docs/validation.csv and prints a markdown table.

What it measures
----------------
* slope_err_pct   : error of the measured R/B colour-index slope (%)
* aer_slope_in_ci : is the true aerosol colour slope inside the reported 95 % CI?
* aod_in_ci       : is the true AOD(550) inside the CI for the TRUE alpha
                    (i.e. if alpha were known, e.g. from AERONET)?
* alpha_hat, alpha_ci : free 3-channel inversion and its 95 % CI
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from skycolour.fit import retrieve_from_frames  # noqa: E402
from skycolour.physics import aod_from_aerosol_slope, colour_index_slope  # noqa: E402
from skycolour.pipeline import process_folder  # noqa: E402
from skycolour.synthetic import make_sunset_sequence  # noqa: E402

CASES = [(b, a) for a in (0.5, 1.3, 2.0) for b in (0.05, 0.15)]
INPUTS = {  # format, noise, systematic slope error, wavelength sd (um)
    "PNG, calibrated camera": ("png", 0.003, 0.005, 0.0),
    "JPEG, calibrated camera": ("jpg", 0.003, 0.02, 0.0),
    "JPEG, uncalibrated camera": ("jpg", 0.003, 0.02, 0.010),
}


def main():
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for label, (fmt, noise, sys_rel, wl_sd) in INPUTS.items():
            for i, (beta, alpha) in enumerate(CASES):
                d = Path(tmp) / f"{fmt}_{i}"
                if not d.exists():
                    make_sunset_sequence(d, beta, alpha, fmt=fmt, noise_rel=noise, interval_s=90, seed=i)
                frames = process_folder(d, verbose=False)
                s = retrieve_from_frames(frames, n_boot=400, slope_sys_rel=sys_rel,
                                         wavelength_sd_um=wl_sd, seed=i)
                true_k = colour_index_slope(0.610, 0.465, beta, alpha)
                true_a = true_k - s["clean_sky_slope_RB"]
                true_aod = beta * 0.55 ** -alpha
                lo, hi = s["bootstrap"]["ci95"]["aerosol_colour_slope_RB"]
                # AOD for the true alpha, with CI scaled from the aerosol-slope CI.
                aod_hat = aod_from_aerosol_slope(s["aerosol_colour_slope_RB"], alpha)
                aod_lo, aod_hi = (aod_from_aerosol_slope(v, alpha) for v in (lo, hi))
                r = s["retrieval"]
                aci = s["bootstrap"]["ci95"]["alpha"]
                rows.append({
                    "input": label, "beta": beta, "alpha": alpha, "true_aod_550": round(true_aod, 3),
                    "slope_err_pct": round(100 * (s["fit_RB"]["slope"] / true_k - 1), 2),
                    "aer_slope_in_ci": bool(lo <= true_a <= hi),
                    "aod_550_hat_known_alpha": round(aod_hat, 3),
                    "aod_in_ci": bool(aod_lo <= true_aod <= aod_hi),
                    "alpha_hat": round(r.get("alpha", np.nan), 2),
                    "alpha_ci": f"[{aci[0]:.2f}, {aci[1]:.2f}]",
                    "alpha_reliability": r.get("reliability", r.get("status")),
                })
    df = pd.DataFrame(rows)
    out = ROOT / "docs" / "validation.csv"
    out.parent.mkdir(exist_ok=True)
    df.to_csv(out, index=False)

    summary = df.groupby("input", sort=False).agg(
        mean_abs_slope_err_pct=("slope_err_pct", lambda v: round(np.abs(v).mean(), 2)),
        aerosol_slope_ci_coverage=("aer_slope_in_ci", "mean"),
        aod_ci_coverage_known_alpha=("aod_in_ci", "mean"),
        alpha_constrained=("alpha_reliability", lambda v: (v == "constrained").mean()),
    )
    print(df.to_markdown(index=False))
    print()
    print(summary.to_markdown())
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
