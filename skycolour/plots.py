"""Diagnostic plot: colour index vs air mass, measured fit vs clean-sky baseline."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INK, MUTED, MEASURED, BASELINE = "#1f2328", "#9aa0a6", "#c2410c", "#2563eb"


def plot_colour_index(frames: pd.DataFrame, summary: dict, path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharex=True)
    for ax, col, key, clean_key in (
        (axes[0], "CI_RB", "fit_RB", "clean_sky_slope_RB"),
        (axes[1], "CI_GB", "fit_GB", "clean_sky_slope_GB"),
    ):
        ok = frames["qc_ok"].astype(bool) & frames[col].notna()
        bad = ~frames["qc_ok"].astype(bool) & frames[col].notna() & frames["airmass"].notna()
        if bad.any():
            ax.scatter(frames.loc[bad, "airmass"], frames.loc[bad, col], s=14, c=MUTED, label="rejected (QC)")
        lo, hi = summary.get("airmass_range", [0, np.inf])
        inside = ok & frames["airmass"].between(lo, hi)
        ax.scatter(frames.loc[inside, "airmass"], frames.loc[inside, col], s=18, c=INK, label="frames used in fit")
        if (ok & ~inside).any():
            ax.scatter(frames.loc[ok & ~inside, "airmass"], frames.loc[ok & ~inside, col], s=18,
                       facecolors="none", edgecolors=INK, label="outside air-mass window")
        fit = summary.get(key, {})
        if fit.get("status") == "ok":
            m_max = np.nanmax(frames["airmass"].to_numpy(float)) if frames["airmass"].notna().any() else 1.0
            m = np.linspace(0, max(m_max, 1.0) * 1.02, 50)
            ax.plot(m, fit["intercept"] + fit["slope"] * m, c=MEASURED, lw=2,
                    label=f"measured k = {fit['slope']:.4f} ± {fit['slope_se']:.4f}")
            ax.plot(m, fit["intercept"] + summary[clean_key] * m, c=BASELINE, lw=1.5, ls="--",
                    label=f"clean sky (Rayleigh only) k = {summary[clean_key]:.4f}")
            ax.axvspan(lo, hi, color=MUTED, alpha=0.08, lw=0)
        ax.set_xlabel("relative air mass m (–)")
        ax.set_ylabel(f"ln({col[3]}/{col[4]}) (–)")
        ax.set_title(f"{col[3]}/{col[4]} colour index", loc="left", fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8, loc="upper left")

    r = summary.get("retrieval", {})
    if "aerosol_colour_slope_RB" in summary:
        g = summary["aod_550_given_alpha"]
        text = (f"aerosol colour slope (R/B) = {summary['aerosol_colour_slope_RB']:.4f}  "
                f"(+{summary['excess_reddening_RB_pct']:.0f} % vs clean sky)   "
                f"AOD(550) ≈ {g.get('1', np.nan):.2f} if α=1, {g.get('1.5', np.nan):.2f} if α=1.5")
        if r.get("reliability") == "constrained":
            text += f"   free fit: α = {r['alpha']:.2f}, AOD = {r['aod_550']:.3f}"
    else:
        text = f"retrieval: {r.get('status', 'not run')}"
    fig.suptitle(text, x=0.01, ha="left", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)
