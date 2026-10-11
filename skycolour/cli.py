"""Command line: `python -m skycolour process ...` and `python -m skycolour synth ...`."""

from __future__ import annotations

import argparse
import json
import sys

from .fit import DEFAULT_AIRMASS_RANGE
from .physics import STANDARD_PRESSURE_HPA, pressure_from_altitude


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="skycolour",
        description="Retrieve aerosol load from a sunset or sunrise photo sequence "
                    "(colour index vs air mass), or generate a synthetic test sequence.")
    sub = p.add_subparsers(dest="cmd", required=True)

    pr = sub.add_parser("process", help="Retrieve aerosol properties from a folder of photos.")
    pr.add_argument("folder", help="Folder with one sunset/sunrise sequence (JPEG, TIFF, PNG or RAW).")
    pr.add_argument("--out", default="results", help="Output folder (default: results).")
    pr.add_argument("--lat", type=float, help="Latitude, + north (overrides GPS EXIF).")
    pr.add_argument("--lon", type=float, help="Longitude, + east (overrides GPS EXIF).")
    pr.add_argument("--tz", help="IANA time zone of the camera clock, e.g. America/Toronto, Europe/Paris, Asia/Beirut.")
    pr.add_argument("--manifest", help="CSV with columns file,time[,lat,lon] (overrides EXIF).")
    pr.add_argument("--clock-offset-s", type=float, default=0.0,
                    help="Seconds to ADD to camera time to get true time (camera slow -> positive).")
    g = pr.add_mutually_exclusive_group()
    g.add_argument("--pressure-hpa", type=float, help=f"Surface pressure (default {STANDARD_PRESSURE_HPA}).")
    g.add_argument("--altitude-m", type=float, help="Site altitude; converted to standard-atmosphere pressure.")
    pr.add_argument("--airmass-min", type=float, default=DEFAULT_AIRMASS_RANGE[0])
    pr.add_argument("--airmass-max", type=float, default=DEFAULT_AIRMASS_RANGE[1])
    pr.add_argument("--alpha-fixed", type=float, help="Fix the Angstrom exponent and solve beta only.")
    pr.add_argument("--wavelengths", nargs=3, type=float, metavar=("R", "G", "B"),
                    help="Effective channel wavelengths in um (default 0.610 0.535 0.465).")
    pr.add_argument("--roi", nargs=4, type=int, metavar=("X0", "Y0", "X1", "Y1"),
                    help="Only search for the sun inside this pixel box.")
    pr.add_argument("--n-boot", type=int, default=1000)
    pr.add_argument("--slope-sys", type=float,
                    help="Relative systematic slope error (default 0.02 for JPEG/PNG input, 0.005 for RAW).")
    pr.add_argument("--wavelength-sd", type=float, default=0.010,
                    help="Uncertainty (um) of the effective channel wavelengths (default 0.010; "
                         "lower it once your camera is calibrated).")
    pr.add_argument("--quiet", action="store_true")

    sy = sub.add_parser("synth", help="Write a synthetic sunset sequence with known aerosol.")
    sy.add_argument("out", help="Output folder for the images.")
    sy.add_argument("--beta", type=float, default=0.10)
    sy.add_argument("--alpha", type=float, default=1.3)
    sy.add_argument("--lat", type=float, default=45.5017)
    sy.add_argument("--lon", type=float, default=-73.5673)
    sy.add_argument("--date", default="2026-10-09")
    sy.add_argument("--tz", default="America/Toronto")
    sy.add_argument("--interval-s", type=int, default=60)
    sy.add_argument("--noise", type=float, default=0.003)
    sy.add_argument("--fmt", choices=("jpg", "png"), default="jpg",
                    help="png mimics linear RAW-quality data; jpg includes compression artefacts.")
    return p


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    if args.cmd == "synth":
        from .synthetic import make_sunset_sequence
        gt = make_sunset_sequence(args.out, args.beta, args.alpha, args.lat, args.lon, args.date,
                                  args.tz, interval_s=args.interval_s, noise_rel=args.noise, fmt=args.fmt)
        print(f"Wrote {len(gt)} frames to {args.out} (truth: beta={args.beta}, alpha={args.alpha}).")
        return 0

    from .pipeline import run
    pressure = (pressure_from_altitude(args.altitude_m) if args.altitude_m is not None
                else args.pressure_hpa or STANDARD_PRESSURE_HPA)
    wl = dict(zip("RGB", args.wavelengths)) if args.wavelengths else None
    summary = run(
        args.folder, args.out, lat=args.lat, lon=args.lon, tz=args.tz, manifest=args.manifest,
        clock_offset_s=args.clock_offset_s, pressure_hpa=pressure,
        airmass_range=(args.airmass_min, args.airmass_max), alpha_fixed=args.alpha_fixed,
        wavelengths_um=wl, roi=tuple(args.roi) if args.roi else None, n_boot=args.n_boot,
        slope_sys_rel=args.slope_sys,
        wavelength_sd_um=args.wavelength_sd, verbose=not args.quiet,
    )
    print(format_report(summary))
    print(f"Outputs in {args.out}/ (frames.csv, summary.json, colour_index_vs_airmass.png)")
    return 0


def _ci(v, fmt="{:.3f}"):
    try:
        return f"[{fmt.format(v[0])}, {fmt.format(v[1])}]"
    except (TypeError, ValueError, IndexError):
        return "[n/a]"


def format_report(summary: dict) -> str:
    """Plain-text summary, most robust quantity first."""
    lines = [f"\nFrames: {summary['n_frames_qc_ok']}/{summary['n_frames_total']} passed QC"]
    fr = summary.get("fit_RB", {})
    if "aerosol_colour_slope_RB" not in summary:
        lines.append(f"Fit failed: {summary.get('retrieval', {}).get('status')}")
        return "\n".join(lines)
    b = summary["bootstrap"]
    lines += [
        f"R/B slope vs air mass : {fr['slope']:.4f} ± {fr['slope_se']:.4f}  (R² = {fr['r2']:.3f}, n = {fr['n']})",
        f"Clean-sky R/B slope   : {summary['clean_sky_slope_RB']:.4f}",
        f"Aerosol colour slope  : {summary['aerosol_colour_slope_RB']:.4f}  95% CI {_ci(b['ci95']['aerosol_colour_slope_RB'], '{:.4f}')}"
        f"  (+{summary['excess_reddening_RB_pct']:.0f} % reddening vs clean sky)",
        "AOD(550 nm) for an assumed Angstrom exponent:",
    ]
    for a, v in summary["aod_550_given_alpha"].items():
        lines.append(f"    alpha = {a:>4}: {v:.3f}  95% CI {_ci(b['ci95_aod_550_given_alpha'][a])}")
    r = summary["retrieval"]
    lines.append(f"Free 3-channel inversion: {r.get('status')}; {r.get('reliability', '')}")
    if str(r.get("status", "")).startswith("ok"):
        ci = b["ci95"]
        lines.append(f"    alpha = {r['alpha']:.2f} {_ci(ci['alpha'], '{:.2f}')},  "
                     f"AOD(550) = {r['aod_550']:.3f} {_ci(ci['aod_550'])}")
    return "\n".join(lines)


if __name__ == "__main__":
    sys.exit(main())
