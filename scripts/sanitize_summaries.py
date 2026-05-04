#!/usr/bin/env python3
"""Batch MAD-any-of sanitization over MM-GBSA summary JSONs.

Walks `outputs/*_calib_s*/mmgbsa_results*/mmgbsa_summary.json` and writes
`mmgbsa_summary_sanitized.json` alongside each original. Non-destructive:
original files are never modified.

Schema contract (per sanitized file):
  * top-level `sanitization_policy` = "mad_any_of_v1"
  * top-level `sanitization` audit block (see utils/mmgbsa_outlier_screen.py)

Usage:
    python3 scripts/sanitize_summaries.py            # full run (writes files)
    python3 scripts/sanitize_summaries.py --dry_run  # report only, no writes
    python3 scripts/sanitize_summaries.py --z 3.5    # custom z-threshold
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import List

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
UTILS_DIR = os.path.join(REPO_ROOT, "utils")
if UTILS_DIR not in sys.path:
    sys.path.insert(0, UTILS_DIR)

from mmgbsa_outlier_screen import (  # noqa: E402
    DEFAULT_Z_THRESHOLD,
    SANITIZATION_POLICY_ID,
    sanitize_mmgbsa_summary,
)

DEFAULT_GLOB = os.path.join(
    REPO_ROOT,
    "outputs",
    "*_calib_s*",
    "mmgbsa_results*",
    "mmgbsa_summary.json",
)
LARGE_SHIFT_KCAL = 5.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--glob",
        default=DEFAULT_GLOB,
        help="Glob pattern for mmgbsa_summary.json inputs.",
    )
    p.add_argument(
        "--z",
        type=float,
        default=DEFAULT_Z_THRESHOLD,
        help=f"|z| threshold (default {DEFAULT_Z_THRESHOLD}).",
    )
    p.add_argument(
        "--dry_run",
        action="store_true",
        help="Do not write sanitized JSONs; only report aggregate stats.",
    )
    return p.parse_args()


def _run_tag(summary_path: str) -> str:
    # outputs/<run_tag>/<mmgbsa_results_dir>/mmgbsa_summary.json
    parent = os.path.dirname(summary_path)  # .../mmgbsa_results[_n25]
    run_dir = os.path.dirname(parent)  # .../<run_tag>
    results_sub = os.path.basename(parent)
    return f"{os.path.basename(run_dir)}::{results_sub}"


def _fmt(v):
    if v is None:
        return "None"
    return f"{v:.3f}"


def main() -> int:
    args = parse_args()
    paths = sorted(glob.glob(args.glob))
    if not paths:
        print(f"[sanitize_summaries] no summaries match glob: {args.glob}")
        return 1

    total_snaps = 0
    total_flagged = 0
    large_shifts: List[str] = []
    per_run_rows: List[tuple] = []
    errors: List[str] = []

    for path in paths:
        tag = _run_tag(path)
        try:
            sanitized = sanitize_mmgbsa_summary(
                path,
                z_threshold=args.z,
                write_sanitized=not args.dry_run,
            )
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{tag}: {type(exc).__name__}: {exc}")
            continue

        audit = sanitized["sanitization"]
        n_orig = audit["n_snap_original"]
        n_flag = audit["n_snap_flagged"]
        mean_orig = audit["mean_dg_original"]
        mean_clean = audit["mean_dg_sanitized"]
        delta = audit["delta_mean_dg"]
        total_snaps += n_orig
        total_flagged += n_flag
        if delta is not None and abs(delta) > LARGE_SHIFT_KCAL:
            large_shifts.append(
                f"{tag}: delta_mean_dg={delta:+.3f} "
                f"(orig={_fmt(mean_orig)} -> clean={_fmt(mean_clean)}, "
                f"{n_flag}/{n_orig} dropped)"
            )
        per_run_rows.append(
            (tag, n_orig, n_flag, mean_orig, mean_clean, delta)
        )

    print("=" * 90)
    print(
        f"MAD-any-of sanitization ({SANITIZATION_POLICY_ID}) "
        f"| z_threshold={args.z} | dry_run={args.dry_run}"
    )
    print(f"Glob: {args.glob}")
    print(f"Summaries scanned: {len(paths)}")
    print(f"Total snaps: {total_snaps}  flagged: {total_flagged}")
    print("=" * 90)

    if large_shifts:
        print(f"\nRuns with |delta_mean_dg| > {LARGE_SHIFT_KCAL:g} kcal/mol "
              f"({len(large_shifts)}):")
        for line in large_shifts:
            print(f"  * {line}")
    else:
        print("\nNo run showed |delta_mean_dg| > "
              f"{LARGE_SHIFT_KCAL:g} kcal/mol.")

    print("\nPer-run summary (tag | n_snap | n_flag | mean_orig | mean_clean | delta):")
    for tag, n_orig, n_flag, m_o, m_c, d in per_run_rows:
        print(
            f"  {tag:<60} {n_orig:>3} {n_flag:>3} "
            f"{_fmt(m_o):>12} {_fmt(m_c):>12} {_fmt(d):>10}"
        )

    if errors:
        print("\nErrors:")
        for e in errors:
            print(f"  ! {e}")
        return 2

    if not args.dry_run:
        print("\nSanitized JSONs written alongside originals as "
              "*_sanitized.json (non-destructive).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
