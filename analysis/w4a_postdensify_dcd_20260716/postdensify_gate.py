#!/usr/bin/env python3
"""Evaluate the frozen W4A post-densify DCD structural gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import median
from typing import Any


ANALYSIS_DIR = Path(__file__).resolve().parent
DEFAULT_SUMMARY = ANALYSIS_DIR / "postdensify_dcd_summary.json"
DEFAULT_JSON_OUT = ANALYSIS_DIR / "postdensify_gate.json"
DEFAULT_REPORT_OUT = ANALYSIS_DIR / "postdensify_gate.md"

SEEDS = ("s7", "s19", "s23", "s42", "s83", "s101", "s127", "s163", "s199", "s251")
HIGH_LEVERAGE_SEED = "s127"
BASELINE_FLOOR = 0.15
AMPLIFICATION_DELTA = 0.10
AMPLIFICATION_RATIO = 1.5
EXPECTED_FRAMES = 750


def frac_lt_0p26(cell: dict[str, Any]) -> float:
    water = cell.get("water_min_nm") or {}
    value = water.get("frac_lt_0p26")
    if value is None:
        raise ValueError(f"{cell.get('cell_id')}: missing frac_lt_0p26")
    return float(value)


def validate_summary(summary: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    if summary.get("status") != "COMPLETE":
        raise ValueError(f"summary status is {summary.get('status')!r}, expected 'COMPLETE'")
    cells = summary.get("cells")
    if not isinstance(cells, list) or len(cells) != 20:
        raise ValueError("gate requires the complete preregistered 20-cell cohort")

    indexed: dict[tuple[str, str], dict[str, Any]] = {}
    for cell in cells:
        arm = cell.get("arm")
        seed = cell.get("seed")
        key = (str(arm), str(seed))
        if key in indexed:
            raise ValueError(f"duplicate cell {key}")
        if arm not in {"carved", "uncarved"} or seed not in SEEDS:
            raise ValueError(f"unexpected cell {key}")
        if cell.get("leg") != "bound" or cell.get("direction") != "dminus":
            raise ValueError(f"{key}: expected bound/dminus")
        if cell.get("n_frames") != EXPECTED_FRAMES:
            raise ValueError(f"{key}: expected {EXPECTED_FRAMES} DCD frames")
        if cell.get("nan_any") is not False:
            raise ValueError(f"{key}: NaN status is not false")
        if cell.get("n_local_water_oxygens", 0) <= 0:
            raise ValueError(f"{key}: no selected local water oxygens")
        frac_lt_0p26(cell)
        indexed[key] = cell

    expected = {(arm, seed) for arm in ("carved", "uncarved") for seed in SEEDS}
    if set(indexed) != expected:
        missing = sorted(expected - set(indexed))
        extra = sorted(set(indexed) - expected)
        raise ValueError(f"cohort mismatch: missing={missing}, extra={extra}")
    return indexed


def evaluate(summary: dict[str, Any]) -> dict[str, Any]:
    cells = validate_summary(summary)
    carved = {seed: frac_lt_0p26(cells[("carved", seed)]) for seed in SEEDS}
    uncarved = {seed: frac_lt_0p26(cells[("uncarved", seed)]) for seed in SEEDS}

    controls = [carved[seed] for seed in SEEDS if seed != HIGH_LEVERAGE_SEED]
    control_median = float(median(controls))
    high = carved[HIGH_LEVERAGE_SEED]
    amplification_delta = high - control_median
    amplification_ratio = high / control_median if control_median > 0 else None

    baseline_persists = control_median >= BASELINE_FLOOR
    amplification_persists = (
        amplification_ratio is not None
        and amplification_delta >= AMPLIFICATION_DELTA
        and amplification_ratio >= AMPLIFICATION_RATIO
    )
    promote_a = baseline_persists or amplification_persists
    paired = {seed: carved[seed] - uncarved[seed] for seed in SEEDS}

    return {
        "schema": "w4a_postdensify_dcd_gate_v1",
        "source_inventory_digest": summary.get("source_inventory_digest"),
        "thresholds": {
            "baseline_floor": BASELINE_FLOOR,
            "amplification_delta": AMPLIFICATION_DELTA,
            "amplification_ratio": AMPLIFICATION_RATIO,
        },
        "high_leverage": {
            "seed": HIGH_LEVERAGE_SEED,
            "frac_lt_0p26": high,
        },
        "carved_controls": {
            "seeds": [seed for seed in SEEDS if seed != HIGH_LEVERAGE_SEED],
            "values": {seed: carved[seed] for seed in SEEDS if seed != HIGH_LEVERAGE_SEED},
            "median": control_median,
        },
        "criteria": {
            "control_baseline_persists": baseline_persists,
            "high_minus_control_median": amplification_delta,
            "high_over_control_median": amplification_ratio,
            "high_leverage_amplification_persists": amplification_persists,
        },
        "secondary_context": {
            "uncarved_median": float(median(uncarved.values())),
            "uncarved_values": uncarved,
            "paired_carved_minus_uncarved": paired,
            "paired_difference_median": float(median(paired.values())),
        },
        "verdict": "PROMOTE_A_CARVE_REDESIGN" if promote_a else "CONTINUE_B_TRANSPORT_PILOT",
        "interpretation_limit": (
            "Short-probe walker/cycle frames are correlated. Values are structural propensity diagnostics, "
            "not independent observations or free-energy estimates."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    criteria = result["criteria"]
    context = result["secondary_context"]
    ratio = criteria["high_over_control_median"]
    ratio_text = "undefined" if ratio is None else f"{ratio:.3f}"
    lines = [
        "# W4A Post-densify DCD Structural Gate",
        "",
        f"Verdict: `{result['verdict']}`",
        "",
        f"- carved s127 frac<0.26: `{result['high_leverage']['frac_lt_0p26']:.3f}`",
        f"- other-nine carved control median: `{result['carved_controls']['median']:.3f}`",
        f"- high minus control median: `{criteria['high_minus_control_median']:.3f}`",
        f"- high/control ratio: `{ratio_text}`",
        f"- control baseline persists: `{criteria['control_baseline_persists']}`",
        f"- high-leverage amplification persists: `{criteria['high_leverage_amplification_persists']}`",
        f"- uncarved median, secondary only: `{context['uncarved_median']:.3f}`",
        f"- paired carved-minus-uncarved median, secondary only: `{context['paired_difference_median']:.3f}`",
        "",
        result["interpretation_limit"],
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON_OUT)
    parser.add_argument("--report-out", type=Path, default=DEFAULT_REPORT_OUT)
    args = parser.parse_args(argv)

    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    result = evaluate(summary)
    encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    args.json_out.write_text(encoded, encoding="utf-8")
    write_report(result, args.report_out)
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
