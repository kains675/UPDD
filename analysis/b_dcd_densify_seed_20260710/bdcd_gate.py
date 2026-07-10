#!/usr/bin/env python
"""DCD structural gate for the B densify + seed expansion campaign."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


DEFAULT_PAIRS = [
    {
        "name": "w23a_bound_dminus",
        "kind": "bound_dminus",
        "high": "w23a_s101_bound_dminus",
        "control": "w23a_ctrl_s42_bound_dminus",
    },
    {
        "name": "w4a_carved_bound_dminus",
        "kind": "bound_dminus",
        "high": "w4a_carved_s127_bound_dminus",
        "control": "w4a_carved_ctrl_s163_bound_dminus",
    },
    {
        "name": "w4a_carved_free_dplus",
        "kind": "free_dplus",
        "high": "w4a_carved_s127_free_dplus",
        "control": "w4a_carved_ctrl_s101_free_dplus",
    },
]


def _load_cells(path: Path) -> dict[str, dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    cells = payload.get("cells")
    if not isinstance(cells, list):
        raise ValueError(f"{path}: expected top-level cells list")
    out: dict[str, dict[str, Any]] = {}
    for cell in cells:
        cid = cell.get("cell_id")
        if not cid:
            raise ValueError(f"{path}: cell without cell_id")
        out[str(cid)] = cell
    return out


def _frac(cell: dict[str, Any]) -> float:
    wm = cell.get("water_min_nm") or {}
    if "frac_lt_0p26" not in wm:
        raise ValueError(f"{cell.get('cell_id')}: missing water_min_nm.frac_lt_0p26")
    return float(wm["frac_lt_0p26"])


def evaluate(
    cells: dict[str, dict[str, Any]],
    *,
    phase: str,
    baseline_floor: float,
    amp_delta: float,
    amp_ratio: float,
) -> dict[str, Any]:
    pairs = []
    post_densify_a_triggers = []

    for pair in DEFAULT_PAIRS:
        high = cells[pair["high"]]
        control = cells[pair["control"]]
        high_frac = _frac(high)
        control_frac = _frac(control)
        delta = high_frac - control_frac
        ratio = high_frac / control_frac if control_frac > 0 else float("inf")
        baseline = pair["kind"] == "bound_dminus" and control_frac >= baseline_floor
        amplified = (
            pair["kind"] == "bound_dminus"
            and delta >= amp_delta
            and ratio >= amp_ratio
        )
        post_trigger = phase == "post_densify" and (baseline or amplified)
        if post_trigger:
            post_densify_a_triggers.append(pair["name"])
        pairs.append(
            {
                "name": pair["name"],
                "kind": pair["kind"],
                "high_cell": pair["high"],
                "control_cell": pair["control"],
                "high_frac_lt_0p26": high_frac,
                "control_frac_lt_0p26": control_frac,
                "delta": delta,
                "ratio": ratio,
                "control_baseline_persists": baseline,
                "high_leverage_amplified": amplified,
                "post_densify_a_trigger": post_trigger,
            }
        )

    if phase == "post_densify" and post_densify_a_triggers:
        verdict = "PROMOTE_A_CARVE_REDESIGN"
    elif phase == "post_densify":
        verdict = "CONTINUE_B_NO_A_TRIGGER"
    else:
        verdict = "PRECHECK_BASELINE_RECORDED"

    return {
        "phase": phase,
        "thresholds": {
            "baseline_floor": baseline_floor,
            "amp_delta": amp_delta,
            "amp_ratio": amp_ratio,
        },
        "pairs": pairs,
        "post_densify_a_triggers": post_densify_a_triggers,
        "verdict": verdict,
        "note": (
            "D2/precheck can record risk only. A-promotion is active only for "
            "phase=post_densify, after the densify+seed expansion campaign."
        ),
    }


def write_report(result: dict[str, Any], path: Path) -> None:
    lines = [
        "# B-DCD Structural Gate Report",
        "",
        f"Phase: `{result['phase']}`",
        f"Verdict: `{result['verdict']}`",
        "",
        "| pair | high frac<0.26 | control frac<0.26 | delta | ratio | baseline | amplified | A trigger |",
        "| --- | ---: | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for pair in result["pairs"]:
        lines.append(
            "| {name} | {hf:.3f} | {cf:.3f} | {delta:.3f} | {ratio:.2f} | {base} | {amp} | {trig} |".format(
                name=pair["name"],
                hf=pair["high_frac_lt_0p26"],
                cf=pair["control_frac_lt_0p26"],
                delta=pair["delta"],
                ratio=pair["ratio"],
                base=pair["control_baseline_persists"],
                amp=pair["high_leverage_amplified"],
                trig=pair["post_densify_a_trigger"],
            )
        )
    lines.extend(["", result["note"], ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument(
        "--phase",
        choices=("precheck", "post_densify"),
        default="precheck",
    )
    parser.add_argument("--baseline-floor", type=float, default=0.15)
    parser.add_argument("--amp-delta", type=float, default=0.10)
    parser.add_argument("--amp-ratio", type=float, default=1.5)
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument("--report-out", type=Path, default=None)
    args = parser.parse_args()

    cells = _load_cells(args.summary)
    result = evaluate(
        cells,
        phase=args.phase,
        baseline_floor=args.baseline_floor,
        amp_delta=args.amp_delta,
        amp_ratio=args.amp_ratio,
    )

    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.json_out:
        args.json_out.write_text(text + "\n", encoding="utf-8")
    if args.report_out:
        write_report(result, args.report_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
