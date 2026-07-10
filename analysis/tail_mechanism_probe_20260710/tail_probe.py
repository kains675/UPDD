#!/usr/bin/env python3
"""Read-only Track B tail-mechanism probe over existing W23A/W4A artifacts."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent

W23A_ROOT = REPO_ROOT / "outputs/_trackb/mdm2_w23a_gateA_20260708"
W4A_C1_ROOT = REPO_ROOT / "outputs/_trackb/w4a_carved_c1_20260709"

TAIL_THRESHOLDS = (120.0, 150.0, 180.0)


@dataclass(frozen=True)
class Cohort:
    name: str
    label: str
    roots_by_leg: dict[str, Path]
    seeds: list[str]
    ddg_by_seed: dict[str, float]
    dgb_by_leg_seed: dict[str, dict[str, float]]
    highlight_seeds: list[str]
    json_path: Path


class SampleBag:
    def __init__(self) -> None:
        self.values: list[float] = []

    def add(self, value: float) -> None:
        self.values.append(value)

    def summary(self) -> dict[str, Any]:
        values = sorted(self.values)
        n = len(values)
        if n == 0:
            return {
                "n": 0,
                "mean": None,
                "sd": None,
                "min": None,
                "p01": None,
                "p05": None,
                "p50": None,
                "p95": None,
                "p99": None,
                "max": None,
                "frac_gt_120": None,
                "frac_gt_150": None,
                "frac_gt_180": None,
                "frac_lt_0": None,
            }
        mean = sum(values) / n
        sd = math.sqrt(sum((v - mean) ** 2 for v in values) / (n - 1)) if n > 1 else 0.0
        return {
            "n": n,
            "mean": mean,
            "sd": sd,
            "min": values[0],
            "p01": percentile(values, 0.01),
            "p05": percentile(values, 0.05),
            "p50": percentile(values, 0.50),
            "p95": percentile(values, 0.95),
            "p99": percentile(values, 0.99),
            "max": values[-1],
            "frac_gt_120": fraction_gt(values, 120.0),
            "frac_gt_150": fraction_gt(values, 150.0),
            "frac_gt_180": fraction_gt(values, 180.0),
            "frac_lt_0": sum(1 for v in values if v < 0.0) / n,
        }


def percentile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        raise ValueError("percentile requires at least one value")
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = q * (len(sorted_values) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return sorted_values[low]
    weight = position - low
    return sorted_values[low] * (1.0 - weight) + sorted_values[high] * weight


def fraction_gt(sorted_values: list[float], threshold: float) -> float:
    return sum(1 for value in sorted_values if value > threshold) / len(sorted_values)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def natural_rep_key(path: Path) -> tuple[int, str]:
    name = path.name
    if name.startswith("rep"):
        try:
            return int(name[3:]), name
        except ValueError:
            pass
    return 10**9, name


def natural_r_key(path: Path) -> tuple[int, str]:
    name = path.parent.name if path.name.endswith(".out") else path.name
    if name.startswith("r"):
        try:
            return int(name[1:]), name
        except ValueError:
            pass
    return 10**9, name


def coerce_seed_map(value: Any, seeds: list[str]) -> dict[str, float]:
    if isinstance(value, dict):
        return {str(seed): float(v) for seed, v in value.items()}
    if isinstance(value, list):
        return {seed: float(value[i]) for i, seed in enumerate(seeds) if i < len(value)}
    return {}


def extract_w23a() -> Cohort:
    json_path = W23A_ROOT / "w23a_ddg.json"
    data = load_json(json_path)
    manifest_paths = sorted(
        (W23A_ROOT / "w23a/bound").glob("rep*/run_manifest.json"),
        key=lambda path: natural_rep_key(path.parent),
    )
    seeds: list[str] = []
    rep_indices: list[int] = []
    for manifest_path in manifest_paths:
        manifest = load_json(manifest_path)
        seeds.append(str(manifest["seed"]))
        rep_indices.append(int(manifest["replicate_index"]))
    per_ddg = data["primary_dgb"]["paired"]["per_rep_ddg"]
    ddg_by_seed = {
        seed: float(per_ddg[rep_idx])
        for seed, rep_idx in zip(seeds, rep_indices)
        if rep_idx < len(per_ddg)
    }
    dgb_by_leg_seed = {
        "bound": coerce_seed_map(data["bound"].get("per_seed_dgb_kcal"), seeds),
        "free": coerce_seed_map(data["free"].get("per_seed_dgb_kcal"), seeds),
    }
    return Cohort(
        name="w23a_gateA",
        label="W23A Gate-A carved",
        roots_by_leg={
            "bound": W23A_ROOT / "w23a/bound",
            "free": W23A_ROOT / "w23a/free",
        },
        seeds=seeds,
        ddg_by_seed=ddg_by_seed,
        dgb_by_leg_seed=dgb_by_leg_seed,
        highlight_seeds=["s101"],
        json_path=json_path,
    )


def extract_w4a(path: Path, name: str, label: str, prefix: str, highlights: list[str]) -> Cohort:
    data = load_json(path)
    pairing_key = "pairing_B_bound_minus_free_carved" if prefix == "carved" else "pairing_B_bound_minus_freeFIXAB"
    bound_key = "bound_carved" if prefix == "carved" else "bound_FIXAB"
    free_key = "free_carved" if prefix == "carved" else "free_FIXAB"
    pairing = data[pairing_key]["primary_dgb"]
    seeds = [str(seed) for seed in pairing["matched_seeds"]]
    ddg_by_seed = dict(zip(seeds, [float(value) for value in pairing["per_rep_ddg"]]))
    return Cohort(
        name=name,
        label=label,
        roots_by_leg={
            "bound": Path(data[bound_key]["out_root"]),
            "free": Path(data[free_key]["out_root"]),
        },
        seeds=seeds,
        ddg_by_seed=ddg_by_seed,
        dgb_by_leg_seed={
            "bound": coerce_seed_map(data[bound_key].get("seed_to_dgb"), seeds),
            "free": coerce_seed_map(data[free_key].get("seed_to_dgb"), seeds),
        },
        highlight_seeds=highlights,
        json_path=path,
    )


def build_cohorts() -> list[Cohort]:
    uncarved_json = W4A_C1_ROOT / "uncarved_reference_dgb.json"
    uncarved_data = load_json(uncarved_json)
    uncarved_pairing = uncarved_data["pairing_B_bound_minus_freeFIXAB"]["primary_dgb"]
    uncarved_seeds = [str(seed) for seed in uncarved_pairing["matched_seeds"]]
    uncarved_ddg = dict(zip(uncarved_seeds, [float(value) for value in uncarved_pairing["per_rep_ddg"]]))
    largest_uncarved = max(uncarved_ddg, key=lambda seed: abs(uncarved_ddg[seed]))
    return [
        extract_w23a(),
        extract_w4a(
            W4A_C1_ROOT / "carved_dgb.json",
            "w4a_c1_carved",
            "W4A C1 carved",
            "carved",
            ["s127"],
        ),
        extract_w4a(
            uncarved_json,
            "w4a_c1_uncarved",
            "W4A C1 uncarved FIXAB",
            "uncarved",
            sorted({"s127", largest_uncarved}),
        ),
    ]


def find_leg_root(base: Path, leg: str) -> Path:
    candidates = [base / "wt" / leg, base / leg, base]
    for candidate in candidates:
        if candidate.exists() and any(candidate.glob("rep*")):
            return candidate
    raise FileNotFoundError(f"could not locate replicate root for {base} leg={leg}")


def direction_out_files(rep_dir: Path, manifest: dict[str, Any], direction_tag: str) -> list[Path]:
    direction_meta = manifest.get("per_direction", {}).get(direction_tag, {})
    if direction_meta.get("subdir"):
        subdir = Path(direction_meta["subdir"])
    else:
        subdir = rep_dir / direction_tag
    if not subdir.exists():
        return []
    files = [
        path
        for path in subdir.glob("r*/trackb*.out")
        if direction_tag in path.name and path.is_file()
    ]
    return sorted(files, key=natural_r_key)


def parse_out_files(paths: list[Path]) -> tuple[SampleBag, dict[int, SampleBag], int]:
    all_values = SampleBag()
    by_state: dict[int, SampleBag] = defaultdict(SampleBag)
    malformed = 0
    for path in paths:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                cols = stripped.split()
                if len(cols) < 10:
                    malformed += 1
                    continue
                try:
                    state = int(float(cols[0]))
                    pert_e = float(cols[9])
                except ValueError:
                    malformed += 1
                    continue
                all_values.add(pert_e)
                by_state[state].add(pert_e)
    return all_values, by_state, malformed


def slim_mixing(manifest: dict[str, Any], direction_tag: str) -> dict[str, Any]:
    meta = manifest.get("per_direction", {}).get(direction_tag, {})
    mixing = meta.get("mixing", {}) or {}
    crossings = mixing.get("adjacent_crossings", {}) or {}
    min_pair = None
    min_crossing = None
    if crossings:
        min_pair, min_crossing = min(crossings.items(), key=lambda item: item[1])
    return {
        "n_states_declared": meta.get("n_states"),
        "n_cycles": meta.get("n_cycles"),
        "gate_passed": mixing.get("gate_passed"),
        "total_round_trips": mixing.get("total_round_trips"),
        "both_ends_visited_count": mixing.get("both_ends_visited_count"),
        "n_adjacent_pairs": mixing.get("n_adjacent_pairs"),
        "n_adjacent_with_crossings": mixing.get("n_adjacent_with_crossings"),
        "min_adjacent_pair": min_pair,
        "min_adjacent_crossing": min_crossing,
        "nan_any": meta.get("nan_any"),
        "apex_state": meta.get("apex_state"),
    }


def summarize_cell(
    cohort: Cohort,
    leg: str,
    rep_dir: Path,
    manifest: dict[str, Any],
    direction_tag: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    seed = str(manifest["seed"])
    out_files = direction_out_files(rep_dir, manifest, direction_tag)
    values, by_state, malformed = parse_out_files(out_files)
    global_summary = values.summary()

    state_rows: list[dict[str, Any]] = []
    for state, bag in sorted(by_state.items()):
        summary = bag.summary()
        state_rows.append(
            {
                "cohort": cohort.name,
                "seed": seed,
                "replicate_index": int(manifest["replicate_index"]),
                "leg": leg,
                "direction": direction_tag,
                "state": state,
                **summary,
            }
        )

    tail_state = None
    if state_rows:
        tail_state = max(state_rows, key=lambda row: none_low(row["p99"]))

    mixing = slim_mixing(manifest, direction_tag)
    cell = {
        "cohort": cohort.name,
        "label": cohort.label,
        "seed": seed,
        "highlight": seed in cohort.highlight_seeds,
        "replicate_index": int(manifest["replicate_index"]),
        "leg": leg,
        "direction": direction_tag,
        "out_files": len(out_files),
        "malformed_lines": malformed,
        "observed_states": sorted(by_state),
        "summary": global_summary,
        "tail_state_by_p99": tail_state,
        "mixing": mixing,
        "root": str(rep_dir),
    }
    return cell, state_rows


def none_low(value: Any) -> float:
    if value is None:
        return -math.inf
    return float(value)


def none_high(value: Any) -> float:
    if value is None:
        return math.inf
    return float(value)


def analyze_cohort(cohort: Cohort) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    cells: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    for leg, base in cohort.roots_by_leg.items():
        leg_root = find_leg_root(base, leg)
        for rep_dir in sorted(leg_root.glob("rep*"), key=natural_rep_key):
            manifest_path = rep_dir / "run_manifest.json"
            if not manifest_path.exists():
                raise FileNotFoundError(f"missing manifest: {manifest_path}")
            manifest = load_json(manifest_path)
            for direction_tag in ("dplus", "dminus"):
                cell, rows = summarize_cell(cohort, leg, rep_dir, manifest, direction_tag)
                cells.append(cell)
                state_rows.extend(rows)
    return cells, state_rows


def build_seed_summaries(cohort: Cohort, cells: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_seed: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for cell in cells:
        if cell["cohort"] == cohort.name:
            by_seed[cell["seed"]].append(cell)

    summaries: list[dict[str, Any]] = []
    for seed in cohort.seeds:
        seed_cells = by_seed.get(seed, [])
        max_p99_cell = max(seed_cells, key=lambda cell: none_low(cell["summary"]["p99"])) if seed_cells else None
        max_frac_cell = (
            max(seed_cells, key=lambda cell: none_low(cell["summary"]["frac_gt_150"]))
            if seed_cells
            else None
        )
        min_cross_cell = (
            min(seed_cells, key=lambda cell: none_high(cell["mixing"]["min_adjacent_crossing"]))
            if seed_cells
            else None
        )
        summaries.append(
            {
                "cohort": cohort.name,
                "seed": seed,
                "highlight": seed in cohort.highlight_seeds,
                "ddg": cohort.ddg_by_seed.get(seed),
                "abs_ddg": abs(cohort.ddg_by_seed[seed]) if seed in cohort.ddg_by_seed else None,
                "bound_dgb": cohort.dgb_by_leg_seed.get("bound", {}).get(seed),
                "free_dgb": cohort.dgb_by_leg_seed.get("free", {}).get(seed),
                "max_p99_cell": cell_ref(max_p99_cell),
                "max_pertE_p99": max_p99_cell["summary"]["p99"] if max_p99_cell else None,
                "max_frac_gt150_cell": cell_ref(max_frac_cell),
                "max_frac_gt150": max_frac_cell["summary"]["frac_gt_150"] if max_frac_cell else None,
                "min_adjacent_cell": cell_ref(min_cross_cell),
                "min_adjacent_crossing": min_cross_cell["mixing"]["min_adjacent_crossing"] if min_cross_cell else None,
                "any_gate_failed": any(cell["mixing"]["gate_passed"] is False for cell in seed_cells),
                "total_round_trips": sum(
                    int(cell["mixing"]["total_round_trips"] or 0) for cell in seed_cells
                ),
                "n_cells": len(seed_cells),
            }
        )
    return summaries


def cell_ref(cell: dict[str, Any] | None) -> str | None:
    if cell is None:
        return None
    tail = cell.get("tail_state_by_p99") or {}
    state = tail.get("state")
    state_part = f":state{state}" if state is not None else ""
    return f"{cell['leg']}/{cell['direction']}{state_part}"


def pearson(rows: list[dict[str, Any]], key_x: str, key_y: str) -> float | None:
    pairs = [
        (float(row[key_x]), float(row[key_y]))
        for row in rows
        if row.get(key_x) is not None and row.get(key_y) is not None
    ]
    if len(pairs) < 3:
        return None
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    sx = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    sy = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if sx == 0.0 or sy == 0.0:
        return None
    return sum((x - mean_x) * (y - mean_y) for x, y in pairs) / (sx * sy)


def median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def cohort_diagnostics(cohort: Cohort, seed_rows: list[dict[str, Any]]) -> dict[str, Any]:
    top_abs = max(seed_rows, key=lambda row: none_low(row["abs_ddg"]))
    top_p99 = max(seed_rows, key=lambda row: none_low(row["max_pertE_p99"]))
    top_frac = max(seed_rows, key=lambda row: none_low(row["max_frac_gt150"]))
    min_cross = min(seed_rows, key=lambda row: none_high(row["min_adjacent_crossing"]))
    controls = [row for row in seed_rows if row["seed"] not in cohort.highlight_seeds]
    highlights = [row for row in seed_rows if row["seed"] in cohort.highlight_seeds]
    control_medians = {
        "abs_ddg": median([float(row["abs_ddg"]) for row in controls if row["abs_ddg"] is not None]),
        "max_pertE_p99": median(
            [float(row["max_pertE_p99"]) for row in controls if row["max_pertE_p99"] is not None]
        ),
        "max_frac_gt150": median(
            [float(row["max_frac_gt150"]) for row in controls if row["max_frac_gt150"] is not None]
        ),
        "min_adjacent_crossing": median(
            [
                float(row["min_adjacent_crossing"])
                for row in controls
                if row["min_adjacent_crossing"] is not None
            ]
        ),
    }
    highlight_deltas = []
    for row in highlights:
        highlight_deltas.append(
            {
                "seed": row["seed"],
                "abs_ddg_minus_control_median": delta(row["abs_ddg"], control_medians["abs_ddg"]),
                "p99_minus_control_median": delta(row["max_pertE_p99"], control_medians["max_pertE_p99"]),
                "frac_gt150_minus_control_median": delta(
                    row["max_frac_gt150"], control_medians["max_frac_gt150"]
                ),
                "min_crossing_minus_control_median": delta(
                    row["min_adjacent_crossing"], control_medians["min_adjacent_crossing"]
                ),
            }
        )
    return {
        "top_abs_ddg_seed": top_abs["seed"],
        "top_abs_ddg": top_abs["abs_ddg"],
        "top_p99_seed": top_p99["seed"],
        "top_p99": top_p99["max_pertE_p99"],
        "top_frac_gt150_seed": top_frac["seed"],
        "top_frac_gt150": top_frac["max_frac_gt150"],
        "min_adjacent_seed": min_cross["seed"],
        "min_adjacent_crossing": min_cross["min_adjacent_crossing"],
        "highlight_seeds": cohort.highlight_seeds,
        "control_medians": control_medians,
        "highlight_minus_control_median": highlight_deltas,
        "pearson_absddg_vs_p99": pearson(seed_rows, "abs_ddg", "max_pertE_p99"),
        "pearson_absddg_vs_frac_gt150": pearson(seed_rows, "abs_ddg", "max_frac_gt150"),
        "pearson_absddg_vs_min_crossing": pearson(seed_rows, "abs_ddg", "min_adjacent_crossing"),
    }


def delta(value: Any, baseline: Any) -> float | None:
    if value is None or baseline is None:
        return None
    return float(value) - float(baseline)


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def flatten_cell_for_csv(cell: dict[str, Any]) -> dict[str, Any]:
    summary = cell["summary"]
    mixing = cell["mixing"]
    tail = cell.get("tail_state_by_p99") or {}
    return {
        "cohort": cell["cohort"],
        "seed": cell["seed"],
        "highlight": cell["highlight"],
        "replicate_index": cell["replicate_index"],
        "leg": cell["leg"],
        "direction": cell["direction"],
        "out_files": cell["out_files"],
        "malformed_lines": cell["malformed_lines"],
        "n": summary["n"],
        "mean": summary["mean"],
        "sd": summary["sd"],
        "p50": summary["p50"],
        "p95": summary["p95"],
        "p99": summary["p99"],
        "max": summary["max"],
        "frac_gt_120": summary["frac_gt_120"],
        "frac_gt_150": summary["frac_gt_150"],
        "frac_gt_180": summary["frac_gt_180"],
        "tail_state": tail.get("state"),
        "tail_state_p99": tail.get("p99"),
        "tail_state_frac_gt150": tail.get("frac_gt_150"),
        "n_states_declared": mixing["n_states_declared"],
        "n_cycles": mixing["n_cycles"],
        "gate_passed": mixing["gate_passed"],
        "total_round_trips": mixing["total_round_trips"],
        "min_adjacent_pair": mixing["min_adjacent_pair"],
        "min_adjacent_crossing": mixing["min_adjacent_crossing"],
        "nan_any": mixing["nan_any"],
        "root": cell["root"],
    }


def top_state_rows(state_rows: list[dict[str, Any]], limit: int = 40) -> list[dict[str, Any]]:
    return sorted(state_rows, key=lambda row: none_low(row["p99"]), reverse=True)[:limit]


def fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def md_table(headers: list[str], rows: list[list[Any]]) -> str:
    safe_headers = [str(header).replace("|", "\\|") for header in headers]
    lines = ["| " + " | ".join(safe_headers) + " |"]
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        safe_row = [fmt(value).replace("|", "\\|") for value in row]
        lines.append("| " + " | ".join(safe_row) + " |")
    return "\n".join(lines)


def write_report(
    cohorts: list[Cohort],
    seed_summaries: dict[str, list[dict[str, Any]]],
    diagnostics: dict[str, dict[str, Any]],
    state_rows: list[dict[str, Any]],
) -> None:
    report = OUT_DIR / "tail_probe_report.md"
    lines: list[str] = [
        "# Tail-Mechanism Probe Report",
        "",
        "Date: 2026-07-10",
        "",
        "Read-only analysis over existing Track B `.out`, manifest, and JSON artifacts.",
        "",
        "## Cohort Diagnostics",
        "",
    ]

    rows = []
    for cohort in cohorts:
        diag = diagnostics[cohort.name]
        rows.append(
            [
                cohort.name,
                ",".join(cohort.highlight_seeds),
                diag["top_abs_ddg_seed"],
                diag["top_abs_ddg"],
                diag["top_p99_seed"],
                diag["top_p99"],
                diag["top_frac_gt150_seed"],
                diag["top_frac_gt150"],
                diag["min_adjacent_seed"],
                diag["min_adjacent_crossing"],
            ]
        )
    lines.append(
        md_table(
            [
                "cohort",
                "highlight",
                "top |ddG| seed",
                "|ddG|",
                "top p99 seed",
                "p99",
                "top frac>150 seed",
                "frac>150",
                "min crossing seed",
                "min crossing",
            ],
            rows,
        )
    )

    for cohort in cohorts:
        lines.extend(["", f"## {cohort.name}", ""])
        rows = []
        for row in seed_summaries[cohort.name]:
            rows.append(
                [
                    row["seed"],
                    "*" if row["highlight"] else "",
                    row["ddg"],
                    row["bound_dgb"],
                    row["free_dgb"],
                    row["max_p99_cell"],
                    row["max_pertE_p99"],
                    row["max_frac_gt150_cell"],
                    row["max_frac_gt150"],
                    row["min_adjacent_cell"],
                    row["min_adjacent_crossing"],
                    row["total_round_trips"],
                    row["any_gate_failed"],
                ]
            )
        lines.append(
            md_table(
                [
                    "seed",
                    "flag",
                    "ddG",
                    "bound dgb",
                    "free dgb",
                    "max p99 cell",
                    "p99",
                    "max frac cell",
                    "frac>150",
                    "min crossing cell",
                    "min crossing",
                    "RT sum",
                    "gate fail",
                ],
                rows,
            )
        )
        diag = diagnostics[cohort.name]
        lines.extend(
            [
                "",
                "Correlation diagnostics:",
                "",
                md_table(
                    ["metric", "value"],
                    [
                        ["Pearson |ddG| vs max p99", diag["pearson_absddg_vs_p99"]],
                        ["Pearson |ddG| vs max frac>150", diag["pearson_absddg_vs_frac_gt150"]],
                        ["Pearson |ddG| vs min crossing", diag["pearson_absddg_vs_min_crossing"]],
                    ],
                ),
            ]
        )
        if diag["highlight_minus_control_median"]:
            lines.extend(["", "Highlight minus non-highlight median:", ""])
            rows = []
            for delta_row in diag["highlight_minus_control_median"]:
                rows.append(
                    [
                        delta_row["seed"],
                        delta_row["abs_ddg_minus_control_median"],
                        delta_row["p99_minus_control_median"],
                        delta_row["frac_gt150_minus_control_median"],
                        delta_row["min_crossing_minus_control_median"],
                    ]
                )
            lines.append(
                md_table(
                    ["seed", "abs ddG", "p99", "frac>150", "min crossing"],
                    rows,
                )
            )

    lines.extend(["", "## Top State-Localized pertE Tails", ""])
    rows = []
    for row in top_state_rows(state_rows, limit=30):
        rows.append(
            [
                row["cohort"],
                row["seed"],
                row["leg"],
                row["direction"],
                row["state"],
                row["p99"],
                row["frac_gt_150"],
                row["mean"],
                row["n"],
            ]
        )
    lines.append(
        md_table(
            ["cohort", "seed", "leg", "dir", "state", "p99", "frac>150", "mean", "n"],
            rows,
        )
    )

    report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    cohorts = build_cohorts()
    all_cells: list[dict[str, Any]] = []
    all_state_rows: list[dict[str, Any]] = []
    seed_summaries: dict[str, list[dict[str, Any]]] = {}
    diagnostics: dict[str, dict[str, Any]] = {}

    for cohort in cohorts:
        cells, state_rows = analyze_cohort(cohort)
        all_cells.extend(cells)
        all_state_rows.extend(state_rows)
        seed_rows = build_seed_summaries(cohort, cells)
        seed_summaries[cohort.name] = seed_rows
        diagnostics[cohort.name] = cohort_diagnostics(cohort, seed_rows)

    write_csv(
        OUT_DIR / "cell_summary.csv",
        [flatten_cell_for_csv(cell) for cell in all_cells],
        [
            "cohort",
            "seed",
            "highlight",
            "replicate_index",
            "leg",
            "direction",
            "out_files",
            "malformed_lines",
            "n",
            "mean",
            "sd",
            "p50",
            "p95",
            "p99",
            "max",
            "frac_gt_120",
            "frac_gt_150",
            "frac_gt_180",
            "tail_state",
            "tail_state_p99",
            "tail_state_frac_gt150",
            "n_states_declared",
            "n_cycles",
            "gate_passed",
            "total_round_trips",
            "min_adjacent_pair",
            "min_adjacent_crossing",
            "nan_any",
            "root",
        ],
    )

    flat_seed_rows = [row for rows in seed_summaries.values() for row in rows]
    write_csv(
        OUT_DIR / "seed_summary.csv",
        flat_seed_rows,
        [
            "cohort",
            "seed",
            "highlight",
            "ddg",
            "abs_ddg",
            "bound_dgb",
            "free_dgb",
            "max_p99_cell",
            "max_pertE_p99",
            "max_frac_gt150_cell",
            "max_frac_gt150",
            "min_adjacent_cell",
            "min_adjacent_crossing",
            "any_gate_failed",
            "total_round_trips",
            "n_cells",
        ],
    )

    write_csv(
        OUT_DIR / "state_tail_summary.csv",
        all_state_rows,
        [
            "cohort",
            "seed",
            "replicate_index",
            "leg",
            "direction",
            "state",
            "n",
            "mean",
            "sd",
            "min",
            "p01",
            "p05",
            "p50",
            "p95",
            "p99",
            "max",
            "frac_gt_120",
            "frac_gt_150",
            "frac_gt_180",
            "frac_lt_0",
        ],
    )

    summary = {
        "analysis": "tail_mechanism_probe_20260710",
        "repo_root": str(REPO_ROOT),
        "inputs": {
            cohort.name: {
                "json_path": str(cohort.json_path),
                "roots_by_leg": {leg: str(path) for leg, path in cohort.roots_by_leg.items()},
            }
            for cohort in cohorts
        },
        "cohorts": {
            cohort.name: {
                "label": cohort.label,
                "seeds": cohort.seeds,
                "highlight_seeds": cohort.highlight_seeds,
                "ddg_by_seed": cohort.ddg_by_seed,
                "seed_summaries": seed_summaries[cohort.name],
                "diagnostics": diagnostics[cohort.name],
            }
            for cohort in cohorts
        },
        "n_cells": len(all_cells),
        "n_state_rows": len(all_state_rows),
        "top_state_rows_by_p99": top_state_rows(all_state_rows, limit=40),
    }
    (OUT_DIR / "tail_probe_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(cohorts, seed_summaries, diagnostics, all_state_rows)


if __name__ == "__main__":
    main()
