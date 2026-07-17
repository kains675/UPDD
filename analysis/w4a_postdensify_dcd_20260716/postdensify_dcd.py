#!/usr/bin/env python3
"""Run the preregistered W4A post-densify bound/dminus DCD sidecar."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path
from statistics import median
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utils.control_center.control_token import PAUSE_EXIT_CODE, boundary_pause_requested

ANALYSIS_DIR = Path(__file__).resolve().parent
PREREGISTRATION = ANALYSIS_DIR / "PREREGISTRATION.md"
SOURCE_INVENTORY = ANALYSIS_DIR / "source_inventory.json"
SUMMARY_JSON = ANALYSIS_DIR / "postdensify_dcd_summary.json"
SUMMARY_MD = ANALYSIS_DIR / "postdensify_dcd_report.md"

SOURCE_ROOT = REPO_ROOT / "outputs/_trackb/w4a_bdcd_h18_20260710_rerun1"
OUT_ROOT = REPO_ROOT / "outputs/_trackb/w4a_postdensify_dcd_20260716"
D2_PATH = ANALYSIS_DIR.parent / "d2_tail_dcd_probe_20260710/d2_tail_dcd_probe.py"

SEEDS = ("s7", "s19", "s23", "s42", "s83", "s101", "s127", "s163", "s199", "s251")
ARMS = {
    "carved": {"source_dir": "carved_bound", "carve_void_waters": True},
    "uncarved": {"source_dir": "uncarved_bound", "carve_void_waters": False},
}
HIGH_LEVERAGE_SEED = "s127"
LAMBDA1_RAMPDOWN = (0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
LAMBDA2_RAMPUP = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
EXPECTED_DMINUS_L1 = tuple(reversed(LAMBDA1_RAMPDOWN)) + (0.0,) * 8
EXPECTED_DMINUS_L2 = (0.5,) * 8 + tuple(reversed(LAMBDA2_RAMPUP[:-1]))
EXPECTED_SOURCE_WALKERS = 15
EXPECTED_SOURCE_ROWS = 800
EXPECTED_SOURCE_COLUMNS = 11
EXPECTED_DCD_FRAMES = 15 * 50


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_d2():
    d2 = load_module("postdensify_d2_core", D2_PATH)
    d2.ANALYSIS_DIR = ANALYSIS_DIR
    d2.OUT_ROOT = OUT_ROOT
    return d2


def relative(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_record(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {
        "path": relative(path),
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def parse_cntl_vector(path: Path, key: str) -> tuple[float, ...]:
    prefix = key + " ="
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith(prefix):
            raw = line.split("=", 1)[1].strip()
            value = ast.literal_eval(raw)
            return tuple(float(item.strip()) for item in value.split(","))
    raise ValueError(f"{path}: missing {key}")


def parse_numeric_out(path: Path, expected_rows: int) -> dict[str, Any]:
    rows = 0
    states: set[int] = set()
    for line_number, line in enumerate(path.open("r", encoding="utf-8", errors="replace"), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) != EXPECTED_SOURCE_COLUMNS:
            raise ValueError(
                f"{path}:{line_number}: expected {EXPECTED_SOURCE_COLUMNS} columns, got {len(fields)}"
            )
        values = [float(field) for field in fields]
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"{path}:{line_number}: non-finite value")
        state = int(values[0])
        if values[0] != state or not 0 <= state < EXPECTED_SOURCE_WALKERS:
            raise ValueError(f"{path}:{line_number}: invalid state {values[0]}")
        rows += 1
        states.add(state)
    if rows != expected_rows:
        raise ValueError(f"{path}: expected {expected_rows} rows, got {rows}")
    return {"rows": rows, "states": sorted(states)}


def inspect_pdb(path: Path) -> dict[str, int]:
    n_atoms = 0
    n_site_heavy = 0
    n_waters = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            n_atoms += 1
            atom_name = line[12:16].strip().upper()
            resname = line[17:20].strip().upper()
            chain = line[21].strip()
            resid = line[22:26].strip()
            element = line[76:78].strip().upper()
            if resname in {"HOH", "WAT", "SOL"} and atom_name in {"O", "OW", "OH2"}:
                n_waters += 1
            if chain == "B" and resid == "4" and resname in {"TRP", "ALA"}:
                is_hydrogen = element == "H" or (not element and atom_name.startswith("H"))
                if not is_hydrogen:
                    n_site_heavy += 1
    if n_site_heavy == 0:
        raise ValueError(f"{path}: no chain B residue 4 site-heavy atoms")
    return {"n_atoms": n_atoms, "n_waters": n_waters, "n_site_heavy": n_site_heavy}


def make_cells(d2) -> tuple[Any, ...]:
    cells = []
    for arm, arm_config in ARMS.items():
        for replicate_index, seed in enumerate(SEEDS):
            cells.append(
                d2.Cell(
                    cell_id=f"w4a_{arm}_{seed}_bound_dminus",
                    system=f"w4a_{arm}",
                    source_rep_dir=(
                        f"outputs/_trackb/w4a_bdcd_h18_20260710_rerun1/"
                        f"{arm_config['source_dir']}/wt/bound/rep{replicate_index}"
                    ),
                    endpoint="wt",
                    leg="bound",
                    replicate_index=replicate_index,
                    seed=seed,
                    direction="dminus",
                    site_chain="B",
                    site_resid="4",
                    lambda1_rampdown=LAMBDA1_RAMPDOWN,
                    lambda2_rampup=LAMBDA2_RAMPUP,
                    control=arm == "carved" and seed != HIGH_LEVERAGE_SEED,
                )
            )
    return tuple(cells)


def validate_arm_preregistration(arm: str) -> dict[str, Any]:
    arm_config = ARMS[arm]
    path = SOURCE_ROOT / arm_config["source_dir"] / "pre_registration.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    config = payload.get("config") or {}
    expected = {
        "leg": "bound",
        "endpoints": ["wt"],
        "seeds": list(SEEDS),
        "directions": ["dplus", "dminus"],
        "n_cycles": EXPECTED_SOURCE_ROWS,
        "md_steps_per_cycle": 250,
        "construction": "twocopy",
        "mutation_spec": "w4a_trp_ala_res4",
        "lambda1_rampdown": list(LAMBDA1_RAMPDOWN),
        "lambda2_rampup": list(LAMBDA2_RAMPUP),
        "carve_void_waters": arm_config["carve_void_waters"],
        "mintimeid": 600,
    }
    for key, value in expected.items():
        if config.get(key) != value:
            raise ValueError(f"{path}: {key}={config.get(key)!r}, expected {value!r}")
    return {"config": expected, "artifact": artifact_record(path)}


def validate_source_cell(cell: Any) -> dict[str, Any]:
    rep_dir = REPO_ROOT / cell.source_rep_dir
    manifest_path = rep_dir / "run_manifest.json"
    xml_path = rep_dir / "inplace_rbfe_bound_sys.xml"
    pdb_path = rep_dir / "inplace_rbfe_bound.pdb"
    cntl_path = rep_dir / "dminus/trackb_dminus_asyncre.cntl"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    direction = (manifest.get("per_direction") or {}).get("dminus") or {}

    expected_manifest = {
        "seed": cell.seed,
        "replicate_index": cell.replicate_index,
        "endpoint": "wt",
        "leg": "bound",
        "construction": "twocopy",
    }
    for key, value in expected_manifest.items():
        if manifest.get(key) != value:
            raise ValueError(f"{manifest_path}: {key}={manifest.get(key)!r}, expected {value!r}")
    for key, value in {"n_states": 15, "n_cycles": 800, "md_steps_per_cycle": 250, "nan_any": False}.items():
        if direction.get(key) != value:
            raise ValueError(f"{manifest_path}: dminus.{key}={direction.get(key)!r}, expected {value!r}")

    lambda1 = parse_cntl_vector(cntl_path, "LAMBDA1")
    lambda2 = parse_cntl_vector(cntl_path, "LAMBDA2")
    if lambda1 != EXPECTED_DMINUS_L1 or lambda2 != EXPECTED_DMINUS_L2:
        raise ValueError(f"{cntl_path}: H18 dminus schedule mismatch")

    observed_states: set[int] = set()
    walker_rows: dict[str, int] = {}
    for walker_index in range(EXPECTED_SOURCE_WALKERS):
        walker = f"r{walker_index}"
        out_path = rep_dir / "dminus" / walker / "trackb_dminus.out"
        parsed = parse_numeric_out(out_path, EXPECTED_SOURCE_ROWS)
        walker_rows[walker] = parsed["rows"]
        observed_states.update(parsed["states"])
    if observed_states != set(range(EXPECTED_SOURCE_WALKERS)):
        raise ValueError(f"{rep_dir}: source state coverage is {sorted(observed_states)}")

    return {
        "cell_id": cell.cell_id,
        "arm": cell.system.removeprefix("w4a_"),
        "seed": cell.seed,
        "replicate_index": cell.replicate_index,
        "source_rep_dir": cell.source_rep_dir,
        "artifacts": {
            "run_manifest": artifact_record(manifest_path),
            "serialized_xml": artifact_record(xml_path),
            "serialized_pdb": artifact_record(pdb_path),
            "dminus_cntl": artifact_record(cntl_path),
        },
        "pdb": inspect_pdb(pdb_path),
        "source_direction": {
            "n_walkers": EXPECTED_SOURCE_WALKERS,
            "rows_per_walker": walker_rows,
            "states": sorted(observed_states),
            "n_columns": EXPECTED_SOURCE_COLUMNS,
            "finite": True,
            "nan_any_manifest": direction["nan_any"],
        },
    }


def inventory_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("inventory_digest", None)
    encoded = json.dumps(canonical, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_source_inventory(d2) -> dict[str, Any]:
    cells = make_cells(d2)
    payload = {
        "schema": "w4a_postdensify_dcd_source_inventory_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "source_root": relative(SOURCE_ROOT),
        "output_root": relative(OUT_ROOT),
        "analysis_preregistration": artifact_record(PREREGISTRATION),
        "arms": {arm: validate_arm_preregistration(arm) for arm in ARMS},
        "cells": [validate_source_cell(cell) for cell in cells],
    }

    by_key = {(cell["arm"], cell["seed"]): cell for cell in payload["cells"]}
    for seed in SEEDS:
        carved = by_key[("carved", seed)]["pdb"]
        uncarved = by_key[("uncarved", seed)]["pdb"]
        atom_delta = uncarved["n_atoms"] - carved["n_atoms"]
        water_delta = uncarved["n_waters"] - carved["n_waters"]
        if atom_delta <= 0 or atom_delta != 3 * water_delta:
            raise ValueError(
                f"{seed}: carve declaration/artifact mismatch: atom delta={atom_delta}, water delta={water_delta}"
            )

    payload["inventory_digest"] = inventory_digest(payload)
    return payload


def verify_source_inventory(d2) -> dict[str, Any]:
    frozen = load_frozen_inventory()
    current = build_source_inventory(d2)
    if frozen["inventory_digest"] != current["inventory_digest"]:
        raise ValueError(
            "source/preregistration drift after inventory freeze: "
            f"stored={frozen['inventory_digest']} current={current['inventory_digest']}"
        )
    return frozen


def load_frozen_inventory(expected_digest: str | None = None) -> dict[str, Any]:
    if not SOURCE_INVENTORY.is_file():
        raise FileNotFoundError(
            f"missing {SOURCE_INVENTORY}; run --prepare-inventory before launch"
        )
    frozen = json.loads(SOURCE_INVENTORY.read_text(encoding="utf-8"))
    if frozen.get("inventory_digest") != inventory_digest(frozen):
        raise ValueError(f"{SOURCE_INVENTORY}: stored inventory digest mismatch")
    if expected_digest is not None and frozen["inventory_digest"] != expected_digest:
        raise ValueError(
            "worker inventory digest mismatch: "
            f"expected={expected_digest} stored={frozen['inventory_digest']}"
        )
    return frozen


def verify_cell_source_inventory(cell: Any, inventory: dict[str, Any]) -> None:
    arm = cell.system.removeprefix("w4a_")
    frozen_arm = (inventory.get("arms") or {}).get(arm)
    if frozen_arm is None:
        raise ValueError(f"{cell.cell_id}: arm {arm!r} absent from frozen inventory")
    current_arm = validate_arm_preregistration(arm)
    if current_arm != frozen_arm:
        raise ValueError(f"{cell.cell_id}: arm preregistration drift after inventory freeze")

    matching = [row for row in inventory.get("cells", []) if row.get("cell_id") == cell.cell_id]
    if len(matching) != 1:
        raise ValueError(f"{cell.cell_id}: expected one frozen cell record, got {len(matching)}")
    current = validate_source_cell(cell)
    if current != matching[0]:
        raise ValueError(f"{cell.cell_id}: source artifacts drifted after inventory freeze")


def cell_rep_dir(d2, cell: Any) -> Path:
    return d2.cell_rep_dir(cell)


def manifest_path(d2, cell: Any) -> Path:
    return cell_rep_dir(d2, cell) / "d2_cell_manifest.json"


def validate_completed_manifest(d2, cell: Any, manifest: dict[str, Any], inventory: dict[str, Any]) -> None:
    cell_data = manifest.get("cell") or {}
    for key, value in {
        "cell_id": cell.cell_id,
        "system": cell.system,
        "seed": cell.seed,
        "replicate_index": cell.replicate_index,
        "leg": "bound",
        "direction": "dminus",
    }.items():
        if cell_data.get(key) != value:
            raise ValueError(f"{cell.cell_id}: completed manifest {key} mismatch")
    context = manifest.get("postdensify_context") or {}
    if context.get("source_inventory_digest") != inventory["inventory_digest"]:
        raise ValueError(f"{cell.cell_id}: completed manifest inventory digest mismatch")
    run = manifest.get("run") or {}
    if run.get("n_cycles") != 50 or run.get("md_steps_per_cycle") != 250 or run.get("nan_any") is not False:
        raise ValueError(f"{cell.cell_id}: completed run protocol or NaN mismatch")
    schedule = manifest.get("schedule") or {}
    if schedule.get("n_states") != 15:
        raise ValueError(f"{cell.cell_id}: completed schedule is not 15-state H18")
    analysis = manifest.get("dcd_analysis") or {}
    if analysis.get("n_frames") != EXPECTED_DCD_FRAMES:
        raise ValueError(f"{cell.cell_id}: expected {EXPECTED_DCD_FRAMES} DCD frames")
    walkers = analysis.get("by_walker") or []
    if len(walkers) != 15 or any(row.get("n_frames") != 50 for row in walkers):
        raise ValueError(f"{cell.cell_id}: DCD walker/frame mismatch")
    selection = manifest.get("selection") or {}
    if selection.get("n_site_heavy_atoms", 0) <= 0 or selection.get("n_local_water_oxygens", 0) <= 0:
        raise ValueError(f"{cell.cell_id}: empty site-heavy or local-water selection")

    direction_dir = cell_rep_dir(d2, cell) / "dminus"
    dcd_paths = sorted((direction_dir / "dcd").glob("r*/trackb_dminus.dcd"))
    if len(dcd_paths) != 15 or any(path.stat().st_size == 0 for path in dcd_paths):
        raise ValueError(f"{cell.cell_id}: DCD file count/size mismatch")
    for walker_index in range(15):
        parse_numeric_out(direction_dir / f"r{walker_index}/trackb_dminus.out", 50)


def load_completed(d2, cell: Any, inventory: dict[str, Any]) -> dict[str, Any] | None:
    path = manifest_path(d2, cell)
    if not path.is_file():
        return None
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    if not isinstance(manifest, dict) or "postdensify_context" not in manifest:
        return None
    validate_completed_manifest(d2, cell, manifest, inventory)
    return manifest


def archive_partial_cell(cell: Any) -> Path | None:
    cell_root = OUT_ROOT / cell.cell_id
    if not cell_root.exists():
        return None
    stamp = time.strftime("%Y%m%d_%H%M%S")
    archive_root = OUT_ROOT / "_archive" / f"{stamp}_{cell.cell_id}"
    archive_root.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(cell_root), str(archive_root))
    atomic_write_json(
        archive_root / "archive_reason.json",
        {
            "cell_id": cell.cell_id,
            "archived": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "reason": "interrupted sidecar cell without a completed manifest",
        },
    )
    return archive_root


def add_context(manifest: dict[str, Any], inventory: dict[str, Any]) -> dict[str, Any]:
    manifest["postdensify_context"] = {
        "schema": "w4a_postdensify_dcd_cell_v1",
        "source_inventory_digest": inventory["inventory_digest"],
        "analysis_preregistration_sha256": inventory["analysis_preregistration"]["sha256"],
        "source_inventory_path": relative(SOURCE_INVENTORY),
    }
    return manifest


def summary_row(manifest: dict[str, Any]) -> dict[str, Any]:
    cell = manifest["cell"]
    analysis = manifest["dcd_analysis"]
    selection = manifest["selection"]
    return {
        "cell_id": cell["cell_id"],
        "arm": cell["system"].removeprefix("w4a_"),
        "seed": cell["seed"],
        "replicate_index": cell["replicate_index"],
        "high_leverage": cell["system"] == "w4a_carved" and cell["seed"] == HIGH_LEVERAGE_SEED,
        "control": bool(cell["control"]),
        "leg": cell["leg"],
        "direction": cell["direction"],
        "n_local_water_residues": selection["n_local_water_residues"],
        "n_local_water_oxygens": selection["n_local_water_oxygens"],
        "n_selected_atoms": selection["n_selected_atoms"],
        "n_frames": analysis["n_frames"],
        "water_min_nm": analysis["water_min_site_heavy_nm"],
        "site_rmsd_nm": analysis["site_heavy_rmsd_nm"],
        "mixing": manifest["mixing"],
        "nan_any": manifest["run"]["nan_any"],
        "elapsed_s": manifest["run"]["elapsed_s"],
    }


def fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return "NA"
    return f"{float(value):.{digits}f}"


def write_summary(manifests: list[dict[str, Any]], cells: tuple[Any, ...], inventory: dict[str, Any]) -> dict[str, Any]:
    rows_by_id = {manifest["cell"]["cell_id"]: summary_row(manifest) for manifest in manifests}
    rows = [rows_by_id[cell.cell_id] for cell in cells if cell.cell_id in rows_by_id]
    pending = [cell.cell_id for cell in cells if cell.cell_id not in rows_by_id]
    elapsed = [float(row["elapsed_s"]) for row in rows]
    eta_s = median(elapsed) * len(pending) if elapsed else None
    payload = {
        "schema": "w4a_postdensify_dcd_summary_v1",
        "analysis": "w4a_postdensify_dcd_20260716",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "COMPLETE" if not pending else "RUNNING",
        "source_inventory_digest": inventory["inventory_digest"],
        "out_root": relative(OUT_ROOT),
        "n_expected": len(cells),
        "n_completed": len(rows),
        "pending_cell_ids": pending,
        "median_elapsed_s": median(elapsed) if elapsed else None,
        "eta_s": eta_s,
        "cells": rows,
    }
    atomic_write_json(SUMMARY_JSON, payload)
    atomic_write_json(OUT_ROOT / "progress.json", payload)

    lines = [
        "# W4A Post-densify DCD Progress",
        "",
        f"Status: `{payload['status']}`",
        f"Completed: `{payload['n_completed']}/{payload['n_expected']}`",
        f"Estimated remaining seconds: `{fmt(payload['eta_s'], 0)}`",
        "",
        "| arm | seed | high | control | waters | frames | frac<0.26 | frac<0.35 | RMSD p95 nm | elapsed min |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        water = row["water_min_nm"]
        lines.append(
            "| {arm} | {seed} | {high} | {control} | {waters} | {frames} | {f026} | {f035} | {rmsd} | {elapsed} |".format(
                arm=row["arm"],
                seed=row["seed"],
                high=row["high_leverage"],
                control=row["control"],
                waters=row["n_local_water_residues"],
                frames=row["n_frames"],
                f026=fmt(water["frac_lt_0p26"]),
                f035=fmt(water["frac_lt_0p35"]),
                rmsd=fmt(row["site_rmsd_nm"]["p95"]),
                elapsed=fmt(row["elapsed_s"] / 60.0, 1),
            )
        )
    if pending:
        lines.extend(["", "Pending: " + ", ".join(f"`{cell_id}`" for cell_id in pending)])
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return payload


def write_run_state(status: str, **extra: Any) -> None:
    payload = {"status": status, "updated": time.strftime("%Y-%m-%d %H:%M:%S %z"), **extra}
    atomic_write_json(OUT_ROOT / "run_state.json", payload)


def cell_by_id(cells: tuple[Any, ...], cell_id: str) -> Any:
    matching = [cell for cell in cells if cell.cell_id == cell_id]
    if len(matching) != 1:
        raise ValueError(f"expected one cell named {cell_id!r}, got {len(matching)}")
    return matching[0]


def cell_worker_command(cell: Any, platform: str, inventory_digest_value: str) -> list[str]:
    return [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker-cell",
        cell.cell_id,
        "--expected-inventory-digest",
        inventory_digest_value,
        "--platform",
        platform,
    ]


def launch_cell_worker(cell: Any, platform: str, inventory_digest_value: str) -> None:
    command = cell_worker_command(cell, platform, inventory_digest_value)
    result = subprocess.run(command, cwd=REPO_ROOT, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"{cell.cell_id}: worker exited with code {result.returncode}")


def run_cell_worker(
    d2,
    cells: tuple[Any, ...],
    inventory: dict[str, Any],
    cell_id: str,
    platform: str,
) -> None:
    cell = cell_by_id(cells, cell_id)
    completed = load_completed(d2, cell, inventory)
    if completed is not None:
        print(f"[POST-DCD-WORKER] already complete {cell.cell_id}", flush=True)
        return
    if (OUT_ROOT / cell.cell_id).exists():
        raise RuntimeError(
            f"{cell.cell_id}: partial output exists; parent must archive it before worker launch"
        )

    verify_cell_source_inventory(cell, inventory)
    prod = d2.load_prod()
    rbfe = prod._load_rbfe()
    print(f"[POST-DCD-WORKER] running {cell.cell_id} (seed={cell.seed})", flush=True)
    manifest = d2.run_cell(prod, rbfe, cell, platform=platform)
    manifest = add_context(manifest, inventory)
    atomic_write_json(manifest_path(d2, cell), manifest)
    validate_completed_manifest(d2, cell, manifest, inventory)
    print(f"[POST-DCD-WORKER] complete {cell.cell_id}", flush=True)


def _pause_at_cell_boundary(manifests: list[dict[str, Any]], cells: tuple[Any, ...]) -> bool:
    progress = {
        "boundary": "dcd_cell",
        "completed": len(manifests),
        "expected": len(cells),
    }
    if not boundary_pause_requested(progress):
        return False
    write_run_state(
        "PAUSED",
        n_completed=len(manifests),
        n_expected=len(cells),
        boundary="dcd_cell",
    )
    print("[POST-DCD] cooperative pause acknowledged at cell boundary", flush=True)
    return True


def run_all(d2, cells: tuple[Any, ...], inventory: dict[str, Any], platform: str) -> bool:
    manifests: list[dict[str, Any]] = []

    for cell in cells:
        completed = load_completed(d2, cell, inventory)
        if completed is not None:
            manifests.append(completed)
    write_summary(manifests, cells, inventory)
    if _pause_at_cell_boundary(manifests, cells):
        return True

    for cell in cells:
        if any(manifest["cell"]["cell_id"] == cell.cell_id for manifest in manifests):
            print(f"[POST-DCD] resume-skip {cell.cell_id}", flush=True)
            continue
        archive = archive_partial_cell(cell)
        if archive is not None:
            print(f"[POST-DCD] archived partial cell at {archive}", flush=True)
        write_run_state("RUNNING", current_cell=cell.cell_id, n_completed=len(manifests), n_expected=len(cells))
        print(f"[POST-DCD] running {cell.cell_id} (seed={cell.seed})", flush=True)
        try:
            launch_cell_worker(cell, platform, inventory["inventory_digest"])
            manifest = load_completed(d2, cell, inventory)
            if manifest is None:
                raise RuntimeError(f"{cell.cell_id}: worker exited successfully without a completed manifest")
        except Exception as exc:
            write_run_state(
                "FAILED",
                current_cell=cell.cell_id,
                n_completed=len(manifests),
                n_expected=len(cells),
                error=f"{type(exc).__name__}: {exc}",
            )
            raise
        manifests.append(manifest)
        write_summary(manifests, cells, inventory)
        if _pause_at_cell_boundary(manifests, cells):
            return True

    write_run_state("COMPLETE", n_completed=len(manifests), n_expected=len(cells))
    return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-inventory", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--platform", default="CUDA")
    parser.add_argument("--worker-cell", help=argparse.SUPPRESS)
    parser.add_argument("--expected-inventory-digest", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    d2 = load_d2()
    cells = make_cells(d2)
    if args.worker_cell:
        if args.prepare_inventory or args.dry_run or args.analyze_only:
            parser.error("--worker-cell cannot be combined with a public execution mode")
        if not args.expected_inventory_digest:
            parser.error("--worker-cell requires --expected-inventory-digest")
        inventory = load_frozen_inventory(args.expected_inventory_digest)
        run_cell_worker(d2, cells, inventory, args.worker_cell, args.platform)
        return 0
    if args.expected_inventory_digest:
        parser.error("--expected-inventory-digest is valid only with --worker-cell")

    if args.prepare_inventory:
        inventory = build_source_inventory(d2)
        atomic_write_json(SOURCE_INVENTORY, inventory)
        print(json.dumps({"n_cells": len(inventory["cells"]), "inventory_digest": inventory["inventory_digest"]}))
        return 0

    inventory = verify_source_inventory(d2)
    if args.dry_run:
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "n_cells": len(cells),
                    "inventory_digest": inventory["inventory_digest"],
                    "cells": [cell.cell_id for cell in cells],
                },
                indent=2,
            )
        )
        return 0

    if args.analyze_only:
        manifests = []
        for cell in cells:
            completed = load_completed(d2, cell, inventory)
            if completed is not None:
                manifests.append(completed)
        write_summary(manifests, cells, inventory)
        return 0

    paused = run_all(d2, cells, inventory, args.platform)
    return PAUSE_EXIT_CODE if paused else 0


if __name__ == "__main__":
    raise SystemExit(main())
