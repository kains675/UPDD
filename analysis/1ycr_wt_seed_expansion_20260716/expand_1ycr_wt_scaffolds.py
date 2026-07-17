#!/usr/bin/env python3
"""Prepare and run the preregistered five-seed 1YCR WT scaffold expansion."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from statistics import median
from typing import Any

import mdtraj as md
import numpy as np
import openmm
from openmm.app import PDBFile


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utils.control_center.control_token import PAUSE_EXIT_CODE, boundary_pause_requested

ANALYSIS_DIR = Path(__file__).resolve().parent
PREREGISTRATION = ANALYSIS_DIR / "PREREGISTRATION.md"
SOURCE_INVENTORY = ANALYSIS_DIR / "source_inventory.json"
SUMMARY_JSON = ANALYSIS_DIR / "scaffold_expansion_summary.json"
SUMMARY_MD = ANALYSIS_DIR / "scaffold_expansion_report.md"
RUNTIME_ROOT = REPO_ROOT / "outputs/_trackb/1ycr_wt_seed_expansion_20260716"

QMMM_PYTHON = "/home/san/miniconda3/envs/qmmm/bin/python"
MD_SCRIPT = REPO_ROOT / "utils/run_restrained_md.py"
SOURCE_INPUT = REPO_ROOT / "outputs/1YCR_WT_calib_s7/_md_input/1YCR_WT.pdb"
EXPECTED_MD_SCRIPT_SHA256 = "b27bc5794cf7bcd61682a67586077024c3d3e15b9a5acc4763718c4a22b5ce10"
EXPECTED_INPUT_SHA256 = "20f87791637e71a189e5e641f3ec5b807a845aca60bca6f1942ac879cda00cf7"
HISTORICAL_COMMIT = "3d2a10c"
HISTORICAL_MD_BLOB = "af0c673571d906f0e616730646da5b0f8830e26c"
CURRENT_MD_BLOB = "fd3249e9c1feff9571097f140561d70a438ce202"

EXISTING_SEEDS = (7, 19, 23, 42, 101)
NEW_SEEDS = (83, 127, 163, 199, 251)
N_STEPS = 2_500_000
EXPECTED_REPORT_ROWS = 100
EXPECTED_FINAL_STEP = 2_525_000
EXPECTED_DCD_FRAMES = 500
DCD_ARCHIVE_ROOT = (
    Path("/media/san/ExpDATA/UPDD_proj_Backup/outputs_dcd_archive_20260623/outputs")
)


def relative(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {"path": relative(path), "size": path.stat().st_size, "sha256": sha256_file(path)}


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def run_checked(command: list[str]) -> str:
    result = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
    )
    return result.stdout.strip()


def parse_md_log(path: Path) -> dict[str, Any]:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) != 7:
                raise ValueError(f"{path}:{line_number}: expected 7 tab-separated fields")
            try:
                progress = float(fields[0].rstrip("%"))
                values = [progress] + [float(field) for field in fields[1:6]]
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: non-numeric reporter row") from exc
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"{path}:{line_number}: non-finite reporter value")
            rows.append(values)

    if len(rows) != EXPECTED_REPORT_ROWS:
        raise ValueError(f"{path}: expected {EXPECTED_REPORT_ROWS} reporter rows, got {len(rows)}")
    if int(rows[-1][1]) != EXPECTED_FINAL_STEP:
        raise ValueError(f"{path}: final step {rows[-1][1]} != {EXPECTED_FINAL_STEP}")
    temperatures = [row[4] for row in rows]
    volumes = [row[5] for row in rows]
    if not all(250.0 <= value <= 350.0 for value in temperatures):
        raise ValueError(f"{path}: temperature outside 250..350 K")
    if not all(value > 0.0 for value in volumes):
        raise ValueError(f"{path}: non-positive box volume")
    return {
        "n_rows": len(rows),
        "first_step": int(rows[0][1]),
        "final_step": int(rows[-1][1]),
        "temperature_min_K": min(temperatures),
        "temperature_max_K": max(temperatures),
        "volume_min_nm3": min(volumes),
        "volume_max_nm3": max(volumes),
    }


def inspect_pdb(path: Path) -> dict[str, Any]:
    pdb = PDBFile(str(path))
    positions = pdb.positions.value_in_unit(openmm.unit.nanometer)
    xyz = np.asarray([[value.x, value.y, value.z] for value in positions], dtype=float)
    if xyz.size == 0 or not np.all(np.isfinite(xyz)):
        raise ValueError(f"{path}: empty or non-finite PDB coordinates")
    chains = sorted({chain.id for chain in pdb.topology.chains()})
    if "A" not in chains or "B" not in chains:
        raise ValueError(f"{path}: expected chains A/B, got {chains}")
    return {
        "n_atoms": int(xyz.shape[0]),
        "n_residues": sum(1 for _ in pdb.topology.residues()),
        "chains": chains,
        "coordinate_min_nm": float(np.min(xyz)),
        "coordinate_max_nm": float(np.max(xyz)),
    }


def inspect_dcd(path: Path) -> dict[str, int]:
    with md.formats.DCDTrajectoryFile(str(path), mode="r") as handle:
        n_frames = len(handle)
        xyz, _, _ = handle.read(n_frames=1)
    if n_frames != EXPECTED_DCD_FRAMES:
        raise ValueError(f"{path}: expected {EXPECTED_DCD_FRAMES} DCD frames, got {n_frames}")
    return {"n_frames": n_frames, "n_atoms": int(xyz.shape[1])}


def validate_outputs(output_root: Path, stdout_log: Path) -> dict[str, Any]:
    mdresult = output_root / "mdresult"
    final_pdb = mdresult / "1YCR_WT_final.pdb"
    md_log = mdresult / "1YCR_WT_md.log"
    dcd = mdresult / "1YCR_WT_restrained.dcd"
    for path in (final_pdb, md_log, dcd, stdout_log):
        if not path.is_file():
            raise FileNotFoundError(path)

    stdout = stdout_log.read_text(encoding="utf-8", errors="replace")
    complete_events = []
    for line in stdout.splitlines():
        if not line.startswith("[DIAG] "):
            continue
        try:
            event = json.loads(line[len("[DIAG] ") :])
        except json.JSONDecodeError as exc:
            raise ValueError(f"{stdout_log}: malformed DIAG JSON") from exc
        if event.get("event") == "md_complete":
            complete_events.append(event)
    if len(complete_events) != 1 or complete_events[0].get("status") != "SUCCESS_WT":
        raise ValueError(f"{stdout_log}: expected one SUCCESS_WT md_complete event, got {complete_events}")
    if "성공 (야생형 MD) : 1 개" not in stdout or "실패 (폭발 등) : 0 개" not in stdout:
        raise ValueError(f"{stdout_log}: WT batch summary did not report one success and zero failures")
    forbidden = ("Traceback", "CUDA_ERROR", "PARTIAL_SUCCESS", '"status": "FAIL"', "[EXPLODED]")
    hits = [pattern for pattern in forbidden if pattern in stdout]
    if hits:
        raise ValueError(f"{stdout_log}: forbidden failure markers {hits}")
    exploded = sorted(str(path) for path in output_root.rglob("*EXPLODED*"))
    partial = sorted(str(path) for path in output_root.rglob("*partial*"))
    if exploded or partial:
        raise ValueError(f"{output_root}: exploded={exploded}, partial={partial}")

    pdb_info = inspect_pdb(final_pdb)
    dcd_info = inspect_dcd(dcd)
    if dcd_info["n_atoms"] != pdb_info["n_atoms"]:
        raise ValueError(
            f"{output_root}: DCD atoms {dcd_info['n_atoms']} != final PDB atoms {pdb_info['n_atoms']}"
        )
    return {
        "md_log": parse_md_log(md_log),
        "pdb": pdb_info,
        "dcd": dcd_info,
        "artifacts": {
            "final_pdb": artifact(final_pdb),
            "md_log": artifact(md_log),
            "dcd": artifact(dcd),
            "stdout_log": artifact(stdout_log),
        },
        "md_complete_event": complete_events[0],
    }


def validate_existing_seed(seed: int) -> dict[str, Any]:
    root = REPO_ROOT / f"outputs/1YCR_WT_calib_s{seed}"
    input_path = root / "_md_input/1YCR_WT.pdb"
    final_pdb = root / "mdresult/1YCR_WT_final.pdb"
    md_log = root / "mdresult/1YCR_WT_md.log"
    dcd = DCD_ARCHIVE_ROOT / f"1YCR_WT_calib_s{seed}/mdresult/1YCR_WT_restrained.dcd"
    if sha256_file(input_path) != EXPECTED_INPUT_SHA256:
        raise ValueError(f"s{seed}: existing input hash mismatch")
    pdb_info = inspect_pdb(final_pdb)
    dcd_info = inspect_dcd(dcd)
    if pdb_info["n_atoms"] != dcd_info["n_atoms"]:
        raise ValueError(f"s{seed}: archived DCD/final PDB atom mismatch")
    return {
        "seed": seed,
        "input": artifact(input_path),
        "final_pdb": artifact(final_pdb),
        "md_log": {"artifact": artifact(md_log), "validation": parse_md_log(md_log)},
        "archived_dcd": {"artifact": artifact(dcd), "validation": dcd_info},
        "pdb_validation": pdb_info,
    }


def inventory_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("inventory_digest", None)
    encoded = json.dumps(canonical, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_source_inventory() -> dict[str, Any]:
    if sha256_file(MD_SCRIPT) != EXPECTED_MD_SCRIPT_SHA256:
        raise ValueError("current run_restrained_md.py SHA-256 drift")
    if sha256_file(SOURCE_INPUT) != EXPECTED_INPUT_SHA256:
        raise ValueError("1YCR WT common input SHA-256 drift")
    if run_checked(["git", "hash-object", relative(MD_SCRIPT)]) != CURRENT_MD_BLOB:
        raise ValueError("current run_restrained_md.py git blob drift")
    if run_checked(["git", "rev-parse", f"{HISTORICAL_COMMIT}:utils/run_restrained_md.py"]) != HISTORICAL_MD_BLOB:
        raise ValueError("historical run_restrained_md.py blob drift")

    try:
        import pdbfixer

        pdbfixer_path = pdbfixer.__file__
    except ImportError:
        pdbfixer_path = None
    payload = {
        "schema": "1ycr_wt_seed_expansion_source_inventory_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "preregistration": artifact(PREREGISTRATION),
        "runner": artifact(Path(__file__).resolve()),
        "md_generator": {
            "artifact": artifact(MD_SCRIPT),
            "current_blob": CURRENT_MD_BLOB,
            "historical_commit": HISTORICAL_COMMIT,
            "historical_blob": HISTORICAL_MD_BLOB,
        },
        "common_input": artifact(SOURCE_INPUT),
        "environment": {
            "python": sys.version.split()[0],
            "executable": sys.executable,
            "openmm": openmm.version.version,
            "numpy": np.__version__,
            "mdtraj": md.__version__,
            "pdbfixer_module": pdbfixer_path,
        },
        "existing_scaffolds": [validate_existing_seed(seed) for seed in EXISTING_SEEDS],
        "new_seeds": list(NEW_SEEDS),
        "protocol": {
            "steps": N_STEPS,
            "topology": "linear",
            "ncaa_label": "none",
            "platform": "CUDA",
            "dt_fs": 2.0,
            "device": "0",
            "dcd_interval_override": None,
        },
    }
    payload["inventory_digest"] = inventory_digest(payload)
    return payload


def verify_source_inventory() -> dict[str, Any]:
    if not SOURCE_INVENTORY.is_file():
        raise FileNotFoundError(f"missing {SOURCE_INVENTORY}; run --prepare-inventory")
    frozen = json.loads(SOURCE_INVENTORY.read_text(encoding="utf-8"))
    if frozen.get("inventory_digest") != inventory_digest(frozen):
        raise ValueError("stored source inventory digest mismatch")
    current = build_source_inventory()
    if frozen["inventory_digest"] != current["inventory_digest"]:
        raise ValueError(
            "source/code/preregistration drift after freeze: "
            f"stored={frozen['inventory_digest']} current={current['inventory_digest']}"
        )
    return frozen


def output_root(seed: int) -> Path:
    return REPO_ROOT / f"outputs/1YCR_WT_calib_s{seed}"


def completed_manifest_path(seed: int) -> Path:
    return output_root(seed) / "scaffold_manifest.json"


def validate_completed_seed(seed: int, inventory: dict[str, Any]) -> dict[str, Any] | None:
    path = completed_manifest_path(seed)
    if not path.is_file():
        return None
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("status") != "COMPLETE" or manifest.get("seed") != seed:
        raise ValueError(f"s{seed}: invalid completed scaffold manifest")
    if manifest.get("source_inventory_digest") != inventory["inventory_digest"]:
        raise ValueError(f"s{seed}: completed scaffold inventory digest mismatch")
    input_path = output_root(seed) / "_md_input/1YCR_WT.pdb"
    if sha256_file(input_path) != EXPECTED_INPUT_SHA256:
        raise ValueError(f"s{seed}: staged input hash mismatch")
    current_validation = validate_outputs(output_root(seed), output_root(seed) / "scaffold_run.log")
    if current_validation["artifacts"] != manifest.get("validation", {}).get("artifacts"):
        raise ValueError(f"s{seed}: completed output artifact drift")
    return manifest


def archive_partial(seed: int) -> Path | None:
    root = output_root(seed)
    if not root.exists():
        return None
    if (root / "mdresult/1YCR_WT_final.pdb").exists():
        raise ValueError(f"s{seed}: ambiguous root has final PDB without a valid completed manifest")
    archive_root = (
        REPO_ROOT
        / "outputs/_archive/1ycr_wt_seed_expansion_20260716"
        / f"{time.strftime('%Y%m%d_%H%M%S')}_s{seed}"
    )
    archive_root.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(root), str(archive_root))
    atomic_write_json(
        archive_root / "archive_reason.json",
        {
            "seed": seed,
            "archived": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "reason": "interrupted scaffold root without final PDB/completed manifest",
        },
    )
    return archive_root


def assert_gpu_idle() -> None:
    result = subprocess.run(
        [
            "nvidia-smi",
            "--query-compute-apps=pid,process_name,used_memory",
            "--format=csv,noheader,nounits",
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
    )
    if result.stdout.strip():
        raise RuntimeError(f"GPU compute process already active:\n{result.stdout.strip()}")


def command_for_seed(seed: int) -> list[str]:
    root = output_root(seed)
    return [
        QMMM_PYTHON,
        str(MD_SCRIPT),
        "--inputdir",
        str(root / "_md_input"),
        "--outputdir",
        str(root / "mdresult"),
        "--steps",
        str(N_STEPS),
        "--topology",
        "linear",
        "--ncaa_label",
        "none",
        "--seed",
        str(seed),
        "--platform",
        "CUDA",
        "--dt_fs",
        "2.0",
    ]


def write_runtime_state(status: str, **extra: Any) -> None:
    atomic_write_json(
        RUNTIME_ROOT / "run_state.json",
        {"status": status, "updated": time.strftime("%Y-%m-%d %H:%M:%S %z"), **extra},
    )


def write_summary(manifests: list[dict[str, Any]]) -> dict[str, Any]:
    by_seed = {manifest["seed"]: manifest for manifest in manifests}
    pending = [seed for seed in NEW_SEEDS if seed not in by_seed]
    elapsed = [float(manifest["elapsed_s"]) for manifest in manifests]
    payload = {
        "schema": "1ycr_wt_seed_expansion_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "COMPLETE" if not pending else "RUNNING",
        "n_expected": len(NEW_SEEDS),
        "n_completed": len(manifests),
        "completed_seeds": [seed for seed in NEW_SEEDS if seed in by_seed],
        "pending_seeds": pending,
        "median_elapsed_s": median(elapsed) if elapsed else None,
        "eta_s": median(elapsed) * len(pending) if elapsed else None,
        "cells": [by_seed[seed] for seed in NEW_SEEDS if seed in by_seed],
    }
    atomic_write_json(SUMMARY_JSON, payload)
    atomic_write_json(RUNTIME_ROOT / "progress.json", payload)
    lines = [
        "# 1YCR WT Scaffold Expansion",
        "",
        f"Status: `{payload['status']}`",
        f"Completed: `{payload['n_completed']}/{payload['n_expected']}`",
        f"Pending seeds: `{','.join(str(seed) for seed in pending) if pending else 'none'}`",
        f"ETA seconds: `{payload['eta_s'] if payload['eta_s'] is not None else 'NA'}`",
        "",
        "| seed | elapsed min | atoms | DCD frames | temp K min/max | final step |",
        "| ---: | ---: | ---: | ---: | --- | ---: |",
    ]
    for seed in NEW_SEEDS:
        if seed not in by_seed:
            continue
        manifest = by_seed[seed]
        validation = manifest["validation"]
        md_log = validation["md_log"]
        lines.append(
            f"| {seed} | {manifest['elapsed_s'] / 60.0:.1f} | {validation['pdb']['n_atoms']} | "
            f"{validation['dcd']['n_frames']} | {md_log['temperature_min_K']:.1f}/"
            f"{md_log['temperature_max_K']:.1f} | {md_log['final_step']} |"
        )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return payload


def run_seed(seed: int, inventory: dict[str, Any]) -> dict[str, Any]:
    archive = archive_partial(seed)
    if archive is not None:
        print(f"[1YCR] archived partial s{seed} at {archive}", flush=True)
    root = output_root(seed)
    input_dir = root / "_md_input"
    input_dir.mkdir(parents=True, exist_ok=True)
    (root / "mdresult").mkdir(parents=True, exist_ok=True)
    staged_input = input_dir / "1YCR_WT.pdb"
    shutil.copy2(SOURCE_INPUT, staged_input)
    if sha256_file(staged_input) != EXPECTED_INPUT_SHA256:
        raise ValueError(f"s{seed}: staged common input hash mismatch")

    assert_gpu_idle()
    command = command_for_seed(seed)
    stdout_log = root / "scaffold_run.log"
    env = os.environ.copy()
    env["UPDD_MD_CUDA_DEVICE"] = "0"
    env.pop("UPDD_MD_DCD_INTERVAL", None)
    started = time.time()
    with stdout_log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
        )
    elapsed = time.time() - started
    if result.returncode != 0:
        raise RuntimeError(f"s{seed}: child exit {result.returncode}")
    validation = validate_outputs(root, stdout_log)
    manifest = {
        "schema": "1ycr_wt_scaffold_manifest_v1",
        "status": "COMPLETE",
        "seed": seed,
        "source_inventory_digest": inventory["inventory_digest"],
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "elapsed_s": elapsed,
        "command": command,
        "environment": {
            "UPDD_MD_CUDA_DEVICE": "0",
            "UPDD_MD_DCD_INTERVAL": None,
        },
        "staged_input": artifact(staged_input),
        "validation": validation,
    }
    atomic_write_json(completed_manifest_path(seed), manifest)
    return manifest


def _pause_at_seed_boundary(manifests: list[dict[str, Any]]) -> bool:
    progress = {
        "boundary": "scaffold_seed",
        "completed": len(manifests),
        "expected": len(NEW_SEEDS),
    }
    if not boundary_pause_requested(progress):
        return False
    write_runtime_state(
        "PAUSED",
        n_completed=len(manifests),
        n_expected=len(NEW_SEEDS),
        boundary="scaffold_seed",
    )
    print("[1YCR] cooperative pause acknowledged at seed boundary", flush=True)
    return True


def run_all(inventory: dict[str, Any]) -> bool:
    manifests = []
    for seed in NEW_SEEDS:
        completed = validate_completed_seed(seed, inventory)
        if completed is not None:
            manifests.append(completed)
    write_summary(manifests)
    if _pause_at_seed_boundary(manifests):
        return True

    for seed in NEW_SEEDS:
        if any(manifest["seed"] == seed for manifest in manifests):
            print(f"[1YCR] resume-skip s{seed}", flush=True)
            continue
        write_runtime_state("RUNNING", current_seed=seed, n_completed=len(manifests), n_expected=len(NEW_SEEDS))
        print(f"[1YCR] running s{seed}", flush=True)
        try:
            manifest = run_seed(seed, inventory)
        except Exception as exc:
            write_runtime_state(
                "FAILED",
                current_seed=seed,
                n_completed=len(manifests),
                n_expected=len(NEW_SEEDS),
                error=f"{type(exc).__name__}: {exc}",
            )
            raise
        manifests.append(manifest)
        write_summary(manifests)
        if _pause_at_seed_boundary(manifests):
            return True
    write_runtime_state("COMPLETE", n_completed=len(manifests), n_expected=len(NEW_SEEDS))
    return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-inventory", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args(argv)

    if args.prepare_inventory:
        inventory = build_source_inventory()
        atomic_write_json(SOURCE_INVENTORY, inventory)
        print(json.dumps({"inventory_digest": inventory["inventory_digest"], "n_existing": 5, "n_new": 5}))
        return 0

    inventory = verify_source_inventory()
    if args.dry_run:
        collisions = [seed for seed in NEW_SEEDS if output_root(seed).exists()]
        if collisions:
            raise ValueError(f"new canonical roots already exist: {collisions}")
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "inventory_digest": inventory["inventory_digest"],
                    "new_seeds": list(NEW_SEEDS),
                    "commands": {str(seed): command_for_seed(seed) for seed in NEW_SEEDS},
                },
                indent=2,
            )
        )
        return 0

    if args.analyze_only:
        manifests = []
        for seed in NEW_SEEDS:
            completed = validate_completed_seed(seed, inventory)
            if completed is not None:
                manifests.append(completed)
        write_summary(manifests)
        return 0

    paused = run_all(inventory)
    return PAUSE_EXIT_CODE if paused else 0


if __name__ == "__main__":
    raise SystemExit(main())
