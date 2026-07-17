#!/usr/bin/env python3
"""Replay rejected P1SD trajectories for raw ATM cap-activation diagnostics."""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
UTILS_DIR = REPO_ROOT / "utils"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))

from analysis.dynamic_ghost_apex_bridge_p0v2_20260717 import (  # noqa: E402
    run_p0v2_reference_audit as p0v2,
)
from analysis.dynamic_ghost_union_bridge_p1sd_20260717 import (  # noqa: E402
    run_p1sd as p1sd,
)
import trackb_apex_bridge as bridge  # noqa: E402
import trackb_dynamic_ghost as dg  # noqa: E402


ANALYSIS_DIR = Path(__file__).resolve().parent
PREREGISTRATION = ANALYSIS_DIR / "PREREGISTRATION.md"
PROTOCOL_PATH = ANALYSIS_DIR / "protocol_r1.json"
FROZEN_MANIFEST = ANALYSIS_DIR / "FROZEN_MANIFEST_R1.json"
RUNNER_PATH = Path(__file__).resolve()

P1_DIR = REPO_ROOT / "analysis/dynamic_ghost_union_bridge_p1sd_20260717"
P1_INVENTORY = P1_DIR / "P1SD_PREBUILD_INVENTORY.json"
P1_NORMALIZED_SUMMARY = P1_DIR / "P1SD_NORMALIZED_SOURCE_SUMMARY.json"
P1_PARTIAL_SUMMARY = P1_DIR / "P1SD_PARTIAL_SUMMARY.json"
P1_PATH_DIAGNOSIS = P1_DIR / "P1SD_PATH_DIAGNOSIS.json"

PREBUILD_INVENTORY = ANALYSIS_DIR / "P1CAP_PREBUILD_INVENTORY_R1.json"
KEEPER_PREFLIGHT = ANALYSIS_DIR / "P1CAP_KEEPER_PREFLIGHT_R1.json"
REPLAY_ROOT = ANALYSIS_DIR / "raw_gap_replay_r1"
RUN_STATE = ANALYSIS_DIR / "P1CAP_RUN_STATE_R1.json"
RUNNER_SUMMARY = ANALYSIS_DIR / "P1CAP_RUNNER_SUMMARY_R1.json"
RUNNER_REPORT = ANALYSIS_DIR / "P1CAP_RUNNER_SUMMARY_R1.md"

EXPECTED_ATM_PYTHON = Path("/home/san/miniconda3/envs/atm/bin/python")
EXPECTED_STATUS = "FROZEN_AFTER_R0_DCD_UNIT_FAILURE_BEFORE_R1_OUTPUT"

relative = p1sd.relative
sha256_file = p1sd.sha256_file
canonical_digest = p1sd.canonical_digest
artifact_record = p1sd.artifact_record
atomic_write_json = p1sd.atomic_write_json
atomic_write_text = p1sd.atomic_write_text


def _hash_map(base: Path, mapping: dict[str, str], label: str) -> None:
    for name, expected in mapping.items():
        path = base / name
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(f"{label} drift: {path} {observed} != {expected}")


def _require_atm_python() -> None:
    if Path(sys.executable).resolve() != EXPECTED_ATM_PYTHON.resolve():
        raise RuntimeError(
            f"P1CAP requires {EXPECTED_ATM_PYTHON}, got {sys.executable}"
        )


def load_protocol_and_verify_freeze() -> dict[str, Any]:
    frozen = json.loads(FROZEN_MANIFEST.read_text(encoding="utf-8"))
    if frozen.get("status") != EXPECTED_STATUS:
        raise ValueError("P1CAP frozen manifest status drifted")
    _hash_map(ANALYSIS_DIR, frozen.get("files") or {}, "P1CAP frozen file")
    _hash_map(REPO_ROOT, frozen.get("parent_evidence") or {}, "P1CAP parent evidence")

    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    if protocol.get("status") != EXPECTED_STATUS:
        raise ValueError("P1CAP protocol status drifted")
    expected = protocol.get("expected") or {}
    if expected != {
        "cells": 4,
        "windows": 16,
        "frames_per_window": 50,
        "frames": 800,
        "source_sample_rows": 800,
    }:
        raise ValueError("P1CAP expected cohort contract drifted")
    cohort = protocol.get("cohort") or []
    if [row.get("cell_id") for row in cohort] != [
        "w4a_union_s101_bound",
        "w4a_union_s127_bound",
        "w4a_union_s163_bound",
        "w4a_union_s101_free",
    ]:
        raise ValueError("P1CAP cohort or order drifted")
    if sum(len(row.get("xi") or []) for row in cohort) != 16:
        raise ValueError("P1CAP window count drifted")
    evaluation = protocol.get("evaluation") or {}
    if evaluation.get("md_steps") != 0 or evaluation.get("minimization_steps") != 0:
        raise ValueError("P1CAP must remain a zero-dynamics replay")
    if evaluation.get("worker_unit") != "cell":
        raise ValueError("P1CAP worker isolation contract drifted")
    if evaluation.get("one_context_per_worker") is not True:
        raise ValueError("P1CAP Context contract drifted")
    if (
        evaluation.get("dcd_reader") != "mdtraj.formats.DCDTrajectoryFile"
        or evaluation.get("dcd_distance_unit") != "angstroms"
        or evaluation.get("dcd_to_nm_scale") != 0.1
    ):
        raise ValueError("P1CAP R1 DCD unit contract drifted")
    softcore = protocol.get("softcore") or {}
    if softcore != {
        "uoffset_kj_mol": 0.0,
        "ubcore_kj_mol": 418.4,
        "umax_kj_mol": 836.8,
        "acore": 0.0625,
        "classification_guard_kj_mol": 5.0,
        "linear_rule": "abs(delta_raw)<=ubcore-guard",
        "ambiguous_rule": "ubcore-guard<abs(delta_raw)<ubcore+guard",
        "cap_active_rule": "abs(delta_raw)>=ubcore+guard",
    }:
        raise ValueError("P1CAP soft-core contract drifted")
    if any((protocol.get("claims") or {}).values()):
        raise ValueError("P1CAP prohibited claim was enabled")

    parent_inventory = p1sd.load_inventory()
    partial = json.loads(P1_PARTIAL_SUMMARY.read_text(encoding="utf-8"))
    path = json.loads(P1_PATH_DIAGNOSIS.read_text(encoding="utf-8"))
    declared_parent = protocol.get("parent") or {}
    if parent_inventory.get("inventory_digest") != declared_parent.get(
        "inventory_digest"
    ):
        raise ValueError("P1CAP parent inventory digest drifted")
    if partial.get("status") != declared_parent.get("partial_status"):
        raise ValueError("P1CAP parent partial status drifted")
    if path.get("path_verdict") != declared_parent.get("path_verdict"):
        raise ValueError("P1CAP parent PATH verdict drifted")
    execution = partial.get("execution") or {}
    if (
        execution.get("completed_window_results") != 16
        or execution.get("sample_rows") != 800
        or execution.get("dcd_frames") != 800
    ):
        raise ValueError("P1CAP parent completed-artifact counts drifted")
    if declared_parent.get("pilot_samples_reusable") is not False:
        raise ValueError("P1CAP parent-sample reuse prohibition drifted")
    return protocol


def inventory_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("inventory_digest", None)
    return canonical_digest(canonical)


def keeper_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("keeper_digest", None)
    return canonical_digest(canonical)


def summary_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("summary_digest", None)
    return canonical_digest(canonical)


def _implementation_records() -> dict[str, dict[str, Any]]:
    paths = {
        "runner": RUNNER_PATH,
        "p1sd_runner": P1_DIR / "run_p1sd.py",
        "p0v2_runner": (
            REPO_ROOT
            / "analysis/dynamic_ghost_apex_bridge_p0v2_20260717"
            / "run_p0v2_reference_audit.py"
        ),
        "apex_bridge": REPO_ROOT / "utils/trackb_apex_bridge.py",
        "dynamic_ghost": REPO_ROOT / "utils/trackb_dynamic_ghost.py",
    }
    return {name: artifact_record(path) for name, path in paths.items()}


def _cohort(protocol: dict[str, Any]) -> list[dict[str, Any]]:
    return list(protocol["cohort"])


def _cell_row(inventory: dict[str, Any], cell_id: str) -> dict[str, Any]:
    rows = [row for row in inventory["cells"] if row["cell_id"] == cell_id]
    if len(rows) != 1:
        raise ValueError(f"expected one P1CAP cell row for {cell_id}, found {len(rows)}")
    return rows[0]


def _verify_artifact(record: dict[str, Any]) -> Path:
    path = REPO_ROOT / record["path"]
    if not path.is_file():
        raise FileNotFoundError(path)
    size = path.stat().st_size
    if size != int(record["size"]):
        raise ValueError(f"artifact size drift: {path} {size} != {record['size']}")
    observed = sha256_file(path)
    if observed != record["sha256"]:
        raise ValueError(
            f"artifact hash drift: {path} {observed} != {record['sha256']}"
        )
    return path


def _assert_declared_artifact(record: dict[str, Any], declared: dict[str, Any]) -> None:
    if record != declared:
        raise ValueError(
            f"artifact declaration drift: {record.get('path')} != {declared.get('path')}"
        )


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        return list(csv.DictReader(handle))


def _dcd_metadata(path: Path) -> dict[str, Any]:
    from mdtraj.formats import DCDTrajectoryFile

    with DCDTrajectoryFile(str(path), mode="r") as reader:
        n_frames = len(reader)
        distance_unit = str(reader.distance_unit)
        xyz, lengths, angles = reader.read(n_frames=1)
    if xyz.shape[0] != 1:
        raise ValueError(f"could not read first DCD frame: {path}")
    return {
        "frames": n_frames,
        "particles": int(xyz.shape[1]),
        "distance_unit": distance_unit,
        "first_lengths_raw": [float(value) for value in lengths[0]],
        "first_angles_deg": [float(value) for value in angles[0]],
    }


def build_inventory() -> dict[str, Any]:
    _require_atm_python()
    protocol = load_protocol_and_verify_freeze()
    if PREBUILD_INVENTORY.exists():
        raise FileExistsError(PREBUILD_INVENTORY)
    parent_inventory = p1sd.load_inventory()
    cells = []
    total_windows = 0
    total_frames = 0
    for spec in _cohort(protocol):
        cell_id = spec["cell_id"]
        parent_cell = p1sd.frozen_cell(parent_inventory, cell_id)
        normalized_dir = P1_DIR / "normalized_sources" / cell_id
        normalization_path = normalized_dir / "normalization_result.json"
        normalization = json.loads(normalization_path.read_text(encoding="utf-8"))
        if (
            normalization.get("status") != "PASS"
            or normalization.get("inventory_digest")
            != parent_inventory["inventory_digest"]
            or normalization.get("cell_id") != cell_id
        ):
            raise ValueError(f"{cell_id}: invalid P1SD normalization result")
        source = {
            "normalization_result": artifact_record(normalization_path),
            "bridge_system": artifact_record(
                normalized_dir / "normalized_bridge_system.xml"
            ),
            "topology_pdb": artifact_record(normalized_dir / "normalized_source.pdb"),
        }
        _assert_declared_artifact(
            source["bridge_system"], normalization["artifacts"]["bridge_system"]
        )
        _assert_declared_artifact(
            source["topology_pdb"], normalization["artifacts"]["inspection_pdb"]
        )

        windows = []
        for xi, expected_status in zip(
            spec["xi"], spec["expected_window_status"], strict=True
        ):
            window_dir = P1_DIR / "sampled_schedule_discovery" / cell_id / p1sd.xi_label(xi)
            result_path = window_dir / "window_result.json"
            result = json.loads(result_path.read_text(encoding="utf-8"))
            if (
                result.get("cell_id") != cell_id
                or float(result.get("xi")) != float(xi)
                or result.get("status") != expected_status
                or result.get("inventory_digest") != parent_inventory["inventory_digest"]
                or result.get("n_rows") != protocol["expected"]["frames_per_window"]
                or result.get("dcd_frames") != protocol["expected"]["frames_per_window"]
            ):
                raise ValueError(f"{cell_id} xi={xi:.2f}: invalid P1SD window result")
            artifacts = {
                "window_result": artifact_record(result_path),
                "samples": artifact_record(window_dir / "samples.csv"),
                "trajectory": artifact_record(window_dir / "trajectory.dcd"),
                "final_positions": artifact_record(window_dir / "final_positions.npy"),
            }
            for name in ("samples", "trajectory", "final_positions"):
                _assert_declared_artifact(artifacts[name], result["artifacts"][name])
            sample_rows = _read_csv_rows(window_dir / "samples.csv")
            dcd = _dcd_metadata(window_dir / "trajectory.dcd")
            if (
                len(sample_rows) != 50
                or dcd["frames"] != 50
                or dcd["distance_unit"]
                != protocol["evaluation"]["dcd_distance_unit"]
            ):
                raise ValueError(f"{cell_id} xi={xi:.2f}: source row/frame mismatch")
            windows.append(
                {
                    "xi": float(xi),
                    "expected_status": expected_status,
                    "source_rows": len(sample_rows),
                    "dcd_frames": dcd["frames"],
                    "dcd_particles": dcd["particles"],
                    "dcd_distance_unit": dcd["distance_unit"],
                    "dcd_first_lengths_raw": dcd["first_lengths_raw"],
                    "dcd_first_angles_deg": dcd["first_angles_deg"],
                    "artifacts": artifacts,
                }
            )
            total_windows += 1
            total_frames += dcd["frames"]
        cells.append(
            {
                "cell_id": cell_id,
                "seed": parent_cell["seed"],
                "leg": spec["leg"],
                "particle_count": normalization["ghost_contract"]["n_particles"],
                "parent_box": parent_cell["parent_box"],
                "state": parent_cell["states"]["apex_dplus"],
                "source": source,
                "windows": windows,
            }
        )
    if total_windows != 16 or total_frames != 800:
        raise ValueError("P1CAP inventory cohort total drifted")

    import mdtraj
    import openmm

    payload = {
        "schema": "updd_raw_gap_cap_activation_p1cap_inventory_v1_r1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "FROZEN_BEFORE_P1CAP_R1_OUTPUT",
        "python": sys.executable,
        "runtime": {
            "openmm": openmm.__version__,
            "mdtraj": mdtraj.__version__,
            "numpy": np.__version__,
        },
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "frozen_manifest_sha256": sha256_file(FROZEN_MANIFEST),
        "implementation": _implementation_records(),
        "parent_inventory_digest": parent_inventory["inventory_digest"],
        "counts": {
            "cells": len(cells),
            "windows": total_windows,
            "frames": total_frames,
        },
        "cells": cells,
    }
    payload["inventory_digest"] = inventory_digest(payload)
    atomic_write_json(PREBUILD_INVENTORY, payload)
    return payload


def _iter_input_artifacts(
    inventory: dict[str, Any],
) -> Iterable[dict[str, Any]]:
    for cell in inventory["cells"]:
        yield from cell["source"].values()
        for window in cell["windows"]:
            yield from window["artifacts"].values()


def load_inventory(*, verify_inputs: bool = False) -> dict[str, Any]:
    protocol = load_protocol_and_verify_freeze()
    if not PREBUILD_INVENTORY.is_file():
        raise FileNotFoundError(f"missing {PREBUILD_INVENTORY}; run --build-inventory")
    inventory = json.loads(PREBUILD_INVENTORY.read_text(encoding="utf-8"))
    if inventory.get("status") != "FROZEN_BEFORE_P1CAP_R1_OUTPUT":
        raise ValueError("P1CAP inventory status drifted")
    if inventory_digest(inventory) != inventory.get("inventory_digest"):
        raise ValueError("P1CAP inventory digest drifted")
    if inventory.get("implementation") != _implementation_records():
        raise ValueError("P1CAP implementation drift after inventory freeze")
    if inventory.get("protocol_sha256") != sha256_file(PROTOCOL_PATH):
        raise ValueError("P1CAP protocol drift after inventory freeze")
    if inventory.get("frozen_manifest_sha256") != sha256_file(FROZEN_MANIFEST):
        raise ValueError("P1CAP frozen manifest drift after inventory freeze")
    if inventory.get("counts") != {
        "cells": protocol["expected"]["cells"],
        "windows": protocol["expected"]["windows"],
        "frames": protocol["expected"]["frames"],
    }:
        raise ValueError("P1CAP inventory counts drifted")
    if verify_inputs:
        for record in _iter_input_artifacts(inventory):
            _verify_artifact(record)
    return inventory


def run_keeper_preflight() -> dict[str, Any]:
    _require_atm_python()
    protocol = load_protocol_and_verify_freeze()
    if KEEPER_PREFLIGHT.exists():
        raise FileExistsError(KEEPER_PREFLIGHT)
    if REPLAY_ROOT.exists() or RUNNER_SUMMARY.exists():
        raise RuntimeError("P1CAP output exists before KEEPER preflight")
    inventory = load_inventory(verify_inputs=True)

    import openmm as mm
    from openmm.app import PDBFile

    cell_checks = []
    total_windows = 0
    total_frames = 0
    for cell in inventory["cells"]:
        cell_id = cell["cell_id"]
        source = cell["source"]
        system_path = REPO_ROOT / source["bridge_system"]["path"]
        pdb_path = REPO_ROOT / source["topology_pdb"]["path"]
        system = mm.XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
        bridge_contract = bridge.inspect_apex_bridge(system)
        pdb = PDBFile(str(pdb_path))
        system_particles = int(system.getNumParticles())
        topology_particles = sum(1 for _atom in pdb.topology.atoms())
        box = np.asarray(
            [
                [vector.x, vector.y, vector.z]
                for vector in system.getDefaultPeriodicBoxVectors()
            ],
            dtype=float,
        )
        target_box = np.asarray(cell["parent_box"]["vectors_nm"], dtype=float)
        window_checks = []
        for window in cell["windows"]:
            samples_path = REPO_ROOT / window["artifacts"]["samples"]["path"]
            trajectory_path = REPO_ROOT / window["artifacts"]["trajectory"]["path"]
            final_path = REPO_ROOT / window["artifacts"]["final_positions"]["path"]
            sample_rows = _read_csv_rows(samples_path)
            dcd = _dcd_metadata(trajectory_path)
            final_positions = np.load(final_path, mmap_mode="r", allow_pickle=False)
            converted_box = _box_vectors(
                np.asarray(dcd["first_lengths_raw"], dtype=float)
                * protocol["evaluation"]["dcd_to_nm_scale"],
                np.asarray(dcd["first_angles_deg"], dtype=float),
            )
            checks = {
                "source_rows_exact": len(sample_rows)
                == protocol["expected"]["frames_per_window"],
                "dcd_frames_exact": dcd["frames"]
                == protocol["expected"]["frames_per_window"],
                "dcd_particles_exact": dcd["particles"] == system_particles,
                "dcd_unit_exact": dcd["distance_unit"]
                == protocol["evaluation"]["dcd_distance_unit"],
                "converted_dcd_box_exact": np.allclose(
                    converted_box,
                    target_box,
                    atol=protocol["validation"]["dcd_box_absolute_tolerance_nm"],
                    rtol=0.0,
                ),
                "final_positions_shape_exact": tuple(final_positions.shape)
                == (system_particles, 3),
            }
            window_checks.append(
                {
                    "xi": window["xi"],
                    "checks": checks,
                    "passed": all(checks.values()),
                }
            )
            total_windows += 1
            total_frames += dcd["frames"]
            del final_positions
        state = cell["state"]
        state_softcore = {
            "uoffset_kj_mol": 0.0,
            "ubcore_kj_mol": state["ubcore_kcal_mol"] * p0v2.KCAL_TO_KJ,
            "umax_kj_mol": state["umax_kcal_mol"] * p0v2.KCAL_TO_KJ,
            "acore": state["acore"],
        }
        expected_softcore = {
            key: protocol["softcore"][key]
            for key in (
                "uoffset_kj_mol",
                "ubcore_kj_mol",
                "umax_kj_mol",
                "acore",
            )
        }
        checks = {
            "normalization_pass": json.loads(
                (REPO_ROOT / source["normalization_result"]["path"]).read_text(
                    encoding="utf-8"
                )
            ).get("status")
            == "PASS",
            "system_particles_exact": system_particles == cell["particle_count"],
            "topology_particles_exact": topology_particles == system_particles,
            "bridge_contract_exact": bridge_contract["energy_function"]
            == bridge.BRIDGE_ATM_EXPRESSION,
            "bridge_derivative_declared": bridge.BRIDGE_PARAMETER
            in bridge_contract["energy_parameter_derivatives"],
            "fixed_parent_box": np.allclose(box, target_box, atol=1e-10, rtol=0.0),
            "softcore_parameters_exact": all(
                math.isclose(
                    float(state_softcore[name]),
                    float(expected_softcore[name]),
                    abs_tol=1e-9,
                    rel_tol=0.0,
                )
                for name in expected_softcore
            ),
            "windows_exact": len(window_checks) == len(cell["windows"]),
            "window_artifacts_valid": all(row["passed"] for row in window_checks),
        }
        cell_checks.append(
            {
                "cell_id": cell_id,
                "system_particles": system_particles,
                "topology_particles": topology_particles,
                "softcore": state_softcore,
                "windows": window_checks,
                "checks": checks,
                "passed": all(checks.values()),
            }
        )
        del pdb, system
        gc.collect()
    overall_checks = {
        "cells_exact": len(cell_checks) == protocol["expected"]["cells"],
        "windows_exact": total_windows == protocol["expected"]["windows"],
        "frames_exact": total_frames == protocol["expected"]["frames"],
        "all_cells_pass": all(row["passed"] for row in cell_checks),
        "context_count_zero": True,
    }
    payload = {
        "schema": "updd_raw_gap_cap_activation_p1cap_keeper_v1_r1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "PASS" if all(overall_checks.values()) else "BLOCK",
        "inventory_digest": inventory["inventory_digest"],
        "context_count": 0,
        "counts": {
            "cells": len(cell_checks),
            "windows": total_windows,
            "frames": total_frames,
        },
        "cell_checks": cell_checks,
        "checks": overall_checks,
    }
    payload["keeper_digest"] = keeper_digest(payload)
    atomic_write_json(KEEPER_PREFLIGHT, payload)
    if payload["status"] != "PASS":
        raise RuntimeError("P1CAP KEEPER BLOCK")
    return payload


def load_keeper(inventory: dict[str, Any]) -> dict[str, Any]:
    if not KEEPER_PREFLIGHT.is_file():
        raise FileNotFoundError(
            f"missing {KEEPER_PREFLIGHT}; run --keeper-preflight"
        )
    keeper = json.loads(KEEPER_PREFLIGHT.read_text(encoding="utf-8"))
    if keeper_digest(keeper) != keeper.get("keeper_digest"):
        raise ValueError("P1CAP KEEPER digest drifted")
    if keeper.get("status") != "PASS" or keeper.get("context_count") != 0:
        raise ValueError("P1CAP KEEPER is not a context-free PASS")
    if keeper.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError("P1CAP KEEPER inventory drifted")
    return keeper


def softcore_value(
    value_kj_mol: float,
    *,
    ubcore_kj_mol: float,
    umax_kj_mol: float,
    acore: float,
) -> float:
    value = float(value_kj_mol)
    ubcore = float(ubcore_kj_mol)
    umax = float(umax_kj_mol)
    a = float(acore)
    if not (a > 0.0 and umax > ubcore):
        raise ValueError("invalid P1CAP soft-core parameters")
    if value <= ubcore:
        return value
    y = (value - ubcore) / (umax - ubcore)
    z = 1.0 + 2.0 * (y / a) + 2.0 * (y / a) ** 2
    f = (z**a - 1.0) / (z**a + 1.0)
    return (umax - ubcore) * f + ubcore


def analytic_bridge_slope(
    delta_raw_kj_mol: float,
    *,
    uoffset_kj_mol: float,
    ubcore_kj_mol: float,
    umax_kj_mol: float,
    acore: float,
) -> float:
    delta_raw = float(delta_raw_kj_mol)
    adjusted = delta_raw - float(uoffset_kj_mol)
    plus = softcore_value(
        adjusted,
        ubcore_kj_mol=ubcore_kj_mol,
        umax_kj_mol=umax_kj_mol,
        acore=acore,
    )
    minus = softcore_value(
        -adjusted,
        ubcore_kj_mol=ubcore_kj_mol,
        umax_kj_mol=umax_kj_mol,
        acore=acore,
    )
    return delta_raw + 0.5 * (minus - plus)


def classify_cap(
    delta_raw_kj_mol: float,
    *,
    ubcore_kj_mol: float,
    guard_kj_mol: float,
) -> str:
    magnitude = abs(float(delta_raw_kj_mol))
    lower = float(ubcore_kj_mol) - float(guard_kj_mol)
    upper = float(ubcore_kj_mol) + float(guard_kj_mol)
    if magnitude <= lower:
        return "linear"
    if magnitude >= upper:
        return "cap_active"
    return "ambiguous"


def directional_cap(delta_raw_kj_mol: float, classification: str) -> str:
    if classification != "cap_active":
        return "none"
    return "dplus" if float(delta_raw_kj_mol) > 0.0 else "dminus"


def slope_tolerance(
    analytic_kj_mol: float, observed_kj_mol: float, validation: dict[str, Any]
) -> float:
    floor = float(validation["analytic_slope_absolute_floor_kj_mol"])
    relative = float(validation["analytic_slope_relative"])
    return max(
        floor,
        relative
        * max(1.0, abs(float(analytic_kj_mol)), abs(float(observed_kj_mol))),
    )


def strict_classification_agrees(first: str, second: str) -> bool:
    if "ambiguous" in {first, second}:
        return True
    return first == second


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write empty P1CAP CSV")
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def _box_vectors(
    lengths_nm: np.ndarray,
    angles_deg: np.ndarray,
) -> np.ndarray:
    from mdtraj.utils.unitcell import lengths_and_angles_to_box_vectors

    vectors = lengths_and_angles_to_box_vectors(
        np.asarray([lengths_nm[0]], dtype=float),
        np.asarray([lengths_nm[1]], dtype=float),
        np.asarray([lengths_nm[2]], dtype=float),
        np.asarray([angles_deg[0]], dtype=float),
        np.asarray([angles_deg[1]], dtype=float),
        np.asarray([angles_deg[2]], dtype=float),
    )
    return np.asarray([vector[0] for vector in vectors], dtype=float)


def _set_frame(
    context: Any,
    xyz_nm: np.ndarray,
    lengths_nm: np.ndarray,
    angles_deg: np.ndarray,
) -> np.ndarray:
    import openmm as mm
    from openmm import unit

    box = _box_vectors(lengths_nm, angles_deg)
    context.setPeriodicBoxVectors(
        *(mm.Vec3(*vector) * unit.nanometer for vector in box)
    )
    context.setPositions(np.asarray(xyz_nm, dtype=float) * unit.nanometer)
    return box


def _evaluate_frame(
    context: Any,
    atm_force: Any,
    softcore: dict[str, float],
    validation: dict[str, Any],
) -> dict[str, Any]:
    from openmm import unit

    raw_u1, raw_u0, _bias = atm_force.getPerturbationEnergy(context)
    u1 = float(raw_u1.value_in_unit(unit.kilojoule_per_mole))
    u0 = float(raw_u0.value_in_unit(unit.kilojoule_per_mole))
    delta = u1 - u0
    state = context.getState(getEnergy=True, getParameterDerivatives=True)
    energy = float(
        state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    )
    derivatives = {
        str(name): float(value)
        for name, value in dict(state.getEnergyParameterDerivatives()).items()
    }
    observed_slope = derivatives[bridge.BRIDGE_PARAMETER]
    analytic_slope = analytic_bridge_slope(
        delta,
        uoffset_kj_mol=softcore["uoffset_kj_mol"],
        ubcore_kj_mol=softcore["ubcore_kj_mol"],
        umax_kj_mol=softcore["umax_kj_mol"],
        acore=softcore["acore"],
    )
    allowed = slope_tolerance(analytic_slope, observed_slope, validation)
    error = abs(observed_slope - analytic_slope)
    classification = classify_cap(
        delta,
        ubcore_kj_mol=softcore["ubcore_kj_mol"],
        guard_kj_mol=softcore["classification_guard_kj_mol"],
    )
    values = (u0, u1, delta, energy, observed_slope, analytic_slope, error)
    return {
        "u0_kj_mol": u0,
        "u1_kj_mol": u1,
        "delta_raw_kj_mol": delta,
        "abs_delta_raw_kj_mol": abs(delta),
        "distance_from_ubcore_kj_mol": abs(delta)
        - softcore["ubcore_kj_mol"],
        "classification": classification,
        "directional_cap": directional_cap(delta, classification),
        "bridge_energy_kj_mol": energy,
        "context_slope_kj_mol": observed_slope,
        "analytic_slope_kj_mol": analytic_slope,
        "algebra_abs_error_kj_mol": error,
        "algebra_allowed_error_kj_mol": allowed,
        "algebra_pass": error <= allowed,
        "finite": all(math.isfinite(value) for value in values),
    }


def _window_output_dir(cell_id: str, xi: float) -> Path:
    return REPLAY_ROOT / cell_id / p1sd.xi_label(xi)


def _cell_result_path(cell_id: str) -> Path:
    return REPLAY_ROOT / cell_id / "cell_result.json"


def _cell_log_path(cell_id: str) -> Path:
    return REPLAY_ROOT / f"{cell_id}.worker.log"


def _write_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%d %H:%M:%S %z")
    with path.open("a", encoding="ascii") as handle:
        handle.write(f"{stamp} {message}\n")


def _window_summary(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    deltas = np.asarray([float(row["delta_raw_kj_mol"]) for row in rows])
    classes = [str(row["classification"]) for row in rows]
    directional = [str(row["directional_cap"]) for row in rows]
    quantile_values = np.quantile(deltas, [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0])
    return {
        "n_frames": len(rows),
        "delta_raw_kj_mol": {
            "min": float(quantile_values[0]),
            "q05": float(quantile_values[1]),
            "q25": float(quantile_values[2]),
            "median": float(quantile_values[3]),
            "q75": float(quantile_values[4]),
            "q95": float(quantile_values[5]),
            "max": float(quantile_values[6]),
            "mean": float(np.mean(deltas)),
            "std": float(np.std(deltas, ddof=1)) if len(deltas) > 1 else 0.0,
        },
        "classification_counts": {
            name: classes.count(name)
            for name in ("linear", "ambiguous", "cap_active")
        },
        "directional_cap_counts": {
            name: directional.count(name) for name in ("dplus", "dminus", "none")
        },
        "max_algebra_abs_error_kj_mol": max(
            float(row["algebra_abs_error_kj_mol"]) for row in rows
        ),
        "original_nonzero_slope_count": sum(
            bool(row["original_slope_nonzero"]) for row in rows
        ),
    }


def run_cell_worker(cell_id: str) -> dict[str, Any]:
    _require_atm_python()
    protocol = load_protocol_and_verify_freeze()
    inventory = load_inventory()
    load_keeper(inventory)
    cell = _cell_row(inventory, cell_id)
    for record in cell["source"].values():
        _verify_artifact(record)
    for window in cell["windows"]:
        for record in window["artifacts"].values():
            _verify_artifact(record)
    output_root = REPLAY_ROOT / cell_id
    if output_root.exists():
        raise FileExistsError(output_root)
    output_root.mkdir(parents=True)
    log_path = _cell_log_path(cell_id)
    started = time.monotonic()
    _write_log(log_path, f"START {cell_id} platform=CUDA md_steps=0")

    import openmm as mm
    from mdtraj.formats import DCDTrajectoryFile
    from openmm import unit

    source_system_path = REPO_ROOT / cell["source"]["bridge_system"]["path"]
    system = mm.XmlSerializer.deserialize(
        source_system_path.read_text(encoding="utf-8")
    )
    if system.getNumParticles() != cell["particle_count"]:
        raise ValueError(f"{cell_id}: System particle count drifted")
    _atm_index, atm_force = dg.find_atm_force(system)
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    platform = mm.Platform.getPlatformByName("CUDA")
    context = mm.Context(
        system,
        integrator,
        platform,
        {
            "DeviceIndex": str(protocol["evaluation"]["device_index"]),
            "Precision": str(protocol["evaluation"]["precision"]),
        },
    )
    p0v2._set_source_parameters(
        context, cell["state"], protocol["evaluation"]["ghost_g"]
    )
    context.setParameter(
        bridge.BRIDGE_PARAMETER,
        float(protocol["evaluation"]["bridge_xi_for_derivative"]),
    )
    softcore = {
        "uoffset_kj_mol": float(context.getParameter("UOffset")),
        "ubcore_kj_mol": float(context.getParameter("Ubcore")),
        "umax_kj_mol": float(context.getParameter("Umax")),
        "acore": float(context.getParameter("Acore")),
        "classification_guard_kj_mol": float(
            protocol["softcore"]["classification_guard_kj_mol"]
        ),
    }
    for name in ("uoffset_kj_mol", "ubcore_kj_mol", "umax_kj_mol", "acore"):
        if not math.isclose(
            softcore[name],
            float(protocol["softcore"][name]),
            abs_tol=1e-9,
            rel_tol=0.0,
        ):
            raise ValueError(
                f"{cell_id}: Context parameter {name} drifted "
                f"{softcore[name]} != {protocol['softcore'][name]}"
            )
    _write_log(log_path, f"CONTEXT_READY {cell_id} particles={cell['particle_count']}")

    window_results = []
    for window in cell["windows"]:
        xi = float(window["xi"])
        _write_log(log_path, f"WINDOW_START xi={xi:.2f}")
        samples_path = REPO_ROOT / window["artifacts"]["samples"]["path"]
        trajectory_path = REPO_ROOT / window["artifacts"]["trajectory"]["path"]
        final_path = REPO_ROOT / window["artifacts"]["final_positions"]["path"]
        source_rows = _read_csv_rows(samples_path)
        frame_rows: list[dict[str, Any]] = []
        fixed_box = True
        last_lengths_nm = None
        last_angles = None
        target_box = np.asarray(cell["parent_box"]["vectors_nm"], dtype=float)
        dcd_to_nm = float(protocol["evaluation"]["dcd_to_nm_scale"])
        with DCDTrajectoryFile(str(trajectory_path), mode="r") as reader:
            if len(reader) != len(source_rows):
                raise ValueError(f"{cell_id} xi={xi:.2f}: DCD/CSV count mismatch")
            if (
                str(reader.distance_unit)
                != protocol["evaluation"]["dcd_distance_unit"]
                or str(reader.distance_unit) != window["dcd_distance_unit"]
            ):
                raise ValueError(
                    f"{cell_id} xi={xi:.2f}: DCD distance unit drifted "
                    f"{reader.distance_unit!s}"
                )
            for frame_index, source_row in enumerate(source_rows, start=1):
                xyz, lengths, angles = reader.read(n_frames=1)
                if xyz.shape != (1, cell["particle_count"], 3):
                    raise ValueError(
                        f"{cell_id} xi={xi:.2f}: invalid DCD frame shape {xyz.shape}"
                    )
                xyz_nm = np.asarray(xyz[0], dtype=float) * dcd_to_nm
                lengths_nm = np.asarray(lengths[0], dtype=float) * dcd_to_nm
                box = _set_frame(context, xyz_nm, lengths_nm, angles[0])
                fixed_box = fixed_box and bool(
                    np.allclose(
                        box,
                        target_box,
                        atol=protocol["validation"][
                            "dcd_box_absolute_tolerance_nm"
                        ],
                        rtol=0.0,
                    )
                )
                evaluated = _evaluate_frame(
                    context,
                    atm_force,
                    softcore,
                    protocol["validation"],
                )
                original_slope = float(source_row["bridge_slope_kj_mol"])
                original_nonzero = (
                    abs(original_slope)
                    > protocol["validation"]["recorded_slope_nonzero_threshold_kj_mol"]
                )
                expected_original_nonzero = (
                    evaluated["classification"] == "cap_active"
                )
                source_consistency = (
                    "ambiguous"
                    if evaluated["classification"] == "ambiguous"
                    else (
                        "match"
                        if original_nonzero == expected_original_nonzero
                        else "mismatch"
                    )
                )
                frame_rows.append(
                    {
                        "cell_id": cell_id,
                        "leg": cell["leg"],
                        "xi": f"{xi:.2f}",
                        "frame": frame_index,
                        "cycle": int(source_row["cycle"]),
                        "step": int(source_row["step"]),
                        **evaluated,
                        "original_slope_kj_mol": original_slope,
                        "original_slope_nonzero": original_nonzero,
                        "original_slope_replay_abs_delta_kj_mol": abs(
                            original_slope - evaluated["context_slope_kj_mol"]
                        ),
                        "original_slope_class_consistency": source_consistency,
                    }
                )
                last_lengths_nm = lengths_nm
                last_angles = np.asarray(angles[0], dtype=float)
        if last_lengths_nm is None or last_angles is None:
            raise ValueError(f"{cell_id} xi={xi:.2f}: empty trajectory")

        final_positions = np.load(final_path, allow_pickle=False)
        _set_frame(context, final_positions, last_lengths_nm, last_angles)
        final_exact = _evaluate_frame(
            context, atm_force, softcore, protocol["validation"]
        )
        final_dcd = frame_rows[-1]
        final_class_agreement = strict_classification_agrees(
            str(final_dcd["classification"]), str(final_exact["classification"])
        )
        checks = {
            "rows_exact": len(frame_rows)
            == protocol["expected"]["frames_per_window"],
            "finite": all(bool(row["finite"]) for row in frame_rows)
            and bool(final_exact["finite"]),
            "algebra": all(bool(row["algebra_pass"]) for row in frame_rows)
            and bool(final_exact["algebra_pass"]),
            "fixed_parent_box": fixed_box,
            "final_float64_classification": final_class_agreement,
            "source_cycles_ordered": [int(row["cycle"]) for row in frame_rows]
            == list(range(1, protocol["expected"]["frames_per_window"] + 1)),
        }
        output_dir = _window_output_dir(cell_id, xi)
        frames_path = output_dir / "frames.csv"
        _write_csv(frames_path, frame_rows)
        result = {
            "schema": "updd_raw_gap_cap_activation_p1cap_window_v1_r1",
            "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "status": "P1CAP_WINDOW_PASS" if all(checks.values()) else "P1CAP_WINDOW_FAIL",
            "cell_id": cell_id,
            "seed": cell["seed"],
            "leg": cell["leg"],
            "xi": xi,
            "inventory_digest": inventory["inventory_digest"],
            "platform": {
                "name": "CUDA",
                "openmm_version": mm.__version__,
                "properties": {
                    "DeviceIndex": str(protocol["evaluation"]["device_index"]),
                    "Precision": str(protocol["evaluation"]["precision"]),
                },
            },
            "md_steps": 0,
            "softcore": softcore,
            "summary": _window_summary(frame_rows),
            "final_float64_check": {
                "dcd_delta_raw_kj_mol": final_dcd["delta_raw_kj_mol"],
                "float64_delta_raw_kj_mol": final_exact["delta_raw_kj_mol"],
                "abs_delta_kj_mol": abs(
                    final_dcd["delta_raw_kj_mol"]
                    - final_exact["delta_raw_kj_mol"]
                ),
                "dcd_classification": final_dcd["classification"],
                "float64_classification": final_exact["classification"],
                "strict_classification_agreement": final_class_agreement,
                "float64_algebra_pass": final_exact["algebra_pass"],
            },
            "checks": checks,
            "artifacts": {"frames": artifact_record(frames_path)},
            "claims": protocol["claims"],
        }
        result_path = output_dir / "window_result.json"
        atomic_write_json(result_path, result)
        result["artifacts"]["result"] = artifact_record(result_path)
        window_results.append(result)
        _write_log(log_path, f"WINDOW_{'PASS' if all(checks.values()) else 'FAIL'} xi={xi:.2f}")
        del final_positions

    del context, integrator, atm_force, system
    gc.collect()
    checks = {
        "windows_exact": len(window_results) == len(cell["windows"]),
        "all_windows_pass": all(
            row["status"] == "P1CAP_WINDOW_PASS" for row in window_results
        ),
        "frames_exact": sum(row["summary"]["n_frames"] for row in window_results)
        == len(cell["windows"]) * protocol["expected"]["frames_per_window"],
        "one_context": True,
        "md_steps_zero": True,
    }
    payload = {
        "schema": "updd_raw_gap_cap_activation_p1cap_cell_v1_r1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "P1CAP_CELL_PASS" if all(checks.values()) else "P1CAP_CELL_FAIL",
        "cell_id": cell_id,
        "seed": cell["seed"],
        "leg": cell["leg"],
        "inventory_digest": inventory["inventory_digest"],
        "context_count": 1,
        "md_steps": 0,
        "elapsed_s": time.monotonic() - started,
        "softcore": softcore,
        "checks": checks,
        "windows": [
            {
                "xi": row["xi"],
                "status": row["status"],
                "summary": row["summary"],
                "final_float64_check": row["final_float64_check"],
                "checks": row["checks"],
                "artifacts": row["artifacts"],
            }
            for row in window_results
        ],
        "claims": protocol["claims"],
    }
    atomic_write_json(_cell_result_path(cell_id), payload)
    _write_log(log_path, payload["status"])
    if payload["status"] != "P1CAP_CELL_PASS":
        raise RuntimeError(f"{cell_id}: P1CAP cell failed")
    return payload


def _numeric_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    bool_fields = {"algebra_pass", "finite", "original_slope_nonzero"}
    int_fields = {"frame", "cycle", "step"}
    string_fields = {
        "cell_id",
        "leg",
        "classification",
        "directional_cap",
        "original_slope_class_consistency",
    }
    for path in sorted(REPLAY_ROOT.glob("*/xi_*/frames.csv")):
        for raw in _read_csv_rows(path):
            parsed: dict[str, Any] = {}
            for key, value in raw.items():
                if key in string_fields:
                    parsed[key] = value
                elif key in bool_fields:
                    parsed[key] = value == "True"
                elif key in int_fields:
                    parsed[key] = int(value)
                elif key == "xi":
                    parsed[key] = float(value)
                else:
                    parsed[key] = float(value)
            rows.append(parsed)
    return rows


def summarize() -> dict[str, Any]:
    protocol = load_protocol_and_verify_freeze()
    inventory = load_inventory()
    load_keeper(inventory)
    cell_results = []
    for spec in _cohort(protocol):
        path = _cell_result_path(spec["cell_id"])
        if not path.is_file():
            raise FileNotFoundError(path)
        result = json.loads(path.read_text(encoding="utf-8"))
        if result.get("inventory_digest") != inventory["inventory_digest"]:
            raise ValueError(f"{spec['cell_id']}: P1CAP cell inventory drifted")
        cell_results.append(result)
    rows = _numeric_rows()
    deltas = np.asarray([row["delta_raw_kj_mol"] for row in rows], dtype=float)
    classes = [row["classification"] for row in rows]
    directions = [row["directional_cap"] for row in rows]
    quantiles = np.quantile(
        deltas, [0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0]
    )
    contingency: dict[str, int] = {}
    for row in rows:
        original = "nonzero" if row["original_slope_nonzero"] else "zero"
        key = f"{row['classification']}|original_{original}"
        contingency[key] = contingency.get(key, 0) + 1
    top = sorted(rows, key=lambda row: abs(row["delta_raw_kj_mol"]), reverse=True)[:20]
    all_window_results = [
        {
            "cell_id": cell["cell_id"],
            "leg": cell["leg"],
            **window,
        }
        for cell in cell_results
        for window in cell["windows"]
    ]
    checks = {
        "cells_exact": len(cell_results) == protocol["expected"]["cells"],
        "windows_exact": len(all_window_results) == protocol["expected"]["windows"],
        "frames_exact": len(rows) == protocol["expected"]["frames"],
        "all_cells_pass": all(
            row["status"] == "P1CAP_CELL_PASS" for row in cell_results
        ),
        "all_windows_pass": all(
            row["status"] == "P1CAP_WINDOW_PASS" for row in all_window_results
        ),
        "all_finite": all(row["finite"] for row in rows),
        "all_algebra_pass": all(row["algebra_pass"] for row in rows),
        "all_final_classifications_pass": all(
            row["final_float64_check"]["strict_classification_agreement"]
            for row in all_window_results
        ),
        "md_steps_zero": all(row["md_steps"] == 0 for row in cell_results),
    }
    payload = {
        "schema": "updd_raw_gap_cap_activation_p1cap_runner_summary_v1_r1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "P1CAP_RUNNER_PASS" if all(checks.values()) else "P1CAP_RUNNER_FAIL",
        "inventory_digest": inventory["inventory_digest"],
        "keeper_digest": json.loads(
            KEEPER_PREFLIGHT.read_text(encoding="utf-8")
        )["keeper_digest"],
        "counts": {
            "cells": len(cell_results),
            "windows": len(all_window_results),
            "frames": len(rows),
            "linear": classes.count("linear"),
            "ambiguous": classes.count("ambiguous"),
            "cap_active": classes.count("cap_active"),
            "dplus_active": directions.count("dplus"),
            "dminus_active": directions.count("dminus"),
            "original_nonzero_slope": sum(
                row["original_slope_nonzero"] for row in rows
            ),
            "original_class_mismatch": sum(
                row["original_slope_class_consistency"] == "mismatch"
                for row in rows
            ),
        },
        "fractions": {
            "linear": classes.count("linear") / len(rows),
            "ambiguous": classes.count("ambiguous") / len(rows),
            "cap_active": classes.count("cap_active") / len(rows),
        },
        "delta_raw_kj_mol": {
            "min": float(quantiles[0]),
            "q01": float(quantiles[1]),
            "q05": float(quantiles[2]),
            "q25": float(quantiles[3]),
            "median": float(quantiles[4]),
            "q75": float(quantiles[5]),
            "q95": float(quantiles[6]),
            "q99": float(quantiles[7]),
            "max": float(quantiles[8]),
            "mean": float(np.mean(deltas)),
            "std": float(np.std(deltas, ddof=1)),
        },
        "max_algebra_abs_error_kj_mol": max(
            row["algebra_abs_error_kj_mol"] for row in rows
        ),
        "max_original_slope_replay_abs_delta_kj_mol": max(
            row["original_slope_replay_abs_delta_kj_mol"] for row in rows
        ),
        "contingency": contingency,
        "cells": [
            {
                "cell_id": row["cell_id"],
                "leg": row["leg"],
                "status": row["status"],
                "elapsed_s": row["elapsed_s"],
                "windows": row["windows"],
            }
            for row in cell_results
        ],
        "highest_absolute_gap_frames": [
            {
                "cell_id": row["cell_id"],
                "leg": row["leg"],
                "xi": row["xi"],
                "frame": row["frame"],
                "cycle": row["cycle"],
                "step": row["step"],
                "delta_raw_kj_mol": row["delta_raw_kj_mol"],
                "classification": row["classification"],
                "directional_cap": row["directional_cap"],
                "original_slope_kj_mol": row["original_slope_kj_mol"],
                "context_slope_kj_mol": row["context_slope_kj_mol"],
            }
            for row in top
        ],
        "checks": checks,
        "claims": protocol["claims"],
    }
    payload["summary_digest"] = summary_digest(payload)
    atomic_write_json(RUNNER_SUMMARY, payload)

    count = payload["counts"]
    lines = [
        "# P1CAP R1 Runner Summary",
        "",
        f"Status: `{payload['status']}`",
        "",
        "P1CAP replayed the 16 completed rejected P1SD windows with zero MD steps.",
        "No free energy or production-use claim is made.",
        "",
        "## Classification",
        "",
        "| class | frames | fraction |",
        "|---|---:|---:|",
        f"| linear | {count['linear']} | {payload['fractions']['linear']:.6f} |",
        f"| ambiguous | {count['ambiguous']} | {payload['fractions']['ambiguous']:.6f} |",
        f"| cap active | {count['cap_active']} | {payload['fractions']['cap_active']:.6f} |",
        "",
        f"dplus active: `{count['dplus_active']}`; "
        f"dminus active: `{count['dminus_active']}`; "
        f"original nonzero slopes: `{count['original_nonzero_slope']}`.",
        "",
        "## Windows",
        "",
        "| cell | xi | status | linear | ambiguous | active |",
        "|---|---:|---|---:|---:|---:|",
    ]
    for row in all_window_results:
        counts = row["summary"]["classification_counts"]
        lines.append(
            f"| {row['cell_id']} | {row['xi']:.2f} | {row['status']} | "
            f"{counts['linear']} | {counts['ambiguous']} | {counts['cap_active']} |"
        )
    lines.extend(
        [
            "",
            "PATH interpretation is recorded separately after the formal Runner gate.",
            "",
        ]
    )
    atomic_write_text(RUNNER_REPORT, "\n".join(lines))
    return payload


def _write_run_state(
    inventory: dict[str, Any],
    *,
    status: str,
    completed: Sequence[str],
    current: str | None,
    error: str | None = None,
) -> None:
    atomic_write_json(
        RUN_STATE,
        {
            "schema": "updd_raw_gap_cap_activation_p1cap_run_state_v1_r1",
            "updated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "status": status,
            "inventory_digest": inventory["inventory_digest"],
            "completed_cells": list(completed),
            "current_cell": current,
            "error": error,
        },
    )


def run_all() -> dict[str, Any]:
    _require_atm_python()
    protocol = load_protocol_and_verify_freeze()
    inventory = load_inventory()
    load_keeper(inventory)
    if REPLAY_ROOT.exists() or RUN_STATE.exists() or RUNNER_SUMMARY.exists():
        raise RuntimeError("P1CAP output already exists; preserve or archive it")
    completed: list[str] = []
    _write_run_state(
        inventory, status="RUNNING", completed=completed, current=None
    )
    for spec in _cohort(protocol):
        cell_id = spec["cell_id"]
        _write_run_state(
            inventory, status="RUNNING", completed=completed, current=cell_id
        )
        result = subprocess.run(
            [str(EXPECTED_ATM_PYTHON), str(RUNNER_PATH), "--worker", cell_id],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.stdout:
            print(result.stdout, end="")
        if result.returncode != 0:
            error = result.stderr.strip() or f"worker return code {result.returncode}"
            _write_run_state(
                inventory,
                status="RUNNER_FAIL_STOPPED",
                completed=completed,
                current=cell_id,
                error=error,
            )
            raise RuntimeError(f"{cell_id}: P1CAP worker failed\n{error}")
        completed.append(cell_id)
    summary = summarize()
    final_status = (
        "P1CAP_COMPLETE"
        if summary["status"] == "P1CAP_RUNNER_PASS"
        else "RUNNER_FAIL_STOPPED"
    )
    _write_run_state(
        inventory,
        status=final_status,
        completed=completed,
        current=None,
        error=None if final_status == "P1CAP_COMPLETE" else summary["status"],
    )
    if final_status != "P1CAP_COMPLETE":
        raise RuntimeError("P1CAP formal Runner gate failed")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--build-inventory", action="store_true")
    modes.add_argument("--keeper-preflight", action="store_true")
    modes.add_argument("--worker", metavar="CELL_ID")
    modes.add_argument("--run", action="store_true")
    modes.add_argument("--summarize", action="store_true")
    args = parser.parse_args()
    try:
        if args.build_inventory:
            result = build_inventory()
            print(
                f"P1CAP_INVENTORY_PASS digest={result['inventory_digest']} "
                f"frames={result['counts']['frames']}"
            )
        elif args.keeper_preflight:
            result = run_keeper_preflight()
            print(
                f"P1CAP_KEEPER_{result['status']} "
                f"digest={result['keeper_digest']}"
            )
        elif args.worker:
            result = run_cell_worker(args.worker)
            print(
                f"{result['status']} cell={result['cell_id']} "
                f"elapsed_s={result['elapsed_s']:.1f}"
            )
        elif args.run:
            result = run_all()
            print(
                f"{result['status']} frames={result['counts']['frames']} "
                f"active={result['counts']['cap_active']}"
            )
        else:
            result = summarize()
            print(
                f"{result['status']} frames={result['counts']['frames']} "
                f"active={result['counts']['cap_active']}"
            )
        return 0
    except Exception:
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
