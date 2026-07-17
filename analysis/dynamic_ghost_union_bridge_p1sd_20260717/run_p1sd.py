#!/usr/bin/env python3
"""Run the frozen union-solvated apex-bridge P1 schedule-discovery pilot."""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import time
import traceback
from pathlib import Path
from statistics import median
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
from analysis.dynamic_ghost_union_bridge_p0v3r2_20260717 import (  # noqa: E402
    run_p0v3r2_reference_audit as p0,
)
import trackb_apex_bridge as bridge  # noqa: E402
import trackb_dynamic_ghost as dg  # noqa: E402


ANALYSIS_DIR = Path(__file__).resolve().parent
PREREGISTRATION = ANALYSIS_DIR / "PREREGISTRATION.md"
PROTOCOL_PATH = ANALYSIS_DIR / "protocol_r2.json"
FROZEN_MANIFEST = ANALYSIS_DIR / "FROZEN_MANIFEST_R2.json"
RUNNER_PATH = Path(__file__).resolve()

PARENT_DIR = REPO_ROOT / "analysis/dynamic_ghost_union_bridge_p0v3r2_20260717"
PARENT_INVENTORY = PARENT_DIR / "P0V3R2_PREBUILD_INVENTORY.json"
PARENT_SOURCE_SUMMARY = PARENT_DIR / "P0V3R2_SOURCE_SUMMARY.json"
PARENT_REFERENCE_SUMMARY = PARENT_DIR / "P0V3R2_REFERENCE_SUMMARY.json"
PARENT_PATH = PARENT_DIR / "P0V3R2_PATH_DIAGNOSIS.json"

PREBUILD_INVENTORY = ANALYSIS_DIR / "P1SD_PREBUILD_INVENTORY.json"
NORMALIZED_ROOT = ANALYSIS_DIR / "normalized_sources"
NORMALIZED_SUMMARY = ANALYSIS_DIR / "P1SD_NORMALIZED_SOURCE_SUMMARY.json"
KEEPER_PREFLIGHT = ANALYSIS_DIR / "P1SD_KEEPER_PREFLIGHT.json"
REFERENCE_ROOT = ANALYSIS_DIR / "reference_preflight"
REFERENCE_SUMMARY = ANALYSIS_DIR / "P1SD_REFERENCE_SUMMARY.json"
SAMPLE_ROOT = ANALYSIS_DIR / "sampled_schedule_discovery"
SAMPLE_SUMMARY = ANALYSIS_DIR / "P1SD_SCHEDULE_SUMMARY.json"
SAMPLE_REPORT = ANALYSIS_DIR / "P1SD_SCHEDULE_SUMMARY.md"
RUN_STATE = ANALYSIS_DIR / "P1SD_RUN_STATE.json"

EXPECTED_ATM_PYTHON = Path("/home/san/miniconda3/envs/atm/bin/python")
SOLVENT_NAMES = frozenset({"HOH", "WAT", "NA", "CL", "Na+", "Cl-"})
WATER_NAMES = frozenset({"HOH", "WAT"})

relative = p0v2.relative
sha256_file = p0v2.sha256_file
sha256_text = p0v2.sha256_text
canonical_digest = p0v2.canonical_digest
artifact_record = p0v2.artifact_record
atomic_write_json = p0v2.atomic_write_json
atomic_write_text = p0v2.atomic_write_text


def _hash_map(base: Path, mapping: dict[str, str], label: str) -> None:
    for name, expected in mapping.items():
        path = base / name
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(f"{label} drift: {path} {observed} != {expected}")


def load_protocol_and_verify_freeze() -> dict[str, Any]:
    frozen = json.loads(FROZEN_MANIFEST.read_text(encoding="utf-8"))
    expected_status = (
        "FROZEN_AFTER_ENDPOINT_SLOPE_REVIEW_BEFORE_P1SD_OFFICIAL_OUTPUT"
    )
    if frozen.get("status") != expected_status:
        raise ValueError("P1SD frozen manifest status drifted")
    _hash_map(ANALYSIS_DIR, frozen.get("files") or {}, "P1SD frozen file")
    _hash_map(REPO_ROOT, frozen.get("parent_evidence") or {}, "P1SD parent evidence")
    _hash_map(REPO_ROOT, frozen.get("trigger_evidence") or {}, "P1SD R1 trigger evidence")

    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    if protocol.get("status") != expected_status:
        raise ValueError("P1SD protocol status drifted")
    if protocol.get("cohort") != [cell.cell_id for cell in p0.make_cells()]:
        raise ValueError("P1SD cohort or order drifted")
    sampling = protocol.get("sampling") or {}
    if sampling.get("xi") != [0.0, 0.25, 0.5, 0.75, 1.0]:
        raise ValueError("P1SD initial bridge grid drifted")
    if sampling.get("barostat") is not False:
        raise ValueError("P1SD must remain fixed-volume NVT")
    if sampling.get("one_live_context") is not True:
        raise ValueError("P1SD one-live-Context contract drifted")
    if sampling.get("final_step") != 14500:
        raise ValueError("P1SD step contract drifted")
    if (protocol.get("mbar") or {}).get("pilot_samples_reusable_in_production") is not False:
        raise ValueError("P1SD pilot-sample reuse prohibition drifted")
    readback = (protocol.get("cross_state_matrix") or {}).get(
        "explicit_readback_tolerance"
    ) or {}
    if readback != {
        "absolute_floor_kj_mol": 1e-5,
        "relative_to_energy_scale": 1e-9,
        "absolute_ceiling_kj_mol": 1e-3,
        "formula": "min(ceiling,max(floor,relative*max(1,abs(observed),abs(expected))))",
    }:
        raise ValueError("P1SD R1 readback tolerance drifted")
    slope = (protocol.get("cross_state_matrix") or {}).get("slope_evaluation")
    if slope != {
        "parameter": bridge.BRIDGE_PARAMETER,
        "probe_xi": 0.5,
        "every_saved_sample": True,
        "restore_generating_xi_before_next_step": True,
        "generating_xi_derivative_is_diagnostic_only": True,
    }:
        raise ValueError("P1SD R2 slope-evaluation contract drifted")

    parent_inventory = json.loads(PARENT_INVENTORY.read_text(encoding="utf-8"))
    parent_sources = json.loads(PARENT_SOURCE_SUMMARY.read_text(encoding="utf-8"))
    parent_reference = json.loads(PARENT_REFERENCE_SUMMARY.read_text(encoding="utf-8"))
    parent_path = json.loads(PARENT_PATH.read_text(encoding="utf-8"))
    declared = protocol.get("parent") or {}
    if parent_inventory.get("inventory_digest") != declared.get("inventory_digest"):
        raise ValueError("P1SD parent inventory digest drifted")
    if parent_sources.get("source_summary_digest") != declared.get("source_summary_digest"):
        raise ValueError("P1SD parent source digest drifted")
    if parent_reference.get("status") != declared.get("reference_status"):
        raise ValueError("P1SD parent Reference status drifted")
    if parent_path.get("path_verdict") != declared.get("path_verdict"):
        raise ValueError("P1SD parent PATH verdict drifted")
    if parent_sources.get("status") != "PASS" or parent_sources.get("n_completed") != 6:
        raise ValueError("P1SD parent source summary is not a complete PASS")
    if parent_reference.get("status") != "P0V3R2_PASS":
        raise ValueError("P1SD parent Reference summary is not P0V3R2_PASS")
    return protocol


def seed_unit_key(cell_id: str, xi: float, role: str) -> str:
    if role not in {"integrator", "velocity"}:
        raise ValueError(f"invalid P1SD seed role {role!r}")
    return f"P1SD|{cell_id}|xi={float(xi):.2f}|{role}"


def derive_seed(cell_id: str, xi: float, role: str) -> int:
    key = seed_unit_key(cell_id, xi, role)
    value = int(hashlib.sha256(key.encode("ascii")).hexdigest()[:8], 16) & 0x7FFFFFFF
    return value or 1


def inventory_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("inventory_digest", None)
    return canonical_digest(canonical)


def normalized_summary_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("normalized_summary_digest", None)
    return canonical_digest(canonical)


def _implementation_records() -> dict[str, dict[str, Any]]:
    return {
        "p1sd_runner": artifact_record(RUNNER_PATH),
        "p0v3r2_runner": artifact_record(p0.RUNNER_PATH),
        "setup_module": artifact_record(p0.SETUP_IMPLEMENTATION),
        "serializer_module": artifact_record(p0.SERIALIZER_IMPLEMENTATION),
        "ghost_module": artifact_record(p0.GHOST_IMPLEMENTATION),
        "bridge_module": artifact_record(p0.BRIDGE_IMPLEMENTATION),
    }


def _parent_rows() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    inventory = json.loads(PARENT_INVENTORY.read_text(encoding="utf-8"))
    source_summary = json.loads(PARENT_SOURCE_SUMMARY.read_text(encoding="utf-8"))
    reference_summary = json.loads(PARENT_REFERENCE_SUMMARY.read_text(encoding="utf-8"))
    return inventory, source_summary, reference_summary


def _parent_cell_maps() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    inventory, source_summary, reference_summary = _parent_rows()
    inventory_map = {row["cell_id"]: row for row in inventory["cells"]}
    source_map = {row["cell"]["cell_id"]: row for row in source_summary["cells"]}
    reference_map = {row["cell_id"]: row for row in reference_summary["cells"]}
    expected = {cell.cell_id for cell in p0.make_cells()}
    if set(inventory_map) != expected or set(source_map) != expected or set(reference_map) != expected:
        raise ValueError("P1SD parent six-cell maps drifted")
    return inventory_map, source_map, reference_map


def build_prebuild_inventory() -> dict[str, Any]:
    protocol = load_protocol_and_verify_freeze()
    if Path(sys.executable).resolve() != EXPECTED_ATM_PYTHON.resolve():
        raise ValueError(f"P1SD inventory requires {EXPECTED_ATM_PYTHON}, got {sys.executable}")
    import openmm
    import pymbar
    import scipy

    parent_inventory, parent_sources, parent_reference = _parent_cell_maps()
    cells = []
    for cell in p0.make_cells():
        source = parent_sources[cell.cell_id]
        reference = parent_reference[cell.cell_id]
        frozen = parent_inventory[cell.cell_id]
        if source.get("status") != "PASS" or reference.get("status") != "P0V3R2_PASS":
            raise ValueError(f"{cell.cell_id}: parent cell is not a complete PASS")
        source_artifacts = source.get("artifacts") or {}
        bridge_artifact = reference.get("bridge_system")
        for record in list(source_artifacts.values()) + [bridge_artifact]:
            path = REPO_ROOT / record["path"]
            if artifact_record(path) != record:
                raise ValueError(f"{cell.cell_id}: parent artifact missing or drifted: {path}")
        cells.append(
            {
                "cell_id": cell.cell_id,
                "seed": cell.seed,
                "leg": cell.leg,
                "replicate_index": cell.replicate_index,
                "parent_box": frozen["parent_contract"]["box"],
                "parent_contract": frozen["parent_contract"],
                "states": frozen["states"],
                "source_result": source["result"],
                "source_artifacts": source_artifacts,
                "bridge_artifact": bridge_artifact,
                "parent_cell_digest": frozen["cell_digest"],
            }
        )
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_inventory_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "FROZEN_BEFORE_P1SD_OUTPUT",
        "python": sys.executable,
        "runtime": {
            "openmm": openmm.__version__,
            "pymbar": pymbar.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "frozen_manifest_sha256": sha256_file(FROZEN_MANIFEST),
        "implementation": _implementation_records(),
        "parent_inventory_digest": protocol["parent"]["inventory_digest"],
        "parent_source_summary_digest": protocol["parent"]["source_summary_digest"],
        "xi": protocol["sampling"]["xi"],
        "cells": cells,
    }
    payload["inventory_digest"] = inventory_digest(payload)
    atomic_write_json(PREBUILD_INVENTORY, payload)
    return payload


def load_inventory() -> dict[str, Any]:
    load_protocol_and_verify_freeze()
    if not PREBUILD_INVENTORY.is_file():
        raise FileNotFoundError(f"missing {PREBUILD_INVENTORY}; run --build-inventory")
    inventory = json.loads(PREBUILD_INVENTORY.read_text(encoding="utf-8"))
    if inventory_digest(inventory) != inventory.get("inventory_digest"):
        raise ValueError("P1SD inventory digest drifted")
    if inventory.get("implementation") != _implementation_records():
        raise ValueError("P1SD implementation drift after inventory freeze")
    if inventory.get("protocol_sha256") != sha256_file(PROTOCOL_PATH):
        raise ValueError("P1SD protocol drift after inventory freeze")
    if inventory.get("frozen_manifest_sha256") != sha256_file(FROZEN_MANIFEST):
        raise ValueError("P1SD frozen manifest drift after inventory freeze")
    return inventory


def frozen_cell(inventory: dict[str, Any], cell_id: str) -> dict[str, Any]:
    rows = [row for row in inventory["cells"] if row["cell_id"] == cell_id]
    if len(rows) != 1:
        raise ValueError(f"expected one P1SD inventory row for {cell_id}, found {len(rows)}")
    return rows[0]


def normalized_cell_dir(cell_id: str) -> Path:
    return NORMALIZED_ROOT / cell_id


def normalized_result_path(cell_id: str) -> Path:
    return normalized_cell_dir(cell_id) / "normalization_result.json"


def normalized_positions_path(cell_id: str) -> Path:
    return normalized_cell_dir(cell_id) / "positions_float64.npy"


def normalized_pdb_path(cell_id: str) -> Path:
    return normalized_cell_dir(cell_id) / "normalized_source.pdb"


def normalized_source_xml_path(cell_id: str) -> Path:
    return normalized_cell_dir(cell_id) / "normalized_source_system.xml"


def normalized_source_ghost_xml_path(cell_id: str) -> Path:
    return normalized_cell_dir(cell_id) / "normalized_source_ghost_system.xml"


def normalized_bridge_xml_path(cell_id: str) -> Path:
    return normalized_cell_dir(cell_id) / "normalized_bridge_system.xml"


def _atomic_save_npy(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp.npy")
    np.save(temp, np.asarray(array, dtype=np.float64), allow_pickle=False)
    os.replace(temp, path)


def _atomic_write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_bytes(payload)
    os.replace(temp, path)


def remap_rigid_residue_anchors(
    positions_nm: np.ndarray,
    old_box_nm: np.ndarray,
    new_box_nm: np.ndarray,
    residue_groups: Sequence[tuple[int, Sequence[int]]],
) -> np.ndarray:
    """Map residue anchors fractionally while preserving each residue internally."""
    positions = np.asarray(positions_nm, dtype=np.float64)
    old_box = np.asarray(old_box_nm, dtype=np.float64)
    new_box = np.asarray(new_box_nm, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must have shape [N,3]")
    if old_box.shape != (3, 3) or new_box.shape != (3, 3):
        raise ValueError("box matrices must have shape [3,3]")
    if np.linalg.det(old_box) <= 0.0 or np.linalg.det(new_box) <= 0.0:
        raise ValueError("box matrices must be right-handed")
    output = positions.copy()
    inverse_old = np.linalg.inv(old_box)
    seen: set[int] = set()
    for anchor, indices_raw in residue_groups:
        indices = [int(index) for index in indices_raw]
        if not indices or int(anchor) not in indices:
            raise ValueError("every remapped residue needs its anchor in the index set")
        if any(index in seen for index in indices):
            raise ValueError("remapped residue groups overlap")
        seen.update(indices)
        fractional = positions[int(anchor)] @ inverse_old
        fractional -= np.floor(fractional)
        target = fractional @ new_box
        output[indices] += target - positions[int(anchor)]
    return output


def _box_matrix_nm(topology: Any) -> np.ndarray:
    from openmm import unit

    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        raise ValueError("topology has no periodic box")
    return np.asarray([v.value_in_unit(unit.nanometer) for v in vectors], dtype=float)


def _set_box(topology: Any, system: Any, box_nm: np.ndarray) -> None:
    from openmm import Vec3, unit

    vectors = tuple(Vec3(*row) for row in np.asarray(box_nm, dtype=float))
    topology.setPeriodicBoxVectors(vectors * unit.nanometer)
    system.setDefaultPeriodicBoxVectors(*(vector * unit.nanometer for vector in vectors))


def _force_xml_digests(system: Any) -> list[str]:
    import openmm as mm

    return [
        sha256_text(mm.XmlSerializer.serialize(system.getForce(index)))
        for index in range(system.getNumForces())
    ]


def _residue_groups_and_indices(topology: Any) -> tuple[list[tuple[int, list[int]]], list[int], list[int]]:
    groups: list[tuple[int, list[int]]] = []
    solute_indices: list[int] = []
    water_oxygen_indices: list[int] = []
    for residue in topology.residues():
        atoms = list(residue.atoms())
        indices = [atom.index for atom in atoms]
        if residue.name in WATER_NAMES:
            oxygens = [
                atom.index
                for atom in atoms
                if atom.element is not None and atom.element.symbol == "O"
            ]
            if len(oxygens) != 1:
                raise ValueError(f"water residue {residue} has {len(oxygens)} oxygens")
            groups.append((oxygens[0], indices))
            water_oxygen_indices.append(oxygens[0])
        elif residue.name in SOLVENT_NAMES:
            if len(atoms) != 1:
                raise ValueError(f"ion residue {residue} is not monatomic")
            groups.append((atoms[0].index, indices))
        else:
            solute_indices.extend(indices)
    return groups, solute_indices, water_oxygen_indices


def _max_internal_vector_delta(
    before: np.ndarray,
    after: np.ndarray,
    groups: Iterable[tuple[int, Sequence[int]]],
) -> float:
    maximum = 0.0
    for anchor, indices in groups:
        idx = np.asarray(list(indices), dtype=int)
        before_rel = before[idx] - before[int(anchor)]
        after_rel = after[idx] - after[int(anchor)]
        maximum = max(maximum, float(np.max(np.abs(before_rel - after_rel))))
    return maximum


def _periodic_minimum_distance(
    query_nm: np.ndarray, target_nm: np.ndarray, box_nm: np.ndarray
) -> float:
    from scipy.spatial import cKDTree

    box = np.asarray(box_nm, dtype=float)
    if not np.allclose(box, np.diag(np.diag(box)), atol=1e-12, rtol=0.0):
        raise ValueError("P1SD periodic KD-tree gate requires an orthogonal box")
    lengths = np.diag(box)
    query = np.mod(np.asarray(query_nm, dtype=float), lengths)
    target = np.mod(np.asarray(target_nm, dtype=float), lengths)
    distances, _indices = cKDTree(target, boxsize=lengths).query(query, k=1, workers=-1)
    return float(np.min(distances))


def _write_inspection_pdb(path: Path, topology: Any, positions_nm: np.ndarray) -> None:
    from openmm import unit
    from openmm.app import PDBFile

    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="ascii") as handle:
        PDBFile.writeFile(topology, positions_nm * unit.nanometer, handle, keepIds=True)
    os.replace(temp, path)


def _validate_parent_artifacts(cell_record: dict[str, Any]) -> None:
    records = list((cell_record.get("source_artifacts") or {}).values())
    records.append(cell_record["bridge_artifact"])
    records.append(cell_record["source_result"])
    for record in records:
        path = REPO_ROOT / record["path"]
        if artifact_record(path) != record:
            raise ValueError(f"{cell_record['cell_id']}: parent artifact drift: {path}")


def normalize_source_worker(cell_id: str, inventory: dict[str, Any]) -> dict[str, Any]:
    import openmm as mm
    from openmm import unit
    from openmm.app import PDBFile

    protocol, _base, bridge_protocol, ghost_protocol = p0.load_protocol_and_verify_freeze()
    del protocol, _base
    frozen = frozen_cell(inventory, cell_id)
    _validate_parent_artifacts(frozen)
    output_dir = normalized_cell_dir(cell_id)
    if output_dir.exists():
        raise RuntimeError(f"{cell_id}: normalized output already exists")
    output_dir.mkdir(parents=True)
    started = time.monotonic()

    source_xml_path = REPO_ROOT / frozen["source_artifacts"]["serialized_xml"]["path"]
    source_pdb_path = REPO_ROOT / frozen["source_artifacts"]["serialized_pdb"]["path"]
    source_xml = source_xml_path.read_text(encoding="utf-8")
    system = mm.XmlSerializer.deserialize(source_xml)
    pdb = PDBFile(str(source_pdb_path))
    before = np.asarray(pdb.positions.value_in_unit(unit.nanometer), dtype=np.float64)
    old_box = _box_matrix_nm(pdb.topology)
    target_box = np.asarray(frozen["parent_box"]["vectors_nm"], dtype=float)
    groups, solute_indices, water_oxygen_indices = _residue_groups_and_indices(pdb.topology)
    after = remap_rigid_residue_anchors(before, old_box, target_box, groups)

    solute_delta = float(np.max(np.abs(after[solute_indices] - before[solute_indices])))
    internal_delta = _max_internal_vector_delta(before, after, groups)
    _set_box(pdb.topology, system, target_box)
    source_serialized = mm.XmlSerializer.serialize(system)

    source_ghost = mm.XmlSerializer.deserialize(source_serialized)
    ghost_report = dg.build_dynamic_ghost_force(source_ghost, pdb.topology, ghost_protocol)
    source_ghost_serialized = mm.XmlSerializer.serialize(source_ghost)

    bridged = mm.XmlSerializer.deserialize(source_ghost_serialized)
    bridge_report = bridge.build_apex_bridge(bridged, bridge_protocol)
    bridge_serialized = mm.XmlSerializer.serialize(bridged)

    ring_geometry = p0._ring_water_geometry(
        system, pdb.topology, after * unit.nanometer
    )
    atoms = list(pdb.topology.atoms())
    solute_heavy = [
        index
        for index in solute_indices
        if atoms[index].element is not None and atoms[index].element.symbol != "H"
    ]
    all_solute_water_min = _periodic_minimum_distance(
        after[solute_heavy], after[water_oxygen_indices], target_box
    )
    source_counts = p0.ats._solvent_counts(pdb.topology)
    source_charge = p0._net_charge_e(system)
    normalization = load_protocol_and_verify_freeze()["parent_box_normalization"]
    tolerances = normalization["tolerances"]
    ghost_contract_valid = bool(
        ghost_report["n_particles_before"] == ghost_report["n_particles_after"]
        and ghost_report["n_forces_after"] == ghost_report["n_forces_before"] + 1
        and ghost_report["force_index"] > ghost_report["atm_force_index"]
        and len(ghost_report["selection"]["ring_indices"]) == 9
        and ghost_report["selection"]["n_water_oxygens"] == source_counts["water"]
        and ghost_report["energy_parameter_derivatives"]
        == [dg.GHOST_GLOBAL_PARAMETER]
    )
    checks = {
        "particle_count_exact": system.getNumParticles() == frozen["parent_contract"]["particles"],
        "topology_count_exact": pdb.topology.getNumAtoms() == frozen["parent_contract"]["topology_atoms"],
        "solvent_counts_exact": source_counts == frozen["parent_contract"]["solvent_counts"],
        "net_charge_exact": abs(source_charge - frozen["parent_contract"]["net_charge_e"]) <= 1e-8,
        "solute_positions_exact": solute_delta <= tolerances["solute_position_max_abs_nm"],
        "solvent_internal_vectors_exact": internal_delta <= tolerances["solvent_internal_vector_max_abs_nm"],
        "box_exact": float(np.max(np.abs(_box_matrix_nm(pdb.topology) - target_box))) <= tolerances["box_vector_max_abs_nm"],
        "ring_geometry": ring_geometry["passed"],
        "all_solute_water_geometry": all_solute_water_min >= tolerances["all_solute_heavy_water_o_min_nm"],
        "ghost_contract": ghost_contract_valid,
        "bridge_contract": all(bridge_report["checks"].values()),
        "source_force_xml_unchanged_except_box": _force_xml_digests(system)
        == _force_xml_digests(mm.XmlSerializer.deserialize(source_xml)),
    }

    _atomic_save_npy(normalized_positions_path(cell_id), after)
    _write_inspection_pdb(normalized_pdb_path(cell_id), pdb.topology, after)
    atomic_write_text(normalized_source_xml_path(cell_id), source_serialized)
    atomic_write_text(normalized_source_ghost_xml_path(cell_id), source_ghost_serialized)
    atomic_write_text(normalized_bridge_xml_path(cell_id), bridge_serialized)
    artifacts = {
        "positions": artifact_record(normalized_positions_path(cell_id)),
        "inspection_pdb": artifact_record(normalized_pdb_path(cell_id)),
        "source_system": artifact_record(normalized_source_xml_path(cell_id)),
        "source_ghost_system": artifact_record(normalized_source_ghost_xml_path(cell_id)),
        "bridge_system": artifact_record(normalized_bridge_xml_path(cell_id)),
    }
    result = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_normalization_cell_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "PASS" if all(checks.values()) else "REJECT_PARENT_BOX_SOURCE",
        "cell_id": cell_id,
        "seed": frozen["seed"],
        "leg": frozen["leg"],
        "inventory_digest": inventory["inventory_digest"],
        "parent_cell_digest": frozen["parent_cell_digest"],
        "parent_artifacts": {
            "source_result": frozen["source_result"],
            "source_artifacts": frozen["source_artifacts"],
            "bridge_artifact": frozen["bridge_artifact"],
        },
        "old_box_vectors_nm": old_box.tolist(),
        "target_box_vectors_nm": target_box.tolist(),
        "volume_relative_delta_before": float(np.linalg.det(old_box) / np.linalg.det(target_box) - 1.0),
        "mapping": normalization["mapping"],
        "counts": source_counts,
        "net_charge_e": source_charge,
        "solute_position_max_abs_delta_nm": solute_delta,
        "solvent_internal_vector_max_abs_delta_nm": internal_delta,
        "ring_geometry": ring_geometry,
        "all_solute_heavy_water_o_min_nm": all_solute_water_min,
        "ghost_contract": p0.parent_p0.compact_force_report(ghost_report),
        "bridge_contract": p0v2._compact_bridge_report(bridge_report),
        "checks": checks,
        "artifacts": artifacts,
        "contexts_created": 0,
        "md_steps": 0,
        "minimization_steps": 0,
        "gpu_used": False,
        "elapsed_s": time.monotonic() - started,
    }
    atomic_write_json(normalized_result_path(cell_id), result)
    return result


def _run_child(args: list[str], log_path: Path) -> subprocess.CompletedProcess[str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        return subprocess.run(
            [str(EXPECTED_ATM_PYTHON), str(RUNNER_PATH), *args],
            cwd=REPO_ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )


def _archive_path(path: Path, reason: str) -> Path | None:
    if not path.exists():
        return None
    stamp = time.strftime("%Y%m%d_%H%M%S")
    target = path.parent / "_archive" / f"{stamp}_{path.name}"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(path), str(target))
    atomic_write_json(
        target / "archive_reason.json",
        {"archived": time.strftime("%Y-%m-%d %H:%M:%S %z"), "reason": reason},
    )
    return target


def validate_normalized_result(cell_id: str, inventory: dict[str, Any]) -> dict[str, Any] | None:
    path = normalized_result_path(cell_id)
    if not path.is_file():
        return None
    result = json.loads(path.read_text(encoding="utf-8"))
    if result.get("status") not in {"PASS", "REJECT_PARENT_BOX_SOURCE"}:
        return None
    if result.get("inventory_digest") != inventory["inventory_digest"]:
        return None
    if result.get("cell_id") != cell_id:
        return None
    for record in (result.get("artifacts") or {}).values():
        artifact_path = REPO_ROOT / record["path"]
        if not artifact_path.is_file() or artifact_record(artifact_path) != record:
            return None
    if result["status"] == "PASS" and not all((result.get("checks") or {}).values()):
        return None
    return result


def normalize_all_sources() -> dict[str, Any]:
    inventory = load_inventory()
    rows = []
    for cell in p0.make_cells():
        existing = validate_normalized_result(cell.cell_id, inventory)
        if existing is not None:
            rows.append(existing)
            if existing["status"] != "PASS":
                break
            continue
        if normalized_cell_dir(cell.cell_id).exists():
            _archive_path(normalized_cell_dir(cell.cell_id), "invalid or partial normalized source")
        completed = _run_child(
            [
                "--worker-normalize",
                cell.cell_id,
                "--expected-inventory-digest",
                inventory["inventory_digest"],
            ],
            NORMALIZED_ROOT / f"{cell.cell_id}.worker.log",
        )
        if completed.returncode != 0:
            raise RuntimeError(f"{cell.cell_id}: normalization worker exited {completed.returncode}")
        result = validate_normalized_result(cell.cell_id, inventory)
        if result is None:
            raise RuntimeError(f"{cell.cell_id}: normalization result is invalid")
        rows.append(result)
        if result["status"] != "PASS":
            break
    status = "PASS" if len(rows) == 6 and all(row["status"] == "PASS" for row in rows) else "REJECT_PARENT_BOX_SOURCE"
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_normalization_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "inventory_digest": inventory["inventory_digest"],
        "n_expected": 6,
        "n_completed": len(rows),
        "pending_cell_ids": [cell.cell_id for cell in p0.make_cells()][len(rows):],
        "cells": [
            {
                "cell_id": row["cell_id"],
                "status": row["status"],
                "result": artifact_record(normalized_result_path(row["cell_id"])),
                "artifacts": row["artifacts"],
                "volume_relative_delta_before": row["volume_relative_delta_before"],
                "stored_min_nm": row["ring_geometry"]["states"]["stored"]["min_ring_water_oxygen_distance_nm"],
                "transformed_min_nm": row["ring_geometry"]["states"]["transformed"]["min_ring_water_oxygen_distance_nm"],
                "all_solute_water_min_nm": row["all_solute_heavy_water_o_min_nm"],
            }
            for row in rows
        ],
        "context_count": 0,
        "gpu_used": False,
        "md_steps": 0,
        "minimization_steps": 0,
    }
    payload["normalized_summary_digest"] = normalized_summary_digest(payload)
    atomic_write_json(NORMALIZED_SUMMARY, payload)
    return payload


def load_normalized_summary(inventory: dict[str, Any]) -> dict[str, Any]:
    if not NORMALIZED_SUMMARY.is_file():
        raise FileNotFoundError(f"missing {NORMALIZED_SUMMARY}; run --normalize-sources")
    summary = json.loads(NORMALIZED_SUMMARY.read_text(encoding="utf-8"))
    if normalized_summary_digest(summary) != summary.get("normalized_summary_digest"):
        raise ValueError("P1SD normalized summary digest drifted")
    if summary.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError("P1SD normalized summary inventory drifted")
    if summary.get("status") != "PASS" or summary.get("n_completed") != 6:
        raise ValueError("P1SD normalized source summary is not a complete PASS")
    for cell in p0.make_cells():
        result = validate_normalized_result(cell.cell_id, inventory)
        if result is None or result.get("status") != "PASS":
            raise ValueError(f"{cell.cell_id}: normalized source is not a valid PASS")
    return summary


def run_keeper_preflight() -> dict[str, Any]:
    inventory = load_inventory()
    summary = load_normalized_summary(inventory)
    cells = []
    for cell in p0.make_cells():
        result = validate_normalized_result(cell.cell_id, inventory)
        assert result is not None
        checks = {
            "normalization_pass": result["status"] == "PASS",
            "all_declared_checks": all(result["checks"].values()),
            "contexts_zero": result["contexts_created"] == 0,
            "gpu_false": result["gpu_used"] is False,
            "sampling_zero": result["md_steps"] == 0 and result["minimization_steps"] == 0,
            "parent_box_exact": np.allclose(
                np.asarray(result["target_box_vectors_nm"]),
                np.asarray(frozen_cell(inventory, cell.cell_id)["parent_box"]["vectors_nm"]),
                atol=1e-10,
                rtol=0.0,
            ),
        }
        cells.append({"cell_id": cell.cell_id, "checks": checks, "passed": all(checks.values())})
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_keeper_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "PASS" if all(row["passed"] for row in cells) else "BLOCK",
        "inventory_digest": inventory["inventory_digest"],
        "normalized_summary_digest": summary["normalized_summary_digest"],
        "n_cells": len(cells),
        "cells": cells,
        "context_count": 0,
        "gpu_used": False,
        "md_steps": 0,
        "minimization_steps": 0,
        "free_energy_estimated": False,
    }
    atomic_write_json(KEEPER_PREFLIGHT, payload)
    return payload


def load_keeper(inventory: dict[str, Any], normalized: dict[str, Any]) -> dict[str, Any]:
    if not KEEPER_PREFLIGHT.is_file():
        raise FileNotFoundError(f"missing {KEEPER_PREFLIGHT}; run --keeper-preflight")
    keeper = json.loads(KEEPER_PREFLIGHT.read_text(encoding="utf-8"))
    if keeper.get("status") != "PASS" or keeper.get("context_count") != 0:
        raise ValueError("P1SD KEEPER is not a context-free PASS")
    if keeper.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError("P1SD KEEPER inventory drifted")
    if keeper.get("normalized_summary_digest") != normalized["normalized_summary_digest"]:
        raise ValueError("P1SD KEEPER normalized summary drifted")
    return keeper


def reference_cell_dir(cell_id: str) -> Path:
    return REFERENCE_ROOT / cell_id


def reference_result_path(cell_id: str) -> Path:
    return reference_cell_dir(cell_id) / "reference_result.json"


def _load_normalized_artifacts(
    cell_id: str, inventory: dict[str, Any]
) -> tuple[dict[str, Any], Any, np.ndarray]:
    from openmm.app import PDBFile

    result = validate_normalized_result(cell_id, inventory)
    if result is None or result.get("status") != "PASS":
        raise ValueError(f"{cell_id}: normalized source is not a valid PASS")
    pdb = PDBFile(str(normalized_pdb_path(cell_id)))
    positions = np.load(normalized_positions_path(cell_id), allow_pickle=False)
    if positions.dtype != np.float64 or positions.shape != (pdb.topology.getNumAtoms(), 3):
        raise ValueError(f"{cell_id}: authoritative positions shape/dtype drifted")
    return result, pdb, positions


def reference_worker(cell_id: str, inventory: dict[str, Any]) -> dict[str, Any]:
    import openmm as mm
    from openmm import unit

    output_dir = reference_cell_dir(cell_id)
    if output_dir.exists():
        raise RuntimeError(f"{cell_id}: Reference output already exists")
    output_dir.mkdir(parents=True)
    started = time.monotonic()
    normalized, pdb, positions_nm = _load_normalized_artifacts(cell_id, inventory)
    frozen = frozen_cell(inventory, cell_id)
    states = frozen["states"]
    positions = positions_nm * unit.nanometer

    source_ghost = mm.XmlSerializer.deserialize(
        normalized_source_ghost_xml_path(cell_id).read_text(encoding="utf-8")
    )
    source_values = p0v2.evaluate_source_apexes(source_ghost, positions, states)
    raw = p0._evaluate_raw_endpoints(source_ghost, positions, states["apex_dplus"])
    max_source_force = max(
        float(np.max(np.abs(value["forces"]))) for value in source_values.values()
    )
    raw_guard = {
        "raw_u0_abs_below_1e8": raw["finite"] and abs(raw["u0_kj_mol"]) < 1e8,
        "raw_u1_abs_below_1e8": raw["finite"] and abs(raw["u1_kj_mol"]) < 1e8,
        "max_force_component_below_1e8": math.isfinite(max_source_force)
        and max_source_force < 1e8,
    }
    raw_guard["passed"] = all(raw_guard.values())

    bridged = mm.XmlSerializer.deserialize(
        normalized_bridge_xml_path(cell_id).read_text(encoding="utf-8")
    )
    bridge_contract = bridge.inspect_apex_bridge(bridged)
    ghost_contract = dg.inspect_dynamic_ghost_force(bridged)
    bridge_values, inactive = p0v2.evaluate_bridge(
        bridged,
        positions,
        states["apex_dplus"],
        [0.0, 0.5, 1.0],
        check_inactive_globals=True,
    )
    plus = source_values["apex_dplus"]
    minus = source_values["apex_dminus"]
    gap_energy = minus["energy_kj_mol"] - plus["energy_kj_mol"]
    gap_forces = minus["forces"] - plus["forces"]
    endpoint_plus = p0v2.parity_record(
        bridge_values["0.00"]["energy_kj_mol"],
        bridge_values["0.00"]["forces"],
        plus["energy_kj_mol"],
        plus["forces"],
    )
    endpoint_minus = p0v2.parity_record(
        bridge_values["1.00"]["energy_kj_mol"],
        bridge_values["1.00"]["forces"],
        minus["energy_kj_mol"],
        minus["forces"],
    )
    expected_mid_energy = plus["energy_kj_mol"] + 0.5 * gap_energy
    expected_mid_forces = plus["forces"] + 0.5 * gap_forces
    midpoint = p0v2.parity_record(
        bridge_values["0.50"]["energy_kj_mol"],
        bridge_values["0.50"]["forces"],
        expected_mid_energy,
        expected_mid_forces,
    )
    derivative = p0v2.derivative_record(
        bridge_values["0.50"]["derivatives_kj_mol"][bridge.BRIDGE_PARAMETER],
        gap_energy,
    )
    assert inactive is not None
    inactive_parity = p0v2.parity_record(
        inactive["changed"]["energy_kj_mol"],
        inactive["changed"]["forces"],
        inactive["baseline"]["energy_kj_mol"],
        inactive["baseline"]["forces"],
    )
    ghost_energies = [value["ghost_energy_kj_mol"] for value in bridge_values.values()]
    ghost_derivatives = [
        value["derivatives_kj_mol"][dg.GHOST_GLOBAL_PARAMETER]
        for value in bridge_values.values()
    ]
    ghost_invariance = {
        "max_ghost_energy_delta_kj_mol": max(ghost_energies) - min(ghost_energies),
        "max_dvdg_delta_kj_mol": max(ghost_derivatives) - min(ghost_derivatives),
    }
    ghost_invariance["passed"] = (
        abs(ghost_invariance["max_ghost_energy_delta_kj_mol"]) <= 1e-5
        and abs(ghost_invariance["max_dvdg_delta_kj_mol"]) <= 1e-5
    )
    finite = all(value["finite"] for value in source_values.values()) and all(
        value["finite"] for value in bridge_values.values()
    )
    checks = {
        "raw_guard": raw_guard["passed"],
        "finite": finite,
        "endpoint_plus": endpoint_plus["passed"],
        "endpoint_minus": endpoint_minus["passed"],
        "midpoint_linear_identity": midpoint["passed"],
        "derivative_identity": derivative["passed"],
        "inactive_globals": inactive_parity["passed"],
        "ghost_invariance": ghost_invariance["passed"],
        "bridge_parameter": bridge_contract["energy_parameter_derivatives"]
        == [bridge.BRIDGE_PARAMETER],
        "ghost_after_atm": ghost_contract["force_index"] > bridge_contract["atm_index"],
        "fixed_parent_box": np.allclose(
            np.asarray(normalized["target_box_vectors_nm"]),
            np.asarray(frozen["parent_box"]["vectors_nm"]),
            atol=1e-10,
            rtol=0.0,
        ),
    }
    status = "P1SD_REFERENCE_PASS" if all(checks.values()) else (
        "REJECT_PARENT_BOX_SOURCE" if not raw_guard["passed"] else "REJECT_BRIDGE"
    )
    result = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_reference_cell_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "cell_id": cell_id,
        "seed": frozen["seed"],
        "leg": frozen["leg"],
        "inventory_digest": inventory["inventory_digest"],
        "normalized_result": artifact_record(normalized_result_path(cell_id)),
        "platform": "Reference",
        "gpu_used": False,
        "md_steps": 0,
        "minimization_steps": 0,
        "free_energy_estimated": False,
        "raw_endpoints": raw,
        "max_source_force_component_kj_mol_nm": max_source_force,
        "raw_guard": raw_guard,
        "source_evaluations": {
            name: p0v2.evaluation_record(value) for name, value in source_values.items()
        },
        "bridge_evaluations": {
            name: p0v2.evaluation_record(value) for name, value in bridge_values.items()
        },
        "endpoint_plus": endpoint_plus,
        "endpoint_minus": endpoint_minus,
        "midpoint": midpoint,
        "derivative": derivative,
        "inactive_parity": inactive_parity,
        "ghost_invariance": ghost_invariance,
        "checks": checks,
        "elapsed_s": time.monotonic() - started,
    }
    del source_ghost, bridged, source_values, bridge_values
    gc.collect()
    atomic_write_json(reference_result_path(cell_id), result)
    return result


def validate_reference_result(cell_id: str, inventory: dict[str, Any]) -> dict[str, Any] | None:
    path = reference_result_path(cell_id)
    if not path.is_file():
        return None
    result = json.loads(path.read_text(encoding="utf-8"))
    if result.get("inventory_digest") != inventory["inventory_digest"]:
        return None
    if result.get("cell_id") != cell_id:
        return None
    allowed = {"P1SD_REFERENCE_PASS", "REJECT_PARENT_BOX_SOURCE", "REJECT_BRIDGE"}
    if result.get("status") not in allowed:
        return None
    if result["status"] == "P1SD_REFERENCE_PASS" and not all(result["checks"].values()):
        return None
    return result


def run_reference_all() -> dict[str, Any]:
    inventory = load_inventory()
    normalized = load_normalized_summary(inventory)
    load_keeper(inventory, normalized)
    rows = []
    for cell in p0.make_cells():
        existing = validate_reference_result(cell.cell_id, inventory)
        if existing is not None:
            rows.append(existing)
            if existing["status"] != "P1SD_REFERENCE_PASS":
                break
            continue
        output_dir = reference_cell_dir(cell.cell_id)
        if output_dir.exists():
            _archive_path(output_dir, "invalid or partial P1SD Reference output")
        completed = _run_child(
            [
                "--worker-reference",
                cell.cell_id,
                "--expected-inventory-digest",
                inventory["inventory_digest"],
            ],
            REFERENCE_ROOT / f"{cell.cell_id}.worker.log",
        )
        if completed.returncode != 0:
            raise RuntimeError(f"{cell.cell_id}: Reference worker exited {completed.returncode}")
        result = validate_reference_result(cell.cell_id, inventory)
        if result is None:
            raise RuntimeError(f"{cell.cell_id}: invalid Reference result")
        rows.append(result)
        if result["status"] != "P1SD_REFERENCE_PASS":
            break
    status = (
        "P1SD_REFERENCE_PASS"
        if len(rows) == 6 and all(row["status"] == "P1SD_REFERENCE_PASS" for row in rows)
        else rows[-1]["status"] if rows else "P1SD_REFERENCE_INTERRUPTED"
    )
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_reference_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "inventory_digest": inventory["inventory_digest"],
        "normalized_summary_digest": normalized["normalized_summary_digest"],
        "n_expected": 6,
        "n_completed": len(rows),
        "pending_cell_ids": [cell.cell_id for cell in p0.make_cells()][len(rows):],
        "cells": [
            {
                "cell_id": row["cell_id"],
                "status": row["status"],
                "raw_u0_kj_mol": row["raw_endpoints"]["u0_kj_mol"],
                "raw_u1_kj_mol": row["raw_endpoints"]["u1_kj_mol"],
                "max_source_force_component_kj_mol_nm": row["max_source_force_component_kj_mol_nm"],
                "result": artifact_record(reference_result_path(row["cell_id"])),
            }
            for row in rows
        ],
        "platform": "Reference",
        "gpu_used": False,
        "md_steps": 0,
        "minimization_steps": 0,
        "free_energy_estimated": False,
    }
    atomic_write_json(REFERENCE_SUMMARY, payload)
    return payload


def load_reference_summary(inventory: dict[str, Any]) -> dict[str, Any]:
    if not REFERENCE_SUMMARY.is_file():
        raise FileNotFoundError(f"missing {REFERENCE_SUMMARY}; run --run-reference")
    summary = json.loads(REFERENCE_SUMMARY.read_text(encoding="utf-8"))
    if summary.get("status") != "P1SD_REFERENCE_PASS" or summary.get("n_completed") != 6:
        raise ValueError("P1SD Reference summary is not a complete PASS")
    if summary.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError("P1SD Reference summary inventory drifted")
    for cell in p0.make_cells():
        result = validate_reference_result(cell.cell_id, inventory)
        if result is None or result.get("status") != "P1SD_REFERENCE_PASS":
            raise ValueError(f"{cell.cell_id}: Reference result is not a valid PASS")
    return summary


def source_prep_dir(cell_id: str) -> Path:
    return SAMPLE_ROOT / cell_id / "source_prep"


def source_prep_result_path(cell_id: str) -> Path:
    return source_prep_dir(cell_id) / "source_prep_result.json"


def source_minimized_positions_path(cell_id: str) -> Path:
    return source_prep_dir(cell_id) / "source_minimized_positions.npy"


def source_minimized_pdb_path(cell_id: str) -> Path:
    return source_prep_dir(cell_id) / "source_minimized.pdb"


def xi_label(xi: float) -> str:
    return f"xi_{float(xi):.2f}".replace(".", "p")


def window_dir(cell_id: str, xi: float) -> Path:
    return SAMPLE_ROOT / cell_id / xi_label(xi)


def window_result_path(cell_id: str, xi: float) -> Path:
    return window_dir(cell_id, xi) / "window_result.json"


def _cuda_properties(protocol: dict[str, Any]) -> dict[str, str]:
    sampling = protocol["sampling"]
    return {
        "DeviceIndex": str(sampling["device_index"]),
        "Precision": str(sampling["precision"]),
    }


def _state_energy_force(context: Any, *, positions: bool = False) -> dict[str, Any]:
    from openmm import unit

    state = context.getState(
        getEnergy=True,
        getForces=True,
        getPositions=positions,
        getParameterDerivatives=True,
    )
    force = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=float,
    )
    payload = {
        "energy_kj_mol": float(
            state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        ),
        "max_force_component_kj_mol_nm": float(np.max(np.abs(force))),
        "derivatives_kj_mol": {
            str(name): float(value)
            for name, value in dict(state.getEnergyParameterDerivatives()).items()
        },
    }
    payload["finite"] = bool(
        math.isfinite(payload["energy_kj_mol"])
        and math.isfinite(payload["max_force_component_kj_mol_nm"])
        and np.isfinite(force).all()
        and all(math.isfinite(value) for value in payload["derivatives_kj_mol"].values())
    )
    if positions:
        payload["positions"] = state.getPositions(asNumpy=True)
    return payload


def source_prep_worker(cell_id: str, inventory: dict[str, Any]) -> dict[str, Any]:
    import openmm as mm
    from openmm import unit
    from openmm.app import PDBFile, Simulation

    protocol = load_protocol_and_verify_freeze()
    frozen = frozen_cell(inventory, cell_id)
    _normalized, pdb, positions_nm = _load_normalized_artifacts(cell_id, inventory)
    output_dir = source_prep_dir(cell_id)
    if output_dir.exists():
        raise RuntimeError(f"{cell_id}: source-prep output already exists")
    output_dir.mkdir(parents=True)
    started = time.monotonic()
    system = mm.XmlSerializer.deserialize(
        normalized_source_ghost_xml_path(cell_id).read_text(encoding="utf-8")
    )
    sampling = protocol["sampling"]
    integrator = mm.LangevinMiddleIntegrator(
        sampling["temperature_K"] * unit.kelvin,
        sampling["friction_per_ps"] / unit.picosecond,
        sampling["timestep_fs"] * unit.femtosecond,
    )
    platform = mm.Platform.getPlatformByName("CUDA")
    simulation = Simulation(
        pdb.topology, system, integrator, platform, _cuda_properties(protocol)
    )
    simulation.context.setPositions(positions_nm * unit.nanometer)
    p0v2._set_source_parameters(simulation.context, frozen["states"]["u0_dplus"], 0.0)
    before = _state_energy_force(simulation.context)
    minimization = sampling["source_minimization"]
    mm.LocalEnergyMinimizer.minimize(
        simulation.context,
        tolerance=minimization["tolerance_kj_mol_nm"]
        * unit.kilojoule_per_mole
        / unit.nanometer,
        maxIterations=int(minimization["max_iterations"]),
    )
    after = _state_energy_force(simulation.context, positions=True)
    final_positions_nm = np.asarray(
        after.pop("positions").value_in_unit(unit.nanometer), dtype=np.float64
    )
    geometry = p0._ring_water_geometry(
        system, pdb.topology, final_positions_nm * unit.nanometer
    )
    final_box = np.asarray(
        [
            vector.value_in_unit(unit.nanometer)
            for vector in simulation.context.getState().getPeriodicBoxVectors()
        ],
        dtype=float,
    )
    target_box = np.asarray(frozen["parent_box"]["vectors_nm"], dtype=float)
    checks = {
        "before_finite": before["finite"],
        "before_gross_guard": before["max_force_component_kj_mol_nm"] < 1e8,
        "after_finite": after["finite"],
        "after_gross_guard": after["max_force_component_kj_mol_nm"] < 1e8,
        "fixed_parent_box": np.allclose(final_box, target_box, atol=1e-10, rtol=0.0),
        "ring_gross_geometry": all(
            row["min_ring_water_oxygen_distance_nm"] >= 0.1
            for row in geometry["states"].values()
        ),
    }
    _atomic_save_npy(source_minimized_positions_path(cell_id), final_positions_nm)
    _write_inspection_pdb(source_minimized_pdb_path(cell_id), pdb.topology, final_positions_nm)
    artifacts = {
        "positions": artifact_record(source_minimized_positions_path(cell_id)),
        "inspection_pdb": artifact_record(source_minimized_pdb_path(cell_id)),
    }
    result = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_source_prep_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "SOURCE_PREP_PASS" if all(checks.values()) else "SOURCE_PREP_FAIL",
        "cell_id": cell_id,
        "inventory_digest": inventory["inventory_digest"],
        "platform": {
            "name": platform.getName(),
            "openmm_version": mm.__version__,
            "properties": _cuda_properties(protocol),
        },
        "state": "physical_u0",
        "ghost_g": 0.0,
        "minimization": minimization,
        "before": before,
        "after": after,
        "geometry": geometry,
        "checks": checks,
        "artifacts": artifacts,
        "md_steps": 0,
        "free_energy_estimated": False,
        "elapsed_s": time.monotonic() - started,
    }
    del simulation, integrator, system
    gc.collect()
    atomic_write_json(source_prep_result_path(cell_id), result)
    return result


def _validate_artifact_map(records: dict[str, dict[str, Any]]) -> bool:
    for record in records.values():
        path = REPO_ROOT / record["path"]
        if not path.is_file() or artifact_record(path) != record:
            return False
    return True


def validate_source_prep(cell_id: str, inventory: dict[str, Any]) -> dict[str, Any] | None:
    path = source_prep_result_path(cell_id)
    if not path.is_file():
        return None
    result = json.loads(path.read_text(encoding="utf-8"))
    if result.get("status") not in {"SOURCE_PREP_PASS", "SOURCE_PREP_FAIL"}:
        return None
    if result.get("cell_id") != cell_id or result.get("inventory_digest") != inventory["inventory_digest"]:
        return None
    if not _validate_artifact_map(result.get("artifacts") or {}):
        return None
    if result["status"] == "SOURCE_PREP_PASS" and not all(result["checks"].values()):
        return None
    return result


def _degrees_of_freedom(system: Any) -> int:
    import openmm as mm

    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    if any(
        isinstance(system.getForce(index), mm.CMMotionRemover)
        for index in range(system.getNumForces())
    ):
        dof -= 3
    if dof <= 0:
        raise ValueError("non-positive system degrees of freedom")
    return dof


def _cross_state_readback(
    context: Any,
    generating_xi: float,
    energy_kj_mol: float,
    derivative_kj_mol: float,
    xi_values: Sequence[float],
    tolerance_config: dict[str, Any],
) -> dict[str, Any]:
    from openmm import unit

    rows = []
    for target_xi in xi_values:
        context.setParameter(bridge.BRIDGE_PARAMETER, float(target_xi))
        observed = float(
            context.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )
        expected = energy_kj_mol + (float(target_xi) - generating_xi) * derivative_kj_mol
        absolute_delta = abs(observed - expected)
        allowed = readback_tolerance_kj_mol(observed, expected, tolerance_config)
        rows.append(
            {
                "target_xi": float(target_xi),
                "observed_kj_mol": observed,
                "expected_kj_mol": expected,
                "abs_delta_kj_mol": absolute_delta,
                "allowed_error_kj_mol": allowed,
                "passed": absolute_delta <= allowed,
            }
        )
    context.setParameter(bridge.BRIDGE_PARAMETER, float(generating_xi))
    return {
        "rows": rows,
        "max_abs_delta_kj_mol": max(row["abs_delta_kj_mol"] for row in rows),
        "max_allowed_error_kj_mol": max(row["allowed_error_kj_mol"] for row in rows),
        "passed": all(row["passed"] for row in rows),
    }


def _bridge_slope_probe(
    context: Any, generating_xi: float, probe_xi: float
) -> dict[str, Any]:
    from openmm import unit

    context.setParameter(bridge.BRIDGE_PARAMETER, float(probe_xi))
    try:
        state = context.getState(getEnergy=True, getParameterDerivatives=True)
        energy = float(
            state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        )
        derivatives = {
            str(name): float(value)
            for name, value in dict(state.getEnergyParameterDerivatives()).items()
        }
        slope = derivatives[bridge.BRIDGE_PARAMETER]
    finally:
        context.setParameter(bridge.BRIDGE_PARAMETER, float(generating_xi))
    restored = float(context.getParameter(bridge.BRIDGE_PARAMETER))
    return {
        "probe_xi": float(probe_xi),
        "probe_energy_kj_mol": energy,
        "slope_kj_mol": slope,
        "finite": math.isfinite(energy) and math.isfinite(slope),
        "generating_xi_restored": restored == float(generating_xi),
    }


def readback_tolerance_kj_mol(
    observed: float, expected: float, config: dict[str, Any]
) -> float:
    floor = float(config["absolute_floor_kj_mol"])
    relative = float(config["relative_to_energy_scale"])
    ceiling = float(config["absolute_ceiling_kj_mol"])
    if not (0.0 < floor <= ceiling and relative > 0.0):
        raise ValueError("invalid P1SD readback tolerance configuration")
    energy_scale = max(1.0, abs(float(observed)), abs(float(expected)))
    return min(ceiling, max(floor, relative * energy_scale))


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty P1SD sample table")
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def _read_openmm_dcd_header(path: Path) -> dict[str, int]:
    with path.open("rb") as handle:
        header = handle.read(20)
    if len(header) != 20:
        raise ValueError(f"truncated DCD header: {path}")
    record_size, magic, frame_count, first_step, interval = struct.unpack(
        "<i4s3i", header
    )
    if record_size != 84 or magic != b"CORD":
        raise ValueError(f"invalid OpenMM DCD header: {path}")
    return {
        "frame_count": frame_count,
        "first_step": first_step,
        "interval": interval,
    }


def window_worker(cell_id: str, xi: float, inventory: dict[str, Any]) -> dict[str, Any]:
    import openmm as mm
    from openmm import app, unit
    from openmm.app import PDBFile, Simulation

    protocol = load_protocol_and_verify_freeze()
    xi_values = [float(value) for value in protocol["sampling"]["xi"]]
    if float(xi) not in xi_values:
        raise ValueError(f"{cell_id}: xi {xi} is outside the frozen grid")
    frozen = frozen_cell(inventory, cell_id)
    source_prep = validate_source_prep(cell_id, inventory)
    if source_prep is None or source_prep.get("status") != "SOURCE_PREP_PASS":
        raise ValueError(f"{cell_id}: source prep is not a valid PASS")
    output_dir = window_dir(cell_id, xi)
    if output_dir.exists():
        raise RuntimeError(f"{cell_id} xi={xi:.2f}: window output already exists")
    output_dir.mkdir(parents=True)
    started = time.monotonic()

    pdb = PDBFile(str(normalized_pdb_path(cell_id)))
    start_positions_nm = np.load(source_minimized_positions_path(cell_id), allow_pickle=False)
    system = mm.XmlSerializer.deserialize(
        normalized_bridge_xml_path(cell_id).read_text(encoding="utf-8")
    )
    sampling = protocol["sampling"]
    integrator_seed = derive_seed(cell_id, xi, "integrator")
    velocity_seed = derive_seed(cell_id, xi, "velocity")
    integrator = mm.LangevinMiddleIntegrator(
        sampling["temperature_K"] * unit.kelvin,
        sampling["friction_per_ps"] / unit.picosecond,
        sampling["timestep_fs"] * unit.femtosecond,
    )
    integrator.setRandomNumberSeed(integrator_seed)
    platform = mm.Platform.getPlatformByName("CUDA")
    simulation = Simulation(
        pdb.topology, system, integrator, platform, _cuda_properties(protocol)
    )
    simulation.context.setPositions(start_positions_nm * unit.nanometer)
    p0v2._set_source_parameters(
        simulation.context, frozen["states"]["apex_dplus"], sampling["ghost_g"]
    )
    simulation.context.setParameter(bridge.BRIDGE_PARAMETER, float(xi))
    before = _state_energy_force(simulation.context)
    window_min = sampling["window_minimization"]
    mm.LocalEnergyMinimizer.minimize(
        simulation.context,
        tolerance=window_min["tolerance_kj_mol_nm"]
        * unit.kilojoule_per_mole
        / unit.nanometer,
        maxIterations=int(window_min["max_iterations"]),
    )
    after_min = _state_energy_force(simulation.context)
    simulation.context.setVelocitiesToTemperature(
        sampling["temperature_K"] * unit.kelvin, velocity_seed
    )
    simulation.step(int(sampling["warmup_steps"]))

    table_path = output_dir / "samples.csv"
    dcd_path = output_dir / "trajectory.dcd"
    dcd_temp = dcd_path.with_name(dcd_path.name + ".tmp")
    checkpoint_path = output_dir / "final.chk"
    final_positions_path = output_dir / "final_positions.npy"
    dof = _degrees_of_freedom(system)
    gas_constant = 0.00831446261815324
    dcd_first_step = int(sampling["warmup_steps"]) + int(
        sampling["steps_per_cycle"]
    )
    sample_rows: list[dict[str, Any]] = []
    explicit_readbacks: dict[str, Any] = {}
    slope_probes: list[dict[str, Any]] = []
    slope_probe_xi = float(
        protocol["cross_state_matrix"]["slope_evaluation"]["probe_xi"]
    )
    with dcd_temp.open("wb") as dcd_handle:
        dcd = app.DCDFile(
            dcd_handle,
            pdb.topology,
            sampling["timestep_fs"] * unit.femtosecond,
            firstStep=dcd_first_step,
            interval=int(sampling["steps_per_cycle"]),
        )
        for cycle in range(1, int(sampling["cycles"]) + 1):
            simulation.step(int(sampling["steps_per_cycle"]))
            state = simulation.context.getState(
                getEnergy=True,
                getPositions=True,
                getParameterDerivatives=True,
            )
            potential = float(
                state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            )
            kinetic = float(
                state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            )
            derivatives = {
                str(name): float(value)
                for name, value in dict(state.getEnergyParameterDerivatives()).items()
            }
            generating_dvdxi = derivatives[bridge.BRIDGE_PARAMETER]
            slope_probe = _bridge_slope_probe(
                simulation.context, float(xi), slope_probe_xi
            )
            bridge_slope = slope_probe["slope_kj_mol"]
            slope_probes.append(slope_probe)
            temperature = 2.0 * kinetic / (dof * gas_constant)
            positions = state.getPositions(asNumpy=True)
            geometry = p0._ring_water_geometry(system, pdb.topology, positions)
            box = state.getPeriodicBoxVectors(asNumpy=True)
            dcd.writeModel(positions, periodicBoxVectors=box)
            sample_rows.append(
                {
                    "cycle": cycle,
                    "step": int(sampling["warmup_steps"])
                    + cycle * int(sampling["steps_per_cycle"]),
                    "xi": f"{float(xi):.2f}",
                    "potential_kj_mol": f"{potential:.12g}",
                    "kinetic_kj_mol": f"{kinetic:.12g}",
                    "temperature_K": f"{temperature:.12g}",
                    "generating_dV_dBridgeXi_kj_mol": f"{generating_dvdxi:.12g}",
                    "bridge_slope_kj_mol": f"{bridge_slope:.12g}",
                    "slope_probe_energy_kj_mol": f"{slope_probe['probe_energy_kj_mol']:.12g}",
                    "stored_ring_water_min_nm": f"{geometry['states']['stored']['min_ring_water_oxygen_distance_nm']:.12g}",
                    "transformed_ring_water_min_nm": f"{geometry['states']['transformed']['min_ring_water_oxygen_distance_nm']:.12g}",
                }
            )
            if cycle in {1, int(sampling["cycles"])}:
                explicit_readbacks[str(cycle)] = _cross_state_readback(
                    simulation.context,
                    float(xi),
                    potential,
                    bridge_slope,
                    xi_values,
                    protocol["cross_state_matrix"]["explicit_readback_tolerance"],
                )
    os.replace(dcd_temp, dcd_path)
    _write_csv(table_path, sample_rows)
    final = _state_energy_force(simulation.context, positions=True)
    final_positions_nm = np.asarray(
        final.pop("positions").value_in_unit(unit.nanometer), dtype=np.float64
    )
    _atomic_save_npy(final_positions_path, final_positions_nm)
    _atomic_write_bytes(checkpoint_path, simulation.context.createCheckpoint())
    final_box = np.asarray(
        [
            vector.value_in_unit(unit.nanometer)
            for vector in simulation.context.getState().getPeriodicBoxVectors()
        ],
        dtype=float,
    )
    target_box = np.asarray(frozen["parent_box"]["vectors_nm"], dtype=float)
    del simulation, integrator, system
    gc.collect()

    from mdtraj.formats import DCDTrajectoryFile

    with DCDTrajectoryFile(str(dcd_path), mode="r") as dcd_reader:
        dcd_frames = len(dcd_reader)
    dcd_header = _read_openmm_dcd_header(dcd_path)
    temperatures = np.asarray([float(row["temperature_K"]) for row in sample_rows])
    stored_minima = np.asarray([float(row["stored_ring_water_min_nm"]) for row in sample_rows])
    transformed_minima = np.asarray(
        [float(row["transformed_ring_water_min_nm"]) for row in sample_rows]
    )
    finite_rows = all(
        all(math.isfinite(float(row[name])) for name in (
            "potential_kj_mol",
            "kinetic_kj_mol",
            "temperature_K",
            "generating_dV_dBridgeXi_kj_mol",
            "bridge_slope_kj_mol",
            "slope_probe_energy_kj_mol",
            "stored_ring_water_min_nm",
            "transformed_ring_water_min_nm",
        ))
        for row in sample_rows
    )
    gates = protocol["execution_gates"]
    checks = {
        "before_finite": before["finite"],
        "after_min_finite": after_min["finite"],
        "rows_exact": len(sample_rows) == sampling["sample_rows"],
        "dcd_frames_exact": dcd_frames == sampling["dcd_frames"],
        "dcd_header_frames_exact": dcd_header["frame_count"]
        == sampling["dcd_frames"],
        "dcd_first_step_exact": dcd_header["first_step"]
        == int(sample_rows[0]["step"]),
        "dcd_interval_exact": dcd_header["interval"]
        == int(sampling["steps_per_cycle"]),
        "final_step_exact": int(sample_rows[-1]["step"]) == sampling["final_step"],
        "finite_rows": finite_rows,
        "slope_probes_finite": all(row["finite"] for row in slope_probes),
        "generating_xi_restored": all(
            row["generating_xi_restored"] for row in slope_probes
        ),
        "temperature_mean": gates["temperature_mean_min_K"]
        <= float(np.mean(temperatures))
        <= gates["temperature_mean_max_K"],
        "temperature_max": float(np.max(temperatures)) <= gates["temperature_sample_max_K"],
        "ring_gross_geometry": min(float(np.min(stored_minima)), float(np.min(transformed_minima)))
        >= gates["ring_water_gross_min_nm"],
        "final_finite": final["finite"],
        "final_force_guard": final["max_force_component_kj_mol_nm"]
        < gates["final_force_component_abs_max_kj_mol_nm"],
        "fixed_parent_box": np.allclose(final_box, target_box, atol=1e-10, rtol=0.0),
        "cross_state_first": explicit_readbacks["1"]["passed"],
        "cross_state_last": explicit_readbacks[str(sampling["cycles"])]["passed"],
    }
    artifacts = {
        "samples": artifact_record(table_path),
        "trajectory": artifact_record(dcd_path),
        "checkpoint": artifact_record(checkpoint_path),
        "final_positions": artifact_record(final_positions_path),
    }
    result = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_window_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "WINDOW_PASS" if all(checks.values()) else "WINDOW_FAIL",
        "cell_id": cell_id,
        "seed": frozen["seed"],
        "leg": frozen["leg"],
        "xi": float(xi),
        "inventory_digest": inventory["inventory_digest"],
        "platform": {
            "name": "CUDA",
            "openmm_version": mm.__version__,
            "properties": _cuda_properties(protocol),
        },
        "seeds": {
            "integrator": {
                "unit_key": seed_unit_key(cell_id, xi, "integrator"),
                "value": integrator_seed,
            },
            "velocity": {
                "unit_key": seed_unit_key(cell_id, xi, "velocity"),
                "value": velocity_seed,
            },
        },
        "protocol": {
            "window_minimization": window_min,
            "warmup_steps": sampling["warmup_steps"],
            "cycles": sampling["cycles"],
            "steps_per_cycle": sampling["steps_per_cycle"],
            "timestep_fs": sampling["timestep_fs"],
            "ghost_g": sampling["ghost_g"],
            "barostat": sampling["barostat"],
        },
        "before": before,
        "after_minimization": after_min,
        "final": final,
        "n_rows": len(sample_rows),
        "dcd_frames": dcd_frames,
        "dcd_header": dcd_header,
        "slope_probe": {
            "xi": slope_probe_xi,
            "n_queries": len(slope_probes),
            "all_finite": all(row["finite"] for row in slope_probes),
            "all_generating_xi_restored": all(
                row["generating_xi_restored"] for row in slope_probes
            ),
        },
        "final_step": int(sample_rows[-1]["step"]),
        "temperature": {
            "mean_K": float(np.mean(temperatures)),
            "min_K": float(np.min(temperatures)),
            "max_K": float(np.max(temperatures)),
        },
        "structure": {
            "stored_min_nm": float(np.min(stored_minima)),
            "transformed_min_nm": float(np.min(transformed_minima)),
            "stored_frac_lt_0p26": float(np.mean(stored_minima < 0.26)),
            "transformed_frac_lt_0p26": float(np.mean(transformed_minima < 0.26)),
        },
        "cross_state_readbacks": explicit_readbacks,
        "checks": checks,
        "artifacts": artifacts,
        "free_energy_estimated": False,
        "elapsed_s": time.monotonic() - started,
    }
    atomic_write_json(window_result_path(cell_id, xi), result)
    return result


def _read_sample_rows(path: Path) -> list[dict[str, float]]:
    rows = []
    with path.open(newline="", encoding="ascii") as handle:
        for raw in csv.DictReader(handle):
            rows.append({name: float(value) for name, value in raw.items()})
    return rows


def validate_window(cell_id: str, xi: float, inventory: dict[str, Any]) -> dict[str, Any] | None:
    path = window_result_path(cell_id, xi)
    if not path.is_file():
        return None
    result = json.loads(path.read_text(encoding="utf-8"))
    if result.get("status") not in {"WINDOW_PASS", "WINDOW_FAIL"}:
        return None
    if result.get("cell_id") != cell_id or result.get("xi") != float(xi):
        return None
    if result.get("inventory_digest") != inventory["inventory_digest"]:
        return None
    if not _validate_artifact_map(result.get("artifacts") or {}):
        return None
    if result["status"] == "WINDOW_PASS" and not all(result["checks"].values()):
        return None
    rows = _read_sample_rows(REPO_ROOT / result["artifacts"]["samples"]["path"])
    if len(rows) != result.get("n_rows"):
        return None
    return result


def cross_state_reduced_potentials(
    sample_rows_by_xi: Sequence[tuple[float, Sequence[dict[str, float]]]],
    xi_values: Sequence[float],
    temperature_K: float,
) -> tuple[np.ndarray, np.ndarray]:
    beta = 1.0 / (0.00831446261815324 * float(temperature_K))
    matrices = []
    n_k = []
    targets = np.asarray(xi_values, dtype=float)
    for generating_xi, rows in sample_rows_by_xi:
        energy = np.asarray([row["potential_kj_mol"] for row in rows], dtype=float)
        derivative = np.asarray([row["bridge_slope_kj_mol"] for row in rows], dtype=float)
        matrix = energy[:, None] + (targets[None, :] - float(generating_xi)) * derivative[:, None]
        matrices.append(beta * matrix)
        n_k.append(len(rows))
    u_nk = np.concatenate(matrices, axis=0)
    u_nk -= u_nk[:, :1]
    return u_nk.T, np.asarray(n_k, dtype=int)


def _run_mbar(u_kn: np.ndarray, n_k: np.ndarray, protocol: dict[str, Any]) -> dict[str, Any]:
    from pymbar import MBAR

    config = protocol["mbar"]
    try:
        mbar = MBAR(
            u_kn,
            n_k,
            maximum_iterations=config["maximum_iterations"],
            relative_tolerance=config["relative_tolerance"],
            solver_protocol=config["solver_protocol"],
        )
    except TypeError:
        mbar = MBAR(u_kn, n_k)
    overlap = np.asarray(mbar.compute_overlap()["matrix"], dtype=float)
    if not np.isfinite(overlap).all():
        raise ValueError("P1SD MBAR overlap matrix is non-finite")
    if float(np.max(np.diag(overlap))) > config["overlap_diagonal_max"]:
        raise ValueError("P1SD MBAR overlap matrix has a non-physical diagonal")
    row_error = float(np.max(np.abs(np.sum(overlap, axis=1) - 1.0)))
    if row_error > config["overlap_row_sum_tolerance"]:
        raise ValueError("P1SD MBAR overlap matrix has a non-physical row sum")
    fe = mbar.compute_free_energy_differences()
    return {
        "overlap_matrix": overlap.tolist(),
        "overlap_row_sum_max_error": row_error,
        "delta_f": np.asarray(fe["Delta_f"], dtype=float).tolist(),
        "d_delta_f": np.asarray(fe["dDelta_f"], dtype=float).tolist(),
    }


def propose_densified_grid(
    xi_values: Sequence[float], interval_min_overlap: Sequence[float], target: float = 0.10
) -> list[float]:
    xis = [float(value) for value in xi_values]
    if len(interval_min_overlap) != len(xis) - 1:
        raise ValueError("overlap interval count does not match xi grid")
    proposed = set(xis)
    for index, overlap in enumerate(interval_min_overlap):
        if not math.isfinite(float(overlap)) or float(overlap) < target:
            proposed.add(0.5 * (xis[index] + xis[index + 1]))
    return sorted(proposed)


def aggregate_cell(cell_id: str, inventory: dict[str, Any]) -> dict[str, Any]:
    protocol = load_protocol_and_verify_freeze()
    xis = [float(value) for value in protocol["sampling"]["xi"]]
    windows = []
    rows_by_xi = []
    for xi in xis:
        result = validate_window(cell_id, xi, inventory)
        if result is None or result.get("status") != "WINDOW_PASS":
            raise ValueError(f"{cell_id} xi={xi:.2f}: window is not a valid PASS")
        rows = _read_sample_rows(REPO_ROOT / result["artifacts"]["samples"]["path"])
        rows_by_xi.append((xi, rows))
        windows.append(result)
    u_kn, n_k = cross_state_reduced_potentials(
        rows_by_xi, xis, protocol["sampling"]["temperature_K"]
    )
    mbar = _run_mbar(u_kn, n_k, protocol)
    matrix = np.asarray(mbar["overlap_matrix"], dtype=float)
    adjacent = [float(matrix[index, index + 1]) for index in range(len(xis) - 1)]
    beta = 1.0 / (
        0.00831446261815324 * float(protocol["sampling"]["temperature_K"])
    )
    delta_f = np.asarray(mbar["delta_f"], dtype=float)
    d_delta_f = np.asarray(mbar["d_delta_f"], dtype=float)
    derivative_g = {}
    try:
        from pymbar import timeseries

        for xi, rows in rows_by_xi:
            values = np.asarray([row["bridge_slope_kj_mol"] for row in rows])
            derivative_g[f"{xi:.2f}"] = float(timeseries.statistical_inefficiency(values))
    except Exception as exc:
        derivative_g = {"status": f"unavailable: {type(exc).__name__}: {exc}"}
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_cell_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "P1SD_CELL_PASS",
        "cell_id": cell_id,
        "inventory_digest": inventory["inventory_digest"],
        "xi": xis,
        "n_k": n_k.tolist(),
        "mbar": mbar,
        "adjacent_overlap": adjacent,
        "min_adjacent_overlap": min(adjacent),
        "exploratory_delta_g_1_minus_0_kj_mol": float(delta_f[0, -1] / beta),
        "exploratory_sigma_kj_mol": float(d_delta_f[0, -1] / beta),
        "exploratory_only_not_production": True,
        "derivative_statistical_inefficiency": derivative_g,
        "structure": {
            "stored_min_nm": min(row["structure"]["stored_min_nm"] for row in windows),
            "transformed_min_nm": min(
                row["structure"]["transformed_min_nm"] for row in windows
            ),
            "stored_frac_lt_0p26_max": max(
                row["structure"]["stored_frac_lt_0p26"] for row in windows
            ),
            "transformed_frac_lt_0p26_max": max(
                row["structure"]["transformed_frac_lt_0p26"] for row in windows
            ),
        },
        "windows": [
            {
                "xi": row["xi"],
                "result": artifact_record(window_result_path(cell_id, row["xi"])),
            }
            for row in windows
        ],
    }
    path = SAMPLE_ROOT / cell_id / "cell_summary.json"
    atomic_write_json(path, payload)
    return payload


def _write_run_state(
    inventory: dict[str, Any], status: str, completed: Sequence[str], active: str | None
) -> None:
    atomic_write_json(
        RUN_STATE,
        {
            "schema": "updd_dynamic_ghost_union_bridge_p1sd_run_state_v1",
            "updated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "status": status,
            "inventory_digest": inventory["inventory_digest"],
            "completed_units": list(completed),
            "active_unit": active,
        },
    )


def run_sampled_schedule() -> dict[str, Any]:
    inventory = load_inventory()
    normalized = load_normalized_summary(inventory)
    load_keeper(inventory, normalized)
    load_reference_summary(inventory)
    protocol = load_protocol_and_verify_freeze()
    xis = [float(value) for value in protocol["sampling"]["xi"]]
    completed_units: list[str] = []
    cell_summaries = []
    for cell in p0.make_cells():
        source_prep = validate_source_prep(cell.cell_id, inventory)
        if source_prep is None:
            prep_dir = source_prep_dir(cell.cell_id)
            if prep_dir.exists():
                _archive_path(prep_dir, "invalid or partial P1SD source prep")
            unit_name = f"{cell.cell_id}/source_prep"
            _write_run_state(inventory, "RUNNING", completed_units, unit_name)
            completed = _run_child(
                [
                    "--worker-source-prep",
                    cell.cell_id,
                    "--expected-inventory-digest",
                    inventory["inventory_digest"],
                ],
                SAMPLE_ROOT / cell.cell_id / "source_prep.worker.log",
            )
            if completed.returncode != 0:
                _write_run_state(inventory, "FAILED", completed_units, unit_name)
                raise RuntimeError(f"{unit_name}: worker exited {completed.returncode}")
            source_prep = validate_source_prep(cell.cell_id, inventory)
        if source_prep is None or source_prep.get("status") != "SOURCE_PREP_PASS":
            raise RuntimeError(f"{cell.cell_id}: source prep failed")
        completed_units.append(f"{cell.cell_id}/source_prep")

        for xi in xis:
            existing = validate_window(cell.cell_id, xi, inventory)
            unit_name = f"{cell.cell_id}/{xi_label(xi)}"
            if existing is None:
                target = window_dir(cell.cell_id, xi)
                if target.exists():
                    _archive_path(target, "invalid or partial P1SD window")
                _write_run_state(inventory, "RUNNING", completed_units, unit_name)
                completed = _run_child(
                    [
                        "--worker-window",
                        cell.cell_id,
                        "--xi",
                        f"{xi:.2f}",
                        "--expected-inventory-digest",
                        inventory["inventory_digest"],
                    ],
                    SAMPLE_ROOT / cell.cell_id / f"{xi_label(xi)}.worker.log",
                )
                if completed.returncode != 0:
                    _write_run_state(inventory, "FAILED", completed_units, unit_name)
                    raise RuntimeError(f"{unit_name}: worker exited {completed.returncode}")
                existing = validate_window(cell.cell_id, xi, inventory)
            if existing is None or existing.get("status") != "WINDOW_PASS":
                _write_run_state(inventory, "FAILED", completed_units, unit_name)
                raise RuntimeError(f"{unit_name}: window failed")
            completed_units.append(unit_name)
        cell_summaries.append(aggregate_cell(cell.cell_id, inventory))

    interval_mins = [
        min(row["adjacent_overlap"][index] for row in cell_summaries)
        for index in range(len(xis) - 1)
    ]
    severe = any(value < protocol["overlap_decision"]["severe_bottleneck_below"] for value in interval_mins)
    needs_densify = any(value < protocol["overlap_decision"]["adjacent_overlap_target"] for value in interval_mins)
    if severe:
        status = "SEVERE_BOTTLENECK"
    elif needs_densify:
        status = "DENSIFY_REQUIRED"
    else:
        status = "GRID_CANDIDATE"
    proposed = propose_densified_grid(
        xis, interval_mins, protocol["overlap_decision"]["adjacent_overlap_target"]
    )
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p1sd_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "inventory_digest": inventory["inventory_digest"],
        "platform": "CUDA",
        "gpu_device_index": protocol["sampling"]["device_index"],
        "precision": protocol["sampling"]["precision"],
        "barostat": False,
        "initial_xi": xis,
        "interval_min_overlap_over_six_cells": interval_mins,
        "proposed_xi": proposed,
        "automatic_next_launch": False,
        "production_free_energy_estimated": False,
        "pilot_samples_reusable_in_production": False,
        "round_trips_measured": False,
        "cells": [
            {
                "cell_id": row["cell_id"],
                "status": row["status"],
                "adjacent_overlap": row["adjacent_overlap"],
                "min_adjacent_overlap": row["min_adjacent_overlap"],
                "exploratory_delta_g_1_minus_0_kj_mol": row[
                    "exploratory_delta_g_1_minus_0_kj_mol"
                ],
                "exploratory_sigma_kj_mol": row["exploratory_sigma_kj_mol"],
                "structure": row["structure"],
                "result": artifact_record(SAMPLE_ROOT / row["cell_id"] / "cell_summary.json"),
            }
            for row in cell_summaries
        ],
    }
    atomic_write_json(SAMPLE_SUMMARY, payload)
    lines = [
        "# P1SD Schedule Discovery Summary",
        "",
        f"Status: `{status}`",
        "",
        "Pilot samples and exploratory MBAR values are not production estimates.",
        "",
        "| cell | min adjacent O | exploratory dG (kJ/mol) | sigma (kJ/mol) |",
        "|---|---:|---:|---:|",
    ]
    for row in payload["cells"]:
        lines.append(
            f"| {row['cell_id']} | {row['min_adjacent_overlap']:.6f} | "
            f"{row['exploratory_delta_g_1_minus_0_kj_mol']:.6f} | "
            f"{row['exploratory_sigma_kj_mol']:.6f} |"
        )
    lines.extend(
        [
            "",
            "Interval minima: " + ", ".join(f"{value:.6f}" for value in interval_mins),
            "",
            "Proposed xi: " + ", ".join(f"{value:.3f}" for value in proposed),
            "",
            "No automatic densified or production launch is authorized.",
        ]
    )
    atomic_write_text(SAMPLE_REPORT, "\n".join(lines) + "\n")
    _write_run_state(inventory, "COMPLETE", completed_units, None)
    return payload


def _parse_cli(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_mutually_exclusive_group(required=True)
    stages.add_argument("--build-inventory", action="store_true")
    stages.add_argument("--normalize-sources", action="store_true")
    stages.add_argument("--keeper-preflight", action="store_true")
    stages.add_argument("--run-reference", action="store_true")
    stages.add_argument("--run", action="store_true")
    stages.add_argument("--worker-normalize")
    stages.add_argument("--worker-reference")
    stages.add_argument("--worker-source-prep")
    stages.add_argument("--worker-window")
    parser.add_argument("--xi", type=float)
    parser.add_argument("--expected-inventory-digest")
    return parser.parse_args(argv)


def main() -> int:
    args = _parse_cli()
    if args.build_inventory:
        print(json.dumps(build_prebuild_inventory(), indent=2))
        return 0
    if args.normalize_sources:
        print(json.dumps(normalize_all_sources(), indent=2))
        return 0
    if args.keeper_preflight:
        print(json.dumps(run_keeper_preflight(), indent=2))
        return 0
    if args.run_reference:
        print(json.dumps(run_reference_all(), indent=2))
        return 0
    if args.run:
        print(json.dumps(run_sampled_schedule(), indent=2))
        return 0
    if args.worker_normalize:
        inventory = load_inventory()
        if args.expected_inventory_digest != inventory["inventory_digest"]:
            raise ValueError("normalization worker inventory digest mismatch")
        result = normalize_source_worker(args.worker_normalize, inventory)
        print(json.dumps({"cell_id": args.worker_normalize, "status": result["status"]}))
        return 0 if result["status"] == "PASS" else 2
    if args.worker_reference:
        inventory = load_inventory()
        if args.expected_inventory_digest != inventory["inventory_digest"]:
            raise ValueError("Reference worker inventory digest mismatch")
        result = reference_worker(args.worker_reference, inventory)
        print(json.dumps({"cell_id": args.worker_reference, "status": result["status"]}))
        return 0 if result["status"] == "P1SD_REFERENCE_PASS" else 2
    if args.worker_source_prep:
        inventory = load_inventory()
        if args.expected_inventory_digest != inventory["inventory_digest"]:
            raise ValueError("source-prep worker inventory digest mismatch")
        result = source_prep_worker(args.worker_source_prep, inventory)
        print(json.dumps({"cell_id": args.worker_source_prep, "status": result["status"]}))
        return 0 if result["status"] == "SOURCE_PREP_PASS" else 2
    if args.worker_window:
        if args.xi is None:
            raise ValueError("--worker-window requires --xi")
        inventory = load_inventory()
        if args.expected_inventory_digest != inventory["inventory_digest"]:
            raise ValueError("window worker inventory digest mismatch")
        result = window_worker(args.worker_window, args.xi, inventory)
        print(
            json.dumps(
                {
                    "cell_id": args.worker_window,
                    "xi": args.xi,
                    "status": result["status"],
                }
            )
        )
        return 0 if result["status"] == "WINDOW_PASS" else 2
    raise AssertionError("unreachable P1SD CLI stage")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as exc:
        traceback.print_exc()
        print(json.dumps({"status": "ERROR", "error": f"{type(exc).__name__}: {exc}"}))
        raise SystemExit(1)
