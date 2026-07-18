#!/usr/bin/env python3
"""Audit an ATM-nested transformed-coordinate dynamic ghost on Reference."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path
from statistics import median
from typing import Any, Iterable

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
UTILS_DIR = REPO_ROOT / "utils"
if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))

import trackb_dynamic_ghost as dg  # noqa: E402


ANALYSIS_DIR = Path(__file__).resolve().parent
PREREGISTRATION = ANALYSIS_DIR / "PREREGISTRATION.md"
SCIVAL_VERDICT = ANALYSIS_DIR / "SCIVAL_VERDICT.md"
PROTOCOL_PATH = ANALYSIS_DIR / "protocol.json"
PROTOCOL_AMENDMENT_PATH = ANALYSIS_DIR / "protocol_amendment_r1.json"
PROTOCOL_AMENDMENT_R2_PATH = ANALYSIS_DIR / "protocol_amendment_r2.json"
FROZEN_MANIFEST = ANALYSIS_DIR / "FROZEN_MANIFEST_R2.json"
RUNNER_PATH = Path(__file__).resolve()

P1_DIR = REPO_ROOT / "analysis/dynamic_ghost_union_bridge_p1sd_20260717"
P1_INVENTORY = P1_DIR / "P1SD_PREBUILD_INVENTORY.json"
P1_NORMALIZED_SUMMARY = P1_DIR / "P1SD_NORMALIZED_SOURCE_SUMMARY.json"
P1_PATH = P1_DIR / "P1SD_PATH_DIAGNOSIS.json"
P1CAP_PATH = (
    REPO_ROOT
    / "analysis/raw_gap_cap_activation_p1cap_20260717"
    / "P1CAP_PATH_DIAGNOSIS.json"
)
NORMALIZED_ROOT = P1_DIR / "normalized_sources"

PREBUILD_INVENTORY = ANALYSIS_DIR / "P0_PREBUILD_INVENTORY.json"
KEEPER_ROOT = ANALYSIS_DIR / "keeper_cells"
KEEPER_PREFLIGHT = ANALYSIS_DIR / "P0_KEEPER_PREFLIGHT.json"
REFERENCE_ROOT = ANALYSIS_DIR / "reference_audit"
RUN_STATE = ANALYSIS_DIR / "P0_RUN_STATE.json"
RUNNER_SUMMARY = ANALYSIS_DIR / "P0_RUNNER_SUMMARY.json"
RUNNER_REPORT = ANALYSIS_DIR / "P0_RUNNER_SUMMARY.md"

EXPECTED_ATM_PYTHON = Path("/home/san/miniconda3/envs/atm/bin/python")
EXPECTED_STATUS = (
    "FROZEN_AFTER_KEEPER_SMOKE_COORDINATE_PRECISION_DISCOVERY_BEFORE_"
    "OFFICIAL_OUTPUT"
)
KCAL_TO_KJ = 4.184


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_digest(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def row_digest(rows: Iterable[Any]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(
            json.dumps(
                row,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ).encode("utf-8")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def artifact_record(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {
        "path": relative(path),
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def verify_artifact(record: dict[str, Any]) -> Path:
    path = REPO_ROOT / record["path"]
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size != int(record["size"]):
        raise ValueError(f"artifact size drift: {path}")
    if sha256_file(path) != record["sha256"]:
        raise ValueError(f"artifact hash drift: {path}")
    return path


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def atomic_write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def _require_atm_python() -> None:
    if Path(sys.executable).resolve() != EXPECTED_ATM_PYTHON.resolve():
        raise RuntimeError(
            f"P0 requires {EXPECTED_ATM_PYTHON}, got {sys.executable}"
        )


def _hash_map(base: Path, mapping: dict[str, str], label: str) -> None:
    for name, expected in mapping.items():
        path = base / name
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(
                f"{label} drift: {path} {observed} != {expected}"
            )


def _expected_cells() -> list[str]:
    return [
        "w4a_union_s101_bound",
        "w4a_union_s127_bound",
        "w4a_union_s163_bound",
        "w4a_union_s101_free",
        "w4a_union_s127_free",
        "w4a_union_s163_free",
    ]


def load_protocol_and_verify_freeze() -> dict[str, Any]:
    frozen = json.loads(FROZEN_MANIFEST.read_text(encoding="utf-8"))
    if frozen.get("status") != EXPECTED_STATUS:
        raise ValueError("P0 frozen-manifest status drifted")
    _hash_map(ANALYSIS_DIR, frozen.get("files") or {}, "P0 frozen file")
    _hash_map(REPO_ROOT, frozen.get("parent_evidence") or {}, "P0 parent evidence")

    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    amendment = json.loads(
        PROTOCOL_AMENDMENT_PATH.read_text(encoding="utf-8")
    )
    amendment_r2 = json.loads(
        PROTOCOL_AMENDMENT_R2_PATH.read_text(encoding="utf-8")
    )
    if (
        protocol.get("status") != "FROZEN_BEFORE_IMPLEMENTATION_OR_OUTPUT"
        or sha256_file(PROTOCOL_PATH)
        != amendment.get("base_protocol_sha256")
        or amendment.get("status")
        != "FROZEN_AFTER_DEV_SMOKE_BASELINE_DISCOVERY_BEFORE_OFFICIAL_OUTPUT"
        or amendment.get("scival_verdict") != "CONDITIONAL_APPROVE"
        or amendment.get("development_smoke_scope")
        != {
            "cell_id": "w4a_union_s101_free",
            "official_output": False,
            "files_written": 0,
            "finding": "endpoint_specific_source_ghost_baseline",
        }
        or amendment.get("gate_semantics")
        != {
            "absolute_g_difference": "diagnostic_only",
            "probe_threshold_input": (
                "probe_g_difference_minus_source_g_difference"
            ),
            "numeric_threshold_changes": 0,
        }
        or any((amendment.get("claims") or {}).values())
    ):
        raise ValueError("P0 R1 amendment declaration drifted")
    if (
        amendment_r2.get("status") != EXPECTED_STATUS
        or amendment_r2.get("scival_verdict") != "CONDITIONAL_APPROVE"
        or amendment_r2.get("base_protocol_sha256")
        != sha256_file(PROTOCOL_PATH)
        or amendment_r2.get("r1_amendment_sha256")
        != sha256_file(PROTOCOL_AMENDMENT_PATH)
        or amendment_r2.get("r1_manifest_sha256")
        != sha256_file(ANALYSIS_DIR / "FROZEN_MANIFEST_R1.json")
        or amendment_r2.get("development_smoke_scope")
        != {
            "cell_id": "w4a_union_s101_free",
            "context_count": 0,
            "official_output": False,
            "finding": (
                "normalized_coordinate_precision_exceeds_r0_offset_tolerance"
            ),
        }
        or amendment_r2.get("gate_semantics")
        != {
            "hamiltonian_changes": 0,
            "scientific_threshold_changes": 0,
            "exact_common_transform_indices_required": True,
            "ring_u0_displacement_allowed": False,
            "water_displacement_allowed": False,
        }
        or any((amendment_r2.get("claims") or {}).values())
    ):
        raise ValueError("P0 R2 amendment declaration drifted")
    protocol["status"] = EXPECTED_STATUS
    protocol["scival_verdict"] = amendment_r2["scival_verdict"]
    protocol["probe"].update(amendment["probe_updates"])
    protocol["transformation"].update(
        amendment_r2["transformation_updates"]
    )
    if (
        protocol.get("status") != EXPECTED_STATUS
        or protocol.get("scival_verdict") != "CONDITIONAL_APPROVE"
    ):
        raise ValueError("P0 status or SciVal verdict drifted")
    dg.validate_atm_nested_ghost_protocol(protocol)
    system = protocol.get("system") or {}
    if system != {
        "mutation_spec": "w4a_trp_ala_res4",
        "cells": _expected_cells(),
        "source_root": (
            "analysis/dynamic_ghost_union_bridge_p1sd_20260717/"
            "normalized_sources"
        ),
        "source_files": [
            "normalized_source_system.xml",
            "normalized_source.pdb",
            "positions_float64.npy",
            "normalization_result.json",
        ],
        "platform": "Reference",
        "md_steps": 0,
        "minimization_steps": 0,
        "gpu_allowed": False,
        "worker_unit": "cell",
        "one_context_per_system_variant": True,
    }:
        raise ValueError("P0 System declaration drifted")
    probe = protocol.get("probe") or {}
    if probe != {
        "target_atom_name": "NE1",
        "outward_vector": "stored_ring_centroid_to_stored_NE1",
        "distance_fraction_of_target_wca_cutoff": 0.85,
        "water_residue_atom_count": 3,
        "water_selection": (
            "ghost_inactive_both_endpoints_then_maximum_H_centroid_"
            "outward_dot_then_lowest_oxygen_index"
        ),
        "coordinate_operation": "rigid_translation_only",
        "coordinate_sets": ["stored_probe", "transformed_probe"],
        "sampling_reuse_allowed": False,
        "source_clearance_nm": 0.42,
        "isolation_metric": (
            "endpoint_matched_source_baseline_double_difference"
        ),
    }:
        raise ValueError("P0 probe declaration drifted")
    gates = protocol.get("gates") or {}
    if gates != {
        "energy_parity_abs_kj_mol": 1.0e-5,
        "force_parity_max_component_kj_mol_nm": 1.0e-5,
        "source_isolated_energy_abs_max_kj_mol": 1000.0,
        "source_isolated_force_max_component_kj_mol_nm": 100000.0,
        "probe_intended_energy_min_kj_mol": 0.1,
        "probe_opposite_energy_abs_max_kj_mol": 1.0e-6,
        "probe_intended_radial_force_min_kj_mol_nm": 1.0,
        "probe_opposite_force_max_component_kj_mol_nm": 1.0e-5,
        "probe_force_balance_max_component_kj_mol_nm": 1.0e-5,
        "finite_dvdg_required": True,
        "all_six_cells_required": True,
    }:
        raise ValueError("P0 gate declaration drifted")
    if any((protocol.get("claims") or {}).values()):
        raise ValueError("P0 prohibited scientific claim was enabled")

    parent_inventory = json.loads(P1_INVENTORY.read_text(encoding="utf-8"))
    normalized = json.loads(P1_NORMALIZED_SUMMARY.read_text(encoding="utf-8"))
    p1_path = json.loads(P1_PATH.read_text(encoding="utf-8"))
    p1cap_path = json.loads(P1CAP_PATH.read_text(encoding="utf-8"))
    if (
        parent_inventory.get("status") != "FROZEN_BEFORE_P1SD_OUTPUT"
        or normalized.get("status") != "PASS"
        or normalized.get("n_completed") != 6
        or normalized.get("n_expected") != 6
        or normalized.get("inventory_digest")
        != parent_inventory.get("inventory_digest")
    ):
        raise ValueError("P1SD normalized parent declaration drifted")
    if (
        p1_path.get("path_verdict")
        != "FAIL_P1SD_APEX_DEGENERACY_AND_TRANSFORMED_VOID_PENETRATION"
        or p1cap_path.get("path_verdict")
        != "CAP_TAIL_PRESENT_DIAGNOSTIC_ONLY"
    ):
        raise ValueError("P0 trigger PATH evidence drifted")
    return protocol


def _implementation_records() -> dict[str, dict[str, Any]]:
    return {
        "runner": artifact_record(RUNNER_PATH),
        "dynamic_ghost": artifact_record(
            REPO_ROOT / "utils/trackb_dynamic_ghost.py"
        ),
    }


def _summary_cell(
    normalized_summary: dict[str, Any], cell_id: str
) -> dict[str, Any]:
    rows = [
        row
        for row in normalized_summary.get("cells", [])
        if row.get("cell_id") == cell_id
    ]
    if len(rows) != 1:
        raise ValueError(
            f"expected one normalized summary row for {cell_id}, found {len(rows)}"
        )
    return rows[0]


def _parent_cell(parent_inventory: dict[str, Any], cell_id: str) -> dict[str, Any]:
    rows = [
        row
        for row in parent_inventory.get("cells", [])
        if row.get("cell_id") == cell_id
    ]
    if len(rows) != 1:
        raise ValueError(
            f"expected one P1SD inventory row for {cell_id}, found {len(rows)}"
        )
    return rows[0]


def inventory_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("inventory_digest", None)
    return canonical_digest(canonical)


def build_inventory() -> dict[str, Any]:
    _require_atm_python()
    protocol = load_protocol_and_verify_freeze()
    parent_inventory = json.loads(P1_INVENTORY.read_text(encoding="utf-8"))
    normalized_summary = json.loads(
        P1_NORMALIZED_SUMMARY.read_text(encoding="utf-8")
    )

    cells: list[dict[str, Any]] = []
    for cell_id in protocol["system"]["cells"]:
        parent = _parent_cell(parent_inventory, cell_id)
        summary = _summary_cell(normalized_summary, cell_id)
        normalized_dir = NORMALIZED_ROOT / cell_id
        normalization_path = normalized_dir / "normalization_result.json"
        normalization = json.loads(
            normalization_path.read_text(encoding="utf-8")
        )
        if (
            normalization.get("status") != "PASS"
            or normalization.get("cell_id") != cell_id
            or normalization.get("inventory_digest")
            != parent_inventory["inventory_digest"]
        ):
            raise ValueError(f"{cell_id}: normalization declaration drifted")
        artifacts = {
            "source_system": artifact_record(
                normalized_dir / "normalized_source_system.xml"
            ),
            "topology_pdb": artifact_record(
                normalized_dir / "normalized_source.pdb"
            ),
            "positions": artifact_record(
                normalized_dir / "positions_float64.npy"
            ),
            "normalization_result": artifact_record(normalization_path),
        }
        declared = normalization.get("artifacts") or {}
        expected_pairs = {
            "source_system": "source_system",
            "topology_pdb": "inspection_pdb",
            "positions": "positions",
        }
        for local_name, parent_name in expected_pairs.items():
            if artifacts[local_name] != declared.get(parent_name):
                raise ValueError(
                    f"{cell_id}: normalization artifact {local_name} drifted"
                )
            if artifacts[local_name] != summary["artifacts"].get(parent_name):
                raise ValueError(
                    f"{cell_id}: normalized summary artifact {local_name} drifted"
                )
        if artifacts["normalization_result"] != summary.get("result"):
            raise ValueError(
                f"{cell_id}: normalization-result artifact drifted"
            )
        cells.append(
            {
                "cell_id": cell_id,
                "seed": parent["seed"],
                "leg": parent["leg"],
                "states": {
                    "u0_dplus": parent["states"]["u0_dplus"],
                    "u1_dminus": parent["states"]["u1_dminus"],
                },
                "parent_box": parent["parent_box"],
                "parent_cell_digest": canonical_digest(parent),
                "normalization_status": normalization["status"],
                "artifacts": artifacts,
            }
        )

    import openmm

    payload = {
        "schema": "updd_atm_nested_dynamic_ghost_p0_inventory_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "FROZEN_BEFORE_P0_OUTPUT",
        "python": sys.executable,
        "runtime": {
            "openmm": openmm.__version__,
            "numpy": np.__version__,
        },
        "protocol": {
            "base": artifact_record(PROTOCOL_PATH),
            "amendment_r1": artifact_record(PROTOCOL_AMENDMENT_PATH),
            "amendment_r2": artifact_record(PROTOCOL_AMENDMENT_R2_PATH),
        },
        "frozen_manifest": artifact_record(FROZEN_MANIFEST),
        "parent_evidence": {
            "p1sd_inventory": artifact_record(P1_INVENTORY),
            "p1sd_normalized_summary": artifact_record(P1_NORMALIZED_SUMMARY),
            "p1sd_path": artifact_record(P1_PATH),
            "p1cap_path": artifact_record(P1CAP_PATH),
        },
        "parent_inventory_digest": parent_inventory["inventory_digest"],
        "implementation": _implementation_records(),
        "counts": {"cells": len(cells), "bound": 3, "free": 3},
        "cells": cells,
    }
    payload["inventory_digest"] = inventory_digest(payload)
    return payload


def load_inventory(expected_digest: str | None = None) -> dict[str, Any]:
    if not PREBUILD_INVENTORY.is_file():
        raise FileNotFoundError(
            f"missing {PREBUILD_INVENTORY}; run --prepare-inventory first"
        )
    inventory = json.loads(PREBUILD_INVENTORY.read_text(encoding="utf-8"))
    if inventory.get("inventory_digest") != inventory_digest(inventory):
        raise ValueError("stored P0 inventory digest is invalid")
    if (
        expected_digest is not None
        and inventory["inventory_digest"] != expected_digest
    ):
        raise ValueError("worker inventory digest mismatch")
    return inventory


def verify_inventory() -> dict[str, Any]:
    frozen = load_inventory()
    current = build_inventory()
    if current["inventory_digest"] != frozen["inventory_digest"]:
        raise ValueError(
            "P0 source/code drift after inventory freeze: "
            f"{current['inventory_digest']} != {frozen['inventory_digest']}"
        )
    return frozen


def frozen_cell(inventory: dict[str, Any], cell_id: str) -> dict[str, Any]:
    rows = [
        row for row in inventory.get("cells", []) if row.get("cell_id") == cell_id
    ]
    if len(rows) != 1:
        raise ValueError(
            f"expected one frozen P0 cell {cell_id}, found {len(rows)}"
        )
    return rows[0]


def _verify_worker_inputs(
    inventory: dict[str, Any], cell_id: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    protocol = load_protocol_and_verify_freeze()
    if inventory.get("implementation") != _implementation_records():
        raise ValueError(f"{cell_id}: implementation drift after inventory freeze")
    cell = frozen_cell(inventory, cell_id)
    for record in cell["artifacts"].values():
        verify_artifact(record)
    return protocol, cell


def _quantity_float(value: Any, target_unit: Any) -> float:
    if hasattr(value, "value_in_unit"):
        return float(value.value_in_unit(target_unit))
    return float(value)


def _vector_nm(value: Any) -> list[float]:
    from openmm import unit

    vector = value.value_in_unit(unit.nanometer)
    return [float(vector[index]) for index in range(3)]


def _transformation_rows(atm_force: Any) -> Iterable[list[Any]]:
    import openmm as mm

    for index in range(atm_force.getNumParticles()):
        transformation = atm_force.getParticleTransformation(index)
        if isinstance(transformation, mm.FixedDisplacement):
            yield [
                index,
                "FixedDisplacement",
                _vector_nm(transformation.getFixedDisplacement0()),
                _vector_nm(transformation.getFixedDisplacement1()),
            ]
        elif isinstance(transformation, mm.ParticleOffsetDisplacement):
            yield [
                index,
                "ParticleOffsetDisplacement",
                int(transformation.getDestinationParticle0()),
                int(transformation.getOriginParticle0()),
                int(transformation.getDestinationParticle1()),
                int(transformation.getOriginParticle1()),
            ]
        else:
            raise dg.DynamicGhostError(
                f"particle {index} has unsupported ATM transformation "
                f"{type(transformation).__name__}"
            )


def _system_declaration(system: Any) -> dict[str, Any]:
    import openmm as mm
    from openmm import unit

    atm_index, atm_force = dg.find_atm_force(system)
    masses = row_digest(
        [
            index,
            _quantity_float(
                system.getParticleMass(index),
                unit.dalton,
            ),
        ]
        for index in range(system.getNumParticles())
    )
    constraints = row_digest(
        [
            index,
            int(system.getConstraintParameters(index)[0]),
            int(system.getConstraintParameters(index)[1]),
            _quantity_float(
                system.getConstraintParameters(index)[2],
                unit.nanometer,
            ),
        ]
        for index in range(system.getNumConstraints())
    )
    virtual_indices = [
        index
        for index in range(system.getNumParticles())
        if system.isVirtualSite(index)
    ]
    top_level_other: list[dict[str, Any]] = []
    for index in range(system.getNumForces()):
        if index == atm_index:
            continue
        force_xml = mm.XmlSerializer.serialize(system.getForce(index))
        top_level_other.append(
            {
                "index": index,
                "name": system.getForce(index).getName(),
                "sha256": hashlib.sha256(force_xml.encode("utf-8")).hexdigest(),
            }
        )
    nested: list[dict[str, Any]] = []
    for index in range(atm_force.getNumForces()):
        force = atm_force.getForce(index)
        force_xml = mm.XmlSerializer.serialize(force)
        nested.append(
            {
                "index": index,
                "name": force.getName(),
                "sha256": hashlib.sha256(force_xml.encode("utf-8")).hexdigest(),
            }
        )
    globals_ = [
        [
            atm_force.getGlobalParameterName(index),
            float(atm_force.getGlobalParameterDefaultValue(index)),
        ]
        for index in range(atm_force.getNumGlobalParameters())
    ]
    derivatives = [
        atm_force.getEnergyParameterDerivativeName(index)
        for index in range(atm_force.getNumEnergyParameterDerivatives())
    ]
    box = [_vector_nm(value) for value in system.getDefaultPeriodicBoxVectors()]
    return {
        "n_particles": int(system.getNumParticles()),
        "n_constraints": int(system.getNumConstraints()),
        "n_virtual_sites": len(virtual_indices),
        "virtual_site_indices_sha256": canonical_digest(virtual_indices),
        "masses_sha256": masses,
        "constraints_sha256": constraints,
        "box_vectors_nm": box,
        "n_top_level_forces": int(system.getNumForces()),
        "top_level_non_atm": top_level_other,
        "atm": {
            "index": int(atm_index),
            "name": atm_force.getName(),
            "force_group": int(atm_force.getForceGroup()),
            "energy_function": atm_force.getEnergyFunction(),
            "global_parameters": globals_,
            "energy_parameter_derivatives": derivatives,
            "n_particles": int(atm_force.getNumParticles()),
            "transformations_sha256": row_digest(
                _transformation_rows(atm_force)
            ),
            "nested_forces": nested,
        },
    }


def compact_force_report(report: dict[str, Any]) -> dict[str, Any]:
    compact = json.loads(json.dumps(report))
    selection = compact.get("selection") or {}
    water_indices = selection.pop("water_oxygen_indices_sha_input", [])
    selection["water_oxygen_indices_sha256"] = canonical_digest(water_indices)
    selection["water_oxygen_index_first"] = (
        water_indices[0] if water_indices else None
    )
    selection["water_oxygen_index_last"] = (
        water_indices[-1] if water_indices else None
    )
    for group in compact.get("interaction_groups") or []:
        second = group.pop("second", [])
        group["second_count"] = len(second)
        group["second_sha256"] = canonical_digest(second)
    return compact


def keeper_cell_dir(cell_id: str) -> Path:
    return KEEPER_ROOT / cell_id


def keeper_cell_result_path(cell_id: str) -> Path:
    return keeper_cell_dir(cell_id) / "keeper_result.json"


def keeper_nested_xml_path(cell_id: str) -> Path:
    return keeper_cell_dir(cell_id) / "atm_nested_ghost_system.xml"


def _write_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S %z')} {message}"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def _keeper_invariants(
    before: dict[str, Any],
    after: dict[str, Any],
) -> dict[str, bool]:
    before_atm = before["atm"]
    after_atm = after["atm"]
    return {
        "particle_count_unchanged": (
            before["n_particles"] == after["n_particles"]
        ),
        "masses_unchanged": (
            before["masses_sha256"] == after["masses_sha256"]
        ),
        "constraints_unchanged": (
            before["n_constraints"] == after["n_constraints"]
            and before["constraints_sha256"] == after["constraints_sha256"]
        ),
        "virtual_sites_unchanged": (
            before["n_virtual_sites"] == after["n_virtual_sites"]
            and before["virtual_site_indices_sha256"]
            == after["virtual_site_indices_sha256"]
        ),
        "box_unchanged": before["box_vectors_nm"] == after["box_vectors_nm"],
        "top_level_forces_unchanged": (
            before["n_top_level_forces"] == after["n_top_level_forces"]
            and before["top_level_non_atm"] == after["top_level_non_atm"]
        ),
        "atm_core_unchanged": all(
            before_atm[key] == after_atm[key]
            for key in (
                "index",
                "name",
                "force_group",
                "energy_function",
                "global_parameters",
                "energy_parameter_derivatives",
                "n_particles",
                "transformations_sha256",
            )
        ),
        "existing_nested_forces_unchanged": (
            after_atm["nested_forces"][: len(before_atm["nested_forces"])]
            == before_atm["nested_forces"]
        ),
        "one_nested_force_appended": (
            len(after_atm["nested_forces"])
            == len(before_atm["nested_forces"]) + 1
            and after_atm["nested_forces"][-1]["name"]
            == dg.ATM_NESTED_GHOST_FORCE_NAME
        ),
    }


def run_keeper_cell(
    inventory: dict[str, Any], cell_id: str
) -> dict[str, Any]:
    import openmm as mm
    from openmm.app import PDBFile

    protocol, cell = _verify_worker_inputs(inventory, cell_id)
    output_dir = keeper_cell_dir(cell_id)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True)
    log_path = output_dir / "keeper.log"
    started = time.monotonic()
    _write_log(log_path, f"START {cell_id} context_count=0")

    source_xml_path = verify_artifact(cell["artifacts"]["source_system"])
    pdb_path = verify_artifact(cell["artifacts"]["topology_pdb"])
    positions_path = verify_artifact(cell["artifacts"]["positions"])
    system = mm.XmlSerializer.deserialize(
        source_xml_path.read_text(encoding="utf-8")
    )
    pdb = PDBFile(str(pdb_path))
    positions_nm = np.load(positions_path)
    if positions_nm.shape != (system.getNumParticles(), 3):
        raise dg.DynamicGhostError(
            f"{cell_id}: position shape {positions_nm.shape} does not match System"
        )
    if pdb.topology.getNumAtoms() != system.getNumParticles():
        raise dg.DynamicGhostError(f"{cell_id}: topology particle-count drift")

    _write_log(log_path, "snapshot source declarations")
    before = _system_declaration(system)
    if any(
        row["name"] in {dg.GHOST_FORCE_NAME, dg.ATM_NESTED_GHOST_FORCE_NAME}
        for row in before["top_level_non_atm"]
        + before["atm"]["nested_forces"]
    ):
        raise dg.DynamicGhostError(f"{cell_id}: source already contains a ghost")
    source_transformations = dg.inspect_atm_coordinate_transformations(
        system,
        pdb.topology,
        positions_nm=positions_nm,
        protocol=protocol,
    )

    _write_log(log_path, "append one ATM-nested ghost and recheck declarations")
    build_report = dg.build_atm_nested_dynamic_ghost_force(
        system,
        pdb.topology,
        protocol,
    )
    after = _system_declaration(system)
    invariant_checks = _keeper_invariants(before, after)
    invariant_checks["charge_unchanged"] = (
        build_report["net_charge_delta_e"] == 0.0
    )
    invariant_checks["builder_scope"] = (
        build_report["n_top_level_forces_before"]
        == build_report["n_top_level_forces_after"]
        and build_report["n_atm_nested_forces_after"]
        == build_report["n_atm_nested_forces_before"] + 1
    )
    if not all(invariant_checks.values()):
        raise dg.DynamicGhostError(
            f"{cell_id}: KEEPER invariant failure: {invariant_checks}"
        )

    _write_log(log_path, "serialize and inspect round-trip without a Context")
    nested_xml = mm.XmlSerializer.serialize(system)
    nested_path = keeper_nested_xml_path(cell_id)
    atomic_write_text(nested_path, nested_xml)
    nested_artifact = artifact_record(nested_path)
    del system
    gc.collect()

    roundtrip = mm.XmlSerializer.deserialize(nested_xml)
    roundtrip_report = dg.inspect_atm_nested_dynamic_ghost_force(roundtrip)
    roundtrip_transformations = dg.inspect_atm_coordinate_transformations(
        roundtrip,
        pdb.topology,
        positions_nm=positions_nm,
        protocol=protocol,
    )
    roundtrip_declaration = _system_declaration(roundtrip)
    serialization_checks = {
        "force_contract_roundtrip": (
            compact_force_report(roundtrip_report)
            == {
                key: value
                for key, value in compact_force_report(build_report).items()
                if key
                in {
                    "atm_force_index",
                    "nested_force_index",
                    "force_name",
                    "force_group",
                    "energy_expression",
                    "global_parameters",
                    "per_particle_parameters",
                    "energy_parameter_derivatives",
                    "n_particles",
                    "nonbonded_method",
                    "cutoff_nm",
                    "long_range_correction",
                    "switching_function",
                    "interaction_groups",
                }
            }
        ),
        "transformations_roundtrip": (
            source_transformations == roundtrip_transformations
        ),
        "system_declaration_roundtrip": after == roundtrip_declaration,
        "serialized_artifact_hash": (
            sha256_file(nested_path) == nested_artifact["sha256"]
        ),
    }
    if not all(serialization_checks.values()):
        raise dg.DynamicGhostError(
            f"{cell_id}: serialization checks failed: {serialization_checks}"
        )

    result = {
        "schema": "updd_atm_nested_dynamic_ghost_p0_keeper_cell_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "PASS",
        "cell_id": cell_id,
        "seed": cell["seed"],
        "leg": cell["leg"],
        "inventory_digest": inventory["inventory_digest"],
        "context_count": 0,
        "md_steps": 0,
        "minimization_steps": 0,
        "gpu_used": False,
        "source_declaration": before,
        "nested_declaration": after,
        "transformations": source_transformations,
        "force_contract": compact_force_report(build_report),
        "invariant_checks": invariant_checks,
        "serialization_checks": serialization_checks,
        "nested_system": nested_artifact,
        "elapsed_s": time.monotonic() - started,
    }
    atomic_write_json(keeper_cell_result_path(cell_id), result)
    _write_log(log_path, f"COMPLETE {cell_id} status=PASS")
    del roundtrip, pdb, positions_nm, nested_xml
    gc.collect()
    return result


def _archive_directory(path: Path, reason: str) -> Path | None:
    if not path.exists():
        return None
    stamp = time.strftime("%Y%m%d_%H%M%S")
    target = path.parent / "_archive" / f"{stamp}_{path.name}"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(path), str(target))
    atomic_write_json(
        target / "archive_reason.json",
        {
            "archived": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "reason": reason,
        },
    )
    return target


def _keeper_command(cell_id: str, digest: str) -> list[str]:
    return [
        str(EXPECTED_ATM_PYTHON),
        str(RUNNER_PATH),
        "--keeper-cell",
        cell_id,
        "--expected-inventory-digest",
        digest,
    ]


def run_keeper_preflight(inventory: dict[str, Any]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for cell in inventory["cells"]:
        cell_id = cell["cell_id"]
        existing = keeper_cell_result_path(cell_id)
        if existing.is_file():
            row = json.loads(existing.read_text(encoding="utf-8"))
            if (
                row.get("status") == "PASS"
                and row.get("inventory_digest") == inventory["inventory_digest"]
                and verify_artifact(row["nested_system"])
            ):
                rows.append(row)
                continue
        _archive_directory(
            keeper_cell_dir(cell_id),
            "incomplete or drifted KEEPER cell",
        )
        completed = subprocess.run(
            _keeper_command(cell_id, inventory["inventory_digest"]),
            cwd=REPO_ROOT,
            check=False,
        )
        if completed.returncode != 0 or not existing.is_file():
            raise RuntimeError(
                f"{cell_id}: KEEPER worker exited {completed.returncode}"
            )
        row = json.loads(existing.read_text(encoding="utf-8"))
        if (
            row.get("status") != "PASS"
            or row.get("inventory_digest") != inventory["inventory_digest"]
        ):
            raise RuntimeError(f"{cell_id}: KEEPER worker did not pass")
        verify_artifact(row["nested_system"])
        rows.append(row)

    payload = {
        "schema": "updd_atm_nested_dynamic_ghost_p0_keeper_preflight_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "PASS",
        "inventory_digest": inventory["inventory_digest"],
        "context_count": 0,
        "md_steps": 0,
        "minimization_steps": 0,
        "gpu_used": False,
        "n_expected": len(inventory["cells"]),
        "n_completed": len(rows),
        "cells": [
            {
                "cell_id": row["cell_id"],
                "seed": row["seed"],
                "leg": row["leg"],
                "status": row["status"],
                "elapsed_s": row["elapsed_s"],
                "keeper_result": artifact_record(
                    keeper_cell_result_path(row["cell_id"])
                ),
                "nested_system": row["nested_system"],
            }
            for row in rows
        ],
    }
    payload["keeper_digest"] = canonical_digest(
        {
            key: value
            for key, value in payload.items()
            if key not in {"generated", "keeper_digest"}
        }
    )
    atomic_write_json(KEEPER_PREFLIGHT, payload)
    return payload


def validate_keeper(inventory: dict[str, Any]) -> dict[str, Any]:
    if not KEEPER_PREFLIGHT.is_file():
        raise FileNotFoundError(
            f"missing {KEEPER_PREFLIGHT}; run --keeper-preflight first"
        )
    keeper = json.loads(KEEPER_PREFLIGHT.read_text(encoding="utf-8"))
    expected_digest = canonical_digest(
        {
            key: value
            for key, value in keeper.items()
            if key not in {"generated", "keeper_digest"}
        }
    )
    if (
        keeper.get("status") != "PASS"
        or keeper.get("context_count") != 0
        or keeper.get("n_completed") != len(inventory["cells"])
        or keeper.get("inventory_digest") != inventory["inventory_digest"]
        or keeper.get("keeper_digest") != expected_digest
    ):
        raise ValueError("P0 KEEPER preflight is invalid")
    for row in keeper["cells"]:
        verify_artifact(row["keeper_result"])
        verify_artifact(row["nested_system"])
    return keeper


def _set_state(
    context: Any,
    state: dict[str, Any],
    *,
    ghost_g: float | None,
) -> None:
    context.setParameter("Lambda1", state["lambda1"])
    context.setParameter("Lambda2", state["lambda2"])
    context.setParameter("Alpha", state["alpha_kcal_inv"] / KCAL_TO_KJ)
    context.setParameter("Uh", state["uh_kcal_mol"] * KCAL_TO_KJ)
    context.setParameter("W0", state["w0_kcal_mol"] * KCAL_TO_KJ)
    context.setParameter("Umax", state["umax_kcal_mol"] * KCAL_TO_KJ)
    context.setParameter("Ubcore", state["ubcore_kcal_mol"] * KCAL_TO_KJ)
    context.setParameter("Acore", state["acore"])
    context.setParameter("Direction", float(state["direction"]))
    context.setParameter("UOffset", 0.0)
    if ghost_g is not None:
        context.setParameter(dg.GHOST_GLOBAL_PARAMETER, float(ghost_g))


def _read_context(
    context: Any,
    atm_force: Any,
    *,
    require_dvdg: bool,
) -> dict[str, Any]:
    from openmm import unit

    started = time.monotonic()
    state = context.getState(
        getEnergy=True,
        getForces=True,
        getParameterDerivatives=require_dvdg,
    )
    energy = float(
        state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    )
    forces = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=float,
    )
    derivatives = (
        {
            str(name): float(value)
            for name, value in dict(
                state.getEnergyParameterDerivatives()
            ).items()
        }
        if require_dvdg
        else {}
    )
    u1, u0, atm_energy = atm_force.getPerturbationEnergy(context)
    raw = {
        "u0_kj_mol": float(u0.value_in_unit(unit.kilojoule_per_mole)),
        "u1_kj_mol": float(u1.value_in_unit(unit.kilojoule_per_mole)),
        "atm_energy_kj_mol": float(
            atm_energy.value_in_unit(unit.kilojoule_per_mole)
        ),
    }
    finite = bool(
        math.isfinite(energy)
        and np.isfinite(forces).all()
        and all(math.isfinite(value) for value in raw.values())
        and (
            not require_dvdg
            or (
                dg.GHOST_GLOBAL_PARAMETER in derivatives
                and math.isfinite(
                    derivatives[dg.GHOST_GLOBAL_PARAMETER]
                )
            )
        )
    )
    return {
        "energy_kj_mol": energy,
        "forces": forces,
        "derivatives_kj_mol": derivatives,
        "raw": raw,
        "finite": finite,
        "elapsed_s": time.monotonic() - started,
    }


def _evaluation_summary(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "energy_kj_mol": value["energy_kj_mol"],
        "max_abs_force_component_kj_mol_nm": float(
            np.max(np.abs(value["forces"]))
        ),
        "dVdg_kj_mol": value["derivatives_kj_mol"].get(
            dg.GHOST_GLOBAL_PARAMETER
        ),
        "raw": value["raw"],
        "finite": value["finite"],
        "elapsed_s": value["elapsed_s"],
    }


def _minimum_image(vector: np.ndarray, box: np.ndarray) -> np.ndarray:
    fractional = vector @ np.linalg.inv(box)
    fractional -= np.round(fractional)
    return fractional @ box


def _transformed_ring_positions(
    system: Any,
    topology: Any,
    positions_nm: np.ndarray,
) -> tuple[dg.GhostSelection, np.ndarray, np.ndarray]:
    import openmm as mm

    selection = dg.select_ghost_atoms(topology)
    _atm_index, atm_force = dg.find_atm_force(system)
    stored = positions_nm[np.asarray(selection.ring_indices, dtype=int)]
    transformed: list[np.ndarray] = []
    for index in selection.ring_indices:
        transformation = atm_force.getParticleTransformation(int(index))
        if not isinstance(transformation, mm.ParticleOffsetDisplacement):
            raise dg.DynamicGhostError(
                f"ring atom {index} lacks ParticleOffsetDisplacement"
            )
        transformed.append(
            positions_nm[index]
            + positions_nm[transformation.getDestinationParticle1()]
            - positions_nm[transformation.getOriginParticle1()]
        )
    return selection, stored, np.asarray(transformed, dtype=float)


def _zero_fixed_displacement(atm_force: Any, index: int) -> bool:
    import openmm as mm

    transformation = atm_force.getParticleTransformation(int(index))
    if not isinstance(transformation, mm.FixedDisplacement):
        return False
    return (
        _vector_nm(transformation.getFixedDisplacement0()) == [0.0, 0.0, 0.0]
        and _vector_nm(transformation.getFixedDisplacement1())
        == [0.0, 0.0, 0.0]
    )


def build_probe_coordinates(
    system: Any,
    topology: Any,
    positions_nm: np.ndarray,
    protocol: dict[str, Any],
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    selection, stored, transformed = _transformed_ring_positions(
        system,
        topology,
        positions_nm,
    )
    box = np.asarray(
        [_vector_nm(value) for value in system.getDefaultPeriodicBoxVectors()],
        dtype=float,
    )
    target_name = protocol["probe"]["target_atom_name"]
    target_offset = [
        index
        for index, row in enumerate(selection.ring_atoms)
        if row["name"] == target_name
    ]
    if len(target_offset) != 1:
        raise dg.DynamicGhostError("probe target NE1 selection drifted")
    target_offset = target_offset[0]
    stored_target = stored[target_offset]
    unwrapped = np.asarray(
        [
            stored_target + _minimum_image(row - stored_target, box)
            for row in stored
        ],
        dtype=float,
    )
    outward = stored_target - np.mean(unwrapped, axis=0)
    outward_norm = float(np.linalg.norm(outward))
    if not math.isfinite(outward_norm) or outward_norm <= 0.0:
        raise dg.DynamicGhostError("invalid ring-centroid outward vector")
    outward /= outward_norm

    _atm_index, atm_force = dg.find_atm_force(system)
    source_clearance_nm = float(protocol["probe"]["source_clearance_nm"])
    candidates: list[
        tuple[
            float,
            int,
            tuple[int, ...],
            np.ndarray,
            float,
            float,
        ]
    ] = []
    for residue in topology.residues():
        if residue.name not in dg.WATER_RESIDUE_NAMES:
            continue
        atoms = tuple(residue.atoms())
        if len(atoms) != protocol["probe"]["water_residue_atom_count"]:
            continue
        oxygens = [
            atom
            for atom in atoms
            if atom.name in dg.WATER_OXYGEN_NAMES
        ]
        hydrogens = [
            atom
            for atom in atoms
            if atom.element is not None and atom.element.symbol == "H"
        ]
        if len(oxygens) != 1 or len(hydrogens) != 2:
            continue
        indices = tuple(int(atom.index) for atom in atoms)
        if not all(
            _zero_fixed_displacement(atm_force, index) for index in indices
        ):
            continue
        oxygen_index = int(oxygens[0].index)
        oxygen = positions_nm[oxygen_index]
        stored_clearance = float(
            min(
                np.linalg.norm(_minimum_image(oxygen - row, box))
                for row in stored
            )
        )
        transformed_clearance = float(
            min(
                np.linalg.norm(_minimum_image(oxygen - row, box))
                for row in transformed
            )
        )
        if (
            stored_clearance < source_clearance_nm
            or transformed_clearance < source_clearance_nm
        ):
            continue
        h_vectors = [
            _minimum_image(positions_nm[atom.index] - oxygen, box)
            for atom in hydrogens
        ]
        bisector = np.mean(np.asarray(h_vectors), axis=0)
        norm = float(np.linalg.norm(bisector))
        if not math.isfinite(norm) or norm <= 0.0:
            continue
        bisector /= norm
        candidates.append(
            (
                float(np.dot(bisector, outward)),
                oxygen_index,
                indices,
                bisector,
                stored_clearance,
                transformed_clearance,
            )
        )
    if not candidates:
        raise dg.DynamicGhostError("no three-atom zero-displacement probe water")
    (
        orientation_dot,
        oxygen_index,
        water_indices,
        bisector,
        selected_stored_clearance,
        selected_transformed_clearance,
    ) = max(
        candidates,
        key=lambda row: (row[0], -row[1]),
    )

    cutoff = float(
        protocol["ghost"]["reference_pair_parameters"]["NE1_water_O"][
            "wca_cutoff_nm"
        ]
    )
    probe_distance = (
        float(protocol["probe"]["distance_fraction_of_target_wca_cutoff"])
        * cutoff
    )
    target_positions = {
        "stored_probe": stored[target_offset] + probe_distance * outward,
        "transformed_probe": (
            transformed[target_offset] + probe_distance * outward
        ),
    }
    coordinate_sets: dict[str, np.ndarray] = {}
    for name, target in target_positions.items():
        coordinates = np.array(positions_nm, copy=True)
        translation = target - positions_nm[oxygen_index]
        coordinates[np.asarray(water_indices, dtype=int)] += translation
        coordinate_sets[name] = coordinates

    opposite_distances = {
        "stored_probe_to_transformed_ring_min_nm": float(
            np.min(
                np.linalg.norm(
                    np.asarray(
                        [
                            _minimum_image(
                                target_positions["stored_probe"] - row,
                                box,
                            )
                            for row in transformed
                        ]
                    ),
                    axis=1,
                )
            )
        ),
        "transformed_probe_to_stored_ring_min_nm": float(
            np.min(
                np.linalg.norm(
                    np.asarray(
                        [
                            _minimum_image(
                                target_positions["transformed_probe"] - row,
                                box,
                            )
                            for row in stored
                        ]
                    ),
                    axis=1,
                )
            )
        ),
    }
    report = {
        "target_atom_name": target_name,
        "target_atom_index": int(selection.ring_indices[target_offset]),
        "outward_unit_vector": [float(value) for value in outward],
        "probe_distance_nm": probe_distance,
        "probe_distance_fraction_of_wca_cutoff": float(
            protocol["probe"]["distance_fraction_of_target_wca_cutoff"]
        ),
        "target_wca_cutoff_nm": cutoff,
        "water_oxygen_index": oxygen_index,
        "water_atom_indices": list(water_indices),
        "water_orientation_dot": orientation_dot,
        "water_bisector_unit_vector": [float(value) for value in bisector],
        "candidate_water_count": len(candidates),
        "source_clearance_threshold_nm": source_clearance_nm,
        "selected_source_stored_clearance_nm": selected_stored_clearance,
        "selected_source_transformed_clearance_nm": (
            selected_transformed_clearance
        ),
        "isolation_metric": protocol["probe"]["isolation_metric"],
        "opposite_site_distances": opposite_distances,
    }
    return coordinate_sets, report


def _evaluate_original(
    system: Any,
    positions_nm: np.ndarray,
    states: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    import openmm as mm
    from openmm import unit

    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(
        system,
        integrator,
        mm.Platform.getPlatformByName("Reference"),
    )
    context.setPositions(positions_nm * unit.nanometer)
    _atm_index, atm_force = dg.find_atm_force(system)
    results: dict[str, dict[str, Any]] = {}
    for name in ("u0_dplus", "u1_dminus"):
        _set_state(context, states[name], ghost_g=None)
        results[name] = _read_context(
            context,
            atm_force,
            require_dvdg=False,
        )
    del context, integrator
    gc.collect()
    return results


def _evaluate_nested(
    system: Any,
    topology: Any,
    positions_nm: np.ndarray,
    states: dict[str, Any],
    protocol: dict[str, Any],
) -> tuple[
    dict[str, dict[str, dict[str, Any]]],
    dict[str, dict[str, dict[str, Any]]],
    dict[str, Any],
]:
    import openmm as mm
    from openmm import unit

    coordinate_sets, probe_report = build_probe_coordinates(
        system,
        topology,
        positions_nm,
        protocol,
    )
    all_coordinates = {"source": positions_nm, **coordinate_sets}
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(
        system,
        integrator,
        mm.Platform.getPlatformByName("Reference"),
    )
    _atm_index, atm_force = dg.find_atm_force(system)
    results: dict[str, dict[str, dict[str, Any]]] = {}
    for coordinate_name, coordinates in all_coordinates.items():
        context.setPositions(coordinates * unit.nanometer)
        coordinate_results: dict[str, dict[str, Any]] = {}
        for state_name in ("u0_dplus", "u1_dminus"):
            state_results: dict[str, Any] = {}
            for g in (0.0, 1.0):
                _set_state(context, states[state_name], ghost_g=g)
                state_results[f"g{int(g)}"] = _read_context(
                    context,
                    atm_force,
                    require_dvdg=True,
                )
            coordinate_results[state_name] = state_results
        results[coordinate_name] = coordinate_results
    summarized = {
        coordinate: {
            state_name: {
                g: _evaluation_summary(value)
                for g, value in g_rows.items()
            }
            for state_name, g_rows in state_rows.items()
        }
        for coordinate, state_rows in results.items()
    }
    del context, integrator
    gc.collect()
    return results, summarized, probe_report


def _parity_record(
    observed: dict[str, Any],
    expected: dict[str, Any],
    gates: dict[str, Any],
) -> dict[str, Any]:
    energy_delta = abs(
        float(observed["energy_kj_mol"]) - float(expected["energy_kj_mol"])
    )
    force_delta = float(
        np.max(np.abs(observed["forces"] - expected["forces"]))
    )
    raw_u0_delta = abs(
        observed["raw"]["u0_kj_mol"] - expected["raw"]["u0_kj_mol"]
    )
    raw_u1_delta = abs(
        observed["raw"]["u1_kj_mol"] - expected["raw"]["u1_kj_mol"]
    )
    passed = bool(
        observed["finite"]
        and expected["finite"]
        and energy_delta <= gates["energy_parity_abs_kj_mol"]
        and force_delta
        <= gates["force_parity_max_component_kj_mol_nm"]
        and raw_u0_delta <= gates["energy_parity_abs_kj_mol"]
        and raw_u1_delta <= gates["energy_parity_abs_kj_mol"]
    )
    return {
        "abs_energy_delta_kj_mol": energy_delta,
        "max_abs_force_component_delta_kj_mol_nm": force_delta,
        "abs_raw_u0_delta_kj_mol": raw_u0_delta,
        "abs_raw_u1_delta_kj_mol": raw_u1_delta,
        "passed": passed,
    }


def _isolated_record(
    g0: dict[str, Any],
    g1: dict[str, Any],
) -> dict[str, Any]:
    forces = g1["forces"] - g0["forces"]
    return {
        "total_energy_delta_kj_mol": (
            g1["energy_kj_mol"] - g0["energy_kj_mol"]
        ),
        "raw_u0_delta_kj_mol": (
            g1["raw"]["u0_kj_mol"] - g0["raw"]["u0_kj_mol"]
        ),
        "raw_u1_delta_kj_mol": (
            g1["raw"]["u1_kj_mol"] - g0["raw"]["u1_kj_mol"]
        ),
        "max_abs_force_component_kj_mol_nm": float(
            np.max(np.abs(forces))
        ),
        "force_balance_vector_kj_mol_nm": [
            float(value) for value in np.sum(forces, axis=0)
        ],
        "forces": forces,
        "finite": bool(
            g0["finite"]
            and g1["finite"]
            and np.isfinite(forces).all()
        ),
    }


def _baseline_corrected_isolated(
    probe: dict[str, Any],
    source: dict[str, Any],
) -> dict[str, Any]:
    forces = probe["forces"] - source["forces"]
    return {
        "total_energy_delta_kj_mol": (
            probe["total_energy_delta_kj_mol"]
            - source["total_energy_delta_kj_mol"]
        ),
        "raw_u0_delta_kj_mol": (
            probe["raw_u0_delta_kj_mol"]
            - source["raw_u0_delta_kj_mol"]
        ),
        "raw_u1_delta_kj_mol": (
            probe["raw_u1_delta_kj_mol"]
            - source["raw_u1_delta_kj_mol"]
        ),
        "max_abs_force_component_kj_mol_nm": float(
            np.max(np.abs(forces))
        ),
        "force_balance_vector_kj_mol_nm": [
            float(value) for value in np.sum(forces, axis=0)
        ],
        "forces": forces,
        "finite": bool(
            probe["finite"]
            and source["finite"]
            and np.isfinite(forces).all()
        ),
    }


def _source_gate(
    nested: dict[str, dict[str, dict[str, Any]]],
    gates: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, bool]]:
    records: dict[str, Any] = {}
    checks: dict[str, bool] = {}
    for state_name in ("u0_dplus", "u1_dminus"):
        isolated = _isolated_record(
            nested["source"][state_name]["g0"],
            nested["source"][state_name]["g1"],
        )
        primary = (
            isolated["raw_u0_delta_kj_mol"]
            if state_name == "u0_dplus"
            else isolated["raw_u1_delta_kj_mol"]
        )
        record = {
            key: value
            for key, value in isolated.items()
            if key != "forces"
        }
        record["primary_endpoint_isolated_energy_kj_mol"] = primary
        records[state_name] = record
        checks[f"{state_name}_finite"] = isolated["finite"]
        checks[f"{state_name}_energy_bounded"] = (
            abs(primary)
            <= gates["source_isolated_energy_abs_max_kj_mol"]
        )
        checks[f"{state_name}_force_bounded"] = (
            isolated["max_abs_force_component_kj_mol_nm"]
            <= gates["source_isolated_force_max_component_kj_mol_nm"]
        )
    return records, checks


def _probe_gate(
    nested: dict[str, dict[str, dict[str, Any]]],
    probe_report: dict[str, Any],
    gates: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, bool]]:
    outward = np.asarray(probe_report["outward_unit_vector"], dtype=float)
    oxygen = int(probe_report["water_oxygen_index"])
    records: dict[str, Any] = {}
    checks: dict[str, bool] = {}
    endpoint_map = {
        "stored_probe": ("u0_dplus", "u1_dminus", "raw_u0_delta_kj_mol"),
        "transformed_probe": (
            "u1_dminus",
            "u0_dplus",
            "raw_u1_delta_kj_mol",
        ),
    }
    for coordinate, (
        intended_state,
        opposite_state,
        intended_raw_key,
    ) in endpoint_map.items():
        absolute_intended = _isolated_record(
            nested[coordinate][intended_state]["g0"],
            nested[coordinate][intended_state]["g1"],
        )
        absolute_opposite = _isolated_record(
            nested[coordinate][opposite_state]["g0"],
            nested[coordinate][opposite_state]["g1"],
        )
        source_intended = _isolated_record(
            nested["source"][intended_state]["g0"],
            nested["source"][intended_state]["g1"],
        )
        source_opposite = _isolated_record(
            nested["source"][opposite_state]["g0"],
            nested["source"][opposite_state]["g1"],
        )
        intended = _baseline_corrected_isolated(
            absolute_intended,
            source_intended,
        )
        opposite = _baseline_corrected_isolated(
            absolute_opposite,
            source_opposite,
        )
        opposite_raw_key = (
            "raw_u1_delta_kj_mol"
            if intended_raw_key == "raw_u0_delta_kj_mol"
            else "raw_u0_delta_kj_mol"
        )
        intended_energy = float(intended[intended_raw_key])
        opposite_energy = float(opposite[opposite_raw_key])
        radial_force = float(np.dot(intended["forces"][oxygen], outward))
        balance_max = float(
            np.max(np.abs(intended["force_balance_vector_kj_mol_nm"]))
        )
        intended_total_consistency = abs(
            intended["total_energy_delta_kj_mol"] - intended_energy
        )
        opposite_total_consistency = abs(
            opposite["total_energy_delta_kj_mol"] - opposite_energy
        )
        records[coordinate] = {
            "intended_state": intended_state,
            "opposite_state": opposite_state,
            "intended_isolated_energy_kj_mol": intended_energy,
            "opposite_isolated_energy_kj_mol": opposite_energy,
            "intended_radial_oxygen_force_kj_mol_nm": radial_force,
            "intended_max_abs_force_component_kj_mol_nm": intended[
                "max_abs_force_component_kj_mol_nm"
            ],
            "opposite_max_abs_force_component_kj_mol_nm": opposite[
                "max_abs_force_component_kj_mol_nm"
            ],
            "intended_force_balance_vector_kj_mol_nm": intended[
                "force_balance_vector_kj_mol_nm"
            ],
            "intended_force_balance_max_component_kj_mol_nm": balance_max,
            "intended_total_raw_consistency_kj_mol": (
                intended_total_consistency
            ),
            "opposite_total_raw_consistency_kj_mol": (
                opposite_total_consistency
            ),
            "absolute_probe_g_difference": {
                "intended_energy_kj_mol": absolute_intended[
                    intended_raw_key
                ],
                "opposite_energy_kj_mol": absolute_opposite[
                    opposite_raw_key
                ],
            },
            "source_g_difference_baseline": {
                "intended_energy_kj_mol": source_intended[
                    intended_raw_key
                ],
                "opposite_energy_kj_mol": source_opposite[
                    opposite_raw_key
                ],
            },
            "threshold_metric": (
                "endpoint_matched_source_baseline_double_difference"
            ),
            "finite": intended["finite"] and opposite["finite"],
        }
        checks[f"{coordinate}_finite"] = (
            intended["finite"] and opposite["finite"]
        )
        checks[f"{coordinate}_intended_energy"] = (
            intended_energy >= gates["probe_intended_energy_min_kj_mol"]
        )
        checks[f"{coordinate}_opposite_energy"] = (
            abs(opposite_energy)
            <= gates["probe_opposite_energy_abs_max_kj_mol"]
        )
        checks[f"{coordinate}_radial_force"] = (
            radial_force
            >= gates["probe_intended_radial_force_min_kj_mol_nm"]
        )
        checks[f"{coordinate}_opposite_force"] = (
            opposite["max_abs_force_component_kj_mol_nm"]
            <= gates["probe_opposite_force_max_component_kj_mol_nm"]
        )
        checks[f"{coordinate}_force_balance"] = (
            balance_max
            <= gates["probe_force_balance_max_component_kj_mol_nm"]
        )
        checks[f"{coordinate}_endpoint_identity"] = (
            intended_total_consistency
            <= gates["energy_parity_abs_kj_mol"]
            and opposite_total_consistency
            <= gates["energy_parity_abs_kj_mol"]
        )
    return records, checks


def reference_cell_dir(cell_id: str) -> Path:
    return REFERENCE_ROOT / cell_id


def reference_cell_result_path(cell_id: str) -> Path:
    return reference_cell_dir(cell_id) / "cell_result.json"


def run_reference_cell(
    inventory: dict[str, Any],
    keeper: dict[str, Any],
    cell_id: str,
) -> dict[str, Any]:
    import openmm as mm
    from openmm.app import PDBFile

    protocol, cell = _verify_worker_inputs(inventory, cell_id)
    output_dir = reference_cell_dir(cell_id)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True)
    log_path = output_dir / "worker.log"
    started = time.monotonic()
    _write_log(
        log_path,
        f"START {cell_id} platform=Reference md=0 minimization=0",
    )

    keeper_rows = [
        row for row in keeper["cells"] if row["cell_id"] == cell_id
    ]
    if len(keeper_rows) != 1:
        raise ValueError(f"{cell_id}: missing KEEPER row")
    nested_xml_path = verify_artifact(keeper_rows[0]["nested_system"])
    source_xml_path = verify_artifact(cell["artifacts"]["source_system"])
    topology_path = verify_artifact(cell["artifacts"]["topology_pdb"])
    positions_path = verify_artifact(cell["artifacts"]["positions"])
    positions_nm = np.load(positions_path)
    pdb = PDBFile(str(topology_path))

    _write_log(log_path, "evaluate unmodified source endpoints")
    original_system = mm.XmlSerializer.deserialize(
        source_xml_path.read_text(encoding="utf-8")
    )
    original = _evaluate_original(
        original_system,
        positions_nm,
        cell["states"],
    )
    del original_system
    gc.collect()

    _write_log(log_path, "evaluate serialized nested source and two probes")
    nested_system = mm.XmlSerializer.deserialize(
        nested_xml_path.read_text(encoding="utf-8")
    )
    dg.inspect_atm_nested_dynamic_ghost_force(nested_system)
    dg.inspect_atm_coordinate_transformations(
        nested_system,
        pdb.topology,
        positions_nm=positions_nm,
        protocol=protocol,
    )
    nested, evaluations, probe_report = _evaluate_nested(
        nested_system,
        pdb.topology,
        positions_nm,
        cell["states"],
        protocol,
    )
    del nested_system
    gc.collect()

    gates = protocol["gates"]
    parity = {
        state_name: _parity_record(
            nested["source"][state_name]["g0"],
            original[state_name],
            gates,
        )
        for state_name in ("u0_dplus", "u1_dminus")
    }
    source_isolation, source_checks = _source_gate(nested, gates)
    probes, probe_checks = _probe_gate(nested, probe_report, gates)
    finite_dvdg = all(
        value["dVdg_kj_mol"] is not None
        and math.isfinite(float(value["dVdg_kj_mol"]))
        for coordinate in evaluations.values()
        for state in coordinate.values()
        for value in state.values()
    )
    checks = {
        "all_g0_parities": all(row["passed"] for row in parity.values()),
        "finite_dvdg_all_nested_evaluations": finite_dvdg,
        **source_checks,
        **probe_checks,
    }
    passed = all(checks.values())

    result = {
        "schema": "updd_atm_nested_dynamic_ghost_p0_reference_cell_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "P0_PASS" if passed else "REJECT_GHOST",
        "cell_id": cell_id,
        "seed": cell["seed"],
        "leg": cell["leg"],
        "inventory_digest": inventory["inventory_digest"],
        "keeper_digest": keeper["keeper_digest"],
        "platform": {
            "name": "Reference",
            "openmm": mm.__version__,
            "python": sys.executable,
            "gpu_used": False,
            "md_steps": 0,
            "minimization_steps": 0,
            "system_variants": 2,
            "max_simultaneous_contexts": 1,
        },
        "source_artifacts": cell["artifacts"],
        "nested_system": keeper_rows[0]["nested_system"],
        "states": cell["states"],
        "original_evaluations": {
            name: _evaluation_summary(value)
            for name, value in original.items()
        },
        "nested_evaluations": evaluations,
        "g0_parity": parity,
        "source_isolation": source_isolation,
        "probe_construction": probe_report,
        "probe_isolation": probes,
        "checks": checks,
        "elapsed_s": time.monotonic() - started,
    }
    atomic_write_json(reference_cell_result_path(cell_id), result)
    _write_log(log_path, f"COMPLETE {cell_id} status={result['status']}")
    del original, nested, positions_nm, pdb
    gc.collect()
    return result


def _validate_completed_reference(
    inventory: dict[str, Any],
    keeper: dict[str, Any],
    cell_id: str,
) -> dict[str, Any] | None:
    path = reference_cell_result_path(cell_id)
    if not path.is_file():
        return None
    try:
        result = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    cell = frozen_cell(inventory, cell_id)
    if (
        result.get("cell_id") != cell_id
        or result.get("seed") != cell["seed"]
        or result.get("leg") != cell["leg"]
        or result.get("inventory_digest") != inventory["inventory_digest"]
        or result.get("keeper_digest") != keeper["keeper_digest"]
        or result.get("status") not in {"P0_PASS", "REJECT_GHOST"}
    ):
        return None
    return result


def _reference_command(cell_id: str, digest: str) -> list[str]:
    return [
        str(EXPECTED_ATM_PYTHON),
        str(RUNNER_PATH),
        "--worker-cell",
        cell_id,
        "--expected-inventory-digest",
        digest,
    ]


def _write_run_state(status: str, **extra: Any) -> None:
    atomic_write_json(
        RUN_STATE,
        {
            "status": status,
            "updated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "platform": "Reference",
            "md_steps": 0,
            "minimization_steps": 0,
            "gpu_used": False,
            **extra,
        },
    )


def write_summary(
    results: list[dict[str, Any]],
    inventory: dict[str, Any],
    keeper: dict[str, Any],
) -> dict[str, Any]:
    by_id = {row["cell_id"]: row for row in results}
    pending = [
        row["cell_id"]
        for row in inventory["cells"]
        if row["cell_id"] not in by_id
    ]
    rejected = [
        row["cell_id"]
        for row in results
        if row["status"] == "REJECT_GHOST"
    ]
    if rejected:
        status = "REJECT_GHOST"
    elif not pending and len(results) == len(inventory["cells"]):
        status = "P0_PASS"
    else:
        status = "RUNNING"
    elapsed = [float(row["elapsed_s"]) for row in results]
    payload = {
        "schema": "updd_atm_nested_dynamic_ghost_p0_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "inventory_digest": inventory["inventory_digest"],
        "keeper_digest": keeper["keeper_digest"],
        "platform": "Reference",
        "md_steps": 0,
        "minimization_steps": 0,
        "gpu_used": False,
        "n_expected": len(inventory["cells"]),
        "n_completed": len(results),
        "pending_cell_ids": pending,
        "rejected_cell_ids": rejected,
        "median_elapsed_s": median(elapsed) if elapsed else None,
        "eta_s": (
            median(elapsed) * len(pending)
            if elapsed and not rejected
            else None
        ),
        "cells": [
            {
                "cell_id": row["cell_id"],
                "seed": row["seed"],
                "leg": row["leg"],
                "status": row["status"],
                "elapsed_s": row["elapsed_s"],
                "stored_probe_energy_kj_mol": row["probe_isolation"][
                    "stored_probe"
                ]["intended_isolated_energy_kj_mol"],
                "transformed_probe_energy_kj_mol": row["probe_isolation"][
                    "transformed_probe"
                ]["intended_isolated_energy_kj_mol"],
                "max_opposite_energy_abs_kj_mol": max(
                    abs(
                        row["probe_isolation"]["stored_probe"][
                            "opposite_isolated_energy_kj_mol"
                        ]
                    ),
                    abs(
                        row["probe_isolation"]["transformed_probe"][
                            "opposite_isolated_energy_kj_mol"
                        ]
                    ),
                ),
                "result": artifact_record(
                    reference_cell_result_path(row["cell_id"])
                ),
            }
            for row in results
        ],
    }
    payload["summary_digest"] = canonical_digest(
        {
            key: value
            for key, value in payload.items()
            if key not in {"generated", "summary_digest"}
        }
    )
    atomic_write_json(RUNNER_SUMMARY, payload)

    lines = [
        "# ATM-Nested Dynamic Ghost P0 Reference Summary",
        "",
        f"Status: `{status}`",
        f"Completed: `{len(results)}/{len(inventory['cells'])}`",
        "Platform: `Reference`",
        "MD/minimization steps: `0/0`",
        "",
        "| leg | seed | status | stored probe E | transformed probe E | max opposite | elapsed min |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in payload["cells"]:
        lines.append(
            "| {leg} | {seed} | {status} | {stored:.6g} | "
            "{transformed:.6g} | {opposite:.3g} | {elapsed:.2f} |".format(
                leg=row["leg"],
                seed=row["seed"],
                status=row["status"],
                stored=row["stored_probe_energy_kj_mol"],
                transformed=row["transformed_probe_energy_kj_mol"],
                opposite=row["max_opposite_energy_abs_kj_mol"],
                elapsed=row["elapsed_s"] / 60.0,
            )
        )
    if pending:
        lines.extend(["", "Pending: " + ", ".join(f"`{x}`" for x in pending)])
    if rejected:
        lines.extend(
            ["", "Rejected: " + ", ".join(f"`{x}`" for x in rejected)]
        )
    lines.extend(
        [
            "",
            "This zero-dynamics audit does not estimate a free energy, establish "
            "delta-G neutrality, or authorize production.",
        ]
    )
    atomic_write_text(RUNNER_REPORT, "\n".join(lines) + "\n")
    return payload


def run_all(inventory: dict[str, Any], keeper: dict[str, Any]) -> int:
    results: list[dict[str, Any]] = []
    for cell in inventory["cells"]:
        completed = _validate_completed_reference(
            inventory,
            keeper,
            cell["cell_id"],
        )
        if completed is not None:
            results.append(completed)
    summary = write_summary(results, inventory, keeper)
    if summary["status"] == "REJECT_GHOST":
        _write_run_state(
            "REJECT_GHOST",
            n_completed=len(results),
            n_expected=len(inventory["cells"]),
        )
        return 2

    for cell in inventory["cells"]:
        cell_id = cell["cell_id"]
        if cell_id in {row["cell_id"] for row in results}:
            continue
        _archive_directory(
            reference_cell_dir(cell_id),
            "incomplete or drifted Reference P0 cell",
        )
        _write_run_state(
            "RUNNING",
            current_cell=cell_id,
            n_completed=len(results),
            n_expected=len(inventory["cells"]),
        )
        print(f"[P0] launch {cell_id} on Reference", flush=True)
        completed = subprocess.run(
            _reference_command(cell_id, inventory["inventory_digest"]),
            cwd=REPO_ROOT,
            check=False,
        )
        result = _validate_completed_reference(
            inventory,
            keeper,
            cell_id,
        )
        if result is None:
            _write_run_state(
                "INDETERMINATE",
                current_cell=cell_id,
                return_code=completed.returncode,
                error="worker produced no valid completed result",
            )
            raise RuntimeError(
                f"{cell_id}: worker exit {completed.returncode} "
                "without a valid result"
            )
        results.append(result)
        summary = write_summary(results, inventory, keeper)
        print(
            f"[P0] complete {cell_id}; "
            f"{summary['n_completed']}/{summary['n_expected']} "
            f"status={result['status']}",
            flush=True,
        )
        if result["status"] == "REJECT_GHOST" or completed.returncode == 2:
            _write_run_state(
                "REJECT_GHOST",
                current_cell=cell_id,
                n_completed=len(results),
                n_expected=len(inventory["cells"]),
            )
            return 2
        if completed.returncode != 0:
            _write_run_state(
                "INDETERMINATE",
                current_cell=cell_id,
                return_code=completed.returncode,
            )
            raise RuntimeError(
                f"{cell_id}: worker exited {completed.returncode}"
            )

    _write_run_state(
        "P0_PASS",
        n_completed=len(results),
        n_expected=len(inventory["cells"]),
    )
    return 0


def _worker_failure(
    root: Path,
    cell_id: str,
    exc: Exception,
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    atomic_write_json(
        root / "failure.json",
        {
            "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "cell_id": cell_id,
            "status": (
                "REJECT_GHOST"
                if isinstance(exc, dg.DynamicGhostError)
                else "INDETERMINATE"
            ),
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
        },
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-inventory", action="store_true")
    parser.add_argument("--keeper-preflight", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--keeper-cell", help=argparse.SUPPRESS)
    parser.add_argument("--worker-cell", help=argparse.SUPPRESS)
    parser.add_argument("--expected-inventory-digest", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    _require_atm_python()
    worker_modes = int(bool(args.keeper_cell)) + int(bool(args.worker_cell))
    public_modes = sum(
        bool(value)
        for value in (
            args.prepare_inventory,
            args.keeper_preflight,
            args.dry_run,
            args.analyze_only,
        )
    )
    if worker_modes > 1 or public_modes > 1 or (worker_modes and public_modes):
        parser.error("choose exactly one worker or public execution mode")
    if worker_modes and not args.expected_inventory_digest:
        parser.error("worker mode requires --expected-inventory-digest")
    if not worker_modes and args.expected_inventory_digest:
        parser.error("--expected-inventory-digest is worker-only")

    if args.prepare_inventory:
        if PREBUILD_INVENTORY.exists():
            raise FileExistsError(PREBUILD_INVENTORY)
        inventory = build_inventory()
        atomic_write_json(PREBUILD_INVENTORY, inventory)
        print(
            json.dumps(
                {
                    "status": "PREPARED",
                    "n_cells": len(inventory["cells"]),
                    "inventory_digest": inventory["inventory_digest"],
                }
            )
        )
        return 0

    if args.keeper_cell:
        inventory = load_inventory(args.expected_inventory_digest)
        try:
            run_keeper_cell(inventory, args.keeper_cell)
        except Exception as exc:
            _worker_failure(keeper_cell_dir(args.keeper_cell), args.keeper_cell, exc)
            raise
        return 0

    if args.worker_cell:
        inventory = load_inventory(args.expected_inventory_digest)
        keeper = validate_keeper(inventory)
        try:
            result = run_reference_cell(
                inventory,
                keeper,
                args.worker_cell,
            )
        except Exception as exc:
            _worker_failure(
                reference_cell_dir(args.worker_cell),
                args.worker_cell,
                exc,
            )
            raise
        return 0 if result["status"] == "P0_PASS" else 2

    inventory = verify_inventory()
    if args.keeper_preflight:
        result = run_keeper_preflight(inventory)
        print(
            json.dumps(
                {
                    "status": result["status"],
                    "context_count": result["context_count"],
                    "n_completed": result["n_completed"],
                    "keeper_digest": result["keeper_digest"],
                },
                indent=2,
            )
        )
        return 0
    keeper = validate_keeper(inventory)
    if args.dry_run:
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "platform": "Reference",
                    "md_steps": 0,
                    "minimization_steps": 0,
                    "gpu_used": False,
                    "inventory_digest": inventory["inventory_digest"],
                    "keeper_digest": keeper["keeper_digest"],
                    "cells": [row["cell_id"] for row in inventory["cells"]],
                },
                indent=2,
            )
        )
        return 0
    if args.analyze_only:
        results = []
        for cell in inventory["cells"]:
            result = _validate_completed_reference(
                inventory,
                keeper,
                cell["cell_id"],
            )
            if result is not None:
                results.append(result)
        write_summary(results, inventory, keeper)
        return 0
    return run_all(inventory, keeper)


if __name__ == "__main__":
    raise SystemExit(main())
