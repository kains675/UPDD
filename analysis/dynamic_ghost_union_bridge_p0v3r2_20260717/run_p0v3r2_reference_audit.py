#!/usr/bin/env python3
"""Build and audit frozen W4A union-solvated apex-bridge P0v3r2 cells."""

from __future__ import annotations

import argparse
import gc
import io
import json
import math
import shutil
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any

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
from analysis.dynamic_ghost_excluded_volume_20260716 import (  # noqa: E402
    run_p0_reference_audit as parent_p0,
)
from utils import atm_trackB_inplace_rbfe as rbfe  # noqa: E402
from utils import atm_trackB_setup as ats  # noqa: E402
import trackb_apex_bridge as bridge  # noqa: E402
import trackb_dynamic_ghost as dg  # noqa: E402


ANALYSIS_DIR = Path(__file__).resolve().parent
PREREGISTRATION = ANALYSIS_DIR / "PREREGISTRATION.md"
PROTOCOL_PATH = ANALYSIS_DIR / "protocol.json"
FROZEN_MANIFEST = ANALYSIS_DIR / "FROZEN_MANIFEST.json"
BASE_ANALYSIS_DIR = REPO_ROOT / "analysis/dynamic_ghost_union_bridge_p0v3_20260717"
BASE_PROTOCOL_PATH = BASE_ANALYSIS_DIR / "protocol.json"
R1_PROTOCOL_PATH = (
    REPO_ROOT / "analysis/dynamic_ghost_union_bridge_p0v3r1_20260717/protocol.json"
)

PREBUILD_INVENTORY = ANALYSIS_DIR / "P0V3R2_PREBUILD_INVENTORY.json"
SOURCE_ROOT = ANALYSIS_DIR / "u0_union_sources"
SOURCE_SUMMARY_JSON = ANALYSIS_DIR / "P0V3R2_SOURCE_SUMMARY.json"
SOURCE_SUMMARY_MD = ANALYSIS_DIR / "P0V3R2_SOURCE_SUMMARY.md"
KEEPER_PREFLIGHT = ANALYSIS_DIR / "P0V3R2_KEEPER_PREFLIGHT.json"
REFERENCE_ROOT = ANALYSIS_DIR / "p0v3r2_reference_audit"
REFERENCE_SUMMARY_JSON = ANALYSIS_DIR / "P0V3R2_REFERENCE_SUMMARY.json"
REFERENCE_SUMMARY_MD = ANALYSIS_DIR / "P0V3R2_REFERENCE_SUMMARY.md"
RUN_STATE = ANALYSIS_DIR / "P0V3R2_RUN_STATE.json"

SETUP_IMPLEMENTATION = REPO_ROOT / "utils/atm_trackB_setup.py"
SERIALIZER_IMPLEMENTATION = REPO_ROOT / "utils/atm_trackB_inplace_rbfe.py"
GHOST_IMPLEMENTATION = REPO_ROOT / "utils/trackb_dynamic_ghost.py"
BRIDGE_IMPLEMENTATION = REPO_ROOT / "utils/trackb_apex_bridge.py"
RUNNER_PATH = Path(__file__).resolve()

ENERGY_TOLERANCE = 1.0e-5
DERIVATIVE_TOLERANCE = 1.0e-5
RAW_ENERGY_MAX = 1.0e8
RAW_FORCE_MAX = 1.0e8
GEOMETRY_MIN_NM = 0.26
PLACEHOLDER_TARGET_TOLERANCE_NM = 1.0e-6
CHARGE_TOLERANCE_E = 1.0e-8
EXPECTED_ATM_PYTHON = Path("/home/san/miniconda3/envs/atm/bin/python")


relative = p0v2.relative
sha256_file = p0v2.sha256_file
sha256_text = p0v2.sha256_text
canonical_digest = p0v2.canonical_digest
artifact_record = p0v2.artifact_record
atomic_write_json = p0v2.atomic_write_json
atomic_write_text = p0v2.atomic_write_text


@dataclass(frozen=True)
class UnionCell:
    cell_id: str
    seed: str
    leg: str
    replicate_index: int
    config: dict[str, int]
    parent_cell: Any

    @property
    def source_dir(self) -> Path:
        return SOURCE_ROOT / self.cell_id

    @property
    def xml_path(self) -> Path:
        return self.source_dir / f"inplace_rbfe_{self.leg}_sys.xml"

    @property
    def pdb_path(self) -> Path:
        return self.source_dir / f"inplace_rbfe_{self.leg}.pdb"

    @property
    def build_result_path(self) -> Path:
        return self.source_dir / "source_build_result.json"


def inventory_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("inventory_digest", None)
    return canonical_digest(canonical)


def source_summary_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("source_summary_digest", None)
    return canonical_digest(canonical)


def keeper_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("keeper_digest", None)
    return canonical_digest(canonical)


def _hash_map(base: Path, mapping: dict[str, str], label: str) -> None:
    for name, expected in mapping.items():
        path = base / name
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(f"{label} drift: {path} {observed} != {expected}")


def load_protocol_and_verify_freeze(
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    frozen = json.loads(FROZEN_MANIFEST.read_text(encoding="utf-8"))
    if frozen.get("status") != "FROZEN_BEFORE_P0V3R2_IMPLEMENTATION_OR_OUTPUT":
        raise ValueError("P0v3r2 frozen manifest status drifted")
    _hash_map(ANALYSIS_DIR, frozen.get("files") or {}, "P0v3r2 frozen file")
    _hash_map(REPO_ROOT, frozen.get("inherited_design") or {}, "inherited design")
    _hash_map(REPO_ROOT, frozen.get("trigger_evidence") or {}, "trigger evidence")

    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    base_protocol = json.loads(BASE_PROTOCOL_PATH.read_text(encoding="utf-8"))
    r1_protocol = json.loads(R1_PROTOCOL_PATH.read_text(encoding="utf-8"))
    inherited = protocol.get("inherits") or {}
    if inherited.get("protocol_sha256") != sha256_file(R1_PROTOCOL_PATH):
        raise ValueError("P0v3r2 inherited P0v3r1 protocol drifted")
    if protocol.get("scival_verdict") != "CONDITIONAL_APPROVE_P0V3R2_ONLY":
        raise ValueError("P0v3r2 SciVal verdict drifted")
    if (r1_protocol.get("override") or {}).get("static_total_box_density_gate") != "REPORT_ONLY":
        raise ValueError("P0v3r2 density override drifted")
    if (r1_protocol.get("override") or {}).get("equilibrium_density_claim_allowed") is not False:
        raise ValueError("P0v3r2 equilibrium-density prohibition drifted")
    appearing_h = protocol.get("appearing_h") or {}
    if appearing_h.get("mutation_spec_object") != "w4a_trp_ala_res4":
        raise ValueError("P0v3r2 appearing-H mutation declaration drifted")
    if appearing_h.get("retry_k") != 5:
        raise ValueError("P0v3r2 appearing-H retry limit drifted")
    seed_rows = appearing_h.get("cells") or {}
    for cell in make_cells(base_protocol):
        unit_key = f"{cell.seed}|{cell.leg}|w4a_trp_ala_res4"
        expected = [ats._derive_appearing_h_seed(unit_key, i) for i in range(5)]
        if seed_rows.get(cell.cell_id) != expected:
            raise ValueError(f"{cell.cell_id}: frozen appearing-H seed trail drifted")
    rng_scope = protocol.get("rng_scope") or {}
    if rng_scope != {
        "union_addsolvent_seed_recorded": True,
        "union_addsolvent_seed_consumed": True,
        "union_addsolvent_global_state_restored": True,
        "whole_builder_global_state_invariant_gate": False,
        "cell_subprocess_isolation": True,
    }:
        raise ValueError("P0v3r2 RNG scope drifted")

    bridge_protocol, ghost_protocol = p0v2.load_protocol_and_verify_freeze()
    bridge.validate_bridge_protocol(bridge_protocol)
    dg.validate_ghost_protocol(ghost_protocol)
    return protocol, base_protocol, bridge_protocol, ghost_protocol


def make_cells(base_protocol: dict[str, Any] | None = None) -> tuple[UnionCell, ...]:
    if base_protocol is None:
        base_protocol = json.loads(BASE_PROTOCOL_PATH.read_text(encoding="utf-8"))
    configs = ((base_protocol.get("union_solvation") or {}).get("cells") or {})
    parent_cells = {
        (cell.seed, cell.leg): cell for cell in parent_p0.make_cells()
    }
    rows: list[UnionCell] = []
    for cell_id, raw in configs.items():
        seed = str(raw["seed"])
        leg = str(raw["leg"])
        expected_id = f"w4a_union_{seed}_{leg}"
        if cell_id != expected_id:
            raise ValueError(f"union cell identity drift: {cell_id} != {expected_id}")
        parent_cell = parent_cells.get((seed, leg))
        if parent_cell is None:
            raise ValueError(f"{cell_id}: no frozen uncarved parent cell")
        if int(raw["replicate_index"]) != int(parent_cell.replicate_index):
            raise ValueError(f"{cell_id}: parent replicate index drift")
        config = {
            name: int(raw[name])
            for name in ("num_added", "water", "na", "cl", "rng_seed", "parent_particles")
        }
        rows.append(
            UnionCell(
                cell_id=cell_id,
                seed=seed,
                leg=leg,
                replicate_index=int(raw["replicate_index"]),
                config=config,
                parent_cell=parent_cell,
            )
        )
    if len(rows) != 6:
        raise ValueError(f"expected six frozen union cells, found {len(rows)}")
    return tuple(rows)


def cell_by_id(cell_id: str, base_protocol: dict[str, Any] | None = None) -> UnionCell:
    matches = [cell for cell in make_cells(base_protocol) if cell.cell_id == cell_id]
    if len(matches) != 1:
        raise ValueError(f"expected one union cell {cell_id!r}, found {len(matches)}")
    return matches[0]


def _implementation_records() -> dict[str, dict[str, Any]]:
    return {
        "setup_module": artifact_record(SETUP_IMPLEMENTATION),
        "serializer_module": artifact_record(SERIALIZER_IMPLEMENTATION),
        "ghost_module": artifact_record(GHOST_IMPLEMENTATION),
        "bridge_module": artifact_record(BRIDGE_IMPLEMENTATION),
        "p0v3r2_runner": artifact_record(RUNNER_PATH),
    }


def _external_artifact_record(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return {
        "path": str(path.resolve()),
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _forcefield_records() -> dict[str, dict[str, Any]]:
    import openmm.app.forcefield as forcefield_module
    import xml.etree.ElementTree as element_tree

    data_roots = [Path(root) for root in forcefield_module._getDataDirectories()]
    records: dict[str, dict[str, Any]] = {}
    pending = list(ats.FF_FILES)
    while pending:
        name = pending.pop(0)
        if name in records:
            continue
        matches = [
            Path(root) / name
            for root in data_roots
            if (Path(root) / name).is_file()
        ]
        if len(matches) != 1:
            raise ValueError(f"force-field input {name!r} resolved to {matches}")
        path = matches[0]
        records[name] = _external_artifact_record(path)
        root = element_tree.parse(path).getroot()
        pending.extend(
            include.attrib["file"]
            for include in root.findall(".//Include")
            if include.attrib.get("file")
        )
    return dict(sorted(records.items()))


def _runtime_module_records() -> dict[str, dict[str, Any]]:
    import openmm.app.forcefield as forcefield_module
    import openmm.app.modeller as modeller_module
    import pdbfixer.pdbfixer as pdbfixer_module

    return {
        "openmm_forcefield": _external_artifact_record(
            Path(forcefield_module.__file__)
        ),
        "openmm_modeller": _external_artifact_record(Path(modeller_module.__file__)),
        "pdbfixer": _external_artifact_record(Path(pdbfixer_module.__file__)),
    }


def _positions_nm(positions: Any) -> np.ndarray:
    from openmm import unit

    if unit.is_quantity(positions):
        positions = positions.value_in_unit(unit.nanometer)
    rows = []
    for value in positions:
        if unit.is_quantity(value):
            value = value.value_in_unit(unit.nanometer)
        rows.append([float(value[0]), float(value[1]), float(value[2])])
    return np.asarray(rows, dtype=float)


def _box_matrix_nm(topology: Any) -> np.ndarray:
    from openmm import unit

    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        raise ValueError("periodic box vectors are absent")
    rows = []
    for vector in vectors:
        if unit.is_quantity(vector):
            vector = vector.value_in_unit(unit.nanometer)
        rows.append([float(vector[0]), float(vector[1]), float(vector[2])])
    matrix = np.asarray(rows, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError(f"periodic box shape is {matrix.shape}, expected (3, 3)")
    return matrix


def _box_contract(topology: Any) -> dict[str, Any]:
    matrix = _box_matrix_nm(topology)
    lengths = np.linalg.norm(matrix, axis=1)
    volume = float(np.linalg.det(matrix))
    finite = bool(np.isfinite(matrix).all())
    right_handed = bool(math.isfinite(volume) and volume > 0.0)
    orthogonal = bool(
        np.allclose(matrix @ matrix.T, np.diag(lengths ** 2), rtol=0.0, atol=1.0e-7)
    )
    cubic = bool(np.max(np.abs(lengths - lengths[0])) <= 1.0e-6)
    return {
        "vectors_nm": matrix.tolist(),
        "lengths_nm": lengths.tolist(),
        "volume_nm3": volume,
        "finite": finite,
        "right_handed": right_handed,
        "orthogonal": orthogonal,
        "cubic": cubic,
        "passed": bool(finite and right_handed and orthogonal and cubic),
    }


def _net_charge_e(system: Any) -> float:
    from openmm import unit

    _atm_index, atm_force = dg.find_atm_force(system)
    _nested_index, nonbonded = dg.resolve_nested_nonbonded(atm_force)
    return float(
        sum(
            nonbonded.getParticleParameters(index)[0].value_in_unit(
                unit.elementary_charge
            )
            for index in range(nonbonded.getNumParticles())
        )
    )


def _solute_signature(topology: Any) -> dict[str, Any]:
    rows = []
    for residue in topology.residues():
        if residue.name in ats._SOLVENT_RESNAMES:
            continue
        rows.append(
            {
                "chain": residue.chain.id,
                "residue_id": str(residue.id),
                "residue_name": residue.name,
                "atoms": [
                    {
                        "name": atom.name,
                        "element": atom.element.symbol if atom.element is not None else None,
                    }
                    for atom in residue.atoms()
                ],
            }
        )
    return {
        "n_residues": len(rows),
        "n_atoms": sum(len(row["atoms"]) for row in rows),
        "digest": canonical_digest(rows),
    }


def _parent_source_contract(cell: UnionCell) -> dict[str, Any]:
    import openmm as mm
    from openmm.app import PDBFile

    system = mm.XmlSerializer.deserialize(cell.parent_cell.xml_path.read_text(encoding="utf-8"))
    pdb = PDBFile(str(cell.parent_cell.pdb_path))
    counts = ats._solvent_counts(pdb.topology)
    contract = {
        "particles": int(system.getNumParticles()),
        "topology_atoms": int(pdb.topology.getNumAtoms()),
        "solvent_counts": counts,
        "net_charge_e": _net_charge_e(system),
        "box": _box_contract(pdb.topology),
        "solute_signature": _solute_signature(pdb.topology),
    }
    del system, pdb
    gc.collect()
    expected_counts = {
        "water": cell.config["water"],
        "na": cell.config["na"],
        "cl": cell.config["cl"],
        "total": cell.config["num_added"],
    }
    if contract["particles"] != cell.config["parent_particles"]:
        raise ValueError(f"{cell.cell_id}: frozen parent particle count drift")
    if contract["topology_atoms"] != cell.config["parent_particles"]:
        raise ValueError(f"{cell.cell_id}: frozen parent topology/System mismatch")
    if counts != expected_counts:
        raise ValueError(f"{cell.cell_id}: frozen parent solvent counts drift")
    if not contract["box"]["passed"]:
        raise ValueError(f"{cell.cell_id}: frozen parent box contract failed")
    return contract


def _scaffold_record(seed: str) -> dict[str, Any]:
    inputs = ats.resolve_leg_inputs(seed)
    path = Path(inputs["final"]["wt"])
    return artifact_record(path)


def _parent_record(
    cell: UnionCell, ghost_protocol: dict[str, Any]
) -> dict[str, Any]:
    source = parent_p0.validate_source_cell(cell.parent_cell, ghost_protocol)
    payload = {
        "cell_id": cell.cell_id,
        "seed": cell.seed,
        "leg": cell.leg,
        "replicate_index": cell.replicate_index,
        "config": cell.config,
        "scaffold": _scaffold_record(cell.seed),
        "parent_source": source,
        "parent_contract": _parent_source_contract(cell),
        "states": source["states"],
    }
    payload["cell_digest"] = canonical_digest(payload)
    return payload


def build_prebuild_inventory() -> dict[str, Any]:
    import openmm as mm

    protocol, base_protocol, bridge_protocol, ghost_protocol = (
        load_protocol_and_verify_freeze()
    )
    frozen = json.loads(FROZEN_MANIFEST.read_text(encoding="utf-8"))
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p0v3r2_prebuild_inventory_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "FROZEN_BEFORE_OFFICIAL_UNION_SOURCE_OUTPUT",
        "platform": "Reference",
        "python": sys.executable,
        "openmm_version": mm.__version__,
        "md_allowed": False,
        "minimization_allowed": False,
        "gpu_allowed": False,
        "free_energy_estimation_allowed": False,
        "equilibrium_density_claim_allowed": False,
        "bridge_states": (base_protocol.get("apex_bridge") or {}).get("audit_states"),
        "ghost_g": (base_protocol.get("dynamic_ghost") or {}).get("bridge_g"),
        "frozen_manifest": artifact_record(FROZEN_MANIFEST),
        "frozen_files": {
            name: artifact_record(ANALYSIS_DIR / name)
            for name in sorted((frozen.get("files") or {}))
        },
        "implementation": _implementation_records(),
        "forcefields": _forcefield_records(),
        "runtime_modules": _runtime_module_records(),
        "protocol_digests": {
            "revision": canonical_digest(protocol),
            "base": canonical_digest(base_protocol),
            "bridge": canonical_digest(bridge_protocol),
            "ghost": canonical_digest(ghost_protocol),
        },
        "cells": [
            _parent_record(cell, ghost_protocol) for cell in make_cells(base_protocol)
        ],
    }
    payload["inventory_digest"] = inventory_digest(payload)
    return payload


def load_prebuild_inventory(expected_digest: str | None = None) -> dict[str, Any]:
    if not PREBUILD_INVENTORY.is_file():
        raise FileNotFoundError(
            f"missing {PREBUILD_INVENTORY}; run --prepare-inventory first"
        )
    payload = json.loads(PREBUILD_INVENTORY.read_text(encoding="utf-8"))
    if payload.get("inventory_digest") != inventory_digest(payload):
        raise ValueError("stored P0v3r2 prebuild inventory digest is invalid")
    if expected_digest is not None and payload["inventory_digest"] != expected_digest:
        raise ValueError("worker prebuild inventory digest mismatch")
    return payload


def verify_prebuild_inventory() -> dict[str, Any]:
    frozen = load_prebuild_inventory()
    current = build_prebuild_inventory()
    if frozen["inventory_digest"] != current["inventory_digest"]:
        raise ValueError(
            "P0v3r2 source/code drift after prebuild freeze: "
            f"{frozen['inventory_digest']} != {current['inventory_digest']}"
        )
    return frozen


def frozen_cell_record(inventory: dict[str, Any], cell: UnionCell) -> dict[str, Any]:
    rows = [row for row in inventory.get("cells", []) if row.get("cell_id") == cell.cell_id]
    if len(rows) != 1:
        raise ValueError(f"{cell.cell_id}: missing or duplicate inventory row")
    return rows[0]


def verify_worker_inputs(
    inventory: dict[str, Any], cell: UnionCell
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    protocol, base_protocol, bridge_protocol, ghost_protocol = (
        load_protocol_and_verify_freeze()
    )
    if _implementation_records() != inventory.get("implementation"):
        raise ValueError(f"{cell.cell_id}: implementation drift after inventory freeze")
    if _forcefield_records() != inventory.get("forcefields"):
        raise ValueError(f"{cell.cell_id}: force-field drift after inventory freeze")
    if _runtime_module_records() != inventory.get("runtime_modules"):
        raise ValueError(f"{cell.cell_id}: runtime module drift after inventory freeze")
    frozen_cell = frozen_cell_record(inventory, cell)
    current_parent = _parent_record(cell, ghost_protocol)
    if current_parent != frozen_cell:
        raise ValueError(f"{cell.cell_id}: parent/scaffold input drift")
    return protocol, base_protocol, bridge_protocol, ghost_protocol, frozen_cell


def _minimum_image_deltas(deltas: np.ndarray, box: np.ndarray) -> np.ndarray:
    inverse = np.linalg.inv(box)
    fractional = deltas @ inverse
    fractional -= np.round(fractional)
    return fractional @ box


def _transformed_ring_positions(
    system: Any, topology: Any, positions_nm: np.ndarray
) -> tuple[Any, np.ndarray, np.ndarray]:
    import openmm as mm

    selection = dg.select_ghost_atoms(topology)
    _atm_index, atm_force = dg.find_atm_force(system)
    stored = positions_nm[np.asarray(selection.ring_indices, dtype=int)]
    transformed_rows = []
    for index in selection.ring_indices:
        transformation = atm_force.getParticleTransformation(int(index))
        if not isinstance(transformation, mm.ParticleOffsetDisplacement):
            raise ValueError(
                f"ring atom {index} has {type(transformation).__name__}, "
                "expected ParticleOffsetDisplacement"
            )
        destination = int(transformation.getDestinationParticle1())
        origin = int(transformation.getOriginParticle1())
        transformed_rows.append(
            positions_nm[int(index)] + positions_nm[destination] - positions_nm[origin]
        )
    return selection, stored, np.asarray(transformed_rows, dtype=float)


def _ring_water_geometry(
    system: Any, topology: Any, positions: Any
) -> dict[str, Any]:
    positions_nm = _positions_nm(positions)
    box = _box_matrix_nm(topology)
    selection, stored, transformed = _transformed_ring_positions(
        system, topology, positions_nm
    )
    water = positions_nm[np.asarray(selection.water_oxygen_indices, dtype=int)]
    rows = {}
    for label, ring in (("stored", stored), ("transformed", transformed)):
        deltas = ring[:, None, :] - water[None, :, :]
        deltas = _minimum_image_deltas(deltas, box)
        distances = np.linalg.norm(deltas, axis=-1)
        minimum = float(np.min(distances))
        count_below = int(np.count_nonzero(distances < GEOMETRY_MIN_NM))
        rows[label] = {
            "min_ring_water_oxygen_distance_nm": minimum,
            "pairs_below_0_26_nm": count_below,
            "passed": bool(minimum >= GEOMETRY_MIN_NM and count_below == 0),
        }
    return {
        "n_ring_atoms": len(selection.ring_indices),
        "n_water_oxygens": len(selection.water_oxygen_indices),
        "states": rows,
        "passed": all(row["passed"] for row in rows.values()),
    }


def _placeholder_parameter_audit(
    system: Any,
    topology: Any,
    positions: Any,
    report: dict[str, Any],
    *,
    target_tolerance_nm: float = PLACEHOLDER_TARGET_TOLERANCE_NM,
) -> dict[str, Any]:
    import openmm as mm
    from openmm import app, unit

    declarations = report.get("placeholder_declarations") or []
    xml = ats._union_placeholder_forcefield_xml(declarations)
    xml_hash_matches = (
        sha256_text(xml) == report.get("placeholder_forcefield_xml_sha256")
    )
    forcefield = app.ForceField(io.StringIO(xml))
    dummy_topology = app.Topology()
    chain = dummy_topology.addChain("U")
    residue_templates = {}
    for row in declarations:
        residue = dummy_topology.addResidue(row["placeholder_residue"], chain)
        dummy_topology.addAtom(
            "DU", app.Element.getBySymbol(row["element"]), residue
        )
        residue_templates[residue] = row["placeholder_residue"]
    dummy_system = forcefield.createSystem(
        dummy_topology,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        residueTemplates=residue_templates,
        removeCMMotion=False,
    )
    dummy_nb = ats._nonbonded_force(dummy_system)
    _atm_index, atm_force = dg.find_atm_force(system)
    _nested_index, source_nb = dg.resolve_nested_nonbonded(atm_force)
    atoms = {atom.index: atom for atom in topology.atoms()}
    positions_nm = _positions_nm(positions)
    _selection, _stored, transformed = _transformed_ring_positions(
        system, topology, positions_nm
    )
    ring_indices = list(dg.select_ghost_atoms(topology).ring_indices)
    transformed_by_index = {
        index: transformed[position] for position, index in enumerate(ring_indices)
    }
    rows = []
    for placeholder_index, declaration in enumerate(declarations):
        source_index = int(declaration["source_atom_index"])
        source_charge, source_sigma, source_epsilon = source_nb.getParticleParameters(
            source_index
        )
        dummy_charge, dummy_sigma, dummy_epsilon = dummy_nb.getParticleParameters(
            placeholder_index
        )
        source_mass = system.getParticleMass(source_index).value_in_unit(unit.dalton)
        dummy_mass = dummy_system.getParticleMass(placeholder_index).value_in_unit(
            unit.dalton
        )
        target = np.asarray(declaration["solvated_target_position_nm"], dtype=float)
        transformed_target = transformed_by_index.get(source_index)
        target_delta = (
            math.inf
            if transformed_target is None
            else float(np.max(np.abs(target - transformed_target)))
        )
        atom = atoms[source_index]
        checks = {
            "source_name": atom.name == declaration["source_atom_name"],
            "source_element": (
                atom.element is not None
                and atom.element.symbol == declaration["element"]
            ),
            "source_mass": float(source_mass) == float(declaration["mass_da"]),
            "source_charge": float(
                source_charge.value_in_unit(unit.elementary_charge)
            ) == float(declaration["source_charge_e"]),
            "source_sigma": float(
                source_sigma.value_in_unit(unit.nanometer)
            ) == float(declaration["sigma_nm"]),
            "source_epsilon": float(
                source_epsilon.value_in_unit(unit.kilojoule_per_mole)
            ) == float(declaration["epsilon_kj_mol"]),
            "placeholder_mass_exact": float(dummy_mass) == float(source_mass),
            "placeholder_charge_zero": float(
                dummy_charge.value_in_unit(unit.elementary_charge)
            ) == 0.0,
            "placeholder_sigma_exact": float(
                dummy_sigma.value_in_unit(unit.nanometer)
            ) == float(source_sigma.value_in_unit(unit.nanometer)),
            "placeholder_epsilon_exact": float(
                dummy_epsilon.value_in_unit(unit.kilojoule_per_mole)
            ) == float(source_epsilon.value_in_unit(unit.kilojoule_per_mole)),
            "target_matches_atm_transform": target_delta <= target_tolerance_nm,
        }
        rows.append(
            {
                "placeholder_index": placeholder_index,
                "source_atom_index": source_index,
                "source_atom_name": atom.name,
                "target_max_abs_delta_nm": target_delta,
                "checks": checks,
                "passed": all(checks.values()),
            }
        )
    structural_checks = {
        "xml_hash_matches": xml_hash_matches,
        "exact_nine_placeholders": len(declarations) == 9,
        "one_force_only": dummy_system.getNumForces() == 1,
        "nonbonded_only": dummy_system.getNumForces() == 1,
        "zero_constraints": dummy_system.getNumConstraints() == 0,
        "zero_virtual_sites": all(
            not dummy_system.isVirtualSite(index)
            for index in range(dummy_system.getNumParticles())
        ),
    }
    return {
        "placeholder_count": len(declarations),
        "structural_checks": structural_checks,
        "atoms": rows,
        "max_target_delta_nm": max(
            (row["target_max_abs_delta_nm"] for row in rows), default=math.inf
        ),
        "target_tolerance_nm": target_tolerance_nm,
        "passed": bool(
            all(structural_checks.values()) and all(row["passed"] for row in rows)
        ),
    }


def _union_source_audit(
    cell: UnionCell,
    serialized: dict[str, Any],
    frozen_cell: dict[str, Any],
) -> dict[str, Any]:
    build = serialized["_build"]
    fused = build["fused_build"]
    system = fused["system"]
    modeller = fused["modeller"]
    report = serialized.get("union_solvation_report") or {}
    expected_counts = {
        "water": cell.config["water"],
        "na": cell.config["na"],
        "cl": cell.config["cl"],
        "total": cell.config["num_added"],
    }
    counts = ats._solvent_counts(modeller.topology)
    parent_contract = frozen_cell["parent_contract"]
    box = _box_contract(modeller.topology)
    geometry = _ring_water_geometry(system, modeller.topology, modeller.positions)
    placeholder = _placeholder_parameter_audit(
        system, modeller.topology, modeller.positions, report
    )
    final_residue_names = {residue.name for residue in modeller.topology.residues()}
    placeholder_names = {
        row["placeholder_residue"]
        for row in report.get("placeholder_declarations") or []
    }
    net_charge = _net_charge_e(system)
    parent_volume = float(parent_contract["box"]["volume_nm3"])
    density_parent = cell.config["num_added"] / parent_volume
    density_union = cell.config["num_added"] / float(box["volume_nm3"])
    density_diagnostic = {
        "claim_status": "REPORT_ONLY_NOT_EQUILIBRIUM_DENSITY",
        "parent_total_solvent_per_nm3": density_parent,
        "union_total_solvent_per_nm3": density_union,
        "relative_delta": density_union / density_parent - 1.0,
        "hard_gate": False,
        "future_npt_preregistration_required": True,
    }
    separation = serialized.get("separation") or {}
    checks = {
        "system_topology_lockstep": (
            system.getNumParticles() == modeller.topology.getNumAtoms()
        ),
        "parent_particle_count_exact": (
            system.getNumParticles() == cell.config["parent_particles"]
        ),
        "solvent_counts_exact": counts == expected_counts,
        "report_counts_exact": report.get("final_solvent_counts") == expected_counts,
        "placeholder_report_exact_nine": report.get("placeholder_count") == 9,
        "placeholder_report_removed": report.get("placeholders_removed") is True,
        "placeholder_names_absent": not (placeholder_names & final_residue_names),
        "placeholder_parameters_and_targets": placeholder["passed"],
        "rng_seed_recorded": (report.get("config") or {}).get("rng_seed")
        == cell.config["rng_seed"],
        "rng_consumed": report.get("seeded_rng_state_changed") is True,
        "global_rng_restored": report.get("global_rng_state_restored") is True,
        "static_carve_absent": (
            serialized.get("carve_report") is None
            and report.get("carve_void_waters") is False
        ),
        "net_charge_matches_parent": abs(
            net_charge - float(parent_contract["net_charge_e"])
        ) <= CHARGE_TOLERANCE_E,
        "solute_order_unchanged": (
            _solute_signature(modeller.topology) == parent_contract["solute_signature"]
        ),
        "box_contract": box["passed"],
        "stored_and_transformed_geometry": geometry["passed"],
        "canonical_twocopy_separation": separation.get("passed") is True,
    }
    return {
        "status": "PASS" if all(checks.values()) else "REJECT_UNION_SOURCE",
        "checks": checks,
        "counts": counts,
        "net_charge_e": net_charge,
        "parent_net_charge_e": parent_contract["net_charge_e"],
        "box": box,
        "density_diagnostic": density_diagnostic,
        "geometry": geometry,
        "placeholder": placeholder,
        "separation": separation,
    }


def write_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S %z')} {message}"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def _appearing_h_retry_audit(
    cell: UnionCell,
    serialized: dict[str, Any],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    report = (serialized.get("_build") or {}).get("appearing_h_retry") or {}
    unit_key = f"{cell.seed}|{cell.leg}|w4a_trp_ala_res4"
    frozen_seeds = list((protocol.get("appearing_h") or {})["cells"][cell.cell_id])
    attempts_used = int(report.get("attempts_used", 0))
    expected_tried = frozen_seeds[:attempts_used]
    checks = {
        "unit_key_exact": report.get("unit_key") == unit_key,
        "retry_k_exact": report.get("retry_k") == 5,
        "attempts_in_range": 1 <= attempts_used <= 5,
        "seed_trail_exact": report.get("seeds_tried") == expected_tried,
        "accepted_seed_exact": (
            bool(expected_tried) and report.get("seed_used") == expected_tried[-1]
        ),
    }
    return {
        "unit_key": unit_key,
        "frozen_attempt_seeds": frozen_seeds,
        "report": report,
        "checks": checks,
        "passed": all(checks.values()),
    }


def run_source_worker(cell: UnionCell, inventory: dict[str, Any]) -> dict[str, Any]:
    if cell.source_dir.exists():
        raise RuntimeError(f"{cell.cell_id}: source output directory already exists")
    cell.source_dir.mkdir(parents=True)
    log_path = cell.source_dir / "source_build.log"
    started = time.monotonic()
    write_log(log_path, f"START {cell.cell_id} isolated union source build")
    protocol, _base, _bridge_protocol, _ghost_protocol, frozen_cell = (
        verify_worker_inputs(inventory, cell)
    )
    try:
        union_config = {
            name: cell.config[name]
            for name in ("num_added", "water", "na", "cl", "rng_seed")
        }
        write_log(log_path, "build canonical two-copy source with temporary placeholders")
        mutation_spec = rbfe.ats.resolve_mutation_spec("w4a_trp_ala_res4")
        serialized = rbfe.serialize_inplace_rbfe_system(
            leg=cell.leg,
            out_dir=str(cell.source_dir),
            seed=cell.seed,
            solvate=True,
            constraints=None,
            tag=cell.leg,
            construction="twocopy",
            displacement_nm=4.0,
            mutation_spec=mutation_spec,
            auto_search_displacement=False,
            appearing_h_retry_k=5,
            carve_void_waters=False,
            union_solvation=union_config,
        )
        write_log(log_path, "audit exact counts, placeholder parameters, and both images")
        audit = _union_source_audit(cell, serialized, frozen_cell)
        appearing_h_audit = _appearing_h_retry_audit(cell, serialized, protocol)
        audit["checks"]["deterministic_appearing_h_retry"] = appearing_h_audit[
            "passed"
        ]
        audit["appearing_h_retry"] = appearing_h_audit
        if not appearing_h_audit["passed"]:
            audit["status"] = "REJECT_UNION_SOURCE"
        artifacts = {
            "serialized_xml": artifact_record(cell.xml_path),
            "serialized_pdb": artifact_record(cell.pdb_path),
        }
        result = {
            "schema": "updd_dynamic_ghost_union_bridge_p0v3r2_source_cell_v1",
            "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "status": audit["status"],
            "cell": {
                "cell_id": cell.cell_id,
                "seed": cell.seed,
                "leg": cell.leg,
                "replicate_index": cell.replicate_index,
            },
            "inventory_digest": inventory["inventory_digest"],
            "frozen_cell_digest": frozen_cell["cell_digest"],
            "config": cell.config,
            "artifacts": artifacts,
            "construction": {
                "construction": serialized["construction"],
                "mutation_spec": "w4a_trp_ala_res4",
                "displacement_vector_nm": serialized["displacement_vector_nm"],
                "displacement_mode": serialized["displacement_mode"],
                "n_copy1": serialized["n_copy1"],
                "union_solvation_report": serialized["union_solvation_report"],
                "appearing_h_retry": serialized["_build"]["appearing_h_retry"],
                "carve_report": serialized["carve_report"],
                "swap": serialized["swap"],
            },
            "audit": audit,
            "platform": {
                "contexts_created": 0,
                "md_steps": 0,
                "minimization_steps": 0,
                "gpu_used": False,
                "free_energy_estimated": False,
            },
            "elapsed_s": time.monotonic() - started,
        }
    except Exception as exc:  # source construction failure is a frozen reject outcome
        write_log(log_path, f"REJECT {type(exc).__name__}: {exc}")
        result = {
            "schema": "updd_dynamic_ghost_union_bridge_p0v3r2_source_cell_v1",
            "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "status": "REJECT_UNION_SOURCE",
            "cell": {
                "cell_id": cell.cell_id,
                "seed": cell.seed,
                "leg": cell.leg,
                "replicate_index": cell.replicate_index,
            },
            "inventory_digest": inventory["inventory_digest"],
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
            "platform": {
                "contexts_created": 0,
                "md_steps": 0,
                "minimization_steps": 0,
                "gpu_used": False,
                "free_energy_estimated": False,
            },
            "elapsed_s": time.monotonic() - started,
        }
    atomic_write_json(cell.build_result_path, result)
    write_log(log_path, f"COMPLETE {cell.cell_id} status={result['status']}")
    return result


def _archive_partial(root: Path, archive_root: Path, reason: str) -> Path | None:
    if not root.exists():
        return None
    stamp = time.strftime("%Y%m%d_%H%M%S")
    target = archive_root / f"{stamp}_{root.name}"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(root), str(target))
    atomic_write_json(
        target / "archive_reason.json",
        {
            "archived": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "reason": reason,
        },
    )
    return target


def validate_source_result(
    cell: UnionCell, inventory: dict[str, Any]
) -> dict[str, Any] | None:
    if not cell.build_result_path.is_file():
        return None
    try:
        result = json.loads(cell.build_result_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    expected_cell = {
        "cell_id": cell.cell_id,
        "seed": cell.seed,
        "leg": cell.leg,
        "replicate_index": cell.replicate_index,
    }
    if result.get("cell") != expected_cell:
        raise ValueError(f"{cell.cell_id}: source result identity mismatch")
    if result.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError(f"{cell.cell_id}: source result inventory mismatch")
    if result.get("status") not in {"PASS", "REJECT_UNION_SOURCE"}:
        raise ValueError(f"{cell.cell_id}: invalid source result status")
    if result["status"] == "PASS":
        if result.get("frozen_cell_digest") != frozen_cell_record(
            inventory, cell
        )["cell_digest"]:
            raise ValueError(f"{cell.cell_id}: frozen cell digest mismatch")
        for record in (result.get("artifacts") or {}).values():
            path = REPO_ROOT / record["path"]
            if not path.is_file() or sha256_file(path) != record["sha256"]:
                raise ValueError(f"{cell.cell_id}: source artifact missing or drifted")
        if not all((result.get("audit") or {}).get("checks", {}).values()):
            raise ValueError(f"{cell.cell_id}: PASS has a false source check")
    return result


def _write_source_summary(
    results: list[dict[str, Any]], inventory: dict[str, Any]
) -> dict[str, Any]:
    cells = make_cells()
    by_id = {row["cell"]["cell_id"]: row for row in results}
    pending = [cell.cell_id for cell in cells if cell.cell_id not in by_id]
    rejected = [
        cell_id
        for cell_id, row in by_id.items()
        if row.get("status") == "REJECT_UNION_SOURCE"
    ]
    if rejected:
        status = "REJECT_UNION_SOURCE"
    elif not pending and len(results) == len(cells):
        status = "PASS"
    else:
        status = "RUNNING"
    elapsed = [float(row.get("elapsed_s", 0.0)) for row in results]
    median_elapsed = median(elapsed) if elapsed else None
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p0v3r2_source_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "inventory_digest": inventory["inventory_digest"],
        "context_count": 0,
        "md_steps": 0,
        "minimization_steps": 0,
        "gpu_used": False,
        "free_energy_estimated": False,
        "n_expected": len(cells),
        "n_completed": len(results),
        "pending_cell_ids": pending,
        "rejected_cell_ids": rejected,
        "median_elapsed_s": median_elapsed,
        "eta_s": median_elapsed * len(pending) if median_elapsed is not None else None,
        "cells": [
            {
                "cell": row["cell"],
                "status": row["status"],
                "elapsed_s": row.get("elapsed_s"),
                "result": artifact_record(
                    cell_by_id(row["cell"]["cell_id"]).build_result_path
                ),
                "artifacts": row.get("artifacts"),
            }
            for row in results
        ],
    }
    payload["source_summary_digest"] = source_summary_digest(payload)
    atomic_write_json(SOURCE_SUMMARY_JSON, payload)
    lines = [
        "# Union Source P0v3r2 Summary",
        "",
        f"Generated: {payload['generated']}",
        "",
        f"Status: **{status}**",
        "",
        "All builds use isolated subprocesses, zero OpenMM Contexts, zero MD or "
        "minimization steps, and no GPU or free-energy estimator.",
        "",
        "| cell | status | elapsed (s) |",
        "|---|---|---:|",
    ]
    for row in payload["cells"]:
        lines.append(
            f"| {row['cell']['cell_id']} | {row['status']} | "
            f"{float(row['elapsed_s'] or 0.0):.3f} |"
        )
    if pending:
        lines.extend(["", "Pending: " + ", ".join(f"`{item}`" for item in pending)])
    if rejected:
        lines.extend(["", "Rejected: " + ", ".join(f"`{item}`" for item in rejected)])
    atomic_write_text(SOURCE_SUMMARY_MD, "\n".join(lines) + "\n")
    return payload


def run_source_parent(inventory: dict[str, Any]) -> dict[str, Any]:
    SOURCE_ROOT.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for cell in make_cells():
        completed = validate_source_result(cell, inventory)
        if completed is not None:
            results.append(completed)
            _write_source_summary(results, inventory)
            if completed["status"] == "REJECT_UNION_SOURCE":
                return _write_source_summary(results, inventory)
            continue
        if cell.source_dir.exists():
            _archive_partial(
                cell.source_dir,
                SOURCE_ROOT / "_archive",
                "incomplete or invalid P0v3r2 source build",
            )
        command = [
            sys.executable,
            str(RUNNER_PATH),
            "--worker-build-cell",
            cell.cell_id,
            "--expected-inventory-digest",
            inventory["inventory_digest"],
        ]
        returncode = subprocess.run(command, cwd=REPO_ROOT, check=False).returncode
        if returncode != 0:
            _write_source_summary(results, inventory)
            raise RuntimeError(
                f"source worker {cell.cell_id} exited with status {returncode}"
            )
        completed = validate_source_result(cell, inventory)
        if completed is None:
            raise RuntimeError(f"{cell.cell_id}: source worker produced no result")
        results.append(completed)
        _write_source_summary(results, inventory)
        if completed["status"] == "REJECT_UNION_SOURCE":
            return _write_source_summary(results, inventory)
    return _write_source_summary(results, inventory)


def load_source_summary_shallow(inventory: dict[str, Any]) -> dict[str, Any]:
    if not SOURCE_SUMMARY_JSON.is_file():
        raise FileNotFoundError(
            f"missing {SOURCE_SUMMARY_JSON}; run --build-sources first"
        )
    payload = json.loads(SOURCE_SUMMARY_JSON.read_text(encoding="utf-8"))
    if payload.get("source_summary_digest") != source_summary_digest(payload):
        raise ValueError("stored P0v3r2 source summary digest is invalid")
    if payload.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError("P0v3r2 source summary inventory drifted")
    if payload.get("status") != "PASS" or payload.get("n_completed") != 6:
        raise ValueError("P0v3r2 source summary is not a complete PASS")
    return payload


def load_source_summary(inventory: dict[str, Any]) -> dict[str, Any]:
    payload = load_source_summary_shallow(inventory)
    summary_rows = {
        row["cell"]["cell_id"]: row for row in payload.get("cells", [])
    }
    for cell in make_cells():
        result = validate_source_result(cell, inventory)
        if result is None or result.get("status") != "PASS":
            raise ValueError(f"{cell.cell_id}: source result is not a valid PASS")
        summary_row = summary_rows.get(cell.cell_id)
        if summary_row is None:
            raise ValueError(f"{cell.cell_id}: source summary row is absent")
        if summary_row.get("status") != result["status"]:
            raise ValueError(f"{cell.cell_id}: source summary status drifted")
        if summary_row.get("result") != artifact_record(cell.build_result_path):
            raise ValueError(f"{cell.cell_id}: source summary result hash drifted")
        if summary_row.get("artifacts") != result.get("artifacts"):
            raise ValueError(f"{cell.cell_id}: source summary artifacts drifted")
    return payload


def run_keeper_preflight(inventory: dict[str, Any]) -> dict[str, Any]:
    """Re-audit all final sources and mutation contracts without a Context."""
    import openmm as mm
    from openmm.app import PDBFile

    source_summary = load_source_summary(inventory)
    _protocol, _base, bridge_protocol, ghost_protocol = (
        load_protocol_and_verify_freeze()
    )
    rows: list[dict[str, Any]] = []
    for cell in make_cells():
        result = validate_source_result(cell, inventory)
        if result is None or result.get("status") != "PASS":
            raise ValueError(f"{cell.cell_id}: KEEPER source is not a PASS")
        frozen_cell = frozen_cell_record(inventory, cell)
        system = mm.XmlSerializer.deserialize(cell.xml_path.read_text(encoding="utf-8"))
        pdb = PDBFile(str(cell.pdb_path))
        expected_counts = {
            "water": cell.config["water"],
            "na": cell.config["na"],
            "cl": cell.config["cl"],
            "total": cell.config["num_added"],
        }
        counts = ats._solvent_counts(pdb.topology)
        parent_contract = frozen_cell["parent_contract"]
        box = _box_contract(pdb.topology)
        geometry = _ring_water_geometry(system, pdb.topology, pdb.positions)
        union_report = result["construction"]["union_solvation_report"]
        # The exact 1e-6 nm target gate was evaluated in memory before PDB
        # serialization. PDB coordinates are independently checked here at their
        # finite decimal precision while the hashed build result preserves the
        # exact pre-serialization result.
        placeholder_roundtrip = _placeholder_parameter_audit(
            system,
            pdb.topology,
            pdb.positions,
            union_report,
            target_tolerance_nm=2.0e-4,
        )
        final_residue_names = {residue.name for residue in pdb.topology.residues()}
        placeholder_names = {
            row["placeholder_residue"]
            for row in union_report.get("placeholder_declarations") or []
        }
        ghost_report = dg.build_dynamic_ghost_force(system, pdb.topology, ghost_protocol)
        source_contract = bridge.inspect_atm_contract(system)
        bridge_report = bridge.build_apex_bridge(system, bridge_protocol)
        ghost_after = dg.inspect_dynamic_ghost_force(system)
        bridge_after = bridge.inspect_apex_bridge(system)
        checks = {
            "source_result_hash_locked": all(
                (result.get("audit") or {}).get("checks", {}).values()
            ),
            "exact_in_memory_target_gate_preserved": (
                float(result["audit"]["placeholder"]["max_target_delta_nm"])
                <= PLACEHOLDER_TARGET_TOLERANCE_NM
            ),
            "serialized_placeholder_parameter_roundtrip": placeholder_roundtrip[
                "passed"
            ],
            "placeholder_names_absent": not (placeholder_names & final_residue_names),
            "system_topology_lockstep": (
                system.getNumParticles() == pdb.topology.getNumAtoms()
            ),
            "particle_count_exact": (
                system.getNumParticles() == cell.config["parent_particles"]
            ),
            "solvent_counts_exact": counts == expected_counts,
            "charge_matches_parent": abs(
                _net_charge_e(system) - float(parent_contract["net_charge_e"])
            ) <= CHARGE_TOLERANCE_E,
            "solute_order_unchanged": (
                _solute_signature(pdb.topology) == parent_contract["solute_signature"]
            ),
            "box_contract": box["passed"],
            "stored_and_transformed_geometry": geometry["passed"],
            "exact_nine_ring_atoms": ghost_report["selection"]["n_ring_atoms"] == 9,
            "all_water_oxygens": (
                ghost_report["selection"]["n_water_oxygens"]
                == ghost_report["selection"]["n_water_residues"]
                and ghost_report["selection"]["n_water_oxygens"] > 0
            ),
            "ghost_particle_and_charge_invariant": (
                ghost_report["n_particles_before"] == ghost_report["n_particles_after"]
                and ghost_report["net_charge_delta_e"] == 0.0
            ),
            "bridge_mutation_scope_exact": all(bridge_report["checks"].values()),
            "source_contract_preserved_by_bridge": (
                source_contract == bridge_report["before"]
            ),
            "ghost_after_atm": ghost_after["force_index"] > bridge_after["atm_index"],
            "ghost_force_group_31": ghost_after["force_group"] == dg.GHOST_FORCE_GROUP,
            "bridge_parameter_consumed": (
                bridge.BRIDGE_PARAMETER in bridge_after["energy_function"]
                and bridge_after["energy_parameter_derivatives"]
                == [bridge.BRIDGE_PARAMETER]
            ),
            "context_free": True,
        }
        if not all(checks.values()):
            raise ValueError(f"{cell.cell_id}: KEEPER checks failed: {checks}")
        rows.append(
            {
                "cell_id": cell.cell_id,
                "seed": cell.seed,
                "leg": cell.leg,
                "source_result": artifact_record(cell.build_result_path),
                "checks": checks,
                "counts": counts,
                "box": box,
                "geometry": geometry,
                "placeholder_roundtrip": placeholder_roundtrip,
                "source_atm_contract": source_contract,
                "ghost_contract": parent_p0.compact_force_report(ghost_report),
                "bridge_contract": p0v2._compact_bridge_report(bridge_report),
            }
        )
        del system, pdb, ghost_report, bridge_report
        gc.collect()

    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p0v3r2_keeper_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "PASS",
        "inventory_digest": inventory["inventory_digest"],
        "source_summary_digest": source_summary["source_summary_digest"],
        "context_count": 0,
        "md_steps": 0,
        "minimization_steps": 0,
        "gpu_used": False,
        "free_energy_estimated": False,
        "n_cells": len(rows),
        "cells": rows,
    }
    payload["keeper_digest"] = keeper_digest(payload)
    atomic_write_json(KEEPER_PREFLIGHT, payload)
    return payload


def validate_keeper_preflight(
    inventory: dict[str, Any], source_summary: dict[str, Any]
) -> dict[str, Any]:
    if not KEEPER_PREFLIGHT.is_file():
        raise FileNotFoundError(
            f"missing {KEEPER_PREFLIGHT}; run --keeper-preflight first"
        )
    payload = json.loads(KEEPER_PREFLIGHT.read_text(encoding="utf-8"))
    if payload.get("keeper_digest") != keeper_digest(payload):
        raise ValueError("P0v3r2 KEEPER digest is invalid")
    if payload.get("status") != "PASS" or payload.get("context_count") != 0:
        raise ValueError("P0v3r2 KEEPER preflight is not a context-free PASS")
    if payload.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError("P0v3r2 KEEPER inventory digest drifted")
    if payload.get("source_summary_digest") != source_summary["source_summary_digest"]:
        raise ValueError("P0v3r2 KEEPER source summary digest drifted")
    if payload.get("n_cells") != len(make_cells()):
        raise ValueError("P0v3r2 KEEPER cell count drifted")
    return payload


def reference_cell_dir(cell: UnionCell) -> Path:
    return REFERENCE_ROOT / cell.cell_id


def reference_result_path(cell: UnionCell) -> Path:
    return reference_cell_dir(cell) / "cell_result.json"


def reference_bridge_xml_path(cell: UnionCell) -> Path:
    return reference_cell_dir(cell) / "apex_bridge_system.xml"


def _evaluate_raw_endpoints(
    system: Any,
    positions: Any,
    state: dict[str, Any],
) -> dict[str, Any]:
    import openmm as mm
    from openmm import unit

    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    p0v2._set_source_parameters(context, state, 0.5)
    _atm_index, atm_force = dg.find_atm_force(system)
    raw_u1, raw_u0, _bias = atm_force.getPerturbationEnergy(context)
    payload = {
        "u0_kj_mol": float(raw_u0.value_in_unit(unit.kilojoule_per_mole)),
        "u1_kj_mol": float(raw_u1.value_in_unit(unit.kilojoule_per_mole)),
    }
    payload["delta_u1_minus_u0_kj_mol"] = payload["u1_kj_mol"] - payload["u0_kj_mol"]
    payload["finite"] = all(math.isfinite(value) for value in payload.values())
    del context, integrator
    gc.collect()
    return payload


def run_reference_worker(
    cell: UnionCell, inventory: dict[str, Any], source_summary: dict[str, Any]
) -> dict[str, Any]:
    import openmm as mm
    from openmm.app import PDBFile

    output_dir = reference_cell_dir(cell)
    if output_dir.exists():
        raise RuntimeError(f"{cell.cell_id}: Reference output directory already exists")
    output_dir.mkdir(parents=True)
    log_path = output_dir / "worker.log"
    started = time.monotonic()
    write_log(log_path, f"START {cell.cell_id} platform=Reference md_steps=0")

    _protocol, _base, bridge_protocol, ghost_protocol, frozen_cell = (
        verify_worker_inputs(inventory, cell)
    )
    validated_source_summary = load_source_summary_shallow(inventory)
    if validated_source_summary["source_summary_digest"] != source_summary[
        "source_summary_digest"
    ]:
        raise ValueError(f"{cell.cell_id}: source summary drift in worker")
    source_result = validate_source_result(cell, inventory)
    if source_result is None or source_result.get("status") != "PASS":
        raise ValueError(f"{cell.cell_id}: Reference source is not a frozen PASS")

    source_xml = cell.xml_path.read_text(encoding="utf-8")
    pdb = PDBFile(str(cell.pdb_path))
    states = frozen_cell["states"]
    xi_values = [float(value) for value in inventory["bridge_states"]]

    write_log(log_path, "evaluate dynamic-ghost source apexes and raw u0/u1")
    source_system = mm.XmlSerializer.deserialize(source_xml)
    source_ghost_report = dg.build_dynamic_ghost_force(
        source_system, pdb.topology, ghost_protocol
    )
    source_values = p0v2.evaluate_source_apexes(source_system, pdb.positions, states)
    raw_endpoints = _evaluate_raw_endpoints(
        source_system, pdb.positions, states["apex_dplus"]
    )
    max_source_force = max(
        float(np.max(np.abs(value["forces"]))) for value in source_values.values()
    )
    raw_guard = {
        "raw_u0_abs_below_1e8": (
            raw_endpoints["finite"] and abs(raw_endpoints["u0_kj_mol"]) < RAW_ENERGY_MAX
        ),
        "raw_u1_abs_below_1e8": (
            raw_endpoints["finite"] and abs(raw_endpoints["u1_kj_mol"]) < RAW_ENERGY_MAX
        ),
        "max_force_component_below_1e8": (
            math.isfinite(max_source_force) and max_source_force < RAW_FORCE_MAX
        ),
    }
    raw_guard["passed"] = all(raw_guard.values())
    if not raw_guard["passed"]:
        result = {
            "schema": "updd_dynamic_ghost_union_bridge_p0v3r2_cell_v1",
            "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "status": "REJECT_UNION_SOURCE",
            "cell": {
                "cell_id": cell.cell_id,
                "seed": cell.seed,
                "leg": cell.leg,
                "replicate_index": cell.replicate_index,
            },
            "platform": {
                "name": "Reference",
                "openmm_version": mm.__version__,
                "python": sys.executable,
                "md_steps": 0,
                "minimization_steps": 0,
                "gpu_used": False,
                "free_energy_estimated": False,
            },
            "inventory_digest": inventory["inventory_digest"],
            "source_summary_digest": source_summary["source_summary_digest"],
            "source_cell_digest": frozen_cell["cell_digest"],
            "source_result": artifact_record(cell.build_result_path),
            "source_artifacts": source_result["artifacts"],
            "implementation": inventory["implementation"],
            "bridge_system": None,
            "bridge_contract": None,
            "ghost_contract": parent_p0.compact_force_report(source_ghost_report),
            "states": {
                "source_apex_dplus": states["apex_dplus"],
                "source_apex_dminus": states["apex_dminus"],
                "bridge_xi": [],
                "ghost_g": 0.5,
            },
            "raw_endpoints": raw_endpoints,
            "raw_guard": raw_guard,
            "evaluations": {
                "source": {
                    name: p0v2.evaluation_record(value)
                    for name, value in source_values.items()
                },
                "bridge": {},
                "roundtrip": {},
            },
            "endpoint_gap": None,
            "bridge_checks": None,
            "stop_reason": "raw union-source gross-clash guard failed before bridge",
            "elapsed_s": time.monotonic() - started,
        }
        del source_system, source_ghost_report
        gc.collect()
        atomic_write_json(reference_result_path(cell), result)
        write_log(log_path, f"COMPLETE {cell.cell_id} status=REJECT_UNION_SOURCE")
        return result
    del source_system, source_ghost_report
    gc.collect()

    write_log(log_path, "append unchanged ghost and exact apex bridge")
    bridge_system = mm.XmlSerializer.deserialize(source_xml)
    ghost_report = dg.build_dynamic_ghost_force(
        bridge_system, pdb.topology, ghost_protocol
    )
    bridge_report = bridge.build_apex_bridge(bridge_system, bridge_protocol)
    ghost_before_serialization = dg.inspect_dynamic_ghost_force(bridge_system)

    serialized = mm.XmlSerializer.serialize(bridge_system)
    bridge_xml_path = reference_bridge_xml_path(cell)
    atomic_write_text(bridge_xml_path, serialized)
    bridge_xml_hash = sha256_text(serialized)

    write_log(log_path, "evaluate bridge xi states and inactive globals")
    bridge_values, inactive_values = p0v2.evaluate_bridge(
        bridge_system,
        pdb.positions,
        states["apex_dplus"],
        xi_values,
        check_inactive_globals=True,
    )
    del bridge_system
    gc.collect()

    write_log(log_path, "deserialize bridge and evaluate roundtrip parity")
    roundtrip_system = mm.XmlSerializer.deserialize(serialized)
    bridge_roundtrip_contract = bridge.inspect_apex_bridge(roundtrip_system)
    ghost_after_serialization = dg.inspect_dynamic_ghost_force(roundtrip_system)
    roundtrip_values, _unused = p0v2.evaluate_bridge(
        roundtrip_system,
        pdb.positions,
        states["apex_dplus"],
        [0.0, 0.5, 1.0],
        check_inactive_globals=False,
    )
    del roundtrip_system, serialized, source_xml
    gc.collect()

    plus = source_values["apex_dplus"]
    minus = source_values["apex_dminus"]
    gap_energy = minus["energy_kj_mol"] - plus["energy_kj_mol"]
    gap_forces = minus["forces"] - plus["forces"]
    endpoint_parities = {
        "xi0_vs_source_dplus": p0v2.parity_record(
            bridge_values["0.00"]["energy_kj_mol"],
            bridge_values["0.00"]["forces"],
            plus["energy_kj_mol"],
            plus["forces"],
        ),
        "xi1_vs_source_dminus": p0v2.parity_record(
            bridge_values["1.00"]["energy_kj_mol"],
            bridge_values["1.00"]["forces"],
            minus["energy_kj_mol"],
            minus["forces"],
        ),
    }
    interior_parities: dict[str, dict[str, Any]] = {}
    derivative_parities: dict[str, dict[str, Any]] = {}
    for xi in (0.25, 0.5, 0.75):
        key = f"{xi:.2f}"
        expected_energy = (1.0 - xi) * plus["energy_kj_mol"] + xi * minus[
            "energy_kj_mol"
        ]
        expected_forces = (1.0 - xi) * plus["forces"] + xi * minus["forces"]
        interior_parities[key] = p0v2.parity_record(
            bridge_values[key]["energy_kj_mol"],
            bridge_values[key]["forces"],
            expected_energy,
            expected_forces,
        )
        derivative_parities[key] = p0v2.derivative_record(
            bridge_values[key]["derivatives_kj_mol"][bridge.BRIDGE_PARAMETER],
            gap_energy,
        )

    ghost_rows = [plus, minus] + [bridge_values[f"{xi:.2f}"] for xi in xi_values]
    ghost_energy_reference = ghost_rows[0]["ghost_energy_kj_mol"]
    dvdg_reference = ghost_rows[0]["derivatives_kj_mol"][dg.GHOST_GLOBAL_PARAMETER]
    ghost_invariance = {
        "max_abs_ghost_energy_delta_kj_mol": max(
            abs(row["ghost_energy_kj_mol"] - ghost_energy_reference)
            for row in ghost_rows
        ),
        "max_abs_dVdg_delta_kj_mol": max(
            abs(
                row["derivatives_kj_mol"][dg.GHOST_GLOBAL_PARAMETER]
                - dvdg_reference
            )
            for row in ghost_rows
        ),
    }
    ghost_invariance["passed"] = bool(
        ghost_invariance["max_abs_ghost_energy_delta_kj_mol"] <= ENERGY_TOLERANCE
        and ghost_invariance["max_abs_dVdg_delta_kj_mol"] <= DERIVATIVE_TOLERANCE
    )

    if inactive_values is None:
        raise RuntimeError("inactive-global audit was not evaluated")
    inactive_parity = p0v2.parity_record(
        inactive_values["changed"]["energy_kj_mol"],
        inactive_values["changed"]["forces"],
        inactive_values["baseline"]["energy_kj_mol"],
        inactive_values["baseline"]["forces"],
    )
    roundtrip_parities = {
        key: p0v2.parity_record(
            roundtrip_values[key]["energy_kj_mol"],
            roundtrip_values[key]["forces"],
            bridge_values[key]["energy_kj_mol"],
            bridge_values[key]["forces"],
        )
        for key in ("0.00", "0.50", "1.00")
    }
    finite_all = all(
        row["finite"]
        for row in list(source_values.values())
        + list(bridge_values.values())
        + list(roundtrip_values.values())
        + [inactive_values["baseline"], inactive_values["changed"]]
    )
    endpoint_gap = {
        "source_hminus_minus_hplus_kj_mol": gap_energy,
        "bridge_xi1_minus_xi0_kj_mol": (
            bridge_values["1.00"]["energy_kj_mol"]
            - bridge_values["0.00"]["energy_kj_mol"]
        ),
        "max_abs_force_component_gap_kj_mol_nm": float(
            np.max(np.abs(gap_forces))
        ),
        "is_pass_fail_threshold": False,
    }
    endpoint_gap["abs_gap_identity_delta_kj_mol"] = abs(
        endpoint_gap["bridge_xi1_minus_xi0_kj_mol"] - gap_energy
    )

    bridge_checks = {
        "bridge_mutation_scope": all(bridge_report["checks"].values()),
        "dynamic_ghost_contract_roundtrip": (
            ghost_before_serialization == ghost_after_serialization
        ),
        "bridge_contract_roundtrip": (
            bridge_roundtrip_contract == bridge_report["after"]
        ),
        "serialized_xml_hash_matches_disk": (
            sha256_file(bridge_xml_path) == bridge_xml_hash
        ),
        "finite_energy_forces_derivatives": finite_all,
        "endpoint_parity": all(row["passed"] for row in endpoint_parities.values()),
        "interior_linear_identity": all(
            row["passed"] for row in interior_parities.values()
        ),
        "interior_derivative_identity": all(
            row["passed"] for row in derivative_parities.values()
        ),
        "ghost_xi_invariance": ghost_invariance["passed"],
        "inactive_legacy_global_invariance": inactive_parity["passed"],
        "serialization_energy_force_parity": all(
            row["passed"] for row in roundtrip_parities.values()
        ),
        "endpoint_gap_identity": (
            endpoint_gap["abs_gap_identity_delta_kj_mol"] <= ENERGY_TOLERANCE
        ),
        "execution_cpu_reference_only": True,
        "zero_md_and_minimization_steps": True,
        "free_energy_not_estimated": True,
    }
    if not raw_guard["passed"]:
        status = "REJECT_UNION_SOURCE"
    elif all(bridge_checks.values()):
        status = "P0V3R2_PASS"
    else:
        status = "REJECT_BRIDGE"
    result = {
        "schema": "updd_dynamic_ghost_union_bridge_p0v3r2_cell_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "cell": {
            "cell_id": cell.cell_id,
            "seed": cell.seed,
            "leg": cell.leg,
            "replicate_index": cell.replicate_index,
        },
        "platform": {
            "name": "Reference",
            "openmm_version": mm.__version__,
            "python": sys.executable,
            "md_steps": 0,
            "minimization_steps": 0,
            "gpu_used": False,
            "free_energy_estimated": False,
        },
        "inventory_digest": inventory["inventory_digest"],
        "source_summary_digest": source_summary["source_summary_digest"],
        "source_cell_digest": frozen_cell["cell_digest"],
        "source_result": artifact_record(cell.build_result_path),
        "source_artifacts": source_result["artifacts"],
        "implementation": inventory["implementation"],
        "bridge_system": artifact_record(bridge_xml_path),
        "bridge_contract": p0v2._compact_bridge_report(bridge_report),
        "ghost_contract": parent_p0.compact_force_report(ghost_report),
        "states": {
            "source_apex_dplus": states["apex_dplus"],
            "source_apex_dminus": states["apex_dminus"],
            "bridge_xi": xi_values,
            "ghost_g": 0.5,
        },
        "raw_endpoints": raw_endpoints,
        "raw_guard": raw_guard,
        "evaluations": {
            "source": {
                name: p0v2.evaluation_record(value)
                for name, value in source_values.items()
            },
            "bridge": {
                name: p0v2.evaluation_record(value)
                for name, value in bridge_values.items()
            },
            "roundtrip": {
                name: p0v2.evaluation_record(value)
                for name, value in roundtrip_values.items()
            },
        },
        "endpoint_gap": endpoint_gap,
        "endpoint_parities": endpoint_parities,
        "interior_parities": interior_parities,
        "derivative_parities": derivative_parities,
        "ghost_invariance": ghost_invariance,
        "inactive_global_parity": {
            **inactive_parity,
            "replacements": inactive_values["replacements"],
        },
        "roundtrip_parities": roundtrip_parities,
        "bridge_checks": bridge_checks,
        "elapsed_s": time.monotonic() - started,
    }
    atomic_write_json(reference_result_path(cell), result)
    write_log(log_path, f"COMPLETE {cell.cell_id} status={status}")
    return result


def validate_reference_result(
    cell: UnionCell,
    inventory: dict[str, Any],
    source_summary: dict[str, Any],
) -> dict[str, Any] | None:
    path = reference_result_path(cell)
    if not path.is_file():
        return None
    try:
        result = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    expected_cell = {
        "cell_id": cell.cell_id,
        "seed": cell.seed,
        "leg": cell.leg,
        "replicate_index": cell.replicate_index,
    }
    if result.get("cell") != expected_cell:
        raise ValueError(f"{cell.cell_id}: Reference result identity mismatch")
    if result.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError(f"{cell.cell_id}: Reference inventory mismatch")
    if result.get("source_summary_digest") != source_summary["source_summary_digest"]:
        raise ValueError(f"{cell.cell_id}: Reference source-summary mismatch")
    if result.get("source_cell_digest") != frozen_cell_record(
        inventory, cell
    )["cell_digest"]:
        raise ValueError(f"{cell.cell_id}: Reference source cell mismatch")
    if result.get("source_result") != artifact_record(cell.build_result_path):
        raise ValueError(f"{cell.cell_id}: Reference source-result hash mismatch")
    if result.get("status") not in {
        "P0V3R2_PASS",
        "REJECT_UNION_SOURCE",
        "REJECT_BRIDGE",
    }:
        raise ValueError(f"{cell.cell_id}: invalid Reference result status")
    raw_passed = (result.get("raw_guard") or {}).get("passed") is True
    if result["status"] == "REJECT_UNION_SOURCE" and raw_passed:
        raise ValueError(f"{cell.cell_id}: source reject has a passing raw guard")
    if result["status"] == "REJECT_BRIDGE":
        if not raw_passed or all((result.get("bridge_checks") or {}).values()):
            raise ValueError(f"{cell.cell_id}: bridge reject classification drifted")
    if result["status"] == "P0V3R2_PASS":
        if not raw_passed or not all((result.get("bridge_checks") or {}).values()):
            raise ValueError(f"{cell.cell_id}: P0V3R2_PASS has a false gate")
    bridge_system = result.get("bridge_system")
    if result["status"] != "REJECT_UNION_SOURCE":
        if not isinstance(bridge_system, dict):
            raise ValueError(f"{cell.cell_id}: bridge XML declaration is absent")
        bridge_path = REPO_ROOT / bridge_system.get("path", "")
        if not bridge_path.is_file() or sha256_file(bridge_path) != bridge_system.get(
            "sha256"
        ):
            raise ValueError(f"{cell.cell_id}: bridge XML missing or drifted")
    return result


def _write_reference_summary(
    results: list[dict[str, Any]],
    inventory: dict[str, Any],
    source_summary: dict[str, Any],
) -> dict[str, Any]:
    cells = make_cells()
    by_id = {row["cell"]["cell_id"]: row for row in results}
    pending = [cell.cell_id for cell in cells if cell.cell_id not in by_id]
    source_rejected = [
        cell_id
        for cell_id, row in by_id.items()
        if row.get("status") == "REJECT_UNION_SOURCE"
    ]
    bridge_rejected = [
        cell_id
        for cell_id, row in by_id.items()
        if row.get("status") == "REJECT_BRIDGE"
    ]
    if source_rejected:
        status = "REJECT_UNION_SOURCE"
    elif bridge_rejected:
        status = "REJECT_BRIDGE"
    elif not pending and len(results) == len(cells):
        status = "P0V3R2_PASS"
    else:
        status = "RUNNING"
    elapsed = [float(row["elapsed_s"]) for row in results]
    median_elapsed = median(elapsed) if elapsed else None
    payload = {
        "schema": "updd_dynamic_ghost_union_bridge_p0v3r2_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "platform": "Reference",
        "gpu_used": False,
        "md_steps": 0,
        "minimization_steps": 0,
        "free_energy_estimated": False,
        "inventory_digest": inventory["inventory_digest"],
        "source_summary_digest": source_summary["source_summary_digest"],
        "n_expected": len(cells),
        "n_completed": len(results),
        "pending_cell_ids": pending,
        "source_rejected_cell_ids": source_rejected,
        "bridge_rejected_cell_ids": bridge_rejected,
        "median_elapsed_s": median_elapsed,
        "eta_s": median_elapsed * len(pending) if median_elapsed is not None else None,
        "cells": [
            {
                "cell_id": row["cell"]["cell_id"],
                "seed": row["cell"]["seed"],
                "leg": row["cell"]["leg"],
                "status": row["status"],
                "elapsed_s": row["elapsed_s"],
                "raw_u0_kj_mol": row["raw_endpoints"]["u0_kj_mol"],
                "raw_u1_kj_mol": row["raw_endpoints"]["u1_kj_mol"],
                "endpoint_gap_kj_mol": (
                    row["endpoint_gap"]["source_hminus_minus_hplus_kj_mol"]
                    if row.get("endpoint_gap") is not None
                    else None
                ),
                "max_force_gap_kj_mol_nm": (
                    row["endpoint_gap"]["max_abs_force_component_gap_kj_mol_nm"]
                    if row.get("endpoint_gap") is not None
                    else None
                ),
                "result": artifact_record(
                    reference_result_path(cell_by_id(row["cell"]["cell_id"]))
                ),
                "bridge_system": row.get("bridge_system"),
            }
            for row in results
        ],
    }
    atomic_write_json(REFERENCE_SUMMARY_JSON, payload)
    lines = [
        "# Union Apex Bridge P0v3r2 Reference Summary",
        "",
        f"Generated: {payload['generated']}",
        "",
        f"Status: **{status}**",
        "",
        "P0v3r2 used OpenMM Reference with zero MD/minimization steps and no GPU. "
        "It did not estimate a free energy or claim equilibrium density.",
        "",
        "| cell | status | raw u0 | raw u1 | bridge endpoint gap | elapsed (s) |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in payload["cells"]:
        gap = (
            f"{row['endpoint_gap_kj_mol']:.12f}"
            if row["endpoint_gap_kj_mol"] is not None
            else "NA"
        )
        lines.append(
            f"| {row['cell_id']} | {row['status']} | {row['raw_u0_kj_mol']:.6f} | "
            f"{row['raw_u1_kj_mol']:.6f} | {gap} | "
            f"{row['elapsed_s']:.3f} |"
        )
    if pending:
        lines.extend(["", "Pending: " + ", ".join(f"`{item}`" for item in pending)])
    if source_rejected:
        lines.extend(
            ["", "Union-source reject: " + ", ".join(f"`{item}`" for item in source_rejected)]
        )
    if bridge_rejected:
        lines.extend(
            ["", "Bridge reject: " + ", ".join(f"`{item}`" for item in bridge_rejected)]
        )
    atomic_write_text(REFERENCE_SUMMARY_MD, "\n".join(lines) + "\n")
    return payload


def run_reference_parent(
    inventory: dict[str, Any], source_summary: dict[str, Any]
) -> dict[str, Any]:
    validate_keeper_preflight(inventory, source_summary)
    REFERENCE_ROOT.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for cell in make_cells():
        completed = validate_reference_result(cell, inventory, source_summary)
        if completed is not None:
            results.append(completed)
            _write_reference_summary(results, inventory, source_summary)
            if completed["status"] != "P0V3R2_PASS":
                return _write_reference_summary(results, inventory, source_summary)
            continue
        output_dir = reference_cell_dir(cell)
        if output_dir.exists():
            _archive_partial(
                output_dir,
                REFERENCE_ROOT / "_archive",
                "incomplete or invalid P0v3r2 Reference cell",
            )
        command = [
            sys.executable,
            str(RUNNER_PATH),
            "--worker-reference-cell",
            cell.cell_id,
            "--expected-inventory-digest",
            inventory["inventory_digest"],
            "--expected-source-summary-digest",
            source_summary["source_summary_digest"],
        ]
        returncode = subprocess.run(command, cwd=REPO_ROOT, check=False).returncode
        if returncode != 0:
            _write_reference_summary(results, inventory, source_summary)
            raise RuntimeError(
                f"Reference worker {cell.cell_id} exited with status {returncode}"
            )
        completed = validate_reference_result(cell, inventory, source_summary)
        if completed is None:
            raise RuntimeError(f"{cell.cell_id}: Reference worker produced no result")
        results.append(completed)
        _write_reference_summary(results, inventory, source_summary)
        if completed["status"] != "P0V3R2_PASS":
            return _write_reference_summary(results, inventory, source_summary)
    return _write_reference_summary(results, inventory, source_summary)


def _assert_atm_interpreter() -> None:
    if Path(sys.executable).resolve() != EXPECTED_ATM_PYTHON.resolve():
        raise RuntimeError(
            "P0v3r2 must run with the Track B ATM interpreter: "
            f"{EXPECTED_ATM_PYTHON}; observed {sys.executable}"
        )


def _write_run_state(status: str, **extra: Any) -> None:
    atomic_write_json(
        RUN_STATE,
        {
            "status": status,
            "updated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "python": sys.executable,
            **extra,
        },
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--verify-freeze", action="store_true")
    group.add_argument("--prepare-inventory", action="store_true")
    group.add_argument("--build-sources", action="store_true")
    group.add_argument("--keeper-preflight", action="store_true")
    group.add_argument("--run", action="store_true")
    group.add_argument("--summarize", action="store_true")
    group.add_argument("--worker-build-cell")
    group.add_argument("--worker-reference-cell")
    parser.add_argument("--expected-inventory-digest")
    parser.add_argument("--expected-source-summary-digest")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        _assert_atm_interpreter()
        if args.verify_freeze:
            _protocol, base, _bridge_protocol, _ghost_protocol = (
                load_protocol_and_verify_freeze()
            )
            print(
                json.dumps(
                    {
                        "status": "PASS",
                        "cells": [cell.cell_id for cell in make_cells(base)],
                        "python": sys.executable,
                    },
                    indent=2,
                )
            )
            return 0
        if args.prepare_inventory:
            payload = build_prebuild_inventory()
            if PREBUILD_INVENTORY.exists():
                existing = load_prebuild_inventory()
                if existing["inventory_digest"] != payload["inventory_digest"]:
                    raise ValueError(
                        "refusing to overwrite a drifted P0v3r2 prebuild inventory"
                    )
            else:
                atomic_write_json(PREBUILD_INVENTORY, payload)
            _write_run_state(
                "PREBUILD_INVENTORY_FROZEN",
                inventory_digest=payload["inventory_digest"],
            )
            print(json.dumps(payload, indent=2))
            return 0

        inventory = load_prebuild_inventory(args.expected_inventory_digest)
        if args.worker_build_cell:
            cell = cell_by_id(args.worker_build_cell)
            result = run_source_worker(cell, inventory)
            print(json.dumps({"cell": cell.cell_id, "status": result["status"]}))
            return 0
        if args.worker_reference_cell:
            source_summary = load_source_summary_shallow(inventory)
            if (
                args.expected_source_summary_digest is None
                or source_summary["source_summary_digest"]
                != args.expected_source_summary_digest
            ):
                raise ValueError("Reference worker source summary digest mismatch")
            cell = cell_by_id(args.worker_reference_cell)
            result = run_reference_worker(cell, inventory, source_summary)
            print(json.dumps({"cell": cell.cell_id, "status": result["status"]}))
            return 0

        inventory = verify_prebuild_inventory()
        if args.build_sources:
            _write_run_state(
                "BUILDING_UNION_SOURCES",
                inventory_digest=inventory["inventory_digest"],
            )
            summary = run_source_parent(inventory)
            _write_run_state(
                summary["status"],
                stage="U0",
                inventory_digest=inventory["inventory_digest"],
                source_summary_digest=summary["source_summary_digest"],
            )
            print(json.dumps(summary, indent=2))
            return 0
        if args.keeper_preflight:
            _write_run_state(
                "KEEPER_PREFLIGHT",
                inventory_digest=inventory["inventory_digest"],
            )
            payload = run_keeper_preflight(inventory)
            _write_run_state(
                "KEEPER_PASS",
                inventory_digest=inventory["inventory_digest"],
                source_summary_digest=payload["source_summary_digest"],
            )
            print(json.dumps(payload, indent=2))
            return 0

        source_summary = load_source_summary(inventory)
        if args.run:
            _write_run_state(
                "REFERENCE_RUNNING",
                inventory_digest=inventory["inventory_digest"],
                source_summary_digest=source_summary["source_summary_digest"],
            )
            summary = run_reference_parent(inventory, source_summary)
            _write_run_state(
                summary["status"],
                stage="P0v3r2",
                inventory_digest=inventory["inventory_digest"],
                source_summary_digest=source_summary["source_summary_digest"],
            )
            print(json.dumps(summary, indent=2))
            return 0
        if args.summarize:
            results = []
            for cell in make_cells():
                result = validate_reference_result(cell, inventory, source_summary)
                if result is not None:
                    results.append(result)
            summary = _write_reference_summary(results, inventory, source_summary)
            print(json.dumps(summary, indent=2))
            return 0
        raise RuntimeError("unhandled P0v3r2 command")
    except Exception as exc:
        _write_run_state(
            "P0V3R2_INTERRUPTED",
            error_type=type(exc).__name__,
            error_message=str(exc),
        )
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
