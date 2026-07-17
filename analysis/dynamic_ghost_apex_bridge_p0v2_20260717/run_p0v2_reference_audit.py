#!/usr/bin/env python3
"""Run the frozen apex-bridge P0v2 audit on OpenMM Reference."""

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
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
UTILS_DIR = REPO_ROOT / "utils"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))

from analysis.dynamic_ghost_excluded_volume_20260716 import (  # noqa: E402
    run_p0_reference_audit as parent_p0,
)
import trackb_apex_bridge as bridge  # noqa: E402
import trackb_dynamic_ghost as dg  # noqa: E402


ANALYSIS_DIR = Path(__file__).resolve().parent
PREREGISTRATION = ANALYSIS_DIR / "PREREGISTRATION.md"
PROTOCOL_PATH = ANALYSIS_DIR / "protocol.json"
FROZEN_MANIFEST = ANALYSIS_DIR / "FROZEN_MANIFEST.json"

P0V2_INVENTORY = ANALYSIS_DIR / "P0V2_SOURCE_INVENTORY.json"
P0V2_KEEPER_PREFLIGHT = ANALYSIS_DIR / "P0V2_KEEPER_PREFLIGHT.json"
P0V2_ROOT = ANALYSIS_DIR / "p0v2_reference_audit"
P0V2_SUMMARY_JSON = ANALYSIS_DIR / "P0V2_REFERENCE_SUMMARY.json"
P0V2_SUMMARY_MD = ANALYSIS_DIR / "P0V2_REFERENCE_SUMMARY.md"

BRIDGE_IMPLEMENTATION = REPO_ROOT / "utils/trackb_apex_bridge.py"
GHOST_IMPLEMENTATION = REPO_ROOT / "utils/trackb_dynamic_ghost.py"
RUNNER_PATH = Path(__file__).resolve()

ENERGY_TOLERANCE = 1.0e-5
FORCE_TOLERANCE = 1.0e-5
DERIVATIVE_TOLERANCE = 1.0e-5
KCAL_TO_KJ = 4.184

Cell = parent_p0.Cell
make_cells = parent_p0.make_cells


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_digest(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


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


def load_protocol_and_verify_freeze() -> tuple[dict[str, Any], dict[str, Any]]:
    frozen = json.loads(FROZEN_MANIFEST.read_text(encoding="utf-8"))
    expected_status = "FROZEN_BEFORE_BRIDGE_IMPLEMENTATION_OR_P0V2_OUTPUT"
    if frozen.get("status") != expected_status:
        raise ValueError("P0v2 frozen manifest status drifted")
    for name, expected_hash in (frozen.get("files") or {}).items():
        path = ANALYSIS_DIR / name
        observed = sha256_file(path)
        if observed != expected_hash:
            raise ValueError(
                f"frozen P0v2 preregistration drift: {path} "
                f"{observed} != {expected_hash}"
            )
    for path_text, expected_hash in (frozen.get("trigger_evidence") or {}).items():
        path = REPO_ROOT / path_text
        observed = sha256_file(path)
        if observed != expected_hash:
            raise ValueError(
                f"P0v2 trigger evidence drift: {path} "
                f"{observed} != {expected_hash}"
            )

    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    bridge.validate_bridge_protocol(protocol)
    trigger = protocol.get("trigger") or {}
    parent_summary = json.loads(
        (REPO_ROOT / "analysis/dynamic_ghost_excluded_volume_20260716/"
         "P0_REFERENCE_SUMMARY.json").read_text(encoding="utf-8")
    )
    if parent_summary.get("status") != "REJECT_GHOST":
        raise ValueError("parent P0 is no longer the frozen REJECT_GHOST trigger")
    if parent_summary.get("inventory_digest") != trigger.get("parent_inventory_digest"):
        raise ValueError("parent P0 inventory digest drifted")
    rejected = parent_summary.get("rejected_cell_ids") or []
    if rejected != [trigger.get("failed_cell")]:
        raise ValueError("parent P0 failed-cell identity drifted")

    parent_protocol, _parent_frozen = parent_p0.load_protocol_and_verify_freeze()
    dg.validate_ghost_protocol(parent_protocol)
    return protocol, parent_protocol


def inventory_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("inventory_digest", None)
    return canonical_digest(canonical)


def build_source_inventory() -> dict[str, Any]:
    protocol, parent_protocol = load_protocol_and_verify_freeze()
    frozen = json.loads(FROZEN_MANIFEST.read_text(encoding="utf-8"))
    payload = {
        "schema": "updd_dynamic_ghost_apex_bridge_p0v2_inventory_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "platform": "Reference",
        "md_allowed": False,
        "minimization_allowed": False,
        "gpu_allowed": False,
        "free_energy_estimation_allowed": False,
        "bridge_states": (protocol.get("bridge") or {}).get("audit_states"),
        "ghost_g": (protocol.get("bridge") or {}).get("ghost_g"),
        "frozen_manifest": artifact_record(FROZEN_MANIFEST),
        "frozen_files": {
            name: artifact_record(ANALYSIS_DIR / name)
            for name in sorted((frozen.get("files") or {}))
        },
        "trigger_evidence": {
            path_text: artifact_record(REPO_ROOT / path_text)
            for path_text in sorted((frozen.get("trigger_evidence") or {}))
        },
        "implementation": {
            "bridge_module": artifact_record(BRIDGE_IMPLEMENTATION),
            "ghost_module": artifact_record(GHOST_IMPLEMENTATION),
            "p0v2_runner": artifact_record(RUNNER_PATH),
        },
        "source_preregistrations": {
            leg: parent_p0.validate_leg_preregistration(leg)
            for leg in parent_p0.LEGS
        },
        "cells": [
            parent_p0.validate_source_cell(cell, parent_protocol)
            for cell in make_cells()
        ],
    }
    payload["inventory_digest"] = inventory_digest(payload)
    return payload


def load_frozen_inventory(expected_digest: str | None = None) -> dict[str, Any]:
    if not P0V2_INVENTORY.is_file():
        raise FileNotFoundError(
            f"missing {P0V2_INVENTORY}; run --prepare-inventory first"
        )
    inventory = json.loads(P0V2_INVENTORY.read_text(encoding="utf-8"))
    if inventory.get("inventory_digest") != inventory_digest(inventory):
        raise ValueError("stored P0v2 inventory digest is invalid")
    if expected_digest is not None and inventory["inventory_digest"] != expected_digest:
        raise ValueError(
            f"worker inventory mismatch: {inventory['inventory_digest']} "
            f"!= {expected_digest}"
        )
    return inventory


def verify_source_inventory() -> dict[str, Any]:
    frozen = load_frozen_inventory()
    current = build_source_inventory()
    if frozen["inventory_digest"] != current["inventory_digest"]:
        raise ValueError(
            "P0v2 source/code drift after inventory freeze: "
            f"{frozen['inventory_digest']} != {current['inventory_digest']}"
        )
    return frozen


def frozen_cell_record(inventory: dict[str, Any], cell: Cell) -> dict[str, Any]:
    rows = [row for row in inventory.get("cells", []) if row.get("cell_id") == cell.cell_id]
    if len(rows) != 1:
        raise ValueError(f"{cell.cell_id}: missing or duplicate frozen inventory row")
    return rows[0]


def cell_by_id(cell_id: str) -> Cell:
    rows = [cell for cell in make_cells() if cell.cell_id == cell_id]
    if len(rows) != 1:
        raise ValueError(f"expected one P0v2 cell {cell_id!r}, found {len(rows)}")
    return rows[0]


def verify_worker_inputs(
    inventory: dict[str, Any], cell: Cell
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    protocol, parent_protocol = load_protocol_and_verify_freeze()
    expected_code = inventory.get("implementation") or {}
    current_code = {
        "bridge_module": artifact_record(BRIDGE_IMPLEMENTATION),
        "ghost_module": artifact_record(GHOST_IMPLEMENTATION),
        "p0v2_runner": artifact_record(RUNNER_PATH),
    }
    if current_code != expected_code:
        raise ValueError(f"{cell.cell_id}: implementation drift after inventory freeze")
    current_cell = parent_p0.validate_source_cell(cell, parent_protocol)
    frozen_cell = frozen_cell_record(inventory, cell)
    if current_cell != frozen_cell:
        raise ValueError(f"{cell.cell_id}: source cell drift after inventory freeze")
    return protocol, parent_protocol, frozen_cell


def _compact_bridge_report(report: dict[str, Any]) -> dict[str, Any]:
    return {
        "atm_index": report["atm_index"],
        "checks": report["checks"],
        "before": report["before"],
        "after": report["after"],
    }


def run_keeper_preflight(inventory: dict[str, Any]) -> dict[str, Any]:
    """Validate all six source-to-bridge declarations without a Context."""
    import openmm as mm
    from openmm.app import PDBFile

    rows: list[dict[str, Any]] = []
    for cell in make_cells():
        protocol, parent_protocol, frozen_cell = verify_worker_inputs(inventory, cell)
        system = mm.XmlSerializer.deserialize(cell.xml_path.read_text(encoding="utf-8"))
        pdb = PDBFile(str(cell.pdb_path))
        ghost_report = dg.build_dynamic_ghost_force(system, pdb.topology, parent_protocol)
        bridge_report = bridge.build_apex_bridge(system, protocol)
        ghost_after = dg.inspect_dynamic_ghost_force(system)
        bridge_after = bridge.inspect_apex_bridge(system)
        checks = {
            "bridge_mutation_scope_exact": all(bridge_report["checks"].values()),
            "exact_nine_ring_atoms": ghost_report["selection"]["n_ring_atoms"] == 9,
            "all_water_oxygens": (
                ghost_report["selection"]["n_water_oxygens"]
                == ghost_report["selection"]["n_water_residues"]
                and ghost_report["selection"]["n_water_oxygens"] > 0
            ),
            "atom_count_invariant": (
                ghost_report["n_particles_before"] == ghost_report["n_particles_after"]
            ),
            "charge_invariant": ghost_report["net_charge_delta_e"] == 0.0,
            "ghost_after_atm": (
                ghost_after["force_index"] > bridge_after["atm_index"]
            ),
            "ghost_group_31": ghost_after["force_group"] == dg.GHOST_FORCE_GROUP,
            "bridge_parameter_consumed": (
                bridge.BRIDGE_PARAMETER in bridge_after["energy_function"]
                and bridge_after["energy_parameter_derivatives"]
                == [bridge.BRIDGE_PARAMETER]
            ),
        }
        if not all(checks.values()):
            raise bridge.ApexBridgeError(
                f"{cell.cell_id}: KEEPER declaration checks failed: {checks}"
            )
        rows.append(
            {
                "cell_id": cell.cell_id,
                "seed": cell.seed,
                "leg": cell.leg,
                "source_cell_digest": frozen_cell["cell_digest"],
                "checks": checks,
                "bridge_contract": _compact_bridge_report(bridge_report),
                "ghost_contract": parent_p0.compact_force_report(ghost_report),
            }
        )
        del system, pdb, ghost_report, bridge_report, bridge_after, ghost_after
        gc.collect()

    payload = {
        "schema": "updd_dynamic_ghost_apex_bridge_keeper_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "PASS",
        "context_count": 0,
        "md_steps": 0,
        "gpu_used": False,
        "inventory_digest": inventory["inventory_digest"],
        "n_cells": len(rows),
        "cells": rows,
    }
    atomic_write_json(P0V2_KEEPER_PREFLIGHT, payload)
    return payload


def validate_keeper_preflight(inventory: dict[str, Any]) -> dict[str, Any]:
    if not P0V2_KEEPER_PREFLIGHT.is_file():
        raise FileNotFoundError(
            f"missing {P0V2_KEEPER_PREFLIGHT}; run --keeper-preflight first"
        )
    payload = json.loads(P0V2_KEEPER_PREFLIGHT.read_text(encoding="utf-8"))
    if payload.get("status") != "PASS" or payload.get("context_count") != 0:
        raise ValueError("P0v2 KEEPER preflight is not a context-free PASS")
    if payload.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError("P0v2 KEEPER preflight inventory digest drifted")
    if payload.get("n_cells") != len(make_cells()):
        raise ValueError("P0v2 KEEPER preflight cell count drifted")
    return payload


def cell_output_dir(cell: Cell) -> Path:
    return P0V2_ROOT / cell.cell_id


def cell_result_path(cell: Cell) -> Path:
    return cell_output_dir(cell) / "cell_result.json"


def cell_bridge_xml_path(cell: Cell) -> Path:
    return cell_output_dir(cell) / "apex_bridge_system.xml"


def archive_partial_cell(cell: Cell) -> Path | None:
    root = cell_output_dir(cell)
    if not root.exists():
        return None
    stamp = time.strftime("%Y%m%d_%H%M%S")
    target = P0V2_ROOT / "_archive" / f"{stamp}_{cell.cell_id}"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(root), str(target))
    atomic_write_json(
        target / "archive_reason.json",
        {
            "cell_id": cell.cell_id,
            "archived": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "reason": "incomplete or invalid P0v2 cell output",
        },
    )
    return target


def _set_source_parameters(context: Any, state: dict[str, Any], ghost_g: float) -> None:
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
    context.setParameter(dg.GHOST_GLOBAL_PARAMETER, float(ghost_g))


def _read_context(context: Any, *, expect_bridge: bool) -> dict[str, Any]:
    import numpy as np
    from openmm import unit

    state = context.getState(
        getEnergy=True,
        getForces=True,
        getParameterDerivatives=True,
    )
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    forces = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=float,
    )
    derivatives = {
        str(name): float(value)
        for name, value in dict(state.getEnergyParameterDerivatives()).items()
    }
    ghost_energy = float(
        context.getState(
            getEnergy=True,
            groups=1 << dg.GHOST_FORCE_GROUP,
        ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    )
    required = [dg.GHOST_GLOBAL_PARAMETER]
    if expect_bridge:
        required.append(bridge.BRIDGE_PARAMETER)
    finite = bool(
        math.isfinite(energy)
        and np.isfinite(forces).all()
        and math.isfinite(ghost_energy)
        and all(name in derivatives and math.isfinite(derivatives[name]) for name in required)
    )
    return {
        "energy_kj_mol": energy,
        "forces": forces,
        "ghost_energy_kj_mol": ghost_energy,
        "derivatives_kj_mol": derivatives,
        "finite": finite,
    }


def evaluate_source_apexes(
    system: Any,
    positions: Any,
    states: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    import openmm as mm
    from openmm import unit

    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    results: dict[str, dict[str, Any]] = {}
    for name in ("apex_dplus", "apex_dminus"):
        _set_source_parameters(context, states[name], 0.5)
        results[name] = _read_context(context, expect_bridge=False)
    del context, integrator
    gc.collect()
    return results


def evaluate_bridge(
    system: Any,
    positions: Any,
    state: dict[str, Any],
    xi_values: list[float],
    *,
    check_inactive_globals: bool,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any] | None]:
    import openmm as mm
    from openmm import unit

    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    _set_source_parameters(context, state, 0.5)
    results: dict[str, dict[str, Any]] = {}
    for xi in xi_values:
        context.setParameter(bridge.BRIDGE_PARAMETER, float(xi))
        results[f"{xi:.2f}"] = _read_context(context, expect_bridge=True)

    inactive = None
    if check_inactive_globals:
        context.setParameter(bridge.BRIDGE_PARAMETER, 0.5)
        baseline = _read_context(context, expect_bridge=True)
        replacements = {
            "Lambda1": 0.123,
            "Lambda2": 0.876,
            "Alpha": 0.321,
            "Uh": 12.3,
            "W0": -4.5,
            "Direction": -1.0,
        }
        for name, value in replacements.items():
            context.setParameter(name, value)
        changed = _read_context(context, expect_bridge=True)
        inactive = {
            "baseline": baseline,
            "changed": changed,
            "replacements": replacements,
        }
    del context, integrator
    gc.collect()
    return results, inactive


def parity_record(
    observed_energy: float,
    observed_forces: Any,
    expected_energy: float,
    expected_forces: Any,
) -> dict[str, Any]:
    import numpy as np

    energy_delta = abs(float(observed_energy) - float(expected_energy))
    force_delta = float(np.max(np.abs(observed_forces - expected_forces)))
    return {
        "abs_energy_delta_kj_mol": energy_delta,
        "max_abs_force_component_delta_kj_mol_nm": force_delta,
        "energy_tolerance_kj_mol": ENERGY_TOLERANCE,
        "force_tolerance_kj_mol_nm": FORCE_TOLERANCE,
        "passed": bool(
            math.isfinite(energy_delta)
            and math.isfinite(force_delta)
            and energy_delta <= ENERGY_TOLERANCE
            and force_delta <= FORCE_TOLERANCE
        ),
    }


def derivative_record(observed: float, expected: float) -> dict[str, Any]:
    delta = abs(float(observed) - float(expected))
    return {
        "observed_kj_mol": float(observed),
        "expected_kj_mol": float(expected),
        "abs_delta_kj_mol": delta,
        "tolerance_kj_mol": DERIVATIVE_TOLERANCE,
        "passed": bool(math.isfinite(delta) and delta <= DERIVATIVE_TOLERANCE),
    }


def evaluation_record(value: dict[str, Any]) -> dict[str, Any]:
    import numpy as np

    return {
        "energy_kj_mol": value["energy_kj_mol"],
        "max_abs_force_component_kj_mol_nm": float(np.max(np.abs(value["forces"]))),
        "ghost_energy_kj_mol": value["ghost_energy_kj_mol"],
        "derivatives_kj_mol": value["derivatives_kj_mol"],
        "finite": value["finite"],
    }


def write_worker_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S %z')} {message}"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def run_worker(cell: Cell, inventory: dict[str, Any]) -> dict[str, Any]:
    import numpy as np
    import openmm as mm
    from openmm.app import PDBFile

    output_dir = cell_output_dir(cell)
    if output_dir.exists():
        raise RuntimeError(
            f"{cell.cell_id}: output exists; parent must validate or archive it"
        )
    output_dir.mkdir(parents=True)
    log_path = output_dir / "worker.log"
    write_worker_log(log_path, f"START {cell.cell_id} platform=Reference md_steps=0")
    started = time.monotonic()

    protocol, parent_protocol, frozen_cell = verify_worker_inputs(inventory, cell)
    source_xml = cell.xml_path.read_text(encoding="utf-8")
    pdb = PDBFile(str(cell.pdb_path))
    states = frozen_cell["states"]
    xi_values = [float(value) for value in inventory["bridge_states"]]

    write_worker_log(log_path, "build unchanged dynamic ghost on source apex system")
    source_system = mm.XmlSerializer.deserialize(source_xml)
    source_ghost_report = dg.build_dynamic_ghost_force(
        source_system, pdb.topology, parent_protocol
    )
    source_values = evaluate_source_apexes(source_system, pdb.positions, states)
    del source_system, source_ghost_report
    gc.collect()

    write_worker_log(log_path, "build bridge and audit exact mutation scope")
    bridge_system = mm.XmlSerializer.deserialize(source_xml)
    ghost_report = dg.build_dynamic_ghost_force(
        bridge_system, pdb.topology, parent_protocol
    )
    bridge_report = bridge.build_apex_bridge(bridge_system, protocol)
    ghost_before_serialization = dg.inspect_dynamic_ghost_force(bridge_system)

    write_worker_log(log_path, "serialize bridge system")
    serialized = mm.XmlSerializer.serialize(bridge_system)
    bridge_xml_path = cell_bridge_xml_path(cell)
    atomic_write_text(bridge_xml_path, serialized)
    bridge_xml_hash = sha256_text(serialized)

    write_worker_log(log_path, "evaluate bridge xi states and inactive globals")
    bridge_values, inactive_values = evaluate_bridge(
        bridge_system,
        pdb.positions,
        states["apex_dplus"],
        xi_values,
        check_inactive_globals=True,
    )
    del bridge_system
    gc.collect()

    write_worker_log(log_path, "deserialize and evaluate endpoint/midpoint parity")
    roundtrip_system = mm.XmlSerializer.deserialize(serialized)
    bridge_roundtrip_contract = bridge.inspect_apex_bridge(roundtrip_system)
    ghost_after_serialization = dg.inspect_dynamic_ghost_force(roundtrip_system)
    roundtrip_values, _unused = evaluate_bridge(
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
        "xi0_vs_source_dplus": parity_record(
            bridge_values["0.00"]["energy_kj_mol"],
            bridge_values["0.00"]["forces"],
            plus["energy_kj_mol"],
            plus["forces"],
        ),
        "xi1_vs_source_dminus": parity_record(
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
        expected_energy = (1.0 - xi) * plus["energy_kj_mol"] + xi * minus["energy_kj_mol"]
        expected_forces = (1.0 - xi) * plus["forces"] + xi * minus["forces"]
        interior_parities[key] = parity_record(
            bridge_values[key]["energy_kj_mol"],
            bridge_values[key]["forces"],
            expected_energy,
            expected_forces,
        )
        derivative_parities[key] = derivative_record(
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
            abs(row["derivatives_kj_mol"][dg.GHOST_GLOBAL_PARAMETER] - dvdg_reference)
            for row in ghost_rows
        ),
    }
    ghost_invariance["passed"] = bool(
        ghost_invariance["max_abs_ghost_energy_delta_kj_mol"] <= ENERGY_TOLERANCE
        and ghost_invariance["max_abs_dVdg_delta_kj_mol"] <= DERIVATIVE_TOLERANCE
    )

    if inactive_values is None:
        raise RuntimeError("inactive-global audit was not evaluated")
    inactive_parity = parity_record(
        inactive_values["changed"]["energy_kj_mol"],
        inactive_values["changed"]["forces"],
        inactive_values["baseline"]["energy_kj_mol"],
        inactive_values["baseline"]["forces"],
    )
    roundtrip_parities = {
        key: parity_record(
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
        "max_abs_force_component_gap_kj_mol_nm": float(np.max(np.abs(gap_forces))),
        "is_pass_fail_threshold": False,
    }
    endpoint_gap["abs_gap_identity_delta_kj_mol"] = abs(
        endpoint_gap["bridge_xi1_minus_xi0_kj_mol"] - gap_energy
    )

    roundtrip_contract_matches = bridge_roundtrip_contract == bridge_report["after"]
    all_checks = {
        "bridge_mutation_scope": all(bridge_report["checks"].values()),
        "dynamic_ghost_contract_roundtrip": (
            ghost_before_serialization == ghost_after_serialization
        ),
        "bridge_contract_roundtrip": roundtrip_contract_matches,
        "serialized_xml_hash_matches_disk": sha256_file(bridge_xml_path) == bridge_xml_hash,
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
    passed = all(all_checks.values())
    result = {
        "schema": "updd_dynamic_ghost_apex_bridge_p0v2_cell_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "P0V2_PASS" if passed else "REJECT_BRIDGE",
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
        "source_cell_digest": frozen_cell["cell_digest"],
        "source_artifacts": frozen_cell["artifacts"],
        "implementation": inventory["implementation"],
        "bridge_system": {
            "path": relative(bridge_xml_path),
            "size": bridge_xml_path.stat().st_size,
            "sha256": bridge_xml_hash,
        },
        "bridge_contract": _compact_bridge_report(bridge_report),
        "ghost_contract": parent_p0.compact_force_report(ghost_report),
        "states": {
            "source_apex_dplus": states["apex_dplus"],
            "source_apex_dminus": states["apex_dminus"],
            "bridge_xi": xi_values,
            "ghost_g": 0.5,
        },
        "evaluations": {
            "source": {
                name: evaluation_record(value) for name, value in source_values.items()
            },
            "bridge": {
                name: evaluation_record(value) for name, value in bridge_values.items()
            },
            "roundtrip": {
                name: evaluation_record(value) for name, value in roundtrip_values.items()
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
        "checks": all_checks,
        "elapsed_s": time.monotonic() - started,
    }
    atomic_write_json(cell_result_path(cell), result)
    write_worker_log(log_path, f"COMPLETE {cell.cell_id} status={result['status']}")
    return result


def validate_completed_cell(
    cell: Cell, inventory: dict[str, Any]
) -> dict[str, Any] | None:
    path = cell_result_path(cell)
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
        raise ValueError(f"{cell.cell_id}: result identity mismatch")
    if result.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError(f"{cell.cell_id}: result inventory mismatch")
    if result.get("source_cell_digest") != frozen_cell_record(inventory, cell)["cell_digest"]:
        raise ValueError(f"{cell.cell_id}: result source digest mismatch")
    bridge_system = result.get("bridge_system") or {}
    system_path = REPO_ROOT / bridge_system.get("path", "")
    if not system_path.is_file() or sha256_file(system_path) != bridge_system.get("sha256"):
        raise ValueError(f"{cell.cell_id}: bridge XML missing or drifted")
    if result.get("status") not in {"P0V2_PASS", "REJECT_BRIDGE"}:
        raise ValueError(f"{cell.cell_id}: invalid completed status")
    return result


def worker_command(cell: Cell, inventory_digest_value: str) -> list[str]:
    return [
        sys.executable,
        str(RUNNER_PATH),
        "--worker-cell",
        cell.cell_id,
        "--expected-inventory-digest",
        inventory_digest_value,
    ]


def launch_worker(cell: Cell, inventory_digest_value: str) -> int:
    return subprocess.run(
        worker_command(cell, inventory_digest_value),
        cwd=REPO_ROOT,
        check=False,
    ).returncode


def write_run_state(status: str, **extra: Any) -> None:
    atomic_write_json(
        P0V2_ROOT / "run_state.json",
        {
            "status": status,
            "updated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            **extra,
        },
    )


def write_summary(
    results: list[dict[str, Any]], inventory: dict[str, Any]
) -> dict[str, Any]:
    cells = make_cells()
    by_id = {row["cell"]["cell_id"]: row for row in results}
    pending = [cell.cell_id for cell in cells if cell.cell_id not in by_id]
    rejected = [
        cell_id
        for cell_id, row in by_id.items()
        if row.get("status") == "REJECT_BRIDGE"
    ]
    elapsed = [float(row["elapsed_s"]) for row in results]
    if rejected:
        status = "REJECT_BRIDGE"
    elif not pending and len(results) == len(cells):
        status = "P0V2_PASS"
    else:
        status = "RUNNING"
    median_elapsed = median(elapsed) if elapsed else None
    payload = {
        "schema": "updd_dynamic_ghost_apex_bridge_p0v2_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "platform": "Reference",
        "gpu_used": False,
        "md_steps": 0,
        "minimization_steps": 0,
        "free_energy_estimated": False,
        "inventory_digest": inventory["inventory_digest"],
        "n_expected": len(cells),
        "n_completed": len(results),
        "pending_cell_ids": pending,
        "rejected_cell_ids": rejected,
        "median_elapsed_s": median_elapsed,
        "eta_s": median_elapsed * len(pending) if median_elapsed is not None else None,
        "cells": [
            {
                "cell_id": row["cell"]["cell_id"],
                "seed": row["cell"]["seed"],
                "leg": row["cell"]["leg"],
                "status": row["status"],
                "elapsed_s": row["elapsed_s"],
                "endpoint_gap_kj_mol": row["endpoint_gap"][
                    "source_hminus_minus_hplus_kj_mol"
                ],
                "max_force_gap_kj_mol_nm": row["endpoint_gap"][
                    "max_abs_force_component_gap_kj_mol_nm"
                ],
                "result": artifact_record(
                    cell_result_path(cell_by_id(row["cell"]["cell_id"]))
                ),
                "bridge_system": row["bridge_system"],
            }
            for row in results
        ],
    }
    atomic_write_json(P0V2_SUMMARY_JSON, payload)

    lines = [
        "# Apex Bridge P0v2 Reference Summary",
        "",
        f"Generated: {payload['generated']}",
        "",
        f"Status: **{status}**",
        "",
        "P0v2 used OpenMM Reference with zero MD/minimization steps and no GPU. "
        "It did not estimate a free energy.",
        "",
        f"Inventory digest: `{inventory['inventory_digest']}`",
        "",
        "| cell | status | endpoint gap (kJ/mol) | max force gap (kJ/mol/nm) | elapsed (s) |",
        "|---|---|---:|---:|---:|",
    ]
    for row in payload["cells"]:
        lines.append(
            f"| {row['cell_id']} | {row['status']} | "
            f"{row['endpoint_gap_kj_mol']:.12f} | "
            f"{row['max_force_gap_kj_mol_nm']:.12f} | {row['elapsed_s']:.3f} |"
        )
    if pending:
        lines.extend(["", "Pending: " + ", ".join(f"`{item}`" for item in pending)])
    if rejected:
        lines.extend(["", "Rejected: " + ", ".join(f"`{item}`" for item in rejected)])
    atomic_write_text(P0V2_SUMMARY_MD, "\n".join(lines) + "\n")
    return payload


def run_parent(inventory: dict[str, Any]) -> dict[str, Any]:
    validate_keeper_preflight(inventory)
    P0V2_ROOT.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for cell in make_cells():
        completed = validate_completed_cell(cell, inventory)
        if completed is not None:
            results.append(completed)
            write_summary(results, inventory)
            if completed["status"] == "REJECT_BRIDGE":
                write_run_state(
                    "REJECT_BRIDGE",
                    failed_cell=cell.cell_id,
                    n_completed=len(results),
                )
                return write_summary(results, inventory)
            continue

        if cell_output_dir(cell).exists():
            archive_partial_cell(cell)
        write_run_state("RUNNING", current_cell=cell.cell_id, n_completed=len(results))
        returncode = launch_worker(cell, inventory["inventory_digest"])
        if returncode != 0:
            write_run_state(
                "P0V2_INTERRUPTED",
                failed_cell=cell.cell_id,
                worker_returncode=returncode,
                n_completed=len(results),
            )
            write_summary(results, inventory)
            raise RuntimeError(
                f"P0v2 worker {cell.cell_id} exited with status {returncode}"
            )
        completed = validate_completed_cell(cell, inventory)
        if completed is None:
            raise RuntimeError(f"{cell.cell_id}: worker produced no valid result")
        results.append(completed)
        write_summary(results, inventory)
        if completed["status"] == "REJECT_BRIDGE":
            write_run_state(
                "REJECT_BRIDGE",
                failed_cell=cell.cell_id,
                n_completed=len(results),
            )
            return write_summary(results, inventory)

    write_run_state("P0V2_PASS", n_completed=len(results))
    return write_summary(results, inventory)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--prepare-inventory", action="store_true")
    group.add_argument("--keeper-preflight", action="store_true")
    group.add_argument("--run", action="store_true")
    group.add_argument("--summarize", action="store_true")
    group.add_argument("--worker-cell")
    parser.add_argument("--expected-inventory-digest")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.prepare_inventory:
            if P0V2_INVENTORY.exists():
                raise FileExistsError(
                    f"{P0V2_INVENTORY} already exists; frozen inventory is immutable"
                )
            inventory = build_source_inventory()
            atomic_write_json(P0V2_INVENTORY, inventory)
            print(inventory["inventory_digest"])
            return 0
        if args.keeper_preflight:
            inventory = verify_source_inventory()
            if P0V2_KEEPER_PREFLIGHT.exists():
                raise FileExistsError(
                    f"{P0V2_KEEPER_PREFLIGHT} already exists; preflight is immutable"
                )
            report = run_keeper_preflight(inventory)
            print(report["status"])
            return 0
        if args.worker_cell:
            if not args.expected_inventory_digest:
                raise ValueError("--worker-cell requires --expected-inventory-digest")
            inventory = load_frozen_inventory(args.expected_inventory_digest)
            validate_keeper_preflight(inventory)
            cell = cell_by_id(args.worker_cell)
            result = run_worker(cell, inventory)
            return 0 if result["status"] in {"P0V2_PASS", "REJECT_BRIDGE"} else 2
        if args.run:
            inventory = verify_source_inventory()
            summary = run_parent(inventory)
            print(summary["status"])
            return 0 if summary["status"] == "P0V2_PASS" else 2
        if args.summarize:
            inventory = verify_source_inventory()
            results = []
            for cell in make_cells():
                completed = validate_completed_cell(cell, inventory)
                if completed is not None:
                    results.append(completed)
            summary = write_summary(results, inventory)
            print(summary["status"])
            return 0
    except Exception:
        traceback.print_exc()
        return 1
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
