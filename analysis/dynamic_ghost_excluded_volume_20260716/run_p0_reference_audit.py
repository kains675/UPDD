#!/usr/bin/env python3
"""Run the frozen dynamic-ghost P0 Hamiltonian audit on OpenMM Reference."""

from __future__ import annotations

import argparse
import ast
import gc
import hashlib
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


REPO_ROOT = Path(__file__).resolve().parents[2]
UTILS_DIR = REPO_ROOT / "utils"
if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))

import trackb_dynamic_ghost as dg


ANALYSIS_DIR = Path(__file__).resolve().parent
PREREGISTRATION = ANALYSIS_DIR / "PREREGISTRATION.md"
PROTOCOL_PATH = ANALYSIS_DIR / "protocol.json"
FROZEN_MANIFEST = ANALYSIS_DIR / "FROZEN_MANIFEST.json"
SOURCE_INVENTORY_PATH = ANALYSIS_DIR.parent / "w4a_postdensify_dcd_20260716/source_inventory.json"

P0_INVENTORY = ANALYSIS_DIR / "P0_SOURCE_INVENTORY.json"
P0_KEEPER_PREFLIGHT = ANALYSIS_DIR / "P0_KEEPER_PREFLIGHT.json"
P0_ROOT = ANALYSIS_DIR / "p0_reference_audit"
P0_SUMMARY_JSON = ANALYSIS_DIR / "P0_REFERENCE_SUMMARY.json"
P0_SUMMARY_MD = ANALYSIS_DIR / "P0_REFERENCE_SUMMARY.md"

SOURCE_ROOT = REPO_ROOT / "outputs/_trackb/w4a_bdcd_h18_20260710_rerun1"
IMPLEMENTATION_PATH = REPO_ROOT / "utils/trackb_dynamic_ghost.py"
RUNNER_PATH = Path(__file__).resolve()

PILOT_SEEDS = ("s101", "s127", "s163")
ALL_SOURCE_SEEDS = (
    "s7",
    "s19",
    "s23",
    "s42",
    "s83",
    "s101",
    "s127",
    "s163",
    "s199",
    "s251",
)
SEED_TO_REPLICATE = {seed: index for index, seed in enumerate(ALL_SOURCE_SEEDS)}
LEGS = ("bound", "free")

ENERGY_TOLERANCE = 1.0e-5
FORCE_TOLERANCE = 1.0e-5
MINIMIZATION_MAX_ITERATIONS = 1
KCAL_TO_KJ = 4.184


@dataclass(frozen=True)
class Cell:
    cell_id: str
    seed: str
    leg: str
    replicate_index: int
    source_rep_dir: Path

    @property
    def xml_path(self) -> Path:
        return self.source_rep_dir / f"inplace_rbfe_{self.leg}_sys.xml"

    @property
    def pdb_path(self) -> Path:
        return self.source_rep_dir / f"inplace_rbfe_{self.leg}.pdb"

    @property
    def manifest_path(self) -> Path:
        return self.source_rep_dir / "run_manifest.json"

    def cntl_path(self, direction: str) -> Path:
        return self.source_rep_dir / direction / f"trackb_{direction}_asyncre.cntl"


def make_cells() -> tuple[Cell, ...]:
    cells: list[Cell] = []
    for leg in LEGS:
        source_lane = SOURCE_ROOT / f"uncarved_{leg}" / "wt" / leg
        for seed in PILOT_SEEDS:
            replicate_index = SEED_TO_REPLICATE[seed]
            cells.append(
                Cell(
                    cell_id=f"w4a_uncarved_{seed}_{leg}",
                    seed=seed,
                    leg=leg,
                    replicate_index=replicate_index,
                    source_rep_dir=source_lane / f"rep{replicate_index}",
                )
            )
    return tuple(cells)


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
    if frozen.get("status") != "FROZEN_BEFORE_IMPLEMENTATION_OR_OUTPUT":
        raise ValueError("frozen manifest status drifted")
    for name, expected_hash in (frozen.get("files") or {}).items():
        path = ANALYSIS_DIR / name
        observed = sha256_file(path)
        if observed != expected_hash:
            raise ValueError(
                f"frozen preregistration drift: {path} {observed} != {expected_hash}"
            )
    for path_text, expected_hash in (frozen.get("trigger_evidence") or {}).items():
        path = REPO_ROOT / path_text
        observed = sha256_file(path)
        if observed != expected_hash:
            raise ValueError(
                f"frozen trigger evidence drift: {path} {observed} != {expected_hash}"
            )

    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    dg.validate_ghost_protocol(protocol)
    p0 = ((protocol.get("stages") or {}).get("P0") or {})
    if p0 != {"platform": "Reference", "md_allowed": False, "legs": ["bound", "free"]}:
        raise ValueError(f"frozen P0 declaration drifted: {p0!r}")
    if protocol.get("system", {}).get("seeds") != list(PILOT_SEEDS):
        raise ValueError("frozen P0 seed cohort drifted")

    source_inventory = json.loads(SOURCE_INVENTORY_PATH.read_text(encoding="utf-8"))
    expected_digest = (protocol.get("evidence") or {}).get("source_inventory_digest")
    if source_inventory.get("inventory_digest") != expected_digest:
        raise ValueError(
            "trigger source inventory digest mismatch: "
            f"{source_inventory.get('inventory_digest')} != {expected_digest}"
        )
    return protocol, frozen


def parse_cntl_value(path: Path, key: str) -> Any:
    prefix = key + " ="
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith(prefix):
            return ast.literal_eval(line.split("=", 1)[1].strip())
    raise ValueError(f"{path}: missing {key}")


def parse_vector(path: Path, key: str, cast=float) -> tuple[Any, ...]:
    raw = parse_cntl_value(path, key)
    if not isinstance(raw, str):
        raise ValueError(f"{path}: {key} must be a quoted comma-separated vector")
    return tuple(cast(item.strip()) for item in raw.split(","))


def parse_schedule(path: Path) -> dict[str, Any]:
    schedule = {
        "lambda1": parse_vector(path, "LAMBDA1", float),
        "lambda2": parse_vector(path, "LAMBDA2", float),
        "directions": parse_vector(path, "DIRECTION", int),
        "intermediate": parse_vector(path, "INTERMEDIATE", int),
        "alpha": parse_vector(path, "ALPHA", float),
        "u0": parse_vector(path, "U0", float),
        "w0": parse_vector(path, "W0COEFF", float),
        "umax": float(parse_cntl_value(path, "UMAX")),
        "ubcore": float(parse_cntl_value(path, "UBCORE")),
        "acore": float(parse_cntl_value(path, "ACORE")),
    }
    lengths = {
        len(value)
        for key, value in schedule.items()
        if key in {"lambda1", "lambda2", "directions", "intermediate", "alpha", "u0", "w0"}
    }
    if lengths != {15}:
        raise ValueError(f"{path}: H18 schedule vector lengths are {sorted(lengths)}")
    return schedule


def validate_h18_schedules(dplus: dict[str, Any], dminus: dict[str, Any], protocol: dict[str, Any]) -> None:
    h18 = protocol.get("h18") or {}
    ramp1 = tuple(float(value) for value in h18.get("lambda1_rampdown") or [])
    ramp2 = tuple(float(value) for value in h18.get("lambda2_rampup") or [])
    expected = {
        "dplus_lambda1": (0.0,) * 8 + ramp1,
        "dplus_lambda2": ramp2 + (0.5,) * 7,
        "dminus_lambda1": tuple(reversed(ramp1)) + (0.0,) * 8,
        "dminus_lambda2": (0.5,) * 8 + tuple(reversed(ramp2[:-1])),
    }
    observed = {
        "dplus_lambda1": dplus["lambda1"],
        "dplus_lambda2": dplus["lambda2"],
        "dminus_lambda1": dminus["lambda1"],
        "dminus_lambda2": dminus["lambda2"],
    }
    if observed != expected:
        raise ValueError("source H18 lambda schedule differs from frozen protocol")
    if dplus["directions"] != (1,) * 15 or dminus["directions"] != (-1,) * 15:
        raise ValueError("source H18 direction vectors drifted")
    for schedule in (dplus, dminus):
        if schedule["alpha"] != (0.1,) * 15:
            raise ValueError("source H18 ALPHA drifted")
        if schedule["u0"] != (110.0,) * 15 or schedule["w0"] != (0.0,) * 15:
            raise ValueError("source H18 U0/W0 drifted")
        if (schedule["umax"], schedule["ubcore"], schedule["acore"]) != (200.0, 100.0, 0.0625):
            raise ValueError("source H18 soft-core constants drifted")


def state_at(schedule: dict[str, Any], index: int, label: str) -> dict[str, Any]:
    return {
        "label": label,
        "lambda1": float(schedule["lambda1"][index]),
        "lambda2": float(schedule["lambda2"][index]),
        "direction": int(schedule["directions"][index]),
        "intermediate": int(schedule["intermediate"][index]),
        "alpha_kcal_inv": float(schedule["alpha"][index]),
        "uh_kcal_mol": float(schedule["u0"][index]),
        "w0_kcal_mol": float(schedule["w0"][index]),
        "umax_kcal_mol": float(schedule["umax"]),
        "ubcore_kcal_mol": float(schedule["ubcore"]),
        "acore": float(schedule["acore"]),
    }


def validate_leg_preregistration(leg: str) -> dict[str, Any]:
    path = SOURCE_ROOT / f"uncarved_{leg}" / "pre_registration.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    config = payload.get("config") or {}
    expected = {
        "leg": leg,
        "endpoints": ["wt"],
        "seeds": list(ALL_SOURCE_SEEDS),
        "directions": ["dplus", "dminus"],
        "n_cycles": 800,
        "md_steps_per_cycle": 250,
        "timestep_fs": 1.0,
        "construction": "twocopy",
        "mutation_spec": "w4a_trp_ala_res4",
        "carve_void_waters": False,
        "mintimeid": 600,
    }
    for key, value in expected.items():
        if config.get(key) != value:
            raise ValueError(
                f"{path}: {key}={config.get(key)!r}, expected {value!r}"
            )
    return {"declaration": expected, "artifact": artifact_record(path)}


def validate_source_cell(cell: Cell, protocol: dict[str, Any]) -> dict[str, Any]:
    manifest = json.loads(cell.manifest_path.read_text(encoding="utf-8"))
    expected = {
        "seed": cell.seed,
        "replicate_index": cell.replicate_index,
        "endpoint": "wt",
        "leg": cell.leg,
        "construction": "twocopy",
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(
                f"{cell.manifest_path}: {key}={manifest.get(key)!r}, expected {value!r}"
            )
    serialize = manifest.get("serialize") or {}
    if Path(serialize.get("sys_xml_path", "")).resolve() != cell.xml_path.resolve():
        raise ValueError(f"{cell.cell_id}: serialized XML path mismatch")
    if Path(serialize.get("pdb_path", "")).resolve() != cell.pdb_path.resolve():
        raise ValueError(f"{cell.cell_id}: serialized PDB path mismatch")

    dplus_path = cell.cntl_path("dplus")
    dminus_path = cell.cntl_path("dminus")
    dplus = parse_schedule(dplus_path)
    dminus = parse_schedule(dminus_path)
    validate_h18_schedules(dplus, dminus, protocol)
    artifacts = {
        "serialized_xml": artifact_record(cell.xml_path),
        "serialized_pdb": artifact_record(cell.pdb_path),
        "run_manifest": artifact_record(cell.manifest_path),
        "dplus_cntl": artifact_record(dplus_path),
        "dminus_cntl": artifact_record(dminus_path),
    }
    reference_source = ((protocol.get("ghost") or {}).get("parameter_source") or {})
    if cell.seed == "s127" and cell.leg == "bound":
        if artifacts["serialized_xml"]["sha256"] != reference_source.get("serialized_xml_sha256"):
            raise ValueError("s127 bound parameter-source XML hash drifted")
        if artifacts["serialized_pdb"]["sha256"] != reference_source.get("serialized_pdb_sha256"):
            raise ValueError("s127 bound parameter-source PDB hash drifted")
        if artifacts["run_manifest"]["sha256"] != reference_source.get("run_manifest_sha256"):
            raise ValueError("s127 bound parameter-source manifest hash drifted")

    states = {
        "u0_dplus": state_at(dplus, 0, "u0_dplus"),
        "apex_dplus": state_at(dplus, -1, "apex_dplus"),
        "apex_dminus": state_at(dminus, 0, "apex_dminus"),
        "u1_dminus": state_at(dminus, -1, "u1_dminus"),
    }
    expected_g = {
        "u0_dplus": 0.0,
        "apex_dplus": 0.5,
        "apex_dminus": 0.5,
        "u1_dminus": 1.0,
    }
    observed_g = {
        "u0_dplus": states["u0_dplus"]["lambda2"],
        "apex_dplus": states["apex_dplus"]["lambda2"],
        "apex_dminus": 1.0 - states["apex_dminus"]["lambda2"],
        "u1_dminus": 1.0 - states["u1_dminus"]["lambda2"],
    }
    if observed_g != expected_g:
        raise ValueError(f"{cell.cell_id}: main-path ghost coupling drifted: {observed_g}")

    record = {
        "cell_id": cell.cell_id,
        "seed": cell.seed,
        "leg": cell.leg,
        "replicate_index": cell.replicate_index,
        "source_rep_dir": relative(cell.source_rep_dir),
        "manifest_declaration": expected,
        "artifacts": artifacts,
        "states": states,
        "main_path_g": observed_g,
    }
    record["cell_digest"] = canonical_digest(record)
    return record


def inventory_digest(payload: dict[str, Any]) -> str:
    canonical = dict(payload)
    canonical.pop("generated", None)
    canonical.pop("inventory_digest", None)
    return canonical_digest(canonical)


def build_source_inventory() -> dict[str, Any]:
    protocol, frozen = load_protocol_and_verify_freeze()
    payload = {
        "schema": "updd_dynamic_ghost_p0_source_inventory_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "platform": "Reference",
        "md_allowed": False,
        "gpu_allowed": False,
        "minimization_max_iterations": MINIMIZATION_MAX_ITERATIONS,
        "frozen_manifest": artifact_record(FROZEN_MANIFEST),
        "frozen_files": {
            name: artifact_record(ANALYSIS_DIR / name)
            for name in sorted((frozen.get("files") or {}))
        },
        "trigger_evidence": {
            path_text: artifact_record(REPO_ROOT / path_text)
            for path_text in sorted((frozen.get("trigger_evidence") or {}))
        },
        "trigger_source_inventory": artifact_record(SOURCE_INVENTORY_PATH),
        "implementation": {
            "force_module": artifact_record(IMPLEMENTATION_PATH),
            "p0_runner": artifact_record(RUNNER_PATH),
        },
        "source_preregistrations": {
            leg: validate_leg_preregistration(leg) for leg in LEGS
        },
        "cells": [validate_source_cell(cell, protocol) for cell in make_cells()],
    }
    payload["inventory_digest"] = inventory_digest(payload)
    return payload


def load_frozen_inventory(expected_digest: str | None = None) -> dict[str, Any]:
    if not P0_INVENTORY.is_file():
        raise FileNotFoundError(
            f"missing {P0_INVENTORY}; run --prepare-inventory before P0"
        )
    inventory = json.loads(P0_INVENTORY.read_text(encoding="utf-8"))
    if inventory.get("inventory_digest") != inventory_digest(inventory):
        raise ValueError("stored P0 source inventory digest is invalid")
    if expected_digest is not None and inventory["inventory_digest"] != expected_digest:
        raise ValueError(
            f"worker inventory mismatch: {inventory['inventory_digest']} != {expected_digest}"
        )
    return inventory


def verify_source_inventory() -> dict[str, Any]:
    frozen = load_frozen_inventory()
    current = build_source_inventory()
    if frozen["inventory_digest"] != current["inventory_digest"]:
        raise ValueError(
            "P0 source/code drift after inventory freeze: "
            f"{frozen['inventory_digest']} != {current['inventory_digest']}"
        )
    return frozen


def cell_by_id(cell_id: str) -> Cell:
    matches = [cell for cell in make_cells() if cell.cell_id == cell_id]
    if len(matches) != 1:
        raise ValueError(f"expected one P0 cell {cell_id!r}, found {len(matches)}")
    return matches[0]


def frozen_cell_record(inventory: dict[str, Any], cell: Cell) -> dict[str, Any]:
    matches = [row for row in inventory.get("cells", []) if row.get("cell_id") == cell.cell_id]
    if len(matches) != 1:
        raise ValueError(f"{cell.cell_id}: missing or duplicate frozen cell inventory")
    return matches[0]


def verify_worker_inputs(inventory: dict[str, Any], cell: Cell) -> tuple[dict[str, Any], dict[str, Any]]:
    protocol, _frozen = load_protocol_and_verify_freeze()
    expected_code = inventory.get("implementation") or {}
    current_code = {
        "force_module": artifact_record(IMPLEMENTATION_PATH),
        "p0_runner": artifact_record(RUNNER_PATH),
    }
    if current_code != expected_code:
        raise ValueError(f"{cell.cell_id}: implementation code drift after inventory freeze")
    current_cell = validate_source_cell(cell, protocol)
    frozen_cell = frozen_cell_record(inventory, cell)
    if current_cell != frozen_cell:
        raise ValueError(f"{cell.cell_id}: source cell drift after inventory freeze")
    return protocol, frozen_cell


def run_keeper_preflight(inventory: dict[str, Any]) -> dict[str, Any]:
    """Validate all six force declarations before any OpenMM Context exists."""
    import openmm as mm
    from openmm.app import PDBFile

    rows: list[dict[str, Any]] = []
    for cell in make_cells():
        protocol, frozen_cell = verify_worker_inputs(inventory, cell)
        system = mm.XmlSerializer.deserialize(cell.xml_path.read_text(encoding="utf-8"))
        pdb = PDBFile(str(cell.pdb_path))
        report = dg.build_dynamic_ghost_force(system, pdb.topology, protocol)
        compact = compact_force_report(report)
        checks = {
            "exact_nine_ring_atoms": compact["selection"]["n_ring_atoms"] == 9,
            "all_water_oxygens": (
                compact["selection"]["n_water_oxygens"]
                == compact["selection"]["n_water_residues"]
                and compact["selection"]["n_water_oxygens"] > 0
            ),
            "atom_count_invariant": compact["n_particles_before"] == compact["n_particles_after"],
            "charge_invariant": compact["net_charge_delta_e"] == 0.0,
            "top_level_after_atm": compact["force_index"] > compact["atm_force_index"],
            "dedicated_force_group": compact["force_group"] == dg.GHOST_FORCE_GROUP,
            "source_pair_parameters_match": True,
        }
        if not all(checks.values()):
            raise dg.DynamicGhostError(
                f"{cell.cell_id}: KEEPER declaration checks failed: {checks}"
            )
        rows.append(
            {
                "cell_id": cell.cell_id,
                "seed": cell.seed,
                "leg": cell.leg,
                "source_cell_digest": frozen_cell["cell_digest"],
                "checks": checks,
                "force_contract": compact,
            }
        )
        del system, pdb, report
        gc.collect()

    payload = {
        "schema": "updd_dynamic_ghost_p0_keeper_preflight_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "PASS",
        "context_count": 0,
        "md_steps": 0,
        "gpu_used": False,
        "inventory_digest": inventory["inventory_digest"],
        "n_cells": len(rows),
        "cells": rows,
    }
    atomic_write_json(P0_KEEPER_PREFLIGHT, payload)
    return payload


def validate_keeper_preflight(inventory: dict[str, Any]) -> dict[str, Any]:
    if not P0_KEEPER_PREFLIGHT.is_file():
        raise FileNotFoundError(
            f"missing {P0_KEEPER_PREFLIGHT}; run --keeper-preflight before P0"
        )
    payload = json.loads(P0_KEEPER_PREFLIGHT.read_text(encoding="utf-8"))
    if payload.get("status") != "PASS" or payload.get("context_count") != 0:
        raise ValueError("P0 KEEPER preflight is not a context-free PASS")
    if payload.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError("P0 KEEPER preflight inventory digest drifted")
    if payload.get("n_cells") != len(make_cells()):
        raise ValueError("P0 KEEPER preflight cell count drifted")
    return payload


def cell_output_dir(cell: Cell) -> Path:
    return P0_ROOT / cell.cell_id


def cell_result_path(cell: Cell) -> Path:
    return cell_output_dir(cell) / "cell_result.json"


def cell_ghost_xml_path(cell: Cell) -> Path:
    return cell_output_dir(cell) / "dynamic_ghost_system.xml"


def archive_partial_cell(cell: Cell) -> Path | None:
    root = cell_output_dir(cell)
    if not root.exists():
        return None
    stamp = time.strftime("%Y%m%d_%H%M%S")
    target = P0_ROOT / "_archive" / f"{stamp}_{cell.cell_id}"
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(root), str(target))
    atomic_write_json(
        target / "archive_reason.json",
        {
            "cell_id": cell.cell_id,
            "archived": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            "reason": "incomplete or invalid P0 cell output",
        },
    )
    return target


def compact_force_report(report: dict[str, Any]) -> dict[str, Any]:
    compact = json.loads(json.dumps(report))
    selection = compact.get("selection") or {}
    water_indices = selection.pop("water_oxygen_indices_sha_input", [])
    selection["water_oxygen_indices_sha256"] = canonical_digest(water_indices)
    selection["water_oxygen_index_first"] = water_indices[0] if water_indices else None
    selection["water_oxygen_index_last"] = water_indices[-1] if water_indices else None
    for group in compact.get("interaction_groups") or []:
        second = group.pop("second", [])
        group["second_count"] = len(second)
        group["second_sha256"] = canonical_digest(second)
    return compact


def _set_context_state(context: Any, state: dict[str, Any], ghost_g: float | None) -> None:
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


def evaluate_reference(
    system: Any,
    positions: Any,
    atm_state: dict[str, Any],
    *,
    ghost_g: float | None,
    minimize_iterations: int = 0,
) -> dict[str, Any]:
    import numpy as np
    import openmm as mm
    from openmm import unit

    started = time.monotonic()
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(
        system,
        integrator,
        mm.Platform.getPlatformByName("Reference"),
    )
    context.setPositions(positions)
    _set_context_state(context, atm_state, ghost_g)
    if minimize_iterations:
        mm.LocalEnergyMinimizer.minimize(
            context,
            tolerance=0.0,
            maxIterations=int(minimize_iterations),
        )
    state = context.getState(
        getEnergy=True,
        getForces=True,
        getPositions=bool(minimize_iterations),
        getParameterDerivatives=ghost_g is not None,
    )
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    forces = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=float,
    )
    output_positions = None
    if minimize_iterations:
        output_positions = np.asarray(
            state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
            dtype=float,
        )
    derivatives = dict(state.getEnergyParameterDerivatives()) if ghost_g is not None else {}
    ghost_energy = None
    if ghost_g is not None:
        ghost_energy = float(
            context.getState(
                getEnergy=True,
                groups=1 << dg.GHOST_FORCE_GROUP,
            ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        )
    elapsed = time.monotonic() - started
    result = {
        "energy_kj_mol": energy,
        "forces": forces,
        "positions_nm": output_positions,
        "ghost_energy_kj_mol": ghost_energy,
        "dVdg_kj_mol": (
            float(derivatives.get(dg.GHOST_GLOBAL_PARAMETER))
            if ghost_g is not None
            else None
        ),
        "elapsed_s": elapsed,
        "finite": bool(
            math.isfinite(energy)
            and np.isfinite(forces).all()
            and (ghost_energy is None or math.isfinite(ghost_energy))
            and (
                ghost_g is None
                or (
                    dg.GHOST_GLOBAL_PARAMETER in derivatives
                    and math.isfinite(float(derivatives[dg.GHOST_GLOBAL_PARAMETER]))
                )
            )
        ),
    }
    del context, integrator, state
    gc.collect()
    return result


def evaluation_record(value: dict[str, Any]) -> dict[str, Any]:
    import numpy as np

    forces = value["forces"]
    return {
        "energy_kj_mol": value["energy_kj_mol"],
        "max_abs_force_component_kj_mol_nm": float(np.max(np.abs(forces))),
        "ghost_energy_kj_mol": value["ghost_energy_kj_mol"],
        "dVdg_kj_mol": value["dVdg_kj_mol"],
        "elapsed_s": value["elapsed_s"],
        "finite": value["finite"],
    }


def parity_record(first: dict[str, Any], second: dict[str, Any]) -> dict[str, Any]:
    import numpy as np

    energy_delta = abs(first["energy_kj_mol"] - second["energy_kj_mol"])
    force_delta = float(np.max(np.abs(first["forces"] - second["forces"])))
    result = {
        "abs_energy_delta_kj_mol": energy_delta,
        "max_abs_force_component_delta_kj_mol_nm": force_delta,
        "energy_tolerance_kj_mol": ENERGY_TOLERANCE,
        "force_tolerance_kj_mol_nm": FORCE_TOLERANCE,
    }
    result["passed"] = bool(
        math.isfinite(energy_delta)
        and math.isfinite(force_delta)
        and energy_delta <= ENERGY_TOLERANCE
        and force_delta <= FORCE_TOLERANCE
    )
    return result


def minimized_parity_record(first: dict[str, Any], second: dict[str, Any]) -> dict[str, Any]:
    import numpy as np

    result = parity_record(first, second)
    result["max_abs_position_component_delta_nm"] = float(
        np.max(np.abs(first["positions_nm"] - second["positions_nm"]))
    )
    result["max_iterations"] = MINIMIZATION_MAX_ITERATIONS
    return result


def group_identity_record(ghost: dict[str, Any], physical: dict[str, Any]) -> dict[str, Any]:
    delta = (
        ghost["energy_kj_mol"]
        - physical["energy_kj_mol"]
        - float(ghost["ghost_energy_kj_mol"])
    )
    return {
        "total_minus_physical_minus_group_kj_mol": delta,
        "tolerance_kj_mol": ENERGY_TOLERANCE,
        "passed": bool(math.isfinite(delta) and abs(delta) <= ENERGY_TOLERANCE),
    }


def write_worker_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S %z')} {message}"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def run_worker(cell: Cell, inventory: dict[str, Any]) -> dict[str, Any]:
    import openmm as mm
    from openmm.app import PDBFile

    output_dir = cell_output_dir(cell)
    if output_dir.exists():
        raise RuntimeError(
            f"{cell.cell_id}: output exists; parent must validate or archive before launch"
        )
    output_dir.mkdir(parents=True)
    log_path = output_dir / "worker.log"
    write_worker_log(log_path, f"START {cell.cell_id} platform=Reference md_allowed=false")
    started = time.monotonic()

    protocol, frozen_cell = verify_worker_inputs(inventory, cell)
    source_xml = cell.xml_path.read_text(encoding="utf-8")
    pdb = PDBFile(str(cell.pdb_path))
    states = frozen_cell["states"]

    write_worker_log(log_path, "revalidate force scope and source parameters before Context creation")
    ghost_system = mm.XmlSerializer.deserialize(source_xml)
    force_report = dg.build_dynamic_ghost_force(ghost_system, pdb.topology, protocol)
    inspected_before = dg.inspect_dynamic_ghost_force(ghost_system)

    write_worker_log(log_path, "evaluate original u0/u1 and one-step u0 minimization")
    original_system = mm.XmlSerializer.deserialize(source_xml)
    original_u0 = evaluate_reference(
        original_system, pdb.positions, states["u0_dplus"], ghost_g=None
    )
    original_u1 = evaluate_reference(
        original_system, pdb.positions, states["u1_dminus"], ghost_g=None
    )
    original_min = evaluate_reference(
        original_system,
        pdb.positions,
        states["u0_dplus"],
        ghost_g=None,
        minimize_iterations=MINIMIZATION_MAX_ITERATIONS,
    )
    del original_system
    gc.collect()

    write_worker_log(log_path, "serialize dynamic ghost and verify round-trip")
    ghost_xml = mm.XmlSerializer.serialize(ghost_system)
    ghost_xml_path = cell_ghost_xml_path(cell)
    atomic_write_text(ghost_xml_path, ghost_xml)
    ghost_xml_hash = sha256_text(ghost_xml)

    write_worker_log(log_path, "evaluate ghost u0/u1/apex states and one-step u0 minimization")
    ghost_u0 = evaluate_reference(
        ghost_system, pdb.positions, states["u0_dplus"], ghost_g=0.0
    )
    correction_u1_g0 = evaluate_reference(
        ghost_system, pdb.positions, states["u1_dminus"], ghost_g=0.0
    )
    main_u1_g1 = evaluate_reference(
        ghost_system, pdb.positions, states["u1_dminus"], ghost_g=1.0
    )
    apex_dplus = evaluate_reference(
        ghost_system, pdb.positions, states["apex_dplus"], ghost_g=0.5
    )
    apex_dminus = evaluate_reference(
        ghost_system, pdb.positions, states["apex_dminus"], ghost_g=0.5
    )
    ghost_min = evaluate_reference(
        ghost_system,
        pdb.positions,
        states["u0_dplus"],
        ghost_g=0.0,
        minimize_iterations=MINIMIZATION_MAX_ITERATIONS,
    )
    del ghost_system
    gc.collect()

    write_worker_log(log_path, "evaluate serialized correction u1,g=1 parity")
    roundtrip_system = mm.XmlSerializer.deserialize(ghost_xml)
    inspected_after = dg.inspect_dynamic_ghost_force(roundtrip_system)
    correction_u1_g1 = evaluate_reference(
        roundtrip_system, pdb.positions, states["u1_dminus"], ghost_g=1.0
    )
    del roundtrip_system, ghost_xml, source_xml
    gc.collect()

    parities = {
        "u0_g0_vs_unmodified": parity_record(ghost_u0, original_u0),
        "correction_u1_g0_vs_unmodified": parity_record(correction_u1_g0, original_u1),
        "main_u1_g1_vs_serialized_correction_u1_g1": parity_record(
            main_u1_g1, correction_u1_g1
        ),
        "shared_apex_dplus_vs_dminus": parity_record(apex_dplus, apex_dminus),
        "minimized_u0_g0_vs_unmodified": minimized_parity_record(
            ghost_min, original_min
        ),
    }
    group_identity = {
        "u0_g0": group_identity_record(ghost_u0, original_u0),
        "u1_g0": group_identity_record(correction_u1_g0, original_u1),
        "u1_g1": group_identity_record(main_u1_g1, original_u1),
    }
    evaluations = {
        "original_u0": evaluation_record(original_u0),
        "original_u1": evaluation_record(original_u1),
        "original_minimized_u0": evaluation_record(original_min),
        "ghost_u0_g0": evaluation_record(ghost_u0),
        "correction_u1_g0": evaluation_record(correction_u1_g0),
        "main_u1_g1": evaluation_record(main_u1_g1),
        "apex_dplus_g0p5": evaluation_record(apex_dplus),
        "apex_dminus_g0p5": evaluation_record(apex_dminus),
        "ghost_minimized_u0_g0": evaluation_record(ghost_min),
        "serialized_correction_u1_g1": evaluation_record(correction_u1_g1),
    }
    finite_all = all(row["finite"] for row in evaluations.values())
    derivatives_finite = all(
        row["dVdg_kj_mol"] is not None and math.isfinite(row["dVdg_kj_mol"])
        for name, row in evaluations.items()
        if name.startswith(("ghost_", "correction_", "main_", "apex_", "serialized_"))
    )
    force_contract_roundtrip = inspected_before == inspected_after
    all_checks = {
        "force_scope_and_source_parameters": True,
        "atom_count_and_charge_invariant": bool(
            force_report["n_particles_before"] == force_report["n_particles_after"]
            and force_report["net_charge_delta_e"] == 0.0
        ),
        "serialization_contract_roundtrip": force_contract_roundtrip,
        "serialized_xml_hash_matches_disk": sha256_file(ghost_xml_path) == ghost_xml_hash,
        "finite_energy_forces_and_dVdg": bool(finite_all and derivatives_finite),
        "all_hamiltonian_parities": all(row["passed"] for row in parities.values()),
        "ghost_group_identity": all(row["passed"] for row in group_identity.values()),
    }
    passed = all(all_checks.values())
    result = {
        "schema": "updd_dynamic_ghost_p0_cell_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": "P0_PASS" if passed else "REJECT_GHOST",
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
            "gpu_used": False,
            "minimization_max_iterations": MINIMIZATION_MAX_ITERATIONS,
        },
        "inventory_digest": inventory["inventory_digest"],
        "source_cell_digest": frozen_cell["cell_digest"],
        "source_artifacts": frozen_cell["artifacts"],
        "implementation": inventory["implementation"],
        "ghost_system": {
            "path": relative(ghost_xml_path),
            "size": ghost_xml_path.stat().st_size,
            "sha256": ghost_xml_hash,
        },
        "force_contract": compact_force_report(force_report),
        "states": states,
        "evaluations": evaluations,
        "parities": parities,
        "ghost_group_identity": group_identity,
        "checks": all_checks,
        "elapsed_s": time.monotonic() - started,
    }
    atomic_write_json(cell_result_path(cell), result)
    write_worker_log(log_path, f"COMPLETE {cell.cell_id} status={result['status']}")
    return result


def validate_completed_cell(cell: Cell, inventory: dict[str, Any]) -> dict[str, Any] | None:
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
        raise ValueError(f"{cell.cell_id}: completed result identity mismatch")
    if result.get("inventory_digest") != inventory["inventory_digest"]:
        raise ValueError(f"{cell.cell_id}: completed result inventory mismatch")
    if result.get("source_cell_digest") != frozen_cell_record(inventory, cell)["cell_digest"]:
        raise ValueError(f"{cell.cell_id}: completed result source digest mismatch")
    ghost = result.get("ghost_system") or {}
    ghost_path = REPO_ROOT / ghost.get("path", "")
    if not ghost_path.is_file() or sha256_file(ghost_path) != ghost.get("sha256"):
        raise ValueError(f"{cell.cell_id}: ghost XML artifact missing or drifted")
    if result.get("status") not in {"P0_PASS", "REJECT_GHOST"}:
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
        P0_ROOT / "run_state.json",
        {
            "status": status,
            "updated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
            **extra,
        },
    )


def write_summary(results: list[dict[str, Any]], inventory: dict[str, Any]) -> dict[str, Any]:
    cells = make_cells()
    by_id = {row["cell"]["cell_id"]: row for row in results}
    pending = [cell.cell_id for cell in cells if cell.cell_id not in by_id]
    rejected = [
        cell_id for cell_id, row in by_id.items() if row.get("status") == "REJECT_GHOST"
    ]
    elapsed = [float(row["elapsed_s"]) for row in results]
    if rejected:
        status = "REJECT_GHOST"
    elif not pending and len(results) == len(cells):
        status = "P0_PASS"
    else:
        status = "RUNNING"
    payload = {
        "schema": "updd_dynamic_ghost_p0_summary_v1",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "status": status,
        "platform": "Reference",
        "md_steps": 0,
        "gpu_used": False,
        "inventory_digest": inventory["inventory_digest"],
        "n_expected": len(cells),
        "n_completed": len(results),
        "pending_cell_ids": pending,
        "rejected_cell_ids": rejected,
        "median_elapsed_s": median(elapsed) if elapsed else None,
        "eta_s": median(elapsed) * len(pending) if elapsed and not rejected else None,
        "cells": [
            {
                "cell_id": row["cell"]["cell_id"],
                "seed": row["cell"]["seed"],
                "leg": row["cell"]["leg"],
                "status": row["status"],
                "elapsed_s": row["elapsed_s"],
                "max_energy_delta_kj_mol": max(
                    value["abs_energy_delta_kj_mol"]
                    for value in row["parities"].values()
                ),
                "max_force_delta_kj_mol_nm": max(
                    value["max_abs_force_component_delta_kj_mol_nm"]
                    for value in row["parities"].values()
                ),
                "result": artifact_record(cell_result_path(cell_by_id(row["cell"]["cell_id"]))),
                "ghost_system": row["ghost_system"],
            }
            for row in results
        ],
    }
    atomic_write_json(P0_SUMMARY_JSON, payload)
    lines = [
        "# Dynamic Ghost P0 Reference Audit",
        "",
        f"Status: `{status}`",
        f"Completed: `{len(results)}/{len(cells)}`",
        "Platform: `Reference`",
        "MD steps: `0`",
        "GPU used: `false`",
        "",
        "| leg | seed | status | max dE (kJ/mol) | max dF component (kJ/mol/nm) | elapsed min |",
        "| --- | --- | --- | ---: | ---: | ---: |",
    ]
    for row in payload["cells"]:
        lines.append(
            "| {leg} | {seed} | {status} | {de:.6g} | {df:.6g} | {minutes:.2f} |".format(
                leg=row["leg"],
                seed=row["seed"],
                status=row["status"],
                de=row["max_energy_delta_kj_mol"],
                df=row["max_force_delta_kj_mol_nm"],
                minutes=row["elapsed_s"] / 60.0,
            )
        )
    if pending:
        lines.extend(["", "Pending: " + ", ".join(f"`{value}`" for value in pending)])
    if rejected:
        lines.extend(["", "Rejected: " + ", ".join(f"`{value}`" for value in rejected)])
    atomic_write_text(P0_SUMMARY_MD, "\n".join(lines) + "\n")
    return payload


def run_all(inventory: dict[str, Any]) -> int:
    validate_keeper_preflight(inventory)
    results: list[dict[str, Any]] = []
    for cell in make_cells():
        completed = validate_completed_cell(cell, inventory)
        if completed is not None:
            results.append(completed)
    initial_summary = write_summary(results, inventory)
    if initial_summary["status"] == "REJECT_GHOST":
        write_run_state(
            "REJECT_GHOST",
            n_completed=len(results),
            n_expected=len(make_cells()),
            reason="existing completed cell carries REJECT_GHOST",
        )
        return 2

    for cell in make_cells():
        if any(row["cell"]["cell_id"] == cell.cell_id for row in results):
            continue
        archived = archive_partial_cell(cell)
        if archived is not None:
            print(f"[P0] archived partial cell at {archived}", flush=True)
        write_run_state(
            "RUNNING",
            current_cell=cell.cell_id,
            n_completed=len(results),
            n_expected=len(make_cells()),
        )
        print(f"[P0] launch {cell.cell_id} on Reference", flush=True)
        return_code = launch_worker(cell, inventory["inventory_digest"])
        completed = validate_completed_cell(cell, inventory)
        if completed is None:
            write_run_state(
                "FAILED",
                current_cell=cell.cell_id,
                return_code=return_code,
                error="worker produced no valid completed result",
            )
            raise RuntimeError(
                f"{cell.cell_id}: worker exit {return_code} without a valid result"
            )
        results.append(completed)
        summary = write_summary(results, inventory)
        if completed["status"] == "REJECT_GHOST" or return_code == 2:
            write_run_state(
                "REJECT_GHOST",
                current_cell=cell.cell_id,
                n_completed=len(results),
                n_expected=len(make_cells()),
            )
            return 2
        if return_code != 0:
            write_run_state(
                "FAILED",
                current_cell=cell.cell_id,
                return_code=return_code,
            )
            raise RuntimeError(f"{cell.cell_id}: worker exited with {return_code}")
        print(
            f"[P0] complete {cell.cell_id}; {summary['n_completed']}/{summary['n_expected']}",
            flush=True,
        )

    write_run_state(
        "P0_PASS",
        n_completed=len(results),
        n_expected=len(make_cells()),
    )
    return 0


def run_worker_entry(cell_id: str, expected_digest: str) -> int:
    cell = cell_by_id(cell_id)
    inventory = load_frozen_inventory(expected_digest)
    try:
        result = run_worker(cell, inventory)
    except Exception as exc:
        output_dir = cell_output_dir(cell)
        output_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_json(
            output_dir / "failure.json",
            {
                "cell_id": cell.cell_id,
                "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
                "status": "REJECT_GHOST" if isinstance(exc, dg.DynamicGhostError) else "FAILED",
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            },
        )
        raise
    return 0 if result["status"] == "P0_PASS" else 2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-inventory", action="store_true")
    parser.add_argument("--keeper-preflight", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--worker-cell", help=argparse.SUPPRESS)
    parser.add_argument("--expected-inventory-digest", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    public_modes = sum(
        bool(value)
        for value in (
            args.prepare_inventory,
            args.keeper_preflight,
            args.dry_run,
            args.analyze_only,
        )
    )
    if public_modes > 1:
        parser.error("choose only one public execution mode")

    if args.worker_cell:
        if args.prepare_inventory or args.keeper_preflight or args.dry_run or args.analyze_only:
            parser.error("--worker-cell cannot be combined with a public mode")
        if not args.expected_inventory_digest:
            parser.error("--worker-cell requires --expected-inventory-digest")
        return run_worker_entry(args.worker_cell, args.expected_inventory_digest)
    if args.expected_inventory_digest:
        parser.error("--expected-inventory-digest is worker-only")

    if args.prepare_inventory:
        inventory = build_source_inventory()
        atomic_write_json(P0_INVENTORY, inventory)
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

    inventory = verify_source_inventory()
    if args.keeper_preflight:
        result = run_keeper_preflight(inventory)
        print(
            json.dumps(
                {
                    "status": result["status"],
                    "context_count": result["context_count"],
                    "n_cells": result["n_cells"],
                    "inventory_digest": result["inventory_digest"],
                },
                indent=2,
            )
        )
        return 0
    if args.dry_run:
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "platform": "Reference",
                    "md_allowed": False,
                    "gpu_allowed": False,
                    "n_cells": len(inventory["cells"]),
                    "inventory_digest": inventory["inventory_digest"],
                    "cells": [row["cell_id"] for row in inventory["cells"]],
                },
                indent=2,
            )
        )
        return 0
    if args.analyze_only:
        results = []
        for cell in make_cells():
            completed = validate_completed_cell(cell, inventory)
            if completed is not None:
                results.append(completed)
        write_summary(results, inventory)
        return 0
    return run_all(inventory)


if __name__ == "__main__":
    raise SystemExit(main())
