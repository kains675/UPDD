#!/usr/bin/env python3
"""D2 short DCD-on structural probe for Track B tail mechanisms.

This is an analysis-only wrapper. It copies source serialized boxes into a
separate output root, runs one short single-direction ladder per preregistered
cell, writes DCDs for non-solvent atoms plus a local water shell, and summarizes
site-water distances plus site-heavy RMSD. It does not run UWHAM.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import shutil
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = Path(__file__).resolve().parent
OUT_ROOT = REPO_ROOT / "outputs/_trackb/d2_tail_dcd_probe_20260710"

ATM_PYTHON = "/home/san/miniconda3/envs/atm/bin/python"
N_CYCLES = 50
MD_STEPS_PER_CYCLE = 250
DCD_STRIDE_CYCLES = 1
LOCAL_WATER_SHELL_NM = 0.60
WATER_INTRUSION_NM = 0.26
WATER_CONTACT_NM = 0.35
JOBNAME = "trackb"
SOLVENT_RESNAMES = {"HOH", "WAT", "SOL", "NA", "CL", "NA+", "CL-", "K", "K+"}
WATER_RESNAMES = {"HOH", "WAT", "SOL"}

W23A_L2 = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
W23A_L1 = [0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
W4A_L1 = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]


@dataclass(frozen=True)
class Cell:
    cell_id: str
    system: str
    source_rep_dir: str
    endpoint: str
    leg: str
    replicate_index: int
    seed: str
    direction: str
    site_chain: str
    site_resid: str
    lambda1_rampdown: tuple[float, ...] | None
    lambda2_rampup: tuple[float, ...] | None
    control: bool = False


CELLS: tuple[Cell, ...] = (
    Cell(
        cell_id="w23a_s101_bound_dplus",
        system="w23a_gateA",
        source_rep_dir="outputs/_trackb/mdm2_w23a_gateA_20260708/w23a/bound/rep4",
        endpoint="w23a",
        leg="bound",
        replicate_index=4,
        seed="s101",
        direction="dplus",
        site_chain="B",
        site_resid="7",
        lambda1_rampdown=tuple(W23A_L1),
        lambda2_rampup=tuple(W23A_L2),
    ),
    Cell(
        cell_id="w23a_s101_bound_dminus",
        system="w23a_gateA",
        source_rep_dir="outputs/_trackb/mdm2_w23a_gateA_20260708/w23a/bound/rep4",
        endpoint="w23a",
        leg="bound",
        replicate_index=4,
        seed="s101",
        direction="dminus",
        site_chain="B",
        site_resid="7",
        lambda1_rampdown=tuple(W23A_L1),
        lambda2_rampup=tuple(W23A_L2),
    ),
    Cell(
        cell_id="w23a_s101_free_dplus",
        system="w23a_gateA",
        source_rep_dir="outputs/_trackb/mdm2_w23a_gateA_20260708/w23a/free/rep4",
        endpoint="w23a",
        leg="free",
        replicate_index=4,
        seed="s101",
        direction="dplus",
        site_chain="B",
        site_resid="7",
        lambda1_rampdown=tuple(W23A_L1),
        lambda2_rampup=tuple(W23A_L2),
    ),
    Cell(
        cell_id="w23a_ctrl_s42_bound_dminus",
        system="w23a_gateA",
        source_rep_dir="outputs/_trackb/mdm2_w23a_gateA_20260708/w23a/bound/rep3",
        endpoint="w23a",
        leg="bound",
        replicate_index=3,
        seed="s42",
        direction="dminus",
        site_chain="B",
        site_resid="7",
        lambda1_rampdown=tuple(W23A_L1),
        lambda2_rampup=tuple(W23A_L2),
        control=True,
    ),
    Cell(
        cell_id="w4a_carved_s127_free_dplus",
        system="w4a_c1_carved",
        source_rep_dir=(
            "outputs/_trackb/w4a_carved_c1_20260709/"
            "twocopy_w4a_free_FIXAB_carved/wt/free/rep3"
        ),
        endpoint="wt",
        leg="free",
        replicate_index=3,
        seed="s127",
        direction="dplus",
        site_chain="B",
        site_resid="4",
        lambda1_rampdown=tuple(W4A_L1),
        lambda2_rampup=None,
    ),
    Cell(
        cell_id="w4a_carved_s127_bound_dminus",
        system="w4a_c1_carved",
        source_rep_dir=(
            "outputs/_trackb/w4a_carved_c1_20260709/"
            "twocopy_w4a_bound_FIXAB_carved/wt/bound/rep3"
        ),
        endpoint="wt",
        leg="bound",
        replicate_index=3,
        seed="s127",
        direction="dminus",
        site_chain="B",
        site_resid="4",
        lambda1_rampdown=tuple(W4A_L1),
        lambda2_rampup=None,
    ),
    Cell(
        cell_id="w4a_carved_ctrl_s101_free_dplus",
        system="w4a_c1_carved",
        source_rep_dir=(
            "outputs/_trackb/w4a_carved_c1_20260709/"
            "twocopy_w4a_free_FIXAB_carved/wt/free/rep2"
        ),
        endpoint="wt",
        leg="free",
        replicate_index=2,
        seed="s101",
        direction="dplus",
        site_chain="B",
        site_resid="4",
        lambda1_rampdown=tuple(W4A_L1),
        lambda2_rampup=None,
        control=True,
    ),
    Cell(
        cell_id="w4a_carved_ctrl_s163_bound_dminus",
        system="w4a_c1_carved",
        source_rep_dir=(
            "outputs/_trackb/w4a_carved_c1_20260709/"
            "twocopy_w4a_bound_FIXAB_carved/wt/bound/rep4"
        ),
        endpoint="wt",
        leg="bound",
        replicate_index=4,
        seed="s163",
        direction="dminus",
        site_chain="B",
        site_resid="4",
        lambda1_rampdown=tuple(W4A_L1),
        lambda2_rampup=None,
        control=True,
    ),
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_prod():
    return load_module("trackb_inplace_rbfe_production", REPO_ROOT / "scripts/trackb_inplace_rbfe_production.py")


def cell_by_id(cell_id: str) -> Cell:
    for cell in CELLS:
        if cell.cell_id == cell_id:
            return cell
    raise KeyError(f"unknown D2 cell {cell_id!r}")


def selected_cells(raw: str | None) -> list[Cell]:
    if not raw:
        return list(CELLS)
    ids = [item.strip() for item in raw.split(",") if item.strip()]
    return [cell_by_id(cell_id) for cell_id in ids]


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def source_manifest(cell: Cell) -> dict[str, Any]:
    path = REPO_ROOT / cell.source_rep_dir / "run_manifest.json"
    with path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if str(manifest.get("seed")) != cell.seed:
        raise RuntimeError(f"{cell.cell_id}: source seed mismatch {manifest.get('seed')} != {cell.seed}")
    if int(manifest.get("replicate_index")) != int(cell.replicate_index):
        raise RuntimeError(
            f"{cell.cell_id}: source replicate mismatch {manifest.get('replicate_index')} != {cell.replicate_index}"
        )
    if cell.direction not in manifest.get("per_direction", {}):
        raise RuntimeError(f"{cell.cell_id}: source direction {cell.direction} absent from manifest")
    return manifest


def cell_rep_dir(cell: Cell) -> Path:
    return OUT_ROOT / cell.cell_id / cell.endpoint / cell.leg / f"rep{cell.replicate_index}"


def stage_box(cell: Cell) -> dict[str, str]:
    source = REPO_ROOT / cell.source_rep_dir
    dest = cell_rep_dir(cell)
    dest.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    for suffix in (f"inplace_rbfe_{cell.leg}_sys.xml", f"inplace_rbfe_{cell.leg}.pdb"):
        src = source / suffix
        if not src.is_file():
            raise FileNotFoundError(f"{cell.cell_id}: missing source box file {src}")
        dst = dest / suffix
        if dst.exists():
            if dst.stat().st_size != src.stat().st_size:
                raise RuntimeError(f"{cell.cell_id}: staged box collision with different size: {dst}")
        else:
            shutil.copy2(src, dst)
        paths[suffix] = str(dst)
    stamp = dest / f".seed_{cell.seed}"
    stamp.write_text(cell.seed + "\n", encoding="utf-8")
    return paths


def quantity_positions_nm(positions) -> np.ndarray:
    import openmm.unit as unit

    arr = positions.value_in_unit(unit.nanometer)
    return np.asarray([[v.x, v.y, v.z] for v in arr], dtype=float)


def box_lengths_nm(topology) -> np.ndarray | None:
    import openmm.unit as unit

    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        return None
    vals = vectors.value_in_unit(unit.nanometer)
    mat = np.asarray([[v.x, v.y, v.z] for v in vals], dtype=float)
    if np.allclose(mat, np.diag(np.diag(mat)), atol=1e-6):
        return np.diag(mat)
    return None


def min_distances_nm(
    a: np.ndarray, b: np.ndarray, lengths: np.ndarray | None = None
) -> np.ndarray:
    delta = a[:, None, :] - b[None, :, :]
    if lengths is not None:
        delta -= lengths * np.round(delta / lengths)
    return np.sqrt(np.sum(delta * delta, axis=2))


def is_hydrogen(atom) -> bool:
    if atom.element is None:
        return atom.name.upper().startswith("H")
    return atom.element.symbol.upper() == "H"


def select_atoms(prod, rbfe, cell: Cell, loaded: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    topology = loaded["topology"]
    positions = quantity_positions_nm(loaded["positions"])
    lengths = box_lengths_nm(topology)
    non_solvent: set[int] = set()
    site_heavy: list[int] = []
    water_residues: dict[int, list[int]] = {}
    water_oxygens: dict[int, int] = {}
    site_residues: list[dict[str, Any]] = []

    for atom in topology.atoms():
        resname = atom.residue.name.upper()
        if resname not in SOLVENT_RESNAMES:
            non_solvent.add(int(atom.index))
        if (
            atom.residue.chain.id == cell.site_chain
            and str(atom.residue.id) == str(cell.site_resid)
            and resname in {"TRP", "ALA", "MTR", "VAL", "ILE", "GLY"}
        ):
            if not is_hydrogen(atom):
                site_heavy.append(int(atom.index))
            key = {
                "residue_index": atom.residue.index,
                "chain": atom.residue.chain.id,
                "resid": atom.residue.id,
                "resname": atom.residue.name,
            }
            if key not in site_residues:
                site_residues.append(key)
        if resname in WATER_RESNAMES:
            water_residues.setdefault(atom.residue.index, []).append(int(atom.index))
            if atom.name.upper() in {"O", "OW", "OH2"}:
                water_oxygens[atom.residue.index] = int(atom.index)

    if not site_heavy:
        raise RuntimeError(f"{cell.cell_id}: no site heavy atoms selected for B:{cell.site_resid}")
    site_xyz = positions[site_heavy]
    local_water_residue_ids: list[int] = []
    for residx, atom_indices in water_residues.items():
        xyz = positions[atom_indices]
        if float(np.min(min_distances_nm(xyz, site_xyz, lengths))) <= LOCAL_WATER_SHELL_NM:
            local_water_residue_ids.append(residx)

    local_water_atoms: set[int] = set()
    local_water_oxygen_indices: list[int] = []
    for residx in local_water_residue_ids:
        local_water_atoms.update(water_residues[residx])
        if residx in water_oxygens:
            local_water_oxygen_indices.append(water_oxygens[residx])

    selected = sorted(non_solvent | local_water_atoms)
    orig_to_subset = {orig: idx for idx, orig in enumerate(selected)}
    site_heavy_subset = [orig_to_subset[i] for i in site_heavy if i in orig_to_subset]
    water_oxygen_subset = [orig_to_subset[i] for i in local_water_oxygen_indices if i in orig_to_subset]

    sub_topology = rbfe._subset_topology(topology, selected)
    import openmm.unit as unit
    from openmm import Vec3

    subset_xyz_nm = quantity_positions_nm(loaded["positions"])[selected]
    subset_positions = [Vec3(float(x), float(y), float(z)) for x, y, z in subset_xyz_nm] * unit.nanometer
    import openmm.app as app

    subset_pdb = out_dir / "dcd_subset_topology.pdb"
    with subset_pdb.open("w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(sub_topology, subset_positions, handle)

    meta = {
        "cell_id": cell.cell_id,
        "selection_mode": "non_solvent_plus_initial_local_water_shell",
        "local_water_shell_nm": LOCAL_WATER_SHELL_NM,
        "n_total_atoms": int(loaded["n_atoms"]),
        "n_non_solvent_atoms": len(non_solvent),
        "n_site_heavy_atoms": len(site_heavy),
        "n_local_water_residues": len(local_water_residue_ids),
        "n_local_water_atoms": len(local_water_atoms),
        "n_local_water_oxygens": len(local_water_oxygen_indices),
        "n_selected_atoms": len(selected),
        "site_residues": site_residues,
        "selected_indices_orig": selected,
        "site_heavy_orig": site_heavy,
        "water_oxygen_orig": local_water_oxygen_indices,
        "site_heavy_subset": site_heavy_subset,
        "water_oxygen_subset": water_oxygen_subset,
        "subset_topology_pdb": str(subset_pdb),
    }
    (out_dir / "dcd_selection.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return meta


def build_schedule(prod, rbfe, cell: Cell) -> dict[str, Any]:
    return prod._build_single_direction_schedule(
        rbfe,
        construction="twocopy",
        direction=prod.DIRECTION_OF_TAG[cell.direction],
        n_windows_half=6,
        softcore_band=2,
        n_apex_bridge=0,
        apex_band=0.5,
        lambda1_rampdown=list(cell.lambda1_rampdown) if cell.lambda1_rampdown else None,
        lambda2_rampup=list(cell.lambda2_rampup) if cell.lambda2_rampup else None,
    )


def summarize_values(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0, "min": None, "p01": None, "p50": None, "p95": None, "max": None}
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(arr.size),
        "min": float(np.min(arr)),
        "p01": float(np.quantile(arr, 0.01)),
        "p50": float(np.quantile(arr, 0.50)),
        "p95": float(np.quantile(arr, 0.95)),
        "max": float(np.max(arr)),
    }


def parse_out_states(path: Path) -> list[int]:
    states: list[int] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            cols = stripped.split()
            if len(cols) >= 1:
                states.append(int(float(cols[0])))
    return states


def analyze_cell_dcd(cell: Cell, direction_dir: Path, selection: dict[str, Any]) -> dict[str, Any]:
    import mdtraj as md

    topology_pdb = Path(selection["subset_topology_pdb"])
    site = list(selection["site_heavy_subset"])
    waters = list(selection["water_oxygen_subset"])
    water_minima: list[float] = []
    site_rmsd: list[float] = []
    frame_mismatch: list[str] = []
    by_walker: list[dict[str, Any]] = []
    n_frames = 0

    for dcd_path in sorted((direction_dir / "dcd").glob("r*/trackb_*.dcd")):
        walker = dcd_path.parent.name
        out_path = direction_dir / walker / f"trackb_{cell.direction}.out"
        states = parse_out_states(out_path)
        traj = md.load(str(dcd_path), top=str(topology_pdb))
        n_frames += int(traj.n_frames)
        if traj.n_frames != len(states):
            frame_mismatch.append(f"{walker}: frames={traj.n_frames} out_rows={len(states)}")
        xyz = traj.xyz

        walker_water_min: list[float] = []
        if waters:
            site_xyz = xyz[:, site, :]
            water_xyz = xyz[:, waters, :]
            diff = site_xyz[:, :, None, :] - water_xyz[:, None, :, :]
            dist = np.sqrt(np.sum(diff * diff, axis=3))
            mins = np.min(dist, axis=(1, 2))
            water_minima.extend(float(v) for v in mins)
            walker_water_min = [float(v) for v in mins]

        site_slice = traj.atom_slice(site)
        rmsd = md.rmsd(site_slice, site_slice, frame=0)
        site_rmsd.extend(float(v) for v in rmsd)

        by_walker.append(
            {
                "walker": walker,
                "n_frames": int(traj.n_frames),
                "states_first10": states[:10],
                "states_last10": states[-10:],
                "water_min_nm": summarize_values(walker_water_min),
                "site_rmsd_nm": summarize_values([float(v) for v in rmsd]),
            }
        )

    if frame_mismatch:
        raise RuntimeError(f"{cell.cell_id}: DCD/.out frame mismatch: {frame_mismatch}")

    water_summary = summarize_values(water_minima)
    if water_minima:
        water_summary["frac_lt_0p26"] = float(np.mean(np.asarray(water_minima) < WATER_INTRUSION_NM))
        water_summary["frac_lt_0p35"] = float(np.mean(np.asarray(water_minima) < WATER_CONTACT_NM))
    else:
        water_summary["frac_lt_0p26"] = None
        water_summary["frac_lt_0p35"] = None

    return {
        "cell_id": cell.cell_id,
        "n_frames": n_frames,
        "water_min_site_heavy_nm": water_summary,
        "site_heavy_rmsd_nm": summarize_values(site_rmsd),
        "by_walker": by_walker,
    }


def run_cell(prod, rbfe, cell: Cell, *, platform: str) -> dict[str, Any]:
    driver = prod._load_driver()
    source = source_manifest(cell)
    staged = stage_box(cell)
    rep_dir = cell_rep_dir(cell)
    direction_dir = rep_dir / cell.direction
    direction_dir.mkdir(parents=True, exist_ok=True)

    xml = rep_dir / f"inplace_rbfe_{cell.leg}_sys.xml"
    pdb = rep_dir / f"inplace_rbfe_{cell.leg}.pdb"
    loaded = rbfe.load_serialized_system(str(xml), str(pdb))
    schedule = build_schedule(prod, rbfe, cell)
    selection = select_atoms(prod, rbfe, cell, loaded, direction_dir)

    log_path = direction_dir / f"{JOBNAME}_{cell.direction}_driver.log"
    cntl_path = direction_dir / f"{JOBNAME}_{cell.direction}_asyncre.cntl"
    prod._write_direction_cntl(str(cntl_path), schedule, cell.direction, MD_STEPS_PER_CYCLE, 1.0)
    dcd_dir = direction_dir / "dcd"
    rng_seed = 20260613 + 1000 * (cell.replicate_index + 1)

    t0 = time.time()
    ladder = rbfe.InplaceRbfeLadder(
        loaded["system"],
        loaded["positions"],
        schedule,
        platform_name=platform,
        temperature_K=schedule["temperature_K"],
        timestep_fs=1.0,
        log_path=str(log_path),
        seed=rng_seed,
        minimize_iters=500,
        backward_equil_steps=500,
        out_dir=str(direction_dir),
        out_basename=f"{JOBNAME}_{cell.direction}",
        staged_min=True,
        dcd_dir=str(dcd_dir),
        dcd_topology=loaded["topology"],
        dcd_atom_indices=selection["selected_indices_orig"],
        dcd_stride_cycles=DCD_STRIDE_CYCLES,
    )
    per_cycle: list[dict[str, Any]] = []
    nan_any = False
    nan_states: set[int] = set()
    apex_state = int(schedule["n_states"]) - 1
    try:
        for _ in range(N_CYCLES):
            info = ladder.run_cycle(md_steps=MD_STEPS_PER_CYCLE)
            nan_any = nan_any or bool(info["nan_seen"])
            nan_states.update(int(s) for s in info.get("nan_states", []))
            rs = info["replica_state"]
            per_cycle.append(
                {
                    "cycle": int(info["cycle"]),
                    "n_accepted": int(info["n_accepted"]),
                    "n_pairs": int(info["n_pairs"]),
                    "nan_seen": bool(info["nan_seen"]),
                    "nan_states": [int(s) for s in info.get("nan_states", [])],
                    "apex_occupied": apex_state in rs,
                    "endpoint_occupied": 0 in rs,
                }
            )
            if info["nan_seen"]:
                raise RuntimeError(f"{cell.cell_id}: NaN states {sorted(nan_states)} at cycle {info['cycle']}")
    finally:
        ladder.close()
    elapsed = time.time() - t0

    tr = driver.parse_state_transitions_from_log(str(log_path), warmup_cycles=2)
    mix = driver.check_atm_mixing(
        str(log_path), schedule_K=schedule["n_states"], warmup_cycles=2, min_crossings=1
    )
    adj = tr.get("adjacent_crossings", {}) or {}
    dcd_analysis = analyze_cell_dcd(cell, direction_dir, selection)

    manifest = {
        "cell": asdict(cell),
        "source_manifest": {
            "path": str(REPO_ROOT / cell.source_rep_dir / "run_manifest.json"),
            "seed": source.get("seed"),
            "replicate_index": source.get("replicate_index"),
            "direction_source_summary": source.get("per_direction", {}).get(cell.direction, {}),
        },
        "staged_box": staged,
        "run": {
            "n_cycles": N_CYCLES,
            "md_steps_per_cycle": MD_STEPS_PER_CYCLE,
            "dcd_stride_cycles": DCD_STRIDE_CYCLES,
            "platform": platform,
            "rng_seed": rng_seed,
            "elapsed_s": elapsed,
            "nan_any": nan_any,
            "nan_states": sorted(nan_states),
            "apex_occupancy_cycles": sum(1 for c in per_cycle if c["apex_occupied"]),
            "endpoint_occupancy_cycles": sum(1 for c in per_cycle if c["endpoint_occupied"]),
        },
        "schedule": {
            "n_states": schedule["n_states"],
            "schedule_kind": schedule.get("schedule_kind"),
            "lambda1_rampdown": list(cell.lambda1_rampdown) if cell.lambda1_rampdown else None,
            "lambda2_rampup": list(cell.lambda2_rampup) if cell.lambda2_rampup else None,
        },
        "mixing": {
            "total_round_trips": mix.get("total_round_trips"),
            "both_ends_visited_count": mix.get("both_ends_visited_count"),
            "walls": mix.get("walls"),
            "gate_passed": mix.get("passed"),
            "n_adjacent_pairs": int(schedule["n_states"]) - 1,
            "n_adjacent_with_crossings": sum(1 for v in adj.values() if v > 0),
            "adjacent_crossings": adj,
        },
        "selection": selection,
        "dcd_analysis": dcd_analysis,
        "per_cycle": per_cycle,
    }
    (rep_dir / "d2_cell_manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n", encoding="utf-8")
    return manifest


def write_summary(manifests: list[dict[str, Any]]) -> None:
    summary = {
        "analysis": "d2_tail_dcd_probe_20260710",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "out_root": str(OUT_ROOT),
        "n_cells": len(manifests),
        "cells": [
            {
                "cell_id": m["cell"]["cell_id"],
                "system": m["cell"]["system"],
                "seed": m["cell"]["seed"],
                "control": m["cell"]["control"],
                "leg": m["cell"]["leg"],
                "direction": m["cell"]["direction"],
                "n_local_water_residues": m["selection"]["n_local_water_residues"],
                "n_selected_atoms": m["selection"]["n_selected_atoms"],
                "n_frames": m["dcd_analysis"]["n_frames"],
                "water_min_nm": m["dcd_analysis"]["water_min_site_heavy_nm"],
                "site_rmsd_nm": m["dcd_analysis"]["site_heavy_rmsd_nm"],
                "mixing": m["mixing"],
                "elapsed_s": m["run"]["elapsed_s"],
            }
            for m in manifests
        ],
    }
    (ANALYSIS_DIR / "d2_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# D2 Tail DCD Probe Report",
        "",
        f"Generated: {summary['generated']}",
        "",
        "| cell | ctrl | waters | frames | min water nm | frac<0.26 | frac<0.35 | site RMSD p95 nm | min crossing | RT |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in summary["cells"]:
        water = row["water_min_nm"]
        rmsd = row["site_rmsd_nm"]
        crossings = row["mixing"]["adjacent_crossings"] or {}
        min_cross = min(crossings.values()) if crossings else None
        lines.append(
            "| {cell_id} | {control} | {waters} | {frames} | {wmin} | {f026} | {f035} | {r95} | {minc} | {rt} |".format(
                cell_id=row["cell_id"],
                control="yes" if row["control"] else "",
                waters=row["n_local_water_residues"],
                frames=row["n_frames"],
                wmin=fmt(water["min"]),
                f026=fmt(water["frac_lt_0p26"]),
                f035=fmt(water["frac_lt_0p35"]),
                r95=fmt(rmsd["p95"]),
                minc=fmt(min_cross, 0),
                rt=fmt(row["mixing"]["total_round_trips"], 0),
            )
        )
    lines.append("")
    lines.append("Interpretation is deferred to `PATH_D2_DIAGNOSIS.md` after Runner verification.")
    (ANALYSIS_DIR / "d2_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value)
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return str(value)


def dry_run(cells: list[Cell]) -> None:
    rows = []
    for cell in cells:
        src = REPO_ROOT / cell.source_rep_dir
        rows.append(
            {
                "cell_id": cell.cell_id,
                "source_exists": src.is_dir(),
                "xml": (src / f"inplace_rbfe_{cell.leg}_sys.xml").is_file(),
                "pdb": (src / f"inplace_rbfe_{cell.leg}.pdb").is_file(),
                "manifest": (src / "run_manifest.json").is_file(),
                "dest": str(cell_rep_dir(cell)),
            }
        )
    print(json.dumps(rows, indent=2))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", default=None, help="comma-separated cell IDs; default all preregistered cells")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--platform", default="CUDA")
    args = parser.parse_args(argv)

    cells = selected_cells(args.cells)
    if args.dry_run:
        dry_run(cells)
        return 0

    prod = load_prod()
    rbfe = prod._load_rbfe()
    manifests: list[dict[str, Any]] = []
    for cell in cells:
        if args.analyze_only:
            manifest_path = cell_rep_dir(cell) / "d2_cell_manifest.json"
            with manifest_path.open("r", encoding="utf-8") as handle:
                manifests.append(json.load(handle))
        else:
            print(f"[D2] running {cell.cell_id} ({cell.leg}/{cell.direction}, seed={cell.seed})", flush=True)
            manifests.append(run_cell(prod, rbfe, cell, platform=args.platform))
            write_summary(manifests)
    write_summary(manifests)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
