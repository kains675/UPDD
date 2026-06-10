#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B production blocker Q2 — Free-leg standalone smoke gate.

Production blocker Q2 (a): hybrid dual-topology requires free-leg
standalone verification because the
5070Ti smoke previously covered only the bound leg. This script:

1. Extracts the cyclic free peptide (chain B) from both Cp4 and WT MD final
   structures (H-complete, no PDBFixer needed for the ncAA path).
2. Builds the free-leg system independently for each endpoint, attaching
   ``ATMForce`` at ``lambda=0.5`` (verdict-specified midpoint).
3. Runs 40-50 step (default 50) MD at lambda=0.5 on the host 5070 Ti.
4. PASS criteria (verdict-defined):
   - Disulfide SG-SG std < 0.1 nm over a 5-step rolling window
     (cyclic_ss signature, per the topology-correction fix).
   - Endpoint XML diff (WT vs MTR): atom-count parity except for the
     residue-4 dual-topology swap (HE1 out, CM+HM1/2/3 in → MTR has +3 atoms).
   - All PE / perturbation energy finite, no NaN, no integrator crash.
5. Writes a JSON report under
   ``outputs/_trackb/freeleg_smoke_<seed>_<timestamp>.json``.

Engine reuse (anti-fragmentation): every build call goes through the existing
``utils/atm_trackB_setup`` API (no duplicated FF stack). The disulfide
detector is the same code path used by the bound-leg smoke.

NOTE: this is the smoke gate — NOT the production ATS free-leg FEP. Production
launch is gated by integrity + execution verification after both blockers
clear.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Optional, List, Dict, Tuple, Any

import numpy as np
import openmm as mm
import openmm.unit as unit
from openmm import app

_PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "utils"))

from atm_trackB_setup import (  # noqa: E402
    HYBRID_MTR_XML,
    ALCH_RESNUM, ALCH_COMMON_ATOM, ALCH_WT_ONLY, ALCH_MTR_ONLY,
    DISULFIDE_MAX_NM,
    resolve_leg_inputs,
    prepare_free_peptide_from_final,
    build_leg_system,
    attach_atm_force,
    identify_alchemical_atoms,
    verify_charge_axis,
)


# --- SG-SG instrumentation -----------------------------------------------
def _sg_indices(modeller, binder_chain: str = "B") -> List[int]:
    out: List[int] = []
    for chain in modeller.topology.chains():
        if chain.id != binder_chain:
            continue
        for res in chain.residues():
            if res.name not in ("CYS", "CYX"):
                continue
            for atom in res.atoms():
                if atom.name == "SG":
                    out.append(atom.index)
    return out


def _measure_sg_distance_nm(positions, sg_pair: Tuple[int, int]) -> float:
    p1 = np.asarray(positions[sg_pair[0]].value_in_unit(unit.nanometer))
    p2 = np.asarray(positions[sg_pair[1]].value_in_unit(unit.nanometer))
    return float(np.linalg.norm(p1 - p2))


# --- Endpoint XML diff ---------------------------------------------------
def _residue4_signature(topology: app.Topology, binder_chain: str = "B",
                        resnum: int = ALCH_RESNUM) -> Dict[str, Any]:
    atoms: List[Dict[str, Any]] = []
    res_name = None
    for chain in topology.chains():
        if chain.id != binder_chain:
            continue
        for res in chain.residues():
            try:
                rn = int(res.id)
            except (TypeError, ValueError):
                continue
            if rn != resnum:
                continue
            res_name = res.name
            for atom in res.atoms():
                atoms.append({
                    "name": atom.name,
                    "element": atom.element.symbol if atom.element else None,
                })
    return {"residue_name": res_name,
            "n_atoms": len(atoms),
            "atom_names": sorted(a["name"] for a in atoms)}


def _topology_signature(topology: app.Topology,
                        binder_chain: str = "B") -> Dict[str, Any]:
    n_atoms = topology.getNumAtoms()
    n_bonds = topology.getNumBonds()
    n_resi = sum(1 for c in topology.chains() if c.id == binder_chain
                 for _ in c.residues())
    return {"n_atoms_total": n_atoms,
            "n_bonds_total": n_bonds,
            "n_binder_residues": n_resi}


def _diff_topology_signatures(wt_sig: Dict[str, Any],
                              mtr_sig: Dict[str, Any]) -> Dict[str, Any]:
    """Compute parity of WT vs MTR free-leg topology.

    Expected diff: only residue-4 atoms differ (HE1 out, CM+HM1-3 in).
    n_atoms_total(MTR) - n_atoms_total(WT) = +3 (one H out, three H in,
    one C in, one H out net: −HE1 +CM +HM1 +HM2 +HM3 = +4 - 1 = +3).
    """
    delta_atoms = mtr_sig["n_atoms_total"] - wt_sig["n_atoms_total"]
    delta_resi  = mtr_sig["n_binder_residues"] - wt_sig["n_binder_residues"]
    expected_delta_atoms = (len(ALCH_MTR_ONLY) - len(ALCH_WT_ONLY))  # +3
    return {
        "delta_n_atoms_total": delta_atoms,
        "delta_n_binder_residues": delta_resi,
        "expected_delta_atoms": expected_delta_atoms,
        "delta_atoms_match_expected": (delta_atoms == expected_delta_atoms),
        "delta_residues_zero": (delta_resi == 0),
    }


def _diff_residue4_signatures(wt_r4: Dict[str, Any],
                              mtr_r4: Dict[str, Any]) -> Dict[str, Any]:
    wt_names = set(wt_r4["atom_names"])
    mtr_names = set(mtr_r4["atom_names"])
    only_wt = sorted(wt_names - mtr_names)
    only_mtr = sorted(mtr_names - wt_names)
    common = sorted(wt_names & mtr_names)
    return {
        "wt_residue": wt_r4["residue_name"],
        "mtr_residue": mtr_r4["residue_name"],
        "only_wt_atoms": only_wt,
        "only_mtr_atoms": only_mtr,
        "n_common_atoms": len(common),
        "ne1_in_common": (ALCH_COMMON_ATOM in common),
        "wt_only_match_expected": (only_wt == sorted(ALCH_WT_ONLY)),
        "mtr_only_match_expected": (only_mtr == sorted(ALCH_MTR_ONLY)),
    }


# --- Smoke run for one endpoint -----------------------------------------
def smoke_one_leg(
    final_pdb: str,
    endpoint: str,
    seed: str,
    work_dir: str,
    n_steps: int = 50,
    lambda_value: float = 0.5,
    binder_chain: str = "B",
    platform_name: str = "CUDA",
    hydrogens_xml: Optional[str] = None,
) -> Dict[str, Any]:
    """Run one free-leg smoke (one endpoint) and report SG-SG + finite + diff.

    The free-peptide PDB is extracted from the H-complete MD ``final.pdb``,
    which preserves the ncAA template atoms exactly (no PDBFixer ncAA path).
    """
    os.makedirs(work_dir, exist_ok=True)
    free_pdb = os.path.join(work_dir, f"{endpoint}_freeleg_{seed}.pdb")
    prepare_free_peptide_from_final(final_pdb, free_pdb,
                                    binder_chain=binder_chain)

    # Build the system from the H-complete free peptide.
    build = build_leg_system(
        free_pdb, leg=f"free_{endpoint}", binder_chain=binder_chain,
        # Solvate=False for smoke (NoCutoff, no PME finite-size effects in the
        # midpoint check). add_hydrogens=False — input is H-complete from MD.
        solvate=False,
        add_hydrogens=False,
        hydrogens_xml=hydrogens_xml,
    )
    system = build["system"]
    modeller = build["modeller"]
    top_sig = _topology_signature(modeller.topology, binder_chain=binder_chain)
    r4_sig = _residue4_signature(modeller.topology,
                                 binder_chain=binder_chain,
                                 resnum=ALCH_RESNUM)

    # Identify SG-SG pair (cyclic_ss).
    sgs = _sg_indices(modeller, binder_chain=binder_chain)
    if len(sgs) < 2:
        return {
            "endpoint": endpoint, "pass": False,
            "error": f"insufficient SG atoms (n={len(sgs)})",
            "topology_signature": top_sig,
            "residue4_signature": r4_sig,
        }
    sg_pair = (sgs[0], sgs[1])

    # Attach ATMForce at lambda=0.5 (verdict midpoint), zero displacement (the
    # dual-topology smoke checks per-leg integrity, not transfer).
    alch = build["alchemical_atoms"]
    displaced = alch["common"] + alch["wt_only"] + alch["mtr_only"]
    attach_atm_force(system, displacement_nm=(0.0, 0.0, 0.0),
                     displaced_atoms=displaced,
                     lambda1=lambda_value, lambda2=lambda_value)

    integrator = mm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin, 1.0 / unit.picosecond, 0.001 * unit.picoseconds
    )
    platform = mm.Platform.getPlatformByName(platform_name)
    sim = app.Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=200)

    # Establish a 5-step rolling baseline at lambda=0.5 to measure SG-SG std.
    rolling_sg: List[float] = []
    rolling_pe: List[float] = []
    pert_log: List[Dict[str, float]] = []
    sim.context.setParameter("Lambda1", lambda_value)
    sim.context.setParameter("Lambda2", lambda_value)

    t0 = time.time()
    crashed = False
    crash_step: Optional[int] = None
    for step in range(n_steps):
        try:
            sim.step(1)
        except Exception as exc:  # OpenMM raises OpenMMException
            crashed = True
            crash_step = step
            crash_msg = str(exc)
            break
        state = sim.context.getState(getEnergy=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        sg = _measure_sg_distance_nm(state.getPositions(), sg_pair)
        rolling_pe.append(pe)
        rolling_sg.append(sg)
        atm = [system.getForce(i) for i in range(system.getNumForces())
               if isinstance(system.getForce(i), mm.ATMForce)][0]
        pert = atm.getPerturbationEnergy(sim.context)
        u1 = (pert[0].value_in_unit(unit.kilocalorie_per_mole)
              if hasattr(pert[0], "value_in_unit") else float(pert[0]))
        u0 = (pert[1].value_in_unit(unit.kilocalorie_per_mole)
              if hasattr(pert[1], "value_in_unit") else float(pert[1]))
        pert_log.append({"step": step, "u1_kcal": u1, "u0_kcal": u0,
                         "delta_kcal": u1 - u0, "pe_kcal": pe, "sg_nm": sg})
    dt = time.time() - t0

    if crashed:
        return {
            "endpoint": endpoint, "pass": False,
            "error": f"integrator crash at step {crash_step}: {crash_msg}",
            "topology_signature": top_sig,
            "residue4_signature": r4_sig,
            "runtime_s": dt,
        }

    sg_arr = np.array(rolling_sg)
    pe_arr = np.array(rolling_pe)
    # Verdict acceptance: SG-SG std < 0.1 nm over a 5-step rolling window.
    if len(sg_arr) >= 5:
        rolling_std = np.array(
            [sg_arr[i:i + 5].std() for i in range(len(sg_arr) - 4)])
        sg_std_max = float(rolling_std.max())
    else:
        sg_std_max = float(sg_arr.std())
    sg_pass = sg_std_max < 0.1  # nm

    all_finite = (bool(np.isfinite(pe_arr).all())
                  and bool(np.isfinite(sg_arr).all())
                  and bool(np.isfinite([p["delta_kcal"] for p in pert_log]).all()))

    return {
        "endpoint": endpoint,
        "seed": seed,
        "free_pdb": free_pdb,
        "lambda": lambda_value,
        "n_steps": n_steps,
        "runtime_s": dt,
        "topology_signature": top_sig,
        "residue4_signature": r4_sig,
        "sg_pair_indices": [int(sg_pair[0]), int(sg_pair[1])],
        "sg_distance_initial_nm": float(rolling_sg[0]) if rolling_sg else None,
        "sg_distance_mean_nm": float(sg_arr.mean()),
        "sg_distance_std_nm": float(sg_arr.std()),
        "sg_rolling_std_max_nm": sg_std_max,
        "sg_acceptance_threshold_nm": 0.1,
        "sg_pass": sg_pass,
        "pe_initial_kcal": float(rolling_pe[0]) if rolling_pe else None,
        "pe_final_kcal": float(rolling_pe[-1]) if rolling_pe else None,
        "pe_min_kcal": float(pe_arr.min()),
        "pe_max_kcal": float(pe_arr.max()),
        "all_finite": all_finite,
        "pass": bool(sg_pass and all_finite),
        "pert_log_first5": pert_log[:5],
        "pert_log_last5": pert_log[-5:],
        "regime": "ranking_only",
    }


# --- Driver ---------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(
        description="Track B Q2 free-leg standalone smoke gate"
    )
    p.add_argument("--seed", default="s7",
                   help="prepared seed under outputs/2QKI_{Cp4_hybrid,WT}_calib_<seed>")
    p.add_argument("--n-steps", type=int, default=50,
                   help="MD steps at lambda=0.5 per endpoint (verdict: 40-50)")
    p.add_argument("--lambda", dest="lam", type=float, default=0.5,
                   help="alchemical midpoint")
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--work-dir", default="outputs/_trackb/freeleg_smoke_work",
                   help="scratch dir for extracted free-peptide PDBs")
    p.add_argument("--report-dir", default="outputs/_trackb",
                   help="JSON report destination directory")
    args = p.parse_args()

    inputs = resolve_leg_inputs(args.seed)
    cp4_final = inputs["final"]["cp4"]
    wt_final = inputs["final"]["wt"]
    hydrogens_xml = inputs.get("hydrogens_xml")
    if not (cp4_final and os.path.isfile(cp4_final)):
        print(f"ERROR: missing Cp4 final.pdb for seed {args.seed}: {cp4_final}",
              file=sys.stderr)
        return 2
    if not (wt_final and os.path.isfile(wt_final)):
        print(f"ERROR: missing WT final.pdb for seed {args.seed}: {wt_final}",
              file=sys.stderr)
        return 2

    # Charge-axis pre-flight: Q1 acceptance must hold.
    axis = verify_charge_axis()
    print("Charge axis (post-Q1 rescale):")
    print(f"  scope = {axis['scope']}, all_within_tol = "
          f"{axis['all_residues_within_tol']}, max|Σq| = "
          f"{axis['max_abs_sigma_q_e']:.2e} e")
    for r in axis["per_residue"]:
        print(f"    {r['residue']:5s}: Σq = {r['sigma_q']:+.2e}  "
              f"NE1 = {r['ne1_charge']}  N = {r['backbone_n_charge']}")
    if not axis["all_residues_within_tol"]:
        print("ERROR: Q1 charge axis FAIL — aborting Q2.", file=sys.stderr)
        return 2
    print()

    # Run both endpoints.
    work = os.path.join(_PROJ_ROOT, args.work_dir)
    os.makedirs(work, exist_ok=True)
    print(f"Free-leg smoke (host {args.platform}): {args.n_steps} step at "
          f"lambda={args.lam}, seed={args.seed}")
    print()

    print("=== Endpoint: MTR (Cp4) ===")
    mtr_result = smoke_one_leg(
        cp4_final, endpoint="mtr", seed=args.seed,
        work_dir=work, n_steps=args.n_steps, lambda_value=args.lam,
        platform_name=args.platform, hydrogens_xml=hydrogens_xml,
    )
    print(json.dumps({k: v for k, v in mtr_result.items()
                      if k not in ("pert_log_first5", "pert_log_last5")},
                     indent=2))
    print()

    print("=== Endpoint: WT (Trp) ===")
    wt_result = smoke_one_leg(
        wt_final, endpoint="wt", seed=args.seed,
        work_dir=work, n_steps=args.n_steps, lambda_value=args.lam,
        platform_name=args.platform, hydrogens_xml=hydrogens_xml,
    )
    print(json.dumps({k: v for k, v in wt_result.items()
                      if k not in ("pert_log_first5", "pert_log_last5")},
                     indent=2))
    print()

    # Endpoint XML diff (atom/bond parity).
    print("=== Endpoint topology diff (WT vs MTR) ===")
    if wt_result.get("topology_signature") and mtr_result.get("topology_signature"):
        top_diff = _diff_topology_signatures(
            wt_result["topology_signature"],
            mtr_result["topology_signature"]
        )
        print(json.dumps(top_diff, indent=2))
    else:
        top_diff = {"error": "missing topology_signature"}

    print()
    print("=== Residue-4 atom diff (WT vs MTR) ===")
    if wt_result.get("residue4_signature") and mtr_result.get("residue4_signature"):
        r4_diff = _diff_residue4_signatures(
            wt_result["residue4_signature"],
            mtr_result["residue4_signature"]
        )
        print(json.dumps(r4_diff, indent=2))
    else:
        r4_diff = {"error": "missing residue4_signature"}

    # PASS gate
    overall_pass = bool(
        wt_result.get("pass") and mtr_result.get("pass")
        and top_diff.get("delta_atoms_match_expected", False)
        and top_diff.get("delta_residues_zero", False)
        and r4_diff.get("ne1_in_common", False)
        and r4_diff.get("wt_only_match_expected", False)
        and r4_diff.get("mtr_only_match_expected", False)
    )

    report = {
        "track": "B",
        "blocker": "Q2_free_leg_standalone_smoke",
        "method_ref": "production blockers cleared, 2026-05-30",
        "seed": args.seed,
        "lambda": args.lam,
        "n_steps": args.n_steps,
        "platform": args.platform,
        "charge_axis": axis,
        "endpoint_mtr": mtr_result,
        "endpoint_wt": wt_result,
        "topology_diff": top_diff,
        "residue4_diff": r4_diff,
        "overall_pass": overall_pass,
        "regime": "ranking_only",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }

    report_dir = os.path.join(_PROJ_ROOT, args.report_dir)
    os.makedirs(report_dir, exist_ok=True)
    report_path = os.path.join(
        report_dir, f"freeleg_smoke_{args.seed}_{time.strftime('%Y%m%dT%H%M%S')}.json")
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)

    print()
    print(f"=== OVERALL ===  {'PASS' if overall_pass else 'FAIL'}")
    print(f"Report: {report_path}")
    return 0 if overall_pass else 1


if __name__ == "__main__":
    sys.exit(main())
