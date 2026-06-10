#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""DEPRECATED — Track B v1 production launcher (architecturally invalid).

**2026-05-31 DEPRECATION (archived, not deleted)**: This launcher uses a single
OpenMM Context per replica and flips ``Direction`` mid-trajectory at
λ=0.5. The ATMForce reference potential is
``select(step(Direction), u0, u1)`` — a step-function discontinuity in the
base potential when Direction flips while the same Context's force cache /
constraint state is reused without re-initialization. v1 produced NaN at r5
on the d=+1→d=-1 transition (3 reproducible crashes 2026-05-30T23:57 →
2026-05-31T00:55); engineering Phase 1 (1 fs warmup, 5000-iter
re-minimization) FAILED to mitigate.

Analysis classifies this as an architectural blunder and prescribes Option 2:
upstream ``atom_openmm.async_re`` with separate replica walkers per
direction (Gallicchio 2021 JCIM 61:5424 standard).

**USE INSTEAD**: ``scripts/trackb_production_v2_asyncre.py``.

v1 archived under
``outputs/_trackb/_archive/v1_architectural_failed_20260531/`` with manifest
README — the 11 r0..r10 forward (d=+1) ``.out`` files are SCRAP for
production analysis (UWHAM bidirectional estimator requires reciprocal d=-1
data, Tan 2012 JCP 136:144102; Mey 2020 LiveCoMS §5.2 forbids
forward-only).

================================================================================
Original docstring preserved for historical context only:

Track B production launcher — ATM/ATS alchemical MTR↔Trp ΔΔG_bind.

Q1 + Q2 production blockers cleared 2026-05-30 (smoke gate PASS on host
5070Ti).

# Design — per-endpoint ABFE → paired ΔΔG_bind

ATS implementation per the existing infrastructure builds **one hybrid system
per endpoint** (Cp4 / WT) where:

* The MTR ncAA hybrid XML (params/MTR_gaff2_hybrid.xml, Option-β/regime-2,
  per-residue Σq |≤5e-4| e after Q1 rescale, NE1 frozen at −0.3418) provides
  the *charge regime* in both endpoints.
* Each endpoint's residue-4 sidechain is the **actual present sidechain**
  (MTR atoms in Cp4, TRP atoms in WT) — the alchemical region in the
  per-endpoint system is the sidechain itself (the smoke verified +3 atom
  parity between endpoints).
* Two legs per endpoint: bound (2QKI complex) and free (cyclic peptide
  extracted, cyclic_ss disulfide preserved).

ΔΔG cycle (verbatim per ``utils/atm_trackB_setup`` module docstring):

    ΔG_bind(endpoint) = ΔG_alch(bound, λ=0→1) − ΔG_alch(free, λ=0→1)
    ΔΔG_bind_MTR↔Trp  = ΔG_bind(Cp4) − ΔG_bind(WT)

Each leg's alchemical scan uses the AToM 22-state symmetric schedule
(11 forward + 11 backward, intermediate at λ=0.5) so that uwham analyzes
the standard ATM trajectory format natively.

# Hardware contract (ADR-0002 sequential allocation)

* CUDA_VISIBLE_DEVICES forced to 0 (host 5070Ti) at launch — V100 stays with
  Track A.
* Single GPU per λ window, sequential through the schedule. Track A is
  presently running on the VM (192.168.122.155); this launcher does not
  touch that path.

# Output schema (ranking-only)

Per-leg per-replicate ``r{i}/trackb.out`` files in the AToM-OpenMM standard
11-column format::

    stateid temperature direction lambda1 lambda2 alpha u0 w0 potE pertE bias

Analysis via ``atom_openmm.uwham.calculate_uwham()`` → per-leg ΔG_alch +
SE → ΔG_bind(endpoint) → ΔΔG_bind_paired. The final JSON
``ddint_kcal_alchemical.json`` records μ, σ_btwn (across 3 replicates), SE,
z_SE, 95% CI, the verbatim launch command, env snapshot, deterministic seed
list, and pre-registered 4 outcomes (sign-stable / σ drift / |μ| drift /
sign-flip).

# Pre-registered outcomes (pre-launch checklist #3)

Multi-outcome pre-registration pattern.

1. **sign-stable**: ΔΔG sign matches Magotti SSOT direction
   (MTR less-favorable than WT, ΔΔG > 0), σ_btwn ≤ 0.5 kcal/mol — primary
   ranking-positive outcome.
2. **σ drift**: σ_btwn > 0.7 kcal/mol (>140% of expected upper bound) —
   replicate count expansion required before ranking claim.
3. **|μ| drift**: |ΔΔG| outside [0.5, 5.0] kcal/mol — methodology
   sensitivity review (cyclic_ss force-field bias check).
4. **sign-flip**: ΔΔG opposite to Magotti SSOT direction — escalate for
   charge regime / sampling regime re-verification.

Ranking-only: NO absolute Magotti comparison. Only sign and rank are
interpreted.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import random
import shlex
import shutil
import subprocess
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
    resolve_leg_inputs,
    prepare_free_peptide_from_final,
    build_leg_system,
    identify_alchemical_atoms,
    verify_charge_axis,
)

# ---------------------------------------------------------------------------
# Canonical AToM-OpenMM 22-state schedule (verdict pre-launch §1).
#
# 11 forward (direction=+1, lambda 0→0.5) + 11 backward (direction=−1,
# lambda 0.5→0). intermediate=1 at the two λ=0.5 midpoints. w0=1 at the
# midpoints (per uwham default). alpha=0.10 (kcal/mol)^-1. u0=110 kcal/mol.
# Identical defaults to atom_openmm.uwham.calculate_uwham() initialiser.
# ---------------------------------------------------------------------------
LAMBDA_FWD = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
LAMBDA_BWD = [0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
DIRECTIONS = [1] * 11 + [-1] * 11
LAMBDAS_1   = LAMBDA_FWD + LAMBDA_BWD
LAMBDAS_2   = LAMBDA_FWD + LAMBDA_BWD
INTERMD     = [0] * 10 + [1, 1] + [0] * 10           # midpoint flag
W0          = [0.0] * 10 + [1.0, 1.0] + [0.0] * 10   # midpoint bias
ALPHA       = [0.10] * 22  # 1/(kcal/mol)
U0          = [110.0] * 22  # kcal/mol
TEMP_K      = 300.0
N_STATES    = len(DIRECTIONS)
assert N_STATES == 22, "Schedule must be 22 states (11+11 symmetric)"


# ---------------------------------------------------------------------------
# ATMForce wiring — per-endpoint ABFE.
#
# For ABFE the displacement vector decouples the *ligand* (the binder chain)
# from the receptor. At λ=0 the ligand is fully coupled (potential = u0);
# at λ=1 it is displaced (potential = u1, scaled by 1−lambda1/2 outside the
# pocket via the soft-core potential). The 22-state symmetric schedule
# integrates ΔG via the AToM ilogistic soft-core.
# ---------------------------------------------------------------------------
DISPLACEMENT_NM_DEFAULT = (2.5, 0.0, 0.0)  # 2.5 nm displacement (verdict B-C5 spirit)

# Soft-core cap constants — must match attach_atm_force_production defaults.
# uwham's _bias_fcn expects perturbation energies already capped to
# [UB_KCAL, UMAX_KCAL] via the AtomUtils softCorePertE smooth bound.
UMAX_KCAL = 200.0
UBCORE_KCAL = 100.0
ACORE = 0.062500


def attach_atm_force_production(
    system: mm.System,
    displacement_nm: Tuple[float, float, float],
    displaced_atoms: List[int],
    lambda1: float = 0.0,
    lambda2: float = 0.0,
    alpha_per_kcal: float = 0.10,
    u0_kcal: float = 110.0,
    w0_kcal: float = 0.0,
    direction: float = 1.0,
    umax_kcal: float = UMAX_KCAL,
    ubcore_kcal: float = UBCORE_KCAL,
    acore: float = ACORE,
) -> int:
    """Add ATMForce with full 9-parameter ABFE wiring.

    Identical to ``atm_trackB_setup.attach_atm_force`` but exposes
    ``direction`` (required for the symmetric backward half of the 22-state
    schedule) and uses the AToM canonical u0=110 / alpha=0.10 defaults
    matched to ``atom_openmm.uwham`` initialiser.
    """
    import copy
    kcal = unit.kilocalorie_per_mole
    kj = unit.kilojoule_per_mole

    kcal_per_kj = (1.0 * kcal).value_in_unit(kj)  # = 4.184
    alpha = (alpha_per_kcal / kcal_per_kj) / kj
    uh = (u0_kcal * kcal).value_in_unit(kj) * kj
    w0 = (w0_kcal * kcal).value_in_unit(kj) * kj
    umax = (umax_kcal * kcal).value_in_unit(kj) * kj
    ubcore = (ubcore_kcal * kcal).value_in_unit(kj) * kj

    atm = mm.ATMForce(lambda1, lambda2, alpha, uh, w0, umax, ubcore, acore,
                      float(direction))

    move_types = (mm.NonbondedForce, mm.HarmonicBondForce,
                  mm.HarmonicAngleForce, mm.PeriodicTorsionForce)
    to_move = [i for i in range(system.getNumForces())
               if isinstance(system.getForce(i), move_types)]
    for i in to_move:
        atm.addForce(copy.copy(system.getForce(i)))
    for i in sorted(to_move, reverse=True):
        system.removeForce(i)

    disp = mm.Vec3(*displacement_nm) * unit.nanometer
    zero = mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer
    dset = set(displaced_atoms or [])
    for idx_p in range(system.getNumParticles()):
        atm.addParticle(disp if idx_p in dset else zero)

    return system.addForce(atm)


# ---------------------------------------------------------------------------
# Ligand displaced-atoms (entire binder chain, for ABFE).
# ---------------------------------------------------------------------------
def collect_binder_atoms(topology: app.Topology,
                         binder_chain: str = "B") -> List[int]:
    """Return atom indices of the entire binder chain (the displaced ligand).

    For ABFE the binder chain is moved out of the binding pocket by the
    displacement vector; the receptor/solvent atoms stay put.
    """
    idxs: List[int] = []
    for chain in topology.chains():
        if chain.id != binder_chain:
            continue
        for res in chain.residues():
            for atom in res.atoms():
                idxs.append(atom.index)
    return idxs


# ---------------------------------------------------------------------------
# AToM-OpenMM .out format writer (matches ommreplica.OMMReplicaATM.save_out).
# ---------------------------------------------------------------------------
def write_out_line(fh, stateid: int, temperature_K: float, direction: float,
                   lambda1: float, lambda2: float, alpha_per_kj: float,
                   uh_kj: float, w0_kj: float, pot_energy_kj: float,
                   pert_energy_kj: float, bias_energy_kj: float) -> None:
    """Write one line in the canonical 11-column AToM .out format.

    Units written are kcal/mol (energies) and 1/(kcal/mol) (alpha), matching
    ``atom_openmm.uwham.calculate_uwham`` parser expectations.
    """
    kj = unit.kilojoule_per_mole.conversion_factor_to(unit.kilocalorie_per_mole)  # 1/4.184
    line = ("%d %f %f %f %f %f %f %f %f %f %f\n" % (
        stateid,
        temperature_K,
        direction,
        lambda1,
        lambda2,
        alpha_per_kj / kj,                    # 1/(kcal/mol)
        uh_kj * kj,                           # kcal/mol
        w0_kj * kj,                           # kcal/mol
        pot_energy_kj * kj,                   # kcal/mol
        pert_energy_kj * kj,                  # kcal/mol
        bias_energy_kj * kj,                  # kcal/mol
    ))
    fh.write(line)


def soft_core_pert_e(u: float, umax: float, ub: float, a: float) -> float:
    """Apply the AToM soft-core cap to a raw u1-u0 (kJ/mol).

    Verbatim port of ``atom_openmm.utils.AtomUtils.softCorePertE``. The cap
    bounds ``u`` to ``[ub, umax]`` smoothly when ``u > ub``. uwham's
    ``_bias_fcn`` expects the **capped** perturbation energy (this is what
    ``ommworker`` writes to the .out file). Without the cap, displaced
    shadow ligands with severe overlap produce 1e15-scale raw u1-u0 values
    that propagate as infinities through the uwham softplus.
    """
    if u <= ub:
        return u
    gu = (u - ub) / (a * (umax - ub))   # y/alpha
    zeta = 1.0 + 2.0 * gu * (gu + 1.0)
    zetap = zeta ** a
    return (umax - ub) * (zetap - 1.0) / (zetap + 1.0) + ub


def softplus_bias_kj(lambda1: float, lambda2: float, alpha_per_kj: float,
                     uh_kj: float, w0_kj: float, pert_kj: float) -> float:
    """Compute the AToM soft-core bias energy (kJ/mol) at one sample.

    Mirrors ``openmm_async_re._softplus`` exactly so the .out file's
    bias_energy column is computed in code (matches what the live ATMForce
    applies). The integrator already incorporates the bias; this is the
    *audit log* of what the bias was.
    """
    import math
    softplus = lambda2 * pert_kj + w0_kj
    if alpha_per_kj > 0:
        ee = 1.0 + math.exp(-alpha_per_kj * (pert_kj - uh_kj))
        softplus += ((lambda2 - lambda1) / alpha_per_kj) * math.log(ee)
    return softplus


# ---------------------------------------------------------------------------
# Production run for one (endpoint, leg, replicate).
# ---------------------------------------------------------------------------
def run_one_replicate(
    pdb_path: str,
    endpoint: str,
    leg: str,
    replicate: int,
    seed: int,
    out_dir: str,
    n_lambda: int = 22,
    equil_ps: float = 500.0,                  # 0.5 ns equilibration per replicate
    prod_ns_per_window: float = 5.0,          # 5 ns per λ window
    sample_interval_ps: float = 50.0,         # sample every 50 ps
    minimize_iter: int = 5000,
    timestep_fs: float = 2.0,
    binder_chain: str = "B",
    platform_name: str = "CUDA",
    displacement_nm: Tuple[float, float, float] = DISPLACEMENT_NM_DEFAULT,
    hydrogens_xml: Optional[str] = None,
    solvate: bool = True,
    jobname: str = "trackb",
    progress_cb=None,
    state_ids: Optional[List[int]] = None,
    warmup_ps: float = 1.0,
    warmup_timestep_fs: float = 1.0,
    direction_flip_warmup_ps: float = 0.1,
    direction_flip_warmup_fs: float = 1.0,
    prod_timestep_fs_after_flip: Optional[float] = None,
    resume_from_state: int = 0,
    resume_skip_equilibration: bool = False,
) -> Dict[str, Any]:
    """Run one replicate of one leg of one endpoint — the production unit.

    Walks the 22-state schedule sequentially (no replica exchange — single-GPU
    sequential per the hardware contract). Each window:
      * sets Lambda1/Lambda2/Direction/Alpha/Uh/W0 via context parameters
      * runs equilibration steps (first window only) or transition steps
      * samples potE/pertE at sample_interval_ps intervals for prod_ns_per_window
    Writes ``r{replicate}/{jobname}.out`` in the canonical 11-column format.
    """
    # AToM layout: ``r{state_id}/{jobname}.out`` (one file per state, replicates
    # appended). ``calculate_uwham`` expects this exact path. Per-replicate log
    # goes to ``_logs/replicate{r}.log`` at out_dir.
    log_dir = os.path.join(out_dir, "_logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"replicate{replicate}_seed{seed}.log")
    # Per-state directories are created on demand inside the schedule loop.

    t_start = time.time()

    # Build the system (solvated PME) — once per replicate (cached topology).
    # Hydrogen-build policy: free leg uses MD final.pdb (H-complete), so
    # skip addHydrogens to avoid mis-placing the MTR methyl set. Bound leg
    # uses _md_input/2QKI_*.pdb (pre-MD, H-incomplete for ncAA), so
    # addHydrogens with the MTR_hydrogens.xml is required.
    add_h = (leg == "bound")
    build = build_leg_system(
        pdb_path, leg=f"{leg}_{endpoint}", binder_chain=binder_chain,
        solvate=solvate,
        add_hydrogens=add_h,
        hydrogens_xml=hydrogens_xml,
    )
    system = build["system"]
    modeller = build["modeller"]

    displaced = collect_binder_atoms(modeller.topology,
                                     binder_chain=binder_chain)
    if not displaced:
        raise RuntimeError(
            f"No displaced atoms collected (binder_chain={binder_chain})"
        )

    # Attach ATMForce with the initial λ=0, direction=+1 state.
    attach_atm_force_production(
        system,
        displacement_nm=displacement_nm,
        displaced_atoms=displaced,
        lambda1=LAMBDAS_1[0], lambda2=LAMBDAS_2[0],
        alpha_per_kcal=ALPHA[0], u0_kcal=U0[0], w0_kcal=W0[0],
        direction=float(DIRECTIONS[0]),
    )

    # Locate the ATMForce + its global-parameter handles.
    atm = next(system.getForce(i) for i in range(system.getNumForces())
               if isinstance(system.getForce(i), mm.ATMForce))

    # NVT only (no barostat). MonteCarloBarostat interacts badly with
    # ATMForce's displaced-particle u1 evaluation: the displaced shadow
    # ligand at u1 has huge effective energy (overlaps mid-box) so the
    # MC accept/reject step computes wrong volumes and triggers NaN.
    # ATM ABFE protocols use NVT on a pre-equilibrated box (the prep PDB
    # already comes from a fully relaxed NPT MD run).
    #
    # Warmup integrator (1 fs default) — gentle initial equilibration
    # absorbs the post-minimization energy gradient at small timestep
    # before switching to the production timestep. Bound legs at 2fs
    # without warmup NaN'd at the first integrator step regardless of
    # minimization length when displacement >= 5 nm.
    warmup_integrator = mm.LangevinMiddleIntegrator(
        TEMP_K * unit.kelvin,
        1.0 / unit.picosecond,
        warmup_timestep_fs * unit.femtoseconds,
    )
    warmup_integrator.setRandomNumberSeed(seed)

    platform = mm.Platform.getPlatformByName(platform_name)
    sim = app.Simulation(modeller.topology, system, warmup_integrator, platform)
    sim.context.setPositions(modeller.positions)
    sim.context.setVelocitiesToTemperature(TEMP_K * unit.kelvin, seed)

    with open(log_path, "w") as logfh:
        def log(msg: str) -> None:
            ts = time.strftime("%Y-%m-%dT%H:%M:%S")
            logfh.write(f"[{ts}] {msg}\n")
            logfh.flush()

        log(f"endpoint={endpoint} leg={leg} replicate={replicate} seed={seed}")
        log(f"pdb_path={pdb_path}")
        log(f"n_atoms={modeller.topology.getNumAtoms()}  "
            f"n_displaced={len(displaced)}  "
            f"displacement_nm={displacement_nm}")

        log("Minimization start")
        sim.minimizeEnergy(maxIterations=minimize_iter)
        log("Minimization done")

        # Warmup at gentle timestep BEFORE switching to production timestep.
        # Sets the initial Lambda1/Lambda2/Direction to state 0 so the
        # ATMForce sees the same configuration the production loop will.
        sim.context.setParameter("Lambda1", LAMBDAS_1[0])
        sim.context.setParameter("Lambda2", LAMBDAS_2[0])
        try:
            sim.context.setParameter("Direction", float(DIRECTIONS[0]))
        except Exception:
            pass
        warmup_steps = int(round(warmup_ps * 1000.0 / warmup_timestep_fs))
        if warmup_steps > 0:
            log(f"Warmup {warmup_ps:.2f} ps at {warmup_timestep_fs:.1f} fs ({warmup_steps} step) ...")
            sim.step(warmup_steps)
            log("Warmup done")

        # Migrate to production integrator/Simulation, transferring state.
        state_pos_vel = sim.context.getState(getPositions=True, getVelocities=True)
        prod_integrator = mm.LangevinMiddleIntegrator(
            TEMP_K * unit.kelvin,
            1.0 / unit.picosecond,
            timestep_fs * unit.femtoseconds,
        )
        prod_integrator.setRandomNumberSeed(seed + 7)
        sim = app.Simulation(modeller.topology, system, prod_integrator, platform)
        sim.context.setPositions(state_pos_vel.getPositions())
        sim.context.setVelocities(state_pos_vel.getVelocities())
        sim.context.setParameter("Lambda1", LAMBDAS_1[0])
        sim.context.setParameter("Lambda2", LAMBDAS_2[0])
        try:
            sim.context.setParameter("Direction", float(DIRECTIONS[0]))
        except Exception:
            pass

        kj = unit.kilojoule_per_mole.conversion_factor_to(unit.kilocalorie_per_mole)
        steps_per_sample = int(round(sample_interval_ps * 1000.0 / timestep_fs))
        prod_steps = int(round(prod_ns_per_window * 1e6 / timestep_fs))
        equil_steps = int(round(equil_ps * 1000.0 / timestep_fs))
        samples_per_window = max(1, prod_steps // steps_per_sample)

        log(f"timestep_fs={timestep_fs} steps_per_sample={steps_per_sample} "
            f"samples_per_window={samples_per_window} prod_steps={prod_steps} "
            f"equil_steps={equil_steps}")

        out_paths: List[str] = []
        try:
            # Initial equilibration at λ=0, direction=+1 (first state).
            sim.context.setParameter("Lambda1", LAMBDAS_1[0])
            sim.context.setParameter("Lambda2", LAMBDAS_2[0])
            # Direction is set via ATMForce — the canonical global-parameter
            # name is "Direction" (ommreplica writes self.atmforce.Direction()).
            try:
                sim.context.setParameter("Direction", float(DIRECTIONS[0]))
            except Exception:
                # Some OpenMM builds use a different name; fall through.
                pass

            if resume_skip_equilibration:
                log("Equilibration SKIPPED (--resume-skip-equilibration)")
            else:
                log("Equilibration start")
                sim.step(equil_steps)
                log("Equilibration done")

            schedule_ids = (list(state_ids) if state_ids is not None
                            else list(range(N_STATES)))

            # Resume support: skip states < resume_from_state. We still need
            # to walk through them at brief equilibration so the context is
            # populated with the correct (Lambda1/2, Direction, w0) before
            # the resume state. Without the walk-through, the resumed state
            # would inherit State 0 parameters from equilibration.
            fast_fwd_steps = max(1, int(round(direction_flip_warmup_ps
                                              * 1000.0 / max(direction_flip_warmup_fs, 1e-9))))
            if resume_from_state > 0:
                log(f"Resume mode: fast-forward through states 0..{resume_from_state - 1} "
                    f"at {direction_flip_warmup_fs:.1f} fs ({fast_fwd_steps} step each), "
                    f"no samples written.")
                # Migrate to a 1fs integrator for the fast-forward to avoid
                # NaN at any direction flip within the skipped range.
                ff_state = sim.context.getState(getPositions=True,
                                                getVelocities=True)
                ff_integrator = mm.LangevinMiddleIntegrator(
                    TEMP_K * unit.kelvin, 1.0 / unit.picosecond,
                    direction_flip_warmup_fs * unit.femtoseconds,
                )
                ff_integrator.setRandomNumberSeed(seed + 13)
                sim = app.Simulation(modeller.topology, system, ff_integrator,
                                     platform)
                sim.context.setPositions(ff_state.getPositions())
                sim.context.setVelocities(ff_state.getVelocities())
                kcal_ff = unit.kilocalorie_per_mole
                kj_unit_ff = unit.kilojoule_per_mole
                kcal_per_kj_ff = (1.0 * kcal_ff).value_in_unit(kj_unit_ff)
                for ff_id in range(resume_from_state):
                    sim.context.setParameter("Lambda1", LAMBDAS_1[ff_id])
                    sim.context.setParameter("Lambda2", LAMBDAS_2[ff_id])
                    try:
                        sim.context.setParameter("Direction",
                                                 float(DIRECTIONS[ff_id]))
                        sim.context.setParameter(
                            "Alpha",
                            ALPHA[ff_id] / kcal_per_kj_ff,
                        )
                        sim.context.setParameter(
                            "Uh", U0[ff_id] * kcal_per_kj_ff
                        )
                        sim.context.setParameter(
                            "W0", W0[ff_id] * kcal_per_kj_ff
                        )
                    except Exception as exc:
                        log(f"WARN: fast-fwd parameter set fallback at "
                            f"State {ff_id} ({exc})")
                    sim.step(fast_fwd_steps)
                log(f"Fast-forward done; resuming production at State "
                    f"{resume_from_state}")

            prev_direction = (float(DIRECTIONS[resume_from_state - 1])
                              if resume_from_state > 0 else None)
            # in_post_flip_region tracks whether we have entered the d=−1
            # half of the symmetric schedule. Once True, the production
            # integrator uses prod_timestep_fs_after_flip (when set) for
            # the remainder of the leg. Initialized False; flipped True at
            # the first direction-change event OR when resuming directly
            # into a d=−1 state.
            in_post_flip_region = False

            for stateid in schedule_ids:
                if stateid < resume_from_state:
                    # Already walked through during fast-forward; skip
                    # production sampling and proceed.
                    continue
                l1 = LAMBDAS_1[stateid]
                l2 = LAMBDAS_2[stateid]
                d = float(DIRECTIONS[stateid])
                a_kcal = ALPHA[stateid]
                uh_kcal = U0[stateid]
                w0_kcal = W0[stateid]

                # Convert to OpenMM-internal units (alpha is in 1/kJ when set
                # via ATMForce; ATMForce parameters use kJ/mol units).
                kcal = unit.kilocalorie_per_mole
                kj_unit = unit.kilojoule_per_mole
                kcal_per_kj = (1.0 * kcal).value_in_unit(kj_unit)
                alpha_per_kj = (a_kcal / kcal_per_kj)
                uh_kj = uh_kcal * kcal_per_kj
                w0_kj = w0_kcal * kcal_per_kj

                sim.context.setParameter("Lambda1", l1)
                sim.context.setParameter("Lambda2", l2)
                try:
                    sim.context.setParameter("Direction", d)
                    sim.context.setParameter("Alpha", alpha_per_kj)
                    sim.context.setParameter("Uh", uh_kj)
                    sim.context.setParameter("W0", w0_kj)
                except Exception as exc:
                    log(f"WARN: parameter set fallback ({exc})")

                # Direction-flip warmup: when the ATMForce Direction parameter
                # changes sign between adjacent states (the +1→−1 mid-point
                # transition of the symmetric 22-state schedule, States 10→11),
                # the potential gradient inverts. A 2fs production timestep
                # cannot absorb that step-function and the integrator NaN's at
                # the first step (observed 2026-05-30 r11 crash, empty .out).
                # Migrate transiently to a 1fs integrator, run
                # direction_flip_warmup_ps at the *new* state's parameters
                # (already set above), then migrate back to the production
                # integrator. State (positions/velocities) is preserved across
                # the migration.
                # Determine the production timestep to use for THIS state.
                # If a flip is happening now (or already happened in a prior
                # state of this leg), and prod_timestep_fs_after_flip is set,
                # use the smaller timestep for the rest of the schedule.
                # Otherwise fall back to the global production timestep.
                flip_now = (prev_direction is not None
                            and d != prev_direction)
                if flip_now:
                    in_post_flip_region = True

                this_state_prod_fs = (
                    prod_timestep_fs_after_flip
                    if (in_post_flip_region
                        and prod_timestep_fs_after_flip is not None)
                    else timestep_fs
                )

                if (prev_direction is not None and d != prev_direction
                        and direction_flip_warmup_ps > 0):
                    flip_steps = max(1, int(round(
                        direction_flip_warmup_ps * 1000.0
                        / max(direction_flip_warmup_fs, 1e-9)
                    )))
                    log(f"  Direction flip {prev_direction:+.0f} → {d:+.0f}: "
                        f"warmup {direction_flip_warmup_ps:.2f} ps at "
                        f"{direction_flip_warmup_fs:.1f} fs "
                        f"({flip_steps} step) before sampling. "
                        f"Post-flip prod timestep = {this_state_prod_fs:.1f} fs.")
                    flip_state = sim.context.getState(getPositions=True,
                                                     getVelocities=True)
                    flip_integrator = mm.LangevinMiddleIntegrator(
                        TEMP_K * unit.kelvin,
                        1.0 / unit.picosecond,
                        direction_flip_warmup_fs * unit.femtoseconds,
                    )
                    flip_integrator.setRandomNumberSeed(seed + 23 + stateid)
                    sim = app.Simulation(modeller.topology, system,
                                         flip_integrator, platform)
                    sim.context.setPositions(flip_state.getPositions())
                    sim.context.setVelocities(flip_state.getVelocities())
                    sim.context.setParameter("Lambda1", l1)
                    sim.context.setParameter("Lambda2", l2)
                    try:
                        sim.context.setParameter("Direction", d)
                        sim.context.setParameter("Alpha", alpha_per_kj)
                        sim.context.setParameter("Uh", uh_kj)
                        sim.context.setParameter("W0", w0_kj)
                    except Exception:
                        pass
                    sim.step(flip_steps)
                    # Migrate to the post-flip production integrator. If
                    # prod_timestep_fs_after_flip is set, this is typically
                    # smaller than the pre-flip 2fs (default 1fs for stability
                    # in the d=−1 half of the symmetric schedule).
                    post_state = sim.context.getState(getPositions=True,
                                                     getVelocities=True)
                    prod_integrator2 = mm.LangevinMiddleIntegrator(
                        TEMP_K * unit.kelvin,
                        1.0 / unit.picosecond,
                        this_state_prod_fs * unit.femtoseconds,
                    )
                    prod_integrator2.setRandomNumberSeed(seed + 41 + stateid)
                    sim = app.Simulation(modeller.topology, system,
                                         prod_integrator2, platform)
                    sim.context.setPositions(post_state.getPositions())
                    sim.context.setVelocities(post_state.getVelocities())
                    sim.context.setParameter("Lambda1", l1)
                    sim.context.setParameter("Lambda2", l2)
                    try:
                        sim.context.setParameter("Direction", d)
                        sim.context.setParameter("Alpha", alpha_per_kj)
                        sim.context.setParameter("Uh", uh_kj)
                        sim.context.setParameter("W0", w0_kj)
                    except Exception:
                        pass
                    # Re-locate ATMForce on the new context's system (same
                    # System object — ATMForce reference still valid).

                # On resume into a d=−1 half (no flip event but we entered
                # mid-schedule), still apply the after-flip timestep override
                # by migrating to the small-step integrator without a warmup.
                elif (prev_direction is None and d < 0
                      and prod_timestep_fs_after_flip is not None
                      and abs(this_state_prod_fs - timestep_fs) > 1e-9):
                    log(f"  Resume into d={d:+.0f} half: switching production "
                        f"timestep {timestep_fs:.1f} fs → "
                        f"{this_state_prod_fs:.1f} fs (stability override).")
                    post_state = sim.context.getState(getPositions=True,
                                                     getVelocities=True)
                    prod_integrator_resume = mm.LangevinMiddleIntegrator(
                        TEMP_K * unit.kelvin,
                        1.0 / unit.picosecond,
                        this_state_prod_fs * unit.femtoseconds,
                    )
                    prod_integrator_resume.setRandomNumberSeed(seed + 41 + stateid)
                    sim = app.Simulation(modeller.topology, system,
                                         prod_integrator_resume, platform)
                    sim.context.setPositions(post_state.getPositions())
                    sim.context.setVelocities(post_state.getVelocities())
                    sim.context.setParameter("Lambda1", l1)
                    sim.context.setParameter("Lambda2", l2)
                    try:
                        sim.context.setParameter("Direction", d)
                        sim.context.setParameter("Alpha", alpha_per_kj)
                        sim.context.setParameter("Uh", uh_kj)
                        sim.context.setParameter("W0", w0_kj)
                    except Exception:
                        pass
                    in_post_flip_region = True

                prev_direction = d

                log(f"State {stateid:02d}: λ1={l1:.3f} λ2={l2:.3f} d={d:+.0f} "
                    f"w0={w0_kcal:.2f} α={a_kcal:.3f}")

                # Per-state output file (AToM canonical layout). Replicates
                # append to the same file so calculate_uwham sees all samples
                # from every replicate at this state.
                state_dir = os.path.join(out_dir, f"r{stateid}")
                os.makedirs(state_dir, exist_ok=True)
                state_path = os.path.join(state_dir, f"{jobname}.out")
                state_fh = open(state_path, "a")
                out_paths.append(state_path)

                # Per-state steps_per_sample: when the timestep was reduced
                # via prod_timestep_fs_after_flip, we need more steps to cover
                # the same sample_interval_ps. samples_per_window stays the
                # same so each state contributes the same number of samples
                # to uwham regardless of integrator timestep.
                this_steps_per_sample = int(round(
                    sample_interval_ps * 1000.0
                    / max(this_state_prod_fs, 1e-9)
                ))

                # Drain the per-window samples: equilibrate-then-sample loop.
                # We do not re-equilibrate between adjacent states (warm
                # crossover is the standard ATM practice).
                for sample_id in range(samples_per_window):
                    sim.step(this_steps_per_sample)
                    state = sim.context.getState(getEnergy=True)
                    pot_kj = state.getPotentialEnergy().value_in_unit(kj_unit)
                    pert_pair = atm.getPerturbationEnergy(sim.context)
                    u1_kj = (pert_pair[0].value_in_unit(kj_unit)
                             if hasattr(pert_pair[0], "value_in_unit")
                             else float(pert_pair[0]))
                    u0_kj_sample = (pert_pair[1].value_in_unit(kj_unit)
                                    if hasattr(pert_pair[1], "value_in_unit")
                                    else float(pert_pair[1]))
                    pert_kj_raw = u1_kj - u0_kj_sample
                    # Apply soft-core cap (matches ommworker.py:415-417 +
                    # AtomUtils.softCorePertE) so uwham's _bias_fcn sees a
                    # bounded pertE. Without this, displaced overlaps produce
                    # 1e15-scale raw pertE that infinitises the softplus.
                    umax_kj = UMAX_KCAL * kcal_per_kj
                    ub_kj = UBCORE_KCAL * kcal_per_kj
                    pert_kj = soft_core_pert_e(
                        pert_kj_raw, umax_kj, ub_kj, ACORE
                    )
                    bias_kj = softplus_bias_kj(
                        l1, l2, alpha_per_kj, uh_kj, w0_kj, pert_kj
                    )

                    write_out_line(
                        state_fh,
                        stateid=stateid,
                        temperature_K=TEMP_K,
                        direction=d,
                        lambda1=l1, lambda2=l2,
                        alpha_per_kj=alpha_per_kj,
                        uh_kj=uh_kj, w0_kj=w0_kj,
                        pot_energy_kj=pot_kj,
                        pert_energy_kj=pert_kj,
                        bias_energy_kj=bias_kj,
                    )
                    state_fh.flush()

                    if progress_cb is not None:
                        progress_cb(stateid, sample_id, samples_per_window)
                state_fh.close()
        finally:
            pass

        log(f"Done in {time.time() - t_start:.1f} s")

    return {
        "endpoint": endpoint,
        "leg": leg,
        "replicate": replicate,
        "seed": seed,
        "pdb_path": pdb_path,
        "out_paths": out_paths,
        "log_path": log_path,
        "n_atoms_total": modeller.topology.getNumAtoms(),
        "n_displaced": len(displaced),
        "displacement_nm": list(displacement_nm),
        "runtime_s": time.time() - t_start,
        "n_states": N_STATES,
        "samples_per_window": samples_per_window,
        "n_total_samples": samples_per_window * N_STATES,
        "prod_ns_per_window": prod_ns_per_window,
        "equil_ps": equil_ps,
        "sample_interval_ps": sample_interval_ps,
        "timestep_fs": timestep_fs,
    }


# ---------------------------------------------------------------------------
# uwham analysis driver (per leg).
# ---------------------------------------------------------------------------
def analyze_leg_with_uwham(
    rundir: str, jobname: str = "trackb",
    mintimeid: Optional[int] = None, maxtimeid: Optional[int] = None,
) -> Dict[str, Any]:
    """Run ``atom_openmm.uwham.calculate_uwham`` on one leg directory.

    rundir must contain ``r0/{jobname}.out`` … ``r{N-1}/{jobname}.out`` — the
    canonical AToM layout.

    The free energies returned by uwham use the symmetric 22-state schedule
    defaults (lambda1/lambda2 = LAMBDA_FWD+LAMBDA_BWD, intermd as above).
    """
    from atom_openmm import uwham
    out = uwham.calculate_uwham(
        rundir=rundir, jobname=jobname,
        mintimeid=mintimeid, maxtimeid=maxtimeid,
        intermd=INTERMD,
        lambda1=LAMBDAS_1, lambda2=LAMBDAS_2,
        alpha=ALPHA, u0=U0, w0=W0,
    )
    return out


# ---------------------------------------------------------------------------
# Top-level orchestration.
# ---------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Track B production launcher — ATM/ATS alchemical "
            "MTR↔Trp ΔΔG_bind (per-endpoint ABFE → paired)."
        )
    )
    p.add_argument("--seed-tag", default="s7",
                   help="prepared seed under outputs/2QKI_{Cp4_hybrid,WT}_calib_<seed>")
    p.add_argument("--out-root", default="outputs/_trackb/production",
                   help="output root under project")
    p.add_argument("--endpoints", default="cp4,wt",
                   help="comma-list of endpoints to run (cp4 / wt)")
    p.add_argument("--legs", default="bound,free",
                   help="comma-list of legs to run (bound / free)")
    p.add_argument("--n-replicates", type=int, default=3,
                   help="independent velocity-seed replicates per leg")
    p.add_argument("--prod-ns", type=float, default=5.0,
                   help="production ns per λ window")
    p.add_argument("--equil-ps", type=float, default=500.0,
                   help="equilibration ps per replicate (before state 0)")
    p.add_argument("--sample-ps", type=float, default=50.0,
                   help="sample interval ps")
    p.add_argument("--minimize-iter", type=int, default=5000,
                   help="energy-minimization max iterations")
    p.add_argument("--timestep-fs", type=float, default=2.0,
                   help="MD timestep (fs)")
    p.add_argument("--displacement-nm", default="2.5,0.0,0.0",
                   help="ABFE displacement vector (nm)")
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--cuda-device", default="0",
                   help="forced CUDA_VISIBLE_DEVICES (host 5070Ti = 0)")
    p.add_argument("--seed-list",
                   help="comma-list of N_replicate*N_endpoint*N_leg ints "
                        "(default: deterministic ascending from --base-seed)")
    p.add_argument("--base-seed", type=int, default=20260530,
                   help="base seed for deterministic seed-list generation")
    p.add_argument("--smoke-mode", action="store_true",
                   help="smoke pre-flight: 1 replicate, 2 states (0, midpoint), "
                        "100 step production, 50 step equil — gates launch")
    p.add_argument("--analyze-only", action="store_true",
                   help="skip MD; rerun uwham analysis on existing .out files")
    p.add_argument("--jobname", default="trackb",
                   help="basename of per-replicate .out / .log files")
    p.add_argument("--resume-from-state", type=int, default=0,
                   help="Resume MD at this state id (0..21). States below are "
                        "fast-forwarded at 1 fs with no samples written so the "
                        "context arrives at the resume state with the correct "
                        "Lambda/Direction. Use to restart after NaN crashes. "
                        "Default 0 = full schedule.")
    p.add_argument("--resume-skip-equilibration", action="store_true",
                   help="Skip the State-0 equilibration entirely on a resume "
                        "(use when the system was already well-equilibrated "
                        "in the prior run; saves ~6 min per replicate).")
    p.add_argument("--direction-flip-warmup-ps", type=float, default=0.1,
                   help="Per-state warmup (ps) at small timestep when ATMForce "
                        "Direction parameter flips sign between adjacent states "
                        "(the +1→−1 midpoint transition of the 22-state "
                        "symmetric schedule). Absorbs the potential-gradient "
                        "discontinuity that would otherwise NaN a 2fs step. "
                        "0 disables (legacy v1 behavior, NaN risk at State 11).")
    p.add_argument("--direction-flip-warmup-fs", type=float, default=1.0,
                   help="Timestep (fs) used during the direction-flip warmup.")
    p.add_argument("--prod-timestep-fs-after-flip", type=float, default=None,
                   help="Production timestep (fs) for States AFTER the first "
                        "Direction-parameter sign flip. When None, falls back "
                        "to --timestep-fs (legacy v1 behavior, 2fs throughout). "
                        "Recommended 1.0 for the d=−1 half of the 22-state "
                        "symmetric schedule when 2fs NaNs immediately after "
                        "the warmup (as observed at State 11 on 2026-05-30 "
                        "even with 0.1ps direction-flip warmup engaged). 2× "
                        "wall-time cost on the second half; correctness >> speed.")
    p.add_argument("--resume-endpoints", default=None,
                   help="Comma-list of endpoints to *restrict* the resume run to "
                        "(e.g. 'cp4' if only cp4/bound was running when the crash "
                        "happened). Default None = same as --endpoints.")
    p.add_argument("--resume-legs", default=None,
                   help="Comma-list of legs to *restrict* the resume run to. "
                        "Default None = same as --legs.")
    args = p.parse_args()

    # ---- CUDA device pin (hardware contract) -----------------------
    os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda_device

    # ---- Pre-flight: charge axis (Q1 gate) -----------------------------
    axis = verify_charge_axis()
    if not axis["all_residues_within_tol"]:
        print(f"ERROR: charge axis FAIL (max|Σq|={axis['max_abs_sigma_q_e']:.2e})",
              file=sys.stderr)
        return 2

    # ---- Resolve leg inputs --------------------------------------------
    inputs = resolve_leg_inputs(args.seed_tag)
    cp4_bound_pdb = inputs["bound"]["cp4"]
    wt_bound_pdb = inputs["bound"]["wt"]
    cp4_final = inputs["final"]["cp4"]
    wt_final = inputs["final"]["wt"]
    hydrogens_xml = inputs.get("hydrogens_xml")

    if not os.path.isfile(cp4_bound_pdb):
        print(f"ERROR: missing Cp4 bound PDB {cp4_bound_pdb}", file=sys.stderr)
        return 2
    if not os.path.isfile(wt_bound_pdb):
        print(f"ERROR: missing WT bound PDB {wt_bound_pdb}", file=sys.stderr)
        return 2
    if not (cp4_final and os.path.isfile(cp4_final)):
        print(f"ERROR: missing Cp4 final PDB {cp4_final}", file=sys.stderr)
        return 2
    if not (wt_final and os.path.isfile(wt_final)):
        print(f"ERROR: missing WT final PDB {wt_final}", file=sys.stderr)
        return 2

    # ---- Prepare free-leg PDBs (one per endpoint) ----------------------
    out_root = os.path.join(_PROJ_ROOT, args.out_root)
    os.makedirs(out_root, exist_ok=True)
    free_work = os.path.join(out_root, "_free_pdbs")
    os.makedirs(free_work, exist_ok=True)

    cp4_free_pdb = os.path.join(free_work, "cp4_freeleg.pdb")
    wt_free_pdb = os.path.join(free_work, "wt_freeleg.pdb")
    if not os.path.isfile(cp4_free_pdb):
        prepare_free_peptide_from_final(cp4_final, cp4_free_pdb)
    if not os.path.isfile(wt_free_pdb):
        prepare_free_peptide_from_final(wt_final, wt_free_pdb)

    endpoint_pdb = {
        ("cp4", "bound"): cp4_bound_pdb,
        ("cp4", "free"):  cp4_free_pdb,
        ("wt",  "bound"): wt_bound_pdb,
        ("wt",  "free"):  wt_free_pdb,
    }

    endpoints = [e.strip() for e in args.endpoints.split(",") if e.strip()]
    legs = [l.strip() for l in args.legs.split(",") if l.strip()]
    for e in endpoints:
        assert e in ("cp4", "wt"), f"unknown endpoint {e}"
    for l in legs:
        assert l in ("bound", "free"), f"unknown leg {l}"

    # --resume-endpoints / --resume-legs restrict the MD execution loop only.
    # The full endpoints / legs set is still used for downstream uwham analysis
    # so existing per-leg .out files remain consumable; if a leg is missing the
    # analyze block will record an error per pre-existing logic.
    if args.resume_endpoints:
        md_endpoints = [e.strip() for e in args.resume_endpoints.split(",")
                        if e.strip()]
        for e in md_endpoints:
            assert e in endpoints, (
                f"--resume-endpoints {e} not in --endpoints {endpoints}"
            )
    else:
        md_endpoints = list(endpoints)
    if args.resume_legs:
        md_legs = [l.strip() for l in args.resume_legs.split(",") if l.strip()]
        for l in md_legs:
            assert l in legs, (
                f"--resume-legs {l} not in --legs {legs}"
            )
    else:
        md_legs = list(legs)

    # ---- Deterministic seed list ---------------------------------------
    n_units = len(endpoints) * len(legs) * args.n_replicates
    if args.seed_list:
        seeds = [int(x) for x in args.seed_list.split(",")]
        if len(seeds) != n_units:
            print(f"ERROR: --seed-list requires {n_units} ints (got {len(seeds)})",
                  file=sys.stderr)
            return 2
    else:
        seeds = [args.base_seed + i for i in range(n_units)]

    displacement_nm = tuple(float(x) for x in args.displacement_nm.split(","))
    assert len(displacement_nm) == 3, "--displacement-nm must be x,y,z"

    # ---- Smoke vs production knobs -------------------------------------
    if args.smoke_mode:
        n_replicates_eff = 1
        prod_ns_eff = 0.001  # ~1 ps per window
        equil_ps_eff = 0.5
        sample_ps_eff = 0.1
        # CLI --minimize-iter overrides; smoke default 2000 (must be enough
        # to relax bound leg ~123k atoms — bound 2fs requires it).
        minimize_iter_eff = (args.minimize_iter
                             if args.minimize_iter != 5000 else 2000)
        warmup_ps_eff = 0.1  # 0.1 ps × 1fs = 100 step warmup
        # Replace 22-state schedule with 2-state (head + midpoint) for smoke.
        smoke_state_ids = [0, 10]  # λ=0 fwd, λ=0.5 fwd-midpoint
    else:
        n_replicates_eff = args.n_replicates
        prod_ns_eff = args.prod_ns
        equil_ps_eff = args.equil_ps
        sample_ps_eff = args.sample_ps
        minimize_iter_eff = args.minimize_iter
        warmup_ps_eff = 1.0  # 1 ps × 1fs = 1000 step warmup
        smoke_state_ids = None

    # Recompute n_units with smoke replicate count.
    if args.smoke_mode:
        n_units = len(endpoints) * len(legs) * n_replicates_eff
        seeds = seeds[:n_units]

    # ---- Build verbatim launch command + env snapshot ------------------
    launch_command = "CUDA_VISIBLE_DEVICES=" + shlex.quote(args.cuda_device) \
        + " /home/san/miniconda3/envs/atm/bin/python " \
        + " ".join(shlex.quote(a) for a in sys.argv)

    try:
        _openmm_ver = mm.version.full_version
    except Exception:
        _openmm_ver = mm.Platform.getOpenMMVersion()
    env_snapshot = {
        "python": sys.executable,
        "openmm": _openmm_ver,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "cwd": os.getcwd(),
        "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
        "atom_openmm_path": None,
    }
    try:
        import atom_openmm as _ao
        env_snapshot["atom_openmm_path"] = _ao.__file__
    except ImportError:
        pass

    pre_register = {
        "schema": "ranking_only_v1",
        "magotti_ssot_direction": "MTR_less_favorable_than_WT (DDG > 0)",
        "magotti_kcal_mol": "-1.4 (ITC) / -3.0 (SPR) — NOT for absolute comparison (ranking-only)",
        "outcomes": {
            "1_sign_stable": {
                "criterion": "DDG > 0 AND sigma_btwn <= 0.5",
                "interpretation": "primary ranking-positive outcome",
            },
            "2_sigma_drift": {
                "criterion": "sigma_btwn > 0.7 (>140% expected upper bound)",
                "interpretation": "replicate count expansion before ranking claim",
            },
            "3_magnitude_drift": {
                "criterion": "|DDG| outside [0.5, 5.0]",
                "interpretation": "methodology sensitivity review (cyclic_ss FF bias check)",
            },
            "4_sign_flip": {
                "criterion": "DDG < 0",
                "interpretation": "escalate (charge regime / sampling regime re-verification)",
            },
        },
        "regime": "ranking_only",
    }

    run_metadata = {
        "track": "B",
        "method": "ATM ABFE per-endpoint with paired DDG_bind",
        "seed_tag": args.seed_tag,
        "endpoints": endpoints,
        "legs": legs,
        "md_endpoints_this_invocation": md_endpoints,
        "md_legs_this_invocation": md_legs,
        "resume_from_state": int(args.resume_from_state),
        "resume_skip_equilibration": bool(args.resume_skip_equilibration),
        "direction_flip_warmup_ps": float(args.direction_flip_warmup_ps),
        "direction_flip_warmup_fs": float(args.direction_flip_warmup_fs),
        "prod_timestep_fs_after_flip": (
            float(args.prod_timestep_fs_after_flip)
            if args.prod_timestep_fs_after_flip is not None else None
        ),
        "n_replicates": n_replicates_eff,
        "schedule": {
            "n_states": N_STATES,
            "lambda1": LAMBDAS_1,
            "lambda2": LAMBDAS_2,
            "directions": DIRECTIONS,
            "intermd": INTERMD,
            "w0_kcal": W0,
            "alpha_per_kcal": ALPHA,
            "u0_kcal": U0,
            "temperature_K": TEMP_K,
        },
        "prod_ns_per_window": prod_ns_eff,
        "equil_ps_per_replicate": equil_ps_eff,
        "sample_interval_ps": sample_ps_eff,
        "minimize_iter": minimize_iter_eff,
        "timestep_fs": args.timestep_fs,
        "displacement_nm": list(displacement_nm),
        "displacement_target": "binder_chain (ABFE — entire ligand decoupled)",
        "platform": args.platform,
        "cuda_device": args.cuda_device,
        "smoke_mode": bool(args.smoke_mode),
        "smoke_state_ids": smoke_state_ids,
        "charge_axis": axis,
        "seeds": seeds,
        "launch_command": launch_command,
        "env_snapshot": env_snapshot,
        "pre_register": pre_register,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "regime": "ranking_only",
        "method_ref": "production blockers cleared, 2026-05-30",
    }

    # On a resume (resume_from_state>0 or resume_endpoints set), preserve the
    # original v1 run_metadata.json (it carries the original launch_command and
    # full schedule reference). Resume metadata goes to a timestamped sidecar.
    is_resume_invocation = (
        args.resume_from_state > 0
        or bool(args.resume_endpoints)
        or bool(args.resume_legs)
        or args.resume_skip_equilibration
    )
    base_meta = os.path.join(out_root, "run_metadata.json")
    if is_resume_invocation and os.path.isfile(base_meta):
        ts_tag = time.strftime("%Y%m%dT%H%M%S")
        meta_path = os.path.join(out_root, f"run_metadata_resume_{ts_tag}.json")
    else:
        meta_path = base_meta
    with open(meta_path, "w") as fh:
        json.dump(run_metadata, fh, indent=2)
    print(f"Run metadata: {meta_path}"
          + (" (resume sidecar)" if meta_path != base_meta else ""))
    print(f"Launch command (verbatim): {launch_command}")
    print()

    # ---- Execute per (endpoint, leg, replicate) ------------------------
    results: List[Dict[str, Any]] = []
    if not args.analyze_only:
        # Resume-aware seed indexing: keep the deterministic seed mapping the
        # same as the original launch (endpoint × leg × replicate ordering),
        # but only execute (md_endpoints × md_legs) pairs. This guarantees
        # that resumed cp4/bound r0 gets the SAME seed it had in the v1 run.
        seed_idx = {}
        for i, e in enumerate(endpoints):
            for j, l in enumerate(legs):
                for r in range(n_replicates_eff):
                    seed_idx[(e, l, r)] = (
                        i * len(legs) * n_replicates_eff
                        + j * n_replicates_eff
                        + r
                    )
        for endpoint in md_endpoints:
            for leg in md_legs:
                leg_dir = os.path.join(out_root, endpoint, leg)
                os.makedirs(leg_dir, exist_ok=True)
                pdb_path = endpoint_pdb[(endpoint, leg)]
                for replicate in range(n_replicates_eff):
                    seed = seeds[seed_idx[(endpoint, leg, replicate)]]
                    print(f"\n=== {endpoint} / {leg} / r{replicate} "
                          f"(seed={seed}) ===")
                    res = run_one_replicate(
                        pdb_path=pdb_path,
                        endpoint=endpoint,
                        leg=leg,
                        replicate=replicate,
                        seed=seed,
                        out_dir=leg_dir,
                        equil_ps=equil_ps_eff,
                        prod_ns_per_window=prod_ns_eff,
                        sample_interval_ps=sample_ps_eff,
                        minimize_iter=minimize_iter_eff,
                        timestep_fs=args.timestep_fs,
                        platform_name=args.platform,
                        displacement_nm=displacement_nm,
                        hydrogens_xml=hydrogens_xml,
                        solvate=(leg == "bound"),
                        jobname=args.jobname,
                        state_ids=smoke_state_ids,
                        warmup_ps=warmup_ps_eff,
                        direction_flip_warmup_ps=args.direction_flip_warmup_ps,
                        direction_flip_warmup_fs=args.direction_flip_warmup_fs,
                        prod_timestep_fs_after_flip=args.prod_timestep_fs_after_flip,
                        resume_from_state=args.resume_from_state,
                        resume_skip_equilibration=args.resume_skip_equilibration,
                    )
                    results.append(res)

    # ---- Analyze: per-leg uwham → per-endpoint ΔG_bind → paired ΔΔG ----
    leg_dgs: Dict[Tuple[str, str], Dict[str, Any]] = {}
    analysis_errors: List[str] = []
    if args.smoke_mode:
        print("\n[smoke-mode] skipping uwham analysis (incomplete schedule).")
        # Inspect each replicate's per-state .out files for finiteness.
        smoke_finite_all = True
        smoke_summary = []
        for r in results:
            for p in r["out_paths"]:
                if not os.path.isfile(p):
                    smoke_finite_all = False
                    smoke_summary.append(f"MISSING {p}")
                    continue
                with open(p) as fh:
                    lines = fh.readlines()
                if not lines:
                    smoke_finite_all = False
                    smoke_summary.append(f"EMPTY {p}")
                    continue
                # Column 8 (0-indexed) = potE, 9 = pertE, 10 = bias
                bad = False
                for line in lines:
                    parts = line.split()
                    if len(parts) < 11:
                        bad = True; break
                    try:
                        pot, pert, bias = float(parts[8]), float(parts[9]), float(parts[10])
                    except ValueError:
                        bad = True; break
                    if not (np.isfinite(pot) and np.isfinite(pert) and np.isfinite(bias)):
                        bad = True; break
                if bad:
                    smoke_finite_all = False
                    smoke_summary.append(f"NONFINITE {p}")
                else:
                    smoke_summary.append(f"OK ({len(lines)} samples) {p}")
        print(f"[smoke-mode] finite check: {'PASS' if smoke_finite_all else 'FAIL'}")
        for s in smoke_summary:
            print(f"  {s}")
    for endpoint in endpoints:
        for leg in legs:
            if args.smoke_mode:
                continue
            leg_dir = os.path.join(out_root, endpoint, leg)
            try:
                out = analyze_leg_with_uwham(leg_dir, jobname=args.jobname)
                # uwham returns dict with "ze" (free energies, last entry is
                # the total ΔG between baseline and the target state). The
                # "result" sub-dict carries per-state ze + ve. Report the
                # global ΔG and SE at the highest λ.
                ze = np.asarray(out["ze"])
                ve = np.asarray(out["ve"])
                # ΔG between state[-1] and state[0] (kcal/mol) — uwham
                # uses ze in reduced units (β U), so divide by β at 300K.
                beta = 1.0 / (0.001986209 * TEMP_K)  # 1/(kcal/mol)
                # uwham returns ze in kT units; convert to kcal/mol.
                dg_kcal = float((ze[-1] - ze[0]) / beta)
                se_kcal = float(np.sqrt(ve[-1] + ve[0]) / beta)
                leg_dgs[(endpoint, leg)] = {
                    "endpoint": endpoint,
                    "leg": leg,
                    "dg_kcal": dg_kcal,
                    "se_kcal": se_kcal,
                    "ze": ze.tolist(),
                    "ve": ve.tolist(),
                }
                print(f"uwham {endpoint}/{leg}: ΔG = {dg_kcal:+.3f} ± "
                      f"{se_kcal:.3f} kcal/mol")
            except Exception as exc:
                msg = f"uwham {endpoint}/{leg} FAIL: {exc}"
                print(msg, file=sys.stderr)
                analysis_errors.append(msg)
                leg_dgs[(endpoint, leg)] = {
                    "endpoint": endpoint, "leg": leg,
                    "dg_kcal": None, "se_kcal": None,
                    "error": str(exc),
                }

    # ΔG_bind(endpoint) = ΔG_alch(bound) − ΔG_alch(free); requires both legs.
    endpoint_dg_bind: Dict[str, Dict[str, Any]] = {}
    for endpoint in endpoints:
        bound = leg_dgs.get((endpoint, "bound"))
        free = leg_dgs.get((endpoint, "free"))
        if (bound and free and bound.get("dg_kcal") is not None
                and free.get("dg_kcal") is not None):
            dg = bound["dg_kcal"] - free["dg_kcal"]
            se = float(np.sqrt(bound["se_kcal"]**2 + free["se_kcal"]**2))
            endpoint_dg_bind[endpoint] = {
                "dg_bind_kcal": dg, "se_kcal": se,
                "dg_alch_bound": bound["dg_kcal"],
                "dg_alch_free": free["dg_kcal"],
            }
        else:
            endpoint_dg_bind[endpoint] = {
                "dg_bind_kcal": None, "se_kcal": None,
                "error": "missing bound or free leg uwham",
            }

    # ΔΔG_bind = ΔG_bind(Cp4) − ΔG_bind(WT).
    cp4_dg = endpoint_dg_bind.get("cp4", {}).get("dg_bind_kcal")
    wt_dg = endpoint_dg_bind.get("wt", {}).get("dg_bind_kcal")
    ddg = None
    ddg_se = None
    if cp4_dg is not None and wt_dg is not None:
        cp4_se = endpoint_dg_bind["cp4"]["se_kcal"]
        wt_se = endpoint_dg_bind["wt"]["se_kcal"]
        ddg = cp4_dg - wt_dg
        ddg_se = float(np.sqrt(cp4_se**2 + wt_se**2))

    final = {
        "track": "B",
        "scope": "ATM ABFE per-endpoint, paired DDG",
        "ddg_bind_kcal": ddg,
        "ddg_se_kcal": ddg_se,
        "endpoint_dg_bind_kcal": endpoint_dg_bind,
        "leg_dg_alch_kcal": {
            f"{e}/{l}": v for (e, l), v in leg_dgs.items()
        },
        "replicates_run": [
            {
                "endpoint": r["endpoint"], "leg": r["leg"],
                "replicate": r["replicate"], "seed": r["seed"],
                "runtime_s": r["runtime_s"], "n_atoms": r["n_atoms_total"],
            } for r in results
        ],
        "n_replicates_per_leg": n_replicates_eff,
        "analysis_errors": analysis_errors,
        "metadata_path": meta_path,
        "regime": "ranking_only",
        "method_ref": "production blockers cleared, 2026-05-30",
        "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }

    final_path = os.path.join(out_root, "ddint_kcal_alchemical.json")
    with open(final_path, "w") as fh:
        json.dump(final, fh, indent=2)
    print(f"\nFinal result: {final_path}")
    if ddg is not None and ddg_se is not None:
        print(f"ΔΔG_bind (Cp4 − WT) = {ddg:+.3f} ± {ddg_se:.3f} kcal/mol "
              f"[regime: ranking_only]")
        # Pre-registered outcome match (informational only — ranking-only).
        if ddg > 0 and ddg_se <= 0.5:
            print("Outcome match: (1) sign-stable")
        elif ddg_se > 0.7:
            print("Outcome match: (2) sigma_drift")
        elif abs(ddg) < 0.5 or abs(ddg) > 5.0:
            print("Outcome match: (3) magnitude_drift")
        elif ddg < 0:
            print("Outcome match: (4) sign_flip — escalation required")
    elif not args.smoke_mode:
        print("ΔΔG_bind: NOT computed (analysis incomplete)")

    if args.smoke_mode:
        return 0 if smoke_finite_all else 1
    return 0 if ddg is not None else 1


if __name__ == "__main__":
    sys.exit(main())
