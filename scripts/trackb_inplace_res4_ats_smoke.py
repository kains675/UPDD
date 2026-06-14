#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B — in-place residue-4 ATS Tier-1 ASSERT smoke (v0.8 estimator de-risk).

Design rationale: C1-C8 coordinate-swap criteria + the Q5
endpoint-equivalence addition. This is a v0.8 PREDICTION test
(R-18) of whether the in-place residue-4 fused dual-topology construction
eliminates the whole-binder coverage NaN — NOT a converged DDG_bind, NOT a
campaign result, NOT a Magotti/absolute comparison (R-11 ranking-only).

The default ``--swap-mode genuine`` builds the REAL HE1<->methyl coupling
toggle: the methyl stays coupled (zero-offset NE1-attach) and HE1's nonbonded
coupling is removed at the u1 endpoint by translating it into bulk against a
massless dummy NE1 reference (its internal angle terms are decoupled so the
perturbation is the coupling change, NOT geometry strain). The resulting
|u1-u0| is a physically plausible ones-to-tens of kcal/mol at the TRUE mutation
geometry — NOT the legacy ``--swap-mode var_park`` parked clash (~1e4 kcal/mol,
retained for diagnostics only).

What this proves (Tier-1, 0-step energy decomposition at Lambda1=Lambda2=0.5,
Direction=+1, on the PME-solvated FREE leg) — every check is an ASSERT, not a
print (C6):
  (a) u0/u1/usc/E_ATM are all FINITE  -> the whole-binder coverage overflow
      (Object A) is defeated for the in-place ~4-5-atom perturbation;
  (b) u1 != u0                        -> the displacement-0 null op (Object B)
      is defeated; the swap is genuinely wired;
  (c) max per-atom force is bounded   -> finite-at-one-frame is not a near-
      singular term that NaNs under integration;
  (d) usc <= Umax (REAL soft-core cap recomputation) AND |u1-u0| within the
      physically-plausible bound for a HE1<->methyl edit -> the genuine swap is
      a real coupling change, not a parked-clash / bond-strain artifact;
  (e) ENDPOINT-EQUIVALENCE            -> the fused REFERENCE energy (u0, the
      un-displaced WT-coupled state, i.e. the TRUE mutation geometry — not the
      decoupled u1 geometry) reproduces the pure single-topology MTR-only energy
      within tolerance (the only cheap check that the box is physically
      equivalent, not merely finite — C6e / Q5).

What this does NOT prove (R-18, honest scope):
  - The eventual DDG_bind is correct (needs the full lambda ladder + overlap +
    convergence — future scope).
  - That common-charge continuity holds with the PRODUCTION (un-harmonized)
    hybrid MTR XML. By default this smoke runs WITHOUT charge harmonization, so
    it will surface the MC1 charge-discontinuity outcome (pre-registered S1
    outcome (ii)) — a REAL finding, not a pass. Pass --harmonize-common-charges
    to run the MECHANICAL validation path (isolates swap/finiteness/endpoint-
    equivalence from the charge gap; NOT production-valid for a quantitative
    DDG).

Pre-registered outcomes (S1):
  (i)   finite + u1!=u0 + endpoint-equiv pass -> proceed Tier-2
  (ii)  finite but endpoint-mismatch          -> MC1 charge-continuity bug
  (iii) finite but u1~u0                       -> swap not wired
  (iv)  NaN/overflow                           -> seed or common-atom-order bug

Runs openmm-only in the ``qmmm`` env (no atom_openmm needed — it lives only in
env ``atm``). The 5070Ti (CUDA) is the target platform; falls back to the
Reference platform if CUDA is unavailable (slower but identical asserts).

DOI references (anti-AI-trail clean; DOIs allowed):
  - Gallicchio 2021 J Chem Theory Comput, DOI 10.1021/acs.jctc.1c00753
  - Azimi et al. 2022 J Chem Inf Model 62(2):309, DOI 10.1021/acs.jcim.1c01129
  - Boresch 2003 (relative cycle); Mey et al. 2020 LiveCoMS
    DOI 10.33011/livecoms.2.1.18378
"""

import argparse
import json
import os
import sys

import numpy as np

import openmm as mm
import openmm.unit as unit

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)
_UTILS = os.path.join(_PROJ, "utils")
if _UTILS not in sys.path:
    sys.path.insert(0, _UTILS)

import atm_trackB_setup as ats  # noqa: E402


# Tier-1 assertion thresholds (C6 spec).
MAX_FORCE_KCAL_PER_MOL_A = 1.0e5      # (c) bounded max per-atom force
U1_NE_U0_MIN_KCAL = 1.0e-3            # (b) the perturbation must be non-trivial
# (d) |u1-u0| must be PHYSICALLY PLAUSIBLE for a HE1<->methyl side-chain edit:
# ones-to-tens of kcal/mol, NOT a parked-clash 1e4+. A genuine HE1 nonbonded
# decoupling is ~10 kcal/mol; we allow generous headroom but firmly reject the
# var_park parked-clash regime (~1.5e4 kcal/mol).
U1_MINUS_U0_PLAUSIBLE_MAX_KCAL = 1.0e3
ENDPOINT_EQUIV_TOL_KCAL = 10.0        # (e) u0 vs MTR-only energy tolerance

# Tier-2 short-dynamics stability thresholds (finite-UNDER-INTEGRATION). Tier-1
# proved finite-at-one-frame (0-step); Tier-2 proves the genuine fused box
# survives short MD at Lambda=0.5/Direction=+1 with a 1 fs UNCONSTRAINED
# alchemical-H integrator (no SHAKE on the appearing/disappearing methyl/HE1 H),
# no HMR. Per the deep-dig S2/S3 spec.
TIER2_DEFAULT_STEPS = 300             # 200-500 step window (S2)
TIER2_CHECKPOINT_EVERY = 50           # PE/KE/T sampled every 50 steps
TIER2_TEMPERATURE_K = 300.0           # Langevin bath temperature
TIER2_FRICTION_PER_PS = 1.0           # 1/ps friction
# Sane temperature band. The UPPER edge is the load-bearing RUNAWAY guard (an
# exploding box co-occurs with T -> inf, Fmax -> inf, NaN). The LOWER edge guards
# against a frozen/sunk box but is applied ONLY after an equilibration window:
# with 1/ps friction the thermalization timescale is ~1 ps (~1000 steps at 1 fs),
# so over a 200-500 step window a small / unsolvated box legitimately sits below
# 300 K (configurational PE released during relaxation has not re-thermalized).
# Judging instantaneous T before equilibration is not physically meaningful, so
# the first TIER2_T_EQUIL_SKIP checkpoints get a finiteness+upper-bound check
# only. The PME-SOLVATED box (thousands of water DOF) thermalizes near 300 K.
TIER2_T_BAND_K = (100.0, 500.0)       # (lower applied post-equilibration only)
TIER2_T_EQUIL_SKIP = 2                # checkpoints to skip for the LOWER bound
TIER2_MAX_FORCE_KCAL_PER_MOL_A = 1.0e5  # bounded max per-atom force throughout
# S3 integrator matrix: the R3 choice (1 fs, UNCONSTRAINED) is the PRIMARY and
# MUST pass. The other cells are EVIDENCE — it is fine/expected if a 2 fs step on
# an appearing/disappearing H is less stable. Short window (S3 default 150).
TIER2_MATRIX_STEPS = 150
TIER2_INTEGRATOR_MATRIX = (
    # (label, dt_fs, constraints) -- primary first.
    ("dt1fs_unconstrained", 1.0, "none"),
    ("dt2fs_unconstrained", 2.0, "none"),
    ("dt2fs_hbonds", 2.0, "hbonds"),
)


def _softcore_usc(u_kcal, umax_kcal, ubcore_kcal, acore):
    """Reproduce the upstream ATM soft-core cap usc(u) (ommsystem.py:783-789).

    For |u| <= Ubcore the perturbation passes through unchanged; above Ubcore it
    is compressed toward Umax via the fsc softplus. Used to make the C6d cap
    check REAL (the capped usc by construction cannot exceed Umax) instead of a
    tautology. Direction=+1, UOffset=0 in this smoke, so u is the signed raw
    perturbation u1-u0.
    """
    u = float(u_kcal)
    if acore <= 0.0:
        return u
    if u <= ubcore_kcal:
        return u
    y = (u - ubcore_kcal) / (umax_kcal - ubcore_kcal)
    z = 1.0 + 2.0 * (y / acore) + 2.0 * (y / acore) ** 2
    fsc = (z ** acore - 1.0) / (z ** acore + 1.0)
    return (umax_kcal - ubcore_kcal) * fsc + ubcore_kcal


def _make_context(system, positions, platform_name):
    """Build a context, preferring CUDA (5070Ti), falling back to Reference."""
    integ = mm.VerletIntegrator(0.001 * unit.picoseconds)
    try:
        plat = mm.Platform.getPlatformByName(platform_name)
        ctx = mm.Context(system, integ, plat)
    except Exception:
        plat = mm.Platform.getPlatformByName("Reference")
        ctx = mm.Context(system, integ, plat)
        platform_name = "Reference"
    ctx.setPositions(positions)
    return ctx, platform_name


def _prepare_leg_structure(final_pdb, out_pdb, leg, binder_chain):
    """Leg-aware structure prep: free => peptide-only; bound => receptor+peptide.
    Mirrors the prep build_inplace_res4_fused_system uses internally so the
    endpoint-equivalence reference is built on the SAME structure as the fused
    box (C1 — no divergent prep path)."""
    if leg == "bound":
        return ats.prepare_bound_complex_from_final(final_pdb, out_pdb, binder_chain)
    return ats.prepare_free_peptide_from_final(final_pdb, out_pdb, binder_chain)


def _pure_mtr_only_energy(seed, binder_chain, solvate, harmonize_common_charges,
                          leg="free"):
    """Build the pure single-topology MTR-only system (for the requested leg) and
    return its potential energy (the endpoint-equivalence reference for C6e).

    When ``harmonize_common_charges`` is set, the SAME common-charge override the
    fused build applies is applied here so the comparison is apples-to-apples
    (the fused u0's common core carries the harmonized charges).
    """
    li = ats.resolve_leg_inputs(seed)
    import tempfile
    td = tempfile.mkdtemp(prefix="ats_mtr_ref_%s_" % (leg,))
    mtr_struct = os.path.join(td, "cp4_%s.pdb" % (leg,))
    _prepare_leg_structure(li["final"]["cp4"], mtr_struct, leg, binder_chain)
    mtr = ats.build_leg_system(
        mtr_struct, leg=leg, binder_chain=binder_chain, solvate=solvate,
        add_hydrogens=False, hydrogens_xml=li["hydrogens_xml"])
    if harmonize_common_charges:
        wt_struct = os.path.join(td, "wt_%s.pdb" % (leg,))
        _prepare_leg_structure(li["final"]["wt"], wt_struct, leg, binder_chain)
        wt = ats.build_leg_system(
            wt_struct, leg=leg, binder_chain=binder_chain, solvate=False,
            add_hydrogens=False, hydrogens_xml=li["hydrogens_xml"])
        cmap = ats._build_common_index_map(wt, mtr, binder_chain)
        ats._harmonize_common_charges_to_wt(wt, mtr, cmap)
    ctx, _ = _make_context(mtr["system"], mtr["modeller"].positions, "Reference")
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole)
    return float(e)


def _endpoint_equivalence_unsolvated(seed, binder_chain, harmonize_common_charges,
                                     var_park_nm, lam, swap_mode="genuine",
                                     genuine_decouple_nm=1.2, leg="free"):
    """Compute the rigorous endpoint-equivalence pair (UNSOLVATED).

    Returns (fused_u0_unsolvated, mtr_only_energy). Both are built protein-only
    (no solvent variance), so u0 of the fused box and the pure MTR-only energy
    are directly comparable: u0 (the WT-coupled / reference endpoint) should
    equal the MTR-only energy plus HE1's residual nonbonded coupling, i.e.
    within a small tolerance. This is the only CHEAP correctness check (C6e/Q5),
    and it is evaluated at the TRUE physical mutation geometry — u0 is the
    REFERENCE (un-displaced) state, so the parked-clash / decoupled u1 geometry
    never enters this gate. A solvated cross-box comparison is confounded by
    independent water placement, so the check is always run unsolvated.
    """
    fused = ats.build_inplace_res4_fused_system(
        leg=leg, seed=seed, binder_chain=binder_chain, solvate=False,
        harmonize_common_charges=harmonize_common_charges,
        var_park_nm=var_park_nm, swap_mode=swap_mode,
        genuine_decouple_nm=genuine_decouple_nm)
    if fused.get("outcome") == "mc1_charge_discontinuity":
        # No fused box attached; endpoint-equivalence is moot in this branch.
        return float("nan"), float("nan")
    fb = fused["fused_build"]
    ctx, _ = _make_context(fb["system"], fb["modeller"].positions, "Reference")
    ctx.setParameter("Lambda1", 0.0)
    ctx.setParameter("Lambda2", 0.0)
    ctx.setParameter("Direction", 1.0)
    atm = next(fb["system"].getForce(i)
               for i in range(fb["system"].getNumForces())
               if isinstance(fb["system"].getForce(i), mm.ATMForce))
    u0 = atm.getPerturbationEnergy(ctx)[1].value_in_unit(
        unit.kilocalorie_per_mole)
    mtr_e = _pure_mtr_only_energy(
        seed, binder_chain, solvate=False,
        harmonize_common_charges=harmonize_common_charges, leg=leg)
    return float(u0), float(mtr_e)


def run_tier1_smoke(seed="s7", binder_chain="B", solvate=True,
                    harmonize_common_charges=False,
                    var_park_nm=None,
                    swap_mode="genuine", genuine_decouple_nm=1.2,
                    lam=0.5, platform_name="CUDA", leg="free"):
    """Build the in-place res-4 fused box and run the Tier-1 ASSERT decomposition.

    ``swap_mode="genuine"`` (default) produces a real, bond/angle-strain-free
    HE1<->methyl perturbation (the methyl stays coupled; HE1's nonbonded coupling
    is removed by translating it into bulk at u1) — a physically plausible
    ones-to-tens of kcal/mol |u1-u0| at the TRUE mutation geometry. The legacy
    ``swap_mode="var_park"`` (or any non-zero ``var_park_nm``) selects the
    DIAGNOSTIC parked-clash path (elevated, bond-strain-dominated |u1-u0|).

    Returns a structured result dict. Raises AssertionError on any C6 violation
    (so an executor sees a hard FAIL), EXCEPT for the pre-registered MC1
    charge-discontinuity outcome (ii), which is surfaced as a structured result
    with outcome="mc1_charge_discontinuity" (a real finding, not a crash).
    """
    build = ats.build_inplace_res4_fused_system(
        leg=leg, seed=seed, binder_chain=binder_chain, solvate=solvate,
        harmonize_common_charges=harmonize_common_charges,
        var_park_nm=var_park_nm, swap_mode=swap_mode,
        genuine_decouple_nm=genuine_decouple_nm)

    if build.get("outcome") == "mc1_charge_discontinuity":
        # Pre-registered outcome (ii): surface, do not crash. This IS the de-risk
        # signal that the production hybrid XML needs common-charge harmonization.
        return {
            "outcome": "mc1_charge_discontinuity",
            "regime": "ranking_only",
            "prediction_test": True,
            "mc1": build["mc1_param_continuity"],
            "note": build["note"],
        }

    fused = build["fused_build"]
    system = fused["system"]
    positions = fused["modeller"].positions

    ctx, platform_used = _make_context(system, positions, platform_name)
    atm = next(system.getForce(i) for i in range(system.getNumForces())
               if isinstance(system.getForce(i), mm.ATMForce))

    ctx.setParameter("Lambda1", lam)
    ctx.setParameter("Lambda2", lam)
    ctx.setParameter("Direction", 1.0)

    state = ctx.getState(getEnergy=True, getForces=True)
    e_atm = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    pert = atm.getPerturbationEnergy(ctx)
    u1 = pert[0].value_in_unit(unit.kilocalorie_per_mole)
    u0 = pert[1].value_in_unit(unit.kilocalorie_per_mole)
    energy_term = (pert[2].value_in_unit(unit.kilocalorie_per_mole)
                   if hasattr(pert[2], "value_in_unit") else float(pert[2]))

    forces = state.getForces().value_in_unit(
        unit.kilocalorie_per_mole / unit.angstrom)
    fmax = max(float(np.linalg.norm(f)) for f in forces)

    # |u1-u0| is the raw perturbation; usc is its soft-core-capped image, which by
    # construction cannot exceed Umax. Compute usc from the upstream soft-core
    # formula (UOffset=0, Direction=+1 here) so the C6d cap check is REAL, not a
    # tautology: usc = u for |u|<=Ubcore; above Ubcore it is compressed toward
    # Umax via the fsc softplus. Reproduces ommsystem.py:783-789 exactly.
    umax_kcal = ats.ATS_UMAX_KCAL
    raw_pert = u1 - u0
    usc_kcal = _softcore_usc(raw_pert, umax_kcal,
                             ats.ATS_UBCORE_KCAL, ats.ATS_ACORE)

    # --- C6 ASSERTS (hard) ---
    assert np.isfinite(e_atm), "C6a FAIL: E_ATM is not finite (NaN/inf)."
    assert np.isfinite(u0), "C6a FAIL: u0 is not finite."
    assert np.isfinite(u1), "C6a FAIL: u1 is not finite."
    assert np.isfinite(energy_term), "C6a FAIL: ATM energy term is not finite."
    assert abs(raw_pert) > U1_NE_U0_MIN_KCAL, (
        "C6b FAIL: u1 ~ u0 (|u1-u0|=%.3e <= %.3e) — the swap is a null op "
        "(pre-registered outcome iii)." % (abs(raw_pert), U1_NE_U0_MIN_KCAL))
    assert fmax < MAX_FORCE_KCAL_PER_MOL_A, (
        "C6c FAIL: max per-atom force %.1f kcal/mol/A exceeds %.1f — near-singular "
        "term, would NaN under integration." % (fmax, MAX_FORCE_KCAL_PER_MOL_A))
    # C6d (REAL assert, not `or True`): the soft-core-capped usc must respect Umax,
    # AND the raw |u1-u0| must be physically plausible for a HE1<->methyl edit
    # (ones-to-tens of kcal/mol). The latter is what distinguishes the GENUINE
    # swap (~10 kcal/mol) from the var_park parked clash (~1.5e4 kcal/mol).
    assert usc_kcal <= umax_kcal + 1e-6, (
        "C6d FAIL: usc %.3f kcal/mol exceeds Umax %.3f — soft-core cap not "
        "respected." % (usc_kcal, umax_kcal))
    assert abs(raw_pert) <= U1_MINUS_U0_PLAUSIBLE_MAX_KCAL, (
        "C6d FAIL: |u1-u0|=%.1f kcal/mol exceeds the plausible bound %.1f for a "
        "HE1<->methyl edit — this is a parked-clash / bond-strain artifact, NOT a "
        "genuine relative perturbation."
        % (abs(raw_pert), U1_MINUS_U0_PLAUSIBLE_MAX_KCAL))

    # --- C6e ENDPOINT-EQUIVALENCE: u0 reproduces the MTR-only energy ---
    # ALWAYS evaluated on the UNSOLVATED paired build (same protein coords, no
    # solvent variance). Comparing two INDEPENDENTLY-solvated boxes is NOT
    # apples-to-apples (each addSolvent draws a different water count/placement,
    # so the solvent term alone differs by >> tol and would confound C6e). The
    # unsolvated protein-only comparison is the rigorous endpoint check; the
    # solvated PME box is what proves finiteness/null-op (C6a-c) under PME.
    fused_u0_unsolv, mtr_only_e = _endpoint_equivalence_unsolvated(
        seed, binder_chain, harmonize_common_charges, var_park_nm, lam,
        swap_mode=swap_mode, genuine_decouple_nm=genuine_decouple_nm, leg=leg)
    endpoint_gap = fused_u0_unsolv - mtr_only_e
    endpoint_equiv = abs(endpoint_gap) <= ENDPOINT_EQUIV_TOL_KCAL

    result = {
        "outcome": ("tier1_pass" if endpoint_equiv else "tier1_endpoint_mismatch"),
        "regime": "ranking_only",
        "prediction_test": True,
        "leg": leg,
        "platform": platform_used,
        "solvated": solvate,
        "harmonize_common_charges": harmonize_common_charges,
        "swap_mode": swap_mode,
        "genuine_decouple_nm": (genuine_decouple_nm
                                if swap_mode == "genuine" else None),
        "genuine_decouple_dir": build.get("genuine_decouple_dir"),
        "var_park_nm": list(var_park_nm) if var_park_nm else None,
        "lambda": lam,
        "direction": 1.0,
        "n_particles": system.getNumParticles(),
        "energies_kcal": {
            "E_ATM": e_atm, "u0": u0, "u1": u1, "u1_minus_u0": raw_pert,
            "usc": usc_kcal, "atm_energy_term": energy_term,
        },
        "max_force_kcal_per_mol_A": fmax,
        "umax_kcal": umax_kcal,
        "usc_kcal": usc_kcal,
        "usc_le_umax": bool(usc_kcal <= umax_kcal + 1e-6),
        "u1_minus_u0_plausible": bool(
            abs(raw_pert) <= U1_MINUS_U0_PLAUSIBLE_MAX_KCAL),
        "endpoint_equivalence": {
            "basis": "unsolvated_paired (rigorous same-coord comparison)",
            "fused_u0_kcal": fused_u0_unsolv,
            "solvated_run_u0_kcal": u0,
            "mtr_only_energy_kcal": mtr_only_e,
            "gap_kcal": endpoint_gap,
            "tol_kcal": ENDPOINT_EQUIV_TOL_KCAL,
            "passed": bool(endpoint_equiv),
        },
        "swap": build["swap"],
        "mc1": build["mc1_param_continuity"],
        "mc2": build["mc2_methyl_bonded"],
        "mc3": build["mc3_disulfide"],
        "seed_assert": build["seed_assert"],
        "note": ("v0.8 PREDICTION test (R-18); ranking-only (R-11); NOT a "
                 "converged DDG_bind, NOT a Magotti/absolute comparison."),
    }
    # C6e is FLAG (structured), not a hard crash: endpoint mismatch is the
    # pre-registered outcome (ii) and points at a residual coupling/charge issue,
    # which is exactly the de-risk information we want surfaced, not hidden.
    return result


# ===========================================================================
# Tier-2: finite-UNDER-INTEGRATION (short dynamics stability at Lambda=0.5).
#
# Tier-1 proved finite-at-one-frame (0-step decomposition). Tier-2 proves the
# genuine fused box SURVIVES short MD at Lambda1=Lambda2=0.5 / Direction=+1
# without NaN / explosion. Per the R3 coordswap choice + the deep-dig S2/S3
# spec: 200-500 steps, 1 fs timestep, UNCONSTRAINED alchemical H (no SHAKE on
# the appearing/disappearing methyl/HE1 H), standard Langevin (LangevinMiddle),
# NO HMR. The constraint set is honored at the System level (createSystem), NOT
# the integrator: OpenMM SHAKE lives on the System, so an unconstrained run
# REBUILDS the genuine box with constraints=None (additive build-module param).
# ===========================================================================

def _constraints_enum(name):
    """Map a constraints label to the OpenMM createSystem constraint object.

    'none'   -> None              (no SHAKE; the R3 unconstrained-alch-H build)
    'hbonds' -> openmm.app.HBonds (SHAKE on every X-H, incl. the alch H)
    """
    key = (name or "none").lower()
    if key in ("none", "off", "unconstrained"):
        return None
    if key in ("hbonds", "h-bonds", "hbond"):
        return mm.app.HBonds
    raise ValueError("unknown constraints label %r (use none|hbonds)" % (name,))


def _build_genuine_fused_for_tier2(seed, binder_chain, solvate,
                                   harmonize_common_charges, genuine_decouple_nm,
                                   constraints_label, leg="free"):
    """Build the genuine fused box for a Tier-2 run with the requested
    constraint set, reusing build_inplace_res4_fused_system VERBATIM (C1).

    Returns (build_dict, constraints_enum). If the on-disk RBFE XML common core
    is discontinuous AND harmonization is off, returns the structured MC1
    outcome (no fused box) so the caller surfaces it honestly rather than
    integrating a charge-discontinuous box.
    """
    cons = _constraints_enum(constraints_label)
    build = ats.build_inplace_res4_fused_system(
        leg=leg, seed=seed, binder_chain=binder_chain, solvate=solvate,
        harmonize_common_charges=harmonize_common_charges,
        swap_mode="genuine", genuine_decouple_nm=genuine_decouple_nm,
        constraints=cons)
    return build, cons


def _run_short_dynamics(system, positions, dt_fs, n_steps, checkpoint_every,
                        platform_name, lam=0.5, minimize_iters=200):
    """Integrate the genuine fused box for n_steps at Lambda1=Lambda2=lam,
    Direction=+1, with a LangevinMiddleIntegrator (300 K, 1/ps) at dt_fs (no
    HMR — masses untouched), checkpointing PE/KE/T/Fmax every checkpoint_every.

    Returns a structured trajectory dict. Detects NaN/inf at every checkpoint
    and SHORT-CIRCUITS (records the failing step + the energy/force trace) — an
    honest failure is recorded, not forced to a pass. Raises nothing for a NaN;
    the caller decides PASS/FAIL from the returned ``survived`` flag.
    """
    integ = mm.LangevinMiddleIntegrator(
        TIER2_TEMPERATURE_K * unit.kelvin,
        TIER2_FRICTION_PER_PS / unit.picosecond,
        dt_fs * 0.001 * unit.picoseconds)
    try:
        plat = mm.Platform.getPlatformByName(platform_name)
        ctx = mm.Context(system, integ, plat)
    except Exception:
        plat = mm.Platform.getPlatformByName("Reference")
        ctx = mm.Context(system, integ, plat)
        platform_name = "Reference"
    ctx.setPositions(positions)
    ctx.setParameter("Lambda1", lam)
    ctx.setParameter("Lambda2", lam)
    ctx.setParameter("Direction", 1.0)

    # Light minimization to relieve the post-injection HE1 placement (the seed is
    # already gated in Tier-1; this only removes the worst local strain so the
    # FIRST integration step is not dominated by a transient). Kept small so it
    # does not mask a genuine dynamics pathology.
    if minimize_iters and minimize_iters > 0:
        ctx.setVelocitiesToTemperature(TIER2_TEMPERATURE_K * unit.kelvin)
        mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=minimize_iters)
    ctx.setVelocitiesToTemperature(TIER2_TEMPERATURE_K * unit.kelvin)

    def _checkpoint(step):
        st = ctx.getState(getEnergy=True, getForces=True)
        pe = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        ke = st.getKineticEnergy().value_in_unit(unit.kilocalorie_per_mole)
        forces = st.getForces().value_in_unit(
            unit.kilocalorie_per_mole / unit.angstrom)
        fmax = max(float(np.linalg.norm(f)) for f in forces)
        # Instantaneous temperature from KE: T = 2*KE/(ndf*kB). ndf via the
        # integrator-aware DOF count (constraints reduce ndf); approximate with
        # 3N - n_constraints - 3 (COM). We read ndf from the system robustly.
        ndf = _system_ndf(system)
        kb = 0.0019872041  # kcal/mol/K
        temp = (2.0 * ke / (ndf * kb)) if ndf > 0 else float("nan")
        finite = bool(np.isfinite(pe) and np.isfinite(ke) and np.isfinite(fmax))
        return {
            "step": step, "pe_kcal": pe, "ke_kcal": ke,
            "temp_K": temp, "fmax_kcal_per_mol_A": fmax, "finite": finite,
        }

    traj = [_checkpoint(0)]
    survived = traj[0]["finite"]
    failed_at = None if survived else 0
    nan_mode = None if survived else "nonfinite_at_start"

    step = 0
    while step < n_steps and survived:
        nblock = min(checkpoint_every, n_steps - step)
        try:
            integ.step(nblock)
        except Exception as exc:  # OpenMM raises on "Particle coordinate is NaN"
            survived = False
            failed_at = step + nblock
            nan_mode = "integrator_exception: %s" % (str(exc)[:200],)
            break
        step += nblock
        cp = _checkpoint(step)
        traj.append(cp)
        if not cp["finite"]:
            survived = False
            failed_at = step
            nan_mode = "nonfinite_checkpoint"
            break

    return {
        "dt_fs": dt_fs,
        "n_steps_requested": n_steps,
        "n_steps_completed": step,
        "checkpoint_every": checkpoint_every,
        "platform": platform_name,
        "trajectory": traj,
        "survived": bool(survived),
        "failed_at_step": failed_at,
        "nan_mode": nan_mode,
    }


def _system_ndf(system):
    """Degrees of freedom = 3*N_massive - n_constraints - 3 (COM removed if a
    CMMotionRemover is present). Massless particles (the dummy NE1 reference)
    contribute no DOF."""
    n_massive = 0
    for i in range(system.getNumParticles()):
        m = system.getParticleMass(i).value_in_unit(unit.dalton)
        if m > 0.0:
            n_massive += 1
    n_constr = system.getNumConstraints()
    has_cmm = any(isinstance(system.getForce(i), mm.CMMotionRemover)
                  for i in range(system.getNumForces()))
    ndf = 3 * n_massive - n_constr - (3 if has_cmm else 0)
    return ndf


def _assert_tier2_stability(run, t_band=TIER2_T_BAND_K,
                            fmax_max=TIER2_MAX_FORCE_KCAL_PER_MOL_A,
                            equil_skip=TIER2_T_EQUIL_SKIP):
    """Hard asserts (not prints) on a primary Tier-2 run. Raises AssertionError
    on any violation so an executor sees a FAIL. The integrator-matrix non-
    primary cells are EVIDENCE (not asserted) — only the primary calls this.

    Every checkpoint must be finite, bounded in Fmax, and BELOW the runaway
    upper temperature edge (the load-bearing explosion guard). The LOWER
    temperature edge is enforced only AFTER ``equil_skip`` checkpoints (the
    thermalization transient — see the TIER2_T_BAND_K rationale)."""
    assert run["survived"], (
        "TIER2 FAIL: did NOT complete %d steps — %s at step %s. The genuine "
        "fused box exploded/NaN'd under integration (finite-at-one-frame did "
        "not survive dynamics)." % (run["n_steps_requested"], run["nan_mode"],
                                    run["failed_at_step"]))
    assert run["n_steps_completed"] == run["n_steps_requested"], (
        "TIER2 FAIL: completed %d/%d steps."
        % (run["n_steps_completed"], run["n_steps_requested"]))
    for cidx, cp in enumerate(run["trajectory"]):
        assert cp["finite"], (
            "TIER2 FAIL: non-finite PE/KE/Fmax at step %d." % cp["step"])
        assert cp["fmax_kcal_per_mol_A"] < fmax_max, (
            "TIER2 FAIL: max per-atom force %.1f kcal/mol/A at step %d exceeds "
            "%.1f." % (cp["fmax_kcal_per_mol_A"], cp["step"], fmax_max))
        # Upper (runaway) bound: ALWAYS enforced — a runaway box explodes upward.
        assert cp["temp_K"] <= t_band[1], (
            "TIER2 FAIL: temperature %.1f K at step %d exceeds the runaway upper "
            "bound %.0f K." % (cp["temp_K"], cp["step"], t_band[1]))
        # Lower bound: only after the equilibration window (transient undershoot).
        if cidx >= equil_skip:
            assert cp["temp_K"] >= t_band[0], (
                "TIER2 FAIL: temperature %.1f K at step %d below the sane lower "
                "bound %.0f K (post-equilibration frozen/sunk box)."
                % (cp["temp_K"], cp["step"], t_band[0]))


def run_tier2_stability(seed="s7", binder_chain="B", solvate=True,
                        harmonize_common_charges=False, genuine_decouple_nm=1.2,
                        steps=TIER2_DEFAULT_STEPS, dt_fs=1.0,
                        constraints_label="none",
                        checkpoint_every=TIER2_CHECKPOINT_EVERY,
                        platform_name="CUDA",
                        run_integrator_matrix=True, lam=0.5, leg="free"):
    """Tier-2 finite-under-integration runner (the primary R3 cell + the S3
    integrator-evidence matrix).

    Primary cell = (dt=1 fs, constraints=none) — the R3 choice; it MUST pass
    (hard asserts). The integrator matrix runs a short window for each
    {1fs/none, 2fs/none, 2fs/hbonds} and RECORDS which survive vs NaN, so the R3
    choice is evidence-backed (it is fine/expected that a 2 fs step on an
    appearing/disappearing H is less stable — that IS the evidence).

    Returns a structured result dict. Raises AssertionError on a PRIMARY failure
    (so an executor sees FAIL) EXCEPT the pre-registered MC1 outcome (surfaced
    structured). An honest Tier-2 failure (box does not survive 200-500 steps)
    is recorded with the failure mode + step + the energy/force trace — NOT
    forced to a pass (R-18). A genuine failure would indicate the single-NE1-ref
    decouple geometry has a dynamics pathology the 0-step Tier-1 cannot see,
    pointing at needing the full two-copy overlay before Tier-2 can pass.
    """
    build, cons = _build_genuine_fused_for_tier2(
        seed, binder_chain, solvate, harmonize_common_charges,
        genuine_decouple_nm, constraints_label, leg=leg)
    if build.get("outcome") == "mc1_charge_discontinuity":
        return {
            "tier": "tier2",
            "outcome": "mc1_charge_discontinuity",
            "regime": "ranking_only",
            "prediction_test": True,
            "mc1": build["mc1_param_continuity"],
            "note": build["note"],
        }

    fused = build["fused_build"]
    system = fused["system"]
    positions = fused["modeller"].positions

    # --- PRIMARY run (R3: dt=1 fs, constraints=none, no HMR) ---
    primary = _run_short_dynamics(
        system, positions, dt_fs=dt_fs, n_steps=steps,
        checkpoint_every=checkpoint_every, platform_name=platform_name, lam=lam)

    # --- S3 integrator-evidence matrix (each cell on a FRESH box; a Context
    #     mutates positions/velocities, so reuse would confound the cells) ---
    matrix = []
    if run_integrator_matrix:
        for label, m_dt, m_cons in TIER2_INTEGRATOR_MATRIX:
            mbuild, _ = _build_genuine_fused_for_tier2(
                seed, binder_chain, solvate, harmonize_common_charges,
                genuine_decouple_nm, m_cons, leg=leg)
            if mbuild.get("outcome") == "mc1_charge_discontinuity":
                matrix.append({"cell": label, "dt_fs": m_dt,
                               "constraints": m_cons,
                               "outcome": "mc1_charge_discontinuity"})
                continue
            mfused = mbuild["fused_build"]
            mrun = _run_short_dynamics(
                mfused["system"], mfused["modeller"].positions,
                dt_fs=m_dt, n_steps=TIER2_MATRIX_STEPS,
                checkpoint_every=checkpoint_every, platform_name=platform_name,
                lam=lam)
            matrix.append({
                "cell": label, "dt_fs": m_dt, "constraints": m_cons,
                "survived": mrun["survived"],
                "n_steps_completed": mrun["n_steps_completed"],
                "n_steps_requested": mrun["n_steps_requested"],
                "failed_at_step": mrun["failed_at_step"],
                "nan_mode": mrun["nan_mode"],
                "final_pe_kcal": mrun["trajectory"][-1]["pe_kcal"],
                "final_temp_K": mrun["trajectory"][-1]["temp_K"],
                "max_fmax_kcal_per_mol_A": max(
                    cp["fmax_kcal_per_mol_A"] for cp in mrun["trajectory"]),
            })

    # --- PRIMARY hard asserts (after the matrix so the evidence is captured even
    #     if the primary later fails the assert) ---
    primary_pass = False
    assert_error = None
    try:
        _assert_tier2_stability(primary)
        primary_pass = True
    except AssertionError as exc:
        assert_error = str(exc)

    pe_series = [cp["pe_kcal"] for cp in primary["trajectory"]]
    t_series = [cp["temp_K"] for cp in primary["trajectory"]]
    result = {
        "tier": "tier2",
        "outcome": ("tier2_pass" if primary_pass else "tier2_fail"),
        "regime": "ranking_only",
        "prediction_test": True,
        "leg": leg,
        "platform": primary["platform"],
        "solvated": solvate,
        "harmonize_common_charges": harmonize_common_charges,
        "swap_mode": "genuine",
        "genuine_decouple_nm": genuine_decouple_nm,
        "genuine_decouple_dir": build.get("genuine_decouple_dir"),
        "lambda": lam,
        "direction": 1.0,
        "n_particles": system.getNumParticles(),
        "primary": {
            "dt_fs": dt_fs,
            "constraints": constraints_label,
            "hmr": False,
            "integrator": "LangevinMiddle 300K 1/ps",
            "n_steps_requested": steps,
            "n_steps_completed": primary["n_steps_completed"],
            "survived": primary["survived"],
            "failed_at_step": primary["failed_at_step"],
            "nan_mode": primary["nan_mode"],
            "pe_kcal_series": pe_series,
            "temp_K_series": t_series,
            "pe_kcal_first": pe_series[0],
            "pe_kcal_last": pe_series[-1],
            "max_fmax_kcal_per_mol_A": max(
                cp["fmax_kcal_per_mol_A"] for cp in primary["trajectory"]),
            "trajectory": primary["trajectory"],
        },
        "integrator_matrix": matrix,
        "mc1": build["mc1_param_continuity"],
        "common_charges_harmonized": build.get("common_charges_harmonized"),
        "thresholds": {
            "max_force_kcal_per_mol_A": TIER2_MAX_FORCE_KCAL_PER_MOL_A,
            "temp_band_K": list(TIER2_T_BAND_K),
        },
        "note": ("v0.8 PREDICTION test (R-18); ranking-only (R-11). Tier-2 = "
                 "finite-UNDER-INTEGRATION (short dynamics), NOT a converged "
                 "DDG_bind."),
    }
    if assert_error is not None:
        result["assertion_error"] = assert_error
    return result


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Track B in-place residue-4 ATS Tier-1 ASSERT smoke "
                    "(v0.8 PREDICTION test; ranking-only).")
    p.add_argument("--seed", default="s7")
    p.add_argument("--binder-chain", default="B")
    p.add_argument("--leg", choices=["free", "bound"], default="free",
                   help="RBFE cycle leg. 'free' (default) = solvated cyclic "
                        "peptide alone. 'bound' = receptor (C3c) + peptide in the "
                        "equilibrated BOUND pose, re-solvated, NO displacement "
                        "(RBFE — the binder stays bound). The DDG_bind = "
                        "DG_mut(bound) - DG_mut(free) requires both legs.")
    p.add_argument("--no-solvate", action="store_true",
                   help="Build unsolvated (NoCutoff) instead of the PME C6 "
                        "target (cheap CPU path; the C6 target is solvated).")
    p.add_argument("--harmonize-common-charges", action="store_true",
                   help="DIAGNOSTIC: force MTR common charges to WT so the "
                        "common core is continuous (isolates the mechanical "
                        "validation from the MC1 charge gap; NOT production "
                        "valid for a quantitative DDG).")
    p.add_argument("--swap-mode", choices=["genuine", "var_park", "null"],
                   default="genuine",
                   help="Var-region swap construction. 'genuine' (default) = "
                        "the real HE1<->methyl coupling toggle (physically "
                        "plausible |u1-u0|, true mutation geometry). 'var_park' "
                        "= the DIAGNOSTIC parked-clash (requires --var-park-nm). "
                        "'null' = the degenerate single-shared-core null op.")
    p.add_argument("--genuine-decouple-nm", type=float, default=1.2,
                   help="Genuine-mode HE1 decouple distance into bulk (nm); the "
                        "dummy NE1-reference offset. Default 1.2 (HE1 fully out "
                        "of the core's interaction range).")
    p.add_argument("--var-park-nm", type=float, nargs=3, default=None,
                   metavar=("X", "Y", "Z"),
                   help="DIAGNOSTIC park vector (the ~4-atom methyl) for "
                        "--swap-mode var_park. Elevated/parked-clash |u1-u0|; "
                        "not a genuine perturbation.")
    p.add_argument("--lambda", dest="lam", type=float, default=0.5)
    p.add_argument("--platform", default="CUDA")
    # --- Tier-2 (finite-under-integration short dynamics) ---
    p.add_argument("--tier2", action="store_true",
                   help="Run the Tier-2 short-dynamics stability runner "
                        "(finite-UNDER-INTEGRATION) instead of the Tier-1 0-step "
                        "decomposition. Builds the genuine fused box, integrates "
                        "N steps at Lambda=0.5/Direction=+1 with a 1 fs "
                        "UNCONSTRAINED-alch-H LangevinMiddle integrator (no HMR), "
                        "and runs the S3 integrator-evidence matrix.")
    p.add_argument("--tier2-steps", type=int, default=TIER2_DEFAULT_STEPS,
                   help="Tier-2 primary-run step count (200-500; default %d)."
                        % TIER2_DEFAULT_STEPS)
    p.add_argument("--tier2-dt-fs", type=float, default=1.0,
                   help="Tier-2 primary-run timestep in fs (R3 default 1.0; the "
                        "unconstrained-alch-H build needs the small step).")
    p.add_argument("--tier2-constraints", choices=["none", "hbonds"],
                   default="none",
                   help="Tier-2 primary-run constraint set. 'none' (R3 default) "
                        "= NO SHAKE on the appearing/disappearing alch H. "
                        "'hbonds' = SHAKE on every X-H (diagnostic).")
    p.add_argument("--tier2-checkpoint-every", type=int,
                   default=TIER2_CHECKPOINT_EVERY,
                   help="Tier-2 checkpoint interval in steps (default %d)."
                        % TIER2_CHECKPOINT_EVERY)
    p.add_argument("--tier2-no-matrix", action="store_true",
                   help="Skip the S3 integrator-evidence matrix (run only the "
                        "primary cell).")
    p.add_argument("--json-out", default=None,
                   help="Write the structured result JSON to this path "
                        "(default: outputs/_trackb/inplace_res4_ats_smoke/"
                        "tier{1,2}_<...>_<seed>.json).")
    args = p.parse_args(argv)

    # --- Tier-2 dispatch (separate path; the Tier-1 flow below is untouched) ---
    if args.tier2:
        json_out = args.json_out
        if json_out is None:
            out_dir = os.path.join(_PROJ, "outputs", "_trackb",
                                   "inplace_res4_ats_smoke")
            os.makedirs(out_dir, exist_ok=True)
            json_out = os.path.join(
                out_dir, "tier2_genuine_%s_%s.json" % (args.leg, args.seed))
        result = run_tier2_stability(
            seed=args.seed, binder_chain=args.binder_chain,
            solvate=not args.no_solvate,
            harmonize_common_charges=args.harmonize_common_charges,
            genuine_decouple_nm=args.genuine_decouple_nm,
            steps=args.tier2_steps, dt_fs=args.tier2_dt_fs,
            constraints_label=args.tier2_constraints,
            checkpoint_every=args.tier2_checkpoint_every,
            platform_name=args.platform,
            run_integrator_matrix=not args.tier2_no_matrix, lam=args.lam,
            leg=args.leg)
        print(json.dumps(result, indent=2, default=str))
        with open(json_out, "w") as fh:
            json.dump(result, fh, indent=2, default=str)
        # Exit codes: 0 = tier2_pass; 6 = mc1 outcome; 8 = tier2_fail.
        if result.get("outcome") == "tier2_pass":
            return 0
        if result.get("outcome") == "mc1_charge_discontinuity":
            return 6
        return 8

    var_park = tuple(args.var_park_nm) if args.var_park_nm is not None else None
    if var_park is not None and not any(abs(c) > 0 for c in var_park):
        var_park = None

    json_out = args.json_out
    if json_out is None:
        out_dir = os.path.join(_PROJ, "outputs", "_trackb",
                               "inplace_res4_ats_smoke")
        os.makedirs(out_dir, exist_ok=True)
        json_out = os.path.join(
            out_dir, "tier1_%s_%s_%s.json" % (args.swap_mode, args.leg,
                                              args.seed))

    try:
        result = run_tier1_smoke(
            seed=args.seed, binder_chain=args.binder_chain,
            solvate=not args.no_solvate,
            harmonize_common_charges=args.harmonize_common_charges,
            var_park_nm=var_park, swap_mode=args.swap_mode,
            genuine_decouple_nm=args.genuine_decouple_nm,
            lam=args.lam, platform_name=args.platform, leg=args.leg)
    except AssertionError as exc:
        result = {"outcome": "tier1_fail", "assertion_error": str(exc),
                  "regime": "ranking_only", "prediction_test": True,
                  "leg": args.leg, "swap_mode": args.swap_mode}
        print(json.dumps(result, indent=2))
        with open(json_out, "w") as fh:
            json.dump(result, fh, indent=2)
        return 5

    print(json.dumps(result, indent=2, default=str))
    with open(json_out, "w") as fh:
        json.dump(result, fh, indent=2, default=str)

    # Exit codes: 0 = tier1_pass; 6 = mc1 outcome; 7 = endpoint mismatch.
    if result.get("outcome") == "tier1_pass":
        return 0
    if result.get("outcome") == "mc1_charge_discontinuity":
        return 6
    return 7


if __name__ == "__main__":
    sys.exit(main())
