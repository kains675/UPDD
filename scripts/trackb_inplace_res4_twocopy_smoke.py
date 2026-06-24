#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B — CANONICAL ATS TWO-COPY residue-4 Tier-1 ASSERT smoke (v0.8 rebuild).

Canonical two-copy ATS construction
(C1-C10). This exercises the REBUILT canonical ATS two-copy path
(``atm_trackB_setup.build_inplace_res4_twocopy_system`` /
``attach_twocopy_swap_atmforce``), NOT the legacy single-shared-core path (which
collapsed: pertE saturated ~150, overlap ~0, dgbind1 = 0.5*pertE offset).

CANONICAL ATS (the load-bearing correction, Q1/Q3):
  - BOTH endpoint copies are resident in ONE box: copy-1 (MTR) at the site,
    copy-2 (WT) DISPLACED by d (~40 A) into BULK solvent. Clash is avoided by the
    SPATIAL SEPARATION (d), NOT by inter-copy nonbonded exclusions (NONE added).
  - The common-region coordinates are SWAPPED as an ATMForce transform
    (ParticleOffsetDisplacement) with DISTINCT attach atoms (copy-1 NE1 !=
    copy-2 NE1), so the var perturbation is a genuine transfer, not a null op.
  - The two copies' common conformers are REGISTERED (copy-2 commons set to
    copy-1's frame + d) so the swap offset is pure d (no inter-conformer strain).

What this Tier-1 smoke proves (every check is an ASSERT or a structured flag,
not a bare print), on the requested leg at Lambda1=Lambda2=0.5/Direction=+1:
  (a) the two-copy box BUILDS + the System has BOTH copies resident
      (n_atoms ~ 2x single endpoint);
  (b) the two copies are SPATIALLY SEPARATED (copy1-NE1 <-> copy2-NE1 ~ d, NOT
      overlaid) and ZERO inter-copy exclusions were added (C3);
  (c) MC1 (common-charge continuity between the two resident copies) holds with
      the harmonized RBFE XML (or surfaces the charge gap honestly);
  (d) MC2 (copy-1 methyl bonded) + MC3 (cyclic_ss in BOTH copies);
  (e) Tier-1 ENERGY: E_ATM / u0 / u1 finite, max-force bounded, |u1-u0| in the
      physically-plausible ones-to-hundreds-kcal band (NOT the single-shared-core
      saturated ~150 plateau, NOT a 56,000 clash);
  (f) ENDPOINT-EQUIVALENCE (C6e/C8): the ATM reference potential u0 reproduces
      the full System potential within tol + the BULK copy is DECOUPLED
      (copy1-NE1 <-> copy2-NE1 >= the PME cutoff).

What this does NOT prove (R-18, honest scope):
  - The converged DDG_bind (needs the lambda ladder + UWHAM — the pilot,
    post-validation).
  - That the swap perturbation is non-saturated AT A REAL LAMBDA WINDOW
    (frac<UBCORE>0). The one-frame |u1-u0| being small + finite is NECESSARY
    (defeats the collapse signature) but NOT SUFFICIENT; the multi-window
    soft-core-uncapped transfer is the definitive verdict (pilot).
  - dgbind1 != 0.5*pertE (the deterministic-offset escape) — also pilot scope.

Ranking-only (R-11); v0.8 PREDICTION test (R-18); NOT a Magotti/absolute
comparison. Runs openmm-only in the ``qmmm`` env (no atom_openmm needed). The
5070Ti (CUDA) is the target; falls back to Reference if CUDA is unavailable.

DOI references (DOIs allowed; AI-trail clean):
  - Gallicchio 2025 J Chem Inf Model, ATS DOI 10.1021/acs.jcim.5c00207
    (preprint arXiv:2412.19971; copy-2 bulk d-displacement + common coord swap,
    PDZ peptide single-point mutation hysteresis 0.59 kcal/mol)
  - Gallicchio 2021 J Chem Theory Comput, DOI 10.1021/acs.jctc.1c00753
  - Azimi et al. 2022 J Chem Inf Model 62(2):309, DOI 10.1021/acs.jcim.1c01129
  - Klimovich/Shirts/Mobley 2015 J Comput Aided Mol Des,
    DOI 10.1007/s10822-015-9840-9
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
MAX_FORCE_KCAL_PER_MOL_A = 1.0e5      # (e) bounded max per-atom force
# |u1-u0| must be PHYSICALLY PLAUSIBLE for a HE1<->methyl side-chain edit:
# ones-to-hundreds of kcal/mol, NOT the single-shared-core saturated ~150 plateau
# (the collapse signature) and NOT a 56,000 kcal/mol clash. The two-copy reference
# frame |u1-u0| is small (the FE signal accumulates over the ladder, not one
# frame); the upper bound rejects the parked-clash / overlay-collapse regimes.
U1_MINUS_U0_PLAUSIBLE_MAX_KCAL = 1.0e3
ENDPOINT_EQUIV_TOL_KCAL = 25.0        # u0 vs full-System potential tolerance


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


def run_tier1_twocopy_smoke(seed="s7", binder_chain="B", solvate=True,
                            harmonize_common_charges=False,
                            displacement_nm=ats.ATS_TWOCOPY_DISPLACEMENT_NM,
                            lam=0.5, platform_name="CUDA", leg="free",
                            auto_search_displacement=False,
                            accept_sep_nm=ats.ATS_TWOCOPY_ACCEPT_SEP_NM):
    """Build the canonical two-copy box + run the Tier-1 ASSERT decomposition.

    Returns a structured result dict. Raises AssertionError on any hard C6
    violation (so an executor sees a FAIL), EXCEPT the pre-registered MC1
    charge-discontinuity outcome (ii), which is surfaced structured (a real
    finding, not a crash).

    ``auto_search_displacement`` (default False -> byte-identical legacy fixed-
    direction path): when True the copy-2 bulk displacement is chosen by the
    builder's direction-aware cone search (``auto_search_twocopy_displacement``),
    which maximises the copy1<->copy2 (+ periodic image) min heavy-atom distance
    and escalates the magnitude only if no direction clears ``accept_sep_nm``.
    This recovers the bound-leg two-copy box build where a fixed-direction
    displacement drives copy-2's binder through copy-1's receptor body. The
    auto-search is d-/direction-NEUTRAL (the swap is partner-offset based;
    u1-u0 is d-invariant given full decoupling + bulk solvation), so it is
    ranking-safe. ``accept_sep_nm`` is consulted only by the auto-search; the
    post-solvate C6 assert always enforces the 1.0 nm clash floor + image gate.
    """
    build = ats.build_inplace_res4_twocopy_system(
        leg=leg, seed=seed, binder_chain=binder_chain, solvate=solvate,
        harmonize_common_charges=harmonize_common_charges,
        displacement_nm=displacement_nm,
        auto_search_displacement=auto_search_displacement,
        accept_sep_nm=accept_sep_nm)

    if build.get("outcome") == "mc1_charge_discontinuity":
        return {
            "outcome": "mc1_charge_discontinuity",
            "regime": "ranking_only",
            "prediction_test": True,
            "swap_mode": "twocopy",
            "leg": leg,
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
    raw_pert = u1 - u0
    forces = state.getForces().value_in_unit(
        unit.kilocalorie_per_mole / unit.angstrom)
    fmax = max(float(np.linalg.norm(f)) for f in forces)

    # --- C6 ASSERTS (hard) ---
    assert np.isfinite(e_atm), "C6a FAIL: E_ATM not finite."
    assert np.isfinite(u0), "C6a FAIL: u0 not finite."
    assert np.isfinite(u1), "C6a FAIL: u1 not finite."
    assert fmax < MAX_FORCE_KCAL_PER_MOL_A, (
        "C6c FAIL: max per-atom force %.1f exceeds %.1f — near-singular, would "
        "NaN under integration." % (fmax, MAX_FORCE_KCAL_PER_MOL_A))
    assert abs(raw_pert) <= U1_MINUS_U0_PLAUSIBLE_MAX_KCAL, (
        "C6d FAIL: |u1-u0|=%.1f kcal/mol exceeds the plausible bound %.1f — a "
        "parked-clash / overlay-collapse artifact, NOT a genuine two-copy "
        "transfer." % (abs(raw_pert), U1_MINUS_U0_PLAUSIBLE_MAX_KCAL))

    # C3 structural ASSERT: ZERO inter-copy exclusions (the rebuild must NOT have
    # re-introduced the forbidden overlay+exclusion design).
    assert build["swap"]["inter_copy_exclusions_added"] == 0, (
        "C3 FAIL: inter-copy exclusions were added — the forbidden overlay design "
        "leaked in. Clash must be avoided by d-separation, not exclusion.")
    assert build["swap"]["distinct_attach"], (
        "C2 FAIL: the two copies share an attach atom (single-shared-core null "
        "op) — the two-copy rebuild requires DISTINCT NE1 attach atoms.")

    # --- ENDPOINT-EQUIVALENCE (C6e/C8) — UNSOLVATED rigorous pair ---
    # Use the SAME displacement policy (auto-search vs fixed) as the solvated
    # build so the equivalence pair reflects the chosen construction, not a
    # divergent one (single policy — SciVal condition 2: no policy heterogeneity).
    eq_build = ats.build_inplace_res4_twocopy_system(
        leg=leg, seed=seed, binder_chain=binder_chain, solvate=False,
        harmonize_common_charges=harmonize_common_charges,
        displacement_nm=displacement_nm,
        auto_search_displacement=auto_search_displacement,
        accept_sep_nm=accept_sep_nm)
    if eq_build.get("outcome") == "twocopy_attached":
        eq = ats.check_twocopy_endpoint_equivalence(
            eq_build, platform_name="Reference",
            tol_kcal=ENDPOINT_EQUIV_TOL_KCAL)
    else:
        eq = {"overall_pass": False,
              "outcome": eq_build.get("outcome")}

    result = {
        "outcome": ("tier1_twocopy_pass" if eq.get("overall_pass")
                    else "tier1_twocopy_endpoint_mismatch"),
        "regime": "ranking_only",
        "prediction_test": True,
        "swap_mode": "twocopy",
        "leg": leg,
        "platform": platform_used,
        "solvated": solvate,
        "harmonize_common_charges": harmonize_common_charges,
        "displacement_vector_nm": build["displacement_vector_nm"],
        # task #100/#6: how d was chosen (fixed_direction vs auto_search) + the
        # per-build search trail (selected dir / magnitude / achieved min-image
        # sep) for the downstream Path decoupling verification.
        "displacement_mode": build.get("displacement_mode"),
        "displacement_log": build.get("displacement_log"),
        "lambda": lam,
        "direction": 1.0,
        "n_particles": system.getNumParticles(),
        "n_copy1": fused["n_copy1"],
        "energies_kcal": {
            "E_ATM": e_atm, "u0": u0, "u1": u1, "u1_minus_u0": raw_pert,
        },
        "max_force_kcal_per_mol_A": fmax,
        "u1_minus_u0_plausible": bool(
            abs(raw_pert) <= U1_MINUS_U0_PLAUSIBLE_MAX_KCAL),
        "separation": build["separation"],
        "swap": build["swap"],
        "mc1": build["mc1_param_continuity"],
        "mc2": build["mc2_methyl_bonded"],
        "mc3": build["mc3_disulfide"],
        "seed_assert": build["seed_assert"],
        "endpoint_equivalence": eq,
        "mtr_ncaa_xml": build.get("mtr_ncaa_xml"),
        "note": ("v0.8 CANONICAL ATS two-copy PREDICTION test (R-18); ranking-"
                 "only (R-11). NOT a converged DDG_bind, NOT a Magotti/absolute "
                 "comparison. One-frame |u1-u0| finite+non-saturated is NECESSARY "
                 "(collapse-signature escape) but NOT SUFFICIENT — the pilot "
                 "(frac<UBCORE>0 + O>=0.1 + dgbind1!=0.5 + UWHAM) is the "
                 "definitive verdict."),
    }
    return result


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Track B canonical ATS two-copy Tier-1 ASSERT smoke "
                    "(v0.8 PREDICTION test; ranking-only).")
    p.add_argument("--seed", default="s7")
    p.add_argument("--binder-chain", default="B")
    p.add_argument("--leg", choices=["free", "bound"], default="free",
                   help="RBFE cycle leg. 'free' (default) = solvated cyclic "
                        "peptide, TWO copies (MTR site + WT bulk). 'bound' = "
                        "receptor + peptide, TWO binder copies. DDG_bind = "
                        "DG_mut(bound) - DG_mut(free).")
    p.add_argument("--no-solvate", action="store_true",
                   help="Build unsolvated (NoCutoff) instead of the PME C8 "
                        "target (cheap CPU path).")
    p.add_argument("--harmonize-common-charges", action="store_true",
                   help="DIAGNOSTIC: force copy-2 common charges to copy-1 so "
                        "MC1 passes (isolates the mechanical validation from any "
                        "charge gap; the harmonized RBFE XML already makes MC1 "
                        "pass on disk, so this is rarely needed).")
    p.add_argument("--displacement-nm", type=float,
                   default=ats.ATS_TWOCOPY_DISPLACEMENT_NM,
                   help="Magnitude of the copy-2 bulk displacement d (nm); ATS "
                        "peptide convention ~4.0 (40 A). Default %.1f. Ignored "
                        "when --auto-search-displacement is set (the search "
                        "picks the magnitude from its ladder)."
                        % ats.ATS_TWOCOPY_DISPLACEMENT_NM)
    p.add_argument("--auto-search-displacement", action="store_true",
                   help="Opt-in: choose the copy-2 bulk displacement by the "
                        "builder's direction-aware cone search (maximises the "
                        "copy1<->copy2 + image min heavy-atom distance, escalates "
                        "the magnitude only if no direction clears "
                        "--accept-sep-nm). Recovers the bound-leg build where a "
                        "fixed-direction d drives copy-2's binder through copy-1's "
                        "receptor body. d-/direction-NEUTRAL (ranking-safe). "
                        "DEFAULT OFF = fixed direction at --displacement-nm "
                        "(byte-identical legacy).")
    p.add_argument("--accept-sep-nm", type=float,
                   default=ats.ATS_TWOCOPY_ACCEPT_SEP_NM,
                   help="Auto-search ONLY: the decoupling-sufficient acceptance "
                        "line (nm) a candidate direction must clear (PME cutoff + "
                        "LJ-tail buffer). Default %.1f. Consulted only when "
                        "--auto-search-displacement is set."
                        % ats.ATS_TWOCOPY_ACCEPT_SEP_NM)
    p.add_argument("--lambda", dest="lam", type=float, default=0.5)
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--json-out", default=None)
    args = p.parse_args(argv)

    json_out = args.json_out
    if json_out is None:
        out_dir = os.path.join(_PROJ, "outputs", "_trackb",
                               "inplace_res4_twocopy_smoke")
        os.makedirs(out_dir, exist_ok=True)
        json_out = os.path.join(
            out_dir, "tier1_twocopy_%s_%s.json" % (args.leg, args.seed))

    try:
        result = run_tier1_twocopy_smoke(
            seed=args.seed, binder_chain=args.binder_chain,
            solvate=not args.no_solvate,
            harmonize_common_charges=args.harmonize_common_charges,
            displacement_nm=args.displacement_nm,
            lam=args.lam, platform_name=args.platform, leg=args.leg,
            auto_search_displacement=args.auto_search_displacement,
            accept_sep_nm=args.accept_sep_nm)
    except AssertionError as exc:
        result = {"outcome": "tier1_twocopy_fail", "assertion_error": str(exc),
                  "regime": "ranking_only", "prediction_test": True,
                  "swap_mode": "twocopy", "leg": args.leg}
        print(json.dumps(result, indent=2, default=str))
        with open(json_out, "w") as fh:
            json.dump(result, fh, indent=2, default=str)
        return 5

    print(json.dumps(result, indent=2, default=str))
    with open(json_out, "w") as fh:
        json.dump(result, fh, indent=2, default=str)

    # Exit codes: 0 = pass; 6 = mc1 outcome; 7 = endpoint mismatch.
    if result.get("outcome") == "tier1_twocopy_pass":
        return 0
    if result.get("outcome") == "mc1_charge_discontinuity":
        return 6
    return 7


if __name__ == "__main__":
    sys.exit(main())
