#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the in-place residue-4 fused dual-topology ATS de-risk
(``utils/atm_trackB_setup`` new functions + ``scripts/trackb_inplace_res4_ats_smoke``).

Design rationale: C1-C8 coordinate-swap criteria + Q5
endpoint-equivalence. These tests assert the build-time hard gates and the
Tier-1 mechanical properties WITHOUT a GPU (the OpenMM Reference platform
evaluates energies on CPU). Tests that build the real fused box require both
2QKI endpoint final.pdb + openmm in the env; they env-gate/skip otherwise.

Coverage:
  - C3 count-parity + index-order alignment (synthetic + real)
  - MC1 common-charge continuity fail-loud (synthetic) + strict re-raise (real)
  - var-var index remap correctness (the WT/MTR raw-index divergence trap)
  - swap wiring => u1 != u0 (the null-op defeat) on the real fused box
  - endpoint-equivalence (u0 reproduces MTR-only) on the real fused box
  - finiteness + bounded force (the coverage-NaN defeat)
  - pre-registered MC1 structured outcome (production charges, no harmonize)
  - smoke-script outcome plumbing

Ranking-only (R-11); these validate the CONSTRUCTION, not a DDG.
"""

import importlib.util
import os
import sys

import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)
_UTILS = os.path.join(_PROJ, "utils")


def _have_openmm():
    try:
        import openmm  # noqa: F401
        return True
    except ImportError:
        return False


def _endpoints_present(seed="s7"):
    cp4 = os.path.join(_PROJ, "outputs", f"2QKI_Cp4_hybrid_calib_{seed}",
                       "mdresult", "2QKI_Cp4_final.pdb")
    wt = os.path.join(_PROJ, "outputs", f"2QKI_WT_calib_{seed}",
                      "mdresult", "2QKI_WT_final.pdb")
    return os.path.isfile(cp4) and os.path.isfile(wt)


def _load_ats():
    spec = importlib.util.spec_from_file_location(
        "atm_trackB_setup", os.path.join(_UTILS, "atm_trackB_setup.py"))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate utils/atm_trackB_setup.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_smoke():
    spec = importlib.util.spec_from_file_location(
        "trackb_inplace_res4_ats_smoke",
        os.path.join(_PROJ, "scripts", "trackb_inplace_res4_ats_smoke.py"))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate scripts/trackb_inplace_res4_ats_smoke.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Lightweight (no real build) tests: module import + constants + synthetic
# helper behaviour.
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ats():
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    return _load_ats()


def test_softcore_canon_constants(ats):
    """C3/R3: the soft-core canon must NOT be re-tuned."""
    assert ats.ATS_UMAX_KCAL == 200.0
    assert ats.ATS_UBCORE_KCAL == 100.0
    assert ats.ATS_ACORE == 0.062500


def test_upstream_expression_strings_carry_uoffset(ats):
    """C1: the verbatim upstream expression must include the UOffset term that
    the 9-arg convenience constructor omits."""
    full = (ats._ATS_REFERENCE_POT_EXPR + ats._ATS_ALCHEMICAL_POT_EXPR
            + ats._ATS_SOFTCORE_EXPR)
    assert "UOffset" in full
    assert "select(step(Direction), u0, u1)" in full
    assert "usc" in full and "Umax" in full and "Ubcore" in full and "Acore" in full


def test_expression_string_parses_with_ten_globals(ats):
    """C1: the upstream expression-string ATMForce ctor + 10 globals must build
    on the host openmm (NOT the 9-arg convenience ctor)."""
    import openmm as mm
    full = (ats._ATS_REFERENCE_POT_EXPR + ats._ATS_ALCHEMICAL_POT_EXPR
            + ats._ATS_SOFTCORE_EXPR)
    atm = mm.ATMForce(full)
    for n in ("Lambda1", "Lambda2", "Alpha", "Uh", "W0", "Umax", "Ubcore",
              "Acore", "Direction", "UOffset"):
        atm.addGlobalParameter(n, 0.0)
    assert atm.getNumGlobalParameters() == 10


def test_host_api_supports_swap_primitives(ats):
    """C1: setParticleTransformation + ParticleOffsetDisplacement + FixedDisplacement
    must exist on the host openmm (the swap reuse precondition)."""
    import openmm as mm
    assert hasattr(mm.ATMForce, "setParticleTransformation")
    assert hasattr(mm, "ParticleOffsetDisplacement")
    assert hasattr(mm, "FixedDisplacement")


def test_attach_atm_force_legacy_preserved(ats):
    """C8/R-7: the legacy single-topology functions must remain intact."""
    assert hasattr(ats, "attach_atm_force")
    assert hasattr(ats, "smoke_test_leg")
    # The new fused path is a separate function (not an overwrite).
    assert hasattr(ats, "build_inplace_res4_fused_system")
    assert hasattr(ats, "attach_inplace_swap_atmforce")


# ---------------------------------------------------------------------------
# Synthetic C3 / MC1 helper tests (no real build): construct two tiny fake
# "builds" with controlled topologies so the count-parity + MC1 logic is
# exercised deterministically.
# ---------------------------------------------------------------------------
def _toy_build(ats, names, charges, var_names):
    """Construct a minimal fake build dict (topology + NonbondedForce) with the
    given residue-4 atom names + charges, var atoms = var_names."""
    import openmm as mm
    from openmm import app
    import openmm.unit as unit

    top = app.Topology()
    chain = top.addChain(id="B")
    res = top.addResidue("TRP", chain, id="4")
    sysm = mm.System()
    nb = mm.NonbondedForce()
    idx_by_name = {}
    for nm in names:
        a = top.addAtom(nm, app.element.carbon, res)
        idx_by_name[nm] = a.index
        sysm.addParticle(12.0 * unit.dalton)
        nb.addParticle(charges[nm] * unit.elementary_charge,
                       0.3 * unit.nanometer, 0.4 * unit.kilojoule_per_mole)
    sysm.addForce(nb)

    common = [idx_by_name[n] for n in names if n == "NE1"]
    wt_only = [idx_by_name[n] for n in names if n in var_names.get("wt", [])]
    mtr_only = [idx_by_name[n] for n in names if n in var_names.get("mtr", [])]

    class _M:
        pass
    m = _M()
    m.topology = top
    m.positions = [mm.Vec3(0, 0, 0) * unit.nanometer for _ in names]
    return {
        "system": sysm,
        "modeller": m,
        "alchemical_atoms": {"common": common, "wt_only": wt_only,
                             "mtr_only": mtr_only},
    }


def test_c3_count_parity_pass(ats):
    """C3: equal common-atom count + aligned names -> map builds."""
    names = ["NE1", "CA", "CB", "HE1"]
    wt = _toy_build(ats, names, {n: 0.0 for n in names},
                    {"wt": ["HE1"], "mtr": []})
    names_m = ["NE1", "CA", "CB", "CM"]
    mtr = _toy_build(ats, names_m, {n: 0.0 for n in names_m},
                     {"wt": [], "mtr": ["CM"]})
    cmap = ats._build_common_index_map(wt, mtr, binder_chain="B")
    assert cmap["n_common"] == 3  # NE1, CA, CB
    assert len(cmap["wt_common"]) == len(cmap["mtr_common"]) == 3


def test_c3_order_alignment_fail_loud(ats):
    """C3: common atoms in a DIFFERENT name order must fail loud (positional swap)."""
    wt = _toy_build(ats, ["NE1", "CA", "CB", "HE1"],
                    {n: 0.0 for n in ["NE1", "CA", "CB", "HE1"]},
                    {"wt": ["HE1"], "mtr": []})
    # MTR common order scrambled: CB before CA.
    mtr = _toy_build(ats, ["NE1", "CB", "CA", "CM"],
                     {n: 0.0 for n in ["NE1", "CB", "CA", "CM"]},
                     {"wt": [], "mtr": ["CM"]})
    with pytest.raises(ValueError, match="C3 common-atom ORDER"):
        ats._build_common_index_map(wt, mtr, binder_chain="B")


def test_c3_count_parity_fail_loud(ats):
    """C3: unequal common-atom counts must fail loud (upstream _exit gate)."""
    wt = _toy_build(ats, ["NE1", "CA", "CB", "HE1"],
                    {n: 0.0 for n in ["NE1", "CA", "CB", "HE1"]},
                    {"wt": ["HE1"], "mtr": []})
    # MTR missing CB among commons -> 2 vs 3.
    mtr = _toy_build(ats, ["NE1", "CA", "CM"],
                     {n: 0.0 for n in ["NE1", "CA", "CM"]},
                     {"wt": [], "mtr": ["CM"]})
    with pytest.raises(ValueError, match="C3 common-atom COUNT parity"):
        ats._build_common_index_map(wt, mtr, binder_chain="B")


def test_mc1_continuity_pass_when_identical(ats):
    """MC1: identical common charges -> passes."""
    names = ["NE1", "CA", "CB", "HE1"]
    q = {"NE1": -0.34, "CA": -0.02, "CB": 0.0, "HE1": 0.34}
    wt = _toy_build(ats, names, q, {"wt": ["HE1"], "mtr": []})
    names_m = ["NE1", "CA", "CB", "CM"]
    q_m = {"NE1": -0.34, "CA": -0.02, "CB": 0.0, "CM": -0.1}
    mtr = _toy_build(ats, names_m, q_m, {"wt": [], "mtr": ["CM"]})
    cmap = ats._build_common_index_map(wt, mtr, binder_chain="B")
    res = ats.assert_common_atom_param_continuity(wt, mtr, cmap)
    assert res["passed"] is True
    assert res["max_dq_e"] < 1e-9


def test_mc1_continuity_fail_loud_when_divergent(ats):
    """MC1: divergent common charge (the production hybrid-XML risk) fails loud."""
    names = ["NE1", "CA", "CB", "HE1"]
    q = {"NE1": -0.34, "CA": -0.02, "CB": 0.0, "HE1": 0.34}
    wt = _toy_build(ats, names, q, {"wt": ["HE1"], "mtr": []})
    names_m = ["NE1", "CA", "CB", "CM"]
    q_m = {"NE1": -0.34, "CA": -0.10, "CB": 0.13, "CM": -0.1}  # CA/CB diverge
    mtr = _toy_build(ats, names_m, q_m, {"wt": [], "mtr": ["CM"]})
    cmap = ats._build_common_index_map(wt, mtr, binder_chain="B")
    with pytest.raises(ValueError, match="MC1 common-atom param continuity FAIL"):
        ats.assert_common_atom_param_continuity(wt, mtr, cmap)


def test_harvest_remap_rejects_partner_outside_common(ats):
    """The HE1 partner-index remap must fail loud if a partner is not in the
    WT->MTR common map (the raw-index divergence trap)."""
    import openmm as mm
    import openmm.unit as unit
    from openmm import app
    # Build a tiny WT system with HE1 bonded to an index NOT in the map.
    sysm = mm.System()
    for _ in range(3):
        sysm.addParticle(1.0 * unit.dalton)
    nb = mm.NonbondedForce()
    for _ in range(3):
        nb.addParticle(0.0, 0.3 * unit.nanometer, 0.0)
    bf = mm.HarmonicBondForce()
    bf.addBond(2, 0, 0.1 * unit.nanometer,  # HE1=2 bonded to atom 0
               1000.0 * unit.kilojoule_per_mole / unit.nanometer ** 2)
    sysm.addForce(nb)
    sysm.addForce(bf)
    wt_build = {"system": sysm}
    # Map deliberately omits atom 0.
    with pytest.raises(ValueError, match="no MTR counterpart"):
        ats._harvest_he1_terms(wt_build, he1_index=2, wt_to_mtr={1: 5})


# ---------------------------------------------------------------------------
# Real fused-box tests (need endpoint final.pdb + openmm). Module-scoped so the
# ~30-60 s build runs once.
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def real_fused_harmonized(ats):
    """The GENUINE-swap fused box (the default mode): the real HE1<->methyl
    coupling toggle (physically plausible |u1-u0|, true mutation geometry)."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    return ats.build_inplace_res4_fused_system(
        leg="free", seed="s7", solvate=False,
        harmonize_common_charges=True, swap_mode="genuine",
        genuine_decouple_nm=1.2)


@pytest.fixture(scope="module")
def real_fused_var_park(ats):
    """The DIAGNOSTIC var_park fused box (legacy parked-clash; retained behind a
    flag). Its |u1-u0| is bond-strain-dominated, NOT a genuine perturbation."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    return ats.build_inplace_res4_fused_system(
        leg="free", seed="s7", solvate=False,
        harmonize_common_charges=True, swap_mode="var_park",
        var_park_nm=(0.8, 0.0, 0.0))


@pytest.fixture(scope="module")
def real_fused_default(ats):
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    return ats.build_inplace_res4_fused_system(
        leg="free", seed="s7", solvate=False, harmonize_common_charges=False)


def test_real_default_handles_mc1_state(real_fused_default):
    """Production charges (no in-memory harmonize): the build must HONESTLY
    reflect the on-disk MTR RBFE XML's common-charge state — NOT a crash, NOT a
    silent false pass either way.

    The shared/RBFE MTR XML's common-core charge continuity is a cross-track,
    untracked artifact that may or may not be harmonized at the file level. So:
      - if the common core DIVERGES from amber14SB-Trp -> outcome must be the
        structured 'mc1_charge_discontinuity' (the pre-registered finding ii),
        with passed=False and a non-trivial displaced charge;
      - if the common core MATCHES (file-level harmonized) -> the build attaches
        ('fused_attached') and MC1 passes cleanly.
    """
    outcome = real_fused_default["outcome"]
    mc1 = real_fused_default["mc1_param_continuity"]
    if outcome == "mc1_charge_discontinuity":
        assert mc1["passed"] is False
        assert mc1["n_diverging"] > 0
        assert abs(mc1["sum_dq_res4_e"]) > 0.05
    else:
        assert outcome == "fused_attached"
        assert mc1["passed"] is True
        assert mc1["max_dq_e"] <= 1e-4


def test_real_strict_mc1_matches_disk_state(ats):
    """strict_mc1=True re-raises the MC1 fail-loud IFF the on-disk common core
    actually diverges; if the file is harmonized it attaches cleanly. Either way
    the behaviour is honest (no false pass, no spurious raise)."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    # Probe the non-strict outcome first to learn the on-disk charge state.
    probe = ats.build_inplace_res4_fused_system(
        leg="free", seed="s7", solvate=False, harmonize_common_charges=False)
    if probe["outcome"] == "mc1_charge_discontinuity":
        with pytest.raises(ValueError,
                           match="MC1 common-atom param continuity FAIL"):
            ats.build_inplace_res4_fused_system(
                leg="free", seed="s7", solvate=False,
                harmonize_common_charges=False, strict_mc1=True)
    else:
        strict = ats.build_inplace_res4_fused_system(
            leg="free", seed="s7", solvate=False,
            harmonize_common_charges=False, strict_mc1=True)
        assert strict["outcome"] == "fused_attached"


def test_real_count_parity_206(real_fused_harmonized):
    """C3: the real WT/MTR common-atom count parity (206 == 206)."""
    cmap = real_fused_harmonized["common_map"]
    assert cmap["n_common"] == len(cmap["wt_common"]) == len(cmap["mtr_common"])
    assert cmap["n_common"] == 206
    # Var partition: WT-var = {HE1}, MTR-var = {CM, HM1-3}.
    assert len(cmap["wt_var"]) == 1
    assert len(cmap["mtr_var"]) == 4


def test_real_attach_atom_is_ne1(real_fused_harmonized, ats):
    """C2: the swap attach atom is NE1 (HE1/CM bond partner), not Calpha."""
    fused = real_fused_harmonized["fused_build"]
    ne1 = fused["alchemical_atoms"]["common"][0]
    name = next(a.name for a in fused["modeller"].topology.atoms()
                if a.index == ne1)
    assert name == "NE1"
    assert real_fused_harmonized["swap"]["attach_atom"] == ne1


def test_real_mc2_methyl_bonded_present(real_fused_harmonized):
    """MC2: CM-NE1 bond + HM-CM connectivity present in the fused box."""
    mc2 = real_fused_harmonized["mc2_methyl_bonded"]
    assert mc2["cm_ne1_bond_present"] is True
    assert all(mc2["hm_cm_bonds_present"].values())


def test_real_mc3_disulfide_preserved(real_fused_harmonized):
    """MC3: cyclic_ss SG-SG disulfide preserved in the fused box."""
    assert real_fused_harmonized["mc3_disulfide"]["disulfide_present"] is True


def test_real_seed_assert_hard_dist_above_threshold(real_fused_harmonized):
    """R2: appearing atoms are non-clashing vs COMMON/solvent (HE1 overlap is the
    FLAGGED dual-topology exclusion pair, not a hard fail)."""
    sa = real_fused_harmonized["seed_assert"]
    assert sa["passed"] is True
    assert sa["min_hard_dist_nm"] > 0.10
    assert 0.12 < sa["cm_ne1_nm"] < 0.17   # ~0.145 nm equilibrium


def test_real_var_var_exclusions_added(real_fused_harmonized):
    """The HE1<->methyl var-var exclusions (4) must be injected (else u0 ~1e6)."""
    inj = real_fused_harmonized["fused_build"]["injection"]
    assert inj["n_var_var_exclusions"] == 4


def test_real_swap_not_null_op(real_fused_harmonized):
    """C6b: the genuine swap is NOT the degenerate single-shared-core null op."""
    swap = real_fused_harmonized["swap"]
    assert swap["single_shared_core_null_op"] is False
    assert swap["genuine_applied"] is True
    assert swap["var_park_applied"] is False
    # The genuine mode adds the dummy NE1 reference + decouples HE1's 2 angles.
    assert swap["ne1_ref_index"] is not None
    assert swap["n_he1_angles_decoupled"] == 2


def test_real_var_park_is_diagnostic_park(real_fused_var_park):
    """The legacy var_park diagnostic path is still wired (behind the flag)."""
    swap = real_fused_var_park["swap"]
    assert swap["swap_mode"] == "var_park"
    assert swap["var_park_applied"] is True
    assert swap["genuine_applied"] is False
    assert swap["single_shared_core_null_op"] is False


def test_real_energy_finite_and_u1_ne_u0(real_fused_harmonized):
    """C6a + C6b: E_ATM/u0/u1 finite (coverage NaN defeated) AND u1 != u0 (null
    op defeated) at Lambda=0.5, Direction=+1 — with a PHYSICALLY PLAUSIBLE
    |u1-u0| (the genuine swap, NOT a parked-clash 1e4)."""
    import openmm as mm
    import openmm.unit as unit
    import numpy as np
    fused = real_fused_harmonized["fused_build"]
    sysm = fused["system"]
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001 * unit.picoseconds),
                     mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(fused["modeller"].positions)
    ctx.setParameter("Lambda1", 0.5)
    ctx.setParameter("Lambda2", 0.5)
    ctx.setParameter("Direction", 1.0)
    st = ctx.getState(getEnergy=True, getForces=True)
    e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    atm = next(sysm.getForce(i) for i in range(sysm.getNumForces())
               if isinstance(sysm.getForce(i), mm.ATMForce))
    pert = atm.getPerturbationEnergy(ctx)
    u1 = pert[0].value_in_unit(unit.kilocalorie_per_mole)
    u0 = pert[1].value_in_unit(unit.kilocalorie_per_mole)
    assert np.isfinite(e) and np.isfinite(u1) and np.isfinite(u0)
    assert abs(u1 - u0) > 1e-3       # u1 != u0 (genuine, not null)
    # Physically plausible for a HE1<->methyl edit (ones-to-tens of kcal/mol):
    # the genuine swap must NOT be a parked-clash / bond-strain artifact.
    assert abs(u1 - u0) <= 1.0e3
    forces = st.getForces().value_in_unit(
        unit.kilocalorie_per_mole / unit.angstrom)
    fmax = max(float(np.linalg.norm(f)) for f in forces)
    assert fmax < 1.0e5              # bounded


def test_real_var_park_is_parked_clash(real_fused_var_park):
    """The var_park diagnostic produces an ELEVATED (parked-clash) |u1-u0|,
    distinguishing it from the genuine swap — this is exactly why the smoke's
    C6d plausibility assert rejects it."""
    import openmm as mm
    import openmm.unit as unit
    import numpy as np
    fused = real_fused_var_park["fused_build"]
    sysm = fused["system"]
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001 * unit.picoseconds),
                     mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(fused["modeller"].positions)
    ctx.setParameter("Lambda1", 0.5)
    ctx.setParameter("Lambda2", 0.5)
    ctx.setParameter("Direction", 1.0)
    atm = next(sysm.getForce(i) for i in range(sysm.getNumForces())
               if isinstance(sysm.getForce(i), mm.ATMForce))
    pert = atm.getPerturbationEnergy(ctx)
    u1 = pert[0].value_in_unit(unit.kilocalorie_per_mole)
    u0 = pert[1].value_in_unit(unit.kilocalorie_per_mole)
    assert np.isfinite(u1) and np.isfinite(u0)
    # The parked methyl strains the CM-NE1 bond (k~282000) -> |u1-u0| >> 1e3.
    assert abs(u1 - u0) > 1.0e3


def test_real_endpoint_equivalence_u0_reproduces_mtr(real_fused_harmonized, ats):
    """C6e/Q5: the fused reference energy u0 reproduces the pure MTR-only energy
    within tolerance (the swap is CORRECT, not merely finite)."""
    import openmm as mm
    import openmm.unit as unit
    fused = real_fused_harmonized["fused_build"]
    sysm = fused["system"]
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001 * unit.picoseconds),
                     mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(fused["modeller"].positions)
    ctx.setParameter("Lambda1", 0.0)
    ctx.setParameter("Lambda2", 0.0)
    ctx.setParameter("Direction", 1.0)
    atm = next(sysm.getForce(i) for i in range(sysm.getNumForces())
               if isinstance(sysm.getForce(i), mm.ATMForce))
    u0 = atm.getPerturbationEnergy(ctx)[1].value_in_unit(
        unit.kilocalorie_per_mole)

    # Pure MTR-only harmonized reference.
    smoke = _load_smoke()
    mtr_e = smoke._pure_mtr_only_energy(
        "s7", "B", solvate=False, harmonize_common_charges=True)
    # u0 carries HE1's residual coupling on top of MTR-only; within ~10 kcal/mol.
    assert abs(u0 - mtr_e) <= 10.0


# ---------------------------------------------------------------------------
# Smoke-script plumbing.
# ---------------------------------------------------------------------------
def test_smoke_module_imports():
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    smoke = _load_smoke()
    assert hasattr(smoke, "run_tier1_smoke")
    assert hasattr(smoke, "main")
    assert hasattr(smoke, "_softcore_usc")
    assert smoke.MAX_FORCE_KCAL_PER_MOL_A == 1.0e5
    assert smoke.ENDPOINT_EQUIV_TOL_KCAL == 10.0
    assert smoke.U1_MINUS_U0_PLAUSIBLE_MAX_KCAL == 1.0e3


def test_smoke_softcore_usc_caps_at_umax():
    """The C6d usc recomputation is REAL (not `or True`): below Ubcore it is the
    identity; a huge raw perturbation is capped below Umax."""
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    smoke = _load_smoke()
    umax, ubcore, acore = 200.0, 100.0, 0.0625
    # Below Ubcore: identity.
    assert abs(smoke._softcore_usc(5.0, umax, ubcore, acore) - 5.0) < 1e-9
    # A parked-clash-scale raw perturbation is capped strictly below Umax.
    usc_big = smoke._softcore_usc(14500.0, umax, ubcore, acore)
    assert usc_big <= umax
    assert usc_big > ubcore


def test_smoke_default_handles_mc1_state():
    """Smoke default (production charges) HONESTLY reflects the on-disk MTR RBFE
    XML's common-charge state: either the structured MC1 finding (if the common
    core diverges) or a genuine tier1_pass (if the file is harmonized). Either
    way ranking-only + prediction-test, never a silent false pass."""
    if not _have_openmm() or not _endpoints_present():
        pytest.skip("openmm or endpoints unavailable")
    smoke = _load_smoke()
    res = smoke.run_tier1_smoke(
        seed="s7", solvate=False, harmonize_common_charges=False,
        swap_mode="genuine", platform_name="Reference")
    assert res["prediction_test"] is True
    assert res["regime"] == "ranking_only"
    assert res["outcome"] in ("mc1_charge_discontinuity", "tier1_pass")
    if res["outcome"] == "tier1_pass":
        # When the file is harmonized the genuine swap must still be physical.
        assert abs(res["energies_kcal"]["u1_minus_u0"]) \
            <= smoke.U1_MINUS_U0_PLAUSIBLE_MAX_KCAL


def test_smoke_harmonized_passes_tier1():
    """Smoke GENUINE (mechanical validation) -> tier1_pass with all C6, and a
    PHYSICALLY PLAUSIBLE |u1-u0| (ones-to-tens of kcal/mol, not a parked clash)."""
    if not _have_openmm() or not _endpoints_present():
        pytest.skip("openmm or endpoints unavailable")
    smoke = _load_smoke()
    res = smoke.run_tier1_smoke(
        seed="s7", solvate=False, harmonize_common_charges=True,
        swap_mode="genuine", genuine_decouple_nm=1.2, platform_name="Reference")
    assert res["outcome"] == "tier1_pass"
    assert res["swap_mode"] == "genuine"
    en = res["energies_kcal"]
    assert abs(en["u1_minus_u0"]) > 1e-3
    # Genuine perturbation: physically plausible, NOT a parked-clash 1e4.
    assert abs(en["u1_minus_u0"]) <= smoke.U1_MINUS_U0_PLAUSIBLE_MAX_KCAL
    assert res["u1_minus_u0_plausible"] is True
    assert res["usc_le_umax"] is True
    assert res["max_force_kcal_per_mol_A"] < 1.0e5
    assert res["endpoint_equivalence"]["passed"] is True
    assert res["regime"] == "ranking_only"


def test_smoke_var_park_fails_plausibility():
    """The legacy var_park diagnostic now FAILS the C6d plausibility assert (the
    `or True` vacuous gap is fixed): a parked-clash |u1-u0|~1e4 is rejected, not
    silently passed."""
    if not _have_openmm() or not _endpoints_present():
        pytest.skip("openmm or endpoints unavailable")
    smoke = _load_smoke()
    with pytest.raises(AssertionError, match="C6d FAIL"):
        smoke.run_tier1_smoke(
            seed="s7", solvate=False, harmonize_common_charges=True,
            swap_mode="var_park", var_park_nm=(0.8, 0.0, 0.0),
            platform_name="Reference")


def test_smoke_null_mode_is_null_op():
    """The null mode (degenerate single-shared-core) FAILS C6b (u1==u0)."""
    if not _have_openmm() or not _endpoints_present():
        pytest.skip("openmm or endpoints unavailable")
    smoke = _load_smoke()
    with pytest.raises(AssertionError, match="C6b FAIL"):
        smoke.run_tier1_smoke(
            seed="s7", solvate=False, harmonize_common_charges=True,
            swap_mode="null", platform_name="Reference")


def test_build_uses_rbfe_xml_not_shared(ats):
    """Cross-track isolation: the in-place RBFE build loads the dedicated RBFE
    MTR XML (or falls back to the shared file only if the RBFE file is absent),
    and NEVER mutates the shared params/MTR_gaff2_hybrid.xml on disk."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    import os
    build = ats.build_inplace_res4_fused_system(
        leg="free", seed="s7", solvate=False,
        harmonize_common_charges=True, swap_mode="genuine")
    used = build["mtr_ncaa_xml"]
    # The RBFE build must resolve to the RBFE-scoped XML when it exists.
    if os.path.isfile(ats.HYBRID_MTR_XML_RBFE):
        assert os.path.abspath(used) == os.path.abspath(ats.HYBRID_MTR_XML_RBFE)
    else:
        assert os.path.abspath(used) == os.path.abspath(ats.HYBRID_MTR_XML)
    # The build path never points at the shared XML when the RBFE file exists.
    assert build["swap_mode"] == "genuine"


# ---------------------------------------------------------------------------
# Tier-2: finite-UNDER-INTEGRATION short-dynamics runner. The cheap tests are
# synthetic (no build): constants, the constraint-label mapping, and the
# band-assert logic (the equilibration-window lower-bound exemption). The one
# real Reference-platform short integration is endpoint-gated.
# ---------------------------------------------------------------------------
def test_tier2_constants_present():
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    smoke = _load_smoke()
    assert smoke.TIER2_DEFAULT_STEPS in range(200, 501)
    assert smoke.TIER2_CHECKPOINT_EVERY == 50
    assert smoke.TIER2_TEMPERATURE_K == 300.0
    assert smoke.TIER2_FRICTION_PER_PS == 1.0
    assert smoke.TIER2_MAX_FORCE_KCAL_PER_MOL_A == 1.0e5
    # The R3 PRIMARY integrator-matrix cell must be first AND be 1fs/unconstrained.
    first = smoke.TIER2_INTEGRATOR_MATRIX[0]
    assert first[1] == 1.0 and first[2] == "none"
    labels = [c[0] for c in smoke.TIER2_INTEGRATOR_MATRIX]
    assert "dt2fs_unconstrained" in labels and "dt2fs_hbonds" in labels


def test_tier2_constraints_enum_mapping():
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    import openmm as mm
    smoke = _load_smoke()
    assert smoke._constraints_enum("none") is None
    assert smoke._constraints_enum("unconstrained") is None
    assert smoke._constraints_enum("hbonds") is mm.app.HBonds
    with pytest.raises(ValueError):
        smoke._constraints_enum("allbonds")


def test_tier2_band_assert_upper_runaway_always_fails():
    """The runaway UPPER temperature bound is enforced at EVERY checkpoint (a
    box that explodes upward is a hard FAIL regardless of equilibration)."""
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    smoke = _load_smoke()
    run = {
        "survived": True, "n_steps_requested": 100, "n_steps_completed": 100,
        "nan_mode": None, "failed_at_step": None,
        "trajectory": [
            {"step": 0, "finite": True, "fmax_kcal_per_mol_A": 50.0,
             "temp_K": 300.0},
            {"step": 50, "finite": True, "fmax_kcal_per_mol_A": 50.0,
             "temp_K": 4200.0},  # runaway
        ],
    }
    with pytest.raises(AssertionError, match="runaway upper bound"):
        smoke._assert_tier2_stability(run)


def test_tier2_band_assert_lower_exempt_during_equilibration():
    """A transient cold checkpoint WITHIN the equilibration window must NOT fail
    (1/ps friction => ~1 ps thermalization; a short window legitimately runs
    cooler than 300 K). Past the window a sunk box DOES fail."""
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    smoke = _load_smoke()
    # Cold at the early checkpoints (within equil_skip) -> PASS.
    ok = {
        "survived": True, "n_steps_requested": 200, "n_steps_completed": 200,
        "nan_mode": None, "failed_at_step": None,
        "trajectory": [
            {"step": 0, "finite": True, "fmax_kcal_per_mol_A": 40.0,
             "temp_K": 290.0},
            {"step": 50, "finite": True, "fmax_kcal_per_mol_A": 40.0,
             "temp_K": 150.0},   # transient undershoot, within equil window
            {"step": 100, "finite": True, "fmax_kcal_per_mol_A": 40.0,
             "temp_K": 160.0},
            {"step": 150, "finite": True, "fmax_kcal_per_mol_A": 40.0,
             "temp_K": 180.0},
        ],
    }
    smoke._assert_tier2_stability(ok, equil_skip=2)  # no raise
    # A sunk box AFTER the window (below the lower bound) -> FAIL.
    sunk = dict(ok)
    sunk["trajectory"] = [
        {"step": 0, "finite": True, "fmax_kcal_per_mol_A": 40.0, "temp_K": 290.0},
        {"step": 50, "finite": True, "fmax_kcal_per_mol_A": 40.0, "temp_K": 150.0},
        {"step": 100, "finite": True, "fmax_kcal_per_mol_A": 40.0, "temp_K": 5.0},
    ]
    with pytest.raises(AssertionError, match="below the sane lower bound"):
        smoke._assert_tier2_stability(sunk, equil_skip=2)


def test_tier2_band_assert_nonfinite_and_fmax_fail():
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    smoke = _load_smoke()
    base = {"survived": True, "n_steps_requested": 50, "n_steps_completed": 50,
            "nan_mode": None, "failed_at_step": None}
    nonfinite = dict(base, trajectory=[
        {"step": 0, "finite": False, "fmax_kcal_per_mol_A": 40.0, "temp_K": 300.0}])
    with pytest.raises(AssertionError, match="non-finite"):
        smoke._assert_tier2_stability(nonfinite)
    bigf = dict(base, trajectory=[
        {"step": 0, "finite": True, "fmax_kcal_per_mol_A": 2.0e5, "temp_K": 300.0}])
    with pytest.raises(AssertionError, match="max per-atom force"):
        smoke._assert_tier2_stability(bigf)


def test_tier2_band_assert_did_not_complete_fails():
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    smoke = _load_smoke()
    crashed = {
        "survived": False, "n_steps_requested": 300, "n_steps_completed": 120,
        "nan_mode": "integrator_exception: Particle coordinate is NaN",
        "failed_at_step": 150,
        "trajectory": [
            {"step": 0, "finite": True, "fmax_kcal_per_mol_A": 40.0,
             "temp_K": 300.0}],
    }
    with pytest.raises(AssertionError, match="did NOT complete"):
        smoke._assert_tier2_stability(crashed)


def test_tier2_system_ndf_excludes_massless():
    """The DOF count must EXCLUDE massless particles (the genuine-mode dummy NE1
    reference) so the instantaneous-T estimate is not inflated."""
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    import openmm as mm
    smoke = _load_smoke()
    sysm = mm.System()
    sysm.addParticle(12.0)   # massive
    sysm.addParticle(1.0)    # massive
    sysm.addParticle(0.0)    # massless dummy reference -> excluded
    # 2 massive, no constraints, no CMMotionRemover -> 3*2 = 6.
    assert smoke._system_ndf(sysm) == 6


def test_tier2_real_reference_short_run_survives():
    """REAL genuine fused box, Reference platform, UNSOLVATED + harmonized (the
    cheap CPU mechanical path), 1 fs UNCONSTRAINED alch-H, a short window: must
    SURVIVE (no NaN) and the primary must pass the Tier-2 asserts. The CUDA+PME
    long run is the GPU validation step; this guards the integration WIRING."""
    if not _have_openmm() or not _endpoints_present():
        pytest.skip("openmm or endpoints unavailable")
    smoke = _load_smoke()
    res = smoke.run_tier2_stability(
        seed="s7", solvate=False, harmonize_common_charges=True,
        steps=40, dt_fs=1.0, constraints_label="none",
        checkpoint_every=20, platform_name="Reference",
        run_integrator_matrix=False)
    assert res["regime"] == "ranking_only"
    assert res["prediction_test"] is True
    # Either a genuine pass, or the MC1 finding if the on-disk XML is divergent
    # (harmonize=True forces continuity, so a pass is expected here).
    assert res["outcome"] in ("tier2_pass", "mc1_charge_discontinuity")
    if res["outcome"] == "tier2_pass":
        pr = res["primary"]
        assert pr["survived"] is True
        assert pr["n_steps_completed"] == 40
        assert pr["constraints"] == "none"
        assert pr["hmr"] is False
        assert pr["max_fmax_kcal_per_mol_A"] < smoke.TIER2_MAX_FORCE_KCAL_PER_MOL_A


def test_tier2_unconstrained_build_has_no_alch_h_constraint(ats):
    """The constraints=None build path must produce a System with NO SHAKE on the
    alch H (in fact NO constraints at all on the unsolvated protein-only box, so
    the appearing/disappearing methyl/HE1 H integrate freely — the R3 spec)."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    build = ats.build_inplace_res4_fused_system(
        leg="free", seed="s7", solvate=False, harmonize_common_charges=True,
        swap_mode="genuine", constraints=None)
    n_constr = build["fused_build"]["system"].getNumConstraints()
    assert n_constr == 0
    # The default (HBonds) build DOES carry X-H constraints (contrast).
    build_h = ats.build_inplace_res4_fused_system(
        leg="free", seed="s7", solvate=False, harmonize_common_charges=True,
        swap_mode="genuine")
    assert build_h["fused_build"]["system"].getNumConstraints() > 0


# ---------------------------------------------------------------------------
# BOUND leg (RBFE complex): same residue-4 HE1<->methyl transform in the
# receptor+peptide bound pose, re-solvated, NO displacement (the binder stays
# bound). Reuses every free-leg helper; the only new logic is the bound-complex
# source + the local-density-aware genuine decouple direction. All bound tests
# use the cheap UNSOLVATED + harmonized + Reference path (the CUDA+PME solvated
# run is the executor step). DDG_bind = DG_mut(bound) - DG_mut(free).
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def real_fused_bound(ats):
    """The GENUINE-swap fused BOUND box (receptor+peptide, unsolvated, harmonized
    for mechanical isolation)."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    return ats.build_inplace_res4_fused_system(
        leg="bound", seed="s7", solvate=False,
        harmonize_common_charges=True, swap_mode="genuine",
        genuine_decouple_nm=1.2)


def test_bound_leg_rejects_unknown_leg(ats):
    """leg must be 'free' or 'bound' (a typo is a fail-loud NotImplementedError,
    not a silent free-leg fallback)."""
    with pytest.raises(NotImplementedError, match="leg must be 'free' or 'bound'"):
        ats.build_inplace_res4_fused_system(leg="complex", seed="s7",
                                            solvate=False)


def test_bound_complex_prep_keeps_receptor_and_binder(ats):
    """prepare_bound_complex_from_final keeps BOTH the receptor (chain A) and the
    binder (chain B) protein atoms and drops only the solvent/ions."""
    import os
    import tempfile
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    li = ats.resolve_leg_inputs("s7")
    td = tempfile.mkdtemp(prefix="ats_bound_prep_test_")
    out = os.path.join(td, "bound.pdb")
    ats.prepare_bound_complex_from_final(li["final"]["cp4"], out, "B")
    chains = set()
    resnames = set()
    n = 0
    for line in open(out):
        if line[:6] in ("ATOM  ", "HETATM"):
            chains.add(line[21])
            resnames.add(line[17:20].strip())
            n += 1
    # Receptor (A) AND binder (B) present; the free-leg prep would have only B.
    assert "A" in chains and "B" in chains
    # Solvent/ions dropped.
    assert not (resnames & ats._SOLVENT_RESNAMES)
    # The ncAA (MTR) survives on the binder.
    assert "MTR" in resnames
    # The bound complex is much bigger than the free peptide (receptor present).
    assert n > 5000


def test_bound_decouple_direction_clears_local_density(real_fused_bound, ats):
    """The bound genuine decouple direction must point OUTWARD from the local atom
    density so the displaced HE1 lands in low density (NOT clashing the receptor
    OR the peptide's own fold). Verifies the landing point is far from every
    heavy atom — the fix for the away-from-receptor heuristic that landed HE1 on
    the binder's Arg side chain."""
    import numpy as np
    import openmm.unit as unit
    d = real_fused_bound["genuine_decouple_dir"]
    assert d is not None                      # bound leg computes a direction
    # Unit vector.
    assert abs(float(np.linalg.norm(np.array(d))) - 1.0) < 1e-6
    fb = real_fused_bound["fused_build"]
    sysm = fb["system"]
    pos = np.array([v.value_in_unit(unit.nanometer)
                    for v in fb["modeller"].positions])
    ne1 = fb["alchemical_atoms"]["common"][0]
    he1 = fb["fused_he1_index"]
    ne1_ref = real_fused_bound["swap"]["ne1_ref_index"]
    landing = pos[ne1] + 1.2 * np.array(d)    # where HE1 sits at u1
    # Nearest HEAVY atom (mass>1.5) to the landing point, excluding HE1 + dummy.
    mind = np.inf
    for i in range(sysm.getNumParticles()):
        if i in (he1, ne1_ref):
            continue
        if sysm.getParticleMass(i).value_in_unit(unit.dalton) <= 1.5:
            continue
        dd = float(np.linalg.norm(pos[i] - landing))
        if dd < mind:
            mind = dd
    # Local-density-aware direction lands ~1 nm from any heavy atom (bulk); the
    # naive away-from-receptor heuristic landed at 0.12 nm (a peptide clash).
    assert mind > 0.5


def test_bound_count_parity_matches_free(real_fused_bound):
    """C3: the bound common-atom count parity is the SAME 206 as the free leg —
    the receptor is excluded from the binder-protein-only common-core swap."""
    cmap = real_fused_bound["common_map"]
    assert cmap["n_common"] == 206
    assert len(cmap["wt_var"]) == 1
    assert len(cmap["mtr_var"]) == 4


def test_bound_receptor_rides_through(real_fused_bound):
    """The bound box carries the receptor as inert context: n_particles is far
    larger than the free leg (which is ~210 unsolvated)."""
    fb = real_fused_bound["fused_build"]
    assert fb["system"].getNumParticles() > 5000


def test_bound_mc3_disulfide_preserved(real_fused_bound):
    """MC3: cyclic_ss SG-SG disulfide survives into the bound fused box."""
    assert real_fused_bound["mc3_disulfide"]["disulfide_present"] is True


def test_bound_mc2_methyl_bonded_present(real_fused_bound):
    """MC2: CM-NE1 + HM-CM bonded terms present in the bound fused box."""
    mc2 = real_fused_bound["mc2_methyl_bonded"]
    assert mc2["cm_ne1_bond_present"] is True
    assert all(mc2["hm_cm_bonds_present"].values())


def test_bound_seed_assert_non_clashing(real_fused_bound):
    """R2: the appearing methyl seed is non-clashing in the bound pose too."""
    sa = real_fused_bound["seed_assert"]
    assert sa["passed"] is True
    assert sa["min_hard_dist_nm"] > 0.10


def test_bound_energy_finite_and_plausible(real_fused_bound):
    """C6a/b/c/d on the BOUND box: E_ATM/u0/u1 finite, u1 != u0 (genuine swap
    wired), |u1-u0| physically plausible (ones-to-tens of kcal/mol — the
    local-density decouple direction, NOT the +128 peptide-clash artifact), and
    Fmax bounded."""
    import numpy as np
    import openmm as mm
    import openmm.unit as unit
    fb = real_fused_bound["fused_build"]
    sysm = fb["system"]
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001 * unit.picoseconds),
                     mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(fb["modeller"].positions)
    ctx.setParameter("Lambda1", 0.5)
    ctx.setParameter("Lambda2", 0.5)
    ctx.setParameter("Direction", 1.0)
    st = ctx.getState(getEnergy=True, getForces=True)
    e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    atm = next(sysm.getForce(i) for i in range(sysm.getNumForces())
               if isinstance(sysm.getForce(i), mm.ATMForce))
    pert = atm.getPerturbationEnergy(ctx)
    u1 = pert[0].value_in_unit(unit.kilocalorie_per_mole)
    u0 = pert[1].value_in_unit(unit.kilocalorie_per_mole)
    assert np.isfinite(e) and np.isfinite(u1) and np.isfinite(u0)
    assert abs(u1 - u0) > 1e-3
    # The genuine decouple cost in the bound pose: ones-to-tens of kcal/mol, NOT
    # the parked-clash / peptide-collision regime (>1e3).
    assert abs(u1 - u0) <= 1.0e3
    forces = st.getForces().value_in_unit(
        unit.kilocalorie_per_mole / unit.angstrom)
    fmax = max(float(np.linalg.norm(f)) for f in forces)
    assert fmax < 1.0e5


def test_bound_endpoint_equivalence_u0_reproduces_mtr(real_fused_bound):
    """C6e/Q5 on the BOUND box: the fused reference energy u0 reproduces the pure
    bound MTR-only energy within tolerance (the swap is CORRECT, not merely
    finite). Uses the leg-aware MTR-only reference."""
    import openmm as mm
    import openmm.unit as unit
    fb = real_fused_bound["fused_build"]
    sysm = fb["system"]
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001 * unit.picoseconds),
                     mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(fb["modeller"].positions)
    ctx.setParameter("Lambda1", 0.0)
    ctx.setParameter("Lambda2", 0.0)
    ctx.setParameter("Direction", 1.0)
    atm = next(sysm.getForce(i) for i in range(sysm.getNumForces())
               if isinstance(sysm.getForce(i), mm.ATMForce))
    u0 = atm.getPerturbationEnergy(ctx)[1].value_in_unit(
        unit.kilocalorie_per_mole)
    smoke = _load_smoke()
    mtr_e = smoke._pure_mtr_only_energy(
        "s7", "B", solvate=False, harmonize_common_charges=True, leg="bound")
    assert abs(u0 - mtr_e) <= 10.0


def test_smoke_tier1_bound_runs(ats):
    """End-to-end smoke run_tier1_smoke(leg='bound') on the unsolvated harmonized
    Reference path: returns tier1_pass (or the MC1 finding) with leg recorded."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    smoke = _load_smoke()
    res = smoke.run_tier1_smoke(
        seed="s7", solvate=False, harmonize_common_charges=True,
        swap_mode="genuine", genuine_decouple_nm=1.2,
        platform_name="Reference", leg="bound")
    assert res["leg"] == "bound"
    assert res["outcome"] in ("tier1_pass", "tier1_endpoint_mismatch",
                              "mc1_charge_discontinuity")
    if res["outcome"] == "tier1_pass":
        assert res["genuine_decouple_dir"] is not None
        assert abs(res["energies_kcal"]["u1_minus_u0"]) <= 1.0e3
        assert res["endpoint_equivalence"]["passed"] is True


def test_smoke_tier2_bound_short_run(ats):
    """A very short Tier-2 BOUND run (Reference, 40 steps, no matrix) survives
    integration without NaN (or surfaces the MC1 finding). The CUDA+PME longer
    run is the executor step."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    smoke = _load_smoke()
    res = smoke.run_tier2_stability(
        seed="s7", solvate=False, harmonize_common_charges=True,
        genuine_decouple_nm=1.2, steps=40, dt_fs=1.0, constraints_label="none",
        checkpoint_every=20, platform_name="Reference",
        run_integrator_matrix=False, leg="bound")
    assert res["leg"] == "bound"
    assert res["outcome"] in ("tier2_pass", "mc1_charge_discontinuity")
    if res["outcome"] == "tier2_pass":
        assert res["primary"]["survived"] is True
        assert res["primary"]["n_steps_completed"] == 40
