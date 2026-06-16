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
    # Bound tol 30 (the free leg uses 10): in the contacting bound pose the WT
    # partner atom HE1 carries a genuine HE1<->receptor nonbonded residual on top
    # of the MTR-only energy -- the bound-only binding contribution, absent in the
    # free leg (no receptor). Single-H single-shell contact is a tens-of-kcal
    # band; 30 keeps ~1.5x margin over the observed ~20 while still flagging a
    # collapse / wrong-swap-wiring (those are ~10x away, z >> 3).
    assert abs(u0 - mtr_e) <= 30.0


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


# ===========================================================================
# CANONICAL ATS TWO-COPY rebuild tests
# (canonical two-copy ATS construction; C1-C10 design criteria).
#
# These validate the NEW two-copy path is canonical (copy-2 d-displaced into bulk,
# common-coord swap, DISTINCT attach atoms, ZERO inter-copy exclusions) and that
# the LEGACY single-shared-core path is preserved byte-for-byte (R-7/C9).
# ===========================================================================
def test_twocopy_functions_present(ats):
    """C9: the new two-copy path exists AND the legacy single-shared-core path is
    preserved (R-7 — densify_pilot and other callers still use the old path)."""
    for fn in ("build_inplace_res4_twocopy_system",
               "attach_twocopy_swap_atmforce",
               "_build_twocopy_index_map",
               "assert_twocopy_common_param_continuity",
               "assert_twocopy_separation",
               "assert_twocopy_methyl_bonded",
               "assert_twocopy_disulfides",
               "assert_twocopy_seed",
               "compute_twocopy_displacement_vector",
               "auto_search_twocopy_displacement",
               "_candidate_directions",
               "_copy_solute_heavy_positions",
               "_box_lengths_nm_from_vectors",
               "_min_image_min_distance_nm",
               "check_twocopy_endpoint_equivalence",
               "_register_copy2_common_to_copy1"):
        assert hasattr(ats, fn), "two-copy fn missing: %s" % fn
    # LEGACY path preserved (R-7/C9).
    assert hasattr(ats, "build_inplace_res4_fused_system")
    assert hasattr(ats, "attach_inplace_swap_atmforce")


def test_twocopy_displacement_constant(ats):
    """C2: the ATS peptide displacement convention (40 A = 4.0 nm)."""
    assert ats.ATS_TWOCOPY_DISPLACEMENT_NM == 4.0


def test_twocopy_shares_softcore_canon(ats):
    """C7: the two-copy path reuses the SAME soft-core canon (NOT re-tuned)."""
    # attach_twocopy_swap_atmforce defaults must be the module canon.
    import inspect
    sig = inspect.signature(ats.attach_twocopy_swap_atmforce)
    assert sig.parameters["umax_kcal"].default == ats.ATS_UMAX_KCAL
    assert sig.parameters["ubcore_kcal"].default == ats.ATS_UBCORE_KCAL
    assert sig.parameters["acore"].default == ats.ATS_ACORE


def test_twocopy_displace_positions_translates(ats):
    """_displace_copy_positions translates by the d-vector exactly."""
    import openmm as mm
    import openmm.unit as unit
    pos = [mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer,
           mm.Vec3(1.0, 2.0, 3.0) * unit.nanometer]
    out = ats._displace_copy_positions(pos, (4.0, 0.0, 0.0))
    a = out[0].value_in_unit(unit.nanometer)
    b = out[1].value_in_unit(unit.nanometer)
    assert abs(a[0] - 4.0) < 1e-9 and abs(a[1]) < 1e-9
    assert abs(b[0] - 5.0) < 1e-9 and abs(b[1] - 2.0) < 1e-9


# --- Real two-copy build (need endpoint final.pdb + openmm). Module-scoped. ---
@pytest.fixture(scope="module")
def real_twocopy_unsolv(ats):
    """The canonical two-copy box, UNSOLVATED (cheap CPU), with the on-disk
    harmonized RBFE XML (production charges; MC1 passes file-level)."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    return ats.build_inplace_res4_twocopy_system(
        leg="free", seed="s7", solvate=False,
        harmonize_common_charges=False)


def test_twocopy_outcome_attached(real_twocopy_unsolv):
    """The two-copy build attaches (or honestly surfaces the MC1 finding)."""
    assert real_twocopy_unsolv["outcome"] in (
        "twocopy_attached", "mc1_charge_discontinuity")


def test_twocopy_both_copies_resident(real_twocopy_unsolv):
    """C2: both endpoint copies are resident (n_atoms ~ 2x single endpoint)."""
    if real_twocopy_unsolv["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome — no attached box to inspect")
    fused = real_twocopy_unsolv["fused_build"]
    # copy-1 = 210 (MTR), copy-2 = 207 (WT) -> 417 unsolvated.
    assert fused["n_copy1"] == 210
    assert fused["n_atoms"] > 400
    assert fused["n_atoms"] == fused["system"].getNumParticles()


def test_twocopy_distinct_attach(real_twocopy_unsolv):
    """C2: the two copies have DISTINCT NE1 attach atoms (NOT a shared core)."""
    if real_twocopy_unsolv["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome")
    cmap = real_twocopy_unsolv["common_map"]
    assert cmap["copy1_attach"] != cmap["copy2_attach"]
    assert real_twocopy_unsolv["swap"]["distinct_attach"] is True


def test_twocopy_zero_inter_copy_exclusions(real_twocopy_unsolv):
    """C3 (LOAD-BEARING): ZERO inter-copy exclusions — clash avoided by the
    d-separation, NOT the forbidden overlay+exclusion design."""
    if real_twocopy_unsolv["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome")
    assert real_twocopy_unsolv["swap"]["inter_copy_exclusions_added"] == 0


def test_twocopy_spatial_separation(real_twocopy_unsolv, ats):
    """C2/C6: the two copies are spatially separated (NE1<->NE1 ~ d, NOT
    overlaid). With conformer registration the separation is exactly d."""
    if real_twocopy_unsolv["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome")
    sep = real_twocopy_unsolv["separation"]["ne1_ne1_sep_nm"]
    # Registered commons => offset is pure d (~4.0 nm). Always >> the PME cutoff.
    assert sep >= 1.0
    assert abs(sep - ats.ATS_TWOCOPY_DISPLACEMENT_NM) < 0.5


def test_twocopy_common_count_parity(real_twocopy_unsolv):
    """C4: the two copies' common-atom count parity (206 == 206)."""
    if real_twocopy_unsolv["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome")
    cmap = real_twocopy_unsolv["common_map"]
    assert cmap["n_common"] == 206
    assert len(cmap["copy1_common"]) == len(cmap["copy2_common"]) == 206
    assert len(cmap["copy1_var"]) == 4   # MTR {CM, HM1-3}
    assert len(cmap["copy2_var"]) == 1   # WT {HE1}


def test_twocopy_two_disulfides(real_twocopy_unsolv):
    """MC3: cyclic_ss preserved in BOTH copies (two disulfides)."""
    if real_twocopy_unsolv["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome")
    assert real_twocopy_unsolv["mc3_disulfide"]["n_disulfides"] == 2


def test_twocopy_mc1_continuity(real_twocopy_unsolv):
    """MC1: with the harmonized RBFE MTR XML the common core is continuous. The MTR
    fixture loads the on-disk harmonized RBFE XML, so the two-copy box attaches with
    a clean (per-atom continuous) common core."""
    assert real_twocopy_unsolv["outcome"] == "twocopy_attached"
    mc1 = real_twocopy_unsolv["mc1_param_continuity"]
    # Harmonized MTR core => per-atom continuous (the strict assert dict is kept,
    # carrying max_dq_e). If a future XML diverged per-atom, MC1 is now reporting-
    # only (the build would still attach with native charges + net sanity), but the
    # on-disk MTR RBFE XML is harmonized so this asserts the continuous path.
    assert mc1["passed"] is True
    assert mc1["max_dq_e"] <= 1e-4


def test_twocopy_endpoint_equivalence_and_nonsaturation(ats):
    """C6e/C8: u0 reproduces the full System potential (endpoint-equivalence) +
    the bulk copy is decoupled + |u1-u0| is NOT the single-shared-core saturated
    ~150 plateau (the collapse-signature escape)."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    build = ats.build_inplace_res4_twocopy_system(
        leg="free", seed="s7", solvate=False, harmonize_common_charges=True)
    if build["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome")
    eq = ats.check_twocopy_endpoint_equivalence(
        build, platform_name="Reference")
    assert eq["endpoint_equivalence"]["passed"] is True
    assert eq["bulk_copy_decoupled"]["passed"] is True
    assert eq["perturbation_regime"]["finite"] is True
    # The collapse signature was a |u1-u0| saturated near ~150 (Umax band); the
    # two-copy reference-frame perturbation must NOT be in that plateau.
    assert eq["perturbation_regime"]["saturated_plateau_flag"] is False
    assert eq["overall_pass"] is True


# ---------------------------------------------------------------------------
# Canonical-fix tests: ONE shared receptor + binder-only copy-2 (bound leg) +
# the whole-solute separation guard that catches receptor-receptor
# interpenetration the NE1<->NE1 gate is blind to.
# ---------------------------------------------------------------------------
def _toy_twocopy_separation_build(ats, copy2_shift_nm):
    """Minimal fake two-copy fused_build + cmap for the solute-separation guard.

    Two solute residues (2 heavy atoms each) + a couple of shared water atoms.
    copy-1 sits at the origin; copy-2 is translated by ``copy2_shift_nm`` on x.
    Returns ``(fused_build, cmap)`` exercising assert_twocopy_separation with no
    real build (deterministic). The NE1 attach atoms are deliberately placed far
    apart (4 nm on x) regardless of the shift so the NE1<->NE1 gate ALWAYS
    passes — the solute-solute gate is the one under test.
    """
    import openmm as mm
    from openmm import app
    import openmm.unit as unit

    top = app.Topology()
    chain = top.addChain(id="B")
    res1 = top.addResidue("TRP", chain, id="4")    # copy-1 solute
    res2 = top.addResidue("TRP", chain, id="104")  # copy-2 solute
    wchain = top.addChain(id="W")
    wres = top.addResidue("HOH", wchain, id="900")  # shared bath

    # copy-1 solute: NE1 (attach) + CA, at the origin.
    a_ne1_c1 = top.addAtom("NE1", app.element.nitrogen, res1)
    a_ca_c1 = top.addAtom("CA", app.element.carbon, res1)
    n_copy1 = top.getNumAtoms()  # boundary AFTER copy-1's atoms
    # copy-2 solute: NE1 (attach) + CA, shifted by copy2_shift_nm on x.
    a_ne1_c2 = top.addAtom("NE1", app.element.nitrogen, res2)
    a_ca_c2 = top.addAtom("CA", app.element.carbon, res2)
    # shared solvent: one O between the copies (must NOT trip the guard).
    a_o = top.addAtom("O", app.element.oxygen, wres)
    top.addAtom("H1", app.element.hydrogen, wres)  # H -> dropped by guard

    s = copy2_shift_nm
    positions = [
        mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer,        # NE1 copy-1
        mm.Vec3(0.1, 0.0, 0.0) * unit.nanometer,        # CA copy-1
        mm.Vec3(4.0 + s, 0.0, 0.0) * unit.nanometer,    # NE1 copy-2 (far on x)
        mm.Vec3(0.1 + s, 0.0, 0.0) * unit.nanometer,    # CA copy-2 (near c1 CA)
        mm.Vec3(2.0, 0.0, 0.0) * unit.nanometer,        # shared water O
        mm.Vec3(2.0, 0.1, 0.0) * unit.nanometer,        # shared water H
    ]

    class _M:
        pass
    m = _M()
    m.topology = top
    m.positions = positions
    fused = {"modeller": m}
    cmap = {
        "copy1_attach": a_ne1_c1.index,
        "copy2_attach": a_ne1_c2.index,
        "n_copy1": n_copy1,
    }
    return fused, cmap


def test_twocopy_separation_solute_gate_passes_when_displaced(ats):
    """The augmented C6 guard PASSES when the two solute copies are >= 1 nm apart
    (copy-2 displaced) and reports the new solute-solute key. Shared solvent
    between the copies is excluded."""
    fused, cmap = _toy_twocopy_separation_build(ats, copy2_shift_nm=4.0)
    res = ats.assert_twocopy_separation(fused, cmap, min_sep_nm=1.0)
    assert res["passed"] is True
    assert res["solute_solute_min_sep_nm"] >= 1.0
    # Hydrogen + solvent excluded: 2 heavy solute atoms per copy.
    assert res["n_copy1_solute_heavy"] == 2
    assert res["n_copy2_solute_heavy"] == 2


def test_twocopy_separation_solute_gate_catches_interpenetration(ats):
    """The augmented C6 guard RAISES on receptor/solute interpenetration that the
    NE1<->NE1 gate cannot see: the two NE1 attach atoms are 4 nm apart (NE1 gate
    passes), but the solute BODIES overlap (copy-2 not displaced) -> the
    solute-solute gate fails loud (the false-green the fix is designed to block)."""
    fused, cmap = _toy_twocopy_separation_build(ats, copy2_shift_nm=0.0)
    # Sanity: the NE1<->NE1 distance alone is still >= 1 nm (gate 1 would pass).
    import numpy as np
    import openmm.unit as unit
    pos = np.array([v.value_in_unit(unit.nanometer)
                    for v in fused["modeller"].positions])
    ne1_sep = float(np.linalg.norm(pos[cmap["copy1_attach"]]
                                   - pos[cmap["copy2_attach"]]))
    assert ne1_sep >= 1.0
    with pytest.raises(ValueError, match="SOLUTE separation FAIL"):
        ats.assert_twocopy_separation(fused, cmap, min_sep_nm=1.0)


@pytest.fixture(scope="module")
def real_twocopy_bound(ats):
    """The canonical two-copy BOUND box, UNSOLVATED, harmonized so MC1 passes:
    copy-1 = receptor + MTR binder, copy-2 = WT binder ONLY (receptor dropped)."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    return ats.build_inplace_res4_twocopy_system(
        leg="bound", seed="s7", solvate=False,
        harmonize_common_charges=True)


def test_twocopy_bound_copy2_is_binder_only(real_twocopy_bound, ats):
    """Canonical fix (a): on the BOUND leg copy-2 carries the binder ONLY (no
    receptor). copy-1 holds the shared receptor (large); copy-2 is the small
    binder-only copy whose atom count matches the free WT peptide (~207)."""
    if real_twocopy_bound["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome — no attached box to inspect")
    fused = real_twocopy_bound["fused_build"]
    n_copy1 = fused["n_copy1"]
    n_copy2 = fused["n_atoms"] - n_copy1
    # copy-1 = receptor + MTR binder -> thousands of atoms (receptor present).
    assert n_copy1 > 5000
    # copy-2 = binder ONLY -> ~207 atoms (receptor dropped); MUST NOT be a second
    # full bound complex (>5000). This is the interpenetration root cause removed.
    assert n_copy2 < 1000
    # copy-2 matches the free-leg WT binder-only count (canonical: single receptor
    # + duplicated binder). The free copy-2 (WT) is 207 atoms unsolvated.
    assert n_copy2 == 207


def test_twocopy_bound_solute_separation_reports_both_gates(real_twocopy_bound):
    """The BOUND two-copy box passes BOTH separation gates and the augmented
    result carries the whole-solute min distance (>= the PME floor)."""
    if real_twocopy_bound["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome")
    sep = real_twocopy_bound["separation"]
    assert sep["passed"] is True
    assert sep["ne1_ne1_sep_nm"] >= 1.0
    assert sep["solute_solute_min_sep_nm"] >= 1.0
    # copy-1 solute (receptor + MTR binder) is large; copy-2 solute (binder) small.
    assert sep["n_copy1_solute_heavy"] > sep["n_copy2_solute_heavy"]


def test_twocopy_free_build_unchanged_reports_solute_gate(real_twocopy_unsolv):
    """Regression (c): the FREE build still attaches and now also reports the new
    solute-solute separation key (the augmented guard runs on the free leg too,
    with no behavioural regression — both copies are already binder-only)."""
    if real_twocopy_unsolv["outcome"] != "twocopy_attached":
        pytest.skip("MC1 outcome")
    sep = real_twocopy_unsolv["separation"]
    assert sep["passed"] is True
    assert sep["ne1_ne1_sep_nm"] >= 1.0
    assert "solute_solute_min_sep_nm" in sep
    assert sep["solute_solute_min_sep_nm"] >= 1.0


# ===========================================================================
# MutationSpec generalization (C5 — res-4 MTR<->Trp default + res-3 Val<->Ile V3I)
# ===========================================================================
def test_mutation_specs_registered(ats):
    """The default + V3I specs are registered and resolvable by name."""
    assert "mtr_trp_res4" in ats.MUTATION_SPECS
    assert "v3i_val_ile_res3" in ats.MUTATION_SPECS
    assert ats.resolve_mutation_spec(None) is ats.MUTATION_MTR_TRP_RES4
    assert ats.resolve_mutation_spec("v3i_val_ile_res3") is ats.MUTATION_VAL_ILE_RES3
    # A MutationSpec instance round-trips unchanged.
    assert ats.resolve_mutation_spec(ats.MUTATION_VAL_ILE_RES3) \
        is ats.MUTATION_VAL_ILE_RES3


def test_mutation_spec_resolve_unknown_raises(ats):
    """An unknown spec name fails loud (no silent default)."""
    with pytest.raises(ValueError):
        ats.resolve_mutation_spec("not_a_spec")
    with pytest.raises(TypeError):
        ats.resolve_mutation_spec(42)


def test_default_spec_byte_identical_to_legacy_constants(ats):
    """The DEFAULT spec reproduces the legacy ALCH_* res-4 MTR<->Trp partition
    (byte-identical: stateB=MTR appears, stateA=WT-Trp disappears)."""
    ms = ats.MUTATION_MTR_TRP_RES4
    assert ms.resnum == ats.ALCH_RESNUM == 4
    assert ms.common_attach_atom == ats.ALCH_COMMON_ATOM == "NE1"
    # stateA disappearing == legacy wt_only ; stateB appearing == legacy mtr_only.
    assert list(ms.stateA_only_atoms) == list(ats.ALCH_WT_ONLY) == ["HE1"]
    assert list(ms.stateB_only_atoms) == list(ats.ALCH_MTR_ONLY) \
        == ["CM", "HM1", "HM2", "HM3"]
    assert ms.hybrid_xml is None         # MTR XML resolved at build time
    assert ms.bonded_heavy_appearing == "CM"
    assert ms.appearing_h_prefix == "HM"


def test_v3i_spec_amber_atom_partition(ats):
    """V3I res-3 Val<->Ile partition matches the amber14 VAL/ILE templates:
    common attach = CG1; disappearing = Val HG11; appearing = Ile CD1 + HD11-13.
    Canonical (no ncAA XML); charge-/parity-neutral (Val and Ile are both 0)."""
    ms = ats.MUTATION_VAL_ILE_RES3
    assert ms.resnum == 3
    assert ms.common_attach_atom == "CG1"
    assert ms.stateA_resname == "VAL"
    assert ms.stateB_resname == "ILE"
    assert list(ms.stateA_only_atoms) == ["HG11"]               # Val gamma-H
    assert list(ms.stateB_only_atoms) == ["CD1", "HD11", "HD12", "HD13"]  # Ile CH3
    assert ms.hybrid_xml is None                                # canonical amber14
    assert ms.bonded_heavy_appearing == "CD1"
    assert ms.appearing_h_prefix == "HD"


class _FakeAtom:
    def __init__(self, name, index):
        self.name = name
        self.index = index


class _FakeResidue:
    def __init__(self, rid, name, atoms):
        self.id = rid
        self.name = name
        self._atoms = atoms

    def atoms(self):
        return iter(self._atoms)


class _FakeChain:
    def __init__(self, cid, residues):
        self.id = cid
        self._res = residues

    def residues(self):
        return iter(self._res)


class _FakeTopology:
    def __init__(self, chains):
        self._chains = chains

    def chains(self):
        return iter(self._chains)


def _val_ile_topology():
    """A synthetic chain-B residue-3 with BOTH Val and Ile var atoms resident
    (a dual-topology box analog), plus the common core, so the partition can be
    classified by name. Indices are arbitrary but distinct."""
    names = [
        # common core (shared Val/Ile)
        "N", "H", "CA", "HA", "C", "O", "CB", "HB",
        "CG1", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23",
        # Val-only (disappearing)
        "HG11",
        # Ile-only (appearing)
        "CD1", "HD11", "HD12", "HD13",
    ]
    atoms = [_FakeAtom(n, i) for i, n in enumerate(names)]
    res = _FakeResidue("3", "VAL", atoms)
    return _FakeTopology([_FakeChain("B", [res])]), {a.name: a.index for a in atoms}


def test_identify_alchemical_atoms_v3i_partition(ats):
    """identify_alchemical_atoms with the V3I spec classifies CG1 as common,
    HG11 as disappearing (wt_only slot), CD1+HD11-13 as appearing (mtr_only slot)
    on a synthetic Val/Ile residue (no OpenMM build needed)."""
    top, idx = _val_ile_topology()
    part = ats.identify_alchemical_atoms(
        top, binder_chain="B", spec="v3i_val_ile_res3")
    assert part["common"] == [idx["CG1"]]
    assert part["wt_only"] == [idx["HG11"]]
    assert sorted(part["mtr_only"]) == sorted(
        [idx["CD1"], idx["HD11"], idx["HD12"], idx["HD13"]])


def test_identify_alchemical_atoms_default_spec_unchanged(ats):
    """identify_alchemical_atoms with spec=None on a synthetic res-4 MTR/Trp
    residue still classifies NE1/HE1/CM+HM1-3 (legacy partition, byte-identical)."""
    names = ["N", "H", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2",
             "CE3", "CZ2", "CZ3", "CH2",
             "HE1",                       # WT-only (disappears)
             "CM", "HM1", "HM2", "HM3"]   # MTR-only (appears)
    atoms = [_FakeAtom(n, i) for i, n in enumerate(names)]
    res = _FakeResidue("4", "MTR", atoms)
    top = _FakeTopology([_FakeChain("B", [res])])
    idx = {a.name: a.index for a in atoms}
    part = ats.identify_alchemical_atoms(top, binder_chain="B")  # spec=None
    assert part["common"] == [idx["NE1"]]
    assert part["wt_only"] == [idx["HE1"]]
    assert sorted(part["mtr_only"]) == sorted(
        [idx["CM"], idx["HM1"], idx["HM2"], idx["HM3"]])


def test_build_twocopy_signature_has_spec(ats):
    """The two-copy builder accepts the spec parameter (default None = legacy)."""
    import inspect
    sig = inspect.signature(ats.build_inplace_res4_twocopy_system)
    assert "spec" in sig.parameters
    assert sig.parameters["spec"].default is None


def test_twocopy_bound_branch_uses_binder_only_prep_in_source(ats):
    """Source-level invariant: the MTR bound branch of
    build_inplace_res4_twocopy_system builds copy-1 (MTR site) with the bound-
    complex prep and copy-2 (WT bulk) with the binder-only prep
    (prepare_free_peptide_from_final), NOT a second bound complex. Guards against a
    regression that re-introduces the duplicated full receptor (the
    interpenetration / NaN root cause).

    Whitespace-insensitive: the prep calls were line-wrapped when the builder
    gained the canonical (V3I) spec branch, so the assertions normalize internal
    whitespace before matching the call signatures.
    """
    import inspect
    import re
    src = inspect.getsource(ats.build_inplace_res4_twocopy_system)
    norm = re.sub(r"\s+", " ", src)
    # copy-1 (MTR) site copy on the bound leg uses the bound-complex prep.
    assert 'prepare_bound_complex_from_final( li["final"]["cp4"]' in norm
    # copy-2 (WT) on the bound leg uses the binder-only prep (receptor dropped).
    assert 'prepare_free_peptide_from_final( li["final"]["wt"]' in norm


# --- V3I (res-3 Val<->Ile) real two-copy build (need WT s7 final.pdb + openmm) ---
def _wt_s7_present():
    wt = os.path.join(_PROJ, "outputs", "2QKI_WT_calib_s7",
                      "mdresult", "2QKI_WT_final.pdb")
    return os.path.isfile(wt)


@pytest.fixture(scope="module")
def real_v3i_unsolv_harmonized(ats):
    """The V3I (res-3 Val<->Ile) canonical two-copy box, UNSOLVATED, with the
    diagnostic charge harmonization ON (isolates the mechanical mapping/swap/C6
    from the amber14 VAL!=ILE common-charge gap — see the MC1-finding test)."""
    if not _wt_s7_present():
        pytest.skip("2QKI WT s7 final.pdb not present")
    return ats.build_inplace_res4_twocopy_system(
        leg="free", seed="s7", solvate=False,
        harmonize_common_charges=True, spec="v3i_val_ile_res3")


def test_v3i_twocopy_mapping_and_separation(real_v3i_unsolv_harmonized):
    """V3I build smoke: common-core/swap-atom mapping + C6 separation + the
    canonical atom partition (attach=CG1, Ile CD1+HD11-13 appear, Val HG11
    disappears)."""
    b = real_v3i_unsolv_harmonized
    assert b["outcome"] == "twocopy_attached"
    fused = b["fused_build"]
    cmap = b["common_map"]
    name = {a.index: a.name for a in fused["modeller"].topology.atoms()}
    # Distinct CG1 attach atoms (NOT a shared core).
    assert cmap["copy1_attach"] != cmap["copy2_attach"]
    assert name[cmap["copy1_attach"]] == "CG1"
    assert name[cmap["copy2_attach"]] == "CG1"
    # Appearing = Ile CD1 + HD11-13; disappearing = Val HG11.
    assert sorted(name[i] for i in cmap["copy1_var"]) == [
        "CD1", "HD11", "HD12", "HD13"]
    assert sorted(name[i] for i in cmap["copy2_var"]) == ["HG11"]
    # C6 separation passes (clash floor + the registered offset is ~d).
    sep = b["separation"]
    assert sep["ne1_ne1_sep_nm"] > 1.0
    assert sep["solute_solute_min_sep_nm"] > 1.0
    # C3: zero inter-copy exclusions; MC2 (CD1 bonded) + MC3 (2 disulfides) pass.
    assert b["swap"]["inter_copy_exclusions_added"] == 0
    assert b["mc2_methyl_bonded"]["passed"] is True
    assert b["mc3_disulfide"]["n_disulfides"] == 2
    assert b["seed_assert"]["passed"] is True


def test_v3i_twocopy_mc1_finding_production(ats):
    """V3I PRODUCTION (non-harmonized) proceeds with NATIVE amber14 charges. MC1 is
    REPORTING-ONLY in the canonical two-copy box (
    ): the ATS swap is
    a coordinate-only transform, so per-atom common-charge divergence between the
    Val copy and the Ile copy is benign (u1-u0 includes the per-copy charge
    difference correctly). The build therefore ATTACHES with native charges; the
    amber14 VAL!=ILE common-charge divergence is surfaced as a non-blocking report.
    The retained C2 NET-charge sanity gate still holds (Val<->Ile is non-charge-
    changing => Sigma dq ~ 0)."""
    if not _wt_s7_present():
        pytest.skip("2QKI WT s7 final.pdb not present")
    b = ats.build_inplace_res4_twocopy_system(
        leg="free", seed="s7", solvate=False,
        harmonize_common_charges=False, spec="v3i_val_ile_res3")
    # NEW: native charges -> attach (no longer the mc1_charge_discontinuity fail).
    assert b["outcome"] == "twocopy_attached"
    # No in-memory harmonization happened (native amber VAL/ILE charges retained).
    assert b["common_charges_harmonized"] is False
    mc1 = b["mc1_param_continuity"]
    # Per-atom divergence is reported (amber14 VAL/ILE common atoms differ) but is
    # NON-BLOCKING (reporting-only).
    assert mc1["per_atom_continuous"] is False
    assert mc1["n_diverging"] > 0
    # C2 NET-charge sanity holds: the FULL alch-residue net charge (common + var)
    # is conserved across copies (Val neutral -> Ile neutral; both totals 0.0).
    assert mc1["net_sanity_ok"] is True
    assert abs(mc1["full_resmut_net_diff_e"]) <= mc1["net_dq_tol_e"]
    assert abs(mc1["full_resmut_net_copy1_e"]) < 1e-3   # Ile residue net ~ 0
    assert abs(mc1["full_resmut_net_copy2_e"]) < 1e-3   # Val residue net ~ 0
    # The COMMON subset alone carries a non-zero net (amber14 VAL!=ILE on shared
    # atoms ~ +0.089 e) — this is BENIGN (var atoms compensate), not the gate.
    assert abs(mc1["sum_dq_all_common_e"]) > 1e-3
    assert abs(mc1["sum_dq_res4_e"]) > 1e-3   # back-compat alias (res-3 common dq)


# --- C2 net-charge sanity (synthetic; no openmm endpoint build needed) ---------
def _toy_twocopy_charge_build(ats, names, q_copy1, q_copy2,
                              var1_q=(), var2_q=(), resnum=4):
    """Build a synthetic MERGED System (single NonbondedForce) + cmap + copy1_build
    that mimics the two-copy box's residue layout, so the C2 net-charge sanity
    logic in ``_summarize_twocopy_charge_divergence`` can be exercised without the
    heavy real endpoint build. Layout (matching the merge convention):
      copy-1 common atoms  -> [0, n)            (residue ``resnum``)
      copy-2 common atoms  -> [n, 2n)           (residue ``resnum``)
      copy-1 variable atoms-> [2n, 2n+v1)       (appearing-state atoms)
      copy-2 variable atoms-> [2n+v1, 2n+v1+v2) (disappearing-state atoms)
    The C2 gate is on the FULL mutated-residue net charge (common + var) per copy,
    so the variable-atom charges (``var1_q``/``var2_q``) complete the residue.
    """
    import openmm as mm
    from openmm import app
    import openmm.unit as unit
    n = len(names)
    system = mm.System()
    nb = mm.NonbondedForce()
    all_q = list(q_copy1) + list(q_copy2) + list(var1_q) + list(var2_q)
    for q in all_q:
        system.addParticle(1.0 * unit.dalton)
        nb.addParticle(q * unit.elementary_charge,
                       0.3 * unit.nanometer, 0.0 * unit.kilojoule_per_mole)
    system.addForce(nb)
    # copy1_build needs only modeller.topology.atoms() (index -> residue.id, name)
    # for the per-residue itemisation; copy-1 commons occupy indices [0, n).
    top = app.Topology()
    chain = top.addChain(id="B")
    res = top.addResidue("TRP", chain, id=str(resnum))
    for name in names:
        top.addAtom(name, app.element.carbon, res)

    class _M:
        pass
    m = _M()
    m.topology = top
    base = 2 * n
    cmap = {
        "copy1_common": list(range(n)),
        "copy2_common": list(range(n, 2 * n)),
        "copy1_var": list(range(base, base + len(var1_q))),
        "copy2_var": list(range(base + len(var1_q),
                                base + len(var1_q) + len(var2_q))),
        "n_common": n,
    }
    return system, cmap, {"modeller": m}


def test_c2_net_sanity_flags_nonzero_residue_net(ats):
    """C2 (retained HARD gate): an ARTIFICIAL build defect where the FULL mutated-
    residue net charge differs between the copies (the residue total is NOT
    conserved) must be flagged net_sanity_ok=False. This is the build-bug signal
    the two-copy MC1 raise keys on — per-atom (and common-subset) divergence is
    benign, a non-conserved RESIDUE TOTAL is not."""
    names = ["CA", "CB", "CG", "CD"]
    # copy-1 residue total (common + var) = 0.00; copy-2 = +0.50 (defect).
    q1 = [-0.10, 0.05, 0.00, 0.05]          # common sum = 0.00
    q2 = [-0.10, 0.05, 0.00, 0.05]          # common sum = 0.00
    var1 = (0.0,)                            # copy-1 var total = 0.00
    var2 = (0.5,)                            # copy-2 var total = +0.50 -> defect
    system, cmap, c1b = _toy_twocopy_charge_build(
        ats, names, q1, q2, var1_q=var1, var2_q=var2)
    rep = ats._summarize_twocopy_charge_divergence(system, cmap, c1b, resnum=4)
    assert rep["net_sanity_ok"] is False
    assert abs(rep["full_resmut_net_diff_e"] + 0.5) < 1e-6   # copy1-copy2 = -0.5
    assert abs(rep["full_resmut_net_copy1_e"]) < 1e-6
    assert abs(rep["full_resmut_net_copy2_e"] - 0.5) < 1e-6


def test_c2_net_sanity_passes_val_ile_analog(ats):
    """C2: the V3I analog — the COMMON subset carries a non-zero net dq (amber14
    VAL!=ILE on shared atoms) but the variable atoms carry the EXACT complementary
    charge, so the FULL residue net is conserved -> net_sanity_ok=True even though
    the common subset and per-atom charges diverge. The canonical two-copy build
    proceeds with native charges in exactly this case (REPORTING-ONLY)."""
    names = ["CA", "CB", "CG", "CD"]
    # Common subset diverges by +0.30 net (copy1 - copy2); the disappearing var on
    # copy-2 carries +0.30 to compensate so BOTH residue totals are 0.00.
    q1 = [-0.10, 0.05, 0.10, -0.05]         # common sum copy1 = 0.00
    q2 = [-0.20, -0.05, 0.05, -0.10]        # common sum copy2 = -0.30
    var1 = (0.0,)                           # copy-1 (appearing) var total = 0.00
    var2 = (0.30,)                          # copy-2 (disappearing) var = +0.30
    system, cmap, c1b = _toy_twocopy_charge_build(
        ats, names, q1, q2, var1_q=var1, var2_q=var2)
    rep = ats._summarize_twocopy_charge_divergence(system, cmap, c1b, resnum=4)
    assert rep["net_sanity_ok"] is True            # FULL residue net conserved
    assert rep["per_atom_continuous"] is False     # common atoms still diverge
    assert rep["n_diverging"] >= 1
    assert abs(rep["full_resmut_net_diff_e"]) <= rep["net_dq_tol_e"]
    assert abs(rep["sum_dq_all_common_e"] - 0.30) < 1e-6   # common-subset net != 0


def test_c2_net_sanity_passes_when_identical(ats):
    """Fully continuous common core (no var atoms) -> net_sanity_ok AND
    per_atom_continuous."""
    names = ["CA", "CB", "CG", "CD"]
    q = [-0.10, 0.05, 0.10, -0.05]
    system, cmap, c1b = _toy_twocopy_charge_build(ats, names, q, q)
    rep = ats._summarize_twocopy_charge_divergence(system, cmap, c1b, resnum=4)
    assert rep["net_sanity_ok"] is True
    assert rep["per_atom_continuous"] is True
    assert rep["passed"] is True
    assert rep["n_diverging"] == 0


def test_v3i_prepare_mutated_binder_atoms(ats):
    """prepare_mutated_binder_from_final mutates res-3 VAL->ILE with the correct
    amber14 atom set (HG11 removed, CD1+HD11-13 added) and stamps the binder
    chain id so the downstream build resolves it."""
    if not _wt_s7_present():
        pytest.skip("2QKI WT s7 final.pdb not present")
    import tempfile
    import openmm.app as app
    li = ats.resolve_leg_inputs("s7")
    tmp = tempfile.mkdtemp(prefix="v3i_test_")
    out = os.path.join(tmp, "ile.pdb")
    ats.prepare_mutated_binder_from_final(
        li["final"]["wt"], out, 3, "VAL", "ILE", "B", add_hydrogens=True)
    pdb = app.PDBFile(out)
    # Binder chain id preserved as 'B' (PDBFile.writeFile re-letters; we re-stamp).
    assert [c.id for c in pdb.topology.chains()] == ["B"]
    res3 = None
    for res in pdb.topology.residues():
        if str(res.id) == "3":
            res3 = res
            break
    assert res3 is not None and res3.name == "ILE"
    names = {a.name for a in res3.atoms()}
    assert {"CD1", "HD11", "HD12", "HD13"} <= names   # Ile gamma-CH3 present
    assert "HG11" not in names                          # Val gamma-H removed
    assert "CG1" in names                               # common attach retained


def test_twocopy_no_inter_copy_exclusions_in_source(ats):
    """C3 source-level invariant: the attach function reports ZERO inter-copy
    exclusions by construction (the wiring adds none — mirrors the upstream
    add_common_var_atoms_to_atmforce which adds none)."""
    import inspect
    src = inspect.getsource(ats.attach_twocopy_swap_atmforce)
    # No addException call in the swap-wiring function.
    assert "addException" not in src
    assert '"inter_copy_exclusions_added": 0' in src


def test_legacy_single_core_path_unchanged(real_fused_harmonized):
    """R-7/C9: the LEGACY single-shared-core path still produces its fused box
    (the new two-copy path did not break it). densify_pilot depends on this."""
    assert real_fused_harmonized["outcome"] == "fused_attached"
    assert real_fused_harmonized["swap_mode"] == "genuine"
    assert real_fused_harmonized["swap"]["swap_mode"] == "genuine"


# --- Two-copy smoke-script plumbing ---
def _load_twocopy_smoke():
    spec = importlib.util.spec_from_file_location(
        "trackb_inplace_res4_twocopy_smoke",
        os.path.join(_PROJ, "scripts", "trackb_inplace_res4_twocopy_smoke.py"))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate scripts/trackb_inplace_res4_twocopy_smoke.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_twocopy_smoke_module_imports():
    if not _have_openmm():
        pytest.skip("openmm not importable")
    smoke = _load_twocopy_smoke()
    assert hasattr(smoke, "run_tier1_twocopy_smoke")
    assert hasattr(smoke, "main")


def test_twocopy_smoke_tier1_free_unsolvated(ats):
    """The two-copy Tier-1 smoke runs end-to-end (Reference, unsolvated) and
    reports a structured outcome (pass / endpoint-mismatch / MC1)."""
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    smoke = _load_twocopy_smoke()
    res = smoke.run_tier1_twocopy_smoke(
        seed="s7", solvate=False, harmonize_common_charges=True,
        lam=0.5, platform_name="Reference", leg="free")
    assert res["swap_mode"] == "twocopy"
    assert res["outcome"] in (
        "tier1_twocopy_pass", "tier1_twocopy_endpoint_mismatch",
        "mc1_charge_discontinuity")
    if res["outcome"] == "tier1_twocopy_pass":
        assert res["swap"]["inter_copy_exclusions_added"] == 0
        assert res["swap"]["distinct_attach"] is True
        assert abs(res["energies_kcal"]["u1_minus_u0"]) <= 1.0e3
        assert res["endpoint_equivalence"]["overall_pass"] is True


# ---------------------------------------------------------------------------
# Bound-complex PBC re-imaging (the cycle-0 / NaN root cause). An MD final.pdb
# written with PBC unwrapping can leave the binder a whole box vector away from
# the receptor; prepare_bound_complex_from_final must minimum-image the binder
# back into receptor contact (pose-preserving) and fail-fast on a genuinely
# broken endpoint. These are pure-coordinate tests (no openmm/GPU); a synthetic
# 2-chain PDB is written so they run in any env. A separate real-data test
# (s7 cp4) is gated on the endpoint being present.
# ---------------------------------------------------------------------------
_REIMAGE_BOX = 30.0


def _write_synthetic_bound_pdb(path, binder_offset=(0.0, 0.0, 0.0),
                               box=_REIMAGE_BOX):
    """Write a minimal 2-chain bound complex: a small receptor cluster (chain A)
    near the origin and a binder cluster (chain B) ~3 A away, optionally
    translated by ``binder_offset`` (used to fake a PBC-image displacement).
    Includes a CRYST1 record. Heavy atoms only (element column set)."""
    rec = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.0, 1.5, 0.0), (0.0, 0.0, 1.5)]
    # Binder cluster ~3 A from the nearest receptor atom (a valid bound contact).
    bnd = [(3.0 + dx, dy, dz) for (dx, dy, dz) in
           ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))]
    ox, oy, oz = binder_offset
    bnd = [(x + ox, y + oy, z + oz) for (x, y, z) in bnd]
    lines = ["CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n"
             % (box, box, box)]
    serial = 1
    for (x, y, z) in rec:
        lines.append("ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                     "           C  \n" % (serial, serial, x, y, z))
        serial += 1
    lines.append("TER\n")
    for (x, y, z) in bnd:
        lines.append("ATOM  %5d  CA  ALA B%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                     "           C  \n" % (serial, serial, x, y, z))
        serial += 1
    lines.append("TER\nEND\n")
    with open(path, "w") as fh:
        fh.writelines(lines)
    return path


def _ab_min_heavy_dist(pdb):
    import math
    a = []
    b = []
    for line in open(pdb):
        if line[:6] in ("ATOM  ", "HETATM"):
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if line[21] == "A":
                a.append(xyz)
            elif line[21] == "B":
                b.append(xyz)
    best = float("inf")
    for (ax, ay, az) in a:
        for (bx, by, bz) in b:
            d = math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
            if d < best:
                best = d
    return best


def test_bound_reimage_preserves_contacting_complex(ats):
    """A normal complex (receptor-binder in contact) is unchanged by re-imaging:
    no box shift, contact distance preserved."""
    import os
    import tempfile
    td = tempfile.mkdtemp(prefix="ats_reimage_normal_")
    src = _write_synthetic_bound_pdb(os.path.join(td, "in.pdb"))
    before = _ab_min_heavy_dist(src)
    out = os.path.join(td, "out.pdb")
    ats.prepare_bound_complex_from_final(src, out, "B", "A")
    after = _ab_min_heavy_dist(out)
    # Contact preserved to floating-point precision (no spurious shift).
    assert abs(after - before) < 1e-3
    assert after < ats._BOUND_CONTACT_MAX_A


def test_bound_reimage_restores_pbc_split_binder(ats):
    """A binder displaced by one whole box vector (PBC-unwrapped) is brought back
    into receptor contact by the minimum-image shift — the cp4 cycle-0 fix."""
    import os
    import tempfile
    td = tempfile.mkdtemp(prefix="ats_reimage_split_")
    # Push the binder one box vector along +z (an integer image displacement).
    src = _write_synthetic_bound_pdb(
        os.path.join(td, "in.pdb"), binder_offset=(0.0, 0.0, _REIMAGE_BOX))
    assert _ab_min_heavy_dist(src) > ats._BOUND_CONTACT_MAX_A   # broken on input
    out = os.path.join(td, "out.pdb")
    ats.prepare_bound_complex_from_final(src, out, "B", "A")
    after = _ab_min_heavy_dist(out)
    # Re-imaged back to the original ~3 A contact.
    assert after < ats._BOUND_CONTACT_MAX_A
    assert after < 5.0


def test_bound_reimage_guard_raises_on_broken_endpoint(ats):
    """A genuinely separated endpoint (binder NOT one box image away — it sits
    mid-box, beyond contact, with no integer shift that restores contact) trips
    the fail-fast guard instead of emitting a non-bound complex."""
    import os
    import tempfile
    td = tempfile.mkdtemp(prefix="ats_reimage_broken_")
    # Offset by HALF a box: min-image rounds to 0 shift, so it stays separated.
    src = _write_synthetic_bound_pdb(
        os.path.join(td, "in.pdb"),
        binder_offset=(0.0, 0.0, _REIMAGE_BOX / 2.0))
    out = os.path.join(td, "out.pdb")
    with pytest.raises(ValueError, match="min heavy-atom distance"):
        ats.prepare_bound_complex_from_final(src, out, "B", "A")


def test_bound_reimage_missing_cryst1_raises(ats):
    """A final.pdb without a CRYST1 record cannot be re-imaged (PBC box unknown)
    and is a hard error — never a silent pass."""
    import os
    import tempfile
    td = tempfile.mkdtemp(prefix="ats_reimage_nobox_")
    src = _write_synthetic_bound_pdb(os.path.join(td, "in.pdb"))
    # Strip the CRYST1 line.
    kept = [ln for ln in open(src) if ln[:6] != "CRYST1"]
    nobox = os.path.join(td, "nobox.pdb")
    with open(nobox, "w") as fh:
        fh.writelines(kept)
    out = os.path.join(td, "out.pdb")
    with pytest.raises(ValueError, match="no CRYST1"):
        ats.prepare_bound_complex_from_final(nobox, out, "B", "A")


def test_bound_reimage_real_cp4_endpoint(ats):
    """On the real s7 cp4 final.pdb (PBC-unwrapped: binder ~44 A from receptor),
    re-imaging restores receptor-binder contact (< guard threshold) — the actual
    failing endpoint that crashed bound asyncre at cycle 0."""
    import os
    import tempfile
    if not _endpoints_present():
        pytest.skip("2QKI endpoint final.pdb not present (s7)")
    li = ats.resolve_leg_inputs("s7")
    td = tempfile.mkdtemp(prefix="ats_reimage_cp4_")
    out = os.path.join(td, "cp4_bound.pdb")
    ats.prepare_bound_complex_from_final(li["final"]["cp4"], out, "B", "A")
    after = _ab_min_heavy_dist(out)
    assert after < ats._BOUND_CONTACT_MAX_A
    # Real bound contact lands at ~2.6 A after the one-image (z) shift.
    assert after < 5.0


# ===========================================================================
# Task #100: direction-aware displacement auto-search + PBC minimum-image C6 gate.
# All synthetic-coordinate / pure-numeric — no real build, no GPU.
# ===========================================================================
def test_min_image_distance_raw_matches_euclidean(ats):
    """Without a box, _min_image_min_distance_nm is the raw Euclidean min dist."""
    import numpy as np
    c1 = np.array([[0.0, 0.0, 0.0]])
    c2 = np.array([[3.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    d = ats._min_image_min_distance_nm(c1, c2, None)
    assert abs(d - 3.0) < 1e-9


def test_min_image_distance_detects_wrap(ats):
    """A copy-2 atom placed far in raw coords but whose nearest periodic IMAGE is
    close to copy-1 is measured at the small image distance (the Q5/C3 wrap)."""
    import numpy as np
    box = np.array([10.0, 10.0, 10.0])
    c1 = np.array([[0.0, 0.0, 0.0]])
    # raw 9.5 nm away on x, but the -10 image lands it at -0.5 nm => 0.5 nm.
    c2 = np.array([[9.5, 0.0, 0.0]])
    raw = ats._min_image_min_distance_nm(c1, c2, None)
    img = ats._min_image_min_distance_nm(c1, c2, box)
    assert abs(raw - 9.5) < 1e-9
    assert abs(img - 0.5) < 1e-9


def test_box_lengths_from_vectors_parses_quantity_and_array(ats):
    """_box_lengths_nm_from_vectors reads OpenMM Quantity Vec3 rows and bare nm
    arrays; returns None for a missing / malformed box."""
    import openmm as mm
    import openmm.unit as unit
    bv = [mm.Vec3(4.0, 0.0, 0.0) * unit.nanometer,
          mm.Vec3(0.0, 5.0, 0.0) * unit.nanometer,
          mm.Vec3(0.0, 0.0, 6.0) * unit.nanometer]
    L = ats._box_lengths_nm_from_vectors(bv)
    assert L is not None
    assert abs(L[0] - 4.0) < 1e-9 and abs(L[1] - 5.0) < 1e-9 and abs(L[2] - 6.0) < 1e-9
    # Bare nested list (nm).
    L2 = ats._box_lengths_nm_from_vectors([[2.0, 0, 0], [0, 3.0, 0], [0, 0, 4.0]])
    assert abs(L2[2] - 4.0) < 1e-9
    # None / degenerate box -> None.
    assert ats._box_lengths_nm_from_vectors(None) is None
    assert ats._box_lengths_nm_from_vectors(
        [[0.0, 0, 0], [0, 0.0, 0], [0, 0, 0.0]]) is None


def test_candidate_directions_base_first_and_unit(ats):
    """_candidate_directions returns the base direction FIRST then a cone ring; all
    are unit vectors; the cone members are tilted ~cone_deg off the base."""
    import numpy as np
    base = (1.0, 0.0, 0.0)
    cands = ats._candidate_directions(base, n_candidates=9, cone_deg=60.0)
    assert len(cands) == 9
    # base first (reproduces the legacy direction when it already clears).
    assert np.allclose(cands[0], base, atol=1e-9)
    for c in cands:
        assert abs(float(np.linalg.norm(c)) - 1.0) < 1e-9
    # Ring members make ~60 deg with the base (cos 60 = 0.5).
    for c in cands[1:]:
        cosang = float(np.dot(c, base))
        assert abs(cosang - 0.5) < 1e-6


def test_separation_image_gate_catches_wrap(ats):
    """assert_twocopy_separation's NEW periodic-image gate fails a build whose
    copy-2 wraps back near copy-1 even though the raw distance passed."""
    import numpy as np
    import openmm as mm
    import openmm.unit as unit
    # Reuse the toy builder but push copy-2 far enough that its image wraps.
    fused, cmap = _toy_twocopy_separation_build(ats, copy2_shift_nm=4.0)
    # raw distances are fine (>= 1 nm). Now supply a SMALL box so copy-2's image
    # wraps close to copy-1: copy-2 solute sits near x ~ 4.1 nm; a 5 nm box images
    # it to ~ -0.9 nm => image distance < 1 nm.
    box = [mm.Vec3(5.0, 0.0, 0.0) * unit.nanometer,
           mm.Vec3(0.0, 5.0, 0.0) * unit.nanometer,
           mm.Vec3(0.0, 0.0, 5.0) * unit.nanometer]
    # Raw-only still passes (box=None).
    res_raw = ats.assert_twocopy_separation(fused, cmap, min_sep_nm=1.0)
    assert res_raw["passed"] is True
    assert res_raw["image_solute_min_sep_nm"] is None
    # With the small box the image gate fires.
    with pytest.raises(ValueError, match="PERIODIC-IMAGE separation FAIL"):
        ats.assert_twocopy_separation(fused, cmap, min_sep_nm=1.0, box_vectors=box)


def test_separation_acceptance_line_rejects_below_threshold(ats):
    """The elevated decoupling-sufficient acceptance line rejects a build that
    clears the 1.0 nm clash floor but not the 1.5 nm acceptance line; the legacy
    default (accept_sep_nm=None) leaves the 1.0 nm-only behaviour byte-identical."""
    # Build copies separated by ~1.2 nm solute-solute (clears 1.0, not 1.5).
    fused, cmap = _toy_twocopy_separation_build(ats, copy2_shift_nm=1.1)
    # Legacy default: passes at the 1.0 nm floor.
    res = ats.assert_twocopy_separation(fused, cmap, min_sep_nm=1.0)
    assert res["passed"] is True
    assert res["solute_solute_min_sep_nm"] >= 1.0
    assert res["solute_solute_min_sep_nm"] < 1.5
    # Elevated acceptance line: fails (decoupling-insufficient).
    with pytest.raises(ValueError, match="ACCEPTANCE separation FAIL"):
        ats.assert_twocopy_separation(
            fused, cmap, min_sep_nm=1.0, accept_sep_nm=1.5)


# --- Synthetic per-copy builds for the auto-search (no real build / GPU). ---
def _toy_copy_build(ats, heavy_xyz_nm, ne1_xyz_nm=None, resid="4"):
    """A minimal per-copy build dict (modeller w/ topology+positions) carrying a
    handful of solute heavy atoms (+ NE1) for the displacement auto-search. The
    auto-search reads only solute heavy positions, so a few atoms suffice."""
    import numpy as np
    import openmm as mm
    from openmm import app
    import openmm.unit as unit
    top = app.Topology()
    chain = top.addChain(id="B")
    res = top.addResidue("TRP", chain, id=resid)
    positions = []
    if ne1_xyz_nm is None:
        ne1_xyz_nm = heavy_xyz_nm[0]
    top.addAtom("NE1", app.element.nitrogen, res)
    positions.append(mm.Vec3(*ne1_xyz_nm) * unit.nanometer)
    for i, xyz in enumerate(heavy_xyz_nm):
        top.addAtom("C%d" % i, app.element.carbon, res)
        positions.append(mm.Vec3(*xyz) * unit.nanometer)

    class _M:
        pass
    m = _M()
    m.topology = top
    m.positions = positions
    return {"modeller": m}


def test_autosearch_accepts_base_direction_when_clear(ats):
    """When the base outward direction already clears the acceptance line at the
    first magnitude, the auto-search accepts it (candidate_index 0) without
    escalating the magnitude — reproducing the legacy choice."""
    import numpy as np
    # copy-1 receptor-like cloud centred near origin; copy-2 small binder near it.
    c1 = _toy_copy_build(
        ats, [(0.0, 0.0, 0.0), (0.3, 0.0, 0.0), (0.0, 0.3, 0.0), (0.0, 0.0, 0.3)])
    c2 = _toy_copy_build(ats, [(0.1, 0.0, 0.0), (0.2, 0.1, 0.0)])
    out = ats.auto_search_twocopy_displacement(
        c1, c2, accept_sep_nm=1.5,
        magnitudes_nm=(4.0, 5.5, 7.0), n_candidates=9)
    assert out["accepted"] is True
    assert out["magnitude_nm"] == 4.0          # first magnitude sufficed
    assert out["achieved_raw_min_nm"] >= 1.5
    assert out["achieved_image_min_nm"] >= 1.5
    # The realised vector has the chosen magnitude.
    v = np.array(out["displacement_vector_nm"])
    assert abs(float(np.linalg.norm(v)) - 4.0) < 1e-6


def test_autosearch_maximises_distance_picks_best_candidate(ats):
    """The auto-search reports the candidate with the largest (image-aware) min
    distance and a full trail; every record is finite and scored."""
    c1 = _toy_copy_build(
        ats, [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (0.0, 0.5, 0.0)])
    c2 = _toy_copy_build(ats, [(0.1, 0.1, 0.0), (0.2, 0.0, 0.1)])
    out = ats.auto_search_twocopy_displacement(
        c1, c2, accept_sep_nm=1.5, magnitudes_nm=(4.0,), n_candidates=12)
    # The chosen score is the max over all candidates at the accepted magnitude.
    mag_recs = [r for r in out["trail"] if r["magnitude_nm"] == 4.0]
    best_score = max(r["score_nm"] for r in mag_recs)
    assert out["accepted"] is True
    assert abs(min(out["achieved_raw_min_nm"], out["achieved_image_min_nm"])
               - best_score) < 1e-6
    for r in out["trail"]:
        assert r["raw_min_nm"] > 0.0 and r["image_min_nm"] > 0.0


def test_autosearch_escalates_magnitude_when_small_d_insufficient(ats):
    """When the smallest magnitude leaves the raw distance below the acceptance
    line, the search escalates to a larger magnitude that clears it — the automated
    d=6->8->10 recovery (with adequate padding so no periodic wrap intervenes)."""
    c1 = _toy_copy_build(
        ats, [(0.0, 0.0, 0.0), (0.3, 0.0, 0.0), (0.0, 0.3, 0.0)])
    c2 = _toy_copy_build(ats, [(0.1, 0.0, 0.0), (0.2, 0.1, 0.0)])
    # 1.0 nm is below the 1.5 nm acceptance; 4.0 nm clears it. padding default-ish.
    out = ats.auto_search_twocopy_displacement(
        c1, c2, accept_sep_nm=1.5, magnitudes_nm=(1.0, 4.0),
        n_candidates=9, padding_nm=1.2)
    assert out["accepted"] is True
    assert out["magnitude_nm"] == 4.0          # escalated past the 1.0 nm try
    assert out["n_magnitudes_tried"] == 2
    assert out["achieved_raw_min_nm"] >= 1.5
    assert out["achieved_image_min_nm"] >= 1.5


def test_autosearch_raises_when_exhausted(ats):
    """When no magnitude clears the acceptance line (every tried d leaves the raw
    distance below it), the search raises a deterministic final failure naming the
    best achieved distance (caller -> Path escalate, no silent pass)."""
    c1 = _toy_copy_build(
        ats, [(0.0, 0.0, 0.0), (0.3, 0.0, 0.0), (0.0, 0.3, 0.0)])
    c2 = _toy_copy_build(ats, [(0.1, 0.0, 0.0)])
    with pytest.raises(ValueError, match="no candidate direction/magnitude cleared"):
        # All magnitudes (~1.0 nm) sit below the 1.5 nm acceptance line; adequate
        # padding so the failure is genuine raw insufficiency, not a wrap artifact.
        ats.auto_search_twocopy_displacement(
            c1, c2, accept_sep_nm=1.5, magnitudes_nm=(1.0, 1.0, 1.0),
            n_candidates=6, padding_nm=1.2)


def test_autosearch_empty_copy_raises(ats):
    """A copy with no solute heavy atoms is a hard error (cannot score)."""
    c1 = _toy_copy_build(ats, [(0.0, 0.0, 0.0)])
    # copy-2 with an explicitly empty heavy list (only the NE1 N is heavy; remove
    # it by passing a build whose topology has only a hydrogen).
    import openmm as mm
    from openmm import app
    import openmm.unit as unit
    top = app.Topology()
    ch = top.addChain(id="B")
    r = top.addResidue("HOH", ch, id="900")   # solvent -> excluded
    top.addAtom("O", app.element.oxygen, r)

    class _M:
        pass
    m = _M()
    m.topology = top
    m.positions = [mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer]
    c2 = {"modeller": m}
    with pytest.raises(ValueError, match="no .*solute heavy atoms"):
        ats.auto_search_twocopy_displacement(c1, c2, accept_sep_nm=1.5)
