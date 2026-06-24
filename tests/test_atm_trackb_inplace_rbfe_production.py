# -*- coding: utf-8 -*-
"""Tests for the in-place residue-4 RBFE per-direction PRODUCTION launcher.

Covers the launcher's pure / logic surface WITHOUT a real OpenMM build or GPU:
  - the C8 SIGN-critical bound-leg decouple-direction gate;
  - the matched-seed availability gate;
  - the C7 Welch-Satterthwaite quadrature SE + z_SE sign rule;
  - the COMBINED + single-direction cntl writers (UWHAM-consumable / SSOT);
  - the C11 pre-registration of the 4 outcomes;
  - the soft-core usc helper (the .out pertE must be the CAPPED perturbation,
    not the raw u1-u0 that would overflow the WHAM solve).

Uses the importlib spec-loader pattern (no sys.modules pollution; the launcher +
bridge are openmm-only and the qmmm env has no atom_openmm).
"""

import importlib.util
import json
import math
import os
import sys

import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)
_UTILS = os.path.join(_PROJ, "utils")
_SCRIPTS = os.path.join(_PROJ, "scripts")


def _load(name, relpath):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_PROJ, relpath))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate %s" % relpath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _have_openmm_units():
    try:
        import openmm.unit  # noqa: F401
        import openmm.vec3  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.fixture(scope="module")
def prod():
    return _load("trackb_inplace_rbfe_production",
                 "scripts/trackb_inplace_rbfe_production.py")


@pytest.fixture(scope="module")
def rbfe():
    if _UTILS not in sys.path:
        sys.path.insert(0, _UTILS)
    return _load("atm_trackB_inplace_rbfe", "utils/atm_trackB_inplace_rbfe.py")


# --------------------------- C8 decouple-direction gate --------------------
def test_c8_free_leg_none_is_ok(prod):
    r = prod.gate_decouple_direction("free", None, raise_on_fail=True)
    assert r["passed"] is True
    assert "free_leg" in r["reason"]


def test_c8_bound_leg_none_fails_loud(prod):
    with pytest.raises(RuntimeError):
        prod.gate_decouple_direction("bound", None, raise_on_fail=True)
    # non-raising form reports the violation.
    r = prod.gate_decouple_direction("bound", None, raise_on_fail=False)
    assert r["passed"] is False and "C8 VIOLATION" in r["reason"]


def test_c8_bound_leg_unit_vector_passes(prod):
    # An outward unit vector (the kind compute_decouple_direction returns).
    v = (0.0, 0.0, 1.0)
    r = prod.gate_decouple_direction("bound", v, raise_on_fail=True)
    assert r["passed"] is True
    assert abs(r["magnitude"] - 1.0) < 1e-6


def test_c8_bound_leg_non_unit_fails(prod):
    with pytest.raises(RuntimeError):
        prod.gate_decouple_direction("bound", (0.0, 0.0, 2.0), raise_on_fail=True)


def test_c8_bound_leg_non_finite_fails(prod):
    with pytest.raises(RuntimeError):
        prod.gate_decouple_direction("bound", (float("nan"), 0.0, 1.0),
                                     raise_on_fail=True)


def test_c8_bound_leg_wrong_length_fails(prod):
    with pytest.raises(RuntimeError):
        prod.gate_decouple_direction("bound", (0.0, 1.0), raise_on_fail=True)


# ---------------- C5 bound two-copy receptor-contact pre-flight gate ----------
def _load_ats_mod():
    if _UTILS not in sys.path:
        sys.path.insert(0, _UTILS)
    return _load("atm_trackB_setup", "utils/atm_trackB_setup.py")


class _CGAtom:
    def __init__(self, name, index, chain_id, resname, symbol="C"):
        self.name = name
        self.index = index
        self.element = type("_El", (), {"symbol": symbol})()
        chain = type("_Ch", (), {"id": chain_id})()
        self.residue = type("_Res", (), {"name": resname, "chain": chain})()


class _CGTopology:
    def __init__(self, atoms):
        self._atoms = atoms

    def atoms(self):
        return iter(self._atoms)


class _CGModeller:
    def __init__(self, topology, positions):
        self.topology = topology
        self.positions = positions


class _CGSystem:
    def __init__(self, box_nm):
        self._box = box_nm

    def getDefaultPeriodicBoxVectors(self):
        return self._box


def _synthetic_bound_fused(separation_nm):
    """Build a synthetic two-copy BOUND ``fused`` dict: copy-1 = a receptor heavy
    (chain A) + a binder heavy (chain B), copy-2 = a displaced binder heavy placed
    ``separation_nm`` from the receptor along +x. A big orthorhombic box so no
    wrap. Positions are OpenMM nm Quantities (the gate value_in_units them)."""
    import openmm.unit as unit
    import openmm.vec3 as _v3
    nm = unit.nanometer
    # copy-1: receptor heavy at origin (chain A), binder heavy at +0.3 (chain B).
    a_rec = _CGAtom("CA", 0, "A", "ALA", "C")
    a_bnd = _CGAtom("CB", 1, "B", "LEU", "C")
    # copy-2: displaced binder heavy at +separation along x.
    c2_bnd = _CGAtom("CB", 2, "B", "LEU", "C")
    top = _CGTopology([a_rec, a_bnd, c2_bnd])
    pos = [
        _v3.Vec3(0.0, 0.0, 0.0) * nm,
        _v3.Vec3(0.3, 0.0, 0.0) * nm,
        _v3.Vec3(float(separation_nm), 0.0, 0.0) * nm,
    ]
    box = [_v3.Vec3(50.0, 0.0, 0.0) * nm,
           _v3.Vec3(0.0, 50.0, 0.0) * nm,
           _v3.Vec3(0.0, 0.0, 50.0) * nm]
    return {
        "modeller": _CGModeller(top, pos),
        "system": _CGSystem(box),
        "n_copy1": 2,          # copy-1 atoms are indices [0, 2)
    }


@pytest.mark.skipif(not _have_openmm_units(),
                    reason="openmm not importable in this env")
def test_c5_contact_gate_not_applicable_for_free(prod):
    """The contact gate is bound two-copy ONLY: free / single-core pass through
    (no measurement, no HALT)."""
    ats = _load_ats_mod()
    r = prod.gate_receptor_contact("free", "twocopy", None, ats,
                                   raise_on_fail=True)
    assert r["passed"] is True and "not applicable" in r["reason"]
    r2 = prod.gate_receptor_contact("bound", "single_core", None, ats,
                                    raise_on_fail=True)
    assert r2["passed"] is True and "not applicable" in r2["reason"]


@pytest.mark.skipif(not _have_openmm_units(),
                    reason="openmm not importable in this env")
def test_c5_contact_gate_accept_when_decoupled(prod):
    """ACCEPT: the displaced binder is 5 nm from the receptor -> 0 contacts AND
    min-image >= 1 nm -> the gate passes (genuine decouple)."""
    ats = _load_ats_mod()
    fused = _synthetic_bound_fused(separation_nm=5.0)
    r = prod.gate_receptor_contact("bound", "twocopy", fused, ats,
                                   endpoint="wt", seed="s7", raise_on_fail=True)
    assert r["passed"] is True
    m = r["measurement"]
    assert m["n_contacts"] == 0
    assert m["min_image_dist_nm"] >= prod.RECEPTOR_DECOUPLE_MIN_NM
    assert m["decoupled"] is True


@pytest.mark.skipif(not _have_openmm_units(),
                    reason="openmm not importable in this env")
def test_c5_contact_gate_halt_when_in_contact(prod):
    """HALT (fail-loud): the displaced binder is 0.2 nm from the receptor -> a
    contact < 0.45 nm -> the gate RAISES (no silent skip)."""
    ats = _load_ats_mod()
    fused = _synthetic_bound_fused(separation_nm=0.2)
    with pytest.raises(RuntimeError) as exc:
        prod.gate_receptor_contact("bound", "twocopy", fused, ats,
                                   endpoint="wt", seed="s163",
                                   raise_on_fail=True)
    assert "C5 VIOLATION" in str(exc.value)
    # The non-raising form reports the violation honestly.
    r = prod.gate_receptor_contact("bound", "twocopy", fused, ats,
                                   endpoint="wt", seed="s163",
                                   raise_on_fail=False)
    assert r["passed"] is False
    assert r["measurement"]["decoupled"] is False
    assert r["measurement"]["n_contacts"] >= 1


@pytest.mark.skipif(not _have_openmm_units(),
                    reason="openmm not importable in this env")
def test_c5_contact_gate_halt_when_under_displaced(prod):
    """HALT: the displaced binder clears the 0.45 nm contact cutoff but sits at
    0.7 nm < the 1.0 nm decouple floor (still pocket-coupled) -> RAISES."""
    ats = _load_ats_mod()
    fused = _synthetic_bound_fused(separation_nm=0.7)
    with pytest.raises(RuntimeError):
        prod.gate_receptor_contact("bound", "twocopy", fused, ats,
                                   endpoint="wt", seed="s199",
                                   raise_on_fail=True)


@pytest.mark.skipif(not _have_openmm_units(),
                    reason="openmm not importable in this env")
def test_c5_contact_gate_missing_build_fails_loud(prod):
    """A bound two-copy gate with no live build dict (e.g. box-reuse path) fails
    loud rather than silently passing."""
    ats = _load_ats_mod()
    with pytest.raises(RuntimeError):
        prod.gate_receptor_contact("bound", "twocopy", None, ats,
                                   raise_on_fail=True)


# --------------------------- seed-availability gate ------------------------
def test_seed_gate_default_seeds_present(prod):
    # The default matched seeds (s7, s101, s127) must exist for both endpoints
    # in this repo — that is the basis of the paired ddG. If this fails the
    # cohort drifted and the default must be updated (a real, not silent, gap).
    r = prod.gate_seed_availability(list(prod.ENDPOINTS), list(prod.DEFAULT_SEEDS))
    assert r["passed"] is True, r["missing"]


def test_seed_gate_missing_seed_fails(prod):
    r = prod.gate_seed_availability(["cp4", "wt"], ["s_does_not_exist_999"])
    assert r["passed"] is False
    assert len(r["missing"]) == 2  # both endpoints missing this seed


# --------------------------- C7 quadrature / z_SE --------------------------
def test_quadrature_se_and_z(prod):
    q = prod._welch_satterthwaite_quadrature(10.0, 1.0, 7.0, 1.0)
    assert abs(q["ddg"] - 3.0) < 1e-12
    assert abs(q["se_quadrature"] - math.sqrt(2.0)) < 1e-9
    assert abs(q["z_se"] - 3.0 / math.sqrt(2.0)) < 1e-9


def test_quadrature_zero_se_returns_none_z(prod):
    q = prod._welch_satterthwaite_quadrature(1.0, 0.0, 1.0, 0.0)
    assert q["z_se"] is None


def test_z_se_threshold_is_three(prod):
    # The pre-registration sign rule: z_SE >= 3 AND n >= 3 (R-18 honest).
    assert prod.Z_SE_SIGN_THRESHOLD == 3.0
    assert prod.MIN_REPLICATES_FOR_SIGN == 3


# --------------------------- cntl writers ----------------------------------
def test_combined_cntl_has_both_directions(prod, rbfe, tmp_path):
    sched = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                   single_direction=None)
    cntl = str(tmp_path / "trackb_asyncre.cntl")
    prod._write_combined_cntl(cntl, sched, "free", 250, 1.0)
    text = open(cntl).read()
    assert "MODE = 'INPLACE_RBFE'" in text
    # The combined cntl is the SSOT: DIRECTION must carry BOTH a +1 and a -1
    # block (UWHAM requires both legs); LIGAND_ATOMS / DISPLACEMENT must NOT
    # appear (that is the ABFE whole-binder path, not the in-place swap).
    dir_line = [ln for ln in text.splitlines() if ln.startswith("DIRECTION")][0]
    assert "1" in dir_line and "-1" in dir_line
    assert "LIGAND_ATOMS" not in text and "DISPLACEMENT" not in text


def test_direction_cntl_single_direction_only(prod, rbfe, tmp_path):
    sched = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                   single_direction="forward")
    cntl = str(tmp_path / "trackb_dplus_asyncre.cntl")
    prod._write_direction_cntl(cntl, sched, "dplus", 250, 1.0)
    text = open(cntl).read()
    assert "BASENAME = 'trackb_dplus'" in text
    dir_line = [ln for ln in text.splitlines() if ln.startswith("DIRECTION")][0]
    # forward standalone -> only +1, no -1.
    assert "-1" not in dir_line


# --------------------------- C11 pre-registration --------------------------
def test_pre_registration_has_four_outcomes(prod, tmp_path):
    path = prod.write_pre_registration(str(tmp_path), {"leg": "free"})
    import json
    payload = json.load(open(path))
    ids = [o["id"] for o in payload["outcomes"]]
    assert len(ids) == 4
    assert "O1_favorable_sign_resolved" in ids
    assert "O3_sign_undetermined" in ids
    assert "O4_non_convergent" in ids
    # The sign-claim policy must state z_SE >= 3 AND n >= 3 (no silent assume).
    assert "z_SE >= 3" in payload["sign_claim_policy"]


# --------------------------- soft-core usc helper --------------------------
def test_softcore_usc_caps_large_perturbation(rbfe):
    # A raw u1-u0 of ~50,000 kcal/mol (the single-shared-core clash) must be
    # capped to < umax. The .out pertE column carries usc, NOT the raw value —
    # otherwise UWHAM's _bias_fcn (no cap of its own) overflows the WHAM solve.
    kcal = 4.184
    umax_kj = 200.0 * kcal
    ubcore_kj = 100.0 * kcal
    acore = 0.0625
    big_pert_kj = 50000.0 * kcal
    base, usc = rbfe._atm_softcore_components_kj(
        big_pert_kj, -17000.0 * kcal, 1.0, umax_kj, ubcore_kj, acore)
    usc_kcal = usc / kcal
    assert usc_kcal < 200.0, usc_kcal
    assert usc_kcal > 100.0  # above ubcore (it was a large perturbation)


def test_softcore_usc_small_perturbation_passthrough(rbfe):
    kcal = 4.184
    umax_kj = 200.0 * kcal
    ubcore_kj = 100.0 * kcal
    acore = 0.0625
    # A small (< ubcore) perturbation passes through uncapped.
    small_kj = 5.0 * kcal
    base, usc = rbfe._atm_softcore_components_kj(
        small_kj, 0.0, 1.0, umax_kj, ubcore_kj, acore)
    assert abs(usc / kcal - 5.0) < 1e-9


def test_softcore_usc_backward_sign(rbfe):
    # Backward (Direction<0) flips the sign of the perturbation + bases on u1.
    kcal = 4.184
    umax_kj = 200.0 * kcal
    ubcore_kj = 100.0 * kcal
    base_f, usc_f = rbfe._atm_softcore_components_kj(
        10.0 * kcal, 0.0, 1.0, umax_kj, ubcore_kj, 0.0625)
    base_b, usc_b = rbfe._atm_softcore_components_kj(
        10.0 * kcal, 0.0, -1.0, umax_kj, ubcore_kj, 0.0625)
    # forward base = u0 = 0; backward base = u1 = u0 + pert.
    assert abs(base_f - 0.0) < 1e-9
    assert abs(base_b - 10.0 * kcal) < 1e-6
    # usc sign flips.
    assert usc_f > 0 and usc_b < 0


def test_state_energy_uses_softcore_components(rbfe):
    # _atm_state_energy_kj must equal base + bias_fcn(usc) so the .out potE +
    # capped pertE reconstruct in UWHAM (e0 = potE - bias_fcn(usc) = base).
    kcal = 4.184
    umax_kj = 200.0 * kcal
    ubcore_kj = 100.0 * kcal
    acore = 0.0625
    pert_kj = 8.0 * kcal
    u0_kj = -100.0 * kcal
    lambda1, lambda2 = 0.3, 0.4
    alpha_per_kj = 0.10 / kcal
    uh_kj = 30.0 * kcal
    w0_kj = 0.0
    base, usc = rbfe._atm_softcore_components_kj(
        pert_kj, u0_kj, 1.0, umax_kj, ubcore_kj, acore)
    U = rbfe._atm_state_energy_kj(
        pert_kj, u0_kj, lambda1, lambda2, alpha_per_kj, uh_kj, w0_kj,
        1.0, umax_kj, ubcore_kj, acore)
    # Reconstruct bias_fcn(usc) in kJ and check U == base + bias.
    bias = (((lambda2 - lambda1) / alpha_per_kj)
            * math.log(1.0 + math.exp(-alpha_per_kj * (usc - uh_kj)))
            + lambda2 * usc + w0_kj)
    assert abs(U - (base + bias)) < 1e-6


# --------------------------- parse helpers ---------------------------------
def test_parse_directions_default_both(prod):
    assert prod._parse_directions(None) == ["dplus", "dminus"]
    assert prod._parse_directions("dplus") == ["dplus"]
    with pytest.raises(ValueError):
        prod._parse_directions("dbogus")


def test_parse_seeds_default_and_custom(prod):
    assert prod._parse_seeds(None) == list(prod.DEFAULT_SEEDS)
    assert prod._parse_seeds("sA,sB") == ["sA", "sB"]
    with pytest.raises(ValueError):
        prod._parse_seeds(",")


def test_max_concurrent_units_packs_to_vram(prod):
    # 14 GiB free, 2 GiB/unit, 1.5 GiB headroom -> floor(12.5/2) = 6 units.
    assert prod._max_concurrent_units(14.0, 2.0) == 6
    # bound: 30 GiB free, 6.5 GiB/unit -> floor(28.5/6.5) = 4.
    assert prod._max_concurrent_units(30.0, 6.5) == 4
    # never below 1.
    assert prod._max_concurrent_units(1.0, 6.5) == 1


# --------------------------- box reuse (--reuse-serialized) ----------------
def test_reuse_serialized_flag_in_parser(prod):
    p = prod.build_arg_parser()
    a = p.parse_args(["--reuse-serialized", "--no-archive-existing",
                      "--leg", "bound", "--directions", "dminus"])
    assert a.reuse_serialized is True
    assert a.no_archive_existing is True
    # default off (legacy behaviour unchanged).
    b = p.parse_args(["--leg", "bound"])
    assert b.reuse_serialized is False


def test_reuse_serialized_requires_no_archive_existing(prod):
    # The safety gate: reuse WITHOUT --no-archive-existing exits 2 before any
    # run (otherwise the existing box_A + producing-direction outputs would be
    # archived/moved away by run_one_replicate's R-7 archive step).
    rc = prod.main(["--reuse-serialized", "--leg", "bound",
                    "--directions", "dminus", "--dry-run"])
    assert rc == 2


def test_reuse_serialized_with_no_archive_dry_run_ok(prod, tmp_path):
    # With --no-archive-existing the gate passes; --dry-run runs nothing.
    rc = prod.main(["--reuse-serialized", "--no-archive-existing",
                    "--leg", "bound", "--endpoints", "cp4",
                    "--directions", "dminus", "--dry-run",
                    "--out-root", str(tmp_path)])
    assert rc == 0


# ---------------- loaded-box decouple recompute (box reuse C8) -------------
def test_loaded_build_adapter_exposes_modeller_and_system(rbfe):
    # The adapter must answer build["modeller"] (subscript) and build.get(
    # "system") the way ats.compute_decouple_direction reads a build dict.
    sentinel_sys = object()
    sentinel_topo = object()
    sentinel_pos = object()
    ad = rbfe._LoadedBuildAdapter({
        "system": sentinel_sys, "topology": sentinel_topo,
        "positions": sentinel_pos})
    assert ad["modeller"].topology is sentinel_topo
    assert ad["modeller"].positions is sentinel_pos
    assert ad.get("system") is sentinel_sys
    assert ad["system"] is sentinel_sys
    # an unrelated key is absent.
    assert ad.get("missing") is None
    with pytest.raises(KeyError):
        _ = ad["missing"]


def test_recompute_decouple_dir_none_when_ne1_absent(rbfe):
    # If NE1 cannot be resolved (no chain B / no res 4), the recompute returns
    # None (the caller's C8 gate accepts None ONLY for the free leg).
    import openmm as mm
    from openmm.app import Topology, element

    topo = Topology()
    chain = topo.addChain(id="A")            # not the binder chain "B"
    res = topo.addResidue("ALA", chain, id="1")
    topo.addAtom("CA", element.carbon, res)
    system = mm.System()
    system.addParticle(12.0)
    import openmm.unit as unit
    positions = [mm.Vec3(0.0, 0.0, 0.0)] * 1 * unit.nanometer
    loaded = {"system": system, "topology": topo, "positions": positions}
    assert rbfe.recompute_decouple_direction_from_loaded(
        loaded, binder_chain="B") is None


# ---------------- two-copy construction wiring (opt-in) --------------------
def test_twocopy_flag_in_parser(prod):
    p = prod.build_arg_parser()
    a = p.parse_args(["--twocopy", "--displacement-nm", "4.0"])
    assert a.twocopy is True
    assert a.displacement_nm == 4.0
    # default OFF (single_core legacy, byte-identical).
    b = p.parse_args(["--leg", "free"])
    assert b.twocopy is False
    assert b.displacement_nm is None


def test_single_direction_schedule_construction_dispatch(prod, rbfe):
    # single_core -> the validated single-shared-core ladder (2*n states).
    sc = prod._build_single_direction_schedule(
        rbfe, construction="single_core", direction="forward",
        n_windows_half=6, softcore_band=2, n_apex_bridge=0, apex_band=0.5)
    assert sc.get("schedule_kind", "rbfe_ladder") != "ats_standard"
    assert sc["n_states"] == 6
    # twocopy -> the canonical ATS standard schedule (2*n - 1 states, Uh=110).
    tc = prod._build_single_direction_schedule(
        rbfe, construction="twocopy", direction="forward",
        n_windows_half=6, softcore_band=2, n_apex_bridge=0, apex_band=0.5)
    assert tc["schedule_kind"] == "ats_standard"
    assert tc["n_states"] == 11
    assert tc["u0"][0] == 110.0
    assert set(tc["directions"]) == {1}


def test_combined_schedule_twocopy_has_both_directions(prod, rbfe):
    cs = prod._build_combined_schedule(
        rbfe, construction="twocopy",
        n_windows_half=6, softcore_band=2, n_apex_bridge=0, apex_band=0.5)
    # 11 forward (+1) + 11 backward (-1) = 22 states; both blocks present so the
    # merge can derive (total=22, fwd=11) from the DIRECTION column.
    assert cs["n_states"] == 22
    assert cs["directions"].count(1) == 11
    assert cs["directions"].count(-1) == 11
    assert cs["schedule_kind"] == "ats_standard"


def test_combined_cntl_twocopy_merge_count_derivation(prod, rbfe, tmp_path):
    # The combined cntl for the two-copy ATS schedule must let the merge derive
    # (total, fwd) from its DIRECTION column (22/11).
    cs = prod._build_combined_schedule(
        rbfe, construction="twocopy",
        n_windows_half=6, softcore_band=2, n_apex_bridge=0, apex_band=0.5)
    cntl = str(tmp_path / "trackb_asyncre.cntl")
    prod._write_combined_cntl(cntl, cs, "free", 250, 1.0)
    driver = _load("trackb_per_direction_production",
                   "scripts/trackb_per_direction_production.py")
    total, fwd = driver._derive_state_counts_from_cntl(cntl)
    assert total == 22
    assert fwd == 11


def test_twocopy_dry_run_surfaces_ats_schedule(prod, tmp_path, capsys):
    rc = prod.main(["--twocopy", "--leg", "free", "--endpoints", "cp4",
                    "--directions", "dplus", "--seeds", "s7", "--dry-run",
                    "--out-root", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    plan = json.loads(out[out.index("{"):])["plan"]
    assert plan["construction"] == "twocopy"
    assert plan["schedule_kind"] == "ats_standard"
    assert plan["n_lambda_per_leg"] == 11
    assert plan["u0_kcal"] == 110.0


def test_default_dry_run_is_single_core(prod, tmp_path, capsys):
    rc = prod.main(["--leg", "free", "--endpoints", "cp4",
                    "--directions", "dplus", "--seeds", "s7", "--dry-run",
                    "--out-root", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    plan = json.loads(out[out.index("{"):])["plan"]
    assert plan["construction"] == "single_core"
    # the single-core plan does NOT carry the two-copy ATS extras.
    assert "schedule_kind" not in plan


# ---------------- auto-search displacement wiring (task #100/#114) ----------
def test_auto_search_flags_in_parser(prod):
    p = prod.build_arg_parser()
    a = p.parse_args(["--twocopy", "--auto-search-displacement",
                      "--accept-sep-nm", "2.0"])
    assert a.auto_search_displacement is True
    assert a.accept_sep_nm == 2.0
    # default OFF + the default acceptance line mirrors the engine constant.
    b = p.parse_args(["--leg", "free"])
    assert b.auto_search_displacement is False
    assert b.accept_sep_nm == prod._ATS_ACCEPT_SEP_NM_DEFAULT


def test_accept_sep_default_locksteps_with_engine_constant(prod):
    # The launcher's argparse default MUST equal the engine SSOT constant (the
    # comment promises lockstep; a drift is a wiring error).
    if _UTILS not in sys.path:
        sys.path.insert(0, _UTILS)
    ats = _load("atm_trackB_setup", "utils/atm_trackB_setup.py")
    assert prod._ATS_ACCEPT_SEP_NM_DEFAULT == ats.ATS_TWOCOPY_ACCEPT_SEP_NM


def test_auto_search_requires_twocopy(prod, tmp_path):
    # --auto-search-displacement on the single-core path is a wiring error: fail
    # loud (rc 2), never silently ignore (the single-core box has no copy-2 bulk
    # displacement to search).
    rc = prod.main(["--auto-search-displacement", "--leg", "free",
                    "--endpoints", "cp4", "--seeds", "s7",
                    "--directions", "dplus", "--out-root", str(tmp_path)])
    assert rc == 2


def test_auto_search_dry_run_registers_policy(prod, tmp_path, capsys):
    # With the search ON the C11 pre-registration records the SINGLE displacement
    # policy (anti-HARKing + Keeper-auditable) and the plan shows the mode.
    rc = prod.main(["--twocopy", "--auto-search-displacement",
                    "--accept-sep-nm", "1.5", "--leg", "bound",
                    "--endpoints", "cp4", "--directions", "dplus",
                    "--seeds", "s7", "--dry-run", "--out-root", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    payload = json.loads(out[out.index("{"):])
    ds = payload["config"]["displacement_search"]
    assert ds["auto_search_displacement"] is True
    assert ds["accept_sep_nm"] == 1.5
    assert "identically" in ds["applies_to"]
    assert payload["plan"]["displacement_mode"] == "auto_search"
    assert payload["plan"]["accept_sep_nm"] == 1.5
    # The on-disk pre_registration.json carries the same policy.
    prereg = json.load(open(os.path.join(str(tmp_path), "pre_registration.json")))
    assert (prereg["config"]["displacement_search"]["auto_search_displacement"]
            is True)


def test_default_off_registers_fixed_direction(prod, tmp_path, capsys):
    # Default OFF => the prereg policy fields are None (byte-identical legacy
    # fixed-direction semantics) and the plan shows fixed_direction.
    rc = prod.main(["--twocopy", "--leg", "bound", "--endpoints", "cp4",
                    "--directions", "dplus", "--seeds", "s7", "--dry-run",
                    "--out-root", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    payload = json.loads(out[out.index("{"):])
    ds = payload["config"]["displacement_search"]
    assert ds["auto_search_displacement"] is False
    assert ds["accept_sep_nm"] is None
    assert payload["plan"]["displacement_mode"] == "fixed_direction"
    assert "accept_sep_nm" not in payload["plan"]


# ---------------- git provenance stamping (run_manifest) -------------------
def test_git_provenance_keys_and_types(prod):
    # The helper must always return both keys with the right types, regardless
    # of whether git resolves (this repo IS a git repo, so commit resolves).
    prov = prod._git_provenance()
    assert set(prov.keys()) == {"git_commit", "git_dirty"}
    assert isinstance(prov["git_commit"], str)
    assert isinstance(prov["git_dirty"], bool)
    # In a real git checkout the commit is a 40-hex sha (or the 'unknown'
    # sentinel if git is somehow unavailable) — never empty.
    assert prov["git_commit"]
    if prov["git_commit"] != "unknown":
        assert len(prov["git_commit"]) == 40
        assert all(c in "0123456789abcdef" for c in prov["git_commit"])


def test_git_provenance_swallows_failures(prod, monkeypatch):
    # A subprocess explosion must degrade to ('unknown', False) and NEVER raise
    # — the manifest is metadata, it can never break a launch.
    import subprocess as _sp

    def _boom(*a, **k):
        raise OSError("git not found")

    monkeypatch.setattr(_sp, "check_output", _boom)
    prov = prod._git_provenance()
    assert prov == {"git_commit": "unknown", "git_dirty": False}


# ---------------------------------------------------------------------------
# FIX-A re-seeding (W4A bound-leg ladder mixing) — launcher flag plumbing.
# DEFAULT OFF must be byte/behaviour-identical; the flags must thread all the way
# to the InplaceRbfeLadder kwargs + be recorded to the run_manifest config.
# ---------------------------------------------------------------------------
def test_reseed_cli_flags_exist_with_off_defaults(prod):
    p = prod.build_arg_parser()
    ns = p.parse_args(["--leg", "bound"])
    assert ns.reseed_perm_seed is None
    assert ns.reseed_endpoint is False
    assert ns.reseed_endpoint_band_lambda2_max == 0.25
    assert ns.reseed_endpoint_equil_steps == 2000


def test_reseed_cli_flags_parse_on(prod):
    p = prod.build_arg_parser()
    ns = p.parse_args([
        "--leg", "bound", "--twocopy", "--mutation", "w4a_trp_ala_res4",
        "--reseed-perm-seed", "20260618", "--reseed-endpoint",
        "--reseed-endpoint-band-lambda2-max", "0.25",
        "--reseed-endpoint-equil-steps", "2000"])
    assert ns.reseed_perm_seed == 20260618
    assert ns.reseed_endpoint is True
    assert ns.reseed_endpoint_band_lambda2_max == 0.25
    assert ns.reseed_endpoint_equil_steps == 2000


def test_reseed_threads_through_run_signatures(prod):
    """The re-seed opt-ins must thread through the whole call chain (run_leg ->
    run_one_replicate -> run_one_direction) with the SAME default-OFF defaults, so
    a default launch is byte-identical and an opt-in launch reaches the ladder."""
    import inspect
    for fn in (prod.run_leg, prod.run_one_replicate, prod.run_one_direction):
        sig = inspect.signature(fn)
        assert sig.parameters["reseed_perm_seed"].default is None
        assert sig.parameters["reseed_endpoint"].default is False
        assert sig.parameters["reseed_endpoint_band_lambda2_max"].default == 0.25
        assert sig.parameters["reseed_endpoint_equil_steps"].default == 2000


def test_reseed_off_config_records_defaults(prod, tmp_path):
    """A default (no-reseed) dry-run records reseed OFF in the pre-registration
    config (reproducibility provenance, default-OFF visible)."""
    out_root = str(tmp_path / "off")
    rc = prod.main(["--leg", "bound", "--twocopy", "--mutation",
                    "w4a_trp_ala_res4", "--dry-run", "--out-root", out_root])
    assert rc == 0
    prereg = json.load(open(os.path.join(out_root, "pre_registration.json")))
    cfg = prereg["config"]
    assert cfg["reseed_perm_seed"] is None
    assert cfg["reseed_endpoint"] is False


def test_reseed_on_config_records_seed(prod, tmp_path):
    """A --reseed-perm-seed / --reseed-endpoint dry-run records the LOGGED seed +
    endpoint flag in the pre-registration config (reproducibility)."""
    out_root = str(tmp_path / "on")
    rc = prod.main(["--leg", "bound", "--twocopy", "--mutation",
                    "w4a_trp_ala_res4", "--reseed-perm-seed", "20260618",
                    "--reseed-endpoint", "--lambda1-rampdown",
                    "0.05,0.1,0.2,0.3,0.4,0.5", "--dry-run",
                    "--out-root", out_root])
    assert rc == 0
    prereg = json.load(open(os.path.join(out_root, "pre_registration.json")))
    cfg = prereg["config"]
    assert cfg["reseed_perm_seed"] == 20260618
    assert cfg["reseed_endpoint"] is True


def test_reseed_pool_cmd_propagates_flags(prod):
    """The pool cmd builder must propagate the re-seed flags to worker subprocesses
    when ON, and OMIT them when OFF (so default workers are byte-identical)."""
    import inspect
    src = inspect.getsource(prod.run_pool_local)
    assert "--reseed-perm-seed" in src
    assert "--reseed-endpoint" in src
    assert "reseed_perm_seed" in src and "reseed_endpoint" in src


# ---------------------------------------------------------------------------
# FIX-C deep-λ2 decouple-tail densification: --lambda2-rampdown (engine kwarg
# lambda2_rampup). Leg-UP λ2 axis, independent + composable with the leg-DOWN
# --lambda1-rampdown. Two-copy ONLY (single_core has no soft-core leg-up).
# ---------------------------------------------------------------------------
def test_lambda2_rampdown_flag_in_parser(prod):
    p = prod.build_arg_parser()
    a = p.parse_args(["--twocopy",
                      "--lambda2-rampdown", "0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5"])
    assert a.lambda2_rampdown == "0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5"
    # default OFF (canonical uniform leg-up, byte-identical).
    b = p.parse_args(["--leg", "free"])
    assert b.lambda2_rampdown is None


def test_parse_lambda2_rampdown_floats(prod):
    assert prod._parse_lambda2_rampdown(None) is None
    assert prod._parse_lambda2_rampdown("") is None
    assert prod._parse_lambda2_rampdown(
        "0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5") == [
            0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    with pytest.raises(ValueError):
        prod._parse_lambda2_rampdown("0.0,abc,0.5")


def test_lambda2_rampdown_threads_through_run_signatures(prod):
    """lambda2_rampup must thread through the whole call chain (run_leg ->
    run_one_replicate -> run_one_direction) with the SAME default-OFF default."""
    import inspect
    for fn in (prod.run_leg, prod.run_one_replicate, prod.run_one_direction):
        sig = inspect.signature(fn)
        assert sig.parameters["lambda2_rampup"].default is None
    # the two schedule builders also accept it.
    for fn in (prod._build_single_direction_schedule,
               prod._build_combined_schedule):
        sig = inspect.signature(fn)
        assert sig.parameters["lambda2_rampup"].default is None


def test_lambda2_rampdown_single_direction_densifies_legup(prod, rbfe):
    """twocopy single-direction schedule densifies the leg-up deep-λ2 tail."""
    tc = prod._build_single_direction_schedule(
        rbfe, construction="twocopy", direction="forward",
        n_windows_half=6, softcore_band=2, n_apex_bridge=0, apex_band=0.5,
        lambda2_rampup=[0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
    assert tc["schedule_kind"] == "ats_standard"
    # 8 leg-up + 5 uniform leg-down = 13 states.
    assert tc["n_states"] == 13
    assert 0.05 in tc["lambdas_2"]
    assert 0.15 in tc["lambdas_2"]
    assert tc["u0"][0] == 110.0


def test_lambda2_rampdown_combined_has_both_directions(prod, rbfe):
    cs = prod._build_combined_schedule(
        rbfe, construction="twocopy",
        n_windows_half=6, softcore_band=2, n_apex_bridge=0, apex_band=0.5,
        lambda2_rampup=[0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
    # 13 forward (+1) + 13 backward (-1) = 26 states.
    assert cs["n_states"] == 26
    assert cs["directions"].count(1) == 13
    assert cs["directions"].count(-1) == 13


def test_lambda2_rampdown_composable_with_lambda1_rampdown(prod, rbfe):
    """Both axes together: independent densification, single shared apex (14)."""
    tc = prod._build_single_direction_schedule(
        rbfe, construction="twocopy", direction="forward",
        n_windows_half=6, softcore_band=2, n_apex_bridge=0, apex_band=0.5,
        lambda2_rampup=[0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5],
        lambda1_rampdown=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    assert tc["n_states"] == 14
    apex_hits = [k for k in range(tc["n_states"])
                 if tc["lambdas_1"][k] == 0.0 and tc["lambdas_2"][k] == 0.5]
    assert apex_hits == [7]  # exactly one shared apex


def test_lambda2_rampdown_two_copy_only_rejected_for_single_core(prod, tmp_path):
    """--lambda2-rampdown without --twocopy fails loud (single_core has no
    soft-core leg-up); exit 2 (mirror of --mutation gate)."""
    rc = prod.main(["--leg", "free", "--endpoints", "cp4", "--seeds", "s7",
                    "--lambda2-rampdown", "0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5",
                    "--dry-run", "--out-root", str(tmp_path)])
    assert rc == 2


def test_lambda2_rampdown_dry_run_surfaces_densified_schedule(prod, tmp_path,
                                                              capsys):
    rc = prod.main(["--twocopy", "--leg", "bound", "--mutation",
                    "w4a_trp_ala_res4", "--endpoints", "cp4",
                    "--directions", "dplus", "--seeds", "s7",
                    "--lambda2-rampdown", "0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5",
                    "--dry-run", "--out-root", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    plan = json.loads(out[out.index("{"):])["plan"]
    assert plan["construction"] == "twocopy"
    assert plan["n_lambda_per_leg"] == 13
    assert plan["lambda2_rampup"] == [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    # the deep-λ2 bridges are present in the leg-up tail.
    assert 0.05 in plan["lambdas_2"]
    assert 0.15 in plan["lambdas_2"]


def test_lambda2_rampdown_off_config_records_none(prod, tmp_path):
    """A default (no-lambda2) dry-run records lambda2_rampup None in the
    pre-registration config (provenance, default-OFF visible)."""
    out_root = str(tmp_path / "off")
    rc = prod.main(["--leg", "bound", "--twocopy", "--mutation",
                    "w4a_trp_ala_res4", "--dry-run", "--out-root", out_root])
    assert rc == 0
    prereg = json.load(open(os.path.join(out_root, "pre_registration.json")))
    assert prereg["config"]["lambda2_rampup"] is None


def test_lambda2_rampdown_on_config_records_knots(prod, tmp_path):
    out_root = str(tmp_path / "on")
    rc = prod.main(["--leg", "bound", "--twocopy", "--mutation",
                    "w4a_trp_ala_res4", "--lambda2-rampdown",
                    "0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5", "--dry-run",
                    "--out-root", out_root])
    assert rc == 0
    prereg = json.load(open(os.path.join(out_root, "pre_registration.json")))
    assert prereg["config"]["lambda2_rampup"] == [
        0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]


def test_lambda2_rampdown_pool_cmd_propagates_flag(prod):
    """The pool cmd builder propagates --lambda2-rampdown to worker subprocesses
    (re-emitted from the dispatcher's parsed knot list)."""
    import inspect
    src = inspect.getsource(prod.run_pool_local)
    assert "--lambda2-rampdown" in src
    assert "lambda2_rampup" in src
