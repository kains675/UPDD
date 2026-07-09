# -*- coding: utf-8 -*-
"""Tests for the barnase Ile96->Ala folding-thermocycle ADDITIVE system-prep layer.

Covers the four additive engine changes (frozen FE core / soft-core / UWHAM / H18
ladder ALL untouched — those remain covered by the existing suites):
  1. MutationSpec.chained_group_certified + the acyclic_connected_group shape;
  2. make_ile_ala_mutation_spec SSOT factory (correct atom partition);
  3. resolve_fold_leg_inputs single-scaffold override (schema + fail-loud);
  4. assert_twocopy_disulfides cysteine-conditional guard (skip when 0 CYS,
     unchanged >=2 gate when cysteines present);
  5. the leg_inputs two-copy-only wiring guard on serialize_inplace_rbfe_system.

Byte-identical regression: the existing named specs (MTR/V3I/A9G/W4A) classify
UNCHANGED (the new flag defaults False).

Fast: no real OpenMM build (the full build + all C4/MC1/C6/MC2/MC3/R2 asserts +
endpoint-equivalence are exercised by scripts/trackb_fold_thermocycle_smoke.py).
"""

import dataclasses
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


def _load_ats():
    spec = importlib.util.spec_from_file_location(
        "atm_trackB_setup", os.path.join(_UTILS, "atm_trackB_setup.py"))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate utils/atm_trackB_setup.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_rbfe():
    spec = importlib.util.spec_from_file_location(
        "atm_trackB_inplace_rbfe",
        os.path.join(_UTILS, "atm_trackB_inplace_rbfe.py"))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate utils/atm_trackB_inplace_rbfe.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# 1) shape classification
# ---------------------------------------------------------------------------
def test_ile_ala_shape_is_acyclic_connected_group():
    ats = _load_ats()
    s = ats.make_ile_ala_mutation_spec(96)
    assert s.shape == "acyclic_connected_group"


def test_ile_ala_uncertified_is_unsupported():
    """Default flag (chained_group_certified=False) == pre-patch behaviour."""
    ats = _load_ats()
    s = dataclasses.replace(ats.make_ile_ala_mutation_spec(96),
                            chained_group_certified=False)
    assert s.shape == "unsupported"


def test_chained_flag_never_certifies_a_ring():
    """A ring spec (non-empty ring_closure_bonds) is NOT rescued by the chained
    flag — the ring gate is FIRST and requires connected_group_certified."""
    ats = _load_ats()
    s = dataclasses.replace(
        ats.make_ile_ala_mutation_spec(96),
        ring_closure_bonds=(("CG1", "CD1"),),   # a spurious ring on the var group
        chained_group_certified=True, connected_group_certified=False)
    assert s.shape == "unsupported"


def test_existing_specs_shape_byte_identical():
    """The new flag must not perturb any existing named spec's classification."""
    ats = _load_ats()
    assert ats.MUTATION_MTR_TRP_RES4.shape == "appearing_heavy"
    assert ats.MUTATION_VAL_ILE_RES3.shape == "appearing_heavy"
    assert ats.MUTATION_ALA_GLY_RES9.shape == "disappearing_heavy"
    assert ats.MUTATION_TRP_ALA_RES4.shape == "single_attach_connected_group"
    # default flag on all registry specs
    for spec in ats.MUTATION_SPECS.values():
        assert spec.chained_group_certified is False


# ---------------------------------------------------------------------------
# 2) make_ile_ala_mutation_spec SSOT factory
# ---------------------------------------------------------------------------
def test_ile_ala_factory_partition():
    ats = _load_ats()
    s = ats.make_ile_ala_mutation_spec(96, name="ile96ala")
    assert s.name == "ile96ala"
    assert s.resnum == 96
    assert s.common_attach_atom == "CB"
    assert s.stateA_resname == "ILE" and s.stateB_resname == "ALA"
    # disappearing: HB + the 3 side-chain heavies + their H's (12 atoms).
    assert set(s.stateA_only_atoms) == {
        "HB", "CG1", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23",
        "CD1", "HD11", "HD12", "HD13"}
    # appearing: the two extra ALA CB hydrogens (no appearing heavy).
    assert set(s.stateB_only_atoms) == {"HB1", "HB2", "HB3"}
    # heavy counts drive the shape: 3 disappearing heavies, 0 appearing.
    assert len(s._heavy_names(s.stateA_only_atoms)) == 3
    assert len(s._heavy_names(s.stateB_only_atoms)) == 0
    assert s.ring_closure_bonds == ()
    assert s.chained_group_certified is True
    assert s.multiheavy_star_certified is False
    assert s.connected_group_certified is False
    # net-charge-neutral canonical amber (no ncAA XML).
    assert s.hybrid_xml is None


def test_ile_ala_factory_default_name():
    ats = _load_ats()
    assert ats.make_ile_ala_mutation_spec(2).name == "ile_ala_res2"


def test_ile_ala_factory_resolves_via_resolve_mutation_spec():
    ats = _load_ats()
    s = ats.make_ile_ala_mutation_spec(96)
    assert ats.resolve_mutation_spec(s) is s   # instance returned unchanged


# ---------------------------------------------------------------------------
# 3) resolve_fold_leg_inputs
# ---------------------------------------------------------------------------
def test_resolve_fold_leg_inputs_schema(tmp_path):
    ats = _load_ats()
    scaffold = tmp_path / "scaffold.pdb"
    scaffold.write_text("END\n")
    li = ats.resolve_fold_leg_inputs(str(scaffold))
    assert set(li.keys()) == {"bound", "free", "final", "hydrogens_xml"}
    # BOTH slots point at the single real scaffold (cp4 only satisfies the guard).
    assert li["final"]["wt"] == str(scaffold)
    assert li["final"]["cp4"] == str(scaffold)
    assert li["hydrogens_xml"] is None


def test_resolve_fold_leg_inputs_hydrogens_xml_passthrough(tmp_path):
    ats = _load_ats()
    scaffold = tmp_path / "s.pdb"
    scaffold.write_text("END\n")
    li = ats.resolve_fold_leg_inputs(str(scaffold), hydrogens_xml="/x/h.xml")
    assert li["hydrogens_xml"] == "/x/h.xml"


def test_resolve_fold_leg_inputs_missing_file_raises():
    ats = _load_ats()
    with pytest.raises(FileNotFoundError):
        ats.resolve_fold_leg_inputs("/no/such/scaffold.pdb")


# ---------------------------------------------------------------------------
# 4) assert_twocopy_disulfides cysteine-conditional guard
# ---------------------------------------------------------------------------
def _make_fused_build_with_residue(resname, atom_names, disulfides):
    """Minimal fused_build stub carrying only a topology (+ disulfides). The
    0-CYS skip path returns before touching the System; the CYS-present <2 path
    raises before touching it — so system=None is sufficient for both."""
    from openmm import app

    top = app.Topology()
    chain = top.addChain("A")
    res = top.addResidue(resname, chain)
    for an in atom_names:
        el = app.element.sulfur if an == "SG" else app.element.carbon
        top.addAtom(an, el, res)

    class _Modeller:
        def __init__(self, topology):
            self.topology = topology

    return {"modeller": _Modeller(top), "system": None, "disulfides": disulfides}


@pytest.mark.skipif(not _have_openmm(), reason="needs openmm")
def test_mc3_skips_disulfide_free_scaffold():
    ats = _load_ats()
    fb = _make_fused_build_with_residue("ALA", ["N", "CA", "CB", "C", "O"], [])
    out = ats.assert_twocopy_disulfides(fb)
    assert out["passed"] is True
    assert out["n_disulfides"] == 0
    assert out["skipped_reason"] == "disulfide_free_scaffold_no_cysteines"


@pytest.mark.skipif(not _have_openmm(), reason="needs openmm")
def test_mc3_raises_when_cysteine_present_but_underdetected():
    ats = _load_ats()
    # A CYS SG is present but zero disulfides detected -> the >=2 gate must still
    # fire (byte-identical to the pre-patch behaviour for disulfide-bearing boxes).
    fb = _make_fused_build_with_residue("CYS", ["N", "CA", "CB", "SG", "C", "O"], [])
    with pytest.raises(ValueError, match="MC3 two-copy FAIL"):
        ats.assert_twocopy_disulfides(fb)


# ---------------------------------------------------------------------------
# 5) leg_inputs wiring guard (two-copy only)
# ---------------------------------------------------------------------------
def test_leg_inputs_rejected_on_single_core():
    rbfe = _load_rbfe()
    with pytest.raises(ValueError, match="leg_inputs is only supported"):
        rbfe.serialize_inplace_rbfe_system(
            leg="free", construction="single_core",
            leg_inputs={"final": {"wt": "/x.pdb", "cp4": "/x.pdb"}})


# ---------------------------------------------------------------------------
# 6) P3-#116 FIX2: deterministic + bounded-retry appearing-H placement.
#    PLACEMENT-ONLY — these tests exercise the seed derivation, the RNG seeding
#    primitive, and the bounded retry / fail-loud CONTROL FLOW with the real build
#    STUBBED OUT (no OpenMM build), so the FE core is provably never entered.
# ---------------------------------------------------------------------------
import hashlib as _hashlib   # noqa: E402
import random as _random     # noqa: E402


def test_derive_appearing_h_seed_deterministic_and_hashlib():
    """The per-unit seed is (a) reproducible for a fixed (unit_key, attempt), (b)
    equal to an INDEPENDENT hashlib.sha256 computation (NOT the process-salted builtin
    hash — that is the whole point: build-ORDER independence + cross-process
    reproducibility), and (c) in numpy's [0, 2**32) range."""
    ats = _load_ats()
    key = "s7|free|ile_ala_folded"
    s0a = ats._derive_appearing_h_seed(key, 0)
    s0b = ats._derive_appearing_h_seed(key, 0)
    assert s0a == s0b                         # reproducible
    # Independent re-derivation with the same recipe (sha256, first 8 hex, 31-bit).
    want = int(_hashlib.sha256(("%s#%d" % (key, 0)).encode("utf-8")).hexdigest()[:8],
               16) & 0x7FFFFFFF
    assert s0a == want                        # sha256, NOT builtin hash
    assert 0 <= s0a < 2 ** 32                 # valid np.random.seed range


def test_derive_appearing_h_seed_increments_per_attempt():
    """Successive attempts yield DISTINCT seeds (the retry re-draws, it does not
    re-try the identical rotamer)."""
    ats = _load_ats()
    key = "s127|free|ile_ala_folded"
    seeds = [ats._derive_appearing_h_seed(key, a) for a in range(5)]
    assert len(set(seeds)) == 5


def test_seed_appearing_h_placement_makes_rng_reproducible():
    """Seeding the process RNG makes the addHydrogens-consumed random.random()
    sequence reproducible (this is what removes the build-ORDER dependence)."""
    ats = _load_ats()
    ats._seed_appearing_h_placement(20260703)
    a = [_random.random() for _ in range(4)]
    ats._seed_appearing_h_placement(20260703)
    b = [_random.random() for _ in range(4)]
    assert a == b


def test_is_r2_seed_failure_matches_only_r2():
    ats = _load_ats()
    r2 = ValueError(
        "R2 two-copy seed min-dist FAIL: copy-1 appearing var group too close "
        "to a copy-1 common/solvent atom (0.0505 nm <= 0.1000 nm) — ('HB1','HD23').")
    other = ValueError("MC1 two-copy NET-charge sanity FAIL: ...")
    assert ats._is_r2_seed_failure(r2) is True
    assert ats._is_r2_seed_failure(other) is False


def _stub_r2_error():
    return ValueError(
        "R2 two-copy seed min-dist FAIL: copy-1 appearing var group too close "
        "to a copy-1 common/solvent atom (0.0505 nm <= 0.1000 nm) — ('HB1','HD23').")


def test_r2_retry_passes_first_attempt(monkeypatch):
    """A healthy first draw returns after ONE build, seeded deterministically, with
    the retry trail attached (no re-placement)."""
    ats = _load_ats()
    calls = []

    def _stub_build(**kwargs):
        calls.append(kwargs.get("appearing_h_seed"))
        return {"outcome": "twocopy_attached"}

    monkeypatch.setattr(ats, "build_inplace_res4_twocopy_system", _stub_build)
    out = ats.build_inplace_res4_twocopy_system_r2_retry(
        unit_key="s7|free|ile_ala_folded", retry_k=5, leg="free")
    assert len(calls) == 1
    assert calls[0] == ats._derive_appearing_h_seed("s7|free|ile_ala_folded", 0)
    trail = out["appearing_h_retry"]
    assert trail["attempts_used"] == 1
    assert trail["retry_k"] == 5
    assert trail["seed_used"] == calls[0]
    assert trail["seeds_tried"] == [calls[0]]


def test_r2_retry_recovers_after_bad_draw(monkeypatch):
    """A bad first draw (R2 FAIL) is re-placed with an INCREMENTED deterministic
    seed and recovered on the next attempt."""
    ats = _load_ats()
    seen = []

    def _stub_build(**kwargs):
        seen.append(kwargs.get("appearing_h_seed"))
        if len(seen) == 1:
            raise _stub_r2_error()          # first draw is the pathological rotamer
        return {"outcome": "twocopy_attached"}

    monkeypatch.setattr(ats, "build_inplace_res4_twocopy_system", _stub_build)
    out = ats.build_inplace_res4_twocopy_system_r2_retry(
        unit_key="s127|free|ile_ala_folded", retry_k=5, leg="free")
    assert len(seen) == 2
    key = "s127|free|ile_ala_folded"
    assert seen == [ats._derive_appearing_h_seed(key, 0),
                    ats._derive_appearing_h_seed(key, 1)]
    assert seen[0] != seen[1]               # seed incremented (fresh draw)
    assert out["appearing_h_retry"]["attempts_used"] == 2


def test_r2_retry_k_exhaustion_fails_loud(monkeypatch):
    """V3: every draw fails R2 -> the wrapper raises RuntimeError (fail-loud) after
    EXACTLY K attempts. The 0.10 nm threshold is NEVER relaxed and NO clashing build
    is returned."""
    ats = _load_ats()
    n = {"calls": 0}

    def _stub_build(**kwargs):
        n["calls"] += 1
        raise _stub_r2_error()

    monkeypatch.setattr(ats, "build_inplace_res4_twocopy_system", _stub_build)
    with pytest.raises(RuntimeError, match="FAILED on ALL 5 deterministic"):
        ats.build_inplace_res4_twocopy_system_r2_retry(
            unit_key="sX|free|ile_ala_folded", retry_k=5, leg="free")
    assert n["calls"] == 5                   # bounded — no infinite retry


def test_r2_retry_non_r2_error_is_not_swallowed(monkeypatch):
    """A NON-R2 build error (a genuine defect) is re-raised immediately — NOT
    retried away, NOT converted into the K-exhaustion RuntimeError."""
    ats = _load_ats()
    n = {"calls": 0}

    def _stub_build(**kwargs):
        n["calls"] += 1
        raise ValueError("MC1 two-copy NET-charge sanity FAIL: genuine defect")

    monkeypatch.setattr(ats, "build_inplace_res4_twocopy_system", _stub_build)
    with pytest.raises(ValueError, match="MC1 two-copy NET-charge sanity FAIL"):
        ats.build_inplace_res4_twocopy_system_r2_retry(
            unit_key="sX|free|ile_ala_folded", retry_k=5, leg="free")
    assert n["calls"] == 1                   # failed fast, no retry


def test_r2_retry_rejects_bad_k():
    ats = _load_ats()
    with pytest.raises(ValueError, match="retry_k must be >= 1"):
        ats.build_inplace_res4_twocopy_system_r2_retry(
            unit_key="s7|free|x", retry_k=0, leg="free")


def test_r2_retry_rejects_explicit_appearing_h_seed():
    """appearing_h_seed is derived per-attempt; passing it explicitly is a wiring
    error (the caller must not pin the seed and defeat the retry increment)."""
    ats = _load_ats()
    with pytest.raises(ValueError, match="appearing_h_seed is"):
        ats.build_inplace_res4_twocopy_system_r2_retry(
            unit_key="s7|free|x", retry_k=5, leg="free", appearing_h_seed=1)


class _RetryPath(Exception):
    pass


class _BasePath(Exception):
    pass


def test_serialize_routes_retry_only_when_opt_in(monkeypatch):
    """V2 opt-in gating: appearing_h_retry_k set -> the retry wrapper; None (default)
    -> the base build (byte-identical legacy path). Also verifies the per-unit key
    derivation (seed|leg|mutation-name)."""
    rbfe = _load_rbfe()
    captured = {}

    def _spy_retry(**kwargs):
        captured.update(kwargs)
        raise _RetryPath()

    def _spy_base(**kwargs):
        raise _BasePath()

    monkeypatch.setattr(rbfe.ats, "build_inplace_res4_twocopy_system_r2_retry",
                        _spy_retry)
    monkeypatch.setattr(rbfe.ats, "build_inplace_res4_twocopy_system", _spy_base)

    # opt-in -> retry wrapper, with the derived unit_key + retry_k.
    with pytest.raises(_RetryPath):
        rbfe._serialize_twocopy_system(
            leg="free", out_dir="/tmp/x", tag="free", seed="s7",
            binder_chain="B", solvate=True, harmonize_common_charges=False,
            displacement_nm=4.0, mtr_ncaa_xml=None, constraints=None,
            mutation_spec=None, appearing_h_retry_k=5)
    assert captured["unit_key"] == "s7|free|default"
    assert captured["retry_k"] == 5

    # default (None) -> base build, byte-identical legacy path (no retry entered).
    with pytest.raises(_BasePath):
        rbfe._serialize_twocopy_system(
            leg="free", out_dir="/tmp/x", tag="free", seed="s7",
            binder_chain="B", solvate=True, harmonize_common_charges=False,
            displacement_nm=4.0, mtr_ncaa_xml=None, constraints=None,
            mutation_spec=None, appearing_h_retry_k=None)


def test_appearing_h_retry_k_rejected_on_single_core():
    rbfe = _load_rbfe()
    with pytest.raises(ValueError, match="appearing_h_retry_k is only supported"):
        rbfe.serialize_inplace_rbfe_system(
            leg="free", construction="single_core", appearing_h_retry_k=5)
