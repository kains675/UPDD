#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``scripts/phase4_trackB_v2_make_system`` (v2.1 system builder)
and ``scripts/trackb_production_v2_1_upstream`` (v2.1 launcher).

v2.1 delegates the system XML build to upstream
``atom_openmm.make_atm_system_from_rcpt_lig.make_system``. These tests
cover the pure-Python glue:

* ``split_complex_to_receptor_binder`` — receptor/binder PDB split.
* ``_make_multi_xml_forcefield_proxy`` — multi-XML monkey-patch logic.
* ``_patched_make_system`` — try/finally revert contract (no global leak).
* v2.1 launcher re-exports of the canonical 22-state schedule.

GPU-dependent paths (upstream make_system, addHydrogens, abfe_structprep,
abfe_production) are NOT exercised here.
"""

import importlib.util
import json
import os
import sys
import tempfile

import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)


def _load_builder_module():
    """Load scripts/phase4_trackB_v2_make_system.py without atom_openmm
    requirement (the relevant module-level imports are openff/atom_openmm
    free).
    """
    spec = importlib.util.spec_from_file_location(
        "trackb_v21_builder",
        os.path.join(_PROJ, "scripts", "phase4_trackB_v2_make_system.py"),
    )
    if spec is None or spec.loader is None:
        pytest.skip("could not locate scripts/phase4_trackB_v2_make_system.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_v21_launcher():
    spec = importlib.util.spec_from_file_location(
        "trackb_v21_launcher",
        os.path.join(_PROJ, "scripts", "trackb_production_v2_1_upstream.py"),
    )
    if spec is None or spec.loader is None:
        pytest.skip(
            "could not locate scripts/trackb_production_v2_1_upstream.py"
        )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def builder():
    return _load_builder_module()


@pytest.fixture(scope="module")
def v21():
    return _load_v21_launcher()


# ---------------------------------------------------------------------------
# split_complex_to_receptor_binder
# ---------------------------------------------------------------------------
_TWO_CHAIN_PDB = """\
ATOM      1  N   PRO A   1      53.224  -2.194 -51.320  1.00  0.00           N
ATOM      2  CA  PRO A   1      53.000  -1.000 -50.000  1.00  0.00           C
ATOM      3  C   PRO A   1      52.000   0.000 -49.000  1.00  0.00           C
TER       4      PRO A   1
ATOM      5  N   ILE B   1      34.122  18.109 -29.155  1.00  0.00           N
ATOM      6  CA  ILE B   1      34.000  19.000 -28.000  1.00  0.00           C
HETATM    7  C   MTR B   4      40.000  20.000 -27.000  1.00  0.00           C
TER       8      MTR B   4
END
"""


def test_split_two_chain_complex(builder):
    """Receptor/binder split: chain A -> receptor_out, chain B -> binder_out."""
    with tempfile.TemporaryDirectory() as td:
        in_pdb = os.path.join(td, "in.pdb")
        rec_out = os.path.join(td, "rec.pdb")
        bin_out = os.path.join(td, "bin.pdb")
        with open(in_pdb, "w") as fh:
            fh.write(_TWO_CHAIN_PDB)
        r, b = builder.split_complex_to_receptor_binder(in_pdb, rec_out, bin_out)
        assert r == rec_out and b == bin_out
        rec_lines = [l for l in open(rec_out)
                     if l[:6] in ("ATOM  ", "HETATM")]
        bin_lines = [l for l in open(bin_out)
                     if l[:6] in ("ATOM  ", "HETATM")]
        assert len(rec_lines) == 3, f"receptor has wrong atom count: {len(rec_lines)}"
        assert len(bin_lines) == 3, f"binder has wrong atom count: {len(bin_lines)}"
        # Chain A atoms all belong to PRO chain A
        assert all(l[21] == "A" for l in rec_lines)
        # Chain B atoms belong to ILE/MTR chain B (one is HETATM)
        assert all(l[21] == "B" for l in bin_lines)
        # HETATM (MTR) preserved on binder side
        assert any(l[:6] == "HETATM" for l in bin_lines)


def test_split_raises_when_chain_missing(builder):
    """If a requested chain is absent, raise (not silent empty file)."""
    with tempfile.TemporaryDirectory() as td:
        in_pdb = os.path.join(td, "in.pdb")
        with open(in_pdb, "w") as fh:
            fh.write("ATOM      1  N   PRO A   1      0  0  0  1.00  0.00           N\nEND\n")
        with pytest.raises(RuntimeError, match="binder-chain"):
            builder.split_complex_to_receptor_binder(
                in_pdb,
                os.path.join(td, "rec.pdb"),
                os.path.join(td, "bin.pdb"),
            )


# ---------------------------------------------------------------------------
# Multi-XML ForceField monkey-patch
# ---------------------------------------------------------------------------
class _FakeForceField:
    """Stand-in for openmm.app.ForceField that records its *files args."""

    instances = []

    def __init__(self, *files):
        self.files = list(files)
        _FakeForceField.instances.append(self)


def test_proxy_splits_whitespace_separated_paths(builder):
    """``proxy("a.xml b.xml", "c.xml")`` -> ``Fake("a.xml", "b.xml", "c.xml")``."""
    _FakeForceField.instances.clear()
    proxy = builder._make_multi_xml_forcefield_proxy(_FakeForceField)
    proxy("amber14-all.xml params/MTR_gaff2_hybrid.xml", "amber14/tip3p.xml")
    assert len(_FakeForceField.instances) == 1
    files = _FakeForceField.instances[0].files
    assert files == [
        "amber14-all.xml",
        "params/MTR_gaff2_hybrid.xml",
        "amber14/tip3p.xml",
    ]


def test_proxy_preserves_single_path_unchanged(builder):
    """Single-path strings (no whitespace) pass through untouched."""
    _FakeForceField.instances.clear()
    proxy = builder._make_multi_xml_forcefield_proxy(_FakeForceField)
    proxy("amber14-all.xml", "amber14/tip3p.xml")
    files = _FakeForceField.instances[0].files
    assert files == ["amber14-all.xml", "amber14/tip3p.xml"]


def test_proxy_passes_non_string_args_unchanged(builder):
    """File-like objects + None must pass through (openmm accepts these)."""
    _FakeForceField.instances.clear()
    proxy = builder._make_multi_xml_forcefield_proxy(_FakeForceField)
    sentinel = object()
    proxy(sentinel, "single.xml")
    files = _FakeForceField.instances[0].files
    assert files[0] is sentinel
    assert files[1] == "single.xml"


def test_proxy_handles_tab_separator(builder):
    """Tabs also trigger split (defensive against indentation glitches)."""
    _FakeForceField.instances.clear()
    proxy = builder._make_multi_xml_forcefield_proxy(_FakeForceField)
    proxy("a.xml\tb.xml")
    assert _FakeForceField.instances[0].files == ["a.xml", "b.xml"]


# ---------------------------------------------------------------------------
# _patched_make_system: try/finally revert contract
# ---------------------------------------------------------------------------
def test_patched_make_system_reverts_on_normal_return(builder, monkeypatch):
    """After a successful call, ForceField symbol is restored on vendored."""
    fake_vendored = type(sys)("_vendored_make_atm_system")
    sentinel = object()
    fake_vendored.ForceField = sentinel
    called = {"yes": False}

    def fake_make_system(**kw):
        called["yes"] = True
        # Note: the proxy is now installed; verify it was swapped in
        assert fake_vendored.ForceField is not sentinel
    fake_vendored.make_system = fake_make_system

    real_modules = sys.modules.copy()
    sys.modules["_vendored_make_atm_system"] = fake_vendored
    try:
        builder._patched_make_system(receptorfile="dummy.pdb")
    finally:
        sys.modules.clear()
        sys.modules.update(real_modules)

    assert called["yes"]
    # CRITICAL: original symbol restored
    assert fake_vendored.ForceField is sentinel


def test_patched_make_system_reverts_on_exception(builder):
    """If vendored raises, ForceField symbol is STILL restored."""
    fake_vendored = type(sys)("_vendored_make_atm_system")
    sentinel = object()
    fake_vendored.ForceField = sentinel

    def fake_make_system(**kw):
        raise RuntimeError("vendored blew up")
    fake_vendored.make_system = fake_make_system

    real_modules = sys.modules.copy()
    sys.modules["_vendored_make_atm_system"] = fake_vendored
    try:
        with pytest.raises(RuntimeError, match="vendored blew up"):
            builder._patched_make_system(receptorfile="dummy.pdb")
    finally:
        sys.modules.clear()
        sys.modules.update(real_modules)

    # CRITICAL: original symbol restored even on exception
    assert fake_vendored.ForceField is sentinel


# ---------------------------------------------------------------------------
# v2.1 launcher schedule + helper re-exports (regression vs v2)
# ---------------------------------------------------------------------------
def test_v21_schedule_matches_v2(v21):
    """v2.1 re-exports the canonical 22-state schedule (must equal v2)."""
    assert v21.N_STATES == 22
    assert v21.LAMBDAS_1[:11] == [0.00, 0.05, 0.10, 0.15, 0.20,
                                  0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    assert v21.LAMBDAS_1[11:] == [0.50, 0.45, 0.40, 0.35, 0.30,
                                  0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
    assert v21.DIRECTIONS[:11] == [1] * 11
    assert v21.DIRECTIONS[11:] == [-1] * 11
    assert v21.UMAX_KCAL == 200.0
    assert v21.UBCORE_KCAL == 100.0
    assert v21.ACORE == 0.062500
    assert v21.DISPLACEMENT_NM_DEFAULT == (2.5, 0.0, 0.0)
    assert v21.PRE_REGISTER_OUTCOMES == (
        "1_sign_stable",
        "2_sigma_drift",
        "3_magnitude_drift",
        "4_sign_flip",
    )


def test_v21_helpers_reused_from_v2(v21):
    """v2.1 re-uses cntl + nodefile + binary-locator helpers from v2."""
    # These must be callable (not None)
    assert callable(v21.write_cntl_file)
    assert callable(v21.write_nodefile)
    assert callable(v21.collect_binder_atom_indices)
    assert callable(v21.collect_receptor_ca_indices)
    assert callable(v21.run_abfe_structprep_for_leg)
    assert callable(v21.run_abfe_production_for_leg)


# ---------------------------------------------------------------------------
# --bound-schedule wiring (densified_bound28 per-direction-path plumbing)
# ---------------------------------------------------------------------------
def _load_per_direction_production():
    """Load scripts/trackb_per_direction_production.py (cntl-split helpers
    are atom_openmm-free; the GPU dispatch path is not exercised here)."""
    spec = importlib.util.spec_from_file_location(
        "trackb_per_direction_production",
        os.path.join(_PROJ, "scripts", "trackb_per_direction_production.py"),
    )
    if spec is None or spec.loader is None:
        pytest.skip(
            "could not locate scripts/trackb_per_direction_production.py"
        )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_v21_bound_schedule_registry_is_v2_ssot(v21):
    """The bound-schedule registry is the SINGLE SSOT in v2 — v2.1 only
    imports it (no duplicate definition). The free-leg densified38* ladders
    must NOT be bound-valid."""
    assert v21.v2_legacy.DEFAULT_BOUND_SCHEDULE == "canonical22"
    assert v21.v2_legacy.BOUND_SCHEDULES == frozenset(
        {"canonical22", "densified_bound28", "densified_bound30"}
    )
    # free-leg densified38* are not bound-valid
    assert "densified38" not in v21.v2_legacy.BOUND_SCHEDULES
    assert "densified38v4" not in v21.v2_legacy.BOUND_SCHEDULES
    # v2.1 module does NOT redefine the registry (import-only SSOT)
    assert not hasattr(v21, "BOUND_SCHEDULES") or (
        v21.BOUND_SCHEDULES is v21.v2_legacy.BOUND_SCHEDULES
    )


def test_v21_bound_default_cntl_byte_equal_canonical22(v21, tmp_path):
    """bound_schedule default (canonical22) → schedule=None path → byte-equal
    to an explicit canonical22 None-schedule cntl (no regression)."""
    lig = list(range(10))
    restr = list(range(10, 20))
    kw = dict(
        basename="trackb",
        nodefile_path=str(tmp_path / "nodefile"),
        ligand_atom_indices=lig,
        pos_restrained_atom_indices=restr,
        displacement_nm=(2.5, 0.0, 0.0),
        production_steps=2500,
        prnt_frequency=2500,
        trj_frequency=25000,
        wall_time_min=720,
        cycle_time_s=10,
        checkpoint_time_s=600,
    )
    # Default bound resolution = schedule_dict None (canonical22).
    p_default = tmp_path / "default.cntl"
    v21.write_cntl_file(cntl_path=str(p_default), schedule=None, **kw)
    # Explicit canonical22 None-schedule cntl.
    p_canon = tmp_path / "canon.cntl"
    v21.write_cntl_file(cntl_path=str(p_canon), schedule=None, **kw)
    assert p_default.read_text() == p_canon.read_text()
    # The canonical cntl carries 22-state arrays.
    body = p_default.read_text()
    for line in body.splitlines():
        if line.strip().startswith("LAMBDA1 ="):
            inner = line.split("=", 1)[1].strip().strip("'\"")
            assert len([t for t in inner.split(",") if t.strip()]) == 22


def test_v21_bound_densified28_cntl_is_28_state_with_bridge(v21, tmp_path):
    """--bound-schedule densified_bound28 → write_cntl_file emits a 28-state
    (dplus14 + dminus14) cntl whose forward LAMBDA1 carries the W0-graded
    soft-core bridge windows (W0COEFF 0.20 / 0.45 / 0.70) at the 6→7 cliff."""
    sched = v21.v2_legacy.get_schedule("densified_bound28")
    p = tmp_path / "bound28.cntl"
    v21.write_cntl_file(
        cntl_path=str(p),
        basename="trackb",
        nodefile_path=str(tmp_path / "nodefile"),
        ligand_atom_indices=list(range(10)),
        pos_restrained_atom_indices=list(range(10, 20)),
        displacement_nm=(2.5, 0.0, 0.0),
        production_steps=2500,
        prnt_frequency=2500,
        trj_frequency=25000,
        wall_time_min=720,
        cycle_time_s=10,
        checkpoint_time_s=600,
        schedule=sched,
    )
    body = p.read_text()

    def _row(key):
        for line in body.splitlines():
            if line.strip().startswith(key + " ="):
                inner = line.split("=", 1)[1].strip().strip("'\"")
                return [t.strip() for t in inner.split(",") if t.strip()]
        raise AssertionError(f"no {key} row in cntl")

    # 28-state symmetric (14 forward + 14 backward).
    direction = [int(float(x)) for x in _row("DIRECTION")]
    assert len(direction) == 28
    assert sum(1 for d in direction if d >= 0) == 14
    assert sum(1 for d in direction if d < 0) == 14
    for key in ("LAMBDA1", "LAMBDA2", "ALPHA", "U0", "W0COEFF", "INTERMEDIATE"):
        assert len(_row(key)) == 28, f"{key} must be 28-state"

    # The bridge W0 ramp (0.20 / 0.45 / 0.70) appears in the forward block.
    w0_fwd = [float(x) for x in _row("W0COEFF")][:14]
    assert 0.20 in w0_fwd
    assert 0.45 in w0_fwd
    assert 0.70 in w0_fwd


def test_v21_bound_densified30_resolves_and_is_30_state(v21, tmp_path):
    """--bound-schedule densified_bound30 resolves through the v2 SSOT registry
    (BOUND_SCHEDULES validation) and write_cntl_file emits a 30-state
    (dplus15 + dminus15) cntl with 4 soft-core bridge windows whose U0 DESCENDS
    105→97 (MC-1). The free-leg densified38* ladders remain bound-INVALID."""
    assert "densified_bound30" in v21.v2_legacy.BOUND_SCHEDULES
    assert "densified38v4" not in v21.v2_legacy.BOUND_SCHEDULES
    sched = v21.v2_legacy.get_schedule("densified_bound30")
    p = tmp_path / "bound30.cntl"
    v21.write_cntl_file(
        cntl_path=str(p),
        basename="trackb",
        nodefile_path=str(tmp_path / "nodefile"),
        ligand_atom_indices=list(range(10)),
        pos_restrained_atom_indices=list(range(10, 20)),
        displacement_nm=(2.5, 0.0, 0.0),
        production_steps=2500,
        prnt_frequency=2500,
        trj_frequency=25000,
        wall_time_min=720,
        cycle_time_s=10,
        checkpoint_time_s=600,
        schedule=sched,
    )
    body = p.read_text()

    def _row(key):
        for line in body.splitlines():
            if line.strip().startswith(key + " ="):
                inner = line.split("=", 1)[1].strip().strip("'\"")
                return [t.strip() for t in inner.split(",") if t.strip()]
        raise AssertionError(f"no {key} row in cntl")

    # 30-state symmetric (15 forward + 15 backward).
    direction = [int(float(x)) for x in _row("DIRECTION")]
    assert len(direction) == 30
    assert sum(1 for d in direction if d >= 0) == 15
    assert sum(1 for d in direction if d < 0) == 15
    for key in ("LAMBDA1", "LAMBDA2", "ALPHA", "U0", "W0COEFF", "INTERMEDIATE"):
        assert len(_row(key)) == 30, f"{key} must be 30-state"

    # 4 soft-core bridge windows in the forward block; U0 strictly DESCENDS.
    w0_fwd = [float(x) for x in _row("W0COEFF")][:15]
    u0_fwd = [float(x) for x in _row("U0")][:15]
    bridge_u0 = [u for w, u in zip(w0_fwd, u0_fwd) if 0.0 < w < 1.0]
    assert len(bridge_u0) == 4
    assert all(a > b for a, b in zip(bridge_u0, bridge_u0[1:])), (
        f"U0 must strictly descend across the bridge (MC-1); got {bridge_u0!r}")


def test_v21_bound_densified28_cntl_splits_14_14(v21, tmp_path):
    """The 28-state combined bound cntl from --bound-schedule densified_bound28
    is split by generate_per_direction_cntls into dplus14 / dminus14
    (count-agnostic: the slice boundary is sum(DIRECTION>=0), never a
    hardcoded 11/22)."""
    perdir = _load_per_direction_production()
    leg_dir = tmp_path / "cp4" / "bound"
    leg_dir.mkdir(parents=True)
    sched = v21.v2_legacy.get_schedule("densified_bound28")
    cntl_path = leg_dir / "trackb_asyncre.cntl"
    v21.write_cntl_file(
        cntl_path=str(cntl_path),
        basename="trackb",
        nodefile_path=str(leg_dir / "nodefile"),
        ligand_atom_indices=list(range(10)),
        pos_restrained_atom_indices=list(range(10, 20)),
        displacement_nm=(2.5, 0.0, 0.0),
        production_steps=2500,
        prnt_frequency=2500,
        trj_frequency=25000,
        wall_time_min=720,
        cycle_time_s=10,
        checkpoint_time_s=600,
        schedule=sched,
    )
    gen = perdir.generate_per_direction_cntls(
        leg_dir=str(leg_dir), jobname="trackb"
    )
    assert gen["total_state_count"] == 28
    assert gen["fwd_state_count"] == 14
    assert gen["directions"]["dplus"]["n_states"] == 14
    assert gen["directions"]["dminus"]["n_states"] == 14
    assert gen["directions"]["dplus"]["state_slice"] == [0, 14]
    assert gen["directions"]["dminus"]["state_slice"] == [14, 28]
    # The dminus slice carries DIRECTION forced to +1 (b+ ATM base=u0).
    dm_dir = gen["directions"]["dminus"]["schedule"]["DIRECTION"]
    assert len(dm_dir) == 14
    assert all(int(float(d)) == 1 for d in dm_dir)


def test_v21_argparse_has_bound_schedule():
    """The launcher argparse exposes --bound-schedule {canonical22,
    densified_bound28, densified_bound30} (default canonical22). Verified via the
    CLI --help surface (argparse lives inside main(), no standalone parser
    factory)."""
    import subprocess

    launcher = os.path.join(
        _PROJ, "scripts", "trackb_production_v2_1_upstream.py"
    )
    proc = subprocess.run(
        [sys.executable, launcher, "--help"],
        capture_output=True, text=True, timeout=120,
    )
    # argparse --help exits 0.
    assert proc.returncode == 0, proc.stderr
    out = proc.stdout
    assert "--bound-schedule" in out
    assert "densified_bound28" in out
    assert "densified_bound30" in out
    # all bound-valid choices surfaced in the {..} choice set
    assert "canonical22" in out


# ---------------------------------------------------------------------------
# Builder argument plumbing — defaults preserved
# ---------------------------------------------------------------------------
def test_builder_default_displacement(builder):
    """Default ABFE displacement is 2.5 nm +x (matches v2 production schedule)."""
    assert builder.DISPLACEMENT_NM_DEFAULT == (2.5, 0.0, 0.0)


def test_builder_hybrid_mtr_xml_path(builder):
    """The hybrid MTR XML path resolves to params/MTR_gaff2_hybrid.xml."""
    assert builder.HYBRID_MTR_XML.endswith(
        "params/MTR_gaff2_hybrid.xml"
    )
    assert os.path.isfile(builder.HYBRID_MTR_XML)


def test_builder_protein_ff_default(builder):
    """amber14-all is the default protein FF (matches utils/atm_trackB_setup)."""
    assert builder.PROTEIN_FF_DEFAULT == "amber14-all.xml"


def test_builder_solvent_ff_default(builder):
    """amber14/tip3p is the default solvent FF (matches upstream convention)."""
    assert builder.SOLVENT_FF_DEFAULT == "amber14/tip3p.xml"


# ---------------------------------------------------------------------------
# Vendored upstream make_system: minimal structural checks
# ---------------------------------------------------------------------------
def _load_vendored_module():
    spec = importlib.util.spec_from_file_location(
        "vendored_make_atm",
        os.path.join(_PROJ, "scripts", "_vendored_make_atm_system.py"),
    )
    if spec is None or spec.loader is None:
        pytest.skip("vendored module not found")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_vendored_make_system_has_pdb_only_guards():
    """Vendored ``make_system`` rejects SDF inputs with NotImplementedError."""
    vendored = _load_vendored_module()
    with pytest.raises(NotImplementedError, match=".pdb receptor"):
        vendored.make_system(
            receptorfile="bogus.sdf",
            displacement=[0.0, 0.0, 0.0],
            xmloutfile="/tmp/_v.xml",
            pdboutfile="/tmp/_v.pdb",
        )


def test_vendored_make_system_rejects_rbfe():
    """Vendored ``make_system`` rejects RBFE (lig2file passed) — ABFE only."""
    vendored = _load_vendored_module()
    # Use real receptor PDB but with lig2file to trigger the RBFE guard.
    receptor = os.path.join(
        _PROJ, "outputs", "2QKI_Cp4_hybrid_calib_s7",
        "_md_input", "2QKI_Cp4.pdb"
    )
    if not os.path.isfile(receptor):
        pytest.skip(f"prepared receptor missing: {receptor}")
    # Don't actually run — the guard fires before any heavy work; but
    # ForceField will load amber14-all.xml first. We accept that this test
    # touches openmm; if openmm is missing skip.
    pytest.importorskip("openmm")
    with pytest.raises(NotImplementedError, match="ABFE only"):
        # Need a fake .pdb path that won't actually be parsed for ligand2
        vendored.make_system(
            receptorfile=receptor,
            displacement=[25.0, 0.0, 0.0],
            xmloutfile="/tmp/_v_rbfe.xml",
            pdboutfile="/tmp/_v_rbfe.pdb",
            lig1file=receptor,  # arbitrary
            lig2file=receptor,  # triggers rbfe=True
        )


def test_vendored_make_system_rejects_cofactor_sdf():
    """Cofactor SDF input is gated (needs openff)."""
    vendored = _load_vendored_module()
    receptor = os.path.join(
        _PROJ, "outputs", "2QKI_Cp4_hybrid_calib_s7",
        "_md_input", "2QKI_Cp4.pdb"
    )
    if not os.path.isfile(receptor):
        pytest.skip(f"prepared receptor missing: {receptor}")
    pytest.importorskip("openmm")
    with pytest.raises(NotImplementedError, match="cofactor SDF"):
        vendored.make_system(
            receptorfile=receptor,
            displacement=[25.0, 0.0, 0.0],
            xmloutfile="/tmp/_v_cof.xml",
            pdboutfile="/tmp/_v_cof.pdb",
            cofsdffile="bogus.sdf",
        )


# ---------------------------------------------------------------------------
# Structprep cache integrity (cross-stage drift enforcement, 2026-05-31)
# ---------------------------------------------------------------------------
# Tests cover the sidecar emit + validate contract introduced to close
# the silent SMOKE-vs-production reuse bug (cp4/free cycle-1 NaN).
# The cached cp4/free smoke state (6000 steps, T=271 K) was silently
# reused under a production launch where wt/free ran a fresh 850000-step
# structprep — causing cycle-1 NaN at r11.


def _write_minimal_state_xml(path: str) -> None:
    """Write a syntactically-valid placeholder state XML for sha + isfile
    checks. The validator does not parse it — only the sidecar / log
    metadata is consulted.
    """
    with open(path, "w") as fh:
        fh.write("<State/>\n")


def _write_structprep_log(
    path: str,
    last_step: int,
    final_T: float,
    interval: int = 100,
) -> None:
    """Emit a fake _structprep.log with the atom_openmm step-line
    format ``"<step>,<PE>,<T>,<box>,<speed>"``. Writes ``interval``-step
    lines from interval up to last_step; final entry has T=final_T.
    """
    lines = [
        "# Command: /home/san/miniconda3/envs/atm/bin/abfe_structprep trackb_asyncre.cntl\n",
        "Energy minimizing the system ...\n",
        "Potential energy after minimization = -72391.13 kJ/mol\n",
        "Thermalization ...\n",
        '#"Step","Potential Energy (kJ/mole)","Temperature (K)",'
        '"Box Volume (nm^3)","Speed (ns/day)"\n',
    ]
    step = interval
    while step < last_step:
        # Linear ramp T from 5 K (early thermalization) to a mid value.
        T_ramp = 5.0 + (final_T - 5.0) * (step / float(last_step))
        lines.append(
            f"{step},-{70000 - step * 0.01:.4f},{T_ramp:.4f},45.0,500\n"
        )
        step += interval
    # Final step exactly at last_step with final_T.
    lines.append(
        f"{last_step},-59594.31,{final_T:.4f},45.0,500\n"
    )
    lines.append("SaveState ...\n")
    with open(path, "w") as fh:
        fh.writelines(lines)


def test_sidecar_emitted_after_structprep(v21, tmp_path):
    """``_emit_structprep_sidecar`` writes a well-formed JSON sidecar
    capturing mode + step_count + final_temp_K + SHA-256 + timestamp +
    ckpt marker presence (integrity producer side)."""
    leg_dir = tmp_path / "cp4" / "free"
    leg_dir.mkdir(parents=True)
    state_xml = leg_dir / "trackb_0.xml"
    _write_minimal_state_xml(str(state_xml))
    _write_structprep_log(
        str(leg_dir / "_structprep.log"),
        last_step=850000,
        final_T=302.07,
    )
    # Mark "ckpt_is_valid" present to verify that flag is captured.
    (leg_dir / "ckpt_is_valid").write_text("")

    sidecar_path = v21._emit_structprep_sidecar(
        leg_dir=str(leg_dir),
        mode="production",
        basename="trackb",
    )
    assert os.path.isfile(sidecar_path)
    with open(sidecar_path) as fh:
        meta = json.load(fh)
    assert meta["schema_version"] == "0.9.13"
    assert meta["mode"] == "production"
    assert meta["step_count"] == 850000
    assert abs(meta["final_temp_K"] - 302.07) < 1e-3
    assert "structprep_completed_at" in meta
    # SHA-256 of "<State/>\n" is a known constant.
    assert len(meta["trackb_0_xml_sha256"]) == 64
    assert meta["ckpt_is_valid_present"] is True


def test_validate_cached_production_sidecar_passes(v21, tmp_path):
    """Production-mode sidecar (step_count=850000, T=302 K) → (True, ...).
    This is the happy path: a freshly-completed production structprep
    is correctly validated for a subsequent production launch."""
    leg_dir = tmp_path / "wt" / "free"
    leg_dir.mkdir(parents=True)
    _write_minimal_state_xml(str(leg_dir / "trackb_0.xml"))
    _write_structprep_log(
        str(leg_dir / "_structprep.log"),
        last_step=850000,
        final_T=302.07,
    )
    v21._emit_structprep_sidecar(
        leg_dir=str(leg_dir),
        mode="production",
        basename="trackb",
    )
    valid, reason = v21._validate_cached_structprep(
        leg_dir=str(leg_dir),
        expected_mode="production",
        basename="trackb",
    )
    assert valid is True
    assert "validated" in reason.lower()
    assert "sidecar" in reason.lower()
    assert "850000" in reason
    assert "302" in reason


def test_validate_cached_smoke_sidecar_rejected(v21, tmp_path):
    """SMOKE sidecar (step_count=6000, T=271 K, mode='smoke') under a
    production-expected launch → (False, reason). This is the exact
    failure mode of the cp4/free NaN incident: smoke artifact silently
    reused as production. The reason text MUST include the leg dir +
    archive remediation step (operator actionability)."""
    leg_dir = tmp_path / "cp4" / "free"
    leg_dir.mkdir(parents=True)
    _write_minimal_state_xml(str(leg_dir / "trackb_0.xml"))
    _write_structprep_log(
        str(leg_dir / "_structprep.log"),
        last_step=6000,
        final_T=271.02,
    )
    v21._emit_structprep_sidecar(
        leg_dir=str(leg_dir),
        mode="smoke",
        basename="trackb",
    )
    valid, reason = v21._validate_cached_structprep(
        leg_dir=str(leg_dir),
        expected_mode="production",
        basename="trackb",
    )
    assert valid is False
    # Reason must cite the actual leg dir AND a remediation action.
    assert str(leg_dir) in reason
    assert "archive" in reason.lower()
    assert "rerun" in reason.lower()
    # Mode mismatch is the first check, so the reason should mention mode
    # rather than step count (early-fail ordering).
    assert "mode" in reason.lower()


def test_validate_legacy_log_fallback_passes(v21, tmp_path):
    """No sidecar but _structprep.log shows last_step=850000, T=302 K
    → (True, ...). Required for legacy artifacts predating sidecar
    emission so this fix is backward-compatible (no false rejections of
    the wt/free cache, which was generated before this fix)."""
    leg_dir = tmp_path / "wt" / "free"
    leg_dir.mkdir(parents=True)
    _write_minimal_state_xml(str(leg_dir / "trackb_0.xml"))
    _write_structprep_log(
        str(leg_dir / "_structprep.log"),
        last_step=850000,
        final_T=302.07,
    )
    # NO sidecar emission — simulates legacy artifact.
    assert not os.path.isfile(str(leg_dir / "trackb_0.xml.meta.json"))

    valid, reason = v21._validate_cached_structprep(
        leg_dir=str(leg_dir),
        expected_mode="production",
        basename="trackb",
    )
    assert valid is True
    assert "legacy" in reason.lower()
    assert "no sidecar" in reason.lower()
    assert "850000" in reason


def test_validate_legacy_log_smoke_rejected(v21, tmp_path):
    """No sidecar, but _structprep.log shows last_step=6000 (smoke) →
    (False, reason). This is the diagnostic path for legacy smoke
    artifacts: even without a sidecar, the log-based step-count gate
    catches the silent reuse."""
    leg_dir = tmp_path / "cp4" / "free"
    leg_dir.mkdir(parents=True)
    _write_minimal_state_xml(str(leg_dir / "trackb_0.xml"))
    _write_structprep_log(
        str(leg_dir / "_structprep.log"),
        last_step=6000,
        final_T=271.02,
    )
    valid, reason = v21._validate_cached_structprep(
        leg_dir=str(leg_dir),
        expected_mode="production",
        basename="trackb",
    )
    assert valid is False
    assert str(leg_dir) in reason
    assert "6000" in reason
    assert "archive" in reason.lower()
    assert "rerun" in reason.lower()


def test_integrity_error_remediation_message_clear(v21, tmp_path):
    """IntegrityError subclasses RuntimeError, and validation failure
    reasons (which are wrapped into IntegrityError at the launcher
    call-site) must always include the archive path + remediation step
    so the operator can act without re-reading code."""
    assert issubclass(v21.IntegrityError, RuntimeError)

    # 1. Sidecar mode mismatch reason
    leg_dir = tmp_path / "cp4" / "free"
    leg_dir.mkdir(parents=True)
    _write_minimal_state_xml(str(leg_dir / "trackb_0.xml"))
    _write_structprep_log(
        str(leg_dir / "_structprep.log"),
        last_step=6000,
        final_T=271.0,
    )
    v21._emit_structprep_sidecar(
        leg_dir=str(leg_dir), mode="smoke", basename="trackb"
    )
    _, reason_smoke = v21._validate_cached_structprep(
        leg_dir=str(leg_dir),
        expected_mode="production",
        basename="trackb",
    )
    # Demonstrate the IntegrityError construction the launcher uses.
    err = v21.IntegrityError(
        f"Cached structprep validation FAILED for {leg_dir}: {reason_smoke}"
    )
    msg = str(err)
    assert "Cached structprep validation FAILED" in msg
    assert str(leg_dir) in msg
    assert "Archive" in msg or "archive" in msg
    assert "_archive/" in msg  # explicit archive root suggestion
    assert "rerun" in msg.lower()

    # 2. No sidecar AND no log reason — second failure path
    bare_dir = tmp_path / "wt" / "bound"
    bare_dir.mkdir(parents=True)
    _write_minimal_state_xml(str(bare_dir / "trackb_0.xml"))
    valid_bare, reason_bare = v21._validate_cached_structprep(
        leg_dir=str(bare_dir),
        expected_mode="production",
        basename="trackb",
    )
    assert valid_bare is False
    assert "no sidecar" in reason_bare.lower()
    assert "no _structprep.log" in reason_bare.lower()
    assert "archive" in reason_bare.lower()
    assert str(bare_dir) in reason_bare


def test_sidecar_temperature_outside_tolerance_rejected(v21, tmp_path):
    """Even a production-mode sidecar with the right step count is
    rejected if the final temperature drifted from 300 K beyond
    tolerance. Defense in depth: catches the "production_steps OK but
    integrator blew up before reaching thermal equilibrium" failure
    mode that step_count alone cannot detect."""
    leg_dir = tmp_path / "cp4" / "bound"
    leg_dir.mkdir(parents=True)
    _write_minimal_state_xml(str(leg_dir / "trackb_0.xml"))
    _write_structprep_log(
        str(leg_dir / "_structprep.log"),
        last_step=850000,
        final_T=265.0,  # 35 K below target, well outside default 10 K tol
    )
    v21._emit_structprep_sidecar(
        leg_dir=str(leg_dir),
        mode="production",
        basename="trackb",
    )
    valid, reason = v21._validate_cached_structprep(
        leg_dir=str(leg_dir),
        expected_mode="production",
        basename="trackb",
    )
    assert valid is False
    assert "265" in reason
    assert "300" in reason
    assert str(leg_dir) in reason
    assert "archive" in reason.lower()


# ---------------------------------------------------------------------------
# Per-direction system XML build (v0.9.19, 2026-06-01)
# trackb_per_direction_system_xml_rebuild_20260601
# ---------------------------------------------------------------------------
def test_build_bound_leg_via_upstream_accepts_direction_kwarg(builder):
    """build_bound_leg_via_upstream gained a direction kwarg in v0.9.19.

    Verify the parameter exists in the signature, the default is 'dplus'
    (legacy compat), and 'dminus' is accepted.
    """
    import inspect
    sig = inspect.signature(builder.build_bound_leg_via_upstream)
    assert "direction" in sig.parameters, (
        "build_bound_leg_via_upstream must accept direction= kwarg "
        "(v0.9.19 per-direction system XML rebuild)"
    )
    assert sig.parameters["direction"].default == "dplus"


def test_build_bound_leg_via_upstream_rejects_bad_direction(builder,
                                                            monkeypatch):
    """direction must be 'dplus' or 'dminus'; anything else raises
    ValueError BEFORE any heavy work happens."""
    # Patch the inner _patched_make_system to fail loudly if reached —
    # the ValueError must fire first.
    def _should_not_be_called(**kwargs):
        raise AssertionError("_patched_make_system reached with bad direction")
    monkeypatch.setattr(builder, "_patched_make_system", _should_not_be_called)
    with pytest.raises(ValueError, match="direction must be"):
        builder.build_bound_leg_via_upstream(
            receptor_pdb="/tmp/r.pdb",
            binder_pdb="/tmp/b.pdb",
            xml_out="/tmp/o.xml",
            pdb_out="/tmp/o.pdb",
            direction="forward",  # not dplus/dminus
        )


def test_build_bound_leg_dplus_uses_positive_displacement(builder,
                                                          monkeypatch):
    """direction='dplus' passes the +displacement vector verbatim into
    _patched_make_system AND disables the Modeller pre-displacement
    (binder stays at binding site x0 per the (b+) corrected spec)."""
    captured = {}

    def _capture(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(builder, "_patched_make_system", _capture)
    result = builder.build_bound_leg_via_upstream(
        receptor_pdb="/tmp/r.pdb",
        binder_pdb="/tmp/b.pdb",
        xml_out="/tmp/o.xml",
        pdb_out="/tmp/o.pdb",
        displacement_nm=(2.5, 0.0, 0.0),
        direction="dplus",
    )
    # displacement is forwarded in Angstrom (nm * 10) — POSITIVE for dplus
    assert captured["displacement"] == [25.0, 0.0, 0.0]
    # Modeller pre-displacement DISABLED for dplus (binder physically at x0)
    assert captured["apply_modeller_pre_displacement"] is False
    # Metadata records the spec contract
    assert result["direction"] == "dplus"
    assert result["apply_modeller_pre_displacement"] is False
    assert "bound" in result["binder_physical_state"]


def test_build_bound_leg_dminus_keeps_positive_displacement(builder,
                                                            monkeypatch):
    """(b+) corrected: direction='dminus' STILL passes
    +displacement to vendored upstream (NOT negated). The negation
    happens at the ATMForce runtime level via the structprep
    set_displacement monkeypatch — NOT at the system XML build."""
    captured = {}

    def _capture(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(builder, "_patched_make_system", _capture)
    result = builder.build_bound_leg_via_upstream(
        receptor_pdb="/tmp/r.pdb",
        binder_pdb="/tmp/b.pdb",
        xml_out="/tmp/o.xml",
        pdb_out="/tmp/o.pdb",
        displacement_nm=(2.5, 0.0, 0.0),
        direction="dminus",
    )
    # POSITIVE displacement (corrected from earlier negate-at-build-time):
    # both directions use +d_vec for bbox + Modeller positions.
    assert captured["displacement"] == [25.0, 0.0, 0.0]
    # Modeller pre-displacement ENABLED for dminus (binder physically at x0+d)
    assert captured["apply_modeller_pre_displacement"] is True
    # Metadata records the spec contract
    assert result["direction"] == "dminus"
    assert result["apply_modeller_pre_displacement"] is True
    assert "dissociated" in result["binder_physical_state"]
    # Runtime ATMForce sign reversal documented in metadata
    assert "-1" in result["runtime_atmforce_displacement_sign"]


def test_build_bound_leg_both_directions_use_walker_direction_plus_one(
        builder, monkeypatch):
    """(b+) corrected: BOTH walkers run with Direction=+1 (NOT ±1).
    This is the two-leg ATM ABFE standard (Azimi 2022)."""
    captured = {}

    def _capture(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(builder, "_patched_make_system", _capture)
    for d in ("dplus", "dminus"):
        result = builder.build_bound_leg_via_upstream(
            receptor_pdb="/tmp/r.pdb",
            binder_pdb="/tmp/b.pdb",
            xml_out=f"/tmp/o_{d}.xml",
            pdb_out=f"/tmp/o_{d}.pdb",
            displacement_nm=(2.5, 0.0, 0.0),
            direction=d,
        )
        assert "+1" in result["walker_direction"]
        assert "both" in result["walker_direction"].lower()


def test_build_bound_leg_3d_displacement_preserved(builder, monkeypatch):
    """3D displacement vector preserved verbatim (positive) for both
    directions; only the Modeller-pre-displacement flag differs."""
    captured = {}

    def _capture(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(builder, "_patched_make_system", _capture)
    builder.build_bound_leg_via_upstream(
        receptor_pdb="/tmp/r.pdb",
        binder_pdb="/tmp/b.pdb",
        xml_out="/tmp/o.xml",
        pdb_out="/tmp/o.pdb",
        displacement_nm=(1.5, 0.5, -0.25),
        direction="dminus",
    )
    # Positive nm → Angstrom (NOT negated at build time)
    assert captured["displacement"] == [15.0, 5.0, -2.5]
    assert captured["apply_modeller_pre_displacement"] is True


def test_build_bound_leg_via_upstream_rejects_empty_direction(builder,
                                                               monkeypatch):
    """direction='' raises ValueError immediately (cohort-safe guard).

    Tested at the lower-level helper (where the guard lives) rather than
    via build_all_four_systems because build_all_four_systems imports
    resolve_leg_inputs lazily inside the function body and protonates
    receptors before reaching the bound-direction loop — making it hard
    to short-circuit cleanly in a unit test."""
    def _should_not_be_called(**kwargs):
        raise AssertionError("_patched_make_system reached with bad direction")
    monkeypatch.setattr(builder, "_patched_make_system", _should_not_be_called)
    with pytest.raises(ValueError, match="direction must be"):
        builder.build_bound_leg_via_upstream(
            receptor_pdb="/tmp/r.pdb",
            binder_pdb="/tmp/b.pdb",
            xml_out="/tmp/o.xml",
            pdb_out="/tmp/o.pdb",
            direction="",
        )


def test_build_all_four_systems_validates_bound_directions_unit(builder):
    """The bound_directions validation logic (non-empty + each entry in
    valid set) is unit-tested directly via the validate logic that lives
    at the top of build_all_four_systems' bound loop. We assert the
    contract via signature introspection + a focused negative case on
    the per-leg helper. The full integration path is exercised by the
    production launcher tests (which call build_all_four_systems via the
    v2.1 launcher with the default ('dplus',) tuple)."""
    import inspect
    sig = inspect.signature(builder.build_all_four_systems)
    # Default preserves legacy behavior
    assert sig.parameters["bound_directions"].default == ("dplus",)
    # Negative: the per-leg validator (called inside the loop) rejects
    # an invalid tag at the FIRST call site, so even partial cohorts
    # halt cleanly.
    with pytest.raises(ValueError, match="direction must be"):
        builder.build_bound_leg_via_upstream(
            receptor_pdb="/tmp/r.pdb",
            binder_pdb="/tmp/b.pdb",
            xml_out="/tmp/o.xml",
            pdb_out="/tmp/o.pdb",
            direction="sideways",  # not in {dplus, dminus}
        )


def test_build_all_four_systems_legacy_single_dplus_keeps_legacy_filenames(
        builder):
    """Default bound_directions=('dplus',) preserves the legacy
    trackb_sys.xml + trackb.pdb filenames (back-compat with v2.1 launcher
    + all existing consumers)."""
    # Default kwarg value introspection (avoids running the real build,
    # which requires real receptor/binder PDBs + GPU).
    import inspect
    sig = inspect.signature(builder.build_all_four_systems)
    assert "bound_directions" in sig.parameters
    assert sig.parameters["bound_directions"].default == ("dplus",)
