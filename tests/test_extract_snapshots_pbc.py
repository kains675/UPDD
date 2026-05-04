"""Tests for utils/extract_snapshots.py L387 PBC defect patch (v54)
and CONECT serial-number wraparound fix (v55).

Phase β CATR-ESCALATE finding (2026-04-28):
    Pre-v54 ``save_snapshots()`` invoked only ``_unwrap_binder_manual()``
    (centroid-shift), which did not repair within-chain peptide bonds when
    the binder chain straddled a periodic boundary at the HETATM/ATOM tag
    junction (e.g. VAL5-C ↔ MTR6-N for 7TL8_MTR6). Result was broken
    peptide bond (~113 Å) in saved PDBs.

v54 fix: layer ``_unwrap_frame_pbc()`` (which calls
``mdtraj.Trajectory.image_molecules(anchor_molecules=..., make_whole=True)``)
BEFORE ``_unwrap_binder_manual()`` so within-chain bonds are made whole
prior to the centroid shift.

Phase β #84 CATR finding (2026-04-29):
    Post-v54 large-system snapshots (>99,999 atoms — 2QKI_Cp4 ~123k,
    7TL8_MTR6_s19/s42 ~144k) failed MM-PBSA at OpenMM's
    ``ForceField.createSystem`` with "No template found for residue ...
    (MTR)" because mdtraj save_pdb wraps the 5-digit serial column at
    99,999 → wrapped water atoms collided with ncAA atoms in
    ``_inject_ncaa_conect`` 's serial-keyed lookup, producing 0 CONECT
    records.

v55 fix: ``_inject_ncaa_conect`` switches to atom-index based
disambiguation (Part A) and re-serializes the entire PDB with
hex-extended serials when total atoms > 99,999 (Part B), so OpenMM's
``atomByNumber`` map is collision-free.

Test scope:
    T1 — 7TL8_MTR6_s7 regression: known BROKEN trajectory must produce
         INTACT (< 5 Å) peptide bond after patch.
    T2 — 7TL8_MTR6_s19 cross-validation: known INTACT trajectory must
         remain INTACT (no regression).
    T3 — 1EBP_MTR13_s7 Phase α regression: verified-INTACT must remain.
    T4 — Synthetic edge case: 2-residue binder straddling periodic
         boundary on a small box.
    T5 — Layered unwrap order: ``_unwrap_frame_pbc`` called BEFORE
         ``_unwrap_binder_manual`` per frame (and both called once).
    T6 — ncAA CONECT injection preserved: CONECT records present in
         saved PDB after patched ``save_snapshots``.
    T7 — Small system regression (1EBP_MTR13 ~62k atoms): v55 patch
         preserves CONECT count and byte-content for under-wraparound
         systems.
    T8 — Large system Cp4 (~123k atoms): pre-v55 produced 0 ncAA CONECT
         records (only 4 mdtraj-default peptide CONECTs); post-v55 must
         produce ≥10 valid CONECT records via re-serialization.
    T9 — Largest system 7TL8_MTR6_s19 (~144k atoms): same as T8 for
         144k-atom system; cross-validate with already-INTACT s7.
    T10 — Boundary edge case at exactly 99,999 / 100,001 atoms.
    T11 — Round-trip via OpenMM PDBFile: patched PDB MTR residue must
          have expected internal bonds parseable back.
    T12 — ForceField template matching: real Cp4 + 7TL8_MTR6_s19
          patched snapshots must pass ``createSystem`` (production
          smoke test for the bug).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

# Ensure utils/ on path (mirrors tests/conftest.py).
_REPO_ROOT = Path(__file__).resolve().parents[1]
_UTILS_DIR = _REPO_ROOT / "utils"
if str(_UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(_UTILS_DIR))


def _restore_openmm_for_mdtraj():
    """Restore real ``openmm`` modules if a prior test (e.g.
    ``test_cofactor_injection.py``) stubbed them — those stubs lack
    ``openmm.unit.Quantity`` which mdtraj requires when loading DCDs.

    Pre-existing test isolation issue (unrelated to this patch); we
    self-heal here so this test file is robust to test ordering.
    """
    omm = sys.modules.get("openmm")
    omm_unit = sys.modules.get("openmm.unit")
    if omm_unit is not None and not hasattr(omm_unit, "Quantity"):
        # Stubbed — drop and re-import the real packages.
        for k in list(sys.modules):
            if k == "openmm" or k.startswith("openmm."):
                del sys.modules[k]
        try:
            import openmm  # noqa: F401
            import openmm.unit  # noqa: F401
            import openmm.app  # noqa: F401
        except ImportError:
            pass


_restore_openmm_for_mdtraj()
import mdtraj as md  # noqa: E402

from extract_snapshots import save_snapshots  # noqa: E402
import extract_snapshots as es  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

OUTPUTS_ROOT = _REPO_ROOT / "outputs"

# Classification threshold (Å). Mirrors verify_pbc_integrity.INTACT_MAX.
INTACT_MAX_ANGSTROM = 5.0

# Extraction parameters (mirror production: K-Means, n_snapshots=25).
N_SNAPSHOTS = 25


def _system_paths(system_dir: Path):
    """Return (dcd_path, top_pdb) for a calib system, or (None, None)."""
    md_root = system_dir / "mdresult"
    if not md_root.is_dir():
        return None, None
    dcds = sorted(md_root.glob("*_restrained.dcd"))
    tops = sorted(md_root.glob("*_final.pdb"))
    if not dcds or not tops:
        return None, None
    return dcds[0], tops[0]


def _load_and_strip(dcd_path: Path, top_pdb: Path):
    """Load DCD with topology and strip solvent (matches production path)."""
    _restore_openmm_for_mdtraj()
    traj = md.load(str(dcd_path), top=str(top_pdb))
    # Use the same selection rule production uses (no solvent / counter-ions).
    traj_protein = es.strip_solvent(traj)
    return traj_protein


def _select_25_frames(n_frames: int):
    """Pick 25 frames spaced as a deterministic stand-in for K-Means.

    Production uses K-Means; for the regression test we need a deterministic
    selection that reliably touches the BROKEN frames in 7TL8_MTR6_s7. An
    even spread covers the trajectory more uniformly and matches the
    *number* of snapshots the production pipeline writes.
    """
    if n_frames <= N_SNAPSHOTS:
        return list(range(n_frames))
    step = max(1, n_frames // N_SNAPSHOTS)
    selected = [i * step for i in range(N_SNAPSHOTS)]
    selected = [min(f, n_frames - 1) for f in selected]
    return sorted(set(selected))


def _find_chain_b(topology):
    """Return the chain whose chain ID is 'B' (binder), or None."""
    for ch in topology.chains:
        cid = getattr(ch, "chain_id", None)
        if cid == "B":
            return ch
    # Fallback: second chain
    chains = list(topology.chains)
    return chains[1] if len(chains) >= 2 else None


def _find_ncaa_in_chain(chain):
    """Return the first non-standard residue (e.g. MTR) in chain."""
    standard = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX",
    }
    solvent = {"HOH", "WAT", "NA", "CL", "K", "MG", "CA", "ZN"}
    for r in chain.residues:
        nm = r.name.upper()
        if nm in standard or nm in solvent:
            continue
        return r
    return None


def _measure_pdb_peptide_distances(pdb_paths):
    """For each PDB, locate (binder ncAA prev_C, ncAA N) atoms and return
    the list of distances in Å.

    Skips a PDB if either atom is absent (e.g. ncAA at N-term).
    """
    _restore_openmm_for_mdtraj()
    distances = []
    for p in pdb_paths:
        try:
            t = md.load(str(p))
        except Exception:
            continue
        chain = _find_chain_b(t.topology)
        if chain is None:
            continue
        ncaa = _find_ncaa_in_chain(chain)
        if ncaa is None:
            continue
        # Find prev residue in same chain.
        prev = None
        for r in chain.residues:
            if r.index == ncaa.index:
                break
            prev = r
        if prev is None:
            continue
        prev_C = next((a for a in prev.atoms if a.name == "C"), None)
        ncaa_N = next((a for a in ncaa.atoms if a.name == "N"), None)
        if prev_C is None or ncaa_N is None:
            continue
        d_nm = np.linalg.norm(
            t.xyz[0, prev_C.index, :] - t.xyz[0, ncaa_N.index, :]
        )
        distances.append(float(d_nm * 10.0))
    return distances


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def s7_7tl8_traj():
    sysdir = OUTPUTS_ROOT / "7TL8_MTR6_calib_s7"
    dcd, top = _system_paths(sysdir)
    if dcd is None:
        pytest.skip("7TL8_MTR6_s7 trajectory not present")
    traj = _load_and_strip(dcd, top)
    return traj


@pytest.fixture(scope="module")
def s19_7tl8_traj():
    sysdir = OUTPUTS_ROOT / "7TL8_MTR6_calib_s19"
    dcd, top = _system_paths(sysdir)
    if dcd is None:
        pytest.skip("7TL8_MTR6_s19 trajectory not present")
    traj = _load_and_strip(dcd, top)
    return traj


@pytest.fixture(scope="module")
def s7_1ebp_traj():
    sysdir = OUTPUTS_ROOT / "1EBP_MTR13_calib_s7"
    dcd, top = _system_paths(sysdir)
    if dcd is None:
        pytest.skip("1EBP_MTR13_s7 trajectory not present")
    traj = _load_and_strip(dcd, top)
    return traj


# ---------------------------------------------------------------------------
# T1 — 7TL8_MTR6_s7 BROKEN→INTACT regression
# ---------------------------------------------------------------------------

def test_t1_7tl8_s7_broken_to_intact(s7_7tl8_traj, tmp_path):
    """Patched save_snapshots() must produce all-INTACT peptide bonds for the
    canonical 'smoking gun' trajectory."""
    selected = _select_25_frames(s7_7tl8_traj.n_frames)
    saved = save_snapshots(
        s7_7tl8_traj, selected, str(tmp_path), "7TL8_MTR6"
    )
    assert len(saved) == len(selected), "save count mismatch"
    distances = _measure_pdb_peptide_distances(saved)
    assert distances, "no measurable peptide-bond distances (test setup bug)"

    median = float(np.median(distances))
    max_d = float(np.max(distances))
    n_intact = sum(1 for d in distances if d < INTACT_MAX_ANGSTROM)
    print(
        f"\n[T1 7TL8_MTR6_s7] frames={len(distances)} median={median:.2f} Å "
        f"max={max_d:.2f} Å INTACT={n_intact}/{len(distances)}"
    )
    # Hard assertion: every frame must be < 5 Å.
    assert max_d < INTACT_MAX_ANGSTROM, (
        f"BROKEN frames remain post-patch: max={max_d:.2f} Å "
        f"(median {median:.2f}, expected < 5 Å). v54 patch ineffective."
    )
    assert median < INTACT_MAX_ANGSTROM


# ---------------------------------------------------------------------------
# T2 — 7TL8_MTR6_s19 INTACT preservation
# ---------------------------------------------------------------------------

def test_t2_7tl8_s19_intact_preserved(s19_7tl8_traj, tmp_path):
    """Already-INTACT trajectory must remain INTACT (no regression)."""
    selected = _select_25_frames(s19_7tl8_traj.n_frames)
    saved = save_snapshots(
        s19_7tl8_traj, selected, str(tmp_path), "7TL8_MTR6"
    )
    distances = _measure_pdb_peptide_distances(saved)
    assert distances, "no measurable peptide-bond distances"

    median = float(np.median(distances))
    max_d = float(np.max(distances))
    n_intact = sum(1 for d in distances if d < INTACT_MAX_ANGSTROM)
    print(
        f"\n[T2 7TL8_MTR6_s19] frames={len(distances)} median={median:.2f} Å "
        f"max={max_d:.2f} Å INTACT={n_intact}/{len(distances)}"
    )
    assert max_d < INTACT_MAX_ANGSTROM
    assert median < 2.0, f"unexpected median for INTACT system: {median:.2f} Å"


# ---------------------------------------------------------------------------
# T3 — 1EBP_MTR13_s7 Phase α preservation
# ---------------------------------------------------------------------------

def test_t3_1ebp_s7_phase_alpha_preserved(s7_1ebp_traj, tmp_path):
    """Phase α verdict §3.3 VALID system must stay INTACT."""
    selected = _select_25_frames(s7_1ebp_traj.n_frames)
    saved = save_snapshots(
        s7_1ebp_traj, selected, str(tmp_path), "1EBP_MTR13"
    )
    distances = _measure_pdb_peptide_distances(saved)
    assert distances, "no measurable peptide-bond distances"

    median = float(np.median(distances))
    max_d = float(np.max(distances))
    n_intact = sum(1 for d in distances if d < INTACT_MAX_ANGSTROM)
    print(
        f"\n[T3 1EBP_MTR13_s7] frames={len(distances)} median={median:.2f} Å "
        f"max={max_d:.2f} Å INTACT={n_intact}/{len(distances)}"
    )
    assert max_d < INTACT_MAX_ANGSTROM
    assert median < 2.5, f"unexpected median for INTACT system: {median:.2f} Å"


# ---------------------------------------------------------------------------
# T4 — Synthetic edge case at periodic boundary
# ---------------------------------------------------------------------------

def test_t4_synthetic_periodic_boundary(tmp_path):
    """Construct a tiny 2-residue binder spanning a periodic boundary and
    confirm that image_molecules + manual unwrap together produce an intact
    bond after save.

    Box: 2 nm cubic. Residue 1 carbonyl C at x=0.10 nm (near origin),
    residue 2 amide N at x=1.90 nm (opposite face). Direct distance is
    1.80 nm; after wrap they should be ~0.20 nm apart (still in lattice
    sense).
    """
    # Build minimum 2-residue ALA dipeptide via mdtraj's Topology API.
    top = md.Topology()
    chain = top.add_chain()
    # Residue 1 (ALA): N, CA, C, O
    r1 = top.add_residue("ALA", chain)
    n1 = top.add_atom("N", md.element.nitrogen, r1)
    ca1 = top.add_atom("CA", md.element.carbon, r1)
    c1 = top.add_atom("C", md.element.carbon, r1)
    o1 = top.add_atom("O", md.element.oxygen, r1)
    # Residue 2 (ALA): N, CA, C, O
    r2 = top.add_residue("ALA", chain)
    n2 = top.add_atom("N", md.element.nitrogen, r2)
    ca2 = top.add_atom("CA", md.element.carbon, r2)
    c2 = top.add_atom("C", md.element.carbon, r2)
    o2 = top.add_atom("O", md.element.oxygen, r2)
    # Bonds — including the peptide bond C1–N2 across the box boundary.
    top.add_bond(n1, ca1)
    top.add_bond(ca1, c1)
    top.add_bond(c1, o1)
    top.add_bond(c1, n2)         # peptide bond (key for image_molecules)
    top.add_bond(n2, ca2)
    top.add_bond(ca2, c2)
    top.add_bond(c2, o2)

    # Coordinates (in nm). Place residue 1 at x≈0.0-0.2; residue 2 at
    # x≈1.7-1.9 — peptide bond C1↔N2 spans 1.8 nm, but on a 2 nm box the
    # minimum-image distance is 0.2 nm.
    xyz = np.array([[
        [0.05, 1.0, 1.0],   # N1
        [0.10, 1.0, 1.0],   # CA1
        [0.15, 1.0, 1.0],   # C1
        [0.20, 1.0, 1.0],   # O1
        [1.95, 1.0, 1.0],   # N2  (wrapped to far face)
        [1.90, 1.0, 1.0],   # CA2
        [1.85, 1.0, 1.0],   # C2
        [1.80, 1.0, 1.0],   # O2
    ]], dtype=np.float32)
    box = np.array([[2.0, 2.0, 2.0]], dtype=np.float32)
    angles = np.array([[90.0, 90.0, 90.0]], dtype=np.float32)
    traj = md.Trajectory(xyz, top, time=[0.0],
                         unitcell_lengths=box, unitcell_angles=angles)

    # Pre-patch the CONECT injector to a no-op (no ncAA registry match so it's
    # already a no-op, but skip lazy import noise on edge-case test).
    with patch.object(es, "_inject_ncaa_conect", lambda _p: None):
        saved = save_snapshots(traj, [0], str(tmp_path), "synthetic")
    assert len(saved) == 1
    out = md.load(str(saved[0]))
    # Peptide bond C1 ↔ N2: atoms at indices 2 and 4 in original ordering;
    # mdtraj.Trajectory.save_pdb may renumber, so locate by name+resseq.
    c1_atom = next(a for a in out.topology.atoms
                   if a.name == "C" and a.residue.resSeq == r1.resSeq)
    n2_atom = next(a for a in out.topology.atoms
                   if a.name == "N" and a.residue.resSeq == r2.resSeq)
    d_nm = float(np.linalg.norm(
        out.xyz[0, c1_atom.index, :] - out.xyz[0, n2_atom.index, :]
    ))
    d_ang = d_nm * 10.0
    print(
        f"\n[T4 synthetic] box=2nm pre-distance=18.0 Å post-distance={d_ang:.2f} Å"
    )
    # Post-image_molecules the C1↔N2 distance should be the minimum-image
    # value (~2 Å), not the 18 Å boundary-cross distance.
    assert d_ang < INTACT_MAX_ANGSTROM, (
        f"image_molecules failed to bridge periodic boundary: {d_ang:.2f} Å"
    )


# ---------------------------------------------------------------------------
# T5 — Layered unwrap order
# ---------------------------------------------------------------------------

def test_t5_unwrap_call_order(tmp_path):
    """Verify save_snapshots calls _unwrap_frame_pbc BEFORE
    _unwrap_binder_manual, and exactly once each per frame.

    Implementation note: ``_unwrap_frame_pbc`` itself calls
    ``_unwrap_binder_manual`` as an internal fallback — so the *outer*
    (top-level) call sequence we're verifying is:
    ``save_snapshots`` → ``_unwrap_frame_pbc`` → ``_unwrap_binder_manual``
    in that order. We patch only the names referenced inside
    ``save_snapshots`` (module-level lookups) by tracking the order of
    invocation.
    """
    # Construct a minimal 2-residue trajectory (reuses T4 builder logic).
    top = md.Topology()
    chain = top.add_chain()
    r1 = top.add_residue("ALA", chain)
    n1 = top.add_atom("N", md.element.nitrogen, r1)
    ca1 = top.add_atom("CA", md.element.carbon, r1)
    c1 = top.add_atom("C", md.element.carbon, r1)
    o1 = top.add_atom("O", md.element.oxygen, r1)
    r2 = top.add_residue("ALA", chain)
    n2 = top.add_atom("N", md.element.nitrogen, r2)
    ca2 = top.add_atom("CA", md.element.carbon, r2)
    c2 = top.add_atom("C", md.element.carbon, r2)
    o2 = top.add_atom("O", md.element.oxygen, r2)
    for a, b in [(n1, ca1), (ca1, c1), (c1, o1), (c1, n2), (n2, ca2),
                 (ca2, c2), (c2, o2)]:
        top.add_bond(a, b)
    xyz = np.zeros((2, 8, 3), dtype=np.float32)
    xyz[:, 1, 0] = 0.10
    xyz[:, 4, 0] = 0.20
    box = np.array([[2.0, 2.0, 2.0]] * 2, dtype=np.float32)
    angles = np.array([[90.0, 90.0, 90.0]] * 2, dtype=np.float32)
    traj = md.Trajectory(xyz, top, time=[0.0, 1.0],
                         unitcell_lengths=box, unitcell_angles=angles)

    call_log = []

    def _wrap_pbc(frame, *args, **kwargs):
        call_log.append("_unwrap_frame_pbc")
        return frame

    def _wrap_manual(frame, *args, **kwargs):
        call_log.append("_unwrap_binder_manual")

    with patch.object(es, "_unwrap_frame_pbc", _wrap_pbc), \
         patch.object(es, "_unwrap_binder_manual", _wrap_manual), \
         patch.object(es, "_inject_ncaa_conect", lambda _p: None):
        save_snapshots(traj, [0, 1], str(tmp_path), "ordertest")

    # Two frames, two outer calls each: pbc then manual.
    assert call_log == [
        "_unwrap_frame_pbc", "_unwrap_binder_manual",
        "_unwrap_frame_pbc", "_unwrap_binder_manual",
    ], f"unexpected outer call order: {call_log}"


# ---------------------------------------------------------------------------
# T6 — ncAA CONECT injection preserved
# ---------------------------------------------------------------------------

def test_t6_ncaa_conect_preserved(s7_7tl8_traj, tmp_path):
    """First saved 7TL8 PDB must contain CONECT records (peptide + sidechain
    bonds for MTR), confirming _inject_ncaa_conect is still chained."""
    saved = save_snapshots(s7_7tl8_traj, [0], str(tmp_path), "7TL8_MTR6")
    assert saved
    text = Path(saved[0]).read_text()
    # CONECT lines emitted post-save by _inject_ncaa_conect for ncAA.
    n_conect = sum(
        1 for line in text.splitlines() if line.startswith("CONECT")
    )
    print(f"\n[T6 ncAA CONECT] count={n_conect} in {Path(saved[0]).name}")
    assert n_conect > 0, "no CONECT records emitted (ncAA injector regression)"


# ===========================================================================
# v55 CONECT serial-number wraparound fix
# ===========================================================================
#
# These tests exercise ``_inject_ncaa_conect`` directly on existing
# postl387-patch snapshot PDBs (READ-ONLY copy to tmp_path). They do NOT
# regenerate snapshots — that is a separate dispatch.
# ---------------------------------------------------------------------------

import shutil  # noqa: E402

# Path to a Cp4 snapshot PDB known to have >99,999 atoms (the bug case).
CP4_SNAP_PDB = (
    OUTPUTS_ROOT
    / "2QKI_Cp4_calib_s7"
    / "snapshots_n25_postl387_patch"
    / "2QKI_Cp4_restrained_trimmed_snap01_f0.pdb"
)
TL8_S19_SNAP_GLOB = (
    OUTPUTS_ROOT / "7TL8_MTR6_calib_s19" / "snapshots_n25_postl387_patch"
)
EBP_SNAP_PDB = (
    OUTPUTS_ROOT
    / "1EBP_MTR13_calib_s7"
    / "snapshots_n25_postl387_patch"
    / "1EBP_MTR13_snap01_f2.pdb"
)


def _count_conect(pdb_path):
    n = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("CONECT"):
                n += 1
    return n


def _count_atoms(pdb_path):
    n = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                n += 1
    return n


# ---------------------------------------------------------------------------
# T7 — Small-system regression (1EBP_MTR13, ~62k atoms, under wraparound)
# ---------------------------------------------------------------------------

def test_t7_small_system_regression(tmp_path):
    """1EBP_MTR13_s7 snapshot has ~62k atoms (under 99,999), so v55 patch
    must be a no-op for re-serialization. CONECT count must equal pre-v55
    (no regression)."""
    if not EBP_SNAP_PDB.is_file():
        pytest.skip(f"1EBP_MTR13_s7 sample PDB not found: {EBP_SNAP_PDB}")
    dst = tmp_path / "1EBP_v55.pdb"
    shutil.copy(str(EBP_SNAP_PDB), str(dst))
    pre_conect = _count_conect(dst)
    pre_atoms = _count_atoms(dst)

    es._inject_ncaa_conect(str(dst))

    post_conect = _count_conect(dst)
    post_atoms = _count_atoms(dst)
    print(
        f"\n[T7 1EBP_MTR13_s7] atoms={pre_atoms} pre_conect={pre_conect} "
        f"post_conect={post_conect}"
    )

    # The source file already had CONECT injected (it's a postl387 snapshot
    # which was processed by an earlier _inject_ncaa_conect run). Re-running
    # on the already-injected file should be idempotent: CONECT count must
    # not decrease, and we expect ≥ 10 (MTR + flanking).
    assert post_conect >= 10, (
        f"unexpected CONECT count {post_conect} (expected ≥ 10) — "
        f"v55 small-system regression"
    )
    assert post_conect == pre_conect, (
        f"CONECT count changed on re-injection: {pre_conect} → {post_conect} "
        f"(expected idempotent for under-wraparound systems)"
    )
    assert post_atoms == pre_atoms, "atom count must be preserved"


# ---------------------------------------------------------------------------
# T8 — Cp4 (~123k atoms, >99,999 — was failing pre-v55)
# ---------------------------------------------------------------------------

def test_t8_cp4_large_system(tmp_path):
    """2QKI_Cp4_s7 (~123k atoms) had wrapped serials in source PDB →
    pre-v55 _inject_ncaa_conect emitted 0 ncAA CONECT records (only 4
    mdtraj-defaults). Post-v55 must produce ≥ 10 (MTR + flanking)."""
    if not CP4_SNAP_PDB.is_file():
        pytest.skip(f"2QKI_Cp4_s7 sample PDB not found: {CP4_SNAP_PDB}")
    dst = tmp_path / "Cp4_v55.pdb"
    shutil.copy(str(CP4_SNAP_PDB), str(dst))
    pre_conect = _count_conect(dst)
    pre_atoms = _count_atoms(dst)

    assert pre_atoms > 99999, f"test fixture atom count: {pre_atoms}"

    es._inject_ncaa_conect(str(dst))

    post_conect = _count_conect(dst)
    post_atoms = _count_atoms(dst)
    print(
        f"\n[T8 2QKI_Cp4_s7] atoms={pre_atoms} pre_conect={pre_conect} "
        f"post_conect={post_conect} (Δ=+{post_conect - pre_conect})"
    )

    # Pre-v55 fixture state: 4 default mdtraj-emitted peptide CONECTs only.
    # Post-v55: should add ~28 MTR-related CONECTs.
    assert post_conect >= 10, (
        f"v55 patch did not emit CONECT records for MTR in Cp4 "
        f"(post={post_conect}, expected ≥ 10)"
    )
    # Atom count preserved (re-serialization renumbers but does not drop atoms).
    assert post_atoms == pre_atoms, "atom count must be preserved"


# ---------------------------------------------------------------------------
# T9 — 7TL8_MTR6_s19 (~144k atoms, largest)
# ---------------------------------------------------------------------------

def test_t9_7tl8_s19_largest(tmp_path):
    """7TL8_MTR6_s19 (~144k atoms) is the largest of the affected systems.
    Same expectations as T8."""
    if not TL8_S19_SNAP_GLOB.is_dir():
        pytest.skip(f"7TL8_MTR6_s19 dir not found: {TL8_S19_SNAP_GLOB}")
    pdbs = sorted(TL8_S19_SNAP_GLOB.glob("*.pdb"))
    if not pdbs:
        pytest.skip("no PDBs in 7TL8_MTR6_s19 snapshots dir")

    dst = tmp_path / "7TL8_s19_v55.pdb"
    shutil.copy(str(pdbs[0]), str(dst))
    pre_conect = _count_conect(dst)
    pre_atoms = _count_atoms(dst)

    assert pre_atoms > 99999, f"test fixture atom count: {pre_atoms}"

    es._inject_ncaa_conect(str(dst))

    post_conect = _count_conect(dst)
    post_atoms = _count_atoms(dst)
    print(
        f"\n[T9 7TL8_MTR6_s19] atoms={pre_atoms} pre_conect={pre_conect} "
        f"post_conect={post_conect} (Δ=+{post_conect - pre_conect})"
    )
    assert post_conect >= 10, (
        f"v55 patch did not emit CONECT records for MTR in 7TL8_s19 "
        f"(post={post_conect}, expected ≥ 10)"
    )
    assert post_atoms == pre_atoms


# ---------------------------------------------------------------------------
# T10 — Synthetic boundary edge case at 99,999 / 100,001 atoms
# ---------------------------------------------------------------------------

def _pdb_atom_line(serial_str, atom_name, resname, chain, resnum, elem):
    """Build a properly-formatted PDB ATOM line per the format spec.

    Cols (1-based):
        1-6   record name 'ATOM  '
        7-11  serial (5)
        13-16 atom name (4) — right-padded for 1-letter element, else padded
        17    alt loc
        18-20 resname (3)
        22    chain
        23-26 resseq (4)
        31-38 x  39-46 y  47-54 z (each 8.3f)
        55-60 occupancy (6.2f)
        61-66 tempfactor (6.2f)
        77-78 element (2)
    """
    # Atom-name placement: 1-letter elements occupy col 14, name in 14-16.
    # Pre-pend space when the first char is the element symbol.
    if len(atom_name) >= 4:
        name_field = atom_name[:4]
    else:
        # 1-3 char names: right-pad after a leading space (typical for H/N/C/O).
        name_field = (" " + atom_name).ljust(4)
    return (
        "ATOM  "                                # 1-6
        + f"{serial_str:>5}"                    # 7-11
        + " "                                   # 12
        + name_field                            # 13-16
        + " "                                   # 17 alt loc
        + f"{resname:>3}"                       # 18-20
        + " "                                   # 21
        + f"{chain:1}"                          # 22
        + f"{resnum:>4d}"                       # 23-26
        + "    "                                # 27-30
        + f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}"     # 31-54
        + f"{1.00:6.2f}{0.00:6.2f}"             # 55-66
        + "          "                          # 67-76
        + f"{elem:>2}"                          # 77-78
        + "\n"
    )


def _build_synthetic_pdb(n_atoms: int, ncaa_at_index: int, out_path: Path):
    """Construct a synthetic PDB with n_atoms total, with an MTR ncAA at
    ``ncaa_at_index``-th atom position (1-based). Filler atoms are HOH
    waters. Coordinates are zero (geometry irrelevant for CONECT injection).

    Serial numbers wrap at 99,999 → 0, mimicking mdtraj save_pdb behavior.
    """
    mtr_atom_names = [
        "N", "H", "CA", "HA", "C", "O", "CB", "HB2", "HB3",
        "CG", "CD1", "HD1", "CD2", "NE1", "CM", "HM1", "HM2", "HM3",
        "CE2", "CE3", "HE3", "CZ2", "HZ2", "CZ3", "HZ3", "CH2", "HH2",
    ]  # 27 atoms (matches MTR sample seen in production)
    n_mtr = len(mtr_atom_names)
    if ncaa_at_index < 1 or ncaa_at_index + n_mtr - 1 > n_atoms:
        raise ValueError(
            f"ncaa_at_index={ncaa_at_index} won't fit {n_mtr} atoms in {n_atoms}"
        )

    def _ser_str(i: int) -> str:
        # mdtraj wraps at 99999 → 0 → 1 → ...
        if i <= 99999:
            return f"{i:5d}"
        return f"{(i - 100000) % 100000:5d}"

    def _water_atom(local_i, water_resnum):
        """Return (atom_name, element) for a HOH atom at sub-index 0/1/2."""
        return [("O", "O"), ("H1", "H"), ("H2", "H")][local_i % 3]

    lines = []
    resnum_water = 0
    # Pre-MTR waters (atoms 1 .. ncaa_at_index - 1)
    for i in range(1, ncaa_at_index):
        local = (i - 1) % 3
        if local == 0:
            resnum_water += 1
        name, elem = _water_atom(local, resnum_water)
        lines.append(_pdb_atom_line(
            _ser_str(i), name, "HOH", "A",
            (resnum_water - 1) % 9999 + 1, elem,
        ))
    # MTR atoms — chain B, resnum 1
    for j, atom_name in enumerate(mtr_atom_names):
        i = ncaa_at_index + j
        if atom_name[0] == "H":
            elem = "H"
        elif atom_name[0] == "N":
            elem = "N"
        elif atom_name[0] == "O":
            elem = "O"
        else:
            elem = "C"
        lines.append(_pdb_atom_line(
            _ser_str(i), atom_name, "MTR", "B", 1, elem,
        ))
    # Post-MTR waters
    for i in range(ncaa_at_index + n_mtr, n_atoms + 1):
        local = (i - ncaa_at_index - n_mtr) % 3
        if local == 0:
            resnum_water += 1
        name, elem = _water_atom(local, resnum_water)
        lines.append(_pdb_atom_line(
            _ser_str(i), name, "HOH", "A",
            (resnum_water - 1) % 9999 + 1, elem,
        ))
    lines.append("END\n")
    out_path.write_text("".join(lines))


def test_t10_boundary_99999(tmp_path):
    """At exactly 99,999 atoms (boundary inclusive): no wraparound. Patch
    must produce ≥ 10 CONECT for MTR.

    At 100,001 atoms (just over): wraparound triggered. Patch must still
    produce ≥ 10 CONECT for MTR.
    """
    for label, n_atoms in [("99999", 99999), ("100001", 100001)]:
        pdb = tmp_path / f"synthetic_{label}.pdb"
        # Place MTR early (at atom 100) — well before any wrap
        _build_synthetic_pdb(n_atoms, ncaa_at_index=100, out_path=pdb)
        n_atoms_actual = _count_atoms(pdb)
        es._inject_ncaa_conect(str(pdb))
        n_conect = _count_conect(pdb)
        n_atoms_post = _count_atoms(pdb)
        print(
            f"\n[T10 boundary={label}] atoms={n_atoms_actual} "
            f"post_conect={n_conect}"
        )
        assert n_atoms_post == n_atoms_actual, (
            f"boundary={label}: atom count changed {n_atoms_actual} → {n_atoms_post}"
        )
        assert n_conect >= 10, (
            f"boundary={label}: only {n_conect} CONECT records emitted "
            f"(expected ≥ 10 for MTR)"
        )


# ---------------------------------------------------------------------------
# T11 — Round-trip via OpenMM PDBFile (replaces mdtraj round-trip — mdtraj
# does not parse hex-extended serials but OpenMM does; the pipeline
# downstream uses OpenMM, so OpenMM round-trip is the relevant one).
# ---------------------------------------------------------------------------

def test_t11_openmm_pdbfile_roundtrip(tmp_path):
    """Patched Cp4 snapshot must be readable by OpenMM PDBFile and the MTR
    residue must carry expected internal bonds.

    Expected MTR internal bonds (from amber14 + extension): ≥ 25
    (parent TRP heavy + parent H + extension CM + 3 methyl Hs).
    """
    if not CP4_SNAP_PDB.is_file():
        pytest.skip(f"2QKI_Cp4_s7 sample PDB not found: {CP4_SNAP_PDB}")
    _restore_openmm_for_mdtraj()
    from openmm.app import PDBFile  # noqa: WPS433

    dst = tmp_path / "Cp4_v55_t11.pdb"
    shutil.copy(str(CP4_SNAP_PDB), str(dst))
    es._inject_ncaa_conect(str(dst))

    pdb = PDBFile(str(dst))
    mtr_atoms = set()
    for r in pdb.topology.residues():
        if r.name == "MTR":
            for a in r.atoms():
                mtr_atoms.add(a.index)

    assert mtr_atoms, "MTR residue not found in patched PDB"

    n_internal = 0
    n_external = 0
    for b in pdb.topology.bonds():
        a1, a2 = b[0], b[1]
        in1 = a1.index in mtr_atoms
        in2 = a2.index in mtr_atoms
        if in1 and in2:
            n_internal += 1
        elif in1 or in2:
            n_external += 1

    print(
        f"\n[T11 OpenMM round-trip] MTR atoms={len(mtr_atoms)} "
        f"internal={n_internal} external={n_external}"
    )
    # MTR has 27 atoms post-strip; expect ≥ 25 internal bonds (covering all
    # heavy + H bonds) and ≥ 1 external (peptide bond to neighbor).
    assert n_internal >= 25, (
        f"insufficient MTR internal bonds: {n_internal} (expected ≥ 25)"
    )
    assert n_external >= 1, (
        f"missing MTR external (peptide) bond: {n_external}"
    )


# ---------------------------------------------------------------------------
# T12 — ForceField template matching (production smoke test)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "label,snap_pdb,params_dir",
    [
        (
            "Cp4",
            CP4_SNAP_PDB,
            OUTPUTS_ROOT / "2QKI_Cp4_calib_s7" / "params",
        ),
        (
            "7TL8_s19",
            None,  # resolved at runtime via glob
            OUTPUTS_ROOT / "7TL8_MTR6_calib_s19" / "params",
        ),
    ],
)
def test_t12_forcefield_create_system(tmp_path, label, snap_pdb, params_dir):
    """Patched snapshot + ncAA XMLs must allow ForceField.createSystem to
    succeed (no "No template found" errors). This is THE production
    smoke test for the wraparound bug.
    """
    if not params_dir.is_dir():
        pytest.skip(f"params dir not found: {params_dir}")
    if snap_pdb is None:
        pdbs = sorted(TL8_S19_SNAP_GLOB.glob("*.pdb"))
        if not pdbs:
            pytest.skip("no PDBs in 7TL8_MTR6_s19 snapshots dir")
        snap_pdb = pdbs[0]
    if not snap_pdb.is_file():
        pytest.skip(f"snapshot not found: {snap_pdb}")

    _restore_openmm_for_mdtraj()
    from openmm.app import (  # noqa: WPS433
        PDBFile, ForceField, NoCutoff, HBonds,
    )

    dst = tmp_path / f"{label}_v55_t12.pdb"
    shutil.copy(str(snap_pdb), str(dst))
    es._inject_ncaa_conect(str(dst))

    xmls = sorted(str(p) for p in params_dir.glob("*.xml"))
    pdb = PDBFile(str(dst))
    ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml", *xmls)

    unmatched, _ = ff.generateTemplatesForUnmatchedResidues(pdb.topology)
    mtr_unmatched = [r for r in unmatched if r.name == "MTR"]
    n_atoms = pdb.topology.getNumAtoms()
    print(
        f"\n[T12 {label}] atoms={n_atoms} MTR_unmatched={len(mtr_unmatched)}"
    )
    assert len(mtr_unmatched) == 0, (
        f"{label}: MTR template not matched after v55 patch — "
        f"{len(mtr_unmatched)} unmatched (was 1 pre-v55)"
    )

    # createSystem should now succeed (no "No template found" for MTR).
    try:
        sys_obj = ff.createSystem(
            pdb.topology,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
            rigidWater=False,
        )
    except Exception as e:
        pytest.fail(
            f"{label}: createSystem failed after v55 patch: "
            f"{type(e).__name__}: {e}"
        )
    assert sys_obj.getNumParticles() == n_atoms, (
        f"{label}: System particle count mismatch"
    )
    print(f"[T12 {label}] createSystem OK: {sys_obj.getNumParticles()} particles")


# ===========================================================================
# v56 — PR-21 Auto-Diagnostic Framework v0.7 (Level 1: Intra-residue Bond Guard)
# ===========================================================================
#
# Trigger event: 7TL8_MTR6_s55 outlier audit (#95, 2026-04-30).
# snap05_f63 had MTR6_N atom imaged across the periodic box (Y = -4.70 Å)
# while H/CA/C/O/sidechain stayed at Y ≈ 109 Å, producing intra-residue
# N↔CA = 113.15 Å while prev_VAL5_C↔MTR6_N stayed at 1.40 Å. The v54
# inter-residue check passes; the intra-residue defect drives ΔG = +686
# kcal/mol single-frame outlier in MM-PBSA.
#
# Test scope:
#     T13 — snap05_f63 forensic detection: load preserved s55 snap05_f63.pdb,
#           run intra_residue_bond_check, expect MTR6 N↔CA flagged broken.
#     T14 — Healthy 1EBP_MTR13_s7 baseline: 0 broken bonds.
#     T15 — Configurable threshold: synthetic ~6 Å bond → BROKEN at thr=5,
#           OK at thr=10.
#     T16 — save_snapshots integration: healthy mini-trajectory, 0 warnings,
#           all snapshots saved.
# ---------------------------------------------------------------------------

# Path to the preserved s55 snap05_f63 forensic PDB.
S55_FORENSIC_PDB = (
    OUTPUTS_ROOT
    / "7TL8_MTR6_calib_s55"
    / "snapshots_n25_postl387_patch_v2"
    / "7TL8_MTR6_snap05_f63.pdb"
)
# Healthy reference snapshot — known INTACT 1EBP_MTR13_s7
EBP_HEALTHY_PDB = (
    OUTPUTS_ROOT
    / "1EBP_MTR13_calib_s7"
    / "snapshots_n25_postl387_patch"
    / "1EBP_MTR13_snap01_f2.pdb"
)


# ---------------------------------------------------------------------------
# T13 — snap05_f63 forensic regression
# ---------------------------------------------------------------------------

def test_t13_snap05_f63_forensic_intra_residue_defect(tmp_path):
    """Load preserved 7TL8_MTR6_s55 snap05_f63 PDB and run
    ``intra_residue_bond_check`` — expect detection of MTR resid 5 (mdtraj
    0-based; resseq 6) N↔CA bond stretched to ~113 Å.

    This is the canonical PR-21 trigger case (#95 audit) and serves as the
    forensic regression: the helper MUST detect this defect.
    """
    if not S55_FORENSIC_PDB.is_file():
        pytest.skip(f"s55 snap05_f63 PDB not present: {S55_FORENSIC_PDB}")
    _restore_openmm_for_mdtraj()
    traj = md.load(str(S55_FORENSIC_PDB))
    broken = es.intra_residue_bond_check(
        traj,
        threshold_angstrom=es.INTRA_RESIDUE_BOND_THRESHOLD_ANGSTROM,
    )
    print(
        f"\n[T13 snap05_f63] intra-residue broken bonds detected: {len(broken)}"
    )
    for rec in broken:
        print(
            f"   frame={rec['frame']} {rec['residue_name']}{rec['residue_index']} "
            f"{rec['atom_a']}↔{rec['atom_b']} = {rec['distance']:.2f} Å"
        )
    assert broken, (
        "Intra-residue bond check FAILED to detect MTR6 N↔CA = 113 Å "
        "defect in snap05_f63 — algorithm bug; PR-21 forensic regression."
    )
    # The smoking-gun defect: MTR (resseq 6 / index 5) N↔CA at ~113 Å +
    # N↔H at ~114 Å. Verify at least one MTR record at > 100 Å.
    mtr_recs = [
        r for r in broken
        if r["residue_name"] == "MTR" and r["distance"] > 100.0
    ]
    assert mtr_recs, (
        f"expected ≥1 MTR intra-residue bond > 100 Å, got {[r for r in broken]}"
    )
    # Specifically the N↔CA bond should be present (canonical 1.46 Å,
    # observed ≈ 113 Å).
    n_ca_recs = [
        r for r in mtr_recs
        if {r["atom_a"], r["atom_b"]} == {"N", "CA"}
    ]
    assert n_ca_recs, (
        f"expected MTR N↔CA bond in broken list; got {mtr_recs}"
    )
    n_ca_dist = n_ca_recs[0]["distance"]
    assert 100.0 < n_ca_dist < 120.0, (
        f"unexpected N↔CA distance: {n_ca_dist:.2f} Å "
        f"(expected ~113 Å per #95 audit)"
    )


# ---------------------------------------------------------------------------
# T14 — Healthy 1EBP_MTR13_s7 baseline
# ---------------------------------------------------------------------------

def test_t14_healthy_1ebp_baseline_no_broken_bonds(tmp_path):
    """Healthy 1EBP_MTR13_s7 snapshot must produce zero broken intra-residue
    bonds (no PBC defect, all covalent distances < 5 Å)."""
    if not EBP_HEALTHY_PDB.is_file():
        pytest.skip(f"1EBP_MTR13_s7 healthy PDB not present: {EBP_HEALTHY_PDB}")
    _restore_openmm_for_mdtraj()
    traj = md.load(str(EBP_HEALTHY_PDB))
    broken = es.intra_residue_bond_check(traj, threshold_angstrom=5.0)
    print(
        f"\n[T14 1EBP_MTR13_s7 healthy] broken bonds: {len(broken)} "
        f"(expected 0)"
    )
    assert broken == [], (
        f"healthy snapshot reported {len(broken)} broken intra-residue "
        f"bonds (expected 0): {broken[:3]}"
    )


# ---------------------------------------------------------------------------
# T15 — Configurable threshold (synthetic 6 Å bond)
# ---------------------------------------------------------------------------

def test_t15_configurable_threshold_synthetic(tmp_path):
    """Construct a synthetic single-residue trajectory with one intra-residue
    bond stretched to 6.0 Å. At ``threshold=5.0`` it must be detected as
    broken; at ``threshold=10.0`` it must NOT be reported.
    """
    top = md.Topology()
    chain = top.add_chain()
    r1 = top.add_residue("ALA", chain)
    n_atom = top.add_atom("N", md.element.nitrogen, r1)
    ca_atom = top.add_atom("CA", md.element.carbon, r1)
    c_atom = top.add_atom("C", md.element.carbon, r1)
    # N–CA bond (key) + CA–C bond (normal)
    top.add_bond(n_atom, ca_atom)
    top.add_bond(ca_atom, c_atom)

    # Coordinates in nm. N at origin, CA at 0.6 nm (= 6.0 Å) on x.
    # CA-C is normal (0.15 nm = 1.5 Å).
    xyz = np.array([[
        [0.00, 0.00, 0.00],   # N
        [0.60, 0.00, 0.00],   # CA  (6.0 Å from N)
        [0.75, 0.00, 0.00],   # C   (1.5 Å from CA)
    ]], dtype=np.float32)
    box = np.array([[10.0, 10.0, 10.0]], dtype=np.float32)
    angles = np.array([[90.0, 90.0, 90.0]], dtype=np.float32)
    traj = md.Trajectory(xyz, top, time=[0.0],
                         unitcell_lengths=box, unitcell_angles=angles)

    broken_at_5 = es.intra_residue_bond_check(traj, threshold_angstrom=5.0)
    broken_at_10 = es.intra_residue_bond_check(traj, threshold_angstrom=10.0)
    print(
        f"\n[T15 synthetic 6 Å bond] thr=5.0: {len(broken_at_5)} broken; "
        f"thr=10.0: {len(broken_at_10)} broken"
    )
    assert len(broken_at_5) == 1, (
        f"thr=5.0: expected 1 broken (N-CA at 6.0 Å), got {broken_at_5}"
    )
    rec = broken_at_5[0]
    assert {rec["atom_a"], rec["atom_b"]} == {"N", "CA"}
    assert rec["distance"] == pytest.approx(6.0, abs=0.01), (
        f"expected ~6.0 Å, got {rec['distance']:.3f}"
    )
    assert broken_at_10 == [], (
        f"thr=10.0: expected 0 broken (6.0 < 10.0), got {broken_at_10}"
    )


# ---------------------------------------------------------------------------
# T16 — save_snapshots integration: healthy mini-trajectory
# ---------------------------------------------------------------------------

def test_t16_save_snapshots_integration_no_warnings(tmp_path, capsys):
    """Run ``save_snapshots`` on a healthy 2-residue mini-trajectory (no
    PBC defect) and verify:
        - no PR-21 WARN messages emitted to stderr
        - all snapshots written to disk
    """
    # Build a 2-residue ALA-ALA trajectory with intact backbone bonds.
    top = md.Topology()
    chain = top.add_chain()
    r1 = top.add_residue("ALA", chain)
    n1 = top.add_atom("N", md.element.nitrogen, r1)
    ca1 = top.add_atom("CA", md.element.carbon, r1)
    c1 = top.add_atom("C", md.element.carbon, r1)
    o1 = top.add_atom("O", md.element.oxygen, r1)
    r2 = top.add_residue("ALA", chain)
    n2 = top.add_atom("N", md.element.nitrogen, r2)
    ca2 = top.add_atom("CA", md.element.carbon, r2)
    c2 = top.add_atom("C", md.element.carbon, r2)
    o2 = top.add_atom("O", md.element.oxygen, r2)
    for a, b in [(n1, ca1), (ca1, c1), (c1, o1), (c1, n2),
                 (n2, ca2), (ca2, c2), (c2, o2)]:
        top.add_bond(a, b)

    # Healthy coordinates: every bond ~1.3-1.5 Å (well below 5 Å threshold)
    # Place residues 0.4 nm apart along x with realistic backbone geometry.
    xyz = np.array([[
        [0.00, 0.00, 0.00],   # N1
        [0.146, 0.00, 0.00],  # CA1 (1.46 Å from N1)
        [0.298, 0.00, 0.00],  # C1  (1.52 Å from CA1)
        [0.421, 0.00, 0.00],  # O1  (1.23 Å from C1)
        [0.431, 0.10, 0.00],  # N2 (within ~1.4 Å of C1; bond C1-N2)
        [0.577, 0.10, 0.00],  # CA2
        [0.729, 0.10, 0.00],  # C2
        [0.852, 0.10, 0.00],  # O2
    ]], dtype=np.float32)
    box = np.array([[5.0, 5.0, 5.0]], dtype=np.float32)
    angles = np.array([[90.0, 90.0, 90.0]], dtype=np.float32)
    # Two-frame trajectory (give save_snapshots more than 1 frame to write)
    xyz2 = np.concatenate([xyz, xyz], axis=0)
    box2 = np.concatenate([box, box], axis=0)
    angles2 = np.concatenate([angles, angles], axis=0)
    traj = md.Trajectory(xyz2, top, time=[0.0, 1.0],
                         unitcell_lengths=box2, unitcell_angles=angles2)

    # Patch the ncAA injector to a no-op (no MTR/etc here; avoid registry noise)
    with patch.object(es, "_inject_ncaa_conect", lambda _p: None):
        saved = es.save_snapshots(traj, [0, 1], str(tmp_path), "T16_healthy")

    captured = capsys.readouterr()
    pr21_warns = [
        line for line in captured.err.splitlines() if "[PR-21" in line
    ]
    print(
        f"\n[T16 healthy save_snapshots] saved={len(saved)} "
        f"PR-21 warnings={len(pr21_warns)}"
    )
    assert len(saved) == 2, f"expected 2 saved, got {len(saved)}"
    for p in saved:
        assert os.path.isfile(p), f"snapshot not written: {p}"
    assert pr21_warns == [], (
        f"unexpected PR-21 warnings on healthy trajectory: {pr21_warns}"
    )
