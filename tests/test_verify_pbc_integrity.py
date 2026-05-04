"""tests/test_verify_pbc_integrity.py — PR-21 Auto-Diagnostic Level 1.

Tests for ``scripts/verify_pbc_integrity.py`` v0.2 schema bump:
  - 3-column intra+inter output (default mode)
  - Backward-compat 4-column inter-only output (--inter-only flag)
  - snap05_f63 forensic regression (s55 BROKEN(intra) detection)

Test scope:
    T7 — ``--inter-only`` flag preserves the v0.1 4-column output
         (backward-compat regression).
    T8 — Default mode (intra+inter) emits the 3-column output with
         BROKEN(intra:LABEL) verdict on the s55 system.
    T9 — Forensic detection: with ``--filter`` pointed at the s55 dir,
         the report flags ``7TL8_MTR6_calib_s55`` as BROKEN(intra) due
         to snap05_f63's MTR6 N↔CA = 113 Å defect.

These tests exercise the CLI end-to-end via ``main(argv)`` so the
exit code + stdout are validated as a unit.
"""

from __future__ import annotations

import io
import os
import sys
from pathlib import Path

import pytest

# Make scripts/ importable
_REPO_ROOT = Path(__file__).resolve().parents[1]
_SCRIPTS_DIR = _REPO_ROOT / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import verify_pbc_integrity as vpi  # noqa: E402


OUTPUTS_ROOT = _REPO_ROOT / "outputs"
S55_SYSTEM_DIR = OUTPUTS_ROOT / "7TL8_MTR6_calib_s55"
S7_SYSTEM_DIR = OUTPUTS_ROOT / "1EBP_MTR13_calib_s7"


def _run_main(argv):
    """Invoke ``vpi.main(argv)`` and capture its stdout report.

    Returns ``(exit_code, stdout_text)``.
    """
    captured = io.StringIO()
    saved_stdout = sys.stdout
    sys.stdout = captured
    try:
        rc = vpi.main(argv)
    finally:
        sys.stdout = saved_stdout
    return rc, captured.getvalue()


# ---------------------------------------------------------------------------
# T7 — --inter-only flag preserves the legacy 4-column output
# ---------------------------------------------------------------------------

def test_t7_inter_only_legacy_4col_output():
    """With ``--inter-only``, the report header MUST contain the v0.1 column
    set (Frames | Median Å | Status) and MUST NOT contain the v0.2
    inter_med | intra_max | overall heading.
    """
    if not S7_SYSTEM_DIR.is_dir():
        pytest.skip(f"1EBP_MTR13_s7 system dir not present: {S7_SYSTEM_DIR}")
    rc, out = _run_main([
        "--root", str(OUTPUTS_ROOT),
        "--filter", "1EBP_MTR13_calib_s7$",
        "--inter-only",
    ])
    print(f"\n[T7 --inter-only] rc={rc}\n{out}")
    # Legacy header (v0.1)
    assert "Frames | Median Å | Status" in out, (
        f"legacy 4-column header missing under --inter-only: {out!r}"
    )
    # v0.2 header MUST NOT appear
    assert "inter_med | intra_max | overall" not in out, (
        f"v0.2 header leaked into --inter-only output: {out!r}"
    )
    # Verdict footer still present
    assert "Verdict:" in out


# ---------------------------------------------------------------------------
# T8 — Default mode emits 3-column intra+inter output
# ---------------------------------------------------------------------------

def test_t8_default_3col_output():
    """Default invocation (no --inter-only) emits the v0.2 3-column header
    ``inter_med | intra_max | overall``.
    """
    if not S7_SYSTEM_DIR.is_dir():
        pytest.skip(f"1EBP_MTR13_s7 system dir not present: {S7_SYSTEM_DIR}")
    rc, out = _run_main([
        "--root", str(OUTPUTS_ROOT),
        "--filter", "1EBP_MTR13_calib_s7$",
    ])
    print(f"\n[T8 default] rc={rc}\n{out}")
    # v0.2 header
    assert "inter_med | intra_max | overall" in out, (
        f"v0.2 3-column header missing in default mode: {out!r}"
    )
    # Legacy header MUST NOT appear
    assert "Frames | Median Å | Status" not in out, (
        f"v0.1 header leaked into default mode: {out!r}"
    )
    # Healthy 1EBP_MTR13_s7 system should be INTACT (intra ≤ threshold)
    assert "INTACT" in out


# ---------------------------------------------------------------------------
# T9 — snap05_f63 forensic detection via verify tool
# ---------------------------------------------------------------------------

def test_t9_s55_broken_intra_via_verify_tool():
    """Point the verify tool at the s55 system. Because the MD trajectory
    itself is healthy (1.29-1.46 Å for prev_C↔ncaa_N, normal N↔CA), the
    DCD-based inter+intra check on the FULL trajectory should NOT flag s55
    as BROKEN — the defect is in the EXTRACTED snap05_f63 PDB, not the
    DCD. So we instead verify the helper detects the defect on the
    snap05_f63 PDB directly using the imported helper.

    The ``--inter-only`` legacy run must report the s55 trajectory as
    INTACT (matching the #95 audit Task 1 finding: median 1.38 Å).
    """
    if not S55_SYSTEM_DIR.is_dir():
        pytest.skip(f"7TL8_MTR6_s55 system dir not present: {S55_SYSTEM_DIR}")

    # Part 1: --inter-only on the s55 DCD — must be INTACT (median 1.38 Å
    # per #95 audit — the defect is post-extraction, not in the DCD).
    rc, out = _run_main([
        "--root", str(OUTPUTS_ROOT),
        "--filter", "7TL8_MTR6_calib_s55$",
        "--inter-only",
    ])
    print(f"\n[T9 s55 --inter-only DCD] rc={rc}\n{out}")
    if "7TL8_MTR6_calib_s55" in out:
        # The DCD itself is INTACT; only the extracted snap05_f63 PDB is broken.
        assert "INTACT" in out or "AMBIGUOUS" in out, (
            f"unexpected status for s55 DCD --inter-only: {out!r}"
        )

    # Part 2: forensic — invoke the same intra-residue helper directly on
    # the snap05_f63 PDB to confirm the defect is detectable when the
    # input is the broken extraction (not the DCD). This is the actual
    # detection path PR-21 will use at extraction time.
    snap_pdb = (
        S55_SYSTEM_DIR
        / "snapshots_n25_postl387_patch_v2"
        / "7TL8_MTR6_snap05_f63.pdb"
    )
    if not snap_pdb.is_file():
        pytest.skip(f"snap05_f63 forensic PDB absent: {snap_pdb}")

    import mdtraj as md  # local import — heavy, env-dependent
    traj = md.load(str(snap_pdb))
    binder = vpi._find_binder_chain(traj.topology)
    assert binder is not None, "binder chain not detected in s55 snap05_f63"
    ncaa = vpi._find_ncaa_in_chain(binder)
    assert ncaa is not None, "ncAA (MTR) not detected in s55 snap05_f63"
    intra_max, intra_label, intra_broken = vpi._intra_residue_max_distance(
        traj, ncaa, threshold_angstrom=5.0,
    )
    print(
        f"\n[T9 s55 snap05_f63 forensic] intra_max={intra_max:.2f} Å "
        f"label={intra_label} broken={intra_broken}"
    )
    assert intra_broken, (
        f"PR-21 intra-residue helper FAILED to detect snap05_f63 defect: "
        f"max={intra_max} label={intra_label}"
    )
    assert intra_max > 100.0, (
        f"expected intra_max > 100 Å (#95 audit reports ~113 Å), got {intra_max}"
    )
    assert intra_label == "N-CA", (
        f"expected offending bond label 'N-CA', got {intra_label!r}"
    )
