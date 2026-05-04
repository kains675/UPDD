"""tests/test_updd_orchestrator_r17_wiring.py

Hermetic orchestrator-level smoke for the R-17 Stage-4 wiring inside
``UPDD.py``. These tests confirm that:

  1. The wiring helpers exist at the expected names and can be imported.
  2. UPDD.py passes ``--target_id`` down to each subprocess that owns a
     gate (G1 preprocess, G5 snapshot extraction).
  3. ``_reinsert_cofactors_into_af2_outputs`` (G2 dispatcher) is present,
     gracefully no-ops for cofactor-absent targets, and writes a backup
     for cofactor-present targets.
  4. ``_run_cofactor_parameterization`` (G3 dispatcher) can be called.
  5. The UPDD.py source references each G-gate entry point (static grep
     guard — regressions in the wiring are caught as source-level drift
     rather than only at runtime).

UPDD.py itself is NOT executed (it would pull in ColabFold / pyscf /
openmm). The tests only import top-level callables from the module after
providing harmless fallbacks for the heavy dashboard deps that the
module already tolerates via ImportError.

Python 3.8+ compatible. Matches the stub discipline used by
test_cofactor_injection.py and test_r17_cofactor_preservation.py.
"""
from __future__ import annotations

import importlib
import json
import os
import shutil
import sys
import types
from pathlib import Path
from typing import Dict, List

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
UPDD_PY_PATH = REPO_ROOT / "UPDD.py"


# ==========================================
# Helpers — import UPDD.py without running main()
# ==========================================
def _ensure_repo_on_path() -> None:
    """UPDD.py imports from utils/ via sys.path.append(UTILS_DIR). The
    project root must also be importable so ``import UPDD`` resolves."""
    repo = str(REPO_ROOT)
    if repo not in sys.path:
        sys.path.insert(0, repo)


def _import_updd_module():
    """Return the UPDD module, imported once per test session. Safe because
    UPDD.py guards its main() behind __name__ == "__main__" (implicit —
    main() is only called via the bottom-of-file ``if __name__``)."""
    _ensure_repo_on_path()
    # updd_dashboard / updd_config / updd_state may be unavailable; UPDD.py
    # already degrades to legacy mode via ImportError. No stubbing needed.
    if "UPDD" in sys.modules:
        return sys.modules["UPDD"]
    # UPDD.py is a top-level script — import by module name.
    return importlib.import_module("UPDD")


# ==========================================
# G-gate presence checks — static grep + import probe
# ==========================================
def test_updd_py_has_g1_target_id_passthrough():
    """G1: run_target_preprocessing accepts a target_id kwarg and forwards
    it to preprocess_target.py as ``--target_id``."""
    src = UPDD_PY_PATH.read_text(encoding="utf-8")
    assert "def run_target_preprocessing(" in src
    assert 'target_id' in src
    # The subprocess arg list must actually contain ``--target_id``.
    assert '"--target_id"' in src, "run_target_preprocessing should forward --target_id to preprocess_target.py"


def test_updd_py_has_g2_reinsert_helper():
    """G2: orchestrator-side dispatcher exists and calls into
    utils/reinsert_cofactors_post_af2.py."""
    src = UPDD_PY_PATH.read_text(encoding="utf-8")
    assert "def _reinsert_cofactors_into_af2_outputs(" in src
    assert "from reinsert_cofactors_post_af2 import reinsert_cofactors" in src


def test_updd_py_has_g3_parameterize_dispatch():
    """G3: Step 8 ncAA parameterization wraps a call into
    utils/parameterize_cofactor.ensure_cofactor_params."""
    src = UPDD_PY_PATH.read_text(encoding="utf-8")
    assert "def _run_cofactor_parameterization(" in src
    assert "from parameterize_cofactor import ensure_cofactor_params" in src


def test_updd_py_has_g5_snapshot_target_id():
    """G5: extract_snapshots.py is invoked with ``--target_id`` when the
    card is present."""
    src = UPDD_PY_PATH.read_text(encoding="utf-8")
    # _run_md_and_snapshots must thread target_id into snap_args.
    assert "def _run_md_and_snapshots(" in src
    assert 'snap_args += ["--target_id", target_id]' in src, (
        "extract_snapshots.py invocation must include --target_id for G5"
    )


def test_updd_py_imports_cleanly():
    """UPDD.py can be imported without executing main(). This catches
    NameError / ImportError / SyntaxError introduced by the Stage 4
    orchestrator wiring (e.g. missing sys.path setup for utils modules)."""
    mod = _import_updd_module()
    # Orchestrator-level helpers must be callable.
    assert callable(getattr(mod, "run_target_preprocessing", None))
    assert callable(getattr(mod, "_reinsert_cofactors_into_af2_outputs", None))
    assert callable(getattr(mod, "_run_cofactor_parameterization", None))
    assert callable(getattr(mod, "_run_md_and_snapshots", None))


# ==========================================
# G2 dispatcher behavior — cofactor-present path
# ==========================================
def _make_af2_rank_pdb(path: Path, chain: str, n_residues: int = 25) -> None:
    """Write a minimal ATOM-only PDB with a target chain of n_residues
    (Cα positions only — matches the shape reinsert_cofactors needs)."""
    lines: List[str] = []
    serial = 1
    for i in range(n_residues):
        x = float(i) * 3.8
        y = 0.0
        z = 0.0
        lines.append(
            f"ATOM  {serial:5d}  CA  ALA {chain}{i + 1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )
        serial += 1
    lines.append("END\n")
    path.write_text("".join(lines))


def _make_crystal_pdb(path: Path, chain: str, n_residues: int = 25,
                     include_hetatm: bool = True) -> None:
    """Same Cα layout as the AF2 stub plus declared HETATM rows (GNP 4
    atoms + MG 1 atom) so the Kabsch path has cofactors to transplant."""
    lines: List[str] = []
    serial = 1
    for i in range(n_residues):
        x = float(i) * 3.8
        y = 0.0
        z = 0.0
        lines.append(
            f"ATOM  {serial:5d}  CA  ALA {chain}{i + 1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        )
        serial += 1
    if include_hetatm:
        for i, name in enumerate(["P1", "O2", "N3", "C4"]):
            x = 5.0 + 0.5 * i
            lines.append(
                f"HETATM{serial:5d}  {name:<3s} GNP {chain} 201    "
                f"{x:8.3f}{5.0:8.3f}{5.0:8.3f}  1.00  0.00           {name[0]:>2s}\n"
            )
            serial += 1
        lines.append(
            f"HETATM{serial:5d}  MG  MG  {chain} 202    "
            f"{6.0:8.3f}{6.0:8.3f}{6.0:8.3f}  1.00  0.00           MG\n"
        )
    lines.append("END\n")
    path.write_text("".join(lines))


def test_g2_dispatcher_noop_when_target_id_missing(tmp_path):
    """No target_id → no backup directory, no crash."""
    mod = _import_updd_module()
    af2 = tmp_path / "af2_results"
    af2.mkdir()
    _make_af2_rank_pdb(af2 / "design1_rank_001.pdb", chain="A")
    report = mod._reinsert_cofactors_into_af2_outputs(
        str(af2), crystal_pdb="", target_id="", target_chain="A",
    )
    assert report["status"] == "no_target_id"
    assert not (af2 / "_pre_reinsert_backup").exists()


def test_g2_dispatcher_noop_when_card_missing(tmp_path):
    """Empty cards_dir → no_card + no backup dir."""
    mod = _import_updd_module()
    af2 = tmp_path / "af2_results"
    af2.mkdir()
    _make_af2_rank_pdb(af2 / "design1_rank_001.pdb", chain="A")
    # target_id provided but no card file in the repo's target_cards for
    # the name "NONEXISTENT_CARD_XYZ".
    report = mod._reinsert_cofactors_into_af2_outputs(
        str(af2), crystal_pdb="", target_id="NONEXISTENT_CARD_XYZ",
        target_chain="A",
    )
    assert report["status"] == "no_card"
    assert not (af2 / "_pre_reinsert_backup").exists()


def test_g2_dispatcher_processes_6wgn_card(tmp_path, monkeypatch):
    """Full G2 round-trip using the real 6WGN.json target_card:
    - AF2 PDB has only ATOM rows (no HETATM).
    - Crystal PDB has matching ATOM + HETATM (GNP + MG).
    - Dispatcher should write a backup + inject cofactors in-place.
    """
    mod = _import_updd_module()

    af2 = tmp_path / "af2_results"
    af2.mkdir()
    af2_pdb = af2 / "design1_rank_001.pdb"
    _make_af2_rank_pdb(af2_pdb, chain="A", n_residues=30)

    crystal = tmp_path / "work_clean.pdb"
    _make_crystal_pdb(crystal, chain="A", n_residues=30, include_hetatm=True)

    report = mod._reinsert_cofactors_into_af2_outputs(
        str(af2), crystal_pdb=str(crystal),
        target_id="6WGN", target_chain="A",
    )
    assert report["status"] == "ok", f"unexpected status: {report}"
    assert report["n_processed"] == 1
    # backup must exist and equal the original (no HETATM).
    backup = af2 / "_pre_reinsert_backup" / "design1_rank_001.pdb"
    assert backup.exists()
    backup_text = backup.read_text()
    assert "HETATM" not in backup_text
    # In-place PDB now carries cofactor HETATMs.
    rewritten = af2_pdb.read_text()
    assert "GNP" in rewritten
    assert " MG " in rewritten or "MG  " in rewritten
    # Report JSON written.
    rep_json = af2 / "_pre_reinsert_backup" / "reinsert_report.json"
    assert rep_json.exists()
    data = json.loads(rep_json.read_text())
    assert data["target_id"] == "6WGN"
    assert data["details"][0]["status"] == "ok"


def test_g2_dispatcher_idempotent_on_second_call(tmp_path):
    """If the backup already exists (resume path), second call reports
    skipped_cached and does not overwrite the rewritten PDB."""
    mod = _import_updd_module()

    af2 = tmp_path / "af2_results"
    af2.mkdir()
    af2_pdb = af2 / "design1_rank_001.pdb"
    _make_af2_rank_pdb(af2_pdb, chain="A", n_residues=30)
    crystal = tmp_path / "work_clean.pdb"
    _make_crystal_pdb(crystal, chain="A", n_residues=30, include_hetatm=True)

    first = mod._reinsert_cofactors_into_af2_outputs(
        str(af2), crystal_pdb=str(crystal),
        target_id="6WGN", target_chain="A",
    )
    assert first["status"] == "ok"
    first_content = af2_pdb.read_text()

    second = mod._reinsert_cofactors_into_af2_outputs(
        str(af2), crystal_pdb=str(crystal),
        target_id="6WGN", target_chain="A",
    )
    assert second["status"] in ("skip_cached", "ok")
    assert second["n_skipped"] >= 1
    # Content unchanged on second invocation.
    assert af2_pdb.read_text() == first_content


# ==========================================
# G3 dispatcher behavior
# ==========================================
def test_g3_dispatcher_noop_when_target_id_missing():
    mod = _import_updd_module()
    report = mod._run_cofactor_parameterization("")
    assert report == {"status": "noop", "reason": "no_target_id"}


def test_g3_dispatcher_noop_when_card_missing():
    mod = _import_updd_module()
    report = mod._run_cofactor_parameterization("NO_SUCH_CARD_XYZ")
    assert report == {"status": "noop", "reason": "no_card"}


def test_g3_dispatcher_ok_for_6wgn(tmp_path, monkeypatch):
    """G3 against the real 6WGN.json card. Mg²⁺ resolves as ion_builtin
    (Li-Merz 12-6-4); GNP has declared ff_parameters paths that may or
    may not exist in the working tree. Acceptable terminal statuses are
    ``ok`` (cache_hit / ion_builtin / regenerated) or
    ``warn_missing`` (GNP missing cache, antechamber opt-in off).
    """
    mod = _import_updd_module()
    # Ensure antechamber opt-in is OFF so the test is deterministic.
    monkeypatch.delenv("UPDD_ALLOW_ANTECHAMBER", raising=False)
    report = mod._run_cofactor_parameterization("6WGN")
    assert report["status"] in ("ok", "warn_missing"), report


# ==========================================
# Regression — static wiring presence
# ==========================================
def test_updd_py_wires_g1_in_main_flow():
    src = UPDD_PY_PATH.read_text(encoding="utf-8")
    # Main flow must actually call run_target_preprocessing with target_id.
    assert "run_target_preprocessing(" in src
    assert "target_id=base_name" in src, (
        "main() must invoke run_target_preprocessing(..., target_id=base_name) for G1"
    )


def test_updd_py_wires_g2_step_6_6():
    src = UPDD_PY_PATH.read_text(encoding="utf-8")
    assert "_reinsert_cofactors_into_af2_outputs(" in src
    assert "Step 6.6" in src, "Step 6.6 G2 print_step header expected"


def test_updd_py_wires_g3_in_parameterization():
    src = UPDD_PY_PATH.read_text(encoding="utf-8")
    assert "_run_cofactor_parameterization" in src


def test_updd_py_wires_g5_via_md_snapshot_helper():
    src = UPDD_PY_PATH.read_text(encoding="utf-8")
    assert "target_id=base_name" in src
    assert 'snap_args += ["--target_id", target_id]' in src
