"""tests/test_updd_py_cofactor_absent.py

Regression safety for cofactor-absent targets (e.g. Compstatin / C3 —
target_card with empty ``cofactor_residues``). The R-17 Stage-4 wiring
in ``UPDD.py`` must be a transparent no-op for any target whose card
declares no required cofactors; this file guards that invariant so
future targets ported into UPDD do not regress.

All tests are hermetic — UPDD.py is imported (not executed) and the
G1..G5 dispatchers are invoked directly with a synthesized cofactor-less
target_card placed under a tmp cards_dir. No subprocess is launched.
"""
from __future__ import annotations

import importlib
import json
import os
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def _ensure_repo_on_path() -> None:
    p = str(REPO_ROOT)
    if p not in sys.path:
        sys.path.insert(0, p)


def _import_updd_module():
    _ensure_repo_on_path()
    if "UPDD" in sys.modules:
        return sys.modules["UPDD"]
    import importlib as _il
    return _il.import_module("UPDD")


def _write_cofactor_absent_card(cards_dir: Path, target_id: str) -> Path:
    """Minimal v0.6.5 card with empty cofactor_residues — mimics
    Compstatin / generic cofactor-less targets."""
    card = {
        "schema_version": "0.6.5",
        "target_id": target_id,
        "target_description": "cofactor-absent test target",
        "target_chain": "A",
        "binder_chain": "B",
        "qm_partitioning": "sidechain_only",
        "target_contact_residues": [10, 20, 30],
        "binder_net_charge": 0,
        "target_iso_net_charge": 0,
        "target_iso_charge_rationale": "test",
        "numbering_convention": {"source": "MD", "md_to_crystal_offset": 0,
                                 "note": "identity mapping"},
        "cofactor_residues": [],
        "partition_rules": {
            "binder_residues_mode": "whole",
            "target_residues_mode": "whole",
            "glycine_handling": "whole_residue",
            "proline_handling": "whole_residue",
            "link_atom_model": "senn_thiel_2009_hcap",
            "cb_link_bond_length_A": 1.09,
        },
        "qm_atom_budget": {"max_n_qm": 600, "expected_n_qm_estimate": 400,
                           "expected_n_qm_range": [300, 500],
                           "enforcement": "hard_fail", "note": "test"},
        "provenance": {"pocket_identification_method": "test",
                       "contact_notes": "test", "review_date": "2026-04-19"},
    }
    path = cards_dir / f"{target_id}.json"
    path.write_text(json.dumps(card))
    return path


# ==========================================
# G2 dispatcher — empty cofactor_residues must no-op
# ==========================================
def test_g2_noop_for_cofactor_absent_card(tmp_path, monkeypatch):
    mod = _import_updd_module()

    cards_dir = tmp_path / "target_cards"
    cards_dir.mkdir()
    target_id = "COFLESS_A"
    _write_cofactor_absent_card(cards_dir, target_id)
    # Point UPDD_DIR at tmp_path so target_cards/ resolves to our fake.
    monkeypatch.setattr(mod, "UPDD_DIR", str(tmp_path))

    af2 = tmp_path / "af2_results"
    af2.mkdir()
    pdb = af2 / "design1_rank_001.pdb"
    pdb.write_text("ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00  0.00           C\nEND\n")

    report = mod._reinsert_cofactors_into_af2_outputs(
        str(af2), crystal_pdb="", target_id=target_id, target_chain="A",
    )
    assert report["status"] == "no_required_cofactors"
    assert report["n_processed"] == 0
    assert not (af2 / "_pre_reinsert_backup").exists()


# ==========================================
# G3 dispatcher — empty cofactor_residues must no-op
# ==========================================
def test_g3_noop_for_cofactor_absent_card(tmp_path, monkeypatch):
    mod = _import_updd_module()

    cards_dir = tmp_path / "target_cards"
    cards_dir.mkdir()
    target_id = "COFLESS_B"
    _write_cofactor_absent_card(cards_dir, target_id)
    monkeypatch.setattr(mod, "UPDD_DIR", str(tmp_path))

    report = mod._run_cofactor_parameterization(target_id)
    assert report["status"] == "noop"
    assert report.get("reason") == "no_required_cofactors"


# ==========================================
# G1 wrapper — run_target_preprocessing argv shape
# ==========================================
def test_run_target_preprocessing_omits_target_id_for_missing_card(
    tmp_path, monkeypatch
):
    """G1 passthrough must *not* inject ``--target_id`` when no card
    exists — preserves legacy behavior for ad-hoc targets."""
    mod = _import_updd_module()
    monkeypatch.setattr(mod, "UPDD_DIR", str(tmp_path))

    captured = {}

    def _fake_run_conda(script_path, args, description, env_key):
        captured["script_path"] = script_path
        captured["args"] = list(args)
        captured["description"] = description
        captured["env_key"] = env_key

    monkeypatch.setattr(mod, "run_conda_command", _fake_run_conda)
    target_pdb = tmp_path / "ad_hoc.pdb"
    target_pdb.write_text("ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00  0.00           C\nEND\n")

    out = mod.run_target_preprocessing(
        str(target_pdb), keep_mode="auto", chain_in="A", keep_res="",
        target_id="NO_SUCH_CARD_XYZ",
    )
    assert "--target_id" not in captured["args"]
    assert out.endswith("_clean.pdb")


def test_run_target_preprocessing_injects_target_id_when_card_present(
    tmp_path, monkeypatch
):
    """When the card exists, --target_id must be threaded into the argv."""
    mod = _import_updd_module()
    monkeypatch.setattr(mod, "UPDD_DIR", str(tmp_path))
    cards_dir = tmp_path / "target_cards"
    cards_dir.mkdir()
    _write_cofactor_absent_card(cards_dir, "COFLESS_C")

    captured = {}

    def _fake_run_conda(script_path, args, description, env_key):
        captured["args"] = list(args)

    monkeypatch.setattr(mod, "run_conda_command", _fake_run_conda)
    target_pdb = tmp_path / "foo.pdb"
    target_pdb.write_text("ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00  0.00           C\nEND\n")

    mod.run_target_preprocessing(
        str(target_pdb), keep_mode="auto", chain_in="A", keep_res="",
        target_id="COFLESS_C",
    )
    assert "--target_id" in captured["args"]
    idx = captured["args"].index("--target_id")
    assert captured["args"][idx + 1] == "COFLESS_C"


# ==========================================
# _run_parameterization — target_id kwarg accepted
# ==========================================
def test_run_parameterization_accepts_target_id_kwarg():
    """_run_parameterization must accept ``target_id`` without TypeError,
    even for wild-type (ncaa_def=None) flows."""
    mod = _import_updd_module()
    # Wild-type mode + no target_id → returns "" without touching disk.
    out = mod._run_parameterization(
        struct_dir="/nonexistent",
        param_dir="/nonexistent",
        ncaa_def=None,
        target_id="",
    )
    assert out == ""


# ==========================================
# _run_md_and_snapshots signature accepts target_id
# ==========================================
def test_run_md_and_snapshots_signature_has_target_id():
    """Regression guard: the helper must expose ``target_id`` kwarg so
    Step 10 snapshot extraction can forward --target_id for G5."""
    import inspect
    mod = _import_updd_module()
    sig = inspect.signature(mod._run_md_and_snapshots)
    assert "target_id" in sig.parameters
    assert sig.parameters["target_id"].default == ""
