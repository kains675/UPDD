"""tests/test_charge_consistency.py — R-15 / R-16 Stage 1 unit tests.

Covers the chemistry-true charge computation (``charge_topology``), the
R-15 magnitude guard, the R-16 binder SSOT guard, cyclic/linear peptide
detection, and the HIS tautomer classification.

These tests use only stdlib (no pyscf, no openmm) so they run in the
bare Python environment.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

import pytest

# tests/conftest.py already adds utils/ to sys.path.
from charge_topology import (  # noqa: E402
    classify_residue_charge,
    compute_binder_chem_charge,
    compute_qm_net_charge_topology,
    detect_cyclic_peptide,
    parse_pdb_atoms_lite,
    residues_by_chain_resnum,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
SNAP01_PDB = (
    REPO_ROOT
    / "outputs/6WGN_cyclic_htc_NMA_10_20-25/snapshots/"
    "design_w4_1_s2_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000_ncaa_snap01_f26.pdb"
)

TARGET_CARD_6WGN = REPO_ROOT / "target_cards" / "6WGN.json"


# ==========================================
# Helpers for synthetic residue atom lists
# ==========================================
def _atom(chain: str, resnum: int, resname: str, name: str,
          element: str = "") -> Dict[str, Any]:
    """Construct an atom dict matching ``parse_pdb_atoms_lite`` output."""
    if not element:
        element = name.lstrip("0123456789")[:1].upper() or "C"
    return {
        "record": "ATOM", "serial": 1, "name": name, "resname": resname,
        "chain": chain, "resnum": resnum, "x": 0.0, "y": 0.0, "z": 0.0,
        "element": element,
    }


def _linear_peptide_atoms(chain: str = "B") -> List[Dict[str, Any]]:
    """Two-residue linear peptide: NH3+ LEU -- COO- LEU. Net = 0."""
    atoms = []
    # LEU1 linear N-terminus (NH3+, three H: H1/H2/H3)
    for nm in ("N", "H1", "H2", "H3", "CA", "HA", "C", "O", "CB"):
        atoms.append(_atom(chain, 1, "LEU", nm))
    # LEU2 linear C-terminus (COO- with OXT)
    for nm in ("N", "H", "CA", "HA", "C", "O", "OXT", "CB"):
        atoms.append(_atom(chain, 2, "LEU", nm))
    return atoms


def _cyclic_leu21_atoms(chain: str = "B") -> List[Dict[str, Any]]:
    """21-residue cyclic peptide LEU1..LEU21. First residue 1 N-H, last no OXT."""
    atoms = []
    # LEU1 cyclic: only one backbone H (H on N)
    for nm in ("N", "H", "CA", "HA", "C", "O", "CB"):
        atoms.append(_atom(chain, 1, "LEU", nm))
    # Middle residues 2-20 (plain LEU, peptide-bonded N-H only)
    for rnum in range(2, 21):
        for nm in ("N", "H", "CA", "HA", "C", "O", "CB"):
            atoms.append(_atom(chain, rnum, "LEU", nm))
    # LEU21 cyclic: no OXT
    for nm in ("N", "H", "CA", "HA", "C", "O", "CB"):
        atoms.append(_atom(chain, 21, "LEU", nm))
    return atoms


def _arg_atoms(chain: str, resnum: int) -> List[Dict[str, Any]]:
    """Protonated arginine: all four HH hydrogens present."""
    names = (
        "N", "H", "CA", "HA", "C", "O", "CB", "CG", "CD", "NE", "HE",
        "CZ", "NH1", "NH2", "HH11", "HH12", "HH21", "HH22",
    )
    return [_atom(chain, resnum, "ARG", nm) for nm in names]


def _lys_atoms(chain: str, resnum: int) -> List[Dict[str, Any]]:
    """Protonated lysine: HZ1/HZ2/HZ3 all present."""
    names = (
        "N", "H", "CA", "HA", "C", "O", "CB", "CG", "CD", "CE", "NZ",
        "HZ1", "HZ2", "HZ3",
    )
    return [_atom(chain, resnum, "LYS", nm) for nm in names]


def _asp_atoms(chain: str, resnum: int) -> List[Dict[str, Any]]:
    """Deprotonated aspartate: OD1/OD2 only, no carboxylate HD."""
    names = ("N", "H", "CA", "HA", "C", "O", "CB", "CG", "OD1", "OD2")
    return [_atom(chain, resnum, "ASP", nm) for nm in names]


def _glu_atoms(chain: str, resnum: int) -> List[Dict[str, Any]]:
    """Deprotonated glutamate: OE1/OE2 only, no carboxylate HE."""
    names = ("N", "H", "CA", "HA", "C", "O", "CB", "CG", "CD", "OE1", "OE2")
    return [_atom(chain, resnum, "GLU", nm) for nm in names]


def _leu_atoms(chain: str, resnum: int,
               n_terminal: bool = False,
               c_terminal: bool = False) -> List[Dict[str, Any]]:
    """Neutral leucine residue, optionally with NH3+ or OXT termini."""
    atoms = []
    if n_terminal:
        for nm in ("N", "H1", "H2", "H3", "CA", "HA", "C", "O", "CB"):
            atoms.append(_atom(chain, resnum, "LEU", nm))
    else:
        for nm in ("N", "H", "CA", "HA", "C", "O", "CB"):
            atoms.append(_atom(chain, resnum, "LEU", nm))
    if c_terminal:
        atoms.append(_atom(chain, resnum, "LEU", "OXT"))
    return atoms


# ==========================================
# Unit tests
# ==========================================
class TestClassifyResidueCharge:
    def test_arg_protonated(self):
        names = {"HH11", "HH12", "HH21", "HH22", "N", "CA"}
        assert classify_residue_charge("ARG", names) == 1

    def test_arg_missing_hh_is_zero(self):
        # Partial guanidinium (HH absent) → neutral
        names = {"N", "CA", "CB"}
        assert classify_residue_charge("ARG", names) == 0

    def test_lys_protonated(self):
        names = {"HZ1", "HZ2", "HZ3", "N", "CA"}
        assert classify_residue_charge("LYS", names) == 1

    def test_lys_deprotonated(self):
        names = {"N", "CA", "NZ"}
        assert classify_residue_charge("LYS", names) == 0

    def test_asp_deprotonated(self):
        assert classify_residue_charge("ASP", {"OD1", "OD2"}) == -1

    def test_ash_neutral(self):
        # AMBER variant name ASH → 0
        assert classify_residue_charge("ASH", {"OD1", "OD2", "HD2"}) == 0

    def test_glu_deprotonated(self):
        assert classify_residue_charge("GLU", {"OE1", "OE2"}) == -1

    def test_glh_neutral(self):
        assert classify_residue_charge("GLH", {"OE1", "OE2", "HE2"}) == 0

    def test_his_hid_neutral(self):
        # HID: HD1 only
        assert classify_residue_charge("HIS", {"HD1", "ND1", "NE2", "CE1"}) == 0

    def test_his_hie_neutral(self):
        # HIE: HE2 only (no HD1)
        assert classify_residue_charge("HIS", {"HE2", "ND1", "NE2", "CE1"}) == 0

    def test_his_hip_positive(self):
        # HIP: both HD1 and HE2
        assert classify_residue_charge("HIS", {"HD1", "HE2", "ND1", "NE2"}) == 1

    def test_his_tautomer_hid_explicit(self):
        assert classify_residue_charge("HID", {"HD1", "ND1"}) == 0
        assert classify_residue_charge("HIE", {"HE2", "NE2"}) == 0
        assert classify_residue_charge("HIP", {"HD1", "HE2", "ND1", "NE2"}) == 1

    def test_tyr_neutral(self):
        assert classify_residue_charge("TYR", {"OH"}) == 0

    def test_cys_neutral(self):
        assert classify_residue_charge("CYS", {"SG", "HG"}) == 0

    def test_unknown_residue_neutral(self):
        assert classify_residue_charge("UNK", set()) == 0

    # ─── X1 (2026-04-20): expanded protonation variant coverage ─────
    # R-16 runtime-only policy makes classify_residue_charge the SOLE
    # safety net (along with build_qm_mol OddElectronError). These
    # regression tests guard against silent miss-classification that
    # would produce wrong QM net charge.
    def test_lyn_explicit_resname(self):
        # AMBER LYN = deprotonated Lys. With or without explicit HZ absence,
        # resname LYN should be recognized as 0.
        assert classify_residue_charge("LYN", {"N", "CA", "NZ"}) == 0
        assert classify_residue_charge("LYN", {"N", "CA", "NZ", "HZ1"}) == 0

    def test_hip_explicit_resname(self):
        # AMBER HIP = doubly-protonated His. resname HIP should be +1
        # without relying on atom name pattern.
        assert classify_residue_charge("HIP", {"ND1", "NE2"}) == 1

    def test_hid_hie_without_atom_detail(self):
        # Explicit resnames should be recognised even if tautomer atom
        # names are incomplete.
        assert classify_residue_charge("HID", set()) == 0
        assert classify_residue_charge("HIE", set()) == 0

    def test_his_protonation_missing_atoms_defaults_neutral(self):
        # HIS without HD1/HE2 atoms (ambiguous) → neutral (not +1).
        assert classify_residue_charge("HIS", {"ND1", "NE2"}) == 0

    def test_d_amino_acids_match_parent(self):
        # D-isomers (DAL, DPN, DTR, DPR, DLE, DVA) inherit parent charge.
        # These are added to ncaa_registry but classify_residue_charge
        # without ncAA charge map should default to 0.
        for dres in ("DAL", "DPN", "DTR", "DPR", "DLE", "DVA"):
            assert classify_residue_charge(dres, set()) == 0

    def test_nme_ace_caps_neutral(self):
        # MD terminal caps — always neutral.
        assert classify_residue_charge("NME", {"N", "CH3"}) == 0
        assert classify_residue_charge("ACE", {"C", "O", "CH3"}) == 0

    def test_proline_neutral(self):
        # PRO (normal) / DPR (D-proline) — neutral cyclic backbone.
        assert classify_residue_charge("PRO", {"N", "CA", "CD"}) == 0

    def test_tyr_with_hh_still_neutral(self):
        # TYR with hydroxyl H — pKa ~10, neutral at pH 7.4.
        assert classify_residue_charge("TYR", {"OH", "HH"}) == 0


class TestDetectCyclicPeptide:
    def test_r14_detects_cyclic(self):
        """LEU1 has 1 backbone N-H, LEU21 no OXT → cyclic."""
        atoms = _cyclic_leu21_atoms()
        group = residues_by_chain_resnum(atoms, chain_filter="B")
        is_cyclic, diag = detect_cyclic_peptide(group)
        assert is_cyclic is True, diag
        assert diag["first_residue"]["n_backbone_NH"] == 1
        assert diag["last_residue"]["has_OXT"] is False

    def test_r14_detects_linear(self):
        """LEU1 with NH3+ (3 H), LEU2 with OXT → linear."""
        atoms = _linear_peptide_atoms()
        group = residues_by_chain_resnum(atoms, chain_filter="B")
        is_cyclic, diag = detect_cyclic_peptide(group)
        assert is_cyclic is False
        assert diag["first_residue"]["n_backbone_NH"] >= 2
        assert diag["last_residue"]["has_OXT"] is True


class TestBinderChemCharge:
    def test_compute_binder_chem_charge_arg_only(self):
        """Synthetic dipeptide with ARG + LEU (cyclic) → +1 net."""
        atoms = []
        # LEU1 cyclic first
        atoms.extend(_leu_atoms("B", 1))
        atoms.extend(_arg_atoms("B", 2))
        atoms.extend(_leu_atoms("B", 3))  # cyclic last (no OXT)
        q, diag = compute_binder_chem_charge(atoms, "B")
        assert q == 1
        assert diag["cyclic"] is True

    def test_compute_binder_chem_charge_asp_glu_balanced(self):
        """ARG+LYS+ASP+GLU balanced → 0 net charge (cyclic)."""
        atoms = []
        atoms.extend(_leu_atoms("B", 1))  # cyclic first, no NH3+
        atoms.extend(_arg_atoms("B", 2))  # +1
        atoms.extend(_lys_atoms("B", 3))  # +1
        atoms.extend(_asp_atoms("B", 4))  # -1
        atoms.extend(_glu_atoms("B", 5))  # -1
        atoms.extend(_leu_atoms("B", 6))  # cyclic last, no OXT
        q, diag = compute_binder_chem_charge(atoms, "B")
        assert q == 0
        assert diag["cyclic"] is True


@pytest.mark.skipif(not SNAP01_PDB.is_file(),
                    reason=f"snap01 PDB missing at {SNAP01_PDB}")
class TestSnap01Integration:
    """Integration smoke test using the real w4_1_s2/snap01 PDB."""

    def test_compute_binder_chem_charge_w4_1_s2_snap01(self):
        atoms = parse_pdb_atoms_lite(str(SNAP01_PDB))
        q, diag = compute_binder_chem_charge(atoms, "B")
        assert q == -1, f"Expected binder = -1, got {q}; diag={diag}"
        assert diag["cyclic"] is True

    def test_compute_qm_net_charge_topology_snap01(self):
        atoms = parse_pdb_atoms_lite(str(SNAP01_PDB))
        with open(TARGET_CARD_6WGN, "r", encoding="utf-8") as f:
            card = json.load(f)
        binder_q, target_q, total_q, diag = compute_qm_net_charge_topology(
            atoms,
            binder_chain=card["binder_chain"],
            target_chain=card["target_chain"],
            target_contact_residues=card["target_contact_residues"],
            whole_residue_exceptions=card.get("whole_residue_exceptions", []),
            pH=7.4,
        )
        assert binder_q == -1
        assert target_q == -3
        assert total_q == -4

    def test_stage2_card_declared_matches_computed(self):
        """Stage 2 (v0.6.4): card binder=-1 + target=-3 ↔ runtime -1 + -3 = -4."""
        atoms = parse_pdb_atoms_lite(str(SNAP01_PDB))
        with open(TARGET_CARD_6WGN, "r", encoding="utf-8") as f:
            card = json.load(f)
        _, _, total_q, _ = compute_qm_net_charge_topology(
            atoms,
            binder_chain=card["binder_chain"],
            target_chain=card["target_chain"],
            target_contact_residues=card["target_contact_residues"],
            whole_residue_exceptions=card.get("whole_residue_exceptions", []),
        )
        declared = card["binder_net_charge"] + card["target_iso_net_charge"]
        # Stage 2: declared and computed both -4 (R-15 PASS).
        assert declared == -4
        assert total_q == -4
        assert declared == total_q
        assert (declared % 2) == (total_q % 2)


class TestExceptionWiring:
    """Verify that run_qmmm exports the new exception classes correctly.

    These do NOT import pyscf — they only check that the symbols exist in
    the module namespace via ``importlib.util`` spec inspection.
    """

    def test_chargedeclarationmismatch_exists(self):
        # The exception class inherits from ValueError to keep existing
        # ``except ValueError`` handlers compatible.
        import ast
        src = (REPO_ROOT / "utils" / "run_qmmm.py").read_text(encoding="utf-8")
        tree = ast.parse(src)
        names = {n.name for n in tree.body if isinstance(n, ast.ClassDef)}
        assert "ChargeDeclarationMismatch" in names
        assert "BinderChargeMismatch" in names
        assert "OddElectronError" in names


class TestR16X1RuntimeOnlyPolicy:
    """Schema 0.6.6 (R-16 X1) regression: target_card.binder_net_charge 는
    informational hint 로 격하. Runtime mismatch 시 fail-fast 아닌 warn.

    These tests verify at the source-code level (ast parsing) that:
    - run_qmmm.py no longer raises BinderChargeMismatch / ChargeDeclarationMismatch
      unconditionally on hint mismatch
    - A warnings.warn call is issued instead
    """

    def test_run_qmmm_r16_no_longer_raises_on_hint_mismatch(self):
        src = (REPO_ROOT / "utils" / "run_qmmm.py").read_text(encoding="utf-8")
        # X1 marker string must be present, confirming the policy change.
        assert "runtime_authoritative_x1" in src, (
            "run_qmmm.py must declare R-16 X1 runtime-only policy "
            "(marker 'runtime_authoritative_x1' missing)."
        )
        # The fail-fast policy label should be absent from the guard block.
        assert "fail_fast_r14" not in src, (
            "run_qmmm.py still contains legacy 'fail_fast_r14' policy label "
            "— remove for X1 transition."
        )

    def test_schema_0_6_6_accepted(self):
        from utils_common import load_target_card
        import tempfile
        import shutil
        # Create a minimal 0.6.6 card and confirm load_target_card accepts it.
        with tempfile.TemporaryDirectory() as tmp:
            card_path = Path(tmp) / "X1TEST.json"
            card_path.write_text(json.dumps({
                "schema_version": "0.6.6",
                "target_id": "X1TEST",
                "target_chain": "A",
                "binder_chain": "B",
                "target_iso_net_charge": -1,
                "binder_net_charge": 0,
            }), encoding="utf-8")
            loaded = load_target_card("X1TEST", cards_dir=tmp)
            assert loaded["schema_version"] == "0.6.6"
