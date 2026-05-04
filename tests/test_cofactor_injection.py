"""tests/test_cofactor_injection.py — Interim cofactor runtime injection
diagnostic (SciVal 8차 verdict §7.3).

Tests for ``run_qmmm._inject_cofactor_point_charges_from_crystal`` — the
gated probe that appends cofactor atoms from the crystal PDB as extra MM
point charges when ``UPDD_COFACTOR_INJECT_DIAGNOSTIC=1``.

This module avoids the real pyscf/openmm imports inside ``run_qmmm`` by
stubbing the relevant submodules in ``sys.modules`` *before* importing
``run_qmmm``. The injection function itself only depends on
``numpy``, ``os``, ``re``, and ``utils_common.parse_pdb_atom_line``.
"""
from __future__ import annotations

import os
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
TARGET_CARD_6WGN_PATH = REPO_ROOT / "target_cards" / "6WGN.json"
CRYSTAL_PDB_6WGN = REPO_ROOT / "target" / "6WGN.pdb"


# ==========================================
# Stub heavy deps so `run_qmmm` imports without pyscf/openmm
# ==========================================
def _install_import_stubs() -> None:
    """Insert minimal stub modules for ``pyscf`` and ``openmm`` so we can
    import ``utils.run_qmmm`` without the real heavy dependencies.

    These stubs satisfy only the *import-time* references in ``run_qmmm``;
    the injection function we test does not call any pyscf/openmm API.
    """
    if "pyscf" not in sys.modules:
        pyscf = types.ModuleType("pyscf")
        pyscf.__config__ = types.SimpleNamespace(min_ao_blksize=64)
        pyscf.gto = types.ModuleType("pyscf.gto")
        pyscf.dft = types.ModuleType("pyscf.dft")
        pyscf.lib = types.ModuleType("pyscf.lib")
        pyscf.lib.num_threads = lambda *_a, **_k: None
        pyscf.qmmm = types.ModuleType("pyscf.qmmm")
        pyscf.qmmm.itrf = types.ModuleType("pyscf.qmmm.itrf")
        sys.modules["pyscf"] = pyscf
        sys.modules["pyscf.gto"] = pyscf.gto
        sys.modules["pyscf.dft"] = pyscf.dft
        sys.modules["pyscf.lib"] = pyscf.lib
        sys.modules["pyscf.qmmm"] = pyscf.qmmm
        sys.modules["pyscf.qmmm.itrf"] = pyscf.qmmm.itrf

    if "openmm" not in sys.modules:
        openmm = types.ModuleType("openmm")
        openmm.unit = types.ModuleType("openmm.unit")
        openmm.app = types.ModuleType("openmm.app")
        openmm.app.ForceField = type("ForceField", (), {})
        openmm.app.PDBFile = type("PDBFile", (), {})
        sys.modules["openmm"] = openmm
        sys.modules["openmm.unit"] = openmm.unit
        sys.modules["openmm.app"] = openmm.app


_install_import_stubs()
# After stubbing: import run_qmmm lazily. conftest.py already puts
# utils/ on sys.path.
import run_qmmm  # noqa: E402


# ==========================================
# Fixtures / helpers
# ==========================================
TARGET_CARD_6WGN = {
    "target_id": "6WGN",
    "cofactor_residues": [
        {"resname": "GNP", "chain": "A", "resnum": 201,
         "treatment": "mm_point_charge", "charge": -4},
        {"resname": "MG",  "chain": "A", "resnum": 202,
         "treatment": "mm_point_charge", "charge": 2},
    ],
}


@pytest.fixture
def unset_diagnostic_env(monkeypatch):
    monkeypatch.delenv("UPDD_COFACTOR_INJECT_DIAGNOSTIC", raising=False)
    yield


@pytest.fixture
def diagnostic_env_on(monkeypatch):
    monkeypatch.setenv("UPDD_COFACTOR_INJECT_DIAGNOSTIC", "1")
    yield


@pytest.fixture
def empty_mm_inputs():
    """Return ``(mm_atoms_in, mm_charges_in)`` both empty — the injection
    function must handle an empty MM arrangement gracefully."""
    import numpy as np
    return [], np.zeros((0, 4), dtype=float)


# ==========================================
# Tests
# ==========================================
class TestInjectNoopWhenEnvUnset:
    def test_noop_when_env_absent(self, unset_diagnostic_env, empty_mm_inputs):
        """With the env var absent, function must be a strict no-op regardless
        of other inputs and must not emit injection diagnostics."""
        mm_atoms_in, mm_charges_in = empty_mm_inputs
        out_atoms, out_charges, meta = (
            run_qmmm._inject_cofactor_point_charges_from_crystal(
                target_card=TARGET_CARD_6WGN,
                mm_atoms_in=mm_atoms_in,
                mm_charges_in=mm_charges_in,
                crystal_pdb_path=str(CRYSTAL_PDB_6WGN),
                enabled_via_env_var=False,
            )
        )
        assert out_atoms is mm_atoms_in
        assert out_charges is mm_charges_in
        assert meta["enabled"] is False
        assert meta["injected_atoms"] == []
        assert meta["charge_sum"] is None
        assert meta["source_pdb"] is None

    def test_noop_with_target_card_none(self, empty_mm_inputs):
        """Even with the env var True, a missing target_card must return
        unchanged inputs (safe fallback, no crash)."""
        mm_atoms_in, mm_charges_in = empty_mm_inputs
        out_atoms, out_charges, meta = (
            run_qmmm._inject_cofactor_point_charges_from_crystal(
                target_card=None,
                mm_atoms_in=mm_atoms_in,
                mm_charges_in=mm_charges_in,
                crystal_pdb_path=str(CRYSTAL_PDB_6WGN),
                enabled_via_env_var=True,
            )
        )
        assert out_atoms is mm_atoms_in
        assert out_charges is mm_charges_in
        assert meta["fallback_reason"] == "target_card_none"


class TestInjectReadsCrystalGnpMg:
    @pytest.mark.skipif(not CRYSTAL_PDB_6WGN.exists(),
                        reason="target/6WGN.pdb not available in this checkout")
    def test_gnp_mg_appended_with_33_atoms_and_charge_sum_minus_two(
        self, diagnostic_env_on, empty_mm_inputs,
    ):
        """With env var on, the real crystal PDB must contribute 32 GNP
        heavy atoms + 1 MG atom = 33 extra MM charges; charge sum = -2."""
        import numpy as np
        mm_atoms_in, mm_charges_in = empty_mm_inputs
        out_atoms, out_charges, meta = (
            run_qmmm._inject_cofactor_point_charges_from_crystal(
                target_card=TARGET_CARD_6WGN,
                mm_atoms_in=mm_atoms_in,
                mm_charges_in=mm_charges_in,
                crystal_pdb_path=str(CRYSTAL_PDB_6WGN),
                enabled_via_env_var=True,
            )
        )
        assert len(out_atoms) == 33, f"expected 33 atoms, got {len(out_atoms)}"
        assert out_charges.shape == (33, 4)

        # Charge sum ≈ -2 (within float tolerance — GNP -4 distributed across
        # 32 heavy atoms + MG +2 on 1 atom).
        charge_sum = float(np.sum(out_charges[:, 3]))
        assert abs(charge_sum - (-2.0)) < 1e-9, (
            f"expected net -2, got {charge_sum}"
        )

        resnames = {a["resname"] for a in out_atoms}
        assert resnames == {"GNP", "MG"}, (
            f"unexpected resnames: {resnames}"
        )

        # Per-residue breakdown diagnostics.
        assert meta["source_pdb"] == str(CRYSTAL_PDB_6WGN)
        assert meta["fallback_reason"] is None
        assert meta["charge_sum"] == pytest.approx(-2.0, abs=1e-9)
        breakdown_names = {b["resname"] for b in meta["injected_residues_breakdown"]}
        assert breakdown_names == {"GNP", "MG"}

        # GNP: 32 heavy atoms at -4/32 each.
        gnp_entry = next(
            b for b in meta["injected_residues_breakdown"] if b["resname"] == "GNP"
        )
        assert gnp_entry["n_heavy_atoms"] == 32
        assert gnp_entry["charge_total"] == pytest.approx(-4.0)
        assert gnp_entry["per_atom_charge"] == pytest.approx(-4.0 / 32.0)

        # MG: 1 atom at +2.
        mg_entry = next(
            b for b in meta["injected_residues_breakdown"] if b["resname"] == "MG"
        )
        assert mg_entry["n_heavy_atoms"] == 1
        assert mg_entry["charge_total"] == pytest.approx(2.0)
        assert mg_entry["per_atom_charge"] == pytest.approx(2.0)


class TestInjectGracefulFallbackMissingCrystal:
    def test_missing_crystal_path_is_no_op(self, diagnostic_env_on, empty_mm_inputs, tmp_path):
        """A non-existent crystal PDB must return inputs unchanged and flag
        fallback_reason=crystal_pdb_missing (no crash)."""
        mm_atoms_in, mm_charges_in = empty_mm_inputs
        bogus_path = str(tmp_path / "does_not_exist.pdb")
        out_atoms, out_charges, meta = (
            run_qmmm._inject_cofactor_point_charges_from_crystal(
                target_card=TARGET_CARD_6WGN,
                mm_atoms_in=mm_atoms_in,
                mm_charges_in=mm_charges_in,
                crystal_pdb_path=bogus_path,
                enabled_via_env_var=True,
            )
        )
        assert out_atoms is mm_atoms_in
        assert out_charges is mm_charges_in
        assert meta["fallback_reason"] == "crystal_pdb_missing"
        assert meta["injected_atoms"] == []
        assert meta["charge_sum"] is None

    def test_none_crystal_path_is_no_op(self, diagnostic_env_on, empty_mm_inputs):
        mm_atoms_in, mm_charges_in = empty_mm_inputs
        out_atoms, out_charges, meta = (
            run_qmmm._inject_cofactor_point_charges_from_crystal(
                target_card=TARGET_CARD_6WGN,
                mm_atoms_in=mm_atoms_in,
                mm_charges_in=mm_charges_in,
                crystal_pdb_path=None,
                enabled_via_env_var=True,
            )
        )
        assert out_atoms is mm_atoms_in
        assert out_charges is mm_charges_in
        assert meta["fallback_reason"] == "crystal_pdb_missing"


class TestInjectGracefulFallbackCofactorNotInPdb:
    def test_synthetic_pdb_without_declared_cofactor_noop(
        self, diagnostic_env_on, empty_mm_inputs, tmp_path,
    ):
        """Crystal PDB exists but contains no declared cofactor residues —
        function must fall back gracefully, not crash."""
        synthetic_pdb = tmp_path / "synthetic_no_cofactor.pdb"
        # Build a minimal protein-only PDB (chain A, 1 residue ALA).
        synthetic_pdb.write_text(
            "HEADER    SYNTHETIC\n"
            "ATOM      1  N   ALA A   1      "
            "  0.000   0.000   0.000  1.00  0.00           N  \n"
            "ATOM      2  CA  ALA A   1      "
            "  1.500   0.000   0.000  1.00  0.00           C  \n"
            "ATOM      3  C   ALA A   1      "
            "  2.000   1.400   0.000  1.00  0.00           C  \n"
            "ATOM      4  O   ALA A   1      "
            "  1.300   2.400   0.000  1.00  0.00           O  \n"
            "END\n"
        )
        mm_atoms_in, mm_charges_in = empty_mm_inputs
        out_atoms, out_charges, meta = (
            run_qmmm._inject_cofactor_point_charges_from_crystal(
                target_card=TARGET_CARD_6WGN,
                mm_atoms_in=mm_atoms_in,
                mm_charges_in=mm_charges_in,
                crystal_pdb_path=str(synthetic_pdb),
                enabled_via_env_var=True,
            )
        )
        assert out_atoms is mm_atoms_in
        assert out_charges is mm_charges_in
        assert meta["fallback_reason"] == "no_declared_cofactor_matched"
        assert meta["injected_atoms"] == []
        assert meta["charge_sum"] is None

    def test_target_card_without_cofactor_residues_noop(
        self, diagnostic_env_on, empty_mm_inputs,
    ):
        """If target_card lacks ``cofactor_residues`` (or it is empty) we must
        fall back gracefully rather than attempt PDB parsing."""
        mm_atoms_in, mm_charges_in = empty_mm_inputs
        bare_card = {"target_id": "6WGN", "cofactor_residues": []}
        out_atoms, out_charges, meta = (
            run_qmmm._inject_cofactor_point_charges_from_crystal(
                target_card=bare_card,
                mm_atoms_in=mm_atoms_in,
                mm_charges_in=mm_charges_in,
                crystal_pdb_path=str(CRYSTAL_PDB_6WGN),
                enabled_via_env_var=True,
            )
        )
        assert out_atoms is mm_atoms_in
        assert out_charges is mm_charges_in
        assert meta["fallback_reason"] == "no_cofactor_residues_in_target_card"
