#!/usr/bin/env python
"""
utils/reinsert_cofactors_post_af2.py
------------------------------------
[R-17 G2] Post-AF2 cofactor reinsertion gate.

Problem
~~~~~~~
AlphaFold2 (ColabFold / AlphaFold-Multimer) predicts protein structures from
sequence only and does NOT include HETATM records. For KRAS-style targets
whose QM/MM calculations require structural cofactors (GNP, Mg²⁺) to be
*present* in the MD snapshot (R-17 SSOT), those cofactors must be transplanted
from the reference crystal PDB onto the AF2-predicted complex.

Algorithm
~~~~~~~~~
1. Parse AF2 output PDB → extract Cα coordinates of the target chain.
2. Parse crystal PDB → extract Cα coordinates of the same target chain over
   a shared residue range.
3. Compute the Kabsch rotation + translation that maps crystal Cα → AF2 Cα
   (least-squares superposition).
4. Apply the transform to every HETATM that matches a declared cofactor
   (resname + chain + resnum triplet from ``target_card.cofactor_residues``
   with ``required=True``).
5. Append the transformed HETATM lines to the AF2 PDB, preserving the
   original header / footer records.
6. Validate RMSD of the superposition — if > 2.0 Å, raise
   ``CofactorReinsertionError`` (unsafe transplant).

Dependencies
~~~~~~~~~~~~
Uses only ``numpy`` and the stdlib. Intentionally avoids Bio.PDB so the
module can be unit-tested in the same stub-heavy environment as G6.

Pipeline integration
~~~~~~~~~~~~~~~~~~~~
``reinsert_cofactors(af2_pdb, crystal_pdb, target_card, out_pdb)`` is the
one entry point. Callers (e.g. the orchestrator between ColabFold and
restrained MD) pass the target_card so the resname/chain/resnum triplet
matching is driven by the SSOT.

Python 3.8+ compatible.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

import numpy as np

from utils_common import parse_pdb_atom_line
from cofactor_errors import CofactorReinsertionError


# Maximum acceptable Cα RMSD between crystal and AF2 target chain after
# Kabsch alignment (Å). Above this threshold the transplant is considered
# geometrically unsafe (wrong reference frame / major fold deviation).
DEFAULT_RMSD_THRESHOLD_A = 2.0

# Minimum number of Cα pairs needed for a trustworthy Kabsch alignment.
# Below this, the fit is ill-conditioned (20 chosen to cover partial
# structures while excluding fragments / peptides).
MIN_CA_PAIRS = 20


# ==========================================
# Kabsch alignment (numpy-only)
# ==========================================
def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Least-squares rigid alignment of P onto Q.

    Args:
        P: (N, 3) source coords (moved).
        Q: (N, 3) target coords (reference).

    Returns:
        (rotation_3x3, translation_3, rmsd_scalar). ``P @ R + t ≈ Q``.
    """
    centroid_P = P.mean(axis=0)
    centroid_Q = Q.mean(axis=0)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # SVD on the covariance matrix Pᵀ Q.
    H = P_centered.T @ Q_centered
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = centroid_Q - centroid_P @ R

    P_aligned = P_centered @ R
    rmsd = float(np.sqrt(((P_aligned - Q_centered) ** 2).sum() / max(len(P), 1)))
    return R, t, rmsd


def _extract_ca_by_chain_resnum(pdb_path: str, chain: str) -> Dict[int, np.ndarray]:
    """Return {resnum: (3,) Cα coord} for every Cα in the given chain."""
    result: Dict[int, np.ndarray] = {}
    if not os.path.exists(pdb_path):
        return result
    target_chain = (chain or "").strip().upper()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parsed = parse_pdb_atom_line(line)
            if parsed is None or parsed["record"] != "ATOM":
                continue
            if parsed["name"].strip() != "CA":
                continue
            if (parsed["chain"] or "").strip().upper() != target_chain:
                continue
            result[int(parsed["resnum"])] = np.asarray(
                [parsed["x"], parsed["y"], parsed["z"]], dtype=float
            )
    return result


def _read_hetatm_lines(pdb_path: str) -> List[str]:
    """Return raw HETATM lines from the PDB, preserving file order."""
    out: List[str] = []
    if not os.path.exists(pdb_path):
        return out
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("HETATM"):
                out.append(line)
    return out


# ==========================================
# Public API
# ==========================================
def reinsert_cofactors(
    af2_pdb: str,
    crystal_pdb: str,
    target_card: dict,
    out_pdb: str,
    target_chain_override: Optional[str] = None,
    rmsd_threshold_a: float = DEFAULT_RMSD_THRESHOLD_A,
) -> dict:
    """[R-17 G2] Transplant declared cofactors from crystal PDB into AF2 output.

    Args:
        af2_pdb: Path to AF2-predicted complex PDB (no HETATM).
        crystal_pdb: Path to crystal reference PDB (contains HETATM).
        target_card: Parsed target_card dict with ``cofactor_residues``
            (v0.6.5; required fields assumed present after auto-upgrade).
        out_pdb: Output path for the augmented PDB.
        target_chain_override: If provided, use this chain letter instead of
            ``target_card.target_chain`` (useful when AF2 renamed chains).
        rmsd_threshold_a: Maximum acceptable Cα RMSD after Kabsch.

    Returns:
        Diagnostic dict with keys: ``rmsd_angstrom``, ``n_ca_pairs``,
        ``inserted_cofactors`` (list), ``output_pdb``.

    Raises:
        CofactorReinsertionError: Crystal PDB missing, insufficient Cα pairs,
            RMSD above threshold, or declared cofactor missing from crystal.
    """
    target_chain = target_chain_override or target_card.get("target_chain", "A")

    # Collect declared required cofactor triplets.
    declared: List[dict] = []
    for entry in (target_card.get("cofactor_residues") or []):
        if not isinstance(entry, dict) or not entry.get("required"):
            continue
        declared.append(entry)
    if not declared:
        # No required cofactors → just copy AF2 to out_pdb.
        if not os.path.exists(af2_pdb):
            raise CofactorReinsertionError(
                f"R-17 G2: AF2 PDB not found: {af2_pdb}",
                gate="G2", detail={"af2_pdb": af2_pdb},
            )
        with open(af2_pdb, "r", encoding="utf-8", errors="ignore") as src, \
             open(out_pdb, "w", encoding="utf-8") as dst:
            dst.write(src.read())
        return {
            "rmsd_angstrom": 0.0, "n_ca_pairs": 0,
            "inserted_cofactors": [], "output_pdb": out_pdb,
        }

    if not os.path.exists(crystal_pdb):
        raise CofactorReinsertionError(
            f"R-17 G2: crystal PDB not found: {crystal_pdb}. "
            f"Required to supply declared cofactors {[e.get('resname') for e in declared]}.",
            gate="G2", detail={"crystal_pdb": crystal_pdb},
        )
    if not os.path.exists(af2_pdb):
        raise CofactorReinsertionError(
            f"R-17 G2: AF2 PDB not found: {af2_pdb}.",
            gate="G2", detail={"af2_pdb": af2_pdb},
        )

    # 1. Kabsch alignment over target-chain Cα pairs.
    ca_af2 = _extract_ca_by_chain_resnum(af2_pdb, target_chain)
    ca_crys = _extract_ca_by_chain_resnum(crystal_pdb, target_chain)
    shared_resnums = sorted(set(ca_af2.keys()) & set(ca_crys.keys()))
    if len(shared_resnums) < MIN_CA_PAIRS:
        raise CofactorReinsertionError(
            f"R-17 G2: only {len(shared_resnums)} shared Cα pairs between AF2 "
            f"and crystal on chain '{target_chain}' (need ≥ {MIN_CA_PAIRS}). "
            f"Superposition would be ill-conditioned.",
            gate="G2", n_ca_reference=len(shared_resnums),
            detail={"shared_resnums_sample": shared_resnums[:10]},
        )

    P = np.vstack([ca_crys[r] for r in shared_resnums])
    Q = np.vstack([ca_af2[r] for r in shared_resnums])
    R, t, rmsd = _kabsch(P, Q)

    if rmsd > rmsd_threshold_a:
        raise CofactorReinsertionError(
            f"R-17 G2: Kabsch RMSD {rmsd:.3f} Å exceeds threshold "
            f"{rmsd_threshold_a:.1f} Å over {len(shared_resnums)} Cα pairs "
            f"(target chain '{target_chain}'). AF2 prediction diverges too "
            f"far from crystal — cofactor transplant is unsafe.",
            gate="G2", rmsd_angstrom=rmsd,
            n_ca_reference=len(shared_resnums),
        )

    # 2. Transform declared cofactor HETATM lines.
    crystal_hetatm_lines = _read_hetatm_lines(crystal_pdb)
    declared_keys = set()
    for e in declared:
        key = (
            str(e.get("resname", "")).strip().upper(),
            str(e.get("chain", "")).strip().upper(),
            int(e.get("resnum")) if e.get("resnum") is not None else None,
        )
        declared_keys.add(key)

    inserted_records: List[str] = []
    inserted_meta: List[dict] = []
    found_keys = set()
    for line in crystal_hetatm_lines:
        parsed = parse_pdb_atom_line(line)
        if parsed is None:
            continue
        key = (
            parsed["resname"].strip().upper(),
            (parsed["chain"] or "").strip().upper(),
            int(parsed["resnum"]),
        )
        if key not in declared_keys:
            continue
        coord = np.asarray([parsed["x"], parsed["y"], parsed["z"]], dtype=float)
        new_coord = coord @ R + t  # single-point: matches Kabsch formula.
        new_line = (
            line[:30]
            + "{:8.3f}{:8.3f}{:8.3f}".format(new_coord[0], new_coord[1], new_coord[2])
            + line[54:]
        )
        inserted_records.append(new_line)
        inserted_meta.append({
            "resname": key[0], "chain": key[1], "resnum": key[2],
            "atom_name": parsed["name"],
            "orig": [float(coord[0]), float(coord[1]), float(coord[2])],
            "transformed": [float(new_coord[0]), float(new_coord[1]), float(new_coord[2])],
        })
        found_keys.add(key)

    missing_keys = declared_keys - found_keys
    if missing_keys:
        raise CofactorReinsertionError(
            "R-17 G2: declared cofactor(s) not present in crystal PDB: "
            "{}. Available HETATM resnames: {}.".format(
                sorted(missing_keys),
                sorted({parse_pdb_atom_line(l)["resname"].strip().upper()
                        for l in crystal_hetatm_lines
                        if parse_pdb_atom_line(l) is not None}),
            ),
            gate="G2", detail={"missing_keys": [list(k) for k in sorted(missing_keys)]},
        )

    # 3. Compose output: AF2 records + inserted HETATMs before END.
    with open(af2_pdb, "r", encoding="utf-8", errors="ignore") as f:
        af2_lines = f.readlines()

    # Split on END / ENDMDL so we can slip HETATM records in before the
    # terminator. If no END record is found, append at EOF.
    split_idx = len(af2_lines)
    for i, line in enumerate(af2_lines):
        if line.startswith("END"):
            split_idx = i
            break
    body = af2_lines[:split_idx]
    tail = af2_lines[split_idx:]

    with open(out_pdb, "w", encoding="utf-8") as dst:
        dst.writelines(body)
        if body and not body[-1].endswith("\n"):
            dst.write("\n")
        dst.writelines(inserted_records)
        dst.writelines(tail)
        if not tail:
            dst.write("END\n")

    return {
        "rmsd_angstrom": rmsd,
        "n_ca_pairs": len(shared_resnums),
        "inserted_cofactors": inserted_meta,
        "output_pdb": out_pdb,
    }


__all__ = ["reinsert_cofactors", "DEFAULT_RMSD_THRESHOLD_A", "MIN_CA_PAIRS"]
