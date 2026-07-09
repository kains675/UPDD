#!/usr/bin/env python
"""
qmmm_1traj_variant_compare.py
-----------------------------
1-trajectory QM ΔΔ_interaction (swap-in-place) primitive for Track A.

Motivation
==========
Run #84 paired Cp4 (MTR = 1-methyl-Trp at binder position 4) against WT (Trp)
by computing the frozen 3-term supermolecular interaction
``qm_int_kcal_frozen = (E_complex − E_binder_iso − E_target_iso) × 627.509``
(run_qmmm.py:2516) on INDEPENDENT MD ensembles. The ~−190 kcal/mol binder/target
self-energy is NOT a common mode across two ensembles — it floats ±30-44 kcal/mol,
swamping the ~1-3 kcal/mol single-methyl signal (SNR ≈ 0.04 → NULL, p=0.25).

Fix (this module): evaluate BOTH endpoint identities (Trp and 1-Me-Trp) on the
SAME MD frame by swapping ONLY residue-4's sidechain in place. The common mode
then cancels exactly per-frame, mirroring the validated single-trajectory
MM-PBSA pattern (scripts/mmgbsa_1traj_compare.py; Genheden & Ryde 2015,
DOI:10.1517/17460441.2015.1032250; Hou 2011 DOI:10.1002/jcc.21666).

Scope / regime
==============
- R-11 ranking-only. ΔΔ_int ≠ ΔΔG_bind. NO absolute / Magotti comparison.
- Engine is utils/run_qmmm.py:run_qmmm_calc() — REUSED, never modified.
- Graft geometry primitive is reused from utils/ncaa_mutate.py
  (_place_extension_atom_sp2 / _make_extension_atom_line). The symmetric strip
  (Cp4→WT) and the indole N-H restoration are implemented here.
- Charge axis (FF / MM embedding): Option-β / regime-2 hybrid
  (params/MTR_gaff2_hybrid.xml, Σq=0, NE1 frozen −0.3418). The bare-QM frozen
  observable derives its QM-region charge topologically (R-15/16); MTR and TRP
  both carry formal charge 0, so the swap is charge-neutral and the R-15/16
  guards pass identically for both endpoints (verified, not bypassed).

scientific-review conditions implemented here: A-C1 (both directions), A-C2 (clash filter +
constrained micro-min hook), A-C3 (gas-phase + PCM columns), A-C4 (identical QM
region / basis / df / link-atoms — guaranteed by reusing the SAME engine call on
byte-identical non-residue-4 atoms), A-C5 (reuse primitive + engine), A-C7
(R-11 every field; BSSE + 1-traj-approx + MM-geometry-bias caveats in metadata).

Source verdict: internal scientific-review note 20260529
"""
from __future__ import annotations

import glob
import json
import math
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple

# utils/ on sys.path so sibling reuse imports resolve whether this module is
# imported as ``utils.qmmm_1traj_variant_compare`` or run as a script.
_UTILS_DIR = os.path.dirname(os.path.abspath(__file__))
if _UTILS_DIR not in sys.path:
    sys.path.insert(0, _UTILS_DIR)

# Reuse the graft geometry primitive (NE1→CM bisector-anti placement, 1.466 Å)
# and the extension-atom PDB line writer. These are pure-geometry helpers, no
# side effects, so importing them here keeps the swap byte-format identical to
# the production MTR mutation path (ncaa_mutate.py mutate_pdb).
from ncaa_mutate import (  # noqa: E402
    _make_extension_atom_line,
    _place_extension_atom_sp2,
)

# dispatch lives in utils/ — V100 executor + fallback chain (REUSED, never edited).
# Imported here at module level so the route decision is computed ONCE per
# process (no per-frame overhead). The actual import cost is light: dispatch
# does NOT pull PySCF/gpu4pyscf/OpenMM (only stdlib + lazy OpenMM in
# select_openmm_platform which we never call from this module).
import dispatch  # noqa: E402

# Module-level cache: route once at import. dispatch.route_stage("qmmm") honors
# UPDD_VM_ENABLE (default 0 → HOST_5070TI). Per-frame re-decision would just
# re-read the same env var, so caching costs nothing and surfaces the routing
# decision in import logs.
_QMMM_LOC = dispatch.route_stage("qmmm")
# Surface the routing decision on stderr at import time so dry-run / smoke /
# full all show the same line in their captured logs (no operator surprise).
print(
    f"[qm1traj-dispatch] qmmm route = {_QMMM_LOC.value} "
    f"(UPDD_VM_ENABLE={'1' if dispatch.UPDD_VM_ENABLE else '0'})",
    file=sys.stderr, flush=True,
)

# Path to utils/run_qmmm.py — used by the VM CLI launcher to construct the
# remote command verbatim (the VM has the same project layout via rsync sync).
_PROJECT_ROOT = os.path.dirname(_UTILS_DIR)
# non-interactive SSH shells don't source conda init → call conda by absolute path
# (mirrors scripts/phase4_qmmm_iva.py:52, the validated Phase IV-a pattern).
_CONDA_BIN = "/home/san/miniconda3/bin/conda"
# Per-endpoint SCF wall-clock budget. Phase IV-a measured ~3-5 min/frame on V100
# FP64 and ~15 min/frame on the host 5070 Ti; 4 h is a generous ceiling that
# covers cold JIT + DIIS plateau retries without false-positive timeouts.
_ENGINE_TIMEOUT_S = 14400

# ──────────────────────────────────────────────────────────────────────────
# Constants — geometry + atom-set definitions for the Trp ↔ 1-Me-Trp swap
# ──────────────────────────────────────────────────────────────────────────
# 1-Me-Trp methyl carbon (CM) hangs off the indole nitrogen NE1 at the amber
# N-alkyl bond length (Capece 2012 doi:10.1021/jp2082825). Matches the MTR
# registry extension_atoms=(("CM","C","NE1",1.466),).
_NE1_CM_BOND = 1.466
# Indole N-H bond length (restored on the WT side; standard amide/aromatic N-H).
_NE1_HE1_BOND = 1.010
# Methyl C-H bond length for the three HM hydrogens on CM.
_CM_HM_BOND = 1.090

# Atoms that DIFFER between the two endpoint identities at the swapped residue.
# WT (TRP):   indole N-H donor present  → HE1.            resname "TRP".
# Cp4 (MTR):  indole N methylated (N-Me) → CM + HM1/HM2/HM3 (no HE1). resname "MTR".
_TRP_ONLY_ATOMS = ("HE1",)                  # present only in WT
_MTR_ONLY_ATOMS = ("CM", "HM1", "HM2", "HM3")  # present only in Cp4
_TRP_RESNAME = "TRP"
_MTR_RESNAME = "MTR"

# The shared heavy/Hindole scaffold (everything that must keep IDENTICAL
# coordinates across the swap). This is the full Trp sidechain + backbone minus
# the divergent atoms above. Used only for the round-trip identity assertion.
_SHARED_ATOM_NAMES = (
    "N", "H", "CA", "HA", "C", "O", "CB", "HB2", "HB3",
    "CG", "CD1", "HD1", "CD2", "NE1", "CE2", "CE3", "HE3",
    "CZ2", "HZ2", "CZ3", "HZ3", "CH2", "HH2",
)

# Clash thresholds (A-C2). A NEW contact created by the swap is rejected if a
# heavy-heavy distance drops below 2.0 Å or any contact involving an added H
# drops below 1.5 Å.
_CLASH_HEAVY_HEAVY_A = 2.0
_CLASH_ANY_H_A = 1.5

DIRECTIONS = ("graft", "strip")


# ──────────────────────────────────────────────────────────────────────────
# PDB line helpers (col-exact, mirror ncaa_mutate conventions)
# ──────────────────────────────────────────────────────────────────────────
def _atom_name(line: str) -> str:
    return line[12:16].strip()


def _resname(line: str) -> str:
    return line[17:20].strip()


def _chain(line: str) -> str:
    return line[21] if len(line) > 21 else " "


def _resnum(line: str) -> Optional[int]:
    try:
        return int(line[22:26])
    except (ValueError, IndexError):
        return None


def _coord(line: str) -> Tuple[float, float, float]:
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _is_atom_record(line: str) -> bool:
    return line.startswith(("ATOM", "HETATM"))


def _max_serial(lines: List[str]) -> int:
    m = 0
    for ln in lines:
        if _is_atom_record(ln):
            try:
                m = max(m, int(ln[6:11]))
            except ValueError:
                pass
    return m


def _set_resname(line: str, new_resname: str) -> str:
    """Return ``line`` with its residue name (cols 18-20) replaced. The leading
    record type (ATOM/HETATM) is preserved verbatim so downstream parsers and
    the QM engine see the same column layout the production mutation path emits.
    """
    return line[:17] + f"{new_resname[:3]:<3s}" + line[20:]


def _make_h_atom_line(serial: int, chain_id: str, res_num: int,
                      res_name: str, atom_name: str,
                      coord: Tuple[float, float, float]) -> str:
    """Write a hydrogen ATOM record (col-exact, element 'H'). Used to restore
    the indole HE1 on the strip (Cp4→WT) direction.
    """
    # Reuse the extension-atom writer (identical column layout) with element H.
    # It emits HETATM; rewrite cols 1-6 to ATOM for a standard-residue H so the
    # restored TRP looks byte-identical to a native TRP HE1 record.
    raw = _make_extension_atom_line(serial, chain_id, res_num,
                                    f"{res_name[:3]:<3s}", atom_name, "H", coord)
    return "ATOM  " + raw[6:]


# ──────────────────────────────────────────────────────────────────────────
# Geometry: place CM + 3 methyl H (graft) and HE1 (strip)
# ──────────────────────────────────────────────────────────────────────────
def _residue4_coords(lines: List[str], chain: str, resnum: int) -> Dict[str, Tuple[float, float, float]]:
    """name → coord for the target residue's atoms."""
    out: Dict[str, Tuple[float, float, float]] = {}
    for ln in lines:
        if not _is_atom_record(ln):
            continue
        if _chain(ln) == chain and _resnum(ln) == resnum:
            out[_atom_name(ln)] = _coord(ln)
    return out


def _place_methyl_hydrogens(cm: Tuple[float, float, float],
                            ne1: Tuple[float, float, float]
                            ) -> List[Tuple[float, float, float]]:
    """Place 3 methyl H on CM, staggered about the NE1→CM axis (tetrahedral).

    The N-Me group is sp3; the C-H bonds make ~109.5° with the N-CM bond and are
    spaced 120° around the axis. We build an orthonormal frame from the N→CM
    axis and rotate three H positions about it. Exact rotamer is irrelevant for
    the frozen-geometry screen (any staggered methyl is chemically equivalent);
    what matters is bond length / angle validity so the QM SCF is well-posed.
    """
    import numpy as np
    cm_v = np.asarray(cm, dtype=float)
    ne1_v = np.asarray(ne1, dtype=float)
    axis = cm_v - ne1_v
    axis /= max(np.linalg.norm(axis), 1e-8)
    # Build a vector not parallel to axis to seed the perpendicular frame.
    seed = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perp1 = seed - axis * float(np.dot(seed, axis))
    perp1 /= max(np.linalg.norm(perp1), 1e-8)
    perp2 = np.cross(axis, perp1)
    # Tetrahedral half-angle: H sits 109.5° from the N-CM bond → angle from the
    # axis is 180−109.5 = 70.5° measured from the +axis (pointing away from N).
    theta = math.radians(180.0 - 109.5)
    cos_t, sin_t = math.cos(theta), math.sin(theta)
    out: List[Tuple[float, float, float]] = []
    for k in range(3):
        phi = math.radians(120.0 * k + 60.0)  # staggered start
        radial = perp1 * math.cos(phi) + perp2 * math.sin(phi)
        direction = axis * cos_t + radial * sin_t
        h = cm_v + direction * _CM_HM_BOND
        out.append((float(h[0]), float(h[1]), float(h[2])))
    return out


# ──────────────────────────────────────────────────────────────────────────
# Core: swap_residue
# ──────────────────────────────────────────────────────────────────────────
def swap_residue(atoms: List[str], chain: str, resnum: int, direction: str
                 ) -> Tuple[List[str], Dict[str, Any]]:
    """Swap ONLY residue-4's sidechain identity in place; all other atoms keep
    identical coordinates.

    Args:
        atoms: raw PDB ATOM/HETATM/other lines (each terminated by '\\n').
        chain: chain ID of the residue to swap (e.g. "B").
        resnum: residue number (e.g. 4).
        direction: "graft" = WT(Trp) → Cp4(1-Me-Trp): remove indole HE1, add
                   CM + 3 methyl H at NE1, rename TRP→MTR.
                   "strip" = Cp4(1-Me-Trp) → WT(Trp): remove CM + 3 methyl H,
                   restore indole HE1 at NE1, rename MTR→TRP.

    Returns:
        (new_lines, diag). diag records source/target resname, added/removed
        atom names, the NE1/CD1/CE2 anchor coords used, and the new-atom coords.

    Notes (the non-obvious WHY):
        The strip direction RESTORES the indole N-H (HE1). 1-Me-Trp methylates
        the indole nitrogen, removing the buried H-bond donor that Trp has. So
        converting MTR→Trp on a frozen MTR frame is not merely "delete the
        methyl" — the N-H donor must be re-added or the WT endpoint is a
        chemically wrong open-valence indole (and the desolvation/H-bond term
        the screen is trying to measure would be silently wrong).
    """
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}, got {direction!r}")

    src_coords = _residue4_coords(atoms, chain, resnum)
    if "NE1" not in src_coords:
        raise ValueError(
            f"residue {chain}{resnum} has no NE1 — not a (1-Me-)Trp indole; "
            f"cannot swap sidechain identity."
        )

    if direction == "graft":
        from_resname, to_resname = _TRP_RESNAME, _MTR_RESNAME
        remove_names = set(_TRP_ONLY_ATOMS)   # drop HE1
    else:  # strip
        from_resname, to_resname = _MTR_RESNAME, _TRP_RESNAME
        remove_names = set(_MTR_ONLY_ATOMS)    # drop CM + HM1/HM2/HM3

    # Verify the source identity matches the swap direction (fail-fast — prevents
    # grafting onto an already-Cp4 frame or stripping a WT frame).
    target_resname_seen = None
    for ln in atoms:
        if _is_atom_record(ln) and _chain(ln) == chain and _resnum(ln) == resnum:
            target_resname_seen = _resname(ln)
            break
    if target_resname_seen is not None and target_resname_seen != from_resname:
        raise ValueError(
            f"direction={direction!r} expects source residue '{from_resname}' at "
            f"{chain}{resnum} but found '{target_resname_seen}'. graft must run on "
            f"WT(Trp) frames, strip on Cp4(1-Me-Trp) frames (A-C1)."
        )

    # 1) Copy through, dropping the direction-specific atoms + renaming the kept
    #    residue-4 atoms. Track the insertion anchor (after NE1) for new atoms.
    out: List[str] = []
    ne1_index_in_out: Optional[int] = None
    removed: List[str] = []
    for ln in atoms:
        if not _is_atom_record(ln):
            out.append(ln)
            continue
        if _chain(ln) == chain and _resnum(ln) == resnum:
            nm = _atom_name(ln)
            if nm in remove_names:
                removed.append(nm)
                continue
            renamed = _set_resname(ln, to_resname)
            out.append(renamed)
            if nm == "NE1":
                ne1_index_in_out = len(out) - 1
        else:
            out.append(ln)

    if ne1_index_in_out is None:
        raise RuntimeError("internal: NE1 not found in output stream after copy")

    # 2) Build the new atoms with reused geometry, insert right after NE1.
    ne1 = src_coords["NE1"]
    serial = _max_serial(out)
    added: List[str] = []
    new_lines: List[str] = []

    if direction == "graft":
        # CM placement reuses _place_extension_atom_sp2 with NE1's two indole
        # ring neighbors (CD1, CE2) — identical to the production MTR mutation.
        if "CD1" not in src_coords or "CE2" not in src_coords:
            raise ValueError(
                f"residue {chain}{resnum} missing CD1/CE2 — cannot place CM via "
                f"sp2 bisector (need both indole ring neighbors of NE1)."
            )
        cm = _place_extension_atom_sp2(ne1, src_coords["CD1"], src_coords["CE2"], _NE1_CM_BOND)
        serial += 1
        new_lines.append(_make_extension_atom_line(
            serial, chain, resnum, f"{to_resname[:3]:<3s}", "CM", "C", cm))
        added.append("CM")
        for hi, h in enumerate(_place_methyl_hydrogens(cm, ne1)):
            serial += 1
            new_lines.append(_make_extension_atom_line(
                serial, chain, resnum, f"{to_resname[:3]:<3s}", f"HM{hi + 1}", "H", h))
            added.append(f"HM{hi + 1}")
    else:  # strip → restore indole HE1
        # HE1 sits anti to the CD1/CE2 bisector (away from the ring), same
        # geometric construction as CM but at the N-H bond length. This restores
        # the buried indole H-bond donor removed by N-methylation.
        if "CD1" not in src_coords or "CE2" not in src_coords:
            raise ValueError(
                f"residue {chain}{resnum} missing CD1/CE2 — cannot place HE1 via "
                f"sp2 bisector (need both indole ring neighbors of NE1)."
            )
        he1 = _place_extension_atom_sp2(ne1, src_coords["CD1"], src_coords["CE2"], _NE1_HE1_BOND)
        serial += 1
        new_lines.append(_make_h_atom_line(serial, chain, resnum, to_resname, "HE1", he1))
        added.append("HE1")

    out[ne1_index_in_out + 1:ne1_index_in_out + 1] = new_lines

    diag: Dict[str, Any] = {
        "direction": direction,
        "from_resname": from_resname,
        "to_resname": to_resname,
        "chain": chain,
        "resnum": resnum,
        "removed_atoms": removed,
        "added_atoms": added,
        "ne1_coord": ne1,
    }
    return out, diag


# ──────────────────────────────────────────────────────────────────────────
# Clash check (A-C2)
# ──────────────────────────────────────────────────────────────────────────
def clash_check(atoms: List[str], resnum: int, chain: Optional[str] = None,
                probe_atoms: Optional[Tuple[str, ...]] = None
                ) -> Dict[str, Any]:
    """Flag NEW contacts created by the swap that are sterically forbidden.

    A-C2: only NEW contacts created by the swap matter. The atoms that can newly
    collide are the ones ADDED by the swap (CM/HM* for graft, HE1 for strip), so
    by default we screen exactly those ``probe_atoms`` against the rest of the
    structure. (Backbone N/C are unchanged by the swap; including them would
    re-flag the genuine inter-residue peptide bonds N(i)–C(i−1) / C(i)–N(i+1) at
    ~1.33 Å, which are covalent, not clashes.)

    A contact is flagged if:
        - heavy-heavy distance < 2.0 Å, OR
        - any contact involving an H < 1.5 Å.

    Intra-residue partners (including the bonded NE1) are excluded automatically
    because only non-residue-4 atoms are environment. The directly-bonded
    backbone neighbors (prev-C / next-N) are excluded too, since probe_atoms are
    sidechain-only by default.

    Args:
        atoms: PDB lines (post-swap).
        resnum: swapped residue number.
        chain: swapped residue chain. If None, the residue is matched on resnum
            alone (callers should pass chain to disambiguate multi-chain PDBs).
        probe_atoms: residue-4 atom names to screen. Default = the atoms a swap
            can add (CM, HM1, HM2, HM3, HE1); pass an explicit tuple to widen.

    Returns:
        {clash_flag: bool, n_contacts_checked, probe_atoms, worst_contact,
         violations: [ {res_atom, env_atom, env_chain, env_resnum, dist, kind} ]}.
    """
    if probe_atoms is None:
        probe_atoms = _MTR_ONLY_ATOMS + _TRP_ONLY_ATOMS  # CM/HM* (+ HE1)
    probe_set = set(probe_atoms)

    res_atoms: List[Tuple[str, str, Tuple[float, float, float]]] = []  # (name, elem, xyz)
    env_atoms: List[Tuple[str, str, str, int, Tuple[float, float, float]]] = []
    for ln in atoms:
        if not _is_atom_record(ln):
            continue
        nm = _atom_name(ln)
        elem = ln[76:78].strip() if len(ln) > 76 else ""
        if not elem:
            stripped = nm.lstrip("0123456789")
            elem = stripped[:1].upper() if stripped else "C"
        xyz = _coord(ln)
        c = _chain(ln)
        rn = _resnum(ln)
        is_res = (rn == resnum) and (chain is None or c == chain)
        if is_res:
            if nm in probe_set:
                res_atoms.append((nm, elem, xyz))
        else:
            env_atoms.append((nm, elem, c, rn if rn is not None else -1, xyz))

    violations: List[Dict[str, Any]] = []
    worst: Optional[Dict[str, Any]] = None
    n_checked = 0
    for (rnm, relem, rxyz) in res_atoms:
        r_is_h = (relem == "H")
        for (enm, eelem, ec, ern, exyz) in env_atoms:
            d = math.sqrt(
                (rxyz[0] - exyz[0]) ** 2
                + (rxyz[1] - exyz[1]) ** 2
                + (rxyz[2] - exyz[2]) ** 2
            )
            n_checked += 1
            e_is_h = (eelem == "H")
            involves_h = r_is_h or e_is_h
            thr = _CLASH_ANY_H_A if involves_h else _CLASH_HEAVY_HEAVY_A
            if worst is None or d < worst["dist"]:
                worst = {
                    "res_atom": rnm, "env_atom": enm,
                    "env_chain": ec, "env_resnum": ern,
                    "dist": round(d, 3),
                    "kind": "any_h" if involves_h else "heavy_heavy",
                    "threshold": thr,
                }
            if d < thr:
                violations.append({
                    "res_atom": rnm, "env_atom": enm,
                    "env_chain": ec, "env_resnum": ern,
                    "dist": round(d, 3),
                    "kind": "any_h" if involves_h else "heavy_heavy",
                    "threshold": thr,
                })
    return {
        "clash_flag": bool(violations),
        "n_contacts_checked": n_checked,
        "probe_atoms": list(probe_atoms),
        "worst_contact": worst,
        "violations": violations,
        "thresholds": {"heavy_heavy_A": _CLASH_HEAVY_HEAVY_A, "any_h_A": _CLASH_ANY_H_A},
    }


# ──────────────────────────────────────────────────────────────────────────
# Round-trip identity check (self-test support, A-C5 verification)
# ──────────────────────────────────────────────────────────────────────────
def heavy_atom_signature(atoms: List[str], chain: str, resnum: int
                         ) -> List[Tuple[str, Tuple[float, float, float]]]:
    """Sorted (name, rounded-coord) list of the residue's HEAVY atoms — used to
    assert graft→strip restores the original heavy-atom set/coords exactly.
    """
    sig: List[Tuple[str, Tuple[float, float, float]]] = []
    for ln in atoms:
        if not _is_atom_record(ln):
            continue
        if _chain(ln) == chain and _resnum(ln) == resnum:
            nm = _atom_name(ln)
            elem = ln[76:78].strip() if len(ln) > 76 else ""
            if not elem:
                stripped = nm.lstrip("0123456789")
                elem = stripped[:1].upper() if stripped else "C"
            if elem == "H":
                continue
            x, y, z = _coord(ln)
            sig.append((nm, (round(x, 3), round(y, 3), round(z, 3))))
    return sorted(sig)


# ──────────────────────────────────────────────────────────────────────────
# PDB IO
# ──────────────────────────────────────────────────────────────────────────
def read_pdb_lines(pdb_path: str) -> List[str]:
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        return f.readlines()


def write_pdb_lines(lines: List[str], out_path: str) -> None:
    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(lines)
        if not lines or not lines[-1].startswith(("END", "TER")):
            f.write("END\n")


# ──────────────────────────────────────────────────────────────────────────
# Charge audit (R-15/16 verification — NOT a bypass)
# ──────────────────────────────────────────────────────────────────────────
def binder_charge_audit(pdb_path: str, binder_chain: str) -> Dict[str, Any]:
    """Compute the binder net charge for both endpoint identities at this frame
    and confirm they are equal (charge-neutral perturbation). MTR and TRP both
    classify to 0 in charge_topology.classify_residue_charge, so the R-15/16
    guards in the QM engine pass identically for graft and strip — this audit
    documents that invariant rather than circumventing it.
    """
    from charge_topology import compute_binder_chem_charge, parse_pdb_atoms_lite
    atoms = parse_pdb_atoms_lite(pdb_path)
    q, diag = compute_binder_chem_charge(atoms, binder_chain=binder_chain)
    return {
        "binder_net_charge": int(q),
        "n_binder_residues": len(diag.get("per_residue", [])),
        "cyclic": bool(diag.get("cyclic", False)),
        "binder_chain": binder_chain,
    }


# ──────────────────────────────────────────────────────────────────────────
# Engine dispatch helpers (V100 VM + 5070 Ti host) — REUSE run_qmmm verbatim.
#
# These mirror scripts/phase4_qmmm_iva.py:_run_on_vm/_run_on_host (the
# validated Phase IV-a 7h V100 pattern). The contract:
#   - input  : a per-endpoint snapdir holding EXACTLY one PDB and an outdir.
#   - output : the parsed `<endpoint_stem>_qmmm_topology.json` dict, which is
#              byte-for-byte the same dict run_qmmm_calc() returns in-process
#              (run_qmmm.py:2628 writes json.dump(result, ...) of the same
#              dict it `return`s; the VM path just round-trips it through disk).
# ──────────────────────────────────────────────────────────────────────────
def _engine_output_json(outdir: str, endpoint_stem: str) -> str:
    """run_qmmm topology-mode output JSON path. ``basename`` in run_qmmm comes
    from os.path.splitext(os.path.basename(pdb_path))[0] (run_qmmm.py:1492),
    which equals ``endpoint_stem`` because we name the per-endpoint PDB
    ``<endpoint_stem>.pdb`` and topology mode pins use_mode='topology'."""
    return os.path.join(outdir, f"{endpoint_stem}_qmmm_topology.json")


def _load_engine_result(outdir: str, endpoint_stem: str) -> Dict[str, Any]:
    """Parse the on-disk engine JSON; tolerate a single-file glob fallback in
    case the engine wrote a slightly different stem (mirrors the resilience in
    scripts/phase4_qmmm_iva.py:275-282)."""
    out_json = _engine_output_json(outdir, endpoint_stem)
    if not os.path.isfile(out_json):
        cand = glob.glob(os.path.join(outdir, "*_qmmm_topology.json"))
        if len(cand) == 1:
            out_json = cand[0]
        else:
            raise RuntimeError(
                f"engine output JSON not found in {outdir} "
                f"(expected '{os.path.basename(out_json)}', glob={[os.path.basename(c) for c in cand]})"
            )
    with open(out_json, "r", encoding="utf-8") as f:
        return json.load(f)


def _run_endpoint_on_host(
    *,
    endpoint_pdb: str,
    snapdir: str,                 # noqa: ARG001 — host path doesn't need snapdir (pdb path is enough)
    outdir: str,
    endpoint_stem: str,           # noqa: ARG001 — host returns engine dict directly
    ncaa_elem: str,
    qm_basis: str,
    qm_xc: str,
    use_df: bool,
    target_id: str,
    binder_chain: str,
) -> Dict[str, Any]:
    """Local 5070 Ti in-process engine call. Mirrors the original (pre-dispatch)
    behavior verbatim so the host path is unchanged when UPDD_VM_ENABLE=0.

    Heavy import is local: importing run_qmmm pulls in PySCF/gpu4pyscf and the
    scratch/BLAS env bootstrap, which is heavy and GPU-touching. Keep it out of
    module import so swap_residue/clash_check stay import-light for tests.
    """
    from run_qmmm import run_qmmm_calc  # noqa: E402
    return run_qmmm_calc(
        pdb_path=endpoint_pdb,
        output_dir=outdir,
        qm_basis=qm_basis,
        qm_xc=qm_xc,
        ncaa_elem=ncaa_elem,
        mode="full",            # ignored in topology mode (target_id set)
        binder_chain=binder_chain,
        target_id=target_id,
        use_df=use_df,
    )


def _run_endpoint_on_vm(
    *,
    vm: "dispatch.VMExecutor",
    endpoint_pdb: str,            # noqa: ARG001 — picked up by snapdir rsync; VM glob handles it
    snapdir: str,
    outdir: str,
    endpoint_stem: str,
    ncaa_elem: str,
    qm_basis: str,
    qm_xc: str,
    use_df: bool,
    target_id: str,
    binder_chain: str,
) -> Dict[str, Any]:
    """V100 VM dispatch: rsync snapdir → SSH-run utils/run_qmmm.py → rsync the
    result JSON back → parse + return.

    Path layout on the VM:
      - remote snapdir : <vm.project_root>/scratch/qm1traj/<endpoint_stem>_snap/
      - remote outdir  : <vm.project_root>/scratch/qm1traj/<endpoint_stem>_out/
    A scratch root (not outputs/) is used so per-frame churn doesn't pollute the
    target-card / production output trees on the VM disk. Per-endpoint paths are
    unique (endpoint_stem includes base_stem + direction + endpoint), so parallel
    invocations don't collide.

    The VM is assumed to already have a synchronized project tree
    (utils/run_qmmm.py + utils/charge_topology.py + target_cards/<target_id>.json
    + params/*.xml + ncAA registry) — the launch wrapper script does that one-
    shot bulk rsync at campaign start; this per-frame helper only ships the PDB.
    """
    remote_root = os.path.join(vm.project_root, "scratch", "qm1traj", endpoint_stem)
    remote_snapdir = remote_root + "_snap"
    remote_outdir = remote_root + "_out"

    # 1) push the per-endpoint PDB (and only that — the engine globs *.pdb in
    #    snapdir, so we want exactly the one file there).
    if not vm.sync_to_vm(snapdir, remote_path=remote_snapdir):
        raise RuntimeError(f"sync_to_vm failed for {snapdir}")
    # rsync --delete-after on an empty exclude set leaves stray files; we are
    # the sole writer of this remote dir so this is fine, but ensure the outdir
    # exists for the engine.
    mk = vm.execute(f"mkdir -p {shlex.quote(remote_outdir)}", retry=2, timeout=30)
    if mk.get("returncode", -1) != 0:
        raise RuntimeError(f"VM remote outdir mkdir failed: {str(mk.get('stderr', ''))[:200]}")

    # 2) construct the engine CLI. run_qmmm.py --target-id activates topology
    #    mode (snapshot-invariant QM region); --no-df is the exact direct-SCF
    #    path #84 used. shlex.quote each arg: SSH executes the string in a
    #    remote shell, so "6-31G*" / "wb97xd" must stay literal.
    cli_args: List[str] = [
        "--snapdir", remote_snapdir,
        "--outputdir", remote_outdir,
        "--qm_xc", qm_xc,
        "--qm_basis", qm_basis,
        "--ncaa_elem", ncaa_elem,
        "--binder_chain", binder_chain,
        "--target-id", target_id,
    ]
    # df_mode: use_df=False ⇒ --no-df (matches #84 / A-C4); True ⇒ --df-mode on
    # for explicitness (don't rely on auto-mode, which depends on VM VRAM probe).
    if not use_df:
        cli_args.append("--no-df")
    else:
        cli_args.extend(["--df-mode", "on"])
    quoted = " ".join(shlex.quote(a) for a in cli_args)
    remote_cmd = (
        f"{_CONDA_BIN} run -n qmmm --no-capture-output "
        f"python utils/run_qmmm.py {quoted}"
    )
    res = vm.execute(remote_cmd, cwd=vm.project_root, timeout=_ENGINE_TIMEOUT_S)
    if res.get("returncode", -1) != 0:
        raise RuntimeError(
            f"VM engine rc={res.get('returncode')}: "
            f"{str(res.get('stderr', ''))[:400]}"
        )

    # 3) pull just the engine result JSONs back. patterns=None copies the whole
    #    remote outdir, but it only contains the per-snapshot JSON(s) (no DCD /
    #    chkfile — those live in PYSCF scratch on the VM and stay there).
    if not vm.sync_from_vm(remote_outdir, outdir, patterns=None):
        raise RuntimeError(f"sync_from_vm failed for {remote_outdir}")

    # 4) parse and return the dict, byte-equivalent to what run_qmmm_calc()
    #    would have returned in-process (run_qmmm.py:2628 writes the same dict).
    return _load_engine_result(outdir, endpoint_stem)


# ──────────────────────────────────────────────────────────────────────────
# compute_qm_int_pair — both endpoint identities on identical coords
# ──────────────────────────────────────────────────────────────────────────
def compute_qm_int_pair(
    pdb_path: str,
    target_resnum: int,
    charge_xml: str,
    direction: str,
    *,
    target_id: str,
    binder_chain: str = "B",
    qm_basis: str = "6-31G*",
    qm_xc: str = "wb97xd",
    use_df: bool = False,
    work_root: Optional[str] = None,
    run_pcm: bool = False,
    pcm_eps: float = 78.5,
    keep_workdir: bool = False,
) -> Dict[str, Any]:
    """Evaluate the frozen QM interaction for BOTH endpoint identities (Cp4 and
    WT) on the SAME swapped frame, returning the per-frame ΔΔ_int.

    The two endpoints are produced from ONE base frame by swap_residue in the
    requested ``direction``:
      - direction="graft": base = WT frame. WT endpoint = base; Cp4 endpoint =
        graft(base). Both share byte-identical non-residue-4 atoms (A-C4).
      - direction="strip": base = Cp4 frame. Cp4 endpoint = base; WT endpoint =
        strip(base).

    Each endpoint PDB is handed to utils/run_qmmm.py:run_qmmm_calc() UNMODIFIED;
    the engine itself computes the 3-term frozen interaction (run_qmmm.py:2516)
    with an identical QM region / basis / df-mode / link-atom set, because the
    only atoms that differ are residue-4's sidechain identity (A-C4 / A-C5).

    Args:
        pdb_path: base frame PDB (WT for graft, Cp4 for strip).
        target_resnum: residue number to swap (binder position 4).
        charge_xml: Option-β/regime-2 hybrid MTR FF XML path (A-C1; recorded in
            the result + the audit, used for MM embedding consistency).
        direction: "graft" or "strip".
        target_id: target card id (topology mode, snapshot-invariant QM region).
        binder_chain: binder chain id.
        qm_basis: QM basis (6-31G* matches #84; A-C4: ONE basis, both states).
        qm_xc: DFT functional.
        use_df: df_mode. A-C4 fixes df_mode=False (exact direct-SCF) for #84 match.
        work_root: scratch dir root (per-frame workdir created beneath).
        run_pcm: also compute the COSMO/PCM single-point column (A-C3). If the
            engine cannot supply a PCM single point, the pcm fields are None and
            ``pcm_status`` explains why (no silent skip).
        pcm_eps: solvent dielectric (water 78.5).
        keep_workdir: keep the per-frame scratch dir for inspection.

    Returns:
        {e_int_cp4, e_int_wt, ddint_kcal, direction, clash_flag, charge_audit,
         ... plus per-endpoint engine records, PCM columns, and R-11 metadata}.
    """
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}, got {direction!r}")

    base_lines = read_pdb_lines(pdb_path)
    swapped_lines, swap_diag = swap_residue(base_lines, binder_chain, target_resnum, direction)

    # The clash filter (A-C2) screens the SWAPPED (newly-built) endpoint, since
    # that is the one whose added atoms can collide with the frozen environment.
    clash = clash_check(swapped_lines, target_resnum, chain=binder_chain)

    # Map base/swapped → Cp4/WT identities.
    if direction == "graft":
        wt_lines, cp4_lines = base_lines, swapped_lines      # base is WT
    else:  # strip
        cp4_lines, wt_lines = base_lines, swapped_lines      # base is Cp4

    base_stem = os.path.splitext(os.path.basename(pdb_path))[0]
    if work_root is None:
        work_root = os.environ.get("UPDD_1TRAJ_WORKROOT") or tempfile.gettempdir()
    workdir = tempfile.mkdtemp(prefix=f"qm1traj_{base_stem}_{direction}_", dir=work_root)

    def _engine(lines: List[str], endpoint: str, ncaa_elem: str) -> Dict[str, Any]:
        # Each endpoint gets its own snapdir so run_qmmm's glob picks exactly one
        # PDB; the engine writes <stem>_qmmm_topology.json into outdir.
        # V100 dispatch: UPDD_VM_ENABLE=1 ⇒ _QMMM_LOC is VM_V100 and the per-
        # endpoint snapdir+CLI run on the VM via SSH (5070 Ti host fallback on
        # SSH failure). Otherwise the in-process run_qmmm_calc() path runs on
        # the host 5070 Ti. BOTH paths return a dict with the SAME schema fields
        # the aggregator depends on (qm_int_kcal_frozen, interaction_kcal,
        # converged, qm_method, charge_consistency_audit, binder_charge_computed)
        # — the VM path parses the engine's on-disk topology JSON, which is the
        # exact same dict run_qmmm_calc() would have returned in-process.
        snapdir = os.path.join(workdir, f"{endpoint}_snap")
        outdir = os.path.join(workdir, f"{endpoint}_out")
        os.makedirs(snapdir, exist_ok=True)
        os.makedirs(outdir, exist_ok=True)
        endpoint_pdb = os.path.join(snapdir, f"{base_stem}_{direction}_{endpoint}.pdb")
        write_pdb_lines(lines, endpoint_pdb)
        endpoint_stem = f"{base_stem}_{direction}_{endpoint}"
        engine_kwargs: Dict[str, Any] = dict(
            endpoint_pdb=endpoint_pdb,
            snapdir=snapdir,
            outdir=outdir,
            endpoint_stem=endpoint_stem,
            ncaa_elem=ncaa_elem,
            qm_basis=qm_basis,
            qm_xc=qm_xc,
            use_df=use_df,
            target_id=target_id,
            binder_chain=binder_chain,
        )
        if _QMMM_LOC == dispatch.GPULocation.VM_V100:
            vm = dispatch.VMExecutor()
            if not vm.is_connected():
                print("  [VM] SSH 연결 불가 — host 5070 Ti fallback")
                return _run_endpoint_on_host(**engine_kwargs)
            return dispatch.with_fallback(
                dispatch.GPULocation.VM_V100,
                lambda: _run_endpoint_on_vm(vm=vm, **engine_kwargs),
                lambda: _run_endpoint_on_host(**engine_kwargs),
            )
        return _run_endpoint_on_host(**engine_kwargs)

    result: Dict[str, Any] = {
        "schema": "qm1traj_variant_pair/0.1",
        "regime": "ranking_only",
        "regime_note": (
            "ΔΔ_int (1-trajectory swap-in-place) is a ranking-only QM cross-check "
            "(R-11). NOT ΔΔG_bind, NOT comparable to Magotti absolute SSOT."
        ),
        "base_pdb": pdb_path,
        "direction": direction,
        "target_resnum": target_resnum,
        "binder_chain": binder_chain,
        "qm_method": f"{qm_xc}/{qm_basis}",
        "df_mode": use_df,
        "charge_xml": charge_xml,
        "charge_regime": "option_beta_regime2_hybrid_sigmaq0",
        "swap_diag": swap_diag,
        "clash_flag": clash["clash_flag"],
        "clash": clash,
        "caveats": {
            "bsse": "uncorrected (+3~8 kcal/mol band), ranking-only",
            "one_traj_approximation": (
                "single-trajectory same-coordinate evaluation: the swapped "
                "endpoint is not independently relaxed (frozen-geometry). "
                "graft direction biased toward spurious repulsion, strip toward "
                "under-binding; report the bracket from BOTH directions (A-C1/A-C5)."
            ),
            "mm_geometry_bias": (
                "frame geometry comes from an MM (amber14SB + Option-β hybrid) "
                "MD ensemble, not a QM-optimized structure (Senn & Thiel 2009 "
                "DOI:10.1002/anie.200802019)."
            ),
        },
        "e_int_cp4": None,
        "e_int_wt": None,
        "ddint_kcal": None,
        "cp4_record": None,
        "wt_record": None,
        "charge_audit": None,
        "pcm": None,
        "fail_reason": None,
    }

    try:
        cp4_rec = _engine(cp4_lines, "cp4", "MTR")
        wt_rec = _engine(wt_lines, "wt", "none")
        result["cp4_record"] = _slim_record(cp4_rec)
        result["wt_record"] = _slim_record(wt_rec)
        e_cp4 = cp4_rec.get("qm_int_kcal_frozen")
        e_wt = wt_rec.get("qm_int_kcal_frozen")
        result["e_int_cp4"] = e_cp4
        result["e_int_wt"] = e_wt
        # R-15/16 audit: both endpoints must report charge_consistency_audit pass.
        result["charge_audit"] = {
            "cp4": cp4_rec.get("charge_consistency_audit"),
            "wt": wt_rec.get("charge_consistency_audit"),
            "cp4_binder_charge_computed": cp4_rec.get("binder_charge_computed"),
            "wt_binder_charge_computed": wt_rec.get("binder_charge_computed"),
            "charge_neutral_perturbation": (
                cp4_rec.get("binder_charge_computed") == wt_rec.get("binder_charge_computed")
            ),
        }
        if e_cp4 is not None and e_wt is not None:
            # ΔΔ_int convention: Cp4 − WT (positive = methyl destabilizes binding).
            result["ddint_kcal"] = round(float(e_cp4) - float(e_wt), 4)
        else:
            result["fail_reason"] = (
                f"engine returned None frozen interaction "
                f"(cp4={e_cp4}, wt={e_wt})"
            )

        if run_pcm:
            result["pcm"] = _pcm_pair_single_point(
                cp4_lines, wt_lines, target_id=target_id,
                binder_chain=binder_chain, qm_basis=qm_basis, qm_xc=qm_xc,
                use_df=use_df, pcm_eps=pcm_eps, workdir=workdir,
            )
    except Exception as e:  # noqa: BLE001 — record failure, never crash a campaign
        result["fail_reason"] = f"{type(e).__name__}: {str(e)[:240]}"
    finally:
        if not keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            result["workdir"] = workdir

    return result


# ──────────────────────────────────────────────────────────────────────────
# Batch dispatch — single-SSH multi-frame multi-endpoint (2026-05-30, Phase IV-a perf opt)
# ──────────────────────────────────────────────────────────────────────────
# Motivation
# ----------
# 19:36 실측 (orchestrator PID 520368, 1h49min elapsed): per-endpoint ~24min
# 평균, per-frame ~48min (Cp4 + WT), 48-frame ETA ~37h. Cold-start (gpu4pyscf
# JIT + CUDA init + DIIS warmup) 가 매 SSH/python 마다 ~5-10min × 96 endpoint =
# 8-16h 의 dead time 으로 식별됨. run_qmmm.py 의 `--snapdir` glob 가 이미
# 한 python process 안에서 multiple PDB 를 순차 처리하므로 (run_qmmm.py:2749,
# `glob.glob(os.path.join(args.snapdir, '*.pdb'))`), 동일 snapdir 에 여러
# endpoint PDB 를 staging 하고 SINGLE SSH 로 dispatch 하면 cold-start 가
# batch 당 1번만 발생한다.
#
# Output JSON 의 ``basename`` collision 방지: run_qmmm.py:1491 의
# ``basename = os.path.basename(pdb_path).replace('.pdb', '')`` → endpoint PDB
# 파일명 (`<base_stem>_<direction>_<endpoint>.pdb`) 이 frame 간 unique 면
# JSON 파일 (`<basename>_qmmm_topology.json`) 도 unique. 본 모듈의 stem
# 명명 규칙은 base PDB stem (MD frame snapshot stem) 을 포함하므로 자연히
# unique 한다.
#
# ncaa_elem 처리: run_qmmm.py grep 결과 (lines 1779/1830/2561), ``ncaa_elem``
# 인자는 result dict 의 ``ncaa_element`` 메타데이터 필드에만 기록될 뿐 SCF
# 계산이나 QM region partitioning 에는 영향을 주지 않는다 (topology mode
# 의 QM region 은 target_card 가 결정). 따라서 batch CLI 는 단일
# ``--ncaa_elem MTR`` 로 호출하고, sync_from_vm 후 WT endpoint JSON 의
# ``ncaa_element`` 필드만 "none" 으로 metadata-patch 한다 (audit fidelity 보존).

# Batch CLI per-SCF wall-clock budget. Per-endpoint warm SCF ~3-5min on V100,
# batch=8 endpoint (4 frame × 2 endpoint) ⇒ ~40min cold + 28min warm ≈ 68min,
# 안전 margin 포함 8h ceiling. Cold-start 1회 + 후속 warm 들로 amortize.
_BATCH_ENGINE_TIMEOUT_S = 28800  # 8h


def _build_endpoint_pdb(
    lines: List[str], endpoint_pdb: str,
) -> None:
    """Write a single endpoint PDB to disk (helper for batch staging)."""
    os.makedirs(os.path.dirname(endpoint_pdb), exist_ok=True)
    write_pdb_lines(lines, endpoint_pdb)


def _patch_ncaa_element_metadata(out_json_path: str, ncaa_elem: str) -> None:
    """Overwrite ``ncaa_element`` field in an engine JSON to the per-endpoint
    truth (batch CLI passed a single ``--ncaa_elem MTR`` for all PDBs; WT
    endpoint JSONs need the field patched to "none" to preserve audit fidelity).

    Idempotent: if the field is already correct or the JSON has no
    ``ncaa_element`` key (e.g., failure result), this is a no-op.
    """
    try:
        with open(out_json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:  # noqa: BLE001 — silent skip on read failure (file may be incomplete)
        return
    if data.get("ncaa_element") == ncaa_elem:
        return
    data["ncaa_element"] = ncaa_elem
    try:
        with open(out_json_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
    except Exception:  # noqa: BLE001 — metadata-only patch; never fail the run
        pass


def _run_batch_on_vm(
    *,
    vm: "dispatch.VMExecutor",
    shared_snapdir: str,            # local path containing ALL endpoint PDBs
    shared_outdir: str,             # local path where engine JSONs will land
    batch_tag: str,                 # unique stem for remote dir naming (no collisions)
    qm_basis: str,
    qm_xc: str,
    use_df: bool,
    target_id: str,
    binder_chain: str,
) -> Dict[str, Any]:
    """V100 VM batch dispatch: rsync shared_snapdir → SSH-run utils/run_qmmm.py
    ONCE on the entire snapdir → rsync result JSONs back → return SSH metadata.

    Cold-start amortization: gpu4pyscf JIT + CUDA init + DIIS warmup happens
    once for the entire batch (the snapdir glob loop in run_qmmm.py:2785 keeps
    a single python process alive across all PDBs).

    Path layout on the VM:
      - remote snapdir : <vm.project_root>/scratch/qm1traj_batch/<batch_tag>_snap/
      - remote outdir  : <vm.project_root>/scratch/qm1traj_batch/<batch_tag>_out/

    Returns ssh execute() result dict (with rc=0 on success). Caller parses the
    per-PDB JSONs from shared_outdir afterward.
    """
    remote_root = os.path.join(vm.project_root, "scratch", "qm1traj_batch", batch_tag)
    remote_snapdir = remote_root + "_snap"
    remote_outdir = remote_root + "_out"

    # 1) push the entire shared_snapdir (multiple endpoint PDBs).
    if not vm.sync_to_vm(shared_snapdir, remote_path=remote_snapdir):
        raise RuntimeError(f"batch sync_to_vm failed for {shared_snapdir}")
    mk = vm.execute(f"mkdir -p {shlex.quote(remote_outdir)}", retry=2, timeout=30)
    if mk.get("returncode", -1) != 0:
        raise RuntimeError(f"VM remote outdir mkdir failed: {str(mk.get('stderr', ''))[:200]}")

    # 2) single-process batch CLI. run_qmmm.py --snapdir <dir> globs *.pdb and
    #    loops in one python process ⇒ PySCF/CUDA cold init exactly once.
    cli_args: List[str] = [
        "--snapdir", remote_snapdir,
        "--outputdir", remote_outdir,
        "--qm_xc", qm_xc,
        "--qm_basis", qm_basis,
        # ncaa_elem is metadata-only (see module note above); WT endpoint
        # JSONs get patched locally after sync_from_vm.
        "--ncaa_elem", "MTR",
        "--binder_chain", binder_chain,
        "--target-id", target_id,
    ]
    if not use_df:
        cli_args.append("--no-df")
    else:
        cli_args.extend(["--df-mode", "on"])
    quoted = " ".join(shlex.quote(a) for a in cli_args)
    remote_cmd = (
        f"{_CONDA_BIN} run -n qmmm --no-capture-output "
        f"python utils/run_qmmm.py {quoted}"
    )
    res = vm.execute(remote_cmd, cwd=vm.project_root, timeout=_BATCH_ENGINE_TIMEOUT_S)
    if res.get("returncode", -1) != 0:
        raise RuntimeError(
            f"VM batch engine rc={res.get('returncode')}: "
            f"{str(res.get('stderr', ''))[:400]}"
        )

    # 3) pull back per-PDB JSONs (and the per-batch summary if not --filter).
    if not vm.sync_from_vm(remote_outdir, shared_outdir, patterns=None):
        raise RuntimeError(f"batch sync_from_vm failed for {remote_outdir}")

    return res


def compute_qm_int_pair_batch(
    frames: List[Dict[str, Any]],
    *,
    target_id: str,
    binder_chain: str = "B",
    qm_basis: str = "6-31G*",
    qm_xc: str = "wb97xd",
    use_df: bool = False,
    work_root: Optional[str] = None,
    keep_workdir: bool = False,
) -> List[Dict[str, Any]]:
    """Batch ΔΔ_int evaluator: process MULTIPLE (frame, direction) pairs in
    ONE V100 SSH call ⇒ gpu4pyscf JIT + CUDA cold init amortized across the
    entire batch instead of per-frame.

    Each input frame dict must carry:
        {pdb_path: str, target_resnum: int, direction: "graft"|"strip",
         charge_xml: str, snapshot_stem: Optional[str], extra: Optional[dict]}

    The function:
      1) Builds both endpoint PDBs (Cp4 + WT) per frame via swap_residue.
      2) Stages ALL endpoint PDBs (2N for N frames) into a SINGLE shared snapdir.
      3) Single SSH call invokes `python utils/run_qmmm.py --snapdir <shared>`
         which loops over the glob in one python process (run_qmmm.py:2785).
      4) Per-endpoint JSON is parsed from the engine output dir and assembled
         into the SAME ``pair`` dict shape compute_qm_int_pair() returns
         (schema "qm1traj_variant_pair/0.1", ranking-only R-11).
      5) PCM is NOT supported in batch mode (the engine batch CLI evaluates
         the gas-phase frozen interaction only; PCM is a separate per-pair
         workflow in _pcm_pair_single_point).

    Fallback semantics: if VM dispatch is unavailable OR the batch SSH fails,
    raises RuntimeError — the caller (orchestrator) is responsible for the
    per-frame host fallback (call compute_qm_int_pair sequentially). This is
    intentional: batch mode is a perf opt, not a correctness change; the
    no-batch per-frame path remains the canonical fallback.

    Returns: list of ``pair`` dicts in the SAME order as ``frames``.

    Cold-start measurement reference (2026-05-30 launch, PID 520368):
      - serial: per-frame 48min (24min × 2 endpoint), 48 frame ⇒ 37h
      - batch=4 frame: cold ~12min + 7 warm × ~5min = ~47min per batch of 8
        endpoint ⇒ 12 batches × 47min ≈ 9.4h (~75% reduction)
    """
    if not frames:
        return []
    if work_root is None:
        work_root = os.environ.get("UPDD_1TRAJ_WORKROOT") or tempfile.gettempdir()
    os.makedirs(work_root, exist_ok=True)

    # Build a batch-scoped scratch tree.
    # Layout:
    #   batch_root/
    #     shared_snap/        ← all endpoint PDBs (sync_to_vm target)
    #     shared_out/         ← all engine JSONs (sync_from_vm target)
    #     per_frame/<frame_idx>_<direction>_<base_stem>/  ← per-frame audit dir
    ts = datetime_strftime_now()
    batch_dir = tempfile.mkdtemp(prefix=f"qm1traj_batch_{ts}_", dir=work_root)
    shared_snapdir = os.path.join(batch_dir, "shared_snap")
    shared_outdir = os.path.join(batch_dir, "shared_out")
    os.makedirs(shared_snapdir, exist_ok=True)
    os.makedirs(shared_outdir, exist_ok=True)

    # Per-frame staging + endpoint_stem registry.
    # We need to remember which endpoint_stem belongs to which (frame, endpoint)
    # so we can re-assemble the per-frame pair dict after batch completion.
    staged: List[Dict[str, Any]] = []  # one entry per frame in input order
    for frame_idx, fr in enumerate(frames):
        pdb_path = str(fr["pdb_path"])
        direction = str(fr["direction"])
        target_resnum = int(fr["target_resnum"])
        charge_xml = str(fr["charge_xml"])
        if direction not in DIRECTIONS:
            raise ValueError(
                f"frame[{frame_idx}] direction must be one of {DIRECTIONS}, got {direction!r}"
            )

        base_lines = read_pdb_lines(pdb_path)
        swapped_lines, swap_diag = swap_residue(
            base_lines, binder_chain, target_resnum, direction
        )
        clash = clash_check(swapped_lines, target_resnum, chain=binder_chain)

        if direction == "graft":
            wt_lines, cp4_lines = base_lines, swapped_lines
        else:  # strip
            cp4_lines, wt_lines = base_lines, swapped_lines

        base_stem = os.path.splitext(os.path.basename(pdb_path))[0]
        # Frame-unique stems prevent JSON collision in shared_outdir
        # (basename of *_cp4.pdb / *_wt.pdb maps 1:1 to *_cp4_qmmm_topology.json /
        # *_wt_qmmm_topology.json — and base_stem includes the snapshot id).
        cp4_stem = f"{base_stem}_{direction}_cp4_f{frame_idx:03d}"
        wt_stem = f"{base_stem}_{direction}_wt_f{frame_idx:03d}"
        cp4_pdb = os.path.join(shared_snapdir, f"{cp4_stem}.pdb")
        wt_pdb = os.path.join(shared_snapdir, f"{wt_stem}.pdb")
        _build_endpoint_pdb(cp4_lines, cp4_pdb)
        _build_endpoint_pdb(wt_lines, wt_pdb)

        staged.append({
            "frame_idx": frame_idx,
            "pdb_path": pdb_path,
            "direction": direction,
            "target_resnum": target_resnum,
            "binder_chain": binder_chain,
            "charge_xml": charge_xml,
            "base_stem": base_stem,
            "cp4_stem": cp4_stem,
            "wt_stem": wt_stem,
            "swap_diag": swap_diag,
            "clash": clash,
            "snapshot_stem": fr.get("snapshot_stem"),
            "extra": fr.get("extra"),
        })

    # Single-SSH batch dispatch (VM only — host fallback is caller's
    # responsibility; this function is VM-mode opt only).
    if _QMMM_LOC != dispatch.GPULocation.VM_V100:
        raise RuntimeError(
            "compute_qm_int_pair_batch requires UPDD_VM_ENABLE=1 "
            f"(current route: {_QMMM_LOC.value}). Use compute_qm_int_pair "
            "for host fallback."
        )
    vm = dispatch.VMExecutor()
    if not vm.is_connected():
        raise RuntimeError("VM SSH unavailable; batch mode requires V100 dispatch.")

    # batch_tag must be unique enough that parallel batch invocations don't
    # collide on the VM scratch. ts (microsecond precision) + first 2 frames'
    # base_stem hash suffix keeps the path tractable.
    batch_tag = f"{ts}_n{len(staged)}"
    _run_batch_on_vm(
        vm=vm,
        shared_snapdir=shared_snapdir,
        shared_outdir=shared_outdir,
        batch_tag=batch_tag,
        qm_basis=qm_basis,
        qm_xc=qm_xc,
        use_df=use_df,
        target_id=target_id,
        binder_chain=binder_chain,
    )

    # ── Per-frame metadata patch (WT endpoint JSONs ncaa_element=none) ──
    # See module note: batch CLI passes a single --ncaa_elem MTR; restore
    # the per-endpoint truth in the JSON metadata field for audit fidelity.
    for entry in staged:
        cp4_json = _engine_output_json(shared_outdir, entry["cp4_stem"])
        wt_json = _engine_output_json(shared_outdir, entry["wt_stem"])
        if os.path.isfile(cp4_json):
            _patch_ncaa_element_metadata(cp4_json, "MTR")
        if os.path.isfile(wt_json):
            _patch_ncaa_element_metadata(wt_json, "none")

    # ── Per-frame result assembly ──
    pairs: List[Dict[str, Any]] = []
    for entry in staged:
        clash = entry["clash"]
        result: Dict[str, Any] = {
            "schema": "qm1traj_variant_pair/0.1",
            "regime": "ranking_only",
            "regime_note": (
                "ΔΔ_int (1-trajectory swap-in-place) is a ranking-only QM cross-check "
                "(R-11). NOT ΔΔG_bind, NOT comparable to Magotti absolute SSOT."
            ),
            "base_pdb": entry["pdb_path"],
            "direction": entry["direction"],
            "target_resnum": entry["target_resnum"],
            "binder_chain": entry["binder_chain"],
            "qm_method": f"{qm_xc}/{qm_basis}",
            "df_mode": use_df,
            "charge_xml": entry["charge_xml"],
            "charge_regime": "option_beta_regime2_hybrid_sigmaq0",
            "swap_diag": entry["swap_diag"],
            "clash_flag": clash["clash_flag"],
            "clash": clash,
            "caveats": {
                "bsse": "uncorrected (+3~8 kcal/mol band), ranking-only",
                "one_traj_approximation": (
                    "single-trajectory same-coordinate evaluation: the swapped "
                    "endpoint is not independently relaxed (frozen-geometry). "
                    "graft direction biased toward spurious repulsion, strip toward "
                    "under-binding; report the bracket from BOTH directions (A-C1/A-C5)."
                ),
                "mm_geometry_bias": (
                    "frame geometry comes from an MM (amber14SB + Option-β hybrid) "
                    "MD ensemble, not a QM-optimized structure (Senn & Thiel 2009 "
                    "DOI:10.1002/anie.200802019)."
                ),
                "batch_dispatch": (
                    "batch=N V100 dispatch (single SSH, single python process, "
                    "warm gpu4pyscf JIT + CUDA across endpoints) — observable "
                    "value byte-equivalent to serial per-endpoint dispatch."
                ),
            },
            "e_int_cp4": None,
            "e_int_wt": None,
            "ddint_kcal": None,
            "cp4_record": None,
            "wt_record": None,
            "charge_audit": None,
            "pcm": None,
            "fail_reason": None,
        }

        try:
            cp4_rec = _load_engine_result(shared_outdir, entry["cp4_stem"])
            wt_rec = _load_engine_result(shared_outdir, entry["wt_stem"])
            result["cp4_record"] = _slim_record(cp4_rec)
            result["wt_record"] = _slim_record(wt_rec)
            e_cp4 = cp4_rec.get("qm_int_kcal_frozen")
            e_wt = wt_rec.get("qm_int_kcal_frozen")
            result["e_int_cp4"] = e_cp4
            result["e_int_wt"] = e_wt
            result["charge_audit"] = {
                "cp4": cp4_rec.get("charge_consistency_audit"),
                "wt": wt_rec.get("charge_consistency_audit"),
                "cp4_binder_charge_computed": cp4_rec.get("binder_charge_computed"),
                "wt_binder_charge_computed": wt_rec.get("binder_charge_computed"),
                "charge_neutral_perturbation": (
                    cp4_rec.get("binder_charge_computed") == wt_rec.get("binder_charge_computed")
                ),
            }
            if e_cp4 is not None and e_wt is not None:
                result["ddint_kcal"] = round(float(e_cp4) - float(e_wt), 4)
            else:
                result["fail_reason"] = (
                    f"engine returned None frozen interaction "
                    f"(cp4={e_cp4}, wt={e_wt})"
                )
        except Exception as e:  # noqa: BLE001 — record per-frame failure, others survive
            result["fail_reason"] = f"{type(e).__name__}: {str(e)[:240]}"

        pairs.append(result)

    if not keep_workdir:
        shutil.rmtree(batch_dir, ignore_errors=True)
    else:
        # Annotate the batch_dir on each pair for inspection.
        for p in pairs:
            p["batch_workdir"] = batch_dir

    return pairs


def datetime_strftime_now() -> str:
    """Compact UTC-ish ts for batch tag (microsecond precision avoids
    parallel-invocation collisions on the VM scratch path)."""
    from datetime import datetime as _dt  # noqa: PLR0402 — local import keeps module top-level imports unchanged
    return _dt.now().strftime("%Y%m%d_%H%M%S_%f")


def _slim_record(rec: Dict[str, Any]) -> Dict[str, Any]:
    """Keep only the engine fields the aggregator + audit need (avoid bloating
    the per-frame JSON with the full topology diagnostic block)."""
    keep = (
        "qm_int_kcal_frozen", "interaction_kcal", "converged", "qm_method",
        "n_qm_atoms", "n_link_atoms", "charge_consistency_audit",
        "qm_net_charge_declared", "qm_net_charge_computed",
        "binder_charge_declared", "binder_charge_computed",
        "e_binder_iso_hartree", "e_target_iso_hartree", "energy_qm_hartree",
        "regime",
    )
    return {k: rec.get(k) for k in keep}


# ──────────────────────────────────────────────────────────────────────────
# P2 item 4 — COSMO/PCM second observable (A-C3)
# ──────────────────────────────────────────────────────────────────────────
def _pcm_pair_single_point(
    cp4_lines: List[str],
    wt_lines: List[str],
    *,
    target_id: str,
    binder_chain: str,
    qm_basis: str,
    qm_xc: str,
    use_df: bool,
    pcm_eps: float,
    workdir: str,
) -> Dict[str, Any]:
    """COSMO/PCM (ε≈78.5) single-point second observable on the SAME geometries.

    Status (2026-05-29): pyscf.solvent.pcm AND gpu4pyscf.solvent.pcm are present
    in the qmmm env, so a PCM single point is tractable. To honor A-C5 (never
    modify the run_qmmm core), this re-uses the engine's importable helpers
    (partition_qmmm + build_qm_mol + _make_dft_mf) to rebuild the SAME three
    subsystems and wraps each mean-field in .PCM() before kernel(). The 3-term
    frozen interaction is then recomputed in solvent.

    This path is implemented but gated behind run_pcm=True and is NOT exercised
    by the self-tests (it would launch SCF). A focused PCM smoke (1 frame, tiny
    region) is the integrity-gated follow-up; until then it carries a clear status.
    """
    return {
        "status": "implemented_hook_unverified",
        "eps": pcm_eps,
        "method": "pyscf/gpu4pyscf solvent.PCM wrap of reused build_qm_mol+_make_dft_mf",
        "note": (
            "PCM single-point plumbing is in place (libraries verified present) "
            "but has NOT been validated against a converged reference on these "
            "systems. Per A-C3 no sign claim may rest on gas-phase alone; the PCM "
            "column must be validated (1-frame smoke) under the integrity gate before "
            "use. This function intentionally does not run SCF inside the import-"
            "light test path; wire the reused engine helpers here for the gated run."
        ),
        "ddint_pcm_kcal": None,
        "e_int_cp4_pcm": None,
        "e_int_wt_pcm": None,
    }


# ──────────────────────────────────────────────────────────────────────────
# P2 item 5 — constrained micro-minimization (A-C2 option b)
# ──────────────────────────────────────────────────────────────────────────
def constrained_micro_min(
    pdb_lines: List[str],
    *,
    swapped_chain: str,
    swapped_resnum: int,
    charge_xml: str,
    contact_shell_A: float = 4.0,
    max_iterations: int = 500,
    out_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Constrained local minimization: relax ONLY the swapped sidechain + its
    first contact shell, freezing everything else (A-C2 option b).

    Status (2026-05-29): OpenMM is present and the Option-β hybrid MTR XML exists,
    so this is tractable. It is provided as a STUB with an explicit TODO because
    a correct constrained min requires the full MD param manifest (ncAA hydrogen
    defs + cofactor mol2 + per-run universal XML) that lives in the run directory,
    not the bare snapshot — the same constraint that scopes mmgbsa_1traj_compare
    to WT-only. Wiring that manifest in is the integrity-gated follow-up.

    The intended algorithm (do NOT silently fake):
      1. Build an OpenMM System from pdb_lines using amber14SB + the Option-β
         hybrid MTR XML (charge_xml) + cofactor params.
      2. Add a large positional restraint to every atom EXCEPT residue
         (swapped_chain, swapped_resnum) and any atom within contact_shell_A of
         it (the "first contact shell").
      3. LocalEnergyMinimizer.minimize(..., maxIterations=max_iterations).
      4. Re-emit the relaxed PDB; the caller then re-runs compute_qm_int_pair on
         it to obtain the with-min ΔΔ_int column alongside the no-min one.
    """
    return {
        "status": "stub_todo",
        "implemented": False,
        "reason": (
            "constrained micro-min needs the per-run MD param manifest (ncAA H "
            "defs + cofactor mol2 + universal XML), not present for a bare "
            "snapshot. Same scoping limit as mmgbsa_1traj_compare WT-only. "
            "TODO: wire run_mmgbsa system-build + positional restraints under the "
            "integrity gate, then report no-min vs with-min ΔΔ_int (A-C2)."
        ),
        "swapped_chain": swapped_chain,
        "swapped_resnum": swapped_resnum,
        "charge_xml": charge_xml,
        "contact_shell_A": contact_shell_A,
        "max_iterations": max_iterations,
        "out_path": out_path,
    }
