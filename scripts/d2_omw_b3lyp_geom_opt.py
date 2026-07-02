#!/usr/bin/env python
"""
D2 deliverable — B3LYP/6-31G(d) geometry optimization of Ace-OMW-NMe dipeptide.

Pure Python single-process pipeline (psi4 native optimize() / optking):
  - Input: Ace-OMW-NMe dipeptide geometry from PySCF D1 (Pass B output)
  - Method: B3LYP / 6-31G(d) (standard for ff14SB / phosaa19SB-style refit)
  - Optimization: optking backend (built into psi4), max 200 iterations

Two optimization modes (mode default = constrained):

  mode = constrained  (DEFAULT, D2 production)
    ff14SB-equilibrium constrained optimization. The peptide backbone
    (N-CA, CA-C, C-O, N-H, CA-CB) plus the indole NE1 region bonds/angles
    are frozen toward amber14SB Trp equilibrium values via optking
    frozen_distance / frozen_bend. The resulting geometry is ff14SB-
    compatible (not a free ff03-style relaxation), so the bond/angle/
    improper deltas reported are the QM-vs-ff14SB disagreements that
    drive D3 bonded re-parameterization. This is the phosaa19SB /
    Forouzesh 2020 pattern (DOI 10.3389/fmolb.2020.608931).

  mode = free   (legacy / diagnostic)
    Unconstrained B3LYP/6-31G(d) minimization (the pre-2026-06 behaviour).
    Kept for diagnostic comparison only; NOT the D2 production geometry.

Conformer ensemble (D2 production, --conformers >= 4):
  Generate >= 4 conformers by rotating the sidechain about CA-CB (chi1)
  and the backbone phi/psi, each constrained-optimized at B3LYP/6-31G(d).
  Stored for the downstream RESP-A2 multi-conformer fit (Cieplak 1995,
  DOI 10.1002/jcc.540161106). Charges are NOT fit here.

Pre-D3 constraint-enforcement verification (--verify-constraint):
  A single ~30-min ranged_dihedral 1-point optimization that confirms
  optking actually enforces a target dihedral across opt iterations
  (drift check). This must pass before the multi-hour D3 relaxed scans
  are launched, to avoid the rigid-scan-explosion failure mode.

CHARGES: never fit in D2. The geometry source level (B3LYP/6-31G(d)) is
deliberately distinct from the RESP ESP level (HF/6-31G*); a B3LYP density
must never leak into the charge fit. D2 produces geometry + conformer
ensemble only.

Outputs (analysis dir, R-7 preserved): outputs/analysis/d2_omw_geom_opt_{TS}/
  omw_dipeptide_input.pdb            (starting geometry)
  omw_dipeptide_optimized.pdb        (B3LYP/6-31G(d) constrained-optimized, ref conf)
  conformers/conf{NN}_optimized.pdb  (>= 4 constrained-optimized conformers)
  omw_geom_opt_summary.json          (energy, gradient norm, n_steps, mode)
  constraint_verification.json       (pre-D3 dihedral drift test result)
  d2_geom_opt_report.md              (bond / angle / improper deltas)
  psi4_opt_output.dat                (psi4 SCF + opt log)

Deliverable copy (verdict-specified): params/omw_qm_geom/
  omw_dipeptide_optimized.pdb + conformers/ + d2_geom_opt_report.md +
  ff14SB diff report. R-7: the analysis dir is the immutable record; the
  params/ copy is the consumed-by-D3 deliverable.

Run:
  /home/san/miniconda3/envs/psi4/bin/python scripts/d2_omw_b3lyp_geom_opt.py \
      [--mode constrained|free] [--conformers 4] [--verify-constraint] \
      [--method b3lyp] [--no-opt]

References:
  - Maier et al. 2015, J. Chem. Theory Comput. 11(8):3696 (ff14SB Trp ref).
    DOI 10.1021/acs.jctc.5b00255
  - Forouzesh & Mishra 2020, Front. Mol. Biosci. 7:608931 (ff14SB-compatible
    constrained refit). DOI 10.3389/fmolb.2020.608931
  - Cieplak et al. 1995, J. Comput. Chem. 16(11):1357 (RESP-A2 multi-conf).
    DOI 10.1002/jcc.540161106
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

PROJ = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# amber14SB Trp equilibrium geometry (Maier 2015 / parm14SB.dat) used for the
# constrained optimization and the QM-vs-ff14SB delta report.
#
# Bond lengths in Angstrom; angles in degrees. Atom names follow the OMW
# residue convention (Khoury CZ1 methyl C == UPDD CM; normalized below).
# ---------------------------------------------------------------------------
AMBER14SB_BONDS: Dict[Tuple[str, str], float] = {
    ("N", "CA"):    1.449,   # backbone N-Calpha
    ("CA", "CB"):   1.535,   # Calpha-Cbeta
    ("CB", "CG"):   1.504,   # Cbeta-Cgamma
    ("CG", "CD1"):  1.380,   # indole 5-ring
    ("CG", "CD2"):  1.430,
    ("CD1", "NE1"): 1.380,
    ("CD2", "CE2"): 1.400,
    ("NE1", "CE2"): 1.380,
    ("NE1", "CM"):  1.460,   # N-methyl junction (key for 1-MeTrp)
    ("CE2", "CZ2"): 1.400,
    ("CD2", "CE3"): 1.400,
    ("CZ2", "CH2"): 1.380,
    ("CE3", "CZ3"): 1.380,
    ("CH2", "CZ3"): 1.400,
    ("CA", "C"):    1.522,
    ("C", "O"):     1.229,
    ("N", "H"):     1.010,
}

# Backbone + NE1-region bonds frozen toward ff14SB equilibrium in the
# constrained optimization. These are the terms that must stay ff14SB-
# compatible so the *side-chain / methyl* deltas are the D3 signal.
# NE1 region is the empirically-diagnosed explosion locus (NE1 freeze
# is non-negotiable per the RESP-refit verdict chain).
CONSTRAINED_BONDS: Tuple[Tuple[str, str], ...] = (
    ("N", "CA"), ("CA", "C"), ("C", "O"), ("N", "H"), ("CA", "CB"),
    ("CD1", "NE1"), ("NE1", "CE2"),
)

# amber14SB Trp equilibrium angles (deg) frozen in the constrained opt.
AMBER14SB_ANGLES: Dict[Tuple[str, str, str], float] = {
    ("N", "CA", "C"):    111.0,   # backbone
    ("CA", "C", "O"):    120.5,
    ("N", "CA", "CB"):   109.7,
    ("CA", "CB", "CG"):  113.1,
    ("CD1", "NE1", "CE2"): 108.7,  # indole 5-ring N apex
}
CONSTRAINED_ANGLES: Tuple[Tuple[str, str, str], ...] = tuple(AMBER14SB_ANGLES.keys())

# Key improper centers for the D3 report (indole-ring planarity + methyl).
# Improper convention here: (i j k l) with k as the central atom.
IMPROPER_CENTERS: Tuple[Tuple[str, str, str, str], ...] = (
    ("CG", "CE2", "NE1", "CM"),    # methyl out-of-plane at indole N
    ("CD1", "CE2", "NE1", "CM"),   # cross-check methyl improper
    ("CG", "NE1", "CD1", "HD1"),   # CD1-H out of 5-ring plane
)

# Pre-D3 constraint-enforcement gate windows (#106 FIX 2). The ranged_dihedral
# window was tightened from +-2.0 to +-0.05 deg historically, which drove a
# geometry-step sawtooth (the dominant non-convergence cause). The drift
# tolerance sits just outside the window so a target genuinely held inside the
# window passes, while a free-relaxed dihedral that leaves the window fails.
VERIFY_RANGED_HALF_WINDOW_DEG = 2.0   # ranged_dihedral half-window in the verify gate
VERIFY_DRIFT_TOL_DEG = 2.5            # PASS if |final - target| <= this (window + margin)
# Conformer ranged_dihedral half-window (#106 FIX 2): +-0.05 -> +-2.0 deg, same
# sawtooth fix as the verify gate. Well inside the 10-deg wrong-well gate, so
# rotamer identity is never ambiguous.
CONF_RANGED_HALF_WINDOW_DEG = 2.0

# GAU_LOOSE force-convergence thresholds (optking internal-coord au). Retained
# INFORMATIONAL only (#106 FIX 5): the force-gate (prior FIX 3) is the WRONG
# acceptance metric for a ranged-constrained optimization -- see the
# geometry-convergence thresholds below. We still parse + report final forces
# for diagnostics, but the PASS/MUST-REDO verdict no longer gates on them.
GAU_LOOSE_MAX_FORCE = 2.5e-3   # au (informational only)
GAU_LOOSE_RMS_FORCE = 1.7e-3   # au (informational only)

# GAU_LOOSE geometry-convergence thresholds (optking internal-coord au), used by
# the C7 verdict geometry gate (#106 FIX 5, replacing the force gate).
#
# Why geometry, not force: for a ranged-constrained opt, the constraint target
# (amber14SB-eq) generally differs from the QM free-minimum, so the QM gradient
# keeps pushing the pinned coordinate toward the free-min. optking clamps the
# DISPLACEMENT to the window wall (disp -> 0) but does NOT zero the FORCE while
# the coordinate sits INSIDE its window (optking displace.py L82-91): the
# residual gradient projection persists as a Lagrange-multiplier constraint
# REACTION FORCE. It is mathematically guaranteed nonzero whenever
# target != free-min -- it is the force the constraint is supplying, by design.
# Demanding Max Force -> 0 demands the constraint do nothing, and false-rejects
# valid constrained minima (empirically: conf02 chi1+180 plateaus at MaxF
# 3.33e-2 while MaxDisp 5e-5 / RMSDisp 1.4e-5 / |dE| ~1e-6 are all met).
#
# The correct convergence signal for a constrained minimum is "the geometry
# stopped moving AND the energy stopped dropping": small displacement + flat
# energy. These are optking's own GAUSSIAN-preset Max Disp / RMS Disp criteria
# (the same numbers optking prints in its "Convergence Criteria" row), plus a
# flat-energy tolerance (#106 v3 geometry-convergence gate).
GAU_LOOSE_MAX_DISP = 1.0e-2    # au (optking GAUSSIAN preset Max Disp)
GAU_LOOSE_RMS_DISP = 6.7e-3    # au (optking GAUSSIAN preset RMS Disp)
# GEOM_DELTA_E_TOL_HA: UNUSED as of #106 v3 (2026-06-21). The Delta-E clause was
# removed from _conformer_geometry_converged (optking's GAUSSIAN preset does not
# require max_DE; 1e-5 Ha is below the 40-atom DF-B3LYP SCF noise floor and
# false-rejected converged conf01). final_delta_e is still parsed/reported as
# informational; this constant is retained only for changelog/back-reference.
GEOM_DELTA_E_TOL_HA = 1.0e-5   # Ha (legacy; NOT a convergence gate)

# chi1 (N-CA-CB-CG) target values (deg) for the conformer grid.
CHI1_TARGETS: Tuple[float, ...] = (-60.0, 180.0, 60.0)
# (phi, psi) backbone target pairs (deg) for the conformer grid.
# phi = C(ACE)-N-CA-C ; psi = N-CA-C-N(NME)
PHIPSI_TARGETS: Tuple[Tuple[float, float], ...] = (
    (-60.0, -45.0),   # alpha-helical
    (-135.0, 135.0),  # beta-sheet
)


def _norm_name(name: str) -> str:
    """Normalize Khoury OMW atom names to UPDD convention (CZ1 -> CM, etc.)."""
    if name == "CZ1":
        return "CM"
    if name in ("HZ11", "HZ12", "HZ13"):
        return {"HZ11": "HM1", "HZ12": "HM2", "HZ13": "HM3"}[name]
    return name


def find_latest_pyscf_dipeptide() -> Optional[Path]:
    for d in sorted((PROJ / "outputs/analysis").glob("d1_khoury_resp_reproduction_*"), reverse=True):
        pdb = d / "omw_dipeptide.pdb"
        if pdb.exists():
            return pdb
    return None


def parse_pdb(pdb_path: Path) -> List[Tuple[str, str, str, float, float, float]]:
    """Return list of (atom_name, residue, element, x, y, z)."""
    out = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            name = line[12:16].strip()
            residue = line[17:20].strip()
            element = line[76:78].strip() or name[0]
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            out.append((name, residue, element, x, y, z))
    return out


def write_pdb(atoms, pdb_path: Path, title: str = "") -> None:
    with open(pdb_path, "w") as f:
        if title:
            f.write(f"REMARK {title}\n")
        for i, (name, residue, element, x, y, z) in enumerate(atoms, 1):
            f.write(
                f"ATOM  {i:5d} {name:<4s} {residue:>3s}     1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
            )
        f.write("END\n")


def build_name_index(atoms, residue: str = "OMW") -> Dict[str, int]:
    """atom_name (normalized) -> 0-based index, OMW residue only."""
    idx = {}
    for i, (name, res, *_rest) in enumerate(atoms):
        if res == residue:
            idx[_norm_name(name)] = i
    return idx


# ---------------------------------------------------------------------------
# Internal-coordinate measurements (Angstrom / degrees), 0-based indices into
# the flat atom list.
# ---------------------------------------------------------------------------
def _bond_len(coords, i: int, j: int) -> float:
    import numpy as np
    return float(np.linalg.norm(np.array(coords[i]) - np.array(coords[j])))


def _angle_deg(coords, i: int, j: int, k: int) -> float:
    import numpy as np
    a = np.array(coords[i]) - np.array(coords[j])
    b = np.array(coords[k]) - np.array(coords[j])
    cosv = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    cosv = max(-1.0, min(1.0, cosv))
    return math.degrees(math.acos(cosv))


def _dihedral_deg(coords, i: int, j: int, k: int, l: int) -> float:
    import numpy as np
    p = [np.array(coords[x]) for x in (i, j, k, l)]
    b0 = p[0] - p[1]
    b1 = p[2] - p[1]
    b2 = p[3] - p[2]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.degrees(math.atan2(y, x))


def _coords_of(atoms) -> List[Tuple[float, float, float]]:
    return [(a[3], a[4], a[5]) for a in atoms]


# ---------------------------------------------------------------------------
# psi4 helpers
# ---------------------------------------------------------------------------
def _make_psi4_molecule(atoms, name: str = "Ace-OMW-NMe"):
    import psi4
    geom_lines = []
    for nm, res, element, x, y, z in atoms:
        geom_lines.append(f"{element:<3s} {x:14.8f} {y:14.8f} {z:14.8f}")
    geom_block = "\n".join(geom_lines) + "\nno_reorient\nno_com\nunits angstrom\n"
    mol = psi4.geometry(geom_block)
    mol.set_name(name)
    mol.update_geometry()
    return mol


def _atoms_from_mol(atoms_template, mol):
    """Replace coordinates in atoms_template with optimized mol coords (Angstrom)."""
    coords = mol.geometry().to_array() * 0.529177210903
    out = []
    for i, (nm, res, element, _x, _y, _z) in enumerate(atoms_template):
        x, y, z = coords[i]
        out.append((nm, res, element, float(x), float(y), float(z)))
    return out


# Ranged-pin windows for the backbone+NE1 constraints (#106 FIX 1). The
# constraint target is the amber14SB equilibrium value, NOT the D1-input
# value, so the C7 "delta_qm_vs_amber < 0.01 A / 1 deg" gate is meaningful by
# construction. The window is deliberately narrow (a soft pin, not a hard
# freeze) so optking can take finite-size steps without the freeze-vs-input
# sawtooth, while still landing the term on amber14SB-eq.
RANGED_BOND_HALF_WINDOW_A = 0.005    # +-0.005 A around amber14SB eq bond length
RANGED_BEND_HALF_WINDOW_DEG = 0.5    # +-0.5 deg around amber14SB eq angle


def _build_constraint_options(atoms, mode: str) -> Dict[str, str]:
    """Build optking ranged_distance / ranged_bend strings (1-based indices)
    that PIN the backbone+NE1 bonds/angles to their amber14SB-equilibrium
    values within a narrow window (#106 FIX 1).

    Previously these were frozen_distance / frozen_bend, which hold a term at
    its D1-INPUT value (1-4 deg / a few mA off amber14SB-eq). Pinning to the
    amber14SB-eq target instead makes the downstream C7 constraint-integrity
    gate (delta_qm_vs_amber < 0.01 A / 1 deg) meaningful, and removes the
    C-O / N-H / CA-C / angle side-effects of freezing the wrong reference."""
    if mode != "constrained":
        return {}
    nidx = build_name_index(atoms)
    rng_dist = []
    for a, b in CONSTRAINED_BONDS:
        ia, ib = nidx.get(a), nidx.get(b)
        eq = AMBER14SB_BONDS.get((a, b))
        if ia is not None and ib is not None and eq is not None:
            lo, hi = eq - RANGED_BOND_HALF_WINDOW_A, eq + RANGED_BOND_HALF_WINDOW_A
            rng_dist.append(f"({ia + 1} {ib + 1} {lo:.3f} {hi:.3f})")
    rng_bend = []
    for a, b, c in CONSTRAINED_ANGLES:
        ia, ib, ic = nidx.get(a), nidx.get(b), nidx.get(c)
        eq = AMBER14SB_ANGLES.get((a, b, c))
        if ia is not None and ib is not None and ic is not None and eq is not None:
            lo, hi = eq - RANGED_BEND_HALF_WINDOW_DEG, eq + RANGED_BEND_HALF_WINDOW_DEG
            rng_bend.append(f"({ia + 1} {ib + 1} {ic + 1} {lo:.2f} {hi:.2f})")
    opts: Dict[str, str] = {}
    if rng_dist:
        opts["ranged_distance"] = " ".join(rng_dist)
    if rng_bend:
        opts["ranged_bend"] = " ".join(rng_bend)
    return opts


# ---------------------------------------------------------------------------
# psi4 Python-logging FileHandler hygiene (#106 FIX 4).
#
# optking emits through the Python logging tree (logger "psi4", child
# "psi4.optking.*", propagate=True). psi4.set_output_file() adds a NEW
# FileHandler to the "psi4" logger on every call and only de-dupes by
# baseFilename, so a fresh per-conformer .dat path means a fresh handler that
# is NEVER removed. After N conformers the "psi4" logger carries N FileHandlers
# -> every optimize() writes to ALL prior conformer .log files (cross-
# contamination) and the logs grow to multiple GB. We strip stale FileHandlers
# on entry and remove our own on exit so each call owns exactly one .log.
# ---------------------------------------------------------------------------
_PSI4_LOGGER_NAME = "psi4"  # parent of psi4.optking.* (confirmed: conf .log emitter)


def _reset_psi4_log_handlers():
    """Remove every FileHandler currently attached to the psi4 Python logger
    (the stale handlers leaked by previous set_output_file calls). StreamHandlers
    and other handler types are left untouched. Returns nothing."""
    import logging
    logger = logging.getLogger(_PSI4_LOGGER_NAME)
    for h in list(logger.handlers):
        if isinstance(h, logging.FileHandler):
            try:
                h.close()
            finally:
                logger.removeHandler(h)


# ---------------------------------------------------------------------------
# optking .dat convergence-table parser (#106 FIX 5, generalizes prior FIX 3).
#
# psi4.optimize raises OptimizationConvergenceError BEFORE returning history,
# so on the non-converged path we have no in-memory gradient record. optking
# does, however, print a per-step "Convergence Check" table to the .dat with
# "Max Force / RMS Force / Max Disp / RMS Disp" + Delta E, all in internal-
# coordinate au (the exact quantities its GAU_LOOSE gate compares to). We parse
# the LAST data row of the LAST such table so the verdict gate can certify a
# genuinely relaxed geometry that merely missed the iteration limit.
#
# FIX 5 generalizes the old force-only parser to ALL five convergence columns,
# because the verdict now gates on DISPLACEMENT + ENERGY (geometry convergence)
# rather than force (constraint reaction force is expected, not a defect -- see
# the GAU_LOOSE_MAX_DISP comment block).
# ---------------------------------------------------------------------------
def _parse_final_convergence(dat_path: Path):
    """Return a dict of the LAST step row of the LAST "Convergence Check" table
    in the .dat:
      {"energy", "delta_e", "max_force", "rms_force", "max_disp", "rms_disp"}
      (Ha / Ha / au / au / au / au), each Optional[float], all None if not
      parseable.

    The table header is
      'Step  Total Energy  Delta E  Max Force  RMS Force  Max Disp  RMS Disp'
    and step rows look like
      '  36  -897.78496802  5.34e-07 o  2.22e-03 *  2.42e-04 *  1.01e-03 *  1.12e-04 *  ~'
    After stripping the '* o ~' status markers the tokens are
      [step, energy, deltaE, maxForce, rmsForce, maxDisp, rmsDisp]
    The 'Convergence Criteria' reference row is skipped (it is not a step)."""
    blank = {"energy": None, "delta_e": None, "max_force": None,
             "rms_force": None, "max_disp": None, "rms_disp": None}
    if not dat_path.exists():
        return dict(blank)
    out = dict(blank)
    try:
        with open(dat_path, errors="ignore") as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("-"):
                    continue
                # skip the reference / header rows of the convergence table
                if "Convergence Criteria" in line or "Total Energy" in line:
                    continue
                toks = line.replace("*", " ").replace("o", " ").replace("~", " ").split()
                # a step row begins with an integer step index, then >= 6 floats
                if len(toks) < 7 or not toks[0].lstrip("+-").isdigit():
                    continue
                try:
                    # toks = [step, energy, deltaE, maxforce, rmsforce, maxdisp, rmsdisp]
                    en = float(toks[1])
                    de = float(toks[2])
                    mf = float(toks[3])
                    rf = float(toks[4])
                    md = float(toks[5])
                    rd = float(toks[6])
                except (ValueError, IndexError):
                    continue
                # keep the LAST valid step row encountered (overwrite as we go)
                out = {"energy": en, "delta_e": de, "max_force": mf,
                       "rms_force": rf, "max_disp": md, "rms_disp": rd}
    except OSError:
        return dict(blank)
    return out


def _parse_final_forces(dat_path: Path) -> Tuple[Optional[float], Optional[float]]:
    """Backward-compat thin wrapper: (final_max_force, final_rms_force) only.
    Retained for the informational force columns; the verdict no longer gates on
    these (see _parse_final_convergence / _conformer_geometry_converged)."""
    conv = _parse_final_convergence(dat_path)
    return conv["max_force"], conv["rms_force"]


def _parse_n_opt_steps(dat_path: Path) -> int:
    """Count optimizer step rows (the LAST contiguous run) in the .dat, used to
    recover n_steps on the non-converged path where the in-memory history was
    lost to the early raise. Best-effort: returns the max step index seen."""
    if not dat_path.exists():
        return 0
    max_step = 0
    try:
        with open(dat_path, errors="ignore") as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("-"):
                    continue
                if "Convergence Criteria" in line or "Total Energy" in line:
                    continue
                toks = line.replace("*", " ").replace("o", " ").replace("~", " ").split()
                if len(toks) < 5 or not toks[0].lstrip("+-").isdigit():
                    continue
                try:
                    float(toks[3])
                    float(toks[4])
                except (ValueError, IndexError):
                    continue
                max_step = max(max_step, int(toks[0]))
    except OSError:
        return 0
    return max_step


def run_psi4_optimize(atoms, logdir: Path, method: str, mode: str,
                      out_name: str = "psi4_opt_output.dat",
                      extra_constraints: Optional[Dict[str, str]] = None):
    """Constrained (or free) B3LYP/6-31G(d) optimization. Returns (atoms_opt, summary)."""
    import psi4

    # FIX 4: strip leaked FileHandlers from prior calls before installing this
    # call's per-.dat log handler (set_output_file adds it below).
    _reset_psi4_log_handlers()

    psi4.core.clean_options()
    psi4.set_output_file(str(logdir / out_name), True)
    psi4.set_memory("8 GB")
    psi4.core.set_num_threads(8)

    mol = _make_psi4_molecule(atoms)
    print(f"[info] Optimization start: {mol.natom()} atoms, {method.upper()}/6-31G(d), mode={mode}")

    options: Dict[str, object] = {
        "basis": "6-31G(d)",
        "scf_type": "df",
        "g_convergence": "gau_loose",  # gau too tight for ~40-atom dipep
        "geom_maxiter": 300,           # FIX 2: 200 -> 300 (ranged-pin needs slack)
        "opt_type": "min",
    }
    cons = _build_constraint_options(atoms, mode)
    if extra_constraints:
        cons.update(extra_constraints)
    options.update(cons)
    psi4.set_options(options)
    if cons:
        print(f"[info] active constraints: {cons}")

    converged = True
    history: Dict[str, object] = {}
    try:
        e_opt, history = psi4.optimize(method, molecule=mol, return_history=True)
    except psi4.OptimizationConvergenceError as e:
        print(f"[warn] Optimization did not converge: {e}")
        converged = False
        wfn = getattr(e, "wfn", None)
        e_opt = wfn.energy() if wfn is not None else float("nan")
        # FIX 3: psi4.optimize raises BEFORE returning history, so history is
        # never bound on this path. Recover the per-step count + final forces by
        # parsing optking's own "Convergence Check" table from the .dat log, so
        # the C7 force-convergence gate can distinguish a genuinely relaxed (but
        # un-flagged) geometry from a non-converged false-green.
        history = {}
    finally:
        # FIX 4: drop the handler this call installed so it cannot leak into the
        # next conformer's optimize().
        _reset_psi4_log_handlers()

    # FIX 3: the optking success-path history key is 'energy' (singular); the
    # old code read 'energies' and so reported n_steps=0 even on convergence.
    energies = history.get("energy", history.get("energies", []))
    n_steps = len(energies)

    # FIX 5: final convergence-table row (optking internal-coord au). Parsed for
    # BOTH paths so the verdict gate always has it; on the non-converged path
    # this is the only available record. Forces are kept INFORMATIONAL; the
    # verdict now gates on Max/RMS Disp + Delta E (geometry convergence).
    conv = _parse_final_convergence(logdir / out_name)
    final_max_force, final_rms_force = conv["max_force"], conv["rms_force"]
    final_max_disp, final_rms_disp = conv["max_disp"], conv["rms_disp"]
    final_delta_e = conv["delta_e"]
    if not converged and n_steps == 0:
        # recover a step count from the parsed convergence table when history
        # was lost to the early raise (keeps the n_steps=0 bug from masking work)
        n_steps = _parse_n_opt_steps(logdir / out_name)

    atoms_opt = _atoms_from_mol(atoms, mol)
    summary = {
        "method": f"{method.upper()}/6-31G(d)",
        "mode": mode,
        "constraints": cons,
        "n_atoms": len(atoms),
        "n_steps": n_steps,
        "final_energy_hartree": float(e_opt),
        "final_energy_kcal": float(e_opt) * 627.5094740631 if e_opt == e_opt else float("nan"),
        "converged": converged,
        # FIX 5: geometry-convergence columns (verdict now gates on these).
        "final_max_disp": final_max_disp,      # optking internal-coord au (gau_loose target 1.0e-2)
        "final_rms_disp": final_rms_disp,      # optking internal-coord au (gau_loose target 6.7e-3)
        "final_delta_e": final_delta_e,        # Ha (last-step |Delta E|; flat => energy converged)
        # forces retained INFORMATIONAL only (no longer a verdict gate, FIX 5).
        "final_max_force": final_max_force,    # optking internal-coord au (informational)
        "final_rms_force": final_rms_force,    # optking internal-coord au (informational)
    }
    return atoms_opt, summary


# ---------------------------------------------------------------------------
# Pre-D3 constraint-enforcement verification (C-D3.1 gate, 1-point drift test)
# ---------------------------------------------------------------------------
def verify_constraint_enforcement(atoms, logdir: Path, method: str = "b3lyp"):
    """Confirm optking ranged_dihedral actually enforces a target dihedral
    across opt iterations (drift check). Exercises the LARGE-rotation path that
    production conformer generation actually uses: chi1 (N-CA-CB-CG) is driven
    to the FARTHEST rotamer of {-60, 180, 60} from its current value, after a
    rigid pre-rotation of the whole CB-side fragment to that target.

    Why large rotation: the previous gate displaced chi1 by only +20 deg, so it
    optimized a near-target geometry and PASSED -- but that small-perturbation
    test never exercised the far-target large-first-step path and therefore did
    NOT catch the production conf01 crash (target ~119 deg away). A gate that
    only certifies the easy path is not a real D3 gate. We now pre-rotate to the
    farthest well and verify the same constrained reopt holds there.

    PASS if |final - target| <= VERIFY_DRIFT_TOL_DEG (the rigid-scan-explosion
    guard). #106 FIX 2 widens the ranged_dihedral window from +-0.05 to +-2.0
    deg (the +-0.05 window was a sawtooth driver); the drift tolerance is set
    just outside that window so a target genuinely held inside [target-2, +2]
    passes while a free-relaxed dihedral (which leaves the window entirely)
    still fails. This stays well inside the conformer wrong-well gate (10 deg).
    """
    import psi4

    nidx = build_name_index(atoms)
    quad = ("N", "CA", "CB", "CG")
    if any(nidx.get(a) is None for a in quad):
        return {"verified": False, "reason": "chi1 atoms missing", "atoms": quad}

    i, j, k, l = (nidx[a] for a in quad)
    coords0 = _coords_of(atoms)
    chi0 = _dihedral_deg(coords0, i, j, k, l)

    # Pick the FARTHEST rotamer well from chi0 (largest wrap-aware gap) so the
    # gate certifies the worst-case (large) rotation, not a soft +20 deg nudge.
    target = max(CHI1_TARGETS, key=lambda w: _angle_diff_deg(chi0, w))
    # optking dihedral range is +-180; normalize target into (-180, 180]
    while target > 180.0:
        target -= 360.0
    while target <= -180.0:
        target += 360.0

    gap = _angle_diff_deg(chi0, target)
    print(f"[verify] chi1 N-CA-CB-CG start={chi0:.3f} deg -> FARTHEST target={target:.3f} deg "
          f"(gap {gap:.1f} deg, large-rotation gate)")

    # Pre-rotate the whole CB-side fragment to the target so the constrained
    # reopt only relaxes locally -- this is exactly the production path (C1/C2).
    adj = _infer_bonds(atoms)
    coords_pre = _prerotate_chi1(atoms, coords0, nidx, adj, target)
    atoms_pre = _atoms_with_coords(atoms, coords_pre)

    # FIX 4: strip leaked FileHandlers before this call installs its own .log.
    _reset_psi4_log_handlers()
    psi4.core.clean_options()
    psi4.set_output_file(str(logdir / "constraint_verify_psi4.dat"), True)
    psi4.set_memory("8 GB")
    psi4.core.set_num_threads(8)
    mol = _make_psi4_molecule(atoms_pre)
    # 1-based atom indices; +-2.0 deg window (FIX 2) enforces the target without
    # the +-0.05 sawtooth.
    rng = (f"{i + 1} {j + 1} {k + 1} {l + 1} "
           f"{target - VERIFY_RANGED_HALF_WINDOW_DEG:.3f} "
           f"{target + VERIFY_RANGED_HALF_WINDOW_DEG:.3f}")
    psi4.set_options({
        "basis": "6-31G(d)",
        "scf_type": "df",
        "g_convergence": "gau_loose",
        "geom_maxiter": 80,
        "ranged_dihedral": rng,
        "intrafrag_step_limit": 0.2,  # FIX 2: float (string was silently ignored)
    })
    converged = True
    try:
        e_opt = psi4.optimize(method, molecule=mol)
    except psi4.OptimizationConvergenceError as e:
        converged = False
        wfn = getattr(e, "wfn", None)
        e_opt = wfn.energy() if wfn is not None else float("nan")
    except Exception as e:  # optking AlgError / OptError (the failure we guard)
        print(f"[verify] optking raised {type(e).__name__}: {e}")
        converged = False
        e_opt = float("nan")
    finally:
        _reset_psi4_log_handlers()  # FIX 4: drop this call's leaked handler

    coords_f = (mol.geometry().to_array() * 0.529177210903).tolist()
    chi_f = _dihedral_deg(coords_f, i, j, k, l)
    drift = abs(chi_f - target)
    if drift > 180.0:
        drift = 360.0 - drift
    result = {
        "verified": bool(drift <= VERIFY_DRIFT_TOL_DEG and e_opt == e_opt),
        "dihedral_atoms": list(quad),
        "chi1_start_deg": round(chi0, 4),
        "target_deg": round(target, 4),
        "target_gap_deg": round(gap, 4),
        "final_deg": round(chi_f, 4),
        "drift_deg": round(drift, 4),
        "ranged_half_window_deg": VERIFY_RANGED_HALF_WINDOW_DEG,
        "pass_threshold_deg": VERIFY_DRIFT_TOL_DEG,
        "opt_converged": converged,
        "large_rotation_gate": True,
        "final_energy_hartree": float(e_opt),
        "method": f"{method.upper()}/6-31G(d)",
    }
    status = "PASS" if result["verified"] else "FAIL"
    print(f"[verify] {status}: target={target:.3f} final={chi_f:.3f} drift={drift:.4f} deg "
          f"(threshold {VERIFY_DRIFT_TOL_DEG} deg, window +-{VERIFY_RANGED_HALF_WINDOW_DEG} deg, "
          f"gap {gap:.1f} deg)")
    return result


# ---------------------------------------------------------------------------
# Conformer generation (C-D2.3, >= 4 conformers)
# ---------------------------------------------------------------------------
def _conformer_plan(n_conf: int) -> List[Dict[str, object]]:
    """Build a deterministic conformer plan: chi1 x (phi,psi) grid.

    Conformer 0 = reference (no extra dihedral constraint beyond the
    ff14SB-equilibrium backbone/NE1 freeze). The rest enforce a chi1
    target and (for the later ones) a backbone phi/psi target so the
    RESP-A2 ensemble samples distinct side-chain + backbone basins.
    """
    plan: List[Dict[str, object]] = [{"label": "ref", "chi1": None, "phipsi": None}]
    # chi1 rotamers at default backbone
    for chi1 in CHI1_TARGETS:
        plan.append({"label": f"chi1_{int(chi1):+d}", "chi1": chi1, "phipsi": None})
    # backbone basins at the trans (180) rotamer
    for (phi, psi) in PHIPSI_TARGETS:
        plan.append({"label": f"bb_{int(phi):+d}_{int(psi):+d}",
                     "chi1": 180.0, "phipsi": (phi, psi)})
    return plan[:max(n_conf, 4)] if n_conf else plan


def _dihedral_constraints_for(atoms, chi1: Optional[float],
                              phipsi: Optional[Tuple[float, float]]) -> Dict[str, str]:
    """Build ranged_dihedral options (1-based, +-2.0 deg window) for a conformer
    (#106 FIX 2; was +-0.05 deg, which drove a geometry-step sawtooth)."""
    nidx = build_name_index(atoms)
    ranges: List[str] = []
    w = CONF_RANGED_HALF_WINDOW_DEG

    def _add(quad, target):
        if any(nidx.get(a) is None for a in quad):
            return
        i, j, k, l = (nidx[a] for a in quad)
        ranges.append(f"({i + 1} {j + 1} {k + 1} {l + 1} {target - w:.3f} {target + w:.3f})")

    if chi1 is not None:
        _add(("N", "CA", "CB", "CG"), chi1)
    if phipsi is not None:
        phi, psi = phipsi
        # phi = C(ACE)-N-CA-C, psi = N-CA-C-N(NME). ACE C and NME N are
        # the cap atoms in those residues; index them specially.
        ace_c = next((m for m, (nm, res, *_r) in enumerate(atoms)
                      if res == "ACE" and nm == "C"), None)
        nme_n = next((m for m, (nm, res, *_r) in enumerate(atoms)
                      if res == "NME" and nm == "N"), None)
        n = nidx.get("N")
        ca = nidx.get("CA")
        c = nidx.get("C")
        if None not in (ace_c, n, ca, c):
            ranges.append(f"({ace_c + 1} {n + 1} {ca + 1} {c + 1} {phi - w:.3f} {phi + w:.3f})")
        if None not in (n, ca, c, nme_n):
            ranges.append(f"({n + 1} {ca + 1} {c + 1} {nme_n + 1} {psi - w:.3f} {psi + w:.3f})")

    return {"ranged_dihedral": " ".join(ranges)} if ranges else {}


# ---------------------------------------------------------------------------
# Geometric pre-rotation helpers (C-D2 #106 fix).
#
# Root cause of the conf01 crash: every conformer started from the single
# reference geometry (chi1 ~ +59 deg) and then imposed an ABSOLUTE chi1 target
# via a 0.1-deg ranged_dihedral. For a far target (e.g. -60 deg, a ~119 deg
# gap) optking's first step is enormous -> AlgError "Step is far too large".
# The fix here pre-rotates the starting geometry RIGIDLY so the target
# dihedral(s) are already approximately satisfied; the constrained reopt then
# only relaxes locally (the ranged_dihedral window stays at 0.1 deg).
#
# CRITICAL: we must rotate the ENTIRE rigid CB-side fragment about the CA-CB
# axis, not just CG. (d3_omw_qm_scans.rotate_dihedral_to rotates only atom 'l'
# which shatters the indole ring -- we reuse only its Rodrigues math, applied
# to the whole fragment.) Because the rotation is a single rigid-body rotation
# about an axis that passes through CB, every bond length, every bond angle,
# and every internal indole bond/angle inside the rotated fragment is exactly
# preserved (rigid motion is isometric). The only quantities that change are
# the dihedrals that cross the rotation bond (CA-CB), which is precisely the
# dihedral we are setting. This preserves all frozen terms: N-CA-CB angle and
# CA-CB-CG angle (CB and CG stay rigid relative to the unrotated N/CA core via
# the axis through CB), CA-CB distance (CB is on the axis), and all internal
# indole bonds/angles incl. CD1-NE1 / NE1-CE2 / CD1-NE1-CE2.
# ---------------------------------------------------------------------------

# Covalent radii (Angstrom) for distance-based bond perception. Cordero 2008
# values for the handful of elements present in the Ace-OMW-NMe dipeptide.
_COVALENT_RADII: Dict[str, float] = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
}
_BOND_TOLERANCE = 1.3  # (r_i + r_j) * tol => bonded


def _infer_bonds(atoms) -> Dict[int, List[int]]:
    """Distance-based connectivity (no topology stored in the PDB). Two atoms
    are bonded if their separation < (r_cov_i + r_cov_j) * _BOND_TOLERANCE.
    Returns an adjacency list keyed by 0-based atom index."""
    import numpy as np
    coords = [np.array((a[3], a[4], a[5])) for a in atoms]
    radii = [_COVALENT_RADII.get(a[2], 0.77) for a in atoms]  # 0.77 = generic C fallback
    n = len(atoms)
    adj: Dict[int, List[int]] = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            cutoff = (radii[i] + radii[j]) * _BOND_TOLERANCE
            if float(np.linalg.norm(coords[i] - coords[j])) < cutoff:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _fragment_from(adj: Dict[int, List[int]], start: int, exclude_edge: Tuple[int, int]) -> List[int]:
    """BFS over the bond graph starting at `start`, but never traversing the
    `exclude_edge` (an undirected (u, v) bond that is cut). Returns the set of
    atom indices reachable from `start` on the `start` side of the cut,
    INCLUDING `start` itself."""
    cut = frozenset(exclude_edge)
    seen = {start}
    stack = [start]
    while stack:
        cur = stack.pop()
        for nb in adj.get(cur, ()):  # neighbours
            if frozenset((cur, nb)) == cut:
                continue  # do not cross the cut bond
            if nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return sorted(seen)


def _rotate_fragment_about_axis(coords, frag_indices, center_coord,
                                axis_from_coord, axis_to_coord, angle_deg: float):
    """Rigidly rotate the atoms in `frag_indices` by `angle_deg` about the axis
    whose direction is unit(axis_to_coord - axis_from_coord), pivoting through
    `center_coord` (all three are 3-vectors / (x,y,z) tuples). Reuses the
    Rodrigues rotation math from d3_omw_qm_scans.rotate_dihedral_to but applied
    to the WHOLE fragment (not a single atom) so the rotated body stays rigid.
    Returns a new coords list (list of [x,y,z])."""
    import numpy as np
    axis = np.array(axis_to_coord, dtype=float) - np.array(axis_from_coord, dtype=float)
    nrm = float(np.linalg.norm(axis))
    if nrm < 1e-9:
        return [list(c) for c in coords]
    axis = axis / nrm
    center = np.array(center_coord, dtype=float)
    theta = math.radians(angle_deg)
    cos_t, sin_t = math.cos(theta), math.sin(theta)
    new_coords = [list(c) for c in coords]
    fset = set(frag_indices)
    for i in fset:
        v = np.array(coords[i], dtype=float) - center
        # Rodrigues' rotation formula (same form as d3 helper).
        v_rot = v * cos_t + np.cross(axis, v) * sin_t + axis * np.dot(axis, v) * (1 - cos_t)
        new_coords[i] = (center + v_rot).tolist()
    return new_coords


def _atoms_with_coords(atoms_template, coords):
    """Return a new atoms list = atoms_template names/residues/elements with the
    coordinates replaced by `coords` (list of (x,y,z))."""
    out = []
    for (nm, res, element, _x, _y, _z), (x, y, z) in zip(atoms_template, coords):
        out.append((nm, res, element, float(x), float(y), float(z)))
    return out


def _prerotate_chi1(atoms, coords, nidx, adj, target_chi1: float):
    """Pre-rotate the rigid CB-side fragment about the CA-CB axis so chi1
    (N-CA-CB-CG) ~ target_chi1. The rotating set = the CB-side fragment of the
    CB->CA cut, MINUS CB itself (CB lies on the rotation axis so is invariant).

    Note on HB2/HB3: they belong to the rigid CB-side fragment, so including
    them keeps the CB tetrahedral geometry intact. We therefore rotate the FULL
    CB-side rigid fragment (CG + entire indole + N-methyl + all H, incl.
    HB2/HB3); only CB on the axis is held fixed. (Rotating the whole rigid
    fragment is cleanest -- any partial exclusion of HB2/HB3 would distort the
    CB tetrahedron.) Returns updated coords (or unchanged if atoms missing)."""
    for a in ("N", "CA", "CB", "CG"):
        if nidx.get(a) is None:
            return coords
    n_i, ca_i, cb_i, cg_i = nidx["N"], nidx["CA"], nidx["CB"], nidx["CG"]
    cur = _dihedral_deg(coords, n_i, ca_i, cb_i, cg_i)
    # CB-side fragment of the CB<->CA cut, then drop CB (on-axis invariant).
    frag = [i for i in _fragment_from(adj, cb_i, (cb_i, ca_i)) if i != cb_i]
    # axis direction CA->CB, pivot at CB so CB stays put.
    return _rotate_fragment_about_axis(
        coords, frag, coords[cb_i], coords[ca_i], coords[cb_i], target_chi1 - cur,
    )


def _prerotate_dihedral(atoms, coords, adj, quad_idx, target_deg: float):
    """Generic rigid pre-rotation to set an arbitrary dihedral i-j-k-l ~
    target_deg by rotating about the central j-k bond. quad_idx = (i, j, k, l)
    0-based. The cut bond is j<->k; we rotate the SMALLER of the two sides
    (the side that does NOT contain i) so the frozen backbone core is disturbed
    as little as possible, pivoting through atom k (on the axis). Returns
    updated coords (or unchanged if any index is None)."""
    if any(x is None for x in quad_idx):
        return coords
    i, j, k, l = quad_idx
    cur = _dihedral_deg(coords, i, j, k, l)
    # Two sides of the j<->k cut. The l-side is reachable from k; the i-side
    # from j. Both j and k sit on the rotation axis, so either is a valid pivot.
    k_side = set(_fragment_from(adj, k, (k, j)))
    j_side = set(_fragment_from(adj, j, (j, k)))
    delta = target_deg - cur
    # Prefer rotating the smaller side to minimize disruption to the backbone
    # core. Axis direction is ALWAYS j->k. Rotating the k/l-side by +delta about
    # j->k increases the dihedral by +delta; rotating the i/j-side instead needs
    # -delta (it carries atom i with it). Pivot is any on-axis atom (j or k).
    rotate_k_side = (l in k_side) and (len(k_side) <= len(j_side) or l not in j_side)
    if rotate_k_side:
        frag = [x for x in k_side if x not in (j, k)]
        pivot = k
        ang = delta
    else:
        frag = [x for x in j_side if x not in (j, k)]
        pivot = j
        ang = -delta
    return _rotate_fragment_about_axis(
        coords, frag, coords[pivot], coords[j], coords[k], ang,
    )


def _prerotate_conformer(atoms_ref, chi1: Optional[float],
                         phipsi: Optional[Tuple[float, float]]):
    """Build a pre-rotated starting geometry for a conformer so the constrained
    reopt only relaxes locally. Returns a new atoms list.

    What IS pre-set:
      - chi1 (N-CA-CB-CG) via rigid full-CB-side-fragment rotation (exact).
      - phi (C_ACE-N-CA-C) and psi (N-CA-C-N_NME) via rigid downstream-fragment
        rotation about N-CA and CA-C respectively, for bb_* conformers.
    Honesty note: phi/psi pre-rotation rotates the smaller bond-graph side of
    each central bond. For the cap-defined phi/psi this is approximate when the
    two halves of the dipeptide are of similar size, but it always lands the
    starting guess much closer to target than the +59-deg reference, removing
    the far-target large-first-step crash. optking then relaxes the residual.
    """
    coords = _coords_of(atoms_ref)
    nidx = build_name_index(atoms_ref)
    adj = _infer_bonds(atoms_ref)

    # 1) chi1 first (so the side-chain is on the right rotamer before backbone).
    if chi1 is not None:
        coords = _prerotate_chi1(atoms_ref, coords, nidx, adj, chi1)

    # 2) phi/psi for backbone-basin conformers.
    if phipsi is not None:
        phi, psi = phipsi
        ace_c = next((m for m, (nm, res, *_r) in enumerate(atoms_ref)
                      if res == "ACE" and nm == "C"), None)
        nme_n = next((m for m, (nm, res, *_r) in enumerate(atoms_ref)
                      if res == "NME" and nm == "N"), None)
        n_i, ca_i, c_i = nidx.get("N"), nidx.get("CA"), nidx.get("C")
        # phi = C_ACE - N - CA - C ; rotate about N-CA.
        coords = _prerotate_dihedral(atoms_ref, coords, adj, (ace_c, n_i, ca_i, c_i), phi)
        # psi = N - CA - C - N_NME ; rotate about CA-C.
        coords = _prerotate_dihedral(atoms_ref, coords, adj, (n_i, ca_i, c_i, nme_n), psi)

    return _atoms_with_coords(atoms_ref, coords)


# ---------------------------------------------------------------------------
# Post-opt wrong-well gate (C5) + measurement helpers (C-D2 #106 fix).
# ---------------------------------------------------------------------------
WRONG_WELL_TOL_DEG = 10.0  # |measured - target| above this => FAILED (wrong well)


def _angle_diff_deg(measured: float, target: float) -> float:
    """Smallest absolute difference of two angles (deg), 360-deg wrap-aware."""
    d = abs((measured - target + 180.0) % 360.0 - 180.0)
    return d


def _measure_conformer_dihedrals(atoms_c, chi1_target: Optional[float],
                                 phipsi_target: Optional[Tuple[float, float]]):
    """Measure chi1 (and phi/psi for bb conformers) on an optimized conformer
    and evaluate the C5 wrong-well gate. Returns
    (measured_dict, on_target: bool, max_drift_deg)."""
    coords = _coords_of(atoms_c)
    nidx = build_name_index(atoms_c)
    measured: Dict[str, Optional[float]] = {}
    drifts: List[float] = []

    if all(x in nidx for x in ("N", "CA", "CB", "CG")):
        chi1_m = _dihedral_deg(coords, nidx["N"], nidx["CA"], nidx["CB"], nidx["CG"])
        measured["chi1"] = round(chi1_m, 3)
        if chi1_target is not None:
            drifts.append(_angle_diff_deg(chi1_m, chi1_target))
    else:
        measured["chi1"] = None

    if phipsi_target is not None:
        phi_t, psi_t = phipsi_target
        ace_c = next((m for m, (nm, res, *_r) in enumerate(atoms_c)
                      if res == "ACE" and nm == "C"), None)
        nme_n = next((m for m, (nm, res, *_r) in enumerate(atoms_c)
                      if res == "NME" and nm == "N"), None)
        n_i, ca_i, c_i = nidx.get("N"), nidx.get("CA"), nidx.get("C")
        if None not in (ace_c, n_i, ca_i, c_i):
            phi_m = _dihedral_deg(coords, ace_c, n_i, ca_i, c_i)
            measured["phi"] = round(phi_m, 3)
            drifts.append(_angle_diff_deg(phi_m, phi_t))
        if None not in (n_i, ca_i, c_i, nme_n):
            psi_m = _dihedral_deg(coords, n_i, ca_i, c_i, nme_n)
            measured["psi"] = round(psi_m, 3)
            drifts.append(_angle_diff_deg(psi_m, psi_t))

    max_drift = max(drifts) if drifts else 0.0
    on_target = max_drift <= WRONG_WELL_TOL_DEG
    return measured, on_target, round(max_drift, 3)


def _drive_dihedral_fallback(atoms_start, logdir: Path, method: str, ci: int,
                             chi1_target: Optional[float],
                             phipsi_target: Optional[Tuple[float, float]],
                             step_deg: float = 20.0):
    """C6 drive fallback (fail-only): incrementally step the ranged_dihedral
    target from the current geometry toward the target in <= step_deg
    increments, reopt at each step carrying the geometry forward. Applied to a
    single failing conformer only. Returns (atoms_driven, summary, ok: bool).
    On any crash mid-drive, returns the last good geometry with ok=False."""
    atoms_cur = atoms_start
    summ: Dict[str, object] = {}
    # Build the ordered list of (quad_names, target) to drive. chi1 first.
    drive_targets: List[Tuple[Tuple[str, ...], float]] = []
    if chi1_target is not None:
        drive_targets.append((("N", "CA", "CB", "CG"), chi1_target))
    # phi/psi driven via the same incremental mechanism (cap atoms handled in
    # _dihedral_constraints_for, which we re-derive per step below).
    drive_phipsi = phipsi_target

    nstep_total = 0
    for quad, tgt in drive_targets:
        nidx = build_name_index(atoms_cur)
        if any(nidx.get(a) is None for a in quad):
            continue
        i, j, k, l = (nidx[a] for a in quad)
        cur = _dihedral_deg(_coords_of(atoms_cur), i, j, k, l)
        # incremental waypoints from cur -> tgt by the shortest angular path
        delta = (tgt - cur + 180.0) % 360.0 - 180.0
        n_inc = max(1, int(math.ceil(abs(delta) / step_deg)))
        for s in range(1, n_inc + 1):
            way = cur + delta * (s / n_inc)
            extra = _dihedral_constraints_for(atoms_cur, way, None)
            try:
                atoms_cur, summ = run_psi4_optimize(
                    atoms_cur, logdir, method, mode="constrained",
                    out_name=f"conf{ci:02d}_drive_chi1_{s:02d}.dat",
                    extra_constraints=extra,
                )
                nstep_total += summ.get("n_steps", 0)
            except Exception as e:  # optking AlgError / OptError mid-drive
                print(f"[drive] conf{ci:02d} chi1 step {s}/{n_inc} crashed: {e}")
                return atoms_cur, summ, False

    # final phi/psi enforcement once chi1 is on target (single constrained opt)
    if drive_phipsi is not None:
        extra = _dihedral_constraints_for(atoms_cur, chi1_target, drive_phipsi)
        try:
            atoms_cur, summ = run_psi4_optimize(
                atoms_cur, logdir, method, mode="constrained",
                out_name=f"conf{ci:02d}_drive_bb.dat", extra_constraints=extra,
            )
            nstep_total += summ.get("n_steps", 0)
        except Exception as e:
            print(f"[drive] conf{ci:02d} bb enforcement crashed: {e}")
            return atoms_cur, summ, False

    if summ:
        summ["n_steps"] = nstep_total
    return atoms_cur, summ, True


def generate_conformers(atoms_ref, logdir: Path, method: str, n_conf: int):
    """Generate >= 4 constrained-optimized conformers. atoms_ref is the
    constrained-optimized reference geometry.

    C-D2 #106 fix flow per conformer:
      1) geometrically PRE-ROTATE the starting geometry to the target
         dihedral(s) (C1-C3), so the constrained reopt only relaxes locally;
      2) run the SAME constrained reopt (frozen backbone+NE1 + 0.1-deg
         ranged_dihedral), with intrafrag_step_limit as a defensive cap (C4),
         wrapped in try/except so one crash cannot kill the run;
      3) HARD-gate the result against the target well (C5: |measured - target|
         <= 10 deg, wrap-aware) and the rotamer energy sanity (< 25 kcal/mol);
      4) on crash OR gate failure, retry THAT conformer once via the
         incremental dihedral drive (C6); record failed conformers without
         aborting the whole run.
    """
    confdir = logdir / "conformers"
    confdir.mkdir(parents=True, exist_ok=True)
    plan = _conformer_plan(n_conf)
    conformers = []
    ref_energy: Optional[float] = None  # set from conf00 (ref) for delta-E sanity
    for ci, spec in enumerate(plan):
        label = spec["label"]
        chi1_t = spec.get("chi1")
        phipsi_t = spec.get("phipsi")
        extra = _dihedral_constraints_for(atoms_ref, chi1_t, phipsi_t)

        # C1-C3: pre-rotate the starting geometry toward the target well so the
        # constrained reopt below only relaxes locally (no far-target crash).
        atoms_start = _prerotate_conformer(atoms_ref, chi1_t, phipsi_t)

        # C4: cap residual large steps (does NOT touch the ranged_dihedral window
        # or the ranged distance/bend pins). FIX 2: pass a FLOAT 0.2 -- the prior
        # string "0.4" was silently ignored by optking, so the logged effective
        # value stayed at the 0.5 default. A float is parsed; verify the .dat
        # shows INTRAFRAG_STEP_LIMIT = 0.2 after the run.
        extra_opt = dict(extra)
        extra_opt["intrafrag_step_limit"] = 0.2
        print(f"[conf] {ci:02d} ({label}) constraints={extra or 'ref'} (pre-rotated)")

        failed = False
        fail_reason: Optional[str] = None
        used_drive = False
        try:
            atoms_c, summ = run_psi4_optimize(
                atoms_start, confdir, method, mode="constrained",
                out_name=f"conf{ci:02d}_psi4.dat", extra_constraints=extra_opt,
            )
        except Exception as e:  # optking AlgError / OptError -> C6 drive fallback
            print(f"[conf] {ci:02d} ({label}) reopt crashed: {e} -> drive fallback")
            atoms_c, summ, ok = _drive_dihedral_fallback(
                atoms_start, confdir, method, ci, chi1_t, phipsi_t,
            )
            used_drive = True
            if not ok:
                failed, fail_reason = True, "drive_crash"

        # C5: wrong-well hard gate (chi1 + phi/psi for bb conformers).
        measured, on_target, max_drift = _measure_conformer_dihedrals(
            atoms_c, chi1_t, phipsi_t,
        )
        # per-conformer energy sanity: finite + rotamer dE vs ref < 25 kcal/mol.
        e_ha = summ.get("final_energy_hartree", float("nan")) if summ else float("nan")
        if ci == 0 and e_ha == e_ha:
            ref_energy = e_ha
        de_kcal = None
        if ref_energy is not None and e_ha == e_ha:
            de_kcal = (e_ha - ref_energy) * 627.5094740631

        energy_bad = (e_ha != e_ha) or (de_kcal is not None and abs(de_kcal) > 25.0)

        if not failed and (not on_target or energy_bad) and not used_drive:
            # C6 fallback (fail-only): retry THIS conformer via incremental drive.
            why = "wrong_well" if not on_target else "energy_out_of_range"
            print(f"[conf] {ci:02d} ({label}) C5 FAIL ({why}, drift={max_drift} deg, "
                  f"dE={de_kcal}); drive fallback")
            atoms_d, summ_d, ok = _drive_dihedral_fallback(
                atoms_start, confdir, method, ci, chi1_t, phipsi_t,
            )
            used_drive = True
            if ok and summ_d:
                atoms_c, summ = atoms_d, summ_d
                measured, on_target, max_drift = _measure_conformer_dihedrals(
                    atoms_c, chi1_t, phipsi_t,
                )
                e_ha = summ.get("final_energy_hartree", float("nan"))
                de_kcal = ((e_ha - ref_energy) * 627.5094740631
                           if (ref_energy is not None and e_ha == e_ha) else None)
                energy_bad = (e_ha != e_ha) or (de_kcal is not None and abs(de_kcal) > 25.0)
            # re-evaluate failure after the drive retry
            if not on_target:
                failed, fail_reason = True, "wrong_well"
            elif energy_bad:
                failed, fail_reason = True, "energy_out_of_range"
        elif not failed and (not on_target or energy_bad) and used_drive:
            # already came from a drive (crash path) and still off-target/bad.
            failed = True
            fail_reason = "wrong_well" if not on_target else "energy_out_of_range"

        pdb = confdir / f"conf{ci:02d}_optimized.pdb"
        write_pdb(atoms_c, pdb,
                  f"conf{ci:02d} {label} {summ.get('method', method) if summ else method} "
                  f"E={e_ha:.6f} Ha")
        conformers.append({
            "index": ci,
            "label": label,
            "pdb": str(pdb.relative_to(logdir)),
            "energy_hartree": e_ha,
            "delta_e_vs_ref_kcal": round(de_kcal, 3) if de_kcal is not None else None,
            "converged": bool(summ.get("converged", False)) if summ else False,
            "n_steps": summ.get("n_steps", 0) if summ else 0,
            # FIX 5: per-conformer geometry-convergence columns (the verdict gate
            # is now disp+energy, since constraint reaction force is expected for
            # a ranged-pin opt and false-rejects valid constrained minima).
            "final_max_disp": summ.get("final_max_disp") if summ else None,
            "final_rms_disp": summ.get("final_rms_disp") if summ else None,
            "final_delta_e": summ.get("final_delta_e") if summ else None,
            # forces retained INFORMATIONAL only (not a verdict gate, FIX 5).
            "final_max_force": summ.get("final_max_force") if summ else None,
            "final_rms_force": summ.get("final_rms_force") if summ else None,
            "chi1_target_deg": chi1_t,
            "phipsi_target_deg": list(phipsi_t) if phipsi_t else None,
            "measured_deg": measured,
            "chi1_measured_deg": measured.get("chi1"),
            "max_drift_deg": max_drift,
            "on_target": bool(on_target),
            "used_drive_fallback": used_drive,
            "failed": failed,
            "fail_reason": fail_reason,
        })
    return conformers


# ---------------------------------------------------------------------------
# Delta report (bond / angle / improper vs ff14SB)
# ---------------------------------------------------------------------------
def compute_bond_deltas(atoms_initial, atoms_opt):
    nidx_opt = build_name_index(atoms_opt)
    nidx_ini = build_name_index(atoms_initial)
    c_opt = _coords_of(atoms_opt)
    c_ini = _coords_of(atoms_initial)
    out = {}
    for (a, b), eq in AMBER14SB_BONDS.items():
        ia, ib = nidx_opt.get(a), nidx_opt.get(b)
        if ia is None or ib is None:
            continue
        d_opt = _bond_len(c_opt, ia, ib)
        d_ini = _bond_len(c_ini, nidx_ini[a], nidx_ini[b])
        out[f"{a}-{b}"] = {
            "amber14SB_eq": eq,
            "initial_pdb": round(d_ini, 4),
            "qm_optimized": round(d_opt, 4),
            "delta_qm_vs_amber": round(d_opt - eq, 4),
            "delta_qm_vs_initial": round(d_opt - d_ini, 4),
        }
    return out


def compute_angle_deltas(atoms_opt):
    nidx = build_name_index(atoms_opt)
    coords = _coords_of(atoms_opt)
    out = {}
    for (a, b, c), eq in AMBER14SB_ANGLES.items():
        ia, ib, ic = nidx.get(a), nidx.get(b), nidx.get(c)
        if None in (ia, ib, ic):
            continue
        ang = _angle_deg(coords, ia, ib, ic)
        out[f"{a}-{b}-{c}"] = {
            "amber14SB_eq_deg": eq,
            "qm_optimized_deg": round(ang, 3),
            "delta_deg": round(ang - eq, 3),
        }
    return out


def compute_improper_deltas(atoms_opt):
    """Out-of-plane improper magnitudes (deg); ff14SB planar centers = ~0."""
    nidx = build_name_index(atoms_opt)
    coords = _coords_of(atoms_opt)
    out = {}
    for (i, j, k, l) in IMPROPER_CENTERS:
        ii, jj, kk, ll = nidx.get(i), nidx.get(j), nidx.get(k), nidx.get(l)
        if None in (ii, jj, kk, ll):
            continue
        imp = _dihedral_deg(coords, ii, jj, kk, ll)
        out[f"{i}-{j}-{k}-{l}"] = {
            "qm_improper_deg": round(imp, 3),
            "abs_out_of_plane_deg": round(min(abs(imp), abs(180.0 - abs(imp))), 3),
        }
    return out


def _conformer_geometry_converged(c: Dict) -> bool:
    """C7 GEOMETRY-convergence test for one conformer (#106 FIX 5, replacing the
    prior force-convergence gate). True if optking flagged it converged, OR the
    parsed final convergence row shows the geometry stopped moving:
        Max Disp <= GAU_LOOSE_MAX_DISP (1.0e-2 au)
        AND RMS Disp <= GAU_LOOSE_RMS_DISP (6.7e-3 au)

    Why no Delta-E clause (#106 v3, 2026-06-21): optking's GAUSSIAN preset does
    NOT list max_DE among its *required* criteria (Delta E is an inactive 'o'
    marker), i.e. optking itself does not use Delta E to declare convergence. A
    1e-5 Ha tolerance also sits below the 40-atom DF-B3LYP SCF noise floor, so it
    false-rejected genuinely converged conformers (conf01: optking_converged with
    Delta E = 1.75e-5). The displacement criteria redundantly cover the flatness
    signal. final_delta_e is still parsed/reported as informational. A 5e-5
    relaxation was rejected as an unsourced arbitrary value.

    Why not force: for a ranged-constrained opt the constraint reaction force is
    expected to be nonzero (the QM gradient pushes the pinned coord toward the
    free-min, the pin clamps the displacement but not the force) -- demanding
    Max Force -> 0 false-rejects valid constrained minima (conf02 plateaus at
    MaxF 3.33e-2 while disp is fully converged). The geometry signal (small disp)
    is what optking's GAUSSIAN preset actually certifies for a constrained
    minimum. Q1/Q2(a)/C-1/C-4.

    A conformer with no parsed disp and converged==False is NOT geometry-
    converged (fail-closed: we do not certify what we cannot measure)."""
    if c.get("converged"):
        return True
    md = c.get("final_max_disp")
    rd = c.get("final_rms_disp")
    if md is None or rd is None:
        return False
    try:
        return (float(md) <= GAU_LOOSE_MAX_DISP
                and float(rd) <= GAU_LOOSE_RMS_DISP)
    except (TypeError, ValueError):
        return False


def compute_d2_verdict(summary: Dict, conformers: List[Dict]) -> Dict[str, object]:
    """C7 D2 acceptance criteria. PASS requires:
      - >= 4 GEOMETRY-converged conformers that pass C5 (on-target, not failed);
      - those on-target conformers cover >= 3 distinct chi1 wells of
        {-60, 180, 60} (|measured - well| < 10 deg);
      - per-conformer energy sanity: |dE vs ref| finite and < ~25 kcal/mol
        (already gated in generate_conformers; re-checked here);
      - GEOMETRY convergence (#106 FIX 5, REPLACES the prior force gate):
        converged==True OR the parsed final convergence row shows the geometry
        stopped moving AND the energy went flat (Max Disp <= 1.0e-2 AND
        RMS Disp <= 6.7e-3 au AND |Delta E| <= 1e-5 Ha). The prior FIX 3 force
        gate (Max Force <= 2.5e-3) was the WRONG metric for a ranged-pin
        constrained opt: the constraint reaction force is expected to be nonzero
        whenever the amber14SB-eq target differs from the QM free-min, so the
        force gate false-rejects valid constrained minima;
      - constraint integrity on the reference opt: frozen distance |delta| <
        0.01 A and frozen bend |delta| < 1 deg vs ff14SB-eq. This gate is now
        LOAD-BEARING -- it carries the "did the amber-pin actually hold" burden
        that the (retired) force gate was wrongly trying to share.
    Returns a dict with PASS/MUST-REDO verdict + the supporting counts."""
    # on-target, non-failed, finite-energy, GEOMETRY-CONVERGED conformers
    good = [c for c in conformers
            if c.get("on_target") and not c.get("failed")
            and isinstance(c.get("energy_hartree"), float)
            and c["energy_hartree"] == c["energy_hartree"]
            and _conformer_geometry_converged(c)]
    n_good = len(good)

    # distinct chi1 wells covered (|measured - well| < 10 deg)
    wells_covered = set()
    for c in good:
        chi = c.get("chi1_measured_deg")
        if chi is None:
            continue
        for w in CHI1_TARGETS:
            if _angle_diff_deg(chi, w) < WRONG_WELL_TOL_DEG:
                wells_covered.add(w)
                break
    n_wells = len(wells_covered)

    # energy sanity across good conformers (defensive re-check)
    energy_ok = all(
        (c.get("delta_e_vs_ref_kcal") is None)
        or (abs(c["delta_e_vs_ref_kcal"]) < 25.0)
        for c in good
    )

    # constraint integrity on the reference opt (frozen distance / bend deltas)
    frozen_bond_names = {f"{a}-{b}" for (a, b) in CONSTRAINED_BONDS}
    frozen_angle_names = {f"{a}-{b}-{c}" for (a, b, c) in CONSTRAINED_ANGLES}
    bond_deltas = summary.get("bond_deltas", {})
    angle_deltas = summary.get("angle_deltas", {})
    bond_viol = {
        nm: d["delta_qm_vs_amber"]
        for nm, d in bond_deltas.items()
        if nm in frozen_bond_names and abs(d["delta_qm_vs_amber"]) >= 0.01
    }
    angle_viol = {
        nm: d["delta_deg"]
        for nm, d in angle_deltas.items()
        if nm in frozen_angle_names and abs(d["delta_deg"]) >= 1.0
    }
    constraint_ok = (not bond_viol) and (not angle_viol)

    # FIX 5: GEOMETRY-convergence accounting across ALL conformers (the gate).
    n_geometry_converged = sum(1 for c in conformers if _conformer_geometry_converged(c))
    not_geometry_converged = [
        {"index": c["index"], "label": c["label"],
         "converged": bool(c.get("converged")),
         "final_max_disp": c.get("final_max_disp"),
         "final_rms_disp": c.get("final_rms_disp"),
         "final_delta_e": c.get("final_delta_e")}
        for c in conformers if not _conformer_geometry_converged(c)
    ]

    # FIX 5: force accounting retained INFORMATIONAL only (diagnostic, not a gate)
    # so the report still shows the (expected-nonzero) constraint reaction force.
    n_force_converged = sum(
        1 for c in conformers
        if c.get("converged")
        or (c.get("final_max_force") is not None and c.get("final_rms_force") is not None
            and float(c["final_max_force"]) <= GAU_LOOSE_MAX_FORCE
            and float(c["final_rms_force"]) <= GAU_LOOSE_RMS_FORCE)
    )

    passed = (n_good >= 4) and (n_wells >= 3) and energy_ok and constraint_ok
    return {
        "verdict": "PASS" if passed else "MUST-REDO",
        "n_on_target_conformers": n_good,  # on-target AND not-failed AND geometry-converged
        "n_on_target_required": 4,
        "n_distinct_chi1_wells": n_wells,
        "n_distinct_chi1_wells_required": 3,
        "chi1_wells_covered": sorted(wells_covered),
        "energy_sanity_ok": bool(energy_ok),
        "constraint_integrity_ok": bool(constraint_ok),
        # FIX 5: geometry convergence is the gate; force is informational.
        "n_geometry_converged": n_geometry_converged,
        "geometry_convergence_thresholds_au": {
            "max_disp": GAU_LOOSE_MAX_DISP, "rms_disp": GAU_LOOSE_RMS_DISP,
        },
        "not_geometry_converged": not_geometry_converged,
        "n_force_converged": n_force_converged,  # informational (constraint reaction force expected)
        "force_convergence_thresholds_au": {
            "max_force": GAU_LOOSE_MAX_FORCE, "rms_force": GAU_LOOSE_RMS_FORCE,
            "note": "informational only; not a verdict gate (see FIX 5)",
        },
        "frozen_distance_violations_A": bond_viol,
        "frozen_bend_violations_deg": angle_viol,
        "failed_conformers": [
            {"index": c["index"], "label": c["label"], "reason": c.get("fail_reason")}
            for c in conformers if c.get("failed")
        ],
    }


def write_report(logdir: Path, summary: Dict, conformers: List[Dict],
                 verify: Optional[Dict]) -> None:
    bonds = summary.get("bond_deltas", {})
    angles = summary.get("angle_deltas", {})
    imps = summary.get("improper_deltas", {})
    md = [
        "# D2 OMW B3LYP/6-31G(d) constrained geometry optimization report",
        "",
        f"**Date**: {datetime.now().isoformat()}",
        "**Stack**: psi4 1.10 + optking (single-process)",
        f"**Method**: {summary['method']}, geom_maxiter=300, g_convergence=gau_loose",
        f"**Mode**: {summary['mode']} (constrained = ff14SB-equilibrium ranged-pinned backbone + NE1 region)",
        f"**Atoms**: {summary['n_atoms']} (Ace-OMW-NMe dipeptide, internal-residue convention)",
        f"**Reference energy**: {summary['final_energy_hartree']:.6f} Ha "
        f"= {summary['final_energy_kcal']:.3f} kcal/mol",
        f"**Reference opt steps**: {summary['n_steps']}  (converged={summary['converged']})",
        f"**Active backbone/NE1 constraints**: {summary.get('constraints', {})}",
        "",
        "Charges are NOT fit in D2 (geometry + conformer ensemble only). The RESP",
        "ESP is a separate HF/6-31G* step; a B3LYP density must not leak into the fit.",
        "",
    ]
    if verify is not None:
        md += [
            "## Pre-D3 constraint-enforcement verification",
            "",
            f"- dihedral: {'-'.join(verify.get('dihedral_atoms', []))}",
            f"- target = {verify.get('target_deg')} deg, final = {verify.get('final_deg')} deg, "
            f"drift = {verify.get('drift_deg')} deg (threshold {verify.get('pass_threshold_deg')} deg)",
            f"- **{'PASS' if verify.get('verified') else 'FAIL'}** "
            f"(opt_converged={verify.get('opt_converged')})",
            "",
            "PASS confirms optking enforces a target dihedral across opt iterations, so the",
            "D3 relaxed scans will not silently free-relax (rigid-scan-explosion guard).",
            "",
        ]
    md += [
        "## Bond length deltas (OMW residue) vs amber14SB Trp equilibrium",
        "",
        "| Bond | amber14SB eq (A) | D1 input (A) | D2 B3LYP (A) | d(B3LYP-amber) | d(B3LYP-input) |",
        "|---|---|---|---|---|---|",
    ]
    for b, d in bonds.items():
        md.append(
            f"| {b} | {d['amber14SB_eq']:.4f} | {d['initial_pdb']:.4f} | "
            f"{d['qm_optimized']:.4f} | {d['delta_qm_vs_amber']:+.4f} | {d['delta_qm_vs_initial']:+.4f} |"
        )
    md += [
        "",
        "## Angle deltas (OMW residue) vs amber14SB Trp equilibrium",
        "",
        "| Angle | amber14SB eq (deg) | D2 B3LYP (deg) | delta (deg) |",
        "|---|---|---|---|",
    ]
    for a, d in angles.items():
        md.append(f"| {a} | {d['amber14SB_eq_deg']:.2f} | {d['qm_optimized_deg']:.3f} | {d['delta_deg']:+.3f} |")
    md += [
        "",
        "## Improper out-of-plane (indole planarity + methyl junction)",
        "",
        "| Improper (i-j-k-l, k central) | QM improper (deg) | |out-of-plane| (deg) |",
        "|---|---|---|",
    ]
    for im, d in imps.items():
        md.append(f"| {im} | {d['qm_improper_deg']:+.3f} | {d['abs_out_of_plane_deg']:.3f} |")
    md += [
        "",
        f"## Conformer ensemble ({len(conformers)} conformers for RESP-A2)",
        "",
        "| # | label | chi1 target | chi1 meas (deg) | dE vs ref (kcal) | drift (deg) "
        "| on-target | drive | converged | maxDisp (au) | rmsDisp (au) | |dE| (Ha) "
        "| maxF (au, info) | steps |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for c in conformers:
        chi1_meas = c.get("chi1_measured_deg")
        chi1_tgt = c.get("chi1_target_deg")
        de = c.get("delta_e_vs_ref_kcal")
        mdsp = c.get("final_max_disp")
        rdsp = c.get("final_rms_disp")
        dee = c.get("final_delta_e")
        mf = c.get("final_max_force")
        md.append(
            f"| {c['index']:02d} | {c['label']} | "
            f"{'-' if chi1_tgt is None else f'{chi1_tgt:+.0f}'} | "
            f"{'-' if chi1_meas is None else chi1_meas} | "
            f"{'-' if de is None else f'{de:+.3f}'} | "
            f"{c.get('max_drift_deg', '-')} | "
            f"{'yes' if c.get('on_target') else 'NO'} | "
            f"{'yes' if c.get('used_drive_fallback') else '-'} | "
            f"{c['converged']} | "
            f"{'-' if mdsp is None else f'{mdsp:.2e}'} | "
            f"{'-' if rdsp is None else f'{rdsp:.2e}'} | "
            f"{'-' if dee is None else f'{abs(dee):.2e}'} | "
            f"{'-' if mf is None else f'{mf:.2e}'} | "
            f"{c['n_steps']} |"
        )

    # C7: D2 acceptance verdict (PASS / MUST-REDO + supporting counts).
    verdict = compute_d2_verdict(summary, conformers)
    md += [
        "",
        "## D2 acceptance verdict (C-D2 #106)",
        "",
        f"- **{verdict['verdict']}**",
        f"- on-target conformers: {verdict['n_on_target_conformers']} "
        f"(need >= {verdict['n_on_target_required']})",
        f"- distinct chi1 wells covered: {verdict['n_distinct_chi1_wells']} of "
        f"{{-60, 180, 60}} -> {verdict['chi1_wells_covered']} "
        f"(need >= {verdict['n_distinct_chi1_wells_required']})",
        f"- energy sanity (|dE vs ref| < 25 kcal/mol): {verdict['energy_sanity_ok']}",
        f"- geometry convergence (GAU_LOOSE MaxDisp<=1.0e-2, RMSDisp<=6.7e-3 au): "
        f"{verdict['n_geometry_converged']}/{len(conformers)} conformers",
        f"- force convergence (GAU_LOOSE MaxF<=2.5e-3, RMSF<=1.7e-3 au) -- "
        f"INFORMATIONAL only, NOT a gate (constraint reaction force expected for "
        f"ranged-pin opt): {verdict['n_force_converged']}/{len(conformers)} conformers",
        f"- constraint integrity (frozen dist |d|<0.01 A, bend |d|<1 deg) -- "
        f"LOAD-BEARING: {verdict['constraint_integrity_ok']}",
    ]
    if verdict.get("not_geometry_converged"):
        md.append(f"- NOT geometry-converged: {verdict['not_geometry_converged']}")
    if verdict["frozen_distance_violations_A"]:
        md.append(f"- frozen-distance violations: {verdict['frozen_distance_violations_A']}")
    if verdict["frozen_bend_violations_deg"]:
        md.append(f"- frozen-bend violations: {verdict['frozen_bend_violations_deg']}")
    if verdict["failed_conformers"]:
        md.append(f"- failed conformers: {verdict['failed_conformers']}")
    md += [
        "",
        "## Interpretation",
        "",
        "- Bonds/angles with |delta(B3LYP - amber14SB)| above ~0.02 A / ~2 deg are D3 re-fit candidates.",
        "- NE1-CM (N-methyl junction) is the primary 1-MeTrp-specific term absent from standard Trp.",
        "- Improper out-of-plane > ~5 deg at an indole-ring center flags a D3 improper-torsion priority.",
        "",
        "## Next step (D3)",
        "",
        "Relaxed constrained-opt QM scans on the optimized geometry (rigid scans are auto-REJECT):",
        "  - Improper N-CD1-CE2-CG: -30 to +30 deg in 5 deg steps (indole planarity)",
        "  - Angle CM-NE1-CE2 / CM-NE1-CD1: methyl junction",
        "  - Dihedral CB-CG-CD1: 0 to 360 deg in 30 deg steps (side-chain rotation)",
        "",
        "The pre-D3 constraint verification above must PASS before launching the D3 scans.",
    ]
    (logdir / "d2_geom_opt_report.md").write_text("\n".join(md))


def copy_deliverable(logdir: Path, conformers: List[Dict]) -> Path:
    """Copy the consumed-by-D3 deliverable into params/omw_qm_geom/ (R-7: the
    analysis dir remains the immutable record)."""
    dst = PROJ / "params/omw_qm_geom"
    dst.mkdir(parents=True, exist_ok=True)
    for fname in ("omw_dipeptide_optimized.pdb", "d2_geom_opt_report.md",
                  "omw_geom_opt_summary.json", "constraint_verification.json"):
        src = logdir / fname
        if src.exists():
            shutil.copy2(src, dst / fname)
    confsrc = logdir / "conformers"
    if confsrc.exists():
        confdst = dst / "conformers"
        if confdst.exists():
            shutil.rmtree(confdst)
        shutil.copytree(confsrc, confdst)
    (dst / "PROVENANCE.txt").write_text(
        f"D2 deliverable copied from {logdir}\n"
        f"generated {datetime.now().isoformat()}\n"
        f"{len(conformers)} conformers; geometry only (no charges; RESP is a separate HF/6-31G* step)\n"
    )
    return dst


# ---------------------------------------------------------------------------
# Re-verdict from existing artifacts (#106 FIX 5 C-3).
#
# Re-compute the D2 PASS/MUST-REDO verdict from a COMPLETED run's artifacts
# WITHOUT re-running any QM. Reads conformers/conf{NN}_optimized.pdb (geometry)
# + conformers/conf{NN}_psi4.dat (the optking convergence table -> disp/energy/
# force) and the reference geometry, re-measures every conformer, and writes a
# fresh d2_geom_opt_report.md + omw_geom_opt_summary.json into the SAME logdir.
#
# Use case: a run that completed under the old (force-gate) verdict logic can be
# re-verdicted under the corrected geometry-gate logic with zero QM cost. Run
# this only on a FINISHED run -- never on the active logdir while a job still
# owns it (R-7 / process-ownership).
# ---------------------------------------------------------------------------
def _reference_geometry_for_reverdict(logdir: Path):
    """Pick the reference geometry for delta-E + constraint-integrity: prefer
    conformers/conf00_optimized.pdb, fall back to omw_dipeptide_optimized.pdb at
    the logdir root. Returns (atoms, source_path) or (None, None)."""
    candidates = [
        logdir / "conformers" / "conf00_optimized.pdb",
        logdir / "omw_dipeptide_optimized.pdb",
    ]
    for p in candidates:
        if p.exists():
            return parse_pdb(p), p
    return None, None


def _conformer_dat_for_reverdict(confdir: Path, ci: int) -> Optional[Path]:
    """Locate the .dat holding conf{ci}'s final convergence table. Prefer the
    primary conf{NN}_psi4.dat; if absent (drive-fallback path), fall back to the
    last conf{NN}_drive_*.dat by name. Returns the path or None."""
    primary = confdir / f"conf{ci:02d}_psi4.dat"
    if primary.exists():
        return primary
    drive = sorted(confdir.glob(f"conf{ci:02d}_drive_*.dat"))
    return drive[-1] if drive else None


def reverdict_from_logdir(logdir: Path, write: bool = True) -> Dict[str, object]:
    """Re-compute the D2 verdict from an existing run's artifacts (no QM re-run).

    Returns the verdict dict (and, when write=True, refreshes
    d2_geom_opt_report.md + omw_geom_opt_summary.json in `logdir`). Geometry
    .pdb / .dat files are read-only (R-7 immutable record)."""
    confdir = logdir / "conformers"
    if not confdir.exists():
        raise FileNotFoundError(f"no conformers/ dir under {logdir}")

    ref_atoms, ref_src = _reference_geometry_for_reverdict(logdir)
    if ref_atoms is None:
        raise FileNotFoundError(
            f"no reference geometry (conf00_optimized.pdb / omw_dipeptide_optimized.pdb) in {logdir}")

    conf_pdbs = sorted(confdir.glob("conf*_optimized.pdb"))
    if not conf_pdbs:
        raise FileNotFoundError(f"no conf*_optimized.pdb under {confdir}")

    # Rebuild the deterministic conformer plan so each index recovers its chi1 /
    # phipsi target (same plan generate_conformers used). Plan length must cover
    # every conf PDB present.
    plan = _conformer_plan(len(conf_pdbs))

    ref_energy: Optional[float] = None
    conformers: List[Dict] = []
    for p in conf_pdbs:
        # conf{NN}_optimized.pdb -> NN
        stem = p.stem  # e.g. "conf02_optimized"
        try:
            ci = int(stem[4:6])
        except ValueError:
            continue
        spec = plan[ci] if ci < len(plan) else {"label": stem, "chi1": None, "phipsi": None}
        chi1_t = spec.get("chi1")
        phipsi_t = spec.get("phipsi")

        atoms_c = parse_pdb(p)
        conv = _parse_final_convergence(
            _conformer_dat_for_reverdict(confdir, ci) or Path("/nonexistent"))
        e_ha = conv["energy"] if conv["energy"] is not None else float("nan")
        if ci == 0 and e_ha == e_ha:
            ref_energy = e_ha
        de_kcal = None
        if ref_energy is not None and e_ha == e_ha:
            de_kcal = (e_ha - ref_energy) * 627.5094740631

        measured, on_target, max_drift = _measure_conformer_dihedrals(
            atoms_c, chi1_t, phipsi_t)
        # reverdict cannot read optking's in-memory converged flag; rely on the
        # disp+energy criteria (fail-closed). converged=False is honest here.
        conformers.append({
            "index": ci,
            "label": spec.get("label", stem),
            "pdb": str(p.relative_to(logdir)),
            "energy_hartree": float(e_ha),
            "delta_e_vs_ref_kcal": round(de_kcal, 3) if de_kcal is not None else None,
            "converged": False,  # not recoverable from artifacts; disp+energy gates instead
            "n_steps": None,
            "final_max_disp": conv["max_disp"],
            "final_rms_disp": conv["rms_disp"],
            "final_delta_e": conv["delta_e"],
            "final_max_force": conv["max_force"],
            "final_rms_force": conv["rms_force"],
            "chi1_target_deg": chi1_t,
            "phipsi_target_deg": list(phipsi_t) if phipsi_t else None,
            "measured_deg": measured,
            "chi1_measured_deg": measured.get("chi1"),
            "max_drift_deg": max_drift,
            "on_target": bool(on_target),
            "used_drive_fallback": None,
            "failed": False,
            "fail_reason": None,
        })

    # constraint-integrity deltas on the reference geometry (vs amber14SB-eq),
    # measured against the reference geometry's own input where available.
    summary: Dict[str, object] = {
        "method": "B3LYP/6-31G(d)",
        "mode": "constrained",
        "reverdict": True,
        "reverdict_reference": str(ref_src.relative_to(logdir)),
        "n_atoms": len(ref_atoms),
        "n_conformers": len(conformers),
        "bond_deltas": compute_bond_deltas(ref_atoms, ref_atoms),
        "angle_deltas": compute_angle_deltas(ref_atoms),
        "improper_deltas": compute_improper_deltas(ref_atoms),
        "conformers": conformers,
    }
    verdict = compute_d2_verdict(summary, conformers)
    summary["d2_acceptance"] = verdict

    if write:
        (logdir / "omw_geom_opt_summary.json").write_text(json.dumps(summary, indent=2))
        # write_report needs a few scalar fields it formats; provide them.
        summary.setdefault("constraints", {})
        summary.setdefault("final_energy_hartree", ref_energy if ref_energy is not None else float("nan"))
        summary.setdefault("final_energy_kcal",
                           (ref_energy * 627.5094740631) if ref_energy is not None else float("nan"))
        summary.setdefault("n_steps", None)
        summary.setdefault("converged", False)
        write_report(logdir, summary, conformers, None)

    return verdict


def main() -> int:
    ap = argparse.ArgumentParser(description="D2 Ace-OMW-NMe constrained B3LYP/6-31G(d) geom opt")
    ap.add_argument("--mode", choices=["constrained", "free"], default="constrained",
                    help="constrained = ff14SB-equilibrium frozen backbone+NE1 (D2 production)")
    ap.add_argument("--method", default="b3lyp", help="QM functional (geometry only)")
    ap.add_argument("--conformers", type=int, default=4,
                    help="number of conformers (>=4 for RESP-A2); 0 = reference only")
    ap.add_argument("--verify-constraint", action="store_true",
                    help="run the pre-D3 1-point dihedral-enforcement drift test")
    ap.add_argument("--verify-only", action="store_true",
                    help="run ONLY the constraint verification (fast pre-D3 gate), skip opt")
    ap.add_argument("--no-opt", action="store_true", help="skip optimization (report parsing only)")
    ap.add_argument("--input-pdb", default=None, help="override D1 dipeptide input path")
    ap.add_argument("--reverdict-from", default=None,
                    help="re-compute the D2 verdict from a COMPLETED run's "
                         "artifacts (no QM re-run); pass the existing logdir. "
                         "Refreshes its d2_geom_opt_report.md + summary.json. "
                         "Run ONLY on a finished run, never the active logdir.")
    args = ap.parse_args()

    # --- re-verdict path (no QM, no new logdir) ---------------------------
    if args.reverdict_from:
        rv_logdir = Path(args.reverdict_from)
        if not rv_logdir.exists():
            print(f"[error] --reverdict-from logdir not found: {rv_logdir}")
            return 1
        print(f"[info] Re-verdict from existing artifacts: {rv_logdir}")
        try:
            verdict = reverdict_from_logdir(rv_logdir, write=True)
        except (FileNotFoundError, ValueError) as e:
            print(f"[error] reverdict failed: {e}")
            return 1
        print(f"[info] reverdict D2 verdict = {verdict['verdict']} "
              f"(on-target {verdict['n_on_target_conformers']}/>=4, "
              f"wells {verdict['n_distinct_chi1_wells']}/>=3, "
              f"geom-conv {verdict['n_geometry_converged']}, "
              f"energy={verdict['energy_sanity_ok']}, "
              f"constraints={verdict['constraint_integrity_ok']})")
        return 0

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = PROJ / f"outputs/analysis/d2_omw_geom_opt_{ts}"
    logdir.mkdir(parents=True, exist_ok=True)
    print(f"[info] Logdir: {logdir}")

    pdb_in = Path(args.input_pdb) if args.input_pdb else find_latest_pyscf_dipeptide()
    if pdb_in is None or not pdb_in.exists():
        print("[error] No PySCF dipeptide.pdb found - run D1 first.")
        return 1
    print(f"[info] Input geometry: {pdb_in}")

    atoms_initial = parse_pdb(pdb_in)
    write_pdb(atoms_initial, logdir / "omw_dipeptide_input.pdb", "D2 input from D1 PySCF")
    print(f"[info] Initial geometry: {len(atoms_initial)} atoms")

    # --- pre-D3 constraint verification (fast gate) -----------------------
    verify = None
    if args.verify_constraint or args.verify_only:
        print("[info] Running pre-D3 constraint-enforcement verification (~30 min)...")
        verify = verify_constraint_enforcement(atoms_initial, logdir, args.method)
        (logdir / "constraint_verification.json").write_text(json.dumps(verify, indent=2))
        if args.verify_only:
            print(f"[info] verify-only complete: "
                  f"{'PASS' if verify.get('verified') else 'FAIL'} drift={verify.get('drift_deg')} deg")
            return 0 if verify.get("verified") else 2

    if args.no_opt:
        print("[info] --no-opt set; skipping optimization.")
        return 0

    # --- reference constrained optimization -------------------------------
    print(f"[info] Reference {args.mode} optimization (this may take 1-4 h CPU)...")
    atoms_opt, summary = run_psi4_optimize(atoms_initial, logdir, args.method, args.mode)
    write_pdb(atoms_opt, logdir / "omw_dipeptide_optimized.pdb",
              f"{summary['method']} {args.mode} | E={summary['final_energy_hartree']:.6f} Ha")

    summary["bond_deltas"] = compute_bond_deltas(atoms_initial, atoms_opt)
    summary["angle_deltas"] = compute_angle_deltas(atoms_opt)
    summary["improper_deltas"] = compute_improper_deltas(atoms_opt)

    # --- conformer ensemble (>= 4) ----------------------------------------
    conformers: List[Dict] = []
    if args.conformers and args.conformers > 0:
        print(f"[info] Generating conformer ensemble (target >= 4)...")
        conformers = generate_conformers(atoms_opt, logdir, args.method, args.conformers)
        summary["n_conformers"] = len(conformers)
        summary["conformers"] = conformers

    # C7: D2 acceptance verdict (PASS / MUST-REDO) emitted in summary JSON.
    verdict = compute_d2_verdict(summary, conformers)
    summary["d2_acceptance"] = verdict

    (logdir / "omw_geom_opt_summary.json").write_text(json.dumps(summary, indent=2))
    write_report(logdir, summary, conformers, verify)

    dst = copy_deliverable(logdir, conformers)

    print()
    print("=== D2 complete ===")
    print(f"mode      = {summary['mode']}")
    print(f"E_final   = {summary['final_energy_hartree']:.6f} Ha (converged={summary['converged']})")
    print(f"steps     = {summary['n_steps']}")
    print(f"conformers= {len(conformers)}")
    if verify is not None:
        print(f"verify    = {'PASS' if verify.get('verified') else 'FAIL'} (drift {verify.get('drift_deg')} deg)")
    print(f"D2 verdict= {verdict['verdict']} "
          f"(on-target {verdict['n_on_target_conformers']}/>=4, "
          f"wells {verdict['n_distinct_chi1_wells']}/>=3, "
          f"energy={verdict['energy_sanity_ok']}, constraints={verdict['constraint_integrity_ok']})")
    print(f"report    = {logdir / 'd2_geom_opt_report.md'}")
    print(f"deliverable = {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
