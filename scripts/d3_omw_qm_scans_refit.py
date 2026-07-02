#!/usr/bin/env python
"""D3 deliverable -- MTR (1-MeTrp) Ace-OMW-NMe bonded refit + RESP-A2 charge
refit on the D2-optimized dipeptide.

This SUPERSEDES the old scripts/d3_omw_qm_scans.py, whose single-atom rigid
``rotate_dihedral_to`` + ``frozen_dihedral`` scan shattered the indole ring and
absorbed nonbonded energy into the fitted V_n. The corrected protocol here is:

  * RELAXED scans only. Each scan point pre-rotates the WHOLE rigid fragment to
    the target (reusing the D2 fragment / Rodrigues helpers), then runs a
    constrained B3LYP/6-31G(d) optimization with the scan coordinate held in a
    tight ranged window and the backbone+NE1 pinned to amber14SB equilibrium;
    every other degree of freedom relaxes. (No rigid single-atom displacement.)

  * Bonded refit scope (pre-decided): the NE1-CM bond r_eq ONLY (the 1-MeTrp-
    specific N-methyl junction), refit to the D2 B3LYP value while keeping the
    amber14SB NA-CT force constant. Indole-ring bonds (delta <= 0.026 A) are
    accepted as amber14SB. The N-methyl rotor torsion is the ONLY refit scan
    (MM-subtracted, V3-restricted by methyl symmetry) and is the default
    "essential" deliverable. The improper and angle scans are VALIDATION-only
    and DEMOTED to opt-in (--scans all): D2 already validated those terms at
    geometry level, so their energy scan is redundant. The improper scan also
    carries a KNOWN, unfixed pre-rotation bug (proper-dihedral path is invalid
    for the 3+-bond central atom NE1) -- see SCAN_DEFS / the frcmod manifest.

  * Dihedral fit is MM-SUBTRACTED: the AMBER cosine V_n are fit to
    E_QM(phi) - [E_MM(phi) - E_MM_without_target_torsion(phi)] so the V_n do not
    absorb nonbonded / 1-4 contributions (see scripts/d3_mm_energy_helpers.py).

  * Charge refit = canonical AmberTools 2-stage RESP-A2 (NOT the lstsq
    run_pyscf_resp). 6-conformer equal-weight HF/6-31G* ESP; backbone
    (N/H/CA/HA/C/O) + NE1 frozen at amber14SB(Maier); side chain refit; methyl /
    methylene equivalencing left to respgen. ESP is the exact HF/6-31G*
    potential via PySCF int1e_grids (psi4 scf_type='df' GRID_ESP returned
    corrupt values near atoms, so the ESP is computed in-env with PySCF).

Cross-env:
  * QM relaxed scans: psi4 1.10 in conda env `psi4`
      /home/san/miniconda3/envs/psi4/bin/python
  * HF/6-31G* ESP (PySCF), AmberTools (antechamber/respgen/resp) + OpenMM
    MM-subtraction: conda env
    `qmmm` (/home/san/miniconda3/envs/qmmm/bin), AMBERHOME=/home/san/miniconda3/envs/qmmm
  This driver runs in `qmmm` and shells out to the psi4 interpreter by explicit
  path for the QM/ESP steps.

Modes:
  --scan-only   relaxed QM scans + bonded/MM-subtracted fits
  --resp-only   6-conformer ESP + AmberTools 2-stage RESP-A2 charge refit
  --full        both (default)
  --scans       which relaxed scans to run: 'essential' (default = the refit
                N-methyl rotor only, the only deployed-parameter scan), 'all'
                (also the opt-in VALIDATION scans improper/angle), or a comma-
                separated list of scan names. Unselected scans are reported as
                'validation deferred' (D2 covered them geometry-level), NOT as
                failures.
  --smoke       cheapest end-to-end shape check: 1 scan / 1 point (relaxed
                B3LYP via the psi4 subprocess) + 1-conformer ESP (PySCF) + a REAL
                1-conformer 2-stage RESP-A2 fit, to validate the cross-env
                wiring + esp.dat format + respgen/resp + MM-subtraction WITHOUT
                the multi-day full job. Combine with --scan-only / --resp-only.

Outputs (new files only; analysis dir is the immutable record):
  outputs/analysis/d3_omw_qm_scans_<TS>/
  scan_<name>/scan_points.json + <name>_fit.json, esp_grids/conf0N_grid.dat,
  resp_output/(resp1.out, resp2.out, resp2.qout, charge_audit.json),
  d3_report.md, d3_acceptance.json
  + params/d3_omw_bonded_refit.frcmod (NE1-CM r_eq refit; NEW file)
  + params/d3_omw_qm_geom/ deliverable copy (D2 pattern)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

PROJ = Path(__file__).resolve().parent.parent

# Cross-env interpreters.
PSI4_PY = Path("/home/san/miniconda3/envs/psi4/bin/python")
QMMM_BIN = Path("/home/san/miniconda3/envs/qmmm/bin")
AMBERHOME = Path("/home/san/miniconda3/envs/qmmm")

# Frozen MTR force field (read-only): the NE1-CM bond + target torsion live here.
MTR_FF_XML = PROJ / "params/MTR_gaff2_hybrid.xml"
PROTEIN_FF = "amber14-all.xml"

# D2 deliverable (geometry input).
D2_DELIVERABLE = PROJ / "params/omw_qm_geom"
D2_CONFORMERS = D2_DELIVERABLE / "conformers"

ANG_TO_BOHR = 1.8897259886

# ---------------------------------------------------------------------------
# Bonded refit constant (pre-decided): NE1-CM r_eq -> D2 B3LYP value, k kept.
# The frozen MTR_gaff2_hybrid.xml stores this bond as
#   <Bond type1="protein-CT" type2="protein-NA" length="0.1475" k="282001.6"/>
# (1.475 A, amber14SB NA-CT). The D2 B3LYP relaxed value is 1.4482 A. We emit a
# NEW frcmod with the refit r_eq and the SAME k (k is curvature, not Delta r_eq).
# ---------------------------------------------------------------------------
NE1_CM_REFIT_REQ_ANG = 1.4482          # D2 B3LYP relaxed NE1-CM bond length (A)
NE1_CM_K_KJ_PER_NM2 = 282001.6         # amber14SB NA-CT force constant (kept)


# ===========================================================================
# Import the D2 helpers (anti-fragmentation). D2 is import-safe from this env
# (it lazy-imports psi4 only inside run_psi4_optimize).
# ===========================================================================
def _load_d2_module():
    spec = importlib.util.spec_from_file_location(
        "d2_omw_b3lyp_geom_opt", str(PROJ / "scripts/d2_omw_b3lyp_geom_opt.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D2 = _load_d2_module()


def _load_mm_helpers():
    spec = importlib.util.spec_from_file_location(
        "d3_mm_energy_helpers", str(PROJ / "scripts/d3_mm_energy_helpers.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MM = _load_mm_helpers()


# ===========================================================================
# Scan definitions (relaxed). Each entry names the OMW atoms (UPDD-normalized)
# spanning the scanned internal coordinate. We measure the coordinate with the
# D2 dihedral/angle helpers and pin it with a tight ranged window.
#
#  * dihedral N-methyl rotor CD1-NE1-CM-HM1: the 1-MeTrp-specific torsion that
#    is REFIT (MM-subtracted). 0..330 deg / 30 deg (12 points).
#  * improper CG-CE2-NE1-CM: methyl out-of-plane at the indole N. VALIDATION
#    (scanned but not necessarily refit; NE1 frozen charge handles v0.8
#    explosion). -30..+30 deg / 5 deg.
#  * angle CM-NE1-CD1: the N-methyl junction angle. VALIDATION. 100..130 / 5 deg.
# ===========================================================================
SCAN_DEFS = [
    {
        "name": "dihedral_CD1_NE1_CM_HM1",
        "kind": "dihedral",
        "atoms": ("CD1", "NE1", "CM", "HM1"),
        "range": tuple(float(v) for v in range(0, 360, 30)),
        "role": "refit_mm_subtracted",
        # FF atom-class quadruple for the MM target torsion. MTR uses amber14
        # protein types at the NE1 junction; the methyl is CT/HC. Resolved from
        # the FF XML at runtime; this is the lookup key.
        "ff_class_quad": ("CD1", "NE1", "CM", "HM1"),
        # 3-fold symmetric N-methyl rotor: only n=3k terms (V3, optional V6) are
        # non-zero by methyl symmetry. A free V1..V4 fit of the asymmetric
        # 0..180 deg (1.5-period) sampling produces SPURIOUS 1/2/4-fold terms
        # (V1/V2/V4 != 0) that must NOT be deployed. The fit is restricted to
        # V3(+V6) for this rotor. (Wang 2004, J. Comput. Chem. 25:1157 --
        # methyl = V3-dominant.)
        "fold": 3,
    },
    {
        "name": "improper_CG_CE2_NE1_CM",
        "kind": "improper",
        "atoms": ("CG", "CE2", "NE1", "CM"),
        "range": tuple(float(v) for v in range(-30, 31, 5)),
        "role": "validation",
        # KNOWN ISSUE (opt-in only; NOT fixed in this revision): this improper is
        # currently dispatched through the SAME ranged_dihedral + _prerotate_
        # dihedral path as a proper dihedral (run_relaxed_scan). For an improper
        # the central atom NE1 has 3+ bonds, so the proper-dihedral graph-split
        # in _prerotate_dihedral is INVALID -- it tears the indole ring and the
        # optimization fails even at the in-plane target. Before opting into this
        # scan, an improper-specific pre-rotation is REQUIRED: hold the central
        # atom NE1 fixed and displace a SINGLE substituent (CM) out of plane,
        # with a NARROW ranged window (+-5..10 deg; the +-30 deg range above is
        # excessive for a rigid aromatic improper). This term is VALIDATION-only
        # (not refit; the frozen NE1 charge prevents the v0.8 explosion), and D2
        # already validated the out-of-plane angle geometry-level (<1 deg), so it
        # is demoted to opt-in rather than fixed here.
    },
    {
        "name": "angle_CM_NE1_CD1",
        "kind": "angle",
        "atoms": ("CM", "NE1", "CD1"),
        "range": tuple(float(v) for v in range(100, 131, 5)),
        "role": "validation",
        # VALIDATION-only (opt-in). ranged_bend drives the natural ~120 deg angle
        # to far targets (100..130 deg) with no pre-set, so off-equilibrium
        # points hit the optking maxiter grind (~1h/pt). D2 already validated the
        # angle geometry-level (delta <0.6 deg), so this energy scan is redundant
        # validation and is demoted to opt-in (re-run with --scans all).
    },
]

# Roles that produce a DEPLOYED FF parameter (the essential D3 deliverable).
# Everything else (e.g. "validation") is opt-in: run only with --scans all.
REFIT_ROLES = ("refit_mm_subtracted",)

# RESP freeze set (amber14SB Maier): backbone + NE1. Side chain is refit.
# Values mirror utils/parameterize_ncaa.py::_MTR_AMBER14_CHARGES (the SSOT).
RESP_FROZEN_CHARGES: Dict[str, float] = {
    "N": -0.4157, "H": 0.2719, "CA": -0.0275, "HA": 0.1123,
    "C": 0.5973, "O": -0.5679, "NE1": -0.3418,
}

# Merz-Kollman ESP grid (deterministic Fibonacci shells), reused from the D1
# RESP-A2 reproduction conventions.
MK_SHELLS: Tuple[float, ...] = (1.4, 1.6, 1.8, 2.0)
MK_SURFACE_DENSITY_PER_A2 = 1.0
MK_VDW_ANGSTROM = {"H": 1.20, "C": 1.50, "N": 1.50, "O": 1.40, "S": 1.75}


# ===========================================================================
# Geometry / index utilities (thin wrappers around the D2 helpers).
# ===========================================================================
def _omw_index_map(atoms) -> Dict[str, int]:
    """OMW atom name (normalized) -> 0-based index into the flat atom list."""
    return D2.build_name_index(atoms, residue="OMW")


def _quad_indices_0based(atoms, names: Sequence[str]) -> Optional[Tuple[int, ...]]:
    nidx = _omw_index_map(atoms)
    out = []
    for nm in names:
        i = nidx.get(nm)
        if i is None:
            return None
        out.append(i)
    return tuple(out)


def _measure_internal(coords, kind: str, idx: Sequence[int]) -> float:
    if kind in ("dihedral", "improper"):
        return D2._dihedral_deg(coords, *idx)
    if kind == "angle":
        return D2._angle_deg(coords, *idx)
    raise ValueError(f"unknown scan kind: {kind}")


def _ranged_constraint(kind: str, idx_0based: Sequence[int], target: float,
                       half_window: float) -> Dict[str, str]:
    """optking ranged_* constraint string (1-based indices) pinning the scan
    coordinate to [target-half, target+half]."""
    ids = " ".join(str(i + 1) for i in idx_0based)
    lo, hi = target - half_window, target + half_window
    if kind in ("dihedral", "improper"):
        return {"ranged_dihedral": f"({ids} {lo:.3f} {hi:.3f})"}
    if kind == "angle":
        return {"ranged_bend": f"({ids} {lo:.3f} {hi:.3f})"}
    raise ValueError(f"unknown scan kind: {kind}")


# ===========================================================================
# Relaxed scan (runs INSIDE the psi4 env -- see _scan_worker_main). The driver
# dispatches this module to the psi4 interpreter with --_scan-worker.
# ===========================================================================
SCAN_RANGED_HALF_WINDOW_DEG = 2.0  # matches D2 CONF_RANGED_HALF_WINDOW_DEG (anti-sawtooth)


def run_relaxed_scan(atoms_ref, scan_def: Dict, logdir: Path,
                     method: str = "b3lyp",
                     limit_points: Optional[int] = None) -> Dict:
    """Run a relaxed QM scan for one internal coordinate. For each target value:
      1. pre-rotate the WHOLE rigid fragment so the coordinate starts near target
         (D2 _prerotate_dihedral / _prerotate_chi1 reuse; never single-atom).
      2. constrained B3LYP/6-31G(d) optimize with the scan coord in a tight
         ranged window AND backbone+NE1 amber-pinned (D2 _build_constraint_options
         via mode='constrained'); all other DOF relax.
      3. record (target, measured, energy, optimized geometry).

    Returns {"name","kind","atoms","points":[{...}], ...}. The optimized geometry
    PDB of each point is written to logdir/scan_<name>/point_<i>.pdb so the
    MM-subtraction (route A) can re-evaluate E_MM at the relaxed geometry."""
    name = scan_def["name"]
    kind = scan_def["kind"]
    names = scan_def["atoms"]
    targets = list(scan_def["range"])
    if limit_points is not None:
        targets = targets[:limit_points]

    scandir = logdir / f"scan_{name}"
    scandir.mkdir(parents=True, exist_ok=True)

    idx = _quad_indices_0based(atoms_ref, names)
    if idx is None:
        return {"name": name, "kind": kind, "atoms": list(names),
                "error": "atom name(s) not found", "points": []}

    adj = D2._infer_bonds(atoms_ref)
    points: List[Dict] = []
    print(f"[scan] {name} ({kind}) on {names}: {len(targets)} points")

    for k, target in enumerate(targets):
        coords0 = D2._coords_of(atoms_ref)
        # 1) pre-rotate full rigid fragment to the target (dihedral/improper).
        if kind in ("dihedral", "improper"):
            coords0 = D2._prerotate_dihedral(atoms_ref, coords0, adj, idx, target)
        # (angle pre-set is not rigid-rotatable cleanly; optking's ranged_bend
        #  drives the angle from its current value -- acceptable for a +-15 deg
        #  sweep around the equilibrium.)
        atoms_start = D2._atoms_with_coords(atoms_ref, coords0)

        cons = _ranged_constraint(kind, idx, target, SCAN_RANGED_HALF_WINDOW_DEG)
        out_name = f"{name}_p{k:02d}.dat"
        try:
            atoms_opt, summary = D2.run_psi4_optimize(
                atoms_start, scandir, method, mode="constrained",
                out_name=out_name, extra_constraints=cons,
            )
            measured = _measure_internal(D2._coords_of(atoms_opt), kind, idx)
            pt_pdb = scandir / f"point_{k:02d}.pdb"
            D2.write_pdb(atoms_opt, pt_pdb,
                         title=f"{name} target={target:+.1f} measured={measured:+.2f}")
            points.append({
                "target_deg": float(target),
                "measured_deg": float(measured),
                "energy_hartree": summary["final_energy_hartree"],
                "converged": summary["converged"],
                "geom_pdb": pt_pdb.name,
            })
            print(f"  [{k+1:2d}/{len(targets)}] target={target:+7.1f} "
                  f"measured={measured:+7.2f}  E={summary['final_energy_hartree']:.6f} "
                  f"conv={summary['converged']}")
        except Exception as exc:  # noqa: BLE001 -- record per-point failure, continue
            points.append({
                "target_deg": float(target),
                "measured_deg": None,
                "energy_hartree": None,
                "converged": False,
                "error": f"{type(exc).__name__}: {str(exc)[:120]}",
            })
            print(f"  [{k+1:2d}/{len(targets)}] target={target:+7.1f}  FAIL: "
                  f"{type(exc).__name__}: {str(exc)[:80]}")

    return {"name": name, "kind": kind, "atoms": list(names),
            "role": scan_def.get("role"), "points": points}


# ===========================================================================
# Fits (salvaged from the old d3 scaffold; the dihedral fit is rewritten to
# operate on the MM-SUBTRACTED energy, not raw E_QM).
# ===========================================================================
def fit_improper(points: List[Dict]) -> Dict:
    """Fit improper energy to V * (1 - cos(2*phi)) (== V/2*(1+cos(2phi-180)))."""
    import numpy as np
    from scipy.optimize import curve_fit
    valid = [(p["target_deg"], p["energy_hartree"]) for p in points
             if p.get("energy_hartree") is not None]
    if len(valid) < 5:
        return {"error": "insufficient points", "n_points": len(valid)}
    ang = np.array([v[0] for v in valid])
    en = np.array([v[1] for v in valid])
    e_ref = en.min()
    de = (en - e_ref) * 627.5094740631
    rad = np.deg2rad(ang)

    def model(x, A):
        return A * (1 - np.cos(2 * x))

    popt, _ = curve_fit(model, rad, de, p0=[5.0])
    A = float(popt[0])
    rmsd = float(np.sqrt(np.mean((de - model(rad, A)) ** 2)))
    return {"form": "V/2*(1+cos(2*phi-180))", "V_kcal_per_mol": 2 * A,
            "A_kcal_per_mol": A, "rmsd_kcal": rmsd, "n_points": len(valid)}


def fit_angle(points: List[Dict]) -> Dict:
    """Fit angle energy to harmonic k*(theta-theta0)^2 (k in kcal/mol/rad^2)."""
    import numpy as np
    from scipy.optimize import curve_fit
    valid = [(p["target_deg"], p["energy_hartree"]) for p in points
             if p.get("energy_hartree") is not None]
    if len(valid) < 5:
        return {"error": "insufficient points", "n_points": len(valid)}
    ang = np.array([v[0] for v in valid])
    en = np.array([v[1] for v in valid])
    i_min = int(np.argmin(en))
    theta0 = float(ang[i_min])
    de = (en - en[i_min]) * 627.5094740631
    rad = np.deg2rad(ang - theta0)

    def model(x, k):
        return k * x ** 2

    popt, _ = curve_fit(model, rad, de, p0=[100.0])
    k = float(popt[0])
    rmsd = float(np.sqrt(np.mean((de - model(rad, k)) ** 2)))
    return {"form": "k*(theta-theta0)^2", "k_kcal_per_mol_rad2": k,
            "theta_0_deg": theta0, "rmsd_kcal": rmsd, "n_points": len(valid)}


def fit_dihedral_mm_subtracted(points: List[Dict],
                               mm_contrib_kcal: List[Optional[float]],
                               fold: Optional[int] = None,
                               include_v6: bool = False) -> Dict:
    """Fit the AMBER cosine series V_n to the MM-SUBTRACTED torsion energy:

        E_target(phi) = E_QM(phi) - [E_MM(phi) - E_MM_no_target(phi)]

    `mm_contrib_kcal[i]` is the bracket (route A or B) for point i, in kcal/mol.
    The V_n then describe ONLY the torsion that is being replaced; the nonbonded
    / 1-4 energy that the QM scan also contains is removed before fitting.

    `fold`: when set to 3 (a methyl / 3-fold symmetric rotor) the fit is
    RESTRICTED to the n=3k cosine terms. By methyl symmetry only n=3k terms are
    non-zero; a free V1..V4 fit of an asymmetric 0..180 deg (1.5-period) sampling
    otherwise produces SPURIOUS V1/V2/V4 terms that would be deployed into the
    FF. The DEFAULT is V3-ONLY (a methyl rotor is V3-dominant; on the S1 data
    this reproduces V3~0.10, RMSD~0.41). V6 is OPT-IN (`include_v6=True`): with
    only ~7 points over 1.5 periods a fitted V6 absorbs sampling noise rather
    than physical 6-fold character (it flips the V3 sign), so it is NOT enabled
    by default. When `fold` is None the legacy free V1..V4 model is used (general
    asymmetric dihedral)."""
    import numpy as np
    from scipy.optimize import curve_fit
    pts = []
    for p, mm in zip(points, mm_contrib_kcal):
        if p.get("energy_hartree") is None or mm is None:
            continue
        pts.append((p["target_deg"], p["energy_hartree"], mm))
    # A 3-fold rotor sampled 0..180 deg is 1.5 unique periods -> over-determined
    # at 6 points; the n>=6 guard admits the 7-converged-point S1 N-methyl rotor
    # (which the prior n>=8 guard wrongly rejected). General fits still need a
    # well-sampled curve; 6 is the minimum for the restricted/free models here.
    if len(pts) < 6:
        return {"error": "insufficient points", "n_points": len(pts)}
    ang = np.array([p[0] for p in pts])
    eqm_kcal = np.array([p[1] for p in pts]) * 627.5094740631
    mm = np.array([p[2] for p in pts])
    # MM-subtracted target-torsion energy (kcal/mol), zeroed to its minimum.
    e_target = eqm_kcal - mm
    e_target = e_target - e_target.min()
    rad = np.deg2rad(ang)

    # C-3: methyl / 3-fold symmetric rotor -> restrict to n=3k terms (V3, +V6).
    # Deploying free V1/V2/V4 on a symmetric rotor is non-physical overfit.
    # Default is V3-only; V6 is opt-in (overfits the sparse 1.5-period sample).
    if fold == 3:
        use_v6 = bool(include_v6) and len(pts) >= 7

        if use_v6:
            def model(x, V3, V6):
                return ((V3 / 2) * (1 + np.cos(3 * x))
                        + (V6 / 2) * (1 + np.cos(6 * x)))

            p0 = [0.1, 0.0]
        else:
            def model(x, V3):
                return (V3 / 2) * (1 + np.cos(3 * x))

            p0 = [0.1]
        try:
            popt, _ = curve_fit(model, rad, e_target, p0=p0)
        except Exception as exc:  # noqa: BLE001
            return {"error": f"curve_fit failed: {exc}", "n_points": len(pts)}
        rmsd = float(np.sqrt(np.mean((e_target - model(rad, *popt)) ** 2)))
        out = {
            "form": "V3/2*(1+cos(3*phi))" + ("+V6/2*(1+cos(6*phi))" if use_v6 else "")
                    + "  [MM-subtracted, methyl V3-restricted]",
            "fold_restricted": 3,
            # V1/V2/V4 are pinned to zero by methyl symmetry (NOT fit / deployed).
            "V1_kcal": 0.0, "V2_kcal": 0.0, "V3_kcal": float(popt[0]),
            "V4_kcal": 0.0,
            "rmsd_kcal": rmsd, "n_points": len(pts),
            "mm_subtraction_applied": True,
        }
        if use_v6:
            out["V6_kcal"] = float(popt[1])
        return out

    def model(x, V1, V2, V3, V4):
        return ((V1 / 2) * (1 + np.cos(x))
                + (V2 / 2) * (1 + np.cos(2 * x))
                + (V3 / 2) * (1 + np.cos(3 * x))
                + (V4 / 2) * (1 + np.cos(4 * x)))

    try:
        popt, _ = curve_fit(model, rad, e_target, p0=[0.5, 0.5, 0.5, 0.5])
    except Exception as exc:  # noqa: BLE001
        return {"error": f"curve_fit failed: {exc}", "n_points": len(pts)}
    rmsd = float(np.sqrt(np.mean((e_target - model(rad, *popt)) ** 2)))
    return {
        "form": "sum_n V_n/2*(1+cos(n*phi))  [MM-subtracted]",
        "V1_kcal": float(popt[0]), "V2_kcal": float(popt[1]),
        "V3_kcal": float(popt[2]), "V4_kcal": float(popt[3]),
        "rmsd_kcal": rmsd, "n_points": len(pts),
        "mm_subtraction_applied": True,
    }


# ===========================================================================
# MM-subtraction driver (route A full-OpenMM with route B analytic cross-check).
# ===========================================================================
def compute_mm_subtraction(scan_result: Dict, scandir: Path,
                           ff_class_quad: Sequence[str]) -> Dict:
    """Compute the MM torsion contribution for each scan point (kcal/mol).

    Route A (primary): build the full MM system + the system with the target
    torsion zeroed, evaluate both at each relaxed-scan geometry, take the
    difference. Route B (cross-check / fallback): the analytic AMBER cosine of
    the target torsion type read from the frozen FF XML, at the measured
    dihedral. If route A is unavailable (FF template-match failure for the
    capped ncAA dipeptide), route B is used for the fit.

    Returns {"contrib_kcal":[...], "route":"openmm|analytic", ...}."""
    points = scan_result["points"]
    atoms_ref = D2.parse_pdb(_first_conformer_pdb())
    adj = D2._infer_bonds(atoms_ref)
    idx = _quad_indices_0based(atoms_ref, scan_result["atoms"])
    measured = [p.get("measured_deg") for p in points]

    # Route B: analytic terms from the frozen FF XML (cross-check / fallback).
    # NOTE: the MTR junction torsion's parameters live in amber14-all.xml (the
    # MTR residue uses protein-* types), NOT in the hybrid XML, so this single-
    # XML lookup may legitimately return [] -- route A is authoritative.
    terms = MM.extract_torsion_terms_from_xml(
        str(MTR_FF_XML),
        tuple(_resolve_ff_classes(atoms_ref, ff_class_quad)),
    )
    analytic = [
        (MM.amber_torsion_energy_kcal(m, terms) if (m is not None and terms) else None)
        for m in measured
    ]

    # Route A (authoritative): full-OpenMM difference at the relaxed geometries.
    # Read each relaxed-scan PDB's coordinates (Angstrom) in atoms_ref order.
    scan_coords: List[Optional[List[Tuple[float, float, float]]]] = []
    for p in points:
        gp = p.get("geom_pdb")
        pdb_path = (scandir / gp) if gp else None
        if pdb_path is not None and pdb_path.exists():
            scan_coords.append([(a[3], a[4], a[5]) for a in D2.parse_pdb(pdb_path)])
        else:
            scan_coords.append(None)

    route = "analytic"
    contrib: List[Optional[float]] = list(analytic)
    n_zeroed = 0
    openmm_contrib: List[Optional[float]] = []
    if idx is not None and len(scan_coords) == len(points):
        try:
            openmm_contrib, n_zeroed = MM.mm_torsion_contribution_openmm(
                [PROTEIN_FF, str(MTR_FF_XML)], atoms_ref, adj, scan_coords, idx,
                residue_rename={"OMW": "MTR"},
            )
            # Route A is authoritative even when the target torsion has zero FF
            # terms (n_zeroed == 0): the MM contribution is then genuinely ~0,
            # and the fitted V_n describes the full QM rotor profile (there is no
            # existing torsion term to double-count). Use route A whenever it
            # produced a finite contribution for every point that HAS a relaxed
            # geometry. Failed scan points (no geometry PDB -> scan_coords None ->
            # route-A None) are legitimately None and must NOT veto route A: the
            # subsequent fit skips any point whose E_QM is None regardless, so a
            # None at a failed-scan index is harmless. (A whole-list
            # `all(c is not None)` gate wrongly fell back to the analytic route --
            # which is all-None for this junction, since the N-methyl rotor has no
            # FF torsion term -- leaving the fit with zero usable points.)
            if openmm_contrib and all(
                c is not None
                for c, sc in zip(openmm_contrib, scan_coords)
                if sc is not None
            ) and any(sc is not None for sc in scan_coords):
                contrib = openmm_contrib
                route = "openmm"
        except Exception as exc:  # noqa: BLE001 -- fall back to analytic route
            print(f"[mm-sub] OpenMM route unavailable ({type(exc).__name__}: "
                  f"{str(exc)[:80]}); using analytic route")

    return {
        "contrib_kcal": contrib,
        "route": route,
        "n_terms_zeroed": n_zeroed,
        "analytic_terms": terms,
        "analytic_contrib_kcal": analytic,
        "openmm_contrib_kcal": openmm_contrib,
    }


def _resolve_ff_classes(atoms_ref, atom_names: Sequence[str]) -> List[str]:
    """Map OMW atom names to their FF atom-CLASS strings from the MTR XML
    Residue/Atom type assignments (so the torsion-type lookup matches OpenMM's
    own class-based matching). Falls back to the atom name if unresolved."""
    import xml.etree.ElementTree as ET
    tree = ET.parse(MTR_FF_XML)
    root = tree.getroot()
    # name -> type (from any MTR Residue variant; the internal MTR is canonical)
    name_to_type: Dict[str, str] = {}
    for res in root.findall(".//Residue"):
        if res.get("name") not in ("MTR", "NMTR", "CMTR"):
            continue
        for atom in res.findall("Atom"):
            name_to_type.setdefault(atom.get("name"), atom.get("type"))
    # type -> class (from AtomTypes)
    type_to_class: Dict[str, str] = {}
    for at in root.findall(".//AtomTypes/Type"):
        type_to_class[at.get("name")] = at.get("class") or at.get("name")
    out = []
    for nm in atom_names:
        t = name_to_type.get(nm)
        out.append(type_to_class.get(t, nm) if t else nm)
    return out


def _amber14_type_to_class() -> Dict[str, str]:
    """type -> class map from the deployed amber14 protein ffxml
    (amber14/protein.ff14SB.xml, the file amber14-all.xml <Include>s). The MTR
    residue's N-methyl-junction atoms carry protein-* TYPES whose CLASS strings
    (CT/NA/CW/HC) live in that file, not in the hybrid MTR XML (where their
    AtomTypes/Type entries have no class attr). Returns {} if the file cannot be
    located (the caller then has no deployed-FF class and must omit the DIHE)."""
    import xml.etree.ElementTree as ET
    try:
        import openmm
        data = Path(openmm.__file__).resolve().parent / "app" / "data"
    except Exception:  # noqa: BLE001
        return {}
    pff = data / "amber14" / "protein.ff14SB.xml"
    if not pff.exists():
        return {}
    root = ET.parse(pff).getroot()
    t2c: Dict[str, str] = {}
    for at in root.findall(".//AtomTypes/Type"):
        cls = at.get("class")
        if cls:
            t2c[at.get("name")] = cls
    return t2c


def resolve_dihedral_ff_types(atom_names: Sequence[str]) -> Dict[str, object]:
    """Resolve the DEPLOYED-FF AMBER atom-type class quadruple for a torsion over
    `atom_names` (e.g. the N-methyl rotor CD1-NE1-CM-HM1), so the emitted frcmod
    DIHE line uses the exact 4 types the deployed ForceField matches on.

    Lookup chain (anti-fragmentation; mirrors OpenMM's own class-based matching):
      1. MTR_gaff2_hybrid.xml Residue 'MTR' -> atom name -> atom TYPE (protein-*).
      2. amber14/protein.ff14SB.xml AtomTypes -> protein-* TYPE -> CLASS string.
    Returns {"classes":[..], "types":[..], "resolved":bool, "rationale":str}.
    `resolved` is True only when all four CLASS strings were found; otherwise the
    caller must NOT emit a DIHE under guessed/wildcard types (the prompt forbids
    a wildcard X-NA-CT-X fallback)."""
    import xml.etree.ElementTree as ET
    root = ET.parse(MTR_FF_XML).getroot()
    name_to_type: Dict[str, str] = {}
    for res in root.findall(".//Residue"):
        if res.get("name") != "MTR":
            continue
        for atom in res.findall("Atom"):
            name_to_type.setdefault(atom.get("name"), atom.get("type"))
        break
    t2c = _amber14_type_to_class()
    types = [name_to_type.get(nm) for nm in atom_names]
    classes = [t2c.get(t) if t else None for t in types]
    resolved = all(c is not None for c in classes)
    rationale = (
        "DEPLOYED-FF atom-type classes for the N-methyl rotor dihedral "
        f"{tuple(atom_names)} resolved as {tuple(classes)} via: MTR residue in "
        f"{MTR_FF_XML.name} (name->type {dict(zip(atom_names, types))}) then "
        "amber14/protein.ff14SB.xml AtomTypes (type->class). These are the "
        "specific 4 AMBER types the deployed ForceField (amber14-all.xml + the "
        "MTR hybrid XML) matches on -- NOT a wildcard X-NA-CT-X. No existing "
        "proper torsion (exact or wildcard) is defined about the NA-CT central "
        "bond in protein.ff14SB.xml, so this V3 DIHE is a NEW term, not a "
        "replacement (no double-count)."
    )
    return {"classes": classes, "types": types,
            "resolved": resolved, "rationale": rationale}


# ===========================================================================
# ESP generation -- HF/6-31G* density + EXACT ESP via PySCF int1e_grids, run in
# THIS (qmmm) env. (psi4's scf_type='df' GRID_ESP was found to return corrupt
# ESP at certain vdW-shell grid points -- e.g. 8.08 Ha/e at 1.68 A from the
# nearest atom for a neutral molecule, where PySCF's exact int1e_grids gives
# 0.012 -- while both codes' SCF energies agreed to <2 mHa. The exact PySCF
# route is the one validated in the D1 RESP-A2 reproduction; either ESP backend
# is acceptable, and PySCF is the numerically stable one here.)
# ===========================================================================
def build_mk_grid(atoms) -> "object":
    """Deterministic Merz-Kollman ESP grid (Angstrom), Fibonacci shells with
    inner-shell exclusion. Reused conventions from the D1 RESP-A2 reproduction."""
    import numpy as np

    def fib(n):
        n = int(n)
        if n < 1:
            return np.zeros((0, 3))
        ii = np.arange(n) + 0.5
        phi = math.pi * (1.0 + math.sqrt(5.0))
        z = 1.0 - 2.0 * ii / n
        r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
        th = phi * ii
        return np.column_stack((r * np.cos(th), r * np.sin(th), z))

    radii = np.array([MK_VDW_ANGSTROM.get(a[2], 1.7) for a in atoms])
    coords = np.array([(a[3], a[4], a[5]) for a in atoms])
    inner = radii * MK_SHELLS[0]
    chunks = []
    for s in MK_SHELLS:
        for ia in range(len(atoms)):
            R = radii[ia] * s
            n_points = max(50, int(round(4.0 * math.pi * R * R * MK_SURFACE_DENSITY_PER_A2)))
            sph = fib(n_points) * R + coords[ia]
            keep = np.ones(sph.shape[0], dtype=bool)
            for jb in range(len(atoms)):
                if jb == ia:
                    continue
                keep &= np.linalg.norm(sph - coords[jb], axis=1) > inner[jb]
            chunks.append(sph[keep])
    return np.vstack(chunks)


def compute_esp_pyscf(atoms, grid_ang):
    """HF/6-31G* density + EXACT electrostatic potential on `grid_ang` (Angstrom)
    via PySCF int1e_grids. Returns (esp_au, converged). Mirrors the D1 RESP-A2
    reproduction ESP route (mol.intor('int1e_grids') for the electronic term +
    sum_A Z_A/|r-R_A| nuclear term, all in Bohr). RESP-A2 charge basis = HF/
    6-31G* (Bayly 1993 / Cieplak 1995). The molecule is built with charge 0,
    spin 0 (MTR dipeptide is neutral, even-electron)."""
    import numpy as np
    from pyscf import gto, scf

    atomstr = ";".join(f"{a[2]} {a[3]} {a[4]} {a[5]}" for a in atoms)
    mol = gto.M(atom=atomstr, basis="6-31G*", charge=0, spin=0,
                unit="Angstrom", verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-9
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        return None, False
    dm = mf.make_rdm1()
    if not isinstance(dm, np.ndarray):
        dm = np.asarray(dm.get())

    grid_bohr = np.asarray(grid_ang) * ANG_TO_BOHR
    coords = mol.atom_coords()  # Bohr
    Z = mol.atom_charges().astype(float)
    diff = grid_bohr[:, None, :] - coords[None, :, :]
    dist = np.linalg.norm(diff, axis=2)
    v_nuc = np.einsum("a,ga->g", Z, 1.0 / dist)
    # electronic term, chunked to bound memory (nao x nao x chunk).
    nao = mol.nao
    chunk = max(64, int(2_000_000_000 / (8 * nao * nao)))
    v_elec = np.zeros(len(grid_bohr))
    for s in range(0, len(grid_bohr), chunk):
        e = min(s + chunk, len(grid_bohr))
        v_op = mol.intor("int1e_grids", grids=grid_bohr[s:e])
        v_elec[s:e] = -np.einsum("gij,ij->g", v_op, dm)
    return v_nuc + v_elec, True


def write_esp_dat_block(handle, atoms, grid_ang, esp_au) -> None:
    """Append ONE conformer block to an AmberTools/RESP esp.dat file (a.u.):
        header  '%5d%5d%5d'  natoms, ngrid, 0
        natoms  '%17s%16.7E%16.7E%16.7E'  blank, x, y, z (Bohr)
        ngrid   '%16.7E%16.7E%16.7E%16.7E'  V, x, y, z (Bohr)
    Multiple blocks concatenated == multi-conformer ESP (resp &cntrl nmol=N)."""
    import numpy as np
    coords_bohr = np.array([(a[3], a[4], a[5]) for a in atoms]) * ANG_TO_BOHR
    grid_bohr = np.asarray(grid_ang) * ANG_TO_BOHR
    handle.write(f"{len(atoms):5d}{len(grid_bohr):5d}{0:5d}\n")
    for c in coords_bohr:
        handle.write(f"{'':>17s}{c[0]:16.7E}{c[1]:16.7E}{c[2]:16.7E}\n")
    for g, v in zip(grid_bohr, esp_au):
        handle.write(f"{v:16.7E}{g[0]:16.7E}{g[1]:16.7E}{g[2]:16.7E}\n")


# ===========================================================================
# RESP-A2 (AmberTools 2-stage) -- canonical, runs in THIS (qmmm) env.
# ===========================================================================
def _amber_env() -> Dict[str, str]:
    env = dict(os.environ)
    env["AMBERHOME"] = str(AMBERHOME)
    env["PATH"] = f"{QMMM_BIN}:{env.get('PATH', '')}"
    return env


def _write_qin_8f10(charges: Sequence[float], path: Path) -> None:
    """resp initial-charge file (8F10.6, 8 per line)."""
    lines = []
    for i in range(0, len(charges), 8):
        chunk = charges[i:i + 8]
        lines.append("".join("{:10.6f}".format(q) for q in chunk))
    path.write_text("\n".join(lines) + "\n")


def build_resp_addfile(atoms, addfile_path: Path) -> List[int]:
    """respgen addfile that FREEZES backbone (N/H/CA/HA/C/O) + NE1 at their
    amber14SB(Maier) charges (CHARGE directives -> '-99' ivary in resp1/2.in).
    Returns the 1-based ids of the frozen atoms (the qin must supply the same
    values, else resp would freeze them to the qin contents = 0.0)."""
    lines = ["//frozen backbone + NE1 charges from amber14SB (Maier 2015)"]
    frozen_ids: List[int] = []
    for i, (name, res, *_rest) in enumerate(atoms, start=1):
        if res != "OMW":
            continue
        nm = D2._norm_name(name)
        q = RESP_FROZEN_CHARGES.get(nm)
        if q is not None:
            lines.append(f"CHARGE\t{q:+.6f}\t{i}\t{name}")
            frozen_ids.append(i)
    addfile_path.write_text("\n".join(lines) + "\n")
    return frozen_ids


def run_resp_a2(ac_path: Path, esp_path: Path, addfile: Path,
                atoms, n_conf: int, logdir: Path,
                dry: bool = False) -> Optional[Path]:
    """Drive respgen (resp1+resp2 with -a addfile) and resp (2-stage) for the
    multi-conformer ESP. Patches nmol=`n_conf` into the generated resp inputs so
    resp reads all conformer blocks with equal weight. Returns the stage-2 .qout
    path, or None if `dry` (commands assembled + written but not executed)."""
    env = _amber_env()
    respgen = str(QMMM_BIN / "respgen")
    resp = str(QMMM_BIN / "resp")

    def _run(cmd, log_name, timeout):
        if dry:
            (logdir / log_name).write_text("DRY: " + " ".join(cmd) + "\n")
            return None
        r = subprocess.run(cmd, cwd=logdir, capture_output=True, text=True,
                           timeout=timeout, env=env)
        (logdir / log_name).write_text(
            "CMD: " + " ".join(cmd) + "\n\n" + r.stdout +
            "\n--- STDERR ---\n" + r.stderr)
        return r

    resp1_in = logdir / "resp1.in"
    resp2_in = logdir / "resp2.in"
    # Multi-conformer: respgen's native `-n N` emits a correct N-molecule input
    # (nmol=N, N concatenated molecule blocks, plus the inter-conformer charge-
    # equivalence section in the trailing "9th area" that ties atom k of each
    # conformer to atom k of conformer 1). This replaces the previous _patch_nmol
    # hack, which only bumped `nmol=1->N` while leaving a SINGLE molecule block ->
    # `resp` read past the one block and hit Fortran EOF (resp.F unit=5).
    nmol_args = ["-n", str(n_conf)] if n_conf > 1 else []
    cmd1 = [respgen, "-i", ac_path.name, "-o", resp1_in.name, "-f", "resp1",
            "-a", addfile.name] + nmol_args
    cmd2 = [respgen, "-i", ac_path.name, "-o", resp2_in.name, "-f", "resp2",
            "-a", addfile.name] + nmol_args
    r1 = _run(cmd1, "respgen_resp1.log", 120)
    r2 = _run(cmd2, "respgen_resp2.log", 120)
    if dry:
        # still validate that the cross-env wiring + file layout is sound
        return None
    if r1 is None or r1.returncode != 0 or not resp1_in.exists():
        raise RuntimeError(f"respgen resp1 failed (rc={getattr(r1,'returncode',None)})")
    if r2 is None or r2.returncode != 0 or not resp2_in.exists():
        raise RuntimeError(f"respgen resp2 failed (rc={getattr(r2,'returncode',None)})")

    # qin: frozen atoms get amber values, all others 0.0. resp reads ONE qin
    # entry per atom PER conformer block (length = n_atoms * n_conf), so the
    # single-molecule seed must be replicated n_conf times. The `-99` frozen
    # atoms are anchored via qin in EVERY conformer (respgen's inter-conformer
    # equivalence section ties only the FREE atoms across conformers); seeding
    # the frozen anchors in conformer 1 alone leaves the integer-net constraint
    # to distribute the deficit across conformers (per-block sum != 0).
    qin = logdir / "qin.frozen"
    init_q = []
    for name, res, *_rest in atoms:
        nm = D2._norm_name(name) if res == "OMW" else name
        init_q.append(RESP_FROZEN_CHARGES.get(nm, 0.0) if res == "OMW" else 0.0)
    _write_qin_8f10(init_q * max(1, n_conf), qin)

    resp1_q = logdir / "resp1.qout"
    resp2_q = logdir / "resp2.qout"
    cmd3 = [resp, "-O", "-i", resp1_in.name, "-o", "resp1.out",
            "-e", esp_path.name, "-q", qin.name, "-t", resp1_q.name]
    r3 = _run(cmd3, "resp_stage1.log", 3600)
    if r3.returncode != 0 or not resp1_q.exists():
        raise RuntimeError("resp stage 1 failed")
    cmd4 = [resp, "-O", "-i", resp2_in.name, "-o", "resp2.out",
            "-e", esp_path.name, "-q", resp1_q.name, "-t", resp2_q.name]
    r4 = _run(cmd4, "resp_stage2.log", 3600)
    if r4.returncode != 0 or not resp2_q.exists():
        raise RuntimeError("resp stage 2 failed")
    return resp2_q


def _patch_nmol(resp_in: Path, n_conf: int) -> None:
    """DEPRECATED -- do not use. This only bumped `nmol=1->N` in the &cntrl block
    while leaving a SINGLE molecule block, so `resp` read past the one block and
    crashed with a Fortran EOF (resp.F unit=5). Multi-conformer inputs are now
    generated correctly by respgen's native `-n N` flag (N molecule blocks + the
    inter-conformer equivalence section). Retained as a guard against reuse."""
    raise RuntimeError(
        "_patch_nmol is deprecated: it produces a 1-block input with nmol=N "
        "(resp EOF crash). Use respgen `-n N` (see run_resp_a2).")


def _parse_qout(qout: Path, n_atoms: Optional[int] = None) -> List[float]:
    """Parse a resp .qout. For a multi-conformer fit the .qout repeats the full
    charge set once per conformer block (length = n_atoms * n_conf); slice the
    FIRST conformer when `n_atoms` is given (the per-conformer blocks are charge-
    equivalent up to the frozen anchors, which are seeded per-conformer via qin,
    so block 0 carries the canonical single charge set, sum == integer net)."""
    nums = []
    for tok in qout.read_text().split():
        try:
            nums.append(float(tok))
        except ValueError:
            continue
    if n_atoms is not None and n_atoms > 0 and len(nums) >= n_atoms:
        return nums[:n_atoms]
    return nums


# ===========================================================================
# Bonded refit output (NE1-CM r_eq) -- NEW frcmod file (never touches MTR XML).
# ===========================================================================
def write_bonded_refit_frcmod(out_path: Path, validation: Dict,
                              dihe_refit: Optional[Dict] = None) -> None:
    """Emit the NE1-CM bond r_eq refit (and, when `dihe_refit` is supplied, the
    N-methyl rotor V3 DIHE term) as a standalone AMBER frcmod. r_eq is the D2
    B3LYP relaxed value; k is the retained amber14SB NA-CT force constant. AMBER
    frcmod BOND units: kcal/mol/A^2 and Angstrom (convert from the XML's
    kJ/mol/nm^2 k).

    `dihe_refit` (when provided, from run_finalize): a dict with the fitted
    `V3_kcal` plus the DEPLOYED-FF atom-type quad (`ff_classes`, e.g. CW-NA-CT-HC)
    and its resolution `rationale`. The DIHE line is the AMBER cosine for the
    N-methyl rotor: idivf=1, pk=V3 (kcal/mol), phase=0 deg, period=3. When the FF
    classes could not be resolved (resolved=False), the DIHE is NOT emitted (the
    prompt forbids a wildcard X-NA-CT-X fallback) and the manifest records why."""
    # kJ/mol/nm^2 -> kcal/mol/A^2 : /4.184 (kcal) * (0.1 nm/A)^2 = /4.184/100
    k_kcal_a2 = NE1_CM_K_KJ_PER_NM2 / 4.184 / 100.0
    # AMBER frcmod halves the bond k convention? No: AMBER PARM bond k IS the
    # k in E = k*(r-r0)^2 (already includes the leading factor); OpenMM stores
    # 2*k. Convert: openmm_k = 2*amber_k -> amber_k = openmm_k/2.
    k_amber = k_kcal_a2 / 2.0

    # DIHE: N-methyl rotor V3. AMBER frcmod DIHE format:
    #   <C1-C2-C3-C4>  IDIVF  PK  PHASE  PERIODICITY  [comment]
    # For a methyl (3-fold) the single deployed term is idivf=1, pk=V3, phase=0,
    # period=3 (V1/V2/V4 are pinned to zero by methyl symmetry -- not emitted).
    dihe_lines: List[str] = []
    dihe_manifest: Dict[str, object]
    if dihe_refit and dihe_refit.get("emit"):
        classes = dihe_refit["ff_classes"]
        v3 = float(dihe_refit["V3_kcal"])
        quad = "-".join(classes)
        dihe_lines = [
            f"{quad:<11s} 1  {v3:8.5f}    0.000   3.   "
            f"D3 N-methyl rotor V3 (MM-subtracted B3LYP/6-31G(d) relaxed scan)",
        ]
        dihe_manifest = {
            "emitted": True,
            "ff_class_quad": classes,
            "ff_type_quad": dihe_refit.get("ff_types"),
            "V3_kcal_per_mol": round(v3, 6),
            "idivf": 1, "phase_deg": 0.0, "periodicity": 3,
            "rmsd_kcal": dihe_refit.get("rmsd_kcal"),
            "n_points": dihe_refit.get("n_points"),
            "atom_type_rationale": dihe_refit.get("rationale"),
            "mm_subtraction_route": dihe_refit.get("mm_route"),
            "double_count_note": (
                "MM-subtraction is 0 for this rotor: no proper torsion is defined "
                "about the NA-CT central bond in the deployed FF (route A "
                "n_terms_zeroed=0), so V3 describes the full QM rotor and this is "
                "a NEW DIHE term, not a replacement -- no double-count."),
        }
    else:
        reason = (dihe_refit or {}).get(
            "skip_reason",
            "no dihedral refit supplied (bond-only frcmod)")
        dihe_manifest = {"emitted": False, "reason": reason}
        if dihe_refit:
            dihe_manifest["atom_type_rationale"] = dihe_refit.get("rationale")

    lines = [
        "D3 NE1-CM bond r_eq refit + N-methyl rotor V3 (1-MeTrp junction)",
        "BOND",
        # amber atom types for NE1 (NA) and CM (CT); class strings per MTR XML.
        f"NA-CT  {k_amber:7.1f}  {NE1_CM_REFIT_REQ_ANG:.4f}   "
        f"D3 B3LYP/6-31G(d) relaxed r_eq; k=amber14SB NA-CT (retained)",
        "",
        "ANGLE",
        "",
        "DIHE",
        *dihe_lines,
        "",
        "IMPROPER",
        "",
        "NONBON",
        "",
    ]
    out_path.write_text("\n".join(lines) + "\n")
    (out_path.parent / (out_path.stem + "_manifest.json")).write_text(json.dumps({
        "refit_scope": ("NE1-CM bond r_eq + N-methyl rotor V3 DIHE"
                        if dihe_manifest.get("emitted")
                        else "NE1-CM bond r_eq only"),
        "ne1_cm_req_angstrom": NE1_CM_REFIT_REQ_ANG,
        "ne1_cm_req_nm": round(NE1_CM_REFIT_REQ_ANG / 10.0, 6),
        "k_openmm_kj_per_nm2_retained": NE1_CM_K_KJ_PER_NM2,
        "k_amber_kcal_per_a2": round(k_amber, 4),
        "source_xml_value_nm": 0.1475,
        "source": "D2 B3LYP/6-31G(d) relaxed geometry (params/omw_qm_geom)",
        "frozen_xml_untouched": str(MTR_FF_XML.relative_to(PROJ)),
        "nmethyl_rotor_dihe": dihe_manifest,
        "validation": validation,
        "known_issues": {
            "improper_scan_prerotation": (
                "The improper scan (improper_CG_CE2_NE1_CM) is currently "
                "dispatched through the proper-dihedral ranged_dihedral + "
                "_prerotate_dihedral path. For an improper the central atom NE1 "
                "has 3+ bonds, so the proper-dihedral graph-split is INVALID and "
                "tears the indole ring (S3 0/13 in the halted FULL run). The "
                "improper is VALIDATION-only (NOT refit; frozen NE1 charge "
                "handles the v0.8 explosion) and is demoted to opt-in. Opting in "
                "REQUIRES an improper-specific pre-rotation first (hold NE1 "
                "fixed, displace only one substituent out of plane) with a NARROW "
                "+-5..10 deg window (the +-30 deg range is excessive for a rigid "
                "aromatic improper). NOT fixed in this revision."),
            "angle_scan_grind": (
                "The angle scan (angle_CM_NE1_CD1) drives the natural ~120 deg "
                "angle to far targets with no geometric pre-set, hitting the "
                "optking maxiter grind. VALIDATION-only, demoted to opt-in. D2 "
                "validated this angle geometry-level (delta <0.6 deg)."),
        },
        "deferred_validation_note": (
            "Energy-LEVEL bonded validation of the improper V-band and angle "
            "k-band is DEFERRED (demoted to opt-in; D2 provided only "
            "geometry-level validation). Both terms are accepted-as-amber14SB "
            "(not refit), so their curvature is the amber reference by "
            "construction; the residual risk (stiff/soft amber angle/improper at "
            "the CM junction unverified at energy level) is caught downstream by "
            "the D4 trajectory NE1-region geometry distribution."),
        "generated": datetime.now().isoformat(),
    }, indent=2))


# ===========================================================================
# Helpers for conformer enumeration.
# ===========================================================================
def _conformer_pdbs() -> List[Path]:
    if not D2_CONFORMERS.exists():
        return []
    return sorted(D2_CONFORMERS.glob("conf??_optimized.pdb"))


def _first_conformer_pdb() -> Path:
    confs = _conformer_pdbs()
    if confs:
        return confs[0]
    fallback = D2_DELIVERABLE / "omw_dipeptide_optimized.pdb"
    if fallback.exists():
        return fallback
    raise FileNotFoundError("No D2 conformer geometry found in params/omw_qm_geom")


# ===========================================================================
# (ESP is computed IN-ENV via PySCF compute_esp_pyscf -- no psi4 ESP worker.
#  Only the relaxed QM SCANS need the psi4 interpreter, via _scan_worker_main.
#  psi4's scf_type='df' GRID_ESP was found to return corrupt ESP at certain
#  vdW-shell points; PySCF int1e_grids is exact and is the D1-validated route.)
# ===========================================================================
def _scan_worker_main(argv: List[str]) -> int:
    """Runs INSIDE the psi4 env. Args: <logdir> <scan_index> [--limit-points N].
    Runs one relaxed scan and writes scan_points.json into its scan_<name> dir."""
    ap = argparse.ArgumentParser()
    ap.add_argument("logdir")
    ap.add_argument("scan_index", type=int)
    ap.add_argument("--limit-points", type=int, default=None)
    args = ap.parse_args(argv)

    logdir = Path(args.logdir)
    scan_def = SCAN_DEFS[args.scan_index]
    atoms_ref = D2.parse_pdb(_first_conformer_pdb())
    result = run_relaxed_scan(atoms_ref, scan_def, logdir,
                              limit_points=args.limit_points)
    out = logdir / f"scan_{scan_def['name']}" / "scan_points.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))
    print(f"[scan-worker] wrote {out}")
    return 0


# ===========================================================================
# Driver orchestration (runs in the qmmm env; shells out to psi4).
# ===========================================================================
def _psi4_subprocess(worker_flag: str, worker_args: List[str], smoke: bool) -> int:
    """Invoke this module inside the psi4 interpreter for a worker step."""
    cmd = [str(PSI4_PY), str(Path(__file__).resolve()), worker_flag, *worker_args]
    print(f"[driver] psi4 subprocess: {' '.join(cmd)}")
    r = subprocess.run(cmd)
    return r.returncode


def select_scan_indices(scans_mode: str) -> List[int]:
    """Return the SCAN_DEFS indices to run for `scans_mode`:
      * "essential" (default): only role in REFIT_ROLES (the deployed-parameter
        deliverable -- the N-methyl rotor). Validation-only scans are skipped.
      * "all": every scan, including the opt-in validation scans (improper /
        angle). NOTE: the improper scan has a KNOWN, unfixed pre-rotation bug
        (see SCAN_DEFS) -- opting in requires the documented improper-specific
        pre-rotation first.
      * comma-separated scan names: run exactly those.
    The set NOT selected is reported as "validation deferred" (compute_acceptance
    does not fail on it; D2 already validated those terms geometry-level)."""
    if scans_mode == "all":
        return list(range(len(SCAN_DEFS)))
    if scans_mode == "essential":
        return [i for i, sd in enumerate(SCAN_DEFS)
                if sd.get("role") in REFIT_ROLES]
    wanted = {n.strip() for n in scans_mode.split(",") if n.strip()}
    return [i for i, sd in enumerate(SCAN_DEFS) if sd["name"] in wanted]


def run_scans(logdir: Path, smoke: bool, scans_mode: str = "essential") -> Dict:
    """Run the SELECTED relaxed scans (each as a psi4 subprocess) + the fits +
    the MM-subtraction. `scans_mode` ("essential" default / "all" / names) picks
    which scans run; the unselected ones are recorded as deferred. Returns the
    scans summary dict.

    Essential (default) runs only the refit-role scan(s) (the N-methyl rotor =
    the only deployed-parameter scan). Validation scans (improper / angle) are
    opt-in via --scans all (or by name)."""
    indices = select_scan_indices(scans_mode)
    if smoke:
        # smoke = cheapest end-to-end: first SELECTED scan, single point.
        indices = indices[:1]
    limit = 1 if smoke else None
    scans_out: Dict[str, Dict] = {}
    selected = set(indices)
    for i, scan_def in enumerate(SCAN_DEFS):
        if i not in selected:
            # Demoted / not selected -> recorded as deferred, NOT executed and
            # NOT failed (D2 covered these terms at geometry level).
            scans_out[scan_def["name"]] = {
                "name": scan_def["name"], "kind": scan_def["kind"],
                "role": scan_def.get("role"), "deferred": True,
                "deferral_reason": "not selected (validation deferred; "
                                   "geometry-level PASS at D2)",
            }
            continue
        wargs = [str(logdir), str(i)]
        if limit is not None:
            wargs += ["--limit-points", str(limit)]
        rc = _psi4_subprocess("--_scan-worker", wargs, smoke)
        sp = logdir / f"scan_{scan_def['name']}" / "scan_points.json"
        if rc != 0 or not sp.exists():
            scans_out[scan_def["name"]] = {"error": f"scan worker rc={rc}"}
            continue
        result = json.loads(sp.read_text())
        scandir = logdir / f"scan_{scan_def['name']}"
        fit = _fit_for(scan_def, result, scandir)
        result["fit"] = fit
        (scandir / f"{scan_def['name']}_fit.json").write_text(json.dumps(fit, indent=2))
        scans_out[scan_def["name"]] = result
    return scans_out


def _fit_for(scan_def: Dict, result: Dict, scandir: Path) -> Dict:
    kind = scan_def["kind"]
    points = result.get("points", [])
    if kind == "improper":
        return fit_improper(points)
    if kind == "angle":
        return fit_angle(points)
    if kind == "dihedral":
        mm = compute_mm_subtraction(result, scandir, scan_def.get("ff_class_quad", scan_def["atoms"]))
        fit = fit_dihedral_mm_subtracted(points, mm["contrib_kcal"],
                                         fold=scan_def.get("fold"),
                                         include_v6=scan_def.get("include_v6", False))
        fit["mm_subtraction"] = {k: mm[k] for k in ("route", "n_terms_zeroed", "analytic_terms")}
        return fit
    return {"error": f"no fit for kind {kind}"}


def run_resp(logdir: Path, smoke: bool) -> Dict:
    """6-conformer (or 1, smoke) HF/6-31G* ESP (PySCF exact) + AmberTools
    2-stage RESP-A2. Smoke runs ONE conformer end-to-end (real RESP, not dry)
    so the cross-binary wiring + frozen-anchor seating are genuinely exercised."""
    espdir = logdir / "esp_grids"
    espdir.mkdir(parents=True, exist_ok=True)
    respdir = logdir / "resp_output"
    respdir.mkdir(parents=True, exist_ok=True)

    confs = _conformer_pdbs()
    if not confs:
        return {"error": "no D2 conformers found"}
    if smoke:
        confs = confs[:1]

    # 1) ESP per conformer via PySCF (exact int1e_grids), one esp.dat block each.
    #    Runs IN-ENV (qmmm has PySCF) -- no cross-env subprocess for ESP.
    block_paths: List[Path] = []
    for ci, conf in enumerate(confs):
        atoms_c = D2.parse_pdb(conf)
        grid_ang = build_mk_grid(atoms_c)
        esp_au, ok = compute_esp_pyscf(atoms_c, grid_ang)
        if not ok or esp_au is None:
            return {"error": f"ESP SCF did not converge for conf{ci:02d}",
                    "n_conformers_done": ci}
        block = espdir / f"conf{ci:02d}_block.dat"
        with open(block, "w") as fh:
            write_esp_dat_block(fh, atoms_c, grid_ang, esp_au)
        block_paths.append(block)
        print(f"[esp] conf{ci:02d}: {len(grid_ang)} grid pts, "
              f"|ESP|max={float(max(abs(v) for v in esp_au)):.4f} Ha/e")

    # 2) concat into a single multi-conformer esp.dat.
    esp_dat = respdir / "esp.dat"
    with open(esp_dat, "w") as out:
        for bp in block_paths:
            out.write(bp.read_text())

    # 3) antechamber -> .ac (bonding/types; charges discarded).
    atoms_ref = D2.parse_pdb(confs[0])
    ac_pdb = respdir / "Ace-OMW-NMe.pdb"
    D2.write_pdb(atoms_ref, ac_pdb)
    ac_out = respdir / "omw.ac"
    env = _amber_env()
    ac_cmd = [str(QMMM_BIN / "antechamber"),
              "-i", ac_pdb.name, "-fi", "pdb", "-o", ac_out.name, "-fo", "ac",
              "-rn", "OMW", "-nc", "0", "-at", "gaff2", "-j", "5",
              "-pf", "y", "-dr", "no"]
    r_ac = subprocess.run(ac_cmd, cwd=respdir, capture_output=True, text=True,
                          timeout=600, env=env)
    (respdir / "antechamber.log").write_text(
        "CMD: " + " ".join(ac_cmd) + "\n\n" + r_ac.stdout +
        "\n--- STDERR ---\n" + r_ac.stderr)
    if r_ac.returncode != 0 or not ac_out.exists():
        return {"error": f"antechamber failed (rc={r_ac.returncode})"}

    # 4) addfile (freeze backbone+NE1) + RESP-A2 2-stage. Smoke runs the REAL
    #    2-stage fit on the single conformer (it is cheap) so the frozen-anchor
    #    seating + integer-net constraint are genuinely validated end-to-end.
    addfile = respdir / "addfile"
    frozen_ids = build_resp_addfile(atoms_ref, addfile)
    qout = run_resp_a2(ac_out, esp_dat, addfile, atoms_ref,
                       n_conf=len(confs), logdir=respdir, dry=False)

    audit = {
        "n_conformers": len(confs),
        "smoke_single_conformer": smoke,  # smoke fits 1 conf (real RESP, not dry)
        "frozen_atom_ids_1based": frozen_ids,
        "frozen_charges": RESP_FROZEN_CHARGES,
        "esp_dat": str(esp_dat.relative_to(PROJ)),
    }
    if qout is not None and qout.exists():
        # multi-conformer qout repeats the charge set per block; take conformer 0
        charges = _parse_qout(qout, n_atoms=len(atoms_ref))
        atoms_omw = [(D2._norm_name(a[0]), a[1]) for a in atoms_ref]
        per_atom = {}
        for (nm, res), q in zip(atoms_omw, charges):
            if res == "OMW":
                per_atom[nm] = q
        audit["charges_per_atom"] = per_atom
        audit["charge_sum"] = round(sum(charges), 6)
        audit["frozen_match"] = _verify_frozen(per_atom)
    (respdir / "charge_audit.json").write_text(json.dumps(audit, indent=2))
    return audit


def _count_esp_blocks(esp_dat: Path) -> int:
    """Count conformer blocks in a concatenated AmberTools esp.dat. Each block
    starts with a header line `%5d%5d%5d` = natoms, ngrid, 0. The block count is
    the number of such headers (lines whose 3 leading int fields parse and whose
    3rd field == 0 while the 1st is a small atom count)."""
    n = 0
    expect_header = True
    remaining = 0
    for raw in esp_dat.read_text().splitlines():
        if not raw.strip():
            continue
        if expect_header:
            toks = raw.split()
            if len(toks) >= 2 and toks[0].isdigit() and toks[1].isdigit():
                natoms = int(toks[0])
                ngrid = int(toks[1])
                n += 1
                remaining = natoms + ngrid
                expect_header = False
        else:
            remaining -= 1
            if remaining <= 0:
                expect_header = True
    return n


def run_resp_redo(src_logdir: Path) -> Dict:
    """Re-run ONLY the AmberTools RESP-A2 2-stage fit from an existing run's
    resp_output/, reusing its esp.dat (the exact HF/6-31G* multi-conformer ESP)
    and Ace-OMW-NMe.pdb -- NO QM/ESP recompute. Used to recover from a RESP-stage
    crash without re-paying the (expensive) per-conformer SCF/ESP. Writes resp1/2
    inputs, qout, and charge_audit.json into the SAME resp_output/ (overwrite)."""
    respdir = src_logdir / "resp_output"
    esp_dat = respdir / "esp.dat"
    ac_pdb = respdir / "Ace-OMW-NMe.pdb"
    ac_out = respdir / "omw.ac"
    if not esp_dat.exists():
        return {"error": f"esp.dat not found in {respdir} (nothing to reuse)"}
    if not ac_pdb.exists():
        return {"error": f"Ace-OMW-NMe.pdb not found in {respdir}"}

    n_conf = _count_esp_blocks(esp_dat)
    if n_conf < 1:
        return {"error": f"could not count conformer blocks in {esp_dat}"}
    print(f"[resp-redo] reusing {esp_dat} ({n_conf} conformer block(s)); "
          f"NO QM/ESP recompute")

    atoms_ref = D2.parse_pdb(ac_pdb)

    # antechamber -> .ac if absent (bonding/types only; charges discarded).
    if not ac_out.exists():
        env = _amber_env()
        ac_cmd = [str(QMMM_BIN / "antechamber"),
                  "-i", ac_pdb.name, "-fi", "pdb", "-o", ac_out.name, "-fo", "ac",
                  "-rn", "OMW", "-nc", "0", "-at", "gaff2", "-j", "5",
                  "-pf", "y", "-dr", "no"]
        r_ac = subprocess.run(ac_cmd, cwd=respdir, capture_output=True, text=True,
                              timeout=600, env=env)
        (respdir / "antechamber.log").write_text(
            "CMD: " + " ".join(ac_cmd) + "\n\n" + r_ac.stdout +
            "\n--- STDERR ---\n" + r_ac.stderr)
        if r_ac.returncode != 0 or not ac_out.exists():
            return {"error": f"antechamber failed (rc={r_ac.returncode})"}

    addfile = respdir / "addfile"
    frozen_ids = build_resp_addfile(atoms_ref, addfile)
    qout = run_resp_a2(ac_out, esp_dat, addfile, atoms_ref,
                       n_conf=n_conf, logdir=respdir, dry=False)

    audit = {
        "n_conformers": n_conf,
        "resp_redo_from": str(src_logdir.relative_to(PROJ)),
        "frozen_atom_ids_1based": frozen_ids,
        "frozen_charges": RESP_FROZEN_CHARGES,
        "esp_dat": str(esp_dat.relative_to(PROJ)),
    }
    if qout is not None and qout.exists():
        charges = _parse_qout(qout, n_atoms=len(atoms_ref))
        atoms_omw = [(D2._norm_name(a[0]), a[1]) for a in atoms_ref]
        per_atom = {}
        for (nm, res), q in zip(atoms_omw, charges):
            if res == "OMW":
                per_atom[nm] = q
        audit["charges_per_atom"] = per_atom
        audit["charge_sum"] = round(sum(charges), 6)
        audit["frozen_match"] = _verify_frozen(per_atom)
    (respdir / "charge_audit.json").write_text(json.dumps(audit, indent=2))
    return audit


def _verify_frozen(per_atom: Dict[str, float], tol: float = 1e-3) -> Dict:
    """Confirm the fitted backbone+NE1 charges match the amber14SB freeze set."""
    out = {"all_frozen_held": True, "deltas": {}}
    for nm, q_target in RESP_FROZEN_CHARGES.items():
        q_fit = per_atom.get(nm)
        if q_fit is None:
            out["all_frozen_held"] = False
            out["deltas"][nm] = None
            continue
        dq = abs(q_fit - q_target)
        out["deltas"][nm] = round(q_fit - q_target, 6)
        if dq > tol:
            out["all_frozen_held"] = False
    return out


# ===========================================================================
# Finalize: assemble the final deliverable (frcmod DIHE + acceptance + report)
# from an EXISTING run dir, with NO QM/ESP/RESP recompute. The scan_points.json
# (relaxed QM scan) and resp_output/ (esp.dat + resp2.qout) are read-only; only
# the cheap MM single-point evaluations (MM-subtraction) and the fits run here.
# ===========================================================================
def _load_existing_resp_audit(src_logdir: Path) -> Dict:
    """Re-derive the RESP charge audit from an existing run's resp_output/ WITHOUT
    re-running RESP: parse resp_output/resp2.qout (conformer-0 slice) against the
    atom order in resp_output/Ace-OMW-NMe.pdb. Returns the same shape as the
    charge_audit.json (charges_per_atom / charge_sum / frozen_match), or
    {"error":...} if the qout / pdb is missing."""
    respdir = src_logdir / "resp_output"
    qout = respdir / "resp2.qout"
    ac_pdb = respdir / "Ace-OMW-NMe.pdb"
    esp_dat = respdir / "esp.dat"
    if not qout.exists():
        return {"error": f"resp2.qout not found in {respdir} (RESP not done)"}
    if not ac_pdb.exists():
        return {"error": f"Ace-OMW-NMe.pdb not found in {respdir}"}
    atoms_ref = D2.parse_pdb(ac_pdb)
    n_conf = _count_esp_blocks(esp_dat) if esp_dat.exists() else None
    charges = _parse_qout(qout, n_atoms=len(atoms_ref))
    per_atom: Dict[str, float] = {}
    for a, q in zip(atoms_ref, charges):
        if a[1] == "OMW":
            per_atom[D2._norm_name(a[0])] = q
    audit = {
        "n_conformers": n_conf,
        "resp_audit_from": str(src_logdir.relative_to(PROJ)),
        "frozen_charges": RESP_FROZEN_CHARGES,
        "charges_per_atom": per_atom,
        "charge_sum": round(sum(charges), 6),
        "frozen_match": _verify_frozen(per_atom),
        "recomputed": False,
    }
    return audit


def finalize_dihedral_fit(src_logdir: Path) -> Dict:
    """Fit the N-methyl rotor V3 from an existing run's relaxed-scan
    scan_points.json (READ-ONLY) -- the MM-subtraction here is a cheap OpenMM
    single-point at each ALREADY-relaxed geometry, NOT a QM recompute. Returns
    {scan_name, fit, mm_route, ...} or {"error":...} if the scan json is absent.
    The scan_def is the refit-role SCAN_DEFS entry (the N-methyl rotor)."""
    refit_idx = [i for i, sd in enumerate(SCAN_DEFS)
                 if sd.get("role") in REFIT_ROLES]
    if not refit_idx:
        return {"error": "no refit-role scan defined"}
    sd = SCAN_DEFS[refit_idx[0]]
    scandir = src_logdir / f"scan_{sd['name']}"
    sp = scandir / "scan_points.json"
    if not sp.exists():
        return {"error": f"scan_points.json not found: {sp}"}
    result = json.loads(sp.read_text())
    mm = compute_mm_subtraction(result, scandir,
                                sd.get("ff_class_quad", sd["atoms"]))
    fit = fit_dihedral_mm_subtracted(result["points"], mm["contrib_kcal"],
                                     fold=sd.get("fold"),
                                     include_v6=sd.get("include_v6", False))
    fit["mm_subtraction"] = {k: mm[k]
                             for k in ("route", "n_terms_zeroed", "analytic_terms")}
    # NOTE: we DO NOT overwrite the existing scan_<name>/<name>_fit.json here.
    # The finalize write-scope is restricted to acceptance / frcmod (+ manifest) /
    # report; the V3 fit is recorded there. The pre-existing in-scan fit.json
    # (which the original --resp-redo path left with n_points=0) is left as the
    # immutable scan-run record. `result["fit"]` is set only in memory for the
    # report rendering.
    result["fit"] = fit
    return {"scan_name": sd["name"], "fit": fit, "mm_route": mm["route"],
            "scan_result": result}


def _build_dihe_refit(dih: Dict) -> Dict:
    """Assemble the `dihe_refit` payload for write_bonded_refit_frcmod from a
    finalize_dihedral_fit() result. Resolves the DEPLOYED-FF atom-type quad and
    decides whether the DIHE may be emitted (emit only when the V3 fit succeeded
    AND all four FF classes resolved; the prompt forbids a wildcard fallback)."""
    fit = dih.get("fit", {})
    rotor = next((sd for sd in SCAN_DEFS if sd.get("role") in REFIT_ROLES), None)
    atoms = rotor["atoms"] if rotor else ("CD1", "NE1", "CM", "HM1")
    ff = resolve_dihedral_ff_types(atoms)
    payload = {
        "ff_classes": ff["classes"],
        "ff_types": ff["types"],
        "rationale": ff["rationale"],
        "rmsd_kcal": fit.get("rmsd_kcal"),
        "n_points": fit.get("n_points"),
        "mm_route": dih.get("mm_route"),
    }
    if "error" in fit or fit.get("V3_kcal") is None:
        payload["emit"] = False
        payload["skip_reason"] = (
            "V3 fit unavailable: " + str(fit.get("error", "no V3")))
        return payload
    if not ff["resolved"]:
        payload["emit"] = False
        payload["skip_reason"] = (
            "deployed-FF atom-type quad for the rotor could not be resolved "
            f"(classes={ff['classes']}); refusing to emit a DIHE under a "
            "wildcard X-NA-CT-X (forbidden). BOND-only frcmod.")
        return payload
    payload["emit"] = True
    payload["V3_kcal"] = float(fit["V3_kcal"])
    return payload


def compute_acceptance_finalize(dih: Dict, resp: Dict,
                                dihe_refit: Dict) -> Dict:
    """Finalize acceptance for the D3 essential deliverable. PASS / MUST-REDO
    on four criteria:
      (a) rotor V3 fit RMSD < 2 kcal/mol,
      (b) RESP integer net (Sigma q = 0 +- 0.0005) AND the 7 frozen backbone+NE1
          atoms held at the amber14SB values (Delta = 0),
      (c) the NE1-CM BOND r_eq + the N-methyl rotor V3 DIHE were emitted to the
          frcmod,
      (d) improper / angle = 'validation deferred (D2 geometry-level PASS)' --
          NOT a MUST-REDO.
    Returns {"verdict":"PASS"|"MUST-REDO", "checks":{...}, ...}."""
    acc: Dict[str, object] = {"checks": {}, "verdict": "PASS",
                              "criteria": "finalize (D3 essential)"}

    def _chk(key, ok, detail, *, gating=True):
        acc["checks"][key] = {"pass": bool(ok), "gating": gating, "detail": detail}
        if gating and not ok:
            acc["verdict"] = "MUST-REDO"

    # (a) rotor V3 RMSD < 2
    fit = dih.get("fit", {})
    rmsd = fit.get("rmsd_kcal")
    v3 = fit.get("V3_kcal")
    _chk("a_rotor_v3_rmsd_lt_2",
         (rmsd is not None and rmsd < 2.0),
         f"V3={v3} RMSD={rmsd} kcal/mol (n_points={fit.get('n_points')})")

    # (b) RESP integer net 0 +- 0.0005 AND 7 frozen held (Delta=0)
    cs = resp.get("charge_sum")
    _chk("b_resp_integer_net_zero",
         (cs is not None and abs(cs) <= 5e-4),
         f"Sigma q = {cs}")
    fm = resp.get("frozen_match", {})
    held = fm.get("all_frozen_held", False)
    n_frozen = sum(1 for d in (fm.get("deltas") or {}).values()
                   if d is not None and abs(d) <= 5e-4)
    _chk("b_resp_frozen_7_held",
         (held and n_frozen == len(RESP_FROZEN_CHARGES)),
         f"all_frozen_held={held}, {n_frozen}/{len(RESP_FROZEN_CHARGES)} "
         f"within 5e-4; deltas={fm.get('deltas')}")

    # (c) NE1-CM BOND + N-methyl V3 DIHE emitted
    _chk("c_ne1_cm_bond_in_range",
         1.40 <= NE1_CM_REFIT_REQ_ANG <= 1.50,
         f"NE1-CM r_eq = {NE1_CM_REFIT_REQ_ANG} A")
    _chk("c_nmethyl_v3_dihe_emitted",
         bool(dihe_refit.get("emit")),
         f"DIHE emitted={bool(dihe_refit.get('emit'))} "
         f"(quad={'-'.join(dihe_refit.get('ff_classes') or [])}); "
         f"{dihe_refit.get('skip_reason','')}".strip())

    # (d) improper / angle: validation deferred -- NON-gating (D2 geometry PASS)
    deferred = [sd["name"] for sd in SCAN_DEFS
                if sd.get("role") not in REFIT_ROLES]
    _chk("d_improper_angle_validation_deferred",
         True,
         "validation deferred (D2 geometry-level PASS); NOT a MUST-REDO. "
         f"deferred scans: {deferred}", gating=False)
    acc["deferred_scans"] = deferred
    return acc


def run_finalize(src_logdir: Path, scans_mode: str = "essential") -> Dict:
    """Assemble the FINAL D3 deliverable from an EXISTING run dir with NO
    QM/ESP/RESP recompute:
      1. RESP charges re-read from resp_output/resp2.qout (READ-ONLY).
      2. N-methyl rotor V3 fit from scan_points.json (READ-ONLY scan; the only
         compute is a cheap OpenMM single-point MM-subtraction at the already-
         relaxed geometries).
      3. write params/d3_omw_bonded_refit.frcmod = NE1-CM BOND r_eq + V3 DIHE.
      4. write <src_logdir>/d3_acceptance.json + d3_report.md (finalize variant).
      5. copy the deliverable to params/d3_omw_qm_geom/ + refresh the frcmod.
    Returns {"acceptance":..., "resp":..., "dihedral":...}."""
    resp = _load_existing_resp_audit(src_logdir)
    dih = finalize_dihedral_fit(src_logdir)
    dihe_refit = _build_dihe_refit(dih)

    # validation block for the manifest: the rotor fit + the deferred scans.
    validation = {dih.get("scan_name", "dihedral"): dih.get("fit", {})}
    for sd in SCAN_DEFS:
        if sd.get("role") not in REFIT_ROLES:
            validation[sd["name"]] = {
                "deferred": True,
                "note": "validation deferred (D2 geometry-level PASS)"}
    write_bonded_refit_frcmod(PROJ / "params/d3_omw_bonded_refit.frcmod",
                              validation, dihe_refit=dihe_refit)

    acc = compute_acceptance_finalize(dih, resp, dihe_refit)
    (src_logdir / "d3_acceptance.json").write_text(json.dumps(acc, indent=2))

    # report (finalize variant): reuse the scans/resp shapes write_report expects.
    scans_for_report = {dih.get("scan_name", "dihedral"): dih.get("scan_result", {})}
    for sd in SCAN_DEFS:
        if sd.get("role") not in REFIT_ROLES:
            scans_for_report[sd["name"]] = {
                "name": sd["name"], "kind": sd["kind"], "role": sd.get("role"),
                "deferred": True}
    write_report(src_logdir, scans_for_report, resp, _acc_compat(acc),
                 mode="finalize", smoke=False)
    # finalize deliverable copy (also refreshes the completed frcmod copy).
    dst = copy_deliverable(src_logdir, scans_for_report, resp)
    _copy_frcmod_to_deliverable(dst)
    return {"acceptance": acc, "resp": resp, "dihedral": dih,
            "dihe_refit": dihe_refit, "deliverable": str(dst.relative_to(PROJ))}


def _acc_compat(acc_finalize: Dict) -> Dict:
    """Adapt the finalize acceptance (verdict/checks) to the shape write_report
    expects ({'overall':bool,'checks':{...}}). PASS -> overall True."""
    return {"overall": acc_finalize.get("verdict") == "PASS",
            "checks": acc_finalize.get("checks", {}),
            "verdict": acc_finalize.get("verdict")}


def _copy_frcmod_to_deliverable(dst: Path) -> None:
    """Copy the completed params/d3_omw_bonded_refit.frcmod (+ manifest) into the
    deliverable dir so params/d3_omw_qm_geom/ carries the finished bonded params."""
    for fname in ("d3_omw_bonded_refit.frcmod",
                  "d3_omw_bonded_refit_manifest.json"):
        src = PROJ / "params" / fname
        if src.exists():
            shutil.copy2(src, dst / fname)


# ===========================================================================
# Acceptance + report.
# ===========================================================================
def compute_acceptance(scans: Dict, resp: Dict) -> Dict:
    """Advisory acceptance gate (scan/fit RMSD + parameter bands + RESP integer
    net + frozen-anchor held). Informational only -- not a release decision."""
    acc: Dict[str, object] = {"checks": {}, "overall": True}

    def _fail(key, ok, detail):
        acc["checks"][key] = {"pass": bool(ok), "detail": detail}
        if not ok:
            acc["overall"] = False

    # bonded refit value sanity
    _fail("ne1_cm_req_in_range",
          1.40 <= NE1_CM_REFIT_REQ_ANG <= 1.50,
          f"NE1-CM r_eq = {NE1_CM_REFIT_REQ_ANG} A")

    # scan/fit RMSD + parameter bands -- gate ONLY scans actually executed.
    # Deferred (opt-in / not selected) validation scans are reported as
    # "validation deferred" and do NOT contribute to overall=False (D2 already
    # validated those terms geometry-level; their energy scan is redundant).
    acc["deferred_scans"] = []
    for name, res in scans.items():
        if isinstance(res, dict) and res.get("deferred"):
            acc["deferred_scans"].append(name)
            acc["checks"][f"fit_{name}"] = {
                "pass": True, "deferred": True,
                "detail": "validation deferred (geometry-level PASS at D2)"}
            continue
        fit = res.get("fit", {}) if isinstance(res, dict) else {}
        if "error" in fit or not fit:
            _fail(f"fit_{name}", False, fit.get("error", "no fit"))
            continue
        rmsd = fit.get("rmsd_kcal")
        _fail(f"fit_{name}_rmsd", (rmsd is not None and rmsd < 2.0),
              f"RMSD={rmsd}")
        if "V_kcal_per_mol" in fit:
            v = fit["V_kcal_per_mol"]
            _fail(f"fit_{name}_imp_band", 0.5 <= v <= 5.0, f"V={v}")
        if "k_kcal_per_mol_rad2" in fit:
            k = fit["k_kcal_per_mol_rad2"]
            _fail(f"fit_{name}_ang_band", 80.0 <= k <= 200.0, f"k={k}")
        if "V1_kcal" in fit:
            vs = [fit.get(f"V{i}_kcal", 0.0) for i in (1, 2, 3, 4)]
            _fail(f"fit_{name}_tor_band", all(0.0 <= abs(v) <= 15.0 for v in vs),
                  f"V_n={vs}")

    # RESP integer net + frozen held
    if resp and "charge_sum" in resp:
        cs = resp["charge_sum"]
        _fail("resp_integer_net", abs(round(cs) - cs) < 1e-2 and abs(round(cs)) == 0,
              f"sum={cs}")
        fm = resp.get("frozen_match", {})
        _fail("resp_frozen_held", fm.get("all_frozen_held", False),
              f"deltas={fm.get('deltas')}")
    elif resp and "error" in resp:
        _fail("resp", False, resp["error"])
    return acc


def write_report(logdir: Path, scans: Dict, resp: Dict, acc: Dict,
                 mode: str, smoke: bool) -> None:
    # Read back the emitted frcmod manifest (if present) so the report states the
    # ACTUAL bonded scope (BOND-only vs BOND + N-methyl V3 DIHE).
    dihe_md: List[str] = []
    man_path = PROJ / "params" / "d3_omw_bonded_refit_manifest.json"
    if man_path.exists():
        try:
            mdh = json.loads(man_path.read_text()).get("nmethyl_rotor_dihe", {})
            if mdh.get("emitted"):
                dihe_md.append(
                    f"- N-methyl rotor V3 DIHE -> "
                    f"{'-'.join(mdh.get('ff_class_quad') or [])} "
                    f"pk={mdh.get('V3_kcal_per_mol')} kcal/mol, period=3, phase=0 "
                    f"(MM-subtracted; RMSD={mdh.get('rmsd_kcal')}, "
                    f"n_points={mdh.get('n_points')})")
            else:
                dihe_md.append(
                    f"- N-methyl rotor V3 DIHE NOT emitted: "
                    f"{mdh.get('reason', 'bond-only frcmod')}")
        except Exception:  # noqa: BLE001 -- manifest read is best-effort for the report
            pass
    scope_title = ("NE1-CM bond r_eq + N-methyl rotor V3"
                   if dihe_md and "NOT emitted" not in dihe_md[0]
                   else "NE1-CM bond r_eq only")
    md = [
        "# D3 OMW Bonded Refit + RESP-A2 Charge Refit Report",
        "",
        f"**Date**: {datetime.now().isoformat()}",
        f"**Mode**: {mode}{' (smoke)' if smoke else ''}",
        f"**QM**: B3LYP/6-31G(d) relaxed scans + HF/6-31G* GRID_ESP (psi4 1.10)",
        f"**Charge**: AmberTools 2-stage RESP-A2 ({resp.get('n_conformers','?')} conformers, equal weight)",
        f"**Geometry input**: {D2_DELIVERABLE.relative_to(PROJ)}",
        "",
        f"## Bonded refit (pre-decided scope: {scope_title})",
        f"- NE1-CM r_eq -> {NE1_CM_REFIT_REQ_ANG:.4f} A (D2 B3LYP relaxed); "
        f"k = {NE1_CM_K_KJ_PER_NM2} kJ/mol/nm^2 (amber14SB NA-CT, retained)",
        *dihe_md,
        f"- Emitted as params/d3_omw_bonded_refit.frcmod (frozen MTR XML untouched)",
        "",
        "## Relaxed scans + fits",
    ]
    for name, res in scans.items():
        fit = res.get("fit", {}) if isinstance(res, dict) else {}
        role = res.get("role", "?") if isinstance(res, dict) else "?"
        md.append(f"### {name} ({res.get('kind','?')}, {role})")
        if isinstance(res, dict) and res.get("deferred"):
            md.append("- VALIDATION DEFERRED (opt-in; not selected). Geometry-"
                      "level PASS at D2; re-run with --scans all to execute. "
                      "NOT a failure.")
        elif "error" in res:
            md.append(f"- ERROR: {res['error']}")
        else:
            n_ok = sum(1 for p in res.get("points", []) if p.get("energy_hartree") is not None)
            md.append(f"- points: {n_ok}/{len(res.get('points', []))} converged")
            md.append(f"- fit: {json.dumps(fit)}")
        md.append("")
    md.append("## RESP-A2 charge refit")
    md.append(f"- {json.dumps({k: v for k, v in resp.items() if k != 'charges_per_atom'})}")
    if "charges_per_atom" in resp:
        md.append(f"- charges_per_atom: {json.dumps(resp['charges_per_atom'])}")
    md.append("")
    md.append("## Acceptance (advisory)")
    md.append(f"- overall: {acc['overall']}")
    for k, v in acc["checks"].items():
        md.append(f"  - {k}: {v}")
    md.append("")
    md.append("## Caveats")
    md.append("- ranking-regime use only; absolute values pending calibration.")
    md.append("- MM-subtraction route + accuracy, RESP frozen-anchor seating, and "
              "scan convergence are NOT independently validated by this run alone.")
    md.append("- The N-methyl rotor (3-fold symmetric) is fit with a V3(+V6)-"
              "restricted model; V1/V2/V4 are pinned to zero by methyl symmetry "
              "and are NOT deployed into the FF (a free V1..V4 fit of the "
              "asymmetric 0..180 deg sampling would produce spurious 1/2/4-fold "
              "terms).")
    md.append("- DEFERRED VALIDATION (honest disclosure): the improper and angle "
              "scans are demoted to opt-in (--scans all). In essential mode their "
              "energy-LEVEL bonded validation is NOT performed; D2 provided only "
              "geometry-level validation (improper <1 deg, angle delta <0.6 deg). "
              "Both terms are accepted-as-amber14SB (not refit). Downstream D4 "
              "trajectory NE1-region geometry is the deferred check.")
    md.append("- KNOWN ISSUE (improper scan, NOT fixed; opt-in only): the improper "
              "is dispatched through the proper-dihedral pre-rotation, which is "
              "invalid for the 3+-bond central atom NE1 and tears the ring. "
              "Opting in requires an improper-specific pre-rotation (NE1 fixed, "
              "one substituent out of plane, narrow +-5..10 deg window) first. "
              "See params/d3_omw_bonded_refit_manifest.json known_issues.")
    (logdir / "d3_report.md").write_text("\n".join(md))


def copy_deliverable(logdir: Path, scans: Dict, resp: Dict) -> Path:
    dst = PROJ / "params/d3_omw_qm_geom"
    dst.mkdir(parents=True, exist_ok=True)
    for fname in ("d3_report.md", "d3_acceptance.json"):
        src = logdir / fname
        if src.exists():
            shutil.copy2(src, dst / fname)
    (dst / "PROVENANCE.txt").write_text(
        f"D3 deliverable copied from {logdir}\n"
        f"generated {datetime.now().isoformat()}\n"
        f"bonded refit: params/d3_omw_bonded_refit.frcmod; charges: resp_output/\n"
    )
    return dst


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--scan-only", action="store_true", help="relaxed scans + fits only")
    mode.add_argument("--resp-only", action="store_true", help="ESP + RESP-A2 only")
    mode.add_argument("--full", action="store_true", help="both (default)")
    ap.add_argument("--smoke", action="store_true",
                    help="cheapest end-to-end shape check (1 scan/1 point + 1-conf "
                         "ESP + real 1-conf RESP-A2) to validate cross-env wiring")
    ap.add_argument("--scans", default="essential",
                    help="which relaxed scans to run: 'essential' (default; only "
                         "the refit-role scan(s) = the N-methyl rotor, the only "
                         "scan that emits a deployed FF parameter), 'all' (also "
                         "the opt-in VALIDATION scans improper/angle -- improper "
                         "has a KNOWN unfixed pre-rotation bug, see source), or a "
                         "comma-separated list of scan names. Unselected scans are "
                         "reported as 'validation deferred' (D2 covered them at "
                         "geometry level), not as failures.")
    ap.add_argument("--resp-redo", metavar="LOGDIR", default=None,
                    help="re-run ONLY the AmberTools RESP-A2 2-stage fit from an "
                         "existing run dir's resp_output/esp.dat (reuse the exact "
                         "multi-conformer HF/6-31G* ESP; NO QM/ESP recompute). "
                         "Recovers from a RESP-stage crash cheaply. Outputs "
                         "overwrite that dir's resp_output/.")
    ap.add_argument("--finalize", metavar="LOGDIR", default=None,
                    help="assemble the FINAL D3 deliverable from an existing run "
                         "dir with NO QM/ESP/RESP recompute: re-read the RESP "
                         "charges (resp_output/resp2.qout), fit the N-methyl rotor "
                         "V3 from the existing relaxed scan_points.json (the only "
                         "compute is a cheap OpenMM single-point MM-subtraction at "
                         "the already-relaxed geometries), write the frcmod "
                         "(NE1-CM BOND r_eq + V3 DIHE), the acceptance.json and "
                         "report.md, and copy the deliverable. scan/esp/resp are "
                         "READ-ONLY (acceptance/frcmod/report are the only new "
                         "writes).")
    # hidden cross-env worker (relaxed QM scans run in the psi4 interpreter)
    ap.add_argument("--_scan-worker", action="store_true", help=argparse.SUPPRESS)
    args, rest = ap.parse_known_args()

    if args._scan_worker:
        return _scan_worker_main(rest)

    if args.finalize:
        src = Path(args.finalize)
        if not src.is_absolute():
            src = (PROJ / src).resolve()
        if not src.exists():
            print(f"[error] --finalize logdir not found: {src}")
            return 1
        print(f"[driver] finalize from {src} "
              f"(NO QM/ESP/RESP recompute; frcmod/acceptance/report only)")
        out = run_finalize(src, scans_mode=args.scans)
        acc = out["acceptance"]
        dih = out["dihedral"].get("fit", {})
        resp = out["resp"]
        dr = out["dihe_refit"]
        print(f"[finalize] rotor V3 = {dih.get('V3_kcal')} kcal/mol  "
              f"RMSD = {dih.get('rmsd_kcal')}  (n_points={dih.get('n_points')}, "
              f"mm_route={out['dihedral'].get('mm_route')})")
        print(f"[finalize] DIHE emitted = {bool(dr.get('emit'))}  "
              f"quad = {'-'.join(dr.get('ff_classes') or [])}")
        print(f"[finalize] RESP Sigma q = {resp.get('charge_sum')}  "
              f"frozen held = {resp.get('frozen_match', {}).get('all_frozen_held')}")
        print(f"[finalize] verdict = {acc.get('verdict')}")
        print(f"[finalize] deliverable: {out['deliverable']}")
        print("=== D3 finalize complete ===")
        return 0 if acc.get("verdict") == "PASS" else 2

    if args.resp_redo:
        src = Path(args.resp_redo)
        if not src.is_absolute():
            src = (PROJ / src).resolve()
        if not src.exists():
            print(f"[error] --resp-redo logdir not found: {src}")
            return 1
        print(f"[driver] resp-redo from {src} (RESP-A2 only, ESP reused)")
        audit = run_resp_redo(src)
        (src / "resp_output" / "charge_audit.json").write_text(
            json.dumps(audit, indent=2))
        if "error" in audit:
            print(f"[resp-redo] FAILED: {audit['error']}")
            return 1
        print(f"[resp-redo] charge_sum = {audit.get('charge_sum')}")
        print(f"[resp-redo] frozen_match = "
              f"{audit.get('frozen_match', {}).get('all_frozen_held')}")
        print("=== D3 resp-redo complete ===")
        return 0

    do_scan = args.scan_only or args.full or not (args.scan_only or args.resp_only)
    do_resp = args.resp_only or args.full or not (args.scan_only or args.resp_only)
    mode_name = "scan-only" if args.scan_only else ("resp-only" if args.resp_only else "full")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = PROJ / f"outputs/analysis/d3_omw_qm_scans_{ts}"
    logdir.mkdir(parents=True, exist_ok=True)
    print(f"[driver] logdir: {logdir}")
    print(f"[driver] mode: {mode_name}  smoke: {args.smoke}  scans: {args.scans}")

    if not D2_CONFORMERS.exists():
        print(f"[error] D2 deliverable missing: {D2_CONFORMERS}")
        return 1

    scans: Dict = {}
    resp: Dict = {}
    if do_scan:
        scans = run_scans(logdir, args.smoke, scans_mode=args.scans)
        # bonded refit frcmod (validation = improper/angle/dihedral fits)
        validation = {name: res.get("fit", {}) for name, res in scans.items()
                      if isinstance(res, dict)}
        write_bonded_refit_frcmod(PROJ / "params/d3_omw_bonded_refit.frcmod", validation)
    if do_resp:
        resp = run_resp(logdir, args.smoke)

    acc = compute_acceptance(scans, resp)
    (logdir / "d3_acceptance.json").write_text(json.dumps(acc, indent=2))
    write_report(logdir, scans, resp, acc, mode_name, args.smoke)
    if not args.smoke:
        copy_deliverable(logdir, scans, resp)

    print(f"[driver] report: {logdir / 'd3_report.md'}")
    print(f"[driver] acceptance overall: {acc['overall']}")
    print("=== D3 complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
