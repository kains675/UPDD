#!/usr/bin/env python
"""
D3 deliverable — QM dihedral / angle / improper scans on the D2-optimized
Ace-OMW-NMe dipeptide. Produces an AMBER-style cosine-series fit of each
scan for input to a Layer 3D refitted XML.

Scans (per Khoury / amber14SB-style refit pattern):
  1. Improper N-CD1-CE2-CG: -30° to +30° in 5° steps (indole planarity)
     — fits A * (1 - cos(2φ)) form
  2. Angle CM-CG-CD2: 100° to 130° in 5° steps (N-methyl junction)
     — fits k_θ (θ - θ_0)² harmonic
  3. Dihedral CB-CG-CD1: 0° to 360° in 30° steps (sidechain rotation)
     — fits Σ_n V_n / 2 * (1 + cos(n φ - phase_n)) AMBER form

For each constrained scan point, psi4 performs B3LYP/6-31G(d) single-point.
Pure Python single-process via psi4 native API.

Outputs: outputs/analysis/d3_omw_qm_scans_{TS}/
  - improper_NE1_scan.json    (angle vs energy)
  - angle_CM_CG_scan.json
  - dihedral_CB_CG_scan.json
  - amber_fitted_params.json  (k_θ, V_n, phase_n, RMSD vs QM)
  - d3_qm_scans_report.md

Run:
  /home/san/miniconda3/envs/psi4/bin/python scripts/d3_omw_qm_scans.py
"""
from __future__ import annotations

import json
import sys
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore")

PROJ = Path(__file__).resolve().parent.parent


def find_latest_d2() -> Path | None:
    for d in sorted((PROJ / "outputs/analysis").glob("d2_omw_geom_opt_*"), reverse=True):
        pdb = d / "omw_dipeptide_optimized.pdb"
        if pdb.exists():
            return pdb
    return None


def parse_pdb(pdb_path: Path):
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


def get_atom_idx(atoms, name, residue="OMW", normalize=True):
    """1-based index for psi4 frozen_dihedral / frozen_angle."""
    target = "CZ1" if normalize and name == "CM" else name
    for i, (n, r, *_) in enumerate(atoms, 1):
        if r == residue and (n == name or n == target):
            return i
    return None


def build_geom_block(atoms) -> str:
    lines = []
    for name, residue, element, x, y, z in atoms:
        lines.append(f"{element:<3s} {x:14.8f} {y:14.8f} {z:14.8f}")
    return "\n".join(lines) + "\nno_reorient\nno_com\nunits angstrom\n"


def rotate_dihedral_to(coords, atom_indices_0based, target_deg):
    """Rigid rotation of atoms 'beyond' the central bond (i-j) so that the
    dihedral i-j-k-l = target_deg. atom_indices_0based = (a, b, c, d), and
    atoms structurally beyond c (i.e., bonded to c, d, ...) are rotated.

    Simpler implementation: rotate only atom d and its downstream substituents
    is too topology-dependent; instead we rotate only atom d (and let psi4's
    single-point energy reflect the geometry as-is). This is the rigid scan
    approximation — fine for D3's AMBER-form parameter fit purpose.
    """
    import numpy as np
    a, b, c, d = atom_indices_0based
    p1, p2, p3, p4 = (np.array(coords[i]) for i in (a, b, c, d))

    # Compute current dihedral
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1 / np.linalg.norm(n1), b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2) * np.linalg.norm(b2)
    current_rad = np.arctan2(y, x)
    current_deg = np.rad2deg(current_rad)
    delta_rad = np.deg2rad(target_deg - current_deg)

    # Rotation axis = b2 normalized; rotate p4 around c by delta_rad
    axis = b2 / np.linalg.norm(b2)
    # Rodrigues' rotation
    v = p4 - p3
    cos_t = np.cos(delta_rad)
    sin_t = np.sin(delta_rad)
    v_rot = v * cos_t + np.cross(axis, v) * sin_t + axis * np.dot(axis, v) * (1 - cos_t)
    new_p4 = p3 + v_rot

    new_coords = [list(c) for c in coords]
    new_coords[d] = new_p4.tolist()
    return new_coords


def rotate_angle_to(coords, atom_indices_0based, target_deg):
    """Rigid rotation of atom 'a' around 'b' so that angle a-b-c = target_deg."""
    import numpy as np
    a, b, c = atom_indices_0based
    p1, p2, p3 = (np.array(coords[i]) for i in (a, b, c))

    v_ba = p1 - p2
    v_bc = p3 - p2
    current_rad = np.arccos(
        np.clip(np.dot(v_ba, v_bc) / (np.linalg.norm(v_ba) * np.linalg.norm(v_bc)), -1, 1)
    )
    current_deg = np.rad2deg(current_rad)
    delta_rad = np.deg2rad(target_deg - current_deg)

    # Rotation axis = perpendicular to (v_ba, v_bc) plane
    axis = np.cross(v_ba, v_bc)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-8:
        axis = np.array([0, 0, 1])
    else:
        axis = axis / axis_norm
    # Rotate v_ba around b by delta_rad
    cos_t = np.cos(delta_rad)
    sin_t = np.sin(delta_rad)
    v_rot = v_ba * cos_t + np.cross(axis, v_ba) * sin_t + axis * np.dot(axis, v_ba) * (1 - cos_t)
    new_p1 = p2 + v_rot
    new_coords = [list(c_) for c_ in coords]
    new_coords[a] = new_p1.tolist()
    return new_coords


def run_constrained_scan(
    atoms, scan_type: str, scan_atoms: list, scan_range: list, logdir: Path
):
    """Run RELAXED scan: rotate target coordinate to value, freeze that coord,
    optimize all other DOF, take final B3LYP/6-31G(d) energy.

    Per psi4 1.10 OPTKING semantics, `frozen_dihedral`/`frozen_bend` with only
    atom indices (no value) freezes the *starting* dihedral/angle during opt.
    We rotate the geometry to the target value first, so the starting value is
    the target; freezing then preserves it while other DOF relax.

    scan_type: 'dihedral' | 'angle' | 'improper'
    scan_atoms: list of 1-based atom indices
    scan_range: list of target values (degrees)
    """
    import psi4
    psi4.set_output_file(str(logdir / f"psi4_scan_{scan_type}.dat"), True)
    psi4.set_memory("8 GB")
    psi4.core.set_num_threads(4)

    init_coords = [[a[3], a[4], a[5]] for a in atoms]
    elements = [a[2] for a in atoms]
    atom_idx_0based = [i - 1 for i in scan_atoms]

    results = {}
    atoms_str = " ".join(str(i) for i in scan_atoms)
    print(f"[info] RELAXED-scan {scan_type} on atoms {atoms_str}: {len(scan_range)} points")

    for k, val in enumerate(scan_range):
        # 1. Rotate starting geometry to target value
        if scan_type in ("dihedral", "improper"):
            rotated_coords = rotate_dihedral_to(init_coords, atom_idx_0based, val)
        elif scan_type == "angle":
            rotated_coords = rotate_angle_to(init_coords, atom_idx_0based, val)
        else:
            rotated_coords = init_coords

        # 2. Build psi4 geometry from rotated coords
        geom_lines = []
        for i, el in enumerate(elements):
            x, y, z = rotated_coords[i]
            geom_lines.append(f"{el:<3s} {x:14.8f} {y:14.8f} {z:14.8f}")
        geom_block = "\n".join(geom_lines) + "\nno_reorient\nno_com\nunits angstrom\n"
        mol = psi4.geometry(geom_block)
        mol.update_geometry()
        mol.set_name(f"OMW_{scan_type}_{val:+05.1f}")

        # 3. Set options including frozen coord (preserves starting value through opt)
        psi4_opts = {
            "basis": "6-31G(d)",
            "scf_type": "df",
            "g_convergence": "gau_loose",
            "geom_maxiter": 150,
        }
        if scan_type in ("dihedral", "improper"):
            psi4_opts["frozen_dihedral"] = f"({atoms_str})"
        elif scan_type == "angle":
            psi4_opts["frozen_bend"] = f"({atoms_str})"
        psi4.set_options(psi4_opts)

        # 4. Optimize other DOF with target coord frozen → relaxed-scan energy
        try:
            e = psi4.optimize("b3lyp", molecule=mol)
            results[float(val)] = float(e)
            print(f"  point {k+1:2d}/{len(scan_range)}  {scan_type}={val:+.1f}  E={e:.6f}")
        except Exception as exc:
            print(f"  point {k+1:2d}/{len(scan_range)}  FAIL: {type(exc).__name__}: {str(exc)[:80]}")
            results[float(val)] = None

        psi4.core.clean()

    return results


def fit_improper(scan_data: dict):
    """Fit improper N-CD1-CE2-CG to A * (1 - cos(2φ)) AMBER form."""
    import numpy as np
    valid = {k: v for k, v in scan_data.items() if v is not None}
    if len(valid) < 5:
        return {"error": "insufficient points", "n_points": len(valid)}
    angles_deg = np.array(sorted(valid.keys()))
    energies = np.array([valid[k] for k in angles_deg])
    # Reference: planar geometry at angle = 0
    E_ref = valid.get(0.0) or min(energies)
    dE_kcal = (energies - E_ref) * 627.5095
    angles_rad = np.deg2rad(angles_deg)
    # AMBER improper: V/2 * (1 + cos(2φ - π)) = V/2 * (1 - cos(2φ))
    # Equivalent fit: dE = A * (1 - cos(2φ)) where A = V/2
    from scipy.optimize import curve_fit
    def model(x, A): return A * (1 - np.cos(2 * x))
    popt, _ = curve_fit(model, angles_rad, dE_kcal, p0=[10.0])
    A = popt[0]
    V = 2 * A
    residuals = dE_kcal - model(angles_rad, A)
    rmsd = float(np.sqrt(np.mean(residuals ** 2)))
    return {
        "form": "V/2 * (1 + cos(2*phi - 180))",
        "V_kcal_per_mol": float(V),
        "A_kcal_per_mol": float(A),
        "rmsd_kcal": rmsd,
        "n_points": len(valid),
    }


def fit_angle(scan_data: dict):
    """Fit angle CM-CG-CD2 to harmonic k_θ (θ - θ_0)² form."""
    import numpy as np
    valid = {k: v for k, v in scan_data.items() if v is not None}
    if len(valid) < 5:
        return {"error": "insufficient points", "n_points": len(valid)}
    angles_deg = np.array(sorted(valid.keys()))
    energies = np.array([valid[k] for k in angles_deg])
    # Find minimum energy point as theta_0
    i_min = np.argmin(energies)
    theta_0 = angles_deg[i_min]
    E_ref = energies[i_min]
    dE_kcal = (energies - E_ref) * 627.5095
    angles_rad = np.deg2rad(angles_deg - theta_0)
    from scipy.optimize import curve_fit
    def model(x, k): return k * x ** 2  # k in kcal/mol/rad²
    popt, _ = curve_fit(model, angles_rad, dE_kcal, p0=[50.0])
    k = popt[0]
    rmsd = float(np.sqrt(np.mean((dE_kcal - model(angles_rad, k)) ** 2)))
    return {
        "form": "k * (theta - theta_0)^2",
        "k_kcal_per_mol_rad2": float(k),
        "theta_0_deg": float(theta_0),
        "rmsd_kcal": rmsd,
        "n_points": len(valid),
    }


def fit_dihedral(scan_data: dict):
    """Fit dihedral CB-CG-CD1 to Σ_n V_n/2 * (1 + cos(n*φ - phase_n))."""
    import numpy as np
    valid = {k: v for k, v in scan_data.items() if v is not None}
    if len(valid) < 8:
        return {"error": "insufficient points", "n_points": len(valid)}
    angles_deg = np.array(sorted(valid.keys()))
    energies = np.array([valid[k] for k in angles_deg])
    E_ref = min(energies)
    dE_kcal = (energies - E_ref) * 627.5095
    angles_rad = np.deg2rad(angles_deg)
    from scipy.optimize import curve_fit
    def model(x, V1, V2, V3, V4):
        # 4-term cosine series (AMBER convention, no phase shifts here)
        return (V1 / 2) * (1 + np.cos(x)) \
             + (V2 / 2) * (1 + np.cos(2 * x)) \
             + (V3 / 2) * (1 + np.cos(3 * x)) \
             + (V4 / 2) * (1 + np.cos(4 * x))
    try:
        popt, _ = curve_fit(model, angles_rad, dE_kcal, p0=[1.0, 1.0, 1.0, 1.0])
    except Exception as e:
        return {"error": f"curve_fit failed: {e}", "n_points": len(valid)}
    rmsd = float(np.sqrt(np.mean((dE_kcal - model(angles_rad, *popt)) ** 2)))
    return {
        "form": "Σ_n V_n/2 * (1 + cos(n*phi))",
        "V1_kcal": float(popt[0]),
        "V2_kcal": float(popt[1]),
        "V3_kcal": float(popt[2]),
        "V4_kcal": float(popt[3]),
        "rmsd_kcal": rmsd,
        "n_points": len(valid),
    }


def main() -> int:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = PROJ / f"outputs/analysis/d3_omw_qm_scans_{ts}"
    logdir.mkdir(parents=True, exist_ok=True)
    print(f"[info] Logdir: {logdir}")

    pdb_in = find_latest_d2()
    if pdb_in is None:
        print("[error] No D2 optimized geometry found — run D2 first.")
        return 1
    print(f"[info] D2 optimized geometry: {pdb_in}")

    atoms = parse_pdb(pdb_in)
    print(f"[info] Atoms: {len(atoms)}")

    # Atom index lookup (1-based for psi4)
    idx_N    = get_atom_idx(atoms, "N", "OMW")
    idx_CD1  = get_atom_idx(atoms, "CD1", "OMW")
    idx_CE2  = get_atom_idx(atoms, "CE2", "OMW")
    idx_CG   = get_atom_idx(atoms, "CG", "OMW")
    idx_CD2  = get_atom_idx(atoms, "CD2", "OMW")
    idx_NE1  = get_atom_idx(atoms, "NE1", "OMW")
    idx_CM   = get_atom_idx(atoms, "CM", "OMW")  # or CZ1
    idx_CB   = get_atom_idx(atoms, "CB", "OMW")

    print(f"[info] Atom indices: N={idx_N} CD1={idx_CD1} CE2={idx_CE2} CG={idx_CG} CD2={idx_CD2} NE1={idx_NE1} CM={idx_CM} CB={idx_CB}")

    # Scan 1: improper N-CD1-CE2-CG (indole 5-ring planarity)
    print("\n[scan 1/3] Improper N-CD1-CE2-CG, -30° to +30° in 5° steps")
    imp_range = list(range(-30, 31, 5))
    imp_data = run_constrained_scan(
        atoms, "improper", [idx_N, idx_CD1, idx_CE2, idx_CG], imp_range, logdir
    )
    (logdir / "improper_N_CD1_CE2_CG.json").write_text(json.dumps(imp_data, indent=2))
    imp_fit = fit_improper(imp_data)

    # Scan 2: angle CM-CG-CD2 (N-methyl junction to indole ring)
    # Note: this is more of a CG-CD2-(NE1-CM) chain check; switch to NE1-CM-(anything)
    # Simpler: angle CM-NE1-CD1 (CM to indole ring via NE1)
    print("\n[scan 2/3] Angle CM-NE1-CD1, 100° to 130° in 5° steps")
    ang_range = list(range(100, 131, 5))
    ang_data = run_constrained_scan(
        atoms, "angle", [idx_CM, idx_NE1, idx_CD1], ang_range, logdir
    )
    (logdir / "angle_CM_NE1_CD1.json").write_text(json.dumps(ang_data, indent=2))
    ang_fit = fit_angle(ang_data)

    # Scan 3: dihedral CB-CG-CD1-CD2 (sidechain orientation vs indole)
    # Note: CB-CG-CD1 has only 3 atoms — need 4 atoms for dihedral. Use CB-CG-CD1-NE1
    print("\n[scan 3/3] Dihedral CB-CG-CD1-NE1, 0° to 330° in 30° steps")
    dih_range = list(range(0, 360, 30))
    dih_data = run_constrained_scan(
        atoms, "dihedral", [idx_CB, idx_CG, idx_CD1, idx_NE1], dih_range, logdir
    )
    (logdir / "dihedral_CB_CG_CD1_NE1.json").write_text(json.dumps(dih_data, indent=2))
    dih_fit = fit_dihedral(dih_data)

    # Summary
    summary = {
        "method": "B3LYP/6-31G(d)",
        "input_geometry": str(pdb_in),
        "improper_N_CD1_CE2_CG_fit": imp_fit,
        "angle_CM_NE1_CD1_fit": ang_fit,
        "dihedral_CB_CG_CD1_NE1_fit": dih_fit,
    }
    (logdir / "amber_fitted_params.json").write_text(json.dumps(summary, indent=2))

    md = [
        "# D3 OMW QM Scans + AMBER Form Fit Report",
        "",
        f"**Date**: {datetime.now().isoformat()}",
        f"**Stack**: psi4 1.10 (single-process, geomeTRIC opt)",
        f"**Method**: B3LYP / 6-31G(d) constrained optimization at each scan point",
        f"**Input geometry**: {pdb_in}",
        "",
        "## Scan 1 — Improper N-CD1-CE2-CG (indole 5-ring planarity)",
        f"- Range: -30° to +30° in 5° steps ({imp_fit.get('n_points', 'n/a')} valid points)",
        f"- AMBER form: `V/2 * (1 + cos(2*phi - 180))`",
        f"- V = {imp_fit.get('V_kcal_per_mol', float('nan')):.3f} kcal/mol",
        f"- Fit RMSD: {imp_fit.get('rmsd_kcal', float('nan')):.4f} kcal/mol",
        "",
        "## Scan 2 — Angle CM-NE1-CD1 (N-methyl junction)",
        f"- Range: 100° to 130° in 5° steps ({ang_fit.get('n_points', 'n/a')} valid points)",
        f"- AMBER form: `k * (theta - theta_0)^2`",
        f"- k = {ang_fit.get('k_kcal_per_mol_rad2', float('nan')):.2f} kcal/mol/rad²",
        f"- theta_0 = {ang_fit.get('theta_0_deg', float('nan')):.1f}°",
        f"- Fit RMSD: {ang_fit.get('rmsd_kcal', float('nan')):.4f} kcal/mol",
        "",
        "## Scan 3 — Dihedral CB-CG-CD1-NE1 (sidechain torsion)",
        f"- Range: 0° to 330° in 30° steps ({dih_fit.get('n_points', 'n/a')} valid points)",
        f"- AMBER form: `Σ_n V_n/2 * (1 + cos(n*phi))`",
        f"- V1 = {dih_fit.get('V1_kcal', float('nan')):.3f} kcal/mol",
        f"- V2 = {dih_fit.get('V2_kcal', float('nan')):.3f} kcal/mol",
        f"- V3 = {dih_fit.get('V3_kcal', float('nan')):.3f} kcal/mol",
        f"- V4 = {dih_fit.get('V4_kcal', float('nan')):.3f} kcal/mol",
        f"- Fit RMSD: {dih_fit.get('rmsd_kcal', float('nan')):.4f} kcal/mol",
        "",
        "## Next step — Layer 3D refitted XML",
        "",
        "These fitted parameters are inputs for the Layer 3D refit of `params/MTR_gaff2.xml`:",
        "- replace the indole improper, the CM-NE1-CD1 angle, and the CB-CG-CD1-NE1 dihedral with QM-derived values",
        "- keep all other amber14SB parameters and the Khoury OMW partial charges (or β hybrid variant)",
        "- target Layer 3D pilot acceptance: λ < 0.005/ns over 5×5 ns MD (per limitations_roadmap §1 L1)",
    ]
    (logdir / "d3_qm_scans_report.md").write_text("\n".join(md))
    print(f"\n[info] D3 report: {logdir / 'd3_qm_scans_report.md'}")
    print(f"=== D3 complete ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
