#!/usr/bin/env python
"""
D1 Psi4 stack Khoury OMW RESP-A2 reproduction via psi4numpy resp_driver.

Pure Python single-process pipeline (no qcfractal / external scripts):
  - psi4 HF/6-31G* single-point on Ace-OMW-NMe dipeptide
  - Merz-Kollman ESP grid (vdw_scale=[1.4, 1.6, 1.8, 2.0], density=1.0/Å²)
  - RESP-A2 2-stage fit (a=0.0005 / 0.001, b=0.1, methyl 3-equivalent)

Reference: psi4numpy/One-Electron-Property/Restrained-Electrostatic-Potential
(Alenaizan, BSD-3-Clause). resp_driver.py / resp_helper.py / espfit.py copied
into scripts/psi4numpy_resp/.

Outputs: outputs/analysis/d1_psi4_resp_reproduction_{TS}/
  - omw_dipeptide.pdb  (geometry reused from PySCF Coder Pass B)
  - omw_psi4_charges.json  (27 OMW atoms, UPDD naming)
  - d1_psi4_reproduction_report.md  (vs Khoury published)

Run:
  /home/san/miniconda3/envs/psi4/bin/python scripts/d1_khoury_psi4_resp.py
"""
from __future__ import annotations

import json
import shutil
import sys
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

PROJ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ / "scripts" / "psi4numpy_resp"))

# Khoury 2014 OMW published RESP-A2 charges (UPDD naming: CZ1→CM, HZ1*→HM*).
KHOURY_OMW_CHARGES = {
    "N":   -0.280943, "H":   0.241023, "CA": -0.124259, "HA":  0.120930,
    "C":    0.575115, "O":  -0.571038, "CB": -0.157536, "HB2": 0.104569,
    "HB3":  0.104569, "CG":  -0.126422, "CD2": 0.083375, "CE2": 0.119586,
    "NE1":  0.077670, "CD1": -0.252655, "HD1": 0.204947, "CM": -0.239062,
    "HM1":  0.112200, "HM2":  0.112200, "HM3":  0.112200, "CZ2": -0.241414,
    "HZ2":  0.146305, "CH2": -0.139177, "HH2":  0.136546, "CZ3": -0.206367,
    "HZ3":  0.143250, "CE3": -0.211658, "HE3":  0.156047,
}


def find_pyscf_dipeptide_pdb() -> Path | None:
    for d in sorted((PROJ / "outputs/analysis").glob("d1_khoury_resp_reproduction_*"), reverse=True):
        pdb = d / "omw_dipeptide.pdb"
        if pdb.exists():
            return pdb
    return None


def parse_pdb(pdb_path: Path) -> list[tuple[str, str, str, float, float, float]]:
    """Parse ATOM/HETATM records. Returns [(atom_name, residue, element, x, y, z), ...]"""
    out = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            name = line[12:16].strip()
            residue = line[17:20].strip()
            element = line[76:78].strip() or name[0]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            out.append((name, residue, element, x, y, z))
    return out


def build_psi4_geom(atoms):
    """Build psi4 Molecule from parsed PDB atoms."""
    import psi4
    geom_lines = []
    for atom in atoms:
        name, residue, element, x, y, z = atom
        geom_lines.append(f"{element:<3s} {x:14.8f} {y:14.8f} {z:14.8f}")
    geom_block = "\n".join(geom_lines) + "\nno_reorient\nno_com\nunits angstrom\n"
    mol = psi4.geometry(geom_block)
    mol.set_name("Ace-OMW-NMe")
    mol.update_geometry()
    return mol


def identify_omw_methyl_atoms(atoms):
    """Find CM (N-methyl carbon of OMW) and HM1/HM2/HM3 atom indices (0-based).
    Filter by residue == 'OMW' to avoid Ace/NMe atom-name collision.
    """
    cm_idx = None
    hm_indices = []
    for i, atom in enumerate(atoms):
        name, residue = atom[0], atom[1]
        if residue != "OMW":
            continue
        if name in ("CM", "CZ1"):
            cm_idx = i
        elif name in ("HM1", "HZ11"):
            hm_indices.append(i)
        elif name in ("HM2", "HZ12"):
            hm_indices.append(i)
        elif name in ("HM3", "HZ13"):
            hm_indices.append(i)
    return cm_idx, sorted(hm_indices)


def identify_omw_methylene_atoms(atoms):
    hb = []
    for i, atom in enumerate(atoms):
        name, residue = atom[0], atom[1]
        if residue == "OMW" and name in ("HB2", "HB3"):
            hb.append(i)
    return sorted(hb)


def identify_cap_atoms(atoms: list[tuple[str, ...]]) -> list[int]:
    """ACE/NME cap atom indices to constrain to neutral sum in stage 2."""
    cap_names = {
        # Ace cap (typical RDKit/PSI4 ordering may vary, identify by name pattern)
        "CH3", "CT", "C1", "C2", "C3",  # Ace methyl C
        "H1", "H2", "H3",                # Ace methyl H (or NME methyl H)
        # NME cap atoms not always identifiable by name alone
    }
    cap = []
    for i, (name, *_) in enumerate(atoms):
        if name in cap_names:
            cap.append(i)
    return sorted(cap)


def run_resp_fit(atoms, logdir: Path) -> dict[str, float]:
    """Run psi4numpy resp_driver 2-stage RESP-A2 fit."""
    import psi4
    import numpy as np
    import resp_driver

    psi4.set_output_file(str(logdir / "psi4_resp.dat"), True)
    psi4.set_memory("4 GB")
    psi4.core.set_num_threads(4)

    mol = build_psi4_geom(atoms)
    n = mol.natom()
    print(f"[info] Molecule built: {n} atoms")

    # Stage 1: weak restraint, no constraints
    stage1_options = {
        "N_VDW_LAYERS":      4,
        "VDW_SCALE_FACTOR":  1.4,
        "VDW_INCREMENT":     0.2,
        "VDW_POINT_DENSITY": 1.0,
        "BASIS_ESP":         "6-31g*",
        "METHOD_ESP":        "hf",
        "resp_a":            0.0005,
        "RESP_B":            0.1,
        "RESTRAINT":         True,
    }

    print("[info] Stage 1 RESP-A2 fit (HF/6-31G*, MK 4-layer, a=0.0005)...")
    charges1 = resp_driver.resp([mol], [stage1_options])
    stage1_esp = charges1[0][0]
    stage1_resp = charges1[0][1]
    print(f"[info]   stage 1 charges sum: {sum(stage1_resp):+.4f}")

    # Stage 2: constraint_charge fixes everything except methyl + methylene H;
    # constraint_group enforces methyl 3-H equivalence + methylene 2-H equivalence.
    cm_idx, hm_indices = identify_omw_methyl_atoms(atoms)
    hb_indices = identify_omw_methylene_atoms(atoms)

    print(f"[info]   Stage 2 equivalence groups:")
    print(f"   methyl HM (1-based): {[i+1 for i in hm_indices]}")
    print(f"   methylene HB (1-based): {[i+1 for i in hb_indices]}")

    stage2_options = dict(stage1_options)
    stage2_options["resp_a"] = 0.001
    # RESP-A2 standard stage 2: ALL atoms free re-fit under stronger restraint, with
    # equivalence groups (methyl 3-H + methylene 2-H) enforced. NO constraint_charge
    # (don't freeze stage 1 values — that defeats the strong-restraint refit purpose).
    constraint_group = []
    if len(hm_indices) >= 2:
        constraint_group.append([i + 1 for i in hm_indices])
    if len(hb_indices) >= 2:
        constraint_group.append([i + 1 for i in hb_indices])
    stage2_options["constraint_group"] = constraint_group

    print("[info] Stage 2 RESP-A2 fit (a=0.001, equivalence constraints)...")
    charges2 = resp_driver.resp([mol], [stage2_options])
    stage2_resp = charges2[0][1]
    print(f"[info]   stage 2 charges sum: {sum(stage2_resp):+.4f}")

    # Map indices → atom names (OMW residue only, with Khoury → UPDD remap on CZ1/HZ1*)
    KH_TO_UPDD = {"CZ1": "CM", "HZ11": "HM1", "HZ12": "HM2", "HZ13": "HM3"}
    final = {}
    for i, atom in enumerate(atoms):
        name, residue = atom[0], atom[1]
        if residue != "OMW":
            continue
        updd_name = KH_TO_UPDD.get(name, name)
        if updd_name in KHOURY_OMW_CHARGES:
            final[updd_name] = float(stage2_resp[i])
    return final


def compare_to_khoury(charges: dict[str, float], report_path: Path) -> tuple[str, float, float]:
    rmsd_sq = 0.0
    max_abs = 0.0
    rows = []
    n = 0
    for name, khoury_q in KHOURY_OMW_CHARGES.items():
        psi4_q = charges.get(name)
        if psi4_q is None:
            rows.append((name, khoury_q, None, None))
            continue
        delta = psi4_q - khoury_q
        rmsd_sq += delta ** 2
        max_abs = max(max_abs, abs(delta))
        rows.append((name, khoury_q, psi4_q, delta))
        n += 1
    rmsd = (rmsd_sq / n) ** 0.5 if n else float("nan")
    verdict = "PASS" if max_abs < 0.005 else "FAIL"

    md = [
        "# D1 Psi4 (psi4numpy resp_driver) Khoury OMW RESP-A2 Reproduction Report",
        "",
        f"**Date**: {datetime.now().isoformat()}",
        f"**Stack**: psi4 1.10 + psi4numpy resp_driver (Alenaizan, BSD-3-Clause)",
        f"**Protocol**: HF/6-31G* + MK ESP (4 layers @ 1.4/1.6/1.8/2.0 × vdW, density 1.0/Å²) + RESP-A2 2-stage (a=0.0005/0.001, b=0.1, methyl + methylene equivalence)",
        f"**Reference**: Khoury 2014 ACS Synth Biol DOI 10.1021/sb400168u (ffncaa.in L778-830)",
        "",
        f"## Verdict: **{verdict}**",
        f"- max |Δq| = {max_abs:.4f} e (threshold 0.005 e per atom)",
        f"- RMSD     = {rmsd:.4f} e",
        f"- Compared = {n} / {len(KHOURY_OMW_CHARGES)} atoms",
        "",
        "## Per-atom comparison",
        "",
        "| Atom | Khoury (e) | Psi4 (e) | Δq (e) | |Δq| > 0.005? |",
        "|---|---|---|---|---|",
    ]
    for name, kq, pq, dq in rows:
        pq_s = f"{pq:+.6f}" if pq is not None else "n/a"
        dq_s = f"{dq:+.4f}" if dq is not None else "n/a"
        flag = "**FAIL**" if (dq is not None and abs(dq) > 0.005) else "ok"
        md.append(f"| {name} | {kq:+.6f} | {pq_s} | {dq_s} | {flag} |")

    report_path.write_text("\n".join(md))
    return verdict, max_abs, rmsd


def main() -> int:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logdir = PROJ / f"outputs/analysis/d1_psi4_resp_reproduction_{ts}"
    logdir.mkdir(parents=True, exist_ok=True)
    print(f"[info] Logdir: {logdir}")

    pdb = find_pyscf_dipeptide_pdb()
    if pdb is None:
        print("[error] No PySCF dipeptide.pdb found — run PySCF D1 first.")
        return 1
    print(f"[info] Geometry reused: {pdb}")
    dipep_pdb = logdir / "omw_dipeptide.pdb"
    shutil.copy(pdb, dipep_pdb)

    atoms = parse_pdb(dipep_pdb)
    print(f"[info] Atoms parsed: {len(atoms)}")

    charges = run_resp_fit(atoms, logdir)

    (logdir / "omw_psi4_charges.json").write_text(json.dumps(charges, indent=2))
    print(f"[info] Charges saved: {logdir / 'omw_psi4_charges.json'}")

    verdict, max_abs, rmsd = compare_to_khoury(
        charges, logdir / "d1_psi4_reproduction_report.md"
    )
    print()
    print(f"=== D1 Psi4 reproduction (psi4numpy): {verdict} ===")
    print(f"max |Δq| = {max_abs:.4f} e")
    print(f"RMSD     = {rmsd:.4f} e")
    print(f"Report   = {logdir / 'd1_psi4_reproduction_report.md'}")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
