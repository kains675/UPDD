#!/usr/bin/env python
"""
utils/run_benchmark_2qki.py
---------------------------
Benchmark runner: Compstatin 4W9A + C3c (PDB 2QKI) full pipeline validation.

Purpose
~~~~~~~
SciVal Verdict 20 benchmark — run UPDD's MD + snapshot + MM-GBSA chain on a
literature-characterized cyclic peptide (compstatin 4W9A, Kd 152 nM, Magotti
2009 SPR) to verify whether the protocol returns physically reasonable
ΔG_bind (literature = −9.30 kcal/mol).

If UPDD's end-to-end result is reasonable (−5 to −20 kcal/mol) → protocol
is sound, KRAS-GNP ΔG inflation is target-specific. If result is
non-physical (<-100 or worse) → protocol-level bug, must fix before KRAS
re-entry.

Pipeline (abbreviated vs. UPDD.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. MD: 2QKI_clean.pdb (already has H, ions stripped in preprocess) →
   TIP3P 0.15 M NaCl (extracellular, C3 is plasma target), cyclic_ss
   topology, dt 2fs, 2ns production (abbreviated for benchmark speed).
2. Snapshot: 5 frames uniform from production.
3. MM-GBSA: each snapshot → complex/receptor/ligand split → 500-iter
   minimize → ΔG.
4. Compare mean ⟨ΔG⟩ to literature −9.30 kcal/mol.

Skipped
~~~~~~~
- Steps 1-8 (RFdiff / MPNN / AF2 / ncAA mutation / parameterize): no
  design needed, crystal structure has native compstatin binder.
- Step 11 (QM/MM): not part of MM-GBSA validation chain.
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
UPDD_ROOT = os.path.abspath(os.path.join(HERE, os.pardir))


def run_or_die(cmd, label):
    print(f"\n{'=' * 60}\n{label}\n{'=' * 60}")
    print("$ " + " ".join(cmd))
    r = subprocess.run(cmd, check=False)
    if r.returncode != 0:
        sys.exit(f"[!] {label} failed (rc={r.returncode})")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--steps", type=int, default=1_000_000,
                    help="MD production steps (2 ns default at 2fs; set "
                         "5_000_000 for 10 ns). Benchmark speed vs accuracy.")
    ap.add_argument("--n_snapshots", type=int, default=5)
    ap.add_argument("--outputdir", default=os.path.join(UPDD_ROOT, "outputs", "2QKI_benchmark"))
    ap.add_argument("--pdb", default=os.path.join(UPDD_ROOT, "target", "2QKI_clean.pdb"))
    args = ap.parse_args()

    outdir = os.path.abspath(args.outputdir)
    md_dir = os.path.join(outdir, "mdresult")
    snap_dir = os.path.join(outdir, "snapshots")
    mmgbsa_dir = os.path.join(outdir, "mmgbsa_results")
    for d in (md_dir, snap_dir, mmgbsa_dir):
        os.makedirs(d, exist_ok=True)

    # Stage PDB into an input dir (run_restrained_md expects inputdir/*.pdb).
    md_input_dir = os.path.join(outdir, "_md_input")
    os.makedirs(md_input_dir, exist_ok=True)
    staged_pdb = os.path.join(md_input_dir, "2QKI_compstatin_4W9A.pdb")
    shutil.copy2(args.pdb, staged_pdb)
    print(f"[benchmark] staged {args.pdb} → {staged_pdb}")

    # Step 1: MD (cyclic_ss, extracellular Na+/Cl- — target_card.subcellular_localization=extracellular)
    md_cmd = [
        "/home/san/miniconda3/envs/md_simulation/bin/python",
        os.path.join(HERE, "run_restrained_md.py"),
        "--inputdir", md_input_dir,
        "--outputdir", md_dir,
        "--steps", str(args.steps),
        "--topology", "cyclic_ss",
        "--binder_chain", "B",
        "--graph_policy", "strict",
        "--target_id", "2QKI",
        "--ncaa_label", "none", "--ncaa_code", "",
        "--dt_fs", "2.0",
        "--platform", "CUDA",
    ]
    run_or_die(md_cmd, "Step 9: Restrained MD (Compstatin 4W9A + C3c)")

    # Step 2: snapshot extract
    snap_cmd = [
        "/home/san/miniconda3/envs/md_simulation/bin/python",
        os.path.join(HERE, "extract_snapshots.py"),
        "--md_dir", md_dir,
        "--outputdir", snap_dir,
        "--n_snapshots", str(args.n_snapshots),
        "--binder_chain", "B",
        "--target_id", "2QKI",
    ]
    run_or_die(snap_cmd, "Step 10: Snapshot Extraction")

    # Step 3: MM-GBSA on snapshots
    mmgbsa_cmd = [
        "/home/san/miniconda3/envs/md_simulation/bin/python",
        os.path.join(HERE, "run_mmgbsa.py"),
        "--md_dir", snap_dir,
        "--outputdir", mmgbsa_dir,
        "--ncaa_elem", "none",
        "--receptor_chain", "A",
        "--binder_chain", "B",
        "--target_id", "2QKI",
    ]
    run_or_die(mmgbsa_cmd, "Step 12: MM-GBSA ΔG calculation")

    # Aggregate results
    print(f"\n{'=' * 60}\nBENCHMARK RESULT\n{'=' * 60}")
    results = []
    for p in sorted(glob.glob(os.path.join(mmgbsa_dir, "*.json"))):
        try:
            with open(p, encoding="utf-8") as f:
                d = json.load(f)
            if "delta_g_kcal" in d:
                results.append(d)
        except (OSError, json.JSONDecodeError):
            continue

    if not results:
        print("[!] no MM-GBSA JSONs produced — check individual step logs")
        sys.exit(1)

    dgs = [r["delta_g_kcal"] for r in results]
    mean = sum(dgs) / len(dgs)
    stdev = (sum((x - mean) ** 2 for x in dgs) / max(len(dgs), 1)) ** 0.5
    print(f"  snapshots:       {len(dgs)}")
    print(f"  ΔG range:        {min(dgs):+.2f} to {max(dgs):+.2f} kcal/mol")
    print(f"  ⟨ΔG⟩ mean:       {mean:+.2f} kcal/mol")
    print(f"  stdev:           {stdev:.2f} kcal/mol")
    print(f"  literature ΔG:   -9.30 kcal/mol (Magotti 2009 SPR, 4W9A Kd=152 nM)")
    print(f"  discrepancy:     {mean - (-9.30):+.2f} kcal/mol")
    if -20 <= mean <= -5:
        print(f"  ✅ PASS: protocol gives physically reasonable ΔG")
    elif mean < -100:
        print(f"  🔴 FAIL: non-physical magnitude, protocol-level bug confirmed")
    else:
        print(f"  🟡 WARN: outside strict physical range but not catastrophic")


if __name__ == "__main__":
    main()
