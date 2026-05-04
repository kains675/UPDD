#!/usr/bin/env python
"""
Test: re-image the existing snapshot PDBs with image_molecules and re-measure
binder integrity. If image_molecules+make_whole brings VAL5-C..MTR-N to ~1.3 Å,
that confirms the σ_w=38 was a snapshot-extraction PBC artifact, NOT a real
trajectory issue.

For each of the 25 snapshots:
  - Read PDB
  - Apply image_molecules with binder chain as one molecule
  - Measure VAL5-C..MTR-N and MTR-C..TRP7-N before and after
  - Save corrected snapshots to a sibling directory `snapshots_n25_reimaged/`
"""
import json
import os
from pathlib import Path
import numpy as np
import mdtraj as md

ROOT = Path('/home/san/UPDD_proj')
SYS_DIR = ROOT / 'outputs' / '7TL8_MTR6_calib_s7'
SNAP_DIR = SYS_DIR / 'snapshots_n25'
EXT_DIR = SYS_DIR / 'mdresult' / 'extension_audit_20260428'
EXT_DIR.mkdir(parents=True, exist_ok=True)
REIMAGED_DIR = SYS_DIR / 'snapshots_n25_reimaged'
REIMAGED_DIR.mkdir(parents=True, exist_ok=True)

records = []
for pdb in sorted(SNAP_DIR.glob('7TL8_MTR6_snap*.pdb')):
    t = md.load(str(pdb))
    top = t.topology

    # Identify binder/receptor by MTR
    mtr = next((r for r in top.residues if r.name == 'MTR'), None)
    if not mtr:
        continue
    binder_chain_idx = mtr.chain.index
    PROT_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU',
               'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','HID','HIE','HIP','MTR'}
    recv_chain_idx = next(ch.index for ch in top.chains
                           if ch.index != binder_chain_idx
                           and len([r for r in ch.residues if r.name in PROT_AA]) > 50)

    # Indices
    val5 = next(r for r in mtr.chain.residues if r.resSeq == 5 and r.name == 'VAL')
    trp7 = next(r for r in mtr.chain.residues if r.resSeq == 7 and r.name == 'TRP')
    val5_C = next(a for a in val5.atoms if a.name == 'C')
    mtr_N = next(a for a in mtr.atoms if a.name == 'N')
    mtr_C = next(a for a in mtr.atoms if a.name == 'C')
    trp7_N = next(a for a in trp7.atoms if a.name == 'N')

    # Pre
    d_pre_v_m = float(np.linalg.norm(t.xyz[0, val5_C.index, :] - t.xyz[0, mtr_N.index, :]) * 10)
    d_pre_m_t = float(np.linalg.norm(t.xyz[0, mtr_C.index, :] - t.xyz[0, trp7_N.index, :]) * 10)

    # Re-image
    try:
        recv_atoms = set(top.chain(recv_chain_idx).atoms)
        binder_atoms = set(top.chain(binder_chain_idx).atoms)
        t2 = t.image_molecules(inplace=False, anchor_molecules=[recv_atoms],
                                other_molecules=[binder_atoms], make_whole=True)
        d_post_v_m = float(np.linalg.norm(t2.xyz[0, val5_C.index, :] - t2.xyz[0, mtr_N.index, :]) * 10)
        d_post_m_t = float(np.linalg.norm(t2.xyz[0, mtr_C.index, :] - t2.xyz[0, trp7_N.index, :]) * 10)
        # Save re-imaged PDB
        out_path = REIMAGED_DIR / pdb.name
        t2.save_pdb(str(out_path))
    except Exception as e:
        d_post_v_m = float('nan')
        d_post_m_t = float('nan')
        out_path = None
        print(f"  {pdb.name}: image_molecules failed: {e}")

    rec = {
        'snap': pdb.stem,
        'pre_VAL5C_MTRN_A': d_pre_v_m,
        'pre_MTRC_TRP7N_A': d_pre_m_t,
        'post_VAL5C_MTRN_A': d_post_v_m,
        'post_MTRC_TRP7N_A': d_post_m_t,
        'reimaged_pdb': str(out_path) if out_path else None,
    }
    records.append(rec)
    print(f"  {pdb.stem}: pre {d_pre_v_m:.2f}/{d_pre_m_t:.2f}  -> post {d_post_v_m:.2f}/{d_post_m_t:.2f} A")

# Summary
print("\n=== SUMMARY ===")
pre_vmn = np.array([r['pre_VAL5C_MTRN_A'] for r in records])
post_vmn = np.array([r['post_VAL5C_MTRN_A'] for r in records])
print(f"VAL5-C..MTR-N before image_molecules: median={np.median(pre_vmn):.2f}, max={pre_vmn.max():.2f} A")
print(f"VAL5-C..MTR-N after image_molecules:  median={np.median(post_vmn):.2f}, max={post_vmn.max():.2f} A")
print(f"Snapshots with bond intact (post < 2 A): {int((post_vmn < 2.0).sum())}/{len(records)}")

with open(EXT_DIR / 'snapshot_imaging_test.json', 'w') as f:
    json.dump({
        'n': len(records),
        'pre_VAL5C_MTRN_A_median': float(np.median(pre_vmn)),
        'post_VAL5C_MTRN_A_median': float(np.median(post_vmn)),
        'n_intact_after_imaging': int((post_vmn < 2.0).sum()),
        'records': records,
    }, f, indent=2)
print(f"\n=> Saved -> {EXT_DIR/'snapshot_imaging_test.json'}")
print(f"   Re-imaged snapshots → {REIMAGED_DIR}")
