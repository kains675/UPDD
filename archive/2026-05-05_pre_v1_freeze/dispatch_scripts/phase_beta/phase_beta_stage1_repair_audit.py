#!/usr/bin/env python
"""
Phase β Stage 1-Repair — 7TL8_MTR6_s7 σ_w=38.33 mechanism audit.

Inputs:
  - outputs/7TL8_MTR6_calib_s7/mdresult/7TL8_MTR6_restrained.dcd  (865 MB, 500 frames)
  - outputs/7TL8_MTR6_calib_s7/mdresult/7TL8_MTR6_final.pdb       (topology)
  - outputs/7TL8_MTR6_calib_s7/snapshots_n25/snap{NN}.pdb         (25 PBSA snapshots)
  - outputs/7TL8_MTR6_calib_s7/mmpbsa_results_fix4_postpatch/snap_NN_*.json  (per-snap energies)

Outputs (under outputs/7TL8_MTR6_calib_s7/mdresult/extension_audit_20260428/):
  - audit_summary.json   (all metrics)
  - rmsd_backbone_binder.csv      (per-frame Cα/N/C/O RMSD on chain B)
  - rmsd_mtr_sidechain.csv        (per-frame MTR side-chain heavy atom RMSD)
  - rmsf_per_residue_binder.csv   (per-residue RMSF on chain B)
  - rmsf_mtr_atoms.csv            (per-atom RMSF for MTR ring/methyl atoms)
  - mtr_dihedrals.csv             (χ1, χ2 per frame for rotamer clustering)
  - per_snap_energy_pose.csv      (joins snapshot frame index to ΔG_PBSA)

Mechanism classification:
  A) Cyclic-peptide / backbone instability: backbone RMSD high, drift over time.
  B) MTR rotamer flipping: MTR side-chain RMSF high, multimodal χ1/χ2.
  C) Binding-pose drift: binder COM drift relative to receptor (chain A) frame.
  D) FF artifact: energy outlier snapshots align with structural anomalies (high
     |F| or torsion strain) rather than equilibrium fluctuations.

Run from project root:
  /home/san/miniconda3/envs/qmmm/bin/python scripts/phase_beta_stage1_repair_audit.py
"""

import os
import json
import math
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import mdtraj as md

# ---------- Paths ----------
ROOT = Path('/home/san/UPDD_proj')
SYS_DIR = ROOT / 'outputs' / '7TL8_MTR6_calib_s7'
DCD = SYS_DIR / 'mdresult' / '7TL8_MTR6_restrained.dcd'
TOP_PDB = SYS_DIR / 'mdresult' / '7TL8_MTR6_final.pdb'
EXT_DIR = SYS_DIR / 'mdresult' / 'extension_audit_20260428'
EXT_DIR.mkdir(parents=True, exist_ok=True)
PBSA_DIR = SYS_DIR / 'mmpbsa_results_fix4_postpatch'

print(f"[audit] ROOT={ROOT}")
print(f"[audit] DCD={DCD} ({DCD.stat().st_size / 1e6:.1f} MB)")
print(f"[audit] TOP={TOP_PDB}")
print(f"[audit] EXT_DIR={EXT_DIR}")

# ---------- Load trajectory ----------
print("[audit] loading trajectory ...")
traj = md.load_dcd(str(DCD), top=str(TOP_PDB))
print(f"[audit] n_frames={traj.n_frames}, n_atoms={traj.n_atoms}")
top = traj.topology

# ---------- Identify chains / MTR ----------
# Chain B = binder, includes MTR at residue 6 (per pre-flight). Chain A = receptor.
# In mdtraj's topology after solvation, ions/water are added — chain mapping changes.
# We rely on resnames + residue order: peptide chains (A protein, B with MTR).
binder_resnames = ['TYR', 'GLN', 'VAL', 'THR', 'VAL', 'MTR', 'TRP', 'ALA',
                   'PRO', 'TRP', 'GLU', 'ASP', 'CYS']

# Identify MTR residue (must exist exactly once in protein-like residues).
mtr_residues = [r for r in top.residues if r.name == 'MTR']
print(f"[audit] MTR residues found: {len(mtr_residues)}  -> "
      f"{[(r.chain.index, r.index, r.resSeq) for r in mtr_residues]}")
assert len(mtr_residues) == 1, "Expected exactly one MTR residue"
mtr_res = mtr_residues[0]

# Identify chain B = MTR's chain. Use its chain index for binder atom selection.
binder_chain_idx = mtr_res.chain.index
binder_residues = list(mtr_res.chain.residues)
print(f"[audit] binder chain idx={binder_chain_idx}, n_residues={len(binder_residues)}")
print(f"[audit] binder resnames: {[r.name for r in binder_residues]}")

# Chain A = receptor (protein chain that is NOT the binder chain and is amino-acid).
PROT_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU',
           'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','HID','HIE','HIP','MTR'}
receptor_chain_idx = None
for ch in top.chains:
    if ch.index == binder_chain_idx:
        continue
    rs = [r for r in ch.residues if r.name in PROT_AA]
    if len(rs) > 50:
        receptor_chain_idx = ch.index
        break
print(f"[audit] receptor chain idx={receptor_chain_idx}")
receptor_residues = [r for r in top.chain(receptor_chain_idx).residues if r.name in PROT_AA]
print(f"[audit] receptor n_residues={len(receptor_residues)}")

# ---------- Atom selections ----------
# Backbone heavy atoms of binder chain (Cα, C, N, O) — for RMSD/RMSF.
binder_bb = top.select(
    f"chainid {binder_chain_idx} and (name CA or name C or name N or name O)"
)
print(f"[audit] binder backbone atoms: {len(binder_bb)}")

# Receptor backbone (for alignment frame).
recv_bb = top.select(
    f"chainid {receptor_chain_idx} and (name CA or name C or name N or name O)"
)
print(f"[audit] receptor backbone atoms: {len(recv_bb)}")

# All-protein heavy atoms (binder + receptor) for global alignment.
prot_chains_idxs = (binder_chain_idx, receptor_chain_idx)
prot_heavy = top.select(
    f"(chainid {binder_chain_idx} or chainid {receptor_chain_idx}) and not symbol H"
)
print(f"[audit] all-protein heavy atoms: {len(prot_heavy)}")

# MTR side-chain heavy atoms — MTR atoms minus N, CA, C, O (backbone) and H*.
mtr_atom_names = sorted({a.name for a in mtr_res.atoms})
print(f"[audit] MTR atom names ({len(mtr_atom_names)}): {mtr_atom_names}")
mtr_sc_idx = []
mtr_sc_names = []
for a in mtr_res.atoms:
    if a.element.symbol == 'H':
        continue
    if a.name in {'N', 'CA', 'C', 'O', 'OXT'}:
        continue
    mtr_sc_idx.append(a.index)
    mtr_sc_names.append(a.name)
mtr_sc_idx = np.array(mtr_sc_idx, dtype=int)
print(f"[audit] MTR side-chain heavy atoms: {len(mtr_sc_idx)}: {mtr_sc_names}")

# ---------- Pre-compute binder atom indices used for PBC unwrap ----------
binder_heavy = top.select(f"chainid {binder_chain_idx} and not symbol H")
binder_chain_atoms_all = np.array(
    [a.index for a in top.chain(binder_chain_idx).atoms], dtype=int
)
print(f"[audit] binder heavy atoms (for COM): {len(binder_heavy)}")
print(f"[audit] binder chain all atoms (for unwrap): {len(binder_chain_atoms_all)}")

# ---------- Alignment ----------
# (1) image_molecules makes molecules whole and centers anchor (receptor)
print("[audit] image_molecules (anchor=receptor) ...")
recv_atom_set = set(top.chain(receptor_chain_idx).atoms)
binder_atom_set = set(top.chain(binder_chain_idx).atoms)
try:
    traj = traj.image_molecules(
        inplace=False,
        anchor_molecules=[recv_atom_set],
        other_molecules=[binder_atom_set],
        make_whole=True,
    )
    print("[audit] image_molecules ok")
except Exception as exc:
    print(f"[audit][warn] image_molecules failed ({exc})")

# (2) Per-frame PBC-correction: shift binder chain so its COM lies in the
#     receptor's primary image (delta_COM reduced to nearest image).
print("[audit] PBC nearest-image correction on binder COM ...")
box_lengths = traj.unitcell_lengths
n_shifted = 0
for i in range(traj.n_frames):
    L = box_lengths[i]
    if not np.all(L > 0):
        continue
    recv_com = traj.xyz[i, recv_bb, :].mean(axis=0)
    bind_com = traj.xyz[i, binder_heavy, :].mean(axis=0)
    delta = bind_com - recv_com
    shift = np.round(delta / L) * L
    if np.any(shift != 0.0):
        traj.xyz[i, binder_chain_atoms_all, :] -= shift
        n_shifted += 1
print(f"[audit] PBC-shifted frames: {n_shifted}/{traj.n_frames}")

# (3) Pick a "canonical bound" reference frame: one whose binder COM is
#     closest to the median binder COM (so we don't anchor to an outlier).
print("[audit] selecting bound-pose reference frame ...")
all_bind_com = np.array([traj.xyz[i, binder_heavy, :].mean(axis=0)
                         for i in range(traj.n_frames)])
all_recv_com = np.array([traj.xyz[i, recv_bb, :].mean(axis=0)
                         for i in range(traj.n_frames)])
delta_com = all_bind_com - all_recv_com
median_delta = np.median(delta_com, axis=0)
distances = np.linalg.norm(delta_com - median_delta, axis=1)
ref_frame = int(np.argmin(distances))
print(f"[audit] ref_frame={ref_frame}  median_delta_COM={median_delta*10} A")

# (4) Align all frames to that reference by receptor backbone
print(f"[audit] aligning trajectory by receptor backbone (frame {ref_frame} ref) ...")
traj.superpose(traj, frame=ref_frame, atom_indices=recv_bb)

# ---------- RMSD: backbone of binder (ref_frame) ----------
print("[audit] computing binder backbone RMSD ...")
ref_xyz = traj.xyz[ref_frame, binder_bb, :]
rmsd_bb_binder = np.zeros(traj.n_frames, dtype=float)
for i in range(traj.n_frames):
    diff = traj.xyz[i, binder_bb, :] - ref_xyz
    rmsd_bb_binder[i] = math.sqrt(float(np.mean(np.sum(diff * diff, axis=1)))) * 10.0

# ---------- RMSD: MTR side chain (ref_frame) ----------
print("[audit] computing MTR side-chain heavy-atom RMSD ...")
ref_sc = traj.xyz[ref_frame, mtr_sc_idx, :]
rmsd_mtr_sc = np.zeros(traj.n_frames, dtype=float)
for i in range(traj.n_frames):
    diff = traj.xyz[i, mtr_sc_idx, :] - ref_sc
    rmsd_mtr_sc[i] = math.sqrt(float(np.mean(np.sum(diff * diff, axis=1)))) * 10.0

# ---------- Binder COM drift relative to receptor frame ----------
print("[audit] computing binder COM drift ...")
masses = np.array([top.atom(int(i)).element.mass for i in binder_heavy], dtype=float)

def com(xyz_frame, atom_idxs, masses_):
    return (masses_[:, None] * xyz_frame[atom_idxs, :]).sum(axis=0) / masses_.sum()

ref_com = com(traj.xyz[ref_frame], binder_heavy, masses)
binder_com_drift = np.zeros(traj.n_frames, dtype=float)
for i in range(traj.n_frames):
    c = com(traj.xyz[i], binder_heavy, masses)
    binder_com_drift[i] = float(np.linalg.norm(c - ref_com)) * 10.0

# ---------- Per-residue RMSF (binder) ----------
print("[audit] computing per-residue binder RMSF ...")
mean_xyz = traj.xyz.mean(axis=0)  # (n_atoms, 3)
diffs = traj.xyz - mean_xyz[None, :, :]  # (n_frames, n_atoms, 3)
sq = (diffs * diffs).sum(axis=2)  # (n_frames, n_atoms)
rmsf_atom = np.sqrt(sq.mean(axis=0)) * 10.0  # Å, shape (n_atoms,)

per_res_rmsf = []
for r in binder_residues:
    heavy = [a.index for a in r.atoms if a.element.symbol != 'H']
    if not heavy:
        per_res_rmsf.append((r.resSeq, r.name, float('nan'), 0))
        continue
    val = float(rmsf_atom[heavy].mean())
    per_res_rmsf.append((r.resSeq, r.name, val, len(heavy)))

# ---------- Per-atom RMSF for MTR ----------
mtr_atom_rmsf = []
for a in mtr_res.atoms:
    if a.element.symbol == 'H':
        continue
    mtr_atom_rmsf.append((a.name, float(rmsf_atom[a.index])))

# ---------- MTR χ1, χ2 dihedrals ----------
print("[audit] computing MTR chi1, chi2 ...")
# MTR (modified Trp) atom layout, derived from name list above:
# Backbone N, CA, C, O. Side-chain: CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2.
# CM = N1 methyl carbon (on NE1 nitrogen).
# χ1 = N - CA - CB - CG
# χ2 = CA - CB - CG - CD1
def name2idx(res, name):
    for a in res.atoms:
        if a.name == name:
            return a.index
    raise KeyError(f"Atom {name} not found in {res.name} {res.resSeq}")

chi1_atoms = [name2idx(mtr_res, n) for n in ('N', 'CA', 'CB', 'CG')]
try:
    chi2_atoms = [name2idx(mtr_res, n) for n in ('CA', 'CB', 'CG', 'CD1')]
except KeyError:
    # fallback: CD2
    chi2_atoms = [name2idx(mtr_res, n) for n in ('CA', 'CB', 'CG', 'CD2')]

chi1 = md.compute_dihedrals(traj, np.array([chi1_atoms]))[:, 0] * (180.0 / math.pi)
chi2 = md.compute_dihedrals(traj, np.array([chi2_atoms]))[:, 0] * (180.0 / math.pi)

# Simple grid clustering on (χ1, χ2): bin into 60° cells
def grid_cluster(c1, c2, bin_deg=60.0):
    b1 = np.floor((c1 + 180.0) / bin_deg).astype(int)
    b2 = np.floor((c2 + 180.0) / bin_deg).astype(int)
    return list(zip(b1, b2))

cells = grid_cluster(chi1, chi2)
cell_count = defaultdict(int)
for c in cells:
    cell_count[c] += 1
top_cells = sorted(cell_count.items(), key=lambda kv: -kv[1])

# Transition events: count when consecutive frames belong to different cells.
n_trans = 0
for i in range(1, len(cells)):
    if cells[i] != cells[i - 1]:
        n_trans += 1

# ---------- Snapshot ↔ frame mapping ----------
# Snapshot filenames embed frame index: 7TL8_MTR6_snapNN_fFFF.pdb
import re
snap_dir = SYS_DIR / 'snapshots_n25'
snap_pdbs = sorted(snap_dir.glob('7TL8_MTR6_snap*_f*.pdb'))
print(f"[audit] {len(snap_pdbs)} snapshot PDBs found")

snap_to_frame = {}
fname_re = re.compile(r'_snap(\d+)_f(\d+)\.pdb$')
for pdb in snap_pdbs:
    m = fname_re.search(pdb.name)
    if m:
        snap_idx = int(m.group(1))
        frame_idx = int(m.group(2))
        snap_to_frame[pdb.name] = frame_idx

# ---------- Per-snap energy join ----------
energy_records = []
for pdb in snap_pdbs:
    name = pdb.stem  # e.g., 7TL8_MTR6_snap01_f9
    cands = list(PBSA_DIR.glob(f"{name}*.json"))
    e_val = None
    e_path = None
    for c in cands:
        try:
            with open(c) as f:
                d = json.load(f)
            for k in ('delta_g_kcal', 'delta_G', 'dG', 'mean_dg', 'binding_int_kcal',
                      'pbsa_dG_kcal', 'binding_kcal'):
                if k in d:
                    e_val = float(d[k])
                    e_path = str(c)
                    break
            if e_val is None and 'energy' in d:
                e_val = float(d['energy'])
                e_path = str(c)
            if e_val is not None:
                break
        except Exception:
            continue
    frame_idx = snap_to_frame.get(pdb.name)
    energy_records.append({
        'snap': name,
        'frame': int(frame_idx) if frame_idx is not None else None,
        'pbsa_dG_kcal': e_val,
        'json_path': e_path,
    })

# ---------- Write CSVs ----------
def write_csv(path, header, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)

# RMSD backbone binder
write_csv(EXT_DIR / 'rmsd_backbone_binder.csv',
          ['frame', 'time_ps', 'rmsd_bb_binder_A', 'binder_com_drift_A'],
          [(i, float(traj.time[i]) if i < len(traj.time) else i,
            float(rmsd_bb_binder[i]), float(binder_com_drift[i]))
           for i in range(traj.n_frames)])

# RMSD MTR side chain
write_csv(EXT_DIR / 'rmsd_mtr_sidechain.csv',
          ['frame', 'time_ps', 'rmsd_mtr_sc_A'],
          [(i, float(traj.time[i]) if i < len(traj.time) else i,
            float(rmsd_mtr_sc[i]))
           for i in range(traj.n_frames)])

# RMSF per-residue
write_csv(EXT_DIR / 'rmsf_per_residue_binder.csv',
          ['resSeq', 'resname', 'rmsf_A', 'n_heavy_atoms'],
          per_res_rmsf)

# RMSF MTR atoms
write_csv(EXT_DIR / 'rmsf_mtr_atoms.csv',
          ['atom_name', 'rmsf_A'],
          mtr_atom_rmsf)

# MTR chi1/chi2
write_csv(EXT_DIR / 'mtr_dihedrals.csv',
          ['frame', 'time_ps', 'chi1_deg', 'chi2_deg'],
          [(i, float(traj.time[i]) if i < len(traj.time) else i,
            float(chi1[i]), float(chi2[i]))
           for i in range(traj.n_frames)])

# Per-snap energy + frame
write_csv(EXT_DIR / 'per_snap_energy_pose.csv',
          ['snap', 'frame', 'pbsa_dG_kcal', 'rmsd_bb_binder_A',
           'rmsd_mtr_sc_A', 'binder_com_drift_A', 'chi1_deg', 'chi2_deg'],
          [(r['snap'], r['frame'], r['pbsa_dG_kcal'],
            float(rmsd_bb_binder[r['frame']]) if r['frame'] is not None else None,
            float(rmsd_mtr_sc[r['frame']]) if r['frame'] is not None else None,
            float(binder_com_drift[r['frame']]) if r['frame'] is not None else None,
            float(chi1[r['frame']]) if r['frame'] is not None else None,
            float(chi2[r['frame']]) if r['frame'] is not None else None)
           for r in energy_records])

# ---------- Aggregate summary ----------
def basic_stats(arr):
    a = np.asarray(arr, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {'n': 0}
    return {
        'n': int(a.size),
        'mean': float(a.mean()),
        'median': float(np.median(a)),
        'std': float(a.std(ddof=1)) if a.size > 1 else 0.0,
        'min': float(a.min()),
        'max': float(a.max()),
        'p95': float(np.percentile(a, 95.0)),
    }

summary = {
    'system': '7TL8_MTR6_calib_s7',
    'trajectory': str(DCD),
    'topology_pdb': str(TOP_PDB),
    'n_frames': int(traj.n_frames),
    'n_atoms': int(traj.n_atoms),
    'pbc_shifted_frames': int(n_shifted),
    'reference_frame': int(ref_frame),
    'binder_chain_idx': int(binder_chain_idx),
    'receptor_chain_idx': int(receptor_chain_idx),
    'mtr_residue': {
        'name': mtr_res.name,
        'chain_idx': int(binder_chain_idx),
        'resSeq': int(mtr_res.resSeq),
        'res_index': int(mtr_res.index),
        'atom_names': mtr_atom_names,
    },
    'rmsd_backbone_binder_A': basic_stats(rmsd_bb_binder),
    'rmsd_mtr_sidechain_A': basic_stats(rmsd_mtr_sc),
    'binder_com_drift_A': basic_stats(binder_com_drift),
    'rmsf_per_residue_binder_A': [
        {'resSeq': r[0], 'resname': r[1], 'rmsf_A': r[2], 'n_heavy': r[3]}
        for r in per_res_rmsf
    ],
    'rmsf_mtr_atoms_A': dict(mtr_atom_rmsf),
    'mtr_chi1_deg': basic_stats(chi1),
    'mtr_chi2_deg': basic_stats(chi2),
    'mtr_chi_grid_clusters': {
        'bin_deg': 60,
        'n_cells_occupied': len(cell_count),
        'top_cells': [{'cell': list(c), 'count': n} for c, n in top_cells[:10]],
        'n_transitions': int(n_trans),
        'transition_rate_per_frame': float(n_trans) / max(1, traj.n_frames - 1),
    },
    'per_snap_energies': energy_records,
}

# ---------- Mechanism classification ----------
# Heuristic thresholds, cited in the report:
#   Backbone RMSD median > 3 Å  -> (A) signal
#   MTR side-chain RMSF average > 1.5 Å  -> (B) signal
#   Binder COM drift median > 4 Å -> (C) signal
#   ChiX cluster transitions > 30 (over 500 frames, i.e. 6%/frame)  -> (B) signal

mtr_sc_avg_rmsf = float(np.mean([v for _, v in mtr_atom_rmsf]))

flags = {
    'A_backbone_instability': summary['rmsd_backbone_binder_A']['median'] > 3.0,
    'B_mtr_rotamer_flipping': (mtr_sc_avg_rmsf > 1.5) or (n_trans > 30),
    'C_binding_pose_drift': summary['binder_com_drift_A']['median'] > 4.0,
    'D_ff_artifact': False,  # populated below if PBSA energies look "spike-like"
}

# (D) check: does the σ_w come from a few high-energy outliers?
energies = [r['pbsa_dG_kcal'] for r in energy_records if r['pbsa_dG_kcal'] is not None]
if len(energies) > 5:
    e_arr = np.asarray(energies, dtype=float)
    e_med = float(np.median(e_arr))
    e_p95 = float(np.percentile(e_arr, 95.0))
    e_max = float(np.max(e_arr))
    n_outliers = int(((e_arr - e_med) > 50.0).sum())  # >50 kcal above median
    summary['energy_outlier_check'] = {
        'median': e_med, 'p95': e_p95, 'max': e_max,
        'n_outliers_gt50_above_median': n_outliers,
    }
    # FF artifact heuristic: ≥3 of 25 snapshots are >+50 kcal above median ⇒ outlier-driven σ_w
    flags['D_ff_artifact'] = n_outliers >= 3

summary['mechanism_flags'] = flags
summary['mtr_sc_avg_rmsf_A'] = mtr_sc_avg_rmsf

# Primary mechanism = strongest signal
# Combine into a single classification (A/B/C/D, ranked by evidence strength)
score = {
    'A': summary['rmsd_backbone_binder_A']['median'] / 3.0,
    'B': max(mtr_sc_avg_rmsf / 1.5, n_trans / 30.0),
    'C': summary['binder_com_drift_A']['median'] / 4.0,
    'D': (summary.get('energy_outlier_check', {}).get('n_outliers_gt50_above_median', 0)) / 3.0,
}
primary = max(score.items(), key=lambda kv: kv[1])
summary['mechanism_classification'] = {
    'scores_normalized': score,
    'primary': primary[0],
    'primary_score': primary[1],
}

def _to_native(o):
    """Convert numpy scalars/arrays to Python natives for JSON."""
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"Type {type(o).__name__} not serializable")

with open(EXT_DIR / 'audit_summary.json', 'w') as f:
    json.dump(summary, f, indent=2, default=_to_native)
print(f"[audit] summary written -> {EXT_DIR / 'audit_summary.json'}")

# ---------- Console summary ----------
print("\n=== AUDIT SUMMARY (pre-extension, on existing 500-frame trajectory) ===")
print(f"  binder backbone RMSD (Å):       median={summary['rmsd_backbone_binder_A']['median']:.2f}  "
      f"p95={summary['rmsd_backbone_binder_A']['p95']:.2f}  "
      f"max={summary['rmsd_backbone_binder_A']['max']:.2f}")
print(f"  MTR side-chain RMSD (Å):        median={summary['rmsd_mtr_sidechain_A']['median']:.2f}  "
      f"p95={summary['rmsd_mtr_sidechain_A']['p95']:.2f}  "
      f"max={summary['rmsd_mtr_sidechain_A']['max']:.2f}")
print(f"  binder COM drift (Å):           median={summary['binder_com_drift_A']['median']:.2f}  "
      f"max={summary['binder_com_drift_A']['max']:.2f}")
print(f"  MTR side-chain mean RMSF (Å):   {mtr_sc_avg_rmsf:.2f}")
print(f"  MTR χ-grid: {summary['mtr_chi_grid_clusters']['n_cells_occupied']} cells, "
      f"{n_trans} transitions, rate={summary['mtr_chi_grid_clusters']['transition_rate_per_frame']:.3f}/frame")
print(f"  Energy outliers (>median+50): {summary.get('energy_outlier_check', {}).get('n_outliers_gt50_above_median', 'NA')}")
print(f"  Flags: {summary['mechanism_flags']}")
print(f"  Primary mechanism: {summary['mechanism_classification']['primary']}  "
      f"(score={summary['mechanism_classification']['primary_score']:.2f})")
print("\nDone.")
