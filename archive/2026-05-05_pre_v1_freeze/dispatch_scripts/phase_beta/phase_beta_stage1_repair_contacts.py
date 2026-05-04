#!/usr/bin/env python
"""
MTR-pocket contact analysis: identify which receptor residues sit within 4 Å of MTR
heavy atoms across the trajectory, and how that contact pattern correlates with
PBSA σ_w outliers.
"""
import os
import json
import csv
from pathlib import Path
import numpy as np
import mdtraj as md

ROOT = Path('/home/san/UPDD_proj')
SYS_DIR = ROOT / 'outputs' / '7TL8_MTR6_calib_s7'
DCD = SYS_DIR / 'mdresult' / '7TL8_MTR6_restrained.dcd'
TOP_PDB = SYS_DIR / 'mdresult' / '7TL8_MTR6_final.pdb'
EXT_DIR = SYS_DIR / 'mdresult' / 'extension_audit_20260428'

print("[contacts] loading trajectory ...")
traj = md.load_dcd(str(DCD), top=str(TOP_PDB))
top = traj.topology

# image_molecules + canonical bound ref
mtr_res = next(r for r in top.residues if r.name == 'MTR')
binder_chain_idx = mtr_res.chain.index
PROT_AA = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP MTR'.split())
recv_chain_idx = next(ch.index for ch in top.chains if ch.index != binder_chain_idx
                     and len([r for r in ch.residues if r.name in PROT_AA])>50)
recv_atoms = set(top.chain(recv_chain_idx).atoms)
binder_atoms_set = set(top.chain(binder_chain_idx).atoms)

print("[contacts] image_molecules ...")
traj = traj.image_molecules(inplace=False, anchor_molecules=[recv_atoms],
                             other_molecules=[binder_atoms_set], make_whole=True)

# Frames bounded mask
recv_bb = top.select(f'chainid {recv_chain_idx} and (name CA or name C or name N or name O)')
binder_heavy = top.select(f'chainid {binder_chain_idx} and not symbol H')
mtr_heavy = np.array([a.index for a in mtr_res.atoms if a.element.symbol != 'H'], dtype=int)

# Per-frame binder COM offset
all_bind_com = np.array([traj.xyz[i, binder_heavy, :].mean(0) for i in range(traj.n_frames)])
all_recv_com = np.array([traj.xyz[i, recv_bb, :].mean(0) for i in range(traj.n_frames)])
delta = all_bind_com - all_recv_com
median_delta = np.median(delta, axis=0)
distances = np.linalg.norm(delta - median_delta, axis=1)
ref_frame = int(np.argmin(distances))
traj.superpose(traj, frame=ref_frame, atom_indices=recv_bb)
print(f"[contacts] ref_frame={ref_frame}")

# Recompute COM after alignment, identify bound vs drift
all_bind_com = np.array([traj.xyz[i, binder_heavy, :].mean(0) for i in range(traj.n_frames)])
ref_com = traj.xyz[ref_frame, binder_heavy, :].mean(0)
binder_com_drift = np.linalg.norm(all_bind_com - ref_com, axis=1) * 10.0  # Å
bound_mask = binder_com_drift < 8.0
print(f"[contacts] bound frames: {bound_mask.sum()}/{traj.n_frames}")

# Receptor heavy atoms (per residue)
recv_residues = [r for r in top.chain(recv_chain_idx).residues if r.name in PROT_AA]

# Compute per-frame minimum distance from MTR heavy atoms to each receptor residue's heavy atoms
# Use only bound frames for this analysis.
print("[contacts] computing min-distance MTR <-> each recv residue (bound frames) ...")
bound_frames = np.where(bound_mask)[0]
contact_freq = np.zeros(len(recv_residues), dtype=float)
mean_min_dist = np.zeros(len(recv_residues), dtype=float)
THRESH_NM = 0.4  # 4 Å
for ridx, r in enumerate(recv_residues):
    ratoms = np.array([a.index for a in r.atoms if a.element.symbol != 'H'], dtype=int)
    if ratoms.size == 0:
        continue
    # Per-frame min distance MTR<->this residue
    mins = np.zeros(len(bound_frames), dtype=float)
    for k, i in enumerate(bound_frames):
        # broadcast difference
        d = traj.xyz[i, mtr_heavy[:, None], :] - traj.xyz[i, ratoms[None, :], :]
        # d shape: (n_mtr, n_r, 3)
        mins[k] = np.sqrt((d * d).sum(axis=2)).min()
    contact_freq[ridx] = float((mins < THRESH_NM).mean())
    mean_min_dist[ridx] = float(mins.mean())

# Top 15 contacts
order = np.argsort(-contact_freq)[:15]
print("\nTop 15 receptor residues in contact with MTR (bound subset)")
print(f"{'rank':<5}{'resSeq':<8}{'resname':<8}{'freq':<10}{'mean_min_d_A'}")
for rank, ridx in enumerate(order, start=1):
    r = recv_residues[ridx]
    print(f"{rank:<5}{r.resSeq:<8}{r.name:<8}{contact_freq[ridx]:<10.2f}{mean_min_dist[ridx]*10:.2f}")

# Save
out = []
for ridx, r in enumerate(recv_residues):
    if contact_freq[ridx] > 0.05:
        out.append({
            'resSeq': int(r.resSeq), 'resname': r.name,
            'contact_freq_bound': float(contact_freq[ridx]),
            'mean_min_dist_A_bound': float(mean_min_dist[ridx] * 10),
        })
with open(EXT_DIR / 'mtr_pocket_contacts.json', 'w') as f:
    json.dump({
        'threshold_nm': THRESH_NM,
        'n_bound_frames': int(bound_mask.sum()),
        'ref_frame': int(ref_frame),
        'contacts_bound': out,
    }, f, indent=2)
print(f"\n=> Saved -> {EXT_DIR/'mtr_pocket_contacts.json'}")
