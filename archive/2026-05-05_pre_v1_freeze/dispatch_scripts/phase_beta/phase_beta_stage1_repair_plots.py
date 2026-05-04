#!/usr/bin/env python
"""Phase β Stage 1-Repair audit plots."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

EXT_DIR = Path('/home/san/UPDD_proj/outputs/7TL8_MTR6_calib_s7/mdresult/extension_audit_20260428')

# 1. Per-frame backbone RMSD + binder COM drift
times = []
bb = []
com = []
with open(EXT_DIR / 'rmsd_backbone_binder.csv') as f:
    for r in csv.DictReader(f):
        times.append(float(r['time_ps']))
        bb.append(float(r['rmsd_bb_binder_A']))
        com.append(float(r['binder_com_drift_A']))
times = np.array(times)
bb = np.array(bb)
com = np.array(com)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
ax1.plot(times, bb, lw=0.7, color='C0')
ax1.set_ylabel('Binder backbone RMSD (Å)')
ax1.axhline(15, color='red', ls='--', lw=0.5, label='15 Å bound/drift cutoff')
ax1.set_ylim(0, 80)
ax1.legend()
ax1.set_title('7TL8_MTR6_s7 — backbone RMSD + binder COM drift\n'
              '(jumps to ~70 Å / ~40 Å signal periodic-image hopping at VAL5/MTR6 junction)')

ax2.plot(times, com, lw=0.7, color='C1')
ax2.set_xlabel('Time (ps)')
ax2.set_ylabel('Binder COM drift (Å)')
ax2.axhline(8, color='red', ls='--', lw=0.5, label='8 Å bound/drift cutoff')
ax2.set_ylim(0, 50)
ax2.legend()
fig.tight_layout()
fig.savefig(EXT_DIR / 'plot_rmsd_com.png', dpi=120)
print(f"saved {EXT_DIR/'plot_rmsd_com.png'}")
plt.close()

# 2. MTR side-chain RMSD + chi1/chi2
times2 = []
mtr_rmsd = []
with open(EXT_DIR / 'rmsd_mtr_sidechain.csv') as f:
    for r in csv.DictReader(f):
        times2.append(float(r['time_ps']))
        mtr_rmsd.append(float(r['rmsd_mtr_sc_A']))
times2 = np.array(times2)
mtr_rmsd = np.array(mtr_rmsd)

times3 = []
chi1 = []
chi2 = []
with open(EXT_DIR / 'mtr_dihedrals.csv') as f:
    for r in csv.DictReader(f):
        times3.append(float(r['time_ps']))
        chi1.append(float(r['chi1_deg']))
        chi2.append(float(r['chi2_deg']))
times3 = np.array(times3)
chi1 = np.array(chi1)
chi2 = np.array(chi2)

fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
axes[0].plot(times2, mtr_rmsd, lw=0.7, color='C2')
axes[0].set_ylabel('MTR side-chain RMSD (Å)')
axes[0].set_title('7TL8_MTR6_s7 — MTR side-chain RMSD + χ₁ + χ₂\n'
                  '(stable single rotamer well: χ₁ ≈ −68° ± 17°, χ₂ ≈ +102° ± 10°)')
axes[1].plot(times3, chi1, lw=0.7, color='C3')
axes[1].set_ylabel('χ₁ (deg)')
axes[1].axhline(-60, color='gray', ls=':', lw=0.4)
axes[1].axhline(-90, color='gray', ls=':', lw=0.4)
axes[1].set_ylim(-180, 180)
axes[2].plot(times3, chi2, lw=0.7, color='C4')
axes[2].set_ylabel('χ₂ (deg)')
axes[2].set_xlabel('Time (ps)')
axes[2].axhline(60, color='gray', ls=':', lw=0.4)
axes[2].axhline(120, color='gray', ls=':', lw=0.4)
axes[2].set_ylim(-180, 180)
fig.tight_layout()
fig.savefig(EXT_DIR / 'plot_mtr_dynamics.png', dpi=120)
print(f"saved {EXT_DIR/'plot_mtr_dynamics.png'}")
plt.close()

# 3. χ₁ vs χ₂ scatter
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(chi1, chi2, s=8, alpha=0.5, c=times3, cmap='viridis')
ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
ax.set_xlabel('χ₁ (deg)')
ax.set_ylabel('χ₂ (deg)')
ax.set_title('7TL8_MTR6_s7 MTR — χ₁/χ₂ scatter (color = time)')
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(EXT_DIR / 'plot_chi1_chi2_scatter.png', dpi=120)
print(f"saved {EXT_DIR/'plot_chi1_chi2_scatter.png'}")
plt.close()

# 4. Per-residue RMSF binder
res_seq = []
res_name = []
rmsf = []
with open(EXT_DIR / 'rmsf_per_residue_binder.csv') as f:
    for r in csv.DictReader(f):
        res_seq.append(int(r['resSeq']))
        res_name.append(r['resname'])
        rmsf.append(float(r['rmsf_A']))

fig, ax = plt.subplots(figsize=(8, 4))
bars = ax.bar(res_seq, rmsf,
              color=['darkred' if v > 5 else 'C0' for v in rmsf])
ax.set_xticks(res_seq)
ax.set_xticklabels([f"{n}\n{s}" for s, n in zip(res_seq, res_name)], fontsize=8)
ax.set_ylabel('Heavy-atom RMSF (Å)')
ax.set_title('7TL8_MTR6_s7 — per-residue RMSF (binder chain)\n'
             'Residues 1–5 RMSF ≈ 40 Å = box/2 (PBC image-hopping signature)')
ax.set_yscale('symlog', linthresh=2)
fig.tight_layout()
fig.savefig(EXT_DIR / 'plot_per_residue_rmsf.png', dpi=120)
print(f"saved {EXT_DIR/'plot_per_residue_rmsf.png'}")
plt.close()

# 5. Per-snap dG vs imaging state
snaps = []
dgs = []
bbs = []
coms = []
with open(EXT_DIR / 'per_snap_energy_pose.csv') as f:
    for r in csv.DictReader(f):
        if r['pbsa_dG_kcal'] in ('', 'None'):
            continue
        snaps.append(int(r['snap'].split('snap')[1].split('_')[0]))
        dgs.append(float(r['pbsa_dG_kcal']))
        bbs.append(float(r['rmsd_bb_binder_A']))
        coms.append(float(r['binder_com_drift_A']))

snaps = np.array(snaps); dgs = np.array(dgs); bbs = np.array(bbs); coms = np.array(coms)
bound = (coms < 8) & (bbs < 15)

fig, ax = plt.subplots(figsize=(9, 4.5))
ax.bar(snaps[bound], dgs[bound], color='C0', label='bound (n=18)')
ax.bar(snaps[~bound], dgs[~bound], color='C3', label='periodic-image artifact (n=7)')
ax.axhline(0, color='black', lw=0.5)
ax.set_xlabel('Snapshot index')
ax.set_ylabel('PBSA ΔG (kcal/mol)')
ax.set_title('7TL8_MTR6_s7 — per-snap ΔG_PBSA, bound vs periodic-image classification\n'
             'σ_w(ALL)=38.33 kcal; σ(PB,bound)=90.61 kcal — driven by HETATM/ATOM PBC split, NOT real fluctuation')
ax.legend()
fig.tight_layout()
fig.savefig(EXT_DIR / 'plot_per_snap_dg.png', dpi=120)
print(f"saved {EXT_DIR/'plot_per_snap_dg.png'}")
plt.close()

print('Done.')
