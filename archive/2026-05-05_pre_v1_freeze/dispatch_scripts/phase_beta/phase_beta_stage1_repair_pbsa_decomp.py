#!/usr/bin/env python
"""
PBSA energy decomposition + pocket-contact analysis for 7TL8_MTR6_s7's σ_w=38.33.

Reads the per-snap mmpbsa JSONs and decomposes ΔG into PB + nonpolar + dispersion
contributions, identifies which term drives the high-σ outliers (snap04 +98,
snap05 +92, snap11 +78, snap21 +121, etc.), and reports MTR-pocket contact
distances on the bound-only subset.
"""
import os
import json
import csv
import math
from pathlib import Path
import numpy as np
import mdtraj as md

ROOT = Path('/home/san/UPDD_proj')
SYS_DIR = ROOT / 'outputs' / '7TL8_MTR6_calib_s7'
PBSA_DIR = SYS_DIR / 'mmpbsa_results_fix4_postpatch'
EXT_DIR = SYS_DIR / 'mdresult' / 'extension_audit_20260428'
EXT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- Read PBSA snap JSONs ----------
records = []
for j in sorted(PBSA_DIR.glob('7TL8_MTR6_snap*_f*.json')):
    if j.name.endswith('mmpbsa_summary.json'):
        continue
    d = json.load(open(j))
    rec = {
        'snap': d.get('snapshot'),
        'pdb_path': d.get('pdb_path'),
        'frame': int(j.stem.split('_f')[1].split('_')[0]),
        'dG': d.get('delta_g_kcal'),
        'pb': d.get('delta_epb_kcal'),
        'np': d.get('delta_enpolar_kcal'),
        'disp': d.get('delta_edisper_kcal'),
        'fav': d.get('favorable'),
    }
    records.append(rec)

# ---------- Load audit CSV (frame -> bound/drift class) ----------
pose = {}
with open(EXT_DIR / 'per_snap_energy_pose.csv') as f:
    for r in csv.DictReader(f):
        pose[r['snap']] = r

# ---------- Decompose σ ----------
def stats(arr):
    a = np.asarray(arr, dtype=float)
    return {
        'mean': float(a.mean()),
        'median': float(np.median(a)),
        'std': float(a.std(ddof=1)) if a.size > 1 else 0.0,
        'min': float(a.min()),
        'max': float(a.max()),
    }

# Partition by bound vs drift (using audit cutoffs)
def bound_or_drift(snap_name):
    p = pose.get(snap_name)
    if not p:
        return 'unknown'
    com = float(p['binder_com_drift_A'])
    bb = float(p['rmsd_bb_binder_A'])
    return 'BOUND' if (com < 8.0 and bb < 15.0) else 'DRIFT'

bound = [r for r in records if bound_or_drift(r['snap']) == 'BOUND']
drift = [r for r in records if bound_or_drift(r['snap']) == 'DRIFT']

print(f"BOUND: n={len(bound)}, DRIFT: n={len(drift)}")
for label, group in [('ALL', records), ('BOUND', bound), ('DRIFT', drift)]:
    if not group:
        continue
    print(f"\n=== {label} (n={len(group)}) ===")
    print(f"  dG    : {stats([r['dG'] for r in group])}")
    print(f"  PB    : {stats([r['pb'] for r in group])}")
    print(f"  NP    : {stats([r['np'] for r in group])}")
    print(f"  Disp  : {stats([r['disp'] for r in group])}")

# ---------- Identify which term carries σ_w ----------
all_dG = np.array([r['dG'] for r in records])
all_pb = np.array([r['pb'] for r in records])
all_np = np.array([r['np'] for r in records])
all_disp = np.array([r['disp'] for r in records])

print("\n=== σ component decomposition (ALL n=25) ===")
print(f"  σ(dG)   = {all_dG.std(ddof=1):.2f} kcal")
print(f"  σ(PB)   = {all_pb.std(ddof=1):.2f} kcal")
print(f"  σ(NP)   = {all_np.std(ddof=1):.2f} kcal")
print(f"  σ(Disp) = {all_disp.std(ddof=1):.2f} kcal")

bound_dG = np.array([r['dG'] for r in bound])
bound_pb = np.array([r['pb'] for r in bound])
bound_np = np.array([r['np'] for r in bound])
bound_disp = np.array([r['disp'] for r in bound])
print(f"\n=== σ component decomposition (BOUND n={len(bound)}) ===")
print(f"  σ(dG)   = {bound_dG.std(ddof=1):.2f} kcal")
print(f"  σ(PB)   = {bound_pb.std(ddof=1):.2f} kcal")
print(f"  σ(NP)   = {bound_np.std(ddof=1):.2f} kcal")
print(f"  σ(Disp) = {bound_disp.std(ddof=1):.2f} kcal")

# ---------- Output JSON ----------
out = {
    'all': {
        'n': len(records),
        'sigma_dG': float(all_dG.std(ddof=1)),
        'sigma_PB': float(all_pb.std(ddof=1)),
        'sigma_NP': float(all_np.std(ddof=1)),
        'sigma_Disp': float(all_disp.std(ddof=1)),
        'mean_dG': float(all_dG.mean()),
    },
    'bound': {
        'n': len(bound),
        'sigma_dG': float(bound_dG.std(ddof=1)) if len(bound) > 1 else 0.0,
        'sigma_PB': float(bound_pb.std(ddof=1)) if len(bound) > 1 else 0.0,
        'sigma_NP': float(bound_np.std(ddof=1)) if len(bound) > 1 else 0.0,
        'sigma_Disp': float(bound_disp.std(ddof=1)) if len(bound) > 1 else 0.0,
        'mean_dG': float(bound_dG.mean()) if len(bound_dG) else None,
    },
    'drift': {
        'n': len(drift),
        'sigma_dG': float(np.std([r['dG'] for r in drift], ddof=1)) if len(drift) > 1 else 0.0,
        'mean_dG': float(np.mean([r['dG'] for r in drift])) if drift else None,
    },
    'top_outliers_by_dG': sorted(records, key=lambda r: -r['dG'])[:5],
    'classification': {},
}

# Classification heuristic:
# If σ(PB) ≫ σ(NP) and σ(Disp) ⇒ electrostatic-pocket flicker (typical for charged-Σq net 0
#   but locally dipole-laden MTR + Glu/Asp neighbors)
# If σ(PB) ~ σ(Disp) with both > 30 ⇒ dispersion + electrostatic co-fluctuation, likely a
#   pocket-volume oscillation
# In Phase α §3 MTR analyses on 1EBP, σ_w ~ 4 with σ(PB) ~ 3 (dominant). Here σ(PB) is
#   what we should sanity-check.

out['classification']['sigma_PB_dominant'] = (
    out['bound']['sigma_PB'] > 0.7 * out['bound']['sigma_dG']
)
out['classification']['sigma_NP_dominant'] = (
    out['bound']['sigma_NP'] > 0.7 * out['bound']['sigma_dG']
)
out['classification']['sigma_Disp_dominant'] = (
    out['bound']['sigma_Disp'] > 0.7 * out['bound']['sigma_dG']
)

with open(EXT_DIR / 'pbsa_decomposition.json', 'w') as f:
    json.dump(out, f, indent=2, default=lambda o: float(o) if hasattr(o, 'item') else str(o))
print(f"\n=> Wrote {EXT_DIR/'pbsa_decomposition.json'}")

# ---------- Top outliers ----------
print("\n=== Top 5 |dG| snaps ===")
for r in sorted(records, key=lambda r: -r['dG'])[:5]:
    print(f"  {r['snap']:<25} dG={r['dG']:>+8.2f}  PB={r['pb']:>+9.2f}  "
          f"NP={r['np']:>+8.2f}  Disp={r['disp']:>+8.2f}  [{bound_or_drift(r['snap'])}]")
