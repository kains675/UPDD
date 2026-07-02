#!/usr/bin/env python3
"""Binder–target heavy-atom contact / interface-engagement analysis (Paper 1 v3, §4.3.2).

Re-computes and PERSISTS the contact-engagement metric that §4.3.2 of the v3
manuscript asserts (patch-off ≈99% engaged, γ ≈20%, WT ≈100%; mean contacts
71±18 / 12±11 / 92±15). The original ADR-0009 full-ensemble computation was a
one-off whose per-snapshot output was never saved to disk; this script makes it
reproducible and versioned.

Metric (locked):
  - binder  = chain B heavy atoms (the 13-residue Cp4 / WT peptide, incl. MTR ncAA).
  - target  = chain A heavy atoms (the C3c receptor; res 1..631).
  - water (chain C, HOH) and ions (chain D, NA/CL) are EXCLUDED.
  - heavy   = PDB element column (77-78) not in {H, D}.
  - contact(snapshot) = number of (binder_heavy, target_heavy) atom pairs with
                        Euclidean distance < 4.0 Angstrom (scipy cKDTree.count_neighbors).
  - engaged(snapshot) = contact >= 20  (the ADR-0009 engagement threshold).

Snapshot set: the n=25 MM-PBSA production ensemble (snapshots_n25_postl387_patch_v2/),
i.e. the SAME snapshots that fed the postl387_v2 branched_ddg numbers.

Cohorts:
  - patch-off Cp4  : outputs/2QKI_Cp4_calib_s*        (Σq=-0.187 GAFF2-default; audit baseline)
  - gamma Cp4      : outputs/2QKI_Cp4_gamma_calib_s*  (Σq=0 Khoury full-RESP / γ production charges)
  - WT             : outputs/2QKI_WT_calib_s*

Output: outputs/analysis/contact_engagement_<stamp>/contact_engagement.json (+ .md summary).
Run with the qmmm env python (numpy + scipy):
  /home/san/miniconda3/envs/qmmm/bin/python scripts/contact_engagement_analysis.py
"""
from __future__ import annotations

import glob
import json
import os
import sys
from statistics import mean, pstdev

import numpy as np
from scipy.spatial import cKDTree

REPO = "/home/san/UPDD_proj"
SNAP_SUBDIR = "snapshots_n25_postl387_patch_v2"
CUTOFF = 4.0          # Angstrom
ENGAGED_MIN = 20      # contacts >= this => engaged
BINDER_CHAIN = "B"
TARGET_CHAIN = "A"

COHORTS = {
    "patchoff_Cp4": "2QKI_Cp4_calib_s*",
    "gamma_Cp4": "2QKI_Cp4_gamma_calib_s*",
    "WT": "2QKI_WT_calib_s*",
    # D3-refit MTR FF pose-recovery cohort (dt=1fs). The glob captures both the
    # fresh d5pose seeds (2QKI_Cp4_d3refit_d5pose_s*_dt1fs) and the reused dt=1fs
    # control (2QKI_Cp4_d3refit_localize_s211_dt1fs) for n>=4, while excluding the
    # dt=2fs localize dirs (..._dt2fs) which crashed and are not band-comparable.
    "d3refit_Cp4": "2QKI_Cp4_d3refit_*_dt1fs",
    # WT positive control (D5 NO-GO specificity, P3 #107) — strongest binder run
    # through the SAME dt=1fs / 25ns protocol + metric. WT engaged => the protocol
    # can report engagement => the d3refit NO-GO is charge-specific; WT disengaged
    # => dt=1fs/pipeline suspect. NOTE: WT ran --ncaa none = UNRESTRAINED (anchor
    # asymmetry vs the res-4-anchored d3refit cohort; see the launcher header).
    "WT_dt1fs": "2QKI_WT_d5posctrl_s*",
    # patch-off MTR PRIMARY positive control (D5 NO-GO specificity, P3 #107) — the
    # same 2QKI_Cp4 (MTR res-4) system run through the SAME dt=1fs / 25ns protocol +
    # metric, with the patch-off MTR FF (calib_s101 MTR_gaff2.xml: NE1=-0.3418 frozen,
    # CM=0.0487, Σq=0) installed instead of the d3refit FF. Anchor-SYMMETRIC vs the
    # d3refit cohort (both MTR res-4 backbone anchored), so the ONLY variable is the
    # sidechain RESP-A2 refit charge -> apples-to-apples single-variable isolation.
    "patchoff_dt1fs": "2QKI_Cp4_patchoff_d5posctrl_s*",
}


def parse_pdb_heavy(path: str):
    """Return (binder_xyz Nx3, target_xyz Mx3) heavy-atom coords for chains B / A."""
    binder, target = [], []
    with open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            elem = line[76:78].strip()
            if elem in ("H", "D"):
                continue
            ch = line[21]
            if ch == BINDER_CHAIN:
                dst = binder
            elif ch == TARGET_CHAIN:
                dst = target
            else:
                continue  # water (C) / ions (D) excluded
            try:
                dst.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    return np.asarray(binder, float), np.asarray(target, float)


def snapshot_contacts(path: str) -> int:
    b, t = parse_pdb_heavy(path)
    if b.size == 0 or t.size == 0:
        return -1  # malformed -> flag
    # count of (binder, target) atom pairs within CUTOFF
    return int(cKDTree(b).count_neighbors(cKDTree(t), CUTOFF))


def seed_token(seed_dir: str, cohort_glob: str) -> str:
    base = os.path.basename(seed_dir)
    return base.split("_calib_")[-1]


def main() -> int:
    stamp = sys.argv[1] if len(sys.argv) > 1 else "20260612"
    out_dir = os.path.join(REPO, "outputs", "analysis", f"contact_engagement_{stamp}")
    os.makedirs(out_dir, exist_ok=True)

    report = {
        "schema": "contact_engagement/0.1",
        "stamp": stamp,
        "method": {
            "binder": "chain B heavy atoms (13-mer peptide incl. MTR)",
            "target": "chain A heavy atoms (C3c receptor, res 1..631); water(C)/ions(D) excluded",
            "heavy": "PDB element col not in {H,D}",
            "contact": "count of (binder_heavy, target_heavy) atom pairs < 4.0 Angstrom (cKDTree.count_neighbors)",
            "engaged_threshold": ENGAGED_MIN,
            "snapshot_set": SNAP_SUBDIR,
        },
        "cohorts": {},
    }

    for cohort, pat in COHORTS.items():
        seed_dirs = sorted(glob.glob(os.path.join(REPO, "outputs", pat)))
        seeds_out = {}
        pooled = []
        for sd in seed_dirs:
            snaps = sorted(glob.glob(os.path.join(sd, SNAP_SUBDIR, "*.pdb")))
            if not snaps:
                continue
            counts = []
            for sp in snaps:
                c = snapshot_contacts(sp)
                if c < 0:
                    print(f"  [WARN] malformed/empty selection: {sp}", file=sys.stderr)
                    continue
                counts.append(c)
            if not counts:
                continue
            tok = seed_token(sd, pat)
            n_eng = sum(1 for c in counts if c >= ENGAGED_MIN)
            seeds_out[tok] = {
                "n_snap": len(counts),
                "mean_contacts": round(mean(counts), 2),
                "sd_contacts": round(pstdev(counts), 2) if len(counts) > 1 else 0.0,
                "n_engaged": n_eng,
                "frac_engaged": round(n_eng / len(counts), 4),
                "per_snapshot_contacts": counts,
            }
            pooled.extend(counts)
        n_eng_pool = sum(1 for c in pooled if c >= ENGAGED_MIN)
        report["cohorts"][cohort] = {
            "n_seeds": len(seeds_out),
            "n_snapshots": len(pooled),
            "pooled_mean_contacts": round(mean(pooled), 2) if pooled else None,
            "pooled_sd_contacts": round(pstdev(pooled), 2) if len(pooled) > 1 else None,
            "pooled_frac_engaged": round(n_eng_pool / len(pooled), 4) if pooled else None,
            "pooled_pct_engaged": round(100 * n_eng_pool / len(pooled), 1) if pooled else None,
            "seeds": seeds_out,
        }

    jpath = os.path.join(out_dir, "contact_engagement.json")
    with open(jpath, "w") as fh:
        json.dump(report, fh, indent=2)

    # readable summary + comparison to the manuscript / ADR-0009 assertions
    lines = ["# Binder-target contact engagement (Paper 1 v3 §4.3.2 re-computation)\n",
             f"stamp={stamp}  cutoff={CUTOFF} A  engaged>= {ENGAGED_MIN} contacts  set={SNAP_SUBDIR}\n",
             "| cohort | seeds | snaps | mean±SD contacts | %engaged |",
             "|---|---|---|---|---|"]
    claim = {"patchoff_Cp4": "71±18 / 99%", "gamma_Cp4": "12±11 / 20%", "WT": "92±15 / 100%"}
    for cohort, d in report["cohorts"].items():
        lines.append(f"| {cohort} | {d['n_seeds']} | {d['n_snapshots']} | "
                     f"{d['pooled_mean_contacts']}±{d['pooled_sd_contacts']} | {d['pooled_pct_engaged']}% |  "
                     f"(ADR-0009 claim: {claim.get(cohort,'?')})")
    summary = "\n".join(lines) + "\n"
    with open(os.path.join(out_dir, "contact_engagement_summary.md"), "w") as fh:
        fh.write(summary)

    print(summary)
    print(f"[written] {jpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
