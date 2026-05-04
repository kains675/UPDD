#!/bin/bash
# Phase β Stage 1 — 7TL8_MTR6 + 2QKH_MTR25 Post-patch PBSA (multi-system v0.7.1 verification)
#
# Reduced scope (data inventory deviation):
#   7TL8_MTR6: only s7 has DCD (s19/s42 NVT-exploded, no trajectory)
#   2QKH_MTR25: all 3 DCDs exist but partial_md.pdb is empty stub (production MD
#               exploded at 0% before any topology PDB written) → unextractable
#   7TL8_WT: s7/s19/s42 all have snapshots_n25 ready
#   2QKH_WT: s7/s19/s42 all have snapshots_n25 ready
#
# Therefore this script runs:
#   - 7TL8_MTR6_s7      (post-patch ncAA, single-seed)
#   - 7TL8_WT_{s7,s19,s42} (control, 3 seeds for σ_btwn estimate)
#   - 2QKH_WT_{s7,s19,s42}  (control reference, 3 seeds)
#   (2QKH_MTR25 SKIPPED — pre-flight unextractable; documented in deliverable)
#
# Protocol matches Phase α Step 6/7 v2: Fix 4 (1-trajectory PBSA), 25 snapshots,
# patched MTR XMLs (Σq=0 verified), seeds {7,19,42}.
# Output dir: <run>/mmpbsa_results_fix4_postpatch/

set -u
PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

export UPDD_MMGBSA_PLATFORM=CUDA

LOG="/home/san/UPDD_proj/outputs/analysis/phase_beta_stage1_pbsa_20260428.log"
mkdir -p /home/san/UPDD_proj/outputs/analysis
echo "==========================================" | tee -a "$LOG"
echo "PHASE β STAGE 1 PBSA START: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "Inventory deviation: see scripts/phase_beta_stage1_pbsa.sh header" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"

run_pbsa () {
    local label="$1"
    local snap_dir="$2"
    local outdir="$3"
    local target_id="$4"
    local ncaa="$5"  # MTR or none
    if [ -f "${outdir}/mmpbsa_summary.json" ]; then
        echo "SKIP ${label} (already done)" | tee -a "$LOG"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        echo "MISS ${label} no pdb in ${snap_dir}" | tee -a "$LOG"
        return 1
    fi
    local start=$(date +%s)
    echo "[START ${label}] $(date --iso-8601=seconds) snap_dir=${snap_dir}" | tee -a "$LOG"
    UPDD_MMGBSA_PLATFORM=CUDA $PY scripts/run_mmpbsa.py \
      --md_dir "${snap_dir}" \
      --outputdir "${outdir}" \
      --ncaa_elem "${ncaa}" --receptor_chain A --binder_chain B \
      --target_id "${target_id}" \
      --protocol 1traj \
      >> "$LOG" 2>&1
    local rc=$?
    local end=$(date +%s)
    local elapsed=$((end-start))
    echo "[DONE ${label}] rc=${rc} elapsed=${elapsed}s" | tee -a "$LOG"
}

# 7TL8_MTR6 (ncAA, post-patch MTR XML) — only s7 viable
run_pbsa MTR6_s7   outputs/7TL8_MTR6_calib_s7/snapshots_n25  outputs/7TL8_MTR6_calib_s7/mmpbsa_results_fix4_postpatch  7TL8 MTR

# 2QKH_MTR25 — SKIPPED (all 3 unextractable)
echo "[SKIP 2QKH_MTR25_*] all DCDs unextractable (partial_md.pdb empty stub)" | tee -a "$LOG"

# 7TL8_WT (control, 3 seeds)
for s in s7 s19 s42; do
    run_pbsa WT7TL8_${s} outputs/7TL8_WT_calib_${s}/snapshots_n25 outputs/7TL8_WT_calib_${s}/mmpbsa_results_fix4_postpatch 7TL8 none
done

# 2QKH_WT (reference control, 3 seeds)
for s in s7 s19 s42; do
    run_pbsa WT2QKH_${s} outputs/2QKH_WT_calib_${s}/snapshots_n25 outputs/2QKH_WT_calib_${s}/mmpbsa_results_fix4_postpatch 2QKH none
done

echo "==========================================" | tee -a "$LOG"
echo "PHASE β STAGE 1 PBSA END: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
