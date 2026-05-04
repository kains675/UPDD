#!/bin/bash
# Track 1: 1EBP_MTR13 Fix 4 PBSA Multi-System Validation
# Cross-references Track B Cp4 PBSA (X.A UNCONDITIONAL) — does PBSA generalize across MTR scaffolds?
# Snapshots already extracted (n25 each for s7/s19/s42, snapshots/ for s23/s101); MD reuse.
set -u
PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

export UPDD_MMGBSA_PLATFORM=CUDA

LOG="/home/san/UPDD_proj/outputs/analysis/track_b_1ebp_mtr13_pbsa_20260427.log"
mkdir -p /home/san/UPDD_proj/outputs/analysis
echo "==========================================" | tee -a "$LOG"
echo "TRACK 1 1EBP_MTR13 PBSA START: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"

run_pbsa () {
    local label="$1"
    local snap_dir="$2"
    local outdir="$3"
    local target_id="$4"
    if [ -f "${outdir}/mmpbsa_summary.json" ]; then
        echo "SKIP ${label} (already done)" | tee -a "$LOG"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        echo "MISS ${label} no pdb in ${snap_dir}" | tee -a "$LOG"
        return 1
    fi
    local start=$(date +%s)
    echo "[START ${label}] $(date --iso-8601=seconds)" | tee -a "$LOG"
    UPDD_MMGBSA_PLATFORM=CUDA $PY scripts/run_mmpbsa.py \
      --md_dir "${snap_dir}" \
      --outputdir "${outdir}" \
      --ncaa_elem MTR --receptor_chain A --binder_chain B --target_id "${target_id}" \
      --protocol 1traj \
      >> "$LOG" 2>&1
    local rc=$?
    local end=$(date +%s)
    local elapsed=$((end-start))
    echo "[DONE ${label}] rc=${rc} elapsed=${elapsed}s" | tee -a "$LOG"
}

# ncAA replicates
run_pbsa MTR13_s7   outputs/1EBP_MTR13_calib_s7/snapshots_n25   outputs/1EBP_MTR13_calib_s7/mmpbsa_results_fix4   1EBP
run_pbsa MTR13_s19  outputs/1EBP_MTR13_calib_s19/snapshots_n25  outputs/1EBP_MTR13_calib_s19/mmpbsa_results_fix4  1EBP
run_pbsa MTR13_s42  outputs/1EBP_MTR13_calib_s42/snapshots_n25  outputs/1EBP_MTR13_calib_s42/mmpbsa_results_fix4  1EBP
run_pbsa MTR13_s23  outputs/1EBP_MTR13_calib_s23/snapshots      outputs/1EBP_MTR13_calib_s23/mmpbsa_results_fix4  1EBP
run_pbsa MTR13_s101 outputs/1EBP_MTR13_calib_s101/snapshots     outputs/1EBP_MTR13_calib_s101/mmpbsa_results_fix4 1EBP

# WT comparator (at minimum s42; if more WT seeds exist, run them too)
for s in s7 s19 s42 s23 s101; do
    snap="outputs/1EBP_WT_calib_${s}/snapshots_n25"
    [ -d "$snap" ] || snap="outputs/1EBP_WT_calib_${s}/snapshots"
    if ls "${snap}"/*.pdb >/dev/null 2>&1; then
        run_pbsa "WT_${s}" "${snap}" "outputs/1EBP_WT_calib_${s}/mmpbsa_results_fix4" 1EBP
    fi
done

echo "==========================================" | tee -a "$LOG"
echo "TRACK 1 1EBP_MTR13 PBSA END: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
