#!/bin/bash
# Phase α Step 7: 1EBP_MTR13 PBSA Post-Patch Re-run + 1EBP_WT control re-run
# Same protocol as Track 1; only difference is patched MTR XMLs (WT unaffected, control sanity).
# Output to *_postpatch/ to preserve Track 1 baseline (X.C: ΔΔG=+2.91, z_SE=1.88).
set -u
PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

export UPDD_MMGBSA_PLATFORM=CUDA

LOG="/home/san/UPDD_proj/outputs/analysis/track_b_1ebp_mtr13_postpatch_pbsa_20260427.log"
mkdir -p /home/san/UPDD_proj/outputs/analysis
echo "==========================================" | tee -a "$LOG"
echo "STEP 7 1EBP_MTR13 POST-PATCH PBSA START: $(date --iso-8601=seconds)" | tee -a "$LOG"
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

# 1EBP_MTR13 (5 seeds, post-patch MTR XMLs)
run_pbsa MTR13_s7   outputs/1EBP_MTR13_calib_s7/snapshots_n25   outputs/1EBP_MTR13_calib_s7/mmpbsa_results_fix4_postpatch   1EBP
run_pbsa MTR13_s19  outputs/1EBP_MTR13_calib_s19/snapshots_n25  outputs/1EBP_MTR13_calib_s19/mmpbsa_results_fix4_postpatch  1EBP
run_pbsa MTR13_s42  outputs/1EBP_MTR13_calib_s42/snapshots_n25  outputs/1EBP_MTR13_calib_s42/mmpbsa_results_fix4_postpatch  1EBP
run_pbsa MTR13_s23  outputs/1EBP_MTR13_calib_s23/snapshots      outputs/1EBP_MTR13_calib_s23/mmpbsa_results_fix4_postpatch  1EBP
run_pbsa MTR13_s101 outputs/1EBP_MTR13_calib_s101/snapshots     outputs/1EBP_MTR13_calib_s101/mmpbsa_results_fix4_postpatch 1EBP

# 1EBP_WT (5 seeds, control — no MTR, parameterization unchanged; expected near-identical to Track 1 WT)
for s in s7 s19 s42 s23 s101; do
    snap="outputs/1EBP_WT_calib_${s}/snapshots_n25"
    [ -d "$snap" ] || snap="outputs/1EBP_WT_calib_${s}/snapshots"
    if ls "${snap}"/*.pdb >/dev/null 2>&1; then
        run_pbsa "WT_${s}" "${snap}" "outputs/1EBP_WT_calib_${s}/mmpbsa_results_fix4_postpatch" 1EBP
    fi
done

echo "==========================================" | tee -a "$LOG"
echo "STEP 7 1EBP_MTR13 POST-PATCH PBSA END: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
