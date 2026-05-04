#!/bin/bash
# Phase α Step 6: Cp4 Fix 4 PBSA Post-Patch Re-run
# Same protocol as Track B Cp4 Phase 2; only difference is patched MTR XMLs (Σq=0, Δq_N=0).
# Output to *_postpatch/ to preserve Track B baseline (X.A UNCONDITIONAL: ΔΔG=+3.92, σ_w=3.62).
set -u
PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

export UPDD_MMGBSA_PLATFORM=CUDA
# Phase α default-ON behavior is in effect (no env override needed)

LOG="/home/san/UPDD_proj/outputs/analysis/track_b_cp4_postpatch_pbsa_20260427.log"
mkdir -p /home/san/UPDD_proj/outputs/analysis
echo "==========================================" | tee -a "$LOG"
echo "STEP 6 Cp4 POST-PATCH PBSA START: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"

run_pbsa () {
    local label="$1"
    local snap_dir="$2"
    local outdir="$3"
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
      --ncaa_elem MTR --receptor_chain A --binder_chain B --target_id 2QKI \
      --protocol 1traj \
      >> "$LOG" 2>&1
    local rc=$?
    local end=$(date +%s)
    local elapsed=$((end-start))
    echo "[DONE ${label}] rc=${rc} elapsed=${elapsed}s" | tee -a "$LOG"
}

run_pbsa Cp4_s7   outputs/2QKI_Cp4_calib_s7/snapshots_n25         outputs/2QKI_Cp4_calib_s7/mmpbsa_results_fix4_postpatch
run_pbsa Cp4_s42  outputs/2QKI_Cp4_calib_s42/snapshots_n25        outputs/2QKI_Cp4_calib_s42/mmpbsa_results_fix4_postpatch
run_pbsa Cp4_s23  outputs/2QKI_Cp4_calib_s23/snapshots            outputs/2QKI_Cp4_calib_s23/mmpbsa_results_fix4_postpatch
run_pbsa Cp4_s101 outputs/2QKI_Cp4_calib_s101/snapshots           outputs/2QKI_Cp4_calib_s101/mmpbsa_results_fix4_postpatch
# s19 used reseed55 in Track B Phase 2 due to MD failure on original s19; reuse same seed for parity
if [ -d outputs/2QKI_Cp4_calib_s19_reseed55/snapshots ]; then
  run_pbsa Cp4_s19_reseed55 outputs/2QKI_Cp4_calib_s19_reseed55/snapshots outputs/2QKI_Cp4_calib_s19_reseed55/mmpbsa_results_fix4_postpatch
else
  run_pbsa Cp4_s19 outputs/2QKI_Cp4_calib_s19/snapshots_n25 outputs/2QKI_Cp4_calib_s19/mmpbsa_results_fix4_postpatch
fi

echo "==========================================" | tee -a "$LOG"
echo "STEP 6 Cp4 POST-PATCH PBSA END: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
