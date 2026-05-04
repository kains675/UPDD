#!/bin/bash
# Track B Phase 2: Cp4 PBSA with fallback fix (5 seeds x 25 snap)
set -u
PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

export UPDD_MMGBSA_PLATFORM=CUDA

echo "=========================================="
echo "TRACK B PHASE 2 Cp4 PBSA BATCH START: $(date --iso-8601=seconds)"
echo "=========================================="

run_cp4_mmpbsa() {
    local seed_label="$1"   # 's7', 's42', 's23', 's101', 's19_reseed55'
    local snap_dir="$2"
    local outdir="$3"
    [ -f "${outdir}/mmpbsa_summary.json" ] && { echo "SKIP ${seed_label} (already done)"; return 0; }
    local start=$(date +%s)
    echo "[START Cp4_${seed_label}] $(date --iso-8601=seconds)"
    UPDD_MMGBSA_PLATFORM=CUDA $PY scripts/run_mmpbsa.py \
      --md_dir "${snap_dir}" \
      --outputdir "${outdir}" \
      --ncaa_elem MTR --receptor_chain A --binder_chain B --target_id 2QKI \
      --protocol 1traj \
      > "/tmp/track_b_cp4_phase2_${seed_label}.log" 2>&1
    local rc=$?
    local end=$(date +%s)
    local elapsed=$((end-start))
    echo "[DONE Cp4_${seed_label}] rc=${rc} elapsed=${elapsed}s"
}

run_cp4_mmpbsa s7            outputs/2QKI_Cp4_calib_s7/snapshots_n25         outputs/2QKI_Cp4_calib_s7/mmpbsa_results_fix4
run_cp4_mmpbsa s42           outputs/2QKI_Cp4_calib_s42/snapshots_n25        outputs/2QKI_Cp4_calib_s42/mmpbsa_results_fix4
run_cp4_mmpbsa s23           outputs/2QKI_Cp4_calib_s23/snapshots            outputs/2QKI_Cp4_calib_s23/mmpbsa_results_fix4
run_cp4_mmpbsa s101          outputs/2QKI_Cp4_calib_s101/snapshots           outputs/2QKI_Cp4_calib_s101/mmpbsa_results_fix4
run_cp4_mmpbsa s19_reseed55  outputs/2QKI_Cp4_calib_s19_reseed55/snapshots   outputs/2QKI_Cp4_calib_s19_reseed55/mmpbsa_results_fix4

echo "=========================================="
echo "TRACK B PHASE 2 Cp4 PBSA BATCH END: $(date --iso-8601=seconds)"
echo "=========================================="
