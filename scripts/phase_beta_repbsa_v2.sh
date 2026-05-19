#!/bin/bash
# ============================================================================
# Phase β re-PBSA v2 — 10 systems with CONECT v55 patch
# ============================================================================
# Inputs : outputs/<sys>/snapshots_n25_postl387_patch_v2/  (CONECT v55)
# Outputs: outputs/<sys>/mmpbsa_results_postl387_v2/       (NEW; preserves
#          mmpbsa_results_postl387/ from #84 per R-7).
#
# Idempotent: skips entries whose mmpbsa_summary.json already exists.
# Pattern: mirrors scripts/phase_beta_mmpbsa_rerun.sh.
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

export UPDD_MMGBSA_PLATFORM=CUDA

# Log location can be overridden via T1_LOG env var (T1 Phase 1 runs use a
# dated dispatch log; default keeps the historical Phase β log path).
LOG="${T1_LOG:-${PROJ}/outputs/analysis/phase_beta_repbsa_v2_20260429.log}"
mkdir -p "$(dirname "$LOG")"

echo "==========================================" | tee -a "$LOG"
echo "PHASE β RE-PBSA v2 START: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "host=$(hostname) gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"

run_pbsa () {
    local label="$1"
    local snap_dir="$2"
    local outdir="$3"
    local target_id="$4"
    local ncaa_elem="$5"

    if [ -f "${outdir}/mmpbsa_summary.json" ]; then
        echo "SKIP ${label} (already done: ${outdir}/mmpbsa_summary.json)" | tee -a "$LOG"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        echo "MISS ${label} no pdb in ${snap_dir}" | tee -a "$LOG"
        return 1
    fi
    mkdir -p "${outdir}"
    local start
    start=$(date +%s)
    echo "[START ${label}] $(date --iso-8601=seconds) snap=${snap_dir} out=${outdir} target=${target_id} ncaa=${ncaa_elem}" | tee -a "$LOG"
    UPDD_MMGBSA_PLATFORM=CUDA "$PY" scripts/run_mmpbsa.py \
        --md_dir "${snap_dir}" \
        --outputdir "${outdir}" \
        --ncaa_elem "${ncaa_elem}" \
        --receptor_chain A --binder_chain B \
        --target_id "${target_id}" \
        --protocol 1traj \
        >> "$LOG" 2>&1
    local rc=$?
    local end
    end=$(date +%s)
    local elapsed=$((end - start))
    echo "[DONE ${label}] rc=${rc} elapsed=${elapsed}s" | tee -a "$LOG"
}

# Convenience wrapper: outputs/<sys>/snapshots_n25_postl387_patch_v2/
#                  →  outputs/<sys>/mmpbsa_results_postl387_v2/
run_sys () {
    local sys="$1"
    local target_id="$2"
    local ncaa_elem="$3"
    run_pbsa "${sys}" \
             "outputs/${sys}/snapshots_n25_postl387_patch_v2" \
             "outputs/${sys}/mmpbsa_results_postl387_v2" \
             "${target_id}" "${ncaa_elem}"
}

# v1-path wrapper (Case-A systems where atom_count < 99,999 → v55 patch
# numerically equivalent; aggregation stays under the v1 directory naming
# for symmetry with reps that pre-date the v55 patch).
run_sys_v1 () {
    local sys="$1"
    local target_id="$2"
    local ncaa_elem="$3"
    run_pbsa "${sys}" \
             "outputs/${sys}/snapshots_n25_postl387_patch" \
             "outputs/${sys}/mmpbsa_results_postl387" \
             "${target_id}" "${ncaa_elem}"
}

# ==========================================================================
# 10 systems — Phase β re-PBSA v2 (CONECT v55 unblocks #84 failures)
# ==========================================================================
echo "------------------------------------------" | tee -a "$LOG"
echo "Phase β re-PBSA v2 — 10 systems" | tee -a "$LOG"
echo "------------------------------------------" | tee -a "$LOG"

# 2QKI_Cp4 (8 seeds; Cp4 binder uses MTR for 4(1MeW))
run_sys 2QKI_Cp4_calib_s7            2QKI MTR
run_sys 2QKI_Cp4_calib_s19_reseed55  2QKI MTR
run_sys 2QKI_Cp4_calib_s23           2QKI MTR
run_sys 2QKI_Cp4_calib_s42           2QKI MTR
run_sys 2QKI_Cp4_calib_s83           2QKI MTR
run_sys 2QKI_Cp4_calib_s101          2QKI MTR
run_sys 2QKI_Cp4_calib_s163          2QKI MTR
run_sys 2QKI_Cp4_calib_s251          2QKI MTR

# 7TL8_MTR6 (2 INTACT seeds; s7 already done in #84 with v55-incompatible snap dir;
#            v55 path proved CONECT INTACT for s7 Cp4 too — but s7 in #84 used
#            the prior dir, not v2; aggregation will use #84's s7 unchanged.)
run_sys 7TL8_MTR6_calib_s19 7TL8 MTR
run_sys 7TL8_MTR6_calib_s42 7TL8 MTR

# ==========================================================================
# T1 N→30 expansion — Phase 1 (2026-05-18)
# ==========================================================================
# Skipped entries (existing summary) → no-op. New work runs and writes
# alongside existing reps for joint aggregation.
echo "------------------------------------------" | tee -a "$LOG"
echo "T1 N→30 expansion — Phase 1 (2026-05-18)" | tee -a "$LOG"
echo "------------------------------------------" | tee -a "$LOG"

# Phase 1A — Cp4 MD-ready seeds matching the published baseline cohort
# (GAFF2-unpatched MD per 2026-05-18 cohort audit). The other 4 seeds
# (s433/s523/s773/s991) carry ff14SB-patched MD and were archived to
# /media/san/ExpDATA/UPDD_proj_Backup/cp4_patched_cohort_20260518/ for
# separate downstream analysis. With these 3 reps the Cp4 cohort grows
# from n=8 (baseline) to n=11; 4 additional fresh MD launched separately
# under UPDD_MTR_AMBER14_PATCH=0 will bring it to n=15.
run_sys 2QKI_Cp4_calib_s127 2QKI MTR
run_sys 2QKI_Cp4_calib_s19  2QKI MTR
run_sys 2QKI_Cp4_calib_s199 2QKI MTR

# Phase 1B — 2QKH_MTR25 in-flight restart only (9 seeds, snapshots intact).
# Case A confirmed (45,358-45,828 atoms < 99,999) → v1 path. Brings the
# cohort from n=15 (already done) to n=24. The 3 MD-only seeds
# (s479/s997/s1549) were DEFERRED: their DCDs were generated by OpenMM 8.4
# producing 45,358-atom solvation, but no `_final.pdb` was written before
# the 2026-05-17 shutdown. The existing OpenMM-8.2 s7 _final.pdb has
# 45,828 atoms (different waterbox population), so the topology-proxy
# approach used for Cp4 GAFF2 seeds does not transfer here. Recovery
# requires either reconstructing topology from the per-seed
# `mdresult/*_universal.xml`, or re-running these 3 MDs from scratch —
# both are post-Phase-1 follow-ups.
run_sys_v1 2QKH_MTR25_calib_s317  2QKH MTR
run_sys_v1 2QKH_MTR25_calib_s367  2QKH MTR
run_sys_v1 2QKH_MTR25_calib_s421  2QKH MTR
run_sys_v1 2QKH_MTR25_calib_s827  2QKH MTR
run_sys_v1 2QKH_MTR25_calib_s881  2QKH MTR
run_sys_v1 2QKH_MTR25_calib_s937  2QKH MTR
run_sys_v1 2QKH_MTR25_calib_s1361 2QKH MTR
run_sys_v1 2QKH_MTR25_calib_s1423 2QKH MTR
run_sys_v1 2QKH_MTR25_calib_s1483 2QKH MTR

# Phase 1C — 2QKI_WT pilot (Cp4 reference branch) under postl387_v2.
# Two seeds matched to Cp4 baseline (s7, s23) for protocol validation
# before committing to full N=30 in Phase 3.
run_sys 2QKI_WT_calib_s7  2QKI none
run_sys 2QKI_WT_calib_s23 2QKI none

echo "==========================================" | tee -a "$LOG"
echo "PHASE β RE-PBSA v2 END: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
