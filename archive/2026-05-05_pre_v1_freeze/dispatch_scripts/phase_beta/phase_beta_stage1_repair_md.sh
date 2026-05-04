#!/bin/bash
# Phase β Stage 1-Repair — MD reseed of 5 failed MTR-bearing seeds
#
# Scope:
#   7TL8_MTR6_calib_s19 — NVT-exploded (Particle coordinate is NaN), no DCD
#   7TL8_MTR6_calib_s42 — NVT-exploded, no DCD
#   2QKH_MTR25_calib_s7  — Production-MD exploded at step 0, partial_md.pdb is 119-byte stub
#   2QKH_MTR25_calib_s19 — same failure mode
#   2QKH_MTR25_calib_s42 — same failure mode
#
# Topology:
#   7TL8_MTR6   — linear (binder chain B, 13-mer with single CYS13)
#   2QKH_MTR25  — linear (GIP 32-mer)
#
# Pass 1: original seed
# Pass 1b: on NVT/MD failure, reseed with seed+13 (canonical UPDD precedent)
# Steps: 2,500,000 × 2 fs = 5 ns (matches successful 7TL8_MTR6_s7 reference)
# Platform: CUDA, dt=2.0 fs
# ncAA: MTR (post Phase α + Option C, Σq=0 on Residue/MTR/NMTR templates)
#
# Sequential execution (single GPU, no contention).

set -u
PY=/home/san/miniconda3/envs/md_simulation/bin/python
export PATH=/home/san/miniconda3/envs/md_simulation/bin:$PATH
export UPDD_MTR_AMBER14_PATCH=1
cd /home/san/UPDD_proj

LAUNCH_TS=$(date +%s)
GPU_BUDGET_H=14
STEPS=2500000   # 5 ns at 2 fs
N_SNAP=25

LOGDIR=/home/san/UPDD_proj/outputs/analysis/phase_beta_stage1_repair_run_logs
mkdir -p "$LOGDIR"

log_status() {
    echo "[stage1_repair $(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/dispatch.log"
}

check_budget() {
    local now=$(date +%s)
    local elapsed_h=$(( (now - LAUNCH_TS) / 3600 ))
    if [ "$elapsed_h" -ge "$GPU_BUDGET_H" ]; then
        log_status "GPU budget ${GPU_BUDGET_H}h exceeded (elapsed=${elapsed_h}h) — HALT"
        return 1
    fi
    return 0
}

run_one_md() {
    local outdir="$1"; local target_id="$2"; local seed="$3"; local tag="$4"
    local manifest="${outdir}/_md_input/${target_id##*_}_params_manifest.json"
    # use canonical params/ manifest (LOCAL xml_path), NOT _md_input/ manifest (/tmp pointer)
    local manifest_canonical="${outdir}/params/MTR_params_manifest.json"
    local logf="${LOGDIR}/${tag}.log"

    log_status "=== ${tag} seed=${seed} 5ns linear MTR MD ==="

    UPDD_MMGBSA_PLATFORM=CUDA $PY utils/run_restrained_md.py \
        --inputdir "${outdir}/_md_input" \
        --outputdir "${outdir}/mdresult" \
        --params_manifest "$manifest_canonical" \
        --steps "$STEPS" \
        --topology linear \
        --binder_chain B \
        --graph_policy strict \
        --target_id "$target_id" \
        --dt_fs 2.0 \
        --platform CUDA \
        --seed "$seed" \
        --ncaa_label MTR \
        --ncaa_code MTR \
        > "$logf" 2>&1

    if grep -qE "성공 \((ncAA|야생형) MD\) : 1 개|부분 성공 \(Recovery\) : 1 개|PARTIAL_SUCCESS|최종 구조 저장" "$logf"; then
        return 0
    fi
    return 1
}

run_one_full() {
    # outdir, target_id, seed_orig, tag_stem
    local outdir="$1"; local target_id="$2"; local seed_orig="$3"; local stem="$4"
    local dcd="${outdir}/mdresult/${stem}_restrained.dcd"

    if [ -f "$dcd" ] && [ "$(stat -c %s "$dcd")" -gt "100000000" ]; then
        log_status "SKIP ${stem}_s${seed_orig} (DCD exists, $(du -h "$dcd" | cut -f1))"
        return 0
    fi

    check_budget || return 2

    # Clear stale partial artifacts before retry
    rm -f "${outdir}/mdresult/"*EXPLODED* "${outdir}/mdresult/"*_partial_* \
          "${outdir}/mdresult/"*_final.pdb "${outdir}/mdresult/"*_md.log 2>/dev/null

    # Pass 1: original seed
    run_one_md "$outdir" "$target_id" "$seed_orig" "${stem}_s${seed_orig}_pass1" \
        && return 0

    log_status "WARN ${stem}_s${seed_orig} pass1 failed -> seed=$((seed_orig+13)) retry"
    rm -f "${outdir}/mdresult/"*EXPLODED* "${outdir}/mdresult/"*_partial_* \
          "${outdir}/mdresult/"*_final.pdb "${outdir}/mdresult/"*_md.log 2>/dev/null
    check_budget || return 2

    local reseed=$((seed_orig+13))
    run_one_md "$outdir" "$target_id" "$reseed" "${stem}_s${seed_orig}_pass1b_reseed${reseed}" \
        && return 0

    log_status "FAIL ${stem}_s${seed_orig} both passes exploded"
    return 1
}

log_status "========================================================"
log_status "Phase β Stage 1-Repair MD reseed dispatch"
log_status "  GPU budget: ${GPU_BUDGET_H} h"
log_status "  Steps:      ${STEPS} (5 ns at 2 fs)"
log_status "  Env:        UPDD_MTR_AMBER14_PATCH=${UPDD_MTR_AMBER14_PATCH}"
log_status "  Platform:   CUDA, md_simulation env (OpenMM 8.2)"
log_status "========================================================"

# 7TL8_MTR6 reseeds (s19, s42) — NVT-exploded
run_one_full outputs/7TL8_MTR6_calib_s19  7TL8  19  7TL8_MTR6
run_one_full outputs/7TL8_MTR6_calib_s42  7TL8  42  7TL8_MTR6

# 2QKH_MTR25 reseeds (s7, s19, s42) — Production-MD exploded at step 0
run_one_full outputs/2QKH_MTR25_calib_s7   2QKH   7  2QKH_MTR25
run_one_full outputs/2QKH_MTR25_calib_s19  2QKH  19  2QKH_MTR25
run_one_full outputs/2QKH_MTR25_calib_s42  2QKH  42  2QKH_MTR25

log_status "========================================================"
log_status "MD reseed phase complete"
log_status "========================================================"

# Per-system DCD verification summary
$PY << 'PYEOF'
import os, json
sys_specs = [
    ('7TL8_MTR6_s19', '/home/san/UPDD_proj/outputs/7TL8_MTR6_calib_s19/mdresult/7TL8_MTR6_restrained.dcd',
     '/home/san/UPDD_proj/outputs/7TL8_MTR6_calib_s19/mdresult/7TL8_MTR6_final.pdb'),
    ('7TL8_MTR6_s42', '/home/san/UPDD_proj/outputs/7TL8_MTR6_calib_s42/mdresult/7TL8_MTR6_restrained.dcd',
     '/home/san/UPDD_proj/outputs/7TL8_MTR6_calib_s42/mdresult/7TL8_MTR6_final.pdb'),
    ('2QKH_MTR25_s7',  '/home/san/UPDD_proj/outputs/2QKH_MTR25_calib_s7/mdresult/2QKH_MTR25_restrained.dcd',
     '/home/san/UPDD_proj/outputs/2QKH_MTR25_calib_s7/mdresult/2QKH_MTR25_final.pdb'),
    ('2QKH_MTR25_s19', '/home/san/UPDD_proj/outputs/2QKH_MTR25_calib_s19/mdresult/2QKH_MTR25_restrained.dcd',
     '/home/san/UPDD_proj/outputs/2QKH_MTR25_calib_s19/mdresult/2QKH_MTR25_final.pdb'),
    ('2QKH_MTR25_s42', '/home/san/UPDD_proj/outputs/2QKH_MTR25_calib_s42/mdresult/2QKH_MTR25_final.pdb',
     '/home/san/UPDD_proj/outputs/2QKH_MTR25_calib_s42/mdresult/2QKH_MTR25_final.pdb'),
]
report = []
for tag, dcd, top in sys_specs:
    rec = {'tag': tag, 'dcd': dcd, 'exists': os.path.exists(dcd)}
    if rec['exists']:
        rec['size_mb'] = round(os.path.getsize(dcd)/1e6, 1)
    report.append(rec)
print(json.dumps(report, indent=2))
with open('/home/san/UPDD_proj/outputs/analysis/phase_beta_stage1_repair_md_status.json', 'w') as f:
    json.dump(report, f, indent=2)
PYEOF
log_status "MD status JSON saved"
