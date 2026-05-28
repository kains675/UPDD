#!/usr/bin/env bash
# ============================================================================
# v09_gamma_rescue.sh — UPDD.py Pass 1b (reseed) + Pass 2 (1fs) fallback 포팅
# ----------------------------------------------------------------------------
# Standalone v09_gamma_production.sh 는 reseed/1fs fallback 미탑재 (full UPDD.py
# L1740-1857 에만 존재) → 폭발 seed 가 그냥 skip 됨. 이 script 가 그 fallback 을
# 포팅하여 exploded/incomplete Cp4 γ seed 를 rescue.
#
# Fallback ladder (UPDD.py 와 동일):
#   Pass 1b: seed+13 (prime offset) 으로 2fs 재시도 (Langevin RNG stream 변경 →
#            bad-draw NaN 회복. UPDD.py L1782: new_seed = seed + 13)
#   Pass 2:  여전히 폭발 시 dt=1fs 재시도 (cyclic_ss 는 UPDD.py 에서 1fs gated
#            아니지만, rescue 에서는 안전망으로 포함 — reseed 도 실패 시 최후수단)
#
# R-7 준수: 폭발 artifact 는 삭제 X, mdresult/_archive/<ts>_pre_reseed/ 로 이동.
#
# Device: host 5070 Ti (MD 에서 V100 보다 2× 빠름, Finding A). UPDD_RESCUE_GPU 로 override.
#
# 사용법:
#   bash scripts/v09_gamma_rescue.sh                 # 모든 exploded/incomplete seed auto-rescue
#   bash scripts/v09_gamma_rescue.sh s23             # 특정 seed 만
#   UPDD_RESCUE_GPU=0 bash scripts/v09_gamma_rescue.sh s23
# ============================================================================
set -u

PROJ=/home/san/UPDD_proj
cd "$PROJ"
source /home/san/miniconda3/etc/profile.d/conda.sh

PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
RESCUE_GPU="${UPDD_RESCUE_GPU:-0}"   # host 5070 Ti default (MD 는 host 가 빠름)

# Production cohort: label → original seed int
declare -A SEED_INT=(
    [s7]=7 [s23]=23 [s42]=42 [s83]=83
    [s101]=101 [s163]=163 [s251]=251 [s19_reseed55]=55
)

LOGDIR=$(ls -td "${PROJ}"/outputs/analysis/v09_gamma_production_* 2>/dev/null | head -1)
[ -z "$LOGDIR" ] && LOGDIR="${PROJ}/outputs/analysis/v09_gamma_rescue_$(date +%Y%m%d_%H%M%S)" && mkdir -p "$LOGDIR"
RESCUE_LOG="${LOGDIR}/rescue_dispatch.log"

log() { printf '[v09_rescue %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$RESCUE_LOG"; }

# 대상 결정: 인자 있으면 그 seed, 없으면 exploded/incomplete 전수
if [ "$#" -ge 1 ]; then
    TARGETS=("$@")
else
    TARGETS=()
    for label in "${!SEED_INT[@]}"; do
        seedir="${PROJ}/outputs/2QKI_Cp4_gamma_calib_${label}"
        # final.pdb 없고 (EXPLODED 마커 있거나 mdresult 존재) → rescue 대상
        if [ ! -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ] && [ -d "${seedir}/mdresult" ]; then
            TARGETS+=("$label")
        fi
    done
fi

log "================================================================"
log "v0.9 γ rescue START — targets: ${TARGETS[*]:-none}"
log "device: GPU ${RESCUE_GPU} (host 5070 Ti default — MD 는 host 가 V100 2× 빠름)"
log "================================================================"

[ "${#TARGETS[@]}" -eq 0 ] && { log "rescue 대상 없음 (모든 seed 완료 또는 미시작). 종료."; exit 0; }

archive_exploded () {
    # R-7: 폭발 artifact 를 _archive/<ts>_pre_<tag>/ 로 이동 (삭제 X)
    local seedir="$1"; local tag="$2"
    local mdr="${seedir}/mdresult"
    local adir="${mdr}/_archive/$(date +%Y%m%d_%H%M%S)_pre_${tag}"
    mkdir -p "$adir"
    # UPDD.py L1763-1769 의 archive pattern
    for pat in _EXPLODED_dt2fs.log _EXPLODED_dt1fs.log _EXPLODED_MD.log \
               _EXPLODED_NVT.log _EXPLODED_COLD_NVT.log _EXPLODED_MIN.log \
               _EXPLODED_CYCLIC_MIN.log _EXPLODED_CYCLIC_RELAX.log \
               _partial_md.pdb _partial_nvt.pdb _final.pdb _md.log _restrained.dcd; do
        for f in "${mdr}"/*"${pat}"; do
            [ -e "$f" ] && mv "$f" "$adir/" 2>/dev/null
        done
    done
    log "    archived exploded artifacts → ${adir#$PROJ/}"
}

run_md () {
    # 단일 seed MD 실행 (host or override GPU)
    local seedir="$1"; local seed="$2"; local dt="$3"; local label="$4"; local pass="$5"
    log "    $label: MD launch (Pass $pass, dt=${dt}fs, seed=${seed}, gpu=${RESCUE_GPU})"
    conda activate md_simulation
    CUDA_VISIBLE_DEVICES="${RESCUE_GPU}" \
        UPDD_NCAA_AMBER14_PATCH=0 \
        UPDD_MMGBSA_PLATFORM=CUDA \
        "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${seedir}/_md_input" \
            --outputdir "${seedir}/mdresult" \
            --params_manifest "${seedir}/params/MTR_params_manifest.json" \
            --steps 2500000 \
            --topology cyclic_ss \
            --binder_chain B \
            --graph_policy strict \
            --target_id 2QKI \
            --dt_fs "$dt" \
            --platform CUDA \
            --seed "$seed" \
            --ncaa_label MTR \
            --ncaa_code MTR \
            > "${LOGDIR}/rescue_md_${label}_pass${pass}.log" 2>&1
    local rc=$?
    conda deactivate
    return $rc
}

for label in "${TARGETS[@]}"; do
    seed_int="${SEED_INT[$label]:-}"
    [ -z "$seed_int" ] && { log "  $label: unknown seed label, skip"; continue; }
    seedir="${PROJ}/outputs/2QKI_Cp4_gamma_calib_${label}"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $label: SKIP (이미 final.pdb 존재)"
        continue
    fi
    if [ ! -d "${seedir}/_md_input" ] || [ ! -f "${seedir}/params/MTR_params_manifest.json" ]; then
        log "  $label: SKIP (setup 미비 — _md_input/params 부재, host launcher 가 아직 도달 안함)"
        continue
    fi

    log "  $label: rescue 시작 (original seed=${seed_int})"

    # ── Pass 1b: reseed (seed + 13) at 2fs (UPDD.py L1782) ──
    archive_exploded "$seedir" "reseed"
    new_seed=$((seed_int + 13))
    run_md "$seedir" "$new_seed" "2.0" "$label" "1b-reseed"
    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $label: ✓ Pass 1b reseed (seed=${new_seed}, 2fs) SUCCESS"
        continue
    fi
    log "  $label: Pass 1b reseed 실패 — 1fs 최후수단"

    # ── Pass 2: dt=1fs (UPDD.py L1850, rescue 안전망) ──
    archive_exploded "$seedir" "dt1fs"
    run_md "$seedir" "$new_seed" "1.0" "$label" "2-dt1fs"
    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $label: ✓ Pass 2 (seed=${new_seed}, 1fs) SUCCESS"
    else
        log "  $label: ✗ rescue 최종 실패 (2fs reseed + 1fs 모두 폭발) — seed 영구 제외 후보"
    fi
done

log "================================================================"
log "γ rescue 완료. host launcher 의 Stage A summary 가 rescued final.pdb 카운트."
log "================================================================"
