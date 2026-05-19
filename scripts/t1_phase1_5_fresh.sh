#!/bin/bash
# ============================================================================
# T1 Phase 1.5 — Cp4 fresh GAFF2 + WT v2 expansion (2026-05-18)
# ============================================================================
# Scope:
#   Stage A: 4 Cp4 fresh GAFF2-MD seeds (s53, s89, s113, s173)
#            with Option A+B remediation: per-seed UNPATCHED params from
#            outputs/_archive/pre_amb14_patch_20260427/ + grep gate
#            verification BEFORE MD launch.
#   Stage B: 2 WT fresh MD seeds (s127, s199) — matched-integer with Cp4
#            GAFF2 new cohort (Phase 1 s127/s19/s199 + Phase 1.5 fresh 4).
#            No MTR patch concern (WT has no ncAA).
#   Stage C: WT v2 reextract for 6 existing-DCD seeds (s19, s42, s83, s101,
#            s163, s251) + 2 fresh (Stage B output) + Cp4 fresh 4 (Stage A
#            output). Plus restore PATCHED params/MTR_gaff2.xml on Cp4 fresh
#            seeds for PBSA scoring stage parity with baseline 11 cohort.
#   Stage D: 4-lane parallel MMPBSA on 4 Cp4 + 8 WT = 12 work units.
#   Stage E: Branched ΔΔG aggregation via utils/branched_ddg.py.
#
# Wall (revised after preflight scope c review):
#   Stage A 4 × 11h = 44h GPU (sequential MD)
#   Stage B 2 × 11h = 22h GPU
#   Stage C 8 × 30min = 4h CPU+light-GPU
#   Stage D 12 / 4-lane × ~7min = ~21 min wall (4 parallel)
#   Stage E aggregation ~5 min
#   TOTAL ~70-72h wall, single GPU.
#
# Cohort consistency (R-4 + R-17 + Path A remediation):
#   Cp4: MD-time = UNPATCHED (matches baseline 11 MD), PBSA-time = PATCHED
#        (matches baseline 11 scoring). Toggled per-seed via XML swap.
#   WT:  no MTR → no patch concern, automatic cohort consistency.
#
# Template: archive/2026-05-05_pre_v1_freeze/dispatch_scripts/system_specific/7tl8_mtr6_n8_expand.sh
# ============================================================================

set -u

# ----- env -----
PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh

# UPDD_MTR_AMBER14_PATCH=0 (defensive) — preflight noted this is inert at MD
# time but harmless. Cohort consistency is enforced via per-seed XML swap, not
# this env var.
export UPDD_MTR_AMBER14_PATCH=0
export UPDD_MMGBSA_PLATFORM=CPU  # CPU-forced to avoid /dev/nvidia0 contention with MD; affinity below isolates further

# ============================================================================
# CPU affinity 자동 분배 (2026-05-19 hardening)
# ----------------------------------------------------------------------------
# 9800X3D 8C/16T 같은 8-physical 시스템에서 MD (GPU launcher, 단일 thread bound)
# 와 MMPBSA (sander multi-thread CPU compute) 가 CPU thread pool 경쟁으로 인해
# GPU kernel dispatch 지연 + MD throughput 3× 저하 + WT MMPBSA stall 발생.
# 해결: 1 physical core 를 MD 전용 + 나머지 (N-1) 을 MMPBSA 에 할당.
# AMD/Intel SMT 패턴: physical core N = logical {N, N+PHYS_CORES}.
# 8 phys 초과 시 free scheduling 으로 fallback (충분한 자원).
# ============================================================================
PHYS_CORES=$(lscpu -p 2>/dev/null | grep -v "^#" | awk -F, '{print $2}' | sort -nu | wc -l)
LOGICAL_CORES=$(nproc)
if [ "${PHYS_CORES:-0}" -ge 4 ] && [ "${PHYS_CORES:-0}" -le 8 ]; then
    MD_AFFINITY="0,${PHYS_CORES}"
    MMPBSA_AFFINITY="1-$((PHYS_CORES - 1)),$((PHYS_CORES + 1))-$((LOGICAL_CORES - 1))"
    USE_AFFINITY=1
    MD_PREFIX="taskset -c ${MD_AFFINITY}"
    MMPBSA_PREFIX="taskset -c ${MMPBSA_AFFINITY}"
else
    USE_AFFINITY=0
    MD_PREFIX=""
    MMPBSA_PREFIX=""
fi

LAUNCH_TS=$(date +%s)
TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/t1_phase1_5_${TS_TAG}"
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"

log () {
    printf '[t1_phase1_5 %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "T1 Phase 1.5 START"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
log "logdir=$LOGDIR"
log "phys_cores=$PHYS_CORES  logical=$LOGICAL_CORES  affinity=${USE_AFFINITY}"
if [ "$USE_AFFINITY" -eq 1 ]; then
    log "  MD affinity:     ${MD_AFFINITY} (1 phys core)"
    log "  MMPBSA affinity: ${MMPBSA_AFFINITY} (${PHYS_CORES} - 1 = $((PHYS_CORES - 1)) phys cores)"
else
    log "  affinity off (free scheduling — system has ${PHYS_CORES} phys cores)"
fi
log "================================================================"

# Seed lists (preflight verified all integers free):
CP4_FRESH=(53 89 113 173)
WT_FRESH=(127 199)
WT_EXISTING=(19 42 83 101 163 251)

# ============================================================================
# Stage A — Cp4 fresh GAFF2 MD (4 seeds) with GATE 2 remediation
# ============================================================================
log "----------------------------------------------------------------"
log "Stage A: Cp4 fresh GAFF2-MD ${#CP4_FRESH[@]} seeds: ${CP4_FRESH[*]}"
log "----------------------------------------------------------------"

CP4_REF=s7
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_${CP4_REF}"
CP4_ARCHIVE_REF="${PROJ}/outputs/_archive/pre_amb14_patch_20260427/2QKI_Cp4_calib_${CP4_REF}"

# Pre-flight: verify archive UNPATCHED reference exists
if [ ! -f "${CP4_ARCHIVE_REF}/params/MTR_gaff2.xml" ]; then
    abort "Archive UNPATCHED MTR_gaff2.xml missing: ${CP4_ARCHIVE_REF}/params/MTR_gaff2.xml"
fi

for seed_int in "${CP4_FRESH[@]}"; do
    seed="s${seed_int}"
    seedir="${PROJ}/outputs/2QKI_Cp4_calib_${seed}"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  $seed: SKIP (DCD exists $(du -h ${seedir}/mdresult/2QKI_Cp4_restrained.dcd | cut -f1))"
        continue
    fi

    log "  $seed: setup dir"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    # Mirror non-MTR input files
    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && cp -n "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    # Step A1 — Install UNPATCHED params/ from archive (q_N=-0.8938)
    log "  $seed: copy UNPATCHED params from archive"
    cp -n "${CP4_ARCHIVE_REF}/params/MTR_gaff2.xml" "${seedir}/params/"
    [ -f "${CP4_ARCHIVE_REF}/params/MTR_hydrogens.xml" ] && cp -n "${CP4_ARCHIVE_REF}/params/MTR_hydrogens.xml" "${seedir}/params/"
    [ -f "${CP4_ARCHIVE_REF}/params/MTR_params_manifest.json" ] && cp -n "${CP4_ARCHIVE_REF}/params/MTR_params_manifest.json" "${seedir}/params/"

    # Update manifest xml_path to point locally
    if [ -f "${seedir}/params/MTR_params_manifest.json" ]; then
        $PY -c "
import json, sys
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"
    fi

    # Step A2 — GATE-B: verify q_N is UNPATCHED before MD launch
    xml="${seedir}/params/MTR_gaff2.xml"
    qN=$(grep -m1 '<Atom name="N" ' "$xml" 2>/dev/null | grep -oP 'charge="\K[^"]+')
    case "$qN" in
        -0.8938|-0.8937999978048875)
            log "  $seed: GATE-B PASS q_N=$qN (UNPATCHED ✓)"
            ;;
        *)
            abort "$seed: GATE-B FAIL q_N=$qN (expected -0.8938 UNPATCHED); halting Phase 1.5"
            ;;
    esac

    # Step A3 — MD launch (env defensive)
    log "  $seed: MD launch (5ns, cyclic_ss, seed=${seed_int})  affinity=[${MD_AFFINITY:-free}]"
    conda activate md_simulation
    UPDD_MTR_AMBER14_PATCH=0 UPDD_MMGBSA_PLATFORM=CUDA \
        $MD_PREFIX "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${seedir}/_md_input" \
            --outputdir "${seedir}/mdresult" \
            --params_manifest "${seedir}/params/MTR_params_manifest.json" \
            --steps 2500000 \
            --topology cyclic_ss \
            --binder_chain B \
            --graph_policy strict \
            --target_id 2QKI \
            --dt_fs 2.0 \
            --platform CUDA \
            --seed "$seed_int" \
            --ncaa_label MTR \
            --ncaa_code MTR \
            > "${LOGDIR}/md_cp4_${seed}.log" 2>&1
    rc=$?
    conda deactivate

    if [ "$rc" -ne 0 ] || [ ! -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  $seed: MD FAIL rc=$rc (see ${LOGDIR}/md_cp4_${seed}.log)"
        continue
    fi
    log "  $seed: MD DONE"

    # Step A4 — Restore PATCHED params/ for PBSA scoring stage parity
    log "  $seed: install PATCHED params from current ${CP4_REF_DIR}/params/ (PBSA-side cohort match)"
    cp -f "${CP4_REF_DIR}/params/MTR_gaff2.xml" "${seedir}/params/"
    [ -f "${CP4_REF_DIR}/params/MTR_hydrogens.xml" ] && cp -f "${CP4_REF_DIR}/params/MTR_hydrogens.xml" "${seedir}/params/"
    [ -f "${CP4_REF_DIR}/params/MTR_params_manifest.json" ] && cp -f "${CP4_REF_DIR}/params/MTR_params_manifest.json" "${seedir}/params/"
    if [ -f "${seedir}/params/MTR_params_manifest.json" ]; then
        $PY -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"
    fi

    # Step A5 — Verify PATCHED q_N for PBSA stage
    qN=$(grep -m1 '<Atom name="N" ' "${seedir}/params/MTR_gaff2.xml" 2>/dev/null | grep -oP 'charge="\K[^"]+')
    case "$qN" in
        -0.4157)
            log "  $seed: PBSA-side q_N=$qN (PATCHED ✓ matches baseline)"
            ;;
        *)
            log "  $seed: WARN PBSA-side q_N=$qN (expected -0.4157 PATCHED)"
            ;;
    esac
done

log "----------------------------------------------------------------"
log "Stage A complete. Cp4 fresh MD: ${#CP4_FRESH[@]} seeds processed."

# ============================================================================
# Stage B — WT fresh MD (2 seeds, no MTR patch concern)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage B: WT fresh MD ${#WT_FRESH[@]} seeds: ${WT_FRESH[*]}"
log "----------------------------------------------------------------"

WT_REF=s7
WT_REF_DIR="${PROJ}/outputs/2QKI_WT_calib_${WT_REF}"

for seed_int in "${WT_FRESH[@]}"; do
    seed="s${seed_int}"
    seedir="${PROJ}/outputs/2QKI_WT_calib_${seed}"

    if [ -f "${seedir}/mdresult/2QKI_WT_restrained.dcd" ]; then
        log "  $seed: SKIP (DCD exists $(du -h ${seedir}/mdresult/2QKI_WT_restrained.dcd | cut -f1))"
        continue
    fi

    log "  $seed: setup dir + mirror from ${WT_REF}"
    mkdir -p "${seedir}/_md_input" "${seedir}/mdresult"
    # WT has no params/ — mirror _md_input/
    if [ -d "${WT_REF_DIR}/_md_input" ]; then
        cp -rn "${WT_REF_DIR}/_md_input/." "${seedir}/_md_input/"
    fi

    log "  $seed: MD launch (5ns, cyclic_ss, seed=${seed_int})  affinity=[${MD_AFFINITY:-free}]"
    conda activate md_simulation
    UPDD_MMGBSA_PLATFORM=CUDA \
        $MD_PREFIX "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${seedir}/_md_input" \
            --outputdir "${seedir}/mdresult" \
            --steps 2500000 \
            --topology cyclic_ss \
            --binder_chain B \
            --graph_policy strict \
            --target_id 2QKI \
            --dt_fs 2.0 \
            --platform CUDA \
            --seed "$seed_int" \
            > "${LOGDIR}/md_wt_${seed}.log" 2>&1
    rc=$?
    conda deactivate

    if [ "$rc" -ne 0 ] || [ ! -f "${seedir}/mdresult/2QKI_WT_restrained.dcd" ]; then
        log "  $seed: MD FAIL rc=$rc"
        continue
    fi
    log "  $seed: MD DONE"
done

log "----------------------------------------------------------------"
log "Stage B complete. WT fresh MD: ${#WT_FRESH[@]} seeds processed."

# ============================================================================
# Stage C — Snapshot reextract (v2 path for all)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage C: snapshot reextract v2 path (Cp4 4 fresh + WT 8 total)"
log "----------------------------------------------------------------"
conda activate qmmm

# Build temporary TARGETS list via Python-as-glue. Use --filter on the existing
# reextract.py since editing its TARGETS array would touch a different concern.
for system_seed in "2QKI_Cp4_calib_s53" "2QKI_Cp4_calib_s89" "2QKI_Cp4_calib_s113" "2QKI_Cp4_calib_s173" \
                   "2QKI_WT_calib_s127" "2QKI_WT_calib_s199" \
                   "2QKI_WT_calib_s19" "2QKI_WT_calib_s42" "2QKI_WT_calib_s83" \
                   "2QKI_WT_calib_s101" "2QKI_WT_calib_s163" "2QKI_WT_calib_s251"; do
    target="${PROJ}/outputs/${system_seed}/snapshots_n25_postl387_patch_v2"
    if [ -d "$target" ] && [ "$(ls $target/*.pdb 2>/dev/null | wc -l)" -ge 25 ]; then
        continue  # already extracted
    fi
    log "  reextract ${system_seed}  affinity=[${MMPBSA_AFFINITY:-free}]"
    CUDA_VISIBLE_DEVICES="" $MMPBSA_PREFIX $PY scripts/phase_beta_repbsa_v2_reextract.py \
        --filter "${system_seed}" \
        --out-json "${LOGDIR}/reextract_${system_seed}.json" \
        >> "${LOGDIR}/reextract.log" 2>&1
done

log "Stage C complete."

# ============================================================================
# Stage D — 4-lane parallel MMPBSA (12 work units)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage D: MMPBSA parallel sweep (4 lanes)"
log "----------------------------------------------------------------"
LANES=4
T1_LOG="${LOGDIR}/mmpbsa.log"

run_one_pbsa () {
    local label="$1" snap_dir="$2" out_dir="$3" target_id="$4" ncaa_elem="$5"
    if [ -f "${out_dir}/mmpbsa_summary.json" ]; then
        printf '[%s SKIP %s]\n' "$(date +%H:%M:%S)" "$label" >> "$T1_LOG"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        printf '[%s MISS %s] no pdb\n' "$(date +%H:%M:%S)" "$label" >> "$T1_LOG"
        return 1
    fi
    mkdir -p "$out_dir"
    local seed_log="${LOGDIR}/mmpbsa_${label}.log"
    local start; start=$(date +%s)
    printf '[%s START %s]\n' "$(date +%H:%M:%S)" "$label" >> "$T1_LOG"
    CUDA_VISIBLE_DEVICES="" UPDD_MMGBSA_PLATFORM=CPU $MMPBSA_PREFIX "$PY" scripts/run_mmpbsa.py \
        --md_dir "$snap_dir" --outputdir "$out_dir" \
        --ncaa_elem "$ncaa_elem" --receptor_chain A --binder_chain B \
        --target_id "$target_id" --protocol 1traj \
        > "$seed_log" 2>&1
    local rc=$?
    local end; end=$(date +%s)
    printf '[%s DONE  %s] rc=%d elapsed=%ds\n' "$(date +%H:%M:%S)" "$label" "$rc" "$((end-start))" >> "$T1_LOG"
}

WORK=(
    "Cp4_s53|outputs/2QKI_Cp4_calib_s53/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s53/mmpbsa_results_postl387_v2|2QKI|MTR"
    "Cp4_s89|outputs/2QKI_Cp4_calib_s89/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s89/mmpbsa_results_postl387_v2|2QKI|MTR"
    "Cp4_s113|outputs/2QKI_Cp4_calib_s113/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s113/mmpbsa_results_postl387_v2|2QKI|MTR"
    "Cp4_s173|outputs/2QKI_Cp4_calib_s173/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s173/mmpbsa_results_postl387_v2|2QKI|MTR"
    "WT_s19|outputs/2QKI_WT_calib_s19/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s19/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s42|outputs/2QKI_WT_calib_s42/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s42/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s83|outputs/2QKI_WT_calib_s83/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s83/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s101|outputs/2QKI_WT_calib_s101/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s101/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s127|outputs/2QKI_WT_calib_s127/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s127/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s163|outputs/2QKI_WT_calib_s163/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s163/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s199|outputs/2QKI_WT_calib_s199/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s199/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s251|outputs/2QKI_WT_calib_s251/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s251/mmpbsa_results_postl387_v2|2QKI|none"
)

running=0
for entry in "${WORK[@]}"; do
    IFS='|' read -r label snap out tgt ncaa <<< "$entry"
    run_one_pbsa "$label" "$snap" "$out" "$tgt" "$ncaa" &
    running=$((running + 1))
    if [ "$running" -ge "$LANES" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait

log "Stage D complete."

# ============================================================================
# Stage E — Branched ΔΔG aggregation
# ============================================================================
log "----------------------------------------------------------------"
log "Stage E: Branched ΔΔG aggregation via utils/branched_ddg.py"
log "----------------------------------------------------------------"
$PY -m utils.branched_ddg \
    --wt "outputs/2QKI_WT_calib_*/mmpbsa_results_postl387_v2" \
    --variant "outputs/2QKI_Cp4_calib_*/mmpbsa_results_postl387_v2" \
    --output "${LOGDIR}/branched_ddg/" \
    >> "$DISPATCH_LOG" 2>&1 || log "  WARN branched_ddg returned non-zero"

# ============================================================================
# Summary
# ============================================================================
ELAPSED=$(( $(date +%s) - LAUNCH_TS ))
ELAPSED_H=$(echo "scale=2; $ELAPSED / 3600" | bc)
log "================================================================"
log "T1 Phase 1.5 COMPLETE  elapsed=${ELAPSED}s (~${ELAPSED_H}h)"
log "  Stage A Cp4 fresh MD:    ${#CP4_FRESH[@]} seeds"
log "  Stage B WT fresh MD:     ${#WT_FRESH[@]} seeds"
log "  Stage C Snapshot extract: 12 dirs"
log "  Stage D MMPBSA:          12 work units (4 lane)"
log "  Stage E Branched ΔΔG:    outputs/analysis/t1_phase1_5_${TS_TAG}/branched_ddg/"
log ""
log "Next: review branched_ddg result for Phase 1.5 final tier classification."
log "================================================================"
