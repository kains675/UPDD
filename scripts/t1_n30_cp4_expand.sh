#!/bin/bash
# ============================================================================
# T1 #1 — N→30 expansion dispatcher (Option β, V100 sequential, 2026-05-22)
# ============================================================================
# Source: byte-identical clone of scripts/t1_phase1_5_fresh.sh per
#         verdict_t1_n30_expansion_20260522.md §C2. Only the seed arrays
#         and Stage C+D WORK list construction differ; the MD-time UNPATCHED
#         + PBSA-time PATCHED 2-layer protocol (lines 144-178 of source) is
#         preserved verbatim.
#
# Pre-launch conditions satisfied (verdict §C1-C5):
#   C1: Option β PRIMARY — Cp4 + WT both N→30 symmetric.
#   C2: Byte-identical clone of t1_phase1_5_fresh.sh; this header + seed
#       arrays + Stage C+D WORK construction differ; algorithmic core
#       preserved.
#   C3: Multi-outcome pre-registration —
#       outputs/paper1/v3_pending/t1_n30_pre_registration.md
#   C4: Deterministic numpy RNG seed selection —
#       outputs/analysis/t1_n30_seed_selection/t1_n30_seed_manifest.json
#       (SHA-256 34df5975e6958ca2f885d14bcf93ff2fc307495ce2761a8e4194d654c66ca027,
#        SeedSequence entropy 20260522).
#   C5: V100 sequential allocation — DO NOT launch concurrent with v0.9 γ
#       production lane on the same V100. Wait for v09_gamma_production.sh
#       to complete (or stop) before launching.
#
# Scope:
#   Stage A: 17 Cp4 fresh GAFF2-MD seeds (CP4_FRESH below) with 2-layer
#            MD-UNPATCHED + PBSA-PATCHED + GATE-B verification.
#   Stage B: 20 WT fresh MD seeds (WT_FRESH below). No MTR → no patch.
#   Stage C: postl387_v2 snapshot reextract on 37 new + 23 existing-verify
#            dirs (Phase 1.5 cohort dirs short-circuit if already extracted).
#   Stage D: 4-lane parallel MMPBSA on 30 Cp4 + 30 WT = 60 work units.
#   Stage E: Branched ΔΔG aggregation (N=30 paired).
#
# Wall envelope (V100 sequential after γ production):
#   Stage A 17 × 11h = ~187h V100
#   Stage B 20 × 11h = ~220h V100
#   Stage C 37 × 30min = ~18h CPU+light-GPU
#   Stage D 60 / 4-lane × ~7min = ~110min
#   Stage E ~5 min
#   TOTAL ~430h (~18 days V100). Interim N=20 checkpoint per
#   pre-registration §4 (after 7 Cp4 + 10 WT fresh complete).
#
# Cohort consistency (R-4 + R-17 + Path A remediation, identical to source):
#   Cp4: MD-time = UNPATCHED, PBSA-time = PATCHED. Per-seed XML swap.
#   WT:  no MTR → automatic cohort consistency.
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

# Seed lists (deterministic numpy RNG, SeedSequence entropy=20260522,
# manifest SHA-256 34df5975e6958ca2f885d14bcf93ff2fc307495ce2761a8e4194d654c66ca027
# — outputs/analysis/t1_n30_seed_selection/t1_n30_seed_manifest.json):
CP4_FRESH=(31 91 152 167 217 288 296 365 522 579 618 626 767 790 793 841 955)
WT_FRESH=(26 112 278 292 316 387 422 426 482 499 501 507 613 666 687 688 743 866 900 909)
# Phase 1.5 existing cohort (Cp4 13 effective; s113/s173 LOAD_FAILED excluded;
# integer 55 corresponds to directory s19_reseed55):
CP4_EXISTING=(7 19 23 42 53 55 83 89 101 127 163 199 251)
# WT existing cohort (n=10):
WT_EXISTING=(7 19 23 42 83 101 127 163 199 251)

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
log "Stage C: snapshot reextract postl387_v2 (37 fresh + 23 existing-verify)"
log "----------------------------------------------------------------"
conda activate qmmm

# Build reextract target list from arrays. Already-extracted dirs (≥25 pdbs)
# are skipped at runtime — Phase 1.5 cohort dirs short-circuit.
REEXTRACT_TARGETS=()
for s in "${CP4_FRESH[@]}"; do REEXTRACT_TARGETS+=("2QKI_Cp4_calib_s${s}"); done
for s in "${WT_FRESH[@]}"; do REEXTRACT_TARGETS+=("2QKI_WT_calib_s${s}"); done
for s in "${CP4_EXISTING[@]}"; do
    if [ "$s" -eq 55 ]; then
        REEXTRACT_TARGETS+=("2QKI_Cp4_calib_s19_reseed55")
    else
        REEXTRACT_TARGETS+=("2QKI_Cp4_calib_s${s}")
    fi
done
for s in "${WT_EXISTING[@]}"; do REEXTRACT_TARGETS+=("2QKI_WT_calib_s${s}"); done

for system_seed in "${REEXTRACT_TARGETS[@]}"; do
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

WORK=()
# Cp4 fresh 17 + Cp4 existing 13 = 30 Cp4 work units (MTR ncaa)
for s in "${CP4_FRESH[@]}" "${CP4_EXISTING[@]}"; do
    if [ "$s" -eq 55 ]; then
        WORK+=("Cp4_s19_reseed55|outputs/2QKI_Cp4_calib_s19_reseed55/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s19_reseed55/mmpbsa_results_postl387_v2|2QKI|MTR")
    else
        WORK+=("Cp4_s${s}|outputs/2QKI_Cp4_calib_s${s}/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s${s}/mmpbsa_results_postl387_v2|2QKI|MTR")
    fi
done
# WT fresh 20 + WT existing 10 = 30 WT work units (no ncaa)
for s in "${WT_FRESH[@]}" "${WT_EXISTING[@]}"; do
    WORK+=("WT_s${s}|outputs/2QKI_WT_calib_s${s}/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s${s}/mmpbsa_results_postl387_v2|2QKI|none")
done
log "Stage D: ${#WORK[@]} work units constructed (Cp4 30 + WT 30)"

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
