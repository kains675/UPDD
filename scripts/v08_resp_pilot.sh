#!/bin/bash
# ============================================================================
# v0.8 Khoury RESP migration pilot — 2QKI_Cp4 ↔ WT first quantitative test
# ============================================================================
# Per tier-1 cross-system sweep §Priority 1: validate Khoury full RESP MTR
# (= OMW in Khoury naming) recovers the literature sign for the Cp4 ↔ WT paired
# ΔΔG. Three success criteria (v3 §4.3):
#   (i)   sign recovery: ΔΔG_RESP < 0 (Magotti SSOT-consistent)
#   (ii)  magnitude convergence: |ΔΔG_RESP − [−1.4, −3.0]| < 3 kcal/mol
#   (iii) σ_btwn stability: ~1× current 6.76 (range ~5–8)
#
# Three informative failure modes (v3 §4.3):
#   - sign recovered + magnitude > 3 kcal/mol off → solvation modeling bottleneck
#   - sign not recovered → broader force-field weakness (dihedral / vdW)
#   - σ_btwn jump > 2× → cohort instability
#
# Stages:
#   A: 8 new Cp4 fresh MD seeds (s17, s29, s47, s61, s73, s79, s97, s131)
#      with params/MTR_gaff2_resp.xml (Khoury OMW RESP-A2 charges, Σq=0)
#      MD-time topology cyclic_ss, dt=2fs, 5ns each.
#   B: 8 Cp4 RESP seed snapshot reextract via phase_beta_repbsa_v2_reextract.
#   C: 8 Cp4 RESP seed MMPBSA postl387_v2 (CUDA mode 4-lane parallel).
#   D: branched_ddg aggregation — RESP Cp4 cohort (n=8) vs existing WT (n=10).
#
# Cohort dir naming: outputs/2QKI_Cp4_resp_calib_s{seed}/ (distinct from
# 2QKI_Cp4_calib_s* legacy cohort).
#
# Compute envelope: 8 seeds × 50min MD = 7h sequential (or ~2h 4-lane GPU
# parallel if VRAM allows; default sequential for safety). + ~30min reextract
# + ~40min MMPBSA 4-lane = ~8-10h total wall.
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh

# CPU affinity 자동 분배 (t1_phase1_5_fresh.sh 와 동일 logic)
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

# v0.8 RESP charges path
RESP_XML="${PROJ}/params/MTR_gaff2_resp.xml"
RESP_HYDROGENS="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_hydrogens.xml"
RESP_MANIFEST_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"  # 기존 baseline 로 _md_input 복제

LAUNCH_TS=$(date +%s)
TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/v08_resp_pilot_${TS_TAG}"
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"

log () {
    printf '[v08_resp_pilot %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "v0.8 Khoury RESP migration pilot START"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
log "logdir=$LOGDIR"
log "RESP XML: $RESP_XML"
log "phys_cores=$PHYS_CORES  logical=$LOGICAL_CORES  affinity=${USE_AFFINITY}"
log "================================================================"

# Pre-flight: verify RESP XML exists + Σq close to zero
if [ ! -f "$RESP_XML" ]; then
    abort "RESP XML not found: $RESP_XML — run scripts/build_khoury_mtr_resp_xml.py first."
fi

$PY -c "
import xml.etree.ElementTree as ET
t = ET.parse('$RESP_XML')
res = next(r for r in t.findall('.//Residue') if r.get('name') == 'MTR')
total = sum(float(a.get('charge')) for a in res.findall('Atom'))
print(f'Σq verification: {total:+.2e}')
assert abs(total) < 1e-5, f'Σq not neutral: {total}'
print('PASS')
" || abort "RESP XML Σq verification failed"

# Seed list (8 new integers, no collision with legacy 2QKI_Cp4_calib_*)
CP4_RESP_SEEDS=(17 29 47 61 73 79 97 131)

# ============================================================================
# Stage A — 8 Cp4 RESP MD seeds (sequential, single-GPU)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage A: 8 Cp4 fresh RESP-MD seeds: ${CP4_RESP_SEEDS[*]}"
log "----------------------------------------------------------------"

for seed_int in "${CP4_RESP_SEEDS[@]}"; do
    seed="s${seed_int}"
    seedir="${PROJ}/outputs/2QKI_Cp4_resp_calib_${seed}"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  $seed: SKIP (DCD exists $(du -h ${seedir}/mdresult/2QKI_Cp4_restrained.dcd | cut -f1))"
        continue
    fi

    log "  $seed: setup dir"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    # _md_input 복제 (PDB 같음, charges 만 다름)
    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && cp -n "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    # RESP XML install (cohort 식별 위해 local copy)
    cp -f "$RESP_XML" "${seedir}/params/MTR_gaff2.xml"  # 이름은 MTR_gaff2 유지 (run_restrained_md.py 의 manifest 의존성)
    cp -n "$RESP_HYDROGENS" "${seedir}/params/"
    cp -n "$RESP_MANIFEST_SRC" "${seedir}/params/"

    # Manifest 의 xml_path 를 local 로 update
    $PY -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
m['charge_source'] = 'Khoury 2014 OMW RESP-A2 (DOI: 10.1021/sb400168u)'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"

    # Verify q_N (Khoury OMW value)
    qN=$(grep -m1 '<Atom name="N" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    log "  $seed: RESP q_N=$qN (Khoury target: −0.280943; amber14SB-patched baseline: −0.4157)"

    log "  $seed: MD launch (5ns, cyclic_ss, RESP charges, seed=${seed_int})"
    conda activate md_simulation
    UPDD_NCAA_AMBER14_PATCH=0 UPDD_MMGBSA_PLATFORM=CUDA \
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
            > "${LOGDIR}/md_resp_${seed}.log" 2>&1
    rc=$?
    conda deactivate

    if [ "$rc" -ne 0 ] || [ ! -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  $seed: MD FAIL rc=$rc (see ${LOGDIR}/md_resp_${seed}.log)"
        continue
    fi
    log "  $seed: MD DONE"
done

log "Stage A complete. ${#CP4_RESP_SEEDS[@]} RESP-MD seeds processed."

# ============================================================================
# Stage B — Snapshot reextract for new RESP seeds
# ============================================================================
log "----------------------------------------------------------------"
log "Stage B: snapshot reextract (8 RESP seeds)"
log "----------------------------------------------------------------"
conda activate qmmm

for seed_int in "${CP4_RESP_SEEDS[@]}"; do
    seed="s${seed_int}"
    SYS="2QKI_Cp4_resp_calib_${seed}"
    target="${PROJ}/outputs/${SYS}/snapshots_n25_postl387_patch_v2"
    if [ -d "$target" ] && [ "$(ls $target/*.pdb 2>/dev/null | wc -l)" -ge 25 ]; then
        log "  reextract ${SYS} SKIP (already $(ls $target/*.pdb | wc -l) PDBs)"
        continue
    fi
    log "  reextract ${SYS}  affinity=[${MMPBSA_AFFINITY:-free}]"
    CUDA_VISIBLE_DEVICES="" $MMPBSA_PREFIX $PY scripts/phase_beta_repbsa_v2_reextract.py \
        --filter "${SYS}" \
        --out-json "${LOGDIR}/reextract_${SYS}.json" \
        >> "${LOGDIR}/reextract.log" 2>&1
done

log "Stage B complete."

# ============================================================================
# Stage C — MMPBSA 4-lane parallel (CUDA mode per Phase 1.5 lesson)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage C: MMPBSA 4-lane CUDA mode"
log "----------------------------------------------------------------"

LANES=4

run_one_pbsa () {
    local label="$1" snap_dir="$2" out_dir="$3"
    if [ -f "${out_dir}/mmpbsa_summary.json" ]; then
        log "  $label: SKIP (summary exists)"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        log "  $label: MISS no pdb"
        return 1
    fi
    log "  $label: START"
    local START=$(date +%s)
    UPDD_MMGBSA_PLATFORM=CUDA $MMPBSA_PREFIX "$PY" -u scripts/run_mmpbsa.py \
        --md_dir "$snap_dir" --outputdir "$out_dir" \
        --ncaa_elem MTR --receptor_chain A --binder_chain B \
        --target_id 2QKI --protocol 1traj \
        > "${LOGDIR}/mmpbsa_${label}.log" 2>&1
    local RC=$?
    local END=$(date +%s)
    log "  $label: DONE rc=$RC elapsed=$((END-START))s"
}

running=0
for seed_int in "${CP4_RESP_SEEDS[@]}"; do
    seed="s${seed_int}"
    SNAP="outputs/2QKI_Cp4_resp_calib_${seed}/snapshots_n25_postl387_patch_v2"
    OUT="outputs/2QKI_Cp4_resp_calib_${seed}/mmpbsa_results_postl387_v2"
    run_one_pbsa "RESP_Cp4_${seed}" "$SNAP" "$OUT" &
    running=$((running + 1))
    if [ "$running" -ge "$LANES" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait

log "Stage C complete."

# ============================================================================
# Stage D — Branched ΔΔG aggregation (RESP Cp4 vs existing WT)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage D: branched_ddg aggregation (RESP Cp4 vs existing WT n=10)"
log "----------------------------------------------------------------"
"$PY" -m utils.branched_ddg \
    --wt "outputs/2QKI_WT_calib_*/mmpbsa_results_postl387_v2" \
    --variant "outputs/2QKI_Cp4_resp_calib_*/mmpbsa_results_postl387_v2" \
    --output "${LOGDIR}/branched_ddg/" \
    >> "$DISPATCH_LOG" 2>&1 && log "Stage D DONE" || log "Stage D WARN non-zero"

# Summary
log "================================================================"
log "v0.8 RESP pilot COMPLETE"
log "FINAL branched_ddg.json: ${LOGDIR}/branched_ddg/branched_ddg.json"
log ""
log "Final state per seed:"
for seed_int in "${CP4_RESP_SEEDS[@]}"; do
    seed="s${seed_int}"
    out="${PROJ}/outputs/2QKI_Cp4_resp_calib_${seed}/mmpbsa_results_postl387_v2"
    if [ -f "${out}/mmpbsa_summary.json" ]; then
        log "  RESP_Cp4_$seed: ✓ done"
    else
        log "  RESP_Cp4_$seed: ✗ failed"
    fi
done
log "================================================================"
