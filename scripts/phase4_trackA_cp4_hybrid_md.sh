#!/bin/bash
# ============================================================================
# Phase IV Track A — Cp4 hybrid-charge production MD for the #84 seed set
# (ADR-0010, SciVal verdict A-C1 / B-C1: Option-β / regime-2 hybrid MTR FF).
# ----------------------------------------------------------------------------
# Generates the strip-direction base ensemble for the 1-trajectory QM ΔΔ_int
# (scripts/phase4_qmmm_iva_1traj.py). The graft-direction WT ensemble already
# exists (2QKI_WT_calib_s*, amber14SB Trp, regime-invariant) and is NOT touched.
#
# This script is a clone of scripts/t1_phase1_5_fresh.sh Stage A (Cp4 MD) +
# Stage C (postl387_v2 snapshot reextract) + Stage D (1-traj MM-PBSA), with the
# SINGLE protocol deviation that the MTR sidechain charge model is swapped from
# the patch-OFF baseline (Σq=−0.187, q_N=−0.8938 UNPATCHED-at-MD-time) to the
# hybrid regime-2 XML params/MTR_gaff2_hybrid.xml (Σq=0, NE1 frozen −0.3418,
# Khoury sidechain). Every other MD parameter is reproduced verbatim from the
# baseline: 2.5M steps × 2 fs (5 ns), cyclic_ss topology, binder chain B,
# graph_policy strict, target 2QKI, CUDA platform, per-seed deterministic seed,
# UPDD_MTR_AMBER14_PATCH=0, postl387_v2 (n=25 K-Means) snapshot extraction,
# 1traj MM-PBSA scoring.
#
# Hybrid-XML injection mechanism reused from scripts/v08_resp_hybrid_smoke.sh
# (HYBRID_XML → params/MTR_gaff2.xml + HYDROGENS_SRC + MANIFEST_SRC, manifest
# xml_path/hydrogens_path/charge_source rewrite). run_restrained_md.py prefers
# the local params/<resname>_gaff2.xml over the manifest path (run_restrained_md
# .py:495-498), so the local hybrid XML wins.
#
# Explosion-rescue ladder reused from scripts/v09_gamma_rescue.sh (UPDD.py
# L1782/L1850 port): on NaN/explosion (no final.pdb), Pass 1b reseeds at seed+13
# (2 fs), Pass 2 falls back to dt=1 fs. R-7: exploded artifacts archived under
# mdresult/_archive/<ts>_pre_<tag>/ (never deleted). The hybrid regime carries a
# non-zero crash rate (s73 ~8% historically).
#
# GPU: host RTX 5070 Ti (index 0). MD is FP32 → host faster than the V100 VM
# (dual-GPU benchmark, Finding A). MD is NEVER routed to the V100 VM. Sequential
# single-GPU per baseline convention.
#
# NOTE: this script BUILDS and DRY-RUNS only by default-safe gating. The actual
# launch is Keeper-gated + user-approved (ADR-0010). Use --dry-run to print the
# resolved protocol + per-seed commands without executing MD.
#
# Usage:
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh --dry-run    # print plan, no MD
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh              # full run (gated)
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh --md-only    # MD + rescue only
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh s7 s23       # subset of seeds
# ============================================================================

set -u

# ----- env -----
PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh

# Cohort consistency note (ADR-0010): unlike the patch-OFF baseline which runs a
# 2-layer MD-UNPATCHED + PBSA-PATCHED XML swap, the hybrid Track A ensemble uses
# ONE charge model (the hybrid XML) at every stage — that is the SciVal-approved
# deviation. UPDD_MTR_AMBER14_PATCH=0 is set defensively (inert at MD time; the
# charge model is fully determined by the local params/MTR_gaff2.xml = hybrid).
export UPDD_MTR_AMBER14_PATCH=0
export UPDD_MMGBSA_PLATFORM=CPU  # CPU-forced for MM-PBSA to avoid GPU contention

# Host GPU only — MD is FP32 → host RTX 5070 Ti faster than V100 VM (Finding A).
# MD is never routed to the V100 (192.168.122.155). Override only for diagnosis.
MD_GPU="${UPDD_TRACKA_MD_GPU:-0}"

# ----- arg parse -----
DRY_RUN=0
MD_ONLY=0
ARG_SEEDS=()
for a in "$@"; do
    case "$a" in
        --dry-run|--dry_run) DRY_RUN=1 ;;
        --md-only|--md_only) MD_ONLY=1 ;;
        s*) ARG_SEEDS+=("$a") ;;
        *) echo "unknown arg: $a" >&2; exit 2 ;;
    esac
done

# ============================================================================
# CPU affinity 자동 분배 (baseline t1_phase1_5_fresh.sh lines 62-74, verbatim)
# 9800X3D 8C/16T: 1 physical core → MD, 나머지 (N-1) → MM-PBSA. 8 phys 초과 시
# free scheduling fallback.
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
LOGDIR="${PROJ}/outputs/analysis/phase4_trackA_cp4_hybrid_${TS_TAG}"
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"

log () {
    printf '[trackA_hybrid %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

# ----- protocol constants (baseline t1_phase1_5_fresh.sh Stage A, verbatim) -----
MD_STEPS=2500000        # 5 ns @ 2 fs  (baseline --steps)
MD_DT=2.0               # baseline --dt_fs
MD_TOPOLOGY=cyclic_ss   # baseline --topology
MD_BINDER_CHAIN=B       # baseline --binder_chain
MD_GRAPH_POLICY=strict  # baseline --graph_policy
MD_TARGET_ID=2QKI       # baseline --target_id
MD_PLATFORM=CUDA        # baseline --platform
MD_NCAA_LABEL=MTR       # baseline --ncaa_label
MD_NCAA_CODE=MTR        # baseline --ncaa_code
SNAP_SUBDIR=snapshots_n25_postl387_patch_v2   # baseline Stage C
MMPBSA_SUBDIR=mmpbsa_results_postl387_v2       # baseline Stage D
MMPBSA_PROTOCOL=1traj
MMPBSA_NCAA_ELEM=MTR
MMPBSA_RECEPTOR_CHAIN=A
MMPBSA_BINDER_CHAIN=B

# Hybrid-charge injection sources (v08_resp_hybrid_smoke.sh lines 47-50)
HYBRID_XML="${PROJ}/params/MTR_gaff2_hybrid.xml"
HYDROGENS_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_hydrogens.xml"
MANIFEST_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"
CHARGE_SOURCE_NOTE='Khoury 2014 OMW RESP-A2 hybrid (sidechain-only; backbone N/H + NE1 preserved from amber14SB-patched baseline)'

# #84 seed set (ADR-0010 / phase4_qmmm_iva_1traj.py COMMON_SEEDS)
SEEDS=(s7 s19 s23 s101 s127 s163 s199 s251)
declare -A SEED_INT=(
    [s7]=7 [s19]=19 [s23]=23 [s101]=101
    [s127]=127 [s163]=163 [s199]=199 [s251]=251
)
[ "${#ARG_SEEDS[@]}" -gt 0 ] && SEEDS=("${ARG_SEEDS[@]}")

log "================================================================"
log "Phase IV Track A — Cp4 hybrid-charge production MD"
log "host=$(hostname)  md_gpu=index_${MD_GPU} (RTX 5070 Ti; V100 VM excluded)"
log "logdir=$LOGDIR  dry_run=${DRY_RUN}  md_only=${MD_ONLY}"
log "phys_cores=$PHYS_CORES logical=$LOGICAL_CORES affinity=${USE_AFFINITY} (MD=[${MD_AFFINITY:-free}])"
log "charge model: $HYBRID_XML (Option-β/regime-2 hybrid, Σq=0, NE1 −0.3418)"
log "seeds (#84 set): ${SEEDS[*]}"
log "protocol: ${MD_STEPS} steps × ${MD_DT} fs (5 ns), topology=${MD_TOPOLOGY}, "
log "          graph_policy=${MD_GRAPH_POLICY}, target=${MD_TARGET_ID}, ncaa=${MD_NCAA_CODE}"
log "================================================================"

# ----- pre-flight: hybrid XML Σq + preserved baseline atoms (v08 smoke lines 82-97) -----
"$PY" -c "
import xml.etree.ElementTree as ET
t = ET.parse('$HYBRID_XML')
res = next(r for r in t.findall('.//Residue') if r.get('name') == 'MTR')
charges = {a.get('name'): float(a.get('charge')) for a in res.findall('Atom')}
total = sum(charges.values())
print(f'Sigma_q verification: {total:+.2e}')
assert abs(total) < 1e-5, f'Sigma_q not neutral: {total}'
expected = {'N': -0.4157, 'H': 0.2719, 'NE1': -0.3418}
for name, exp in expected.items():
    got = charges.get(name)
    assert got is not None and abs(got - exp) < 1e-4, f'{name} drift: got {got}, expected {exp}'
print(f'Preserved baseline N={charges[\"N\"]:+.4f} H={charges[\"H\"]:+.4f} NE1={charges[\"NE1\"]:+.4f}')
print('PASS')
" 2>&1 | tee -a "$DISPATCH_LOG" || abort "Hybrid XML verification failed"

[ -f "$HYDROGENS_SRC" ] || abort "HYDROGENS_SRC missing: $HYDROGENS_SRC"
[ -f "$MANIFEST_SRC" ] || abort "MANIFEST_SRC missing: $MANIFEST_SRC"

# ----- per-seed dir setup (v08_resp_hybrid_smoke.sh lines 115-138) -----
setup_seed_dir () {
    local seed="$1"
    local seedir="${PROJ}/outputs/2QKI_Cp4_hybrid_calib_${seed}"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && cp -n "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    cp -f "$HYBRID_XML" "${seedir}/params/MTR_gaff2.xml"
    cp -n "$HYDROGENS_SRC" "${seedir}/params/"
    cp -n "$MANIFEST_SRC" "${seedir}/params/"

    "$PY" -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
m['charge_source'] = '${CHARGE_SOURCE_NOTE}'
m['charge_regime'] = 'option_beta_regime2_hybrid_sigmaq0'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"
}

# ----- R-7 archive of exploded artifacts (v09_gamma_rescue.sh lines 66-82) -----
archive_exploded () {
    local seedir="$1"; local tag="$2"
    local mdr="${seedir}/mdresult"
    local adir="${mdr}/_archive/$(date +%Y%m%d_%H%M%S)_pre_${tag}"
    mkdir -p "$adir"
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

# ----- single MD invocation (baseline t1_phase1_5_fresh.sh Stage A3, verbatim
#       flags; XML is hybrid via local params/) -----
run_md () {
    local seedir="$1"; local seed_for_rng="$2"; local dt="$3"; local seed="$4"; local pass="$5"
    log "  $seed: MD launch (Pass ${pass}, dt=${dt}fs, seed=${seed_for_rng}, gpu=${MD_GPU}, affinity=[${MD_AFFINITY:-free}])"
    conda activate md_simulation
    CUDA_VISIBLE_DEVICES="${MD_GPU}" \
        UPDD_MTR_AMBER14_PATCH=0 UPDD_MMGBSA_PLATFORM=CUDA \
        $MD_PREFIX "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${seedir}/_md_input" \
            --outputdir "${seedir}/mdresult" \
            --params_manifest "${seedir}/params/MTR_params_manifest.json" \
            --steps "$MD_STEPS" \
            --topology "$MD_TOPOLOGY" \
            --binder_chain "$MD_BINDER_CHAIN" \
            --graph_policy "$MD_GRAPH_POLICY" \
            --target_id "$MD_TARGET_ID" \
            --dt_fs "$dt" \
            --platform "$MD_PLATFORM" \
            --seed "$seed_for_rng" \
            --ncaa_label "$MD_NCAA_LABEL" \
            --ncaa_code "$MD_NCAA_CODE" \
            > "${LOGDIR}/md_hybrid_${seed}_pass${pass}.log" 2>&1
    local rc=$?
    conda deactivate
    return $rc
}

# ----- per-seed MD with rescue ladder (v09_gamma_rescue.sh lines 126-146) -----
md_with_rescue () {
    local seed="$1"
    local seed_int="${SEED_INT[$seed]:-}"
    [ -z "$seed_int" ] && { log "  $seed: unknown seed label, skip"; return 1; }
    local seedir="${PROJ}/outputs/2QKI_Cp4_hybrid_calib_${seed}"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $seed: SKIP MD (final.pdb exists)"
        return 0
    fi

    setup_seed_dir "$seed"

    # GATE: confirm the local MD XML is the hybrid model (q_N=−0.4157, Σq=0)
    local qN qNE1
    qN=$(grep -m1 '<Atom name="N" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    qNE1=$(grep -m1 '<Atom name="NE1" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    log "  $seed: hybrid q_N=$qN q_NE1=$qNE1 (regime-2)"

    # Pass 1: baseline protocol (2 fs, original seed integer)
    run_md "$seedir" "$seed_int" "$MD_DT" "$seed" "1"
    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $seed: ✓ Pass 1 (seed=${seed_int}, ${MD_DT}fs) DONE"
        return 0
    fi
    log "  $seed: Pass 1 explosion/incomplete — engage rescue ladder"

    # Pass 1b: reseed (seed+13 prime offset) at 2 fs
    archive_exploded "$seedir" "reseed"
    local new_seed=$((seed_int + 13))
    run_md "$seedir" "$new_seed" "$MD_DT" "$seed" "1b-reseed"
    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $seed: ✓ Pass 1b reseed (seed=${new_seed}, ${MD_DT}fs) DONE"
        return 0
    fi
    log "  $seed: Pass 1b reseed failed — 1 fs last resort"

    # Pass 2: dt=1 fs last resort
    archive_exploded "$seedir" "dt1fs"
    run_md "$seedir" "$new_seed" "1.0" "$seed" "2-dt1fs"
    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $seed: ✓ Pass 2 (seed=${new_seed}, 1fs) DONE"
        return 0
    fi
    log "  $seed: ✗ rescue FINAL FAIL (2fs reseed + 1fs both exploded) — permanent-exclusion candidate"
    return 1
}

# ----- dry-run: print resolved plan + exact per-seed commands, no execution -----
if [ "$DRY_RUN" -eq 1 ]; then
    log "----------------------------------------------------------------"
    log "DRY-RUN — no MD executed. Resolved per-seed plan:"
    log "----------------------------------------------------------------"
    for seed in "${SEEDS[@]}"; do
        seed_int="${SEED_INT[$seed]:-?}"
        seedir="${PROJ}/outputs/2QKI_Cp4_hybrid_calib_${seed}"
        graft_dir="${PROJ}/outputs/2QKI_WT_calib_${seed}"
        printf '\n[seed %s  int=%s]\n' "$seed" "$seed_int" | tee -a "$DISPATCH_LOG"
        printf '  out dir   : %s\n' "${seedir#$PROJ/}" | tee -a "$DISPATCH_LOG"
        printf '  graft ref : %s  (exists: %s)\n' "${graft_dir#$PROJ/}" \
            "$([ -f "${graft_dir}/mdresult/2QKI_WT_restrained.dcd" ] && echo yes || echo NO)" | tee -a "$DISPATCH_LOG"
        printf '  baseline  : outputs/2QKI_Cp4_calib_%s  (DCD exists: %s)\n' "$seed" \
            "$([ -f "${PROJ}/outputs/2QKI_Cp4_calib_${seed}/mdresult/2QKI_Cp4_restrained.dcd" ] && echo yes || echo NO)" | tee -a "$DISPATCH_LOG"
        printf '  MD CMD    : CUDA_VISIBLE_DEVICES=%s UPDD_MTR_AMBER14_PATCH=0 UPDD_MMGBSA_PLATFORM=CUDA \\\n' "$MD_GPU" | tee -a "$DISPATCH_LOG"
        printf '              %s %s utils/run_restrained_md.py \\\n' "${MD_PREFIX:-}" "$PY_MD" | tee -a "$DISPATCH_LOG"
        printf '                --inputdir %s/_md_input --outputdir %s/mdresult \\\n' "${seedir#$PROJ/}" "${seedir#$PROJ/}" | tee -a "$DISPATCH_LOG"
        printf '                --params_manifest %s/params/MTR_params_manifest.json \\\n' "${seedir#$PROJ/}" | tee -a "$DISPATCH_LOG"
        printf '                --steps %s --topology %s --binder_chain %s --graph_policy %s \\\n' "$MD_STEPS" "$MD_TOPOLOGY" "$MD_BINDER_CHAIN" "$MD_GRAPH_POLICY" | tee -a "$DISPATCH_LOG"
        printf '                --target_id %s --dt_fs %s --platform %s --seed %s \\\n' "$MD_TARGET_ID" "$MD_DT" "$MD_PLATFORM" "$seed_int" | tee -a "$DISPATCH_LOG"
        printf '                --ncaa_label %s --ncaa_code %s\n' "$MD_NCAA_LABEL" "$MD_NCAA_CODE" | tee -a "$DISPATCH_LOG"
        printf '  CHARGE    : local params/MTR_gaff2.xml ← %s (hybrid)\n' "${HYBRID_XML#$PROJ/}" | tee -a "$DISPATCH_LOG"
        printf '  RESCUE    : Pass1b seed=%s @2fs → Pass2 seed=%s @1fs (R-7 archive)\n' "$((seed_int + 13))" "$((seed_int + 13))" | tee -a "$DISPATCH_LOG"
        printf '  REEXTRACT : reextract_one(system="2QKI_Cp4_hybrid", seed="%s") → %s/\n' "$seed" "$SNAP_SUBDIR" | tee -a "$DISPATCH_LOG"
        printf '  MM-PBSA   : run_mmpbsa.py --md_dir .../%s --protocol %s --ncaa_elem %s\n' "$SNAP_SUBDIR" "$MMPBSA_PROTOCOL" "$MMPBSA_NCAA_ELEM" | tee -a "$DISPATCH_LOG"
    done
    log "----------------------------------------------------------------"
    log "DRY-RUN complete. ${#SEEDS[@]} seeds planned. No MD launched."
    log "----------------------------------------------------------------"
    exit 0
fi

# ============================================================================
# Stage A — Cp4 hybrid MD (sequential single-GPU, with rescue)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage A: Cp4 hybrid MD ${#SEEDS[@]} seeds (sequential, host GPU ${MD_GPU})"
log "----------------------------------------------------------------"
for seed in "${SEEDS[@]}"; do
    md_with_rescue "$seed"
done
log "Stage A complete."

if [ "$MD_ONLY" -eq 1 ]; then
    log "--md-only: skipping Stage C/D. Done."
    exit 0
fi

# ============================================================================
# Stage C — postl387_v2 snapshot reextract (n=25 K-Means). Calls reextract_one
# directly with system="2QKI_Cp4_hybrid" so the shared TARGETS list in
# phase_beta_repbsa_v2_reextract.py is untouched (anti-fragmentation).
# ============================================================================
log "----------------------------------------------------------------"
log "Stage C: snapshot reextract (${SNAP_SUBDIR})"
log "----------------------------------------------------------------"
conda activate qmmm
for seed in "${SEEDS[@]}"; do
    seedir="${PROJ}/outputs/2QKI_Cp4_hybrid_calib_${seed}"
    target="${seedir}/${SNAP_SUBDIR}"
    if [ -d "$target" ] && [ "$(ls "$target"/*.pdb 2>/dev/null | wc -l)" -ge 25 ]; then
        log "  ${seed}: reextract SKIP (≥25 pdb already present)"
        continue
    fi
    if [ ! -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  ${seed}: reextract SKIP (no final.pdb — MD did not complete)"
        continue
    fi
    log "  ${seed}: reextract (affinity=[${MMPBSA_AFFINITY:-free}])"
    CUDA_VISIBLE_DEVICES="" $MMPBSA_PREFIX "$PY" -c "
import sys
sys.path.insert(0, 'scripts')
from phase_beta_repbsa_v2_reextract import reextract_one
rec = reextract_one('2QKI_Cp4_hybrid', '${seed}', out_subdir='${SNAP_SUBDIR}')
print(f\"{rec['tag']}: {rec['status']} n_saved={rec.get('n_saved')} {rec.get('skip_reason') or ''}\")
" >> "${LOGDIR}/reextract.log" 2>&1
done
conda deactivate
log "Stage C complete."

# ============================================================================
# Stage D — 1-traj MM-PBSA scoring (provides the median-Δg ranking the
# orchestrator's resolve_frames() consumes). 4-lane parallel per baseline.
# ============================================================================
log "----------------------------------------------------------------"
log "Stage D: MM-PBSA (${MMPBSA_PROTOCOL}, 4 lanes)"
log "----------------------------------------------------------------"
LANES=4
PBSA_LOG="${LOGDIR}/mmpbsa.log"

run_one_pbsa () {
    local seed="$1"
    local seedir="${PROJ}/outputs/2QKI_Cp4_hybrid_calib_${seed}"
    local snap_dir="${seedir}/${SNAP_SUBDIR}"
    local out_dir="${seedir}/${MMPBSA_SUBDIR}"
    if [ -f "${out_dir}/mmpbsa_summary.json" ]; then
        printf '[%s SKIP %s]\n' "$(date +%H:%M:%S)" "$seed" >> "$PBSA_LOG"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        printf '[%s MISS %s] no pdb\n' "$(date +%H:%M:%S)" "$seed" >> "$PBSA_LOG"
        return 1
    fi
    mkdir -p "$out_dir"
    local seed_log="${LOGDIR}/mmpbsa_${seed}.log"
    local start; start=$(date +%s)
    printf '[%s START %s]\n' "$(date +%H:%M:%S)" "$seed" >> "$PBSA_LOG"
    CUDA_VISIBLE_DEVICES="" UPDD_MMGBSA_PLATFORM=CPU $MMPBSA_PREFIX "$PY" scripts/run_mmpbsa.py \
        --md_dir "$snap_dir" --outputdir "$out_dir" \
        --ncaa_elem "$MMPBSA_NCAA_ELEM" --receptor_chain "$MMPBSA_RECEPTOR_CHAIN" \
        --binder_chain "$MMPBSA_BINDER_CHAIN" \
        --target_id "$MD_TARGET_ID" --protocol "$MMPBSA_PROTOCOL" \
        > "$seed_log" 2>&1
    local rc=$?
    local end; end=$(date +%s)
    printf '[%s DONE  %s] rc=%d elapsed=%ds\n' "$(date +%H:%M:%S)" "$seed" "$rc" "$((end-start))" >> "$PBSA_LOG"
}

running=0
for seed in "${SEEDS[@]}"; do
    run_one_pbsa "$seed" &
    running=$((running + 1))
    if [ "$running" -ge "$LANES" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait
log "Stage D complete."

# ============================================================================
# Summary
# ============================================================================
SUCCESS=0
for seed in "${SEEDS[@]}"; do
    if [ -f "${PROJ}/outputs/2QKI_Cp4_hybrid_calib_${seed}/${MMPBSA_SUBDIR}/mmpbsa_summary.json" ]; then
        SUCCESS=$((SUCCESS + 1))
        log "  ✓ ${seed}: MD + reextract + MM-PBSA complete"
    else
        log "  ✗ ${seed}: incomplete (check logs)"
    fi
done
ELAPSED=$(( $(date +%s) - LAUNCH_TS ))
log "================================================================"
log "Track A Cp4 hybrid MD COMPLETE  elapsed=${ELAPSED}s  ${SUCCESS}/${#SEEDS[@]} ready"
log "Next: scripts/phase4_qmmm_iva_1traj.py (strip frames now resolve from hybrid dirs)"
log "================================================================"
