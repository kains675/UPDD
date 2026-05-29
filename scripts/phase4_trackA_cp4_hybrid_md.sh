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
# 2-LANE DISPATCH (added): the seed batch is split across two concurrent lanes —
# the host RTX 5070 Ti (local) and the VM Tesla V100 (san@192.168.122.155, VFIO
# passthrough, host-invisible). Lanes run DISJOINT seeds → disjoint output dirs
# (2QKI_Cp4_hybrid_calib_s{N}) → zero contention, no locks. Within a lane seeds
# run sequentially; the two lanes run in parallel (wait for both). Default split
# for the #84 batch: host = s7 s19 s23 s101 s127 (5), vm = s163 s199 s251 (3) —
# balanced for the ~1.3× host:V100 MD ratio (target makespan ~4h vs ~6.7h serial).
#
# The seed split, system name, variant tag, and charge XML are all parameters
# (--host-seeds / --vm-seeds / --system / --variant / --charge-xml) so a future
# batch can pass any disjoint split. Disjointness is asserted at startup (abort
# on overlap); the union is reported.
#
# V100 lane transport reuses utils/dispatch.py::VMExecutor (rsync + SSH + retry)
# — the SAME transport the QM orchestrator uses; no new transport is invented.
# Per VM seed: (a) VMExecutor.sync_to_vm pushes the seed prerequisites + the
# hybrid params XML to the VM, (b) VMExecutor.execute SSH-runs THIS SAME script
# on the VM in --vm-worker mode (so the identical setup/MD/rescue/Stage C/Stage D
# logic runs there — no MD logic duplicated), pinned to the V100 with absolute
# conda paths (non-interactive SSH lacks the conda PATH), (c) VMExecutor.sync_
# from_vm rsyncs the result dir (snapshot subdir + mmpbsa_summary.json) back.
# The VM has full AmberTools (MMPBSA.py/ante-MMPBSA.py/tleap/sander in env qmmm,
# probed read-only) so Stage D MM-PBSA runs ON the VM — no DCD sync-back-for-host.
#
# NOTE: this script BUILDS and DRY-RUNS only by default-safe gating. The actual
# launch is Keeper-gated + user-approved (ADR-0010). Use --dry-run to print the
# resolved protocol + per-lane seed assignment + per-seed host commands + the VM
# rsync/SSH commands without executing MD.
#
# Usage:
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh --dry-run    # print 2-lane plan, no MD
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh              # full 2-lane run (gated)
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh --md-only    # MD + rescue only (both lanes)
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh \
#        --host-seeds "s7 s19 s23 s101 s127" --vm-seeds "s163 s199 s251"
#   bash scripts/phase4_trackA_cp4_hybrid_md.sh --vm-seeds "" s7 s23   # host-only subset
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

# Host-lane GPU index (RTX 5070 Ti = 0). Override only for diagnosis.
MD_GPU="${UPDD_TRACKA_MD_GPU:-0}"

# ----- V100 VM lane config (transport via utils/dispatch.py::VMExecutor) -----
# These mirror dispatch.py's UPDD_VM_* env defaults so the bash lane and the
# Python VMExecutor agree on the same target. The VM's V100 is index 0 inside
# the VM (probed: single Tesla V100-PCIE-32GB). conda is NOT on the non-login
# SSH PATH → absolute path required (this bit us before).
VM_SSH_TARGET="${UPDD_VM_SSH_TARGET:-san@192.168.122.155}"
VM_PROJECT_ROOT="${UPDD_VM_PROJECT_ROOT:-/home/san/UPDD_proj}"
VM_CONDA="${UPDD_VM_CONDA:-/home/san/miniconda3/bin/conda}"
VM_MD_GPU="${UPDD_VM_MD_GPU:-0}"   # V100 index inside the VM

# ----- arg parse -----
DRY_RUN=0
MD_ONLY=0
VM_WORKER=0                # internal: set when THIS script runs on the VM for one lane
ARG_SEEDS=()               # positional seed subset (host-only convenience; restricts the batch)
HOST_SEEDS_ARG=""          # --host-seeds "s7 s19 ..."  (unset → default split)
VM_SEEDS_ARG=""            # --vm-seeds "s163 ..."      (unset → default split; "" → host-only)
VM_SEEDS_SET=0             # 1 once --vm-seeds is seen (so an explicit "" means host-only)
# Reusability parameters (a future batch overrides system/variant/charge model)
SYSTEM_BASE="${UPDD_TRACKA_SYSTEM:-2QKI_Cp4}"        # base system (matches Cp4 ref + reextract system stem)
VARIANT_TAG="${UPDD_TRACKA_VARIANT:-hybrid}"         # variant suffix → 2QKI_Cp4_hybrid_calib_s{N}
CHARGE_XML_ARG=""                                    # --charge-xml PATH (default params/MTR_gaff2_hybrid.xml below)
while [ "$#" -gt 0 ]; do
    case "$1" in
        --dry-run|--dry_run) DRY_RUN=1 ;;
        --md-only|--md_only) MD_ONLY=1 ;;
        --vm-worker|--vm_worker) VM_WORKER=1 ;;
        --host-seeds|--host_seeds) HOST_SEEDS_ARG="${2:-}"; shift ;;
        --vm-seeds|--vm_seeds)     VM_SEEDS_ARG="${2:-}"; VM_SEEDS_SET=1; shift ;;
        --system) SYSTEM_BASE="${2:-}"; shift ;;
        --variant) VARIANT_TAG="${2:-}"; shift ;;
        --charge-xml|--charge_xml) CHARGE_XML_ARG="${2:-}"; shift ;;
        s*) ARG_SEEDS+=("$1") ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
    shift
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

# Derived system stems (parameterized for reusability)
#   SYSTEM_BASE   = "2QKI_Cp4"            (ref-dir + reextract system stem)
#   VARIANT_TAG   = "hybrid"             → SYSTEM_TAG = "2QKI_Cp4_hybrid"
#   per-seed out  = outputs/2QKI_Cp4_hybrid_calib_s{N}
SYSTEM_TAG="${SYSTEM_BASE}_${VARIANT_TAG}"
REEXTRACT_SYSTEM="$SYSTEM_TAG"          # reextract_one(system=...) stem
REF_SEED="${UPDD_TRACKA_REF_SEED:-s7}"  # baseline seed dir supplying _md_input + param sources

# Hybrid-charge injection sources (v08_resp_hybrid_smoke.sh lines 47-50)
HYBRID_XML="${CHARGE_XML_ARG:-${PROJ}/params/MTR_gaff2_hybrid.xml}"
HYDROGENS_SRC="${PROJ}/outputs/${SYSTEM_BASE}_calib_${REF_SEED}/params/MTR_hydrogens.xml"
MANIFEST_SRC="${PROJ}/outputs/${SYSTEM_BASE}_calib_${REF_SEED}/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/${SYSTEM_BASE}_calib_${REF_SEED}"
CHARGE_SOURCE_NOTE='Khoury 2014 OMW RESP-A2 hybrid (sidechain-only; backbone N/H + NE1 preserved from amber14SB-patched baseline)'

# #84 seed set (ADR-0010 / phase4_qmmm_iva_1traj.py COMMON_SEEDS)
declare -A SEED_INT=(
    [s7]=7 [s19]=19 [s23]=23 [s101]=101
    [s127]=127 [s163]=163 [s199]=199 [s251]=251
)

# ----- 2-lane seed assignment (default split balanced for ~1.3× host:V100) -----
# Default for the #84 batch: host=5 seeds, vm=3 seeds. Overridable via
# --host-seeds / --vm-seeds. A positional seed subset (ARG_SEEDS) restricts the
# host lane and disables the VM lane (host-only convenience path).
DEFAULT_HOST_SEEDS="s7 s19 s23 s101 s127"
DEFAULT_VM_SEEDS="s163 s199 s251"

if [ "${#ARG_SEEDS[@]}" -gt 0 ]; then
    # positional subset → host-only (no VM lane), VM seeds empty
    HOST_SEEDS=("${ARG_SEEDS[@]}")
    VM_SEEDS=()
else
    # shellcheck disable=SC2206
    HOST_SEEDS=(${HOST_SEEDS_ARG:-$DEFAULT_HOST_SEEDS})
    if [ "$VM_SEEDS_SET" -eq 1 ]; then
        # explicit --vm-seeds (may be "" to force host-only)
        # shellcheck disable=SC2206
        VM_SEEDS=(${VM_SEEDS_ARG})
    else
        # shellcheck disable=SC2206
        VM_SEEDS=(${DEFAULT_VM_SEEDS})
    fi
fi
# When this script runs as the VM worker it owns only the VM seeds, executed
# LOCALLY on the V100 (no nested VM lane).
if [ "$VM_WORKER" -eq 1 ]; then
    HOST_SEEDS=("${VM_SEEDS[@]}")
    VM_SEEDS=()
    MD_GPU="$VM_MD_GPU"
fi
# Full batch = union of the two lanes (for reporting + completeness assert)
ALL_SEEDS=("${HOST_SEEDS[@]}" "${VM_SEEDS[@]}")

# ----- disjointness + known-seed asserts (abort on overlap) -----
assert_lane_seeds () {
    local -A seen=()
    local s overlap=""
    for s in "${HOST_SEEDS[@]}"; do
        [ -z "$s" ] && continue
        [ -n "${SEED_INT[$s]:-}" ] || abort "unknown seed '$s' (not in SEED_INT registry)"
        seen[$s]=host
    done
    for s in "${VM_SEEDS[@]}"; do
        [ -z "$s" ] && continue
        [ -n "${SEED_INT[$s]:-}" ] || abort "unknown seed '$s' (not in SEED_INT registry)"
        if [ "${seen[$s]:-}" = "host" ]; then
            overlap="${overlap} $s"
        fi
        seen[$s]=vm
    done
    if [ -n "$overlap" ]; then
        abort "host-lane and vm-lane seeds OVERLAP →${overlap} (lanes MUST be disjoint to avoid output-dir contention)"
    fi
    if [ "${#HOST_SEEDS[@]}" -eq 0 ] && [ "${#VM_SEEDS[@]}" -eq 0 ]; then
        abort "no seeds assigned to either lane"
    fi
}

assert_lane_seeds

log "================================================================"
log "Phase IV Track A — ${SYSTEM_TAG} production MD (2-lane dispatch)"
log "host=$(hostname)  vm_worker=${VM_WORKER}"
log "logdir=$LOGDIR  dry_run=${DRY_RUN}  md_only=${MD_ONLY}"
log "phys_cores=$PHYS_CORES logical=$LOGICAL_CORES affinity=${USE_AFFINITY} (MD=[${MD_AFFINITY:-free}])"
log "charge model: $HYBRID_XML (Option-β/regime-2 hybrid, Σq=0, NE1 −0.3418)"
log "system_tag=${SYSTEM_TAG}  ref_seed=${REF_SEED}  reextract_system=${REEXTRACT_SYSTEM}"
if [ "$VM_WORKER" -eq 1 ]; then
    log "LANE=vm-worker (running ON the V100 VM)  local_gpu=index_${MD_GPU} (Tesla V100)"
    log "  seeds (this worker): ${HOST_SEEDS[*]:-none}"
else
    log "LANE SPLIT (disjoint, union = full batch):"
    log "  host lane (RTX 5070 Ti idx ${MD_GPU}): ${HOST_SEEDS[*]:-none}  (n=${#HOST_SEEDS[@]})"
    log "  vm   lane (V100 ${VM_SSH_TARGET}): ${VM_SEEDS[*]:-none}  (n=${#VM_SEEDS[@]})"
    log "  full batch (union): ${ALL_SEEDS[*]:-none}  (n=${#ALL_SEEDS[@]})"
fi
log "protocol: ${MD_STEPS} steps × ${MD_DT} fs (5 ns), topology=${MD_TOPOLOGY}, "
log "          graph_policy=${MD_GRAPH_POLICY}, target=${MD_TARGET_ID}, ncaa=${MD_NCAA_CODE}"
log "================================================================"

# ----- pre-flight: hybrid XML Σq + preserved baseline atoms (v08 smoke lines 82-97) -----
# Reads the PROJECT-level hybrid XML + REF_SEED baseline param sources, which live
# OUTSIDE the per-seed dir the VM transport pushes. In --vm-worker mode the worker
# trusts the already-staged-and-pushed per-seed params/ (host ran this pre-flight +
# setup_seed_dir BEFORE pushing), so the worker skips these host-only source reads;
# its only charge gate is the in-worker grep-gate on the delivered per-seed
# params/MTR_gaff2.xml (md_with_rescue, below).
if [ "$VM_WORKER" -eq 0 ]; then
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
else
    log "vm-worker: skip host-only hybrid-XML pre-flight + REF_SEED source checks"
    log "           (trusting the pushed per-seed params/; charge gate = in-worker grep)"
fi

# ----- per-seed dir setup (v08_resp_hybrid_smoke.sh lines 115-138) -----
# On the host this stages the per-seed dir from PROJECT-level + REF_SEED sources
# (hybrid XML → params/MTR_gaff2.xml, hydrogens/manifest copy, manifest rewrite),
# fully populating _md_input/ + params/ BEFORE the VM push. On the VM worker those
# source paths are NOT transport-delivered, so the worker skips the source-copies
# and uses the already-pushed per-seed layout verbatim (only ensuring the dirs
# exist for the MD/reextract/MM-PBSA outputs).
setup_seed_dir () {
    local seed="$1"
    local seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    if [ "$VM_WORKER" -eq 1 ]; then
        # Per-seed _md_input/ + params/{MTR_gaff2.xml,MTR_hydrogens.xml,
        # MTR_params_manifest.json} were staged by the host and pushed by
        # sync_to_vm; do NOT re-read un-pushed host/REF_SEED sources here.
        return 0
    fi

    for f in "${SYSTEM_BASE}.pdb" "${SYSTEM_BASE}_renum.pdb"; do
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
    local seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"

    if [ -f "${seedir}/mdresult/${SYSTEM_BASE}_final.pdb" ]; then
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
    if [ -f "${seedir}/mdresult/${SYSTEM_BASE}_final.pdb" ]; then
        log "  $seed: ✓ Pass 1 (seed=${seed_int}, ${MD_DT}fs) DONE"
        return 0
    fi
    log "  $seed: Pass 1 explosion/incomplete — engage rescue ladder"

    # Pass 1b: reseed (seed+13 prime offset) at 2 fs
    archive_exploded "$seedir" "reseed"
    local new_seed=$((seed_int + 13))
    run_md "$seedir" "$new_seed" "$MD_DT" "$seed" "1b-reseed"
    if [ -f "${seedir}/mdresult/${SYSTEM_BASE}_final.pdb" ]; then
        log "  $seed: ✓ Pass 1b reseed (seed=${new_seed}, ${MD_DT}fs) DONE"
        return 0
    fi
    log "  $seed: Pass 1b reseed failed — 1 fs last resort"

    # Pass 2: dt=1 fs last resort
    archive_exploded "$seedir" "dt1fs"
    run_md "$seedir" "$new_seed" "1.0" "$seed" "2-dt1fs"
    if [ -f "${seedir}/mdresult/${SYSTEM_BASE}_final.pdb" ]; then
        log "  $seed: ✓ Pass 2 (seed=${new_seed}, 1fs) DONE"
        return 0
    fi
    log "  $seed: ✗ rescue FINAL FAIL (2fs reseed + 1fs both exploded) — permanent-exclusion candidate"
    return 1
}

# ============================================================================
# V100 VM lane — transport via utils/dispatch.py::VMExecutor (NOT a new transport)
# ----------------------------------------------------------------------------
# Per VM seed the lane does three VMExecutor calls:
#   (1) sync_to_vm  — push the seed prerequisites + the fresh hybrid XML
#   (2) execute     — SSH-run THIS script in --vm-worker mode on the V100
#                     (identical setup/MD/rescue/Stage C/Stage D logic; the VM
#                      has full AmberTools so Stage D runs there)
#   (3) sync_from_vm— pull the result dir (snapshot subdir + mmpbsa subdir);
#                     the heavy DCD stays on the VM (orchestrator never needs it)
# The remote --vm-worker command pins CUDA_VISIBLE_DEVICES to the V100 and uses
# the absolute conda path (non-interactive SSH lacks the conda PATH).
# ============================================================================

# Remote conda env to run the --vm-worker pipeline in. The worker shells into
# md_simulation/qmmm internally per stage (conda activate inside the script), so
# the OUTER env only needs to provide bash + the project's python entrypoints;
# qmmm carries AmberTools (Stage D) and the openmm import is via md_simulation.
REMOTE_RUN_ENV="${UPDD_VM_RUN_ENV:-qmmm}"

# (0) one-time VM bootstrap — push THIS launcher to the VM scripts/ dir.
# VMExecutor.execute SSH-references the VM-side path `scripts/<this script>`
# (dispatch.py builds `ssh <target> "cd <root> && ... bash scripts/<name> ..."`),
# so the launcher file must physically exist on the VM. The per-seed sync_to_vm
# pushes ONLY outputs/<tag>/, never scripts/; this closes that gap. With the
# --vm-worker Option-B guards (host-only pre-flight + setup_seed_dir source-copies
# skipped), the launcher is the ONLY un-per-seed prerequisite the worker needs:
# the byte-identical run_restrained_md.py / phase_beta_repbsa_v2_reextract.py /
# run_mmpbsa.py already live on the rsync-deployed VM project tree.
vm_bootstrap_script () {
    local script_path; script_path=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
    local remote_scripts="${VM_PROJECT_ROOT}/scripts"
    "$PY" - "$script_path" "$remote_scripts" <<'PYEOF' 2>&1 | tee -a "$DISPATCH_LOG"
import os, sys
from utils.dispatch import VMExecutor
local_script, remote_scripts = sys.argv[1], sys.argv[2]
vm = VMExecutor()
vm.execute(f"mkdir -p {remote_scripts}", retry=2, timeout=30)
# Push the single launcher file (a file, not a dir) into the VM scripts/ dir.
import subprocess
target = f"{vm.ssh_target}:{remote_scripts}/{os.path.basename(local_script)}"
r = subprocess.run(
    ["rsync", "-az",
     "-e", "ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR",
     local_script, target],
    capture_output=True, text=True, timeout=120,
)
ok = r.returncode == 0
print(f"  [VM bootstrap] launcher -> {target}: {'OK' if ok else 'FAIL'}")
if not ok:
    sys.stderr.write((r.stderr or "")[:300])
sys.exit(0 if ok else 1)
PYEOF
    return "${PIPESTATUS[0]}"
}

# Remote --vm-worker command string (also printed verbatim by --dry-run).
vm_worker_cmd () {
    local seed="$1"
    local maybe_md_only=""
    [ "$MD_ONLY" -eq 1 ] && maybe_md_only="--md-only"
    printf 'CUDA_VISIBLE_DEVICES=%s UPDD_MD_CUDA_DEVICE=%s %s run -n %s bash scripts/%s --vm-worker %s --system %s --variant %s --vm-seeds %s' \
        "$VM_MD_GPU" "$VM_MD_GPU" "$VM_CONDA" "$REMOTE_RUN_ENV" "$(basename "$0")" \
        "$maybe_md_only" "$SYSTEM_BASE" "$VARIANT_TAG" "$(printf '%q' "$seed")"
}

# (1) push prerequisites + fresh hybrid XML via VMExecutor.sync_to_vm
vm_push_seed () {
    local seed="$1"
    local seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"
    # Stage the inputs locally first (reuse setup_seed_dir so the VM gets the
    # SAME prepared layout: _md_input PDBs + params/ with the hybrid XML injected
    # + manifest rewritten). Then sync only that small dir (no DCD yet).
    setup_seed_dir "$seed"
    "$PY" - "$seedir" "${VM_PROJECT_ROOT}/outputs/${SYSTEM_TAG}_calib_${seed}" <<'PYEOF' 2>&1 | tee -a "$DISPATCH_LOG"
import sys
from utils.dispatch import VMExecutor
local_dir, remote_dir = sys.argv[1], sys.argv[2]
vm = VMExecutor()
# Override default excludes: we WANT _md_input + params (small) on the VM; only
# guard against any accidental large artifacts (DCD/NC/checkpoints).
ok = vm.sync_to_vm(
    local_dir, remote_path=remote_dir,
    excludes=["*.dcd", "*.nc", "*.h5", "*.chk", "__pycache__/", "_archive/", "mdresult/"],
    timeout=600,
)
print(f"  [VM push] {local_dir} -> {remote_dir}: {'OK' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
PYEOF
    return "${PIPESTATUS[0]}"
}

# (2) SSH-run the --vm-worker pipeline on the V100 via VMExecutor.execute
vm_run_seed () {
    local seed="$1"
    "$PY" - "$seed" "$(vm_worker_cmd "$seed")" <<'PYEOF' 2>&1 | tee -a "$DISPATCH_LOG"
import os, sys
from utils.dispatch import VMExecutor
seed, cmd = sys.argv[1], sys.argv[2]
vm = VMExecutor()
# 5 ns MD + reextract + MM-PBSA on the V100 → allow a wide ceiling (5 h).
res = vm.execute(cmd, cwd=os.environ["VM_PROJECT_ROOT"], timeout=18000, retry=2)
sys.stdout.write(res.get("stdout", "") or "")
sys.stderr.write(res.get("stderr", "") or "")
print(f"  [VM run {seed}] rc={res['returncode']} elapsed={res['elapsed_s']:.0f}s retried={res['retried']}")
sys.exit(0 if res["returncode"] == 0 else 1)
PYEOF
    return "${PIPESTATUS[0]}"
}

# (3) pull result dir (snapshot subdir + mmpbsa subdir) via VMExecutor.sync_from_vm
vm_pull_seed () {
    local seed="$1"
    local seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"
    "$PY" - "${VM_PROJECT_ROOT}/outputs/${SYSTEM_TAG}_calib_${seed}" "$seedir" \
            "$SNAP_SUBDIR" "$MMPBSA_SUBDIR" <<'PYEOF' 2>&1 | tee -a "$DISPATCH_LOG"
import sys
from utils.dispatch import VMExecutor
remote_dir, local_dir, snap_subdir, mmpbsa_subdir = sys.argv[1:5]
vm = VMExecutor()
ok = vm.sync_from_vm(
    remote_dir, local_dir,
    patterns=[snap_subdir + "/", mmpbsa_subdir + "/"],
    timeout=900,
)
print(f"  [VM pull] {remote_dir} -> {local_dir} ({snap_subdir}/,{mmpbsa_subdir}/): {'OK' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
PYEOF
    return "${PIPESTATUS[0]}"
}

# Whole VM lane: each seed sequentially push→run→pull, idempotent skip.
run_vm_lane () {
    [ "${#VM_SEEDS[@]}" -eq 0 ] && { log "VM lane: no seeds assigned, skip."; return 0; }
    export VM_PROJECT_ROOT
    log "----------------------------------------------------------------"
    log "VM lane (V100 ${VM_SSH_TARGET}): ${VM_SEEDS[*]} (sequential push→run→pull)"
    log "----------------------------------------------------------------"
    # (0) one-time: deploy THIS launcher to the VM so the --vm-worker SSH command
    # (which references scripts/<name> on the VM) can actually start. Abort the
    # whole VM lane on failure rather than silently no-op every seed.
    log "  VM bootstrap: deploy launcher to ${VM_PROJECT_ROOT}/scripts/"
    vm_bootstrap_script || { log "  VM bootstrap FAILED — VM lane cannot start, abort lane."; return 1; }
    local seed seedir
    for seed in "${VM_SEEDS[@]}"; do
        seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"
        if [ -f "${seedir}/${MMPBSA_SUBDIR}/mmpbsa_summary.json" ]; then
            log "  VM ${seed}: SKIP (mmpbsa_summary.json already present locally)"
            continue
        fi
        log "  VM ${seed}: (1/3) push prerequisites + hybrid XML"
        vm_push_seed "$seed" || { log "  VM ${seed}: push FAIL — skip"; continue; }
        log "  VM ${seed}: (2/3) run --vm-worker pipeline on V100"
        vm_run_seed "$seed"  || { log "  VM ${seed}: remote run FAIL — attempting result pull anyway"; }
        log "  VM ${seed}: (3/3) pull result dir"
        vm_pull_seed "$seed" || log "  VM ${seed}: pull FAIL"
    done
    log "VM lane complete."
}

# ----- dry-run: print resolved 2-lane plan + exact host commands + VM transport -----
print_host_seed_plan () {
    local seed="$1"
    local seed_int="${SEED_INT[$seed]:-?}"
    local seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"
    local graft_dir="${PROJ}/outputs/2QKI_WT_calib_${seed}"
    printf '\n  [seed %s  int=%s]\n' "$seed" "$seed_int" | tee -a "$DISPATCH_LOG"
    printf '    out dir   : %s\n' "${seedir#$PROJ/}" | tee -a "$DISPATCH_LOG"
    printf '    graft ref : %s  (exists: %s)\n' "${graft_dir#$PROJ/}" \
        "$([ -f "${graft_dir}/mdresult/2QKI_WT_restrained.dcd" ] && echo yes || echo NO)" | tee -a "$DISPATCH_LOG"
    printf '    baseline  : outputs/%s_calib_%s  (DCD exists: %s)\n' "$SYSTEM_BASE" "$seed" \
        "$([ -f "${PROJ}/outputs/${SYSTEM_BASE}_calib_${seed}/mdresult/${SYSTEM_BASE}_restrained.dcd" ] && echo yes || echo NO)" | tee -a "$DISPATCH_LOG"
    printf '    MD CMD    : CUDA_VISIBLE_DEVICES=%s UPDD_MTR_AMBER14_PATCH=0 UPDD_MMGBSA_PLATFORM=CUDA \\\n' "$MD_GPU" | tee -a "$DISPATCH_LOG"
    printf '                %s %s utils/run_restrained_md.py \\\n' "${MD_PREFIX:-}" "$PY_MD" | tee -a "$DISPATCH_LOG"
    printf '                  --inputdir %s/_md_input --outputdir %s/mdresult \\\n' "${seedir#$PROJ/}" "${seedir#$PROJ/}" | tee -a "$DISPATCH_LOG"
    printf '                  --params_manifest %s/params/MTR_params_manifest.json \\\n' "${seedir#$PROJ/}" | tee -a "$DISPATCH_LOG"
    printf '                  --steps %s --topology %s --binder_chain %s --graph_policy %s \\\n' "$MD_STEPS" "$MD_TOPOLOGY" "$MD_BINDER_CHAIN" "$MD_GRAPH_POLICY" | tee -a "$DISPATCH_LOG"
    printf '                  --target_id %s --dt_fs %s --platform %s --seed %s \\\n' "$MD_TARGET_ID" "$MD_DT" "$MD_PLATFORM" "$seed_int" | tee -a "$DISPATCH_LOG"
    printf '                  --ncaa_label %s --ncaa_code %s\n' "$MD_NCAA_LABEL" "$MD_NCAA_CODE" | tee -a "$DISPATCH_LOG"
    printf '    CHARGE    : local params/MTR_gaff2.xml <- %s (hybrid)\n' "${HYBRID_XML#$PROJ/}" | tee -a "$DISPATCH_LOG"
    printf '    RESCUE    : Pass1b seed=%s @2fs -> Pass2 seed=%s @1fs (R-7 archive)\n' "$((seed_int + 13))" "$((seed_int + 13))" | tee -a "$DISPATCH_LOG"
    printf '    REEXTRACT : reextract_one(system="%s", seed="%s") -> %s/\n' "$REEXTRACT_SYSTEM" "$seed" "$SNAP_SUBDIR" | tee -a "$DISPATCH_LOG"
    printf '    MM-PBSA   : run_mmpbsa.py --md_dir .../%s --protocol %s --ncaa_elem %s\n' "$SNAP_SUBDIR" "$MMPBSA_PROTOCOL" "$MMPBSA_NCAA_ELEM" | tee -a "$DISPATCH_LOG"
}

print_vm_seed_plan () {
    local seed="$1"
    local seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"
    local rdir="${VM_PROJECT_ROOT}/outputs/${SYSTEM_TAG}_calib_${seed}"
    printf '\n  [seed %s  int=%s]  -> V100 (%s)\n' "$seed" "${SEED_INT[$seed]:-?}" "$VM_SSH_TARGET" | tee -a "$DISPATCH_LOG"
    printf '    (1) HOST pre-stage (setup_seed_dir): %s/_md_input/*.pdb + params/{MTR_gaff2.xml(hybrid),MTR_hydrogens.xml,MTR_params_manifest.json}\n' "${seedir#$PROJ/}" | tee -a "$DISPATCH_LOG"
    printf '        VMExecutor.sync_to_vm  : %s/  ->  %s:%s/\n' "${seedir#$PROJ/}" "$VM_SSH_TARGET" "$rdir" | tee -a "$DISPATCH_LOG"
    printf '        (excludes: *.dcd *.nc *.h5 *.chk mdresult/ ; pushes the COMPLETE per-seed _md_input/ + params/)\n' | tee -a "$DISPATCH_LOG"
    printf '    (2) VMExecutor.execute     : ssh %s (cwd=%s):\n' "$VM_SSH_TARGET" "$VM_PROJECT_ROOT" | tee -a "$DISPATCH_LOG"
    printf '          %s\n' "$(vm_worker_cmd "$seed")" | tee -a "$DISPATCH_LOG"
    printf '        (V100 worker SKIPS host-only pre-flight + setup_seed_dir source-copies;\n' | tee -a "$DISPATCH_LOG"
    printf '         charge gate = in-worker grep on the PUSHED params/MTR_gaff2.xml; reads ONLY\n' | tee -a "$DISPATCH_LOG"
    printf '         transport-delivered per-seed dir + VM-resident run_restrained_md/reextract/run_mmpbsa)\n' | tee -a "$DISPATCH_LOG"
    printf '    (3) VMExecutor.sync_from_vm: %s:%s/  ->  %s/\n' "$VM_SSH_TARGET" "$rdir" "${seedir#$PROJ/}" | tee -a "$DISPATCH_LOG"
    printf '        (patterns: %s/ , %s/ ; DCD stays on VM)\n' "$SNAP_SUBDIR" "$MMPBSA_SUBDIR" | tee -a "$DISPATCH_LOG"
}

if [ "$DRY_RUN" -eq 1 ]; then
    log "----------------------------------------------------------------"
    log "DRY-RUN — no MD executed. Resolved 2-lane plan:"
    log "----------------------------------------------------------------"
    # disjointness confirmation (assert already passed above; restate explicitly)
    log "DISJOINTNESS CONFIRMED: host ∩ vm = ∅ ; host ∪ vm = full batch"
    log "  host lane (${#HOST_SEEDS[@]}): ${HOST_SEEDS[*]:-none}"
    log "  vm   lane (${#VM_SEEDS[@]}): ${VM_SEEDS[*]:-none}"
    log "  union   (${#ALL_SEEDS[@]}): ${ALL_SEEDS[*]:-none}"

    log "................................................................"
    log "HOST LANE (RTX 5070 Ti idx ${MD_GPU}) — sequential per-seed commands:"
    if [ "${#HOST_SEEDS[@]}" -eq 0 ]; then
        log "  (none)"
    else
        for seed in "${HOST_SEEDS[@]}"; do print_host_seed_plan "$seed"; done
    fi

    log "................................................................"
    log "VM LANE (V100 ${VM_SSH_TARGET}) — VMExecutor transport per seed:"
    if [ "${#VM_SEEDS[@]}" -eq 0 ]; then
        log "  (none)"
    else
        export VM_PROJECT_ROOT
        log "  (0) one-time VM bootstrap (before per-seed loop):"
        log "      rsync $(basename "$0")  ->  ${VM_SSH_TARGET}:${VM_PROJECT_ROOT}/scripts/"
        log "      (the --vm-worker SSH command references scripts/$(basename "$0") on the VM;"
        log "       hybrid XML + REF_SEED sources are NOT pushed — worker no longer reads them)"
        for seed in "${VM_SEEDS[@]}"; do print_vm_seed_plan "$seed"; done
    fi

    log "----------------------------------------------------------------"
    log "DRY-RUN complete. host=${#HOST_SEEDS[@]} vm=${#VM_SEEDS[@]} (total ${#ALL_SEEDS[@]}). No MD launched."
    log "Lanes run CONCURRENTLY (both backgrounded, wait for both)."
    log "----------------------------------------------------------------"
    exit 0
fi

PBSA_LOG="${LOGDIR}/mmpbsa.log"

# ----- Stage D worker (one MM-PBSA per seed; 4-lane parallel within the lane) --
run_one_pbsa () {
    local seed="$1"
    local seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"
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

# ============================================================================
# Host lane — Stage A (MD+rescue) → Stage C (reextract) → Stage D (MM-PBSA) on
# the seeds assigned to THIS lane (HOST_SEEDS). The SAME function body runs on
# the host (RTX 5070 Ti) and — when invoked over SSH as --vm-worker — on the
# V100 (where HOST_SEEDS has been set to the VM seeds, MD_GPU to the V100 index).
# ============================================================================
run_host_lane () {
    local lane_label="$1"
    [ "${#HOST_SEEDS[@]}" -eq 0 ] && { log "${lane_label}: no seeds assigned, skip."; return 0; }

    # ----- Stage A — MD (sequential single-GPU, with rescue) -----
    log "----------------------------------------------------------------"
    log "${lane_label} Stage A: ${SYSTEM_TAG} MD ${#HOST_SEEDS[@]} seeds (sequential, GPU ${MD_GPU})"
    log "----------------------------------------------------------------"
    local seed
    for seed in "${HOST_SEEDS[@]}"; do
        md_with_rescue "$seed"
    done
    log "${lane_label} Stage A complete."

    if [ "$MD_ONLY" -eq 1 ]; then
        log "${lane_label} --md-only: skipping Stage C/D."
        return 0
    fi

    # ----- Stage C — postl387_v2 snapshot reextract (n=25 K-Means). Calls
    #       reextract_one directly with system=${REEXTRACT_SYSTEM} so the shared
    #       TARGETS list in phase_beta_repbsa_v2_reextract.py is untouched. -----
    log "----------------------------------------------------------------"
    log "${lane_label} Stage C: snapshot reextract (${SNAP_SUBDIR})"
    log "----------------------------------------------------------------"
    conda activate qmmm
    for seed in "${HOST_SEEDS[@]}"; do
        local seedir="${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}"
        local target="${seedir}/${SNAP_SUBDIR}"
        if [ -d "$target" ] && [ "$(ls "$target"/*.pdb 2>/dev/null | wc -l)" -ge 25 ]; then
            log "  ${seed}: reextract SKIP (≥25 pdb already present)"
            continue
        fi
        if [ ! -f "${seedir}/mdresult/${SYSTEM_BASE}_final.pdb" ]; then
            log "  ${seed}: reextract SKIP (no final.pdb — MD did not complete)"
            continue
        fi
        log "  ${seed}: reextract (affinity=[${MMPBSA_AFFINITY:-free}])"
        CUDA_VISIBLE_DEVICES="" $MMPBSA_PREFIX "$PY" -c "
import sys
sys.path.insert(0, 'scripts')
from phase_beta_repbsa_v2_reextract import reextract_one
rec = reextract_one('${REEXTRACT_SYSTEM}', '${seed}', out_subdir='${SNAP_SUBDIR}')
print(f\"{rec['tag']}: {rec['status']} n_saved={rec.get('n_saved')} {rec.get('skip_reason') or ''}\")
" >> "${LOGDIR}/reextract.log" 2>&1
    done
    conda deactivate
    log "${lane_label} Stage C complete."

    # ----- Stage D — 1-traj MM-PBSA scoring (median-Δg ranking; 4-lane parallel) -
    log "----------------------------------------------------------------"
    log "${lane_label} Stage D: MM-PBSA (${MMPBSA_PROTOCOL}, 4 lanes)"
    log "----------------------------------------------------------------"
    local running=0
    for seed in "${HOST_SEEDS[@]}"; do
        run_one_pbsa "$seed" &
        running=$((running + 1))
        if [ "$running" -ge 4 ]; then
            wait -n
            running=$((running - 1))
        fi
    done
    wait
    log "${lane_label} Stage D complete."
}

# ============================================================================
# 2-lane dispatch — host lane + V100 VM lane run CONCURRENTLY (both backgrounded,
# wait for both). Disjoint seeds → disjoint 2QKI_Cp4_hybrid_calib_s{N} dirs → no
# lock/race. When --vm-worker, run only the (local) host lane on the V100.
# ============================================================================
if [ "$VM_WORKER" -eq 1 ]; then
    # Running ON the V100 over SSH — execute the local lane only (no nested VM lane).
    run_host_lane "vm-worker"
else
    log "================================================================"
    log "Launching 2 lanes CONCURRENTLY: host (${#HOST_SEEDS[@]} seeds) ∥ vm (${#VM_SEEDS[@]} seeds)"
    log "================================================================"
    HOST_LANE_LOG="${LOGDIR}/host_lane.log"
    VM_LANE_LOG="${LOGDIR}/vm_lane.log"

    # Host lane (backgrounded). Its own log; dispatch.log still aggregates.
    ( run_host_lane "host-lane" ) > >(tee -a "$HOST_LANE_LOG") 2>&1 &
    HOST_PID=$!

    # VM lane (backgrounded). Skipped cleanly if no VM seeds.
    ( run_vm_lane ) > >(tee -a "$VM_LANE_LOG") 2>&1 &
    VM_PID=$!

    HOST_RC=0; VM_RC=0
    wait "$HOST_PID" || HOST_RC=$?
    wait "$VM_PID"   || VM_RC=$?
    log "Both lanes joined. host_rc=${HOST_RC} vm_rc=${VM_RC}"
fi

# ============================================================================
# Summary — report over the seeds owned by THIS invocation. For the main run
# that is the full union (host ∪ vm, all synced back); for --vm-worker it is the
# worker's local seeds.
# ============================================================================
if [ "$VM_WORKER" -eq 1 ]; then
    SUMMARY_SEEDS=("${HOST_SEEDS[@]}")   # worker owns the VM seeds (set into HOST_SEEDS)
else
    SUMMARY_SEEDS=("${ALL_SEEDS[@]}")
fi
SUCCESS=0
for seed in "${SUMMARY_SEEDS[@]}"; do
    [ -z "$seed" ] && continue
    if [ -f "${PROJ}/outputs/${SYSTEM_TAG}_calib_${seed}/${MMPBSA_SUBDIR}/mmpbsa_summary.json" ]; then
        SUCCESS=$((SUCCESS + 1))
        log "  ✓ ${seed}: MD + reextract + MM-PBSA complete"
    else
        log "  ✗ ${seed}: incomplete (check logs)"
    fi
done
ELAPSED=$(( $(date +%s) - LAUNCH_TS ))
log "================================================================"
log "Track A ${SYSTEM_TAG} MD COMPLETE  elapsed=${ELAPSED}s  ${SUCCESS}/${#SUMMARY_SEEDS[@]} ready"
log "Next: scripts/phase4_qmmm_iva_1traj.py (strip frames now resolve from ${VARIANT_TAG} dirs)"
log "================================================================"
