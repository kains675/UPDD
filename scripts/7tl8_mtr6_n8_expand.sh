#!/bin/bash
# ============================================================================
# 7TL8_MTR6 Sampling Expansion N=3 → N=8 — Convergence Index Validation
# ============================================================================
# Goal: Validate the user's σ_btwn-as-classifier hypothesis by extending the
#       7TL8_MTR6 family from N=3 (s7/s19/s42, σ_btwn=5.547, CI=0.70) to N=8
#       with 5 new MD seeds (s55/s73/s89/s107/s131; prime-gap spacing).
#
# Pipeline per seed:
#   1. Mirror outputs/7TL8_MTR6_calib_s7/{_md_input,params}/ into NEW dir
#      outputs/7TL8_MTR6_calib_s<NEW>/
#   2. MD: utils/run_restrained_md.py linear / chain B / 5 ns @ 2 fs
#      (md_simulation env, mirrors phase_beta_stage1_repair_md.sh pattern).
#   3. Snapshot extract: invoke save_snapshots() (v54 PBC + v55 CONECT) into
#      snapshots_n25_postl387_patch_v2/ (consistent with #90).
#   4. PBSA verify: assert all 25 PDBs prev_C↔ncAA_N < 5 Å, CONECT > 10.
#   5. MM-PBSA: scripts/run_mmpbsa.py --protocol 1traj into
#      mmpbsa_results_postl387_v2/ (qmmm env, mirrors #90).
#   6. Aggregate: append seed result + recompute σ_btwn into
#      outputs/analysis/phase_beta_repbsa_v2_aggregate.json (extend, do not
#      break v0.7 format).
#
# Idempotent: skips MD if mdresult/<basename>_restrained.dcd > 100MB exists;
#             skips snapshots if 25 PDBs already in v2 dir;
#             skips PBSA if mmpbsa_summary.json exists.
#
# Safety:
#   - DO NOT touch outputs/7TL8_MTR6_calib_s{7,19,42}/ (R-7 preservation).
#   - DO NOT modify any source code (extract_snapshots.py, run_restrained_md.py
#     are immutable; reference only).
#   - GPU contention check: bail if other compute apps detected
#     (UPDD memory: GPU exclusive during user compute).
#
# Per-seed wall ~11-13h GPU; 5 seeds sequential ~55-65h (target ~70h budget,
# 5-day hard limit 120h).
# ============================================================================

set -u

PROJ=/home/san/UPDD_proj
cd "$PROJ"

# Conda env paths (mirrors phase_beta_stage1_repair_md.sh + phase_beta_repbsa_v2.sh)
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PY_PBSA=/home/san/miniconda3/envs/qmmm/bin/python

source /home/san/miniconda3/etc/profile.d/conda.sh

# Production MD env vars (mirrors stage1_repair_md.sh)
export UPDD_MTR_AMBER14_PATCH=1
export UPDD_MMGBSA_PLATFORM=CUDA

LAUNCH_TS=$(date +%s)
GPU_BUDGET_H=120          # 5-day hard limit
TARGET_BUDGET_H=70        # ~3-day target
STEPS=2500000             # 5 ns at 2 fs (matches s7 successful reference)
N_SNAP=25
TARGET_ID=7TL8
NCAA_LABEL=MTR
NCAA_CODE=MTR
BASENAME=7TL8_MTR6
REF_SEED=s7               # Canonical input/params source

# 5 NEW seeds (prime-gap spacing, all distinct from existing s7/s19/s42)
NEW_SEEDS=(55 73 89 107 131)

LOGDIR=${PROJ}/outputs/analysis/7tl8_mtr6_n8_expand_logs
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"
AGGREGATE_JSON="${PROJ}/outputs/analysis/phase_beta_repbsa_v2_aggregate.json"

log() {
    echo "[7tl8_n8_expand $(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$DISPATCH_LOG"
}

check_budget() {
    local now=$(date +%s)
    local elapsed_s=$(( now - LAUNCH_TS ))
    local elapsed_h=$(( elapsed_s / 3600 ))
    if [ "$elapsed_h" -ge "$GPU_BUDGET_H" ]; then
        log "[CATR-ESCALATE] GPU budget ${GPU_BUDGET_H}h exceeded (elapsed=${elapsed_h}h) — HALT"
        return 1
    fi
    if [ "$elapsed_h" -ge "$TARGET_BUDGET_H" ]; then
        log "WARN target budget ${TARGET_BUDGET_H}h exceeded (elapsed=${elapsed_h}h) — continuing under hard cap"
    fi
    return 0
}

check_gpu_exclusive() {
    # GPU exclusive during user compute — bail if a non-self compute proc holds GPU.
    # We track our launched MD/PBSA child PIDs and tolerate them. For external
    # SCF/MD/etc started by the user, halt and CATR-escalate.
    local self_pgid
    self_pgid=$(ps -o pgid= -p $$ 2>/dev/null | tr -d ' ')
    local apps
    apps=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null)
    local n_external=0
    while IFS= read -r pid; do
        [ -z "$pid" ] && continue
        local pid_pgid
        pid_pgid=$(ps -o pgid= -p "$pid" 2>/dev/null | tr -d ' ')
        if [ -n "$pid_pgid" ] && [ "$pid_pgid" != "$self_pgid" ]; then
            n_external=$((n_external + 1))
        fi
    done <<< "$apps"
    if [ "$n_external" -gt 0 ]; then
        log "WARN GPU has ${n_external} external compute apps; pausing 60s"
        sleep 60
        # Re-check
        apps=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null)
        n_external=0
        while IFS= read -r pid; do
            [ -z "$pid" ] && continue
            local pid_pgid
            pid_pgid=$(ps -o pgid= -p "$pid" 2>/dev/null | tr -d ' ')
            if [ -n "$pid_pgid" ] && [ "$pid_pgid" != "$self_pgid" ]; then
                n_external=$((n_external + 1))
            fi
        done <<< "$apps"
        if [ "$n_external" -gt 0 ]; then
            log "[CATR-ESCALATE] GPU still has ${n_external} external compute apps — HALT"
            return 1
        fi
    fi
    return 0
}

setup_seed_dir() {
    local seed_int="$1"
    local outdir="${PROJ}/outputs/7TL8_MTR6_calib_s${seed_int}"
    local refdir="${PROJ}/outputs/7TL8_MTR6_calib_${REF_SEED}"

    if [ -d "${outdir}/_md_input" ] && [ -d "${outdir}/params" ]; then
        log "  setup s${seed_int}: dirs already present"
        return 0
    fi

    log "  setup s${seed_int}: mirroring _md_input/ + params/ from ${REF_SEED}"
    mkdir -p "${outdir}/_md_input" "${outdir}/params" "${outdir}/mdresult"
    # Copy MD input (input PDB + ncaa params used by run_restrained_md.py)
    cp -n "${refdir}/_md_input/${BASENAME}.pdb" "${outdir}/_md_input/" || true
    cp -n "${refdir}/_md_input/${BASENAME}_renum.pdb" "${outdir}/_md_input/" || true
    cp -n "${refdir}/_md_input/MTR_gaff2.xml" "${outdir}/_md_input/" || true
    cp -n "${refdir}/_md_input/MTR_hydrogens.xml" "${outdir}/_md_input/" || true
    cp -n "${refdir}/_md_input/MTR_params_manifest.json" "${outdir}/_md_input/" || true
    # Copy canonical params/ (xml_path points LOCALLY)
    cp -n "${refdir}/params/MTR_gaff2.xml" "${outdir}/params/" || true
    cp -n "${refdir}/params/MTR_hydrogens.xml" "${outdir}/params/" || true
    # Re-write params/MTR_params_manifest.json with this seed's local path
    "$PY_PBSA" - <<PYEOF
import json
src = "${refdir}/params/MTR_params_manifest.json"
dst = "${outdir}/params/MTR_params_manifest.json"
with open(src) as f:
    m = json.load(f)
m["xml_path"] = "${outdir}/params/MTR_gaff2.xml"
m["hydrogens_path"] = "${outdir}/params/MTR_hydrogens.xml"
with open(dst, "w") as f:
    json.dump(m, f, indent=2)
PYEOF
    return 0
}

run_md_for_seed() {
    local seed_int="$1"
    local outdir="${PROJ}/outputs/7TL8_MTR6_calib_s${seed_int}"
    local manifest="${outdir}/params/MTR_params_manifest.json"
    local dcd="${outdir}/mdresult/${BASENAME}_restrained.dcd"
    local logf="${LOGDIR}/md_s${seed_int}.log"

    if [ -f "$dcd" ] && [ "$(stat -c %s "$dcd")" -gt 100000000 ]; then
        log "  MD s${seed_int}: SKIP (DCD exists, $(du -h "$dcd" | cut -f1))"
        return 0
    fi

    check_budget || return 2
    check_gpu_exclusive || return 2

    log "  MD s${seed_int}: launching 5ns linear MTR (env=md_simulation)"
    conda activate md_simulation
    UPDD_MMGBSA_PLATFORM=CUDA "$PY_MD" utils/run_restrained_md.py \
        --inputdir "${outdir}/_md_input" \
        --outputdir "${outdir}/mdresult" \
        --params_manifest "$manifest" \
        --steps "$STEPS" \
        --topology linear \
        --binder_chain B \
        --graph_policy strict \
        --target_id "$TARGET_ID" \
        --dt_fs 2.0 \
        --platform CUDA \
        --seed "$seed_int" \
        --ncaa_label "$NCAA_LABEL" \
        --ncaa_code "$NCAA_CODE" \
        > "$logf" 2>&1
    local rc=$?
    conda deactivate

    if [ "$rc" -ne 0 ] || ! [ -f "$dcd" ]; then
        log "  MD s${seed_int}: FAIL rc=${rc} (NVT explosion or other) — see ${logf}"
        # Try seed+13 reseed (canonical UPDD precedent from stage1_repair)
        local reseed=$((seed_int + 13))
        log "  MD s${seed_int}: pass1b reseed=${reseed}"
        rm -f "${outdir}/mdresult/"*EXPLODED* "${outdir}/mdresult/"*_partial_* \
              "${outdir}/mdresult/"*_final.pdb "${outdir}/mdresult/"*_md.log 2>/dev/null
        local logf2="${LOGDIR}/md_s${seed_int}_pass1b_reseed${reseed}.log"
        conda activate md_simulation
        UPDD_MMGBSA_PLATFORM=CUDA "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${outdir}/_md_input" \
            --outputdir "${outdir}/mdresult" \
            --params_manifest "$manifest" \
            --steps "$STEPS" \
            --topology linear \
            --binder_chain B \
            --graph_policy strict \
            --target_id "$TARGET_ID" \
            --dt_fs 2.0 \
            --platform CUDA \
            --seed "$reseed" \
            --ncaa_label "$NCAA_LABEL" \
            --ncaa_code "$NCAA_CODE" \
            > "$logf2" 2>&1
        rc=$?
        conda deactivate
        if [ "$rc" -ne 0 ] || ! [ -f "$dcd" ]; then
            log "  MD s${seed_int}: FAIL after pass1b — escalating"
            return 1
        fi
        log "  MD s${seed_int}: PASS via reseed ${reseed}"
    else
        log "  MD s${seed_int}: PASS (DCD $(du -h "$dcd" | cut -f1))"
    fi
    return 0
}

extract_snapshots_for_seed() {
    local seed_int="$1"
    local outdir="${PROJ}/outputs/7TL8_MTR6_calib_s${seed_int}"
    local snapdir="${outdir}/snapshots_n25_postl387_patch_v2"
    local n_existing
    n_existing=$(ls "${snapdir}"/*_snap*_f*.pdb 2>/dev/null | wc -l)
    if [ "$n_existing" -ge "$N_SNAP" ]; then
        log "  extract s${seed_int}: SKIP (n_existing=${n_existing})"
        return 0
    fi
    log "  extract s${seed_int}: launching save_snapshots (v54 PBC + v55 CONECT) inline"
    local logf="${LOGDIR}/extract_s${seed_int}.log"
    # The phase_beta_repbsa_v2_reextract.py TARGETS list does not include the
    # new seeds, so use the inline save_snapshots invocation directly.
    conda activate qmmm
    "$PY_PBSA" - <<PYEOF > "$logf" 2>&1
import sys, os, warnings, logging, json
warnings.filterwarnings("ignore")
logging.getLogger("mdtraj").setLevel(logging.ERROR)
logging.getLogger("openmm").setLevel(logging.ERROR)
sys.path.insert(0, "${PROJ}")
sys.path.insert(0, "${PROJ}/utils")
import mdtraj as md
from utils.extract_snapshots import cluster_and_select, save_snapshots

outdir = "${outdir}"
snapdir = os.path.join(outdir, "snapshots_n25_postl387_patch_v2")
os.makedirs(snapdir, exist_ok=True)
dcd = os.path.join(outdir, "mdresult", "${BASENAME}_restrained.dcd")
top_candidates = [
    os.path.join(outdir, "mdresult", "${BASENAME}_final.pdb"),
    os.path.join(outdir, "_md_input", "${BASENAME}_renum.pdb"),
    os.path.join(outdir, "_md_input", "${BASENAME}.pdb"),
]
top = next((c for c in top_candidates if os.path.isfile(c)), None)
if top is None:
    raise SystemExit("no topology")
traj = md.load(dcd, top=top)
frames = cluster_and_select(traj, ${N_SNAP}, binder_chain="B")
saved = save_snapshots(traj, sorted(frames), snapdir, "${BASENAME}")
print(f"saved={len(saved)} frames={sorted(frames)}")
PYEOF
    local rc=$?
    conda deactivate
    if [ "$rc" -ne 0 ]; then
        log "  extract s${seed_int}: FAIL rc=${rc} — see ${logf}"
        return 1
    fi
    n_after=$(ls "${snapdir}"/*_snap*_f*.pdb 2>/dev/null | wc -l)
    log "  extract s${seed_int}: PASS (n_pdb=${n_after})"
    return 0
}

verify_snapshots_for_seed() {
    local seed_int="$1"
    local outdir="${PROJ}/outputs/7TL8_MTR6_calib_s${seed_int}"
    local snapdir="${outdir}/snapshots_n25_postl387_patch_v2"
    local logf="${LOGDIR}/verify_s${seed_int}.log"
    log "  verify s${seed_int}: prev_C↔ncAA_N + CONECT integrity"
    conda activate qmmm
    "$PY_PBSA" - <<PYEOF > "$logf" 2>&1
import sys, os, glob, json
import numpy as np
sys.path.insert(0, "${PROJ}")
import mdtraj as md
import warnings, logging
warnings.filterwarnings("ignore")
logging.getLogger("mdtraj").setLevel(logging.ERROR)

snapdir = "${snapdir}"
pdbs = sorted(glob.glob(os.path.join(snapdir, "*_snap*_f*.pdb")))
n_pdb = len(pdbs)
n_intact = 0
n_conect_ok = 0
worst_dist = 0.0
worst_pdb = ""
for p in pdbs:
    n_c = sum(1 for line in open(p) if line.startswith("CONECT"))
    if n_c > 10:
        n_conect_ok += 1
    try:
        t = md.load(p)
        top = t.topology
        # find MTR (ncAA) N atom and prev residue C atom
        for ch in top.chains:
            residues = list(ch.residues)
            for i, r in enumerate(residues):
                if r.name == "MTR" and i > 0:
                    prev = residues[i-1]
                    try:
                        c_idx = next(a.index for a in prev.atoms if a.name == "C")
                        n_idx = next(a.index for a in r.atoms if a.name == "N")
                        d = np.linalg.norm(t.xyz[0, c_idx] - t.xyz[0, n_idx]) * 10.0  # nm→Å
                        if d < 5.0:
                            n_intact += 1
                        if d > worst_dist:
                            worst_dist = d
                            worst_pdb = os.path.basename(p)
                    except StopIteration:
                        pass
                    break
    except Exception as e:
        print(f"  WARN parse {p}: {e}")
print(f"n_pdb={n_pdb} n_intact_le5A={n_intact} n_conect_gt10={n_conect_ok} worst={worst_dist:.2f}A ({worst_pdb})")
status = "OK" if (n_intact >= 24 and n_conect_ok >= 24) else "WARN"
print(f"verdict: {status}")
with open("${LOGDIR}/verify_s${seed_int}.json", "w") as f:
    json.dump({"n_pdb": n_pdb, "n_intact": n_intact, "n_conect_ok": n_conect_ok,
               "worst_dist_A": worst_dist, "worst_pdb": worst_pdb,
               "status": status}, f, indent=2)
PYEOF
    local rc=$?
    conda deactivate
    if [ "$rc" -ne 0 ]; then
        log "  verify s${seed_int}: WARN rc=${rc}"
    fi
    return 0
}

run_pbsa_for_seed() {
    local seed_int="$1"
    local outdir="${PROJ}/outputs/7TL8_MTR6_calib_s${seed_int}"
    local snapdir="${outdir}/snapshots_n25_postl387_patch_v2"
    local pbsadir="${outdir}/mmpbsa_results_postl387_v2"
    local summary="${pbsadir}/mmpbsa_summary.json"
    local logf="${LOGDIR}/pbsa_s${seed_int}.log"
    if [ -f "$summary" ]; then
        log "  pbsa s${seed_int}: SKIP (mmpbsa_summary.json exists)"
        return 0
    fi
    if ! ls "${snapdir}"/*.pdb >/dev/null 2>&1; then
        log "  pbsa s${seed_int}: FAIL (no snapshots)"
        return 1
    fi
    check_budget || return 2
    check_gpu_exclusive || return 2

    mkdir -p "$pbsadir"
    log "  pbsa s${seed_int}: launching MM-PBSA 1traj (env=qmmm)"
    conda activate qmmm
    UPDD_MMGBSA_PLATFORM=CUDA "$PY_PBSA" "${PROJ}/scripts/run_mmpbsa.py" \
        --md_dir "$snapdir" \
        --outputdir "$pbsadir" \
        --ncaa_elem MTR \
        --receptor_chain A --binder_chain B \
        --target_id "$TARGET_ID" \
        --protocol 1traj \
        > "$logf" 2>&1
    local rc=$?
    conda deactivate
    if [ "$rc" -ne 0 ] || [ ! -f "$summary" ]; then
        log "  pbsa s${seed_int}: FAIL rc=${rc} — see ${logf}"
        return 1
    fi
    log "  pbsa s${seed_int}: PASS"
    return 0
}

aggregate_seed() {
    # Append/refresh seed in aggregate JSON, recompute σ_btwn for 7TL8_MTR6.
    local seed_int="$1"
    local outdir="${PROJ}/outputs/7TL8_MTR6_calib_s${seed_int}"
    local summary="${outdir}/mmpbsa_results_postl387_v2/mmpbsa_summary.json"
    if [ ! -f "$summary" ]; then
        log "  aggregate s${seed_int}: SKIP (no summary)"
        return 1
    fi
    log "  aggregate s${seed_int}: updating ${AGGREGATE_JSON}"
    "$PY_PBSA" - <<PYEOF
import json, os, statistics
agg_path = "${AGGREGATE_JSON}"
summary_path = "${summary}"
seed = "s${seed_int}"
with open(agg_path) as f:
    agg = json.load(f)
with open(summary_path) as f:
    summ = json.load(f)
# extract per-snapshot ΔG values from the summary
# (mirrors phase_beta_repbsa_v2_aggregate.py logic)
def mu_sigma(values):
    if not values:
        return 0.0, 0.0
    n = len(values)
    m = sum(values)/n
    if n < 2:
        return m, 0.0
    s = (sum((x-m)**2 for x in values) / (n-1))**0.5
    return m, s
# Try multiple keys for the per-snapshot list
records = summ.get("records") or summ.get("entries") or summ.get("snapshots") or []
delta_gs = []
for r in records:
    if not isinstance(r, dict):
        continue
    v = r.get("delta_g") or r.get("dG") or r.get("binding_kcal") or r.get("ΔG_total")
    if v is None:
        # nested?
        for k in ("mmpbsa", "results", "result"):
            if k in r and isinstance(r[k], dict):
                v = r[k].get("delta_g") or r[k].get("dG")
                if v is not None:
                    break
    if v is not None:
        try:
            fv = float(v)
            if fv == fv and abs(fv) < 1e6:  # NaN/Inf filter
                delta_gs.append(fv)
        except (TypeError, ValueError):
            pass
# Fallback: read individual snapshot JSONs
if not delta_gs:
    pbsadir = os.path.dirname(summary_path)
    import glob
    for sp in sorted(glob.glob(os.path.join(pbsadir, "*_mmpbsa.json"))):
        try:
            with open(sp) as f:
                sj = json.load(f)
            v = sj.get("delta_g") or sj.get("dG") or sj.get("binding_kcal")
            if v is None and "result" in sj:
                v = sj["result"].get("delta_g") or sj["result"].get("dG")
            if v is None and "mmpbsa" in sj:
                v = sj["mmpbsa"].get("delta_g") or sj["mmpbsa"].get("dG")
            if v is not None:
                fv = float(v)
                if fv == fv and abs(fv) < 1e6:
                    delta_gs.append(fv)
        except Exception:
            pass
n = len(delta_gs)
mu, sw = mu_sigma(delta_gs)
n_pos = sum(1 for x in delta_gs if x > 0)
median = statistics.median(delta_gs) if delta_gs else 0.0
mn = min(delta_gs) if delta_gs else 0.0
mx = max(delta_gs) if delta_gs else 0.0
print(f"seed={seed} n={n} mu={mu:.4f} sigma_w={sw:.4f} median={median:.4f} n_pos={n_pos}")

fam = agg["families"].get("7TL8_MTR6")
if fam is None:
    fam = {"family": "7TL8_MTR6", "n_seeds": 0, "rows": []}
    agg["families"]["7TL8_MTR6"] = fam
# Update or append seed
existing = next((r for r in fam["rows"] if r.get("seed") == seed), None)
row = {
    "seed": seed,
    "status": "OK",
    "source_dir": "mmpbsa_results_postl387_v2",
    "n": n,
    "mu": mu,
    "sigma_w": sw,
    "median": median,
    "n_pos": n_pos,
    "min": mn,
    "max": mx,
}
if existing is not None:
    fam["rows"] = [r if r.get("seed") != seed else row for r in fam["rows"]]
else:
    fam["rows"].append(row)
fam["n_seeds"] = len(fam["rows"])
mus = [r["mu"] for r in fam["rows"]]
sws = [r["sigma_w"] for r in fam["rows"]]
M = sum(mus)/len(mus)
sigma_btwn = (sum((m-M)**2 for m in mus) / (len(mus)-1))**0.5 if len(mus) > 1 else 0.0
sigma_w_med = statistics.median(sws) if sws else 0.0
SE = sigma_btwn / (len(mus)**0.5) if mus else 0.0
CI95 = SE * 1.96
z_SE = abs(M) / SE if SE > 0 else 0.0
if z_SE >= 3:
    tier = "Tier-1"
elif z_SE >= 2:
    tier = "Tier-2"
else:
    tier = "Tier-3"
fam["M"] = M
fam["sigma_btwn"] = sigma_btwn
fam["sigma_w_med"] = sigma_w_med
fam["SE"] = SE
fam["CI95"] = CI95
fam["z_SE"] = z_SE
fam["tier"] = tier
print(f"7TL8_MTR6 family @ N={len(mus)}: ⟨ΔG⟩={M:.3f} σ_btwn={sigma_btwn:.3f} CI={sigma_btwn/abs(M):.3f} z_SE={z_SE:.2f} tier={tier}")
with open(agg_path, "w") as f:
    json.dump(agg, f, indent=2)
print("aggregate JSON updated")
PYEOF
    return $?
}

# ============================================================================
# Main loop
# ============================================================================
log "========================================================"
log "7TL8_MTR6 N=3→N=8 sampling expansion dispatch"
log "  GPU budget hard:   ${GPU_BUDGET_H} h (5 days)"
log "  GPU budget target: ${TARGET_BUDGET_H} h (~3 days)"
log "  Steps per seed:    ${STEPS} (5 ns @ 2 fs)"
log "  Snapshots per:     ${N_SNAP}"
log "  New seeds:         ${NEW_SEEDS[@]}"
log "  Aggregate JSON:    ${AGGREGATE_JSON}"
log "========================================================"

for seed_int in "${NEW_SEEDS[@]}"; do
    log "==== seed s${seed_int} START ===="
    setup_seed_dir "$seed_int" || { log "[CATR-ESCALATE] setup s${seed_int} FAIL"; continue; }
    run_md_for_seed "$seed_int"
    rc=$?
    if [ "$rc" -eq 2 ]; then
        log "[CATR-ESCALATE] budget/GPU halt at MD s${seed_int} — stopping remaining"
        break
    elif [ "$rc" -ne 0 ]; then
        log "[CATR-ESCALATE] MD s${seed_int} FAIL — skipping to next seed"
        continue
    fi
    extract_snapshots_for_seed "$seed_int" || { log "extract s${seed_int} FAIL — skip"; continue; }
    verify_snapshots_for_seed "$seed_int"
    run_pbsa_for_seed "$seed_int"
    rc=$?
    if [ "$rc" -eq 2 ]; then
        log "[CATR-ESCALATE] budget/GPU halt at PBSA s${seed_int} — stopping remaining"
        break
    elif [ "$rc" -ne 0 ]; then
        log "[CATR-ESCALATE] PBSA s${seed_int} FAIL — skipping aggregation"
        continue
    fi
    aggregate_seed "$seed_int" || log "aggregate s${seed_int} WARN (data preserved)"
    log "==== seed s${seed_int} DONE ===="
done

log "========================================================"
log "All new seeds processed; updating final summary"
log "========================================================"

# Write a small status JSON for cron to poll
"$PY_PBSA" - <<PYEOF
import json, os
status = {"timestamp": "$(date --iso-8601=seconds)", "seeds": []}
for seed in [55, 73, 89, 107, 131]:
    summ = f"/home/san/UPDD_proj/outputs/7TL8_MTR6_calib_s{seed}/mmpbsa_results_postl387_v2/mmpbsa_summary.json"
    dcd = f"/home/san/UPDD_proj/outputs/7TL8_MTR6_calib_s{seed}/mdresult/${BASENAME}_restrained.dcd"
    snap_dir = f"/home/san/UPDD_proj/outputs/7TL8_MTR6_calib_s{seed}/snapshots_n25_postl387_patch_v2"
    n_pdb = len([f for f in os.listdir(snap_dir) if f.endswith('.pdb')]) if os.path.isdir(snap_dir) else 0
    status["seeds"].append({
        "seed": seed,
        "dcd_exists": os.path.isfile(dcd),
        "dcd_mb": round(os.path.getsize(dcd)/1e6, 1) if os.path.isfile(dcd) else 0,
        "n_pdb": n_pdb,
        "pbsa_done": os.path.isfile(summ),
    })
with open("${LOGDIR}/dispatch_status.json", "w") as f:
    json.dump(status, f, indent=2)
print(json.dumps(status, indent=2))
PYEOF

log "dispatch complete"
log "========================================================"
