#!/bin/bash
# ============================================================================
# WT positive control (NO-GO specificity check) — dt=1fs restrained-protocol
# pose-engagement reproduction on the STRONGEST binder (2QKI_WT).
# (P3 #107, follow-up to the D5 d3refit-MTR NO-GO: 4/4 seeds 0% engaged.)
#
# Purpose (specificity positive control ONLY): the D5 pose-recovery gate scored
# the D3-refit MTR force field as NO-GO (fully disengaged by 25ns, all seeds).
# That NO-GO can have two distinct roots that the D5 run alone cannot separate:
#   (a) MTR-charge-specific  — the refit MTR FF genuinely fails to hold the pose;
#   (b) pipeline / dt=1fs artefact — the dt=1fs restrained 25ns protocol itself
#       (or the contact metric) disengages ANY binder, even a strong one.
# This launcher runs the strongest binder (WT) through the SAME dt=1fs / 25ns
# protocol + the SAME contact-engagement metric, so the WT outcome discriminates:
#   WT engaged  => the protocol CAN report engagement => the D5 NO-GO is NOT a
#                  dt=1fs/pipeline artefact => d3refit NO-GO is charge-specific.
#   WT disengage => the dt=1fs protocol disengages even a strong binder => the
#                  dt=1fs artefact hypothesis takes priority; the d3refit
#                  "MTR-specific" conclusion is SUSPENDED pending diagnosis.
# It is a GEOMETRIC go/no-go on the binding pose — NOT a free-energy / affinity
# computation.
#
# ------------------------------------------------------------------------
# !! ANCHOR ASYMMETRY DISCLOSURE (read before interpreting — load-bearing) !!
# ------------------------------------------------------------------------
# The D5 d3refit cohort restrains the MTR (res-4) backbone N/CA/C at k=100
# kJ/mol/nm^2 because run_restrained_md.py keys its position restraint off the
# ncAA residue name (ctx.xml_res_name). The WT path uses --ncaa none, for which
# run_restrained_md.py NEVER sets xml_res_name (it stays None) and therefore
# the restraint block `if ctx.xml_res_name:` (build stage) is SKIPPED entirely:
# a WT run is, by construction of the current MD engine, a FULLY UNRESTRAINED
# free MD (n_restrained = 0). There is NO CLI flag on run_restrained_md.py that
# can anchor an arbitrary (res-4 Trp4) backbone for a WT run; the only anchor
# mechanism is the ncAA-name-keyed one. So this WT control is:
#   * APPLES-TO-APPLES on: system family (2QKI), topology (cyclic_ss),
#     timestep (dt=1fs), length (25ns), platform, contact metric / bands.
#   * NOT apples-to-apples on: the res-4 backbone anchor. D5 res-4 is anchored
#     (k=100); WT res-4 (Trp4) is FREE (unrestrained, like the whole solute).
# Consequence for interpretation: a WT control that runs FREE is the WEAKER
# (more conservative) test of "does the protocol disengage a strong binder?":
#   - If even FREE WT stays engaged -> engagement is robust w/o any anchor ->
#     the D5 disengagement is NOT a generic protocol artefact (the d3refit run
#     even had the HELPFUL res-4 anchor and still disengaged) -> strengthens
#     the charge-specific reading. GREEN is then conservative / well-founded.
#   - If FREE WT disengages -> ambiguous: could be the missing anchor OR a
#     dt=1fs artefact; cannot cleanly blame dt=1fs alone. RED must be reported
#     as "WT-free disengaged" and escalated (an anchored-WT rerun would need an
#     MD-engine change + a separate science-review sign-off; out of scope here).
# This asymmetry is a property of the existing run_restrained_md.py engine, NOT
# introduced by this launcher. It is surfaced here, in the driver log, and in
# the handback so the decision is read with the correct caveat (R-18).
# ------------------------------------------------------------------------
#
# 3-OUTCOME PRE-REGISTER (decided BEFORE any data, echoed to driver log):
#   GREEN  (pipeline valid): pooled mean contacts WT-like (~92, the ADR-0009 WT
#          band) / high frac_engaged AND sustained-to-25ns (no late-window
#          disengagement). => the dt=1fs / 25ns protocol + contact metric CAN
#          report engagement for a strong binder => the D5 d3refit NO-GO is NOT
#          a dt=1fs/pipeline artefact => d3refit NO-GO is charge-specific.
#   RED    (dt=1fs/pipeline suspect): WT also disengages (pooled mean in the
#          gamma band ~12 / frac_engaged < ~0.20 / disengaged by 25ns). => even
#          a strong binder disengages under this protocol => the dt=1fs artefact
#          hypothesis takes priority; the d3refit "MTR-specific" conclusion is
#          SUSPENDED pending a dt=1fs / anchor-asymmetry diagnosis. (Report as
#          "WT-free disengaged" given the anchor caveat above.)
#   AMBIGUOUS (intermediate): seeds split (one engaged, one disengaged) OR
#          pooled frac_engaged 0.2~0.7 OR time-series shows late disengagement.
#          Report the distribution, grow n, NO magnitude claim.
#
# FORBIDDEN (R-11/R-18): any magnitude / sign / lambda / binding-affinity /
# ddG claim. Contact count is a geometric observable only. A dt=1fs crash here
# is logged loudly and HALTED (not silently rescued) — but note WT at dt=1fs is
# the canonical amber FF and is not expected to carry the H1 MTR-improper risk.
#
# Metric (PRIMARY go/no-go DECISION, locked, identical to D5):
#   n=25 snapshot extraction (cluster_and_select + CONECT-v55 save_snapshots,
#   ADR-0009 SSOT pipeline) -> snapshots_n25_postl387_patch_v2/ ->
#   contact_engagement_analysis.py (locked cutoff 4.0A / n_engaged>=20 /
#   chainB<->chainA). Directly band-comparable to the FIXED WT/gamma/patchoff
#   ADR-0009 bands and to the D5 d3refit cohort.
#
# Runs (2 fresh dt=1fs WT seeds, sequential on a single GPU = 5070 Ti;
# concurrent execution forbidden = contention; V100 not attached):
#   run  seed  dt_fs  --steps    DCD_INTERVAL  length  outdir suffix
#   1    411   1.0    25000000   2000 (=2ps)    25ns    2QKI_WT_d5posctrl_s411_dt1fs
#   2    433   1.0    25000000   2000 (=2ps)    25ns    2QKI_WT_d5posctrl_s433_dt1fs
# n=2 = WT is the strongest binder; 2 fresh dt=1fs seeds are sufficient for the
# discriminating control (a clean engaged/disengaged readout). Length = 25ns
# (= 25M steps at dt=1fs) to MATCH the D5 cohort length (band-comparability) and
# to capture slow / progressive disengagement (a short window would FALSELY
# score a late-disengaging seed as engaged).
#
# WT input provenance: mirrored from the validated WT calib cohort
# outputs/2QKI_WT_calib_s7/_md_input/ (2QKI_WT.pdb + 2QKI_WT_renum.pdb) — the
# SAME inputs that fed the ADR-0009 WT band. The MD invocation reproduces the
# validated WT calib launch (scripts/t1_n30_cp4_expand.sh Stage B) verbatim
# except for dt (2.0->1.0), length (2.5M->25M steps) and the dense DCD interval.
# WT carries NO MTR — the frozen params/MTR_gaff2_hybrid.xml /
# MTR_gaff2_layer3d_d3refit.xml are NEVER read or touched by this launcher.
#
# Dirs: NEW per-seed cohort dirs outputs/2QKI_WT_d5posctrl_s{seed}_dt1fs. The
# existing 2QKI_WT_calib_* dirs are PRESERVED and never clobbered (R-7).
#
# Honest tally: each run is judged by the presence of an _EXPLODED_ marker (NOT
# by exit code "DONE" and NOT by _final.pdb, which is also written on a
# recoverable crash). Poisson rule-of-three on crash-rate (N_CRASH/2 is a fact,
# not "0%").
#
# Usage:
#   scripts/wt_positive_control_dt1fs.sh                        # 2 runs, sequential
#   CUDA_VISIBLE_DEVICES=0 scripts/wt_positive_control_dt1fs.sh # pin GPU
#   DRY_RUN=1 scripts/wt_positive_control_dt1fs.sh              # print plan, no MD
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ" || { echo "cannot cd $PROJ"; exit 2; }

source /home/san/miniconda3/etc/profile.d/conda.sh

# ----- knobs (env-overridable) -----
DRY_RUN="${DRY_RUN:-0}"
# Dense DCD interval for the time-series / sustained-to-25ns diagnostic
# (2ps/frame, matches the D5 d3refit cohort). Override via UPDD_MD_DCD_INTERVAL.
DCD_INTERVAL="${UPDD_MD_DCD_INTERVAL:-2000}"

# ----- CPU affinity (d5_pose_recovery_production.sh / d4 logic) -----
PHYS_CORES=$(lscpu -p 2>/dev/null | grep -v "^#" | awk -F, '{print $2}' | sort -nu | wc -l)
LOGICAL_CORES=$(nproc)
if [ "${PHYS_CORES:-0}" -ge 4 ] && [ "${PHYS_CORES:-0}" -le 8 ]; then
    MD_AFFINITY="0,${PHYS_CORES}"
    MD_PREFIX="taskset -c ${MD_AFFINITY}"
else
    MD_PREFIX=""
fi

# ----- WT input provenance (validated WT calib cohort; NO MTR, NO params) -----
WT_REF_DIR="${PROJ}/outputs/2QKI_WT_calib_s7"

# ----- snapshot extraction SSOT (ADR-0009 / section 4.3.2) -----
SNAP_SUBDIR="snapshots_n25_postl387_patch_v2"
N_SNAPSHOTS=25

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/wt_positive_control_${TS_TAG}"
mkdir -p "$LOGDIR"
mkdir -p "${PROJ}/logs"
DISPATCH_LOG="${LOGDIR}/dispatch.log"
DRIVER_LOG="${PROJ}/logs/wt_positive_control.log"

log () {
    printf '[wt_posctrl %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG" "$DRIVER_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "WT positive control (D5 NO-GO specificity) START — 2QKI_WT cyclic_ss, dt=1fs, 25ns"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
log "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<all>}"
log "logdir=$LOGDIR  driver_log=$DRIVER_LOG"
log "WT input provenance: ${WT_REF_DIR}/_md_input (2QKI_WT.pdb + 2QKI_WT_renum.pdb; ADR-0009 WT band inputs)"
log "2 fresh dt1fs WT seeds: 411 | 433  (25ns each = 25M steps; dense DCD=${DCD_INTERVAL} step)"
log "Goal: WT contact engagement vs FIXED ADR-0009 bands -> WT(92+/-15)/patchoff(71+/-18)/gamma(12+/-11)"
log "----------------------------------------------------------------"
log "ANCHOR ASYMMETRY (load-bearing caveat): D5 d3refit anchors MTR res-4 backbone (k=100);"
log "  WT (--ncaa none) is FULLY UNRESTRAINED (run_restrained_md.py keys restraint off ncAA name,"
log "  so WT n_restrained=0). NOT apples-to-apples on the res-4 anchor; apples-to-apples on"
log "  system/topology/dt/length/metric. FREE-WT is the more conservative engagement test (see header)."
log "----------------------------------------------------------------"
log "PRE-REGISTER (3 outcomes, decided before data):"
log "  GREEN  (pipeline valid): pooled mean WT-like (~92) + high frac_engaged + sustained-to-25ns"
log "         -> protocol+metric report engagement -> D5 d3refit NO-GO is charge-specific"
log "  RED    (dt1fs/pipeline suspect): WT also disengages (gamma-band ~12 / frac_engaged <~0.20 / disengaged-by-25ns)"
log "         -> dt1fs/anchor-asymmetry artefact priority -> d3refit 'MTR-specific' conclusion SUSPENDED"
log "  AMBIGUOUS: seed-split or pooled frac_engaged 0.2~0.7 -> report distribution, grow n"
log "FORBIDDEN (R-11/R-18): magnitude / sign / lambda / binding-affinity / ddG claims; contact count is geometric only"
log "frozen params/MTR_gaff2_hybrid.xml & MTR_gaff2_layer3d_d3refit.xml: NOT read/touched (WT has no MTR)"
log "================================================================"

# ----- Pre-flight: WT input must exist; no MTR XML involved -----
if [ ! -f "${WT_REF_DIR}/_md_input/2QKI_WT.pdb" ]; then
    abort "WT input not found: ${WT_REF_DIR}/_md_input/2QKI_WT.pdb -- cannot mirror WT calib inputs."
fi
log "Pre-flight PASS: WT calib inputs present (no ncAA params / no MTR XML required for WT)."

# ----- per-seed snapshot extraction (PRIMARY metric pipeline, ADR-0009 SSOT) -----
# Identical SSOT extraction to d5_pose_recovery_production.sh: cluster_and_select(
# traj, 25, binder_chain="B") + CONECT-v55 save_snapshots ->
# snapshots_n25_postl387_patch_v2/. The reextract scripts are hardwired to the
# *_calib_* naming and cannot target a _d5posctrl_ dir, so the SSOT extractor
# functions are driven directly here (same functions, same parameters =
# band-comparable to both the WT and the D5 d3refit cohorts).
extract_snapshots_one () {
    local seedir="$1"
    local suffix="$2"
    local dcd="${seedir}/mdresult/2QKI_WT_restrained.dcd"
    local top="${seedir}/mdresult/2QKI_WT_final.pdb"
    local snapdir="${seedir}/${SNAP_SUBDIR}"

    if [ ! -f "$dcd" ]; then
        log "  [${suffix}] EXTRACT SKIP: no DCD at $dcd"
        return 1
    fi
    if [ ! -f "$top" ]; then
        log "  [${suffix}] EXTRACT SKIP: no topology (final.pdb) at $top"
        return 1
    fi

    log "  [${suffix}] extract n=${N_SNAPSHOTS} snapshots -> ${SNAP_SUBDIR}/ (SSOT cluster_and_select + CONECT-v55)"
    SEEDIR="$seedir" DCD="$dcd" TOP="$top" SNAPDIR="$snapdir" \
        NSNAP="$N_SNAPSHOTS" $PY - <<'PYEXTRACT'
import logging, os, sys, warnings
warnings.filterwarnings("ignore")
logging.getLogger("mdtraj").setLevel(logging.ERROR)
logging.getLogger("openmm").setLevel(logging.ERROR)
sys.path.insert(0, "/home/san/UPDD_proj")
import mdtraj as md
try:
    md.set_logger_level(logging.ERROR)
except Exception:
    pass
from utils.extract_snapshots import cluster_and_select, save_snapshots

dcd = os.environ["DCD"]
top = os.environ["TOP"]
snapdir = os.environ["SNAPDIR"]
nsnap = int(os.environ["NSNAP"])
basename = os.path.basename(top).replace("_final.pdb", "")

traj = md.load(dcd, top=top)
frames = cluster_and_select(traj, nsnap, binder_chain="B")
frames = sorted(set(int(f) for f in frames))[:nsnap]
os.makedirs(snapdir, exist_ok=True)
saved = save_snapshots(traj, frames, snapdir, basename)
n_conect = []
import glob as _glob
for p in sorted(_glob.glob(os.path.join(snapdir, "*_snap*_f*.pdb"))):
    n_conect.append(sum(1 for line in open(p) if line.startswith("CONECT")))
cmin = min(n_conect) if n_conect else 0
print("EXTRACT_OK n_saved=%d frames=%d CONECT_min=%d" % (len(saved), len(frames), cmin))
PYEXTRACT
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        log "  [${suffix}] EXTRACT FAIL (rc=$rc) -- see above"
        return 1
    fi

    # PBC-integrity gate on the saved snapshots (ADR-0009 identical:
    # verify_pbc_postpatch peptide-bond distance distribution). Non-fatal:
    # report status; BROKEN snapshots are flagged but do not abort the launcher.
    log "  [${suffix}] verify_pbc_postpatch on ${SNAP_SUBDIR}/ (peptide-bond integrity)"
    SNAPDIR="$snapdir" $PY - <<'PYPBC'
import os, sys
sys.path.insert(0, "/home/san/UPDD_proj/scripts")
from pathlib import Path
import verify_pbc_postpatch as vpp
snapdir = Path(os.environ["SNAPDIR"])
row = vpp.measure_system_dir(snapdir)
print("PBC %s: n=%s median=%s min/max=%s/%s" % (
    row.status, row.n_pdb, row.median_ang, row.min_ang, row.max_ang))
PYPBC
    return 0
}

# ----- one-seed launcher helper: (seed_int, dt_fs, steps, dcd_interval, suffix) -----
# Mirrors d5_pose_recovery_production.sh::run_one_seed (fresh dir, input mirror,
# post-MD snapshot extraction + PBC gate) but for the WT system: NO params
# manifest, NO ncAA XML install, NO manifest patch, --ncaa none (unrestrained).
run_one_seed () {
    local seed_int="$1"
    local dt_fs="$2"
    local steps="$3"
    local dcd_interval="$4"
    local suffix="$5"
    local seedir="${PROJ}/outputs/2QKI_WT_${suffix}"

    log "  [${suffix}] setup dir $seedir (seed=${seed_int}, dt=${dt_fs}fs, steps=${steps}, DCD=${dcd_interval} step)"
    mkdir -p "${seedir}/_md_input" "${seedir}/mdresult"

    # WT has no params/ — mirror only the _md_input PDBs (verbatim t1_n30 Stage B).
    if [ -d "${WT_REF_DIR}/_md_input" ]; then
        cp -rn "${WT_REF_DIR}/_md_input/." "${seedir}/_md_input/"
    fi
    if [ ! -f "${seedir}/_md_input/2QKI_WT.pdb" ]; then
        log "  [${suffix}] input mirror FAIL (no 2QKI_WT.pdb)"; return 1
    fi

    local ns
    ns=$(awk "BEGIN{print ${steps}*${dt_fs}/1e6}")
    if [ "$DRY_RUN" = "1" ]; then
        log "  [${suffix}] DRY_RUN -- would launch ${ns}ns WT MD (steps=${steps}, dt=${dt_fs}fs, DCD=${dcd_interval} step, --ncaa none unrestrained) + n=${N_SNAPSHOTS} snapshot extraction, not launching"
        return 0
    fi

    log "  [${suffix}] MD launch (${ns}ns physical, cyclic_ss, dt=${dt_fs}fs, dense DCD=${dcd_interval} step, seed=${seed_int}, --ncaa none = UNRESTRAINED WT)"
    conda activate md_simulation
    UPDD_MMGBSA_PLATFORM=CUDA UPDD_MD_DCD_INTERVAL="$dcd_interval" \
        $MD_PREFIX "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${seedir}/_md_input" \
            --outputdir "${seedir}/mdresult" \
            --steps "$steps" \
            --topology cyclic_ss \
            --binder_chain B \
            --graph_policy strict \
            --target_id 2QKI \
            --dt_fs "$dt_fs" \
            --platform CUDA \
            --seed "$seed_int" \
            > "${LOGDIR}/md_${suffix}.log" 2>&1
    local rc=$?
    conda deactivate

    log "  [${suffix}] MD finished rc=$rc (see ${LOGDIR}/md_${suffix}.log)"

    # honest per-seed crash check (by _EXPLODED_ marker, NOT by rc / final.pdb).
    local expl
    expl=$(ls "${seedir}/mdresult/"*_EXPLODED_* 2>/dev/null | head -1)
    if [ -n "$expl" ]; then
        log "  [${suffix}] CRASH: explosion marker $(basename "$expl"). Skipping snapshot extraction for this seed."
        return 0
    fi

    # PRIMARY metric: extract n=25 snapshots (only on a clean 25ns trajectory).
    extract_snapshots_one "$seedir" "$suffix"
    return 0
}

log "----------------------------------------------------------------"
log "Stage A: 2 fresh dt1fs WT seeds (sequential -- single GPU, concurrency forbidden)"
log "----------------------------------------------------------------"

# run 1: s411 dt=1fs, 25ns (25M step), DCD 2000 step = 2ps
run_one_seed 411 1.0 25000000 "$DCD_INTERVAL" "d5posctrl_s411_dt1fs"
# run 2: s433 dt=1fs, 25ns (25M step), DCD 2000 step = 2ps
run_one_seed 433 1.0 25000000 "$DCD_INTERVAL" "d5posctrl_s433_dt1fs"

# ----- honest tally: crash judged by _EXPLODED_ marker, NOT exit code / final.pdb -----
log "================================================================"
log "WT positive control Stage A done -- crash tally (by _EXPLODED_ marker):"
SEEDS="d5posctrl_s411_dt1fs d5posctrl_s433_dt1fs"
N_CRASH=0
N_COMPLETE=0
N_INCONCLUSIVE=0
N_SNAP_OK=0
for suffix in $SEEDS; do
    seedir="${PROJ}/outputs/2QKI_WT_${suffix}"
    expl=$(ls "${seedir}/mdresult/"*_EXPLODED_* 2>/dev/null | head -1)
    final="${seedir}/mdresult/2QKI_WT_final.pdb"
    n_snap=$(ls "${seedir}/${SNAP_SUBDIR}/"*_snap*_f*.pdb 2>/dev/null | wc -l)
    if [ -n "$expl" ]; then
        N_CRASH=$((N_CRASH + 1))
        log "  CRASH  ${suffix}: explosion marker $(basename "$expl")"
    elif [ -f "$final" ]; then
        N_COMPLETE=$((N_COMPLETE + 1))
        if [ "$n_snap" -ge "$N_SNAPSHOTS" ]; then
            N_SNAP_OK=$((N_SNAP_OK + 1))
            log "  COMPLETE ${suffix}: 25ns reached, no explosion marker, ${n_snap} snapshots extracted"
        else
            log "  COMPLETE ${suffix}: 25ns reached, no explosion marker, but only ${n_snap} snapshots (expected ${N_SNAPSHOTS}) -- check extraction"
        fi
    else
        N_INCONCLUSIVE=$((N_INCONCLUSIVE + 1))
        log "  INCONCLUSIVE ${suffix}: neither marker nor final.pdb (check ${LOGDIR}/md_${suffix}.log)"
    fi
done

log "----------------------------------------------------------------"
log "Tally: ${N_CRASH} crashed / ${N_COMPLETE} completed / ${N_INCONCLUSIVE} inconclusive (of 2 fresh WT seeds)"
log "Snapshots extracted (n>=${N_SNAPSHOTS}) on ${N_SNAP_OK} fresh seeds."
log "NOTE: crash-rate is NOT '0' even if N_CRASH=0 (Poisson rule-of-three); report 'N_CRASH/2' as a fact only."
log "----------------------------------------------------------------"
log "NEXT (analysis, run separately after this launcher):"
log "  PRIMARY (go/no-go DECISION):"
log "    $PY scripts/contact_engagement_analysis.py <stamp>"
log "    -> compare WT_dt1fs cohort pooled mean contacts / frac_engaged vs FIXED bands"
log "       WT 92+/-15 / patchoff 71+/-18 / gamma 12+/-11 AND vs the D5 d3refit_Cp4 cohort"
log "  DECISION: GREEN(WT-like + sustained-25ns => pipeline valid, d3refit NO-GO charge-specific) /"
log "    RED(WT-free disengaged => dt1fs/anchor-asymmetry suspect, d3refit 'MTR-specific' SUSPENDED) /"
log "    AMBIGUOUS(seed-split / frac 0.2~0.7). NO magnitude claim (R-11)."
log "  CAVEAT (R-18): WT ran UNRESTRAINED (anchor asymmetry, see header). FREE-WT engaged = conservative"
log "    support for the charge-specific reading; FREE-WT disengaged cannot cleanly isolate dt1fs (anchor"
log "    confound) -> an anchored-WT rerun would need an MD-engine change + a separate science-review sign-off (out of scope)."
log "================================================================"
