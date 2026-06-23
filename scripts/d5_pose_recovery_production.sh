#!/bin/bash
# ============================================================================
# D5 pose-recovery gate — D3-refit MTR FF binding-pose engagement test
# (P3 #106, follow-up to D4b localize/H1-H2 diagnostic).
#
# Purpose (pose-recovery gate ONLY): re-run 2QKI_Cp4 cyclic_ss bound MD with
# the D3-refit MTR force field (params/MTR_gaff2_layer3d_d3refit.xml) at a
# dt=1fs BASELINE timestep, then measure whether the binder-target contact
# engagement lands in the WT band (92+/-15), the patch-off band (71+/-18), or
# the gamma band (12+/-11) established in ADR-0009 / manuscript section 4.3.2.
# This is a geometric go/no-go DECISION on whether the refit FF retains the
# binding pose -- it is NOT a free-energy / binding-affinity computation.
#
# 3-OUTCOME PRE-REGISTER (decided BEFORE any data, echoed to driver log):
#   NO-GO  (most-likely, STOP): pooled mean contacts in the gamma band
#          (~12+/-11) / frac_engaged ~0.20 / fully disengaged by 25ns. The
#          refit FF does not restore the pose -> pose-sampling root confirmed,
#          no free-energy claim. NOT a D5 failure -- a publishable axis
#          separation (partial-engagement, manuscript X.B class).
#   PASS   (bonus, -> GATE-2 RBFE): pooled mean patchoff-like (~71 / 99%
#          engaged) OR WT-like (~92 / 100%) AND sustained-to-25ns (no
#          late-window disengagement in the dense time-series).
#   AMBIGUOUS (multi-basin): seeds split (some engaged, some disengaged) OR
#          pooled frac_engaged 0.2~0.7 OR time-series shows some seeds
#          disengaging late. Report the distribution, grow n, NO magnitude
#          claim, partial-engagement (X.B) framing.
#
# FORBIDDEN (R-11/R-18): any magnitude / sign / lambda / binding-affinity
# claim. Contact count is a geometric observable only. A dt=1fs crash here is
# the diagnostic signal that the H1 (MTR-improper, timestep-independent FF
# defect) hypothesis re-activates -- do NOT silently rescue; log it loudly and
# HALT escalation (dt=1fs baseline, NOT a rescue-after-2fs run).
#
# Metric (two roles, locked):
#   PRIMARY (go/no-go DECISION): n=25 snapshot extraction
#     (cluster_and_select + CONECT-v55 save_snapshots, ADR-0009 SSOT pipeline)
#     -> snapshots_n25_postl387_patch_v2/ -> contact_engagement_analysis.py
#     (locked cutoff 4.0A / n_engaged>=20 / chainB<->chainA). Directly
#     band-comparable to the FIXED WT/gamma/patchoff ADR-0009 bands.
#   SECONDARY (onset / sustained diagnostic): dense DCD (2ps/frame) PBC
#     min-image contact time-series (the D4b-style trajectory analysis). Used
#     for the WHEN / sustained-to-25ns narrative, NOT band-comparable to the
#     SSOT (it skips the cluster-select n=25 + imaging pipeline). Run
#     separately after this launcher (the dense DCD is the artefact it needs).
#
# Runs (3 fresh dt=1fs seeds, sequential on a single GPU = 5070 Ti; concurrent
# execution forbidden = contention; V100 not attached):
#   run  seed  dt_fs  --steps    DCD_INTERVAL  length  outdir suffix
#   1    311   1.0    25000000   2000 (=2ps)    25ns    _d5pose_s311_dt1fs
#   2    347   1.0    25000000   2000 (=2ps)    25ns    _d5pose_s347_dt1fs
#   3    379   1.0    25000000   2000 (=2ps)    25ns    _d5pose_s379_dt1fs
# Length = 25ns (= 25M steps at dt=1fs) to MATCH the ADR-0009 cohort length
# (band-comparability) and to capture the slow / progressive disengagement
# observed in the existing s211_dt1fs control (engaged early, fully disengaged
# only by ~25ns; a short window would FALSELY score it engaged).
#
# REUSE: the existing dt=1fs control outputs/2QKI_Cp4_d3refit_localize_s211_dt1fs
# (18GB dense DCD, 25ns) is the 4th member of the n>=4 ensemble. This launcher
# does NOT touch it -- s211 is read directly by the analysis step.
#
# FF path: the FROZEN refit XML params/MTR_gaff2_layer3d_d3refit.xml is
# installed read-only as the per-seed cohort params/MTR_gaff2.xml exactly like
# d4b_localize_dense_dcd.sh / d4_stability_smoke.sh (anti-fragmentation). Only
# dt and the cohort dir differ from the calib cohort -- the FF is identical
# (the precondition that makes the contact metric band-comparable).
#
# Dirs: NEW per-seed cohort dirs outputs/2QKI_Cp4_d3refit_d5pose_s{seed}_dt1fs.
# The previously created _localize_ / _calib_ dirs are PRESERVED and never
# clobbered (R-7).
#
# Honest tally: each run is judged by the presence of an _EXPLODED_ marker (NOT
# by exit code "DONE" and NOT by _final.pdb, which is also written on a
# recoverable crash). A dt=1fs crash is the H1 re-activation signal.
#
# Usage:
#   scripts/d5_pose_recovery_production.sh                        # 3 runs, sequential
#   CUDA_VISIBLE_DEVICES=0 scripts/d5_pose_recovery_production.sh # pin GPU
#   DRY_RUN=1 scripts/d5_pose_recovery_production.sh              # print plan, no MD
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ" || { echo "cannot cd $PROJ"; exit 2; }

source /home/san/miniconda3/etc/profile.d/conda.sh

# ----- knobs (env-overridable) -----
DRY_RUN="${DRY_RUN:-0}"
# Dense DCD interval for the SECONDARY time-series (2ps/frame, matches the
# existing s211_dt1fs control). Override via UPDD_MD_DCD_INTERVAL if needed.
DCD_INTERVAL="${UPDD_MD_DCD_INTERVAL:-2000}"

# ----- CPU affinity (d4_stability_smoke.sh / d4b logic) -----
PHYS_CORES=$(lscpu -p 2>/dev/null | grep -v "^#" | awk -F, '{print $2}' | sort -nu | wc -l)
LOGICAL_CORES=$(nproc)
if [ "${PHYS_CORES:-0}" -ge 4 ] && [ "${PHYS_CORES:-0}" -le 8 ]; then
    MD_AFFINITY="0,${PHYS_CORES}"
    MD_PREFIX="taskset -c ${MD_AFFINITY}"
else
    MD_PREFIX=""
fi

# ----- refit FF path (frozen, read-only) -----
REFIT_XML="${PROJ}/params/MTR_gaff2_layer3d_d3refit.xml"
HYDROGENS_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_hydrogens.xml"
MANIFEST_SRC_TEMPLATE="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"

# ----- snapshot extraction SSOT (ADR-0009 / section 4.3.2) -----
SNAP_SUBDIR="snapshots_n25_postl387_patch_v2"
N_SNAPSHOTS=25

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/d5_pose_recovery_${TS_TAG}"
mkdir -p "$LOGDIR"
mkdir -p "${PROJ}/logs"
DISPATCH_LOG="${LOGDIR}/dispatch.log"
DRIVER_LOG="${PROJ}/logs/d5_pose_recovery.log"

log () {
    printf '[d5_pose %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG" "$DRIVER_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "D5 pose-recovery gate START (D3 refit MTR FF, 2QKI_Cp4 cyclic_ss, dt=1fs baseline)"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
log "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<all>}"
log "logdir=$LOGDIR  driver_log=$DRIVER_LOG"
log "refit XML (frozen): $REFIT_XML"
log "3 fresh dt1fs seeds: 311 | 347 | 379  (25ns each = 25M steps; dense DCD=${DCD_INTERVAL} step)"
log "REUSE (not launched here): outputs/2QKI_Cp4_d3refit_localize_s211_dt1fs (existing 25ns dt1fs control) = ensemble member 4 (n>=4)"
log "Goal: contact engagement vs FIXED ADR-0009 bands -> WT(92+/-15)/patchoff(71+/-18)/gamma(12+/-11)"
log "----------------------------------------------------------------"
log "PRE-REGISTER (3 outcomes, decided before data):"
log "  NO-GO  (most-likely): pooled mean ~gamma-band (~12) / frac_engaged ~0.20 / disengaged by 25ns -> refit does not restore pose; NOT a failure (publishable X.B axis separation)"
log "  PASS   (bonus):       pooled mean patchoff-like (~71/99%) or WT-like (~92/100%) AND sustained-to-25ns -> GATE-2 RBFE"
log "  AMBIGUOUS (multi-basin): seed-split or pooled frac_engaged 0.2~0.7 -> report distribution, grow n, partial-engagement (X.B)"
log "FORBIDDEN (R-11/R-18): magnitude / sign / lambda / binding-affinity claims; contact count is a geometric observable only"
log "dt=1fs is a BASELINE (NOT a rescue-after-2fs run). dt=1fs crash -> H1 (MTR-improper) re-activation signal -> HALT escalation, no silent rescue"
log "================================================================"

# ----- Pre-flight: refit XML must exist + Sigma q == 0 + frozen anchors intact
# + NE1-CM r_eq + N-methyl rotor (identical verification to d4b/d4_stability_smoke) -----
if [ ! -f "$REFIT_XML" ]; then
    abort "refit XML not found: $REFIT_XML -- run scripts/d4_build_refit_xml.py first."
fi

$PY -c "
import sys
import xml.etree.ElementTree as ET
t = ET.parse('$REFIT_XML')
res = next(r for r in t.findall('.//Residue') if r.get('name') == 'MTR')
charges = {a.get('name'): float(a.get('charge')) for a in res.findall('Atom')}
total = sum(charges.values())
print('Sigma_q (MTR) = %+.2e' % total)
assert abs(total) < 1e-5, 'Sigma_q not neutral: %r' % total
# frozen anchors (amber14SB Maier)
expected = {'N': -0.4157, 'H': 0.2719, 'CA': -0.0275, 'HA': 0.1123,
            'C': 0.5973, 'O': -0.5679, 'NE1': -0.3418}
for name, exp in expected.items():
    got = charges.get(name)
    assert got is not None and abs(got - exp) < 1e-4, \
        '%s drift: got %r expected %r' % (name, got, exp)
# NE1-CM bond r_eq + new rotor Proper present
bond = next((b for b in t.findall('.//HarmonicBondForce/Bond')
             if b.get('type1') == 'protein-CT' and b.get('type2') == 'protein-NA'), None)
assert bond is not None and abs(float(bond.get('length')) - 0.14482) < 1e-6, \
    'NE1-CM bond r_eq != 0.14482: %r' % (bond.get('length') if bond is not None else None)
rotor = next((p for p in t.findall('.//PeriodicTorsionForce/Proper')
              if (p.get('type1'), p.get('type2'), p.get('type3'), p.get('type4'))
              == ('protein-CW', 'protein-NA', 'protein-CT', 'protein-HC')), None)
assert rotor is not None, 'N-methyl rotor Proper missing'
assert abs(float(rotor.get('k1')) - 0.461282) < 1e-3, \
    'rotor k1 != 0.461282: %r' % rotor.get('k1')
assert int(rotor.get('periodicity1')) == 3 and abs(float(rotor.get('phase1'))) < 1e-9
print('Pre-flight PASS: Sigma q=%+.2e, 7 frozen anchors intact, NE1-CM r_eq=0.14482, rotor k1=%.6f'
      % (total, float(rotor.get('k1'))))
" || abort "refit XML pre-flight verification failed"

# ----- per-seed snapshot extraction (PRIMARY metric pipeline, ADR-0009 SSOT) -----
# Mirrors phase_beta_repbsa_v2_reextract.py exactly: cluster_and_select(traj,
# 25, binder_chain="B") + CONECT-v55 save_snapshots -> snapshots_n25_postl387_patch_v2/.
# The reextract / verify_pbc_postpatch scripts are hardwired to the *_calib_*
# dir naming and cannot target a _d5pose_ dir, so the SSOT extractor functions
# are driven directly here (same functions, same parameters = band-comparable).
extract_snapshots_one () {
    local seedir="$1"
    local suffix="$2"
    local dcd="${seedir}/mdresult/2QKI_Cp4_restrained.dcd"
    local top="${seedir}/mdresult/2QKI_Cp4_final.pdb"
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
    # verify_pbc_postpatch peptide-bond distance distribution). This is the
    # downstream consumer artefact that the contact metric actually reads.
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
# non-fatal: report status; BROKEN snapshots are flagged but do not abort the
# launcher (the analysis step records the integrity status alongside contacts).
PYPBC
    return 0
}

# ----- one-seed launcher helper: (seed_int, dt_fs, steps, dcd_interval, suffix) -----
# Mirrors d4b_localize_dense_dcd.sh::run_one_arm per-seed setup (fresh dir,
# frozen XML install, manifest patch, charge echo) with a fresh d5pose cohort
# dir and a post-MD snapshot-extraction + PBC-integrity step.
run_one_seed () {
    local seed_int="$1"
    local dt_fs="$2"
    local steps="$3"
    local dcd_interval="$4"
    local suffix="$5"
    local seedir="${PROJ}/outputs/2QKI_Cp4_d3refit_${suffix}"

    log "  [${suffix}] setup dir $seedir (seed=${seed_int}, dt=${dt_fs}fs, steps=${steps}, DCD=${dcd_interval} step)"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && [ ! -f "${seedir}/_md_input/${f}" ] && \
            cp "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    # install frozen refit XML as cohort MTR_gaff2.xml (read-only source; the
    # FROZEN template + refit source are never modified -- cp only).
    cp -f "$REFIT_XML" "${seedir}/params/MTR_gaff2.xml"
    [ -f "$HYDROGENS_SRC" ] && [ ! -f "${seedir}/params/$(basename "$HYDROGENS_SRC")" ] && \
        cp "$HYDROGENS_SRC" "${seedir}/params/"
    [ -f "$MANIFEST_SRC_TEMPLATE" ] && [ ! -f "${seedir}/params/$(basename "$MANIFEST_SRC_TEMPLATE")" ] && \
        cp "$MANIFEST_SRC_TEMPLATE" "${seedir}/params/"

    $PY -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
m['charge_source'] = 'D3 RESP-A2 (NE1/backbone frozen amber14SB) + NA-CT r_eq 1.4482A + N-methyl rotor V3; refit MTR_gaff2_layer3d_d3refit.xml'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
" || { log "  [${suffix}] manifest patch FAIL"; return 1; }

    local qN qNE1 qCM
    qN=$(grep -m1 '<Atom name="N" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    qNE1=$(grep -m1 '<Atom name="NE1" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    qCM=$(grep -m1 '<Atom name="CM" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    log "  [${suffix}] refit q_N=$qN (frozen) q_NE1=$qNE1 (frozen) q_CM=$qCM (RESP-A2)"

    local ns
    ns=$(awk "BEGIN{print ${steps}*${dt_fs}/1e6}")
    if [ "$DRY_RUN" = "1" ]; then
        log "  [${suffix}] DRY_RUN -- would launch ${ns}ns MD (steps=${steps}, dt=${dt_fs}fs, DCD=${dcd_interval} step) + n=${N_SNAPSHOTS} snapshot extraction, not launching"
        return 0
    fi

    log "  [${suffix}] MD launch (${ns}ns physical, cyclic_ss, dt=${dt_fs}fs baseline, dense DCD=${dcd_interval} step, seed=${seed_int})"
    conda activate md_simulation
    UPDD_NCAA_AMBER14_PATCH=0 UPDD_MMGBSA_PLATFORM=CUDA UPDD_MD_DCD_INTERVAL="$dcd_interval" \
        $MD_PREFIX "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${seedir}/_md_input" \
            --outputdir "${seedir}/mdresult" \
            --params_manifest "${seedir}/params/MTR_params_manifest.json" \
            --steps "$steps" \
            --topology cyclic_ss \
            --binder_chain B \
            --graph_policy strict \
            --target_id 2QKI \
            --dt_fs "$dt_fs" \
            --platform CUDA \
            --seed "$seed_int" \
            --ncaa_label MTR \
            --ncaa_code MTR \
            > "${LOGDIR}/md_${suffix}.log" 2>&1
    local rc=$?
    conda deactivate

    log "  [${suffix}] MD finished rc=$rc (see ${LOGDIR}/md_${suffix}.log)"

    # honest per-seed crash check (by _EXPLODED_ marker, NOT by rc / final.pdb).
    # On a recoverable crash both the marker AND _final.pdb are written, so a
    # marker presence is authoritative. dt=1fs crash = H1 re-activation signal.
    local expl
    expl=$(ls "${seedir}/mdresult/"*_EXPLODED_* 2>/dev/null | head -1)
    if [ -n "$expl" ]; then
        log "  [${suffix}] CRASH: explosion marker $(basename "$expl") (dt=1fs baseline -> H1 re-activation signal). Skipping snapshot extraction for this seed."
        return 0
    fi

    # PRIMARY metric: extract n=25 snapshots (only on a clean 25ns trajectory).
    extract_snapshots_one "$seedir" "$suffix"
    return 0
}

log "----------------------------------------------------------------"
log "Stage A: 3 fresh dt1fs seeds (sequential -- single GPU, concurrency forbidden)"
log "----------------------------------------------------------------"

# run 1: s311 dt=1fs, 25ns (25M step), DCD 2000 step = 2ps
run_one_seed 311 1.0 25000000 "$DCD_INTERVAL" "d5pose_s311_dt1fs"
# run 2: s347 dt=1fs, 25ns (25M step), DCD 2000 step = 2ps
run_one_seed 347 1.0 25000000 "$DCD_INTERVAL" "d5pose_s347_dt1fs"
# run 3: s379 dt=1fs, 25ns (25M step), DCD 2000 step = 2ps
run_one_seed 379 1.0 25000000 "$DCD_INTERVAL" "d5pose_s379_dt1fs"

# ----- honest tally: crash judged by _EXPLODED_ marker, NOT exit code / final.pdb -----
log "================================================================"
log "D5 pose-recovery Stage A done -- crash tally (by _EXPLODED_ marker):"
SEEDS="d5pose_s311_dt1fs d5pose_s347_dt1fs d5pose_s379_dt1fs"
N_CRASH=0
N_COMPLETE=0
N_INCONCLUSIVE=0
N_SNAP_OK=0
for suffix in $SEEDS; do
    seedir="${PROJ}/outputs/2QKI_Cp4_d3refit_${suffix}"
    expl=$(ls "${seedir}/mdresult/"*_EXPLODED_* 2>/dev/null | head -1)
    final="${seedir}/mdresult/2QKI_Cp4_final.pdb"
    n_snap=$(ls "${seedir}/${SNAP_SUBDIR}/"*_snap*_f*.pdb 2>/dev/null | wc -l)
    if [ -n "$expl" ]; then
        N_CRASH=$((N_CRASH + 1))
        log "  CRASH  ${suffix}: explosion marker $(basename "$expl") (dt=1fs -> H1 re-activation)"
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

# ----- reuse member: extract s211_dt1fs snapshots (idempotent, R-7 safe) -----
# The existing dt=1fs control outputs/2QKI_Cp4_d3refit_localize_s211_dt1fs is
# ensemble member 4 (n>=4). It has the 25ns dense DCD but its snapshots were
# never extracted. Extract them the SAME SSOT way so the d3refit_Cp4 cohort
# (glob 2QKI_Cp4_d3refit_*_dt1fs in contact_engagement_analysis.py) includes it.
# This only WRITES the (currently absent) snapshot subdir -- the trajectory and
# all other s211 artefacts are read-only (never clobbered).
S211_DIR="${PROJ}/outputs/2QKI_Cp4_d3refit_localize_s211_dt1fs"
if [ -d "$S211_DIR" ]; then
    s211_snaps=$(ls "${S211_DIR}/${SNAP_SUBDIR}/"*_snap*_f*.pdb 2>/dev/null | wc -l)
    if [ "$s211_snaps" -ge "$N_SNAPSHOTS" ]; then
        log "  [reuse s211_dt1fs] ${s211_snaps} snapshots already present -- skip extraction"
    elif [ "$DRY_RUN" = "1" ]; then
        log "  [reuse s211_dt1fs] DRY_RUN -- would extract n=${N_SNAPSHOTS} snapshots (none present), not extracting"
    else
        log "  [reuse s211_dt1fs] extracting n=${N_SNAPSHOTS} snapshots (ensemble member 4)"
        extract_snapshots_one "$S211_DIR" "reuse_s211_dt1fs"
    fi
else
    log "  [reuse s211_dt1fs] WARN dir absent ($S211_DIR) -- d3refit_Cp4 cohort will be fresh-seeds only"
fi

log "----------------------------------------------------------------"
log "Tally: ${N_CRASH} crashed / ${N_COMPLETE} completed / ${N_INCONCLUSIVE} inconclusive (of 3 fresh seeds; +1 reused s211_dt1fs = n>=4 total)"
log "Snapshots extracted (n>=${N_SNAPSHOTS}) on ${N_SNAP_OK} fresh seeds."
log "NOTE: crash-rate is NOT '0' even if N_CRASH=0 (Poisson rule-of-three); report 'N_CRASH/3' as a fact only."
log "----------------------------------------------------------------"
log "NEXT (analysis, run separately after this launcher):"
log "  PRIMARY (go/no-go DECISION): extract s211_dt1fs snapshots the same way if absent, then:"
log "    $PY scripts/contact_engagement_analysis.py <stamp>"
log "    -> compare d3refit_Cp4 pooled mean contacts / frac_engaged vs FIXED bands"
log "       WT 92+/-15 / patchoff 71+/-18 / gamma 12+/-11 (PRIMARY DECISION)"
log "  SECONDARY (onset / sustained-to-25ns diagnostic): dense-DCD PBC min-image"
log "    contact time-series (D4b-style) on each 2ps/frame DCD -- NOT band-comparable,"
log "    used for the WHEN narrative + sustained-to-25ns false-PASS guard."
log "  DECISION: NO-GO(gamma-band/disengaged-by-25ns, most-likely) / PASS(patchoff-like"
log "    AND sustained-25ns) / AMBIGUOUS(seed-split / frac 0.2~0.7). NO magnitude claim (R-11)."
log "  NO-GO is NOT a failure (R-18): pose-root confirmed = publishable axis separation (X.B)."
log "================================================================"
