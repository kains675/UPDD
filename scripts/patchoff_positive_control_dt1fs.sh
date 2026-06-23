#!/bin/bash
# ============================================================================
# patch-off MTR positive control (D5 NO-GO specificity, PRIMARY) — dt=1fs
# pose-engagement reproduction on the SAME 2QKI_Cp4 (MTR res-4) system with the
# PATCH-OFF MTR force field installed instead of the D3-refit FF.
# (P3 #107, the PRIMARY positive control for the D5 d3refit-MTR NO-GO.)
#
# Purpose (specificity positive control ONLY): the D5 pose-recovery gate scored
# the D3-refit MTR force field as NO-GO (fully disengaged by 25ns, all seeds).
# That NO-GO can have two distinct roots that the D5 run alone cannot separate:
#   (a) sidechain-RESP-refit-specific — the refit MTR sidechain charges genuinely
#       fail to hold the pose;
#   (b) pipeline / dt=1fs artefact — the dt=1fs restrained 25ns protocol itself
#       (or the contact metric) disengages ANY binder, even a strong one.
# This launcher runs the patch-off MTR FF (the validated calib FF that engaged
# 72% / 99.7% in ADR-0009) through the SAME dt=1fs / 25ns protocol + the SAME
# contact-engagement metric, so the patch-off outcome discriminates:
#   patchoff engaged  => the protocol+metric CAN report engagement at dt=1fs =>
#                  the D5 NO-GO is NOT a dt=1fs/pipeline artefact => the d3refit
#                  NO-GO is sidechain-RESP-refit-specific.
#   patchoff disengage => the dt=1fs protocol disengages even the validated
#                  patch-off FF => the dt=1fs artefact hypothesis takes priority;
#                  the d3refit "refit-specific" conclusion is SUSPENDED (#106
#                  SUSPEND) pending a dt=1fs diagnosis.
# It is a GEOMETRIC go/no-go on the binding pose — NOT a free-energy / affinity
# computation.
#
# ------------------------------------------------------------------------
# !! ANCHOR SYMMETRY (why this is the PRIMARY control — load-bearing) !!
# ------------------------------------------------------------------------
# Both the D5 d3refit cohort AND this patch-off cohort carry the MTR ncAA at
# res-4, so run_restrained_md.py keys its res-4 backbone position restraint off
# the ncAA residue name (ctx.xml_res_name = "MTR") at k=100 kJ/mol/nm^2 in BOTH
# cohorts (--ncaa_label MTR --ncaa_code MTR). The restraint build block
# `if ctx.xml_res_name:` therefore fires identically for both -> n_restrained is
# the SAME res-4-anchored set in both arms. This is ANCHOR-SYMMETRIC, unlike the
# secondary WT control (--ncaa none -> FULLY UNRESTRAINED free MD). Consequently
# this patch-off control is fully apples-to-apples:
#   * SAME on: system family (2QKI_Cp4, MTR res-4), topology (cyclic_ss),
#     res-4 backbone anchor (k=100), timestep (dt=1fs), length (25ns), platform,
#     contact metric / bands.
#   * DIFFERS ONLY on: the MTR sidechain charges (patch-off CM=0.0487 vs d3refit
#     CM=0.109137 + the rest of the D3 RESP-A2 sidechain refit). The frozen
#     backbone+NE1 anchors (N=-0.4157, NE1=-0.3418, ...) are BIT-IDENTICAL in
#     both XMLs and Σq=0 in both. => single-variable isolation = sidechain refit.
# Because the variable is isolated to the sidechain RESP-A2 refit charge, a
# patch-off-engaged readout cleanly attributes the d3refit NO-GO to that refit.
# ------------------------------------------------------------------------
#
# 3-OUTCOME PRE-REGISTER (decided BEFORE any data, echoed to driver log):
#   GREEN  (pipeline valid): pooled mean contacts patch-off-like (~71±18, the
#          ADR-0009 patch-off band ~72 / 99.7% engaged) / high frac_engaged AND
#          sustained-to-25ns (no late-window disengagement). => the dt=1fs / 25ns
#          protocol + contact metric CAN report engagement for the validated
#          patch-off FF => the D5 d3refit NO-GO is NOT a dt=1fs/pipeline artefact
#          AND dt=1fs is OK => the d3refit NO-GO is sidechain-RESP-refit-specific
#          (CONFIRMED).
#   RED    (dt=1fs/pipeline suspect): patch-off also disengages (pooled mean in
#          the gamma band ~12 / frac_engaged < ~0.20 / disengaged by 25ns). =>
#          even the validated patch-off FF disengages under this protocol => the
#          dt=1fs artefact hypothesis takes priority; #106 is SUSPENDED pending a
#          dt=1fs diagnosis. (The d3refit "refit-specific" conclusion is held.)
#   AMBIGUOUS (intermediate): seeds split (one engaged, one disengaged) OR
#          pooled frac_engaged 0.2~0.7 OR time-series shows late disengagement.
#          Report the distribution, grow n, NO magnitude claim.
#
# FORBIDDEN (R-11/R-18): any magnitude / sign / lambda / binding-affinity /
# ddG claim. Contact count is a geometric observable only. A dt=1fs crash here is
# logged loudly and HALTED (not silently rescued).
#
# Metric (PRIMARY go/no-go DECISION, locked, identical to D5):
#   n=25 snapshot extraction (cluster_and_select + CONECT-v55 save_snapshots,
#   ADR-0009 SSOT pipeline) -> snapshots_n25_postl387_patch_v2/ ->
#   contact_engagement_analysis.py (locked cutoff 4.0A / n_engaged>=20 /
#   chainB<->chainA). Directly band-comparable to the FIXED WT/gamma/patchoff
#   ADR-0009 bands and to the D5 d3refit cohort.
#
# Runs (2 fresh dt=1fs patch-off seeds, sequential on a single GPU = 5070 Ti;
# concurrent execution forbidden = contention; V100 not attached):
#   run  seed  dt_fs  --steps    DCD_INTERVAL  length  outdir suffix
#   1    443   1.0    25000000   2000 (=2ps)    25ns    2QKI_Cp4_patchoff_d5posctrl_s443_dt1fs
#   2    457   1.0    25000000   2000 (=2ps)    25ns    2QKI_Cp4_patchoff_d5posctrl_s457_dt1fs
# n=2 = the patch-off FF is the validated strong-engagement FF; 2 fresh dt=1fs
# seeds are sufficient for the discriminating control (a clean engaged/disengaged
# readout). Seeds 443/457 are confirmed unused across all outputs (the secondary
# WT control already uses 411/433; these are kept distinct to avoid cohort
# cross-contamination). Length = 25ns (= 25M steps at dt=1fs) to MATCH the D5
# cohort length (band-comparability) and to capture slow / progressive
# disengagement (a short window would FALSELY score a late-disengaging seed as
# engaged).
#
# FF path (DECISIVE — read before launch): the patch-off MTR XML
# outputs/2QKI_Cp4_calib_s101/params/MTR_gaff2.xml (NE1=-0.3418 frozen,
# CM=0.0487, Σq=0) is installed read-only as the per-seed cohort
# params/MTR_gaff2.xml. The FROZEN refit XML params/MTR_gaff2_layer3d_d3refit.xml
# (CM=0.109137) is NEVER installed here — installing it would make this control
# identical to D5 and scientifically meaningless. The frozen
# params/MTR_gaff2_hybrid.xml is likewise never read or touched.
#
# Cp4 input provenance: mirrored from the SAME validated Cp4 calib cohort that
# fed the D5 cohort (outputs/2QKI_Cp4_calib_s101/_md_input: 2QKI_Cp4.pdb +
# 2QKI_Cp4_renum.pdb). The MD invocation reproduces the D5 d5pose launch verbatim
# (dt=1fs, 25ns, cyclic_ss, MTR res-4 anchor, UPDD_NCAA_AMBER14_PATCH=0) except
# for the installed MTR_gaff2.xml (patch-off, not d3refit) and the cohort dir.
#
# Dirs: NEW per-seed cohort dirs outputs/2QKI_Cp4_patchoff_d5posctrl_s{seed}_dt1fs.
# The existing 2QKI_Cp4_calib_* and 2QKI_Cp4_d3refit_* dirs are PRESERVED and
# never clobbered (R-7).
#
# Honest tally: each run is judged by the presence of an _EXPLODED_ marker (NOT
# by exit code "DONE" and NOT by _final.pdb, which is also written on a
# recoverable crash). Poisson rule-of-three on crash-rate (N_CRASH/2 is a fact,
# not "0%").
#
# Usage:
#   scripts/patchoff_positive_control_dt1fs.sh                        # 2 runs, sequential
#   CUDA_VISIBLE_DEVICES=0 scripts/patchoff_positive_control_dt1fs.sh # pin GPU
#   DRY_RUN=1 scripts/patchoff_positive_control_dt1fs.sh              # print plan, no MD
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

# ----- patch-off MTR FF path (DECISIVE: patch-off, NOT d3refit; read-only) -----
# Installed as the per-seed cohort params/MTR_gaff2.xml. The d3refit XML must
# NEVER be installed here (it would make this control identical to D5).
PATCHOFF_XML="${PROJ}/outputs/2QKI_Cp4_calib_s101/params/MTR_gaff2.xml"
HYDROGENS_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s101/params/MTR_hydrogens.xml"
MANIFEST_SRC_TEMPLATE="${PROJ}/outputs/2QKI_Cp4_calib_s101/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s101"
# the FROZEN d3refit XML — referenced ONLY to assert it is NOT what we install.
D3REFIT_XML="${PROJ}/params/MTR_gaff2_layer3d_d3refit.xml"

# ----- snapshot extraction SSOT (ADR-0009 / section 4.3.2) -----
SNAP_SUBDIR="snapshots_n25_postl387_patch_v2"
N_SNAPSHOTS=25

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/patchoff_positive_control_${TS_TAG}"
mkdir -p "$LOGDIR"
mkdir -p "${PROJ}/logs"
DISPATCH_LOG="${LOGDIR}/dispatch.log"
DRIVER_LOG="${PROJ}/logs/patchoff_positive_control.log"

log () {
    printf '[patchoff_posctrl %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG" "$DRIVER_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "patch-off MTR positive control (D5 NO-GO specificity, PRIMARY) START — 2QKI_Cp4 MTR res-4 cyclic_ss, dt=1fs, 25ns"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
log "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<all>}"
log "logdir=$LOGDIR  driver_log=$DRIVER_LOG"
log "patch-off MTR XML (INSTALLED): $PATCHOFF_XML"
log "d3refit MTR XML (NEVER installed here): $D3REFIT_XML"
log "Cp4 input provenance: ${CP4_REF_DIR}/_md_input (2QKI_Cp4.pdb + 2QKI_Cp4_renum.pdb; SAME as D5 cohort)"
log "2 fresh dt1fs patch-off seeds: 443 | 457  (25ns each = 25M steps; dense DCD=${DCD_INTERVAL} step)"
log "Goal: patch-off contact engagement vs FIXED ADR-0009 bands -> WT(92+/-15)/patchoff(71+/-18)/gamma(12+/-11)"
log "----------------------------------------------------------------"
log "ANCHOR SYMMETRY (why PRIMARY): MTR res-4 backbone anchored (k=100) in BOTH this cohort AND the"
log "  D5 d3refit cohort (--ncaa_label MTR -> ctx.xml_res_name=MTR -> restraint fires identically)."
log "  Apples-to-apples on system/topology/anchor/dt/length/metric; the ONLY variable is the MTR"
log "  sidechain RESP-A2 refit charge (patch-off CM=0.0487 vs d3refit CM=0.109137; backbone+NE1 anchors"
log "  bit-identical, Σq=0 in both) -> single-variable isolation = sidechain refit."
log "----------------------------------------------------------------"
log "PRE-REGISTER (3 outcomes, decided before data):"
log "  GREEN  (pipeline valid): pooled mean patch-off-like (~71/99.7%) + high frac_engaged + sustained-to-25ns"
log "         -> protocol+metric report engagement at dt=1fs -> D5 d3refit NO-GO is sidechain-RESP-refit-specific"
log "  RED    (dt1fs/pipeline suspect): patch-off also disengages (gamma-band ~12 / frac_engaged <~0.20 / disengaged-by-25ns)"
log "         -> dt1fs artefact priority -> #106 SUSPENDED pending dt1fs diagnosis"
log "  AMBIGUOUS: seed-split or pooled frac_engaged 0.2~0.7 -> report distribution, grow n"
log "FORBIDDEN (R-11/R-18): magnitude / sign / lambda / binding-affinity / ddG claims; contact count is geometric only"
log "frozen params/MTR_gaff2_hybrid.xml & MTR_gaff2_layer3d_d3refit.xml: NOT installed (patch-off XML is the FF here)"
log "================================================================"

# ----- Pre-flight: patch-off XML must exist + Sigma q == 0 + frozen anchors
# intact + CM == 0.0487 (patch-off, NOT d3refit's 0.109137). This both verifies
# integrity AND guards against accidentally installing the d3refit XML. -----
if [ ! -f "$PATCHOFF_XML" ]; then
    abort "patch-off XML not found: $PATCHOFF_XML -- expected the validated 2QKI_Cp4_calib_s101 MTR FF."
fi

$PY -c "
import sys
import xml.etree.ElementTree as ET
t = ET.parse('$PATCHOFF_XML')
res = next(r for r in t.findall('.//Residue') if r.get('name') == 'MTR')
charges = {a.get('name'): float(a.get('charge')) for a in res.findall('Atom')}
total = sum(charges.values())
print('Sigma_q (MTR) = %+.2e' % total)
assert abs(total) < 1e-5, 'Sigma_q not neutral: %r' % total
# frozen anchors (amber14SB Maier) -- bit-identical to the d3refit cohort.
expected = {'N': -0.4157, 'H': 0.2719, 'CA': -0.0275, 'HA': 0.1123,
            'C': 0.5973, 'O': -0.5679, 'NE1': -0.3418}
for name, exp in expected.items():
    got = charges.get(name)
    assert got is not None and abs(got - exp) < 1e-4, \
        '%s drift: got %r expected %r' % (name, got, exp)
# CM is the SINGLE-VARIABLE discriminator: patch-off CM must be 0.0487
# (NOT the d3refit 0.109137). This is the control's whole point.
cm = charges.get('CM')
assert cm is not None and abs(cm - 0.0487) < 1e-4, \
    'CM != patch-off 0.0487 (got %r). If this is ~0.1091 the d3refit XML was installed by mistake -- ABORT.' % cm
print('Pre-flight PASS: Sigma q=%+.2e, 7 frozen anchors intact, CM=%.4f (patch-off, NOT d3refit 0.109137)'
      % (total, cm))
" || abort "patch-off XML pre-flight verification failed"

# ----- per-seed snapshot extraction (PRIMARY metric pipeline, ADR-0009 SSOT) -----
# Mirrors phase_beta_repbsa_v2_reextract.py exactly: cluster_and_select(traj,
# 25, binder_chain="B") + CONECT-v55 save_snapshots -> snapshots_n25_postl387_patch_v2/.
# The reextract / verify_pbc_postpatch scripts are hardwired to the *_calib_*
# dir naming and cannot target a _patchoff_ dir, so the SSOT extractor functions
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
# Mirrors d5_pose_recovery_production.sh::run_one_seed (fresh dir, frozen XML
# install, manifest patch, charge echo, MTR res-4 anchor, post-MD snapshot
# extraction + PBC gate) but installs the PATCH-OFF MTR XML (NOT d3refit) and
# uses a fresh patchoff cohort dir.
run_one_seed () {
    local seed_int="$1"
    local dt_fs="$2"
    local steps="$3"
    local dcd_interval="$4"
    local suffix="$5"
    local seedir="${PROJ}/outputs/2QKI_Cp4_patchoff_${suffix}"

    log "  [${suffix}] setup dir $seedir (seed=${seed_int}, dt=${dt_fs}fs, steps=${steps}, DCD=${dcd_interval} step)"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && [ ! -f "${seedir}/_md_input/${f}" ] && \
            cp "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    # install PATCH-OFF MTR XML as cohort MTR_gaff2.xml (read-only source; the
    # patch-off calib XML is never modified -- cp only). NOT the d3refit XML.
    cp -f "$PATCHOFF_XML" "${seedir}/params/MTR_gaff2.xml"
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
m['charge_source'] = 'patch-off MTR FF (calib_s101 MTR_gaff2.xml: NE1/backbone frozen amber14SB, CM=0.0487 GAFF2-default-derived, Sigma_q=0); positive control for the D5 d3refit NO-GO'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
" || { log "  [${suffix}] manifest patch FAIL"; return 1; }

    local qN qNE1 qCM
    qN=$(grep -m1 '<Atom name="N" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    qNE1=$(grep -m1 '<Atom name="NE1" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    qCM=$(grep -m1 '<Atom name="CM" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    log "  [${suffix}] installed q_N=$qN (frozen) q_NE1=$qNE1 (frozen) q_CM=$qCM (patch-off; d3refit would be 0.109137)"

    local ns
    ns=$(awk "BEGIN{print ${steps}*${dt_fs}/1e6}")
    if [ "$DRY_RUN" = "1" ]; then
        log "  [${suffix}] DRY_RUN -- would launch ${ns}ns MD (steps=${steps}, dt=${dt_fs}fs, DCD=${dcd_interval} step, MTR res-4 anchor) + n=${N_SNAPSHOTS} snapshot extraction, not launching"
        return 0
    fi

    log "  [${suffix}] MD launch (${ns}ns physical, cyclic_ss, dt=${dt_fs}fs, dense DCD=${dcd_interval} step, seed=${seed_int}, --ncaa_label MTR = MTR res-4 anchored)"
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
log "Stage A: 2 fresh dt1fs patch-off seeds (sequential -- single GPU, concurrency forbidden)"
log "----------------------------------------------------------------"

# run 1: s443 dt=1fs, 25ns (25M step), DCD 2000 step = 2ps
run_one_seed 443 1.0 25000000 "$DCD_INTERVAL" "d5posctrl_s443_dt1fs"
# run 2: s457 dt=1fs, 25ns (25M step), DCD 2000 step = 2ps
run_one_seed 457 1.0 25000000 "$DCD_INTERVAL" "d5posctrl_s457_dt1fs"

# ----- honest tally: crash judged by _EXPLODED_ marker, NOT exit code / final.pdb -----
log "================================================================"
log "patch-off positive control Stage A done -- crash tally (by _EXPLODED_ marker):"
SEEDS="d5posctrl_s443_dt1fs d5posctrl_s457_dt1fs"
N_CRASH=0
N_COMPLETE=0
N_INCONCLUSIVE=0
N_SNAP_OK=0
for suffix in $SEEDS; do
    seedir="${PROJ}/outputs/2QKI_Cp4_patchoff_${suffix}"
    expl=$(ls "${seedir}/mdresult/"*_EXPLODED_* 2>/dev/null | head -1)
    final="${seedir}/mdresult/2QKI_Cp4_final.pdb"
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
log "Tally: ${N_CRASH} crashed / ${N_COMPLETE} completed / ${N_INCONCLUSIVE} inconclusive (of 2 fresh patch-off seeds)"
log "Snapshots extracted (n>=${N_SNAPSHOTS}) on ${N_SNAP_OK} fresh seeds."
log "NOTE: crash-rate is NOT '0' even if N_CRASH=0 (Poisson rule-of-three); report 'N_CRASH/2' as a fact only."
log "----------------------------------------------------------------"
log "NEXT (analysis, run separately after this launcher):"
log "  PRIMARY (go/no-go DECISION):"
log "    $PY scripts/contact_engagement_analysis.py <stamp>"
log "    -> compare patchoff_dt1fs cohort pooled mean contacts / frac_engaged vs FIXED bands"
log "       WT 92+/-15 / patchoff 71+/-18 / gamma 12+/-11 AND vs the D5 d3refit_Cp4 cohort"
log "  DECISION: GREEN(patch-off-like + sustained-25ns => pipeline valid AND dt1fs OK, d3refit NO-GO sidechain-RESP-refit-specific) /"
log "    RED(patch-off disengaged => dt1fs/pipeline suspect, #106 SUSPENDED) /"
log "    AMBIGUOUS(seed-split / frac 0.2~0.7). NO magnitude claim (R-11)."
log "  ANCHOR SYMMETRY (R-18): this PRIMARY control is anchor-symmetric to D5 (MTR res-4 anchored in both),"
log "    so the ONLY variable is the MTR sidechain RESP-A2 refit charge -> clean single-variable attribution."
log "================================================================"
