#!/usr/bin/env bash
# =============================================================================
# Track B in-place residue-4 RBFE — BOUND leg OOM recovery (reboot-survivable).
#
# Purpose: after a host reboot (RAM upgrade attempt), recover the cp4 bound
# leg dminus replicates that were lost to a dplus->dminus transition OOM on the
# VM, WITHOUT destroying the 14h cp4 dplus .out files that already completed.
#
# This script is IDEMPOTENT and SAFE BY DEFAULT:
#   - With NO flag it only INSPECTS + PLANS (dry-run): host-RAM scenario
#     detection, data-integrity check, dplus backup, and the exact rerun
#     command it WOULD run. It changes NOTHING on the VM and launches NOTHING.
#   - Only with --go does it (1) reallocate VM memory per scenario and
#     (2) launch the cp4 dminus rerun (pool, --no-archive-existing) detached.
#
# The companion human runbook (zero-context readable) is:
#   outputs/_trackb/REBOOT_RECOVERY_RUNBOOK.md
#
# Safety contract for the irreplaceable cp4 dplus data (400 cycles x 6 windows
# x 3 reps):
#   1. The rerun uses pool mode with --directions dminus only + --endpoints cp4
#      only. The dplus subdir is NEVER touched (only the dminus subdir of each
#      rep is rewritten by run_one_direction).
#   2. The rerun passes --no-archive-existing so the launcher does NOT
#      shutil.move the whole rep dir (which WOULD include dplus) to _archive/.
#   3. BEFORE any rerun this script copies (R-7) every cp4 dplus .out into a
#      separate backup dir on the VM. The copy is verified before --go proceeds.
#   4. AFTER the rerun launches, the operator (and STEP 6 of the runbook) must
#      re-verify cp4 dplus is still 400x6x3.
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Config (matches the original bound worker args VERBATIM so the ladder is
# identical — do NOT change these without re-validating; ladder identity is
# required for the dplus already on disk to remain mergeable with the new
# dminus).
# ---------------------------------------------------------------------------
VM_SSH="san@192.168.122.155"
VM_DOMAIN="v100-vm"
VM_PY="/home/san/miniconda3/envs/atm/bin/python3.11"
PROJ="/home/san/UPDD_proj"
OUT_ROOT="outputs/_trackb/inplace_rbfe_prod_bound"        # relative to PROJ (VM cwd)
OUT_ROOT_ABS="${PROJ}/${OUT_ROOT}"
LAUNCHER="scripts/trackb_inplace_rbfe_production.py"

# Ladder / run params (original bound worker args).
N_WINDOWS_HALF=6
SOFTCORE_BAND=2
N_CYCLES=400
MD_STEPS_PER_CYCLE=250
PLATFORM="CUDA"
TIMESTEP_FS=1.0
MINIMIZE_ITERS=500
BACKWARD_EQUIL_STEPS=500
GENUINE_DECOUPLE_NM=1.2
BINDER_CHAIN="B"
DEVICE_INDEX=0
MINTIMEID=100

# Recovery target = cp4 dminus (the lost direction). dplus is reused in place.
RECOVER_ENDPOINT="cp4"
RECOVER_DIRECTION="dminus"
RECOVER_SEEDS="s7,s101,s127"
EXPECTED_CYCLES=400          # a direction r0 .out is COMPLETE only at >= this
EXPECTED_WINDOWS=6           # r0..r5 per direction (n-windows-half)
N_REPS=3                     # rep0,rep1,rep2

# VM memory targets (KiB) + concurrency per scenario.
SCEN_A_VM_KIB=50331648       # 48 GiB
SCEN_A_CONCURRENCY=3
SCEN_B_VM_KIB=41943040       # 40 GiB
SCEN_B_CONCURRENCY=2
HOST_SCENARIO_A_MIN_GIB=88   # host total RAM (GiB) threshold for scenario A

DPLUS_BACKUP_DIR="${OUT_ROOT}/_dplus_backup_${RECOVER_ENDPOINT}"   # VM-relative (R-7)

GO=0
for arg in "$@"; do
    case "$arg" in
        --go) GO=1 ;;
        -h|--help)
            echo "Usage: $0 [--go]"
            echo "  (no flag) inspect + plan only (dry-run, changes nothing)"
            echo "  --go      reallocate VM memory + launch cp4 dminus rerun"
            exit 0 ;;
        *) echo "ERROR: unknown arg '$arg' (use --go or no arg)"; exit 2 ;;
    esac
done

log()  { printf '[recovery] %s\n' "$*"; }
warn() { printf '[recovery][WARN] %s\n' "$*" >&2; }
die()  { printf '[recovery][FATAL] %s\n' "$*" >&2; exit 1; }

ssh_vm() { ssh -o ConnectTimeout=15 -o BatchMode=yes "$VM_SSH" "$@"; }

# ---------------------------------------------------------------------------
# STEP 0 — VM up?
# ---------------------------------------------------------------------------
log "STEP 0: VM domain state"
VM_STATE="$(virsh domstate "$VM_DOMAIN" 2>/dev/null || echo unknown)"
log "  virsh domstate ${VM_DOMAIN} = ${VM_STATE}"
if [ "$VM_STATE" != "running" ]; then
    warn "VM not running. With --go this script will NOT auto-start it (memory"
    warn "must be set while shutoff first). See runbook STEP 0/STEP 2."
fi

# ---------------------------------------------------------------------------
# STEP 1 — host RAM -> scenario
# ---------------------------------------------------------------------------
log "STEP 1: host RAM scenario detection"
HOST_TOTAL_MIB="$(free -m | awk 'NR==2{print $2}')"
HOST_TOTAL_GIB=$(( HOST_TOTAL_MIB / 1024 ))
log "  host total RAM = ${HOST_TOTAL_GIB} GiB (${HOST_TOTAL_MIB} MiB)"
if [ "$HOST_TOTAL_GIB" -ge "$HOST_SCENARIO_A_MIN_GIB" ]; then
    SCENARIO="A"
    VM_KIB="$SCEN_A_VM_KIB"
    CONCURRENCY="$SCEN_A_CONCURRENCY"
    log "  -> SCENARIO A (host >= ${HOST_SCENARIO_A_MIN_GIB} GiB): VM 48 GiB, concurrency 3"
else
    SCENARIO="B"
    VM_KIB="$SCEN_B_VM_KIB"
    CONCURRENCY="$SCEN_B_CONCURRENCY"
    log "  -> SCENARIO B (host < ${HOST_SCENARIO_A_MIN_GIB} GiB): VM 40 GiB, concurrency 2"
fi
VM_GIB=$(( VM_KIB / 1024 / 1024 ))

# ---------------------------------------------------------------------------
# STEP 3 — data integrity (run before backup so we know what we are protecting).
# A direction is COMPLETE for a rep only if r0..r{W-1} each have >= EXPECTED_CYCLES
# lines. EMPTY/short .out files (the OOM trap: files present but 0 bytes) count
# as INCOMPLETE.
# ---------------------------------------------------------------------------
log "STEP 3: data-integrity check (line counts, not file presence)"

# Remote probe BODY is a fully-quoted heredoc (zero local expansion). The host
# values it needs are injected as a remote env prelude (PROBE_ENV) prepended to
# the body, so nothing in the body is mangled by local shell expansion.
PROBE_ENV="OUT='${OUT_ROOT_ABS}'; EXP=${EXPECTED_CYCLES}; W=${EXPECTED_WINDOWS}; NREPS=${N_REPS};"
REMOTE_PROBE_BODY=$(cat <<'REMOTE'
probe_dir() {  # ep leg rep tag
  ep=$1; leg=$2; rep=$3; tag=$4
  ok=1
  w=0
  while [ "$w" -lt "$W" ]; do
    f="${OUT}/${ep}/${leg}/rep${rep}/${tag}/r${w}/trackb_${tag}.out"
    n=0
    [ -f "$f" ] && n=$(wc -l < "$f" 2>/dev/null || echo 0)
    [ "$n" -lt "$EXP" ] && ok=0
    w=$((w+1))
  done
  if [ "$ok" -eq 1 ]; then echo "${ep} ${leg} rep${rep} ${tag} complete";
  else echo "${ep} ${leg} rep${rep} ${tag} INCOMPLETE"; fi
}
r=0
while [ "$r" -lt "$NREPS" ]; do
  probe_dir cp4 bound "$r" dplus
  probe_dir cp4 bound "$r" dminus
  probe_dir wt  bound "$r" dplus
  probe_dir wt  bound "$r" dminus
  r=$((r+1))
done
REMOTE
)
REMOTE_PROBE="${PROBE_ENV}
${REMOTE_PROBE_BODY}"

if [ "$VM_STATE" = "running" ]; then
    INTEGRITY="$(ssh_vm "bash -s" <<<"$REMOTE_PROBE" || true)"
    echo "$INTEGRITY" | sed 's/^/  /'

    # cp4 dplus MUST be complete on all reps or we abort (cannot lose it).
    CP4_DPLUS_BAD="$(echo "$INTEGRITY" | grep -c '^cp4 bound rep[0-9] dplus INCOMPLETE' || true)"
    if [ "${CP4_DPLUS_BAD:-0}" -ne 0 ]; then
        die "cp4 dplus is INCOMPLETE on ${CP4_DPLUS_BAD} rep(s). STOP. Do NOT rerun; investigate (the irreplaceable dplus is the whole point of recovery)."
    fi
    log "  cp4 dplus complete on all ${N_REPS} reps (safe to proceed)."

    # Reps that still need a rerun (any direction INCOMPLETE). For info; the
    # cp4 dminus rerun is the explicit target. wt INCOMPLETE reps are listed so
    # the operator can extend the rerun (see runbook STEP 3 note).
    log "  reps with an INCOMPLETE direction (rerun candidates):"
    echo "$INTEGRITY" | grep 'INCOMPLETE' | sed 's/^/    /' || log "    (none — nothing to recover?)"
else
    warn "VM not running — cannot probe data. Start the VM (STEP 0) then re-run this script."
    INTEGRITY=""
fi

# ---------------------------------------------------------------------------
# STEP 4 — cp4 dplus backup (R-7 belt-and-suspenders). Copy (not move) every
# cp4 dplus .out into a separate VM dir. Verified before any rerun.
# ---------------------------------------------------------------------------
log "STEP 4: cp4 dplus backup plan -> ${DPLUS_BACKUP_DIR} (VM-side, copy/R-7)"
BACKUP_ENV="OUT='${OUT_ROOT_ABS}'; BK='${PROJ}/${DPLUS_BACKUP_DIR}'; NREPS=${N_REPS};"
BACKUP_BODY=$(cat <<'REMOTE'
set -e
mkdir -p "$BK"
copied=0
r=0
while [ "$r" -lt "$NREPS" ]; do
  src="${OUT}/cp4/bound/rep${r}/dplus"
  dst="${BK}/rep${r}/dplus"
  if [ -d "$src" ]; then
    mkdir -p "$dst"
    cp -an "$src/." "$dst/" 2>/dev/null || cp -rn "$src/." "$dst/"
    n=$(find "$dst" -name 'trackb_dplus.out' | wc -l)
    copied=$((copied+n))
  fi
  r=$((r+1))
done
echo "backup_out_files=$copied"
REMOTE
)
BACKUP_CMD="${BACKUP_ENV}
${BACKUP_BODY}"

# ---------------------------------------------------------------------------
# STEP 5 — the rerun command (cp4 dminus only, pool, no-archive-existing).
# Built here so dry-run can PRINT it exactly and --go can RUN it verbatim.
# ---------------------------------------------------------------------------
RERUN_LOG="${OUT_ROOT}/cp4_dminus_recovery.$(date +%Y%m%d_%H%M%S).log"
RERUN_CMD=( "$VM_PY" "$LAUNCHER"
    --leg bound
    --endpoints "$RECOVER_ENDPOINT"
    --directions "$RECOVER_DIRECTION"
    --seeds "$RECOVER_SEEDS"
    --pool
    --max-concurrent "$CONCURRENCY"
    --no-archive-existing
    --n-windows-half "$N_WINDOWS_HALF"
    --softcore-band "$SOFTCORE_BAND"
    --n-cycles "$N_CYCLES"
    --md-steps-per-cycle "$MD_STEPS_PER_CYCLE"
    --platform "$PLATFORM"
    --timestep-fs "$TIMESTEP_FS"
    --minimize-iters "$MINIMIZE_ITERS"
    --backward-equil-steps "$BACKWARD_EQUIL_STEPS"
    --genuine-decouple-nm "$GENUINE_DECOUPLE_NM"
    --binder-chain "$BINDER_CHAIN"
    --device-index "$DEVICE_INDEX"
    --mintimeid "$MINTIMEID"
    --out-root "$OUT_ROOT" )

# A single shell string for the detached nohup (cwd = PROJ on the VM).
RERUN_STR="cd ${PROJ} && nohup ${RERUN_CMD[*]} > ${OUT_ROOT_ABS%/*}/$(basename "$RERUN_LOG") 2>&1 < /dev/null &"

log "STEP 5: rerun plan (scenario ${SCENARIO})"
log "  VM memory target: ${VM_GIB} GiB (${VM_KIB} KiB)   pool concurrency: ${CONCURRENCY}"
log "  rerun (verbatim, detached on VM):"
printf '    %s\n' "${RERUN_CMD[*]}"

# ---------------------------------------------------------------------------
# DRY-RUN gate. Without --go, stop here. Nothing has been changed or launched.
# ---------------------------------------------------------------------------
if [ "$GO" -ne 1 ]; then
    log "DRY-RUN complete. No VM memory change, no backup write, no launch."
    log "Re-run with --go (after confirming scenario + integrity above) to execute."
    exit 0
fi

# =============================== --go path ==================================
log "--go: executing recovery (scenario ${SCENARIO})"

[ "$VM_STATE" = "running" ] || die "VM is not running; start it first (see runbook STEP 0) — refusing to act blind."
[ "${CP4_DPLUS_BAD:-1}" -eq 0 ] || die "cp4 dplus integrity not confirmed; refusing to rerun."

# --- STEP 4 (execute): backup cp4 dplus and verify count before touching VM mem.
log "STEP 4 (execute): backing up cp4 dplus .out files on the VM"
BK_RESULT="$(ssh_vm "bash -s" <<<"$BACKUP_CMD")"
echo "$BK_RESULT" | sed 's/^/  /'
BK_N="$(echo "$BK_RESULT" | sed -n 's/^backup_out_files=//p')"
EXPECT_BK=$(( N_REPS * EXPECTED_WINDOWS ))   # 3 reps x 6 windows = 18
if [ "${BK_N:-0}" -lt "$EXPECT_BK" ]; then
    die "dplus backup incomplete (${BK_N:-0}/${EXPECT_BK} .out files). STOP before any rerun."
fi
log "  dplus backup verified: ${BK_N}/${EXPECT_BK} .out files copied (R-7)."

# --- STEP 2 (execute): VM memory reallocation. setmaxmem requires shutoff.
CUR_MAX_KIB="$(virsh dominfo "$VM_DOMAIN" 2>/dev/null | awk -F: '/Max memory/{gsub(/[^0-9]/,"",$2);print $2}')"
log "STEP 2 (execute): VM memory ${CUR_MAX_KIB:-?} KiB -> ${VM_KIB} KiB"
if [ "${CUR_MAX_KIB:-0}" != "$VM_KIB" ]; then
    log "  shutting down VM to set --config maxmem (no hotplug slot present)..."
    virsh shutdown "$VM_DOMAIN" || warn "virsh shutdown returned non-zero"
    # Wait (bounded) for shutoff.
    for _ in $(seq 1 60); do
        st="$(virsh domstate "$VM_DOMAIN" 2>/dev/null || echo unknown)"
        [ "$st" = "shut off" ] && break
        sleep 5
    done
    st="$(virsh domstate "$VM_DOMAIN" 2>/dev/null || echo unknown)"
    [ "$st" = "shut off" ] || die "VM did not reach 'shut off' (state=${st}). Resolve manually (runbook STEP 2)."
    virsh setmaxmem "$VM_DOMAIN" "$VM_KIB" --config
    virsh setmem    "$VM_DOMAIN" "$VM_KIB" --config
    virsh start     "$VM_DOMAIN"
    log "  VM restarted; waiting for SSH..."
    for _ in $(seq 1 60); do
        if ssh_vm "true" 2>/dev/null; then break; fi
        sleep 5
    done
    ssh_vm "true" 2>/dev/null || die "VM SSH not reachable after restart. Resolve manually."
    log "  VM SSH reachable. New Max memory: $(virsh dominfo "$VM_DOMAIN" | awk -F: '/Max memory/{print $2}')"
else
    log "  VM already at target memory; no shutdown needed."
fi

# --- STEP 5 (execute): launch cp4 dminus rerun detached on the VM.
log "STEP 5 (execute): launching cp4 ${RECOVER_DIRECTION} rerun (detached) on VM"
ssh_vm "$RERUN_STR"
sleep 3
log "  launched. Driver/worker processes on VM:"
ssh_vm "pgrep -af 'trackb_inplace_rbfe_production' || echo '    (none yet — check the log)'" | sed 's/^/    /'
log "  rerun log (on VM): ${OUT_ROOT}/$(basename "$RERUN_LOG")"

# --- post-launch dplus non-destruction verify (sentinel).
log "POST-LAUNCH: verifying cp4 dplus still 400x6x3 (non-destruction sentinel)"
POST="$(ssh_vm "bash -s" <<<"$REMOTE_PROBE" || true)"
echo "$POST" | grep '^cp4 bound rep[0-9] dplus' | sed 's/^/  /'
POST_BAD="$(echo "$POST" | grep -c '^cp4 bound rep[0-9] dplus INCOMPLETE' || true)"
if [ "${POST_BAD:-0}" -ne 0 ]; then
    die "cp4 dplus became INCOMPLETE after launch (${POST_BAD} reps)! Restore from ${DPLUS_BACKUP_DIR} immediately (runbook STEP 4 restore)."
fi
log "  cp4 dplus intact after launch."

log "DONE. Monitor the rerun, then follow runbook STEP 6 (merge -> rsync -> UWHAM -> ddint_bound -> ddG_bind sign)."
