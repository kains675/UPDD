#!/bin/bash
# ============================================================================
# D1 auto-chain watcher: PySCF Coder completion → Psi4 D1 → 3-stack cross-validation
#
# Polls every 60 s for the PySCF reproduction's `omw_pyscf_charges.json` (the
# Coder-produced output). Once detected, sequentially:
#   1. Runs the Psi4 stack reproduction (`scripts/d1_khoury_psi4_resp.py`)
#   2. Runs the 3-stack cross-validation (`scripts/d1_qmstack_crossvalidation.py`)
#   3. Creates `MORNING_REVIEW_D1_RESULT.md` symlink for quick access
#
# Timeout: 8 h. If PySCF output never arrives, proceeds to Psi4 standalone
# (RDKit-built geometry fallback in the psi4 script).
#
# Run:
#   nohup bash scripts/d1_auto_chain.sh > /tmp/d1_chain_nohup.log 2>&1 &
# ============================================================================

set -u

PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGFILE="/tmp/d1_chain_${TS_TAG}.log"

log () {
    printf '[d1_chain %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$LOGFILE"
}

log "==============================================================="
log "D1 auto-chain watcher start"
log "Logfile: $LOGFILE"
log "Polling cadence: 60 s | Timeout: 8 h | Stack chain: PySCF → Psi4 → 3-way"
log "==============================================================="

# Phase 1 — Wait for PySCF Coder output (omw_pyscf_charges.json)
MAX_WAIT_S=$((8 * 3600))
START_TS=$(date +%s)
PYSCF_CHARGES=""
PYSCF_REPORT=""

while true; do
    NOW=$(date +%s)
    ELAPSED=$((NOW - START_TS))

    if [ "$ELAPSED" -ge "$MAX_WAIT_S" ]; then
        log "TIMEOUT after $((MAX_WAIT_S / 3600)) h — proceeding to Psi4 standalone (no PySCF reuse)"
        break
    fi

    # Look for latest PySCF output dir with charges json + report
    PYSCF_DIR=$(ls -td "$PROJ"/outputs/analysis/d1_khoury_resp_reproduction_* 2>/dev/null | head -1)
    if [ -n "$PYSCF_DIR" ] && [ -d "$PYSCF_DIR" ]; then
        # Try standard json name first, fallback to dat / npy
        for cand in omw_pyscf_charges.json pyscf_resp_charges.json omw_charges.json final_charges.json; do
            if [ -f "$PYSCF_DIR/$cand" ]; then
                PYSCF_CHARGES="$PYSCF_DIR/$cand"
                break
            fi
        done
        # Report file (any *.md with "report" in name)
        PYSCF_REPORT=$(ls "$PYSCF_DIR"/*report*.md 2>/dev/null | head -1)
    fi

    if [ -n "$PYSCF_CHARGES" ]; then
        log "PySCF reproduction DETECTED after ${ELAPSED}s (~$((ELAPSED / 60)) min)"
        log "  charges: $PYSCF_CHARGES"
        [ -n "$PYSCF_REPORT" ] && log "  report:  $PYSCF_REPORT"
        break
    fi

    # Periodic status log every 10 minutes
    if [ $((ELAPSED % 600)) -lt 60 ] && [ "$ELAPSED" -gt 60 ]; then
        log "  [waiting] ${ELAPSED}s elapsed, no PySCF output yet"
    fi

    sleep 60
done

# Phase 2 — Psi4 D1 reproduction
log "==============================================================="
log "Phase 2: Psi4 D1 RESP-A2 reproduction"
log "==============================================================="
conda activate psi4
python "$PROJ/scripts/d1_khoury_psi4_resp.py" 2>&1 | tee -a "$LOGFILE"
PSI4_RC=${PIPESTATUS[0]}
conda deactivate
log "Psi4 D1 exit code: $PSI4_RC"

# Phase 3 — 3-stack cross-validation
log "==============================================================="
log "Phase 3: 3-stack cross-validation (PySCF vs Psi4 vs Khoury)"
log "==============================================================="
conda activate qmmm
python "$PROJ/scripts/d1_qmstack_crossvalidation.py" 2>&1 | tee -a "$LOGFILE"
CROSSVAL_RC=${PIPESTATUS[0]}
conda deactivate
log "Cross-validation exit code: $CROSSVAL_RC"

# Phase 4 — Morning review symlink
log "==============================================================="
log "Phase 4: Morning review preview"
log "==============================================================="

LATEST_CROSSVAL=$(ls -td "$PROJ"/outputs/analysis/d1_qmstack_crossvalidation_*/ 2>/dev/null | head -1)
LATEST_PSI4=$(ls -td "$PROJ"/outputs/analysis/d1_psi4_resp_reproduction_*/ 2>/dev/null | head -1)
LATEST_PYSCF=$(ls -td "$PROJ"/outputs/analysis/d1_khoury_resp_reproduction_*/ 2>/dev/null | head -1)

PREVIEW="$PROJ/MORNING_REVIEW_D1_RESULT.md"
{
    echo "# D1 morning review — auto-generated $(date '+%Y-%m-%d %H:%M:%S')"
    echo
    echo "## Sources"
    echo
    echo "- PySCF reproduction dir: ${LATEST_PYSCF:-NOT FOUND}"
    echo "- Psi4 reproduction dir:  ${LATEST_PSI4:-NOT FOUND}"
    echo "- Cross-validation dir:   ${LATEST_CROSSVAL:-NOT FOUND}"
    echo
    echo "## Exit codes"
    echo
    echo "- Psi4 D1: $PSI4_RC ($([ $PSI4_RC -eq 0 ] && echo OK || echo NOT-PASS))"
    echo "- Cross-validation: $CROSSVAL_RC ($([ $CROSSVAL_RC -eq 0 ] && echo OK || echo NOT-PASS))"
    echo
    if [ -n "$LATEST_CROSSVAL" ] && [ -f "$LATEST_CROSSVAL/comparison.md" ]; then
        echo "## 3-stack cross-validation report"
        echo
        cat "$LATEST_CROSSVAL/comparison.md"
    elif [ -n "$LATEST_PSI4" ] && [ -f "$LATEST_PSI4/d1_psi4_reproduction_report.md" ]; then
        echo "## Psi4 reproduction report (cross-validation unavailable)"
        echo
        cat "$LATEST_PSI4/d1_psi4_reproduction_report.md"
    elif [ -n "$LATEST_PYSCF" ] && [ -n "$PYSCF_REPORT" ]; then
        echo "## PySCF reproduction report (Psi4 + cross-validation unavailable)"
        echo
        cat "$PYSCF_REPORT"
    else
        echo "## No reports available — investigate logs:"
        echo "- Watcher log: $LOGFILE"
    fi
    echo
    echo "## Watcher log"
    echo
    echo "Full chain log at: $LOGFILE"
} > "$PREVIEW"

log "MORNING REVIEW preview: $PREVIEW"
log "  $(wc -l "$PREVIEW" | awk '{print $1}') lines"
log "D1 auto-chain COMPLETE. Total wall time: $(($(date +%s) - START_TS)) s"
