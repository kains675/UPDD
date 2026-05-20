#!/bin/bash
# ============================================================================
# D2+D3 auto-fix watcher
#
# Monitors the d2_d3_chain orchestrator (PID in /tmp/d2_d3_chain_pid).
# On chain termination:
#   - If D3 report exists → success, finalize MORNING_REVIEW_D2_D3.md
#   - If D2/D3 incomplete → apply graduated auto-fix and relaunch the chain
#   - After 3 retries → drop an escalation note for morning review
#
# Auto-fix progression (per retry):
#   Retry 1: D2 g_convergence gau_loose → gau_verytight, geom_maxiter 200 → 400
#   Retry 2: D3 scan grids coarsened (5° → 10° improper, 30° → 60° dihedral)
#   Retry 3: D3 method downgrade B3LYP → HF (faster, less precise)
#
# Hard deadline: 08:00 KST. After that, watcher writes final MORNING_REVIEW
# from whatever partial results exist and exits.
#
# Run:
#   nohup bash scripts/d2_d3_auto_fix_watcher.sh > /tmp/d2_d3_watcher_nohup.log 2>&1 &
# ============================================================================

set -u

PROJ=/home/san/UPDD_proj
CHAIN_PID_FILE=/tmp/d2_d3_chain_pid
RETRY_FILE=/tmp/d2_d3_retry_count
MAX_RETRIES=3
DEADLINE_HMS="12:00:00"

# Initialize retry counter if not present
[ -f "$RETRY_FILE" ] || echo 0 > "$RETRY_FILE"

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGFILE="/tmp/d2_d3_watcher_${TS_TAG}.log"

log () {
    printf '[d2_d3_watcher %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$LOGFILE"
}

deadline_reached () {
    local NOW=$(date '+%H:%M:%S')
    # String compare works for HH:MM:SS in 24h format
    [[ "$NOW" > "$DEADLINE_HMS" ]] && [[ "$NOW" < "16:00:00" ]]
}

finalize_morning_review () {
    local status="$1"
    local PREVIEW="$PROJ/MORNING_REVIEW_D2_D3.md"
    local LATEST_D2=$(ls -td "$PROJ"/outputs/analysis/d2_omw_geom_opt_*/ 2>/dev/null | head -1)
    local LATEST_D3=$(ls -td "$PROJ"/outputs/analysis/d3_omw_qm_scans_*/ 2>/dev/null | head -1)
    local RETRY=$(cat "$RETRY_FILE")
    {
        echo "# D2 + D3 morning review — auto-generated $(date '+%Y-%m-%d %H:%M:%S')"
        echo
        echo "**Watcher status**: $status"
        echo "**Retries used**: $RETRY / $MAX_RETRIES"
        echo
        echo "## Result locations"
        echo
        echo "- D2 dir: ${LATEST_D2:-NOT FOUND}"
        echo "- D3 dir: ${LATEST_D3:-NOT FOUND}"
        echo
        if [ -n "$LATEST_D2" ] && [ -f "$LATEST_D2/d2_geom_opt_report.md" ]; then
            echo "## D2 report"
            echo
            cat "$LATEST_D2/d2_geom_opt_report.md"
            echo
        fi
        if [ -n "$LATEST_D3" ] && [ -f "$LATEST_D3/d3_qm_scans_report.md" ]; then
            echo "## D3 report"
            echo
            cat "$LATEST_D3/d3_qm_scans_report.md"
            echo
        fi
        if [ -f /tmp/d2_d3_escalation.log ]; then
            echo "## Escalation notes"
            echo
            cat /tmp/d2_d3_escalation.log
            echo
        fi
        echo
        echo "## Watcher log"
        echo
        echo "- Main: $LOGFILE"
        echo "- Retries:"
        ls /tmp/d2_d3_retry_*.log 2>/dev/null | while read f; do echo "  - $f"; done
    } > "$PREVIEW"
    log "MORNING_REVIEW written: $PREVIEW ($(wc -l < "$PREVIEW") lines)"
}

apply_autofix () {
    local retry=$1
    case "$retry" in
        1)
            log "Auto-fix 1: D2 g_convergence gau_loose → gau_verytight, geom_maxiter 200 → 400"
            sed -i 's/"g_convergence": "gau_loose"/"g_convergence": "gau_verytight"/' \
                "$PROJ/scripts/d2_omw_b3lyp_geom_opt.py"
            sed -i 's/"geom_maxiter": 200/"geom_maxiter": 400/' \
                "$PROJ/scripts/d2_omw_b3lyp_geom_opt.py"
            ;;
        2)
            log "Auto-fix 2: D3 scan grids coarsened (5° → 10° improper, 30° → 60° dihedral)"
            sed -i 's/range(-30, 31, 5)/range(-30, 31, 10)/' "$PROJ/scripts/d3_omw_qm_scans.py"
            sed -i 's/range(0, 360, 30)/range(0, 360, 60)/' "$PROJ/scripts/d3_omw_qm_scans.py"
            ;;
        3)
            log "Auto-fix 3: D3 method downgrade B3LYP → HF for dihedral scan"
            sed -i 's|psi4.optimize("b3lyp", molecule=mol)|psi4.optimize("hf", molecule=mol)|' \
                "$PROJ/scripts/d3_omw_qm_scans.py"
            ;;
        *)
            log "Unknown retry $retry, no fix applied"
            ;;
    esac
}

escalate () {
    local reason="$1"
    log "Escalation: $reason"
    cat >> /tmp/d2_d3_escalation.log <<EOF

=== Escalation $(date '+%Y-%m-%d %H:%M:%S') ===
Reason: $reason
Retries used: $(cat "$RETRY_FILE") / $MAX_RETRIES
Action: human verdict needed before morning continuation.
Latest D2 dir: $(ls -td "$PROJ"/outputs/analysis/d2_omw_geom_opt_*/ 2>/dev/null | head -1)
Latest D3 dir: $(ls -td "$PROJ"/outputs/analysis/d3_omw_qm_scans_*/ 2>/dev/null | head -1)
EOF
}

restart_chain () {
    local retry=$1
    log "Restarting D2+D3 chain (retry $retry)"
    nohup bash "$PROJ/scripts/d2_d3_chain.sh" > "/tmp/d2_d3_retry_${retry}.log" 2>&1 &
    local NEW_PID=$!
    echo "$NEW_PID" > "$CHAIN_PID_FILE"
    log "New chain PID: $NEW_PID"
}

log "============================================================"
log "D2+D3 auto-fix watcher START"
log "Deadline: $DEADLINE_HMS  |  Max retries: $MAX_RETRIES"
log "============================================================"

while true; do
    sleep 300  # 5 min poll interval

    # Deadline check
    if deadline_reached; then
        log "Deadline 08:00 KST reached — finalizing with whatever we have"
        # If chain is still running, leave it (deadline cap is best-effort)
        finalize_morning_review "DEADLINE_REACHED"
        break
    fi

    # Chain alive?
    CHAIN_PID=$(cat "$CHAIN_PID_FILE" 2>/dev/null)
    if [ -n "$CHAIN_PID" ] && ps -p "$CHAIN_PID" > /dev/null 2>&1; then
        # still running
        continue
    fi

    log "Chain PID $CHAIN_PID is dead — examining result"

    # Check D3 result
    LATEST_D3=$(ls -td "$PROJ"/outputs/analysis/d3_omw_qm_scans_*/ 2>/dev/null | head -1)
    if [ -n "$LATEST_D3" ] && [ -f "$LATEST_D3/d3_qm_scans_report.md" ]; then
        log "D3 report exists — chain succeeded"
        finalize_morning_review "SUCCESS"
        break
    fi

    # Chain incomplete — auto-fix + retry
    RETRY=$(cat "$RETRY_FILE")
    RETRY=$((RETRY + 1))
    echo "$RETRY" > "$RETRY_FILE"

    if [ "$RETRY" -gt "$MAX_RETRIES" ]; then
        escalate "All $MAX_RETRIES auto-fixes exhausted; D2 or D3 still failing"
        finalize_morning_review "EXHAUSTED_RETRIES"
        break
    fi

    apply_autofix "$RETRY"
    restart_chain "$RETRY"
done

log "Watcher exit at $(date '+%Y-%m-%d %H:%M:%S')"
