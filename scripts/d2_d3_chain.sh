#!/bin/bash
# ============================================================================
# D2 + D3 chain for v0.9 Layer 3D pilot.
#
# Sequential pipeline:
#   D2 — Ace-OMW-NMe B3LYP/6-31G(d) geometry optimization (~1-4 h CPU)
#   D3 — Constrained dihedral/angle/improper scans on D2-opt geom (~5-15 h CPU)
#   D2+D3 morning preview symlink update
#
# Run:
#   nohup bash scripts/d2_d3_chain.sh > /tmp/d2_d3_nohup.log 2>&1 &
# ============================================================================

set -u

PROJ=/home/san/UPDD_proj
cd "$PROJ"
source /home/san/miniconda3/etc/profile.d/conda.sh

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGFILE="/tmp/d2_d3_chain_${TS_TAG}.log"

log () {
    printf '[d2_d3_chain %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$LOGFILE"
}

log "============================================================"
log "D2 + D3 chain START"
log "============================================================"

conda activate psi4

# Phase D2 — B3LYP/6-31G(d) geometry optimization
log "Phase D2: B3LYP/6-31G(d) geometry optimization (~1-4 h)"
START_D2=$(date +%s)
python "$PROJ/scripts/d2_omw_b3lyp_geom_opt.py" 2>&1 | tee -a "$LOGFILE"
D2_RC=${PIPESTATUS[0]}
ELAPSED_D2=$(($(date +%s) - START_D2))
log "D2 exit code: $D2_RC (wall: ${ELAPSED_D2}s = $((ELAPSED_D2 / 60)) min)"

if [ "$D2_RC" -ne 0 ]; then
    log "D2 FAILED — skipping D3"
    exit 1
fi

# Phase D3 — QM scans on D2-optimized geometry
log "Phase D3: QM dihedral/angle/improper scans (~5-15 h, 3 scans × multi-point)"
START_D3=$(date +%s)
python "$PROJ/scripts/d3_omw_qm_scans.py" 2>&1 | tee -a "$LOGFILE"
D3_RC=${PIPESTATUS[0]}
ELAPSED_D3=$(($(date +%s) - START_D3))
log "D3 exit code: $D3_RC (wall: ${ELAPSED_D3}s = $((ELAPSED_D3 / 60)) min)"

conda deactivate

# Phase morning preview
log "============================================================"
log "Updating MORNING_REVIEW_D2_D3.md"
log "============================================================"

LATEST_D2=$(ls -td "$PROJ"/outputs/analysis/d2_omw_geom_opt_*/ 2>/dev/null | head -1)
LATEST_D3=$(ls -td "$PROJ"/outputs/analysis/d3_omw_qm_scans_*/ 2>/dev/null | head -1)

PREVIEW="$PROJ/MORNING_REVIEW_D2_D3.md"
{
    echo "# D2 + D3 morning review — auto-generated $(date '+%Y-%m-%d %H:%M:%S')"
    echo
    echo "## Summary"
    echo
    echo "| Phase | Wall time | Exit code | Result dir |"
    echo "|---|---|---|---|"
    echo "| D2 B3LYP geom opt | $((ELAPSED_D2 / 60)) min | $D2_RC | ${LATEST_D2:-NOT FOUND} |"
    echo "| D3 QM scans       | $((ELAPSED_D3 / 60)) min | $D3_RC | ${LATEST_D3:-NOT FOUND} |"
    echo
    if [ -n "$LATEST_D2" ] && [ -f "$LATEST_D2/d2_geom_opt_report.md" ]; then
        echo "## D2 Geometry optimization report"
        echo
        cat "$LATEST_D2/d2_geom_opt_report.md"
        echo
    fi
    if [ -n "$LATEST_D3" ] && [ -f "$LATEST_D3/d3_qm_scans_report.md" ]; then
        echo "## D3 QM Scans report"
        echo
        cat "$LATEST_D3/d3_qm_scans_report.md"
        echo
    fi
    echo
    echo "## Chain log"
    echo
    echo "Full log: $LOGFILE"
} > "$PREVIEW"

log "MORNING_REVIEW: $PREVIEW ($(wc -l "$PREVIEW" | awk '{print $1}') lines)"
log "D2+D3 chain COMPLETE. Total wall: $((($(date +%s) - START_D2) / 60)) min"
