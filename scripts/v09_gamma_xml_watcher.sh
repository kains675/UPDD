#!/bin/bash
# ============================================================================
# γ XML completion watcher + auto-launch smoke test
#
# Polls every 60 s for `params/MTR_gaff2_layer3d_gamma.xml` (XML build output).
# Once detected + valid (Σq=0 cross-check), launches the γ smoke test via
# `scripts/v09_gamma_smoke.sh` in nohup mode.
#
# Hard timeout: 4 h (= 12:30 KST if started 08:30); if XML never appears,
# logs an escalation note and exits.
#
# Run:
#   nohup bash scripts/v09_gamma_xml_watcher.sh > /tmp/v09_gamma_watcher_nohup.log 2>&1 &
# ============================================================================

set -u

PROJ=/home/san/UPDD_proj
GAMMA_XML="${PROJ}/params/MTR_gaff2_layer3d_gamma.xml"

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGFILE="/tmp/v09_gamma_watcher_${TS_TAG}.log"
MAX_WAIT_S=$((4 * 3600))   # 4 h

log () {
    printf '[v09_gamma_watcher %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$LOGFILE"
}

log "============================================================"
log "γ XML completion watcher START"
log "Watching: $GAMMA_XML"
log "Timeout: $((MAX_WAIT_S / 3600)) h"
log "============================================================"

START_TS=$(date +%s)

while true; do
    NOW=$(date +%s)
    ELAPSED=$((NOW - START_TS))

    if [ "$ELAPSED" -ge "$MAX_WAIT_S" ]; then
        log "TIMEOUT after $((ELAPSED / 60)) min — γ XML never appeared. Exiting."
        exit 1
    fi

    if [ -f "$GAMMA_XML" ]; then
        log "γ XML detected after $((ELAPSED / 60)) min"
        SIGMA_Q=$(/home/san/miniconda3/envs/qmmm/bin/python -c "
import xml.etree.ElementTree as ET
try:
    t = ET.parse('$GAMMA_XML')
    res = next((r for r in t.findall('.//Residue') if r.get('name') == 'MTR'), None)
    if res is None:
        print('NO_MTR_RESIDUE')
    else:
        total = sum(float(a.get('charge', 0)) for a in res.findall('Atom'))
        print(f'{total:+.6f}')
except Exception as e:
    print(f'PARSE_ERROR:{e}')
")
        log "γ XML Σq = $SIGMA_Q"
        if [[ "$SIGMA_Q" =~ ^[+-]?[0-9]+\.[0-9]+$ ]]; then
            ABS=$(awk -v s="$SIGMA_Q" 'BEGIN { print (s < 0 ? -s : s) }')
            if awk -v a="$ABS" 'BEGIN { exit !(a < 0.00001) }'; then
                log "Σq verification PASS — launching γ smoke test"
                break
            else
                log "Σq verification FAIL ($ABS >= 1e-5) — sleep 60s, recheck"
            fi
        else
            log "Σq parse anomaly ($SIGMA_Q) — sleep 60s, recheck"
        fi
    fi

    sleep 60
done

# Launch γ smoke test
log "============================================================"
log "Launching γ smoke test (nohup)..."
log "============================================================"
nohup bash "$PROJ/scripts/v09_gamma_smoke.sh" > /tmp/v09_gamma_smoke_nohup.log 2>&1 &
GAMMA_SMOKE_PID=$!
echo "$GAMMA_SMOKE_PID" > /tmp/v09_gamma_smoke_pid
log "γ smoke PID: $GAMMA_SMOKE_PID"
log "γ smoke log: /tmp/v09_gamma_smoke_nohup.log"
log "Smoke ETA: ~3 h (4 seeds × 40 min sequential)"
log "============================================================"
log "Watcher exit (smoke handoff complete)"
