#!/usr/bin/env bash
set -euo pipefail

PROJ=/home/san/UPDD_proj
DCD_ROOT=${PROJ}/outputs/_trackb/w4a_postdensify_dcd_20260716
DCD_STATE=${DCD_ROOT}/run_state.json
DCD_SESSION=w4a_postdcd_20260716
GATE_SCRIPT=${PROJ}/analysis/w4a_postdensify_dcd_20260716/postdensify_gate.py
GATE_JSON=${PROJ}/analysis/w4a_postdensify_dcd_20260716/postdensify_gate.json
GATE_LOCK=${DCD_ROOT}/postdensify_gate.lock
RUNNER=${PROJ}/analysis/1ycr_wt_seed_expansion_20260716/expand_1ycr_wt_scaffolds.py
PYTHON=/home/san/miniconda3/envs/qmmm/bin/python

while true; do
  state=$(jq -r '.status // "UNKNOWN"' "${DCD_STATE}" 2>/dev/null || echo UNKNOWN)
  if [ "${state}" = "FAILED" ]; then
    echo "[1YCR-WAITER] DCD runner failed; scaffold queue will not launch."
    exit 1
  fi
  if [ "${state}" = "COMPLETE" ]; then
    break
  fi
  sleep 300
done

while tmux has-session -t "${DCD_SESSION}" 2>/dev/null; do
  sleep 10
done

exec 9>"${GATE_LOCK}"
flock 9
if [ ! -f "${GATE_JSON}" ]; then
  echo "[1YCR-WAITER] evaluating frozen post-densify DCD gate"
  "${PYTHON}" "${GATE_SCRIPT}"
fi
flock -u 9

echo "[1YCR-WAITER] DCD gate complete; launching independent 1YCR WT seed expansion"
cd "${PROJ}"
exec "${PYTHON}" "${RUNNER}"
