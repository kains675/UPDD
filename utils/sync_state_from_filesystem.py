#!/usr/bin/env python
"""
utils/sync_state_from_filesystem.py
-----------------------------------
UPDD 파이프라인이 ``--no-dashboard`` 로 구동되어 ``state.json`` 이 쓰이지
않는 동안에도 ``--dashboard-only`` 옵저버가 현재 진행도를 볼 수 있도록,
실제 산출물 파일을 스캔하여 ``state.json`` 을 한 번 생성/갱신한다.

사용:
    python utils/sync_state_from_filesystem.py \
        --project outputs/6WGN_cyclic_htc_NMA_12_12-16 \
        --log-dir log/6WGN

자동 감지:
    python utils/sync_state_from_filesystem.py --target 6WGN

watch 와 결합하여 주기 갱신:
    watch -n 10 'python utils/sync_state_from_filesystem.py --target 6WGN'

FIX 이후 (PipelineState 가 항상 쓰기) 에는 rerun 하면 자동 갱신되므로 본 유틸
은 "이미 no-dashboard 로 시작된 세션" 용 hotfix.
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import sys
import time
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
try:
    from updd_state import PipelineState, PIPELINE_STEPS  # type: ignore
    from updd_state import check_outputs_exist  # type: ignore
except ImportError as e:
    sys.exit(f"[!] updd_state import 실패: {e}")

UPDD_ROOT = os.path.abspath(os.path.join(HERE, os.pardir))


def load_inputs(project_dir: str) -> dict:
    p = os.path.join(project_dir, "updd_status.json")
    if not os.path.exists(p):
        return {}
    try:
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}


def count_md_progress(project_dir: str) -> tuple[int, int, int]:
    """Return (n_final, n_exploded, n_input)."""
    md = os.path.join(project_dir, "mdresult")
    ncaa = os.path.join(project_dir, "ncaa_mutants")
    n_final = len(glob.glob(os.path.join(md, "*_final.pdb")))
    n_exp = len(glob.glob(os.path.join(md, "*_EXPLODED_*.log")))
    # input PDB count (exclude _renum / _minimized / _partial)
    n_input = 0
    for p in glob.glob(os.path.join(ncaa, "**", "*.pdb"), recursive=True):
        if not any(x in p for x in ("_renum.pdb", "_minimized.pdb",
                                     "_partial_")):
            n_input += 1
    return n_final, n_exp, n_input


def find_current_md_design(project_dir: str) -> str:
    """Return the basename (stem) of the currently-active MD design, or ''."""
    md = os.path.join(project_dir, "mdresult")
    logs = sorted(
        glob.glob(os.path.join(md, "*_md.log")),
        key=lambda p: os.path.getmtime(p), reverse=True,
    )
    if not logs:
        return ""
    # Fresh if modified within 2 minutes
    if time.time() - os.path.getmtime(logs[0]) < 120:
        return os.path.basename(logs[0]).replace("_md.log", "")
    return ""


STEP_OUTPUT_PROBES = {
    # (step_name, check_key_for_check_outputs_exist_or_empty)
    "preprocess": "preprocess",
    "rfdiff": "rfdiffusion",
    "proteinmpnn": "proteinmpnn",
    "af2": "af2",
    "ncaa_mutation": "ncaa_mutation",
    "parameterize": "parameterize",
    "md": "md",
    "snapshots": "snapshots",
    "qmmm": "qmmm",
    "mmgbsa": "mmgbsa",
    "rank": "",  # no cheap probe, only mark complete from completed_steps
}


def sync(project_dir: str, log_dir: str) -> None:
    inputs = load_inputs(project_dir)
    project_name = os.path.basename(project_dir.rstrip("/"))

    state = PipelineState(log_dir=log_dir, project_name=project_name)
    completed = inputs.get("completed_steps", []) if isinstance(inputs, dict) else []
    completed_set = {str(c).lower() for c in completed}

    # Legacy → new name map (mirrors _STATE_STEP_MAP in UPDD.py)
    legacy_map = {
        "preprocess": ["preprocess"],
        "rfdiff": ["rfdiff"],
        "rfdiffusion": ["rfdiff"],
        "proteinmpnn": ["proteinmpnn"],
        "mpnn": ["proteinmpnn"],
        "af2": ["af2"],
        "alphafold": ["af2"],
        "ncaa_mutation": ["ncaa_mutation"],
        "parameterize": ["parameterize"],
        "md": ["md"],
        "md_snapshots": ["md", "snapshots"],
        "snapshots": ["snapshots"],
        "qmmm": ["qmmm"],
        "mmgbsa": ["mmgbsa"],
        "rank": ["rank"],
    }

    # Mark completed steps (skip=True) — first from updd_status.json hints,
    # then re-verify via filesystem probes for all known steps.
    marked = set()
    for leg, targets in legacy_map.items():
        if leg in completed_set:
            for name in targets:
                probe = STEP_OUTPUT_PROBES.get(name, "")
                if probe and check_outputs_exist(probe, inputs, project_dir):
                    state.complete_step(name, skipped=True)
                    marked.add(name)
                elif not probe:
                    state.complete_step(name, skipped=True)
                    marked.add(name)

    # Filesystem-based backfill for steps that have output artifacts but
    # weren't in completed_steps (e.g. parameterize → params/, af2 →
    # af2_results/). preprocess always counts as done (input step).
    if "preprocess" not in marked:
        state.complete_step("preprocess", skipped=True)
    for name in ("rfdiff", "proteinmpnn", "af2", "ncaa_mutation",
                 "parameterize"):
        if name in marked:
            continue
        probe = STEP_OUTPUT_PROBES.get(name, "")
        if probe and check_outputs_exist(probe, inputs, project_dir):
            state.complete_step(name, skipped=True)
            marked.add(name)

    # MD: derive running progress from filesystem
    n_final, n_exp, n_input = count_md_progress(project_dir)
    if n_input > 0 and (n_final + n_exp) > 0:
        current = find_current_md_design(project_dir)
        if n_final + n_exp < n_input and state.steps["md"].status.value != "complete":
            # currently running
            state.start_step("md", total=n_input, workers=1,
                             detail=f"{n_final} SUCCESS / {n_exp} EXPLODED "
                                    f"/ {n_input} total"
                                    + (f" — running: {current[:30]}..." if current else ""))
            state.update_step("md", completed=n_final)
        elif n_final == n_input:
            state.complete_step("md")

    # Additional hints for snapshots
    snap_dir = os.path.join(project_dir, "snapshots")
    if os.path.isdir(snap_dir):
        snap_pdbs = glob.glob(os.path.join(snap_dir, "**", "*.pdb"),
                              recursive=True)
        if snap_pdbs:
            state.complete_step("snapshots", skipped=False)

    # Emit a sync marker event
    state.emit_event(
        "info",
        f"🔄 state synced from filesystem "
        f"({n_final}✅/{n_exp}❌/{n_input} designs, "
        f"project={project_name})",
    )
    print(f"✅ state.json synced → {state.state_path}")
    print(f"   project  : {project_name}")
    print(f"   progress : {state.total_progress*100:.1f}%")
    print(f"   md       : {n_final} success / {n_exp} exploded / {n_input} total")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--project",
                    help="outputs/<project_folder> 경로 (생략 시 --target 으로 자동 탐색)")
    ap.add_argument("--log-dir",
                    help="log/<target>/ 경로 (생략 시 --target 으로 자동 탐색)")
    ap.add_argument("--target",
                    help="target_id (예: 6WGN) — project / log-dir 자동 탐색")
    args = ap.parse_args()

    project_dir = args.project
    log_dir = args.log_dir

    if args.target:
        if not project_dir:
            cands = [d for d in glob.glob(
                os.path.join(UPDD_ROOT, "outputs", f"{args.target}_*"))
                if os.path.isdir(d)]
            if not cands:
                sys.exit(f"[!] '{args.target}*' project folder 없음")
            project_dir = max(cands, key=os.path.getmtime)
            print(f"[auto] project = {project_dir}")
        if not log_dir:
            log_dir = os.path.join(UPDD_ROOT, "log", args.target)
            print(f"[auto] log_dir = {log_dir}")

    if not project_dir or not log_dir:
        ap.error("--project + --log-dir 지정 또는 --target 필요")

    if not os.path.isdir(project_dir):
        sys.exit(f"[!] project 디렉토리 없음: {project_dir}")
    os.makedirs(log_dir, exist_ok=True)

    sync(project_dir, log_dir)


if __name__ == "__main__":
    main()
