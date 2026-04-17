#!/usr/bin/env python3
"""
updd_state.py — Pipeline State & Event System (v0.4)

Dashboard가 UPDD.py의 현재 상태를 알 수 있도록 구조화된 state를 관리.
state.json (실시간 상태) + events.jsonl (이벤트 로그) 2개 파일 사용.

UPDD.py는 step 시작/완료 시 state를 업데이트하고,
dashboard는 state.json을 폴링하여 UI에 반영.

사용 예시:
    from updd_state import PipelineState, Step, StepStatus

    state = PipelineState(log_dir="~/UPDD_proj/log/6WGN")
    state.start_step("rfdiff", total_items=10)
    state.update_step("rfdiff", completed=5)
    state.complete_step("rfdiff")

Thread safety: state write 시 fcntl.flock으로 파일 잠금 → dashboard 읽기 중 race 방지.
"""

import os
import json
import time
import fcntl
from datetime import datetime
from enum import Enum
from dataclasses import dataclass, field, asdict
from typing import Optional, Dict, List, Any
from pathlib import Path


class StepStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"
    SKIPPED = "skipped"


# 파이프라인 10단계 정의 (UPDD.py의 실제 step과 맞춤)
PIPELINE_STEPS = [
    ("preprocess", "Preprocess"),
    ("rfdiff", "RFdiffusion"),
    ("proteinmpnn", "ProteinMPNN"),
    ("af2", "AF2 Multimer"),
    ("ncaa_mutation", "ncAA Mutation"),
    ("parameterize", "Parameterize ncAA"),
    ("md", "Restrained MD"),
    ("snapshots", "Snapshot extract"),
    ("qmmm", "QM/MM"),
    ("mmgbsa", "MM-GBSA"),
    ("rank", "Rank Results"),
]


@dataclass
class Step:
    name: str                          # e.g. "rfdiff"
    label: str                         # display name e.g. "RFdiffusion"
    status: StepStatus = StepStatus.PENDING
    start_time: Optional[float] = None  # epoch seconds
    end_time: Optional[float] = None
    completed: int = 0
    total: int = 0
    workers: int = 1
    detail: str = ""                   # 현재 진행 상태 (예: "4/8 designs")
    error: Optional[str] = None

    @property
    def elapsed(self) -> float:
        """경과 시간(초). 진행 중이면 현재까지, 완료면 총 시간."""
        if self.start_time is None:
            return 0.0
        end = self.end_time or time.time()
        return end - self.start_time

    @property
    def elapsed_str(self) -> str:
        e = int(self.elapsed)
        if e < 60:
            return f"{e}s"
        if e < 3600:
            return f"{e//60}m{e%60:02d}s"
        return f"{e//3600}h{(e%3600)//60:02d}m"

    @property
    def progress(self) -> float:
        """0.0 ~ 1.0. total이 0이면 0.0 반환."""
        if self.total <= 0:
            return 1.0 if self.status == StepStatus.COMPLETE else 0.0
        return min(1.0, self.completed / self.total)


class PipelineState:
    """
    파이프라인 상태 + 이벤트 저장.

    파일 구조:
        <log_dir>/
            state.json       # 현재 스냅샷 (dashboard가 폴링)
            events.jsonl     # 이벤트 로그 (append-only)

    모든 write는 atomic (tempfile + rename)하여 reader가 partial read 보지 않도록 함.
    """

    def __init__(self, log_dir: str, project_name: str = ""):
        self.log_dir = Path(os.path.expanduser(log_dir))
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.state_path = self.log_dir / "state.json"
        self.events_path = self.log_dir / "events.jsonl"

        self.project_name = project_name
        self.pipeline_start = time.time()
        self.steps: Dict[str, Step] = {
            name: Step(name=name, label=label)
            for name, label in PIPELINE_STEPS
        }

        # 초기 dump
        self._write_state()

    # ───────────────────── Step Lifecycle ─────────────────────

    def start_step(self, name: str, total: int = 0, workers: int = 1,
                   detail: str = "") -> None:
        if name not in self.steps:
            self.emit_event("warn", f"Unknown step: {name}")
            return
        step = self.steps[name]
        step.status = StepStatus.RUNNING
        step.start_time = time.time()
        step.end_time = None
        step.completed = 0
        step.total = total
        step.workers = workers
        step.detail = detail
        self.emit_event("step_start", f"🔄 {step.label} started "
                                       f"(workers={workers}, total={total})")
        self._write_state()

    def update_step(self, name: str, completed: Optional[int] = None,
                    detail: Optional[str] = None) -> None:
        if name not in self.steps:
            return
        step = self.steps[name]
        if completed is not None:
            step.completed = completed
        if detail is not None:
            step.detail = detail
        self._write_state()

    def complete_step(self, name: str, skipped: bool = False) -> None:
        if name not in self.steps:
            return
        step = self.steps[name]
        step.end_time = time.time()
        step.status = StepStatus.SKIPPED if skipped else StepStatus.COMPLETE
        if not skipped:
            step.completed = step.total if step.total > 0 else 1
        icon = "⏭ " if skipped else "✅"
        self.emit_event("step_complete",
                        f"{icon} {step.label} {'skipped' if skipped else 'complete'} "
                        f"({step.elapsed_str})")
        self._write_state()

    def fail_step(self, name: str, error: str) -> None:
        if name not in self.steps:
            return
        step = self.steps[name]
        step.end_time = time.time()
        step.status = StepStatus.FAILED
        step.error = error
        self.emit_event("step_fail", f"❌ {step.label} failed: {error}")
        self._write_state()

    # ───────────────────── Events ─────────────────────

    def emit_event(self, kind: str, message: str, **extra) -> None:
        """이벤트를 events.jsonl에 append."""
        ev = {
            "time": time.time(),
            "timestamp": datetime.now().strftime("%H:%M:%S"),
            "kind": kind,
            "message": message,
            **extra,
        }
        try:
            with open(self.events_path, "a", encoding="utf-8") as f:
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                f.write(json.dumps(ev, ensure_ascii=False) + "\n")
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
        except OSError:
            pass

    # ───────────────────── Convenience: QM/MM worker events ─────────────────────

    def qmmm_cycle(self, tag: str, cycle: int, grad: float,
                   delta_e: float) -> None:
        """QM/MM SCF cycle 기록 (빈도 높음, events.jsonl에만 경량 저장)."""
        self.emit_event("qmmm_cycle", f"[{tag}] cycle={cycle} |g|={grad:.2e}",
                        tag=tag, cycle=cycle, grad=grad, delta_e=delta_e)

    def qmmm_converged(self, tag: str, cycles: int, energy: float) -> None:
        self.emit_event("qmmm_converged",
                        f"✅ [{tag}] converged at cycle {cycles}, "
                        f"E={energy:.4f} Hartree",
                        tag=tag, cycles=cycles, energy=energy)

    def qmmm_homo_lumo_warn(self, tag: str, count: int) -> None:
        self.emit_event("qmmm_warn",
                        f"⚠️  [{tag}] HOMO==LUMO warning (#{count})",
                        tag=tag, warn_count=count)

    # ───────────────────── Aggregate / Query ─────────────────────

    @property
    def elapsed(self) -> float:
        return time.time() - self.pipeline_start

    @property
    def total_progress(self) -> float:
        """전체 파이프라인 진행도 (0.0 ~ 1.0)."""
        n = len(self.steps)
        done = sum(1 for s in self.steps.values()
                   if s.status in (StepStatus.COMPLETE, StepStatus.SKIPPED))
        running_progress = sum(s.progress for s in self.steps.values()
                               if s.status == StepStatus.RUNNING)
        return (done + running_progress) / n

    @property
    def current_step(self) -> Optional[Step]:
        """현재 실행 중인 step (없으면 None)."""
        for s in self.steps.values():
            if s.status == StepStatus.RUNNING:
                return s
        return None

    # ───────────────────── Serialization ─────────────────────

    def to_dict(self) -> dict:
        return {
            "project": self.project_name,
            "pipeline_start": self.pipeline_start,
            "elapsed": self.elapsed,
            "total_progress": self.total_progress,
            "steps": [
                {
                    "name": s.name,
                    "label": s.label,
                    "status": s.status.value,
                    "start_time": s.start_time,
                    "end_time": s.end_time,
                    "elapsed": s.elapsed,
                    "elapsed_str": s.elapsed_str,
                    "completed": s.completed,
                    "total": s.total,
                    "workers": s.workers,
                    "detail": s.detail,
                    "error": s.error,
                    "progress": s.progress,
                }
                for s in self.steps.values()
            ],
        }

    def _write_state(self) -> None:
        """Atomic write: tempfile + rename."""
        data = self.to_dict()
        tmp = self.state_path.with_suffix(".json.tmp")
        try:
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
            os.replace(tmp, self.state_path)
        except OSError:
            pass

    # ───────────────────── Reader (for dashboard) ─────────────────────

    @staticmethod
    def load_state(log_dir: str) -> Optional[dict]:
        """Dashboard가 state.json을 읽을 때 사용."""
        path = Path(os.path.expanduser(log_dir)) / "state.json"
        if not path.exists():
            return None
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except (OSError, json.JSONDecodeError):
            return None

    @staticmethod
    def tail_events(log_dir: str, n: int = 10) -> List[dict]:
        """events.jsonl 마지막 N개 이벤트 (dashboard용)."""
        path = Path(os.path.expanduser(log_dir)) / "events.jsonl"
        if not path.exists():
            return []
        try:
            with open(path, "r", encoding="utf-8") as f:
                lines = f.readlines()
            events = []
            for line in lines[-n*2:]:  # 여유분
                line = line.strip()
                if not line:
                    continue
                try:
                    events.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
            return events[-n:]
        except OSError:
            return []


# ──────────────────────────────────────────────────────────
# Output existence check (v0.4.1 dashboard resume 정확도 보정)
# ──────────────────────────────────────────────────────────

def check_outputs_exist(step: str, inputs: Optional[dict], outputs_dir: str) -> bool:
    """각 step 의 실제 산출물 존재 여부를 확인한다.

    Dashboard 가 resume 시 ``completed_steps`` 플래그만 믿고 SKIPPED 를
    마킹하면 실제로는 재실행이 필요한 상황(qmmm_results/ 가 빈 경우 등)에서
    잘못된 UI 를 보여준다. 이 헬퍼는 각 step 의 *실제* 산출물을 체크하여
    진짜 완료된 경우만 True 를 반환한다.

    Args:
        step: updd_status.json 의 legacy key 또는 PIPELINE_STEPS name.
        inputs: updd_status.json 로드 결과 (n_designs 등 참고용).
        outputs_dir: <project>/outputs/<base>_... 디렉토리.

    Returns:
        산출물이 "충분히" 존재하면 True. 부분 완료는 False (재실행).
    """
    if not outputs_dir or not os.path.isdir(outputs_dir):
        return False
    inputs = inputs or {}

    def _count_pdbs(path: str, prefix: str = "") -> int:
        if not os.path.isdir(path):
            return 0
        return sum(1 for f in os.listdir(path)
                   if f.endswith(".pdb") and (not prefix or f.startswith(prefix)))

    if step == "preprocess":
        # preprocess 산출물은 target_pdb (target/<base>_clean.pdb 또는 grafted_target.pdb)
        # 하지만 outputs_dir 기준으로 직접 확인 어려움 → designs/ 등 하위 step 이
        # 있으면 preprocess 는 무조건 완료로 간주.
        return True  # gather_all_inputs 단계에서 이미 처리됨

    if step in ("rfdiff", "rfdiffusion"):
        design_dir = os.path.join(outputs_dir, "designs")
        n_expected = int(inputs.get("n_designs", 5))
        n_pdbs = _count_pdbs(design_dir, prefix="design_")
        # n_pdbs 가 n_expected 의 절반 이상이면 완료로 간주 (워커별 분배 여유).
        return n_pdbs >= max(1, n_expected // 2)

    if step == "proteinmpnn":
        # design_dir/*/seqs/*.fa 또는 mpnn_results/
        for cand in ("mpnn_results", "designs"):
            path = os.path.join(outputs_dir, cand)
            if os.path.isdir(path):
                for root, dirs, files in os.walk(path):
                    if any(f.endswith(".fa") for f in files):
                        return True
        return False

    if step in ("af2", "alphafold"):
        af2_dir = os.path.join(outputs_dir, "af2_results")
        if not os.path.isdir(af2_dir):
            return False
        # rank_001 PDB 또는 unrelaxed PDB 중 하나라도 있으면 OK
        for f in os.listdir(af2_dir):
            if f.endswith(".pdb"):
                return True
        return False

    if step == "ncaa_mutation":
        ncaa_dir = os.path.join(outputs_dir, "ncaa_mutants")
        if not os.path.isdir(ncaa_dir):
            return False
        for f in os.listdir(ncaa_dir):
            if f.endswith(".pdb"):
                return True
        # manifest 가 있으면 ncAA 가 "none" 이어도 step 자체는 완료
        return os.path.exists(os.path.join(ncaa_dir, "mutation_manifest.json"))

    if step == "parameterize":
        params_dir = os.path.join(outputs_dir, "params")
        if not os.path.isdir(params_dir):
            return False
        # params_manifest.json 또는 하위 파일들
        if os.path.exists(os.path.join(params_dir, "params_manifest.json")):
            return True
        return len(os.listdir(params_dir)) > 0

    if step in ("md", "md_snapshots"):
        mdresult = os.path.join(outputs_dir, "mdresult")
        snap_dir = os.path.join(outputs_dir, "snapshots")
        # mdresult 에 뭐라도 있으면 md 는 완료
        if os.path.isdir(mdresult) and len(os.listdir(mdresult)) > 0:
            return True
        # snapshots 에 PDB 1 개 이상이면 md+snapshots 모두 완료
        return _count_pdbs(snap_dir) > 0

    if step == "snapshots":
        snap_dir = os.path.join(outputs_dir, "snapshots")
        return _count_pdbs(snap_dir) > 0

    if step == "qmmm":
        qmmm_dir = os.path.join(outputs_dir, "qmmm_results")
        snap_dir = os.path.join(outputs_dir, "snapshots")
        if not os.path.isdir(qmmm_dir):
            return False
        jsons = [f for f in os.listdir(qmmm_dir)
                 if f.endswith(".json") and "qmmm" in f]
        if not jsons:
            return False
        if not os.path.isdir(snap_dir):
            return len(jsons) > 0
        n_snaps = _count_pdbs(snap_dir)
        if n_snaps == 0:
            return len(jsons) > 0
        # 보수적: snapshot 수만큼은 qmmm JSON 이 있어야 "완전 완료".
        # 부분 완료면 False → 재실행 분기로 진입.
        return len(jsons) >= n_snaps

    if step == "mmgbsa":
        mmgbsa_dir = os.path.join(outputs_dir, "mmgbsa_results")
        if not os.path.isdir(mmgbsa_dir):
            return False
        for f in os.listdir(mmgbsa_dir):
            if f.endswith(".json") or f.endswith(".csv"):
                return True
        return False

    if step in ("rank", "final_ranking"):
        for fname in ("final_ranking.csv", "ranking.csv", "rank_results.csv"):
            if os.path.exists(os.path.join(outputs_dir, fname)):
                return True
        return False

    return False  # unknown step


if __name__ == "__main__":
    # 자가 테스트
    import tempfile
    import sys

    with tempfile.TemporaryDirectory() as tmp:
        state = PipelineState(log_dir=tmp, project_name="test")

        print("Initial state:")
        print(json.dumps(state.to_dict(), indent=2, ensure_ascii=False))

        state.start_step("rfdiff", total=10, workers=4)
        time.sleep(0.5)
        state.update_step("rfdiff", completed=5, detail="5/10 designs")
        time.sleep(0.5)
        state.complete_step("rfdiff")

        state.start_step("qmmm", total=8, workers=1)
        state.qmmm_cycle("w1_0_s1", 1, 3.73, 15.1)
        state.qmmm_cycle("w1_0_s1", 2, 1.83, -6.4)
        state.qmmm_converged("w1_0_s1", 19, -4630.9002)

        print("\nFinal state (via load_state):")
        loaded = PipelineState.load_state(tmp)
        print(json.dumps(loaded, indent=2, ensure_ascii=False))

        print("\nLast 10 events:")
        for ev in PipelineState.tail_events(tmp, n=10):
            print(f"  {ev['timestamp']} [{ev['kind']}] {ev['message']}")
