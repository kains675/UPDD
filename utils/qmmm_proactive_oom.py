#!/usr/bin/env python3
"""
qmmm_proactive_oom.py — Proactive OOM Prevention with Convergence-Aware Priority (v0.3.3)

QM/MM 계산 중 VRAM 사용률이 임계치를 초과하면,
수렴에 가장 먼 워커를 선제 종료하여 수렴 직전인 워커를 보호한다.

v0.3.3 개선사항:
  - rich 기반 btop 스타일 고정 화면 대시보드
  - 상단 Header: 현재 UPDD 단계, 시작 시각, 경과 시간
  - VRAM 패널: 큰 게이지 + 상세
  - Workers 테이블: 실시간 cycle/|g|/남은 cycle
  - Events 로그: 최근 10개 이벤트 (HOMO warn, 수렴, kill 등)

사용법:
    python qmmm_proactive_oom.py                          # 기본 (VRAM 90%, 15초 간격)
    python qmmm_proactive_oom.py --vram-threshold 85      # 85%에서 트리거
    python qmmm_proactive_oom.py --dry-run                # kill 안 함 (모니터링만)

요구사항:
    pip install rich  (qmmm conda env에서)

Author: UPDD Team
Version: 0.3.3
"""

import os
import re
import sys
import time
import signal
import argparse
import subprocess
import numpy as np
from datetime import datetime, timedelta
from collections import deque

try:
    from rich.console import Console
    from rich.live import Live
    from rich.panel import Panel
    from rich.table import Table
    from rich.layout import Layout
    from rich.text import Text
    from rich.align import Align
    from rich import box
except ImportError:
    print("[!] rich 라이브러리 필요. 설치:")
    print("    conda run -n qmmm pip install rich")
    sys.exit(1)


# ──────────────────────────────────────────────────────────
# 1. VRAM 모니터링
# ──────────────────────────────────────────────────────────

def get_vram_usage():
    """pynvml 또는 nvidia-smi로 VRAM 사용률(%)과 사용량(MB)을 반환."""
    try:
        import pynvml
        pynvml.nvmlInit()
        handle = pynvml.nvmlDeviceGetHandleByIndex(0)
        info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        used_mb = info.used / (1024 ** 2)
        total_mb = info.total / (1024 ** 2)
        pct = info.used / info.total * 100
        pynvml.nvmlShutdown()
        return pct, used_mb, total_mb
    except ImportError:
        pass

    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.used,memory.total",
             "--format=csv,noheader,nounits"],
            text=True
        ).strip()
        used, total = [float(x) for x in out.split(",")]
        return used / total * 100, used, total
    except Exception:
        return None, None, None


# ──────────────────────────────────────────────────────────
# 2. QM/MM 워커 프로세스 탐색
# ──────────────────────────────────────────────────────────

def find_qmmm_workers():
    """실행 중인 run_qmmm.py 프로세스를 찾아 PID와 design_id를 반환.

    conda wrapper 제외, 실제 python 실행만 기록. tag로 dedup.
    """
    workers = []
    seen_tags = set()

    try:
        out = subprocess.check_output(
            ["ps", "aux"], text=True, stderr=subprocess.DEVNULL
        )
        for line in out.strip().split("\n"):
            if "run_qmmm.py" not in line or "grep" in line:
                continue

            if "conda run" in line:
                continue

            parts = line.split()
            pid = int(parts[1])

            design_id = "unknown"
            if "--filter" in line:
                idx = line.index("--filter")
                match = re.search(r"--filter\s+(\S+)", line[idx:])
                if match:
                    design_id = match.group(1)

            tag = design_id.replace("design_", "") if design_id.startswith("design_") else design_id

            if tag in seen_tags:
                continue
            seen_tags.add(tag)

            workers.append({"pid": pid, "design_id": design_id, "tag": tag})
    except Exception:
        pass
    return workers


def kill_worker_by_tag(tag, active_workers):
    """tag에 해당하는 워커를 SIGTERM으로 종료."""
    for w in active_workers:
        if w["tag"] == tag:
            try:
                os.kill(w["pid"], signal.SIGTERM)
                return True
            except ProcessLookupError:
                return False
            except PermissionError:
                return False
    return False


# ──────────────────────────────────────────────────────────
# 3. 로그 파일 자동 탐색
# ──────────────────────────────────────────────────────────

def find_latest_log():
    """~/UPDD_proj/log/**/qmmm_live.log 중 가장 최근 것."""
    import glob
    home = os.path.expanduser("~")
    candidates = glob.glob(os.path.join(home, "UPDD_proj", "log", "**", "qmmm_live.log"),
                           recursive=True)
    candidates += glob.glob(os.path.join(home, "UPDD_proj", "log", "**", "updd_*.log"),
                            recursive=True)
    if not candidates:
        return None
    return max(candidates, key=os.path.getmtime)


# ──────────────────────────────────────────────────────────
# 4. 수렴 추적
# ──────────────────────────────────────────────────────────

class ConvergenceTracker:
    """워커별 SCF 수렴 이력을 추적하고 kill 우선순위를 계산."""

    CYCLE_RE = re.compile(
        r"\[(\w+)\]\s+cycle=\s*(\d+)\s+E=\s*(\S+)\s+delta_E=\s*(\S+)"
        r"\s+\|g\|=\s*(\S+)\s+\|ddm\|=\s*(\S+)"
    )
    CONVERGED_RE = re.compile(r"\[(\w+)\]\s+converged SCF energy\s*=\s*([\d.eE+-]+)")
    HOMO_LUMO_WARN_RE = re.compile(r"\[(\w+)\]\s+WARN:\s+HOMO\s+.*==\s+LUMO")
    DESIGN_START_RE = re.compile(r"\[(\w+)\]\s+계산:\s*(design_\S+)")
    PASS_RE = re.compile(r"QM/MM Pass\s+(\d+)")
    MODE_RE = re.compile(r"\[(\w+)\].*(?:DFT계산|QM/MM Calculation)")
    MODE_QM_ONLY_RE = re.compile(r"\[(\w+)\].*QM 단독계산")
    TARGET_GRAD = 1e-4

    def __init__(self, event_log: deque):
        self.workers = {}
        self.log_position = 0
        self.event_log = event_log
        self.current_pass = 1
        self.current_designs = {}
        # [v0.3.4 fix #3] 전체 SCF 수렴 횟수 누적 (per-worker converged 플래그는
        # cycle=1 감지 시 리셋되므로, 세션 누적은 별도 카운터로 추적한다.)
        self.total_convergence_count = 0

    def update_from_log(self, log_path):
        if not os.path.exists(log_path):
            return

        try:
            with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
                f.seek(self.log_position)
                new_lines = f.readlines()
                self.log_position = f.tell()
        except IOError:
            return

        for line in new_lines:
            m = self.PASS_RE.search(line)
            if m:
                self.current_pass = int(m.group(1))
                continue

            m = self.DESIGN_START_RE.search(line)
            if m:
                tag = m.group(1)
                design = m.group(2)
                self.current_designs[tag] = design
                continue

            m = self.MODE_QM_ONLY_RE.search(line)
            if m:
                tag = m.group(1)
                if tag in self.workers:
                    self.workers[tag]["current_mode"] = "QM-only"
                else:
                    self.workers[tag] = {
                        "cycles": [],
                        "converged": False,
                        "homo_lumo_warns": 0,
                        "last_update": time.time(),
                        "current_mode": "QM-only",
                    }
                continue

            m = self.MODE_RE.search(line)
            if m:
                tag = m.group(1)
                if tag in self.workers:
                    self.workers[tag]["current_mode"] = "QM/MM"
                else:
                    self.workers[tag] = {
                        "cycles": [],
                        "converged": False,
                        "homo_lumo_warns": 0,
                        "last_update": time.time(),
                        "current_mode": "QM/MM",
                    }
                continue

            m = self.CYCLE_RE.search(line)
            if m:
                tag = m.group(1)
                cycle = int(m.group(2))
                delta_e = float(m.group(4))
                grad = float(m.group(5))

                if tag not in self.workers:
                    self.workers[tag] = {
                        "cycles": [],
                        "converged": False,
                        "homo_lumo_warns": 0,
                        "last_update": time.time(),
                        "current_mode": "\u2014",
                    }

                if cycle == 1 and (self.workers[tag]["converged"] or
                                   len(self.workers[tag]["cycles"]) > 0):
                    self.workers[tag]["converged"] = False
                    self.workers[tag]["cycles"] = []
                    self.workers[tag]["homo_lumo_warns"] = 0

                ddm = float(m.group(6))
                self.workers[tag]["cycles"].append((cycle, grad, delta_e, ddm))
                self.workers[tag]["last_update"] = time.time()
                continue

            m = self.CONVERGED_RE.search(line)
            if m:
                tag = m.group(1)
                if tag in self.workers:
                    if not self.workers[tag]["converged"]:
                        ts = datetime.now().strftime("%H:%M:%S")
                        ncycles = len(self.workers[tag]["cycles"])
                        self.event_log.append(
                            (ts, "✅", f"[{tag}] converged (cycles={ncycles})")
                        )
                        # [v0.3.4 fix #3] 새 SCF 수렴 1회 → 누적 카운트 증가.
                        # fix #1 에 의해 converged 플래그가 cycle=1 시 리셋되므로
                        # 이 분기(False→True 전이) 는 각 SCF 별 최초 1회만 진입한다.
                        self.total_convergence_count += 1
                    self.workers[tag]["converged"] = True
                continue

            m = self.HOMO_LUMO_WARN_RE.search(line)
            if m:
                tag = m.group(1)
                if tag in self.workers:
                    self.workers[tag]["homo_lumo_warns"] += 1
                    ts = datetime.now().strftime("%H:%M:%S")
                    warns = self.workers[tag]["homo_lumo_warns"]
                    self.event_log.append(
                        (ts, "⚠️ ", f"[{tag}] HOMO==LUMO warning (#{warns})")
                    )

    def estimate_remaining(self, tag):
        if tag not in self.workers:
            return float("inf")

        w = self.workers[tag]
        if w["converged"]:
            return 0

        cycles = w["cycles"]
        if len(cycles) < 3:
            return float("inf")

        recent = cycles[-10:]
        x = np.array([c[0] for c in recent], dtype=float)
        y = np.array([np.log10(max(c[1], 1e-15)) for c in recent])

        valid = np.isfinite(y)
        if valid.sum() < 2:
            return float("inf")

        x, y = x[valid], y[valid]

        try:
            a, b = np.polyfit(x, y, 1)
        except (np.linalg.LinAlgError, ValueError):
            return float("inf")

        y_std = np.std(y)
        if y_std < 0.05:
            return float("inf")

        if a >= -1e-3:
            return float("inf")

        target_y = np.log10(self.TARGET_GRAD)
        current_y = a * x[-1] + b

        if current_y <= target_y:
            return 0

        remaining = (target_y - current_y) / a

        if remaining > 500:
            return float("inf")

        return max(0, remaining)

    def classify_worker(self, tag):
        if tag not in self.workers:
            return "UNKNOWN"

        w = self.workers[tag]
        if w["converged"]:
            return "CONVERGED"

        cycles = w["cycles"]
        if len(cycles) < 10:
            return "UNKNOWN"

        all_grads = [c[1] for c in cycles]
        recent_grads = all_grads[-5:]
        current_g = recent_grads[-1]

        if current_g < 1e-2:
            return "CONVERGING"

        if current_g > 10 and current_g > recent_grads[0]:
            return "DIVERGING"

        total = len(cycles)
        if total > 30:
            warn_ratio = w["homo_lumo_warns"] / total
            if warn_ratio > 0.5:
                return "DIVERGING"

        if len(cycles) > 50 and current_g > 1.0:
            return "STALLED"

        return "CONVERGING"

    def get_kill_priority(self):
        priorities = []
        for tag, w in self.workers.items():
            if w["converged"]:
                continue

            remaining = self.estimate_remaining(tag)
            classification = self.classify_worker(tag)

            if classification == "DIVERGING":
                score = 1e9
            elif classification == "STALLED":
                score = 1e6
            elif classification == "UNKNOWN":
                score = 0
            else:
                score = remaining if remaining != float("inf") else 1e5

            priorities.append((tag, remaining, classification, score))

        priorities.sort(key=lambda x: x[3], reverse=True)
        return priorities


# ──────────────────────────────────────────────────────────
# 5. btop 스타일 대시보드 렌더링
# ──────────────────────────────────────────────────────────

def render_vram_bar(pct, width=40):
    """유니코드 블록으로 VRAM 바 렌더."""
    filled = int(pct / 100 * width)
    bar = "█" * filled + "░" * (width - filled)

    if pct < 60:
        color = "green"
    elif pct < 80:
        color = "yellow"
    else:
        color = "red"

    return f"[{color}]{bar}[/{color}]"


def build_header_panel(start_time, log_path, tracker, active_workers, dry_run, version="v0.3.3"):
    """상단 Header 패널."""
    now = datetime.now()
    elapsed = now - start_time
    elapsed_str = str(elapsed).split(".")[0]

    n_converged = sum(1 for w in tracker.workers.values() if w["converged"])
    n_total = len(tracker.workers) if tracker.workers else 0
    active_tags = [w["tag"] for w in active_workers]

    if active_tags:
        step_str = (f"QM/MM Pass {tracker.current_pass} │ "
                    f"Active: {', '.join(active_tags)} │ "
                    f"Converged: {n_converged}/{n_total}")
    else:
        step_str = f"대기 중 │ 총 {n_total}개 워커 기록됨"

    text = Text()
    text.append(f"⏱  Started: {start_time.strftime('%H:%M:%S')}    ", style="cyan")
    text.append(f"Elapsed: {elapsed_str}    ", style="bold cyan")
    text.append(f"Now: {now.strftime('%H:%M:%S')}\n", style="cyan")
    text.append(f"📍 {step_str}\n", style="white")
    text.append(f"📄 {log_path}", style="dim")
    if dry_run:
        text.append("    ", style="dim")
        text.append("[DRY-RUN]", style="bold yellow on black")

    return Panel(
        text,
        title=f"[bold cyan]🛡️  QM/MM Watchdog {version}[/]",
        title_align="left",
        border_style="cyan",
        box=box.ROUNDED,
    )


def build_vram_panel(pct, used, total, threshold, n_active, dry_run):
    """VRAM 패널."""
    bar = render_vram_bar(pct, width=40)

    if pct < 60:
        icon, color = "🟢", "green"
    elif pct < threshold:
        icon, color = "🟡", "yellow"
    else:
        icon, color = "🔴", "red"

    text = Text()
    text.append(f"{icon} ", style=color)
    text.append_text(Text.from_markup(bar))
    text.append("\n")
    text.append(f"   {pct:>5.1f}%  ({used:>7,.0f} / {total:>7,.0f} MB)   ",
                style=f"bold {color}")
    text.append(f"Threshold: {threshold:.0f}%   ", style="dim")
    text.append(f"Active: {n_active}   ", style="dim")
    mode_str = "dry-run" if dry_run else "live kill"
    mode_color = "yellow" if dry_run else "red"
    text.append(f"Mode: ", style="dim")
    text.append(mode_str, style=mode_color)

    return Panel(
        text,
        title="[bold]VRAM[/]",
        title_align="left",
        border_style=color,
        box=box.ROUNDED,
    )


def build_workers_table(tracker, active_workers):
    """Workers 테이블."""
    table = Table(
        show_header=True,
        header_style="bold white",
        box=box.SIMPLE,
        expand=True,
        padding=(0, 1),
    )
    table.add_column("●", width=2)
    table.add_column("TAG", style="cyan", min_width=10)
    table.add_column("MODE", style="magenta", width=9)
    table.add_column("STATUS", min_width=14)
    table.add_column("CYCLE", justify="right", min_width=5)
    table.add_column("|g|", justify="right", min_width=11)
    table.add_column("|ddm|", justify="right", min_width=9)
    table.add_column("Δ E", justify="right", min_width=11)
    table.add_column("REMAIN", justify="right", min_width=8)

    active_tags = {w["tag"] for w in active_workers}

    all_tags = set(tracker.workers.keys()) | active_tags
    if not all_tags:
        table.add_row("", "[dim]no workers yet[/]", "", "", "", "", "", "", "")
        return Panel(table, title="[bold]Workers[/]", title_align="left",
                     border_style="blue", box=box.ROUNDED)

    def sort_key(tag):
        w = tracker.workers.get(tag, {})
        if w.get("converged", False):
            return (3, tag)
        cls = tracker.classify_worker(tag)
        order = {"CONVERGING": 0, "UNKNOWN": 1, "STALLED": 2, "DIVERGING": 0}.get(cls, 1)
        return (order, tag)

    for tag in sorted(all_tags, key=sort_key):
        w = tracker.workers.get(tag, {"cycles": [], "converged": False})
        cycles = w.get("cycles", [])
        converged = w.get("converged", False)

        alive = "[green]●[/]" if tag in active_tags else "[dim]○[/]"

        mode_str = w.get("current_mode", "\u2014")

        if cycles:
            last = cycles[-1]
            last_cycle = last[0]
            last_grad = last[1]
            last_delta = last[2]
            last_ddm = last[3] if len(last) >= 4 else 0.0
            cycle_str = f"c{last_cycle}"
            grad_str = f"{last_grad:.2e}"
            ddm_str = f"{last_ddm:.2e}" if last_ddm else "\u2014"
            delta_str = f"{last_delta:+.1e}"
        else:
            cycle_str = "-"
            grad_str = "-"
            ddm_str = "\u2014"
            delta_str = "-"

        if converged:
            status_str = "[green]✅ CONVERGED[/]"
            remain_str = "[green]done[/]"
        else:
            cls = tracker.classify_worker(tag)
            remaining = tracker.estimate_remaining(tag)
            if cls == "CONVERGING":
                status_str = "[green]📉 CONVERGING[/]"
            elif cls == "DIVERGING":
                status_str = "[red]💥 DIVERGING[/]"
            elif cls == "STALLED":
                status_str = "[yellow]⏸  STALLED[/]"
            else:
                status_str = "[dim]❓ UNKNOWN[/]"

            if remaining == float("inf"):
                remain_str = "[red]∞[/]"
            elif remaining == 0:
                remain_str = "[green]≈0[/]"
            else:
                remain_str = f"~{remaining:.0f}"

        table.add_row(alive, tag, mode_str, status_str, cycle_str, grad_str, ddm_str, delta_str, remain_str)

    return Panel(
        table,
        title="[bold]Workers[/]",
        title_align="left",
        border_style="blue",
        box=box.ROUNDED,
    )


def build_events_panel(event_log, kill_history):
    """최근 이벤트 패널."""
    text = Text()

    events = list(event_log)[-10:]

    if not events:
        text.append("  (no events yet)", style="dim")
    else:
        for ts, icon, msg in reversed(events):
            text.append(f"  {ts}  {icon}  ", style="dim")
            text.append(f"{msg}\n")

    return Panel(
        text,
        title="[bold]Events[/]",
        title_align="left",
        border_style="magenta",
        box=box.ROUNDED,
    )


def build_dashboard(tracker, vram_info, active_workers, start_time, log_path,
                    threshold, dry_run, event_log, kill_history):
    """전체 대시보드 Layout."""
    layout = Layout()

    pct, used, total = vram_info

    layout.split(
        Layout(build_header_panel(start_time, log_path, tracker, active_workers, dry_run),
               size=6, name="header"),
        Layout(build_vram_panel(pct, used, total, threshold, len(active_workers), dry_run),
               size=5, name="vram"),
        Layout(build_workers_table(tracker, active_workers),
               name="workers"),
        Layout(build_events_panel(event_log, kill_history),
               size=14, name="events"),
    )

    return layout


# ──────────────────────────────────────────────────────────
# 6. 메인 루프
# ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Proactive OOM Prevention — btop 스타일 대시보드",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--vram-threshold", type=float, default=90.0)
    parser.add_argument("--interval", type=int, default=15,
                        help="갱신 간격 (초, 기본: 15)")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--log", type=str, default=None)
    parser.add_argument("--max-kills", type=int, default=3)
    args = parser.parse_args()

    log_path = args.log or find_latest_log()
    if not log_path:
        print("[!] UPDD 로그 파일을 찾을 수 없습니다.")
        print("    --log 옵션으로 직접 지정하거나, UPDD.py를 먼저 실행하세요.")
        sys.exit(1)

    console = Console()
    event_log = deque(maxlen=50)
    tracker = ConvergenceTracker(event_log)
    kill_count = 0
    kill_history = []
    start_time = datetime.now()

    # [v0.3.4 fix #2] no-worker grace period — 모든 worker 가 사라진 시점부터
    # 일정 시간(default 120s) 이 지나도록 재등장이 없으면 watchdog 자연 종료.
    # 이전 converged-기반 exit 은 fix #1 이 없을 때 "첫 SCF 수렴 → 즉시 종료"
    # 문제를 야기했다. pass 1→2 전환처럼 worker 가 잠시 사라졌다가 재등장할
    # 수 있으므로 단일 순간 폴링(active_workers empty) 만으로 종료하지 않는다.
    no_worker_since = None  # type: ignore[assignment]
    GRACE_PERIOD_SEC = 120

    event_log.append((start_time.strftime("%H:%M:%S"), "🛡️ ",
                      f"Watchdog started (threshold={args.vram_threshold}%, "
                      f"dry_run={args.dry_run})"))

    try:
        with Live(console=console, refresh_per_second=4, screen=True) as live:
            while True:
                tracker.update_from_log(log_path)

                vram_pct, vram_used, vram_total = get_vram_usage()
                if vram_pct is None:
                    vram_pct, vram_used, vram_total = 0.0, 0.0, 16384.0

                active_workers = find_qmmm_workers()
                active_tags = {w["tag"] for w in active_workers}

                layout = build_dashboard(
                    tracker,
                    (vram_pct, vram_used, vram_total),
                    active_workers,
                    start_time,
                    log_path,
                    args.vram_threshold,
                    args.dry_run,
                    event_log,
                    kill_history,
                )
                live.update(layout)

                if vram_pct >= args.vram_threshold:
                    priorities = tracker.get_kill_priority()
                    kill_candidates = [
                        (tag, remaining, cls, score)
                        for tag, remaining, cls, score in priorities
                        if tag in active_tags and cls != "CONVERGED"
                    ]

                    if kill_candidates and kill_count < args.max_kills:
                        target_tag, target_remaining, target_class, _ = kill_candidates[0]
                        ts = datetime.now().strftime("%H:%M:%S")

                        if args.dry_run:
                            event_log.append((ts, "🏷️ ",
                                            f"[DRY] would kill [{target_tag}] "
                                            f"({target_class}, VRAM {vram_pct:.1f}%)"))
                        else:
                            success = kill_worker_by_tag(target_tag, active_workers)
                            if success:
                                kill_count += 1
                                kill_history.append({
                                    "time": ts,
                                    "tag": target_tag,
                                    "classification": target_class,
                                    "vram_pct": vram_pct,
                                    "remaining_est": target_remaining,
                                })
                                event_log.append((ts, "🚨",
                                                f"KILLED [{target_tag}] "
                                                f"({target_class}, VRAM {vram_pct:.1f}%)"))

                # [v0.3.4 fix #2] Active worker 부재가 GRACE_PERIOD_SEC 이상
                # 지속될 때만 watchdog 종료. 이전 구현은 per-SCF convergence
                # 플래그만 보고 2 초 내 종료되는 버그가 있었다 (fix #1 설명 참고).
                if not active_workers:
                    if no_worker_since is None:
                        no_worker_since = time.time()
                        ts = datetime.now().strftime("%H:%M:%S")
                        event_log.append(
                            (ts, "⏸ ",
                             f"No active workers — grace period "
                             f"({GRACE_PERIOD_SEC}s) started"),
                        )
                    elif time.time() - no_worker_since > GRACE_PERIOD_SEC:
                        ts = datetime.now().strftime("%H:%M:%S")
                        event_log.append(
                            (ts, "🏁",
                             f"No workers for {GRACE_PERIOD_SEC}s — exiting"),
                        )
                        live.update(build_dashboard(
                            tracker, (vram_pct, vram_used, vram_total),
                            active_workers, start_time, log_path,
                            args.vram_threshold, args.dry_run,
                            event_log, kill_history,
                        ))
                        time.sleep(2)
                        break
                else:
                    # Worker 재등장 → grace period 리셋 (pass 전환 등 정상 케이스).
                    if no_worker_since is not None:
                        no_worker_since = None
                        ts = datetime.now().strftime("%H:%M:%S")
                        event_log.append(
                            (ts, "▶ ", "Workers resumed — grace period reset"),
                        )

                time.sleep(args.interval)

    except KeyboardInterrupt:
        pass

    # [v0.3.4 fix #3] Session Summary — "Converged" 카운트의 의미를 명확화.
    # fix #1 로 per-worker converged 플래그는 매 SCF 마다 리셋되므로,
    # 기존 방식은 "종료 시점에 마지막 SCF 가 converged 상태인 worker 수" 만
    # 세어 오해를 유발했다. 세션 누적 수렴 횟수는 tracker.total_convergence_count
    # 로 별도 추적한다.
    _n_still_active = sum(
        1 for w in tracker.workers.values() if not w.get("converged", False)
    )
    console.print()
    console.print(Panel(
        f"[bold]Session Summary[/]\n\n"
        f"  Total kills: {kill_count}\n"
        f"  Workers tracked: {len(tracker.workers)}\n"
        f"  Total SCF convergences: {tracker.total_convergence_count}\n"
        f"  In-flight at exit: {_n_still_active}\n"
        f"  Elapsed: {datetime.now() - start_time}",
        title="🛡️  Watchdog Ended",
        border_style="cyan",
    ))

    if kill_history:
        table = Table(title="Kill History", box=box.SIMPLE)
        table.add_column("Time")
        table.add_column("Tag", style="cyan")
        table.add_column("Class")
        table.add_column("VRAM %")
        table.add_column("Remain est")
        for k in kill_history:
            table.add_row(
                k["time"], k["tag"], k["classification"],
                f"{k['vram_pct']:.1f}",
                f"~{k['remaining_est']:.0f}" if k['remaining_est'] != float('inf') else "∞",
            )
        console.print(table)


if __name__ == "__main__":
    main()
