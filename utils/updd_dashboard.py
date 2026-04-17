#!/usr/bin/env python3
"""
updd_dashboard.py — UPDD Pipeline Dashboard (v0.4, btop-style)

단일 명령 `python UPDD.py` 실행 시 이 대시보드가 메인 스레드에 뜨고,
실제 pipeline은 백그라운드 스레드에서 실행됩니다.

화면 구성 (4개 패널):
  ┌─ Header ─────────────────────────────────────┐
  │ Project name, elapsed, total progress bar   │
  └──────────────────────────────────────────────┘
  ┌─ Pipeline Steps ─────────────────────────────┐
  │ ✅ ⏭  🔄 ⏳ 각 step의 상태와 경과시간        │
  └──────────────────────────────────────────────┘
  ┌─ Current Step Detail ────────────────────────┐
  │ VRAM gauge + active workers (QM/MM 시)       │
  └──────────────────────────────────────────────┘
  ┌─ Events (최근 10개) ─────────────────────────┐
  │ 수렴/warn/kill 이벤트 실시간                │
  └──────────────────────────────────────────────┘

단독 실행도 가능:
    python utils/updd_dashboard.py --log-dir ~/UPDD_proj/log/6WGN
    (UPDD.py가 돌고 있는 동안 외부에서 관찰만 하기)

UPDD.py에서 사용:
    from updd_dashboard import run_dashboard_with
    run_dashboard_with(pipeline_fn=my_pipeline, log_dir=log_dir)
"""

import os
import sys
import time
import signal
import argparse
import threading
import subprocess
from pathlib import Path
from datetime import datetime
from collections import deque

try:
    from rich.console import Console
    from rich.live import Live
    from rich.panel import Panel
    from rich.table import Table
    from rich.layout import Layout
    from rich.text import Text
    from rich.progress import Progress, BarColumn, TextColumn, TaskProgressColumn
    from rich import box
except ImportError:
    print("[!] rich 라이브러리 필요: pip install rich")
    sys.exit(1)

# 같은 디렉토리의 updd_state 임포트
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from updd_state import PipelineState, StepStatus, PIPELINE_STEPS


# ──────────────────────────────────────────────────────────
# 1. VRAM 모니터링 (v0.3.3 watchdog에서 가져옴)
# ──────────────────────────────────────────────────────────

def get_vram():
    """Return (pct, used_mb, total_mb) or (None, None, None)."""
    try:
        import pynvml
        pynvml.nvmlInit()
        h = pynvml.nvmlDeviceGetHandleByIndex(0)
        info = pynvml.nvmlDeviceGetMemoryInfo(h)
        pynvml.nvmlShutdown()
        return (info.used / info.total * 100,
                info.used / 1024**2,
                info.total / 1024**2)
    except Exception:
        pass
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.used,memory.total",
             "--format=csv,noheader,nounits"],
            text=True,
        ).strip()
        used, total = [float(x) for x in out.split(",")]
        return used / total * 100, used, total
    except Exception:
        return None, None, None


def find_active_qmmm_workers():
    """Return list of {"pid", "tag"} for active run_qmmm.py processes."""
    workers = []
    try:
        pids = os.listdir("/proc")
    except OSError:
        return workers

    seen = set()
    for pid_str in pids:
        if not pid_str.isdigit():
            continue
        pid = int(pid_str)
        try:
            with open(f"/proc/{pid}/cmdline", "rb") as f:
                cmdline = f.read()
        except OSError:
            continue
        if b"run_qmmm.py" not in cmdline:
            continue
        args = [a for a in cmdline.decode("utf-8", errors="replace").split("\0") if a]
        if not args:
            continue
        exe = os.path.basename(args[0])
        if not exe.startswith("python"):
            continue  # skip conda wrapper

        tag = "unknown"
        for i, a in enumerate(args):
            if a == "--filter" and i + 1 < len(args):
                tag = args[i + 1].replace("design_", "")
                break
        if tag in seen:
            continue
        seen.add(tag)
        workers.append({"pid": pid, "tag": tag})
    return workers


# ──────────────────────────────────────────────────────────
# 2. Panel Builders
# ──────────────────────────────────────────────────────────

def build_header_panel(state_dict: dict, log_dir: str) -> Panel:
    elapsed = state_dict.get("elapsed", 0.0)
    hrs, rem = divmod(int(elapsed), 3600)
    mins, secs = divmod(rem, 60)
    elapsed_str = f"{hrs:02d}:{mins:02d}:{secs:02d}"

    project = state_dict.get("project", "(unknown)")
    total_pct = state_dict.get("total_progress", 0.0) * 100

    # 진행률 바 (수동)
    bar_width = 40
    filled = int(total_pct / 100 * bar_width)
    bar = "█" * filled + "░" * (bar_width - filled)

    n_steps = len(state_dict.get("steps", []))
    n_done = sum(1 for s in state_dict.get("steps", [])
                 if s["status"] in ("complete", "skipped"))

    text = Text()
    text.append("📦 ", style="cyan")
    text.append(f"{project}\n", style="bold cyan")
    text.append(f"⏱  Elapsed: ", style="dim")
    text.append(f"{elapsed_str}    ", style="bold")
    text.append(f"🏃 Now: ", style="dim")
    text.append(f"{datetime.now().strftime('%H:%M:%S')}\n", style="white")

    # progress bar — Text.append() 는 markup 을 파싱하지 않으므로
    # bar 문자열을 style 인자로 직접 적용한다 ([green]...[/green] leakage 방지).
    if total_pct < 40:
        color = "yellow"
    elif total_pct < 80:
        color = "cyan"
    else:
        color = "green"
    text.append("📊 ")
    text.append(bar, style=color)
    text.append(f"  {total_pct:>5.1f}%  ({n_done}/{n_steps} steps)")

    return Panel(
        text,
        title="[bold cyan]🧬  UPDD Pipeline Dashboard  v0.4[/]",
        title_align="left",
        border_style="cyan",
        box=box.ROUNDED,
    )


def build_steps_panel(state_dict: dict) -> Panel:
    steps = state_dict.get("steps", [])
    table = Table(show_header=False, box=box.SIMPLE, expand=True,
                  padding=(0, 1))
    table.add_column(width=3)   # 아이콘
    table.add_column(width=3)   # 번호
    table.add_column(min_width=18)  # 이름
    table.add_column(justify="right", min_width=10)  # elapsed
    table.add_column(min_width=20)  # detail

    icon_map = {
        "pending":  "[dim]⏳[/]",
        "running":  "[yellow]🔄[/]",
        "complete": "[green]✅[/]",
        "skipped":  "[dim]⏭ [/]",
        "failed":   "[red]❌[/]",
    }

    for i, s in enumerate(steps, 1):
        status = s["status"]
        icon = icon_map.get(status, "?")
        num = f"{i}."

        # name 스타일
        if status == "running":
            name = f"[bold yellow]{s['label']}[/]"
        elif status == "complete":
            name = f"[green]{s['label']}[/]"
        elif status == "failed":
            name = f"[red]{s['label']}[/]"
        else:
            name = f"[dim]{s['label']}[/]"

        elapsed_str = s.get("elapsed_str", "") if status != "pending" else "-"

        # detail
        if status == "running":
            if s["total"] > 0:
                detail = f"{s['completed']}/{s['total']}  (workers={s['workers']})"
            else:
                detail = s.get("detail") or "running..."
            detail = f"[cyan]{detail}[/]"
        elif status == "complete":
            detail = s.get("detail", "") or "done"
            detail = f"[dim]{detail}[/]"
        elif status == "skipped":
            detail = "[dim]already complete[/]"
        elif status == "failed":
            detail = f"[red]{s.get('error', 'error')}[/]"
        else:
            detail = ""

        table.add_row(icon, num, name, elapsed_str, detail)

    return Panel(
        table,
        title="[bold]Pipeline Steps[/]",
        title_align="left",
        border_style="blue",
        box=box.ROUNDED,
    )


def build_current_step_panel(state_dict: dict) -> Panel:
    """현재 step 상세 — QM/MM이면 VRAM + workers, 아니면 간단히."""
    current = None
    for s in state_dict.get("steps", []):
        if s["status"] == "running":
            current = s
            break

    text = Text()

    if not current:
        text.append("(대기 중 — 실행 중 step 없음)", style="dim")
        return Panel(text, title="[bold]Current Step[/]",
                     title_align="left", border_style="dim", box=box.ROUNDED)

    # VRAM (항상 표시, GPU step일 때 특히 유용)
    vram_pct, vram_used, vram_total = get_vram()
    if vram_pct is not None:
        bar_w = 36
        filled = int(vram_pct / 100 * bar_w)
        vram_bar = "█" * filled + "░" * (bar_w - filled)
        if vram_pct < 60:
            vcolor, vicon = "green", "🟢"
        elif vram_pct < 85:
            vcolor, vicon = "yellow", "🟡"
        else:
            vcolor, vicon = "red", "🔴"
        text.append(f"💻 VRAM  {vicon} ", style="")
        text.append(vram_bar, style=vcolor)
        text.append("  ")
        text.append(f"{vram_pct:>5.1f}%  "
                    f"({vram_used:>5.0f} / {vram_total:>5.0f} MB)\n\n",
                    style=f"bold {vcolor}")

    # 현재 step 정보
    text.append(f"🔄 ", style="yellow")
    text.append(f"{current['label']}  ", style="bold yellow")
    text.append(f"(workers={current['workers']}, "
                f"elapsed={current['elapsed_str']})\n", style="dim")
    if current.get("detail"):
        text.append(f"   {current['detail']}\n", style="cyan")

    # QM/MM step이면 활성 워커 표시
    if current["name"] == "qmmm":
        active = find_active_qmmm_workers()
        if active:
            text.append("\n")
            for w in active:
                text.append(f"   ● ", style="green")
                text.append(f"[{w['tag']}]", style="bold cyan")
                text.append(f"  pid={w['pid']}\n", style="dim")
            # 활성 워커 수도 detail 로 안내 (state.workers 와 교차 확인 가능)
            text.append(f"   (active subprocesses: {len(active)})\n", style="dim")
        else:
            text.append("\n   [dim]대기 중... (워커 시작 전)[/]\n", style="")

    return Panel(
        text,
        title="[bold]Current Step Details[/]",
        title_align="left",
        border_style="yellow" if current else "dim",
        box=box.ROUNDED,
    )


def build_events_panel(events: list) -> Panel:
    text = Text()
    if not events:
        text.append("  (이벤트 없음)", style="dim")
    else:
        # 최신이 위
        for ev in reversed(events[-10:]):
            ts = ev.get("timestamp", "")
            kind = ev.get("kind", "")
            msg = ev.get("message", "")

            # kind 별 색상
            color_map = {
                "step_start": "yellow",
                "step_complete": "green",
                "step_fail": "red",
                "qmmm_converged": "green",
                "qmmm_warn": "yellow",
                "warn": "yellow",
                "error": "red",
            }
            color = color_map.get(kind, "white")

            text.append(f"  {ts}  ", style="dim")
            text.append(f"{msg}\n", style=color)

    return Panel(
        text,
        title="[bold]Events[/] (최근 10개)",
        title_align="left",
        border_style="magenta",
        box=box.ROUNDED,
    )


# ──────────────────────────────────────────────────────────
# 3. Dashboard main loop
# ──────────────────────────────────────────────────────────

def build_footer_panel() -> Panel:
    """정적 단축키 힌트 (Phase 3 키입력은 v0.4.1에서 보강 예정)."""
    text = Text()
    text.append("  Press ", style="dim")
    text.append("Ctrl+C", style="bold yellow")
    text.append(" to stop pipeline  ", style="dim")
    text.append(" · ", style="dim")
    text.append("UPDD v0.4", style="bold cyan")
    return Panel(text, border_style="dim", box=box.SIMPLE)


def build_layout(state_dict: dict, events: list, log_dir: str) -> Layout:
    layout = Layout()
    layout.split(
        Layout(build_header_panel(state_dict, log_dir), size=6, name="header"),
        Layout(build_steps_panel(state_dict), size=15, name="steps"),
        Layout(build_current_step_panel(state_dict), name="current"),
        Layout(build_events_panel(events), size=14, name="events"),
        Layout(build_footer_panel(), size=3, name="footer"),
    )
    return layout


def run_dashboard(log_dir: str, stop_event: threading.Event = None,
                  refresh_hz: int = 4, console: Console = None):
    """
    Dashboard 메인 루프. stop_event가 set되면 종료.

    단독 실행 시 Ctrl+C로 종료.
    UPDD.py와 함께 쓸 때는 pipeline thread 완료 후 stop_event.set() 호출.

    Args:
        console: 외부에서 캡처한 Console 객체 (UPDD.py가 stdout 리다이렉트하기
                 전에 `Console(file=sys.__stdout__, force_terminal=True)` 로
                 캡처해 전달하면, dashboard 출력이 pipeline.log 로 새지 않는다.
                 None 이면 원본 stdout 에 바인딩된 Console 을 생성한다.
    """
    if console is None:
        console = Console(file=sys.__stdout__, force_terminal=True)
    log_dir = os.path.expanduser(log_dir)

    # 시작 시 state 파일이 생성되길 잠깐 기다림
    wait_time = 0
    while not os.path.exists(os.path.join(log_dir, "state.json")):
        if stop_event and stop_event.is_set():
            return
        time.sleep(0.2)
        wait_time += 0.2
        if wait_time > 10:
            console.print(f"[yellow]state.json이 생성되지 않음: {log_dir}[/]")
            return

    try:
        with Live(console=console, refresh_per_second=refresh_hz,
                  screen=True, vertical_overflow="visible") as live:
            while True:
                if stop_event and stop_event.is_set():
                    break

                state_dict = PipelineState.load_state(log_dir)
                if state_dict is None:
                    state_dict = {"steps": [], "elapsed": 0.0,
                                  "total_progress": 0.0, "project": "(loading)"}

                events = PipelineState.tail_events(log_dir, n=10)
                layout = build_layout(state_dict, events, log_dir)
                live.update(layout)

                time.sleep(1.0 / refresh_hz)
    except KeyboardInterrupt:
        pass


def run_dashboard_in_thread(log_dir: str, console: Console = None) -> tuple:
    """
    대시보드를 백그라운드 스레드에서 실행 (UPDD.py 통합용).

    Args:
        console: 외부에서 미리 캡처한 Console (stdout 리다이렉트 방지용).
                 run_dashboard() 의 docstring 참조.

    Returns:
        (thread, stop_event) - UPDD.py 종료 시 stop_event.set() 필요
    """
    stop_event = threading.Event()
    t = threading.Thread(
        target=run_dashboard,
        args=(log_dir, stop_event, 4, console),
        daemon=True,  # UPDD 종료 시 강제 종료
    )
    t.start()
    return t, stop_event


# ──────────────────────────────────────────────────────────
# 4. CLI (단독 실행)
# ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="UPDD Pipeline Dashboard (standalone mode)")
    parser.add_argument("--log-dir", type=str, default=None,
                        help="로그 디렉토리 (자동 탐색 시 생략)")
    parser.add_argument("--refresh-hz", type=int, default=4,
                        help="갱신 주파수 (Hz, 기본 4)")
    args = parser.parse_args()

    log_dir = args.log_dir
    if not log_dir:
        # 자동 탐색: ~/UPDD_proj/log/<project>/
        import glob
        home = os.path.expanduser("~")
        candidates = glob.glob(os.path.join(home, "UPDD_proj", "log", "*"))
        candidates = [c for c in candidates if os.path.isdir(c)]
        if not candidates:
            print("[!] 로그 디렉토리를 찾을 수 없음. --log-dir 지정 필요.")
            sys.exit(1)
        log_dir = max(candidates, key=os.path.getmtime)
        print(f"자동 탐색된 로그 디렉토리: {log_dir}")
        time.sleep(1)

    run_dashboard(log_dir=log_dir, refresh_hz=args.refresh_hz)


if __name__ == "__main__":
    main()
