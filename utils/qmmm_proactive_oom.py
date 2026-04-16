#!/usr/bin/env python3
"""
qmmm_proactive_oom.py — Proactive OOM Prevention with Convergence-Aware Priority

QM/MM 계산 중 VRAM 사용률이 임계치를 초과하면,
수렴에 가장 먼 워커를 선제 종료하여 수렴 직전인 워커를 보호한다.

사용법:
    # UPDD.py가 QM/MM을 실행 중일 때, 별도 터미널(또는 tmux 창)에서:
    python qmmm_proactive_oom.py

    # 옵션:
    python qmmm_proactive_oom.py --vram-threshold 88 --interval 15 --dry-run

    # tmux 분할 창에서 모니터링:
    # 창 0: python UPDD.py          (QM/MM 실행)
    # 창 1: python qmmm_proactive_oom.py  (OOM 감시)

v0.3에서 UPDD.py에 통합 예정. 현재는 독립 실행 스크립트로 운용.

Author: UPDD Team
Version: 0.1 (standalone watchdog)
"""

import os
import re
import sys
import time
import signal
import argparse
import subprocess
import numpy as np
from datetime import datetime
from collections import defaultdict

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

    # fallback: nvidia-smi
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.used,memory.total",
             "--format=csv,noheader,nounits"],
            text=True
        ).strip()
        used, total = [float(x) for x in out.split(",")]
        return used / total * 100, used, total
    except Exception as e:
        print(f"  [!] VRAM 조회 실패: {e}")
        return None, None, None


# ──────────────────────────────────────────────────────────
# 2. QM/MM 워커 프로세스 탐색
# ──────────────────────────────────────────────────────────

def find_qmmm_workers():
    """실행 중인 run_qmmm.py 프로세스를 찾아 PID와 design_id를 반환.

    Returns:
        list[dict]: [{"pid": int, "design_id": str, "tag": str}, ...]
    """
    workers = []
    try:
        # ps aux에서 run_qmmm.py 프로세스 탐색
        out = subprocess.check_output(
            ["ps", "aux"], text=True, stderr=subprocess.DEVNULL
        )
        for line in out.strip().split("\n"):
            if "run_qmmm.py" not in line or "grep" in line:
                continue
            parts = line.split()
            pid = int(parts[1])

            # --filter 인자에서 design_id 추출
            design_id = "unknown"
            if "--filter" in line:
                idx = line.index("--filter")
                # --filter 뒤의 값 추출
                match = re.search(r"--filter\s+(\S+)", line[idx:])
                if match:
                    design_id = match.group(1)

            # 짧은 태그
            tag = design_id.replace("design_", "").split("_unrelaxed")[0]
            if len(tag) > 10:
                tag = tag[:10]

            workers.append({"pid": pid, "design_id": design_id, "tag": tag})
    except Exception as e:
        print(f"  [!] 프로세스 탐색 실패: {e}")

    return workers


# ──────────────────────────────────────────────────────────
# 3. SCF 수렴 추적 (로그 파싱)
# ──────────────────────────────────────────────────────────

class ConvergenceTracker:
    """UPDD.py 로그를 실시간 파싱하여 워커별 SCF 수렴 상태를 추적."""

    # 로그 패턴: [w4_1_s2] cycle= 31 E= -6379.52100  delta_E= -6.46e-06  |g|= 0.000997  |ddm|= 0.056
    CYCLE_RE = re.compile(
        r"\[(\w+)\]\s+cycle=\s*(\d+)\s+E=\s*([\d.eE+-]+)\s+"
        r"delta_E=\s*([\d.eE+-]+)\s+\|g\|=\s*([\d.eE+-]+)\s+\|ddm\|=\s*([\d.eE+-]+)"
    )

    # 수렴 완료: [w4_1_s2] converged SCF energy = ...
    CONVERGED_RE = re.compile(r"\[(\w+)\]\s+converged SCF energy")

    # 발산/실패: [w4_1_s1] WARN: HOMO ... == LUMO
    HOMO_LUMO_WARN_RE = re.compile(r"\[(\w+)\]\s+WARN:\s+HOMO\s+.*==\s+LUMO")

    # QM-MM 결과: [w4_1_s2]   QM-MM 상호작용: -1.9303 kcal/mol
    RESULT_RE = re.compile(r"\[(\w+)\].*QM-MM 상호작용:\s*([\d.eE+-]+)\s*kcal/mol")

    TARGET_GRAD = 1e-4  # |g| 수렴 목표

    def __init__(self):
        self.workers = {}  # tag → {"cycles": [(cycle, grad, delta_e)], "converged": bool, ...}
        self.log_position = 0  # 로그 파일 읽기 위치

    def update_from_log(self, log_path):
        """로그 파일의 새로운 라인을 파싱하여 수렴 상태를 업데이트."""
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
            # cycle 데이터
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
                    }
                self.workers[tag]["cycles"].append((cycle, grad, delta_e))
                self.workers[tag]["last_update"] = time.time()
                continue

            # 수렴 완료
            m = self.CONVERGED_RE.search(line)
            if m:
                tag = m.group(1)
                if tag in self.workers:
                    self.workers[tag]["converged"] = True
                continue

            # HOMO==LUMO 경고 (발산 신호)
            m = self.HOMO_LUMO_WARN_RE.search(line)
            if m:
                tag = m.group(1)
                if tag in self.workers:
                    self.workers[tag]["homo_lumo_warns"] += 1

    def estimate_remaining(self, tag):
        """log₁₀(|g|) 선형 회귀로 수렴까지 남은 cycle 수를 추정.

        Returns:
            float: 추정 남은 cycle 수.
                   양수 = 수렴 중, 음수/inf = 발산 또는 추정 불가.
        """
        if tag not in self.workers:
            return float("inf")

        w = self.workers[tag]
        if w["converged"]:
            return 0

        cycles = w["cycles"]
        if len(cycles) < 3:
            return float("inf")  # 데이터 부족

        # 최근 10개 cycle로 추세 파악
        recent = cycles[-10:]
        x = np.array([c[0] for c in recent], dtype=float)  # cycle number
        y = np.array([np.log10(max(c[1], 1e-15)) for c in recent])  # log₁₀(|g|)

        # |g|가 0이거나 음수인 비정상 값 필터
        valid = np.isfinite(y)
        if valid.sum() < 2:
            return float("inf")

        x, y = x[valid], y[valid]

        # 선형 회귀: y = a*x + b
        try:
            a, b = np.polyfit(x, y, 1)
        except (np.linalg.LinAlgError, ValueError):
            return float("inf")

        # 기울기가 양수면 발산 중
        if a >= 0:
            return float("inf")

        # 목표: log₁₀(0.0001) = -4
        target_y = np.log10(self.TARGET_GRAD)
        current_y = a * x[-1] + b

        if current_y <= target_y:
            return 0  # 이미 수렴 근처

        remaining = (target_y - current_y) / a
        return max(0, remaining)

    def classify_worker(self, tag):
        """워커를 분류: CONVERGED, CONVERGING, STALLED, DIVERGING, UNKNOWN."""
        if tag not in self.workers:
            return "UNKNOWN"

        w = self.workers[tag]
        if w["converged"]:
            return "CONVERGED"

        cycles = w["cycles"]
        if len(cycles) < 3:
            return "UNKNOWN"

        # 최근 5개 |g| 추세
        recent_grads = [c[1] for c in cycles[-5:]]

        # 발산 판정: |g| > 10이고 증가 추세
        if recent_grads[-1] > 10 and recent_grads[-1] > recent_grads[0]:
            return "DIVERGING"

        # HOMO==LUMO 경고가 3회 이상
        if w["homo_lumo_warns"] >= 3:
            return "DIVERGING"

        # 정체 판정: |g|가 줄지 않고 진동
        if len(cycles) > 50 and recent_grads[-1] > 1.0:
            return "STALLED"

        return "CONVERGING"

    def get_kill_priority(self):
        """kill 우선순위 리스트를 반환 (높은 순서 = 먼저 kill).

        Returns:
            list[tuple]: [(tag, remaining_cycles, classification), ...] 정렬됨
        """
        priorities = []
        for tag, w in self.workers.items():
            if w["converged"]:
                continue  # 수렴 완료된 건 제외

            remaining = self.estimate_remaining(tag)
            classification = self.classify_worker(tag)

            # 우선순위 점수 (높을수록 먼저 kill)
            if classification == "DIVERGING":
                score = 999999
            elif classification == "STALLED":
                score = 99999
            elif remaining == float("inf"):
                score = 88888
            else:
                score = remaining

            priorities.append((tag, remaining, classification, score))

        # 점수 높은 순 (= 가장 먼저 kill할 후보)
        priorities.sort(key=lambda x: -x[3])
        return priorities


# ──────────────────────────────────────────────────────────
# 4. 선제 종료 실행
# ──────────────────────────────────────────────────────────

def kill_worker_by_tag(tag, workers_info):
    """태그에 해당하는 워커 프로세스를 SIGTERM으로 종료.

    Args:
        tag: 워커 태그 (예: "w4_1_s1")
        workers_info: find_qmmm_workers()의 결과

    Returns:
        bool: 종료 성공 여부
    """
    for w in workers_info:
        if w["tag"] == tag:
            pid = w["pid"]
            try:
                os.kill(pid, signal.SIGTERM)
                print(f"  🔪 [KILL] PID {pid} ({tag}) — SIGTERM 전송")
                return True
            except ProcessLookupError:
                print(f"  [!] PID {pid} ({tag}) — 이미 종료됨")
                return False
            except PermissionError:
                print(f"  [!] PID {pid} ({tag}) — 권한 부족 (sudo 필요?)")
                return False
    print(f"  [!] 태그 '{tag}'에 해당하는 워커를 찾을 수 없음")
    return False


# ──────────────────────────────────────────────────────────
# 5. 메인 감시 루프
# ──────────────────────────────────────────────────────────

def find_latest_log():
    """가장 최근 UPDD 로그 파일을 찾는다."""
    log_dirs = [
        os.path.expanduser("~/UPDD_proj/log/6WGN"),
        os.path.expanduser("~/UPDD_proj/log"),
    ]
    latest = None
    latest_mtime = 0

    for log_dir in log_dirs:
        if not os.path.isdir(log_dir):
            continue
        for f in os.listdir(log_dir):
            if f.startswith("updd_") and f.endswith(".log"):
                path = os.path.join(log_dir, f)
                mtime = os.path.getmtime(path)
                if mtime > latest_mtime:
                    latest_mtime = mtime
                    latest = path

    return latest


def main():
    parser = argparse.ArgumentParser(
        description="QM/MM Proactive OOM Prevention — 수렴 인지 선제 종료",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
사용 예시:
  python qmmm_proactive_oom.py                    # 기본 설정 (VRAM 90%, 30초 간격)
  python qmmm_proactive_oom.py --vram-threshold 85  # 85%에서 트리거
  python qmmm_proactive_oom.py --dry-run             # 실제 kill 없이 모니터링만
  python qmmm_proactive_oom.py --log /path/to/updd.log  # 특정 로그 파일 지정
        """
    )
    parser.add_argument("--vram-threshold", type=float, default=90.0,
                        help="VRAM 사용률 임계치 (%%, 기본: 90)")
    parser.add_argument("--interval", type=int, default=30,
                        help="모니터링 간격 (초, 기본: 30)")
    parser.add_argument("--dry-run", action="store_true",
                        help="실제 kill 없이 모니터링만 수행")
    parser.add_argument("--log", type=str, default=None,
                        help="UPDD 로그 파일 경로 (자동 탐색 시 생략)")
    parser.add_argument("--max-kills", type=int, default=3,
                        help="한 세션에서 최대 kill 횟수 (기본: 3)")

    args = parser.parse_args()

    # 로그 파일 탐색
    log_path = args.log or find_latest_log()
    if not log_path:
        print("[!] UPDD 로그 파일을 찾을 수 없습니다.")
        print("    --log 옵션으로 직접 지정하거나, UPDD.py를 먼저 실행하세요.")
        sys.exit(1)

    print("=" * 60)
    print(" 🛡️  QM/MM Proactive OOM Prevention Watchdog")
    print("=" * 60)
    print(f"  로그 파일    : {log_path}")
    print(f"  VRAM 임계치  : {args.vram_threshold:.0f}%")
    print(f"  모니터링 간격: {args.interval}초")
    print(f"  Dry-run 모드 : {'✅ ON (kill 안 함)' if args.dry_run else '❌ OFF (실제 kill)'}")
    print(f"  최대 kill 수 : {args.max_kills}")
    print("=" * 60)

    tracker = ConvergenceTracker()
    kill_count = 0
    kill_history = []

    try:
        while True:
            # ── A. 로그 파싱 ──
            tracker.update_from_log(log_path)

            # ── B. VRAM 체크 ──
            vram_pct, vram_used, vram_total = get_vram_usage()
            if vram_pct is None:
                time.sleep(args.interval)
                continue

            # ── C. 현재 상태 출력 ──
            now = datetime.now().strftime("%H:%M:%S")
            vram_bar_len = int(vram_pct / 5)
            vram_bar = "█" * vram_bar_len + "░" * (20 - vram_bar_len)
            vram_color = "🟢" if vram_pct < 70 else ("🟡" if vram_pct < args.vram_threshold else "🔴")

            print(f"\n[{now}] VRAM {vram_color} {vram_bar} {vram_pct:.1f}% "
                  f"({vram_used:.0f}/{vram_total:.0f} MB)")

            # 워커별 수렴 상태
            active_workers = find_qmmm_workers()
            active_tags = {w["tag"] for w in active_workers}
            priorities = tracker.get_kill_priority()

            if priorities:
                print(f"  {'TAG':<12} {'STATUS':<12} {'Cycle':<8} {'|g|':<12} "
                      f"{'Remaining':<12} {'Kill순위'}")
                print(f"  {'─'*12} {'─'*12} {'─'*8} {'─'*12} {'─'*12} {'─'*8}")

                for i, (tag, remaining, classification, score) in enumerate(priorities):
                    w = tracker.workers.get(tag, {})
                    cycles = w.get("cycles", [])
                    last_cycle = cycles[-1][0] if cycles else 0
                    last_grad = cycles[-1][1] if cycles else 0

                    # 활성 상태 표시
                    alive = "🟢" if tag in active_tags else "⚪"

                    # 남은 cycle 표시
                    if remaining == float("inf"):
                        rem_str = "∞"
                    elif remaining == 0:
                        rem_str = "수렴!"
                    else:
                        rem_str = f"~{remaining:.0f} cyc"

                    # classification 이모지
                    cls_emoji = {
                        "CONVERGING": "📉",
                        "CONVERGED": "✅",
                        "STALLED": "⏸️",
                        "DIVERGING": "💥",
                        "UNKNOWN": "❓",
                    }.get(classification, "❓")

                    print(f"  {alive}{tag:<11} {cls_emoji}{classification:<11} "
                          f"c{last_cycle:<7} {last_grad:<12.6f} {rem_str:<12} #{i+1}")

            # ── D. 선제 종료 판단 ──
            if vram_pct >= args.vram_threshold and priorities:
                # kill 대상: 활성 워커 중 kill 우선순위 최상위
                kill_candidates = [
                    (tag, remaining, classification, score)
                    for tag, remaining, classification, score in priorities
                    if tag in active_tags and classification != "CONVERGED"
                ]

                if kill_candidates and kill_count < args.max_kills:
                    target_tag = kill_candidates[0][0]
                    target_remaining = kill_candidates[0][1]
                    target_class = kill_candidates[0][2]

                    print(f"\n  🚨 VRAM {vram_pct:.1f}% ≥ {args.vram_threshold:.0f}% — "
                          f"선제 종료 트리거!")
                    print(f"  🎯 대상: [{target_tag}] "
                          f"({target_class}, ~{target_remaining:.0f} cycles 남음)")

                    if args.dry_run:
                        print(f"  🏷️  [DRY-RUN] kill 건너뜀 — 실제 모드에서는 종료됩니다")
                    else:
                        success = kill_worker_by_tag(target_tag, active_workers)
                        if success:
                            kill_count += 1
                            kill_history.append({
                                "time": now,
                                "tag": target_tag,
                                "classification": target_class,
                                "vram_pct": vram_pct,
                                "remaining_est": target_remaining,
                            })
                            print(f"  ✅ 종료 성공 (누적 {kill_count}/{args.max_kills})")

                            # 수렴 직전 워커 보호 성공 알림
                            protected = [
                                (t, r) for t, r, c, s in kill_candidates[1:]
                                if c == "CONVERGING" and r < 50
                            ]
                            if protected:
                                names = ", ".join(t for t, _ in protected)
                                print(f"  🛡️  보호된 워커: {names}")

                elif kill_count >= args.max_kills:
                    print(f"\n  ⚠️  VRAM {vram_pct:.1f}% 초과했지만 "
                          f"최대 kill 횟수({args.max_kills}) 도달 — 대기")

            # ── E. 모든 워커 수렴 완료 시 종료 ──
            if active_workers and all(
                tracker.workers.get(w["tag"], {}).get("converged", False)
                for w in active_workers
            ):
                print("\n  🎉 모든 활성 워커 수렴 완료! 감시 종료.")
                break

            if not active_workers and tracker.workers:
                print("\n  ⚪ 활성 QM/MM 워커 없음 — 대기 중...")

            time.sleep(args.interval)

    except KeyboardInterrupt:
        print("\n\n  🛑 수동 중단")

    # ── 종료 요약 ──
    print("\n" + "=" * 60)
    print(" 📊 Watchdog 세션 요약")
    print("=" * 60)
    print(f"  총 kill 횟수: {kill_count}")
    if kill_history:
        for k in kill_history:
            print(f"    {k['time']} [{k['tag']}] {k['classification']} "
                  f"(VRAM {k['vram_pct']:.1f}%, ~{k['remaining_est']:.0f} cyc 남음)")

    # 최종 수렴 상태
    print(f"\n  워커 최종 상태:")
    for tag, w in sorted(tracker.workers.items()):
        cycles = w["cycles"]
        status = "✅ 수렴" if w["converged"] else tracker.classify_worker(tag)
        last_grad = cycles[-1][1] if cycles else "N/A"
        print(f"    [{tag}] {status} — {len(cycles)} cycles, last |g|={last_grad}")


if __name__ == "__main__":
    main()
