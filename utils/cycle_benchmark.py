#!/usr/bin/env python3
"""
cycle_benchmark.py — QM/MM SCF cycle 간격 측정 도구

로그 파일을 실시간 tail하면서 각 worker의 cycle 간 시간 간격을 측정.
Ctrl+C 시 per-worker + 전체 통계를 출력.

사용:
    python cycle_benchmark.py
    python cycle_benchmark.py --log /path/to/qmmm_live.log
    python cycle_benchmark.py --quiet  # 요약만 출력

실험 프로토콜:
    1. UPDD_QMMM_MAX_WORKERS=2 python UPDD.py
    2. 별도 창에서 이 스크립트 실행
    3. 5-10 cycle 관찰 후 Ctrl+C
    4. 워커 수 바꿔서 반복
"""

import os
import re
import sys
import time
import subprocess
import argparse
from datetime import datetime
from collections import defaultdict
import statistics

CYCLE_RE = re.compile(r"\[(\w+)\]\s+cycle=\s*(\d+)\s+E=")
CONVERGED_RE = re.compile(r"\[(\w+)\]\s+converged SCF energy")


def find_log():
    """기본 로그 경로 탐색."""
    home = os.path.expanduser("~")
    import glob
    candidates = glob.glob(
        os.path.join(home, "UPDD_proj", "log", "**", "qmmm_live.log"),
        recursive=True,
    )
    return max(candidates, key=os.path.getmtime) if candidates else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=str, default=None)
    parser.add_argument("--quiet", action="store_true",
                        help="요약 통계만 출력 (실시간 라인 숨김)")
    args = parser.parse_args()

    log_path = args.log or find_log()
    if not log_path or not os.path.exists(log_path):
        print(f"[!] 로그 없음: {log_path}")
        print("    UPDD.py가 먼저 실행 중이어야 합니다.")
        sys.exit(1)

    last_times = {}               # tag -> last cycle timestamp
    intervals = defaultdict(list) # tag -> [(cycle, interval_s)]
    converged = set()
    session_start = time.time()

    print(f"📊 Cycle Benchmark — {os.path.basename(log_path)}")
    print(f"   Start: {datetime.now().strftime('%H:%M:%S')}")
    print(f"   (Ctrl+C 로 종료 → 통계 출력)\n")
    if not args.quiet:
        print(f"{'Time':<10} {'Worker':<12} {'Cycle':<7} {'Interval':<10}")
        print("-" * 45)

    # -n 0 은 로그 끝부터 (과거 로그 제외)
    proc = subprocess.Popen(
        ["tail", "-F", "-n", "0", log_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )

    try:
        for line in proc.stdout:
            # cycle 라인
            m = CYCLE_RE.search(line)
            if m:
                tag = m.group(1)
                cycle = int(m.group(2))
                now = time.time()

                if tag in last_times:
                    dt = now - last_times[tag]
                    intervals[tag].append((cycle, dt))
                    if not args.quiet:
                        ts = datetime.now().strftime("%H:%M:%S")
                        print(f"{ts}   {tag:<12} c{cycle:<6} +{dt:>6.1f}s")
                else:
                    if not args.quiet:
                        ts = datetime.now().strftime("%H:%M:%S")
                        print(f"{ts}   {tag:<12} c{cycle:<6} (first)")

                last_times[tag] = now
                continue

            # 수렴 완료
            m = CONVERGED_RE.search(line)
            if m:
                tag = m.group(1)
                if tag not in converged:
                    converged.add(tag)
                    if not args.quiet:
                        ts = datetime.now().strftime("%H:%M:%S")
                        n = len(intervals.get(tag, []))
                        print(f"{ts}   {tag:<12} ✅ converged ({n} intervals recorded)")

    except KeyboardInterrupt:
        pass
    finally:
        proc.terminate()

    # ──────────── 통계 출력 ────────────
    session_elapsed = time.time() - session_start

    print("\n" + "=" * 60)
    print(f" 📊 Benchmark Summary  (session: {session_elapsed:.0f}s)")
    print("=" * 60)

    if not intervals:
        print("  (no cycle data collected)")
        return

    # per-worker
    print(f"\n  {'Worker':<12} {'N':<5} {'Avg':<8} {'Median':<8} "
          f"{'Min':<8} {'Max':<8} {'Stdev':<8}")
    print(f"  {'-'*12} {'-'*5} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")

    all_durations = []
    for tag in sorted(intervals):
        durs = [d for _, d in intervals[tag]]
        if not durs:
            continue
        all_durations.extend(durs)
        avg = statistics.mean(durs)
        med = statistics.median(durs)
        stdev = statistics.stdev(durs) if len(durs) > 1 else 0.0
        mark = "✅" if tag in converged else "  "
        print(f"  {mark}{tag:<10} {len(durs):<5} "
              f"{avg:<8.1f} {med:<8.1f} {min(durs):<8.1f} "
              f"{max(durs):<8.1f} {stdev:<8.2f}")

    # 전체
    print(f"\n  {'─'*12} 전체 ─")
    print(f"  총 cycle interval: {len(all_durations)}")
    print(f"  평균 cycle 시간: {statistics.mean(all_durations):.2f}s")
    print(f"  중앙값:          {statistics.median(all_durations):.2f}s")
    print(f"  Min / Max:       {min(all_durations):.1f}s / {max(all_durations):.1f}s")

    # throughput
    if all_durations:
        cycles_per_min = 60 / statistics.mean(all_durations) * len(intervals)
        print(f"\n  💡 throughput: ~{cycles_per_min:.1f} cycles/min "
              f"(across {len(intervals)} workers)")

    print("\n  이 결과를 다른 워커 수와 비교하세요:")
    print("     UPDD_QMMM_MAX_WORKERS=2 / 3 / 4 python UPDD.py")


if __name__ == "__main__":
    main()
