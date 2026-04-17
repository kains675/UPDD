#!/usr/bin/env python3
"""
updd_config.py — Step-independent worker configuration (v0.4)

v0.3까지는 parallel_workers가 모든 step에서 공유됨 (RFdiffusion, QM/MM, MM-GBSA 등).
이것이 "오늘 발견된 문제": QM/MM만 worker 수를 바꾸고 싶어도 RFdiffusion resume 로직이
재발동되어 무한 재생성을 일으킴.

v0.4 해결책: step별 독립 설정.
  workers: {
    rfdiff: 4,    # CPU/GPU 여유 있음
    af2:    1,    # AF2는 GPU 풀사용, 1이 최적
    qmmm:   1,    # RTX 5070 Ti 실측: 1이 최적 (1.62 cyc/min)
    mmgbsa: 4,    # CPU 병렬 안전
    md:     1     # OpenMM GPU
  }

Backward compat:
  - updd_status.json에 workers 키 없으면 parallel_workers로 fallback
  - 기존 프로젝트 그대로 작동
  - 새 프로젝트는 workers 사용

환경변수 override:
  UPDD_RFDIFF_WORKERS=2 python UPDD.py
  UPDD_QMMM_WORKERS=2 python UPDD.py
"""

import os
from typing import Dict, Optional

# Step 이름 → 기본값 (RTX 5070 Ti 16GB 기준)
DEFAULT_WORKERS: Dict[str, int] = {
    "rfdiff": 4,
    "af2": 1,
    "qmmm": 1,       # 오늘 실측으로 확정된 최적값
    "mmgbsa": 4,
    "md": 1,
}

# 환경변수 이름 매핑
ENV_VAR_MAP: Dict[str, str] = {
    "rfdiff": "UPDD_RFDIFF_WORKERS",
    "af2":    "UPDD_AF2_WORKERS",
    "qmmm":   "UPDD_QMMM_WORKERS",
    "mmgbsa": "UPDD_MMGBSA_WORKERS",
    "md":     "UPDD_MD_WORKERS",
}

# v0.3.3 시절 env var (deprecated but still honored)
LEGACY_ENV_VAR = "UPDD_QMMM_MAX_WORKERS"


def get_workers(step: str, inputs: Optional[dict] = None) -> int:
    """
    주어진 step에 대한 worker 수를 반환.

    우선순위 (높은 순):
      1. 환경변수 UPDD_{STEP}_WORKERS
      2. inputs["workers"][step] (v0.4 스키마)
      3. inputs["parallel_workers"] (v0.3 호환)
      4. DEFAULT_WORKERS[step]

    Args:
        step: "rfdiff", "af2", "qmmm", "mmgbsa", "md" 중 하나
        inputs: updd_status.json으로 로드된 dict (optional)

    Returns:
        worker 수 (1 이상의 정수)

    Example:
        >>> get_workers("qmmm", {"workers": {"qmmm": 2}})
        2
        >>> get_workers("rfdiff", {"parallel_workers": 4})  # legacy
        4
        >>> get_workers("md")  # no config
        1
    """
    if step not in DEFAULT_WORKERS:
        raise ValueError(f"Unknown step: {step}. "
                         f"Valid: {list(DEFAULT_WORKERS.keys())}")

    # 1. 환경변수 (최우선)
    env_var = ENV_VAR_MAP[step]
    env_val = os.environ.get(env_var)
    if env_val:
        try:
            return max(1, int(env_val))
        except ValueError:
            pass  # 무효한 값이면 무시하고 다음 우선순위로

    # 1.5. Legacy env var (v0.3.3 호환)
    if step == "qmmm":
        legacy = os.environ.get(LEGACY_ENV_VAR)
        if legacy:
            try:
                return max(1, int(legacy))
            except ValueError:
                pass

    # 2. inputs["workers"][step] (v0.4 스키마)
    if inputs and isinstance(inputs.get("workers"), dict):
        workers_dict = inputs["workers"]
        if step in workers_dict:
            try:
                return max(1, int(workers_dict[step]))
            except (ValueError, TypeError):
                pass

    # 3. inputs["parallel_workers"] (v0.3 fallback)
    if inputs and "parallel_workers" in inputs:
        try:
            return max(1, int(inputs["parallel_workers"]))
        except (ValueError, TypeError):
            pass

    # 4. 기본값
    return DEFAULT_WORKERS[step]


def resolve_all(inputs: Optional[dict] = None) -> Dict[str, int]:
    """모든 step의 worker 수를 한 번에 resolve."""
    return {step: get_workers(step, inputs) for step in DEFAULT_WORKERS}


def summary_line(inputs: Optional[dict] = None) -> str:
    """dashboard header용 한 줄 요약."""
    w = resolve_all(inputs)
    return " | ".join(f"{k}={v}" for k, v in w.items())


if __name__ == "__main__":
    # 자가 테스트
    print("=== updd_config.py self-test ===\n")

    print("1. 기본값 (no config):")
    for step in DEFAULT_WORKERS:
        print(f"   {step}: {get_workers(step)}")

    print("\n2. v0.3 legacy (parallel_workers=2):")
    legacy = {"parallel_workers": 2}
    for step in DEFAULT_WORKERS:
        print(f"   {step}: {get_workers(step, legacy)}")

    print("\n3. v0.4 스키마 (mixed):")
    v04 = {"workers": {"rfdiff": 4, "qmmm": 1}, "parallel_workers": 2}
    for step in DEFAULT_WORKERS:
        print(f"   {step}: {get_workers(step, v04)}  "
              f"# rfdiff/qmmm from workers, others fallback to parallel_workers")

    print("\n4. Summary:")
    print(f"   {summary_line(v04)}")
