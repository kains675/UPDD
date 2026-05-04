"""[v3 13-3] rank_results_qmmm.minmax_normalize 동작 검증.

검증 범위:
    - 모든 값이 동일하면 0.5 균일 배열
    - invert=False: 최솟값=0.0, 최댓값=1.0
    - invert=True:  최솟값=1.0, 최댓값=0.0 (낮을수록 좋은 지표 정규화)
    - NaN 입력 처리 (np.nanmin/np.nanmax 사용 → NaN 항목 제외)
"""
import math
import numpy as np
import pytest

from rank_results_qmmm import minmax_normalize


def test_identical_values():
    result = minmax_normalize([5.0, 5.0, 5.0])
    assert all(v == 0.5 for v in result)


def test_normal_normalization():
    result = minmax_normalize([1.0, 2.0, 3.0], invert=False)
    assert math.isclose(result[0], 0.0)
    assert math.isclose(result[1], 0.5)
    assert math.isclose(result[2], 1.0)


def test_invert():
    # 낮을수록 좋은 지표는 invert=True 로 반전
    result = minmax_normalize([1.0, 2.0, 3.0], invert=True)
    assert math.isclose(result[0], 1.0)  # 최솟값이 1.0
    assert math.isclose(result[1], 0.5)
    assert math.isclose(result[2], 0.0)  # 최댓값이 0.0


def test_nan_passthrough():
    # NaN 항목은 출력에서도 NaN 으로 보존되어야 한다 (호출부가 가중치 재분배)
    arr = minmax_normalize([1.0, float("nan"), 3.0], invert=False)
    assert math.isclose(arr[0], 0.0)
    assert math.isnan(arr[1])
    assert math.isclose(arr[2], 1.0)
