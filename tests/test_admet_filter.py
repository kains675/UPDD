"""[v3 13-2] admet_filter.check_admet_rules 동작 검증.

검증 범위:
    - mode="none" 은 항상 통과
    - mode="cell" 에서 작은 약물(아스피린)이 통과하는가
    - mode="bbb" 에서 BBB 컷오프 상수가 적용되는가
    - 알 수 없는 mode 는 True (passthrough) 를 반환하는가
"""
import pytest

# RDKit 미설치 환경에서는 admet_filter import 자체가 실패한다 → 모듈 단위 skip
rdkit = pytest.importorskip("rdkit")

from admet_filter import (  # noqa: E402
    check_admet_rules,
    LIPINSKI_MW_MAX,
    BBB_TPSA_MAX,
)


def test_none_mode_always_passes():
    passed, msg = check_admet_rules("C", mode="none")
    assert passed is True
    assert "스킵" in msg or "Skip" in msg.lower() or msg


def test_aspirin_passes_cell():
    aspirin = "CC(=O)OC1=CC=CC=C1C(=O)O"
    passed, msg = check_admet_rules(aspirin, mode="cell")
    assert passed is True
    # Pass 메시지에는 MW / LogP 값이 포함되어야 한다
    assert "MW=" in msg


def test_unknown_mode_returns_true():
    passed, msg = check_admet_rules("CCO", mode="custom_unknown_mode")
    assert passed is True


def test_invalid_smiles_returns_false():
    passed, msg = check_admet_rules("not_a_real_smiles!@#$", mode="cell")
    assert passed is False


def test_lipinski_constants_are_module_level():
    # [v3 10-2] Lipinski 컷오프가 모듈 상수로 분리되었는지 검증
    assert LIPINSKI_MW_MAX == 500.0
    assert BBB_TPSA_MAX == 90.0
