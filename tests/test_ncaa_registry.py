"""[v3 13-1] ncaa_registry 의 외부 API (resolve_ncaa_definition) 동작 검증.

검증 범위:
    - 알려진 code (NMA) 가 정확한 NCAADef 객체로 해석되는가
    - "none" / "None" 등 wild-type 별칭이 None 으로 해석되는가
    - 알 수 없는 code 가 ValueError 를 raise 하는가
    - import 시 _build_and_validate_registry 가 충돌 없이 통과하는가
    - import 시 _validate_smiles_if_rdkit_available 이 (RDKit 사용 가능 환경에서) 통과하는가
"""
import pytest

from ncaa_registry import (
    NCAA_REGISTRY_DATA,
    resolve_ncaa_definition,
    _build_and_validate_registry,
    get_all_ncaa_labels,
)


def test_resolve_known_code():
    result = resolve_ncaa_definition("NMA")
    assert result is not None
    assert result.code == "NMA"
    assert result.element == "C"
    # case-insensitive lookup도 통과해야 한다
    assert resolve_ncaa_definition("nma").code == "NMA"


def test_resolve_none_returns_none():
    assert resolve_ncaa_definition("none") is None
    assert resolve_ncaa_definition("None") is None
    assert resolve_ncaa_definition("") is None
    assert resolve_ncaa_definition("wild-type") is None


def test_unknown_raises():
    with pytest.raises(ValueError):
        resolve_ncaa_definition("INVALID_CODE_XYZ")


def test_registry_idempotent_validation():
    # 임포트 시점에 _build_and_validate_registry() 는 이미 1회 실행되었지만,
    # 재호출 시에도 ValueError 없이 통과해야 한다 (재진입 안전성).
    _build_and_validate_registry()


def test_all_labels_unique():
    labels = get_all_ncaa_labels()
    assert len(labels) == len(set(labels)), "ncAA 라벨이 중복됨"
    assert len(labels) == len(NCAA_REGISTRY_DATA)
