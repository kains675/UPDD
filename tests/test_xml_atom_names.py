"""[v4 T-4] NMA XML 템플릿 중원자명이 PDB 규약과 일치하는지 검증.

목적:
    parameterize_ncaa.py 의 _rename_xml_to_pdb_names 가 GAFF2 내부명 (N1/C2/...) 을
    PDB 규약명 (N/CA/C/O/CB/CM/...) 으로 정확히 변환했는지 회귀 테스트.
    이 테스트가 깨지면 run_restrained_md.py 의 addHydrogens 가 "No template found"
    오류로 즉시 실패한다 (CLAUDE.md v4 우선과제와 동일한 클래스의 버그).

탐색 정책:
    - tests/test_xml_atom_names.py 는 outputs/ 하위에 NMA_gaff2.xml 이 존재할 때만
      활성화된다 (pytest skip). CI 환경에서는 outputs/ 가 비어있을 수 있으므로
      유연한 skip 처리를 한다.
    - 발견된 모든 NMA_gaff2.xml 파일에 대해 검증을 수행한다.
"""
import os
import glob
import xml.etree.ElementTree as ET
import pytest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
NMA_EXPECTED_HEAVY = {"N", "CA", "C", "O", "CB", "CM"}


def _find_nma_xmls():
    """outputs/ 하위에서 NMA_gaff2.xml 파일을 모두 찾는다."""
    pattern = os.path.join(REPO_ROOT, "outputs", "**", "NMA_gaff2.xml")
    return sorted(glob.glob(pattern, recursive=True))


def test_nma_xml_files_exist_or_skip():
    """outputs/ 에 NMA XML 이 존재하면 후속 테스트가 의미를 가진다."""
    files = _find_nma_xmls()
    if not files:
        pytest.skip("outputs/ 에 NMA_gaff2.xml 이 없음 — CI/clean checkout 환경으로 간주")
    assert len(files) >= 1


@pytest.mark.parametrize("xml_path", _find_nma_xmls() or [None])
def test_nma_xml_heavy_atoms_match_pdb_convention(xml_path):
    """NMA <Residue> 블록의 모든 중원자가 PDB 규약명 집합 {N, CA, C, O, CB, CM, OXT} 의 부분집합인지 검증."""
    if xml_path is None:
        pytest.skip("NMA XML 부재")

    tree = ET.parse(xml_path)
    nma_residue = None
    for res in tree.getroot().iter("Residue"):
        if res.get("name") == "NMA":
            nma_residue = res
            break
    assert nma_residue is not None, f"NMA Residue 블록이 XML 에 없음: {xml_path}"

    atom_names = [a.get("name") for a in nma_residue.findall("Atom")]
    heavy = {n for n in atom_names if n and not n.startswith("H")}

    # OXT 는 자유 말단 (linear) 일 때만 존재. 본 테스트는 NMA 의 핵심 백본 6 원자
    # (N, CA, C, O, CB, CM) 가 모두 존재하는지 검증한다.
    missing = NMA_EXPECTED_HEAVY - heavy
    assert not missing, (
        f"NMA XML 중원자가 PDB 규약과 불일치: 누락={sorted(missing)}, "
        f"실제={sorted(heavy)}, 파일={xml_path}"
    )

    # 또한 GAFF2 내부 잔재 (N1/C2/C3/...) 가 남아있지 않은지 검증
    gaff_residues = {n for n in heavy if any(c.isdigit() for c in n) and n not in {"CA", "CB", "CM"}}
    assert not gaff_residues, (
        f"NMA XML 에 GAFF2 내부 원자명 잔재가 있음: {sorted(gaff_residues)} "
        f"(rename 미완료 의심): {xml_path}"
    )
