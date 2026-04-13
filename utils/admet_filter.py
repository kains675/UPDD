#!/usr/bin/env python
"""
admet_filter.py
---------------
RDKit을 활용하여 화합물의 물리화학적 특성(Descriptors)을 계산하고,
세포막(Cell Membrane) 및 혈뇌장벽(BBB) 투과성을 초고속으로 평가하는 모듈입니다.
"""

from typing import Tuple, Union, Any
from rdkit import Chem
from rdkit.Chem import Descriptors

# ==========================================
# 모듈 상수 — 컷오프 분리
# ==========================================
# [Refactor] Lipinski Rule of Five 및 BBB(혈뇌장벽) 컷오프를 모듈 상수로 분리하여
# 값의 출처(문헌)를 명확히 하고 운영자가 한 곳에서 조정할 수 있게 한다.
LIPINSKI_MW_MAX     = 1200.0  # Da, bRo5 펩타이드 기준 (Doak et al. JCIM 2016)
                              # FDA 승인 환형 펩타이드: Cyclosporine A 1202 Da, Vancomycin 1449 Da
                              # 기존 Lipinski 500 Da 는 경구 소분자 기준으로, ncAA 환형 펩타이드
                              # (MW 800-2000 Da) 에는 부적합하여 상향 조정.
LIPINSKI_LOGP_MAX   = 5.0     # -, Lipinski 1997
LIPINSKI_HBD_MAX    = 10      # -, bRo5 기준 (Doak et al. JCIM 2016)
                              # 펩타이드는 다수의 NH 를 가지므로 기존 5 에서 상향.
LIPINSKI_HBA_MAX    = 10      # -, Lipinski 1997
LIPINSKI_VIOLATIONS = 1       # 허용 위반 수 (Lipinski는 1개 위반까지 허용)

BBB_MW_MAX    = 400.0  # Da, CNS MPO simplified
BBB_LOGP_MIN  = 2.0    # -, CNS MPO simplified
BBB_LOGP_MAX  = 5.0    # -
BBB_HBD_MAX   = 3      # -
BBB_TPSA_MAX  = 90.0   # Å^2


def check_admet_rules(smiles_or_mol: Union[str, Any], mode: str = "none") -> Tuple[bool, str]:
    """
    분자의 물성을 계산하고 투과성을 평가하는 톨게이트 함수
    - mode: "none" (스킵), "cell" (세포막/리핀스키), "bbb" (혈뇌장벽)
    """
    if mode == "none":
        return True, "ADMET 필터 스킵됨"

    # 1. 분자 객체 준비 (SMILES 텍스트를 RDKit 분자로 변환)
    mol = Chem.MolFromSmiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else smiles_or_mol
    if not mol:
        return False, "분자 구조 인식 실패"

    # 2. 핵심 물성(Descriptors) 계산
    mw = Descriptors.MolWt(mol)              # 분자량
    logp = Descriptors.MolLogP(mol)          # 지용성
    hbd = Descriptors.NumHDonors(mol)        # 수소 결합 주개
    hba = Descriptors.NumHAcceptors(mol)     # 수소 결합 받개
    tpsa = Descriptors.TPSA(mol)             # 극성 표면적

    # 3. 모드별 컷오프(Cut-off) 심사
    if mode == "cell":
        # [리핀스키 5법칙] - 일반 세포막 투과 조건 (모듈 상수 LIPINSKI_*)
        violations = 0
        if mw   > LIPINSKI_MW_MAX:   violations += 1
        if logp > LIPINSKI_LOGP_MAX: violations += 1
        if hbd  > LIPINSKI_HBD_MAX:  violations += 1
        if hba  > LIPINSKI_HBA_MAX:  violations += 1

        if violations <= LIPINSKI_VIOLATIONS:
            return True, f"Pass (세포막): MW={mw:.1f}, LogP={logp:.2f}, 위반={violations}개"
        else:
            return False, f"Fail (세포막): 위반 {violations}개 (MW={mw:.1f}, LogP={logp:.2f})"

    elif mode == "bbb":
        # [혈뇌장벽(BBB) 투과 조건] - CNS MPO 기반 간소화 룰 (모듈 상수 BBB_*)
        if mw > BBB_MW_MAX or logp > BBB_LOGP_MAX or logp < BBB_LOGP_MIN or hbd > BBB_HBD_MAX or tpsa > BBB_TPSA_MAX:
            return False, f"Fail (BBB): 뇌 통과 불가 (MW={mw:.1f}, TPSA={tpsa:.1f}, LogP={logp:.2f})"
        else:
            return True, f"Pass (BBB): 뇌 통과 가능! (MW={mw:.1f}, TPSA={tpsa:.1f}, LogP={logp:.2f})"

    return True, f"알 수 없는 모드: {mode}"


# ==========================================
# 단독 테스트용 블록
# ==========================================
if __name__ == "__main__":
    print("=== ADMET 필터 테스트 ===")
    
    # 테스트 1: 아스피린 (작고 훌륭한 약물)
    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    print("\n[테스트 분자: 아스피린]")
    is_pass_cell, msg_cell = check_admet_rules(aspirin_smiles, mode="cell")
    is_pass_bbb, msg_bbb = check_admet_rules(aspirin_smiles, mode="bbb")
    print(f"세포막 테스트: {msg_cell}")
    print(f"BBB 테스트  : {msg_bbb}")

    # 테스트 2: 가상의 거대한 펩타이드 모방체 (무겁고 수소결합이 많음)
    giant_smiles = "CC(C)CC(NC(=O)C(CC1=CC=CC=C1)NC(=O)C(CC(=O)O)NC(=O)C(CO)NC(=O)C(C)NC(=O)C(CC(=O)N)NC(=O)C(CC1=CN=CN1)N)C(=O)O"
    print("\n[테스트 분자: 무거운 펩타이드 사슬]")
    is_pass_cell2, msg_cell2 = check_admet_rules(giant_smiles, mode="cell")
    is_pass_bbb2, msg_bbb2 = check_admet_rules(giant_smiles, mode="bbb")
    print(f"세포막 테스트: {msg_cell2}")
    print(f"BBB 테스트  : {msg_bbb2}")
