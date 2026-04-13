#!/usr/bin/env python
"""
ncaa_registry.py
----------------
UPDD 파이프라인 전체에서 공유되는 ncAA (비천연 아미노산) 단일 진실 공급원(Single Source of Truth).
모든 하위 스크립트(Mutation, Parameterization, MD)는 이 레지스트리의 정의(Contract)를 따릅니다.

v2
- [VALIDATION] 임포트 시점에 모든 code, resname, alias의 중복(Collision)을 검사하는 validate_registry() 도입
- [SERIALIZATION] 하위 스크립트의 manifest 생성을 돕기 위한 to_manifest_dict() 메서드 추가
- [CHEM_META] mutation 단계의 하드코딩을 없애기 위한 keep_atoms 속성 추가
- [CHEM_META] QM/MM 및 전하 상태(Protonation)에 대한 명시적 설명(charge_model_note) 추가
- [PERFORMANCE] O(1) 탐색 및 충돌 방지를 위한 build_alias_index() 캐싱 구조 도입
"""

from dataclasses import dataclass, asdict
from typing import Tuple, Optional, Set

@dataclass(frozen=True)
class NCAADef:
    label: str               # 사용자 표시명 (예: N-Methylalanine)
    code: str                # 내부 Canonical Code (예: NMA)
    pdb_resname: str         # PDB 3-letter 잔기명 (예: NMA)
    xml_resname: str         # XML 템플릿 잔기명 (예: NMA)
    element: str             # 기본 중심 원소 (MD/일반용, 예: C, Si)
    formal_charge: int       # Parameterization 계산용 정수 전하 (제공된 SMILES 양성자화 상태 기준)
    charge_model_note: str   # 전하 상태 및 양성자화(Protonation) 기준에 대한 화학적 메모
    smiles_free: str         # 양 말단 자유산 SMILES (RESP 전하 계산용)
    smiles_capped: str       # 양 말단 Capped SMILES (파라미터화 보조용)
    mutation_type: str       # Mutate 스크립트용 보존/치환 규칙 키
    keep_atoms: Tuple[str, ...] # Mutate 시 절제하지 않고 남길 기본 백본 원자 목록
    # [Convention] aliases는 모두 소문자(lowercase)로 등록한다.
    # _build_and_validate_registry()가 모든 alias 비교를 .lower() 정규화를 거쳐
    # 수행하므로, 대소문자 혼용은 충돌 검출 및 O(1) 탐색에서 의도치 않은 누락을
    #유발할 수 있다. 새 NCAADef 항목을 추가할 때도 이 규약을 준수한다.
    aliases: Tuple[str, ...] # 사용자 CLI 입력을 매핑하기 위한 약어 및 별칭 모음 (lowercase 권장)

    def to_manifest_dict(self) -> dict:
        """하위 스크립트가 JSON Manifest에 기록하기 쉽도록 직렬화된 딕셔너리를 반환합니다."""
        return asdict(self)

# ==========================================
# 공통 상수 (기본 펩타이드 뼈대)
# ==========================================
STANDARD_BACKBONE = ("N", "CA", "C", "O", "CB")
N_METHYL_BACKBONE = ("N", "CA", "C", "O", "CB", "CM") # N-Methyl 탄소 포함

# ==========================================
# 통합 레지스트리 데이터베이스
# ==========================================
NCAA_REGISTRY_DATA = [
    NCAADef(
        label="N-Methylalanine",
        code="NMA",
        pdb_resname="NMA",
        xml_resname="NMA",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone",
        smiles_free="CN[C@@H](C)C(=O)O",
        smiles_capped="CC(=O)N(C)[C@@H](C)C(=O)NC",
        mutation_type="N-M",
        keep_atoms=N_METHYL_BACKBONE,
        aliases=("nma", "n-me-ala", "n-methylalanine", "NMe-Ala"),
    ),
    NCAADef(
        label="Trimethylsilylalanine",
        code="TMS",
        pdb_resname="TMS",
        xml_resname="TMS",
        element="Si",
        formal_charge=0,
        charge_model_note="Neutral backbone",
        smiles_free="N[C@@H](C[Si](C)(C)C)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](C[Si](C)(C)C)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("tms", "tmsala", "tms-ala", "trimethylsilylalanine"),
    ),
    NCAADef(
        label="Silaproline",
        code="SIP",
        pdb_resname="SIP",
        xml_resname="SIP",
        element="Si",
        formal_charge=0,
        charge_model_note="Neutral cyclic backbone",
        smiles_free="O=C(O)C1C[Si](C)(C)CCN1",
        smiles_capped="CC(=O)N1CC[Si](C)(C)C1C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("sip", "silaproline", "sila-pro"),
    ),
    NCAADef(
        label="Dimethylphenylsilylalanine",
        code="DPS",
        pdb_resname="DPS",
        xml_resname="DPS",
        element="Si",
        formal_charge=0,
        charge_model_note="Neutral backbone",
        smiles_free="N[C@@H](C[Si](C)(C)c1ccccc1)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](C[Si](C)(C)c1ccccc1)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dps", "dimethylphenylsilylalanine"),
    ),
    NCAADef(
        label="4-Fluorophenylalanine",
        code="PFF",
        pdb_resname="PFF",
        xml_resname="PFF",
        element="F",
        formal_charge=0,
        charge_model_note="Neutral backbone",
        smiles_free="N[C@@H](Cc1ccc(F)cc1)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](Cc1ccc(F)cc1)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("pff", "4-fluorophenylalanine", "4-f-phe"),
    ),
    NCAADef(
        label="4-Chlorophenylalanine",
        code="CLP",
        pdb_resname="CLP",
        xml_resname="CLP",
        element="Cl",
        formal_charge=0,
        charge_model_note="Neutral backbone",
        smiles_free="N[C@@H](Cc1ccc(Cl)cc1)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](Cc1ccc(Cl)cc1)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("clp", "4-chlorophenylalanine", "4-cl-phe"),
    ),
    NCAADef(
        label="4-Bromophenylalanine",
        code="BRP",
        pdb_resname="BRP",
        xml_resname="BRP",
        element="Br",
        formal_charge=0,
        charge_model_note="Neutral backbone",
        smiles_free="N[C@@H](Cc1ccc(Br)cc1)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](Cc1ccc(Br)cc1)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("brp", "4-bromophenylalanine", "4-br-phe"),
    ),
    NCAADef(
        label="Aminoisobutyric acid (Aib)",
        code="AIB",
        pdb_resname="AIB",
        xml_resname="AIB",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone",
        smiles_free="NC(C)(C)C(=O)O",
        smiles_capped="CC(=O)NC(C)(C)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("aib", "aminoisobutyric acid"),
    ),
    NCAADef(
        label="D-Phenylalanine",
        code="DPN",
        pdb_resname="DPN",
        xml_resname="DPN",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone (D-isomer)",
        smiles_free="N[C@H](Cc1ccccc1)C(=O)O",
        smiles_capped="CC(=O)N[C@H](Cc1ccccc1)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dpn", "d-phenylalanine", "d-phe"),
    ),
    NCAADef(
        label="D-Alanine",
        code="DAL",
        pdb_resname="DAL",
        xml_resname="DAL",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone (D-isomer)",
        smiles_free="N[C@H](C)C(=O)O",
        smiles_capped="CC(=O)N[C@H](C)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dal", "d-alanine", "d-ala"),
    ),
    NCAADef(
        label="Propargylglycine",
        code="PRA",
        pdb_resname="PRA",
        xml_resname="PRA",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone (Alkyne)",
        smiles_free="N[C@@H](CC#C)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](CC#C)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("pra", "ppg", "propargylglycine"),
    ),
    NCAADef(
        label="Azidohomoalanine",
        code="AHA",
        pdb_resname="AHA",
        xml_resname="AHA",
        element="N",
        formal_charge=0,
        charge_model_note="Neutral backbone; pendant azide (-N=N+=N-) carries internal +1/-1 formal charges with net 0",
        smiles_free="N[C@@H](CCN=[N+]=[N-])C(=O)O",
        smiles_capped="CC(=O)N[C@@H](CCN=[N+]=[N-])C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("aha", "azidohomoalanine"),
    ),
    NCAADef(
        label="Phosphoserine",
        code="SEP",
        pdb_resname="SEP",
        xml_resname="SEP",
        element="P",
        formal_charge=0,
        charge_model_note="Neutral state P(=O)(O)O provided via SMILES. (Note: Physiological pH usually requires -2 charge)",
        smiles_free="N[C@@H](COP(=O)(O)O)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](COP(=O)(O)O)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("sep", "phosphoserine", "p-ser"),
    ),
]

# ==========================================
# 무결성 검증 및 인덱싱 (Init 시점 실행)
# ==========================================
_ALIAS_INDEX = {}
_NCAA_REGISTRY_MAP = {}

def _build_and_validate_registry():
    """
    레지스트리 데이터의 무결성(이름, 코드, PDB명, Alias 간의 충돌 여부)을 검사하고,
    O(1) 속도로 탐색할 수 있는 글로벌 소문자 인덱스를 생성합니다.
    """
    global _ALIAS_INDEX, _NCAA_REGISTRY_MAP
    
    seen_codes = set()
    seen_pdb_names = set()
    seen_aliases = set()

    # [Diagnostics] 충돌 발생 시 출처(이전 등록 항목의 code)를 함께 보고하기 위한 보조 인덱스
    code_owner: dict = {}
    pdb_owner: dict = {}
    alias_owner: dict = {}

    for spec in NCAA_REGISTRY_DATA:
        # 1. Code 중복 검사
        code_upper = spec.code.upper()
        if code_upper in seen_codes:
            raise ValueError(
                f"[Registry Error] 중복된 ncAA Code 발견: '{code_upper}' "
                f"(신규: '{spec.code}', 기존 등록자: '{code_owner.get(code_upper)}')"
            )
        seen_codes.add(code_upper)
        code_owner[code_upper] = spec.code
        _NCAA_REGISTRY_MAP[code_upper] = spec

        # 2. PDB Resname 중복 검사
        pdb_upper = spec.pdb_resname.upper()
        if pdb_upper in seen_pdb_names:
            raise ValueError(
                f"[Registry Error] 중복된 PDB Resname 발견: '{pdb_upper}' "
                f"(신규 code: '{spec.code}', 기존 등록자 code: '{pdb_owner.get(pdb_upper)}')"
            )
        seen_pdb_names.add(pdb_upper)
        pdb_owner[pdb_upper] = spec.code

        # 3. Alias 충돌 검사 및 인덱싱
        names_to_index = {
            spec.label.lower(),
            spec.code.lower(),
            spec.pdb_resname.lower(),
            spec.xml_resname.lower()
        }
        names_to_index.update([a.lower() for a in spec.aliases])

        for name in names_to_index:
            if name in seen_aliases:
                raise ValueError(
                    f"[Registry Error] 중복된 식별자/Alias 발견: '{name}' "
                    f"(신규 code: '{spec.code}', 기존 등록자 code: '{alias_owner.get(name)}')"
                )
            seen_aliases.add(name)
            alias_owner[name] = spec.code
            _ALIAS_INDEX[name] = spec

# 모듈 로드 시 즉시 검증 및 빌드 (Fail-Fast)
_build_and_validate_registry()


# ==========================================
# SMILES 유효성 검증 (CLAUDE.md v3 8-1)
# ==========================================
# [Chemistry] RDKit이 사용 가능한 환경에서는 모든 NCAA 항목의 smiles_free 와
# smiles_capped 필드를 RDKit Chem.MolFromSmiles 로 파싱하여 유효성을 검증한다.
# 파싱이 실패하는 SMILES 가 발견되면 즉시 ValueError 를 발생시켜 (Fail-Fast)
# 손상된 데이터가 하위 단계 (parameterize_ncaa, run_restrained_md) 로 전파되는
# 것을 차단한다. RDKit 미설치 환경에서는 검증을 건너뛰고 정상 진입을 허용한다.
def _validate_smiles_if_rdkit_available() -> None:
    try:
        from rdkit import Chem
        from rdkit import RDLogger
    except ImportError:
        return  # RDKit 미설치: 본 검증은 RDKit 의존 환경에서만 활성화

    # RDKit 의 SMILES 파서 경고를 일시적으로 억제 (검증 자체는 그대로 수행)
    RDLogger.DisableLog("rdApp.warning")
    try:
        for spec in NCAA_REGISTRY_DATA:
            for field_name in ("smiles_free", "smiles_capped"):
                smi = getattr(spec, field_name)
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    raise ValueError(
                        f"[Registry] SMILES 파싱 실패 — code='{spec.code}' "
                        f"field='{field_name}' (RDKit MolFromSmiles 가 None 을 반환)."
                    )
    finally:
        RDLogger.EnableLog("rdApp.warning")


# 임포트 시 SMILES 유효성 검증 (RDKit 사용 가능 시).
_validate_smiles_if_rdkit_available()


# ==========================================
# 외부 공개 API
# ==========================================
def get_all_ncaa_labels() -> list:
    """UI 메뉴 표시를 위한 모든 ncAA 라벨 리스트를 반환합니다."""
    return [spec.label for spec in NCAA_REGISTRY_DATA]

def resolve_ncaa_definition(selected_name: str) -> Optional[NCAADef]:
    """
    사용자 입력 문자열을 받아 인덱스에서 정확한 NCAADef 객체를 찾아 반환합니다.
    O(1) 탐색을 보장하며, 알 수 없는 이름일 경우 명확한 에러를 발생시킵니다.
    """
    if not selected_name or selected_name.lower().strip() in ("none", "none (wild-type md)", "wild-type"):
        return None
        
    query = selected_name.lower().strip()
    spec = _ALIAS_INDEX.get(query)
    
    if spec is not None:
        return spec
            
    supported = get_all_ncaa_labels()
    raise ValueError(f"[!] 알 수 없는 ncAA 이름입니다: '{selected_name}'.\n지원되는 이름 목록: {supported}")
