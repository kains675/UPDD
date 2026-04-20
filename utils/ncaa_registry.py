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

from dataclasses import dataclass, asdict, field
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
    # [v0.6.6] Strategy A (hybrid amber14SB + GAFF2) sidechain preservation
    # 핵심: parent_residue 가 지정되면 source PDB 의 parent sidechain 전체를
    # 보존 (backbone-only strip 방지), extension_atoms 는 기하학적으로 부여.
    # parameterize_ncaa 의 _rename_xml_to_pdb_names 는 이 값을 사용하여
    # amber-compatible 네이밍으로 rename 한다.
    # Ref: Capece 2012 (doi:10.1021/jp2082825), Khoury 2014 (doi:10.1021/sb400168u).
    parent_residue: Optional[str] = None  # 예: "TRP" for MTR, "SER" for SEP
    # extension_atoms: (atom_name, element, attach_atom, bond_length_A) 튜플 목록
    # 예: MTR 의 CM 은 ("CM", "C", "NE1", 1.466) — NE1 에 1.466Å 로 부착된 sp3 C
    extension_atoms: Tuple[Tuple[str, str, str, float], ...] = ()

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
    NCAADef(
        label="1-Methyltryptophan",
        code="MTR",
        pdb_resname="MTR",
        xml_resname="MTR",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone. Indole N1 is methylated (N-Me), "
                          "removing the indole N-H H-bond donor present in Trp. "
                          "Side chain remains neutral.",
        smiles_free="Cn1cc(C[C@H](N)C(=O)O)c2ccccc21",
        smiles_capped="CC(=O)N[C@@H](Cc1cn(C)c2ccccc12)C(=O)NC",
        mutation_type="STANDARD",
        # [v0.6.6] Full Trp sidechain preserved from source residue per Strategy A
        # (Capece 2012 doi:10.1021/jp2082825 — hybrid amber14SB Trp scaffold + GAFF2 methyl).
        # CM (N1-methyl) appended via extension_atoms geometry at NE1.
        keep_atoms=("N", "CA", "C", "O", "CB", "CG", "CD1", "NE1", "CE2",
                    "CD2", "CE3", "CZ2", "CZ3", "CH2", "CM"),
        parent_residue="TRP",
        extension_atoms=(("CM", "C", "NE1", 1.466),),  # NE1-CM amber N-alkyl bond
        aliases=("mtr", "1-me-trp", "1-methyltryptophan", "1-methyl-tryptophan",
                 "n1-methyl-trp", "1-me-tryptophan"),
    ),

    # ==========================================
    # Tier-1 expansion (verdict_ncaa_registry_expansion_20260419)
    # Category A — N-Methyl Backbone (7 entries)
    # ==========================================
    NCAADef(
        label="N-Methyl-Leucine",
        code="NML",
        pdb_resname="MLE",
        xml_resname="MLE",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral N-methyl backbone (Nα-Me). Keep CM atom for Nα-methyl; side chain unchanged.",
        smiles_free="CN[C@@H](CC(C)C)C(=O)O",
        smiles_capped="CC(=O)N(C)[C@@H](CC(C)C)C(=O)NC",
        mutation_type="N-M",
        keep_atoms=N_METHYL_BACKBONE,
        aliases=("nml", "mle", "n-me-leu", "nme-leu", "n-methyl-leu", "n-methyl-leucine"),
    ),
    NCAADef(
        label="N-Methyl-Valine",
        code="NMV",
        pdb_resname="MVA",
        xml_resname="MVA",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral N-methyl backbone.",
        smiles_free="CN[C@@H](C(C)C)C(=O)O",
        smiles_capped="CC(=O)N(C)[C@@H](C(C)C)C(=O)NC",
        mutation_type="N-M",
        keep_atoms=N_METHYL_BACKBONE,
        aliases=("nmv", "mva", "n-me-val", "nme-val", "n-methyl-val", "n-methyl-valine"),
    ),
    NCAADef(
        label="N-Methyl-Phenylalanine",
        code="MEA",
        pdb_resname="MEA",
        xml_resname="MEA",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral N-methyl backbone.",
        smiles_free="CN[C@@H](Cc1ccccc1)C(=O)O",
        smiles_capped="CC(=O)N(C)[C@@H](Cc1ccccc1)C(=O)NC",
        mutation_type="N-M",
        keep_atoms=N_METHYL_BACKBONE,
        aliases=("mea", "nmf", "n-me-phe", "nme-phe", "n-methyl-phe", "n-methyl-phenylalanine"),
    ),
    NCAADef(
        label="N-Methyl-Lysine (Nα-methyl)",
        code="NMK",
        pdb_resname="NMK",
        xml_resname="NMK",
        element="N",
        formal_charge=1,
        charge_model_note="Nα-methyl backbone + free Nε (side-chain pKa ~10.5). pH 7: protonated ammonium → +1. "
                          "Distinct from Nε-methyl-Lys (MLZ).",
        smiles_free="CN[C@@H](CCCCN)C(=O)O",
        smiles_capped="CC(=O)N(C)[C@@H](CCCC[NH3+])C(=O)NC",
        mutation_type="N-M",
        keep_atoms=N_METHYL_BACKBONE,
        aliases=("nmk", "n-me-lys", "nme-lys", "n-alpha-me-lys", "n-methyl-lysine"),
    ),
    NCAADef(
        label="N-Methyl-Arginine (Nα-methyl)",
        code="NMR",
        pdb_resname="NMR",
        xml_resname="NMR",
        element="N",
        formal_charge=1,
        charge_model_note="Nα-methyl backbone + free guanidinium (pKa ~12.5). pH 7: +1. Distinct from Nω-Me-Arg (AGM).",
        smiles_free="CN[C@@H](CCCNC(=N)N)C(=O)O",
        smiles_capped="CC(=O)N(C)[C@@H](CCCNC(=[NH2+])N)C(=O)NC",
        mutation_type="N-M",
        keep_atoms=N_METHYL_BACKBONE,
        aliases=("nmr-arg", "n-me-arg", "nme-arg", "n-alpha-me-arg", "n-methyl-arginine"),
    ),
    NCAADef(
        label="Sarcosine (N-methyl-Glycine)",
        code="SAR",
        pdb_resname="SAR",
        xml_resname="SAR",
        element="C",
        formal_charge=0,
        charge_model_note="N-methyl-Gly. Neutral. Gly base has no CB — keep_atoms omits CB.",
        smiles_free="CNCC(=O)O",
        smiles_capped="CC(=O)N(C)CC(=O)NC",
        mutation_type="N-M",
        keep_atoms=("N", "CA", "C", "O", "CM"),
        aliases=("sar", "sarcosine", "n-me-gly", "nme-gly", "n-methyl-glycine"),
    ),
    NCAADef(
        label="N-Methyl-Glutamine",
        code="NMQ",
        pdb_resname="MEQ",
        xml_resname="MEQ",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral N-methyl backbone; neutral amide side chain.",
        smiles_free="CN[C@@H](CCC(=O)N)C(=O)O",
        smiles_capped="CC(=O)N(C)[C@@H](CCC(=O)N)C(=O)NC",
        mutation_type="N-M",
        keep_atoms=N_METHYL_BACKBONE,
        aliases=("nmq", "meq", "n-me-gln", "nme-gln", "n-methyl-gln", "n-methyl-glutamine"),
    ),

    # ==========================================
    # Category B — D-Amino Acids (5 entries)
    # ==========================================
    NCAADef(
        label="D-Valine",
        code="DVA",
        pdb_resname="DVA",
        xml_resname="DVA",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone (D-isomer).",
        smiles_free="N[C@H](C(C)C)C(=O)O",
        smiles_capped="CC(=O)N[C@H](C(C)C)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dva", "d-val", "d-valine"),
    ),
    NCAADef(
        label="D-Leucine",
        code="DLE",
        pdb_resname="DLE",
        xml_resname="DLE",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone (D-isomer).",
        smiles_free="N[C@H](CC(C)C)C(=O)O",
        smiles_capped="CC(=O)N[C@H](CC(C)C)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dle", "d-leu", "d-leucine"),
    ),
    NCAADef(
        label="D-Tryptophan",
        code="DTR",
        pdb_resname="DTR",
        xml_resname="DTR",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral backbone (D-isomer); indole N-H retained.",
        smiles_free="N[C@H](Cc1c[nH]c2ccccc12)C(=O)O",
        smiles_capped="CC(=O)N[C@H](Cc1c[nH]c2ccccc12)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dtr", "d-trp", "d-tryptophan"),
    ),
    NCAADef(
        label="D-Arginine",
        code="DAR",
        pdb_resname="DAR",
        xml_resname="DAR",
        element="N",
        formal_charge=1,
        charge_model_note="Guanidinium +1 (pKa ~12.5). D-isomer. R-16 guard: signature matches ARG HH11-HH22 pattern.",
        smiles_free="N[C@H](CCCNC(=N)N)C(=O)O",
        smiles_capped="CC(=O)N[C@H](CCCNC(=[NH2+])N)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dar", "d-arg", "d-arginine"),
    ),
    NCAADef(
        label="D-Proline",
        code="DPR",
        pdb_resname="DPR",
        xml_resname="DPR",
        element="C",
        formal_charge=0,
        charge_model_note="Neutral cyclic backbone (D-Pro). β-turn inducer.",
        smiles_free="OC(=O)[C@H]1CCCN1",
        smiles_capped="CC(=O)N1CCC[C@H]1C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dpr", "d-pro", "d-proline"),
    ),

    # ==========================================
    # Category C — Cyclization-Friendly Side Chains (7 entries)
    # ==========================================
    NCAADef(
        label="Ornithine",
        code="ORN",
        pdb_resname="ORN",
        xml_resname="ORN",
        element="N",
        formal_charge=1,
        charge_model_note="Side-chain δ-amine (pKa ~10.8) → pH 7 ammonium +1. Shorter than Lys by one CH2. "
                          "R-16 signature: terminal NH3+ pattern analogous to LYS HZ.",
        smiles_free="N[C@@H](CCCN)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](CCC[NH3+])C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("orn", "ornithine", "l-ornithine"),
    ),
    NCAADef(
        label="2,4-Diaminobutyrate",
        code="DAB",
        pdb_resname="DAB",
        xml_resname="DAB",
        element="N",
        formal_charge=1,
        charge_model_note="γ-amine side chain (pKa ~10) → pH 7 ammonium +1. Two CH2 shorter than Lys.",
        smiles_free="N[C@@H](CCN)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](CC[NH3+])C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("dab", "2,4-diaminobutyrate", "dbu"),
    ),
    NCAADef(
        label="Homoarginine",
        code="HAR",
        pdb_resname="HAR",
        xml_resname="HAR",
        element="N",
        formal_charge=1,
        charge_model_note="Arg + one extra CH2 on side chain. Guanidinium +1 at pH 7.",
        smiles_free="N[C@@H](CCCCNC(=N)N)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](CCCCNC(=[NH2+])N)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("har", "homoarginine", "h-arg"),
    ),
    NCAADef(
        label="Norleucine",
        code="NLE",
        pdb_resname="NLE",
        xml_resname="NLE",
        element="C",
        formal_charge=0,
        charge_model_note="Unbranched n-butyl side chain. Leu isostere without δ-branch. Met-oxidation-free.",
        smiles_free="N[C@@H](CCCCC)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](CCCCC)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("nle", "norleucine", "nor-leu"),
    ),
    NCAADef(
        label="4-Hydroxyproline (trans)",
        code="HYP",
        pdb_resname="HYP",
        xml_resname="HYP",
        element="C",
        formal_charge=0,
        charge_model_note="Proline + 4-OH (trans). Neutral. Collagen-mimic rigid turn.",
        smiles_free="O[C@@H]1C[C@H](C(=O)O)NC1",
        smiles_capped="CC(=O)N1C[C@H](O)C[C@@H]1C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("hyp", "4-hydroxyproline", "trans-hyp"),
    ),
    NCAADef(
        label="Tic (Tetrahydroisoquinoline-3-COOH)",
        code="TIC",
        pdb_resname="TIC",
        xml_resname="TIC",
        element="C",
        formal_charge=0,
        charge_model_note="Fused aromatic-aliphatic rigid bicyclic Phe analog. Neutral.",
        smiles_free="OC(=O)[C@@H]1Cc2ccccc2CN1",
        smiles_capped="CC(=O)N1Cc2ccccc2C[C@@H]1C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("tic", "tic-3", "tetrahydroisoquinoline-3-carboxylate"),
    ),
    NCAADef(
        label="Cyclohexylalanine",
        code="CHA",
        pdb_resname="CHA",
        xml_resname="CHA",
        element="C",
        formal_charge=0,
        charge_model_note="Fully-saturated Phe isostere (phenyl → cyclohexyl). Hydrophobic, achiral side-chain.",
        smiles_free="N[C@@H](CC1CCCCC1)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](CC1CCCCC1)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("cha", "cyclohexylalanine", "chx-ala"),
    ),

    # ==========================================
    # Category E — Phospho PTMs (2 entries)
    # Note: formal_charge reflects pH 7 dominant state (-2). SEP is kept at 0
    # for backwards compatibility per existing R-15 rationale; TPO/PTR adopt
    # the explicit -2 SSOT. See verdict §5.1 for the follow-up SEP migration.
    # ==========================================
    NCAADef(
        label="Phosphothreonine",
        code="TPO",
        pdb_resname="TPO",
        xml_resname="TPO",
        element="P",
        formal_charge=-2,
        charge_model_note="Phosphate pKa1 ~2, pKa2 ~6.5 → pH 7 dominant -2. R-15 target_iso_net_charge must "
                          "include -2 per TPO. smiles_free provided protonated for RESP geometry; "
                          "formal_charge declares the total pH-7 state.",
        smiles_free="N[C@@H]([C@@H](C)OP(=O)(O)O)C(=O)O",
        smiles_capped="CC(=O)N[C@@H]([C@@H](C)OP(=O)(O)O)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("tpo", "phosphothreonine", "p-thr", "pt"),
    ),
    NCAADef(
        label="Phosphotyrosine",
        code="PTR",
        pdb_resname="PTR",
        xml_resname="PTR",
        element="P",
        formal_charge=-2,
        charge_model_note="Aryl phosphate pKa1 ~1, pKa2 ~6.0 → pH 7 -2. R-15 dual verification: "
                          "target_card computed charge must include -2 for each PTR residue.",
        smiles_free="N[C@@H](Cc1ccc(OP(=O)(O)O)cc1)C(=O)O",
        smiles_capped="CC(=O)N[C@@H](Cc1ccc(OP(=O)(O)O)cc1)C(=O)NC",
        mutation_type="STANDARD",
        keep_atoms=STANDARD_BACKBONE,
        aliases=("ptr", "phosphotyrosine", "p-tyr", "py"),
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


def get_ncaa_charge_map() -> dict:
    """Return ``{pdb_resname_upper: formal_charge}`` for every registered ncAA.

    Consumed by ``utils/charge_topology.py::compute_binder_chem_charge`` so
    R-16 guards pick up charge-carrying ncAAs (ORN/DAB/HAR/NMK/NMR/DAR/TPO/PTR
    …) automatically — registry is the single source of truth.
    """
    out = {}
    for spec in NCAA_REGISTRY_DATA:
        out[spec.pdb_resname.upper()] = int(spec.formal_charge)
        out[spec.xml_resname.upper()] = int(spec.formal_charge)
        out[spec.code.upper()] = int(spec.formal_charge)
    return out

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
