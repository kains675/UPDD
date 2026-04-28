#!/usr/bin/env python
"""
run_mmgbsa.py
-------------
OpenMM + PySCF RESP 전하를 이용한 MM-GBSA ΔG 결합 에너지 계산
- ΔG_bind = G_complex - G_receptor - G_ligand
- GB 모델: GBn2 (implicit solvent)
- ncAA Si/Br/Cl/F 원소의 GB 반경 자동 Override
사용 환경: qmmm (conda)
"""

import os
import sys
import glob
import json
import argparse
from typing import Optional, Dict
import numpy as np

from openmm import app, unit
import openmm as mm
from openmm.app import PDBFile, ForceField, Simulation
from openmm.app import HBonds  # [v39 FIX — Bug 2] GBn2 enum 더 이상 createSystem 에 전달하지 않음

# [v3 6-2] 보조인자 식별자(SSoT) 를 utils_common 에서 임포트
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils_common import KNOWN_COFACTORS, parse_pdb_atom_line  # noqa: E402

# [R-17 G6 FIX 2026-04-20] MM-GBSA implicit-solvent createSystem 은 GNP/ATP
# 등 비표준 cofactor 에 대해 OpenMM ForceField 의 residue template 을 필요로
# 한다. MD 쪽에서 쓰는 register_cofactor_templates / inject_cofactor_bonds
# 헬퍼를 재활용하여 동일 패턴으로 해결한다. 모듈 글로벌 변수는 multiprocessing
# 워커가 자신의 프로세스 공간에만 설정 → race condition 없음.
_COFACTOR_MOL2_PATHS: list = []

# [Option β, 2026-04-20 — SciVal Verdict 19 + ultrareview]
# Cyclic peptide HTC / SS bond 는 MD 에서 이미 검증된 helper 함수로 처리.
# 이전 write_temp_pdb 의 +5Å translation + 합성 OXT hack (F2, dominant ΔG
# bias 원인) 은 **제거**되었음. 올바른 순서: addHydrogens → trim termini →
# detect HTC pair (post-trim) → addBond. `residueTemplates` override 없이
# amber14 의 atom-set + external-bond-count 매칭이 자동으로 mid-chain GLU
# 템플릿을 선택함 (ultrareview 확증).
try:
    from run_restrained_md import (  # noqa: E402
        detect_head_to_tail_pair,
        detect_disulfide_pair,
        _trim_termini_for_htc,
        commit_bond_if_missing,
    )
    _CYCLIC_HELPERS_AVAILABLE = True
except ImportError:
    _CYCLIC_HELPERS_AVAILABLE = False


# ==========================================
# GB 반경 Override 테이블
# ==========================================
GB_RADII_OVERRIDE = {
    "Si": 2.10,
    "Br": 1.85,
    "Cl": 1.75,
    "F":  1.47,
    "P":  1.80,
    # [SciVal Verdict 18, 2026-04-20] mbondi2-style extension for divalent
    # metal ions (Mg/Zn/Ca/Mn/Fe/Cu/Co/Ni). GBn2 (Nguyen 2013) does not
    # define Born radii for metal cations; these are literature-backed
    # values (Li-Merz 2014 + Bondi 1964 VDW + SciVal Verdict 18 §4).
    # Pool-symmetric in ΔG_bind (metal lives only in receptor/complex
    # pools), so minor radius error cancels in rank ordering.
    "Mg": 1.18,
    "Zn": 1.09,
    "Ca": 1.37,
    "Mn": 1.13,
    "Fe": 1.08,
    "Cu": 1.00,
    "Co": 1.03,
    "Ni": 1.00,
}

GB_SCALE_OVERRIDE = {
    "Si": 0.80,
    "Br": 0.80,
    "Cl": 0.80,
    "F":  0.80,
    "P":  0.86,
    # [SciVal Verdict 18] Divalent metal ions — GB scale 0.85 (uniform).
    "Mg": 0.85, "Zn": 0.85, "Ca": 0.85, "Mn": 0.85,
    "Fe": 0.85, "Cu": 0.85, "Co": 0.85, "Ni": 0.85,
}



# ==========================================
# PDB 분리 유틸 — Complex / Receptor / Ligand
# ==========================================
# [v39 FIX 2026-04-19] implicit/gbn2.xml 은 물·카운터이온 템플릿을 포함하지
# 않으므로, MM-GBSA 에 투입하기 전에 반드시 제거해야 한다. explicit water 와
# 단순 카운터이온은 implicit 용매 모델이 연속 매질로 대체한다.
# [SciVal Verdict 17/18 ROLLBACK 2026-04-20] 이전 v39 FIX 의 "구조적 금속
# 이온 제거" 조항은 **literature 위반 (Gohlke-Case 2004, Panteva 2015,
# Kazemi 2021 모두 Mg 유지)** + F1 root cause (verdict_mmgbsa_rootcause)
# 로 확인되어 철회한다. Li-Merz 12-6-4 ion XML (amber14/tip3p_HFE_multivalent.xml)
# 이 OpenMM 기본 탑재 → structural metal 유지, GB radius 는
# GB_RADII_OVERRIDE 의 mbondi2 extension 이 담당.
# [v52 FIX] MD 말단 캡핑 잔기. PDBFixer가 MD 안정화를 위해 추가한 것이므로
# MM-GBSA implicit solvent 계산 전에 제거한다. 이 잔기들은 실제 단백질의
# 일부가 아니며, implicit/gbn2.xml에 템플릿이 없어 에러를 유발한다.
_CAP_RESNAMES = {"NME", "ACE", "NHE"}

_MMGBSA_STRIP_RESNAMES = {
    "HOH", "WAT", "TIP", "TIP3", "SOL",       # water
    "NA", "CL", "K", "BR", "IOD",              # monovalent counterions
    # [P0-2 ultrareview 2026-04-20] Crystallographic artifacts (buffer /
    # cryoprotectant / precipitant) — KNOWN_COFACTORS 에 들어있지만 implicit
    # GBn2 FF 에 template 없음. Receptor pool 로 routing 되면 createSystem 이
    # "No template for GOL/EDO/..." 로 silent fail → snapshot drop.
    "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "ACT", "IMD",
    # Structural divalents (MG/ZN/CA/MN/FE/CU/CO/NI/FE2/FE3) DELIBERATELY
    # NOT stripped — SciVal Verdict 18 winner (Li-Merz 12-6-4 via
    # amber14/tip3p_HFE_multivalent.xml) provides their templates.
}


def split_complex(pdb_path, receptor_chain="A", binder_chain="B"):
    """
    PDB에서 Complex / Receptor / Ligand 를 분리한다.

    [v3 6-2] HETATM 잔기 중 KNOWN_COFACTORS (GDP, GTP, NAD, HEM, MG, ZN, ...)
    에 등재된 보조인자/이온은 receptor 측에 포함되어야 한다. 이는 binding
    pocket 의 핵심 화학 환경을 보존하기 위한 조치이다.

    [v39 FIX] implicit/gbn2.xml 이 물·이온 템플릿을 포함하지 않으므로
    water, counterion, 단순 금속 이온 (MG, CL, NA 등) 은 세 풀 모두에서
    제거한다. 구조적 보조인자 (GDP, ATP, HEM 등) 는 receptor 에 유지한다.

    분류 우선순위:
        0. resname ∈ _MMGBSA_STRIP_RESNAMES → 제거 (implicit FF 호환)
        1. record == TER/END           → complex 만 기록 후 continue
        2. chain == receptor_chain     → receptor (단백질 본체)
        3. resname ∈ KNOWN_COFACTORS   → receptor (보조인자, 체인 무관)
        4. chain == binder_chain       → ligand (펩타이드)
        5. chain == ' ' or rec=HETATM  → ligand 폴백 (미등록 ncAA 등)

    Returns: (complex_lines, receptor_lines, ligand_lines)
    """
    complex_lines  = []
    receptor_lines = []
    ligand_lines   = []
    n_blank_chain = 0
    n_cofactor_to_receptor = 0
    n_stripped = 0
    # [v0.6.7] Track serial → bucket mapping so CONECT records survive the
    # split. Without this, external peptide bonds for ncAA residues (e.g.,
    # VAL3-MTR4 in Cp4) are dropped and OpenMM falls back to residue-name
    # bond inference, which fails at the ncAA boundary because MTR is not
    # a standard residue.
    serial_bucket: dict = {}  # serial -> set of buckets {"complex","receptor","ligand"}
    conect_raw_lines: list = []

    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "CONECT":
                conect_raw_lines.append(line)
                continue
            if rec not in ("ATOM", "HETATM", "TER", "END"):
                continue

            if rec in ("TER", "END"):
                complex_lines.append(line)
                continue

            chain = line[21] if len(line) > 21 else " "
            resname = line[17:20].strip() if len(line) >= 20 else ""

            # [v39] implicit FF 호환: water/ion/metal 제거
            if resname in _MMGBSA_STRIP_RESNAMES:
                n_stripped += 1
                continue

            # [v52 FIX] MD 말단 캡(NME/ACE) 제거 — implicit FF 템플릿 미지원
            if resname in _CAP_RESNAMES:
                n_stripped += 1
                continue

            complex_lines.append(line)
            try:
                atom_serial = int(line[6:11])
            except (ValueError, IndexError):
                atom_serial = None
            if atom_serial is not None:
                serial_bucket.setdefault(atom_serial, set()).add("complex")

            if chain == receptor_chain:
                receptor_lines.append(line)
                if atom_serial is not None:
                    serial_bucket[atom_serial].add("receptor")
            elif rec == "HETATM" and resname in KNOWN_COFACTORS:
                receptor_lines.append(line)
                n_cofactor_to_receptor += 1
                if atom_serial is not None:
                    serial_bucket[atom_serial].add("receptor")
            elif chain == binder_chain:
                ligand_lines.append(line)
                if atom_serial is not None:
                    serial_bucket[atom_serial].add("ligand")
            elif chain == " " or rec == "HETATM":
                ligand_lines.append(line)
                if atom_serial is not None:
                    serial_bucket[atom_serial].add("ligand")
                if chain == " ":
                    n_blank_chain += 1

    # Filter CONECT lines per bucket: emit only those where ALL referenced
    # serials are in that bucket (otherwise the bond would point at a removed
    # atom and PDBFile would fail to parse).
    def _filter_conect(bucket_name: str) -> list:
        out = []
        for cl in conect_raw_lines:
            parts = cl[6:].split()
            try:
                serials = [int(p) for p in parts]
            except ValueError:
                continue
            if not serials:
                continue
            if all(bucket_name in serial_bucket.get(s, set()) for s in serials):
                out.append(cl)
        return out

    complex_conect  = _filter_conect("complex")
    receptor_conect = _filter_conect("receptor")
    ligand_conect   = _filter_conect("ligand")
    complex_lines.extend(complex_conect)
    receptor_lines.extend(receptor_conect)
    ligand_lines.extend(ligand_conect)

    if n_stripped > 0:
        print(f"  [split_complex] implicit GB 호환: water/ion/metal {n_stripped}개 원자 제거")
    if n_cofactor_to_receptor > 0:
        print(f"  [split_complex] 보조인자 {n_cofactor_to_receptor}개 원자를 receptor 풀에 포함 (KNOWN_COFACTORS).")
    if n_blank_chain > 0:
        print(f"  [!] split_complex 경고: chain ID가 공백인 원자 {n_blank_chain}개를 ligand로 분류했습니다.")
    if conect_raw_lines:
        print(f"  [split_complex] CONECT: complex={len(complex_conect)}, receptor={len(receptor_conect)}, ligand={len(ligand_conect)} (of {len(conect_raw_lines)} total)")

    return complex_lines, receptor_lines, ligand_lines


def write_temp_pdb(lines, path):
    """임시 PDB 를 작성한다.

    [Option β 2026-04-20] 이전 v44 FIX 의 cyclic HTC disruption hack
    (+5 Å x-translation + fabricated OXT) 는 F2 root cause 로 확인되어 철회.
    [v0.6.7] CONECT 제거도 철회 — parent-based ncAA (MTR, SEP, PTR, ...) 는
    비표준 잔기 경계에서 peptide-bond CONECT 가 필요 (VAL3.C ↔ MTR4.N).
    split_complex 가 이미 bucket 별 CONECT 를 정확히 분배하므로 그대로 기록.
    """
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
        if not lines or not lines[-1].startswith("END"):
            f.write("END\n")


def _apply_cyclic_bonds_mmgbsa(modeller, binder_chain: str = ""):
    """[Option β] MM-GBSA Complex / Ligand sub-system 에 cyclic peptide
    HTC 또는 SS bond 를 주입. Receptor sub-system (binder_chain 미포함)
    에서는 no-op. addHydrogens 직후, createSystem 직전에 호출.

    Args:
        modeller: OpenMM Modeller (already addHydrogens 된 상태)
        binder_chain: binder chain id ("B" 기본). 빈 문자열이면 no-op
            (receptor-only PDB 에서 skip 용).

    Returns:
        (modeller, cyclic_info_dict)
    """
    info = {"htc_formed": False, "ss_formed": False, "binder_chain": binder_chain}
    if not binder_chain or not _CYCLIC_HELPERS_AVAILABLE:
        return modeller, info

    # 1. HTC (head-to-tail) 감지 → trim 말단 cap → re-detect → addBond
    htc_pair = detect_head_to_tail_pair(modeller, binder_chain)
    if htc_pair:
        modeller = _trim_termini_for_htc(modeller, binder_chain)
        htc_pair = detect_head_to_tail_pair(modeller, binder_chain)
        if htc_pair:
            is_new = commit_bond_if_missing(modeller.topology, *htc_pair)
            if is_new:
                print(f"  [cyclic_htc] binder chain {binder_chain} HTC peptide "
                      f"bond 주입: {htc_pair[0].residue.name}{htc_pair[0].residue.id}"
                      f".{htc_pair[0].name} → {htc_pair[1].residue.name}"
                      f"{htc_pair[1].residue.id}.{htc_pair[1].name}")
                info["htc_formed"] = True

    # 2. Disulfide bridge (cyclic_ss) 감지
    try:
        ss_pair = detect_disulfide_pair(modeller, binder_chain)
    except RuntimeError:
        ss_pair = None  # ambiguous disulfide — MD 에서 error, MM-GBSA 는 skip
    if ss_pair:
        is_new = commit_bond_if_missing(modeller.topology, *ss_pair)
        if is_new:
            print(f"  [cyclic_ss] binder chain {binder_chain} SS bond 주입")
            info["ss_formed"] = True

    return modeller, info


# ==========================================
# GB 반경 Override 적용
# ==========================================
def apply_gb_radius_override(system, topology, ncaa_elem,
                             si_radius: float = 2.10,
                             custom_radii: Optional[Dict[str, float]] = None):
    """[v3 6-1] GBSAOBCForce / CustomGBForce 의 ncAA 이종원소 GB 반경·스케일을 override 한다.

    Args:
        system: OpenMM System
        topology: OpenMM Topology
        ncaa_elem: 진단/로깅용 ncAA 핵심 원소 (예: "Si"). None 또는 "none" 이면
            ncAA 미사용 상태로 간주하나, 모듈 기본 GB_RADII_OVERRIDE 는 항상 적용된다.
        si_radius: [Backwards-Compat] Si 의 GB 반경 (Å). custom_radii 가 None 일 때만
            사용되며, ncaa_elem 이 활성 상태인 경우 override_map["Si"] 에 주입된다.
        custom_radii: 임의 원소별 GB 반경 dict (예: {"Si": 2.10, "Br": 1.90}). 지정 시
            si_radius 를 우회하여 모든 원소를 한 번에 override 할 수 있다.

    Returns:
        OpenMM System (in-place 수정 후 그대로 반환)
    """
    override_map = dict(GB_RADII_OVERRIDE)
    if custom_radii:
        override_map.update(custom_radii)
    elif ncaa_elem and ncaa_elem != "none":
        override_map["Si"] = si_radius  # 하위 호환: Si 단일 인자 경로

    # [P0-1 FIX, ultrareview 2026-04-20] GBn2 (implicit/gbn2.xml) uses
    # OpenMM's CustomGBForce with per-particle slots:
    #   params[0] = charge (e)
    #   params[1] = offset_radius (Born_radius − 0.009 Å  == 0.0009 nm)
    #   params[2] = scaled_offset_radius (scale * offset_radius)
    #   params[3] = radius_index_for_Ni (lookup table index)
    # Writing our Born radius directly into params[1] was OFF-BY-0.009 Å
    # (small). More critically, writing scale into params[2] without the
    # `scale × offset_radius` product silently zeroed the atom's GB pair
    # screening. This patch handles CustomGBForce (isinstance check)
    # explicitly with the correct offset formula.
    _GBN2_OFFSET_NM = 0.0009  # Nguyen 2013 offset
    n_overridden = 0
    for force in system.getForces():
        force_type = type(force).__name__
        if "GB" not in force_type and "GBSA" not in force_type:
            continue
        is_custom_gb = isinstance(force, mm.CustomGBForce)

        for atom in topology.atoms():
            elem = atom.element.symbol if atom.element else "C"
            if elem not in override_map:
                continue
            idx = atom.index
            r_nm = override_map[elem] * 0.1  # Å → nm
            scale = GB_SCALE_OVERRIDE.get(elem, 0.80)
            if is_custom_gb:
                try:
                    params = list(force.getParticleParameters(idx))
                    if len(params) >= 3:
                        offset_r = max(r_nm - _GBN2_OFFSET_NM, 0.01)
                        params[1] = offset_r
                        params[2] = scale * offset_r
                        force.setParticleParameters(idx, params)
                        n_overridden += 1
                    else:
                        print(f"  [!] CustomGBForce unexpected param count "
                              f"{len(params)} for atom {idx} ({elem})")
                except Exception as _cgb_err:
                    print(f"  [!] CustomGBForce param set failed for atom "
                          f"{idx} ({elem}): {_cgb_err.__class__.__name__}")
            else:
                try:
                    # GBSAOBCForce signature: (charge, radius, scale)
                    charge, _r, _s = force.getParticleParameters(idx)
                    force.setParticleParameters(
                        idx, charge, r_nm * unit.nanometer, scale,
                    )
                    n_overridden += 1
                except Exception as _obc_err:
                    print(f"  [!] GBSAOBCForce param set failed for atom "
                          f"{idx} ({elem}): {_obc_err.__class__.__name__}")

    if n_overridden > 0:
        print(f"  GB 반경 Override: {n_overridden}개 원자 ({ncaa_elem} {override_map.get(ncaa_elem,'?')} Å)")
    return system


# ==========================================
# 단일 PDB → OpenMM 에너지 계산
# ==========================================
def calc_energy(pdb_path, ff, ncaa_elem, si_radius, binder_chain: str = ""):
    """
    PDB → OpenMM implicit solvent 에너지 (kcal/mol)

    Args:
        pdb_path: input PDB path
        ff: OpenMM ForceField (cofactor + metal ion + GBn2 templates registered)
        ncaa_elem: ncAA key element (e.g., "C")
        si_radius: Si GB radius Å (legacy)
        binder_chain: [Option β] Complex / Ligand sub-PDB 의 binder chain id.
            비어 있으면 cyclic bond 주입 skip (receptor-only path).
    """
    # [v55 FIX] NME/ACE 캡 잔기를 PDB 로드 전에 제거.
    # split_complex에서 이미 필터하지만, 방어적으로 여기서도 필터.
    # 원본 PDB를 in-place 수정하지 않고 임시 파일을 사용한다.
    load_path = pdb_path
    _cap_tmp = None
    try:
        with open(pdb_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        has_cap = any(
            l[17:20].strip() in _CAP_RESNAMES
            for l in lines if l.startswith(("ATOM", "HETATM")) and len(l) >= 20
        )
        if has_cap:
            filtered = [l for l in lines if not (
                l.startswith(("ATOM", "HETATM")) and len(l) >= 20
                and l[17:20].strip() in _CAP_RESNAMES
            )]
            _cap_tmp = pdb_path + ".nocap.tmp"
            with open(_cap_tmp, "w", encoding="utf-8") as f:
                f.writelines(filtered)
            load_path = _cap_tmp
            print(f"    [Cap Filter] {os.path.basename(pdb_path)}에서 NME/ACE 캡 제거 (임시 파일)")
    except Exception as _cap_err:
        print(f"    [!] Cap 필터링 실패(무시): {_cap_err}")

    try:
        pdb = PDBFile(load_path)
    except Exception as e:
        print(f"  [!] PDB 로드 실패 ({os.path.basename(pdb_path)}): {e}")
        return None
    finally:
        if _cap_tmp and os.path.exists(_cap_tmp):
            try:
                os.remove(_cap_tmp)
            except OSError:
                pass

    try:
        # [v44 FIX] Modeller.addHydrogens 로 누락된 말단 수소 원자를 보충한다.
        # cyclic HTC bond 가 형성된 펩타이드는 표준 말단기 (N-ter NH3+, C-ter COO-)
        # 가 없어서 amber14-all.xml 템플릿과 불일치한다. addHydrogens 는 누락된
        # 말단 원자를 자동으로 추가하여 이 불일치를 해소한다.
        modeller = app.Modeller(pdb.topology, pdb.positions)
        # [R-17 G6 FIX] GNP / ATP 등 비표준 cofactor 는 mol2 connectivity 가
        # 없으면 GAFFTemplateGenerator 의 graph-isomorphism match 실패 →
        # "No template found". MD 와 동일 패턴으로 topology 에 결합 주입.
        if _COFACTOR_MOL2_PATHS:
            try:
                from run_restrained_md import inject_cofactor_bonds  # type: ignore
                inject_cofactor_bonds(modeller, _COFACTOR_MOL2_PATHS)
            except ImportError:
                pass

        # [Option β 2026-04-20] Cyclic HTC / SS bond 주입을 **addHydrogens
        # 이전** 에 수행. MD snapshot 은 이미 H 원자 완비 상태로 저장되므로
        # addHydrogens 는 불필요 (더구나 cyclic 경계에서 addHydrogens 의
        # internal createSystem 이 positional N/C-term detection 으로 인해
        # template match 실패함, ultrareview 확증).
        if binder_chain:
            modeller, _cyclic_info = _apply_cyclic_bonds_mmgbsa(modeller, binder_chain)
            # SKIP addHydrogens — cyclic 경계 template match 실패 방지.
            # MD snap 에 모든 H 가 이미 있어 no-op 에 가까움.
        else:
            # Receptor sub-PDB: cyclic 무관, addHydrogens 안전하게 호출 가능
            # (혹시 snap 에 결측 H 가 있을 때 보완).
            modeller.addHydrogens(ff)

        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod  = app.NoCutoff,
            constraints      = HBonds,
            soluteDielectric = 1.0,
            solventDielectric= 78.5,
        )
    except Exception as e:
        print(f"  [!] 시스템 생성 실패: {e}")
        return None

    # GB 반경 Override
    system = apply_gb_radius_override(system, modeller.topology, ncaa_elem, si_radius)

    integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
    # [2026-04-20] Platform 선택: UPDD_MMGBSA_PLATFORM 환경변수 우선,
    # 기본값은 CUDA fallback CPU. CPU 만 쓰려면 UPDD_MMGBSA_PLATFORM=CPU
    # 로 명시. 대형 복합체 (2QKI ~10k 원자) 에서 CPU 는 실용적이지 않음.
    _pref = os.environ.get("UPDD_MMGBSA_PLATFORM", "CUDA")
    platform = None
    for name in (_pref, "CUDA", "OpenCL", "CPU"):
        try:
            platform = mm.Platform.getPlatformByName(name)
            break
        except Exception:
            continue
    if platform is None:
        platform = mm.Platform.getPlatformByName("CPU")
    sim = Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)

    # [P1 ultrareview 2026-04-20] 500 iter / tol=10 은 large protein
    # (≥10k atom) receptor 에 완전히 불충분. 2QKI C3c (9896 원자) benchmark
    # 에서 500 iter 미완료 PE = +76,600 kcal/mol → 5000 iter/tol=1 → -20,768
    # kcal/mol (정상). 문헌 (Hou 2011, Genheden 2015) 기준 tolerance 0.1-1.0
    # kJ/mol/nm 필요. 환경변수 UPDD_MMGBSA_MIN_ITER / UPDD_MMGBSA_MIN_TOL
    # 로 override 가능.
    _min_iter = int(os.environ.get("UPDD_MMGBSA_MIN_ITER", "5000"))
    _min_tol = float(os.environ.get("UPDD_MMGBSA_MIN_TOL", "1.0"))
    sim.minimizeEnergy(maxIterations=_min_iter, tolerance=_min_tol)

    state = sim.context.getState(getEnergy=True)
    e_kj  = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    e_kcal = e_kj / 4.184
    return e_kcal


# ==========================================
# [v4 S-1] Interaction Entropy (Duan et al. JCTC 2016)
# ==========================================
def calc_interaction_entropy(delta_g_frames, temperature: float = 300.0) -> dict:
    """Duan 2016 Interaction Entropy 로 −TΔS_int 를 추정한다.

    수식:
        ΔE_int^k = E_int^k − <E_int>
        −TΔS = β^(-1) · ln <exp(β · ΔE_int^k)>
    여기서 β = 1/(k_B · T) 이고, 합산은 NPT 궤적의 모든 (또는 충분한) 스냅샷에 대해 수행된다.

    구현 노트:
        - exp(β·ΔE) 가 큰 양수에서 폭발하므로, log-sum-exp 트릭으로 안정화한다.
        - 본 추정의 신뢰 구간은 N (스냅샷 수) 의 √N 에 비례하므로 N≥5 권장.
        - 입력 단위는 kcal/mol, 출력도 kcal/mol.
        - T=300 K 에서 k_B·T ≈ 0.5961 kcal/mol.

    Args:
        delta_g_frames: 각 스냅샷의 ΔG_bind (kcal/mol) 리스트
        temperature: 절대온도 (K), 기본 300

    Returns:
        dict:
            - mean_delta_g_kcal      : ΔG 평균
            - std_delta_g_kcal       : ΔG 표준편차
            - minus_TdS_kcal         : −TΔS 추정값 (kcal/mol)
            - delta_g_corrected_kcal : 평균 ΔG + (−TΔS) (Duan 보정 후 자유에너지)
            - n_frames               : 사용된 프레임 수
            - reliable               : N≥5 인 경우 True
    """
    import math as _math
    n = len(delta_g_frames)
    if n == 0:
        return {
            "mean_delta_g_kcal": None, "std_delta_g_kcal": None,
            "minus_TdS_kcal": None, "delta_g_corrected_kcal": None,
            "n_frames": 0, "reliable": False,
        }

    arr = np.asarray(delta_g_frames, dtype=float)
    mean_g = float(np.mean(arr))
    std_g  = float(np.std(arr, ddof=1)) if n > 1 else 0.0

    # k_B in kcal/(mol·K) = 1.987204259e-3
    KB = 1.987204259e-3
    beta = 1.0 / (KB * temperature)

    deviations = arr - mean_g  # ΔE_int^k
    # log-sum-exp 안정화: ln <exp(βΔE)> = ln (Σ exp(βΔE_k)) − ln N
    bx = beta * deviations
    bx_max = float(np.max(bx))
    log_mean_exp = bx_max + _math.log(float(np.sum(np.exp(bx - bx_max))) / n)
    minus_TdS = (1.0 / beta) * log_mean_exp  # kcal/mol

    return {
        "mean_delta_g_kcal":      mean_g,
        "std_delta_g_kcal":       std_g,
        "minus_TdS_kcal":         minus_TdS,
        "delta_g_corrected_kcal": mean_g + minus_TdS,
        "n_frames":               n,
        "reliable":               n >= 5,
        "temperature_K":          temperature,
    }


# ==========================================
# MM-GBSA 계산 메인
# ==========================================
def calc_mmgbsa(pdb_path, output_dir, ff, ncaa_elem, si_radius, receptor_chain="A", binder_chain="B"):
    """단일 스냅샷에 대한 ΔG_bind 계산"""

    basename = os.path.basename(pdb_path).replace(".pdb", "")
    tmp_dir  = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    print(f"\n  계산: {basename}")

    # ── 1. Complex / Receptor / Ligand 분리 ─────────────────
    complex_lines, receptor_lines, ligand_lines = split_complex(pdb_path, receptor_chain, binder_chain)

    if not receptor_lines:
        print(f"  [!] Receptor (chain {receptor_chain}) 없음 — 건너뜀")
        return None
    if not ligand_lines:
        print(f"  [!] Ligand (chain {binder_chain}) 없음 — 건너뜀")
        return None

    cplx_pdb = os.path.join(tmp_dir, f"{basename}_complex.pdb")
    recv_pdb = os.path.join(tmp_dir, f"{basename}_receptor.pdb")
    lig_pdb  = os.path.join(tmp_dir, f"{basename}_ligand.pdb")

    # [Refactor] Principle 4: 임시 PDB 파일 정리는 어떠한 경로(에너지 계산
    # 실패, JSON 저장 실패, OpenMM 예외 등)에서도 누락되지 않도록 try/finally
    # 블록으로 감싼다. 이전 구현은 정상 경로에서만 cleanup을 수행하여 실패한
    # 디자인의 임시 파일이 tmp_dir에 누적되었다.
    try:
        write_temp_pdb(complex_lines,  cplx_pdb)
        write_temp_pdb(receptor_lines, recv_pdb)
        write_temp_pdb(ligand_lines,   lig_pdb)

        # ── 2. 각 파트 에너지 계산 ──────────────────────────────
        # [Option β] binder_chain 을 Complex / Ligand 에 전달 → cyclic HTC
        # peptide bond 가 올바르게 주입됨. Receptor 는 binder 부재 → skip.
        print("  Complex  에너지 계산 중...")
        e_complex  = calc_energy(cplx_pdb, ff, ncaa_elem, si_radius,
                                 binder_chain=binder_chain)

        print("  Receptor 에너지 계산 중...")
        e_receptor = calc_energy(recv_pdb, ff, ncaa_elem="none",
                                 si_radius=si_radius, binder_chain="")

        print("  Ligand   에너지 계산 중...")
        e_ligand   = calc_energy(lig_pdb,  ff, ncaa_elem, si_radius,
                                 binder_chain=binder_chain)

        # ── 3. ΔG 계산 ──────────────────────────────────────────
        if None in (e_complex, e_receptor, e_ligand):
            print("  [!] 에너지 계산 실패 — 건너뜀")
            return None

        delta_g = e_complex - e_receptor - e_ligand

        print(f"  E_complex  : {e_complex:>12.4f} kcal/mol")
        print(f"  E_receptor : {e_receptor:>12.4f} kcal/mol")
        print(f"  E_ligand   : {e_ligand:>12.4f} kcal/mol")
        print(f"  ΔG_bind    : {delta_g:>+12.4f} kcal/mol  {'✅ 결합 유리' if delta_g < 0 else '⚠️  결합 불리'}")

        result = {
            "snapshot":     basename,
            "pdb_path":     pdb_path,
            "ncaa_element": ncaa_elem,
            "si_radius_A":  si_radius,
            "e_complex_kcal":  e_complex,
            "e_receptor_kcal": e_receptor,
            "e_ligand_kcal":   e_ligand,
            "delta_g_kcal":    delta_g,
            "favorable":       delta_g < 0,
        }

        result["protocol_variant"] = "3traj"

        # 개별 결과 저장
        out_json = os.path.join(output_dir, f"{basename}_mmgbsa.json")
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)

        return result
    finally:
        # 임시 파일 정리는 정상/예외 경로 모두에서 보장된다.
        for tmp in [cplx_pdb, recv_pdb, lig_pdb]:
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except OSError as _rm_err:
                    print(f"  [!] 임시 파일 정리 실패(무시): {tmp} ({_rm_err})")


# ==========================================
# [v0.6.8 Fix 4] 1-trajectory MM-GBSA (opt-in via UPDD_MMGBSA_PROTOCOL=1traj)
# ==========================================
# Physical rationale (Genheden & Ryde 2015 §2.1, DOI:10.1517/17460441.2015.1032250):
#   - Minimize ONLY the complex (bound state).
#   - Evaluate E_receptor, E_ligand STATICALLY at the same bound geometry
#     (subsystem-rebuild: new System with amber14 + GBn2, but positions
#      transferred from the minimized complex — no independent minimization).
#   - Bound-state strain cancels symmetrically between E_complex and
#     (E_rec_static + E_lig_static) → variance ↓ + systematic bias ↓.
#
# Empirical backing:
#   - SciVal verdict 2026-04-22: 2QKI_WT (σ₁/σ₃ = 0.111, 9.0× tighter)
#     + 3IOL_WT (σ₁/σ₃ = 0.563, 1.77× tighter) both approved (Case A).
#   - Reference standalone validator: scripts/mmgbsa_1traj_compare.py (513 LOC).
#
# Integration notes:
#   - Co-exists with 3-traj default; dispatcher reads UPDD_MMGBSA_PROTOCOL.
#   - Reuses split_complex() (CONECT bucket preservation) → ncAA boundary
#     bonds survive. 1-traj does not alter F1/F2 behavior; F1 Mg²⁺ preserved
#     via _MMGBSA_STRIP_RESNAMES (unchanged), F2 no +5 Å translate hack.
#   - P1 minimize_iter default 5000, env UPDD_MMGBSA_MIN_ITER override honored.


def split_complex_1traj(complex_pdb, receptor_chain="A", binder_chain="B"):
    """[v0.6.8 Fix 4] Partition complex atoms by chain (receptor/ligand) and
    return three tmp PDB paths operating on the SAME input coordinates.

    Unlike the 3-traj `split_complex`, this helper is intended for the 1-traj
    workflow where receptor / ligand subsystems share the bound geometry of
    the minimized complex. It reuses `split_complex` (same STRIP rules, same
    CONECT bucket filtering, same KNOWN_COFACTORS → receptor routing) to avoid
    duplicating parsing logic. The only structural change is that the CALLER
    evaluates receptor/ligand energies statically (no separate minimization).

    Args:
        complex_pdb: source snapshot PDB path.
        receptor_chain: receptor chain id (default "A").
        binder_chain: binder chain id (default "B").

    Returns:
        (complex_lines, receptor_lines, ligand_lines) — PDB lines in memory.
        Caller writes tmp PDBs via write_temp_pdb (identical to 3-traj path).
    """
    return split_complex(complex_pdb, receptor_chain, binder_chain)


def calc_energy_static(pdb_path, ff, ncaa_elem, si_radius,
                       override_positions=None, binder_chain: str = ""):
    """[v0.6.8 Fix 4] Evaluate a subsystem's potential energy WITHOUT
    minimization. Mirrors `calc_energy` structurally (same Modeller / FF /
    GB radius override / platform selection), but skips `minimizeEnergy()`.

    Used for E_rec_static, E_lig_static in the 1-traj workflow.

    Args:
        pdb_path: subsystem PDB (receptor-only or ligand-only).
        ff: OpenMM ForceField (same as 3-traj).
        ncaa_elem: ncAA key element (passed to apply_gb_radius_override).
        si_radius: legacy Si radius argument.
        override_positions: optional sequence of OpenMM Vec3 positions to
            inject after addHydrogens. If None, uses the PDB's loaded
            positions (legacy behavior — unused in 1-traj, retained for
            downstream flexibility).
        binder_chain: binder chain id for cyclic HTC / SS injection.
            Receptor-only subsystems should pass "" (no binder).

    Returns:
        Potential energy in kcal/mol, or None on failure.
    """
    # Cap filter (identical to calc_energy).
    load_path = pdb_path
    _cap_tmp = None
    try:
        with open(pdb_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        has_cap = any(
            l[17:20].strip() in _CAP_RESNAMES
            for l in lines if l.startswith(("ATOM", "HETATM")) and len(l) >= 20
        )
        if has_cap:
            filtered = [l for l in lines if not (
                l.startswith(("ATOM", "HETATM")) and len(l) >= 20
                and l[17:20].strip() in _CAP_RESNAMES
            )]
            _cap_tmp = pdb_path + ".nocap.tmp"
            with open(_cap_tmp, "w", encoding="utf-8") as f:
                f.writelines(filtered)
            load_path = _cap_tmp
    except Exception as _cap_err:
        print(f"    [!] Cap 필터링 실패(무시): {_cap_err}")

    try:
        pdb = PDBFile(load_path)
    except Exception as e:
        print(f"  [!] PDB 로드 실패 ({os.path.basename(pdb_path)}): {e}")
        return None
    finally:
        if _cap_tmp and os.path.exists(_cap_tmp):
            try:
                os.remove(_cap_tmp)
            except OSError:
                pass

    try:
        modeller = app.Modeller(pdb.topology, pdb.positions)
        if _COFACTOR_MOL2_PATHS:
            try:
                from run_restrained_md import inject_cofactor_bonds  # type: ignore
                inject_cofactor_bonds(modeller, _COFACTOR_MOL2_PATHS)
            except ImportError:
                pass
        if binder_chain:
            modeller, _cyclic_info = _apply_cyclic_bonds_mmgbsa(modeller, binder_chain)
        else:
            modeller.addHydrogens(ff)

        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod  = app.NoCutoff,
            constraints      = HBonds,
            soluteDielectric = 1.0,
            solventDielectric= 78.5,
        )
    except Exception as e:
        print(f"  [!] 시스템 생성 실패: {e}")
        return None

    system = apply_gb_radius_override(system, modeller.topology, ncaa_elem, si_radius)

    integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
    _pref = os.environ.get("UPDD_MMGBSA_PLATFORM", "CUDA")
    platform = None
    for name in (_pref, "CUDA", "OpenCL", "CPU"):
        try:
            platform = mm.Platform.getPlatformByName(name)
            break
        except Exception:
            continue
    if platform is None:
        platform = mm.Platform.getPlatformByName("CPU")
    sim = Simulation(modeller.topology, system, integrator, platform)

    # Optional position injection (for 1-traj bound-state positions). Fallback
    # to modeller.positions when no override is provided.
    if override_positions is not None:
        sim.context.setPositions(override_positions)
    else:
        sim.context.setPositions(modeller.positions)

    # NO minimization — static energy only.
    state = sim.context.getState(getEnergy=True)
    e_kj  = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    e_kcal = e_kj / 4.184
    return e_kcal


def calc_mmgbsa_1traj(pdb_path, output_dir, ff, ncaa_elem, si_radius,
                      receptor_chain="A", binder_chain="B"):
    """[v0.6.8 Fix 4] 1-trajectory MM-GBSA for a single snapshot.

    Algorithm (Genheden & Ryde 2015 §2.1):
      1. Split complex → complex / receptor / ligand tmp PDBs (split_complex).
      2. Build complex System (amber14 + GBn2 + HBonds + NoCutoff).
      3. Minimize complex: UPDD_MMGBSA_MIN_ITER (default 5000) / tol (default 1.0).
      4. Extract minimized positions from the complex context.
      5. Build receptor / ligand subsystems (same FF, only their atoms).
      6. Transfer minimized positions to each subsystem via
         (chain, resSeq, atomName) keying (identical to
         mmgbsa_1traj_compare._transfer_positions).
      7. Evaluate E_rec_static, E_lig_static via calc_energy_static.
      8. ΔG_1traj = E_complex_min − E_rec_static − E_lig_static.

    Result dict keys:
        snapshot, pdb_path, ncaa_element, si_radius_A,
        e_complex_kcal, e_receptor_kcal, e_ligand_kcal, delta_g_kcal,
        favorable, protocol_variant="1traj", minimize_iter, minimize_tol,
        platform.
    """
    basename = os.path.basename(pdb_path).replace(".pdb", "")
    tmp_dir  = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    print(f"\n  계산 [1-traj]: {basename}")

    # ── 1. Complex / Receptor / Ligand 분리 (reuse split_complex) ──
    complex_lines, receptor_lines, ligand_lines = split_complex_1traj(
        pdb_path, receptor_chain, binder_chain,
    )
    if not receptor_lines:
        print(f"  [!] Receptor (chain {receptor_chain}) 없음 — 건너뜀")
        return None
    if not ligand_lines:
        print(f"  [!] Ligand (chain {binder_chain}) 없음 — 건너뜀")
        return None

    cplx_pdb = os.path.join(tmp_dir, f"{basename}_complex.pdb")
    recv_pdb = os.path.join(tmp_dir, f"{basename}_receptor.pdb")
    lig_pdb  = os.path.join(tmp_dir, f"{basename}_ligand.pdb")

    _min_iter = int(os.environ.get("UPDD_MMGBSA_MIN_ITER", "5000"))
    _min_tol  = float(os.environ.get("UPDD_MMGBSA_MIN_TOL", "1.0"))

    try:
        write_temp_pdb(complex_lines,  cplx_pdb)
        write_temp_pdb(receptor_lines, recv_pdb)
        write_temp_pdb(ligand_lines,   lig_pdb)

        # ── 2. Complex 빌드 + 최소화 ───────────────────────────────
        # Cap filter for the complex (mirrors calc_energy).
        load_path = cplx_pdb
        _cap_tmp = None
        try:
            with open(cplx_pdb, "r", encoding="utf-8") as _f:
                _lines = _f.readlines()
            has_cap = any(
                l[17:20].strip() in _CAP_RESNAMES
                for l in _lines if l.startswith(("ATOM", "HETATM")) and len(l) >= 20
            )
            if has_cap:
                filtered = [l for l in _lines if not (
                    l.startswith(("ATOM", "HETATM")) and len(l) >= 20
                    and l[17:20].strip() in _CAP_RESNAMES
                )]
                _cap_tmp = cplx_pdb + ".nocap.tmp"
                with open(_cap_tmp, "w", encoding="utf-8") as _f:
                    _f.writelines(filtered)
                load_path = _cap_tmp
        except Exception as _cap_err:
            print(f"    [!] Cap 필터링 실패(무시): {_cap_err}")

        try:
            pdb = PDBFile(load_path)
        except Exception as e:
            print(f"  [!] Complex PDB 로드 실패: {e}")
            return None
        finally:
            if _cap_tmp and os.path.exists(_cap_tmp):
                try:
                    os.remove(_cap_tmp)
                except OSError:
                    pass

        cplx_modeller = app.Modeller(pdb.topology, pdb.positions)
        if _COFACTOR_MOL2_PATHS:
            try:
                from run_restrained_md import inject_cofactor_bonds  # type: ignore
                inject_cofactor_bonds(cplx_modeller, _COFACTOR_MOL2_PATHS)
            except ImportError:
                pass
        if binder_chain:
            cplx_modeller, _cyclic_info = _apply_cyclic_bonds_mmgbsa(
                cplx_modeller, binder_chain,
            )

        try:
            cplx_system = ff.createSystem(
                cplx_modeller.topology,
                nonbondedMethod  = app.NoCutoff,
                constraints      = HBonds,
                soluteDielectric = 1.0,
                solventDielectric= 78.5,
            )
        except Exception as e:
            print(f"  [!] Complex 시스템 생성 실패: {e}")
            return None
        cplx_system = apply_gb_radius_override(
            cplx_system, cplx_modeller.topology, ncaa_elem, si_radius,
        )

        _pref = os.environ.get("UPDD_MMGBSA_PLATFORM", "CUDA")
        platform = None
        for name in (_pref, "CUDA", "OpenCL", "CPU"):
            try:
                platform = mm.Platform.getPlatformByName(name)
                break
            except Exception:
                continue
        if platform is None:
            platform = mm.Platform.getPlatformByName("CPU")

        integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
        sim = Simulation(cplx_modeller.topology, cplx_system, integrator, platform)
        sim.context.setPositions(cplx_modeller.positions)

        print(f"  [1-traj] Complex 최소화 중 ({_min_iter} iter, tol={_min_tol})...")
        sim.minimizeEnergy(maxIterations=_min_iter, tolerance=_min_tol)
        min_state = sim.context.getState(getEnergy=True, getPositions=True)
        e_complex = min_state.getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole) / 4.184
        cplx_min_positions = min_state.getPositions(asNumpy=False)

        # Build key → complex atom index lookup for position transfer.
        def _atom_key(atom):
            res = atom.residue
            chain_id = res.chain.id if res.chain is not None else " "
            return (chain_id, res.name, str(res.id), atom.name)

        cplx_index_map = {
            _atom_key(a): a.index for a in cplx_modeller.topology.atoms()
        }

        # ── 3. Receptor subsystem static energy at bound geometry ──
        def _load_subsystem(sub_pdb_path, sub_binder_chain):
            """Load, addHydrogens (or cyclic-inject), return (modeller, system).
            Returns (None, None) on failure.
            """
            sub_load_path = sub_pdb_path
            sub_cap_tmp = None
            try:
                with open(sub_pdb_path, "r", encoding="utf-8") as _f:
                    _sub_lines = _f.readlines()
                has_cap = any(
                    l[17:20].strip() in _CAP_RESNAMES
                    for l in _sub_lines if l.startswith(("ATOM", "HETATM")) and len(l) >= 20
                )
                if has_cap:
                    filtered = [l for l in _sub_lines if not (
                        l.startswith(("ATOM", "HETATM")) and len(l) >= 20
                        and l[17:20].strip() in _CAP_RESNAMES
                    )]
                    sub_cap_tmp = sub_pdb_path + ".nocap.tmp"
                    with open(sub_cap_tmp, "w", encoding="utf-8") as _f:
                        _f.writelines(filtered)
                    sub_load_path = sub_cap_tmp
            except Exception:
                pass

            try:
                sub_pdb = PDBFile(sub_load_path)
            except Exception as e:
                print(f"  [!] Subsystem PDB 로드 실패 ({sub_pdb_path}): {e}")
                return None, None
            finally:
                if sub_cap_tmp and os.path.exists(sub_cap_tmp):
                    try:
                        os.remove(sub_cap_tmp)
                    except OSError:
                        pass

            sub_modeller = app.Modeller(sub_pdb.topology, sub_pdb.positions)
            if _COFACTOR_MOL2_PATHS:
                try:
                    from run_restrained_md import inject_cofactor_bonds  # type: ignore
                    inject_cofactor_bonds(sub_modeller, _COFACTOR_MOL2_PATHS)
                except ImportError:
                    pass
            if sub_binder_chain:
                sub_modeller, _ = _apply_cyclic_bonds_mmgbsa(
                    sub_modeller, sub_binder_chain,
                )
            else:
                sub_modeller.addHydrogens(ff)

            try:
                sub_system = ff.createSystem(
                    sub_modeller.topology,
                    nonbondedMethod  = app.NoCutoff,
                    constraints      = HBonds,
                    soluteDielectric = 1.0,
                    solventDielectric= 78.5,
                )
            except Exception as e:
                print(f"  [!] Subsystem 시스템 생성 실패: {e}")
                return None, None
            return sub_modeller, sub_system

        def _transfer_positions(sub_modeller, label):
            sub_positions = list(sub_modeller.positions)
            n_matched = 0
            n_fallback = 0
            for atom in sub_modeller.topology.atoms():
                key = _atom_key(atom)
                if key in cplx_index_map:
                    sub_positions[atom.index] = cplx_min_positions[cplx_index_map[key]]
                    n_matched += 1
                else:
                    n_fallback += 1
            if n_fallback > 0:
                print(f"    [{label}] position transfer: {n_matched} matched, "
                      f"{n_fallback} fallback (addHydrogens edge)")
            return sub_positions

        # Receptor
        recv_modeller, recv_system = _load_subsystem(recv_pdb, sub_binder_chain="")
        if recv_modeller is None:
            return None
        recv_system = apply_gb_radius_override(
            recv_system, recv_modeller.topology, "none", si_radius,
        )
        recv_positions = _transfer_positions(recv_modeller, "receptor")
        recv_integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
        recv_sim = Simulation(
            recv_modeller.topology, recv_system, recv_integrator, platform,
        )
        recv_sim.context.setPositions(recv_positions)
        e_receptor = recv_sim.context.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(unit.kilojoules_per_mole) / 4.184

        # Ligand
        lig_modeller, lig_system = _load_subsystem(lig_pdb, sub_binder_chain=binder_chain)
        if lig_modeller is None:
            return None
        lig_system = apply_gb_radius_override(
            lig_system, lig_modeller.topology, ncaa_elem, si_radius,
        )
        lig_positions = _transfer_positions(lig_modeller, "ligand")
        lig_integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
        lig_sim = Simulation(
            lig_modeller.topology, lig_system, lig_integrator, platform,
        )
        lig_sim.context.setPositions(lig_positions)
        e_ligand = lig_sim.context.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(unit.kilojoules_per_mole) / 4.184

        # ── 4. ΔG 계산 ──────────────────────────────────────────
        delta_g = e_complex - e_receptor - e_ligand

        print(f"  [1-traj] E_complex  : {e_complex:>12.4f} kcal/mol (minimized)")
        print(f"  [1-traj] E_receptor : {e_receptor:>12.4f} kcal/mol (static)")
        print(f"  [1-traj] E_ligand   : {e_ligand:>12.4f} kcal/mol (static)")
        print(f"  [1-traj] ΔG_bind    : {delta_g:>+12.4f} kcal/mol  "
              f"{'✅ 결합 유리' if delta_g < 0 else '⚠️  결합 불리'}")

        result = {
            "snapshot":         basename,
            "pdb_path":         pdb_path,
            "ncaa_element":     ncaa_elem,
            "si_radius_A":      si_radius,
            "e_complex_kcal":   e_complex,
            "e_receptor_kcal":  e_receptor,
            "e_ligand_kcal":    e_ligand,
            "delta_g_kcal":     delta_g,
            "favorable":        delta_g < 0,
            "protocol_variant": "1traj",
            "minimize_iter":    _min_iter,
            "minimize_tol":     _min_tol,
            "platform":         platform.getName(),
        }

        out_json = os.path.join(output_dir, f"{basename}_mmgbsa.json")
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)

        return result
    finally:
        for tmp in [cplx_pdb, recv_pdb, lig_pdb]:
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except OSError as _rm_err:
                    print(f"  [!] 임시 파일 정리 실패(무시): {tmp} ({_rm_err})")


# ==========================================
# [v0.6.8] Protocol dispatcher (reads UPDD_MMGBSA_PROTOCOL)
# ==========================================
_VALID_PROTOCOLS = ("3traj", "1traj")


def _resolve_protocol() -> str:
    """Read UPDD_MMGBSA_PROTOCOL env var. Returns "3traj" (default) or "1traj".
    Raises ValueError on any other value.
    """
    val = os.environ.get("UPDD_MMGBSA_PROTOCOL", "3traj").strip().lower()
    if val not in _VALID_PROTOCOLS:
        raise ValueError(
            f"UPDD_MMGBSA_PROTOCOL={val!r} invalid. "
            f"Expected one of {_VALID_PROTOCOLS}. Unset to use default 3traj."
        )
    return val


def calc_mmgbsa_dispatch(pdb_path, output_dir, ff, ncaa_elem, si_radius,
                         receptor_chain="A", binder_chain="B"):
    """[v0.6.8] Dispatch to 1-traj or 3-traj path based on UPDD_MMGBSA_PROTOCOL.

    - UPDD_MMGBSA_PROTOCOL unset or "3traj": call calc_mmgbsa() (existing).
    - UPDD_MMGBSA_PROTOCOL="1traj": call calc_mmgbsa_1traj() (new).
    - Any other value: ValueError.
    """
    protocol = _resolve_protocol()
    if protocol == "1traj":
        return calc_mmgbsa_1traj(
            pdb_path, output_dir, ff, ncaa_elem, si_radius,
            receptor_chain, binder_chain,
        )
    return calc_mmgbsa(
        pdb_path, output_dir, ff, ncaa_elem, si_radius,
        receptor_chain, binder_chain,
    )


# ==========================================
# CPU 병렬 워커 (v0.2 — multiprocessing.Pool)
# ==========================================
def _mmgbsa_worker(worker_args):
    """multiprocessing.Pool 워커. ForceField 는 pickle 불가능할 수 있어
    경로만 전달받고 워커 내부에서 재구성한다. 한 스냅샷의 실패가 다른
    스냅샷을 중단시키지 않도록 예외를 포착하여 None 을 반환한다.
    """
    (pdb_path, outputdir, ff_files, ncaa_xmls, ncaa_hydrogens,
     ncaa_elem, si_radius, receptor_chain, binder_chain,
     cofactor_mol2_paths) = worker_args
    try:
        # ncAA Hydrogens 는 프로세스별로 한 번만 로드
        for h_xml in ncaa_hydrogens:
            try:
                app.Modeller.loadHydrogenDefinitions(h_xml)
            except Exception:
                pass
        try:
            ff = ForceField(*ff_files, *ncaa_xmls)
        except Exception:
            ff = ForceField(*ff_files)
        # [R-17 G6 FIX] per-process GAFF template registration for cofactors.
        # Also set module-global so calc_energy can inject bonds per PDB.
        global _COFACTOR_MOL2_PATHS
        _COFACTOR_MOL2_PATHS = list(cofactor_mol2_paths or [])
        if _COFACTOR_MOL2_PATHS:
            try:
                from run_restrained_md import register_cofactor_templates  # type: ignore
                register_cofactor_templates(ff, _COFACTOR_MOL2_PATHS)
            except ImportError:
                pass
        return calc_mmgbsa_dispatch(
            pdb_path, outputdir, ff, ncaa_elem, si_radius, receptor_chain, binder_chain
        )
    except Exception as _werr:
        print(f"  [!] 워커 실패 ({os.path.basename(pdb_path)}): {_werr}")
        return None


# ==========================================
# CLI 메인
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description="MM-GBSA ΔG — OpenMM GBn2, ncAA GB 반경 자동 보정"
    )
    parser.add_argument("--md_dir",    required=True,        help="MD 결과 디렉토리 (final PDB)")
    parser.add_argument("--outputdir", required=True,        help="MM-GBSA 결과 출력 디렉토리")
    parser.add_argument("--ncaa_elem", default="none",       help="ncAA 핵심 원소 (예: Si, F)")
    parser.add_argument("--si_radius", type=float, default=2.10, help="Si GB 반경 Å (기본 2.10)")
    parser.add_argument("--receptor_chain", default="A",     help="Receptor chain ID (기본 A)")
    parser.add_argument("--binder_chain", default="B",       help="Binder chain ID (기본 B)")
    parser.add_argument("--target_id", default="",
                        help="[R-17 G6] target_card stem (예: 6WGN). 지정 시 "
                             "GNP/ATP 등 cofactor GAFF 템플릿을 OpenMM ForceField 에 "
                             "등록하여 implicit-solvent createSystem 이 성공하도록 한다.")
    args = parser.parse_args()

    os.makedirs(args.outputdir, exist_ok=True)

    # [v39 FIX — Bug 2] ForceField 로드 (ncAA XML 자동 포함)
    # 기존 `amber14-all.xml + amber14/tip3pfb.xml` 조합은 OpenMM 8.x 에서
    # `createSystem(implicitSolvent=GBn2)` 를 무시하여
    # `The argument 'implicitSolvent' was specified to createSystem() but was never used`
    # 에러를 유발한다 (모든 Complex/Receptor/Ligand 에너지가 silently 실패).
    # 올바른 조합은 `amber14-all.xml + implicit/gbn2.xml` 이며, 이 경우
    # createSystem 은 implicitSolvent 인자 없이 GBn2 파라미터를 자동 적용한다.
    # tip3pfb.xml 은 explicit water 정의로 implicit 계산과 무관하므로 제거한다.
    # [SciVal Verdict 18, 2026-04-20] Li-Merz 12-6-4 divalent ion FF via
    # OpenMM-shipped amber14/tip3p_HFE_multivalent.xml. Provides residue
    # templates + 12-6-4 r⁻⁴ polarization term for Mg/Zn/Ca/Mn/Fe/Cu/Co/Ni.
    # Mandatory for KRAS-GNP where Mg²⁺ screens the β,γ-phosphate charge
    # in the receptor pool (v39 FIX metal-strip rollback per verdict §5.3).
    ff_files  = [
        "amber14-all.xml",
        "amber/tip3p_HFE_multivalent.xml",
        "implicit/gbn2.xml",
    ]
    # [v55 FIX] manifest 기반 ncAA XML 로딩.
    # 기존 glob("*_resp.xml")은 GAFF2 파라미터화(resp_used=false)일 때
    # NMA_gaff2.xml을 찾지 못해 ForceField에 ncAA 템플릿이 누락되었다.
    # params_manifest.json의 xml_path/hydrogens_path를 사용하여
    # 파라미터화 방식에 관계없이 올바른 XML을 로드한다.
    # [Option C v2] Cache indirection deprecation — prefer LOCAL params/<resname>_gaff2.xml
    # over manifest xml_path which may point to /tmp/calib_params/ (volatile).
    # Alias-aware: ncaa_code (e.g., NML) may differ from xml_resname (e.g., MLE);
    # the on-disk filename uses xml_resname. Falls back to manifest path with
    # a warning when local file is missing (legacy compat).
    import warnings
    params_dir = os.path.join(os.path.dirname(args.md_dir), "params")
    ncaa_xmls = []
    ncaa_hydrogens = []
    for mf in sorted(glob.glob(os.path.join(params_dir, "*_params_manifest.json"))):
        try:
            with open(mf, "r", encoding="utf-8") as _mf:
                manifest = json.load(_mf)
            ncaa_code = manifest.get("ncaa_code", "")
            xml_resname = manifest.get("xml_resname", ncaa_code)  # NML→MLE alias
            if not xml_resname:
                continue
            local_xml = os.path.join(params_dir, f"{xml_resname}_gaff2.xml")
            if os.path.exists(local_xml):
                ncaa_xmls.append(local_xml)
            else:
                xml_path = manifest.get("xml_path", "")
                if xml_path and os.path.exists(xml_path):
                    warnings.warn(
                        f"[Option C legacy fallback] {ncaa_code} (xml_resname={xml_resname}): "
                        f"local {local_xml} missing; using manifest xml_path={xml_path}.",
                        stacklevel=2,
                    )
                    ncaa_xmls.append(xml_path)
            # Hydrogens: alias-aware prefer-local-with-fallback
            local_h = os.path.join(params_dir, f"{xml_resname}_hydrogens.xml")
            if os.path.exists(local_h):
                ncaa_hydrogens.append(local_h)
            else:
                h_path = manifest.get("hydrogens_path", "")
                if h_path and os.path.exists(h_path):
                    ncaa_hydrogens.append(h_path)
        except Exception as _me:
            print(f"  [!] manifest 파싱 실패 ({os.path.basename(mf)}): {_me}")
    if ncaa_xmls:
        print(f"  ncAA XML 로드: {[os.path.basename(x) for x in ncaa_xmls]}")
    if ncaa_hydrogens:
        for h_xml in ncaa_hydrogens:
            app.Modeller.loadHydrogenDefinitions(h_xml)
        print(f"  ncAA Hydrogens 로드: {[os.path.basename(x) for x in ncaa_hydrogens]}")
    try:
        ff = ForceField(*ff_files, *ncaa_xmls)
    except Exception as _ff_err:
        print(f"  [!] ncAA XML 로드 실패 ({_ff_err.__class__.__name__}: {_ff_err}). "
              f"표준 amber14 only 로 폴백합니다 — 결과 해석 시 주의!")
        ff = ForceField(*ff_files)

    # [R-17 G6] Load cofactor mol2 paths from target_card (for workers).
    cofactor_mol2_paths = []
    if args.target_id:
        try:
            from run_restrained_md import load_cofactor_ff_parameters  # type: ignore
            from utils_common import load_target_card  # type: ignore
            _cards_dir = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                os.pardir, "target_cards",
            )
            _card = load_target_card(args.target_id, cards_dir=_cards_dir)
            cofactor_mol2_paths = [
                (r, p) for r, p in load_cofactor_ff_parameters(_card)
                if p.endswith(".mol2")
            ]
            if cofactor_mol2_paths:
                print(f"  [R-17 G6] cofactor mol2 감지: "
                      f"{[r for r,_ in cofactor_mol2_paths]} — 워커에 전달")
                # Also register on the main-process ff (for any non-worker paths)
                try:
                    from run_restrained_md import register_cofactor_templates  # type: ignore
                    register_cofactor_templates(ff, cofactor_mol2_paths)
                    global _COFACTOR_MOL2_PATHS
                    _COFACTOR_MOL2_PATHS = list(cofactor_mol2_paths)
                except ImportError:
                    pass
        except (FileNotFoundError, ValueError, ImportError) as _tc_err:
            print(f"  [R-17 G6][WARN] target_card '{args.target_id}' 로드 실패 "
                  f"({_tc_err}) — cofactor FF 등록 생략")

    # final PDB 파일 탐색
    pdb_files = sorted(glob.glob(os.path.join(args.md_dir, "*_final.pdb")))
    if not pdb_files:
        pdb_files = sorted(glob.glob(os.path.join(args.md_dir, "*.pdb")))
    if not pdb_files:
        print(f"[!] PDB 파일 없음: {args.md_dir}")
        return

    print(f"\n[MM-GBSA ΔG Calculation]")
    print(f"  대상 PDB  : {len(pdb_files)}개")
    print(f"  GB 모델   : GBn2")
    print(f"  ncAA 원소 : {args.ncaa_elem}")
    print(f"  Si 반경   : {args.si_radius} Å")

    # [v0.2] CPU 멀티프로세싱 — ForceField 객체는 pickle 불가능할 수 있으므로
    # ff_files/ncaa_xmls/ncaa_hydrogens 경로 리스트를 워커에 전달하여
    # 각 워커가 독립적으로 ForceField 를 재구성한다. 한 스냅샷 실패가 전체를
    # 중단시키지 않도록 try/except 로 감싼다.
    import multiprocessing as _mp

    _worker_args_list = [
        (
            pdb_path, args.outputdir,
            list(ff_files), list(ncaa_xmls), list(ncaa_hydrogens),
            args.ncaa_elem, args.si_radius, args.receptor_chain, args.binder_chain,
            list(cofactor_mol2_paths),
        )
        for pdb_path in pdb_files
    ]

    n_workers = min(4, len(pdb_files)) or 1
    print(f"  병렬 워커 : {n_workers} (multiprocessing.Pool)")

    all_results = []
    try:
        with _mp.Pool(processes=n_workers) as pool:
            for res in pool.imap_unordered(_mmgbsa_worker, _worker_args_list):
                if res is not None:
                    all_results.append(res)
    except Exception as _pool_err:
        print(f"  [!] 병렬 처리 실패 ({_pool_err}). 순차 실행으로 폴백합니다.")
        for wa in _worker_args_list:
            res = _mmgbsa_worker(wa)
            if res is not None:
                all_results.append(res)

    if not all_results:
        print("\n[!] 계산된 결과 없음")
        return

    # ΔG 순 정렬
    all_results.sort(key=lambda r: r["delta_g_kcal"])

    # [v4 S-1] Interaction Entropy (Duan 2016) — 동일 design 의 스냅샷 묶음에 대해
    # ΔG 분포로부터 −TΔS 를 추정하고 보정 ΔG 를 계산한다.
    #
    # 그룹화 키: snapshot 이름의 "_snap" 접미사 이전 prefix (= design id).
    # 동일 design 의 스냅샷이 5개 이상 있으면 reliable=True 로 보고된다.
    from collections import defaultdict as _dd
    by_design = _dd(list)
    for r in all_results:
        snap_name = r.get("snapshot", "")
        design_key = snap_name.split("_snap")[0] if "_snap" in snap_name else snap_name
        by_design[design_key].append(r["delta_g_kcal"])

    interaction_entropy_per_design = {}
    for design_key, dg_list in by_design.items():
        ie_result = calc_interaction_entropy(dg_list, temperature=300.0)
        interaction_entropy_per_design[design_key] = ie_result

    # [v0.6.8 Fix 4] Per-snap protocol_variant already set by calc_mmgbsa /
    # calc_mmgbsa_1traj. Top-level summary field records "1traj" / "3traj" /
    # "mixed" so consumers can short-circuit legacy compat logic.
    _protocols_seen = {r.get("protocol_variant", "3traj") for r in all_results}
    if len(_protocols_seen) == 1:
        _top_protocol = _protocols_seen.pop()
    else:
        _top_protocol = "mixed"

    # 전체 요약 JSON
    summary_path = os.path.join(args.outputdir, "mmgbsa_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "schema_version":   "0.6.8",  # [v0.6.8] bumped from 0.6.7 (Fix 4)
            "protocol_variant": _top_protocol,
            "gb_model":   "GBn2",
            "ncaa_elem":  args.ncaa_elem,
            "si_radius":  args.si_radius,
            "n_calc":     len(all_results),
            "mean_dg":    float(np.mean([r["delta_g_kcal"] for r in all_results])),
            "best_dg":    float(all_results[0]["delta_g_kcal"]),
            "results":    all_results,
            "interaction_entropy_per_design": interaction_entropy_per_design,
            "interaction_entropy_method":     "Duan 2016 JCTC, T=300K",
        }, f, indent=2)

    # [v4 S-1] IE 보정 결과 콘솔 출력 (reliable 그룹만)
    reliable_ie = [(k, v) for k, v in interaction_entropy_per_design.items() if v.get("reliable")]
    if reliable_ie:
        print(f"\n  [Interaction Entropy 보정 (Duan 2016, N≥5)]")
        for design_key, ie in sorted(reliable_ie, key=lambda x: x[1]["delta_g_corrected_kcal"]):
            print(f"    {design_key[:40]:<42} <ΔG>={ie['mean_delta_g_kcal']:+8.3f}  "
                  f"−TΔS={ie['minus_TdS_kcal']:+7.3f}  ΔG_corr={ie['delta_g_corrected_kcal']:+8.3f}  "
                  f"(N={ie['n_frames']})")

    print(f"\n  [ΔG_bind 순위 (낮을수록 강한 결합)]")
    print(f"  {'순위':<4} {'스냅샷':<42} {'ΔG (kcal/mol)':>14}")
    print("  " + "-" * 62)
    for i, r in enumerate(all_results, 1):
        flag = "✅" if r["favorable"] else "⚠️ "
        print(f"  {i:<4} {r['snapshot'][:40]:<42} {r['delta_g_kcal']:>+12.4f}  {flag}")

    print(f"\n  평균 ΔG : {np.mean([r['delta_g_kcal'] for r in all_results]):+.4f} kcal/mol")
    print(f"  최선 ΔG : {all_results[0]['delta_g_kcal']:+.4f} kcal/mol")
    print(f"  요약 저장: {summary_path}")


if __name__ == "__main__":
    main()
