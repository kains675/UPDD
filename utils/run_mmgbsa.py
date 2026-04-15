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


# ==========================================
# GB 반경 Override 테이블
# ==========================================
GB_RADII_OVERRIDE = {
    "Si": 2.10,
    "Br": 1.85,
    "Cl": 1.75,
    "F":  1.47,
    "P":  1.80,
}

GB_SCALE_OVERRIDE = {
    "Si": 0.80,
    "Br": 0.80,
    "Cl": 0.80,
    "F":  0.80,
    "P":  0.86,
}



# ==========================================
# PDB 분리 유틸 — Complex / Receptor / Ligand
# ==========================================
# [v39 FIX] implicit/gbn2.xml 은 물·카운터이온 템플릿을 포함하지 않으므로,
# MM-GBSA 에 투입하기 전에 반드시 제거해야 한다. explicit water 와 단순
# 카운터이온은 implicit 용매 모델이 연속 매질로 대체한다. 구조적 금속
# 이온 (MG, ZN, CA 등) 도 현재 implicit FF 에 템플릿이 없어 제거하되,
# 향후 커스텀 이온 XML 로 복원 가능하도록 별도 카운트를 남긴다.
# [v52 FIX] MD 말단 캡핑 잔기. PDBFixer가 MD 안정화를 위해 추가한 것이므로
# MM-GBSA implicit solvent 계산 전에 제거한다. 이 잔기들은 실제 단백질의
# 일부가 아니며, implicit/gbn2.xml에 템플릿이 없어 에러를 유발한다.
_CAP_RESNAMES = {"NME", "ACE", "NHE"}

_MMGBSA_STRIP_RESNAMES = {
    "HOH", "WAT", "TIP", "TIP3", "SOL",                    # water
    "NA", "CL", "K", "BR", "IOD",                           # counterions
    "MG", "ZN", "CA", "FE", "MN", "CU", "CO", "NI",        # structural metals
    "FE2", "FE3",
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

    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            rec = line[:6].strip()
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

            if chain == receptor_chain:
                receptor_lines.append(line)
            elif rec == "HETATM" and resname in KNOWN_COFACTORS:
                receptor_lines.append(line)
                n_cofactor_to_receptor += 1
            elif chain == binder_chain:
                ligand_lines.append(line)
            elif chain == " " or rec == "HETATM":
                ligand_lines.append(line)
                if chain == " ":
                    n_blank_chain += 1

    if n_stripped > 0:
        print(f"  [split_complex] implicit GB 호환: water/ion/metal {n_stripped}개 원자 제거")
    if n_cofactor_to_receptor > 0:
        print(f"  [split_complex] 보조인자 {n_cofactor_to_receptor}개 원자를 receptor 풀에 포함 (KNOWN_COFACTORS).")
    if n_blank_chain > 0:
        print(f"  [!] split_complex 경고: chain ID가 공백인 원자 {n_blank_chain}개를 ligand로 분류했습니다.")

    return complex_lines, receptor_lines, ligand_lines


def write_temp_pdb(lines, path):
    """임시 PDB 를 작성한다.

    [v44 FIX] CONECT 레코드를 제거하고, cyclic HTC bond 가 감지되면
    말단 잔기 좌표를 분리하여 OpenMM 의 자동 결합 검출이 cyclic bond 를
    형성하지 않도록 한다. 이는 MM-GBSA 에너지 계산에서 cyclic 펩타이드의
    템플릿 불일치를 방지한다 (amide bond 1.33 Å → auto-bond 검출 문턱 이내).
    cyclic bond 에너지 기여분은 complex 와 ligand 양쪽에서 동일하게 취급되므로
    ΔG_bind ranking 에 영향을 주지 않는다.
    """
    # CONECT 제거 + cyclic HTC 해제
    cleaned = []
    atom_lines_by_chain = {}  # chain → list of (idx_in_cleaned, resnum, atom_name, x, y, z)
    for line in lines:
        rec = line[:6].strip()
        if rec == "CONECT":
            continue
        cleaned.append(line)
        if rec in ("ATOM", "HETATM") and len(line) >= 54:
            chain = line[21:22]
            try:
                resnum = int(line[22:26])
                atom_name = line[12:16].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atom_lines_by_chain.setdefault(chain, []).append(
                    (len(cleaned) - 1, resnum, atom_name, x, y, z))
            except ValueError:
                pass

    # 각 chain 에서 cyclic N-C contact (< 2.0 Å) 를 감지하고 해제
    for chain, atoms in atom_lines_by_chain.items():
        if len(atoms) < 4:
            continue
        resnums = sorted({a[1] for a in atoms})
        if len(resnums) < 2:
            continue
        first_res, last_res = resnums[0], resnums[-1]
        n_info = next(((i, x, y, z) for i, rn, an, x, y, z in atoms
                       if rn == first_res and an == "N"), None)
        c_info = next(((i, x, y, z) for i, rn, an, x, y, z in atoms
                       if rn == last_res and an == "C"), None)
        if n_info and c_info:
            _, nx, ny, nz = n_info
            ci, cx, cy, cz = c_info
            dist = ((nx - cx)**2 + (ny - cy)**2 + (nz - cz)**2) ** 0.5
            if dist < 2.0:
                # cyclic HTC bond 감지 → 마지막 잔기 전체를 5 Å 이동 + OXT 원자 추가.
                # 이동만으로는 C-terminal OXT 부재로 CGLU 템플릿 불일치가 발생.
                c_atom = o_atom = ca_atom = None
                for ai, arn, aan, ax, ay, az in atoms:
                    if arn == last_res:
                        old_line = cleaned[ai]
                        cleaned[ai] = old_line[:30] + f"{ax + 5.0:8.3f}{ay:8.3f}{az:8.3f}" + old_line[54:]
                        if aan == "C":
                            c_atom = (ax + 5.0, ay, az)
                        elif aan == "O":
                            o_atom = (ax + 5.0, ay, az)
                        elif aan == "CA":
                            ca_atom = (ax + 5.0, ay, az)

                # OXT 원자 추가 (C-terminal carboxylate 의 두 번째 산소)
                # 위치: C 를 중심으로 O 의 반대편, C 에서 1.25 Å
                if c_atom and o_atom:
                    dx = c_atom[0] - o_atom[0]
                    dy = c_atom[1] - o_atom[1]
                    dz = c_atom[2] - o_atom[2]
                    d_co = max((dx**2 + dy**2 + dz**2)**0.5, 0.01)
                    scale = 1.25 / d_co
                    oxt_x = c_atom[0] + dx * scale
                    oxt_y = c_atom[1] + dy * scale
                    oxt_z = c_atom[2] + dz * scale
                    # 마지막 잔기의 chain/resnum 정보 추출
                    ref_line = None
                    for ai, arn, aan, _, _, _ in atoms:
                        if arn == last_res and aan == "O":
                            ref_line = cleaned[ai]
                            break
                    if ref_line:
                        serial = len([l for l in cleaned if l[:4] in ("ATOM", "HETA")]) + 1
                        # PDB ATOM 80-char 정확한 컬럼 포맷
                        buf = list(" " * 80)
                        buf[0:6]   = list("ATOM  ")
                        buf[6:11]  = list(f"{serial:5d}")
                        buf[12:16] = list(" OXT")
                        buf[17:20] = list(ref_line[17:20])  # resName
                        buf[21:22] = list(ref_line[21:22])  # chainID
                        buf[22:26] = list(ref_line[22:26])  # resSeq
                        buf[30:38] = list(f"{oxt_x:8.3f}")
                        buf[38:46] = list(f"{oxt_y:8.3f}")
                        buf[46:54] = list(f"{oxt_z:8.3f}")
                        buf[54:60] = list("  1.00")
                        buf[60:66] = list("  0.00")
                        buf[76:78] = list(" O")
                        oxt_line = "".join(buf).rstrip() + "\n"
                        # 마지막 잔기의 마지막 ATOM 행 바로 뒤에 삽입
                        last_atom_idx = max(ai for ai, arn, _, _, _, _ in atoms
                                            if arn == last_res)
                        cleaned.insert(last_atom_idx + 1, oxt_line)

    with open(path, "w", encoding="utf-8") as f:
        f.writelines(cleaned)
        if not cleaned or not cleaned[-1].startswith("END"):
            f.write("END\n")


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

    n_overridden = 0
    for force in system.getForces():
        force_type = type(force).__name__
        if "GB" not in force_type and "GBSA" not in force_type:
            continue

        for atom in topology.atoms():
            elem = atom.element.symbol if atom.element else "C"
            if elem in override_map:
                idx = atom.index
                try:
                    # GBSAOBCForce signature: (charge, radius, scale)
                    charge, radius, scale = force.getParticleParameters(idx)
                    new_r = override_map[elem] * unit.angstrom
                    new_s = GB_SCALE_OVERRIDE.get(elem, 0.80)
                    force.setParticleParameters(idx, charge, new_r, new_s)
                    n_overridden += 1
                except Exception as _obc_err:
                    # [v35 AUDIT] Principle 1: CustomGBForce fallback retained —
                    # it is a valid routing for GBn2/CustomGB systems. Principle 5:
                    # final failure path now logs rather than silently passing, so
                    # operators can trace GB-radius override misses that corrupt
                    # ΔG_bind magnitudes for Si/Br/Cl/F-bearing ncAAs.
                    try:
                        params = list(force.getParticleParameters(idx))
                        if len(params) >= 2:
                            params[1] = override_map[elem] * unit.angstrom
                            force.setParticleParameters(idx, params)
                            n_overridden += 1
                    except Exception as _custom_err:
                        print(f"  [!] GB radius override failed for atom idx={idx} elem={elem}: "
                              f"OBC={_obc_err.__class__.__name__}, Custom={_custom_err.__class__.__name__}")

    if n_overridden > 0:
        print(f"  GB 반경 Override: {n_overridden}개 원자 ({ncaa_elem} {override_map.get(ncaa_elem,'?')} Å)")
    return system


# ==========================================
# 단일 PDB → OpenMM 에너지 계산
# ==========================================
def calc_energy(pdb_path, ff, ncaa_elem, si_radius):
    """
    PDB → OpenMM implicit solvent 에너지 (kcal/mol)
    """
    # [v55 FIX] NME/ACE 캡 잔기를 PDB 로드 전에 제거.
    # split_complex에서 이미 필터하지만, 방어적으로 여기서도 필터.
    # [SciVal] 원본 PDB를 in-place 수정하지 않고 임시 파일을 사용한다.
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
    platform   = mm.Platform.getPlatformByName("CPU")
    sim = Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)

    # 에너지 최소화 (구조 충돌 방지)
    sim.minimizeEnergy(maxIterations=500, tolerance=10.0)

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
        print("  Complex  에너지 계산 중...")
        e_complex  = calc_energy(cplx_pdb, ff, ncaa_elem, si_radius)

        print("  Receptor 에너지 계산 중...")
        e_receptor = calc_energy(recv_pdb, ff, ncaa_elem="none", si_radius=si_radius)

        print("  Ligand   에너지 계산 중...")
        e_ligand   = calc_energy(lig_pdb,  ff, ncaa_elem, si_radius)

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
# CPU 병렬 워커 (v0.2 — multiprocessing.Pool)
# ==========================================
def _mmgbsa_worker(worker_args):
    """multiprocessing.Pool 워커. ForceField 는 pickle 불가능할 수 있어
    경로만 전달받고 워커 내부에서 재구성한다. 한 스냅샷의 실패가 다른
    스냅샷을 중단시키지 않도록 예외를 포착하여 None 을 반환한다.
    """
    (pdb_path, outputdir, ff_files, ncaa_xmls, ncaa_hydrogens,
     ncaa_elem, si_radius, receptor_chain, binder_chain) = worker_args
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
        return calc_mmgbsa(
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
    ff_files  = ["amber14-all.xml", "implicit/gbn2.xml"]
    # [v55 FIX] manifest 기반 ncAA XML 로딩.
    # 기존 glob("*_resp.xml")은 GAFF2 파라미터화(resp_used=false)일 때
    # NMA_gaff2.xml을 찾지 못해 ForceField에 ncAA 템플릿이 누락되었다.
    # params_manifest.json의 xml_path/hydrogens_path를 사용하여
    # 파라미터화 방식에 관계없이 올바른 XML을 로드한다.
    params_dir = os.path.join(os.path.dirname(args.md_dir), "params")
    ncaa_xmls = []
    ncaa_hydrogens = []
    for mf in sorted(glob.glob(os.path.join(params_dir, "*_params_manifest.json"))):
        try:
            with open(mf, "r", encoding="utf-8") as _mf:
                manifest = json.load(_mf)
            xml_path = manifest.get("xml_path", "")
            if xml_path and os.path.exists(xml_path):
                ncaa_xmls.append(xml_path)
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

    # 전체 요약 JSON
    summary_path = os.path.join(args.outputdir, "mmgbsa_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
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
