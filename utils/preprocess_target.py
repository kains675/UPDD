#!/usr/bin/env python
"""
preprocess_target.py
--------------------
타겟 PDB 자동 전처리 (PDBFixer 완벽 통합 버전)
- 선택 체인 유지 / 물 제거 / HETATM 자동·수동 보존
- [핵심 수정] 원본 잔기 번호 유지 (Renumbering 삭제 -> 핫스팟 번호 꼬임 원천 차단)
- [핵심 수정] PDBFixer 자동 호출 (유실된 루프/원자 물리적 복구)
"""
import os
import sys
import csv
import argparse
import subprocess
from collections import defaultdict
import numpy as np

# [v3 12] 도메인 상수와 공통 헬퍼는 utils_common(SSoT) 에서 임포트한다.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils_common import (  # noqa: E402
    WATER_RESNAMES,
    COMMON_IONS,
    COMMON_COFACTORS,
    parse_pdb_atom_line,
    atom_distance,
)


def parse_keep_list(text):
    if not text:
        return set()
    return {x.strip().upper() for x in text.split(',') if x.strip()}


def atom_xyz(line):
    """[Backwards-Compat] PDB ATOM/HETATM 줄에서 (x, y, z) 튜플을 반환한다.

    내부적으로 utils_common.parse_pdb_atom_line 을 사용하여 파싱 정합성을
    유지한다. 파싱 실패 시 ValueError 를 raise 한다.
    """
    parsed = parse_pdb_atom_line(line)
    if parsed is None:
        raise ValueError(f"PDB 좌표 파싱 실패: {line[:30]!r}")
    return (parsed["x"], parsed["y"], parsed["z"])


def dist(a, b):
    """[Backwards-Compat Wrapper] utils_common.atom_distance 위임."""
    return atom_distance(a, b)


def _is_close_to_protein(het_coords_list, protein_coords_arr, cutoff=5.0):
    """HETATM 좌표 리스트와 단백질 좌표 배열 간 최소 거리 < cutoff 여부를 벡터 연산으로 판별.

    [v3 10-1] 기존 이중 for-loop 의 O(n_het × n_prot) 비용을 numpy 브로드캐스팅
    으로 O(n_het × n_prot) 단일 메모리 패스로 개선한다 (CPython 인터프리터 오버헤드 제거).

    Args:
        het_coords_list: HETATM 잔기의 (x, y, z) 튜플 리스트
        protein_coords_arr: 단백질 원자 좌표 numpy 배열 (shape (N_prot, 3))
        cutoff: Å 단위 거리 임계값

    Returns:
        bool: 어떤 HETATM 원자라도 단백질 원자와 cutoff 이내의 거리에 있으면 True.
    """
    if len(protein_coords_arr) == 0 or len(het_coords_list) == 0:
        return False
    het_arr = np.asarray(het_coords_list, dtype=float)
    # 각 HETATM 원자에 대해 모든 단백질 원자와의 거리를 계산
    # diff: shape (n_het, n_prot, 3); norm 후 shape (n_het, n_prot)
    diff = het_arr[:, None, :] - protein_coords_arr[None, :, :]
    dists = np.linalg.norm(diff, axis=2)
    return bool(np.any(dists <= cutoff))

def collect_records(path):
    """[v3 3-4] PDB 레코드를 utils_common.parse_pdb_atom_line 으로 일관 파싱한다."""
    protein_lines = []
    protein_atoms = []
    hetatm_lines = []
    passthrough = []

    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('ATOM'):
                protein_lines.append(line)
                parsed = parse_pdb_atom_line(line)
                if parsed is not None:
                    protein_atoms.append({
                        'xyz': (parsed["x"], parsed["y"], parsed["z"]),
                        'chain': parsed["chain"] or line[21],  # chain 컬럼 공백 보존
                        'resid': str(parsed["resnum"]),
                        'resname': parsed["resname"],
                    })
            elif line.startswith('HETATM'):
                hetatm_lines.append(line)
            elif line.startswith(('CONECT', 'MASTER', 'END')):
                continue
            else:
                passthrough.append(line)
    return protein_lines, protein_atoms, hetatm_lines, passthrough


def filter_hetatms(hetatm_lines, protein_atoms, mode, manual_keeps):
    """[v3 10-1] HETATM 분류기.

    핵심 개선:
        - 단백질-HETATM 거리 계산을 numpy 브로드캐스팅으로 벡터화 (이중 루프 제거).
        - 단백질 좌표 배열은 함수 진입 시 한 번만 numpy 변환되어 모든 HETATM
          잔기에 재사용된다 (메모리 지역성 개선).
    """
    kept_lines = []
    kept_keys = set()
    removed_het = []

    # 단백질 좌표를 한 번에 numpy 배열로 변환
    protein_coords_arr = np.asarray([pa['xyz'] for pa in protein_atoms], dtype=float) \
        if protein_atoms else np.empty((0, 3), dtype=float)

    grouped_het = defaultdict(list)
    for line in hetatm_lines:
        chain = line[21]
        resid = line[22:26].strip()
        resname = line[17:20].strip()
        grouped_het[(chain, resid, resname)].append(line)

    if mode == 'none':
        for k, lines in grouped_het.items():
            removed_het.append((k[0], k[1], k[2], "mode_none"))
        return [], set(), removed_het

    for (chain, resid, resname), lines in grouped_het.items():
        # 물 분자 무조건 제거
        if resname in WATER_RESNAMES:
            removed_het.append((chain, resid, resname, "water"))
            continue

        if mode == 'manual':
            if resname in manual_keeps:
                kept_lines.extend(lines)
                kept_keys.add((chain, resid, resname, "manual_keep"))
            else:
                removed_het.append((chain, resid, resname, "manual_exclude"))
            continue

        if mode == 'auto':
            if resname in manual_keeps:
                kept_lines.extend(lines)
                kept_keys.add((chain, resid, resname, "manual_override"))
                continue

            if resname in COMMON_IONS or resname in COMMON_COFACTORS:
                kept_lines.extend(lines)
                kept_keys.add((chain, resid, resname, "common_cofactor_or_ion"))
                continue

            # 단백질과 5Å 이내에 있는지 벡터 거리 계산
            het_coords = []
            for hl in lines:
                parsed = parse_pdb_atom_line(hl)
                if parsed is not None:
                    het_coords.append((parsed["x"], parsed["y"], parsed["z"]))
            is_close = _is_close_to_protein(het_coords, protein_coords_arr, cutoff=5.0)

            if is_close:
                kept_lines.extend(lines)
                kept_keys.add((chain, resid, resname, "distance<=5A"))
            else:
                removed_het.append((chain, resid, resname, "distance>5A"))

    return kept_lines, kept_keys, removed_het

def main():
    parser = argparse.ArgumentParser(description="PDB 전처리 스크립트 (PDBFixer 포함)")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chains", default="A", help="유지할 체인 (예: A,B)")
    parser.add_argument("--hetatm_mode", choices=["auto", "manual", "none"], default="auto")
    parser.add_argument("--keep_hetatms", default="", help="강제 유지할 HETATM 이름 (예: GDP,MG)")
    parser.add_argument("--report", default="", help="리포트 CSV 출력 경로")
    args = parser.parse_args()

    keep_chains = parse_keep_list(args.chains)
    manual_keeps = parse_keep_list(args.keep_hetatms)

    prot_lines, prot_atoms, het_lines, passthrough = collect_records(args.input)

    # 체인 필터링
    if keep_chains:
        prot_lines = [l for l in prot_lines if l[21] in keep_chains]
        prot_atoms = [a for a in prot_atoms if a['chain'] in keep_chains]
        het_lines = [l for l in het_lines if l[21] in keep_chains]

    # HETATM 필터링
    kept_het_lines, kept_het_keys, removed_het = filter_hetatms(
        het_lines, prot_atoms, args.hetatm_mode, manual_keeps
    )

    output_lines = []
    output_lines.extend(prot_lines)
    if output_lines and not output_lines[-1].startswith('TER'):
        output_lines.append('TER\n')
    output_lines.extend(kept_het_lines)
    output_lines.extend(passthrough) # SEQRES 등 헤더 유지 (PDBFixer가 참조함)
    output_lines.append('END\n')

    # 1. 1차 필터링 결과를 임시 파일로 저장
    temp_pdb = args.output + ".tmp"
    with open(temp_pdb, 'w', encoding='utf-8') as f:
        f.writelines(output_lines)

    # 2. PDBFixer 실행 (물리적 간극 수선 및 원자 복구)
    print("\n[+] PDBFixer를 통한 유실된 루프 및 원자 복구 시작...")
    pdbfixer_status = "applied"
    pdbfixer_error  = ""
    try:
        cmd = [
            "pdbfixer", temp_pdb,
            f"--output={args.output}",
            "--add-atoms=all",
            "--add-residues",
            "--keep-heterogens=all"
        ]
        # PDBFixer 호출
        subprocess.run(cmd, check=True)
        print(f"[+] PDBFixer 수술 완료: {args.output}")
        os.remove(temp_pdb)
    except Exception as e:
        pdbfixer_status = "fallback_no_repair"
        pdbfixer_error  = str(e)
        print(f"[!] PDBFixer 실행 실패 (환경 확인 필요): {e}")
        print(f"[!] 복구 없이 기본 필터링만 적용하여 저장합니다.")
        os.rename(temp_pdb, args.output)

    # 리포트 작성
    report = args.report or os.path.splitext(args.output)[0] + '_preprocess_report.csv'
    with open(report, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['section', 'chain', 'resid', 'icode', 'resname', 'value'])
        for chain, resid, resname, reason in kept_het_keys:
            w.writerow(['kept_het', chain, resid, '', resname, reason])
        for chain, resid, resname, reason in removed_het:
            w.writerow(['removed_het', chain, resid, '', resname, reason])
        # [Provenance] PDBFixer 적용 여부를 리포트에 기록하여 후속 분석가가
        # 누락된 루프/원자 복구의 출처를 추적할 수 있게 한다.
        w.writerow(['pdbfixer', '', '', '', '', pdbfixer_status])
        if pdbfixer_error:
            w.writerow(['pdbfixer_error', '', '', '', '', pdbfixer_error])

    print(f"[+] 전처리 리포트 생성 완료: {report}")
    print(f"[+] 최종 출력: {args.output}")

if __name__ == "__main__":
    main()
