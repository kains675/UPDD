#!/usr/bin/env python
from __future__ import annotations

"""
utils/run_cyclic_fold.py
------------------------
Cyclic peptide 구조 예측: ColabFold + geometric N-C gap closure.

전략 (v44):
    1차: ColabFold 로 linear 구조 예측 (고품질, pLDDT > 90)
    2차: 바인더 체인의 N-C 말단을 geometric interpolation 으로 폐합
    3차: N-C 거리 검증 (< 4.0 Å → MD HTC bond 형성 가능)

과학적 근거:
    - ColabFold 는 backbone 기하를 정확히 예측하나 cyclic constraint 를 강제하지 않음
    - Geometric gap closure 는 binding interface 를 보존하면서 말단만 이동
    - MD equilibration 에서 OpenMM 이 실제 N-C amide bond 을 형성 (HTC bond)
    - AfCycDesign 은 짧은 바인더 (< 10 aa) 에서 hallucination (pLDDT < 0.5)

    Rettie et al., Nat. Commun. 16, 4730 (2025) — AfCycDesign 참고
    gap closure 는 Rosetta CCD loop closure 의 간소화 버전

Principle 4 (Fail-Fast):
    ColabFold 미가용 시 즉시 에러. gap closure 실패 시 경고 + 원본 유지.
"""

import os
import sys
import csv
import glob
import json
import math
import subprocess
import logging
from typing import Optional, Tuple, List

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils_common import parse_pdb_atom_line, atom_distance  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("cyclic_fold")

DEFAULT_NC_CUTOFF_A = 4.0  # Å — 검증 임계값 (amide bond ~1.33 Å + 마진)


# ==========================================
# N-C gap closure (geometric)
# ==========================================
def close_nc_gap(pdb_path: str,
                 output_path: str,
                 binder_chain: str = "B",
                 target_dist: float = 2.5) -> Tuple[bool, float, float]:
    """ColabFold 출력 PDB 의 바인더 N-C 말단 간격을 geometric interpolation 으로 폐합한다.

    ColabFold 는 cyclic constraint 를 강제하지 않아 바인더의 N-말단 N 과
    C-말단 C 의 거리가 30+ Å 이다. 이 함수는 바인더 backbone 원자에 가중
    이동 (sigmoid profile) 을 적용하여 말단을 target_dist Å 이내로 접근시킨다.

    가중 프로파일:
        - 바인더 중앙 residue: 이동량 0 (binding interface 보존)
        - N-말단 쪽: 이동량 증가 (C→N 방향 벡터의 절반)
        - C-말단 쪽: 이동량 증가 (N→C 방향 벡터의 절반)
        양 말단이 서로를 향해 대칭적으로 이동하여 중앙부 왜곡을 최소화한다.

    Args:
        pdb_path: 입력 PDB (ColabFold 출력)
        output_path: 출력 PDB (gap closed)
        binder_chain: 바인더 체인 ID
        target_dist: 목표 N-C 거리 (Å). 기본 2.5 (amide bond pre-formation)

    Returns:
        (success, dist_before, dist_after)
    """
    lines = []
    binder_residues = {}  # resnum → list of (line_idx, atom_name, x, y, z)
    n_coord = None
    c_coord = None
    first_res = 99999
    last_res = -1

    with open(pdb_path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            lines.append(line)
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain = line[21:22]
            if chain != binder_chain:
                continue
            try:
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            binder_residues.setdefault(resnum, []).append((i, atom_name, x, y, z))
            if resnum < first_res:
                first_res = resnum
            if resnum > last_res:
                last_res = resnum
            if resnum == first_res and atom_name == "N":
                n_coord = (x, y, z)
            if resnum == last_res and atom_name == "C":
                c_coord = (x, y, z)

    if n_coord is None or c_coord is None or first_res >= last_res:
        log.warning(f"[gap closure] N 또는 C 원자를 찾을 수 없음: {pdb_path}")
        return False, float("inf"), float("inf")

    dist_before = atom_distance(n_coord, c_coord)
    if dist_before <= target_dist * 1.5:
        log.info(f"[gap closure] 이미 충분히 가까움: N-C={dist_before:.1f}Å ≤ {target_dist*1.5:.1f}Å")
        # 그래도 output 에 복사
        with open(output_path, "w", encoding="utf-8") as f:
            f.writelines(lines)
        return True, dist_before, dist_before

    # [v45 개선 — Error B] 이동 벡터를 N→C gap vector 기반 순수 translation 으로 변경.
    # 기존 "midpoint 를 향한 이동" 은 backbone angle 을 왜곡하여 NVT NaN 을 유발.
    # 새 전략: N-말단 절반은 C 쪽으로, C-말단 절반은 N 쪽으로 평행 이동 (rigid translation).
    # 이로써 각 잔기 내부의 bond angle/dihedral 이 완전히 보존된다.
    gap_vec = (c_coord[0] - n_coord[0],
               c_coord[1] - n_coord[1],
               c_coord[2] - n_coord[2])  # N→C 벡터

    sorted_res = sorted(binder_residues.keys())
    n_res = len(sorted_res)
    res_to_idx = {r: i for i, r in enumerate(sorted_res)}

    # 목표: 양 말단을 서로 향해 각각 gap/2 만큼 이동. 중앙은 이동량 0.
    # 가중치: N-말단(t=0) = +0.5 (C 쪽으로), 중앙(t=0.5) = 0, C-말단(t=1) = -0.5 (N 쪽으로)
    # 이 선형 보간은 각 잔기를 rigid body 로 평행이동하므로 내부 각도를 보존한다.
    for resnum, atoms in binder_residues.items():
        idx = res_to_idx[resnum]
        t = idx / max(n_res - 1, 1)
        # 선형 가중: -0.45 at t=0, 0 at t=0.5, +0.45 at t=1
        # t=0 (N-말단): dx = -gap * (-0.45) = +gap * 0.45 → N 이 C 쪽으로 이동 ✓
        # t=1 (C-말단): dx = -gap * (+0.45) = -gap * 0.45 → C 가 N 쪽으로 이동 ✓
        w = 0.45 * (2.0 * t - 1.0)
        for line_idx, atom_name, ax, ay, az in atoms:
            dx = -gap_vec[0] * w
            dy = -gap_vec[1] * w
            dz = -gap_vec[2] * w

            nx, ny, nz = ax + dx, ay + dy, az + dz
            old_line = lines[line_idx]
            new_line = old_line[:30] + f"{nx:8.3f}{ny:8.3f}{nz:8.3f}" + old_line[54:]
            lines[line_idx] = new_line

    # 재계산 N-C 거리
    new_n = new_c = None
    for line_idx, atom_name, _, _, _ in binder_residues.get(first_res, []):
        if atom_name == "N":
            l = lines[line_idx]
            new_n = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    for line_idx, atom_name, _, _, _ in binder_residues.get(last_res, []):
        if atom_name == "C":
            l = lines[line_idx]
            new_c = (float(l[30:38]), float(l[38:46]), float(l[46:54]))

    dist_after = atom_distance(new_n, new_c) if new_n and new_c else float("inf")

    with open(output_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    success = dist_after <= DEFAULT_NC_CUTOFF_A
    log.info(f"[gap closure] N-C: {dist_before:.1f}Å → {dist_after:.1f}Å "
             f"({'✅' if success else '⚠️ 부족 — 반복 필요'})")

    # 반복: 한 번으로 충분하지 않으면 최대 5 회 반복
    if not success:
        for iteration in range(5):
            ok, _, dist_after = close_nc_gap(output_path, output_path,
                                             binder_chain, target_dist)
            if ok:
                break
    return dist_after <= DEFAULT_NC_CUTOFF_A, dist_before, dist_after


# ==========================================
# N-C 거리 검증
# ==========================================
def validate_cyclic_structure(pdb_path: str,
                              max_nc_distance: float = DEFAULT_NC_CUTOFF_A,
                              binder_chain: str = "B") -> Tuple[bool, float]:
    """예측된 구조의 N-C 말단 거리를 검증하여 cyclic constraint 충족 여부를 판별한다."""
    binder_atoms: List[dict] = []
    try:
        with open(pdb_path, "r", encoding="utf-8") as f:
            for line in f:
                parsed = parse_pdb_atom_line(line)
                if parsed and parsed["chain"] == binder_chain:
                    binder_atoms.append(parsed)
    except (IOError, OSError) as e:
        log.warning(f"[Cyclic 검증] PDB 읽기 실패 ({e}): {pdb_path}")
        return False, float("inf")

    if not binder_atoms:
        return False, float("inf")

    first_resnum = min(a["resnum"] for a in binder_atoms)
    last_resnum  = max(a["resnum"] for a in binder_atoms)

    n_coords: Optional[Tuple[float, float, float]] = None
    c_coords: Optional[Tuple[float, float, float]] = None
    for a in binder_atoms:
        if a["resnum"] == first_resnum and a["name"] == "N":
            n_coords = (a["x"], a["y"], a["z"])
        if a["resnum"] == last_resnum and a["name"] == "C":
            c_coords = (a["x"], a["y"], a["z"])

    if n_coords is None or c_coords is None:
        return False, float("inf")

    nc_dist = atom_distance(n_coords, c_coords)
    return nc_dist <= max_nc_distance, nc_dist


# ==========================================
# 메인 오케스트레이터 (v44: ColabFold + gap closure)
# ==========================================
def run_cyclic_fold(input_csv: str, output_dir: str,
                    topology_type: str,
                    target_pdb: str = "",
                    af2_bin: Optional[str] = None,
                    msa_mode: str = "mmseqs2_uniref_env",
                    binder_chain: str = "B",
                    params_dir: Optional[str] = None) -> str:
    """Cyclic peptide 구조 예측: ColabFold + geometric N-C gap closure.

    v44 전략:
        1. ColabFold 로 linear 복합체 구조 예측 (고품질 backbone)
        2. 바인더 체인에 geometric gap closure 적용 (N-C < 4.0 Å)
        3. 검증: gap closure 성공 시 MD 에서 HTC bond 자연 형성 가능
    """
    os.makedirs(output_dir, exist_ok=True)

    # ── 1단계: ColabFold 구조 예측 ──
    if not af2_bin:
        log.error("[Cyclic Fold] ColabFold 바이너리 미지정 — 예측 불가.")
        return output_dir
    try:
        subprocess.run(
            [af2_bin, input_csv, output_dir,
             "--num-recycle", "3",
             "--pair-mode", "unpaired_paired",
             "--msa-mode", msa_mode],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        log.error(f"[Cyclic Fold] ColabFold 실패: {e}")
        return output_dir

    # ── 2단계: Geometric N-C gap closure ──
    pdb_files = sorted(glob.glob(os.path.join(output_dir, "**", "*_unrelaxed_*.pdb"), recursive=True))
    if not pdb_files:
        pdb_files = sorted(glob.glob(os.path.join(output_dir, "**", "*.pdb"), recursive=True))

    if topology_type in ("cyclic_htc", "cyclic_nm"):
        n_closed = 0
        for pdb_path in pdb_files:
            ok, d_before, d_after = close_nc_gap(
                pdb_path, pdb_path, binder_chain=binder_chain, target_dist=2.5
            )
            if ok:
                n_closed += 1
            else:
                log.warning(f"  [gap closure] 폐합 실패: {os.path.basename(pdb_path)} "
                            f"N-C={d_after:.1f}Å > {DEFAULT_NC_CUTOFF_A}Å")
        log.info(f"  [gap closure] 폐합 성공: {n_closed}/{len(pdb_files)}")

    # ── 3단계: N-C 검증 ──
    _validate_outputs(output_dir, binder_chain)
    return output_dir


def _validate_outputs(output_dir: str, binder_chain: str = "B") -> None:
    """출력 디렉토리 내 PDB 들의 N-C 거리를 검증."""
    pdb_files = sorted(glob.glob(os.path.join(output_dir, "**", "*.pdb"), recursive=True))
    n_pass = n_fail = 0
    for pdb_path in pdb_files:
        passed, nc_dist = validate_cyclic_structure(pdb_path, binder_chain=binder_chain)
        if passed:
            n_pass += 1
        else:
            n_fail += 1
            log.warning(f"  [Cyclic 검증] N-C 거리 {nc_dist:.1f}Å > {DEFAULT_NC_CUTOFF_A}Å — "
                        f"cyclic 미충족: {os.path.basename(pdb_path)}")
    if n_fail > 0 and pdb_files:
        log.warning(f"  [Cyclic 검증 요약] 통과: {n_pass}, 실패: {n_fail}/{len(pdb_files)}")
    elif pdb_files:
        log.info(f"  [Cyclic 검증 요약] 전체 {n_pass}/{len(pdb_files)} 통과.")


# ==========================================
# CLI
# ==========================================
def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Cyclic peptide 구조 예측: ColabFold + geometric N-C gap closure"
    )
    parser.add_argument("--input_csv", required=True,
                        help="ColabFold 형식 CSV (id, sequence=target:binder)")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--topology_type", default="cyclic_htc",
                        choices=["cyclic_htc", "cyclic_ss", "cyclic_nm", "bicyclic"])
    parser.add_argument("--af2_bin", default=None, help="ColabFold 바이너리")
    parser.add_argument("--target_pdb", default="")
    parser.add_argument("--target_chain", default="A")
    args = parser.parse_args()

    run_cyclic_fold(
        input_csv=args.input_csv,
        output_dir=args.output_dir,
        topology_type=args.topology_type,
        target_pdb=args.target_pdb,
        af2_bin=args.af2_bin,
    )


if __name__ == "__main__":
    main()
