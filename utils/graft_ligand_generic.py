"""
graft_ligand_generic.py
-----------------------
Grafts HETATM ligand records from a reference complex onto a freshly-prepared
target protein by CA-based superposition.

[v35 AUDIT] Principle 7: chain-aware. Previously the CA superposition hardcoded
`chainid 0` for BOTH target and reference, silently producing a nonsensical
alignment whenever the target protein chain was not the first chain parsed
by mdtraj (e.g., multi-chain receptors, or reference files where a ligand
chain appeared before the protein chain). --target_chain and --ref_chain
expose the chain letter used for alignment; defaults preserve the historical
behavior so existing UPDD invocations remain valid.
"""

import os
import sys
import argparse
import mdtraj as md

# [Refactor] 공통 helper로 위임. utils 디렉토리를 sys.path에 등록하여
# CLI 모드(`python utils/graft_ligand_generic.py`)에서도 동일하게 import 가능.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils_common import resolve_chainid_by_letter as _resolve_chainid_common  # noqa: E402


def _resolve_chainid_by_letter(topology, chain_letter):
    """[Backwards-Compat Wrapper] utils_common.resolve_chainid_by_letter 위임."""
    return _resolve_chainid_common(topology, chain_letter)


def _select_protein_ca(traj, chain_letter, label):
    # Prefer an exact chain-letter match; if the letter cannot be resolved or
    # yields no CA atoms, fall back to the legacy "first chain" heuristic so
    # pre-existing runs continue to work without regression.
    chainid = _resolve_chainid_by_letter(traj.topology, chain_letter)
    if chainid is not None:
        sel = traj.topology.select(f"chainid {chainid} and name CA")
        if len(sel) > 0:
            print(f"  [{label}] 체인 '{chain_letter}' → chainid {chainid} (CA {len(sel)}개)")
            return sel
        print(f"  [{label}] 체인 '{chain_letter}' 해석 성공했으나 CA 원자 없음 — chainid 0 폴백")
    sel = traj.topology.select("chainid 0 and name CA")
    print(f"  [{label}] 기본 chainid 0 사용 (CA {len(sel)}개)")
    return sel


def graft_ligand(target_pdb, ref_pdb, output_pdb, target_chain="", ref_chain=""):
    print(f"🧬 Grafting from {ref_pdb} to {target_pdb}...")

    t_target = md.load(target_pdb)
    t_ref    = md.load(ref_pdb)

    ca_target = _select_protein_ca(t_target, target_chain, "Target")
    ca_ref    = _select_protein_ca(t_ref,    ref_chain,    "Reference")

    if len(ca_target) == 0 or len(ca_ref) == 0:
        raise RuntimeError(
            f"Grafting 불가: CA 원자 집합이 비어 있습니다 (target={len(ca_target)}, ref={len(ca_ref)}). "
            "chain 식별자가 PDB 파일의 실제 체인과 일치하는지 확인해 주십시오."
        )

    # [Quality] CA 원자 수가 비대칭이면 정렬에 사용되는 잔기가 잘리므로
    # superposition의 RMSD가 왜곡될 수 있다. 운영자에게 명시적으로 경고한다.
    if len(ca_target) != len(ca_ref):
        diff = abs(len(ca_target) - len(ca_ref))
        print(f"  [!] 경고: target ({len(ca_target)}) vs reference ({len(ca_ref)}) "
              f"CA 원자 수 불일치 (차이 {diff}개) — 짧은 쪽 길이로 절단하여 정렬합니다.")

    min_len = min(len(ca_target), len(ca_ref))

    t_target.superpose(
        t_ref,
        atom_indices=ca_target[:min_len],
        ref_atom_indices=ca_ref[:min_len]
    )

    # [v25 FIX retained] temp 파일을 output과 같은 디렉토리에 생성 (권한 오류 방지)
    temp_pdb = os.path.join(os.path.dirname(os.path.abspath(output_pdb)), "temp_aligned.pdb")
    t_target.save(temp_pdb)

    # [v35 AUDIT] Principle 7/v27: encoding="utf-8" for all text I/O.
    with open(ref_pdb, 'r', encoding="utf-8") as f:
        ligands = [l for l in f if l.startswith("HETATM") and "HOH" not in l]
    with open(temp_pdb, 'r', encoding="utf-8") as f:
        target_lines = [l for l in f if l.startswith("ATOM") or l.startswith("TER")]

    if not ligands:
        print("⚠️ [Warning] Reference PDB에 리간드(HETATM)가 없습니다. Grafting 없이 저장합니다.")

    with open(output_pdb, 'w', encoding="utf-8") as f:
        f.writelines(target_lines + ligands + ["END\n"])

    os.remove(temp_pdb)
    print(f"✅ Saved: {output_pdb}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True, help="타겟 단백질 PDB 경로")
    parser.add_argument("--ref",    required=True, help="리간드를 포함한 Reference PDB 경로")
    parser.add_argument("--output", required=True, help="출력 PDB 경로")
    parser.add_argument("--target_chain", default="", help="CA 정렬에 사용할 타겟 체인 문자 (비우면 첫 번째 체인)")
    parser.add_argument("--ref_chain",    default="", help="CA 정렬에 사용할 레퍼런스 체인 문자 (비우면 첫 번째 체인)")
    args = parser.parse_args()
    graft_ligand(args.target, args.ref, args.output, target_chain=args.target_chain, ref_chain=args.ref_chain)
