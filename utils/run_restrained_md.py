#!/usr/bin/env python
"""
run_restrained_md.py
--------------------
ncAA 치환된 구조를 OpenMM으로 Position-Restrained MD 수행한다.
- ncAA 잔기 백본(N/CA/C)은 강하게 구속하고 측쇄·용매·단백질만 적응시킨다.
- CUDA → OpenCL → CPU 순으로 플랫폼을 자동 선택한다.
- 사용 환경: qmmm (conda)

[Architecture Summary — v4 ~ v28]
v4–v9에서 토폴로지 패치/말단 보강/Cyclic 결합 탐색 및 PDBFixer 통합을 점진적으로
완성하였고, v10에서 "Topology First Protocol"(PDBFixer → 고리/말단 처리 → XML 패치 →
addHydrogens) 파이프라인 순서를 확정하였다. v11–v13은 WT/ncAA Rename 범위 정밀화,
다중 ncAA Hard Fail, Atom-index 기반 Restraint 변경, 병렬 XML Race condition 방지
등을 통해 Strict Production 등급으로 안정화하였다. v14–v17은 --ncaa none(WT) 경로
정상화, |q|>2.0 RESP Fail-Fast, Backbone 기반 ncAA 식별, --target_resid/--seed/
--dispersion CLI, 폭발 시 단계별 Partial PDB 저장 등 운영 신뢰성을 강화하였다.
v18–v23은 NetworkX Anchored Graph Isomorphism, 2-Stage Matching(Topology→거리 기반
Fallback), 추정 엣지 Provenance 로깅, HTC/SS Cyclization 트랜잭션화 및 메타데이터
반환을 도입하였다. v24–v28은 HTC/SS 시맨틱 분리, --graph_policy strict 기본값,
PDBFixer NMA→NME 위장(Masking), 모든 파일 입출력 encoding="utf-8" 강제, 파일명
해시 태그 12자리 확장 및 병렬 안정성 강화를 마무리한 단계이다. 단계별 상세 변경
내역은 UPDATE.md 의 각 버전 섹션을 참조한다.
"""

import os
import sys
import json
import glob
import argparse
import re
import shutil
import hashlib
import xml.etree.ElementTree as ET
from collections import defaultdict
import numpy as np
import networkx as nx

import openmm as mm
import openmm.unit as unit
from openmm import app
from openmm.app import (
    PDBFile, ForceField, Modeller, Topology, element,
    PME, HBonds, Simulation, StateDataReporter, DCDReporter
)
from pdbfixer import PDBFixer

# ==========================================
# 모듈 상수 및 유틸리티
# ==========================================
STD_AA_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "HID", "HIE", "HIP"
}

SOLVENT_ION_NAMES = {
    "HOH", "WAT", "NA", "CL", "K", "MG", "CA", "ZN"
}

def get_best_platform(preferred: str = "CUDA"):
    """[v4 H-1] OpenMM 플랫폼 자동 선택 — CUDA (mixed precision) 우선.

    Args:
        preferred: 우선 시도할 플랫폼명. 기본은 CUDA. CLI --platform 인자로 override 가능.

    Returns:
        (Platform, platform_name). 모든 선택 실패 시 (None, "Reference").

    폴백 경고:
        CUDA → OpenCL → CPU 순서로 시도하되, CUDA 실패 시 MD 가 100x 이상 느려지므로
        운영자에게 명시적 경고를 출력한다 (hardware / driver 문제 조기 포착).
    """
    # preferred 가 지정되면 그 플랫폼을 우선 시도하고, 실패 시 표준 폴백 체인
    search_order = []
    if preferred and preferred not in ("", "auto"):
        search_order.append(preferred)
    for name in ["CUDA", "OpenCL", "CPU"]:
        if name not in search_order:
            search_order.append(name)

    for i, name in enumerate(search_order):
        try:
            plat = mm.Platform.getPlatformByName(name)
            speed = plat.getSpeed()
            if i == 0:
                print(f"  [Platform] {name} 선택 (속도 지수 {speed:.0f}x)")
            else:
                fallback_msg = (
                    f"  [Platform] {name} 폴백 선택 (속도 지수 {speed:.0f}x) — "
                    f"상위 후보 실패. GPU 가속이 필요하다면 CUDA 드라이버·cudatoolkit 설치를 확인하세요."
                )
                print(fallback_msg)
                if name == "CPU":
                    print(f"  ⚠️  [Platform WARNING] CPU 플랫폼으로 실행됩니다 — MD 속도가 50-100배 느려집니다!")
            return plat, name
        except Exception as _plat_err:
            print(f"  [Platform] {name} 사용 불가 ({_plat_err.__class__.__name__})")
            continue
    print(f"  ⚠️  [Platform WARNING] CUDA/OpenCL/CPU 모두 사용 불가 → Reference 플랫폼 폴백")
    return None, "Reference"


def get_platform_properties(platform_name):
    """[v4 H-1] 플랫폼별 properties. CUDA 는 mixed precision + DeviceIndex 0 기본."""
    if platform_name == "CUDA":
        return {"DeviceIndex": "0", "Precision": "mixed"}
    elif platform_name == "OpenCL":
        return {"DeviceIndex": "0", "Precision": "mixed"}
    return {}

def is_peptide_like_atomset(atom_names):
    return {"N", "CA", "C"}.issubset(atom_names)

def make_unique_tag(pdb_path, topology, ncaa_code):
    raw_str = f"{os.path.abspath(pdb_path)}_{topology}_{ncaa_code}"
    return hashlib.md5(raw_str.encode()).hexdigest()[:12]

# ==========================================
# RESP XML 로드 및 검증
# ==========================================
def load_ncaa_manifest(manifest_path):
    if not manifest_path or not os.path.exists(manifest_path):
        raise RuntimeError(f"ncAA 모드: 파라미터 매니페스트를 찾을 수 없습니다 ({manifest_path}). 즉시 중단합니다.")
    try:
        # [v35 AUDIT] Principle 7/v27: encoding="utf-8" enforced on manifest I/O to
        # match the producer side in parameterize_ncaa.write_manifest. Without this,
        # a 한글 path segment in manifest_path + a non-UTF-8 locale raised
        # UnicodeDecodeError and v31's load-bearing fail-fast was never reached.
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)
    except Exception as e:
        raise RuntimeError(f"매니페스트 파싱 실패: {e}")

    xml_path = manifest.get("xml_path")
    if not xml_path or not os.path.exists(xml_path):
        raise RuntimeError(f"매니페스트 내 XML 경로 누락 또는 파일 없음: {xml_path}")

    xml_res_name = manifest.get("xml_resname")
    if not xml_res_name:
        raise RuntimeError("매니페스트 내 XML 잔기명(xml_resname) 정보 누락")

    # [v36 CRITICAL FIX-2] Hydrogen Definitions 파일 경로. 부재해도 downstream 은
    # 가능한 모든 경로를 시도하지만, 로딩이 누락되면 OpenMM addHydrogens 가
    # ncAA 잔기에 H 를 추가하지 못해 createSystem 매칭이 실패한다. Principle 4
    # (Fail-Fast): 부재 시 명시 경고를 출력하여 Silent-Fallback 을 방지.
    hydrogens_path = manifest.get("hydrogens_path", "")
    if hydrogens_path and not os.path.exists(hydrogens_path):
        print(f"  [경고] manifest 내 hydrogens_path 가 가리키는 파일이 존재하지 않습니다: {hydrogens_path}")
        hydrogens_path = ""
    if not hydrogens_path:
        print(f"  [경고] manifest 에 hydrogens_path 가 없습니다 — 구버전 캐시일 수 있습니다. 재파라미터화를 권장합니다.")

    print(f"  ncAA XML 파라미터 로드 (Manifest 기반): [{os.path.basename(xml_path)}] (XML 내부 잔기명: {xml_res_name})")
    if hydrogens_path:
        print(f"  ncAA Hydrogen Definitions 로드: [{os.path.basename(hydrogens_path)}]")
    return [xml_path], xml_res_name, hydrogens_path

def find_peptide_like_nonstandard(lines, binder_chain):
    non_std = []
    by_res = defaultdict(set)
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
        if line[21] != binder_chain: continue
        rname = line[17:20].strip()
        if rname in STD_AA_NAMES or rname in SOLVENT_ION_NAMES: continue
        resid = line[22:26].strip()
        atom_name = line[12:16].strip()
        by_res[(rname, resid)].add(atom_name)
        
    for (rname, resid), atoms in by_res.items():
        if is_peptide_like_atomset(atoms):
            non_std.append((rname, resid))
    return non_std

# ==========================================
# Graph Isomorphism Helper Functions
# ==========================================
def prepare_signatures(G):
    for n in G.nodes():
        G.nodes[n]["degree"] = G.degree[n]
        nbr_elems = sorted(G.nodes[m].get("element") or "?" for m in G.neighbors(n))
        G.nodes[n]["nbr_sig"] = tuple(nbr_elems)

    anchors_dict = {d["atom_name"]: n for n, d in G.nodes(data=True) if d.get("anchor")}
    for anchor_type in ["N", "CA", "C"]:
        anchor_node = anchors_dict.get(anchor_type)
        if anchor_node is not None and anchor_node in G:
            lengths = nx.single_source_shortest_path_length(G, anchor_node)
            for n in G.nodes():
                G.nodes[n][f"dist_{anchor_type}"] = lengths.get(n, 999)
        else:
            for n in G.nodes():
                G.nodes[n][f"dist_{anchor_type}"] = 999

def node_match(td, xd):
    if td.get("element") != xd.get("element"): return False
    if td.get("degree") != xd.get("degree"): return False
    if td.get("anchor") != xd.get("anchor"): return False
    if td.get("n_external") != xd.get("n_external"): return False
    if td.get("nbr_sig") != xd.get("nbr_sig"): return False
    if td.get("dist_N") != xd.get("dist_N"): return False
    if td.get("dist_CA") != xd.get("dist_CA"): return False
    if td.get("dist_C") != xd.get("dist_C"): return False
    return True

def graph_signature_table(G):
    rows = []
    for n, d in G.nodes(data=True):
        rows.append({
            "node": d.get("atom_name", str(n)),
            "element": d.get("element"),
            "degree": d.get("degree"),
            "n_ext": d.get("n_external"),
            "nbr_sig": d.get("nbr_sig"),
            "d_N": d.get("dist_N"),
            "d_CA": d.get("dist_CA"),
            "d_C": d.get("dist_C")
        })
    return rows

def explain_graph_mismatch(G_tmpl, G_tgt):
    tmpl_rows = graph_signature_table(G_tmpl)
    tgt_rows = graph_signature_table(G_tgt)
    res = "--- Template Signatures (Top 15) ---\n"
    for r in tmpl_rows[:15]: res += f"  {r}\n"
    res += "--- Target Signatures (Top 15) ---\n"
    for r in tgt_rows[:15]: res += f"  {r}\n"
    return res

def try_graph_match(G_tmpl, G_tgt, required_anchor_pairs):
    matcher = nx.algorithms.isomorphism.GraphMatcher(G_tmpl, G_tgt, node_match=node_match)
    valid_mappings = []
    for mapping in matcher.isomorphisms_iter():
        if all(mapping.get(t) == x for t, x in required_anchor_pairs.items()):
            valid_mappings.append(mapping)
    return valid_mappings

# ==========================================
# Explicit Bond Graph Re-injection
# ==========================================
def inject_xml_bonds(topology, n_xmls, xml_res_name):
    """[v38 Multi-Site] PDBFixer 가 삭제한 ncAA 내부 결합을 XML 템플릿으로부터
    복구한다. 동일 ncAA 가 chain 내 다중 위치에 배치된 경우 (AUTO_SASA:N 의
    정상 동작) **모든 잔기 인스턴스**에 대해 동일 결합 그래프를 주입한다.

    변경 이력:
        v37 이전: ``next(...)`` 로 첫 잔기 1 개에만 bond 주입. 다중 위치 치환
            시 두 번째 이후 잔기의 내부 결합이 PDBFixer 에 의해 삭제된 채로
            ForceField createSystem 에 전달되어 template 매칭 실패.
        v38: 모든 ``r.name == xml_res_name`` 잔기에 대해 반복 적용. 단일 위치
            치환에서는 기존 동작과 동등 (idempotent).
    """
    if not n_xmls or not xml_res_name: return
    import xml.etree.ElementTree as ET

    target_residues = [r for r in topology.residues() if r.name == xml_res_name]
    if not target_residues: return

    tree = ET.parse(n_xmls[0])
    xml_res_node = next((r for r in tree.getroot().iter("Residue") if r.get("name") == xml_res_name), None)
    if not xml_res_node: return

    existing_bonds = {frozenset([b.atom1.index, b.atom2.index]) for b in topology.bonds()}
    xml_bond_pairs = [
        (b.get("atomName1", "").strip(), b.get("atomName2", "").strip())
        for b in xml_res_node.findall("Bond")
    ]

    added_internal_total = 0
    for target_res in target_residues:
        atom_by_name = {a.name: a for a in target_res.atoms()}
        added_for_this_res = 0
        for a1, a2 in xml_bond_pairs:
            if a1 in atom_by_name and a2 in atom_by_name:
                oa1, oa2 = atom_by_name[a1], atom_by_name[a2]
                pair = frozenset([oa1.index, oa2.index])
                if pair not in existing_bonds:
                    topology.addBond(oa1, oa2)
                    existing_bonds.add(pair)
                    added_for_this_res += 1
        added_internal_total += added_for_this_res

    if added_internal_total > 0:
        print(
            f"  [XML Bond Injection] PDBFixer가 삭제한 '{xml_res_name}' 내부 결합 "
            f"{added_internal_total}개 완전 복구 완료 ({len(target_residues)}개 잔기 인스턴스)"
        )

# ==========================================
# Symmetry-Equivalent Mapping Resolver
# ==========================================
def _mapping_ff_signature(mapping, tmpl_name_to_type):
    """단일 graph isomorphism 매핑을 {target_atom_index: xml_atom_type} 서명으로 변환한다.

    Graph isomorphism 이 반환하는 mapping 은 {template_atom_name: target_atom_index}
    형식이다. AMBER ForceField 는 residue template 의 Atom 원소별 ``type``
    속성으로 질량/전하/LJ 파라미터를 주입하므로, 매핑의 물리적 의미는 "어떤
    target 원자에 어떤 XML atom type 이 할당되는가" 로 환원된다.
    """
    return {
        tgt_idx: tmpl_name_to_type.get(tmpl_name)
        for tmpl_name, tgt_idx in mapping.items()
    }

def _are_mappings_forcefield_equivalent(mappings, tmpl_name_to_type):
    """다수 매핑이 동일한 force field 파라미터를 산출하는지 검증한다.

    대칭 동등 원자(homotopic atoms — 예: 메틸기의 HB1/HB2/HB3 세 수소) 의
    순열은 동일 XML atom type 을 공유하므로, 어느 매핑을 선택해도 AMBER
    residue template 이 target 원자에 부여하는 atom type 이 완전히 동일하다.
    본 함수는 각 매핑의 (target_idx → xml_type) 서명을 dict 비교로 검증하여
    이 조건을 확인한다. 휴리스틱(3! 거듭제곱 패턴 등) 에 의존하지 않고
    물리적 파라미터 동등성을 직접 증명한다 (Principle 4: Fail-Fast over
    Silent-Fallback — 동등하지 않으면 즉시 False 반환).

    Returns:
        bool: 모든 매핑이 동일 FF 서명을 가지면 True.
    """
    if len(mappings) <= 1:
        return True
    reference = _mapping_ff_signature(mappings[0], tmpl_name_to_type)
    for m in mappings[1:]:
        if _mapping_ff_signature(m, tmpl_name_to_type) != reference:
            return False
    return True

# ==========================================
# Universal XML Patch
# ==========================================
def apply_universal_xml_patch(mdl, n_xmls, xml_res_name, output_dir, patched_basename, is_cyclic_stage=False, graph_policy="strict"):
    """[v38 Multi-Site] ncAA 잔기의 XML 템플릿을 topology 의 실제 원자 그래프에
    맞게 rename 하여 단일 patched XML 을 생성한다.

    다중 위치 치환 처리 (CLAUDE.md v6):
        ForceField 의 Residue 템플릿은 잔기 **이름** 기준으로 매칭되므로, 동일
        ncAA 가 chain 내 여러 위치에 배치되어도 patched XML 은 한 개만 필요
        하다. 본 함수는 다음을 수행한다:

        1. Topology 내 xml_res_name 에 해당하는 모든 잔기를 수집.
        2. 모든 잔기가 동일한 atom-name 집합을 가지는지 검증 (ncaa_mutate /
           parameterize_ncaa 파이프라인의 불변량 — Principle 7: 데이터 계약).
        3. 대표 잔기 (첫 번째) 로 graph isomorphism 매핑 수행.
        4. 매핑된 atom-name / external bond 카운트를 사용해 XML 재작성.

        변경 이력:
            v37 이전: ``len(target_residues) > 1 → raise`` — 다중 치환 시 즉시 실패.
            v38: 이종 ncAA 가 아닌 동일 ncAA 다중 위치 (AUTO_SASA:N 정상 경로)
                 를 허용. 원자명 불일치는 여전히 Fail-Fast 로 차단 (Principle 4).
    """
    from openmm import unit as _u

    target_residues = [res for res in mdl.topology.residues() if res.name == xml_res_name]
    if not target_residues: raise RuntimeError(f"Topology 내에 '{xml_res_name}' 잔기를 찾을 수 없습니다.")

    # [v38 Multi-Site] 다중 잔기 atom-name 일관성 검증. 본 체크는 parameterize_ncaa
    # 가 단일 XML 을 모든 인스턴스에 공유하는 데이터 계약 (동일 ncaa_code → 동일
    # 원자명) 을 강제한다. 불일치 시 runtime 에서 즉시 실패하여 silent
    # corruption 을 예방한다.
    if len(target_residues) > 1:
        ref_atom_names = frozenset(a.name.strip() for a in target_residues[0].atoms())
        for idx, res in enumerate(target_residues[1:], start=2):
            this_names = frozenset(a.name.strip() for a in res.atoms())
            if this_names != ref_atom_names:
                raise RuntimeError(
                    f"동일 ncAA '{xml_res_name}' 의 다중 인스턴스 간 원자명 불일치 감지: "
                    f"#{1} (ID={target_residues[0].id}): {sorted(ref_atom_names)} vs "
                    f"#{idx} (ID={res.id}): {sorted(this_names)}. "
                    f"ncaa_mutate/parameterize_ncaa 파이프라인의 데이터 계약 위반."
                )
        print(
            f"  [Graph Patch] 다중 위치 치환 감지: '{xml_res_name}' × {len(target_residues)}개소 "
            f"(ID: {[r.id for r in target_residues]}). 원자명 일관성 검증 통과."
        )
    target_res = target_residues[0]

    tree = ET.parse(n_xmls[0])
    root = tree.getroot()
    xml_res_node = next((r for r in root.iter("Residue") if r.get("name") == xml_res_name), None)
    if xml_res_node is None: raise RuntimeError(f"XML 내에 '{xml_res_name}' 정의가 없습니다.")

    # ── 1. 템플릿 그래프 구성 ──
    atomtypes = {t.get("name"): {"element": t.get("element"), "class": t.get("class")} for t in root.findall("AtomTypes/Type")}
    G_tmpl = nx.Graph()
    for atom in xml_res_node.findall("Atom"):
        aname = atom.get("name").strip()
        atype = atom.get("type")
        element_sym = atomtypes.get(atype, {}).get("element")
        G_tmpl.add_node(aname, atom_name=aname, element=element_sym, anchor=(aname in {"N", "CA", "C"}), n_external=0)
    
    for bond in xml_res_node.findall("Bond"):
        G_tmpl.add_edge(bond.get("atomName1").strip(), bond.get("atomName2").strip())
        
    for ext_bond in xml_res_node.findall("ExternalBond"):
        aname = ext_bond.get("atomName").strip()
        if aname in G_tmpl: G_tmpl.nodes[aname]["n_external"] += 1

    # ── 2. 타겟 그래프 구성 (1차 Strict / 2차 Fallback) ──
    target_atoms = list(target_res.atoms())
    target_atom_indices = {a.index for a in target_atoms}
    
    def build_target_graph(use_distance_fallback=False):
        G_tgt = nx.Graph()
        inferred_edges = []
        for atom in target_atoms:
            elem_sym = atom.element.symbol if atom.element else None
            G_tgt.add_node(atom.index, atom_name=atom.name.strip(), element=elem_sym, anchor=(atom.name.strip() in {"N", "CA", "C"}), n_external=0)
        
        existing_bonds = set()
        for bond in mdl.topology.bonds():
            i, j = bond.atom1.index, bond.atom2.index
            existing_bonds.add(frozenset([i, j]))
            if i in target_atom_indices and j in target_atom_indices:
                G_tgt.add_edge(i, j, inferred=False)
            elif i in target_atom_indices and j not in target_atom_indices:
                G_tgt.nodes[i]["n_external"] += 1
            elif j in target_atom_indices and i not in target_atom_indices:
                G_tgt.nodes[j]["n_external"] += 1
        
        if use_distance_fallback:
            COVALENT_RADII = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05, "Cl": 1.02}
            def get_radius(sym): return COVALENT_RADII.get(sym.capitalize(), 0.8)
            pos = np.array(mdl.positions.value_in_unit(_u.angstroms))
            for i in range(len(target_atoms)):
                for j in range(i + 1, len(target_atoms)):
                    a1, a2 = target_atoms[i], target_atoms[j]
                    if frozenset([a1.index, a2.index]) not in existing_bonds:
                        if a1.element is None or a2.element is None: continue
                        dist = np.linalg.norm(pos[a1.index] - pos[a2.index])
                        thresh = (get_radius(a1.element.symbol) + get_radius(a2.element.symbol)) * 1.15
                        if dist <= thresh:
                            G_tgt.add_edge(a1.index, a2.index, inferred=True)
                            inferred_edges.append({
                                "a1_name": a1.name, "a1_idx": a1.index,
                                "a2_name": a2.name, "a2_idx": a2.index,
                                "dist_A": round(float(dist * 10), 3), "thresh_A": round(float(thresh * 10), 3)
                            })
        return G_tgt, inferred_edges

    # ── 3. 앵커 검증 ──
    def get_req_anchors(G_tgt):
        tmpl_anchors_dict = {d["atom_name"]: n for n, d in G_tmpl.nodes(data=True) if d.get("anchor")}
        tgt_anchors_dict = {d["atom_name"]: n for n, d in G_tgt.nodes(data=True) if d.get("anchor")}
        req_pairs = {}
        for t_name in ["N", "CA", "C"]:
            t_node = tmpl_anchors_dict.get(t_name)
            x_node = tgt_anchors_dict.get(t_name)
            if t_node and not x_node: raise RuntimeError(f"타겟 잔기에서 필수 백본 앵커 '{t_name}' 누락")
            if t_node and x_node: req_pairs[t_node] = x_node
        return req_pairs

    # ── 4. 매핑 시도 (Strict -> Relaxed Fallback) ──
    G_tgt_strict, _ = build_target_graph(use_distance_fallback=False)
    prepare_signatures(G_tmpl)
    prepare_signatures(G_tgt_strict)
    req_anchors = get_req_anchors(G_tgt_strict)
    
    valid_mappings = try_graph_match(G_tmpl, G_tgt_strict, req_anchors)
    match_mode = "Strict Topology (100% 신뢰)"
    inferred_edges_final = []

    if not valid_mappings:
        if graph_policy == "strict":
            mismatch_dump = explain_graph_mismatch(G_tmpl, G_tgt_strict)
            raise RuntimeError(f"그래프 매칭 완전 실패 (Strict 모드에서 거리 기반 보완은 엄격히 금지됩니다).\n\n{mismatch_dump}")
        print("  [Graph Match] Strict 매칭 실패. 3D 거리 기반 보완 그래프(Relaxed Fallback)로 재시도합니다.")
        G_tgt_relaxed, inferred_edges_final = build_target_graph(use_distance_fallback=True)
        prepare_signatures(G_tgt_relaxed)
        valid_mappings = try_graph_match(G_tmpl, G_tgt_relaxed, req_anchors)
        match_mode = "Relaxed Distance Fallback (거리 추정 엣지 의존)"
        G_tgt = G_tgt_relaxed
    else:
        G_tgt = G_tgt_strict

    if not valid_mappings:
        mismatch_dump = explain_graph_mismatch(G_tmpl, G_tgt)
        raise RuntimeError(f"그래프 매칭 완전 실패 (Strict 및 Relaxed 모두 동형성 확보 불가).\n\n{mismatch_dump}")

    # [v39] 대칭 동등 매핑(symmetry-equivalent automorphism) 허용
    # ─────────────────────────────────────────────────────────────
    # Graph isomorphism 이 다수 매핑을 반환하는 경우, 그 대부분은 homotopic
    # atoms (메틸기 수소 HM1/HM2/HM3, HB1/HB2/HB3 등) 의 순열로 발생한다.
    # 이들 원자는 동일한 AMBER atom type 과 부분 전하를 공유하므로 어느
    # 매핑을 선택해도 force field 파라미터와 MD 시뮬레이션 결과가 수치적으로
    # 완전히 동일하다 (automorphism under force field equivalence). 따라서
    # 유효 매핑 간 (target_idx → xml_type) 서명이 일치하면 첫 번째 매핑을
    # 채택하고, 불일치하는 비대칭 모호성만 Fail-Fast 로 차단한다.
    # Principle 4 (Fail-Fast over Silent-Fallback) 유지: 진정한 모호성은
    # 여전히 _FAILED_XML_PATCH.log 와 함께 파이프라인을 중단시킨다.
    if len(valid_mappings) > 1:
        tmpl_name_to_type = {
            atom.get("name").strip(): atom.get("type")
            for atom in xml_res_node.findall("Atom")
        }
        if _are_mappings_forcefield_equivalent(valid_mappings, tmpl_name_to_type):
            n = len(valid_mappings)
            import math as _math
            # 메틸기 1개 = 3! = 6, 2개 = 36, 3개 = 216 의 전형적 패턴
            log6 = _math.log(n) / _math.log(6) if n > 0 else 0.0
            est_methyl = round(log6) if abs(log6 - round(log6)) < 1e-6 else None
            methyl_note = f"(메틸기 ≈ {est_methyl}개 추정) " if est_methyl else ""
            print(
                f"  [Graph Patch] 대칭 동등 매핑 {n}개 감지 {methyl_note}— "
                f"force field 파라미터 동등성 검증 통과. 첫 번째 automorphism 채택."
            )
        else:
            raise RuntimeError(
                f"비대칭 모호성: {len(valid_mappings)}개의 유효 매핑이 서로 다른 "
                f"force field 파라미터(XML atom type) 를 산출합니다. 원자명 지정 "
                f"또는 템플릿 수정이 필요합니다. (Fail-Fast)"
            )

    if inferred_edges_final and graph_policy == "strict":
        raise RuntimeError("Relaxed Graph Fallback 금지 정책 위반 (strict 모드에서는 거리 기반 추정 엣지를 절대 허용하지 않습니다.)")

    best_mapping = valid_mappings[0]

    # ── 5. XML 적용 ──
    rename_map = {G_tmpl.nodes[t]["atom_name"]: G_tgt.nodes[x]["atom_name"] for t, x in best_mapping.items()}

    for atom in xml_res_node.findall("Atom"):
        old_name = atom.get("name").strip()
        if old_name in rename_map: atom.set("name", rename_map[old_name])

    for b in xml_res_node.findall("Bond"):
        a1, a2 = b.get("atomName1").strip(), b.get("atomName2").strip()
        b.set("atomName1", rename_map.get(a1, a1))
        b.set("atomName2", rename_map.get(a2, a2))

    for eb in xml_res_node.findall("ExternalBond"): xml_res_node.remove(eb)
    
    ext_bonds = set(G_tgt.nodes[n]["atom_name"] for n, d in G_tgt.nodes(data=True) if d.get("n_external") > 0)
    for ext in ext_bonds:
        ET.SubElement(xml_res_node, "ExternalBond", {"atomName": ext})

    patched_xml = os.path.join(output_dir, f"{patched_basename}_universal.xml")
    tree.write(patched_xml)
    
    patch_info = {"match_mode": match_mode, "inferred_edges": inferred_edges_final}
    stage = "고리형" if is_cyclic_stage else "선형"
    print(f"  [Graph Patch ({stage})] {xml_res_name} 동형성 매핑 성공! [{match_mode}]")
    return [patched_xml], patch_info

# ==========================================
# Cyclic 펩타이드 토폴로지
# ==========================================
def atom_ref(atom):
    return {"chain": atom.residue.chain.id, "resid": atom.residue.id, "resname": atom.residue.name, "atom": atom.name}

def detect_head_to_tail_pair(modeller, binder_chain="B", max_htc_dist_nm=0.40):
    # [v45 FIX — Error A] 임계값 0.30 nm (3.0 Å) → 0.40 nm (4.0 Å).
    # gap closure 출력의 N-C 거리는 1.3-4.0 Å 범위이며, 3-4 Å 은
    # bond 형성 후 OpenMM minimization 으로 자연 수렴하는 거리이다.
    # Engh & Huber (1991): amide bond 1.33 Å, strain at 4 Å ≈ 5 kcal/mol.
    chains = [c for c in modeller.topology.chains() if c.id == binder_chain]
    if not chains: return None
    residues = list(chains[0].residues())
    if not residues: return None

    n_term, c_term = residues[0], residues[-1]
    from openmm import unit as _u
    pos = list(modeller.positions)

    n_atom = next((a for a in n_term.atoms() if a.name == "N"), None)
    if not n_atom: n_atom = next((a for a in n_term.atoms() if a.element and a.element.symbol == "N"), None)

    if not n_atom: return None
    np_pos = np.array(pos[n_atom.index].value_in_unit(_u.nanometers))

    c_atom = next((a for a in c_term.atoms() if a.name == "C"), None)
    if not c_atom:
        carbonyl_candidates = []
        for a in c_term.atoms():
            if a.element and a.element.symbol == "C":
                cp = np.array(pos[a.index].value_in_unit(_u.nanometers))
                dist = np.linalg.norm(cp - np_pos)
                if dist < max_htc_dist_nm:
                    carbonyl_candidates.append((dist, a))
        if carbonyl_candidates:
            carbonyl_candidates.sort(key=lambda x: x[0])
            c_atom = carbonyl_candidates[0][1]

    if n_atom and c_atom:
        dist = np.linalg.norm(np.array(pos[c_atom.index].value_in_unit(_u.nanometers)) - np_pos)
        if dist <= max_htc_dist_nm:
            return (c_atom, n_atom)
        else:
            print(f"  [Warning] HTC 예측 모호함: N-말단 N과 C-말단 C 거리가 {dist*10:.2f} Å 입니다.")
    return None

def detect_disulfide_pair(modeller, binder_chain="B", max_sg_dist_nm=0.24):
    chains = [c for c in modeller.topology.chains() if c.id == binder_chain]
    if not chains: return None
    residues = list(chains[0].residues())
    sg_atoms = [a for r in residues if r.name in ("CYS", "CYX") for a in r.atoms() if a.name == "SG"]
    
    if len(sg_atoms) < 2: return None
    
    from openmm import unit as _u
    pos = list(modeller.positions)
    existing = {frozenset([b.atom1.index, b.atom2.index]) for b in modeller.topology.bonds()}
    
    candidates = []
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            a1, a2 = sg_atoms[i], sg_atoms[j]
            p1 = np.array(pos[a1.index].value_in_unit(_u.nanometers))
            p2 = np.array(pos[a2.index].value_in_unit(_u.nanometers))
            dist = np.linalg.norm(p1 - p2)
            already_bonded = frozenset([a1.index, a2.index]) in existing
            if already_bonded or dist <= max_sg_dist_nm:
                candidates.append((dist, a1, a2, already_bonded))
                
    if not candidates: return None
    
    candidates.sort(key=lambda x: (not x[3], x[0]))
    best = candidates[0]
    
    if len(candidates) > 1 and not best[3]: 
        if abs(candidates[1][0] - candidates[0][0]) < 0.02:
            raise RuntimeError("Disulfide pair가 둘 이상으로 모호합니다 (거리 차이 < 0.02 nm). 명시적 수동 규칙이 필요합니다.")
            
    return (best[1], best[2])

def commit_bond_if_missing(topology, atom1, atom2):
    existing = {frozenset([b.atom1.index, b.atom2.index]) for b in topology.bonds()}
    pair = frozenset([atom1.index, atom2.index])
    if pair not in existing:
        topology.addBond(atom1, atom2)
        return True
    return False

def _trim_termini_for_htc(modeller, binder_chain):
    atoms_to_delete = []
    for chain in modeller.topology.chains():
        if chain.id == binder_chain:
            res_list = list(chain.residues())
            if res_list:
                n_term, c_term = res_list[0], res_list[-1]
                for a in n_term.atoms():
                    if a.name in ["H2", "H3"]: atoms_to_delete.append(a)
                    elif a.name == "H1": a.name = "H"
                for a in c_term.atoms():
                    if a.name in ["OXT", "HXT"]: atoms_to_delete.append(a)
    if atoms_to_delete:
        modeller.delete(atoms_to_delete)
        print(f"  [cyclic_htc] HTC 결합 적용 전 체인 {binder_chain} 말단 잉여 원자(OXT/H2 등) {len(atoms_to_delete)}개 절제 완료")
    return modeller

def add_cyclic_bond_to_topology(modeller, topology_type, binder_chain="B"):
    htc_pair = detect_head_to_tail_pair(modeller, binder_chain)
    ss_pair = detect_disulfide_pair(modeller, binder_chain)
    
    meta = {"is_cyclized": False, "htc_formed": False, "ss_formed": False, "htc_pair": None, "ss_pair": None, "topology_type": topology_type}

    if topology_type in ("cyclic_htc", "cyclic_nm"):
        if not htc_pair: return meta, modeller
        modeller = _trim_termini_for_htc(modeller, binder_chain)
        htc_pair = detect_head_to_tail_pair(modeller, binder_chain)
        if not htc_pair: raise RuntimeError("말단 절제 후 HTC pair 재탐색에 실패했습니다. 구조 결함 의심.")
            
        is_new = commit_bond_if_missing(modeller.topology, *htc_pair)
        status = "추가 완료" if is_new else "이미 존재함"
        print(f"  [{topology_type}] Peptide bond {status}: {htc_pair[0].residue.name}{htc_pair[0].residue.id}.{htc_pair[0].name} → {htc_pair[1].residue.name}{htc_pair[1].residue.id}.{htc_pair[1].name}")
        meta.update({"is_cyclized": True, "htc_formed": True, "htc_pair": (atom_ref(htc_pair[0]), atom_ref(htc_pair[1]))})
        return meta, modeller

    if topology_type == "cyclic_ss":
        if not ss_pair: return meta, modeller
        is_new = commit_bond_if_missing(modeller.topology, *ss_pair)
        status = "추가 완료" if is_new else "이미 존재함"
        print(f"  [cyclic_ss] Disulfide bond {status}: {ss_pair[0].residue.name}{ss_pair[0].residue.id}.SG ↔ {ss_pair[1].residue.name}{ss_pair[1].residue.id}.SG")
        meta.update({"is_cyclized": True, "ss_formed": True, "ss_pair": (atom_ref(ss_pair[0]), atom_ref(ss_pair[1]))})
        return meta, modeller

    if topology_type == "bicyclic":
        if not (htc_pair and ss_pair): return meta, modeller
        modeller = _trim_termini_for_htc(modeller, binder_chain)
        htc_pair = detect_head_to_tail_pair(modeller, binder_chain)
        ss_pair = detect_disulfide_pair(modeller, binder_chain)
        if not (htc_pair and ss_pair): raise RuntimeError("말단 절제 후 Bicyclic pair 재탐색에 실패했습니다.")
            
        commit_bond_if_missing(modeller.topology, *htc_pair)
        commit_bond_if_missing(modeller.topology, *ss_pair)
        print(f"  [bicyclic] HTC 및 SS 결합 트랜잭션 반영 완료")
        meta.update({"is_cyclized": True, "htc_formed": True, "ss_formed": True, "htc_pair": (atom_ref(htc_pair[0]), atom_ref(htc_pair[1])), "ss_pair": (atom_ref(ss_pair[0]), atom_ref(ss_pair[1]))})
        return meta, modeller

    return meta, modeller

# ==========================================
# Peptide Bond 보정 강화
# ==========================================
def is_carbonyl_carbon(residue, c_atom, positions_nm, max_co_nm=0.14):
    if c_atom.element is None or c_atom.element.symbol != "C": return False
    c_pos = np.array(positions_nm[c_atom.index].value_in_unit(unit.nanometers))
    for a in residue.atoms():
        if a.index == c_atom.index: continue
        if a.element and a.element.symbol == "O":
            o_pos = np.array(positions_nm[a.index].value_in_unit(unit.nanometers))
            if np.linalg.norm(c_pos - o_pos) <= max_co_nm:
                return True
    return False

def add_missing_peptide_bonds_safe(modeller, binder_chain="B", max_cn_distance_nm=0.20):
    from openmm import unit as _u
    pos = list(modeller.positions)
    added = 0
    existing = {frozenset([b.atom1.index, b.atom2.index]) for b in modeller.topology.bonds()}

    for chain in modeller.topology.chains():
        if chain.id != binder_chain: continue
        residues = list(chain.residues())
        for r1, r2 in zip(residues[:-1], residues[1:]):
            r1_names = {a.name for a in r1.atoms()}
            r2_names = {a.name for a in r2.atoms()}
            if not (is_peptide_like_atomset(r1_names) and is_peptide_like_atomset(r2_names)): continue

            c_atom = next((a for a in r1.atoms() if a.name == "C"), None)
            n_atom = next((a for a in r2.atoms() if a.name == "N"), None)
            if c_atom is None or n_atom is None: continue
            
            if not is_carbonyl_carbon(r1, c_atom, pos): continue

            pair = frozenset([c_atom.index, n_atom.index])
            if pair in existing: continue

            c_pos = np.array(pos[c_atom.index].value_in_unit(_u.nanometers))
            n_pos = np.array(pos[n_atom.index].value_in_unit(_u.nanometers))
            dist  = np.linalg.norm(c_pos - n_pos)

            if dist <= max_cn_distance_nm:
                modeller.topology.addBond(c_atom, n_atom)
                existing.add(pair)
                print(f"  [peptide bond 보정] 체인 {chain.id}: {r1.name}{r1.id}.C → {r2.name}{r2.id}.N ({dist*10:.3f} Å)")
                added += 1
    return added

def add_missing_oxt(modeller, cyclic_chains=None):
    from openmm.app import Topology as _Topo
    from openmm.app import element as _elem
    from openmm import unit as _u

    if cyclic_chains is None: cyclic_chains = set()

    old_topo = modeller.topology
    old_pos  = list(modeller.positions)

    need_oxt = set()
    for chain in old_topo.chains():
        if chain.id in cyclic_chains: continue
        residues = list(chain.residues())
        if not residues: continue
        c_term = residues[-1]
        names  = {a.name for a in c_term.atoms()}
        if "OXT" in names or "OC2" in names: continue
        atom_names = {a.name: a for a in c_term.atoms()}
        if all(k in atom_names for k in ("C", "CA", "O")):
            need_oxt.add(id(c_term))

    if not need_oxt: return modeller, 0

    new_topo = _Topo()
    if old_topo.getPeriodicBoxVectors() is not None:
        new_topo.setPeriodicBoxVectors(old_topo.getPeriodicBoxVectors())

    old_to_new, new_pos, added = {}, [], 0

    for chain in old_topo.chains():
        new_chain = new_topo.addChain(chain.id)
        for res in chain.residues():
            new_res = new_topo.addResidue(res.name, new_chain, res.id, res.insertionCode)
            for atom in res.atoms():
                new_atom = new_topo.addAtom(atom.name, atom.element, new_res)
                old_to_new[atom.index] = new_atom
                new_pos.append(old_pos[atom.index])

            if id(res) in need_oxt:
                amap   = {a.name: a for a in res.atoms()}
                c_pos  = np.array(old_pos[amap["C"].index].value_in_unit(_u.nanometers))
                ca_pos = np.array(old_pos[amap["CA"].index].value_in_unit(_u.nanometers))
                o_pos  = np.array(old_pos[amap["O"].index].value_in_unit(_u.nanometers))

                v_ca     = ca_pos - c_pos
                norm_ca  = v_ca / np.linalg.norm(v_ca)
                v_o      = o_pos  - c_pos
                norm_o   = v_o  / np.linalg.norm(v_o)
                oxt_dir  = -(norm_ca + norm_o)
                oxt_dir /= np.linalg.norm(oxt_dir)
                c_o_dist = np.linalg.norm(v_o)
                oxt_pos  = c_pos + oxt_dir * c_o_dist

                oxt_atom = new_topo.addAtom("OXT", _elem.oxygen, new_res)
                new_topo.addBond(old_to_new[amap["C"].index], oxt_atom)
                new_pos.append(oxt_pos * _u.nanometers)
                print(f"  [OXT 추가] 체인 {chain.id}: {res.name}{res.id} (C-OXT = {c_o_dist*10:.3f} Å)")
                added += 1

    for bond in old_topo.bonds():
        a1 = old_to_new.get(bond.atom1.index)
        a2 = old_to_new.get(bond.atom2.index)
        if a1 and a2:
            new_topo.addBond(a1, a2)

    return Modeller(new_topo, new_pos), added

def renumber_protein_residues(pdb_path):
    renum_lines, last_key_by_chain, current_id_by_chain = [], {}, {}
    out_path = pdb_path.replace(".pdb", "_renum.pdb")

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                chain   = line[21]
                res_key = (line[22:26], line[26], line[17:20])
                if last_key_by_chain.get(chain) != res_key:
                    current_id_by_chain[chain] = current_id_by_chain.get(chain, 0) + 1
                    last_key_by_chain[chain]   = res_key
                new_resid = f"{current_id_by_chain[chain]:>4d}"
                renum_lines.append(line[:22] + new_resid + line[26:])
            else:
                renum_lines.append(line)

    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(renum_lines)
        if not renum_lines or not renum_lines[-1].strip().startswith("END"): f.write("END\n")
    return out_path

# ==========================================
# Position Restraints
# ==========================================
def add_position_restraints(system, positions, topology, xml_res_name, k_restrain=1000.0, junction_resids=None):
    """[v46] ncAA 잔기의 백본 (N, CA, C) 에 조화 위치 구속을 적용한다.

    동일 ncAA 가 chain 내 여러 위치에 배치된 경우 (AUTO_SASA:N 정상 동작),
    **모든 인스턴스의 백본 원자** 를 단일 CustomExternalForce 에 등록하여
    restraint 그룹을 구성한다. 잔기 간 k 값 차별화는 현재 필요하지 않으므로
    글로벌 파라미터 ``k`` 하나로 동기 스케줄링된다 (Minimization → NVT 워밍업
    → Production 전환 시 일괄 완화).

    변경 이력:
        v37 이전: ``len(target_residues) > 1 → raise`` — 다중 치환 불가.
        v38: 모든 target 잔기의 백본을 루프로 append. 단일 치환 시 기존 동작과
             동등 (idempotent).
        v46: junction_resids 파라미터 추가 — cyclic HTC junction 잔기(첫/마지막
             2잔기)는 restraint 에서 제외하여 angle/dihedral 최적화 허용.
    """
    restraint = mm.CustomExternalForce("k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addGlobalParameter("k", k_restrain)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    target_residues = [res for res in topology.residues() if res.name == xml_res_name]
    if not target_residues: raise RuntimeError("Position Restraint 대상 잔기를 찾지 못했습니다.")

    n_restrained = 0
    BACKBONE_ATOMS = {"N", "CA", "C"}
    per_residue_counts = []

    n_junction_skipped = 0
    for target_res in target_residues:
        rc = 0
        for atom in target_res.atoms():
            if atom.name.strip() not in BACKBONE_ATOMS: continue
            if junction_resids and atom.residue.index in junction_resids:
                n_junction_skipped += 1
                continue
            pos = positions[atom.index]
            restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])
            n_restrained += 1
            rc += 1
        per_residue_counts.append((target_res.id, rc))

    if n_restrained > 0:
        system.addForce(restraint)
        junction_msg = f", junction 제외={n_junction_skipped}" if n_junction_skipped > 0 else ""
        if len(target_residues) == 1:
            print(f"  Position Restraint 준비: {n_restrained}개 백본(N, CA, C) 원자 앵커링 (k_init={k_restrain} kJ/mol/nm²{junction_msg})")
        else:
            print(
                f"  Position Restraint 준비: {xml_res_name} × {len(target_residues)}개소에서 "
                f"총 {n_restrained}개 백본(N, CA, C) 원자 앵커링 "
                f"(per-res: {per_residue_counts}, k_init={k_restrain} kJ/mol/nm²{junction_msg})"
            )
    elif n_junction_skipped > 0:
        print(f"  [v46] Position Restraint: ncAA가 junction 잔기에만 위치하여 restraint 0개 (junction 제외={n_junction_skipped})")
    else:
        raise RuntimeError("Position Restraint 대상 백본 원자(N, CA, C)를 찾지 못했습니다.")

    return system, n_restrained

# ==========================================
# CONECT Appender
# ==========================================
def find_serial(lines, ref):
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")): continue
        if line[21].strip() != ref["chain"]: continue
        if line[22:26].strip() != str(ref["resid"]).strip(): continue
        if line[17:20].strip() != ref["resname"]: continue
        if line[12:16].strip() != ref["atom"]: continue
        return int(line[6:11].strip())
    return None

def append_conect_records(out_pdb, cycle_meta):
    with open(out_pdb, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    conect_lines = []
    
    if cycle_meta.get("htc_formed") and cycle_meta.get("htc_pair"):
        a, b = cycle_meta["htc_pair"]
        s1, s2 = find_serial(lines, a), find_serial(lines, b)
        if s1 and s2:
            conect_lines.append(f"CONECT{s1:5d}{s2:5d}\n")
            conect_lines.append(f"CONECT{s2:5d}{s1:5d}\n")
            print(f"  [Visual] HTC CONECT 추가 완료 (N={s2}, C={s1})")
            
    if cycle_meta.get("ss_formed") and cycle_meta.get("ss_pair"):
        a, b = cycle_meta["ss_pair"]
        s1, s2 = find_serial(lines, a), find_serial(lines, b)
        if s1 and s2:
            conect_lines.append(f"CONECT{s1:5d}{s2:5d}\n")
            conect_lines.append(f"CONECT{s2:5d}{s1:5d}\n")
            print(f"  [Visual] SS CONECT 추가 완료 (SG1={s1}, SG2={s2})")
            
    if conect_lines:
        with open(out_pdb, "w", encoding="utf-8") as f:
            for line in lines:
                if line.strip() != "END": f.write(line)
            for line in conect_lines: f.write(line)
            f.write("END\n")

# ==========================================
# 메인 MD 함수
# ==========================================
# ==========================================
# [v3 2-1] run_md 4 분할 — 헬퍼 함수 (R-1: run_md 정의보다 위에 배치)
# ==========================================
class MDStageError(Exception):
    """run_md 의 4 단계 (prepare → build → equilibrate → produce) 중 하나가
    회복 불가한 사유로 중단될 때 발생시키는 예외.

    Attributes:
        status      : run_md 가 호출자에게 반환할 상태 문자열
                      ("FAIL" | "SUCCESS_WITH_WARNING" 등).
        log_suffix  : 진단용 로그 파일 접미사 (예: "_FAILED_FF.log").
        partial_stage: Partial PDB 저장이 필요한 단계명 ("Minimization" /
                       "NVT" / "Production MD") 또는 None.
        partial_suffix: Partial PDB 파일명 접미사 (단계명과 함께 사용).
    """
    def __init__(self, message, status="FAIL", log_suffix="",
                 partial_stage=None, partial_suffix=""):
        super().__init__(message)
        self.status = status
        self.log_suffix = log_suffix
        self.partial_stage = partial_stage
        self.partial_suffix = partial_suffix


class _MDContext:
    """run_md 4 단계 헬퍼들이 공유하는 가변 상태 컨테이너.

    [v3 2-1] 단일 거대한 함수에 흩어져 있던 지역 변수를 명시적 객체로
    응집하여, 각 단계가 다른 단계의 산물에 명확히 의존하도록 만든다.
    이는 향후 단위 테스트(특정 단계만 격리 실행) 의 토대가 된다.
    """
    def __init__(self, pdb_path, params_manifest, output_dir, n_steps,
                 topology_type, ncaa_label, ncaa_code, hdd_path,
                 binder_chain, seed, disp_corr, target_resid, graph_policy,
                 platform_pref="CUDA", dt_fs=2.0):
        # ── 입력 인자 (불변) ──
        self.pdb_path = pdb_path
        self.params_manifest = params_manifest
        self.output_dir = output_dir
        self.n_steps = n_steps
        self.topology_type = topology_type
        self.ncaa_label = ncaa_label
        self.ncaa_code = ncaa_code
        self.hdd_path = hdd_path
        self.binder_chain = binder_chain
        self.seed = seed
        self.disp_corr = disp_corr
        self.target_resid = target_resid
        self.graph_policy = graph_policy
        self.platform_pref = platform_pref
        self.dt_fs = dt_fs

        # ── 파일 경로 ──
        self.basename = os.path.basename(pdb_path).replace(".pdb", "")
        self.out_dcd  = os.path.join(output_dir, f"{self.basename}_restrained.dcd")
        self.out_pdb  = os.path.join(output_dir, f"{self.basename}_final.pdb")
        self.out_log  = os.path.join(output_dir, f"{self.basename}_md.log")
        self.std_pdb  = None  # renumber 후 prepare_structure 가 채움

        # ── 단계간 산물 ──
        self.modeller = None
        self.ncaa_xmls = []
        self.xml_res_name = None
        self.n_xmls = None
        # [v36] manifest 로부터 추출되는 OpenMM Hydrogens.xml 경로. addHydrogens
        # 단계 전에 Modeller.loadHydrogenDefinitions 로 등록되어야 함.
        self.hydrogens_path = ""
        self.cycle_meta = {}
        self.cyclic_chains = set()
        self.system = None
        self.integrator = None
        self.simulation = None
        self.n_restrained = 0

        # ── 상태 플래그 ──
        self.warnings_found = False
        self.production_partial = False

    def log_status(self, reason, status_type, suffix=""):
        print(f"  [{status_type}] {reason}")
        if status_type == "WARNING":
            self.warnings_found = True
        if suffix:
            with open(self.out_log.replace("_md.log", suffix), "w", encoding="utf-8") as f:
                f.write(f"{reason}\n")

    def save_partial_pdb(self, stage_name, suffix):
        partial_pdb = os.path.join(self.output_dir, f"{self.basename}{suffix}")
        try:
            state = self.simulation.context.getState(getPositions=True)
            with open(partial_pdb, "w", encoding="utf-8") as f:
                PDBFile.writeFile(self.modeller.topology, state.getPositions(), f)
            print(f"  [!] {stage_name} 터지기 전 Partial PDB 저장 완료: {os.path.basename(partial_pdb)}")
        except Exception as save_err:
            print(f"  [!] 좌표계 붕괴로 인해 PDB 강제 저장 실패: {save_err}")


def _diagnose_template_mismatch(topology, xml_res_name, xml_path):
    """[v4 CLAUDE.md 방안 B] addHydrogens "No template found" 오류 발생 시 원자명 불일치 진단.

    OpenMM 이 XML 템플릿을 찾지 못하는 사유는 거의 항상 다음 중 하나이다:
        (a) PDB 내 잔기의 중원자 이름이 XML <Residue> 블록의 <Atom> 이름과 불일치
        (b) PDB 잔기가 원자를 덜/더 가지고 있음 (chain ID 가 잘못 배치되어 원자
            일부가 다른 유령 잔기로 분리된 v4 버그와 같은 현상)

    본 함수는 두 집합을 수집·비교하여 어느 쪽이 문제인지 즉시 식별할 수 있도록
    상세 진단을 stdout 으로 출력한다. 호출 비용이 낮으므로 실패 경로에서만 실행된다.
    """
    pdb_atoms = set()
    pdb_residue_count = 0
    for res in topology.residues():
        if res.name == xml_res_name:
            pdb_residue_count += 1
            pdb_atoms.update({atom.name for atom in res.atoms()})

    xml_atoms = set()
    if xml_path and os.path.exists(xml_path):
        try:
            tree = ET.parse(xml_path)
            for res_node in tree.getroot().iter("Residue"):
                if res_node.get("name") == xml_res_name:
                    xml_atoms = {a.get("name") for a in res_node.findall("Atom")}
                    break
        except ET.ParseError as _xml_err:
            print(f"  [진단] XML 파싱 실패: {xml_path} — {_xml_err}")

    print(f"\n  [진단] Template Mismatch 분석 ({xml_res_name}):")
    print(f"    Topology 내 '{xml_res_name}' 잔기 개수: {pdb_residue_count}")
    if pdb_residue_count == 0:
        print(f"    ⚠️  Topology 에 해당 잔기가 전혀 없음 → XML 로드 경로 재확인 필요")
    print(f"    PDB 원자 (합집합): {sorted(pdb_atoms)}")
    print(f"    XML 템플릿 원자:   {sorted(xml_atoms)}")

    pdb_heavy = {a for a in pdb_atoms if not a.startswith('H')}
    xml_heavy = {a for a in xml_atoms if not a.startswith('H')}
    print(f"    PDB 중원자: {sorted(pdb_heavy)}")
    print(f"    XML 중원자: {sorted(xml_heavy)}")

    only_pdb = pdb_atoms - xml_atoms
    only_xml = xml_atoms - pdb_atoms
    if only_pdb:
        print(f"    PDB 에만 있는 원자: {sorted(only_pdb)}")
    if only_xml:
        print(f"    XML 에만 있는 원자: {sorted(only_xml)}")

    if pdb_heavy != xml_heavy:
        print(f"    ⚠️  중원자 이름 불일치 — 하단 원인 중 하나를 의심:")
        print(f"         · ncaa_mutate._make_cm_atom_line 등의 PDB 포맷 문자열 col 오정렬")
        print(f"         · parameterize_ncaa._rename_xml_to_pdb_names 실패")
        print(f"         · PDB 내 잔기가 chain ID 누락으로 복수 잔기로 분리됨")
    else:
        print(f"    ✅ 중원자 일치 — 수소 원자 정의 또는 결합 그래프 문제 의심")


def prepare_structure(ctx: "_MDContext") -> None:
    """[v3 2-1] PDBFixer 정리, ncAA UNK 마스킹, Topology 재구성, XML 패치, 고리화, OXT 보강.

    수행 단계:
        1) PDB renumber → ctx.std_pdb
        2) WT/ncAA 모드 판별 및 manifest 로드
        3) ncAA 잔기 식별 + UNK 임시 위장 (PDBFixer 보호)
        4) PDBFixer 실행 + UNK → 본래 잔기명 복구
        5) Modeller 생성, inject_xml_bonds
        6) addHydrogens (선형 캡핑)
        7) apply_universal_xml_patch (Graph Isomorphism)
        8) add_missing_peptide_bonds_safe + add_cyclic_bond_to_topology
        9) add_missing_oxt

    상태 변경:
        ctx.modeller, ctx.ncaa_xmls, ctx.xml_res_name, ctx.n_xmls,
        ctx.cycle_meta, ctx.cyclic_chains 를 채운다.

    Raises: MDStageError on any stage failure.
    """
    print(f"\n  처리 중: {ctx.basename}")
    ctx.std_pdb = renumber_protein_residues(ctx.pdb_path)
    unique_tag = make_unique_tag(ctx.pdb_path, ctx.topology_type, ctx.ncaa_code)
    patched_basename = f"{ctx.basename}_{unique_tag}"

    # ── 1. 모드 판별 및 정밀 Rename 방어 (위장/Masking 포함) ──
    with open(ctx.std_pdb, "r", encoding="utf-8") as f:
        lines = f.readlines()

    target_orig_rname = None
    target_orig_resids: set = set()
    if ctx.ncaa_label.lower() == "none":
        print("  [WT Mode] ncAA 라벨이 'none'입니다. 야생형(Wild-type) 단백질 일반 MD를 수행합니다.")
        wt_offenders = find_peptide_like_nonstandard(lines, ctx.binder_chain)
        if wt_offenders:
            raise MDStageError(
                f"WT 모드이나, 입력 구조에 비표준 펩타이드 잔기가 존재합니다: {wt_offenders}. 자동 진행 불가.",
                log_suffix="_FAILED_WT_HAS_NONSTANDARD.log",
            )
    else:
        try:
            ctx.n_xmls, ctx.xml_res_name, ctx.hydrogens_path = load_ncaa_manifest(ctx.params_manifest)
            if ctx.ncaa_code and ctx.ncaa_code.upper() != ctx.xml_res_name.upper():
                raise MDStageError(
                    f"--ncaa_code 입력값('{ctx.ncaa_code}')이 XML 잔기명('{ctx.xml_res_name}')과 다릅니다. 명시적 검증 실패.",
                    log_suffix="_FAILED_NCAA_MISMATCH.log",
                )
        except RuntimeError as e:
            raise MDStageError(str(e), log_suffix="_FAILED_NO_XML.log") from e

        non_std_res_map = defaultdict(set)
        for line in lines:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21] == ctx.binder_chain:
                rname = line[17:20].strip()
                if rname not in STD_AA_NAMES and rname not in SOLVENT_ION_NAMES:
                    resid = line[22:26].strip()
                    atom_name = line[12:16].strip()
                    non_std_res_map[(rname, resid)].add(atom_name)

        ncaa_candidates = []
        if ctx.target_resid:
            matched = [(rname, resid) for (rname, resid) in non_std_res_map.keys() if resid == ctx.target_resid]
            if len(matched) == 0:
                raise MDStageError(
                    f"target_resid={ctx.target_resid} 에 해당하는 비표준 잔기가 없습니다.",
                    log_suffix="_FAILED_MISSING_NCAA.log",
                )
            if len(matched) > 1:
                raise MDStageError(
                    f"target_resid={ctx.target_resid} 에 해당하는 비표준 잔기가 여러 개입니다.",
                    log_suffix="_FAILED_MULTIPLE_HETATM.log",
                )
            rname, resid = matched[0]
            if not is_peptide_like_atomset(non_std_res_map[(rname, resid)]):
                raise MDStageError(
                    f"target_resid={ctx.target_resid} 는 비표준 잔기이지만 펩타이드 뼈대(N,CA,C)가 없습니다.",
                    log_suffix="_FAILED_NON_PEPTIDE.log",
                )
            ncaa_candidates.append((rname, resid))
        else:
            for (rname, resid), atoms in non_std_res_map.items():
                if is_peptide_like_atomset(atoms):
                    ncaa_candidates.append((rname, resid))

        if len(ncaa_candidates) == 0:
            raise MDStageError(
                f"체인 {ctx.binder_chain}에 펩타이드 뼈대(N,CA,C)를 가진 비표준 잔기가 없습니다. 구조 결함 의심.",
                log_suffix="_FAILED_MISSING_NCAA.log",
            )

        # [v38 Multi-Site CLAUDE.md v6] 동일 ncAA 다중 위치 치환은 허용하고,
        # 이종 ncAA 혼재만 차단한다. 원래 방어막 (v37 이전) 은 NMA × 2 와 같은
        # AUTO_SASA:N 의 정상 출력까지 차단하여 다중 치환 파이프라인 자체가
        # 불가능했다. 차단해야 할 진짜 오염 케이스 (예: NMA + TMS 혼재,
        # peptide-backbone 이 없는 ligand 혼입) 는 (a) 이종 잔기 종류 개수
        # 검사와 (b) 기존의 ``is_peptide_like_atomset`` 전제로 이미 구분된다.
        unique_ncaa_types = {rname for rname, _resid in ncaa_candidates}
        if len(unique_ncaa_types) > 1:
            raise MDStageError(
                f"체인 {ctx.binder_chain}에 서로 다른 종류의 비표준 잔기가 혼재합니다: "
                f"{sorted(unique_ncaa_types)}. 이종 ncAA 오염 방지를 위해 중단합니다.",
                log_suffix="_FAILED_MIXED_NCAA.log",
            )
        if len(ncaa_candidates) > 1:
            _ncaa_type = ncaa_candidates[0][0]
            _resid_list = [resid for _, resid in ncaa_candidates]
            print(
                f"  [Multi-Site] 동일 ncAA '{_ncaa_type}'의 다중 위치 치환 감지: "
                f"잔기 ID {_resid_list} ({len(ncaa_candidates)}개) — 허용 (AUTO_SASA:N 정상 경로)"
            )

        # [v38] target_orig_rname 은 단일 타입이 보장되므로 대표값으로 유지한다
        # (unique_ncaa_types 의 유일한 원소). target_orig_resids 는 set 으로
        # 일반화하여 UNK 마스킹 / 복구 단계에서 모든 인스턴스를 추적한다.
        target_orig_rname = ncaa_candidates[0][0]
        target_orig_resids = {resid for _, resid in ncaa_candidates}

        renamed_count = 0
        with open(ctx.std_pdb, "w", encoding="utf-8") as f:
            for line in lines:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21] == ctx.binder_chain:
                    if line[17:20].strip() == target_orig_rname and line[22:26].strip() in target_orig_resids:
                        line = line[:17] + "UNK" + line[20:]
                        renamed_count += 1
                f.write(line)
        print(
            f"  [Rename] 체인 {ctx.binder_chain}의 대상 잔기 '{target_orig_rname}' × "
            f"{len(target_orig_resids)}개소 (ID: {sorted(target_orig_resids)}) "
            f"원자 {renamed_count}개를 PDBFixer 보호용 'UNK'로 임시 위장 완료"
        )

    # ── PDBFixer + UNK 복구 ──
    print("  [MD] PDBFixer로 AF2 구조 정리 중...")
    fixer = PDBFixer(filename=ctx.std_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if ctx.ncaa_label.lower() != "none" and ctx.xml_res_name:
        # [v38 Multi-Site] 마스킹된 모든 resid 를 set 기준으로 복구한다. ID 기반
        # 정밀 매칭이 우선이며, 매칭 실패 시에만 UNK 전체 복구 폴백을 사용한다
        # (fallback 은 단일 ID 대응 실패 케이스의 잔존 안전망 — Principle 4).
        target_resid_strs = {str(r) for r in target_orig_resids}
        restored_ids = set()
        for res in fixer.topology.residues():
            if res.name == "UNK" and str(res.id) in target_resid_strs:
                res.name = ctx.xml_res_name
                restored_ids.add(str(res.id))
        missing_ids = target_resid_strs - restored_ids
        if missing_ids:
            # PDBFixer 가 residue id 를 재할당한 경우에 대비한 fallback. 남은
            # UNK 전체를 ctx.xml_res_name 으로 복구한다. 단, 이 경로는 운영자에게
            # 명시 경고를 남긴다.
            unk_remaining = [res for res in fixer.topology.residues() if res.name == "UNK"]
            for res in unk_remaining:
                res.name = ctx.xml_res_name
                restored_ids.add(str(res.id))
            print(
                f"  [Rename][WARN] 정밀 ID 매칭 실패 (미복구 ID: {sorted(missing_ids)}). "
                f"잔존 UNK {len(unk_remaining)}개를 '{ctx.xml_res_name}'으로 일괄 폴백 복구."
            )
        print(
            f"  [Rename] 임시 위장된 잔기 {len(restored_ids)}개소 "
            f"(ID: {sorted(restored_ids)})를 본래의 '{ctx.xml_res_name}'으로 최종 복구 완료"
        )

    with open(ctx.std_pdb, "w", encoding="utf-8") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    ff_files = ["amber14-all.xml", "amber14/tip3pfb.xml"]
    ctx.modeller = Modeller(fixer.topology, fixer.positions)

    # ── 1.5 ncAA 다중 방어 (v38 Multi-Site 완화) ──
    # [v38] 동일 ncAA 다중 위치 치환은 AUTO_SASA:N 의 정상 출력이므로 허용한다.
    # 0개 케이스만 Fail-Fast 로 차단 (구조 결함 / 이름 복구 실패). 이종 ncAA
    # 혼재는 위 단계의 unique_ncaa_types 검사가 이미 걸러냈다.
    if ctx.xml_res_name:
        target_residues = [r for r in ctx.modeller.topology.residues() if r.name == ctx.xml_res_name]
        if len(target_residues) == 0:
            raise MDStageError(
                f"구조 내에 '{ctx.xml_res_name}' 잔기가 존재하지 않습니다.",
                log_suffix="_FAILED_MISSING_NCAA.log",
            )
        _chain_ids = sorted({r.chain.id for r in target_residues})
        _res_ids = [r.id for r in target_residues]
        print(
            f"  [검증] 대상 ncAA 잔기 확정: 체인 {_chain_ids}의 {ctx.xml_res_name} "
            f"× {len(target_residues)}개소 (ID: {_res_ids})"
        )

    # ── 2. Explicit Bond Graph Re-injection (Topology First) ──
    if ctx.xml_res_name and ctx.n_xmls:
        inject_xml_bonds(ctx.modeller.topology, ctx.n_xmls, ctx.xml_res_name)

    # ── 3. addHydrogens 선형 캡핑 ──
    # [v4 CLAUDE.md 방안 B] ForceField 로드를 파일별로 명시 검증하여 XML 문제를 조기 포착.
    # [v36 CRITICAL FIX-2] ncAA Hydrogen Definitions 등록 — Modeller.addHydrogens
    # 는 _residueHydrogens 데이터베이스에 등록된 잔기명에 대해서만 H 를 추가한다.
    # 표준 아미노산만 기본 등록되어 있으므로, ncAA 가 manifest 에 포함된 경우
    # 반드시 Modeller.loadHydrogenDefinitions() 로 등록해야 한다. 미등록 시
    # NMA 는 H 가 추가되지 않은 채 통과되어 createSystem 의 템플릿 매칭이
    # heavy-only(6) vs heavy+H(14) mismatch 로 실패한다.
    if ctx.hydrogens_path and os.path.exists(ctx.hydrogens_path):
        try:
            Modeller.loadHydrogenDefinitions(ctx.hydrogens_path)
            print(f"  [Hydrogens] {os.path.basename(ctx.hydrogens_path)} 등록 완료 (Modeller._residueHydrogens 갱신)")
        except Exception as _hd_err:
            raise MDStageError(
                f"Hydrogen Definitions 로드 실패 ({ctx.hydrogens_path}): {_hd_err}",
                log_suffix="_FAILED_HYDROGENS_LOAD.log",
            ) from _hd_err
    print("  [MD] OpenMM addHydrogens 실행 (선형 구조 토폴로지 동기화)...")
    try:
        ff_temp = ForceField(*ff_files)
        if ctx.n_xmls:
            for xml_path in ctx.n_xmls:
                try:
                    ff_temp.loadFile(xml_path)
                    print(f"  [ForceField] ncAA XML 로드 성공: {os.path.basename(xml_path)}")
                except (mm.OpenMMException, ValueError, FileNotFoundError) as _load_err:
                    print(f"  [!] ncAA XML 로드 실패: {xml_path} — {_load_err}")
                    _diagnose_template_mismatch(ctx.modeller.topology, ctx.xml_res_name, xml_path)
                    raise MDStageError(
                        f"ncAA XML ForceField 로드 실패 ({os.path.basename(xml_path)}): {_load_err}",
                        log_suffix="_FAILED_XML_LOAD.log",
                    ) from _load_err
        ctx.modeller.addHydrogens(ff_temp)
        print("  [INFO] 수소 추가 완료 (선형 캡핑 및 내부 결합 매핑 성공)")
    except MDStageError:
        # ForceField 로드 단계에서 이미 raise 한 MDStageError 는 그대로 전파
        raise
    except (mm.OpenMMException, ValueError) as e:
        # [v4] "No template found" 계열 오류 발생 시 원자명 불일치 진단 실행
        _diagnose_template_mismatch(
            ctx.modeller.topology, ctx.xml_res_name,
            ctx.n_xmls[0] if ctx.n_xmls else None,
        )
        raise MDStageError(
            f"선형 수소 추가/토폴로지 동기화 OpenMM 오류: {e}",
            log_suffix="_FAILED_ADDH.log",
        ) from e
    except Exception as e:
        _diagnose_template_mismatch(
            ctx.modeller.topology, ctx.xml_res_name,
            ctx.n_xmls[0] if ctx.n_xmls else None,
        )
        raise MDStageError(
            f"선형 수소 추가/토폴로지 동기화 예상치 못한 오류 ({e.__class__.__name__}): {e}",
            log_suffix="_FAILED_ADDH.log",
        ) from e

    if ctx.xml_res_name and ctx.n_xmls:
        # [v38 Multi-Site] addHydrogens 후 모든 ncAA 잔기 인스턴스의 원자 수가
        # XML 템플릿과 일치하는지 검증. 하나라도 불일치면 protonation sync
        # 경고를 기록하며 잔기 ID 별로 차이를 보고한다.
        all_target_res = [r for r in ctx.modeller.topology.residues() if r.name == ctx.xml_res_name]
        if all_target_res:
            tree = ET.parse(ctx.n_xmls[0])
            xml_res_node = next((r for r in tree.getroot().iter("Residue") if r.get("name") == ctx.xml_res_name), None)
            if xml_res_node is not None:
                xml_atom_count = len(xml_res_node.findall("Atom"))
                mismatches = [
                    (r.id, len(list(r.atoms())))
                    for r in all_target_res
                    if len(list(r.atoms())) != xml_atom_count
                ]
                if mismatches:
                    ctx.log_status(
                        f"Modeller.addHydrogens() 복구 후에도 원자 수가 XML({xml_atom_count})과 일치하지 않는 "
                        f"잔기 인스턴스: {mismatches}. 강제로 Graph Matching으로 넘어가면 실패할 수 있습니다.",
                        "WARNING", "_WARNING_PROTONATION_SYNC.log",
                    )

    # ── 4. XML 패치 (Graph Isomorphism) ──
    if ctx.xml_res_name and ctx.n_xmls:
        try:
            ctx.ncaa_xmls, patch_info = apply_universal_xml_patch(
                ctx.modeller, ctx.n_xmls, ctx.xml_res_name, ctx.output_dir,
                patched_basename, is_cyclic_stage=False, graph_policy=ctx.graph_policy,
            )
            if patch_info.get("inferred_edges"):
                msg = f"Relaxed Graph Fallback 사용됨: 추정된 엣지 {len(patch_info['inferred_edges'])}개"
                for edge in patch_info["inferred_edges"]:
                    msg += f"\n    ↳ {edge['a1_name']} ↔ {edge['a2_name']} (dist: {edge['dist_A']}A, thresh: {edge['thresh_A']}A)"
                if ctx.graph_policy == "strict":
                    raise MDStageError(msg, log_suffix="_FAILED_RELAXED_GRAPH_MATCH.log")
                elif ctx.graph_policy == "warn":
                    ctx.log_status(msg, "WARNING", "_WARNING_RELAXED_GRAPH_MATCH.log")
                else:
                    print(f"  [INFO] {msg}")
        except RuntimeError as e:
            raise MDStageError(str(e), log_suffix="_FAILED_XML_PATCH.log") from e
    else:
        ctx.ncaa_xmls = []

    # ── 5. Peptide Bonds & Cyclization ──
    add_missing_peptide_bonds_safe(ctx.modeller, binder_chain=ctx.binder_chain)
    if ctx.topology_type not in ("linear", None, ""):
        try:
            ctx.cycle_meta, ctx.modeller = add_cyclic_bond_to_topology(
                ctx.modeller, ctx.topology_type, binder_chain=ctx.binder_chain
            )
        except RuntimeError as e:
            raise MDStageError(f"고리화 단계 실패: {e}", log_suffix="_FAILED_CYCLIZATION.log") from e
        if ctx.cycle_meta.get("is_cyclized"):
            if ctx.cycle_meta.get("htc_formed"):
                ctx.cyclic_chains.add(ctx.binder_chain)
        else:
            raise MDStageError(
                f"요청한 고리형 토폴로지('{ctx.topology_type}') 형성 실패. 과학적 위험을 방지하기 위해 중단합니다.",
                log_suffix="_FAILED_TOPOLOGY_UNMET.log",
            )

    ctx.modeller, _ = add_missing_oxt(ctx.modeller, cyclic_chains=ctx.cyclic_chains)


def build_openmm_system(ctx: "_MDContext") -> None:
    """[v3 2-1] ForceField 로드, 용매화, System 생성, Position Restraint, Dispersion Correction.

    상태 변경:
        ctx.system, ctx.integrator, ctx.simulation, ctx.n_restrained 채움.

    Raises: MDStageError on ForceField 로드 / Restraint 적용 실패.
    """
    ff_files = ["amber14-all.xml", "amber14/tip3pfb.xml"]
    try:
        ff = ForceField(*ff_files, *ctx.ncaa_xmls)
    except (mm.OpenMMException, ValueError) as e:
        raise MDStageError(
            f"ForceField 로드 OpenMM 오류: {e}", log_suffix="_FAILED_FF.log"
        ) from e
    except Exception as e:
        raise MDStageError(
            f"ForceField 로드 예상치 못한 오류 ({e.__class__.__name__}): {e}",
            log_suffix="_FAILED_FF.log",
        ) from e

    ctx.modeller.addSolvent(ff, model="tip3p", padding=1.2 * unit.nanometers)
    print(f"  원자 수 (용매화 후): {ctx.modeller.topology.getNumAtoms()}")

    ctx.system = ff.createSystem(
        ctx.modeller.topology,
        nonbondedMethod    = PME,
        nonbondedCutoff    = 1.0 * unit.nanometers,
        constraints        = HBonds,
        rigidWater         = True,
        ewaldErrorTolerance= 0.0005,
    )

    ctx.n_restrained = 0
    junction_resids = set()
    if ctx.topology_type in ("cyclic_htc", "cyclic_nm"):
        binder_residues = [r for c in ctx.modeller.topology.chains() if c.id == ctx.binder_chain for r in c.residues()]
        if len(binder_residues) >= 4:
            junction_resids = {
                binder_residues[0].index, binder_residues[1].index,
                binder_residues[-1].index, binder_residues[-2].index,
            }
            junction_names = [f"{r.name}{r.id}" for r in (binder_residues[:2] + binder_residues[-2:])]
            print(f"  [v46 FIX] Junction 잔기 {len(junction_resids)}개 restraint 제외: {junction_names}")
    ctx.junction_resids = junction_resids

    if ctx.xml_res_name:
        try:
            ctx.system, ctx.n_restrained = add_position_restraints(
                ctx.system, ctx.modeller.positions, ctx.modeller.topology,
                ctx.xml_res_name, k_restrain=100.0, junction_resids=junction_resids,
            )
        except RuntimeError as e:
            raise MDStageError(str(e), log_suffix="_FAILED_RESTRAINT.log") from e

    is_backbone_closed = ctx.cycle_meta.get("htc_formed", False)
    if ctx.disp_corr == "off" or (ctx.disp_corr == "auto" and is_backbone_closed):
        for f in ctx.system.getForces():
            if isinstance(f, mm.NonbondedForce):
                f.setUseDispersionCorrection(False)
        print("  [Physics] NonbondedForce Dispersion Correction: 비활성화됨 (OFF)")
    else:
        print("  [Physics] NonbondedForce Dispersion Correction: 활성화됨 (ON)")

    ctx.integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picoseconds, ctx.dt_fs * unit.femtoseconds
    )
    print(f"  [MD] Integrator timestep: {ctx.dt_fs}fs")
    if ctx.seed != -1:
        ctx.integrator.setRandomNumberSeed(ctx.seed)

    platform, plat_name = get_best_platform(preferred=getattr(ctx, "platform_pref", "CUDA"))
    ctx.simulation = app.Simulation(
        ctx.modeller.topology, ctx.system, ctx.integrator,
        platform, get_platform_properties(plat_name),
    )
    ctx.simulation.context.setPositions(ctx.modeller.positions)


def run_cyclic_relaxation(ctx: "_MDContext") -> None:
    """[v46 FIX] HTC junction angle/dihedral 교정을 위한 Cyclization Relaxation.

    HTC bond 형성 직후 junction 잔기(첫/마지막 2잔기)만 자유롭게 두고
    나머지 모든 원자에 강한 position restraint(k=5000)를 적용한 뒤
    energy minimization을 수행한다. 이를 통해 cyclization junction의
    C-N-CA angle과 ω dihedral이 물리적 값으로 수렴된다.

    Raises: MDStageError on relaxation 실패.
    """
    junction_resids = getattr(ctx, "junction_resids", set())
    if ctx.topology_type not in ("cyclic_htc", "cyclic_nm"):
        return

    binder_residues = [r for c in ctx.modeller.topology.chains() if c.id == ctx.binder_chain for r in c.residues()]
    if len(binder_residues) < 4:
        return

    if not junction_resids:
        junction_resids = {
            binder_residues[0].index, binder_residues[1].index,
            binder_residues[-1].index, binder_residues[-2].index,
        }

    print("  [Cyclic Relax] HTC junction 다단계 relaxation 시작...")

    k_schedule = [(5000.0, 5000, "강한 restraint"), (2000.0, 3000, "중간 restraint"), (500.0, 2000, "약한 restraint")]

    for k_val, max_iter, label in k_schedule:
        relax_force = mm.CustomExternalForce("0.5*k_relax*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        relax_force.addPerParticleParameter("k_relax")
        relax_force.addPerParticleParameter("x0")
        relax_force.addPerParticleParameter("y0")
        relax_force.addPerParticleParameter("z0")

        positions = ctx.simulation.context.getState(getPositions=True).getPositions()
        n_relax_restrained = 0
        n_junction_free = 0
        for atom in ctx.modeller.topology.atoms():
            if atom.residue.index in junction_resids:
                n_junction_free += 1
                continue
            pos = positions[atom.index]
            relax_force.addParticle(atom.index, [k_val, pos.x, pos.y, pos.z])
            n_relax_restrained += 1

        force_idx = ctx.system.addForce(relax_force)
        ctx.simulation.context.reinitialize(preserveState=True)

        try:
            ctx.simulation.minimizeEnergy(maxIterations=max_iter)
            state = ctx.simulation.context.getState(getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print(f"    {label} (k={k_val}, iter={max_iter}): PE = {pe:.2f} kJ/mol")
        except (mm.OpenMMException, ValueError) as e:
            ctx.system.removeForce(force_idx)
            ctx.simulation.context.reinitialize(preserveState=True)
            raise MDStageError(
                f"Cyclization Relaxation 폭발: {e}",
                log_suffix="_EXPLODED_CYCLIC_RELAX.log",
                partial_stage="CyclicRelax", partial_suffix="_partial_relax.pdb",
            ) from e

        ctx.system.removeForce(force_idx)
        ctx.simulation.context.reinitialize(preserveState=True)

    print(f"    (restrained={n_relax_restrained}, free={n_junction_free})")
    print("  [Cyclic Relax] Junction relaxation 완료")


def run_equilibration(ctx: "_MDContext", k_stages=(100.0, 500.0, 1000.0)) -> None:
    """[v46] 3 단계 Staged Restraint Minimization + NVT 워밍업.

    Cyclic 구조일 때 NVT 워밍업을 연장한다 (5K→300K, 120ps).

    Raises: MDStageError on Minimization/NVT 폭발 (좌표 NaN 등).
    """
    print("  [Stage 1] Staged Minimization...")
    k_schedule = list(k_stages) if ctx.n_restrained > 0 else [None]
    try:
        for k_val in k_schedule:
            if k_val is not None:
                ctx.simulation.context.setParameter("k", k_val)
            ctx.simulation.minimizeEnergy(maxIterations=1000)
            state = ctx.simulation.context.getState(getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print(f"    k={k_val if k_val is not None else 'N/A'}: PE = {pe:.2f} kJ/mol")
    except (mm.OpenMMException, ValueError) as e:
        raise MDStageError(
            f"Minimization OpenMM 폭발: {e}",
            log_suffix="_EXPLODED_MIN.log",
            partial_stage="Minimization", partial_suffix="_partial_min.pdb",
        ) from e
    except Exception as e:
        raise MDStageError(
            f"Minimization 예상치 못한 오류 ({e.__class__.__name__}): {e}",
            log_suffix="_EXPLODED_MIN.log",
            partial_stage="Minimization", partial_suffix="_partial_min.pdb",
        ) from e

    is_cyclic = ctx.topology_type in ("cyclic_htc", "cyclic_nm")

    if is_cyclic:
        print("  [Stage 1.5] Cyclic 추가 자유 minimization (restraint 완화)...")
        try:
            if ctx.n_restrained > 0:
                ctx.simulation.context.setParameter("k", 10.0)
            ctx.simulation.minimizeEnergy(maxIterations=5000)
            state = ctx.simulation.context.getState(getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print(f"    k=10.0 (약한 restraint): PE = {pe:.2f} kJ/mol")

            if ctx.n_restrained > 0:
                ctx.simulation.context.setParameter("k", 0.0)
            ctx.simulation.minimizeEnergy(maxIterations=5000)
            state = ctx.simulation.context.getState(getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print(f"    k=0.0 (자유 minimization): PE = {pe:.2f} kJ/mol")

            if ctx.n_restrained > 0:
                ctx.simulation.context.setParameter("k", 100.0)
        except (mm.OpenMMException, ValueError) as e:
            raise MDStageError(
                f"Cyclic 추가 minimization 폭발: {e}",
                log_suffix="_EXPLODED_CYCLIC_MIN.log",
                partial_stage="CyclicMin", partial_suffix="_partial_cycmin.pdb",
            ) from e

    if is_cyclic:
        print("  [Stage 2a] Cyclic 초저온 NVT (dt=1fs, 1K → 10K, 20ps)...")
        try:
            orig_dt = ctx.integrator.getStepSize()
            ctx.integrator.setStepSize(1.0 * unit.femtoseconds)
            ctx.simulation.context.setVelocitiesToTemperature(1.0 * unit.kelvin)
            for temp_10x in range(10, 110, 10):
                temp = temp_10x / 10.0
                ctx.integrator.setTemperature(temp * unit.kelvin)
                ctx.simulation.step(2000)
            print(f"    초저온 phase 완료 (10K 도달, dt=1fs)")
            ctx.integrator.setStepSize(orig_dt)

            nvt_start_temp = 10
            nvt_step_size = 5
            nvt_steps_per_temp = 2000
            nvt_label = "10K → 300K, 120ps"
        except (mm.OpenMMException, ValueError) as e:
            raise MDStageError(
                f"Cyclic 초저온 NVT 폭발: {e}",
                log_suffix="_EXPLODED_COLD_NVT.log",
                partial_stage="ColdNVT", partial_suffix="_partial_cold_nvt.pdb",
            ) from e
    else:
        nvt_start_temp = 10
        nvt_step_size = 10
        nvt_steps_per_temp = 1000
        nvt_label = "10K → 300K, 60ps"

    print(f"  [Stage 2] NVT 워밍업 ({nvt_label}, dt=2fs)...")
    try:
        if not is_cyclic:
            ctx.simulation.context.setVelocitiesToTemperature(nvt_start_temp * unit.kelvin)
        for temp in range(nvt_start_temp, 310, nvt_step_size):
            ctx.integrator.setTemperature(temp * unit.kelvin)
            ctx.simulation.step(nvt_steps_per_temp)
    except (mm.OpenMMException, ValueError) as e:
        raise MDStageError(
            f"NVT 워밍업 OpenMM 폭발: {e}",
            log_suffix="_EXPLODED_NVT.log",
            partial_stage="NVT", partial_suffix="_partial_nvt.pdb",
        ) from e
    except Exception as e:
        raise MDStageError(
            f"NVT 워밍업 예상치 못한 오류 ({e.__class__.__name__}): {e}",
            log_suffix="_EXPLODED_NVT.log",
            partial_stage="NVT", partial_suffix="_partial_nvt.pdb",
        ) from e


_CHECKPOINT_INTERVAL = 50000  # 100ps (50000 steps × 2fs) or 50ps at 1fs
_MIN_COMPLETION_RATIO = 0.20  # 20% 이상 완료 시 PARTIAL_SUCCESS


def run_production(ctx: "_MDContext") -> None:
    """[v47] NPT 전환, Reporter 등록, 체크포인트 기반 Production MD, Recovery, PDB 저장.

    체크포인트 간격(_CHECKPOINT_INTERVAL)마다 정상 상태를 기록하고,
    NaN 폭발 시 마지막 정상 상태를 _final.pdb로 저장하여 부분 궤적을 보존한다.
    """
    print("  [Stage 3] NPT 전환...")
    ctx.integrator.setTemperature(300 * unit.kelvin)
    ctx.system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 300 * unit.kelvin))
    ctx.simulation.context.reinitialize(preserveState=True)
    print("  NPT 전환 완료.")

    dt_fs = ctx.dt_fs
    dcd_interval = max(1, ctx.n_steps // 500)
    ctx.simulation.reporters.append(DCDReporter(ctx.out_dcd, dcd_interval))
    ctx.simulation.reporters.append(StateDataReporter(
        ctx.out_log, max(1, ctx.n_steps // 100),
        step=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, volume=True, progress=True,
        remainingTime=True, totalSteps=ctx.n_steps, separator="\t",
    ))
    ctx.simulation.reporters.append(StateDataReporter(
        sys.stdout, max(1, ctx.n_steps // 10),
        step=True, potentialEnergy=True, temperature=True,
        progress=True, totalSteps=ctx.n_steps,
    ))

    ps_total = ctx.n_steps * dt_fs / 1000.0
    print(f"  MD 실행 중 ({ctx.n_steps} steps = {ps_total:.1f} ps, dt={dt_fs}fs)...")

    n_chunks = ctx.n_steps // _CHECKPOINT_INTERVAL
    remainder = ctx.n_steps % _CHECKPOINT_INTERVAL
    completed_steps = 0
    last_good_state = None
    nan_occurred = False
    log_every = max(1, n_chunks // 10)

    for chunk_idx in range(n_chunks):
        try:
            ctx.simulation.step(_CHECKPOINT_INTERVAL)
        except Exception as e:
            nan_occurred = True
            completed_ps = completed_steps * dt_fs / 1000.0
            print(f"  [!] Production MD 폭발 at step {completed_steps} ({completed_ps:.0f} ps): {e}")
            break
        completed_steps = (chunk_idx + 1) * _CHECKPOINT_INTERVAL
        last_good_state = ctx.simulation.context.getState(getPositions=True, getEnergy=True)
        if log_every and (chunk_idx + 1) % log_every == 0:
            pct = completed_steps / ctx.n_steps * 100
            pe = last_good_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            print(f"  [Checkpoint {pct:.0f}%] step={completed_steps}, PE={pe:.0f} kJ/mol")

    if not nan_occurred and remainder > 0:
        try:
            ctx.simulation.step(remainder)
            completed_steps += remainder
            last_good_state = ctx.simulation.context.getState(getPositions=True, getEnergy=True)
        except Exception:
            nan_occurred = True

    completion_ratio = completed_steps / ctx.n_steps if ctx.n_steps > 0 else 0.0

    if nan_occurred:
        exploded_log = os.path.join(ctx.output_dir,
                                    f"{ctx.basename}_EXPLODED_dt{dt_fs:.0f}fs.log")
        with open(exploded_log, "w", encoding="utf-8") as ef:
            ef.write(f"NaN at step {completed_steps}/{ctx.n_steps} "
                     f"({completion_ratio*100:.0f}%), dt={dt_fs}fs\n")

        if last_good_state is not None and completion_ratio >= _MIN_COMPLETION_RATIO:
            with open(ctx.out_pdb, "w", encoding="utf-8") as f:
                PDBFile.writeFile(ctx.modeller.topology, last_good_state.getPositions(), f)
            completed_ps = completed_steps * dt_fs / 1000.0
            print(f"  [Recovery] {completion_ratio*100:.0f}% 완료 ({completed_ps:.0f}ps). "
                  f"마지막 정상 상태를 {os.path.basename(ctx.out_pdb)}에 저장.")
            print(f"  [EXPLODED] 마커 생성: {os.path.basename(exploded_log)}")
            ctx.production_partial = True
        else:
            raise MDStageError(
                f"Production MD 폭발 ({completion_ratio*100:.0f}% 완료, 최소 {_MIN_COMPLETION_RATIO*100:.0f}% 미달)",
                log_suffix="_EXPLODED_MD.log",
                partial_stage="Production MD", partial_suffix="_partial_md.pdb",
            )
    else:
        state = ctx.simulation.context.getState(getPositions=True)
        with open(ctx.out_pdb, "w", encoding="utf-8") as f:
            PDBFile.writeFile(ctx.modeller.topology, state.getPositions(), f)
        ctx.production_partial = False
    print(f"  최종 구조 저장 완료: {os.path.basename(ctx.out_pdb)}")

    if ctx.cycle_meta.get("is_cyclized"):
        append_conect_records(ctx.out_pdb, ctx.cycle_meta)

    if ctx.hdd_path and os.path.exists(ctx.hdd_path):
        dest_dcd = os.path.join(ctx.hdd_path, os.path.basename(ctx.out_dcd))
        try:
            shutil.move(ctx.out_dcd, dest_dcd)
            print(f"  [Auto-Move] {os.path.basename(ctx.out_dcd)} → HDD 이동 완료")
        except (OSError, shutil.Error) as e:
            raise MDStageError(
                f"DCD 이동 실패: {e}",
                status="SUCCESS_WITH_WARNING",
                log_suffix="_WARNING_DCD_MOVE.log",
            ) from e


def run_md(pdb_path, params_manifest, output_dir, n_steps, topology_type, ncaa_label, ncaa_code, hdd_path="", binder_chain="B", seed=42, disp_corr="auto", target_resid="", graph_policy="strict", platform_pref="CUDA", dt_fs=2.0):
    """[v48] 4 단계 헬퍼의 얇은 오케스트레이터.

    Returns:
        str: "SUCCESS" / "SUCCESS_WT" / "SUCCESS_WITH_WARNING"
             / "PARTIAL_SUCCESS" / "SKIP" / "FAIL"
    """
    ctx = _MDContext(
        pdb_path=pdb_path, params_manifest=params_manifest, output_dir=output_dir,
        n_steps=n_steps, topology_type=topology_type, ncaa_label=ncaa_label,
        ncaa_code=ncaa_code, hdd_path=hdd_path, binder_chain=binder_chain,
        seed=seed, disp_corr=disp_corr, target_resid=target_resid,
        graph_policy=graph_policy, platform_pref=platform_pref,
        dt_fs=dt_fs,
    )

    exploded_markers = glob.glob(os.path.join(output_dir, f"{ctx.basename}_EXPLODED_*.log"))
    if os.path.exists(ctx.out_pdb) and not exploded_markers:
        ctx.log_status(f"이미 완료된 MD 결과 존재: {ctx.basename}", "SKIP", "_SKIPPED_ALREADY_DONE.log")
        return "SKIP"

    try:
        prepare_structure(ctx)
        build_openmm_system(ctx)
        run_cyclic_relaxation(ctx)
        run_equilibration(ctx)
        run_production(ctx)
    except MDStageError as stage_err:
        ctx.log_status(str(stage_err), "WARNING" if stage_err.status == "SUCCESS_WITH_WARNING" else "FAIL", stage_err.log_suffix)
        if stage_err.partial_stage and ctx.simulation is not None:
            ctx.save_partial_pdb(stage_err.partial_stage, stage_err.partial_suffix)
        if stage_err.status == "SUCCESS_WITH_WARNING":
            return "SUCCESS_WITH_WARNING"
        return "FAIL"

    if getattr(ctx, "production_partial", False):
        return "PARTIAL_SUCCESS"

    if ctx.ncaa_label.lower() == "none":
        return "SUCCESS_WT" if not ctx.warnings_found else "SUCCESS_WITH_WARNING"
    return "SUCCESS" if not ctx.warnings_found else "SUCCESS_WITH_WARNING"


# ==========================================
# CLI 메인
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Restrained MD — ncAA 치환 구조 및 야생형(WT) 단백질 일반 MD (v28)")
    parser.add_argument("--inputdir",  required=True)
    parser.add_argument("--params_manifest",  required=False, default="", help="파라미터 매니페스트 (WT 모드일 경우 생략 가능)")
    parser.add_argument("--outputdir", required=True)
    parser.add_argument("--steps",     type=int, default=50000)
    parser.add_argument("--topology",  default="linear", help="linear, cyclic_htc(또는 cyclic_nm), cyclic_ss, bicyclic")
    parser.add_argument("--ncaa_label", default="none", help="사용자 지정 ncAA 라벨명 (none=야생형 일반 MD 수행)")
    parser.add_argument("--ncaa_code",  default="", help="XML 내부의 실제 PDB 잔기명. --ncaa_label 과 다를 때 엄격한 검증을 위해 사용됨.")
    parser.add_argument("--binder_chain", type=str, default="B", help="Binder 체인 ID (기본: B)")
    parser.add_argument("--target_resid", type=str, default="", help="ncAA가 위치한 명시적 Residue ID (예: 15). 생략 시 자동 탐색")
    parser.add_argument("--hdd_path",  type=str, default="")
    parser.add_argument("--seed",      type=int, default=42, help="Langevin Integrator 난수 시드 (-1: 무작위)")
    parser.add_argument("--dispersion", choices=["auto", "on", "off"], default="auto", help="NonbondedForce 분산 보정 (auto: HTC 고리형 펩타이드일 때만 off)")
    parser.add_argument("--graph_policy", choices=["strict", "warn", "relaxed"], default="strict", help="[v28 기본값 상향] Graph matching 시 거리 기반 추정 엣지 허용 정책")
    parser.add_argument("--platform", choices=["CUDA", "OpenCL", "CPU", "auto"], default="CUDA",
                        help="[v4 H-1] OpenMM 플랫폼 우선순위 (기본 CUDA mixed precision). auto 는 v3 의 기본 폴백 체인.")
    parser.add_argument("--dt_fs", type=float, default=2.0,
                        help="[v48] Production timestep (fs). 기본 2.0, cyclic 재시도 시 1.0")
    args = parser.parse_args()

    os.makedirs(args.outputdir, exist_ok=True)

    pdb_files = [
        p for p in glob.glob(os.path.join(args.inputdir, "**", "*.pdb"), recursive=True)
        if not any(x in p for x in ("_renum.pdb", "_minimized.pdb", "_partial_"))
    ]

    if not pdb_files:
        print(f"[!] PDB 파일 없음: {args.inputdir}")
        return

    actual_ncaa_code = args.ncaa_code.strip()

    print(f"\n[Restrained MD v28 - The Absolute Integrity]")
    print(f"  대상 PDB  : {len(pdb_files)}개")
    print(f"  스텝 수   : {args.steps} ({args.steps * 2 / 1000:.1f} ps)")
    print(f"  토폴로지  : {args.topology}")
    print(f"  ncAA 라벨 : {args.ncaa_label} (Code: {actual_ncaa_code if actual_ncaa_code else '(미지정: XML 이름 검증 생략)'})")
    print(f"  바인더    : 체인 {args.binder_chain}")
    print(f"  난수 시드 : {'무작위' if args.seed == -1 else args.seed}")
    print(f"  분산 보정 : {args.dispersion}")
    print(f"  그래프정책: {args.graph_policy}")
    print(f"  플랫폼    : {args.platform}")
    print(f"  timestep  : {args.dt_fs}fs")

    success_cnt, wt_cnt, warn_cnt, skip_cnt, fail_cnt, partial_cnt = 0, 0, 0, 0, 0, 0

    for pdb_path in pdb_files:
        status = run_md(
            pdb_path     = pdb_path,
            params_manifest = args.params_manifest,
            output_dir   = args.outputdir,
            n_steps      = args.steps,
            topology_type= args.topology,
            ncaa_label   = args.ncaa_label,
            ncaa_code    = actual_ncaa_code,
            hdd_path     = args.hdd_path,
            binder_chain = args.binder_chain,
            seed         = args.seed,
            disp_corr    = args.dispersion,
            target_resid = args.target_resid,
            graph_policy = args.graph_policy,
            platform_pref= args.platform,
            dt_fs        = args.dt_fs,
        )
        if status == "SUCCESS": success_cnt += 1
        elif status == "SUCCESS_WT": wt_cnt += 1
        elif status == "SUCCESS_WITH_WARNING": warn_cnt += 1
        elif status == "PARTIAL_SUCCESS": partial_cnt += 1
        elif status == "SKIP": skip_cnt += 1
        elif status == "FAIL": fail_cnt += 1

    print(f"\n========================================")
    print(f"  [MD 배치 완료 통계]")
    print(f"  - 총 입력 파일 : {len(pdb_files)} 개")
    print(f"  - 성공 (ncAA MD) : {success_cnt} 개")
    if wt_cnt > 0:
        print(f"  - 성공 (야생형 MD) : {wt_cnt} 개")
    if warn_cnt > 0:
        print(f"  - 성공 (경고 동반) : {warn_cnt} 개 (로그 확인 요망)")
    if partial_cnt > 0:
        print(f"  - 부분 성공 (Recovery) : {partial_cnt} 개 (부분 궤적 유효)")
    print(f"  - 스킵 (완료됨 등): {skip_cnt} 개")
    print(f"  - 실패 (폭발 등) : {fail_cnt} 개")
    print(f"========================================")

if __name__ == "__main__":
    main()
