#!/usr/bin/env python
"""
parameterize_ncaa.py
--------------------
SMILES 기반 ncAA GAFF2+RESP 파라미터화.
UPDD의 ncaa_registry 정보를 받아 Manifest.json을 출력하고 하드코딩된 변수를 제거합니다.
"""

import os
import sys
import re
import argparse
import subprocess
import logging
import json
import copy
import xml.etree.ElementTree as ET
import numpy as np

# utils 디렉토리 내부에서 직접 registry 참조
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from ncaa_registry import resolve_ncaa_definition
except ImportError:
    sys.exit("[!] ncaa_registry.py 파일을 찾을 수 없습니다. utils 폴더에 위치해야 합니다.")

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("parameterize")

try:
    from pyscf import gto, scf, lib
    lib.num_threads(os.cpu_count() or 8)
    HAS_PYSCF = True
except ImportError:
    HAS_PYSCF = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    log.error("[RDKit] 미설치 — parameterize_ncaa 실행 불가")
    sys.exit(1)

try:
    import parmed as pmd
    HAS_PARMED = True
except ImportError:
    HAS_PARMED = False

GB_RADII = {
    "C": 1.70, "N": 1.55, "O": 1.52, "H": 1.20, "S": 1.80, "P": 1.80,
    "F": 1.47, "Cl": 1.75, "Br": 1.85, "Si": 2.10,
}

def _emit_hydrogen_definitions(xml_path: str, res_name: str, hyd_path: str) -> bool:
    """[v36 CRITICAL FIX-2] OpenMM ``Modeller.addHydrogens()`` 가 ncAA 잔기에
    수소를 추가할 수 있도록, 잔기 템플릿으로부터 Hydrogen Definitions XML 을
    유도하여 별도 파일로 출력한다.

    배경:
        OpenMM 의 ``Modeller.addHydrogens()`` 는 ForceField 의 잔기 템플릿이
        아닌, 별도의 ``Modeller._residueHydrogens`` 데이터베이스
        (`Hydrogens.xml`) 를 참조하여 어느 부모 heavy atom 에 어떤 H 를
        추가할지 결정한다. 표준 아미노산 (ALA, GLY, …) 만 기본 등록되어
        있으며, ncAA 잔기명은 등록되어 있지 않다. 등록되지 않은 잔기는
        ``addHydrogens`` 가 H 를 추가하지 않고 그대로 통과시키므로, 이후
        ``ForceField.createSystem()`` 의 템플릿 매칭이 (heavy-only 6 원자 vs
        heavy+H 14 원자) mismatch 로 실패한다 ("No template found ... is
        similar to ASP" 에러).

    근거:
        OpenMM 5.x ~ 7.x 의 Modeller.addHydrogens 구현은 다음 조건을 갖는다
        (lib/openmm/app/modeller.py 참조):

            if residue.name in Modeller._residueHydrogens or ...:
                <add hydrogens>

        즉, ncAA 잔기명을 ``_residueHydrogens`` 에 등록하기 전에는 H 가
        추가되지 않는다. ``Modeller.loadHydrogenDefinitions(file)`` 가 이를
        등록하는 공식 진입점이므로, 본 함수는 그에 맞는 XML 을 생성한다.

    데이터 계약 (Principle 7):
        - 출력 파일: ``{res_name}_hydrogens.xml`` (parameterize_ncaa 출력
          디렉토리 내) — manifest 의 ``hydrogens_path`` 키로 노출
        - 형식 (OpenMM Hydrogens.xml 표준):

            <Residues>
              <Residue name="NMA">
                <H name="HN" parent="N"/>
                <H name="HA" parent="CA"/>
                ...
              </Residue>
            </Residues>

        - 본 함수는 ``terminal`` 속성을 부여하지 않는다. 즉, 모든 H 를
          internal (비-terminal) 로 등록한다. ncAA 가 chain terminal 에
          위치하는 경우의 cap 처리는 별도 메커니즘 (apply_universal_xml_patch
          또는 사전 cyclization) 에 위임한다 — Principle 5 (영구 해법) 를
          위해 단순한 가정에 의존하지 않는다.

    Args:
        xml_path: ParmEd → _rename → _postprocess 가 완료된 잔기 템플릿 XML
        res_name: 잔기명 (예: "NMA")
        hyd_path: 출력 Hydrogen Definitions XML 경로

    Returns:
        bool: 정상 출력 시 True, 잔기 블록을 찾지 못해 no-op 인 경우 False.
    """
    if not os.path.exists(xml_path):
        log.warning(f"[Hydrogens] 입력 XML 부재: {xml_path}")
        return False

    tree = ET.parse(xml_path)
    root = tree.getroot()

    # [v5 해법 B] 단일 XML 안의 internal/N-term/C-term 3개 변형 잔기 모두에
    # 대응하는 H 정의를 출력한다. _generate_topology_variants 가 N{RES}/C{RES}
    # 블록을 추가하므로, addHydrogens 가 잔기명을 통해 어느 변형이 들어와도
    # H 를 부착할 수 있어야 한다. 변형 블록이 부재하면 해당 entry 만 생략한다
    # (Principle 7: 의존성 일관성 — 변형 미생성 환경에서도 무결성 유지).
    target_names = [res_name, f"N{res_name}", f"C{res_name}"]

    # AtomTypes 색인 — element 추론용 (XML 전체 단일 인스턴스)
    atomtype_to_elem = {
        t.get("name"): t.get("element")
        for t in root.findall("AtomTypes/Type")
    }

    out_root = ET.Element("Residues")
    emitted_count = 0
    base_h_count = 0

    for tname in target_names:
        res_node = next(
            (r for r in root.iter("Residue") if r.get("name") == tname),
            None,
        )
        if res_node is None:
            continue

        # 원자명 → element 매핑 (잔기 내부)
        name_to_elem: dict = {}
        for atom in res_node.findall("Atom"):
            aname = atom.get("name")
            atype = atom.get("type")
            name_to_elem[aname] = atomtype_to_elem.get(atype)

        # 결합으로부터 H → 부모 heavy atom 매핑 추출
        h_to_parent: dict = {}
        for bond in res_node.findall("Bond"):
            a1 = bond.get("atomName1")
            a2 = bond.get("atomName2")
            e1 = name_to_elem.get(a1)
            e2 = name_to_elem.get(a2)
            if e1 == "H" and e2 != "H":
                h_to_parent[a1] = a2
            elif e2 == "H" and e1 != "H":
                h_to_parent[a2] = a1

        if not h_to_parent:
            log.warning(f"[Hydrogens] {tname}: 추출된 H 결합 없음 — 해당 변형 entry 생략")
            continue

        out_res = ET.SubElement(out_root, "Residue", attrib={"name": tname})
        # 안정적 출력 순서: 부모 heavy atom 명 기준 알파벳 → H 이름
        for h_name in sorted(h_to_parent.keys(), key=lambda x: (h_to_parent[x], x)):
            ET.SubElement(out_res, "H", attrib={
                "name": h_name,
                "parent": h_to_parent[h_name],
            })
        emitted_count += 1
        if tname == res_name:
            base_h_count = len(h_to_parent)

    if emitted_count == 0:
        log.warning(f"[Hydrogens] {res_name}: 잔기 블록 미발견 — Hydrogens.xml 생성 생략")
        return False

    out_tree = ET.ElementTree(out_root)
    out_tree.write(hyd_path, xml_declaration=True, encoding="utf-8")
    log.info(
        f"[Hydrogens] {res_name} Hydrogen Definitions 작성 완료 "
        f"(변형={emitted_count}, base H={base_h_count}): {hyd_path}"
    )
    return True


def write_manifest(output_dir, args, ncaa_def, xml_out, resp_used, mol2, frcmod, hydrogens_path=""):
    # [v35 AUDIT] Principle 4/7: xml_resname MUST be sourced from the canonical
    # ncaa_registry (ncaa_def.xml_resname), NOT from args.ncaa_code. The Antigravity
    # Architecture Fix in UPDATE.md explicitly mandates xml_resname as the single
    # source of truth — prior code wrote args.ncaa_code here, so any registry entry
    # whose xml_resname diverges from its CLI code produced a manifest whose
    # xml_resname disagreed with the actual XML <Residue name="..."> block. The
    # downstream graph-isomorphism gate in run_restrained_md.py then loaded the
    # wrong residue name and Fail-Fast with a cryptic "template not found" error.
    # [v36] hydrogens_path 추가 — Modeller.addHydrogens 가 ncAA 에 H 를 추가할
    # 수 있도록 _emit_hydrogen_definitions 가 생성한 OpenMM Hydrogens.xml 의
    # 절대 경로를 manifest 에 노출한다 (Principle 7: 데이터 계약 일관성).
    manifest = {
        "ncaa_code": args.ncaa_code,
        "xml_resname": ncaa_def.xml_resname,
        "formal_charge": ncaa_def.formal_charge,
        "smiles_used": ncaa_def.smiles_free,
        "resp_used": resp_used,
        "xml_path": xml_out,
        "mol2_path": mol2,
        "frcmod_path": frcmod,
        "hydrogens_path": hydrogens_path,
        "status": "SUCCESS"
    }
    # [v35 AUDIT] Principle 7/v27: encoding="utf-8" enforced for manifest I/O —
    # output_dir may contain 한글 path segments on the operator's workstation.
    manifest_path = os.path.join(output_dir, f"{args.ncaa_code}_params_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    return manifest_path

def smiles_to_mol(smiles, output_mol, name="ncAA"):
    from rdkit.Chem import SDWriter
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
    # [v35 AUDIT] Principle 5: Silent 'except Exception: pass' replaced with a
    # logged warning. The scientific fallback (ETKDG-embedded geometry without
    # MMFF refinement) is preserved — only the silent swallow is eliminated, so
    # the operator now learns *why* RESP charges may be noisier on this ncAA.
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as _mmff_err:
        log.warning(f"[MMFF] optimization failed for {name}: {_mmff_err}. Falling back to raw ETKDG geometry.")
    with SDWriter(output_mol) as w: w.write(mol)
    return mol

def _mol_to_xyz(mol_path, xyz_path):
    from rdkit.Chem import SDMolSupplier
    suppl = SDMolSupplier(mol_path, removeHs=False)
    mol = next(iter(suppl))
    conf = mol.GetConformer()
    with open(xyz_path, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\nPySCF RESP\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol():2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")

def run_antechamber(mol_path, work_dir, res_name, charge=0):
    mol2 = os.path.join(work_dir, f"{res_name}.mol2")
    frcmod = os.path.join(work_dir, f"{res_name}.frcmod")
    cmd = ["antechamber", "-i", mol_path, "-fi", "mdl", "-o", mol2, "-fo", "mol2", "-c", "bcc", "-nc", str(charge), "-at", "gaff2", "-rn", res_name, "-dr", "no", "-s", "2"]
    subprocess.run(cmd, check=True, capture_output=True, cwd=work_dir)
    subprocess.run(["parmchk2", "-i", mol2, "-f", "mol2", "-o", frcmod, "-s", "gaff2"], check=True, capture_output=True)
    return mol2, frcmod

def _compute_electron_esp(mol, dm, grid_pts: np.ndarray) -> np.ndarray:
    # [v35 AUDIT] Principle 1: This fallback is a *scientifically valid*
    # physics routing — newer PySCF builds expose `int1e_grids` as a fast path,
    # while older builds require the `int3c2e` aux-e2 integral via `df.incore`.
    # Per Principle 1 the fallback must NOT be deleted; per Principle 5 we now
    # log the routing decision so operators can audit which integral code path
    # produced the ESP grid used in RESP fitting.
    try:
        j3c = mol.intor('int1e_grids', grids=grid_pts)
        return np.einsum('pij,ij->p', j3c, dm)
    except Exception as _int1e_err:
        log.info(f"[ESP] int1e_grids unavailable ({_int1e_err.__class__.__name__}); routing through df.incore int3c2e fallback.")
        from pyscf import df
        from pyscf.gto import fakemol_for_charges
        nao, n_pts = mol.nao, len(grid_pts)
        esp_elec = np.zeros(n_pts)
        for start in range(0, n_pts, 500):
            end = min(start + 500, n_pts)
            fmol = fakemol_for_charges(grid_pts[start:end])
            ints = df.incore.aux_e2(mol, fmol, 'int3c2e', aosym='s1', comp=1)
            if ints.ndim == 2:
                ints = ints.reshape(nao, nao, -1)
            esp_elec[start:end] = np.einsum('ij,ijn->n', dm, ints)
        return esp_elec

def _run_hf_calculation(xyz_path, charge):
    """[v3 5-2] HF/6-31G* SCF 계산을 수행하고 (mol, dm) 튜플을 반환한다.

    Args:
        xyz_path: PySCF 가 읽을 수 있는 XYZ 형식의 좌표 파일 경로
        charge:   분자 전체 정수 전하 (parameterize_ncaa 의 ncaa_def.formal_charge)

    Returns:
        Optional[Tuple[gto.Mole, np.ndarray]]: 수렴 시 (mol, density matrix),
        SCF 가 수렴하지 못하면 None.
    """
    atoms_str = ""
    with open(xyz_path, encoding="utf-8") as f:
        for line in f.readlines()[2:]:
            p = line.split()
            if len(p) == 4:
                atoms_str += f"{p[0]} {p[1]} {p[2]} {p[3]};"
    mol = gto.M(atom=atoms_str.rstrip(";"), basis="6-31G*", charge=charge, spin=0, unit="Angstrom", verbose=0)
    mf = scf.RHF(mol)
    mf.max_cycle, mf.conv_tol = 200, 1e-9
    mf.kernel()
    if not mf.converged:
        return None
    return mol, mf.make_rdm1()


def _build_esp_grid(mol, dm):
    """[v3 5-2] 원자별 다중 동심구 격자 위에서 정전기 전위(ESP) 를 계산한다.

    원자핵의 양전하 기여와 전자 밀도의 음전하 기여를 모두 합산하여, RESP
    피팅에 사용되는 grid 정전기 전위 배열을 반환한다.

    Args:
        mol: PySCF gto.Mole 객체 (수렴된 SCF 의 산물)
        dm:  HF density matrix (mol 과 동일한 좌표계)

    Returns:
        Tuple[np.ndarray, np.ndarray]: (grid_pts, esp) — grid 좌표(Bohr) 와
        해당 격자 위의 정전기 전위(원자 단위).
    """
    n_atoms = mol.natm
    grid_pts = []
    np.random.seed(42)
    for i in range(n_atoms):
        coord = mol.atom_coord(i)
        r_vdw = GB_RADII.get(mol.atom_symbol(i), 1.7) * 1.8897259886
        for scale in [1.4, 1.6, 1.8, 2.0]:
            for _ in range(100):
                phi, ct = np.random.uniform(0, 2 * np.pi), np.random.uniform(-1, 1)
                st = np.sqrt(1 - ct ** 2)
                grid_pts.append(coord + r_vdw * scale * np.array([st * np.cos(phi), st * np.sin(phi), ct]))

    grid_pts = np.array(grid_pts)
    esp = np.zeros(len(grid_pts))
    for i in range(n_atoms):
        d = np.linalg.norm(grid_pts - mol.atom_coord(i), axis=1)
        d = np.where(d < 1e-10, 1e-10, d)
        esp += mol.atom_charge(i) / d
    esp -= _compute_electron_esp(mol, dm, grid_pts)
    return grid_pts, esp


def _fit_resp_charges(esp, grid_pts, mol, charge):
    """[v3 5-2] 라그랑주 승수법으로 RESP 전하를 최소제곱 피팅한다.

    제약: sum(q_i) = total_charge (분자 전체 전하 보존)
    목적: || A · q − b ||^2 최소화 (ESP 재현)

    Args:
        esp:      _build_esp_grid 가 산출한 ESP 배열
        grid_pts: 동일 격자 좌표 (Bohr)
        mol:      PySCF gto.Mole 객체
        charge:   분자 전체 정수 전하

    Returns:
        np.ndarray: 원자별 RESP 전하 (n_atoms,) — 마지막 라그랑주 승수는 제외됨.
    """
    n_atoms = mol.natm
    A = np.zeros((n_atoms, n_atoms))
    b = np.zeros(n_atoms)
    for idx, pt in enumerate(grid_pts):
        di = np.array([max(np.linalg.norm(pt - mol.atom_coord(i)), 1e-10) for i in range(n_atoms)])
        A += np.outer(1 / di, 1 / di)
        b += esp[idx] / di

    A_ext = np.zeros((n_atoms + 1, n_atoms + 1))
    b_ext = np.zeros(n_atoms + 1)
    A_ext[:n_atoms, :n_atoms] = A
    A_ext[n_atoms, :n_atoms] = 1.0
    A_ext[:n_atoms, n_atoms] = 1.0
    b_ext[:n_atoms] = b
    b_ext[n_atoms] = charge
    return np.linalg.lstsq(A_ext, b_ext, rcond=None)[0][:n_atoms]


def run_pyscf_resp(xyz_path, charge=0):
    """[v3 5-2] PySCF 를 이용한 RESP 전하 도출 — 3단계 헬퍼의 얇은 오케스트레이터.

    1) _run_hf_calculation: HF/6-31G* SCF 수렴
    2) _build_esp_grid:     원자별 다중 동심구 ESP 격자 구축
    3) _fit_resp_charges:   라그랑주 승수법 최소제곱 RESP 피팅

    추가 안전장치:
        - PySCF 미설치 환경: 즉시 None 반환 (downstream AM1-BCC 폴백)
        - SCF 미수렴: None 반환
        - |q|_max > 2.0 e: charge-sanity Fail-Fast (downstream AM1-BCC 폴백)
        - 예기치 못한 예외: 경고 로깅 후 None 반환
    """
    if not HAS_PYSCF:
        log.warning("[RESP] PySCF unavailable — skipping RESP charge derivation; downstream will fall back to AM1-BCC via antechamber.")
        return None
    try:
        hf_result = _run_hf_calculation(xyz_path, charge)
        if hf_result is None:
            return None
        mol, dm = hf_result

        grid_pts, esp = _build_esp_grid(mol, dm)
        charges = _fit_resp_charges(esp, grid_pts, mol, charge)

        # [v35 AUDIT] Principle 5: log charge-sanity Fail-Fast with the offending
        # magnitude so the operator can trace RESP divergence on this ncAA
        # (mirrors the run_restrained_md.py v15 |q| > 2.0 threshold).
        max_q = float(np.max(np.abs(charges)))
        if max_q > 2.0:
            log.warning(f"[RESP] rejected: peak |q|={max_q:.3f} e > 2.0 threshold — downstream will fall back to AM1-BCC.")
            return None
        return charges, mol
    except Exception as _resp_err:
        log.warning(f"[RESP] PySCF routine crashed ({_resp_err.__class__.__name__}: {_resp_err}). Falling back to AM1-BCC.")
        return None

def _inject_resp_into_mol2(mol2_path, resp_result):
    if not HAS_PARMED or resp_result is None:
        return False
    resp_charges, _ = resp_result
    mol2_obj = pmd.load_file(mol2_path)
    atoms = list(list(mol2_obj.to_library().values())[0].atoms) if isinstance(mol2_obj, pmd.modeller.ResidueTemplateContainer) else list(mol2_obj.residues[0].atoms) if hasattr(mol2_obj, "residues") else list(mol2_obj.atoms)
    if len(resp_charges) != len(atoms):
        return False
    for i, atom in enumerate(atoms):
        atom.charge = float(resp_charges[i])
    mol2_obj.save(mol2_path, overwrite=True)
    return True

def _run_tleap(mol2_path, frcmod_path, work_dir, res_name):
    # [v35 AUDIT] Principle 4: tleap invocation now uses check=True so a failed
    # leap.in (invalid mol2, missing GAFF2 parameters, radical spin mismatch)
    # triggers Fail-Fast with the captured stdout/stderr rather than silently
    # producing an empty/corrupt prmtop that downstream parmed loads and
    # emits malformed OpenMM XML from.
    prmtop = os.path.join(work_dir, f"{res_name}.prmtop")
    inpcrd = os.path.join(work_dir, f"{res_name}.inpcrd")
    leap_in = os.path.join(work_dir, "leap.in")
    with open(leap_in, "w", encoding="utf-8") as f:
        f.write(f"source leaprc.gaff2\nloadamberparams {frcmod_path}\nmol = loadmol2 {mol2_path}\ncheck mol\nsaveamberparm mol {prmtop} {inpcrd}\nquit\n")
    try:
        subprocess.run(["tleap", "-f", leap_in], check=True, capture_output=True, text=True, cwd=work_dir)
    except subprocess.CalledProcessError as _leap_err:
        log.error(f"[tleap] Non-zero exit for residue {res_name}. stdout:\n{_leap_err.stdout}\nstderr:\n{_leap_err.stderr}")
        raise RuntimeError(f"tleap failed for {res_name}: {_leap_err}") from _leap_err
    if not (os.path.exists(prmtop) and os.path.exists(inpcrd)):
        raise RuntimeError(f"tleap completed but prmtop/inpcrd missing for {res_name}: {prmtop}, {inpcrd}")
    return prmtop, inpcrd

def _build_pdb_name_map(atoms, ncaa_mutation_type) -> dict:
    """[v3 5-1] mol2 atoms 리스트를 분석하여 GAFF2 내부 원자명 → PDB 규약명 매핑 dict 를 생성한다.

    백본 N/CA/C/O/CB/CM 식별 및 측쇄 H 원자명 추론을 수행한다. mutation type
    이 "N-M" 이 아니거나, 핵심 anchor (N) 식별이 실패하면 빈 dict 를 반환한다.

    Args:
        atoms: parmed Atom 객체 리스트
        ncaa_mutation_type: ncaa_def.mutation_type

    Returns:
        dict: {gaff2_name: pdb_name} 매핑. 빈 dict 는 rename 비대상 의미.
    """
    if ncaa_mutation_type != "N-M":
        return {}

    heavy = [a for a in atoms if a.atomic_number > 1]
    name_map: dict = {}

    def heavy_neighbors(atom):
        return [
            (b.atom1 if b.atom2 is atom else b.atom2)
            for b in atom.bonds
            if (b.atom1 if b.atom2 is atom else b.atom2).atomic_number > 1
        ]

    n_atom = next((a for a in heavy if a.atomic_number == 7), None)
    if not n_atom:
        return {}
    name_map[n_atom.name] = "N"

    ca_atom, c_atom = None, None
    for c in [a for a in heavy_neighbors(n_atom) if a.atomic_number == 6]:
        if len(heavy_neighbors(c)) <= 1:
            name_map[c.name] = "CM"
        else:
            name_map[c.name] = "CA"
            ca_atom = c

    if ca_atom is not None:
        for nbr in heavy_neighbors(ca_atom):
            if nbr is n_atom or nbr.name in name_map:
                continue
            if nbr.atomic_number == 6:
                if any(x.atomic_number == 8 for x in heavy_neighbors(nbr)):
                    name_map[nbr.name] = "C"
                    c_atom = nbr
                else:
                    name_map[nbr.name] = "CB"

    if c_atom is not None:
        for o in heavy_neighbors(c_atom):
            if o.atomic_number != 8:
                continue
            has_h = any((b.atom1 if b.atom2 is o else b.atom2).atomic_number == 1 for b in o.bonds)
            name_map[o.name] = "OXT" if has_h else "O"

    # 수소 원자 명명: 부모 heavy atom 의 PDB 이름에 따라 결정
    H_NAME_TABLE = {"CM": ["HM1", "HM2", "HM3"], "CB": ["HB1", "HB2", "HB3"]}
    h_idx: dict = {}
    for atom in atoms:
        if atom.atomic_number != 1:
            continue
        parents = [b.atom1 if b.atom2 is atom else b.atom2 for b in atom.bonds]
        if not parents:
            continue
        pname = name_map.get(parents[0].name, parents[0].name)
        if pname in H_NAME_TABLE:
            idx = h_idx.get(pname, 0)
            tbl = H_NAME_TABLE[pname]
            name_map[atom.name] = tbl[idx] if idx < len(tbl) else f"HX{idx}"
            h_idx[pname] = idx + 1
        elif pname == "N":
            name_map[atom.name] = "HN"
        elif pname == "CA":
            name_map[atom.name] = "HA"
        elif pname == "OXT":
            name_map[atom.name] = "HXT"

    return name_map


# [v0.6.6 Strategy A] MTR amber14SB charge patch table.
# Source: /openmm/app/data/amber14/protein.ff14SB.xml <Residue name="TRP"> (verified
# sum = 0.0000). Reference: Maier 2015 doi:10.1021/acs.jctc.5b00255.
# HE1 (+0.3412) is removed in MTR (methylated); its charge is redistributed to
# CM + 3×HM with total = +0.3412 so that MTR net charge preserves TRP's 0.0000.
# NE1 kept at ff14SB value (−0.3418) — local dipole preserved per the RESP
# group-equivalence principle (Bayly et al. 1993 §III.B
# doi:10.1021/j100142a004; Cieplak et al. 1995 §II.B
# doi:10.1002/jcc.540161106). The earlier inline citation to Capece 2012
# (doi:10.1021/jp2082825) was a misattribution: that paper is a heme
# dioxygenase reaction-mechanism study (IDO QM/MM), not RESP methodology.
# PDF re-acquisition was attempted 2026-04-28 (see
# outputs/analysis/capece_2012_acquisition_log_20260428.md) and the
# attribution was re-traced to the actual RESP source pair above.
_MTR_AMBER14_CHARGES = {
    "N":  -0.4157, "H":   0.2719,
    "CA": -0.0275, "HA":  0.1123,
    "CB": -0.0050, "HB2": 0.0339, "HB3": 0.0339,
    "CG": -0.1415,
    "CD1":-0.1638, "HD1": 0.2062,
    "NE1":-0.3418,
    "CE2": 0.1380,
    "CD2": 0.1243,
    "CE3":-0.2387, "HE3": 0.1700,
    "CZ3":-0.1972, "HZ3": 0.1447,
    "CH2":-0.1134, "HH2": 0.1417,
    "CZ2":-0.2601, "HZ2": 0.1572,
    "C":   0.5973, "O":  -0.5679,
    # Extension atoms (CM + HMs) sum = +0.3412 (replaces removed HE1)
    # Split: CM carries residual (0.3412 - 3*HM_q), HMs equal.
    # HM_q = 0.0975 (consistent with Khoury et al. 2014 Forcefield_PTM RESP-A2
    # fit for N-methyl side-chain protons, doi:10.1021/sb400168u; cf. Maier
    # 2015 ff14SB aliphatic-CH₃ HC range 0.07–0.11 e). Local RESP would
    # refine to ±0.05 e. Earlier "Capece 2012 indicative" attribution
    # corrected 2026-04-28 — see acquisition log
    # outputs/analysis/capece_2012_acquisition_log_20260428.md.
    "CM":  0.0487,
    "HM1": 0.0975, "HM2": 0.0975, "HM3": 0.0975,
}


# [v0.6.6 Strategy A full] amber14 atom type patch table for MTR.
# amber14 types (from amber14-all.xml AtomTypes block) — replaces GAFF2 types so
# bonded params (bond/angle/dihedral) at MTR-neighbor junction come from amber14,
# avoiding Khoury 2013 pure-GAFF2 backbone bias.
# CM + HMs use amber14 CT/HC (sp3 C methyl on aromatic N) so all bonds are amber14.
# Ref: Maier 2015 (doi:10.1021/acs.jctc.5b00255) TRP <Residue> ff14SB atom types.
_MTR_AMBER14_TYPES = {
    "N":   "protein-N",   "H":   "protein-H",
    "CA":  "protein-CX",  "HA":  "protein-H1",
    "CB":  "protein-CT",  "HB2": "protein-HC", "HB3": "protein-HC",
    "CG":  "protein-C*",
    "CD1": "protein-CW",  "HD1": "protein-H4",
    "NE1": "protein-NA",
    "CE2": "protein-CN",
    "CD2": "protein-CB",
    "CE3": "protein-CA",  "HE3": "protein-HA",
    "CZ3": "protein-CA",  "HZ3": "protein-HA",
    "CH2": "protein-CA",  "HH2": "protein-HA",
    "CZ2": "protein-CA",  "HZ2": "protein-HA",
    "C":   "protein-C",   "O":   "protein-O",
    # CM/HMs: use amber14 CT/HC so NE1(NA)-CM(CT) bond param is found in amber14.
    # amber14 has NA-CT (nucleic acid: ribose-base N-C bond) and CT-HC (alkyl).
    "CM":  "protein-CT",
    "HM1": "protein-HC",  "HM2": "protein-HC", "HM3": "protein-HC",
}


def _strip_spurious_ncaa_hydrogens(xml_path: str, ncaa_def) -> bool:
    """[v0.6.7] parent_residue 기반 ncAA 의 XML 에서 amber14 parent template 보다
    **많은** H 를 가진 heavy atom 의 남은 H 를 제거한다.

    Free-acid smiles_free 의 NH2 (2 H on backbone N) 는 parent internal form
    (amber14 TRP 등) 의 NH (1 H) 기대치를 초과한다. Walker rename 은 첫 H 만
    표준 이름 부여하고 나머지 H 는 GAFF 이름 (H8, H9 등) 으로 남음. 이 스퍼리어스
    H 는 downstream addHydrogens + MM-GBSA template match 에서 snap atom set 과
    XML template atom set 의 불일치를 초래한다 (Cp4 MTR H9 버그).

    Scan rules:
      - Heavy atom 의 expected H 수 = parent_h_name_table(parent)[heavy] 의 len
        + extension_atoms 의 해당 atom 의 H 예상 (CM: 3, OG/OH: 1 등)
      - 실제 H 수가 expected 초과이면 **이름이 GAFF-pattern** (Hdigit-only) 인
        H 부터 순차 제거.
    """
    parent_residue = getattr(ncaa_def, "parent_residue", None)
    if not parent_residue or not os.path.exists(xml_path):
        return False
    try:
        from parent_topology import parent_h_name_table  # type: ignore
    except Exception:
        return False

    import re as _re
    tree = ET.parse(xml_path)
    root = tree.getroot()
    modified = False
    pdb_resname = getattr(ncaa_def, "pdb_resname", "") or ""
    variant_names = {pdb_resname, f"N{pdb_resname}", f"C{pdb_resname}"}
    h_table = parent_h_name_table(parent_residue)
    extension_atoms = getattr(ncaa_def, "extension_atoms", ()) or ()

    for res in root.iter("Residue"):
        if res.get("name") not in variant_names:
            continue
        # Build atom_name -> atom_element map + H-to-heavy bond adjacency
        atoms_by_name = {a.get("name"): a for a in res.findall("Atom")}
        h_to_heavy: Dict[str, str] = {}
        heavy_to_hs: Dict[str, List[str]] = {}
        for b in res.findall("Bond"):
            a1, a2 = b.get("atomName1"), b.get("atomName2")
            if a1 is None or a2 is None:
                continue
            # H atoms start with 'H' (convention), heavy atoms don't
            is_h1 = a1.startswith("H")
            is_h2 = a2.startswith("H")
            if is_h1 and not is_h2:
                h_to_heavy[a1] = a2
                heavy_to_hs.setdefault(a2, []).append(a1)
            elif is_h2 and not is_h1:
                h_to_heavy[a2] = a1
                heavy_to_hs.setdefault(a1, []).append(a2)

        # Expected H count per heavy, from parent template + extensions.
        expected_h_count: Dict[str, int] = {}
        for heavy, h_names in h_table.items():
            expected_h_count[heavy] = len(h_names)
        # Extensions: heuristic — suffix length 1 (e.g., CM) → 3 H's (methyl);
        # suffix >1 (e.g., O1P) → 1 H (hydroxyl), unless element is not a
        # hydrogen-acceptor (F/Cl/Br → 0 H).
        for (ext_name, ext_elem, _attach, _bl) in extension_atoms:
            if ext_elem in ("F", "Cl", "Br", "I"):
                expected_h_count[ext_name] = 0
                continue
            if not ext_name:
                continue
            suffix = ext_name[1:]
            expected_h_count[ext_name] = 3 if len(suffix) == 1 else 1

        # For each heavy atom with excess H, drop GAFF-pattern-named H's first,
        # then any H beyond expected count.
        gaff_name_re = _re.compile(r"^H\d+[a-z]?$")  # H1, H12, H9b, etc.
        atoms_to_remove: List[str] = []
        for heavy, hs in heavy_to_hs.items():
            limit = expected_h_count.get(heavy)
            if limit is None:
                # Unknown heavy — leave intact
                continue
            excess = len(hs) - limit
            if excess <= 0:
                continue
            # Rank removal: GAFF-style first, then lexicographic
            gaff_style = [h for h in hs if gaff_name_re.match(h)]
            std_style = [h for h in hs if not gaff_name_re.match(h)]
            to_drop = (gaff_style + std_style)[:excess]
            atoms_to_remove.extend(to_drop)

        if not atoms_to_remove:
            continue

        drop_set = set(atoms_to_remove)
        for atom in list(res.findall("Atom")):
            if atom.get("name") in drop_set:
                res.remove(atom)
                modified = True
        for bond in list(res.findall("Bond")):
            if bond.get("atomName1") in drop_set or bond.get("atomName2") in drop_set:
                res.remove(bond)
        log.info(
            f"[spurious H strip] {res.get('name')}: removed {sorted(drop_set)} "
            f"(parent={parent_residue} over-H cleanup)"
        )

    if modified:
        tree.write(xml_path)
    return modified


def _apply_mtr_amber14_charge_patch(xml_path: str) -> bool:
    """[v0.6.6 Strategy A] Replace MTR atom charges with amber14SB TRP values.

    Fixes pure-GAFF2 backbone bias (Khoury 2013 §2.3): GAFF2 assigns terminal-amine
    charges (e.g., N=−0.89) to backbone N, causing NVT electrostatic explosion when
    placed adjacent to amber14 internal amides (−0.42).

    Side effect: removes spurious HX1 atom (2nd backbone H from capped SMILES's
    free-amine N-terminus, invalid for internal residue).

    Args:
        xml_path: MTR_gaff2.xml path (in-place modification).

    Returns:
        bool: True if any atom was patched.
    """
    if not os.path.exists(xml_path):
        return False
    tree = ET.parse(xml_path)
    root = tree.getroot()
    modified = False

    for res in root.iter("Residue"):
        if res.get("name") not in ("MTR", "NMTR", "CMTR"):
            continue
        # Remove spurious HX1 atom + its bonds (before charge patch)
        for atom in list(res.findall("Atom")):
            if atom.get("name") == "HX1":
                res.remove(atom)
                modified = True
                for bond in list(res.findall("Bond")):
                    if "HX1" in (bond.get("atomName1"), bond.get("atomName2")):
                        res.remove(bond)
        # Patch charges + types (Strategy A full)
        atom_names = set()
        for atom in res.findall("Atom"):
            name = atom.get("name")
            atom_names.add(name)
            if name in _MTR_AMBER14_CHARGES:
                old_q = atom.get("charge")
                new_q = _MTR_AMBER14_CHARGES[name]
                if abs(float(old_q) - new_q) > 1e-5:
                    atom.set("charge", f"{new_q:.4f}")
                    modified = True
            if name in _MTR_AMBER14_TYPES:
                old_type = atom.get("type")
                new_type = _MTR_AMBER14_TYPES[name]
                if old_type != new_type:
                    atom.set("type", new_type)
                    modified = True
        # Sanity: compute net charge sum for MTR variant
        total = sum(float(a.get("charge", 0)) for a in res.findall("Atom"))
        log.info(
            f"[amber14 patch] {res.get('name')}: net charge = {total:+.4f} "
            f"(expected ~0 for internal; HXT adds +0.4 for CMTR)"
        )

    # [v0.6.6 Strategy A] Inject cross-force-field bridge bonded params at CM-NE1 junction.
    # amber14 ff14SB has no protein-CT↔protein-NA bond (TRP NE1 is only bonded to
    # CW/CN in native Trp). Methylated NE1 (MTR) requires CT-NA bond param for
    # Modeller.createSystem. We use DNA.OL15 CT-N* values (same chemical motif:
    # sp3 methyl C to aromatic N-alkyl in nucleic acid methylation).
    # Ref: DNA.OL15 CT-N* (amber14/DNA.OL15.xml) — 5-methylcytosine methyl model.
    if modified:
        # Find or create HarmonicBondForce block
        bond_force = root.find("HarmonicBondForce")
        if bond_force is None:
            bond_force = ET.SubElement(root, "HarmonicBondForce")
        existing_bond_keys = {
            frozenset((b.get("type1"), b.get("type2")))
            for b in bond_force.findall("Bond")
        }
        bridge_bond = ("protein-CT", "protein-NA")
        if frozenset(bridge_bond) not in existing_bond_keys:
            ET.SubElement(bond_force, "Bond", {
                "type1": "protein-CT", "type2": "protein-NA",
                "length": "0.1475", "k": "282001.6",
            })
        # HarmonicAngleForce: HC-CT-NA, CT-NA-CW, CT-NA-CN
        angle_force = root.find("HarmonicAngleForce")
        if angle_force is None:
            angle_force = ET.SubElement(root, "HarmonicAngleForce")
        existing_angle_keys = {
            (a.get("type1"), a.get("type2"), a.get("type3"))
            for a in angle_force.findall("Angle")
        }
        bridge_angles = [
            # From DNA.OL15 H1-CT-N* (same tetrahedral sp3 geometry):
            ("protein-HC", "protein-CT", "protein-NA", "1.9111355309", "418.40"),
            # From DNA.OL15 CB-N*-CT: 5-ring junction (CW equivalent):
            ("protein-CT", "protein-NA", "protein-CW", "2.1956241990", "585.76"),
            # From DNA.OL15 CM-N*-CT: 5-6 ring junction (CN equivalent):
            ("protein-CT", "protein-NA", "protein-CN", "2.1153390534", "585.76"),
        ]
        for t1, t2, t3, angle_rad, k in bridge_angles:
            key = (t1, t2, t3)
            key_rev = (t3, t2, t1)
            if key not in existing_angle_keys and key_rev not in existing_angle_keys:
                ET.SubElement(angle_force, "Angle", {
                    "type1": t1, "type2": t2, "type3": t3,
                    "angle": angle_rad, "k": k,
                })
        tree.write(xml_path)
    return modified


def _build_pdb_name_map_trp_bootstrap(atoms) -> dict:
    """[v0.6.6 Strategy A] MTR (1-Me-Trp) 전용 GAFF2 → amber14 PDB 원자명 매핑.

    구현 배경 (Khoury 2014 §2.3 doi:10.1021/sb400168u — Forcefield_PTM
    ncAA atom-naming map; cf. Maier 2015 amber14 TRP scaffold
    doi:10.1021/acs.jctc.5b00255). 이전 Capece 2012 (jp2082825) 인용은
    오인용으로 2026-04-28 정정됨 — 자세한 내용은
    outputs/analysis/capece_2012_acquisition_log_20260428.md 참조.
    - MTR 은 amber14 TRP scaffold 에 N1-methyl 확장. 모든 backbone + indole
      heavy atom 은 amber14 TRP 의 표준 원자명을 갖고, CM 은 NE1 에 부착된
      methyl C. 이 함수는 GAFF2 mol2 의 atom graph 를 걸어 amber14 호환
      이름을 산출한다.
    - ncAA_mutate 에서 source TRP 의 sidechain 을 보존했고, CM 을 NE1 의
      ring-바깥쪽 bisector 에 배치했으므로, PDB 측은 amber14 네이밍 완료 상태.
      본 함수는 XML 측도 동일 네이밍으로 매핑해 OpenMM 템플릿 매칭 성공 보장.

    Args:
        atoms: parmed Atom 객체 리스트 (capped MTR mol2)

    Returns:
        dict: {gaff2_name: pdb_name} 매핑. 빈 dict → 매핑 실패 (수용 불가).
    """
    heavy = [a for a in atoms if a.atomic_number > 1]
    name_map: dict = {}

    def heavy_neighbors(atom):
        return [
            (b.atom1 if b.atom2 is atom else b.atom2)
            for b in atom.bonds
            if (b.atom1 if b.atom2 is atom else b.atom2).atomic_number > 1
        ]

    # 1. Backbone N: CA 이웃이 3개 heavy neighbor (N, C, CB) 를 갖고,
    #    그 중 하나는 O 를 이웃으로 갖는 (C=O) 를 찾는다.
    backbone_n = None
    for n in [a for a in heavy if a.atomic_number == 7]:
        for nb in heavy_neighbors(n):
            if nb.atomic_number != 6:
                continue
            nb_heavy = heavy_neighbors(nb)
            if len(nb_heavy) < 3:
                continue
            has_carbonyl = any(
                x.atomic_number == 6 and any(y.atomic_number == 8 for y in heavy_neighbors(x))
                for x in nb_heavy
            )
            if has_carbonyl:
                backbone_n = n
                name_map[n.name] = "N"
                name_map[nb.name] = "CA"
                ca_atom = nb
                break
        if backbone_n is not None:
            break
    if backbone_n is None:
        return {}

    ca_atom = next((a for a in heavy if a.name == name_map.get(next((n for n, v in name_map.items() if v == "CA"), ""))), None)
    # Re-resolve ca_atom via name_map
    ca_atom = None
    for a in heavy:
        if name_map.get(a.name) == "CA":
            ca_atom = a
            break
    if ca_atom is None:
        return {}

    # 2. CA 에서 C (carbonyl), CB 분리
    c_atom = cb_atom = None
    for nbr in heavy_neighbors(ca_atom):
        if nbr is backbone_n or nbr.name in name_map:
            continue
        if nbr.atomic_number != 6:
            continue
        if any(x.atomic_number == 8 for x in heavy_neighbors(nbr)):
            name_map[nbr.name] = "C"
            c_atom = nbr
        else:
            name_map[nbr.name] = "CB"
            cb_atom = nbr
    if c_atom is None or cb_atom is None:
        return {}

    # 3. C → O, OXT (free acid 용) 분리
    for o in heavy_neighbors(c_atom):
        if o.atomic_number != 8:
            continue
        has_h = any((b.atom1 if b.atom2 is o else b.atom2).atomic_number == 1 for b in o.bonds)
        name_map[o.name] = "OXT" if has_h else "O"

    # 4. CB → CG (sp2 aromatic C with 3 heavy neighbors)
    cg_atom = None
    for nbr in heavy_neighbors(cb_atom):
        if nbr is ca_atom or nbr.name in name_map:
            continue
        if nbr.atomic_number == 6 and len(heavy_neighbors(nbr)) >= 3:
            name_map[nbr.name] = "CG"
            cg_atom = nbr
            break
    if cg_atom is None:
        return {}

    # 5. CG 의 2개 aromatic C neighbor (CB 제외): CD1 (→ NE1 N), CD2 (→ aromatic C 링 junction)
    cd1_atom = cd2_atom = None
    ne1_atom = None
    for nbr in heavy_neighbors(cg_atom):
        if nbr is cb_atom or nbr.name in name_map:
            continue
        if nbr.atomic_number != 6:
            continue
        nbr_heavy_N = [x for x in heavy_neighbors(nbr) if x.atomic_number == 7]
        if nbr_heavy_N:
            cd1_atom = nbr
            ne1_atom = nbr_heavy_N[0]
            name_map[nbr.name] = "CD1"
            name_map[ne1_atom.name] = "NE1"
        else:
            cd2_atom = nbr
            name_map[nbr.name] = "CD2"
    if cd1_atom is None or cd2_atom is None or ne1_atom is None:
        return {}

    # 6. NE1 의 나머지 heavy neighbor: CE2 (5-ring 완성) + CM (methyl — MTR 전용)
    ce2_atom = cm_atom = None
    for nbr in heavy_neighbors(ne1_atom):
        if nbr is cd1_atom or nbr.name in name_map:
            continue
        if nbr.atomic_number != 6:
            continue
        nbr_heavy = heavy_neighbors(nbr)
        if len(nbr_heavy) == 1:
            # only NE1 as heavy neighbor → sp3 CH3 methyl
            cm_atom = nbr
            name_map[nbr.name] = "CM"
        else:
            # aromatic ring junction
            ce2_atom = nbr
            name_map[nbr.name] = "CE2"
    if ce2_atom is None:
        return {}

    # 7. CE2 → CZ2 (aromatic CH in 6-ring, neighbor not in {NE1, CD2})
    for nbr in heavy_neighbors(ce2_atom):
        if nbr is ne1_atom or nbr is cd2_atom or nbr.name in name_map:
            continue
        if nbr.atomic_number == 6:
            name_map[nbr.name] = "CZ2"

    # 8. CD2 → CE3 (aromatic CH in 6-ring, neighbor not in {CG, CE2})
    ce3_atom = None
    for nbr in heavy_neighbors(cd2_atom):
        if nbr is cg_atom or nbr is ce2_atom or nbr.name in name_map:
            continue
        if nbr.atomic_number == 6:
            name_map[nbr.name] = "CE3"
            ce3_atom = nbr

    # 9. CE3 → CZ3
    cz3_atom = None
    if ce3_atom is not None:
        for nbr in heavy_neighbors(ce3_atom):
            if nbr is cd2_atom or nbr.name in name_map:
                continue
            if nbr.atomic_number == 6:
                name_map[nbr.name] = "CZ3"
                cz3_atom = nbr

    # 10. CZ3 → CH2
    if cz3_atom is not None:
        for nbr in heavy_neighbors(cz3_atom):
            if nbr is ce3_atom or nbr.name in name_map:
                continue
            if nbr.atomic_number == 6:
                name_map[nbr.name] = "CH2"

    # 11. H 원자 명명: 부모 heavy atom PDB 이름 기반.
    H_NAME_TABLE = {
        "N":   ["H"],                 # amide H (backbone)
        "CA":  ["HA"],
        "CB":  ["HB2", "HB3"],
        "CD1": ["HD1"],
        "CM":  ["HM1", "HM2", "HM3"],
        "CE3": ["HE3"],
        "CZ2": ["HZ2"],
        "CZ3": ["HZ3"],
        "CH2": ["HH2"],
        "OXT": ["HXT"],
    }
    h_idx: dict = {}
    for atom in atoms:
        if atom.atomic_number != 1:
            continue
        parents = [b.atom1 if b.atom2 is atom else b.atom2 for b in atom.bonds]
        if not parents:
            continue
        pname = name_map.get(parents[0].name)
        if pname in H_NAME_TABLE:
            idx = h_idx.get(pname, 0)
            tbl = H_NAME_TABLE[pname]
            name_map[atom.name] = tbl[idx] if idx < len(tbl) else f"HX{idx}"
            h_idx[pname] = idx + 1

    return name_map


def _apply_name_map_to_xml(xml_path: str, name_map: dict) -> bool:
    """[v3 5-1] XML 파일 내 원자명을 name_map 에 따라 일괄 치환한다.

    Args:
        xml_path: 대상 OpenMM XML 파일 경로 (in-place 수정)
        name_map: {gaff2_name: pdb_name} 매핑. 빈 dict 면 즉시 False 반환.

    Returns:
        bool: 실제 파일 수정이 수행되었으면 True, 아니면 False.
    """
    if not name_map:
        return False
    with open(xml_path, encoding="utf-8") as f:
        content = f.read()
    for old, new in sorted(name_map.items(), key=lambda x: -len(x[0])):
        content = re.sub(
            r'((?:name|atomName1|atomName2|atomName)\s*=\s*")' + re.escape(old) + r'(")',
            r'\g<1>' + new + r'\g<2>',
            content,
        )
    with open(xml_path, "w", encoding="utf-8") as f:
        f.write(content)
    return True


def _rename_xml_to_pdb_names(xml_path: str, mol2_path: str, ncaa_code: str,
                             ncaa_mutation_type: str, ncaa_def=None) -> bool:
    """[v3 5-1 | v0.6.6/v0.6.7 ext] _build_pdb_name_map + _apply_name_map_to_xml 의 얇은 오케스트레이터.

    Dispatch 우선순위:
    1. ncaa_def.parent_residue 지정 시 → generic amber14 template-driven walker
       (utils/parent_topology.build_pdb_name_map_via_parent). 모든 canonical
       amino acid parent 자동 지원.
    2. mutation_type == "N-M" (sidechain 없는 ncAA, 예: NMA) → legacy walker.
    3. 그 외 → no-op (rename skip).
    """
    if not HAS_PARMED:
        return False
    parent_residue = getattr(ncaa_def, "parent_residue", None) if ncaa_def else None
    extension_defs = getattr(ncaa_def, "extension_atoms", ()) if ncaa_def else ()
    if ncaa_mutation_type != "N-M" and parent_residue is None:
        return False
    try:
        mol2_obj = pmd.load_file(mol2_path)
        atoms = list(list(mol2_obj.to_library().values())[0].atoms) if isinstance(mol2_obj, pmd.modeller.ResidueTemplateContainer) else list(mol2_obj.residues[0].atoms) if hasattr(mol2_obj, "residues") else list(mol2_obj.atoms)
        if parent_residue:
            # v0.6.7: generic parent walker (any canonical amino acid parent)
            from parent_topology import build_pdb_name_map_via_parent
            name_map = build_pdb_name_map_via_parent(
                atoms, parent_residue, extension_defs or (),
            )
            if not name_map:
                log.warning(
                    f"[XML rename] generic parent walker failed for {ncaa_code} "
                    f"(parent={parent_residue}). Falling back to legacy N-M walker."
                )
                name_map = _build_pdb_name_map(atoms, ncaa_mutation_type)
        else:
            name_map = _build_pdb_name_map(atoms, ncaa_mutation_type)
        if not name_map:
            return False
        return _apply_name_map_to_xml(xml_path, name_map)
    except Exception as _rename_err:
        # [v35 AUDIT] Principle 5: the GAFF2→PDB atom-name rename is required
        # for the anchored graph-isomorphism gate in run_restrained_md.py to
        # find N/CA/C anchors. A silent False return caused downstream strict
        # matching to fail with a cryptic "anchor not found" message; now the
        # exception is logged so operators can trace rename pathology.
        log.warning(f"[XML rename] failed for {ncaa_code}: {_rename_err.__class__.__name__}: {_rename_err}")
        return False

def _generate_topology_variants(xml_path: str, res_name: str) -> bool:
    """[v5 해법 B] 단일 ncAA XML 안에 internal/N-term/C-term 세 변형 잔기
    템플릿을 동시에 정의한다 (계층적 ncAA 방어 전략의 두 번째 방어선).

    배경 (CLAUDE.md v5 §해법 B):
        해법 A 가 AUTO_SASA 후보에서 chain 말단을 제외하여 ncAA 의 말단 배치를
        예방하지만, 사용자가 명시 위치(예: ``--residues B5,B12``) 로 말단 잔기를
        지정하거나 동적 분석에서 ncAA 가 chain 끝에 배치되는 경우 OpenMM
        ``addHydrogens`` 가 적절한 terminal 템플릿을 찾지 못해 즉시 실패한다.
        본 함수는 amber14 의 표준 아미노산 처리 규약 (LYS / NLYS / CLYS) 을 모방
        하여, ParmEd 가 생성한 자유산(free-acid) 잔기 블록으로부터 두 개의
        terminal 변형 잔기 블록을 deep-copy 로 파생시켜 동일 XML 에 추가한다:

            {RES}  : Internal residue
                       - OXT/HXT 제거
                       - <ExternalBond atomName="N"/>
                       - <ExternalBond atomName="C"/>
                       - 본 함수가 아닌 ``_postprocess_xml_for_internal_residue`` 가 처리.
            N{RES} : N-terminal residue
                       - OXT/HXT 제거 (C 가 다음 잔기의 N 과 결합)
                       - <ExternalBond atomName="C"/> 만 존재
                       - 본 함수가 deep-copy 로 생성.
            C{RES} : C-terminal residue
                       - OXT/HXT 유지 (자유산 형태로 노출)
                       - <ExternalBond atomName="N"/> 만 존재
                       - 본 함수가 deep-copy 로 생성.

    동작 순서 (Principle 7 — 의존성 순서):
        본 함수는 반드시 ``_rename_xml_to_pdb_names`` 직후, ``_postprocess_xml_
        for_internal_residue`` *직전* 에 호출되어야 한다. 호출 시점에 원본
        ``{RES}`` 블록은 자유산 형태(OXT/HXT 포함, ExternalBond 부재) 이며,
        이 상태에서 deep-copy 해야 C-terminal 변형이 OXT/HXT 를 보존할 수 있다.
        후속 ``_postprocess_xml_for_internal_residue`` 는 ``{RES}`` 블록만
        필터링하여 수정하므로, 새로 추가된 ``N{RES}``/``C{RES}`` 는 영향받지
        않는다.

    데이터 신뢰성 (Principle 1):
        본 함수는 GAFF2 force-field 파라미터(전하, 결합, 각도) 를 일체 변경
        하지 않고 잔기 블록의 위상학적 (topological) 정의만 추가한다. 즉,
        ncAA 가 어느 위치에 배치되어도 동일한 RESP 전하 + GAFF2 파라미터 가
        적용되어 데이터 신뢰성이 유지된다 (변형 간 에너지 일관성).

    Idempotent:
        ``N{RES}`` 또는 ``C{RES}`` 블록이 이미 존재하면 해당 변형 생성을
        건너뛴다. 캐시 재실행 시 중복 생성을 방지한다.

    Args:
        xml_path: ParmEd → rename 직후 OpenMM ForceField XML 경로 (in-place 수정)
        res_name: 원본 잔기명 (예: "NMA")

    Returns:
        bool: 변형 잔기를 1개 이상 추가했으면 True, no-op 이면 False.
    """
    if not os.path.exists(xml_path):
        log.warning(f"[Topology Variants] 파일이 존재하지 않음: {xml_path}")
        return False

    REMOVE_NAMES = {"OXT", "HXT"}
    tree = ET.parse(xml_path)
    root = tree.getroot()

    residues_node = root.find("Residues")
    if residues_node is None:
        log.warning(f"[Topology Variants] <Residues> 컨테이너 미발견: {xml_path}")
        return False

    # 원본 internal/free-acid 잔기 블록 식별
    src_node = next(
        (r for r in residues_node.findall("Residue") if r.get("name") == res_name),
        None,
    )
    if src_node is None:
        log.warning(f"[Topology Variants] 원본 잔기 블록 미발견: name={res_name}")
        return False

    existing_names = {r.get("name") for r in residues_node.findall("Residue")}
    added_variants = []

    # ── N-terminal 변형: OXT/HXT 제거 + ExternalBond C 만 ──
    n_variant_name = f"N{res_name}"
    if n_variant_name not in existing_names:
        n_node = copy.deepcopy(src_node)
        n_node.set("name", n_variant_name)

        for atom in list(n_node.findall("Atom")):
            if atom.get("name") in REMOVE_NAMES:
                n_node.remove(atom)
        for bond in list(n_node.findall("Bond")):
            if bond.get("atomName1", "") in REMOVE_NAMES or bond.get("atomName2", "") in REMOVE_NAMES:
                n_node.remove(bond)
        for eb in list(n_node.findall("ExternalBond")):
            n_node.remove(eb)
        # N-terminal: C 만 다음 잔기와 결합
        present_atoms = {a.get("name") for a in n_node.findall("Atom")}
        if "C" in present_atoms:
            ET.SubElement(n_node, "ExternalBond", {"atomName": "C"})

        residues_node.append(n_node)
        added_variants.append(n_variant_name)

    # ── C-terminal 변형: OXT/HXT 유지 + ExternalBond N 만 ──
    c_variant_name = f"C{res_name}"
    if c_variant_name not in existing_names:
        c_node = copy.deepcopy(src_node)
        c_node.set("name", c_variant_name)

        # OXT/HXT 보존 (자유산 형태). 기존 ExternalBond 만 제거.
        for eb in list(c_node.findall("ExternalBond")):
            c_node.remove(eb)
        # C-terminal: N 만 이전 잔기와 결합
        present_atoms = {a.get("name") for a in c_node.findall("Atom")}
        if "N" in present_atoms:
            ET.SubElement(c_node, "ExternalBond", {"atomName": "N"})

        residues_node.append(c_node)
        added_variants.append(c_variant_name)

    if added_variants:
        tree.write(xml_path, xml_declaration=True, encoding="utf-8")
        log.info(
            f"[Topology Variants] {res_name}: 추가 변형 잔기 생성 = {added_variants}"
        )
        return True

    log.info(f"[Topology Variants] {res_name}: 모든 변형이 이미 존재 — no-op")
    return False


def _postprocess_xml_for_internal_residue(xml_path: str, res_name: str) -> bool:
    """[v36 CRITICAL FIX] 자유산(free acid) 형태의 GAFF2 XML 템플릿을 펩타이드 내부
    잔기(internal residue) 템플릿으로 변환한다.

    근본 원인 (CLAUDE.md 확정 진단):
        ncaa_registry 의 ``smiles_free`` 는 의도적으로 양 말단이 자유산 형태이며
        (RESP 전하 정확도 확보를 위해 -COOH 말단을 갖는다), 이 SMILES 로부터
        antechamber → tleap → ParmEd → OpenMM XML 변환을 거치면 결과 XML 의
        ``<Residue>`` 블록에 카르복실 산소(OXT)와 카르복실 수소(HXT)가 자유산
        잔기로 그대로 보존된다. 반면 ``ncaa_mutate.py`` 는 ncAA 를 펩타이드
        체인 *내부* 잔기로 삽입하므로 PDB 구조에는 OXT/HXT 가 존재하지 않는다.

        결과: XML 중원자 7 (N, CA, C, O, CB, CM, OXT) ≠ PDB 중원자 6 → OpenMM
        ``Modeller.addHydrogens()`` 의 템플릿 매칭이 실패하고 ``No template
        found for residue 167 (NMA). The set of atoms is similar to ASP``
        오류로 즉시 중단된다 (FAILED_ADDH.log 다수 사례 확인).

    해결 (Principle 1: 데이터 신뢰성 향상):
        1) ``<Residue name=res_name>`` 블록에서 OXT/HXT ``<Atom>`` 노드 제거
        2) OXT/HXT 를 참조하는 ``<Bond atomName1/atomName2>`` 노드 제거
        3) N(head) 및 C(tail) 에 ``<ExternalBond atomName="..."/>`` 가 누락된
           경우 보강 (이미 존재하면 idempotent 하게 no-op)

        본 함수는 잔기명에 의존하지 않으며, OXT/HXT 가 처음부터 없으면 ExternalBond
        보강만 수행하므로 모든 ncAA 에 안전하게 적용 가능하다 (Principle 7:
        의존성 일관성). 단, GAFF2→PDB 원자명 rename 이 선행되어 있어야 하므로
        호출부는 반드시 ``_rename_xml_to_pdb_names`` 직후에 호출한다.

        Principle 4 (Fail-Fast over Silent-Fallback): XML 파싱 실패는 명시적
        예외로 전파되며, 호출부의 ``frcmod_to_openmm_xml`` 가 캐치하지 않으면
        전체 파라미터화가 즉시 실패하여 운영자가 손상된 XML 을 즉각 인지할 수
        있다.

    Args:
        xml_path: ParmEd 가 생성한 OpenMM ForceField XML 의 경로 (in-place 수정)
        res_name: 후처리 대상 잔기명 (예: "NMA")

    Returns:
        bool: 실제 수정이 발생했으면 True, 변경이 없었으면 False.
    """
    if not os.path.exists(xml_path):
        log.warning(f"[XML 후처리] 파일이 존재하지 않음: {xml_path}")
        return False

    REMOVE_NAMES = {"OXT", "HXT"}
    tree = ET.parse(xml_path)
    root = tree.getroot()
    modified = False

    for res_node in root.iter("Residue"):
        if res_node.get("name") != res_name:
            continue

        # ── 1. 자유산 전용 원자(OXT, HXT) 제거 ──
        atoms_removed = []
        for atom in list(res_node.findall("Atom")):
            if atom.get("name") in REMOVE_NAMES:
                res_node.remove(atom)
                atoms_removed.append(atom.get("name"))

        # ── 2. 제거된 원자를 참조하는 결합 제거 (OpenMM XML 은 atomName1/2 사용) ──
        bonds_removed = 0
        for bond in list(res_node.findall("Bond")):
            a1 = bond.get("atomName1", "")
            a2 = bond.get("atomName2", "")
            if a1 in REMOVE_NAMES or a2 in REMOVE_NAMES:
                res_node.remove(bond)
                bonds_removed += 1

        # ── 3. N(head)/C(tail) ExternalBond 보강 (idempotent) ──
        existing_ext = {eb.get("atomName") for eb in res_node.findall("ExternalBond")}
        present_atoms = {a.get("name") for a in res_node.findall("Atom")}
        ext_added = []
        for anchor in ("N", "C"):
            if anchor in present_atoms and anchor not in existing_ext:
                ET.SubElement(res_node, "ExternalBond", {"atomName": anchor})
                ext_added.append(anchor)

        if atoms_removed or bonds_removed or ext_added:
            modified = True
            log.info(
                f"[XML 후처리] {res_name}: removed atoms={atoms_removed}, "
                f"removed bonds={bonds_removed}, added ExternalBond={ext_added}"
            )

    if modified:
        # encoding 명시하여 한글 path/메타데이터 안전성 확보 (Principle 7/v27 규약)
        tree.write(xml_path, xml_declaration=True, encoding="utf-8")
    return modified


def frcmod_to_openmm_xml(mol2_path, frcmod_path, output_xml, res_name, ncaa_mutation_type, resp_result=None, ncaa_def=None):
    if not HAS_PARMED:
        return {"HarmonicBondForce"}
    if resp_result:
        _inject_resp_into_mol2(mol2_path, resp_result)
    prmtop, inpcrd = _run_tleap(mol2_path, frcmod_path, os.path.dirname(mol2_path), res_name)
    struct = pmd.load_file(prmtop, inpcrd)
    ff = pmd.openmm.OpenMMParameterSet.from_structure(struct)
    tmpl = pmd.modeller.ResidueTemplate.from_residue(struct.residues[0])
    tmpl.name, tmpl.head, tmpl.tail = res_name, None, None
    # [v0.6.6 FIX] parent_residue 기반 ncAA (MTR 등) 는 residue 에 backbone 외
    # sidechain N 이 추가로 존재 (MTR: NE1). 단순 "첫 N" 선택은 NE1 을 head 로
    # 잘못 picking 하여 spurious ExternalBond 발생. parent-bootstrap rename map
    # 을 먼저 계산하여 backbone N/C 를 지정한다.
    parent_residue = getattr(ncaa_def, "parent_residue", None) if ncaa_def else None
    extension_defs = getattr(ncaa_def, "extension_atoms", ()) if ncaa_def else ()
    pre_name_map = {}
    if parent_residue:
        try:
            from parent_topology import build_pdb_name_map_via_parent
            pre_name_map = build_pdb_name_map_via_parent(
                list(tmpl.atoms), parent_residue, extension_defs or (),
            )
        except Exception as _e:
            log.warning(f"[head/tail] generic parent walker failed: {_e}")
            pre_name_map = {}
    for atom in tmpl.atoms:
        mapped = pre_name_map.get(atom.name)
        if mapped == "N":
            tmpl.head = atom
        elif mapped == "C":
            tmpl.tail = atom
    if tmpl.head is None and tmpl.tail is None:
        # Fall back to legacy first-N/carbonyl-C heuristic (non-bootstrap ncAAs)
        for atom in tmpl.atoms:
            if atom.atomic_number == 7 and tmpl.head is None:
                tmpl.head = atom
            elif atom.atomic_number == 6 and tmpl.tail is None and any(
                (b.atom1 if b.atom2 is atom else b.atom2).atomic_number == 8 for b in atom.bonds
            ):
                tmpl.tail = atom
    ff.residues[res_name] = tmpl
    ff.write(output_xml)

    # [v4 CLAUDE.md 방안 A] rename 실패를 silent 로 넘기지 말고 명시적으로 경고.
    # 실패 시 XML 내부 원자명 목록을 로그로 출력하여 downstream addHydrogens 실패
    # 진단을 돕는다. STANDARD mutation type 은 현재 rename 대상이 아니므로 (N-M 에만
    # 수행) 이 경로에서도 원자명 로그만 남겨 운영자가 XML 상태를 확인할 수 있게 한다.
    renamed = _rename_xml_to_pdb_names(output_xml, mol2_path, res_name, ncaa_mutation_type, ncaa_def=ncaa_def)
    if not renamed:
        log.warning(
            f"[XML Rename] GAFF2→PDB 원자명 변환 스킵 또는 실패 "
            f"(ncaa_code='{res_name}', mutation_type='{ncaa_mutation_type}'). "
            f"addHydrogens 단계에서 template 매칭 실패가 발생하면 XML 원자명을 점검하세요."
        )
        _log_xml_atom_names(output_xml, res_name)

    # [v5 해법 B] rename 직후, postprocess 직전 — 자유산 잔기 블록을 deep-copy
    # 하여 N-terminal/C-terminal 변형을 추가한다. C-terminal 변형은 OXT/HXT 가
    # 보존된 상태에서 복사되어야 하므로 반드시 postprocess 보다 먼저 실행한다.
    _generate_topology_variants(output_xml, res_name)

    # [v36 CRITICAL FIX] rename 직후 자유산→내부 잔기 후처리를 수행한다. rename 이
    # 성공한 경우에만 OXT/HXT 가 PDB 규약명으로 식별 가능하므로, rename skip 경로
    # (현재 STANDARD 타입) 에서는 후처리도 no-op 가 된다. 이 순서 의존성은
    # _postprocess_xml_for_internal_residue 의 docstring 에 명시되어 있다.
    # 본 호출은 {RES} 잔기명만 필터링하므로 위에서 추가된 N{RES}/C{RES} 는
    # 영향받지 않는다.
    _postprocess_xml_for_internal_residue(output_xml, res_name)

    # NOTE: amber14SB TRP type+charge patch is applied AFTER _emit_hydrogen_definitions
    # in main() — running it here would break hydrogen inference because
    # _emit_hydrogen_definitions resolves element via <AtomTypes> which doesn't contain
    # amber14 'protein-*' types. See main() for the correct ordering.

    return set()


def _log_xml_atom_names(xml_path: str, res_name: str) -> None:
    """[v4 CLAUDE.md 방안 A] XML 내 잔기 템플릿의 원자명을 진단 로그로 출력한다.

    rename 실패 경로 또는 미지원 mutation_type 경로에서 호출되어, 운영자가 XML
    원자명이 PDB 규약 (N/CA/C/O/CB/CM/...) 과 일치하는지 즉시 확인할 수 있게 한다.
    """
    try:
        tree = ET.parse(xml_path)
        for res_node in tree.getroot().iter("Residue"):
            if res_node.get("name") == res_name:
                atoms = [a.get("name") for a in res_node.findall("Atom")]
                heavy = [a for a in atoms if a and not a.startswith("H")]
                hs = [a for a in atoms if a and a.startswith("H")]
                log.warning(f"  [XML 진단] {res_name} heavy atoms: {heavy}")
                log.warning(f"  [XML 진단] {res_name} hydrogens:   {hs}")
                return
        log.warning(f"  [XML 진단] {res_name} 잔기 블록을 XML 에서 찾지 못함: {xml_path}")
    except (ET.ParseError, FileNotFoundError) as _err:
        log.warning(f"  [XML 진단] XML 읽기 실패 ({_err.__class__.__name__}): {_err}")

def _is_cache_valid(xml_out: str, manifest_out: str) -> bool:
    """파라미터화 캐시(XML + manifest 페어)의 유효성을 검사한다.

    XML과 manifest JSON이 모두 존재할 때에만 캐시가 유효하다고 판단한다.
    XML만 존재하는 부분 캐시는 stale/broken으로 간주하여 False를 반환하고,
    호출부가 재파라미터화를 수행하도록 위임한다 (run_restrained_md.py 의
    load_ncaa_manifest() Fail-Fast 경로와의 데이터 계약 일관성 보장).
    """
    if os.path.exists(xml_out) and os.path.exists(manifest_out):
        return True
    if os.path.exists(xml_out) and not os.path.exists(manifest_out):
        log.warning(f"[캐시 불완전] XML은 존재하나 manifest 누락 → 재파라미터화 강제: {xml_out}")
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputdir", required=True)
    parser.add_argument("--outputdir", required=True)
    parser.add_argument("--ncaa_code", required=True)
    args = parser.parse_args()

    os.makedirs(args.outputdir, exist_ok=True)

    ncaa_def = resolve_ncaa_definition(args.ncaa_code)
    if ncaa_def is None:
        log.error(f"[!] 알 수 없는 ncAA Code: {args.ncaa_code}")
        sys.exit(1)

    res_name = ncaa_def.xml_resname
    work_dir = os.path.join(args.outputdir, f"_work_{res_name}")
    os.makedirs(work_dir, exist_ok=True)

    xml_out = os.path.join(args.outputdir, f"{res_name}_gaff2.xml")
    manifest_out = os.path.join(args.outputdir, f"{args.ncaa_code}_params_manifest.json")
    # [v35 AUDIT] Principle 7: cache-hit must ALSO verify the manifest exists.
    # 검증 로직은 _is_cache_valid()로 분리되어 main()의 본 흐름을 단순화한다.
    if _is_cache_valid(xml_out, manifest_out):
        log.info(f"[캐시 발견] 기존 XML + manifest 재사용: {xml_out}")
        return

    mol_free = os.path.join(work_dir, f"{res_name}_free.mol")
    smiles_to_mol(ncaa_def.smiles_free, mol_free, res_name)

    mol2, frcmod = run_antechamber(mol_free, work_dir, res_name, ncaa_def.formal_charge)
    resp_result = None
    if HAS_PYSCF:
        xyz_free = os.path.join(work_dir, f"{res_name}_free.xyz")
        _mol_to_xyz(mol_free, xyz_free)
        resp_result = run_pyscf_resp(xyz_free, ncaa_def.formal_charge)

    frcmod_to_openmm_xml(mol2, frcmod, xml_out, res_name, ncaa_def.mutation_type, resp_result, ncaa_def=ncaa_def)

    # [v0.6.7] parent_residue ncAA 의 XML 에서 free-acid SMILES 기반 GAFF 생성물의
    # 스퍼리어스 H (backbone N 의 2번째 H, 또는 rename 안 된 GAFF-이름 H) 를 제거.
    # 이 cleanup 은 opt-in amber14 patch 와 무관하게 항상 실행되어야 함 —
    # downstream MM-GBSA template match 에서 snap atom set 과 XML atom set 의
    # 1:1 매칭을 강제.
    _strip_spurious_ncaa_hydrogens(xml_out, ncaa_def)

    # [v36 CRITICAL FIX-2] Hydrogen Definitions XML 생성 — addHydrogens 가
    # ncAA 잔기에 수소를 추가하기 위해 필요. 본 파일이 부재하면 downstream
    # run_restrained_md.py 의 addHydrogens 단계가 6-원자 NMA → 14-원자 템플릿
    # mismatch 로 createSystem 에서 실패한다.
    hyd_out = os.path.join(args.outputdir, f"{res_name}_hydrogens.xml")
    _emit_hydrogen_definitions(xml_out, res_name, hyd_out)

    # [v0.6.6 Strategy A] amber14SB TRP type+charge patch for MTR (parent_residue=TRP).
    # 핵심: hydrogen definitions 생성 뒤에 실행해야 함 (_emit_hydrogen_definitions
    # 는 <AtomTypes> 블록으로 element 추론하는데 amber14 'protein-*' type 은
    # MTR_gaff2.xml 에 정의되지 않아 lookup 실패 → H entries 누락).
    # Fixes Khoury 2013 §2.3 pure-GAFF2 backbone bias.
    # Refs: Maier 2015 (ff14SB TRP charge set, doi:10.1021/acs.jctc.5b00255);
    # Khoury 2013 / 2014 (Forcefield_PTM ncAA RESP-A2 protocol,
    # doi:10.1021/ct400556v / doi:10.1021/sb400168u); Bayly 1993
    # (RESP charge-equivalence, doi:10.1021/j100142a004); Cieplak 1995
    # (multi-conformer RESP for biopolymers, doi:10.1002/jcc.540161106).
    # Earlier "Capece 2012 doi:10.1021/jp2082825" inline ref was a
    # misattribution (heme-dioxygenase reaction-mechanism paper, not FF
    # methodology) — corrected 2026-04-28; see acquisition log
    # outputs/analysis/capece_2012_acquisition_log_20260428.md.
    parent_residue = getattr(ncaa_def, "parent_residue", None)
    # [v0.6.6 Strategy A; comment refresh 2026-05-11] amber14 type+charge patch is gated by
    # UPDD_MTR_AMBER14_PATCH; default value is "1" (patch ON), so the Trp-derived MTR class
    # XMLs (24/24 production set) ship with q_N = -0.4157 e (amber14SB Trp reference) and
    # |Σq| ≤ 1×10⁻⁶ e by default. Setting UPDD_MTR_AMBER14_PATCH=0 falls back to pure GAFF2
    # backbone-N (q_N = -0.8938 e, residual Σq = -0.187 e), retained as the audit baseline
    # in outputs/_archive/pre_amb14_patch_20260427/.
    if parent_residue == "TRP" and os.environ.get("UPDD_MTR_AMBER14_PATCH", "1") != "0":
        patched = _apply_mtr_amber14_charge_patch(xml_out)
        if patched:
            log.info(f"[amber14 patch] {res_name}: amber14SB TRP type+charge overlay applied (opt-in)")
        # hydrogens.xml 도 HX1 제거 (MTR template 에서 제거된 spurious 2nd amide H)
        if os.path.exists(hyd_out):
            hyd_tree = ET.parse(hyd_out)
            hyd_root = hyd_tree.getroot()
            removed = 0
            for res in hyd_root.findall("Residue"):
                for h in list(res.findall("H")):
                    if h.get("name") == "HX1":
                        res.remove(h)
                        removed += 1
            if removed:
                hyd_tree.write(hyd_out, xml_declaration=True, encoding="utf-8")
                log.info(f"[amber14 patch] {res_name}_hydrogens.xml: HX1 × {removed} 제거")

    write_manifest(
        args.outputdir, args, ncaa_def, xml_out,
        bool(resp_result), mol2, frcmod,
        hydrogens_path=hyd_out if os.path.exists(hyd_out) else "",
    )

if __name__ == "__main__":
    main()
