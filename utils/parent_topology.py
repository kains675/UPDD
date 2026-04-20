"""Generic amber14 parent-residue topology walker for ncAA parameterization.

Replaces hardcoded per-parent walkers (_build_pdb_name_map_trp_bootstrap etc.)
with a single `build_pdb_name_map_via_parent(atoms, parent_resname, extension_atoms)`
that works for any canonical amino acid declared as `parent_residue` in NCAADef.

Approach:
  1. Parse amber14 protein.ff14SB.xml for the parent residue's internal template
     (excludes OXT/HXT — only internal form used for ncAA substitution).
  2. Build heavy-atom graph from parent template (amber atom names as node IDs).
  3. Identify ACE/NME caps + extension atoms in GAFF2 mol2 via structural
     fingerprint (no hardcoded names).
  4. Build stripped GAFF2 heavy-atom subgraph (caps + extensions removed).
  5. VF2 graph isomorphism match (parent vs stripped-GAFF) with element equality.
  6. Return mapping {gaff_name: amber_name} + H naming via parent template H's.

Derived helpers also expose:
  - `parent_heavy_bonds(parent_resname)` — list of (amber_name, amber_name) for
    CONECT records emission in ncaa_mutate.
  - `parent_neighbor_hint(parent_resname, attach_name)` — 2 heavy neighbors of
    an attach atom, for sp2 bisector extension placement.
  - `parent_h_name_table(parent_resname)` — {heavy_name: [h_names...]} from
    amber14 template, for hydrogen naming in XML rename.

All data read from `/openmm/app/data/amber14/protein.ff14SB.xml` verbatim —
no hardcoded atom names, ensuring any future amber update propagates.

Ref:
  - Maier 2015 ff14SB (doi:10.1021/acs.jctc.5b00255)
  - Capece 2012 Strategy A hybrid pattern (doi:10.1021/jp2082825)
"""
from __future__ import annotations
import os
import functools
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Set, Tuple

try:
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher
    HAS_NX = True
except ImportError:
    HAS_NX = False

_ELEMENT_BY_Z = {
    1: "H", 6: "C", 7: "N", 8: "O", 9: "F",
    15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I",
}


def _amber14_ff14SB_path() -> str:
    """Locate amber14/protein.ff14SB.xml in the OpenMM install."""
    import openmm
    return os.path.join(
        os.path.dirname(openmm.__file__),
        "app", "data", "amber14", "protein.ff14SB.xml",
    )


@functools.lru_cache(maxsize=32)
def _load_parent_template(parent_resname: str) -> Optional[Tuple[Dict[str, str], List[Tuple[str, str]]]]:
    """Parse amber14 ff14SB.xml for a canonical residue template.

    Returns (atoms: {name → element}, bonds: [(name1, name2), ...]) for the
    internal form (OXT/HXT excluded). None if residue not found.
    """
    ff_path = _amber14_ff14SB_path()
    if not os.path.exists(ff_path):
        return None
    tree = ET.parse(ff_path)
    root = tree.getroot()
    type_to_elem = {t.get("name"): t.get("element") for t in root.findall("AtomTypes/Type")}
    for r in root.iter("Residue"):
        if r.get("name") != parent_resname:
            continue
        atoms = {}
        for a in r.findall("Atom"):
            name = a.get("name")
            # Exclude terminal-only atoms from the internal form
            if name in ("OXT", "HXT"):
                continue
            elem = type_to_elem.get(a.get("type"))
            atoms[name] = elem
        bonds = []
        for b in r.findall("Bond"):
            a1, a2 = b.get("atomName1"), b.get("atomName2")
            if a1 in atoms and a2 in atoms:
                bonds.append((a1, a2))
        return atoms, bonds
    return None


@functools.lru_cache(maxsize=32)
def parent_heavy_bonds(parent_resname: str) -> List[Tuple[str, str]]:
    """Return list of (heavy_name1, heavy_name2) bonds for the parent residue.

    Used by ncaa_mutate to emit CONECT records for non-standard residues so
    run_restrained_md.py graph_policy=strict finds explicit bonds.
    """
    tmpl = _load_parent_template(parent_resname)
    if tmpl is None:
        return []
    atoms, bonds = tmpl
    heavy_bonds = []
    for a1, a2 in bonds:
        if atoms.get(a1) != "H" and atoms.get(a2) != "H":
            heavy_bonds.append((a1, a2))
    return heavy_bonds


@functools.lru_cache(maxsize=128)
def parent_neighbor_hint(parent_resname: str, attach_name: str) -> Tuple[str, ...]:
    """Return up to 2 heavy neighbors of attach_name in parent_resname template.

    Used by ncaa_mutate for sp2 bisector extension placement. Sp2 aromatic N
    (e.g., TRP NE1) has 2 aromatic C neighbors — their bisector (opposite side)
    is the natural placement for a methyl extension (Capece 2012 §2).
    """
    tmpl = _load_parent_template(parent_resname)
    if tmpl is None:
        return ()
    atoms, bonds = tmpl
    neighbors = []
    for a1, a2 in bonds:
        if a1 == attach_name and atoms.get(a2) != "H":
            neighbors.append(a2)
        elif a2 == attach_name and atoms.get(a1) != "H":
            neighbors.append(a1)
    return tuple(neighbors[:2])


@functools.lru_cache(maxsize=32)
def parent_h_name_table(parent_resname: str) -> Dict[str, Tuple[str, ...]]:
    """Return {heavy_name: (h_name1, h_name2, ...)} from amber14 template.

    Groups the H atoms declared in amber14 ff14SB.xml by their bonded heavy
    atom. Used for H naming after graph-match of heavy atoms.
    """
    tmpl = _load_parent_template(parent_resname)
    if tmpl is None:
        return {}
    atoms, bonds = tmpl
    table: Dict[str, List[str]] = {}
    for a1, a2 in bonds:
        e1, e2 = atoms.get(a1), atoms.get(a2)
        if e1 == "H" and e2 != "H":
            table.setdefault(a2, []).append(a1)
        elif e2 == "H" and e1 != "H":
            table.setdefault(a1, []).append(a2)
    return {k: tuple(v) for k, v in table.items()}


def _heavy_neighbors(atom):
    """parmed helper: list of heavy-atom neighbors."""
    return [
        (b.atom1 if b.atom2 is atom else b.atom2)
        for b in atom.bonds
        if (b.atom1 if b.atom2 is atom else b.atom2).atomic_number > 1
    ]


def _identify_caps_and_backbone(atoms) -> Tuple[Set[str], Set[str], Optional[object], Optional[object]]:
    """Identify terminal cap atoms to exclude from the backbone heavy graph.

    Handles three forms:
      1. Capped (ACE-ncAA-NME): ACE 3 heavy + NME 2 heavy
      2. Free-acid (H2N-ncAA-COOH): N-terminal NH2 (kept) + C-terminal OXT (excluded)
      3. Mixed: N-cap ACE + C free-acid OXT, or N-terminal free amine + NME

    Heuristic: find backbone N via CA candidate (sp3 C with ≥3 heavy neighbors
    including C=O) + optional ACE cap. Then find backbone C via CA traversal
    and detect its terminal form: if it has 2 O neighbors (one =O, one -OH)
    → free-acid (mark OXT for exclusion). If it has 1 O + 1 N with methyl
    attached → NME cap. If only 1 O and nothing else → already internal (no cap).
    """
    heavy = [a for a in atoms if a.atomic_number > 1]
    backbone_n = None
    ace_c = None
    ca_atom = None

    # Pass 1: find CA first (sp3 C with ≥3 heavy neighbors where ≥2 are C and
    # at least one connected C has an O neighbor via backbone carbonyl).
    for c in [a for a in heavy if a.atomic_number == 6]:
        c_nbrs = _heavy_neighbors(c)
        if len(c_nbrs) < 3:
            continue
        # CA must have at least 1 N neighbor (backbone N or free amine) and
        # at least 1 C neighbor that is a carbonyl C (has =O).
        n_nbrs = [x for x in c_nbrs if x.atomic_number == 7]
        co_nbrs = [
            x for x in c_nbrs
            if x.atomic_number == 6 and any(y.atomic_number == 8 for y in _heavy_neighbors(x))
        ]
        if not n_nbrs or not co_nbrs:
            continue
        # Take this CA
        ca_atom = c
        backbone_n = n_nbrs[0]
        # Identify if backbone_n has an ACE cap (neighbor C with =O and methyl)
        for c2 in _heavy_neighbors(backbone_n):
            if c2 is ca_atom or c2.atomic_number != 6:
                continue
            c2_inner = _heavy_neighbors(c2)
            has_o = any(x.atomic_number == 8 for x in c2_inner)
            has_methyl = any(
                x.atomic_number == 6 and len(_heavy_neighbors(x)) == 1
                for x in c2_inner
            )
            if has_o and has_methyl:
                ace_c = c2
                break
        break

    ace_atoms: Set[str] = set()
    nme_atoms: Set[str] = set()
    backbone_c = None

    if ace_c is not None:
        ace_atoms.add(ace_c.name)
        for x in _heavy_neighbors(ace_c):
            if x is backbone_n:
                continue
            ace_atoms.add(x.name)

    # Find backbone C (CA's carbonyl neighbor)
    if ca_atom is not None:
        for c in _heavy_neighbors(ca_atom):
            if c is backbone_n or c.atomic_number != 6:
                continue
            c_inner = _heavy_neighbors(c)
            if any(x.atomic_number == 8 for x in c_inner):
                backbone_c = c
                break

    # Identify C-terminal cap: either NME (amide N with methyl) or free-acid (OXT).
    # In free-acid form, backbone C has 2 O neighbors: =O (internal) and -OH (OXT).
    # Strip only the hydroxyl-bearing oxygen (OXT) — keep the carbonyl =O.
    if backbone_c is not None:
        c_inner = _heavy_neighbors(backbone_c)
        o_neighbors = [x for x in c_inner if x.atomic_number == 8]
        n_neighbors = [x for x in c_inner if x.atomic_number == 7]
        # NME cap branch
        for n in n_neighbors:
            if n is backbone_n:
                continue
            nme_atoms.add(n.name)
            for x in _heavy_neighbors(n):
                if x is backbone_c:
                    continue
                nme_atoms.add(x.name)
            break
        # Free-acid branch: exclude OXT (the O with H attached, i.e., hydroxyl)
        if len(o_neighbors) >= 2:
            # Pick the O whose bonds include an H
            for o in o_neighbors:
                has_h = any(
                    (b.atom1 if b.atom2 is o else b.atom2).atomic_number == 1
                    for b in o.bonds
                )
                if has_h:
                    nme_atoms.add(o.name)
                    break

    return ace_atoms, nme_atoms, backbone_n, backbone_c


def _parent_heavy_graph(parent_resname: str):
    """Build networkx graph of parent residue heavy atoms."""
    if not HAS_NX:
        return None
    tmpl = _load_parent_template(parent_resname)
    if tmpl is None:
        return None
    atoms, bonds = tmpl
    g = nx.Graph()
    for name, elem in atoms.items():
        if elem == "H":
            continue
        g.add_node(name, atom_name=name, element=elem)
    for a1, a2 in bonds:
        if atoms.get(a1) != "H" and atoms.get(a2) != "H":
            g.add_edge(a1, a2)
    return g


def _gaff_heavy_subgraph(atoms, extension_atom_names: Set[str], cap_atom_names: Set[str]):
    """Build networkx graph of ncAA GAFF2 heavy atoms, stripping caps + extensions."""
    if not HAS_NX:
        return None
    excluded = cap_atom_names | set(extension_atom_names)
    g = nx.Graph()
    for a in atoms:
        if a.atomic_number <= 1 or a.name in excluded:
            continue
        elem = _ELEMENT_BY_Z.get(a.atomic_number, "?")
        g.add_node(a.name, atom_name=a.name, element=elem)
    for a in atoms:
        if a.atomic_number <= 1 or a.name in excluded:
            continue
        for b in a.bonds:
            other = b.atom2 if b.atom1 is a else b.atom1
            if other.atomic_number <= 1 or other.name in excluded:
                continue
            if a.name == other.name:
                continue
            g.add_edge(a.name, other.name)
    return g


def _find_extensions_structurally(atoms, parent_resname: str, cap_atoms: Set[str],
                                   extension_defs: Tuple[Tuple[str, str, str, float], ...]):
    """Identify GAFF2 atom names corresponding to declared extension atoms.

    Uses iterative "remove + check isomorphism" strategy: try excluding each
    combination of candidate methyl-like atoms and check if the remaining
    heavy graph matches the parent template. The winning combination is the
    extension set.

    Args:
        atoms: parmed Atom list.
        parent_resname: amber14 parent code.
        cap_atoms: set of ACE/NME cap atom names (to exclude from candidates).
        extension_defs: tuple of (name, element, attach_amber_name, bond_A).

    Returns:
        dict {gaff_name: declared_pdb_name} for each matched extension atom.
        Empty dict if no combination produces a valid match.
    """
    if not extension_defs or not HAS_NX:
        return {}
    n_ext = len(extension_defs)
    # Candidates: heavy atoms with exactly 1 heavy neighbor (methyl-like terminal),
    # not in caps, whose element matches one of the declared extension elements.
    ext_elements = [_elem.upper() for (_, _elem, *_rest) in extension_defs]
    candidates = []
    for a in atoms:
        if a.atomic_number <= 1 or a.name in cap_atoms:
            continue
        elem = _ELEMENT_BY_Z.get(a.atomic_number, "?")
        if elem.upper() not in ext_elements:
            continue
        heavy_nbrs = _heavy_neighbors(a)
        if len(heavy_nbrs) != 1:
            continue
        candidates.append(a)

    if len(candidates) < n_ext:
        return {}

    from itertools import combinations
    parent_g = _parent_heavy_graph(parent_resname)
    if parent_g is None:
        return {}

    for combo in combinations(candidates, n_ext):
        ext_names = set(a.name for a in combo)
        gaff_g = _gaff_heavy_subgraph(atoms, ext_names, cap_atoms)
        if gaff_g is None or parent_g.number_of_nodes() != gaff_g.number_of_nodes():
            continue
        gm = GraphMatcher(
            parent_g, gaff_g,
            node_match=lambda p, g: p.get("element") == g.get("element"),
        )
        if not gm.is_isomorphic():
            continue
        # Good combo — now match each GAFF extension atom to its declared name
        # via its attach atom (declared_pdb_name's amber attach atom in parent).
        iso = next(gm.isomorphisms_iter())  # parent → gaff
        amber_to_gaff = dict(iso)
        result: Dict[str, str] = {}
        for ext_atom in combo:
            # ext_atom's single heavy neighbor's gaff_name maps back to amber
            attach_gaff = _heavy_neighbors(ext_atom)[0].name
            attach_amber = None
            for amber_n, gaff_n in amber_to_gaff.items():
                if gaff_n == attach_gaff:
                    attach_amber = amber_n
                    break
            # Find which extension_def has this attach_amber
            for (ext_name, ext_elem, attach_name, _bond) in extension_defs:
                if attach_name == attach_amber and ext_name not in result.values():
                    result[ext_atom.name] = ext_name
                    break
        if len(result) == n_ext:
            return result
    return {}


def build_pdb_name_map_via_parent(
    atoms,
    parent_resname: str,
    extension_defs: Tuple[Tuple[str, str, str, float], ...] = (),
) -> Dict[str, str]:
    """Generic parent-residue topology walker.

    Maps GAFF2 mol2 atom names → amber14 PDB convention names for a ncAA whose
    scaffold is a canonical amino acid (parent_residue), with optional
    extension atoms (e.g., MTR's CM methyl on NE1) identified structurally.

    Args:
        atoms: parmed Atom list from the capped ncAA mol2 (ACE-ncAA-NME).
        parent_resname: amber14 canonical residue code (e.g., "LEU", "PHE").
        extension_defs: NCAADef.extension_atoms tuple of
            (name, element, attach_amber_name, bond_length_A).

    Returns:
        {gaff_name: amber_name} mapping including heavy + H atoms + extensions.
        Empty dict on match failure.
    """
    if not HAS_NX:
        return {}

    ace_atoms, nme_atoms, _n, _c = _identify_caps_and_backbone(atoms)
    cap_atoms = ace_atoms | nme_atoms

    # Structurally identify extension atoms by trying removal combinations
    ext_map = _find_extensions_structurally(atoms, parent_resname, cap_atoms, extension_defs)
    extension_gaff_names = set(ext_map.keys())

    parent_g = _parent_heavy_graph(parent_resname)
    if parent_g is None:
        return {}
    gaff_g = _gaff_heavy_subgraph(atoms, extension_gaff_names, cap_atoms)
    if gaff_g is None:
        return {}

    if parent_g.number_of_nodes() != gaff_g.number_of_nodes():
        return {}

    gm = GraphMatcher(
        parent_g, gaff_g,
        node_match=lambda p, g: p.get("element") == g.get("element"),
    )
    if not gm.is_isomorphic():
        return {}
    mapping = next(gm.isomorphisms_iter())  # {parent_name → gaff_name}

    name_map: Dict[str, str] = {}
    for parent_name, gaff_name in mapping.items():
        name_map[gaff_name] = parent_name
    # Extension atoms: keep their declared PDB names
    for gaff_name, pdb_name in ext_map.items():
        name_map[gaff_name] = pdb_name
    # Free-acid OXT: detected in nme_atoms for free-acid form (the O with H).
    # Rename it to OXT so downstream _postprocess_xml_for_internal_residue finds
    # it by name for removal + ExternalBond emission.
    if nme_atoms and not ace_atoms:
        # Free-acid form (no ACE): the single atom in nme_atoms is OXT
        for oxt_gaff_name in nme_atoms:
            # Verify it's an O (not NME N or methyl C)
            for a in atoms:
                if a.name == oxt_gaff_name and a.atomic_number == 8:
                    name_map[oxt_gaff_name] = "OXT"
                    # Its H is HXT
                    for b in a.bonds:
                        other = b.atom2 if b.atom1 is a else b.atom1
                        if other.atomic_number == 1:
                            name_map[other.name] = "HXT"
                    break

    # Hydrogens via amber14 template H-to-parent table
    h_table = parent_h_name_table(parent_resname)
    # Extension H naming convention: H{ext_name_letter}1, H{letter}2, ...
    # e.g., CM → HM1/HM2/HM3 (MTR pattern, Capece 2012)
    ext_h_tables: Dict[str, Tuple[str, ...]] = {}
    for gaff_name, pdb_name in ext_map.items():
        # If extension is CM, H's are HM1/HM2/HM3
        # If extension is a non-C single atom (e.g., F, O), no H
        letter = pdb_name[-1] if pdb_name else ""
        ext_h_tables[gaff_name] = tuple(f"H{letter}{i+1}" for i in range(3))

    gaff_to_amber = dict(name_map)
    h_idx: Dict[str, int] = {}
    for atom in atoms:
        if atom.atomic_number != 1:
            continue
        parents = [
            (b.atom1 if b.atom2 is atom else b.atom2) for b in atom.bonds
        ]
        if not parents:
            continue
        pname_gaff = parents[0].name
        # Skip H's attached to caps (ACE CH3, NME NH+CH3)
        if pname_gaff in cap_atoms:
            continue
        # H's on extension atoms (e.g., HM1-3 on CM): use ext_h_tables
        if pname_gaff in ext_h_tables:
            pname_amber_key = ext_map[pname_gaff]
            h_names = ext_h_tables[pname_gaff]
        else:
            pname_amber_key = gaff_to_amber.get(pname_gaff)
            if pname_amber_key is None:
                continue
            h_names = h_table.get(pname_amber_key, ())
        if not h_names:
            continue
        idx = h_idx.get(pname_amber_key, 0)
        if idx < len(h_names):
            name_map[atom.name] = h_names[idx]
            h_idx[pname_amber_key] = idx + 1

    return name_map


__all__ = [
    "build_pdb_name_map_via_parent",
    "parent_heavy_bonds",
    "parent_neighbor_hint",
    "parent_h_name_table",
]
