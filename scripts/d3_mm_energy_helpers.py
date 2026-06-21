#!/usr/bin/env python
"""D3 MM-energy helpers (OpenMM, qmmm env).

Supports the MM-subtraction step of the D3 dihedral refit: the QM relaxed-scan
energy E_QM(phi) contains BOTH the target torsion and every nonbonded /
1-4 / coupled-bonded contribution that changes as the molecule relaxes around
phi. Fitting a cosine series straight onto E_QM(phi) would let the fitted V_n
silently ABSORB those nonbonded contributions. The standard AMBER torsion-
refit protocol therefore fits

    E_target_torsion(phi) = E_QM(phi) - [ E_MM(phi) - E_MM_no_target(phi) ]

where the bracket is the MM energy attributable to the single torsion type that
is being refit, evaluated at the scan-point geometry. At a FIXED geometry,
removing exactly one PeriodicTorsionForce term changes the total MM energy by
exactly that term's value -- every other force (nonbonded, 1-4, bonds, angles,
other torsions) is identical between the two systems and cancels in the
difference. So the bracket isolates the MM contribution of the term we are
replacing with the QM-derived V_n.

This module provides two routes to that bracket, which cross-validate:

  (A) full-OpenMM route (`mm_torsion_contribution_openmm`): build the System
      from the supplied force field, then build a SECOND System with the
      matching PeriodicTorsionForce term(s) zeroed, and take the energy
      difference at each scan-point geometry. This is the literal
      "E_MM - E_MM_without_target_torsion" the protocol asks for.

  (B) analytic route (`mm_torsion_contribution_analytic`): evaluate the AMBER
      cosine series of the target torsion type directly from its per-term
      (periodicity, phase, k) parameters at the measured dihedral. Used as a
      fallback when force-field template matching for a capped ncAA dipeptide
      is unavailable, and as an independent check on route (A).

NOTE: this module NEVER edits the frozen MTR force-field XML. It only reads it
(via openmm.app.ForceField) and, for route (A), clones the resulting System in
memory.

Run env: conda env `qmmm` (/home/san/miniconda3/envs/qmmm/bin/python, OpenMM 8.4).
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Analytic AMBER torsion energy (route B + cross-check for route A).
# ---------------------------------------------------------------------------
def amber_torsion_energy_kcal(phi_deg: float,
                              terms: Sequence[Tuple[int, float, float]]) -> float:
    """AMBER periodic-torsion energy at dihedral `phi_deg` (degrees) for a list
    of cosine terms.

    `terms` is a list of (periodicity_n, phase_deg, k_kcal_per_mol). The AMBER
    functional form is  E = sum_n k_n * (1 + cos(n*phi - phase_n))  with k in
    kcal/mol (this is the OpenMM/AMBER convention where k already absorbs the
    1/2; i.e. OpenMM `PeriodicTorsionForce` energy = k*(1+cos(n*theta-phase))).

    Returns energy in kcal/mol.
    """
    phi = math.radians(phi_deg)
    e = 0.0
    for n, phase_deg, k in terms:
        e += float(k) * (1.0 + math.cos(n * phi - math.radians(phase_deg)))
    return e


def mm_torsion_contribution_analytic(phi_deg_list: Sequence[float],
                                     terms: Sequence[Tuple[int, float, float]]
                                     ) -> List[float]:
    """Route B: MM torsion contribution (kcal/mol) at each scan-point dihedral,
    computed analytically from the target torsion type's cosine terms.

    `terms` come from the force-field XML PeriodicTorsionForce entry for the
    torsion being refit (see `extract_torsion_terms_from_xml`)."""
    return [amber_torsion_energy_kcal(p, terms) for p in phi_deg_list]


# ---------------------------------------------------------------------------
# Read the target torsion's cosine terms straight from the FF XML (read-only).
# ---------------------------------------------------------------------------
def extract_torsion_terms_from_xml(
    xml_path: str,
    class_quad: Tuple[str, str, str, str],
) -> List[Tuple[int, float, float]]:
    """Read the PeriodicTorsionForce entry matching the atom-CLASS quadruple
    `class_quad` (e.g. ('ca','na','c3','hc') for CD1-NE1-CM-HM) from the FF XML
    and return its cosine terms as (periodicity, phase_deg, k_kcal_per_mol).

    AMBER/OpenMM XML stores k in kJ/mol and phase in radians; we convert to
    kcal/mol and degrees so the analytic energy matches the kcal-mol QM scan.
    Matching is order-and-reverse-insensitive on the class quadruple, which is
    how OpenMM's Proper torsion matching itself behaves. Returns [] if no entry
    matches (caller decides whether that is acceptable)."""
    import xml.etree.ElementTree as ET

    KJ_PER_KCAL = 4.184
    want = tuple(class_quad)
    want_rev = tuple(reversed(class_quad))
    tree = ET.parse(xml_path)
    root = tree.getroot()
    out: List[Tuple[int, float, float]] = []
    for ptf in root.findall(".//PeriodicTorsionForce"):
        for proper in ptf.findall("Proper"):
            quad = (
                proper.get("class1") or proper.get("type1"),
                proper.get("class2") or proper.get("type2"),
                proper.get("class3") or proper.get("type3"),
                proper.get("class4") or proper.get("type4"),
            )
            if quad != want and quad != want_rev:
                continue
            # collect periodicity{i}/phase{i}/k{i} term sets
            i = 1
            while proper.get(f"periodicity{i}") is not None:
                n = int(proper.get(f"periodicity{i}"))
                phase_rad = float(proper.get(f"phase{i}"))
                k_kj = float(proper.get(f"k{i}"))
                out.append((n, math.degrees(phase_rad), k_kj / KJ_PER_KCAL))
                i += 1
            if out:
                return out
    return out


# ---------------------------------------------------------------------------
# Full-OpenMM route (route A): E_MM - E_MM_without_target_torsion.
# ---------------------------------------------------------------------------
def _zero_matching_torsions(system, atom_quad_0based: Tuple[int, int, int, int]) -> int:
    """In-place: zero every PeriodicTorsionForce term in `system` that acts on
    the (unordered) atom quadruple `atom_quad_0based`. Returns the number of
    terms zeroed. Zeroing k (not deleting the term) keeps every force-group /
    index layout identical so the two systems differ ONLY by this torsion."""
    import openmm

    want = frozenset(atom_quad_0based)
    n_zeroed = 0
    for force in system.getForces():
        if not isinstance(force, openmm.PeriodicTorsionForce):
            continue
        for ti in range(force.getNumTorsions()):
            p1, p2, p3, p4, periodicity, phase, k = force.getTorsionParameters(ti)
            if frozenset((p1, p2, p3, p4)) == want:
                force.setTorsionParameters(ti, p1, p2, p3, p4, periodicity, phase, 0.0)
                n_zeroed += 1
    return n_zeroed


def _single_point_mm_energy_kcal(system, topology, positions_nm) -> float:
    """Evaluate the potential energy (kcal/mol) of `system` at `positions_nm`
    (list of (x,y,z) in nanometres) on the Reference platform."""
    import openmm
    from openmm import unit

    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    platform = openmm.Platform.getPlatformByName("Reference")
    context = openmm.Context(system, integrator, platform)
    try:
        context.setPositions([openmm.Vec3(*p) * unit.nanometer for p in positions_nm])
        state = context.getState(getEnergy=True)
        e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    finally:
        del context
        del integrator
    return e_kj / 4.184


def build_normalized_topology(atoms,
                              adjacency: Dict[int, List[int]],
                              residue_rename: Optional[Dict[str, str]] = None):
    """Build an OpenMM Topology for the capped ncAA dipeptide from a D2-style
    atom list (name, residue, element, x, y, z) plus a distance-inferred
    adjacency (0-based, undirected) so EVERY bond -- intra-residue and the two
    peptide bonds -- is present. Without explicit bonds OpenMM's ForceField
    cannot template-match a non-standard residue.

    `residue_rename` maps the PDB residue name to the FF template name (e.g.
    {"OMW": "MTR"}); the FF residue template is matched by name, and the PDB
    stores the ncAA under its own name. Residues are numbered sequentially in
    one chain so the backbone is contiguous. Returns (topology, positions_nm)."""
    from openmm import app, unit, Vec3

    residue_rename = residue_rename or {}
    topology = app.Topology()
    chain = topology.addChain("A")
    elem = app.element

    omm_atoms = []
    last_res_key = None
    cur_res = None
    for (name, resname, element, _x, _y, _z) in atoms:
        rn = residue_rename.get(resname, resname)
        key = (resname, rn)
        if key != last_res_key:
            cur_res = topology.addResidue(rn, chain)
            last_res_key = key
        try:
            e = elem.get_by_symbol(element)
        except KeyError:
            e = elem.get_by_symbol(name[0])
        omm_atoms.append(topology.addAtom(name, e, cur_res))

    # all bonds from the distance adjacency (covers intra-residue + peptide).
    seen = set()
    for i, nbrs in adjacency.items():
        for j in nbrs:
            key = frozenset((i, j))
            if i != j and key not in seen:
                topology.addBond(omm_atoms[i], omm_atoms[j])
                seen.add(key)

    positions = [Vec3(a[3], a[4], a[5]) * 0.1 * unit.nanometer for a in atoms]
    return topology, positions


def build_mm_systems(forcefield_xmls: Sequence[str],
                     topology,
                     atom_quad_0based: Tuple[int, int, int, int]):
    """Build the (full, target-torsion-zeroed) OpenMM System pair from a
    pre-built normalized `topology`. The two systems differ ONLY by the target
    torsion. Returns (system_full, system_no_target, n_zeroed).

    Raises whatever openmm.app.ForceField raises on a template-matching failure
    -- the driver catches that and falls back to the analytic route."""
    import copy
    from openmm import app

    forcefield = app.ForceField(*forcefield_xmls)
    system_full = forcefield.createSystem(
        topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False,
    )
    system_no_target = copy.deepcopy(system_full)
    n_zeroed = _zero_matching_torsions(system_no_target, atom_quad_0based)
    return system_full, system_no_target, n_zeroed


def mm_torsion_contribution_openmm(
    forcefield_xmls: Sequence[str],
    atoms_ref,
    adjacency: Dict[int, List[int]],
    scan_coords_ang: Sequence[Sequence[Tuple[float, float, float]]],
    atom_quad_0based: Tuple[int, int, int, int],
    residue_rename: Optional[Dict[str, str]] = None,
) -> Tuple[List[Optional[float]], int]:
    """Route A: for each scan-point coordinate set, return
    E_MM(full) - E_MM(no-target) in kcal/mol -- the MM energy of the single
    target torsion at that relaxed geometry. Returns (contributions, n_zeroed).

    A normalized Topology (residue renamed to the FF template name, every bond
    present) is built ONCE from `atoms_ref` + `adjacency`; the two Systems differ
    ONLY by the zeroed target torsion. `scan_coords_ang` holds the per-point
    coordinates (Angstrom) in the same atom order as `atoms_ref`. A None entry
    means that point's coordinates were missing."""
    if not scan_coords_ang:
        return [], 0

    topology, _pos0 = build_normalized_topology(atoms_ref, adjacency, residue_rename)
    system_full, system_no_target, n_zeroed = build_mm_systems(
        forcefield_xmls, topology, atom_quad_0based,
    )
    contributions: List[Optional[float]] = []
    for coords_ang in scan_coords_ang:
        if coords_ang is None:
            contributions.append(None)
            continue
        try:
            pos_nm = [(c[0] * 0.1, c[1] * 0.1, c[2] * 0.1) for c in coords_ang]
            e_full = _single_point_mm_energy_kcal(system_full, topology, pos_nm)
            e_no = _single_point_mm_energy_kcal(system_no_target, topology, pos_nm)
            contributions.append(e_full - e_no)
        except Exception:  # noqa: BLE001 -- per-point failure recorded as None
            contributions.append(None)
    return contributions, n_zeroed


# ---------------------------------------------------------------------------
# Self-test (cheap, no QM, no full force field): the analytic energy must equal
# the OpenMM single-torsion energy for a hand-built 4-atom system.
# ---------------------------------------------------------------------------
def _selftest() -> int:
    """Cross-check route B (analytic) against a 4-atom OpenMM PeriodicTorsion
    over a dihedral sweep. Returns 0 on pass, 1 on failure. Uses only OpenMM +
    stdlib; no force-field XML and no QM."""
    import openmm
    from openmm import unit

    terms = [(1, 0.0, 0.8), (2, 180.0, 0.2), (3, 0.0, 0.15)]  # kcal/mol

    system = openmm.System()
    for _ in range(4):
        system.addParticle(12.0)
    ptf = openmm.PeriodicTorsionForce()
    for n, phase_deg, k in terms:
        ptf.addTorsion(0, 1, 2, 3, n,
                       math.radians(phase_deg) * unit.radian,
                       k * 4.184 * unit.kilojoule_per_mole)
    system.addForce(ptf)

    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    context = openmm.Context(system, integrator,
                             openmm.Platform.getPlatformByName("Reference"))

    max_err = 0.0
    for phi_deg in range(-180, 181, 15):
        phi = math.radians(phi_deg)
        # place 4 atoms so the 0-1-2-3 dihedral == phi (standard construction)
        positions = [
            openmm.Vec3(1.0, 1.0, 0.0),
            openmm.Vec3(1.0, 0.0, 0.0),
            openmm.Vec3(0.0, 0.0, 0.0),
            openmm.Vec3(0.0, math.cos(phi), math.sin(phi)),
        ]
        context.setPositions([p * unit.nanometer for p in positions])
        e_kj = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)
        e_omm = e_kj / 4.184
        e_ana = amber_torsion_energy_kcal(phi_deg, terms)
        max_err = max(max_err, abs(e_omm - e_ana))
    del context
    del integrator

    ok = max_err < 1e-6
    print(f"[selftest] analytic-vs-OpenMM torsion max |Δ| = {max_err:.3e} kcal/mol "
          f"-> {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(_selftest())
