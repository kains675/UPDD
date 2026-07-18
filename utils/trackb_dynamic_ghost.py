"""Dynamic water-only excluded-volume ghost for Track B two-copy systems.

This module implements only the force construction and declaration checks
frozen by ``analysis/dynamic_ghost_excluded_volume_20260716``. It does not run
MD, estimate a free energy, or authorize a production launch.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable


GHOST_FORCE_NAME = "UPDDDynamicGhostWCA"
ATM_NESTED_GHOST_FORCE_NAME = "UPDDATMNestedDynamicGhostWCA"
GHOST_GLOBAL_PARAMETER = "g"
GHOST_ALPHA_PARAMETER = "alpha_sc"
GHOST_FORCE_GROUP = 31

RING_ATOM_NAMES = (
    "CG",
    "CD1",
    "CD2",
    "NE1",
    "CE2",
    "CE3",
    "CZ2",
    "CZ3",
    "CH2",
)
WATER_RESIDUE_NAMES = frozenset({"HOH", "WAT", "SOL"})
WATER_OXYGEN_NAMES = frozenset({"O", "OW", "OH2"})

EXPECTED_MUTATION_SPEC = "w4a_trp_ala_res4"
EXPECTED_BINDER_CHAIN = "B"
EXPECTED_RESIDUE_ID = "4"
EXPECTED_RESIDUE_NAME = "TRP"

ALPHA_SC = 0.5
CUTOFF_NM = 0.40
SOURCE_PARAMETER_TOLERANCE = 1.0e-12

ENERGY_EXPRESSION = (
    "g*step(rc-rsc)*(4*epsilon*((sigma/rsc)^12-(sigma/rsc)^6)+epsilon);"
    "rsc=(r^6+alpha_sc*(1-g)*sigma^6)^(1.0/6.0);"
    "rc=2^(1.0/6.0)*sigma;"
    "sigma=0.5*(sigma1+sigma2);"
    "epsilon=sqrt(epsilon1*epsilon2)"
)


class DynamicGhostError(RuntimeError):
    """A frozen dynamic-ghost declaration or integrity gate failed."""


@dataclass(frozen=True)
class GhostSelection:
    ring_indices: tuple[int, ...]
    water_oxygen_indices: tuple[int, ...]
    ring_atoms: tuple[dict[str, Any], ...]
    n_water_residues: int
    water_name_pairs: tuple[tuple[str, str], ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "ring_indices": list(self.ring_indices),
            "water_oxygen_indices_sha_input": list(self.water_oxygen_indices),
            "n_ring_atoms": len(self.ring_indices),
            "n_water_oxygens": len(self.water_oxygen_indices),
            "n_water_residues": self.n_water_residues,
            "ring_atoms": list(self.ring_atoms),
            "water_name_pairs": [list(pair) for pair in self.water_name_pairs],
        }


def _as_float(value: Any, target_unit: Any) -> float:
    if hasattr(value, "value_in_unit"):
        return float(value.value_in_unit(target_unit))
    return float(value)


def _vector_nm(value: Any) -> tuple[float, float, float]:
    from openmm import unit

    vector = value.value_in_unit(unit.nanometer) if hasattr(value, "value_in_unit") else value
    return tuple(float(vector[index]) for index in range(3))


def _assert_close(observed: float, expected: float, label: str) -> None:
    if not math.isclose(
        float(observed),
        float(expected),
        rel_tol=0.0,
        abs_tol=SOURCE_PARAMETER_TOLERANCE,
    ):
        raise DynamicGhostError(
            f"{label}: observed {observed:.16g}, expected {expected:.16g} "
            f"within {SOURCE_PARAMETER_TOLERANCE:g}"
        )


def find_atm_force(system: Any) -> tuple[int, Any]:
    import openmm as mm

    matches = [
        (index, system.getForce(index))
        for index in range(system.getNumForces())
        if isinstance(system.getForce(index), mm.ATMForce)
    ]
    if len(matches) != 1:
        raise DynamicGhostError(
            f"expected exactly one top-level ATMForce, found {len(matches)}"
        )
    return matches[0]


def resolve_nested_nonbonded(atm_force: Any) -> tuple[int, Any]:
    """Return a typed copy of the canonical NonbondedForce inside ATMForce."""
    import openmm as mm

    matches: list[tuple[int, Any]] = []
    for index in range(atm_force.getNumForces()):
        nested = atm_force.getForce(index)
        if isinstance(nested, mm.NonbondedForce):
            typed = nested
        else:
            typed = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(nested))
        if isinstance(typed, mm.NonbondedForce):
            matches.append((index, typed))
    if len(matches) != 1:
        raise DynamicGhostError(
            "expected exactly one NonbondedForce nested inside ATMForce, "
            f"found {len(matches)}"
        )
    return matches[0]


def _typed_force(force: Any) -> Any:
    """Return the concrete OpenMM force type for an ATM nested-force proxy."""
    import openmm as mm

    return mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(force))


def select_ghost_atoms(
    topology: Any,
    *,
    binder_chain: str = EXPECTED_BINDER_CHAIN,
    residue_id: str = EXPECTED_RESIDUE_ID,
) -> GhostSelection:
    """Select the frozen disappearing TRP ring and all explicit-water oxygens."""
    ring_name_set = set(RING_ATOM_NAMES)
    candidate_residues: list[Any] = []
    for residue in topology.residues():
        if (
            residue.chain.id == binder_chain
            and str(residue.id) == str(residue_id)
            and residue.name == EXPECTED_RESIDUE_NAME
        ):
            candidate_residues.append(residue)
    if len(candidate_residues) != 1:
        raise DynamicGhostError(
            "expected exactly one disappearing TRP at "
            f"chain {binder_chain} residue {residue_id}, found "
            f"{len(candidate_residues)}"
        )

    residue = candidate_residues[0]
    by_name: dict[str, Any] = {}
    for atom in residue.atoms():
        if atom.name in ring_name_set:
            if atom.name in by_name:
                raise DynamicGhostError(
                    f"duplicate ring atom name {atom.name} in target TRP"
                )
            by_name[atom.name] = atom
    if set(by_name) != ring_name_set:
        missing = sorted(ring_name_set - set(by_name))
        extra = sorted(set(by_name) - ring_name_set)
        raise DynamicGhostError(
            f"target TRP ring selection mismatch: missing={missing}, extra={extra}"
        )

    ring_atoms = tuple(by_name[name] for name in RING_ATOM_NAMES)
    ring_records = tuple(
        {
            "index": int(atom.index),
            "name": atom.name,
            "residue_name": atom.residue.name,
            "residue_id": str(atom.residue.id),
            "chain_id": atom.residue.chain.id,
            "element": atom.element.symbol if atom.element is not None else None,
        }
        for atom in ring_atoms
    )
    if any(record["element"] == "H" for record in ring_records):
        raise DynamicGhostError("ring selection contains a hydrogen")

    water_indices: list[int] = []
    water_residues = 0
    water_name_pairs: set[tuple[str, str]] = set()
    for water in topology.residues():
        if water.name not in WATER_RESIDUE_NAMES:
            continue
        water_residues += 1
        oxygens = [
            atom
            for atom in water.atoms()
            if atom.name in WATER_OXYGEN_NAMES
        ]
        if len(oxygens) != 1:
            raise DynamicGhostError(
                f"water residue {water.name}{water.id} has {len(oxygens)} "
                "declared oxygen atoms; expected exactly one"
            )
        oxygen = oxygens[0]
        if oxygen.element is not None and oxygen.element.symbol != "O":
            raise DynamicGhostError(
                f"water selection {water.name}{water.id}:{oxygen.name} is not oxygen"
            )
        water_indices.append(int(oxygen.index))
        water_name_pairs.add((water.name, oxygen.name))

    if not water_indices or len(water_indices) != water_residues:
        raise DynamicGhostError(
            f"water oxygen selection mismatch: residues={water_residues}, "
            f"oxygens={len(water_indices)}"
        )
    if set(int(atom.index) for atom in ring_atoms) & set(water_indices):
        raise DynamicGhostError("ring and water interaction groups overlap")

    return GhostSelection(
        ring_indices=tuple(int(atom.index) for atom in ring_atoms),
        water_oxygen_indices=tuple(water_indices),
        ring_atoms=ring_records,
        n_water_residues=water_residues,
        water_name_pairs=tuple(sorted(water_name_pairs)),
    )


def _particle_parameters(nonbonded: Any, index: int) -> tuple[float, float, float]:
    from openmm import unit

    charge, sigma, epsilon = nonbonded.getParticleParameters(int(index))
    return (
        _as_float(charge, unit.elementary_charge),
        _as_float(sigma, unit.nanometer),
        _as_float(epsilon, unit.kilojoule_per_mole),
    )


def _summarize_source_parameters(
    nonbonded: Any,
    selection: GhostSelection,
    reference_pair_parameters: dict[str, Any],
) -> dict[str, Any]:
    water_parameters = [
        _particle_parameters(nonbonded, index)
        for index in selection.water_oxygen_indices
    ]
    water_charge, water_sigma, water_epsilon = water_parameters[0]
    for position, values in enumerate(water_parameters[1:], 1):
        _assert_close(values[0], water_charge, f"water oxygen {position} charge")
        _assert_close(values[1], water_sigma, f"water oxygen {position} sigma")
        _assert_close(values[2], water_epsilon, f"water oxygen {position} epsilon")

    ring_rows: list[dict[str, Any]] = []
    class_rows: dict[str, list[dict[str, float]]] = {
        "aromatic_C_water_O": [],
        "NE1_water_O": [],
    }
    for atom_record, index in zip(selection.ring_atoms, selection.ring_indices):
        charge, sigma, epsilon = _particle_parameters(nonbonded, index)
        mixed_sigma = 0.5 * (sigma + water_sigma)
        mixed_epsilon = math.sqrt(epsilon * water_epsilon)
        cutoff = (2.0 ** (1.0 / 6.0)) * mixed_sigma
        pair_class = "NE1_water_O" if atom_record["name"] == "NE1" else "aromatic_C_water_O"
        row = {
            **atom_record,
            "charge_e": charge,
            "sigma_nm": sigma,
            "epsilon_kj_mol": epsilon,
            "pair_class": pair_class,
            "mixed_sigma_nm": mixed_sigma,
            "mixed_epsilon_kj_mol": mixed_epsilon,
            "wca_cutoff_nm": cutoff,
        }
        ring_rows.append(row)
        class_rows[pair_class].append(
            {
                "sigma_nm": mixed_sigma,
                "epsilon_kj_mol": mixed_epsilon,
                "wca_cutoff_nm": cutoff,
            }
        )

    if len(class_rows["aromatic_C_water_O"]) != 8 or len(class_rows["NE1_water_O"]) != 1:
        raise DynamicGhostError(
            "ring pair classes must contain eight aromatic carbons and one NE1"
        )

    pair_classes: dict[str, dict[str, float]] = {}
    for name, rows in class_rows.items():
        first = rows[0]
        for row in rows[1:]:
            for key, value in first.items():
                _assert_close(row[key], value, f"{name} {key} homogeneity")
        expected = reference_pair_parameters.get(name)
        if not isinstance(expected, dict):
            raise DynamicGhostError(f"missing frozen reference pair class {name}")
        for key in ("sigma_nm", "epsilon_kj_mol", "wca_cutoff_nm"):
            _assert_close(first[key], float(expected[key]), f"{name} {key}")
        pair_classes[name] = dict(first)

    return {
        "water_oxygen": {
            "charge_e": water_charge,
            "sigma_nm": water_sigma,
            "epsilon_kj_mol": water_epsilon,
            "n_atoms": len(water_parameters),
        },
        "ring_atoms": ring_rows,
        "pair_classes": pair_classes,
        "source_parameter_tolerance": SOURCE_PARAMETER_TOLERANCE,
    }


def validate_ghost_protocol(protocol: dict[str, Any]) -> dict[str, Any]:
    ghost = protocol.get("ghost") or {}
    expected = {
        "implementation": "top_level_CustomNonbondedForce_interaction_group",
        "coordinate_policy": "stored_coordinates_of_disappearing_trp_atoms",
        "ring_atom_names": list(RING_ATOM_NAMES),
        "water_residue_names": sorted(WATER_RESIDUE_NAMES),
        "water_oxygen_names": sorted(WATER_OXYGEN_NAMES),
        "electrostatics": False,
        "attraction": False,
        "water_deletion": False,
        "added_particles": 0,
        "potential": "softcore_WCA",
        "alpha_sc": ALPHA_SC,
        "mixing_rule": "Lorentz-Berthelot",
        "nonbonded_method": "CutoffPeriodic",
        "cutoff_nm": CUTOFF_NM,
        "long_range_correction": False,
        "force_group_required": True,
        "energy_parameter_derivative_required": True,
    }
    for key, value in expected.items():
        observed = ghost.get(key)
        if key in {"water_residue_names", "water_oxygen_names"}:
            observed = sorted(observed or [])
        if observed != value:
            raise DynamicGhostError(
                f"frozen ghost protocol mismatch for {key}: {observed!r} != {value!r}"
            )
    system = protocol.get("system") or {}
    if system.get("mutation_spec") != EXPECTED_MUTATION_SPEC:
        raise DynamicGhostError("dynamic ghost is frozen only for w4a_trp_ala_res4")
    return ghost


def validate_atm_nested_ghost_protocol(
    protocol: dict[str, Any],
) -> dict[str, Any]:
    """Validate the frozen transformed-coordinate ATM-nested ghost contract."""
    ghost = protocol.get("ghost") or {}
    expected = {
        "implementation": "ATMForce_nested_CustomNonbondedForce_interaction_group",
        "coordinate_policy": "u0_stored_u1_ATM_transformed_coordinates",
        "force_name": ATM_NESTED_GHOST_FORCE_NAME,
        "ring_atom_names": list(RING_ATOM_NAMES),
        "water_residue_names": sorted(WATER_RESIDUE_NAMES),
        "water_oxygen_names": sorted(WATER_OXYGEN_NAMES),
        "electrostatics": False,
        "attraction": False,
        "water_deletion": False,
        "added_particles": 0,
        "potential": "softcore_WCA",
        "alpha_sc": ALPHA_SC,
        "mixing_rule": "Lorentz-Berthelot",
        "nonbonded_method": "CutoffPeriodic",
        "cutoff_nm": CUTOFF_NM,
        "long_range_correction": False,
        "switching_function": False,
        "top_level_force_count_change": 0,
        "atm_nested_force_count_change": 1,
        "energy_parameter_derivative_required": True,
    }
    for key, value in expected.items():
        observed = ghost.get(key)
        if key in {"water_residue_names", "water_oxygen_names"}:
            observed = sorted(observed or [])
        if observed != value:
            raise DynamicGhostError(
                f"frozen ATM-nested ghost protocol mismatch for {key}: "
                f"{observed!r} != {value!r}"
            )
    system = protocol.get("system") or {}
    if system.get("mutation_spec") != EXPECTED_MUTATION_SPEC:
        raise DynamicGhostError(
            "ATM-nested dynamic ghost is frozen only for w4a_trp_ala_res4"
        )
    return ghost


def build_dynamic_ghost_force(
    system: Any,
    topology: Any,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    """Append the frozen ghost as a separate top-level CustomNonbondedForce."""
    import openmm as mm
    from openmm import unit

    ghost_config = validate_ghost_protocol(protocol)
    n_particles_before = int(system.getNumParticles())
    n_forces_before = int(system.getNumForces())
    atm_index, atm_force = find_atm_force(system)
    nested_index, nonbonded = resolve_nested_nonbonded(atm_force)
    if nonbonded.getNumParticles() != n_particles_before:
        raise DynamicGhostError(
            "nested NonbondedForce particle count does not match System"
        )
    if topology.getNumAtoms() != n_particles_before:
        raise DynamicGhostError(
            f"topology/System count mismatch: {topology.getNumAtoms()} != "
            f"{n_particles_before}"
        )
    for index in range(n_forces_before):
        if system.getForce(index).getName() == GHOST_FORCE_NAME:
            raise DynamicGhostError("system already contains a dynamic ghost force")
    if any(system.getForce(index).getForceGroup() == GHOST_FORCE_GROUP for index in range(n_forces_before)):
        raise DynamicGhostError(
            f"force group {GHOST_FORCE_GROUP} is already used by a top-level force"
        )

    selection = select_ghost_atoms(topology)
    source_parameters = _summarize_source_parameters(
        nonbonded,
        selection,
        ghost_config.get("reference_pair_parameters") or {},
    )

    force = mm.CustomNonbondedForce(ENERGY_EXPRESSION)
    force.setName(GHOST_FORCE_NAME)
    force.addGlobalParameter(GHOST_GLOBAL_PARAMETER, 0.0)
    force.addGlobalParameter(GHOST_ALPHA_PARAMETER, ALPHA_SC)
    force.addPerParticleParameter("sigma")
    force.addPerParticleParameter("epsilon")

    net_charge = 0.0
    for index in range(n_particles_before):
        charge, sigma, epsilon = _particle_parameters(nonbonded, index)
        net_charge += charge
        force.addParticle([sigma * unit.nanometer, epsilon * unit.kilojoule_per_mole])

    force.addInteractionGroup(
        set(selection.ring_indices),
        set(selection.water_oxygen_indices),
    )
    force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    force.setCutoffDistance(CUTOFF_NM * unit.nanometer)
    force.setUseLongRangeCorrection(False)
    force.setUseSwitchingFunction(False)
    force.addEnergyParameterDerivative(GHOST_GLOBAL_PARAMETER)
    force.setForceGroup(GHOST_FORCE_GROUP)
    ghost_index = int(system.addForce(force))

    if ghost_index <= atm_index or ghost_index != n_forces_before:
        raise DynamicGhostError(
            "dynamic ghost was not appended after the canonical ATMForce"
        )
    if system.getNumParticles() != n_particles_before:
        raise DynamicGhostError("dynamic ghost changed the System particle count")

    report = inspect_dynamic_ghost_force(system)
    report.update(
        {
            "atm_force_index": atm_index,
            "nested_nonbonded_index": nested_index,
            "n_particles_before": n_particles_before,
            "n_particles_after": int(system.getNumParticles()),
            "n_forces_before": n_forces_before,
            "n_forces_after": int(system.getNumForces()),
            "net_charge_before_e": net_charge,
            "net_charge_after_e": net_charge,
            "net_charge_delta_e": 0.0,
            "selection": selection.as_dict(),
            "source_parameters": source_parameters,
        }
    )
    return report


def inspect_dynamic_ghost_force(system: Any) -> dict[str, Any]:
    import openmm as mm
    from openmm import unit

    matches = [
        (index, system.getForce(index))
        for index in range(system.getNumForces())
        if isinstance(system.getForce(index), mm.CustomNonbondedForce)
        and system.getForce(index).getName() == GHOST_FORCE_NAME
    ]
    if len(matches) != 1:
        raise DynamicGhostError(
            f"expected exactly one {GHOST_FORCE_NAME}, found {len(matches)}"
        )
    index, force = matches[0]
    global_parameters = {
        force.getGlobalParameterName(i): float(force.getGlobalParameterDefaultValue(i))
        for i in range(force.getNumGlobalParameters())
    }
    per_particle_parameters = [
        force.getPerParticleParameterName(i)
        for i in range(force.getNumPerParticleParameters())
    ]
    derivatives = [
        force.getEnergyParameterDerivativeName(i)
        for i in range(force.getNumEnergyParameterDerivatives())
    ]
    interaction_groups: list[dict[str, Any]] = []
    for i in range(force.getNumInteractionGroups()):
        first, second = force.getInteractionGroupParameters(i)
        interaction_groups.append(
            {
                "first": sorted(int(value) for value in first),
                "second": sorted(int(value) for value in second),
            }
        )

    expected_globals = {
        GHOST_GLOBAL_PARAMETER: 0.0,
        GHOST_ALPHA_PARAMETER: ALPHA_SC,
    }
    if global_parameters != expected_globals:
        raise DynamicGhostError(
            f"ghost global parameters drifted: {global_parameters!r}"
        )
    if per_particle_parameters != ["sigma", "epsilon"]:
        raise DynamicGhostError(
            f"ghost per-particle parameters drifted: {per_particle_parameters!r}"
        )
    if derivatives != [GHOST_GLOBAL_PARAMETER]:
        raise DynamicGhostError(
            f"ghost derivative declaration drifted: {derivatives!r}"
        )
    if force.getEnergyFunction() != ENERGY_EXPRESSION:
        raise DynamicGhostError("ghost energy expression drifted")
    if force.getNonbondedMethod() != mm.CustomNonbondedForce.CutoffPeriodic:
        raise DynamicGhostError("ghost nonbonded method is not CutoffPeriodic")
    cutoff_nm = _as_float(force.getCutoffDistance(), unit.nanometer)
    _assert_close(cutoff_nm, CUTOFF_NM, "ghost cutoff")
    if force.getUseLongRangeCorrection() or force.getUseSwitchingFunction():
        raise DynamicGhostError("ghost LRC/switching must be disabled")
    if force.getForceGroup() != GHOST_FORCE_GROUP:
        raise DynamicGhostError("ghost force group drifted")
    if len(interaction_groups) != 1:
        raise DynamicGhostError("ghost must contain exactly one interaction group")

    return {
        "force_index": int(index),
        "force_name": force.getName(),
        "force_group": int(force.getForceGroup()),
        "energy_expression": force.getEnergyFunction(),
        "global_parameters": global_parameters,
        "per_particle_parameters": per_particle_parameters,
        "energy_parameter_derivatives": derivatives,
        "n_particles": int(force.getNumParticles()),
        "nonbonded_method": "CutoffPeriodic",
        "cutoff_nm": cutoff_nm,
        "long_range_correction": bool(force.getUseLongRangeCorrection()),
        "switching_function": bool(force.getUseSwitchingFunction()),
        "interaction_groups": interaction_groups,
    }


def inspect_atm_coordinate_transformations(
    system: Any,
    topology: Any,
    *,
    positions_nm: Any | None = None,
    protocol: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate the ring and water endpoint transformations used by the ghost."""
    import openmm as mm

    selection = select_ghost_atoms(topology)
    atm_index, atm_force = find_atm_force(system)

    ring_rows: list[dict[str, Any]] = []
    u1_pairs: set[tuple[int, int]] = set()
    for atom, index in zip(selection.ring_atoms, selection.ring_indices):
        transformation = atm_force.getParticleTransformation(int(index))
        if not isinstance(transformation, mm.ParticleOffsetDisplacement):
            raise DynamicGhostError(
                f"ring atom {index} uses {type(transformation).__name__}; "
                "expected ParticleOffsetDisplacement"
            )
        destination1 = int(transformation.getDestinationParticle1())
        origin1 = int(transformation.getOriginParticle1())
        destination0 = int(transformation.getDestinationParticle0())
        origin0 = int(transformation.getOriginParticle0())
        if destination1 < 0 or origin1 < 0:
            raise DynamicGhostError(f"ring atom {index} has an invalid u1 offset")
        if destination0 != -1 or origin0 != -1:
            raise DynamicGhostError(f"ring atom {index} unexpectedly displaces u0")
        u1_pairs.add((destination1, origin1))
        ring_rows.append(
            {
                "index": int(index),
                "name": atom["name"],
                "type": "ParticleOffsetDisplacement",
                "destination1": destination1,
                "origin1": origin1,
                "destination0": destination0,
                "origin0": origin0,
            }
        )
    if len(u1_pairs) != 1:
        raise DynamicGhostError(
            f"ring atoms do not share one u1 particle offset: {sorted(u1_pairs)}"
        )

    zero = (0.0, 0.0, 0.0)
    for index in selection.water_oxygen_indices:
        transformation = atm_force.getParticleTransformation(int(index))
        if not isinstance(transformation, mm.FixedDisplacement):
            raise DynamicGhostError(
                f"water oxygen {index} uses {type(transformation).__name__}; "
                "expected FixedDisplacement"
            )
        displacement0 = _vector_nm(transformation.getFixedDisplacement0())
        displacement1 = _vector_nm(transformation.getFixedDisplacement1())
        if displacement0 != zero or displacement1 != zero:
            raise DynamicGhostError(
                f"water oxygen {index} has a nonzero endpoint displacement"
            )

    destination1, origin1 = next(iter(u1_pairs))
    offset_nm = None
    offset_norm_nm = None
    if positions_nm is not None:
        destination = positions_nm[destination1]
        origin = positions_nm[origin1]
        offset_nm = tuple(float(destination[i] - origin[i]) for i in range(3))
        offset_norm_nm = math.sqrt(sum(value * value for value in offset_nm))

    if protocol is not None:
        transformation = protocol.get("transformation") or {}
        expected = {
            "ring_type": "ParticleOffsetDisplacement",
            "ring_common_u1_offset_required": True,
            "ring_u0_displacement": False,
            "expected_u1_offset_norm_nm": 4.0,
            "water_type": "FixedDisplacement",
            "water_zero_u0_u1_required": True,
        }
        if any(transformation.get(key) != value for key, value in expected.items()):
            raise DynamicGhostError(
                f"frozen ATM transformation protocol drifted: {transformation!r}"
            )
        offset_tolerance = transformation.get("offset_norm_tolerance_nm")
        if (
            not isinstance(offset_tolerance, (int, float))
            or not math.isfinite(float(offset_tolerance))
            or not 0.0 < float(offset_tolerance) <= 1.0e-3
        ):
            raise DynamicGhostError(
                "ATM transformation offset tolerance must be in (0, 1e-3] nm"
            )
        if set(transformation) != {*expected, "offset_norm_tolerance_nm"}:
            raise DynamicGhostError(
                f"frozen ATM transformation fields drifted: {transformation!r}"
            )
        if offset_norm_nm is not None and not math.isclose(
            offset_norm_nm,
            float(transformation["expected_u1_offset_norm_nm"]),
            rel_tol=0.0,
            abs_tol=float(offset_tolerance),
        ):
            raise DynamicGhostError(
                f"ring u1 offset norm {offset_norm_nm:.16g} nm differs from "
                f"{transformation['expected_u1_offset_norm_nm']:.16g} nm"
            )

    return {
        "atm_force_index": int(atm_index),
        "ring": ring_rows,
        "common_u1_offset": {
            "destination1": destination1,
            "origin1": origin1,
            "vector_nm": list(offset_nm) if offset_nm is not None else None,
            "norm_nm": offset_norm_nm,
        },
        "n_ring_atoms": len(selection.ring_indices),
        "n_water_oxygens": len(selection.water_oxygen_indices),
        "water_oxygen_transform": "zero_FixedDisplacement",
    }


def build_atm_nested_dynamic_ghost_force(
    system: Any,
    topology: Any,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    """Append the frozen WCA ghost inside the canonical ATMForce."""
    import openmm as mm
    from openmm import unit

    ghost_config = validate_atm_nested_ghost_protocol(protocol)
    n_particles_before = int(system.getNumParticles())
    n_top_level_before = int(system.getNumForces())
    atm_index, atm_force = find_atm_force(system)
    n_nested_before = int(atm_force.getNumForces())
    canonical_nonbonded_index, nonbonded = resolve_nested_nonbonded(atm_force)
    if nonbonded.getNumParticles() != n_particles_before:
        raise DynamicGhostError(
            "nested NonbondedForce particle count does not match System"
        )
    if topology.getNumAtoms() != n_particles_before:
        raise DynamicGhostError(
            f"topology/System count mismatch: {topology.getNumAtoms()} != "
            f"{n_particles_before}"
        )

    forbidden_names = {GHOST_FORCE_NAME, ATM_NESTED_GHOST_FORCE_NAME}
    for index in range(n_top_level_before):
        if system.getForce(index).getName() in forbidden_names:
            raise DynamicGhostError(
                "system already contains a top-level dynamic ghost force"
            )
    for index in range(n_nested_before):
        if atm_force.getForce(index).getName() in forbidden_names:
            raise DynamicGhostError(
                "ATMForce already contains a dynamic ghost force"
            )

    selection = select_ghost_atoms(topology)
    transformation_report = inspect_atm_coordinate_transformations(
        system,
        topology,
        protocol=protocol,
    )
    source_parameters = _summarize_source_parameters(
        nonbonded,
        selection,
        ghost_config.get("reference_pair_parameters") or {},
    )

    force = mm.CustomNonbondedForce(ENERGY_EXPRESSION)
    force.setName(ATM_NESTED_GHOST_FORCE_NAME)
    force.addGlobalParameter(GHOST_GLOBAL_PARAMETER, 0.0)
    force.addGlobalParameter(GHOST_ALPHA_PARAMETER, ALPHA_SC)
    force.addPerParticleParameter("sigma")
    force.addPerParticleParameter("epsilon")

    net_charge = 0.0
    for index in range(n_particles_before):
        charge, sigma, epsilon = _particle_parameters(nonbonded, index)
        net_charge += charge
        force.addParticle(
            [sigma * unit.nanometer, epsilon * unit.kilojoule_per_mole]
        )
    force.addInteractionGroup(
        set(selection.ring_indices),
        set(selection.water_oxygen_indices),
    )
    force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    force.setCutoffDistance(CUTOFF_NM * unit.nanometer)
    force.setUseLongRangeCorrection(False)
    force.setUseSwitchingFunction(False)
    force.addEnergyParameterDerivative(GHOST_GLOBAL_PARAMETER)
    force.setForceGroup(0)

    nested_ghost_index = int(atm_force.addForce(force))
    if nested_ghost_index != n_nested_before:
        raise DynamicGhostError(
            "ATM-nested dynamic ghost was not appended after existing forces"
        )
    if system.getNumForces() != n_top_level_before:
        raise DynamicGhostError("ATM-nested ghost changed top-level force count")
    if system.getNumParticles() != n_particles_before:
        raise DynamicGhostError("ATM-nested ghost changed System particle count")

    report = inspect_atm_nested_dynamic_ghost_force(system)
    report.update(
        {
            "canonical_nonbonded_index": canonical_nonbonded_index,
            "n_particles_before": n_particles_before,
            "n_particles_after": int(system.getNumParticles()),
            "n_top_level_forces_before": n_top_level_before,
            "n_top_level_forces_after": int(system.getNumForces()),
            "n_atm_nested_forces_before": n_nested_before,
            "n_atm_nested_forces_after": int(atm_force.getNumForces()),
            "net_charge_before_e": net_charge,
            "net_charge_after_e": net_charge,
            "net_charge_delta_e": 0.0,
            "selection": selection.as_dict(),
            "source_parameters": source_parameters,
            "transformations": transformation_report,
        }
    )
    return report


def inspect_atm_nested_dynamic_ghost_force(system: Any) -> dict[str, Any]:
    """Inspect and validate the concrete ghost force nested in ATMForce."""
    import openmm as mm
    from openmm import unit

    atm_index, atm_force = find_atm_force(system)
    named: list[tuple[int, Any]] = []
    for index in range(atm_force.getNumForces()):
        nested = atm_force.getForce(index)
        if nested.getName() == ATM_NESTED_GHOST_FORCE_NAME:
            named.append((index, _typed_force(nested)))
    if len(named) != 1:
        raise DynamicGhostError(
            f"expected exactly one nested {ATM_NESTED_GHOST_FORCE_NAME}, "
            f"found {len(named)}"
        )
    nested_index, force = named[0]
    if not isinstance(force, mm.CustomNonbondedForce):
        raise DynamicGhostError(
            f"nested {ATM_NESTED_GHOST_FORCE_NAME} is "
            f"{type(force).__name__}, expected CustomNonbondedForce"
        )

    global_parameters = {
        force.getGlobalParameterName(i): float(
            force.getGlobalParameterDefaultValue(i)
        )
        for i in range(force.getNumGlobalParameters())
    }
    per_particle_parameters = [
        force.getPerParticleParameterName(i)
        for i in range(force.getNumPerParticleParameters())
    ]
    derivatives = [
        force.getEnergyParameterDerivativeName(i)
        for i in range(force.getNumEnergyParameterDerivatives())
    ]
    interaction_groups: list[dict[str, Any]] = []
    for index in range(force.getNumInteractionGroups()):
        first, second = force.getInteractionGroupParameters(index)
        interaction_groups.append(
            {
                "first": sorted(int(value) for value in first),
                "second": sorted(int(value) for value in second),
            }
        )

    expected_globals = {
        GHOST_GLOBAL_PARAMETER: 0.0,
        GHOST_ALPHA_PARAMETER: ALPHA_SC,
    }
    if global_parameters != expected_globals:
        raise DynamicGhostError(
            f"nested ghost global parameters drifted: {global_parameters!r}"
        )
    if per_particle_parameters != ["sigma", "epsilon"]:
        raise DynamicGhostError(
            "nested ghost per-particle parameters drifted: "
            f"{per_particle_parameters!r}"
        )
    if derivatives != [GHOST_GLOBAL_PARAMETER]:
        raise DynamicGhostError(
            f"nested ghost derivative declaration drifted: {derivatives!r}"
        )
    if force.getEnergyFunction() != ENERGY_EXPRESSION:
        raise DynamicGhostError("nested ghost energy expression drifted")
    if force.getNonbondedMethod() != mm.CustomNonbondedForce.CutoffPeriodic:
        raise DynamicGhostError(
            "nested ghost nonbonded method is not CutoffPeriodic"
        )
    cutoff_nm = _as_float(force.getCutoffDistance(), unit.nanometer)
    _assert_close(cutoff_nm, CUTOFF_NM, "nested ghost cutoff")
    if force.getUseLongRangeCorrection() or force.getUseSwitchingFunction():
        raise DynamicGhostError("nested ghost LRC/switching must be disabled")
    if force.getForceGroup() != 0:
        raise DynamicGhostError("nested ghost force group must remain zero")
    if force.getNumParticles() != system.getNumParticles():
        raise DynamicGhostError(
            "nested ghost particle count does not match the System"
        )
    if len(interaction_groups) != 1:
        raise DynamicGhostError(
            "nested ghost must contain exactly one interaction group"
        )

    return {
        "atm_force_index": int(atm_index),
        "nested_force_index": int(nested_index),
        "force_name": force.getName(),
        "force_group": int(force.getForceGroup()),
        "energy_expression": force.getEnergyFunction(),
        "global_parameters": global_parameters,
        "per_particle_parameters": per_particle_parameters,
        "energy_parameter_derivatives": derivatives,
        "n_particles": int(force.getNumParticles()),
        "nonbonded_method": "CutoffPeriodic",
        "cutoff_nm": cutoff_nm,
        "long_range_correction": bool(force.getUseLongRangeCorrection()),
        "switching_function": bool(force.getUseSwitchingFunction()),
        "interaction_groups": interaction_groups,
    }


def mixed_pair_parameters(
    ring_rows: Iterable[dict[str, Any]],
    water_row: dict[str, Any],
) -> list[dict[str, float]]:
    """Return mixed parameters in ring order for analytic audit calculations."""
    rows: list[dict[str, float]] = []
    for ring in ring_rows:
        sigma = 0.5 * (float(ring["sigma_nm"]) + float(water_row["sigma_nm"]))
        epsilon = math.sqrt(
            float(ring["epsilon_kj_mol"])
            * float(water_row["epsilon_kj_mol"])
        )
        rows.append({"sigma_nm": sigma, "epsilon_kj_mol": epsilon})
    return rows
