"""Reference-platform tests for the frozen Track B dynamic ghost."""

from __future__ import annotations

import copy
import math

import numpy as np
import pytest

import trackb_dynamic_ghost as dg


def _openmm():
    return pytest.importorskip("openmm")


def _protocol() -> dict:
    return {
        "system": {"mutation_spec": "w4a_trp_ala_res4"},
        "ghost": {
            "implementation": "top_level_CustomNonbondedForce_interaction_group",
            "coordinate_policy": "stored_coordinates_of_disappearing_trp_atoms",
            "ring_atom_names": list(dg.RING_ATOM_NAMES),
            "water_residue_names": ["HOH", "WAT", "SOL"],
            "water_oxygen_names": ["O", "OW", "OH2"],
            "electrostatics": False,
            "attraction": False,
            "water_deletion": False,
            "added_particles": 0,
            "potential": "softcore_WCA",
            "alpha_sc": 0.5,
            "mixing_rule": "Lorentz-Berthelot",
            "nonbonded_method": "CutoffPeriodic",
            "cutoff_nm": 0.4,
            "long_range_correction": False,
            "force_group_required": True,
            "energy_parameter_derivative_required": True,
            "reference_pair_parameters": {
                "aromatic_C_water_O": {
                    "sigma_nm": 0.32,
                    "epsilon_kj_mol": 0.6,
                    "wca_cutoff_nm": (2.0 ** (1.0 / 6.0)) * 0.32,
                },
                "NE1_water_O": {
                    "sigma_nm": 0.31,
                    "epsilon_kj_mol": 0.8,
                    "wca_cutoff_nm": (2.0 ** (1.0 / 6.0)) * 0.31,
                },
            },
        },
    }


def _nested_protocol() -> dict:
    protocol = _protocol()
    ghost = protocol["ghost"]
    ghost.update(
        {
            "implementation": "ATMForce_nested_CustomNonbondedForce_interaction_group",
            "coordinate_policy": "u0_stored_u1_ATM_transformed_coordinates",
            "force_name": dg.ATM_NESTED_GHOST_FORCE_NAME,
            "switching_function": False,
            "top_level_force_count_change": 0,
            "atm_nested_force_count_change": 1,
        }
    )
    ghost.pop("force_group_required")
    protocol["transformation"] = {
        "ring_type": "ParticleOffsetDisplacement",
        "ring_common_u1_offset_required": True,
        "ring_u0_displacement": False,
        "expected_u1_offset_norm_nm": 4.0,
        "offset_norm_tolerance_nm": 1.0e-8,
        "water_type": "FixedDisplacement",
        "water_zero_u0_u1_required": True,
    }
    return protocol


def _fixture_system():
    mm = _openmm()
    from openmm import app, unit

    topology = app.Topology()
    chain = topology.addChain("B")
    ala = topology.addResidue("ALA", chain, "4")
    topology.addAtom("CB", app.element.carbon, ala)
    trp = topology.addResidue("TRP", chain, "4")
    for name in dg.RING_ATOM_NAMES:
        element = app.element.nitrogen if name == "NE1" else app.element.carbon
        topology.addAtom(name, element, trp)
    other = topology.addResidue("GLY", chain, "8")
    topology.addAtom("O", app.element.oxygen, other)

    water_chain = topology.addChain("W")
    for resid, resname, oxygen_name in (("1", "HOH", "O"), ("2", "WAT", "OW")):
        water = topology.addResidue(resname, water_chain, resid)
        topology.addAtom(oxygen_name, app.element.oxygen, water)
        topology.addAtom("H1", app.element.hydrogen, water)
        topology.addAtom("H2", app.element.hydrogen, water)

    system = mm.System()
    nonbonded = mm.NonbondedForce()
    nonbonded.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    positions = []
    for atom in topology.atoms():
        system.addParticle(1.0)
        if atom.residue.name == "TRP" and atom.name in dg.RING_ATOM_NAMES:
            sigma = 0.32 if atom.name == "NE1" else 0.34
            epsilon = 1.6 if atom.name == "NE1" else 0.9
        elif atom.residue.name in dg.WATER_RESIDUE_NAMES and atom.name in dg.WATER_OXYGEN_NAMES:
            sigma, epsilon = 0.30, 0.40
        else:
            sigma, epsilon = 0.25, 0.0
        nonbonded.addParticle(0.0, sigma, epsilon)
        positions.append(mm.Vec3(1.5, 1.5, 1.5))

    # Keep the synthetic canonical potential finite and identically zero so the
    # test isolates the added ghost interaction.
    for first in range(system.getNumParticles()):
        for second in range(first + 1, system.getNumParticles()):
            nonbonded.addException(first, second, 0.0, 1.0, 0.0)

    selection = dg.select_ghost_atoms(topology)
    for offset, atom_index in enumerate(selection.ring_indices):
        positions[atom_index] = mm.Vec3(0.5, 0.5 + 0.01 * offset, 0.5)
    positions[selection.water_oxygen_indices[0]] = mm.Vec3(0.78, 0.5, 0.5)
    positions[selection.water_oxygen_indices[1]] = mm.Vec3(1.5, 1.5, 1.5)

    atm = mm.ATMForce("u0")
    atm.addForce(copy.copy(nonbonded))
    for _ in range(system.getNumParticles()):
        atm.addParticle()
    system.addForce(atm)
    system.setDefaultPeriodicBoxVectors(
        mm.Vec3(2.0, 0.0, 0.0),
        mm.Vec3(0.0, 2.0, 0.0),
        mm.Vec3(0.0, 0.0, 2.0),
    )
    return system, topology, positions * unit.nanometer, selection


def _nested_fixture_system():
    mm = _openmm()
    from openmm import unit

    system, topology, positions, selection = _fixture_system()
    _atm_index, atm = dg.find_atm_force(system)
    atm.setEnergyFunction("(1-select)*u0+select*u1")
    atm.addGlobalParameter("select", 0.0)

    atoms = list(topology.atoms())
    origin = next(atom.index for atom in atoms if atom.name == "CB")
    destination = next(
        atom.index
        for atom in atoms
        if atom.residue.name == "GLY" and atom.name == "O"
    )
    position_array = np.asarray(
        positions.value_in_unit(unit.nanometer),
        dtype=float,
    )
    position_array[origin] = (0.2, 0.2, 0.2)
    position_array[destination] = (1.0, 0.2, 0.2)
    for index in selection.ring_indices:
        atm.setParticleTransformation(
            int(index),
            mm.ParticleOffsetDisplacement(destination, origin),
        )
    return (
        system,
        topology,
        position_array * unit.nanometer,
        selection,
        np.array([0.8, 0.0, 0.0]),
    )


def _state(system, positions, *, g=None):
    mm = _openmm()
    from openmm import unit

    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(
        system,
        integrator,
        mm.Platform.getPlatformByName("Reference"),
    )
    context.setPositions(positions)
    if g is not None:
        context.setParameter(dg.GHOST_GLOBAL_PARAMETER, float(g))
    state = context.getState(
        getEnergy=True,
        getForces=True,
        getParameterDerivatives=g is not None,
    )
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=float,
    )
    derivatives = dict(state.getEnergyParameterDerivatives()) if g is not None else {}
    ghost_energy = None
    if g is not None:
        ghost_energy = context.getState(
            getEnergy=True,
            groups=1 << dg.GHOST_FORCE_GROUP,
        ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    del context, integrator
    return energy, forces, derivatives, ghost_energy


def _nested_state(system, positions, *, endpoint, g):
    mm = _openmm()
    from openmm import unit

    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(
        system,
        integrator,
        mm.Platform.getPlatformByName("Reference"),
    )
    context.setPositions(positions)
    context.setParameter("select", float(endpoint))
    context.setParameter(dg.GHOST_GLOBAL_PARAMETER, float(g))
    state = context.getState(
        getEnergy=True,
        getForces=True,
        getParameterDerivatives=True,
    )
    energy = float(
        state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    )
    forces = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=float,
    )
    derivatives = {
        str(name): float(value)
        for name, value in dict(state.getEnergyParameterDerivatives()).items()
    }
    _atm_index, atm = dg.find_atm_force(system)
    u1, u0, _alpha = atm.getPerturbationEnergy(context)
    raw = {
        "u0": float(u0.value_in_unit(unit.kilojoule_per_mole)),
        "u1": float(u1.value_in_unit(unit.kilojoule_per_mole)),
    }
    del context, integrator
    return energy, forces, derivatives, raw


def test_selection_is_exact_ring_and_all_water_oxygens():
    _system, _topology, _positions, selection = _fixture_system()

    assert [row["name"] for row in selection.ring_atoms] == list(dg.RING_ATOM_NAMES)
    assert len(selection.ring_indices) == 9
    assert len(selection.water_oxygen_indices) == 2
    assert selection.n_water_residues == 2
    assert selection.water_name_pairs == (("HOH", "O"), ("WAT", "OW"))
    assert not set(selection.ring_indices) & set(selection.water_oxygen_indices)


def test_builder_preserves_particles_charge_and_top_level_scope():
    system, topology, _positions, selection = _fixture_system()
    n_particles = system.getNumParticles()
    report = dg.build_dynamic_ghost_force(system, topology, _protocol())

    assert report["n_particles_before"] == report["n_particles_after"] == n_particles
    assert report["net_charge_delta_e"] == 0.0
    assert report["force_index"] > report["atm_force_index"]
    assert report["force_group"] == 31
    assert report["energy_parameter_derivatives"] == ["g"]
    assert report["long_range_correction"] is False
    assert report["switching_function"] is False
    group = report["interaction_groups"][0]
    assert group["first"] == sorted(selection.ring_indices)
    assert group["second"] == sorted(selection.water_oxygen_indices)


def test_g_zero_is_exact_energy_force_parity_and_g_one_is_repulsive():
    mm = _openmm()
    system, topology, positions, _selection = _fixture_system()
    original = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))
    dg.build_dynamic_ghost_force(system, topology, _protocol())

    e_original, f_original, _d0, _g0 = _state(original, positions)
    e_zero, f_zero, derivatives_zero, ghost_zero = _state(system, positions, g=0.0)
    e_one, _f_one, derivatives_one, ghost_one = _state(system, positions, g=1.0)

    assert e_zero == pytest.approx(e_original, abs=1.0e-10)
    assert np.max(np.abs(f_zero - f_original)) < 1.0e-10
    assert ghost_zero == pytest.approx(0.0, abs=1.0e-12)
    assert e_one > e_zero
    assert ghost_one is not None and ghost_one > 0.0
    assert e_one - e_zero == pytest.approx(ghost_one, abs=1.0e-9)
    assert math.isfinite(float(derivatives_zero["g"]))
    assert math.isfinite(float(derivatives_one["g"]))


def test_serialization_roundtrip_preserves_force_contract_and_energy():
    mm = _openmm()
    system, topology, positions, _selection = _fixture_system()
    dg.build_dynamic_ghost_force(system, topology, _protocol())
    roundtrip = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))

    before = dg.inspect_dynamic_ghost_force(system)
    after = dg.inspect_dynamic_ghost_force(roundtrip)
    e_before, f_before, _d_before, g_before = _state(system, positions, g=1.0)
    e_after, f_after, _d_after, g_after = _state(roundtrip, positions, g=1.0)

    assert after == before
    assert e_after == pytest.approx(e_before, abs=1.0e-10)
    assert g_after == pytest.approx(g_before, abs=1.0e-10)
    assert np.max(np.abs(f_after - f_before)) < 1.0e-10


def test_protocol_or_target_selection_drift_fails_loudly():
    system, topology, _positions, _selection = _fixture_system()
    bad_protocol = _protocol()
    bad_protocol["ghost"]["attraction"] = True
    with pytest.raises(dg.DynamicGhostError, match="attraction"):
        dg.build_dynamic_ghost_force(system, topology, bad_protocol)

    with pytest.raises(dg.DynamicGhostError, match="exactly one disappearing TRP"):
        dg.select_ghost_atoms(topology, binder_chain="Z")


def test_atm_nested_builder_preserves_top_level_scope_and_transformations():
    system, topology, positions, selection, offset = _nested_fixture_system()
    _atm_index, atm = dg.find_atm_force(system)
    n_top_level = system.getNumForces()
    n_nested = atm.getNumForces()

    report = dg.build_atm_nested_dynamic_ghost_force(
        system,
        topology,
        _nested_protocol(),
    )
    transformations = dg.inspect_atm_coordinate_transformations(
        system,
        topology,
        positions_nm=np.asarray(
            positions.value_in_unit(_openmm().unit.nanometer),
            dtype=float,
        ),
    )

    assert report["n_top_level_forces_before"] == n_top_level
    assert report["n_top_level_forces_after"] == n_top_level
    assert report["n_atm_nested_forces_before"] == n_nested
    assert report["n_atm_nested_forces_after"] == n_nested + 1
    assert report["nested_force_index"] == n_nested
    assert report["force_name"] == dg.ATM_NESTED_GHOST_FORCE_NAME
    assert report["force_group"] == 0
    assert report["net_charge_delta_e"] == 0.0
    assert transformations["n_ring_atoms"] == 9
    assert transformations["n_water_oxygens"] == 2
    assert transformations["common_u1_offset"]["norm_nm"] == pytest.approx(
        np.linalg.norm(offset),
        abs=1.0e-12,
    )
    group = report["interaction_groups"][0]
    assert group["first"] == sorted(selection.ring_indices)
    assert group["second"] == sorted(selection.water_oxygen_indices)


def test_atm_nested_ghost_tracks_stored_and_transformed_endpoint_coordinates():
    mm = _openmm()
    from openmm import unit

    system, topology, positions, selection, offset = _nested_fixture_system()
    dg.build_atm_nested_dynamic_ghost_force(system, topology, _nested_protocol())

    stored_g0 = _nested_state(system, positions, endpoint=0.0, g=0.0)
    stored_u0 = _nested_state(system, positions, endpoint=0.0, g=1.0)
    stored_u1 = _nested_state(system, positions, endpoint=1.0, g=1.0)
    assert stored_g0[0] == pytest.approx(0.0, abs=1.0e-12)
    assert stored_u0[0] > 0.1
    assert stored_u1[0] == pytest.approx(0.0, abs=1.0e-12)
    assert stored_u0[3]["u0"] > 0.1
    assert stored_u0[3]["u1"] == pytest.approx(0.0, abs=1.0e-12)
    assert math.isfinite(stored_u0[2][dg.GHOST_GLOBAL_PARAMETER])

    transformed_positions = np.asarray(
        positions.value_in_unit(unit.nanometer),
        dtype=float,
    )
    transformed_positions[selection.water_oxygen_indices[0]] += offset
    transformed_positions *= unit.nanometer
    transformed_u0 = _nested_state(
        system,
        transformed_positions,
        endpoint=0.0,
        g=1.0,
    )
    transformed_u1 = _nested_state(
        system,
        transformed_positions,
        endpoint=1.0,
        g=1.0,
    )
    assert transformed_u0[0] == pytest.approx(0.0, abs=1.0e-12)
    assert transformed_u1[0] > 0.1
    assert transformed_u1[3]["u0"] == pytest.approx(0.0, abs=1.0e-12)
    assert transformed_u1[3]["u1"] > 0.1


def test_atm_nested_serialization_preserves_contract_and_endpoint_energy():
    mm = _openmm()
    system, topology, positions, _selection, _offset = _nested_fixture_system()
    dg.build_atm_nested_dynamic_ghost_force(system, topology, _nested_protocol())
    roundtrip = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))

    before = dg.inspect_atm_nested_dynamic_ghost_force(system)
    after = dg.inspect_atm_nested_dynamic_ghost_force(roundtrip)
    expected = _nested_state(system, positions, endpoint=0.0, g=1.0)
    observed = _nested_state(roundtrip, positions, endpoint=0.0, g=1.0)

    assert after == before
    assert observed[0] == pytest.approx(expected[0], abs=1.0e-10)
    assert observed[3] == pytest.approx(expected[3], abs=1.0e-10)
    assert np.max(np.abs(observed[1] - expected[1])) < 1.0e-10
