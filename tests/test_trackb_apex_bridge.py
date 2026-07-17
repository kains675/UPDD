"""Reference tests for the explicit Track B apex Hamiltonian bridge."""

from __future__ import annotations

import copy

import numpy as np
import pytest

import trackb_apex_bridge as bridge


def _openmm():
    return pytest.importorskip("openmm")


def _protocol() -> dict:
    return {
        "status": "FROZEN_BEFORE_BRIDGE_IMPLEMENTATION_OR_P0V2_OUTPUT",
        "bridge": {
            "parameter": "BridgeXi",
            "parameter_default": 0.0,
            "domain": [0.0, 1.0],
            "audit_states": [0.0, 0.25, 0.5, 0.75, 1.0],
            "derivative_audit_states": [0.25, 0.5, 0.75],
            "ghost_g": 0.5,
            "interpolation": "(1-BridgeXi)*Hplus+BridgeXi*Hminus",
            "inactive_legacy_globals": list(bridge.INACTIVE_SOURCE_GLOBALS),
            "active_source_globals": list(bridge.ACTIVE_SOURCE_GLOBALS),
        },
        "sampled_bridge": {"sampling_launch_authorized": False},
        "stages": {
            "P0v2": {
                "platform": "Reference",
                "cells": 6,
                "legs": ["bound", "free"],
                "md_allowed": False,
                "minimization_allowed": False,
                "gpu_allowed": False,
                "free_energy_estimation_allowed": False,
                "cell_subprocess_required": True,
            }
        },
    }


def _source_system():
    mm = _openmm()
    from openmm import unit

    system = mm.System()
    system.addParticle(12.0)
    nested = mm.CustomExternalForce("0.5*k*x^2")
    nested.addGlobalParameter("k", 12000.0)
    nested.addParticle(0, [])

    atm = mm.ATMForce(bridge.SOURCE_ATM_EXPRESSION)
    atm.addForce(copy.copy(nested))
    atm.addParticle(mm.Vec3(1.0, 0.0, 0.0), mm.Vec3(0.0, 0.0, 0.0))
    defaults = {
        "Lambda1": 0.5,
        "Lambda2": 0.5,
        "Alpha": 0.1,
        "Uh": 0.0,
        "W0": 0.0,
        "Umax": 836.8,
        "Ubcore": 418.4,
        "Acore": 0.0625,
        "Direction": 1.0,
        "UOffset": 0.0,
    }
    for name in bridge.EXPECTED_SOURCE_GLOBALS:
        atm.addGlobalParameter(name, defaults[name])
    system.addForce(atm)
    return system, [mm.Vec3(0.1, 0.0, 0.0)] * unit.nanometer


def _evaluate(system, positions, *, direction=None, xi=None):
    mm = _openmm()
    from openmm import unit

    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    if direction is not None:
        context.setParameter("Direction", float(direction))
    if xi is not None:
        context.setParameter(bridge.BRIDGE_PARAMETER, float(xi))
    state = context.getState(
        getEnergy=True,
        getForces=True,
        getParameterDerivatives=xi is not None,
    )
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    forces = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=float,
    )
    derivatives = dict(state.getEnergyParameterDerivatives()) if xi is not None else {}
    del context, integrator
    return energy, forces, derivatives


def test_builder_changes_only_declared_atm_fields():
    system, _positions = _source_system()

    report = bridge.build_apex_bridge(system, _protocol())

    assert all(report["checks"].values())
    assert report["after"]["atm_name"] == bridge.BRIDGE_FORCE_NAME
    assert report["after"]["energy_function"] == bridge.BRIDGE_ATM_EXPRESSION
    assert report["after"]["energy_parameter_derivatives"] == ["BridgeXi"]
    assert report["before"]["nested_forces"] == report["after"]["nested_forces"]
    assert report["before"]["particle_contract"] == report["after"]["particle_contract"]


def test_endpoints_interior_forces_and_derivative_are_exact():
    mm = _openmm()
    source, positions = _source_system()
    plus_energy, plus_forces, _ = _evaluate(source, positions, direction=1.0)
    minus_energy, minus_forces, _ = _evaluate(source, positions, direction=-1.0)
    assert abs(minus_energy - plus_energy) > 1.0

    bridged = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(source))
    bridge.build_apex_bridge(bridged, _protocol())
    gap = minus_energy - plus_energy
    for xi in (0.0, 0.25, 0.5, 0.75, 1.0):
        energy, forces, derivatives = _evaluate(bridged, positions, xi=xi)
        expected_energy = (1.0 - xi) * plus_energy + xi * minus_energy
        expected_forces = (1.0 - xi) * plus_forces + xi * minus_forces
        assert energy == pytest.approx(expected_energy, abs=1.0e-9)
        assert np.max(np.abs(forces - expected_forces)) < 1.0e-8
        if 0.0 < xi < 1.0:
            assert float(derivatives["BridgeXi"]) == pytest.approx(gap, abs=1.0e-9)


def test_serialization_roundtrip_preserves_bridge_contract_and_values():
    mm = _openmm()
    system, positions = _source_system()
    bridge.build_apex_bridge(system, _protocol())
    roundtrip = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))

    assert bridge.inspect_apex_bridge(roundtrip) == bridge.inspect_apex_bridge(system)
    before = _evaluate(system, positions, xi=0.5)
    after = _evaluate(roundtrip, positions, xi=0.5)
    assert after[0] == pytest.approx(before[0], abs=1.0e-10)
    assert np.max(np.abs(after[1] - before[1])) < 1.0e-10
    assert float(after[2]["BridgeXi"]) == pytest.approx(
        float(before[2]["BridgeXi"]), abs=1.0e-10
    )


def test_source_or_protocol_drift_fails_loudly():
    system, _positions = _source_system()
    _, atm = bridge.find_atm_force(system)
    atm.setEnergyFunction("u0")
    with pytest.raises(bridge.ApexBridgeError, match="source ATM energy expression"):
        bridge.build_apex_bridge(system, _protocol())

    system, _positions = _source_system()
    bad_protocol = _protocol()
    bad_protocol["bridge"]["ghost_g"] = 0.4
    with pytest.raises(bridge.ApexBridgeError, match="ghost_g"):
        bridge.build_apex_bridge(system, bad_protocol)

