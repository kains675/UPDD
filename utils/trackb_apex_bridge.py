"""Explicit Hamiltonian bridge between the two Track B ATM apex states.

The bridge is limited to the CPU-only P0v2 audit preregistered under
``analysis/dynamic_ghost_apex_bridge_p0v2_20260717``.  It constructs the
Hamiltonian but does not sample it or estimate a free energy.
"""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any


BRIDGE_FORCE_NAME = "UPDDApexBridgeATM"
BRIDGE_PARAMETER = "BridgeXi"
BRIDGE_PARAMETER_DEFAULT = 0.0

SOURCE_ATM_EXPRESSION = (
    "select(step(Direction), u0, u1) + "
    "select(Lambda2-Lambda1 , "
    "((Lambda2-Lambda1)/Alpha)*log(1+exp(-Alpha*(usc-Uh))) + "
    "Lambda2*usc + W0, Lambda2*usc + W0);"
    "usc = select(Acore, select(step(u-Ubcore), "
    "(Umax-Ubcore)*fsc+Ubcore, u), u);"
    "fsc = (z^Acore-1)/(z^Acore+1);"
    "z = 1 + 2*(y/Acore) + 2*(y/Acore)^2;"
    "y = (u-Ubcore)/(Umax-Ubcore);"
    "u = select(step(Direction), 1, -1)*(u1-(u0 + UOffset))"
)

# The outer selects preserve the exact endpoint branches while the interior
# branch is the declared linear interpolation.
BRIDGE_ATM_EXPRESSION = (
    "select(BridgeXi,"
    "select(1-BridgeXi,(1-BridgeXi)*Hplus+BridgeXi*Hminus,Hminus),"
    "Hplus);"
    "Hplus=u0+0.5*uscplus;"
    "Hminus=u1+0.5*uscminus;"
    "uscplus=select(Acore,select(step(upos-Ubcore),"
    "(Umax-Ubcore)*fplus+Ubcore,upos),upos);"
    "fplus=(zplus^Acore-1)/(zplus^Acore+1);"
    "zplus=1+2*(yplus/Acore)+2*(yplus/Acore)^2;"
    "yplus=(upos-Ubcore)/(Umax-Ubcore);"
    "uscminus=select(Acore,select(step(uneg-Ubcore),"
    "(Umax-Ubcore)*fminus+Ubcore,uneg),uneg);"
    "fminus=(zneg^Acore-1)/(zneg^Acore+1);"
    "zneg=1+2*(yneg/Acore)+2*(yneg/Acore)^2;"
    "yneg=(uneg-Ubcore)/(Umax-Ubcore);"
    "upos=u1-(u0+UOffset);"
    "uneg=-1*(u1-(u0+UOffset))"
)

EXPECTED_SOURCE_GLOBALS = (
    "Lambda1",
    "Lambda2",
    "Alpha",
    "Uh",
    "W0",
    "Umax",
    "Ubcore",
    "Acore",
    "Direction",
    "UOffset",
)
ACTIVE_SOURCE_GLOBALS = ("Umax", "Ubcore", "Acore", "UOffset")
INACTIVE_SOURCE_GLOBALS = (
    "Lambda1",
    "Lambda2",
    "Alpha",
    "Uh",
    "W0",
    "Direction",
)


class ApexBridgeError(RuntimeError):
    """A frozen apex-bridge declaration or integrity gate failed."""


def _vec3_values(value: Any) -> list[float]:
    from openmm import unit

    if hasattr(value, "value_in_unit"):
        value = value.value_in_unit(unit.nanometer)
    return [float(value[0]), float(value[1]), float(value[2])]


def find_atm_force(system: Any) -> tuple[int, Any]:
    import openmm as mm

    matches = [
        (index, system.getForce(index))
        for index in range(system.getNumForces())
        if isinstance(system.getForce(index), mm.ATMForce)
    ]
    if len(matches) != 1:
        raise ApexBridgeError(
            f"expected exactly one top-level ATMForce, found {len(matches)}"
        )
    return matches[0]


def _global_parameters(atm_force: Any) -> list[dict[str, Any]]:
    return [
        {
            "name": atm_force.getGlobalParameterName(index),
            "default": float(atm_force.getGlobalParameterDefaultValue(index)),
        }
        for index in range(atm_force.getNumGlobalParameters())
    ]


def _energy_derivatives(atm_force: Any) -> list[str]:
    return [
        atm_force.getEnergyParameterDerivativeName(index)
        for index in range(atm_force.getNumEnergyParameterDerivatives())
    ]


def _nested_force_hashes(atm_force: Any) -> list[dict[str, Any]]:
    import openmm as mm

    rows = []
    for index in range(atm_force.getNumForces()):
        force = atm_force.getForce(index)
        xml = mm.XmlSerializer.serialize(force)
        rows.append(
            {
                "index": index,
                "type": type(force).__name__,
                "name": force.getName(),
                "sha256": hashlib.sha256(xml.encode("utf-8")).hexdigest(),
            }
        )
    return rows


def _particle_contract(atm_force: Any) -> dict[str, Any]:
    digest = hashlib.sha256()
    type_counts: dict[str, int] = {}
    for index in range(atm_force.getNumParticles()):
        transformation = atm_force.getParticleTransformation(index)
        transform_type = type(transformation).__name__
        type_counts[transform_type] = type_counts.get(transform_type, 0) + 1
        if transform_type == "FixedDisplacement":
            displacement1, displacement0 = atm_force.getParticleParameters(index)
            transform = {
                "type": transform_type,
                "fixed0": _vec3_values(transformation.getFixedDisplacement0()),
                "fixed1": _vec3_values(transformation.getFixedDisplacement1()),
            }
            particle_parameters = {
                "displacement1_nm": _vec3_values(displacement1),
                "displacement0_nm": _vec3_values(displacement0),
            }
        elif transform_type == "ParticleOffsetDisplacement":
            transform = {
                "type": transform_type,
                "origin0": int(transformation.getOriginParticle0()),
                "origin1": int(transformation.getOriginParticle1()),
                "destination0": int(transformation.getDestinationParticle0()),
                "destination1": int(transformation.getDestinationParticle1()),
            }
            particle_parameters = None
        else:
            raise ApexBridgeError(
                f"unsupported ATM particle transformation {transform_type!r}"
            )
        record = {
            "index": index,
            "particle_parameters": particle_parameters,
            "transformation": transform,
        }
        digest.update(
            json.dumps(
                record,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ).encode("utf-8")
        )
        digest.update(b"\n")
    return {
        "count": int(atm_force.getNumParticles()),
        "transformation_type_counts": dict(sorted(type_counts.items())),
        "sha256": digest.hexdigest(),
    }


def _other_top_level_force_hashes(system: Any, atm_index: int) -> list[dict[str, Any]]:
    import openmm as mm

    rows = []
    for index in range(system.getNumForces()):
        if index == atm_index:
            continue
        force = system.getForce(index)
        xml = mm.XmlSerializer.serialize(force)
        rows.append(
            {
                "index": index,
                "type": type(force).__name__,
                "name": force.getName(),
                "force_group": int(force.getForceGroup()),
                "sha256": hashlib.sha256(xml.encode("utf-8")).hexdigest(),
            }
        )
    return rows


def inspect_atm_contract(system: Any) -> dict[str, Any]:
    """Return an exact contract for bridge-permitted and invariant fields."""
    atm_index, atm_force = find_atm_force(system)
    return {
        "system_particles": int(system.getNumParticles()),
        "system_forces": int(system.getNumForces()),
        "atm_index": int(atm_index),
        "atm_name": atm_force.getName(),
        "atm_force_group": int(atm_force.getForceGroup()),
        "energy_function": atm_force.getEnergyFunction(),
        "global_parameters": _global_parameters(atm_force),
        "energy_parameter_derivatives": _energy_derivatives(atm_force),
        "nested_forces": _nested_force_hashes(atm_force),
        "particle_contract": _particle_contract(atm_force),
        "other_top_level_forces": _other_top_level_force_hashes(system, atm_index),
    }


def validate_bridge_protocol(protocol: dict[str, Any]) -> dict[str, Any]:
    if protocol.get("status") != "FROZEN_BEFORE_BRIDGE_IMPLEMENTATION_OR_P0V2_OUTPUT":
        raise ApexBridgeError("apex-bridge protocol is not in the frozen state")
    bridge = protocol.get("bridge") or {}
    expected = {
        "parameter": BRIDGE_PARAMETER,
        "parameter_default": BRIDGE_PARAMETER_DEFAULT,
        "domain": [0.0, 1.0],
        "audit_states": [0.0, 0.25, 0.5, 0.75, 1.0],
        "derivative_audit_states": [0.25, 0.5, 0.75],
        "ghost_g": 0.5,
        "interpolation": "(1-BridgeXi)*Hplus+BridgeXi*Hminus",
        "inactive_legacy_globals": list(INACTIVE_SOURCE_GLOBALS),
        "active_source_globals": list(ACTIVE_SOURCE_GLOBALS),
    }
    for key, value in expected.items():
        if bridge.get(key) != value:
            raise ApexBridgeError(
                f"frozen bridge protocol mismatch for {key}: "
                f"{bridge.get(key)!r} != {value!r}"
            )
    stage = (protocol.get("stages") or {}).get("P0v2") or {}
    required_stage = {
        "platform": "Reference",
        "cells": 6,
        "legs": ["bound", "free"],
        "md_allowed": False,
        "minimization_allowed": False,
        "gpu_allowed": False,
        "free_energy_estimation_allowed": False,
        "cell_subprocess_required": True,
    }
    if stage != required_stage:
        raise ApexBridgeError(f"frozen P0v2 stage drifted: {stage!r}")
    if (protocol.get("sampled_bridge") or {}).get("sampling_launch_authorized"):
        raise ApexBridgeError("P0v2 must not authorize sampled-bridge launch")
    return bridge


def _validate_source_atm(atm_force: Any) -> None:
    if atm_force.getEnergyFunction() != SOURCE_ATM_EXPRESSION:
        raise ApexBridgeError("source ATM energy expression drifted")
    globals_before = _global_parameters(atm_force)
    names = tuple(row["name"] for row in globals_before)
    if names != EXPECTED_SOURCE_GLOBALS:
        raise ApexBridgeError(f"source ATM globals drifted: {names!r}")
    if _energy_derivatives(atm_force):
        raise ApexBridgeError("source ATMForce unexpectedly declares derivatives")
    if BRIDGE_PARAMETER in names:
        raise ApexBridgeError("source ATMForce already contains BridgeXi")
    constants = {row["name"]: row["default"] for row in globals_before}
    expected_constants = {
        "UOffset": 0.0,
        "Umax": 836.8,
        "Ubcore": 418.4,
        "Acore": 0.0625,
    }
    for name, expected in expected_constants.items():
        if not math.isclose(constants[name], expected, rel_tol=0.0, abs_tol=1.0e-10):
            raise ApexBridgeError(
                f"source ATM {name}={constants[name]!r}, expected {expected!r}"
            )


def build_apex_bridge(system: Any, protocol: dict[str, Any]) -> dict[str, Any]:
    """Mutate one source ATMForce into the preregistered apex bridge."""
    validate_bridge_protocol(protocol)
    atm_index, atm_force = find_atm_force(system)
    _validate_source_atm(atm_force)
    before = inspect_atm_contract(system)

    atm_force.setEnergyFunction(BRIDGE_ATM_EXPRESSION)
    atm_force.setName(BRIDGE_FORCE_NAME)
    atm_force.addGlobalParameter(BRIDGE_PARAMETER, BRIDGE_PARAMETER_DEFAULT)
    atm_force.addEnergyParameterDerivative(BRIDGE_PARAMETER)

    after = inspect_atm_contract(system)
    expected_globals = before["global_parameters"] + [
        {"name": BRIDGE_PARAMETER, "default": BRIDGE_PARAMETER_DEFAULT}
    ]
    checks = {
        "system_particles_invariant": (
            before["system_particles"] == after["system_particles"]
        ),
        "system_forces_invariant": before["system_forces"] == after["system_forces"],
        "atm_index_invariant": before["atm_index"] == after["atm_index"] == atm_index,
        "atm_force_group_invariant": (
            before["atm_force_group"] == after["atm_force_group"]
        ),
        "source_expression_exact": before["energy_function"] == SOURCE_ATM_EXPRESSION,
        "bridge_expression_exact": after["energy_function"] == BRIDGE_ATM_EXPRESSION,
        "bridge_name_exact": after["atm_name"] == BRIDGE_FORCE_NAME,
        "globals_append_only": after["global_parameters"] == expected_globals,
        "bridge_derivative_only": (
            after["energy_parameter_derivatives"] == [BRIDGE_PARAMETER]
        ),
        "nested_forces_invariant": before["nested_forces"] == after["nested_forces"],
        "particle_transformations_invariant": (
            before["particle_contract"] == after["particle_contract"]
        ),
        "other_top_level_forces_invariant": (
            before["other_top_level_forces"] == after["other_top_level_forces"]
        ),
    }
    if not all(checks.values()):
        raise ApexBridgeError(f"apex bridge mutated forbidden fields: {checks}")
    return {
        "atm_index": int(atm_index),
        "checks": checks,
        "before": before,
        "after": after,
    }


def inspect_apex_bridge(system: Any) -> dict[str, Any]:
    """Validate and report a constructed or deserialized apex bridge."""
    _atm_index, atm_force = find_atm_force(system)
    contract = inspect_atm_contract(system)
    globals_after = contract["global_parameters"]
    names = tuple(row["name"] for row in globals_after)
    if names != EXPECTED_SOURCE_GLOBALS + (BRIDGE_PARAMETER,):
        raise ApexBridgeError(f"bridge ATM globals drifted: {names!r}")
    if globals_after[-1]["default"] != BRIDGE_PARAMETER_DEFAULT:
        raise ApexBridgeError("BridgeXi default drifted")
    if contract["atm_name"] != BRIDGE_FORCE_NAME:
        raise ApexBridgeError("bridge ATM name drifted")
    if contract["energy_function"] != BRIDGE_ATM_EXPRESSION:
        raise ApexBridgeError("bridge ATM expression drifted")
    if contract["energy_parameter_derivatives"] != [BRIDGE_PARAMETER]:
        raise ApexBridgeError("bridge derivative declaration drifted")
    return contract
