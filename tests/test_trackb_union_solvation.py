#!/usr/bin/env python3
"""Focused tests for count-preserving two-copy union solvation."""

from __future__ import annotations

import random

import pytest


pytest.importorskip("openmm")

import openmm as mm
import openmm.unit as unit
from openmm import app

from utils import atm_trackB_setup as ats
from utils import atm_trackB_inplace_rbfe as rbfe


class _SyntheticSpec:
    resnum = 4
    stateA_resname = "SRC"
    stateA_only_atoms = ("CG",)

    @staticmethod
    def _heavy_names(names):
        return tuple(names)


def _write_synthetic_forcefield(path):
    path.write_text(
        """<ForceField>
  <AtomTypes>
    <Type name="APP-C" class="APP-C" element="C" mass="12.01078"/>
    <Type name="APP-O" class="APP-O" element="O" mass="15.99943"/>
    <Type name="SRC-C" class="SRC-C" element="C" mass="12.01078"/>
    <Type name="SRC-N" class="SRC-N" element="N" mass="14.00672"/>
  </AtomTypes>
  <Residues>
    <Residue name="APP">
      <Atom name="CA" type="APP-C" charge="0.0"/>
      <Atom name="OX" type="APP-O" charge="0.0"/>
      <Bond atomName1="CA" atomName2="OX"/>
    </Residue>
    <Residue name="SRC">
      <Atom name="CG" type="SRC-C" charge="0.0"/>
      <Atom name="NX" type="SRC-N" charge="0.0"/>
      <Bond atomName1="CG" atomName2="NX"/>
    </Residue>
  </Residues>
  <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">
    <UseAttributeFromResidue name="charge"/>
    <Atom type="APP-C" sigma="0.3399669508423535" epsilon="0.4577296"/>
    <Atom type="APP-O" sigma="0.2959921901149463" epsilon="0.87864"/>
    <Atom type="SRC-C" sigma="0.3399669508423535" epsilon="0.4577296"/>
    <Atom type="SRC-N" sigma="0.325" epsilon="0.71128"/>
  </NonbondedForce>
  <HarmonicBondForce>
    <Bond class1="SRC-C" class2="SRC-N" length="0.145" k="300000.0"/>
    <Bond class1="APP-C" class2="APP-O" length="0.145" k="300000.0"/>
  </HarmonicBondForce>
</ForceField>
""",
        encoding="utf-8",
    )


def _synthetic_merged():
    topology = app.Topology()
    chain = topology.addChain("B")
    appearing = topology.addResidue("APP", chain, "4")
    ca = topology.addAtom("CA", app.element.carbon, appearing)
    ox = topology.addAtom("OX", app.element.oxygen, appearing)
    topology.addBond(ca, ox)
    disappearing = topology.addResidue("SRC", chain, "4")
    cg = topology.addAtom("CG", app.element.carbon, disappearing)
    nx = topology.addAtom("NX", app.element.nitrogen, disappearing)
    topology.addBond(cg, nx)
    positions = [
        mm.Vec3(0.0, 0.0, 0.0),
        mm.Vec3(0.0, 0.145, 0.0),
        mm.Vec3(1.0, 0.0, 0.0),
        mm.Vec3(1.145, 0.0, 0.0),
    ]
    return app.Modeller(topology, positions * unit.nanometer)


def test_union_config_is_exact_and_integer_typed():
    valid = {"num_added": 10, "water": 8, "na": 1, "cl": 1, "rng_seed": 7}
    assert ats._validate_union_config(valid) == valid

    with pytest.raises(ats.UnionSolvationError, match="keys"):
        ats._validate_union_config({**valid, "extra": 1})
    with pytest.raises(ats.UnionSolvationError, match="must be an integer"):
        ats._validate_union_config({**valid, "water": 8.0})
    with pytest.raises(ats.UnionSolvationError, match="do not sum"):
        ats._validate_union_config({**valid, "water": 7})


def test_union_and_static_carve_are_mutually_exclusive():
    with pytest.raises(ats.UnionSolvationError, match="mutually exclusive"):
        ats.build_inplace_res4_twocopy_system(
            leg="free",
            seed="s101",
            solvate=True,
            carve_void_waters=True,
            union_solvation={
                "num_added": 10,
                "water": 8,
                "na": 1,
                "cl": 1,
                "rng_seed": 7,
            },
        )


def test_serializer_forwards_union_config_only_to_twocopy(monkeypatch, tmp_path):
    config = {"num_added": 10, "water": 8, "na": 1, "cl": 1, "rng_seed": 7}
    captured = {}

    def fake_serialize(**kwargs):
        captured.update(kwargs)
        return {"union_solvation_report": {"config": kwargs["union_solvation"]}}

    monkeypatch.setattr(rbfe, "_serialize_twocopy_system", fake_serialize)
    result = rbfe.serialize_inplace_rbfe_system(
        leg="free",
        out_dir=str(tmp_path / "twocopy"),
        construction="twocopy",
        union_solvation=config,
    )
    assert captured["union_solvation"] is config
    assert result["union_solvation_report"]["config"] is config

    with pytest.raises(ValueError, match="only supported with construction='twocopy'"):
        rbfe.serialize_inplace_rbfe_system(
            leg="free",
            out_dir=str(tmp_path / "single_core"),
            construction="single_core",
            union_solvation=config,
        )


def test_union_solvation_removes_placeholders_and_restores_rng(tmp_path):
    forcefield_path = tmp_path / "synthetic.xml"
    _write_synthetic_forcefield(forcefield_path)
    forcefield_inputs = [str(forcefield_path), "amber14/tip3pfb.xml"]
    canonical = app.ForceField(*forcefield_inputs)
    merged = _synthetic_merged()
    state_before = random.getstate()

    report = ats._count_preserving_union_solvate(
        merged,
        n_copy1=2,
        dvec=(1.0, 0.0, 0.0),
        spec=_SyntheticSpec(),
        binder_chain="B",
        ff_inputs=forcefield_inputs,
        canonical_forcefield=canonical,
        config={
            "num_added": 256,
            "water": 254,
            "na": 1,
            "cl": 1,
            "rng_seed": 20260717,
        },
    )

    assert random.getstate() == state_before
    assert report["placeholder_count"] == 1
    assert report["max_placeholder_target_delta_nm"] <= 1.0e-6
    assert report["seeded_rng_state_changed"] is True
    assert report["global_rng_state_restored"] is True
    assert report["placeholders_removed"] is True
    assert report["final_solvent_counts"] == {
        "water": 254,
        "na": 1,
        "cl": 1,
        "total": 256,
    }
    assert merged.topology.getNumAtoms() == 4 + 254 * 3 + 2
    assert all(not residue.name.startswith("U") for residue in merged.topology.residues())

    final_system = canonical.createSystem(
        merged.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=None,
        rigidWater=True,
    )
    assert final_system.getNumParticles() == merged.topology.getNumAtoms()
