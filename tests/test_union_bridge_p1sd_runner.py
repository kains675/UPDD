#!/usr/bin/env python3
"""Static and pure-numeric gates for the frozen P1 schedule-discovery runner."""

from __future__ import annotations

import hashlib
import struct

import numpy as np
import pytest

from analysis.dynamic_ghost_union_bridge_p1sd_20260717 import run_p1sd as runner


def test_freeze_is_six_cell_fixed_volume_schedule_discovery():
    protocol = runner.load_protocol_and_verify_freeze()
    assert protocol["schema"] == "updd_dynamic_ghost_union_bridge_p1sd_v1_r2"
    assert protocol["cross_state_matrix"]["slope_evaluation"]["probe_xi"] == 0.5
    assert protocol["cohort"] == [
        "w4a_union_s101_bound",
        "w4a_union_s127_bound",
        "w4a_union_s163_bound",
        "w4a_union_s101_free",
        "w4a_union_s127_free",
        "w4a_union_s163_free",
    ]
    assert protocol["sampling"]["xi"] == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert protocol["sampling"]["barostat"] is False
    assert protocol["sampling"]["one_live_context"] is True
    assert protocol["sampling"]["final_step"] == 14500
    assert protocol["mbar"]["production_estimate"] is False
    assert protocol["mbar"]["pilot_samples_reusable_in_production"] is False


def test_seed_derivation_is_exact_and_role_scoped():
    cell_id = "w4a_union_s127_bound"
    xi = 0.25
    for role in ("integrator", "velocity"):
        key = f"P1SD|{cell_id}|xi=0.25|{role}"
        expected = int(hashlib.sha256(key.encode("ascii")).hexdigest()[:8], 16)
        expected &= 0x7FFFFFFF
        assert runner.seed_unit_key(cell_id, xi, role) == key
        assert runner.derive_seed(cell_id, xi, role) == (expected or 1)
    assert runner.derive_seed(cell_id, xi, "integrator") != runner.derive_seed(
        cell_id, xi, "velocity"
    )
    with pytest.raises(ValueError):
        runner.derive_seed(cell_id, xi, "other")


def test_rigid_anchor_remap_preserves_solute_and_internal_vectors():
    positions = np.asarray(
        [
            [1.0, 1.0, 1.0],
            [8.0, 5.0, 5.0],
            [8.1, 5.0, 5.0],
            [8.0, 5.1, 5.0],
            [2.0, 7.0, 4.0],
        ],
        dtype=float,
    )
    old_box = np.diag([10.0, 10.0, 10.0])
    new_box = np.diag([9.0, 9.0, 9.0])
    remapped = runner.remap_rigid_residue_anchors(
        positions, old_box, new_box, [(1, [1, 2, 3]), (4, [4])]
    )
    assert np.array_equal(remapped[0], positions[0])
    assert np.allclose(remapped[1], [7.2, 4.5, 4.5], atol=1e-15)
    assert np.allclose(remapped[4], [1.8, 6.3, 3.6], atol=1e-15)
    assert np.allclose(remapped[[2, 3]] - remapped[1], positions[[2, 3]] - positions[1])


def test_cross_state_matrix_uses_exact_linear_bridge_identity():
    rows = [
        (
            0.0,
            [
                {"potential_kj_mol": 10.0, "bridge_slope_kj_mol": 4.0},
                {"potential_kj_mol": 11.0, "bridge_slope_kj_mol": -2.0},
            ],
        ),
        (
            1.0,
            [
                {"potential_kj_mol": 20.0, "bridge_slope_kj_mol": 8.0},
            ],
        ),
    ]
    xi = [0.0, 0.5, 1.0]
    u_kn, n_k = runner.cross_state_reduced_potentials(rows, xi, 300.0)
    beta = 1.0 / (0.00831446261815324 * 300.0)
    assert n_k.tolist() == [2, 1]
    assert u_kn.shape == (3, 3)
    assert np.allclose(u_kn[:, 0], beta * np.asarray([0.0, 2.0, 4.0]))
    assert np.allclose(u_kn[:, 1], beta * np.asarray([0.0, -1.0, -2.0]))
    assert np.allclose(u_kn[:, 2], beta * np.asarray([0.0, 4.0, 8.0]))


def test_midpoint_rule_densifies_only_failed_intervals():
    initial = [0.0, 0.25, 0.5, 0.75, 1.0]
    assert runner.propose_densified_grid(initial, [0.2, 0.09, 0.1, 0.01]) == [
        0.0,
        0.25,
        0.375,
        0.5,
        0.75,
        0.875,
        1.0,
    ]
    assert runner.propose_densified_grid(initial, [0.1, 0.2, 0.3, 0.4]) == initial
    with pytest.raises(ValueError):
        runner.propose_densified_grid(initial, [0.1])


def test_openmm_dcd_header_contract(tmp_path):
    path = tmp_path / "trajectory.dcd"
    path.write_bytes(struct.pack("<i4s3i", 84, b"CORD", 50, 2250, 250))
    assert runner._read_openmm_dcd_header(path) == {
        "frame_count": 50,
        "first_step": 2250,
        "interval": 250,
    }
    path.write_bytes(b"short")
    with pytest.raises(ValueError):
        runner._read_openmm_dcd_header(path)


def test_mixed_precision_readback_tolerance_is_bounded():
    config = runner.load_protocol_and_verify_freeze()["cross_state_matrix"][
        "explicit_readback_tolerance"
    ]
    assert runner.readback_tolerance_kj_mol(-5.0, -5.0, config) == 1e-5
    assert runner.readback_tolerance_kj_mol(-500000.0, -500000.0, config) == 5e-4
    assert runner.readback_tolerance_kj_mol(-2000000.0, -2000000.0, config) == 1e-3
    with pytest.raises(ValueError):
        runner.readback_tolerance_kj_mol(
            0.0,
            0.0,
            {
                "absolute_floor_kj_mol": 1e-2,
                "relative_to_energy_scale": 1e-9,
                "absolute_ceiling_kj_mol": 1e-3,
            },
        )


def test_cli_stages_are_mutually_exclusive():
    assert runner._parse_cli(["--build-inventory"]).build_inventory is True
    assert runner._parse_cli(["--worker-window", "x", "--xi", "0.5"]).xi == 0.5
    with pytest.raises(SystemExit):
        runner._parse_cli(["--build-inventory", "--run"])
