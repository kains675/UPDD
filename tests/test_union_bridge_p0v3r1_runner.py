#!/usr/bin/env python3
"""Static gates for the frozen union/bridge P0v3r1 runner."""

from __future__ import annotations

from analysis.dynamic_ghost_union_bridge_p0v3r1_20260717 import (
    run_p0v3r1_reference_audit as runner,
)


def test_revised_freeze_resolves_exact_six_cells():
    protocol, base, _bridge, _ghost = runner.load_protocol_and_verify_freeze()
    cells = runner.make_cells(base)
    assert protocol["override"]["static_total_box_density_gate"] == "REPORT_ONLY"
    assert protocol["override"]["equilibrium_density_claim_allowed"] is False
    assert len(cells) == 6
    assert [cell.cell_id for cell in cells] == [
        "w4a_union_s101_bound",
        "w4a_union_s127_bound",
        "w4a_union_s163_bound",
        "w4a_union_s101_free",
        "w4a_union_s127_free",
        "w4a_union_s163_free",
    ]


def test_cell_configs_preserve_component_and_parent_counts():
    for cell in runner.make_cells():
        assert cell.config["water"] + cell.config["na"] + cell.config["cl"] == cell.config[
            "num_added"
        ]
        assert cell.config["parent_particles"] > cell.config["num_added"]
        assert cell.parent_cell.seed == cell.seed
        assert cell.parent_cell.leg == cell.leg
        assert cell.parent_cell.replicate_index == cell.replicate_index


def test_inventory_digests_ignore_only_generation_metadata():
    first = {"schema": "x", "generated": "one", "value": 3}
    second = {"schema": "x", "generated": "two", "value": 3}
    assert runner.inventory_digest(first) == runner.inventory_digest(second)

    first["inventory_digest"] = "old"
    second["inventory_digest"] = "new"
    assert runner.inventory_digest(first) == runner.inventory_digest(second)


def test_cli_stages_are_mutually_exclusive():
    parsed = runner.parse_args(["--build-sources"])
    assert parsed.build_sources is True
    assert parsed.run is False
