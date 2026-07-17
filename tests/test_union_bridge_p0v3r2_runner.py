#!/usr/bin/env python3
"""Static gates for the frozen union/bridge P0v3r2 runner."""

from __future__ import annotations

from analysis.dynamic_ghost_union_bridge_p0v3r2_20260717 import (
    run_p0v3r2_reference_audit as runner,
)


def test_revised_freeze_resolves_exact_six_cells():
    protocol, base, _bridge, _ghost = runner.load_protocol_and_verify_freeze()
    cells = runner.make_cells(base)
    assert protocol["scival_verdict"] == "CONDITIONAL_APPROVE_P0V3R2_ONLY"
    assert protocol["appearing_h"]["retry_k"] == 5
    assert protocol["rng_scope"]["whole_builder_global_state_invariant_gate"] is False
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


def test_deterministic_appearing_h_trail_is_a_hard_source_gate():
    protocol, base, _bridge, _ghost = runner.load_protocol_and_verify_freeze()
    cell = runner.cell_by_id("w4a_union_s101_bound", base)
    frozen = protocol["appearing_h"]["cells"][cell.cell_id]
    serialized = {
        "_build": {
            "appearing_h_retry": {
                "unit_key": "s101|bound|w4a_trp_ala_res4",
                "retry_k": 5,
                "attempts_used": 1,
                "seeds_tried": frozen[:1],
                "seed_used": frozen[0],
            }
        }
    }
    assert runner._appearing_h_retry_audit(cell, serialized, protocol)["passed"]

    serialized["_build"]["appearing_h_retry"]["unit_key"] = "wrong"
    assert not runner._appearing_h_retry_audit(cell, serialized, protocol)["passed"]

def test_cli_stages_are_mutually_exclusive():
    parsed = runner.parse_args(["--build-sources"])
    assert parsed.build_sources is True
    assert parsed.run is False
