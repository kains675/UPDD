"""Contract tests for the frozen dynamic-ghost P0 Reference runner."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = ROOT / "analysis/dynamic_ghost_excluded_volume_20260716/run_p0_reference_audit.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("dynamic_ghost_p0_runner", RUNNER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_cell_matrix_is_exact_frozen_cohort():
    runner = load_runner()

    cells = runner.make_cells()

    assert [(cell.leg, cell.seed, cell.replicate_index) for cell in cells] == [
        ("bound", "s101", 5),
        ("bound", "s127", 6),
        ("bound", "s163", 7),
        ("free", "s101", 5),
        ("free", "s127", 6),
        ("free", "s163", 7),
    ]


def test_worker_command_is_current_interpreter_and_has_no_platform_override():
    runner = load_runner()
    cell = runner.make_cells()[0]

    command = runner.worker_command(cell, "digest-123")

    assert command == [
        sys.executable,
        str(RUNNER_PATH),
        "--worker-cell",
        cell.cell_id,
        "--expected-inventory-digest",
        "digest-123",
    ]
    assert "CUDA" not in command
    assert "--platform" not in command


def test_real_source_schedules_resolve_frozen_endpoints_and_apexes():
    runner = load_runner()
    protocol, _frozen = runner.load_protocol_and_verify_freeze()
    cell = runner.cell_by_id("w4a_uncarved_s127_free")

    record = runner.validate_source_cell(cell, protocol)

    assert record["main_path_g"] == {
        "u0_dplus": 0.0,
        "apex_dplus": 0.5,
        "apex_dminus": 0.5,
        "u1_dminus": 1.0,
    }
    assert record["states"]["u0_dplus"]["direction"] == 1
    assert record["states"]["u1_dminus"]["direction"] == -1
    assert record["states"]["apex_dplus"]["lambda1"] == 0.5
    assert record["states"]["apex_dminus"]["lambda2"] == 0.5


def test_compact_force_report_keeps_exact_selection_digests():
    runner = load_runner()
    report = {
        "selection": {"water_oxygen_indices_sha_input": [10, 20, 30]},
        "interaction_groups": [{"first": [1, 2], "second": [10, 20, 30]}],
    }

    compact = runner.compact_force_report(report)

    assert "water_oxygen_indices_sha_input" not in compact["selection"]
    assert compact["selection"]["water_oxygen_index_first"] == 10
    assert compact["selection"]["water_oxygen_index_last"] == 30
    assert compact["interaction_groups"][0]["second_count"] == 3
    assert "second" not in compact["interaction_groups"][0]


def test_run_all_stops_on_existing_reject(monkeypatch):
    runner = load_runner()
    cell = runner.make_cells()[0]
    rejected = {"cell": {"cell_id": cell.cell_id}, "status": "REJECT_GHOST"}
    states = []

    monkeypatch.setattr(runner, "validate_keeper_preflight", lambda inventory: {"status": "PASS"})
    monkeypatch.setattr(
        runner,
        "validate_completed_cell",
        lambda candidate, inventory: rejected if candidate == cell else None,
    )
    monkeypatch.setattr(
        runner,
        "write_summary",
        lambda results, inventory: {"status": "REJECT_GHOST"},
    )
    monkeypatch.setattr(
        runner,
        "write_run_state",
        lambda status, **extra: states.append((status, extra)),
    )
    monkeypatch.setattr(
        runner,
        "launch_worker",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("must not launch")),
    )

    assert runner.run_all({"inventory_digest": "digest"}) == 2
    assert states[0][0] == "REJECT_GHOST"
