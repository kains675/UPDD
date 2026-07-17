"""Contract tests for the frozen apex-bridge P0v2 runner."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = (
    ROOT
    / "analysis/dynamic_ghost_apex_bridge_p0v2_20260717"
    / "run_p0v2_reference_audit.py"
)


def load_runner():
    spec = importlib.util.spec_from_file_location("apex_bridge_p0v2_runner", RUNNER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_cell_matrix_and_execution_boundary_are_exact():
    runner = load_runner()
    protocol, _parent_protocol = runner.load_protocol_and_verify_freeze()

    assert [
        (cell.leg, cell.seed, cell.replicate_index) for cell in runner.make_cells()
    ] == [
        ("bound", "s101", 5),
        ("bound", "s127", 6),
        ("bound", "s163", 7),
        ("free", "s101", 5),
        ("free", "s127", 6),
        ("free", "s163", 7),
    ]
    assert protocol["bridge"]["audit_states"] == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert protocol["stages"]["P0v2"] == {
        "platform": "Reference",
        "cells": 6,
        "legs": ["bound", "free"],
        "md_allowed": False,
        "minimization_allowed": False,
        "gpu_allowed": False,
        "free_energy_estimation_allowed": False,
        "cell_subprocess_required": True,
    }


def test_worker_command_uses_current_interpreter_without_gpu_override():
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


def test_source_cells_resolve_distinct_existing_apex_directions():
    runner = load_runner()
    _protocol, parent_protocol = runner.load_protocol_and_verify_freeze()
    cell = runner.cell_by_id("w4a_uncarved_s127_free")

    record = runner.parent_p0.validate_source_cell(cell, parent_protocol)

    assert record["states"]["apex_dplus"]["direction"] == 1
    assert record["states"]["apex_dminus"]["direction"] == -1
    assert record["states"]["apex_dplus"]["lambda1"] == 0.5
    assert record["states"]["apex_dminus"]["lambda2"] == 0.5


def test_parent_stops_on_existing_bridge_reject(monkeypatch):
    runner = load_runner()
    cell = runner.make_cells()[0]
    rejected = {
        "cell": {"cell_id": cell.cell_id},
        "status": "REJECT_BRIDGE",
    }
    states = []

    monkeypatch.setattr(
        runner, "validate_keeper_preflight", lambda inventory: {"status": "PASS"}
    )
    monkeypatch.setattr(
        runner,
        "validate_completed_cell",
        lambda candidate, inventory: rejected if candidate == cell else None,
    )
    monkeypatch.setattr(
        runner,
        "write_summary",
        lambda results, inventory: {"status": "REJECT_BRIDGE"},
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

    assert runner.run_parent({"inventory_digest": "digest"}) == {
        "status": "REJECT_BRIDGE"
    }
    assert states[0][0] == "REJECT_BRIDGE"

