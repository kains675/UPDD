"""Tests for W4A post-densify DCD per-cell process isolation."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = ROOT / "analysis/w4a_postdensify_dcd_20260716/postdensify_dcd.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("w4a_postdensify_dcd_runner", RUNNER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_worker_command_uses_current_interpreter_and_frozen_digest():
    runner = load_runner()
    cell = SimpleNamespace(cell_id="w4a_carved_s7_bound_dminus")

    command = runner.cell_worker_command(cell, "CUDA", "digest-123")

    assert command == [
        runner.sys.executable,
        str(RUNNER_PATH.resolve()),
        "--worker-cell",
        cell.cell_id,
        "--expected-inventory-digest",
        "digest-123",
        "--platform",
        "CUDA",
    ]


def test_worker_revalidates_source_and_runs_exactly_one_cell(monkeypatch, tmp_path):
    runner = load_runner()
    runner.OUT_ROOT = tmp_path / "out"
    cell = SimpleNamespace(cell_id="w4a_carved_s7_bound_dminus", seed="s7")
    inventory = {"inventory_digest": "digest-123"}
    events = []
    written = {}

    class FakeProd:
        def _load_rbfe(self):
            events.append("load-rbfe")
            return "rbfe"

    class FakeD2:
        def load_prod(self):
            events.append("load-prod")
            return FakeProd()

        def run_cell(self, prod, rbfe, selected_cell, *, platform):
            events.append(("run-cell", selected_cell.cell_id, rbfe, platform))
            return {"cell": {"cell_id": selected_cell.cell_id}}

    monkeypatch.setattr(runner, "load_completed", lambda *args: None)
    monkeypatch.setattr(
        runner,
        "verify_cell_source_inventory",
        lambda selected_cell, frozen: events.append(("verify", selected_cell.cell_id, frozen["inventory_digest"])),
    )
    monkeypatch.setattr(
        runner,
        "add_context",
        lambda manifest, frozen: {**manifest, "digest": frozen["inventory_digest"]},
    )
    monkeypatch.setattr(runner, "manifest_path", lambda *args: tmp_path / "manifest.json")
    monkeypatch.setattr(runner, "atomic_write_json", lambda path, payload: written.update({path: payload}))
    monkeypatch.setattr(
        runner,
        "validate_completed_manifest",
        lambda d2, selected_cell, manifest, frozen: events.append(("validate", selected_cell.cell_id)),
    )

    runner.run_cell_worker(FakeD2(), (cell,), inventory, cell.cell_id, "CUDA")

    assert events == [
        ("verify", cell.cell_id, "digest-123"),
        "load-prod",
        "load-rbfe",
        ("run-cell", cell.cell_id, "rbfe", "CUDA"),
        ("validate", cell.cell_id),
    ]
    assert written[tmp_path / "manifest.json"]["digest"] == "digest-123"


def test_parent_launches_worker_without_loading_openmm_stack(monkeypatch):
    runner = load_runner()
    cell = SimpleNamespace(cell_id="w4a_carved_s7_bound_dminus", seed="s7")
    inventory = {"inventory_digest": "digest-123"}
    completed_manifest = {"cell": {"cell_id": cell.cell_id}}
    state = {"launched": False}
    run_states = []
    summaries = []

    def fake_load_completed(d2, selected_cell, frozen):
        return completed_manifest if state["launched"] else None

    def fake_launch(selected_cell, platform, digest):
        assert selected_cell is cell
        assert platform == "CUDA"
        assert digest == "digest-123"
        state["launched"] = True

    monkeypatch.setattr(runner, "load_completed", fake_load_completed)
    monkeypatch.setattr(runner, "launch_cell_worker", fake_launch)
    monkeypatch.setattr(runner, "archive_partial_cell", lambda selected_cell: None)
    monkeypatch.setattr(runner, "write_run_state", lambda status, **extra: run_states.append((status, extra)))
    monkeypatch.setattr(
        runner,
        "write_summary",
        lambda manifests, cells, frozen: summaries.append([row["cell"]["cell_id"] for row in manifests]),
    )

    runner.run_all(object(), (cell,), inventory, "CUDA")

    assert state["launched"] is True
    assert summaries == [[], [cell.cell_id]]
    assert [status for status, _ in run_states] == ["RUNNING", "COMPLETE"]


def test_parent_stops_without_retry_when_worker_fails(monkeypatch):
    runner = load_runner()
    cell = SimpleNamespace(cell_id="w4a_carved_s7_bound_dminus", seed="s7")
    inventory = {"inventory_digest": "digest-123"}
    run_states = []
    launch_count = 0

    def fail_launch(*args):
        nonlocal launch_count
        launch_count += 1
        raise RuntimeError("worker exited with code 9")

    monkeypatch.setattr(runner, "load_completed", lambda *args: None)
    monkeypatch.setattr(runner, "launch_cell_worker", fail_launch)
    monkeypatch.setattr(runner, "archive_partial_cell", lambda selected_cell: None)
    monkeypatch.setattr(runner, "write_summary", lambda *args: None)
    monkeypatch.setattr(runner, "write_run_state", lambda status, **extra: run_states.append((status, extra)))

    with pytest.raises(RuntimeError, match="worker exited with code 9"):
        runner.run_all(object(), (cell,), inventory, "CUDA")

    assert launch_count == 1
    assert [status for status, _ in run_states] == ["RUNNING", "FAILED"]
    assert "worker exited with code 9" in run_states[-1][1]["error"]


def test_pause_is_acknowledged_only_at_cell_boundary(monkeypatch):
    runner = load_runner()
    seen = []
    monkeypatch.setattr(
        runner,
        "boundary_pause_requested",
        lambda progress: seen.append(("token", progress)) or True,
    )
    monkeypatch.setattr(
        runner,
        "write_run_state",
        lambda status, **extra: seen.append((status, extra)),
    )

    paused = runner._pause_at_cell_boundary([{"cell": {"cell_id": "done"}}], (1, 2, 3))

    assert paused is True
    assert seen[0] == ("token", {"boundary": "dcd_cell", "completed": 1, "expected": 3})
    assert seen[1] == (
        "PAUSED",
        {"n_completed": 1, "n_expected": 3, "boundary": "dcd_cell"},
    )


def test_worker_cell_source_guard_detects_drift(monkeypatch):
    runner = load_runner()
    cell = SimpleNamespace(cell_id="w4a_carved_s7_bound_dminus", system="w4a_carved")
    arm_record = {"config": {"carve_void_waters": True}}
    source_record = {"cell_id": cell.cell_id, "seed": "s7"}
    inventory = {"arms": {"carved": arm_record}, "cells": [source_record]}

    monkeypatch.setattr(runner, "validate_arm_preregistration", lambda arm: arm_record)
    monkeypatch.setattr(runner, "validate_source_cell", lambda selected_cell: source_record)
    runner.verify_cell_source_inventory(cell, inventory)

    monkeypatch.setattr(
        runner,
        "validate_source_cell",
        lambda selected_cell: {"cell_id": selected_cell.cell_id, "seed": "s19"},
    )
    with pytest.raises(ValueError, match="source artifacts drifted"):
        runner.verify_cell_source_inventory(cell, inventory)


@pytest.mark.parametrize("content", ["{", '{"cell": {"cell_id": "partial"}}'])
def test_unaccepted_worker_manifest_is_treated_as_partial(monkeypatch, tmp_path, content):
    runner = load_runner()
    path = tmp_path / "d2_cell_manifest.json"
    path.write_text(content, encoding="utf-8")
    monkeypatch.setattr(runner, "manifest_path", lambda *args: path)
    monkeypatch.setattr(
        runner,
        "validate_completed_manifest",
        lambda *args: pytest.fail("an unaccepted worker manifest must not be validated as complete"),
    )

    assert runner.load_completed(object(), object(), {"inventory_digest": "digest-123"}) is None
