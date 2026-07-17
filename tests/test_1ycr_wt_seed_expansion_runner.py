"""Control-boundary tests for the frozen 1YCR WT seed-expansion runner."""

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = ROOT / "analysis/1ycr_wt_seed_expansion_20260716/expand_1ycr_wt_scaffolds.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("expand_1ycr_wt_scaffolds_runner", RUNNER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pause_is_acknowledged_only_at_seed_boundary(monkeypatch):
    runner = load_runner()
    seen = []
    monkeypatch.setattr(
        runner,
        "boundary_pause_requested",
        lambda progress: seen.append(("token", progress)) or True,
    )
    monkeypatch.setattr(
        runner,
        "write_runtime_state",
        lambda status, **extra: seen.append((status, extra)),
    )

    paused = runner._pause_at_seed_boundary([{"seed": 83}, {"seed": 127}])

    assert paused is True
    assert seen[0] == (
        "token",
        {"boundary": "scaffold_seed", "completed": 2, "expected": 5},
    )
    assert seen[1] == (
        "PAUSED",
        {"n_completed": 2, "n_expected": 5, "boundary": "scaffold_seed"},
    )
