"""Tests for the frozen W4A post-densify DCD cohort gate."""

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
GATE_PATH = ROOT / "analysis/w4a_postdensify_dcd_20260716/postdensify_gate.py"


def load_gate():
    spec = importlib.util.spec_from_file_location("w4a_postdensify_gate", GATE_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_summary(*, control: float, high: float, uncarved: float = 0.05):
    gate = load_gate()
    cells = []
    for arm in ("carved", "uncarved"):
        for seed in gate.SEEDS:
            if arm == "carved":
                value = high if seed == gate.HIGH_LEVERAGE_SEED else control
            else:
                value = uncarved
            cells.append(
                {
                    "cell_id": f"w4a_{arm}_{seed}_bound_dminus",
                    "arm": arm,
                    "seed": seed,
                    "leg": "bound",
                    "direction": "dminus",
                    "n_frames": 750,
                    "nan_any": False,
                    "n_local_water_oxygens": 4,
                    "water_min_nm": {"frac_lt_0p26": value},
                }
            )
    return {"status": "COMPLETE", "source_inventory_digest": "test", "cells": cells}


def test_control_baseline_promotes_a():
    gate = load_gate()
    result = gate.evaluate(make_summary(control=0.15, high=0.15))
    assert result["criteria"]["control_baseline_persists"] is True
    assert result["verdict"] == "PROMOTE_A_CARVE_REDESIGN"


def test_high_leverage_amplification_promotes_a_below_baseline():
    gate = load_gate()
    result = gate.evaluate(make_summary(control=0.10, high=0.20))
    assert result["criteria"]["control_baseline_persists"] is False
    assert result["criteria"]["high_leverage_amplification_persists"] is True
    assert result["verdict"] == "PROMOTE_A_CARVE_REDESIGN"


def test_no_trigger_continues_b_transport_pilot():
    gate = load_gate()
    result = gate.evaluate(make_summary(control=0.10, high=0.14))
    assert result["criteria"]["high_leverage_amplification_persists"] is False
    assert result["verdict"] == "CONTINUE_B_TRANSPORT_PILOT"


def test_zero_control_median_does_not_define_ratio_trigger():
    gate = load_gate()
    result = gate.evaluate(make_summary(control=0.0, high=0.20))
    assert result["criteria"]["high_over_control_median"] is None
    assert result["criteria"]["high_leverage_amplification_persists"] is False
    assert result["verdict"] == "CONTINUE_B_TRANSPORT_PILOT"


def test_incomplete_summary_is_rejected():
    gate = load_gate()
    summary = make_summary(control=0.1, high=0.2)
    summary["status"] = "RUNNING"
    with pytest.raises(ValueError, match="expected 'COMPLETE'"):
        gate.evaluate(summary)
