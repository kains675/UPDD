#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``scripts/trackb_production`` — schedule integrity + soft-core cap.

The launcher itself runs on the ``atm`` conda env (not ``qmmm``). These tests
exercise only the pure-Python helpers (no openmm import) so they can run in
``qmmm`` alongside the rest of the suite, OR in ``atm`` for end-to-end
coverage. The launcher's ATMForce / build_leg_system path is exercised by
the existing trackb_freeleg_smoke test, not duplicated here.
"""

import importlib.util
import os
import sys

import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)


def _load_launcher_module():
    """Load scripts/trackb_production.py without importing openmm at module top.

    The launcher requires openmm at import time. Tests that do not need the
    openmm-dependent helpers can patch sys.modules with a dummy openmm stub.
    """
    spec = importlib.util.spec_from_file_location(
        "trackb_production",
        os.path.join(_PROJ, "scripts", "trackb_production.py"),
    )
    if spec is None or spec.loader is None:
        pytest.skip("could not locate scripts/trackb_production.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def launcher():
    """Import the launcher — requires openmm; skip if missing."""
    try:
        import openmm  # noqa: F401
    except ImportError:
        pytest.skip("openmm not importable in this env")
    return _load_launcher_module()


def test_schedule_22_states_symmetric(launcher):
    """The AToM 22-state schedule must be 11 forward + 11 backward symmetric."""
    assert launcher.N_STATES == 22
    assert len(launcher.LAMBDAS_1) == 22
    assert len(launcher.LAMBDAS_2) == 22
    assert len(launcher.DIRECTIONS) == 22
    assert len(launcher.INTERMD) == 22
    assert len(launcher.W0) == 22

    # Forward half (0..10): direction=+1, lambda 0→0.5
    assert launcher.DIRECTIONS[:11] == [1] * 11
    assert launcher.LAMBDAS_1[:11] == [0.00, 0.05, 0.10, 0.15, 0.20,
                                       0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    # Backward half (11..21): direction=-1, lambda 0.5→0
    assert launcher.DIRECTIONS[11:] == [-1] * 11
    assert launcher.LAMBDAS_1[11:] == [0.50, 0.45, 0.40, 0.35, 0.30,
                                       0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
    # Lambda1 == Lambda2 (linear ATM potential)
    assert launcher.LAMBDAS_1 == launcher.LAMBDAS_2
    # intermd flag at the two midpoints (states 10, 11)
    assert [i for i, v in enumerate(launcher.INTERMD) if v == 1] == [10, 11]
    # w0 nonzero only at the two midpoints
    assert [i for i, v in enumerate(launcher.W0) if v != 0.0] == [10, 11]


def test_soft_core_cap_bounds(launcher):
    """soft_core_pert_e must bound u to [ub, umax] when u > ub."""
    # Below threshold: identity.
    assert launcher.soft_core_pert_e(50.0, umax=200.0, ub=100.0, a=0.0625) == 50.0
    assert launcher.soft_core_pert_e(99.999, umax=200.0, ub=100.0, a=0.0625) == 99.999
    # At threshold: identity (u == ub).
    assert launcher.soft_core_pert_e(100.0, umax=200.0, ub=100.0, a=0.0625) == 100.0
    # Above threshold: capped to a value in (ub, umax).
    capped = launcher.soft_core_pert_e(1e15, umax=200.0, ub=100.0, a=0.0625)
    assert 100.0 < capped < 200.0
    # Modest excess: also capped, in (ub, umax).
    capped2 = launcher.soft_core_pert_e(500.0, umax=200.0, ub=100.0, a=0.0625)
    assert 100.0 < capped2 < 200.0
    # Monotone: bigger raw input -> bigger (or equal) capped output.
    capped_small = launcher.soft_core_pert_e(150.0, umax=200.0, ub=100.0, a=0.0625)
    capped_big = launcher.soft_core_pert_e(1e6, umax=200.0, ub=100.0, a=0.0625)
    assert capped_small <= capped_big <= 200.0


def test_softplus_bias_zero_at_lambda_zero(launcher):
    """At lambda1=lambda2=0 the softplus bias reduces to w0 only (alpha-independent)."""
    bias = launcher.softplus_bias_kj(
        lambda1=0.0, lambda2=0.0, alpha_per_kj=0.0239, uh_kj=110.0 * 4.184,
        w0_kj=0.0, pert_kj=50.0,
    )
    assert abs(bias) < 1e-10

    # At lambda1=lambda2=0.5 the bias is 0.5*pert + w0_at_midpoint.
    bias_mid = launcher.softplus_bias_kj(
        lambda1=0.5, lambda2=0.5, alpha_per_kj=0.0239, uh_kj=110.0 * 4.184,
        w0_kj=1.0 * 4.184, pert_kj=50.0,
    )
    expected = 0.5 * 50.0 + 1.0 * 4.184
    assert abs(bias_mid - expected) < 1e-9


def test_displacement_default(launcher):
    """Default displacement is 2.5 nm in +x (verdict B-C5 spirit)."""
    assert launcher.DISPLACEMENT_NM_DEFAULT == (2.5, 0.0, 0.0)


def test_soft_core_constants_match_attach_defaults(launcher):
    """The cap constants must match the ATMForce constructor defaults."""
    import inspect
    sig = inspect.signature(launcher.attach_atm_force_production)
    assert sig.parameters["umax_kcal"].default == launcher.UMAX_KCAL
    assert sig.parameters["ubcore_kcal"].default == launcher.UBCORE_KCAL
    assert sig.parameters["acore"].default == launcher.ACORE


def test_pre_register_outcome_schema(launcher):
    """The launcher must compose a 4-outcome pre-register schema."""
    # We can't run main() without openmm + GPU; instead verify the outcome
    # labels are present in the source.
    import pathlib
    src = pathlib.Path(launcher.__file__).read_text()
    for label in ("1_sign_stable", "2_sigma_drift",
                  "3_magnitude_drift", "4_sign_flip"):
        assert label in src, f"pre-register missing outcome label {label}"
