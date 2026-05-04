"""tests/test_schema_regression.py — Phase α Step 4 / Track 4 schema regression.

Validates that ``utils.audit_charge_consistency.schema_regression_check``
catches the Cron #34 outlier pattern (mmgbsa_summary.ncaa_element="none"
while target_card.binder.ncaa_code is populated), correctly passes the
normal aligned case, returns INCONCLUSIVE on legacy schema (missing
ncaa_element), and handles the MLE↔NML alias pair.

Schema regression rationale: pipeline summary writers historically
emitted the ncAA element identity inconsistently. The Cron #34 outlier
(ncaa_element="none" while binder is MTR) caused a downstream ranking
to treat that snapshot as wild-type, polluting the calibration mean.
This regression test prevents reintroduction of that pattern.

Reference: outputs/analysis/option_y_replicate_axis_20260427.md
(ncaa_element vs ncaa_code assertion).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from audit_charge_consistency import (
    SchemaRegressionError,
    schema_regression_check,
)


# ---------------------------------------------------------------------------
# Fixtures — minimal summary / target_card JSONs written to tmp_path
# ---------------------------------------------------------------------------


def _write_json(p: Path, obj: dict) -> Path:
    p.write_text(json.dumps(obj), encoding="utf-8")
    return p


@pytest.fixture
def make_summary(tmp_path):
    counter = {"i": 0}

    def _make(payload: dict) -> str:
        counter["i"] += 1
        path = tmp_path / f"mmgbsa_summary_{counter['i']}.json"
        return str(_write_json(path, payload))

    return _make


@pytest.fixture
def make_card(tmp_path):
    counter = {"i": 0}

    def _make(payload: dict) -> str:
        counter["i"] += 1
        path = tmp_path / f"target_card_{counter['i']}.json"
        return str(_write_json(path, payload))

    return _make


# ---------------------------------------------------------------------------
# Test 1: Cron #34 outlier — ncaa_element="none", card.binder.ncaa_code="MTR"
# Expectation: MUST FAIL (status='fail').
# ---------------------------------------------------------------------------


def test_cron_34_outlier_pattern_fails(make_summary, make_card):
    summary = make_summary({
        "snapshot": "snap01",
        "ncaa_element": "none",
        "deltaG_kcal_per_mol": -8.7,
    })
    card = make_card({
        "target_id": "2QKI",
        "binder": {"ncaa_code": "MTR", "sequence": "ICVVQDWGHHRCT"},
    })

    result = schema_regression_check(summary, card)
    assert result["status"] == "fail", (
        f"Cron #34 outlier (ncaa_element='none' vs binder.ncaa_code='MTR') "
        f"must trigger FAIL; got {result}"
    )
    assert result["expected"] == "MTR"
    assert result["actual"] == "none"
    assert "Cron #34" in result["reason"] or "disagrees" in result["reason"]


def test_cron_34_outlier_via_exception_wrapper(make_summary, make_card):
    """Verify the SchemaRegressionError class can wrap a fail result."""
    summary = make_summary({"snapshot": "snap01", "ncaa_element": "none"})
    card = make_card({"binder": {"ncaa_code": "MTR"}})

    result = schema_regression_check(summary, card)
    assert result["status"] == "fail"

    # Caller pattern: hard-fail by raising on a non-pass status.
    if result["status"] == "fail":
        with pytest.raises(SchemaRegressionError):
            raise SchemaRegressionError(result["reason"])


# ---------------------------------------------------------------------------
# Test 2: Normal aligned pattern — ncaa_element="MTR", card.ncaa_code="MTR"
# Expectation: MUST PASS (status='pass').
# ---------------------------------------------------------------------------


def test_normal_aligned_pattern_passes(make_summary, make_card):
    summary = make_summary({
        "snapshot": "snap01",
        "ncaa_element": "MTR",
        "deltaG_kcal_per_mol": -8.7,
    })
    card = make_card({
        "target_id": "2QKI",
        "binder": {"ncaa_code": "MTR"},
    })

    result = schema_regression_check(summary, card)
    assert result["status"] == "pass", (
        f"Aligned MTR↔MTR must PASS; got {result}"
    )
    assert result["expected"] == "MTR"
    assert result["actual"] == "MTR"


# ---------------------------------------------------------------------------
# Test 3: Missing ncaa_element field → INCONCLUSIVE (warning, not block).
# ---------------------------------------------------------------------------


def test_missing_ncaa_element_inconclusive(make_summary, make_card):
    summary = make_summary({
        "snapshot": "snap01",
        # No ncaa_element / ncaa_code field at all.
        "deltaG_kcal_per_mol": -8.7,
    })
    card = make_card({"binder": {"ncaa_code": "MTR"}})

    result = schema_regression_check(summary, card)
    assert result["status"] == "inconclusive", (
        f"Missing ncaa_element must yield INCONCLUSIVE (warning, not block); "
        f"got {result}"
    )
    assert result["expected"] == "MTR"
    assert result["actual"] is None
    assert "missing" in result["reason"].lower()


# ---------------------------------------------------------------------------
# Test 4: MLE/NML alias — handles bidirectional alias mapping.
# ---------------------------------------------------------------------------


def test_mle_nml_alias_passes_default_aliases(make_summary, make_card):
    """Default _DEFAULT_SCHEMA_ALIASES maps MLE↔NML."""
    summary = make_summary({"snapshot": "snap01", "ncaa_element": "MLE"})
    card = make_card({"binder": {"ncaa_code": "NML"}})

    result = schema_regression_check(summary, card)
    assert result["status"] == "pass", (
        f"MLE↔NML alias must PASS via default aliases; got {result}"
    )
    assert result["expected"] == "NML"
    assert result["actual"] == "MLE"
    assert "alias" in result["reason"].lower()


def test_nml_mle_alias_reverse_passes(make_summary, make_card):
    """Bidirectional: NML in summary, MLE in card."""
    summary = make_summary({"snapshot": "snap01", "ncaa_element": "NML"})
    card = make_card({"binder": {"ncaa_code": "MLE"}})

    result = schema_regression_check(summary, card)
    assert result["status"] == "pass"


def test_alias_with_explicit_dict(make_summary, make_card):
    """Caller can override allow_aliases (e.g. for custom mappings)."""
    summary = make_summary({"snapshot": "snap01", "ncaa_element": "FOO"})
    card = make_card({"binder": {"ncaa_code": "BAR"}})

    # Without an alias mapping, FOO vs BAR must FAIL.
    r1 = schema_regression_check(summary, card, allow_aliases={})
    assert r1["status"] == "fail"

    # With explicit FOO↔BAR alias, must PASS.
    r2 = schema_regression_check(summary, card, allow_aliases={"FOO": "BAR", "BAR": "FOO"})
    assert r2["status"] == "pass"
