from __future__ import annotations

import math

import numpy as np
import pytest

from analysis.raw_gap_cap_activation_p1cap_20260717 import run_p1cap as p1cap


SOFTCORE = {
    "uoffset_kj_mol": 0.0,
    "ubcore_kj_mol": 418.4,
    "umax_kj_mol": 836.8,
    "acore": 0.0625,
}


@pytest.mark.parametrize("delta", [-418.4, -100.0, 0.0, 100.0, 418.4])
def test_bridge_slope_is_zero_in_linear_branch(delta: float) -> None:
    assert p1cap.analytic_bridge_slope(delta, **SOFTCORE) == pytest.approx(
        0.0, abs=1e-12
    )


def test_bridge_slope_is_directional_when_cap_is_active() -> None:
    positive = p1cap.analytic_bridge_slope(1000.0, **SOFTCORE)
    negative = p1cap.analytic_bridge_slope(-1000.0, **SOFTCORE)
    assert positive > 0.0
    assert negative < 0.0
    assert positive == pytest.approx(-negative, rel=1e-12)


@pytest.mark.parametrize(
    ("delta", "expected"),
    [
        (0.0, "linear"),
        (413.4, "linear"),
        (413.400001, "ambiguous"),
        (418.4, "ambiguous"),
        (423.399999, "ambiguous"),
        (423.4, "cap_active"),
        (-423.4, "cap_active"),
    ],
)
def test_cap_classification_guard(delta: float, expected: str) -> None:
    assert (
        p1cap.classify_cap(
            delta, ubcore_kj_mol=418.4, guard_kj_mol=5.0
        )
        == expected
    )


def test_softcore_is_continuous_at_ubcore() -> None:
    at = p1cap.softcore_value(
        418.4,
        ubcore_kj_mol=418.4,
        umax_kj_mol=836.8,
        acore=0.0625,
    )
    above = p1cap.softcore_value(
        math.nextafter(418.4, math.inf),
        ubcore_kj_mol=418.4,
        umax_kj_mol=836.8,
        acore=0.0625,
    )
    assert at == pytest.approx(418.4)
    assert above == pytest.approx(at, abs=1e-9)


def test_ambiguous_final_classification_does_not_force_mismatch() -> None:
    assert p1cap.strict_classification_agrees("linear", "ambiguous")
    assert p1cap.strict_classification_agrees("ambiguous", "cap_active")
    assert not p1cap.strict_classification_agrees("linear", "cap_active")


def test_r1_angstrom_box_is_explicitly_converted_to_nm() -> None:
    raw_lengths_angstrom = np.asarray([145.795, 145.795, 145.795])
    angles_deg = np.asarray([90.0, 90.0, 90.0])
    box_nm = p1cap._box_vectors(raw_lengths_angstrom * 0.1, angles_deg)
    assert box_nm == pytest.approx(np.diag([14.5795, 14.5795, 14.5795]))


def test_r1_protocol_freezes_low_level_dcd_unit_contract() -> None:
    protocol = p1cap.load_protocol_and_verify_freeze()
    assert protocol["evaluation"]["dcd_distance_unit"] == "angstroms"
    assert protocol["evaluation"]["dcd_to_nm_scale"] == 0.1
