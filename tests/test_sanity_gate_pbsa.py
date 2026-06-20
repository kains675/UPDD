"""tests/test_sanity_gate_pbsa.py — numerical sign-validity sanity gate.

The gate is variant-pair (Mode B) catastrophic-cancellation only:

    r1 = abs(ΔΔnet) / (abs(ΔΔgas) + abs(ΔΔsolv))
    sign_ill_conditioned_cancellation = (r1 < c)        # c = 0.15

The per-snap / per-system endpoint level has NO gate (cancellation there is
universal and undiscriminating); ``apply_pbsa_sanity_gate`` only LOGS the
endpoint cancellation ratio.

Covers:
  (a) Cp4 canonical ΔΔ (r1 ≈ 0.069)            -> flagged True
  (b) WT-vs-WT null ΔΔ (large net, small r1)   -> flagged True (artifact)
  (c) cutoff-c boundary behaviour (0.1 / 0.15 / 0.2)
  (d) pair undefined (single-system) -> per-snap gate = N/A, r_A logged only
  (e) additive regression: pre-existing per-snap keys unchanged
  (f) synthetic large, well-resolved ΔΔ (r1 ≫ 0.2) -> passes [SYNTHETIC ONLY]
  (g) DELTA G gas / DELTA G solv parsing from a sample .dat fixture
"""
import importlib.util as _ilu
import os
import sys
from pathlib import Path

import pytest

# conftest already puts utils/ on sys.path; import the module directly.
_REPO_ROOT = Path(__file__).resolve().parents[1]
_UTILS = _REPO_ROOT / "utils"
if str(_UTILS) not in sys.path:
    sys.path.insert(0, str(_UTILS))

from sanity_gate_pbsa import (  # noqa: E402
    SANITY_GATE_CANCEL_C,
    apply_pbsa_pair_sanity_gate,
    apply_pbsa_pair_sanity_gate_from_summaries,
    apply_pbsa_sanity_gate,
)


# --------------------------------------------------------------------
# (a) Cp4 canonical ΔΔ -> flagged (catastrophic cancellation)
# --------------------------------------------------------------------
def test_a_cp4_canonical_dd_flagged():
    # Cp4(s251) - WT(s42), the canonical Cp4 case from the raw cohort:
    #   ΔΔgas = +187.20, ΔΔsolv = -162.97, ΔΔnet = +24.23
    #   r1 = 24.23 / (187.20 + 162.97) = 0.0692  < c=0.15  -> flagged.
    out = apply_pbsa_pair_sanity_gate(
        dd_gas=187.20, dd_solv=-162.97, dd_net=24.23
    )
    assert out["dd_cancellation_r1"] == pytest.approx(0.0692, abs=1e-3)
    assert out["sign_ill_conditioned_cancellation"] is True
    assert out["sanity_gate_reason"] is not None
    # r2 companion = abs(net)/max(abs(gas), abs(solv)) = 24.23/187.20.
    assert out["dd_cancellation_r2"] == pytest.approx(24.23 / 187.20, abs=1e-4)
    # Phenomenological tag only: numerical conditioning, no FF attribution.
    reason = out["sanity_gate_reason"].lower()
    assert "ill-conditioned" in reason
    assert "not a force-field attribution" in reason
    # default cutoff echoed
    assert out["sanity_gate_cancel_c"] == pytest.approx(0.15)


# --------------------------------------------------------------------
# (b) WT-vs-WT null ΔΔ (true ΔΔ = 0) -> flagged True (correctly: artifact)
# --------------------------------------------------------------------
def test_b_wt_null_pair_flagged_as_artifact():
    # Same chemistry (true ΔΔ = 0) but seed differs: ΔΔnet ≈ +9.59 with
    # r1 ≈ 0.0875 (from the raw WT_s101 x WT_s163 null pair). A large net
    # with a small cancellation ratio is exactly an ill-conditioned
    # artifact, so the gate SHOULD fire here. (ΔΔgas/ΔΔsolv chosen to
    # reproduce r1 = 9.59 / (60.0 + 49.6) = 0.0875.)
    out = apply_pbsa_pair_sanity_gate(
        dd_gas=60.0, dd_solv=-49.6, dd_net=9.59
    )
    assert out["dd_cancellation_r1"] == pytest.approx(0.0875, abs=1e-3)
    assert out["sign_ill_conditioned_cancellation"] is True


def test_b_wt_null_tiny_net_flagged():
    # The other null pair (WT_s42 x WT_s83): ΔΔnet ≈ -0.16, r1 ≈ 0.0065,
    # deep in cancellation -> flagged.
    out = apply_pbsa_pair_sanity_gate(
        dd_gas=15.0, dd_solv=-9.7, dd_net=-0.16
    )
    assert out["dd_cancellation_r1"] == pytest.approx(0.16 / 24.7, abs=1e-3)
    assert out["sign_ill_conditioned_cancellation"] is True


# --------------------------------------------------------------------
# (c) cutoff-c boundary behaviour
# --------------------------------------------------------------------
def test_c_boundary_default_is_015():
    assert SANITY_GATE_CANCEL_C == pytest.approx(0.15)


def test_c_ratio_just_below_c_fires():
    # r1 = 14.0 / (100.0 + 0.0) = 0.14 < 0.15  -> fires.
    out = apply_pbsa_pair_sanity_gate(dd_gas=100.0, dd_solv=0.0, dd_net=14.0)
    assert out["dd_cancellation_r1"] == pytest.approx(0.14, abs=1e-9)
    assert out["sign_ill_conditioned_cancellation"] is True


def test_c_ratio_at_and_above_c_does_not_fire():
    # r1 == c exactly (strict <): 15.0 / 100.0 = 0.15 -> NOT flagged.
    out_eq = apply_pbsa_pair_sanity_gate(dd_gas=100.0, dd_solv=0.0, dd_net=15.0)
    assert out_eq["dd_cancellation_r1"] == pytest.approx(0.15, abs=1e-9)
    assert out_eq["sign_ill_conditioned_cancellation"] is False
    # r1 = 0.20 > c -> not flagged.
    out_gt = apply_pbsa_pair_sanity_gate(dd_gas=100.0, dd_solv=0.0, dd_net=20.0)
    assert out_gt["dd_cancellation_r1"] == pytest.approx(0.20, abs=1e-9)
    assert out_gt["sign_ill_conditioned_cancellation"] is False


def test_c_custom_cutoff_overrides_default():
    # The same ΔΔ (r1 = 0.069) does not fire at a stricter c = 0.05, and
    # the chosen cutoff is echoed in the output.
    out = apply_pbsa_pair_sanity_gate(
        dd_gas=187.20, dd_solv=-162.97, dd_net=24.23, c=0.05
    )
    assert out["sign_ill_conditioned_cancellation"] is False
    assert out["sanity_gate_cancel_c"] == pytest.approx(0.05)
    # at c = 0.1 it would fire (0.069 < 0.1) and at c = 0.2 also.
    assert apply_pbsa_pair_sanity_gate(
        dd_gas=187.20, dd_solv=-162.97, dd_net=24.23, c=0.1
    )["sign_ill_conditioned_cancellation"] is True
    assert apply_pbsa_pair_sanity_gate(
        dd_gas=187.20, dd_solv=-162.97, dd_net=24.23, c=0.2
    )["sign_ill_conditioned_cancellation"] is True


# --------------------------------------------------------------------
# (d) pair undefined (single system) -> per-snap gate = N/A; r_A logged only
# --------------------------------------------------------------------
def test_d_per_snap_logs_ratio_no_flag():
    # Real Cp4 s251 per-snap sample: gas -90.65 / solv +118.31 / net +27.66.
    # r_A = 27.66 / (90.65 + 118.31) = 0.1324, logged only. The per-snap
    # gate raises NO flag (sign_invalid_gas_dominated is the retired Mode A
    # key, always None now).
    res = {
        "delta_g_kcal":      27.66,
        "delta_g_gas_kcal":  -90.65,
        "delta_g_solv_kcal": 118.31,
    }
    out = apply_pbsa_sanity_gate(res)
    assert out["endpoint_cancellation_ratio"] == pytest.approx(
        27.66 / (90.65 + 118.31), abs=1e-6
    )
    assert out["sign_invalid_gas_dominated"] is None
    assert out["sanity_gate_reason"] is None


def test_d_summaries_pair_undefined_raises():
    # A single cohort summary lacking the partner cannot define a ΔΔ pair;
    # the from-summaries wrapper raises rather than silently fabricating one.
    summary_a = {"mean_gas": -90.0, "mean_solv": 118.0, "mean_dg": 27.0}
    summary_b = {"mean_dg": -18.0}  # no gas/solv -> undefined
    with pytest.raises(ValueError):
        apply_pbsa_pair_sanity_gate_from_summaries(summary_a, summary_b)


def test_d_summaries_pair_defined_flags():
    # Two full cohort summaries (Cp4 vs WT) reproduce the canonical ΔΔ and
    # flag it. ΔΔgas = -31.83-(-219.03)=187.20, ΔΔsolv=37.85-200.82=-162.97,
    # ΔΔnet = 6.02-(-18.21)=24.23 -> r1 = 0.0692 -> flagged.
    summary_cp4 = {"mean_gas": -31.83, "mean_solv": 37.85, "mean_dg": 6.02}
    summary_wt = {"mean_gas": -219.03, "mean_solv": 200.82, "mean_dg": -18.21}
    out = apply_pbsa_pair_sanity_gate_from_summaries(summary_cp4, summary_wt)
    assert out["dd_gas_kcal"] == pytest.approx(187.20, abs=1e-2)
    assert out["dd_solv_kcal"] == pytest.approx(-162.97, abs=1e-2)
    assert out["dd_net_kcal"] == pytest.approx(24.23, abs=1e-2)
    assert out["dd_cancellation_r1"] == pytest.approx(0.0692, abs=1e-3)
    assert out["sign_ill_conditioned_cancellation"] is True


# --------------------------------------------------------------------
# (e) additive regression: existing per-snap keys unchanged
# --------------------------------------------------------------------
def test_e_per_snap_is_additive_only():
    res = {
        "snapshot":           "demo_snap",
        "delta_g_kcal":       7.2729,
        "delta_epb_kcal":     1.0720,
        "delta_enpolar_kcal": -8.9565,
        "delta_edisper_kcal": 21.9772,
        "delta_g_gas_kcal":   -6.8198,
        "delta_g_solv_kcal":  14.0927,
        "favorable":          False,
    }
    before = dict(res)
    out = apply_pbsa_sanity_gate(res)
    # original keys byte-identical
    for key, val in before.items():
        assert out[key] == val
    # exactly the additive per-snap keys
    added = set(out) - set(before)
    assert added == {
        "endpoint_cancellation_ratio",
        "sign_invalid_gas_dominated",
        "sanity_gate_reason",
    }


def test_e_per_snap_missing_decomposition_noop():
    res = {"delta_g_kcal": -5.0}
    out = apply_pbsa_sanity_gate(res)
    assert out["endpoint_cancellation_ratio"] is None
    assert out["sign_invalid_gas_dominated"] is None
    assert out["sanity_gate_reason"] is None


def test_e_pair_gate_keys_additive():
    out = apply_pbsa_pair_sanity_gate(dd_gas=10.0, dd_solv=-9.0, dd_net=1.0)
    assert set(out) == {
        "dd_gas_kcal",
        "dd_solv_kcal",
        "dd_net_kcal",
        "dd_cancellation_r1",
        "dd_cancellation_r2",
        "sanity_gate_cancel_c",
        "sign_ill_conditioned_cancellation",
        "sanity_gate_reason",
    }


def test_e_pair_zero_dd_magnitude_no_crash():
    # Both ΔΔgas and ΔΔsolv vanish -> no cancellation magnitude to assess,
    # ratios None, flag stays False (cannot judge).
    out = apply_pbsa_pair_sanity_gate(dd_gas=0.0, dd_solv=0.0, dd_net=0.0)
    assert out["dd_cancellation_r1"] is None
    assert out["dd_cancellation_r2"] is None
    assert out["sign_ill_conditioned_cancellation"] is False


# --------------------------------------------------------------------
# (f) synthetic large, well-resolved ΔΔ passes -- SYNTHETIC FIXTURE ONLY
# --------------------------------------------------------------------
def test_f_synthetic_large_effect_passes():
    """A large, well-resolved ΔΔ (net not a small residual) must pass.

    IMPORTANT — SYNTHETIC FIXTURE: the 2QKI Cp4/WT cohorts are entirely in
    the cancellation regime and contain NO real variant pair with r1 >> 0.2.
    These numbers are invented to exercise the pass-through branch; they are
    NOT measured data. The upper-side separation of c = 0.15 cannot be
    verified on 2QKI alone and should be re-calibrated on another system
    (e.g. 7TL8) when a large-effect reference exists.
    """
    # r1 = 50.0 / (40.0 + 30.0) = 0.714 >> 0.15 -> not flagged.
    out = apply_pbsa_pair_sanity_gate(dd_gas=-40.0, dd_solv=30.0, dd_net=-50.0)
    assert out["dd_cancellation_r1"] == pytest.approx(50.0 / 70.0, abs=1e-6)
    assert out["sign_ill_conditioned_cancellation"] is False
    assert out["sanity_gate_reason"] is None


# --------------------------------------------------------------------
# (g) parse DELTA G gas / DELTA G solv from a sample .dat fixture
# --------------------------------------------------------------------
_SAMPLE_DAT = """\
| Run on ...

Complex:
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
VDWAALS                  -5300.6530                0.0000              0.0000
EEL                     -44696.9901                0.0000              0.0000

Differences (Complex - Receptor - Ligand):
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
VDWAALS                    -11.8562                0.0000              0.0000
EEL                          5.0364                0.0000              0.0000
EPB                          1.0720                0.0000              0.0000
ENPOLAR                     -8.9565                0.0000              0.0000
EDISPER                     21.9772                0.0000              0.0000

DELTA G gas                 -6.8198                0.0000              0.0000
DELTA G solv                14.0927                0.0000              0.0000

DELTA TOTAL                  7.2729                0.0000              0.0000
"""


def _load_run_mmpbsa():
    """Load scripts/run_mmpbsa.py by filespec under a stable name so its
    private parser is testable without importing the package."""
    path = _REPO_ROOT / "scripts" / "run_mmpbsa.py"
    spec = _ilu.spec_from_file_location("scripts_run_mmpbsa_test", path)
    assert spec is not None and spec.loader is not None
    mod = _ilu.module_from_spec(spec)
    sys.modules["scripts_run_mmpbsa_test"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def test_g_parse_gas_solv(tmp_path):
    mod = _load_run_mmpbsa()
    dat = tmp_path / "FINAL_RESULTS_MMPBSA.dat"
    dat.write_text(_SAMPLE_DAT, encoding="utf-8")
    parsed = mod._parse_mmpbsa_output(str(dat))
    assert parsed is not None
    assert parsed["delta_g_gas_kcal"] == pytest.approx(-6.8198, abs=1e-4)
    assert parsed["delta_g_solv_kcal"] == pytest.approx(14.0927, abs=1e-4)
    assert parsed["delta_g_kcal"] == pytest.approx(7.2729, abs=1e-4)
    # Cross-check: GGAS must equal VDWAALS + EEL (gas-phase consistency).
    assert parsed["delta_g_gas_kcal"] == pytest.approx(-11.8562 + 5.0364,
                                                       abs=1e-3)


def test_g_parse_missing_gas_solv_returns_none(tmp_path):
    # A legacy .dat with no GGAS/GSOLV lines still parses TOTAL/EPB and
    # leaves the new keys None (additive, no crash).
    mod = _load_run_mmpbsa()
    legacy = "\n".join(
        ln for ln in _SAMPLE_DAT.splitlines()
        if not ln.startswith("DELTA G gas")
        and not ln.startswith("DELTA G solv")
    )
    dat = tmp_path / "FINAL_RESULTS_MMPBSA.dat"
    dat.write_text(legacy + "\n", encoding="utf-8")
    parsed = mod._parse_mmpbsa_output(str(dat))
    assert parsed is not None
    assert parsed["delta_g_gas_kcal"] is None
    assert parsed["delta_g_solv_kcal"] is None
    assert parsed["delta_g_kcal"] == pytest.approx(7.2729, abs=1e-4)
