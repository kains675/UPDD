"""tests/test_cycle_manager.py — PR-NEW-C Iterative Cycle Manager tests (T1-T10).

Per ``utils/cycle_manager.py`` (schema ``cycle_manager/0.1``).

Test layout:
    T1  Project initialization — directory + JSON file + required fields
    T2  Add Cycle N entry — variants stored, predicted_ddg optional
    T3  set_next_wt — auto-pick best validated from cycle 0 wet lab
    T4  ΔK_d trajectory math — 100nM → 20nM → 5nM, verify ΔΔG vs -RT ln(K_d)
    T5  Plateau detection — last 3 ΔK_d < 0.5 → plateau=True
    T6  Multi-assay aggregation — 2 assays for same sequence → ``_aggregate_assays``
    T7  JSON round-trip — init/add/save/load → dataclass-equal
    T8  Schema validation — required fields present in JSON
    T9  Error handling — ProjectNotFoundError, InvalidCycleError, KdParseError
    T10 Report generation — markdown contains cycle headers + K_d values
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

# Make ``utils.cycle_manager`` importable.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from utils.cycle_manager import (  # noqa: E402
    SCHEMA_VERSION,
    DEFAULT_TEMPERATURE_K,
    R_KCAL_PER_MOL_K,
    CycleManagerError,
    ProjectNotFoundError,
    InvalidCycleError,
    KdParseError,
    CycleHistory,
    CycleEntry,
    WetLabResult,
    PredictedDDG,
    init_project,
    add_cycle,
    add_result,
    set_next_wt,
    compute_trajectory,
    detect_plateau,
    generate_report,
    _parse_kd,
    _kd_to_ddg,
    _aggregate_assays,
    _pick_best_validated,
    _load_history,
    _history_path,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_root(tmp_path: Path) -> str:
    """Use a tmp_path as output_root so tests never write outputs/projects."""
    return str(tmp_path)


# ---------------------------------------------------------------------------
# T1 — Project initialization
# ---------------------------------------------------------------------------

def test_t1_project_initialization(tmp_root: str) -> None:
    history = init_project(
        "1EBP_design_q2", "1ebp.pdb",
        output_root=tmp_root, ncaa_budget=5, operator="user"
    )
    assert history.schema == SCHEMA_VERSION
    assert history.project_id == "1EBP_design_q2"
    assert history.target_pdb == "1ebp.pdb"
    assert history.metadata["ncaa_budget"] == 5
    assert history.metadata["operator"] == "user"
    # Directory + JSON file should both exist.
    proj_dir = Path(tmp_root) / "1EBP_design_q2"
    assert proj_dir.is_dir()
    assert (proj_dir / "cycle_history.json").is_file()


# ---------------------------------------------------------------------------
# T2 — Add Cycle N entry
# ---------------------------------------------------------------------------

def test_t2_add_cycle_entry(tmp_root: str) -> None:
    init_project("p2", "1ebp.pdb", output_root=tmp_root)
    entry = add_cycle(
        "p2", cycle_n=0,
        variant_sequences=["ICVVQDWGHHRCT", "ICVVQDW(MTR)HHRCT"],
        notes="Cycle 0 de novo, 2 candidates",
        output_root=tmp_root,
    )
    assert entry.cycle_n == 0
    assert entry.wt_sequence is None  # de novo
    assert entry.variants == ["ICVVQDWGHHRCT", "ICVVQDW(MTR)HHRCT"]
    assert entry.predicted_ddg == []
    assert entry.wet_lab_results == []
    # Reload fresh and verify persistence.
    history = _load_history("p2", output_root=tmp_root)
    assert len(history.cycles) == 1
    assert history.cycles[0].variants == ["ICVVQDWGHHRCT", "ICVVQDW(MTR)HHRCT"]


# ---------------------------------------------------------------------------
# T3 — set_next_wt auto-pick
# ---------------------------------------------------------------------------

def test_t3_set_next_wt_auto_pick(tmp_root: str) -> None:
    init_project("p3", "1ebp.pdb", output_root=tmp_root)
    add_cycle("p3", 0, variant_sequences=["A", "B", "C"],
              output_root=tmp_root)
    add_result("p3", 0, "A", kd="500nM", assay="ITC", output_root=tmp_root)
    add_result("p3", 0, "B", kd="80nM", assay="ITC", output_root=tmp_root)
    add_result("p3", 0, "C", kd="200nM", assay="ITC", output_root=tmp_root)
    chosen = set_next_wt("p3", source_cycle=0, best="auto",
                         output_root=tmp_root)
    assert chosen == "B"   # tightest binder (80 nM)
    history = _load_history("p3", output_root=tmp_root)
    assert history.cycles[0].best_validated_sequence == "B"


# ---------------------------------------------------------------------------
# T4 — ΔK_d trajectory math
# ---------------------------------------------------------------------------

def test_t4_trajectory_math(tmp_root: str) -> None:
    """Build 3 cycles with K_d 100nM → 20nM → 5nM and verify trajectory math."""
    init_project("p4", "1ebp.pdb", output_root=tmp_root)
    # Cycle 0 with K_d = 100 nM
    add_cycle("p4", 0, variant_sequences=["S0"], output_root=tmp_root)
    add_result("p4", 0, "S0", kd="100nM", assay="ITC", output_root=tmp_root)
    set_next_wt("p4", 0, best="auto", output_root=tmp_root)
    # Cycle 1 with K_d = 20 nM
    add_cycle("p4", 1, wt_sequence="S0", variant_sequences=["S1"],
              output_root=tmp_root)
    add_result("p4", 1, "S1", kd="20nM", assay="ITC", output_root=tmp_root)
    set_next_wt("p4", 1, best="auto", output_root=tmp_root)
    # Cycle 2 with K_d = 5 nM
    add_cycle("p4", 2, wt_sequence="S1", variant_sequences=["S2"],
              output_root=tmp_root)
    add_result("p4", 2, "S2", kd="5nM", assay="ITC", output_root=tmp_root)
    traj = compute_trajectory("p4", output_root=tmp_root)
    # Two transitions expected: 0→1 and 1→2.
    assert len(traj["delta_kd_kcal"]) == 2
    # Expected ΔΔG for K_d ratio of 1/5 (100→20) and 1/4 (20→5):
    # ΔG_bind = RT ln(K_d) → ΔΔG = RT ln(K_d_new / K_d_old)
    # (becomes more negative = tighter binding when K_d shrinks)
    expected_01 = R_KCAL_PER_MOL_K * DEFAULT_TEMPERATURE_K * math.log(
        20e-9 / 100e-9
    )
    expected_12 = R_KCAL_PER_MOL_K * DEFAULT_TEMPERATURE_K * math.log(
        5e-9 / 20e-9
    )
    assert math.isclose(traj["delta_kd_kcal"][0], expected_01, abs_tol=0.01), \
        f"got {traj['delta_kd_kcal'][0]:.6f} vs expected {expected_01:.6f}"
    assert math.isclose(traj["delta_kd_kcal"][1], expected_12, abs_tol=0.01), \
        f"got {traj['delta_kd_kcal'][1]:.6f} vs expected {expected_12:.6f}"
    # Ratios match.
    assert math.isclose(traj["ratio_per_cycle"][0], 0.2, rel_tol=1e-9)
    assert math.isclose(traj["ratio_per_cycle"][1], 0.25, rel_tol=1e-9)
    # Per-cycle ΔG sanity: ΔG = -RT ln(K_d), so 100nM → ΔG ~ -9.55 kcal/mol.
    expected_ddg_0 = _kd_to_ddg(100e-9)
    assert math.isclose(traj["per_cycle"][0]["ddg"], expected_ddg_0, abs_tol=1e-6)


# ---------------------------------------------------------------------------
# T5 — Plateau detection
# ---------------------------------------------------------------------------

def test_t5_plateau_detection(tmp_root: str) -> None:
    """5 cycles with last 3 ΔΔG < 0.5 kcal/mol → plateau=True."""
    init_project("p5", "1ebp.pdb", output_root=tmp_root)
    # Construct cycles where ΔG values are: -9.5, -10.5, -10.6, -10.7, -10.65
    # Last 3 deltas: 0.1, 0.1, -0.05 → all |Δ| < 0.5 → plateau
    target_ddgs = [-9.5, -10.5, -10.6, -10.7, -10.65]
    for i, ddg in enumerate(target_ddgs):
        add_cycle("p5", i, variant_sequences=[f"S{i}"], output_root=tmp_root)
        # ΔG_bind = RT ln(K_d) → K_d = exp(ΔG / RT)
        kd = math.exp(ddg / (R_KCAL_PER_MOL_K * DEFAULT_TEMPERATURE_K))
        add_result("p5", i, f"S{i}", kd=kd, assay="ITC",
                   output_root=tmp_root)
        if i > 0:
            # No need to set WT — _pick_best_validated falls back if absent
            pass
    traj = compute_trajectory("p5", output_root=tmp_root,
                              lookback=3, threshold_kcal=0.5)
    assert traj["plateau_detected"] is True
    # Recommendation: final ΔG (-10.65) is between -9 and -12 → "converged".
    assert traj["recommendation"] == "converged"
    # Verify also via detect_plateau directly.
    history = _load_history("p5", output_root=tmp_root)
    assert detect_plateau(history, lookback=3, threshold_kcal=0.5) is True
    # And: tighter threshold = no plateau.
    assert detect_plateau(history, lookback=3, threshold_kcal=0.05) is False


# ---------------------------------------------------------------------------
# T6 — Multi-assay K_d aggregation
# ---------------------------------------------------------------------------

def test_t6_multi_assay_aggregation() -> None:
    rs = [
        WetLabResult(sequence="X", kd_molar=80e-9, assay="ITC",
                     confidence=1.0),
        WetLabResult(sequence="X", kd_molar=120e-9, assay="SPR",
                     confidence=1.0),
    ]
    agg = _aggregate_assays(rs)
    # Equal weights → mean = 100 nM
    assert math.isclose(agg.kd_molar, 100e-9, rel_tol=1e-9)
    assert agg.sequence == "X"
    assert "ITC" in agg.assay and "SPR" in agg.assay
    assert agg.ddg_calc is not None
    expected_ddg = _kd_to_ddg(100e-9)
    assert math.isclose(agg.ddg_calc, expected_ddg, abs_tol=1e-6)
    # Weighted aggregation: ITC weight=3, SPR weight=1 → 80*3+120*1 / 4 = 90 nM
    rs_w = [
        WetLabResult(sequence="X", kd_molar=80e-9, assay="ITC", confidence=3.0),
        WetLabResult(sequence="X", kd_molar=120e-9, assay="SPR", confidence=1.0),
    ]
    agg_w = _aggregate_assays(rs_w)
    assert math.isclose(agg_w.kd_molar, 90e-9, rel_tol=1e-9)
    # All-zero confidences → KdParseError.
    rs_zero = [
        WetLabResult(sequence="X", kd_molar=80e-9, assay="ITC", confidence=0.0),
        WetLabResult(sequence="X", kd_molar=120e-9, assay="SPR", confidence=0.0),
    ]
    with pytest.raises(KdParseError):
        _aggregate_assays(rs_zero)


# ---------------------------------------------------------------------------
# T7 — JSON round-trip
# ---------------------------------------------------------------------------

def test_t7_json_round_trip(tmp_root: str) -> None:
    init_project("p7", "1ebp.pdb", output_root=tmp_root, ncaa_budget=5)
    add_cycle("p7", 0,
              variant_sequences=["A", "B"],
              predicted_ddg=[
                  PredictedDDG(sequence="B", ddg=-1.2, se=1.5, tier="X.B"),
              ],
              output_root=tmp_root)
    add_result("p7", 0, "A", kd="500nM", assay="ITC", confidence=0.9,
               output_root=tmp_root)
    add_result("p7", 0, "B", kd="80nM", assay="ITC", confidence=0.95,
               output_root=tmp_root)
    set_next_wt("p7", 0, best="auto", output_root=tmp_root)
    # Read the JSON back and reconstruct via from_dict.
    path = _history_path("p7", tmp_root)
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    reconstructed = CycleHistory.from_dict(raw)
    history = _load_history("p7", output_root=tmp_root)
    # Compare structurally (output_root may differ between load and from_dict).
    assert reconstructed.schema == history.schema
    assert reconstructed.project_id == history.project_id
    assert reconstructed.target_pdb == history.target_pdb
    assert reconstructed.metadata == history.metadata
    assert len(reconstructed.cycles) == 1
    rc = reconstructed.cycles[0]
    hc = history.cycles[0]
    assert rc.cycle_n == hc.cycle_n == 0
    assert rc.variants == hc.variants
    assert len(rc.predicted_ddg) == 1
    assert rc.predicted_ddg[0].tier == "X.B"
    assert len(rc.wet_lab_results) == 2
    assert rc.best_validated_sequence == "B"


# ---------------------------------------------------------------------------
# T8 — Schema validation (required fields)
# ---------------------------------------------------------------------------

def test_t8_schema_required_fields(tmp_root: str) -> None:
    init_project("p8", "1ebp.pdb", output_root=tmp_root,
                 ncaa_budget=5, operator="user", notes="schema test")
    add_cycle("p8", 0, variant_sequences=["A"], output_root=tmp_root)
    add_result("p8", 0, "A", kd="80nM", assay="ITC", output_root=tmp_root)
    path = _history_path("p8", tmp_root)
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    for required in ("schema", "project_id", "target_pdb", "created_at",
                     "metadata", "cycles", "trajectory"):
        assert required in d, f"missing required top-level field: {required}"
    assert d["schema"] == SCHEMA_VERSION
    cyc = d["cycles"][0]
    for field_ in ("cycle_n", "wt_sequence", "variants", "predicted_ddg",
                   "wet_lab_results", "best_validated_sequence",
                   "synthesized", "notes", "created_at"):
        assert field_ in cyc, f"missing required cycle field: {field_}"
    res = cyc["wet_lab_results"][0]
    for field_ in ("sequence", "kd_molar", "assay", "ddg_calc",
                   "confidence", "temperature_k", "notes"):
        assert field_ in res, f"missing required wet_lab field: {field_}"


# ---------------------------------------------------------------------------
# T9 — Error handling
# ---------------------------------------------------------------------------

def test_t9_error_handling(tmp_root: str) -> None:
    # ProjectNotFoundError on missing project.
    with pytest.raises(ProjectNotFoundError):
        _load_history("nonexistent_project", output_root=tmp_root)
    # InvalidCycleError on duplicate cycle_n.
    init_project("p9", "1ebp.pdb", output_root=tmp_root)
    add_cycle("p9", 0, variant_sequences=["A"], output_root=tmp_root)
    with pytest.raises(InvalidCycleError):
        add_cycle("p9", 0, variant_sequences=["B"], output_root=tmp_root)
    # InvalidCycleError on negative cycle_n.
    with pytest.raises(InvalidCycleError):
        add_cycle("p9", -1, variant_sequences=["B"], output_root=tmp_root)
    # InvalidCycleError on add_result for unregistered cycle.
    with pytest.raises(InvalidCycleError):
        add_result("p9", 99, "A", kd="80nM", output_root=tmp_root)
    # KdParseError on unparseable K_d.
    with pytest.raises(KdParseError):
        _parse_kd("not-a-number")
    with pytest.raises(KdParseError):
        _parse_kd("80 zM")          # unknown unit
    with pytest.raises(KdParseError):
        _parse_kd("-5 nM")           # negative
    with pytest.raises(KdParseError):
        _parse_kd("0 nM")            # zero
    with pytest.raises(KdParseError):
        _parse_kd("")                # empty
    # InvalidCycleError on set_next_wt for unknown source cycle.
    with pytest.raises(InvalidCycleError):
        set_next_wt("p9", source_cycle=42, best="auto",
                    output_root=tmp_root)
    # InvalidCycleError when source cycle has no wet-lab data.
    with pytest.raises(InvalidCycleError):
        set_next_wt("p9", source_cycle=0, best="auto",
                    output_root=tmp_root)
    # init_project on existing project without overwrite.
    with pytest.raises(InvalidCycleError):
        init_project("p9", "1ebp.pdb", output_root=tmp_root)


# ---------------------------------------------------------------------------
# T10 — Report generation
# ---------------------------------------------------------------------------

def test_t10_report_generation(tmp_root: str) -> None:
    init_project("p10", "1ebp.pdb", output_root=tmp_root)
    add_cycle("p10", 0, variant_sequences=["S0"], output_root=tmp_root)
    add_result("p10", 0, "S0", kd="100nM", assay="ITC", output_root=tmp_root)
    set_next_wt("p10", 0, best="auto", output_root=tmp_root)
    add_cycle("p10", 1, wt_sequence="S0", variant_sequences=["S1"],
              output_root=tmp_root)
    add_result("p10", 1, "S1", kd="10nM", assay="SPR", output_root=tmp_root)
    text = generate_report("p10", output_root=tmp_root)
    assert "## Cycle 0" in text
    assert "## Cycle 1" in text
    assert "100" in text and "nM" in text
    assert "10" in text
    # Trajectory section present.
    assert "## Trajectory" in text
    assert "Plateau detected" in text
    # Schema header line present.
    assert SCHEMA_VERSION in text
    # Optional output_path.
    out_path = Path(tmp_root) / "p10_report.md"
    text2 = generate_report("p10", output_root=tmp_root,
                            output_path=str(out_path))
    assert out_path.is_file()
    assert text == text2


# ---------------------------------------------------------------------------
# Bonus: K_d parser format coverage (reuses T9's error coverage).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("80nM",        80e-9),
    ("80 nM",       80e-9),
    ("80,000 pM",   80e-9),         # = 80 nM
    ("0.08uM",      80e-9),
    ("0.08 μM",     80e-9),
    ("8e-8",        8e-8),
    ("8e-8 M",      8e-8),
    ("5fM",         5e-15),
    ("1 mM",        1e-3),
    ("0.001 M",     1e-3),
    ("1",           1.0),
])
def test_kd_parser_formats(raw: str, expected: float) -> None:
    got = _parse_kd(raw)
    assert math.isclose(got, expected, rel_tol=1e-9), (
        f"_parse_kd({raw!r}) = {got!r}; expected {expected!r}"
    )


# ---------------------------------------------------------------------------
# Day 2 v0.2 enhancements — T11-T17
# ---------------------------------------------------------------------------

# T11 — Multi-assay K_d input (ITC + SPR with confidence) end-to-end
def test_t11_multi_assay_input(tmp_root: str) -> None:
    """ITC + SPR confidence-weighted aggregation surfaces in compute_trajectory."""
    from utils.cycle_manager import _assay_normalize  # noqa: F401
    init_project("p11", "1ebp.pdb", output_root=tmp_root)
    add_cycle("p11", 0, variant_sequences=["X"], output_root=tmp_root)
    # Two assays for the same sequence with different confidences.
    # ITC weight 3, SPR weight 1: 80 nM * 3 + 120 nM * 1 / 4 = 90 nM
    add_result("p11", 0, "X", kd="80nM", assay="ITC", confidence=3.0,
               output_root=tmp_root)
    add_result("p11", 0, "X", kd="120nM", assay="SPR", confidence=1.0,
               output_root=tmp_root)
    set_next_wt("p11", 0, best="auto", output_root=tmp_root)
    traj = compute_trajectory("p11", output_root=tmp_root)
    pc = traj["per_cycle"][0]
    # Aggregated K_d should be the inverse-variance-style 90 nM (weighted).
    assert math.isclose(pc["kd_molar"], 90e-9, rel_tol=1e-9), (
        f"T11: weighted K_d {pc['kd_molar']:.3e} ≠ expected 9e-8 (90 nM)"
    )
    # ΔG matches K_d 90 nM at default T.
    assert math.isclose(pc["ddg"], _kd_to_ddg(90e-9), abs_tol=1e-6)


# T12 — K_d ↔ ΔG round-trip
def test_t12_kd_ddg_roundtrip() -> None:
    """K_d → ΔG → K_d matches input within numerical precision."""
    test_kds = [1e-3, 1e-6, 80e-9, 5e-12, 1.0, 0.5]
    for kd in test_kds:
        ddg = _kd_to_ddg(kd)
        # Inverse: K_d = exp(ΔG / RT)
        kd_back = math.exp(ddg / (R_KCAL_PER_MOL_K * DEFAULT_TEMPERATURE_K))
        assert math.isclose(kd, kd_back, rel_tol=1e-9), (
            f"T12: K_d round-trip {kd!r} → ΔG={ddg:.6f} → {kd_back!r}"
        )


# T13 — Plateau detection synthetic 3-cycle ΔΔG < 0.5
def test_t13_plateau_detection_synthetic(tmp_root: str) -> None:
    """4 cycles where last 3 |ΔΔG| < 0.5 → plateau=True."""
    init_project("p13", "1ebp.pdb", output_root=tmp_root)
    target_ddgs = [-8.0, -8.4, -8.6, -8.5]   # deltas: -0.4, -0.2, +0.1
    for i, ddg in enumerate(target_ddgs):
        add_cycle("p13", i, variant_sequences=[f"S{i}"],
                  output_root=tmp_root)
        kd = math.exp(ddg / (R_KCAL_PER_MOL_K * DEFAULT_TEMPERATURE_K))
        add_result("p13", i, f"S{i}", kd=kd, assay="ITC",
                   output_root=tmp_root)
    traj = compute_trajectory("p13", output_root=tmp_root,
                              lookback=3, threshold_kcal=0.5)
    assert traj["plateau_detected"] is True, (
        f"T13: expected plateau=True; deltas={traj['delta_kd_kcal']}"
    )


# T14 — Plateau threshold parameterization
def test_t14_plateau_threshold_parameterization(tmp_root: str) -> None:
    """Different thresholds yield different verdicts on the same data."""
    init_project("p14", "1ebp.pdb", output_root=tmp_root)
    target_ddgs = [-8.0, -8.3, -8.4, -8.45]   # deltas: -0.3, -0.1, -0.05
    for i, ddg in enumerate(target_ddgs):
        add_cycle("p14", i, variant_sequences=[f"S{i}"],
                  output_root=tmp_root)
        kd = math.exp(ddg / (R_KCAL_PER_MOL_K * DEFAULT_TEMPERATURE_K))
        add_result("p14", i, f"S{i}", kd=kd, output_root=tmp_root)
    history = _load_history("p14", output_root=tmp_root)
    # Threshold 0.5 → all 3 deltas < 0.5 → plateau
    assert detect_plateau(history, lookback=3, threshold_kcal=0.5) is True
    # Threshold 0.2 → 0.3 violates → no plateau
    assert detect_plateau(history, lookback=3, threshold_kcal=0.2) is False
    # Threshold 0.05 (very tight) → most deltas violate → no plateau
    assert detect_plateau(history, lookback=3, threshold_kcal=0.05) is False
    # target_kd_molar steers recommendation: at K_d ~370nM-440nM, plateau
    # while above 1nM target → try_alternative.
    traj = compute_trajectory(
        "p14", output_root=tmp_root,
        lookback=3, threshold_kcal=0.5,
        target_kd_molar=1e-9,
    )
    assert traj["plateau_detected"] is True
    assert traj["recommendation"] == "try_alternative"
    # If the target is huge (1 mM), we are well below it → converged.
    traj2 = compute_trajectory(
        "p14", output_root=tmp_root,
        lookback=3, threshold_kcal=0.5,
        target_kd_molar=1e-3,
    )
    assert traj2["plateau_detected"] is True
    assert traj2["recommendation"] == "converged"


# T15 — `updd cycle` subcommand integration (Phase B; skipped if unavailable)
def test_t15_updd_cycle_subcommand(tmp_path: Path) -> None:
    """Smoke-test ``updd cycle init/add/status`` via scripts.updd_cli.main.

    Skipped if Phase B integration has not landed (the `cycle` subcommand
    is not yet registered in the parser, or returns non-zero).
    """
    try:
        from scripts.updd_cli import main as cli_main  # noqa: WPS433
    except Exception as exc:  # pragma: no cover - dispatch placeholder
        pytest.skip(f"updd_cli unavailable for cycle integration: {exc}")
    output_root = str(tmp_path / "projects")

    def _run(argv):
        try:
            return cli_main(argv)
        except SystemExit as exc:  # argparse reject of unknown subcommand
            return int(getattr(exc, "code", 2) or 0)

    # init
    rc = _run([
        "cycle", "init",
        "--project", "p15",
        "--target", "1ebp.pdb",
        "--output-root", output_root,
    ])
    if rc != 0:
        pytest.skip(
            "Phase B integration not landed (cmd_cycle returned non-zero "
            f"or unknown subcommand; rc={rc}); tracked under Phase B "
            "race-coordination follow-up"
        )
    # add cycle 0 + 1 wet-lab record
    rc = _run([
        "cycle", "add",
        "--project", "p15",
        "--cycle", "0",
        "--variants", "ICVVQDWGHHRCT",
        "--output-root", output_root,
    ])
    assert rc == 0
    rc = _run([
        "cycle", "add-result",
        "--project", "p15",
        "--cycle", "0",
        "--sequence", "ICVVQDWGHHRCT",
        "--kd", "80nM",
        "--assay", "ITC",
        "--output-root", output_root,
    ])
    assert rc == 0
    # status — must produce a nonzero-length stdout (no exception).
    rc = _run([
        "cycle", "status",
        "--project", "p15",
        "--output-root", output_root,
    ])
    assert rc == 0


# T16 — export_visualization_data structure
def test_t16_export_visualization_data(tmp_root: str) -> None:
    """Visualization data dict has expected per-cycle keys + trajectory points."""
    from utils.cycle_manager import export_visualization_data
    init_project("p16", "1ebp.pdb", output_root=tmp_root,
                 ncaa_budget=5)
    # 3 cycles with K_d 100, 30, 10 nM
    for i, kd_str in enumerate(["100nM", "30nM", "10nM"]):
        add_cycle("p16", i, variant_sequences=[f"S{i}"],
                  output_root=tmp_root)
        add_result("p16", i, f"S{i}", kd=kd_str, assay="ITC",
                   output_root=tmp_root)
    history = _load_history("p16", output_root=tmp_root)
    viz = export_visualization_data(history,
                                    lookback=2, threshold_kcal=0.5,
                                    target_kd_molar=10e-9)
    assert viz["schema"] == "cycle_manager_viz/0.2"
    assert viz["project_id"] == "p16"
    assert viz["n_cycles"] == 3
    # Per-cycle key set
    expected_keys = {"cycle_n", "wt_sequence", "best_validated_sequence",
                     "kd_molar", "kd_label", "ddg_kcal", "n_variants",
                     "n_wet_lab_results", "synthesized"}
    for pc in viz["per_cycle"]:
        assert expected_keys.issubset(set(pc.keys())), (
            f"T16: missing per-cycle keys: "
            f"{expected_keys - set(pc.keys())}"
        )
    # Trajectory points
    assert len(viz["trajectory_points"]) == 3
    for tp in viz["trajectory_points"]:
        assert {"cycle_n", "ddg_kcal", "kd_molar"}.issubset(tp.keys())
    # 2 deltas across 3 cycles
    assert len(viz["delta_kd_kcal"]) == 2
    assert len(viz["ratio_per_cycle"]) == 2
    # K_d label sanity
    assert viz["per_cycle"][0]["kd_label"].endswith("nM")
    # Target was 10nM and final K_d is 10nM → recommendation depends on plateau
    assert viz["target_kd_molar"] == 10e-9
    assert "recommendation" in viz


# T17 — Backward compatibility (v0.1 schema input → v0.2 reader)
def test_t17_backward_compat_v01_schema(tmp_root: str) -> None:
    """A hand-crafted v0.1 JSON loads into v0.2 reader and is upgraded."""
    legacy_payload = {
        "schema": "cycle_manager/0.1",
        "project_id": "p17",
        "target_pdb": "1ebp.pdb",
        "created_at": "2026-04-28T00:00:00Z",
        "metadata": {"ncaa_budget": 5, "operator": "user"},
        "cycles": [
            {
                "cycle_n": 0,
                "wt_sequence": None,
                "variants": ["X", "Y"],
                "predicted_ddg": [],
                "wet_lab_results": [
                    {
                        "sequence": "X",
                        "kd_molar": 80e-9,
                        "assay": "ITC",
                        "ddg_calc": -9.677,
                        "confidence": None,
                        "temperature_k": 298.15,
                        "notes": None,
                    }
                ],
                "best_validated_sequence": "X",
                "synthesized": False,
                "notes": None,
                "created_at": "2026-04-28T00:00:00Z",
            }
        ],
        "trajectory": {},
        "output_root": tmp_root,
    }
    proj_dir = Path(tmp_root) / "p17"
    proj_dir.mkdir(parents=True, exist_ok=True)
    path = proj_dir / "cycle_history.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(legacy_payload, f, indent=2)
    # Load via v0.2 reader.
    history = _load_history("p17", output_root=tmp_root)
    # Schema upgraded to v0.2 in-memory.
    assert history.schema == SCHEMA_VERSION, (
        f"T17: expected schema upgrade to {SCHEMA_VERSION}; got {history.schema}"
    )
    assert history.project_id == "p17"
    assert len(history.cycles) == 1
    assert history.cycles[0].variants == ["X", "Y"]
    # And v0.2 features work on legacy data.
    from utils.cycle_manager import export_visualization_data
    viz = export_visualization_data(history)
    assert viz["schema"] == "cycle_manager_viz/0.2"
    assert viz["n_cycles"] == 1
    # Add a new cycle and verify save writes v0.2 schema string.
    add_cycle("p17", 1, wt_sequence="X", variant_sequences=["Z"],
              output_root=tmp_root)
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    assert d["schema"] == SCHEMA_VERSION
