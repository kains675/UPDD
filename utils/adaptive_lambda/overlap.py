"""Adaptive λ-scheduling DIAGNOSE module — per-adjacent-pair overlap.

Reads an existing per-direction ATM async_re pilot's ``.out`` files (GPU-FREE —
``.out`` text only, no OpenMM, no GPU), groups the perturbation-energy (col 9)
by state (col 0), and computes the per-adjacent-pair Bhattacharyya coefficient
(the decision PROXY) plus, when pymbar is importable, the MBAR overlap
matrix nearest-neighbor off-diagonal (the PRIMARY gate when available).

``.out`` format (ommreplica.save_out, 11 whitespace-separated floats per row):
    col0=stateid col1=temperature col2=direction col3=lambda1 col4=lambda2
    col5=alpha col6=u0 col7=w0 col8=potE col9=pertE col10=trash

In async replica exchange each ``r{i}/`` is a WALKER, not a state; its rows carry
whichever state the walker currently occupies (col0). So pertE-by-state is built
by aggregating rows ACROSS ALL ``r*`` walker dirs and grouping by col0
(re-aggregated pertE per state across ALL dplus replicas).

Warmup: the first ``warmup_cycles`` rows of EACH walker are discarded (in
async_re one row == one cycle per walker), using the G1 setting
(warmup 5 cycles discarded), which reproduces dplus 8→9 BC≈0.131 and
dminus 9→10 BC≈0.364.

Bhattacharyya gate metric (Gaussian/parametric form — robust at N~20, no
histogram binning artifact):
    BC = exp(-[0.25·ln(0.25·(v1/v2 + v2/v1 + 2)) + 0.25·(m1-m2)²/(v1+v2)])
Histogram-intersection is BANNED as a gate (N~20 false-positive-prone per the
verdict); it is provided ONLY as an explicitly-labeled advisory diagnostic.

Python 3.8+; numpy hard dep; pymbar optional (auto).
"""

import glob
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# ``.out`` column indices (ommreplica.save_out / uwham postprocess convention).
_COL_STATEID = 0
_COL_PERTE = 9
_N_OUT_COLS = 11


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------
@dataclass
class PairOverlap:
    """Overlap diagnostic for one adjacent state pair (i, i+1) within a run."""

    i: int
    j: int
    bc: float
    n_i: int
    n_j: int
    mean_i: float
    mean_j: float
    sd_i: float
    sd_j: float
    # MBAR nearest-neighbor off-diagonal O_{i,j} when available (else None).
    mbar_o: Optional[float] = None
    # Advisory-only histogram intersection (NEVER a gate).
    hist_intersection_advisory: Optional[float] = None

    def as_dict(self) -> Dict[str, Any]:
        return {
            "i": self.i,
            "j": self.j,
            "bc": self.bc,
            "n_i": self.n_i,
            "n_j": self.n_j,
            "mean_i": self.mean_i,
            "mean_j": self.mean_j,
            "sd_i": self.sd_i,
            "sd_j": self.sd_j,
            "mbar_o": self.mbar_o,
            "hist_intersection_advisory": self.hist_intersection_advisory,
        }


@dataclass
class OverlapReport:
    """Per-direction overlap report for a Track B per-direction pilot leg."""

    run_dir: str
    jobname: str
    warmup_cycles: int
    # direction tag ("dplus"/"dminus") -> list[PairOverlap] (adjacent pairs)
    per_direction: Dict[str, List[PairOverlap]] = field(default_factory=dict)
    # direction tag -> optional MBAR overlap matrix (None if pymbar absent)
    mbar_matrices: Dict[str, Optional[np.ndarray]] = field(default_factory=dict)
    # direction tag -> {state: n_samples}
    sample_counts: Dict[str, Dict[int, int]] = field(default_factory=dict)
    pymbar_available: bool = False

    def worst_bc(self) -> Optional[Tuple[str, PairOverlap]]:
        """Return (direction, PairOverlap) of the globally worst-overlapping
        adjacent pair (min BC across both directions), or None if empty."""
        worst: Optional[Tuple[str, PairOverlap]] = None
        for direction, pairs in self.per_direction.items():
            for p in pairs:
                if worst is None or p.bc < worst[1].bc:
                    worst = (direction, p)
        return worst

    def as_dict(self) -> Dict[str, Any]:
        return {
            "run_dir": self.run_dir,
            "jobname": self.jobname,
            "warmup_cycles": self.warmup_cycles,
            "pymbar_available": self.pymbar_available,
            "per_direction": {
                d: [p.as_dict() for p in pairs]
                for d, pairs in self.per_direction.items()
            },
            "sample_counts": {
                d: {int(k): int(v) for k, v in counts.items()}
                for d, counts in self.sample_counts.items()
            },
        }


# ---------------------------------------------------------------------------
# Extraction — pertE by state from per-direction .out files (GPU-free)
# ---------------------------------------------------------------------------
def extract_pertE_by_state(
    run_dir: str,
    jobname: str,
    direction: str,
    warmup_cycles: int = 5,
) -> Dict[int, np.ndarray]:
    """Read ``run_dir/<direction>/r*/{jobname}_{direction}.out`` and return a
    ``{state_id: pertE_array}`` dict.

    For each walker dir ``r*`` the first ``warmup_cycles`` rows are discarded
    (one row == one async_re cycle per walker), then rows are grouped by their
    col-0 stateid; col-9 is the perturbation energy. Aggregates across ALL
    walker dirs (async_re: a walker's state changes per cycle, so pertE-by-state
    requires the union over walkers — verdict re-derivation).

    ``direction`` is the tag ("dplus" / "dminus"); the per-direction subdir is
    ``run_dir/<direction>`` and the file basename is ``{jobname}_{direction}.out``.
    Returns an empty dict if no walker .out files are found (caller decides
    whether that is a refuse condition).
    """
    if direction not in ("dplus", "dminus"):
        raise ValueError(
            f"direction must be 'dplus' or 'dminus', got {direction!r}"
        )
    if warmup_cycles < 0:
        raise ValueError(f"warmup_cycles must be >= 0, got {warmup_cycles}")

    subdir = os.path.join(run_dir, direction)
    out_basename = f"{jobname}_{direction}.out"
    walker_dirs = sorted(
        glob.glob(os.path.join(subdir, "r*")),
        key=_natural_r_key,
    )

    by_state: Dict[int, List[float]] = {}
    for wdir in walker_dirs:
        if not os.path.isdir(wdir):
            continue
        out_path = os.path.join(wdir, out_basename)
        if not os.path.isfile(out_path):
            continue
        rows = _read_out_rows(out_path)
        # Discard the first warmup_cycles rows of THIS walker.
        rows = rows[warmup_cycles:]
        for sid, pe in rows:
            by_state.setdefault(sid, []).append(pe)

    return {sid: np.asarray(vals, dtype=float)
            for sid, vals in by_state.items()}


def _natural_r_key(path: str) -> Tuple[int, str]:
    """Sort ``r0, r1, r2, ..., r10`` numerically (not lexicographically)."""
    base = os.path.basename(path)
    if base.startswith("r") and base[1:].isdigit():
        return (int(base[1:]), base)
    return (10 ** 9, base)


def _read_out_rows(out_path: str) -> List[Tuple[int, float]]:
    """Read an ATM ``.out`` file → list of (stateid, pertE). Non-numeric / short
    lines are skipped defensively (abfe_production .out has no header, but this
    keeps the reader robust to blank trailing lines)."""
    rows: List[Tuple[int, float]] = []
    with open(out_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < _N_OUT_COLS:
                continue
            try:
                sid = int(float(parts[_COL_STATEID]))
                pe = float(parts[_COL_PERTE])
            except (ValueError, IndexError):
                continue
            rows.append((sid, pe))
    return rows


# ---------------------------------------------------------------------------
# Bhattacharyya coefficient (Gaussian) — the decision proxy / gate
# ---------------------------------------------------------------------------
def bhattacharyya_coefficient(a: np.ndarray, b: np.ndarray) -> float:
    """Gaussian (parametric) Bhattacharyya coefficient between two pertE samples.

        BC = exp(-[0.25·ln(0.25·(v1/v2 + v2/v1 + 2)) + 0.25·(m1-m2)²/(v1+v2)])

    where (m1, v1) and (m2, v2) are the sample mean / variance of ``a`` and
    ``b``. This is the closed-form BC for two univariate Gaussians (the
    Bhattacharyya distance D_B, BC = exp(-D_B)). Robust at N~20 (no histogram
    binning). Variance is the sample variance (ddof=1) to match the manual
    re-derivation; degenerate / zero variances are floored to a tiny epsilon so
    BC stays finite.

    Returns BC in (0, 1]. Raises ValueError if either sample has < 2 points
    (variance undefined).
    """
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if a.size < 2 or b.size < 2:
        raise ValueError(
            f"bhattacharyya_coefficient needs >= 2 samples each; got "
            f"n_a={a.size} n_b={b.size}"
        )
    m1, m2 = float(a.mean()), float(b.mean())
    v1, v2 = float(a.var(ddof=1)), float(b.var(ddof=1))
    eps = 1e-12
    if v1 <= 0.0:
        v1 = eps
    if v2 <= 0.0:
        v2 = eps
    term_var = 0.25 * np.log(0.25 * (v1 / v2 + v2 / v1 + 2.0))
    term_mean = 0.25 * (m1 - m2) ** 2 / (v1 + v2)
    return float(np.exp(-(term_var + term_mean)))


def histogram_intersection_advisory(
    a: np.ndarray,
    b: np.ndarray,
    bins: int = 23,
) -> float:
    """ADVISORY-ONLY histogram intersection of two pertE samples.

    ⚠️ BANNED as a gate (verdict: N~20 false-positive-prone — it overcalled 4
    pairs vs Bhattacharyya's 1 in the densified38 pilot). Provided ONLY as an
    explicitly-labeled advisory diagnostic for human inspection; NO judge / gate
    code may consume this value. The default 23 bins matches the verdict's
    histogram diagnostic.

    Returns the shared-area fraction in [0, 1] over a common binning of the
    pooled range.
    """
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if a.size == 0 or b.size == 0:
        return 0.0
    lo = float(min(a.min(), b.min()))
    hi = float(max(a.max(), b.max()))
    if hi <= lo:
        return 1.0
    edges = np.linspace(lo, hi, bins + 1)
    ha, _ = np.histogram(a, bins=edges, density=False)
    hb, _ = np.histogram(b, bins=edges, density=False)
    pa = ha / ha.sum() if ha.sum() > 0 else ha
    pb = hb / hb.sum() if hb.sum() > 0 else hb
    return float(np.minimum(pa, pb).sum())


# ---------------------------------------------------------------------------
# Adjacent overlaps from a pertE-by-state dict
# ---------------------------------------------------------------------------
def adjacent_overlaps(
    pertE_by_state: Dict[int, np.ndarray],
    include_advisory: bool = False,
) -> List[PairOverlap]:
    """Compute adjacent-pair Bhattacharyya overlaps for a single direction.

    States are taken in sorted order; for each contiguous pair (i, i+1) where
    BOTH states are present with >= 2 samples a ``PairOverlap`` is emitted.
    Pairs where either state has < 2 samples are SKIPPED (BC undefined) — the
    caller can detect a gap by the missing pair index.

    ``include_advisory=True`` additionally fills the (banned-as-gate) histogram
    intersection field for human inspection only.
    """
    states = sorted(pertE_by_state.keys())
    pairs: List[PairOverlap] = []
    for k in range(len(states) - 1):
        i, j = states[k], states[k + 1]
        if j != i + 1:
            # Non-contiguous (a state had no samples at all) — skip; the gap is
            # visible to the caller via the missing pair.
            continue
        a = pertE_by_state[i]
        b = pertE_by_state[j]
        if a.size < 2 or b.size < 2:
            continue
        bc = bhattacharyya_coefficient(a, b)
        adv = (
            histogram_intersection_advisory(a, b)
            if include_advisory else None
        )
        pairs.append(PairOverlap(
            i=i, j=j, bc=bc,
            n_i=int(a.size), n_j=int(b.size),
            mean_i=float(a.mean()), mean_j=float(b.mean()),
            sd_i=float(a.std(ddof=1)), sd_j=float(b.std(ddof=1)),
            hist_intersection_advisory=adv,
        ))
    return pairs


# ---------------------------------------------------------------------------
# MBAR overlap matrix (optional — pymbar auto-detected)
# ---------------------------------------------------------------------------
def mbar_overlap_matrix(
    pertE_by_state: Dict[int, np.ndarray],
    schedule: Optional[Dict[str, Any]] = None,
    temperature_K: float = 300.0,
) -> Optional[np.ndarray]:
    """Compute the MBAR overlap matrix from per-state pertE samples IF pymbar is
    importable; otherwise return None (BC remains the gate — the MVP must NOT
    hard-depend on pymbar).

    The reduced-potential matrix u_kn[k, n] = β · λ-bias(sample n evaluated at
    state k) is approximated for the MVP from the perturbation energy and the
    per-state ``lambdas`` coupling (linear: u = λ·pertE). This is the SAME
    linear-region surrogate used here (the soft-core bias is more
    involved; in the MVP MBAR is ADVISORY and Bhattacharyya is the gate, so a
    linear-λ surrogate is acceptable). When ``schedule`` is None or lacks a
    ``lambdas`` array, a uniform λ ramp over the present states is assumed.

    Returns the (n_states × n_states) overlap matrix, or None if pymbar is
    unavailable or the MBAR solve fails (the failure is swallowed — MBAR is
    optional; the caller falls back to BC).
    """
    try:
        from pymbar import MBAR
    except Exception:
        return None

    states = sorted(pertE_by_state.keys())
    n_states = len(states)
    if n_states < 2:
        return None
    samples = [pertE_by_state[s] for s in states]
    if any(s.size < 1 for s in samples):
        return None

    # β in (kcal/mol)^-1 — pertE is in the .out's native energy unit; for an
    # overlap-matrix diagnostic the absolute scale only rescales β, which
    # cancels in the overlap structure. We use kcal/mol·K convention.
    kB_kcal = 0.0019872041  # Boltzmann constant in kcal/(mol·K)
    beta = 1.0 / (kB_kcal * temperature_K)

    if schedule is not None and "lambdas" in schedule:
        lam_full = [float(x) for x in schedule["lambdas"]]
        # Map by state index when lengths agree; else uniform fallback.
        if len(lam_full) >= max(states) + 1:
            lambdas = [lam_full[s] for s in states]
        else:
            lambdas = list(np.linspace(0.0, 1.0, n_states))
    else:
        lambdas = list(np.linspace(0.0, 1.0, n_states))

    N_k = np.array([s.size for s in samples], dtype=int)
    n_total = int(N_k.sum())
    # Concatenate all samples; u_kn[k, n] = β · λ_k · pertE_n (linear surrogate).
    all_pertE = np.concatenate(samples) if n_total > 0 else np.array([])
    u_kn = np.zeros((n_states, n_total), dtype=float)
    for k in range(n_states):
        u_kn[k, :] = beta * lambdas[k] * all_pertE

    try:
        mbar = MBAR(u_kn, N_k)
        result = mbar.compute_overlap()
        # pymbar 4.x returns a dict with key "matrix"; 3.x returns a tuple.
        if isinstance(result, dict):
            matrix = np.asarray(result["matrix"], dtype=float)
        else:
            # (eigenvalues, eigenvectors, matrix) tuple form
            matrix = np.asarray(result[-1], dtype=float)
        return matrix
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Top-level estimate — per-direction adjacent BC + optional MBAR (GPU-free)
# ---------------------------------------------------------------------------
def estimate_overlap(
    run_dir: str,
    jobname: str,
    schedule: Optional[Dict[str, Any]] = None,
    use_mbar: str = "auto",
    warmup_cycles: int = 5,
    directions: Tuple[str, ...] = ("dplus", "dminus"),
    include_advisory: bool = False,
    temperature_K: float = 300.0,
) -> OverlapReport:
    """Build an ``OverlapReport`` for a per-direction pilot leg (GPU-FREE).

    Reads ONLY the per-direction ``.out`` files under
    ``run_dir/<direction>/r*/{jobname}_{direction}.out`` — no OpenMM import, no
    GPU, no launch. For each present direction it computes the adjacent-pair
    Bhattacharyya overlaps and (when ``use_mbar`` allows AND pymbar is importable)
    the MBAR overlap matrix.

    ``use_mbar``:
      * ``"auto"`` (default) — compute MBAR if pymbar importable, else skip.
      * ``"never"`` — never attempt MBAR (BC only).
      * ``"require"`` — attempt MBAR; if pymbar absent, raise RuntimeError
        (for callers that explicitly want the primary gate).

    BC always populated (the gate proxy); MBAR populated only when available.
    """
    if use_mbar not in ("auto", "never", "require"):
        raise ValueError(
            f"use_mbar must be 'auto'/'never'/'require', got {use_mbar!r}"
        )

    pymbar_ok = _pymbar_importable()
    if use_mbar == "require" and not pymbar_ok:
        raise RuntimeError(
            "use_mbar='require' but pymbar is not importable in this env"
        )

    report = OverlapReport(
        run_dir=run_dir,
        jobname=jobname,
        warmup_cycles=warmup_cycles,
        pymbar_available=pymbar_ok,
    )

    for direction in directions:
        by_state = extract_pertE_by_state(
            run_dir, jobname, direction, warmup_cycles=warmup_cycles,
        )
        if not by_state:
            continue
        pairs = adjacent_overlaps(by_state, include_advisory=include_advisory)

        matrix: Optional[np.ndarray] = None
        if use_mbar in ("auto", "require") and pymbar_ok:
            matrix = mbar_overlap_matrix(
                by_state, schedule=schedule, temperature_K=temperature_K,
            )
            if matrix is not None:
                _attach_mbar_to_pairs(pairs, by_state, matrix)

        report.per_direction[direction] = pairs
        report.mbar_matrices[direction] = matrix
        report.sample_counts[direction] = {
            int(s): int(arr.size) for s, arr in by_state.items()
        }

    return report


def _pymbar_importable() -> bool:
    try:
        import pymbar  # noqa: F401
        return True
    except Exception:
        return False


def _attach_mbar_to_pairs(
    pairs: List[PairOverlap],
    by_state: Dict[int, np.ndarray],
    matrix: np.ndarray,
) -> None:
    """Fill each PairOverlap.mbar_o with the matrix nearest-neighbor off-diagonal
    O_{row, row+1}, where rows index the SORTED present states. The matrix is
    built over sorted(by_state); pair (i, j=i+1) maps to consecutive matrix rows.
    """
    states = sorted(by_state.keys())
    pos = {s: k for k, s in enumerate(states)}
    n = matrix.shape[0]
    for p in pairs:
        ri = pos.get(p.i)
        rj = pos.get(p.j)
        if ri is None or rj is None:
            continue
        if 0 <= ri < n and 0 <= rj < n:
            p.mbar_o = float(matrix[ri, rj])
