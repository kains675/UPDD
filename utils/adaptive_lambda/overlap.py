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
#   col0=stateid col1=temperature col2=direction col3=lambda1 col4=lambda2
#   col5=alpha col6=u0 col7=w0 col8=potE col9=pertE col10=trash
_COL_STATEID = 0
_COL_LAMBDA1 = 3
_COL_LAMBDA2 = 4
_COL_ALPHA = 5
_COL_U0 = 6
_COL_W0 = 7
_COL_POTE = 8
_COL_PERTE = 9
_N_OUT_COLS = 11

# UWHAM β convention — EXACTLY the constant atom_openmm.uwham / the postprocess
# (``trackb_uwham_postprocess.py::_calculate_uwham_multi_intermediate``) uses:
# ``bet = 1.0 / (0.001986209 * tempt)``. Bit-equivalence (BC4) demands
# this exact kB, NOT a rounded 0.0019872041 — a different β silently shifts
# every reduced potential.
_UWHAM_KB_KCAL = 0.001986209

# Field order for the full-row sample arrays (build_softcore_neg_pot consumes
# them by name via this index map).
_SAMPLE_FIELDS = ("lambda1", "lambda2", "alpha", "u0", "w0", "potE", "pertE")


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


def extract_samples_by_state(
    run_dir: str,
    jobname: str,
    direction: str,
    warmup_cycles: int = 5,
) -> Dict[int, np.ndarray]:
    """Like ``extract_pertE_by_state`` but returns the FULL per-sample soft-core
    columns: ``{state_id: array[n_samples, 7]}`` where the 7 columns are, in
    order, (lambda1, lambda2, alpha, u0, w0, potE, pertE) = ``_SAMPLE_FIELDS``.

    Same across-walker aggregation + per-walker warmup discard as
    ``extract_pertE_by_state`` (one row == one async_re cycle; the first
    ``warmup_cycles`` rows of EACH walker are dropped). Needed by the TRUE
    soft-core reduced-potential reconstruction (P6) which requires each
    sample's OWN sampled bias params + potE (cols 3-8), not just pertE (col 9)."""
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

    by_state: Dict[int, List[Tuple[float, ...]]] = {}
    for wdir in walker_dirs:
        if not os.path.isdir(wdir):
            continue
        out_path = os.path.join(wdir, out_basename)
        if not os.path.isfile(out_path):
            continue
        rows = _read_out_rows_full(out_path)
        # Discard the first warmup_cycles rows of THIS walker.
        rows = rows[warmup_cycles:]
        for sid, sample in rows:
            by_state.setdefault(sid, []).append(sample)

    return {
        sid: np.asarray(vals, dtype=float).reshape(-1, len(_SAMPLE_FIELDS))
        for sid, vals in by_state.items()
    }


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


def _read_out_rows_full(out_path: str) -> List[Tuple[int, Tuple[float, ...]]]:
    """Read an ATM ``.out`` file → list of ``(stateid, sample_tuple)`` where
    ``sample_tuple`` carries the per-sample soft-core columns in ``_SAMPLE_FIELDS``
    order: (lambda1, lambda2, alpha, u0, w0, potE, pertE).

    ``_read_out_rows`` keeps only (stateid, pertE) — the linear-λ surrogate's
    inputs. The TRUE soft-core reduced potential (P6 / BC2) additionally
    needs the walker's OWN sampled (λ1,λ2,α,u0,w0,potE) to subtract the per-row
    bias (``e0 = potE − _bias_fcn(...)``), so this reader keeps cols 3-8 too.
    Non-numeric / short lines are skipped defensively (matches _read_out_rows)."""
    rows: List[Tuple[int, Tuple[float, ...]]] = []
    with open(out_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < _N_OUT_COLS:
                continue
            try:
                sid = int(float(parts[_COL_STATEID]))
                sample = (
                    float(parts[_COL_LAMBDA1]),
                    float(parts[_COL_LAMBDA2]),
                    float(parts[_COL_ALPHA]),
                    float(parts[_COL_U0]),
                    float(parts[_COL_W0]),
                    float(parts[_COL_POTE]),
                    float(parts[_COL_PERTE]),
                )
            except (ValueError, IndexError):
                continue
            rows.append((sid, sample))
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
# TRUE ATM soft-core reduced potential (P6 — reuses atom_openmm.uwham)
# ---------------------------------------------------------------------------
def _import_uwham_bias():
    """Lazily import the ATM soft-core bias functions from ``atom_openmm.uwham``
    (the SAME helpers the postprocess UWHAM estimator trusts — R-12
    reuse, NOT a re-implementation of the soft-core math).

    Returns ``(_bias_fcn, _npot_fcn)`` or ``None`` if ``atom_openmm`` is not
    importable in this env (it lives only in the ``atm`` conda env, not
    ``qmmm``). Callers MUST treat ``None`` as "true soft-core unavailable" and
    degrade gracefully (return None matrix) — they must NOT fall back to the
    deleted linear-λ surrogate, because a surrogate O is the physics-WRONG
    operator in the soft-core region (P6 BC5)."""
    try:
        from atom_openmm.uwham import _bias_fcn, _npot_fcn  # type: ignore
        return _bias_fcn, _npot_fcn
    except Exception:
        return None


def build_softcore_neg_pot(
    samples_by_state: Dict[int, np.ndarray],
    target_states: List[int],
    schedule: Dict[str, Any],
    temperature_K: float = 300.0,
) -> np.ndarray:
    """Reconstruct the TRUE ATM soft-core NEGATIVE reduced-potential matrix
    ``neg_pot[N, K]`` — bit-equivalent (BC1, atol≤1e-9) to the matrix
    ``trackb_uwham_postprocess.py::_calculate_uwham_multi_intermediate`` feeds
    ``_uwham_r``.

    Procedure (mirrors postprocess line-for-line, reusing the SAME
    ``atom_openmm.uwham._bias_fcn`` / ``_npot_fcn``):
      1. Concatenate the per-sample rows of ``samples_by_state`` in ascending
         state-id order (the matrix-ROW order). Each row carries its OWN sampled
         (λ1,λ2,α,u0,w0,potE,pertE) (``_SAMPLE_FIELDS``).
      2. ``e0[n] = potE[n] − _bias_fcn(pertE[n], λ1[n], λ2[n], α[n], u0[n], w0[n])``
         using the WALKER'S OWN sampled bias params (postprocess L342-346).
      3. For each TARGET column state ``k`` in ``target_states`` (the matrix-
         COLUMN order), ``neg_pot[:, k] = _npot_fcn(e0, pertE, β, λ1_s[k],
         λ2_s[k], α_s[k], u0_s[k], w0_s[k])`` using the per-target-state SCHEDULE
         SSOT params (BC2 — NOT walker rows), indexed by GLOBAL state id.
      β uses the EXACT postprocess kB (``_UWHAM_KB_KCAL``), BC4.

    ``samples_by_state`` :  {global_state_id: array[n_i, 7]} (``_SAMPLE_FIELDS``).
    ``target_states``     :  ordered list of GLOBAL state ids = matrix columns
                             (caller supplies the leg/adjacency ordering; the
                             postprocess leg-1 order is ``range(leg1istate+1)``).
    ``schedule``          :  per-GLOBAL-state SSOT arrays (``lambda1`` / ``lambda2``
                             / ``alpha`` / ``u0`` / ``w0``), e.g. the
                             ``_parse_cntl_schedule`` dict. Indexed by state id.

    Raises ``RuntimeError`` if ``atom_openmm.uwham`` is unavailable (the soft-core
    math MUST come from the trusted helper — no surrogate fallback)."""
    bias = _import_uwham_bias()
    if bias is None:
        raise RuntimeError(
            "atom_openmm.uwham not importable — the TRUE soft-core reduced "
            "potential requires the trusted _bias_fcn/_npot_fcn (atm env). "
            "No linear-λ surrogate fallback (P6: surrogate O is the "
            "physics-wrong operator in the soft-core region)."
        )
    _bias_fcn, _npot_fcn = bias

    # Resolve per-state SSOT arrays. The postprocess SSOT (_parse_cntl_schedule)
    # uses 'lambda1'/'lambda2'; the asyncre get_schedule dict uses
    # 'lambdas_1'/'lambdas_2'. Accept either (same physical column).
    def _sched_arr(*names: str) -> np.ndarray:
        for nm in names:
            if nm in schedule and schedule[nm] is not None:
                return np.asarray(schedule[nm], dtype=float)
        raise KeyError(
            f"schedule must expose one of {names} for the true soft-core "
            f"reconstruction (got keys {sorted(schedule)})"
        )

    lam1_s = _sched_arr("lambda1", "lambdas_1")
    lam2_s = _sched_arr("lambda2", "lambdas_2")
    alpha_s = _sched_arr("alpha")
    u0_s = _sched_arr("u0")
    w0_s = _sched_arr("w0")

    # Row order = samples concatenated in ascending state-id (matches
    # pandas concat of r0..rN read in ascending stateid; for the per-direction
    # diagnostic any stable order works since overlap is row-permutation
    # invariant, but ascending matches the postprocess leg partition).
    row_states = sorted(samples_by_state.keys())
    blocks = [np.asarray(samples_by_state[s], dtype=float).reshape(
        -1, len(_SAMPLE_FIELDS)) for s in row_states]
    if not blocks:
        return np.zeros((0, len(target_states)), dtype=float)
    data = np.concatenate(blocks, axis=0)
    f = {name: idx for idx, name in enumerate(_SAMPLE_FIELDS)}
    pertE = data[:, f["pertE"]]
    potE = data[:, f["potE"]]

    beta = 1.0 / (_UWHAM_KB_KCAL * temperature_K)

    # e0 = potE − bias(OWN sampled params) — per-sample SCALAR bias subtraction.
    # MUST loop per-sample exactly like the postprocess (L341-346): upstream
    # ``_bias_fcn`` has ``if alpha > 0`` which is NOT array-safe, so the OWN-row
    # bias (where α/u0/w0 vary per sample) cannot be vectorized over rows. The
    # per-TARGET-state column build below CAN vectorize (scalar schedule params).
    e0 = potE.copy()
    lam1_own = data[:, f["lambda1"]]
    lam2_own = data[:, f["lambda2"]]
    alpha_own = data[:, f["alpha"]]
    u0_own = data[:, f["u0"]]
    w0_own = data[:, f["w0"]]
    n = e0.shape[0]
    for i in range(n):
        e0[i] -= _bias_fcn(
            pertE[i], lam1_own[i], lam2_own[i],
            alpha_own[i], u0_own[i], w0_own[i],
        )

    neg_pot = np.zeros((n, len(target_states)), dtype=float)
    for col, st in enumerate(target_states):
        neg_pot[:, col] = _npot_fcn(
            e0, pertE, beta,
            lam1_s[st], lam2_s[st], alpha_s[st],
            u0_s[st], w0_s[st],
        )
    return neg_pot


# ---------------------------------------------------------------------------
# MBAR overlap matrix (optional — pymbar auto-detected)
# ---------------------------------------------------------------------------
def mbar_overlap_matrix(
    samples_by_state: Dict[int, np.ndarray],
    schedule: Optional[Dict[str, Any]] = None,
    temperature_K: float = 300.0,
) -> Optional[np.ndarray]:
    """Compute the MBAR overlap matrix on the TRUE ATM soft-core reduced
    potential IF both pymbar AND atom_openmm.uwham are importable; otherwise
    return None (BC remains the always-available gate proxy — this must NOT
    hard-depend on either).

    The reduced potential is the genuine soft-core ``neg_pot`` built by
    ``build_softcore_neg_pot`` (reusing ``atom_openmm.uwham._bias_fcn`` /
    ``_npot_fcn``), NOT the deleted linear-λ surrogate ``u = λ·pertE`` (verdict
    P6: omitting the soft-core bias precisely where it defines the states is the
    WRONG operator, not an approximation). ``u_kn = −neg_pot.T`` matches the
    postprocess P11 convention (``trackb_uwham_postprocess.py::_compute_mbar_
    overlap_matrix``: neg_pot is the NEGATIVE reduced potential, MBAR wants the
    POSITIVE).

    ``samples_by_state`` is the FULL-column dict from ``extract_samples_by_state``
    ({state: array[n,7]}, ``_SAMPLE_FIELDS`` order). ``schedule`` must expose the
    per-GLOBAL-state ``lambda1``/``lambda2``/``alpha``/``u0``/``w0`` SSOT arrays.

    Returns the (n_states × n_states) overlap matrix over the present states in
    sorted order, or None if pymbar/atom_openmm is unavailable, ``schedule`` is
    missing the soft-core arrays, or the MBAR solve fails (swallowed — MBAR is
    optional; the caller falls back to BC). It is UNWIRED from every production
    gate (BC3) — failure here never blocks anything.

    Solver: MBAR is constructed with ``solver_protocol="robust"`` +
    ``maximum_iterations=20000`` + ``relative_tolerance=1e-8``. The pymbar 4.2
    DEFAULT solver is a first-order self-consistent/hybr root-finder that STALLS
    on the huge-dynamic-range two-copy ATS reduced potential (u_kn span ~2261
    reduced units) → "No solution to within tolerance" yet returns a non-physical
    overlap matrix (diagonal > 1, row-sum ≠ 1) WITHOUT raising. The convex
    UWHAM/MBAR objective (Tan 2012; Shirts-Chodera 2008) has a UNIQUE minimum, so
    the robust solver does not change the answer — it reaches the same fixed point
    the default fails to find (config B==C bit-identical confirms uniqueness).
    A convergence GUARD (diagonal ≤ 1+tol AND |row-sum − 1| ≤ tol) rejects the
    non-physical default-solver garbage as ``None`` (never emits a clean 0.0,
    which would be misread as "measured, no overlap" → phantom λ-cliff densify).
    """
    try:
        from pymbar import MBAR
    except Exception:
        return None

    states = sorted(samples_by_state.keys())
    n_states = len(states)
    if n_states < 2:
        return None

    if schedule is None:
        # No SSOT soft-core arrays → cannot build the TRUE operator. Refuse to
        # emit a surrogate number (P6 BC5); BC remains the gate.
        return None

    # Per-state sample counts. ``samples_by_state`` must be the FULL-column dict
    # ({state: array[n,7]}); a malformed shape (e.g. a bare 1-D pertE dict from
    # the OLD surrogate API) degrades gracefully to None — never raises, never
    # surrogate-falls-back.
    try:
        counts = [np.asarray(samples_by_state[s]).reshape(
            -1, len(_SAMPLE_FIELDS)).shape[0] for s in states]
    except (ValueError, TypeError):
        return None
    if any(c < 1 for c in counts):
        return None

    try:
        neg_pot = build_softcore_neg_pot(
            samples_by_state, target_states=states,
            schedule=schedule, temperature_K=temperature_K,
        )
    except Exception:
        # atom_openmm absent / schedule incomplete → true soft-core
        # unavailable. Swallow (MBAR optional); never surrogate-fallback.
        return None

    N_k = np.array(counts, dtype=int)
    # neg_pot is the NEGATIVE reduced potential; MBAR consumes the POSITIVE
    # reduced potential, so u_kn = -neg_pot.T (postprocess P11 convention).
    u_kn = -neg_pot.T

    try:
        # Robust solver kwargs (pymbar 4.2). The default first-order solver
        # stalls on the high-dynamic-range two-copy u_kn and returns a
        # non-physical matrix without raising; the robust trust-region path
        # reaches the unique convex minimum. pymbar 3.x lacks these kwargs →
        # graceful TypeError fallback to the default ctor (still guarded below).
        try:
            mbar = MBAR(
                u_kn, N_k,
                maximum_iterations=20000,
                relative_tolerance=1e-8,
                solver_protocol="robust",
            )
        except TypeError:
            mbar = MBAR(u_kn, N_k)
        result = mbar.compute_overlap()
        # pymbar 4.x returns a dict with key "matrix"; 3.x returns a tuple.
        if isinstance(result, dict):
            matrix = np.asarray(result["matrix"], dtype=float)
        else:
            # (eigenvalues, eigenvectors, matrix) tuple form
            matrix = np.asarray(result[-1], dtype=float)
        # Convergence guard: a converged overlap matrix is (near) doubly
        # stochastic — every diagonal ≤ 1 and every row sums to 1. A stalled
        # default solver emits garbage (diag ≫ 1, row-sum ≫ 1) WITHOUT raising;
        # reject it as non-convergence (None) so the caller BC-fallbacks rather
        # than reading a phantom 0.0 adjacent overlap.
        guard_tol = 1e-3
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            return None
        if not np.all(np.isfinite(matrix)):
            return None
        if float(np.diag(matrix).max()) > 1.0 + guard_tol:
            return None
        row_sums = matrix.sum(axis=1)
        if float(np.max(np.abs(row_sums - 1.0))) > guard_tol:
            return None
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
            # TRUE soft-core path: needs the full per-sample columns (cols 3-8 +
            # pertE), not just pertE. If atom_openmm.uwham / the soft-core
            # schedule arrays are unavailable, mbar_overlap_matrix returns None
            # and BC (already computed above) remains the gate proxy — NO
            # linear-λ surrogate fallback (P6).
            samples_by_state = extract_samples_by_state(
                run_dir, jobname, direction, warmup_cycles=warmup_cycles,
            )
            matrix = mbar_overlap_matrix(
                samples_by_state, schedule=schedule,
                temperature_K=temperature_K,
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
