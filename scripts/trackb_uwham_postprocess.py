"""Track B Plan A UWHAM post-processing — free-leg ΔΔG_int aggregation.

Consumes the asyncre ``r{0..N-1}/trackb.out`` sample data emitted by
``trackb_production_v2_asyncre.py`` / ``trackb_production_v2_1_upstream.py``
and produces ``ddint_free.json`` (ranking-only).

State-count-agnostic (2026-06-05, λ-densify campaign): the free leg may
be analyzed either as a COMBINED leg OR (REVISED architecture) after a
PER-DIRECTION split + merge — same path as the bound legs. It may carry the
canonical 22-state schedule, the DEPRECATED densified 34-state schedule
(forensic only — Factor-B broken) OR the REVISED densified 38-state schedule
(revised ladder fix, 2026-06-05).
The expected replica/state count + the λ ladder are DERIVED from the leg's own
``trackb_asyncre.cntl`` — NEVER a hardcoded 22. The BOUND leg stays 22
(= 11+11 per-direction split); the REVISED free leg is 38 (= 19+19
per-direction split). Both are analyzed by
``trackb_per_direction_production.merge_per_direction_outputs`` + this module's
leg path after the merge.

WARNING — C.1 UWHAM per-state SSOT (BLOCKING). The densified38
schedule RAMPS the per-state ALPHA (0.10→0.25) and U0 (110→82). Upstream
``atom_openmm.uwham.calculate_uwham`` DEFAULTS alpha/u0/λ/w0 to 22-STATE
HARDCODED values — analyzing a 38-state ramped run with those defaults would
SILENTLY BIAS ΔG (no NaN; a silent-wrong-state failure mode). This module therefore PARSES the
per-state (LAMBDA1/LAMBDA2/ALPHA/U0/W0COEFF) arrays from each leg's PRODUCTION
cntl (the single source of truth shared by simulation + analysis) and passes
them EXPLICITLY to the estimator. A missing/unparseable cntl is a HARD ERROR
(default ``require_cntl_schedule=True``); the ``--allow-legacy-defaults``
escape hatch is forensic-only and flags the result UNVERIFIED.

Schedule semantics consumed here:
  * canonical22 — LAMBDAS ascending 0→0.5 then descending 0.5→0; INTERMEDIATE
    at the two λ=0.5 midpoints (states 10 & 11); DIRECTION +1 then -1.
  * densified34 (DEPRECATED) — 10 linear forward states + 7-window ilogistic
    anneal ladder (INTERMEDIATE=1 across all 7) + whole-tuple-reversed backward
    half. INTERMEDIATE==1 spans states 10..16 (fwd) and 17..23 (bwd).
  * densified38 (REVISED, production) — 10 linear forward states + 9-window
    ilogistic anneal ladder (INTERMEDIATE=1 across all 9, per-state α/U0 ramp)
    + whole-tuple-reversed backward half. INTERMEDIATE==1 spans states 10..18
    (fwd) and 19..27 (bwd). DIRECTION +1 (states 0..18) then -1 (states 19..37).

``atom_openmm.uwham.calculate_uwham`` partitions the forward / backward legs
using ``leg1istate = where(intermd==1)[0][0]`` and
``leg2istate = where(intermd==1)[0][1]`` — i.e. the FIRST TWO intermediate
indices. That is correct for the canonical 22-state schedule (intermediates at
10, 11) but WRONG for the densified 34-state schedule (intermediates at
10..16 and 17..23 → upstream would split the forward leg at state 10 and
discard the entire ladder). For the densified schedule this module therefore
uses a vendored, source-equivalent ``_calculate_uwham_multi_intermediate``
that derives the leg boundaries from the DIRECTION column instead
(leg1istate = last forward state, leg2istate = first backward state). The
22-state path still calls the upstream ``calculate_uwham`` unchanged.

``atom_openmm.uwham.calculate_uwham`` returns
  (dgb, ddgb, dgbind1, dgbind2, samples_per_replica)
where
  dgbind1 = forward-leg ΔG (state 0 → last forward state, DIRECTION +1)
  dgbind2 = backward-leg ΔG (last state → first backward state, DIRECTION -1)
  dgb     = dgbind1 - dgbind2  (closed-cycle leg-internal closure check)

ΔΔG_int_free  =  dgbind1_cp4  −  dgbind1_wt           (ranking)
err(ΔΔG)      =  sqrt(err_cp4² + err_wt²)             (single-run diagnostic)

When multiple independent-seed replicates are supplied (``--cp4-leg`` /
``--wt-leg`` repeated, or ``--cp4-replicate-legs`` / ``--wt-replicate-legs``
comma-lists), the PRIMARY reported uncertainty is σ_btwn (inter-replicate
std / SEM of dgbind1 across replicates), NOT the single-run UWHAM analytic
error (which is autocorrelated and underestimates, 2026-06-05). The
single-run analytic error + per-replicate block-bootstrap CI are retained as
diagnostics only.

USAGE (single-run, backward compatible):
    python scripts/trackb_uwham_postprocess.py \\
        --cp4-leg outputs/_trackb/production_v2_1/cp4/free \\
        --wt-leg  outputs/_trackb/production_v2_1/wt/free  \\
        --output  outputs/_trackb/production_v2_1/ddint_free.json

USAGE (multi-replicate, σ_btwn primary):
    python scripts/trackb_uwham_postprocess.py \\
        --cp4-replicate-legs out/cp4/free_s1,out/cp4/free_s2,out/cp4/free_s3 \\
        --wt-replicate-legs  out/wt/free_s1,out/wt/free_s2,out/wt/free_s3  \\
        --output  out/ddint_free.json

``--mintimeid`` defaults to 100 (= equilibration discard;
confirmed cp4_free dgbind1 89.58→91.97 and wt closure 1.43→0.12 at
mintimeid=100). Pass ``--mintimeid 0`` to include all cycles (legacy "all"
behavior; uwham treats 0 the same as None internally because timeid starts at
1, so the explicit None sentinel is preserved via ``--mintimeid -1``).

Ranking regime only. No absolute ΔG claims. Magotti SSOT comparison is
SIGN ONLY, never magnitude. Free leg alone CANNOT give ΔΔG_bind — that
requires bound-leg production as well (ΔΔG_bind = ΔΔG_int_bound −
ΔΔG_int_free per Gallicchio ATM/ABFE thermodynamic cycle).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
import time
from typing import Any, Dict, List, Optional, Tuple


# Equilibration-discard default (2026-06-05). 100 cycles = 500 ps
# warmup discard at 5 ps/cycle. mintimeid=None sentinel (include all) is still
# reachable via --mintimeid -1.
DEFAULT_MINTIMEID = 100

# Block-bootstrap defaults.
DEFAULT_N_BLOCKS = 5
DEFAULT_N_BOOTSTRAP = 1000


def _import_uwham():
    """Lazy import of atom_openmm.uwham (requires conda env `atm`)."""
    try:
        from atom_openmm import uwham  # noqa: F401
        return uwham
    except ImportError as e:
        sys.stderr.write(
            "ERROR: atom_openmm.uwham not importable. Activate the `atm` "
            "conda env first:\n"
            "    conda activate atm\n"
            f"Original error: {e}\n"
        )
        sys.exit(2)


def _statistical_inefficiency():
    """Return pymbar's statistical-inefficiency callable, or None.

    pymbar 4.x exposes ``pymbar.timeseries.statistical_inefficiency``
    (snake_case); pymbar 3.x exposes ``statisticalInefficiency`` (camelCase).
    Probe both; return None if pymbar is absent so block-count auto-tuning
    degrades gracefully to the fixed default.
    """
    try:
        from pymbar import timeseries  # type: ignore
    except Exception:
        return None
    for name in ("statistical_inefficiency", "statisticalInefficiency"):
        fn = getattr(timeseries, name, None)
        if callable(fn):
            return fn
    return None


# ---------------------------------------------------------------------------
# State-count derivation (Path 2026-06-05 — NOT hardcoded 22).
# ---------------------------------------------------------------------------
def _parse_cntl_schedule(cntl_path: str) -> Optional[Dict[str, List[float]]]:
    """Parse the per-state schedule arrays from a Track B ``.cntl``.

    Returns a dict with int/float lists for keys ``directions`` / ``intermd``
    / ``lambda1`` / ``lambda2`` / ``alpha`` / ``u0`` / ``w0`` / ``n_states``,
    or None if the cntl is absent / malformed. Used to derive the expected
    replica count + (for the densified schedule) the explicit UWHAM override
    arrays.
    """
    if not os.path.isfile(cntl_path):
        return None
    keymap = {
        "DIRECTION": "directions",
        "INTERMEDIATE": "intermd",
        "LAMBDA1": "lambda1",
        "LAMBDA2": "lambda2",
        "ALPHA": "alpha",
        "U0": "u0",
        "W0COEFF": "w0",
    }
    out: Dict[str, List[float]] = {}
    with open(cntl_path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or "=" not in stripped:
                continue
            key, _, value = stripped.partition("=")
            key = key.strip()
            if key not in keymap:
                continue
            body = value.strip().strip("'\"")
            try:
                vals = [float(tok.strip()) for tok in body.split(",")
                        if tok.strip() != ""]
            except ValueError:
                continue
            out[keymap[key]] = vals
    if "intermd" not in out:
        return None
    n = len(out["intermd"])
    out["n_states"] = n
    # Sanity: every parsed per-state array must agree in length.
    for k, v in out.items():
        if k == "n_states":
            continue
        if len(v) != n:
            return None
    return out


def _derive_expected_replicas(leg_dir: str, jobname: str) -> int:
    """Derive the expected replica/state count for a COMBINED leg.

    Priority:
      1. Parse ``leg_dir/{jobname}_asyncre.cntl`` INTERMEDIATE length.
      2. Fall back to counting ``leg_dir/r*/`` directories with a
         ``{jobname}.out`` file.

    Never returns a hardcoded literal — the free leg may be 22 (canonical) or
    34 (densified). Raises FileNotFoundError if neither source resolves.
    """
    cntl_path = os.path.join(leg_dir, jobname + "_asyncre.cntl")
    sched = _parse_cntl_schedule(cntl_path)
    if sched is not None:
        return int(sched["n_states"])
    # Fallback: count r{i}/ dirs that hold a .out file.
    rdirs = sorted(
        glob.glob(os.path.join(leg_dir, "r*", f"{jobname}.out"))
    )
    if rdirs:
        # r-index set must be contiguous 0..N-1; count the highest index + 1.
        idxs = []
        for p in rdirs:
            m = re.search(r"r(\d+)" + re.escape(os.sep) + re.escape(jobname),
                          p)
            if m:
                idxs.append(int(m.group(1)))
        if idxs:
            return max(idxs) + 1
    raise FileNotFoundError(
        f"cannot derive replica count for {leg_dir}: no cntl with "
        f"INTERMEDIATE and no r*/{jobname}.out files found"
    )


# ---------------------------------------------------------------------------
# Vendored multi-intermediate UWHAM (densified schedule support).
# Source-equivalent to atom_openmm.uwham.calculate_uwham EXCEPT the leg
# boundary derivation: upstream uses the FIRST TWO intermediate indices
# (correct only for a 2-intermediate 22-state schedule). For a schedule with
# >2 intermediate states (densified34: 14 intermediates), the forward/backward
# split must come from the DIRECTION column: leg1istate = last DIRECTION==+1
# state, leg2istate = first DIRECTION==-1 state.
# ---------------------------------------------------------------------------
def _calculate_uwham_multi_intermediate(
    rundir: str,
    jobname: str,
    schedule: Dict[str, List[float]],
    mintimeid: Optional[int],
    maxtimeid: Optional[int],
) -> Tuple[float, float, float, float, int]:
    """Drop-in replacement for ``uwham.calculate_uwham`` that supports a
    schedule with more than two INTERMEDIATE states.

    Mirrors the upstream UWHAM math line-for-line (``_uwham_r``, ``_bias_fcn``,
    ``_npot_fcn`` are reused unchanged from atom_openmm.uwham) but partitions
    the legs by the DIRECTION column. Returns the same 5-tuple
    (dgb, ddgb, dgbind1, dgbind2, samples_per_replica).
    """
    import numpy as np
    import pandas as pd

    uwham = _import_uwham()

    intermd = np.array(schedule["intermd"])
    lambda1 = np.array(schedule["lambda1"])
    lambda2 = np.array(schedule["lambda2"])
    alpha = np.array(schedule["alpha"])
    u0 = np.array(schedule["u0"])
    w0 = np.array(schedule["w0"])
    directions = np.array(schedule["directions"])
    nstates = len(intermd)

    tempt = np.array([300])
    bet = 1.0 / (0.001986209 * tempt)

    # Leg boundaries from DIRECTION (not the first-two-intermediate heuristic).
    # Forward leg = DIRECTION >= 0; backward leg = DIRECTION < 0.
    fwd_idx = np.where(directions >= 0)[0]
    bwd_idx = np.where(directions < 0)[0]
    if len(fwd_idx) == 0 or len(bwd_idx) == 0:
        raise ValueError(
            "schedule must contain both forward (DIRECTION>=0) and backward "
            "(DIRECTION<0) states for the bidirectional UWHAM cycle"
        )
    leg1istate = int(fwd_idx[-1])   # last forward state (λ=0.5 turn point)
    leg2istate = int(bwd_idx[0])    # first backward state

    columns = [
        "stateid", "temperature", "direction", "lambda1", "lambda2",
        "alpha", "u0", "w0", "potE", "pertE", "trash",
    ]
    datafiles = [
        os.path.join(rundir, f"r{i}", f"{jobname}.out") for i in range(nstates)
    ]
    dfs = []
    for file in datafiles:
        df = pd.read_csv(file, sep=r"\s+", header=None, names=columns,
                         index_col=False)
        df["timeid"] = np.arange(1, len(df) + 1)
        dfs.append(df)
    data = pd.concat(dfs, ignore_index=True)
    data["bet"] = 1.0 / (0.001986209 * data["temperature"])

    timemask = np.ones(len(data), dtype=bool)
    if mintimeid is not None:
        timemask = timemask & (data["timeid"] >= mintimeid)
    if maxtimeid is not None:
        timemask = timemask & (data["timeid"] <= maxtimeid)
    if not np.any(timemask):
        raise ValueError(
            f"No data found in the given time range. Time values range "
            f"between {data['timeid'].min()} and {data['timeid'].max()} "
            f"and the user requested data between {mintimeid} and {maxtimeid}."
        )
    nsamples = len(data[timemask])
    samplesperreplica = nsamples // nstates

    # ---- Leg 1 (forward, states 0..leg1istate) ----
    data1 = data[timemask & (data["stateid"] <= leg1istate)]
    mtempt = len(bet)
    leg1stateids = np.arange(leg1istate + 1)
    mlam = len(leg1stateids)
    m = mlam * mtempt
    N = len(data1)
    e0 = data1["potE"].values.copy()
    for i in range(N):
        e0[i] -= uwham._bias_fcn(
            data1["pertE"].iloc[i], data1["lambda1"].iloc[i],
            data1["lambda2"].iloc[i], data1["alpha"].iloc[i],
            data1["u0"].iloc[i], data1["w0"].iloc[i],
        )
    neg_pot = np.zeros((N, m))
    sid = 0
    for be in leg1stateids:
        for te in range(mtempt):
            neg_pot[:, sid] = uwham._npot_fcn(
                e0, data1["pertE"].values, bet[te],
                lambda1[be], lambda2[be], alpha[be], u0[be], w0[be],
            )
            sid += 1
    statelabels = data1["stateid"].values.astype(int) + 1
    out = uwham._uwham_r(label=statelabels, logQ=neg_pot,
                         ufactormax=1, ufactormin=1)
    ze = np.array(out["ze"]).reshape(mtempt, mlam)
    ve = np.array(out["ve"]).reshape(mtempt, mlam)
    dgbind1 = (-ze[:, -1] / bet) - (-ze[:, 0] / bet)
    ddgbind1 = np.sqrt(ve[:, -1] + ve[:, 0]) / bet

    # ---- Leg 2 (backward, states leg2istate..nstates-1, reversed) ----
    data1 = data[timemask & (data["stateid"] >= leg2istate)]
    leg2stateids = np.arange(leg2istate, nstates)[::-1]
    mlam = len(leg2stateids)
    m = mlam * mtempt
    N = len(data1)
    e0 = data1["potE"].copy()
    for i in range(N):
        e0.iloc[i] -= uwham._bias_fcn(
            data1["pertE"].iloc[i], data1["lambda1"].iloc[i],
            data1["lambda2"].iloc[i], data1["alpha"].iloc[i],
            data1["u0"].iloc[i], data1["w0"].iloc[i],
        )
    neg_pot = np.zeros((N, m))
    sid = 0
    for be in leg2stateids:
        for te in range(mtempt):
            neg_pot[:, sid] = uwham._npot_fcn(
                e0=e0, epert=data1["pertE"], bet=bet[te],
                lam1=lambda1[be], lam2=lambda2[be], alpha=alpha[be],
                u0=u0[be], w0=w0[be],
            )
            sid += 1
    statelabels = nstates - data1["stateid"]
    out = uwham._uwham_r(label=statelabels, logQ=neg_pot,
                         ufactormax=1, ufactormin=1)
    ze = out["ze"].reshape(mtempt, mlam)
    ve = out["ve"].reshape(mtempt, mlam)
    dgbind2 = (-ze[:, -1] / bet) - (-ze[:, 0] / bet)
    ddgbind2 = np.sqrt(ve[:, -1] + ve[:, 0]) / bet

    dgb = dgbind1 - dgbind2
    ddgb = np.sqrt(ddgbind2 * ddgbind2 + ddgbind1 * ddgbind1)
    return (float(dgb[0]), float(ddgb[0]), float(dgbind1[0]),
            float(dgbind2[0]), int(samplesperreplica))


# ---------------------------------------------------------------------------
# Block bootstrap (per-replicate convergence diagnostic).
# ---------------------------------------------------------------------------
def _auto_block_count(
    leg_dir: str,
    jobname: str,
    n_states: int,
    default_blocks: int,
) -> int:
    """Auto-tune the block count from per-replica autocorrelation.

    Uses pymbar's statistical-inefficiency (g) on the forward-leg replica's
    pertE series to size the block count so each block holds >~ g samples.
    Falls back to ``default_blocks`` when pymbar is absent or the series is
    too short. Bounded to [2, 10] to keep block estimates well-conditioned.
    """
    stat_ineff = _statistical_inefficiency()
    if stat_ineff is None:
        return default_blocks
    import numpy as np
    # Use the first forward replica's pertE column (col index 9) as the
    # representative timeseries. Single representative series is enough to
    # size the block count (autocorr is dominated by the gap states; a
    # conservative g is fine).
    rep_out = os.path.join(leg_dir, "r0", f"{jobname}.out")
    if not os.path.isfile(rep_out):
        return default_blocks
    try:
        arr = np.loadtxt(rep_out, usecols=(9,))
    except Exception:
        return default_blocks
    if arr.ndim == 0 or arr.size < 10:
        return default_blocks
    try:
        g = float(stat_ineff(arr))
    except Exception:
        return default_blocks
    if g <= 0:
        return default_blocks
    n_independent = max(1.0, arr.size / g)
    n_blocks = int(min(10, max(2, round(n_independent / 2.0))))
    return n_blocks


def _block_bootstrap_dgbind1(
    leg_dir: str,
    jobname: str,
    schedule: Optional[Dict[str, List[float]]],
    mintimeid: Optional[int],
    maxtimeid: Optional[int],
    n_states: int,
    n_blocks: int,
    n_bootstrap: int,
    rng_seed: int,
    analyze_fn,
) -> Dict[str, Any]:
    """Per-replicate block-bootstrap CI of dgbind1 (convergence diagnostic).

    Splits the production cycles into ``n_blocks`` contiguous time-blocks,
    runs the leg analysis on each block window (block b = cycles [lo_b, hi_b]
    via mintimeid/maxtimeid), then bootstrap-resamples the block-level
    dgbind1 values ``n_bootstrap`` times to form a CI. Returns the bootstrap
    mean / std / 2.5–97.5 percentile CI + the per-block dgbind1 list.

    ``analyze_fn(leg_dir, jobname, schedule, lo, hi, n_states)`` must return a
    5-tuple (dgb, ddgb, dgbind1, dgbind2, samples). Blocks that raise (too few
    samples) are skipped; if <2 usable blocks survive the diagnostic returns
    ``{"status": "insufficient_blocks", ...}``.
    """
    import numpy as np

    # Determine the cycle range present in r0 (post-equilibration if mintimeid).
    rep_out = os.path.join(leg_dir, "r0", f"{jobname}.out")
    if not os.path.isfile(rep_out):
        return {"status": "no_replica_out", "leg_dir": leg_dir}
    with open(rep_out) as fh:
        n_cycles = sum(1 for ln in fh if ln.strip())
    lo_start = (mintimeid if (mintimeid is not None and mintimeid > 0) else 1)
    hi_end = (maxtimeid if maxtimeid is not None else n_cycles)
    if hi_end <= lo_start:
        return {"status": "empty_range", "lo": lo_start, "hi": hi_end}
    span = hi_end - lo_start + 1
    if span < n_blocks:
        n_blocks = max(2, span)
    edges = np.linspace(lo_start, hi_end + 1, n_blocks + 1).astype(int)

    block_dgbind1: List[float] = []
    for b in range(n_blocks):
        lo = int(edges[b])
        hi = int(edges[b + 1]) - 1
        if hi < lo:
            continue
        try:
            _, _, dgb1, _, _ = analyze_fn(leg_dir, jobname, schedule, lo, hi,
                                          n_states)
        except Exception:
            continue
        if dgb1 == dgb1:  # not NaN
            block_dgbind1.append(float(dgb1))

    if len(block_dgbind1) < 2:
        return {
            "status": "insufficient_blocks",
            "n_blocks_requested": n_blocks,
            "n_blocks_usable": len(block_dgbind1),
            "block_dgbind1": block_dgbind1,
        }

    blocks = np.array(block_dgbind1)
    rng = np.random.default_rng(rng_seed)
    boot_means = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        sample = rng.choice(blocks, size=len(blocks), replace=True)
        boot_means[i] = sample.mean()
    return {
        "status": "ok",
        "n_blocks": len(block_dgbind1),
        "n_bootstrap": n_bootstrap,
        "block_dgbind1": block_dgbind1,
        "block_mean_dgbind1": float(blocks.mean()),
        "block_std_dgbind1": float(blocks.std(ddof=1)),
        "bootstrap_mean_dgbind1": float(boot_means.mean()),
        "bootstrap_std_dgbind1": float(boot_means.std(ddof=1)),
        "bootstrap_ci95_dgbind1": [
            float(np.percentile(boot_means, 2.5)),
            float(np.percentile(boot_means, 97.5)),
        ],
        "rng_seed": rng_seed,
    }


# ---------------------------------------------------------------------------
# Per-leg analysis (schedule-aware dispatch + block bootstrap).
# ---------------------------------------------------------------------------
def _leg_analysis_window(
    leg_dir: str,
    jobname: str,
    schedule: Optional[Dict[str, List[float]]],
    mintimeid: Optional[int],
    maxtimeid: Optional[int],
    n_states: int,
) -> Tuple[float, float, float, float, int]:
    """Dispatch to upstream calculate_uwham (canonical) or the vendored
    multi-intermediate variant (densified) for one analysis window.

    WARNING — C.1 UWHAM per-state SSOT (BLOCKING). ``schedule``
    MUST be the per-state arrays parsed from the leg's PRODUCTION cntl (the
    single source of truth shared by simulation and analysis). It is passed
    EXPLICITLY to ``calculate_uwham`` (or the vendored multi-intermediate
    variant) so the estimator re-evaluates samples with the SAME
    (λ1/λ2/α/u0/w0) bias function the simulation used.

    ``schedule`` is NEVER None here: upstream ``calculate_uwham`` defaults
    α/u0/λ/w0 to 22-STATE HARDCODED values (alpha=0.10, u0=110, 22-state λ
    arrays). A 38-state ramped run (α 0.10→0.25, u0 110→82) analyzed with those
    defaults would be SILENTLY BIASED (no NaN — a silent-wrong
    pattern). The caller (``analyze_one_leg``) raises before reaching this
    function if the cntl cannot be parsed, so a missing/malformed cntl can
    never fall through to upstream defaults.

    The densified (>2 INTERMEDIATE) path uses the vendored
    multi-intermediate variant (DIRECTION-based leg partition); the canonical /
    custom (≤2 INTERMEDIATE) path uses upstream calculate_uwham WITH the cntl
    arrays as explicit overrides.
    """
    # C.1 guard FIRST (before any atom_openmm import) so the silent-bias
    # rejection fires regardless of whether the atm env is active.
    if schedule is None:
        raise ValueError(
            "C.1 UWHAM per-state SSOT violation: _leg_analysis_window requires "
            "an explicit schedule parsed from the production cntl (got None). "
            "Falling back to upstream 22-state α/u0/λ/w0 defaults would "
            "SILENTLY BIAS a ramped (densified38) run — per-state SSOT "
            "BLOCKING condition. Ensure the leg cntl is present + parseable."
        )
    uwham = _import_uwham()
    if sum(int(x) for x in schedule["intermd"]) > 2:
        return _calculate_uwham_multi_intermediate(
            rundir=leg_dir, jobname=jobname, schedule=schedule,
            mintimeid=mintimeid, maxtimeid=maxtimeid,
        )
    dgb, ddgb, dgbind1, dgbind2, samples = uwham.calculate_uwham(
        rundir=leg_dir, jobname=jobname,
        mintimeid=mintimeid, maxtimeid=maxtimeid,
        intermd=schedule["intermd"], lambda1=schedule["lambda1"],
        lambda2=schedule["lambda2"], alpha=schedule["alpha"],
        u0=schedule["u0"], w0=schedule["w0"],
    )
    return (float(dgb), float(ddgb), float(dgbind1), float(dgbind2),
            int(samples))


def analyze_one_leg(
    leg_dir: str,
    jobname: str,
    mintimeid: Optional[int],
    maxtimeid: Optional[int],
    block_bootstrap: bool = True,
    n_blocks: int = DEFAULT_N_BLOCKS,
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    auto_block_count: bool = True,
    rng_seed: int = 20260605,
    require_cntl_schedule: bool = True,
) -> Dict[str, Any]:
    """Run UWHAM on a single COMBINED leg directory (state-count-agnostic).

    Derives the expected replica/state count from the leg cntl; PARSES the
    per-state (λ1/λ2/α/u0/w0) arrays from that PRODUCTION cntl and passes them
    EXPLICITLY to the estimator (C.1 UWHAM per-state SSOT —
    cntl is the single source of truth shared by simulation + analysis).
    Selects upstream calculate_uwham (canonical 22-state / custom ≤2
    intermediate, WITH cntl overrides) or the vendored multi-intermediate
    variant (densified34/38, DIRECTION-based leg partition). Returns the five
    UWHAM scalars + provenance + (optionally) the per-replicate block-bootstrap
    CI.

    WARNING — C.1 (BLOCKING): when ``require_cntl_schedule`` is True (default) a
    missing/unparseable cntl is a HARD ERROR — upstream calculate_uwham would
    otherwise silently default α/u0/λ/w0 to 22-STATE values and SILENTLY BIAS
    a ramped (densified38) run (no NaN; silent-wrong-state risk). Set False ONLY for
    forensic re-analysis of a legacy 22-state run whose cntl is genuinely
    absent (the result is then UNVERIFIED and flagged in ``schedule_kind``).
    """
    if not os.path.isdir(leg_dir):
        raise FileNotFoundError(f"leg_dir not found: {leg_dir}")

    # Derive state count + parse the explicit per-state schedule from the cntl
    # (C.1 SSOT). The schedule arrays are byte-identical to what the simulation
    # used because both read the same cntl.
    expected_replicas = _derive_expected_replicas(leg_dir, jobname)
    cntl_path = os.path.join(leg_dir, jobname + "_asyncre.cntl")
    schedule = _parse_cntl_schedule(cntl_path)

    if schedule is None and require_cntl_schedule:
        raise FileNotFoundError(
            f"C.1 UWHAM per-state SSOT: cannot parse the production cntl "
            f"{cntl_path} (per-state λ1/λ2/α/u0/w0 arrays). Analyzing without "
            f"it would fall back to upstream 22-state α/u0 defaults and "
            f"SILENTLY BIAS a ramped (densified38) run (per-state SSOT "
            f"BLOCKING). Stage the cntl alongside r*/ or pass "
            f"require_cntl_schedule=False ONLY for a legacy 22-state forensic "
            f"re-analysis."
        )

    missing = [
        i for i in range(expected_replicas)
        if not os.path.isfile(os.path.join(leg_dir, f"r{i}", f"{jobname}.out"))
    ]
    if missing:
        raise FileNotFoundError(
            f"missing replica .out files in {leg_dir}: "
            f"r{missing} (expected 0..{expected_replicas - 1})"
        )

    n_intermd = (sum(int(x) for x in schedule["intermd"])
                 if schedule is not None else 2)
    if schedule is None:
        # require_cntl_schedule=False legacy path — UNVERIFIED (upstream
        # 22-state defaults). Flagged so the JSON consumer can see it.
        schedule_kind = f"UNVERIFIED_legacy_defaults_{expected_replicas}state"
    elif expected_replicas == 38 and n_intermd > 2:
        schedule_kind = "densified38"
    elif expected_replicas == 34 and n_intermd > 2:
        schedule_kind = "densified34_DEPRECATED"
    elif expected_replicas == 22:
        schedule_kind = "canonical22"
    else:
        schedule_kind = f"custom_{expected_replicas}state"

    # Choose the analysis callable. With a parsed schedule → the C.1-strict
    # window (explicit cntl overrides). Legacy fallback (schedule None +
    # require_cntl_schedule=False) → upstream defaults, UNVERIFIED.
    if schedule is not None:
        analyze_fn = _leg_analysis_window
    else:
        analyze_fn = _leg_analysis_window_legacy_defaults

    t0 = time.time()
    dgb, ddgb, dgbind1, dgbind2, samples = analyze_fn(
        leg_dir=leg_dir, jobname=jobname, schedule=schedule,
        mintimeid=mintimeid, maxtimeid=maxtimeid, n_states=expected_replicas,
    )
    wall_s = time.time() - t0

    result: Dict[str, Any] = {
        "leg_dir": leg_dir,
        "jobname": jobname,
        "mintimeid": mintimeid,
        "maxtimeid": maxtimeid,
        "n_states": expected_replicas,
        "n_intermediate_states": n_intermd,
        "schedule_kind": schedule_kind,
        "samples_per_replica": int(samples),
        "dgb_kcal": float(dgb),         # leg-internal closure (dgbind1 - dgbind2)
        "ddgb_kcal": float(ddgb),       # propagated single-run analytic err
        "dgbind1_kcal": float(dgbind1), # forward leg (state 0 → last fwd state)
        "dgbind2_kcal": float(dgbind2), # backward leg
        "wall_s": round(wall_s, 2),
    }

    if block_bootstrap:
        blocks = n_blocks
        if auto_block_count:
            blocks = _auto_block_count(leg_dir, jobname, expected_replicas,
                                       n_blocks)
        result["block_bootstrap"] = _block_bootstrap_dgbind1(
            leg_dir=leg_dir, jobname=jobname, schedule=schedule,
            mintimeid=mintimeid, maxtimeid=maxtimeid,
            n_states=expected_replicas, n_blocks=blocks,
            n_bootstrap=n_bootstrap, rng_seed=rng_seed,
            analyze_fn=analyze_fn,
        )
    return result


def _leg_analysis_window_legacy_defaults(
    leg_dir: str,
    jobname: str,
    schedule: Optional[Dict[str, List[float]]],
    mintimeid: Optional[int],
    maxtimeid: Optional[int],
    n_states: int,
) -> Tuple[float, float, float, float, int]:
    """LEGACY-ONLY analysis path: upstream calculate_uwham with NO per-state
    overrides (22-state α/u0/λ/w0 defaults).

    Reachable ONLY when ``analyze_one_leg(..., require_cntl_schedule=False)``
    is called for a legacy 22-state run whose cntl is genuinely absent. The
    result is UNVERIFIED (not C.1-SSOT) and ``schedule_kind`` flags it. NEVER
    used for densified runs — those have >2 intermediate states which upstream
    defaults cannot represent. ``schedule`` is ignored (always None here).
    """
    uwham = _import_uwham()
    dgb, ddgb, dgbind1, dgbind2, samples = uwham.calculate_uwham(
        rundir=leg_dir, jobname=jobname,
        mintimeid=mintimeid, maxtimeid=maxtimeid,
    )
    return (float(dgb), float(ddgb), float(dgbind1), float(dgbind2),
            int(samples))


def analyze_replicate_set(
    leg_dirs: List[str],
    jobname: str,
    mintimeid: Optional[int],
    maxtimeid: Optional[int],
    block_bootstrap: bool = True,
    n_blocks: int = DEFAULT_N_BLOCKS,
    n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
    auto_block_count: bool = True,
    require_cntl_schedule: bool = True,
) -> Dict[str, Any]:
    """Analyze one endpoint's set of independent-seed replicate legs.

    Each replicate leg is analyzed via ``analyze_one_leg`` (C.1 SSOT: per-state
    arrays read from each leg's own production cntl). σ_btwn is the
    inter-replicate std / SEM of dgbind1 (PRIMARY error). The
    single-run analytic error and per-replicate block-bootstrap are retained
    as diagnostics. With a single leg, σ_btwn is undefined (None) and the
    analytic / block-bootstrap error is the only available diagnostic.

    Returns a dict with ``replicates`` (per-leg results), ``mean_dgbind1`` /
    ``sigma_btwn_dgbind1`` / ``sem_dgbind1`` and provenance.
    """
    import numpy as np

    reps: List[Dict[str, Any]] = []
    for i, leg_dir in enumerate(leg_dirs):
        r = analyze_one_leg(
            leg_dir=leg_dir, jobname=jobname,
            mintimeid=mintimeid, maxtimeid=maxtimeid,
            block_bootstrap=block_bootstrap, n_blocks=n_blocks,
            n_bootstrap=n_bootstrap, auto_block_count=auto_block_count,
            rng_seed=20260605 + i,  # per-replicate deterministic bootstrap seed
            require_cntl_schedule=require_cntl_schedule,
        )
        reps.append(r)

    dgbind1_vals = np.array([r["dgbind1_kcal"] for r in reps])
    n = len(dgbind1_vals)
    mean_dgbind1 = float(dgbind1_vals.mean())
    if n >= 2:
        sigma_btwn = float(dgbind1_vals.std(ddof=1))   # inter-replicate std
        sem = float(sigma_btwn / (n ** 0.5))           # SEM
    else:
        sigma_btwn = None
        sem = None

    return {
        "n_replicates": n,
        "replicate_leg_dirs": list(leg_dirs),
        "replicates": reps,
        "mean_dgbind1_kcal": mean_dgbind1,
        "sigma_btwn_dgbind1_kcal": sigma_btwn,   # PRIMARY error (inter-rep std)
        "sem_dgbind1_kcal": sem,                 # inter-rep SEM
        # Single-run analytic error of replicate-0 (diagnostic only).
        "single_run_analytic_ddgb_kcal_diag": reps[0]["ddgb_kcal"],
    }


def _resolve_legs(single: Optional[str],
                  replicates_csv: Optional[str]) -> List[str]:
    """Resolve the leg-dir list from either --<ep>-leg (single) or
    --<ep>-replicate-legs (comma list). Returns absolute paths.

    Exactly one of the two must be provided per endpoint (the caller asserts).
    """
    if replicates_csv:
        return [os.path.abspath(p.strip())
                for p in replicates_csv.split(",") if p.strip()]
    if single:
        return [os.path.abspath(single)]
    return []


def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Track B Plan A UWHAM post-processing: free-leg ΔΔG_int "
            "(ranking-only; state-count-agnostic 22/34; σ_btwn primary)."
        )
    )
    # Single-run (backward compatible) leg dirs.
    p.add_argument(
        "--cp4-leg",
        help="cp4 endpoint free-leg dir (single run; r0/..rN-1/trackb.out)",
    )
    p.add_argument(
        "--wt-leg",
        help="wt endpoint free-leg dir (single run)",
    )
    # Multi-replicate (σ_btwn primary).
    p.add_argument(
        "--cp4-replicate-legs",
        help="comma-list of cp4 replicate leg dirs (independent seeds)",
    )
    p.add_argument(
        "--wt-replicate-legs",
        help="comma-list of wt replicate leg dirs (independent seeds)",
    )
    p.add_argument(
        "--output",
        required=True,
        help="output JSON path (ddint_free.json)",
    )
    p.add_argument(
        "--jobname",
        default="trackb",
        help="basename of per-replica .out files (default: trackb)",
    )
    p.add_argument(
        "--mintimeid",
        type=int,
        default=DEFAULT_MINTIMEID,
        help=(
            "lower bound (inclusive) timeid for the sample window — discards "
            "early equilibration cycles. Default: 100 (= 500 ps warmup "
            "discard; 2026-06-05). Pass -1 for the None sentinel "
            "(include all cycles)."
        ),
    )
    p.add_argument(
        "--maxtimeid",
        type=int,
        default=None,
        help="upper bound (inclusive) timeid for sample window. Default: end.",
    )
    p.add_argument(
        "--no-block-bootstrap",
        action="store_true",
        help="skip the per-replicate block-bootstrap CI diagnostic.",
    )
    p.add_argument(
        "--n-blocks",
        type=int,
        default=DEFAULT_N_BLOCKS,
        help=f"block count for block-bootstrap (default {DEFAULT_N_BLOCKS}).",
    )
    p.add_argument(
        "--n-bootstrap",
        type=int,
        default=DEFAULT_N_BOOTSTRAP,
        help=f"bootstrap resamples (default {DEFAULT_N_BOOTSTRAP}).",
    )
    p.add_argument(
        "--no-auto-block-count",
        action="store_true",
        help="disable pymbar statistical-inefficiency block-count auto-tuning.",
    )
    p.add_argument(
        "--allow-legacy-defaults",
        action="store_true",
        help=(
            "C.1 ESCAPE HATCH (forensic only): permit analyzing a leg whose "
            "production cntl is genuinely absent, falling back to upstream "
            "22-state α/u0/λ/w0 defaults. The result is UNVERIFIED (flagged "
            "in schedule_kind). NEVER use for densified38 — it would silently "
            "bias ΔG. Default OFF: a missing/unparseable cntl is a hard error "
            "(C.1 per-state SSOT, BLOCKING)."
        ),
    )
    p.add_argument(
        "--protocol-tag",
        default="atm-openmm async_re b+ corrected (v2.1 upstream-vendored)",
        help="protocol provenance string written into the JSON",
    )
    p.add_argument(
        "--notes",
        default=(
            "Free leg only. Bound leg structprep b+ FULL PASS (4-leg × 2-dir "
            "audit clean 2026-06-01 04:34 KST); bound-leg production pending "
            "user explicit confirmation (quadruple-gate)."
        ),
        help="free-form note string for the JSON",
    )
    args = p.parse_args()

    # mintimeid sentinel: -1 → None (include all cycles).
    mintimeid: Optional[int] = (None if args.mintimeid is not None
                                and args.mintimeid < 0 else args.mintimeid)

    cp4_legs = _resolve_legs(args.cp4_leg, args.cp4_replicate_legs)
    wt_legs = _resolve_legs(args.wt_leg, args.wt_replicate_legs)
    if not cp4_legs:
        sys.stderr.write("ERROR: provide --cp4-leg or --cp4-replicate-legs\n")
        return 2
    if not wt_legs:
        sys.stderr.write("ERROR: provide --wt-leg or --wt-replicate-legs\n")
        return 2

    block_bootstrap = not args.no_block_bootstrap
    auto_block = not args.no_auto_block_count
    require_cntl = not args.allow_legacy_defaults

    print(f"[uwham] analyzing cp4 ({len(cp4_legs)} replicate(s))")
    cp4 = analyze_replicate_set(
        leg_dirs=cp4_legs, jobname=args.jobname,
        mintimeid=mintimeid, maxtimeid=args.maxtimeid,
        block_bootstrap=block_bootstrap, n_blocks=args.n_blocks,
        n_bootstrap=args.n_bootstrap, auto_block_count=auto_block,
        require_cntl_schedule=require_cntl,
    )
    print(f"[uwham] analyzing wt ({len(wt_legs)} replicate(s))")
    wt = analyze_replicate_set(
        leg_dirs=wt_legs, jobname=args.jobname,
        mintimeid=mintimeid, maxtimeid=args.maxtimeid,
        block_bootstrap=block_bootstrap, n_blocks=args.n_blocks,
        n_bootstrap=args.n_bootstrap, auto_block_count=auto_block,
        require_cntl_schedule=require_cntl,
    )

    # ΔΔG_int_free using the forward leg (dgbind1). PRIMARY error = σ_btwn
    # (inter-replicate) when ≥2 replicates per endpoint; otherwise fall back to
    # the single-run analytic error (diagnostic) so the single-run CLI still
    # reports a number (with a clear flag).
    ddint_free = cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"]

    cp4_err = cp4["sigma_btwn_dgbind1_kcal"]
    wt_err = wt["sigma_btwn_dgbind1_kcal"]
    error_basis = "sigma_btwn"
    if cp4_err is None or wt_err is None:
        # Single-run fallback — analytic error of replicate 0 (diagnostic).
        cp4_err = cp4["replicates"][0]["ddgb_kcal"]
        wt_err = wt["replicates"][0]["ddgb_kcal"]
        error_basis = "single_run_analytic_diag"
    ddint_err = (cp4_err ** 2 + wt_err ** 2) ** 0.5

    # Closure check per endpoint (replicate-0 representative).
    closure_cp4 = abs(cp4["replicates"][0]["dgb_kcal"])
    closure_wt = abs(wt["replicates"][0]["dgb_kcal"])

    payload: Dict[str, Any] = {
        "schema_version": "trackb_ddint_free_v2",
        "regime": "ranking_only_R11",
        "protocol": args.protocol_tag,
        "mintimeid": mintimeid,
        "maxtimeid": args.maxtimeid,
        "error_basis": error_basis,
        "cp4": cp4,
        "wt": wt,
        "ddint_free_kcal": float(ddint_free),
        "ddint_err_kcal": float(ddint_err),
        "closure_cp4_kcal_abs": float(closure_cp4),
        "closure_wt_kcal_abs": float(closure_wt),
        "interpretation_note": (
            "ddint_free_kcal = cp4.mean_dgbind1 - wt.mean_dgbind1 (forward "
            "+1 leg). PRIMARY error = σ_btwn (inter-replicate std) when ≥2 "
            "replicates per endpoint; single-run analytic error is a DIAGNOSTIC "
            "fallback (autocorrelated, underestimate — Path 2026-06-05). "
            "Closure check = |dgb| per endpoint (forward minus backward; ideal "
            "~0 for free leg). Free leg ΔΔG_int ALONE does NOT give ΔΔG_bind — "
            "pair with bound-leg ΔΔG_int (ΔΔG_bind = ΔΔG_int_bound − "
            "ΔΔG_int_free, Gallicchio ATM/ABFE cycle)."
        ),
        "magotti_ssot_comparison_policy": (
            "ranking-only: SIGN comparison with Magotti 2009 anti-C3b "
            "compstatin ITC/SPR allowed (Magotti ΔΔG_bind Cp4-vs-WT ≈ -1.4 "
            "ITC / -3.0 SPR kcal/mol, Cp4 better). MAGNITUDE comparison "
            "FORBIDDEN — UPDD Track B has not passed v0.7 calibration."
        ),
        "free_leg_warning": (
            "Free leg only. Bound leg production pending. Free-leg ΔΔG_int "
            "alone is NOT publishable; it indexes only gas-to-solvated "
            "transfer + intra-binder conformational entropy of the displaced "
            "ligand and CANNOT be compared to experimental ΔΔG_bind."
        ),
        "lambda_densify_note": (
            "Path λ-densify spec 2026-06-05: the free leg crossover at "
            "λ 0.45→0.50 collapses ~191 kcal/mol within one λ step → adjacent "
            "states have ZERO overlap. The densified34 schedule "
            "(--free-schedule densified34) inserts a 7-window ilogistic anneal "
            "ladder; analysis is state-count-agnostic (22 or 34) + uses the "
            "DIRECTION-based leg partition for the >2-intermediate ladder."
        ),
        "notes": args.notes,
    }

    out_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    # ---- console summary ----
    def _fmt_err(v: Optional[float]) -> str:
        return f"{v:.3f}" if v is not None else "n/a"

    print()
    print("=" * 70)
    print(f"ΔΔG_int_free  =  {ddint_free:+.3f} ± {ddint_err:.3f} kcal/mol "
          f"(ranking-only; error_basis={error_basis})")
    print(f"   cp4: mean_dgbind1={cp4['mean_dgbind1_kcal']:+.3f}  "
          f"σ_btwn={_fmt_err(cp4['sigma_btwn_dgbind1_kcal'])}  "
          f"closure |dgb|={closure_cp4:.3f}  "
          f"n_rep={cp4['n_replicates']}  "
          f"states={cp4['replicates'][0]['n_states']}")
    print(f"   wt:  mean_dgbind1={wt['mean_dgbind1_kcal']:+.3f}  "
          f"σ_btwn={_fmt_err(wt['sigma_btwn_dgbind1_kcal'])}  "
          f"closure |dgb|={closure_wt:.3f}  "
          f"n_rep={wt['n_replicates']}  "
          f"states={wt['replicates'][0]['n_states']}")
    print("=" * 70)
    print(f"written: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
