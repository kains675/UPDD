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
# State-count derivation (2026-06-05 — NOT hardcoded 22).
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
    neg_pot_sink: Optional[Dict[str, Any]] = None,
) -> Tuple[float, float, float, float, int]:
    """Drop-in replacement for ``uwham.calculate_uwham`` that supports a
    schedule with more than two INTERMEDIATE states.

    Mirrors the upstream UWHAM math line-for-line (``_uwham_r``, ``_bias_fcn``,
    ``_npot_fcn`` are reused unchanged from atom_openmm.uwham) but partitions
    the legs by the DIRECTION column. Returns the same 5-tuple
    (dgb, ddgb, dgbind1, dgbind2, samples_per_replica).

    ``neg_pot_sink`` (P11 MBAR overlap QC, REPORTING-ONLY): an optional dict
    captured by reference. When provided it is populated, AFTER the leg-1 /
    leg-2 ``neg_pot`` matrices are built for the UWHAM estimator, with copies of
    those SSOT matrices + the per-target-state sample counts (keys
    ``leg1_neg_pot`` / ``leg1_N_k`` / ``leg2_neg_pot`` / ``leg2_N_k``). This is
    a PURE READ of matrices that already exist in the estimator's scope — it
    triggers NO additional ``_npot_fcn`` evaluation and changes NONE of the
    returned scalars (dgb/ddgb/dgbind1/dgbind2/samples). The downstream overlap
    QC consumes these copies; the estimator math is untouched.
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

    # P11 (REPORTING-ONLY): hand the leg-1 SSOT neg_pot + per-target-state
    # sample counts to the overlap-QC sink. Single-T (mtempt==1) so the column
    # order IS leg1stateids. N_k counts come from the state labels of data1.
    if neg_pot_sink is not None:
        _counts1 = data1["stateid"].value_counts()
        neg_pot_sink["leg1_neg_pot"] = neg_pot.copy()
        neg_pot_sink["leg1_N_k"] = [int(_counts1.get(int(s), 0))
                                    for s in leg1stateids]

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

    # P11 (REPORTING-ONLY): leg-2 SSOT neg_pot + counts. Column order IS
    # leg2stateids (reversed-state order: leg2istate..nstates-1 then [::-1]);
    # adjacency in this column order is the backward-leg λ adjacency.
    if neg_pot_sink is not None:
        _counts2 = data1["stateid"].value_counts()
        neg_pot_sink["leg2_neg_pot"] = neg_pot.copy()
        neg_pot_sink["leg2_N_k"] = [int(_counts2.get(int(s), 0))
                                    for s in leg2stateids]

    dgb = dgbind1 - dgbind2
    ddgb = np.sqrt(ddgbind2 * ddgbind2 + ddgbind1 * ddgbind1)
    return (float(dgb[0]), float(ddgb[0]), float(dgbind1[0]),
            float(dgbind2[0]), int(samplesperreplica))


# ---------------------------------------------------------------------------
# P11 — MBAR overlap-matrix QC (REPORTING-ONLY; numbers-neutral).
#
# Computes pymbar.MBAR(...).compute_overlap() (Shirts & Chodera 2008, DOI
# 10.1063/1.2978177) on the SAME per-target-state reduced-potential matrix the
# UWHAM estimator already built (the captured neg_pot, NOT the linear-λ
# surrogate in utils/adaptive_lambda/overlap.py). Reports per-adjacent
# O_{i,i+1} + the matrix min off-diagonal + a WARNING string listing adjacent
# pairs below the FLAG threshold.
#
# HARD SCOPE LOCK (overlap-QC reporting-only policy):
#   * This NEVER rejects / recomputes / rescales / skips any ΔG, error bar,
#     ze/ve/dgbind1/dgbind2, σ_btwn, bootstrap CI, or paired SEM. It GATES
#     NOTHING. It only adds a descriptive dict key.
#   * The overlap value O is NEVER folded into any uncertainty.
#   * Threshold O < 0.03 is a FLAG (report), not a gate. 0.10 would false-alarm
#     dense-by-design ladders; 0.03 flags genuine single-step collapse (e.g.
#     the λ0.45→0.50 ~191 kcal/mol zero-overlap class) only.
#   * Histogram-intersection is BANNED — overlap is MBAR W_nk·W_nkᵀ ONLY.
#   * Graceful degrade: pymbar / atom_openmm may be absent (qmmm env) — any
#     failure yields {"status": "unavailable", "reason": ...} and NOTHING
#     raises, so the existing numbers + the qmmm pytest path are unaffected.
# DOI: Klimovich, Shirts, Mobley 2015, J. Comput. Aided Mol. Des. 29(5):397-411
#      (DOI 10.1007/s10822-015-9840-9) — MBAR overlap matrix as the convergence
#      diagnostic; O_{i,i+1} >~ 0.03 "estimable", >~ 0.10 "well-determined".
# ---------------------------------------------------------------------------
OVERLAP_FLAG_THRESHOLD = 0.03      # below this an adjacent pair is FLAGGED
OVERLAP_WELL_DETERMINED = 0.10     # informational tier label only (not a flag)


# ---------------------------------------------------------------------------
# P5 — per-replicate cycle-closure DISTRUST flag (REPORTING-ONLY; ASYMMETRIC).
#
# Surfaces the EXISTING per-leg closure |dgb| = |dgbind1 - dgbind2| (already
# in result["dgb_kcal"]; forward-minus-backward UWHAM-leg disagreement of ONE
# per-direction system) as a descriptive DISTRUST label. It reads a number
# that already exists; it computes no new estimator value and changes no ΔG.
#
# ASYMMETRIC BY CONSTRUCTION (the defining property):
#   * A large |closure| sets distrust=True (a human-review prompt).
#   * A small |closure| sets distrust=False, which asserts NOTHING about
#     correctness — small closure is necessary-not-sufficient (the
#     state-0-collapse legs of 2026-06-03 had clean cheap-checks too). There
#     is NO "trust" / "pass" / "green" key anywhere in this report.
#
# HARD SCOPE LOCK (closure-distrust reporting-only policy):
#   * MUST NOT BAR/MBAR-combine dplus and dminus — they are NON-mirror
#     per-direction systems ~51 waters apart on different base ensembles;
#     combining injects a water-count/ensemble systematic. (The bidirectional
#     dgbind1/dgbind2 closure WITHIN one per-direction system is legitimate and
#     is exactly what dgb already measures — that is what this labels.)
#   * MUST NOT inject |closure|/2 (or any |dgb| function) into ddint_err or any
#     uncertainty. Closure is a SEPARATE distrust axis (a SYSTEMATIC) reported
#     alongside σ_btwn / paired-SEM (RANDOM) — never folded into σ.
#   * MUST NOT assert TRUST. distrust=False == "no closure red-flag", never
#     "this leg is correct".
#   * MUST NOT alter / re-emit the frozen free-leg numbers — purely additive.
#
# CLOSURE-ASYMMETRY CAVEAT (advisory-for-human-review, NOT a non-convergence verdict):
#   For densified38v4 dminus legs the closure asymmetry is PARTLY STRUCTURAL
#   non-mirror hysteresis (the 25 Å ATM DISPLACEMENT + the intended dminus
#   W0-graded micro-bridges, which are ΔG-UNBIASED) — NOT necessarily
#   non-convergence. A flagged dminus leg therefore reads "review: may be
#   non-mirror hysteresis", not "broken". This flag is a DISTRUST PROMPT, not a
#   pass/fail gate.
#
# Threshold |closure| > 2.0 kcal/mol = distrust. Calibrated to precedent: the
# v2_1 free-leg soft band (closure_wt=1.43 mild / closure_cp4=0.37 clean) is
# PRESERVED below the 2.0 floor (NOT re-classified as distrust), while the
# current campaign (closure_cp4=2.727 / closure_wt=3.424) DOES flag — honest,
# because that free-leg hysteresis is real and the user is owed the signal.
# Graded tier (label, not a gate): clean (<1.0) / mild (1.0-2.0, = v2_1 soft
# band) / distrust (>2.0).
# DOI: Pohorille, Jarzynski, Chipot 2010, J. Phys. Chem. B 114(32):10235-10253
#      (DOI 10.1021/jp102971x) — directional (forward/backward) disagreement as
#      a bias/convergence signature; Klimovich, Shirts, Mobley 2015 (DOI
#      10.1007/s10822-015-9840-9) closure interpretation.
# ---------------------------------------------------------------------------
CLOSURE_DISTRUST_THRESHOLD = 2.0   # |closure| above this -> distrust=True
CLOSURE_MILD_BAND = 1.0            # 1.0-2.0 mild (v2_1 soft-flag band) label

# ---------------------------------------------------------------------------
# A6: calibration-anchor ADVISORY sign-check (REPORTING-ONLY; numbers-neutral).
# The active per-direction convention is ddint = mean(cp4.dgbind1) -
# mean(wt.dgbind1); a Cp4-favorable result (Cp4 binds TIGHTER than WT, i.e.
# ΔΔG_bind < 0) therefore maps to a NEGATIVE ddint. The anchor is the
# experimental sign of ΔΔG_bind for the anti-C3b compstatin Cp4-vs-WT pair.
# R-11 / R-18: this is a SIGN advisory ONLY — magnitude is NOT calibrated, and
# the FREE-leg ddint alone is NOT ΔΔG_bind (the bound leg is required), so a
# free-leg sign check is PROVISIONAL. Anchor: Katragadda & Lambris 2006
# (DOI 10.1021/jm0603419, K_d ~15 nM); ITC/SPR SSOT: Magotti 2009
# (DOI 10.1002/jmr.972, ΔΔG_bind ≈ -1.4 ITC / -3.0 SPR kcal/mol, Cp4 better).
# ---------------------------------------------------------------------------
# sign_status (from paired_difference_stats) -> favorable-direction mapping for
# the favorable_cp4 anchor: ddint < 0 (negative) is CONSISTENT with Cp4 binding
# tighter; ddint > 0 (positive) is DISCORDANT; CI-spans-0 is UNDETERMINED.
ANCHOR_FAVORABLE_CP4_CONSISTENT_SIGN = "negative"

# ---------------------------------------------------------------------------
# A7: ADVISORY extend-recommended verdict (REPORTING-ONLY; numbers-neutral).
# Encodes the human pre-registration C4 rule: an endpoint inter-replicate
# σ_btwn (std of dgbind1) above ~2× the free paired σ (≈ 2.2 kcal/mol) is the
# stop-and-extend signal (n -> n+1 replicates) — the spread is too large for
# the current replicate count to resolve the difference. This ADVISORY merely
# compares the existing σ_btwn (already in analyze_replicate_set output) to the
# target and emits advice; it does NOT change n, does NOT gate, and is NEVER
# folded into any error bar. Convergence-diagnostic lineage: Klimovich, Shirts,
# Mobley 2015 (DOI 10.1007/s10822-015-9840-9).
# ---------------------------------------------------------------------------
EXTEND_SIGMA_BTWN_TARGET_KCAL = 2.2   # pre-reg C4: ~2x free paired σ threshold


def _closure_distrust_from_dgb(dgb, threshold=CLOSURE_DISTRUST_THRESHOLD):
    """Reduce an already-computed leg closure ``dgb`` to the DISTRUST report.

    Pure / numeric (no estimator re-eval). ``dgb`` is the leg-internal
    ``dgbind1 - dgbind2`` (forward minus backward); its absolute value is the
    cycle-closure disagreement. ASYMMETRIC: ``distrust=True`` iff
    ``|dgb| > threshold``; ``distrust=False`` is "no closure red-flag" and
    asserts NOTHING about correctness. NEVER emits a trust / pass / green key.
    """
    closure_abs = abs(float(dgb))
    distrust = closure_abs > float(threshold)
    if closure_abs < CLOSURE_MILD_BAND:
        tier = "clean"
    elif closure_abs <= float(threshold):
        tier = "mild"
    else:
        tier = "distrust"
    warning = None
    if distrust:
        warning = (
            f"CLOSURE_DISTRUST: |closure| = {closure_abs:.3f} kcal/mol > "
            f"{float(threshold):.1f} (forward/backward UWHAM-leg disagreement). "
            f"ADVISORY distrust-only — review this leg's directional "
            f"consistency. For a densified38v4 dminus leg this may be NON-MIRROR "
            f"hysteresis (25 Å DISPLACEMENT + ΔG-unbiased W0 micro-bridges), "
            f"NOT necessarily non-convergence. Does NOT gate / rescale ΔG and is "
            f"NEVER folded into any error bar."
        )
    return {
        "status": "ok",
        "basis": "abs_dgbind1_minus_dgbind2",
        "threshold": float(threshold),
        "closure_kcal_abs": closure_abs,
        "distrust": bool(distrust),   # ASYMMETRIC: True flags; False is NOT "trust"
        "tier": tier,
        "warning": warning,
        "note": (
            "ADVISORY distrust-only; False == no closure red-flag (NOT a trust "
            "assertion). Non-mirror per-direction hysteresis is partly "
            "structural, not necessarily non-convergence — review, don't reject."
        ),
    }


def compute_closure_distrust(result, threshold=CLOSURE_DISTRUST_THRESHOLD):
    """Build the per-leg closure-distrust report from a leg ``result`` dict.

    Reads the EXISTING ``result["dgb_kcal"]`` (leg-internal closure) and labels
    it. Graceful: if ``dgb_kcal`` is absent / non-finite (e.g. a leg whose
    backward dgbind2 was unavailable), returns ``{"status":
    "closure_unavailable", ...}`` and NEVER raises — so the existing numbers
    are unaffected. REPORTING-ONLY, ASYMMETRIC (scope-lock policy): no ΔG /
    error bar is rejected, recomputed, or rescaled; no trust is asserted.
    """
    import math

    if not isinstance(result, dict) or "dgb_kcal" not in result:
        return {
            "status": "closure_unavailable",
            "reason": "no dgb_kcal on leg result (dgbind2/closure unavailable)",
            "basis": "abs_dgbind1_minus_dgbind2",
            "threshold": float(threshold),
        }
    dgb = result["dgb_kcal"]
    try:
        dgb_f = float(dgb)
    except (TypeError, ValueError):
        return {
            "status": "closure_unavailable",
            "reason": f"dgb_kcal not a finite number ({dgb!r})",
            "basis": "abs_dgbind1_minus_dgbind2",
            "threshold": float(threshold),
        }
    if not math.isfinite(dgb_f):
        return {
            "status": "closure_unavailable",
            "reason": f"dgb_kcal non-finite ({dgb_f})",
            "basis": "abs_dgbind1_minus_dgbind2",
            "threshold": float(threshold),
        }
    return _closure_distrust_from_dgb(dgb_f, threshold=threshold)


def _compute_mbar_overlap_matrix(neg_pot, n_k):
    """Run ``pymbar.MBAR(...).compute_overlap()`` on a captured neg_pot matrix.

    ``neg_pot`` is the UWHAM ``[N_samples, K_states]`` NEGATIVE reduced-potential
    matrix (``-β·U_k(x_n)``; the ``_npot_fcn`` output ``_uwham_r`` consumes as
    ``logQ``). pymbar's ``u_kn`` convention is the POSITIVE reduced potential
    ``u_kn[k, n] = β·U_k(x_n)`` shaped ``[K_states, N_samples]`` — so
    ``u_kn = -neg_pot.T``. ``n_k`` is the per-state sample count (len == K,
    sum == N).

    Returns the K×K overlap matrix (numpy array) or raises — the caller wraps
    this in the graceful-degrade try/except so an absent pymbar never crashes.
    """
    import numpy as np
    from pymbar import MBAR  # local import — graceful-degrade if absent

    arr = np.asarray(neg_pot, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"neg_pot must be 2-D, got shape {arr.shape}")
    n_samples, k_states = arr.shape
    nk = np.asarray(n_k, dtype=int)
    if nk.shape != (k_states,):
        raise ValueError(
            f"n_k length {nk.shape} != K_states {k_states} (overlap matrix is "
            f"computed over the per-target-state columns)"
        )
    if int(nk.sum()) != n_samples:
        raise ValueError(
            f"sum(n_k)={int(nk.sum())} != N_samples={n_samples} — the captured "
            f"sample counts must partition the neg_pot rows"
        )
    u_kn = -arr.T  # negative reduced potential -> positive reduced potential
    mbar = MBAR(u_kn, nk)
    overlap = mbar.compute_overlap()
    matrix = np.asarray(overlap["matrix"], dtype=float)
    return matrix


def _overlap_qc_from_matrix(matrix, threshold=OVERLAP_FLAG_THRESHOLD):
    """Reduce a K×K MBAR overlap matrix to the per-leg QC report.

    Pure / numeric: reads adjacent off-diagonals ``O_{i,i+1}``, the matrix min
    off-diagonal, flags adjacent pairs below ``threshold`` (FLAG, not gate),
    and assigns an informational tier label. No estimator value is touched.
    """
    import numpy as np

    m = np.asarray(matrix, dtype=float)
    k = m.shape[0]
    adjacent_O = [float(m[i, i + 1]) for i in range(k - 1)]
    # Min off-diagonal across the whole matrix (diagonal excluded).
    if k >= 2:
        off = m.copy()
        np.fill_diagonal(off, np.inf)
        min_off_diagonal = float(off.min())
        min_adjacent_O = float(min(adjacent_O)) if adjacent_O else None
    else:
        min_off_diagonal = None
        min_adjacent_O = None

    flagged_pairs = [
        {"i": i, "j": i + 1, "O": adjacent_O[i]}
        for i in range(len(adjacent_O))
        if adjacent_O[i] < threshold
    ]
    if min_adjacent_O is None:
        tier = "single_state"
    elif min_adjacent_O < threshold:
        tier = "COLLAPSE"
    elif min_adjacent_O < OVERLAP_WELL_DETERMINED:
        tier = "MARGINAL"
    else:
        tier = "WELL_DETERMINED"

    if flagged_pairs:
        pairs_str = ", ".join(
            f"O_{p['i']},{p['j']}={p['O']:.4f}" for p in flagged_pairs
        )
        warning = (
            f"LOW_OVERLAP: {len(flagged_pairs)} adjacent pair(s) below "
            f"{threshold:.3f} MBAR overlap (single-step phase-space collapse): "
            f"{pairs_str}. Reporting-only (does NOT gate / rescale ΔG)."
        )
    else:
        warning = None

    return {
        "status": "ok",
        "basis": "mbar_compute_overlap",
        "threshold": float(threshold),
        "n_states": int(k),
        "adjacent_O": adjacent_O,
        "min_adjacent_O": min_adjacent_O,
        "min_off_diagonal_O": min_off_diagonal,
        "flagged_pairs": flagged_pairs,
        "tier": tier,
        "warning": warning,
    }


def compute_overlap_qc(neg_pot_sink, threshold=OVERLAP_FLAG_THRESHOLD):
    """Build the per-leg overlap-QC report from a populated ``neg_pot_sink``.

    Wraps the MBAR overlap computation in a graceful-degrade try/except: if
    pymbar is absent, the matrix is unavailable, or anything raises, returns
    ``{"status": "unavailable", "reason": ...}`` and NEVER propagates — so the
    qmmm-env path (no pymbar/atom_openmm guarantee) and every existing number
    are unaffected. REPORTING-ONLY (scope-lock policy): no ΔG / error bar is
    rejected, recomputed, or rescaled here.

    Reports the FORWARD (leg-1) overlap as the primary block plus, when present,
    the backward (leg-2) overlap. ``min_adjacent_O`` aggregates both legs.
    """
    if not neg_pot_sink or "leg1_neg_pot" not in neg_pot_sink:
        return {
            "status": "unavailable",
            "reason": (
                "no captured neg_pot matrix (canonical-22 upstream path does "
                "not expose it, or capture was disabled)"
            ),
            "basis": "mbar_compute_overlap",
            "threshold": float(threshold),
        }

    legs: Dict[str, Any] = {}
    mins: List[float] = []
    flagged: List[Dict[str, Any]] = []
    try:
        for tag, mat_key, nk_key in (
            ("forward", "leg1_neg_pot", "leg1_N_k"),
            ("backward", "leg2_neg_pot", "leg2_N_k"),
        ):
            if mat_key not in neg_pot_sink:
                continue
            matrix = _compute_mbar_overlap_matrix(
                neg_pot_sink[mat_key], neg_pot_sink[nk_key]
            )
            leg_qc = _overlap_qc_from_matrix(matrix, threshold=threshold)
            legs[tag] = leg_qc
            if leg_qc["min_adjacent_O"] is not None:
                mins.append(leg_qc["min_adjacent_O"])
            for p in leg_qc["flagged_pairs"]:
                flagged.append({"leg": tag, **p})
    except Exception as exc:  # noqa: BLE001 — graceful degrade (any failure)
        return {
            "status": "unavailable",
            "reason": f"{type(exc).__name__}: {exc}",
            "basis": "mbar_compute_overlap",
            "threshold": float(threshold),
        }

    if not legs:
        return {
            "status": "unavailable",
            "reason": "no overlap matrix could be computed for any leg",
            "basis": "mbar_compute_overlap",
            "threshold": float(threshold),
        }

    min_adjacent_O = float(min(mins)) if mins else None
    if min_adjacent_O is None:
        tier = "single_state"
    elif min_adjacent_O < threshold:
        tier = "COLLAPSE"
    elif min_adjacent_O < OVERLAP_WELL_DETERMINED:
        tier = "MARGINAL"
    else:
        tier = "WELL_DETERMINED"

    if flagged:
        pairs_str = ", ".join(
            f"{p['leg']} O_{p['i']},{p['j']}={p['O']:.4f}" for p in flagged
        )
        warning = (
            f"LOW_OVERLAP: {len(flagged)} adjacent pair(s) below "
            f"{threshold:.3f} MBAR overlap: {pairs_str}. Reporting-only "
            f"(does NOT gate / rescale ΔG)."
        )
    else:
        warning = None

    return {
        "status": "ok",
        "basis": "mbar_compute_overlap",
        "threshold": float(threshold),
        "min_adjacent_O": min_adjacent_O,
        "tier": tier,
        "flagged_pairs": flagged,
        "warning": warning,
        "legs": legs,
    }


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
    neg_pot_sink: Optional[Dict[str, Any]] = None,
) -> Tuple[float, float, float, float, int]:
    """Dispatch to upstream calculate_uwham (canonical) or the vendored
    multi-intermediate variant (densified) for one analysis window.

    ``neg_pot_sink`` (P11 overlap QC, REPORTING-ONLY): when provided AND the
    densified (>2 INTERMEDIATE) path is taken, it is forwarded to
    ``_calculate_uwham_multi_intermediate`` to capture the SSOT neg_pot
    matrices. The canonical ≤2-INTERMEDIATE upstream path does NOT expose
    neg_pot, so the sink stays empty there and the QC degrades to
    ``status:"unavailable"``. The block-bootstrap callers never pass a sink
    (no capture during the per-block windows).

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
            neg_pot_sink=neg_pot_sink,
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

    # P11 overlap-QC sink (REPORTING-ONLY): only the strict (C.1) window path
    # populates it (the densified multi-intermediate variant exposes neg_pot);
    # the legacy upstream path leaves it empty → QC degrades to "unavailable".
    neg_pot_sink: Optional[Dict[str, Any]] = (
        {} if analyze_fn is _leg_analysis_window else None
    )
    analyze_kwargs: Dict[str, Any] = dict(
        leg_dir=leg_dir, jobname=jobname, schedule=schedule,
        mintimeid=mintimeid, maxtimeid=maxtimeid, n_states=expected_replicas,
    )
    if neg_pot_sink is not None:
        analyze_kwargs["neg_pot_sink"] = neg_pot_sink

    t0 = time.time()
    dgb, ddgb, dgbind1, dgbind2, samples = analyze_fn(**analyze_kwargs)
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

    # P11 MBAR overlap QC (REPORTING-ONLY; graceful-degrade — never raises).
    # Additive dict key; does NOT touch any ΔG / error bar above.
    result["overlap_qc"] = compute_overlap_qc(neg_pot_sink)

    # P5 per-replicate cycle-closure DISTRUST flag (REPORTING-ONLY; ASYMMETRIC;
    # graceful — never raises). Labels the EXISTING result["dgb_kcal"]; adds NO
    # estimator value and changes NO ΔG / error bar. distrust=True is a
    # human-review prompt, distrust=False is NOT a trust assertion.
    result["closure_distrust"] = compute_closure_distrust(result)

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


def paired_difference_stats(
    a_vals: List[float],
    b_vals: List[float],
) -> Dict[str, Any]:
    """Paired-difference statistics over MATCHED-seed replicates.

    Reusable helper for any cp4-vs-wt (or, later, bound-vs-free) ΔΔG combiner
    that pairs replicate i of one endpoint with replicate i of the other (same
    independent seed). Honest paired statistics: the naive independent
    quadrature of the two endpoints' STDEVs ignores the cross-leg correlation
    and conflates STDEV with SEM. With matched seeds, the per-replicate
    difference ΔΔG_i = a_i − b_i has a paired SD whose inter-replicate SEM is
    the correct uncertainty on the mean ΔΔG.

    Inputs ``a_vals`` / ``b_vals`` MUST be index-aligned (a_vals[i] and
    b_vals[i] are the SAME seed) and equal length. Returns:

      * ``n``                — number of matched pairs
      * ``per_rep_ddg``      — [a_i − b_i] (ranking vector)
      * ``mean``             — mean ΔΔG (the POINT ESTIMATE; identical to
                               mean(a) − mean(b))
      * ``paired_sd``        — sample SD of per_rep_ddg (ddof=1)
      * ``sem``              — paired SEM = paired_sd / sqrt(n)
      * ``pearson_r``        — cross-leg Pearson correlation of a vs b
      * ``df``               — degrees of freedom (n − 1)
      * ``t_stat``           — paired t = mean / sem
      * ``t_crit_95``        — two-sided 0.975 t critical value at df
      * ``ci95``             — [mean − t_crit·sem, mean + t_crit·sem]
      * ``sign_status``      — 'undetermined' if ci95 spans 0, else
                               'positive' / 'negative'
      * ``error_basis``      — 'paired_sigma_btwn_t' (provenance tag; the
                               primary error is the paired SEM, NOT a STDEV)

    Pure / read-only — no estimator re-evaluation; operates on already-computed
    per-replicate dgbind1 values. With n < 2 the dispersion fields are None
    (no paired SD / SEM / t available from a single pair).
    """
    import numpy as np

    if len(a_vals) != len(b_vals):
        raise ValueError(
            f"paired_difference_stats: matched-seed vectors must be equal "
            f"length (got len(a)={len(a_vals)} != len(b)={len(b_vals)}). "
            f"Replicate sets must be index-aligned by seed."
        )
    a = np.asarray(a_vals, dtype=float)
    b = np.asarray(b_vals, dtype=float)
    diff = a - b
    n = int(diff.shape[0])
    mean = float(diff.mean()) if n >= 1 else None

    result: Dict[str, Any] = {
        "n": n,
        "per_rep_ddg": [float(x) for x in diff],
        "mean": mean,
        "paired_sd": None,
        "sem": None,
        "pearson_r": None,
        "df": (n - 1) if n >= 1 else None,
        "t_stat": None,
        "t_crit_95": None,
        "ci95": None,
        "sign_status": "undetermined",
        "error_basis": "paired_sigma_btwn_t",
    }
    if n < 2:
        return result

    paired_sd = float(diff.std(ddof=1))
    sem = float(paired_sd / (n ** 0.5))
    # Cross-leg Pearson r (guard against a degenerate constant vector).
    if a.std() > 0 and b.std() > 0:
        pearson_r = float(np.corrcoef(a, b)[0, 1])
    else:
        pearson_r = None
    df = n - 1
    t_crit = _student_t_ppf_975(df)
    t_stat = float(mean / sem) if sem > 0 else None
    if t_crit is not None:
        lo = float(mean - t_crit * sem)
        hi = float(mean + t_crit * sem)
        ci95 = [lo, hi]
        if lo > 0:
            sign_status = "positive"
        elif hi < 0:
            sign_status = "negative"
        else:
            sign_status = "undetermined"
    else:
        ci95 = None
        sign_status = "undetermined"

    result.update({
        "paired_sd": paired_sd,
        "sem": sem,
        "pearson_r": pearson_r,
        "df": df,
        "t_stat": t_stat,
        "t_crit_95": t_crit,
        "ci95": ci95,
        "sign_status": sign_status,
    })
    return result


def _student_t_ppf_975(df: int) -> Optional[float]:
    """Two-sided 0.975 Student-t critical value at ``df`` degrees of freedom.

    Prefers scipy (exact); falls back to a small hard-coded table for the
    low-df regime relevant to a handful of replicates (so the helper has no
    hard scipy dependency). Returns None when df < 1.
    """
    if df < 1:
        return None
    try:
        from scipy import stats  # type: ignore
        return float(stats.t.ppf(0.975, df))
    except Exception:
        # df -> t_{0.975} table (two-sided 95%). Covers the realistic
        # replicate-count range; beyond the table fall back to the z value.
        table = {
            1: 12.7062, 2: 4.3027, 3: 3.1824, 4: 2.7764, 5: 2.5706,
            6: 2.4469, 7: 2.3646, 8: 2.3060, 9: 2.2622, 10: 2.2281,
            11: 2.2010, 12: 2.1788, 15: 2.1314, 20: 2.0860, 30: 2.0423,
        }
        if df in table:
            return table[df]
        if df > 30:
            return 1.9600  # z_{0.975} asymptote
        # Nearest tabulated df below (conservative-ish for the gap rows).
        keys = sorted(k for k in table if k <= df)
        return table[keys[-1]] if keys else None


def compute_sign_vs_anchor(
    paired: Optional[Dict[str, Any]],
    ddint_point: Optional[float] = None,
    expected_sign: str = "favorable_cp4",
    leg: str = "free",
) -> Dict[str, Any]:
    """ADVISORY comparison of the computed ΔΔG sign vs the experimental anchor.

    REPORTING-ONLY (R-11 / R-18): compares the SIGN of the active-convention
    ΔΔG (ddint = mean(cp4.dgbind1) - mean(wt.dgbind1)) to the expected
    favorable direction from the calibration anchor. Emits one of:

      * ``consistent``    — the sign agrees with the favorable direction AND the
                            paired 95% CI excludes 0 (sign is statistically
                            resolved).
      * ``discordant``    — the sign DISAGREES with the favorable direction with
                            the CI excluding 0 (a resolved opposite sign).
      * ``undetermined``  — the paired 95% CI spans 0 (sign not resolved), or no
                            paired statistics are available.

    Does NOT gate, does NOT change any number, NEVER makes a magnitude /
    calibration claim. Crucially, for ``leg == "free"`` the ddint is the
    FREE-leg ΔΔG_int only and is NOT ΔΔG_bind (the bound leg is required for a
    true binding free energy), so the result carries a ``provisional`` flag and
    a caveat: a free-leg sign that is discordant with the binding anchor is NOT
    by itself evidence the design ranks wrong — only the full bound-minus-free
    cycle sign can be compared to the binding anchor.

    Inputs: ``paired`` is the paired_difference_stats dict (or None);
    ``ddint_point`` is the point estimate (``paired['mean']`` when omitted).
    ``expected_sign`` is the anchor's ``expected_ddg_bind_sign`` (only
    ``favorable_cp4`` is wired; any other value -> undetermined with a reason).

    DOI 10.1021/jm0603419 (Katragadda & Lambris 2006 anchor K_d),
    10.1002/jmr.972 (Magotti 2009 ITC/SPR ΔΔG_bind SSOT, sign only).
    """
    consistent_sign = ANCHOR_FAVORABLE_CP4_CONSISTENT_SIGN  # 'negative'
    provisional = (leg == "free")
    base: Dict[str, Any] = {
        "basis": "sign_of_ddint_vs_favorable_direction",
        "expected_ddg_bind_sign": expected_sign,
        "convention": "ddint = mean(cp4.dgbind1) - mean(wt.dgbind1)",
        "favorable_ddint_sign": consistent_sign,
        "leg": leg,
        "provisional": provisional,
        "sign_vs_anchor": "undetermined",
        "computed_sign_status": None,
        "ddint_point_estimate_kcal": (
            float(ddint_point) if ddint_point is not None
            else (float(paired["mean"]) if isinstance(paired, dict)
                  and paired.get("mean") is not None else None)
        ),
        "warning": None,
        "note": (
            "ADVISORY sign-only comparison (R-11 / R-18): magnitude is NOT "
            "calibrated and is NEVER claimed. FREE-leg ddint is NOT ΔΔG_bind — "
            "this free-leg sign check is PROVISIONAL; only the full "
            "bound-minus-free cycle sign may be compared to the binding anchor."
        ),
    }

    if expected_sign != "favorable_cp4":
        base["sign_vs_anchor"] = "undetermined"
        base["reason"] = (
            f"unsupported expected_ddg_bind_sign={expected_sign!r}; only "
            f"'favorable_cp4' is wired for the sign mapping."
        )
        return base

    sign_status = (
        paired.get("sign_status") if isinstance(paired, dict) else None
    )
    base["computed_sign_status"] = sign_status

    if sign_status in (None, "undetermined"):
        # CI spans 0 (or no paired stats) → sign not statistically resolved.
        base["sign_vs_anchor"] = "undetermined"
        base["reason"] = (
            "paired 95% CI spans 0 (sign unresolved)" if sign_status
            == "undetermined" else "no paired statistics available"
        )
        return base

    if sign_status == consistent_sign:
        base["sign_vs_anchor"] = "consistent"
    else:
        base["sign_vs_anchor"] = "discordant"
        base["warning"] = (
            f"SIGN_VS_ANCHOR DISCORDANT: computed ddint sign is "
            f"'{sign_status}' but the favorable direction for the anti-C3b "
            f"compstatin Cp4-vs-WT anchor (Katragadda & Lambris 2006, "
            f"DOI 10.1021/jm0603419; Magotti 2009 ITC/SPR, DOI 10.1002/jmr.972) "
            f"is '{consistent_sign}' (Cp4 binds tighter). ADVISORY sign-only — "
            f"does NOT gate / rescale any number. "
            + (
                "NOTE this is the FREE leg alone, NOT ΔΔG_bind: a free-leg "
                "sign discordance is NOT by itself a ranking failure; the "
                "bound leg is required before comparing to the binding anchor."
                if provisional else
                "Compared against the binding anchor for the full cycle."
            )
        )
    return base


def compute_extend_recommended(
    cp4: Dict[str, Any],
    wt: Dict[str, Any],
    target_kcal: float = EXTEND_SIGMA_BTWN_TARGET_KCAL,
) -> Dict[str, Any]:
    """ADVISORY extend-recommended verdict from the pre-registration C4 rule.

    REPORTING-ONLY (numbers-neutral): reads the EXISTING per-endpoint
    inter-replicate ``sigma_btwn_dgbind1_kcal`` (already computed by
    ``analyze_replicate_set``) and recommends extending the replicate count
    (n -> n+1) when the WORST endpoint σ_btwn exceeds ``target_kcal`` (the
    pre-reg C4 stop-and-extend threshold, ≈ 2× the free paired σ ≈ 2.2
    kcal/mol). Does NOT change n, does NOT gate, is NEVER folded into an error
    bar — it only emits advice for the operator to act on.

    Returns a dict with ``recommended`` (bool), ``current_sigma_btwn_kcal``
    (the max of the two endpoint σ_btwn — the binding constraint), ``target_kcal``,
    ``n_current`` (replicate count), ``reason``, and per-endpoint σ_btwn for
    transparency. When σ_btwn is unavailable (single replicate) the verdict is
    ``recommended: False`` with an explanatory reason (cannot assess spread from
    one replicate — a separate single-run concern, not a C4 extend trigger).

    DOI 10.1007/s10822-015-9840-9 (Klimovich/Shirts/Mobley 2015 convergence
    diagnostics).
    """
    cp4_sigma = cp4.get("sigma_btwn_dgbind1_kcal")
    wt_sigma = wt.get("sigma_btwn_dgbind1_kcal")
    # n is the matched replicate count (min of the two; they are normally equal).
    n_cp4 = cp4.get("n_replicates")
    n_wt = wt.get("n_replicates")
    n_current = (
        min(n_cp4, n_wt) if isinstance(n_cp4, int) and isinstance(n_wt, int)
        else (n_cp4 if isinstance(n_cp4, int) else n_wt)
    )

    sigmas = [s for s in (cp4_sigma, wt_sigma) if s is not None]
    if not sigmas:
        return {
            "basis": "max_endpoint_sigma_btwn_vs_pre_reg_C4_target",
            "recommended": False,
            "current_sigma_btwn_kcal": None,
            "target_kcal": float(target_kcal),
            "n_current": n_current,
            "per_endpoint_sigma_btwn_kcal": {
                "cp4": cp4_sigma, "wt": wt_sigma,
            },
            "reason": (
                "σ_btwn unavailable (single replicate per endpoint) — cannot "
                "assess inter-replicate spread; the pre-reg C4 extend trigger "
                "needs >=2 replicates. Not a C4 recommendation either way."
            ),
            "note": (
                "ADVISORY only — does NOT change n, does NOT gate, NEVER folded "
                "into any error bar."
            ),
        }

    worst = float(max(sigmas))
    recommended = worst > float(target_kcal)
    if recommended:
        reason = (
            f"σ_btwn={worst:.3f} > target {float(target_kcal):.2f} kcal/mol "
            f"(pre-reg C4 ≈ 2× free paired σ): inter-replicate spread is too "
            f"large for n={n_current} to resolve the difference — RECOMMEND "
            f"extending replicates (n -> n+1) and re-assessing."
        )
    else:
        reason = (
            f"σ_btwn={worst:.3f} <= target {float(target_kcal):.2f} kcal/mol: "
            f"inter-replicate spread within the pre-reg C4 budget at "
            f"n={n_current}; no extension recommended on this criterion."
        )
    return {
        "basis": "max_endpoint_sigma_btwn_vs_pre_reg_C4_target",
        "recommended": bool(recommended),
        "current_sigma_btwn_kcal": worst,
        "target_kcal": float(target_kcal),
        "n_current": n_current,
        "per_endpoint_sigma_btwn_kcal": {"cp4": cp4_sigma, "wt": wt_sigma},
        "reason": reason,
        "note": (
            "ADVISORY only — does NOT change n, does NOT gate, NEVER folded "
            "into any error bar. Encodes the human pre-registration C4 "
            "stop-and-extend rule (σ_btwn > 2× free paired σ ≈ 2.2 kcal/mol)."
        ),
    }


def build_v3_payload(
    cp4: Dict[str, Any],
    wt: Dict[str, Any],
    protocol: str,
    mintimeid: Optional[int],
    maxtimeid: Optional[int],
    notes: str,
) -> Dict[str, Any]:
    """Assemble the ``trackb_ddint_free_v3`` payload from analyzed endpoint
    blocks.

    Shared by the live analysis path (``main``) and the ``--reemit-v3-from``
    path (which feeds already-computed cp4/wt blocks from an archived v2
    artifact). The POINT ESTIMATE ``ddint_free_kcal`` is mean(cp4.dgbind1) −
    mean(wt.dgbind1), unchanged from v2. The PRIMARY uncertainty is the
    paired-difference SEM over matched-seed replicates; the v2 STDEV
    quadrature is retained as a clearly-labelled legacy field (units differ:
    SEM vs STDEV — the two are NOT overloaded onto one key).
    """
    ddint_free = cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"]

    cp4_err = cp4["sigma_btwn_dgbind1_kcal"]
    wt_err = wt["sigma_btwn_dgbind1_kcal"]
    error_basis = "sigma_btwn"
    if cp4_err is None or wt_err is None:
        # Single-run fallback — analytic error of replicate 0 (diagnostic).
        cp4_err = cp4["replicates"][0]["ddgb_kcal"]
        wt_err = wt["replicates"][0]["ddgb_kcal"]
        error_basis = "single_run_analytic_diag"
    # LEGACY error: naive independent quadrature of the two endpoint STDEVs.
    # Retained ONLY as a clearly-labelled legacy field — it ignores the
    # cross-leg correlation and conflates STDEV with SEM (units!). The
    # honest paired SEM below is the primary uncertainty.
    ddint_err_quadrature_legacy = (cp4_err ** 2 + wt_err ** 2) ** 0.5

    # PRIMARY error: paired-difference statistics over matched-seed replicates
    # (paired SD / SEM / cross-leg r / t / 95% CI / sign_status). Only valid
    # when both endpoints have the SAME number of matched-seed replicates.
    cp4_dgbind1 = [r["dgbind1_kcal"] for r in cp4["replicates"]]
    wt_dgbind1 = [r["dgbind1_kcal"] for r in wt["replicates"]]
    paired: Optional[Dict[str, Any]] = None
    if len(cp4_dgbind1) == len(wt_dgbind1) and len(cp4_dgbind1) >= 1:
        paired = paired_difference_stats(cp4_dgbind1, wt_dgbind1)

    # Closure check per endpoint (replicate-0 representative).
    closure_cp4 = abs(cp4["replicates"][0]["dgb_kcal"])
    closure_wt = abs(wt["replicates"][0]["dgb_kcal"])

    # error_basis: when matched-seed paired stats are available they are the
    # PRIMARY uncertainty; otherwise keep the upstream σ_btwn / single-run tag.
    if paired is not None and paired.get("sem") is not None:
        payload_error_basis = "paired_sigma_btwn_t"
    else:
        payload_error_basis = error_basis

    # P11 MBAR overlap-QC top-level summary (REPORTING-ONLY; ADDITIVE). Each
    # endpoint's per-replicate result already carries its own ``overlap_qc``
    # (in cp4/wt → replicates[i]); this is a discoverable top-level digest of
    # the replicate-0 representative per endpoint. NEVER feeds any ΔG / error
    # bar — purely descriptive QC metadata. Absent on synthetic blocks that
    # lack a captured matrix → "unavailable".
    overlap_qc_summary = _build_overlap_qc_summary(cp4, wt)

    # P5 closure-quarantine digest (REPORTING-ONLY; ASYMMETRIC; ADDITIVE). Reads
    # the per-replicate closure_distrust already attached to each endpoint's
    # replicate-0 result; lists the flagged endpoints (|closure|>2.0). NEVER
    # feeds any ΔG / error bar; distrust=False is not a trust assertion.
    closure_quarantine = _build_closure_quarantine(cp4, wt)

    # A6 calibration-anchor ADVISORY sign-check (REPORTING-ONLY; numbers-
    # neutral). Compares the SIGN of the computed ddint (active convention) to
    # the favorable_cp4 anchor direction. FREE leg here → provisional sign check
    # (free-leg ddint is NOT ΔΔG_bind). Never gates / changes a number.
    sign_vs_anchor = compute_sign_vs_anchor(
        paired=paired, ddint_point=float(ddint_free),
        expected_sign="favorable_cp4", leg="free",
    )

    # A7 extend-recommended ADVISORY (REPORTING-ONLY; numbers-neutral). Encodes
    # the pre-reg C4 stop-and-extend rule (σ_btwn > 2× free paired σ ≈ 2.2
    # kcal/mol → recommend n -> n+1). Never changes n / gates / touches σ.
    extend_recommended = compute_extend_recommended(cp4, wt)

    return {
        "schema_version": "trackb_ddint_free_v3",
        "regime": "ranking_only_R11",
        "protocol": protocol,
        "mintimeid": mintimeid,
        "maxtimeid": maxtimeid,
        "error_basis": payload_error_basis,
        "cp4": cp4,
        "wt": wt,
        "ddint_free_kcal": float(ddint_free),
        # PRIMARY uncertainty (v3): paired-difference SEM over matched seeds.
        # This is a SEM (units: kcal/mol on the MEAN), NOT a STDEV — it is a
        # distinct field from the legacy quadrature value below; the two are
        # deliberately not overloaded onto one key (the v2 ddint_err_kcal was
        # a STDEV quadrature; silently reusing that key for a SEM would change
        # the meaning of the number without changing its name).
        "ddint_err_paired_sem_kcal": (
            float(paired["sem"]) if paired and paired.get("sem") is not None
            else None
        ),
        # LEGACY (v2): naive independent quadrature of the two endpoint STDEVs
        # — ignores cross-leg correlation (r) + conflates STDEV with SEM.
        # Retained verbatim for back-compat / provenance only.
        "ddint_err_quadrature_legacy_kcal": float(ddint_err_quadrature_legacy),
        # Full paired-difference statistics block (None when endpoints have
        # unequal replicate counts so pairing is undefined).
        "paired_difference_stats": paired,
        "closure_cp4_kcal_abs": float(closure_cp4),
        "closure_wt_kcal_abs": float(closure_wt),
        "interpretation_note": (
            "ddint_free_kcal = cp4.mean_dgbind1 - wt.mean_dgbind1 (forward "
            "+1 leg) — POINT ESTIMATE, unchanged from v2. PRIMARY error "
            "(v3) = ddint_err_paired_sem_kcal: the paired-difference SEM over "
            "MATCHED-seed replicates (paired SD / sqrt(n)), with full stats "
            "in paired_difference_stats (cross-leg Pearson r, paired t(df), "
            "t-based 95% CI, sign_status). The v2 ddint_err_quadrature_legacy_"
            "kcal is a naive independent quadrature of the two endpoint STDEVs "
            "— it IGNORES the cross-leg correlation and conflates STDEV with "
            "SEM; kept only as a labelled legacy field. sign_status is "
            "'undetermined' when the 95% CI spans 0. "
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
            "λ-densify spec 2026-06-05: the free leg crossover at "
            "λ 0.45→0.50 collapses ~191 kcal/mol within one λ step → adjacent "
            "states have ZERO overlap. The densified34 schedule "
            "(--free-schedule densified34) inserts a 7-window ilogistic anneal "
            "ladder; analysis is state-count-agnostic (22 or 34) + uses the "
            "DIRECTION-based leg partition for the >2-intermediate ladder."
        ),
        # P11 MBAR overlap-matrix QC (REPORTING-ONLY; ADDITIVE; numbers-neutral).
        # Top-level digest of the per-endpoint replicate-0 overlap report. The
        # full per-replicate reports live under cp4/wt → replicates[i] →
        # overlap_qc. Threshold O<0.03 is a FLAG, NOT a gate; O is NEVER folded
        # into any error bar (reporting-only scope-lock). DOI 10.1063/
        # 1.2978177 (MBAR), 10.1007/s10822-015-9840-9 (overlap diagnostic).
        "overlap_qc": overlap_qc_summary,
        # P5 closure-quarantine flag (REPORTING-ONLY; ASYMMETRIC; ADDITIVE;
        # numbers-neutral). per_rep_closure + flagged endpoints where
        # |closure|>2.0; distrust-only (NEVER asserts trust). The closure
        # numbers themselves are byte-identical to closure_cp4/wt_kcal_abs —
        # this only LABELS them. Non-mirror per-direction hysteresis is partly
        # structural, not necessarily non-convergence (reporting-only
        # scope-lock). DOI 10.1021/jp102971x (directional bias), 10.1007/
        # s10822-015-9840-9 (closure diagnostic).
        "closure_quarantine": closure_quarantine,
        # A6 calibration-anchor ADVISORY sign-check (REPORTING-ONLY; numbers-
        # neutral). sign_vs_anchor ∈ {consistent, discordant, undetermined}
        # comparing the computed ddint sign (active convention) to the
        # favorable_cp4 anchor. PROVISIONAL on the free leg (ddint is NOT
        # ΔΔG_bind). Sign-only, R-11/R-18 — never a magnitude / calibration
        # claim, never gates. Anchor DOI 10.1021/jm0603419 (Katragadda &
        # Lambris 2006), ITC/SPR SSOT DOI 10.1002/jmr.972 (Magotti 2009).
        "sign_vs_anchor": sign_vs_anchor,
        # A7 extend-recommended ADVISORY (REPORTING-ONLY; numbers-neutral).
        # recommended=True when the worst endpoint σ_btwn exceeds the pre-reg C4
        # target (≈ 2× free paired σ ≈ 2.2 kcal/mol). Advises n -> n+1; NEVER
        # changes n / gates / folds into an error bar. DOI 10.1007/
        # s10822-015-9840-9 (convergence diagnostics).
        "extend_recommended": extend_recommended,
        "notes": notes,
    }


def _build_overlap_qc_summary(
    cp4: Dict[str, Any],
    wt: Dict[str, Any],
) -> Dict[str, Any]:
    """Assemble the top-level overlap-QC digest from each endpoint's
    replicate-0 ``overlap_qc`` report (REPORTING-ONLY; numbers-neutral).

    Reads (does not compute) the per-replicate ``overlap_qc`` already attached
    to each endpoint by ``analyze_one_leg``. When an endpoint block has no
    captured report (e.g. a synthetic test block, or the canonical-22 upstream
    path that does not expose neg_pot), that endpoint's digest is
    ``{"status": "unavailable"}``. The combined ``min_adjacent_O`` is the min
    over whichever endpoints reported a value; ``any_flagged`` is True iff any
    endpoint flagged an adjacent pair below its threshold. This digest NEVER
    feeds a ΔG or error bar.
    """
    def _ep_report(block: Dict[str, Any]) -> Dict[str, Any]:
        reps = block.get("replicates") or []
        if not reps:
            return {"status": "unavailable", "reason": "no replicates"}
        rep0 = reps[0]
        qc = rep0.get("overlap_qc")
        if not isinstance(qc, dict):
            return {"status": "unavailable",
                    "reason": "no overlap_qc on replicate-0"}
        return qc

    cp4_qc = _ep_report(cp4)
    wt_qc = _ep_report(wt)

    mins = [
        qc["min_adjacent_O"] for qc in (cp4_qc, wt_qc)
        if qc.get("status") == "ok" and qc.get("min_adjacent_O") is not None
    ]
    combined_min = float(min(mins)) if mins else None
    any_flagged = any(
        bool(qc.get("flagged_pairs"))
        for qc in (cp4_qc, wt_qc) if qc.get("status") == "ok"
    )
    threshold = next(
        (qc["threshold"] for qc in (cp4_qc, wt_qc)
         if qc.get("threshold") is not None),
        OVERLAP_FLAG_THRESHOLD,
    )
    return {
        "basis": "mbar_compute_overlap",
        "threshold": float(threshold),
        "min_adjacent_O": combined_min,
        "any_flagged": any_flagged,
        "cp4": cp4_qc,
        "wt": wt_qc,
        "note": (
            "REPORTING-ONLY MBAR overlap QC: O<threshold is a FLAG (phase-space "
            "collapse signal), NOT a gate; O is never folded into any error "
            "bar. ddint_free_kcal + all error bars are unchanged by this field."
        ),
    }


def _build_closure_quarantine(
    cp4: Dict[str, Any],
    wt: Dict[str, Any],
    threshold: float = CLOSURE_DISTRUST_THRESHOLD,
) -> Dict[str, Any]:
    """Assemble the top-level closure-quarantine digest from each endpoint's
    replicate-0 ``closure_distrust`` report (REPORTING-ONLY; ASYMMETRIC).

    Reads (does not compute) the per-replicate ``closure_distrust`` already
    attached to each endpoint by ``analyze_one_leg``, plus the existing
    ``dgb_kcal`` per replicate for the ``per_rep_closure`` list. Lists the
    endpoints whose |closure| exceeds ``threshold`` in ``flagged``; emits a
    WARNING string for them. ASYMMETRIC: there is NO trust / pass / green key —
    an endpoint NOT in ``flagged`` is "no closure red-flag", never "trusted".
    When an endpoint's closure is unavailable it is marked
    ``closure_unavailable`` (no crash). This digest NEVER feeds a ΔG or error
    bar; the closure numbers are byte-identical to closure_cp4/wt_kcal_abs.
    """
    per_rep_closure: List[Dict[str, Any]] = []
    flagged: List[str] = []
    unavailable: List[str] = []

    for label, block in (("cp4", cp4), ("wt", wt)):
        reps = block.get("replicates") or []
        for i, rep in enumerate(reps):
            rid = f"{label}/rep{i}"
            cd = rep.get("closure_distrust") if isinstance(rep, dict) else None
            if not isinstance(cd, dict) or cd.get("status") != "ok":
                # Live path attaches result["closure_distrust"]; the v2-archive
                # re-emit path does NOT (it predates the flag) but DOES carry the
                # per-replicate ``dgb_kcal``. Derive the distrust label from that
                # via the EXISTING graceful wrapper (no tier-logic duplication):
                # compute_closure_distrust reads rep["dgb_kcal"] -> ok report or
                # a closure_unavailable status when dgb_kcal is absent/non-finite.
                if isinstance(rep, dict):
                    cd = compute_closure_distrust(rep, threshold=threshold)
            if not isinstance(cd, dict) or cd.get("status") != "ok":
                # closure_unavailable only when BOTH closure_distrust AND a finite
                # dgb_kcal are missing.
                reason = (
                    cd.get("reason") if isinstance(cd, dict)
                    else "no closure_distrust on replicate"
                )
                per_rep_closure.append({
                    "rep": rid,
                    "status": "closure_unavailable",
                    "reason": reason,
                })
                unavailable.append(rid)
                continue
            entry = {
                "rep": rid,
                "status": "ok",
                "closure_kcal_abs": cd["closure_kcal_abs"],
                "distrust": bool(cd["distrust"]),  # ASYMMETRIC (True flags only)
                "tier": cd["tier"],
            }
            per_rep_closure.append(entry)
            if cd["distrust"]:
                flagged.append(rid)

    if flagged:
        flagged_str = ", ".join(
            f"{e['rep']} |closure|={e['closure_kcal_abs']:.3f}"
            for e in per_rep_closure
            if e.get("status") == "ok" and e["rep"] in flagged
        )
        warning = (
            f"CLOSURE_DISTRUST: {len(flagged)} leg(s) with |closure| > "
            f"{float(threshold):.1f} kcal/mol: {flagged_str}. ADVISORY "
            f"distrust-only (NEVER asserts trust) — review directional "
            f"consistency. A flagged densified38v4 dminus leg may be NON-MIRROR "
            f"hysteresis (structural), NOT necessarily non-convergence. Does "
            f"NOT gate / rescale ΔG and is NEVER folded into any error bar."
        )
    else:
        warning = None

    return {
        "basis": "abs_dgbind1_minus_dgbind2",
        "threshold": float(threshold),
        "per_rep_closure": per_rep_closure,
        "flagged": flagged,            # ASYMMETRIC: distrust ids ONLY (no trust list)
        "unavailable": unavailable,
        "warning": warning,
        "note": (
            "advisory distrust-only; non-mirror per-direction hysteresis is "
            "partly structural, not necessarily non-convergence. A flagged "
            "dminus leg reads 'review: may be non-mirror hysteresis', not "
            "'broken'. distrust=False is NOT a trust assertion; closure is a "
            "SEPARATE axis NEVER folded into σ_btwn / paired-SEM. The closure "
            "values are byte-identical to closure_cp4/wt_kcal_abs."
        ),
    }


def reemit_v3_from_v2(
    v2_path: str,
    output_path: str,
    notes: Optional[str] = None,
) -> Dict[str, Any]:
    """Re-derive a v3 payload from an existing v2 artifact's stored blocks.

    Numbers-neutral: the cp4/wt blocks (with per-replicate dgbind1) are taken
    VERBATIM from the v2 JSON, so the point estimate and every per-replicate
    value are byte-identical to what shipped; only the uncertainty
    representation is upgraded (paired SEM + stats added; STDEV quadrature
    demoted to a legacy field). Use when the source .out files require a
    different conda env (atom_openmm) than the analysis runner — this path
    needs neither openmm nor atom_openmm.

    Raises ValueError if the source is not a v2 artifact (guards against
    re-emitting an already-v3 file or an unrelated JSON).
    """
    with open(v2_path) as fh:
        v2 = json.load(fh)
    schema = v2.get("schema_version")
    if schema != "trackb_ddint_free_v2":
        raise ValueError(
            f"--reemit-v3-from expects a trackb_ddint_free_v2 artifact; got "
            f"schema_version={schema!r} at {v2_path}."
        )
    for key in ("cp4", "wt"):
        if key not in v2 or "replicates" not in v2[key]:
            raise ValueError(
                f"v2 artifact {v2_path} missing '{key}.replicates' — cannot "
                f"re-derive paired statistics."
            )
    payload = build_v3_payload(
        cp4=v2["cp4"], wt=v2["wt"],
        protocol=v2.get("protocol", ""),
        mintimeid=v2.get("mintimeid"), maxtimeid=v2.get("maxtimeid"),
        notes=(notes if notes is not None else v2.get("notes", "")),
    )
    payload["reemit_provenance"] = {
        "reemitted_from": os.path.abspath(v2_path),
        "source_schema": schema,
        "method": "blocks_verbatim_paired_stats_added",
    }
    out_abs = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(out_abs), exist_ok=True)
    with open(out_abs, "w") as fh:
        json.dump(payload, fh, indent=2)
    return payload


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
        "--reemit-v3-from",
        default=None,
        metavar="V2_JSON",
        help=(
            "Re-derive a v3 payload from an existing trackb_ddint_free_v2 "
            "artifact's stored cp4/wt blocks (numbers-neutral: point estimate "
            "+ per-replicate values byte-identical; adds paired-difference SEM "
            "+ stats, demotes the STDEV quadrature to a legacy field). Needs "
            "no openmm / atom_openmm — skips the UWHAM re-analysis entirely. "
            "Writes to --output."
        ),
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

    # Re-emit path: re-derive v3 from an archived v2 artifact (no UWHAM,
    # no openmm/atom_openmm). Numbers-neutral — blocks copied verbatim.
    if args.reemit_v3_from is not None:
        payload = reemit_v3_from_v2(
            v2_path=args.reemit_v3_from,
            output_path=args.output,
            notes=(args.notes if args.notes else None),
        )
        out_path = os.path.abspath(args.output)
        paired = payload.get("paired_difference_stats")
        print(f"[reemit] {args.reemit_v3_from} -> {out_path} "
              f"(schema={payload['schema_version']})")
        print(f"   ddint_free={payload['ddint_free_kcal']:+.4f} kcal/mol "
              f"(point estimate unchanged)")
        if paired and paired.get("sem") is not None:
            ci = paired.get("ci95")
            print(f"   paired SEM={paired['sem']:.4f}  SD={paired['paired_sd']:.4f}  "
                  f"r={paired['pearson_r']:+.4f}  "
                  f"t({paired['df']})={paired['t_stat']:.4f}  "
                  f"95%CI=[{ci[0]:+.3f}, {ci[1]:+.3f}]  "
                  f"sign={paired['sign_status']}")
        # A6/A7 ADVISORY prompts on the re-emit path too (REPORTING-ONLY).
        _sa = payload.get("sign_vs_anchor")
        if isinstance(_sa, dict):
            print(f"   sign_vs_anchor={_sa.get('sign_vs_anchor')} "
                  f"(provisional={_sa.get('provisional')})")
            if _sa.get("warning"):
                print(f"   ⚠ {_sa['warning']}")
        _er = payload.get("extend_recommended")
        if isinstance(_er, dict) and _er.get("recommended"):
            print(f"   ⚠ EXTEND_RECOMMENDED: {_er.get('reason')}")
        print(f"written: {out_path}")
        return 0

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

    payload = build_v3_payload(
        cp4=cp4, wt=wt,
        protocol=args.protocol_tag,
        mintimeid=mintimeid, maxtimeid=args.maxtimeid,
        notes=args.notes,
    )
    paired = payload["paired_difference_stats"]
    payload_error_basis = payload["error_basis"]
    ddint_free = payload["ddint_free_kcal"]
    ddint_err_quadrature_legacy = payload["ddint_err_quadrature_legacy_kcal"]
    closure_cp4 = payload["closure_cp4_kcal_abs"]
    closure_wt = payload["closure_wt_kcal_abs"]

    out_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    # ---- console summary ----
    def _fmt_err(v: Optional[float]) -> str:
        return f"{v:.3f}" if v is not None else "n/a"

    print()
    print("=" * 70)
    if paired is not None and paired.get("sem") is not None:
        _err = paired["sem"]
        print(f"ΔΔG_int_free  =  {ddint_free:+.3f} ± {_err:.3f} kcal/mol "
              f"(ranking-only; error_basis={payload_error_basis} [paired SEM]; "
              f"legacy quadrature STDEV={ddint_err_quadrature_legacy:.3f})")
        _ci = paired.get("ci95")
        if _ci is not None:
            print(f"   paired: SD={paired['paired_sd']:.4f}  "
                  f"SEM={paired['sem']:.4f}  r={_fmt_err(paired['pearson_r'])}  "
                  f"t({paired['df']})={_fmt_err(paired['t_stat'])}  "
                  f"95%CI=[{_ci[0]:+.3f}, {_ci[1]:+.3f}]  "
                  f"sign={paired['sign_status']}")
    else:
        print(f"ΔΔG_int_free  =  {ddint_free:+.3f} ± "
              f"{ddint_err_quadrature_legacy:.3f} kcal/mol "
              f"(ranking-only; error_basis={payload_error_basis} "
              f"[quadrature legacy — paired stats unavailable])")
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
    # P5 closure-quarantine WARNING (REPORTING-ONLY; ASYMMETRIC). Distrust-only
    # prompt — does NOT gate / rescale ΔG.
    _cq = payload.get("closure_quarantine")
    if isinstance(_cq, dict) and _cq.get("warning"):
        print(f"   ⚠ {_cq['warning']}")
    # A6 sign-vs-anchor ADVISORY (REPORTING-ONLY; sign-only). Discordant prompt
    # — does NOT gate / rescale ΔG; provisional on the free leg.
    _sa = payload.get("sign_vs_anchor")
    if isinstance(_sa, dict):
        print(f"   sign_vs_anchor={_sa.get('sign_vs_anchor')} "
              f"(computed_sign={_sa.get('computed_sign_status')}, "
              f"provisional={_sa.get('provisional')})")
        if _sa.get("warning"):
            print(f"   ⚠ {_sa['warning']}")
    # A7 extend-recommended ADVISORY (REPORTING-ONLY). Advice prompt — does NOT
    # change n / gate / touch σ.
    _er = payload.get("extend_recommended")
    if isinstance(_er, dict) and _er.get("recommended"):
        print(f"   ⚠ EXTEND_RECOMMENDED: {_er.get('reason')}")
    print("=" * 70)
    print(f"written: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
