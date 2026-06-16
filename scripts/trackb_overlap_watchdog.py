# -*- coding: utf-8 -*-
"""Track B per-direction MBAR adjacent-overlap WATCHDOG (alert-only).

WHY this exists
---------------
The production gate digest (``trackb_inplace_rbfe_production.py::
_collect_gate_digest``) gates on C3/C4/C9 only. The MBAR ADJACENT overlap —
the metric that distinguishes a genuine λ-ladder traversal from an
occupancy-only false-green (occupancy ≠ mixing) or a frozen-plateau
deterministic-offset ΔG — is *reporting-only* there: it never auto-halts. A
human had to read the postprocess output to catch the three historical
false-greens. With ~30 campaigns nobody eyeballs every overlap. This watchdog
re-computes the adjacent overlaps with the ROBUST MBAR harness, pools them
across matched-seed replicates, bootstrap-CIs each adjacent pair, and emits a
machine-readable ALERT (stdout + JSON, ``exit 1``) when a pooled pair is below
the hard floor or its CI contains zero.

WHAT it is NOT (HARKing guard)
------------------------------
This is ALERT + RECOMMENDATION ONLY. It NEVER densifies, NEVER re-runs, NEVER
mutates a schedule. The only recommendation it prints is "insert 1 λ-bridge
window" (λ-spacing lever only — window count is the dominant mixing lever, the
soft-core constants U0/α are secondary and stay frozen). The actual bridge
insertion is a human review decision. There is deliberately NO call to any
production / densify / rerun entry point in this module.

REUSE (no re-implementation)
----------------------------
The overlap math is the EXISTING robust harness committed at df2d9e6:
  - ``utils.adaptive_lambda.overlap.extract_samples_by_state`` — the GPU-free
    per-direction ``.out`` reader (cols 3-9 by state across all walkers).
  - ``utils.adaptive_lambda.overlap.mbar_overlap_matrix`` — builds the TRUE
    soft-core ``neg_pot`` (``atom_openmm.uwham._bias_fcn`` / ``_npot_fcn``) and
    runs MBAR with ``solver_protocol='robust'`` + the convergence guard
    (diagonal ≤ 1+tol, row-sum ≈ 1) that rejects stalled-solver garbage as
    None.
The per-state SSOT schedule arrays come from
``scripts.trackb_uwham_postprocess._parse_cntl_schedule`` (the per-direction
``<jobname>_<direction>_asyncre.cntl``).

CPU-only robust path
--------------------
pymbar 4.2 routes through JAX; the cuBLAS path has an error history on this
host, so this module forces the JAX CPU backend (``JAX_PLATFORM_NAME=cpu`` /
``JAX_PLATFORMS=cpu`` / ``CUDA_VISIBLE_DEVICES=""``) at import time, BEFORE any
pymbar/JAX import. Overlap is a cheap CPU solve; no GPU is wanted here.

ENV
---
Run with the ``atm`` interpreter (atom_openmm + pymbar):
    /home/san/miniconda3/envs/atm/bin/python3.11 scripts/trackb_overlap_watchdog.py ...

Python 3.8+ (typing.Optional / List / Dict; no PEP 604 unions).
"""

from __future__ import annotations

# --- Force the JAX CPU backend BEFORE pymbar/JAX import (cuBLAS error history).
# setdefault so an explicit caller override is still honoured.
import os as _os_for_env

_os_for_env.environ.setdefault("JAX_PLATFORM_NAME", "cpu")
_os_for_env.environ.setdefault("JAX_PLATFORMS", "cpu")
_os_for_env.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

import argparse
import glob
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Defaults (verdict-anchored; see module docstring)
# ---------------------------------------------------------------------------
# HARD floor: a pooled adjacent overlap below this (OR a CI that contains 0) is
# a genuine ladder break → exit 1. 0.03 is the postprocess COLLAPSE tier.
DEFAULT_HARD_FLOOR = 0.03
# THIN flag: 0.03 ≤ pooled < 0.10 is informational (e.g. the leg-switch knee at
# pairs 5-6 is physically thin but stable) → NOT exit 1 (cry-wolf guard).
DEFAULT_THIN_FLAG = 0.10
# Equilibration discard: keep cycle index >= mintimeid (1-based per walker), so
# warmup_cycles = mintimeid - 1 (drop the first mintimeid-1 cycles per walker).
DEFAULT_MINTIMEID = 100
DEFAULT_N_BOOT = 2000
JOBNAME = "trackb"
DIRECTION_TAGS = ("dplus", "dminus")

VERDICT_HEALTHY = "HEALTHY"
VERDICT_ALERT = "ALERT"
PAIR_HARD_ALERT = "HARD_ALERT"
PAIR_THIN_FLAG = "THIN_FLAG"
PAIR_OK = "OK"
PAIR_UNAVAILABLE = "unavailable"


# ---------------------------------------------------------------------------
# Lazy import of the EXISTING robust overlap harness + the postprocess schedule
# parser (kept lazy so ``--help`` / arg parsing works even in an env without
# atom_openmm; the actual compute requires the atm env).
# ---------------------------------------------------------------------------
def _import_overlap_harness():
    """Import the committed robust overlap harness + the postprocess schedule
    parser. Returns ``(overlap_module, parse_cntl_schedule)``.

    Raises ImportError (with an env hint) if the modules are not importable —
    the watchdog's whole job IS the overlap compute, so a missing harness is a
    hard failure (NOT a silent HEALTHY), unlike a per-pair MBAR non-convergence
    which degrades to ``unavailable``."""
    here = os.path.dirname(os.path.abspath(__file__))
    proj = os.path.dirname(here)
    utils_dir = os.path.join(proj, "utils")
    if utils_dir not in sys.path:
        sys.path.insert(0, utils_dir)
    if proj not in sys.path:
        sys.path.insert(0, proj)
    try:
        from adaptive_lambda import overlap as overlap_mod  # type: ignore
    except Exception:  # pragma: no cover - exercised only in non-atm envs
        from utils.adaptive_lambda import overlap as overlap_mod  # type: ignore
    # Reuse the postprocess cntl schedule parser (SSOT) without importing the
    # whole postprocess module's atom_openmm-dependent surface at import time.
    import importlib.util

    pp_path = os.path.join(here, "trackb_uwham_postprocess.py")
    spec = importlib.util.spec_from_file_location(
        "trackb_uwham_postprocess_wd", pp_path)
    if spec is None or spec.loader is None:  # pragma: no cover
        raise ImportError("cannot locate trackb_uwham_postprocess.py")
    pp_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pp_mod)  # type: ignore[union-attr]
    return overlap_mod, pp_mod._parse_cntl_schedule


# ---------------------------------------------------------------------------
# Replicate discovery + per-(endpoint, direction) adjacent-overlap extraction.
# ---------------------------------------------------------------------------
def discover_replicate_dirs(out_root: str, endpoint: str, leg: str
                            ) -> List[str]:
    """Return the sorted ``rep*`` directories for one (endpoint, leg).

    Layout: ``<out_root>/<endpoint>/<leg>/rep{N}/`` (the production launcher's
    per-replicate dir). Sorted numerically by the rep index."""
    leg_root = os.path.join(out_root, endpoint, leg)
    rep_dirs = glob.glob(os.path.join(leg_root, "rep*"))

    def _rep_key(path: str) -> Tuple[int, str]:
        base = os.path.basename(path)
        if base.startswith("rep") and base[3:].isdigit():
            return (int(base[3:]), base)
        return (10 ** 9, base)

    return [d for d in sorted(rep_dirs, key=_rep_key) if os.path.isdir(d)]


def adjacent_overlaps_for_rep(
    overlap_mod,
    parse_cntl_schedule,
    rep_dir: str,
    direction: str,
    mintimeid: int,
    maxtimeid: Optional[int],
) -> Optional[List[Optional[float]]]:
    """Adjacent-pair MBAR overlaps ``[O_{0,1}, O_{1,2}, ...]`` for ONE rep's
    ONE direction, using the EXISTING robust harness.

    Returns ``None`` when the harness could not produce a matrix at all (pymbar
    / atom_openmm absent, schedule missing, or MBAR non-convergence → the
    harness returns None) — the caller records the whole rep/direction as
    unavailable. When a matrix IS returned, the list has ``K-1`` adjacent
    off-diagonals; individual entries are never None (the matrix is complete),
    but the signature keeps Optional for forward-compat.

    Time window: ``mintimeid`` (1-based cycle) maps to the harness'
    ``warmup_cycles = mintimeid - 1`` (drop the first mintimeid-1 cycles per
    walker — equilibration discard). ``maxtimeid`` (upper cap) is applied here
    by trimming the per-state sample arrays, because the harness reader only
    takes a lower warmup.
    """
    subdir = os.path.join(rep_dir, direction)
    if not os.path.isdir(subdir):
        return None
    cntl_path = os.path.join(subdir, "%s_%s_asyncre.cntl" % (JOBNAME, direction))
    schedule = parse_cntl_schedule(cntl_path)
    if schedule is None:
        return None

    warmup = max(0, mintimeid - 1)
    samples_by_state = overlap_mod.extract_samples_by_state(
        rep_dir, JOBNAME, direction, warmup_cycles=warmup)
    if not samples_by_state:
        return None

    # Upper-cap (maxtimeid): the harness reader keeps cols-3..9 rows per state
    # in walker-then-cycle order; after the lower warmup the per-state arrays
    # hold the post-warmup samples. maxtimeid trims to keep at most
    # (maxtimeid - warmup) samples per walker-equivalent. To stay faithful to
    # the per-walker cycle window WITHOUT re-reading, we cap the per-state
    # sample COUNT to (maxtimeid - mintimeid + 1) * n_walkers — but n_walkers is
    # not known per state here, so we apply the cap conservatively per state by
    # truncating to the post-warmup count implied by maxtimeid only when
    # maxtimeid is set AND smaller than the available cycle span. In practice
    # the completed campaigns have a uniform cycle count, so the upper cap is a
    # rarely-used analysis convenience; the lower warmup is the load-bearing one.
    if maxtimeid is not None:
        keep_cycles = maxtimeid - warmup
        if keep_cycles < 1:
            return None
        n_walkers = _count_walkers(subdir, direction)
        if n_walkers >= 1:
            cap = keep_cycles * n_walkers
            samples_by_state = {
                s: arr[:cap] for s, arr in samples_by_state.items()
            }

    schedule_arrays = _schedule_to_arrays(schedule)
    matrix = overlap_mod.mbar_overlap_matrix(
        samples_by_state, schedule=schedule_arrays, temperature_K=300.0)
    if matrix is None:
        return None
    m = np.asarray(matrix, dtype=float)
    k = m.shape[0]
    if k < 2:
        return None
    return [float(m[i, i + 1]) for i in range(k - 1)]


def _count_walkers(subdir: str, direction: str) -> int:
    """Count ``r*`` walker dirs that hold a ``<jobname>_<direction>.out``."""
    out_basename = "%s_%s.out" % (JOBNAME, direction)
    n = 0
    for w in glob.glob(os.path.join(subdir, "r*")):
        if os.path.isdir(w) and os.path.isfile(os.path.join(w, out_basename)):
            n += 1
    return n


def _schedule_to_arrays(schedule: Dict[str, Any]) -> Dict[str, Any]:
    """Project the postprocess cntl-schedule dict onto the soft-core SSOT keys
    ``mbar_overlap_matrix`` consumes (``lambda1``/``lambda2``/``alpha``/``u0``/
    ``w0``), as numpy arrays indexed by GLOBAL (here per-direction) state id."""
    out: Dict[str, Any] = {}
    for key in ("lambda1", "lambda2", "alpha", "u0", "w0"):
        if key in schedule and schedule[key] is not None:
            out[key] = np.asarray(schedule[key], dtype=float)
    return out


# ---------------------------------------------------------------------------
# Pooling across matched-seed replicates + per-pair bootstrap CI.
# ---------------------------------------------------------------------------
def pool_pairs_across_reps(
    per_rep_pairs: List[Optional[List[Optional[float]]]],
) -> Tuple[int, List[List[float]]]:
    """Collect each adjacent pair's overlap value across reps.

    ``per_rep_pairs`` is one ``[O_{0,1}, ...]`` list per rep (or None for an
    unavailable rep). Returns ``(n_pairs, by_pair)`` where ``by_pair[p]`` is the
    list of finite overlap values for adjacent pair ``p`` across all available
    reps. ``n_pairs`` is the max K-1 seen across available reps (0 if none)."""
    available = [p for p in per_rep_pairs if p is not None]
    if not available:
        return 0, []
    n_pairs = max(len(p) for p in available)
    by_pair: List[List[float]] = [[] for _ in range(n_pairs)]
    for p in available:
        for idx, val in enumerate(p):
            if val is None:
                continue
            v = float(val)
            if np.isfinite(v):
                by_pair[idx].append(v)
    return n_pairs, by_pair


def bootstrap_pair_ci(
    values: List[float],
    n_boot: int,
    rng: np.random.Generator,
    ci: float = 0.95,
) -> Tuple[float, float, float]:
    """Bootstrap the rep-pooled mean of one adjacent pair.

    Resamples the per-rep overlap values WITH replacement (rep is the
    resampling unit), recomputes the mean ``n_boot`` times, and returns
    ``(pooled_mean, ci_lo, ci_hi)`` at the central ``ci`` level. With a single
    rep the CI collapses to the point value (degenerate, flagged downstream by
    the n<5 cry-wolf note). Empty ``values`` → ``(nan, nan, nan)``."""
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return float("nan"), float("nan"), float("nan")
    pooled = float(np.mean(arr))
    if arr.size == 1:
        return pooled, pooled, pooled
    lo_q = (1.0 - ci) / 2.0
    hi_q = 1.0 - lo_q
    idx = rng.integers(0, arr.size, size=(n_boot, arr.size))
    boot_means = arr[idx].mean(axis=1)
    ci_lo = float(np.quantile(boot_means, lo_q))
    ci_hi = float(np.quantile(boot_means, hi_q))
    return pooled, ci_lo, ci_hi


def classify_pair(
    pooled: float,
    ci_lo: float,
    hard_floor: float,
    thin_flag: float,
) -> str:
    """Three-tier verdict for one adjacent pair.

    HARD_ALERT : pooled < hard_floor OR ci_lo <= 0 (CI contains 0). Only this
                 tier drives ``exit 1`` (the cry-wolf guard: THIN is advisory).
    THIN_FLAG  : hard_floor <= pooled < thin_flag (physically thin but stable,
                 e.g. the leg-switch knee — informational, never exit 1).
    OK         : pooled >= thin_flag.
    Non-finite pooled → unavailable (handled by the caller)."""
    if not np.isfinite(pooled):
        return PAIR_UNAVAILABLE
    if pooled < hard_floor or (np.isfinite(ci_lo) and ci_lo <= 0.0):
        return PAIR_HARD_ALERT
    if pooled < thin_flag:
        return PAIR_THIN_FLAG
    return PAIR_OK


# ---------------------------------------------------------------------------
# Per-(endpoint, direction) analysis.
# ---------------------------------------------------------------------------
def analyze_endpoint_direction(
    overlap_mod,
    parse_cntl_schedule,
    out_root: str,
    endpoint: str,
    leg: str,
    direction: str,
    mintimeid: int,
    maxtimeid: Optional[int],
    hard_floor: float,
    thin_flag: float,
    n_boot: int,
    rng: np.random.Generator,
) -> Dict[str, Any]:
    """Pool the adjacent overlaps across this (endpoint, direction)'s reps and
    classify each adjacent pair. Returns a JSON-serializable dict."""
    rep_dirs = discover_replicate_dirs(out_root, endpoint, leg)
    per_rep_pairs: List[Optional[List[Optional[float]]]] = []
    per_rep_status: List[str] = []
    for rep_dir in rep_dirs:
        pairs = adjacent_overlaps_for_rep(
            overlap_mod, parse_cntl_schedule, rep_dir, direction,
            mintimeid, maxtimeid)
        per_rep_pairs.append(pairs)
        per_rep_status.append("ok" if pairs is not None else PAIR_UNAVAILABLE)

    n_pairs, by_pair = pool_pairs_across_reps(per_rep_pairs)
    n_reps_available = sum(1 for p in per_rep_pairs if p is not None)

    per_pair: List[Dict[str, Any]] = []
    n_hard = 0
    n_thin = 0
    for p in range(n_pairs):
        vals = by_pair[p]
        if not vals:
            per_pair.append({
                "i": p, "j": p + 1,
                "pooled_overlap": None, "ci_lo": None, "ci_hi": None,
                "n_rep": 0, "verdict": PAIR_UNAVAILABLE,
            })
            continue
        pooled, ci_lo, ci_hi = bootstrap_pair_ci(vals, n_boot, rng)
        verdict = classify_pair(pooled, ci_lo, hard_floor, thin_flag)
        if verdict == PAIR_HARD_ALERT:
            n_hard += 1
        elif verdict == PAIR_THIN_FLAG:
            n_thin += 1
        per_pair.append({
            "i": p, "j": p + 1,
            "pooled_overlap": pooled,
            "ci_lo": ci_lo, "ci_hi": ci_hi,
            "n_rep": len(vals),
            "verdict": verdict,
        })

    return {
        "endpoint": endpoint,
        "direction": direction,
        "leg": leg,
        "n_reps_total": len(rep_dirs),
        "n_reps_available": n_reps_available,
        "per_rep_status": per_rep_status,
        "n_pairs": n_pairs,
        "per_pair": per_pair,
        "n_hard_alert": n_hard,
        "n_thin": n_thin,
    }


# ---------------------------------------------------------------------------
# Top-level run.
# ---------------------------------------------------------------------------
def run_watchdog(
    out_root: str,
    leg: str,
    endpoints: List[str],
    directions: List[str],
    mintimeid: int,
    maxtimeid: Optional[int],
    hard_floor: float,
    thin_flag: float,
    n_boot: int,
    seed: int = 0,
) -> Dict[str, Any]:
    """Analyze every (endpoint, direction), pool/classify, and build the report.

    The report dict is JSON-serializable and carries an ``overall`` verdict
    (ALERT iff any (endpoint, direction) has >= 1 HARD_ALERT pair, else
    HEALTHY). THIN_FLAG never moves ``overall`` (cry-wolf guard)."""
    overlap_mod, parse_cntl_schedule = _import_overlap_harness()
    rng = np.random.default_rng(seed)

    groups: List[Dict[str, Any]] = []
    total_hard = 0
    total_thin = 0
    for endpoint in endpoints:
        for direction in directions:
            grp = analyze_endpoint_direction(
                overlap_mod, parse_cntl_schedule, out_root,
                endpoint, leg, direction, mintimeid, maxtimeid,
                hard_floor, thin_flag, n_boot, rng)
            groups.append(grp)
            total_hard += grp["n_hard_alert"]
            total_thin += grp["n_thin"]

    overall = VERDICT_ALERT if total_hard > 0 else VERDICT_HEALTHY
    return {
        "out_root": out_root,
        "leg": leg,
        "endpoints": list(endpoints),
        "directions": list(directions),
        "mintimeid": mintimeid,
        "maxtimeid": maxtimeid,
        "hard_floor": hard_floor,
        "thin_flag": thin_flag,
        "n_boot": n_boot,
        "groups": groups,
        "n_hard_alert": total_hard,
        "n_thin": total_thin,
        "overall": overall,
    }


# ---------------------------------------------------------------------------
# stdout rendering.
# ---------------------------------------------------------------------------
def _fmt_val(v: Optional[float]) -> str:
    return "n/a" if v is None or not np.isfinite(v) else "%.4f" % v


def render_stdout(report: Dict[str, Any]) -> str:
    """One line per adjacent pair + a final RESULT line.

    The RESULT line is the machine hook the orchestrator greps:
    ``RESULT: ALERT ...`` / ``RESULT: HEALTHY``. PushNotification is the
    orchestrator's job, NOT this module's."""
    lines: List[str] = []
    lines.append(
        "Track B overlap watchdog — leg=%s out_root=%s" % (
            report["leg"], report["out_root"]))
    lines.append(
        "  mintimeid=%s maxtimeid=%s hard_floor=%.3f thin_flag=%.3f n_boot=%d"
        % (report["mintimeid"], report["maxtimeid"], report["hard_floor"],
           report["thin_flag"], report["n_boot"]))
    for grp in report["groups"]:
        lines.append(
            "  [%s / %s] reps=%d/%d  hard=%d thin=%d" % (
                grp["endpoint"], grp["direction"],
                grp["n_reps_available"], grp["n_reps_total"],
                grp["n_hard_alert"], grp["n_thin"]))
        n_avail = grp["n_reps_available"]
        for pair in grp["per_pair"]:
            note = ""
            if pair["verdict"] == PAIR_HARD_ALERT:
                note = "  -> recommend: insert 1 lambda-bridge window (lambda-spacing only)"
            if pair["verdict"] in (PAIR_HARD_ALERT, PAIR_THIN_FLAG) and n_avail < 5:
                note += "  [n<5 reps: bootstrap CI noisy]"
            lines.append(
                "    pair %d-%d  O=%s  CI=[%s, %s]  n_rep=%d  %s%s" % (
                    pair["i"], pair["j"],
                    _fmt_val(pair["pooled_overlap"]),
                    _fmt_val(pair["ci_lo"]), _fmt_val(pair["ci_hi"]),
                    pair["n_rep"], pair["verdict"], note))
    if report["overall"] == VERDICT_ALERT:
        lines.append(
            "RESULT: ALERT  %d hard-alert adjacent pair(s) across %d group(s) "
            "(%d thin-flag, advisory)" % (
                report["n_hard_alert"], len(report["groups"]),
                report["n_thin"]))
    else:
        lines.append(
            "RESULT: HEALTHY  0 hard-alert adjacent pairs (%d thin-flag, "
            "advisory)" % report["n_thin"])
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI.
# ---------------------------------------------------------------------------
def _parse_csv(raw: Optional[str]) -> List[str]:
    if raw is None:
        return []
    return [tok.strip() for tok in raw.split(",") if tok.strip()]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Track B per-direction MBAR adjacent-overlap watchdog "
                    "(alert-only; never densifies or re-runs).")
    p.add_argument("--out-root", required=True,
                   help="campaign out-root (e.g. "
                        "outputs/_trackb/twocopy_bound_n3b)")
    p.add_argument("--leg", default="bound", choices=("bound", "free"),
                   help="leg name (default bound)")
    p.add_argument("--endpoints", default="cp4,wt",
                   help="comma-separated endpoints (default cp4,wt)")
    p.add_argument("--directions", default="dplus,dminus",
                   help="comma-separated directions (default dplus,dminus)")
    p.add_argument("--mintimeid", type=int, default=DEFAULT_MINTIMEID,
                   help="keep cycles >= mintimeid (1-based; equilibration "
                        "discard). Default %d." % DEFAULT_MINTIMEID)
    p.add_argument("--maxtimeid", type=int, default=None,
                   help="optional upper cycle cap (default: all)")
    p.add_argument("--hard-floor", type=float, default=DEFAULT_HARD_FLOOR,
                   help="pooled overlap below this (or CI contains 0) -> "
                        "HARD_ALERT (exit 1). Default %.3f." % DEFAULT_HARD_FLOOR)
    p.add_argument("--thin-flag", type=float, default=DEFAULT_THIN_FLAG,
                   help="pooled overlap below this (>= hard-floor) -> THIN_FLAG "
                        "(advisory, no exit 1). Default %.3f." % DEFAULT_THIN_FLAG)
    p.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT,
                   help="bootstrap resamples for the per-pair CI (default %d)"
                        % DEFAULT_N_BOOT)
    p.add_argument("--seed", type=int, default=0,
                   help="bootstrap RNG seed (default 0; deterministic report)")
    p.add_argument("--json", default=None,
                   help="optional path to write the report JSON")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    endpoints = _parse_csv(args.endpoints)
    directions = _parse_csv(args.directions)
    if not endpoints:
        print("ERROR: no endpoints parsed from --endpoints", file=sys.stderr)
        return 2
    if not directions:
        print("ERROR: no directions parsed from --directions", file=sys.stderr)
        return 2
    for d in directions:
        if d not in DIRECTION_TAGS:
            print("ERROR: unknown direction %r (expected one of %s)"
                  % (d, DIRECTION_TAGS), file=sys.stderr)
            return 2

    maxtimeid = args.maxtimeid if (args.maxtimeid is not None
                                   and args.maxtimeid > 0) else None

    report = run_watchdog(
        out_root=args.out_root, leg=args.leg,
        endpoints=endpoints, directions=directions,
        mintimeid=args.mintimeid, maxtimeid=maxtimeid,
        hard_floor=args.hard_floor, thin_flag=args.thin_flag,
        n_boot=args.n_boot, seed=args.seed)

    print(render_stdout(report))

    if args.json:
        os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=2)

    return 1 if report["overall"] == VERDICT_ALERT else 0


if __name__ == "__main__":
    sys.exit(main())
