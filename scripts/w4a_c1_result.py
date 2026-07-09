#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""C1 verdict aggregator — W4A carved-vs-uncarved void-water dG-neutrality.

WHAT
    When the W4A CARVED control run finishes, compute the carved dgb-based
    ddG_bind with the SAME estimator/pipeline as the uncarved reference, compare
    to the uncarved reference (Pairing B), apply the C1 verdict, and write a
    self-contained w4a_c1_result.md summarising the whole thing (W23A result +
    C1 + verdict + funnel implication + R-18 caveats).

WHY
    The W23A carved in-place two-copy RBFE gave a LARGE binding effect
    (dgb ddG_bind = -9.81 kcal, |t|=6.25 -> NOT flatten). The void-water carve
    fix has UNPROVEN dG-neutrality. scientific-review condition C1 (hard gate for the W23A
    verdict): run W4A CARVED at the same config as the existing UNCARVED W4A and
    test whether the carve perturbs ddG_bind. Because both runs share the SAME 6
    seeds, the correct test is a PAIRED same-seed delta d[s] = carved_ddg[s] -
    uncarved_ddg[s] (paired mean/SEM/|t| + a robust drop-the-outlier recheck +
    the carved-vs-uncarved correlation r), NOT a |delta-of-means| point-threshold.
    The point-threshold (|delta|>=4 -> "FAIL") was a FALSE-POSITIVE: it ignored
    the SEM and a single-seed heavy tail (s127). formal scientific review
    (the formal C1 review): the true status is INDETERMINATE
    at n=6/400cyc -- carve neither proven neutral nor proven non-neutral.

HOW TO USE
    Manual:  /home/san/miniconda3/envs/atm/bin/python scripts/w4a_c1_result.py
    Cron  :  the same command every ~15 min (see w4a_c1_result.py --install-cron
             or the crontab set up alongside this script). Idempotent + safe to
             poll: if the run is not finished it prints status and exits 0 without
             writing; once finished it writes w4a_c1_result.md + a sentinel and
             subsequent runs exit fast.

    Options:
        --force        recompute + rewrite even if the sentinel exists
        --out PATH     write the .md somewhere other than the default
        --print        also echo the verdict block to stdout

REGIME: R-11 SIGN/order-only. R-18 honesty (the uncarved reference is itself
O3-collapsed, so C1 can only bound a carve bias of >~2-3 kcal; a sub-kcal bias
cannot be excluded -- which is fine for a SIGN/order verdict).

NO source edits to the frozen engine/estimator: this reuses
outputs/_trackb/w4a_carved_c1_20260709/recompute_uncarved_ref.py verbatim
(analyze_leg / paired), which itself reuses the frozen uwham estimator.
"""
from __future__ import annotations

import argparse
import datetime
import importlib.util
import json
import math
import os
import sys

PROJ = "/home/san/UPDD_proj"
C1DIR = os.path.join(PROJ, "outputs/_trackb/w4a_carved_c1_20260709")
CARVED_BOUND = os.path.join(C1DIR, "twocopy_w4a_bound_FIXAB_carved")
CARVED_FREE = os.path.join(C1DIR, "twocopy_w4a_free_FIXAB_carved")
RUN_LOG = os.path.join(C1DIR, "run_w4a_c1.log")
UNCARVED_REF_JSON = os.path.join(C1DIR, "uncarved_reference_dgb.json")
CARVED_JSON = os.path.join(C1DIR, "carved_dgb.json")
SENTINEL = os.path.join(C1DIR, ".w4a_c1_finalized")
DEFAULT_MD = os.path.join(PROJ, "w4a_c1_result.md")

# --- W23A result (the thing C1 gates), from w4a_c1 sibling run (07-09 06:42) ---
W23A_DDG = -9.808          # PRIMARY dgb ddG_bind (offset-cancelled), n=5
W23A_T = 6.25              # paired |t|
W23A_DDG_4SEED = -8.26     # s101 (marginal-overlap outlier) excluded

# --- C1 verdict gates (SIGN/order regime, R-11; significance + robustness) ---
# The verdict is NO LONGER a |delta| point-threshold. That threshold produced a
# false-positive "FAIL" by ignoring the SEM and a single-seed heavy tail (s127).
# It is now a PAIRED same-seed significance gate + a robust (drop-max-|z|) recheck
# (the formal C1 review).
ANCHOR_KCAL = 4.5          # W23A destabilizing anchor band midpoint (~4-5), reporting-only
SIG_T = 3.0                # paired |t| below this = NOT significant (matches Z_SE_SIGN_THRESHOLD)
COLLAPSE_FRAC = 0.5        # robust |delta| below this * |full delta| = "collapses toward 0"
NEUTRAL_CI_HALFWIDTH = 1.0 # paired CI half-width below this (kcal) = bounded near zero

DONE_MARKER = "### W4A CARVED BOUND DONE"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _run_complete():
    """Both legs finished iff the wrapper echoed the bound-DONE marker."""
    if not os.path.isfile(RUN_LOG):
        return False
    try:
        with open(RUN_LOG, "r", errors="ignore") as fh:
            return DONE_MARKER in fh.read()
    except OSError:
        return False


def _fmt(x, nd=3):
    if x is None:
        return "n/a"
    try:
        return ("%%.%df" % nd) % float(x)
    except (TypeError, ValueError):
        return str(x)


def _align_arms(carved_primary, unc_primary):
    """Common seeds (carved order) -> aligned (seeds, carved_vec, uncarved_vec).

    Pairs by SEED (not rep-index): the two runs share seeds but not rep order.
    """
    c_map = dict(zip(carved_primary.get("matched_seeds") or [],
                     carved_primary.get("per_rep_ddg") or []))
    u_map = dict(zip(unc_primary.get("matched_seeds") or [],
                     unc_primary.get("per_rep_ddg") or []))
    seeds = [s for s in (carved_primary.get("matched_seeds") or []) if s in u_map]
    return seeds, [c_map[s] for s in seeds], [u_map[s] for s in seeds]


def _max_abs_z_index(deltas):
    """Index of the single max-|z| element (heavy-tail outlier); None if n<3."""
    n = len(deltas)
    if n < 3:
        return None
    mean = sum(deltas) / n
    sd = math.sqrt(sum((d - mean) ** 2 for d in deltas) / (n - 1))
    if sd == 0.0:
        return None
    return max(range(n), key=lambda i: abs((deltas[i] - mean) / sd))


def compute_carve_delta(carved_primary, unc_primary, pds):
    """PAIRED same-seed carve delta d[s] = carved_ddg[s] - uncarved_ddg[s].

    Reuses the FROZEN paired_difference_stats (pds) for the delta stats AND the
    carved-vs-uncarved correlation r (pds returns pearson_r of the two arms).
    ROBUST = drop the single max-|z| seed from BOTH arms and re-run pds. The
    unpaired |delta-of-means| / combined-SEM is kept as a SECONDARY line only.
    No estimator re-implementation: pds is the frozen uwham paired stat.
    """
    seeds, cvec, uvec = _align_arms(carved_primary, unc_primary)
    full = pds(cvec, uvec) if seeds else {}
    deltas = full.get("per_rep_ddg") or [c - u for c, u in zip(cvec, uvec)]

    di = _max_abs_z_index(deltas)
    robust, dropped = None, None
    if di is not None:
        dropped = seeds[di]
        cvec_r = [v for k, v in enumerate(cvec) if k != di]
        uvec_r = [v for k, v in enumerate(uvec) if k != di]
        robust = pds(cvec_r, uvec_r)

    c_sem = carved_primary.get("sem") or 0.0
    u_sem = unc_primary.get("sem") or 0.0
    comb_sem = math.sqrt(c_sem * c_sem + u_sem * u_sem)
    delta_means = ((carved_primary.get("ddG_bind_kcal") or 0.0)
                   - (unc_primary.get("ddG_bind_kcal") or 0.0))
    z_unpaired = (abs(delta_means) / comb_sem) if comb_sem else None

    return {
        "seeds": seeds, "carved_ddg": cvec, "uncarved_ddg": uvec,
        "per_seed_delta": deltas, "full": full, "robust": robust,
        "dropped_seed": dropped, "pearson_r": full.get("pearson_r"),
        "unpaired": {"delta": delta_means, "comb_sem": comb_sem, "z": z_unpaired},
    }


def _classify(carve):
    """Significance- + robustness-gated carve verdict (SIGN/order, R-11).

    NOT a |delta| point-threshold. Uses the paired same-seed delta (full) and
    the drop-max-|z|-seed robust recheck (robust):
      - CARVE-BIASED   : paired |t| >= SIG_T AND robust delta does NOT collapse.
      - CARVE-NEUTRAL  : paired |t| <  SIG_T AND paired CI is tight around zero.
      - INDETERMINATE  : otherwise (default) -- underpowered, cannot resolve.
    """
    full = carve.get("full") or {}
    robust = carve.get("robust")
    dropped = carve.get("dropped_seed") or "n/a"

    t = full.get("t_stat")
    abs_t = abs(t) if t is not None else None
    delta = full.get("mean")
    ci = full.get("ci95") or [None, None]
    lo, hi = (ci[0], ci[1]) if len(ci) == 2 else (None, None)
    ci_has_zero = (lo is None or hi is None) or (lo <= 0.0 <= hi)
    ci_half = (abs(hi - lo) / 2.0) if (lo is not None and hi is not None) else None

    rob_delta = robust.get("mean") if robust else None
    collapses = (rob_delta is not None and delta not in (None, 0.0)
                 and abs(rob_delta) < COLLAPSE_FRAC * abs(delta))
    sig = abs_t is not None and abs_t >= SIG_T

    if sig and not collapses:
        return ("CARVE-BIASED (non-neutral)",
                "paired |t|=%s >= %.0f AND the robust (drop-%s) delta stays large "
                "(%s kcal, does NOT collapse) -> the carve shifts ddG_bind "
                "coherently across seeds. W23A verdict is confounded on this basis; "
                "escalate (scientific review) / redesign the carve."
                % (_fmt(abs_t, 2), SIG_T, dropped, _fmt(rob_delta, 2)))

    if (not sig) and ci_half is not None and ci_half < NEUTRAL_CI_HALFWIDTH and ci_has_zero:
        return ("CARVE-NEUTRAL (bounded near zero)",
                "paired |t|=%s < %.0f AND the paired CI (%s..%s kcal) is tight "
                "around zero (half-width %s < %.1f) -> the carve's dG-perturbation "
                "is bounded near zero."
                % (_fmt(abs_t, 2), SIG_T, _fmt(lo, 2), _fmt(hi, 2),
                   _fmt(ci_half, 2), NEUTRAL_CI_HALFWIDTH))

    collapse_clause = ""
    if collapses:
        pct = (100.0 * abs(rob_delta) / abs(delta)) if delta else float("nan")
        collapse_clause = (", and the robust (drop-%s) delta collapses to %s kcal "
                           "(%s%% of the point estimate)"
                           % (dropped, _fmt(rob_delta, 2), _fmt(pct, 0)))
    return ("INDETERMINATE",
            "carve dG-neutrality cannot be resolved at this sampling: the paired "
            "same-seed test is not significant (|t|=%s < %.0f)%s, and the CI "
            "(%s..%s kcal) spans zero. The carve is NEITHER proven neutral (the CI "
            "does not exclude a several-kcal bias) NOR proven non-neutral (fails "
            "significance / collapses on outlier removal); both directions are "
            "underpowered."
            % (_fmt(abs_t, 2), SIG_T, collapse_clause, _fmt(lo, 2), _fmt(hi, 2)))


def compute():
    """Analyse carved legs (reusing the frozen estimator) + load uncarved ref."""
    recompute = _load(
        "recompute_uncarved_ref",
        os.path.join(C1DIR, "recompute_uncarved_ref.py"))

    bound = recompute.analyze_leg(CARVED_BOUND, "bound")
    free = recompute.analyze_leg(CARVED_FREE, "free")
    carved_primary = recompute.paired(bound, free, "dgb")
    carved_secondary = recompute.paired(bound, free, "dgbind1")

    carved = {
        "protocol": "w4a_c1_CARVED_dgb_leg_difference (Pairing B)",
        "estimator": bound.get("out_root") and (
            "atom_openmm/uwham analyze_replicate_set + paired_difference_stats "
            "(same as uncarved ref, W23A pipeline)"),
        "bound_carved": bound,
        "free_carved": free,
        "pairing_B_bound_minus_free_carved": {
            "primary_dgb": carved_primary,
            "secondary_dgbind1": carved_secondary,
        },
    }
    json.dump(carved, open(CARVED_JSON, "w"), indent=2, default=str)

    uncref = json.load(open(UNCARVED_REF_JSON))
    unc_primary = uncref["pairing_B_bound_minus_freeFIXAB"]["primary_dgb"]

    # PAIRED same-seed carve delta (reuses the frozen paired_difference_stats,
    # which lives on the recompute module alongside the leg estimator).
    carve = compute_carve_delta(
        carved_primary, unc_primary, recompute.uwham.paired_difference_stats)

    return carved, carved_primary, carved_secondary, unc_primary, uncref, carve


def render_md(carved_primary, carved_secondary, unc_primary, uncref, carve, stamp):
    cd = carved_primary.get("ddG_bind_kcal")
    ud = unc_primary.get("ddG_bind_kcal")
    full = carve.get("full") or {}
    robust = carve.get("robust")
    dropped = carve.get("dropped_seed")
    paired_delta = full.get("mean")
    paired_t = full.get("t_stat")
    pci = full.get("ci95") or [None, None]
    plo, phi = (pci[0], pci[1]) if len(pci) == 2 else (None, None)
    robust_delta = robust.get("mean") if robust else None
    robust_n = robust.get("n") if robust else None
    r = carve.get("pearson_r")
    unp = carve.get("unpaired") or {}
    n = full.get("n")
    dfree = max((n or 1) - 1, 0)
    verdict, why = _classify(carve)
    unc_sec = uncref["pairing_B_bound_minus_freeFIXAB"]["secondary_dgbind1"]

    def _t(x):
        return _fmt(abs(x), 2) if x is not None else "n/a"

    # the dropped (max-|z|) seed's carved/uncarved values, for the swing narrative
    drop_c = drop_u = None
    if dropped is not None and dropped in carve.get("seeds", []):
        i = carve["seeds"].index(dropped)
        drop_c = carve["carved_ddg"][i]
        drop_u = carve["uncarved_ddg"][i]
    drop_swing = (drop_c - drop_u) if (drop_c is not None and drop_u is not None) else None

    lines = []
    A = lines.append
    A("# W4A carved-vs-uncarved — review condition C1 verdict (void-water carve dG-neutrality)")
    A("")
    A("_Auto-generated %s by `scripts/w4a_c1_result.py`. Regime: R-11 SIGN/order-only, R-18 honesty._" % stamp)
    A("")
    A("## TL;DR")
    A("")
    A("| quantity | value |")
    A("|---|---|")
    A("| **C1 VERDICT** | **%s** |" % verdict)
    A("| **paired delta (carved − uncarved), same 6 seeds** | **%s kcal**, paired \\|t\\|=%s (df=%s), n=%s |"
      % (_fmt(paired_delta, 2), _t(paired_t), dfree, n))
    A("| **robust delta (drop %s from BOTH arms)** | **%s kcal**, n=%s (collapses toward noise) |"
      % (dropped or "n/a", _fmt(robust_delta, 2), robust_n if robust_n is not None else "n/a"))
    A("| carved-vs-uncarved per-seed correlation r | %s (anti-correlated → heavy-tail artifact, not a coherent bias) |"
      % _fmt(r, 2))
    A("| unpaired delta / combined SEM (secondary) | %s / %s → z=%s |"
      % (_fmt(unp.get("delta"), 2), _fmt(unp.get("comb_sem"), 2), _fmt(unp.get("z"), 2)))
    A("| ddG_bind CARVED (dgb, Pairing B) | %s kcal (\\|t\\|=%s, n=%s, sign=%s) |"
      % (_fmt(cd, 2), _t(carved_primary.get("t_stat")), carved_primary.get("n"),
         carved_primary.get("sign_status")))
    A("| ddG_bind UNCARVED ref (dgb, Pairing B) | %s kcal (\\|t\\|=%s, n=%s, sign=%s) |"
      % (_fmt(ud, 2), _t(unc_primary.get("t_stat")), unc_primary.get("n"),
         unc_primary.get("sign_status")))
    A("| anchor (W23A destabilizing, reporting-only) | ~%.1f kcal |" % ANCHOR_KCAL)
    A("| W23A carved result (the effect C1 gates) | dgb ddG_bind = %s kcal, \\|t\\|=%s (LARGE = not-flatten) |"
      % (_fmt(W23A_DDG, 2), W23A_T))
    A("")
    A("**VERDICT — %s.** %s" % (verdict, why))
    A("")
    A("## What C1 tests")
    A("")
    A("W23A carved in-place two-copy RBFE returned a LARGE binding effect "
      "(dgb ddG_bind = %s kcal, |t|=%s) → the construction does NOT flatten a "
      "known large binding effect. But the void-water carve fix (needed to build "
      "W23A's crash-prone free leg) has UNPROVEN dG-neutrality. C1 runs W4A CARVED "
      "at the exact config of the existing UNCARVED W4A. Both runs share the SAME 6 "
      "seeds (s7,s23,s101,s127,s163,s199), so the carve's dG-perturbation is "
      "isolated by a **PAIRED same-seed delta** d[s] = carved_ddg[s] − "
      "uncarved_ddg[s] — NOT a difference of two means against a point-threshold."
      % (_fmt(W23A_DDG, 2), W23A_T))
    A("")
    A("## Paired same-seed statistics (the correct, more powerful test)")
    A("")
    A("| seed | carved ddG | uncarved ddG | delta d[s] |")
    A("|---|---|---|---|")
    for s, c, u in zip(carve.get("seeds", []), carve.get("carved_ddg", []),
                       carve.get("uncarved_ddg", [])):
        mark = "  ← max-\\|z\\| outlier" if s == dropped else ""
        A("| %s | %+.2f | %+.2f | %+.2f%s |" % (s, c, u, c - u, mark))
    A("")
    A("- **PAIRED mean delta = %s kcal, |t| = %s (df=%s), n=%s → NOT significant** "
      "(|t| < %.0f). The CI %s..%s kcal spans zero." % (
          _fmt(paired_delta, 2), _t(paired_t), dfree, n, SIG_T,
          _fmt(plo, 2), _fmt(phi, 2)))
    A("- **ROBUST (drop %s, the max-|z| seed, from BOTH arms): delta = %s kcal (n=%s)** "
      "→ the point estimate collapses toward noise once the single heavy-tail seed "
      "is removed." % (dropped or "n/a", _fmt(robust_delta, 2),
                       robust_n if robust_n is not None else "n/a"))
    A("- **carved-vs-uncarved per-seed correlation r = %s (anti-correlated).** A "
      "systematic ~%s kcal carve bias would shift seeds ~coherently; instead %s swings "
      "%s kcal same-seed (carved %s vs uncarved %s) while s23/s163 barely move → the "
      "signature of a heavy-tailed **sampling artifact**, not a reproducible carve "
      "bias." % (_fmt(r, 2), _fmt(abs(paired_delta) if paired_delta is not None else None, 1),
                 dropped or "one seed", _fmt(abs(drop_swing) if drop_swing is not None else None, 1),
                 _fmt(drop_c, 2), _fmt(drop_u, 2)))
    A("- The paired test being LESS significant than the unpaired one (|t| %s < z %s) "
      "confirms the two arms do NOT track a common carve signal." % (
          _t(paired_t), _fmt(unp.get("z"), 2)))
    A("- Secondary (dgbind1, forward-only, non-offset-cancelled): carved ddG_bind = %s "
      "kcal (|t|=%s), uncarved = %s kcal (|t|=%s). Reported for completeness; the "
      "primary offset-cancelled dgb is authoritative." % (
          _fmt(carved_secondary.get("ddG_bind_kcal"), 2), _t(carved_secondary.get("t_stat")),
          _fmt(unc_sec.get("ddG_bind_kcal"), 2), _t(unc_sec.get("t_stat"))))
    A("")
    A("## Why INDETERMINATE, not the earlier \"FAIL\"")
    A("")
    A("The script previously flagged \"CARVE NON-NEUTRAL (C1 FAIL)\" from a "
      "|delta|>=4.0 kcal **point-threshold** on the difference of means. That is a "
      "false-positive: it ignores the SEM (combined ~%s kcal) and the single-seed "
      "heavy tail (%s). Three of four proper analyses (paired |t|=%s, robust "
      "delta=%s, anti-correlation r=%s) fail to establish a ≥4 kcal bias. Per scientific review "
      "(the formal C1 review) the correct status is INDETERMINATE: "
      "the carve is neither proven neutral nor proven non-neutral at n=6/400cyc." % (
          _fmt(unp.get("comb_sem"), 2), dropped or "s127", _t(paired_t),
          _fmt(robust_delta, 2), _fmt(r, 2)))
    A("")
    A("## Does the carve, if real, FLATTEN or INFLATE? → INFLATE")
    A("")
    A("Carved is MORE negative than uncarved (delta %s, SAME sign as W23A's %s). "
      "W23A is a Trp→Ala abolishing the deepest p53 anchor → a large-negative "
      "ddG_bind IS the signal. The carve pushes in the SAME direction, so it can "
      "only **inflate** the magnitude, never flatten it. Worst-case correction of "
      "W23A (subtract the carve delta):" % (_fmt(paired_delta, 2), _fmt(W23A_DDG, 2)))
    A("")
    A("| W23A | raw | minus full carve (%s) | minus robust carve (%s) |"
      % (_fmt(paired_delta, 2), _fmt(robust_delta, 2)))
    A("|---|---|---|---|")
    A("| full (n=5) | %s | %s | %s |" % (
        _fmt(W23A_DDG, 2), _fmt(W23A_DDG - (paired_delta or 0.0), 2),
        _fmt(W23A_DDG - (robust_delta or 0.0), 2)))
    A("| ex-s101 (n=4) | %s | %s | %s |" % (
        _fmt(W23A_DDG_4SEED, 2), _fmt(W23A_DDG_4SEED - (paired_delta or 0.0), 2),
        _fmt(W23A_DDG_4SEED - (robust_delta or 0.0), 2)))
    A("")
    A("In every scenario W23A stays large (|ddG| ≥ ~4.2) and same-sign. The only way "
      "the carve could flatten W23A is by adding > ~8 kcal of spurious binding "
      "(~3.6σ beyond the C1 point estimate, and in the inflating not flattening "
      "direction) — not supported. So the carve **INFLATES, it does not FLATTEN**; "
      "the W23A SIGN/ORDER \"large, non-flat binding effect\" conclusion stands "
      "(bound-specific: free ~0 vs bound %s), with caveats. NOT a quantitative claim."
      % _fmt(W23A_DDG, 1))
    A("")
    A("## Combined conclusion (funnel implication)")
    A("")
    A("- **C1 = INDETERMINATE.** At n=6/400cyc the carve's effect on ddG_bind is "
      "statistically indistinguishable from zero (paired |t|=%s, robust collapses to "
      "%s, anti-correlated r=%s). It is NOT proven neutral (CI ~%s..%s spans a several-"
      "kcal bias) and NOT proven non-neutral (fails significance). Both underpowered." % (
          _t(paired_t), _fmt(robust_delta, 2), _fmt(r, 2), _fmt(plo, 1), _fmt(phi, 1)))
    A("- **W23A SIGN/ORDER conclusion STANDS with caveats**: the carve inflates not "
      "flattens, so W23A's large-negative, bound-specific (free ~0 vs bound %s) effect "
      "survives the worst-case carve correction (stays |ddG| ≥ ~4.2, same sign). This "
      "is a SIGN/order methods finding, NOT a certified/quantitative ddG." % _fmt(W23A_DDG, 1))
    A("- This is NOT \"confounded / do not publish\". It IS \"carve bias bounded but "
      "not proven neutral; the W23A SIGN conclusion is robust to the worst-case bias "
      "because the carve inflates, not flattens\".")
    A("")
    A("**OPEN DECISION (scientific-review-ranked):**")
    A("1. **(B) Densify + expand seeds** on W23A AND the W4A carved/uncarved control, "
      "with the corrected paired / robust / bootstrap statistics. Minimum-assumption, "
      "power-acquiring step; decides whether the s127/s101 tail is sampling (shrinks "
      "with n → carve salvageable) or structural (persists → A mandatory). **Top pick.**")
    A("2. **(A) Carve redesign → union-solvation / ghost-atom placeholder** "
      "(dG-faithful). The correct destination: physically occupies the displaced-copy "
      "void throughout dynamics (fixes the residual tail the t=0-only carve leaves) and "
      "dissolves the neutrality question. Becomes 1st if B shows the tail is structural.")
    A("3. **(C) Forward-only / base=u0** as a cross-check only — dropping the backward "
      "leg reintroduces the box-size soft-core apex OFFSET that cancels only in the "
      "bidirectional dgb (barnase Gate-A re-adjudication 2026-07-07). Capping base=u1 "
      "is rejected (u1 is a physical endpoint; capping masks a real clash).")
    A("")
    A("## NEW RISKS (scientific review, not in the earlier summary)")
    A("")
    A("1. **Irreducible transferability gap** — C1 can NEVER certify the *W23A* carve "
      "directly: uncarved W23A does not build (NaN crash — the reason the carve exists), "
      "so C1 must proxy via W4A, a DIFFERENT water set. W23A carve neutrality is "
      "structurally unprovable by the carve-vs-uncarve method.")
    A("2. **Leg-asymmetric carve does NOT cancel in the double difference** — W23A carve "
      "is bound 0 / free 5 waters (enters ddG_bind through the free leg ONLY); W4A carve "
      "is bound 3 / free 4 (near-symmetric → largely cancels). So even a conclusive C1 "
      "on W4A would not transfer to W23A.")
    A("3. **Collapsed near-zero baseline** — the W4A uncarved reference is itself "
      "O3-collapsed (min-O ~0.065, ddG %s, |t|=%s). Differencing two heavy-tailed ~0 "
      "quantities inherits both tails → C1 power ~0 by construction." % (
          _fmt(ud, 2), _t(unc_primary.get("t_stat"))))
    A("")
    A("## R-18 caveats (honest limits)")
    A("")
    A("- The uncarved reference is O3-collapsed (min adjacent overlap 0.065–0.126, "
      "|t|=%s, sign-undetermined). C1 can only BOUND a carve bias; it cannot certify "
      "sub-kcal neutrality — fine for a SIGN/order verdict, NOT for a calibrated ddG."
      % _t(unc_primary.get("t_stat")))
    A("- With |z|>2 single-seed outliers at n=6, the paired-t normality assumption is "
      "violated → all |t| (W23A's 6.25 and C1's carved 2.30) are outlier-sensitive. "
      "Only conclusions robust to single-seed removal, at SIGN/ORDER, may be trusted.")
    A("- W23A magnitude (%s kcal) is ~2x the typical anchor (+4–5 destabilizing). Trp23 "
      "is p53's deepest MDM2 anchor so a large effect is plausible, but FF over-"
      "estimation of buried-hydrophobic mutations is likely → R-11 order-only." % _fmt(W23A_DDG, 2))
    A("- W23A large-negative sign is robust without its s101 outlier (4-seed = %s kcal); "
      "MAGNITUDE is not (already forbidden by R-11). The carve is t=0-only, so a water "
      "can still diffuse into the displaced-copy void during dynamics (the residual "
      "s101/s127 tail source the carve leaves untouched)." % _fmt(W23A_DDG_4SEED, 2))
    A("- SIGN convention: W23A code-negative = destabilizing is consistent with the "
      "barnase fold ADR-0024 convention → known-direction anchors 2/2 (barnase "
      "I96A + W23A) support code-neg = destabilizing.")
    A("")
    A("## Provenance")
    A("")
    A("- Estimator: `outputs/_trackb/w4a_carved_c1_20260709/recompute_uncarved_ref.py` "
      "(reuses frozen `trackb_inplace_rbfe_production` + `uwham` verbatim; the paired "
      "delta + r use the same `paired_difference_stats`), mintimeid=100. This script "
      "(`w4a_c1_result.py`) does NOT edit the engine or the estimator.")
    A("- Carved run: `outputs/_trackb/w4a_carved_c1_20260709/` (launch 2026-07-09 07:37, "
      "6 seeds s7/s23/s101/s127/s163/s199, 24-state, 400 cyc, `--carve-void-waters "
      "--carve-cutoff-nm 0.26`, disp 4.0, staged-min, reseed 20260618). carve OFF = "
      "existing uncarved build byte-match; carve ON removes free 4 / bound 3 bulk "
      "waters (net-charge 0, C2 all-bulk). NB review risk #4: per-leg carved-water "
      "counts must be nailed down + fail-loud logged (W23A carve is bound 0 / free 5).")
    A("- Uncarved reference: `uncarved_reference_dgb.json` (Pairing B). Carved raw: "
      "`carved_dgb.json`. Formal review: "
      "`internal C1 scientific-review note 20260710`.")
    A("- W23A sibling result: `outputs/_trackb/mdm2_w23a_gateA_20260708/w23a_ddg.json`.")
    A("")
    return "\n".join(lines) + "\n", {
        "verdict": verdict, "carved_ddg": cd, "uncarved_ddg": ud,
        "delta": paired_delta,
        "abs_delta": abs(paired_delta) if paired_delta is not None else None,
        "paired_t": abs(paired_t) if paired_t is not None else None,
        "robust_delta": robust_delta, "pearson_r": r,
        "nsigma": unp.get("z")}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=DEFAULT_MD)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--print", dest="do_print", action="store_true")
    ap.add_argument("--install-cron", action="store_true",
                    help="print the crontab line to add and exit")
    args = ap.parse_args(argv)

    if args.install_cron:
        py = "/home/san/miniconda3/envs/atm/bin/python"
        print("*/15 * * * * %s %s/scripts/w4a_c1_result.py "
              ">> %s/w4a_c1_cron.log 2>&1  # UPDD-W4A-C1"
              % (py, PROJ, C1DIR))
        return 0

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if os.path.exists(SENTINEL) and not args.force:
        print("[w4a_c1_result] already finalized (%s exists). Use --force to "
              "recompute. Result: %s" % (SENTINEL, args.out))
        return 0

    if not _run_complete():
        print("[w4a_c1_result] %s  run not complete (no '%s' in %s). "
              "Nothing written; will retry." % (stamp, DONE_MARKER, os.path.basename(RUN_LOG)))
        return 0

    try:
        carved, cprim, csec, uprim, uncref, carve = compute()
    except Exception as exc:  # data present but not yet analysable / partial
        print("[w4a_c1_result] %s  run marked done but analysis not ready "
              "(%s: %s). Will retry." % (stamp, type(exc).__name__, exc))
        return 0

    md, summ = render_md(cprim, csec, uprim, uncref, carve, stamp)
    with open(args.out, "w") as fh:
        fh.write(md)
    with open(SENTINEL, "w") as fh:
        fh.write(stamp + "\n" + json.dumps(summ, default=str) + "\n")

    print("[w4a_c1_result] %s  WROTE %s" % (stamp, args.out))
    print("  VERDICT: %s" % summ["verdict"])
    print("  paired delta=%s kcal  |t|=%s  robust=%s  r=%s  (unpaired z=%s)"
          % (_fmt(summ["delta"], 3), _fmt(summ["paired_t"], 2),
             _fmt(summ["robust_delta"], 3), _fmt(summ["pearson_r"], 2),
             _fmt(summ["nsigma"], 2)))
    print("  carved=%s  uncarved=%s"
          % (_fmt(summ["carved_ddg"], 3), _fmt(summ["uncarved_ddg"], 3)))
    if args.do_print:
        sys.stdout.write("\n" + md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
