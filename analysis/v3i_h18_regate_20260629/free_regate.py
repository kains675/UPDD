#!/usr/bin/env python
"""H18-grade sampling-validity re-gate for the EXISTING V3I FREE leg.

READ-ONLY. No simulation re-run. Applies the H18 nonparametric sampling-validity
gate (hist-OVL floor 0.45) to the already-completed n=6 V3I FREE-leg .out files to
decide whether the existing 12-state/400cyc free leg can be REUSED (paired with the
new bound re-run) or must ALSO be re-run for a valid paired ddG.

This is the FREE-leg twin of regate.py (which re-gated the bound leg). It reuses the
EXACT same overlap kernels + col0-grouping logic; only the out-root and the
seed->rep resolution differ (the free leg's 6 seeds are co-located under one root).

OVERLAP COMPUTATION reused VERBATIM from regate.py / reprobe_gate.py:
  - hist_ovl(), gauss_bhatt(), roundtrips(), and the load() grouping logic
    (pool col9=pertE by col0=state across ALL walker files; dminus state-flip
     NS-1-s so adjacent pairs are physically aligned).

CRITICAL correctness rules enforced (R-18 / project memory):
  - trackb .out columns 0-indexed: col0=state, col3=lambda1, col4=lambda2, col9=pertE.
  - GROUP BY col0=state, NEVER by directory/rep-index. A given rN dir is a WALKER,
    not a ladder-state -> directory grouping fabricates phantom overlaps.
  - hist-OVL is the truer metric; gauss-Bhatt over-estimates for left-skewed
    leg-up distributions (reported for context only, NOT the gate).
  - We do NOT use ADR-0014's OLD weak free gate (overlap >= 0.10, "11/11 pass").
    Floor is hist-OVL >= 0.45 (the H18 floor applied to the bound leg).

Layout (verified): <root>/wt/free/repN/{dplus,dminus}/r{0..NS-1}/trackb_{dplus,dminus}.out
  NS auto-detected from #walker files (=12 here). n_cycles=400.
  (A separate non-directional <root>/wt/free/repN/rN/trackb.out set also exists; it is
   NOT the ADR-0014 directional free ladder and is deliberately ignored.)
"""
import os, sys, glob, json
import numpy as np

OVL_FLOOR = 0.45          # H18 sampling-validity floor (nonparametric hist-OVL)
BHATT_FLOOR = 0.80        # context only (reprobe_gate dual-floor); NOT the H18 gate
OUTDIR = "/home/san/UPDD_proj/analysis/v3i_h18_regate_20260629"
FREE_ROOT = "/home/san/UPDD_proj/outputs/_trackb/twocopy_v3i_free_n6"

# seed -> rep dir. The free leg's 6 reps are co-located under one out-root
# (pool_manifest_free.json: rep0=s7, rep1=s101, rep2=s127, rep3=s23, rep4=s163, rep5=s199).
SEED_PATHS = {
    "s7":   f"{FREE_ROOT}/wt/free/rep0",
    "s101": f"{FREE_ROOT}/wt/free/rep1",
    "s127": f"{FREE_ROOT}/wt/free/rep2",
    "s23":  f"{FREE_ROOT}/wt/free/rep3",
    "s163": f"{FREE_ROOT}/wt/free/rep4",
    "s199": f"{FREE_ROOT}/wt/free/rep5",
}
DIRECTIONS = ("dplus", "dminus")


# ---- overlap kernels reused VERBATIM from regate.py ------------------------
def hist_ovl(x, y, nb=60):
    if len(x) < 5 or len(y) < 5:
        return np.nan
    lo = min(x.min(), y.min()); hi = max(x.max(), y.max())
    if hi - lo < 1e-9:
        return 1.0
    b = np.linspace(lo, hi, nb + 1)
    px, _ = np.histogram(x, bins=b, density=True)
    py, _ = np.histogram(y, bins=b, density=True)
    bw = b[1] - b[0]
    return float(np.sum(np.minimum(px, py)) * bw)

def gauss_bhatt(x, y):
    if len(x) < 3 or len(y) < 3:
        return np.nan
    m1, s1 = x.mean(), x.std(); m2, s2 = y.mean(), y.std()
    if s1 < 1e-9 or s2 < 1e-9:
        return np.nan
    return float(np.sqrt(2 * s1 * s2 / (s1**2 + s2**2)) * np.exp(-(m1 - m2)**2 / (4 * (s1**2 + s2**2))))


# ---- loaders: GROUP BY col0=state (regate.py load() logic) -----------------
def detect_NS(rep_dir, d):
    return len(glob.glob(f"{rep_dir}/{d}/r*/trackb_{d}.out"))

def load(rep_dir, d, NS):
    """Pool col9=pertE by col0=state across ALL walker files (rN = walker, not state).
    dminus states are flipped NS-1-s so adjacent pairs align physically (reprobe_gate)."""
    W = {}                                   # walker index -> col0 state trace
    pe = {s: [] for s in range(NS)}          # physical-state -> pooled pertE
    for w in range(NS):
        f = f"{rep_dir}/{d}/r{w}/trackb_{d}.out"
        if not os.path.exists(f):
            continue
        a = np.loadtxt(f)
        a = a[None, :] if a.ndim == 1 else a
        W[w] = a[:, 0].astype(int)
        for s, p in zip(a[:, 0].astype(int), a[:, 9]):
            ps = s if d == "dplus" else NS - 1 - s
            if 0 <= ps < NS:
                pe[ps].append(p)
    return W, {s: np.array(v) for s, v in pe.items()}

def roundtrips(W, NS):
    """state0 <-> state_max boundary crossings, counted on col0 walker traces."""
    tot = 0
    for st in W.values():
        hits = []
        for s in st:
            if s == 0 or s == NS - 1:
                if not hits or hits[-1] != s:
                    hits.append(s)
        tot += sum(1 for i in range(1, len(hits)) if hits[i] != hits[i - 1]) // 2
    return tot


def gate_one_leg(seed, rep_dir, d):
    NS = detect_NS(rep_dir, d)
    if NS == 0:
        return None
    W, pe = load(rep_dir, d, NS)
    pairs = []
    for j in range(NS - 1):
        o = hist_ovl(pe.get(j, np.array([])), pe.get(j + 1, np.array([])))
        b = gauss_bhatt(pe.get(j, np.array([])), pe.get(j + 1, np.array([])))
        pairs.append({"pair": [j, j + 1], "hist_ovl": o, "gauss_bhatt": b,
                      "n_lo": int(len(pe.get(j, []))), "n_hi": int(len(pe.get(j + 1, [])))})
    valid = [p for p in pairs if not np.isnan(p["hist_ovl"])]
    ovls = [p["hist_ovl"] for p in valid]
    mn = min(ovls) if ovls else float("nan")
    worst = min(valid, key=lambda p: p["hist_ovl"]) if valid else None
    n_below = sum(1 for o in ovls if o < OVL_FLOOR)
    rt = roundtrips(W, NS)
    return {
        "seed": seed, "direction": d, "NS": NS, "n_pairs": len(valid),
        "min_hist_ovl": round(mn, 4) if not np.isnan(mn) else None,
        "worst_pair": f"({worst['pair'][0]},{worst['pair'][1]})" if worst else None,
        "n_pairs_below_045": n_below,
        "pass_045": bool(ovls) and (mn >= OVL_FLOOR),
        "round_trips": int(rt),
        "all_pairs": [{"pair": f"({p['pair'][0]},{p['pair'][1]})",
                       "hist_ovl": round(p["hist_ovl"], 4) if not np.isnan(p["hist_ovl"]) else None,
                       "gauss_bhatt": round(p["gauss_bhatt"], 4) if not np.isnan(p["gauss_bhatt"]) else None,
                       "below_045": (not np.isnan(p["hist_ovl"])) and p["hist_ovl"] < OVL_FLOOR}
                      for p in pairs],
    }


def main():
    print("=" * 90)
    print("V3I FREE-leg H18-grade sampling-validity re-gate (nonparametric hist-OVL)")
    print(f"floor: hist-OVL >= {OVL_FLOOR} for ALL adjacent state pairs (group by col0=state)")
    print(f"data: existing n=6 V3I FREE, lambda1_rampdown 12-state ladder, 400 cyc (ADR-0014)")
    print(f"NOTE: ADR-0014 reported 11/11 >= 0.10 on the OLD weak gate; this applies the NEW 0.45 floor")
    print("=" * 90)
    results = []
    for seed, rep_dir in SEED_PATHS.items():
        for d in DIRECTIONS:
            r = gate_one_leg(seed, rep_dir, d)
            if r is None:
                print(f"\n### {seed}/{d}: NO DATA")
                continue
            results.append(r)
            mn = r["min_hist_ovl"]
            print(f"\n### {seed} / {d}  (NS={r['NS']}, pairs={r['n_pairs']})")
            print(f"  min hist-OVL = {mn} @ {r['worst_pair']}   "
                  f"pairs<{OVL_FLOOR} = {r['n_pairs_below_045']}   "
                  f"round_trips(col0) = {r['round_trips']}   "
                  f"{'PASS' if r['pass_045'] else 'FAIL'}")
            below = [p for p in r["all_pairs"] if p["below_045"]]
            if below:
                print("   below-floor pairs: " +
                      ", ".join(f"{p['pair']}={p['hist_ovl']}" for p in below))

    n_total = len(results)
    n_pass = sum(1 for r in results if r["pass_045"])
    # verdict vocabulary per task: reuse_ok (all pass) / free_needs_rerun (any fail) / mixed
    if n_total == 0:
        verdict = "free_needs_rerun"   # no data == cannot reuse
    elif n_pass == n_total:
        verdict = "reuse_ok"
    elif n_pass == 0:
        verdict = "free_needs_rerun"
    else:
        verdict = "mixed"

    print("\n" + "=" * 90)
    print(f"n_legs_pass_045 / n_legs_total = {n_pass} / {n_total}")
    print(f"OVERALL VERDICT: {verdict}")
    if verdict != "reuse_ok":
        print("CAVEAT: if the free leg fails the 0.45 floor, the re-run MUST include the free")
        print("        leg too (cheap, ~0.9h/seed) so the paired ddG_bind uses a matched-grade free.")
    print("=" * 90)

    out = {
        "schema": "v3i_h18_free_regate_v1",
        "date": "2026-06-29",
        "leg": "free",
        "out_root": FREE_ROOT,
        "method": ("nonparametric histogram overlap (hist_ovl reused verbatim from regate.py / "
                   "apex_probe_gate_20260628/reprobe_gate.py); col9=pertE pooled by col0=state "
                   "across all walker files (dminus state-flip NS-1-s); floor hist-OVL>=0.45; "
                   "gauss-Bhatt reported context-only (overestimates left-skew); "
                   "ADR-0014's old 0.10 free gate NOT used"),
        "data_caveat": ("existing free run = lambda1_rampdown 12-state/400cyc (ADR-0014); ADR-0014 "
                        "reported 11/11 >= 0.10 on the OLD weak gate. A PASS here means the free "
                        "spacing already meets the H18 0.45 floor and the free leg can be REUSED "
                        "paired with the new bound re-run. A FAIL means the free leg must ALSO be "
                        "re-run (cheap, ~0.9h/seed) for a valid paired ddG_bind."),
        "ovl_floor": OVL_FLOOR,
        "n_legs_total": n_total,
        "n_legs_pass_045": n_pass,
        "verdict": verdict,
        "per_leg": results,
    }
    with open(os.path.join(OUTDIR, "free_regate_result.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {os.path.join(OUTDIR, 'free_regate_result.json')}")


if __name__ == "__main__":
    main()
