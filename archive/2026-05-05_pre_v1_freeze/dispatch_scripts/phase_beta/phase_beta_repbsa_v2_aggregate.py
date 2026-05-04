#!/usr/bin/env python
"""scripts/phase_beta_repbsa_v2_aggregate.py — Phase β re-PBSA v2 family aggregation.

Reads ``outputs/<sys>/mmpbsa_results_postl387_v2/mmpbsa_summary.json`` for
the v2 sweep + falls back to ``mmpbsa_results_postl387/`` for #84-completed
systems (e.g. 7TL8_MTR6_s7) when computing per-family aggregates.

Computed per-family quantities follow Definition 3 multi-system σ_btwn:
    seed-mean      mu_i  = mean(ΔG_25 of seed i)
    seed-σ_w       sw_i  = stdev(ΔG_25 of seed i, ddof=1)
    family-⟨ΔG⟩    M     = mean over seeds of mu_i
    family-σ_btwn σ_btw  = stdev over seeds of mu_i, ddof=1
    σ_w_med       med    = median over seeds of sw_i
    SE             SE    = σ_btw / sqrt(n_seeds)
    CI95           ±1.96·SE
    z_SE           |M| / SE
    tier           Tier-1 if z_SE >= 3, Tier-2 if 2 <= z_SE < 3, Tier-3 < 2

Output: one JSON + a markdown table fragment.
"""

from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional

REPO = Path(__file__).resolve().parent.parent

# Per-family seed lists (post-#84 + v2 incorporated).
FAMILIES: Dict[str, Dict] = {
    "2QKI_Cp4": {
        "seeds": ["s7", "s19_reseed55", "s23", "s42", "s83", "s101", "s163", "s251"],
        # v2-only family; #84 had 0/8 (CONECT bug).
        "preferred_dir": "mmpbsa_results_postl387_v2",
        "fallback_dir": None,
    },
    "7TL8_MTR6": {
        # s7 done in #84; s19, s42 done in v2.
        "seeds": ["s7", "s19", "s42"],
        "preferred_dir": "mmpbsa_results_postl387_v2",
        "fallback_dir": "mmpbsa_results_postl387",
    },
    # Other families come from #84 unchanged (for cross-comparison).
    "1EBP_MTR13": {
        "seeds": ["s7", "s19", "s23", "s42", "s101"],
        "preferred_dir": "mmpbsa_results_postl387",
        "fallback_dir": None,
    },
    "2QKH_MTR25": {
        "seeds": ["s7", "s19", "s42"],
        "preferred_dir": "mmpbsa_results_postl387",
        "fallback_dir": None,
    },
    "3IOL_NML20": {
        # s101 excluded — real dissociation per #84.
        "seeds": ["s19", "s23", "s42"],
        "preferred_dir": "mmpbsa_results_postl387",
        "fallback_dir": None,
    },
    "1YCR_NML22": {
        # s19/s42 excluded — real detachment per #84.
        "seeds": ["s7", "s23", "s101"],
        "preferred_dir": "mmpbsa_results_postl387",
        "fallback_dir": None,
    },
}


def _load_summary(family: str, seed: str, prefer: str, fallback: Optional[str]) -> Optional[Dict]:
    sysdir = REPO / "outputs" / f"{family}_calib_{seed}"
    for sub in (prefer, fallback) if fallback else (prefer,):
        if sub is None:
            continue
        p = sysdir / sub / "mmpbsa_summary.json"
        if p.exists():
            try:
                d = json.load(open(p))
                d["__source_dir__"] = sub
                return d
            except Exception as exc:  # noqa: BLE001
                print(f"[WARN] {family}_{seed} {sub} parse fail: {exc}", file=sys.stderr)
    return None


def _seed_stats(summary: Dict) -> Dict:
    dgs = [r["delta_g_kcal"] for r in summary.get("results", [])]
    n = len(dgs)
    mu = sum(dgs) / n if n else float("nan")
    sw = statistics.stdev(dgs) if n >= 2 else float("nan")
    md = statistics.median(dgs) if n else float("nan")
    n_pos = sum(1 for x in dgs if x > 0)
    return {
        "n": n, "mu": mu, "sigma_w": sw, "median": md, "n_pos": n_pos,
        "min": min(dgs) if dgs else float("nan"),
        "max": max(dgs) if dgs else float("nan"),
    }


def _family_aggregate(family: str, cfg: Dict) -> Dict:
    rows: List[Dict] = []
    for seed in cfg["seeds"]:
        s = _load_summary(family, seed, cfg["preferred_dir"], cfg["fallback_dir"])
        if s is None:
            rows.append({"seed": seed, "status": "MISSING"})
            continue
        st = _seed_stats(s)
        rows.append({
            "seed": seed,
            "status": "OK",
            "source_dir": s["__source_dir__"],
            **st,
        })
    valid = [r for r in rows if r.get("status") == "OK" and not math.isnan(r.get("mu", float("nan")))]
    n_seeds = len(valid)
    if n_seeds >= 2:
        mus = [r["mu"] for r in valid]
        sws = [r["sigma_w"] for r in valid]
        M = sum(mus) / len(mus)
        sigma_btw = statistics.stdev(mus)
        sigma_w_med = statistics.median(sws)
        SE = sigma_btw / math.sqrt(n_seeds)
        CI95 = 1.96 * SE
        z_SE = abs(M) / SE if SE > 0 else float("inf")
        tier = "Tier-1" if z_SE >= 3 else ("Tier-2" if z_SE >= 2 else "Tier-3")
    else:
        M = sigma_btw = sigma_w_med = SE = CI95 = z_SE = float("nan")
        tier = "INSUFFICIENT_SEEDS"
    return {
        "family": family,
        "n_seeds": n_seeds,
        "rows": rows,
        "M": M,
        "sigma_btwn": sigma_btw,
        "sigma_w_med": sigma_w_med,
        "SE": SE,
        "CI95": CI95,
        "z_SE": z_SE,
        "tier": tier,
    }


def main():
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/phase_beta_repbsa_v2_aggregate.json")

    families = {f: _family_aggregate(f, cfg) for f, cfg in FAMILIES.items()}

    # Print per-family table
    print("Phase β re-PBSA v2 — Per-family Aggregation")
    print("=" * 110)
    print(f"{'Family':<14s} {'n_seeds':>7s} {'⟨ΔG⟩':>10s} {'σ_btwn':>8s} {'σ_w_med':>8s} {'SE':>7s} {'CI95':>7s} {'z_SE':>6s} {'tier':<10s}")
    print("-" * 110)
    for f, agg in families.items():
        if agg["n_seeds"] >= 2:
            print(f"{f:<14s} {agg['n_seeds']:>7d} {agg['M']:>10.3f} {agg['sigma_btwn']:>8.3f} "
                  f"{agg['sigma_w_med']:>8.3f} {agg['SE']:>7.3f} {agg['CI95']:>7.3f} {agg['z_SE']:>6.2f} {agg['tier']:<10s}")
        else:
            print(f"{f:<14s} {agg['n_seeds']:>7d} (insufficient seeds — {len([r for r in agg['rows'] if r.get('status')=='OK'])} valid)")

    print()
    print("Per-seed detail:")
    print("-" * 110)
    for f, agg in families.items():
        for r in agg["rows"]:
            if r.get("status") == "OK":
                print(f"  {f:<14s} {r['seed']:<14s} src={r['source_dir']:<32s} n={r['n']:>3d} mu={r['mu']:>9.3f}  σ_w={r['sigma_w']:>6.3f}  med={r['median']:>9.3f}  n+={r['n_pos']}")
            else:
                print(f"  {f:<14s} {r['seed']:<14s} {r.get('status', 'MISSING')}")

    out_path.write_text(json.dumps({
        "schema": "phase_beta_repbsa_v2_aggregate/0.1",
        "families": families,
    }, indent=2, default=lambda o: float(o) if isinstance(o, float) else o))
    print(f"\nJSON: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
