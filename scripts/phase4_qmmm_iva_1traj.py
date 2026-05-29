#!/usr/bin/env python3
"""Phase IV Track A — 1-trajectory QM ΔΔ_int (swap-in-place) orchestrator.

Mirrors scripts/phase4_qmmm_iva.py, but replaces the two-ensemble paired design
(which produced the #84 NULL, ΔΔ_int +24.6 ± 55.1, p=0.25) with a single-
trajectory same-coordinate evaluation: on each MD frame, residue-4's sidechain
is swapped in place between Trp (WT) and 1-Me-Trp (Cp4), so the ~−190 kcal/mol
binder/target self-energy cancels EXACTLY per frame (Genheden & Ryde 2015;
Hou 2011 DOI:10.1002/jcc.21666). Both endpoint identities are evaluated by the
SAME utils/run_qmmm.py engine (never modified) on byte-identical non-residue-4
atoms, so the QM region / basis / df-mode / link-atom set are identical (A-C4).

A-C1 (both directions, mandatory):
  - graft direction: base = WT MD frames; graft methyl → Cp4; ΔΔ_int = Cp4 − WT.
  - strip direction: base = Cp4 MD frames; strip methyl → WT; ΔΔ_int = Cp4 − WT.
  Per-frame ΔΔ_int is aggregated separately for each direction; the two means
  must agree in SIGN and bracket the true value. Divergence = 1-traj breakdown
  (non-local perturbation) → escalate (skip to Track B).

A-C2: clash filter MANDATORY (reject new heavy-contact <2.0 Å / H-contact
<1.5 Å on the swapped endpoint), THEN (follow-up) constrained micro-min. This
orchestrator records the no-min clash decision per frame; the with-min column is
the Keeper-gated follow-up (see utils.qmmm_1traj_variant_compare.constrained_micro_min).

A-C3: gas-phase qm_int_kcal_frozen is the primary column; the COSMO/PCM (ε=78.5)
second column is produced on the SAME geometries when --pcm is set (the engine-
reuse PCM hook is wired but unverified — no sign claim rests on gas-phase alone).

Charge axis (cross-track / A-C1): Option-β / regime-2 hybrid MTR FF XML
(params/MTR_gaff2_hybrid.xml, Σq=0, NE1 frozen −0.3418). MTR and Trp both carry
formal charge 0, so the swap is charge-neutral and the R-15/16 guards pass
identically for both endpoints (verified per frame, not bypassed).

Regime: ranking-only (R-11). ΔΔ_int ≠ ΔΔG_bind; NO Magotti absolute comparison.
BSSE +3-8 band + 1-traj-approximation + MM-geometry-bias disclosed in the JSON.

Source verdict: .claude/agent-memory/scival/verdict_trackAB_qm1traj_fep_prevalidation_20260529.md
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from typing import Dict, List, Optional, Tuple

# utils/ on sys.path for dispatch + the swap primitive (reuse, no edits)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "utils"))
from utils import dispatch  # noqa: E402

CONDA_BIN = "/home/san/miniconda3/bin/conda"

# ── Locked experiment constants (mirror phase4_qmmm_iva.py §"Method") ──
SCHEMA_VERSION = "phase4_iv_1traj_swap_ddint_v1"
TARGET_ID = "2QKI"
QM_XC = "wb97xd"
QM_BASIS = "6-31G*"            # A-C4: ONE basis, both endpoint identities (matches #84)
BINDER_CHAIN = "B"
TARGET_RESNUM = 4              # binder position 4 (Trp / 1-Me-Trp)
USE_DF = False                # A-C4: exact direct-SCF (df_mode=False), matches #84
# Option-β / regime-2 hybrid MTR charge XML (A-C1 / cross-track). MM/FF axis only;
# the bare-QM frozen observable derives its QM-region charge topologically.
CHARGE_XML = os.path.join(PROJECT_ROOT, "params", "MTR_gaff2_hybrid.xml")
CHARGE_REGIME = "option_beta_regime2_hybrid_sigmaq0"

# MM-PBSA summary that ranks frames (median-Δg representative rule, A-C4).
MMPBSA_SUMMARY_REL = os.path.join("mmpbsa_results_postl387_v2", "mmpbsa_summary.json")

# ΔΔ_int = Cp4 − WT in BOTH directions (positive = methyl destabilizes binding).
# A-C1: graft base cohort = WT ensemble; strip base cohort = Cp4 ensemble.
DIRECTIONS: Dict[str, Dict[str, object]] = {
    "graft": {
        "base_cohort": "wt",
        "dir_tmpl": "2QKI_WT_calib_{seed}",
        "seeds": ["s7", "s19", "s23", "s101", "s127", "s163", "s199", "s251"],
    },
    "strip": {
        "base_cohort": "cp4",
        # A-C1 / ADR-0010: strip base frames come from the Cp4 hybrid regime-2
        # ensemble (Σq=0, NE1 −0.3418) so the geometry FF axis matches the graft
        # WT ensemble (plain amber14SB Trp, regime-invariant). NOT the patch-OFF
        # 2QKI_Cp4_calib_* ensemble (Σq=−0.187), which would break the cross-track
        # FF axis the verdict's A-C1 protects.
        "dir_tmpl": "2QKI_Cp4_hybrid_calib_{seed}",
        "seeds": ["s7", "s19", "s23", "s101", "s127", "s163", "s199", "s251"],
    },
}
COMMON_SEEDS = ["s7", "s19", "s23", "s101", "s127", "s163", "s199", "s251"]


# ──────────────────────────────────────────────────────────────
# Frame resolution — median-Δg representative + neighbors (≥3/seed)
# ──────────────────────────────────────────────────────────────
def _seed_dir(direction: str, seed: str) -> str:
    return os.path.join(PROJECT_ROOT, "outputs",
                        str(DIRECTIONS[direction]["dir_tmpl"]).format(seed=seed))


def resolve_frames(direction: str, seed: str, n_frames: int) -> Dict[str, object]:
    """Return ``n_frames`` base-frame PDBs for a seed, centered on the median-Δg
    snapshot (deterministic, no cherry-pick: median plus nearest neighbors in the
    Δg-sorted order).

    The base cohort is the WT ensemble for graft and the Cp4 ensemble for strip
    (A-C1). Each frame is later swapped in place to the opposite identity.
    """
    sdir = _seed_dir(direction, seed)
    summ_path = os.path.join(sdir, MMPBSA_SUMMARY_REL)
    info: Dict[str, object] = {
        "direction": direction, "seed": seed, "seed_dir": sdir,
        "summary": summ_path, "frames": [], "resolve_ok": False, "error": None,
    }
    if not os.path.isfile(summ_path):
        info["error"] = f"missing MM-PBSA summary: {os.path.relpath(summ_path, PROJECT_ROOT)}"
        return info
    try:
        with open(summ_path, "r", encoding="utf-8") as f:
            summ = json.load(f)
    except Exception as e:  # noqa: BLE001
        info["error"] = f"summary read failed: {type(e).__name__}: {str(e)[:160]}"
        return info
    results = summ.get("results", [])
    if not results:
        info["error"] = "no MM-PBSA results"
        return info

    ordered = sorted(results, key=lambda r: r["delta_g_kcal"])
    n = len(ordered)
    med_idx = (n - 1) // 2
    half = max(1, n_frames) // 2
    lo = max(0, med_idx - half)
    hi = min(n, lo + n_frames)
    lo = max(0, hi - n_frames)            # re-anchor if we hit the top edge
    chosen = ordered[lo:hi]

    frames: List[Dict[str, object]] = []
    for rec in chosen:
        stem = rec["snapshot"]
        pdb_rel = rec.get("pdb_path")
        pdb_abs = (pdb_rel if os.path.isabs(str(pdb_rel))
                   else os.path.join(PROJECT_ROOT, str(pdb_rel)))
        frames.append({
            "snapshot_stem": stem,
            "delta_g_kcal": float(rec["delta_g_kcal"]),
            "pdb_path": pdb_abs,
            "pdb_exists": os.path.isfile(pdb_abs),
        })
    info["frames"] = frames
    info["n_total_snapshots"] = n
    info["median_index"] = med_idx
    info["resolve_ok"] = bool(frames) and all(fr["pdb_exists"] for fr in frames)
    return info


# ──────────────────────────────────────────────────────────────
# Per-frame compute (reuse the swap primitive + run_qmmm engine)
# ──────────────────────────────────────────────────────────────
def compute_frame(direction: str, seed: str, frame: Dict[str, object],
                  run_pcm: bool, keep_workdir: bool) -> Dict[str, object]:
    """Run one frame's two-state pair via the swap primitive. Heavy import is
    local (pulls PySCF/gpu4pyscf only when actually computing)."""
    from qmmm_1traj_variant_compare import compute_qm_int_pair  # noqa: E402
    pair = compute_qm_int_pair(
        pdb_path=str(frame["pdb_path"]),
        target_resnum=TARGET_RESNUM,
        charge_xml=CHARGE_XML,
        direction=direction,
        target_id=TARGET_ID,
        binder_chain=BINDER_CHAIN,
        qm_basis=QM_BASIS,
        qm_xc=QM_XC,
        use_df=USE_DF,
        run_pcm=run_pcm,
        keep_workdir=keep_workdir,
    )
    rec: Dict[str, object] = {
        "direction": direction,
        "seed": seed,
        "snapshot_stem": frame["snapshot_stem"],
        "base_delta_g_kcal": frame["delta_g_kcal"],
        "ddint_kcal": pair.get("ddint_kcal"),
        "e_int_cp4": pair.get("e_int_cp4"),
        "e_int_wt": pair.get("e_int_wt"),
        "clash_flag": pair.get("clash_flag"),
        "charge_audit": pair.get("charge_audit"),
        "pcm": pair.get("pcm"),
        "fail_reason": pair.get("fail_reason"),
        # A-C2: a clashed no-min frame is excluded from the no-min aggregate;
        # the with-min column is the Keeper-gated follow-up.
        "ok": (pair.get("ddint_kcal") is not None and not pair.get("clash_flag")),
        "ok_including_clashed": pair.get("ddint_kcal") is not None,
    }
    return rec


# ──────────────────────────────────────────────────────────────
# Aggregation — paired stats per direction + direction agreement (A-C1)
# ──────────────────────────────────────────────────────────────
def _direction_stats(records: List[Dict[str, object]], direction: str,
                     include_clashed: bool) -> Dict[str, object]:
    import numpy as np
    from scipy import stats
    key = "ok_including_clashed" if include_clashed else "ok"
    vals = [float(r["ddint_kcal"]) for r in records
            if r["direction"] == direction and r.get(key)]
    n = len(vals)
    out: Dict[str, object] = {
        "direction": direction,
        "n_frames": n,
        "include_clashed": include_clashed,
        "per_frame_ddint_kcal": [round(v, 4) for v in vals],
    }
    if n >= 2:
        arr = np.asarray(vals, dtype=float)
        mean = float(arr.mean())
        sd = float(arr.std(ddof=1))
        se = sd / (n ** 0.5)
        df = n - 1
        tcrit = float(stats.t.ppf(0.975, df))
        tt = stats.ttest_1samp(arr, 0.0)
        out.update({
            "mean_ddint_kcal": round(mean, 4),
            "sd_kcal": round(sd, 4),
            "se_kcal": round(se, 4),
            "df": df,
            "ci95_kcal": [round(mean - tcrit * se, 4), round(mean + tcrit * se, 4)],
            "z_se": round(abs(mean) / se, 4) if se > 0 else None,
            "t_stat": round(float(tt.statistic), 4),
            "p_value": float(tt.pvalue),
        })
    elif n == 1:
        out["mean_ddint_kcal"] = round(vals[0], 4)
        out["note"] = "single frame — no dispersion statistics"
    else:
        out["note"] = "no usable frames"
    return out


def _direction_agreement(graft: Dict[str, object], strip: Dict[str, object]
                         ) -> Dict[str, object]:
    """A-C1 gate: the two directions must agree in SIGN and bracket the true
    value. graft (build methyl onto WT frames) is biased toward spurious
    repulsion (more positive ΔΔ_int); strip (remove methyl on Cp4 frames) toward
    under-binding (more negative). Sign agreement + a finite bracket = 1-traj
    validity; sign divergence = non-local perturbation → escalate.
    """
    mg = graft.get("mean_ddint_kcal")
    ms = strip.get("mean_ddint_kcal")
    res: Dict[str, object] = {
        "graft_mean_ddint_kcal": mg,
        "strip_mean_ddint_kcal": ms,
    }
    if mg is None or ms is None:
        res["verdict"] = "indeterminate"
        res["reason"] = "one or both directions lack a mean (insufficient frames)"
        return res
    mg_f, ms_f = float(mg), float(ms)
    same_sign = (mg_f >= 0) == (ms_f >= 0)
    res["bracket_kcal"] = [round(min(mg_f, ms_f), 4), round(max(mg_f, ms_f), 4)]
    res["bracket_width_kcal"] = round(abs(mg_f - ms_f), 4)
    res["sign_agreement"] = bool(same_sign)
    if same_sign:
        res["verdict"] = "agree_1traj_valid"
        res["reason"] = "both directions agree in sign; true ΔΔ_int lies in the bracket"
    else:
        res["verdict"] = "diverge_escalate"
        res["reason"] = (
            "directions disagree in sign — 1-trajectory approximation breaks down "
            "(non-local perturbation). Per verdict A-C1, escalate to Track B."
        )
    return res


def aggregate(records: List[Dict[str, object]], logdir: str) -> Dict[str, object]:
    graft_nomin = _direction_stats(records, "graft", include_clashed=False)
    strip_nomin = _direction_stats(records, "strip", include_clashed=False)
    graft_all = _direction_stats(records, "graft", include_clashed=True)
    strip_all = _direction_stats(records, "strip", include_clashed=True)
    payload: Dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "framing": "ranking_only_R11",
        "target_id": TARGET_ID,
        "binder_chain": BINDER_CHAIN,
        "swapped_resnum": TARGET_RESNUM,
        "qm_method": f"{QM_XC}/{QM_BASIS}",
        "scf": "exact_direct_no_df",
        "charge_xml": os.path.relpath(CHARGE_XML, PROJECT_ROOT),
        "charge_regime": CHARGE_REGIME,
        "observable": "qm_int_kcal_frozen (1-traj swap-in-place)",
        "observable_definition": (
            "per-frame ΔΔ_int = qm_int_kcal_frozen(Cp4) − qm_int_kcal_frozen(WT), "
            "BOTH evaluated on the SAME swapped MD frame "
            "(E_qm − E_binder_iso − E_target_iso, run_qmmm.py:2516)"
        ),
        "ddint_sign_convention": "Cp4 − WT (positive = methyl destabilizes binding)",
        "caveats": {
            "regime": "ranking-only (R-11); ΔΔ_int ≠ ΔΔG_bind; NO Magotti absolute comparison",
            "bsse": "uncorrected (+3~8 kcal/mol band)",
            "one_traj_approximation": (
                "same-coordinate frozen evaluation; swapped endpoint not "
                "independently relaxed. graft biased to repulsion, strip to "
                "under-binding; the two-direction bracket is the honest interval (A-C1)."
            ),
            "mm_geometry_bias": (
                "frames come from an MM MD ensemble (amber14SB + Option-β hybrid), "
                "not QM-optimized geometries (Senn & Thiel 2009 DOI:10.1002/anie.200802019)."
            ),
            "clash_filter": (
                "no-min aggregate excludes frames whose swap created a new "
                "heavy-contact <2.0 Å or H-contact <1.5 Å (A-C2); the constrained "
                "micro-min (with-min) column is the Keeper-gated follow-up."
            ),
        },
        "per_frame": records,
        "directions": {
            "graft": {"no_min": graft_nomin, "all_frames": graft_all},
            "strip": {"no_min": strip_nomin, "all_frames": strip_all},
        },
        "direction_agreement_no_min": _direction_agreement(graft_nomin, strip_nomin),
    }
    out_path = os.path.join(logdir, "phase4_iva_1traj_ddint.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    return {"path": out_path, "payload": payload}


# ──────────────────────────────────────────────────────────────
# Plan enumeration + logdir
# ──────────────────────────────────────────────────────────────
def enumerate_plan() -> List[Tuple[str, str]]:
    plan: List[Tuple[str, str]] = []
    for direction, cfg in DIRECTIONS.items():
        for seed in cfg["seeds"]:  # type: ignore[union-attr]
            plan.append((direction, seed))
    return plan


def _new_logdir() -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    d = os.path.join(PROJECT_ROOT, "outputs", "analysis", f"phase4_qmmm_iva_1traj_{ts}")
    os.makedirs(d, exist_ok=True)
    return d


# ──────────────────────────────────────────────────────────────
# CLI modes
# ──────────────────────────────────────────────────────────────
def mode_dry_run(n_frames: int) -> int:
    print("=" * 78)
    print("Phase IV Track A — 1-traj QM ΔΔ_int (swap-in-place) — DRY RUN (no execution)")
    print(f"  directions: graft (base=WT frames → build methyl), "
          f"strip (base=Cp4 frames → remove methyl)   [A-C1]")
    print(f"  frames/seed: {n_frames} (median-Δg ± neighbors, deterministic)")
    print(f"  common seeds (n={len(COMMON_SEEDS)}): {COMMON_SEEDS}")
    print(f"  engine lock: {QM_XC}/{QM_BASIS}  df_mode={USE_DF}  "
          f"--target-id {TARGET_ID}  --binder_chain {BINDER_CHAIN}  resnum={TARGET_RESNUM}")
    print(f"  charge axis: {os.path.relpath(CHARGE_XML, PROJECT_ROOT)} "
          f"({CHARGE_REGIME}); exists={os.path.isfile(CHARGE_XML)}")
    print("=" * 78)
    missing = 0
    total_frames = 0
    for direction, seed in enumerate_plan():
        info = resolve_frames(direction, seed, n_frames)
        if not info["resolve_ok"]:
            print(f"  [{direction:5s} {seed:6s}] RESOLVE-FAIL: {info.get('error') or 'frame PDB missing'}")
            missing += 1
            continue
        frames = info["frames"]  # type: ignore[assignment]
        total_frames += len(frames)
        stems = ", ".join(str(fr["snapshot_stem"]).split("_", 1)[-1] for fr in frames)  # type: ignore[union-attr]
        print(f"  [{direction:5s} {seed:6s}] {len(frames)} frames "
              f"(median idx {info.get('median_index')}/{info.get('n_total_snapshots')}): {stems}")
    print("-" * 78)
    print("Aggregation targets:")
    print("  per-direction paired stats (no-min + all-frames), direction-agreement gate (A-C1)")
    print(f"  output: <logdir>/phase4_iva_1traj_ddint.json (schema={SCHEMA_VERSION})")
    print("=" * 78)
    if missing:
        print(f"[!] {missing} seed/direction(s) unresolved — investigate before launch.")
        return 1
    print(f"[OK] all {len(enumerate_plan())} seed/direction plans resolved "
          f"({total_frames} frames total); plan ready for Keeper/Runner gate.")
    return 0


def mode_smoke(direction: str, seed: str, n_frames: int, run_pcm: bool) -> int:
    if direction not in DIRECTIONS:
        print(f"[!] unknown direction '{direction}' (choose: {list(DIRECTIONS)})")
        return 2
    if seed not in DIRECTIONS[direction]["seeds"]:  # type: ignore[operator]
        print(f"[!] seed '{seed}' not in {direction} seed list")
        return 2
    logdir = _new_logdir()
    print(f"[SMOKE] logdir={logdir}  direction={direction} seed={seed} frames={n_frames}")
    info = resolve_frames(direction, seed, n_frames)
    if not info["resolve_ok"]:
        print(f"[!] resolve failed: {info.get('error')}")
        return 1
    records: List[Dict[str, object]] = []
    for frame in info["frames"]:  # type: ignore[union-attr]
        rec = compute_frame(direction, seed, frame, run_pcm=run_pcm, keep_workdir=False)
        records.append(rec)
        print(json.dumps(rec, indent=2, ensure_ascii=False))
    agg = aggregate(records, logdir)
    print(f"  → {agg['path']}")
    return 0 if any(r["ok"] for r in records) else 1


def mode_full(n_frames: int, run_pcm: bool) -> int:
    logdir = _new_logdir()
    print(f"[FULL] logdir={logdir}  (sequential; V100 one-at-a-time per UPDD_VM_ENABLE)")
    # Surface the routing decision (the engine itself does the V100 dispatch via
    # dispatch.route_stage('qmmm') inside run_qmmm; here we report the resolved
    # location so the operator sees where SCF will run).
    loc = dispatch.route_stage("qmmm")
    print(f"  qmmm route: {loc} (UPDD_VM_ENABLE={'1' if dispatch.UPDD_VM_ENABLE else '0'})")
    records: List[Dict[str, object]] = []
    for direction, seed in enumerate_plan():
        info = resolve_frames(direction, seed, n_frames)
        if not info["resolve_ok"]:
            print(f"  ── {direction} {seed} ── RESOLVE-FAIL: {info.get('error')}")
            continue
        for frame in info["frames"]:  # type: ignore[union-attr]
            print(f"  ── {direction} {seed} {frame['snapshot_stem']} ──")
            rec = compute_frame(direction, seed, frame, run_pcm=run_pcm, keep_workdir=False)
            tag = "ok" if rec["ok"] else (
                "CLASH-REJECT" if (rec["ok_including_clashed"] and rec["clash_flag"])
                else f"FAIL({rec['fail_reason']})")
            print(f"     {tag}  ddint={rec['ddint_kcal']} "
                  f"(cp4={rec['e_int_cp4']} wt={rec['e_int_wt']} clash={rec['clash_flag']})")
            records.append(rec)
    agg = aggregate(records, logdir)
    p = agg["payload"]
    print("=" * 78)
    for direction in ("graft", "strip"):
        d = p["directions"][direction]["no_min"]  # type: ignore[index]
        if "mean_ddint_kcal" in d:
            print(f"  {direction:5s} (no-min) ΔΔ_int = {d['mean_ddint_kcal']} "
                  f"± {d.get('sd_kcal')} (SE {d.get('se_kcal')}, z_SE {d.get('z_se')}, "
                  f"CI95 {d.get('ci95_kcal')}, n={d['n_frames']})")
        else:
            print(f"  {direction:5s} (no-min): {d.get('note')}")
    ag = p["direction_agreement_no_min"]
    print(f"  direction agreement: {ag.get('verdict')} — {ag.get('reason')}")
    print(f"  → {agg['path']}")
    print("=" * 78)
    n_ok = sum(1 for r in records if r["ok"])
    return 0 if (records and n_ok > 0) else 1


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Phase IV Track A 1-traj QM ΔΔ_int (swap-in-place) orchestrator (ranking-only, R-11)."
    )
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve frames + print plan; no execution")
    ap.add_argument("--smoke", nargs=2, metavar=("DIRECTION", "SEED"),
                    help="run one seed of one direction (V100 dispatch path check)")
    ap.add_argument("--frames", type=int, default=3,
                    help="frames per seed (median-Δg ± neighbors; default 3, A-C4 ≥3)")
    ap.add_argument("--pcm", action="store_true",
                    help="also compute the COSMO/PCM (ε=78.5) second column (A-C3); "
                         "PCM hook is wired but unverified — see module docstring")
    args = ap.parse_args()

    if args.dry_run:
        return mode_dry_run(args.frames)
    if args.smoke:
        return mode_smoke(args.smoke[0], args.smoke[1], args.frames, args.pcm)
    return mode_full(args.frames, args.pcm)


if __name__ == "__main__":
    sys.exit(main())
