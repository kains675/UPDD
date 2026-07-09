#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B — barnase Ile96->Ala FOLDING-thermocycle production orchestrator.

Gate-A LARGE-effect engine control (governance
analysis/ncaa_funnel_gate_governance_20260702.md §Phase-1; scientific review GO-WITH-COND
verdict_gateA_largeeffect_control_design_20260702). Runs a canonical, single-
scaffold Ile->Ala side-chain morph through the IDENTICAL two-copy ATS in-place
RBFE engine used by the ncAA vs-WT pairs (the FROZEN soft-core canon, the H18
15-state ladder, and the UWHAM estimator are ALL reused verbatim — this
orchestrator is strictly system-prep + dispatch + a LEG-axis subtraction).

THERMOCYCLE (folding-stability double difference):
  ddG_fold(Ile->Ala) = ddG_mut(folded barnase) - ddG_mut(Ac-Ile-NMe)
                     = dgbind1_FOLDED - dgbind1_TRIPEP
Both legs use the CANONICAL all-amber path (is_ncaa_mtr=False) and the FREE leg
engine mechanics (binder-only two-copy — barnase folded is a single chain with no
receptor/binder split, so the whole molecule IS the duplicated 'binder'). The
number is the LEG-AXIS per-seed difference (paired_difference_stats over
dgbind1_folded[s] - dgbind1_tripep[s]) — NOT analyze_leg_cohort's cp4-wt endpoint
subtraction (which is DEGENERATE for a single-scaffold mutation: the builder never
receives a distinct cp4/wt endpoint, so that axis reports ~0 = a false null).

Cohorts (dispatched via the production run_leg over MATCHED velocity seeds):
  FOLDED : scaffold = an equilibrated H-complete single-chain barnase folded
           final.pdb ; endpoint label "folded" ; the minuend.
  TRIPEP : scaffold = an Ac-Ile-NMe (ACE-ILE-NME) equilibrated final.pdb ;
           endpoint label "tripep" ; the subtrahend (the ff14SB FF-bias probe).

Two MutationSpec INSTANCES (Ile->Ala, identical atom partition / soft-core,
different resnum) are constructed here (make_ile_ala_mutation_spec) and passed as
objects — NO registry edit. auto_search_displacement=True is REQUIRED (the 4 nm
default cannot separate two barnase copies; the search escalates the magnitude).

Sign is claimed ONLY when the paired |z| >= 3 AND n >= 3 (ranking-only, R-11/R-18);
the anchor exp ddG_fold(I96A) ~ 3.3-4.0 kcal/mol (destabilizing) is ABOVE the
detection floor. This orchestrator does the DISPATCH + the LEG-axis combine; the
scaffold MD prep is a separate DATA step (out of this code's scope).

Ranking-only (R-11); NOT a converged calibrated ddG. DOI references:
  - Gallicchio 2025 J Chem Inf Model, ATS DOI 10.1021/acs.jcim.5c00207
  - anchor exp ddG_fold(I96A): Kellis/Nyberg/Sali/Fersht 1988 Nature 333:784
    (+4.0 kcal, DOI 10.1038/333784a0) + Serrano/Kellis/.../Fersht 1992 JMB 224:783
    (+3.3 kcal). NB I96V = +1.1 kcal (single -CH2-, SUB-FLOOR) is a DIFFERENT mutation.
  - amber-family FEP baseline (I96A = 3.9 kcal, agrees with exp; shows the FF is NOT
    the bottleneck for this aliphatic cavity): Sun/Veenstra/Kollman 1996 Protein Eng
    9(3):273, DOI 10.1093/protein/9.3.273.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from typing import Any, Dict, List, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)
if _PROJ not in sys.path:
    sys.path.insert(0, _PROJ)


def _load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _load_prod():
    return _load_module(
        "trackb_inplace_rbfe_production",
        os.path.join(_HERE, "trackb_inplace_rbfe_production.py"))


def _load_ats():
    return _load_module(
        "atm_trackB_setup", os.path.join(_PROJ, "utils", "atm_trackB_setup.py"))


def _find_core_ile_resnum(pdb_path: str, prev: str = "LEU",
                          nxt: str = "TYR") -> int:
    """Locate the buried core Ile by the unique ``prev``-ILE-``nxt`` motif
    (barnase = Leu95-Ile96-Tyr97), robust to any residue renumbering in prep."""
    from openmm.app import PDBFile
    rl = [(int(r.id), r.name) for r in PDBFile(pdb_path).topology.residues()]
    hits = [rl[i][0] for i in range(1, len(rl) - 1)
            if rl[i][1] == "ILE" and rl[i - 1][1] == prev and rl[i + 1][1] == nxt]
    if not hits:
        raise RuntimeError(
            "_find_core_ile_resnum: no %s-ILE-%s motif in %s"
            % (prev, nxt, pdb_path))
    return hits[0]


def _dgbind1_vector(uwham, out_root: str, endpoint: str, leg: str,
                    n_replicates: int, mintimeid: Optional[int]) -> Dict[str, Any]:
    """Per-seed dgbind1 for one cohort via analyze_replicate_set (VERBATIM)."""
    prod = _load_prod()
    leg_dirs = [prod._rep_dir(out_root, endpoint, leg, j)
                for j in range(n_replicates)]
    res = uwham.analyze_replicate_set(
        leg_dirs, jobname=prod.JOBNAME, mintimeid=mintimeid, maxtimeid=None,
        require_cntl_schedule=True)
    dgbind1 = [float(r["dgbind1_kcal"]) for r in res["replicates"]]
    return {
        "endpoint": endpoint,
        "leg_dirs": leg_dirs,
        "n_replicates": res["n_replicates"],
        "mean_dgbind1_kcal": res["mean_dgbind1_kcal"],
        "sigma_btwn_dgbind1_kcal": res["sigma_btwn_dgbind1_kcal"],
        "sem_dgbind1_kcal": res["sem_dgbind1_kcal"],
        "per_seed_dgbind1_kcal": dgbind1,
        "frozen_plateau_guard": res.get("frozen_plateau_guard"),
    }


def analyze_fold_thermocycle(out_root: str, *, leg: str, n_replicates: int,
                             mintimeid: Optional[int],
                             output_json: Optional[str] = None) -> Dict[str, Any]:
    """LEG-axis folding double-difference: ddG_fold[s] = dgbind1_folded[s] -
    dgbind1_tripep[s], combined by paired_difference_stats (VERBATIM).

    Deliberately does NOT use analyze_leg_cohort (cp4-vs-wt endpoint subtraction),
    which is degenerate for a single-scaffold mutation (endpoint-axis trap).
    """
    prod = _load_prod()
    uwham = prod._load_uwham()

    folded = _dgbind1_vector(uwham, out_root, "folded", leg, n_replicates, mintimeid)
    tripep = _dgbind1_vector(uwham, out_root, "tripep", leg, n_replicates, mintimeid)

    a = folded["per_seed_dgbind1_kcal"]   # minuend (folded)
    b = tripep["per_seed_dgbind1_kcal"]   # subtrahend (tripeptide reference)
    n_matched = min(len(a), len(b))
    a, b = a[:n_matched], b[:n_matched]
    paired = uwham.paired_difference_stats(a, b)

    # Sign gate (ranking-only): the paired-t sign_status AND n >= MIN_REPLICATES.
    sign_status = paired.get("sign_status", "undetermined")
    sign_reason: Optional[str] = None
    if n_matched < prod.MIN_REPLICATES_FOR_SIGN:
        sign_status = "undetermined"
        sign_reason = ("n=%d < %d matched seeds — sign not claimed (R-18)"
                       % (n_matched, prod.MIN_REPLICATES_FOR_SIGN))
    else:
        t = paired.get("t_stat")
        if t is None or abs(t) < prod.Z_SE_SIGN_THRESHOLD:
            sign_status = "undetermined"
            sign_reason = ("paired |t|=%s < %.0f — sign not claimed"
                           % (("%.2f" % t) if t is not None else "n/a",
                              prod.Z_SE_SIGN_THRESHOLD))
        else:
            sign_reason = ("paired |t|=%.2f >= %.0f AND n=%d >= %d — sign %s"
                           % (t, prod.Z_SE_SIGN_THRESHOLD, n_matched,
                              prod.MIN_REPLICATES_FOR_SIGN, sign_status))

    payload = {
        "protocol": "fold_thermocycle_leg_axis",
        "leg": leg,
        "axis": "ddG_fold[s] = dgbind1_folded[s] - dgbind1_tripep[s]",
        "n_matched_seeds": n_matched,
        "folded": folded,
        "tripep": tripep,
        "ddG_fold_kcal": paired.get("mean"),
        "paired": paired,
        "sign_status": sign_status,
        "sign_reason": sign_reason,
        "regime": "ranking_only",
        "note": ("LEG-axis folding double difference (NOT the degenerate cp4-wt "
                 "endpoint subtraction). Convention: positive ddG_fold = "
                 "DESTABILIZING (matches Fersht ddG_(D-N)); a mutated barnase "
                 "with dgbind1 more positive than the tripeptide reference is "
                 "destabilized. Confirm the sign convention against the anchor "
                 "before any wet-lab claim (R-18)."),
    }
    if output_json:
        with open(output_json, "w") as fh:
            json.dump(payload, fh, indent=2)
    return payload


def _dispatch_cohort(prod, ats, *, out_root: str, endpoint: str,
                     scaffold_final: str, resnum: int, binder_chain: str,
                     seeds: List[str], directions: List[str], args) -> Dict[str, Any]:
    """Dispatch one folding cohort (FOLDED or TRIPEP) through the production
    run_leg (leg='free' two-copy, auto_search_displacement, per-scaffold spec +
    leg_inputs override)."""
    spec = ats.make_ile_ala_mutation_spec(resnum, name="ile_ala_%s" % endpoint)
    if spec.shape != "acyclic_connected_group":
        raise RuntimeError(
            "cohort %s: Ile->Ala spec must classify acyclic_connected_group, got %r"
            % (endpoint, spec.shape))
    leg_inputs = ats.resolve_fold_leg_inputs(scaffold_final)
    # P3-#116 FIX2: deterministic + bounded-retry appearing-H placement (0 => DISABLE,
    # legacy unseeded single-attempt path). The retry seeds the addHydrogens jitter
    # per unit and re-places on the rare R2 seed-clash FAIL, up to K attempts.
    #
    # RESILIENCE SCOPE (R-18 honest): a per-seed K-EXHAUSTION raises loudly and HALTS
    # the campaign (it propagates through run_leg, which dispatches the cohort's seeds
    # in one call). Full drop-that-seed-from-BOTH-cohorts-and-continue resilience is
    # DEFERRED, not silently dropped: the rep dirs are POSITIONAL (rep0..repN indexed
    # off the --seeds order) and the LEG-axis paired analysis pairs folded[s]-tripep[s]
    # positionally, so cleanly excising one seed mid-cohort while keeping folded/tripep
    # paired would touch the positional rep indexing + the analysis path (not a minimal
    # change). This is acceptable because (a) K>=5 makes exhaustion ~3e-7 (effectively
    # never) and (b) K-exhaustion is a REAL scaffold/mutation geometry problem where
    # HALT + escalate is the CORRECT response, not a stochastic outlier to skip past
    # (review condition C2 / post-run review). Prior seeds' completed ladders are already on disk; recover a
    # partial cohort via --analyze-only on the surviving --seeds subset.
    retry_k = args.appearing_h_retry_k if args.appearing_h_retry_k > 0 else None
    # H18 densification knots — reuse the inplace SSOT parsers VERBATIM (no new
    # parsing logic; anti-fragmentation). GOTCHA: the CLI flag --lambda2-rampdown
    # maps to the engine kwarg lambda2_rampup (the leg-up climbs λ2 UP; "-rampdown"
    # is a symmetry misnomer), EXACTLY as trackb_inplace_rbfe_production.main threads
    # it (its args.lambda2_rampdown -> lambda2_rampup -> run_leg). None => uniform.
    lambda1_rampdown = prod._parse_lambda1_rampdown(args.lambda1_rampdown)
    lambda2_rampup = prod._parse_lambda2_rampdown(args.lambda2_rampdown)
    return prod.run_leg(
        out_root=out_root, endpoint=endpoint, leg="free",
        seeds=seeds, directions=directions,
        n_windows_half=args.n_windows_half, softcore_band=args.softcore_band,
        n_apex_bridge=args.n_apex_bridge, apex_band=args.apex_band,
        n_cycles=args.n_cycles, md_steps_per_cycle=args.md_steps_per_cycle,
        platform_name=args.platform, timestep_fs=args.timestep_fs,
        minimize_iters=args.minimize_iters,
        backward_equil_steps=args.backward_equil_steps,
        genuine_decouple_nm=args.genuine_decouple_nm,
        mtr_ncaa_xml=None, binder_chain=binder_chain,
        archive_existing=(not args.no_archive),
        construction="twocopy", auto_search_displacement=True,
        accept_sep_nm=args.accept_sep_nm,
        lambda1_rampdown=lambda1_rampdown, lambda2_rampup=lambda2_rampup,
        mutation_spec=spec, staged_min=args.staged_min,
        leg_inputs=leg_inputs, appearing_h_retry_k=retry_k)


def main(argv: Optional[List[str]] = None) -> int:
    prod = _load_prod()
    ats = _load_ats()

    p = argparse.ArgumentParser(
        description="Barnase Ile96->Ala folding-thermocycle production orchestrator.",
        epilog=(
            "THROUGHPUT NOTE: the inplace launcher's --pool/--max-concurrent 2-way "
            "concurrency is DEFERRED here (not mirrored). The pool re-invokes the "
            "launcher as CLI --worker subprocesses that resolve --mutation by "
            "REGISTRY NAME, but this orchestrator builds each cohort's Ile->Ala "
            "MutationSpec as an in-process OBJECT with no registry entry (by design), "
            "so it cannot be passed over the worker CLI. Cohorts + seeds therefore "
            "run SERIALLY (run_leg loops seeds). Pin the GPU with "
            "CUDA_VISIBLE_DEVICES=<idx>. 2-way concurrency would require a new "
            "worker mode (out of this additive fix's scope)."))
    p.add_argument("--out-root", required=True,
                   help="cohort dispatch root (folded/ + tripep/ subtrees).")
    p.add_argument("--barnase-final", default=None,
                   help="equilibrated H-complete single-chain barnase final.pdb "
                        "(FOLDED leg scaffold).")
    p.add_argument("--barnase-chain", default="A")
    p.add_argument("--barnase-resnum", type=int, default=None,
                   help="core Ile resnum in the barnase scaffold (default: locate "
                        "the Leu-Ile-Tyr motif = crystal Ile96).")
    p.add_argument("--tripep-final", default=None,
                   help="Ac-Ile-NMe (ACE-ILE-NME) equilibrated final.pdb "
                        "(TRIPEP reference leg scaffold).")
    p.add_argument("--tripep-chain", default="A")
    p.add_argument("--tripep-resnum", type=int, default=2,
                   help="ILE resnum in the tripeptide (ACE=1, ILE=2, NME=3).")
    p.add_argument("--seeds", default="s7,s101,s127,s19,s163,s199",
                   help="comma-separated matched velocity seeds (production n=6).")
    p.add_argument("--directions", default="dplus,dminus",
                   help="ladder directions (both -> UWHAM apex stitch).")
    p.add_argument("--lambda1-rampdown", default=None,
                   help="Two-copy densification (leg-DOWN): comma list of explicit "
                        "leg-down λ1 knots (each in (0,0.5], strictly increasing, "
                        "ending at 0.5; the apex λ1=0 state is placed by the leg-up "
                        "phase). Threaded VERBATIM to run_leg as the engine kwarg "
                        "lambda1_rampdown (parsed with the inplace SSOT "
                        "_parse_lambda1_rampdown). Default None => the canonical "
                        "UNIFORM leg-down (NOT the densified H18 ladder). V3I H18 SSOT "
                        "(v3i_h18_launch.sh): 0.025,0.05,0.1,0.2,0.3,0.4,0.5 = 7 "
                        "leg-down states. Interior-λ reshaping = ΔG-unbiased "
                        "(ranking-only, R-11); soft-core canon untouched (C7).")
    p.add_argument("--lambda2-rampdown", default=None,
                   help="Two-copy densification (leg-UP): comma list of explicit "
                        "leg-up λ2 knots (each in [0,0.5], strictly increasing, MUST "
                        "start at 0.0 and end at 0.5). GOTCHA: the flag name mirrors "
                        "--lambda1-rampdown for symmetry, but the leg-up climbs λ2 UP "
                        "— it is threaded to run_leg as the engine kwarg lambda2_rampup "
                        "(NOT lambda2_rampdown; reversing the mapping silently corrupts "
                        "the ladder), EXACTLY as trackb_inplace_rbfe_production maps it. "
                        "Default None => the canonical UNIFORM leg-up (NOT densified). "
                        "V3I H18 SSOT: 0,0.05,0.1,0.15,0.2,0.3,0.4,0.5 = 8 leg-up "
                        "states. Interior-λ reshaping = ΔG-unbiased (ranking-only, R-11).")
    p.add_argument("--n-windows-half", type=int, default=6,
                   help="Leg-up window count for the UNIFORM-fallback ladder only "
                        "(mirrors the inplace SSOT default 6). The DENSIFIED 15-state "
                        "H18 ladder is produced by --lambda2-rampdown (8 leg-up) + "
                        "--lambda1-rampdown (7 leg-down); when BOTH knot lists are given "
                        "the state counts come from the lists and this value only "
                        "satisfies the n>=2 guard (byte-identical whether 6 or 8). "
                        "WITHOUT the rampdown flags the ladder is UNIFORM (the "
                        "known-FAILING H18 grade), NOT the densified 15-state.")
    p.add_argument("--softcore-band", type=int, default=2)
    p.add_argument("--n-apex-bridge", type=int, default=0)
    p.add_argument("--apex-band", type=float, default=0.5)
    p.add_argument("--n-cycles", type=int, default=prod.DEFAULT_N_CYCLES)
    p.add_argument("--md-steps-per-cycle", type=int,
                   default=prod.DEFAULT_MD_STEPS_PER_CYCLE)
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--timestep-fs", type=float, default=1.0)
    p.add_argument("--minimize-iters", type=int, default=500)
    p.add_argument("--backward-equil-steps", type=int, default=500)
    p.add_argument("--genuine-decouple-nm", type=float, default=1.2)
    p.add_argument("--accept-sep-nm", type=float,
                   default=prod._ATS_ACCEPT_SEP_NM_DEFAULT)
    p.add_argument("--staged-min", action="store_true",
                   help="staged minimization (>=5000 iters floor) — recommended "
                        "for the large buried-core box.")
    p.add_argument("--appearing-h-retry-k", type=int, default=5,
                   help="P3-#116 FIX2: bounded R2-retry budget for the deterministic "
                        "appearing-H (Ala CB methyl) placement. The appearing methyl "
                        "H are placed by an UNSEEDED addHydrogens jitter, so a rare "
                        "(~5%%) build draws a pathological rotamer that trips the R2 "
                        "seed guard (build-ORDER dependent, NOT velocity-seed physics; "
                        "crashed the s127 build). This orchestrator seeds the placement "
                        "DETERMINISTICALLY per unit and re-places (incrementing the "
                        "seed) up to K times on an R2 FAIL. K>=5 -> ~3e-7 exhaustion; "
                        "K exhaustion is FAIL-LOUD (a real-geometry escalation, the "
                        "0.10 nm R2 threshold is never relaxed). Set 0 to DISABLE (fall "
                        "back to the legacy unseeded single-attempt placement). "
                        "PLACEMENT-ONLY: the soft-core C7 canon / λ-schedule / cycles / "
                        "templates / box / frozen FE core are byte-identical.")
    p.add_argument("--no-archive", action="store_true",
                   help="do NOT archive an existing rep dir (R-7 archive is default).")
    p.add_argument("--mintimeid", type=int, default=None,
                   help="UWHAM equilibration discard (analyze phase).")
    p.add_argument("--analyze-only", action="store_true",
                   help="skip dispatch; only run the LEG-axis analysis.")
    p.add_argument("--json-out", default=None)
    args = p.parse_args(argv)

    seeds = [s.strip() for s in args.seeds.split(",") if s.strip()]
    directions = [d.strip() for d in args.directions.split(",") if d.strip()]
    os.makedirs(args.out_root, exist_ok=True)

    dispatch: Dict[str, Any] = {}
    if not args.analyze_only:
        if not args.barnase_final or not args.tripep_final:
            print("ERROR: --barnase-final and --tripep-final are required for "
                  "dispatch (or use --analyze-only). Scaffold MD prep is a "
                  "separate data step.")
            return 2
        barn_resnum = (args.barnase_resnum
                       if args.barnase_resnum is not None
                       else _find_core_ile_resnum(args.barnase_final))
        print("[dispatch] FOLDED cohort: scaffold=%s resnum=%d chain=%s seeds=%s"
              % (args.barnase_final, barn_resnum, args.barnase_chain, seeds))
        dispatch["folded"] = _dispatch_cohort(
            prod, ats, out_root=args.out_root, endpoint="folded",
            scaffold_final=args.barnase_final, resnum=barn_resnum,
            binder_chain=args.barnase_chain, seeds=seeds, directions=directions,
            args=args)
        print("[dispatch] TRIPEP cohort: scaffold=%s resnum=%d chain=%s seeds=%s"
              % (args.tripep_final, args.tripep_resnum, args.tripep_chain, seeds))
        dispatch["tripep"] = _dispatch_cohort(
            prod, ats, out_root=args.out_root, endpoint="tripep",
            scaffold_final=args.tripep_final, resnum=args.tripep_resnum,
            binder_chain=args.tripep_chain, seeds=seeds, directions=directions,
            args=args)

    # LEG-axis analysis (only when both directions ran, so the merged .out exists).
    payload: Dict[str, Any] = {"dispatch_done": (not args.analyze_only)}
    if set(directions) == set(prod.DIRECTION_TAGS):
        payload["analysis"] = analyze_fold_thermocycle(
            args.out_root, leg="free", n_replicates=len(seeds),
            mintimeid=args.mintimeid, output_json=args.json_out)
        a = payload["analysis"]
        print("\nddG_fold(Ile->Ala) = %s kcal/mol  (n=%d, sign=%s)"
              % (a.get("ddG_fold_kcal"), a.get("n_matched_seeds"),
                 a.get("sign_status")))
        print("  reason:", a.get("sign_reason"))
    else:
        print("\n[analysis skipped] both directions (dplus,dminus) needed for the "
              "UWHAM apex stitch; got %s" % (directions,))

    if args.json_out and "analysis" not in payload:
        with open(args.json_out, "w") as fh:
            json.dump({"dispatch": {"folded": bool(dispatch.get("folded")),
                                    "tripep": bool(dispatch.get("tripep"))},
                       "generated": time.strftime("%Y-%m-%d %H:%M:%S")}, fh, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
