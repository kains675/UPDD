#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B — in-place residue-4 RBFE λ-ladder asyncre SMOKE (v0.8 bridge step).

This is the SHORT asyncre validation of the in-place RBFE ladder bridge
(``utils/atm_trackB_inplace_rbfe.py``), NOT full production. It:

  1. SERIALIZES the in-place fused System (genuine HE1<->methyl swap, S0-
     harmonized RBFE XML) to the ``inplace_rbfe_<leg>_sys.xml`` /
     ``inplace_rbfe_<leg>.pdb`` format + verifies the load-only round-trip
     (the per-direction-driver-side deserialize contract).
  2. Builds an INITIAL RBFE λ-ladder (~8-16 windows; the count is a PARAMETER,
     validated empirically here — NOT pre-assumed converged).
  3. Runs a FEW asyncre cycles (default 12) of the ladder on the FREE 5070Ti
     (CUDA), in-process, via the thin ``InplaceRbfeLadder`` adapter that NEVER
     rebuilds the ATMForce (only sets its per-state global parameters).
  4. PARSES the run with the per-direction driver's OWN mixing gate
     (``parse_state_transitions_from_log`` + ``check_atm_mixing``), VERBATIM —
     so the window-count validator is reused, not reimplemented.

What this SMOKE answers (R-18, honest scope):
  (a) does the ladder LOAD + RUN without NaN?  (loads? nan_seen?)
  (b) is there ANY replica-exchange mixing — adjacent crossings / both-ends-
      visited starting to populate — i.e. is the ladder wired + do adjacent
      windows overlap at all?  (walls / both_ends_visited_count / round_trips)

What it does NOT prove (NEXT step, not this one):
  - converged ΔΔG_int / ΔΔG_bind (needs full production + UWHAM + n>=3);
  - the final window count (this measures the INITIAL count's mixing signal;
    add windows where overlap is thin — the empirical escalation);
  - that the per-direction ABFE driver can consume the in-place System (it
    CANNOT without its displacement/LIGAND_ATOMS assumptions — see the module
    docstring; this is a real integration finding, hence the in-process adapter).

Runs openmm-only in the ``qmmm`` env (no atom_openmm needed — the mixing gate is
spec-loaded from the driver without triggering its lazy atom_openmm imports).

DOI references:
  - Gallicchio 2021 J Chem Theory Comput, DOI 10.1021/acs.jctc.1c00753
  - Azimi et al. 2022 J Chem Inf Model 62(2):309, DOI 10.1021/acs.jcim.1c01129
  - Mey et al. 2020 LiveCoMS, DOI 10.33011/livecoms.2.1.18378
"""

import argparse
import importlib.util
import json
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)
_UTILS = os.path.join(_PROJ, "utils")
if _UTILS not in sys.path:
    sys.path.insert(0, _UTILS)

import atm_trackB_inplace_rbfe as rbfe  # noqa: E402


def _load_mixing_gate():
    """Spec-load the per-direction driver's mixing-gate functions WITHOUT
    importing the whole module's lazy atom_openmm chain. The mixing functions
    (parse_state_transitions_from_log / check_atm_mixing) are pure-Python log
    parsers — they do not touch atom_openmm — so importing the module is safe in
    the qmmm env (the atom_openmm imports inside that module are all lazy/local).
    """
    spec = importlib.util.spec_from_file_location(
        "trackb_per_direction_production",
        os.path.join(_PROJ, "scripts", "trackb_per_direction_production.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_ladder_smoke(
    leg="free",
    seed="s7",
    binder_chain="B",
    solvate=True,
    harmonize_common_charges=False,
    swap_mode="genuine",
    genuine_decouple_nm=1.2,
    n_windows_half=8,
    softcore_band=2,
    n_apex_bridge=0,
    apex_band=0.5,
    single_direction=None,
    n_cycles=12,
    md_steps_per_cycle=250,
    warmup_cycles=2,
    min_crossings=1,
    platform_name="CUDA",
    timestep_fs=1.0,
    rng_seed=20260613,
    minimize_iters=500,
    backward_equil_steps=500,
    out_dir=None,
):
    """Serialize -> load-check -> build ladder -> short asyncre -> mixing gate.

    ``single_direction`` (None -> the full symmetric two-direction ladder;
    "forward"/"backward" -> a STANDALONE single-direction ladder, the DECISIVE
    fork test for whether the per-direction-separate estimator round-trips
    INTERNALLY on the single-shared-core box) is forwarded verbatim to
    ``build_rbfe_ladder``. For a single-direction ladder the mixing gate's
    ``round_trips`` measure the INTERNAL state-0 <-> apex traversal of one
    standalone direction (NOT the combined-ladder apex Dir-flip handoff).

    Returns the structured result dict (written to JSON by ``main``).
    """
    t0 = time.time()
    if out_dir is None:
        out_dir = os.path.join(_PROJ, "outputs", "_trackb",
                               "inplace_rbfe_ladder_smoke")
    os.makedirs(out_dir, exist_ok=True)

    result = {
        "outcome": None,
        "leg": leg,
        "seed": seed,
        "swap_mode": swap_mode,
        "solvated": solvate,
        "harmonize_common_charges": harmonize_common_charges,
        "single_direction": single_direction,
        "regime": "ranking_only",
        "prediction_test": True,
        "note": ("In-place RBFE λ-ladder asyncre SMOKE (R-18) — a SHORT mixing "
                 "check of an INITIAL window count, NOT a converged ΔΔG. The "
                 "per-direction ABFE driver cannot consume the in-place System "
                 "(displacement/LIGAND_ATOMS assumptions); this runs the ladder "
                 "via the in-process adapter that reuses the driver's mixing "
                 "gate verbatim."),
    }

    # 1) Serialize + load-only round-trip (the driver-side deserialize contract).
    ser = rbfe.serialize_inplace_rbfe_system(
        leg=leg, out_dir=out_dir, seed=seed, binder_chain=binder_chain,
        solvate=solvate, harmonize_common_charges=harmonize_common_charges,
        swap_mode=swap_mode, genuine_decouple_nm=genuine_decouple_nm,
        constraints=None)
    loaded = rbfe.load_serialized_system(ser["sys_xml_path"], ser["pdb_path"])
    result["serialize"] = {
        "sys_xml_path": ser["sys_xml_path"],
        "pdb_path": ser["pdb_path"],
        "n_atoms": ser["n_atoms"],
        "atmforce_index": ser["atmforce_index"],
        "mtr_ncaa_xml": ser["mtr_ncaa_xml"],
        "alchemical_atoms": ser["alchemical_atoms"],
        "deserialize_atmforce_index": loaded["atmforce_index"],
        "deserialize_n_atoms": loaded["n_atoms"],
        "roundtrip_ok": (loaded["n_atoms"] == ser["n_atoms"]
                         and loaded["atmforce_index"] is not None),
    }

    # 2) Build the INITIAL ladder (optionally with an apex bridge, optionally a
    #    STANDALONE single-direction ladder — the decisive fork test).
    schedule = rbfe.build_rbfe_ladder(
        n_windows_half=n_windows_half, softcore_band=softcore_band,
        n_apex_bridge=n_apex_bridge, apex_band=apex_band,
        single_direction=single_direction)
    result["ladder"] = {
        "schedule_name": schedule["schedule_name"],
        "n_states": schedule["n_states"],
        "n_windows_half": schedule["n_windows_half"],
        "n_apex_bridge": schedule.get("n_apex_bridge"),
        "apex_band": schedule.get("apex_band"),
        "single_direction": schedule.get("single_direction"),
        "softcore_band": schedule["softcore_band"],
        "lambdas_1": schedule["lambdas_1"],
        "lambdas_2": schedule["lambdas_2"],
        "directions": schedule["directions"],
        "intermd": schedule["intermd"],
        "w0": schedule["w0"],
        "alpha": schedule["alpha"],
        "u0": schedule["u0"],
        "umax": schedule["umax"],
        "ubcore": schedule["ubcore"],
        "acore": schedule["acore"],
    }

    # Emit the asyncre cntl for the ladder (the per-direction-driver-side input
    # format) so the schedule is reproducible from disk.
    cntl_path = os.path.join(out_dir, "inplace_rbfe_%s_asyncre.cntl" % (leg,))
    _write_ladder_cntl(cntl_path, schedule, ser, leg, md_steps_per_cycle,
                       timestep_fs)
    result["cntl_path"] = cntl_path

    # 3) Run the SHORT asyncre ladder via the in-process adapter (5070Ti CUDA).
    log_path = os.path.join(out_dir, "inplace_rbfe_%s_driver.log" % (leg,))
    ladder = rbfe.InplaceRbfeLadder(
        loaded["system"], loaded["positions"], schedule,
        platform_name=platform_name, temperature_K=schedule["temperature_K"],
        timestep_fs=timestep_fs, log_path=log_path, seed=rng_seed,
        minimize_iters=minimize_iters,
        backward_equil_steps=backward_equil_steps)
    result["platform"] = ladder.platform_name
    result["minimize_iters"] = minimize_iters
    result["backward_equil_steps"] = backward_equil_steps
    result["backward_endpoint_state"] = ladder._backward_endpoint

    nan_any = False
    nan_states_all = set()
    cycle_log = []
    for _ in range(n_cycles):
        info = ladder.run_cycle(md_steps=md_steps_per_cycle)
        nan_any = nan_any or info["nan_seen"]
        nan_states_all.update(info.get("nan_states", []))
        cycle_log.append({
            "cycle": info["cycle"], "n_accepted": info["n_accepted"],
            "n_pairs": info["n_pairs"], "nan_seen": info["nan_seen"],
            "nan_states": info.get("nan_states", []),
        })
    ladder.close()
    result["run"] = {
        "n_cycles": n_cycles,
        "md_steps_per_cycle": md_steps_per_cycle,
        "timestep_fs": timestep_fs,
        "nan_any": nan_any,
        "nan_states": sorted(nan_states_all),
        "driver_log": log_path,
        "per_cycle": cycle_log,
        "total_accepted": sum(c["n_accepted"] for c in cycle_log),
        "total_pairs_attempted": sum(c["n_pairs"] for c in cycle_log),
    }

    # 4) Mixing gate — REUSED VERBATIM from the per-direction driver.
    pdp = _load_mixing_gate()
    tr = pdp.parse_state_transitions_from_log(
        log_path, warmup_cycles=warmup_cycles)
    try:
        mix = pdp.check_atm_mixing(
            log_path, schedule_K=schedule["n_states"],
            warmup_cycles=warmup_cycles, min_crossings=min_crossings)
        mix_passed = mix.get("passed")
        mix_walls = mix.get("walls")
        both_ends = mix.get("both_ends_visited_count")
        total_round_trips = mix.get("total_round_trips")
        mix_reason = mix.get("reason")
    except Exception as exc:  # noqa: BLE001  (gate is advisory in the smoke)
        mix = {"error": str(exc)}
        mix_passed = None
        mix_walls = None
        both_ends = None
        total_round_trips = None
        mix_reason = "mixing gate raised: %s" % (exc,)

    # Count how many adjacent ladder pairs saw ANY crossing post-warmup (the
    # overlap signal: a pair with 0 crossings is a candidate window to densify).
    adj = tr.get("adjacent_crossings", {})
    n_adjacent_with_crossings = sum(1 for v in adj.values() if v > 0)
    n_adjacent_pairs = schedule["n_states"] - 1
    result["mixing"] = {
        "warmup_cycles": warmup_cycles,
        "min_crossings": min_crossings,
        "n_cycles_total": tr.get("n_cycles_total"),
        "n_samples": tr.get("n_samples"),
        "replicas_seen": tr.get("replicas_seen"),
        "adjacent_crossings": adj,
        "n_adjacent_pairs": n_adjacent_pairs,
        "n_adjacent_with_crossings": n_adjacent_with_crossings,
        "round_trips": tr.get("round_trips"),
        "both_ends_visited_count": both_ends,
        "total_round_trips": total_round_trips,
        "gate_passed": mix_passed,
        "gate_walls": mix_walls,
        "gate_reason": mix_reason,
    }

    # Outcome classification (R-18 — honest; the smoke reports the signal, it
    # does not assert convergence). A per-replica NaN on a SUBSET of states is no
    # longer fatal to the run (the adapter catches + flags it), so a partial
    # outcome distinguishes "specific windows unstable" from "whole run NaN'd".
    nan_states = result["run"]["nan_states"]
    has_mixing = n_adjacent_with_crossings > 0 or (both_ends and both_ends > 0)
    if not result["serialize"]["roundtrip_ok"]:
        result["outcome"] = "serialize_roundtrip_fail"
    elif nan_states and len(nan_states) >= schedule["n_states"]:
        # Every state NaN'd — the box is not runnable at any λ (a build/seed bug).
        result["outcome"] = "ladder_nan"
    elif nan_states:
        # SOME ladder window(s) are unstable while others ran + mixed. The
        # honest finding: report WHICH states NaN'd (a schedule-stability signal
        # — those windows need a softer soft-core / smaller step / a bridge), do
        # not hide it and do not let it abort the whole ladder.
        result["outcome"] = ("ladder_partial_nan_with_mixing" if has_mixing
                             else "ladder_partial_nan_no_mixing")
    elif has_mixing:
        # Wired + at least SOME adjacent overlap — the ladder is runnable and
        # windows overlap somewhere. NEXT: more cycles for round-trips +
        # densify any 0-crossing pair (the empirical window-count escalation).
        result["outcome"] = "ladder_mixing_present"
    else:
        # Ran clean, no NaN, but NO adjacent crossings observed in this short
        # window → either too few cycles or windows do not overlap (need more
        # windows). The honest signal: report the zero crossings, do not hide.
        result["outcome"] = "ladder_no_mixing_observed"

    # DECISIVE single-direction fork-test verdict (R-18 — the smoke reports the
    # number, it does not fabricate a pass). For a STANDALONE single-direction
    # ladder the two mixing-gate ends are the two PHYSICAL λ-endpoints of THAT
    # direction (λ=0 endpoint and λ=0.5 apex), so total_round_trips counts the
    # INTERNAL 0<->apex traversal of one standalone direction. The verdict is
    # ONLY emitted for single_direction ladders (the symmetric ladder's round
    # trips would be the combined-ladder apex handoff, a different question).
    if single_direction is not None:
        rt = total_round_trips if total_round_trips is not None else 0
        all_pairs_cross = (n_adjacent_with_crossings == n_adjacent_pairs
                           and n_adjacent_pairs > 0)
        if nan_states:
            verdict = "indeterminate_nan"
            verdict_msg = (
                "Single-direction (%s) ladder has unstable window(s) %s — the "
                "round-trip count is not yet decisive; soften/bridge those "
                "windows + re-smoke." % (single_direction, nan_states))
        elif rt > 0:
            verdict = "per_direction_viable"
            verdict_msg = (
                "DECISIVE: the standalone %s ladder round-trips INTERNALLY "
                "(total_round_trips=%d, both_ends_visited=%s, %d/%d adjacent "
                "pairs cross). The single-shared-core box CAN sample a "
                "converged direction => the apex direction-flip wall was a "
                "COMBINED-ladder artifact; the per-direction-separate estimator "
                "(dplus+dminus as separate asyncre runs, UWHAM-merged at the "
                "shared apex) is VIABLE — NO two-copy overlay needed."
                % (single_direction, rt, both_ends, n_adjacent_with_crossings,
                   n_adjacent_pairs))
        elif all_pairs_cross:
            verdict = "internal_mixing_no_roundtrip_yet"
            verdict_msg = (
                "The standalone %s ladder mixes across ALL %d adjacent pairs "
                "but recorded 0 end-to-end round-trips in this SHORT smoke "
                "(both_ends_visited=%s). Not yet decisive: every boundary is "
                "permeable, so this is an under-sampling (more cycles) signal, "
                "NOT a structural wall — re-run with more cycles to confirm "
                "round_trips > 0." % (single_direction, n_adjacent_pairs,
                                      both_ends))
        else:
            verdict = "two_copy_overlay_required"
            verdict_msg = (
                "The standalone %s ladder does NOT round-trip internally "
                "(total_round_trips=0, both_ends_visited=%s) AND has a 0-"
                "crossing wall (%d/%d adjacent pairs cross; walls=%s). A single "
                "standalone direction cannot sample a converged path on the "
                "single-shared-core box => the two-copy overlay IS required."
                % (single_direction, both_ends, n_adjacent_with_crossings,
                   n_adjacent_pairs, mix_walls))
        result["single_direction_verdict"] = {
            "single_direction": single_direction,
            "verdict": verdict,
            "total_round_trips": rt,
            "both_ends_visited_count": both_ends,
            "n_adjacent_with_crossings": n_adjacent_with_crossings,
            "n_adjacent_pairs": n_adjacent_pairs,
            "gate_walls": mix_walls,
            "message": verdict_msg,
        }

    result["elapsed_s"] = round(time.time() - t0, 1)
    return result


def _write_ladder_cntl(cntl_path, schedule, ser, leg, md_steps, timestep_fs):
    """Emit the in-place RBFE asyncre cntl (schedule + System/topology pointers).

    This is the per-direction-driver-side input format for the in-place ladder.
    It is NOT consumed by upstream OMMSystemABFE (the in-place System carries its
    own ATMForce; OMMSystemABFE would rebuild one). It records the SAME per-state
    arrays the in-process adapter applies, so the schedule is reproducible from
    disk and a future RBFE driver can read it. ``MODE = INPLACE_RBFE`` marks the
    System as ATMForce-pre-attached so a consumer does NOT re-wrap it.
    """
    def _csv(vals):
        return ", ".join(str(v) for v in vals)

    lines = [
        "# Track B in-place residue-4 RBFE asyncre control file",
        "# Generated by trackb_inplace_rbfe_ladder_smoke.py — DO NOT EDIT BY HAND.",
        "# MODE = INPLACE_RBFE: the System XML carries its OWN ATMForce (genuine",
        "# single-shared-core HE1<->methyl swap). A consumer MUST NOT rebuild the",
        "# ATMForce from LIGAND_ATOMS/DISPLACEMENT (that is the ABFE path) — it",
        "# only sets the per-state global parameters below.",
        "",
        "MODE = 'INPLACE_RBFE'",
        "JOB_TRANSPORT = 'LOCAL_OPENMM'",
        "BASENAME = 'inplace_rbfe_%s'" % (leg,),
        "SYSTEM_XML = '%s'" % (os.path.basename(ser["sys_xml_path"]),),
        "TOPOLOGY_PDB = '%s'" % (os.path.basename(ser["pdb_path"]),),
        "",
        "TEMPERATURES = '%s'" % (schedule["temperature_K"],),
        "LAMBDAS =      '%s'" % (_csv(schedule["lambdas"]),),
        "DIRECTION =    '%s'" % (_csv(schedule["directions"]),),
        "INTERMEDIATE = '%s'" % (_csv(schedule["intermd"]),),
        "LAMBDA1 =      '%s'" % (_csv(schedule["lambdas_1"]),),
        "LAMBDA2 =      '%s'" % (_csv(schedule["lambdas_2"]),),
        "ALPHA =        '%s'" % (_csv(schedule["alpha"]),),
        "U0 =           '%s'" % (_csv(schedule["u0"]),),
        "W0COEFF =      '%s'" % (_csv(schedule["w0"]),),
        "",
        "UMAX = %s" % (schedule["umax"],),
        "UBCORE = %s" % (schedule["ubcore"],),
        "ACORE = %s" % (schedule["acore"],),
        "",
        "# No DISPLACEMENT / LIGAND_ATOMS: this is an in-place single-shared-core",
        "# swap (NE1 attach, dummy-NE1-ref decouple), NOT a whole-binder",
        "# displacement decoupling.",
        "PRODUCTION_STEPS = '%d'" % (md_steps,),
        "TIME_STEP = %s" % (timestep_fs / 1000.0,),  # ps
        "",
    ]
    with open(cntl_path, "w") as fh:
        fh.write("\n".join(lines))


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Track B in-place residue-4 RBFE λ-ladder asyncre smoke.")
    p.add_argument("--leg", choices=["free", "bound"], default="free",
                   help="Thermodynamic leg (default free — the now-free 5070Ti).")
    p.add_argument("--seed", default="s7", help="Endpoint MD seed (default s7).")
    p.add_argument("--binder-chain", default="B")
    p.add_argument("--no-solvate", action="store_true",
                   help="Build the unsolvated box (CPU-cheap mechanical test).")
    p.add_argument("--harmonize-common-charges", action="store_true",
                   help="Force common charges to WT (mechanical path). DEFAULT "
                        "off — relies on the on-disk S0-harmonized RBFE XML for "
                        "MC1 continuity (production-valid).")
    p.add_argument("--swap-mode", choices=["genuine", "var_park", "null"],
                   default="genuine")
    p.add_argument("--genuine-decouple-nm", type=float, default=1.2)
    p.add_argument("--n-windows-half", type=int, default=8,
                   help="Forward-half window count (total states = 2x this). "
                        "The count is a PARAMETER validated empirically by this "
                        "smoke's mixing signal — NOT pre-assumed converged "
                        "(default 8 => 16 states).")
    p.add_argument("--softcore-band", type=int, default=2,
                   help="States at the λ=0.5 apex carrying the soft-core anneal.")
    p.add_argument("--n-apex-bridge", type=int, default=0,
                   help="EXTRA soft-core windows on each half just below the "
                        "λ=0.5 apex (finer W0 near saturation; default 0 = no "
                        "bridge). The right count is found empirically by this "
                        "smoke (R-18). NOTE: for the single-shared-core in-place "
                        "box the apex direction-flip wall is NOT soft-core-"
                        "bridgeable (the ATM base u0/u1 is un-softened) — the "
                        "bridge is correct infra but cannot close that wall.")
    p.add_argument("--apex-band", type=float, default=0.5,
                   help="Fraction of the soft-core band (from λ=0.5 inward) the "
                        "apex-bridge windows occupy (default 0.5).")
    p.add_argument("--single-direction", choices=["forward", "backward"],
                   default=None,
                   help="Build + run ONLY one STANDALONE direction half (the "
                        "DECISIVE fork test): 'forward' = the dplus ladder (λ1 "
                        "0->0.5, DIRECTION=+1, n-windows-half states), 'backward' "
                        "= the dminus ladder (reverse, DIRECTION=-1). NO apex "
                        "direction-flip boundary, so the mixing gate's "
                        "round_trips measure the INTERNAL state-0 <-> apex "
                        "traversal of one standalone direction — does the per-"
                        "direction-separate estimator round-trip on the single-"
                        "shared-core box? (default None = full symmetric ladder).")
    p.add_argument("--n-cycles", type=int, default=12,
                   help="Asyncre cycles to run (SHORT smoke; default 12).")
    p.add_argument("--md-steps-per-cycle", type=int, default=250,
                   help="MD steps per replica per cycle (default 250).")
    p.add_argument("--warmup-cycles", type=int, default=2,
                   help="Mixing-gate warmup cycles to skip (default 2).")
    p.add_argument("--min-crossings", type=int, default=1)
    p.add_argument("--platform", default="CUDA",
                   help="OpenMM platform (default CUDA — the 5070Ti).")
    p.add_argument("--timestep-fs", type=float, default=1.0)
    p.add_argument("--rng-seed", type=int, default=20260613)
    p.add_argument("--minimize-iters", type=int, default=500,
                   help="Per-replica minimization iterations at the assigned "
                        "state before integration (default 500; the validated "
                        "Tier-2 path requires this on the fresh PME box).")
    p.add_argument("--backward-equil-steps", type=int, default=500,
                   help="Steps to equilibrate each backward (Direction<0) "
                        "replica at the backward endpoint (u1-decoupled basin) "
                        "BEFORE its assigned state (default 500). The backward "
                        "soft-core anneal-edge NaNs if seeded from the u0 "
                        "geometry; 0 disables (reproduces the partial-NaN "
                        "finding).")
    p.add_argument("--json-out", default=None)
    args = p.parse_args(argv)

    result = run_ladder_smoke(
        leg=args.leg, seed=args.seed, binder_chain=args.binder_chain,
        solvate=not args.no_solvate,
        harmonize_common_charges=args.harmonize_common_charges,
        swap_mode=args.swap_mode, genuine_decouple_nm=args.genuine_decouple_nm,
        n_windows_half=args.n_windows_half, softcore_band=args.softcore_band,
        n_apex_bridge=args.n_apex_bridge, apex_band=args.apex_band,
        single_direction=args.single_direction,
        n_cycles=args.n_cycles, md_steps_per_cycle=args.md_steps_per_cycle,
        warmup_cycles=args.warmup_cycles, min_crossings=args.min_crossings,
        platform_name=args.platform, timestep_fs=args.timestep_fs,
        rng_seed=args.rng_seed, minimize_iters=args.minimize_iters,
        backward_equil_steps=args.backward_equil_steps)

    json_out = args.json_out
    if json_out is None:
        out_dir = os.path.join(_PROJ, "outputs", "_trackb",
                               "inplace_rbfe_ladder_smoke")
        os.makedirs(out_dir, exist_ok=True)
        n_states_total = result.get("ladder", {}).get(
            "n_states", 2 * args.n_windows_half)
        bridge_tag = ("_ab%d" % args.n_apex_bridge) if args.n_apex_bridge else ""
        dir_tag = ("_%s" % args.single_direction) if args.single_direction else ""
        json_out = os.path.join(
            out_dir, "ladder_smoke_%s%s_%dw%s_%s.json"
            % (args.leg, dir_tag, n_states_total, bridge_tag, args.seed))
    print(json.dumps(result, indent=2, default=str))
    with open(json_out, "w") as fh:
        json.dump(result, fh, indent=2, default=str)

    # Exit codes (R-18 — each is a distinct empirical signal, not pass/fail):
    #   0 = mixing present (ladder wired + overlapping; NEXT: convergence)
    #   3 = partial NaN but the rest of the ladder mixed (specific unstable
    #       window(s); NEXT: soften those windows + re-smoke)
    #   4 = ran clean but NO mixing observed (densify / more cycles)
    #   5 = whole-ladder NaN (build/seed bug)
    #   6 = serialize round-trip fail
    outcome = result.get("outcome")
    if outcome == "ladder_mixing_present":
        return 0
    if outcome in ("ladder_partial_nan_with_mixing",
                   "ladder_partial_nan_no_mixing"):
        return 3
    if outcome == "ladder_no_mixing_observed":
        return 4
    if outcome == "ladder_nan":
        return 5
    return 6


if __name__ == "__main__":
    sys.exit(main())
