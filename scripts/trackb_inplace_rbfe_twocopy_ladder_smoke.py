#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B — CANONICAL ATS TWO-COPY RBFE λ-ladder asyncre SMOKE (v0.8 wiring).

This is the SHORT asyncre validation that the CANONICAL ATS TWO-COPY box
(canonical two-copy ATS construction;
wiring-validated) now RUNS THROUGH THE LADDER — i.e. that
the two-copy build is WIRED into serialize + load + ``InplaceRbfeLadder``, NOT
just buildable. It is the wiring counterpart of
``trackb_inplace_rbfe_ladder_smoke.py`` (single-shared-core).

It:
  1. SERIALIZES the two-copy box (``construction='twocopy'``) — copy-1 (MTR) at
     the site + copy-2 (WT) displaced ~40 A into bulk, common-coord swap with
     DISTINCT attach atoms, NO inter-copy exclusions — to the
     ``inplace_rbfe_<leg>_sys.xml`` / ``.pdb`` format + verifies the load-only
     round-trip (the ATMForce survives, larger box than the single-endpoint one).
  2. Builds the CANONICAL ATS standard standalone λ-schedule
     (``build_ats_standard_ladder``; spec C: λ1=0/λ2 climbs then λ2=0.5/λ1
     climbs; n_windows_half=6 => 11 λ/leg; Uh=110, umax/ubcore/acore canon).
  3. Runs a FEW asyncre cycles of THAT ladder via the in-process
     ``InplaceRbfeLadder`` adapter (the same adapter the production launcher
     drives), writing per-walker ``.out`` rows.
  4. Re-affirms the Tier-1 escape on the running box: NaN == 0, the .out files
     were written, and the per-walker raw |u1-u0| is the physically-plausible
     few-kcal/mol two-copy transfer (NOT the single-shared-core saturated ~150
     plateau, NOT a 56,000 clash).

What this SMOKE answers (R-18, honest scope):
  (a) does the two-copy box LOAD + RUN through the ladder without NaN? (loads?
      nan_seen?)
  (b) are the per-walker .out files produced? (the UWHAM-consumable layout)
  (c) is the running-box |u1-u0| still finite + NON-saturated (few kcal/mol)? —
      the collapse-signature escape re-checked under integration, not just one
      static frame.

What it does NOT prove (the PILOT, NOT this — full production run):
  - converged ΔΔG_int / ΔΔG_bind (needs full production + UWHAM + n>=3);
  - frac<UBCORE> > 0 + adjacent MBAR-O >= 0.10 + dgbind1 != 0.5 + UWHAM
    convergence + SIMULTANEOUS endpoint-equivalence (the definitive verdict);
  - the final window count (this measures the INITIAL count's runnability).

Ranking-only (R-11); v0.8 PREDICTION test (R-18); NOT a Magotti/absolute
comparison. Runs openmm-only in the ``qmmm`` env (no atom_openmm needed — the
mixing gate is spec-loaded from the driver without its lazy atom_openmm imports).
The 5070Ti (CUDA) is the target; falls back to Reference if CUDA is unavailable.

DOI references:
  - Gallicchio 2025 J Chem Inf Model, ATS DOI 10.1021/acs.jcim.5c00207
    (preprint arXiv:2412.19971; copy-2 bulk d-displacement + common coord swap)
  - Gallicchio 2021 J Chem Theory Comput, DOI 10.1021/acs.jctc.1c00753
  - Azimi et al. 2022 J Chem Inf Model 62(2):309, DOI 10.1021/acs.jcim.1c01129
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


# The two-copy ATS transfer is a real coordinate transfer in the ABFE-cliff band,
# so the per-walker raw |u1-u0| must be FINITE and physically plausible — the
# upper bound rejects a parked-clash / overlay-collapse artifact, and the
# saturated band (140..201) flags the single-shared-core collapse signature.
U1_MINUS_U0_PLAUSIBLE_MAX_KCAL = 1.0e3
SATURATED_FLOOR_KCAL = 140.0       # single-shared-core collapse plateau floor
SATURATED_CEIL_KCAL = 201.0        # ~Umax (200) + slack


def _load_mixing_gate():
    """Spec-load the per-direction driver's mixing-gate functions WITHOUT
    importing the whole module's lazy atom_openmm chain (the mixing functions are
    pure-Python log parsers; the atom_openmm imports in that module are lazy)."""
    spec = importlib.util.spec_from_file_location(
        "trackb_per_direction_production",
        os.path.join(_PROJ, "scripts", "trackb_per_direction_production.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _read_out_pert_range(out_dir, n_states, basename):
    """Scan the per-walker .out files; return (n_rows, min|pert|, max|pert|).

    The .out row is ``stateid T dir l1 l2 alpha u0 w0 potE pertE 0`` (col 9 =
    pertE = the soft-core-CAPPED usc the UWHAM estimator consumes). For the
    running-box plausibility re-check we want the soft-core-capped magnitude;
    |usc| <= Umax(200), so a max |usc| pinned at ~150 would be the collapse
    plateau and a small |usc| (few kcal/mol) is the genuine two-copy transfer.
    """
    n_rows = 0
    pmin = None
    pmax = None
    for r in range(n_states):
        path = os.path.join(out_dir, "r%d" % (r,), basename + ".out")
        if not os.path.isfile(path):
            continue
        with open(path) as fh:
            for line in fh:
                cols = line.split()
                if len(cols) < 10:
                    continue
                try:
                    pert = abs(float(cols[9]))
                except ValueError:
                    continue
                n_rows += 1
                pmin = pert if pmin is None else min(pmin, pert)
                pmax = pert if pmax is None else max(pmax, pert)
    return n_rows, pmin, pmax


def run_twocopy_ladder_smoke(
    leg="free",
    seed="s7",
    binder_chain="B",
    solvate=True,
    harmonize_common_charges=False,
    displacement_nm=None,
    n_windows_half=6,
    single_direction="forward",
    n_cycles=6,
    md_steps_per_cycle=100,
    warmup_cycles=1,
    min_crossings=1,
    platform_name="CUDA",
    timestep_fs=1.0,
    rng_seed=20260614,
    minimize_iters=500,
    out_dir=None,
):
    """Serialize(twocopy) -> load -> ATS-standard ladder -> short asyncre -> checks.

    Returns the structured result dict (written to JSON by ``main``). Raises
    AssertionError on a hard plausibility violation (running-box |u1-u0|
    saturated / parked-clash) so an executor sees a FAIL; the MC1 charge-
    discontinuity build outcome is surfaced structured (a real finding) instead.
    """
    t0 = time.time()
    if out_dir is None:
        out_dir = os.path.join(_PROJ, "outputs", "_trackb",
                               "inplace_rbfe_twocopy_ladder_smoke")
    os.makedirs(out_dir, exist_ok=True)

    result = {
        "outcome": None,
        "leg": leg,
        "seed": seed,
        "construction": "twocopy",
        "single_direction": single_direction,
        "solvated": solvate,
        "harmonize_common_charges": harmonize_common_charges,
        "regime": "ranking_only",
        "prediction_test": True,
        "note": ("CANONICAL ATS two-copy λ-ladder asyncre WIRING smoke (R-18) — a "
                 "SHORT runnability check that the two-copy box runs through "
                 "serialize -> load -> InplaceRbfeLadder, NOT a converged ΔΔG. "
                 "The full overlap/UWHAM pilot is the full production run."),
    }

    # 1) Serialize the TWO-COPY box + load-only round-trip. A two-copy MC1 charge
    #    discontinuity is a real finding (surfaced structured, not a crash).
    serialize_kwargs = dict(
        leg=leg, out_dir=out_dir, seed=seed, binder_chain=binder_chain,
        solvate=solvate, harmonize_common_charges=harmonize_common_charges,
        constraints=None, construction="twocopy")
    if displacement_nm is not None:
        serialize_kwargs["displacement_nm"] = displacement_nm
    try:
        ser = rbfe.serialize_inplace_rbfe_system(**serialize_kwargs)
    except RuntimeError as exc:
        if "mc1" in str(exc).lower() or "charge discontinuity" in str(exc).lower():
            result["outcome"] = "mc1_charge_discontinuity"
            result["mc1_error"] = str(exc)
            result["note"] += (" [MC1: the two-copy common core is charge-"
                               "discontinuous; supply the harmonized RBFE XML.]")
            result["elapsed_s"] = round(time.time() - t0, 1)
            return result
        raise

    loaded = rbfe.load_serialized_system(ser["sys_xml_path"], ser["pdb_path"])
    result["serialize"] = {
        "sys_xml_path": ser["sys_xml_path"],
        "pdb_path": ser["pdb_path"],
        "n_atoms": ser["n_atoms"],
        "n_copy1": ser.get("n_copy1"),
        "atmforce_index": ser["atmforce_index"],
        "construction": ser.get("construction"),
        "displacement_vector_nm": ser.get("displacement_vector_nm"),
        "genuine_decouple_dir": ser.get("genuine_decouple_dir"),
        "mtr_ncaa_xml": ser.get("mtr_ncaa_xml"),
        "alchemical_atoms": ser.get("alchemical_atoms"),
        "deserialize_atmforce_index": loaded["atmforce_index"],
        "deserialize_n_atoms": loaded["n_atoms"],
        "roundtrip_ok": (loaded["n_atoms"] == ser["n_atoms"]
                         and loaded["atmforce_index"] is not None),
    }

    # 2) Build the CANONICAL ATS standard standalone schedule (spec C).
    schedule = rbfe.build_ats_standard_ladder(
        n_windows_half=n_windows_half, single_direction=single_direction)
    result["ladder"] = {
        "schedule_name": schedule["schedule_name"],
        "schedule_kind": schedule.get("schedule_kind"),
        "n_states": schedule["n_states"],
        "lambdas_1": schedule["lambdas_1"],
        "lambdas_2": schedule["lambdas_2"],
        "directions": schedule["directions"],
        "intermd": schedule["intermd"],
        "u0": schedule["u0"],
        "alpha": schedule["alpha"],
        "w0": schedule["w0"],
        "umax": schedule["umax"],
        "ubcore": schedule["ubcore"],
        "acore": schedule["acore"],
    }

    # 3) Run the SHORT asyncre ladder via the in-process adapter, writing the
    #    per-walker .out tree (the production layout) so the .out plausibility +
    #    generation can be re-checked from disk.
    log_path = os.path.join(out_dir, "twocopy_%s_%s_driver.log"
                            % (leg, single_direction))
    base = "trackb_" + ("dplus" if single_direction == "forward" else "dminus")
    ladder = rbfe.InplaceRbfeLadder(
        loaded["system"], loaded["positions"], schedule,
        platform_name=platform_name, temperature_K=schedule["temperature_K"],
        timestep_fs=timestep_fs, log_path=log_path, seed=rng_seed,
        minimize_iters=minimize_iters, backward_equil_steps=500,
        out_dir=out_dir, out_basename=base)
    result["platform"] = ladder.platform_name

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

    n_rows, pmin, pmax = _read_out_pert_range(out_dir, schedule["n_states"], base)
    saturated = bool(
        pmax is not None and SATURATED_FLOOR_KCAL <= pmax <= SATURATED_CEIL_KCAL)
    result["run"] = {
        "n_cycles": n_cycles,
        "md_steps_per_cycle": md_steps_per_cycle,
        "timestep_fs": timestep_fs,
        "nan_any": nan_any,
        "nan_states": sorted(nan_states_all),
        "driver_log": log_path,
        "out_basename": base,
        "out_rows_written": n_rows,
        "pert_abs_min_kcal": pmin,
        "pert_abs_max_kcal": pmax,
        "pert_saturated_plateau": saturated,
        "per_cycle": cycle_log,
        "total_accepted": sum(c["n_accepted"] for c in cycle_log),
        "total_pairs_attempted": sum(c["n_pairs"] for c in cycle_log),
    }

    # 4) Advisory mixing gate (the SHORT smoke does not assert round-trips — that
    #    is the pilot; this just confirms the log parses + the adapter wired).
    pdp = _load_mixing_gate()
    tr = pdp.parse_state_transitions_from_log(log_path, warmup_cycles=warmup_cycles)
    adj = tr.get("adjacent_crossings", {})
    result["mixing"] = {
        "n_cycles_total": tr.get("n_cycles_total"),
        "n_samples": tr.get("n_samples"),
        "n_adjacent_pairs": schedule["n_states"] - 1,
        "n_adjacent_with_crossings": sum(1 for v in adj.values() if v > 0),
        "adjacent_crossings": adj,
    }

    # --- HARD plausibility ASSERTS on the RUNNING box (R-18 escape re-check) ---
    if not result["serialize"]["roundtrip_ok"]:
        result["outcome"] = "serialize_roundtrip_fail"
        result["elapsed_s"] = round(time.time() - t0, 1)
        return result
    if nan_states_all and len(nan_states_all) >= schedule["n_states"]:
        result["outcome"] = "ladder_nan"
        result["elapsed_s"] = round(time.time() - t0, 1)
        return result

    assert n_rows > 0, (
        "WIRING FAIL: no per-walker .out rows were written — the ladder did not "
        "run or the .out writer is not wired.")
    if pmax is not None:
        assert pmax <= U1_MINUS_U0_PLAUSIBLE_MAX_KCAL, (
            "ESCAPE FAIL: running-box max |usc|=%.1f kcal/mol exceeds the "
            "plausible bound %.1f — a parked-clash / overlay-collapse artifact, "
            "NOT a genuine two-copy transfer." % (pmax, U1_MINUS_U0_PLAUSIBLE_MAX_KCAL))
        assert not saturated, (
            "ESCAPE FAIL: running-box max |usc|=%.1f kcal/mol is pinned in the "
            "soft-core saturated band [%.0f, %.0f] — the single-shared-core "
            "COLLAPSE signature re-appeared in the two-copy box."
            % (pmax, SATURATED_FLOOR_KCAL, SATURATED_CEIL_KCAL))

    if nan_states_all:
        result["outcome"] = "ladder_partial_nan"
    else:
        result["outcome"] = "twocopy_ladder_runs"
    result["elapsed_s"] = round(time.time() - t0, 1)
    return result


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Track B canonical ATS two-copy RBFE λ-ladder asyncre WIRING "
                    "smoke (v0.8 PREDICTION test; ranking-only).")
    p.add_argument("--leg", choices=["free", "bound"], default="free")
    p.add_argument("--seed", default="s7")
    p.add_argument("--binder-chain", default="B")
    p.add_argument("--no-solvate", action="store_true",
                   help="Build unsolvated (NoCutoff) — cheap CPU path.")
    p.add_argument("--harmonize-common-charges", action="store_true",
                   help="DIAGNOSTIC: force copy-2 common charges to copy-1 so MC1 "
                        "passes (the harmonized RBFE XML already passes on disk).")
    p.add_argument("--displacement-nm", type=float, default=None,
                   help="Magnitude of the copy-2 bulk displacement d (nm; default "
                        "ATS_TWOCOPY_DISPLACEMENT_NM ~4.0).")
    p.add_argument("--n-windows-half", type=int, default=6,
                   help="ATS leg-up phase length (n=6 => 11 λ/leg, spec C).")
    p.add_argument("--single-direction", choices=["forward", "backward"],
                   default="forward")
    p.add_argument("--n-cycles", type=int, default=6)
    p.add_argument("--md-steps-per-cycle", type=int, default=100)
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--timestep-fs", type=float, default=1.0)
    p.add_argument("--minimize-iters", type=int, default=500)
    p.add_argument("--rng-seed", type=int, default=20260614)
    p.add_argument("--json-out", default=None)
    args = p.parse_args(argv)

    json_out = args.json_out
    if json_out is None:
        out_dir = os.path.join(_PROJ, "outputs", "_trackb",
                               "inplace_rbfe_twocopy_ladder_smoke")
        os.makedirs(out_dir, exist_ok=True)
        json_out = os.path.join(
            out_dir, "twocopy_ladder_%s_%s_%s.json"
            % (args.leg, args.single_direction, args.seed))

    try:
        result = run_twocopy_ladder_smoke(
            leg=args.leg, seed=args.seed, binder_chain=args.binder_chain,
            solvate=not args.no_solvate,
            harmonize_common_charges=args.harmonize_common_charges,
            displacement_nm=args.displacement_nm,
            n_windows_half=args.n_windows_half,
            single_direction=args.single_direction,
            n_cycles=args.n_cycles, md_steps_per_cycle=args.md_steps_per_cycle,
            platform_name=args.platform, timestep_fs=args.timestep_fs,
            rng_seed=args.rng_seed, minimize_iters=args.minimize_iters)
    except AssertionError as exc:
        result = {"outcome": "twocopy_ladder_escape_fail",
                  "assertion_error": str(exc), "regime": "ranking_only",
                  "prediction_test": True, "construction": "twocopy",
                  "leg": args.leg}
        print(json.dumps(result, indent=2, default=str))
        with open(json_out, "w") as fh:
            json.dump(result, fh, indent=2, default=str)
        return 5

    print(json.dumps(result, indent=2, default=str))
    with open(json_out, "w") as fh:
        json.dump(result, fh, indent=2, default=str)

    # Exit codes: 0 = runs; 5 = escape fail; 6 = mc1; 7 = roundtrip/nan.
    if result.get("outcome") == "twocopy_ladder_runs":
        return 0
    if result.get("outcome") == "mc1_charge_discontinuity":
        return 6
    return 7


if __name__ == "__main__":
    sys.exit(main())
