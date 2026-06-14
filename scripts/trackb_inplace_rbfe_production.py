#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B in-place residue-4 RBFE — per-direction PRODUCTION launcher.

Runs the in-place HE1<->methyl RBFE estimator as a per-direction-SEPARATE
ladder set and yields, after UWHAM, the ranking-only

    ddG_bind  =  ddG_int_bound  -  ddG_int_free      (Gallicchio ATM/ABFE cycle)

The estimator GEOMETRY is identical to how the validated ABFE per-direction
driver runs: for each {leg in free,bound} x {direction in dplus(forward),
dminus(backward)} x {replicate 1..n, matched seeds} this launcher runs a
STANDALONE single-direction RBFE ladder (build_rbfe_ladder(single_direction=...)),
the two directions NEVER exchange, and they are merged POST-HOC by UWHAM at the
shared lambda=0.5 apex lambda-state (merge_per_direction_outputs). The decisive
fork test proved both directions round-trip INTERNALLY on the single-shared-core
box (free fwd round_trips=6, bwd=5; bound fwd round_trips=11 at 50 cycles), so
the per-direction-separate estimator is viable WITHOUT a two-copy overlay.

WHY a bespoke launcher (R-18, honest): the per-direction ABFE driver
(scripts/trackb_per_direction_production.py) is hard-wired to upstream
abfe_production, which deserializes a PLAIN System XML and REBUILDS its own
ATMForce from cntl LIGAND_ATOMS/DISPLACEMENT (whole-binder decoupling). The
in-place fused box is SINGLE-shared-core (the genuine HE1 decouple lives only in
the attached ATMForce), so that driver cannot consume it. This launcher instead
drives the validated in-process adapter (InplaceRbfeLadder) that NEVER rebuilds
the ATMForce — it only sets the per-state global parameters — and writes the
SAME per-walker .out layout the asyncre engine emits, so the EXISTING merge +
UWHAM postprocess path consumes it unchanged. The ABFE displacement path / S0
charges / shared XML / densified arrays are UNTOUCHED.

This is PURELY ADDITIVE (off-campaign): it imports the validated bridge module
+ the UWHAM postprocess + the per-direction merge/mixing-gate VERBATIM. It is a
new launcher script + nothing else.

Ranking-only (R-11). Magotti SSOT comparison is SIGN ONLY, never magnitude. The
SIGN of any ddG is NOT claimed until z_SE >= 3 (C7, Welch-Satterthwaite
quadrature, NOT paired) AND n >= 3 matched seeds (C5). The free leg alone is the
EARLY SIGNAL; ddG_bind requires the bound leg (the multi-day follow-up).

C1-C12 gates (the prerequisites this launcher enforces / records):
  C1  free leg runs in its OWN in-place RBFE box (its own ladder, NOT cross-box).
  C3  apex co-sampling — dplus last + dminus first lambda=0.5 both occupancy>0.
  C4  round-trip >= 1 per replicate (per-direction mixing gate, reused verbatim).
  C5  n >= 3 matched seeds.
  C6  per-adjacent MBAR overlap O >= 0.1 (flag < 0.03) with empirical window
      escalation (cap 2 — densify thin pairs via --n-windows-half escalation).
  C7  z_SE >= 3 (Welch-Satterthwaite quadrature) for any SIGN claim.
  C8  BOUND decouple receptor-away (SIGN-critical) — genuine_decouple_dir is the
      finite outward (low-density) unit vector for the bound leg, fail-loud.
  C9  NaN fail-fast + per-cycle data-integrity logging (round-trips / apex
      occupancy, NOT just occupancy).
  C11 pre-register the 4 outcomes (written to a pre-registration JSON).
  C12 R-7 — completed run dirs are archived, never deleted.

DOI references:
  - Gallicchio 2021 J Chem Theory Comput, DOI 10.1021/acs.jctc.1c00753
  - Azimi et al. 2022 J Chem Inf Model 62(2):309, DOI 10.1021/acs.jcim.1c01129
  - Shirts & Chodera 2008 J Chem Phys 129:124105, DOI 10.1063/1.2978177 (MBAR)
  - Klimovich/Shirts/Mobley 2015, DOI 10.1007/s10822-015-9840-9 (overlap diag)
  - Welch 1947 Biometrika 34:28 (unequal-variance quadrature)
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)
_UTILS = os.path.join(_PROJ, "utils")
_SCRIPTS = os.path.join(_PROJ, "scripts")
for _p in (_UTILS, _SCRIPTS):
    if _p not in sys.path:
        sys.path.insert(0, _p)


# Endpoints + per-direction tags (fixed conventions, mirror the ABFE launcher).
ENDPOINTS = ("cp4", "wt")
DIRECTION_TAGS = ("dplus", "dminus")               # forward, backward
DIRECTION_OF_TAG = {"dplus": "forward", "dminus": "backward"}
JOBNAME = "trackb"

# C5 — minimum matched-seed replicates before a SIGN claim is even attempted.
MIN_REPLICATES_FOR_SIGN = 3
# C7 — z_SE threshold (|mean| / quadrature SE) below which the sign is "undetermined".
Z_SE_SIGN_THRESHOLD = 3.0
# C6 — empirical window-count escalation cap (densify thin pairs at most twice).
WINDOW_ESCALATION_CAP = 2
# C6 — adjacent MBAR-overlap targets: O >= WELL is good, O < FLAG is a collapse.
OVERLAP_WELL = 0.10
OVERLAP_FLAG = 0.03

# Default per-direction matched seeds (n=3). The QM-snapshot seed cohort; each
# also gets a distinct integer velocity seed (replicate index + 1) so the
# Langevin RNG differs per replicate (genuine sigma_btwn). NOT a single seed
# repeated — matched ACROSS the cp4/wt endpoints (paired) but distinct across
# replicates within an endpoint. These three QM-snapshot seeds are present for
# BOTH endpoints (2QKI_Cp4_hybrid_calib_<seed> + 2QKI_WT_calib_<seed> final
# PDBs) so the matched-seed pairing is real (seed-availability is gated
# fail-loud before launch — a missing endpoint PDB is a real blocker, R-18).
DEFAULT_SEEDS = ("s7", "s101", "s127")


def _endpoint_pdb_paths(endpoint: str, seed: str) -> Tuple[str, str]:
    """Return (dir, final_pdb) for an endpoint+seed (mirrors the build's
    resolution: Cp4 -> 2QKI_Cp4_hybrid_calib_<seed>; WT -> 2QKI_WT_calib_<seed>)."""
    if endpoint == "cp4":
        d = os.path.join(_PROJ, "outputs", "2QKI_Cp4_hybrid_calib_%s" % (seed,))
        return d, os.path.join(d, "mdresult", "2QKI_Cp4_final.pdb")
    d = os.path.join(_PROJ, "outputs", "2QKI_WT_calib_%s" % (seed,))
    return d, os.path.join(d, "mdresult", "2QKI_WT_final.pdb")


def gate_seed_availability(endpoints: List[str], seeds: List[str]
                           ) -> Dict[str, Any]:
    """Fail-loud gate: every (endpoint, seed) must have its endpoint final PDB.

    A missing endpoint PDB is a REAL blocker — the matched-seed pairing across
    cp4/wt is the basis of the paired ddG, so a seed present for one endpoint but
    not the other is NOT a valid matched seed (R-18: a real gap, not papered
    over with a single repeated seed)."""
    missing: List[str] = []
    present: List[str] = []
    for seed in seeds:
        for ep in endpoints:
            _d, pdb = _endpoint_pdb_paths(ep, seed)
            if os.path.isfile(pdb):
                present.append("%s/%s" % (ep, seed))
            else:
                missing.append("%s/%s -> %s" % (ep, seed, pdb))
    return {"passed": not missing, "present": present, "missing": missing}

# Default production cycle length. The fork test needed ~50 cycles for the bound
# leg to round-trip; production runs many round-trips per replica. 400 cycles at
# 250 MD steps/cycle = 100k MD steps/replica (~the per-direction structprep
# warmup order), comfortably multiple round-trips.
DEFAULT_N_CYCLES = 400
DEFAULT_MD_STEPS_PER_CYCLE = 250


# ---------------------------------------------------------------------------
# Module loaders (importlib spec-load to avoid sys.modules pollution + so the
# atom_openmm-free qmmm env can import the openmm-only bridge + the mixing gate;
# the UWHAM postprocess lazily imports atom_openmm only inside the estimator).
# ---------------------------------------------------------------------------
def _load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError("cannot locate %s at %s" % (name, path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_rbfe():
    return _load_module(
        "atm_trackB_inplace_rbfe",
        os.path.join(_UTILS, "atm_trackB_inplace_rbfe.py"))


def _load_driver():
    """The per-direction driver — REUSED for merge_per_direction_outputs +
    the mixing gate (parse_state_transitions_from_log / check_atm_mixing). This
    module is NOT modified by this launcher."""
    return _load_module(
        "trackb_per_direction_production",
        os.path.join(_SCRIPTS, "trackb_per_direction_production.py"))


def _load_uwham():
    return _load_module(
        "trackb_uwham_postprocess",
        os.path.join(_SCRIPTS, "trackb_uwham_postprocess.py"))


# ---------------------------------------------------------------------------
# Output layout (mirrors the ABFE per-direction tree so merge + UWHAM consume it
# verbatim):
#
#   <out_root>/<endpoint>/<leg>/rep{j}/
#       trackb_asyncre.cntl                 (COMBINED symmetric schedule cntl —
#                                            the SSOT both merge + UWHAM derive
#                                            the state count + per-state arrays
#                                            from)
#       dplus/  trackb_dplus_asyncre.cntl   (single-direction forward cntl)
#               r0../trackb_dplus.out
#       dminus/ trackb_dminus_asyncre.cntl  (single-direction backward cntl)
#               r0../trackb_dminus.out
#       r0..r{K-1}/trackb.out               (written by merge_per_direction_outputs)
#       inplace_rbfe_<leg>_sys.xml / .pdb   (the serialized in-place System)
#       run_manifest.json                   (per-direction mixing + data-integrity)
# ---------------------------------------------------------------------------
def _rep_dir(out_root: str, endpoint: str, leg: str, replicate_index: int) -> str:
    return os.path.join(out_root, endpoint, leg, "rep%d" % (replicate_index,))


def _csv(vals) -> str:
    return ", ".join(str(v) for v in vals)


def _build_single_direction_schedule(
    rbfe, *, construction: str, direction: str,
    n_windows_half: int, softcore_band: int,
    n_apex_bridge: int, apex_band: float,
    lambda1_rampdown: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Build ONE standalone direction's schedule for the chosen construction.

    ``single_core`` (default) -> ``build_rbfe_ladder(single_direction=...)`` (the
    validated single-shared-core λ1==λ2 linear + apex anneal band ladder).
    ``twocopy`` -> ``build_ats_standard_ladder(single_direction=...)`` (the
    canonical ATS two-phase per-leg schedule, spec C: λ1=0/λ2 climbs then
    λ2=0.5/λ1 climbs, single DIRECTION, Uh=110). The two-copy box's swap is a real
    transfer, so it uses the ATS canon schedule, NOT the single-shared-core ladder.

    ``lambda1_rampdown`` (two-copy ONLY): explicit leg-down λ1 knots that densify
    the leg-switch handoff (e.g. ``[0.05,0.1,0.2,0.3,0.4,0.5]`` = a λ1=0.05 bridge
    window => 12 λ/leg). ``None`` keeps the canonical uniform leg-down (11 λ/leg).
    Ignored for single_core (it has no leg-switch boundary).
    """
    if construction == "twocopy":
        return rbfe.build_ats_standard_ladder(
            n_windows_half=n_windows_half, single_direction=direction,
            lambda1_rampdown=lambda1_rampdown)
    return rbfe.build_rbfe_ladder(
        n_windows_half=n_windows_half, softcore_band=softcore_band,
        n_apex_bridge=n_apex_bridge, apex_band=apex_band,
        single_direction=direction)


def _build_combined_schedule(
    rbfe, *, construction: str,
    n_windows_half: int, softcore_band: int,
    n_apex_bridge: int, apex_band: float,
    lambda1_rampdown: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Build the COMBINED symmetric schedule (the merge + UWHAM per-state SSOT).

    The combined cntl MUST carry a +1 (dplus) block AND a -1 (dminus) block so
    ``merge_per_direction_outputs`` derives (total, fwd) from its DIRECTION column.

    ``single_core`` -> ``build_rbfe_ladder(single_direction=None)`` (the symmetric
    forward+backward ladder, byte-identical to the pre-existing behaviour).
    ``twocopy`` -> the concatenation of the forward + backward ATS standard
    schedules (forward DIRECTION=+1 block then backward DIRECTION=-1 block), the
    ATS-standard analogue of the symmetric ladder. fwd_count = n_states of one
    direction (e.g. 11 for n_windows_half=6), total = 2 * that.
    """
    if construction != "twocopy":
        return rbfe.build_rbfe_ladder(
            n_windows_half=n_windows_half, softcore_band=softcore_band,
            n_apex_bridge=n_apex_bridge, apex_band=apex_band,
            single_direction=None)

    fwd = rbfe.build_ats_standard_ladder(
        n_windows_half=n_windows_half, single_direction="forward",
        lambda1_rampdown=lambda1_rampdown)
    bwd = rbfe.build_ats_standard_ladder(
        n_windows_half=n_windows_half, single_direction="backward",
        lambda1_rampdown=lambda1_rampdown)
    # Concatenate the two standalone halves (forward +1 block, backward -1 block)
    # — the symmetric combined ladder the merge + UWHAM consume. Per-state arrays
    # are joined; scalar canon (umax/ubcore/acore/temperature) is shared.
    combined: Dict[str, Any] = {
        "lambdas": list(fwd["lambdas"]) + list(bwd["lambdas"]),
        "lambdas_1": list(fwd["lambdas_1"]) + list(bwd["lambdas_1"]),
        "lambdas_2": list(fwd["lambdas_2"]) + list(bwd["lambdas_2"]),
        "directions": list(fwd["directions"]) + list(bwd["directions"]),
        "intermd": list(fwd["intermd"]) + list(bwd["intermd"]),
        "alpha": list(fwd["alpha"]) + list(bwd["alpha"]),
        "u0": list(fwd["u0"]) + list(bwd["u0"]),
        "w0": list(fwd["w0"]) + list(bwd["w0"]),
        "umax": fwd["umax"], "ubcore": fwd["ubcore"], "acore": fwd["acore"],
        "n_states": fwd["n_states"] + bwd["n_states"],
        "n_windows_half": fwd["n_states"],
        "schedule_kind": "ats_standard",
        "single_direction": None,
        "temperature_K": fwd["temperature_K"],
        "schedule_name": "ats_standard_combined_%dw" % (
            fwd["n_states"] + bwd["n_states"],),
        "regime": "ranking_only",
    }
    return combined


def _write_combined_cntl(cntl_path: str, schedule: Dict[str, Any], leg: str,
                         md_steps: int, timestep_fs: float) -> None:
    """Write the COMBINED symmetric-ladder cntl (the per-state SSOT).

    This is the file merge_per_direction_outputs derives (total, fwd) from and
    the UWHAM postprocess parses the per-state (LAMBDA1/LAMBDA2/ALPHA/U0/W0COEFF)
    arrays from (C.1 UWHAM per-state SSOT). It MUST be the symmetric ladder
    (forward dplus states + backward dminus states), so the DIRECTION column has
    both a +1 block and a -1 block (UWHAM requires both legs).
    """
    lines = [
        "# Track B in-place residue-4 RBFE COMBINED asyncre control file",
        "# Generated by trackb_inplace_rbfe_production.py — DO NOT EDIT BY HAND.",
        "# MODE = INPLACE_RBFE: the System XML carries its OWN ATMForce (genuine",
        "# single-shared-core HE1<->methyl swap). The per-direction subdir runs",
        "# (dplus/dminus) are SEPARATE standalone ladders merged here by UWHAM at",
        "# the shared lambda=0.5 apex (the two directions never exchange).",
        "",
        "MODE = 'INPLACE_RBFE'",
        "JOB_TRANSPORT = 'LOCAL_OPENMM'",
        "BASENAME = '%s'" % (JOBNAME,),
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
        "PRODUCTION_STEPS = '%d'" % (md_steps,),
        "TIME_STEP = %s" % (timestep_fs / 1000.0,),
        "",
    ]
    with open(cntl_path, "w") as fh:
        fh.write("\n".join(lines))


def _write_direction_cntl(cntl_path: str, schedule: Dict[str, Any], tag: str,
                          md_steps: int, timestep_fs: float) -> None:
    """Write a single-direction (standalone) cntl for one direction subdir.

    Records the per-state arrays of THAT direction only (DIRECTION all +1 for
    dplus / all -1 for dminus). It is a provenance + reproducibility artifact
    for the standalone run; UWHAM never reads it (UWHAM reads the COMBINED cntl
    after merge). The basename is ``trackb_<tag>`` so the .out files match
    ``r*/trackb_<tag>.out`` (the merge source layout).
    """
    base = JOBNAME + "_" + tag
    lines = [
        "# Track B in-place RBFE STANDALONE single-direction (%s) asyncre cntl"
        % (tag,),
        "# Generated by trackb_inplace_rbfe_production.py — DO NOT EDIT BY HAND.",
        "MODE = 'INPLACE_RBFE'",
        "BASENAME = '%s'" % (base,),
        "TEMPERATURES = '%s'" % (schedule["temperature_K"],),
        "LAMBDAS =      '%s'" % (_csv(schedule["lambdas"]),),
        "DIRECTION =    '%s'" % (_csv(schedule["directions"]),),
        "INTERMEDIATE = '%s'" % (_csv(schedule["intermd"]),),
        "LAMBDA1 =      '%s'" % (_csv(schedule["lambdas_1"]),),
        "LAMBDA2 =      '%s'" % (_csv(schedule["lambdas_2"]),),
        "ALPHA =        '%s'" % (_csv(schedule["alpha"]),),
        "U0 =           '%s'" % (_csv(schedule["u0"]),),
        "W0COEFF =      '%s'" % (_csv(schedule["w0"]),),
        "UMAX = %s" % (schedule["umax"],),
        "UBCORE = %s" % (schedule["ubcore"],),
        "ACORE = %s" % (schedule["acore"],),
        "PRODUCTION_STEPS = '%d'" % (md_steps,),
        "TIME_STEP = %s" % (timestep_fs / 1000.0,),
        "",
    ]
    with open(cntl_path, "w") as fh:
        fh.write("\n".join(lines))


# ---------------------------------------------------------------------------
# C8 SIGN-critical bound-leg decouple-direction gate.
# ---------------------------------------------------------------------------
def gate_decouple_direction(leg: str, decouple_dir, *, raise_on_fail: bool = True
                            ) -> Dict[str, Any]:
    """C8 — the BOUND-leg genuine HE1 decouple direction is SIGN-critical.

    For the BOUND leg the dummy NE1 reference (the HE1's u1 landing point) MUST
    be a finite, unit-magnitude outward vector into LOW density (clearing both
    the receptor AND the binder's own fold). ``compute_decouple_direction``
    derives exactly that (NE1 - centroid(local heavy atoms)); a None / zero /
    non-unit vector means the decouple landed in a clash or could not be
    resolved, which would replace the genuine decouple cost with a packing
    artifact and CORRUPT THE SIGN of ddG_int_bound. So the bound leg fails loud
    if the direction is not a clean unit vector.

    For the FREE leg None is CORRECT by design (any direction is bulk for a free
    peptide; the legacy fixed +Z is used) — the gate passes with reason="free_leg".
    """
    report: Dict[str, Any] = {
        "leg": leg,
        "decouple_dir": (list(decouple_dir) if decouple_dir is not None
                         else None),
        "passed": False,
        "reason": None,
    }
    if leg == "free":
        report["passed"] = True
        report["reason"] = (
            "free_leg: decouple direction None is correct (any direction is "
            "bulk solvent for a free peptide; legacy fixed +Z used)")
        return report

    # BOUND leg — SIGN-critical.
    if decouple_dir is None:
        report["reason"] = (
            "C8 VIOLATION (bound leg): genuine_decouple_dir is None — "
            "compute_decouple_direction could not resolve NE1 or its local "
            "shell. A failed/degenerate decouple direction corrupts the SIGN "
            "of ddG_int_bound. Fix the build before launching the bound leg.")
        if raise_on_fail:
            raise RuntimeError(report["reason"])
        return report
    try:
        comps = [float(c) for c in decouple_dir]
    except (TypeError, ValueError):
        report["reason"] = (
            "C8 VIOLATION (bound leg): genuine_decouple_dir is not numeric "
            "(%r)" % (decouple_dir,))
        if raise_on_fail:
            raise RuntimeError(report["reason"])
        return report
    if len(comps) != 3:
        report["reason"] = (
            "C8 VIOLATION (bound leg): genuine_decouple_dir must be a 3-vector, "
            "got length %d" % (len(comps),))
        if raise_on_fail:
            raise RuntimeError(report["reason"])
        return report
    mag = sum(c * c for c in comps) ** 0.5
    report["magnitude"] = mag
    finite = all(c == c and abs(c) != float("inf") for c in comps)
    if not finite:
        report["reason"] = (
            "C8 VIOLATION (bound leg): genuine_decouple_dir has a non-finite "
            "component (%r)" % (comps,))
        if raise_on_fail:
            raise RuntimeError(report["reason"])
        return report
    if abs(mag - 1.0) > 1e-3:
        report["reason"] = (
            "C8 VIOLATION (bound leg): genuine_decouple_dir is not a unit "
            "vector (|v|=%.6f). compute_decouple_direction returns a "
            "normalized outward vector; a non-unit magnitude means a "
            "degenerate / corrupted direction." % (mag,))
        if raise_on_fail:
            raise RuntimeError(report["reason"])
        return report
    report["passed"] = True
    report["reason"] = (
        "C8 OK (bound leg): genuine_decouple_dir is a finite unit outward "
        "(low-density) vector — the decoupled HE1 clears both the receptor and "
        "the binder fold (NE1 - local-centroid). Sign-safe.")
    return report


# ---------------------------------------------------------------------------
# One standalone single-direction ladder run (in-process, the validated adapter).
# ---------------------------------------------------------------------------
def run_one_direction(
    rbfe,
    loaded: Dict[str, Any],
    *,
    leg: str,
    direction_tag: str,
    subdir: str,
    n_windows_half: int,
    softcore_band: int,
    n_apex_bridge: int,
    apex_band: float,
    n_cycles: int,
    md_steps_per_cycle: int,
    platform_name: str,
    timestep_fs: float,
    rng_seed: int,
    minimize_iters: int,
    backward_equil_steps: int,
    construction: str = "single_core",
    lambda1_rampdown: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Run ONE standalone direction ladder and write its per-walker .out tree.

    Writes ``<subdir>/r*/trackb_<tag>.out`` (UWHAM-consumable) + the
    single-direction cntl + a per-direction driver log (for the mixing gate). C9
    NaN fail-fast: if ANY cycle sees a NaN the run RAISES (no silent UWHAM on a
    corrupt sample set) and records which states NaN'd.

    ``construction`` selects the schedule (single_core -> the validated single-
    shared-core ladder; twocopy -> the canonical ATS standard schedule).

    Returns a per-direction manifest (schedule, mixing gate, per-cycle integrity).
    """
    direction = DIRECTION_OF_TAG[direction_tag]
    os.makedirs(subdir, exist_ok=True)

    schedule = _build_single_direction_schedule(
        rbfe, construction=construction, direction=direction,
        n_windows_half=n_windows_half, softcore_band=softcore_band,
        n_apex_bridge=n_apex_bridge, apex_band=apex_band,
        lambda1_rampdown=lambda1_rampdown)
    base = JOBNAME + "_" + direction_tag

    cntl_path = os.path.join(subdir, base + "_asyncre.cntl")
    _write_direction_cntl(cntl_path, schedule, direction_tag,
                          md_steps_per_cycle, timestep_fs)
    log_path = os.path.join(subdir, base + "_driver.log")

    ladder = rbfe.InplaceRbfeLadder(
        loaded["system"], loaded["positions"], schedule,
        platform_name=platform_name, temperature_K=schedule["temperature_K"],
        timestep_fs=timestep_fs, log_path=log_path, seed=rng_seed,
        minimize_iters=minimize_iters,
        backward_equil_steps=backward_equil_steps,
        out_dir=subdir, out_basename=base)

    nan_any = False
    nan_states_all: set = set()
    per_cycle: List[Dict[str, Any]] = []
    apex_state = schedule["n_states"] - 1   # the lambda=0.5 apex end state
    for _ in range(n_cycles):
        info = ladder.run_cycle(md_steps=md_steps_per_cycle)
        nan_any = nan_any or info["nan_seen"]
        nan_states_all.update(info.get("nan_states", []))
        # C9 data-integrity per cycle: record apex occupancy (is the lambda=0.5
        # end state visited?) + acceptance — NOT just whether something ran.
        rs = info["replica_state"]
        per_cycle.append({
            "cycle": info["cycle"],
            "n_accepted": info["n_accepted"],
            "n_pairs": info["n_pairs"],
            "nan_seen": info["nan_seen"],
            "nan_states": info.get("nan_states", []),
            "apex_occupied": apex_state in rs,
            "endpoint_occupied": 0 in rs,
        })
        if info["nan_seen"]:
            # C9 fail-fast: stop the moment a NaN appears (do not keep stacking
            # corrupt samples). The exact unstable states are recorded.
            ladder.close()
            raise RuntimeError(
                "C9 NaN fail-fast: direction %s/%s NaN'd at cycle %d on "
                "state(s) %s. Soften / bridge those windows + relaunch (the "
                "in-place box is sign-sensitive; a hidden NaN would corrupt "
                "ddG)." % (leg, direction_tag, info["cycle"],
                           sorted(set(info.get("nan_states", [])))))
    ladder.close()

    # Mixing gate — REUSED VERBATIM from the per-direction driver.
    driver = _load_driver()
    tr = driver.parse_state_transitions_from_log(log_path, warmup_cycles=2)
    try:
        mix = driver.check_atm_mixing(
            log_path, schedule_K=schedule["n_states"],
            warmup_cycles=2, min_crossings=1)
    except Exception as exc:  # noqa: BLE001 — gate is advisory in the manifest
        mix = {"error": str(exc), "passed": None}

    adj = tr.get("adjacent_crossings", {})
    n_adjacent_with_crossings = sum(1 for v in adj.values() if v > 0)
    n_adjacent_pairs = schedule["n_states"] - 1

    return {
        "leg": leg,
        "direction_tag": direction_tag,
        "direction": direction,
        "subdir": subdir,
        "cntl_path": cntl_path,
        "driver_log": log_path,
        "platform": ladder.platform_name,
        "construction": construction,
        "schedule_kind": schedule.get("schedule_kind", "rbfe_ladder"),
        "n_states": schedule["n_states"],
        "n_windows_half": schedule["n_windows_half"],
        "n_apex_bridge": schedule.get("n_apex_bridge"),
        "softcore_band": schedule["softcore_band"],
        "n_cycles": n_cycles,
        "md_steps_per_cycle": md_steps_per_cycle,
        "nan_any": nan_any,
        "nan_states": sorted(nan_states_all),
        # C4 / C9 mixing + data integrity (round-trips + apex occupancy — NOT
        # just occupancy; occupancy is not mixing).
        "mixing": {
            "total_round_trips": mix.get("total_round_trips"),
            "both_ends_visited_count": mix.get("both_ends_visited_count"),
            "walls": mix.get("walls"),
            "gate_passed": mix.get("passed"),
            "n_adjacent_pairs": n_adjacent_pairs,
            "n_adjacent_with_crossings": n_adjacent_with_crossings,
            "adjacent_crossings": adj,
        },
        "apex_state": apex_state,
        "apex_occupancy_cycles": sum(1 for c in per_cycle
                                     if c["apex_occupied"]),
        "endpoint_occupancy_cycles": sum(1 for c in per_cycle
                                         if c["endpoint_occupied"]),
        "per_cycle": per_cycle,
    }


# ---------------------------------------------------------------------------
# One replicate = serialize the leg box + run dplus + dminus + merge for UWHAM.
# ---------------------------------------------------------------------------
def run_one_replicate(
    rbfe,
    driver,
    *,
    out_root: str,
    endpoint: str,
    leg: str,
    replicate_index: int,
    seed: str,
    directions: List[str],
    n_windows_half: int,
    softcore_band: int,
    n_apex_bridge: int,
    apex_band: float,
    n_cycles: int,
    md_steps_per_cycle: int,
    platform_name: str,
    timestep_fs: float,
    minimize_iters: int,
    backward_equil_steps: int,
    genuine_decouple_nm: float,
    mtr_ncaa_xml: Optional[str],
    binder_chain: str,
    archive_existing: bool,
    reuse_serialized: bool = False,
    construction: str = "single_core",
    displacement_nm: Optional[float] = None,
    lambda1_rampdown: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Run one matched-seed replicate of one (endpoint, leg): the requested
    direction(s) + (when BOTH ran) merge into the combined leg dir UWHAM consumes.

    The per-replicate Langevin RNG seed is distinct (replicate_index+1 offset)
    so the replicates are genuinely independent (sigma_btwn meaningful).

    ``reuse_serialized`` (box reuse): when True AND a saved box
    (``inplace_rbfe_<leg>_sys.xml`` + ``inplace_rbfe_<leg>.pdb``) already exists
    in the rep dir, the box is NOT re-serialized (no re-solvation) — the existing
    box_A is loaded as-is and the C8 decouple direction is RE-DERIVED from that
    box's own geometry. This is the path for re-running a single direction
    (``--directions dminus``) against the SAME box the producing direction
    (dplus) used, so the two directions' ``.out`` stitch in UWHAM. It is meant
    with ``archive_existing=False`` so the existing dplus subdir is preserved;
    only the requested direction's subdir is (re)written.
    """
    rep_dir = _rep_dir(out_root, endpoint, leg, replicate_index)
    # C12 R-7: archive (never delete) an existing rep dir before a fresh run.
    # Box reuse keeps the rep dir in place (do NOT archive — that would move the
    # producing direction's outputs away); the reuse path requires the box, and
    # is meant to be paired with archive_existing=False (see launcher gate).
    if archive_existing and os.path.isdir(rep_dir) and os.listdir(rep_dir):
        ts = time.strftime("%Y%m%d_%H%M%S")
        archive_root = os.path.join(out_root, "_archive",
                                    "%s_%s_rep%d_%s"
                                    % (endpoint, leg, replicate_index, ts))
        os.makedirs(os.path.dirname(archive_root), exist_ok=True)
        shutil.move(rep_dir, archive_root)
    os.makedirs(rep_dir, exist_ok=True)

    # 1) Obtain the in-place fused System for this leg. Two paths:
    #    (a) FRESH (default): serialize a new box (C1: each leg gets its OWN
    #        in-place RBFE box; the build picks free vs bound topology + the
    #        bound decouple direction). The serialized System carries the genuine
    #        ATMForce (so the in-process adapter never rebuilds it).
    #    (b) REUSE (--reuse-serialized, box present): load the EXISTING box_A
    #        without re-solvating, and re-derive the C8 decouple direction from
    #        that box's geometry — so a single-direction re-run samples the SAME
    #        box the producing direction used (UWHAM stitch validity).
    reuse_xml = os.path.join(rep_dir, "inplace_rbfe_%s_sys.xml" % (leg,))
    reuse_pdb = os.path.join(rep_dir, "inplace_rbfe_%s.pdb" % (leg,))
    reused_box = (reuse_serialized
                  and os.path.isfile(reuse_xml) and os.path.isfile(reuse_pdb))
    if reused_box:
        loaded = rbfe.load_serialized_system(reuse_xml, reuse_pdb)
        decouple_dir = rbfe.recompute_decouple_direction_from_loaded(
            loaded, binder_chain=binder_chain,
            decouple_nm=genuine_decouple_nm)
        # C8 SIGN-critical: validate the RE-DERIVED bound-leg direction (fail
        # loud). For the bound leg this must be a finite unit outward vector.
        c8 = gate_decouple_direction(leg, decouple_dir, raise_on_fail=True)
        c8["source"] = "reused_box_recomputed"
        ser = {
            "sys_xml_path": reuse_xml,
            "pdb_path": reuse_pdb,
            "n_atoms": loaded["n_atoms"],
            "genuine_decouple_dir": decouple_dir,
            "reused": True,
        }
    else:
        if reuse_serialized:
            # Asked to reuse but the box is absent — fail loud (do NOT silently
            # fall back to a fresh re-serialize, which would build a DIFFERENT
            # box and break the same-box UWHAM stitch the flag exists to keep).
            raise RuntimeError(
                "run_one_replicate(reuse_serialized=True): the saved box is "
                "missing for %s/%s rep%d — expected\n  %s\n  %s\nA single-"
                "direction reuse run REQUIRES the producing direction's box_A "
                "(re-serializing would build a different box and break the "
                "UWHAM stitch). Stage the box first." % (
                    endpoint, leg, replicate_index, reuse_xml, reuse_pdb))
        serialize_kwargs = dict(
            leg=leg, out_dir=rep_dir, seed=seed, binder_chain=binder_chain,
            solvate=True, harmonize_common_charges=False, swap_mode="genuine",
            genuine_decouple_nm=genuine_decouple_nm, mtr_ncaa_xml=mtr_ncaa_xml,
            constraints=None, tag=leg, construction=construction)
        # displacement_nm is two-copy-only; pass it through only when set so the
        # single-core path signature is unaffected (the rbfe serialize uses the
        # canonical ATS default when None).
        if construction == "twocopy" and displacement_nm is not None:
            serialize_kwargs["displacement_nm"] = displacement_nm
        ser = rbfe.serialize_inplace_rbfe_system(**serialize_kwargs)
        # C8 SIGN-critical: validate the bound-leg decouple direction (fail loud).
        c8 = gate_decouple_direction(leg, ser.get("genuine_decouple_dir"),
                                     raise_on_fail=True)
        c8["source"] = "fresh_serialize"
        loaded = rbfe.load_serialized_system(ser["sys_xml_path"], ser["pdb_path"])

    # 2) Run each requested direction as a STANDALONE ladder. The per-replicate
    #    velocity/exchange RNG seed differs per replicate (genuine independence).
    rng_seed = 20260613 + 1000 * (replicate_index + 1)
    leg_dir = os.path.join(out_root, endpoint, leg, "rep%d" % (replicate_index,))
    per_direction: Dict[str, Any] = {}
    for tag in directions:
        subdir = os.path.join(leg_dir, tag)
        per_direction[tag] = run_one_direction(
            rbfe, loaded,
            leg=leg, direction_tag=tag, subdir=subdir,
            n_windows_half=n_windows_half, softcore_band=softcore_band,
            n_apex_bridge=n_apex_bridge, apex_band=apex_band,
            n_cycles=n_cycles, md_steps_per_cycle=md_steps_per_cycle,
            platform_name=platform_name, timestep_fs=timestep_fs,
            rng_seed=rng_seed, minimize_iters=minimize_iters,
            backward_equil_steps=backward_equil_steps,
            construction=construction, lambda1_rampdown=lambda1_rampdown)

    # 3) Write the COMBINED symmetric cntl (the SSOT for merge + UWHAM) +, when
    #    BOTH directions ran, merge the per-direction outputs into r*/trackb.out.
    merged = None
    if set(directions) == set(DIRECTION_TAGS):
        combined_schedule = _build_combined_schedule(
            rbfe, construction=construction,
            n_windows_half=n_windows_half, softcore_band=softcore_band,
            n_apex_bridge=n_apex_bridge, apex_band=apex_band,
            lambda1_rampdown=lambda1_rampdown)
        combined_cntl = os.path.join(leg_dir, JOBNAME + "_asyncre.cntl")
        _write_combined_cntl(combined_cntl, combined_schedule, leg,
                             md_steps_per_cycle, timestep_fs)
        # merge_per_direction_outputs derives (total, fwd) from the combined cntl
        # DIRECTION column and shifts dminus stateids up by fwd_count.
        merged = driver.merge_per_direction_outputs(
            leg_dir, jobname=JOBNAME)

    # C3 apex co-sampling check: both directions occupied their lambda=0.5 apex.
    apex_cosampled = None
    if set(directions) == set(DIRECTION_TAGS):
        dp = per_direction["dplus"]["apex_occupancy_cycles"]
        dm = per_direction["dminus"]["apex_occupancy_cycles"]
        apex_cosampled = (dp > 0 and dm > 0)

    manifest = {
        "endpoint": endpoint,
        "leg": leg,
        "construction": construction,
        "replicate_index": replicate_index,
        "seed": seed,
        "rep_dir": leg_dir,
        "rng_seed": rng_seed,
        "directions": list(directions),
        "c8_decouple_gate": c8,
        "reused_box": reused_box,
        "serialize": {
            "sys_xml_path": ser["sys_xml_path"],
            "pdb_path": ser["pdb_path"],
            "n_atoms": ser["n_atoms"],
            "genuine_decouple_dir": ser.get("genuine_decouple_dir"),
            "reused": bool(ser.get("reused")),
            "construction": ser.get("construction", construction),
        },
        "per_direction": per_direction,
        "merged": merged,
        "apex_cosampled_C3": apex_cosampled,
    }
    with open(os.path.join(leg_dir, "run_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    return manifest


# ---------------------------------------------------------------------------
# Leg run (n replicates) + the empirical C6 window-count escalation.
# ---------------------------------------------------------------------------
def run_leg(
    *,
    out_root: str,
    endpoint: str,
    leg: str,
    seeds: List[str],
    directions: List[str],
    n_windows_half: int,
    softcore_band: int,
    n_apex_bridge: int,
    apex_band: float,
    n_cycles: int,
    md_steps_per_cycle: int,
    platform_name: str,
    timestep_fs: float,
    minimize_iters: int,
    backward_equil_steps: int,
    genuine_decouple_nm: float,
    mtr_ncaa_xml: Optional[str],
    binder_chain: str,
    archive_existing: bool,
    reuse_serialized: bool = False,
    construction: str = "single_core",
    displacement_nm: Optional[float] = None,
    lambda1_rampdown: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Run all matched-seed replicates of one (endpoint, leg)."""
    rbfe = _load_rbfe()
    driver = _load_driver()
    reps: List[Dict[str, Any]] = []
    for j, seed in enumerate(seeds):
        reps.append(run_one_replicate(
            rbfe, driver,
            out_root=out_root, endpoint=endpoint, leg=leg,
            replicate_index=j, seed=seed, directions=directions,
            n_windows_half=n_windows_half, softcore_band=softcore_band,
            n_apex_bridge=n_apex_bridge, apex_band=apex_band,
            n_cycles=n_cycles, md_steps_per_cycle=md_steps_per_cycle,
            platform_name=platform_name, timestep_fs=timestep_fs,
            minimize_iters=minimize_iters,
            backward_equil_steps=backward_equil_steps,
            genuine_decouple_nm=genuine_decouple_nm,
            mtr_ncaa_xml=mtr_ncaa_xml, binder_chain=binder_chain,
            archive_existing=archive_existing,
            reuse_serialized=reuse_serialized,
            construction=construction, displacement_nm=displacement_nm,
            lambda1_rampdown=lambda1_rampdown))
    return {
        "endpoint": endpoint,
        "leg": leg,
        "construction": construction,
        "n_replicates": len(reps),
        "replicates": reps,
    }


# ---------------------------------------------------------------------------
# C11 pre-registration of the 4 outcomes.
# ---------------------------------------------------------------------------
PRE_REGISTERED_OUTCOMES = [
    {
        "id": "O1_favorable_sign_resolved",
        "criteria": ("ddG_int (free or, after the bound leg, ddG_bind) sign is "
                     "RESOLVED (z_SE >= 3, n >= 3) AND NEGATIVE — consistent with "
                     "the favorable_cp4 Magotti SSOT sign (Cp4 binds tighter)."),
    },
    {
        "id": "O2_discordant_sign_resolved",
        "criteria": ("ddG sign is RESOLVED (z_SE >= 3, n >= 3) but POSITIVE — "
                     "discordant with the favorable_cp4 sign. SIGN-only; this is "
                     "a genuine finding, NOT a magnitude claim."),
    },
    {
        "id": "O3_sign_undetermined",
        "criteria": ("z_SE < 3 OR n < 3 — the sign is NOT resolved. The point "
                     "estimate is reported with the quadrature SE; NO sign claim "
                     "(R-18 honest: undetermined, not assumed favorable)."),
    },
    {
        "id": "O4_non_convergent",
        "criteria": ("A gate FAILED: C9 NaN, C4 zero round-trips on a direction, "
                     "C6 collapse (adjacent O < 0.03) un-fixable within the "
                     "escalation cap, or C3 apex not co-sampled. ddG is NOT "
                     "reported as converged."),
    },
]


def write_pre_registration(out_root: str, args_dict: Dict[str, Any]) -> str:
    """Write the C11 pre-registration JSON (the 4 outcomes + the launch config),
    BEFORE any analysis decides which one obtains."""
    os.makedirs(out_root, exist_ok=True)
    path = os.path.join(out_root, "pre_registration.json")
    payload = {
        "schema": "trackb_inplace_rbfe_prereg_v1",
        "regime": "ranking_only_R11",
        "registered_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "sign_claim_policy": (
            "SIGN is claimed ONLY if z_SE >= %.0f (Welch-Satterthwaite "
            "quadrature, NOT paired) AND n >= %d matched seeds. Otherwise the "
            "outcome is O3 (undetermined). Magotti SSOT comparison is SIGN "
            "ONLY, never magnitude (R-11 / R-18)."
            % (Z_SE_SIGN_THRESHOLD, MIN_REPLICATES_FOR_SIGN)),
        "outcomes": PRE_REGISTERED_OUTCOMES,
        "config": args_dict,
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# UWHAM analysis of a completed leg cohort (n replicates) -> ddG_int + gates.
# ---------------------------------------------------------------------------
def _welch_satterthwaite_quadrature(
    mean_a: float, se_a: float, mean_b: float, se_b: float
) -> Dict[str, Any]:
    """Quadrature SE of (mean_a - mean_b) with independent endpoint SEs (C7).

    NOT paired — the prompt's C7 explicitly requires the quadrature (Welch-
    Satterthwaite) SE for any SIGN claim. ddG = mean_a - mean_b,
    SE = sqrt(se_a^2 + se_b^2), z_SE = |ddG| / SE. The sign is claimed only when
    z_SE >= Z_SE_SIGN_THRESHOLD.
    """
    ddg = mean_a - mean_b
    se = (se_a ** 2 + se_b ** 2) ** 0.5
    z = (abs(ddg) / se) if se > 0 else None
    return {"ddg": ddg, "se_quadrature": se, "z_se": z}


def analyze_leg_cohort(
    out_root: str,
    leg: str,
    *,
    n_replicates: int,
    mintimeid: Optional[int],
    output_json: Optional[str] = None,
) -> Dict[str, Any]:
    """UWHAM analysis of a completed (cp4 + wt) leg cohort -> ddG_int.

    Consumes the merged per-replicate leg dirs via the EXISTING UWHAM
    postprocess (analyze_replicate_set + build_v3_payload). Reports the
    quadrature SE + z_SE (C7) and the per-direction mixing / apex / overlap gate
    status; a SIGN is claimed ONLY when z_SE >= threshold AND n >= 3 (C5/C7).
    """
    uwham = _load_uwham()

    cp4_dirs = [_rep_dir(out_root, "cp4", leg, j) for j in range(n_replicates)]
    wt_dirs = [_rep_dir(out_root, "wt", leg, j) for j in range(n_replicates)]

    cp4 = uwham.analyze_replicate_set(
        cp4_dirs, jobname=JOBNAME, mintimeid=mintimeid, maxtimeid=None,
        require_cntl_schedule=True)
    wt = uwham.analyze_replicate_set(
        wt_dirs, jobname=JOBNAME, mintimeid=mintimeid, maxtimeid=None,
        require_cntl_schedule=True)

    payload = uwham.build_v3_payload(
        cp4=cp4, wt=wt, protocol="inplace_rbfe_%s" % (leg,),
        mintimeid=mintimeid, maxtimeid=None,
        notes=("In-place residue-4 RBFE per-direction-separate estimator "
               "(single-shared-core HE1<->methyl swap). ddint_%s = "
               "mean(cp4.dgbind1) - mean(wt.dgbind1)." % (leg,)))

    # C7 quadrature SE + z_SE (the SIGN gate; NOT the paired SEM). Endpoint SE =
    # sigma_btwn / sqrt(n) (the inter-replicate SEM already in analyze_replicate_set).
    cp4_se = cp4.get("sem_dgbind1_kcal")
    wt_se = wt.get("sem_dgbind1_kcal")
    n_cp4 = cp4.get("n_replicates")
    n_wt = wt.get("n_replicates")
    n_matched = min(n_cp4, n_wt) if isinstance(n_cp4, int) and isinstance(n_wt, int) else 0

    quad: Optional[Dict[str, Any]] = None
    sign_status = "undetermined"
    sign_reason = None
    if cp4_se is not None and wt_se is not None:
        quad = _welch_satterthwaite_quadrature(
            cp4["mean_dgbind1_kcal"], cp4_se,
            wt["mean_dgbind1_kcal"], wt_se)
        z = quad["z_se"]
        if n_matched < MIN_REPLICATES_FOR_SIGN:
            sign_status = "undetermined"
            sign_reason = ("O3: n=%d < %d matched seeds (C5) — sign not claimed"
                           % (n_matched, MIN_REPLICATES_FOR_SIGN))
        elif z is None or z < Z_SE_SIGN_THRESHOLD:
            sign_status = "undetermined"
            sign_reason = ("O3: z_SE=%s < %.0f (C7 quadrature) — sign not claimed"
                           % (("%.2f" % z) if z is not None else "n/a",
                              Z_SE_SIGN_THRESHOLD))
        else:
            sign_status = "negative" if quad["ddg"] < 0 else "positive"
            sign_reason = ("z_SE=%.2f >= %.0f AND n=%d >= %d — sign RESOLVED (%s)"
                           % (z, Z_SE_SIGN_THRESHOLD, n_matched,
                              MIN_REPLICATES_FOR_SIGN, sign_status))
    else:
        sign_reason = ("O3: sigma_btwn unavailable (need >= 2 replicates per "
                       "endpoint for an inter-replicate SE) — sign not claimed")

    # Gate digest (C3 / C4 / C6) read from the per-replicate run manifests.
    gate_digest = _collect_gate_digest(out_root, leg, n_replicates)

    result = {
        "schema": "trackb_inplace_rbfe_ddint_v1",
        "regime": "ranking_only_R11",
        "leg": leg,
        "n_matched_replicates": n_matched,
        "ddint_%s_kcal" % (leg,): payload["ddint_free_kcal"]
            if leg == "free" else (cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"]),
        "ddint_point_estimate_kcal": cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"],
        "cp4_mean_dgbind1_kcal": cp4["mean_dgbind1_kcal"],
        "wt_mean_dgbind1_kcal": wt["mean_dgbind1_kcal"],
        "cp4_sem_kcal": cp4_se,
        "wt_sem_kcal": wt_se,
        # C7 — the SIGN gate uses the QUADRATURE SE / z_SE (NOT the paired SEM;
        # the paired stats remain in the payload as a separate diagnostic).
        "quadrature_C7": quad,
        "sign_status": sign_status,
        "sign_reason": sign_reason,
        "sign_undetermined": sign_status == "undetermined",
        "gates": gate_digest,
        # The full UWHAM v3 payload (paired stats, closure quarantine, overlap
        # QC, sign-vs-anchor advisory) for transparency.
        "uwham_v3_payload": payload,
        "magotti_policy": ("ranking-only: SIGN comparison with Magotti 2009 "
                           "ITC/SPR allowed (Cp4 better) — MAGNITUDE FORBIDDEN "
                           "(pre-v0.7-calibration). Free leg alone is NOT "
                           "ddG_bind."),
    }
    if output_json:
        with open(output_json, "w") as fh:
            json.dump(result, fh, indent=2, default=str)
    return result


def _collect_gate_digest(out_root: str, leg: str, n_replicates: int
                         ) -> Dict[str, Any]:
    """Read the per-replicate run manifests and digest C3/C4/C6 + NaN status."""
    rows: List[Dict[str, Any]] = []
    all_c4_ok = True
    all_c3_ok = True
    any_nan = False
    for endpoint in ENDPOINTS:
        for j in range(n_replicates):
            mpath = os.path.join(_rep_dir(out_root, endpoint, leg, j),
                                 "run_manifest.json")
            if not os.path.isfile(mpath):
                rows.append({"endpoint": endpoint, "rep": j,
                             "status": "manifest_missing"})
                all_c4_ok = False
                continue
            with open(mpath) as fh:
                m = json.load(fh)
            for tag, pd in m.get("per_direction", {}).items():
                rt = pd.get("mixing", {}).get("total_round_trips")
                rt = rt if rt is not None else 0
                nan = pd.get("nan_any")
                any_nan = any_nan or bool(nan)
                if rt < 1:
                    all_c4_ok = False
                rows.append({
                    "endpoint": endpoint, "rep": j, "direction": tag,
                    "total_round_trips_C4": rt,
                    "nan_any_C9": bool(nan),
                    "n_adjacent_with_crossings": pd.get("mixing", {}).get(
                        "n_adjacent_with_crossings"),
                    "n_adjacent_pairs": pd.get("mixing", {}).get(
                        "n_adjacent_pairs"),
                })
            if m.get("apex_cosampled_C3") is False:
                all_c3_ok = False
    return {
        "C3_apex_cosampled_all": all_c3_ok,
        "C4_roundtrip_ge_1_all": all_c4_ok,
        "C9_nan_seen": any_nan,
        "per_run": rows,
    }


# ---------------------------------------------------------------------------
# VRAM-aware concurrent pool (maximize hardware utilization — no idle GPU).
#
# Each POOL UNIT is one (endpoint, leg, replicate): the worker serializes the
# in-place box ONCE then runs dplus + dminus sequentially + merges (the
# expensive System build is amortized over both directions). A unit is launched
# as a DETACHED subprocess re-invoking this script in --worker mode pinned to a
# GPU via CUDA_VISIBLE_DEVICES. The dispatcher packs concurrent units per GPU to
# the GPU's free VRAM (with headroom), so both GPUs stay saturated.
#
# THREAD-PIN (REAL utilization, not thrashing): every worker inherits the BLAS /
# OMP single-thread pin (OPENBLAS/MKL/NUMEXPR/VECLIB = 1) so N concurrent
# OpenMM + numpy workers do not oversubscribe the 16 HW threads. OpenMM's GPU
# kernels are the compute; the CPU side is light per worker.
#
# DUAL-GPU: the free leg (~1-2 GB/unit, 4252-atom box) packs many units on the
# 5070Ti (device 0). The bound leg (~6 GB/unit, 541k-atom box) is V100-class
# (device dispatch + the V100 SSH path) and is BUILT here but stays behind the
# user gate (NOT auto-launched) — the user explicitly gates the multi-day bound
# run after the free ddG_int_free preview lands.
# ---------------------------------------------------------------------------
_FREE_UNIT_VRAM_GB = 2.0       # generous per-free-unit VRAM budget (observed ~1-2)
_BOUND_UNIT_VRAM_GB = 6.5      # per-bound-unit (541k box, ~6 GB + headroom)
_GPU_VRAM_HEADROOM_GB = 1.5    # leave this free on each GPU


def _query_local_gpu_free_gb(device_index: int = 0) -> Optional[float]:
    """Free VRAM (GiB) on a local CUDA device via nvidia-smi, or None."""
    import subprocess
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.free",
             "--format=csv,noheader,nounits", "-i", str(device_index)],
            stderr=subprocess.DEVNULL, timeout=20).decode().strip()
        return float(out.splitlines()[0]) / 1024.0
    except Exception:
        return None


def _max_concurrent_units(free_gb: float, unit_gb: float) -> int:
    """How many units fit in ``free_gb`` (with headroom), at least 1."""
    usable = max(0.0, free_gb - _GPU_VRAM_HEADROOM_GB)
    return max(1, int(usable // unit_gb))


def _worker_env(device_index: int) -> Dict[str, str]:
    """Subprocess env: pin to one GPU + single-thread BLAS/OMP (no thrash)."""
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(device_index)
    env.setdefault("OMP_NUM_THREADS", "2")
    for k in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[k] = "1"
    return env


def run_pool_local(
    *,
    out_root: str,
    leg: str,
    endpoints: List[str],
    seeds: List[str],
    directions: List[str],
    device_index: int,
    max_concurrent: int,
    ladder_args: Dict[str, Any],
    archive_existing: bool,
    reuse_serialized: bool = False,
) -> Dict[str, Any]:
    """Launch all (endpoint, seed) units of one leg as a CONCURRENT local pool.

    Packs up to ``max_concurrent`` detached worker subprocesses on
    ``device_index`` at once; as each finishes the next is started, so the GPU
    stays saturated until every unit completes. Returns the per-unit launch +
    exit manifest. Each worker writes its own per-replicate run_manifest.json +
    log, so the merge + UWHAM path consumes the result regardless of pool
    ordering.
    """
    import subprocess

    units: List[Dict[str, Any]] = []
    for ep in endpoints:
        for j, seed in enumerate(seeds):
            units.append({"endpoint": ep, "replicate_index": j, "seed": seed})

    logdir = os.path.join(out_root, "_pool_logs")
    os.makedirs(logdir, exist_ok=True)
    env = _worker_env(device_index)

    running: List[Tuple[Dict[str, Any], Any, Any]] = []   # (unit, proc, logfh)
    done: List[Dict[str, Any]] = []
    queue = list(units)

    def _launch(unit):
        log_path = os.path.join(
            logdir, "worker_%s_%s_rep%d.log"
            % (unit["endpoint"], leg, unit["replicate_index"]))
        cmd = [
            sys.executable, os.path.abspath(__file__),
            "--worker",
            "--worker-endpoint", unit["endpoint"],
            "--worker-replicate", str(unit["replicate_index"]),
            "--worker-seed", unit["seed"],
            "--leg", leg,
            "--out-root", out_root,
            "--directions", ",".join(directions),
            "--n-windows-half", str(ladder_args["n_windows_half"]),
            "--softcore-band", str(ladder_args["softcore_band"]),
            "--n-apex-bridge", str(ladder_args["n_apex_bridge"]),
            "--apex-band", str(ladder_args["apex_band"]),
            "--n-cycles", str(ladder_args["n_cycles"]),
            "--md-steps-per-cycle", str(ladder_args["md_steps_per_cycle"]),
            "--platform", ladder_args["platform"],
            "--timestep-fs", str(ladder_args["timestep_fs"]),
            "--minimize-iters", str(ladder_args["minimize_iters"]),
            "--backward-equil-steps", str(ladder_args["backward_equil_steps"]),
            "--genuine-decouple-nm", str(ladder_args["genuine_decouple_nm"]),
            "--binder-chain", ladder_args["binder_chain"],
        ]
        if ladder_args.get("mtr_ncaa_xml"):
            cmd += ["--mtr-ncaa-xml", ladder_args["mtr_ncaa_xml"]]
        # Two-copy construction (opt-in) — propagate the flag + displacement so
        # the worker builds the SAME box the dispatcher requested.
        if ladder_args.get("construction") == "twocopy":
            cmd.append("--twocopy")
            if ladder_args.get("displacement_nm") is not None:
                cmd += ["--displacement-nm",
                        str(ladder_args["displacement_nm"])]
            # Leg-down densification knots (two-copy-only): thread the same
            # comma-list the dispatcher parsed so the worker builds the SAME
            # ladder (e.g. the λ1=0.05 leg-switch bridge => 12 λ/leg).
            if ladder_args.get("lambda1_rampdown") is not None:
                cmd += ["--lambda1-rampdown",
                        ",".join(str(x) for x in ladder_args["lambda1_rampdown"])]
        if not archive_existing:
            cmd.append("--no-archive-existing")
        if reuse_serialized:
            cmd.append("--reuse-serialized")
        logfh = open(log_path, "w")
        proc = subprocess.Popen(cmd, stdout=logfh, stderr=subprocess.STDOUT,
                                env=env, cwd=_PROJ)
        unit["pid"] = proc.pid
        unit["log"] = log_path
        unit["device_index"] = device_index
        print("  [pool] launched %s/%s rep%d seed=%s -> PID %d (GPU %d) log=%s"
              % (unit["endpoint"], leg, unit["replicate_index"], unit["seed"],
                 proc.pid, device_index, log_path))
        return (unit, proc, logfh)

    import time as _t
    while queue or running:
        while queue and len(running) < max_concurrent:
            running.append(_launch(queue.pop(0)))
        # Poll the running set.
        still: List[Tuple[Dict[str, Any], Any, Any]] = []
        for (unit, proc, logfh) in running:
            rc = proc.poll()
            if rc is None:
                still.append((unit, proc, logfh))
            else:
                logfh.close()
                unit["exit_code"] = rc
                done.append(unit)
                print("  [pool] finished %s/%s rep%d -> exit %d"
                      % (unit["endpoint"], leg, unit["replicate_index"], rc))
        running = still
        if queue or running:
            _t.sleep(5)

    return {
        "leg": leg, "device_index": device_index,
        "max_concurrent": max_concurrent,
        "n_units": len(units), "units": done,
        "all_ok": all(u.get("exit_code") == 0 for u in done),
    }


# ---------------------------------------------------------------------------
# CLI / dispatch.
# ---------------------------------------------------------------------------
def _parse_seeds(raw: Optional[str]) -> List[str]:
    if not raw:
        return list(DEFAULT_SEEDS)
    seeds = [s.strip() for s in raw.split(",") if s.strip()]
    if not seeds:
        raise ValueError("--seeds parsed to an empty list")
    return seeds


def _parse_directions(raw: Optional[str]) -> List[str]:
    if not raw:
        return list(DIRECTION_TAGS)
    dirs = [d.strip() for d in raw.split(",") if d.strip()]
    for d in dirs:
        if d not in DIRECTION_TAGS:
            raise ValueError("unknown direction %r (expected dplus/dminus)" % (d,))
    return dirs


def _parse_lambda1_rampdown(raw: Optional[str]) -> Optional[List[float]]:
    """Parse the --lambda1-rampdown comma list into floats (fail loud).

    ``None`` / empty -> ``None`` (canonical uniform leg-down). The per-knot
    range / monotonicity / endpoint validation lives in build_ats_standard_ladder
    (the SSOT); this only turns the CLI string into floats.
    """
    if not raw:
        return None
    knots = [s.strip() for s in raw.split(",") if s.strip()]
    if not knots:
        raise ValueError("--lambda1-rampdown parsed to an empty list")
    try:
        return [float(x) for x in knots]
    except ValueError as exc:
        raise ValueError(
            "--lambda1-rampdown must be a comma list of floats, got %r (%s)"
            % (raw, exc))


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Track B in-place residue-4 RBFE per-direction production "
                    "launcher (ranking-only, R-11).")
    p.add_argument("--leg", choices=["free", "bound"], default="free",
                   help="Thermodynamic leg to run (default free — the EARLY "
                        "signal on the 5070Ti; bound is the multi-day follow-up).")
    p.add_argument("--endpoints", default="cp4,wt",
                   help="Comma list of endpoints (default cp4,wt — matched-seed "
                        "paired ddG).")
    p.add_argument("--seeds", default=None,
                   help="Comma list of matched QM-snapshot seeds (default "
                        "s7,s101,s127 => n=3, all present for BOTH endpoints; "
                        "C5 needs n>=3 for a sign claim).")
    p.add_argument("--directions", default=None,
                   help="Comma list of direction tags (default dplus,dminus = "
                        "both, required for the merge + UWHAM cycle).")
    p.add_argument("--out-root", default=None,
                   help="Output root (default outputs/_trackb/inplace_rbfe_prod).")
    p.add_argument("--twocopy", action="store_true",
                   help="Use the CANONICAL ATS TWO-COPY box (copy-1 MTR at site + "
                        "copy-2 WT displaced into bulk, common-coord swap, real "
                        "transfer) instead of the legacy single-shared-core box. "
                        "Selects construction='twocopy' AND the canonical ATS "
                        "standard schedule (spec C: λ1=0/λ2 climbs then "
                        "λ2=0.5/λ1 climbs; n_windows_half=6 => 11 λ/leg). DEFAULT "
                        "OFF = single_core (legacy, byte-identical).")
    p.add_argument("--displacement-nm", type=float, default=None,
                   help="Two-copy ONLY: magnitude of the copy-2 bulk displacement "
                        "d (nm; ATS peptide convention ~4.0 = 40 A). Default None "
                        "=> the canonical ATS_TWOCOPY_DISPLACEMENT_NM.")
    p.add_argument("--lambda1-rampdown", default=None,
                   help="Two-copy ONLY: comma list of explicit leg-down λ1 knots "
                        "to densify the leg-switch handoff (each in (0,0.5], "
                        "strictly increasing, ending at 0.5; the apex λ1=0 state "
                        "is placed by the leg-up phase). E.g. "
                        "0.05,0.1,0.2,0.3,0.4,0.5 inserts a λ1=0.05 bridge window "
                        "=> 12 λ/leg. Default None => the canonical uniform "
                        "leg-down (11 λ/leg). Interior-λ reshaping = ΔG-unbiased "
                        "(ranking-only, R-11); soft-core canon untouched (C7).")
    p.add_argument("--n-windows-half", type=int, default=6,
                   help="Forward-half window count (the fork-validated default "
                        "6; C6 escalation densifies thin pairs up to the cap). "
                        "For --twocopy this is the ATS leg-up phase length "
                        "(n=6 => 11 λ states per leg, the spec C count).")
    p.add_argument("--softcore-band", type=int, default=2)
    p.add_argument("--n-apex-bridge", type=int, default=0)
    p.add_argument("--apex-band", type=float, default=0.5)
    p.add_argument("--n-cycles", type=int, default=DEFAULT_N_CYCLES,
                   help="Production asyncre cycles per direction (default %d; "
                        "fork needed ~50 for the bound leg to round-trip, so "
                        "production runs many round-trips)." % DEFAULT_N_CYCLES)
    p.add_argument("--md-steps-per-cycle", type=int,
                   default=DEFAULT_MD_STEPS_PER_CYCLE)
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--timestep-fs", type=float, default=1.0)
    p.add_argument("--minimize-iters", type=int, default=500)
    p.add_argument("--backward-equil-steps", type=int, default=500)
    p.add_argument("--genuine-decouple-nm", type=float, default=1.2)
    p.add_argument("--mtr-ncaa-xml", default=None,
                   help="Optional S0-harmonized RBFE hybrid XML (MTR).")
    p.add_argument("--binder-chain", default="B")
    p.add_argument("--mintimeid", type=int, default=100,
                   help="UWHAM equilibration-discard (default 100 cycles; "
                        "-1 = include all).")
    p.add_argument("--no-archive-existing", action="store_true",
                   help="Do NOT archive (R-7) an existing rep dir before a "
                        "fresh run (default archives).")
    p.add_argument("--reuse-serialized", action="store_true",
                   help="Box reuse: if a saved box "
                        "(inplace_rbfe_<leg>_sys.xml + .pdb) already exists in a "
                        "rep dir, LOAD it instead of re-serializing (no re-"
                        "solvation), and RE-DERIVE the C8 decouple direction "
                        "from that box. Use with --no-archive-existing and a "
                        "single --directions value to re-run one direction "
                        "against the SAME box the producing direction used "
                        "(UWHAM stitch validity). Fails loud if the box is "
                        "absent (no silent fresh re-serialize).")
    p.add_argument("--analyze-only", action="store_true",
                   help="Skip the run; UWHAM-analyze already-completed leg dirs.")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the plan (legs x directions x replicates) + write "
                        "the pre-registration; run NOTHING.")
    # Concurrent pool (maximize hardware utilization — no idle GPU).
    p.add_argument("--pool", action="store_true",
                   help="Launch all (endpoint, seed) units of the leg as a "
                        "VRAM-aware CONCURRENT pool on one local GPU "
                        "(--device-index), packed to free VRAM. Each unit is a "
                        "detached worker subprocess. UWHAM runs after the pool "
                        "drains.")
    p.add_argument("--device-index", type=int, default=0,
                   help="Local CUDA device index for --pool (default 0 = 5070Ti).")
    p.add_argument("--max-concurrent", type=int, default=0,
                   help="Override pool concurrency (0 = auto from free VRAM).")
    # Worker mode (internal — one (endpoint, leg, replicate) unit; both
    # directions + merge). NOT for direct human use; the pool re-invokes this.
    p.add_argument("--worker", action="store_true",
                   help=argparse.SUPPRESS)
    p.add_argument("--worker-endpoint", default=None, help=argparse.SUPPRESS)
    p.add_argument("--worker-replicate", type=int, default=None,
                   help=argparse.SUPPRESS)
    p.add_argument("--worker-seed", default=None, help=argparse.SUPPRESS)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    out_root = args.out_root or os.path.join(
        _PROJ, "outputs", "_trackb", "inplace_rbfe_prod")
    out_root = os.path.abspath(out_root)
    endpoints = [e.strip() for e in args.endpoints.split(",") if e.strip()]
    seeds = _parse_seeds(args.seeds)
    directions = _parse_directions(args.directions)
    leg = args.leg
    construction = "twocopy" if args.twocopy else "single_core"
    displacement_nm = args.displacement_nm
    lambda1_rampdown = _parse_lambda1_rampdown(args.lambda1_rampdown)

    mintimeid = None if args.mintimeid == -1 else args.mintimeid

    # Box-reuse safety gate: --reuse-serialized MUST be paired with
    # --no-archive-existing. Otherwise run_one_replicate archives (shutil.move)
    # the existing rep dir BEFORE the reuse load could find the box — which both
    # destroys the producing direction's outputs AND removes the box the reuse
    # path needs. Fail loud rather than silently moving box_A + dplus away.
    if args.reuse_serialized and not args.no_archive_existing:
        print("ERROR: --reuse-serialized requires --no-archive-existing "
              "(otherwise the existing rep dir — box_A + the producing "
              "direction's outputs — is archived/moved before the reuse load).",
              file=sys.stderr)
        return 2

    # WORKER MODE: run exactly ONE (endpoint, leg, replicate) unit (both
    # directions + merge) and exit. The pool dispatcher re-invokes this with a
    # GPU pinned via CUDA_VISIBLE_DEVICES. NaN fail-fast (C9) propagates as a
    # non-zero exit so the pool records the failure.
    if args.worker:
        if args.worker_endpoint is None or args.worker_replicate is None \
                or args.worker_seed is None:
            print("ERROR: --worker requires --worker-endpoint / "
                  "--worker-replicate / --worker-seed", file=sys.stderr)
            return 2
        rbfe = _load_rbfe()
        driver = _load_driver()
        t0 = time.time()
        try:
            run_one_replicate(
                rbfe, driver,
                out_root=out_root, endpoint=args.worker_endpoint, leg=leg,
                replicate_index=args.worker_replicate, seed=args.worker_seed,
                directions=directions,
                n_windows_half=args.n_windows_half,
                softcore_band=args.softcore_band,
                n_apex_bridge=args.n_apex_bridge, apex_band=args.apex_band,
                n_cycles=args.n_cycles,
                md_steps_per_cycle=args.md_steps_per_cycle,
                platform_name=args.platform, timestep_fs=args.timestep_fs,
                minimize_iters=args.minimize_iters,
                backward_equil_steps=args.backward_equil_steps,
                genuine_decouple_nm=args.genuine_decouple_nm,
                mtr_ncaa_xml=args.mtr_ncaa_xml, binder_chain=args.binder_chain,
                archive_existing=not args.no_archive_existing,
                reuse_serialized=args.reuse_serialized,
                construction=construction, displacement_nm=displacement_nm,
                lambda1_rampdown=lambda1_rampdown)
        except Exception as exc:  # noqa: BLE001 — surface as non-zero worker exit
            print("WORKER FAILED %s/%s rep%d: %s"
                  % (args.worker_endpoint, leg, args.worker_replicate, exc),
                  file=sys.stderr)
            return 1
        print("WORKER OK %s/%s rep%d seed=%s in %.1f s"
              % (args.worker_endpoint, leg, args.worker_replicate,
                 args.worker_seed, time.time() - t0))
        return 0

    config = {
        "leg": leg, "endpoints": endpoints, "seeds": seeds,
        "directions": directions, "n_windows_half": args.n_windows_half,
        "softcore_band": args.softcore_band, "n_apex_bridge": args.n_apex_bridge,
        "n_cycles": args.n_cycles, "md_steps_per_cycle": args.md_steps_per_cycle,
        "platform": args.platform, "timestep_fs": args.timestep_fs,
        "genuine_decouple_nm": args.genuine_decouple_nm,
        "construction": construction, "displacement_nm": displacement_nm,
        "lambda1_rampdown": lambda1_rampdown,
        "out_root": out_root, "mintimeid": mintimeid,
    }

    # C11 pre-registration (BEFORE the run / analysis).
    prereg_path = write_pre_registration(out_root, config)
    print("C11 pre-registration -> %s" % (prereg_path,))

    if args.dry_run:
        plan = {
            "leg": leg,
            "construction": construction,
            "endpoints": endpoints,
            "seeds (n=%d)" % len(seeds): seeds,
            "directions": directions,
            "ladder": "%d windows/half, softcore_band=%d, apex_bridge=%d"
                      % (args.n_windows_half, args.softcore_band,
                         args.n_apex_bridge),
            "cycles": args.n_cycles,
            "runs": len(endpoints) * len(seeds) * len(directions),
            "out_root": out_root,
        }
        if construction == "twocopy":
            # Surface the two-copy ATS schedule shape (n_windows_half=6 => 11 λ/leg;
            # with --lambda1-rampdown the leg-down is densified => 12 λ/leg).
            rbfe_mod = _load_rbfe()
            sch = rbfe_mod.build_ats_standard_ladder(
                n_windows_half=args.n_windows_half, single_direction="forward",
                lambda1_rampdown=lambda1_rampdown)
            plan["schedule_kind"] = "ats_standard"
            plan["n_lambda_per_leg"] = sch["n_states"]
            plan["lambdas_1"] = sch["lambdas_1"]
            plan["lambdas_2"] = sch["lambdas_2"]
            plan["lambda1_rampdown"] = lambda1_rampdown
            plan["u0_kcal"] = sch["u0"][0]
            plan["displacement_nm"] = (
                displacement_nm if displacement_nm is not None
                else "default (ATS_TWOCOPY_DISPLACEMENT_NM)")
        print(json.dumps({"plan": plan, "config": config}, indent=2,
                         default=str))
        return 0

    if not args.analyze_only:
        # Fail-loud seed-availability gate: every matched (endpoint, seed) must
        # have its endpoint final PDB (a missing one breaks the paired ddG; R-18).
        seed_gate = gate_seed_availability(endpoints, seeds)
        if not seed_gate["passed"]:
            print("SEED-AVAILABILITY GATE FAILED — missing endpoint PDB(s):",
                  file=sys.stderr)
            for m in seed_gate["missing"]:
                print("  " + m, file=sys.stderr)
            return 2

        ladder_args = {
            "n_windows_half": args.n_windows_half,
            "softcore_band": args.softcore_band,
            "n_apex_bridge": args.n_apex_bridge, "apex_band": args.apex_band,
            "n_cycles": args.n_cycles,
            "md_steps_per_cycle": args.md_steps_per_cycle,
            "platform": args.platform, "timestep_fs": args.timestep_fs,
            "minimize_iters": args.minimize_iters,
            "backward_equil_steps": args.backward_equil_steps,
            "genuine_decouple_nm": args.genuine_decouple_nm,
            "mtr_ncaa_xml": args.mtr_ncaa_xml, "binder_chain": args.binder_chain,
            "construction": construction, "displacement_nm": displacement_nm,
            "lambda1_rampdown": lambda1_rampdown,
        }

        if args.pool:
            # VRAM-aware concurrent pool on one local GPU (no idle GPU).
            free_gb = _query_local_gpu_free_gb(args.device_index)
            unit_gb = _BOUND_UNIT_VRAM_GB if leg == "bound" else _FREE_UNIT_VRAM_GB
            if args.max_concurrent > 0:
                max_conc = args.max_concurrent
            elif free_gb is not None:
                max_conc = _max_concurrent_units(free_gb, unit_gb)
            else:
                max_conc = 1
            n_units = len(endpoints) * len(seeds)
            max_conc = min(max_conc, n_units)
            print("[pool/%s] GPU %d free=%s GiB; per-unit~%.1f GiB; "
                  "max_concurrent=%d over %d units"
                  % (leg, args.device_index,
                     ("%.1f" % free_gb) if free_gb is not None else "?",
                     unit_gb, max_conc, n_units))
            t0 = time.time()
            pool_result = run_pool_local(
                out_root=out_root, leg=leg, endpoints=endpoints, seeds=seeds,
                directions=directions, device_index=args.device_index,
                max_concurrent=max_conc, ladder_args=ladder_args,
                archive_existing=not args.no_archive_existing,
                reuse_serialized=args.reuse_serialized)
            print("[pool/%s] drained in %.1f s; all_ok=%s"
                  % (leg, time.time() - t0, pool_result["all_ok"]))
            with open(os.path.join(out_root, "pool_manifest_%s.json" % (leg,)),
                      "w") as fh:
                json.dump(pool_result, fh, indent=2, default=str)
        else:
            for endpoint in endpoints:
                t0 = time.time()
                print("[%s/%s] launching %d replicate(s) x %d direction(s)..."
                      % (endpoint, leg, len(seeds), len(directions)))
                leg_result = run_leg(
                    out_root=out_root, endpoint=endpoint, leg=leg,
                    seeds=seeds, directions=directions,
                    n_windows_half=args.n_windows_half,
                    softcore_band=args.softcore_band,
                    n_apex_bridge=args.n_apex_bridge, apex_band=args.apex_band,
                    n_cycles=args.n_cycles,
                    md_steps_per_cycle=args.md_steps_per_cycle,
                    platform_name=args.platform, timestep_fs=args.timestep_fs,
                    minimize_iters=args.minimize_iters,
                    backward_equil_steps=args.backward_equil_steps,
                    genuine_decouple_nm=args.genuine_decouple_nm,
                    mtr_ncaa_xml=args.mtr_ncaa_xml,
                    binder_chain=args.binder_chain,
                    archive_existing=not args.no_archive_existing,
                    reuse_serialized=args.reuse_serialized,
                    construction=construction, displacement_nm=displacement_nm,
                    lambda1_rampdown=lambda1_rampdown)
                print("[%s/%s] done in %.1f s (%d replicates)"
                      % (endpoint, leg, time.time() - t0,
                         leg_result["n_replicates"]))

    # UWHAM analysis — only when BOTH cp4 + wt + both directions ran (the cycle).
    if set(endpoints) == set(ENDPOINTS) and set(directions) == set(DIRECTION_TAGS):
        analysis_json = os.path.join(out_root, "ddint_%s.json" % (leg,))
        print("[%s] UWHAM analysis -> %s" % (leg, analysis_json))
        result = analyze_leg_cohort(
            out_root, leg, n_replicates=len(seeds), mintimeid=mintimeid,
            output_json=analysis_json)
        ddint = result.get("ddint_point_estimate_kcal")
        quad = result.get("quadrature_C7") or {}
        print("ddint_%s point = %.4f kcal/mol; quadrature SE = %s; z_SE = %s; "
              "sign = %s (%s)"
              % (leg, ddint, quad.get("se_quadrature"), quad.get("z_se"),
                 result.get("sign_status"), result.get("sign_reason")))
    else:
        print("(skipping UWHAM: need cp4+wt + both directions for the cycle)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
