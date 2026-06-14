#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B production v2.1 — upstream make_system + async_re ABFE launcher.

v2 (commit a9442f6) was architecturally invalid because
``utils/atm_trackB_setup.build_leg_system`` placed the binder at the
binding site BEFORE Modeller.addSolvent — INVERTED from upstream's
"binder pre-displaced into bulk, then solvate" pattern. Result:
bound-leg first-step NaN at lambda=0.5 (r11/r12 replicas, evidence at
``_archive/v2_custom_build_partial_20260531/cp4/bound/_production.log``).

v2.1 replaces the custom system builder with
``scripts/phase4_trackB_v2_make_system.build_all_four_systems``, which
wraps the upstream ``atom_openmm.make_atm_system_from_rcpt_lig.make_system``
gold-standard ABFE builder. The cntl / nodefile / structprep /
production glue is reused verbatim from v2 (preserved valuable v2 work
from commit a9442f6).

# Per-leg directory layout (mirrors v2 + upstream contract)

    outputs/_trackb/production_v2_1/<endpoint>/<leg>/
        trackb.pdb            (built by phase4_trackB_v2_make_system)
        trackb_sys.xml        (built by phase4_trackB_v2_make_system)
        trackb_asyncre.cntl   (written by this launcher)
        nodefile              (written by this launcher)
        _structprep.log       (abfe_structprep stdout/stderr)
        _production.log       (abfe_production stdout/stderr)
        r0/ .. r21/           (per-replica walker output)
        trackb_stat.txt       (live replica status)

# Thermodynamic cycle (unchanged)

    DeltaG_bind(endpoint) = DG_alch(bound) - DG_alch(free)
    DeltaDeltaG_bind_Cp4_minus_WT = DG_bind(Cp4) - DG_bind(WT)

# Charge axis (unchanged)

Hybrid MTR XML (``params/MTR_gaff2_hybrid.xml``), per-residue |Sigma q|
<= 5e-4 e, NE1 frozen at -0.3418, cyclic_ss disulfide preserved on both
legs (commit per binder protonation).

# Hardware contract

* host 5070Ti (single GPU, CUDA device 0) for the smoke and initial
  production. Track A V100 untouched until Track A QM completion.
* 22 replica walkers on 1 GPU -> throughput-limited by GPU compute.

# Cross-references

* v2 archive (custom-build NaN): ``outputs/_trackb/_archive/v2_custom_build_partial_20260531/``
* v1 archive (direction-flip NaN): ``outputs/_trackb/_archive/v1_architectural_failed_20260531/``
* System builder: ``scripts/phase4_trackB_v2_make_system.py``
* v2 launcher (deprecated, kept for regression): ``scripts/trackb_production_v2_asyncre.py``
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime
from typing import Optional, List, Dict, Tuple, Any

_PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "utils"))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "scripts"))

# Re-use v2's cntl + nodefile + binary-locator + production helpers
# (preserves the work from v2 commit a9442f6).
import trackb_production_v2_asyncre as v2_legacy  # noqa: E402

# Schedule constants + helpers (re-exported from v2 for test regression)
LAMBDA_FWD = v2_legacy.LAMBDA_FWD
LAMBDA_BWD = v2_legacy.LAMBDA_BWD
DIRECTIONS = v2_legacy.DIRECTIONS
LAMBDAS_1 = v2_legacy.LAMBDAS_1
LAMBDAS_2 = v2_legacy.LAMBDAS_2
INTERMD = v2_legacy.INTERMD
W0 = v2_legacy.W0
ALPHA = v2_legacy.ALPHA
U0 = v2_legacy.U0
TEMP_K = v2_legacy.TEMP_K
N_STATES = v2_legacy.N_STATES
UMAX_KCAL = v2_legacy.UMAX_KCAL
UBCORE_KCAL = v2_legacy.UBCORE_KCAL
ACORE = v2_legacy.ACORE
DISPLACEMENT_NM_DEFAULT = v2_legacy.DISPLACEMENT_NM_DEFAULT
PRE_REGISTER_OUTCOMES = v2_legacy.PRE_REGISTER_OUTCOMES

# Helpers (reused verbatim from v2)
write_cntl_file = v2_legacy.write_cntl_file
write_nodefile = v2_legacy.write_nodefile
collect_binder_atom_indices = v2_legacy.collect_binder_atom_indices
collect_receptor_ca_indices = v2_legacy.collect_receptor_ca_indices
run_abfe_structprep_for_leg = v2_legacy.run_abfe_structprep_for_leg
run_abfe_production_for_leg = v2_legacy.run_abfe_production_for_leg
_find_atm_binary = v2_legacy._find_atm_binary
_atom_openmm_version = v2_legacy._atom_openmm_version


# ---------------------------------------------------------------------------
# Structprep cache integrity (cross-stage drift enforcement, 2026-05-31)
# ---------------------------------------------------------------------------
# Background: the earlier cache check at
# ``scripts/trackb_production_v2_1_upstream.py:594-601``
# silently reused a SMOKE-mode
# structprep artifact (cp4/free 6000 steps, T=271 K) as a production
# starting state, causing cycle-1 NaN at the BWD intermediate replica.
#
# Web evidence (Snakemake / Nextflow / cwltool): mature pipeline managers
# carry per-output metadata sidecars (mode / params / code hash) and
# fail-fast on cache invalidation rather than silently reuse stale state.
# We adopt the same pattern for Track B structprep artifacts.
#
# Schema (``<basename>_0.xml.meta.json``):
#   {
#     "schema_version":            "0.9.13",
#     "mode":                       "smoke" | "production",
#     "step_count":                 int,        # last step recorded in log
#     "final_temp_K":               float,      # last T column of log
#     "structprep_completed_at":    ISO-8601 str,
#     "trackb_0_xml_sha256":        str,        # SHA-256 of the cached state
#     "ckpt_is_valid_present":      bool,
#   }
#
# Validation contract:
#   1. Prefer sidecar (deterministic, includes mode flag).
#   2. Fallback to ``_structprep.log`` parsing (legacy artifacts predating
#      this sidecar emission contract).
#   3. Neither present → IntegrityError (fail-fast).
#
# Auto-deletion is DELIBERATELY rejected (Snakemake ``--cleanup-metadata``
# style operator intervention pattern). The launcher prints the archive
# path and remediation steps; an operator (or a dedicated rerun
# script) must archive the stale artifact before the cohort proceeds.
# ---------------------------------------------------------------------------


SIDECAR_SCHEMA_VERSION = "0.9.13"

# Production-quality thresholds. Used by ``_validate_cached_structprep``
# when the caller does not override.
#
# 100000 steps is the conservative lower bound:
#   - production default = 850000 (550k thermalization + 100k NPT +
#     100k NVT + 100k annealing + 50k equilibration; varies slightly by
#     atom_openmm version).
#   - smoke default = ~6000 (1000 therm + 2000 anneal + 2000 equil +
#     residual).
#   - 100000 sits firmly between the two, with no expected
#     legitimate-production run dipping below it.
PROD_MIN_STRUCTPREP_STEPS = 100000
PROD_TARGET_TEMP_K = 300.0
PROD_TEMP_TOLERANCE_K = 10.0


class IntegrityError(RuntimeError):
    """Cross-stage integrity violation.

    Raised when the declared cache state (``trackb_0.xml`` present →
    ``status="cached"``) does not match the actual artifact metadata
    (e.g., smoke-mode structprep silently reused under a production
    launch).

    The launcher does NOT auto-archive on this error. The operator
    is responsible for archiving the stale artifact + relaunch.
    """
    pass


def _sha256_file(path: str) -> str:
    """SHA-256 hex digest of a file, 1 MiB streaming reads.

    Used to record the cache state hash in the sidecar so subsequent
    launches can detect silent in-place mutation of ``trackb_0.xml``.
    """
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


# atom_openmm step lines have the form "<step>,<PE>,<T>,<box?>,<speed>" —
# always integer step followed by comma + numeric. Match the LAST such
# line in the log (the FINAL recorded step).
_STEP_LINE_RE = re.compile(r"^(\d+),(-?\d+\.?\d*[eE]?-?\d*),(-?\d+\.?\d*[eE]?-?\d*),")


def _parse_structprep_log(log_path: str) -> Tuple[int, float]:
    """Parse ``_structprep.log`` to extract (last_step, final_temp_K).

    Returns (0, 0.0) if no step lines are found (e.g., crashed before any
    integration step). Caller must treat (0, 0.0) as a validation failure.

    The atom_openmm structprep log writes a step line every
    ``prnt_frequency`` (default 100) integration steps in the format
    ``"<step>,<PE>,<T>,<box?>,<speed>"`` (NPT phase includes box,
    NVT/anneal/equil omit it but still have T as the 3rd column).
    """
    last_step = 0
    final_T = 0.0
    if not os.path.isfile(log_path):
        return (0, 0.0)
    with open(log_path, "r", errors="replace") as fh:
        for line in fh:
            m = _STEP_LINE_RE.match(line)
            if m is None:
                continue
            try:
                step = int(m.group(1))
                T = float(m.group(3))
            except (ValueError, IndexError):
                continue
            if step >= last_step:
                last_step = step
                final_T = T
    return (last_step, final_T)


def _emit_structprep_sidecar(
    leg_dir: str,
    mode: str,
    basename: str = "trackb",
    structprep_log_path: Optional[str] = None,
) -> str:
    """Emit ``<basename>_0.xml.meta.json`` after structprep completes.

    Records the protocol mode (smoke / production), the observed final
    step count + temperature, the state-XML SHA-256, and the
    ``ckpt_is_valid`` marker presence. Future launches read this sidecar
    via ``_validate_cached_structprep`` to detect silent mode mismatches.

    Returns the absolute path of the emitted sidecar.

    Raises:
        ValueError: ``mode`` is not "smoke" or "production".
        FileNotFoundError: ``<basename>_0.xml`` is missing in ``leg_dir``.
    """
    if mode not in ("smoke", "production"):
        raise ValueError(f"sidecar mode must be 'smoke' or 'production', got {mode!r}")

    state_xml = os.path.join(leg_dir, f"{basename}_0.xml")
    if not os.path.isfile(state_xml):
        raise FileNotFoundError(
            f"cannot emit sidecar — state XML missing: {state_xml}"
        )

    if structprep_log_path is None:
        structprep_log_path = os.path.join(leg_dir, "_structprep.log")
    last_step, final_T = _parse_structprep_log(structprep_log_path)

    sidecar_path = os.path.join(leg_dir, f"{basename}_0.xml.meta.json")
    meta = {
        "schema_version": SIDECAR_SCHEMA_VERSION,
        "mode": mode,
        "step_count": int(last_step),
        "final_temp_K": float(final_T),
        "structprep_completed_at": datetime.now().isoformat(timespec="seconds"),
        "trackb_0_xml_sha256": _sha256_file(state_xml),
        "ckpt_is_valid_present": os.path.isfile(
            os.path.join(leg_dir, "ckpt_is_valid")
        ),
    }
    with open(sidecar_path, "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")
    return sidecar_path


def _validate_cached_structprep(
    leg_dir: str,
    expected_mode: str = "production",
    expected_min_steps: int = PROD_MIN_STRUCTPREP_STEPS,
    expected_temp_K: float = PROD_TARGET_TEMP_K,
    temp_tolerance_K: float = PROD_TEMP_TOLERANCE_K,
    basename: str = "trackb",
) -> Tuple[bool, str]:
    """Validate that a cached structprep artifact matches the expected
    protocol BEFORE marking it ``"cached"``.

    Returns:
        (True, validation summary) if the artifact is acceptable.
        (False, actionable failure reason) otherwise.

    Validation order:
        1. Sidecar ``<basename>_0.xml.meta.json`` exists → preferred.
        2. Fallback: parse ``_structprep.log`` last_step + final_T.
        3. Neither present → fail.

    Failure reasons explicitly include the leg directory + remediation
    step (archive path suggestion) so the operator can act without
    re-reading code.
    """
    sidecar_path = os.path.join(leg_dir, f"{basename}_0.xml.meta.json")

    if os.path.isfile(sidecar_path):
        try:
            with open(sidecar_path, "r") as fh:
                meta = json.load(fh)
        except (json.JSONDecodeError, OSError) as exc:
            return False, (
                f"Cached structprep sidecar unreadable at {sidecar_path}: {exc}. "
                f"Archive {leg_dir} contents to "
                f"outputs/_trackb/_archive/<basename>_<timestamp>/ "
                f"(preserve, do not delete) and rerun this launcher."
            )

        observed_mode = meta.get("mode")
        if observed_mode != expected_mode:
            return False, (
                f"Cached structprep mode mismatch in {leg_dir}: sidecar "
                f"says mode={observed_mode!r}, current launch expects "
                f"{expected_mode!r}. Archive {leg_dir} contents to "
                f"outputs/_trackb/_archive/<basename>_<timestamp>/ "
                f"(preserve, do not delete) and rerun this launcher to regenerate "
                f"at {expected_mode} quality."
            )

        observed_steps = int(meta.get("step_count", 0))
        if observed_steps < expected_min_steps:
            return False, (
                f"Cached structprep step_count={observed_steps} < "
                f"{expected_min_steps} (smoke-mode artifact) in {leg_dir}. "
                f"Archive {leg_dir} contents to "
                f"outputs/_trackb/_archive/<basename>_<timestamp>/ "
                f"(preserve, do not delete) and rerun this launcher to regenerate "
                f"at production quality."
            )

        observed_T = float(meta.get("final_temp_K", 0.0))
        if abs(observed_T - expected_temp_K) > temp_tolerance_K:
            return False, (
                f"Cached structprep final_T={observed_T:.2f} K differs from "
                f"target {expected_temp_K:.1f} +/- {temp_tolerance_K:.1f} K "
                f"in {leg_dir} (cold / under-equilibrated artifact). "
                f"Archive {leg_dir} contents to "
                f"outputs/_trackb/_archive/<basename>_<timestamp>/ "
                f"(preserve, do not delete) and rerun this launcher."
            )

        return True, (
            f"Cached structprep validated via sidecar "
            f"(mode={observed_mode}, steps={observed_steps}, "
            f"T={observed_T:.2f} K)"
        )

    # Fallback: parse _structprep.log (legacy artifacts pre-dating sidecar).
    log_path = os.path.join(leg_dir, "_structprep.log")
    if not os.path.isfile(log_path):
        return False, (
            f"No sidecar AND no _structprep.log in {leg_dir} — cannot "
            f"validate cache freshness. Archive {leg_dir} contents to "
            f"outputs/_trackb/_archive/<basename>_<timestamp>/ "
            f"(preserve, do not delete) and rerun this launcher."
        )

    last_step, final_T = _parse_structprep_log(log_path)
    if last_step == 0:
        return False, (
            f"Legacy _structprep.log in {leg_dir} contains no parseable "
            f"step lines — cannot validate cache freshness. Archive "
            f"{leg_dir} contents to "
            f"outputs/_trackb/_archive/<basename>_<timestamp>/ "
            f"(preserve, do not delete) and rerun this launcher."
        )
    if last_step < expected_min_steps:
        return False, (
            f"Legacy _structprep.log last_step={last_step} < "
            f"{expected_min_steps} in {leg_dir} (smoke-mode artifact, "
            f"no sidecar). Archive {leg_dir} contents to "
            f"outputs/_trackb/_archive/<basename>_<timestamp>/ "
            f"(preserve, do not delete) and rerun. Future runs will auto-emit "
            f"a sidecar so this fallback path is bypassed."
        )
    if abs(final_T - expected_temp_K) > temp_tolerance_K:
        return False, (
            f"Legacy _structprep.log final_T={final_T:.2f} K differs from "
            f"target {expected_temp_K:.1f} +/- {temp_tolerance_K:.1f} K "
            f"in {leg_dir} (cold / under-equilibrated). Archive {leg_dir} "
            f"contents to "
            f"outputs/_trackb/_archive/<basename>_<timestamp>/ "
            f"(preserve, do not delete) and rerun this launcher."
        )

    return True, (
        f"Legacy cache validated (no sidecar; _structprep.log "
        f"last_step={last_step}, final_T={final_T:.2f} K)"
    )


# ---------------------------------------------------------------------------
# Per-leg orchestration (v2.1 build-path-only differs from v2)
# ---------------------------------------------------------------------------
def setup_one_leg_v21(
    endpoint: str,
    leg: str,
    out_root: str,
    build_meta: Dict[str, Any],
    jobname: str = "trackb",
    binder_chain: str = "B",
    receptor_chain: str = "A",
    displacement_nm: Tuple[float, float, float] = DISPLACEMENT_NM_DEFAULT,
    cuda_devices: Optional[List[int]] = None,
    production_steps: int = 2500,
    prnt_frequency: int = 2500,
    trj_frequency: int = 25000,
    max_samples: int = 1000,
    wall_time_min: int = 720,
    cycle_time_s: int = 10,
    checkpoint_time_s: int = 600,
    smoke: bool = False,
    free_schedule: str = v2_legacy.DEFAULT_FREE_SCHEDULE,
    bound_schedule: str = v2_legacy.DEFAULT_BOUND_SCHEDULE,
) -> Dict[str, Any]:
    """Prepare one (endpoint, leg) for abfe_production. System XML pre-built.

    The upstream-built system XML + topology PDB already live in
    ``out_root/<endpoint>/<leg>/`` (named ``trackb_sys.xml`` and
    ``trackb.pdb``). This function only writes the cntl + nodefile.

    ``free_schedule`` (default ``"canonical22"``) selects the λ ladder for the
    FREE leg only (λ-densify spec, 2026-06-05). ``"densified34"`` swaps the
    free-leg cntl to the 34-state ilogistic anneal ladder.

    ``bound_schedule`` (default ``"canonical22"``) selects the λ ladder for the
    BOUND leg only (densified_bound28 spec 2026-06-11). ``"densified_bound28"``
    writes the combined bound-leg cntl as the 28-state (dplus 14 + dminus 14)
    soft-core-bridge ladder that closes the 6→7 cliff (λ0.30→0.35) BC=0
    zero-overlap; ``"densified_bound30"`` is the 4-window cap-1 escalation of the
    same cliff (30-state = dplus 15 + dminus 15, denser λ across the same span —
    window COUNT is the lever, U0 still DESCENDS 105→97). This combined cntl is
    the schedule SSOT that ``trackb_per_direction_structprep`` +
    ``generate_per_direction_cntls`` then split (count-agnostically) into the
    per-direction cntls (28→14/14, 30→15/15). The free/bound selectors are
    INDEPENDENT (a bound densify does not perturb the free leg and vice-versa).
    The valid bound schedules are ``v2_legacy.BOUND_SCHEDULES`` (the free-leg
    densified38* ladders are NOT bound-valid — they fix the free-leg decoupling
    crossover the receptor-held bound leg does not have).
    """
    leg_dir = os.path.join(out_root, endpoint, leg)
    if not os.path.isdir(leg_dir):
        raise RuntimeError(
            f"Leg dir missing — run phase4_trackB_v2_make_system first: {leg_dir}"
        )
    topology_pdb = os.path.join(leg_dir, jobname + ".pdb")
    system_xml = os.path.join(leg_dir, jobname + "_sys.xml")
    if not os.path.isfile(topology_pdb):
        raise RuntimeError(f"Missing topology PDB: {topology_pdb}")
    if not os.path.isfile(system_xml):
        raise RuntimeError(f"Missing system XML:  {system_xml}")

    # Smoke override: clamp production cycles + structprep step counts.
    thermalization_steps = None
    annealing_steps = None
    equilibration_steps = None
    steps_per_cycle = None
    if smoke:
        production_steps = 50
        prnt_frequency = 50
        trj_frequency = 50
        max_samples = 2
        wall_time_min = 15
        checkpoint_time_s = 60
        thermalization_steps = 1000
        annealing_steps = 2000
        equilibration_steps = 1000
        steps_per_cycle = 100

    # Upstream-built PDB has chain "L" for the binder (upstream renames at
    # make_atm_system_from_rcpt_lig:241-243). Our v2 helpers expect chain
    # "B"; we collect by the upstream-renamed chain id.
    upstream_binder_chain = "L" if leg == "bound" else binder_chain
    ligand_atoms = collect_binder_atom_indices(
        topology_pdb, binder_chain=upstream_binder_chain
    )
    if not ligand_atoms and upstream_binder_chain != binder_chain:
        # Fallback: chain rename may not have happened on PDB load path.
        ligand_atoms = collect_binder_atom_indices(
            topology_pdb, binder_chain=binder_chain
        )
    if not ligand_atoms:
        raise RuntimeError(
            f"No binder-chain atoms found in {topology_pdb}"
        )

    if leg == "bound":
        # Upstream preserves the receptor chain id from the input PDB; our
        # receptor PDB has chain A.
        pos_restrained = collect_receptor_ca_indices(
            topology_pdb, receptor_chain=receptor_chain
        )
        if not pos_restrained:
            raise RuntimeError(
                f"No receptor Calpha atoms found in {topology_pdb} "
                f"(chain {receptor_chain})"
            )
    else:
        pos_restrained = []

    # Resolve the λ schedule. The free + bound selectors are INDEPENDENT:
    # ``free_schedule`` (densified38*) applies to the FREE leg only; the
    # densified38* ladders fix the free-leg decoupling crossover. The BOUND
    # leg uses ``bound_schedule`` (densified_bound28) which closes the bound
    # 6→7 cliff (λ0.30→0.35) BC=0 zero-overlap. Default canonical22 on both
    # legs → schedule_dict=None → byte-equal to the historical cntl. The
    # resulting combined bound cntl is the schedule SSOT that
    # generate_per_direction_cntls splits (count-agnostically) into dplus/dminus.
    schedule_dict: Optional[Dict[str, Any]] = None
    schedule_name = v2_legacy.DEFAULT_FREE_SCHEDULE
    if leg == "free" and free_schedule != v2_legacy.DEFAULT_FREE_SCHEDULE:
        schedule_dict = v2_legacy.get_schedule(free_schedule)
        schedule_name = free_schedule
    elif leg == "bound" and bound_schedule != v2_legacy.DEFAULT_BOUND_SCHEDULE:
        if bound_schedule not in v2_legacy.BOUND_SCHEDULES:
            raise ValueError(
                f"bound_schedule {bound_schedule!r} is not valid for the bound "
                f"leg; valid bound schedules: "
                f"{sorted(v2_legacy.BOUND_SCHEDULES)}. The free-leg densified38* "
                f"ladders fix the free-leg decoupling crossover, which the "
                f"receptor-held bound leg does not have."
            )
        schedule_dict = v2_legacy.get_schedule(bound_schedule)
        schedule_name = bound_schedule

    # Write nodefile + cntl
    nodefile_path = os.path.join(leg_dir, "nodefile")
    write_nodefile(nodefile_path,
                   gpu_indices=(cuda_devices or [0]),
                   platform="CUDA")
    cntl_path = os.path.join(leg_dir, jobname + "_asyncre.cntl")
    write_cntl_file(
        cntl_path=cntl_path,
        basename=jobname,
        nodefile_path=nodefile_path,
        ligand_atom_indices=ligand_atoms,
        pos_restrained_atom_indices=pos_restrained,
        displacement_nm=displacement_nm,
        production_steps=production_steps,
        prnt_frequency=prnt_frequency,
        trj_frequency=trj_frequency,
        max_samples=max_samples,
        wall_time_min=wall_time_min,
        cycle_time_s=cycle_time_s,
        checkpoint_time_s=checkpoint_time_s,
        thermalization_steps=thermalization_steps,
        annealing_steps=annealing_steps,
        equilibration_steps=equilibration_steps,
        steps_per_cycle=steps_per_cycle,
        schedule=schedule_dict,
    )

    return {
        "endpoint": endpoint,
        "leg": leg,
        "leg_dir": leg_dir,
        "cntl_path": cntl_path,
        "nodefile_path": nodefile_path,
        "topology_pdb": topology_pdb,
        "system_xml": system_xml,
        "n_ligand_atoms": len(ligand_atoms),
        "n_pos_restrained": len(pos_restrained),
        "build_meta": build_meta,
        "smoke_mode": smoke,
        "schedule_name": schedule_name,
        "n_states": (schedule_dict["n_states"] if schedule_dict is not None
                     else N_STATES),
    }


# ---------------------------------------------------------------------------
# Smoke gate (bound leg, both directions, lambda=0.5 mid-point)
# ---------------------------------------------------------------------------
def bound_leg_bidirectional_smoke(
    leg_dir: str,
    n_steps: int = 30,
    jobname: str = "trackb",
) -> Dict[str, Any]:
    """Bidirectional NaN gate at lambda=0.5 via OMMSystemABFE.

    Rebuilds the production-equivalent ATMForce-wrapped system from the
    cntl file using upstream ``atom_openmm.ommsystem.OMMSystemABFE``
    (same code path abfe_production uses), then runs ``n_steps`` MD
    steps at lambda1=lambda2=0.5 with direction=+1 and again with
    direction=-1 — both in independent Contexts (no in-place flip;
    matches the v2.1 architectural fix). Acceptance: no NaN, |u0|/|u1|
    finite, max force < 1e4 kcal/mol/A, SG-SG std < 0.1 A (cyclic_ss
    intact).

    The post-structprep equilibrated State (``trackb_0.xml``) provides
    the starting positions/velocities/box. The bare System
    (``trackb_sys.xml``) is re-wrapped via ``OMMSystemABFE.create_system``
    to ensure ATMForce + restraints match production.
    """
    import openmm as mm
    import openmm.unit as unit
    from openmm import app
    import numpy as np
    import math
    import logging

    state_xml = os.path.join(leg_dir, jobname + "_0.xml")
    if not os.path.isfile(state_xml):
        raise RuntimeError(
            f"No post-structprep state XML at {state_xml}; "
            f"run abfe_structprep first."
        )
    sys_xml = os.path.join(leg_dir, jobname + "_sys.xml")
    topology_pdb = os.path.join(leg_dir, jobname + ".pdb")
    cntl_path = os.path.join(leg_dir, jobname + "_asyncre.cntl")
    if not (os.path.isfile(sys_xml) and os.path.isfile(topology_pdb)
            and os.path.isfile(cntl_path)):
        raise RuntimeError(
            f"Missing system XML / topology PDB / cntl in {leg_dir}"
        )

    # Rebuild via OMMSystemABFE using upstream's parse_config (which
    # int-casts the LIGAND_CM_ATOMS / RCPT_CM_ATOMS / POS_RESTRAINED_ATOMS
    # lists so CustomCentroidBondForce.addGroup gets list[int] not
    # list[str]). Then run massage_keywords with restrain_solutes=False to
    # set TIME_STEP without overriding POS_RESTRAINED_ATOMS.
    from atom_openmm.ommsystem import OMMSystemABFE
    from atom_openmm.abfe_structprep import massage_keywords
    from atom_openmm.utils.config import parse_config
    keywords = parse_config(cntl_path)
    massage_keywords(keywords, restrain_solutes=False)
    logger = logging.getLogger("trackb_v21_smoke")
    syswrap = OMMSystemABFE(jobname, keywords, topology_pdb, sys_xml, logger)
    syswrap.create_system()
    wrapped_system = syswrap.system

    # Read PDB topology (for atom enumeration; positions come from State).
    pdb = app.PDBFile(topology_pdb)

    # CYS SG indices for cyclic_ss audit
    sg_idx = []
    for chain in pdb.topology.chains():
        if chain.id not in ("B", "L"):
            continue
        for res in chain.residues():
            if res.name in ("CYS", "CYX"):
                for atom in res.atoms():
                    if atom.name == "SG":
                        sg_idx.append(atom.index)

    out: Dict[str, Any] = {
        "leg_dir": leg_dir,
        "state_xml": state_xml,
        "system_xml": sys_xml,
        "n_steps": n_steps,
        "directions": {},
    }

    for direction_val in (1.0, -1.0):
        # Fresh wrapper per direction (independent walker — matches the
        # v2.1 architectural fix: no in-place direction flip). Re-builds
        # the System AND the production-canonical ATMMTSLangevinIntegrator
        # (multi-timestep, ATMForce in slow group) — using a single-step
        # LangevinMiddleIntegrator here would explode at d=-1 because the
        # soft-core gradient is too stiff for 2 fs without MTS.
        from atom_openmm.ommsystem import OMMSystemABFE
        syswrap_dir = OMMSystemABFE(
            jobname, keywords, topology_pdb, sys_xml, logger
        )
        syswrap_dir.create_system()
        wrapped_system_dir = syswrap_dir.system
        integrator = syswrap_dir.integrator
        platform = mm.Platform.getPlatformByName("CUDA")
        sim = app.Simulation(
            pdb.topology, wrapped_system_dir, integrator, platform
        )
        # Load the equilibrated state (positions, velocities, box).
        with open(state_xml) as fh:
            sim.context.setState(mm.XmlSerializer.deserialize(fh.read()))
        # Set the alchemical mid-point and direction.
        sim.context.setParameter("Lambda1", 0.5)
        sim.context.setParameter("Lambda2", 0.5)
        sim.context.setParameter("Direction", direction_val)
        # Reference the per-direction wrapped system for downstream
        # ATMForce introspection (replaces the outer wrapped_system).
        wrapped_system = wrapped_system_dir

        nans = 0
        crashed_at = None
        max_force = 0.0
        sg_pairs: List[float] = []
        for step in range(n_steps):
            try:
                sim.step(1)
            except Exception as e:
                crashed_at = step
                nans = n_steps - step
                out["directions"][f"d={int(direction_val):+d}"] = {
                    "crashed_at_step": step,
                    "error": repr(e)[:300],
                    "completed_steps": step,
                    "n_nans": nans,
                }
                break
            st = sim.context.getState(
                getEnergy=True, getForces=True, getPositions=True
            )
            pe_kj = st.getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole
            )
            if not math.isfinite(pe_kj):
                nans += 1
                break
            forces = st.getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer
            )
            max_force = max(max_force, float(np.max(np.abs(forces))))
            if len(sg_idx) >= 2:
                pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                sg_pairs.append(
                    float(np.linalg.norm(pos[sg_idx[0]] - pos[sg_idx[1]]))
                )

        if crashed_at is not None:
            continue  # already recorded above

        # Final ATMForce u0/u1 perturbation
        u0_kj = u1_kj = float("nan")
        atm_forces = [wrapped_system.getForce(i)
                      for i in range(wrapped_system.getNumForces())
                      if isinstance(wrapped_system.getForce(i), mm.ATMForce)]
        if atm_forces:
            atm = atm_forces[0]
            try:
                pe = atm.getPerturbationEnergy(sim.context)
                u1_raw, u0_raw = pe[0], pe[1]
                u1_kj = (u1_raw.value_in_unit(unit.kilojoule_per_mole)
                         if hasattr(u1_raw, "value_in_unit") else float(u1_raw))
                u0_kj = (u0_raw.value_in_unit(unit.kilojoule_per_mole)
                         if hasattr(u0_raw, "value_in_unit") else float(u0_raw))
            except Exception:
                pass

        sg_std_nm = float(np.std(sg_pairs)) if sg_pairs else None
        sg_mean_nm = float(np.mean(sg_pairs)) if sg_pairs else None
        out["directions"][f"d={int(direction_val):+d}"] = {
            "completed_steps": n_steps - nans,
            "n_nans": nans,
            "u0_kj_per_mol": u0_kj,
            "u1_kj_per_mol": u1_kj,
            "max_force_kj_per_mol_nm": max_force,
            "max_force_kcal_per_mol_A": max_force / 4.184 * 0.1,
            "sg_sg_mean_nm": sg_mean_nm,
            "sg_sg_std_nm": sg_std_nm,
            "finite_u0_u1": bool(math.isfinite(u0_kj) and math.isfinite(u1_kj)),
        }

    # Acceptance gates
    ok = all(
        d.get("n_nans", 0) == 0
        and d.get("finite_u0_u1", False)
        and d.get("max_force_kcal_per_mol_A", float("inf")) < 1e4
        and (d.get("sg_sg_std_nm") is None or d["sg_sg_std_nm"] < 0.01)
        for d in out["directions"].values()
    )
    out["smoke_pass"] = ok
    out["regime"] = "ranking_only"
    return out


# ---------------------------------------------------------------------------
# Top-level main
# ---------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Track B production v2.1 — upstream make_system + async_re ABFE."
        )
    )
    p.add_argument("--seed-tag", default="s7",
                   help="prepared seed (outputs/2QKI_*_calib_<seed>)")
    p.add_argument("--out-root", default="outputs/_trackb/production_v2_1",
                   help="output root under project")
    p.add_argument("--endpoints", default="cp4,wt",
                   help="comma-list of endpoints (cp4 / wt)")
    p.add_argument("--legs", default="bound,free",
                   help="comma-list of legs (bound / free)")
    p.add_argument("--displacement-nm", default="2.5,0.0,0.0",
                   help="ABFE displacement vector (nm)")
    p.add_argument("--cuda-device", default="0",
                   help="forced CUDA_VISIBLE_DEVICES (host 5070Ti = 0)")
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--jobname", default="trackb",
                   help="basename of per-leg cntl / pdb / xml files")
    p.add_argument("--production-steps", type=int, default=2500,
                   help="MD steps per replica per cycle (default 2500)")
    p.add_argument("--max-samples", type=int, default=1000,
                   help="max cycles per replica (default 1000 ~ 5 ns)")
    p.add_argument("--wall-time-min", type=int, default=720,
                   help="abfe_production wall-time budget per leg (minutes)")
    p.add_argument("--checkpoint-time-s", type=int, default=600)
    p.add_argument("--cycle-time-s", type=int, default=10)
    p.add_argument("--smoke", action="store_true",
                   help="smoke gate: 50 steps x 2 cycles x 1 leg")
    p.add_argument("--free-schedule", default=v2_legacy.DEFAULT_FREE_SCHEDULE,
                   choices=sorted(v2_legacy.SCHEDULES),
                   help=("lambda ladder for the FREE leg only "
                         "(lambda-densify spec, 2026-06-05). 'canonical22' "
                         "(default) = 22-state symmetric; 'densified38' = "
                         "REVISED 38-state ilogistic anneal ladder (softened "
                         "backward W0-peak + per-state alpha/U0 ramp) bridging "
                         "the lambda=0.45->0.50 decoupling crossover "
                         "(PRODUCTION); 'densified34' = DEPRECATED (Factor-B "
                         "broken, forensic only). The bound leg uses "
                         "--bound-schedule (independent selector)."))
    p.add_argument("--bound-schedule",
                   default=v2_legacy.DEFAULT_BOUND_SCHEDULE,
                   choices=sorted(v2_legacy.BOUND_SCHEDULES),
                   help=("lambda ladder for the BOUND leg only "
                         "(densified_bound28 spec 2026-06-11). 'canonical22' "
                         "(default) = 22-state symmetric; 'densified_bound28' = "
                         "28-state (dplus14 + dminus14) soft-core-bridge ladder "
                         "closing the bound 6->7 cliff (lambda0.30->0.35) BC=0 "
                         "zero-overlap; 'densified_bound30' = 4-window cap-1 "
                         "escalation (30-state = dplus15 + dminus15, denser "
                         "lambda across the same cliff span). The combined bound "
                         "cntl is the schedule SSOT that per_direction_structprep "
                         "+ generate_per_direction_cntls split count-agnostically "
                         "(28->14/14, 30->15/15). Independent of --free-schedule."))
    p.add_argument("--skip-build", action="store_true",
                   help="skip system rebuild (use existing trackb.{pdb,sys.xml})")
    p.add_argument("--build-only", action="store_true",
                   help="only run phase4_trackB_v2_make_system, no structprep")
    p.add_argument("--bidir-smoke-only", action="store_true",
                   help="run only the 30-step bound bidirectional smoke gate")
    p.add_argument("--setup-only", action="store_true",
                   help="prepare cntl/system/structprep but do NOT launch production")
    p.add_argument("--skip-structprep", action="store_true",
                   help="assume <basename>_0.xml already exists from prior structprep")
    p.add_argument("--background", action="store_true",
                   help="spawn abfe_production in background (Popen) and return")
    p.add_argument("--no-charge-axis-gate", action="store_true",
                   help="skip the charge_axis verifier pre-flight (debug only)")
    p.add_argument(
        "--bound-directions",
        default="dplus,dminus",
        help=(
            "Comma-separated bound-leg directions to build via "
            "build_all_four_systems. The bound leg REQUIRES distinct "
            "per-direction systems: dplus = bound-state base (binder at "
            "binding-site coords, ~92855 particles), dminus = "
            "dissociated-state base (binder pre-displaced to bulk solvent, "
            "binding-site pocket re-hydrated by addSolvent, ~92804 "
            "particles = -17 waters). The two systems describe physically "
            "DISTINCT solvation environments and are NOT interchangeable. "
            "Default 'dplus,dminus' emits trackb_sys_dplus.xml + "
            "trackb_sys_dminus.xml so both walker directions start "
            "clash-free. Pass 'dplus' only for legacy single-system "
            "combined builds (free leg / forensic). Free legs ignore this "
            "(direction-agnostic single system)."
        ),
    )
    args = p.parse_args()

    # CUDA pin (hardware contract — Track A V100 untouched)
    os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda_device
    cuda_devs = [int(d) for d in args.cuda_device.split(",") if d.strip()]

    # Charge-axis gate
    if not args.no_charge_axis_gate:
        from atm_trackB_setup import verify_charge_axis
        axis = verify_charge_axis()
        if not axis["all_residues_within_tol"]:
            print(f"ERROR: charge axis FAIL "
                  f"(max|Sigma q|={axis['max_abs_sigma_q_e']:.2e})",
                  file=sys.stderr)
            return 2

    # Parse displacement
    displ_parts = [float(x) for x in args.displacement_nm.split(",")]
    if len(displ_parts) != 3:
        print(f"ERROR: --displacement-nm needs 3 floats, got {displ_parts}",
              file=sys.stderr)
        return 2
    displacement_nm = tuple(displ_parts)

    out_root = os.path.join(_PROJ_ROOT, args.out_root)
    os.makedirs(out_root, exist_ok=True)

    endpoints = [e.strip() for e in args.endpoints.split(",") if e.strip()]
    legs = [l.strip() for l in args.legs.split(",") if l.strip()]
    for e in endpoints:
        assert e in ("cp4", "wt"), f"unknown endpoint {e}"
    for l in legs:
        assert l in ("bound", "free"), f"unknown leg {l}"

    # smoke restriction
    if args.smoke and "--endpoints" not in sys.argv and "--legs" not in sys.argv:
        endpoints = ["cp4"]
        legs = ["bound"]  # v2.1 smoke targets the formerly-failing bound leg
        print("# smoke mode: restricted to cp4/bound "
              "(override with --endpoints/--legs)")

    # Parse bound-leg directions. The bound leg requires per-direction
    # systems (dplus = bound-state base / dminus = dissociated-state base,
    # binding-site pocket re-hydrated). A combined-only ('dplus') build
    # leaves the dminus walker starting against the wrong (bound) topology,
    # so the dissociated-state base never exists and the alchemical
    # perturbation is ill-defined.
    bound_directions = tuple(
        d.strip() for d in args.bound_directions.split(",") if d.strip()
    )
    if not bound_directions:
        print("ERROR: --bound-directions must be non-empty "
              "(use 'dplus,dminus' for the bound leg)", file=sys.stderr)
        return 2
    for d in bound_directions:
        if d not in ("dplus", "dminus"):
            print(f"ERROR: --bound-directions entry must be 'dplus' or "
                  f"'dminus', got {d!r}", file=sys.stderr)
            return 2

    # Step 1: build 4 systems via upstream wrapper
    build_meta: Dict[str, Any] = {}
    if not args.skip_build:
        from phase4_trackB_v2_make_system import build_all_four_systems
        print(f"\n=== BUILD 4 systems (upstream make_system) ===")
        print(f"    bound_directions = {list(bound_directions)}")
        t0 = time.time()
        build_meta = build_all_four_systems(
            out_root=out_root,
            seed_tag=args.seed_tag,
            displacement_nm=displacement_nm,
            bound_directions=bound_directions,
        )
        print(f"    build wall: {time.time() - t0:.1f} s")
        build_meta["bound_directions"] = list(bound_directions)
        with open(os.path.join(out_root, "build_metadata.json"), "w") as fh:
            json.dump(build_meta, fh, indent=2)
    else:
        print("--skip-build: assuming pre-built systems exist under out-root")

    if args.build_only:
        print("\n--build-only: stopping after system build.")
        return 0

    # Compose verbatim launch command for audit
    launch_command = " ".join(shlex.quote(a) for a in sys.argv)
    env_snapshot = {
        "host": platform.node(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES", ""),
        "atom_openmm_version": _atom_openmm_version(),
    }
    pre_register = {
        "outcomes": list(PRE_REGISTER_OUTCOMES),
        "sign_ssot_kcal": -1.4,
        "sign_ssot_source": "Magotti 2009 ITC (Compstatin Cp4 vs WT)",
        "sign_ssot_regime": "ranking_only",
        "decision_at": "post-uwham-aggregate-of-all-4-legs",
    }

    # Step 2: per-leg cntl + nodefile + structprep
    setup_results: List[Dict[str, Any]] = []
    structprep_results: List[Dict[str, Any]] = []
    production_results: List[Dict[str, Any]] = []
    procs: List[Dict[str, Any]] = []

    for endpoint in endpoints:
        for leg in legs:
            print(f"\n=== SETUP {endpoint} / {leg} ===")
            leg_info = setup_one_leg_v21(
                endpoint=endpoint,
                leg=leg,
                out_root=out_root,
                build_meta=build_meta.get("systems", {}).get(
                    f"{endpoint}/{leg}", {}
                ),
                jobname=args.jobname,
                displacement_nm=displacement_nm,
                cuda_devices=cuda_devs,
                production_steps=args.production_steps,
                prnt_frequency=args.production_steps,
                trj_frequency=args.production_steps * 10,
                max_samples=args.max_samples,
                wall_time_min=args.wall_time_min,
                cycle_time_s=args.cycle_time_s,
                checkpoint_time_s=args.checkpoint_time_s,
                smoke=args.smoke,
                free_schedule=args.free_schedule,
                bound_schedule=args.bound_schedule,
            )
            setup_results.append(leg_info)
            print(f"    cntl: {leg_info['cntl_path']}")
            print(f"    n_ligand: {leg_info['n_ligand_atoms']}  "
                  f"n_restr: {leg_info['n_pos_restrained']}")

            if not args.skip_structprep:
                ready_marker = os.path.join(
                    leg_info["leg_dir"], args.jobname + "_0.xml"
                )
                # Expected mode = "smoke" iff this launch is smoke (smoke
                # cache is only legitimate when intentionally rerunning a
                # smoke gate); otherwise expected_mode = "production"
                # (default). This is the cross-stage drift gate that
                # closes the silent-smoke-reuse bug (cp4/free NaN).
                expected_mode = "smoke" if args.smoke else "production"
                if os.path.isfile(ready_marker):
                    valid, reason = _validate_cached_structprep(
                        leg_dir=leg_info["leg_dir"],
                        expected_mode=expected_mode,
                        basename=args.jobname,
                    )
                    if valid:
                        print(
                            f"    structprep already done (skip): "
                            f"{ready_marker}"
                        )
                        print(f"    cache validation: {reason}")
                        structprep_results.append({
                            "endpoint": endpoint,
                            "leg": leg,
                            "status": "cached",
                            "cache_validation": reason,
                        })
                        continue  # fall through to next leg, do not
                                  # re-run structprep
                    # Fail-fast. Do NOT auto-archive — operator
                    # responsibility (preserve the stale artifact; a
                    # dedicated rerun script handles regeneration).
                    raise IntegrityError(
                        f"Cached structprep validation FAILED for "
                        f"{leg_info['leg_dir']}: {reason}"
                    )
                print(f"--- STRUCTPREP {endpoint} / {leg} ---")
                rc = run_abfe_structprep_for_leg(leg_info)
                if rc != 0:
                    structprep_results.append({
                        "endpoint": endpoint, "leg": leg, "rc": rc,
                        "log": os.path.join(leg_info["leg_dir"],
                                            "_structprep.log"),
                    })
                    print(
                        f"    structprep FAILED rc={rc} — see "
                        f"_structprep.log",
                        file=sys.stderr,
                    )
                    return 3
                # rc == 0: emit sidecar so future launches can validate
                # this artifact (integrity producer side). Wrap in try so a
                # sidecar emission failure (e.g. missing log) does not
                # invalidate a successful structprep — the legacy log
                # fallback inside _validate_cached_structprep handles
                # the no-sidecar case.
                sidecar_path: Optional[str] = None
                sidecar_emit_error: Optional[str] = None
                try:
                    sidecar_path = _emit_structprep_sidecar(
                        leg_dir=leg_info["leg_dir"],
                        mode=expected_mode,
                        basename=args.jobname,
                    )
                    print(f"    sidecar: {sidecar_path}")
                except Exception as exc:  # noqa: BLE001 (defensive)
                    sidecar_emit_error = repr(exc)[:200]
                    print(
                        f"    WARN: sidecar emission failed: "
                        f"{sidecar_emit_error} "
                        f"(legacy log-parse fallback will be used)",
                        file=sys.stderr,
                    )
                structprep_results.append({
                    "endpoint": endpoint, "leg": leg, "rc": rc,
                    "log": os.path.join(leg_info["leg_dir"],
                                        "_structprep.log"),
                    "sidecar": sidecar_path,
                    "sidecar_mode": expected_mode,
                    "sidecar_emit_error": sidecar_emit_error,
                })

    # Step 3 (optional): 30-step bound bidirectional smoke gate
    if args.bidir_smoke_only or args.smoke:
        bound_legs = [li for li in setup_results if li["leg"] == "bound"]
        smoke_results = []
        for li in bound_legs:
            print(f"\n=== BIDIR SMOKE {li['endpoint']} / {li['leg']} (30 steps) ===")
            try:
                sr = bound_leg_bidirectional_smoke(li["leg_dir"], n_steps=30,
                                                   jobname=args.jobname)
                smoke_results.append(sr)
                print(f"    smoke_pass: {sr.get('smoke_pass', False)}")
                for tag, d in sr.get("directions", {}).items():
                    print(f"      {tag}: u0={d.get('u0_kj_per_mol'):.3e}  "
                          f"u1={d.get('u1_kj_per_mol'):.3e}  "
                          f"steps={d.get('completed_steps')}  "
                          f"max_f={d.get('max_force_kcal_per_mol_A')}")
            except Exception as e:
                print(f"    smoke FAILED: {e}", file=sys.stderr)
                smoke_results.append({"endpoint": li["endpoint"],
                                      "error": repr(e)[:300]})
        smoke_path = os.path.join(out_root, "bidir_smoke_results.json")
        with open(smoke_path, "w") as fh:
            json.dump(smoke_results, fh, indent=2)
        print(f"\nBidir smoke results: {smoke_path}")
        if args.bidir_smoke_only:
            ok = all(s.get("smoke_pass") for s in smoke_results)
            print(f"\nSmoke gate: {'PASS' if ok else 'FAIL'}")
            return 0 if ok else 5

    # Run metadata
    run_metadata = {
        "track": "B",
        "version": "v2.1",
        "method": "ATM/ATS ABFE per-endpoint (async_re, Gallicchio 2021); "
                  "upstream make_atm_system_from_rcpt_lig",
        "launcher": "scripts/trackb_production_v2_1_upstream.py",
        "method_ref": "upstream make_system bound-leg shadow setup, 2026-05-31",
        "v2_archive": (
            "outputs/_trackb/_archive/v2_custom_build_partial_20260531/"
        ),
        "v1_archive": (
            "outputs/_trackb/_archive/v1_architectural_failed_20260531/"
        ),
        "endpoints": endpoints,
        "legs": legs,
        "displacement_nm": list(displacement_nm),
        "cuda_visible_devices": os.environ["CUDA_VISIBLE_DEVICES"],
        "platform": args.platform,
        "schedule": {
            "n_states": N_STATES,
            "lambdas_1": LAMBDAS_1,
            "lambdas_2": LAMBDAS_2,
            "directions": DIRECTIONS,
            "intermd": INTERMD,
            "w0_kcal": W0,
            "alpha_per_kcal": ALPHA,
            "u0_kcal": U0,
            "temperature_K": TEMP_K,
        },
        # Per-leg schedule selectors (lambda-densify spec 2026-06-05 +
        # densified_bound28 spec 2026-06-11). free_schedule applies to the
        # free leg only; bound_schedule applies to the bound leg only — they
        # are INDEPENDENT. densified38*/densified_bound28 record the full
        # per-state arrays for audit of the free/bound asymmetry.
        "free_schedule": args.free_schedule,
        "free_schedule_detail": (
            None if args.free_schedule == v2_legacy.DEFAULT_FREE_SCHEDULE
            else v2_legacy.get_schedule(args.free_schedule)
        ),
        "bound_schedule": args.bound_schedule,
        "bound_schedule_detail": (
            None if args.bound_schedule == v2_legacy.DEFAULT_BOUND_SCHEDULE
            else v2_legacy.get_schedule(args.bound_schedule)
        ),
        "production_steps_per_cycle": (
            50 if args.smoke else args.production_steps
        ),
        "max_samples_per_replica": (
            2 if args.smoke else args.max_samples
        ),
        "wall_time_min_per_leg": (
            15 if args.smoke else args.wall_time_min
        ),
        "smoke_mode": bool(args.smoke),
        "launch_command": launch_command,
        "env_snapshot": env_snapshot,
        "pre_register": pre_register,
        "build_metadata_summary": {
            "n_systems": len(build_meta.get("systems", {})),
            "displacement_nm": build_meta.get("displacement_nm"),
        } if build_meta else {},
        "setup_results": setup_results,
        "structprep_results": structprep_results,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "regime": "ranking_only",
    }
    meta_path = os.path.join(out_root, "run_metadata.json")
    if os.path.isfile(meta_path):
        ts_tag = time.strftime("%Y%m%dT%H%M%S")
        meta_path = os.path.join(out_root, f"run_metadata_{ts_tag}.json")
    with open(meta_path, "w") as fh:
        json.dump(run_metadata, fh, indent=2)
    print(f"\nRun metadata: {meta_path}")
    print(f"Launch command (verbatim): {launch_command}")

    if args.setup_only:
        print("\n--setup-only: production NOT launched. To launch:")
        for leg in setup_results:
            print(f"  abfe_production {leg['cntl_path']}  "
                  f"> {os.path.join(leg['leg_dir'], '_production.log')} 2>&1 &")
        return 0

    # Launch production
    for leg_info in setup_results:
        print(f"\n--- PRODUCTION {leg_info['endpoint']} / {leg_info['leg']} ---")
        log_path = os.path.join(leg_info["leg_dir"], "_production.log")
        if args.background:
            proc = run_abfe_production_for_leg(leg_info,
                                               log_path=log_path,
                                               background=True)
            procs.append({
                "endpoint": leg_info["endpoint"],
                "leg": leg_info["leg"],
                "pid": proc.pid,
                "log": log_path,
                "cntl": leg_info["cntl_path"],
            })
            print(f"    BACKGROUND launched PID={proc.pid}  log={log_path}")
        else:
            rc = run_abfe_production_for_leg(leg_info, log_path=log_path)
            production_results.append({
                "endpoint": leg_info["endpoint"],
                "leg": leg_info["leg"],
                "rc": rc,
                "log": log_path,
            })
            print(f"    DONE rc={rc}  log={log_path}")
            if rc != 0:
                print(f"    production FAILED — see _production.log",
                      file=sys.stderr)
                return 4

    if args.background and procs:
        pid_path = os.path.join(out_root, "background_pids.json")
        with open(pid_path, "w") as fh:
            json.dump({"procs": procs,
                       "launched_at": time.strftime("%Y-%m-%dT%H:%M:%S")},
                      fh, indent=2)
        print(f"\nBackground PIDs: {pid_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
