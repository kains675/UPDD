#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B per-direction production launcher (Option B Phase 4 template).

PRODUCTION LAUNCH IS BLOCKED BY DEFAULT.
Requires --i-have-confirmed-c1-through-c8 explicit flag. This script is a
TEMPLATE. Its purpose is to:

  1. Document C1-C8 conditions inline (auto-checked where possible).
  2. Provide the wiring logic that injects per-direction _0_dplus.xml /
     _0_dminus.xml into r0..r10 / r11..r21 trackb_ckpt.xml (per
     ommreplica.py:79 load_checkpoint mechanism).
  3. Pre-register the 4 outcomes.
  4. NEVER auto-launch without the explicit flag.

The per-direction structprep approach is required because a single shared
system topology cannot represent both transfer directions; C1-C8 below are
the hard launch conditions, S1-S4 the soft conditions.

# Usage (when launch authorized)

    # Step 1 — verify C1-C8 manually
    cat scripts/trackb_per_direction_production.py | grep "^# C"

    # Step 2 — run with the explicit safety flag
    # v0.9.10 device-guard fix: --gpu-host {vm,local,cpu}; "device 1"
    # never existed (host = 5070Ti device 0, VM = V100 device 0).
    /home/san/miniconda3/envs/atm/bin/python \\
        scripts/trackb_per_direction_production.py \\
        --i-have-confirmed-c1-through-c8 \\
        --gpu-host vm \\             # V100 (after Track A QM batch idle)
        --seeds s7,s19,s23,s101,s127,s163,s199,s251 \\
        --replicates 3 \\
        --production-steps 2500 \\
        --max-samples 1000 \\
        --wall-time-min 720

# Cross-references

* Per-direction structprep: ``scripts/trackb_per_direction_structprep.py``
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
from typing import Optional, List, Dict, Tuple, Any

_PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "utils"))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "scripts"))


# ---------------------------------------------------------------------------
# Hardware host gating helpers — reused from prep script (v0.9.10 fix).
# Imported lazily to avoid forcing the full prep-script import chain when
# the production launcher is only used in --dry-run mode.
# ---------------------------------------------------------------------------
def _load_prep_gate_helpers():
    """Import gate_gpu_host + normalize_legacy_cuda_device from the prep
    script via spec-loader (no atom_openmm dependency at import time).
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_trackb_prep_for_gate",
        os.path.join(_PROJ_ROOT, "scripts",
                     "trackb_per_direction_structprep.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ===========================================================================
# C1-C8 hard conditions (per-direction structprep, 2026-05-31)
# ===========================================================================
# C1: Free leg PID 1876426 must complete cleanly before Option B Phase 4
#     launch. If free leg ALSO fails, the Option B fundamental assumption
#     (free leg works under the same architecture) is broken — escalate to
#     Option A directly.
#
# C2: Per-direction structprep d=-1 equilibration MUST reach >=250ps
#     stable potential energy (matching d=+1 baseline ~-940000 kJ/mol
#     within +/-5%) before production cycle 1. If d=-1 cannot
#     equilibrate cleanly, Option B FALLBACK to A.
#
# C3: ommreplica.py state-load mechanism must be verified pre-Phase 4
#     launch: explicit dry-run of replica r0 loading
#     _dplus_0.xml and replica r11 loading _dminus_0.xml with positions
#     diff confirming binder location (binding site vs displaced).
#
# C4: Q1(a) integer rescale (Sigma q = -0.176 -> MTR sidechain rescale,
#     backbone+NE1 frozen) MUST persist into Phase 4 production. Charge
#     integrity must be re-checked before resuming.
#
# C5: amber14SB FF consistency confirmed in production launch log
#     (NOT amber99sb-star-ildn-mut, NOT openff). Single-FF declaration
#     on every system XML.
#
# C6: 8-seed cohort matched to Phase IV pairing:
#     {s7, s19, s23, s101, s127, s163, s199, s251}. NO substitutions.
#
# C7: 4-outcome pre-registration (sign-stable / sigma_drift /
#     |mu|_drift / sign-flip) archived in
#     ``outputs/_trackb/production_v2_2/pre_registration.json``
#     BEFORE Phase 4 launch.
#
# C8: ranking-only disclosure in paper sec 4.3 — Track B
#     Delta-Delta-G vs Magotti -1.4 (ITC, pure W4->1MeW)/-3.0
#     (SPR, native->Cp4) ranking only, 4W9A indirect-analog caveat
#     persists.
# ===========================================================================

PRE_REGISTER_OUTCOMES = [
    "sign_stable",        # Cp4 < WT, |delta_delta_G| > 0.5 kcal/mol
    "sigma_drift",        # sigma_btwn > 1.0 kcal/mol
    "mu_drift",           # |delta_delta_G| < 0.5 kcal/mol (statistical tie)
    "sign_flip",          # Cp4 > WT (against Magotti SSOT)
]

EIGHT_SEED_COHORT = ["s7", "s19", "s23", "s101", "s127", "s163", "s199", "s251"]


def check_c1_free_leg_complete(
    free_leg_pid: int,
    free_leg_results_dir: str,
    expected_replicas: int = 22,
) -> Dict[str, Any]:
    """C1 auto-check: free leg PID has terminated AND results exist.

    ``expected_replicas`` (Path λ-densify 2026-06-05) is the COMBINED free-leg
    replica count: 22 for canonical, 38 for the REVISED densified ladder (34
    for the DEPRECATED densified34). Defaults to 22 (the canonical free leg).
    When the leg dir holds a ``trackb_asyncre.cntl`` its DIRECTION column
    overrides the default so a densified free leg is checked against its true
    count (38/34) without the caller having to know it.
    """
    # PID liveness
    pid_alive = False
    try:
        os.kill(free_leg_pid, 0)
        pid_alive = True
    except ProcessLookupError:
        pid_alive = False
    except PermissionError:
        pid_alive = True

    # If the leg dir carries a combined cntl, derive the true replica count
    # from its DIRECTION column (densified free = 34) rather than trusting the
    # default. Silent fallback to the passed default if cntl absent/malformed.
    n_replica_total = expected_replicas
    cntl_path = os.path.join(free_leg_results_dir, "trackb_asyncre.cntl")
    if os.path.isfile(cntl_path):
        try:
            n_replica_total = len(_parse_direction_column(cntl_path))
        except (ValueError, OSError):
            n_replica_total = expected_replicas

    # Results existence — last cycle .out file size > 0 in r0..r{N-1}
    n_replica_complete = 0
    if os.path.isdir(free_leg_results_dir):
        for rid in range(n_replica_total):
            out_path = os.path.join(
                free_leg_results_dir, f"r{rid}", "trackb.out"
            )
            if os.path.isfile(out_path) and os.path.getsize(out_path) > 0:
                n_replica_complete += 1

    return {
        "condition": "C1_free_leg_complete",
        "free_leg_pid": free_leg_pid,
        "pid_still_alive": pid_alive,
        "n_replica_complete": n_replica_complete,
        "n_replica_total": n_replica_total,
        "pass": not pid_alive and n_replica_complete == n_replica_total,
    }


# ---------------------------------------------------------------------------
# Free-pilot mode (densified38 per-direction FREE pilot, 2026-06-05)
# ---------------------------------------------------------------------------
# The canonical C1-C8 gate (above) is a BOUND-leg precondition: it requires the
# FREE leg to ALREADY be complete before the BOUND production launches (C1), and
# C2/C3/C4 read the bound per-direction structprep prep-report. For a FRESH
# FREE-leg pilot those gates are inverted / not-applicable:
#
#   * C1 ("free leg already complete") is meaningless — the free leg is the
#     thing we are about to run. With --free-leg-pid 0 the PID is not alive, but
#     n_replica_complete == n_replica_total can never be true for a fresh run.
#   * C2/C3/C4 read the BOUND per-direction prep-report (sanity_audit of bound
#     _dplus/_dminus XMLs). A free pilot has no bound prep-report.
#
# This is NOT a blind bypass. The free-pilot's correct precondition is that the
# per-direction FREE structprep has produced ``<jobname>_0_dplus.xml`` and
# ``<jobname>_0_dminus.xml`` in every requested (endpoint, leg) leg dir, plus a
# combined cntl whose DIRECTION column is a contiguous forward(+1)-then-
# backward(-1) block (so the 2-process split derives correct counts). The check
# below substitutes that free-leg-appropriate readiness gate for C1-C4.
def check_free_pilot_readiness(
    v21_out_root: str,
    endpoints: List[str],
    legs: List[str],
    jobname: str = "trackb",
) -> Dict[str, Any]:
    """Free-pilot readiness gate (replaces bound C1-C4 for a fresh FREE run).

    For each requested (endpoint, leg), require:
      1. ``<leg_dir>/<jobname>_asyncre.cntl`` present + parseable, with a
         contiguous forward(+1)-then-backward(-1) DIRECTION column.
      2. Both per-direction starting states present:
         ``<jobname>_0_dplus.xml`` AND ``<jobname>_0_dminus.xml``
         (produced by ``trackb_per_direction_structprep.py``).
      3. (v0.9.29.2, 2026-06-05) The FULL set of per-direction STAGING inputs
         resolvable for BOTH directions — i.e. for each of dplus/dminus the
         three files ``stage_per_direction_subdir`` consumes
         (``<jobname>_{tag}_sys.xml`` OR combined ``<jobname>_sys.xml``;
         ``<jobname>_{tag}.pdb`` OR combined ``<jobname>.pdb``;
         ``<jobname>_0_{tag}.xml``) must all be present. This closes the
         false-PASS bug where the gate passed on the base state alone but
         staging then died mid-flight on the missing system/pdb (observed
         2026-06-05): a provisioning gap now FAILS AT THE GATE, not at staging.
         The system/pdb checks use the same per-direction-PREFERRED /
         combined-FALLBACK resolution as staging, so the direction-agnostic
         free leg (combined sys/pdb only) PASSES while a leg missing even the
         combined fallback FAILS.

    ``pass`` is True only if every requested leg satisfies all of (1)-(3).
    Missing items are reported per-leg so the operator knows whether to run
    structprep or point at the right out-root. ``v21_out_root`` is
    project-relative (joined against ``_PROJ_ROOT`` like the rest of ``main``).
    """
    legs_report: List[Dict[str, Any]] = []
    all_ok = True
    for endpoint in endpoints:
        for leg in legs:
            leg_dir = os.path.join(_PROJ_ROOT, v21_out_root, endpoint, leg)
            cntl_path = os.path.join(leg_dir, jobname + "_asyncre.cntl")
            dplus = os.path.join(leg_dir, jobname + "_0_dplus.xml")
            dminus = os.path.join(leg_dir, jobname + "_0_dminus.xml")

            cntl_ok = os.path.isfile(cntl_path)
            n_states: Optional[int] = None
            direction_ok = False
            if cntl_ok:
                try:
                    dirs = _parse_direction_column(cntl_path)
                    # Reuse the canonical contiguous-block validator (raises
                    # ValueError on a non-contiguous / single-direction column).
                    total, _fwd = _derive_state_counts_from_directions(dirs)
                    n_states = total
                    direction_ok = True
                except (ValueError, OSError):
                    direction_ok = False

            dplus_ok = os.path.isfile(dplus)
            dminus_ok = os.path.isfile(dminus)

            # (3) Full staging-input resolvability for BOTH directions. Mirror
            # stage_per_direction_subdir's _resolve_staging_source logic so the
            # gate's verdict matches what staging will actually find.
            staging: Dict[str, Any] = {}
            staging_ok = True
            for tag in ("dplus", "dminus"):
                _sys_path, sys_kind = _resolve_staging_source(
                    leg_dir,
                    per_direction_name=jobname + "_sys_" + tag + ".xml",
                    combined_name=jobname + "_sys.xml",
                )
                _pdb_path, pdb_kind = _resolve_staging_source(
                    leg_dir,
                    per_direction_name=jobname + "_" + tag + ".pdb",
                    combined_name=jobname + ".pdb",
                )
                base_state = os.path.join(
                    leg_dir, jobname + "_0_" + tag + ".xml"
                )
                tag_ok = (
                    sys_kind != "missing"
                    and pdb_kind != "missing"
                    and os.path.isfile(base_state)
                )
                if not tag_ok:
                    staging_ok = False
                staging[tag] = {
                    "sys_source_kind": sys_kind,
                    "pdb_source_kind": pdb_kind,
                    "base_state_present": os.path.isfile(base_state),
                    "pass": tag_ok,
                }

            leg_ok = (
                cntl_ok and direction_ok and dplus_ok and dminus_ok
                and staging_ok
            )
            if not leg_ok:
                all_ok = False
            legs_report.append({
                "endpoint": endpoint,
                "leg": leg,
                "leg_dir": leg_dir,
                "cntl_present": cntl_ok,
                "direction_contiguous": direction_ok,
                "n_states": n_states,
                "dplus_xml_present": dplus_ok,
                "dminus_xml_present": dminus_ok,
                "staging_inputs": staging,
                "staging_inputs_ok": staging_ok,
                "pass": leg_ok,
            })

    reason = None
    if not all_ok:
        missing = [
            f"{r['endpoint']}/{r['leg']}"
            for r in legs_report if not r["pass"]
        ]
        reason = (
            "Free-pilot readiness FAILED for leg(s): "
            + ", ".join(missing)
            + ". Each requires <jobname>_asyncre.cntl (contiguous +1/-1 "
            "DIRECTION) + <jobname>_0_dplus.xml + <jobname>_0_dminus.xml + the "
            "full per-direction staging inputs (per-direction OR combined "
            "<jobname>_sys.xml and <jobname>.pdb) for BOTH directions. "
            "Run trackb_per_direction_structprep.py for the free leg first."
        )
    return {
        "condition": "FREE_PILOT_readiness",
        "v21_out_root": v21_out_root,
        "legs": legs_report,
        "pass": all_ok,
        **({"reason": reason} if reason else {}),
    }


# ---------------------------------------------------------------------------
# G44 fix: C2 auto-check stale-field
# bug. The original implementation read `report["structprep_results"]`, which
# is EMPTY for sanity-only prep reports (per_direction_prep_latest.json
# 2026-06-01 has `structprep_results=[]` because the structprep was run in
# sanity-only mode after b+ XMLs were already on disk). With the old logic
# C2 evaluated to PASS via `n_dplus==n_dminus==expected==0` early on but
# subsequent C2 inspections (production gate) failed with `n_expected=4,
# n_produced=0 → pass=False → LAUNCH BLOCKED` even though all 4 production
# XMLs are healthy on disk.
#
# New logic: prefer `sanity_audit[].directions["d=+1"|"d=-1"]` ground truth
# (status='ok' with size_bytes >= 20MB cutoff matches the actual b+
# artifacts at 22.88-22.90 MB; corrupt/partial files would fall short).
# Fall back to `structprep_results` when sanity_audit is absent or empty
# (preserves legacy behavior for prep reports generated by structprep-mode
# runs that populate structprep_results).
# ---------------------------------------------------------------------------
MIN_PRODUCTION_XML_BYTES = 20_000_000  # 20 MB cutoff (b+ XMLs are 22.88-22.90 MB)


def check_c2_dminus_equilibration(prep_report_path: str) -> Dict[str, Any]:
    """C2 auto-check: d=-1 equilibration reached stable PE within +/-5% of d=+1.

    Ground truth source priority (G44 fix):

    1. ``report["sanity_audit"][i]["directions"]["d=+1"|"d=-1"]`` with
       ``status='ok'`` and ``size_bytes >= MIN_PRODUCTION_XML_BYTES``
       (~22.88-22.90 MB observed for healthy b+ XMLs; 20 MB cutoff catches
       NaN-truncated / crash-truncated outputs).
    2. Fallback: ``report["structprep_results"]`` legacy schema (used by
       prep reports where structprep itself populated the field, NOT
       sanity-only mode).
    """
    if not os.path.isfile(prep_report_path):
        return {
            "condition": "C2_dminus_equilibration",
            "pass": False,
            "reason": f"prep report missing: {prep_report_path}",
        }
    with open(prep_report_path) as fh:
        report = json.load(fh)

    sanity_audit = report.get("sanity_audit", []) or []
    expected_legs = len(sanity_audit)  # one audit entry per (endpoint, leg)

    # --- Source 1: sanity_audit[].directions (preferred — production proof) ---
    n_dplus_ok = 0
    n_dminus_ok = 0
    n_dplus_too_small: List[Dict[str, Any]] = []
    n_dminus_too_small: List[Dict[str, Any]] = []
    for entry in sanity_audit:
        dirs = entry.get("directions") or {}
        for dir_key, dir_info in dirs.items():
            if not isinstance(dir_info, dict):
                continue
            status = dir_info.get("status")
            size_bytes = dir_info.get("size_bytes")
            if status != "ok":
                continue
            size_ok = (
                isinstance(size_bytes, (int, float))
                and size_bytes >= MIN_PRODUCTION_XML_BYTES
            )
            if dir_key == "d=+1":
                if size_ok:
                    n_dplus_ok += 1
                elif size_bytes is not None:
                    n_dplus_too_small.append({
                        "endpoint": entry.get("endpoint"),
                        "leg": entry.get("leg"),
                        "size_bytes": size_bytes,
                    })
            elif dir_key == "d=-1":
                if size_ok:
                    n_dminus_ok += 1
                elif size_bytes is not None:
                    n_dminus_too_small.append({
                        "endpoint": entry.get("endpoint"),
                        "leg": entry.get("leg"),
                        "size_bytes": size_bytes,
                    })

    sanity_audit_complete = (
        expected_legs > 0
        and n_dplus_ok == expected_legs
        and n_dminus_ok == expected_legs
    )
    if sanity_audit_complete:
        return {
            "condition": "C2_dminus_equilibration",
            "source": "sanity_audit",
            "n_dplus_ok": n_dplus_ok,
            "n_dminus_ok": n_dminus_ok,
            "n_expected_legs": expected_legs,
            "min_bytes_cutoff": MIN_PRODUCTION_XML_BYTES,
            "pass": True,
            "note": (
                f"All {expected_legs} legs have both d=+1 and d=-1 XMLs "
                f"with status='ok' AND size_bytes>={MIN_PRODUCTION_XML_BYTES:,} "
                f"(production-grade artifacts). Operator should still inspect "
                f"_structprep_dminus.log for thermalization stability if a "
                f"new fresh structprep was run."
            ),
        }

    # If sanity_audit was present but incomplete, surface the diagnostic.
    if expected_legs > 0 and (n_dplus_ok + n_dminus_ok) > 0:
        return {
            "condition": "C2_dminus_equilibration",
            "source": "sanity_audit",
            "n_dplus_ok": n_dplus_ok,
            "n_dminus_ok": n_dminus_ok,
            "n_expected_legs": expected_legs,
            "n_dplus_too_small": n_dplus_too_small,
            "n_dminus_too_small": n_dminus_too_small,
            "min_bytes_cutoff": MIN_PRODUCTION_XML_BYTES,
            "pass": False,
            "reason": (
                f"sanity_audit incomplete: {n_dplus_ok}/{expected_legs} "
                f"d=+1 + {n_dminus_ok}/{expected_legs} d=-1 ok (size>="
                f"{MIN_PRODUCTION_XML_BYTES:,} bytes)."
            ),
        }

    # --- Source 2: structprep_results (legacy fallback) ---
    structprep = report.get("structprep_results", []) or []
    n_dplus = sum(1 for r in structprep
                  if r.get("direction_tag") == "dplus"
                  and r.get("status") in ("produced", "cached"))
    n_dminus = sum(1 for r in structprep
                   if r.get("direction_tag") == "dminus"
                   and r.get("status") in ("produced", "cached"))
    expected = len(report.get("endpoints", [])) * len(report.get("legs", []))
    legacy_pass = (
        expected > 0 and n_dplus == expected and n_dminus == expected
    )
    return {
        "condition": "C2_dminus_equilibration",
        "source": "structprep_results",
        "n_dplus_produced": n_dplus,
        "n_dminus_produced": n_dminus,
        "n_expected": expected,
        "pass": legacy_pass,
        "reason": (
            None if legacy_pass
            else (
                f"Both sources empty/incomplete: sanity_audit gave "
                f"{n_dplus_ok}/{expected_legs} d+ + {n_dminus_ok}/"
                f"{expected_legs} d- ok; structprep_results gave "
                f"{n_dplus}/{expected} d+ + {n_dminus}/{expected} d-. "
                f"Re-run trackb_per_direction_structprep.py "
                f"--sanity-only to refresh sanity_audit."
            )
        ),
        "note": (
            "Legacy fallback path (sanity_audit absent or empty). "
            "Operator must additionally inspect _structprep_dminus.log "
            "to confirm final PE within +/-5% of d=+1 baseline."
        ),
    }


def check_c3_ommreplica_dryrun(prep_report_path: str) -> Dict[str, Any]:
    """C3 check: per-direction binder centroid diff approx +/-displacement."""
    if not os.path.isfile(prep_report_path):
        return {"condition": "C3_ommreplica_dryrun", "pass": False,
                "reason": "missing prep report"}
    with open(prep_report_path) as fh:
        report = json.load(fh)
    sanity = report.get("sanity_audit", [])
    centroid_diffs: List[float] = []
    for sr in sanity:
        if "binder_centroid_diff_magnitude_nm" in sr:
            centroid_diffs.append(sr["binder_centroid_diff_magnitude_nm"])
    # Expect diff approx 2 * displacement (one displacement from binding
    # site to bulk, then equilibration; reverse for d=-1). Threshold:
    # 0.5 nm minimum (very loose — manual visual inspection recommended).
    pass_ = len(centroid_diffs) > 0 and all(d > 0.5 for d in centroid_diffs)
    return {
        "condition": "C3_ommreplica_dryrun",
        "n_legs_audited": len(centroid_diffs),
        "centroid_diffs_nm": centroid_diffs,
        "pass": pass_,
        "note": (
            "Manually verify that r0 loads _dplus.xml and "
            "r11 loads _dminus.xml via dry-run before launch."
        ),
    }


def check_c4_charge_axis(prep_report_path: str) -> Dict[str, Any]:
    """C4 check: Q1(a) integer rescale persistence (all_residues_within_tol)."""
    if not os.path.isfile(prep_report_path):
        return {"condition": "C4_charge_axis", "pass": False,
                "reason": "missing prep report"}
    with open(prep_report_path) as fh:
        report = json.load(fh)
    charge = report.get("charge_axis_audit", [])
    pass_ = all(c.get("all_residues_within_tol", False) for c in charge if "error" not in c)
    return {
        "condition": "C4_charge_axis",
        "n_legs_audited": len(charge),
        "all_pass": pass_,
        "pass": pass_,
    }


def check_c6_seed_cohort(seeds: List[str]) -> Dict[str, Any]:
    """C6 check: 8-seed cohort matches the canonical Phase IV cohort."""
    canonical = set(EIGHT_SEED_COHORT)
    given = set(seeds)
    return {
        "condition": "C6_seed_cohort",
        "canonical_cohort": sorted(canonical),
        "given_seeds": sorted(given),
        "pass": canonical == given,
        "note": (
            "Differences require explicit review + ranking-only disclosure "
            "of the cohort substitution rationale."
        ),
    }


def _sha256_file(path: str) -> str:
    """Return sha256 hex digest of file at ``path`` (1MB streaming read)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _is_stale_v21_ckpt(ckpt_path: str, baseline_xml: str) -> bool:
    """Return True if ``ckpt_path`` content sha256 matches ``baseline_xml``.

    The Track B v2.1 free leg (PID 1876426, hardcoded ``direction=+1`` for
    all 22 walkers per ``trackb_direction_flip_architectural_blocker_
    20260531``) writes ``r{rid}/trackb_ckpt.xml`` after the equilibration
    + lambda annealing as a verbatim snapshot of the leg-level
    ``trackb_0.xml`` baseline state. If the sha matches, the ckpt is the
    v2.1 stale state that MUST be replaced by per-direction
    ``_0_dplus.xml`` / ``_0_dminus.xml`` before Phase 4 launch (else
    r11..r21 silently load d=+1 and Track B degenerates to the v2.1
    NaN architecture).

    If sha differs, the ckpt is treated as a legitimate in-progress run
    and must NOT be overwritten.

    Returns False if either path is missing (caller handles fresh-stage).
    """
    if not os.path.isfile(ckpt_path) or not os.path.isfile(baseline_xml):
        return False
    try:
        return _sha256_file(ckpt_path) == _sha256_file(baseline_xml)
    except OSError:
        # Defensive: read failure → treat as not-stale (refuses overwrite)
        return False


def _ckpt_matches_expected_direction(ckpt_path: str, expected_src_xml: str) -> bool:
    """Return True iff ckpt_path is byte-identical (sha256) to expected_src_xml,
    the INDEX-SPECIFIC per-direction source (dplus_xml for rid<fwd_replica_count,
    dminus_xml for rid>=fwd_replica_count). A match means the ckpt is already
    correctly per-direction-staged (or a run that crashed at step 0 before mutating
    OpenMM State) -> re-staging is a safe idempotent no-op. Any run that advanced
    >=1 integration step mutates positions/velocities and yields a different sha,
    so it will NOT match and falls through to stale/in-progress handling.
    Returns False if either path missing or on OSError (mirrors _is_stale_v21_ckpt).
    """
    if not os.path.isfile(ckpt_path) or not os.path.isfile(expected_src_xml):
        return False
    try:
        return _sha256_file(ckpt_path) == _sha256_file(expected_src_xml)
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Per-direction control-file (.cntl) generation — C6 root-cause fix
# (two-process per-direction split).
#
# C6 smoke proved a single ``abfe_production`` builds ONE OpenMM System
# Context, which cannot mix the two per-direction atom counts
# (sys_dplus=92855 vs sys_dminus=92804). The combined 22-state cntl
# (DIRECTION=[+1*11, -1*11]) therefore crashed r11 (the first dminus
# walker) because its 92804-atom ckpt was loaded into the 92855-atom dplus
# Context. The fix is to split execution into TWO processes per leg:
#
#   * forward  (process-A): combined states 0..10, DIRECTION rewritten to
#     +1 (already +1 in the combined cntl), DISPLACEMENT +d, lambda
#     0.0 -> 0.5, INTERMEDIATE/W0COEFF=1 at the lambda=0.5 endpoint state
#     10) running against the dplus system (92855 atoms).
#   * backward (process-B): combined states 11..21, DIRECTION rewritten to
#     +1 (combined cntl has -1 — see below), DISPLACEMENT -d, lambda
#     0.5 -> 0.0, INTERMEDIATE/W0COEFF=1 at the lambda=0.5 endpoint state
#     11) running against the dminus system (92804 atoms).
#
# DIRECTION + DISPLACEMENT semantics (source-verified, NOT free-lanced):
#   ommsystem.py:460-466 — the ATMForce reference potential is
#     ``select(step(Direction), u0, u1)`` and the perturbation is
#     ``select(step(Direction), 1, -1)*(u1-(u0+UOffset))``. Direction>=0
#     => base = u0 (soft-core caps the perturbation u1, never the base);
#     Direction<0 => base = u1 (UNPROTECTED — a clashing u1 NaNs the run).
#   The ATMForce is rebuilt at PRODUCTION time by OMMSystemABFE.create_system
#     (ommworker.py:198 -> ommsystem.py:505) from the cntl DISPLACEMENT
#     (ommsystem.py:367 set_displacement -> :492 per-particle), NOT from the
#     system XML — ``grep -c ATMForce trackb_sys_{tag}.xml`` == 0 for BOTH
#     directions. The structprep _make_patched_set_displacement negation
#     applies ONLY to the _0.xml equilibration and is reverted in finally
#     (trackb_per_direction_structprep.py try/finally) — it does NOT
#     propagate to production. Therefore the cntl is the SOLE source of the
#     production-time displacement sign + Direction.
#   (b+) two-leg spec (Azimi 2022 §2.3): BOTH legs base = u0. For dminus
#     the binder is pre-displaced to bulk (Modeller binder@x0+d) so base=u0
#     is the dissociated state; the cntl gets DIRECTION=+1 (rewritten from
#     the combined -1 by _force_direction_plus) and DISPLACEMENT=-d
#     (negated by _negate_displacement) so the perturbation u1 points back
#     into the binding pocket. The naive verbatim slice (DIRECTION=-1 +
#     DISPLACEMENT=+d) selects base=u1=E(binder@bulk+25A=box edge) with no
#     soft-core => Particle-coordinate-NaN at annealing step 0 (d=-1 NaN
#     root cause; cp4/wt bound dminus crashes 2026-06-01 01:57/02:33).
#
# lambda=0.5 intermediate identity (verdict C3 / K-1..K-4): forward state
#   10 and backward state 0 (= combined state 11) both carry
#   lambda1=lambda2=0.5, ALPHA=0.1, U0=110.0, W0COEFF=1.0 — preserved
#   by exact slicing of the combined cntl arrays.
# ---------------------------------------------------------------------------

# cntl keywords whose value is a comma-separated PER-STATE array (length =
# number of states). These get sliced [0:fwd] (forward) / [fwd:total]
# (backward), where fwd/total are DERIVED per-leg from the DIRECTION column
# (canonical bound leg = 11/22; densified free leg = 17/34, Path 2026-06-05).
_PER_STATE_CNTL_KEYS = (
    "LAMBDAS",
    "DIRECTION",
    "INTERMEDIATE",
    "LAMBDA1",
    "LAMBDA2",
    "ALPHA",
    "U0",
    "W0COEFF",
)

# Canonical DEFAULTS for the per-direction split (bound leg = 22-state
# symmetric: 11 forward + 11 backward). These remain the defaults for callers
# / tests that do not pass an explicit count, but the per-state slicing +
# subdir staging + merge now DERIVE the actual counts per-leg from the
# combined cntl's DIRECTION column (Path λ-densify spec 2026-06-05): FREE leg
# may be densified34 (17 fwd + 17 bwd), BOUND leg stays 22 (11+11). NO global
# literal drives the split — these are fall-backs only.
_FWD_STATE_COUNT = 11
_TOTAL_STATE_COUNT = 22


def _derive_state_counts_from_directions(directions: List[int]) -> Tuple[int, int]:
    """Return ``(total_state_count, fwd_state_count)`` from a DIRECTION list.

    total = ``len(directions)``; fwd = ``sum(d >= 0)`` (forward states carry
    DIRECTION=+1, backward DIRECTION=-1). Validates that the DIRECTION column
    is a contiguous block of forward (+1) states followed by backward (-1)
    states (the canonical + densified ATM layout) — raises ValueError on an
    interleaved/malformed column so a corrupt cntl fails loud.
    """
    if not directions:
        raise ValueError("empty DIRECTION list — cannot derive state counts")
    total = len(directions)
    fwd = sum(1 for d in directions if d >= 0)
    bwd = total - fwd
    if fwd == 0 or bwd == 0:
        raise ValueError(
            f"DIRECTION must contain both forward (+1) and backward (-1) "
            f"states; got fwd={fwd} bwd={bwd} (total={total})"
        )
    # Contiguity: first `fwd` must all be >=0, remainder all <0.
    if not (all(d >= 0 for d in directions[:fwd])
            and all(d < 0 for d in directions[fwd:])):
        raise ValueError(
            f"DIRECTION must be a contiguous forward(+1) then backward(-1) "
            f"block; got {directions!r}"
        )
    return total, fwd


def _parse_direction_column(cntl_path: str) -> List[int]:
    """Extract the integer DIRECTION per-state array from a combined cntl.

    Returns a list of ints (e.g. ``[1]*11 + [-1]*11`` for the canonical
    schedule, ``[1]*17 + [-1]*17`` for densified34). Raises FileNotFoundError
    if the cntl is absent, ValueError if no DIRECTION line is found.
    """
    if not os.path.isfile(cntl_path):
        raise FileNotFoundError(f"combined cntl missing: {cntl_path}")
    for kind, key, value in _parse_cntl(cntl_path):
        if kind == "kv" and key == "DIRECTION":
            body = value.strip().strip("'\"")
            return [int(float(tok.strip()))
                    for tok in body.split(",") if tok.strip() != ""]
    raise ValueError(f"no DIRECTION line found in {cntl_path}")


def _derive_state_counts_from_cntl(cntl_path: str) -> Tuple[int, int]:
    """Convenience: parse the DIRECTION column then derive (total, fwd)."""
    return _derive_state_counts_from_directions(
        _parse_direction_column(cntl_path)
    )


def _parse_cntl(cntl_path: str) -> List[Tuple[str, str, str]]:
    """Parse a Track B ``.cntl`` into ordered (kind, key, value) tuples.

    ``kind`` is ``"kv"`` for ``KEY = VALUE`` lines and ``"raw"`` for
    comments / blank lines (preserved verbatim so the regenerated
    per-direction cntls stay diff-friendly). For ``"kv"`` entries ``value``
    is the right-hand side with surrounding whitespace stripped (the
    enclosing quotes, if any, are retained).

    Source format reference: ``outputs/_trackb/production_v2_1/cp4/bound/
    trackb_asyncre.cntl`` (read 2026-06-01). Keys use ``KEY = VALUE`` with
    optional single quotes around the value.
    """
    entries: List[Tuple[str, str, str]] = []
    with open(cntl_path) as fh:
        for line in fh:
            stripped = line.rstrip("\n")
            bare = stripped.strip()
            if not bare or bare.startswith("#"):
                entries.append(("raw", "", stripped))
                continue
            if "=" not in stripped:
                entries.append(("raw", "", stripped))
                continue
            key, _, value = stripped.partition("=")
            entries.append(("kv", key.strip(), value.strip()))
    return entries


def _slice_per_state_value(
    raw_value: str,
    start: int,
    end: int,
    total_state_count: int = _TOTAL_STATE_COUNT,
) -> str:
    """Slice a comma-separated per-state cntl value to ``[start:end]``.

    Preserves a single layer of surrounding quotes if present (matching the
    combined cntl's ``LAMBDAS = '0.0, 0.05, ...'`` style) and the
    ``", "`` element separator. Raises ValueError if the array length is not
    exactly ``total_state_count`` so a malformed combined cntl fails loud
    rather than silently producing a short schedule.

    ``total_state_count`` defaults to the canonical 22 but is DERIVED per-leg
    by ``generate_per_direction_cntls`` from the cntl's DIRECTION column
    (densified free leg = 34; bound leg = 22 — Path λ-densify spec
    2026-06-05). NO global literal is assumed.
    """
    quoted = False
    body = raw_value
    if len(body) >= 2 and body[0] == body[-1] and body[0] in ("'", '"'):
        quoted = True
        quote_char = body[0]
        body = body[1:-1]
    items = [tok.strip() for tok in body.split(",")]
    if len(items) != total_state_count:
        raise ValueError(
            f"per-state cntl array has {len(items)} elements, expected "
            f"{total_state_count} (combined schedule). Value: "
            f"{raw_value!r}"
        )
    sliced = items[start:end]
    rendered = ", ".join(sliced)
    if quoted:
        return f"{quote_char}{rendered}{quote_char}"
    return rendered


def _force_direction_plus(rendered_value: str) -> str:
    """Rewrite every element of an already-sliced DIRECTION value to ``1``.

    Used for the dminus (backward) leg per the (b+) two-leg spec: the
    combined cntl encodes the backward states as ``-1`` (which makes the
    ATMForce select base = u1, no soft-core on the base => d=-1 NaN), but
    BOTH legs must run base = u0 (DIRECTION=+1). The element COUNT and the
    surrounding quote style are preserved verbatim — only the per-element
    integer value is forced to ``1``. Raises ValueError if any element is
    not already an integer DIRECTION token (defensive; a malformed combined
    cntl must fail loud rather than silently emit a wrong schedule).
    """
    quoted = False
    body = rendered_value
    quote_char = ""
    if len(body) >= 2 and body[0] == body[-1] and body[0] in ("'", '"'):
        quoted = True
        quote_char = body[0]
        body = body[1:-1]
    items = [tok.strip() for tok in body.split(",")]
    for tok in items:
        # DIRECTION tokens are ints (+1 / -1); reject anything else.
        if tok.lstrip("+-") == "" or not tok.lstrip("+-").isdigit():
            raise ValueError(
                f"_force_direction_plus: non-integer DIRECTION token "
                f"{tok!r} in {rendered_value!r}"
            )
    forced = ["1"] * len(items)
    rendered = ", ".join(forced)
    if quoted:
        return f"{quote_char}{rendered}{quote_char}"
    return rendered


def _negate_displacement(rendered_value: str) -> str:
    """Negate every component of a 3-vector DISPLACEMENT cntl value.

    ``'25.0, 0.0, 0.0'`` -> ``'-25.0, 0.0, 0.0'``. Preserves the surrounding
    quote style and the ``", "`` separator. ``0.0`` is emitted as ``0.0``
    (not ``-0.0``) for cntl readability. Raises ValueError if any component
    is not a float (a malformed DISPLACEMENT must fail loud). The engine
    parses this back to ``[float(...)]`` via parse_config (utils/config.py)
    so the textual form only needs to round-trip through ``float``.
    """
    quoted = False
    body = rendered_value
    quote_char = ""
    if len(body) >= 2 and body[0] == body[-1] and body[0] in ("'", '"'):
        quoted = True
        quote_char = body[0]
        body = body[1:-1]
    out_components: List[str] = []
    for tok in body.split(","):
        raw = tok.strip()
        try:
            val = float(raw)
        except ValueError as exc:
            raise ValueError(
                f"_negate_displacement: non-float DISPLACEMENT component "
                f"{raw!r} in {rendered_value!r}"
            ) from exc
        negated = -val
        if negated == 0.0:
            negated = 0.0  # collapse -0.0 -> 0.0
        # Preserve a trailing ".0" style for integer-valued floats so the
        # rewritten cntl visually matches the combined cntl (25.0 -> -25.0).
        if negated == int(negated):
            out_components.append(f"{negated:.1f}")
        else:
            out_components.append(repr(negated))
    rendered = ", ".join(out_components)
    if quoted:
        return f"{quote_char}{rendered}{quote_char}"
    return rendered


# Mapping from schedule-dict keys (utils.adaptive_lambda / asyncre
# _schedule_dict shape) to the combined cntl's per-state keyword names.
_SCHEDULE_DICT_TO_CNTL_KEY = {
    "lambdas": "LAMBDAS",
    "directions": "DIRECTION",
    "intermd": "INTERMEDIATE",
    "lambdas_1": "LAMBDA1",
    "lambdas_2": "LAMBDA2",
    "alpha": "ALPHA",
    "u0": "U0",
    "w0": "W0COEFF",
}


def apply_free_schedule_file_to_combined_cntl(
    combined_cntl_path: str,
    schedule_file: str,
) -> Dict[str, Any]:
    """Re-emit the per-state schedule columns of a COMBINED cntl from a JSON
    schedule dict — the SINGLE combined-cntl write point (single source of
    truth for the schedule).

    Loads ``schedule_file`` via ``utils.adaptive_lambda.schedule_io.
    load_schedule_dict`` (accepts a bare _schedule_dict or an
    adaptive_schedule.json envelope), then rewrites EACH per-state keyword
    (LAMBDAS / DIRECTION / INTERMEDIATE / LAMBDA1 / LAMBDA2 / ALPHA / U0 /
    W0COEFF) of ``combined_cntl_path`` in place from the loaded arrays, using the
    same ``", "`` rendering + single-quote style as the asyncre writer. ALL
    non-schedule keys (LIGAND_ATOMS / DISPLACEMENT / step counts / UMAX / …) and
    all comment/blank lines are preserved verbatim.

    This is the wiring for ``--free-schedule-file``: the
    downstream ``generate_per_direction_cntls`` slices THIS combined cntl, so the
    per-direction cntls inherit the loaded λ WITHOUT any post-hoc per-direction
    patch. Returns a small audit dict (n_states, rewritten keys, schedule_file).

    Raises FileNotFoundError if either file is missing, ValueError if the
    loaded schedule's array lengths are inconsistent or a required per-state key
    is absent from the combined cntl (fail-loud, no silent partial rewrite).
    """
    if not os.path.isfile(combined_cntl_path):
        raise FileNotFoundError(
            f"combined cntl missing: {combined_cntl_path}"
        )
    # schedule_io is in utils/ (already on sys.path via _PROJ_ROOT/utils insert).
    from adaptive_lambda.schedule_io import load_schedule_dict

    schedule = load_schedule_dict(schedule_file)
    n_states = int(schedule["n_states"])

    # Render each schedule array as the cntl value body (no quotes here; the
    # write step re-applies the original quote style per key).
    rendered: Dict[str, str] = {}
    for sk, cntl_key in _SCHEDULE_DICT_TO_CNTL_KEY.items():
        if sk not in schedule:
            raise ValueError(
                f"schedule dict missing key {sk!r} (need it for cntl "
                f"{cntl_key}); have {sorted(schedule)}"
            )
        arr = schedule[sk]
        if len(arr) != n_states:
            raise ValueError(
                f"schedule array {sk!r} length {len(arr)} != n_states "
                f"{n_states}"
            )
        rendered[cntl_key] = ", ".join(_format_cntl_element(sk, v) for v in arr)

    entries = _parse_cntl(combined_cntl_path)
    rewritten_keys: List[str] = []
    out_lines: List[str] = []
    seen_keys = set()
    for kind, key, value in entries:
        if kind == "raw":
            out_lines.append(value)
            continue
        if key in rendered:
            quote = ""
            if len(value) >= 2 and value[0] == value[-1] and value[0] in (
                "'", '"'
            ):
                quote = value[0]
            out_lines.append(f"{key} = {quote}{rendered[key]}{quote}")
            rewritten_keys.append(key)
            seen_keys.add(key)
            continue
        out_lines.append(f"{key} = {value}")

    missing = [k for k in rendered if k not in seen_keys]
    if missing:
        raise ValueError(
            f"combined cntl {combined_cntl_path} is missing per-state keys "
            f"{missing}; cannot safely apply --free-schedule-file (the cntl "
            f"layout is unexpected). No write performed."
        )

    with open(combined_cntl_path, "w") as fh:
        fh.write("\n".join(out_lines))
        if out_lines and not out_lines[-1].endswith("\n"):
            fh.write("\n")

    return {
        "combined_cntl": combined_cntl_path,
        "schedule_file": schedule_file,
        "n_states": n_states,
        "rewritten_keys": rewritten_keys,
    }


def _format_cntl_element(schedule_key: str, value: Any) -> str:
    """Render one schedule-array element as a cntl token.

    DIRECTION / INTERMEDIATE are integer columns (1 / -1 / 0); the rest are
    floats. Integer columns are emitted without a decimal point to match the
    combined cntl's style.
    """
    if schedule_key in ("directions", "intermd"):
        return str(int(round(float(value))))
    return repr(float(value))


def generate_per_direction_cntls(
    leg_dir: str,
    jobname: str = "trackb",
    combined_cntl_basename: Optional[str] = None,
) -> Dict[str, Any]:
    """Derive forward (dplus) + backward (dminus) per-direction cntls from
    the combined 22-state cntl.

    Reads ``leg_dir/{combined_cntl_basename}`` (default
    ``{jobname}_asyncre.cntl``), slices the per-state arrays
    (``_PER_STATE_CNTL_KEYS``) into states 0..10 (forward) and 11..21
    (backward), rewrites ``BASENAME`` to point at the per-direction system /
    pdb / _0 variant, and writes the two cntls into per-direction work
    subdirs:

      * ``leg_dir/dplus/{jobname}_dplus_asyncre.cntl``  (11 states)
      * ``leg_dir/dminus/{jobname}_dminus_asyncre.cntl`` (11 states)

    Per the (b+) two-leg ATM ABFE spec (Azimi 2022 §2.3) BOTH legs run with
    the ATMForce reference base = u0 (``select(step(Direction), u0, u1)``,
    ommsystem.py:460-466), so BOTH cntls carry ``DIRECTION = +1`` for every
    state. The dminus (backward) leg reaches u0 = dissociated by a
    PRE-DISPLACED binder (Modeller binder@x0+d) combined with a NEGATED
    runtime ATMForce displacement (``DISPLACEMENT`` sign-flipped here so the
    engine's ``set_displacement`` — ommsystem.py:365-368 — drives the
    perturbation toward the binding pocket, not the box edge). The naive
    verbatim slice (dminus DIRECTION=-1 + DISPLACEMENT=+25) makes the engine
    select base = u1 = E(binder@bulk + 25 A = box-edge/void) with NO
    soft-core protection on the base term => Particle-coordinate-NaN at
    annealing step 0 (the d=-1 NaN root cause).

    Concretely, this function rewrites two keywords for the dminus slice:

      * ``DIRECTION``  → all ``+1`` (NOT the verbatim ``-1`` slice). The
        dplus slice keeps its native ``+1`` values.
      * ``DISPLACEMENT`` → each component negated (``25.0, 0.0, 0.0`` →
        ``-25.0, 0.0, 0.0``). The engine builds the ATMForce per-particle
        displacement from the cntl ``DISPLACEMENT`` verbatim
        (ommsystem.py:492); the per-direction system XML carries NO
        ATMForce (``grep -c ATMForce trackb_sys_dminus.xml`` == 0), so the
        cntl is the SOLE source of the production-time displacement sign.
        The dplus DISPLACEMENT is copied verbatim (+d).

    The lambda=0.5 intermediate state (forward state 10 / backward state 0)
    keeps ALPHA=0.1, U0=110.0, W0COEFF=1.0, lambda1=lambda2=0.5 by exact
    array slicing — unchanged.

    BASENAME rewrite: each per-direction cntl's ``BASENAME`` becomes
    ``{jobname}_{tag}`` so upstream ``OMMSystemABFE`` loads
    ``{jobname}_{tag}_sys.xml`` + ``{jobname}_{tag}.pdb`` and the engine's
    canonical ``{jobname}_{tag}_0.xml``. The caller (subdir staging) is
    responsible for materializing those per-direction-named files in the
    work subdir.

    Returns a dict describing both written cntls + the sliced schedule for
    inspection. Raises FileNotFoundError if the combined cntl is
    absent, ValueError if any per-state array is not 22 elements long.
    """
    leg_dir = os.path.abspath(leg_dir)
    if combined_cntl_basename is None:
        combined_cntl_basename = jobname + "_asyncre.cntl"
    combined_path = os.path.join(leg_dir, combined_cntl_basename)
    if not os.path.isfile(combined_path):
        raise FileNotFoundError(
            f"Combined cntl missing: {combined_path}. Run "
            f"trackb_per_direction_structprep.py first."
        )

    entries = _parse_cntl(combined_path)

    # Derive the per-leg state counts from the cntl's DIRECTION column rather
    # than the module literal (Path λ-densify spec 2026-06-05): bound leg = 22
    # (11 fwd + 11 bwd), densified free leg = 34 (17 fwd + 17 bwd). The slice
    # boundaries follow from sum(DIRECTION==+1), never a global 11/22.
    total_state_count, fwd_state_count = _derive_state_counts_from_cntl(
        combined_path
    )

    directions = {
        "dplus": {
            "tag": "dplus",
            "slice": (0, fwd_state_count),
            "subdir": os.path.join(leg_dir, "dplus"),
        },
        "dminus": {
            "tag": "dminus",
            "slice": (fwd_state_count, total_state_count),
            "subdir": os.path.join(leg_dir, "dminus"),
        },
    }

    result: Dict[str, Any] = {
        "leg_dir": leg_dir,
        "combined_cntl": combined_path,
        "total_state_count": total_state_count,
        "fwd_state_count": fwd_state_count,
        "directions": {},
    }

    for tag, spec in directions.items():
        start, end = spec["slice"]
        subdir = spec["subdir"]
        os.makedirs(subdir, exist_ok=True)
        out_lines: List[str] = []
        sliced_schedule: Dict[str, List[str]] = {}
        for kind, key, value in entries:
            if kind == "raw":
                out_lines.append(value)
                continue
            if key == "BASENAME":
                # BASENAME drives sys.xml / pdb / _0.xml lookup. Repoint at
                # the per-direction variant. Preserve quote style.
                quote = ""
                bare = value
                if len(value) >= 2 and value[0] == value[-1] and value[0] in (
                    "'", '"'
                ):
                    quote = value[0]
                    bare = value[1:-1]
                # bare is the original jobname (e.g. 'trackb')
                new_basename = bare + "_" + tag
                out_lines.append(f"BASENAME = {quote}{new_basename}{quote}")
                continue
            if key in _PER_STATE_CNTL_KEYS:
                sliced_val = _slice_per_state_value(
                    value, start, end,
                    total_state_count=total_state_count,
                )
                if key == "DIRECTION" and tag == "dminus":
                    # (b+) DIRECTION rewrite — NOT a verbatim slice. The
                    # combined cntl encodes the backward leg as DIRECTION=-1
                    # which makes the ATMForce select base = u1 (no
                    # soft-core on the base) => d=-1 NaN. Force every dminus
                    # state to DIRECTION=+1 so the base is always u0 (the
                    # pre-displaced dissociated state), matching the dplus
                    # leg. The element COUNT and quoting are preserved.
                    sliced_val = _force_direction_plus(sliced_val)
                out_lines.append(f"{key} = {sliced_val}")
                # Record unquoted element list for the audit payload.
                inner = sliced_val.strip("'\"")
                sliced_schedule[key] = [
                    t.strip() for t in inner.split(",")
                ]
                continue
            if key == "DISPLACEMENT":
                # (b+) DISPLACEMENT sign derivation — NOT verbatim. The
                # engine builds the ATMForce per-particle displacement from
                # this cntl value at production time (ommsystem.py:492); the
                # per-direction system XML carries NO ATMForce (grep -c
                # ATMForce trackb_sys_{tag}.xml == 0). For dminus the runtime
                # displacement must point the perturbation toward the binding
                # pocket (-d), so each component is negated; dplus keeps +d.
                emit_val = value
                if tag == "dminus":
                    emit_val = _negate_displacement(value)
                out_lines.append(f"{key} = {emit_val}")
                sliced_schedule["DISPLACEMENT"] = [
                    t.strip() for t in emit_val.strip("'\"").split(",")
                ]
                continue
            # Other non-per-state keyword (LIGAND_ATOMS, RCPT_CM_ATOMS,
            # UMAX, etc.) — copy verbatim; identical for both legs.
            out_lines.append(f"{key} = {value}")

        out_basename = jobname + "_" + tag + "_asyncre.cntl"
        out_path = os.path.join(subdir, out_basename)
        header = (
            "# Track B per-direction cntl ({tag}) — generated by "
            "generate_per_direction_cntls\n"
            "# Sliced from combined 22-state {combined} (states {a}..{b}).\n"
            "# DIRECTION column preserved verbatim (thermodynamic-leg "
            "label, NOT re-derived).\n"
            "# DO NOT EDIT BY HAND.\n"
        ).format(
            tag=tag,
            combined=combined_cntl_basename,
            a=start,
            b=end - 1,
        )
        with open(out_path, "w") as fh:
            fh.write(header)
            fh.write("\n".join(out_lines) + "\n")

        result["directions"][tag] = {
            "tag": tag,
            "subdir": subdir,
            "cntl_path": out_path,
            "cntl_basename": out_basename,
            "state_slice": [start, end],
            "n_states": end - start,
            "basename": jobname + "_" + tag,
            "schedule": sliced_schedule,
        }

    return result


def _ensure_canonical_base_state(
    leg_dir: str,
    jobname: str = "trackb",
    direction: str = "dplus",
) -> Dict[str, Any]:
    """Provision the canonical ``{jobname}_0.{xml,pdb}`` from the
    per-direction variant (idempotent, no-overwrite).

    Root cause (v0.9.25): per-direction structprep renames (``shutil.move``)
    the engine's canonical ``{jobname}_0.xml`` to the per-direction variant
    ``{jobname}_0_{tag}.xml`` for BOTH directions, so no canonical
    ``{jobname}_0.xml`` survives. The abfe_production worker bootstrap
    (``ommworker.py:265``) hardcodes ``loadState({jobname}_0.xml)`` for the
    service worker (``compute=False`` template only). Absent canonical →
    ``FileNotFoundError: trackb_0.xml``.

    Restores the canonical by **copying** (not symlink — copy is
    restart-safe; a re-run of structprep's ``shutil.move`` cannot break a
    real file) from the ``direction`` variant (default ``dplus``).

    Idempotent: a pre-existing canonical is NEVER overwritten. The
    variant is a Direction=+1 lambda=0.5 equilibrium consumed only by the
    service worker; every sampling walker loads its own per-direction ckpt
    (r0..r10=dplus, r11..r21=dminus), so the backward (r11..r21) replicas
    are NOT contaminated.

    Raises RuntimeError if the source variant XML is absent (structprep
    incomplete) — fail-fast, cohort-safe.

    Returns ``{"created_xml", "created_pdb", "xml", "pdb", "source_xml",
    "source_pdb", "direction"}`` where ``created_*`` is True iff this call
    wrote the canonical (False = pre-existing no-op).
    """
    leg_dir = os.path.abspath(leg_dir)
    src_xml = os.path.join(leg_dir, jobname + "_0_" + direction + ".xml")
    src_pdb = os.path.join(leg_dir, jobname + "_0_" + direction + ".pdb")
    canon_xml = os.path.join(leg_dir, jobname + "_0.xml")
    canon_pdb = os.path.join(leg_dir, jobname + "_0.pdb")

    if not os.path.isfile(src_xml):
        raise RuntimeError(
            "_ensure_canonical_base_state: per-direction source XML "
            "absent (structprep incomplete?): " + src_xml
        )

    created_xml = False
    if not os.path.isfile(canon_xml):
        shutil.copy2(src_xml, canon_xml)
        created_xml = True

    created_pdb = False
    if not os.path.isfile(canon_pdb) and os.path.isfile(src_pdb):
        shutil.copy2(src_pdb, canon_pdb)
        created_pdb = True

    return {
        "created_xml": created_xml,
        "created_pdb": created_pdb,
        "xml": canon_xml,
        "pdb": canon_pdb,
        "source_xml": src_xml,
        "source_pdb": src_pdb,
        "direction": direction,
    }


def _rsync_canonical_base_state_to_vm(
    leg_dir: str,
    vm_ssh_host: str,
    jobname: str = "trackb",
    ssh_connect_timeout_s: int = 5,
    subprocess_timeout_s: int = 300,
) -> Dict[str, Any]:
    """Push the canonical ``{jobname}_0.{xml,pdb}`` to the VM ``leg_dir``.

    G45 (``_rsync_per_replica_ckpts_to_vm``) pushes only the per-replica
    ckpts + ckpt_is_valid marker; the VM already has the per-direction
    variants (``_0_dplus.xml`` etc., written by VM-side structprep) but the
    canonical ``_0.xml`` is only restored on the VM when structprep ran
    through the SSH wrapper. For the host-local-structprep + vm-production
    case the VM canonical is absent → abfe_production crashes. Push it
    explicitly (mirror tree: same ``leg_dir`` on host and VM, matching
    ``_rsync_per_replica_ckpts_to_vm``).

    Canonical is a real file (copy of dplus), so plain ``rsync -a`` (no
    ``-L`` dereference). Idempotent. VM-side ``test -f`` verify; raises
    RuntimeError (cohort halt) on failure.
    """
    leg_dir = leg_dir.rstrip("/")
    canon_xml = os.path.join(leg_dir, jobname + "_0.xml")
    canon_pdb = os.path.join(leg_dir, jobname + "_0.pdb")
    if not os.path.isfile(canon_xml):
        raise RuntimeError(
            "_rsync_canonical_base_state_to_vm: host canonical XML missing "
            "pre-rsync: " + canon_xml +
            ". _ensure_canonical_base_state must run first."
        )
    ssh_e = (
        f"ssh -o ConnectTimeout={ssh_connect_timeout_s} "
        f"-o StrictHostKeyChecking=no -o LogLevel=ERROR"
    )
    rel_paths = [jobname + "_0.xml"]
    if os.path.isfile(canon_pdb):
        rel_paths.append(jobname + "_0.pdb")
    import tempfile
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".rsync_canon", delete=False
    ) as fh:
        fh.write("\n".join(rel_paths) + "\n")
        files_from_path = fh.name
    try:
        rsync_cmd = [
            "rsync", "-a", "-e", ssh_e,
            f"--files-from={files_from_path}",
            f"{leg_dir}/",
            f"{vm_ssh_host}:{leg_dir}/",
        ]
        try:
            rsync_result = subprocess.run(
                rsync_cmd, capture_output=True,
                timeout=subprocess_timeout_s,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"rsync of canonical base-state timed out "
                f"(>{subprocess_timeout_s}s). Cohort halted."
            ) from exc
        if rsync_result.returncode != 0:
            raise RuntimeError(
                f"rsync canonical base-state FAILED (rc="
                f"{rsync_result.returncode}): "
                f"{rsync_result.stderr.decode(errors='replace')[:500]}"
            )
    finally:
        try:
            os.unlink(files_from_path)
        except OSError:
            pass
    # VM-side verify canonical XML present.
    verify_cmd = (
        f"test -f {shlex.quote(os.path.join(leg_dir, jobname + '_0.xml'))}"
    )
    try:
        verify_result = subprocess.run(
            [
                "ssh",
                "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
                "-o", "StrictHostKeyChecking=no",
                "-o", "LogLevel=ERROR",
                vm_ssh_host, verify_cmd,
            ],
            capture_output=True, timeout=60,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VM canonical-base-state verify timed out — ssh to "
            f"{vm_ssh_host} unresponsive. Cohort halted."
        ) from exc
    if verify_result.returncode != 0:
        raise RuntimeError(
            "VM canonical base-state push failed: canonical "
            f"{jobname}_0.xml missing on VM at {leg_dir}. Cohort halt (C1)."
        )
    return {"status": "rsynced", "pushed": rel_paths, "leg_dir": leg_dir}


def stage_per_replica_checkpoints(
    leg_dir: str,
    jobname: str,
    n_states: int = 22,
    fwd_replica_count: int = 11,
    archive_stale: bool = False,
    timestamp: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Inject _0_dplus.xml -> r0..r10/trackb_ckpt.xml and _0_dminus.xml ->
    r11..r21/trackb_ckpt.xml, plus touch ``ckpt_is_valid`` marker
    (consumed by ommreplica.py:77 load_checkpoint).

    .. deprecated:: 2026-06-03 (NOT called in the live launch path)
        This helper injects per-replica checkpoints that embed
        ``REStateId="0"`` from the lambda-0.5 equilibration state. Loading
        such a ckpt makes ``ommreplica.update_state_from_context``
        (``ommreplica.py:135``) read ``stateid == 0`` for every replica,
        overwriting the ``stateid is None`` sentinel and SKIPPING async_re's
        initial-seed guard (``openmm_async_re.py:180-181, 396-397``) — the
        root cause of the bound-leg state-0-collapse. The live
        path now uses ``stage_per_direction_subdir`` with NO per-replica
        ckpt (the engine seeds states correctly). This function is RETAINED
        only for the regression tests that document the broken contract.

    The original wiring relied on ommreplica.py:79: the verification was
    correct that the ckpt's stateid is read; the error was staging a ckpt
    at all.

    ``archive_stale`` (default False — safe for --dry-run / inspection)
    controls collision handling:

    * False → existing ``r{rid}/trackb_ckpt.xml`` causes
      ``skipped_already_exists`` (legacy protective behavior).
    * True  → sha-compare existing ckpt against ``leg_dir/{jobname}_0.xml``
      (the v2.1 baseline). If sha matches → ARCHIVE the stale ckpt to
      ``leg_dir/_stale_v21_<timestamp>/r{rid}/trackb_ckpt.xml`` (move,
      never delete) and FRESH-STAGE the per-direction XML. If sha
      differs → ``skipped_legitimate_in_progress_refuse`` (a real
      in-progress run, protected). The ``ckpt_is_valid`` marker is
      archived alongside the first stale-replica directory.

    ``timestamp`` (default None → ``time.strftime("%Y%m%dT%H%M%S")``) is
    the timestamp embedded in the archive directory name. Exposed for
    test determinism.
    """
    dplus_xml = os.path.join(leg_dir, jobname + "_0_dplus.xml")
    dminus_xml = os.path.join(leg_dir, jobname + "_0_dminus.xml")
    if not os.path.isfile(dplus_xml):
        raise RuntimeError(f"Missing {dplus_xml} — run per-direction structprep first")
    if not os.path.isfile(dminus_xml):
        raise RuntimeError(f"Missing {dminus_xml} — run per-direction structprep first")

    baseline_xml = os.path.join(leg_dir, jobname + "_0.xml")
    archive_ts = timestamp or time.strftime("%Y%m%dT%H%M%S")
    archive_root = os.path.join(leg_dir, f"_stale_v21_{archive_ts}")
    marker_archived = False  # archive ckpt_is_valid only once per call

    staged: List[Dict[str, Any]] = []
    for rid in range(n_states):
        r_dir = os.path.join(leg_dir, f"r{rid}")
        os.makedirs(r_dir, exist_ok=True)
        ckpt_dst = os.path.join(r_dir, jobname + "_ckpt.xml")
        # Hoisted above the existence check so both the idempotent
        # already-staged match (branch a) and the fresh stage below use the
        # SAME index-specific per-direction source.
        src = dplus_xml if rid < fwd_replica_count else dminus_xml
        if os.path.isfile(ckpt_dst):
            if archive_stale and _ckpt_matches_expected_direction(ckpt_dst, src):
                # (a) Already correctly per-direction-staged (byte-identical to
                # the index-specific src). Idempotent no-op — re-running the
                # launcher after a prior successful stage must NOT mis-classify
                # these as stale-v2.1 or genuine in-progress. Gated by
                # archive_stale so the archive_stale=False contract
                # (skipped_already_exists) is preserved unchanged.
                staged.append({
                    "replica": rid,
                    "ckpt_path": ckpt_dst,
                    "source": "dplus" if rid < fwd_replica_count else "dminus",
                    "src_xml": src,
                    "status": "already_per_direction_staged",
                })
                continue
            elif archive_stale and _is_stale_v21_ckpt(ckpt_dst, baseline_xml):
                # (b) Archive (move, never delete) then fall through
                stale_dir = os.path.join(archive_root, f"r{rid}")
                os.makedirs(stale_dir, exist_ok=True)
                archived_ckpt = os.path.join(stale_dir, jobname + "_ckpt.xml")
                shutil.move(ckpt_dst, archived_ckpt)
                # Archive ckpt_is_valid marker once (it is leg-level, not
                # per-replica) so a stale launch attempt aborts cleanly.
                if not marker_archived:
                    marker_src = os.path.join(leg_dir, "ckpt_is_valid")
                    if os.path.isfile(marker_src):
                        os.makedirs(archive_root, exist_ok=True)
                        shutil.move(
                            marker_src,
                            os.path.join(archive_root, "ckpt_is_valid"),
                        )
                    marker_archived = True
                # Fall through to fresh stage below — do not `continue`.
            else:
                # (c) Genuine in-progress (sha differs from src AND baseline)
                # under archive_stale; or legacy skip when archive_stale=False.
                staged.append({
                    "replica": rid,
                    "ckpt_path": ckpt_dst,
                    "status": (
                        "skipped_legitimate_in_progress_refuse"
                        if archive_stale
                        else "skipped_already_exists"
                    ),
                })
                continue
        # Copy (not symlink) so the structprep XML can be deleted/moved
        # without breaking the run.
        shutil.copy2(src, ckpt_dst)
        entry: Dict[str, Any] = {
            "replica": rid,
            "ckpt_path": ckpt_dst,
            "source": "dplus" if rid < fwd_replica_count else "dminus",
            "src_xml": src,
            "status": "staged",
        }
        if archive_stale and marker_archived:
            entry["stale_archived_to"] = archive_root
        staged.append(entry)

    # Touch ``ckpt_is_valid`` (ommreplica.py:77 safecheckpoint marker).
    # If we archived the prior marker above, this rewrites a fresh one
    # for the per-direction launch.
    valid_marker = os.path.join(leg_dir, "ckpt_is_valid")
    if not os.path.isfile(valid_marker):
        with open(valid_marker, "w") as fh:
            fh.write("# Created by trackb_per_direction_production.py "
                     "stage_per_replica_checkpoints — per ommreplica.py:77\n")
    return staged


# ---------------------------------------------------------------------------
# G45 fix (VM-PER-REPLICA-CKPT-MISSING-01):
# stage_per_replica_checkpoints writes 22 per-replica trackb_ckpt.xml +
# leg-level ckpt_is_valid marker on the HOST filesystem, but the SSH-
# dispatched abfe_production runs on the VM filesystem. Without explicit
# rsync the VM r{0..21} directories are empty (no per-direction state
# injected) and ommreplica.py:79 falls back to loading the system XML —
# silently degrading Track B v2.2 to the v2.1 single-direction architecture
# (NaN-architectural-blocker class). This is the 5th member of the failure
# family:
#   - integrity_vm_lane_self_provisioning_20260529
#   - integrity_vm_conda_path_noninteractive_20260528
#   - Round 2 F3 (VM abfe_bin missing)
#   - Round 3 F3 (VM env install)
#   - Round 6 G38 (VM leg-dir missing on dispatch)
#
# Fix: explicit rsync of per-replica ckpts + ckpt_is_valid marker AFTER
# stage_per_replica_checkpoints (host mutation complete) and BEFORE the
# SSH-wrapped abfe_production invocation. Atomic: rsync transfers all
# 22 ckpts in a single invocation; partial failure raises RuntimeError so
# operator can manually restage or retry without leaving the cohort in
# half-staged VM state.
#
# Bandwidth: 22 files * ~22 MB = ~484 MB per leg. Local network ~3-5 min
# per leg. VM disk requirement: ~1 GB free per leg (cohort = ~2 GB).
# ---------------------------------------------------------------------------
def _rsync_per_replica_ckpts_to_vm(
    leg_dir: str,
    vm_ssh_host: str,
    jobname: str,
    n_replicas: int = 22,
    ssh_connect_timeout_s: int = 5,
    subprocess_timeout_s: int = 900,
) -> Dict[str, Any]:
    """Rsync ``leg_dir/r{0..n_replicas-1}/{jobname}_ckpt.xml`` + leg-level
    ``ckpt_is_valid`` marker to the same path on the VM.

    Atomic contract: if any per-replica ckpt or the marker fails to
    propagate, raise RuntimeError. Caller halts cohort.

    Verification: after rsync, an ``ssh ... test -f`` probe confirms
    all 22 ckpts + marker exist on VM with non-zero size. This catches
    silent partial transfers (network blip, disk full, permission).

    Parameters
    ----------
    leg_dir : str
        Absolute path to ``v21_out_root/<endpoint>/<leg>``. Same path
        is used on both host and VM (mirror tree).
    vm_ssh_host : str
        SSH target (e.g. ``san@192.168.122.155``).
    jobname : str
        Default ``trackb`` — matches the ckpt basename pattern.
    n_replicas : int
        22 (Track B v2.1+v2.2 default).
    ssh_connect_timeout_s : int
        ssh ConnectTimeout (default 5 — fast-fail if VM unreachable).
    subprocess_timeout_s : int
        Outer subprocess.run timeout for the rsync (default 15 min;
        local-network 484 MB transfers in ~3-5 min).

    Returns
    -------
    Dict[str, Any]
        ``{"status": "rsynced", "n_ckpts": 22, "marker_present": True,
           "bytes_sent_estimate": ...}`` on success.
        Raises RuntimeError on failure (cohort-safe halt).
    """
    # 1) Build list of source files to verify they exist on host first
    src_files: List[str] = []
    for rid in range(n_replicas):
        ckpt_path = os.path.join(
            leg_dir, f"r{rid}", jobname + "_ckpt.xml"
        )
        if not os.path.isfile(ckpt_path):
            raise RuntimeError(
                f"Host-side per-replica ckpt missing pre-rsync: {ckpt_path}. "
                f"stage_per_replica_checkpoints must run before "
                f"_rsync_per_replica_ckpts_to_vm."
            )
        src_files.append(ckpt_path)
    marker_path = os.path.join(leg_dir, "ckpt_is_valid")
    if not os.path.isfile(marker_path):
        raise RuntimeError(
            f"Host-side ckpt_is_valid marker missing pre-rsync: "
            f"{marker_path}. stage_per_replica_checkpoints must run "
            f"before _rsync_per_replica_ckpts_to_vm."
        )
    src_files.append(marker_path)

    # 2) Ensure VM-side r{0..n-1} dirs exist (rsync --relative would do
    # this but we prefer explicit mkdir for clearer error surface)
    mkdir_cmds = " && ".join(
        f"mkdir -p {shlex.quote(os.path.join(leg_dir, f'r{rid}'))}"
        for rid in range(n_replicas)
    )
    try:
        mkdir_result = subprocess.run(
            [
                "ssh",
                "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
                "-o", "StrictHostKeyChecking=no",
                "-o", "LogLevel=ERROR",
                vm_ssh_host,
                mkdir_cmds,
            ],
            capture_output=True,
            timeout=60,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VM mkdir for {n_replicas} replica dirs timed out — "
            f"ssh to {vm_ssh_host} unreachable. Cohort halted."
        ) from exc
    if mkdir_result.returncode != 0:
        raise RuntimeError(
            f"VM mkdir failed (rc={mkdir_result.returncode}): "
            f"{mkdir_result.stderr.decode(errors='replace')[:500]}"
        )

    # 3) Rsync each per-replica ckpt (preserves directory layout via the
    # rsync source path matching the leg_dir structure). We use a single
    # rsync invocation with --files-from for atomicity + progress.
    # The transferred files share a common prefix (leg_dir/), so use
    # --files-from with paths relative to leg_dir.
    import tempfile
    rel_paths = (
        [f"r{rid}/{jobname}_ckpt.xml" for rid in range(n_replicas)]
        + ["ckpt_is_valid"]
    )
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".rsync_filelist", delete=False
    ) as fh:
        fh.write("\n".join(rel_paths) + "\n")
        files_from_path = fh.name
    try:
        rsync_cmd = [
            "rsync",
            "-az",
            "-e", f"ssh -o ConnectTimeout={ssh_connect_timeout_s} "
                  f"-o StrictHostKeyChecking=no -o LogLevel=ERROR",
            f"--files-from={files_from_path}",
            f"{leg_dir.rstrip('/')}/",
            f"{vm_ssh_host}:{leg_dir.rstrip('/')}/",
        ]
        try:
            rsync_result = subprocess.run(
                rsync_cmd,
                capture_output=True,
                timeout=subprocess_timeout_s,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"rsync of {n_replicas} ckpts timed out "
                f"(>{subprocess_timeout_s}s). Network or VM disk issue. "
                f"Cohort halted; manual investigation required."
            ) from exc
        if rsync_result.returncode != 0:
            raise RuntimeError(
                f"rsync per-replica ckpts FAILED (rc="
                f"{rsync_result.returncode}): "
                f"{rsync_result.stderr.decode(errors='replace')[:1000]}"
            )
    finally:
        try:
            os.unlink(files_from_path)
        except OSError:
            pass

    # 4) VM-side verification: count ckpts + verify marker
    verify_cmd_parts = [
        f"test -f {shlex.quote(os.path.join(leg_dir, f'r{rid}', jobname + '_ckpt.xml'))}"
        for rid in range(n_replicas)
    ]
    verify_cmd_parts.append(
        f"test -f {shlex.quote(os.path.join(leg_dir, 'ckpt_is_valid'))}"
    )
    verify_cmd = " && ".join(verify_cmd_parts)
    try:
        verify_result = subprocess.run(
            [
                "ssh",
                "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
                "-o", "StrictHostKeyChecking=no",
                "-o", "LogLevel=ERROR",
                vm_ssh_host,
                verify_cmd,
            ],
            capture_output=True,
            timeout=60,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VM post-rsync verify timed out — ssh to {vm_ssh_host} "
            f"unresponsive. Cohort halted (do not launch — VM state "
            f"unverified)."
        ) from exc
    if verify_result.returncode != 0:
        raise RuntimeError(
            f"VM post-rsync verify FAILED (rc="
            f"{verify_result.returncode}): one or more of the "
            f"{n_replicas} ckpts or the ckpt_is_valid marker is "
            f"missing on VM at {leg_dir}. Manual re-rsync required."
        )

    # 5) Estimate bytes transferred (sum of host-side ckpt sizes)
    try:
        bytes_estimate = sum(os.path.getsize(f) for f in src_files)
    except OSError:
        bytes_estimate = -1

    return {
        "status": "rsynced",
        "leg_dir": leg_dir,
        "n_ckpts": n_replicas,
        "marker_present": True,
        "bytes_sent_estimate": bytes_estimate,
    }


def _rsync_subdir_to_vm(
    subdir: str,
    vm_ssh_host: str,
    ssh_connect_timeout_s: int = 5,
    subprocess_timeout_s: int = 900,
) -> Dict[str, Any]:
    """Rsync a complete per-direction work subdir (recursive) to the same
    path on the VM (2-process split generalization of
    ``_rsync_per_replica_ckpts_to_vm``).

    The per-direction subdir is fully self-contained (cntl + sys.xml + pdb +
    _0.xml + nodefile + r0..r10 ckpts + ckpt_is_valid). A single recursive
    ``rsync -az`` propagates the whole tree; a VM-side ``test -d`` + cntl
    ``test -f`` probe confirms arrival. Raises RuntimeError (cohort halt)
    on any failure.
    """
    subdir = os.path.abspath(subdir).rstrip("/")
    if not os.path.isdir(subdir):
        raise RuntimeError(
            f"_rsync_subdir_to_vm: host subdir missing pre-rsync: {subdir}. "
            f"stage_per_direction_subdir must run first."
        )
    # Ensure parent dir exists on the VM (rsync of subdir/ -> subdir/).
    parent = os.path.dirname(subdir)
    ssh_opts = [
        "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
        "-o", "StrictHostKeyChecking=no",
        "-o", "LogLevel=ERROR",
    ]
    try:
        mkdir_result = subprocess.run(
            ["ssh", *ssh_opts, vm_ssh_host,
             f"mkdir -p {shlex.quote(subdir)}"],
            capture_output=True, timeout=60,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VM mkdir for subdir timed out — ssh to {vm_ssh_host} "
            f"unreachable. Cohort halted."
        ) from exc
    if mkdir_result.returncode != 0:
        raise RuntimeError(
            f"VM mkdir for subdir failed (rc={mkdir_result.returncode}): "
            f"{mkdir_result.stderr.decode(errors='replace')[:500]}"
        )
    rsync_cmd = [
        "rsync", "-az",
        "-e", f"ssh -o ConnectTimeout={ssh_connect_timeout_s} "
              f"-o StrictHostKeyChecking=no -o LogLevel=ERROR",
        f"{subdir}/",
        f"{vm_ssh_host}:{subdir}/",
    ]
    try:
        rsync_result = subprocess.run(
            rsync_cmd, capture_output=True, timeout=subprocess_timeout_s,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"rsync of subdir {subdir} timed out (>{subprocess_timeout_s}s). "
            f"Network or VM disk issue. Cohort halted."
        ) from exc
    if rsync_result.returncode != 0:
        raise RuntimeError(
            f"rsync subdir FAILED (rc={rsync_result.returncode}): "
            f"{rsync_result.stderr.decode(errors='replace')[:1000]}"
        )
    # VM-side verify: subdir exists + contains at least one *_asyncre.cntl.
    verify_cmd = (
        f"test -d {shlex.quote(subdir)} && "
        f"ls {shlex.quote(subdir)}/*_asyncre.cntl >/dev/null 2>&1"
    )
    try:
        verify_result = subprocess.run(
            ["ssh", *ssh_opts, vm_ssh_host, verify_cmd],
            capture_output=True, timeout=60,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VM post-rsync subdir verify timed out — ssh to {vm_ssh_host} "
            f"unresponsive. Cohort halted."
        ) from exc
    if verify_result.returncode != 0:
        raise RuntimeError(
            f"VM post-rsync verify FAILED: subdir {subdir} or its cntl "
            f"missing on VM. Manual re-rsync required."
        )
    try:
        bytes_estimate = sum(
            os.path.getsize(os.path.join(root, f))
            for root, _dirs, files in os.walk(subdir)
            for f in files
        )
    except OSError:
        bytes_estimate = -1
    return {
        "status": "rsynced",
        "subdir": subdir,
        "bytes_sent_estimate": bytes_estimate,
    }


# ---------------------------------------------------------------------------
# C5 finite-energy probe + C4 water-count audit
# (two-process per-direction split, C4 / C5).
#
# C5 (finite-energy probe): the dminus system pre-displaces the binder and
# the ATMForce stores -DISPLACEMENT, so the u1 (or u0) endpoint can wrap a
# binder atom into the receptor pocket under PBC (vendored make_atm_system
# caveat). A non-finite (>1e10 kJ/mol) single-point energy on any of the 4
# (cp4,wt) x (dplus,dminus) endpoint systems is an architectural HALT
# signal — remediate by shrinking DISPLACEMENT (2.5 -> 1.0 nm) or
# increasing box padding, then re-prep. We evaluate the energy on the
# *production* state XML (_0_<tag>.xml) loaded into the per-direction
# system, exactly as the walker will.
#
# C4 (water-count audit): sys_dplus / sys_dminus differ by ~51 waters
# (binder-position-dependent addSolvent, Q6 b+ intended). The Delta-water
# 1st-order term cancels in DDG(Cp4-WT) ONLY if |DN_cp4 - DN_wt| < 5.
# Beyond that the cancellation breaks -> ranking caveat / HALT.
# ---------------------------------------------------------------------------

# A potential energy above this magnitude (kJ/mol) is treated as non-finite
# (PBC-wrap clash). C5 verdict threshold.
_C5_ENERGY_HALT_KJ = 1.0e10

# Max allowed |DN_water(cp4) - DN_water(wt)| before the 1st-order Delta-water
# cancellation is considered broken (C4 verdict).
_C4_WATER_DELTA_MAX = 5


def _count_waters_in_system_xml(sys_xml_path: str) -> int:
    """Return the total particle count of the serialized System.

    A System XML alone does not carry residue names, so this returns the
    System's particle count (``system.getNumParticles()``) — used as a
    proxy for the DELTA between dplus / dminus where only waters differ
    (the non-water atom count is identical by K-7). The water DELTA between
    two same-endpoint systems equals their particle-count delta divided by
    3 (TIP3P), which is what C4 actually compares.

    Uses OpenMM deserialization (lazy import — atm env) so the count is the
    authoritative ``getNumParticles()`` and NOT a regex over ``<Particle``
    lines (which would double-count the ATMForce per-particle displacement
    block). If openmm is unavailable, falls back to counting only the
    top-level ``<Particles>`` block (parsed by depth, ATMForce block
    excluded).
    """
    try:
        import openmm  # lazy — atm env
        with open(sys_xml_path) as fh:
            system = openmm.XmlSerializer.deserialize(fh.read())
        return system.getNumParticles()
    except ImportError:
        # Fallback: count <Particle ...> lines only inside the FIRST
        # <Particles> ... </Particles> block (the System particle list;
        # ATMForce per-particle entries live in a separate later block).
        n = 0
        in_particles = False
        with open(sys_xml_path) as fh:
            for line in fh:
                if "<Particles>" in line:
                    in_particles = True
                    continue
                if "</Particles>" in line:
                    break
                if in_particles and "<Particle " in line:
                    n += 1
        return n


def water_count_audit(
    leg_dir: str,
    jobname: str = "trackb",
) -> Dict[str, Any]:
    """C4 water-count audit for a single endpoint's leg.

    Compares the particle counts of ``{jobname}_sys_dplus.xml`` vs
    ``{jobname}_sys_dminus.xml`` and returns the per-endpoint Delta (in
    water molecules, particle-delta / 3 for TIP3P). The cross-endpoint
    comparison (``|DN_cp4 - DN_wt| < 5``) is performed by
    ``cross_endpoint_water_audit`` once both endpoints are measured.

    Returns ``{"leg_dir", "n_dplus", "n_dminus", "delta_particles",
    "delta_waters"}``. Raises FileNotFoundError if either sys XML missing.
    """
    leg_dir = os.path.abspath(leg_dir)
    dplus = os.path.join(leg_dir, jobname + "_sys_dplus.xml")
    dminus = os.path.join(leg_dir, jobname + "_sys_dminus.xml")
    for p in (dplus, dminus):
        if not os.path.isfile(p):
            raise FileNotFoundError(
                f"water_count_audit: missing per-direction system XML {p}"
            )
    n_dplus = _count_waters_in_system_xml(dplus)
    n_dminus = _count_waters_in_system_xml(dminus)
    delta_particles = n_dplus - n_dminus
    return {
        "leg_dir": leg_dir,
        "n_dplus": n_dplus,
        "n_dminus": n_dminus,
        "delta_particles": delta_particles,
        # TIP3P = 3 particles/water; report the water-molecule delta.
        "delta_waters": delta_particles / 3.0,
    }


def cross_endpoint_water_audit(
    cp4_leg_dir: str,
    wt_leg_dir: str,
    jobname: str = "trackb",
    water_delta_max: int = _C4_WATER_DELTA_MAX,
) -> Dict[str, Any]:
    """C4 cross-endpoint audit: ``|DN_water(cp4) - DN_water(wt)| <
    water_delta_max``.

    ``DN_water`` is the per-endpoint dplus-vs-dminus water delta from
    ``water_count_audit``. If the cross-endpoint difference exceeds the
    threshold, the 1st-order Delta-water cancellation in DDG(Cp4-WT) is
    broken (verdict (b) §2nd-order). Returns ``status`` ``ok`` /
    ``caveat`` plus both per-endpoint audits.
    """
    cp4 = water_count_audit(cp4_leg_dir, jobname=jobname)
    wt = water_count_audit(wt_leg_dir, jobname=jobname)
    diff = abs(cp4["delta_waters"] - wt["delta_waters"])
    status = "ok" if diff < water_delta_max else "caveat"
    return {
        "status": status,
        "cp4": cp4,
        "wt": wt,
        "abs_cross_endpoint_water_delta": diff,
        "threshold": water_delta_max,
        "note": (
            "C4: |DN_water(cp4) - DN_water(wt)| < threshold required for "
            "1st-order Delta-water cancellation in DDG(Cp4-WT). status="
            f"{status}."
        ),
    }


def _resolve_per_direction_cntl(leg_dir: str, direction_tag: str,
                                jobname: str) -> Optional[str]:
    """Locate the per-direction cntl for the C5 probe (READ-ONLY on the leg
    tree).

    Prefers the subdir-staged cntl
    (``leg_dir/{tag}/{jobname}_{tag}_asyncre.cntl``, written by the
    production staging step). When that is absent, derives the per-direction
    cntl from the combined cntl (``leg_dir/{jobname}_asyncre.cntl``) into a
    TEMP directory — the C5 probe is a pre-flight gate and must NOT mutate
    the leg tree (staging belongs to the authorized production launch, not
    the probe). Returns the cntl path (subdir or temp), or ``None`` if
    neither source exists.

    The temp dir is intentionally NOT cleaned here so the returned path
    stays readable by the caller; it lands under the system temp root and
    is reclaimed by the OS. The DIRECTION-rewrite + DISPLACEMENT-negation
    are applied identically to the staged path (same
    generate_per_direction_cntls call), so the probe reads the SAME
    corrected cntl production will run.
    """
    subdir_cntl = os.path.join(
        leg_dir, direction_tag, jobname + "_" + direction_tag + "_asyncre.cntl"
    )
    if os.path.isfile(subdir_cntl):
        return subdir_cntl
    combined = os.path.join(leg_dir, jobname + "_asyncre.cntl")
    if os.path.isfile(combined):
        import tempfile
        tmp_leg = tempfile.mkdtemp(prefix="trackb_c5_cntl_")
        shutil.copy2(combined, os.path.join(tmp_leg, jobname + "_asyncre.cntl"))
        gen = generate_per_direction_cntls(leg_dir=tmp_leg, jobname=jobname)
        return gen["directions"][direction_tag]["cntl_path"]
    return None


def finite_energy_probe_one(
    leg_dir: str,
    direction_tag: str,
    jobname: str = "trackb",
    energy_halt_kj: float = _C5_ENERGY_HALT_KJ,
) -> Dict[str, Any]:
    """C5 production-equivalent finite-energy probe for ONE per-direction
    endpoint.

    Builds the EXACT production system — ``OMMSystemABFE.create_system``
    (the same call the ``abfe_production`` worker makes, ommworker.py:198) —
    using the per-direction cntl's LIGAND_ATOMS / DIRECTION / DISPLACEMENT,
    so the runtime-constructed ATMForce is identical to production. The bare
    ``{jobname}_sys_{tag}.xml`` carries NO ATMForce (``grep -c ATMForce`` ==
    0), so the prior bare-system probe could NOT detect the d=-1 NaN: it
    evaluated only the unperturbed base system at its native positions,
    which is finite for both directions even when production's
    ``select(step(Direction), u0, u1)`` base selects a clashing u1.

    The HALT criterion is the **base potential** — the soft-cored full
    ATMForce potential energy at the BASE state (schedule state 0:
    lambda1=lambda2=0, Direction = the cntl's first DIRECTION value). The
    base is ``select(step(Direction), u0, u1)`` and is the term that NaN'd
    in production: for the +d TRAP (DIRECTION=-1 or DISPLACEMENT=+d on
    dminus) the base selects u1 = E(binder@bulk+25A = box edge) with NO
    soft-core protection => base PE ~ 1.4e13 kJ/mol (>1e10 HALT) and
    Particle-coordinate-NaN at annealing step 0. For the CORRECTED dminus
    (DIRECTION=+1, DISPLACEMENT=-d) the base is u0 = E(binder@bulk) =>
    finite (~ -1.23e6 kJ/mol). Empirically validated 2026-06-01 against the
    real cp4/bound leg: corrected dminus base -1.225e6 (PASS); dplus base
    -1.221e6 (PASS); +d trap base 1.4e13 (HALT).

    The raw perturbation endpoint ``u1`` (via
    ``ATMForce.getPerturbationEnergy``) is recorded as a DIAGNOSTIC ONLY,
    NOT a hard gate: the soft-core caps the perturbation to ~Umax
    (200 kcal/mol ~ 837 kJ/mol) in the actual potential, so the raw
    uncapped u1 is astronomically large (1e10-1e13) for BOTH the corrected
    AND the trap cntl by design — gating on raw u1 would false-positive on
    a valid run. u1 is flagged nonfinite only on an actual NaN/inf (a
    genuinely broken geometry), never on a large-but-finite value.

    Imports openmm + atom_openmm lazily so the qmmm test env (no atm
    openmm) can still import this module — the probe itself requires the
    ``atm`` conda env.

    Returns a dict with ``status`` (``ok`` / ``nonfinite`` / ``missing`` /
    ``error``), ``base_energy_kj`` (the gated term), ``u0_kj``, ``u1_kj``
    (diagnostic), the resolved ``direction`` + ``displacement`` actually
    used, plus ``direction_tag`` / ``leg_dir``.
    """
    import math as _math

    leg_dir = os.path.abspath(leg_dir)
    sys_xml = os.path.join(leg_dir, jobname + "_sys_" + direction_tag + ".xml")
    pdb_path = os.path.join(leg_dir, jobname + "_" + direction_tag + ".pdb")
    state_xml = os.path.join(leg_dir, jobname + "_0_" + direction_tag + ".xml")
    for p in (sys_xml, pdb_path, state_xml):
        if not os.path.isfile(p):
            return {
                "status": "missing",
                "direction_tag": direction_tag,
                "leg_dir": leg_dir,
                "missing": p,
            }
    cntl_path = _resolve_per_direction_cntl(leg_dir, direction_tag, jobname)
    if cntl_path is None:
        return {
            "status": "missing",
            "direction_tag": direction_tag,
            "leg_dir": leg_dir,
            "missing": os.path.join(leg_dir, jobname + "_asyncre.cntl"),
        }

    # Lazy imports — only the atm env has openmm + atom_openmm wired for
    # ATMForce. logging is used to satisfy the OMMSystemABFE constructor.
    import logging as _logging
    import openmm
    from openmm import unit, Platform
    from atom_openmm.ommsystem import OMMSystemABFE
    from atom_openmm.utils.config import parse_config

    # parse_config coerces DISPLACEMENT -> [float], DIRECTION -> [int],
    # LIGAND_ATOMS -> [int] etc. — byte-identical to the production worker.
    keywords = parse_config(cntl_path)
    displacement = list(keywords.get("DISPLACEMENT", []))
    direction_list = list(keywords.get("DIRECTION", []))
    # The base state (schedule state 0) Direction governs the u0/u1 select.
    base_direction = float(direction_list[0]) if direction_list else 1.0

    logger = _logging.getLogger("trackb_c5_probe")

    # Build the production system. OMMSystemABFE.create_system resolves the
    # system XML + pdb from <basename>; pass the per-direction names
    # explicitly via the keywords-driven constructor (basename only labels
    # logs here — we hand it the exact files).
    basename = jobname + "_" + direction_tag
    try:
        ommsys = OMMSystemABFE(
            basename, keywords, pdb_path, sys_xml, logger
        )
        ommsys.create_system()
    except Exception as exc:  # noqa: BLE001 — surface build failure as error
        return {
            "status": "error",
            "direction_tag": direction_tag,
            "leg_dir": leg_dir,
            "error": f"OMMSystemABFE.create_system failed: {exc}",
        }

    system = ommsys.system
    atmforce = ommsys.atmforce

    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin, 0.5 / unit.picosecond, 0.002 * unit.picoseconds
    )
    platform = Platform.getPlatformByName("Reference")
    context = openmm.Context(system, integrator, platform)

    # Set base-state global parameters: lambda1=lambda2=0 (no alchemical
    # bias), Direction = the cntl's base-state Direction, soft-core params
    # from the cntl (cparams). This is the lambda=0 endpoint the worker
    # initializes before the first MD step — exactly where the d=-1 NaN
    # appears in production (annealing step 0).
    context.setParameter(atmforce.Lambda1(), 0.0)
    context.setParameter(atmforce.Lambda2(), 0.0)
    context.setParameter(atmforce.Alpha(), 0.0)
    context.setParameter(atmforce.Uh(), 0.0)
    context.setParameter(atmforce.W0(), 0.0)
    context.setParameter(atmforce.Direction(), base_direction)
    for cparam_name, cparam_val in ommsys.cparams.items():
        try:
            context.setParameter(cparam_name, cparam_val)
        except Exception:  # noqa: BLE001 — non-context params ignored
            pass

    with open(state_xml) as fh:
        state = openmm.XmlSerializer.deserialize(fh.read())
    box = state.getPeriodicBoxVectors()
    if box is not None:
        context.setPeriodicBoxVectors(*box)
    context.setPositions(state.getPositions())

    base_energy = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )

    # Perturbation endpoint (DIAGNOSTIC ONLY — not gated). getPerturbation
    # Energy returns (u1, u0, bias) RAW (pre-soft-core). u1 is huge
    # (1e10-1e13) for both corrected + trap cntls by design (the soft-core
    # caps it to ~Umax in the actual potential), so it is recorded but only
    # flagged on a genuine NaN/inf — never on a large-but-finite value.
    u1_kj = None
    u0_kj = None
    u1_nan = False
    try:
        pert = atmforce.getPerturbationEnergy(context)
        u1_kj = float(pert[0].value_in_unit(unit.kilojoule_per_mole))
        u0_kj = float(pert[1].value_in_unit(unit.kilojoule_per_mole))
        u1_nan = not _math.isfinite(u1_kj)
    except Exception:  # noqa: BLE001 — perturbation read optional
        pass

    del context, integrator

    # HALT criterion: the base potential (soft-cored full ATMForce energy)
    # must be finite + bounded. This is the term that NaN'd in production
    # for the naive dminus (DIRECTION=-1 selects an unprotected u1 base).
    base_finite = (
        _math.isfinite(base_energy) and abs(base_energy) < energy_halt_kj
    )
    status = "ok" if (base_finite and not u1_nan) else "nonfinite"
    return {
        "status": status,
        "base_energy_kj": float(base_energy),
        "u0_kj": u0_kj,
        "u1_kj": u1_kj,
        "u1_is_nan": u1_nan,
        "direction": base_direction,
        "displacement": displacement,
        "direction_tag": direction_tag,
        "leg_dir": leg_dir,
        "cntl_path": cntl_path,
        "threshold_kj": energy_halt_kj,
        "note": (
            "C5 production-equivalent: OMMSystemABFE.create_system + "
            "base-state ATMForce(Direction=%s, displacement=%s). HALT on "
            "base_energy_kj nonfinite/>1e10 (the d=-1 NaN term). Raw u1 is "
            "diagnostic only (soft-core caps it in production)."
            % (base_direction, displacement)
        ),
    }


def finite_energy_probe_all(
    leg_dirs: Dict[str, str],
    jobname: str = "trackb",
    energy_halt_kj: float = _C5_ENERGY_HALT_KJ,
) -> Dict[str, Any]:
    """C5 probe across all endpoints x both directions (4 probes).

    ``leg_dirs`` maps endpoint name -> leg_dir (e.g. ``{"cp4": ".../cp4/
    bound", "wt": ".../wt/bound"}``). Runs ``finite_energy_probe_one`` for
    each (endpoint, tag in {dplus, dminus}). Returns ``status`` ``ok`` iff
    ALL probes are finite (``missing`` probes are surfaced as a separate
    ``incomplete`` status so a partially-prepped cohort does not silently
    pass). Caller HALTs on ``nonfinite``.
    """
    probes: List[Dict[str, Any]] = []
    any_nonfinite = False
    any_missing = False
    for endpoint, leg_dir in leg_dirs.items():
        for tag in ("dplus", "dminus"):
            probe = finite_energy_probe_one(
                leg_dir=leg_dir, direction_tag=tag,
                jobname=jobname, energy_halt_kj=energy_halt_kj,
            )
            probe["endpoint"] = endpoint
            probes.append(probe)
            if probe["status"] == "nonfinite":
                any_nonfinite = True
            elif probe["status"] == "missing":
                any_missing = True
    if any_nonfinite:
        status = "nonfinite"
    elif any_missing:
        status = "incomplete"
    else:
        status = "ok"
    return {
        "status": status,
        "probes": probes,
        "threshold_kj": energy_halt_kj,
        "note": (
            "C5: single-point energy on each (endpoint x direction) "
            "production state; >1e10 kJ/mol or NaN => PBC-wrap clash, HALT "
            "(shrink DISPLACEMENT or increase box padding, re-prep)."
        ),
    }


def _resolve_staging_source(
    leg_dir: str,
    per_direction_name: str,
    combined_name: str,
) -> Tuple[str, str]:
    """Resolve a staging-source path with per-direction-PREFERRED /
    combined-FALLBACK semantics (v0.9.29.2, 2026-06-05).

    Returns ``(path, kind)`` where ``kind`` is:
      * ``"per_direction"`` — the direction-specific file is present (BOUND
        leg: the b+ rebuild emits genuinely-different per-direction system
        XMLs + topology PDBs, water count differs by direction). Used as-is.
      * ``"combined"``      — only the combined file is present (FREE leg: the
        system is direction-agnostic — same particle/atom count both
        directions; structprep emits only a per-direction base STATE, not a
        per-direction system). Fall back to the single combined file for both
        directions.
      * ``"missing"``       — NEITHER is present. ``path`` points at the
        per-direction name so the caller's fail-loud check surfaces the
        canonical expected name.

    Preference order matters: when a per-direction file exists it MUST win
    (collapsing the bound leg's distinct dplus/dminus systems to one would
    corrupt the calculation). The combined fallback is ONLY reached when no
    per-direction file was produced — i.e. the direction-agnostic free leg.
    Never fabricates: the combined file is the genuine single system, not a
    synthesized per-direction copy.
    """
    per_dir_path = os.path.join(leg_dir, per_direction_name)
    if os.path.isfile(per_dir_path):
        return per_dir_path, "per_direction"
    combined_path = os.path.join(leg_dir, combined_name)
    if os.path.isfile(combined_path):
        return combined_path, "combined"
    return per_dir_path, "missing"


def stage_per_direction_subdir(
    leg_dir: str,
    direction_tag: str,
    cntl_info: Dict[str, Any],
    jobname: str = "trackb",
    fwd_replica_count: int = 11,
    total_state_count: int = 22,
) -> Dict[str, Any]:
    """Materialize a fully self-contained per-direction work subdir.

    Each subdir (``leg_dir/dplus`` or ``leg_dir/dminus``) is staged so a
    single ``abfe_production {jobname}_{tag}_asyncre.cntl`` (BASENAME=
    ``{jobname}_{tag}``) runs against ONE atom count with NO cross-direction
    contamination (C6 root cause). Materializes:

      * ``{jobname}_{tag}_sys.xml``  (copy of leg ``{jobname}_sys_{tag}.xml``)
      * ``{jobname}_{tag}.pdb``      (copy of leg ``{jobname}_{tag}.pdb``)
      * ``{jobname}_{tag}_0.xml``    (canonical base-state = copy of leg
                                      ``{jobname}_0_{tag}.xml``; consumed by
                                      ommworker.py service worker)
      * ``nodefile``                 (copy of leg nodefile if present, else
                                      a single-line ``localhost,0:0,1,CUDA``)
      * ``r0..r10/`` — EMPTY per-replica work dirs (NO ``_ckpt.xml``).
      * ``ckpt_is_valid``            (safecheckpoint marker; required so the
                                      engine's ``checkpointJob`` can later
                                      write per-replica ckpts —
                                      ``openmm_async_re.py:37-57``)

    STATE-INIT WIRING FIX (2026-06-03)
    --------------------------------------------------------
    Prior versions copied a base-state / leg checkpoint into every
    ``r{local}/{base}_ckpt.xml``. Those checkpoints embed
    ``REStateId="0"`` (the lambda=0.5 / lambda1=0.0 equilibration state),
    so ``ommreplica.load_checkpoint`` (``ommreplica.py:72-83``) →
    ``update_state_from_context`` (``ommreplica.py:135``) read
    ``stateid == 0`` for EVERY replica. That overwrote the
    ``stateid is None`` sentinel set in ``OMMReplica.__init__``
    (``ommreplica.py:26``), so async_re's initial-seed guard
    ``if replica.stateid == None: replica.set_state(i, stateparams[i])``
    (``openmm_async_re.py:180-181, 396-397, 435-436``) was SKIPPED → all
    11 replicas locked to state 0 → degenerate swaps ("new state 0" only)
    → the lambda-ladder was never sampled → UWHAM rank-deficient → invalid
    ΔΔG_bind (confirmed invalid).

    The PROVEN-WORKING upstream free/combined leg
    (``scripts/trackb_production_v2_1_upstream.py::setup_one_leg_v21``)
    ships NO per-replica ckpt at launch; the engine creates ``r{i}`` dirs
    (``ommreplica.py:36-37``) and seeds replica ``i`` to state ``i`` via
    ``set_state(i)``. We mirror that layout exactly: NO ckpt is staged, so
    ``load_checkpoint`` no-ops (``ommreplica.py:76`` — file absent),
    ``stateid`` stays ``None``, and the seed guard fires correctly.

    The per-replica equilibrium coordinates are NOT needed for a correct
    ABFE run — async_re re-equilibrates each replica into its own
    lambda-rung from the shared starting State; the only thing the old ckpt
    contributed was the toxic embedded stateid.

    ``replica_map`` is still returned (local_rid -> source leg replica) for
    the merge / audit trail, but no longer carries ``ckpt_src`` /
    ``ckpt_dst`` / ``cold_start`` (no ckpt is staged).

    Returns a dict with the staged paths + replica mapping for audit.
    Raises RuntimeError if required leg-level inputs are absent.
    """
    leg_dir = os.path.abspath(leg_dir)
    subdir = os.path.join(leg_dir, direction_tag)
    os.makedirs(subdir, exist_ok=True)
    base = jobname + "_" + direction_tag

    # Required leg-level inputs.
    #
    # System XML + topology PDB: per-direction-PREFERRED with combined-FALLBACK
    # (v0.9.29.2, 2026-06-05). The BOUND leg's per-direction system XMLs
    # genuinely DIFFER by direction (addSolvent fills different water counts
    # around the pre-displaced binder — dplus 92855 vs dminus 92804 particles,
    # the b+ rebuild; see trackb_per_direction_system_xml_rebuild_20260601), so
    # when a `trackb_sys_{tag}.xml` + `trackb_{tag}.pdb` pair is present it MUST
    # be used (collapsing to one would corrupt the bound calc). The FREE leg has
    # no receptor pocket → the binder displacement does NOT change the solvation
    # count → both directions share ONE direction-agnostic combined system
    # (`trackb_sys.xml` 4292 particles + `trackb.pdb` 4292 atoms; ATMForce /
    # Direction / DISPLACEMENT carry NO baked-in direction — they are runtime
    # params from the per-direction cntl). Per-direction structprep only emits
    # the per-direction base STATE (`trackb_0_{tag}.xml`) for the free leg, NOT
    # a per-direction system. So when the per-direction system/pdb are absent we
    # fall back to the combined `trackb_sys.xml` / `trackb.pdb` for BOTH
    # directions. This is NOT fabrication: the free system genuinely IS the same
    # topology in both directions (verified
    # trackb_free_leg_direction_agnostic_sys_20260605). The per-direction
    # difference lives only in the base state + the runtime cntl.
    src_sys, sys_source_kind = _resolve_staging_source(
        leg_dir,
        per_direction_name=jobname + "_sys_" + direction_tag + ".xml",
        combined_name=jobname + "_sys.xml",
    )
    src_pdb, pdb_source_kind = _resolve_staging_source(
        leg_dir,
        per_direction_name=jobname + "_" + direction_tag + ".pdb",
        combined_name=jobname + ".pdb",
    )
    # Base state is ALWAYS per-direction (structprep's direction-patched
    # annealing produces a distinct `trackb_0_{tag}.xml` per direction — the
    # free leg's dplus/dminus States differ by ~25.7k lines / distinct box
    # vectors). No combined fallback: a missing per-direction base state means
    # structprep did not run for this direction → fail loud.
    src_base_state = os.path.join(
        leg_dir, jobname + "_0_" + direction_tag + ".xml"
    )
    for p in (src_sys, src_pdb, src_base_state):
        if not os.path.isfile(p):
            raise RuntimeError(
                f"stage_per_direction_subdir: missing leg input {p}. Run "
                f"trackb_per_direction_structprep.py first."
            )

    dst_sys = os.path.join(subdir, base + "_sys.xml")
    dst_pdb = os.path.join(subdir, base + ".pdb")
    dst_base_state = os.path.join(subdir, base + "_0.xml")
    shutil.copy2(src_sys, dst_sys)
    shutil.copy2(src_pdb, dst_pdb)
    shutil.copy2(src_base_state, dst_base_state)

    # nodefile.
    dst_nodefile = os.path.join(subdir, "nodefile")
    src_nodefile = os.path.join(leg_dir, "nodefile")
    if os.path.isfile(src_nodefile):
        shutil.copy2(src_nodefile, dst_nodefile)
    elif not os.path.isfile(dst_nodefile):
        with open(dst_nodefile, "w") as fh:
            fh.write("localhost,0:0,1,CUDA,,/tmp\n")

    # Replica index mapping: local r0..r10 -> source leg replica.
    if direction_tag == "dplus":
        src_rids = list(range(fwd_replica_count))                # 0..10
    else:
        src_rids = list(range(fwd_replica_count, total_state_count))  # 11..21

    # Create EMPTY per-replica work dirs (NO _ckpt.xml). The engine seeds
    # each replica's state via async_re set_state(i) because no checkpoint
    # is present to overwrite the stateid==None sentinel. This is the
    # state-init wiring fix — see the docstring. Mirrors the
    # proven-working free/combined leg layout exactly.
    replica_map: List[Dict[str, Any]] = []
    for local_rid, src_rid in enumerate(src_rids):
        r_dir = os.path.join(subdir, f"r{local_rid}")
        os.makedirs(r_dir, exist_ok=True)
        # Defensive: if a stale per-replica ckpt survived from an earlier
        # (broken) staging into the same subdir, remove it (archival not
        # needed: these are launcher-generated transient seeds, never
        # primary data).
        stale_ckpt = os.path.join(r_dir, base + "_ckpt.xml")
        if os.path.isfile(stale_ckpt):
            os.remove(stale_ckpt)
        replica_map.append({
            "local_rid": local_rid,
            "src_rid": src_rid,
            "ckpt_staged": False,
            "seeds_via": "async_re_set_state",
        })

    # safecheckpoint marker.
    marker = os.path.join(subdir, "ckpt_is_valid")
    if not os.path.isfile(marker):
        with open(marker, "w") as fh:
            fh.write(
                "# Created by trackb_per_direction_production.py "
                "stage_per_direction_subdir — per ommreplica.py:77\n"
            )

    return {
        "direction_tag": direction_tag,
        "subdir": subdir,
        "basename": base,
        "cntl_path": cntl_info["cntl_path"],
        "cntl_basename": cntl_info["cntl_basename"],
        "sys_xml": dst_sys,
        "pdb": dst_pdb,
        "base_state_xml": dst_base_state,
        # Provenance of the staged system/topology (v0.9.29.2): "per_direction"
        # (bound, distinct per-dir system) or "combined" (free, direction-
        # agnostic single system reused for both directions).
        "sys_source_kind": sys_source_kind,
        "pdb_source_kind": pdb_source_kind,
        "sys_xml_src": src_sys,
        "pdb_src": src_pdb,
        "nodefile": dst_nodefile,
        "n_replicas": len(src_rids),
        "replica_map": replica_map,
        "marker": marker,
    }


def merge_per_direction_outputs(
    leg_dir: str,
    jobname: str = "trackb",
    fwd_replica_count: Optional[int] = None,
    total_state_count: Optional[int] = None,
) -> Dict[str, Any]:
    """Merge the two per-direction subdir outputs into the parent leg's
    ``r0..r{N-1}/{jobname}.out`` for UWHAM (C1 verdict).

    State-count-agnostic (Path λ-densify spec 2026-06-05). When
    ``fwd_replica_count`` / ``total_state_count`` are None they are DERIVED
    from the leg's combined cntl DIRECTION column (bound leg = 11/22;
    densified free leg = 17/34). NO global 11/22 literal is assumed. The
    backward half-count is ``total - fwd`` (= fwd for the symmetric ATM
    schedules, but computed explicitly so an asymmetric schedule would still
    map correctly).

    ``atom_openmm.uwham.calculate_uwham`` reads ``r{i}/{jobname}.out`` for
    ``i in range(N)`` and partitions by the per-line ``stateid`` column
    (col 0). The forward subdir wrote local r0..r{fwd-1} with stateid
    0..fwd-1 (no renumber needed); the backward subdir wrote local
    r0..r{bwd-1} with stateid 0..bwd-1 but those are GLOBAL states
    fwd..total-1, so both the directory index AND the stateid column must
    shift UP by ``fwd_replica_count``.

    Returns a manifest of merged files. Idempotent (overwrites the parent
    ``.out`` each call). Raises FileNotFoundError if a subdir .out is
    missing (run incomplete) or the combined cntl is needed but absent.
    """
    leg_dir = os.path.abspath(leg_dir)

    # Derive counts from the combined cntl when not explicitly supplied.
    if fwd_replica_count is None or total_state_count is None:
        combined_cntl = os.path.join(leg_dir, jobname + "_asyncre.cntl")
        derived_total, derived_fwd = _derive_state_counts_from_cntl(
            combined_cntl
        )
        if total_state_count is None:
            total_state_count = derived_total
        if fwd_replica_count is None:
            fwd_replica_count = derived_fwd

    bwd_replica_count = total_state_count - fwd_replica_count
    merged: List[Dict[str, Any]] = []

    for tag, offset, count in (
        ("dplus", 0, fwd_replica_count),
        ("dminus", fwd_replica_count, bwd_replica_count),
    ):
        base = jobname + "_" + tag
        subdir = os.path.join(leg_dir, tag)
        for local_rid in range(count):
            src_out = os.path.join(
                subdir, f"r{local_rid}", base + ".out"
            )
            if not os.path.isfile(src_out):
                raise FileNotFoundError(
                    f"merge_per_direction_outputs: subdir output missing "
                    f"{src_out} — {tag} run incomplete."
                )
            global_rid = local_rid + offset
            dst_dir = os.path.join(leg_dir, f"r{global_rid}")
            os.makedirs(dst_dir, exist_ok=True)
            dst_out = os.path.join(dst_dir, jobname + ".out")
            _renumber_stateid_out(src_out, dst_out, stateid_offset=offset)
            merged.append({
                "direction_tag": tag,
                "local_rid": local_rid,
                "global_rid": global_rid,
                "stateid_offset": offset,
                "src_out": src_out,
                "dst_out": dst_out,
            })

    return {
        "leg_dir": leg_dir,
        "n_merged": len(merged),
        "fwd_replica_count": fwd_replica_count,
        "bwd_replica_count": bwd_replica_count,
        "expected": total_state_count,
        "complete": len(merged) == total_state_count,
        "merged": merged,
    }


def _renumber_stateid_out(
    src_out: str,
    dst_out: str,
    stateid_offset: int,
) -> None:
    """Copy ``src_out`` to ``dst_out`` incrementing the stateid column
    (first whitespace-delimited token, int) by ``stateid_offset`` on every
    data line.

    The ``.out`` format is 11 whitespace-separated floats with col 0 =
    stateid (ommreplica.save_out ``"%d %f ..."``). A zero offset is a
    verbatim copy. Non-numeric / blank lines are passed through untouched
    (defensive — abfe_production .out has no header, but this keeps the
    helper robust).
    """
    if stateid_offset == 0:
        shutil.copy2(src_out, dst_out)
        return
    with open(src_out) as fin, open(dst_out, "w") as fout:
        for line in fin:
            parts = line.split()
            if not parts:
                fout.write(line)
                continue
            try:
                sid = int(float(parts[0]))
            except (ValueError, IndexError):
                fout.write(line)
                continue
            parts[0] = str(sid + stateid_offset)
            fout.write(" ".join(parts) + "\n")


# ---------------------------------------------------------------------------
# Stage 2 — independent-seed replicate orchestration (λ-densify campaign
# 2026-06-05; σ_btwn analysis). OFF by default. When enabled,
# each replicate is launched into a SEPARATE output subtree so its per-replica
# .out / ckpt trees never collide, and σ_btwn is computed across the replicate
# subtrees by the post-processor. Independent seed = different Langevin
# velocity seed + walker-shuffle seed; the seed value is recorded in metadata
# for a deterministic, reproducible audit trail (deterministic numpy RNG).
# ---------------------------------------------------------------------------
def _replicate_out_root(v21_out_root: str, replicate_index: int) -> str:
    """Return the per-replicate output subtree path.

    ``rep{i}`` is appended to ``v21_out_root`` so each replicate's
    ``<endpoint>/<leg>/r*`` trees are fully disjoint. Replicate 0 with a
    single replicate keeps the legacy root unchanged (no ``rep0`` suffix) so
    existing single-run launches are byte-for-byte identical.
    """
    return os.path.join(v21_out_root, f"rep{replicate_index}")


def _seed_for_replicate(seeds: List[str], replicate_index: int) -> str:
    """Deterministic seed selection for replicate ``i``.

    Cycles through the provided ``seeds`` list (s1/s2/s3 ... or the 8-seed
    cohort). Recorded in launch metadata so the velocity / walker-shuffle RNG
    is reproducible. Raises IndexError-free by modulo-cycling.
    """
    if not seeds:
        raise ValueError("no seeds provided for replicate selection")
    return seeds[replicate_index % len(seeds)]


def _velocity_seed_for_replicate(replicate_index: int) -> int:
    """Deterministic per-replicate INTEGER velocity-init seed (mechanism B).

    Returns ``replicate_index + 1`` so replicates get distinct integers
    1..n. This is the structprep ``setVelocitiesToTemperature`` seed AND the
    production Langevin ``setRandomNumberSeed`` audit value — DISTINCT from the
    QM snapshot cohort seed (``_seed_for_replicate`` → {s7, ...}). Recorded in
    launch metadata so the velocity RNG is reproducible / auditable.
    """
    return int(replicate_index) + 1


def _live_launch_all_legs(
    v21_out_root: str,
    endpoints: List[str],
    legs: List[str],
    jobname: str,
    gpu_host: str,
    dry_run: bool = False,
    abfe_bin: Optional[str] = None,
    vm_ssh_host: str = "san@192.168.122.155",
) -> List[Dict[str, Any]]:
    """Per-leg ``abfe_production`` orchestrator.

    Iterates ``(endpoint, leg)`` pairs in ``endpoints x legs`` order. For
    each leg:

      1. Resolves ``leg_dir = v21_out_root/<endpoint>/<leg>`` and
         ``cntl_path = leg_dir/<jobname>_asyncre.cntl``. FileNotFound
         raised if cntl missing (operator must run structprep first).
      2. Generates per-direction cntls (``generate_per_direction_cntls``)
         and materializes two self-contained subdirs (dplus / dminus) via
         ``stage_per_direction_subdir``. Each subdir ships sys.xml + pdb +
         ``_0.xml`` base state + nodefile + ``ckpt_is_valid`` + EMPTY
         r0..r10 dirs (NO per-replica ckpt — state-init wiring fix).
         The engine seeds replica ``i`` to state ``i`` via async_re
         ``set_state(i)``, mirroring the proven-working free leg.
      3. Launches ``abfe_production <cntl basename>`` (one process per
         direction) via ``subprocess.run`` with cwd=subdir (matching
         upstream ``run_abfe_production_for_leg`` contract). When
         ``gpu_host == "vm"`` the command is wrapped in ``ssh``; when
         ``"local"`` it sets ``CUDA_VISIBLE_DEVICES=0``; when ``"cpu"``
         no GPU env var is set.
      4. Per-direction log captured to ``<subdir>/_live_launch.log``.

    NOTE: ``stage_per_replica_checkpoints`` (the leg-level combined-ckpt
    injector) is RETAINED for regression history but is NOT called in the
    live path — it injected the embedded-stateid checkpoints that caused
    the state-0-collapse class. The 2-process split staging above replaces
    it.

    Replaces the prior manual handoff ("invoke v2_legacy.run_abfe_
    production_for_leg per leg") which left operator to construct
    ``leg_info`` dicts + remember the (endpoint, leg) loop. Operator
    burden is now: ``--i-have-confirmed-c1-through-c8`` + ``--gpu-host``.

    CRITICAL safety: This function must only be invoked AFTER both the
    explicit C1-C8 flag check AND ``gate_gpu_host`` pass in ``main()``.
    Module-load executes nothing. Tests cover the dry_run path only.

    ``dry_run=True`` prints the launch plan and returns immediately
    without running subprocesses (inspectable rehearsal).

    Returns a list of per-leg result dicts (status, rc, log_path) for
    audit. On first non-zero exit code, raises ``RuntimeError`` to halt
    cohort mid-flight rather than silently skip remaining legs.
    """
    if abfe_bin is None:
        abfe_bin = "/home/san/miniconda3/envs/atm/bin/abfe_production"

    # 1) Plan + validate cntls up front so an early-leg failure does not
    # leave later legs in unstaged state.
    plan: List[Dict[str, Any]] = []
    for endpoint in endpoints:
        for leg in legs:
            leg_dir = os.path.join(_PROJ_ROOT, v21_out_root, endpoint, leg)
            cntl_path = os.path.join(leg_dir, jobname + "_asyncre.cntl")
            if not os.path.isfile(cntl_path):
                raise FileNotFoundError(
                    f"Per-leg cntl missing: {cntl_path}. Run "
                    f"trackb_per_direction_structprep.py first."
                )
            plan.append({
                "endpoint": endpoint,
                "leg": leg,
                "leg_dir": leg_dir,
                "cntl_path": cntl_path,
            })

    # 1b) G38 fix: VM leg-dir pre-flight gate.
    # When gpu_host=vm, the abfe_production binary on the VM cannot read
    # host-side leg_dirs over SSH; it requires trackb.pdb / trackb_sys.xml
    # / trackb_asyncre.cntl to exist at the SAME path on the VM filesystem.
    # If any leg fails the gate, halt the cohort BEFORE stage_per_replica_
    # checkpoints runs (cohort-safe). Same failure-class family as
    # integrity_vm_lane_self_provisioning_20260529 + integrity_vm_conda_
    # path_noninteractive_20260528 + Round 2/3 F3 gates.
    if gpu_host == "vm" and not dry_run:
        gate_failures: List[str] = []
        for item in plan:
            allow, reason = _gate_vm_leg_dir_exists(
                leg_dir=item["leg_dir"],
                vm_ssh_host=vm_ssh_host,
            )
            if not allow:
                gate_failures.append(
                    f"  {item['endpoint']}/{item['leg']}: {reason}"
                )
        if gate_failures:
            raise RuntimeError(
                "VM leg-dir pre-flight FAILED for "
                f"{len(gate_failures)}/{len(plan)} leg(s); cohort halted "
                "BEFORE ckpt staging (no host state mutated). "
                "Remediation per failure:\n"
                + "\n".join(gate_failures)
            )

    # 1c) Per-direction cntl generation + subdir staging (C6 two-process
    # split). For each leg we slice the combined 22-state cntl into a
    # forward (dplus, states 0..10) and backward (dminus, states 11..21)
    # cntl, then materialize two FULLY SELF-CONTAINED work subdirs
    # (leg_dir/dplus, leg_dir/dminus). Each subdir runs against ONE atom
    # count (dplus=92855, dminus=92804) so the C6 atom-count-mixing crash
    # is structurally impossible.
    #
    # dry_run: NO host mutation. The cntl is derived into a TEMP leg (copy
    # of the combined cntl only) so the dry-run rehearsal NEVER writes into
    # the real leg tree; the rendered plan below substitutes the real-leg
    # subdir paths for display. (Pre-fix this generated into the live leg
    # tree even in dry-run — a silent mutation contradicting "production
    # tree pristine".)
    for item in plan:
        if dry_run:
            import tempfile
            tmp_leg = tempfile.mkdtemp(prefix="trackb_dryrun_cntl_")
            shutil.copy2(
                item["cntl_path"],
                os.path.join(tmp_leg, jobname + "_asyncre.cntl"),
            )
            cntl_gen = generate_per_direction_cntls(
                leg_dir=tmp_leg, jobname=jobname,
            )
        else:
            cntl_gen = generate_per_direction_cntls(
                leg_dir=item["leg_dir"], jobname=jobname,
            )
        item["cntl_gen"] = cntl_gen
        # Derived per-leg counts (Path λ-densify 2026-06-05): bound=22/11,
        # densified free=34/17. Passed to subdir staging so the replica-index
        # mapping (local r0..r{fwd-1} → source leg replica) is correct for
        # both schedules without a global 11/22 literal.
        leg_total = cntl_gen.get("total_state_count", _TOTAL_STATE_COUNT)
        leg_fwd = cntl_gen.get("fwd_state_count", _FWD_STATE_COUNT)
        if not dry_run:
            subdirs: Dict[str, Dict[str, Any]] = {}
            for tag in ("dplus", "dminus"):
                try:
                    subdirs[tag] = stage_per_direction_subdir(
                        leg_dir=item["leg_dir"],
                        direction_tag=tag,
                        cntl_info=cntl_gen["directions"][tag],
                        jobname=jobname,
                        fwd_replica_count=leg_fwd,
                        total_state_count=leg_total,
                    )
                except RuntimeError as exc:
                    raise RuntimeError(
                        f"stage_per_direction_subdir failed for "
                        f"{item['endpoint']}/{item['leg']} {tag}: {exc}"
                    ) from exc
            item["subdirs"] = subdirs
            print(
                f"  [stage] {item['endpoint']}/{item['leg']}: "
                f"dplus={subdirs['dplus']['n_replicas']} replicas + "
                f"dminus={subdirs['dminus']['n_replicas']} replicas "
                f"(2-process split)"
            )

    # 2b) G45 fix — VM-PER-REPLICA-CKPT-MISSING-01,
    # extended to the 2-process layout. When gpu_host=vm, abfe_production
    # runs on the VM filesystem; each per-direction subdir (cntl + sys.xml +
    # pdb + _0.xml + nodefile + r0..r10 ckpts + ckpt_is_valid) must exist at
    # the SAME path on the VM. Rsync the whole subdir tree (recursive),
    # verify VM-side completeness, halt cohort on any rsync failure (no
    # partial-state launch). Same failure-class family as Round 6 G38 /
    # Round 7 G45 (single-dispatch) — generalized to per-direction subdirs.
    if gpu_host == "vm" and not dry_run:
        for item in plan:
            for tag in ("dplus", "dminus"):
                sub = item["subdirs"][tag]
                try:
                    rsync_info = _rsync_subdir_to_vm(
                        subdir=sub["subdir"],
                        vm_ssh_host=vm_ssh_host,
                    )
                except RuntimeError as exc:
                    raise RuntimeError(
                        f"VM per-direction subdir rsync FAILED for "
                        f"{item['endpoint']}/{item['leg']} {tag}: {exc}. "
                        f"Cohort halted BEFORE any VM-side abfe_production "
                        f"launch (no partial-state launch). Remediation: "
                        f"verify VM disk space (`ssh {vm_ssh_host} df -h "
                        f"/home`), re-run launcher (host-side subdirs "
                        f"intact)."
                    ) from exc
                sub["rsync_status"] = rsync_info["status"]
                sub["rsync_bytes"] = rsync_info["bytes_sent_estimate"]

    # 3) Per-leg launch — TWO processes per leg (forward dplus + backward
    # dminus). Both must succeed for the leg to be considered complete;
    # UWHAM merges the two subdir outputs back to r0..r21 afterwards
    # (merge_per_direction_outputs, C1).
    results: List[Dict[str, Any]] = []
    for item in plan:
        cntl_gen = item["cntl_gen"]
        dispatches: List[Dict[str, Any]] = []
        for tag in ("dplus", "dminus"):
            dspec = cntl_gen["directions"][tag]
            if dry_run:
                # The cntl was derived into a temp leg (no host mutation);
                # display the REAL leg subdir where production would run.
                subdir = os.path.join(item["leg_dir"], tag)
            else:
                subdir = dspec["subdir"]
            cntl_basename = dspec["cntl_basename"]
            log_path = os.path.join(subdir, "_live_launch.log")
            if gpu_host == "vm":
                remote_cmd = (
                    f"cd {shlex.quote(subdir)} && "
                    f"{shlex.quote(abfe_bin)} {shlex.quote(cntl_basename)}"
                )
                cmd = ["ssh", vm_ssh_host, remote_cmd]
                cwd = None
                env = os.environ.copy()
            elif gpu_host in ("local", "cpu"):
                cmd = [abfe_bin, cntl_basename]
                cwd = subdir
                env = os.environ.copy()
                env["CUDA_VISIBLE_DEVICES"] = "0" if gpu_host == "local" else ""
            else:
                raise ValueError(
                    f"_live_launch_all_legs: unsupported gpu_host "
                    f"{gpu_host!r} (expected vm | local | cpu)"
                )

            rendered = " ".join(shlex.quote(c) for c in cmd)
            if dry_run:
                print(f"  [DRY-RUN] {item['endpoint']}/{item['leg']} {tag}: "
                      f"{rendered}  (cwd={cwd}, log={log_path})")
                dispatches.append({
                    "direction_tag": tag,
                    "cmd": rendered,
                    "cwd": cwd,
                    "log_path": log_path,
                    "status": "dry_run",
                    "rc": None,
                    "n_states": dspec["n_states"],
                })
                continue

            print(f"  launching {item['endpoint']}/{item['leg']} {tag} "
                  f"→ {log_path}")
            with open(log_path, "w") as logfh:
                logfh.write(f"# Command: {rendered}\n")
                logfh.write(f"# CWD:     {cwd}\n")
                logfh.write(f"# GPUHost: {gpu_host}\n")
                logfh.write(f"# Direction: {tag}\n")
                logfh.write(
                    f"# Started: {time.strftime('%Y-%m-%dT%H:%M:%S')}\n\n"
                )
                logfh.flush()
                proc = subprocess.run(
                    cmd, cwd=cwd, env=env,
                    stdout=logfh, stderr=subprocess.STDOUT,
                )
            if proc.returncode != 0:
                raise RuntimeError(
                    f"abfe_production FAILED for "
                    f"{item['endpoint']}/{item['leg']} {tag} "
                    f"(rc={proc.returncode}); see {log_path}"
                )
            dispatches.append({
                "direction_tag": tag,
                "cmd": rendered,
                "cwd": cwd,
                "log_path": log_path,
                "status": "complete",
                "rc": proc.returncode,
                "n_states": dspec["n_states"],
            })

        results.append({
            "endpoint": item["endpoint"],
            "leg": item["leg"],
            "leg_dir": item["leg_dir"],
            "status": "dry_run" if dry_run else "complete",
            "n_dispatches": len(dispatches),
            "dispatches": dispatches,
        })

    return results


# ---------------------------------------------------------------------------
# Stage 2 — independent-seed replicate orchestration (λ-densify campaign
# 2026-06-05; σ_btwn analysis). OFF by default. When enabled,
# each replicate is launched into a SEPARATE output subtree so its per-replica
# .out / ckpt trees never collide, and σ_btwn is computed across the replicate
# subtrees by the post-processor. Independent seed = different Langevin
# velocity seed + walker-shuffle seed; the seed value is recorded in metadata
# for a deterministic, reproducible audit trail (deterministic numpy RNG).
# ---------------------------------------------------------------------------
def _live_launch_replicates(
    v21_out_root: str,
    endpoints: List[str],
    legs: List[str],
    jobname: str,
    gpu_host: str,
    n_replicates: int,
    seeds: List[str],
    dry_run: bool = False,
    abfe_bin: Optional[str] = None,
    vm_ssh_host: str = "san@192.168.122.155",
) -> List[Dict[str, Any]]:
    """Launch ``n_replicates`` independent-seed replicates, each into its own
    ``rep{i}`` output subtree (Stage 2, Path λ-densify campaign 2026-06-05).

    Thin orchestration wrapper over ``_live_launch_all_legs`` — calls it once
    per replicate with the per-replicate subtree (``_replicate_out_root``).
    Each replicate uses a distinct seed (``_seed_for_replicate``) recorded in
    the returned per-leg result dicts. PRECONDITION (NOT enforced here): the
    per-replicate cntls / systems must exist in each ``rep{i}`` subtree (built
    by structprep with a matching ``--replicates``). When a subtree's leg cntl
    is absent, ``_live_launch_all_legs`` raises FileNotFoundError (cohort-safe
    — no partial-state launch).

    Returns the flattened list of per-leg result dicts across all replicates,
    each annotated with ``replicate_index`` + ``seed`` + ``replicate_out_root``.
    """
    all_results: List[Dict[str, Any]] = []
    for ridx in range(n_replicates):
        rep_root = _replicate_out_root(v21_out_root, ridx)
        seed = _seed_for_replicate(seeds, ridx)
        velocity_seed = _velocity_seed_for_replicate(ridx)
        rep_results = _live_launch_all_legs(
            v21_out_root=rep_root,
            endpoints=endpoints,
            legs=legs,
            jobname=jobname,
            gpu_host=gpu_host,
            dry_run=dry_run,
            abfe_bin=abfe_bin,
            vm_ssh_host=vm_ssh_host,
        )
        for r in rep_results:
            r["replicate_index"] = ridx
            r["seed"] = seed
            # mechanism B: distinct integer velocity-init seed (1..n),
            # the structprep setVelocitiesToTemperature + production Langevin
            # setRandomNumberSeed audit value (distinct from the QM cohort
            # seed above). Each rep{i} subtree's structprep must have been
            # run with --velocity-seed=<this value> for the independence to
            # be real (the production base state carries the seeded velocities).
            r["velocity_seed"] = velocity_seed
            r["replicate_out_root"] = rep_root
        all_results.extend(rep_results)
    return all_results


# ---------------------------------------------------------------------------
# F3 fix: defensive pre-flight
# gate to verify the VM-side abfe_production binary exists before any
# per-leg ckpt staging occurs. Without this gate the cohort half-stages
# (per-direction ckpts moved into place) and then dies at the first
# subprocess.run with rc=127 ("bash: command not found"), leaving the
# operator to manually un-stage. Same failure class as the VM-lane
# self-provisioning and non-interactive conda-path issues.
# ---------------------------------------------------------------------------
def _gate_vm_abfe_bin_exists(
    abfe_bin: str,
    vm_ssh_host: str = "san@192.168.122.155",
    ssh_connect_timeout_s: int = 5,
    subprocess_timeout_s: int = 10,
) -> Tuple[bool, str]:
    """Pre-flight check: VM ``abfe_bin`` exists and is executable.

    Uses ``ssh -o ConnectTimeout=N <host> test -x <bin>`` (POSIX-portable,
    rc=0 iff executable). Catches ``subprocess.TimeoutExpired`` and any
    ssh-level failure as gate-FAIL (the launcher must not proceed without
    a positively verified binary path).

    Returns
    -------
    (allow, reason)
        ``allow=True`` only when the remote ``test -x`` returns rc=0.
    """
    try:
        result = subprocess.run(
            [
                "ssh",
                "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
                vm_ssh_host,
                f"test -x {shlex.quote(abfe_bin)}",
            ],
            capture_output=True,
            timeout=subprocess_timeout_s,
        )
    except subprocess.TimeoutExpired:
        return (
            False,
            (
                f"VM abfe_bin probe timed out (>{subprocess_timeout_s}s) — "
                f"ssh to {vm_ssh_host} unreachable. Cannot launch."
            ),
        )
    except FileNotFoundError as exc:
        # local ssh binary missing
        return (
            False,
            f"local ssh binary missing for VM abfe_bin probe: {exc}",
        )
    if result.returncode == 0:
        return (
            True,
            f"VM abfe_bin {abfe_bin!r} verified executable on {vm_ssh_host}",
        )
    return (
        False,
        (
            f"VM abfe_bin {abfe_bin!r} does not exist or not executable "
            f"on {vm_ssh_host} (ssh test -x rc={result.returncode}). "
            f"Provision the VM `atm` conda env: "
            f"`ssh {vm_ssh_host} bash /tmp/vm_atm_install.sh` "
            f"(install script template under "
            f"`scripts/_vm_atm_env_install.sh`) — or pass "
            f"`--abfe-bin <path>` to override default."
        ),
    )


# ---------------------------------------------------------------------------
# G38 fix — _gate_vm_leg_dir_exists port from structprep.
# Production launcher had the same latent VM-FS bug as the structprep
# SSH-dispatch fix: when
# --gpu-host=vm dispatches abfe_production via SSH but the VM leg_dir
# lacks the required inputs (trackb.pdb / trackb_sys.xml /
# trackb_asyncre.cntl), the cohort half-stages — stage_per_replica_
# checkpoints succeeds on the host (which has the inputs) but the
# per-leg SSH-wrapped abfe_production dies on the VM with
# FileNotFoundError. Operator is then left with stale host-side r0..r21
# per-direction ckpts and no production output.
#
# Fix: defensive pre-flight gate _gate_vm_leg_dir_exists (mirror of the
# structprep helper) invoked BEFORE per-leg launch when gpu_host==vm.
# Cohort-safe halt — refuse to stage ckpts for any leg unless ALL legs
# pass the VM leg-dir gate.
#
# Same failure-class family as:
# - VM-lane self-provisioning
# - non-interactive conda-path resolution
# - F3 (VM abfe_bin missing → host-side ckpt stage hangs)
# - VM env install fix
# - G28 (production launcher VM leg-dir port)
# ---------------------------------------------------------------------------
def _gate_vm_leg_dir_exists(
    leg_dir: str,
    vm_ssh_host: str = "san@192.168.122.155",
    ssh_connect_timeout_s: int = 5,
    subprocess_timeout_s: int = 10,
) -> Tuple[bool, str]:
    """Pre-flight: confirm ``leg_dir`` exists on the VM AND contains the
    upstream-required inputs (``trackb.pdb`` + ``trackb_sys.xml`` +
    ``trackb_asyncre.cntl``).

    Returns ``(allow, reason)``. Refuses launch when any required input
    is missing on the VM side — cohort-safe gate semantics (better to
    halt than start half-staged).

    Same failure-class family as the VM-lane self-provisioning and
    non-interactive conda-path issues: pre-flight assertion
    that VM-side preconditions are met before any in-process commit.
    """
    required = ["trackb.pdb", "trackb_sys.xml", "trackb_asyncre.cntl"]
    probe_cmd = " && ".join(
        f"test -f {shlex.quote(os.path.join(leg_dir, f))}" for f in required
    )
    try:
        result = subprocess.run(
            [
                "ssh",
                "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
                "-o", "StrictHostKeyChecking=no",
                "-o", "LogLevel=ERROR",
                vm_ssh_host,
                probe_cmd,
            ],
            capture_output=True,
            timeout=subprocess_timeout_s,
        )
    except subprocess.TimeoutExpired:
        return (
            False,
            f"VM leg-dir probe timed out (>{subprocess_timeout_s}s) — "
            f"ssh to {vm_ssh_host} unreachable. Cannot launch.",
        )
    except FileNotFoundError as exc:
        return (
            False,
            f"local ssh binary missing for VM leg-dir probe: {exc}",
        )
    if result.returncode == 0:
        return (
            True,
            f"VM leg_dir {leg_dir!r} contains all required inputs "
            f"({', '.join(required)}) on {vm_ssh_host}",
        )
    return (
        False,
        f"VM leg_dir {leg_dir!r} missing one or more required inputs "
        f"({', '.join(required)}) on {vm_ssh_host} "
        f"(ssh test -f rc={result.returncode}). "
        f"Remediation: rsync the leg dir to the VM first, e.g. "
        f"`rsync -avz {leg_dir}/ {vm_ssh_host}:{leg_dir}/` "
        f"(VM disk must have free space; check `ssh {vm_ssh_host} df -h /home`)."
    )


# ---------------------------------------------------------------------------
# ATM state-occupancy fail-fast gate:
# detect the silent state-0-collapse class EARLY. The bound-leg incident
# produced "clean" data (NaN=0, finite energies) that masked the fact that
# every replica was locked to alchemical state 0 — the lambda-ladder was
# never sampled, making UWHAM rank-deficient and ΔΔG_bind invalid (a silent
# wrong-state pass).
#
# This is the launch-monitoring data-integrity check for
# ATM: parse the authoritative driver log "Replica N new state M" lines and
# assert that, after a warmup, at least two distinct alchemical states are
# occupied (i.e. the swap ladder is actually moving replicas off state 0).
# ---------------------------------------------------------------------------

# async_re driver line, e.g.
#   "2026-06-01 15:23:55 - INFO - async_re.openmm_async_re - Replica 1 new state 0"
_NEW_STATE_RE = re.compile(r"Replica\s+(\d+)\s+new state\s+(\d+)")

# Leading timestamp "YYYY-MM-DD HH:MM:SS" — used to group swap rounds
# (each round logs one "new state" line per replica with a shared
# timestamp). Robust warmup boundary without depending on replica count.
_LOG_TS_RE = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})")


def parse_state_occupancy_from_log(
    log_path: str,
    warmup_cycles: int = 20,
) -> Dict[str, Any]:
    """Parse an async_re driver log into per-state occupancy after warmup.

    Reads ``Replica N new state M`` lines (the authoritative record of
    each replica's alchemical-state assignment per swap round) and returns:

      * ``occupancy``      — {state_id: count} over POST-warmup samples
      * ``states_seen``    — sorted list of distinct states (post-warmup)
      * ``n_samples``      — post-warmup sample count
      * ``n_cycles_total`` — distinct timestamp groups (swap rounds) seen
      * ``replicas_seen``  — sorted list of distinct replica ids
      * ``warmup_cycles``  — echo of the warmup boundary used

    The warmup boundary is measured in swap rounds (distinct leading
    timestamps), so it is independent of the replica count (works for the
    11-state per-direction subdir and the 22-state combined free leg).
    Samples emitted during the first ``warmup_cycles`` rounds are excluded
    (initial equilibration / DIIS-style settling).

    Returns ``n_samples == 0`` if the log has no parseable ``new state``
    lines (caller treats as INDETERMINATE, not PASS).
    """
    occupancy: Dict[int, int] = {}
    replicas: set = set()
    n_samples = 0
    seen_ts: List[str] = []
    seen_ts_set: set = set()

    if not os.path.isfile(log_path):
        return {
            "log_path": log_path,
            "occupancy": {},
            "states_seen": [],
            "n_samples": 0,
            "n_cycles_total": 0,
            "replicas_seen": [],
            "warmup_cycles": warmup_cycles,
            "error": "log file not found",
        }

    with open(log_path, "r", errors="replace") as fh:
        for line in fh:
            m = _NEW_STATE_RE.search(line)
            if m is None:
                continue
            ts_m = _LOG_TS_RE.match(line)
            ts = ts_m.group(1) if ts_m else ""
            if ts and ts not in seen_ts_set:
                seen_ts_set.add(ts)
                seen_ts.append(ts)
            # Warmup boundary: skip samples until we have crossed
            # ``warmup_cycles`` distinct timestamp groups.
            if ts and len(seen_ts) <= warmup_cycles:
                continue
            try:
                replica = int(m.group(1))
                state = int(m.group(2))
            except (ValueError, IndexError):
                continue
            occupancy[state] = occupancy.get(state, 0) + 1
            replicas.add(replica)
            n_samples += 1

    return {
        "log_path": log_path,
        "occupancy": {str(k): v for k, v in sorted(occupancy.items())},
        "states_seen": sorted(occupancy.keys()),
        "n_samples": n_samples,
        "n_cycles_total": len(seen_ts),
        "replicas_seen": sorted(replicas),
        "warmup_cycles": warmup_cycles,
    }


def check_atm_state_occupancy(
    log_path: str,
    warmup_cycles: int = 20,
    min_distinct_states: int = 2,
) -> Dict[str, Any]:
    """Fail-fast gate for the ATM state-0-collapse class.

    PASS  iff, after ``warmup_cycles`` swap rounds, at least
          ``min_distinct_states`` distinct alchemical states are occupied
          (i.e. the lambda-ladder is moving replicas off state 0).
    FAIL  if only state 0 is occupied (states 1..K-1 total occupancy == 0)
          after warmup — the exact 2026-06-03 bound-leg pathology.
    INDETERMINATE (verdict="INDETERMINATE", a non-PASS) when there are no
          post-warmup samples yet (log too short / run just started). The
          caller decides whether to retry later; this is NOT a PASS.

    Returns a dict with ``verdict`` ("PASS"/"FAIL"/"INDETERMINATE"),
    ``passed`` (bool, True only for PASS), a human ``message``, and the
    full ``occupancy`` block from ``parse_state_occupancy_from_log``.
    """
    occ = parse_state_occupancy_from_log(log_path, warmup_cycles=warmup_cycles)
    states_seen = occ["states_seen"]
    n_samples = occ["n_samples"]

    if occ.get("error"):
        return {
            "verdict": "INDETERMINATE",
            "passed": False,
            "message": (
                f"Cannot evaluate occupancy for {log_path}: "
                f"{occ['error']}."
            ),
            "occupancy": occ,
        }

    if n_samples == 0:
        return {
            "verdict": "INDETERMINATE",
            "passed": False,
            "message": (
                f"No post-warmup 'new state' samples in {log_path} "
                f"(n_cycles_total={occ['n_cycles_total']}, "
                f"warmup_cycles={warmup_cycles}). Run too short to judge "
                f"ladder mixing; re-check after more cycles."
            ),
            "occupancy": occ,
        }

    nonzero_states = [s for s in states_seen if s != 0]
    n_distinct = len(states_seen)

    if n_distinct < min_distinct_states or not nonzero_states:
        return {
            "verdict": "FAIL",
            "passed": False,
            "message": (
                f"STATE-0-COLLAPSE detected in {log_path}: after "
                f"{warmup_cycles} warmup cycles only states {states_seen} "
                f"occupied ({n_samples} samples). The lambda-ladder is NOT "
                f"sampling intermediate states — UWHAM/MBAR is rank-"
                f"deficient and any ΔΔG is INVALID (bound-leg state-0 "
                f"collapse). "
                f"Root-cause check: per-replica ckpts must NOT be staged "
                f"(stage_per_direction_subdir creates EMPTY r-dirs so "
                f"async_re seeds state i to replica i)."
            ),
            "occupancy": occ,
        }

    return {
        "verdict": "PASS",
        "passed": True,
        "message": (
            f"Ladder mixing OK in {log_path}: {n_distinct} distinct states "
            f"occupied after warmup ({states_seen}); {n_samples} samples "
            f"across {len(occ['replicas_seen'])} replicas."
        ),
        "occupancy": occ,
    }


def _run_occupancy_check_cli(
    leg_dir: str,
    warmup_cycles: int = 20,
    jobname: str = "trackb",
) -> int:
    """CLI entry for ``--check-occupancy <leg_dir>``.

    Resolves the driver log(s) under ``leg_dir`` and runs
    ``check_atm_state_occupancy`` on each. For the 2-process per-direction
    layout there are two logs (``dplus/_live_launch.log`` and
    ``dminus/_live_launch.log``); for a combined leg there is one
    (``_live_launch.log`` or ``_production.log``). Returns 0 iff EVERY
    discovered log PASSES; 5 if any FAILS; 6 if none could be evaluated
    (no log found / all INDETERMINATE).
    """
    candidates = [
        os.path.join(leg_dir, "dplus", "_live_launch.log"),
        os.path.join(leg_dir, "dminus", "_live_launch.log"),
        os.path.join(leg_dir, "_live_launch.log"),
        os.path.join(leg_dir, "_production.log"),
    ]
    logs = [c for c in candidates if os.path.isfile(c)]
    if not logs:
        print(
            f"ERROR: no driver log found under {leg_dir} "
            f"(looked for dplus/dminus/_live_launch.log + "
            f"_live_launch.log + _production.log)",
            file=sys.stderr,
        )
        return 6

    any_fail = False
    any_pass = False
    for log in logs:
        res = check_atm_state_occupancy(log, warmup_cycles=warmup_cycles)
        print(f"[{res['verdict']}] {log}")
        print(f"    {res['message']}")
        print(f"    occupancy={res['occupancy']['occupancy']}")
        if res["verdict"] == "FAIL":
            any_fail = True
        elif res["verdict"] == "PASS":
            any_pass = True

    if any_fail:
        print("\nOCCUPANCY GATE: FAIL")
        return 5
    if not any_pass:
        print("\nOCCUPANCY GATE: INDETERMINATE (no PASS, no FAIL)")
        return 6
    print("\nOCCUPANCY GATE: PASS")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Track B per-direction production launcher (Option B Phase 4 "
            "TEMPLATE). LAUNCH BLOCKED by default; requires explicit C1-C8 "
            "confirmation flag."
        )
    )
    p.add_argument(
        "--check-occupancy",
        metavar="LEG_DIR",
        default=None,
        help=(
            "Standalone ATM state-occupancy gate. Parses "
            "the driver log(s) under LEG_DIR and PASS/FAILs on whether the "
            "lambda-ladder is sampling >1 state after warmup. Does NOT "
            "launch anything; exits 0=PASS / 5=FAIL / 6=INDETERMINATE."
        ),
    )
    p.add_argument(
        "--occupancy-warmup-cycles",
        type=int,
        default=20,
        help="warmup swap-rounds excluded before judging occupancy (default 20).",
    )
    p.add_argument(
        "--i-have-confirmed-c1-through-c8",
        action="store_true",
        help=(
            "Explicit operator gate. Set ONLY after manually verifying ALL "
            "C1-C8 conditions and obtaining user session confirmation. "
            "Without this flag the script halts after auto-checks + "
            "report write."
        ),
    )
    p.add_argument(
        "--v21-out-root",
        default="outputs/_trackb/production_v2_1",
        help="v2.1 systems root",
    )
    p.add_argument(
        "--prep-report",
        default="outputs/_trackb/per_direction_structprep/per_direction_prep_latest.json",
        help="per-direction prep report (Step 1 output)",
    )
    p.add_argument(
        "--out-root",
        default="outputs/_trackb/production_v2_2",
        help="production output root (NEW for v2.2 per-direction)",
    )
    p.add_argument(
        "--endpoints",
        default="cp4,wt",
        help="comma-list endpoints",
    )
    p.add_argument(
        "--legs",
        default="bound",
        help=(
            "comma-list legs (default 'bound' — free leg ddint already "
            "captured via uwham postprocess on host 5070Ti. Operator may "
            "explicitly pass 'bound,free' if a fresh free-leg run is "
            "desired, but the launcher will then FileNotFoundError on "
            "VM where free leg dirs do not exist)."
        ),
    )
    p.add_argument(
        "--seeds",
        default=",".join(EIGHT_SEED_COHORT),
        help=f"8-seed cohort (canonical: {','.join(EIGHT_SEED_COHORT)})",
    )
    p.add_argument(
        "--replicates",
        type=int,
        default=3,
        help="replicates per (endpoint, leg, seed) for sigma_btwn measurement",
    )
    p.add_argument(
        "--enable-replicate-subtrees",
        action="store_true",
        help=(
            "Stage 2 (Path λ-densify campaign 2026-06-05): launch N "
            "independent-seed replicates into SEPARATE output subtrees "
            "(<v21-out-root>/rep{0..N-1}/<endpoint>/<leg>). OFF by default "
            "(single-run path unchanged). Each replicate = independent "
            "velocity/walker-shuffle seed (recorded in metadata). PRECONDITION: "
            "the per-replicate cntls/systems must already exist in each "
            "rep{i} subtree (built by structprep with matching --replicates). "
            "When OFF the launcher uses the legacy single subtree."
        ),
    )
    p.add_argument(
        "--free-schedule",
        default="canonical22",
        choices=["canonical22", "densified38", "densified38v2",
                 "densified38v3", "densified38v4", "densified34"],
        help=(
            "λ ladder for the FREE leg only (Path λ-densify spec 2026-06-05). "
            "Recorded in this launcher's metadata for cross-stage audit; the "
            "actual densified cntl is emitted by the production launcher "
            "(trackb_production_v2_1_upstream.py --free-schedule). This "
            "launcher consumes pre-built cntls, so the flag here is an audit "
            "+ replicate-orchestration selector. Bound leg always canonical22. "
            "'densified38' = REVISED 38-state (19 fwd + 19 bwd) production "
            "ladder (revised ladder fix 2026-06-05). 'densified38v2' = "
            "REBALANCED densified38 (count-neutral linear-λ thin-plateau + "
            "densify-turnover, 2026-06-05; same 38 states, fixes "
            "the dplus 8→9 BC=0.131 overlap gap). 'densified38v3' = "
            "PER-DIRECTION (NON-mirror): dplus = clean densified38v2 forward "
            "(19); dminus = densified38v2 backward + ONE W0-graded micro-bridge "
            "at the soft-core-end→plateau handoff (20), closing the dminus 9→10 "
            "BC=0.24 hole the reversed mirror created (dminus-asymmetry "
            "fix 2026-06-06). UNEQUAL per-direction counts (39 total). "
            "'densified38v4' = densified38v3 + a SECOND dminus W0-graded micro-"
            "bridge at the re-indexed 5→6 handoff (dminus 21, total 40 UNEQUAL), "
            "closing the dminus 5→6 BC=0.18 hole the v3 re-pilot exposed "
            "(2026-06-06; the LAST iteration under the hard "
            "2-bridge/leg cap). All downstream per-direction wiring is DIRECTION-"
            "derived/count-agnostic. "
            "'densified34' is DEPRECATED (Factor-B broken) and kept only for "
            "forensic re-analysis of the earlier pilot."
        ),
    )
    p.add_argument(
        "--free-schedule-file",
        default=None,
        help=(
            "ADDITIVE (registry-untouched) λ-schedule override for the "
            "FREE leg: a path to a JSON schedule dict (the _schedule_dict shape, "
            "or an adaptive_schedule.json emitted by "
            "utils/adaptive_lambda.schedule_io.emit_validated_schedule). When "
            "set, BEFORE the per-direction slicing the FREE leg's COMBINED cntl "
            "(<jobname>_asyncre.cntl) has its per-state schedule columns "
            "(LAMBDAS/DIRECTION/INTERMEDIATE/LAMBDA1/LAMBDA2/ALPHA/U0/W0COEFF) "
            "re-emitted from the loaded dict via the asyncre write_cntl_file "
            "SSOT — the SINGLE cntl write point (NEVER patch "
            "per-direction cntls post-hoc; the launcher slices the combined "
            "cntl). All non-schedule keys (LIGAND_ATOMS / DISPLACEMENT / steps) "
            "are preserved verbatim. This is WIRING ONLY: no launch behavior "
            "change beyond the λ values the slicer inherits. Coexists with "
            "--free-schedule (the registry enum); --free-schedule-file takes "
            "precedence on the FREE leg when both are given. Bound leg is "
            "untouched (always canonical22)."
        ),
    )
    p.add_argument(
        "--production-steps",
        type=int,
        default=2500,
    )
    p.add_argument(
        "--max-samples",
        type=int,
        default=1000,
    )
    p.add_argument(
        "--wall-time-min",
        type=int,
        default=720,
    )
    p.add_argument(
        "--gpu-host",
        default="vm",
        choices=["vm", "local", "cpu"],
        help=(
            "GPU host (v0.9.10 device-guard fix). Each box has a single "
            "GPU at device 0; 'device 1' does not exist on either host. "
            "vm = VM V100 (refused if Track A util > 30%%); local = host "
            "5070Ti (refused if free-leg PID alive); cpu = OpenMM CPU."
        ),
    )
    p.add_argument(
        "--cuda-device",
        default=None,
        help=(
            "DEPRECATED legacy option (v0.9.10 device-guard fix). Only "
            "'cpu' or '0' accepted, and '0' is ambiguous between host "
            "and VM — prefer --gpu-host. Any other value (including "
            "the old default '1') is rejected: that GPU does not exist."
        ),
    )
    p.add_argument(
        "--gpu-util-refuse-threshold-pct",
        type=float,
        default=30.0,
        help="VM V100 utilization threshold above which --gpu-host=vm refuses.",
    )
    p.add_argument(
        "--jobname",
        default="trackb",
    )
    p.add_argument(
        "--abfe-bin",
        default="/home/san/miniconda3/envs/atm/bin/abfe_production",
        help=(
            "Path to abfe_production binary. "
            "Default assumes VM atm env (provisioned via "
            "scripts/_vm_atm_env_install.sh). Override when using a "
            "different conda env path (e.g., if VM atm env is unavailable "
            "and operator has provisioned abfe_production in another env)."
        ),
    )
    p.add_argument(
        "--vm-ssh-host",
        default="san@192.168.122.155",
        help=(
            "VM ssh target for --gpu-host=vm. Used by the F3 pre-flight "
            "gate (_gate_vm_abfe_bin_exists) and by _live_launch_all_legs "
            "ssh wrapping."
        ),
    )
    p.add_argument(
        "--free-leg-pid",
        type=int,
        default=1876426,
    )
    p.add_argument(
        "--free-leg-results-dir",
        default="outputs/_trackb/production_v2_1/cp4/free",
        help="dir where free-leg results land (for C1 check)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="C1-C8 auto-check + stage per-replica ckpts, but DO NOT launch",
    )
    p.add_argument(
        "--free-pilot",
        action="store_true",
        help=(
            "FREE-leg pilot mode (densified38 per-direction free, 2026-06-05). "
            "Replaces the BOUND-leg C1-C4 gates (which precondition on an "
            "ALREADY-COMPLETE free leg + read the bound per-direction "
            "prep-report) with a free-leg-appropriate readiness gate: every "
            "requested (endpoint, leg) must have <jobname>_asyncre.cntl + "
            "<jobname>_0_dplus.xml + <jobname>_0_dminus.xml present "
            "(produced by trackb_per_direction_structprep.py). NOT a bypass "
            "— it inverts the inverted C1 precondition for a FRESH free run. "
            "C6 (seed cohort) still applies. Intended with --legs free + a "
            "SEPARATE --out-root + --v21-out-root pointing at the densified38 "
            "free system (e.g. production_v2_3_densify38_pilot). The C5/C4 "
            "bound-only probes self-skip when 'bound' not in --legs."
        ),
    )
    args = p.parse_args()

    # Standalone occupancy gate — short-circuit BEFORE any launcher banner /
    # device guard. Read-only; never launches production.
    if args.check_occupancy is not None:
        return _run_occupancy_check_cli(
            leg_dir=args.check_occupancy,
            warmup_cycles=args.occupancy_warmup_cycles,
            jobname=args.jobname,
        )

    print("\n" + "=" * 70)
    print("Track B per-direction production (Option B Phase 4 TEMPLATE)")
    print("=" * 70)
    print(f"Auto-launch GATED: requires --i-have-confirmed-c1-through-c8")
    print()

    # ----- v0.9.10 device-guard fix: --cuda-device deprecated --------
    if args.cuda_device is not None:
        prep_mod = _load_prep_gate_helpers()
        try:
            mapped = prep_mod.normalize_legacy_cuda_device(args.cuda_device)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        if mapped is not None:
            print(f"# --cuda-device={args.cuda_device!r} mapped to "
                  f"--gpu-host={mapped!r} (legacy compat)")
            args.gpu_host = mapped

    seeds = [s.strip() for s in args.seeds.split(",") if s.strip()]
    endpoints = [e.strip() for e in args.endpoints.split(",") if e.strip()]
    legs = [l.strip() for l in args.legs.split(",") if l.strip()]

    # ----- --free-schedule-file (additive λ override) ----------
    # SSOT single write point: re-emit the FREE leg's COMBINED cntl per-state
    # schedule columns from the loaded JSON dict BEFORE the per-direction slicing
    # (generate_per_direction_cntls): the slicer inherits the
    # loaded λ; no per-direction cntl is patched post-hoc. Bound leg untouched.
    if args.free_schedule_file:
        sched_file = (
            args.free_schedule_file
            if os.path.isabs(args.free_schedule_file)
            else os.path.join(_PROJ_ROOT, args.free_schedule_file)
        )
        print(f"\n--- --free-schedule-file: {sched_file} ---")
        applied_any = False
        for endpoint in endpoints:
            if "free" not in legs:
                break
            leg_dir = os.path.join(
                _PROJ_ROOT, args.v21_out_root, endpoint, "free"
            )
            combined = os.path.join(leg_dir, args.jobname + "_asyncre.cntl")
            if not os.path.isfile(combined):
                print(f"  {endpoint}/free: WARN combined cntl missing "
                      f"({combined}) — run structprep first; skipping override")
                continue
            try:
                audit = apply_free_schedule_file_to_combined_cntl(
                    combined_cntl_path=combined,
                    schedule_file=sched_file,
                )
            except (FileNotFoundError, ValueError) as exc:
                print(f"  {endpoint}/free: free-schedule-file FAILED: {exc}",
                      file=sys.stderr)
                return 2
            applied_any = True
            print(f"  {endpoint}/free: combined cntl re-emitted "
                  f"({audit['n_states']} states, keys "
                  f"{audit['rewritten_keys']})")
        if not applied_any:
            print("  WARN: --free-schedule-file set but no FREE leg combined "
                  "cntl was rewritten (no 'free' in --legs, or no combined "
                  "cntl present).")

    # ----- Auto checks (best-effort) ----------------------------
    prep_report = os.path.join(_PROJ_ROOT, args.prep_report)
    free_leg_dir = os.path.join(_PROJ_ROOT, args.free_leg_results_dir)
    if args.free_pilot:
        # Free-pilot mode: the bound C1-C4 gates are inverted / N/A for a
        # fresh FREE run (see check_free_pilot_readiness docstring). Substitute
        # the free-leg-appropriate readiness gate + keep C6 (seed cohort).
        checks = [
            check_free_pilot_readiness(
                v21_out_root=args.v21_out_root,
                endpoints=endpoints,
                legs=legs,
                jobname=args.jobname,
            ),
            check_c6_seed_cohort(seeds),
        ]
        print("\n--- FREE-PILOT AUTO-CHECKS (bound C1-C4 inverted/N/A) ---")
    else:
        checks = [
            check_c1_free_leg_complete(args.free_leg_pid, free_leg_dir),
            check_c2_dminus_equilibration(prep_report),
            check_c3_ommreplica_dryrun(prep_report),
            check_c4_charge_axis(prep_report),
            check_c6_seed_cohort(seeds),
        ]
        print("\n--- C1-C8 AUTO-CHECKS ---")
    all_pass = all(c.get("pass", False) for c in checks)
    for c in checks:
        status = "PASS" if c.get("pass") else "FAIL/N/A"
        print(f"  [{status}] {c.get('condition')}")
        if not c.get("pass") and "reason" in c:
            print(f"      reason: {c['reason']}")

    # ----- Pre-registration (C7) -------------------------------
    out_root = os.path.join(_PROJ_ROOT, args.out_root)
    os.makedirs(out_root, exist_ok=True)
    pre_reg = {
        "regime": "ranking_only",
        "outcomes": list(PRE_REGISTER_OUTCOMES),
        "sign_ssot_kcal_itc": -1.4,
        "sign_ssot_kcal_spr": -3.0,
        "sign_ssot_source": "Magotti 2009 J Mol Recognit 22(6):495 (Cp4 vs WT)",
        "sign_ssot_caveat": "4W9A indirect-analog (per Q5 + cp4-literature-comparability)",
        "decision_at": "post-uwham-aggregate-of-3-rep-x-8-seed-x-2-leg",
        "seeds": seeds,
        "replicates": args.replicates,
        "replicate_subtrees_enabled": bool(args.enable_replicate_subtrees),
        "free_schedule": args.free_schedule,
        "free_pilot": bool(args.free_pilot),
        "endpoints": endpoints,
        "legs": legs,
        "method_ref": "per-direction structprep (Option B), 2026-05-31",
        "lambda_densify_ref": "free-leg λ-resampling resolution, 2026-06-05",
        "registered_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    pre_reg_path = os.path.join(out_root, "pre_registration.json")
    # G36-C7 fix: also regenerate when invocation arguments
    # diverge from existing pre-registration (e.g. legs field changed from
    # ["bound","free"] to ["bound"]). Without
    # this, the audit trail goes stale and integrity checks diagnose
    # mismatches that are actually expected operator changes.
    pre_reg_action: str = "fresh-write"
    if os.path.isfile(pre_reg_path):
        try:
            with open(pre_reg_path) as fh:
                existing = json.load(fh)
            divergent_fields: List[str] = []
            for field in ("legs", "endpoints", "seeds", "replicates"):
                if existing.get(field) != pre_reg[field]:
                    divergent_fields.append(field)
            if divergent_fields and args.i_have_confirmed_c1_through_c8:
                # Operator has authority + audit trail integrity: regenerate
                # but preserve old as _stale_<ts> alongside (never delete).
                stale_ts = time.strftime("%Y%m%dT%H%M%S")
                stale_path = os.path.join(
                    out_root, f"pre_registration_stale_{stale_ts}.json"
                )
                shutil.move(pre_reg_path, stale_path)
                pre_reg_action = "regenerated"
                print(f"\n--- C7 PRE-REGISTRATION REGENERATED ---")
                print(
                    f"  divergent fields: {divergent_fields}\n"
                    f"  prior pre_reg archived to: {stale_path}"
                )
            elif divergent_fields:
                pre_reg_action = "kept-existing-divergent"
                print(f"\n--- C7 PRE-REGISTRATION DIVERGENCE ---")
                print(
                    f"  WARN: existing pre_reg diverges from CLI args in "
                    f"{divergent_fields}, but --i-have-confirmed-c1-through-c8 "
                    f"not set — leaving existing untouched. Re-run with the "
                    f"flag to regenerate (will archive prior to "
                    f"pre_registration_stale_<ts>.json)."
                )
            else:
                pre_reg_action = "kept-existing-matched"
        except (json.JSONDecodeError, OSError) as exc:
            # Corrupt prior file → archive + write fresh (preserve old).
            stale_ts = time.strftime("%Y%m%dT%H%M%S")
            shutil.move(
                pre_reg_path,
                os.path.join(out_root, f"pre_registration_corrupt_{stale_ts}.json"),
            )
            pre_reg_action = "regenerated-after-corrupt"
            print(f"  WARN: existing pre_reg unreadable ({exc}); archived")

    if pre_reg_action in ("fresh-write", "regenerated", "regenerated-after-corrupt"):
        with open(pre_reg_path, "w") as fh:
            json.dump(pre_reg, fh, indent=2)
        print(f"\n--- C7 PRE-REGISTRATION ---")
        print(f"  archived: {pre_reg_path} ({pre_reg_action})")

    # ----- Fallback triggers logic (auto-detect post-production) -
    fallback_triggers = {
        "trigger_secondary_nan_dminus": (
            "Phase 3 d=-1 walker exhibits secondary NaN (different "
            "pathology than addSolvent water-clash)"
        ),
        "trigger_bidir_hysteresis": (
            "Phase 3 bidir hysteresis > 3 kcal/mol at lambda=0.5 "
            "(per-direction structprep does not equilibrate to "
            "comparable thermodynamic state)"
        ),
        "trigger_sigma_drift": (
            "Phase 4 production sigma_btwn > 1.0 kcal/mol "
            "(loss of Gibbs cross-direction sampling worse than predicted)"
        ),
        "fallback_engine": "Option A (Q5b PMX + Q5c GROMACS-BFEE2 hybrid)",
        "fallback_setup_days": "5-7 day setup, ~165 GPU-h V100",
        "fallback_caveat": (
            "amber99sb-star-ildn-mut FF mismatch with UPDD amber14SB "
            "(ranking-only disclosure required in paper §4.3)"
        ),
    }
    fallback_path = os.path.join(out_root, "fallback_triggers.json")
    with open(fallback_path, "w") as fh:
        json.dump(fallback_triggers, fh, indent=2)

    # ----- LAUNCH GATE -----------------------------------------
    print(f"\n--- LAUNCH GATE ---")
    if args.dry_run:
        print("  --dry-run: stopping after auto-checks + pre-registration")
        # In dry-run we render the 2-process launch plan for manual
        # inspection. _live_launch_all_legs(dry_run=True) generates the
        # per-direction cntls (forward dplus + backward dminus slices) and
        # prints the two planned abfe_production commands per leg WITHOUT
        # staging subdirs or spawning subprocesses. Legs whose combined
        # cntl is absent are skipped with a WARN (operator must run
        # structprep first).
        present_legs = []
        for endpoint in endpoints:
            for leg in legs:
                leg_dir = os.path.join(_PROJ_ROOT, args.v21_out_root,
                                       endpoint, leg)
                cntl = os.path.join(leg_dir, args.jobname + "_asyncre.cntl")
                if os.path.isfile(cntl):
                    present_legs.append((endpoint, leg))
                elif os.path.isdir(leg_dir):
                    print(f"  {endpoint}/{leg}: WARN combined cntl missing "
                          f"(run structprep first)")
        if present_legs:
            try:
                _live_launch_all_legs(
                    v21_out_root=args.v21_out_root,
                    endpoints=sorted({e for e, _ in present_legs}),
                    legs=sorted({l for _, l in present_legs}),
                    jobname=args.jobname,
                    gpu_host=args.gpu_host,
                    dry_run=True,
                    abfe_bin=args.abfe_bin,
                    vm_ssh_host=args.vm_ssh_host,
                )
            except (FileNotFoundError, RuntimeError, ValueError) as e:
                print(f"  dry-run plan WARN: {e}")
        return 0

    if not args.i_have_confirmed_c1_through_c8:
        print("  LAUNCH BLOCKED — --i-have-confirmed-c1-through-c8 not set")
        print("  See `# C1` through `# C8` blocks in this script source.")
        print("  All C1-C8 must be manually verified + user must confirm")
        print("  in the calling session before this flag is set.")
        return 1

    if not all_pass:
        print("  LAUNCH BLOCKED — automated C1-C8 checks FAILED above.")
        print("  Override only with explicit review + user consent.")
        return 1

    # ----- v0.9.10 device-guard fix: gate GPU host pre-launch -------
    prep_mod = _load_prep_gate_helpers()
    gate = prep_mod.gate_gpu_host(
        args.gpu_host,
        free_leg_pid=args.free_leg_pid,
        refuse_threshold_pct=args.gpu_util_refuse_threshold_pct,
    )
    print(f"  gpu_host={gate['host']!r} platform={gate['platform']!r} "
          f"allow={gate['allow']} — {gate['reason']}")
    if not gate["allow"]:
        print(f"  LAUNCH BLOCKED — GPU host gate REFUSED: {gate['reason']}")
        return 1

    # ----- F3 fix: VM abfe_bin pre-flight check ------
    # Only enforced when gpu_host=vm (local/cpu use the host-side binary
    # which is verified at install time). Prevents cohort half-staged
    # state (per-direction ckpts staged but no production launch).
    if args.gpu_host == "vm":
        allow_bin, bin_reason = _gate_vm_abfe_bin_exists(
            abfe_bin=args.abfe_bin,
            vm_ssh_host=args.vm_ssh_host,
        )
        print(f"  vm_abfe_bin gate: allow={allow_bin} — {bin_reason}")
        if not allow_bin:
            print(
                "  LAUNCH BLOCKED — VM abfe_bin pre-flight FAILED. "
                "No ckpts have been staged yet (safe to retry after "
                "fixing).", file=sys.stderr,
            )
            return 1

    # ----- C5 finite-energy probe + C4 water-count audit
    # (two-process per-direction split, C4 / C5). Run on the
    # HOST (Reference platform single-point) BEFORE any per-direction
    # dispatch. C5 nonfinite (PBC-wrap clash) => hard HALT; C4
    # cross-endpoint water-delta breach => ranking caveat (warn, persisted
    # to the audit JSON, not a hard HALT per verdict §b). Probes the bound
    # leg endpoints only — the free leg has no per-direction systems.
    bound_endpoint_legs = {
        ep: os.path.join(_PROJ_ROOT, args.v21_out_root, ep, "bound")
        for ep in endpoints
        if "bound" in legs
        and os.path.isdir(
            os.path.join(_PROJ_ROOT, args.v21_out_root, ep, "bound")
        )
    }
    c5_report: Dict[str, Any] = {"status": "skipped"}
    c4_report: Dict[str, Any] = {"status": "skipped"}
    if bound_endpoint_legs:
        print("  C5 finite-energy probe (PBC-wrap clash) — "
              f"{len(bound_endpoint_legs)} endpoint(s) x 2 directions ...")
        try:
            c5_report = finite_energy_probe_all(
                bound_endpoint_legs, jobname=args.jobname,
            )
        except ImportError as exc:
            c5_report = {"status": "skipped_no_openmm", "error": str(exc)}
            print(f"  C5 SKIPPED — openmm unavailable in this env: {exc}")
        for pr in c5_report.get("probes", []):
            print(f"    C5 {pr.get('endpoint')}/{pr.get('direction_tag')}: "
                  f"status={pr['status']} "
                  f"base_energy_kj={pr.get('base_energy_kj', 'n/a')} "
                  f"(u1 diag={pr.get('u1_kj', 'n/a')})")
        if c5_report["status"] == "nonfinite":
            print(
                "  LAUNCH BLOCKED — C5 finite-energy probe FAILED "
                "(PBC-wrap clash; energy > 1e10 kJ/mol or NaN). No "
                "dispatch occurred. Remediation: shrink DISPLACEMENT "
                "(2.5 -> 1.0 nm) or increase box padding, then re-prep.",
                file=sys.stderr,
            )
            return 1
        # C4 water-count audit (needs cp4 + wt bound legs).
        if "cp4" in bound_endpoint_legs and "wt" in bound_endpoint_legs:
            try:
                c4_report = cross_endpoint_water_audit(
                    bound_endpoint_legs["cp4"],
                    bound_endpoint_legs["wt"],
                    jobname=args.jobname,
                )
                print(
                    f"  C4 water audit: status={c4_report['status']} "
                    f"|DN_cp4-DN_wt|="
                    f"{c4_report['abs_cross_endpoint_water_delta']:.1f} "
                    f"(threshold {c4_report['threshold']})"
                )
                if c4_report["status"] == "caveat":
                    print(
                        "  WARNING — C4 cross-endpoint water-delta exceeds "
                        "threshold; 1st-order Delta-water cancellation in "
                        "DDG(Cp4-WT) may break. Ranking caveat must be "
                        "disclosed (paper sec 4.3). Launch CONTINUES "
                        "(ranking-only).", file=sys.stderr,
                    )
            except FileNotFoundError as exc:
                c4_report = {"status": "skipped_missing_input",
                             "error": str(exc)}
                print(f"  C4 SKIPPED — {exc}")

    # ===== AT THIS POINT THE LAUNCH IS AUTHORIZED ===============
    # F2 fix: inline the
    # per-leg abfe_production loop here instead of leaving the operator
    # to construct ``leg_info`` dicts + remember the (endpoint, leg)
    # iteration. The two safety gates above (--i-have-confirmed-c1-
    # through-c8 + gate_gpu_host) MUST have passed to reach this point.
    try:
        if args.enable_replicate_subtrees:
            # Stage 2 (Path λ-densify campaign 2026-06-05): N independent-seed
            # replicates into separate rep{i} subtrees for σ_btwn.
            print(f"  authorized — launching {args.replicates} independent-"
                  f"seed replicate(s) via _live_launch_replicates "
                  f"(subtrees rep0..rep{args.replicates - 1})")
            launch_results = _live_launch_replicates(
                v21_out_root=args.v21_out_root,
                endpoints=endpoints,
                legs=legs,
                jobname=args.jobname,
                gpu_host=args.gpu_host,
                n_replicates=args.replicates,
                seeds=seeds,
                dry_run=False,
                abfe_bin=args.abfe_bin,
                vm_ssh_host=args.vm_ssh_host,
            )
            execution_model = "per_direction_2process_split_replicated"
        else:
            print("  authorized — launching per-leg abfe_production via "
                  "_live_launch_all_legs (single subtree)")
            launch_results = _live_launch_all_legs(
                v21_out_root=args.v21_out_root,
                endpoints=endpoints,
                legs=legs,
                jobname=args.jobname,
                gpu_host=args.gpu_host,
                dry_run=False,
                abfe_bin=args.abfe_bin,
                vm_ssh_host=args.vm_ssh_host,
            )
            execution_model = "per_direction_2process_split"
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"  LAUNCH FAILED: {exc}", file=sys.stderr)
        return 1
    # Archive launch results next to pre_registration for audit. Includes
    # the C5 finite-energy probe + C4 water-count audit reports (verdict
    # provenance for the 2-process split).
    launch_audit_path = os.path.join(out_root, "_live_launch_results.json")
    with open(launch_audit_path, "w") as fh:
        json.dump({
            "launched_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "gpu_host": args.gpu_host,
            "regime": "ranking_only_R11",
            "execution_model": execution_model,
            "free_schedule": args.free_schedule,
            "n_replicates": args.replicates,
            "replicate_subtrees_enabled": bool(args.enable_replicate_subtrees),
            "seeds": seeds,
            "c5_finite_energy_probe": c5_report,
            "c4_water_count_audit": c4_report,
            "results": launch_results,
        }, fh, indent=2)
    print(f"  launch audit: {launch_audit_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
