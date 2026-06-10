#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B production v2 — upstream ``atom_openmm.async_re`` per-endpoint ABFE.

Option 2 (PRIMARY). The v1 launcher (``scripts/trackb_production.py``) was
architecturally invalid because it ran a single OpenMM Context per replica
and flipped the ATMForce ``Direction`` global parameter mid-trajectory at
λ=0.5 (states 10→11). The ATMForce reference potential expression
``select(step(Direction), u0, u1)`` is a step function in the BASE
potential — flipping Direction in a live Context without re-initialization
discontinuously switches u0↔u1 → NaN at the first integrator step.
Engineering Phase 1 (1 fs warmup, 5000-iter re-minimization) FAILED
because the underlying defect is in the launcher topology, not the timestep.

This v2 launcher uses Gallicchio 2021 JCIM 61:5424 standard ``async_re``
(implemented in ``atom_openmm`` v8.4.0 ``abfe_production`` /
``abfe_structprep``): each direction is a SEPARATE replica walker with
its OWN OpenMM Context. Direction is set at Context construction and
never flipped. Bidirectional UWHAM analysis (Tan 2012 JCP 136:144102)
operates on per-replica reciprocal work distributions.

# Thermodynamic cycle (unchanged from v1)

    ΔG_bind(endpoint) = ΔG_alch(bound, λ=0→1) − ΔG_alch(free, λ=0→1)
    ΔΔG_bind_MTR↔Trp  = ΔG_bind(Cp4) − ΔG_bind(WT)

4 ABFE legs per production: cp4/bound, cp4/free, wt/bound, wt/free.

# Charge axis (unchanged)

Hybrid MTR XML (``params/MTR_gaff2_hybrid.xml``, Option-β/regime-2,
per-residue |Σq|≤5e-4 e, NE1 frozen at -0.3418, Khoury 2014 OMW
sidechain). Cyclic SS disulfide (CYS2-CYS12) preserved on both legs.

# Schedule (unchanged 22-state symmetric)

Canonical AToM-OpenMM 22-state schedule: 11 forward (direction=+1,
λ=0→0.5) + 11 backward (direction=-1, λ=0.5→0), intermediate=1 at the
two λ=0.5 midpoints. Identical to upstream
``examples/ABFE/temoa-g1/temoa-g1_asyncre.cntl``.

# Output schema (ranking-only)

Per AToM-OpenMM v8.4.0: ``r0/{jobname}.out`` ... ``r21/{jobname}.out``
where r0..r21 = 22 separate replica walkers (one per (λ, direction)
state). Per-leg directory tree:

    outputs/_trackb/production_v2/<endpoint>/<leg>/
        ├── trackb.cntl              # async_re config
        ├── trackb.pdb               # topology
        ├── trackb_sys.xml           # system
        ├── trackb_min.xml           # post-minimization
        ├── trackb_therm.xml         # post-thermalization
        ├── trackb_npt.xml           # post-NPT
        ├── trackb_equil.xml         # post-equilibration (NVT)
        ├── trackb_mdlambda.xml      # post-λ-annealing (mid-λ pivot)
        ├── trackb_0.xml             # production starting point
        ├── trackb_0.pdb
        ├── r{0..21}/                # per-replica walker output
        │     ├── trackb.out
        │     ├── trackb.xtc
        │     └── trackb_ckpt.xml
        ├── trackb_stat.txt          # live replica status
        └── _logs/                   # async_re log + structprep log

# Hardware contract

* Smoke + initial production: host 5070Ti (single GPU, CUDA device 0).
* Track A V100 stays untouched until Track A completion (~6/1 10:30 KST);
  thereafter, optionally migrate remaining legs to V100 by editing the
  nodefile + relaunching ``abfe_production``.
* 22 replica walkers on 1 GPU → throughput-limited by GPU compute, async
  RE scheduler queues replicas in fast-execution buffer.

# Cross-references

* v1 archive:    ``outputs/_trackb/_archive/v1_architectural_failed_20260531/``
* v1 launcher:   ``scripts/trackb_production.py`` (deprecated, retained
                  for historical post-mortem reference)
* Setup module:  ``utils/atm_trackB_setup.py`` (system XML build path
                  reused unchanged)
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import time
from typing import Optional, List, Dict, Tuple, Any

_PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "utils"))

# ---------------------------------------------------------------------------
# Canonical AToM 22-state schedule (unchanged from v1; required for
# test_atm_trackb_production.py regression coverage).
# ---------------------------------------------------------------------------
LAMBDA_FWD = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
LAMBDA_BWD = [0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
DIRECTIONS = [1] * 11 + [-1] * 11
LAMBDAS_1 = LAMBDA_FWD + LAMBDA_BWD
LAMBDAS_2 = LAMBDA_FWD + LAMBDA_BWD
INTERMD = [0] * 10 + [1, 1] + [0] * 10
W0 = [0.0] * 10 + [1.0, 1.0] + [0.0] * 10
ALPHA = [0.10] * 22
U0 = [110.0] * 22
TEMP_K = 300.0
N_STATES = len(DIRECTIONS)
assert N_STATES == 22, "Schedule must be 22 states (11+11 symmetric)"


# ---------------------------------------------------------------------------
# Densified 34-state free-leg schedule (λ-densify spec 2026-06-05).
#
# WHY this exists (do NOT replace the 22-state default with it): the free-leg
# ATM decoupling crossover at λ 0.45→0.50 collapses ~191 kcal/mol of binding
# ΔU within ONE λ step → adjacent states 9↔10 / 11↔12 have ZERO phase-space
# overlap → UWHAM cannot bridge them (free-leg ΔΔG_int error structurally
# floored ~±0.86 honest; analytic ±0.377 underestimate). The decoupling is
# driven by the W0 / soft-core anneal at the INTERMEDIATE midpoint, NOT by the
# λ value — so plain λ insertion (λ1=λ2, W0=0) lands in the SAME coupled basin
# and bridges nothing. The fix is a NON-LINEAR ilogistic anneal ladder:
# λ1<λ2 (softplus active) + a W0 ramp 0→1 across 7 windows.
#
# Forward half (17 states) = 10 linear states (λ 0.00..0.45, λ1=λ2, W0=0,
# INTERMEDIATE=0; identical to the 22-state forward leg) + 7-window anneal
# ladder (INTERMEDIATE=1, ALPHA=0.10, U0=110.0). Backward half (17 states) is
# the WHOLE-TUPLE reverse of the forward half (per-state arrays reversed; λ1
# stays λ1, λ2 stays λ2 — NOT a λ1↔λ2 swap) with DIRECTION=-1. Total = 34.
#
# This is a physically-grounded STARTING POINT, NOT a guaranteed solution: a
# short pilot MUST re-measure adjacent OVL (target >0.3 every step) before any
# full campaign. If any ladder step still <0.3, add windows there (escalation:
# u0 ramp 200→110 across the ladder). Bound leg is UNCHANGED (22-state) —
# the receptor holds the ligand in the 192-194 plateau across all λ (no
# crossover, no gap); its limiter is forward/backward hysteresis, fixed by
# more sampling/replicates, not λ-densification.
#
# Free leg ONLY. Selected via ``--free-schedule densified34`` (default
# ``canonical22``).
# ---------------------------------------------------------------------------
# Linear forward states (λ1=λ2=λ, W0=0, INTERMEDIATE=0).
DENSE34_LAMBDA_FWD_LINEAR = [
    0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
]
# Anneal ladder windows as (λ1, λ2, W0coeff) — INTERMEDIATE=1, softplus active.
# Exact reproduction of Path's LADDER (7 windows).
DENSE34_LADDER = [
    (0.45, 0.46, 0.20),
    (0.45, 0.47, 0.40),
    (0.46, 0.48, 0.55),
    (0.47, 0.49, 0.70),
    (0.48, 0.495, 0.85),
    (0.49, 0.50, 0.95),
    (0.50, 0.50, 1.00),
]


def _build_densified34_schedule() -> Dict[str, List[Any]]:
    """Construct the 34-state densified free-leg schedule arrays.

    Forward half (17) = 10 linear + 7-window anneal ladder; backward half (17)
    = whole-tuple reverse of forward half with DIRECTION=-1. ALPHA=0.10 and
    U0=110.0 across all 34 states.

    Returns a dict with keys ``lambdas_1`` / ``lambdas_2`` / ``lambdas`` /
    ``directions`` / ``intermd`` / ``w0`` / ``alpha`` / ``u0`` /
    ``n_states``. The ``lambdas`` column tracks λ2 in the ladder (matching
    Path's LAMBDAS array, which equals λ2 in the soft-core windows).
    """
    fwd_l1 = list(DENSE34_LAMBDA_FWD_LINEAR) + [t[0] for t in DENSE34_LADDER]
    fwd_l2 = list(DENSE34_LAMBDA_FWD_LINEAR) + [t[1] for t in DENSE34_LADDER]
    fwd_w0 = [0.0] * 10 + [t[2] for t in DENSE34_LADDER]
    fwd_int = [0] * 10 + [1] * 7
    fwd_lam = list(DENSE34_LAMBDA_FWD_LINEAR) + [t[1] for t in DENSE34_LADDER]

    # Backward half = reverse the forward half as whole per-state tuples (λ1
    # stays λ1, λ2 stays λ2; NOT a λ1↔λ2 swap). DIRECTION = -1.
    bwd_l1 = list(reversed(fwd_l1))
    bwd_l2 = list(reversed(fwd_l2))
    bwd_w0 = list(reversed(fwd_w0))
    bwd_int = list(reversed(fwd_int))
    bwd_lam = list(reversed(fwd_lam))

    lambdas_1 = fwd_l1 + bwd_l1
    lambdas_2 = fwd_l2 + bwd_l2
    lambdas = fwd_lam + bwd_lam
    intermd = fwd_int + bwd_int
    w0 = fwd_w0 + bwd_w0
    directions = [1] * 17 + [-1] * 17
    n_states = len(directions)
    alpha = [0.10] * n_states
    u0 = [110.0] * n_states
    return {
        "lambdas_1": lambdas_1,
        "lambdas_2": lambdas_2,
        "lambdas": lambdas,
        "directions": directions,
        "intermd": intermd,
        "w0": w0,
        "alpha": alpha,
        "u0": u0,
        "n_states": n_states,
    }


_DENSE34 = _build_densified34_schedule()
DENSE34_LAMBDAS_1 = _DENSE34["lambdas_1"]
DENSE34_LAMBDAS_2 = _DENSE34["lambdas_2"]
DENSE34_LAMBDAS = _DENSE34["lambdas"]
DENSE34_DIRECTIONS = _DENSE34["directions"]
DENSE34_INTERMD = _DENSE34["intermd"]
DENSE34_W0 = _DENSE34["w0"]
DENSE34_ALPHA = _DENSE34["alpha"]
DENSE34_U0 = _DENSE34["u0"]
DENSE34_N_STATES = _DENSE34["n_states"]
assert DENSE34_N_STATES == 34, "Densified schedule must be 34 states (17+17)"
assert sum(DENSE34_INTERMD) == 14, "Densified schedule must have 14 INTERMEDIATE windows"


# ---------------------------------------------------------------------------
# REVISED densified 38-state free-leg schedule (revised ladder fix
# 2026-06-05).
#
# WHY densified34 is DEPRECATED (do NOT use for production): the densified34
# smoke built cleanly but abfe_production NaN'd 6× at the BACKWARD W0-peak
# intermediates (replica 17 ×4 + replica 18 ×2; states 17/18 = λ=0.5, W0≈1.0,
# DIR=−1). The discriminator (full 850k-step / 303 K structprep) STILL NaN'd
# the backward intermediates → Factor B (backward-intermediate base
# mis-equilibration) is genuine and densification-amplified, NOT just the cold
# 265 K smoke structprep (Factor A). The 7-window flat-α/U0 ladder was too
# stiff at the backward peak.
#
# The REVISED 38-state ladder softens the backward peak in three ways
# (endpoint-INVARIANT per ATMForce source: usc is independent of
# α/Uh/W0, so ramping them reshapes only the integration PATH / overlap, not
# the endpoint ΔG; Kirkwood 1935 / Zwanzig 1954 path-independence):
#   (a) finer W0 near saturation (ΔW0 step 0.15 far from peak → 0.04 at peak),
#   (b) per-state ALPHA ramp 0.10→0.25 (gentler ilogistic knee → wider pertE
#       distribution → MORE overlap),
#   (c) per-state U0/Uh ramp 110→82 (knee below the capped ~191 plateau).
# UMAX/ACORE/UBCORE stay GLOBAL-fixed (ommsystem.py:445 — cannot ramp per-state
# without rewriting the ATMForce expression; out of scope).
#
# Forward half (19 states) = 10 linear states (λ 0.00..0.45, λ1=λ2, W0=0,
# ALPHA=0.10, U0=110, INTERMEDIATE=0; identical to the 22-state forward leg) +
# 9-window anneal ladder (INTERMEDIATE=1) with per-state (λ1, λ2, W0, ALPHA,
# U0). Backward half (19 states) is the WHOLE-TUPLE reverse of the forward half
# (per-state arrays reversed; λ1 stays λ1, λ2 stays λ2 — NOT a λ1↔λ2 swap) with
# DIRECTION=-1. Total = 38.
#
# This is still a physically-grounded STARTING POINT, NOT a guaranteed
# solution. The re-smoke pilot MUST verify, before any full campaign:
#   (1) NaN-free ≥20 production cycles, (2) adjacent OVL>0.3 at EVERY pair
#   (esp. backward W0-peak 17..23 + 9→10/11→12 crossover), (3) 38-state
#   occupancy>0, (4) endpoint sanity (state-0 pertE ~191 plateau + λ=0.5/W0=1
#   apex decoupled ~0). Escalation if any pair <0.3: (i) add W0/λ window →
#   (ii) extra α ramp → (iii) extra U0 descent. The LOAD-BEARING fix is the
#   per-direction free split (backward equilibrates from a dminus base); α/U0
#   softening is the cushion (Path: "the dminus base IS the fix").
#
# Free leg ONLY. Selected via ``--free-schedule densified38``. The BOUND leg
# always uses canonical22 (receptor holds the ligand → no crossover gap).
# ---------------------------------------------------------------------------
# Forward 9-window ladder (states 10..18) as per-state tuples
# (λ1, λ2, W0coeff, ALPHA, U0). Exact reproduction of Path's REVISED arrays
# (pathology L142-146, ladder slices). Linear forward states 0..9 are
# λ1=λ2=λ, W0=0, ALPHA=0.10, U0=110.0, INTERMEDIATE=0 (same as canonical).
DENSE38_LAMBDA_FWD_LINEAR = [
    0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
]
DENSE38_LADDER = [
    # (lambda1, lambda2, w0coeff, alpha, u0)
    (0.45,  0.46,  0.15, 0.12, 105.0),
    (0.45,  0.47,  0.30, 0.14, 100.0),
    (0.46,  0.48,  0.45, 0.16,  95.0),
    (0.47,  0.49,  0.60, 0.18,  92.0),
    (0.48,  0.495, 0.72, 0.20,  90.0),
    (0.485, 0.498, 0.82, 0.22,  88.0),
    (0.49,  0.499, 0.90, 0.23,  86.0),
    (0.495, 0.5,   0.96, 0.24,  84.0),
    (0.5,   0.5,   1.00, 0.25,  82.0),
]


def _build_densified38_arrays(linear_lambdas: List[float]) -> Dict[str, List[Any]]:
    """Shared 38-state densified free-leg builder core.

    Both ``densified38`` and ``densified38v2`` use this single builder so the
    soft-core ladder (DENSE38_LADDER, 9 windows), the ALPHA ramp (0.10→0.25),
    the U0 ramp (110→82), the W0 ramp (0→1.0), the INTERMEDIATE pattern, and
    the whole-tuple-reverse backward construction stay byte-identical between
    the two schedules — ONLY the forward LINEAR λ array (``linear_lambdas``)
    differs. This avoids a divergent copy-paste implementation:
    densified38v2 is densified38 with a rebalanced linear segment, nothing else.

    Forward half = ``len(linear_lambdas)`` linear states (λ1=λ2=λ, W0=0,
    ALPHA=0.10, U0=110, INTERMEDIATE=0) + the 9-window DENSE38_LADDER
    (INTERMEDIATE=1). Backward half = whole-tuple reverse of forward half with
    DIRECTION=-1 (λ1 stays λ1, λ2 stays λ2 — NOT a λ1↔λ2 swap; all per-state
    arrays incl. ALPHA + U0 reversed so the ladder is symmetric). Total = 38
    states for a 10-state linear segment (10 linear + 9 ladder = 19/direction).

    Returns a dict with keys ``lambdas_1`` / ``lambdas_2`` / ``lambdas`` /
    ``directions`` / ``intermd`` / ``w0`` / ``alpha`` / ``u0`` / ``n_states``.
    The ``lambdas`` column tracks λ2 in the ladder (matching Path's LAMBDAS
    array, which equals λ2 in the soft-core windows).
    """
    n_linear = len(linear_lambdas)
    fwd_l1 = list(linear_lambdas) + [t[0] for t in DENSE38_LADDER]
    fwd_l2 = list(linear_lambdas) + [t[1] for t in DENSE38_LADDER]
    fwd_w0 = [0.0] * n_linear + [t[2] for t in DENSE38_LADDER]
    fwd_alpha = [0.10] * n_linear + [t[3] for t in DENSE38_LADDER]
    fwd_u0 = [110.0] * n_linear + [t[4] for t in DENSE38_LADDER]
    fwd_int = [0] * n_linear + [1] * len(DENSE38_LADDER)
    fwd_lam = list(linear_lambdas) + [t[1] for t in DENSE38_LADDER]

    # Backward half = reverse the forward half as whole per-state tuples (λ1
    # stays λ1, λ2 stays λ2; NOT a λ1↔λ2 swap). DIRECTION = -1. ALL per-state
    # arrays (incl. ALPHA + U0) are reversed so the ladder is symmetric.
    bwd_l1 = list(reversed(fwd_l1))
    bwd_l2 = list(reversed(fwd_l2))
    bwd_w0 = list(reversed(fwd_w0))
    bwd_alpha = list(reversed(fwd_alpha))
    bwd_u0 = list(reversed(fwd_u0))
    bwd_int = list(reversed(fwd_int))
    bwd_lam = list(reversed(fwd_lam))

    fwd_n = n_linear + len(DENSE38_LADDER)
    directions = [1] * fwd_n + [-1] * fwd_n
    return {
        "lambdas_1": fwd_l1 + bwd_l1,
        "lambdas_2": fwd_l2 + bwd_l2,
        "lambdas": fwd_lam + bwd_lam,
        "directions": directions,
        "intermd": fwd_int + bwd_int,
        "w0": fwd_w0 + bwd_w0,
        "alpha": fwd_alpha + bwd_alpha,
        "u0": fwd_u0 + bwd_u0,
        "n_states": len(directions),
    }


def _build_densified38_schedule() -> Dict[str, List[Any]]:
    """Construct the REVISED 38-state densified free-leg schedule arrays.

    Forward half (19) = 10 linear + 9-window anneal ladder; backward half (19)
    = whole-tuple reverse of forward half with DIRECTION=-1. Per-state ALPHA
    (0.10 linear → 0.25 apex) and U0 (110 linear → 82 apex) ramps soften the
    backward W0-peak (revised ladder fix; endpoint-invariant).

    Delegates to the shared ``_build_densified38_arrays`` core with the
    uniform-Δλ=0.05 linear array (DENSE38_LAMBDA_FWD_LINEAR).
    """
    return _build_densified38_arrays(DENSE38_LAMBDA_FWD_LINEAR)


_DENSE38 = _build_densified38_schedule()
DENSE38_LAMBDAS_1 = _DENSE38["lambdas_1"]
DENSE38_LAMBDAS_2 = _DENSE38["lambdas_2"]
DENSE38_LAMBDAS = _DENSE38["lambdas"]
DENSE38_DIRECTIONS = _DENSE38["directions"]
DENSE38_INTERMD = _DENSE38["intermd"]
DENSE38_W0 = _DENSE38["w0"]
DENSE38_ALPHA = _DENSE38["alpha"]
DENSE38_U0 = _DENSE38["u0"]
DENSE38_N_STATES = _DENSE38["n_states"]
assert DENSE38_N_STATES == 38, "REVISED densified schedule must be 38 states (19+19)"
assert sum(DENSE38_INTERMD) == 18, "REVISED densified schedule must have 18 INTERMEDIATE windows"


# ---------------------------------------------------------------------------
# densified38v2 — count-neutral λ-REBALANCE of densified38's LINEAR segment
# (dplus 8→9 overlap fix, 2026-06-05).
#
# WHY: the densified38 free pilot passed NaN-fix + occupancy, but the dplus
# state 8→9 pair (λ=0.40→0.45, both W0=0 / INTERMEDIATE=0 — the PLAIN-LINEAR
# region BEFORE the soft-core ladder, which starts at state 10) had a genuine
# λ-spacing overlap gap (Bhattacharyya BC=0.131, below the well-determined
# floor). Diagnosis: states 0-8 (λ=0.00→0.40) are an OVERCROWDED plateau
# (pertE ~190±3-8, adjacent BC 0.91-0.99 — overlap-saturated) while the entire
# coupled→decoupled flip is crammed into the SINGLE step λ=0.40→0.45. The fix is
# a PURE LINEAR-λ rebalance: thin the saturated plateau (8 windows over
# 0.00-0.35 → 5 windows) and densify the turnover (Δλ=0.02 across 0.40→0.45).
#
# This is RIGOR-NEUTRAL + ranking-safe (ΔG-unbias proof): free energy is
# a state function (Zwanzig 1954 / Kirkwood 1935); both physical endpoints
# (λ=0.00 coupled, λ=0.5/W0=1.0 apex) are LEFT INVARIANT, so relocating interior
# λ states changes ONLY the estimator's overlap/variance, never E[ΔG]. The fix
# touches NEITHER the soft-core ladder (DENSE38_LADDER, 9 windows — byte-
# identical) NOR the α/U0/W0 ramps (those live in the soft-core region where
# α/U0/W0 are inert in the linear segment by construction).
#
# It is applied to the FORWARD linear array ONLY: the backward half is the
# whole-tuple reverse of the forward half (_build_densified38_arrays), so dminus
# AUTOMATICALLY inherits the same rebalanced linear spacing, mirror-imaged
# (Q4 — symmetry is automatic + beneficial; dminus's marginal 9↔10 turnover
# BC=0.364 can only IMPROVE with finer spacing). NO separate backward edit.
#
# densified38 is NOT overwritten — it remains selectable for forensic
# reproducibility. densified38v2 is ADDED alongside it.
#
# Forward linear segment (10 states): 5 thinned plateau (was 8 over 0.00-0.35)
# + 5 turnover at Δλ=0.02 across the 0.40→0.45 flip. λ=0.00 coupled endpoint
# preserved as the first state. W0=0, ALPHA=0.10, U0=110, INTERMEDIATE=0 for
# all 10 (same as densified38's linear segment). 10 linear + 9 ladder = 19
# states/direction = 38 total. Free leg ONLY (bound leg always canonical22).
# ---------------------------------------------------------------------------
DENSE38v2_LAMBDA_FWD_LINEAR = [
    0.00, 0.08, 0.16, 0.24, 0.32,   # 5 plateau (saturated overlap, thinned)
    0.38, 0.40, 0.42, 0.44, 0.45,   # 5 turnover (Δλ=0.02 across the flip)
]


def _build_densified38v2_schedule() -> Dict[str, List[Any]]:
    """Construct the densified38v2 (REBALANCED) 38-state free-leg schedule.

    Identical to densified38 in EVERY array except the forward LINEAR λ values:
    the saturated coupled plateau is thinned and the 0.40→0.45 turnover is
    densified to Δλ=0.02. Delegates to the shared
    ``_build_densified38_arrays`` core with DENSE38v2_LAMBDA_FWD_LINEAR — so the
    9-window soft-core ladder, the ALPHA/U0/W0 ramps, the INTERMEDIATE pattern,
    and the whole-tuple-reverse backward construction stay byte-identical to
    densified38.
    """
    return _build_densified38_arrays(DENSE38v2_LAMBDA_FWD_LINEAR)


_DENSE38v2 = _build_densified38v2_schedule()
DENSE38v2_LAMBDAS_1 = _DENSE38v2["lambdas_1"]
DENSE38v2_LAMBDAS_2 = _DENSE38v2["lambdas_2"]
DENSE38v2_LAMBDAS = _DENSE38v2["lambdas"]
DENSE38v2_DIRECTIONS = _DENSE38v2["directions"]
DENSE38v2_INTERMD = _DENSE38v2["intermd"]
DENSE38v2_W0 = _DENSE38v2["w0"]
DENSE38v2_ALPHA = _DENSE38v2["alpha"]
DENSE38v2_U0 = _DENSE38v2["u0"]
DENSE38v2_N_STATES = _DENSE38v2["n_states"]
assert DENSE38v2_N_STATES == 38, "densified38v2 schedule must be 38 states (19+19)"
assert sum(DENSE38v2_INTERMD) == 18, "densified38v2 schedule must have 18 INTERMEDIATE windows"
assert len(DENSE38v2_LAMBDA_FWD_LINEAR) == 10, "densified38v2 linear segment must be 10 states"
assert DENSE38v2_LAMBDA_FWD_LINEAR[0] == 0.00, "densified38v2 coupled endpoint λ=0.00 must be first"
# Soft-core ladder block must be byte-identical to densified38 (only the linear
# segment differs). States 10-18 forward = the 9-window ladder.
assert DENSE38v2_W0[10:19] == DENSE38_W0[10:19], "densified38v2 ladder W0 must match densified38"
assert DENSE38v2_ALPHA[10:19] == DENSE38_ALPHA[10:19], "densified38v2 ladder ALPHA must match densified38"
assert DENSE38v2_U0[10:19] == DENSE38_U0[10:19], "densified38v2 ladder U0 must match densified38"


# ---------------------------------------------------------------------------
# densified38v3 — PER-DIRECTION (NON-mirror) schedule: dminus gets a W0-graded
# micro-bridge at its soft-core-end → plateau handoff
# (densified38v2 dminus-asymmetry fix, 2026-06-06).
#
# WHY (the empirical finding that BREAKS densified38v2's mirror symmetry):
#   densified38v2 = forward linear rebalance + ``reversed()`` whole-tuple mirror
#   for dminus. The v2 re-pilot showed the dplus 8→9 gap was FIXED (binning-free
#   BC 0.131 → 0.70, all dplus pairs ≥ 0.63) BUT the reversed mirror DEGRADED
#   the dminus 9→10 turnover (BC 0.364 → 0.243 Gauss / 0.37 KDE — a genuine gate
#   FAIL). Root cause (re-derived from raw col9 pertE): the soft-core
#   W0 ladder sits on OPPOSITE sides of the coupled→decoupled flip in the two
#   directions because DISPLACEMENT = +25 (dplus, coupled base) vs −25 (dminus,
#   decoupled base). In dplus the broad transition state (state 9, λ=0.45 W0=0,
#   mean pertE 78) hands off to the apex THROUGH the W0-graded ladder that
#   immediately follows it (state 10 W0=0.15) — the ladder BUFFERS the turnover.
#   In dminus (reversed) the ladder ENDS before the turnover, so the broad
#   transition state (state 9, λ=0.45 W0=0, mean pertE 112, sd 79) abuts the
#   BARE tight W0=0 plateau (state 10, λ=0.44 W0=0, mean pertE 177, sd 43)
#   DIRECTLY, with NO W0-graded intermediate → the 9→10 overlap hole. Note:
#   index reverse-symmetry ≠ overlap-symmetry; each direction's overlap
#   profile must be rebalanced from its OWN data, NOT inherited via reverse().
#
# THE FIX (per-direction, minimal):
#   * dplus = the CLEAN densified38v2 forward arrays, UNCHANGED (dplus already
#     passes every pair; reuses _build_densified38_arrays + DENSE38v2 linear +
#     the shared DENSE38_LADDER). dplus stays byte-identical to densified38v2's
#     dplus — verified by regression test.
#   * dminus = NOT reversed(dplus). It is reversed(dplus) PLUS one EXTRA
#     W0-graded micro-bridge window inserted at the soft-core-end → plateau
#     handoff (between the reversed-mirror broad transition state and the bare
#     tight plateau, i.e. the 9→10 hole). The bridge is a soft-core intermediate
#     so the broad transition hands off THROUGH a W0-graded step instead of
#     abutting the bare plateau. Micro-bridge spec: λ1=λ2=0.45 (the broad
#     transition's λ), W0 ≈ 0.07, with ALPHA/U0 monotone-interpolated between the
#     ladder-end neighbor (W0=0.15: α=0.12, U0=105) and the plateau neighbor
#     (W0=0: α=0.10, U0=110). At W0=0.07 (≈ midway of 0.15) α≈0.11, U0≈107.5.
#     INTERMEDIATE=1 (it is a soft-core / W0>0 window).
#
# COUNT (pre-registered, G7): UNEQUAL — dplus 19, dminus 20. The
#   bridge adds ONE state to dminus only. This is SAFE because EVERY count-
#   dependent code path is DIRECTION-derived + count-agnostic (audited 2026-06-06,
#   see G7 notes): write_cntl_file emits per-state arrays
#   generically (_csv); the per-direction launcher derives slices from the cntl
#   DIRECTION column (_derive_state_counts_from_directions: total=len, fwd=Σ(d≥0),
#   bwd=total−fwd — NO len(dplus)==len(dminus) / zip / hardcoded-19 assumption);
#   generate_per_direction_cntls slices dplus=(0,fwd) / dminus=(fwd,total);
#   stage_per_direction_subdir creates bwd=total−fwd replica dirs;
#   merge_per_direction_outputs offsets dminus by fwd ("computed explicitly so an
#   asymmetric schedule would still map correctly"); the UWHAM post-processor
#   (_calculate_uwham_multi_intermediate) solves each leg SEPARATELY from
#   where(DIRECTION≥0)/where(DIRECTION<0) and combines at the SCALAR level — the
#   two legs are NEVER index-aligned. So unequal counts are the cleaner path: the
#   dplus arrays stay byte-identical to v2 (no redundant no-op window injected to
#   force equality). Endpoints (λ=0 coupled first, W0=1.0/λ=0.5 apex) are IMMUTABLE
#   in BOTH directions → ΔG-unbiased / ranking-safe (Zwanzig 1954 / Kirkwood
#   1935). densified38 / densified38v2 are NOT overwritten.
#
# Free leg ONLY. The BOUND leg always uses canonical22.
# ---------------------------------------------------------------------------
# dminus W0-graded micro-bridge (one per-state tuple, in dminus-local order:
# inserted between the reversed-mirror broad transition state and the bare tight
# plateau). (λ1, λ2, W0coeff, ALPHA, U0). λ matches the broad transition (0.45);
# W0=0.07; α/U0 monotone-interpolated between ladder-end (0.15/0.12/105) and
# plateau (0.0/0.10/110).
DENSE38v3_DMINUS_BRIDGE = (0.45, 0.45, 0.07, 0.11, 107.5)


def _build_densified38v3_arrays() -> Dict[str, List[Any]]:
    """Construct the densified38v3 PER-DIRECTION (NON-mirror) free-leg arrays.

    dplus (forward, 19 states) = the CLEAN densified38v2 forward half, reused
    verbatim from ``_build_densified38_arrays(DENSE38v2_LAMBDA_FWD_LINEAR)`` (the
    forward slice ``[:19]``). dplus is already overlap-balanced (re-pilot: all
    pairs ≥ 0.63) — it is NOT touched here.

    dminus (backward, 20 states) = the densified38v2 backward half (the whole-
    tuple ``reversed()`` of the forward half — 19 states) PLUS one extra
    W0-graded micro-bridge window (``DENSE38v3_DMINUS_BRIDGE``) inserted at the
    soft-core-end → plateau handoff. DIRECTION = -1 for all dminus states.

    This BREAKS densified38v2's forward-only-mirror symmetry deliberately
    (reverse-symmetric INDEX construction does not produce symmetric
    OVERLAP when the soft-core ladder sits on opposite sides of the ±DISPLACEMENT
    flip). dplus and dminus therefore have UNEQUAL per-direction state counts
    (19 vs 20) — SAFE because every downstream code path is DIRECTION-derived /
    count-agnostic (see module block comment, G7 audit).

    Bridge placement: the dminus backward half is
    ``reversed(forward) = [reversed 9-ladder | reversed 10-linear]``. In the
    re-pilot .out data this is dminus-local states 0-8 (reversed ladder, W0
    1.0→0.15) then 9-18 (reversed linear, W0=0). The broad transition is at
    dminus-local state 9 (λ=0.45, W0=0) and the bare tight plateau begins at
    state 10 (λ=0.44, W0=0); the gate-failing bridge is 9→10. The micro-bridge
    is inserted BETWEEN local state 9 and local state 10 so the broad transition
    hands off to the plateau through a W0=0.07 soft-core intermediate. This makes
    dminus 20 states; the insert index is the count of reversed-ladder states
    (= ``len(DENSE38_LADDER)``) + the count of reversed-linear states up to and
    INCLUDING the broad transition state (the last W0=0 state at λ=0.45). Because
    the reversed-linear segment begins at λ=0.45 (the forward linear's last
    value), the broad transition is the FIRST reversed-linear state, so the
    insert index = ``len(DENSE38_LADDER) + 1``.

    Returns a dict with the same keys as ``_build_densified38_arrays``
    (``lambdas_1`` / ``lambdas_2`` / ``lambdas`` / ``directions`` / ``intermd``
    / ``w0`` / ``alpha`` / ``u0`` / ``n_states``). n_states = 39 (19 + 20).
    """
    # dplus = clean densified38v2 forward half (states [:19]).
    v2 = _build_densified38_arrays(DENSE38v2_LAMBDA_FWD_LINEAR)
    n_fwd = len(DENSE38v2_LAMBDA_FWD_LINEAR) + len(DENSE38_LADDER)  # 10 + 9 = 19
    fwd = {k: list(v2[k][:n_fwd]) for k in (
        "lambdas_1", "lambdas_2", "lambdas", "intermd", "w0", "alpha", "u0")}

    # dminus base = the v2 backward half (whole-tuple reverse of forward), states
    # [n_fwd:]. This is 19 states; the bridge insertion makes it 20.
    bwd = {k: list(v2[k][n_fwd:]) for k in (
        "lambdas_1", "lambdas_2", "lambdas", "intermd", "w0", "alpha", "u0")}

    # Insert index: after the reversed ladder (len DENSE38_LADDER) + the broad
    # transition state (the first reversed-linear state at λ=0.45). The bridge
    # sits BETWEEN that broad transition state and the bare tight plateau that
    # immediately follows it (the 9→10 overlap hole).
    insert_idx = len(DENSE38_LADDER) + 1  # = 10 (between dminus-local 9 and 10)
    bl1, bl2, bw0, balpha, bu0 = DENSE38v3_DMINUS_BRIDGE
    bridge = {
        "lambdas_1": bl1,
        "lambdas_2": bl2,
        "lambdas": bl2,        # LAMBDAS tracks λ2 in soft-core windows
        "intermd": 1,          # soft-core / W0>0 → INTERMEDIATE
        "w0": bw0,
        "alpha": balpha,
        "u0": bu0,
    }
    for key in bwd:
        bwd[key].insert(insert_idx, bridge[key])

    n_bwd = n_fwd + 1  # 20
    directions = [1] * n_fwd + [-1] * n_bwd
    return {
        "lambdas_1": fwd["lambdas_1"] + bwd["lambdas_1"],
        "lambdas_2": fwd["lambdas_2"] + bwd["lambdas_2"],
        "lambdas": fwd["lambdas"] + bwd["lambdas"],
        "directions": directions,
        "intermd": fwd["intermd"] + bwd["intermd"],
        "w0": fwd["w0"] + bwd["w0"],
        "alpha": fwd["alpha"] + bwd["alpha"],
        "u0": fwd["u0"] + bwd["u0"],
        "n_states": len(directions),
    }


def _build_densified38v3_schedule() -> Dict[str, List[Any]]:
    """Construct the densified38v3 (PER-DIRECTION, NON-mirror) free-leg schedule.

    dplus = clean densified38v2 forward half (19 states); dminus = densified38v2
    backward half (19 states) + one W0-graded micro-bridge at the soft-core-end
    → plateau handoff (= 20 states). Total = 39 states (UNEQUAL per-direction).
    Delegates the heavy lifting to ``_build_densified38v3_arrays`` (which reuses
    the shared ``_build_densified38_arrays`` core for the dplus side + the shared
    DENSE38_LADDER — no divergent copy-paste ladder).
    """
    return _build_densified38v3_arrays()


_DENSE38v3 = _build_densified38v3_schedule()
DENSE38v3_LAMBDAS_1 = _DENSE38v3["lambdas_1"]
DENSE38v3_LAMBDAS_2 = _DENSE38v3["lambdas_2"]
DENSE38v3_LAMBDAS = _DENSE38v3["lambdas"]
DENSE38v3_DIRECTIONS = _DENSE38v3["directions"]
DENSE38v3_INTERMD = _DENSE38v3["intermd"]
DENSE38v3_W0 = _DENSE38v3["w0"]
DENSE38v3_ALPHA = _DENSE38v3["alpha"]
DENSE38v3_U0 = _DENSE38v3["u0"]
DENSE38v3_N_STATES = _DENSE38v3["n_states"]
# Per-direction counts (UNEQUAL: dplus 19, dminus 20 — pre-registered).
DENSE38v3_N_FWD = sum(1 for d in DENSE38v3_DIRECTIONS if d >= 0)
DENSE38v3_N_BWD = DENSE38v3_N_STATES - DENSE38v3_N_FWD
assert DENSE38v3_N_STATES == 39, "densified38v3 must be 39 states (19 dplus + 20 dminus)"
assert DENSE38v3_N_FWD == 19, "densified38v3 dplus (forward) must be 19 states"
assert DENSE38v3_N_BWD == 20, "densified38v3 dminus (backward) must be 20 states (+1 bridge)"
# DIRECTION column is a contiguous +1 block then -1 block (required by the
# per-direction launcher's _derive_state_counts_from_directions contiguity gate).
assert DENSE38v3_DIRECTIONS[:DENSE38v3_N_FWD] == [1] * DENSE38v3_N_FWD
assert DENSE38v3_DIRECTIONS[DENSE38v3_N_FWD:] == [-1] * DENSE38v3_N_BWD
# dplus (forward [:19]) MUST be byte-identical to densified38v2's dplus — the
# clean side is reused verbatim, NOT rebuilt (do not change dplus).
for _k in ("lambdas_1", "lambdas_2", "lambdas", "intermd", "w0", "alpha", "u0"):
    assert _DENSE38v3[_k][:DENSE38v3_N_FWD] == _DENSE38v2[_k][:19], (
        f"densified38v3 dplus {_k} must equal densified38v2 dplus"
    )
# Endpoints IMMUTABLE both directions (ΔG-unbias / ranking-only):
#   dplus first state = coupled λ=0.00, W0=0; dplus last fwd state = W0=1.0/λ=0.5 apex.
assert DENSE38v3_LAMBDAS_1[0] == 0.00 and DENSE38v3_W0[0] == 0.0, "dplus coupled endpoint"
assert DENSE38v3_W0[DENSE38v3_N_FWD - 1] == 1.0 and DENSE38v3_LAMBDAS_2[DENSE38v3_N_FWD - 1] == 0.5, "dplus apex"
#   dminus first state = W0=1.0/λ=0.5 apex (reversed); dminus last state = coupled λ=0.00.
assert DENSE38v3_W0[DENSE38v3_N_FWD] == 1.0 and DENSE38v3_LAMBDAS_2[DENSE38v3_N_FWD] == 0.5, "dminus apex (first bwd state)"
assert DENSE38v3_LAMBDAS_1[-1] == 0.00 and DENSE38v3_W0[-1] == 0.0, "dminus coupled endpoint (last state)"
# The dminus bridge is the ONLY W0>0 window NOT present in densified38v2's
# backward half (the extra state). Confirm exactly one extra INTERMEDIATE in
# dminus vs the v2 dminus (the bridge carries INTERMEDIATE=1).
assert sum(int(x) for x in DENSE38v3_INTERMD[DENSE38v3_N_FWD:]) == (
    sum(int(x) for x in _DENSE38v2["intermd"][19:]) + 1
), "densified38v3 dminus must have exactly one extra INTERMEDIATE (the bridge)"


# ---------------------------------------------------------------------------
# densified38v4 — PER-DIRECTION (NON-mirror): dminus gets a SECOND W0-graded
# micro-bridge at the (re-indexed) soft-core-ladder-end → plateau handoff
# (dminus second-bridge fix, 2026-06-06).
#
# WHY (the convergent — NOT whack-a-mole — finding):
#   The v3 re-pilot CLOSED the dminus 9→10 hole the reversed mirror created
#   (binning-free BC 0.24 → 0.84). dplus stayed clean (min Gauss 0.443). BUT the
#   v3 dminus exposed ONE more broad(sd>40)→tight(sd<10) abutment, now at the
#   v3-dminus-local 5→6 handoff (re-pilot col9 pertE, warmup-5, 20 replicas:
#   state 5 mean 79.05 / sd 75.30 / p95 196 — BROAD + bimodal, ~30% of samples
#   reach the coupled needle; state 6 mean 190.98 / sd 3.85 — TIGHT plateau;
#   Gauss-BC 0.184 / KDE-BC 0.242, ABOVE the MBAR O≥0.03 estimability floor but
#   BELOW the 0.3 soft target).
#
#   DECISIVE finding: this is NOT a second pathology. The dminus
#   leg has EXACTLY ONE coupled→decoupled turnover, and it is the SAME single
#   turnover as v2/v3 — it RELOCATED by index when the v3 bridge re-indexed the
#   W0 ladder; it did not multiply. The leg-wide broad(sd>40) state count is
#   SHRINKING (v2 {7,8,9}=3 → v3 {5}=1) and each bridge RAISES its target's
#   overlap (9→10: 0.24→0.84). This is the divergence-detector's CONVERGENT
#   signature, not barrier-limited whack-a-mole. One more graded-W0 bridge has
#   ~85% chance of closure (≤2 worst case, finite/bounded).
#
#   ITERATION CAP (the anti-whack-a-mole rule): at most 2 targeted
#   micro-bridge pilots per leg. v4 is the 2nd (count→2). After the v4 re-pilot,
#   regardless of result, STOP pilot-iteration and LAUNCH (with a documented
#   soft-flag if 5→6 is still 0.1–0.3). NO densified38v5.
#
# THE FIX (per-direction, minimal):
#   * dplus = the CLEAN densified38v2/v3 forward arrays, UNCHANGED (byte-equal —
#     dplus passes every pair; reuses _build_densified38v3_arrays' dplus side,
#     which itself reuses _build_densified38_arrays + DENSE38v2 linear). dplus
#     stays byte-identical to densified38v2's AND densified38v3's dplus.
#   * dminus = the densified38v3 dminus (20 states — which ALREADY has the
#     9→10 W0=0.07 bridge) PLUS one MORE W0-graded micro-bridge inserted at the
#     5→6 handoff (between the v3-dminus-local broad transition state 5 —
#     λ1=0.470/λ2=0.490, W0=0.600 — and the tight plateau side state 6 —
#     λ1=0.460/λ2=0.480, W0=0.450). Micro-bridge spec: λ1=0.465, λ2=0.485
#     (interior to the two neighbors' λ), W0=0.17, with ALPHA/U0 monotone-
#     interpolated between state 5 (α=0.18, U0=92) and state 6 (α=0.16, U0=95):
#     α=0.17, U0=93.5 (the midpoints). INTERMEDIATE=1 (soft-core / W0>0 window).
#     This puts a graded intermediate exactly at the soft-core-ladder-end →
#     plateau handoff that currently has none. Construction is IDENTICAL in kind
#     to the v3 9→10 bridge that demonstrably worked (no divergent ladder).
#
# COUNT (pre-registered, G7 re-audit): UNEQUAL — dplus 19,
#   dminus 21 (v3 dminus 20 + this bridge), total 40. SAFE because EVERY count-
#   dependent path is DIRECTION-derived + count-agnostic (re-audited 2026-06-06,
#   same 6 sites verified for v3 — see G7 notes):
#   _derive_state_counts_from_directions (total=len, fwd=Σ(d≥0), bwd=total−fwd —
#   NO len(dplus)==len(dminus) / zip / hardcoded count); generate_per_direction_
#   cntls (slices dplus=(0,fwd) / dminus=(fwd,total)); stage_per_direction_subdir
#   (bwd=total−fwd replica dirs); merge_per_direction_outputs (dminus offset=fwd,
#   "computed explicitly so an asymmetric schedule would still map correctly");
#   _calculate_uwham_multi_intermediate (legs solved separately from
#   where(DIRECTION≥0)/where(DIRECTION<0), combined at the SCALAR level — never
#   index-aligned); write_cntl_file (generic _csv serialization). 40=19/21 flows
#   identically to 39=19/20. Endpoints (λ=0 coupled first, W0=1.0/λ=0.5 apex) are
#   IMMUTABLE both directions → ΔG-unbiased / ranking-safe (Zwanzig 1954 /
#   Kirkwood 1935). densified38 / v2 / v3 / canonical22 NOT overwritten.
#
# Free leg ONLY. The BOUND leg always uses canonical22.
# ---------------------------------------------------------------------------
# dminus 5→6 W0-graded micro-bridge (one per-state tuple; (λ1, λ2, W0coeff,
# ALPHA, U0)). λ1=0.465, λ2=0.485 (interior to the v3-dminus-local 5 / 6
# neighbors); W0=0.17; α/U0 = midpoints of state 5 (0.18/92) and state 6
# (0.16/95).
DENSE38v4_DMINUS_BRIDGE_56 = (0.465, 0.485, 0.17, 0.17, 93.5)


def _build_densified38v4_arrays() -> Dict[str, List[Any]]:
    """Construct the densified38v4 PER-DIRECTION (NON-mirror) free-leg arrays.

    dplus (forward, 19 states) = the CLEAN densified38v2/v3 forward half, reused
    verbatim from ``_build_densified38v3_arrays`` (whose dplus side is itself the
    densified38v2 forward half). dplus is overlap-balanced (re-pilot: all pairs
    ≥ 0.443 Gauss) — it is NOT touched here, and stays byte-identical to v2/v3.

    dminus (backward, 21 states) = the densified38v3 dminus (20 states — which
    ALREADY carries the 9→10 W0=0.07 micro-bridge) PLUS one MORE W0-graded
    micro-bridge (``DENSE38v4_DMINUS_BRIDGE_56``) inserted at the v3-dminus-local
    5→6 soft-core-ladder-end → plateau handoff. DIRECTION = -1 for all dminus.

    Bridge placement (robust — located by NEIGHBOR signature, NOT a magic index):
    the v3 dminus-local state 5 is the broad transition (λ1=0.470, λ2=0.490,
    W0=0.600) and state 6 is the tight plateau side (λ1=0.460, λ2=0.480,
    W0=0.450). The 5→6 abutment is the single broad↔tight handoff counted
    in the v3 re-pilot. We find the index of that exact (l1, l2, W0) state-5
    signature in the v3 dminus block and insert the new bridge immediately AFTER
    it (between state 5 and state 6), so the broad transition hands off to the
    plateau through a W0=0.17 soft-core intermediate. This reuses the v3 dminus
    array as the base (no divergent re-derivation of the v3 9→10 bridge).

    Returns a dict with the same keys as ``_build_densified38_arrays``. n_states
    = 40 (19 dplus + 21 dminus, UNEQUAL).
    """
    # Base = the full v3 per-direction arrays (dplus 19 + dminus 20).
    v3 = _build_densified38v3_arrays()
    n_fwd = DENSE38v3_N_FWD  # 19 (forward, contiguous +1 block)
    keys = ("lambdas_1", "lambdas_2", "lambdas", "intermd", "w0", "alpha", "u0")
    fwd = {k: list(v3[k][:n_fwd]) for k in keys}          # dplus, UNCHANGED
    bwd = {k: list(v3[k][n_fwd:]) for k in keys}          # v3 dminus (20 states)

    # Locate the 5→6 handoff in the v3 dminus block by the state-5 broad-
    # transition signature (λ1=0.470, λ2=0.490, W0=0.600). Insert the new bridge
    # immediately AFTER state 5 (i.e. between state 5 and state 6). Fail loud if
    # the signature is absent or non-unique (guards a future v3 array change from
    # silently mis-placing the bridge — source-verification posture).
    bl1, bl2, bw0, balpha, bu0 = DENSE38v4_DMINUS_BRIDGE_56
    s5_l1, s5_l2, s5_w0 = 0.470, 0.490, 0.600
    s5_positions = [
        j for j in range(len(bwd["w0"]))
        if abs(bwd["lambdas_1"][j] - s5_l1) < 1e-9
        and abs(bwd["lambdas_2"][j] - s5_l2) < 1e-9
        and abs(bwd["w0"][j] - s5_w0) < 1e-9
    ]
    if len(s5_positions) != 1:
        raise ValueError(
            "densified38v4: expected exactly ONE v3-dminus broad-transition "
            "state matching the 5→6 handoff signature "
            f"(λ1={s5_l1}, λ2={s5_l2}, W0={s5_w0}); found {len(s5_positions)} "
            f"at {s5_positions!r}. The v3 dminus array layout changed — re-derive "
            "the bridge placement against the v3 re-pilot .out before proceeding."
        )
    insert_idx = s5_positions[0] + 1  # immediately after state 5 → between 5 and 6
    bridge = {
        "lambdas_1": bl1,
        "lambdas_2": bl2,
        "lambdas": bl2,        # LAMBDAS tracks λ2 in soft-core windows
        "intermd": 1,          # soft-core / W0>0 → INTERMEDIATE
        "w0": bw0,
        "alpha": balpha,
        "u0": bu0,
    }
    for key in keys:
        bwd[key].insert(insert_idx, bridge[key])

    n_bwd = len(bwd["w0"])  # 21 (v3 dminus 20 + this bridge)
    directions = [1] * n_fwd + [-1] * n_bwd
    return {
        "lambdas_1": fwd["lambdas_1"] + bwd["lambdas_1"],
        "lambdas_2": fwd["lambdas_2"] + bwd["lambdas_2"],
        "lambdas": fwd["lambdas"] + bwd["lambdas"],
        "directions": directions,
        "intermd": fwd["intermd"] + bwd["intermd"],
        "w0": fwd["w0"] + bwd["w0"],
        "alpha": fwd["alpha"] + bwd["alpha"],
        "u0": fwd["u0"] + bwd["u0"],
        "n_states": len(directions),
    }


def _build_densified38v4_schedule() -> Dict[str, List[Any]]:
    """Construct the densified38v4 (PER-DIRECTION, NON-mirror) free-leg schedule.

    dplus = clean densified38v2/v3 forward half (19 states, byte-equal); dminus =
    densified38v3 dminus (20 states, incl. the 9→10 W0=0.07 bridge) + one MORE
    W0-graded micro-bridge at the 5→6 handoff (= 21 states). Total = 40 states
    (UNEQUAL per-direction). Delegates to ``_build_densified38v4_arrays`` (which
    reuses ``_build_densified38v3_arrays`` as its base — no divergent re-build of
    the dplus side or the v3 9→10 bridge).
    """
    return _build_densified38v4_arrays()


_DENSE38v4 = _build_densified38v4_schedule()
DENSE38v4_LAMBDAS_1 = _DENSE38v4["lambdas_1"]
DENSE38v4_LAMBDAS_2 = _DENSE38v4["lambdas_2"]
DENSE38v4_LAMBDAS = _DENSE38v4["lambdas"]
DENSE38v4_DIRECTIONS = _DENSE38v4["directions"]
DENSE38v4_INTERMD = _DENSE38v4["intermd"]
DENSE38v4_W0 = _DENSE38v4["w0"]
DENSE38v4_ALPHA = _DENSE38v4["alpha"]
DENSE38v4_U0 = _DENSE38v4["u0"]
DENSE38v4_N_STATES = _DENSE38v4["n_states"]
# Per-direction counts (UNEQUAL: dplus 19, dminus 21 — pre-registered).
DENSE38v4_N_FWD = sum(1 for d in DENSE38v4_DIRECTIONS if d >= 0)
DENSE38v4_N_BWD = DENSE38v4_N_STATES - DENSE38v4_N_FWD
assert DENSE38v4_N_STATES == 40, "densified38v4 must be 40 states (19 dplus + 21 dminus)"
assert DENSE38v4_N_FWD == 19, "densified38v4 dplus (forward) must be 19 states"
assert DENSE38v4_N_BWD == 21, "densified38v4 dminus (backward) must be 21 states (v3 20 + 1 bridge)"
# DIRECTION column is a contiguous +1 block then -1 block (required by the
# per-direction launcher's _derive_state_counts_from_directions contiguity gate).
assert DENSE38v4_DIRECTIONS[:DENSE38v4_N_FWD] == [1] * DENSE38v4_N_FWD
assert DENSE38v4_DIRECTIONS[DENSE38v4_N_FWD:] == [-1] * DENSE38v4_N_BWD
# dplus (forward [:19]) MUST be byte-identical to densified38v2's AND densified38v3's
# dplus — the clean side is reused verbatim, NOT rebuilt (dplus is
# clean, do not touch it").
for _k in ("lambdas_1", "lambdas_2", "lambdas", "intermd", "w0", "alpha", "u0"):
    assert _DENSE38v4[_k][:DENSE38v4_N_FWD] == _DENSE38v2[_k][:19], (
        f"densified38v4 dplus {_k} must equal densified38v2 dplus"
    )
    assert _DENSE38v4[_k][:DENSE38v4_N_FWD] == _DENSE38v3[_k][:19], (
        f"densified38v4 dplus {_k} must equal densified38v3 dplus"
    )
# Endpoints IMMUTABLE both directions (ΔG-unbias / ranking-only):
#   dplus first state = coupled λ=0.00, W0=0; dplus last fwd state = W0=1.0/λ=0.5 apex.
assert DENSE38v4_LAMBDAS_1[0] == 0.00 and DENSE38v4_W0[0] == 0.0, "dplus coupled endpoint"
assert DENSE38v4_W0[DENSE38v4_N_FWD - 1] == 1.0 and DENSE38v4_LAMBDAS_2[DENSE38v4_N_FWD - 1] == 0.5, "dplus apex"
#   dminus first state = W0=1.0/λ=0.5 apex (reversed); dminus last state = coupled λ=0.00.
assert DENSE38v4_W0[DENSE38v4_N_FWD] == 1.0 and DENSE38v4_LAMBDAS_2[DENSE38v4_N_FWD] == 0.5, "dminus apex (first bwd state)"
assert DENSE38v4_LAMBDAS_1[-1] == 0.00 and DENSE38v4_W0[-1] == 0.0, "dminus coupled endpoint (last state)"
# dminus carries EXACTLY one more INTERMEDIATE than v3 dminus (the new 5→6 bridge,
# INTERMEDIATE=1). v3 dminus already had +1 over v2 (the 9→10 bridge); so v4
# dminus = v2 dminus + 2 bridges.
assert sum(int(x) for x in DENSE38v4_INTERMD[DENSE38v4_N_FWD:]) == (
    sum(int(x) for x in _DENSE38v3["intermd"][DENSE38v3_N_FWD:]) + 1
), "densified38v4 dminus must have exactly one more INTERMEDIATE than v3 dminus"
# Both bridges present in v4 dminus: the v3 9→10 W0=0.07 AND the new 5→6 W0=0.17.
assert sum(1 for w in DENSE38v4_W0[DENSE38v4_N_FWD:] if abs(w - 0.07) < 1e-9) == 1, (
    "densified38v4 dminus must keep the v3 9→10 W0=0.07 bridge"
)
assert sum(1 for w in DENSE38v4_W0[DENSE38v4_N_FWD:] if abs(w - 0.17) < 1e-9) == 1, (
    "densified38v4 dminus must add the new 5→6 W0=0.17 bridge"
)


# ---------------------------------------------------------------------------
# Schedule registry — selectable via ``--free-schedule`` (default canonical22).
# Each entry is a dict the cntl writer + run_metadata read uniformly. The
# 22-state default preserves the existing per-direction split path (bound leg
# stays 22) and all existing regression tests (the 22-state path is NOT
# overwritten, the densified ladder is ADDED).
# ---------------------------------------------------------------------------
def _schedule_dict(lambdas_1, lambdas_2, directions, intermd, w0, alpha, u0,
                   lambdas=None) -> Dict[str, Any]:
    """Pack a per-state schedule into the dict shape ``write_cntl_file``
    + ``run_metadata`` consume. ``lambdas`` (the LAMBDAS cntl column)
    defaults to ``lambdas_1`` when not supplied (true for the canonical
    22-state schedule where LAMBDAS == LAMBDA1 == LAMBDA2).
    """
    return {
        "lambdas": list(lambdas if lambdas is not None else lambdas_1),
        "lambdas_1": list(lambdas_1),
        "lambdas_2": list(lambdas_2),
        "directions": list(directions),
        "intermd": list(intermd),
        "w0": list(w0),
        "alpha": list(alpha),
        "u0": list(u0),
        "n_states": len(directions),
    }


CANONICAL22_SCHEDULE = _schedule_dict(
    lambdas_1=LAMBDAS_1, lambdas_2=LAMBDAS_2, directions=DIRECTIONS,
    intermd=INTERMD, w0=W0, alpha=ALPHA, u0=U0, lambdas=LAMBDAS_1,
)
# DEPRECATED — Factor-B broken (NaN'd at backward W0-peak even with full
# structprep). Retained for forensic reproduction of the 2026-06-05 smoke
# failure; do NOT use for production. Use densified38 instead.
DENSIFIED34_SCHEDULE = _schedule_dict(
    lambdas_1=DENSE34_LAMBDAS_1, lambdas_2=DENSE34_LAMBDAS_2,
    directions=DENSE34_DIRECTIONS, intermd=DENSE34_INTERMD, w0=DENSE34_W0,
    alpha=DENSE34_ALPHA, u0=DENSE34_U0, lambdas=DENSE34_LAMBDAS,
)
# REVISED production free-leg ladder (softened backward peak + per-state α/U0
# ramp). Pairs with the per-direction free split (dminus base = Factor-B fix).
DENSIFIED38_SCHEDULE = _schedule_dict(
    lambdas_1=DENSE38_LAMBDAS_1, lambdas_2=DENSE38_LAMBDAS_2,
    directions=DENSE38_DIRECTIONS, intermd=DENSE38_INTERMD, w0=DENSE38_W0,
    alpha=DENSE38_ALPHA, u0=DENSE38_U0, lambdas=DENSE38_LAMBDAS,
)
# REBALANCED production free-leg ladder (count-neutral linear-λ rebalance of
# densified38). Identical soft-core ladder + α/U0/W0 ramps;
# ONLY the forward linear segment is thinned (plateau) + densified (turnover) to
# close the dplus 8→9 BC=0.131 overlap gap. Endpoints invariant → rigor-neutral.
DENSIFIED38v2_SCHEDULE = _schedule_dict(
    lambdas_1=DENSE38v2_LAMBDAS_1, lambdas_2=DENSE38v2_LAMBDAS_2,
    directions=DENSE38v2_DIRECTIONS, intermd=DENSE38v2_INTERMD, w0=DENSE38v2_W0,
    alpha=DENSE38v2_ALPHA, u0=DENSE38v2_U0, lambdas=DENSE38v2_LAMBDAS,
)
# PER-DIRECTION (NON-mirror) free-leg ladder (dminus-asymmetry fix
# 2026-06-06, Q3-b/-d). dplus = clean densified38v2 forward (UNCHANGED); dminus =
# densified38v2 backward + one W0-graded micro-bridge at the soft-core-end →
# plateau handoff (closes the dminus 9→10 BC=0.24 overlap hole the reversed
# mirror created). UNEQUAL per-direction counts (dplus 19, dminus 20 = 39 total).
# Endpoints invariant both directions → rigor-neutral / ranking-safe. densified38v2
# is NOT overwritten.
DENSIFIED38v3_SCHEDULE = _schedule_dict(
    lambdas_1=DENSE38v3_LAMBDAS_1, lambdas_2=DENSE38v3_LAMBDAS_2,
    directions=DENSE38v3_DIRECTIONS, intermd=DENSE38v3_INTERMD, w0=DENSE38v3_W0,
    alpha=DENSE38v3_ALPHA, u0=DENSE38v3_U0, lambdas=DENSE38v3_LAMBDAS,
)
# PER-DIRECTION (NON-mirror) free-leg ladder — SECOND dminus bridge
# whack-a-mole verdict 2026-06-06, Q1/Q3-a). dplus = clean densified38v2/v3
# forward (UNCHANGED, byte-equal); dminus = densified38v3 dminus (which already
# has the 9→10 W0=0.07 bridge) + one MORE W0-graded micro-bridge at the 5→6
# soft-core-ladder-end → plateau handoff (closes the dminus 5→6 BC=0.18 hole the
# v3 re-pilot exposed — the SAME single turnover, re-indexed). UNEQUAL per-
# direction counts (dplus 19, dminus 21 = 40 total). Endpoints invariant both
# directions → rigor-neutral / ranking-safe. This is the LAST schedule iteration
# under the hard iteration cap of 2 bridges/leg. densified38/v2/v3 NOT
# overwritten.
DENSIFIED38v4_SCHEDULE = _schedule_dict(
    lambdas_1=DENSE38v4_LAMBDAS_1, lambdas_2=DENSE38v4_LAMBDAS_2,
    directions=DENSE38v4_DIRECTIONS, intermd=DENSE38v4_INTERMD, w0=DENSE38v4_W0,
    alpha=DENSE38v4_ALPHA, u0=DENSE38v4_U0, lambdas=DENSE38v4_LAMBDAS,
)

# Schedules deprecated for production (Factor-B broken) — kept for forensic
# reproduction but flagged so the launcher can warn if a caller selects one.
DEPRECATED_SCHEDULES = frozenset({"densified34"})

SCHEDULES: Dict[str, Dict[str, Any]] = {
    "canonical22": CANONICAL22_SCHEDULE,
    "densified34": DENSIFIED34_SCHEDULE,     # DEPRECATED (Factor-B broken)
    "densified38": DENSIFIED38_SCHEDULE,     # REVISED production free-leg ladder
    "densified38v2": DENSIFIED38v2_SCHEDULE,  # REBALANCED linear-λ (8→9 OVL fix)
    "densified38v3": DENSIFIED38v3_SCHEDULE,  # PER-DIRECTION (dminus bridge, 9→10 fix)
    "densified38v4": DENSIFIED38v4_SCHEDULE,  # PER-DIRECTION (2nd dminus bridge, 5→6 fix)
}
DEFAULT_FREE_SCHEDULE = "canonical22"


def get_schedule(name: str) -> Dict[str, Any]:
    """Return the named schedule dict (copy). Raises KeyError on unknown name.

    Used by the cntl writer + launcher so ``--free-schedule`` selects between
    the canonical 22-state ladder (default; bound leg always uses this), the
    REVISED densified 38-state free-leg ladder (production), and the DEPRECATED
    densified 34-state ladder (forensic reproduction only — Factor-B
    broken). Selecting a deprecated schedule emits a stderr warning.
    """
    if name not in SCHEDULES:
        raise KeyError(
            f"unknown schedule {name!r}; available: {sorted(SCHEDULES)}"
        )
    if name in DEPRECATED_SCHEDULES:
        sys.stderr.write(
            f"WARNING: schedule {name!r} is DEPRECATED (Factor-B broken: NaN'd "
            f"at the backward W0-peak even with full structprep, Path REVISED "
            f"LADDER FIX 2026-06-05). Use 'densified38' for production. "
            f"Retained for forensic reproduction only.\n"
        )
    sched = SCHEDULES[name]
    return {k: (list(v) if isinstance(v, list) else v)
            for k, v in sched.items()}


# Soft-core cap constants (must match upstream OMMSystemABFE defaults +
# v1 attach_atm_force_production). uwham reads these from the cntl.
UMAX_KCAL = 200.0
UBCORE_KCAL = 100.0
ACORE = 0.062500

# Default ABFE displacement (binder fully decoupled to bulk).
# 2.5 nm matches v1; upstream temoa-g1 uses 2.2 nm. Either works as long as
# the displaced shadow ligand has no PME image overlap with the receptor.
DISPLACEMENT_NM_DEFAULT = (2.5, 0.0, 0.0)


# ---------------------------------------------------------------------------
# Test-compatibility re-exports — preserve the v1 helper API so the existing
# tests/test_atm_trackb_production.py regression suite continues to load this
# module under the same import path used by v1. The launcher does NOT call
# these in the async_re path (upstream handles ATMForce wiring), but the
# helpers are physically equivalent and useful for off-line audit work.
# ---------------------------------------------------------------------------
def attach_atm_force_production(
    system,
    displacement_nm: Tuple[float, float, float],
    displaced_atoms: List[int],
    lambda1: float = 0.0,
    lambda2: float = 0.0,
    alpha_per_kcal: float = 0.10,
    u0_kcal: float = 110.0,
    w0_kcal: float = 0.0,
    direction: float = 1.0,
    umax_kcal: float = UMAX_KCAL,
    ubcore_kcal: float = UBCORE_KCAL,
    acore: float = ACORE,
) -> int:
    """ATMForce attachment (v1-compatible helper, retained for regression).

    Identical to ``scripts/trackb_production.py::attach_atm_force_production``
    so the existing test_soft_core_constants_match_attach_defaults test
    continues to pass against this module.
    """
    import copy
    import openmm as mm
    import openmm.unit as unit

    kcal = unit.kilocalorie_per_mole
    kj = unit.kilojoule_per_mole

    kcal_per_kj = (1.0 * kcal).value_in_unit(kj)  # = 4.184
    alpha = (alpha_per_kcal / kcal_per_kj) / kj
    uh = (u0_kcal * kcal).value_in_unit(kj) * kj
    w0 = (w0_kcal * kcal).value_in_unit(kj) * kj
    umax = (umax_kcal * kcal).value_in_unit(kj) * kj
    ubcore = (ubcore_kcal * kcal).value_in_unit(kj) * kj

    atm = mm.ATMForce(lambda1, lambda2, alpha, uh, w0, umax, ubcore, acore,
                      float(direction))

    move_types = (mm.NonbondedForce, mm.HarmonicBondForce,
                  mm.HarmonicAngleForce, mm.PeriodicTorsionForce)
    to_move = [i for i in range(system.getNumForces())
               if isinstance(system.getForce(i), move_types)]
    for i in to_move:
        atm.addForce(copy.copy(system.getForce(i)))
    for i in sorted(to_move, reverse=True):
        system.removeForce(i)

    disp = mm.Vec3(*displacement_nm) * unit.nanometer
    zero = mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer
    dset = set(displaced_atoms or [])
    for idx_p in range(system.getNumParticles()):
        atm.addParticle(disp if idx_p in dset else zero)

    return system.addForce(atm)


def soft_core_pert_e(u: float, umax: float, ub: float, a: float) -> float:
    """AToM soft-core perturbation-energy cap. Verbatim from v1.

    Bounds ``u`` to ``[ub, umax]`` smoothly when ``u > ub``. uwham's bias
    function expects the capped value (this is what ``ommworker`` writes).
    """
    if u <= ub:
        return u
    gu = (u - ub) / (a * (umax - ub))
    zeta = 1.0 + 2.0 * gu * (gu + 1.0)
    zetap = zeta ** a
    return (umax - ub) * (zetap - 1.0) / (zetap + 1.0) + ub


def softplus_bias_kj(lambda1: float, lambda2: float, alpha_per_kj: float,
                     uh_kj: float, w0_kj: float, pert_kj: float) -> float:
    """AToM softplus bias energy (kJ/mol). Verbatim from v1."""
    import math
    softplus = lambda2 * pert_kj + w0_kj
    if alpha_per_kj > 0:
        ee = 1.0 + math.exp(-alpha_per_kj * (pert_kj - uh_kj))
        softplus += ((lambda2 - lambda1) / alpha_per_kj) * math.log(ee)
    return softplus


# ---------------------------------------------------------------------------
# Pre-register outcome labels (used by test_pre_register_outcome_schema).
# Verbatim from v1.
# ---------------------------------------------------------------------------
PRE_REGISTER_OUTCOMES = (
    "1_sign_stable",     # |μ| sign agrees with Magotti -1.4 SSOT, σ_btwn < 0.5
    "2_sigma_drift",     # sign agrees, σ_btwn 0.5..1.5 (mid-confidence)
    "3_magnitude_drift", # sign agrees, |μ| differs >2× from SSOT (calibration cue)
    "4_sign_flip",       # μ sign disagrees with SSOT (gate failure, escalate)
)


# ---------------------------------------------------------------------------
# Async-RE cntl file generator (canonical 22-state ABFE schedule)
# ---------------------------------------------------------------------------
def _csv(vals: List[Any]) -> str:
    """Render a list as ConfigObj-compatible comma-separated string."""
    return ", ".join(str(v) for v in vals)


def write_cntl_file(
    cntl_path: str,
    basename: str,
    nodefile_path: str,
    ligand_atom_indices: List[int],
    pos_restrained_atom_indices: List[int],
    displacement_nm: Tuple[float, float, float],
    production_steps: int,
    prnt_frequency: int,
    trj_frequency: int,
    wall_time_min: int,
    cycle_time_s: int,
    checkpoint_time_s: int,
    max_samples: Optional[int] = None,
    posre_force_constant: float = 25.0,
    posre_tolerance: float = 0.5,
    cm_kf: float = 25.0,
    cm_tol: float = 5.0,
    friction_coeff: float = 0.5,
    time_step_ps: float = 0.002,
    temperature_K: float = TEMP_K,
    openmm_platform: str = "CUDA",
    verbose: str = "no",
    subjobs_buffer_size: float = 1.0,
    thermalization_steps: Optional[int] = None,
    annealing_steps: Optional[int] = None,
    equilibration_steps: Optional[int] = None,
    steps_per_cycle: Optional[int] = None,
    schedule: Optional[Dict[str, Any]] = None,
) -> None:
    """Write the async_re control file for ``abfe_production``.

    Displacement converted from nm to Å (upstream cntl uses Å). LIGAND_ATOMS
    and LIGAND_CM_ATOMS are set to the binder chain atoms (ABFE: whole binder
    decoupled).

    ``schedule`` (optional, Path λ-densify spec 2026-06-05) selects the
    per-state λ ladder. When ``None`` (default) the canonical module-level
    22-state schedule is emitted (LAMBDAS == LAMBDA1 == LAMBDA2, INTERMEDIATE
    at the two λ=0.5 midpoints) — backward-compatible with every existing
    caller + regression test. When a schedule dict is passed (e.g.
    ``get_schedule("densified34")``), its ``lambdas`` / ``lambdas_1`` /
    ``lambdas_2`` / ``directions`` / ``intermd`` / ``alpha`` / ``u0`` / ``w0``
    arrays are emitted verbatim. The densified34 schedule is FREE-LEG ONLY;
    the bound leg always uses the canonical 22-state schedule.
    """
    displ_A = tuple(d * 10.0 for d in displacement_nm)

    # Resolve the per-state schedule. None → canonical module globals (the
    # LAMBDAS column == LAMBDA1 for the symmetric 22-state default).
    if schedule is None:
        sched_lambdas = LAMBDAS_1
        sched_lambda1 = LAMBDAS_1
        sched_lambda2 = LAMBDAS_2
        sched_directions = DIRECTIONS
        sched_intermd = INTERMD
        sched_alpha = ALPHA
        sched_u0 = U0
        sched_w0 = W0
    else:
        sched_lambdas = schedule["lambdas"]
        sched_lambda1 = schedule["lambdas_1"]
        sched_lambda2 = schedule["lambdas_2"]
        sched_directions = schedule["directions"]
        sched_intermd = schedule["intermd"]
        sched_alpha = schedule["alpha"]
        sched_u0 = schedule["u0"]
        sched_w0 = schedule["w0"]

    lines = [
        "# Track B v2 ABFE control file — generated by trackb_production_v2_asyncre.py",
        "# DO NOT EDIT BY HAND — re-run the launcher to regenerate.",
        "",
        "JOB_TRANSPORT = 'LOCAL_OPENMM'",
        f"BASENAME = '{basename}'",
        f"NODEFILE = '{os.path.basename(nodefile_path)}'",
        "",
        f"TEMPERATURES = '{temperature_K}'",
        f"LAMBDAS =      '{_csv(sched_lambdas)}'",
        f"DIRECTION =    '{_csv(sched_directions)}'",
        f"INTERMEDIATE = '{_csv(sched_intermd)}'",
        f"LAMBDA1 =      '{_csv(sched_lambda1)}'",
        f"LAMBDA2 =      '{_csv(sched_lambda2)}'",
        f"ALPHA =        '{_csv(sched_alpha)}'",
        f"U0 =           '{_csv(sched_u0)}'",
        f"W0COEFF =      '{_csv(sched_w0)}'",
        "",
        f"DISPLACEMENT = '{displ_A[0]}, {displ_A[1]}, {displ_A[2]}'",
        "LIGOFFSET = '0., 0., 0.'",
        "",
        f"WALL_TIME = {wall_time_min}",
        f"CYCLE_TIME = {cycle_time_s}",
        f"CHECKPOINT_TIME = {checkpoint_time_s}",
        f"SUBJOBS_BUFFER_SIZE = '{subjobs_buffer_size}'",
        "",
        f"PRODUCTION_STEPS = '{production_steps}'",
        f"PRNT_FREQUENCY = '{prnt_frequency}'",
        f"TRJ_FREQUENCY = '{trj_frequency}'",
    ]
    if max_samples is not None:
        lines.append(f"MAX_SAMPLES = {max_samples}")

    lines += [
        "",
        f"LIGAND_ATOMS = {_csv(ligand_atom_indices)}",
        f"LIGAND_CM_ATOMS = {_csv(ligand_atom_indices)}",
        # No RCPT_CM_ATOMS or CM_KF/CM_TOL on the free leg: the free-leg PDB
        # has no receptor. We omit these keys entirely when the caller passes
        # an empty pos_restrained list (the free leg sentinel).
    ]
    # POSRE constants are ALWAYS written because abfe_structprep's
    # massage_keywords() auto-populates POS_RESTRAINED_ATOMS from non-solvent
    # atoms during mintherm (MINTHERM_RESTRAIN_SOLUTES default=YES). Without
    # POSRE_FORCE_CONSTANT/POSRE_TOLERANCE, OMMSystem.set_positional_restraints
    # raises TypeError on float(None). The constants are harmless on the free
    # leg: the auto-populated POS_RESTRAINED_ATOMS = all peptide atoms (a
    # reasonable mintherm restraint), and POSRE is restored to the cntl's
    # original value (typically None for free leg) after mintherm completes.
    lines += [
        f"POSRE_FORCE_CONSTANT = {posre_force_constant}",
        f"POSRE_TOLERANCE = {posre_tolerance}",
    ]

    if pos_restrained_atom_indices:
        # Bound leg: CM-CM Vsite restraint anchors displaced binder vs receptor
        # CA atoms; position restraints stabilize the receptor.
        rcpt_cm = _csv(pos_restrained_atom_indices[:80])  # subset for CM
        lines += [
            f"RCPT_CM_ATOMS = {rcpt_cm}",
            f"CM_KF = {cm_kf}",
            f"CM_TOL = {cm_tol}",
            f"POS_RESTRAINED_ATOMS = {_csv(pos_restrained_atom_indices)}",
        ]

    # Optional structprep step counts (THERMALIZATION/ANNEALING/EQUILIBRATION).
    # Upstream defaults are 150000/250000/150000 = ~550 ps total at 2fs.
    # Smoke mode reduces to ~10 ps to keep wall-time under 1 min.
    if thermalization_steps is not None:
        lines.append(f"THERMALIZATION_STEPS = {thermalization_steps}")
    if annealing_steps is not None:
        lines.append(f"ANNEALING_STEPS = {annealing_steps}")
    if equilibration_steps is not None:
        lines.append(f"EQUILIBRATION_STEPS = {equilibration_steps}")
    if steps_per_cycle is not None:
        lines.append(f"STEPS_PER_CYCLE = {steps_per_cycle}")

    lines += [
        "",
        f"UMAX = {UMAX_KCAL}",
        f"ACORE = {ACORE}",
        f"UBCORE = {UBCORE_KCAL}",
        "",
        f"FRICTION_COEFF = {friction_coeff}",
        f"TIME_STEP = {time_step_ps}",
        f"OPENMM_PLATFORM = {openmm_platform}",
        f"VERBOSE = '{verbose}'",
        "",
    ]

    with open(cntl_path, "w") as fh:
        fh.write("\n".join(lines))


def write_nodefile(nodefile_path: str,
                   gpu_indices: List[int],
                   platform: str = "CUDA",
                   tmp_root: str = "/tmp") -> None:
    """Write the async_re nodefile (one row per GPU device).

    Row format: ``<node_name>,<platform_id>:<device_id>,<threads>,<arch>,<user>,<tmp_folder>``
    The middle two ints are platform_id (for CUDA: ignored) and device_id
    (the CUDA index; with CUDA_VISIBLE_DEVICES pinned, always 0).
    """
    if not gpu_indices:
        gpu_indices = [0]
    with open(nodefile_path, "w") as fh:
        for d in gpu_indices:
            fh.write(f"localhost,0:{d},8,{platform},,{tmp_root}\n")


# ---------------------------------------------------------------------------
# System XML / topology PDB build — per leg (delegates to atm_trackB_setup)
# ---------------------------------------------------------------------------
def build_leg_xml(
    pdb_path: str,
    basename: str,
    leg_dir: str,
    add_hydrogens: bool,
    hydrogens_xml: Optional[str],
    binder_chain: str = "B",
) -> Dict[str, Any]:
    """Build ``<basename>.pdb`` + ``<basename>_sys.xml`` for one leg.

    Reuses ``atm_trackB_setup.build_leg_system`` which is the validated
    UPDD path (amber14-all + tip3pfb + MTR_gaff2_hybrid.xml, MTR hydrogen
    definitions registered, ncAA peptide bonds added, cyclic_ss disulfide
    committed). The output XML is fed directly to ``abfe_structprep`` —
    upstream sees a complete System with all bonded/non-bonded forces ready
    for the ATMForce wrapping that ``OMMSystemABFE.create_system`` performs.

    Returns the build dict (atoms, alch indices, disulfide) for run_metadata.
    """
    import openmm
    from openmm.app import PDBFile
    from atm_trackB_setup import build_leg_system

    build = build_leg_system(
        pdb_path,
        leg=os.path.basename(leg_dir),
        binder_chain=binder_chain,
        solvate=True,
        add_hydrogens=add_hydrogens,
        hydrogens_xml=hydrogens_xml,
    )
    system = build["system"]
    modeller = build["modeller"]

    # Set the periodic box from the modeller before serialization (otherwise
    # PME nonbonded barfs at structprep time on a non-periodic-aware system).
    try:
        boxvec = modeller.topology.getPeriodicBoxVectors()
        if boxvec is not None:
            system.setDefaultPeriodicBoxVectors(*boxvec)
    except Exception:
        # Modeller.addSolvent already set the box; if not, leave system alone.
        pass

    os.makedirs(leg_dir, exist_ok=True)
    pdb_out = os.path.join(leg_dir, basename + ".pdb")
    xml_out = os.path.join(leg_dir, basename + "_sys.xml")

    with open(pdb_out, "w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
    with open(xml_out, "w") as fh:
        fh.write(openmm.XmlSerializer.serialize(system))

    return {
        "pdb_path": pdb_out,
        "xml_path": xml_out,
        "n_atoms": build["n_atoms"],
        "alchemical_atoms": build["alchemical_atoms"],
        "disulfide": build["disulfide"],
        "n_peptide_bonds_added": build["n_peptide_bonds_added"],
        "n_internal_bonds_added": build["n_internal_bonds_added"],
        "ff_inputs": build["ff_inputs"],
    }


def collect_binder_atom_indices(pdb_path: str, binder_chain: str = "B") -> List[int]:
    """Collect every atom index of the binder chain from a topology PDB.

    Used to populate LIGAND_ATOMS in the cntl. The indices are global topology
    indices (matching the order of the PDB used by abfe_structprep).
    """
    from openmm.app import PDBFile

    pdb = PDBFile(pdb_path)
    idxs: List[int] = []
    for chain in pdb.topology.chains():
        if chain.id != binder_chain:
            continue
        for atom in chain.atoms():
            idxs.append(atom.index)
    return idxs


def collect_receptor_ca_indices(pdb_path: str,
                                receptor_chain: str = "A") -> List[int]:
    """Collect receptor Cα atom indices (for RCPT_CM_ATOMS + POS_RESTRAINED).

    Returns the indices of the protein backbone Cα atoms of the receptor
    chain; uwham / position-restraint helpers use this set both for the
    CM-CM Vsite restraint anchor and (a subset of) the position restraints.
    """
    from openmm.app import PDBFile

    pdb = PDBFile(pdb_path)
    ca: List[int] = []
    for chain in pdb.topology.chains():
        if chain.id != receptor_chain:
            continue
        for residue in chain.residues():
            for atom in residue.atoms():
                if atom.name == "CA":
                    ca.append(atom.index)
    return ca


# ---------------------------------------------------------------------------
# Resolve inputs (per-endpoint per-leg PDB).
# ---------------------------------------------------------------------------
def resolve_endpoint_pdbs(seed_tag: str,
                          out_root: str) -> Dict[Tuple[str, str], str]:
    """Resolve the prepared (endpoint, leg) → PDB mapping.

    Bound leg uses the per-endpoint ``_md_input`` complex PDB. Free leg
    uses the H-complete ``mdresult/<endpoint>_final.pdb`` extracted to
    binder-only via ``prepare_free_peptide_from_final``. The extraction
    is cached under ``<out_root>/_free_pdbs/<endpoint>_freeleg.pdb`` so
    repeated launches do not re-extract.
    """
    from atm_trackB_setup import resolve_leg_inputs, prepare_free_peptide_from_final

    inputs = resolve_leg_inputs(seed_tag)
    cp4_bound = inputs["bound"]["cp4"]
    wt_bound = inputs["bound"]["wt"]
    cp4_final = inputs["final"]["cp4"]
    wt_final = inputs["final"]["wt"]

    free_work = os.path.join(out_root, "_free_pdbs")
    os.makedirs(free_work, exist_ok=True)
    cp4_free = os.path.join(free_work, "cp4_freeleg.pdb")
    wt_free = os.path.join(free_work, "wt_freeleg.pdb")
    if not os.path.isfile(cp4_free):
        if cp4_final and os.path.isfile(cp4_final):
            prepare_free_peptide_from_final(cp4_final, cp4_free)
        else:
            raise FileNotFoundError(
                f"Cp4 final.pdb missing for free-leg extraction: {cp4_final}"
            )
    if not os.path.isfile(wt_free):
        if wt_final and os.path.isfile(wt_final):
            prepare_free_peptide_from_final(wt_final, wt_free)
        else:
            raise FileNotFoundError(
                f"WT final.pdb missing for free-leg extraction: {wt_final}"
            )

    return {
        ("cp4", "bound"): cp4_bound,
        ("cp4", "free"): cp4_free,
        ("wt", "bound"): wt_bound,
        ("wt", "free"): wt_free,
    }


# ---------------------------------------------------------------------------
# Top-level per-leg orchestration.
# ---------------------------------------------------------------------------
def setup_one_leg(
    endpoint: str,
    leg: str,
    pdb_path: str,
    out_root: str,
    jobname: str = "trackb",
    binder_chain: str = "B",
    receptor_chain: str = "A",
    hydrogens_xml: Optional[str] = None,
    displacement_nm: Tuple[float, float, float] = DISPLACEMENT_NM_DEFAULT,
    cuda_devices: Optional[List[int]] = None,
    production_steps: int = 2500,    # 5 ps per cycle @ 2fs (production-scale)
    prnt_frequency: int = 2500,      # match production_steps
    trj_frequency: int = 25000,      # save trajectory every 10 cycles
    max_samples: int = 1000,         # 1000 cycles ≈ 5 ns per replica
    wall_time_min: int = 720,        # 12 h per leg
    cycle_time_s: int = 10,
    checkpoint_time_s: int = 600,
    smoke: bool = False,
    free_schedule: str = DEFAULT_FREE_SCHEDULE,
) -> Dict[str, Any]:
    """Prepare one (endpoint, leg) directory for ``abfe_production``.

    Builds the system XML + topology PDB, writes the cntl + nodefile,
    runs ``abfe_structprep`` if not already done. Production launch is
    deliberately NOT triggered here — the caller drives it via
    ``run_abfe_production_for_leg``.

    In smoke mode the production cycle counts are clamped to a small
    budget (50 steps per cycle, max_samples=2) so the full pipeline runs
    end-to-end in minutes.

    ``free_schedule`` (default ``"canonical22"``) selects the λ ladder for
    the FREE leg only (Path λ-densify spec 2026-06-05). ``"densified34"``
    swaps the free-leg cntl to the 34-state ilogistic anneal ladder. The
    BOUND leg ALWAYS uses the canonical 22-state schedule regardless of this
    value (the receptor holds the ligand → no decoupling crossover → no gap).
    """
    leg_dir = os.path.join(out_root, endpoint, leg)
    os.makedirs(leg_dir, exist_ok=True)

    # Hydrogen-build policy (verbatim from v1):
    # bound leg uses pre-MD _md_input PDB (H-incomplete for ncAA, needs
    #   addHydrogens + MTR_hydrogens.xml).
    # free leg uses MD final.pdb (H-complete, skip addHydrogens — re-running
    #   it on a fully protonated ncAA mis-handles the methyl set).
    add_h = (leg == "bound")

    # Smoke override: reduce production cycles AND structprep step counts.
    # Default upstream THERMALIZATION/ANNEALING/EQUILIBRATION = 550k steps
    # (~1100 ps @ 2fs ≈ 4 min on 5070Ti) — fine for the first smoke but
    # painful for iteration. Smoke clamps to ~5 ps total structprep + 2
    # production cycles = end-to-end in ~30 s.
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
        # Aggressive structprep clamp — proves the pipeline, not the
        # thermodynamics. Production runs use upstream defaults.
        thermalization_steps = 1000
        annealing_steps = 2000
        equilibration_steps = 1000
        steps_per_cycle = 100

    # Build the system XML + topology PDB.
    build_info = build_leg_xml(
        pdb_path=pdb_path,
        basename=jobname,
        leg_dir=leg_dir,
        add_hydrogens=add_h,
        hydrogens_xml=hydrogens_xml,
        binder_chain=binder_chain,
    )
    topology_pdb = build_info["pdb_path"]

    # Collect atom indices for the cntl.
    ligand_atoms = collect_binder_atom_indices(topology_pdb,
                                               binder_chain=binder_chain)
    if not ligand_atoms:
        raise RuntimeError(
            f"No binder-chain atoms found in {topology_pdb} (chain {binder_chain})"
        )

    if leg == "bound":
        pos_restrained = collect_receptor_ca_indices(
            topology_pdb, receptor_chain=receptor_chain
        )
        if not pos_restrained:
            raise RuntimeError(
                f"No receptor Cα atoms found in {topology_pdb} "
                f"(chain {receptor_chain}) — bound leg needs RCPT_CM_ATOMS"
            )
    else:
        # Free leg: no receptor → no CM restraint, no position restraint.
        # abfe_structprep auto-restrains all non-solvent atoms during mintherm
        # via the MINTHERM_RESTRAIN_SOLUTES path; we let that default take over.
        pos_restrained = []

    # Write the nodefile.
    nodefile_path = os.path.join(leg_dir, "nodefile")
    write_nodefile(nodefile_path,
                   gpu_indices=(cuda_devices or [0]),
                   platform="CUDA")

    # Resolve the λ schedule. densified34 applies to the FREE leg only; the
    # bound leg always uses the canonical 22-state schedule (Path: receptor
    # holds the ligand → no crossover gap → no densification needed).
    schedule_dict: Optional[Dict[str, Any]] = None
    schedule_name = DEFAULT_FREE_SCHEDULE
    if leg == "free" and free_schedule != DEFAULT_FREE_SCHEDULE:
        schedule_dict = get_schedule(free_schedule)
        schedule_name = free_schedule

    # Write the cntl.
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
        "system_xml": build_info["xml_path"],
        "n_atoms": build_info["n_atoms"],
        "n_ligand_atoms": len(ligand_atoms),
        "n_pos_restrained": len(pos_restrained),
        "disulfide": build_info["disulfide"],
        "alchemical_atoms": build_info["alchemical_atoms"],
        "n_peptide_bonds_added": build_info["n_peptide_bonds_added"],
        "n_internal_bonds_added": build_info["n_internal_bonds_added"],
        "smoke_mode": smoke,
        "schedule_name": schedule_name,
        "n_states": (schedule_dict["n_states"] if schedule_dict is not None
                     else N_STATES),
    }


def _find_atm_binary(name: str) -> str:
    """Locate an atom_openmm binary, falling back to the atm conda env.

    Honors PATH first, then probes the canonical ``atm`` env install location.
    Raises if neither finds the binary.
    """
    found = shutil.which(name)
    if found:
        return found
    fallback = f"/home/san/miniconda3/envs/atm/bin/{name}"
    if os.path.isfile(fallback) and os.access(fallback, os.X_OK):
        return fallback
    raise RuntimeError(
        f"{name} not found on PATH and not at {fallback}. "
        f"Activate the atm conda env or install atom_openmm."
    )


def run_abfe_structprep_for_leg(leg_info: Dict[str, Any],
                                log_path: Optional[str] = None) -> int:
    """Invoke ``abfe_structprep`` on one leg's cntl. Returns exit code."""
    cntl = leg_info["cntl_path"]
    binary = _find_atm_binary("abfe_structprep")
    cmd = [binary, os.path.basename(cntl)]
    cwd = os.path.dirname(cntl)
    log_path = log_path or os.path.join(cwd, "_structprep.log")
    with open(log_path, "w") as logfh:
        logfh.write(f"# Command: {' '.join(shlex.quote(c) for c in cmd)}\n")
        logfh.write(f"# CWD:     {cwd}\n")
        logfh.write(f"# Started: {time.strftime('%Y-%m-%dT%H:%M:%S')}\n\n")
        logfh.flush()
        proc = subprocess.run(cmd, cwd=cwd, stdout=logfh, stderr=subprocess.STDOUT)
    return proc.returncode


def run_abfe_production_for_leg(leg_info: Dict[str, Any],
                                log_path: Optional[str] = None,
                                background: bool = False) -> Any:
    """Invoke ``abfe_production`` on one leg's cntl.

    When background=True, returns the spawned subprocess.Popen object (caller
    is responsible for monitoring). When False, returns the exit code after
    abfe_production completes.
    """
    cntl = leg_info["cntl_path"]
    binary = _find_atm_binary("abfe_production")
    cmd = [binary, os.path.basename(cntl)]
    cwd = os.path.dirname(cntl)
    log_path = log_path or os.path.join(cwd, "_production.log")
    logfh = open(log_path, "w")
    logfh.write(f"# Command: {' '.join(shlex.quote(c) for c in cmd)}\n")
    logfh.write(f"# CWD:     {cwd}\n")
    logfh.write(f"# Started: {time.strftime('%Y-%m-%dT%H:%M:%S')}\n\n")
    logfh.flush()
    if background:
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=logfh, stderr=subprocess.STDOUT)
        return proc
    try:
        proc = subprocess.run(cmd, cwd=cwd, stdout=logfh, stderr=subprocess.STDOUT)
        return proc.returncode
    finally:
        logfh.close()


# ---------------------------------------------------------------------------
# Top-level main()
# ---------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Track B production v2 — upstream atom_openmm async_re ABFE "
            "(per-endpoint MTR↔Trp ΔΔG_bind)."
        )
    )
    p.add_argument("--seed-tag", default="s7",
                   help="prepared seed under outputs/2QKI_{Cp4_hybrid,WT}_calib_<seed>")
    p.add_argument("--out-root", default="outputs/_trackb/production_v2",
                   help="output root under project")
    p.add_argument("--endpoints", default="cp4,wt",
                   help="comma-list of endpoints to run (cp4 / wt)")
    p.add_argument("--legs", default="bound,free",
                   help="comma-list of legs to run (bound / free)")
    p.add_argument("--displacement-nm", default="2.5,0.0,0.0",
                   help="ABFE displacement vector (nm)")
    p.add_argument("--cuda-device", default="0",
                   help="forced CUDA_VISIBLE_DEVICES (host 5070Ti = 0)")
    p.add_argument("--platform", default="CUDA")
    p.add_argument("--jobname", default="trackb",
                   help="basename of the per-leg cntl / pdb / xml files")
    p.add_argument("--production-steps", type=int, default=2500,
                   help="MD steps per replica per cycle (default 2500 = 5ps@2fs)")
    p.add_argument("--max-samples", type=int, default=1000,
                   help="max cycles per replica (default 1000 ≈ 5 ns per replica)")
    p.add_argument("--wall-time-min", type=int, default=720,
                   help="abfe_production wall-time budget per leg in minutes")
    p.add_argument("--checkpoint-time-s", type=int, default=600,
                   help="checkpoint interval (s)")
    p.add_argument("--cycle-time-s", type=int, default=10,
                   help="async-RE swap-attempt interval (s)")
    p.add_argument("--smoke", action="store_true",
                   help="smoke gate: 50 steps × 2 cycles × 1 leg "
                        "(host 5070Ti, ~10-15 min)")
    p.add_argument("--free-schedule", default=DEFAULT_FREE_SCHEDULE,
                   choices=sorted(SCHEDULES),
                   help=("λ ladder for the FREE leg only (Path λ-densify "
                         "spec 2026-06-05). 'canonical22' (default) = 22-state "
                         "symmetric; 'densified38' = REVISED 38-state ilogistic "
                         "anneal ladder (softened backward W0-peak + per-state "
                         "α/U0 ramp) bridging the λ=0.45→0.50 decoupling "
                         "crossover (PRODUCTION free-leg ladder); 'densified38v2' "
                         "= REBALANCED densified38 (count-neutral linear-λ "
                         "thin-plateau + densify-turnover, fixes the "
                         "dplus 8→9 BC=0.131 overlap gap); 'densified38v3' = "
                         "PER-DIRECTION (NON-mirror): dplus = clean densified38v2 "
                         "forward (19); dminus = densified38v2 backward + ONE "
                         "W0-graded micro-bridge at the soft-core-end→plateau "
                         "handoff (20), closing the dminus 9→10 BC=0.24 hole "
                         "(dminus-asymmetry fix 2026-06-06; UNEQUAL "
                         "39-state); 'densified38v4' = densified38v3 + a SECOND "
                         "dminus W0-graded micro-bridge at the re-indexed 5→6 "
                         "handoff (dminus 21, UNEQUAL 40-state), closing the "
                         "dminus 5→6 BC=0.18 hole the v3 re-pilot exposed "
                         "(2026-06-06; LAST iteration under "
                         "the hard 2-bridge/leg cap); 'densified34' = "
                         "DEPRECATED (Factor-B broken, forensic only). "
                         "The bound leg always uses canonical22."))
    p.add_argument("--setup-only", action="store_true",
                   help="prepare cntl/system/structprep but DO NOT launch production")
    p.add_argument("--skip-structprep", action="store_true",
                   help="assume <basename>_0.xml already exists from a prior structprep")
    p.add_argument("--background", action="store_true",
                   help="spawn abfe_production in background (Popen) and return")
    p.add_argument("--no-charge-axis-gate", action="store_true",
                   help="skip the charge_axis verifier pre-flight (debug only)")
    args = p.parse_args()

    # CUDA pin (hardware contract — Track A V100 untouched)
    os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda_device
    cuda_devs = [int(d) for d in args.cuda_device.split(",") if d.strip()]

    # Charge-axis gate (PATCH-01 surface)
    if not args.no_charge_axis_gate:
        from atm_trackB_setup import verify_charge_axis
        axis = verify_charge_axis()
        if not axis["all_residues_within_tol"]:
            print(f"ERROR: charge axis FAIL "
                  f"(max|Σq|={axis['max_abs_sigma_q_e']:.2e})",
                  file=sys.stderr)
            return 2

    # Parse displacement vector
    displ_parts = [float(x) for x in args.displacement_nm.split(",")]
    if len(displ_parts) != 3:
        print(f"ERROR: --displacement-nm needs 3 floats, got {displ_parts}",
              file=sys.stderr)
        return 2
    displacement_nm = tuple(displ_parts)

    # Resolve inputs
    out_root = os.path.join(_PROJ_ROOT, args.out_root)
    os.makedirs(out_root, exist_ok=True)
    endpoint_pdb = resolve_endpoint_pdbs(args.seed_tag, out_root)

    # Resolve hydrogens xml (used by bound legs only)
    from atm_trackB_setup import resolve_leg_inputs
    inputs = resolve_leg_inputs(args.seed_tag)
    hydrogens_xml = inputs.get("hydrogens_xml")

    endpoints = [e.strip() for e in args.endpoints.split(",") if e.strip()]
    legs = [l.strip() for l in args.legs.split(",") if l.strip()]
    for e in endpoints:
        assert e in ("cp4", "wt"), f"unknown endpoint {e}"
    for l in legs:
        assert l in ("bound", "free"), f"unknown leg {l}"

    # In smoke mode, restrict to one leg unless caller overrides explicitly.
    if args.smoke and "--endpoints" not in sys.argv and "--legs" not in sys.argv:
        endpoints = ["cp4"]
        legs = ["free"]
        print("# smoke mode: restricted to cp4/free (override with --endpoints/--legs)")

    # Compose verbatim launch command for audit.
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

    # ---- Per-leg orchestration ----------------------------------------
    setup_results: List[Dict[str, Any]] = []
    structprep_results: List[Dict[str, Any]] = []
    production_results: List[Dict[str, Any]] = []
    procs = []

    for endpoint in endpoints:
        for leg in legs:
            pdb_path = endpoint_pdb[(endpoint, leg)]
            print(f"\n=== SETUP   {endpoint} / {leg} ===")
            print(f"    pdb: {pdb_path}")
            leg_info = setup_one_leg(
                endpoint=endpoint,
                leg=leg,
                pdb_path=pdb_path,
                out_root=out_root,
                jobname=args.jobname,
                hydrogens_xml=hydrogens_xml,
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
            )
            setup_results.append(leg_info)
            print(f"    cntl: {leg_info['cntl_path']}")
            print(f"    n_atoms: {leg_info['n_atoms']}  "
                  f"n_ligand: {leg_info['n_ligand_atoms']}  "
                  f"n_restr: {leg_info['n_pos_restrained']}")
            print(f"    disulfide: {leg_info['disulfide']}")

            if not args.skip_structprep:
                # Idempotent: skip if <basename>_0.xml already present.
                ready_marker = os.path.join(
                    leg_info["leg_dir"], args.jobname + "_0.xml"
                )
                if os.path.isfile(ready_marker):
                    print(f"    structprep already done (skip): {ready_marker}")
                    structprep_results.append({"endpoint": endpoint, "leg": leg,
                                               "status": "cached"})
                else:
                    print(f"--- STRUCTPREP {endpoint} / {leg} ---")
                    rc = run_abfe_structprep_for_leg(leg_info)
                    structprep_results.append({"endpoint": endpoint, "leg": leg,
                                               "rc": rc,
                                               "log": os.path.join(
                                                   leg_info["leg_dir"],
                                                   "_structprep.log")})
                    if rc != 0:
                        print(f"    structprep FAILED rc={rc} — see _structprep.log",
                              file=sys.stderr)
                        return 3

    # ---- run_metadata --------------------------------------------------
    run_metadata = {
        "track": "B",
        "method": "ATM/ATS ABFE per-endpoint (async_re, Gallicchio 2021)",
        "launcher": "scripts/trackb_production_v2_asyncre.py",
        "method_ref": "per-replica async_re ABFE (architectural NaN fix), 2026-05-31",
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
        # Per-leg free-leg schedule selector (Path λ-densify spec 2026-06-05).
        # The bound leg ALWAYS uses canonical22; free_schedule applies to the
        # free leg only. When densified34, the full 34-state arrays are
        # recorded for audit of the free=34/bound=22 asymmetry.
        "free_schedule": args.free_schedule,
        "free_schedule_detail": (
            None if args.free_schedule == DEFAULT_FREE_SCHEDULE
            else get_schedule(args.free_schedule)
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
        "setup_results": setup_results,
        "structprep_results": structprep_results,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "regime": "ranking_only",
    }
    meta_path = os.path.join(out_root, "run_metadata.json")
    if os.path.isfile(meta_path):
        # Sidecar: preserve the original run_metadata; new launches get a
        # timestamped sidecar (mirrors v1 resume sidecar convention).
        ts_tag = time.strftime("%Y%m%dT%H%M%S")
        meta_path = os.path.join(out_root,
                                 f"run_metadata_{ts_tag}.json")
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

    # ---- Launch production --------------------------------------------
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
        # Persist the launched PIDs / log paths for the monitor.
        pid_path = os.path.join(out_root, "background_pids.json")
        with open(pid_path, "w") as fh:
            json.dump({"procs": procs, "launched_at": time.strftime(
                "%Y-%m-%dT%H:%M:%S")}, fh, indent=2)
        print(f"\nBackground PIDs: {pid_path}")

    return 0


def _atom_openmm_version() -> str:
    try:
        import atom_openmm
        v = getattr(atom_openmm, "__version__", None)
        if v:
            return v
        from atom_openmm.async_re import __version__ as v2
        return v2
    except Exception:
        return "unknown"


if __name__ == "__main__":
    sys.exit(main())
