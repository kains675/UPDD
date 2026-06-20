#!/usr/bin/env python
"""
sanity_gate_pbsa.py
-------------------
Numerical sanity gate for MM-PBSA ΔG_bind results.

Scope and intent
================
MM-PBSA (and MM-GBSA) are used in this pipeline as a **cheap-triage /
ranking** layer only. They are NOT authoritative for the *sign* of a
binding free energy. The gas-phase term (ΔG_gas = ΔE_vdW + ΔE_EEL) and the
solvation term (ΔG_solv = ΔE_PB + ΔE_nonpolar) carry opposite signs and
nearly cancel, so any single net total ΔG_bind is the small residual of two
large opposing numbers. That is the universal state of a normal binding
endpoint — at the per-snapshot / per-system level the cancellation is
*always* near-total, so an endpoint-level test cannot discriminate a
trustworthy result from an untrustworthy one.

The diagnostic that DOES discriminate lives one level up, at the
**variant-pair difference** ΔΔ (e.g. Cp4 cohort − WT cohort). There the
quantity of interest, ΔΔnet, is itself the small residual of two large
opposing differences ΔΔgas and ΔΔsolv. When that residual is a small
fraction of the magnitude of its parts, a few-percent systematic error in
either part can flip the sign of ΔΔnet, so the ΔΔ sign is numerically
ill-conditioned and must not be trusted for a ranking decision.

This gate measures exactly that conditioning at the ΔΔ level:

    r1 = abs(ΔΔnet) / (abs(ΔΔgas) + abs(ΔΔsolv))
    sign_ill_conditioned = (r1 < SANITY_GATE_CANCEL_C)     # c = 0.15

with a continuous companion ratio r2 = abs(ΔΔnet) / max(abs(ΔΔgas),
abs(ΔΔsolv)) reported alongside. r1 is the reciprocal of the relative
condition number of the subtraction net = gas + solv: when
net = gas + solv, κ ≈ (abs(gas) + abs(solv)) / abs(net) = 1 / r1. r1 < 0.15
⟺ κ > 6.7 ⟺ a 1% error in the inputs amplifies to ~7%+ in the residual —
the standard definition of catastrophic subtractive cancellation
(Higham 2002, §1.7).

What this gate does and does NOT claim
======================================
- It DOES claim: "for this variant pair the ΔΔ sign is numerically
  ill-conditioned and unreliable; route to a more expensive method".
- It DOES NOT claim any root cause. In particular it does NOT assert that
  the force-field gas-phase term is "wrong" — the underlying cause may be
  conformational sampling (e.g. a partially disengaged binding pose),
  which no scoring method can repair. The tag is a *phenomenological*
  description of numerical conditioning, not an attribution.

z is NOT a gate input
=====================
The statistical resolution z = abs(ΔΔnet) / SE_pooled is reported for
context but MUST NOT be used to decide the gate. A same-chemistry null pair
(true ΔΔ = 0) can still show ΔΔnet ≈ +9.6 kcal/mol with z > 4 ("resolved")
because that residual is a *systematic* between-seed cancellation artifact,
not statistical noise; SE captures only within-cohort variance. A
reproducible (high-z) ΔΔnet that is the small difference of two large terms
is still sign-ill-conditioned. Only r1 exposes that conditioning.

Cheap-triage demotion semantics
===============================
A flagged result is not deleted. The flag is a router signal: pairs with an
ill-conditioned ΔΔ sign should be sent to an expensive method (QM/MM, FEP),
while the MM-PBSA energies, components, and JSON are preserved intact.

Regime
======
ranking / cheap-triage only. This gate produces no new magnitude or
absolute ΔG claim; it only tags the reliability of the existing ΔΔ sign.
Magnitudes remain ranking-only and must not be used for absolute
comparison.

Modes
=====
- Mode B (paired ΔΔ across a variant pair) is the ONLY gate. It fires at
  the cohort/summary level where a variant pair is defined
  (net = ΔΔG_bind, gas = ΔΔG_gas, solv = ΔΔG_solv). Use
  ``apply_pbsa_pair_sanity_gate``.
- A per-snap / per-system endpoint has NO gate. The cancellation there is
  universal and undiscriminating, so ``apply_pbsa_sanity_gate`` only logs
  the endpoint ratio r_A = abs(net) / (abs(gas) + abs(solv)) for visibility
  and never raises a flag.

Reference
=========
- Higham N J 2002, *Accuracy and Stability of Numerical Algorithms*, 2nd
  ed., SIAM, §1.7 (subtractive cancellation; condition number
  κ = Σ abs(term) / abs(sum)). ISBN 0-89871-521-0.
- Genheden S, Ryde U 2015, Expert Opin. Drug Discov. 10(5):449-461,
  DOI 10.1517/17460441.2015.1032936 (systematic ≫ statistical error;
  single-traj cancellation limits).
- Hou T et al. 2011, J. Chem. Inf. Model. 51(1):69-82,
  DOI 10.1021/ci100275a (EEL/PB cancellation governs ΔG accuracy).

Limitations
===========
- 1-trajectory only. The 3-trajectory protocol does not cancel the
  intramolecular gas term and would change the ΔΔgas magnitude; introducing
  3traj support is a separate task.
- The cutoff c = 0.15 is calibrated on the 2QKI Cp4/WT cohorts, which are
  entirely in the cancellation regime. A reference variant pair with a
  large, well-resolved ΔΔ (r1 ≫ 0.2) is absent in that system, so the
  upper-side separation of c is not empirically verified here; re-calibrate
  c on another system (e.g. 7TL8) when such data exist.
- MM-GBSA decomposition is not supported here: the OpenMM single-total GBSA
  energy is not separated into gas/solv force-group components, so this gate
  cannot be applied to MM-GBSA results (follow-up).
"""
from typing import Any, Dict, Optional, Tuple


# Catastrophic-cancellation cutoff for the ΔΔ ratio r1. r1 < c flags the
# ΔΔ sign as numerically ill-conditioned. c = 0.15 ⟺ condition number
# κ ≈ 1/r1 > 6.7. Calibrated on the 2QKI Cp4/WT cohorts (canonical Cp4 ΔΔ
# r1 ≈ 0.069; same-chemistry null pairs ≤ ~0.09; over-flag cost is only an
# expensive recompute, never data loss). Exposed as a module constant.
SANITY_GATE_CANCEL_C = 0.15


def _cancellation_ratios(
    dd_gas: float, dd_solv: float, dd_net: float
) -> Tuple[Optional[float], Optional[float]]:
    """Return (r1, r2) for a ΔΔ decomposition.

    r1 = abs(dd_net) / (abs(dd_gas) + abs(dd_solv))  — reciprocal condition
        number of net = gas + solv (small ⟺ ill-conditioned).
    r2 = abs(dd_net) / max(abs(dd_gas), abs(dd_solv)) — companion continuous
        indicator (r2 ≈ 2*r1 in the balanced-cancellation regime).

    Both are None when the denominator is zero (no cancellation magnitude to
    speak of).
    """
    denom1 = abs(dd_gas) + abs(dd_solv)
    denom2 = max(abs(dd_gas), abs(dd_solv))
    r1 = (abs(dd_net) / denom1) if denom1 != 0.0 else None
    r2 = (abs(dd_net) / denom2) if denom2 != 0.0 else None
    return r1, r2


def apply_pbsa_pair_sanity_gate(
    dd_gas: float,
    dd_solv: float,
    dd_net: float,
    c: float = SANITY_GATE_CANCEL_C,
) -> Dict[str, Any]:
    """Mode B gate — tag a variant-pair ΔΔ sign for catastrophic cancellation.

    This is the ONLY sanity gate that raises a flag. It operates on the
    ΔΔ(variant_A − variant_B) decomposition (e.g. Cp4 cohort − WT cohort):

        r1 = abs(dd_net) / (abs(dd_gas) + abs(dd_solv))
        sign_ill_conditioned_cancellation = (r1 < c)        # c = 0.15

    The returned dict has additive ΔΔ-level keys:

    - ``sign_ill_conditioned_cancellation`` (bool): True when r1 < c, i.e.
      ΔΔnet is the small residual of two large opposing ΔΔ terms and its
      sign is numerically ill-conditioned / unreliable.
    - ``dd_cancellation_r1`` (float | None): the r1 ratio (None if the
      ΔΔgas/ΔΔsolv magnitudes are both zero).
    - ``dd_cancellation_r2`` (float | None): companion r2 ratio.
    - ``dd_gas_kcal`` / ``dd_solv_kcal`` / ``dd_net_kcal`` (float): the
      ΔΔ decomposition echoed for the record.
    - ``sanity_gate_cancel_c`` (float): the cutoff used.
    - ``sanity_gate_reason`` (str | None): a phenomenological description
      when flagged; None otherwise.

    This is a sign-reliability tag only (ranking / cheap-triage regime); it
    makes no magnitude or absolute ΔG claim, and asserts NO root cause —
    in particular it does NOT attribute the cancellation to the force-field
    gas term; the root may be conformational sampling (a disengaged pose).

    Note: z = abs(dd_net) / SE_pooled, where it can be computed by the
    caller, is statistical-reproducibility context only and is NOT used by
    this gate (a high-z ΔΔnet can still be sign-ill-conditioned).

    Parameters
    ----------
    dd_gas, dd_solv, dd_net : float
        The ΔΔ(A − B) gas, solvation, and net (ΔΔG_bind) differences.
    c : float
        Cancellation cutoff; r1 < c flags. Defaults to
        ``SANITY_GATE_CANCEL_C`` (0.15).
    """
    r1, r2 = _cancellation_ratios(dd_gas, dd_solv, dd_net)

    out: Dict[str, Any] = {
        "dd_gas_kcal": dd_gas,
        "dd_solv_kcal": dd_solv,
        "dd_net_kcal": dd_net,
        "dd_cancellation_r1": r1,
        "dd_cancellation_r2": r2,
        "sanity_gate_cancel_c": c,
        "sign_ill_conditioned_cancellation": False,
        "sanity_gate_reason": None,
    }

    # No ΔΔ magnitude to assess (both ΔΔgas and ΔΔsolv vanish): cannot judge,
    # leave the flag False.
    if r1 is None:
        return out

    if r1 < c:
        out["sign_ill_conditioned_cancellation"] = True
        out["sanity_gate_reason"] = (
            f"ΔΔnet ({dd_net:+.1f}) is the small residual of ΔΔgas "
            f"({dd_gas:+.1f}) and ΔΔsolv ({dd_solv:+.1f}); cancellation "
            f"ratio r1={r1:.3f} < {c:.2f} (condition number "
            f">{1.0 / c:.1f}); ΔΔ sign numerically ill-conditioned — route "
            f"to expensive method. Phenomenological tag only: this is NOT a "
            f"force-field attribution; the root may be conformational "
            f"sampling (e.g. a disengaged binding pose)."
        )
    return out


def apply_pbsa_pair_sanity_gate_from_summaries(
    summary_a: Dict[str, Any],
    summary_b: Dict[str, Any],
    c: float = SANITY_GATE_CANCEL_C,
    gas_key: str = "mean_gas",
    solv_key: str = "mean_solv",
    net_key: str = "mean_dg",
) -> Dict[str, Any]:
    """Mode B convenience wrapper — compute ΔΔ from two cohort summaries.

    ΔΔ = summary_a − summary_b for the gas / solv / net cohort means, then
    delegate to ``apply_pbsa_pair_sanity_gate``. The cohort-mean keys must be
    present and non-None in both summaries; otherwise a ValueError is raised
    (the gate is undefined without a decomposition for both variants).

    Parameters
    ----------
    summary_a, summary_b : dict
        Cohort summaries for variant A (e.g. Cp4) and B (e.g. WT).
    gas_key, solv_key, net_key : str
        Keys for the cohort-mean gas, solvation, and net ΔG in each summary.
    """
    a_gas = summary_a.get(gas_key)
    a_solv = summary_a.get(solv_key)
    a_net = summary_a.get(net_key)
    b_gas = summary_b.get(gas_key)
    b_solv = summary_b.get(solv_key)
    b_net = summary_b.get(net_key)
    if None in (a_gas, a_solv, a_net, b_gas, b_solv, b_net):
        raise ValueError(
            "apply_pbsa_pair_sanity_gate_from_summaries requires "
            f"non-None '{gas_key}'/'{solv_key}'/'{net_key}' in both "
            "summaries (gate undefined without both decompositions)"
        )
    return apply_pbsa_pair_sanity_gate(
        dd_gas=float(a_gas) - float(b_gas),
        dd_solv=float(a_solv) - float(b_solv),
        dd_net=float(a_net) - float(b_net),
        c=c,
    )


def apply_pbsa_sanity_gate(result: Dict[str, Any]) -> Dict[str, Any]:
    """Per-snap / per-system endpoint logging (NO gate flag).

    A single endpoint is always deep in the cancellation regime
    (abs(net) ≪ abs(gas) + abs(solv)), so an endpoint-level cancellation
    test cannot discriminate trustworthy from untrustworthy results and is
    NOT a gate. This function therefore only LOGS the endpoint cancellation
    ratio for visibility; it never raises ``sign_ill_conditioned_*``.

    Non-destructive: the input dict is mutated in place and also returned,
    but no existing key is altered. Only additive keys are written:

    - ``endpoint_cancellation_ratio`` (float | None): r_A = abs(net) /
      (abs(gas) + abs(solv)) for this endpoint, logged for visibility.
      None when the gas/solv decomposition is missing or both vanish.
    - ``sign_invalid_gas_dominated`` (None): RETIRED Mode A flag, kept for
      schema continuity but always None (the per-snap gate was removed
      because per-snap cancellation is universal and undiscriminating).
    - ``sanity_gate_reason`` (None): no per-snap reason.

    The actual sign gate is ``apply_pbsa_pair_sanity_gate`` at the
    variant-pair ΔΔ level (Mode B). This is a ranking / cheap-triage logging
    helper only; it makes no magnitude or absolute ΔG claim.

    Parameters
    ----------
    result : dict
        May contain ``delta_g_kcal`` (net), ``delta_g_gas_kcal`` and
        ``delta_g_solv_kcal``. With any of those missing/None the ratio is
        logged as None.
    """
    net = result.get("delta_g_kcal")
    gas = result.get("delta_g_gas_kcal")
    solv = result.get("delta_g_solv_kcal")

    # RETIRED Mode A flag: kept for schema continuity, always None.
    result["sign_invalid_gas_dominated"] = None
    result["endpoint_cancellation_ratio"] = None
    result["sanity_gate_reason"] = None

    if net is None or gas is None or solv is None:
        return result

    denom = abs(gas) + abs(solv)
    if denom != 0.0:
        result["endpoint_cancellation_ratio"] = abs(net) / denom

    return result
