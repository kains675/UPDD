"""Adaptive λ-scheduling schedule I/O — registry-dict adapter + validated JSON.

Round-trips the schedule dict shape used by the asyncre launcher's
``_schedule_dict`` / ``get_schedule`` / ``write_cntl_file`` (lambdas / lambdas_1
/ lambdas_2 / directions / intermd / w0 / alpha / u0 / n_states), and emits a
validated ``adaptive_schedule.json`` (the proposed schedule + provenance:
seed, per-iteration overlap table, endpoint hashes, regime="ranking_only").

It also exposes ``load_schedule_dict`` — the loader the additive launcher flag
``--free-schedule-file`` uses to read a JSON schedule and feed it (verbatim) into
the SSOT ``write_cntl_file(schedule=...)``.

Python 3.8+; numpy + json (stdlib).
"""

import hashlib
import json
import os
from typing import Any, Dict, List, Optional

# The canonical per-state array keys (the _schedule_dict shape).
_SCHEDULE_KEYS = (
    "lambdas", "lambdas_1", "lambdas_2", "directions", "intermd",
    "w0", "alpha", "u0",
)
_INT_KEYS = ("directions", "intermd")


# ---------------------------------------------------------------------------
# Registry-dict round-trip
# ---------------------------------------------------------------------------
def to_registry_dict(schedule: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize a schedule dict into the registry shape (plain Python lists,
    ints for directions/intermd, floats elsewhere, n_states present).

    Drops any private ``_provenance`` payload (that travels in the JSON sidecar,
    not the registry-shaped dict). Raises KeyError if a required per-state array
    is missing or ValueError if array lengths disagree.
    """
    out: Dict[str, Any] = {}
    n: Optional[int] = None
    for key in _SCHEDULE_KEYS:
        if key not in schedule:
            if key == "lambdas":
                # lambdas defaults to lambdas_1 (canonical 22-state convention).
                if "lambdas_1" not in schedule:
                    raise KeyError(
                        "schedule missing both 'lambdas' and 'lambdas_1'"
                    )
                vals = list(schedule["lambdas_1"])
            else:
                raise KeyError(f"schedule missing required key {key!r}")
        else:
            vals = list(schedule[key])
        if key in _INT_KEYS:
            out[key] = [int(round(float(v))) for v in vals]
        else:
            out[key] = [float(v) for v in vals]
        if n is None:
            n = len(out[key])
        elif len(out[key]) != n:
            raise ValueError(
                f"schedule array {key!r} length {len(out[key])} != {n}"
            )
    out["n_states"] = int(n if n is not None else 0)
    return out


def from_registry_dict(d: Dict[str, Any]) -> Dict[str, Any]:
    """Inverse of ``to_registry_dict``: validate + return a schedule dict in the
    shape ``write_cntl_file(schedule=...)`` consumes. Idempotent with
    ``to_registry_dict`` (round-trip stable). Raises ValueError on a malformed
    dict (missing arrays / inconsistent lengths) so a corrupt JSON fails loud.
    """
    norm = to_registry_dict(d)
    # to_registry_dict already validates lengths + required keys.
    return norm


# ---------------------------------------------------------------------------
# Load a schedule dict from JSON (the --free-schedule-file loader)
# ---------------------------------------------------------------------------
def load_schedule_dict(path: str) -> Dict[str, Any]:
    """Load a schedule dict from a JSON file and return it in the
    ``write_cntl_file`` shape.

    The JSON may be either:
      * a bare schedule dict (the _schedule_dict shape), or
      * an ``adaptive_schedule.json`` emitted by ``emit_validated_schedule``
        (which nests the schedule under the ``"schedule"`` key) — this loader
        unwraps that automatically.

    Raises FileNotFoundError if absent, ValueError if malformed. This is the
    SOLE loader the launcher's ``--free-schedule-file`` uses; it feeds the
    result straight into the SSOT ``write_cntl_file(schedule=...)`` (no per-state
    array is re-derived).
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"schedule JSON not found: {path}")
    with open(path) as fh:
        raw = json.load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"schedule JSON must be an object, got {type(raw)}")
    # Unwrap an adaptive_schedule.json envelope.
    if "schedule" in raw and isinstance(raw["schedule"], dict):
        raw = raw["schedule"]
    return from_registry_dict(raw)


# ---------------------------------------------------------------------------
# Endpoint hashing (provenance / integrity)
# ---------------------------------------------------------------------------
def _endpoint_hashes(schedule: Dict[str, Any]) -> Dict[str, str]:
    """SHA-256 of the coupled (λ=0, forward first state) + apex (W0==1.0) state
    signatures. These are compared between seed and proposed to confirm the
    endpoints are byte-invariant (the empirical ranking-unbias check)."""
    directions = [float(d) for d in schedule["directions"]]
    fwd = sum(1 for d in directions if d >= 0)

    def sig(idx: int) -> str:
        tup = (
            float(schedule["lambdas_1"][idx]),
            float(schedule["lambdas_2"][idx]),
            float(schedule["w0"][idx]),
            float(schedule["alpha"][idx]),
            float(schedule["u0"][idx]),
            int(round(float(schedule["intermd"][idx]))),
        )
        return hashlib.sha256(repr(tup).encode("utf-8")).hexdigest()

    coupled_idx = 0
    apex_candidates = [i for i in range(fwd) if float(schedule["w0"][i]) == 1.0]
    hashes = {"coupled_lambda0": sig(coupled_idx)}
    if apex_candidates:
        hashes["apex_w0_1"] = sig(apex_candidates[0])
    return hashes


# ---------------------------------------------------------------------------
# Emit validated adaptive_schedule.json
# ---------------------------------------------------------------------------
def emit_validated_schedule(
    result: Dict[str, Any],
    path: str,
    seed: Optional[Dict[str, Any]] = None,
    overlap_table: Optional[List[Dict[str, Any]]] = None,
    extra_provenance: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Write ``adaptive_schedule.json`` = the final schedule dict + provenance.

    ``result`` — the proposed schedule dict (from ``propose_rebalance``); its
        embedded ``_provenance`` (if any) is merged into the JSON's provenance
        block and stripped from the registry-shaped schedule.
    ``seed`` — the seed schedule (its endpoint hashes are recorded for the
        integrity gate).
    ``overlap_table`` — the per-iteration / per-pair overlap rows (e.g.
        ``[p.as_dict() for p in pairs]``) for the audit trail.
    ``extra_provenance`` — any additional caller-supplied provenance (e.g. the
        pilot run_dir, warmup, host).

    Provenance ALWAYS carries ``regime="ranking_only"`` + endpoint hashes
    (seed + proposed) so a downstream consumer can verify the endpoints are
    invariant without re-running the guards. Returns the written envelope dict.
    """
    embedded_prov = {}
    if isinstance(result, dict) and "_provenance" in result:
        embedded_prov = dict(result["_provenance"])
    schedule_clean = {k: v for k, v in result.items() if k != "_provenance"}
    schedule_norm = to_registry_dict(schedule_clean)

    provenance: Dict[str, Any] = {
        "regime": "ranking_only",
        "n_states": schedule_norm["n_states"],
        "proposed_endpoint_hashes": _endpoint_hashes(schedule_norm),
    }
    if seed is not None:
        seed_norm = to_registry_dict(
            {k: v for k, v in seed.items() if k != "_provenance"}
        )
        provenance["seed_endpoint_hashes"] = _endpoint_hashes(seed_norm)
        provenance["seed_n_states"] = seed_norm["n_states"]
        provenance["endpoints_invariant"] = (
            provenance["seed_endpoint_hashes"]
            == provenance["proposed_endpoint_hashes"]
        )
        provenance["seed_schedule"] = seed_norm
    if overlap_table is not None:
        provenance["overlap_table"] = overlap_table
    if embedded_prov:
        provenance["rebalance"] = embedded_prov
    if extra_provenance:
        provenance.update(extra_provenance)

    envelope = {
        "schema": "adaptive_schedule.v1",
        "schedule": schedule_norm,
        "provenance": provenance,
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(envelope, fh, indent=2)
    return envelope
