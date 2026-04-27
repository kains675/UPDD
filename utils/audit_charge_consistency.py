#!/usr/bin/env python
"""utils/audit_charge_consistency.py — Stage 1 historical JSON audit.

Scans a ``qmmm_results/`` directory for ``*_qmmm_topology.json`` files and
classifies each against the chemistry-true charge computed from its
snapshot PDB (:func:`charge_topology.compute_qm_net_charge_topology`).

Classification per SciVal verdict v2 §5:
    PASS          — declared total == computed total (magnitude OK).
    PARITY_ONLY   — parity matches but magnitudes differ (silent-pass risk).
    FAIL          — parity differs (would trip OddElectronError with the
                    correct charge — yet the run converged at wrong charge).
    SKIPPED       — JSON in an ``_archive/`` subdirectory, or missing fields.

CLI (dry-run by default — the ``--apply`` flag is the only way to move
any JSON file):

    python -m utils.audit_charge_consistency \
        --outputs-dir outputs/6WGN_cyclic_htc_NMA_10_20-25/qmmm_results/ \
        --target-card target_cards/6WGN.json \
        --report _audit/charge_consistency_20260419.json

Note: ``--apply`` is intentionally user-gated. The audit default is dry-run
so the operator can review the manifest before any file movement.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Allow running either as ``python -m utils.audit_charge_consistency`` or
# ``python utils/audit_charge_consistency.py``.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from charge_topology import (  # noqa: E402
    compute_qm_net_charge_topology,
    parse_pdb_atoms_lite,
)


# ==========================================
# Constants
# ==========================================
QMMM_JSON_SUFFIX = "_qmmm_topology.json"
ARCHIVE_DIR_NAME = "_archive"
WRONG_CHARGE_BUCKET = "v0.6_wrong_charge"

# Default ncAA code aliases (Track 4 schema regression — Phase α Step 4).
# Catches the MLE↔NML asymmetry: target_card writes the canonical
# registry code (NML) while mmgbsa_summary may report the PDB resname
# (MLE), and vice versa. Allowed pairs are bidirectional.
_DEFAULT_SCHEMA_ALIASES: Dict[str, str] = {
    "MLE": "NML",
    "NML": "MLE",
}


# ==========================================================================
# Track 4 — Schema regression check (Phase α Step 4 / Cron #34 outlier)
# ==========================================================================

class SchemaRegressionError(Exception):
    """Raised when summary.ncaa_element disagrees with target_card.binder.ncaa_code.

    Catches the Cron #34 outlier pattern (mmgbsa_summary.ncaa_element="none"
    while target_card.binder.ncaa_code="MTR"). The mismatch indicates the
    pipeline's summary writer dropped the ncAA identity on the way out
    despite a fully-defined target_card; downstream ranking would
    treat such a run as wild-type by mistake.
    """


def _read_json_safe(path: str) -> Dict[str, Any]:
    """Load a JSON; raise FileNotFoundError or json.JSONDecodeError on error."""
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _resolve_card_ncaa_code(card: Dict[str, Any]) -> Optional[str]:
    """Return the canonical ncAA code from a target_card.

    Looks for ``binder.ncaa_code`` first (Phase α schema), then falls back
    to legacy ``binder.ncaa`` and ``ncaa_code`` top-level fields.
    Returns None if no field is present.
    """
    binder = card.get("binder")
    if isinstance(binder, dict):
        for key in ("ncaa_code", "ncaa", "code"):
            v = binder.get(key)
            if isinstance(v, str) and v:
                return v
    for key in ("ncaa_code", "ncaa"):
        v = card.get(key)
        if isinstance(v, str) and v:
            return v
    return None


def _resolve_summary_ncaa_element(summary: Dict[str, Any]) -> Optional[str]:
    """Return ``summary.ncaa_element`` (or equivalent), preserving 'none'.

    The Cron #34 outlier pattern emits literal string ``"none"`` (lowercase)
    instead of dropping the field; we preserve that as a sentinel value
    so the caller can distinguish "missing" (None) from "explicitly none"
    (string sentinel).
    """
    for key in ("ncaa_element", "ncaa_code", "ncaa"):
        v = summary.get(key)
        if isinstance(v, str):
            return v
    return None


def _aliases_match(a: str, b: str, allow_aliases: Optional[Dict[str, str]]) -> bool:
    """Return True iff a and b are equal or paired in allow_aliases."""
    if a == b:
        return True
    if allow_aliases is None:
        return False
    return allow_aliases.get(a) == b or allow_aliases.get(b) == a


def schema_regression_check(
    summary_path: str,
    target_card_path: str,
    allow_aliases: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Validate that summary.ncaa_element matches target_card.binder.ncaa_code.

    Args:
        summary_path: path to mmgbsa_summary.json (or compatible).
        target_card_path: path to target_cards/<id>.json.
        allow_aliases: optional dict mapping equivalent codes; default uses
            ``_DEFAULT_SCHEMA_ALIASES`` (MLE↔NML).

    Returns:
        dict with keys ``status`` ('pass'|'fail'|'inconclusive'),
        ``expected``, ``actual``, ``reason``.

    Notes:
        Pure check — never raises. Caller wraps with ``SchemaRegressionError``
        when a hard fail is desired. Missing summary fields → INCONCLUSIVE
        (a warning, not a block) because legacy summary JSONs predate the
        ncaa_element field.
    """
    aliases = _DEFAULT_SCHEMA_ALIASES if allow_aliases is None else allow_aliases

    try:
        summary = _read_json_safe(summary_path)
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "status": "inconclusive",
            "expected": None,
            "actual": None,
            "reason": f"summary_read_error: {type(exc).__name__}: {exc}",
        }

    try:
        card = _read_json_safe(target_card_path)
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "status": "inconclusive",
            "expected": None,
            "actual": None,
            "reason": f"target_card_read_error: {type(exc).__name__}: {exc}",
        }

    expected = _resolve_card_ncaa_code(card)
    actual = _resolve_summary_ncaa_element(summary)

    if expected is None:
        return {
            "status": "inconclusive",
            "expected": None,
            "actual": actual,
            "reason": "target_card has no ncaa_code (binder.ncaa_code / ncaa)",
        }

    if actual is None:
        return {
            "status": "inconclusive",
            "expected": expected,
            "actual": None,
            "reason": "summary missing ncaa_element field (legacy schema)",
        }

    # Cron #34 outlier sentinel: literal 'none' string disagrees with any
    # populated card. Treat as FAIL.
    if actual.lower() == "none" and expected.lower() != "none":
        return {
            "status": "fail",
            "expected": expected,
            "actual": actual,
            "reason": (f"summary.ncaa_element={actual!r} disagrees with "
                       f"target_card.binder.ncaa_code={expected!r} "
                       f"(Cron #34 outlier pattern)"),
        }

    if _aliases_match(actual, expected, aliases):
        return {
            "status": "pass",
            "expected": expected,
            "actual": actual,
            "reason": ("identical" if actual == expected
                       else f"alias_match: {actual} ↔ {expected}"),
        }

    return {
        "status": "fail",
        "expected": expected,
        "actual": actual,
        "reason": (f"summary.ncaa_element={actual!r} disagrees with "
                   f"target_card.binder.ncaa_code={expected!r}"),
    }


# ==========================================================================
# Track 4 — PATCH-01 backbone N audit wrapper (Phase α Step 4 / SciVal V2.5)
# ==========================================================================

def _run_patch_01_audit(
    run_xmls: List[str],
    strict_mode: bool = False,
    manifest: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Apply PATCH-01 backbone N gate on each XML.

    Skips non-TRP-derived ncAAs (deferred to N-methyl overlay R&D / PR-6).
    The current `_apply_*_amber14_charge_patch` is gated to parent=TRP only;
    N-methyl class XMLs (NMA/NML/MEA/SAR/...) carry an unpatched GAFF2
    backbone N and would always fail PATCH-01. Marking them as "skipped"
    (not "failed") keeps the gate focused on the population it can speak to.

    Args:
        run_xmls: list of XML paths to audit.
        strict_mode: if True, raise Patch01BackboneNError on first failure.
        manifest: optional dict to populate with patch_01_failures /
            patch_01_skipped / patch_01_passes / patch_01_inconclusive
            entries (mutated in place).

    Returns:
        Summary dict {n_total, n_pass, n_fail, n_skip, n_inconclusive,
        failures, skipped, entries}.

    Raises:
        Patch01BackboneNError: in strict_mode, on the first failing XML.
    """
    # Lazy import — utils_draft is the source of truth and must not be
    # mutated; keeping the import local prevents accidental side-effects
    # on plain `import utils.audit_charge_consistency`.
    repo_root = Path(__file__).resolve().parent.parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from utils_draft.gates.patch_01_backbone_n import (  # noqa: E402
        check_xml,
        Patch01BackboneNError,
        _CODE_PARENT_CLASS,
    )

    if manifest is not None:
        manifest.setdefault("patch_01_failures", [])
        manifest.setdefault("patch_01_skipped", [])
        manifest.setdefault("patch_01_passes", [])
        manifest.setdefault("patch_01_inconclusive", [])

    entries: List[Dict[str, Any]] = []
    n_pass = 0
    n_fail = 0
    n_skip = 0
    n_inconc = 0

    for xml_path in run_xmls:
        p = Path(xml_path)
        # Derive code from filename to determine if scope-eligible.
        m = re.match(r"^(?P<code>[A-Z0-9]{3,4})_gaff2\.xml$", p.name)
        code = m.group("code") if m else None
        parent, cls = _CODE_PARENT_CLASS.get(code or "", (None, None))

        # Out-of-scope skip: anything whose parent is not TRP. This is the
        # PR-6 deferral envelope (N-methyl overlay R&D); the gate exists
        # but cannot speak to non-TRP class until that overlay lands.
        if parent != "TRP":
            entry = {
                "xml": xml_path,
                "code": code,
                "parent": parent,
                "class": cls,
                "status": "skipped",
                "reason": (f"non-TRP parent (parent={parent!r}, class={cls!r}) "
                           f"— deferred to N-methyl overlay R&D / PR-6"),
            }
            entries.append(entry)
            if manifest is not None:
                manifest["patch_01_skipped"].append(entry)
            n_skip += 1
            continue

        res = check_xml(p)
        entry = {
            "xml": xml_path,
            "code": code,
            "parent": parent,
            "class": cls,
            "primary_residue": res.primary_residue,
            "threshold": res.threshold,
            "backbone_n_charge": res.backbone_n_charge,
            "amber14_ref_n": res.amber14_ref_n,
            "delta_q": res.delta_q,
            "passed": bool(res.passed),
            "inconclusive": bool(res.inconclusive),
            "status": ("pass" if res.passed
                       else "inconclusive" if res.inconclusive
                       else "fail"),
            "reason": res.reason,
        }
        entries.append(entry)

        if res.passed:
            n_pass += 1
            if manifest is not None:
                manifest["patch_01_passes"].append(entry)
            continue

        if res.inconclusive:
            n_inconc += 1
            if manifest is not None:
                manifest["patch_01_inconclusive"].append(entry)
            if strict_mode:
                raise Patch01BackboneNError(
                    f"PATCH-01 INCONCLUSIVE (strict_mode): {res.reason} "
                    f"({xml_path})"
                )
            continue

        # Hard fail.
        n_fail += 1
        if manifest is not None:
            manifest["patch_01_failures"].append(entry)
        if strict_mode:
            raise Patch01BackboneNError(
                f"Backbone N charge deviation {abs(res.delta_q):.4f} e > "
                f"{res.threshold:.2f} e threshold for ncAA "
                f"{res.primary_residue or '?'} ({xml_path}). "
                f"Run Phase α patch first."
            )

    return {
        "n_total": len(run_xmls),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "n_skip": n_skip,
        "n_inconclusive": n_inconc,
        "failures": [e for e in entries if e["status"] == "fail"],
        "skipped": [e for e in entries if e["status"] == "skipped"],
        "passes": [e for e in entries if e["status"] == "pass"],
        "inconclusive": [e for e in entries if e["status"] == "inconclusive"],
    }


def _classify(declared_total: int, computed_total: int) -> str:
    """Classify per SciVal verdict v2 §5."""
    if declared_total == computed_total:
        return "PASS"
    if (declared_total % 2) == (computed_total % 2):
        return "PARITY_ONLY"
    return "FAIL"


def _in_archive(path: Path) -> bool:
    """Return True if any component of ``path`` is the archive sentinel."""
    return any(part == ARCHIVE_DIR_NAME for part in path.parts)


def _find_snapshot_pdb(json_obj: Dict[str, Any], fallback_dir: Path) -> Optional[Path]:
    """Locate the snapshot PDB for a qmmm JSON record.

    Preference order:
        1. ``pdb_path`` field (absolute).
        2. Derive from ``snapshot`` field + sibling ``snapshots/`` directory.
    """
    pdb_path = json_obj.get("pdb_path")
    if pdb_path and os.path.isfile(pdb_path):
        return Path(pdb_path)

    snapshot = json_obj.get("snapshot")
    if not snapshot:
        return None

    # Walk up from qmmm_results/ to the project run dir, then look at snapshots/
    candidate = fallback_dir.parent / "snapshots" / f"{snapshot}.pdb"
    if candidate.is_file():
        return candidate
    # Fallback: scan fallback_dir.parent for *.pdb matching the snapshot stem
    for p in fallback_dir.parent.rglob(f"{snapshot}.pdb"):
        return p
    return None


def _load_target_card(target_card_path: Optional[Path]) -> Dict[str, Any]:
    """Load target card JSON, with a permissive default if missing."""
    if target_card_path is None:
        return {}
    try:
        with open(target_card_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"[WARN] could not load target card {target_card_path}: {exc}",
              file=sys.stderr)
        return {}


def audit_one_json(
    json_path: Path,
    target_card: Dict[str, Any],
    ncaa_charge_map: Optional[Dict[str, int]] = None,
) -> Dict[str, Any]:
    """Classify a single qmmm_topology.json file.

    Returns:
        dict entry suitable for inclusion in the audit manifest.
    """
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            rec = json.load(f)
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "json_path": str(json_path),
            "status": "SKIPPED",
            "reason": f"json_read_error: {exc}",
        }

    snapshot = rec.get("snapshot") or json_path.stem.replace("_qmmm_topology", "")

    # Skip failure results (partition_error / odd_electron / etc.) — those
    # already carry a status field. The audit is for *converged* runs that
    # may have silently passed with a wrong declared charge.
    status_in_record = rec.get("status")
    if status_in_record == "FAILED":
        return {
            "json_path": str(json_path),
            "snapshot": snapshot,
            "status": "SKIPPED",
            "reason": (f"original_run_failed: {rec.get('reason')}"),
        }

    pdb_path = _find_snapshot_pdb(rec, json_path.parent)
    if pdb_path is None:
        return {
            "json_path": str(json_path),
            "snapshot": snapshot,
            "status": "SKIPPED",
            "reason": "snapshot_pdb_not_found",
            "pdb_path_declared": rec.get("pdb_path"),
        }

    # Resolve binder/target chains and contact list from target_card if
    # available, else from JSON record.
    binder_chain = target_card.get("binder_chain") or rec.get("binder_chain") or "B"
    target_chain = target_card.get("target_chain") or rec.get("target_chain") or "A"
    contacts = target_card.get("target_contact_residues") or []
    whole_excs = target_card.get("whole_residue_exceptions") or []

    try:
        atoms = parse_pdb_atoms_lite(str(pdb_path))
        binder_q, target_q, total_q, diag = compute_qm_net_charge_topology(
            atoms,
            binder_chain=binder_chain,
            target_chain=target_chain,
            target_contact_residues=contacts,
            whole_residue_exceptions=whole_excs,
            pH=7.4,
            ncaa_charge_map=ncaa_charge_map,
        )
    except Exception as exc:  # noqa: BLE001 — audit must never crash
        return {
            "json_path": str(json_path),
            "snapshot": snapshot,
            "status": "SKIPPED",
            "reason": f"chemistry_compute_error: {exc}",
            "pdb_path": str(pdb_path),
        }

    # Declared charge: prefer explicit JSON fields (Stage 1 adds these). If
    # absent (historical JSON), fall back to target_card defaults.
    declared_total = rec.get("qm_net_charge_declared")
    declared_binder = rec.get("binder_charge_declared")
    if declared_total is None:
        # Legacy path: reconstruct declared total from card fields.
        card_binder = int(target_card.get("binder_net_charge", 0))
        card_target = int(target_card.get("target_iso_net_charge", 0))
        declared_total = card_binder + card_target
        if declared_binder is None:
            declared_binder = card_binder

    declared_total = int(declared_total)
    declared_binder = int(declared_binder if declared_binder is not None else 0)

    parity_ok = (declared_total % 2) == (total_q % 2)
    magnitude_ok = declared_total == total_q
    status = _classify(declared_total, total_q)

    return {
        "json_path": str(json_path),
        "snapshot": snapshot,
        "pdb_path": str(pdb_path),
        "declared": declared_total,
        "binder_declared": declared_binder,
        "binder_chem": int(binder_q),
        "target_chem": int(target_q),
        "total_chem": int(total_q),
        "parity_ok": bool(parity_ok),
        "magnitude_ok": bool(magnitude_ok),
        "status": status,
        "cyclic": bool(diag.get("binder_diag", {}).get("cyclic", False)),
    }


def _collect_json_paths(outputs_dir: Path) -> List[Path]:
    """Find all qmmm topology JSONs under ``outputs_dir`` excluding ``_archive``."""
    paths = []
    for p in sorted(outputs_dir.rglob(f"*{QMMM_JSON_SUFFIX}")):
        if _in_archive(p):
            continue
        paths.append(p)
    return paths


def _collect_param_xmls(outputs_dir: Path) -> List[Path]:
    """Find all *_gaff2.xml ncAA params under ``outputs_dir`` excluding ``_archive``."""
    paths = []
    for p in sorted(outputs_dir.rglob("*_gaff2.xml")):
        if _in_archive(p):
            continue
        if p.parent.name != "params":
            continue
        paths.append(p)
    return paths


def run_audit(
    outputs_dir: Path,
    report_path: Path,
    target_card_path: Optional[Path] = None,
    apply_moves: bool = False,
    run_patch_01: bool = False,
    patch_01_strict: bool = False,
) -> Dict[str, Any]:
    """Run the full audit and write a manifest JSON.

    The base R-15/R-16 audit is unchanged. When ``run_patch_01=True``,
    a non-disruptive PATCH-01 backbone N gate sweep is added on every
    ``*_gaff2.xml`` under ``outputs_dir/<run>/params/``. Non-TRP-class
    XMLs are auto-skipped (PR-6 deferral). The base audit's PASS counts
    are NOT affected by PATCH-01.
    """
    target_card = _load_target_card(target_card_path)
    json_paths = _collect_json_paths(outputs_dir)

    entries: List[Dict[str, Any]] = []
    for p in json_paths:
        entries.append(audit_one_json(p, target_card))

    summary: Dict[str, int] = {"PASS": 0, "PARITY_ONLY": 0, "FAIL": 0, "SKIPPED": 0}
    for e in entries:
        summary[e.get("status", "SKIPPED")] = summary.get(e.get("status", "SKIPPED"), 0) + 1

    manifest = {
        "schema_version": "1.1",
        "outputs_dir": str(outputs_dir),
        "target_card": str(target_card_path) if target_card_path else None,
        "dry_run": not apply_moves,
        "summary": summary,
        "entries": entries,
    }

    # PATCH-01 add-on (Phase α Step 4 — non-disruptive). Failures/skips
    # are recorded in the manifest under separate keys; they do NOT change
    # the R-15/R-16 PASS/PARITY_ONLY/FAIL classification.
    if run_patch_01:
        xmls = [str(x) for x in _collect_param_xmls(outputs_dir)]
        patch01_summary = _run_patch_01_audit(
            xmls,
            strict_mode=patch_01_strict,
            manifest=manifest,
        )
        manifest["patch_01_summary"] = patch01_summary

    # Dry-run plan: entries that *would* be moved.
    plan: List[Dict[str, str]] = []
    archive_root = outputs_dir / ARCHIVE_DIR_NAME / WRONG_CHARGE_BUCKET
    for e in entries:
        if e.get("status") in ("FAIL", "PARITY_ONLY"):
            src = Path(e["json_path"])
            dest = archive_root / src.name
            plan.append({"src": str(src), "dest": str(dest), "status": e["status"]})
    manifest["move_plan"] = plan

    # Materialise the manifest
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    # Execute moves only when --apply is passed.
    moved: List[Dict[str, str]] = []
    if apply_moves and plan:
        archive_root.mkdir(parents=True, exist_ok=True)
        import shutil
        for item in plan:
            try:
                shutil.move(item["src"], item["dest"])
                moved.append(item)
            except OSError as exc:
                print(f"[WARN] move failed {item['src']}: {exc}", file=sys.stderr)
        manifest["moved"] = moved
        with open(report_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)

    return manifest


def _print_summary(manifest: Dict[str, Any]) -> None:
    summary = manifest.get("summary", {})
    total = sum(summary.values())
    print("=" * 70)
    print(f"Audit summary ({'APPLIED' if not manifest.get('dry_run') else 'DRY-RUN'})")
    print(f"  outputs_dir: {manifest.get('outputs_dir')}")
    print(f"  target_card: {manifest.get('target_card')}")
    print(f"  total JSON scanned: {total}")
    for k in ("PASS", "PARITY_ONLY", "FAIL", "SKIPPED"):
        print(f"    {k:>12}: {summary.get(k, 0)}")
    print("=" * 70)
    plan = manifest.get("move_plan", [])
    if plan:
        action = "Moved" if not manifest.get("dry_run") else "Would move"
        print(f"\n{action} to _archive/{WRONG_CHARGE_BUCKET}/:")
        for item in plan:
            print(f"  [{item['status']:>11}] {Path(item['src']).name}")
    else:
        print("\n(no FAIL / PARITY_ONLY entries)")

    entries = manifest.get("entries", [])
    # Show per-entry classification for transparency
    print("\nPer-entry classification:")
    for e in entries:
        st = e.get("status", "?")
        name = Path(e.get("json_path", "?")).name
        extra = ""
        if st in ("PASS", "PARITY_ONLY", "FAIL"):
            extra = (f" decl={e.get('declared')}, chem={e.get('total_chem')} "
                     f"(binder {e.get('binder_chem')}, target {e.get('target_chem')})")
        elif st == "SKIPPED":
            extra = f" reason={e.get('reason')}"
        print(f"  [{st:>11}] {name}{extra}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="UPDD charge consistency audit (Stage 1, R-15/R-16).",
    )
    parser.add_argument(
        "--outputs-dir", required=True,
        help="Directory containing *_qmmm_topology.json results (recursive).",
    )
    parser.add_argument(
        "--target-card", default=None,
        help="target_cards/<target_id>.json to resolve chains & contacts "
             "(default: auto-detect from first JSON).",
    )
    parser.add_argument(
        "--report", required=True,
        help="Output manifest JSON path.",
    )
    parser.add_argument(
        "--apply", action="store_true",
        help="Actually move FAIL/PARITY_ONLY JSONs into _archive/. "
             "Default is dry-run (manifest only).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="[default behaviour] Do not move files; write manifest only.",
    )
    parser.add_argument(
        "--run-patch-01", action="store_true",
        help="Also run the PATCH-01 backbone N gate over <run>/params/*_gaff2.xml. "
             "Non-TRP-derived ncAAs are auto-skipped (PR-6 deferral).",
    )
    parser.add_argument(
        "--patch-01-strict", action="store_true",
        help="With --run-patch-01: raise on first PATCH-01 failure (rc!=0).",
    )
    args = parser.parse_args()

    outputs_dir = Path(args.outputs_dir).resolve()
    report_path = Path(args.report).resolve()
    target_card_path = Path(args.target_card).resolve() if args.target_card else None

    if not outputs_dir.is_dir():
        print(f"[ERROR] outputs dir not found: {outputs_dir}", file=sys.stderr)
        return 2

    # Auto-detect the target_card when not provided by reading target_card_id
    # from the first candidate JSON.
    if target_card_path is None:
        candidates = _collect_json_paths(outputs_dir)
        for p in candidates:
            try:
                with open(p, "r", encoding="utf-8") as f:
                    rec = json.load(f)
                tid = rec.get("target_card_id")
                if tid:
                    guess = Path(__file__).resolve().parent.parent / "target_cards" / f"{tid}.json"
                    if guess.is_file():
                        target_card_path = guess
                        print(f"[INFO] auto-detected target card: {guess}")
                        break
            except (OSError, json.JSONDecodeError):
                continue

    # Resolve dry-run / apply flags. --apply wins; --dry-run is informational.
    apply_moves = bool(args.apply)

    manifest = run_audit(outputs_dir, report_path, target_card_path,
                         apply_moves=apply_moves,
                         run_patch_01=bool(args.run_patch_01),
                         patch_01_strict=bool(args.patch_01_strict))
    _print_summary(manifest)
    print(f"\n[INFO] manifest written: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
