#!/usr/bin/env python
"""
utils/parameterize_cofactor.py
------------------------------
[R-17 G3] Cofactor force-field parameterization gate.

Scope
~~~~~
This module implements the third preservation gate (G3) from SciVal's 8th
verdict (verdict_cofactor_preservation_policy_20260419.md §2.3). For each
required cofactor declared in ``target_card.cofactor_residues`` it:

  1. Resolves ``ff_parameters.frcmod`` / ``ff_parameters.mol2`` paths against
     the project root.
  2. If the files exist → cache-hit → return the abs paths.
  3. If the files do not exist and ``source == "library_insertion"`` with
     ``ff_parameters.method in {"gaff2_resp", "gaff2"}`` → attempt on-demand
     regeneration via antechamber + parmchk2 (graceful if the tool is
     missing — raises ``CofactorParamMissingError`` with a clear remediation
     message).
  4. For ions with ``method == "li_merz_12_6_4"`` (e.g. Mg²⁺) → parameters
     are shipped with AMBER (frcmod.ions1lm_126_iod / frcmod.ions_mg…) so
     no regeneration is needed; report cache-hit with path=None.

Interaction with other gates
~~~~~~~~~~~~~~~~~~~~~~~~~~~
- G2 (post-AF2 reinsertion) produces a PDB with cofactor HETATM records.
- G3 (this module) ensures the FF parameters are present at a deterministic
  location (``params/<resname>/<resname>.frcmod`` etc.) so G4
  (``run_restrained_md``) can call ``ForceField.loadFile(...)``.
- G4 consumes G3's output via ``run_restrained_md.load_cofactor_ff_parameters``.

The module intentionally does NOT invoke antechamber in CI / unit tests —
regeneration is only attempted when ``UPDD_ALLOW_ANTECHAMBER=1`` is set. In
all other cases, an existing cache is the only path to success. This avoids
silently generating parameters with unvalidated charges and preserves
scientific reproducibility.

Python 3.8+ compatible.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from typing import Dict, List, Optional

from cofactor_errors import CofactorParamMissingError


# ==========================================
# Path helpers
# ==========================================
_PROJ_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _resolve_abs(path: Optional[str]) -> Optional[str]:
    """Resolve a possibly-relative parameter path against the project root."""
    if not path:
        return None
    return path if os.path.isabs(path) else os.path.join(_PROJ_ROOT, path)


# ==========================================
# Antechamber wrapper (opt-in via UPDD_ALLOW_ANTECHAMBER=1)
# ==========================================
def _can_run_antechamber() -> bool:
    """Return True iff antechamber is on PATH and user explicitly opted in."""
    if os.environ.get("UPDD_ALLOW_ANTECHAMBER", "0") != "1":
        return False
    return bool(shutil.which("antechamber")) and bool(shutil.which("parmchk2"))


def _regenerate_gaff2_params(
    resname: str,
    source_pdb: str,
    out_dir: str,
    net_charge: int,
    method: str = "bcc",
) -> Dict[str, str]:
    """Run antechamber + parmchk2 to produce ``<resname>.mol2`` and ``.frcmod``.

    Args:
        resname: 3-letter residue name (UPPER).
        source_pdb: PDB with ONLY this residue (extracted upstream).
        out_dir: Directory to write the ``.mol2`` / ``.frcmod`` files into.
        net_charge: Integer total charge for the residue (passed to
            antechamber ``-nc``).
        method: Charge assignment method: ``"bcc"`` (AM1-BCC, default) or
            ``"resp"`` (requires external RESP server — NOT automated here).

    Returns:
        Dict with ``{"mol2": abs_path, "frcmod": abs_path}``.

    Raises:
        CofactorParamMissingError: Any subprocess failure. The caller (G3
            resolver) can decide whether to escalate or fallback.
    """
    if method not in ("bcc", "resp"):
        raise ValueError(f"Unsupported charge method '{method}' (bcc|resp).")
    if method == "resp":
        raise CofactorParamMissingError(
            f"RESP automation not implemented for {resname}. Supply "
            f"pre-computed .mol2/.frcmod via target_card.ff_parameters and "
            f"place them at {out_dir}/.",
            method=method,
        )

    os.makedirs(out_dir, exist_ok=True)
    mol2_path = os.path.join(out_dir, f"{resname}.mol2")
    frcmod_path = os.path.join(out_dir, f"{resname}.frcmod")

    try:
        subprocess.run(
            [
                "antechamber", "-i", source_pdb, "-fi", "pdb",
                "-o", mol2_path, "-fo", "mol2",
                "-c", "bcc", "-nc", str(net_charge),
                "-s", "0", "-rn", resname,
                "-at", "gaff2",
            ],
            check=True, capture_output=True, text=True, cwd=out_dir,
        )
        subprocess.run(
            ["parmchk2", "-i", mol2_path, "-f", "mol2", "-o", frcmod_path,
             "-s", "gaff2"],
            check=True, capture_output=True, text=True, cwd=out_dir,
        )
    except FileNotFoundError as exc:
        raise CofactorParamMissingError(
            f"antechamber/parmchk2 not available on PATH: {exc}",
            method=method, resname=resname, gate="G3",
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise CofactorParamMissingError(
            f"antechamber/parmchk2 failed for {resname}: rc={exc.returncode} "
            f"stderr={(exc.stderr or '')[:240]}",
            method=method, resname=resname, gate="G3",
        ) from exc

    return {"mol2": mol2_path, "frcmod": frcmod_path}


# ==========================================
# Public API — G3 entry point
# ==========================================
def ensure_cofactor_params(
    target_card: dict,
    source_pdbs: Optional[Dict[str, str]] = None,
) -> List[dict]:
    """[R-17 G3] Validate / regenerate force-field parameters for every
    required cofactor declared in ``target_card.cofactor_residues``.

    Args:
        target_card: Parsed target_card dict (must be v0.6.5 or auto-upgraded
            by ``utils_common.load_target_card``).
        source_pdbs: Optional mapping ``{resname: path_to_isolated_pdb}`` used
            as input for antechamber when regeneration is needed. Callers
            (e.g. the pipeline orchestrator) are responsible for producing
            these stub PDBs. If the mapping is missing and regeneration is
            required, ``CofactorParamMissingError`` is raised.

    Returns:
        List of dicts, one per required cofactor entry, with keys:
            - resname (str)
            - status ("cache_hit" | "regenerated" | "ion_builtin" | "missing")
            - method (str, from ff_parameters.method)
            - mol2 (Optional[str])
            - frcmod (Optional[str])

    Raises:
        CofactorParamMissingError: A required cofactor lacks valid
            parameters AND regeneration failed (or was not permitted).
    """
    if source_pdbs is None:
        source_pdbs = {}

    out: List[dict] = []
    for entry in (target_card.get("cofactor_residues") or []):
        if not isinstance(entry, dict) or not entry.get("required"):
            continue

        resname = str(entry.get("resname", "")).strip().upper()
        if not resname:
            continue
        ff_params = entry.get("ff_parameters") or {}
        method = str(ff_params.get("method", "")).strip().lower()

        mol2_abs = _resolve_abs(ff_params.get("mol2"))
        frcmod_abs = _resolve_abs(ff_params.get("frcmod"))

        # Ion models (Li-Merz 12-6-4 etc.) ship with AMBER and do not need a
        # per-target frcmod. Report as ion_builtin so G4 skips loadFile().
        if method in ("li_merz_12_6_4", "li_merz", "amber_ion", "ion_builtin"):
            out.append({
                "resname": resname, "status": "ion_builtin",
                "method": method, "mol2": None, "frcmod": None,
            })
            continue

        # Cache hit: both files exist.
        if mol2_abs and frcmod_abs and os.path.exists(mol2_abs) and os.path.exists(frcmod_abs):
            out.append({
                "resname": resname, "status": "cache_hit",
                "method": method, "mol2": mol2_abs, "frcmod": frcmod_abs,
            })
            continue

        # Regeneration branch (opt-in only).
        source_key = entry.get("source", "pdb_literal")
        if source_key == "library_insertion" and _can_run_antechamber():
            if resname not in source_pdbs or not os.path.exists(source_pdbs[resname]):
                raise CofactorParamMissingError(
                    f"R-17 G3: required cofactor '{resname}' has no cached "
                    f"frcmod/mol2 (looked at {frcmod_abs}, {mol2_abs}) and "
                    f"no source PDB was provided for regeneration.",
                    method=method, resname=resname, gate="G3",
                    param_paths=[p for p in (mol2_abs, frcmod_abs) if p],
                )
            out_dir = os.path.dirname(frcmod_abs) if frcmod_abs else os.path.join(
                _PROJ_ROOT, "params", resname.lower()
            )
            charge = int(entry.get("charge", 0))
            regen = _regenerate_gaff2_params(
                resname, source_pdbs[resname], out_dir, charge, method="bcc",
            )
            out.append({
                "resname": resname, "status": "regenerated",
                "method": method, "mol2": regen["mol2"], "frcmod": regen["frcmod"],
            })
            continue

        # No cache, no regeneration permitted → fail-fast.
        raise CofactorParamMissingError(
            f"R-17 G3: required cofactor '{resname}' lacks force-field "
            f"parameters. Checked: mol2={mol2_abs}, frcmod={frcmod_abs}. "
            f"Set UPDD_ALLOW_ANTECHAMBER=1 + source='library_insertion' to "
            f"regenerate, or ship the .mol2/.frcmod under "
            f"params/{resname.lower()}/.",
            method=method, resname=resname, gate="G3",
            param_paths=[p for p in (mol2_abs, frcmod_abs) if p],
        )

    return out


__all__ = ["ensure_cofactor_params"]
