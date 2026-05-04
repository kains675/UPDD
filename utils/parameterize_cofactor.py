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
    # Filenames are lowercase to match target_card's canonical layout
    # ``params/<resname_lower>/<resname_lower>.{mol2,frcmod}``. The internal
    # residue name inside the files stays uppercase via antechamber ``-rn``.
    _stem = resname.lower()
    mol2_path = os.path.join(out_dir, f"{_stem}.mol2")
    frcmod_path = os.path.join(out_dir, f"{_stem}.frcmod")

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
        # Also emit (a) an SDF and (b) a SYBYL-typed mol2 alongside the GAFF2
        # mol2. The SYBYL mol2 is what OpenFF/RDKit can actually parse —
        # RDKit's Mol2 parser rejects GAFF2 atom type codes like 'p5'/'n7',
        # while SYBYL codes ('P.3'/'N.am') are standard. The SYBYL file
        # preserves the AM1-BCC partial charges from the GAFF2 file, so G4
        # can feed charge-carrying OpenFF Molecules to GAFFTemplateGenerator
        # without needing to re-run AM1-BCC.
        sdf_path = os.path.join(out_dir, f"{_stem}.sdf")
        sybyl_mol2_path = os.path.join(out_dir, f"{_stem}_sybyl.mol2")
        try:
            subprocess.run(
                ["antechamber", "-i", mol2_path, "-fi", "mol2",
                 "-o", sdf_path, "-fo", "sdf",
                 "-at", "gaff2", "-rn", resname, "-s", "0", "-dr", "no"],
                check=True, capture_output=True, text=True, cwd=out_dir,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            pass
        try:
            subprocess.run(
                ["antechamber", "-i", mol2_path, "-fi", "mol2",
                 "-o", sybyl_mol2_path, "-fo", "mol2",
                 "-at", "sybyl", "-rn", resname, "-s", "0", "-dr", "no"],
                check=True, capture_output=True, text=True, cwd=out_dir,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            pass
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
# Method-token classification (SciVal verdict_cofactor_method_catalog)
# ==========================================
# Tokens that resolve to AMBER-shipped ion FF (no per-target mol2/frcmod).
# G4 consumers skip loadFile() for these resnames.
_ION_BUILTIN_TOKENS = frozenset({
    "li_merz_12_6_4",
    "li_merz",
    "li_merz_12_6_4_trivalent",  # Fe3+ (frcmod.ions34lm_1264)
    "li_merz_126",               # Cu+ monovalent transition metal
    "panteva_12_6_4",            # Mg2+/Mn2+/Zn2+ nucleic-acid context
    "joung_cheatham",            # Na/K/Cs/Li/Cl/Br/I/F
    "amber_ion",
    "ion_builtin",
})

# Tokens that require pre-shipped community frcmod/mol2 (Meagher polyphosphate,
# Ryde NAD, Shahrokh heme, Aleksandrov CoA, Kashefolgheta oxoanions, …).
# G3 does NOT regenerate via antechamber for these — user must ship the files
# under params/<resname>/ or set ff_parameters.{mol2,frcmod} to their paths.
_LIBRARY_TOKENS = frozenset({
    "amber_library",
    "amber_library_heme",
    "amber_library_heme_covalent",
    "gaff2_anion",
})

# Tokens eligible for antechamber AM1-BCC auto-regeneration from a source PDB.
_GAFF2_AUTO_TOKENS = frozenset({"gaff2", "gaff2_bcc"})

# Tokens that require user-supplied RESP-fit mol2. G3 verifies the file
# exists (at ff_parameters.resp_mol2 or .mol2) and runs parmchk2 against it.
_GAFF2_RESP_EXTERNAL_TOKENS = frozenset({"gaff2_resp", "gaff2_resp_external"})


def _library_remediation_hint(resname: str, method: str) -> str:
    """One-line hint appended to CofactorParamMissingError messages pointing
    the user at the canonical community library source for ``method``."""
    hints = {
        "amber_library":
            "Ship community frcmod/mol2 (Meagher 2003 polyphosphate / Ryde 1995 NAD / "
            "Aleksandrov 2006 CoA / Bryce DB). See verdict_cofactor_method_catalog §2.",
        "amber_library_heme":
            "Ship Shahrokh 2012 heme library for the declared heme_state; "
            "set ff_parameters.heme_state and supply the matching .frcmod/.mol2.",
        "amber_library_heme_covalent":
            "Run MCPB.py per-target to generate the bonded Fe-S-C frcmod (Li/Merz 2016).",
        "gaff2_anion":
            "Ship Kashefolgheta 2017 oxoanion frcmod/mol2.",
        "gaff2_resp_external":
            "Run Gaussian/ORCA + antechamber -c resp and place the .mol2 at "
            "ff_parameters.resp_mol2 (or .mol2). parmchk2 will derive the frcmod.",
    }
    return hints.get(method, "")


def _parmchk2_from_resp_mol2(
    resname: str, resp_mol2: str, out_dir: str,
) -> Dict[str, str]:
    """Run parmchk2 against a user-supplied RESP .mol2 to derive the frcmod.

    Does NOT touch charges (RESP file is authoritative). Raises
    CofactorParamMissingError if parmchk2 is unavailable or fails.
    """
    os.makedirs(out_dir, exist_ok=True)
    frcmod_path = os.path.join(out_dir, f"{resname.lower()}.frcmod")
    try:
        subprocess.run(
            ["parmchk2", "-i", resp_mol2, "-f", "mol2", "-o", frcmod_path,
             "-s", "gaff2"],
            check=True, capture_output=True, text=True, cwd=out_dir,
        )
    except FileNotFoundError as exc:
        raise CofactorParamMissingError(
            f"parmchk2 not available on PATH for {resname} "
            f"(gaff2_resp_external): {exc}",
            method="gaff2_resp_external", resname=resname, gate="G3",
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise CofactorParamMissingError(
            f"parmchk2 failed for {resname} (gaff2_resp_external): "
            f"rc={exc.returncode} stderr={(exc.stderr or '')[:240]}",
            method="gaff2_resp_external", resname=resname, gate="G3",
        ) from exc
    return {"mol2": os.path.abspath(resp_mol2), "frcmod": frcmod_path}


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
        resp_mol2_abs = _resolve_abs(ff_params.get("resp_mol2"))

        # (1) Ion models ship with AMBER and do not need a per-target frcmod.
        # See verdict_cofactor_method_catalog §2.1 for the token catalog.
        if method in _ION_BUILTIN_TOKENS:
            out.append({
                "resname": resname, "status": "ion_builtin",
                "method": method, "mol2": None, "frcmod": None,
            })
            continue

        # (2) Cache hit: both files exist.
        if mol2_abs and frcmod_abs and os.path.exists(mol2_abs) and os.path.exists(frcmod_abs):
            out.append({
                "resname": resname, "status": "cache_hit",
                "method": method, "mol2": mol2_abs, "frcmod": frcmod_abs,
            })
            continue

        # (3) Community library tokens — user MUST ship frcmod/mol2.
        # No antechamber regen (polyphosphate / heme / oxoanion cannot be
        # reproduced via AM1-BCC). Fail-fast with remediation hint.
        if method in _LIBRARY_TOKENS:
            raise CofactorParamMissingError(
                f"R-17 G3: cofactor '{resname}' uses method='{method}' which "
                f"requires a pre-shipped community library. Checked: "
                f"mol2={mol2_abs}, frcmod={frcmod_abs}. "
                f"{_library_remediation_hint(resname, method)}",
                method=method, resname=resname, gate="G3",
                param_paths=[p for p in (mol2_abs, frcmod_abs) if p],
            )

        # (4) User-supplied RESP .mol2 path → run parmchk2 to derive frcmod.
        if method in _GAFF2_RESP_EXTERNAL_TOKENS and resp_mol2_abs and \
                os.path.exists(resp_mol2_abs):
            out_dir = os.path.dirname(frcmod_abs) if frcmod_abs else os.path.join(
                _PROJ_ROOT, "params", resname.lower()
            )
            regen = _parmchk2_from_resp_mol2(resname, resp_mol2_abs, out_dir)
            out.append({
                "resname": resname, "status": "regenerated_resp_external",
                "method": method, "mol2": regen["mol2"], "frcmod": regen["frcmod"],
            })
            continue

        # (5) AM1-BCC regeneration branch (opt-in only). Accepts both
        # gaff2/gaff2_bcc and gaff2_resp/gaff2_resp_external as fallback when
        # no resp_mol2 is shipped — the caller gets a 'regenerated_bcc_fallback'
        # status to signal the RESP step is deferred.
        source_key = entry.get("source", "pdb_literal")
        if (method in _GAFF2_AUTO_TOKENS or method in _GAFF2_RESP_EXTERNAL_TOKENS) \
                and _can_run_antechamber():
            if resname not in source_pdbs or not os.path.exists(source_pdbs[resname]):
                raise CofactorParamMissingError(
                    f"R-17 G3: required cofactor '{resname}' has no cached "
                    f"frcmod/mol2 (looked at {frcmod_abs}, {mol2_abs}) and "
                    f"no source PDB was provided for regeneration "
                    f"(source='{source_key}', method='{method}'). Pass "
                    f"source_pdbs={{'{resname}': path}} or ship the params "
                    f"manually.",
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
            status = ("regenerated_bcc_fallback"
                      if method in _GAFF2_RESP_EXTERNAL_TOKENS
                      else "regenerated")
            out.append({
                "resname": resname, "status": status,
                "method": method, "mol2": regen["mol2"], "frcmod": regen["frcmod"],
            })
            continue

        # (6) No cache, unknown / unsupported method token → fail-fast.
        raise CofactorParamMissingError(
            f"R-17 G3: required cofactor '{resname}' (method='{method}') "
            f"lacks force-field parameters. Checked: mol2={mol2_abs}, "
            f"frcmod={frcmod_abs}. Supported tokens: ion_builtin="
            f"{sorted(_ION_BUILTIN_TOKENS)}, library={sorted(_LIBRARY_TOKENS)}, "
            f"auto-regen={sorted(_GAFF2_AUTO_TOKENS | _GAFF2_RESP_EXTERNAL_TOKENS)}. "
            f"Set UPDD_ALLOW_ANTECHAMBER=1 (auto-regen path) or ship the "
            f".mol2/.frcmod under params/{resname.lower()}/. "
            f"{_library_remediation_hint(resname, method)}",
            method=method, resname=resname, gate="G3",
            param_paths=[p for p in (mol2_abs, frcmod_abs) if p],
        )

    return out


__all__ = ["ensure_cofactor_params"]
