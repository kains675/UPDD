#!/usr/bin/env python
"""
utils/cofactor_errors.py
------------------------
R-17 (Cofactor Preservation SSOT) error class hierarchy.

SciVal 8th verdict (verdict_cofactor_preservation_policy_20260419.md §4.2, §8.2):
the R-17 enforcement pipeline spans six gates (G1 preprocess ... G6 QM/MM build).
Each gate raises a *specific* subclass of ``PipelineError`` so Keeper's delegate
reports (COFACTOR-01 / PRESERVE-01 / SURVIVE-01) can narrate which stage failed
and how to remediate without parsing free text.

Hierarchy:

    PipelineError (base)
     ├── CofactorError (R-17 base)
     │    ├── CofactorMissingError         — G1 / G5 (atom absence)
     │    ├── CofactorChargeSumMismatch    — G6 (declared-vs-injected charge sum)
     │    ├── CofactorParamMissingError    — G3 / G4 (frcmod/mol2/xml not loadable)
     │    └── CofactorReinsertionError     — G2 (post-AF2 superposition failure)

Design notes
~~~~~~~~~~~~
- Each class carries structured attributes (``resname``, ``chain``, ``resnum``,
  ``gate``, ``declared``, ``computed``, ``detail``) so callers (Keeper) can
  emit deterministic JSON without pattern matching.
- ``__str__`` is informative enough for Runner's log output.
- R-15 ``ChargeDeclarationMismatch`` and R-16 ``BinderChargeMismatch`` stay in
  their existing module (``utils/charge_topology.py``) — R-17 is orthogonal
  (atom presence), not a charge-topology concern.

Python 3.8+ compatible (Optional[str], no PEP 604 unions).
"""
from typing import Any, List, Optional


# ==========================================
# Base — shared with the rest of the pipeline
# ==========================================
class PipelineError(Exception):
    """Base exception for UPDD pipeline stage failures.

    Kept minimal to avoid pulling in stage-specific imports. Any stage-level
    error subclass (R-15/R-16 charge, R-17 cofactor, R-13 parity) inherits
    from this so Runner / Keeper can catch "pipeline-level" failures with a
    single ``except PipelineError`` clause.
    """

    pass


# ==========================================
# R-17 base
# ==========================================
class CofactorError(PipelineError):
    """Base class for all R-17 cofactor preservation failures.

    Subclasses carry structured context for Keeper's gate reports. The
    ``gate`` attribute identifies which of G1..G6 raised the exception so
    upstream remediation is deterministic.
    """

    def __init__(
        self,
        message: str,
        gate: Optional[str] = None,
        resname: Optional[str] = None,
        chain: Optional[str] = None,
        resnum: Optional[int] = None,
        detail: Any = None,
    ):
        super().__init__(message)
        self.gate = gate
        self.resname = resname
        self.chain = chain
        self.resnum = resnum
        self.detail = detail

    def to_dict(self) -> dict:
        """Return a JSON-serializable dict for diagnostic logs."""
        return {
            "error":   self.__class__.__name__,
            "message": str(self),
            "gate":    self.gate,
            "resname": self.resname,
            "chain":   self.chain,
            "resnum":  self.resnum,
            "detail":  self.detail,
        }


# ==========================================
# Gate-specific subclasses
# ==========================================
class CofactorMissingError(CofactorError):
    """Declared (required=True) cofactor not present at a gate checkpoint.

    Raised by:
      - G1 (preprocess_target): declared cofactor absent from crystal PDB /
        removed by chain-filter / dropped by PDBFixer.
      - G5 (extract_snapshots): declared cofactor absent from MD snapshot
        (MD dropped it — usually due to G4 force-field template mismatch or
        G1 bug).
      - G6 (run_qmmm.validate_cofactor_presence): declared cofactor absent
        from the snapshot PDB passed to the QM/MM builder.
    """

    pass


class CofactorChargeSumMismatch(CofactorError):
    """Sum of declared cofactor charges differs from what was injected.

    Raised at G6 when ``sum(cofactor.charge for c in target_card)`` does not
    match the sum actually added to the MM point-charge set (integer charge
    comparison, tolerance < 1e-3 e). This protects against silent under-
    injection (e.g. MG charge+2 declared but only GNP atoms actually
    appended to the MM charges).
    """

    def __init__(
        self,
        message: str,
        declared: Optional[float] = None,
        computed: Optional[float] = None,
        **kwargs: Any,
    ):
        super().__init__(message, **kwargs)
        self.declared = declared
        self.computed = computed

    def to_dict(self) -> dict:
        payload = super().to_dict()
        payload["declared"] = self.declared
        payload["computed"] = self.computed
        return payload


class CofactorParamMissingError(CofactorError):
    """Cofactor force-field parameters (frcmod/mol2/xml) are not available.

    Raised by:
      - G3 (parameterize_cofactor): ``ff_parameters.frcmod`` / ``mol2`` path
        does not resolve *and* we cannot regenerate (antechamber missing,
        ``source != 'library_insertion'``, or RESP failed).
      - G4 (run_restrained_md): OpenMM ForceField.loadFile() fails for a
        declared cofactor's force-field XML/frcmod file.
    """

    def __init__(
        self,
        message: str,
        param_paths: Optional[List[str]] = None,
        method: Optional[str] = None,
        **kwargs: Any,
    ):
        super().__init__(message, **kwargs)
        self.param_paths = list(param_paths) if param_paths else []
        self.method = method

    def to_dict(self) -> dict:
        payload = super().to_dict()
        payload["param_paths"] = self.param_paths
        payload["method"] = self.method
        return payload


class CofactorReinsertionError(CofactorError):
    """Post-AF2 cofactor reinsertion (G2) failed.

    Raised when:
      - Crystal PDB not found / does not contain declared cofactor.
      - Superposition RMSD > 2.0 Å (wrong reference frame; unsafe to
        transplant cofactor coords into AF2 output).
      - Target-chain Cα count in AF2 output < 20 (structurally incomplete
        reference frame for Kabsch).
    """

    def __init__(
        self,
        message: str,
        rmsd_angstrom: Optional[float] = None,
        n_ca_reference: Optional[int] = None,
        **kwargs: Any,
    ):
        super().__init__(message, **kwargs)
        self.rmsd_angstrom = rmsd_angstrom
        self.n_ca_reference = n_ca_reference

    def to_dict(self) -> dict:
        payload = super().to_dict()
        payload["rmsd_angstrom"] = self.rmsd_angstrom
        payload["n_ca_reference"] = self.n_ca_reference
        return payload


__all__ = [
    "PipelineError",
    "CofactorError",
    "CofactorMissingError",
    "CofactorChargeSumMismatch",
    "CofactorParamMissingError",
    "CofactorReinsertionError",
]
