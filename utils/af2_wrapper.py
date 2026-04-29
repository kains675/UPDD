"""utils/af2_wrapper.py — AlphaFold2 / ColabFold subprocess wrapper.

Schema: ``stage12_pipeline/0.1`` (PR-NEW-D Day 3, 2026-04-29).

Stage 2 of the v0.7 de-novo discovery cascade. Predicts the structure of a
ProteinMPNN-emitted sequence and returns quality metrics (pLDDT, pTM, PAE)
that the orchestrator uses to filter weak designs out of subsequent
QM/MM-driven Stages 3-5.

The wrapper accepts ColabFold (``colabfold_batch``) or local AlphaFold2
(``run_alphafold.py``) binaries; mock fallback emits deterministic
quality metrics so the de-novo chain remains testable on CPU-only / no-GPU
checkouts.

Direct invocation::

    python -m utils.af2_wrapper --check
    python -m utils.af2_wrapper --run --sequence GLPDSAVTAAFEKMVERLG --recycles 3
"""

from __future__ import annotations

import argparse
import json
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional


SCHEMA_VERSION = "stage12_pipeline/0.1"
TOOL_NAME = "utils.af2_wrapper"

_ENV_BINARY = "AF2_PATH"
_BINARY_CANDIDATES = ("colabfold_batch", "run_alphafold.py", "alphafold")


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class AF2Error(RuntimeError):
    """Base class for AF2 wrapper errors."""


class AF2NotFoundError(AF2Error):
    """No AF2/ColabFold binary found and mock disabled."""


class AF2InvocationError(AF2Error):
    """Subprocess call failed or output unparseable."""


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class AF2Result:
    """Output of :meth:`AF2Wrapper.run_af2`.

    Attributes
    ----------
    pdb_path : str
        Predicted structure PDB path.
    plddt : float
        Mean pLDDT (0-100).
    plddt_per_residue : list[float]
        Per-residue pLDDT.
    pae_matrix : list[list[float]]
        Predicted Aligned Error matrix (LxL). Empty list when not produced.
    ptm : float | None
        Predicted TM-score (0-1) when available.
    sequence : str
        Input sequence (post-translation if any).
    recycles : int
        Number of AF2 recycles used.
    mock : bool
        True for mock results.
    elapsed_s : float
    schema : str
    binary : str | None
    """
    pdb_path: str
    plddt: float
    plddt_per_residue: List[float]
    pae_matrix: List[List[float]]
    ptm: Optional[float]
    sequence: str
    recycles: int = 3
    mock: bool = False
    elapsed_s: float = 0.0
    schema: str = SCHEMA_VERSION
    binary: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# Wrapper class
# ---------------------------------------------------------------------------

class AF2Wrapper:
    """Subprocess wrapper around ColabFold / local AF2."""

    def __init__(self, binary: Optional[str] = None, *,
                 allow_mock: bool = True,
                 timeout_s: Optional[int] = None) -> None:
        self._explicit_binary = binary
        self.allow_mock = allow_mock
        self.timeout_s = timeout_s

    # ------------------------------------------------------------------
    # Availability
    # ------------------------------------------------------------------

    def _resolve_binary(self) -> Optional[str]:
        if self._explicit_binary:
            return (str(self._explicit_binary)
                    if Path(self._explicit_binary).exists() else None)
        env = os.environ.get(_ENV_BINARY)
        if env and Path(env).exists():
            return env
        for cand in _BINARY_CANDIDATES:
            p = shutil.which(cand)
            if p:
                return p
        return None

    def check_availability(self) -> Dict[str, Any]:
        path = self._resolve_binary()
        if not path:
            return {
                "available": False, "version": None, "path": None,
                "schema": SCHEMA_VERSION, "tool": TOOL_NAME,
            }
        version = "unknown"
        try:
            cp = subprocess.run([path, "--version"], capture_output=True,
                                text=True, timeout=10)
            out = (cp.stdout or "") + (cp.stderr or "")
            if cp.returncode == 0 and out.strip():
                version = out.strip().splitlines()[0][:80]
        except Exception:  # noqa: BLE001
            pass
        return {
            "available": True, "version": version, "path": path,
            "schema": SCHEMA_VERSION, "tool": TOOL_NAME,
        }

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run_af2(self, sequence: str, *,
                msa: Optional[str] = None,
                templates: Optional[List[str]] = None,
                recycles: int = 3,
                output_dir: Optional[str] = None,
                ) -> AF2Result:
        """Predict structure for ``sequence``.

        Parameters
        ----------
        sequence : str
            Protein sequence (1-letter, uppercase). Whitespace stripped.
        msa : str | None
            Optional MSA / a3m file path.
        templates : list[str] | None
            Optional template PDB list.
        recycles : int
            Number of AF2 recycle iterations (default: 3).
        output_dir : str | None
            Output directory.

        Returns
        -------
        :class:`AF2Result`

        Raises
        ------
        AF2NotFoundError
            Binary missing and mock disabled.
        AF2InvocationError
            Subprocess failure.
        ValueError
            Empty / invalid sequence.
        """
        if not sequence:
            raise ValueError("sequence is required")
        seq = "".join(sequence.split()).upper()
        if not seq.isalpha():
            raise ValueError(f"sequence must be alphabetic; got {seq!r}")
        if recycles < 0:
            raise ValueError(f"recycles must be >= 0; got {recycles}")

        path = self._resolve_binary()
        if not path:
            if not self.allow_mock:
                raise AF2NotFoundError(
                    "AF2/ColabFold binary not found in $AF2_PATH or PATH."
                )
            return self._mock_run(sequence=seq, recycles=recycles,
                                  output_dir=output_dir)

        return self._real_run(
            binary=path, sequence=seq, msa=msa,
            templates=templates, recycles=recycles,
            output_dir=output_dir,
        )

    # ------------------------------------------------------------------
    # Real subprocess invocation
    # ------------------------------------------------------------------

    def _real_run(self, *, binary: str, sequence: str,
                  msa: Optional[str], templates: Optional[List[str]],
                  recycles: int,
                  output_dir: Optional[str]) -> AF2Result:
        out = Path(output_dir) if output_dir else Path(
            tempfile.mkdtemp(prefix="af2_"))
        out.mkdir(parents=True, exist_ok=True)

        # Write a single-record FASTA the binary can ingest.
        fa = out / "input.fasta"
        fa.write_text(f">design\n{sequence}\n")

        argv: List[str] = [binary, str(fa), str(out),
                           "--num-recycle", str(int(recycles))]
        if msa:
            argv += ["--msa-mode", "custom", "--custom-msa", msa]
        if templates:
            argv += ["--templates", ",".join(templates)]

        t0 = datetime.now(timezone.utc)
        try:
            cp = subprocess.run(argv, capture_output=True, text=True,
                                timeout=self.timeout_s, check=False)
        except subprocess.TimeoutExpired as exc:
            raise AF2InvocationError(
                f"AF2 timed out after {self.timeout_s}s"
            ) from exc
        elapsed = (datetime.now(timezone.utc) - t0).total_seconds()

        if cp.returncode != 0:
            tail = (cp.stderr or "")[-2048:]
            raise AF2InvocationError(
                f"AF2 exited {cp.returncode}; stderr tail:\n{tail}"
            )

        pdb = self._locate_output_pdb(out)
        if not pdb:
            raise AF2InvocationError(
                f"AF2 completed but no PDB found under {out}"
            )
        # ColabFold writes a JSON of scores adjacent to each PDB.
        plddt_mean, plddt_pr, pae, ptm = self._parse_scores(out, pdb)

        return AF2Result(
            pdb_path=str(pdb), plddt=plddt_mean,
            plddt_per_residue=plddt_pr, pae_matrix=pae,
            ptm=ptm, sequence=sequence, recycles=recycles,
            mock=False, elapsed_s=elapsed, binary=binary,
        )

    @staticmethod
    def _locate_output_pdb(out: Path) -> Optional[Path]:
        # Prefer rank-1 / unrelaxed_rank_001 / model_1 patterns.
        prefixes = ("rank_1", "rank_001", "ranked_0",
                    "unrelaxed_rank_1", "unrelaxed_rank_001",
                    "model_1")
        for pre in prefixes:
            hits = sorted(out.glob(f"{pre}*.pdb"))
            if hits:
                return hits[0]
        hits = sorted(out.rglob("*.pdb"))
        return hits[0] if hits else None

    @staticmethod
    def _parse_scores(out: Path, pdb: Path) -> tuple:
        # ColabFold score JSON co-located with PDB.
        for cand in (pdb.with_suffix(".json"),
                     out / "scores.json",
                     out / "rank_1_scores.json"):
            if cand.exists():
                try:
                    data = json.loads(cand.read_text())
                except Exception:  # noqa: BLE001
                    continue
                pl = data.get("plddt") or data.get("pLDDT") or []
                pl_mean = float(sum(pl) / len(pl)) if pl else 0.0
                pae = data.get("pae") or data.get("predicted_aligned_error") \
                      or []
                ptm = data.get("ptm") or data.get("pTM")
                return (pl_mean, list(pl), list(pae),
                        float(ptm) if ptm is not None else None)
        return (0.0, [], [], None)

    # ------------------------------------------------------------------
    # Mock fallback
    # ------------------------------------------------------------------

    def _mock_run(self, *, sequence: str, recycles: int,
                  output_dir: Optional[str]) -> AF2Result:
        """Synthesise a deterministic AF2 result.

        Strategy:
            * pLDDT mean ∈ [55, 90] — derived from SHA-256 of the sequence.
            * pLDDT per residue: mean ± small noise.
            * pTM ∈ [0.5, 0.9].
            * PAE matrix: empty (skipped for cost; tests check shape contract
              but not numeric values).
            * Writes a tiny placeholder PDB (3 CA polyglycine).
        """
        out = Path(output_dir) if output_dir else Path(
            tempfile.mkdtemp(prefix="af2_mock_"))
        out.mkdir(parents=True, exist_ok=True)

        h = hashlib.sha256(sequence.encode("utf-8")).digest()
        plddt_mean = 55.0 + (h[0] / 255.0) * 35.0     # 55-90
        ptm = 0.5 + (h[1] / 255.0) * 0.4              # 0.5-0.9
        L = len(sequence)
        plddt_pr = [round(plddt_mean + ((h[i % len(h)] - 128) / 128.0) * 2.0,
                          2) for i in range(L)]

        pdb = out / "mock_pred.pdb"
        # 3 CA polyglycine placeholder (matches RFdiffusion mock format).
        pdb.write_text(
            f"REMARK  Mock AF2 prediction L={L} pLDDT={plddt_mean:.1f}\n"
            f"ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 30.00           C\n"
            f"ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00 30.00           C\n"
            f"ATOM      3  CA  GLY A   3       7.600   0.000   0.000  1.00 30.00           C\n"
            f"TER       4      GLY A   3\n"
            f"END\n"
        )
        return AF2Result(
            pdb_path=str(pdb),
            plddt=round(plddt_mean, 2),
            plddt_per_residue=plddt_pr,
            pae_matrix=[],
            ptm=round(ptm, 3),
            sequence=sequence, recycles=recycles,
            mock=True, elapsed_s=0.0, binary=None,
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="utils.af2_wrapper",
        description=("AF2 / ColabFold subprocess wrapper "
                     "(Stage 2 of UPDD de-novo discovery)."),
    )
    p.add_argument("--check", action="store_true")
    p.add_argument("--run", action="store_true")
    p.add_argument("--sequence", default=None)
    p.add_argument("--msa", default=None)
    p.add_argument("--templates", default=None,
                   help="Comma-separated template PDB paths.")
    p.add_argument("--recycles", type=int, default=3)
    p.add_argument("--output-dir", default=None)
    p.add_argument("--no-mock", action="store_true")
    p.add_argument("--binary", default=None)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    w = AF2Wrapper(binary=args.binary, allow_mock=(not args.no_mock))
    if args.check:
        info = w.check_availability()
        print(json.dumps(info, indent=2))
        return 0 if info["available"] else 1
    if args.run:
        if not args.sequence:
            print("--run requires --sequence", file=sys.stderr)
            return 2
        templates = ([t.strip() for t in args.templates.split(",")
                      if t.strip()] if args.templates else None)
        try:
            res = w.run_af2(
                sequence=args.sequence, msa=args.msa,
                templates=templates, recycles=args.recycles,
                output_dir=args.output_dir,
            )
        except AF2Error as exc:
            print(f"[af2] {type(exc).__name__}: {exc}", file=sys.stderr)
            return 3
        # Drop the bulky pae_matrix from the CLI print to keep output sane.
        out = res.to_dict()
        if out.get("pae_matrix"):
            out["pae_matrix"] = f"<{len(out['pae_matrix'])}x... matrix>"
        print(json.dumps(out, indent=2))
        return 0
    print("[af2] no action; use --check or --run", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
