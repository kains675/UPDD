"""utils/proteinmpnn_wrapper.py — ProteinMPNN subprocess wrapper.

Schema: ``stage12_pipeline/0.1`` (PR-NEW-D Day 3, 2026-04-29).

Stage 1.5 of the v0.7 de-novo discovery cascade. Given an RFdiffusion
backbone PDB, propose ``num_sequences`` sequences scored by ProteinMPNN.

The Baker-lab ProteinMPNN entry point is ``protein_mpnn_run.py``; we look
for it via ``$PROTEINMPNN_PATH`` or PATH. Mock fallback emits deterministic
random-looking sequences whose length matches the backbone CA count.

Direct invocation::

    python -m utils.proteinmpnn_wrapper --check
    python -m utils.proteinmpnn_wrapper --run --backbone /path/design_0.pdb \
        --num-sequences 4
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
TOOL_NAME = "utils.proteinmpnn_wrapper"

_ENV_BINARY = "PROTEINMPNN_PATH"
_BINARY_CANDIDATES = ("protein_mpnn_run.py", "proteinmpnn", "protein_mpnn")

# Standard 20-aa alphabet for mock sampling. Order roughly mirrors PDB usage.
_AA_ALPHABET = "ARNDCQEGHILKMFPSTWYV"


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class ProteinMPNNError(RuntimeError):
    """Base error class for the ProteinMPNN wrapper."""


class ProteinMPNNNotFoundError(ProteinMPNNError):
    """Raised when ProteinMPNN binary cannot be located AND mock fallback is
    disabled."""


class ProteinMPNNInvocationError(ProteinMPNNError):
    """Raised when the subprocess returned non-zero or output unparseable."""


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class ProteinMPNNResult:
    """Output of :meth:`ProteinMPNNWrapper.run_mpnn`.

    Attributes
    ----------
    sequences : list[dict]
        ``[{"seq": "MPNN-output-seq", "score": float}, ...]``. Score is the
        ProteinMPNN per-residue NLL (lower = better) for real runs and a
        deterministic mock score in mock mode.
    backbone_pdb : str
        Input backbone path.
    num_sequences : int
        Requested number of sequences.
    backbone_length : int
        Detected backbone length (CA count) — useful when the upstream
        contigs spec is opaque.
    mock : bool
        ``True`` for mock-mode results, ``False`` otherwise.
    elapsed_s : float
        Wall-clock seconds.
    schema : str
        Schema tag.
    binary : str | None
        Resolved binary path (None for mock).
    """
    sequences: List[Dict[str, Any]]
    backbone_pdb: str
    num_sequences: int
    backbone_length: int = 0
    mock: bool = False
    elapsed_s: float = 0.0
    schema: str = SCHEMA_VERSION
    binary: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# Wrapper class
# ---------------------------------------------------------------------------

class ProteinMPNNWrapper:
    """Subprocess wrapper around ``protein_mpnn_run.py``."""

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
            cp = subprocess.run([path, "--help"], capture_output=True,
                                text=True, timeout=10)
            if "ProteinMPNN" in (cp.stdout or "") + (cp.stderr or ""):
                version = "ProteinMPNN-detected"
        except Exception:  # noqa: BLE001
            pass
        return {
            "available": True, "version": version, "path": path,
            "schema": SCHEMA_VERSION, "tool": TOOL_NAME,
        }

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run_mpnn(self, backbone_pdb: str, *,
                 num_sequences: int = 8,
                 sampling_temperature: float = 0.1,
                 chains_to_design: Optional[str] = None,
                 fixed_positions: Optional[List[int]] = None,
                 output_dir: Optional[str] = None,
                 ) -> ProteinMPNNResult:
        """Generate ``num_sequences`` sequences for ``backbone_pdb``.

        Parameters
        ----------
        backbone_pdb : str
            Path to a backbone PDB (RFdiffusion output).
        num_sequences : int
            Sequences to sample (default: 8).
        sampling_temperature : float
            ProteinMPNN sampling temperature.
        chains_to_design : str | None
            Comma-separated chain IDs (default: all-non-target).
        fixed_positions : list[int] | None
            Positions to keep fixed (1-indexed).
        output_dir : str | None
            Where to write intermediate files.

        Returns
        -------
        :class:`ProteinMPNNResult`

        Raises
        ------
        ProteinMPNNNotFoundError
            Binary missing and mock disabled.
        ProteinMPNNInvocationError
            Subprocess failure.
        ValueError
            ``num_sequences <= 0`` or empty backbone path.
        """
        if num_sequences <= 0:
            raise ValueError(f"num_sequences must be >= 1; got {num_sequences}")
        if not backbone_pdb:
            raise ValueError("backbone_pdb is required")
        if not Path(backbone_pdb).exists():
            raise ValueError(f"backbone PDB not found: {backbone_pdb}")

        path = self._resolve_binary()
        if not path:
            if not self.allow_mock:
                raise ProteinMPNNNotFoundError(
                    "ProteinMPNN binary not found in $PROTEINMPNN_PATH or PATH."
                )
            return self._mock_run(backbone_pdb=backbone_pdb,
                                  num_sequences=num_sequences)

        return self._real_run(
            binary=path,
            backbone_pdb=backbone_pdb,
            num_sequences=num_sequences,
            sampling_temperature=sampling_temperature,
            chains_to_design=chains_to_design,
            fixed_positions=fixed_positions,
            output_dir=output_dir,
        )

    # ------------------------------------------------------------------
    # Real subprocess invocation
    # ------------------------------------------------------------------

    def _real_run(self, *, binary: str, backbone_pdb: str,
                  num_sequences: int,
                  sampling_temperature: float,
                  chains_to_design: Optional[str],
                  fixed_positions: Optional[List[int]],
                  output_dir: Optional[str]) -> ProteinMPNNResult:
        out = Path(output_dir) if output_dir else Path(
            tempfile.mkdtemp(prefix="mpnn_"))
        out.mkdir(parents=True, exist_ok=True)

        argv: List[str] = [
            binary,
            "--pdb_path", str(backbone_pdb),
            "--out_folder", str(out),
            "--num_seq_per_target", str(int(num_sequences)),
            "--sampling_temp", str(float(sampling_temperature)),
        ]
        if chains_to_design:
            argv += ["--chains_to_design", chains_to_design]
        if fixed_positions:
            pos = " ".join(str(p) for p in fixed_positions)
            argv += ["--fixed_positions", pos]

        t0 = datetime.now(timezone.utc)
        try:
            cp = subprocess.run(argv, capture_output=True, text=True,
                                timeout=self.timeout_s, check=False)
        except subprocess.TimeoutExpired as exc:
            raise ProteinMPNNInvocationError(
                f"ProteinMPNN timed out after {self.timeout_s}s"
            ) from exc
        elapsed = (datetime.now(timezone.utc) - t0).total_seconds()

        if cp.returncode != 0:
            tail = (cp.stderr or "")[-2048:]
            raise ProteinMPNNInvocationError(
                f"ProteinMPNN exited {cp.returncode}; stderr tail:\n{tail}"
            )

        # Parse FASTA output ProteinMPNN typically emits to
        # ``out/seqs/<backbone_basename>.fa``.
        fa_path = self._locate_output_fasta(out, backbone_pdb)
        sequences = _parse_proteinmpnn_fasta(fa_path) if fa_path else []
        if not sequences:
            raise ProteinMPNNInvocationError(
                f"ProteinMPNN completed but no sequences parsed from {fa_path}"
            )
        return ProteinMPNNResult(
            sequences=sequences[:num_sequences],
            backbone_pdb=str(backbone_pdb),
            num_sequences=num_sequences,
            backbone_length=_count_ca_atoms(backbone_pdb),
            mock=False,
            elapsed_s=elapsed,
            binary=binary,
        )

    @staticmethod
    def _locate_output_fasta(out: Path,
                             backbone_pdb: str) -> Optional[Path]:
        stem = Path(backbone_pdb).stem
        for cand in (out / "seqs" / f"{stem}.fa",
                     out / f"{stem}.fa",
                     out / "seqs" / f"{stem}.fasta"):
            if cand.exists():
                return cand
        # Fallback: first *.fa under out/.
        hits = sorted(out.rglob("*.fa")) + sorted(out.rglob("*.fasta"))
        return hits[0] if hits else None

    # ------------------------------------------------------------------
    # Mock fallback
    # ------------------------------------------------------------------

    def _mock_run(self, *, backbone_pdb: str,
                  num_sequences: int) -> ProteinMPNNResult:
        """Emit ``num_sequences`` deterministic mock sequences whose length
        matches the backbone CA count (or 13, the demo default, when CA
        cannot be parsed).
        """
        L = _count_ca_atoms(backbone_pdb) or 13
        sequences: List[Dict[str, Any]] = []
        # Use SHA-256 of (backbone path + idx) as the per-sequence PRNG seed
        # so the mock is reproducible across runs.
        for i in range(num_sequences):
            seed = hashlib.sha256(
                f"{backbone_pdb}|{i}".encode("utf-8")).digest()
            seq_chars = []
            for j in range(L):
                seq_chars.append(_AA_ALPHABET[seed[j % len(seed)] %
                                              len(_AA_ALPHABET)])
            seq_str = "".join(seq_chars)
            score = round(0.30 + (i * 0.05) + (seed[0] % 50) / 1000.0, 4)
            sequences.append({"seq": seq_str, "score": score, "mock": True})
        return ProteinMPNNResult(
            sequences=sequences,
            backbone_pdb=str(backbone_pdb),
            num_sequences=num_sequences,
            backbone_length=L,
            mock=True,
            elapsed_s=0.0,
            binary=None,
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _count_ca_atoms(pdb_path: str) -> int:
    """Count CA atoms in a PDB file. Returns 0 if file cannot be read."""
    try:
        n = 0
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("ATOM") and " CA " in line[:20]:
                    n += 1
        return n
    except Exception:  # noqa: BLE001
        return 0


def _parse_proteinmpnn_fasta(fa_path: Path) -> List[Dict[str, Any]]:
    """Parse a ProteinMPNN FASTA output. Each record header looks like::

        >T=0.1, score=0.8732, ...

    We keep score from the header when present; fall back to ``None``.
    """
    if not fa_path or not fa_path.exists():
        return []
    out: List[Dict[str, Any]] = []
    cur_header: Optional[str] = None
    cur_seq: List[str] = []
    with open(fa_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_header is not None:
                    out.append(_assemble_record(cur_header, "".join(cur_seq)))
                cur_header = line[1:]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    if cur_header is not None:
        out.append(_assemble_record(cur_header, "".join(cur_seq)))
    # Skip the canonical first record (input native sequence) when present.
    if out and "score=" not in out[0].get("header", ""):
        out = out[1:]
    return out


def _assemble_record(header: str, seq: str) -> Dict[str, Any]:
    score: Optional[float] = None
    for tok in header.split(","):
        tok = tok.strip()
        if tok.startswith("score="):
            try:
                score = float(tok.split("=", 1)[1])
            except ValueError:
                score = None
    return {"seq": seq, "score": score, "header": header}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="utils.proteinmpnn_wrapper",
        description=("ProteinMPNN subprocess wrapper "
                     "(Stage 1.5 of UPDD de-novo discovery)."),
    )
    p.add_argument("--check", action="store_true")
    p.add_argument("--run", action="store_true")
    p.add_argument("--backbone", default=None,
                   help="Input backbone PDB.")
    p.add_argument("--num-sequences", type=int, default=8)
    p.add_argument("--sampling-temperature", type=float, default=0.1)
    p.add_argument("--chains-to-design", default=None)
    p.add_argument("--no-mock", action="store_true")
    p.add_argument("--binary", default=None)
    p.add_argument("--output-dir", default=None)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    w = ProteinMPNNWrapper(binary=args.binary,
                           allow_mock=(not args.no_mock))
    if args.check:
        info = w.check_availability()
        print(json.dumps(info, indent=2))
        return 0 if info["available"] else 1
    if args.run:
        if not args.backbone:
            print("--run requires --backbone", file=sys.stderr)
            return 2
        try:
            res = w.run_mpnn(
                backbone_pdb=args.backbone,
                num_sequences=args.num_sequences,
                sampling_temperature=args.sampling_temperature,
                chains_to_design=args.chains_to_design,
                output_dir=args.output_dir,
            )
        except ProteinMPNNError as exc:
            print(f"[proteinmpnn] {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            return 3
        print(json.dumps(res.to_dict(), indent=2))
        return 0
    print("[proteinmpnn] no action; use --check or --run", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
