"""utils/rfdiffusion_wrapper.py — RFdiffusion subprocess wrapper.

Schema: ``stage12_pipeline/0.1`` (PR-NEW-D Day 3, 2026-04-29).

Thin Python wrapper around the RFdiffusion binary (``run_inference.py`` from
the Baker-lab ``RFdiffusion`` repo). The wrapper:

    * Detects the binary via ``RFDIFFUSION_PATH`` env var or ``which``.
    * Invokes inference via :mod:`subprocess` with a contig string.
    * Falls back to a deterministic mock when the binary is absent so the
      Day-3 :mod:`scripts.updd_cli` ``--de-novo`` chain remains testable on
      CPU-only / no-tools-installed checkouts.

Result type :class:`RFdiffusionResult` is a dataclass with explicit ``mock:
bool`` flag — downstream consumers (Stage 1.5 ProteinMPNN feeder) must check
this flag rather than inferring real-vs-mock from the file content.

Direct invocation::

    python -m utils.rfdiffusion_wrapper --check
    python -m utils.rfdiffusion_wrapper --run --target /path/target.pdb \
        --contigs '[A1-100/0 50-100]' --num-designs 4
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


SCHEMA_VERSION = "stage12_pipeline/0.1"
TOOL_NAME = "utils.rfdiffusion_wrapper"

# Env var override for RFdiffusion binary location.
_ENV_BINARY = "RFDIFFUSION_PATH"
# Default candidate names that ``shutil.which`` will probe.
_BINARY_CANDIDATES = ("rfdiffusion", "run_inference.py")


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class RFdiffusionError(RuntimeError):
    """Base error class for the RFdiffusion wrapper."""


class RFdiffusionNotFoundError(RFdiffusionError):
    """Raised when RFdiffusion binary cannot be located AND mock fallback is
    disabled (``allow_mock=False``)."""


class RFdiffusionInvocationError(RFdiffusionError):
    """Raised when the subprocess returned a non-zero exit code."""


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class RFdiffusionResult:
    """Output of :meth:`RFdiffusionWrapper.run_diffusion`.

    Attributes
    ----------
    backbones : list[str]
        Absolute paths to backbone PDB files (one per design).
    target_pdb : str
        Path to the input target PDB the diffusion was conditioned on.
    contigs : str
        The contig spec passed to RFdiffusion.
    num_designs : int
        Requested design count (== ``len(backbones)`` on success).
    hotspot_residues : list[str] | None
        Optional hotspot list (e.g. ``["A12", "A45"]``) propagated through.
    mock : bool
        ``True`` when the result was synthesised because no RFdiffusion
        binary was available. Real-tool runs set this to ``False``.
    elapsed_s : float
        Wall-clock seconds for the diffusion call.
    schema : str
        Schema tag (``stage12_pipeline/0.1``).
    binary : str | None
        Resolved binary path used for the run (``None`` for mock).
    stderr_tail : str
        Last 2 KB of stderr (real runs only); empty for mock.
    """
    backbones: List[str]
    target_pdb: str
    contigs: str
    num_designs: int
    hotspot_residues: Optional[List[str]] = None
    mock: bool = False
    elapsed_s: float = 0.0
    schema: str = SCHEMA_VERSION
    binary: Optional[str] = None
    stderr_tail: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# Wrapper class
# ---------------------------------------------------------------------------

class RFdiffusionWrapper:
    """Subprocess wrapper around the RFdiffusion ``run_inference.py`` binary.

    Parameters
    ----------
    binary : str | None
        Explicit path to the RFdiffusion binary (overrides env / PATH lookup).
    allow_mock : bool
        When ``True`` (default), :meth:`run_diffusion` falls back to
        synthesising mock backbone PDBs if no binary is found. When
        ``False``, missing binary raises :class:`RFdiffusionNotFoundError`.
    timeout_s : int | None
        Optional subprocess timeout in seconds. ``None`` (default) ⇒ no
        timeout. Set this for CI / smoke-test runs.
    """

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
        """Return a binary path or ``None``. Lookup order:

            1. ``self._explicit_binary`` (constructor argument)
            2. ``$RFDIFFUSION_PATH`` env var
            3. ``shutil.which`` over :data:`_BINARY_CANDIDATES`
        """
        if self._explicit_binary:
            if Path(self._explicit_binary).exists():
                return str(self._explicit_binary)
            return None
        env = os.environ.get(_ENV_BINARY)
        if env and Path(env).exists():
            return env
        for cand in _BINARY_CANDIDATES:
            p = shutil.which(cand)
            if p:
                return p
        return None

    def check_availability(self) -> Dict[str, Any]:
        """Return a dict ``{available, version, path}`` describing the local
        RFdiffusion install.

        ``version`` is reported as ``"unknown"`` when the binary exists but
        ``--version`` is not supported by the upstream script. The function
        never raises — failures are reflected via ``available=False``.
        """
        path = self._resolve_binary()
        if not path:
            return {
                "available": False,
                "version": None,
                "path": None,
                "schema": SCHEMA_VERSION,
                "tool": TOOL_NAME,
            }
        version = "unknown"
        try:
            cp = subprocess.run(
                [path, "--version"],
                capture_output=True, text=True, timeout=10,
            )
            out = (cp.stdout or "") + (cp.stderr or "")
            if cp.returncode == 0 and out.strip():
                version = out.strip().splitlines()[0][:80]
        except Exception:  # noqa: BLE001
            # --version is not part of every upstream distribution; treat
            # absence as "unknown" rather than "unavailable".
            version = "unknown"
        return {
            "available": True,
            "version": version,
            "path": path,
            "schema": SCHEMA_VERSION,
            "tool": TOOL_NAME,
        }

    # ------------------------------------------------------------------
    # Diffusion
    # ------------------------------------------------------------------

    def run_diffusion(self, target_pdb: str, contigs: str,
                      num_designs: int, *,
                      hotspot_residues: Optional[List[str]] = None,
                      output_dir: Optional[str] = None,
                      extra_args: Optional[List[str]] = None,
                      ) -> RFdiffusionResult:
        """Generate ``num_designs`` backbone PDBs conditioned on ``target_pdb``.

        Parameters
        ----------
        target_pdb : str
            Path to the receptor / target PDB.
        contigs : str
            RFdiffusion contig spec (e.g. ``"[A1-100/0 50-100]"``).
        num_designs : int
            Number of backbones to sample. Must be ``>= 1``.
        hotspot_residues : list[str] | None
            Optional hotspot residues (e.g. ``["A12", "A45"]``).
        output_dir : str | None
            Directory for the design PDBs. ``None`` (default) ⇒ a
            tmp directory is created.
        extra_args : list[str] | None
            Additional CLI args appended verbatim.

        Returns
        -------
        :class:`RFdiffusionResult`
            ``mock=True`` when binary missing AND ``self.allow_mock`` is set.

        Raises
        ------
        RFdiffusionNotFoundError
            Binary missing and ``allow_mock`` is ``False``.
        RFdiffusionInvocationError
            Real subprocess call returned non-zero.
        ValueError
            On clearly-invalid arguments (e.g. ``num_designs <= 0``).
        """
        if num_designs <= 0:
            raise ValueError(f"num_designs must be >= 1; got {num_designs}")
        if not target_pdb:
            raise ValueError("target_pdb is required")

        path = self._resolve_binary()
        if not path:
            if not self.allow_mock:
                raise RFdiffusionNotFoundError(
                    "RFdiffusion binary not found in $RFDIFFUSION_PATH or "
                    "PATH. Set RFDIFFUSION_PATH or pass binary=... explicitly."
                )
            return self._mock_run(
                target_pdb=target_pdb, contigs=contigs,
                num_designs=num_designs,
                hotspot_residues=hotspot_residues,
                output_dir=output_dir,
            )

        # Real path — invoke subprocess.
        return self._real_run(
            binary=path, target_pdb=target_pdb, contigs=contigs,
            num_designs=num_designs,
            hotspot_residues=hotspot_residues,
            output_dir=output_dir, extra_args=extra_args,
        )

    # ------------------------------------------------------------------
    # Real subprocess invocation
    # ------------------------------------------------------------------

    def _real_run(self, *, binary: str, target_pdb: str, contigs: str,
                  num_designs: int,
                  hotspot_residues: Optional[List[str]],
                  output_dir: Optional[str],
                  extra_args: Optional[List[str]]) -> RFdiffusionResult:
        out = Path(output_dir) if output_dir else Path(
            tempfile.mkdtemp(prefix="rfdiff_"))
        out.mkdir(parents=True, exist_ok=True)

        argv: List[str] = [binary,
                           f"inference.input_pdb={target_pdb}",
                           f"contigmap.contigs={contigs}",
                           f"inference.num_designs={int(num_designs)}",
                           f"inference.output_prefix={out}/design"]
        if hotspot_residues:
            hot = ",".join(hotspot_residues)
            argv.append(f"ppi.hotspot_res=[{hot}]")
        if extra_args:
            argv.extend(extra_args)

        t0 = datetime.now(timezone.utc)
        try:
            cp = subprocess.run(
                argv, capture_output=True, text=True,
                timeout=self.timeout_s, check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise RFdiffusionInvocationError(
                f"RFdiffusion timed out after {self.timeout_s}s"
            ) from exc
        elapsed = (datetime.now(timezone.utc) - t0).total_seconds()

        if cp.returncode != 0:
            tail = (cp.stderr or "")[-2048:]
            raise RFdiffusionInvocationError(
                f"RFdiffusion exited {cp.returncode}; stderr tail:\n{tail}"
            )

        # Collect produced backbone PDBs (design_*.pdb / design_0.pdb …).
        backbones = sorted(str(p) for p in out.glob("design*.pdb"))
        if not backbones:
            # Some upstream variants emit *.pdb in nested dirs; widen.
            backbones = sorted(str(p) for p in out.rglob("*.pdb"))
        if not backbones:
            raise RFdiffusionInvocationError(
                f"RFdiffusion completed but no PDBs found under {out}"
            )

        return RFdiffusionResult(
            backbones=backbones,
            target_pdb=str(target_pdb),
            contigs=contigs,
            num_designs=num_designs,
            hotspot_residues=hotspot_residues,
            mock=False,
            elapsed_s=elapsed,
            binary=binary,
            stderr_tail=(cp.stderr or "")[-512:],
        )

    # ------------------------------------------------------------------
    # Mock fallback
    # ------------------------------------------------------------------

    def _mock_run(self, *, target_pdb: str, contigs: str,
                  num_designs: int,
                  hotspot_residues: Optional[List[str]],
                  output_dir: Optional[str]) -> RFdiffusionResult:
        """Synthesise ``num_designs`` minimal backbone PDB stubs.

        Each stub is a 3-CA-atom polyglycine placeholder. This is **not** a
        physically meaningful backbone — it exists purely to drive Stage
        1.5 / 2 mock chain in the de-novo discovery test harness.
        """
        out = Path(output_dir) if output_dir else Path(
            tempfile.mkdtemp(prefix="rfdiff_mock_"))
        out.mkdir(parents=True, exist_ok=True)

        backbones: List[str] = []
        for i in range(num_designs):
            p = out / f"mock_design_{i}.pdb"
            p.write_text(_synthesise_minimal_backbone(i, contigs))
            backbones.append(str(p))

        return RFdiffusionResult(
            backbones=backbones,
            target_pdb=str(target_pdb),
            contigs=contigs,
            num_designs=num_designs,
            hotspot_residues=hotspot_residues,
            mock=True,
            elapsed_s=0.0,
            binary=None,
            stderr_tail="",
        )


def _synthesise_minimal_backbone(idx: int, contigs: str) -> str:
    """Return a 3-CA polyglycine PDB stub. Deterministic per idx."""
    # Three CA atoms spaced 3.8 Å — placeholder geometry only.
    z = float(idx) * 1.0
    return (
        f"REMARK  Mock RFdiffusion backbone idx={idx} contigs={contigs}\n"
        f"ATOM      1  CA  GLY A   1       0.000   0.000  {z:7.3f}  1.00 30.00           C\n"
        f"ATOM      2  CA  GLY A   2       3.800   0.000  {z:7.3f}  1.00 30.00           C\n"
        f"ATOM      3  CA  GLY A   3       7.600   0.000  {z:7.3f}  1.00 30.00           C\n"
        f"TER       4      GLY A   3\n"
        f"END\n"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="utils.rfdiffusion_wrapper",
        description=("RFdiffusion subprocess wrapper "
                     "(Stage 1 of UPDD de-novo discovery)."),
    )
    p.add_argument("--check", action="store_true",
                   help="Probe binary availability and exit.")
    p.add_argument("--run", action="store_true",
                   help="Run a diffusion job (requires --target/--contigs).")
    p.add_argument("--target", default=None, help="Target PDB path.")
    p.add_argument("--contigs", default="[A1-100/0 50-100]",
                   help="Contig spec (default: '[A1-100/0 50-100]').")
    p.add_argument("--num-designs", type=int, default=1)
    p.add_argument("--hotspot-residues", default=None,
                   help="Comma-separated, e.g. 'A12,A45'.")
    p.add_argument("--output-dir", default=None)
    p.add_argument("--no-mock", action="store_true",
                   help="Disable mock fallback when binary missing.")
    p.add_argument("--binary", default=None,
                   help="Explicit RFdiffusion binary path.")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    w = RFdiffusionWrapper(binary=args.binary,
                           allow_mock=(not args.no_mock))
    if args.check:
        info = w.check_availability()
        print(json.dumps(info, indent=2))
        return 0 if info["available"] else 1
    if args.run:
        if not args.target:
            print("--run requires --target", file=sys.stderr)
            return 2
        hot = ([h.strip() for h in args.hotspot_residues.split(",") if h.strip()]
               if args.hotspot_residues else None)
        try:
            res = w.run_diffusion(
                target_pdb=args.target, contigs=args.contigs,
                num_designs=args.num_designs,
                hotspot_residues=hot, output_dir=args.output_dir,
            )
        except RFdiffusionError as exc:
            print(f"[rfdiffusion] {type(exc).__name__}: {exc}", file=sys.stderr)
            return 3
        print(json.dumps(res.to_dict(), indent=2))
        return 0
    print("[rfdiffusion] no action; use --check or --run", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
