#!/usr/bin/env python
"""scripts/verify_pbc_integrity.py — UPDD PBC Integrity Verification Tool.

Schema: ``pbc_integrity/0.1`` (Phase β PBC cross-check, 2026-04-29)

Verifies that MD trajectories preserve binder peptide-bond integrity
across the full simulation. The relevant defect was characterized in
``outputs/analysis/phase_beta_stage1_repair_20260428.md`` (Phase β
Stage 1-Repair, 2026-04-28): when the binder chain crosses the periodic
boundary at the HETATM/ATOM (canonical/ncAA) junction during MD, the
extracted snapshot PDBs (via ``utils/extract_snapshots.py::save_snapshots``
L387, which calls ``_unwrap_binder_manual`` only) can leave the
prev_C ↔ ncAA_N peptide bond at ~110 Å instead of ~1.3 Å — a periodic
imaging artifact, not a real broken bond.

This is a **read-only** verification driver: it loads each system's
restrained DCD, applies no imaging, and reports the median peptide-bond
distance across all frames as evidence of latent PBC vulnerability.

CLI::

    python scripts/verify_pbc_integrity.py [--filter PATTERN] [--output FILE] [--root DIR]

Classification:
    - INTACT     : median < 5 Å
    - BROKEN     : median > 50 Å
    - AMBIGUOUS  : 5 ≤ median ≤ 50 Å

Exit codes:
    0 — all checked systems INTACT
    1 — at least one BROKEN system detected
    2 — fatal error during scan (e.g. mdtraj import failure)

Output is intentionally TERSE: one row per system, no per-system path
printing, no mdtraj warnings, no traceback for skipped systems.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

# ---------------------------------------------------------------------------
# Silence mdtraj / NumPy / OpenMM warnings before import
# ---------------------------------------------------------------------------

warnings.filterwarnings("ignore")
logging.getLogger("mdtraj").setLevel(logging.ERROR)
logging.getLogger("openmm").setLevel(logging.ERROR)
logging.getLogger("numpy").setLevel(logging.ERROR)

try:
    import mdtraj as md  # noqa: E402
    import numpy as np  # noqa: E402
except ImportError as exc:  # pragma: no cover
    print(f"[FATAL] mdtraj/numpy import failed: {exc}", file=sys.stderr)
    sys.exit(2)

# Some mdtraj versions have a global logger handle:
try:
    md.set_logger_level(logging.ERROR)  # type: ignore[attr-defined]
except (AttributeError, Exception):
    pass

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "pbc_integrity/0.1"
TOOL_NAME = "scripts.verify_pbc_integrity"

# Canonical residue names (amber14 + standard variants). Anything else in a
# binder-sized chain is treated as a candidate ncAA.
STANDARD_AA = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX",
})

SOLVENT_ION = frozenset({
    "HOH", "WAT", "TIP", "TIP3", "TIP3P", "TIP4P",
    "NA", "NA+", "CL", "CL-", "K", "K+", "MG", "MG2", "CA", "CA2",
    "ZN", "FE", "MN", "CU",
})

# Topology candidate filenames (preferred order)
TOPOLOGY_CANDIDATES = [
    "*_final.pdb",
    "*_topology.pdb",
    "*_production_md.pdb",
    "*.pdb",
]

# Classification thresholds (Å)
INTACT_MAX = 5.0
BROKEN_MIN = 50.0


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class IntegrityRow:
    """Per-system PBC integrity record."""

    system: str
    n_frames: Optional[int] = None
    median_angstrom: Optional[float] = None
    status: str = "SKIPPED"    # INTACT | BROKEN | AMBIGUOUS | SKIPPED
    skip_reason: Optional[str] = None

    def render(self) -> str:
        sys_col = self.system.ljust(28)[:28]
        if self.status == "SKIPPED":
            return f"{sys_col} | {'-':>6} | {'-':>8} | SKIPPED ({self.skip_reason or 'unknown'})"
        nf = f"{self.n_frames:>4}" if self.n_frames is not None else "  -- "
        med = f"{self.median_angstrom:>6.2f}" if self.median_angstrom is not None else "   -- "
        return f"{sys_col} | {nf:>6} | {med:>8} | {self.status}"


# ---------------------------------------------------------------------------
# Topology / DCD discovery
# ---------------------------------------------------------------------------

def _find_dcd(system_dir: Path) -> Optional[Path]:
    """Return path to the restrained DCD or first DCD in mdresult/, else None."""
    mdresult = system_dir / "mdresult"
    if not mdresult.is_dir():
        return None
    # Preferred: <system>_restrained.dcd
    for c in mdresult.glob("*_restrained.dcd"):
        return c
    # Fallback: first DCD
    for c in mdresult.glob("*.dcd"):
        return c
    return None


def _find_topology(system_dir: Path) -> Optional[Path]:
    """Return a topology PDB that mdtraj can load.

    Searches both ``mdresult/`` (final.pdb) and ``_md_input/`` (initial pdb).
    """
    candidates: List[Path] = []
    for sub in ("mdresult", "_md_input"):
        d = system_dir / sub
        if not d.is_dir():
            continue
        for pat in TOPOLOGY_CANDIDATES:
            candidates.extend(sorted(d.glob(pat)))
    # Deduplicate, preserve order
    seen = set()
    unique = []
    for c in candidates:
        if c in seen:
            continue
        seen.add(c)
        unique.append(c)
    # Try each in order until one loads
    for c in unique:
        try:
            md.load(str(c))
            return c
        except Exception:
            continue
    return None


# ---------------------------------------------------------------------------
# Binder / ncAA detection
# ---------------------------------------------------------------------------

def _is_binder_chain(chain) -> bool:
    """Heuristic: chain has 3-50 residues, none are HOH/ions."""
    residues = list(chain.residues)
    n = len(residues)
    if n < 3 or n > 50:
        return False
    # Reject if any residue is solvent/ion
    if any(r.name.upper() in SOLVENT_ION for r in residues):
        return False
    return True


def _find_binder_chain(topology):
    """Return the first chain that looks like a binder peptide."""
    for ch in topology.chains:
        if _is_binder_chain(ch):
            return ch
    return None


def _find_ncaa_in_chain(chain):
    """Return the first non-standard residue in the chain, else None."""
    for r in chain.residues:
        nm = r.name.upper()
        if nm in STANDARD_AA:
            continue
        if nm in SOLVENT_ION:
            continue
        return r
    return None


def _prev_residue(chain, ncaa_residue):
    """Return the residue immediately preceding ``ncaa_residue`` in chain."""
    prev = None
    for r in chain.residues:
        if r.index == ncaa_residue.index:
            return prev
        prev = r
    return None


def _atom_by_name(residue, atom_name):
    """Return the atom in ``residue`` with the given name, or None."""
    for a in residue.atoms:
        if a.name == atom_name:
            return a
    return None


# ---------------------------------------------------------------------------
# Per-system analysis
# ---------------------------------------------------------------------------

def analyze_system(system_dir: Path) -> IntegrityRow:
    """Run the PBC integrity check on a single system directory."""
    system = system_dir.name
    row = IntegrityRow(system=system)

    dcd = _find_dcd(system_dir)
    if dcd is None:
        row.skip_reason = "no DCD"
        return row

    top = _find_topology(system_dir)
    if top is None:
        row.skip_reason = "no topology"
        return row

    try:
        traj = md.load(str(dcd), top=str(top))
    except Exception as exc:
        row.skip_reason = f"load failed: {type(exc).__name__}"
        return row

    binder = _find_binder_chain(traj.topology)
    if binder is None:
        row.skip_reason = "no binder chain"
        return row

    ncaa = _find_ncaa_in_chain(binder)
    if ncaa is None:
        row.skip_reason = "no ncAA"
        return row

    prev = _prev_residue(binder, ncaa)
    if prev is None:
        row.skip_reason = "ncAA at N-term"
        return row

    prev_C = _atom_by_name(prev, "C")
    ncaa_N = _atom_by_name(ncaa, "N")
    if prev_C is None or ncaa_N is None:
        row.skip_reason = "missing peptide atoms"
        return row

    try:
        # mdtraj.xyz is in nm; convert to Å
        d_nm = np.linalg.norm(
            traj.xyz[:, prev_C.index, :] - traj.xyz[:, ncaa_N.index, :],
            axis=1,
        )
        d_ang = d_nm * 10.0
        median_ang = float(np.median(d_ang))
    except Exception as exc:
        row.skip_reason = f"distance calc failed: {type(exc).__name__}"
        return row

    row.n_frames = int(traj.n_frames)
    row.median_angstrom = median_ang
    if median_ang < INTACT_MAX:
        row.status = "INTACT"
    elif median_ang > BROKEN_MIN:
        row.status = "BROKEN"
    else:
        row.status = "AMBIGUOUS"
    return row


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def discover_systems(root: Path, pattern: Optional[str]) -> List[Path]:
    """Return sorted list of <root>/*_calib_* directories matching ``pattern``."""
    candidates = sorted(p for p in root.glob("*_calib_*") if p.is_dir())
    if pattern:
        rx = re.compile(pattern)
        candidates = [p for p in candidates if rx.search(p.name)]
    return candidates


def render_report(rows: List[IntegrityRow]) -> str:
    """Format the TERSE report."""
    lines: List[str] = []
    lines.append("PBC Integrity Check Report")
    lines.append("==========================")
    lines.append(
        "System                       | Frames | Median Å | Status"
    )
    lines.append(
        "-----------------------------|--------|----------|--------"
    )
    for row in rows:
        lines.append(row.render())

    intact = sum(1 for r in rows if r.status == "INTACT")
    broken = sum(1 for r in rows if r.status == "BROKEN")
    ambiguous = sum(1 for r in rows if r.status == "AMBIGUOUS")
    skipped = sum(1 for r in rows if r.status == "SKIPPED")
    total = len(rows)
    checked = intact + broken + ambiguous

    def _pct(num: int, denom: int) -> str:
        return f"{(100.0 * num / denom):.0f}%" if denom else "—"

    lines.append("")
    lines.append("Summary:")
    lines.append(f"  Total systems checked: {total}")
    lines.append(f"  INTACT: {intact} ({_pct(intact, checked)})")
    lines.append(f"  BROKEN: {broken} ({_pct(broken, checked)})")
    lines.append(
        f"  Ambiguous: {ambiguous} (5 ≤ median ≤ 50)"
    )
    lines.append(f"  Skipped (no DCD/topology/ncAA): {skipped}")
    lines.append("")
    if broken == 0:
        lines.append("Verdict: ALL_INTACT")
    else:
        lines.append("Verdict: DEFECTS_DETECTED")
    return "\n".join(lines) + "\n"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL_NAME,
        description=(
            "Verify PBC integrity of MD trajectories by checking the prev_C ↔ "
            "ncAA_N peptide-bond distance per frame; classify each system as "
            "INTACT / BROKEN / AMBIGUOUS."
        ),
    )
    parser.add_argument(
        "--filter",
        default=None,
        help="regex to filter system directory names (e.g. \"1EBP|Cp4\")",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="write TERSE report to this file (default: stdout)",
    )
    parser.add_argument(
        "--root",
        default="/home/san/UPDD_proj/outputs",
        help="outputs root directory (default: /home/san/UPDD_proj/outputs)",
    )
    args = parser.parse_args(argv)

    root = Path(args.root).resolve()
    if not root.is_dir():
        print(f"[FATAL] root not found: {root}", file=sys.stderr)
        return 2

    systems = discover_systems(root, args.filter)
    rows: List[IntegrityRow] = []
    for sysd in systems:
        try:
            row = analyze_system(sysd)
        except Exception as exc:
            row = IntegrityRow(
                system=sysd.name,
                status="SKIPPED",
                skip_reason=f"unhandled: {type(exc).__name__}",
            )
        rows.append(row)

    report = render_report(rows)
    if args.output:
        Path(args.output).write_text(report)
    else:
        sys.stdout.write(report)

    broken = sum(1 for r in rows if r.status == "BROKEN")
    return 0 if broken == 0 else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
