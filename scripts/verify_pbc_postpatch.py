#!/usr/bin/env python
"""scripts/verify_pbc_postpatch.py — verify saved snapshot PDBs.

Reads PDBs from ``snapshots_n25_postl387_patch/`` (or ``snapshots_n25/``
if --pre-patch given) and reports the prev_C ↔ ncAA_N peptide-bond
distance distribution per system. Unlike ``verify_pbc_integrity.py``
(which scans the raw DCD), this verifier validates the *output of
save_snapshots* — exactly the artefact that downstream MM-GBSA / DFT
QM/MM consumes.

Each saved PDB is its own single-frame snapshot, so the median over the
25 PDBs is the proper post-patch integrity statistic.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

warnings.filterwarnings("ignore")
logging.getLogger("mdtraj").setLevel(logging.ERROR)

import mdtraj as md  # noqa: E402
import numpy as np  # noqa: E402

try:
    md.set_logger_level(logging.ERROR)  # type: ignore[attr-defined]
except Exception:
    pass


STANDARD_AA = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX",
    "ACE", "NME", "NMA",  # caps
})

SOLVENT_ION = frozenset({
    "HOH", "WAT", "TIP", "TIP3", "TIP3P", "TIP4P",
    "NA", "NA+", "CL", "CL-", "K", "K+", "MG", "MG2", "CA", "CA2",
    "ZN", "FE", "MN", "CU",
})

INTACT_MAX = 5.0
BROKEN_MIN = 50.0


@dataclass
class SystemRow:
    tag: str
    n_pdb: int = 0
    median_ang: Optional[float] = None
    min_ang: Optional[float] = None
    max_ang: Optional[float] = None
    status: str = "SKIPPED"
    skip_reason: Optional[str] = None
    per_frame: List[float] = field(default_factory=list)

    def render(self) -> str:
        tag = self.tag.ljust(36)[:36]
        if self.status == "SKIPPED":
            return f"{tag} | n=- | med=-      | SKIPPED ({self.skip_reason or '-'})"
        med = f"{self.median_ang:>6.2f}" if self.median_ang is not None else "  -- "
        mn = f"{self.min_ang:>6.2f}" if self.min_ang is not None else "  -- "
        mx = f"{self.max_ang:>6.2f}" if self.max_ang is not None else "  -- "
        return (
            f"{tag} | n={self.n_pdb:>2d} | med={med} | "
            f"min={mn} max={mx} | {self.status}"
        )


def _is_binder_chain(chain) -> bool:
    residues = list(chain.residues)
    n = len(residues)
    if n < 2 or n > 50:
        return False
    if any(r.name.upper() in SOLVENT_ION for r in residues):
        return False
    return True


def _find_binder_chain(top):
    for ch in top.chains:
        if _is_binder_chain(ch):
            return ch
    return None


def _find_ncaa_in_chain(chain):
    for r in chain.residues:
        nm = r.name.upper()
        if nm in STANDARD_AA or nm in SOLVENT_ION:
            continue
        return r
    return None


def _prev_residue(chain, ncaa_residue):
    prev = None
    for r in chain.residues:
        if r.index == ncaa_residue.index:
            return prev
        prev = r
    return None


def _atom_by_name(residue, atom_name: str):
    for a in residue.atoms:
        if a.name == atom_name:
            return a
    return None


def _measure_pdb(pdb_path: Path) -> Optional[float]:
    try:
        frame = md.load(str(pdb_path))
    except Exception:
        return None
    binder = _find_binder_chain(frame.topology)
    if binder is None:
        return None
    ncaa = _find_ncaa_in_chain(binder)
    if ncaa is None:
        return None
    prev = _prev_residue(binder, ncaa)
    if prev is None:
        return None
    pc = _atom_by_name(prev, "C")
    nn = _atom_by_name(ncaa, "N")
    if pc is None or nn is None:
        return None
    try:
        d_nm = float(np.linalg.norm(
            frame.xyz[0, pc.index, :] - frame.xyz[0, nn.index, :]
        ))
    except Exception:
        return None
    return d_nm * 10.0


def measure_system_dir(snap_dir: Path) -> SystemRow:
    row = SystemRow(tag=snap_dir.parent.name)
    if not snap_dir.is_dir():
        row.skip_reason = f"missing {snap_dir.name}"
        return row
    pdbs = sorted(snap_dir.glob("*_snap*_f*.pdb"))
    if not pdbs:
        row.skip_reason = f"no PDBs in {snap_dir.name}"
        return row
    distances: List[float] = []
    for p in pdbs:
        d = _measure_pdb(p)
        if d is not None:
            distances.append(d)
    row.n_pdb = len(distances)
    if not distances:
        row.skip_reason = "no measurable peptide bond"
        return row
    arr = np.array(distances)
    row.median_ang = float(np.median(arr))
    row.min_ang = float(arr.min())
    row.max_ang = float(arr.max())
    row.per_frame = [round(x, 3) for x in distances]
    if row.median_ang < INTACT_MAX:
        row.status = "INTACT"
    elif row.median_ang > BROKEN_MIN:
        row.status = "BROKEN"
    else:
        row.status = "AMBIGUOUS"
    return row


def render_report(rows: List[SystemRow]) -> str:
    lines = [
        "Post-patch PBC Integrity (saved PDB snapshots)",
        "=" * 72,
        "Tag                                  | n  | median Å | min/max     | Status",
        "-" * 72,
    ]
    for r in rows:
        lines.append(r.render())
    intact = sum(1 for r in rows if r.status == "INTACT")
    broken = sum(1 for r in rows if r.status == "BROKEN")
    amb = sum(1 for r in rows if r.status == "AMBIGUOUS")
    skip = sum(1 for r in rows if r.status == "SKIPPED")
    total = len(rows)
    checked = intact + broken + amb
    lines.append("")
    lines.append("Summary:")
    lines.append(f"  total={total}  INTACT={intact}  BROKEN={broken}  AMBIGUOUS={amb}  SKIPPED={skip}")
    if checked:
        lines.append(f"  intact_pct={100.0 * intact / checked:.1f}%")
    lines.append("")
    if broken == 0 and amb == 0 and skip == 0:
        lines.append("Verdict: ALL_INTACT")
    elif broken == 0:
        lines.append("Verdict: NO_BROKEN (some AMBIGUOUS / SKIPPED)")
    else:
        lines.append("Verdict: DEFECTS_REMAIN")
    return "\n".join(lines) + "\n"


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--targets-json",
        default=None,
        help=("Optional JSON listing systems to check. Default: scan "
              "outputs/*_calib_*/snapshots_n25_postl387_patch."),
    )
    p.add_argument(
        "--root", default="/home/san/UPDD_proj/outputs",
    )
    p.add_argument(
        "--prefer-postl387", action="store_true",
        help="If patched dir present, use it; otherwise fallback to snapshots_n25.",
    )
    p.add_argument("--output", default=None)
    p.add_argument("--json-out", default=None)
    p.add_argument(
        "--filter", default=None,
        help="Substring filter on tag.",
    )
    args = p.parse_args(argv)

    root = Path(args.root)
    rows: List[SystemRow] = []

    if args.targets_json:
        data = json.loads(Path(args.targets_json).read_text())
        records = data.get("records", [])
        tags = [r["tag"] for r in records]
    else:
        tags = sorted(d.name for d in root.glob("*_calib_*") if d.is_dir())

    for tag in tags:
        if args.filter and args.filter not in tag:
            continue
        sysdir = root / tag
        patched = sysdir / "snapshots_n25_postl387_patch"
        legacy = sysdir / "snapshots_n25"
        if args.prefer_postl387 and patched.is_dir():
            chosen = patched
        elif patched.is_dir():
            chosen = patched
        else:
            chosen = legacy
        row = measure_system_dir(chosen)
        # Override row.tag to reflect the actual scanned dir
        row.tag = f"{tag}/{chosen.name}"
        rows.append(row)

    report = render_report(rows)
    if args.output:
        Path(args.output).write_text(report)
    sys.stdout.write(report)

    if args.json_out:
        Path(args.json_out).write_text(json.dumps({
            "schema": "pbc_postpatch/0.1",
            "rows": [
                {
                    "tag": r.tag,
                    "n_pdb": r.n_pdb,
                    "median_ang": r.median_ang,
                    "min_ang": r.min_ang,
                    "max_ang": r.max_ang,
                    "status": r.status,
                    "skip_reason": r.skip_reason,
                    "per_frame": r.per_frame,
                }
                for r in rows
            ],
        }, indent=2))

    broken = sum(1 for r in rows if r.status == "BROKEN")
    return 0 if broken == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
