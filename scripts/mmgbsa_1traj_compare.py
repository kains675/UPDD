#!/usr/bin/env python
"""
mmgbsa_1traj_compare.py
-----------------------
1-snapshot 1-trajectory MM-GBSA (Fix 4, Genheden & Ryde 2015 §2.1,
DOI: 10.1517/17460441.2015.1032250) comparison helper.

- Builds ONE complex System (amber14 + GBn2 implicit solvent) and minimizes it.
- Evaluates E_complex, E_receptor, E_ligand from the SAME minimized geometry
  (subsystem-rebuild strategy — no independent minimization of receptor/ligand).
- Compares against the existing 3-trajectory ΔG in the matching
  `outputs/{run}/mmgbsa_results/mmgbsa_summary.json`.

WT only (v1): ncaa_label != "none" raises NotImplementedError because ncAA
parameter wiring (params_manifest + hydrogen defs + cofactor mol2) is tied to
the full MD run directory, not a stand-alone snapshot. Cp4/ncAA variants are
explicitly out of scope for Phase 3 test.

Constraints honored:
    - F1: Mg²⁺ is NOT stripped (WT 2QKI/3IOL snapshots have no Mg²⁺ anyway).
    - F2: No +5 Å translate / fabricated-OXT hack; we never move coordinates.
    - P1: 5000-iter minimum via --minimize_iter default.
    - R-11: ranking-only output, no absolute-ΔG claim.
"""

import argparse
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from openmm import app, unit
import openmm as mm
from openmm.app import ForceField, HBonds, PDBFile, Simulation

# Reuse the canonical FF bundle + strip rules from utils/run_mmgbsa.py so the
# comparison is apples-to-apples (same templates, same Born radii overrides).
_UTILS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "utils"))
if _UTILS_DIR not in sys.path:
    sys.path.insert(0, _UTILS_DIR)
from run_mmgbsa import (  # noqa: E402
    _CAP_RESNAMES,
    _MMGBSA_STRIP_RESNAMES,
    apply_gb_radius_override,
)
from utils_common import KNOWN_COFACTORS  # noqa: E402


# Matches outputs/{run}/snapshots/{snap_stem}.pdb (allows nested runs).
_SNAP_PATH_RE = re.compile(
    r".*/outputs/(?P<run>[^/]+)/snapshots/(?P<stem>[^/]+)\.pdb$"
)


# ==========================================
# PDB parsing — classify atoms by bucket
# ==========================================
def classify_atoms(pdb_path: str, receptor_chain: str, binder_chain: str
                   ) -> Tuple[List[str], List[str], List[str]]:
    """Split raw ATOM/HETATM/TER/END lines into complex / receptor / binder
    buckets using the same rules as utils/run_mmgbsa.split_complex, but without
    writing tmp PDBs (we keep everything in-memory for the single System).

    Returns: (complex_lines, receptor_lines, binder_lines) — all with TER/END
    kept on the complex bucket only, and CONECT lines filtered per bucket.
    """
    complex_lines: List[str] = []
    receptor_lines: List[str] = []
    binder_lines: List[str] = []
    serial_bucket: Dict[int, set] = {}
    conect_raw: List[str] = []

    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "CONECT":
                conect_raw.append(line)
                continue
            if rec not in ("ATOM", "HETATM", "TER", "END"):
                continue
            if rec in ("TER", "END"):
                complex_lines.append(line)
                continue

            chain = line[21] if len(line) > 21 else " "
            resname = line[17:20].strip() if len(line) >= 20 else ""

            if resname in _MMGBSA_STRIP_RESNAMES or resname in _CAP_RESNAMES:
                continue

            complex_lines.append(line)
            try:
                serial = int(line[6:11])
            except (ValueError, IndexError):
                serial = None
            if serial is not None:
                serial_bucket.setdefault(serial, set()).add("complex")

            if chain == receptor_chain:
                receptor_lines.append(line)
                if serial is not None:
                    serial_bucket[serial].add("receptor")
            elif rec == "HETATM" and resname in KNOWN_COFACTORS:
                receptor_lines.append(line)
                if serial is not None:
                    serial_bucket[serial].add("receptor")
            elif chain == binder_chain:
                binder_lines.append(line)
                if serial is not None:
                    serial_bucket[serial].add("binder")
            elif chain == " " or rec == "HETATM":
                binder_lines.append(line)
                if serial is not None:
                    serial_bucket[serial].add("binder")

    def _filter(bucket_name: str) -> List[str]:
        out: List[str] = []
        for cl in conect_raw:
            try:
                serials = [int(p) for p in cl[6:].split()]
            except ValueError:
                continue
            if serials and all(bucket_name in serial_bucket.get(s, set()) for s in serials):
                out.append(cl)
        return out

    complex_lines.extend(_filter("complex"))
    receptor_lines.extend(_filter("receptor"))
    binder_lines.extend(_filter("binder"))
    return complex_lines, receptor_lines, binder_lines


def _write_pdb(lines: List[str], path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
        if not lines or not lines[-1].startswith("END"):
            fh.write("END\n")


# ==========================================
# OpenMM energy helpers
# ==========================================
def _get_platform() -> mm.Platform:
    """Select platform. Env UPDD_1TRAJ_PLATFORM or OPENMM_PLATFORM wins; default
    'Reference' — deterministic CPU for debugging; avoids stealing GPU from the
    running Quick Win 1 calibration job.
    """
    pref = (
        os.environ.get("UPDD_1TRAJ_PLATFORM")
        or os.environ.get("OPENMM_PLATFORM")
        or "Reference"
    )
    for name in (pref, "Reference", "CPU"):
        try:
            return mm.Platform.getPlatformByName(name)
        except Exception:
            continue
    return mm.Platform.getPlatformByName("CPU")


def _build_system(pdb_path: str, ff: ForceField):
    """Load PDB and return (modeller, system) for a single bucket.

    No cyclic-bond injection here: Phase-3 v1 is WT-only, so HTC / SS wiring
    is not needed. If the project later extends this to cyclic WT variants,
    wire in `run_mmgbsa._apply_cyclic_bonds_mmgbsa` the same way.
    """
    pdb = PDBFile(pdb_path)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    # Snapshots coming out of the MM-GBSA pipeline already have explicit H; a
    # defensive addHydrogens keeps the receptor subset safe if any residue
    # lost an H during the split.
    modeller.addHydrogens(ff)
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=HBonds,
        soluteDielectric=1.0,
        solventDielectric=78.5,
    )
    return modeller, system


def _evaluate_energy_kcal(modeller, system, positions, platform) -> float:
    """Build a throwaway Simulation, inject the given positions, return
    potential energy in kcal/mol. No minimization here (the caller decides).
    """
    integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
    sim = Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(positions)
    state = sim.context.getState(getEnergy=True)
    e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    return e_kj / 4.184


# ==========================================
# 1-trajectory MM-GBSA core
# ==========================================
def run_1traj_mmgbsa(snap_pdb: str, outputdir: str, receptor_chain: str,
                     binder_chain: str, ncaa_elem: str, minimize_iter: int,
                     minimize_tol: float) -> Dict[str, Any]:
    """Perform the 1-snapshot 1-trajectory MM-GBSA on `snap_pdb`.

    Strategy: subsystem-rebuild. Split atom lines into complex/receptor/binder
    buckets → write three tmp PDBs → build three Systems with the SAME
    amber14 + GBn2 ForceField. Minimize ONLY the complex (bound state).
    For the receptor and binder subsystems, extract the bound-geometry
    positions of the corresponding atoms directly from the minimized complex
    (by order-preserving correspondence of the PDB lines) and evaluate their
    energies statically.

    This captures interaction + polarization differences while cancelling
    bound-state strain on both sides (cf. Genheden & Ryde 2015 §2.1).
    """
    os.makedirs(outputdir, exist_ok=True)
    tmp_dir = os.path.join(outputdir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    complex_lines, receptor_lines, binder_lines = classify_atoms(
        snap_pdb, receptor_chain, binder_chain
    )
    if not receptor_lines:
        raise RuntimeError(f"Receptor (chain {receptor_chain}) empty in {snap_pdb}")
    if not binder_lines:
        raise RuntimeError(f"Binder (chain {binder_chain}) empty in {snap_pdb}")

    stem = os.path.basename(snap_pdb).replace(".pdb", "")
    cplx_pdb = os.path.join(tmp_dir, f"{stem}_complex.pdb")
    recv_pdb = os.path.join(tmp_dir, f"{stem}_receptor.pdb")
    lig_pdb = os.path.join(tmp_dir, f"{stem}_binder.pdb")
    _write_pdb(complex_lines, cplx_pdb)
    _write_pdb(receptor_lines, recv_pdb)
    _write_pdb(binder_lines, lig_pdb)

    ff_files = [
        "amber14-all.xml",
        "amber/tip3p_HFE_multivalent.xml",  # Li-Merz metal templates (F1 safety)
        "implicit/gbn2.xml",
    ]
    ff = ForceField(*ff_files)
    platform = _get_platform()

    try:
        # --- Complex: build → minimize → record minimized positions + energy
        cplx_modeller, cplx_system = _build_system(cplx_pdb, ff)
        cplx_system = apply_gb_radius_override(
            cplx_system, cplx_modeller.topology, ncaa_elem, si_radius=2.10
        )

        integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
        sim = Simulation(cplx_modeller.topology, cplx_system, integrator, platform)
        sim.context.setPositions(cplx_modeller.positions)
        sim.minimizeEnergy(maxIterations=minimize_iter, tolerance=minimize_tol)

        min_state = sim.context.getState(getEnergy=True, getPositions=True)
        e_complex = min_state.getPotentialEnergy().value_in_unit(
            unit.kilojoules_per_mole) / 4.184
        cplx_min_positions = min_state.getPositions(asNumpy=False)

        # Build a lookup from (chain, resSeq, atomName) → index in complex
        # topology so we can transfer minimized coordinates to receptor/binder
        # subsystems cleanly (atom indices differ because subsystems drop
        # the other chain).
        def _atom_key(atom) -> Tuple[str, str, str, str]:
            res = atom.residue
            chain_id = res.chain.id if res.chain is not None else " "
            return (chain_id, res.name, str(res.id), atom.name)

        cplx_index_map = {_atom_key(a): a.index for a in cplx_modeller.topology.atoms()}

        # --- Receptor subsystem: static energy at bound (minimized) coords
        recv_modeller, recv_system = _build_system(recv_pdb, ff)
        recv_system = apply_gb_radius_override(
            recv_system, recv_modeller.topology, "none", si_radius=2.10
        )
        recv_positions = _transfer_positions(
            recv_modeller, cplx_index_map, cplx_min_positions,
            label=f"{stem} receptor",
        )
        e_receptor = _evaluate_energy_kcal(
            recv_modeller, recv_system, recv_positions, platform
        )

        # --- Binder subsystem: static energy at bound (minimized) coords
        lig_modeller, lig_system = _build_system(lig_pdb, ff)
        lig_system = apply_gb_radius_override(
            lig_system, lig_modeller.topology, ncaa_elem, si_radius=2.10
        )
        lig_positions = _transfer_positions(
            lig_modeller, cplx_index_map, cplx_min_positions,
            label=f"{stem} binder",
        )
        e_ligand = _evaluate_energy_kcal(
            lig_modeller, lig_system, lig_positions, platform
        )

        delta_g = e_complex - e_receptor - e_ligand
        return {
            "snapshot": stem,
            "method": "1-traj (Fix 4, subsystem-rebuild)",
            "minimize_iter": minimize_iter,
            "minimize_tol_kj_per_mol_nm": minimize_tol,
            "platform": platform.getName(),
            "e_complex_kcal": e_complex,
            "e_receptor_kcal": e_receptor,
            "e_ligand_kcal": e_ligand,
            "delta_g_kcal": delta_g,
            "favorable": delta_g < 0,
        }
    finally:
        for tmp in (cplx_pdb, recv_pdb, lig_pdb):
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except OSError:
                    pass


def _transfer_positions(sub_modeller, cplx_index_map: Dict[Tuple[str, str, str, str], int],
                        cplx_positions, label: str):
    """Copy minimized-complex coordinates into the subsystem's atom order.

    addHydrogens may have introduced new H atoms in the subsystem that were
    already present in the complex (no-op) or, in edge cases, genuinely new
    atoms. Missing keys fall back to the subsystem's own loaded positions.
    """
    sub_positions = list(sub_modeller.positions)
    n_matched = 0
    n_fallback = 0
    for atom in sub_modeller.topology.atoms():
        res = atom.residue
        chain_id = res.chain.id if res.chain is not None else " "
        key = (chain_id, res.name, str(res.id), atom.name)
        if key in cplx_index_map:
            sub_positions[atom.index] = cplx_positions[cplx_index_map[key]]
            n_matched += 1
        else:
            n_fallback += 1
    if n_fallback > 0:
        print(f"  [{label}] transfer_positions: {n_matched} matched, "
              f"{n_fallback} fallback (likely addHydrogens edge cases)")
    return sub_positions


# ==========================================
# 3-trajectory baseline lookup
# ==========================================
def lookup_3traj_baseline(snap_pdb: str) -> Optional[Dict[str, Any]]:
    """If `snap_pdb` matches outputs/{run}/snapshots/{stem}.pdb, load the
    sibling `mmgbsa_results/mmgbsa_summary.json` and return the result entry
    whose `snapshot` equals {stem}. Returns None if any step fails.
    """
    abs_path = os.path.abspath(snap_pdb)
    m = _SNAP_PATH_RE.match(abs_path)
    if not m:
        return None
    run = m.group("run")
    stem = m.group("stem")
    idx = abs_path.rfind(f"/outputs/{run}/")
    if idx < 0:
        return None
    run_root = abs_path[: idx + len(f"/outputs/{run}")]
    summary_path = os.path.join(run_root, "mmgbsa_results", "mmgbsa_summary.json")
    if not os.path.exists(summary_path):
        return None
    try:
        with open(summary_path, encoding="utf-8") as fh:
            summary = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None
    for res in summary.get("results", []):
        if res.get("snapshot") == stem:
            return {
                "snapshot": stem,
                "summary_path": summary_path,
                "e_complex_kcal": res.get("e_complex_kcal"),
                "e_receptor_kcal": res.get("e_receptor_kcal"),
                "e_ligand_kcal": res.get("e_ligand_kcal"),
                "delta_g_kcal": res.get("delta_g_kcal"),
            }
    return None


# ==========================================
# CLI
# ==========================================
def _format_cell(v: Optional[float]) -> str:
    return f"{v:>14.4f}" if isinstance(v, (int, float)) else f"{'n/a':>14}"


def _format_delta(a: Optional[float], b: Optional[float]) -> str:
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return f"{a - b:>+14.4f}"
    return f"{'n/a':>14}"


def print_comparison_table(one_traj: Dict[str, Any],
                           three_traj: Optional[Dict[str, Any]]) -> None:
    stem = one_traj["snapshot"]
    print(f"\nSnap: {stem}\n")
    hdr = f"  {'':<12} | {'3-traj (baseline)':>18} | {'1-traj (Fix 4)':>16} | {'Δ(1-3)':>14}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    rows = [
        ("E_complex", "e_complex_kcal"),
        ("E_receptor", "e_receptor_kcal"),
        ("E_ligand", "e_ligand_kcal"),
        ("ΔG", "delta_g_kcal"),
    ]
    for label, key in rows:
        base = three_traj.get(key) if three_traj else None
        new = one_traj.get(key)
        print(f"  {label:<12} | {_format_cell(base):>18} | {_format_cell(new):>16} "
              f"| {_format_delta(new, base):>14}")
    if three_traj is None:
        print("\n  [note] 3-traj baseline not located — outputs/{run}/mmgbsa_results/"
              "mmgbsa_summary.json missing or stem mismatch.")


def main():
    parser = argparse.ArgumentParser(
        description="1-trajectory MM-GBSA (Fix 4) vs 3-trajectory comparison "
                    "for a single snapshot. WT only (v1)."
    )
    parser.add_argument("--snap_pdb", required=True,
                        help="Pre-stripped complex snapshot PDB "
                             "(chain A=receptor, chain B=binder)")
    parser.add_argument("--outputdir", required=True,
                        help="Output directory for the JSON + tmp PDBs")
    parser.add_argument("--receptor_chain", default="A",
                        help="Receptor chain ID (default: A)")
    parser.add_argument("--binder_chain", default="B",
                        help="Binder chain ID (default: B)")
    parser.add_argument("--target_id", default="",
                        help="target_card stem (e.g., 2QKI) — currently "
                             "informational only; WT snaps need no cofactor "
                             "template registration")
    parser.add_argument("--ncaa_elem", default="none",
                        help="ncAA key element passed to GB radius override "
                             "(default: none). Leave 'none' for WT.")
    parser.add_argument("--ncaa_label", default="none",
                        help="ncAA variant label (e.g., 'Cp4'). Must be "
                             "'none' in v1 — ncAA parameterization is "
                             "out of scope for Phase 3 test.")
    parser.add_argument("--minimize_iter", type=int, default=5000,
                        help="Complex minimization iterations (P1: ≥5000)")
    parser.add_argument("--minimize_tol", type=float, default=1.0,
                        help="Minimization tolerance in kJ/mol/nm (default: 1.0)")
    args = parser.parse_args()

    if args.ncaa_label and args.ncaa_label.lower() != "none":
        raise NotImplementedError(
            f"ncaa_label={args.ncaa_label!r} not supported in v1. "
            "Phase 3 1-traj test is WT-only (ncAA parameter wiring requires "
            "the full MD run's params_manifest/hydrogens/mol2 bundle, which "
            "is not reachable from a standalone snapshot). Rerun with "
            "--ncaa_label none on a WT snapshot."
        )

    if not os.path.exists(args.snap_pdb):
        raise FileNotFoundError(f"snap_pdb not found: {args.snap_pdb}")

    os.makedirs(args.outputdir, exist_ok=True)

    print(f"[mmgbsa_1traj_compare] snap   : {args.snap_pdb}")
    print(f"[mmgbsa_1traj_compare] out    : {args.outputdir}")
    print(f"[mmgbsa_1traj_compare] chains : receptor={args.receptor_chain} "
          f"binder={args.binder_chain}")
    print(f"[mmgbsa_1traj_compare] min    : {args.minimize_iter} iter, "
          f"tol={args.minimize_tol} kJ/mol/nm")

    one_traj = run_1traj_mmgbsa(
        snap_pdb=args.snap_pdb,
        outputdir=args.outputdir,
        receptor_chain=args.receptor_chain,
        binder_chain=args.binder_chain,
        ncaa_elem=args.ncaa_elem,
        minimize_iter=args.minimize_iter,
        minimize_tol=args.minimize_tol,
    )
    three_traj = lookup_3traj_baseline(args.snap_pdb)

    print_comparison_table(one_traj, three_traj)

    out_json = os.path.join(args.outputdir, "mmgbsa_1traj_compare.json")
    payload = {
        "input_pdb": os.path.abspath(args.snap_pdb),
        "target_id": args.target_id,
        "ncaa_label": args.ncaa_label,
        "ncaa_elem": args.ncaa_elem,
        "receptor_chain": args.receptor_chain,
        "binder_chain": args.binder_chain,
        "one_traj": one_traj,
        "three_traj_baseline": three_traj,
        "regime": "ranking-only (R-11)",
        "reference": "Genheden & Ryde 2015 §2.1 "
                     "DOI:10.1517/17460441.2015.1032250",
    }
    if isinstance(one_traj.get("delta_g_kcal"), (int, float)) and \
       three_traj and isinstance(three_traj.get("delta_g_kcal"), (int, float)):
        payload["delta_delta_g_1traj_minus_3traj_kcal"] = (
            one_traj["delta_g_kcal"] - three_traj["delta_g_kcal"]
        )
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    print(f"\n  saved: {out_json}")


if __name__ == "__main__":
    main()
