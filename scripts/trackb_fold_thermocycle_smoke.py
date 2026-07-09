#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B — barnase Ile96->Ala FOLDING-thermocycle STAGE-0 ASSERT smoke (CPU).

STAGE-0 of the Gate-A large-effect engine control
(analysis/ncaa_funnel_gate_governance_20260702.md §Phase-1; scientific review GO-WITH-COND
verdict_gateA_largeeffect_control_design_20260702). Exercises the ADDITIVE
system-prep layer that lets the FROZEN two-copy ATS in-place RBFE engine build a
CANONICAL single-scaffold folding control through the IDENTICAL code path
(``atm_trackB_setup.build_inplace_res4_twocopy_system`` / the frozen soft-core
canon / H18 ladder / UWHAM are ALL untouched).

What this STAGE-0 (solvate=False, assert-only, CPU) proves — for BOTH the barnase
folded scaffold AND the Ac-Ile-NMe tripeptide reference, on leg='free',
Lambda1=Lambda2=0, auto_search_displacement=True:

  (C4)  the two-copy box builds (two full copies resident, n_atoms ~ 2x);
  (MC1) full-residue net-charge sanity (Ile/Ala both net-0 -> passes);
  (C6)  the two copies are spatially SEPARATED (auto-search escalates d to clear
        two barnase copies; min-image / raw solute separation >= accept);
  (MC2) the Ile->Ala ACYCLIC single-attach connected group (CG1,CG2 -> CB roots +
        CD1 chained off CG1, no ring) is certified via the NEW
        chained_group_certified branch;
  (MC3) the disulfide guard SKIPS/PASSES (barnase / Ac-Ile-NMe carry zero CYS);
  (R2)  per-copy seed geometry is non-clashing;
  (C6e/C8) endpoint-equivalence (u0 == full E) + bulk-copy decouple.

It ALSO deliberately EXPOSES the two predicted failure modes the patch opens
(patched vs pre-patch contrast, both run POST-patch so the contrast is clean):
  (i)  MC2 shape wall: an Ile->Ala spec WITHOUT chained_group_certified (== the
       pre-patch behaviour, the flag did not exist) classifies as 'unsupported'
       and the build fail-louds — while the certified spec passes.
  (ii) MC3 disulfide wall: the PRE-PATCH assert body (``len(disulfides) < 2 ->
       raise``) run on the built box raises 'detected 0', while the patched
       (cysteine-conditional) assert skips/passes.

What this does NOT prove (R-18, honest scope): the converged ddG_fold (needs the
H18 ladder + UWHAM — the paired production campaign, post-validation); that the
buried-core perturbation keeps overlap at deep lambda (window escalation may be
needed). Ranking-only (R-11); a PREDICTION-path build test, NOT a converged ddG.

The barnase scaffold is prepared here from a fetched WT barnase PDB (default
1BNI) ONLY for the assert-only smoke; the PRODUCTION scaffold is an equilibrated,
CRYST1-bearing, solvated MD final.pdb (a separate DATA step, out of this code's
scope). PDBFixer renumbers the crystal residues sequentially, so the core Ile
(crystal Ile96) is located STRUCTURALLY by the unique Leu-Ile-Tyr motif rather
than a hard-coded number.

DOI references:
  - Gallicchio 2025 J Chem Inf Model, ATS DOI 10.1021/acs.jcim.5c00207
  - Prevost/Serrano 1996 Protein Eng 9(3):273, DOI 10.1093/protein/9.3.273
    (barnase I96A amber-family FEP ddG_fold ~3.9 kcal/mol)
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple

_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)
if _PROJ not in sys.path:
    sys.path.insert(0, _PROJ)


def _load_ats():
    path = os.path.join(_PROJ, "utils", "atm_trackB_setup.py")
    spec = importlib.util.spec_from_file_location("atm_trackB_setup", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------
# Scaffold prep (DATA helpers for the assert-only smoke; NOT the production
# MD prep — production uses an equilibrated solvated CRYST1-bearing final.pdb).
# --------------------------------------------------------------------------
def prep_barnase_single_chain(out_pdb: str, pdbid: str = "1BNI",
                              chain_keep: str = "A") -> str:
    """Fetch a WT barnase PDB, keep one protein chain, drop het/water, complete
    heavy atoms + hydrogens, and stamp the chain id. H-complete single chain."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(pdbid=pdbid)
    keep_idx = None
    for i, ch in enumerate(fixer.topology.chains()):
        names = [r.name for r in ch.residues()]
        if ch.id == chain_keep and any(
                n in ("ILE", "ALA", "LEU", "GLY", "TYR") for n in names):
            keep_idx = i
            break
    if keep_idx is None:
        raise RuntimeError(
            "prep_barnase_single_chain: no protein chain %r in %s"
            % (chain_keep, pdbid))
    fixer.removeChains(
        [i for i, _ in enumerate(fixer.topology.chains()) if i != keep_idx])
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do not insert internal gaps
    fixer.removeHeterogens(False)       # drop waters / ions
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    tmp = out_pdb + ".w.pdb"
    with open(tmp, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    with open(tmp) as s, open(out_pdb, "w") as d:
        for ln in s:
            if ln[:6] in ("ATOM  ", "HETATM"):
                d.write(ln[:21] + chain_keep + ln[22:])
            else:
                d.write(ln)
    return out_pdb


def prep_ac_ile_nme(out_pdb: str, chain: str = "A") -> str:
    """Build an Ac-Ile-NMe (ACE-ILE-NME) tripeptide via tleap (ff14SB templates,
    idealized coords + hydrogens). ILE is residue 2."""
    tleap = None
    for cand in ("/home/san/miniconda3/envs/qmmm/bin/tleap",):
        if os.path.isfile(cand):
            tleap = cand
            break
    if tleap is None:
        tleap = "tleap"
    work = tempfile.mkdtemp(prefix="fold_tripep_")
    raw = os.path.join(work, "tripep_raw.pdb")
    inp = os.path.join(work, "build.in")
    with open(inp, "w") as fh:
        fh.write("source leaprc.protein.ff14SB\n")
        fh.write("mol = sequence { ACE ILE NME }\n")
        fh.write("savepdb mol %s\n" % raw)
        fh.write("quit\n")
    r = subprocess.run([tleap, "-f", inp], capture_output=True, text=True,
                       cwd=work)
    if not os.path.isfile(raw):
        raise RuntimeError("tleap Ac-Ile-NMe build failed (rc=%d):\n%s\n%s"
                           % (r.returncode, r.stdout[-800:], r.stderr[-400:]))
    with open(raw) as s, open(out_pdb, "w") as d:
        for ln in s:
            if ln[:6] in ("ATOM  ", "HETATM"):
                d.write(ln[:21] + chain + ln[22:])
            else:
                d.write(ln)
    return out_pdb


def find_core_ile_resnum(pdb_path: str, prev: str = "LEU",
                         nxt: str = "TYR") -> int:
    """Locate the buried core Ile by the unique ``prev``-ILE-``nxt`` sequence
    motif (barnase Ile96 = Leu95-Ile96-Tyr97), robust to PDBFixer renumbering."""
    from openmm.app import PDBFile
    rl = [(int(r.id), r.name) for r in PDBFile(pdb_path).topology.residues()]
    hits = [rl[i][0] for i in range(1, len(rl) - 1)
            if rl[i][1] == "ILE" and rl[i - 1][1] == prev and rl[i + 1][1] == nxt]
    if not hits:
        raise RuntimeError(
            "find_core_ile_resnum: no %s-ILE-%s motif in %s" % (prev, nxt, pdb_path))
    return hits[0]


# --------------------------------------------------------------------------
# STAGE-0 assert-only build for one scaffold.
# --------------------------------------------------------------------------
def run_stage0_for_scaffold(ats, scaffold_pdb: str, resnum: int,
                            binder_chain: str, tag: str,
                            platform_name: str = "Reference") -> Dict[str, Any]:
    """Build the two-copy box for one scaffold (leg='free', solvate=False,
    auto_search_displacement) + certify every gate + endpoint-equivalence."""
    li = ats.resolve_fold_leg_inputs(scaffold_pdb)
    spec = ats.make_ile_ala_mutation_spec(resnum, name="ile_ala_%s" % tag)
    assert spec.shape == "acyclic_connected_group", \
        "certified Ile->Ala must classify acyclic_connected_group, got %r" % spec.shape

    build = ats.build_inplace_res4_twocopy_system(
        leg="free", binder_chain=binder_chain, solvate=False,
        auto_search_displacement=True, spec=spec, leg_inputs=li)

    assert build["outcome"] == "twocopy_attached", build.get("outcome")
    mc2 = build["mc2_methyl_bonded"]
    mc3 = build["mc3_disulfide"]
    sep = build["separation"]
    seed_a = build["seed_assert"]
    assert mc2["passed"] and mc2["shape"] == "acyclic_connected_group", mc2
    assert mc3["passed"] and mc3.get("skipped_reason") == \
        "disulfide_free_scaffold_no_cysteines", mc3
    assert sep["passed"], sep
    assert seed_a["passed"], seed_a

    eq = ats.check_twocopy_endpoint_equivalence(build, platform_name=platform_name)
    assert eq["overall_pass"], eq

    return {
        "tag": tag,
        "scaffold_pdb": scaffold_pdb,
        "resnum": resnum,
        "binder_chain": binder_chain,
        "shape": spec.shape,
        "n_atoms": build["fused_build"]["n_atoms"],
        "n_copy1": build["fused_build"]["n_copy1"],
        "mc2": {k: mc2.get(k) for k in
                ("shape", "side", "attach_bonded_roots",
                 "n_heavies_certified", "n_subgraph_edges", "passed")},
        "mc3": mc3,
        "separation_solute_min_nm": sep.get("solute_solute_min_sep_nm"),
        "separation_image_solute_min_nm": sep.get("image_solute_min_sep_nm"),
        "r2_seed_passed": seed_a["passed"],
        "endpoint_equivalence": eq.get("endpoint_equivalence"),
        "bulk_copy_decoupled": eq.get("bulk_copy_decoupled"),
        "energies_kcal": eq.get("energies_kcal"),
        "overall_pass": True,
    }


def demonstrate_walls(ats, scaffold_pdb: str, resnum: int,
                      binder_chain: str, tag: str) -> Dict[str, Any]:
    """Reproduce the two PRE-PATCH failure walls the patch opens (post-patch, so
    the contrast is clean): (i) uncertified Ile->Ala shape -> MC2 'unsupported';
    (ii) the pre-patch MC3 body (>=2 unconditional) -> 'detected 0'."""
    import dataclasses

    # (i) MC2 shape wall: default-flag (chained_group_certified=False) == pre-patch.
    spec_unc = dataclasses.replace(
        ats.make_ile_ala_mutation_spec(resnum, name="ile_ala_%s_unc" % tag),
        chained_group_certified=False)
    unc_shape = spec_unc.shape
    li = ats.resolve_fold_leg_inputs(scaffold_pdb)
    mc2_raise = None
    try:
        ats.build_inplace_res4_twocopy_system(
            leg="free", binder_chain=binder_chain, solvate=False,
            auto_search_displacement=True, spec=spec_unc, leg_inputs=li)
    except ValueError as exc:
        mc2_raise = str(exc)

    # (ii) MC3 disulfide wall: build the box (patched), then run the PRE-PATCH body
    #      (>=2 unconditional) on it vs the patched assert.
    spec_ok = ats.make_ile_ala_mutation_spec(resnum, name="ile_ala_%s" % tag)
    build = ats.build_inplace_res4_twocopy_system(
        leg="free", binder_chain=binder_chain, solvate=False,
        auto_search_displacement=True, spec=spec_ok, leg_inputs=li)
    disulfides = build["fused_build"].get("disulfides") or []
    # PRE-PATCH MC3 body (verbatim: the unconditional >=2 guard).
    prepatch_mc3_raise = None
    if len(disulfides) < 2:
        prepatch_mc3_raise = (
            "MC3 two-copy FAIL: expected 2 cyclic_ss disulfides (one per copy), "
            "detected %d." % (len(disulfides),))
    # PATCHED MC3 (cysteine-conditional) on the same box.
    patched_mc3 = ats.assert_twocopy_disulfides(build["fused_build"])

    return {
        "mc2_wall": {
            "uncertified_shape": unc_shape,
            "uncertified_is_unsupported": (unc_shape == "unsupported"),
            "build_raised": mc2_raise is not None,
            "raise_is_unsupported_shape": bool(
                mc2_raise and "unsupported mutation shape" in mc2_raise),
        },
        "mc3_wall": {
            "n_disulfides_detected": len(disulfides),
            "prepatch_would_raise": prepatch_mc3_raise is not None,
            "prepatch_message": prepatch_mc3_raise,
            "patched_passed": patched_mc3.get("passed"),
            "patched_skipped_reason": patched_mc3.get("skipped_reason"),
        },
    }


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Barnase Ile96->Ala folding-thermocycle STAGE-0 assert smoke.")
    p.add_argument("--barnase-pdb", default=None,
                   help="prepared H-complete single-chain barnase PDB (else "
                        "fetch + prep from --barnase-pdbid).")
    p.add_argument("--barnase-pdbid", default="1BNI",
                   help="PDB id to fetch when --barnase-pdb is not supplied.")
    p.add_argument("--tripep-pdb", default=None,
                   help="prepared Ac-Ile-NMe PDB (else build via tleap).")
    p.add_argument("--skip-barnase", action="store_true",
                   help="run the tripeptide scaffold only (no network fetch).")
    p.add_argument("--skip-tripep", action="store_true",
                   help="run the barnase scaffold only (no tleap).")
    p.add_argument("--work-dir", default=None,
                   help="dir for the prepared scaffolds (default: a temp dir).")
    p.add_argument("--platform", default="Reference",
                   help="OpenMM platform for the endpoint-equivalence eval.")
    p.add_argument("--json-out", default=None, help="write the summary JSON here.")
    args = p.parse_args(argv)

    ats = _load_ats()
    work = args.work_dir or tempfile.mkdtemp(prefix="fold_stage0_")
    os.makedirs(work, exist_ok=True)

    summary: Dict[str, Any] = {"scaffolds": [], "walls": {}, "all_pass": False}

    scaffolds: List[Tuple[str, str, int, str]] = []  # (tag, pdb, resnum, chain)

    if not args.skip_barnase:
        barn = args.barnase_pdb
        if barn is None:
            barn = prep_barnase_single_chain(
                os.path.join(work, "barnase_chainA.pdb"),
                pdbid=args.barnase_pdbid)
        core = find_core_ile_resnum(barn)
        scaffolds.append(("barnase", barn, core, "A"))

    if not args.skip_tripep:
        tri = args.tripep_pdb
        if tri is None:
            tri = prep_ac_ile_nme(os.path.join(work, "ac_ile_nme.pdb"))
        scaffolds.append(("tripep", tri, 2, "A"))

    if not scaffolds:
        print("ERROR: nothing to run (both scaffolds skipped).")
        return 2

    ok = True
    for tag, pdb, resnum, chain in scaffolds:
        try:
            res = run_stage0_for_scaffold(
                ats, pdb, resnum, chain, tag, platform_name=args.platform)
            summary["scaffolds"].append(res)
            print("[STAGE-0 %-8s] PASS  n_atoms=%d resnum=%d shape=%s "
                  "sep_min=%.3f nm eq=%s"
                  % (tag, res["n_atoms"], resnum, res["shape"],
                     res["separation_solute_min_nm"], res["overall_pass"]))
        except Exception as exc:  # noqa: BLE001 - smoke surfaces the cause
            ok = False
            summary["scaffolds"].append({"tag": tag, "error": str(exc)})
            print("[STAGE-0 %-8s] FAIL  %s" % (tag, exc))

    # Wall demonstration on the first available scaffold (path is scaffold-neutral).
    wtag, wpdb, wres, wchain = scaffolds[0]
    try:
        walls = demonstrate_walls(ats, wpdb, wres, wchain, wtag)
        summary["walls"] = walls
        w2 = walls["mc2_wall"]
        w3 = walls["mc3_wall"]
        w2_ok = (w2["uncertified_is_unsupported"] and w2["build_raised"]
                 and w2["raise_is_unsupported_shape"])
        w3_ok = (w3["prepatch_would_raise"] and w3["patched_passed"]
                 and w3["patched_skipped_reason"]
                 == "disulfide_free_scaffold_no_cysteines")
        ok = ok and w2_ok and w3_ok
        print("[WALL mc2 shape ] pre-patch(uncertified)=%r -> build raises "
              "unsupported=%s ; patched(certified) builds"
              % (w2["uncertified_shape"], w2["raise_is_unsupported_shape"]))
        print("[WALL mc3 disulf] pre-patch(>=2) would raise 'detected %d' ; "
              "patched skips=%s"
              % (w3["n_disulfides_detected"], w3["patched_skipped_reason"]))
    except Exception as exc:  # noqa: BLE001
        ok = False
        summary["walls"] = {"error": str(exc)}
        print("[WALL           ] FAIL  %s" % exc)

    summary["all_pass"] = ok
    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(summary, fh, indent=2)
    print("\nSTAGE-0 ALL_PASS =", ok)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
