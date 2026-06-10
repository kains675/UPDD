#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B v2.1 — upstream ``atom_openmm.make_atm_system_from_rcpt_lig`` wrapper.

Q2-B PRIMARY: build the ABFE system with the canonical upstream pattern
(binder pre-translated to the displaced position BEFORE Modeller.addSolvent,
so the binding site is empty and water fills it cleanly at solvation time).

Our previous v2 path (``utils/atm_trackB_setup.build_leg_system``) solvated
the binder at the binding site — the INVERTED pattern. Result: bound-leg
first-step NaN at λ=0.5 (the alchemical mid-point) because the binder's
displaced "ghost" position was buried in water → grad(u1) divergent.

This wrapper splits our prepared 2QKI {Cp4,WT} complex PDB into receptor +
binder PDB inputs, registers the MTR ncAA force field (hybrid GAFF2 +
amber14SB-frozen backbone XML + MTR_hydrogens.xml) BEFORE calling upstream,
and lets ``make_atm_system_from_rcpt_lig.make_system`` perform the actual
Modeller.addSolvent + ForceField.createSystem with the upstream's
binder-pre-displaced pattern.

The output (system XML + topology PDB) is consumed verbatim by
``abfe_structprep`` and ``abfe_production`` — no custom ATMForce wrapping
required (upstream does this in ``OMMSystemABFE``).

# Inputs

The 4 leg systems are built from 2 source complexes:

* Cp4: ``outputs/2QKI_Cp4_hybrid_calib_s7/_md_input/2QKI_Cp4.pdb`` (chain A
  = MASP-2 receptor 1-9897, chain B = cyclic compstatin with MTR at
  residue 4 starting at atom 9898)
* WT:  ``outputs/2QKI_WT_calib_s7/_md_input/2QKI_WT.pdb`` (chain A = same
  receptor, chain B = cyclic compstatin with TRP at residue 4)

Each source produces two leg builds:

* ``bound`` = receptor + binder co-solvated (upstream pre-displaces the
  binder into bulk before addSolvent)
* ``free``  = binder only in water (free leg of the thermodynamic cycle)

For the free leg we pass ``receptorinFile=binder.pdb`` and OMIT
``LIG1inFile`` — upstream treats the receptor as the sole non-solvent unit
and solvates around it. No displacement applied (the free-leg cycle has
no binding site).

# MTR force field registration

Upstream's ``ForceField(proteinforcefield, solventforcefield)`` accepts a
SPACE-SEPARATED list in its first arg (Python's ``ForceField`` ctor
unpacks ``*args``). We pass:

    ``--proteinForceField "amber14-all.xml params/MTR_gaff2_hybrid.xml"``

so the hybrid MTR XML is loaded alongside amber14SB. We then load the
MTR_hydrogens.xml via ``Modeller.loadHydrogenDefinitions`` BEFORE
addHydrogens — required so the ncAA N-methyl group (HM1/HM2/HM3 on CM)
positions are placed correctly. Upstream does NOT call addHydrogens for
PDB inputs, so we monkey-patch the call to insert hydrogen-definition
loading at the right point (alternatively, we pre-protonate the PDB and
skip addHydrogens — see ``add_hydrogens_to_pdb`` below).

# Displacement choice

We use 2.5 nm (= 25 Å) +x by default — same as v2 — well above the 10 Å
minimum (compmolbiophys docs) and well above the 22 Å recommended for
binder volume ≲ 1000 Å³. The 5070Ti can handle the larger box.

# Cyclic-SS disulfide

The compstatin macrocycle (CYS2-CYS12 SG-SG) is auto-detected and
committed by ``utils.atm_trackB_setup.detect_disulfide_pair`` on the
output topology — we apply it post-build (after upstream creates the
Modeller-+ ForceField system).

# Cross-references

* Upstream source: ``/home/san/miniconda3/envs/atm/lib/python3.11/site-packages/atom_openmm/make_atm_system_from_rcpt_lig.py``
* Replaced custom path: ``utils/atm_trackB_setup.build_leg_system``
  (archived: ``_archive/trackb_v2_custom_build_20260531/``)
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from typing import Optional, List, Dict, Tuple, Any

import numpy as np

_PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "utils"))


# ---------------------------------------------------------------------------
# Canonical UPDD asset paths (mirror utils/atm_trackB_setup.py)
# ---------------------------------------------------------------------------
HYBRID_MTR_XML = os.path.join(_PROJ_ROOT, "params", "MTR_gaff2_hybrid.xml")
PROTEIN_FF_DEFAULT = "amber14-all.xml"
SOLVENT_FF_DEFAULT = "amber14/tip3p.xml"
DISPLACEMENT_NM_DEFAULT = (2.5, 0.0, 0.0)


# ---------------------------------------------------------------------------
# Source-PDB → receptor / binder split
# ---------------------------------------------------------------------------
def split_complex_to_receptor_binder(
    complex_pdb: str,
    receptor_out: str,
    binder_out: str,
    receptor_chain: str = "A",
    binder_chain: str = "B",
) -> Tuple[str, str]:
    """Split a 2-chain complex PDB into receptor + binder PDB files.

    Upstream ``make_atm_system_from_rcpt_lig.make_system`` takes the
    receptor file and the ligand file separately. Our prepared complex
    PDBs (``outputs/2QKI_*/_md_input/2QKI_*.pdb``) have chain A = receptor
    and chain B = binder; this writes them out as two PDBs. CRYST1 / END
    / MODEL records are dropped (upstream parses each file with
    ``openmm.app.PDBFile`` which infers box later).

    HETATM lines are preserved on whichever chain owns them (the MTR ncAA
    is HETATM on chain B).
    """
    rec_lines: List[str] = []
    bin_lines: List[str] = []
    with open(complex_pdb) as fh:
        for line in fh:
            tag = line[:6]
            if tag in ("ATOM  ", "HETATM"):
                chain = line[21]
                if chain == receptor_chain:
                    rec_lines.append(line)
                elif chain == binder_chain:
                    bin_lines.append(line)
            elif tag == "TER   ":
                # Mirror TER into whichever chain the preceding atom belonged.
                if rec_lines and rec_lines[-1][:6] in ("ATOM  ", "HETATM"):
                    rec_lines.append(line)
                if bin_lines and bin_lines[-1][:6] in ("ATOM  ", "HETATM"):
                    bin_lines.append(line)
    if not rec_lines:
        raise RuntimeError(
            f"No receptor-chain ({receptor_chain}) atoms in {complex_pdb}"
        )
    if not bin_lines:
        raise RuntimeError(
            f"No binder-chain ({binder_chain}) atoms in {complex_pdb}"
        )
    with open(receptor_out, "w") as fh:
        fh.writelines(rec_lines)
        fh.write("END\n")
    with open(binder_out, "w") as fh:
        fh.writelines(bin_lines)
        fh.write("END\n")
    return receptor_out, binder_out


# ---------------------------------------------------------------------------
# Hydrogen pre-build (so upstream sees a fully-protonated PDB)
# ---------------------------------------------------------------------------
def protonate_binder_with_mtr_xml(
    binder_pdb_in: str,
    binder_pdb_out: str,
    mtr_xml: str = HYBRID_MTR_XML,
    hydrogens_xml: Optional[str] = None,
    binder_chain: str = "B",
) -> Dict[str, Any]:
    """Add hydrogens to the binder PDB, registering MTR ncAA definitions first.

    Upstream ``make_atm_system_from_rcpt_lig`` does NOT call ``addHydrogens``
    on PDB inputs — it assumes the PDB is fully protonated. Our prepared
    binder PDBs from ``outputs/2QKI_*/_md_input/`` may be missing the
    N-methyl HM1/HM2/HM3 on MTR's CM atom. This function pre-runs
    ``Modeller.addHydrogens`` with the MTR XML registered, so when upstream
    reads the PDB it sees a complete H set.

    Idempotent: if the binder already has hydrogens on every heavy atom,
    addHydrogens is a no-op.

    Also commits ncAA peptide bonds + cyclic_ss disulfide (the
    OpenMM/PDBFile parser does not auto-infer ATOM↔HETATM peptide bonds at
    MTR junctions, and the cyclic disulfide needs an explicit bond).
    """
    from openmm.app import ForceField, Modeller, PDBFile
    from atm_trackB_setup import (
        add_missing_peptide_bonds_safe,
        commit_bond_if_missing,
        detect_disulfide_pair,
        inject_xml_internal_bonds,
    )

    ff_inputs = [PROTEIN_FF_DEFAULT, mtr_xml]
    ff = ForceField(*ff_inputs)

    pdb = PDBFile(binder_pdb_in)
    modeller = Modeller(pdb.topology, pdb.positions)

    # Inject MTR ncAA internal bonds from the XML (HETATM parser skips them).
    n_internal_added = inject_xml_internal_bonds(
        modeller.topology, [mtr_xml], xml_res_name="MTR"
    )

    # Inter-residue peptide bonds at the HETATM/ATOM junctions.
    n_peptide_added = add_missing_peptide_bonds_safe(
        modeller, binder_chain=binder_chain, max_cn_distance_nm=0.20
    )

    # cyclic_ss disulfide (CYS2-CYS12 SG-SG).
    disulfide = None
    pair = detect_disulfide_pair(modeller, binder_chain=binder_chain)
    if pair is not None:
        added = commit_bond_if_missing(modeller.topology, pair[0], pair[1])
        disulfide = {
            "sg1_res": pair[0].residue.id,
            "sg2_res": pair[1].residue.id,
            "bond_added": added,
        }

    # Register MTR hydrogen definitions if available (needed for HM1-3 on CM).
    if hydrogens_xml and os.path.isfile(hydrogens_xml):
        Modeller.loadHydrogenDefinitions(hydrogens_xml)
    modeller.addHydrogens(ff)

    with open(binder_pdb_out, "w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)

    return {
        "n_atoms": modeller.topology.getNumAtoms(),
        "n_internal_bonds_added": n_internal_added,
        "n_peptide_bonds_added": n_peptide_added,
        "disulfide": disulfide,
        "ff_inputs": ff_inputs,
    }


# ---------------------------------------------------------------------------
# Receptor pre-build (add receptor hydrogens once, reuse across endpoints)
# ---------------------------------------------------------------------------
def protonate_receptor(
    receptor_pdb_in: str,
    receptor_pdb_out: str,
) -> Dict[str, Any]:
    """Add hydrogens to the bare receptor PDB.

    The receptor (MASP-2 catalytic domain in chain A) is a standard
    amber14SB-compatible protein, no ncAA. addHydrogens here ensures all
    backbone/sidechain Hs are placed before upstream sees the file.
    """
    from openmm.app import ForceField, Modeller, PDBFile

    ff = ForceField(PROTEIN_FF_DEFAULT)
    pdb = PDBFile(receptor_pdb_in)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)

    with open(receptor_pdb_out, "w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)

    return {"n_atoms": modeller.topology.getNumAtoms()}


# ---------------------------------------------------------------------------
# ForceField multi-XML monkey-patch
# ---------------------------------------------------------------------------
# Upstream ``make_atm_system_from_rcpt_lig.make_system`` calls
# ``ForceField(proteinforcefield, solventforcefield)`` at module line 133,
# where each arg is a SINGLE XML filename. To register the MTR ncAA hybrid
# XML (which the binder needs) we wrap the ForceField symbol on the
# upstream module with a thin proxy that whitespace-splits multi-XML
# strings into separate positional ``*files`` for openmm.app.ForceField.
# Reverted via try/finally so no global state pollution.
# ---------------------------------------------------------------------------
def _make_multi_xml_forcefield_proxy(real_forcefield_cls):
    """Return a callable that splits whitespace-separated XML paths.

    ``proxy("a.xml b.xml", "c.xml")`` -> ``real_ForceField("a.xml", "b.xml", "c.xml")``.
    Splitting is whitespace-safe because XML filenames in our project
    layout never contain spaces.
    """
    def _proxy(*files):
        expanded: List[str] = []
        for f in files:
            if isinstance(f, str) and (" " in f or "\t" in f):
                expanded.extend(f.split())
            else:
                expanded.append(f)
        return real_forcefield_cls(*expanded)
    return _proxy


def _patched_make_system(**kwargs) -> None:
    """Call the vendored upstream make_system with the multi-XML proxy active.

    Uses ``_vendored_make_atm_system.make_system`` (a verbatim copy of
    upstream's PDB code path, with the openff dependency stripped — the
    ``atm`` conda env intentionally does not have openff). The proxy is
    installed on the vendored module's ForceField symbol for the duration
    of the call, then reverted. Guarantees revert even on upstream raise.

    Forwards ALL kwargs verbatim (including v0.9.19
    ``apply_modeller_pre_displacement`` — see vendored module for details).

    See ``scripts/_vendored_make_atm_system.py`` docstring for the
    rationale and the upstream provenance manifest.
    """
    import _vendored_make_atm_system as vendored
    from openmm.app import ForceField as RealForceField

    original = vendored.ForceField
    vendored.ForceField = _make_multi_xml_forcefield_proxy(RealForceField)
    try:
        vendored.make_system(**kwargs)
    finally:
        vendored.ForceField = original


# ---------------------------------------------------------------------------
# Upstream system build invocation
# ---------------------------------------------------------------------------
def build_bound_leg_via_upstream(
    receptor_pdb: str,
    binder_pdb: str,
    xml_out: str,
    pdb_out: str,
    displacement_nm: Tuple[float, float, float] = DISPLACEMENT_NM_DEFAULT,
    protein_ff: str = PROTEIN_FF_DEFAULT,
    solvent_ff: str = SOLVENT_FF_DEFAULT,
    mtr_xml: str = HYBRID_MTR_XML,
    ionic_strength_M: float = 0.15,
    hmass_amu: float = 1.0,
    direction: str = "dplus",
) -> Dict[str, Any]:
    """Build the bound-leg ABFE system via upstream ``make_system``.

    Upstream's ABFE path:
    1. Modeller.add(receptor) at binding-site coordinates.
    2. Modeller.add(binder) at binding-site coordinates.
    3. Translate binder coords by ``displacement`` (binder → bulk solvent).
    4. Bounding box from receptor + displaced binder positions + 2 nm padding.
    5. Modeller.addSolvent(forceField, boxVectors=..., ionicStrength=...).
       — addSolvent clash-detects the binder atoms, so water fills the
       empty binding site cleanly.
    6. forcefield.createSystem(modeller.topology, PME, 0.9 nm cutoff, ...).

    Per-direction support (2026-06-01 v0.9.19, (b+) corrected spec,
    [[trackb_per_direction_system_xml_rebuild_20260601]]):

    BOTH directions pass ``+displacement_nm`` (positive) to the vendored
    ``make_system``. The difference is the Modeller pre-displacement flag:

    ``direction="dplus"`` — sys_dplus.xml = "bound state base":
      * Modeller binder@x0 (binding site coords) — vendored deviation
        DISABLED via ``apply_modeller_pre_displacement=False``.
      * addSolvent fills +displacement region with water.
      * ATMForce per-particle displacement = +d_vec (from cntl
        ``DISPLACEMENT`` keyword at runtime, default).
      * Walker dynamics base (Direction=+1) = u0 = E(binder@x0) =
        bound state, STABLE.
      * u1 perturbation = E(binder@x0+d_vec) (in water region) —
        soft-cored as PERTURBATION, not as base.

    ``direction="dminus"`` — sys_dminus.xml = "dissociated state base":
      * Modeller binder pre-displaced to x0+d_vec (bulk water region)
        — vendored deviation ENABLED via
        ``apply_modeller_pre_displacement=True``.
      * addSolvent fills binding-site region (x0) with water.
      * ATMForce per-particle displacement = -d_vec (must be set at
        runtime by structprep cntl override; this builder writes the
        system with the SAME +displacement_nm for bbox/Modeller
        consistency — the runtime negate is what makes the alchemical
        perturbation point back to the binding site).
      * Walker dynamics base (Direction=+1) = u0 = E(binder@x0+d_vec) =
        dissociated state, STABLE.
      * u1 perturbation = E(binder@x0+d_vec + (-d_vec)) = E(binder@x0)
        — bound state, soft-cored as PERTURBATION.

    Both walkers use Direction=+1 (NOT ±1). This is the standard
    two-leg ATM ABFE structure (Azimi 2022 §2.3, ATS Khuttan 2024) where
    leg 1 starts bound + perturbs toward unbound, leg 2 starts
    dissociated + perturbs toward bound. They meet at the alchemical
    intermediate λ=0.5; UWHAM combines the two work distributions.

    The runtime ATMForce displacement sign reversal for sys_dminus is
    handled by the structprep / production launcher (cntl
    ``DISPLACEMENT`` override or ``OMMSystemABFE.set_displacement``
    monkeypatch); this builder is concerned ONLY with the system XML
    (Modeller positions + water + bbox).

    Returns the build metadata dict (atoms, displacement, ff stack,
    direction, apply_modeller_pre_displacement flag).
    """
    if direction not in ("dplus", "dminus"):
        raise ValueError(
            f"direction must be 'dplus' or 'dminus', got {direction!r}"
        )
    # (b+) corrected spec: BOTH directions pass +displacement_nm
    # to vendored upstream. The differentiator is the Modeller pre-
    # displacement flag (dplus=OFF, dminus=ON). The ATMForce per-particle
    # displacement sign reversal for dminus walker is a runtime concern
    # (set in cntl by structprep launcher), not a build-time concern.
    apply_modeller_pre_disp = (direction == "dminus")
    effective_displacement_nm = tuple(float(d) for d in displacement_nm)

    # MTR ncAA hybrid XML is appended to amber14SB via the multi-XML monkey-
    # patch (``_patched_make_system`` -> ``_make_multi_xml_forcefield_proxy``).
    proteinforcefield = f"{protein_ff} {mtr_xml}"
    displ_A = [float(d) * 10.0 for d in effective_displacement_nm]  # nm -> A

    _patched_make_system(
        receptorfile=receptor_pdb,
        displacement=displ_A,
        xmloutfile=xml_out,
        pdboutfile=pdb_out,
        lig1file=binder_pdb,
        proteinforcefield=proteinforcefield,
        solventforcefield=solvent_ff,
        ligandforcefield="openff-2.0.0",  # unused (binder is PDB, not SDF)
        ionicstrength=ionic_strength_M,
        hmass=hmass_amu,
        apply_modeller_pre_displacement=apply_modeller_pre_disp,
    )

    return {
        "xml_path": xml_out,
        "pdb_path": pdb_out,
        "displacement_nm": list(effective_displacement_nm),
        "displacement_A": displ_A,
        "direction": direction,
        "apply_modeller_pre_displacement": apply_modeller_pre_disp,
        "binder_physical_state": (
            "bound (binder@x0)" if direction == "dplus"
            else "dissociated (binder@x0+d_vec)"
        ),
        "runtime_atmforce_displacement_sign": (
            "+1 (cntl default)" if direction == "dplus"
            else "-1 (cntl override required at structprep)"
        ),
        "walker_direction": "+1 (both directions, (b+) corrected)",
        "proteinforcefield": proteinforcefield,
        "solventforcefield": solvent_ff,
        "ionic_strength_M": ionic_strength_M,
        "leg": "bound",
    }


def build_free_leg_via_upstream(
    binder_pdb: str,
    xml_out: str,
    pdb_out: str,
    protein_ff: str = PROTEIN_FF_DEFAULT,
    solvent_ff: str = SOLVENT_FF_DEFAULT,
    mtr_xml: str = HYBRID_MTR_XML,
    ionic_strength_M: float = 0.15,
    hmass_amu: float = 1.0,
) -> Dict[str, Any]:
    """Build the free-leg ABFE system (binder only in water).

    Free leg has no receptor and no displacement — the binder IS the
    "receptor" from upstream's perspective. We pass ``receptorinFile=binder``
    and omit ``LIG1inFile``. Upstream then runs addSolvent around the
    binder alone and produces a free-peptide-in-water system.
    """
    proteinforcefield = f"{protein_ff} {mtr_xml}"

    _patched_make_system(
        receptorfile=binder_pdb,
        displacement=[0.0, 0.0, 0.0],     # no displacement on free leg
        xmloutfile=xml_out,
        pdboutfile=pdb_out,
        lig1file=None,                    # no second ligand
        proteinforcefield=proteinforcefield,
        solventforcefield=solvent_ff,
        ligandforcefield="openff-2.0.0",
        ionicstrength=ionic_strength_M,
        hmass=hmass_amu,
    )

    return {
        "xml_path": xml_out,
        "pdb_path": pdb_out,
        "displacement_nm": [0.0, 0.0, 0.0],
        "proteinforcefield": proteinforcefield,
        "solventforcefield": solvent_ff,
        "ionic_strength_M": ionic_strength_M,
        "leg": "free",
    }


# ---------------------------------------------------------------------------
# 4-system orchestrator (cp4/bound, cp4/free, wt/bound, wt/free)
# ---------------------------------------------------------------------------
def build_all_four_systems(
    out_root: str,
    seed_tag: str = "s7",
    displacement_nm: Tuple[float, float, float] = DISPLACEMENT_NM_DEFAULT,
    binder_chain: str = "B",
    receptor_chain: str = "A",
    bound_directions: Tuple[str, ...] = ("dplus",),
) -> Dict[str, Any]:
    """Build cp4/bound, cp4/free, wt/bound, wt/free system XML + topology PDB.

    Layout::

        <out_root>/_inputs_split/
            cp4_receptor.pdb       (chain A from 2QKI_Cp4)
            cp4_binder.pdb         (chain B, MTR-4)
            wt_receptor.pdb        (chain A from 2QKI_WT)
            wt_binder.pdb          (chain B, TRP-4)
            cp4_binder_protonated.pdb
            wt_binder_protonated.pdb
            cp4_receptor_protonated.pdb
            wt_receptor_protonated.pdb
        <out_root>/cp4/bound/trackb.pdb + trackb_sys.xml (legacy compat,
            symlinks to dplus when bound_directions includes dplus)
        <out_root>/cp4/bound/trackb_sys_dplus.xml + trackb_dplus.pdb
            (per-direction artifacts when "dplus" in bound_directions)
        <out_root>/cp4/bound/trackb_sys_dminus.xml + trackb_dminus.pdb
            (per-direction artifacts when "dminus" in bound_directions)
        <out_root>/cp4/free/trackb.{pdb,_sys.xml}  (single, no displacement)
        <out_root>/wt/bound/trackb.{pdb,_sys.xml}
        <out_root>/wt/free/trackb.{pdb,_sys.xml}

    ``bound_directions`` (default ``("dplus",)``): which per-direction
    bound-leg systems to build. Pass ``("dplus", "dminus")`` for the
    bidirectional fix per [[trackb_per_direction_system_xml_rebuild_20260601]].
    ``"dminus"`` builds a binder-pre-displaced-by-NEGATIVE-displacement
    system so r11..r21 walkers (Direction=-1) start clash-free.

    Free legs ignore ``bound_directions`` (no displacement applied — the
    free-leg cycle has no binding site, so direction is moot).

    Backward compat: when ``bound_directions == ("dplus",)`` the legacy
    ``trackb_sys.xml`` + ``trackb.pdb`` file names ARE preserved (as the
    primary write target, not symlinks), so existing consumers that look
    up ``trackb_sys.xml`` keep working unchanged. When the call includes
    ``"dminus"``, ``trackb_sys.xml`` becomes a symlink to
    ``trackb_sys_dplus.xml`` (if dplus also built) or absent (if only
    dminus built — callers must use the per-direction names).

    Returns a nested dict mapping (endpoint, leg) → build metadata. For
    bound legs with multiple directions, metadata contains a
    ``per_direction`` sub-dict keyed by direction tag.
    """
    from atm_trackB_setup import resolve_leg_inputs

    inputs = resolve_leg_inputs(seed_tag)
    cp4_complex = inputs["bound"]["cp4"]
    wt_complex = inputs["bound"]["wt"]
    hydrogens_xml = inputs.get("hydrogens_xml")

    if not os.path.isfile(cp4_complex):
        raise FileNotFoundError(f"Cp4 complex PDB missing: {cp4_complex}")
    if not os.path.isfile(wt_complex):
        raise FileNotFoundError(f"WT complex PDB missing: {wt_complex}")

    inputs_dir = os.path.join(out_root, "_inputs_split")
    os.makedirs(inputs_dir, exist_ok=True)

    # Split complexes → receptor + binder PDBs (raw, no protonation yet).
    cp4_rec_raw = os.path.join(inputs_dir, "cp4_receptor_raw.pdb")
    cp4_bin_raw = os.path.join(inputs_dir, "cp4_binder_raw.pdb")
    wt_rec_raw = os.path.join(inputs_dir, "wt_receptor_raw.pdb")
    wt_bin_raw = os.path.join(inputs_dir, "wt_binder_raw.pdb")
    split_complex_to_receptor_binder(
        cp4_complex, cp4_rec_raw, cp4_bin_raw,
        receptor_chain=receptor_chain, binder_chain=binder_chain,
    )
    split_complex_to_receptor_binder(
        wt_complex, wt_rec_raw, wt_bin_raw,
        receptor_chain=receptor_chain, binder_chain=binder_chain,
    )

    # Protonate the receptor + binder PDBs.
    cp4_rec = os.path.join(inputs_dir, "cp4_receptor.pdb")
    wt_rec = os.path.join(inputs_dir, "wt_receptor.pdb")
    cp4_bin = os.path.join(inputs_dir, "cp4_binder.pdb")
    wt_bin = os.path.join(inputs_dir, "wt_binder.pdb")

    rec_meta_cp4 = protonate_receptor(cp4_rec_raw, cp4_rec)
    rec_meta_wt = protonate_receptor(wt_rec_raw, wt_rec)
    bin_meta_cp4 = protonate_binder_with_mtr_xml(
        cp4_bin_raw, cp4_bin,
        mtr_xml=HYBRID_MTR_XML,
        hydrogens_xml=hydrogens_xml,
        binder_chain=binder_chain,
    )
    bin_meta_wt = protonate_binder_with_mtr_xml(
        wt_bin_raw, wt_bin,
        mtr_xml=HYBRID_MTR_XML,
        hydrogens_xml=hydrogens_xml,
        binder_chain=binder_chain,
    )

    # Validate bound_directions input (fail fast before any build).
    valid_dirs = ("dplus", "dminus")
    for d in bound_directions:
        if d not in valid_dirs:
            raise ValueError(
                f"bound_directions entry must be one of {valid_dirs}, "
                f"got {d!r}"
            )
    if not bound_directions:
        raise ValueError(
            "bound_directions must be non-empty (at minimum ('dplus',))"
        )

    results: Dict[Tuple[str, str], Dict[str, Any]] = {}
    legacy_compat_single_dplus = (tuple(bound_directions) == ("dplus",))

    for endpoint, rec, bin_, rec_meta, bin_meta in (
        ("cp4", cp4_rec, cp4_bin, rec_meta_cp4, bin_meta_cp4),
        ("wt",  wt_rec,  wt_bin,  rec_meta_wt,  bin_meta_wt),
    ):
        # Bound leg — per-direction systems.
        bound_dir = os.path.join(out_root, endpoint, "bound")
        os.makedirs(bound_dir, exist_ok=True)

        per_dir_meta: Dict[str, Dict[str, Any]] = {}
        for direction_tag in bound_directions:
            if legacy_compat_single_dplus:
                # Preserve legacy filenames so existing consumers continue
                # to find trackb_sys.xml + trackb.pdb in place.
                d_xml = os.path.join(bound_dir, "trackb_sys.xml")
                d_pdb = os.path.join(bound_dir, "trackb.pdb")
            else:
                d_xml = os.path.join(
                    bound_dir, f"trackb_sys_{direction_tag}.xml"
                )
                d_pdb = os.path.join(
                    bound_dir, f"trackb_{direction_tag}.pdb"
                )
            d_meta = build_bound_leg_via_upstream(
                receptor_pdb=rec,
                binder_pdb=bin_,
                xml_out=d_xml,
                pdb_out=d_pdb,
                displacement_nm=displacement_nm,
                direction=direction_tag,
            )
            d_meta.update({
                "receptor_n_atoms": rec_meta["n_atoms"],
                "binder_n_atoms": bin_meta["n_atoms"],
                "disulfide": bin_meta["disulfide"],
            })
            per_dir_meta[direction_tag] = d_meta

        # When the bidirectional build is requested, also publish a
        # legacy symlink ``trackb_sys.xml`` -> ``trackb_sys_dplus.xml``
        # (if dplus was built) so legacy consumers (e.g. v2.1 launcher
        # cached-sidecar check) still resolve a file. We do NOT symlink
        # the .pdb (the topology PDB content is identical for dplus
        # and dminus — only the Modeller pre-displacement differs, the
        # atom list / bonds / residue ids are the same — so we symlink
        # trackb.pdb -> trackb_dplus.pdb when both built).
        if (not legacy_compat_single_dplus) and "dplus" in per_dir_meta:
            legacy_sys = os.path.join(bound_dir, "trackb_sys.xml")
            legacy_pdb = os.path.join(bound_dir, "trackb.pdb")
            dplus_sys_basename = "trackb_sys_dplus.xml"
            dplus_pdb_basename = "trackb_dplus.pdb"
            for target_basename, link_path in (
                (dplus_sys_basename, legacy_sys),
                (dplus_pdb_basename, legacy_pdb),
            ):
                try:
                    if os.path.islink(link_path) or os.path.exists(link_path):
                        os.remove(link_path)
                    os.symlink(target_basename, link_path)
                except OSError:
                    # Symlink unsupported (rare on POSIX). Fall back to
                    # copying the file so downstream still sees it.
                    src = os.path.join(bound_dir, target_basename)
                    if os.path.isfile(src):
                        shutil.copy2(src, link_path)

        # Bound-leg metadata in result: legacy single-dplus mode keeps the
        # flat schema (xml_path / pdb_path) for back-compat; multi-direction
        # mode adds per_direction sub-dict.
        if legacy_compat_single_dplus:
            results[(endpoint, "bound")] = per_dir_meta["dplus"]
        else:
            primary_tag = (
                "dplus" if "dplus" in per_dir_meta else next(iter(per_dir_meta))
            )
            results[(endpoint, "bound")] = {
                "leg": "bound",
                "per_direction": per_dir_meta,
                "directions_built": list(per_dir_meta.keys()),
                "primary_direction": primary_tag,
                "receptor_n_atoms": rec_meta["n_atoms"],
                "binder_n_atoms": bin_meta["n_atoms"],
                "disulfide": bin_meta["disulfide"],
                # Surface representative paths so JSON consumers can still
                # find a canonical xml/pdb path (primary direction).
                "xml_path": per_dir_meta[primary_tag]["xml_path"],
                "pdb_path": per_dir_meta[primary_tag]["pdb_path"],
            }

        # Free leg (single, no displacement)
        free_dir = os.path.join(out_root, endpoint, "free")
        os.makedirs(free_dir, exist_ok=True)
        free_xml = os.path.join(free_dir, "trackb_sys.xml")
        free_pdb = os.path.join(free_dir, "trackb.pdb")
        free_meta = build_free_leg_via_upstream(
            binder_pdb=bin_,
            xml_out=free_xml,
            pdb_out=free_pdb,
        )
        free_meta.update({
            "binder_n_atoms": bin_meta["n_atoms"],
            "disulfide": bin_meta["disulfide"],
        })
        results[(endpoint, "free")] = free_meta

    return {
        "out_root": out_root,
        "seed_tag": seed_tag,
        "inputs_dir": inputs_dir,
        "complex_paths": {"cp4": cp4_complex, "wt": wt_complex},
        "hybrid_mtr_xml": HYBRID_MTR_XML,
        "hydrogens_xml": hydrogens_xml,
        "displacement_nm": list(displacement_nm),
        # JSON-friendly results (tuple keys serialize as "endpoint/leg")
        "systems": {f"{e}/{l}": meta for (e, l), meta in results.items()},
    }


# ---------------------------------------------------------------------------
# Main (smoke / production entry)
# ---------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Track B v2.1 ABFE system builder — wraps upstream "
            "atom_openmm.make_atm_system_from_rcpt_lig."
        )
    )
    p.add_argument("--out-root", required=True,
                   help="output root (e.g. outputs/_trackb/production_v2_1)")
    p.add_argument("--seed-tag", default="s7",
                   help="prepared seed (outputs/2QKI_*_calib_<seed>/)")
    p.add_argument("--displacement-nm", default="2.5,0.0,0.0",
                   help="ABFE displacement vector (nm)")
    p.add_argument(
        "--bound-directions",
        default="dplus",
        help=(
            "Comma-separated bound-leg directions to build. Use "
            "'dplus,dminus' for the bidirectional fix per "
            "[[trackb_per_direction_system_xml_rebuild_20260601]] — "
            "produces trackb_sys_dplus.xml + trackb_sys_dminus.xml "
            "(binder physically pre-displaced by +/-displacement in each, "
            "so r0..r10 + r11..r21 walkers are both clash-free at start). "
            "Default 'dplus' preserves the legacy single-direction "
            "trackb_sys.xml filename for back-compat."
        ),
    )
    args = p.parse_args()

    out_root = os.path.abspath(args.out_root)
    os.makedirs(out_root, exist_ok=True)

    displ_parts = [float(x) for x in args.displacement_nm.split(",")]
    if len(displ_parts) != 3:
        print(f"ERROR: --displacement-nm needs 3 floats, got {displ_parts}",
              file=sys.stderr)
        return 2
    displacement_nm = tuple(displ_parts)

    bound_directions = tuple(
        d.strip() for d in args.bound_directions.split(",") if d.strip()
    )
    if not bound_directions:
        print("ERROR: --bound-directions must be non-empty", file=sys.stderr)
        return 2

    t0 = time.time()
    meta = build_all_four_systems(
        out_root=out_root,
        seed_tag=args.seed_tag,
        displacement_nm=displacement_nm,
        bound_directions=bound_directions,
    )
    t1 = time.time()
    meta["wall_time_s"] = t1 - t0
    meta["regime"] = "ranking_only"
    meta["builder"] = "scripts/phase4_trackB_v2_make_system.py"
    meta["bound_directions"] = list(bound_directions)
    meta["method_ref"] = "upstream make_system bound-leg shadow setup, 2026-05-31"
    meta["started_at"] = time.strftime("%Y-%m-%dT%H:%M:%S",
                                       time.gmtime(t0))
    meta["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S",
                                        time.gmtime(t1))

    meta_path = os.path.join(out_root, "build_metadata.json")
    if os.path.isfile(meta_path):
        ts = time.strftime("%Y%m%dT%H%M%S")
        meta_path = os.path.join(out_root, f"build_metadata_{ts}.json")
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"\nBuild metadata: {meta_path}")
    print(f"Wall time:      {t1 - t0:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
