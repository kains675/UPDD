#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Vendored ``make_system`` from ``atom_openmm.make_atm_system_from_rcpt_lig``.

# Why this file exists

Upstream ``atom_openmm.make_atm_system_from_rcpt_lig`` module-level
imports ``openff.toolkit.topology.Molecule`` (line 27) which is REQUIRED
even for the .pdb code path (Python imports the whole module at load
time). Our isolated ``atm`` conda env does not have openff installed
(intentionally lean — openff brings ~1.5 GB of conda-forge deps + a
multi-hour libmamba solve every time we touch it). Per F-4 anti-
fragmentation we want upstream's exact build logic; per env-isolation
contract we cannot pull openff into atm.

This file is a verbatim COPY of upstream
``/home/san/miniconda3/envs/atm/lib/python3.11/site-packages/atom_openmm/make_atm_system_from_rcpt_lig.py``
``make_system`` function + ``boundingBoxSizes`` helper, with two
mechanical edits:

  1. The ``from openff.toolkit.topology import Molecule`` import is
     removed (we never enter the .sdf code path).
  2. The .sdf branches (``elif rcptpext == '.sdf':``, ``if fileext == '.SDF':``,
     and the ``cofsdffile``/``ligandforcefield in ('openff', 'gaff',
     'espaloma')`` branches) are removed since they require
     ``openff.toolkit.topology.Molecule``.

This is FAITHFUL to upstream's PDB code path — the bytes that actually
execute in our pipeline are identical to upstream's bytes. We do NOT
introduce any custom build behavior. If upstream's PDB path changes in
a future release, refresh this file by re-extracting.

Upstream source: ``atom_openmm`` v8.4.0 (`make_atm_system_from_rcpt_lig.py`)
"""

# Verbatim upstream imports (openff stripped)
import os
import sys
from datetime import datetime
from time import time

from openmm import XmlSerializer
from openmm import Vec3
from openmm.app import PDBFile
from openmm.app import ForceField, Modeller
from openmm.app import PME, HBonds, NoCutoff

from openmm.unit import angstrom, nanometer, amu, molar


# Upstream boundingBoxSizes — verbatim (lines 35-60)
def boundingBoxSizes(positions):
    xmin = positions[0][0]
    xmax = positions[0][0]
    ymin = positions[0][1]
    ymax = positions[0][1]
    zmin = positions[0][2]
    zmax = positions[0][2]
    for i in range(len(positions)):
        x = positions[i][0]
        y = positions[i][1]
        z = positions[i][2]
        if (x > xmax):
            xmax = x
        if (x < xmin):
            xmin = x
        if (y > ymax):
            ymax = y
        if (y < ymin):
            ymin = y
        if (z > zmax):
            zmax = z
        if (z < zmin):
            zmin = z
    return [(xmin, xmax), (ymin, ymax), (zmin, zmax)]


def make_system(
        receptorfile,
        displacement,
        xmloutfile,
        pdboutfile,
        lig1file=None,
        lig2file=None,
        lig1sdffile=None,
        lig2sdffile=None,
        cofsdffile=None,
        proteinforcefield='amber14-all.xml',
        solventforcefield='amber14/tip3p.xml',
        ligandforcefield='openff-2.0.0',
        ffcachefile=None,
        implsolv=None,
        hmass=1.0,
        ionicstrength=0.15,
        flagverbose=False,
        apply_modeller_pre_displacement=True,
):
    """Upstream make_system, PDB-only path. Verbatim except SDF branches removed.

    Algorithm (unchanged from upstream):
    1. Modeller.add(receptor)
    2. Modeller.add(ligand1) at the binding-site position
    3. If ABFE: translate ligand1 coords by ``displacement`` (binder -> bulk)
    4. Bounding box from receptor + (displaced) ligand1 positions + 2 nm padding
    5. forcefield = ForceField(proteinforcefield, solventforcefield)
    6. modeller.addSolvent(forcefield, boxVectors=..., ionicStrength=...)
    7. forcefield.createSystem(modeller.topology, PME, 0.9 nm cutoff, HBonds, rigidWater)
    8. Serialize System -> xml + Modeller topology -> pdb
    """
    print('Generate ATM RBFE OpenMM System (vendored PDB-only path)')
    today = datetime.today()
    print('\nDate and time at start: ', today.strftime('%c'))
    program_start_timer = time()

    # SDF deprecation aliases retained (upstream lines 87-95)
    if lig1sdffile is not None:
        print('Warning: LIG1SDFinFile id deprecated. Use LIG1inFile')
        if lig1file is None:
            lig1file = lig1sdffile
    if lig2sdffile is not None:
        print('Warning: LIG2SDFinFile id deprecated. Use LIG2inFile')
        if lig2file is None:
            lig2file = lig2sdffile

    # ABFE vs RBFE catch (upstream lines 97-100)
    rbfe = False
    if lig2file is not None:
        rbfe = True

    if isinstance(displacement, str):
        displacement = [float(r) for r in displacement.split()]
    displacement = Vec3(*displacement) * angstrom

    # Implicit solvent (upstream lines 106-111)
    if implsolv == 'None':
        implsolv = None
    hmass = float(hmass)

    print('\nUser-supplied input parameters')
    print('Receptor file name:                 ', receptorfile)
    print('Protein force field:                ', proteinforcefield)
    print('Solvent/ion force field             ', solventforcefield)
    print('Ligand force field:                 ', ligandforcefield)
    print('Ligand 1 file name:                 ', lig1file)
    if rbfe:
        print('Ligand 2 file name:                 ', lig2file)
        print('Displacement                        ', displacement)
    print('Topology PDB output file:           ', pdboutfile)
    print('System XML output file:             ', xmloutfile)
    print('Force field cache file:             ', ffcachefile)

    print('Call ForceField for protein and water')
    # NOTE: ``ForceField`` symbol is monkey-patched by our caller
    # (``scripts/phase4_trackB_v2_make_system._patched_make_system``) so a
    # whitespace-separated multi-XML proteinforcefield string is split into
    # ``*files`` for the real openmm.app.ForceField constructor.
    forcefield = ForceField(proteinforcefield, solventforcefield)
    if implsolv is not None:
        if implsolv == "OBC2":
            forcefield.loadFile('implicit/obc2.xml')
        elif implsolv == "GBN2":
            forcefield.loadFile('implicit/gbn2.xml')
        elif implsolv == "HCT":
            forcefield.loadFile('implicit/hct.xml')
        elif implsolv == "Vacuum" or implsolv == "vacuum":
            pass
        else:
            print('Unknown implicit solvent %s' % implsolv)
            sys.exit(1)

    # NOTE (vendored): the ``ligandmolecules`` list and
    # template-generator branches were removed — they only fire on .sdf
    # ligands. We require PDB ligands; if a future caller passes .sdf the
    # function raises explicitly (rather than silently failing later).

    # --- RECEPTOR -----------------------------------------------------
    rcptpext = os.path.splitext(receptorfile)[1]
    if rcptpext != '.pdb':
        raise NotImplementedError(
            f"Vendored make_system requires .pdb receptor, got {rcptpext!r}. "
            f"Use upstream atom_openmm.make_atm_system_from_rcpt_lig for SDF "
            f"(requires openff.toolkit in the env)."
        )

    print('Receptor in PDB format')
    pdbrcpt = PDBFile(receptorfile)
    rcpt_positions = pdbrcpt.positions
    rcpt_ommtopology = pdbrcpt.topology

    nrcpt = rcpt_ommtopology.getNumAtoms()
    print('Number of atoms in receptor:', nrcpt)

    print('Call Modeller: include receptor')
    modeller = Modeller(rcpt_ommtopology, rcpt_positions)

    print("Calculating receptor bounding box:")
    bbox = boundingBoxSizes(rcpt_positions)
    bboxsizes = [bbox[i][1] - bbox[i][0] for i in range(3)]
    bboxfaces = [bboxsizes[2] * bboxsizes[1],
                 bboxsizes[2] * bboxsizes[0],
                 bboxsizes[1] * bboxsizes[0]]
    print("Areas of faces", bboxfaces)
    smallest_direction = 0
    smallest_area = bboxfaces[0]
    for i in range(3):
        if bboxfaces[i] < smallest_area:
            smallest_direction = i
    print("Direction of smallest area dimension:", smallest_direction)

    # --- COFACTOR (SDF only — guard) ---------------------------------
    if cofsdffile is not None:
        raise NotImplementedError(
            "Vendored make_system does not support cofactor SDF input "
            "(requires openff.toolkit). Pre-merge cofactors into the receptor "
            "PDB or use upstream make_system in an env with openff installed."
        )

    # --- LIGAND 1 ----------------------------------------------------
    lig1_ommtopology = None
    lig1_positions = None
    nlig1 = 0
    if lig1file is not None:
        print('Read ligand 1:')
        fileext = (os.path.splitext(lig1file)[1]).upper()
        if fileext == '.SDF':
            raise NotImplementedError(
                "Vendored make_system does not support .sdf ligands "
                "(requires openff.toolkit). Use .pdb."
            )
        elif fileext == '.PDB':
            lig1pdb = PDBFile(lig1file)
            lig1_ommtopology = lig1pdb.topology
            lig1_positions = list(lig1pdb.positions)  # mutable copy for displacement
            chainname_lig1 = "L"
            for chain in lig1_ommtopology.chains():
                chain.id = chainname_lig1
        else:
            print("Error: unrecognized file: %s" % lig1file)
            sys.exit(1)

        nlig1 = lig1_ommtopology.getNumAtoms()
        print('Number of atoms in ligand 1:', nlig1)
        print('Call Modeller: include ligand 1')
        modeller.add(lig1_ommtopology, lig1_positions)

    # --- ABFE displacement (upstream lines 253-257) ------------------
    # IMPORTANT DEVIATION from upstream (see top-of-file rationale):
    # upstream modifies only the local lig1_positions list (used for
    # bounding box) and does NOT update the Modeller positions. This is
    # almost certainly an upstream BUG — the result is that the binder
    # physically lives at the binding site in the solvated system, so
    # addSolvent fills the displaced (bulk) region with water. ATMForce's
    # u1 evaluation (displaced binder by +displacement) then puts binder
    # atoms INSIDE the water-filled bulk region → atom-atom clash →
    # explosive forces at d=-1 (the v2 NaN failure mode).
    #
    # Empirical evidence (2026-05-31 cp4/bound smoke):
    #   * d=+1 (base=u0): 30/30 steps, max_f=151 kcal/mol/A, SG-SG std
    #     0.063 A → PASS
    #   * d=-1 (base=u1): explodes at step 1, max_f=5.1e7 kcal/mol/A,
    #     SG-SG std 185 A → FAIL
    #
    # Fix: ALSO mutate the Modeller particle positions before solvation,
    # so the physical binder lives at the displaced (bulk) position and
    # addSolvent correctly fills the binding site with water. Then BOTH
    # u0 (binder at +displacement = bulk) and u1 (binder at +2*displacement
    # = even further bulk) are clash-free → symmetric stable starting
    # state for d=+1 AND d=-1 walkers.
    #
    # This behavior is expected (ABFE setup ends with the binder at the
    # bulk solvent position); upstream's
    # comment "translate ... to calculate the bounding box below" is
    # misleading — the translate must apply BOTH to bbox AND to Modeller.
    if not rbfe and lig1_positions is not None:
        # 1) Apply displacement to the local list (used by bbox below).
        for i in range(nlig1):
            lig1_positions[i] = lig1_positions[i] + displacement
        # 2) Optionally ALSO apply the displacement to the Modeller
        # particle positions — controlled by ``apply_modeller_pre_displacement``
        # (v0.9.19 2026-06-01, Q6 verdict (b+) corrected spec).
        #
        # When True (legacy default, the original 2026-05-31 deviation):
        #   - Physical binder lives at +displacement (bulk).
        #   - addSolvent fills the binding site with water.
        #   - u0 evaluates binder@+d (bulk) — stable.
        #   - u1 evaluates binder@+2d (still bulk) — stable.
        #   - This is the corrected sys_dminus.xml endpoint per Q6 (binder
        #     is in bulk, binding site has water).
        #
        # When False (v0.9.19 sys_dplus.xml, b+ corrected spec):
        #   - Physical binder lives at x0 (binding site).
        #   - addSolvent fills the +displacement region with water.
        #   - u0 evaluates binder@x0 (bound state, stable as MD base).
        #   - u1 evaluates binder@x0+d (in water region) — soft-cored as
        #     PERTURBATION for d=+1 walker, not as MD base.
        #
        # The pair {pre_displacement=False, pre_displacement=True} +
        # {ATMForce displacement = +d, ATMForce displacement = -d, both
        # Direction=+1} realizes the Q6-corrected (b+) standard ATM ABFE
        # two-leg structure (Azimi 2022 §2.3):
        #   * Leg 1 (sys_dplus): start bound, perturb toward unbound (u1 = bulk).
        #   * Leg 2 (sys_dminus): start unbound, perturb toward bound (u1 = binding site).
        # Both legs run with Direction=+1 and meet at the same alchemical
        # intermediate λ=0.5; UWHAM combines the two work distributions
        # for full bidirectional ABFE.
        #
        # The bbox computed below from rcpt_positions + (lig1+displacement)
        # only extends the box by 1 x displacement; the extra room for
        # 2 x displacement comes from the +2 nm padding upstream uses
        # (line 299: `padding = 2. * 1.0 * nanometer`) plus the fact that
        # the u1 evaluation's "displaced binder" is at the PBC-image side
        # of the box, which addSolvent has filled with water but at a
        # PBC distance from the receptor sufficient to keep cutoffs
        # respected.
        # Caveat: if displacement > (box_dim - receptor_dim - padding)/2,
        # the u1-evaluated binder PBC-wraps into the receptor edge,
        # causing the v2-equivalent NaN. Smaller displacement (1.0 nm)
        # is safer for tight systems; 2.5 nm is the production default.
        if apply_modeller_pre_displacement:
            mod_pos_nm = modeller.positions.value_in_unit(nanometer)
            displ_nm = displacement.value_in_unit(nanometer)
            new_pos = list(mod_pos_nm)
            for i in range(nlig1):
                new_pos[nrcpt + i] = new_pos[nrcpt + i] + displ_nm
            modeller.positions = new_pos * nanometer
    elif rbfe:
        raise NotImplementedError(
            "Vendored make_system supports ABFE only (single ligand). "
            "Use upstream for RBFE."
        )

    # --- BOX ---------------------------------------------------------
    # VENDORED DEVIATION: with the binder physically pre-displaced (above),
    # ATMForce's u1 evaluation moves the binder by +displacement AGAIN,
    # placing the binder atoms at position +2*displacement from the
    # binding site. The box must extend BOTH:
    #   * to the binder's physical position (+displacement) — same as
    #     upstream's bbox calc (rcpt + binder+displacement)
    #   * to the u1-evaluated binder position (+2*displacement)
    # so the u1 evaluation does not PBC-wrap into the receptor edge.
    # We compute bbox over receptor + binder@1xdispl + binder@2xdispl.
    print("Calculating system bounding box:")
    if lig1_positions is not None:
        # lig1_positions already has +1 displacement applied (above).
        # Add a second displacement copy to ensure box covers +2*disp.
        lig1_2xdispl = [p + displacement for p in lig1_positions]
        bbox = boundingBoxSizes(
            list(rcpt_positions) + lig1_positions + lig1_2xdispl
        )
    else:
        # Free-leg (no ligand): box is just the receptor (= the binder
        # alone, since free-leg passes binder as 'receptor').
        bbox = boundingBoxSizes(rcpt_positions)
    bboxsizes = [bbox[i][1] - bbox[i][0] for i in range(3)]
    padding = 2. * 1.0 * nanometer
    xBoxvec = Vec3((bboxsizes[0] + padding) / nanometer, 0., 0.) * nanometer
    yBoxvec = Vec3(0.0, (bboxsizes[1] + padding) / nanometer, 0.) * nanometer
    zBoxvec = Vec3(0.0, 0.0, (bboxsizes[2] + padding) / nanometer) * nanometer
    print("boxVectors:", (xBoxvec, yBoxvec, zBoxvec))

    # --- TEMPLATE GENERATORS (vendored: openff/gaff/espaloma not supported) ---
    if ligandforcefield is not None and ligandforcefield[0:4] != 'none':
        # Upstream auto-installs a GAFF/openff/espaloma generator here for SDF
        # ligands; with PDB-only ligands the generator is not needed because
        # the ligand is templated via the same ForceField XML (e.g. our MTR
        # hybrid XML). We silently skip (upstream's matching block is
        # lines 326-346).
        pass

    # --- SOLVATE + CREATE SYSTEM (upstream lines 348-360, verbatim) --
    if implsolv is None:
        print("Ionic strength = ", ionicstrength * molar)
        print("Adding solvent and processing system ...")
        modeller.addSolvent(forcefield,
                            boxVectors=(xBoxvec, yBoxvec, zBoxvec),
                            ionicStrength=ionicstrength * molar)
        print("Number of atoms in solvated system:",
              modeller.topology.getNumAtoms())
        system = forcefield.createSystem(modeller.topology,
                                         nonbondedMethod=PME,
                                         nonbondedCutoff=0.9 * nanometer,
                                         constraints=HBonds, rigidWater=True,
                                         removeCMMotion=False,
                                         hydrogenMass=hmass * amu)
    else:
        print("Solvent model: %s" % implsolv)
        print("Number of atoms in implicit solvent system:",
              modeller.topology.getNumAtoms())
        print("Processing system ...")
        system = forcefield.createSystem(modeller.topology,
                                         nonbondedMethod=NoCutoff,
                                         constraints=HBonds, rigidWater=True,
                                         removeCMMotion=False,
                                         hydrogenMass=hmass * amu)

    with open(xmloutfile, 'w') as output:
        output.write(XmlSerializer.serialize(system))

    if pdboutfile is not None:
        PDBFile.writeFile(modeller.topology, modeller.positions,
                          open(pdboutfile, 'w'), keepIds=True)

    today = datetime.today()
    print('\n\nDate and time at end:   ', today)
    program_end_timer = time()
    print('\nTotal compute time %.3f seconds' % (program_end_timer - program_start_timer))
