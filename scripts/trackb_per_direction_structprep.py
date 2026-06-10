#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B per-direction structprep — Option B (2026-05-31).

Option B is the PRIMARY architectural fix for the v2.1 d=-1 NaN failure mode:
upstream ``abfe_structprep`` hardcodes ``direction=+1`` in ``do_equil``
(``atom_openmm/abfe_structprep.py:340``), so the canonical
``trackb_0.xml`` is an equilibrium at lambda=0.5 / direction=+1 — perfectly
suited to forward-leg replicas (r0..r10) but mis-equilibrated for the
backward-leg replicas (r11..r21) which start at direction=-1.

The Option-B fix: run structprep TWICE per leg, once with direction=+1
(unmodified upstream path → produces ``trackb_0_dplus.xml``) and once with
direction=-1 (monkeypatched ``do_lambda_annealing`` + ``do_equil`` →
produces ``trackb_0_dminus.xml``). Each direction equilibrates to its own
ground truth at lambda=0.5. Production then injects ``_dplus.xml`` into
r0..r10 and ``_dminus.xml`` into r11..r21 via per-replica
``trackb_ckpt.xml`` copies (consumed by ``ommreplica.load_checkpoint`` at
``ommreplica.py:79``).

# PREP-ONLY contract

This script:
  1. PREP: monkeypatched per-direction structprep on already-built systems
     (``trackb.pdb`` + ``trackb_sys.xml`` + ``trackb_asyncre.cntl`` from
     v2.1 ``setup_one_leg_v21``).
  2. SANITY: per-direction atom-count + bond-count + cyclic_ss SG-SG dist +
     per-residue ``|Sigma q|`` audit.
  3. WRITE per-direction starting XMLs ``trackb_0_dplus.xml`` /
     ``trackb_0_dminus.xml`` for dry-run input.
  4. NEVER touches production. NEVER launches anything that writes
     into ``r0..r21/``. Free-leg PID 1876426 + V100 Track A
     untouched.

# Hardware contract (v0.9.10 device-guard fix, 2026-05-31)

* **Actual hardware** (no "device 1" GPU exists anywhere):
  - Host (this box): RTX 5070 Ti 16 GB — single GPU at device 0.
    Currently occupied by Track B free leg PID 1876426.
  - VM (san@192.168.122.155): Tesla V100 32 GB — single GPU at device 0.
    Currently occupied by Track A QM 1-traj batch (util ~88 %).
  - "V100 dual-GPU" in the project narrative = host 5070Ti + VM V100
    treated as a virtual pair; it is NOT two V100 cards.
* Per condition C1: free leg PID 1876426 (5070Ti, host) MUST complete
  before any production launch. The prep script supports two GPU
  hosts, both occupy ``CUDA_VISIBLE_DEVICES=0`` on their respective
  box. Choice via ``--gpu-host``:

  - ``--gpu-host vm`` (default, recommended): ssh dispatch to
    ``san@192.168.122.155`` after V100 Track A QM batch idle
    (auto-checked via read-only ``nvidia-smi --query-gpu=utilization.gpu``;
    util > 30 % → REFUSE).
  - ``--gpu-host local``: run on host 5070Ti after free leg PID
    1876426 has terminated (auto-checked via ``os.kill(pid, 0)``;
    alive → REFUSE).
  - ``--gpu-host cpu``: OpenMM CPU platform (~10× slower, safe for
    audit-only sanity-check; ignores both gates).

* Legacy ``--cuda-device`` option is retained for backward compat but
  ONLY accepts ``0`` or ``cpu``. Any other value (including ``1``)
  raises a clear error — that GPU does not exist on either host. The
  literal device index is forced to ``0`` in both VM and local modes
  because each box has a single GPU.
* Memory: 4 systems × structprep (~6 GB GPU memory per run) — only
  one structprep at a time per script invocation.
* SSH safety: only read-only ``nvidia-smi --query-gpu=...`` queries.
  Never kill, never reset, never touch GUI session.

# Cross-references

* v2.1 launcher (source of v2.1 systems): ``scripts/trackb_production_v2_1_upstream.py``
* Upstream structprep being patched: ``atom_openmm/abfe_structprep.py``
* Upstream ommreplica state-load: ``atom_openmm/ommreplica.py:75-79``
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shlex
import shutil
import subprocess
import sys
import textwrap
import time
from typing import Optional, List, Dict, Tuple, Any

_PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "utils"))
sys.path.insert(0, os.path.join(_PROJ_ROOT, "scripts"))


# ---------------------------------------------------------------------------
# Per-direction system XML selection (v0.9.19 architectural fix, 2026-06-01)
# ---------------------------------------------------------------------------
#
# Per [[trackb_per_direction_system_xml_rebuild_20260601]] the bound-leg
# dminus NaN failure mode is NOT resolvable by direction-patching equil +
# annealing alone — the underlying system XML must also be built with the
# binder pre-displaced by the OPPOSITE displacement vector so that the
# dminus u0/u1 endpoints both land in solvated bulk (instead of inside the
# receptor pocket).
#
# Upstream ``atom_openmm.ommsystem.OMMSystemABFE(basename, ...)`` always
# reads ``<basename>_sys.xml`` regardless of direction. To consume the new
# per-direction system XMLs (``trackb_sys_dplus.xml`` /
# ``trackb_sys_dminus.xml``) we transiently swap the active
# ``trackb_sys.xml`` before invoking ``abfe_structprep`` and restore the
# original afterwards (try / finally), exactly mirroring how the in-process
# monkeypatch of ``do_equil`` is reverted at the end of a run.
#
# Legacy compatibility: if only ``trackb_sys.xml`` exists (single-direction
# build, e.g. legacy v2.1 cohorts predating 2026-06-01), the swap is a
# no-op — the run proceeds against the shared system XML exactly as before.
# ---------------------------------------------------------------------------
def _stage_one_file(
    per_direction_path: str,
    per_direction_basename: str,
    active_path: str,
    backup_basename: str,
    leg_dir: str,
) -> Dict[str, Any]:
    """Helper for the file-level swap (.xml or .pdb).

    Stages ``per_direction_path`` as ``active_path`` via symlink (basename
    target = portable); falls back to copy on filesystems without symlink
    support. Returns dict describing what was staged + how to restore.
    Idempotent when ``active_path`` already resolves to the per-direction
    file.
    """
    entry: Dict[str, Any] = {
        "active_path": active_path,
        "per_direction_path": per_direction_path,
        "swapped": False,
        "backup_path": None,
        "preexisting_symlink_target": None,
    }
    # Idempotency: if active already symlinked to the per-direction file
    # (basename target or absolute equivalent), the swap is a no-op — but
    # we MUST record the existing link target as preexisting so that
    # _restore_one_file recreates it after the run. Without this, restore
    # would just remove the link and leave nothing (observed empirically
    # 2026-06-01 03:39: v0.9.19 wiped the VM-side trackb_sys.xml +
    # trackb.pdb symlinks across the 4 dispatches because the no-op
    # idempotent path left preexisting_symlink_target=None → restore was
    # a destructive no-op).
    if os.path.islink(active_path):
        link_target = os.readlink(active_path)
        if (link_target == per_direction_basename
                or os.path.abspath(link_target)
                == os.path.abspath(per_direction_path)):
            entry["preexisting_symlink_target"] = link_target
            entry["swapped"] = True
            return entry
        entry["preexisting_symlink_target"] = link_target
        os.remove(active_path)
    elif os.path.isfile(active_path):
        backup_path = os.path.join(leg_dir, backup_basename)
        if os.path.exists(backup_path):
            ts = time.strftime("%Y%m%dT%H%M%S")
            backup_path = backup_path + f".{ts}"
        os.rename(active_path, backup_path)
        entry["backup_path"] = backup_path
    try:
        os.symlink(per_direction_basename, active_path)
    except OSError:
        shutil.copy2(per_direction_path, active_path)
    entry["swapped"] = True
    return entry


def _restore_one_file(entry: Dict[str, Any]) -> None:
    """Helper for the file-level restore (.xml or .pdb). Best-effort:
    log restore failures to stderr but never raise."""
    if not entry or not entry.get("swapped"):
        return
    active_path = entry["active_path"]
    backup_path = entry.get("backup_path")
    preexisting_link = entry.get("preexisting_symlink_target")
    try:
        if os.path.islink(active_path) or os.path.isfile(active_path):
            os.remove(active_path)
    except OSError as exc:
        print(
            f"WARN: could not remove staged file {active_path!r}: {exc}",
            file=sys.stderr,
        )
    if preexisting_link is not None:
        try:
            os.symlink(preexisting_link, active_path)
        except OSError as exc:
            print(
                f"WARN: could not restore preexisting symlink "
                f"{active_path!r} -> {preexisting_link!r}: {exc}",
                file=sys.stderr,
            )
    elif backup_path is not None and os.path.exists(backup_path):
        try:
            os.rename(backup_path, active_path)
        except OSError as exc:
            print(
                f"WARN: could not restore backup "
                f"{backup_path!r} -> {active_path!r}: {exc}",
                file=sys.stderr,
            )


def _select_sys_xml_for_direction(
    leg_dir: str,
    jobname: str,
    direction_tag: str,
) -> Dict[str, Any]:
    """Stage BOTH ``<jobname>_sys.xml`` AND ``<jobname>.pdb`` as the
    per-direction variants so upstream ``OMMSystemABFE`` (which loads
    BOTH files) sees a matching atom-count pair.

    v0.9.19.1 bug-fix (2026-06-01 03:30) on v0.9.19: the initial impl
    only swapped the .xml. Per the per-direction system rebuild the
    .pdb topology atom count ALSO differs across directions (addSolvent
    places different water counts based on binder physical position).
    Loading dplus .pdb (92855 atoms) against dminus .xml (92804 atoms)
    raised
        OpenMMException: Called setPositions() on a Context with the
        wrong number of positions
    in ``do_mintherm`` immediately on dispatch start (verified empirically
    cp4/bound dminus 2026-06-01 03:12:45). Fixed by also swapping
    ``<jobname>.pdb`` → ``<jobname>_<tag>.pdb`` for the run.

    Returns dict with:
      * ``swapped`` (bool): True iff BOTH .xml AND .pdb per-direction
        files were found + staged.
      * ``per_direction_path`` (str|None): backwards-compat — the .xml
        path. Existing audit JSON consumers keep working.
      * ``active_sys_path`` (str): canonical .xml path (backwards-compat).
      * ``backup_path`` (str|None): legacy .xml backup path
        (backwards-compat).
      * ``preexisting_symlink_target`` (str|None): legacy .xml preexisting
        symlink target (backwards-compat).
      * ``xml_entry`` (dict): full xml stage info.
      * ``pdb_entry`` (dict): full pdb stage info.

    Cohort-safe: never deletes data; always restorable via
    ``_restore_sys_xml_after_direction``.

    Idempotent: if BOTH active paths already resolve to the requested
    per-direction files, returns ``swapped=True`` without rewriting.

    Fail-fast: if the .xml per-direction file exists but the matching
    .pdb is missing, raises ``RuntimeError`` (would have caused
    atom-count mismatch silent failure otherwise).
    """
    # XML
    xml_per_dir_basename = f"{jobname}_sys_{direction_tag}.xml"
    xml_per_dir_path = os.path.join(leg_dir, xml_per_dir_basename)
    xml_active_path = os.path.join(leg_dir, f"{jobname}_sys.xml")
    xml_backup_basename = f"{jobname}_sys.xml.bak_{direction_tag}"
    # PDB
    pdb_per_dir_basename = f"{jobname}_{direction_tag}.pdb"
    pdb_per_dir_path = os.path.join(leg_dir, pdb_per_dir_basename)
    pdb_active_path = os.path.join(leg_dir, f"{jobname}.pdb")
    pdb_backup_basename = f"{jobname}.pdb.bak_{direction_tag}"

    info: Dict[str, Any] = {
        "swapped": False,
        "per_direction_path": None,
        "active_sys_path": xml_active_path,
        "backup_path": None,
        "preexisting_symlink_target": None,
        "xml_entry": None,
        "pdb_entry": None,
    }

    if not os.path.isfile(xml_per_dir_path):
        # No per-direction XML — legacy single-system flow.
        return info
    if not os.path.isfile(pdb_per_dir_path):
        # XML present but PDB missing — would have caused atom-count
        # mismatch silently. Fail-fast cohort-safe halt.
        raise RuntimeError(
            f"Per-direction XML {xml_per_dir_path!r} exists but matching "
            f"PDB {pdb_per_dir_path!r} is missing. Both must be present "
            f"for the per-direction structprep to work — rebuild via "
            f"scripts/phase4_trackB_v2_make_system.py "
            f"--bound-directions {direction_tag}."
        )

    info["per_direction_path"] = xml_per_dir_path

    xml_entry = _stage_one_file(
        per_direction_path=xml_per_dir_path,
        per_direction_basename=xml_per_dir_basename,
        active_path=xml_active_path,
        backup_basename=xml_backup_basename,
        leg_dir=leg_dir,
    )
    pdb_entry = _stage_one_file(
        per_direction_path=pdb_per_dir_path,
        per_direction_basename=pdb_per_dir_basename,
        active_path=pdb_active_path,
        backup_basename=pdb_backup_basename,
        leg_dir=leg_dir,
    )

    info["xml_entry"] = xml_entry
    info["pdb_entry"] = pdb_entry
    # Surface back-compat keys at top-level (matches v0.9.19 schema).
    info["backup_path"] = xml_entry.get("backup_path")
    info["preexisting_symlink_target"] = xml_entry.get(
        "preexisting_symlink_target"
    )
    info["swapped"] = bool(xml_entry["swapped"] and pdb_entry["swapped"])
    return info


def _make_patched_set_displacement(direction_val: int):
    """Return a OMMSystemABFE.set_displacement clone that negates
    ``self.displ`` for the dminus walker.

    v0.9.19 (b+) corrected spec: the ATMForce per-particle
    displacement vector is per-particle FIXED at addParticle time
    (``ommsystem.py:492``: ``self.atmforce.setParticleParameters(i,
    Vec3(self.displ[0], ...) / nanometer)``). The Direction global
    parameter does NOT flip this vector at runtime. Therefore, to make
    the sys_dminus.xml walker's u1 evaluation reach the bound state, we
    must NEGATE ``self.displ`` BEFORE ``set_atmforce`` runs.

    Upstream ``set_displacement`` (``ommsystem.py:365-371``) simply reads
    ``DISPLACEMENT`` from the cntl keywords (in Angstrom) and stores as
    ``self.displ`` (in Angstrom). The patched clone calls the original
    behavior then negates element-wise when ``direction_val == -1``.

    For ``direction_val == +1`` the patched clone is byte-equivalent to
    upstream (no negation). The patch installs cheaply per-run and is
    reverted in try/finally by the caller.
    """
    from openmm.unit import angstrom

    def set_displacement_patched(self):
        if self.keywords.get('DISPLACEMENT') is None:
            msg = "Error: DISPLACEMENT is required"
            self._exit(msg)
        raw = self.keywords.get('DISPLACEMENT')
        # Upstream stores as (raw * angstrom). For dminus we negate the
        # raw list element-wise BEFORE multiplying — equivalent to
        # negating the displacement Vec3.
        if direction_val == -1:
            raw = [-float(x) for x in raw]
        self.displ = raw * angstrom

    return set_displacement_patched


def _restore_sys_xml_after_direction(info: Dict[str, Any]) -> None:
    """Revert the staging performed by ``_select_sys_xml_for_direction``.

    v0.9.19.1: restores BOTH .xml AND .pdb via per-file ``_restore_one_file``
    helpers when the new info-dict shape is used. Falls back to the legacy
    single-file restore for any in-flight callers that constructed an
    older-shape info dict (no ``xml_entry`` / ``pdb_entry``).

    Best-effort: any restore failure is logged to stderr but does not raise.
    """
    if not info.get("swapped"):
        return
    # v0.9.19.1 path: per-file entries.
    if info.get("xml_entry") is not None:
        _restore_one_file(info["xml_entry"])
    if info.get("pdb_entry") is not None:
        _restore_one_file(info["pdb_entry"])
    if info.get("xml_entry") is not None or info.get("pdb_entry") is not None:
        return

    # Legacy fallback (v0.9.19 info dict without per-file entries).
    active_path = info["active_sys_path"]
    backup_path = info.get("backup_path")
    preexisting_link = info.get("preexisting_symlink_target")
    try:
        if os.path.islink(active_path) or os.path.isfile(active_path):
            os.remove(active_path)
    except OSError as exc:
        print(
            f"WARN: could not remove staged sys.xml {active_path!r}: {exc}",
            file=sys.stderr,
        )
    if preexisting_link is not None:
        try:
            os.symlink(preexisting_link, active_path)
        except OSError as exc:
            print(
                f"WARN: could not restore preexisting symlink "
                f"{active_path!r} -> {preexisting_link!r}: {exc}",
                file=sys.stderr,
            )
    elif backup_path is not None and os.path.exists(backup_path):
        try:
            os.rename(backup_path, active_path)
        except OSError as exc:
            print(
                f"WARN: could not restore sys.xml backup "
                f"{backup_path!r} -> {active_path!r}: {exc}",
                file=sys.stderr,
            )


# ---------------------------------------------------------------------------
# Per-direction monkeypatch of abfe_structprep.do_equil and do_lambda_annealing
# ---------------------------------------------------------------------------
def _make_patched_do_equil(direction_val: int):
    """Return a do_equil clone with patched walker Direction.

    Upstream ``atom_openmm/abfe_structprep.py:do_equil`` hardcodes
    ``direction = 1`` at line 340. The rest of the function loads
    ``trackb_mdlambda.xml``, overrides ATMForce parameters
    (Lambda1=Lambda2=0.5, alpha=0, uh=0, w0coeff=0, Direction=direction),
    runs ``EQUILIBRATION_STEPS`` MD at lambda=0.5, and saves
    ``trackb_0.xml`` + ``trackb_0.pdb``.

    v0.9.19 (b+) corrected spec (2026-06-01):
    Walker Direction is ALWAYS +1 for both legs (two-leg ATM ABFE
    standard, Azimi 2022 §2.3). The ``direction_val`` argument now
    selects the SYSTEM (via _select_sys_xml_for_direction handled by
    the caller) and the ATMForce displacement-vector sign (via the
    OMMSystemABFE.set_displacement monkeypatch handled by the caller).
    The walker Direction parameter set INTO the ATMForce at equilibration
    time is +1 regardless.

    The patched clone:
      * forces walker Direction = +1 (corrected from upstream's hardcoded
        +1 — same value, but now justified rather than
        accidental)
      * is otherwise byte-equivalent to upstream
      * preserves the upstream parameter ladder (alpha=0, uh=0, w0coeff=0,
        umsc=1000, ubcore=500, acore=0.0625) verbatim
    """
    from atom_openmm.abfe_structprep import (
        OMMSystemABFEnoATM,
        set_platform,
    )
    from atom_openmm.ommsystem import OMMSystemABFE
    from sys import stdout
    from openmm.app import Simulation, PDBFile, XTCReporter, StateDataReporter
    from openmm.unit import (
        kelvin, kilojoules_per_mole, kilocalorie_per_mole,
        kilocalories_per_mole, nanometer,
    )
    from openmm import unit

    def do_equil_patched(keywords, logger):
        basename = keywords.get('BASENAME')
        jobname = basename

        pdbtopfile = basename + ".pdb"
        systemfile = basename + "_sys.xml"

        syst = OMMSystemABFE(basename, keywords, pdbtopfile, systemfile, logger)
        syst.create_system()

        (platform, platform_properties) = set_platform(keywords)
        simulation = Simulation(syst.topology, syst.system, syst.integrator,
                                platform, platform_properties)
        simulation.context.setPositions(syst.positions)
        if syst.boxvectors is not None:
            simulation.context.setPeriodicBoxVectors(
                syst.boxvectors[0], syst.boxvectors[1], syst.boxvectors[2]
            )

        print("Using platform %s" % simulation.context.getPlatform().getName())
        syst.barostat.setFrequency(0)

        temp = keywords.get("TEMPERATURES")
        if temp is None:
            temperature = 300.0 * kelvin
        elif isinstance(temp, list):
            temperature = float(temp[0]) * kelvin
        else:
            temperature = float(temp) * kelvin

        lmbd = 0.5
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalorie_per_mole
        uh = 0.0 * kilocalorie_per_mole
        w0coeff = 0.0 * kilocalorie_per_mole
        umsc = 1000.0 * kilocalorie_per_mole
        ubcore = 500.0 * kilocalorie_per_mole
        acore = 0.062500
        # ============================================================
        # SINGLE PATCH POINT (vs upstream do_equil L340 ``direction = 1``)
        # v0.9.19 (b+) corrected: walker Direction = +1 for
        # BOTH legs (sys_dplus + sys_dminus). direction_val argument
        # selects sys.xml + ATMForce displacement sign at the caller
        # level, NOT the walker Direction parameter.
        # ============================================================
        direction = 1  # walker Direction = +1 always (Q6 (b+) corrected)
        uoffset = 0.0 * kilocalorie_per_mole

        print("LoadState ...")
        simulation.loadState(jobname + '_mdlambda.xml')

        simulation.context.setParameter(syst.atmforce.Lambda1(), lambda1)
        simulation.context.setParameter(syst.atmforce.Lambda2(), lambda2)
        simulation.context.setParameter(syst.atmforce.Alpha(),
                                        alpha * kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.Uh(),
                                        uh / kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.W0(),
                                        w0coeff / kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.Umax(),
                                        umsc / kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.Ubcore(),
                                        ubcore / kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.Acore(), acore)
        simulation.context.setParameter(syst.atmforce.Direction(), direction)
        simulation.context.setParameter('UOffset',
                                        uoffset / kilojoules_per_mole)

        epot = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        print("Potential Energy =",
              simulation.context.getState(getEnergy=True).getPotentialEnergy())
        print("Equilibration at lambda = 1/2, walker direction = %d (Q6 b+)" % direction)

        totalSteps = int(keywords.get("EQUILIBRATION_STEPS", 150000))
        steps_per_cycle = int(keywords.get("STEPS_PER_CYCLE", 5000))
        simulation.reporters.append(
            StateDataReporter(stdout, steps_per_cycle, step=True,
                              potentialEnergy=True, temperature=True, speed=True)
        )
        if os.path.exists(jobname + "_0.xtc"):
            os.remove(jobname + "_0.xtc")
        simulation.reporters.append(
            XTCReporter(jobname + "_0.xtc", steps_per_cycle,
                        enforcePeriodicBox=False)
        )
        simulation.step(totalSteps)

        print("SaveState ...")
        simulation.saveState(jobname + "_0.xml")

        positions = simulation.context.getState(getPositions=True).getPositions()
        boxsize = simulation.context.getState().getPeriodicBoxVectors()
        simulation.topology.setPeriodicBoxVectors(boxsize)
        with open(jobname + '_0.pdb', 'w') as output:
            PDBFile.writeFile(simulation.topology, positions, output, keepIds=True)

    return do_equil_patched


def _make_patched_do_lambda_annealing(direction_val: int):
    """Return a do_lambda_annealing clone with hardcoded direction override.

    Upstream ``atom_openmm/abfe_structprep.py:do_lambda_annealing`` (L191)
    runs lambda annealing from 0.0 -> 0.5 with hardcoded direction=+1
    (L229). For direction=-1 path we want the annealing to also reflect
    the d=-1 reference potential so the binder population at lambda=0.5
    is representative of the d=-1 equilibrium (avoids a discontinuous
    Direction switch between annealing and equilibration).
    """
    from atom_openmm.ommsystem import OMMSystemABFE
    from atom_openmm.abfe_structprep import set_platform
    from sys import stdout
    from openmm.app import Simulation, PDBFile, XTCReporter, StateDataReporter
    from openmm.unit import (
        kelvin, kilojoules_per_mole, kilocalorie_per_mole, nanometer,
    )

    def do_lambda_annealing_patched(keywords, logger):
        basename = keywords.get('BASENAME')
        jobname = basename
        pdbtopfile = basename + ".pdb"
        systemfile = basename + "_sys.xml"

        syst = OMMSystemABFE(basename, keywords, pdbtopfile, systemfile, logger)
        syst.create_system()

        (platform, platform_properties) = set_platform(keywords)
        simulation = Simulation(syst.topology, syst.system, syst.integrator,
                                platform, platform_properties)
        simulation.context.setPositions(syst.positions)
        if syst.boxvectors is not None:
            simulation.context.setPeriodicBoxVectors(
                syst.boxvectors[0], syst.boxvectors[1], syst.boxvectors[2]
            )

        print("Using platform %s" % simulation.context.getPlatform().getName())
        syst.barostat.setFrequency(0)

        temp = keywords.get("TEMPERATURES")
        if temp is None:
            temperature = 300.0 * kelvin
        elif isinstance(temp, list):
            temperature = float(temp[0]) * kelvin
        else:
            temperature = float(temp) * kelvin

        lmbd = 0.0
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalorie_per_mole
        uh = 0.0 * kilocalorie_per_mole
        w0coeff = 0.0 * kilocalorie_per_mole
        umsc = 1000.0 * kilocalorie_per_mole
        ubcore = 500.0 * kilocalorie_per_mole
        acore = 0.062500
        # ============================================================
        # SINGLE PATCH POINT (vs upstream do_lambda_annealing L229)
        # v0.9.19 (b+) corrected: walker Direction = +1 for
        # BOTH legs. direction_val selects sys.xml + ATMForce
        # displacement sign at the caller level.
        # ============================================================
        direction = 1  # walker Direction = +1 always (Q6 (b+) corrected)
        uoffset = 0.0 * kilocalorie_per_mole

        print("LoadState ...")
        simulation.loadState(jobname + '_equil.xml')

        simulation.context.setParameter(syst.atmforce.Lambda1(), lambda1)
        simulation.context.setParameter(syst.atmforce.Lambda2(), lambda2)
        simulation.context.setParameter(syst.atmforce.Alpha(),
                                        alpha * kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.Uh(),
                                        uh / kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.W0(),
                                        w0coeff / kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.Umax(),
                                        umsc / kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.Ubcore(),
                                        ubcore / kilojoules_per_mole)
        simulation.context.setParameter(syst.atmforce.Acore(), acore)
        simulation.context.setParameter(syst.atmforce.Direction(), direction)
        simulation.context.setParameter('UOffset',
                                        uoffset / kilojoules_per_mole)

        epot = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        print("Potential Energy =",
              simulation.context.getState(getEnergy=True).getPotentialEnergy())
        print("Annealing to lambda = 1/2 at direction = %d ..." % direction)

        totalSteps = int(keywords.get("ANNEALING_STEPS", 250000))
        steps_per_cycle = int(keywords.get("STEPS_PER_CYCLE", 5000))
        number_of_cycles = int(totalSteps / steps_per_cycle)
        deltalambda = (0.5 - 0.0) / float(number_of_cycles)
        simulation.reporters.append(
            StateDataReporter(stdout, steps_per_cycle, step=True,
                              potentialEnergy=True, temperature=True, speed=True)
        )
        if os.path.exists(jobname + "_mdlambda.xtc"):
            os.remove(jobname + "_mdlambda.xtc")
        simulation.reporters.append(
            XTCReporter(jobname + "_mdlambda.xtc", steps_per_cycle,
                        enforcePeriodicBox=False)
        )

        binding_file = jobname + '_mdlambda.out'
        f = open(binding_file, 'w')

        for i in range(number_of_cycles):
            simulation.step(steps_per_cycle)
            state = simulation.context.getState(getEnergy=True)
            pot_energy = state.getPotentialEnergy()
            (u1, u0, ebias) = syst.atmforce.getPerturbationEnergy(simulation.context)
            umcore = simulation.context.getParameter(syst.atmforce.Umax()) * kilojoules_per_mole
            ubcore_v = simulation.context.getParameter(syst.atmforce.Ubcore()) * kilojoules_per_mole
            acore_v = simulation.context.getParameter(syst.atmforce.Acore())
            cur_direction = simulation.context.getParameter(syst.atmforce.Direction())
            if cur_direction > 0:
                pert_energy = syst.atm_utils.softCorePertE(
                    u1 - (u0 + uoffset), umcore, ubcore_v, acore_v
                )
            else:
                pert_energy = syst.atm_utils.softCorePertE(
                    (u0 + uoffset) - u1, umcore, ubcore_v, acore_v
                )
            l1 = simulation.context.getParameter(syst.atmforce.Lambda1())
            l2 = simulation.context.getParameter(syst.atmforce.Lambda2())
            a = simulation.context.getParameter(syst.atmforce.Alpha()) / kilojoules_per_mole
            umid = simulation.context.getParameter(syst.atmforce.Uh()) * kilojoules_per_mole
            w0 = simulation.context.getParameter(syst.atmforce.W0()) * kilojoules_per_mole
            print("%f %f %f %f %f %f %f %f %f" % (
                temperature / kelvin, lmbd, l1, l2,
                a * kilocalorie_per_mole, umid / kilocalorie_per_mole,
                w0 / kilocalorie_per_mole, pot_energy / kilocalorie_per_mole,
                pert_energy / kilocalorie_per_mole
            ), file=f)
            f.flush()
            lmbd += deltalambda
            lambda1 += deltalambda
            lambda2 += deltalambda
            simulation.context.setParameter(syst.atmforce.Lambda1(), lambda1)
            simulation.context.setParameter(syst.atmforce.Lambda2(), lambda2)
        f.close()

        print("SaveState ...")
        simulation.saveState(jobname + "_mdlambda.xml")
        positions = simulation.context.getState(getPositions=True).getPositions()
        boxsize = simulation.context.getState().getPeriodicBoxVectors()
        simulation.topology.setPeriodicBoxVectors(boxsize)
        with open(jobname + '_mdlambda.pdb', 'w') as output:
            PDBFile.writeFile(simulation.topology, positions, output, keepIds=True)

    return do_lambda_annealing_patched


def _make_patched_do_mintherm(velocity_seed: int):
    """Return a do_mintherm clone that seeds INDEPENDENT initial velocities.

    Velocity-seed mechanism B: each replicate gets fresh,
    INDEPENDENT initial velocities so an ensemble of per-replicate structpreps
    samples distinct initial conditions → σ_btwn (inter-replicate
    reproducibility error) becomes measurable (Wan 2021). Upstream
    ``atom_openmm/abfe_structprep.py:do_mintherm`` (L88-189) never calls
    ``setVelocitiesToTemperature`` — OpenMM auto-assigns velocities from a
    nondeterministic global RNG, so independent runs are NOT reproducibly
    seeded.

    SINGLE PATCH POINT vs upstream: after minimization + the initial
    ``syst.integrator.setTemperature(initial_temperature)`` (upstream L142-143)
    and IMMEDIATELY BEFORE the thermalization temperature-ramp loop (upstream
    L144), inject::

        simulation.context.setVelocitiesToTemperature(initial_temperature,
                                                       velocity_seed)

    The Maxwell-Boltzmann draw at ``initial_temperature`` (50 K default) with
    an explicit integer ``velocity_seed`` makes the per-replicate initial
    velocities deterministic AND independent (distinct seeds → distinct
    draws). The rest of the function is byte-equivalent to upstream
    (minimization, thermalization ramp, NPT, NVT, all saveState calls).

    This patch is direction-AGNOSTIC (do_mintherm produces the
    non-direction-dependent ``_equil.xml``); the direction override lives in
    the do_lambda_annealing / do_equil patches. Composable with them.
    """
    from atom_openmm.abfe_structprep import (
        OMMSystemABFEnoATM,
        set_platform,
    )
    from sys import stdout
    from openmm.app import Simulation, PDBFile, StateDataReporter
    from openmm.unit import kelvin

    def do_mintherm_patched(keywords, logger):
        basename = keywords.get('BASENAME')
        jobname = basename

        pdbtopfile = basename + ".pdb"
        systemfile = basename + "_sys.xml"

        # OpenMM system for minimization, thermalization, NPT, NVT
        # (does not include ATM Force).
        syst = OMMSystemABFEnoATM(basename, keywords, pdbtopfile, systemfile,
                                  logger)
        syst.create_system()

        (platform, platform_properties) = set_platform(keywords)
        simulation = Simulation(syst.topology, syst.system, syst.integrator,
                                platform, platform_properties)
        simulation.context.setPositions(syst.positions)
        if syst.boxvectors is not None:
            simulation.context.setPeriodicBoxVectors(
                syst.boxvectors[0], syst.boxvectors[1], syst.boxvectors[2])
        simulation.context.applyConstraints(0.00001)
        print("Using platform %s" % simulation.context.getPlatform().getName())

        print("Potential energy before minimization =",
              simulation.context.getState(getEnergy=True).getPotentialEnergy())
        print("Energy minimizing the system ...")
        simulation.minimizeEnergy()
        print("Potential energy after minimization =",
              simulation.context.getState(getEnergy=True).getPotentialEnergy())

        # saves minimization checkpoint
        simulation.saveState(jobname + '_min.xml')
        positions = simulation.context.getState(
            getPositions=True).getPositions()
        boxsize = simulation.context.getState().getPeriodicBoxVectors()
        simulation.topology.setPeriodicBoxVectors(boxsize)
        with open(jobname + '_min.pdb', 'w') as output:
            PDBFile.writeFile(simulation.topology, positions, output,
                              keepIds=True)

        print("Thermalization ...")

        totalSteps = int(keywords.get("THERMALIZATION_STEPS", 150000))
        steps_per_cycle = int(keywords.get("STEPS_PER_CYCLE", 5000))
        number_of_cycles = int(totalSteps / steps_per_cycle)
        simulation.reporters.append(
            StateDataReporter(stdout, steps_per_cycle, step=True,
                              potentialEnergy=True, temperature=True,
                              volume=True, speed=True))

        # initial temperature
        initial_temp = keywords.get("INITIAL_TEMPERATURE")
        if initial_temp is None:
            initial_temperature = 50.0 * kelvin
        else:
            initial_temperature = float(initial_temp) * kelvin

        final_temperature = syst.temperature
        delta_temperature = (
            (final_temperature - initial_temperature) / number_of_cycles)

        syst.barostat.setFrequency(0)  # disabled

        # MD with temperature ramp
        temperature = initial_temperature
        syst.integrator.setTemperature(temperature)

        # ============================================================
        # SINGLE PATCH POINT (velocity-seed mechanism B): seed
        # INDEPENDENT initial velocities so an ensemble of per-replicate
        # structpreps samples distinct initial conditions (→ measurable
        # σ_btwn). Distinct integer velocity_seed per replicate → distinct
        # Maxwell-Boltzmann draw at the initial temperature. Upstream
        # do_mintherm does NOT seed velocities (nondeterministic global RNG).
        # ============================================================
        simulation.context.setVelocitiesToTemperature(
            initial_temperature, int(velocity_seed))
        print("Velocity-seed (mechanism B) = %d at initial_temperature = %s"
              % (int(velocity_seed), initial_temperature))

        for i in range(number_of_cycles):
            simulation.step(steps_per_cycle)
            # prepare system for new temperature
            temperature = temperature + delta_temperature
            syst.integrator.setTemperature(temperature)

        # saves thermalized checkpoint
        simulation.saveState(jobname + '_therm.xml')
        positions = simulation.context.getState(
            getPositions=True).getPositions()
        boxsize = simulation.context.getState().getPeriodicBoxVectors()
        simulation.topology.setPeriodicBoxVectors(boxsize)
        with open(jobname + '_therm.pdb', 'w') as output:
            PDBFile.writeFile(simulation.topology, positions, output,
                              keepIds=True)

        print("NPT equilibration ...")
        syst.barostat.setFrequency(25)
        for i in range(number_of_cycles):
            simulation.step(steps_per_cycle)
        simulation.saveState(jobname + '_npt.xml')
        positions = simulation.context.getState(
            getPositions=True).getPositions()
        boxsize = simulation.context.getState().getPeriodicBoxVectors()
        simulation.topology.setPeriodicBoxVectors(boxsize)
        with open(jobname + '_npt.pdb', 'w') as output:
            PDBFile.writeFile(simulation.topology, positions, output,
                              keepIds=True)

        print("NVT equilibration ...")
        syst.barostat.setFrequency(0)  # disabled
        for i in range(number_of_cycles):
            simulation.step(steps_per_cycle)
        simulation.saveState(jobname + '_equil.xml')
        positions = simulation.context.getState(
            getPositions=True).getPositions()
        boxsize = simulation.context.getState().getPeriodicBoxVectors()
        simulation.topology.setPeriodicBoxVectors(boxsize)
        with open(jobname + '_equil.pdb', 'w') as output:
            PDBFile.writeFile(simulation.topology, positions, output,
                              keepIds=True)

    return do_mintherm_patched


# ---------------------------------------------------------------------------
# v0.9.16 VM dispatch helpers
# ---------------------------------------------------------------------------
# The structprep monkeypatch (direction=+/-1 override for do_equil +
# do_lambda_annealing) cannot be applied to a remote process via simple
# ssh-cd-execute (the patch must be loaded in the SAME Python process that
# calls ``abfe_structprep``). The VM-dispatch path therefore generates a
# self-contained Python wrapper script that re-applies the monkeypatch on
# the VM side, then ssh-executes it via the VM's ``atm`` env Python.
#
# Pattern mirrors ``trackb_per_direction_production.py:_live_launch_all_legs``
# L528-545 (``cmd = ["ssh", vm_ssh_host, remote_cmd]``) and the F3
# pre-flight gate pattern (``_gate_vm_abfe_bin_exists``).
# ---------------------------------------------------------------------------

DEFAULT_VM_PYTHON_BIN = "/home/san/miniconda3/envs/atm/bin/python"
# Match the canonical VM_SSH_HOST default declared later in the
# "Hardware host gating" block (line ~702). Duplicated here as a literal
# (not a forward reference) so this helper can be tested in isolation
# before module-level execution reaches the VM_SSH_HOST constant.
_DEFAULT_VM_SSH_HOST = "san@192.168.122.155"


def _gate_vm_leg_dir_exists(
    leg_dir: str,
    vm_ssh_host: str = _DEFAULT_VM_SSH_HOST,
    ssh_connect_timeout_s: int = 5,
    subprocess_timeout_s: int = 10,
) -> Tuple[bool, str]:
    """Pre-flight: confirm ``leg_dir`` exists on the VM AND contains the
    upstream-required inputs (``trackb.pdb`` + ``trackb_sys.xml`` +
    ``trackb_asyncre.cntl``).

    Returns ``(allow, reason)``. Refuses launch when any required input
    is missing on the VM side — cohort-safe gate semantics (better to
    halt than start half-staged).

    Same failure-class family as
    [[integrity-vm-lane-self-provisioning-20260529]] and
    [[integrity-vm-conda-path-noninteractive]]: pre-flight assertion
    that VM-side preconditions are met before any in-process commit.
    """
    required = ["trackb.pdb", "trackb_sys.xml", "trackb_asyncre.cntl"]
    probe_cmd = " && ".join(
        f"test -f {shlex.quote(os.path.join(leg_dir, f))}" for f in required
    )
    try:
        result = subprocess.run(
            [
                "ssh",
                "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
                "-o", "StrictHostKeyChecking=no",
                "-o", "LogLevel=ERROR",
                vm_ssh_host,
                probe_cmd,
            ],
            capture_output=True,
            timeout=subprocess_timeout_s,
        )
    except subprocess.TimeoutExpired:
        return (
            False,
            f"VM leg-dir probe timed out (>{subprocess_timeout_s}s) — "
            f"ssh to {vm_ssh_host} unreachable. Cannot launch.",
        )
    except FileNotFoundError as exc:
        return (
            False,
            f"local ssh binary missing for VM leg-dir probe: {exc}",
        )
    if result.returncode == 0:
        return (
            True,
            f"VM leg_dir {leg_dir!r} contains all required inputs "
            f"({', '.join(required)}) on {vm_ssh_host}",
        )
    return (
        False,
        f"VM leg_dir {leg_dir!r} missing one or more required inputs "
        f"({', '.join(required)}) on {vm_ssh_host} "
        f"(ssh test -f rc={result.returncode}). "
        f"Remediation: rsync the leg dir to the VM first, e.g. "
        f"`rsync -avz {leg_dir}/ {vm_ssh_host}:{leg_dir}/` "
        f"(VM disk must have free space; check `ssh {vm_ssh_host} df -h /home`)."
    )


# Structprep build inputs that must be present on the VM leg_dir before a
# VM-dispatched structprep. The first three are MANDATORY (mirror
# ``_gate_vm_leg_dir_exists``); the per-direction system/topology variants
# are OPTIONAL (present only for legs that use per-direction systems — the
# BOUND leg; the FREE leg is direction-agnostic and ships only the combined
# ``trackb_sys.xml`` + ``trackb.pdb``, so the optional variants are simply
# absent and skipped).
_VM_STRUCTPREP_REQUIRED_INPUTS = (
    "{jobname}.pdb",
    "{jobname}_sys.xml",
    "{jobname}_asyncre.cntl",
)
_VM_STRUCTPREP_OPTIONAL_INPUTS = (
    "{jobname}_sys_dplus.xml",
    "{jobname}_sys_dminus.xml",
    "{jobname}_dplus.pdb",
    "{jobname}_dminus.pdb",
)


def _rsync_leg_inputs_to_vm(
    leg_dir: str,
    vm_ssh_host: str = _DEFAULT_VM_SSH_HOST,
    jobname: str = "trackb",
    ssh_connect_timeout_s: int = 5,
    subprocess_timeout_s: int = 600,
) -> Dict[str, Any]:
    """Provision a leg's structprep BUILD inputs onto the VM ``leg_dir``.

    Task 2 (free-system V100 provisioning): the FREE leg combined system
    (``trackb_sys.xml`` ~4292 particles, direction-agnostic) + densified38v4
    cntl (``trackb_asyncre.cntl`` — the schedule SSOT, materialized host-side
    by the v2_asyncre launcher) + topology (``trackb.pdb``) must be on the VM
    before ``--gpu-host=vm`` structprep can build the per-direction
    ``_0_{dplus,dminus}.xml`` base states there.

    Mirrors the production launcher's
    ``trackb_per_direction_production.py:_rsync_canonical_base_state_to_vm``
    pattern (``rsync -a --files-from`` + ssh-side ``test -f`` verify) so the
    structprep and production lanes provision the VM identically.

    The MANDATORY inputs (``trackb.pdb`` / ``trackb_sys.xml`` /
    ``trackb_asyncre.cntl``) must exist on the HOST or this raises
    ``RuntimeError`` (cohort-safe fail-fast — no partial provision). OPTIONAL
    per-direction variants are pushed only when present on the host (BOUND
    leg ships them; FREE leg does not).

    Why a dedicated helper rather than ``rsync -avz leg_dir/`` wholesale: the
    leg dir can already contain large VM-incompatible artifacts (per-replica
    ckpts, prior ``_0_*.xml``); we push only the deterministic build inputs
    so re-provisioning is cheap and never clobbers VM-side per-direction
    outputs that a prior structprep produced.

    NOTE: schedule code (``get_schedule('densified38v4')``) is NOT required
    VM-side — the structprep VM wrapper imports only ``atom_openmm`` and reads
    the cntl. The cntl carries the 40-state densified38v4 λ ladder verbatim,
    so rsyncing the cntl is the complete schedule-provisioning step.

    Returns ``{"status", "pushed", "skipped_absent", "leg_dir"}``.
    Raises RuntimeError on missing mandatory input / rsync failure / VM
    verify failure (cohort halt).
    """
    leg_dir = leg_dir.rstrip("/")
    required = [t.format(jobname=jobname)
                for t in _VM_STRUCTPREP_REQUIRED_INPUTS]
    optional = [t.format(jobname=jobname)
                for t in _VM_STRUCTPREP_OPTIONAL_INPUTS]

    missing_required = [
        rel for rel in required
        if not os.path.isfile(os.path.join(leg_dir, rel))
    ]
    if missing_required:
        raise RuntimeError(
            "_rsync_leg_inputs_to_vm: host leg_dir missing mandatory "
            f"structprep input(s) {missing_required} under {leg_dir}. "
            "Build the leg (v2_3 free system + densified38v4 cntl) on the "
            "host first; cohort halted (no partial VM provision)."
        )

    rel_paths = list(required)
    skipped_absent = []
    for rel in optional:
        if os.path.isfile(os.path.join(leg_dir, rel)):
            rel_paths.append(rel)
        else:
            skipped_absent.append(rel)

    ssh_e = (
        f"ssh -o ConnectTimeout={ssh_connect_timeout_s} "
        f"-o StrictHostKeyChecking=no -o LogLevel=ERROR"
    )
    # Ensure the VM leg_dir parent path exists before rsync (rsync -a does
    # not mkdir intermediate dirs of the destination root).
    mkdir_cmd = [
        "ssh",
        "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
        "-o", "StrictHostKeyChecking=no",
        "-o", "LogLevel=ERROR",
        vm_ssh_host,
        f"mkdir -p {shlex.quote(leg_dir)}",
    ]
    try:
        mk = subprocess.run(mkdir_cmd, capture_output=True,
                            timeout=60)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VM leg_dir mkdir timed out — ssh to {vm_ssh_host} "
            "unresponsive. Cohort halted."
        ) from exc
    if mk.returncode != 0:
        raise RuntimeError(
            f"VM leg_dir mkdir FAILED (rc={mk.returncode}): "
            f"{mk.stderr.decode(errors='replace')[:300]}"
        )

    import tempfile
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".rsync_leg_inputs", delete=False
    ) as fh:
        fh.write("\n".join(rel_paths) + "\n")
        files_from_path = fh.name
    try:
        rsync_cmd = [
            "rsync", "-a", "-e", ssh_e,
            f"--files-from={files_from_path}",
            f"{leg_dir}/",
            f"{vm_ssh_host}:{leg_dir}/",
        ]
        try:
            rsync_result = subprocess.run(
                rsync_cmd, capture_output=True,
                timeout=subprocess_timeout_s,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"rsync of leg structprep inputs timed out "
                f"(>{subprocess_timeout_s}s). Cohort halted."
            ) from exc
        if rsync_result.returncode != 0:
            raise RuntimeError(
                f"rsync leg structprep inputs FAILED (rc="
                f"{rsync_result.returncode}): "
                f"{rsync_result.stderr.decode(errors='replace')[:500]}"
            )
    finally:
        try:
            os.unlink(files_from_path)
        except OSError:
            pass

    # VM-side verify: all MANDATORY inputs present after the push.
    verify_cmd = " && ".join(
        f"test -f {shlex.quote(os.path.join(leg_dir, rel))}"
        for rel in required
    )
    try:
        verify_result = subprocess.run(
            [
                "ssh",
                "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
                "-o", "StrictHostKeyChecking=no",
                "-o", "LogLevel=ERROR",
                vm_ssh_host, verify_cmd,
            ],
            capture_output=True, timeout=60,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VM leg-input verify timed out — ssh to {vm_ssh_host} "
            "unresponsive. Cohort halted."
        ) from exc
    if verify_result.returncode != 0:
        raise RuntimeError(
            "VM leg-input push failed: one or more mandatory inputs "
            f"({', '.join(required)}) missing on VM at {leg_dir}. "
            "Cohort halt."
        )
    return {
        "status": "rsynced",
        "pushed": rel_paths,
        "skipped_absent": skipped_absent,
        "leg_dir": leg_dir,
    }


# Per-direction structprep OUTPUTS that the VM produces and the host must
# receive back (VM-local-only filesystem; host is local ext4, NOT shared/NFS
# with the VM). ``.xml`` is MANDATORY (the base state the production readiness
# gate + host-side subdir staging require); ``.pdb`` is OPTIONAL (the
# direction-final-frame snapshot — present when abfe_structprep emitted it).
_VM_STRUCTPREP_OUTPUT_TEMPLATES = (
    ("{jobname}_0_{tag}.xml", True),    # mandatory
    ("{jobname}_0_{tag}.pdb", False),   # optional
)


def _rsync_per_direction_outputs_from_vm(
    leg_dir: str,
    direction_tag: str,
    vm_ssh_host: str = _DEFAULT_VM_SSH_HOST,
    jobname: str = "trackb",
    ssh_connect_timeout_s: int = 5,
    subprocess_timeout_s: int = 300,
) -> Dict[str, Any]:
    """Pull ONE direction's structprep base-state outputs from the VM back to
    the host ``leg_dir`` (VM->host; reverse of
    ``_rsync_canonical_base_state_to_vm`` /
    ``_rsync_leg_inputs_to_vm`` which are host->VM).

    The ``--gpu-host=vm`` structprep runs ``abfe_structprep`` on the VM and the
    VM wrapper writes the renamed ``{jobname}_0_{tag}.{xml,pdb}`` IN-PLACE on
    the VM leg_dir. The host filesystem is local ext4 (NOT shared/NFS with the
    VM), so those genuine VM-produced, velocity-seed-distinct base states are
    invisible to the host until pulled back. Both downstream consumers run
    HOST-side and read the HOST leg_dir:

      * the production readiness gate
        ``trackb_per_direction_production.py:check_free_pilot_readiness``
        (requires ``{jobname}_0_dplus.xml`` + ``{jobname}_0_dminus.xml`` in
        the host leg_dir);
      * host-side per-direction subdir staging
        ``stage_per_direction_subdir`` (copies the host
        ``{jobname}_0_{tag}.xml`` into each work subdir; that per-subdir base
        state is what then travels to the VM at production time via
        ``_rsync_subdir_to_vm``, where ``ommworker.py:265`` loads it as the
        per-direction-BASENAME ``{jobname}_{tag}_0.xml``).

    INTEGRITY-CRITICAL (do NOT fake a host file): this helper FAILS LOUD
    (``RuntimeError``) when the mandatory VM-side output is absent, rather than
    creating an empty/placeholder host file. A faked host base state would seed
    the wrong (non-seeded or stale) replicate velocities -> a fake-independent
    replicate -> invalid σ_btwn. The rsync moves the GENUINE
    VM-produced seeded base state; the VM-side ``test -f`` pre-check below is
    the fail-fast gate that guarantees we never silently fabricate.

    Mirrors ``_rsync_leg_inputs_to_vm`` machinery (ssh + ``rsync -a
    --files-from`` + verify) but with src/dst reversed and a VM-side
    pre-check ``ssh ... test -f`` before the pull (so an absent VM output is a
    clear, early RuntimeError rather than a confusing post-rsync host-miss).
    Real files (no symlink) -> plain ``rsync -a`` (no ``-L``). Cohort-safe:
    raises on any failure (mkdir / ssh / rsync / verify / absent VM output).

    Returns ``{"status", "pulled", "skipped_absent", "leg_dir",
    "direction_tag"}``.
    """
    if direction_tag not in ("dplus", "dminus"):
        raise ValueError(
            f"direction_tag must be 'dplus' or 'dminus', got "
            f"{direction_tag!r}"
        )
    leg_dir = leg_dir.rstrip("/")
    templated = [
        (t.format(jobname=jobname, tag=direction_tag), mandatory)
        for (t, mandatory) in _VM_STRUCTPREP_OUTPUT_TEMPLATES
    ]
    mandatory_rels = [rel for rel, m in templated if m]
    optional_rels = [rel for rel, m in templated if not m]

    ssh_opts = [
        "-o", f"ConnectTimeout={ssh_connect_timeout_s}",
        "-o", "StrictHostKeyChecking=no",
        "-o", "LogLevel=ERROR",
    ]

    # VM-side pre-check: every MANDATORY output must exist on the VM BEFORE we
    # pull. This is the fail-loud guard against faking a host base state.
    precheck_cmd = " && ".join(
        f"test -f {shlex.quote(os.path.join(leg_dir, rel))}"
        for rel in mandatory_rels
    )
    try:
        pre = subprocess.run(
            ["ssh", *ssh_opts, vm_ssh_host, precheck_cmd],
            capture_output=True, timeout=60,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VM per-direction output pre-check timed out — ssh to "
            f"{vm_ssh_host} unresponsive. Cohort halted."
        ) from exc
    if pre.returncode != 0:
        raise RuntimeError(
            "_rsync_per_direction_outputs_from_vm: mandatory VM-produced "
            f"output(s) {mandatory_rels} ABSENT on VM at {leg_dir} "
            f"(host={vm_ssh_host}). Refusing to fabricate a host base state "
            "(a faked base state would seed the wrong replicate "
            "velocities -> invalid sigma_btwn). Verify the VM structprep "
            "actually produced the seeded base state, then re-run."
        )

    # Determine which OPTIONAL outputs are present on the VM (pull only those).
    present_optional: List[str] = []
    skipped_absent: List[str] = []
    for rel in optional_rels:
        try:
            chk = subprocess.run(
                ["ssh", *ssh_opts, vm_ssh_host,
                 f"test -f {shlex.quote(os.path.join(leg_dir, rel))}"],
                capture_output=True, timeout=60,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"VM per-direction optional-output probe timed out — ssh to "
                f"{vm_ssh_host} unresponsive. Cohort halted."
            ) from exc
        if chk.returncode == 0:
            present_optional.append(rel)
        else:
            skipped_absent.append(rel)

    rel_paths = list(mandatory_rels) + present_optional

    # Ensure the host leg_dir exists before the pull (the leg_dir was built on
    # the host, so it normally exists; mkdir -p is defensive + idempotent).
    os.makedirs(leg_dir, exist_ok=True)

    ssh_e = (
        f"ssh -o ConnectTimeout={ssh_connect_timeout_s} "
        f"-o StrictHostKeyChecking=no -o LogLevel=ERROR"
    )
    import tempfile
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".rsync_outputs_back", delete=False
    ) as fh:
        fh.write("\n".join(rel_paths) + "\n")
        files_from_path = fh.name
    try:
        rsync_cmd = [
            "rsync", "-a", "-e", ssh_e,
            f"--files-from={files_from_path}",
            f"{vm_ssh_host}:{leg_dir}/",
            f"{leg_dir}/",
        ]
        try:
            rsync_result = subprocess.run(
                rsync_cmd, capture_output=True,
                timeout=subprocess_timeout_s,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"rsync-back of per-direction outputs timed out "
                f"(>{subprocess_timeout_s}s). Cohort halted."
            ) from exc
        if rsync_result.returncode != 0:
            raise RuntimeError(
                f"rsync-back per-direction outputs FAILED (rc="
                f"{rsync_result.returncode}): "
                f"{rsync_result.stderr.decode(errors='replace')[:500]}"
            )
    finally:
        try:
            os.unlink(files_from_path)
        except OSError:
            pass

    # Host-side verify: the mandatory XML must now be a real file on the host.
    for rel in mandatory_rels:
        host_path = os.path.join(leg_dir, rel)
        if not os.path.isfile(host_path):
            raise RuntimeError(
                "VM->host rsync-back reported success but mandatory output "
                f"{host_path} not present on host. The VM-produced base "
                "state was NOT moved; refusing to proceed (no fabrication)."
            )
    return {
        "status": "pulled",
        "pulled": rel_paths,
        "skipped_absent": skipped_absent,
        "leg_dir": leg_dir,
        "direction_tag": direction_tag,
    }


def _generate_vm_structprep_wrapper(
    leg_dir: str,
    direction_val: int,
    jobname: str,
    cntl_basename: str,
    velocity_seed: Optional[int] = None,
) -> str:
    """Generate a self-contained Python wrapper that re-applies the
    direction-patched monkeypatch on the VM side and invokes
    ``abfe_structprep``.

    The wrapper inlines minimal copies of ``_make_patched_do_equil`` and
    ``_make_patched_do_lambda_annealing`` (textually identical to the host
    helpers above, byte-for-byte where it matters — the upstream parameter
    ladder is preserved verbatim). No host project imports are needed on
    VM (the wrapper depends only on ``atom_openmm``).

    ``velocity_seed`` (velocity-seed mechanism B; v0.9.30 VM parity):
    when set, the wrapper additionally inlines a copy of
    ``_make_patched_do_mintherm`` (byte-equivalent to the host helper, single
    patch point = ``setVelocitiesToTemperature(initial_temperature, seed)``
    immediately before the thermalization ramp) and assigns
    ``upstream.do_mintherm = _make_patched_do_mintherm(VELOCITY_SEED)`` in the
    VM-side ``main()`` BEFORE ``abfe_structprep`` runs. The seed is baked into
    the generated source as the compile-time constant ``VELOCITY_SEED = <int>``
    so it travels with the wrapper over ``ssh ... <python> -`` (stdin pipe) and
    is genuinely applied inside the VM-side OpenMM ``Simulation.context`` — the
    VM-produced ``_equil.xml`` (and therefore the per-direction ``_0_*.xml``)
    is seed-distinct, NOT silently identical. When ``velocity_seed is None``
    (default, single-run path) the wrapper does NOT touch ``do_mintherm`` and
    is byte-for-byte identical to the prior behaviour.

    INTEGRITY-CRITICAL: if ``velocity_seed`` were dropped here, the VM
    replicates would share OpenMM's nondeterministic global-RNG velocities →
    fake-independent replicates → invalid σ_btwn. The constant is therefore
    embedded literally and the wrapper's ``main()`` asserts it was applied.

    Returns the wrapper script as a UTF-8 string for the caller to write
    to the VM-side temp path. Wrapper exits non-zero on any failure; the
    host caller maps that to ``status="error"``.

    NOTE: The wrapper writes the renamed ``_0_{direction_tag}.xml`` and
    ``_0_{direction_tag}.pdb`` outputs IN-PLACE inside ``leg_dir`` on the
    VM, matching the host-side dispatch contract. The host caller does
    not need to rsync-back results when the leg_dir is shared via NFS;
    when leg_dir is VM-local-only (common case), the caller must rsync
    the per-direction outputs back to the host before the audit step.
    """
    if direction_val not in (1, -1):
        raise ValueError(
            f"direction_val must be 1 or -1, got {direction_val}"
        )
    if velocity_seed is not None and not isinstance(velocity_seed, int):
        raise ValueError(
            f"velocity_seed must be an int or None, got "
            f"{type(velocity_seed).__name__}"
        )
    direction_tag = "dplus" if direction_val == 1 else "dminus"
    # textwrap.dedent so the wrapper renders left-justified regardless of
    # the caller indentation level.
    wrapper = textwrap.dedent(f'''\
        #!/usr/bin/env python
        """VM-side per-direction structprep wrapper.

        Auto-generated by trackb_per_direction_structprep.py
        (_generate_vm_structprep_wrapper). DO NOT EDIT — regenerate
        from the host script if the patch surface changes.

        v0.9.19 (2026-06-01): wrapper additionally stages
        trackb_sys_<tag>.xml -> trackb_sys.xml (symlink) before invoking
        abfe_structprep so OMMSystemABFE consumes the direction-specific
        system XML (binder pre-displaced into the correct bulk region).
        Restored in finally block. Legacy single-system builds (only
        trackb_sys.xml present) no-op the swap.

        v0.9.30 (2026-06-06): when VELOCITY_SEED is not None the wrapper
        re-applies _make_patched_do_mintherm(VELOCITY_SEED) on the VM side
        so each replicate's initial velocities are seed-INDEPENDENT
        (setVelocitiesToTemperature(T, seed) before the thermalization
        ramp). The seed is a compile-time constant baked into THIS source,
        so it genuinely reaches the VM-side OpenMM Simulation (no silent
        drop). VELOCITY_SEED is None on the single-run path (do_mintherm
        untouched, byte-equivalent to pre-v0.9.30 behaviour).

        Direction: {direction_val:+d} ({direction_tag})
        Leg dir:   {leg_dir}
        Jobname:   {jobname}
        Cntl:      {cntl_basename}
        VeloSeed:  {velocity_seed}
        """
        import os
        import shutil
        import sys
        import time

        LEG_DIR = {leg_dir!r}
        JOBNAME = {jobname!r}
        CNTL_BASENAME = {cntl_basename!r}
        DIRECTION_VAL = {direction_val}
        DIRECTION_TAG = {direction_tag!r}
        VELOCITY_SEED = {velocity_seed!r}


        # ----------------------------------------------------------------
        # Per-direction system XML + topology PDB staging
        # (mirror of host helpers _stage_one_file + _restore_one_file +
        # _select_sys_xml_for_direction + _restore_sys_xml_after_direction
        # — v0.9.19.1 swaps BOTH .xml AND .pdb).
        # ----------------------------------------------------------------
        def _stage_one_file(per_dir_path, per_dir_basename, active_path,
                            backup_basename):
            entry = {{
                "active_path": active_path,
                "per_direction_path": per_dir_path,
                "swapped": False,
                "backup_path": None,
                "preexisting_symlink_target": None,
            }}
            if os.path.islink(active_path):
                link_target = os.readlink(active_path)
                if (link_target == per_dir_basename
                        or os.path.abspath(link_target)
                        == os.path.abspath(per_dir_path)):
                    # No-op idempotent: record preexisting target so
                    # restore recreates the symlink (v0.9.19.2 fix).
                    entry["preexisting_symlink_target"] = link_target
                    entry["swapped"] = True
                    return entry
                entry["preexisting_symlink_target"] = link_target
                os.remove(active_path)
            elif os.path.isfile(active_path):
                backup_path = os.path.join(LEG_DIR, backup_basename)
                if os.path.exists(backup_path):
                    ts = time.strftime("%Y%m%dT%H%M%S")
                    backup_path = backup_path + "." + ts
                os.rename(active_path, backup_path)
                entry["backup_path"] = backup_path
            try:
                os.symlink(per_dir_basename, active_path)
            except OSError:
                shutil.copy2(per_dir_path, active_path)
            entry["swapped"] = True
            return entry


        def _restore_one_file(entry):
            if not entry or not entry.get("swapped"):
                return
            active_path = entry["active_path"]
            backup_path = entry.get("backup_path")
            preexisting_link = entry.get("preexisting_symlink_target")
            try:
                if (os.path.islink(active_path)
                        or os.path.isfile(active_path)):
                    os.remove(active_path)
            except OSError as exc:
                print("WARN: could not remove staged file: " + repr(exc),
                      file=sys.stderr)
            if preexisting_link is not None:
                try:
                    os.symlink(preexisting_link, active_path)
                except OSError as exc:
                    print("WARN: could not restore preexisting link: "
                          + repr(exc), file=sys.stderr)
            elif backup_path is not None and os.path.exists(backup_path):
                try:
                    os.rename(backup_path, active_path)
                except OSError as exc:
                    print("WARN: could not restore backup: "
                          + repr(exc), file=sys.stderr)


        def _select_sys_xml():
            # XML
            xml_per_dir_basename = (JOBNAME + "_sys_" + DIRECTION_TAG
                                    + ".xml")
            xml_per_dir_path = os.path.join(LEG_DIR, xml_per_dir_basename)
            xml_active_path = os.path.join(LEG_DIR, JOBNAME + "_sys.xml")
            xml_backup_basename = (JOBNAME + "_sys.xml.bak_"
                                   + DIRECTION_TAG)
            # PDB
            pdb_per_dir_basename = (JOBNAME + "_" + DIRECTION_TAG + ".pdb")
            pdb_per_dir_path = os.path.join(LEG_DIR, pdb_per_dir_basename)
            pdb_active_path = os.path.join(LEG_DIR, JOBNAME + ".pdb")
            pdb_backup_basename = (JOBNAME + ".pdb.bak_" + DIRECTION_TAG)

            info = {{
                "swapped": False,
                "per_direction_path": None,
                "active_sys_path": xml_active_path,
                "backup_path": None,
                "preexisting_symlink_target": None,
                "xml_entry": None,
                "pdb_entry": None,
            }}
            if not os.path.isfile(xml_per_dir_path):
                return info
            if not os.path.isfile(pdb_per_dir_path):
                raise RuntimeError(
                    "Per-direction XML " + repr(xml_per_dir_path)
                    + " exists but matching PDB "
                    + repr(pdb_per_dir_path)
                    + " is missing. Both must be present for the "
                    "per-direction structprep to work — rebuild via "
                    "scripts/phase4_trackB_v2_make_system.py "
                    "--bound-directions " + DIRECTION_TAG + "."
                )
            info["per_direction_path"] = xml_per_dir_path
            xml_entry = _stage_one_file(xml_per_dir_path,
                                        xml_per_dir_basename,
                                        xml_active_path,
                                        xml_backup_basename)
            pdb_entry = _stage_one_file(pdb_per_dir_path,
                                        pdb_per_dir_basename,
                                        pdb_active_path,
                                        pdb_backup_basename)
            info["xml_entry"] = xml_entry
            info["pdb_entry"] = pdb_entry
            info["backup_path"] = xml_entry.get("backup_path")
            info["preexisting_symlink_target"] = xml_entry.get(
                "preexisting_symlink_target")
            info["swapped"] = bool(
                xml_entry["swapped"] and pdb_entry["swapped"])
            return info


        def _restore_sys_xml(info):
            if not info.get("swapped"):
                return
            if info.get("xml_entry") is not None:
                _restore_one_file(info["xml_entry"])
            if info.get("pdb_entry") is not None:
                _restore_one_file(info["pdb_entry"])

        # ----------------------------------------------------------------
        # Re-implements _make_patched_do_mintherm from host script
        # (v0.9.30 VM velocity-seed parity, mechanism B).
        # SINGLE PATCH POINT vs upstream do_mintherm: inject
        # setVelocitiesToTemperature(initial_temperature, velocity_seed)
        # immediately before the thermalization ramp so this replicate's
        # initial velocities are seed-INDEPENDENT (distinct integer seeds
        # -> distinct Maxwell-Boltzmann draws -> measurable sigma_btwn).
        # Byte-equivalent to the host helper except this VM copy is invoked
        # only when VELOCITY_SEED is not None. Direction-AGNOSTIC
        # (do_mintherm produces the non-direction-dependent _equil.xml).
        # ----------------------------------------------------------------
        def _make_patched_do_mintherm(velocity_seed):
            from atom_openmm.abfe_structprep import (
                OMMSystemABFEnoATM,
                set_platform,
            )
            from sys import stdout
            from openmm.app import Simulation, PDBFile, StateDataReporter
            from openmm.unit import kelvin

            def do_mintherm_patched(keywords, logger):
                basename = keywords.get('BASENAME')
                jobname = basename
                pdbtopfile = basename + ".pdb"
                systemfile = basename + "_sys.xml"
                syst = OMMSystemABFEnoATM(basename, keywords, pdbtopfile,
                                          systemfile, logger)
                syst.create_system()
                (platform, platform_properties) = set_platform(keywords)
                simulation = Simulation(syst.topology, syst.system,
                                        syst.integrator, platform,
                                        platform_properties)
                simulation.context.setPositions(syst.positions)
                if syst.boxvectors is not None:
                    simulation.context.setPeriodicBoxVectors(
                        syst.boxvectors[0], syst.boxvectors[1],
                        syst.boxvectors[2])
                simulation.context.applyConstraints(0.00001)
                print("Using platform %s"
                      % simulation.context.getPlatform().getName())
                print("Potential energy before minimization =",
                      simulation.context.getState(
                          getEnergy=True).getPotentialEnergy())
                print("Energy minimizing the system ...")
                simulation.minimizeEnergy()
                print("Potential energy after minimization =",
                      simulation.context.getState(
                          getEnergy=True).getPotentialEnergy())
                simulation.saveState(jobname + '_min.xml')
                positions = simulation.context.getState(
                    getPositions=True).getPositions()
                boxsize = simulation.context.getState(
                    ).getPeriodicBoxVectors()
                simulation.topology.setPeriodicBoxVectors(boxsize)
                with open(jobname + '_min.pdb', 'w') as output:
                    PDBFile.writeFile(simulation.topology, positions,
                                      output, keepIds=True)
                print("Thermalization ...")
                totalSteps = int(keywords.get("THERMALIZATION_STEPS",
                                              150000))
                steps_per_cycle = int(keywords.get("STEPS_PER_CYCLE", 5000))
                number_of_cycles = int(totalSteps / steps_per_cycle)
                simulation.reporters.append(
                    StateDataReporter(stdout, steps_per_cycle, step=True,
                                      potentialEnergy=True, temperature=True,
                                      volume=True, speed=True))
                initial_temp = keywords.get("INITIAL_TEMPERATURE")
                if initial_temp is None:
                    initial_temperature = 50.0 * kelvin
                else:
                    initial_temperature = float(initial_temp) * kelvin
                final_temperature = syst.temperature
                delta_temperature = (
                    (final_temperature - initial_temperature)
                    / number_of_cycles)
                syst.barostat.setFrequency(0)  # disabled
                temperature = initial_temperature
                syst.integrator.setTemperature(temperature)
                # === SINGLE PATCH POINT (velocity-seed mechanism B): seed
                # INDEPENDENT initial velocities on the VM side so the
                # VM-produced _equil.xml is genuinely seed-distinct. ===
                simulation.context.setVelocitiesToTemperature(
                    initial_temperature, int(velocity_seed))
                print("Velocity-seed (mechanism B, VM) = %d at "
                      "initial_temperature = %s"
                      % (int(velocity_seed), initial_temperature))
                for i in range(number_of_cycles):
                    simulation.step(steps_per_cycle)
                    temperature = temperature + delta_temperature
                    syst.integrator.setTemperature(temperature)
                simulation.saveState(jobname + '_therm.xml')
                positions = simulation.context.getState(
                    getPositions=True).getPositions()
                boxsize = simulation.context.getState(
                    ).getPeriodicBoxVectors()
                simulation.topology.setPeriodicBoxVectors(boxsize)
                with open(jobname + '_therm.pdb', 'w') as output:
                    PDBFile.writeFile(simulation.topology, positions,
                                      output, keepIds=True)
                print("NPT equilibration ...")
                syst.barostat.setFrequency(25)
                for i in range(number_of_cycles):
                    simulation.step(steps_per_cycle)
                simulation.saveState(jobname + '_npt.xml')
                positions = simulation.context.getState(
                    getPositions=True).getPositions()
                boxsize = simulation.context.getState(
                    ).getPeriodicBoxVectors()
                simulation.topology.setPeriodicBoxVectors(boxsize)
                with open(jobname + '_npt.pdb', 'w') as output:
                    PDBFile.writeFile(simulation.topology, positions,
                                      output, keepIds=True)
                print("NVT equilibration ...")
                syst.barostat.setFrequency(0)  # disabled
                for i in range(number_of_cycles):
                    simulation.step(steps_per_cycle)
                simulation.saveState(jobname + '_equil.xml')
                positions = simulation.context.getState(
                    getPositions=True).getPositions()
                boxsize = simulation.context.getState(
                    ).getPeriodicBoxVectors()
                simulation.topology.setPeriodicBoxVectors(boxsize)
                with open(jobname + '_equil.pdb', 'w') as output:
                    PDBFile.writeFile(simulation.topology, positions,
                                      output, keepIds=True)

            return do_mintherm_patched

        # Re-implements _make_patched_do_equil from host script.
        # v0.9.19 (b+) corrected: walker Direction is +1 for
        # BOTH legs. DIRECTION_VAL selects the sys.xml + ATMForce
        # displacement sign (via _select_sys_xml + the
        # _make_patched_set_displacement wrapper below), NOT the walker
        # Direction parameter. (Upstream's abfe_structprep.py:L340
        # hardcoded direction=1; we keep that value but reinterpret it
        # per the corrected ATM two-leg ABFE spec.)
        def _make_patched_do_equil(direction_val):
            from atom_openmm.ommsystem import OMMSystemABFE
            from atom_openmm.abfe_structprep import set_platform
            from sys import stdout
            from openmm.app import (
                Simulation, PDBFile, XTCReporter, StateDataReporter,
            )
            from openmm.unit import (
                kelvin, kilojoules_per_mole, kilocalorie_per_mole,
            )

            def do_equil_patched(keywords, logger):
                basename = keywords.get("BASENAME")
                jobname = basename
                pdbtopfile = basename + ".pdb"
                systemfile = basename + "_sys.xml"

                syst = OMMSystemABFE(basename, keywords, pdbtopfile,
                                     systemfile, logger)
                syst.create_system()

                (platform, platform_properties) = set_platform(keywords)
                simulation = Simulation(syst.topology, syst.system,
                                        syst.integrator, platform,
                                        platform_properties)
                simulation.context.setPositions(syst.positions)
                if syst.boxvectors is not None:
                    simulation.context.setPeriodicBoxVectors(
                        syst.boxvectors[0], syst.boxvectors[1],
                        syst.boxvectors[2],
                    )

                print("Using platform %s"
                      % simulation.context.getPlatform().getName())
                syst.barostat.setFrequency(0)

                temp = keywords.get("TEMPERATURES")
                if temp is None:
                    temperature = 300.0 * kelvin
                elif isinstance(temp, list):
                    temperature = float(temp[0]) * kelvin
                else:
                    temperature = float(temp) * kelvin

                lmbd = 0.5
                lambda1 = lmbd
                lambda2 = lmbd
                alpha = 0.0 / kilocalorie_per_mole
                uh = 0.0 * kilocalorie_per_mole
                w0coeff = 0.0 * kilocalorie_per_mole
                umsc = 1000.0 * kilocalorie_per_mole
                ubcore = 500.0 * kilocalorie_per_mole
                acore = 0.0625
                # v0.9.19 Q6 (b+): walker Direction = +1 always.
                direction = 1
                uoffset = 0.0 * kilocalorie_per_mole

                print("LoadState ...")
                simulation.loadState(jobname + "_mdlambda.xml")
                simulation.context.setParameter(
                    syst.atmforce.Lambda1(), lambda1)
                simulation.context.setParameter(
                    syst.atmforce.Lambda2(), lambda2)
                simulation.context.setParameter(
                    syst.atmforce.Alpha(),
                    alpha * kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.Uh(), uh / kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.W0(), w0coeff / kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.Umax(), umsc / kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.Ubcore(),
                    ubcore / kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.Acore(), acore)
                simulation.context.setParameter(
                    syst.atmforce.Direction(), direction)
                simulation.context.setParameter(
                    "UOffset", uoffset / kilojoules_per_mole)

                print("Equilibration at lambda=1/2, direction=%d ..."
                      % direction)
                totalSteps = int(keywords.get(
                    "EQUILIBRATION_STEPS", 150000))
                steps_per_cycle = int(keywords.get(
                    "STEPS_PER_CYCLE", 5000))
                simulation.reporters.append(
                    StateDataReporter(stdout, steps_per_cycle,
                                      step=True, potentialEnergy=True,
                                      temperature=True, speed=True))
                if os.path.exists(jobname + "_0.xtc"):
                    os.remove(jobname + "_0.xtc")
                simulation.reporters.append(
                    XTCReporter(jobname + "_0.xtc",
                                steps_per_cycle,
                                enforcePeriodicBox=False))
                simulation.step(totalSteps)

                print("SaveState ...")
                simulation.saveState(jobname + "_0.xml")
                positions = simulation.context.getState(
                    getPositions=True).getPositions()
                boxsize = simulation.context.getState(
                    ).getPeriodicBoxVectors()
                simulation.topology.setPeriodicBoxVectors(boxsize)
                with open(jobname + "_0.pdb", "w") as out:
                    PDBFile.writeFile(simulation.topology, positions,
                                      out, keepIds=True)
            return do_equil_patched


        def _make_patched_do_lambda_annealing(direction_val):
            from atom_openmm.ommsystem import OMMSystemABFE
            from atom_openmm.abfe_structprep import set_platform
            from sys import stdout
            from openmm.app import (
                Simulation, PDBFile, XTCReporter, StateDataReporter,
            )
            from openmm.unit import (
                kelvin, kilojoules_per_mole, kilocalorie_per_mole,
            )

            def do_lambda_annealing_patched(keywords, logger):
                basename = keywords.get("BASENAME")
                jobname = basename
                pdbtopfile = basename + ".pdb"
                systemfile = basename + "_sys.xml"

                syst = OMMSystemABFE(basename, keywords, pdbtopfile,
                                     systemfile, logger)
                syst.create_system()

                (platform, platform_properties) = set_platform(keywords)
                simulation = Simulation(syst.topology, syst.system,
                                        syst.integrator, platform,
                                        platform_properties)
                simulation.context.setPositions(syst.positions)
                if syst.boxvectors is not None:
                    simulation.context.setPeriodicBoxVectors(
                        syst.boxvectors[0], syst.boxvectors[1],
                        syst.boxvectors[2],
                    )

                print("Using platform %s"
                      % simulation.context.getPlatform().getName())
                syst.barostat.setFrequency(0)

                temp = keywords.get("TEMPERATURES")
                if temp is None:
                    temperature = 300.0 * kelvin
                elif isinstance(temp, list):
                    temperature = float(temp[0]) * kelvin
                else:
                    temperature = float(temp) * kelvin

                lmbd = 0.0
                lambda1 = lmbd
                lambda2 = lmbd
                alpha = 0.0 / kilocalorie_per_mole
                uh = 0.0 * kilocalorie_per_mole
                w0coeff = 0.0 * kilocalorie_per_mole
                umsc = 1000.0 * kilocalorie_per_mole
                ubcore = 500.0 * kilocalorie_per_mole
                acore = 0.0625
                # v0.9.19 Q6 (b+): walker Direction = +1 always.
                direction = 1
                uoffset = 0.0 * kilocalorie_per_mole

                print("LoadState ...")
                simulation.loadState(jobname + "_equil.xml")
                simulation.context.setParameter(
                    syst.atmforce.Lambda1(), lambda1)
                simulation.context.setParameter(
                    syst.atmforce.Lambda2(), lambda2)
                simulation.context.setParameter(
                    syst.atmforce.Alpha(),
                    alpha * kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.Uh(), uh / kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.W0(), w0coeff / kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.Umax(), umsc / kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.Ubcore(),
                    ubcore / kilojoules_per_mole)
                simulation.context.setParameter(
                    syst.atmforce.Acore(), acore)
                simulation.context.setParameter(
                    syst.atmforce.Direction(), direction)
                simulation.context.setParameter(
                    "UOffset", uoffset / kilojoules_per_mole)

                nsteps = int(keywords.get("ANNEALING_STEPS", 250000))
                steps_per_cycle = int(keywords.get(
                    "STEPS_PER_CYCLE", 5000))
                # Linear lambda ramp 0.0 -> 0.5 over nsteps
                ncycles = nsteps // steps_per_cycle
                lmb_step = 0.5 / max(ncycles, 1)
                simulation.reporters.append(
                    StateDataReporter(stdout, steps_per_cycle,
                                      step=True, potentialEnergy=True,
                                      temperature=True, speed=True))
                if os.path.exists(jobname + "_mdlambda.xtc"):
                    os.remove(jobname + "_mdlambda.xtc")
                simulation.reporters.append(
                    XTCReporter(jobname + "_mdlambda.xtc",
                                steps_per_cycle,
                                enforcePeriodicBox=False))
                lmb = 0.0
                for cyc in range(ncycles):
                    lmb = min(0.5, lmb + lmb_step)
                    simulation.context.setParameter(
                        syst.atmforce.Lambda1(), lmb)
                    simulation.context.setParameter(
                        syst.atmforce.Lambda2(), lmb)
                    simulation.step(steps_per_cycle)

                print("SaveState ...")
                simulation.saveState(jobname + "_mdlambda.xml")
                positions = simulation.context.getState(
                    getPositions=True).getPositions()
                boxsize = simulation.context.getState(
                    ).getPeriodicBoxVectors()
                simulation.topology.setPeriodicBoxVectors(boxsize)
                with open(jobname + "_mdlambda.pdb", "w") as out:
                    PDBFile.writeFile(simulation.topology, positions,
                                      out, keepIds=True)
            return do_lambda_annealing_patched


        def _make_patched_set_displacement(direction_val):
            """v0.9.19 Q6 (b+) negate ATMForce displacement for dminus."""
            from openmm.unit import angstrom

            def set_displacement_patched(self):
                if self.keywords.get("DISPLACEMENT") is None:
                    msg = "Error: DISPLACEMENT is required"
                    self._exit(msg)
                raw = self.keywords.get("DISPLACEMENT")
                if direction_val == -1:
                    raw = [-float(x) for x in raw]
                self.displ = raw * angstrom
            return set_displacement_patched


        def main():
            import atom_openmm.abfe_structprep as upstream
            from atom_openmm.ommsystem import OMMSystemABFE
            os.chdir(LEG_DIR)
            # v0.9.19: stage per-direction sys.xml BEFORE OMMSystemABFE
            # construction inside abfe_structprep, restore in finally.
            sys_xml_info = _select_sys_xml()
            print("# sys.xml stage: swapped=%r per_direction=%r"
                  % (sys_xml_info.get("swapped"),
                     sys_xml_info.get("per_direction_path")))
            upstream.do_equil = _make_patched_do_equil(DIRECTION_VAL)
            upstream.do_lambda_annealing = (
                _make_patched_do_lambda_annealing(DIRECTION_VAL))
            # v0.9.30 VM velocity-seed parity (mechanism B): when a
            # per-replicate VELOCITY_SEED is baked into this wrapper, re-apply
            # _make_patched_do_mintherm on the VM side so the VM-produced
            # _equil.xml (and the per-direction _0_*.xml that anneals from it)
            # is genuinely seed-distinct. INTEGRITY-CRITICAL: this is the line
            # that prevents fake-independent VM replicates (invalid sigma_btwn).
            orig_mintherm = upstream.do_mintherm
            if VELOCITY_SEED is not None:
                upstream.do_mintherm = _make_patched_do_mintherm(
                    int(VELOCITY_SEED))
                if upstream.do_mintherm is orig_mintherm:
                    print("ERROR: velocity-seed patch did not take effect "
                          "(do_mintherm unchanged)", file=sys.stderr)
                    sys.exit(3)
                print("# velocity-seed (mechanism B, VM) applied: seed=%d"
                      % int(VELOCITY_SEED))
            else:
                print("# velocity-seed: None (upstream do_mintherm, "
                      "single-run path)")
            # v0.9.19 Q6 (b+) ATMForce displacement sign reversal for dminus
            orig_set_displ = OMMSystemABFE.set_displacement
            OMMSystemABFE.set_displacement = (
                _make_patched_set_displacement(DIRECTION_VAL))
            print("# ATMForce displacement sign: %s"
                  % ("+1 (dplus)" if DIRECTION_VAL == 1
                     else "-1 (dminus, manually reversed)"))
            t0 = time.time()
            try:
                upstream.abfe_structprep(CNTL_BASENAME)
            finally:
                _restore_sys_xml(sys_xml_info)
                OMMSystemABFE.set_displacement = orig_set_displ
                upstream.do_mintherm = orig_mintherm
            wall = time.time() - t0

            produced_xml = JOBNAME + "_0.xml"
            produced_pdb = JOBNAME + "_0.pdb"
            target_xml = JOBNAME + "_0_" + DIRECTION_TAG + ".xml"
            target_pdb = JOBNAME + "_0_" + DIRECTION_TAG + ".pdb"
            if not os.path.isfile(produced_xml):
                print("ERROR: abfe_structprep did not produce "
                      + produced_xml, file=sys.stderr)
                sys.exit(1)
            shutil.move(produced_xml, target_xml)
            if os.path.isfile(produced_pdb):
                shutil.move(produced_pdb, target_pdb)
            # C1 root-cause (v0.9.25): restore the canonical _0.{{xml,pdb}}
            # on the VM from the dplus variant so the abfe_production worker
            # (ommworker.py:265 loadState(_0.xml)) finds it. Only the dplus
            # dispatch restores it; dminus must NOT overwrite (dplus
            # already did). Copy (not move) keeps the variant intact.
            if DIRECTION_TAG == "dplus":
                _canon_xml = JOBNAME + "_0.xml"
                _dplus_xml = JOBNAME + "_0_dplus.xml"
                if (os.path.isfile(_dplus_xml)
                        and not os.path.isfile(_canon_xml)):
                    shutil.copy2(_dplus_xml, _canon_xml)
                    _canon_pdb = JOBNAME + "_0.pdb"
                    _dplus_pdb = JOBNAME + "_0_dplus.pdb"
                    if (os.path.isfile(_dplus_pdb)
                            and not os.path.isfile(_canon_pdb)):
                        shutil.copy2(_dplus_pdb, _canon_pdb)
            print("OK wrapper produced %s in %.1fs"
                  % (target_xml, wall))
            sys.exit(0)


        if __name__ == "__main__":
            main()
    ''')
    return wrapper


def _run_per_direction_structprep_via_ssh(
    leg_dir: str,
    direction_val: int,
    jobname: str,
    cntl_path: str,
    vm_ssh_host: str,
    vm_python_bin: str,
    log_path: str,
    velocity_seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Inner helper: dispatch ONE direction's structprep to the VM via
    a stdin-piped Python wrapper. Returns dict with status/wall/log.

    The wrapper is generated in-memory by ``_generate_vm_structprep_wrapper``
    and piped to ``ssh VM <python> -`` via stdin so no leftover file is
    left on the VM after the run (avoids tmpfile cleanup race).

    ``velocity_seed`` (v0.9.30 VM parity): forwarded into the generated
    wrapper as the compile-time ``VELOCITY_SEED`` constant so the seed
    travels over the stdin pipe and is genuinely applied inside the VM-side
    OpenMM ``Simulation.context`` (mechanism B). ``None`` = single-run path.
    """
    direction_tag = "dplus" if direction_val == 1 else "dminus"
    cntl_basename = os.path.basename(cntl_path)
    wrapper_src = _generate_vm_structprep_wrapper(
        leg_dir=leg_dir,
        direction_val=direction_val,
        jobname=jobname,
        cntl_basename=cntl_basename,
        velocity_seed=velocity_seed,
    )
    # Command: ssh VM '<python> -' (reads wrapper from stdin)
    ssh_cmd = [
        "ssh",
        "-o", "StrictHostKeyChecking=no",
        "-o", "LogLevel=ERROR",
        vm_ssh_host,
        f"{shlex.quote(vm_python_bin)} -",
    ]
    print(f"  [VM] dispatch {direction_tag} to {vm_ssh_host}: "
          f"{' '.join(shlex.quote(c) for c in ssh_cmd)}")
    t0 = time.time()
    with open(log_path, "w") as logfh:
        logfh.write(f"# VM dispatch: {' '.join(ssh_cmd)}\n")
        logfh.write(f"# leg_dir:     {leg_dir}\n")
        logfh.write(f"# direction:   {direction_val:+d}\n")
        logfh.write(f"# velocity_seed: {velocity_seed}\n")
        logfh.write(f"# started:     {time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
        logfh.write("# --- wrapper src follows (head) ---\n")
        for line in wrapper_src.splitlines()[:8]:
            logfh.write(f"# {line}\n")
        logfh.write("# --- wrapper output ---\n")
        logfh.flush()
        proc = subprocess.run(
            ssh_cmd,
            input=wrapper_src,
            stdout=logfh,
            stderr=subprocess.STDOUT,
            text=True,
        )
    wall = time.time() - t0
    if proc.returncode != 0:
        return {
            "status": "error",
            "rc": proc.returncode,
            "wall_seconds": round(wall, 1),
            "log_path": log_path,
            "error": (
                f"VM structprep failed (rc={proc.returncode}); "
                f"see {log_path}"
            ),
        }
    return {
        "status": "produced",
        "rc": 0,
        "wall_seconds": round(wall, 1),
        "log_path": log_path,
    }


# ---------------------------------------------------------------------------
# Per-direction structprep driver (per leg, per direction)
# ---------------------------------------------------------------------------
def run_per_direction_structprep(
    leg_dir: str,
    direction_val: int,
    jobname: str = "trackb",
    smoke: bool = False,
    gpu_host: str = "local",
    vm_ssh_host: str = "san@192.168.122.155",
    vm_python_bin: str = DEFAULT_VM_PYTHON_BIN,
    velocity_seed: Optional[int] = None,
    rsync_outputs_from_vm: bool = False,
) -> Dict[str, Any]:
    """Run structprep in ``leg_dir`` with hardcoded direction override.

    ``velocity_seed`` (velocity-seed mechanism B): when set, patches
    ``do_mintherm`` to call ``setVelocitiesToTemperature(T, velocity_seed)``
    before the thermalization ramp so this replicate gets INDEPENDENT initial
    velocities (per-replicate distinct integer → distinct Maxwell-Boltzmann
    draw → measurable σ_btwn). When None (default, single-run path) the
    upstream do_mintherm is used unchanged (no velocity seeding) — byte-for-byte
    identical to the prior behavior.

    Replaces ``atom_openmm.abfe_structprep.do_lambda_annealing`` and
    ``do_equil`` with direction-patched clones, then invokes
    ``abfe_structprep(cntl_path)`` exactly as upstream would. The
    resulting ``{jobname}_0.xml`` is renamed to
    ``{jobname}_0_dplus.xml`` (direction=+1) or
    ``{jobname}_0_dminus.xml`` (direction=-1).

    Pre-equilibration steps ``do_mintherm`` (minimization +
    thermalization + NPT + NVT) are NOT direction-dependent and
    produce ``{jobname}_equil.xml`` which is reused as the input to
    the direction-patched annealing.

    ``gpu_host`` (v0.9.16 SSH-dispatch fix):
      * ``"vm"``    — dispatch wrapper script via ssh to ``vm_ssh_host``;
                      monkeypatch is applied in the VM-side Python
                      process (matches production launcher pattern at
                      ``trackb_per_direction_production.py:535``).
      * ``"local"`` — in-process monkeypatch + ``abfe_structprep`` on
                      the host (legacy path, used when host 5070Ti is
                      free).
      * ``"cpu"``   — in-process monkeypatch with CUDA disabled
                      (``CUDA_VISIBLE_DEVICES=""``), OpenMM CPU
                      platform. Slow but safe for audit / smoke.

    Returns dict with output paths + sha256 + size.

    ``rsync_outputs_from_vm`` ((A) VM->host rsync-BACK, v0.9.31): only
    meaningful with ``gpu_host=="vm"``. When True, after the VM wrapper
    succeeds the GENUINE VM-produced base state(s)
    ``{jobname}_0_{tag}.{xml,pdb}`` are pulled BACK to the host leg_dir via
    ``_rsync_per_direction_outputs_from_vm`` (fail-loud if absent on the VM —
    no fabrication). This makes the host leg_dir contain the real seeded
    base state so the HOST-side production readiness gate
    (``check_free_pilot_readiness``) and HOST-side subdir staging
    (``stage_per_direction_subdir``) resolve it.

    NOTE on the VM-side production canonical (part (B) of the V100 free
    campaign gap): the live ``--gpu-host vm`` PRODUCTION path does NOT consume
    a leg-level ``{jobname}_0.xml`` on the VM. The 2-process split
    (``trackb_per_direction_production.py:stage_per_direction_subdir`` +
    ``_rsync_subdir_to_vm``) materializes a self-contained per-direction
    subdir whose BASENAME is ``{jobname}_{tag}`` and whose base state is
    ``{jobname}_{tag}_0.xml`` (a copy of the host leg-level
    ``{jobname}_0_{tag}.xml``), then rsyncs that ENTIRE subdir to the VM. The
    VM production worker (``ommworker.py:265`` loadState ``BASENAME + "_0.xml"``
    == ``{jobname}_{tag}_0.xml``) reads the per-subdir base state, NOT a
    leg-level ``{jobname}_0.xml``. So the VM canonical for production is
    provisioned by the production subdir rsync — the rsync-back of part (A) is
    the load-bearing step (it puts the per-direction base state on the host
    where staging picks it up). The VM wrapper's inline dplus->canonical cp
    (below) is retained ONLY for the legacy single-dispatch consumer and is
    harmless for the 2-process path.

    NOTE: caller must already have validated ``gate_gpu_host`` PASS +
    (for ``vm`` mode) ``_gate_vm_leg_dir_exists`` PASS. This driver
    does not re-gate (caller is ``main()``).
    """
    if direction_val not in (1, -1):
        raise ValueError(f"direction_val must be 1 or -1, got {direction_val}")
    if gpu_host not in ("vm", "local", "cpu"):
        raise ValueError(
            f"gpu_host must be one of vm|local|cpu, got {gpu_host!r}"
        )
    direction_tag = "dplus" if direction_val == 1 else "dminus"

    cntl_path = os.path.join(leg_dir, jobname + "_asyncre.cntl")
    if not os.path.isfile(cntl_path):
        raise RuntimeError(f"Missing cntl: {cntl_path}")
    pdb_path = os.path.join(leg_dir, jobname + ".pdb")
    sys_path = os.path.join(leg_dir, jobname + "_sys.xml")
    if not os.path.isfile(pdb_path) or not os.path.isfile(sys_path):
        raise RuntimeError(
            f"Missing {jobname}.pdb or {jobname}_sys.xml in {leg_dir}; "
            f"run v2.1 build_all_four_systems first"
        )

    target_xml = os.path.join(leg_dir, jobname + f"_0_{direction_tag}.xml")
    target_pdb = os.path.join(leg_dir, jobname + f"_0_{direction_tag}.pdb")
    if os.path.isfile(target_xml):
        st = os.stat(target_xml)
        return {
            "leg_dir": leg_dir,
            "direction": direction_val,
            "status": "cached",
            "target_xml": target_xml,
            "target_pdb": target_pdb if os.path.isfile(target_pdb) else None,
            "size_bytes": st.st_size,
            "sha256": _sha256(target_xml),
        }

    log_path = os.path.join(leg_dir, f"_structprep_{direction_tag}.log")
    print(f"\n=== STRUCTPREP {leg_dir} direction={direction_val:+d} "
          f"gpu_host={gpu_host} ===")

    # Velocity-seed mechanism B (v0.9.30 VM parity): the VM ssh-dispatch
    # wrapper now plumbs the seed into its self-contained do_mintherm
    # re-implementation (_generate_vm_structprep_wrapper bakes it in as the
    # VELOCITY_SEED compile-time constant; the wrapper's main() applies
    # _make_patched_do_mintherm(VELOCITY_SEED) before abfe_structprep). The
    # earlier NotImplementedError silent-drop guard is therefore lifted
    # for the VM path. Defense-in-depth type check remains so a non-int seed
    # fails fast rather than silently disabling seeding via a bad literal.
    if velocity_seed is not None and not isinstance(velocity_seed, int):
        raise ValueError(
            f"velocity_seed must be an int or None, got "
            f"{type(velocity_seed).__name__}"
        )

    if gpu_host == "vm":
        # v0.9.16 fix: dispatch to VM via ssh + stdin-piped wrapper.
        # The VM-side wrapper writes the renamed _0_{tag}.xml + _0_{tag}.pdb
        # in-place; no host-side rename needed when VM and host share the
        # same filesystem (NFS/shared) OR when the operator rsyncs the
        # result back before the audit step.
        # v0.9.30: velocity_seed is forwarded into the wrapper so the
        # VM-produced base state is genuinely seed-distinct (mechanism B).
        ssh_result = _run_per_direction_structprep_via_ssh(
            leg_dir=leg_dir,
            direction_val=direction_val,
            jobname=jobname,
            cntl_path=cntl_path,
            vm_ssh_host=vm_ssh_host,
            vm_python_bin=vm_python_bin,
            log_path=log_path,
            velocity_seed=velocity_seed,
        )
        if ssh_result["status"] != "produced":
            return {
                "leg_dir": leg_dir,
                "direction": direction_val,
                "status": "error",
                "error": ssh_result.get("error", "unknown VM error"),
                "rc": ssh_result.get("rc"),
                "log_path": log_path,
                "wall_seconds": ssh_result.get("wall_seconds"),
                "dispatch": "ssh_vm",
                "vm_ssh_host": vm_ssh_host,
            }
        # (A) VM->host rsync-back (v0.9.31): the host filesystem is local
        # ext4 (NOT shared/NFS with the VM), so the genuine VM-produced,
        # velocity-seed-distinct base state is invisible to the host until
        # pulled back. When ``rsync_outputs_from_vm`` is set, pull the
        # mandatory ``{jobname}_0_{tag}.xml`` (+ optional .pdb) back so the
        # HOST-side readiness gate (check_free_pilot_readiness) and HOST-side
        # subdir staging (stage_per_direction_subdir) find the genuine seeded
        # base state. The helper FAILS LOUD if the VM output is absent —
        # never fabricates a host file (a faked base state would seed
        # the wrong replicate velocities -> invalid sigma_btwn).
        rsync_back_info: Optional[Dict[str, Any]] = None
        if rsync_outputs_from_vm:
            rsync_back_info = _rsync_per_direction_outputs_from_vm(
                leg_dir=leg_dir,
                direction_tag=direction_tag,
                vm_ssh_host=vm_ssh_host,
                jobname=jobname,
            )
        # When leg_dir is shared with VM (NFS/9p), the renamed XML is
        # already present on host. When VM-local-only AND the rsync-back was
        # NOT requested, the operator must rsync it back before this exists
        # locally — surface a clear error (do NOT fabricate).
        if not os.path.isfile(target_xml):
            return {
                "leg_dir": leg_dir,
                "direction": direction_val,
                "status": "error",
                "error": (
                    f"VM produced wrapper succeeded (rc=0) but "
                    f"{target_xml} not visible on host. The VM leg_dir "
                    f"may not be shared (NFS/9p) — re-run with "
                    f"--rsync-outputs-from-vm to pull the per-direction "
                    f"base states back, or rsync manually: `rsync -avz "
                    f"{vm_ssh_host}:{leg_dir}/{jobname}_0_{direction_tag}.* "
                    f"{leg_dir}/`"
                ),
                "log_path": log_path,
                "wall_seconds": ssh_result.get("wall_seconds"),
                "dispatch": "ssh_vm",
                "vm_ssh_host": vm_ssh_host,
            }
        return {
            "leg_dir": leg_dir,
            "direction": direction_val,
            "status": "produced",
            "target_xml": target_xml,
            "target_pdb": target_pdb if os.path.isfile(target_pdb) else None,
            "size_bytes": os.stat(target_xml).st_size,
            "sha256": _sha256(target_xml),
            "wall_seconds": ssh_result.get("wall_seconds"),
            "log_path": log_path,
            "dispatch": "ssh_vm",
            "vm_ssh_host": vm_ssh_host,
            "velocity_seed": velocity_seed,
            "rsync_back": rsync_back_info,
        }

    # ----------------------------------------------------------------
    # local / cpu: in-process monkeypatch + abfe_structprep
    # ----------------------------------------------------------------
    # v0.9.19 (2026-06-01): stage trackb_sys_<tag>.xml as the active
    # trackb_sys.xml BEFORE invoking upstream so OMMSystemABFE picks up the
    # per-direction system (binder physically pre-displaced in the correct
    # bulk-solvent region). Legacy single-system builds (only trackb_sys.xml
    # present) no-op the stage and run exactly as before.
    sys_xml_stage_info = _select_sys_xml_for_direction(
        leg_dir, jobname, direction_tag,
    )

    # Monkeypatch upstream in-process. Save originals so we can restore.
    # v0.9.19 (Q6 (b+) corrected): also patch OMMSystemABFE.set_displacement
    # to negate the ATMForce displacement vector for dminus walkers
    # (per-particle FixedDisplacement does NOT auto-flip with Direction —
    # see _make_patched_set_displacement docstring).
    import atom_openmm.abfe_structprep as upstream
    from atom_openmm.ommsystem import OMMSystemABFE
    orig_annealing = upstream.do_lambda_annealing
    orig_equil = upstream.do_equil
    orig_mintherm = upstream.do_mintherm
    orig_set_displacement = OMMSystemABFE.set_displacement
    upstream.do_lambda_annealing = _make_patched_do_lambda_annealing(direction_val)
    upstream.do_equil = _make_patched_do_equil(direction_val)
    # Velocity-seed mechanism B: seed independent initial velocities when a
    # per-replicate velocity_seed is supplied (else upstream do_mintherm).
    if velocity_seed is not None:
        upstream.do_mintherm = _make_patched_do_mintherm(velocity_seed)
    OMMSystemABFE.set_displacement = _make_patched_set_displacement(direction_val)

    t0 = time.time()
    cwd_save = os.getcwd()
    try:
        # Smoke override via env (consumed by parse_config in cntl).
        # Smoke mode reduces THERMALIZATION_STEPS / ANNEALING_STEPS /
        # EQUILIBRATION_STEPS via cntl file (already set by v2.1 launcher
        # when --smoke is passed). Here we just respect the cntl as-is.
        os.chdir(leg_dir)
        upstream.abfe_structprep(os.path.basename(cntl_path))
        # Rename output to direction-tagged name
        produced_xml = os.path.join(leg_dir, jobname + "_0.xml")
        produced_pdb = os.path.join(leg_dir, jobname + "_0.pdb")
        if not os.path.isfile(produced_xml):
            raise RuntimeError(
                f"abfe_structprep did not produce {produced_xml}"
            )
        # Move (not copy) — second-direction run will overwrite _0.xml
        # but original dplus _0.xml must already be renamed first.
        shutil.move(produced_xml, target_xml)
        if os.path.isfile(produced_pdb):
            shutil.move(produced_pdb, target_pdb)
    finally:
        os.chdir(cwd_save)
        # Restore originals (avoid polluting other in-process invocations).
        upstream.do_lambda_annealing = orig_annealing
        upstream.do_equil = orig_equil
        upstream.do_mintherm = orig_mintherm
        OMMSystemABFE.set_displacement = orig_set_displacement
        # Restore the active sys.xml to its pre-stage state.
        _restore_sys_xml_after_direction(sys_xml_stage_info)

    wall = time.time() - t0
    return {
        "leg_dir": leg_dir,
        "direction": direction_val,
        "status": "produced",
        "target_xml": target_xml,
        "target_pdb": target_pdb if os.path.isfile(target_pdb) else None,
        "size_bytes": os.stat(target_xml).st_size,
        "sha256": _sha256(target_xml),
        "wall_seconds": round(wall, 1),
        "log_path": log_path,
        "dispatch": "in_process",
        "gpu_host": gpu_host,
        "sys_xml_swap_applied": bool(sys_xml_stage_info.get("swapped")),
        "per_direction_sys_xml": sys_xml_stage_info.get("per_direction_path"),
        "velocity_seed": velocity_seed,
    }


# ---------------------------------------------------------------------------
# Sanity audit (dry-run input)
# ---------------------------------------------------------------------------
def audit_per_direction_xml(
    leg_dir: str,
    jobname: str = "trackb",
) -> Dict[str, Any]:
    """Per-direction starting-XML sanity audit.

    For each leg's ``trackb_0_dplus.xml`` and ``trackb_0_dminus.xml``:
      - Load via XmlSerializer.deserialize → State
      - Extract positions, velocities, box
      - Identify binder chain L atoms (vs receptor chain A) from topology
      - Compute binder centroid for each direction → diff should approx
        +/- displacement_nm[0] (validates per-direction equilibration
        sampled different physical regions)
      - CYS SG-SG dist for cyclic_ss check (must be ~ 0.20-0.22 nm,
        std < 0.1 nm across both directions)
    """
    import openmm as mm
    from openmm import unit, XmlSerializer
    from openmm.app import PDBFile
    import numpy as np

    pdb_path = os.path.join(leg_dir, jobname + ".pdb")
    if not os.path.isfile(pdb_path):
        raise RuntimeError(f"Missing topology PDB: {pdb_path}")
    pdb = PDBFile(pdb_path)

    # Identify binder chain ("L" upstream-renamed for bound, "B" for free)
    # Identify CYS SG atom indices for cyclic_ss audit
    binder_atom_indices: List[int] = []
    sg_atom_indices: List[int] = []
    for chain in pdb.topology.chains():
        if chain.id in ("L", "B"):
            for res in chain.residues():
                for atom in res.atoms():
                    binder_atom_indices.append(atom.index)
                if res.name in ("CYS", "CYX"):
                    for atom in res.atoms():
                        if atom.name == "SG":
                            sg_atom_indices.append(atom.index)

    results: Dict[str, Any] = {
        "leg_dir": leg_dir,
        "topology_pdb": pdb_path,
        "n_binder_atoms": len(binder_atom_indices),
        "n_cys_sg": len(sg_atom_indices),
        "directions": {},
    }

    centroids: Dict[str, Optional[List[float]]] = {}
    for direction_tag, direction_val in (("dplus", 1), ("dminus", -1)):
        xml = os.path.join(leg_dir, jobname + f"_0_{direction_tag}.xml")
        if not os.path.isfile(xml):
            results["directions"][f"d={direction_val:+d}"] = {
                "status": "missing",
                "xml": xml,
            }
            centroids[direction_tag] = None
            continue
        with open(xml) as fh:
            state = XmlSerializer.deserialize(fh.read())

        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        # binder centroid
        if binder_atom_indices:
            bp = pos[binder_atom_indices]
            binder_centroid = bp.mean(axis=0).tolist()
        else:
            binder_centroid = None
        # cyclic_ss SG-SG
        sg_dist_nm = None
        if len(sg_atom_indices) >= 2:
            d = pos[sg_atom_indices[0]] - pos[sg_atom_indices[1]]
            import numpy as np
            sg_dist_nm = float(np.linalg.norm(d))

        results["directions"][f"d={direction_val:+d}"] = {
            "status": "ok",
            "xml": xml,
            "size_bytes": os.stat(xml).st_size,
            "sha256": _sha256(xml),
            "box_nm": [[float(x) for x in row] for row in box],
            "binder_centroid_nm": binder_centroid,
            "sg_sg_distance_nm": sg_dist_nm,
        }
        centroids[direction_tag] = binder_centroid

    # Per-direction binder centroid diff (should reflect different equilibria)
    if centroids.get("dplus") and centroids.get("dminus"):
        import numpy as np
        c_plus = np.array(centroids["dplus"])
        c_minus = np.array(centroids["dminus"])
        diff_nm = (c_plus - c_minus).tolist()
        diff_mag_nm = float(np.linalg.norm(c_plus - c_minus))
        results["binder_centroid_diff_nm"] = diff_nm
        results["binder_centroid_diff_magnitude_nm"] = diff_mag_nm

    # Cross-direction SG-SG std (cyclic_ss intact in both directions)
    sg_vals = [
        results["directions"][k].get("sg_sg_distance_nm")
        for k in ("d=+1", "d=-1")
        if results["directions"].get(k, {}).get("sg_sg_distance_nm") is not None
    ]
    if len(sg_vals) == 2:
        import numpy as np
        results["sg_sg_std_nm_across_directions"] = float(np.std(sg_vals))
        results["sg_sg_max_nm"] = float(np.max(sg_vals))
        results["sg_sg_min_nm"] = float(np.min(sg_vals))
        # cyclic_ss intact criterion: 0.18 < SG-SG < 0.25 nm AND std < 0.05 nm
        # Cast to native bool — `np.std(...) < 0.05` returns numpy.bool_
        # which is not JSON-serializable by Python stdlib json (numpy 2.x).
        sg_intact = bool(
            all(0.18 < v < 0.25 for v in sg_vals)
            and (float(np.std(sg_vals)) < 0.05)
        )
        results["cyclic_ss_intact_both_directions"] = sg_intact

    return results


# ---------------------------------------------------------------------------
# Charge audit (Q1(a) integer rescale persistence check)
# ---------------------------------------------------------------------------
def audit_charge_axis_persistence(
    leg_dir: str,
    jobname: str = "trackb",
    tol_e: float = 5e-4,
) -> Dict[str, Any]:
    """Confirm per-residue |Sigma q| <= tol_e after structprep.

    The Q1(a) integer rescale is applied at system XML build time
    (utils/atm_trackB_setup.py + scripts/trackb_rescale_mtr_xml.py).
    Per condition C4 the rescale must persist into Phase 4 production —
    re-verify after structprep by inspecting the bare system XML
    (structprep does NOT rebuild the system; it only equilibrates the
    state).
    """
    import openmm as mm
    from openmm import XmlSerializer, NonbondedForce
    from openmm.app import PDBFile

    sys_xml = os.path.join(leg_dir, jobname + "_sys.xml")
    pdb_path = os.path.join(leg_dir, jobname + ".pdb")
    if not os.path.isfile(sys_xml):
        raise RuntimeError(f"Missing system XML: {sys_xml}")
    with open(sys_xml) as fh:
        system = XmlSerializer.deserialize(fh.read())
    pdb = PDBFile(pdb_path)

    # Per-residue charge sum
    nbforces = [system.getForce(i) for i in range(system.getNumForces())
                if isinstance(system.getForce(i), NonbondedForce)]
    if not nbforces:
        raise RuntimeError("No NonbondedForce in system")
    nb = nbforces[0]

    residue_sigma_q: List[Dict[str, Any]] = []
    max_abs_sigma_q = 0.0
    violating: List[Dict[str, Any]] = []
    for res in pdb.topology.residues():
        if res.name in ("HOH", "WAT", "NA", "CL", "K", "MG", "ZN"):
            continue
        sigma_q = 0.0
        for atom in res.atoms():
            q, sig, eps = nb.getParticleParameters(atom.index)
            try:
                sigma_q += q.value_in_unit(q.unit)
            except AttributeError:
                sigma_q += float(q)
        residue_sigma_q.append({
            "chain": res.chain.id,
            "res_name": res.name,
            "res_id": res.id,
            "sigma_q_e": round(sigma_q, 6),
        })
        if abs(sigma_q - round(sigma_q)) > tol_e:
            violating.append({
                "chain": res.chain.id,
                "res_name": res.name,
                "res_id": res.id,
                "sigma_q_e": round(sigma_q, 6),
                "fractional_part_e": round(sigma_q - round(sigma_q), 6),
            })
        max_abs_sigma_q = max(max_abs_sigma_q, abs(sigma_q - round(sigma_q)))

    return {
        "leg_dir": leg_dir,
        "system_xml": sys_xml,
        "n_residues_checked": len(residue_sigma_q),
        "max_abs_fractional_sigma_q_e": round(max_abs_sigma_q, 6),
        "tolerance_e": tol_e,
        "all_residues_within_tol": max_abs_sigma_q <= tol_e,
        "n_violating": len(violating),
        "violating_residues": violating[:20],  # cap to first 20
    }


def _restore_canonical_from_dplus(
    leg_dir: str, jobname: str = "trackb",
) -> Dict[str, Any]:
    """Restore canonical ``{jobname}_0.{xml,pdb}`` from the dplus variant.

    Idempotent (writes only when the canonical is absent; never
    overwrites an existing canonical — e.g. an in-progress run's restart
    state). Used after per-direction structprep, whose renames remove the
    canonical that the abfe_production worker (``ommworker.py:265
    loadState({jobname}_0.xml)``) hardcodes. Mirrors the launcher backstop
    ``trackb_per_direction_production._ensure_canonical_base_state``.

    canonical=dplus is a Direction=+1 lambda=0.5 equilibrium consumed only
    by the service worker (``compute=False``); every sampling walker loads
    its own per-direction ckpt (r0..r10=dplus, r11..r21=dminus), so the
    backward (r11..r21) replicas are NOT contaminated.

    Returns ``{"created_xml", "created_pdb", "xml", "pdb"}``.
    """
    leg_dir = os.path.abspath(leg_dir)
    src_xml = os.path.join(leg_dir, jobname + "_0_dplus.xml")
    src_pdb = os.path.join(leg_dir, jobname + "_0_dplus.pdb")
    canon_xml = os.path.join(leg_dir, jobname + "_0.xml")
    canon_pdb = os.path.join(leg_dir, jobname + "_0.pdb")
    created_xml = False
    created_pdb = False
    if os.path.isfile(src_xml) and not os.path.isfile(canon_xml):
        shutil.copy2(src_xml, canon_xml)
        created_xml = True
    if os.path.isfile(src_pdb) and not os.path.isfile(canon_pdb):
        shutil.copy2(src_pdb, canon_pdb)
        created_pdb = True
    return {
        "created_xml": created_xml,
        "created_pdb": created_pdb,
        "xml": canon_xml,
        "pdb": canon_pdb,
    }


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Hardware host gating (v0.9.10 device-guard fix, 2026-05-31)
# ---------------------------------------------------------------------------
#
# Replaces the broken "device 1" assumption with a host-aware gate. Each box
# (host + VM) has exactly one GPU at index 0; "device 1" never existed. The
# choice of which GPU is encoded by ``--gpu-host`` (vm | local | cpu).
#
# Both gates are READ-ONLY. SSH issues only ``nvidia-smi --query-gpu=...``
# queries to the VM; never kill, never reset, never touch GUI session.
# ---------------------------------------------------------------------------

VM_SSH_HOST = "san@192.168.122.155"
VM_SSH_OPTS = ["-o", "BatchMode=yes", "-o", "ConnectTimeout=8"]
GPU_UTIL_REFUSE_THRESHOLD_PCT = 30


def _query_vm_gpu_utilization_pct(
    ssh_host: str = VM_SSH_HOST,
    timeout_s: int = 12,
) -> Dict[str, Any]:
    """Read-only nvidia-smi query against the VM. Returns dict with keys
    ``util_pct`` (Optional[float], None on query failure), ``raw`` (raw
    stdout for debugging), and ``error`` (str if any).

    Implementation: ``ssh -o BatchMode=yes -o ConnectTimeout=8 san@<vm>
    'nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits'``.
    Always read-only — never sends a destructive subcommand.
    """
    import subprocess
    cmd = [
        "ssh",
        *VM_SSH_OPTS,
        ssh_host,
        "nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits",
    ]
    try:
        out = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        return {"util_pct": None, "raw": "", "error": "ssh_timeout"}
    except FileNotFoundError:
        return {"util_pct": None, "raw": "", "error": "ssh_not_found"}
    if out.returncode != 0:
        return {
            "util_pct": None, "raw": out.stdout,
            "error": f"ssh_rc={out.returncode}: {out.stderr[:200].strip()}",
        }
    text = out.stdout.strip()
    # nvidia-smi returns one row per GPU; VM has exactly one so first line.
    first = text.splitlines()[0] if text else ""
    try:
        util_pct = float(first.strip())
    except ValueError:
        return {"util_pct": None, "raw": text,
                "error": f"unparseable: {first!r}"}
    return {"util_pct": util_pct, "raw": text, "error": ""}


def _check_local_free_leg_pid_alive(pid: int) -> bool:
    """Returns True if local PID is alive (signal-0 probe). PID <= 0
    skips the check (returns False — caller treats as 'free leg gone').
    """
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        # Exists but not ours (different uid). Treat as alive for safety.
        return True


def gate_gpu_host(
    gpu_host: str,
    free_leg_pid: int,
    refuse_threshold_pct: float = GPU_UTIL_REFUSE_THRESHOLD_PCT,
) -> Dict[str, Any]:
    """Pre-launch gate. Returns dict ``{"allow": bool, "reason": str,
    "platform": str, "cuda_visible_devices": str, "host": str,
    "checks": {...}}``.

    Modes:
      * ``vm``    — refuse if VM V100 util > refuse_threshold_pct
                    (Track A QM batch still active)
      * ``local`` — refuse if free_leg_pid is alive
                    (Track B free leg still active)
      * ``cpu``   — always allow (no GPU)
    """
    gpu_host = (gpu_host or "").strip().lower()
    if gpu_host == "cpu":
        return {
            "allow": True,
            "host": "cpu",
            "platform": "CPU",
            "cuda_visible_devices": "",
            "reason": "cpu mode requested (no GPU)",
            "checks": {},
        }
    if gpu_host == "vm":
        q = _query_vm_gpu_utilization_pct()
        if q["util_pct"] is None:
            return {
                "allow": False, "host": "vm", "platform": "CUDA",
                "cuda_visible_devices": "0",
                "reason": f"VM nvidia-smi probe failed: {q['error']}",
                "checks": {"vm_query": q},
            }
        if q["util_pct"] > refuse_threshold_pct:
            return {
                "allow": False, "host": "vm", "platform": "CUDA",
                "cuda_visible_devices": "0",
                "reason": (
                    f"VM V100 util {q['util_pct']:.0f}% > "
                    f"{refuse_threshold_pct:.0f}% — Track A QM batch "
                    f"still active. Refusing to share GPU."
                ),
                "checks": {"vm_query": q},
            }
        return {
            "allow": True, "host": "vm", "platform": "CUDA",
            "cuda_visible_devices": "0",
            "reason": f"VM V100 util {q['util_pct']:.0f}% <= threshold",
            "checks": {"vm_query": q},
        }
    if gpu_host == "local":
        alive = _check_local_free_leg_pid_alive(free_leg_pid)
        if alive:
            return {
                "allow": False, "host": "local", "platform": "CUDA",
                "cuda_visible_devices": "0",
                "reason": (
                    f"Free-leg PID {free_leg_pid} still alive on host "
                    f"5070Ti (device 0). Refusing to share GPU."
                ),
                "checks": {"free_leg_pid_alive": True},
            }
        return {
            "allow": True, "host": "local", "platform": "CUDA",
            "cuda_visible_devices": "0",
            "reason": (
                f"Free-leg PID {free_leg_pid} terminated; host 5070Ti free"
            ),
            "checks": {"free_leg_pid_alive": False},
        }
    return {
        "allow": False, "host": gpu_host, "platform": None,
        "cuda_visible_devices": None,
        "reason": f"unknown --gpu-host value: {gpu_host!r} "
                  f"(choose vm | local | cpu)",
        "checks": {},
    }


def normalize_legacy_cuda_device(value: Optional[str]) -> Optional[str]:
    """Backward compat: accepts legacy ``--cuda-device`` only when value is
    ``0`` or ``cpu``. Any other value (notably ``1`` — the old broken
    default) raises ValueError with a clear message.

    Returns the equivalent ``--gpu-host`` value (``local`` / ``cpu``) or
    ``None`` if the caller did not pass ``--cuda-device``.
    """
    if value is None or value == "":
        return None
    v = value.strip().lower()
    if v == "cpu":
        return "cpu"
    if v == "0":
        # Cannot disambiguate vm vs local from device index alone; instruct
        # caller to use --gpu-host explicitly.
        raise ValueError(
            "--cuda-device=0 is ambiguous (could be host 5070Ti or VM "
            "V100). Use --gpu-host {vm,local} explicitly."
        )
    raise ValueError(
        f"--cuda-device={value!r} is not a valid GPU on this project. "
        f"Both host (5070Ti) and VM (V100) expose only device 0. "
        f"Use --gpu-host {{vm,local,cpu}} instead. "
        f"(Legacy 'device 1' references in old docs were incorrect.)"
    )


# ---------------------------------------------------------------------------
# Top-level main()
# ---------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Track B per-direction structprep (Option B). "
            "PREP-ONLY: builds _dplus.xml + _dminus.xml per (endpoint, leg). "
            "Does NOT touch production. Free-leg PID 1876426 + V100 untouched."
        )
    )
    p.add_argument(
        "--v21-out-root",
        default="outputs/_trackb/production_v2_1",
        help="root of v2.1 systems (must contain <endpoint>/<leg>/trackb_sys.xml)",
    )
    p.add_argument(
        "--out-root",
        default="outputs/_trackb/per_direction_structprep",
        help="output root for prep artifacts + sanity audit",
    )
    p.add_argument(
        "--endpoints",
        default="cp4,wt",
        help="comma-list of endpoints",
    )
    p.add_argument(
        "--legs",
        default="bound,free",
        help="comma-list of legs",
    )
    p.add_argument(
        "--directions",
        default="dplus,dminus",
        help="comma-list of direction tags (dplus / dminus)",
    )
    p.add_argument(
        "--jobname",
        default="trackb",
        help="basename of per-leg cntl/pdb/xml",
    )
    p.add_argument(
        "--gpu-host",
        default="vm",
        choices=["vm", "local", "cpu"],
        help=(
            "GPU host (v0.9.10 device-guard fix). Each box has a single "
            "GPU at device 0; 'device 1' does not exist on either host. "
            "vm = ssh dispatch to VM V100 (refused if Track A util > "
            "30%%); local = host 5070Ti (refused if free-leg PID alive); "
            "cpu = OpenMM CPU platform (slow, audit-only safe)."
        ),
    )
    p.add_argument(
        "--cuda-device",
        default=None,
        help=(
            "DEPRECATED legacy option (v0.9.10 device-guard fix). Only "
            "'0' or 'cpu' accepted, and '0' is ambiguous between host "
            "and VM — prefer --gpu-host. Any other value (including "
            "the old default '1') is rejected: that GPU does not exist."
        ),
    )
    p.add_argument(
        "--gpu-util-refuse-threshold-pct",
        type=float,
        default=GPU_UTIL_REFUSE_THRESHOLD_PCT,
        help=(
            "VM V100 utilization threshold (%%) above which "
            "--gpu-host=vm is refused (Track A QM active). Default 30."
        ),
    )
    p.add_argument(
        "--skip-structprep",
        action="store_true",
        help="skip the per-direction structprep; only run audit on cached _0_*.xml",
    )
    p.add_argument(
        "--audit-only",
        action="store_true",
        help=(
            "alias of --skip-structprep — only do per-direction sanity audit + "
            "charge axis check (no GPU work, safe to run while free leg occupies 5070Ti)"
        ),
    )
    p.add_argument(
        "--smoke",
        action="store_true",
        help=(
            "smoke mode: respect _asyncre.cntl smoke step settings "
            "(THERMALIZATION_STEPS / ANNEALING_STEPS / EQUILIBRATION_STEPS already "
            "downsized by v2.1 launcher's --smoke; runs ~1-2 min per direction)"
        ),
    )
    p.add_argument(
        "--free-leg-pid",
        type=int,
        default=1876426,
        help="PID of free-leg run that must remain alive (don't touch). 0 = skip check.",
    )
    p.add_argument(
        "--charge-tol-e",
        type=float,
        default=5e-4,
        help="per-residue |Sigma q| tolerance (Q1(a) integer rescale gate)",
    )
    # v0.9.16 SSH-dispatch fix args
    p.add_argument(
        "--vm-ssh-host",
        default=_DEFAULT_VM_SSH_HOST,
        help=(
            "SSH host for --gpu-host=vm dispatch. Default "
            "san@192.168.122.155 (matches production launcher convention)."
        ),
    )
    p.add_argument(
        "--vm-python-bin",
        default=DEFAULT_VM_PYTHON_BIN,
        help=(
            "Python interpreter on VM for the structprep wrapper "
            "(must have atom_openmm + openmm installed; default "
            f"{DEFAULT_VM_PYTHON_BIN})."
        ),
    )
    p.add_argument(
        "--skip-vm-leg-dir-precheck",
        action="store_true",
        help=(
            "DEBUG: skip the VM leg-dir pre-flight gate. Operator confirms "
            "the VM leg_dir has the required inputs (trackb.pdb + _sys.xml "
            "+ _asyncre.cntl). Only used in tests."
        ),
    )
    p.add_argument(
        "--rsync-leg-inputs-to-vm",
        action="store_true",
        help=(
            "Task 2 (free-system V100 provisioning): before the VM leg-dir "
            "pre-flight gate, rsync each leg's structprep BUILD inputs "
            "(trackb.pdb + trackb_sys.xml + trackb_asyncre.cntl; plus the "
            "per-direction sys/pdb variants when present) from the host to "
            "the VM leg_dir. The cntl carries the densified38v4 schedule "
            "verbatim (schedule code is NOT needed VM-side). Only meaningful "
            "with --gpu-host=vm. Cohort-safe: any rsync failure halts before "
            "any structprep dispatch."
        ),
    )
    p.add_argument(
        "--rsync-outputs-from-vm",
        action="store_true",
        help=(
            "(A) VM->host rsync-BACK (v0.9.31): after each VM-dispatched "
            "structprep succeeds (rc=0), pull the GENUINE VM-produced "
            "per-direction base state(s) trackb_0_{dplus,dminus}.{xml,pdb} "
            "BACK to the host leg_dir so the host-side production readiness "
            "gate (check_free_pilot_readiness) and host-side subdir staging "
            "(stage_per_direction_subdir) find them (host FS is local ext4, "
            "NOT shared/NFS with the VM). FAILS LOUD if the VM output is "
            "absent — never fabricates a host base state (a faked base "
            "state would seed the wrong replicate velocities -> invalid "
            "sigma_btwn). Only meaningful with --gpu-host=vm."
        ),
    )
    p.add_argument(
        "--velocity-seed",
        type=int,
        default=None,
        help=(
            "Velocity-seed mechanism B: integer seed for "
            "setVelocitiesToTemperature in do_mintherm so this structprep "
            "draws INDEPENDENT initial velocities (per-replicate distinct "
            "integer 1..n → measurable sigma_btwn). Default None = upstream "
            "do_mintherm (no velocity seeding; single-run path unchanged). "
            "Distinct from the QM snapshot cohort {s7,...}: this is a "
            "velocity-init seed only. v0.9.30: plumbed through ALL dispatch "
            "modes (local | cpu | vm). For --gpu-host=vm the seed is baked "
            "into the ssh-dispatched wrapper as the VELOCITY_SEED constant "
            "and applied on the VM-side OpenMM Simulation (genuine VM "
            "seed-distinct base states; no silent drop)."
        ),
    )
    args = p.parse_args()

    # --------------------------------------------------------------
    # Hardware safety guards (conditions C1, C5; PID protection)
    # v0.9.10 device-guard fix: replace broken "device 1" assumption
    # with host-aware gate (vm | local | cpu).
    # --------------------------------------------------------------
    # Backward-compat: --cuda-device overrides --gpu-host iff it maps
    # cleanly to a {cpu} target; '0' is now ambiguous and raises.
    if args.cuda_device is not None:
        try:
            mapped = normalize_legacy_cuda_device(args.cuda_device)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        if mapped is not None:
            print(
                f"# --cuda-device={args.cuda_device!r} mapped to "
                f"--gpu-host={mapped!r} (legacy compat)"
            )
            args.gpu_host = mapped

    # Free leg PID alive snapshot (also used in report)
    free_leg_alive = _check_local_free_leg_pid_alive(args.free_leg_pid)
    print(f"# free leg PID {args.free_leg_pid} alive: {free_leg_alive}")

    # Host gate
    gate = gate_gpu_host(
        args.gpu_host,
        free_leg_pid=args.free_leg_pid,
        refuse_threshold_pct=args.gpu_util_refuse_threshold_pct,
    )
    print(f"# gpu_host={gate['host']!r} platform={gate['platform']!r} "
          f"allow={gate['allow']} — {gate['reason']}")
    if not gate["allow"]:
        print(f"ERROR: GPU host gate REFUSED: {gate['reason']}",
              file=sys.stderr)
        return 2
    if gate["host"] == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = gate["cuda_visible_devices"]
    platform_name = gate["platform"]
    print(f"# CUDA_VISIBLE_DEVICES: "
          f"{os.environ.get('CUDA_VISIBLE_DEVICES')!r}")

    # --------------------------------------------------------------
    # Setup
    # --------------------------------------------------------------
    v21_root = os.path.join(_PROJ_ROOT, args.v21_out_root)
    out_root = os.path.join(_PROJ_ROOT, args.out_root)
    os.makedirs(out_root, exist_ok=True)

    endpoints = [e.strip() for e in args.endpoints.split(",") if e.strip()]
    legs = [l.strip() for l in args.legs.split(",") if l.strip()]
    directions = [d.strip() for d in args.directions.split(",") if d.strip()]
    for e in endpoints:
        assert e in ("cp4", "wt"), f"unknown endpoint {e}"
    for l in legs:
        assert l in ("bound", "free"), f"unknown leg {l}"
    for d in directions:
        assert d in ("dplus", "dminus"), f"unknown direction tag {d}"
    direction_vals = {"dplus": 1, "dminus": -1}

    # --------------------------------------------------------------
    # Step 1: Per-direction structprep (if not --audit-only)
    # --------------------------------------------------------------
    structprep_results: List[Dict[str, Any]] = []
    vm_leg_dir_gates: List[Dict[str, Any]] = []  # populated only when --gpu-host=vm
    any_fail = False  # only meaningful when --gpu-host=vm cohort gate runs
    if not (args.skip_structprep or args.audit_only):
        # v0.9.30 Task 2: optionally provision the VM leg dirs with the
        # structprep build inputs (free-system V100 provisioning) BEFORE the
        # pre-flight gate, so the gate then passes. Cohort-safe: any rsync
        # failure raises RuntimeError -> halt before any structprep dispatch.
        if (gate["host"] == "vm" and args.rsync_leg_inputs_to_vm):
            print("\n# --- VM leg-input provisioning (Task 2 rsync) ---")
            for endpoint in endpoints:
                for leg in legs:
                    leg_dir = os.path.join(v21_root, endpoint, leg)
                    prov = _rsync_leg_inputs_to_vm(
                        leg_dir=leg_dir,
                        vm_ssh_host=args.vm_ssh_host,
                        jobname=args.jobname,
                    )
                    print(f"  {endpoint}/{leg}: {prov['status']} — "
                          f"pushed {len(prov['pushed'])} file(s)"
                          + (f", skipped absent {prov['skipped_absent']}"
                             if prov['skipped_absent'] else ""))

        # v0.9.16 SSH-dispatch fix: when --gpu-host=vm, run the F3-pattern
        # pre-flight gate per leg BEFORE any direction dispatch. Cohort-safe:
        # if even one leg fails the gate, halt the entire cohort rather than
        # half-stage. Matches the F3 gate semantics
        # (trackb_per_direction_production.py _gate_vm_abfe_bin_exists).
        if gate["host"] == "vm" and not args.skip_vm_leg_dir_precheck:
            print("\n# --- VM leg-dir pre-flight gate (F3-pattern) ---")
            any_fail = False
            for endpoint in endpoints:
                for leg in legs:
                    leg_dir = os.path.join(v21_root, endpoint, leg)
                    allow, reason = _gate_vm_leg_dir_exists(
                        leg_dir=leg_dir,
                        vm_ssh_host=args.vm_ssh_host,
                    )
                    vm_leg_dir_gates.append({
                        "endpoint": endpoint,
                        "leg": leg,
                        "leg_dir": leg_dir,
                        "allow": allow,
                        "reason": reason,
                    })
                    status_tag = "PASS" if allow else "FAIL"
                    print(f"  {endpoint}/{leg}: {status_tag} — {reason}")
                    if not allow:
                        any_fail = True
            if any_fail:
                print(
                    "\nERROR: One or more VM leg-dir pre-flight gates "
                    "FAILED. Halting cohort to avoid half-staged state. "
                    "Remediation: ensure leg_dir contents are present "
                    "on the VM (rsync if VM disk has space).",
                    file=sys.stderr,
                )
                # Still write the report so operator can see which legs
                # passed and which failed.
                structprep_results.append({
                    "status": "cohort_halted_vm_leg_dir_gate_fail",
                    "n_failed_gates": sum(
                        1 for g in vm_leg_dir_gates if not g["allow"]
                    ),
                    "gates": vm_leg_dir_gates,
                })

        if not any(
            r.get("status") == "cohort_halted_vm_leg_dir_gate_fail"
            for r in structprep_results
        ):
            for endpoint in endpoints:
                for leg in legs:
                    leg_dir = os.path.join(v21_root, endpoint, leg)
                    if not os.path.isdir(leg_dir):
                        print(f"WARN: missing v2.1 leg dir: {leg_dir}",
                              file=sys.stderr)
                        continue
                    for direction_tag in directions:
                        direction_val = direction_vals[direction_tag]
                        try:
                            r = run_per_direction_structprep(
                                leg_dir=leg_dir,
                                direction_val=direction_val,
                                jobname=args.jobname,
                                smoke=args.smoke,
                                gpu_host=gate["host"],
                                vm_ssh_host=args.vm_ssh_host,
                                vm_python_bin=args.vm_python_bin,
                                velocity_seed=args.velocity_seed,
                                rsync_outputs_from_vm=(
                                    args.rsync_outputs_from_vm),
                            )
                            r["endpoint"] = endpoint
                            r["leg"] = leg
                            r["direction_tag"] = direction_tag
                            structprep_results.append(r)
                            print(
                                f"   {endpoint}/{leg} d={direction_val:+d}: "
                                f"{r.get('status')} "
                                f"({r.get('size_bytes', 0):,} bytes, "
                                f"wall={r.get('wall_seconds', '-')}s)"
                            )
                        except Exception as e:
                            print(
                                f"ERROR: {endpoint}/{leg} "
                                f"d={direction_val:+d}: {e}",
                                file=sys.stderr,
                            )
                            structprep_results.append({
                                "endpoint": endpoint,
                                "leg": leg,
                                "direction_tag": direction_tag,
                                "status": "error",
                                "error": repr(e)[:300],
                            })
                    # C1 root-cause (v0.9.25): both directions' structprep
                    # shutil.move the upstream canonical {jobname}_0.{xml,pdb}
                    # to per-direction variant names, so no canonical
                    # {jobname}_0.xml survives. The abfe_production worker
                    # (ommworker.py:265) hardcodes loadState({jobname}_0.xml)
                    # for the service worker (compute=False template only —
                    # per-direction base state restore).
                    # Restore the canonical from the dplus variant (copy,
                    # idempotent, no-overwrite). Covers the host-local
                    # path; the SSH path's VM wrapper restores its own VM-side
                    # canonical inline.
                    try:
                        cr = _restore_canonical_from_dplus(
                            leg_dir, jobname=args.jobname,
                        )
                        if cr.get("created_xml"):
                            print(f"   {endpoint}/{leg}: restored canonical "
                                  f"{args.jobname}_0.xml from dplus variant")
                    except Exception as e:
                        print(f"WARN: canonical restore {endpoint}/{leg}: {e}",
                              file=sys.stderr)

    # --------------------------------------------------------------
    # Step 2: Per-direction sanity audit
    # --------------------------------------------------------------
    sanity_results: List[Dict[str, Any]] = []
    charge_results: List[Dict[str, Any]] = []
    for endpoint in endpoints:
        for leg in legs:
            leg_dir = os.path.join(v21_root, endpoint, leg)
            if not os.path.isdir(leg_dir):
                continue
            try:
                sr = audit_per_direction_xml(leg_dir, jobname=args.jobname)
                sr["endpoint"] = endpoint
                sr["leg"] = leg
                sanity_results.append(sr)
            except Exception as e:
                sanity_results.append({
                    "endpoint": endpoint, "leg": leg,
                    "error": repr(e)[:300],
                })
            try:
                cr = audit_charge_axis_persistence(
                    leg_dir, jobname=args.jobname,
                    tol_e=args.charge_tol_e,
                )
                cr["endpoint"] = endpoint
                cr["leg"] = leg
                charge_results.append(cr)
            except Exception as e:
                charge_results.append({
                    "endpoint": endpoint, "leg": leg,
                    "error": repr(e)[:300],
                })

    # --------------------------------------------------------------
    # Aggregate report
    # --------------------------------------------------------------
    report = {
        "regime": "ranking_only",
        "track": "B",
        "method": "Option B per-direction structprep (2026-05-31)",
        "method_ref": "per-direction structprep (Option B), 2026-05-31",
        "v21_out_root": v21_root,
        "out_root": out_root,
        "endpoints": endpoints,
        "legs": legs,
        "directions": directions,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "platform": platform_name,
        "gpu_host": gate["host"],
        "gpu_host_reason": gate["reason"],
        "gpu_host_checks": gate["checks"],
        "free_leg_pid_alive": free_leg_alive,
        "free_leg_pid": args.free_leg_pid,
        "vm_ssh_host": args.vm_ssh_host,
        "vm_python_bin": args.vm_python_bin,
        "vm_leg_dir_gates": vm_leg_dir_gates,
        "rsync_leg_inputs_to_vm": bool(args.rsync_leg_inputs_to_vm),
        "rsync_outputs_from_vm": bool(args.rsync_outputs_from_vm),
        "velocity_seed": args.velocity_seed,
        "structprep_results": structprep_results,
        "sanity_audit": sanity_results,
        "charge_axis_audit": charge_results,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "schema_version": (
            "trackb_per_direction_structprep_v4_vm_rsync_back"
        ),
        "production_launch_blocker": (
            "Production launch BLOCKED — requires (a) free leg PID complete + "
            "results validated, (b) explicit user confirm, (c) "
            "C1-C8 dry-run gate PASS. See scripts/trackb_per_direction_production.py"
        ),
    }
    ts = time.strftime("%Y%m%dT%H%M%S")
    report_path = os.path.join(out_root, f"per_direction_prep_{ts}.json")
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    latest = os.path.join(out_root, "per_direction_prep_latest.json")
    try:
        if os.path.islink(latest) or os.path.exists(latest):
            os.remove(latest)
        os.symlink(os.path.basename(report_path), latest)
    except OSError:
        pass

    print(f"\n=== REPORT ===")
    print(f"  structprep results: {len(structprep_results)}")
    print(f"  sanity audit:       {len(sanity_results)}")
    print(f"  charge axis audit:  {len(charge_results)}")
    print(f"  report: {report_path}")

    # Final acceptance gate summary
    # Charge gate: require non-empty audit set + all PASS (vacuous-true guard)
    charge_audits_with_data = [cr for cr in charge_results if "error" not in cr]
    all_charge_pass = (
        len(charge_audits_with_data) > 0
        and all(cr.get("all_residues_within_tol", False)
                for cr in charge_audits_with_data)
    )
    # Cyclic gate: require non-empty audit set + all PASS + all legs have
    # both _dplus.xml and _dminus.xml produced (no missing-direction false PASS)
    cyclic_audits_with_both_dir = [
        sr for sr in sanity_results
        if "error" not in sr
        and "cyclic_ss_intact_both_directions" in sr
        and all(d.get("status") == "ok"
                for d in sr.get("directions", {}).values())
    ]
    all_cyclic_pass = (
        len(cyclic_audits_with_both_dir) == len(sanity_results)
        and len(cyclic_audits_with_both_dir) > 0
        and all(sr["cyclic_ss_intact_both_directions"]
                for sr in cyclic_audits_with_both_dir)
    )
    # Status label includes N/A when XMLs missing (audit-only mode)
    if len(cyclic_audits_with_both_dir) == 0:
        cyclic_label = "N/A (per-direction XMLs not yet produced)"
    elif len(cyclic_audits_with_both_dir) < len(sanity_results):
        cyclic_label = (
            f"PARTIAL ({len(cyclic_audits_with_both_dir)}/"
            f"{len(sanity_results)} legs)"
        )
    else:
        cyclic_label = "PASS" if all_cyclic_pass else "FAIL"
    print(f"  Q1(a) charge persistence (all legs): "
          f"{'PASS' if all_charge_pass else 'FAIL/N/A'}")
    print(f"  cyclic_ss intact (both directions, all legs): "
          f"{cyclic_label}")
    print(
        f"\n*** PRODUCTION LAUNCH BLOCKED — requires free leg complete + "
        f"user confirm + dry-run gate PASS. See "
        f"scripts/trackb_per_direction_production.py ***"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
