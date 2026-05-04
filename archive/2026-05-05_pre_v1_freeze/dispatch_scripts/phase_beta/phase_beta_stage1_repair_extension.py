#!/usr/bin/env python
"""
2 ns extension on 7TL8_MTR6_s7 final structure.

Strategy:
  - Load 7TL8_MTR6_final.pdb (already-solvated end-state from the 5 ns run).
  - Use the same XML stack (amber14SB protein + tip3pfb water + MTR_gaff2.xml).
  - Minimal equilibration (1000 minimization steps, 5 ps NVT, 5 ps NPT to relax velocities).
  - 2 ns production at dt=2 fs, save DCD every 5000 steps (10 ps → 200 frames).

Output: outputs/7TL8_MTR6_calib_s7/mdresult/extension_audit_20260428/extension_2ns.dcd
"""
import os
import sys
import time
from pathlib import Path

import openmm as mm
import openmm.unit as unit
from openmm import app

ROOT = Path('/home/san/UPDD_proj')
SYS_DIR = ROOT / 'outputs' / '7TL8_MTR6_calib_s7'
EXT_DIR = SYS_DIR / 'mdresult' / 'extension_audit_20260428'
EXT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_PDB = SYS_DIR / 'mdresult' / '7TL8_MTR6_final.pdb'
MTR_FF_XML = SYS_DIR / 'params' / 'MTR_gaff2.xml'
MTR_HYD_XML = SYS_DIR / 'params' / 'MTR_hydrogens.xml'

OUT_DCD = EXT_DIR / 'extension_2ns.dcd'
OUT_LOG = EXT_DIR / 'extension_2ns.log'
OUT_FINAL = EXT_DIR / 'extension_2ns_final.pdb'

# ---------- Parameters ----------
DT_FS = 2.0
TEMP_K = 300.0
N_PROD_STEPS = 1_000_000  # 2 ns at 2 fs
N_NVT_EQ_STEPS = 2500     # 5 ps
N_NPT_EQ_STEPS = 2500     # 5 ps
N_MIN_STEPS = 1000
DCD_INTERVAL = 5000       # 10 ps  → 200 frames over 2 ns
SEED = 314159

print(f"[ext] PDB={INPUT_PDB} ({INPUT_PDB.stat().st_size/1e6:.1f} MB)")
print(f"[ext] MTR FF={MTR_FF_XML}")
print(f"[ext] OUT  ={OUT_DCD}")

# ---------- Load PDB ----------
print("[ext] loading PDB ...")
pdb = app.PDBFile(str(INPUT_PDB))
print(f"[ext] n_atoms={pdb.topology.getNumAtoms()}, "
      f"n_residues={pdb.topology.getNumResidues()}, "
      f"box={pdb.topology.getPeriodicBoxVectors()}")

# Register MTR Hydrogens definitions for Modeller compatibility (not strictly
# needed since the input PDB already has H, but harmless).
app.Modeller.loadHydrogenDefinitions(str(MTR_HYD_XML))

# ---------- ForceField ----------
print("[ext] building ForceField ...")
ff = app.ForceField(
    'amber14-all.xml',
    'amber14/tip3pfb.xml',
    str(MTR_FF_XML),
)

print("[ext] creating system ...")
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
)

# Add MonteCarloBarostat for NPT
print("[ext] adding MonteCarloBarostat (1 atm, 300 K) ...")
system.addForce(mm.MonteCarloBarostat(1.0 * unit.atmospheres,
                                      TEMP_K * unit.kelvin, 25))

# ---------- Integrator + Platform ----------
integrator = mm.LangevinMiddleIntegrator(
    TEMP_K * unit.kelvin,
    1.0 / unit.picosecond,
    DT_FS * unit.femtoseconds,
)
integrator.setRandomNumberSeed(SEED)

plat = mm.Platform.getPlatformByName('CUDA')
props = {'Precision': 'mixed'}
sim = app.Simulation(pdb.topology, system, integrator, plat, props)
sim.context.setPositions(pdb.positions)
print(f"[ext] platform={plat.getName()}")

# ---------- Minimize ----------
print("[ext] minimizing (1000 iter) ...")
t0 = time.time()
sim.minimizeEnergy(maxIterations=N_MIN_STEPS)
print(f"[ext] min done ({time.time()-t0:.1f}s)")

# ---------- Set velocities @ 300 K ----------
sim.context.setVelocitiesToTemperature(TEMP_K * unit.kelvin, SEED)

# ---------- NVT equilibration (5 ps, no barostat for this segment) ----------
# Disable barostat by setting frequency very high (skip)
# Actually leaving barostat on; combined NPT is fine for short eq.
print(f"[ext] equilibration {N_NVT_EQ_STEPS + N_NPT_EQ_STEPS} steps "
      f"({(N_NVT_EQ_STEPS+N_NPT_EQ_STEPS)*DT_FS/1000:.0f} ps) ...")
t0 = time.time()
sim.step(N_NVT_EQ_STEPS + N_NPT_EQ_STEPS)
print(f"[ext] equilibration done ({time.time()-t0:.1f}s)")

# ---------- Production: DCD reporter + State data ----------
print(f"[ext] production: {N_PROD_STEPS} steps "
      f"({N_PROD_STEPS*DT_FS/1000:.0f} ps), DCD every {DCD_INTERVAL} steps "
      f"({DCD_INTERVAL*DT_FS/1000:.0f} ps) → "
      f"{N_PROD_STEPS//DCD_INTERVAL} frames")

sim.reporters.append(app.DCDReporter(str(OUT_DCD), DCD_INTERVAL))
sim.reporters.append(app.StateDataReporter(
    str(OUT_LOG), DCD_INTERVAL, step=True, time=True,
    potentialEnergy=True, kineticEnergy=True, temperature=True,
    volume=True, speed=True, separator='\t',
))

t0 = time.time()
sim.step(N_PROD_STEPS)
prod_time = time.time() - t0
print(f"[ext] production done ({prod_time:.1f}s = {prod_time/60:.1f} min)")

# ---------- Save final state ----------
state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
positions = state.getPositions()
with open(OUT_FINAL, 'w') as f:
    app.PDBFile.writeFile(pdb.topology, positions, f, keepIds=True)
print(f"[ext] final PDB → {OUT_FINAL}")

print(f"[ext] DONE. DCD={OUT_DCD} ({OUT_DCD.stat().st_size/1e6:.1f} MB)")
