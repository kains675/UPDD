#!/usr/bin/env python
"""
sapt_dimer_benchmark.py
=======================
First-principles gas-phase model-dimer benchmark to diagnose the accuracy of the
UPDD force-field (FF) *interaction energy* for ncAA substitutions, by triangulating

        FF  <->  wB97X-D  <->  FNO-CCSD(T)          (all on identical geometry)

and attributing any FF error to a physical SAPT component (elst/exch/ind/disp).

It answers the methods-audit fork "FF-accuracy is orthogonal to sampling": for a
fixed crystal geometry, is the ~1 kcal/mol methylation effect (MTR = 1-Me-Trp vs
WT Trp) resolvable above chemical accuracy AND decomposable by QM, and if so which
FF term is wrong?

scientific design gate (GO-WITH-COND):
    internal scientific-review note 20260703

--------------------------------------------------------------------------------
GEOMETRY SSOT (condition C1)
--------------------------------------------------------------------------------
Crystal `target/2QKI_clean.pdb` (compstatin.C3c, Janssen 2007 JBC 282:29241).
    chain B TRP4 = WT bound pose (contacts intact).
    chain A = C3c target; dominant partner = HIS383 (HID neutral).
MD final frames are FORBIDDEN as geometry source: the MTR4 indole is disengaged
(0 chain-A contacts <=4.5A) in equilibrated snapshots (pose-disengagement pathology).
MTR is built by VERTICAL methylation of the fixed crystal WT heavy frame (only the
added methyl H's are idealized) so the diagnostic isolates the pure electronic
interaction-energy effect of N1-methylation.

--------------------------------------------------------------------------------
FRAGMENTS (condition C2)
--------------------------------------------------------------------------------
Monomer A (the ncAA sidechain analog, built on the crystal Trp4 frame):
    WT  = 3-methylindole (skatole)          : Trp sidechain, Ca-Cb bond H-capped
    MTR = 1,3-dimethylindole                 : + N1-CH3 (vertical, replaces HE1)
    W4A = methane                            : indole removed -> CH4 at Cb (control)
Monomer B (the pocket):
    Tier-1 = capped His383 (Ca-H + imidazole, HID neutral)     ~15 atoms
    Tier-2 = Tier-1 + Thr382 backbone C=O(+Ca) + Gly336 backbone ~27 atoms
    (Thr382 C=O is 3.39A from NE1 = the methyl-attach electrostatic channel.)

--------------------------------------------------------------------------------
METHOD STACK (condition C3), identical geometry, per fragment pair
--------------------------------------------------------------------------------
    1. SAPT0/jun-cc-pVDZ                 -> elst/exch/ind/disp decomposition (BSSE-free)
    2. SAPT2+(3)dMP2/aug-cc-pVDZ         -> refined decomposition + near-CCSD(T) total
    3. FNO-CCSD(T)/aug-cc-pVTZ + CP      -> GOLD total  (Tier-1 only; focal-point fallback)
    4. wB97X-D/def2-TZVP + CP            -> UPDD production DFT layer
    5. FF nonbonded (ff14SB + MTR GAFF2 hybrid) -> the layer under test
DLPNO-CCSD(T) is unavailable in psi4 -> gold = FNO-CCSD(T) (frozen NO + DF).
Counterpoise (Boys-Bernardi) applied to all supermolecular (CCSD(T)/wB97X-D).
SAPT is BSSE-free by construction. FF has no BSSE.

--------------------------------------------------------------------------------
PRE-REGISTERED 3-OUTCOME (condition C4, anti-HARKing) -- see classify_outcome()
--------------------------------------------------------------------------------
(i)  FF-FIXABLE          : |ddE_int^CCSD(T)| >= ~1 kcal/mol (above chem accuracy)
                           AND SAPT attributes ddE to a component
                           AND |ddE^FF - ddE^CCSD(T)| resolvable in that component
                           => "FF term X is wrong and calibratable".
(ii) RANKING-ONLY        : |ddE_int^CCSD(T)| <~ chemical accuracy AND buried in the
                           method spread (CCSD(T) vs SAPT2+(3), basis incompleteness,
                           Tier-1 vs Tier-2) => "not QM-decomposable even in a fixed
                           gas-phase dimer; ranking-only is justified".
(iii)FF STRUCTURAL ERROR : sign(ddE^FF) != sign(ddE^CCSD(T)) OR magnitude off by
                           >2-3x AND SAPT names a physical component FF omits
                           => "FF structural error" (links to MM-PBSA wrong-sign history).
W4A positive control: |ddE_int^CCSD(T)| MUST be large & unfavorable (indole.pocket
contact lost). If it is not, the method/geometry is broken (control FAIL) and the
MTR result cannot be trusted.

--------------------------------------------------------------------------------
R-11 / R-18 CAVEATS (condition C5) -- this is NOT a binding-affinity prediction
--------------------------------------------------------------------------------
The reported quantity is ddE_int: a GAS-PHASE, FIXED-GEOMETRY, 2-BODY interaction
energy difference. It is NOT ddG_bind. Magotti (2009, J Mol Recognit 22:495) ITC
-1.4 / SPR -3.0 kcal/mol are a *different* quantity (solution free energy); any sign
comparison is SUGGESTIVE ONLY, never validation. Model truncations (2-body pocket vs
full multi-residue site, absent solvation, absent geometry relaxation, link/methyl-H
charge approximation) are stated limitations; if the partner fragment is inaccurate
the result is indeterminate.

--------------------------------------------------------------------------------
EXECUTION (condition C4) -- env-split orchestration
--------------------------------------------------------------------------------
This single file runs in three roles (top level imports only numpy + stdlib):
    orchestrator      : build fragments, dispatch sub-jobs, aggregate.  (any env w/ numpy)
    --_qm-driver JSON : psi4 methods.  Re-invoked with the psi4 env python.
    --_ff-driver JSON : OpenMM nonbonded.  Re-invoked with the qmmm env python.
psi4 is CPU-only and fully independent of the GPU (Gate-A MD runs concurrently).
Full method-stack execution (~1-2 days wall) is the execution review's job.

DOIs
----
SAPT levels/basis : Parker, Burns, Sherrill 2014 JCP 140:094106  10.1063/1.4867135
FNO-CCSD(T)       : DePrince & Sherrill 2013 JCTC 9:293          10.1021/ct300780u
Counterpoise      : Boys & Bernardi 1970 Mol Phys 19:553         10.1080/00268977000101561
wB97X-D           : Chai & Head-Gordon 2008 PCCP 10:6615         10.1039/b810189b
Link-atom H-cap   : Senn & Thiel 2009 Angew Chem Int Ed 48:1198  (reused convention)
Magotti (SIGN)    : 2009 J Mol Recognit 22:495                   10.1002/jmr.972
2QKI structure    : Janssen 2007 JBC 282:29241
"""

import os
import sys
import json
import math
import argparse
import subprocess

import numpy as np

# numpy 2.x compat shim: qcelemental/models/molecule.py:369 calls
# np.core.defchararray.title(...), a path removed in numpy >= 2.4 (raises
# "module 'numpy._core' has no attribute 'defchararray'"). It fires when the
# gold FNO-CCSD(T)/wB97X-D drivers build a qcelemental Molecule (SAPT uses the
# psi4-native Molecule and is unaffected). Alias the legacy submodule to the
# public np.char (equivalent element-wise title-case) so the dependency is left
# untouched but the driver works. Runs at import so the QM sub-process inherits it.
try:
    _npcore = getattr(np, "_core", None) or np.core
    if not hasattr(_npcore, "defchararray"):
        import numpy.char as _npchar
        _npcore.defchararray = _npchar
except Exception:
    pass

# --------------------------------------------------------------------------
# Paths / environment interpreters (verdict Q5 / C4)
# --------------------------------------------------------------------------
_THIS = os.path.abspath(__file__)
REPO = os.path.dirname(os.path.dirname(_THIS))

CRYSTAL_PDB = os.path.join(REPO, "target", "2QKI_clean.pdb")
MTR_XML = os.path.join(REPO, "params", "MTR_gaff2_hybrid.xml")

PSI4_PY = "/home/san/miniconda3/envs/psi4/bin/python"
QMMM_PY = "/home/san/miniconda3/envs/qmmm/bin/python"

# psi4 resource limits (verdict Q5): scratch -> dedicated SSD, threads <= 14.
PSI_SCRATCH_DIR = os.environ.get("SAPT_PSI_SCRATCH", "/media/san/San")
PSI_MEMORY = os.environ.get("SAPT_PSI_MEM", "50 GB")  # env-override (2026-07-05: WT/MTR gold crashed at 50GB — SIGSEGV/DSYEV from un-managed LAPACK alloc pushing 60GB host; lower to leave headroom)
PSI_THREADS = 14

HARTREE2KCAL = 627.5094740631  # CODATA-consistent Hartree -> kcal/mol

# --------------------------------------------------------------------------
# Fragment definitions.  Atom names refer to the crystal TRP4 (chain B) and
# the chain-A partner residues.  ff_source picks the residue template used for
# FF charge/LJ lookup; caps (role="cap") are geometric H's whose charge is a
# per-fragment neutralization sink.
# --------------------------------------------------------------------------
BINDER_CHAIN = "B"
TARGET_CHAIN = "A"
NCAA_RESI = "4"          # Trp4 in compstatin numbering (crystal)
HIS_RESI = "383"
THR_RESI = "382"
GLY_RESI = "336"

# TRP sidechain heavy+H names kept for the 3-methylindole (skatole) analog.
_TRP_SIDECHAIN = [
    "CB", "HB2", "HB3", "CG", "CD1", "HD1", "NE1", "HE1",
    "CE2", "CZ2", "HZ2", "CH2", "HH2", "CZ3", "HZ3", "CE3", "HE3", "CD2",
]
# MTR (1-Me-Trp) analog = skatole minus NE1-H (HE1) plus the N-methyl (CM,HM*).
_MTR_SIDECHAIN_REAL = [a for a in _TRP_SIDECHAIN if a != "HE1"]  # HE1 replaced by CM

# Capped His383 (HID): Ca-H + imidazole; backbone N/C -> caps.
_HIS_KEEP = ["CA", "HA", "CB", "HB2", "HB3", "CG",
             "ND1", "HD1", "CD2", "HD2", "CE1", "HE1", "NE2"]
# Thr382 backbone C=O + Ca (Tier-2); Gly336 backbone (Tier-2).
_THR_KEEP = ["CA", "HA", "C", "O"]
_GLY_KEEP = ["CA", "HA2", "HA3", "C", "O"]

# Bond lengths (Angstrom).
BOND_CH = 1.09     # C-H / cap H, Senn & Thiel 2009 (matches run_qmmm._make_link_h)
BOND_NCH3 = 1.47   # aromatic N-CH3


# ==========================================================================
# Geometry primitives
# ==========================================================================
def _cap_h_position(p_keep, p_remove, bond_len=BOND_CH):
    """H-cap coordinate: place an H on the kept atom along the cut bond direction.

    Reuses verbatim the geometry convention of
    ``utils/run_qmmm.py::_make_link_h`` (Senn & Thiel 2009 QM/MM link atom):
        h = p_keep + unit(p_remove - p_keep) * bond_len
    Reimplemented here (not imported) so the fragment builder stays numpy-only
    and env-agnostic; run_qmmm imports PySCF/OpenMM at module load and would
    couple this orchestrator to the qmmm env.  The frozen core is untouched.
    """
    p_keep = np.asarray(p_keep, float)
    p_remove = np.asarray(p_remove, float)
    v = p_remove - p_keep
    n = np.linalg.norm(v)
    if n < 1e-6:
        v = np.array([1.0, 0.0, 0.0])
        n = 1.0
    return p_keep + (v / n) * bond_len


def _place_methyl_hydrogens(c_pos, x_pos, ref_pos, ch=BOND_CH):
    """Idealized sp3 methyl H's on carbon ``c_pos`` bonded to heavy ``x_pos``.

    3 tetrahedral H's (X-C-H = 109.47 deg, 120 deg apart), staggered relative to
    ``ref_pos`` (a neighbor of X) for deterministic placement.  Only the methyl
    H's are idealized; the heavy frame stays at crystal coordinates (verdict Q2).
    """
    c_pos = np.asarray(c_pos, float)
    x_pos = np.asarray(x_pos, float)
    a = c_pos - x_pos
    a /= np.linalg.norm(a)                      # axis X->C ; H's have h.a = +1/3
    ref = np.asarray(ref_pos, float) - c_pos
    ref_perp = ref - np.dot(ref, a) * a
    if np.linalg.norm(ref_perp) < 1e-6:         # ref colinear with axis -> arbitrary
        seed = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(seed, a)) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
        ref_perp = seed - np.dot(seed, a) * a
    u = -ref_perp / np.linalg.norm(ref_perp)    # first H staggered away from ref
    v = np.cross(a, u)
    perp = math.sqrt(8.0) / 3.0                 # sin(arccos(1/3))
    hs = []
    for k in range(3):
        theta = k * (2.0 * math.pi / 3.0)
        direction = (1.0 / 3.0) * a + perp * (math.cos(theta) * u + math.sin(theta) * v)
        hs.append(c_pos + ch * direction)       # |direction| == 1
    return hs


# ==========================================================================
# Minimal PDB reader (matches the scientific review reference scripts' column parse)
# ==========================================================================
def parse_pdb(fn):
    atoms = []
    with open(fn) as f:
        for l in f:
            if l.startswith(("ATOM", "HETATM")):
                atoms.append(dict(
                    name=l[12:16].strip(), resn=l[17:20].strip(), chain=l[21],
                    resi=l[22:26].strip(),
                    x=float(l[30:38]), y=float(l[38:46]), z=float(l[46:54]),
                ))
    return atoms


def _index(atoms, chain, resi):
    """{atom_name: np.array([x,y,z])} for one residue."""
    return {a["name"]: np.array([a["x"], a["y"], a["z"]])
            for a in atoms if a["chain"] == chain and a["resi"] == resi}


def _element_of(name):
    """PDB atom name -> element (aromatic/aliphatic naming; H first)."""
    n = name.strip()
    if not n:
        return "C"
    if n[0].isdigit():
        n = n[1:]
    return "H" if n[:1] == "H" else n[:1].upper()


def _atom(name, xyz, element, ff_source, ff_name, role="real"):
    return {
        "name": name,
        "xyz": [float(xyz[0]), float(xyz[1]), float(xyz[2])],
        "element": element,
        "ff_source": ff_source,   # TRP/MTR/HID/THR/GLY/ALA/CAP
        "ff_name": ff_name,       # name in the source residue template (None for caps)
        "role": role,             # "real" | "cap"
    }


# ==========================================================================
# Monomer A (ncAA sidechain analog) builders
# ==========================================================================
def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v


def build_monomer_A(variant, trp):
    """Return atom list for the ncAA sidechain analog on the crystal Trp4 frame.

    variant: "WT" (3-methylindole), "MTR" (1,3-dimethylindole), "W4A" (methane).
    trp: {atom_name: xyz} of crystal TRP4.
    """
    out = []
    if variant == "W4A":
        # methane at Cb: keep Cb,HB2,HB3 ; cap the Ca and CG bonds -> CH4.
        out.append(_atom("CB", trp["CB"], "C", "ALA", "CB"))
        out.append(_atom("HB2", trp["HB2"], "H", "ALA", "HB2"))
        out.append(_atom("HB3", trp["HB3"], "H", "ALA", "HB3"))
        out.append(_atom("HC1", _cap_h_position(trp["CB"], trp["CA"]), "H", "CAP", None, "cap"))
        out.append(_atom("HC2", _cap_h_position(trp["CB"], trp["CG"]), "H", "CAP", None, "cap"))
        return out

    ff_source = "MTR" if variant == "MTR" else "TRP"
    keep = _MTR_SIDECHAIN_REAL if variant == "MTR" else _TRP_SIDECHAIN
    for nm in keep:
        out.append(_atom(nm, trp[nm], _element_of(nm), ff_source, nm))
    # Ca-Cb cut -> H-cap (skatole methyl completion).
    out.append(_atom("HC1", _cap_h_position(trp["CB"], trp["CA"]), "H", "CAP", None, "cap"))

    if variant == "MTR":
        # Vertical N1-methylation: CM along NE1->HE1 (crystal) at 1.47A; 3 idealized H's.
        cm = np.asarray(trp["NE1"], float) + BOND_NCH3 * _unit(np.asarray(trp["HE1"], float)
                                                              - np.asarray(trp["NE1"], float))
        out.append(_atom("CM", cm, "C", "MTR", "CM"))
        for i, h in enumerate(_place_methyl_hydrogens(cm, trp["NE1"], trp["CD1"]), start=1):
            out.append(_atom("HM%d" % i, h, "H", "MTR", "HM%d" % i))
    return out


# ==========================================================================
# Monomer B (pocket) builders
# ==========================================================================
def _capped_residue(res_atoms, keep_names, ff_source, cut_pairs):
    """Keep ``keep_names`` from one residue; add an H-cap for each (keep, remove) cut."""
    out = []
    for nm in keep_names:
        if nm not in res_atoms:
            continue
        out.append(_atom(nm, res_atoms[nm], _element_of(nm), ff_source, nm))
    for i, (keep, remove) in enumerate(cut_pairs, start=1):
        if keep in res_atoms and remove in res_atoms:
            hc = _cap_h_position(res_atoms[keep], res_atoms[remove])
            out.append(_atom("H%s%d" % (ff_source[:2], i), hc, "H", "CAP", None, "cap"))
    return out


def build_monomer_B(tier, atoms):
    """Pocket fragment.  tier=1 -> capped His383 ; tier=2 -> + Thr382 C=O(+Ca) + Gly336."""
    his = _index(atoms, TARGET_CHAIN, HIS_RESI)
    out = _capped_residue(his, _HIS_KEEP, "HID", [("CA", "N"), ("CA", "C")])
    if tier == 2:
        thr = _index(atoms, TARGET_CHAIN, THR_RESI)
        gly = _index(atoms, TARGET_CHAIN, GLY_RESI)
        out += _capped_residue(thr, _THR_KEEP, "THR", [("CA", "N"), ("CA", "CB")])
        out += _capped_residue(gly, _GLY_KEEP, "GLY", [("CA", "N")])
    return out


# ==========================================================================
# Fragment pair assembly
# ==========================================================================
def build_pair(variant, tier, atoms):
    trp = _index(atoms, BINDER_CHAIN, NCAA_RESI)
    missing = [n for n in ("CA", "CB", "CG", "NE1", "HE1", "HB2", "HB3") if n not in trp]
    if missing:
        raise RuntimeError("crystal TRP%s missing atoms %s in %s"
                           % (NCAA_RESI, missing, CRYSTAL_PDB))
    A = build_monomer_A(variant, trp)
    B = build_monomer_B(tier, atoms)
    return {"variant": variant, "tier": tier, "A": A, "B": B}


def pair_key(variant, tier):
    return "%s_tier%d" % (variant, tier)


# ==========================================================================
# QM driver (runs under the psi4 env python)
# ==========================================================================
QM_METHODS = {
    # name: (psi4 method, basis, mode)   mode = "sapt" | "cp"
    "sapt0": ("sapt0", "jun-cc-pvdz", "sapt"),
    "sapt2p3_dmp2": ("sapt2+(3)dmp2", "aug-cc-pvdz", "sapt"),
    "fno_ccsd_t": ("fno-ccsd(t)", "aug-cc-pvdz", "cp"),      # Tier-1 only (gold); aVTZ->aVDZ 2026-07-03 (aVTZ CCSD(T) SIGSEGV'd, scientifically approved fallback)
    "wb97xd": ("wb97x-d", "def2-tzvp", "cp"),
}
# Tier-2 is SAPT-only (SAPT0 + SAPT2+(3)) per verdict (CCSD(T)/aVTZ too heavy).
TIER2_METHODS = ("sapt0", "sapt2p3_dmp2")


def _geom_block(frag):
    return "\n".join("%2s %18.10f %18.10f %18.10f"
                     % (a["element"], a["xyz"][0], a["xyz"][1], a["xyz"][2])
                     for a in frag)


def run_qm_driver(job_path):
    """Execute the QM method stack for one fragment pair.  Writes <job>.qm.json.

    Robustness: each method is isolated in try/except so one failure (e.g. an
    aVTZ CCSD(T) OOM) does not lose the cheaper results; partial JSON is flushed
    after every method.
    """
    import psi4

    with open(job_path) as f:
        job = json.load(f)

    workdir = job["workdir"]
    scratch = job.get("psi_scratch") or PSI_SCRATCH_DIR
    if os.path.isdir(scratch):
        os.environ.setdefault("PSI_SCRATCH", scratch)
    psi4.set_memory(job.get("memory", PSI_MEMORY))
    psi4.set_num_threads(int(job.get("threads", PSI_THREADS)))
    out_log = os.path.join(workdir, "%s.psi4.log" % pair_key(job["variant"], job["tier"]))
    psi4.core.set_output_file(out_log, False)

    methods = job.get("methods")
    if not methods:
        methods = list(QM_METHODS) if job["tier"] == 1 else list(TIER2_METHODS)

    geom = "0 1\n%s\n--\n0 1\n%s\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n" % (
        _geom_block(job["A"]), _geom_block(job["B"]))

    results = {"variant": job["variant"], "tier": job["tier"],
               "n_atoms_A": len(job["A"]), "n_atoms_B": len(job["B"]), "methods": {}}
    out_json = os.path.splitext(job_path)[0] + ".qm.json"

    # Merge with any prior run's results so a targeted re-run (e.g. only the
    # gold CCSD(T) after a disk-full crash) preserves already-computed methods.
    if os.path.isfile(out_json):
        try:
            _prior = json.load(open(out_json))
            if _prior.get("variant") == job["variant"] and _prior.get("tier") == job["tier"]:
                results["methods"].update(_prior.get("methods", {}))
        except Exception:
            pass

    def _flush():
        with open(out_json, "w") as f:
            json.dump(results, f, indent=2)

    for m in methods:
        if m not in QM_METHODS:
            continue
        method, basis, mode = QM_METHODS[m]
        if m == "fno_ccsd_t" and job["tier"] != 1:
            continue
        try:
            psi4.core.clean()
            psi4.core.clean_options()
            psi4.set_memory(job.get("memory", PSI_MEMORY))
            psi4.set_num_threads(int(job.get("threads", PSI_THREADS)))
            mol = psi4.geometry(geom)
            basis_use = job.get("basis_override", {}).get(m, basis)
            opts = {"basis": basis_use, "freeze_core": "true"}
            if mode == "sapt":
                opts["scf_type"] = "df"
                psi4.set_options(opts)
                psi4.energy(method, molecule=mol)
                entry = {
                    "method": method, "basis": basis_use, "mode": mode,
                    "elst_kcal": psi4.variable("SAPT ELST ENERGY") * HARTREE2KCAL,
                    "exch_kcal": psi4.variable("SAPT EXCH ENERGY") * HARTREE2KCAL,
                    "ind_kcal": psi4.variable("SAPT IND ENERGY") * HARTREE2KCAL,
                    "disp_kcal": psi4.variable("SAPT DISP ENERGY") * HARTREE2KCAL,
                    "total_kcal": psi4.variable("SAPT TOTAL ENERGY") * HARTREE2KCAL,
                }
            else:  # supermolecular counterpoise interaction energy
                opts["scf_type"] = "df"
                # aug-cc-pVDZ + counterpoise ghosts → near-singular overlap
                # (recip cond ~7e-8); default symmetric orthogonalization is
                # numerically unstable → force canonical (drops the near-
                # linearly-dependent functions) for a clean SCF MO basis.
                opts["s_orthogonalization"] = "canonical"
                opts["s_tolerance"] = 1.0e-6
                if "ccsd" in method.lower():
                    # DF-CCSD(T): the conventional (IWL) 2e-integral transform
                    # overflows psi4's 32-bit file indexing at ~11e9 integrals
                    # for the aug-cc-pVDZ dimer → SIGSEGV in the SO-Ints presort
                    # (GitHub psi4 #35; dmesg-confirmed 2026-07-05, memory- AND
                    # canonical-independent — the light monomers passed, only the
                    # heavy dimer overflowed). Density-fitting the MP2 (for the
                    # FNO natural orbitals) and the CC replaces the 4-index sort
                    # with 3-index tensors, removing the overflow entirely.
                    # ~0.1 kcal DF error, cancels in the MTR−WT differential.
                    # Smoke-verified 2026-07-05 (scratchpad/df_ccsdt_smoke.py:
                    # DF+FNO+CP all PASS on a water dimer / aug-cc-pVDZ).
                    opts["mp2_type"] = "df"
                    opts["cc_type"] = "df"
                psi4.set_options(opts)
                e_int = psi4.energy(method, molecule=mol, bsse_type="cp")
                entry = {"method": method, "basis": basis_use, "mode": mode,
                         "total_kcal": float(e_int) * HARTREE2KCAL}
            results["methods"][m] = entry
        except Exception as exc:  # noqa: BLE001 -- isolate per-method failure
            results["methods"][m] = {"error": "%s: %s" % (type(exc).__name__, exc)}
        _flush()
    _flush()
    print("[qm-driver] %s -> %s" % (pair_key(job["variant"], job["tier"]), out_json))
    return out_json


# ==========================================================================
# FF driver (runs under the qmmm env python) -- OpenMM nonbonded interaction
# ==========================================================================
def _parse_ff_templates(protein_ff_xml, mtr_xml):
    """Return (charges, ljtypes, lj) from the production FF XMLs.

    charges[(resname, atomname)] = float
    ljtypes[(resname, atomname)] = amber type string (e.g. 'protein-CT')
    lj[type] = (sigma_nm, epsilon_kJ)
    """
    import xml.etree.ElementTree as ET
    charges, ljtypes, lj = {}, {}, {}
    for xmlf in (protein_ff_xml, mtr_xml):
        root = ET.parse(xmlf).getroot()
        res_root = root.find("Residues")
        if res_root is not None:
            for res in res_root.findall("Residue"):
                rn = res.get("name")
                for at in res.findall("Atom"):
                    key = (rn, at.get("name"))
                    if at.get("charge") is not None:
                        charges[key] = float(at.get("charge"))
                    if at.get("type") is not None:
                        ljtypes[key] = at.get("type")
        nb = root.find("NonbondedForce")
        if nb is not None:
            for at in nb.findall("Atom"):
                t = at.get("type")
                if t is not None and at.get("sigma") is not None:
                    lj[t] = (float(at.get("sigma")), float(at.get("epsilon")))
    return charges, ljtypes, lj


# nonpolar cap H uses the amber14 aliphatic-H LJ type; charge is neutralization sink.
_CAP_LJ_TYPE = "protein-HC"


def _fragment_ff_params(frag, charges, ljtypes, lj):
    """Per-atom (charge, sigma_nm, eps_kJ) with per-fragment neutralization of caps."""
    q, sig, eps = [], [], []
    cap_idx = []
    net_real = 0.0
    for i, a in enumerate(frag):
        if a["role"] == "cap":
            cap_idx.append(i)
            q.append(0.0)
            s, e = lj[_CAP_LJ_TYPE]
            sig.append(s)
            eps.append(e)
        else:
            key = (a["ff_source"], a["ff_name"])
            qa = charges[key]
            net_real += qa
            q.append(qa)
            t = ljtypes[key]
            s, e = lj[t]
            sig.append(s)
            eps.append(e)
    if cap_idx:                       # caps absorb the residual -> neutral fragment
        share = -net_real / len(cap_idx)
        for i in cap_idx:
            q[i] = share
    return q, sig, eps


def _nb_energy(charges, sigmas, epsilons, positions_nm):
    """OpenMM NonbondedForce (NoCutoff, no exclusions) single-point, kJ/mol."""
    import openmm as mm
    system = mm.System()
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    for q, s, e in zip(charges, sigmas, epsilons):
        system.addParticle(1.0)
        nb.addParticle(q, s, e)
    system.addForce(nb)
    integ = mm.VerletIntegrator(1.0)
    ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(positions_nm)
    e = ctx.getState(getEnergy=True).getPotentialEnergy()
    from openmm import unit
    return e.value_in_unit(unit.kilojoule_per_mole)


def run_ff_driver(job_path):
    """FF interaction energy E_int = E_nb(dimer) - E_nb(A) - E_nb(B).

    No exclusions in any of the three systems -> intra-A and intra-B cancel and
    the difference is exactly the intermolecular (cross) nonbonded energy.
    Writes <job>.ff.json.
    """
    with open(job_path) as f:
        job = json.load(f)
    import openmm.app as app
    protein_ff = os.path.join(os.path.dirname(app.__file__),
                              "data", "amber14", "protein.ff14SB.xml")
    charges, ljtypes, lj = _parse_ff_templates(protein_ff, job["mtr_xml"])

    qA, sA, eA = _fragment_ff_params(job["A"], charges, ljtypes, lj)
    qB, sB, eB = _fragment_ff_params(job["B"], charges, ljtypes, lj)
    posA = [[c / 10.0 for c in a["xyz"]] for a in job["A"]]   # A -> nm
    posB = [[c / 10.0 for c in a["xyz"]] for a in job["B"]]

    e_A = _nb_energy(qA, sA, eA, posA)
    e_B = _nb_energy(qB, sB, eB, posB)
    e_AB = _nb_energy(qA + qB, sA + sB, eA + eB, posA + posB)
    kj2kcal = 1.0 / 4.184
    entry = {
        "variant": job["variant"], "tier": job["tier"],
        "e_int_kcal": (e_AB - e_A - e_B) * kj2kcal,
        "e_dimer_kJ": e_AB, "e_A_kJ": e_A, "e_B_kJ": e_B,
        "net_charge_A": round(sum(qA), 6), "net_charge_B": round(sum(qB), 6),
    }
    out_json = os.path.splitext(job_path)[0] + ".ff.json"
    with open(out_json, "w") as f:
        json.dump(entry, f, indent=2)
    print("[ff-driver] %s -> %s" % (pair_key(job["variant"], job["tier"]), out_json))
    return out_json


# ==========================================================================
# Orchestration
# ==========================================================================
def _write_job(workdir, pair, extra=None):
    job = dict(pair)
    job["workdir"] = workdir
    job["mtr_xml"] = MTR_XML
    job["psi_scratch"] = PSI_SCRATCH_DIR
    job["memory"] = PSI_MEMORY
    job["threads"] = PSI_THREADS
    if extra:
        job.update(extra)
    path = os.path.join(workdir, pair_key(pair["variant"], pair["tier"]) + ".job.json")
    with open(path, "w") as f:
        json.dump(job, f, indent=2)
    return path


def _write_xyz(workdir, pair):
    """Human-readable dimer .xyz (fragment A then B) for inspection."""
    all_atoms = pair["A"] + pair["B"]
    path = os.path.join(workdir, pair_key(pair["variant"], pair["tier"]) + ".xyz")
    with open(path, "w") as f:
        f.write("%d\n%s A=%d B=%d\n" % (len(all_atoms), pair_key(pair["variant"], pair["tier"]),
                                        len(pair["A"]), len(pair["B"])))
        for a in all_atoms:
            f.write("%2s %14.6f %14.6f %14.6f\n" % (a["element"], *a["xyz"]))
    return path


def _run_sub(env_py, flag, job_path):
    env = dict(os.environ)
    if os.path.isdir(PSI_SCRATCH_DIR):
        env.setdefault("PSI_SCRATCH", PSI_SCRATCH_DIR)
    r = subprocess.run([env_py, _THIS, flag, job_path], check=False, env=env)
    if r.returncode != 0:
        print("[warn] sub-job rc=%d (%s) — continuing to next pair" % (
            r.returncode, os.path.basename(job_path)))
    return r.returncode


def run_benchmark(args):
    workdir = os.path.abspath(args.workdir)
    os.makedirs(workdir, exist_ok=True)
    atoms = parse_pdb(CRYSTAL_PDB)

    variants = args.variants
    tiers = args.tiers
    pairs = {}
    for v in variants:
        for t in tiers:
            pair = build_pair(v, t, atoms)
            pairs[pair_key(v, t)] = pair
            _write_xyz(workdir, pair)
    _report_fragments(pairs)

    if args.build_only:
        print("\n[build-only] fragments + xyz written to %s" % workdir)
        return

    # QM + FF sub-jobs.
    qm_results, ff_results = {}, {}
    for key, pair in pairs.items():
        job_path = _write_job(workdir, pair,
                              extra={"methods": args.qm_methods} if args.qm_methods else None)
        if not args.ff_only:
            _run_sub(PSI4_PY, "--_qm-driver", job_path)
            qf = os.path.splitext(job_path)[0] + ".qm.json"
            if os.path.exists(qf):
                qm_results[key] = json.load(open(qf))
        if not args.qm_only:
            _run_sub(QMMM_PY, "--_ff-driver", job_path)
            ff = os.path.splitext(job_path)[0] + ".ff.json"
            if os.path.exists(ff):
                ff_results[key] = json.load(open(ff))

    summary = aggregate(pairs, qm_results, ff_results, tiers)
    out = os.path.join(workdir, "sapt_dimer_benchmark_summary.json")
    with open(out, "w") as f:
        json.dump(summary, f, indent=2)
    _print_summary(summary)
    print("\n[done] summary -> %s" % out)


# ==========================================================================
# Aggregation + pre-registered outcome classification
# ==========================================================================
CHEM_ACCURACY = 1.0     # kcal/mol (verdict Q4 threshold "~1 kcal/mol")


def _e_int(res, key, method):
    r = res.get(key, {}).get("methods", {}).get(method)
    if r and "total_kcal" in r:
        return r["total_kcal"]
    return None


def aggregate(pairs, qm_results, ff_results, tiers):
    """ddE_int(variant) = E_int(variant) - E_int(WT), per method + FF, per tier."""
    summary = {
        "geometry_ssot": os.path.relpath(CRYSTAL_PDB, REPO),
        "prereg_outcomes": {
            "i_ff_fixable": "|ddE_int^CCSD(T)| >= ~1 kcal AND SAPT-attributed AND FF deviates in that component",
            "ii_ranking_only": "|ddE_int^CCSD(T)| <~ chem accuracy AND within method spread",
            "iii_ff_structural_error": "sign(ddE^FF) != sign(ddE^CCSD(T)) OR >2-3x AND SAPT names an omitted component",
            "w4a_control": "|ddE_int^CCSD(T)| MUST be large & unfavorable; else control FAIL -> MTR untrustworthy",
        },
        "caveats_R11_R18": "ddE_int = gas-phase fixed-geometry 2-body interaction energy, NOT ddG_bind; "
                           "Magotti -1.4/-3.0 kcal is a different quantity (SIGN suggestive only).",
        "tiers": {},
    }
    methods = list(QM_METHODS)
    for t in tiers:
        wt = pair_key("WT", t)
        tier_block = {"per_variant": {}, "ddE_int": {}, "sapt_decomposition_ddE": {}}
        # absolute E_int per variant
        for v in ("WT", "MTR", "W4A"):
            k = pair_key(v, t)
            if k not in pairs:
                continue
            row = {"ff_kcal": ff_results.get(k, {}).get("e_int_kcal")}
            for m in methods:
                row[m + "_kcal"] = _e_int(qm_results, k, m)
            tier_block["per_variant"][v] = row
        # ddE relative to WT
        wt_ff = ff_results.get(wt, {}).get("e_int_kcal")
        for v in ("MTR", "W4A"):
            k = pair_key(v, t)
            if k not in pairs:
                continue
            dd = {}
            v_ff = ff_results.get(k, {}).get("e_int_kcal")
            dd["ff_kcal"] = None if (v_ff is None or wt_ff is None) else v_ff - wt_ff
            for m in methods:
                ev, ew = _e_int(qm_results, k, m), _e_int(qm_results, wt, m)
                dd[m + "_kcal"] = None if (ev is None or ew is None) else ev - ew
            tier_block["ddE_int"][v] = dd
            # SAPT component ddE (from SAPT0 as primary decomposition)
            tier_block["sapt_decomposition_ddE"][v] = _sapt_dd(qm_results, k, wt)
        tier_block["classification"] = {
            v: classify_outcome(tier_block, v) for v in tier_block["ddE_int"]
        }
        summary["tiers"]["tier%d" % t] = tier_block
    return summary


def _sapt_dd(qm_results, key, wt_key):
    comps = ("elst", "exch", "ind", "disp")
    out = {}
    for level in ("sapt0", "sapt2p3_dmp2"):
        v = qm_results.get(key, {}).get("methods", {}).get(level, {})
        w = qm_results.get(wt_key, {}).get("methods", {}).get(level, {})
        if "elst_kcal" in v and "elst_kcal" in w:
            out[level] = {c + "_kcal": v[c + "_kcal"] - w[c + "_kcal"] for c in comps}
    return out


def classify_outcome(tier_block, variant):
    """Assign the pre-registered outcome (i/ii/iii) + W4A control, or 'incomplete'."""
    dd = tier_block["ddE_int"].get(variant, {})
    gold = dd.get("fno_ccsd_t_kcal")
    silver = dd.get("sapt2p3_dmp2_kcal")
    ff = dd.get("ff_kcal")
    ref = gold if gold is not None else silver
    if ref is None or ff is None:
        return {"outcome": "incomplete", "reason": "missing CCSD(T)/SAPT2+(3) or FF ddE"}
    note = {"ddE_ref_kcal": ref, "ddE_ff_kcal": ff, "ref_method": "fno_ccsd_t" if gold is not None else "sapt2p3_dmp2"}
    if variant == "W4A":
        note["control_pass"] = bool(abs(ref) >= 2.0 * CHEM_ACCURACY)
        note["outcome"] = "control_ok" if note["control_pass"] else "control_fail"
        return note
    same_sign = (ff == 0 and ref == 0) or (ff * ref > 0)
    if (not same_sign) or (abs(ref) > 1e-9 and abs(ff - ref) > 2.0 * abs(ref)):
        note["outcome"] = "iii_ff_structural_error"
    elif abs(ref) < CHEM_ACCURACY:
        note["outcome"] = "ii_ranking_only"
    else:
        note["outcome"] = "i_ff_fixable"
    return note


# ==========================================================================
# Reporting
# ==========================================================================
def _report_fragments(pairs):
    print("=" * 74)
    print(" Fragment inventory (crystal 2QKI, vertical substitution, C1/C2)")
    print("=" * 74)
    for key in sorted(pairs):
        p = pairs[key]
        na, nb = len(p["A"]), len(p["B"])
        print("  %-10s  monomer A = %2d atoms   monomer B = %2d atoms   (dimer %d)"
              % (key, na, nb, na + nb))


def _fmt(x):
    return "   n/a  " if x is None else "%8.3f" % x


def _print_summary(summary):
    print("\n" + "=" * 74)
    print(" ddE_int triangulation  (variant - WT, kcal/mol)  R-11: SIGN/decomp only")
    print("=" * 74)
    for tname, tb in summary["tiers"].items():
        print("\n-- %s --" % tname)
        print("  %-5s %9s %9s %9s %9s %9s"
              % ("var", "FF", "wB97X-D", "SAPT0", "SAPT2+(3)", "CCSD(T)"))
        for v, dd in tb["ddE_int"].items():
            print("  %-5s %s %s %s %s %s" % (
                v, _fmt(dd.get("ff_kcal")), _fmt(dd.get("wb97xd_kcal")),
                _fmt(dd.get("sapt0_kcal")), _fmt(dd.get("sapt2p3_dmp2_kcal")),
                _fmt(dd.get("fno_ccsd_t_kcal"))))
        for v, cl in tb.get("classification", {}).items():
            print("     %-5s -> %s" % (v, cl.get("outcome", "?")))
    print("\n  Pre-registered outcomes (anti-HARKing):")
    for k, txt in summary["prereg_outcomes"].items():
        print("    - %s: %s" % (k, txt))
    print("  %s" % summary["caveats_R11_R18"])


# ==========================================================================
# CLI
# ==========================================================================
def build_parser():
    p = argparse.ArgumentParser(
        description="SAPT0/FNO-CCSD(T)/wB97X-D/FF reference-dimer benchmark (2QKI Trp4).")
    p.add_argument("--workdir", default=os.path.join(REPO, "analysis", "sapt_dimer_benchmark"),
                   help="output directory for jobs/xyz/results")
    p.add_argument("--variants", nargs="+", default=["WT", "MTR", "W4A"],
                   choices=["WT", "MTR", "W4A"])
    p.add_argument("--tiers", nargs="+", type=int, default=[1, 2], choices=[1, 2])
    p.add_argument("--qm-methods", nargs="+", default=None,
                   help="subset of QM method keys (default: all applicable per tier)")
    p.add_argument("--build-only", action="store_true",
                   help="build fragments + xyz only, no QM/FF")
    p.add_argument("--qm-only", action="store_true", help="run QM stack only")
    p.add_argument("--ff-only", action="store_true", help="run FF single-point only")
    p.add_argument("--dry-run", action="store_true",
                   help="build fragments + FF + one fast SAPT0 (sto-3g) smoke on WT tier-1")
    # hidden driver subcommands (self re-invocation)
    p.add_argument("--_qm-driver", dest="qm_driver", default=None, help=argparse.SUPPRESS)
    p.add_argument("--_ff-driver", dest="ff_driver", default=None, help=argparse.SUPPRESS)
    return p


def run_dry_run(args):
    """implementation self-check: fragments + FF + one fast SAPT0 decomposition on WT tier-1."""
    workdir = os.path.abspath(args.workdir)
    os.makedirs(workdir, exist_ok=True)
    atoms = parse_pdb(CRYSTAL_PDB)
    pair = build_pair("WT", 1, atoms)
    _report_fragments({pair_key("WT", 1): pair})
    _write_xyz(workdir, pair)
    job_path = _write_job(workdir, pair,
                          extra={"methods": ["sapt0"],
                                 "basis_override": {"sapt0": "sto-3g"}})
    print("\n[dry-run] FF single-point (OpenMM, qmmm env) ...")
    _run_sub(QMMM_PY, "--_ff-driver", job_path)
    ff = json.load(open(os.path.splitext(job_path)[0] + ".ff.json"))
    print("  FF E_int(WT tier-1) = %.3f kcal/mol  (netA=%s netB=%s)"
          % (ff["e_int_kcal"], ff["net_charge_A"], ff["net_charge_B"]))
    print("\n[dry-run] SAPT0/sto-3g smoke (psi4 env; production uses jun-cc-pVDZ) ...")
    _run_sub(PSI4_PY, "--_qm-driver", job_path)
    qm = json.load(open(os.path.splitext(job_path)[0] + ".qm.json"))
    s0 = qm["methods"].get("sapt0", {})
    if "error" in s0:
        print("  SAPT0 smoke ERROR: %s" % s0["error"])
    else:
        print("  SAPT0 decomposition (kcal/mol): elst=%.3f exch=%.3f ind=%.3f disp=%.3f total=%.3f"
              % (s0["elst_kcal"], s0["exch_kcal"], s0["ind_kcal"], s0["disp_kcal"], s0["total_kcal"]))
    print("\n[dry-run] OK -- plumbing verified. Full stack (jun-cc-pVDZ + FNO-CCSD(T)/aVTZ) is execution review's job.")


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.qm_driver:
        run_qm_driver(args.qm_driver)
        return
    if args.ff_driver:
        run_ff_driver(args.ff_driver)
        return
    if args.dry_run:
        run_dry_run(args)
        return
    run_benchmark(args)


if __name__ == "__main__":
    main()
