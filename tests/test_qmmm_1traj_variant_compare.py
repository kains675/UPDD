"""tests/test_qmmm_1traj_variant_compare.py — Track A 1-traj swap primitive.

Covers utils/qmmm_1traj_variant_compare.py (schema ``qm1traj_variant_pair/0.1``).

Self-test scope (NO production QM — those are Keeper-gated):
    T1: swap_residue graft (Trp→1-Me-Trp): drops HE1, adds CM+HM1/2/3, renames MTR.
    T2: swap_residue strip (1-Me-Trp→Trp): drops CM+HM*, restores HE1, renames TRP.
    T3: graft→strip round-trip restores the original heavy-atom set/coords exactly.
    T4: swap direction guard — graft on an MTR frame (or strip on TRP) fails fast.
    T5: clash_check flags a known graft clash (methyl into a planted neighbor).
    T6: clash_check is clean when no environment atom is near the added methyl.
    T7: CM bond length to NE1 ≈ 1.466 Å; restored HE1 ≈ 1.01 Å (geometry sanity).

Isolation: load the module under test by filespec (importlib) so the test does
not depend on pytest-session sys.path ordering — mirrors tests/test_updd_cli.py.
"""

from __future__ import annotations

import importlib.util as _ilu
import math
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_UTILS_DIR = _REPO_ROOT / "utils"
if str(_UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(_UTILS_DIR))

_MOD_PATH = _UTILS_DIR / "qmmm_1traj_variant_compare.py"
_spec = _ilu.spec_from_file_location("qmmm_1traj_variant_compare_uut", _MOD_PATH)
assert _spec is not None and _spec.loader is not None
_mod = _ilu.module_from_spec(_spec)
sys.modules["qmmm_1traj_variant_compare_uut"] = _mod
_spec.loader.exec_module(_mod)  # type: ignore[union-attr]
m = _mod


# ---------------------------------------------------------------------------
# Fixtures — minimal in-memory PDB residue-4 (Trp / 1-Me-Trp) on chain B
# ---------------------------------------------------------------------------
def _atom(serial, name, resname, chain, resnum, x, y, z, elem):
    buf = list(" " * 80)
    buf[0:6] = list("ATOM  ")
    buf[6:11] = list(f"{serial:5d}")
    # 1-letter element names are right-justified into cols 13-14 (PDB convention).
    nm = name if len(name) >= 4 else f" {name:<3s}"
    buf[12:16] = list(nm[:4])
    buf[17:20] = list(f"{resname[:3]:<3s}")
    buf[21] = chain
    buf[22:26] = list(f"{resnum:4d}")
    buf[30:38] = list(f"{x:8.3f}")
    buf[38:46] = list(f"{y:8.3f}")
    buf[46:54] = list(f"{z:8.3f}")
    buf[54:60] = list("  1.00")
    buf[60:66] = list("  0.00")
    buf[76:78] = list(elem.rjust(2))
    return "".join(buf).rstrip() + "\n"


# A planar-ish indole skeleton sufficient for the swap geometry: we only need
# N/CA/C/O backbone + CB/CG/CD1/CD2/NE1/CE2/CE3/CZ2/CZ3/CH2 heavy atoms plus the
# divergent atoms. Coordinates are a compact synthetic indole (not a real MD
# frame) — the swap is geometric, so exact biophysics is irrelevant here.
def _trp_residue4():
    a = []
    s = 1
    # backbone
    a.append(_atom(s, "N", "TRP", "B", 4, 0.000, 0.000, 0.000, "N")); s += 1
    a.append(_atom(s, "H", "TRP", "B", 4, -0.500, 0.800, 0.000, "H")); s += 1
    a.append(_atom(s, "CA", "TRP", "B", 4, 1.450, 0.000, 0.000, "C")); s += 1
    a.append(_atom(s, "HA", "TRP", "B", 4, 1.800, 1.020, 0.000, "H")); s += 1
    a.append(_atom(s, "C", "TRP", "B", 4, 2.000, -0.700, 1.250, "C")); s += 1
    a.append(_atom(s, "O", "TRP", "B", 4, 3.200, -0.700, 1.500, "O")); s += 1
    # sidechain
    a.append(_atom(s, "CB", "TRP", "B", 4, 2.000, -0.750, -1.250, "C")); s += 1
    a.append(_atom(s, "HB2", "TRP", "B", 4, 1.600, -1.760, -1.300, "H")); s += 1
    a.append(_atom(s, "HB3", "TRP", "B", 4, 3.080, -0.800, -1.200, "H")); s += 1
    a.append(_atom(s, "CG", "TRP", "B", 4, 1.600, 0.000, -2.500, "C")); s += 1
    a.append(_atom(s, "CD1", "TRP", "B", 4, 2.300, 0.100, -3.700, "C")); s += 1
    a.append(_atom(s, "HD1", "TRP", "B", 4, 3.260, -0.380, -3.820, "H")); s += 1
    a.append(_atom(s, "CD2", "TRP", "B", 4, 0.450, 0.850, -2.700, "C")); s += 1
    a.append(_atom(s, "NE1", "TRP", "B", 4, 1.650, 0.950, -4.620, "N")); s += 1
    a.append(_atom(s, "HE1", "TRP", "B", 4, 1.980, 1.180, -5.540, "H")); s += 1
    a.append(_atom(s, "CE2", "TRP", "B", 4, 0.500, 1.430, -4.020, "C")); s += 1
    a.append(_atom(s, "CE3", "TRP", "B", 4, -0.650, 1.150, -1.980, "C")); s += 1
    a.append(_atom(s, "HE3", "TRP", "B", 4, -0.700, 0.780, -0.960, "H")); s += 1
    a.append(_atom(s, "CZ2", "TRP", "B", 4, -0.450, 2.300, -4.560, "C")); s += 1
    a.append(_atom(s, "HZ2", "TRP", "B", 4, -0.400, 2.680, -5.580, "H")); s += 1
    a.append(_atom(s, "CZ3", "TRP", "B", 4, -1.580, 2.020, -2.540, "C")); s += 1
    a.append(_atom(s, "HZ3", "TRP", "B", 4, -2.470, 2.230, -1.950, "H")); s += 1
    a.append(_atom(s, "CH2", "TRP", "B", 4, -1.410, 2.620, -3.800, "C")); s += 1
    a.append(_atom(s, "HH2", "TRP", "B", 4, -2.160, 3.290, -4.220, "H")); s += 1
    a.append("TER\n")
    a.append("END\n")
    return a


def _names_chain_res(lines, chain, resnum):
    out = []
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")) and ln[21] == chain and ln[22:26].strip() == str(resnum):
            out.append(ln[12:16].strip())
    return out


def _coord_of(lines, chain, resnum, name):
    for ln in lines:
        if (ln.startswith(("ATOM", "HETATM")) and ln[21] == chain
                and ln[22:26].strip() == str(resnum) and ln[12:16].strip() == name):
            return (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return None


# ---------------------------------------------------------------------------
# T1 — graft adds the methyl + drops HE1 + renames MTR
# ---------------------------------------------------------------------------
def test_t1_graft_atom_set_and_rename():
    g, diag = m.swap_residue(_trp_residue4(), "B", 4, "graft")
    names = set(_names_chain_res(g, "B", 4))
    assert "HE1" not in names                       # indole donor removed
    assert {"CM", "HM1", "HM2", "HM3"} <= names      # methyl added
    assert diag["to_resname"] == "MTR"
    assert sorted(diag["added_atoms"]) == ["CM", "HM1", "HM2", "HM3"]
    assert diag["removed_atoms"] == ["HE1"]
    # every residue-4 line must read resname MTR now
    for ln in g:
        if ln.startswith(("ATOM", "HETATM")) and ln[21] == "B" and ln[22:26].strip() == "4":
            assert ln[17:20].strip() == "MTR"


# ---------------------------------------------------------------------------
# T2 — strip drops the methyl + restores HE1 + renames TRP
# ---------------------------------------------------------------------------
def test_t2_strip_restores_indole_nh():
    # build a 1-Me-Trp frame first by grafting, then strip it back.
    g, _ = m.swap_residue(_trp_residue4(), "B", 4, "graft")
    s, diag = m.swap_residue(g, "B", 4, "strip")
    names = set(_names_chain_res(s, "B", 4))
    assert "HE1" in names                            # buried donor restored (the WHY)
    assert not ({"CM", "HM1", "HM2", "HM3"} & names)  # methyl gone
    assert diag["to_resname"] == "TRP"
    assert diag["added_atoms"] == ["HE1"]
    assert sorted(diag["removed_atoms"]) == ["CM", "HM1", "HM2", "HM3"]


# ---------------------------------------------------------------------------
# T3 — round-trip heavy-atom identity (graft then strip == original heavies)
# ---------------------------------------------------------------------------
def test_t3_round_trip_heavy_identity():
    orig = _trp_residue4()
    sig0 = m.heavy_atom_signature(orig, "B", 4)
    g, _ = m.swap_residue(orig, "B", 4, "graft")
    s, _ = m.swap_residue(g, "B", 4, "strip")
    sig1 = m.heavy_atom_signature(s, "B", 4)
    assert sig0 == sig1, "graft→strip must restore the original heavy-atom set/coords"
    # the shared scaffold coords must be byte-identical between WT and grafted
    sigg = m.heavy_atom_signature(g, "B", 4)
    shared = {n for n, _ in sig0}
    sig0_d = dict(sig0)
    sigg_d = dict(sigg)
    for n in shared:
        assert sig0_d[n] == sigg_d[n], f"shared atom {n} moved during graft (must be frozen)"


# ---------------------------------------------------------------------------
# T4 — direction guard rejects wrong source identity (A-C1 safety)
# ---------------------------------------------------------------------------
def test_t4_direction_guard_fail_fast():
    g, _ = m.swap_residue(_trp_residue4(), "B", 4, "graft")  # now MTR
    # grafting again (MTR source) must fail — graft expects TRP
    with pytest.raises(ValueError):
        m.swap_residue(g, "B", 4, "graft")
    # stripping the original TRP frame must fail — strip expects MTR
    with pytest.raises(ValueError):
        m.swap_residue(_trp_residue4(), "B", 4, "strip")
    with pytest.raises(ValueError):
        m.swap_residue(_trp_residue4(), "B", 4, "sideways")  # bad direction


# ---------------------------------------------------------------------------
# T5 — clash_check flags a planted clash near the grafted methyl
# ---------------------------------------------------------------------------
def test_t5_clash_check_flags_known_clash():
    g, diag = m.swap_residue(_trp_residue4(), "B", 4, "graft")
    cm = _coord_of(g, "B", 4, "CM")
    assert cm is not None
    # plant a heavy environment atom 1.0 Å from CM on a different residue/chain
    intruder = _atom(999, "OX", "HOH", "C", 200, cm[0] + 1.0, cm[1], cm[2], "O")
    clashed = g + [intruder]
    res = m.clash_check(clashed, 4, chain="B")
    assert res["clash_flag"] is True
    assert any(v["res_atom"] == "CM" and v["kind"] == "heavy_heavy" for v in res["violations"])


# ---------------------------------------------------------------------------
# T6 — clash_check clean when nothing is near the added atoms
# ---------------------------------------------------------------------------
def test_t6_clash_check_clean_case():
    g, _ = m.swap_residue(_trp_residue4(), "B", 4, "graft")
    # a far-away environment atom (20 Å) must not trip the filter
    far = _atom(999, "OX", "HOH", "C", 200, 100.0, 100.0, 100.0, "O")
    res = m.clash_check(g + [far], 4, chain="B")
    assert res["clash_flag"] is False
    assert res["violations"] == []
    # the peptide N/C backbone bonds (if any synthetic neighbor existed) must NOT
    # be probed — probe set is the swap-added sidechain atoms only.
    assert set(res["probe_atoms"]) == {"CM", "HM1", "HM2", "HM3", "HE1"}


# ---------------------------------------------------------------------------
# T7 — geometry sanity: CM–NE1 ≈ 1.466 Å, restored HE1–NE1 ≈ 1.01 Å
# ---------------------------------------------------------------------------
def test_t7_bond_length_sanity():
    g, _ = m.swap_residue(_trp_residue4(), "B", 4, "graft")
    ne1 = _coord_of(g, "B", 4, "NE1")
    cm = _coord_of(g, "B", 4, "CM")
    d_cm = math.dist(ne1, cm)
    assert abs(d_cm - m._NE1_CM_BOND) < 1e-2, f"CM–NE1 = {d_cm:.3f}, expected ~{m._NE1_CM_BOND}"
    # methyl H bond lengths
    for hi in (1, 2, 3):
        h = _coord_of(g, "B", 4, f"HM{hi}")
        assert abs(math.dist(cm, h) - m._CM_HM_BOND) < 1e-2
    # restored HE1 on the strip path
    s, _ = m.swap_residue(g, "B", 4, "strip")
    ne1_s = _coord_of(s, "B", 4, "NE1")
    he1 = _coord_of(s, "B", 4, "HE1")
    d_he1 = math.dist(ne1_s, he1)
    assert abs(d_he1 - m._NE1_HE1_BOND) < 1e-2, f"HE1–NE1 = {d_he1:.3f}, expected ~{m._NE1_HE1_BOND}"
