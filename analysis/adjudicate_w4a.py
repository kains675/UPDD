#!/usr/bin/env python3
"""Adjudicate W4A bound non-mixing: recompute per-bond p_acc two ways."""
import sys, os, math, glob
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "utils"))
from atm_trackB_inplace_rbfe import _atm_softcore_components_kj, _atm_state_energy_kj
import atm_trackB_inplace_rbfe as M

KCAL2KJ = M.RBFE_TEMP_K  # placeholder, will set properly
def kcal2kj(x):
    import openmm.unit as u
    return (x * u.kilocalorie_per_mole).value_in_unit(u.kilojoule_per_mole)

# kT in kcal at 300K
kT_kcal = 0.0019872041 * 300.0
# umax/ubcore/acore canon (from schedule). Read from build to be sure.
# Defaults used in two-copy ladder:
ACORE = M.RBFE_ALPHA_DEFAULT  # alpha is per-state; acore is softcore param, separate
# Need umax,ubcore,acore from the schedule canon
sch = M.build_ats_standard_ladder(n_windows_half=6, single_direction="backward")
UMAX = sch["umax"]; UBCORE = sch["ubcore"]; ACORE = sch["acore"]
print("CANON umax=%.3f ubcore=%.3f acore=%.4f (kcal)" % (UMAX, UBCORE, ACORE), file=sys.stderr)
umax_kj = kcal2kj(UMAX); ubcore_kj = kcal2kj(UBCORE)

def load_dir(d):
    """Return rows[r] = list of (col0_state, lam1, lam2, alpha, u0, w0, potE_kcal, usc_kcal, direction)."""
    rows = {}
    for f in sorted(glob.glob(os.path.join(d, "r*/trackb_*.out"))):
        r = int(os.path.basename(os.path.dirname(f))[1:])
        rr = []
        with open(f) as fh:
            for line in fh:
                p = line.split()
                if len(p) < 11: continue
                rr.append((int(p[0]), float(p[3]), float(p[4]), float(p[5]),
                           float(p[6]), float(p[7]), float(p[8]), float(p[9]), int(p[2])))
        rows[r] = rr
    return rows

def reconstruct_base_kj(potE_kcal, usc_kcal, lam1, lam2, alpha_kcal, uh_kcal, w0_kcal, direction):
    """Invert the write formula: potE = base + hybrid(usc). Return base in kJ."""
    usc_kj = kcal2kj(usc_kcal)
    uh_kj = kcal2kj(uh_kcal)
    w0_kj = kcal2kj(w0_kcal)
    # alpha stored is per-state 'alpha' col; engine uses alpha_per_kj. The ladder's
    # alpha is in 1/kcal? In _atm_state_energy_kj, alpha_per_kj multiplies usc (kj).
    # The .out col5 'alpha' = schedule alpha (0.1). Convert: per-kcal -> per-kj.
    # hybrid uses alpha_per_kj * (usc-uh). We need consistent units. The engine
    # stores alpha in per-kcal; convert to per-kj by dividing by kcal2kj(1).
    alpha_per_kj = alpha_kcal / kcal2kj(1.0)
    if lam2 != lam1 and alpha_per_kj > 0.0:
        hybrid = (((lam2 - lam1) / alpha_per_kj) * math.log(1.0 + math.exp(-alpha_per_kj*(usc_kj - uh_kj)))
                  + lam2*usc_kj + w0_kj)
    else:
        hybrid = lam2*usc_kj + w0_kj
    potE_kj = kcal2kj(potE_kcal)
    return potE_kj - hybrid

def state_energy_from_base(base_kj, usc_kj, lam1, lam2, alpha_kcal, uh_kcal, w0_kcal):
    uh_kj = kcal2kj(uh_kcal); w0_kj = kcal2kj(w0_kcal)
    alpha_per_kj = alpha_kcal / kcal2kj(1.0)
    if lam2 != lam1 and alpha_per_kj > 0.0:
        hybrid = (((lam2 - lam1)/alpha_per_kj)*math.log(1.0+math.exp(-alpha_per_kj*(usc_kj-uh_kj)))
                  + lam2*usc_kj + w0_kj)
    else:
        hybrid = lam2*usc_kj + w0_kj
    return base_kj + hybrid

def per_state_schedule(rows):
    """state -> (lam1,lam2,alpha,u0,w0,direction) from any row at that state."""
    sched = {}
    for r, rr in rows.items():
        for (st, l1, l2, al, u0, w0, pe, usc, dr) in rr:
            if st not in sched:
                sched[st] = (l1, l2, al, u0, w0, dr)
    return sched

def correct_pacc(rows):
    """State-grouped: per cycle k, gather (state, base_kj, usc_kj) from row k of each r-file.
    Then for adjacent states (b,b+1), compute swap delta and Metropolis p."""
    sched = per_state_schedule(rows)
    ncyc = min(len(rr) for rr in rows.values())
    # config[c][state] = (base_kj, usc_kj)
    accept = {}; trials = {}
    kT_kj = kcal2kj(kT_kcal)
    for c in range(ncyc):
        by_state = {}
        for r, rr in rows.items():
            st, l1, l2, al, u0, w0, pe, usc, dr = rr[c]
            base_kj = reconstruct_base_kj(pe, usc, l1, l2, al, u0, w0, dr)
            usc_kj = kcal2kj(usc)
            by_state[st] = (base_kj, usc_kj)
        for b in range(10):
            if b not in by_state or (b+1) not in by_state: continue
            base_lo, usc_lo = by_state[b]
            base_hi, usc_hi = by_state[b+1]
            l1_lo,l2_lo,al_lo,u0_lo,w0_lo,_ = sched[b]
            l1_hi,l2_hi,al_hi,u0_hi,w0_hi,_ = sched[b+1]
            # energy of config-lo at state lo / hi ; config-hi at hi/lo
            e_lo_lo = state_energy_from_base(base_lo, usc_lo, l1_lo,l2_lo,al_lo,u0_lo,w0_lo)
            e_lo_hi = state_energy_from_base(base_lo, usc_lo, l1_hi,l2_hi,al_hi,u0_hi,w0_hi)
            e_hi_hi = state_energy_from_base(base_hi, usc_hi, l1_hi,l2_hi,al_hi,u0_hi,w0_hi)
            e_hi_lo = state_energy_from_base(base_hi, usc_hi, l1_lo,l2_lo,al_lo,u0_lo,w0_lo)
            delta = (e_lo_hi + e_hi_lo) - (e_lo_lo + e_hi_hi)
            p = 1.0 if delta <= 0 else math.exp(-delta/kT_kj)
            accept[b] = accept.get(b,0.0) + p
            trials[b] = trials.get(b,0) + 1
    return {b: accept[b]/trials[b] for b in trials}

def dir_grouped_pacc(rows):
    """WRONG method: treat each r-file as a fixed ladder index, pair adjacent r-files."""
    sched = per_state_schedule(rows)
    ncyc = min(len(rr) for rr in rows.values())
    accept = {}; trials = {}
    kT_kj = kcal2kj(kT_kcal)
    for c in range(ncyc):
        for b in range(10):
            if b not in rows or (b+1) not in rows: continue
            st_lo, l1l,l2l,all_,u0l,w0l,pel,uscl,drl = rows[b][c]
            st_hi, l1h,l2h,alh,u0h,w0h,peh,usch,drh = rows[b+1][c]
            base_lo = reconstruct_base_kj(pel,uscl,l1l,l2l,all_,u0l,w0l,drl)
            base_hi = reconstruct_base_kj(peh,usch,l1h,l2h,alh,u0h,w0h,drh)
            usc_lo_kj = kcal2kj(uscl); usc_hi_kj = kcal2kj(usch)
            # use the DIRECTORY's own lambdas as if r is a fixed state (the artifact)
            e_lo_lo = state_energy_from_base(base_lo, usc_lo_kj, l1l,l2l,all_,u0l,w0l)
            e_lo_hi = state_energy_from_base(base_lo, usc_lo_kj, l1h,l2h,alh,u0h,w0h)
            e_hi_hi = state_energy_from_base(base_hi, usc_hi_kj, l1h,l2h,alh,u0h,w0h)
            e_hi_lo = state_energy_from_base(base_hi, usc_hi_kj, l1l,l2l,all_,u0l,w0l)
            delta = (e_lo_hi + e_hi_lo) - (e_lo_lo + e_hi_hi)
            p = 1.0 if delta <= 0 else math.exp(-delta/kT_kj)
            accept[b] = accept.get(b,0.0)+p
            trials[b] = trials.get(b,0)+1
    return {b: accept[b]/trials[b] for b in trials}

def usc_means(rows):
    sched = per_state_schedule(rows)
    by = {}
    for r, rr in rows.items():
        for (st,l1,l2,al,u0,w0,pe,usc,dr) in rr:
            by.setdefault(st,[]).append(usc)
    return {s: sum(v)/len(v) for s,v in by.items()}

if __name__ == "__main__":
    d = sys.argv[1]
    rows = load_dir(d)
    print("DIR:", d)
    cp = correct_pacc(rows)
    dp = dir_grouped_pacc(rows)
    um = usc_means(rows)
    print("\nstate usc-means:", " ".join("%.1f"%um[s] for s in sorted(um)))
    gaps = [um[s+1]-um[s] for s in range(10)]
    print("consecutive usc gaps:", " ".join("%.1f"%g for g in gaps))
    print("\nbond  p_acc(STATE-correct)  p_acc(DIR-artifact)  usc_gap")
    for b in range(10):
        print("%d/%d   %.4f               %.4f             %.1f" % (b,b+1, cp.get(b,float('nan')), dp.get(b,float('nan')), gaps[b]))
