#!/usr/bin/env python
"""ADR-0020 §A probe gate: densified apex가 marginal cohort zero-RT를 복구했나.
3 cell: s101_wt(fixed), s7_wt(fixed), s199_cp4(auto-search). densified=13 states/dir, apex=state6.
원본(11 state/apex5) marginal zero-RT legs: s101 dplus+dminus, s7 dplus, s199 dminus.
metric: per-walker full round-trip(col0=state, 0<->max traversal). [[reference_trackb_out_index_is_replica]]"""
import os, glob, numpy as np
G="/home/san/UPDD_proj/analysis/apex_probe_gate_20260628"
NS=13; APEX=6  # densified: 2*6+1=13 states, apex middle

def walkers(tag, d):
    """rN/trackb_{d}.out = walker N trace; col0=state."""
    W={}
    for w in range(NS):
        f=f"{G}/{tag}/{d}/r{w}/trackb_{d}.out"
        if not os.path.exists(f): continue
        a=np.loadtxt(f)
        if a.ndim==1: a=a[None,:]
        W[w]=a[:,0].astype(int)
    return W

def roundtrips(W):
    """count 0<->max round trips summed over walkers + spans."""
    tot=0; spans=[]; per=[]
    for w,st in sorted(W.items()):
        spans.append((int(st.min()),int(st.max())))
        hits=[]
        for s in st:
            if s==0 or s==NS-1:
                if not hits or hits[-1]!=s: hits.append(s)
        alt=sum(1 for i in range(1,len(hits)) if hits[i]!=hits[i-1])
        rt=alt//2; tot+=rt; per.append(rt)
    return tot, spans, per

def apex_cross(W):
    """boundary crossings at apex: (APEX-1<->APEX) and (APEX<->APEX+1)."""
    c_lo=0; c_hi=0
    for st in W.values():
        for i in range(1,len(st)):
            a,b=st[i-1],st[i]
            lo,hi=min(a,b),max(a,b)
            if lo==APEX-1 and hi==APEX: c_lo+=1
            if lo==APEX and hi==APEX+1: c_hi+=1
    return c_lo,c_hi

CELLS=[("s101_wt","wt/s101","both"),("s7_wt","wt/s7","dplus"),("s199_cp4","cp4/s199","dminus")]
print("="*78)
print("ADR-0020 §A PROBE GATE — densified apex round-trip 복구 평가")
print("원본 marginal zero-RT: s101(both), s7/dplus, s199/dminus  →  densify 후 RT≥1 ?")
print("="*78)
verdict={}
for tag,label,marg in CELLS:
    print(f"\n### {label}  (marginal leg: {marg})")
    cell_ok=True
    for d in ("dplus","dminus"):
        W=walkers(tag,d)
        if not W: print(f"  {d}: NO DATA"); cell_ok=False; continue
        rt,spans,per=roundtrips(W)
        nfull=sum(1 for lo,hi in spans if lo==0 and hi==NS-1)
        clo,chi=apex_cross(W)
        ncyc=max(len(s) for s in W.values())
        flag = "RT>=1 ✅" if rt>=1 else "RT=0 ❌"
        marg_hit = (marg=="both") or (marg==d)
        tag_m = "  <<MARGINAL" if marg_hit else ""
        print(f"  {d}: total_RT={rt} ({flag})  full-span-walkers={nfull}/{len(W)}  apex_cross[{APEX-1}-{APEX}]={clo} [{APEX}-{APEX+1}]={chi}  ncyc={ncyc}{tag_m}")
        if marg_hit and rt<1: cell_ok=False
    verdict[label]=cell_ok

print("\n"+"="*78)
print("VERDICT (marginal legs RT>=1 복구 여부):")
allok=True
for tag,label,marg in CELLS:
    ok=verdict[label]; allok = allok and ok
    print(f"  {label} (marg {marg}): {'✅ 복구' if ok else '❌ 미복구'}")
print(f"\n  1차 게이트(전 marginal RT>=1): {'🟢 PASS' if allok else '🔴 FAIL'}")
print("  ※ apex MBAR overlap O>=0.10 + frozen 서명 해소는 round-trip PASS 시 별도 확증 단계")
print("="*78)
