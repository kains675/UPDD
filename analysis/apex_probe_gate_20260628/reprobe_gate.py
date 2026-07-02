#!/usr/bin/env python
"""② re-probe gate (nonparametric, cycle-독립 1차 + crossing 2차).
워크플로 정정: gaussian-Bhatt는 leg-up left-skew에 속아 overlap 과대평가 → hist-OVL(비모수)+walker-crossing으로 판정.
NS는 데이터에서 자동(H16=16). dminus fwd=NS-1-i. col0=state, col9=pertE.

사전등록 PASS (anti-HARK, re-probe 실행 前 고정):
 1차(필수, cycle-독립): 전 인접쌍 hist-OVL >= 0.45 AND gaussian-Bhatt >= 0.80 (dual; skew false-pass 차단).
    근거: 실패 probe 최악 hist-OVL~0.25(Bhatt0.60), healthy~0.67+. 0.45 = 모호대(0.36) 위 healthy band 하단.
 2차(보조, 방향): (a) 직전 막힌 boundary(leg-up λ2 0-0.2, leg-down λ1 0.1-0.3)의 walker crossing이 0/near-impassable
    벗어나 healthy boundary 수준 회복. (b) s101 both-dir full round-trip 0->>=1 (count 아닌 0->nonzero 전이; 175cyc).
 FAIL → ADR-0020 §B(DCD structural-OP) 또는 추가 window.
사용: python reprobe_gate.py <out_root_base>  (각 cell: <base>/<tag>/{dplus,dminus}/r*/trackb_*.out)"""
import os, sys, glob, numpy as np
BASE = sys.argv[1] if len(sys.argv)>1 else "/home/san/UPDD_proj/analysis/reprobe_H16_gate"
CELLS=[("s101_wt","wt/s101","both"),("s7_wt","wt/s7","dplus"),("s199_cp4","cp4/s199","dminus")]
OVL_FLOOR=0.45; BHATT_FLOOR=0.80

def detect_NS(tag):
    ns=0
    for d in ("dplus","dminus"):
        rs=glob.glob(f"{BASE}/{tag}/{d}/r*/trackb_{d}.out")
        ns=max(ns, len(rs))
    return ns

def load(tag,d,NS):
    """walker traces(col0) + pertE by fwd state(col9)."""
    W={}; pe={s:[] for s in range(NS)}
    for w in range(NS):
        f=f"{BASE}/{tag}/{d}/r{w}/trackb_{d}.out"
        if not os.path.exists(f): continue
        a=np.loadtxt(f); a=a[None,:] if a.ndim==1 else a
        W[w]=a[:,0].astype(int)
        for s,p in zip(a[:,0].astype(int),a[:,9]):
            if 0<=s<NS: pe[(s if d=="dplus" else NS-1-s)].append(p)
    return W, {s:np.array(v) for s,v in pe.items()}

def hist_ovl(x,y,nb=60):
    if len(x)<5 or len(y)<5: return np.nan
    lo=min(x.min(),y.min()); hi=max(x.max(),y.max())
    if hi-lo<1e-9: return 1.0
    b=np.linspace(lo,hi,nb+1)
    px,_=np.histogram(x,bins=b,density=True); py,_=np.histogram(y,bins=b,density=True)
    bw=b[1]-b[0]; return float(np.sum(np.minimum(px,py))*bw)
def gauss_bhatt(x,y):
    if len(x)<3 or len(y)<3: return np.nan
    m1,s1=x.mean(),x.std(); m2,s2=y.mean(),y.std()
    if s1<1e-9 or s2<1e-9: return np.nan
    return float(np.sqrt(2*s1*s2/(s1**2+s2**2))*np.exp(-(m1-m2)**2/(4*(s1**2+s2**2))))
def roundtrips(W,NS):
    tot=0
    for st in W.values():
        hits=[]
        for s in st:
            if s==0 or s==NS-1:
                if not hits or hits[-1]!=s: hits.append(s)
        tot+=sum(1 for i in range(1,len(hits)) if hits[i]!=hits[i-1])//2
    return tot
def crossings(W,a,b):  # boundary between state a,b (a<b adjacent)
    c=0
    for st in W.values():
        for i in range(1,len(st)):
            lo,hi=min(st[i-1],st[i]),max(st[i-1],st[i])
            if lo==a and hi==b: c+=1
    return c

print("="*84); print("ADR-0021 re-probe gate (nonparametric)"); print(f"floors: hist-OVL>={OVL_FLOOR} AND gauss-Bhatt>={BHATT_FLOOR} (1차); s101 RT 0->>=1 (2차)"); print("="*84)
gate1_ok=True; s101_rt={}
for tag,label,marg in CELLS:
    NS=detect_NS(tag)
    if NS==0: print(f"\n### {label}: NO DATA (미완주?)"); gate1_ok=False; continue
    print(f"\n### {label} (marg {marg}, NS={NS})")
    for d in ("dplus","dminus"):
        W,pe=load(tag,d,NS)
        if not W: print(f"  {d}: NO DATA"); gate1_ok=False; continue
        ovls=[hist_ovl(pe.get(j,np.array([])),pe.get(j+1,np.array([]))) for j in range(NS-1)]
        bhs=[gauss_bhatt(pe.get(j,np.array([])),pe.get(j+1,np.array([]))) for j in range(NS-1)]
        ovls=[o for o in ovls if not np.isnan(o)]; bhs=[b for b in bhs if not np.isnan(b)]
        rt=roundtrips(W,NS)
        if tag=="s101_wt": s101_rt[d]=rt
        mn_o=min(ovls) if ovls else float('nan'); mn_b=min(bhs) if bhs else float('nan')
        jo=int(np.argmin(ovls)) if ovls else -1
        ok = (mn_o>=OVL_FLOOR) and (mn_b>=BHATT_FLOOR)
        gate1_ok = gate1_ok and ok
        print(f"  {d}: min hist-OVL={mn_o:.3f}@({jo},{jo+1}) min gauss-Bhatt={mn_b:.3f}  RT={rt}  1차={'✅' if ok else '❌'}")
print("\n"+"="*84)
print(f"1차(전 인접쌍 OVL>={OVL_FLOOR} & Bhatt>={BHATT_FLOOR}): {'🟢 PASS' if gate1_ok else '🔴 FAIL'}")
if s101_rt: print(f"2차(s101 RT 0->>=1): dplus RT={s101_rt.get('dplus')}, dminus RT={s101_rt.get('dminus')} → {'🟢' if all(v>=1 for v in s101_rt.values()) else '🟡(RT count는 175cyc confound; 0->nonzero 전이가 신호)'}")
print("판정: 1차 PASS → full run(H16, 전 cell auto-search) GO. FAIL → ADR-0020 §B(DCD OP) / 추가 window.")
print("="*84)
