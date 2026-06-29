#!/usr/bin/env python
"""ADR-0020 §A 2차: apex adjacent overlap (cycle-독립) — densified가 floor 0.10 넘겼나.
pertE=col9(idx9). state s의 pertE 분포 = col0==s 행들의 col9 집계. 인접 state Bhattacharyya=overlap proxy.
densified 13-state, apex=state6 → 핵심 인접쌍 5-6,6-7. 원본(11-state apex5)에서 5-6/6-7 MBAR O~0.05, Bhatt~0.68-0.71."""
import os, numpy as np
G="/home/san/UPDD_proj/analysis/apex_probe_gate_20260628"
NS=13; APEX=6

def state_pertE(tag,d):
    """{state: array of pertE(col9) sampled in that state}."""
    pe={s:[] for s in range(NS)}
    for w in range(NS):
        f=f"{G}/{tag}/{d}/r{w}/trackb_{d}.out"
        if not os.path.exists(f): continue
        a=np.loadtxt(f)
        if a.ndim==1: a=a[None,:]
        for s,p in zip(a[:,0].astype(int), a[:,9]):
            if 0<=s<NS: pe[s].append(p)
    return {s:np.array(v) for s,v in pe.items()}

def bhatt(x,y):
    """Bhattacharyya coeff of two gaussians (unequal var)."""
    if len(x)<3 or len(y)<3: return float('nan')
    m1,s1=x.mean(),x.std(); m2,s2=y.mean(),y.std()
    if s1<1e-9 or s2<1e-9: return float('nan')
    return float(np.sqrt(2*s1*s2/(s1**2+s2**2))*np.exp(-(m1-m2)**2/(4*(s1**2+s2**2))))

CELLS=[("s101_wt","wt/s101","both"),("s7_wt","wt/s7","dplus"),("s199_cp4","cp4/s199","dminus")]
ORIG_APEX_BHATT=0.69  # 원본 11-state apex(5-6,6-7) median Bhatt (MBAR O~0.05 동반)
print("="*82)
print("ADR-0020 §A 2차 — apex adjacent overlap (Bhattacharyya, pertE col9) RELATIVE 비교")
print("⚠ floor 0.10은 MBAR O 기준(별도). 여기선 densify로 apex Bhatt가 원본(~%.2f) 대비 ↑ 했나만 본다."%ORIG_APEX_BHATT)
print("densified apex=state6 → 인접쌍 5-6,6-7. (원본 apex Bhatt~0.68-0.71 = 사다리 최저, over-resolved쌍은 0.95)")
print("="*82)
for tag,label,marg in CELLS:
    print(f"\n### {label} (marginal {marg})")
    for d in ("dplus","dminus"):
        pe=state_pertE(tag,d)
        adj={}
        for s in range(NS-1):
            adj[(s,s+1)]=bhatt(pe[s],pe[s+1])
        vals=[(k,v) for k,v in adj.items() if not np.isnan(v)]
        if not vals: print(f"  {d}: NO DATA"); continue
        minpair=min(vals,key=lambda kv:kv[1])
        apex_lo=adj.get((APEX-1,APEX),float('nan')); apex_hi=adj.get((APEX,APEX+1),float('nan'))
        apex_min=min(apex_lo,apex_hi)
        marg_hit=(marg=="both") or (marg==d)
        m="  <<MARGINAL" if marg_hit else ""
        trend = "↑개선" if apex_min>ORIG_APEX_BHATT+0.05 else ("≈동일" if apex_min>ORIG_APEX_BHATT-0.05 else "↓악화")
        print(f"  {d}: apex[5-6]={apex_lo:.3f} apex[6-7]={apex_hi:.3f} (apex-min={apex_min:.3f} vs 원본~{ORIG_APEX_BHATT}: {trend}) | min-adj-ladder={minpair[1]:.3f}@{minpair[0]}{m}")
        if marg_hit:
            print("     full adj Bhatt: "+" ".join(f"{k[0]}-{k[1]}:{v:.2f}" for k,v in sorted(adj.items())))
print("\n"+"="*82)
print("해석: apex-min Bhatt가 원본~0.69보다 뚜렷히 ↑(예 >0.8)면 densify가 apex overlap 개선(round-trip은 cycle부족).")
print("apex-min이 여전히 ~0.69 근처면 densify로도 apex 병목 미해소 → 분기 B(DCD structural-OP).")
print("※ 절대 MBAR O>=0.10 확증은 trackb_uwham_postprocess 필요(full run). 본 Bhatt는 상대 cheap-proxy.")
print("="*82)
