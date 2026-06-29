#!/usr/bin/env python
"""① Thermodynamic-length 윈도 배치 (곡선이 결정). 2 방법:
 (A) SUBDIVIDE-ONLY  : 측정된 각 segment를 D=-ln(BC) 기준 ceil(D/Dmax)개로 쪼갬(선형 λ). 병합 없음
                       → 외삽 위험 0(전부 측정 구간 내 interpolation). over-resolved는 유지(state↑).
 (B) EQUIPART+CAP    : 누적-D 등간격 재배치 BUT max λ-step cap → easy 구간 과병합(거대점프) 금지.
                       over-resolved에서 일부 빼서 net 절감 가능, 거대점프는 cap이 막음.
floor가 윈도수 하한 결정(Dmax=-ln floor). worst-case=cohort 4 leg MIN Bhatt(제일 막힌 seed 기준).
leg-up: λ1=0/λ2 0→0.5 (state0..6). leg-down: λ2=0.5/λ1 0→0.5 (state6..12).
정직: D 가산·segment 내 D 균등분할은 국소 smooth 가정 → re-probe overlap이 실측 검증(②). C7 불변·ΔG-unbiased."""
import os, numpy as np
G="/home/san/UPDD_proj/analysis/apex_probe_gate_20260628"; NS=13
LAM2=[0,0.1,0.2,0.4,0.45,0.475,0.5]; LAM1_PATH=[0.0,0.025,0.05,0.1,0.3,0.4,0.5]

def pertE_by_fwd(tag,d):
    pe={s:[] for s in range(NS)}
    for w in range(NS):
        f=f"{G}/{tag}/{d}/r{w}/trackb_{d}.out"
        if not os.path.exists(f): continue
        a=np.loadtxt(f); a=a[None,:] if a.ndim==1 else a
        for s,p in zip(a[:,0].astype(int),a[:,9]):
            if 0<=s<NS: pe[s if d=="dplus" else NS-1-s].append(p)
    return {s:np.array(v) for s,v in pe.items()}
def bc(x,y):
    if len(x)<3 or len(y)<3: return np.nan
    m1,s1=x.mean(),x.std(); m2,s2=y.mean(),y.std()
    if s1<1e-9 or s2<1e-9: return np.nan
    return float(np.sqrt(2*s1*s2/(s1**2+s2**2))*np.exp(-(m1-m2)**2/(4*(s1**2+s2**2))))
COHORT=[("s101_wt","dplus"),("s101_wt","dminus"),("s7_wt","dplus"),("s199_cp4","dminus")]
data={c:pertE_by_fwd(*c) for c in COHORT}
def seg_worstBC(j):
    vals=[bc(data[c].get(j,np.array([])),data[c].get(j+1,np.array([]))) for c in COHORT]
    vals=[v for v in vals if not np.isnan(v)]; return min(vals) if vals else np.nan

# leg별 (knot λ 리스트, 그 leg의 fwd pair 인덱스)
LEGS={"legup(λ2)":(LAM2,list(range(0,6))), "legdown(λ1)":(LAM1_PATH,list(range(6,12)))}
print("측정 worst-case segment D=-ln(BC):")
for name,(lam,segs) in LEGS.items():
    print(f"  {name}: "+" ".join(f"[{lam[i-segs[0]]:.3f}→{lam[i-segs[0]+1]:.3f}]D={-np.log(seg_worstBC(i)):.3f}" for i in segs))

def subdivide(lam,segs,floor):
    Dmax=-np.log(floor); knots=[lam[0]]
    for k,j in enumerate(segs):
        D=-np.log(seg_worstBC(j)); nsub=max(1,int(np.ceil(D/Dmax-1e-9)))
        a,b=lam[k],lam[k+1]
        for t in range(1,nsub+1): knots.append(round(a+(b-a)*t/nsub,4))
    return knots
def equipart_cap(lam,segs,floor,cap):
    """greedy: λ=0에서 누적 D를 Dmax씩 전진 + λ-step을 cap으로 제한. D(λ)=segment내 선형."""
    Dmax=-np.log(floor); D=[-np.log(seg_worstBC(j)) for j in segs]
    lamA=np.array(lam,float); cum=np.concatenate([[0],np.cumsum(D)]); tot=cum[-1]
    def lam_at(dval):  # 누적 D=dval 인 λ (선형보간)
        return float(np.interp(dval, cum, lamA))
    knots=[float(lamA[0])]; curD=0.0
    while knots[-1] < lamA[-1]-1e-6:
        nextD=min(curD+Dmax, tot)
        nl=lam_at(nextD)
        if nl - knots[-1] > cap:        # cap 위반 → λ-step만큼만 전진(D는 그만큼만)
            nl=round(knots[-1]+cap,4); nextD=float(np.interp(nl,lamA,cum))
        nl=round(nl,4)
        if nl<=knots[-1]+1e-4: nl=round(knots[-1]+cap,4)
        if nl>=lamA[-1]-1e-4: nl=float(lamA[-1])
        knots.append(nl); curD=nextD
        if len(knots)>60: break
    # 중복 제거 + 마지막 0.5
    out=[knots[0]]
    for x in knots[1:]:
        if x>out[-1]+1e-4: out.append(x)
    if out[-1]!=lamA[-1]: out.append(float(lamA[-1]))
    return out

for floor in (0.8,0.85):
    print(f"\n{'='*78}\nfloor={floor} (Dmax=-ln={-np.log(floor):.3f}) — 모든 인접쌍 overlap>={floor} 목표")
    for name,(lam,segs) in LEGS.items():
        A=subdivide(lam,segs,floor); B=equipart_cap(lam,segs,floor,0.12)
        # leg-down은 apex(0) 제외하고 CLI knots
        if "legdown" in name:
            A=[x for x in A if x>1e-6]; B=[x for x in B if x>1e-6]
        if A[-1]!=0.5: A[-1]=0.5
        if B[-1]!=0.5: B[-1]=0.5
        print(f"  {name}:")
        print(f"     (A)subdivide  [{len(A)}st]: {A}")
        print(f"     (B)equip+cap  [{len(B)}st]: {B}")
    # total
    for lbl,fn in (("subdivide",subdivide),("equip+cap",lambda l,s,f:equipart_cap(l,s,f,0.12))):
        up=fn(LAM2,list(range(0,6)),floor); dn=fn(LAM1_PATH,list(range(6,12)),floor)
        dn=[x for x in dn if x>1e-6]
        if up[-1]!=0.5:up[-1]=0.5
        if dn[-1]!=0.5:dn[-1]=0.5
        print(f"   => {lbl}: total/dir = {len(up)+len(dn)} (probe13 대비 {len(up)+len(dn)-13:+d})")
print("\nR-18: D-가산·균등분할은 국소 smooth 가정 → re-probe overlap이 floor 충족 실측 검증. cap=0.12로 거대점프 차단.")
