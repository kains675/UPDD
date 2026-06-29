# apex probe 재설계 grounding (2026-06-28)

## 목표
densified apex probe가 apex overlap은 고쳤으나(Bhatt 0.69→~1.0) **병목이 다른 λ-구간으로 이동**해 round-trip 미복구. 균등 overlap λ-schedule로 재설계 → cheap re-probe → full run.

## 스케줄 구성 SSOT
`utils/atm_trackB_inplace_rbfe.py::build_ats_standard_ladder` (L410-667). 호출: `scripts/trackb_inplace_rbfe_production.py::_build_combined_schedule` (forward+backward 각각 build, UWHAM merge).
- CLI `--lambda2-rampdown` → 함수 `lambda2_rampup` (leg-up λ2 knots; **0.0에서 시작, 0.5에서 끝**, [0,0.5] strictly increasing).
- CLI `--lambda1-rampdown` → 함수 `lambda1_rampdown` (leg-down λ1 knots; **(0,0.5]**, 0.5에서 끝, strictly increasing).
- leg-up: λ1=0, λ2가 knots 따라 0→0.5 (soft-core alchemical region). state 0..len(λ2knots)-1. 마지막 = apex(λ1=0/λ2=0.5).
- leg-down: λ2=0.5, λ1이 knots 따라 →0.5. apex 다음부터. 마지막 = coupled apex(λ1=λ2=0.5).
- total states/dir = len(λ2_rampup) + len(λ1_rampdown). (apex λ1=0/λ2=0.5는 leg-up이 1회 배치, 중복 없음)
- **ΔG-unbiased**: 내부 λ 재배치는 Kirkwood path independence로 ΔG 불변(endpoints + soft-core canon C7 불변). overlap/robustness lever일 뿐 (ranking-only, R-11). **C7 soft-core(umax/ubcore/acore) 불가침**.

## probe(실패) 스케줄 = 13 states/dir
- leg-up λ2 = [0, 0.1, 0.2, 0.4, 0.45, 0.475, 0.5]  (7 states: 0-6, state6=apex)
- leg-down λ1 = [0.025, 0.05, 0.1, 0.3, 0.4, 0.5]  (6 states: 7-12)
- forward state별 (λ1,λ2): s0(0,0) s1(0,0.1) s2(0,0.2) s3(0,0.4) s4(0,0.45) s5(0,0.475) s6(0,0.5=APEX) s7(0.025,0.5) s8(0.05,0.5) s9(0.1,0.5) s10(0.3,0.5) s11(0.4,0.5) s12(0.5,0.5)
- dminus = forward tuple 전체 reverse → dminus state i = forward state (12-i).

## 측정 adjacent overlap (Bhattacharyya, pertE=col9; analysis/apex_probe_gate_20260628/overlap.py)
원본 11-state apex(구 bottleneck) Bhatt~0.69 / MBAR O~0.05. over-resolved 쌍은 ~0.95-1.0.
- s101/dplus: 0-1:0.61 1-2:0.81 2-3:0.97 3-4:1.00 4-5:1.00 5-6:1.00 6-7:0.91 7-8:0.91 8-9:0.84 **9-10:0.60** 10-11:0.95 11-12:0.97
- s101/dminus: 0-1:0.97 1-2:0.93 **2-3:0.61** 3-4:0.86 4-5:0.89 5-6:0.93 6-7:1.00 7-8:1.00 8-9:1.00 9-10:0.95 10-11:0.65 11-12:0.88
- s7/dplus:   0-1:0.66 1-2:0.85 2-3:0.97 3-4:1.00 4-5:1.00 5-6:1.00 6-7:0.91 7-8:0.91 8-9:0.83 **9-10:0.73** 10-11:0.96 11-12:0.96
- s199/dminus:0-1:0.98 1-2:0.93 **2-3:0.74** 3-4:0.86 4-5:0.91 5-6:0.87 6-7:1.00 7-8:1.00 8-9:0.99 9-10:0.98 10-11:0.94 11-12:0.77

## round-trip (gate.py): 원본 marginal zero-RT가 densify 후에도
- s101 dplus RT=0, dminus RT=0 (미복구); s7 dplus RT=0 (미복구); s199 dminus RT=1 (복구).
- (probe=175cyc vs 원본=800cyc → RT 개수는 cycle 비례 confound. overlap이 cycle-독립 판별자.)

## 도출된 병목 귀속 (가설 — 검증 대상)
dminus state i = fwd (12-i)이므로 dminus 병목을 fwd로 환산:
- dminus 2-3 = fwd pair (9,10). dminus 10-11 = fwd pair (1,2). dminus 11-12 = fwd (0,1).
- **#1 병목 = λ1 leg-down 0.1→0.3** (fwd state9→10, Δλ1=0.2): dplus 9-10=0.60 + dminus 2-3=0.61/0.74 양방향 확인. DOMINANT.
- **#2 = λ2 leg-up 0→0.1** (fwd state0→1, Δλ2=0.1): dplus 0-1=0.61/0.66.
- (minor: λ2 0.1→0.2 fwd(1,2) dminus 10-11=0.65; λ1 0.4→0.5 fwd(11,12) s199 dminus 11-12=0.77.)
- ⚠️ **λ2 leg-up 0.2→0.4** (fwd state2→3, Δλ2=0.2)는 Bhatt 0.97로 **멀쩡** — soft-core가 leg-up 큰 λ2 step을 싸게 만듦. **모든 0.2 gap이 병목이 아니다.** (leg-down λ1 step은 비싸고, leg-up λ2 step은 쌈 = 비대칭.)

## 재설계 핵심 질문
1. 병목 귀속(위 가설)이 raw .out으로 맞나? 방향별로 정확히 어느 물리 λ-구간이 under-resolved인가?
2. 어디에 window를 넣고(leg-down λ1 0.1-0.3 우선, leg-up λ2 0-0.1) 어디서 뺄 수 있나(over-resolved apex λ2 0.4/0.45/0.475? λ1 0.025/0.05?)?
3. window 넣으면 병목이 **또 이동(whack-a-mole)**하지 않나? 균등 overlap 달성 후보는?
4. 총 state 수 ↔ GPU cost trade-off (13→몇 개까지 허용? full run 7-GPU-day 기준).
5. 변위는 전 cell `--auto-search-displacement --accept-sep-nm 1.5` 통일(baseline, C6 fix).

## 제약
- C7 soft-core canon 불가침. ΔG-unbiased 유지(내부 λ 재배치만, endpoints 불변). ranking-only(R-11, magnitude/sign 무주장). frozen XML/charge SSOT 불가침. no-AI-trail. SciVal 게이트(R-4) 필수.
