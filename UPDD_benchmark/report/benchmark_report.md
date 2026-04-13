# UPDD Benchmark Report

## 벤치마크 방법론

본 벤치마크는 UPDD 파이프라인의 **통합 성능** (구조 예측 + 스코어링) 을 평가한다.

파이프라인: 타겟 PDB 전처리 → ColabFold 복합체 예측 → MD → QM/MM → MM-GBSA

실험 구조 (PDB) 는 **타겟 체인 전처리에만** 사용되며, 복합체 구조는 ColabFold 가 데이터셋의 바인더 서열로부터 예측한다. 따라서 본 평가에는 AF2 구조 예측 오차가 포함되어 있으며, 스코어링 단독 성능과는 구분된다. 스코어링 단독 성능을 측정하려면 `--use-experimental` 플래그를 사용한다.

### cyclic_htc 토폴로지 주의사항

데이터셋의 cyclic_htc 펩타이드는 ColabFold 로 구조 예측 후 **linear 토폴로지로 MD 를 수행** 하였다. ColabFold 는 cyclic constraint 를 강제하지 않으므로, 예측 구조는 linear 이며, 이 구조에서 head-to-tail cyclic bond 를 post-hoc 으로 형성하는 것은 비물리적이다 (N-C 거리 30+ Å). cyclic_htc 하위그룹의 스코어링 결과는 "linear 구조 기반 스코어링"으로 해석해야 하며, cyclic 구조 기반 스코어링은 AfCycDesign 통합 후 가능하다. cyclic_ss (disulfide) 는 AF2 구조에서도 자연 형성되므로 원래 토폴로지로 MD 를 수행한다.

- Dataset 레코드: 2
- 평가 가능 (valid): 1
- 실패: 0
- 누락 (진행 불가): 2

## Overall Performance

| 지표 | 값 |
|------|----|
| Spearman ρ | +nan (p=nan) |
| Kendall τ  | +nan |
| Pearson r  | +nan |
| ROC-AUC    | None |
| EF 1%      | 0.00 |
| EF 5%      | 0.00 |
| EF 10%     | 0.00 |
| **Tier**   | **N/A** |

해석: 핵심 지표는 Spearman ρ 이며, ΔG 가 낮을수록 pKd 가 높아야 하므로 이상적으로 ρ<0 이다. |ρ|>0.4 = Tier1, |ρ|>0.6 = Tier2, |ρ|>0.75 = Tier3.

## Cross-Validation

- (스킵: 샘플 수 부족 (n=1, folds=5))

## Leave-One-Target-Out

- 타겟 수: 0
- Spearman ρ 평균: +nan

## Topology Subgroup

| Topology | n | Spearman ρ | ROC-AUC | EF 10% |
|----------|---|------------|---------|--------|
| cyclic_htc | 1 | — | — | — |

## Failed / Missing Cases

2 개의 데이터 포인트가 유효 점수를 산출하지 못했습니다.

| peptide_id | status | stage |
|------------|--------|-------|
| keap1_005 | MISSING |  |
| brd4_001 | MISSING |  |