# UPDD v0.6 Refactor — Physical/Chemical Validity Recovery (Rev. 2)

**작성일:** 2026-04-19
**개정:** Rev. 2 (2026-04-19 세션 발견 반영)
**이전 Rev:** Rev. 1 (2026-04-19 초안)
**이전 버전:** v0.5.x (SSOT refactor — SUCCESS/MARGINAL tier classification)
**목표 배포:** v0.6.x → v0.7 Calibration 진입 조건 충족

---

## 📌 Mission Statement

v0.5 refactor 에서 SSOT 구조가 정립되었고, v0.6 에서 **물리/화학적 correctness recovery** 가 진행 중이다. 2026-04-19 세션에서 다수의 optimization 옵션이 시도되었고, 아래 **엄격한 경계 조건** 이 확인되었다:

### 확인된 물리적/기술적 경계 (변경 불가)

1. **RTX 5070 Ti 16GB VRAM** 에서 `n_qm ≥ 466` 의 **DF-J/K 불가능**
   - gpu4pyscf 내부 C++ kernel (libgint, libgvhf_rys) 이 managed memory 우회
   - Application-layer allocator override 로 해결 불가
2. **wB97X-D3 (range-separated hybrid)** 는 SR+LR K 를 두 번 build → 메모리 2×
3. **Direct-SCF mode 에서 cycle 당 130-164초** (n_qm=466 기준)
   - B2/B3/B5 SCF parameter tuning 효과 미미 (164s → 163s)
   - chkfile propagation 도 rks_lowmem 미호환으로 제한적
4. **Design-specific HOMO-LUMO near-degeneracy** 존재 (w4_1_s2 계열)
   - Metastable SCF convergence → unphysical energy 위험
   - Init guess (huckel, minao) 에 민감
5. **Budget 500 atom 하드 제약** — 과학적 근거로 축소 불가
   - 6WGN DxxGQ motif, α3-groove contact residue 는 Zhang 2020 Cell Chem Biol 명시
   - `enforcement="hard_fail"` 하드코딩 → n_qm > 500 design 은 partition_error

### v0.6 의 재정의된 목표

Rev. 1 에서 "15× 가속" 을 기대했으나 **실측에서 거의 불가능함이 확인**. Rev. 2 의 현실적 목표:

- ❌ NOT: DFT QM/MM 를 대규모 가속
- ✅ YES: 물리/화학적 correctness 복구 (charge state consistency, topology integrity)
- ✅ YES: Pre-QM filtering 으로 **실패할 design 을 사전 탈락**
- ✅ YES: MM-GBSA 와 DFT QM/MM 의 **hybrid 활용** 으로 정보량 확보
- ✅ YES: v0.7 Calibration 진입 조건 충족 (최소 3 design × 5 snap converged)

### 하드웨어 업그레이드 경로 명시

현재 설계는 **RTX 5070 Ti 16GB 고정** 전제. 다음 업그레이드 시 자연스럽게 확장:
- **RTX 5090 32GB**: DF auto-mode 가 자동으로 on → 2-3× 가속
- **H100 80GB**: n_qm 1000+ 가능, replicate 확장

업그레이드 전까지는 **현 하드웨어의 물리적 한계 수용** + **영리한 pipeline 설계** 로 진행.

---

## 🎯 UPDD v0.6 Constitution

### Core Principles (Rev. 2)

1. **물리적 의미 우선** (Physical Meaning Supreme)
   파이프라인의 모든 출력 값은 정의된 물리/화학적 의미를 가져야 한다.

2. **Workflow 보존** (Preserve Workflow)
   UPDD 의 파이프라인 단계 순서와 주요 인터페이스는 유지.

3. **Raw 보존 + Derived 추가** (Raw Preserved, Derived Added)
   기존 필드 삭제 금지. 새 의미는 새 필드로.

4. **Fail Fast on Ambiguity**
   물리/화학적 의미가 불명확한 조건에서 explicit failure.

5. **Diagnostic Transparency**
   모든 결과 JSON 은 계산 조건을 명시 (charge state, QM region mode, SCF level).

6. **Ranking-Invariance Regime** (신규, Rev. 2)
   v0.7 Calibration 진입 전까지 **ranking invariance 심사 기준** 적용.
   - 모든 최적화는 design 전체에 일괄 적용 (bias 방지)
   - Absolute accuracy 향상 주장 금지 (ranking 만 validated)
   - Calibration 돌입 시 regime 재선언 (`ranking_invariance_regime.md` 참조)

7. **Hardware-Aware Design** (신규, Rev. 2)
   현재 하드웨어의 물리적 경계를 인정하고 설계.
   - DF 는 VRAM ≥ 20GB 일 때만 자동 on
   - n_qm > 500 design 은 pre-filter 로 탈락
   - 업그레이드 시 코드 수정 없이 확장 (env var + auto-detect)

### Absolute Rules

- **R-1 ~ R-3** (Coder 기본): 변수 순서, 기계적 치환 금지, import 무결성
- **R-4**: SciVal ❌ 판정 시 코드 배포 금지
- **R-5**: Runner 검증 없이 완료 선언 금지
- **R-6**: Path 재진단 전 task 종료 금지
- **R-7**: 탈락 데이터/결과는 `_archive/` 로 이동 (삭제 금지)
- **R-8**: 모든 수정은 `UPDATE.md` 상단에 기록
- **R-9**: Git 커밋/푸시 전 사용자 허가 필수 (커맨드는 제시, 실행 금지)
- **R-10**: 완료된 upgrade 는 `UPGRADE.md` 에 체크 표시
- **R-11** (신규, Rev. 2): **Ranking regime 외 주장 금지**. Absolute ΔG 값의 정량 의미는 v0.7 calibration 완료 전까지 주장 불가. 결과 리포트에 "(ranking-only)" 명시.
- **R-12** (신규, Rev. 2): **설정 속성 효과 소스 검증**. `mf.with_df.max_memory` 같은 "읽히지 않는 속성" 존재. 새 환경변수/config 추가 시 upstream grep 으로 참조 확인 필수.
- **R-13**: 예약 (향후 Long SCF nohup rule 자리).
- **R-14**: 예약.
- **R-15** (신규, 2026-04-19): **Charge Declaration SSOT Enforcement (magnitude + parity + rationale)**.
  `target_card.json` 의 `target_iso_net_charge` 및 `binder_net_charge` 는 QM 실행 시 snapshot PDB 로부터 chemistry-true 값과 **magnitude-level 일치 필수**. 1-bit parity check 외에 추가 magnitude check + rationale auto-regenerate. Mismatch 시 `ChargeDeclarationMismatch` 로 fail-fast; JSON diag 에 chemistry-computed charge + declared charge 양쪽 기록. `target_iso_charge_rationale` 텍스트는 code-generated only (수동 편집 금지). 기존 `OddElectronError` (parity-only) 는 magnitude 와 독립 safety net 으로 유지.
- **R-16** (신규, 2026-04-19): **Binder Chemistry SSOT Enforcement**.
  Binder net charge 는 snapshot PDB 로부터 residue identity + atom name pattern (ARG HH11-HH22, ASP OD-HD, GLU OE-HE, LYS HZ, HIS HD1/HE2 tautomer) 을 기반으로 runtime 계산한다. `target_card.binder_net_charge` 는 "expected hint" 필드일 뿐이며 runtime PDB 가 ground truth. Cyclic 여부는 LEU1 (또는 첫 잔기) N-H count 와 마지막 잔기 OXT 유무로 자동 감지. 불일치 시 `BinderChargeMismatch` 로 fail-fast.

---

## 🤖 5-Agent 구조 (Path 신설, Rev. 2 에서 역할 강화)

| Agent | 역할 | 권한 | 호출 타이밍 |
|---|---|---|---|
| **Archi** | 아키텍트, 전략, 조율 | Read-only, 지시 권한 | 전체 계획, 작업 간 조율 |
| **SciVal** | 사전 과학적 검증 | Read-only, Veto | 코드 변경 **전** 검증 |
| **Coder** | 코드 구현, 수정 | Write | Archi 지시 + SciVal 승인 후 |
| **Runner** | 실행, 로그 수집, 기본 검증 | Bash, Read | Coder 완료 후 |
| **Path** | 결과의 물리/화학/병리학적 진단 | Read-only + Bash (read), Advisory | Runner 완료 후 |

### Rev. 2 에서 강화된 Path 의 역할

Rev. 1 에서 Path 는 "결과 진단" 역할이었다. 2026-04-19 세션에서 드러난 실패 패턴이 **다양하고 복잡**하므로 Rev. 2 에서 Path 의 진단 범위가 확장된다:

**Path 가 진단할 실패 카테고리 (Rev. 2)**:

1. **수렴 실패** (max_cycle 초과)
2. **Metastable 수렴** (converged=True 지만 unphysical energy)
3. **HOMO-LUMO near-degeneracy** (w4_1_s2 계열)
4. **Budget violation** (n_qm > 500)
5. **QM region inconsistency** (topology mode 전 snap 비교)
6. **Charge correction 불일치** (total_protons parity 변화)
7. **SCF init_guess 발산** (huckel 에서 ΔE > 1000 Ha)
8. **Library 레벨 silent failure** (with_df.max_memory 류 무시 속성)

각 카테고리마다 Path 의 **증상 → 원인 → 권고** template 적용.

---

## 🎯 v0.6 작업 우선순위 (Rev. 2, 현실 기반 재정렬)

### 원칙

**"Do less, but ensure correctness"**. Rev. 1 의 6 Phase 를 **4 Phase 로 축소**하고, 실패한 시도를 Phase 에서 제거.

### Phase 1: ✅ **완료** — Topology-based QM Region Selection

Rev. 1 에서 이미 구현됨 (커밋 `14b0295 feat: v0.6.0`).

- Snapshot-invariant n_qm
- `target_card.json` 스키마 확립
- Default mode = `"topology"`

### Phase 2: ✅ **완료** — 3-term Supermolecular Decomposition

Rev. 1 에서 이미 구현됨.

- `qm_int_kcal_frozen` 필드 (3-term 분해)
- `interaction_kcal` 보존 (embedding coupling)

### Phase 3: 🟢 **Stage 2 complete (2026-04-19), awaiting baseline recompute** — Pre-SCF Charge Guard + R-15/R-16 + schema 0.6.4 partitioning

**Stage 1 (merged `4034ee6`, 2026-04-19)**: `compute_qm_net_charge_topology()` + `compute_binder_chem_charge()` 헬퍼, `run_qmmm.py` R-15 magnitude guard + R-16 binder SSOT, `utils/audit_charge_consistency.py` (dry-run default), result JSON schema 확장 (`qm_net_charge_declared/computed`, `binder_charge_declared/computed`, `charge_consistency_audit`), sub-issue #3 L1579 fix.

**Stage 2 (code complete 2026-04-19)**: Option B3 partitioning 재설계.
- `target_cards/6WGN.json` schema 0.6.3 → 0.6.4: `target_residues_mode=whole`, `max_n_qm 500→600`, `whole_residue_exceptions` 제거, `numbering_convention={source: MD, md_to_crystal_offset: 3}` 신규, `target_iso_charge_rationale` auto-generated (MD/crystal 페어 포함), `binder_net_charge=-1` (R-16 hint runtime 일치).
- `utils/run_qmmm.py::partition_qmmm()`: `target_residues_mode` dispatch (whole / sidechain_from_cb). whole mode 에서는 Cα-Cβ cut 없이 peptide-bond 경계만 link H (20 H 예상).
- `utils/run_qmmm.py::run_qmmm_calc()`: `[DIAG] numbering_convention` event, 0.6.3 backward-compat warning, `[DIAG] qm_region` 에 `target_residues_mode` 추가.
- `utils/generate_target_card.py`: `generate_target_iso_rationale()` (MD/crystal pair + breakdown + tautomer notes), `validate_numbering_convention()`, schema 0.6.4 writer, CLI `--md-to-crystal-offset` / `--numbering-source` / `--target-residues-mode`.
- `utils/utils_common.py::load_target_card()`: 0.6.4 schema 지원 (0.6.0-0.6.3 backward-compat 유지).
- `tests/test_stage2_partitioning.py` 신규 (8 tests). `tests/` 전체 44 passed / 1 skipped.
- Partition-only smoke test on w4_1_s2/snap01: `{n_qm: 546, n_link: 20, binder: -1, target: -3, total: -4, declared: -4, r15_passed: true}`.

**Scope (Stage 3, 대기)**: baseline recompute (full QM/MM 20 snap × 4 design) — SciVal precondition 으로 user 승인 후 실행.

**재평가 조건 (Stage 2 완료 판정, 충족)**: R-15/R-16 guard 유지, schema 0.6.4 load 성공, partition_qmmm whole path n_qm 500-640 범위, 모든 테스트 pass. 남은 baseline recompute 완료 시 Phase 3 전체 🟢 로 확정.

### Phase 4 (확장): 🔴 **진행 중** — Pre-QM Budget Enforcement + Charge Parity Filter

**원 문제**: w4_0_s2 (n_qm=578), w1_0_s1 (n_qm=539-561) 가 budget 500 초과 → partition_error 로 **QM/MM 진입 전 무조건 실패**.

**Rev. 2 확장 (2026-04-19)**: Phase 3 발견으로 scope 확장됨. 기존 n_qm filter + **charge parity + chemistry-declared 일치 filter** 추가. AF2 output + target_card 로 `compute_binder_chem_charge()` (R-16) 및 `compute_qm_net_charge_topology()` (R-15) 을 돌려 parity-incompatible / magnitude-mismatch design 을 MD 진입 전 탈락.

**현재 정책**: `enforcement="hard_fail"` (SciVal 승인)

**v0.6 에서 할 일**:

1. **탈락 design 조기 감지** — MD stage 이전에 AF2 output 에서 `estimate_n_qm()` 사전 계산
2. **Charge-aware pre-filter** — 같은 n_qm 이어도 binder identity 에 따라 parity 다름 → 각 design 의 binder residue 분석 → declared vs computed charge 자동 일치 검사 → 불일치 design 은 pre-filter 탈락
3. **Pipeline 레벨 필터** — n_qm > 500 design 은 MD 에서 제외 (MD time 낭비 방지)
4. **Logging 강화** — 탈락 design 목록을 `outputs/*/rejected_designs.json` 에 기록
5. **탈락 이유 분류**:
   - `n_qm_exceeded` (hardware limit)
   - `chemistry_parity_incompatible` (R-15/R-16 neg-mag guard) — 신규
   - `contact_residue_count` (과학적 이슈)
   - `topology_error` (partition 실패)

**수용 기준 (Acceptance Criteria)**:
- MD 진입 전 탈락 design 확인됨
- 탈락 이유가 `rejected_designs.json` 에 기록됨
- Runner 가 "N designs passed budget+charge, M designs rejected" 리포트
- Phase 3 + Phase 4 공동 merged 이전 Phase 5 BLOCK

**구현 위치**: `utils/pre_qm_filter.py` (새 모듈)

### Phase 5 (New): 🟡 **계획** — MM-GBSA + DFT Hybrid Workflow (Phase 3+4 precondition)

**Precondition (2026-04-19 강화)**: Phase 3 (R-15 charge guard) + Phase 4 (charge parity filter) 모두 merged 상태여야 Phase 5 착수 가능. 근거: Representative DFT snap 의 wrong-charge 는 hybrid weighted average 의 DFT term 을 오염 → MM-GBSA term 의 long-range ensemble benefit 을 상쇄.

**목표**: Compute budget 현실을 수용하면서 ensemble 정보량 확보.

**현재 pipeline**:
```
MD 100 ns → 5 snapshot 추출 → DFT QM/MM 5회 → median/stdev
```

**Hybrid pipeline**:
```
MD 100 ns
  → 20 snapshot 추출
  → MM-GBSA on all 20 (값싼 ensemble coverage)
  → Cluster → 2-3 representative
  → DFT QM/MM on representatives only
  → Weighted average: ΔG = w1 · MM-GBSA + w2 · (DFT - MM-GBSA correction)
```

**왜 이 접근이 타당한가**:

1. **MM-GBSA** 는 이미 UPDD 에 존재 (`run_mmgbsa.py`) 하지만 ncAA 특수 물리 부정확
2. **DFT** 는 정확하지만 비쌈 — snapshot 당 55분
3. **Hybrid** 는 MM-GBSA 로 ensemble 주요 부분 cover, DFT 로 ncAA correction 적용
4. **ΔΔG ranking 관점**: MM-GBSA 의 canceling error + DFT 의 ncAA sensitivity = 최적 조합

**DFT compute 절감**:
- 현재: 5 snap/design × 55min = 275 분/design
- Hybrid: 2-3 repr/design × 55min = 110-165 분/design
- **40% 감소** (정확도는 비슷하거나 향상)

**Coder 지시 세부사항**:

1. `run_qmmm.py` 에 representative-selection 옵션 추가
   - `--representatives` flag (default 2-3, `None` → 모든 snap)
2. `utils/cluster_snapshots.py` 신규 모듈 작성
   - Input: 20 snapshot PDB
   - Algorithm: RMSD-based k-medoids (k=2 or 3)
   - Output: representative snap ID 리스트
3. `rank_results.py` 에 hybrid integration
   - MM-GBSA (20 snap) + DFT (2-3 repr) → weighted ΔG
   - Uncertainty: MM-GBSA stderr + DFT variance

**SciVal 검증 필요**:
- Weighted average 방식의 science validity
- MM-GBSA 의 ncAA parameter 정확성 (ff14SB + GAFF 잔존)
- Representative selection 의 ensemble coverage

**수용 기준**:
- Hybrid pipeline 이 Compstatin pilot 에서 실행 성공
- MM-GBSA-only 대비 ΔΔG 개선 확인
- SciVal 승인
- **Phase 3 + Phase 4 merged 이전 BLOCK (2026-04-19 강화)**

### Phase 6 (New): 🟡 **계획** — Level Shift (B10) 완전 구현

**배경**: w4_1_s2 snap04 의 +1416 kcal/mol metastable, snap01 의 cycle 76+ plateau — HOMO-LUMO near-degeneracy 증상.

**SciVal 조건부 승인 (7 조건)**:

1. **Adaptive 점감** (V=0.25 → 0) via 2-stage kernel
2. **모든 design 일괄 적용** (ranking bias 방지)
3. **`conv_tol_grad=1e-5` 명시** (identity gate)
4. **`mf.stability()` 병행** → instability snap `_archive/` 격리
5. **Phase 4 isolation SCF 도 동일 V 적용** (method consistency)
6. **w4_1_s2 snap04 regression test** (+1416 kcal/mol 재현 여부)
7. **JSON diag 필드** (homo_lumo_gap, level_shift_V_final, converged_grad_norm)

**구현 타이밍**: 현재 실행 중인 w4_1_s2/snap01 결과 확인 후 결정
- **수렴 + 물리적** → B10 불필요
- **max_cycle 미수렴** → B10 필수
- **Wrong state converged** → B10 + stability gate 필수

### Deprecated (구현하지 않음, 향후 hardware 업그레이드 후 재고)

- **DF-J/K on 16GB VRAM** (Phase 5 Rev. 1) — 확인된 불가능
- **rks_lowmem.RKS 전환** — QMMMSCF 구조적 비호환
- **Budget 500 → 300 축소** — 과학 convergence 위배
- **Huckel default init_guess** — w4_1_s2 발산 확인

### Opt-in (환경변수로 활성화)

- **B6 huckel init_guess** — `UPDD_SCF_INIT_GUESS=huckel` 명시 시만
- **Unified Memory cascade** — 진단 용도로만, 성능 개선 없음
- **B7 mo_coeff warm-start** — T-Rank-1 smoke test 통과 후 opt-in

---

## 🖥 하드웨어 최적화 Guidelines (Rev. 2)

**Current Target**: Ryzen 9800X3D (8C/16T, 96MB L3) + RTX 5070 Ti 16GB + 32GB RAM

### VRAM 예산 (확정치)

| 항목 | VRAM 예산 |
|---|---|
| OS + display | ~1 GB |
| CUDA context + cupy mempool | ~1 GB |
| DFT electron density + grid | 2-4 GB (n_qm 400-500) |
| Direct-SCF integrals (on-the-fly) | 3-5 GB |
| Runtime buffer | ~1 GB |
| **Total 사용 가능** | **~10-13 GB** |
| **DF CDERI storage** | **14-140 GB (불가)** ❌ |

### CPU 병렬화

- PySCF OpenMP: `lib.num_threads(os.cpu_count())` 유지
- Direct-SCF 의 Fock build 는 8 core 활용 가능
- HDF5 I/O 병목 주의 (chkfile write/read)

### RAM (32GB)

- MM charges + pinned memory: ~4 GB
- 20 snapshot 동시 로드 시: ~6 GB
- Runtime buffer: ~4 GB
- **여유**: ~16 GB (future replicate expansion 용)

### Storage

- 활성 run: SSD (1TB) — chkfile, intermediate
- Archive: HDD (ExpDATA 4TB) — `_archive/`
- Trajectory (DCD 500MB+): 생성 직후 HDD 로 이동

### 환경변수 (현재 구현된 control point)

```bash
# SCF parameters
UPDD_SCF_CONV_TOL=1e-7              # B3 (ranking-safe)
UPDD_SCF_DIRECT_TOL=1e-12           # B5
UPDD_SCF_INIT_GUESS=minao           # B6 (huckel opt-in)
UPDD_DFT_GRIDS_LEVEL=2              # B2

# GPU / Memory
UPDD_MIN_AO_BLKSIZE=64              # B8 (gpu4pyscf import 전)
UPDD_VRAM_POOL_FRACTION=0.65        # B9 (cupy mempool)
UPDD_DF_AUTO_VRAM_THRESHOLD_GB=20   # DF auto threshold
UPDD_QM_MAX_MEMORY=16000            # (레거시, 참조만)

# Deprecated / No-op (참조는 유지)
UPDD_DF_JK_MAX_MEMORY_MB=4000       # gpu4pyscf 가 무시, 명시성 용도
```

---

## 🔁 진단 → 수정 → 재검증 루프 (Rev. 2 updated)

```
┌──────────────────────────────────────────────────────────────────┐
│                                                                  │
│   [새 요청]                                                        │
│      ↓                                                            │
│   ┌──────┐                                                        │
│   │ Archi│ ─── 계획 수립 ──── → SciVal 사전 검증                   │
│   └──────┘                       ↓                                │
│      ↑                        SciVal                              │
│      │                           │ 승인/거부                       │
│      │                           ↓                                │
│      │ 루프                   Coder — 코드 수정 + UPDATE.md        │
│      │ (최대 3회,                ↓                                │
│      │  그 이상은                Runner — 실행 + 기본 검증          │
│      │  user 상의)                ↓                                │
│      │                         Path — 물리/화학/병리학적 진단       │
│      │                           │                                │
│      │      ┌────────────────────┤                                │
│      │      │ 진단 결과                                            │
│      │      ↓                                                      │
│      │   🟢 통과                                                    │
│      │      → user 보고 → git 커맨드 제시 → 허가 대기               │
│      │                                                            │
│      │   🟡 조건부 통과                                              │
│      │      → user 상의 → 계속 진행 or 추가 fix                     │
│      │                                                            │
│      │   🔴 실패                                                    │
│      │      → Archi 로 복귀, refactor 방향 재설계                   │
│      └──── ↑                                                      │
│            (루프 재진입)                                             │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘

Rev. 2 루프 제한:
- Path 가 동일 증상을 3회 연속 🔴 → user 에게 에스컬레이션
- SciVal 이 동일 제안을 2회 이상 거부 → 제안 폐기
- Hardware 물리 한계 도달 → user 와 "수용 or 업그레이드" 결정
```

---

## 🧪 Calibration (v0.7) 진입 조건 (Rev. 2, 엄격화)

### 필수 조건

v0.7 Compstatin/Sa-D3 calibration 에 진입하려면 v0.6 에서 다음을 달성:

1. **최소 3 design × 5 snap = 15 snap 이 `converged=True`**
2. **각 design 의 5 snap 의 stdev < |median| × 0.3** (unified ensemble 기준)
3. **Design 간 ΔΔG_UPDD 순위가 v0.5 archive 와 Kendall τ ≥ 0.9**
4. **모든 converged snap 의 `binding_int_kcal` 이 물리적 범위 내** (Path 진단 🟢)
5. **Pre-QM filter 적용 — budget violation 은 0**
6. **모든 snap 이 `charge_consistency_audit: passed` (R-15/R-16 guard)** — 신규, 2026-04-19
7. **B10 level_shift 가 필요한 design 은 일괄 적용 완료**

### v0.7 에서 추가될 것

1. **WT counterfactual mode** — ncAA vs WT pair 비교
2. **1-Me-Trp registry** — Compstatin 용
3. **Post-hoc calibration layer** — `utils/calibrate_results.py`
4. **Config version 관리** — `calibration/v0.7.0-pilot3.json`

### Rev. 2 의 Scope 축소

Rev. 1 에서 3 pilot pair (Compstatin, Sa-D3, 15-IgBP) 계획. Rev. 2 에서는 **2 pair 축소**:
- **Compstatin (Val4 → Trp4 → 1-MeW4)**: 3 variants, 1 target
- **Sa-D3 (N-Me → ΔN-Me)**: 2 variants, 1 target
- **15-IgBP**: v0.8 로 연기

이유: Current throughput 기준 2 pair × 5 snap × 3 variants = 30 calc ≈ 30시간. 현실적.

---

## 📋 First Task for Archi (Rev. 2)

시작 지점:

1. **현재 실행 중인 w4_1_s2/snap01 결과 확인**
   - 위치: `outputs/6WGN_cyclic_htc_NMA_10_20-25/qmmm_results/*snap01*.json`
   - 수렴 여부, energy 값, 물리적 타당성
2. **Path 에 첫 진단 요청**
   - 입력: 지금까지 completed snap 들 + log
   - 기대 output: v0.6 에서 Phase 4 (Pre-QM filter) / Phase 5 (Hybrid) / Phase 6 (B10) 의 우선순위
3. **Phase 4 (Pre-QM Budget Enforcement) 착수**
   - `utils/pre_qm_filter.py` 스케치
   - `estimate_n_qm()` 함수 위치 결정
   - Runner 가 AF2 output 에서 pre-check 실행
4. **w4_0_s2, w1_0_s1 영구 제외 결정**
   - `_archive/v0.5_budget_exceeded/` 로 이동
   - `rejected_designs.json` 생성
5. **커밋 준비**: 오늘 세션 모든 uncommitted 변경사항
   - utils/run_qmmm.py (번들, cascade, DF auto)
   - utils/qmmm_proactive_oom.py (v0.3.5)
   - utils/cycle_benchmark.py (regex)
   - UPDATE.md (여러 엔트리)
   - SciVal memory files

---

## 🔄 UPDATE.md / UPGRADE.md 관리

### UPDATE.md entry format (Rev. 2)

```markdown
## [v0.6.x] YYYY-MM-DD — Brief Title

**Agent**: Coder (또는 해당 agent)
**Category**: [Bug Fix | Feature | Refactor | Performance | Documentation | Bundle]
**Regime**: ranking-only | absolute-accuracy | methodology

### 변경 사항
- `path/to/file.py::function()`: 변경 내용 (L123-145)

### 이유
과학적/기술적 근거. 연관 SciVal verdict 참조.

### 검증
- [ ] SciVal 승인 (verdict 파일: `sci_memory/xxx.md`)
- [ ] Runner 실행 검증 (regression: N cases pass)
- [ ] Path 재진단 🟢

### Regime impact
- Ranking invariance: [preserved | possibly violated — explain]
- Absolute accuracy: [unchanged | improved | degraded]

### 관련
- UPGRADE.md Phase N 완료 표시
- 관련 실패: RECENTv06.md Section X
```

### UPGRADE.md roadmap (Rev. 2)

```markdown
# UPDD Upgrade Roadmap (Rev. 2)

## v0.6 — Physical/Chemical Validity Recovery (진행 중)

### Phase 1: Topology-based QM Region
- [x] `partition_qmmm()` topology mode
- [x] `target_card.json` schema
- [x] Default mode 변경
- [x] 커밋 14b0295

### Phase 2: 3-term Supermolecular Decomposition
- [x] `qm_int_kcal_frozen` 필드
- [x] `interaction_kcal` 보존
- [x] 커밋 14b0295

### Phase 3: Charge Correction Policy (🟢 Stage 1+2 code complete, baseline recompute 대기)
- [x] R-15/R-16 guard + audit module (Stage 1 완료, 커밋 4034ee6)
- [x] Option B3 partitioning 재설계 (Stage 2, 2026-04-19 코드 완료)
- [x] schema 0.6.4 + numbering_convention + auto rationale (Stage 2)
- [ ] baseline recompute on w4_1_s2 (Stage 3, user 승인 후)
- [ ] 기존 JSON audit 완료 + `_archive/v0.6_wrong_charge/` 격리 (Stage 3)

### Phase 4: Pre-QM Budget Enforcement + Charge Parity Filter (신규)
- [ ] `utils/pre_qm_filter.py` 모듈
- [ ] `estimate_n_qm()` 함수
- [ ] Charge-aware pre-filter (R-15/R-16 재사용)
- [ ] `rejected_designs.json` output
- [ ] MD 진입 전 filtering

### Phase 5: MM-GBSA + DFT Hybrid (신규, Phase 3+4 precondition)
- [ ] `utils/cluster_snapshots.py` 모듈
- [ ] Representative selection (k-medoids RMSD)
- [ ] Weighted average ΔG in `rank_results.py`
- [ ] SciVal 승인 필요

### Phase 6: B10 Level Shift (대기)
- [ ] w4_1_s2/snap01 결과 기반 결정
- [ ] SciVal 7 조건 모두 구현 (condition)

### Bundle B1-B9 Status (2026-04-19)
- [x] B2 grids_level=2 (env var override)
- [x] B3 conv_tol=1e-7 (env var override)
- [x] B5 direct_scf_tol=1e-12 (env var override)
- [x] B6 init_guess=minao (huckel opt-in)
- [x] B8 min_ao_blksize=64 (gpu4pyscf import 전)
- [x] B9 cupy mempool 0.65 (env var override)
- [x] DF auto-mode + threshold env var
- [-] B1 deprecated (rks_lowmem 비호환)
- [ ] B7 conditional (T-Rank-1 smoke test 통과 후)
- [ ] B10 pending (현재 run 결과 대기)

## v0.7 — WT Counterfactual + Calibration

### 진입 조건 (v0.6 에서 달성)
- [ ] 3+ design × 5 snap converged
- [ ] Path 진단 🟢
- [ ] Pre-QM filter 적용
- [ ] R-15/R-16 guard: 모든 snap `charge_consistency_audit: passed`

### Phase 7: WT Counterfactual
- [ ] `variants/{ncaa,wt}/` 디렉토리 구조
- [ ] `assess_wt_comparability()` 3-tier 판정
- [ ] 1-Me-Trp registry addition

### Phase 8: Compstatin Pilot
- [ ] PDB 2QKI target preparation
- [ ] 3 variants: Val4, Trp4, 1-MeW4
- [ ] Calibration fit (2 points)

### Phase 9: Sa-D3 Pilot
- [ ] PDB 7TL8 target preparation
- [ ] 2 variants: N-Me, ΔN-Me
- [ ] Calibration 확장 (3 points)

## v0.8 — Calibration Expansion + Publication Track
- [ ] 15-IgBP pilot
- [ ] MD replicates (3x)
- [ ] Bootstrap CI
- [ ] Paper draft

## Completed ✅

- [x] v0.5 SSOT refactor (2026-04-17)
- [x] v0.5 Bug 1-3 fix (2026-04-17)
- [x] Path agent 도입 (2026-04-19 Rev. 1)
- [x] v0.6.0 QM/MM physical validity (2026-04-18 14b0295)
- [x] SCF bundle B2/B3/B5/B6/B8/B9 (2026-04-19 uncommitted)
- [x] DF auto-mode (2026-04-19 uncommitted)
- [x] Monitoring tools regex fix (2026-04-19 uncommitted)
- [x] CLAUDE.md Rev. 2 반영 (2026-04-19)
- [x] R-15/R-16 Charge guards + audit module Stage 1 (2026-04-19)
```

---

## 🚀 Git Commit & Push Policy (Rev. 2 유지)

### User-Controlled Push

**Push 는 사용자만 수행**. Agent 는 커맨드 제시 후 허가 대기.

### 커밋 구조 제안 (Rev. 2)

오늘 세션의 uncommitted 변경사항을 **구조적으로 분할** 권고:

```bash
# Commit 1: 모니터링 도구 버그 수정 (독립)
git add utils/qmmm_proactive_oom.py utils/cycle_benchmark.py
git commit -m "fix(monitoring): regex patterns for new log prefix format

- qmmm_proactive_oom.py v0.3.3 → v0.3.5
  - 6 regex patches for [tag/snap_stem] format
  - New SNAP column + wall-time estimate
  - Flicker fix (auto_refresh=False + screen=True)
- cycle_benchmark.py: 2 regex patches

Refs: RECENTv06.md Section 1.11"

# Commit 2: SCF 최적화 번들 (B2/B3/B5/B6/B8/B9)
git add utils/run_qmmm.py UPDATE.md
git commit -m "feat(qmmm): SCF optimization bundle B1-B9

Ranking-invariant optimizations with env var overrides:
- B2: grids_level=2 (ORCA DefGrid1)
- B3: conv_tol=1e-7
- B5: direct_scf_tol=1e-12
- B6: init_guess=minao (huckel opt-in, w4_1_s2 divergence fix)
- B8: min_ao_blksize=64 (gpu4pyscf import 전)
- B9: cupy mempool 0.65

Deprecated: B1 (rks_lowmem, QMMMSCF 비호환)
Conditional: B7 (mo_coeff warm-start, T-Rank-1 smoke test 대기)
Pending: B10 (level shift)

SciVal memory: bundle_b1_b9_verdict.md
Refs: RECENTv06.md Section 1.9, 1.10"

# Commit 3: DF auto-mode + VRAM cascade
git add utils/run_qmmm.py
git commit -m "feat(qmmm): DF auto-mode based on VRAM detection

- DF auto on if VRAM >= 20 GiB (threshold env var)
- RTX 5070 Ti 16GB → auto off → direct-SCF
- RTX 5090 32GB → auto on (future-ready)
- --df-mode {auto,on,off}, --no-df alias

VRAM cascade (diagnostic only, no perf gain):
- Tier 0/1/1.5/2: GPU/UM-VRAM/UM-RAM/CPU
- gpu4pyscf 내부 raw cudaMalloc 우회 확인

Refs: RECENTv06.md Section 1.2, 1.4"

# Commit 4: CLAUDE.md Rev. 2
git add CLAUDE.md UPGRADE.md
git commit -m "docs: CLAUDE.md Rev. 2 — 2026-04-19 session findings

- Confirmed hardware limits (16GB VRAM, DF 불가)
- Rev. 2 Phase 재정렬: Pre-QM filter + Hybrid MM-GBSA
- Deprecated list (DF-JK, rks_lowmem, budget 300, huckel default)
- Ranking-invariance regime 공식화
- v0.7 calibration 진입 조건 엄격화 (3 design × 5 snap converged)

Refs: RECENTv06.md (full session log)"

# Commit 5: RECENTv06.md 자체 (optional, archive)
git add RECENTv06.md
git commit -m "docs: archive RECENTv06.md — 2026-04-19 session log"
```

**Push 는 user 가 수동 실행**:
```bash
# user 가 직접:
git push origin <branch>
```

---

## 📞 Agent Invocation 예시 (Rev. 2)

### Archi 호출
```
[Archi 에이전트 호출]

컨텍스트: CLAUDE.md Rev. 2 반영 완료.
RECENTv06.md 참조하여 v0.6 Phase 4 (Pre-QM filter) 착수.

Task:
1. 현재 실행 중인 w4_1_s2/snap01 결과 확인
2. 탈락 design (w4_0_s2, w1_0_s1) _archive/ 이동 계획
3. Phase 4 implementation plan (pre_qm_filter.py 스케치)
```

### Path 호출
```
[Path 에이전트 호출]

컨텍스트: Runner 가 w4_1_s2/snap01 완주 (cycle 76+ plateau,
HOMO-LUMO near-degeneracy 의심).

Task:
1. 최종 JSON 의 converged, energy_total_hartree, qm_int_kcal_frozen 확인
2. 물리적 타당성 진단 (metastable 여부)
3. B10 level_shift 필요 여부 권고
4. 지금까지 완주 snap 들의 v0.7 calibration 진입 조건 평가
```

### SciVal 호출
```
[SciVal 에이전트 호출]

컨텍스트: Phase 4 pre-QM filter 구현 위한 사전 검증.

Task:
1. `estimate_n_qm()` algorithm (topology + contact_residues count) 검증
2. Budget 500 enforcement 가 design space 에 주는 영향 평가
3. 탈락 design 의 재사용 가능성 (예: 다른 target 에 사용 가능?) 판단
```

---

## 💾 Memory Hygiene (Rev. 2 updated)

### Path 의 memory 카테고리

- **pathology**: Recurring symptom patterns (신규)
- **feedback**: Diagnostic rule of thumb
- **project**: Current pipeline state
- **reference**: External benchmarks, literature

### Path 가 기록할 Rev. 2 패턴 예시

1. "HOMO-LUMO near-degeneracy: w4_1_s2 계열. Init guess 민감. B10 level_shift 필요"
2. "n_qm > 500: hard_fail, pre-QM filter 에서 탈락"
3. "charge correction silent: topology mode 후 발동 여부 일관. Phase 3 개입 불필요"
4. "with_df.max_memory: gpu4pyscf 가 무시. 레거시 속성"
5. "def2-SVP + weigend aux basis: PySCF alias 로 def2-universal-jfit 과 동일 naux"

이런 패턴을 축적하면 다음 세션에서 Path 가 즉시 진단 가능.

---

## ⏹ Stop Conditions (Rev. 2 강화)

Agent 루프가 다음 도달 시 즉시 user 에게 에스컬레이션:

1. Path 가 동일 증상 3회 연속 🔴
2. SciVal 이 동일 제안 2회 이상 거부
3. Runner 에서 command-line 에러 (dependency, import, permission)
4. 디스크 공간 < 10 GB
5. `target_card.json` 생성 불가 (수동 curation 필요)
6. Hardware limit (VRAM OOM + CPU fallback 도 실패)
7. **SCF 수렴 실패 design 이 50% 초과** (신규, Rev. 2)
8. **B10 level_shift 후에도 metastable 지속** (신규, Rev. 2)
9. **MD 진입 전 pre-filter 에서 design pool 80% 이상 탈락** (신규, Rev. 2)

---

## 📊 Scope 재조정 — v0.6 → v0.7 → v0.8 (Rev. 2)

### v0.6 (현재, ~2026-04-30 목표)
- Physical validity recovery
- Pre-QM filter
- Hybrid MM-GBSA + DFT pipeline 스케치
- 3-5 design × 5 snap converged 확보

### v0.7 (~2026-05-15 목표)
- WT counterfactual mode
- Compstatin pilot (2 variants: Trp4, 1-MeW4)
- Sa-D3 pilot (2 variants: N-Me, ΔN-Me)
- Post-hoc calibration layer

### v0.8 (~2026-06-15 목표)
- MD replicates (3x)
- 15-IgBP pilot
- Statistical rigor
- Paper draft

### v0.9+ (미래)
- Hardware 업그레이드 후 DF-JK 재활성화
- n_qm 1000+ 가능
- Active ncAA design + synthesis recommendation

---

**마지막 줄 (Rev. 2)**: 이 문서가 진실의 원천이다. Agent 간 충돌 시 이 문서가 우선한다. Hardware 물리 경계를 인정하되, 영리한 design 으로 최대한 활용한다. v0.7 calibration 진입은 v0.6 Phase 4-5 완료 후에만.
