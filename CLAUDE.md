# UPDD v0.6 Refactor — Physical/Chemical Validity Recovery (Rev. 3)

**작성일:** 2026-04-19
**개정:** Rev. 3 (2026-04-19 저녁 세션 — Phase 3 Stage 1 merged, JoltQC eval side-quest 완료)
**이전 Rev:** Rev. 2 (2026-04-19 낮 세션)
**이전 버전 archive:** `archive/CLAUDE_Rev2_20260418.md`
**이전 배포:** v0.5.x (SSOT refactor — SUCCESS/MARGINAL tier classification)
**목표 배포:** v0.6.x → v0.7 Calibration 진입 조건 충족

---

## 🚨 세션 시작 최우선 Action (Rev. 3.1)

**어느 세션이든 첫 턴에서 Archi 는 다음 문서들을 순서대로 읽어야 task 설계를 시작할 수 있다** (자세한 boot sequence + anti-pattern 은 아래 `📋 First Task for Archi § 세션 시작 최초 Step` 참조):

1. `CLAUDE.md` (본 파일) → 2. `UPGRADE.md` → 3. `UPDATE.md` → 4. `UPDDnow.md` → 5. `UPDDfuture.md`

> User 가 "간단한 수정" / "바로 X 하자" 라고 지시해도 Archi 는 boot sequence 완료 후에만 응답. Skip 은 R-10 위반.

---

## 🐍 Environment Setup (Rev. 3 clarification + 2026-05-04 pytest hardening)

### UPDD main project (this document)
- Python env: **conda env `qmmm`** (Python 3.10 + openmm 8.4 + pytest 9.0.3)
- Canonical interpreter: `/home/san/miniconda3/envs/qmmm/bin/python`
- Command pattern: `python utils/run_qmmm.py ...` (with the qmmm env activated)
- Location: `/home/san/UPDD_proj/`

### Pytest invocation (REQUIRED — do NOT use system Python)

```bash
/home/san/miniconda3/envs/qmmm/bin/python -m pytest tests/
```

- **Why**: tests import `openmm`, `mdtraj`, `numpy` from the qmmm env. Running `pytest` via system Python or `/home/san/miniconda3/bin/python` (root miniconda, no openmm) silently fails 13+ tests with `ModuleNotFoundError`.
- **Reference**: `outputs/analysis/test_failure_diagnosis_20260504.md` — 2026-05-04 ULTRAREVIEW that disambiguated the 27-failure pytest pattern into two orthogonal infrastructure issues (wrong interpreter + sys.path pollution), neither of which was a PR-21 regression.

> Test isolation note (2026-05-04): `tests/test_verify_pbc_integrity.py` uses `importlib.util.spec_from_file_location` (not `sys.path.insert`) to avoid pytest-session pollution that previously caused `arrow_menu` ImportError cascade in 14 tests (scripts/updd_cli.py was shadowing utils/updd_cli.py). New tests that load files from `scripts/` should follow the same idiom — see `tests/test_updd_cli.py:38-46` for the reference pattern.

### Related external projects
- JoltQC eval: pixi-managed, `pixi run python ...` (독립 프로젝트)
- Cross-project reference only, document-only interface


## 📋 Rev. 3 Changes Summary

Rev. 2 → Rev. 3 delta (빠른 diff):

### 🔴 Critical Progress

- **Phase 3 Stage 1**: ✅ **merged to origin/main (4034ee6)** — R-15/R-16 guards + audit module 완성, 23/23 tests pass
- **Phase 3 Stage 2**: 🔄 **진행 중** — partitioning redesign (`sidechain_from_cb` → `whole`, budget 500 → 600), schema v0.6.4
- **w4_1_s2/snap01 "v0.6 first valid" status**: **WITHDRAWN** (electronically wrong state, 1 electron deficit, metastable SCF convergence)
  - Archive: `outputs/_archive/v0.6_wrong_charge/` with manifest
  - Chkfile archive: `outputs/_archive/chkfile_benchmarks/` (SCF profiling reference only)

### 🟢 New Structural Additions

- **🚀 QM/MM Acceleration Track** (version-agnostic, 별도 track) — JoltQC integration plan
  - Axis 1 (Numerical): 🟢 DEFINITIVE PASS
  - Axis 2 (Speed): 2-4× speedup confirmed
  - Axis 3 (Physical accuracy): ⏸ DEFERRED pending Phase 3 완료
  - Integration mode: Opt-in first (environment variable)

### 📜 Session Incidents (Phase 3 description 내 inline)

- **2026-04-19 10:32**: tmux crash (Whoopsie ID b45edeed), w4_1_s2/snap02 SCF 중단
- **2026-04-19 15:00**: Parity contradiction discovery (JoltQC eval 중), Phase 3 재개 trigger

### 🆕 Blockers Introduced

- snap02 clean restart: Phase 3 Stage 1+2 merge 까지 block
- 다른 design 확대 (w1_0_s1, w4_0_s2): R-15/R-16 guard + Stage 2 partitioning merge 까지 block
- JoltQC production integration: Phase 3 완료 + external Stage 7 re-run 까지 block

### 📁 Cross-Project Reference (외부, user-managed, document-only)

- `/home/san/joltqc_updd_eval/` — JoltQC evaluation side-quest (independent project)

### 🏷 GitHub Issues Filed

- #XX (P0-critical) Parity-incompatible QM/MM partitioning — main issue
- #YY (P2) target_card whole_residue_exceptions numbering offset — sub-issue
- #ZZ (P4) run_qmmm.py L1579 misleading comment — sub-issue

### 📝 References Added

- SciVal verdict chain: `.claude/agent-memory/scival/verdict_target_card_charge_rationale_20260419*.md` (v1, v2, v3 addendum)

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
5. **Budget 500 atom 하드 제약** — Rev. 3 에서 **600 으로 완화** (Phase 3 Stage 2, B3 partitioning redesign)
   - 6WGN DxxGQ motif, α3-groove contact residue 는 Zhang 2020 Cell Chem Biol 명시
   - `enforcement="hard_fail"` 하드코딩 → n_qm > 600 design 은 partition_error

### v0.6 의 재정의된 목표

Rev. 1 에서 "15× 가속" 을 기대했으나 **실측에서 거의 불가능함이 확인**. Rev. 2 현실적 목표 유지 + **Rev. 3 에서 JoltQC track 추가**:

- ❌ NOT: DFT QM/MM 를 대규모 가속 (단일 path)
- ✅ YES: 물리/화학적 correctness 복구 (charge state consistency, topology integrity) — **R-15/R-16 Stage 1 ✅ merged**
- ✅ YES: Pre-QM filtering 으로 **실패할 design 을 사전 탈락**
- ✅ YES: MM-GBSA 와 DFT QM/MM 의 **hybrid 활용** 으로 정보량 확보
- ✅ YES: v0.7 Calibration 진입 조건 충족 (최소 3 design × 5 snap converged)
- ✅ **NEW (Rev. 3)**: JoltQC 2-4× 가속을 **opt-in** 으로 제공 (Axis 3 검증 후 default 전환 검토)

### 하드웨어 업그레이드 경로 명시

현재 설계는 **RTX 5070 Ti 16GB 고정** 전제. 다음 업그레이드 시 자연스럽게 확장:
- **RTX 5090 32GB**: DF auto-mode 가 자동으로 on → 2-3× 가속
- **H100 80GB**: n_qm 1000+ 가능, replicate 확장

업그레이드 전까지는 **현 하드웨어의 물리적 한계 수용** + **영리한 pipeline 설계** + **JoltQC opt-in acceleration** 으로 진행.

---

## 🎯 UPDD v0.6 Constitution

### Core Principles (Rev. 2 유지)

[Rev. 2 Core Principles 전체 유지 — 변경 없음]

### Absolute Rules

- **R-1 ~ R-3** (Coder 기본): 변수 순서, 기계적 치환 금지, import 무결성
- **R-4**: SciVal ❌ 판정 시 코드 배포 금지
- **R-5**: Runner 검증 없이 완료 선언 금지
- **R-6**: Path 재진단 전 task 종료 금지
- **R-7**: 탈락 데이터/결과는 `_archive/` 로 이동 (삭제 금지)
- **R-8**: 모든 수정은 `UPDATE.md` 상단에 기록 (entry format: Agent/Category/Regime/변경/이유/검증/impact/관련)
- **R-9**: Git 커밋/푸시 전 사용자 허가 필수 (커맨드는 제시, 실행 금지)
- **R-10** (Rev. 3.1 mechanism 명시화): **UPGRADE.md 체크박스 연동 의무**. 작업 시작 전 Archi 는 UPGRADE.md 조회 필수 (precondition + 중복 방지); 작업 완료 후 UPDATE.md 엔트리와 **동일 session 내** UPGRADE.md 체크박스 `[x]` 로 전환 + commit hash + 날짜 기록. Anonymous check (hash/date 생략) 금지. 자세한 루틴은 §"🔄 UPDATE.md / UPGRADE.md 관리" 참조.
- **R-11** (Rev. 2): **Ranking regime 외 주장 금지**. Absolute ΔG 값의 정량 의미는 v0.7 calibration 완료 전까지 주장 불가. 결과 리포트에 "(ranking-only)" 명시.
- **R-12** (Rev. 2): **설정 속성 효과 소스 검증**. `mf.with_df.max_memory` 같은 "읽히지 않는 속성" 존재. 새 환경변수/config 추가 시 upstream grep 으로 참조 확인 필수.
- **R-13**: 예약 (향후 Long SCF nohup rule 자리).
- **R-14**: 예약.
- **R-15** (Rev. 2, 2026-04-19): **Charge Declaration SSOT Enforcement (magnitude + parity + rationale)**.
  `target_card.json` 의 `target_iso_net_charge` 및 `binder_net_charge` 는 QM 실행 시 snapshot PDB 로부터 chemistry-true 값과 **magnitude-level 일치 필수**. 1-bit parity check 외에 추가 magnitude check + rationale auto-regenerate. Mismatch 시 `ChargeDeclarationMismatch` 로 fail-fast; JSON diag 에 chemistry-computed charge + declared charge 양쪽 기록. `target_iso_charge_rationale` 텍스트는 code-generated only (수동 편집 금지). 기존 `OddElectronError` (parity-only) 는 magnitude 와 독립 safety net 으로 유지.
  **Rev. 3 구현 상태**: ✅ Stage 1 merged (4034ee6), 23/23 tests pass, historical audit Stage 3 대기.
- **R-16** (Rev. 2, 2026-04-19): **Binder Chemistry SSOT Enforcement**.
  Binder net charge 는 snapshot PDB 로부터 residue identity + atom name pattern (ARG HH11-HH22, ASP OD-HD, GLU OE-HE, LYS HZ, HIS HD1/HE2 tautomer) 을 기반으로 runtime 계산한다. `target_card.binder_net_charge` 는 "expected hint" 필드일 뿐이며 runtime PDB 가 ground truth. Cyclic 여부는 LEU1 (또는 첫 잔기) N-H count 와 마지막 잔기 OXT 유무로 자동 감지. 불일치 시 `BinderChargeMismatch` 로 fail-fast.
  **Rev. 3 구현 상태**: ✅ Stage 1 merged (4034ee6), `utils/charge_topology.py` + `utils/audit_charge_consistency.py` 완성.
- **R-17** (Rev. 4, 2026-04-19, refined per SciVal 10th verdict §6.3): **Cofactor Preservation SSOT Enforcement**.
  `target_card.cofactor_residues` declares non-residue non-solvent entities whose **identity, topology, formal charges, and force-field parameters MUST be preserved** across the pipeline; their **coordinates and inter-atomic distances MUST be allowed to evolve** under MD thermal sampling. The pipeline must never substitute a post-hoc rigid body, crystal-frozen position, or time-averaged centroid for MD-sampled snapshot coords. Enforced at 6 pipeline gates (G1-G6: preprocess, AF2 reinsertion, parameterize, MD, snapshot extract, QM/MM build). Target-agnostic — driven by target_card declaration, not hardcoded resnames.
  **FORBIDDEN** (SciVal 10th §6.3): (a) post-hoc grafting of crystal cofactor into MD snapshot without full re-equilibration; (b) rigid-body superposition of crystal-frame cofactor onto snapshot frame (KABSCH-based); (c) time-averaged cofactor coords in place of per-snapshot MD coords.
  **Exceptions** (`diagnostic_mode` only): flagged runs may bypass coord preservation for probe purposes but must carry `diagnostic_mode: true` + `not_for_ranking: true` + exclusion from v0.7+ calibration.
  **Exception classes**: `CofactorMissingError` (declared but absent), `CofactorChargeSumMismatch` (sum-of-charges drift), `CofactorParamMissingError` (FF param lookup fail), `CofactorCoordFrameError` (coord-frame mismatch — centroid distance check, RMSD drift, etc.).
  **MM injection**: `treatment: "mm_point_charge"` cofactors MUST be loaded into QM/MM background charge set from **MD-sampled snapshot positions**, not crystal coords. Sum must match declared total.
  **구현 상태**: Stage 4 (R-17 full fix, schema v0.6.5 + 6-gate enforcement + new cofactor_errors module + 15 tests) ✅ **implementation complete 2026-04-19** (66/66 tests pass). Orchestrator wiring into UPDD.py in progress. Empirical validation via full pipeline rerun from Step 4 pending. See `docs/known_issues.md`.

---

## 🤖 6-Agent 구조 (Keeper 신설, Rev. 4 에서 infrastructure 확장)

| Agent | 역할 | 권한 | Veto 클래스 | 호출 타이밍 |
|---|---|---|---|---|
| **Archi** | 아키텍트, 전략, 조율 | Read-only, 지시 권한 | — (조율 권한) | 전체 계획, 작업 간 조율 |
| **SciVal** | 사전 과학적 검증 | Read-only, Veto | **Science (R-4)** | 코드 변경 **전** 검증 |
| **Coder** | 코드 구현, 수정 | Write | — | Archi 지시 + SciVal 승인 후 |
| **Keeper** | Cross-stage integrity 게이트 | Read-only, Advisory Gate | **Integrity (R-17)** | Coder 완료 후 / Runner 실행 **전** |
| **Runner** | 실행, 로그 수집, 기본 검증 | Bash, Read | **Execution (R-5)** | Keeper 🟢 후 |
| **Path** | 결과의 물리/화학/병리학적 진단 | Read-only + Bash (read), Advisory | **Pathology (R-6)** | Runner 완료 후 |

### Agent Veto Taxonomy (Rev. 4)

4개 독립 veto 축: 각 실패 유형은 정확히 하나의 책임자에게 배정됨.

| Veto Class | Owner | Rule | Scope |
|---|---|---|---|
| **Science** | SciVal | R-4 | 과학적 오류 (물리/화학적 타당성, 문헌 일관성) |
| **Execution** | Runner | R-5 | 실행 실패 (output 생성 실패, 검증 실패) |
| **Pathology** | Path | R-6 | converged 결과의 물리/화학적 이상 (metastable, wrong-state) |
| **Integrity** | Keeper | R-17 | 선언-실제 drift (target_card 선언이 artifact 에 반영 안 됨) |

실패를 마주치면 "어떤 class 인가?" 를 먼저 판별 → 해당 agent 에게 delegate. 4 class 는 서로 **orthogonal** (한 이슈가 복수 class 에 걸치는 경우 드물며, 걸치면 primary class 먼저 처리).

### Rev. 4 에서 신설된 Keeper 의 역할

2026-04-19 세션 Stage 3 실패 (GNP/Mg cofactor 가 target_card 에 선언됐으나 MD snapshot 에 부재) 가 SciVal/Path 영역이 아닌 **cross-stage integrity** 문제임이 드러남. Keeper 가 이 gap 을 담당:

- **Pre-flight**: 컴퓨트 launch 전 input integrity 검증 (declared cofactors present, schema version compatible, partition dispatch correct, numbering offset consistent)
- **Post-generation**: 각 pipeline stage 의 artifact 가 downstream consumer 기대치와 일치하는지 검증
- **Schema migration**: 스키마 버전 bump 시 모든 consumer 가 함께 업데이트됐는지 검증
- **Delegation first**: 기존 in-code guard 가 있는 checks (R-15/R-16 runtime guards, `partition_qmmm` budget enforcer) 는 reimplement 하지 않고 delegate. 자체 구현은 guard 없는 checks 만 (numbering offset, env var usage, schema consumer audit).

**Invocation 예시**:
```
[Keeper 에이전트 호출 — pre-SCF gate]

컨텍스트: w4_1_s2/snap01 full SCF launch 직전.

Task:
1. COFACTOR-01 gate: target_card.cofactor_residues 가 snapshot PDB MM charge set 에 주입됐는지
2. CHARGE-01 gate: R-15/R-16 runtime guard 가 예상 magnitude 와 일치할 것인지 (dry-run)
3. NUMBERING-01 gate: MD↔crystal 매핑이 모든 contact residue 에 일관
4. VRAM-01 gate: free VRAM ≥ 예상 peak

Expected: ALL 🟢 → SCF launch 승인. ANY 🔴 → block + escalate.
```

### Invocation ordering (Rev. 4 — Phase 1)

```
Archi → SciVal → Coder → Keeper (pre-execution) → Runner → Path
```

Phase 2 (v0.7+): post-Runner Keeper gate 추가 (artifact integrity for downstream consumer).
Phase 3 (v0.8+): pre-Coder Keeper 선택적 (plan's precondition validation; SciVal pre-check 와 redundant 하지 않게).

Current (Phase 1) = pre-Runner Keeper 단일 point.

### Rev. 2 에서 강화된 Path 의 역할

Rev. 1 에서 Path 는 "결과 진단" 역할이었다. 2026-04-19 세션에서 드러난 실패 패턴이 **다양하고 복잡**하므로 Rev. 2 에서 Path 의 진단 범위가 확장된다:

**Path 가 진단할 실패 카테고리 (Rev. 2 + Rev. 3 추가)**:

1. **수렴 실패** (max_cycle 초과)
2. **Metastable 수렴** (converged=True 지만 unphysical energy)
3. **HOMO-LUMO near-degeneracy** (w4_1_s2 계열)
4. **Budget violation** (n_qm > 600, Rev. 3 update)
5. **QM region inconsistency** (topology mode 전 snap 비교)
6. **Charge correction 불일치** (total_protons parity 변화)
7. **SCF init_guess 발산** (huckel 에서 ΔE > 1000 Ha)
8. **Library 레벨 silent failure** (with_df.max_memory 류 무시 속성)
9. **Charge magnitude mismatch** (R-15/R-16 guard trigger) — **Rev. 3 신규**
10. **Wrong-state convergence** (parity 통과 but magnitude 틀린 SCF) — **Rev. 3 신규**

각 카테고리마다 Path 의 **증상 → 원인 → 권고** template 적용.

### Rev. 3 cross-project 명시

Agent 의 메모리 범위는 **project scope 로 한정**. Cross-project dependency 는 **document-only interface** 로 처리:

- UPDD 의 agents (Archi/SciVal/Coder/Runner/Path) 는 UPDD_proj 범위 내에서만 작동
- 외부 프로젝트 (예: JoltQC eval) 는 사용자가 직접 관리, UPDD agents 는 인지하지 않음
- 외부 프로젝트 결과는 `docs/known_issues.md` 또는 `.claude/agent-memory/scival/` cross-reference 로만 공유

---

## 🎯 v0.6 작업 우선순위 (Rev. 3, Stage 1 merged 반영)

### 원칙

v0.6 의 핵심은 **correctness recovery** + **pre-filtering**. **Performance optimization 은 v0.8 이후 (hardware 업그레이드 + JoltQC integration)**.

### Phase 1: ✅ **완료** — Topology-based QM Region Selection

[Rev. 2 내용 유지]

### Phase 2: ✅ **완료** — 3-term Supermolecular Decomposition

[Rev. 2 내용 유지]

### Phase 3: 🔄 **진행 중** — Pre-SCF Charge Guard + R-15/R-16 Enforcement + Audit Module

**현재 상태 (Rev. 3, 2026-04-19 저녁)**:

- **Stage 1**: ✅ **merged to origin/main (4034ee6)**
  - `compute_qm_net_charge_topology()` + `compute_binder_chem_charge()` 헬퍼 완성
  - `utils/charge_topology.py` 신규 모듈
  - `run_qmmm.py::build_qm_mol()` R-15 guard 삽입, R-16 binder SSOT 검증
  - `utils/audit_charge_consistency.py` 신규 (dry-run default)
  - JSON schema 확장 (`qm_net_charge_declared/computed`, `binder_charge_declared/computed`, `charge_consistency_audit`)
  - 23/23 tests pass
  - Git history scrub 완료 (filter-repo, commit SHA 4034ee6 로 변경)

- **Stage 2**: 🔄 **진행 중**
  - Scope: Option B3 partitioning 재설계 (target_residues_mode=whole, max_n_qm=600)
  - Schema v0.6.4: `numbering_convention` 필드 추가 (md_to_crystal_offset=3 for 6WGN), `whole_residue_exceptions` 제거
  - Auto-generated rationale (수동 curation 금지, R-15 enforcement)
  - Smoke test criteria: n_qm ~580, link atoms 17→5-7
  - Stage 2 completion 시 issue #XX sub-issue #YY (numbering offset) close

- **Stage 3**: ⏸ **v0.7 진입 조건, Stage 1+2 완료 후 착수**
  - 기존 12개 topology-mode JSON audit 완료 (PASS/PARITY_ONLY/FAIL 분류)
  - v0.5 archive 의 historical audit 포함 (JoltQC Stage 7 baseline integrity)
  - Mismatch 결과는 `outputs/_archive/v0.6_wrong_charge/` 또는 `_archive/v0.5_wrong_charge/` 로 격리

### Phase 3 재개 배경 (2026-04-19 Rev. 3 — Session Incidents inline)

**Trigger**: JoltQC evaluation side-quest (외부 프로젝트) 에서 target_card 검증 중 발견.

**발견 경위**:
1. **10:32 — tmux crash** (Whoopsie ID b45edeed-3b8f-11f1-abe6-fa163efa23ad)
   - w4_1_s2/snap02 SCF 중단
   - /tmp/ 에 3 chkfiles orphaned (140MB 각, total 420MB)
   - R-13 예약 자리의 rationale: 이 incident 로 확정 (Long SCF nohup enforcement)

2. **10:32 — 14:30 — 회복 + warm-start 유혹**
   - snap01/02 chkfiles 발견 → "warm-start 로 시간 절약" 유혹
   - SciVal 1차: HDF5 energy discontinuity 확인 (snap01 E_tot = -3233 Ha, snap02 partial = -10559 Ha, Δ ≈ 7326 Ha)
   - **판정**: warm-start ❌ (non-stationary WF, near-degenerate system metastable 위험)
   - **결정**: clean minao restart 수용

3. **14:30 — JoltQC eval side-quest 착수** (clean restart 대기 시간 활용)
   - Stage 0-7 모두 PASS
   - target_card 검증 과정에서 HIS tautomer 체크

4. **15:00 — Parity contradiction 발견** (CRITICAL)
   - Finding 1: topology mode residue list design-invariant (정상)
   - Finding 2: HIS protonation은 tautomer flip 위험 (Bayesian inference: snap01 converged ⟹ HIS92=HID, 실제 PDB atom pattern 으로 확증)
   - Finding 3: target_card rationale literally wrong (ASP9+ASP92+GLU62 — 실제 crystal 은 VAL9/HIS92/SER62)
   - Finding 4: even-parity silent pass (R-12 의 odd_electron_error 는 parity 만 감지, magnitude mismatch silent)
   - Finding 5: Phase 3 재개 trigger (CLAUDE.md Rev. 2 의 "snapshot 간 charge 불일치 감지" 조건 성립)
   - **Finding 6 (CRITICAL)**: **Parity contradiction**
     - Binder Z=1060 (even), Target+link Z=583 (odd), Total Z=1643 (odd)
     - Chemistry-true charge: binder(-1) + target(-3) = -4 → 1647 electrons (ODD) → RKS 불가능
     - Declared charge: -3 → 1646 electrons (EVEN) → RKS 통과 → **silent wrong-state convergence**
     - snap01 `E_total = -10559.88 Ha, converged=True` 는 **electronically wrong by 1 electron**

5. **15:30 — SciVal verdict chain**
   - v1: Phase 3 (Charge Correction Policy) 공식 재개, R-15 제안
   - v2: Option B3 (whole-residue + budget 600) 승인, A/B1/B2/C 기각
   - v3 addendum: Historical disposition — Option i (complete discard) primary, rescue options (Koopmans, thermocycle) 기각 (DFT SIE wrong-charge amplification)

**snap01 결과 처리**:
- JSON → `outputs/_archive/v0.6_wrong_charge/` (manifest 포함)
- Chkfile → `outputs/_archive/chkfile_benchmarks/` (SCF profiling reference only, energy_invalid=true)
- "v0.6 first valid result" 지위 **WITHDRAWN**
- v0.7 진입 counter 에서 제외

**재평가 조건 (Stage 1 완료 판정)**: ✅ 충족
- 기존 JSON audit 완료 (Stage 1 merged 의 audit module)
- R-15/R-16 guard 동작 확인 (regression test 에서 Mismatch 감지)
- CLAUDE.md Phase 3 Stage 1 🟢 완료 표시

**Stage 2 완료 판정 조건**:
- target_card v0.6.4 schema applied (numbering_convention, whole-mode partitioning)
- Partition smoke test pass (n_qm ~580, link atoms 5-7)
- Issue #XX Stage 2 comment + sub-issue #YY close

### Phase 4 (확장): 🔴 **진행 중** — Pre-QM Budget Enforcement + Charge Parity Filter

**원 문제**: w4_0_s2 (n_qm=578), w1_0_s1 (n_qm=539-561) 가 budget 500 초과 → partition_error 로 **QM/MM 진입 전 무조건 실패**.

**Rev. 2 확장 (2026-04-19)**: Phase 3 발견으로 scope 확장됨. 기존 n_qm filter + **charge parity + chemistry-declared 일치 filter** 추가. AF2 output + target_card 로 `compute_binder_chem_charge()` (R-16) 및 `compute_qm_net_charge_topology()` (R-15) 을 돌려 parity-incompatible / magnitude-mismatch design 을 MD 진입 전 탈락.

**Rev. 3 update**: **Budget 500 → 600 완화** (Phase 3 Stage 2 결정). w4_0_s2 는 새 budget 으로 재평가 가능 (578 < 600). w1_0_s1 도 재평가 가능 (539-561 < 600).

**현재 정책**: `enforcement="hard_fail"` (SciVal 승인), threshold=600 (Rev. 3)

**v0.6 에서 할 일**:

1. **탈락 design 조기 감지** — MD stage 이전에 AF2 output 에서 `estimate_n_qm()` 사전 계산 (threshold=600)
2. **Charge-aware pre-filter** — 같은 n_qm 이어도 binder identity 에 따라 parity 다름 → 각 design 의 binder residue 분석 → declared vs computed charge 자동 일치 검사 → 불일치 design 은 pre-filter 탈락
3. **Pipeline 레벨 필터** — n_qm > 600 design 은 MD 에서 제외 (MD time 낭비 방지)
4. **Logging 강화** — 탈락 design 목록을 `outputs/*/rejected_designs.json` 에 기록
5. **탈락 이유 분류**:
   - `n_qm_exceeded` (hardware limit, >600)
   - `chemistry_parity_incompatible` (R-15/R-16 neg-mag guard) — Rev. 2
   - `contact_residue_count` (과학적 이슈)
   - `topology_error` (partition 실패)

**수용 기준 (Acceptance Criteria)**:
- MD 진입 전 탈락 design 확인됨
- 탈락 이유가 `rejected_designs.json` 에 기록됨
- Runner 가 "N designs passed budget+charge, M designs rejected" 리포트
- Phase 3 + Phase 4 공동 merged 이전 Phase 5 BLOCK

**구현 위치**: `utils/pre_qm_filter.py` (새 모듈)

### Phase 5 (New): 🟡 **계획** — MM-GBSA + DFT Hybrid Workflow (Phase 3+4 precondition)

[Rev. 2 내용 유지 — 변경 없음]

### Phase 6 (New): 🟡 **계획** — Level Shift (B10) 완전 구현

[Rev. 2 내용 유지 — 변경 없음]

### Deprecated (구현하지 않음, 향후 hardware 업그레이드 후 재고)

[Rev. 2 내용 유지]

### Opt-in (환경변수로 활성화)

- **B6 huckel init_guess** — `UPDD_SCF_INIT_GUESS=huckel` 명시 시만
- **Unified Memory cascade** — 진단 용도로만, 성능 개선 없음
- **B7 mo_coeff warm-start** — T-Rank-1 smoke test 통과 후 opt-in
- **JoltQC mixed-precision SCF** — `UPDD_USE_JOLTQC=1` opt-in (**Rev. 3 신규**, QM/MM Acceleration Track 참조)

---

## 🚀 QM/MM Acceleration Track (Rev. 3 신규)

**버전 무관 track** — v0.7 calibration 이나 v0.8 publication 과 독립적으로 진행되는 가속 track. Phase number 를 부여하지 않음 (Phase 7-9 는 이미 v0.7 WT/Compstatin/Sa-D3 pilot 에 배정).

### Track 배경

2026-04-19 세션에서 **JoltQC evaluation side-quest** 로 시작. Crash 복구 대기 시간 활용으로 착수된 independent project (`/home/san/joltqc_updd_eval/`). 결과:

- **Axis 1 (Numerical consistency)**: 🟢 **DEFINITIVE PASS**
  - Stage 5 (ncaa_only 52 atoms): dE = 1.47×10⁻⁷ mHa
  - Stage 6 (peptide_only 289 atoms): dE = 7.1×10⁻⁶ mHa
  - Stage 7 (full 466 atoms): dE = 0.139 mHa (gate 0.5 mHa 3.6× margin)
- **Axis 2 (Speed/Efficiency)**: 2-4× speedup confirmed
  - 52 atoms: 3.96×
  - 289 atoms: 2.14×
  - 466 atoms: ~2× (indirect, DF vs direct SCF mixed)
- **Axis 3 (Physical accuracy)**: ⏸ **DEFERRED**
  - Stage 7 reference (UPDD snap01, charge=-3) 자체가 wrong-state (Phase 3 discovery)
  - Post-Phase-3 correct-charge re-run 필요 → Axis 3 verdict 가능

### 주요 scientific observations (FP32 noise floor 실측)

**Near-degenerate SCF convergence pathology** (Stage 7 trajectory):
- Storm (cycle 0-4): initial oscillation
- **Plateau (cycle 5-83, 79 cycles)**: |ΔD| 1-10 범위 oscillation, median 3.47 — FP32 noise floor evidence
- Breakthrough (cycle 84-130): DIIS subspace 가 eventually correct direction 찾음
- Asymptote (cycle 131-151): FP32 noise floor ~10⁻⁵ 에서 oscillate

**Implications**:
- JoltQC 가 FP32 noise + near-degenerate HOMO-LUMO gap 조건에서도 **robust** (eventually converge)
- Cycle cost 증가 +37% (UPDD FP64 111 cycles → JoltQC 152 cycles), 하지만 per-cycle time 이 ~2× 빨라 net speedup
- **Max cycle ≥ 200 권장** (near-degenerate design 의 plateau 허용)

### Integration Roadmap (3-stage, opt-in first)

**Stage 1 (현재 — Rev. 3)**: **Opt-in environment variable**
- `UPDD_USE_JOLTQC=1` 명시 시에만 활성화
- 기본 pipeline 변경 zero (v0.7 calibration 에 영향 없음)
- Power user 가 speed 필요할 때 case-by-case 사용
- Fallback: JoltQC 실패 시 FP64 로 자동 복귀

**Stage 2 (v0.8+, conditional)**: **Hybrid mode**
- Screening tasks: default JoltQC
- Calibration tasks: default FP64
- `run_qmmm.py` 에 `--mode screening|calibration` flag
- Precondition: Axis 3 DEFINITIVE PASS (correct-charge Stage 7 re-run)

**Stage 3 (v0.9+, future)**: **Default acceleration**
- 충분한 production 데이터 축적 후 검토
- Calibration 결과가 JoltQC 와 FP64 의 statistical equivalence 증명 후
- Default FP64 → default JoltQC 전환

### 진입 조건 (Track 전체)

**Required before Stage 1 opt-in deployment**:
- [ ] UPDD Phase 3 Stage 1 merged (✅ completed 4034ee6)
- [ ] UPDD Phase 3 Stage 2 merged (🔄 진행 중)
- [ ] Axis 3 preliminary verdict (Stage 7 re-run after correct-charge)
- [ ] SciVal approval for opt-in integration
- [ ] Fallback mechanism 검증 (JoltQC failure → FP64 automatic)
- [ ] Regression test (기존 converged FP64 result 재현)

**Required before Stage 2 hybrid**:
- [ ] Stage 1 production 데이터 >10 SCF runs
- [ ] At least 2 designs 에서 Stage 7 재현
- [ ] Calibration-relevant accuracy 증명 (ΔΔG ranking invariance)

### Constraints / Known limitations

- **JoltQC 는 DF (density fitting) 미지원** → direct SCF 로 전환 필요 (VRAM 부담 증가 가능)
- **JoltQC 는 UKS 미지원** → open-shell calculations 는 FP64 강제
- **JIT compile overhead**: 첫 실행 시 2-5분 추가 (kernel caching)
- **Near-degenerate overhead**: 약 +37% cycles (but per-cycle 2× 빠름 → net positive)

### Cross-project reference

- **External project**: `/home/san/joltqc_updd_eval/` (user-managed, document-only interface)
- **Re-run plan**: `external/STAGE7_RERUN_PLAN.md` (9 sections, Case A/B/C decision tree)
- **Decision log**: `external/DECISION_LOG.md` (3-axis framework application)
- **Known issues**: `external/KNOWN_ISSUES.md` (bidirectional cross-reference)

### Integration mode 결정 원칙

JoltQC 가 UPDD 에 어떻게 들어가는지는 **사용자가 환경변수로 통제**:

```bash
# Default (Rev. 3 현재): FP64 only
pixi run python run_qmmm.py ...

# Opt-in: JoltQC mixed precision
UPDD_USE_JOLTQC=1 pixi run python run_qmmm.py ...
```

**과학적 rationale**:
- FP64 가 default 인 이유: Calibration 에 영향 주지 않음, v0.7 진입 조건 유지
- Opt-in 인 이유: 사용자가 "지금은 speed 우선" 이라고 명시적 선언 필요
- Rev. 2 의 opt-in pattern (B6 huckel, B7 warm-start) 과 일관

---

## 🖥 하드웨어 최적화 Guidelines (Rev. 2 유지)

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
| **JoltQC mixed-precision (Rev. 3)** | **~8 GB** (Stage 7 실측) |

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
# SCF parameters (Rev. 2)
UPDD_SCF_CONV_TOL=1e-7              # B3 (ranking-safe)
UPDD_SCF_DIRECT_TOL=1e-12           # B5
UPDD_SCF_INIT_GUESS=minao           # B6 (huckel opt-in)
UPDD_DFT_GRIDS_LEVEL=2              # B2

# GPU / Memory (Rev. 2)
UPDD_MIN_AO_BLKSIZE=64              # B8 (gpu4pyscf import 전)
UPDD_VRAM_POOL_FRACTION=0.65        # B9 (cupy mempool)
UPDD_DF_AUTO_VRAM_THRESHOLD_GB=20   # DF auto threshold
UPDD_QM_MAX_MEMORY=16000            # (레거시, 참조만)

# QM/MM Acceleration (Rev. 3 신규)
UPDD_USE_JOLTQC=0                   # Opt-in JoltQC mixed-precision (0=off default, 1=on)

# Deprecated / No-op (참조는 유지)
UPDD_DF_JK_MAX_MEMORY_MB=4000       # gpu4pyscf 가 무시, 명시성 용도
```

---

## 🔁 진단 → 수정 → 재검증 루프 (Rev. 4 updated — Keeper inserted)

Rev. 2 루프에 Keeper 가 Coder ↔ Runner 사이 pre-execution gate 로 삽입됨:

```
[새 요청]
   ↓
Archi — 계획 수립
   ↓
SciVal — 사전 과학적 검증 (Science gate, R-4)
   ↓ 승인
Coder — 코드 수정 + UPDATE.md (R-8)
   ↓
Keeper — cross-stage integrity gate (Integrity gate, R-17) ← Rev. 4 NEW
   ↓ 🟢 pass 시
Runner — 실행 + 기본 검증 (Execution gate, R-5)
   ↓
Path — 물리/화학/병리학적 진단 (Pathology gate, R-6)
   ↓
🟢 통과 → user 보고 → git command 제시 (R-9)
🟡 조건부 통과 → user 상의
🔴 실패 → Archi 로 복귀
```

Keeper 가 🔴 block 한 경우: Runner 실행 안 함, Coder 로 되돌려 integrity 문제 해결 후 re-gate. 이 패턴이 cofactor-absence 류 silent drift 를 pre-SCF 에서 잡음.

Rev. 2 루프 제한 규칙 (Path 3회 연속 🔴 → escalate 등) 전체 유지.

---

## 🧪 Calibration (v0.7) 진입 조건 (Rev. 3 update)

### 필수 조건

v0.7 Compstatin/Sa-D3 calibration 에 진입하려면 v0.6 에서 다음을 달성:

1. **최소 3 design × 5 snap = 15 snap 이 `converged=True`**
   - **Rev. 3 추가**: snap 은 R-15/R-16 guard 통과한 것만 counted
   - **Rev. 3 WITHDRAWN**: w4_1_s2/snap01 (wrong-state) 는 제외
2. **각 design 의 5 snap 의 stdev < |median| × 0.3** (unified ensemble 기준)
3. **Design 간 ΔΔG_UPDD 순위가 v0.5 archive 와 Kendall τ ≥ 0.9**
   - **Rev. 3 caveat**: v0.5 archive 의 historical audit (Stage 3) 완료 후 재평가
4. **모든 converged snap 의 `binding_int_kcal` 이 물리적 범위 내** (Path 진단 🟢)
5. **Pre-QM filter 적용 — budget violation 은 0** (budget 600, Rev. 3 update)
6. **모든 snap 이 `charge_consistency_audit: passed` (R-15/R-16 guard)** — Rev. 2 추가
7. **B10 level_shift 가 필요한 design 은 일괄 적용 완료**
8. **Rev. 3 추가**: JoltQC integration track 은 **v0.7 진입 조건과 무관** (opt-in only, default FP64 유지)

### v0.7 에서 추가될 것

[Rev. 2 내용 유지]

### Rev. 2 의 Scope 축소

[Rev. 2 내용 유지]

---

## 📋 First Task for Archi (Rev. 3)

### 🚨 세션 시작 최초 Step (Rev. 3.1 추가)

**모든 UPDD 세션의 첫 행동 은 Archi 의 document boot sequence.** User 의 task 가 무엇이든, Archi 는 아래 순서대로 읽은 후에만 다른 agent (SciVal/Coder/Runner/Path/Keeper) 호출 혹은 작업 설계를 시작한다:

1. **`CLAUDE.md`** (이 파일) — Rev. 번호 + Core Principles + Absolute Rules + Phase 우선순위 확인
2. **`UPGRADE.md`** — roadmap 체크박스 현황 + v0.6.x 완료 상태 + v0.7 진입 조건
3. **`UPDATE.md`** — 최근 3-5 entries 만 훑어 직전 세션 변경 파악
4. **`UPDDnow.md`** (있다면) — 현재 구현 상태 snapshot
5. **`UPDDfuture.md`** (있다면) — 향후 로드맵 + SciVal pending verdicts

**이 boot sequence 를 건너뛰면** precondition drift, 중복 구현, completed 작업 재시도 등 **세션 시작 1-2 턴 내 구조적 오류** 가 필연적으로 발생. Skip 은 R-10 위반으로 간주.

**Boot sequence 예외 없음**:
- "간단한 수정" 이라도 skip 금지 (간단해 보이는 변경이 완료된 Phase 를 되돌릴 수 있음)
- User 가 "바로 X 하자" 라고 해도 Archi 는 1-2 턴 내 boot sequence 완료 후 응답
- 두 번째 이상 세션 시 CLAUDE.md 가 system-reminder 로 이미 주입된 경우에도 UPGRADE.md 는 반드시 별도 read

### ⚠️ 작업 시작 전 필수 루틴 (Rev. 3.1 추가)

**어떤 작업이든 착수 전, Archi 는 반드시 `UPGRADE.md` 를 먼저 읽고 task 설계에 반영한다.**

- `UPGRADE.md` 는 **진실의 원천 (source of truth)** 이며 현재 진행 중 / 완료 / 계획된 모든 upgrade 항목의 **체크리스트** 를 보유. 로컬 전용 문서 (`.git/info/exclude`, push 되지 않음).
- 작업 시작 시: 관련 Phase / 섹션의 체크박스를 읽고 **중복 여부, 전제 조건, 블로커 의존성** 을 확인. 예: `§A 1-trajectory Fix 4` 진입 전 Phase 3 Stage 3 audit 상태 확인.
- **`UPDATE.md` 작성 직후**: Archi 는 `UPGRADE.md` 의 해당 체크박스를 **반드시 `[x]` 로 전환** 하고 commit hash / 완료일자 기록. Skip 시 다음 세션에서 중복 구현 또는 precondition 누락 발생.
- R-10 (완료된 upgrade 는 `UPGRADE.md` 체크 표시) 의 **메커니즘 명시화**: UPDATE.md 엔트리 → UPGRADE.md 체크 전환이 한 쌍.

**루틴 요약**:
```
[작업 시작] Archi: UPGRADE.md 조회 → 관련 체크박스 + precondition 확인 → 작업 설계
[작업 중]   Coder/Runner/Path/Keeper/SciVal: 실행 + 검증
[작업 완료] Archi: UPDATE.md 상단에 엔트리 추가 (R-8)
            Archi: UPGRADE.md 해당 체크박스 [x] 로 전환 + commit hash / 날짜 기록 (R-10)
[검증]      다음 세션 Archi 첫 invocation: UPGRADE.md 읽으며 이전 세션 진행 확인
```

**Anti-pattern (금지)**:
- UPDATE.md 에만 기록하고 UPGRADE.md 미갱신 → precondition drift + 중복 구현 위험
- UPGRADE.md 를 읽지 않고 새 작업 설계 → 완료된 작업 재구현 또는 블로커 무시
- UPGRADE.md 업데이트 없이 SciVal / Keeper / Coder 에 task 분배

### 시작 지점 (Rev. 3 업데이트):

1. **Phase 3 Stage 2 Coder 완료 리뷰**
   - Coder background 작업 완료 대기
   - `utils/run_qmmm.py::partition_qmmm()` whole-residue path 검증
   - `target_cards/6WGN.json` v0.6.4 schema 확인
   - Smoke test 결과 (n_qm ~580, link atoms 5-7) 확인
2. **Path 에 Stage 1 post-merge 진단 요청**
   - 입력: 4034ee6 (Stage 1 merged) 의 R-15/R-16 guard 가 regression test 에서 정확히 동작하는지
   - 기대 output: Stage 2 방향 confirm or redirect
3. **Stage 2 merge 후 Phase 3 Stage 3 착수**
   - Historical audit execution (user explicit approval 필요, R-7)
   - `outputs/*/qmmm_results/*.json` 일괄 audit
   - PASS/PARITY_ONLY/FAIL 분류 → 격리 이동
4. **Phase 4 Pre-QM filter 업데이트**
   - Budget 500 → 600 반영 (Rev. 3)
   - w4_0_s2, w1_0_s1 재평가 (새 budget 에서 다시 eligible)
5. **GitHub issue tracking**
   - Main issue #XX (P0-critical parity) 진행 상황 업데이트
   - Sub-issue #YY (numbering offset): Stage 2 완료 시 close
   - Sub-issue #ZZ (comment fix): 간단한 cleanup PR 로 close
6. **snap02 clean restart (Stage 1+2 merged 후)**
   - R-15/R-16 guard 통과 확인
   - nohup 실행 (R-13 reserved 의 spirit)
   - v0.7 진입 counter 에 포함

---

## 🔄 UPDATE.md / UPGRADE.md 관리 (Rev. 3.1 확장)

### 두 문서의 역할

- **UPDATE.md** (change log, 로컬 전용): 모든 코드 수정의 **상세 엔트리** 를 상단 역순 시간으로 기록. R-8 의무사항.
  - 포맷: `## [v0.6.x] YYYY-MM-DD — Brief Title` + **Agent / Category / Regime / 변경 사항 / 이유 / 검증 / Regime impact / 관련** 8-field.
  - 문헌 출처 (DOI) 반드시 포함.
  - No Co-Authored-By, no subagent names (feedback memory).

- **UPGRADE.md** (roadmap checklist, 로컬 전용): 모든 Phase / Fix / Milestone 의 **체크박스 진척도**.
  - Rev. 3 현재: v0.6 (대부분 완료) + v0.6.8 (MM-GBSA 1-traj migration) + v0.7 (Calibration) + v0.8 (Expansion) + v0.9+ (Future).
  - 체크박스 `[x]` 전환 시 **commit hash + 완료일자** 병기.
  - Hardware-Aware Design 섹션 + Deprecated 섹션 유지.
  - 새 Phase/Milestone 추가 시 **Archi 가 SciVal 사전 검증 후** 등재.

### Archi 의 **필수 루틴** (Rev. 3.1 추가)

**1) 작업 시작 전**:
   - `UPGRADE.md` 를 최초 조회 (grep / head 도 가능, 전체 read 권장)
   - 관련 Phase / 섹션의 체크박스 status 확인
   - **Precondition dependency 매핑** (예: §A Fix 4 → Phase 3 Stage 3 audit → §B Switch RMSD)
   - Blocker 체크: `UPDDfuture.md` §Q Risk Register + UPDDnow.md §11 미해결 이슈 교차 참조
   - Task 설계 시 중복 방지: 이미 `[x]` 인 항목은 **재구현하지 않음**

**2) 작업 완료 직후** (코드 merge 또는 commit):
   - `UPDATE.md` 상단에 신규 엔트리 작성 (R-8)
   - **동일 session 내** `UPGRADE.md` 의 해당 체크박스 `[x]` 로 전환
   - 새 completed 엔트리는 UPGRADE.md `## Completed ✅` 섹션에도 추가
   - 만약 새 Phase / Milestone 가 등장했다면 UPGRADE.md 에 TBA (To Be Approved) 로 추가 후 Archi + SciVal 승인 대기

**3) 세션 종료 직전**:
   - UPDATE.md 와 UPGRADE.md 의 entries 가 **일관** 한지 확인 (UPDATE.md 에 있는데 UPGRADE.md 미반영 등 drift 방지)
   - 다음 세션을 위한 **handoff note** 를 UPGRADE.md 의 "Rev. N 다음 업데이트 트리거" 섹션에 업데이트

### Anti-pattern (금지)

| ❌ | 이유 | 복구 |
|---|---|---|
| UPDATE.md 만 기록 + UPGRADE.md 스킵 | 다음 세션에서 완료 여부 불명 → 중복 구현 위험 | 즉시 UPGRADE.md `[x]` 전환 |
| UPGRADE.md 읽지 않고 task 설계 | 완료된 작업 재구현 / precondition 누락 | `First Task for Archi` 의 "시작 전 필수 루틴" 재확인 |
| Anonymous 체크박스 `[x]` (commit hash/date 생략) | Audit trail 유실 → Keeper integrity gate 에서 drift 감지 불가 | 체크박스에 `— v0.6.x (커밋 abcd123)` 병기 필수 |
| SciVal 미승인 Phase 를 UPGRADE.md 에 직접 추가 | R-4 위배 + Archi 의 읽기-only 권한 위반 | TBA 로 표기 + SciVal 승인 verdict 파일 링크 |

### Entry format (Rev. 2 유지 + Rev. 3 예시)

**UPDATE.md 엔트리**:

```markdown
## [v0.6.7] 2026-04-21 — Generic amber14 parent walker + 25 ncAA parameterization

**Agent**: Coder (SciVal 🟢 `verdict_ncaa_strategy_20260420.md`)
**Category**: Feature (ncAA pipeline)
**Regime**: ranking-only (R-11)

### 변경 사항
- `utils/parent_topology.py` 신규 (510 lines)
- `utils/parameterize_ncaa.py` Phase 2 generic walker dispatch
- 25 ncAA registry entries (NMA/NML/.../PTR) 에 parent_residue 필드

### 이유
- Rev. 2 의 keep_atoms=STANDARD_BACKBONE 이 MTR 등 sidechain ncAA 에서 template match 실패
- Capece 2012 Strategy A 의 hybrid amber14SB + GAFF2 methyl graft 패턴

### 검증
- MTR 4W9A compstatin 2ns MD SUCCESS (94.71s)
- Cp4 MM-GBSA ⟨ΔG⟩ = -13.63 kcal/mol (2QKI ladder 2-point)
- 25/36 tests pass

### Regime impact
- Ranking-regime safe: Capece 2012 per-residue < 1 kcal/mol bias
- Absolute calibration 은 v0.7 Phase 8 에서 재평가

### 관련
- UPGRADE.md v0.6.6 → v0.6.7 [x] ncAA Universal Pipeline 전체
- UPDDnow §4 / UPDDfuture §M
- Commits: dc46690, 05944fc, 584d3e0, 7075416, 9c3f8a5, a4c1047, 1333e17, f8fe6f9
```

**UPGRADE.md 체크박스 업데이트 예시**:

```diff
- [ ] `utils/parent_topology.py` — amber14 generic walker
+ [x] `utils/parent_topology.py` — amber14 generic walker (v0.6.7, commit `dc46690` 2026-04-21)
```

**Completed 섹션 추가 예시**:

```markdown
## Completed ✅ (Rev. 3 누적)

### v0.6.7 (2026-04-21)
- [x] `dc46690` Generic amber14 parent walker + 14 ncAA parent_residue
- [x] `05944fc` Phospho multi-atom chained extensions + halogens + D-aa
...
```

---

### 연동 도식

```
┌─────────────┐      ┌──────────────┐      ┌──────────────┐
│  Archi 작업  │─────▶│  UPDATE.md   │─────▶│ UPGRADE.md   │
│  시작 (read) │      │  (append)    │      │ (check [x])  │
│             │      │              │      │              │
│ UPGRADE.md  │◀─────│ 다음 세션 조회 │◀─────│ handoff note │
│  조회       │      │              │      │              │
└─────────────┘      └──────────────┘      └──────────────┘
```

이 싸이클이 **깨지면** Phase drift + precondition 누락 + 중복 구현 발생.
이 싸이클이 **유지되면** UPDD 의 v0.6 → v0.7 → v1.0 migration 이 결정론적.

---

## v0.6 — Physical/Chemical Validity Recovery (진행 중)

### Phase 1: Topology-based QM Region
[Rev. 2 내용 유지]

### Phase 2: 3-term Supermolecular Decomposition
[Rev. 2 내용 유지]

### Phase 3: Charge Correction Policy (🔄 진행 중)

**Stage 1**: ✅ merged (4034ee6)
- R-15 Charge Declaration SSOT guard
- R-16 Binder Chemistry SSOT guard
- `utils/charge_topology.py` 신규
- `utils/audit_charge_consistency.py` 신규 (dry-run default)
- 23/23 tests pass

**Stage 2**: 🔄 진행 중
- target_card schema v0.6.4 (numbering_convention, whole-mode partitioning)
- Budget 500 → 600 완화
- Auto-generated rationale

**Stage 3**: ⏸ Stage 2 후 착수
- Historical audit (outputs/*/ + v0.5 archive)

### Phase 4: Pre-QM Budget Enforcement + Charge Parity Filter (신규)
- Budget 500 → 600 (Rev. 3)
- R-15/R-16 pre-filter 통합

### Phase 5: MM-GBSA + DFT Hybrid (신규, Phase 3+4 precondition)
[Rev. 2 내용 유지]

### Phase 6: B10 Level Shift (대기)
[Rev. 2 내용 유지]

### Bundle B1-B9 Status (2026-04-19)
[Rev. 2 내용 유지]

---

## v0.7 — WT Counterfactual + Calibration

### 진입 조건 (v0.6 에서 달성)
- [ ] 3+ design × 5 snap converged
- [ ] Path 진단 🟢
- [ ] Pre-QM filter 적용 (budget 600)
- [ ] R-15/R-16 guard: 모든 snap `charge_consistency_audit: passed`
- [x] **Rev. 3**: R-15/R-16 Stage 1 merged (4034ee6)

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

---

## v0.8 — Calibration Expansion + Publication Track
- [ ] 15-IgBP pilot
- [ ] MD replicates (3x)
- [ ] Bootstrap CI
- [ ] Paper draft
- [ ] **Rev. 3 추가**: QM/MM Acceleration Track Stage 2 검토 (hybrid mode, 조건부)

---

## Phase 8-13 (Reserved — Rev. 2 참조)

Rev. 2 의 Phase 8-13 상세 계획은 `archive/CLAUDE_Rev2_20260418.md` 에 보존. Rev. 3 에서는 변경 없음 — v0.6 진행 완료 + v0.7 pilot 착수 후 Rev. 3.1 또는 Rev. 4 에서 재평가.

---

## Completed ✅

- [x] v0.5 SSOT refactor (2026-04-17)
- [x] v0.5 Bug 1-3 fix (2026-04-17)
- [x] Path agent 도입 (2026-04-19 Rev. 1)
- [x] v0.6.0 QM/MM physical validity (2026-04-18 14b0295)
- [x] SCF bundle B2/B3/B5/B6/B8/B9 (2026-04-19 uncommitted)
- [x] DF auto-mode (2026-04-19 uncommitted)
- [x] Monitoring tools regex fix (2026-04-19 uncommitted)
- [x] CLAUDE.md Rev. 2 반영 (2026-04-19)
- [x] R-15/R-16 Charge guards + audit module Stage 1 (2026-04-19, commit 4034ee6)
- [x] **Rev. 3 추가 (2026-04-19 저녁)**:
  - [x] JoltQC evaluation side-quest Stage 0-7 complete (external project)
  - [x] SciVal verdict chain v1/v2/v3 for charge policy
  - [x] GitHub issue filing (main #XX + sub-issues #YY, #ZZ)
  - [x] snap01 archive to `_archive/v0.6_wrong_charge/` (manifest 포함)
  - [x] Chkfile archive to `_archive/chkfile_benchmarks/` (SCF profiling reference)
  - [x] CLAUDE.md Rev. 3 published

---

## 🚀 Git Commit & Push Policy (Rev. 2 유지)

[Rev. 2 내용 유지 — User-Controlled Push + 커밋 구조 제안]

---

## 📞 Agent Invocation 예시 (Rev. 2 유지)

[Rev. 2 Archi/Path/SciVal 호출 예시 전체 유지]

---

## 💾 Memory Hygiene (Rev. 3 update)

### Path 의 memory 카테고리
[Rev. 2 유지]

### Path 가 기록할 Rev. 2 패턴 예시
[Rev. 2 유지]

### Rev. 3 추가: SciVal verdict chain layer

`.claude/agent-memory/scival/` 에 verdict chain 저장:
- `verdict_target_card_charge_rationale_20260419.md` (v1 — initial diagnosis)
- `verdict_target_card_charge_rationale_20260419_v2.md` (v2 — Option B3 decision)
- `verdict_target_card_charge_rationale_20260419_v3_addendum.md` (v3 — historical disposition)

각 verdict 는 **cross-referenced**: GitHub issue 에서 인용, Path 진단에서 참조, CLAUDE.md 에서 link.

### Cross-project memory isolation

Agent 의 memory 범위는 **project scope 로 한정**. Cross-project memory sharing 은 **document-only** (SciVal verdict 파일 같은 shared document) 경로로만. Agents 는 서로 직접 통신하지 않음.

---

## ⏹ Stop Conditions (Rev. 3 update)

[Rev. 2 내용 유지 + Rev. 3 추가:]

### Rev. 3 추가 Stop Conditions

- **R-15/R-16 guard fail-fast**: 기존 converged result 에서 Mismatch 감지 시 즉시 halt, SciVal 에스컬레이션
- **Partitioning parity violation 재발**: Stage 2 merge 후에도 parity contradiction 발견 시 Stage 2 rollback
- **Cross-project dependency conflict**: 외부 프로젝트 (JoltQC eval 등) 결과가 UPDD 의 계획과 상충하면 user 상의 필요

---

## 📊 Scope 재조정 — v0.6 → v0.7 → v0.8 (Rev. 3)

### v0.6 (현재, ~2026-04-30 목표)
- Phase 3 Stage 1 ✅ merged (2026-04-19)
- Phase 3 Stage 2 🔄 진행 중
- Phase 3 Stage 3 ⏸ Stage 2 후
- Phase 4 확장 진행 중 (budget 600)
- Phase 5-6 계획

### v0.7 (~2026-05-15 목표)
- Compstatin + Sa-D3 pilot
- Calibration 2-3 points
- **Rev. 3**: QM/MM Acceleration Track Stage 1 (opt-in) 검토 시작

### v0.8 (~2026-06-15 목표)
- 15-IgBP pilot
- MD replicates (3x)
- Bootstrap CI + paper draft
- **Rev. 3**: JoltQC hybrid mode 조건부 배포 (Axis 3 DEFINITIVE PASS 시)

### v0.9+ (미래)
- H100 upgrade path
- **Rev. 3**: JoltQC default integration (v0.8 production 데이터 기반)

---

## 📜 Rev. 3 후기

Rev. 2 에서 R-15/R-16 설계를 명시한 지 하루 만에, JoltQC evaluation side-quest 중에 그 guard 의 필요성이 **실증적으로 확인**되었다. w4_1_s2/snap01 의 1-electron deficit 은 silent pass 의 대표 case 이며, parity check + magnitude check 의 이중 방어의 정당성을 보여준다.

이번 incident 는 UPDD 의 correctness 철학의 정당성을 확증한다:
- ✅ Pre-QM filtering 의 중요성 (Phase 4)
- ✅ Runtime verification 의 필수성 (R-15/R-16)
- ✅ Historical audit 의 가치 (Stage 3, v0.5 archive 도 포함)
- ✅ Multi-agent verification (SciVal veto chain 3회)

또한 JoltQC eval 은 UPDD 에게 **외부 관점의 검증 도구**를 제공했다. 향후 유사한 side-quest 가 다른 숨은 bug 를 발견할 수 있도록, cross-project document 채널 (`docs/known_issues.md`, SciVal verdicts) 을 지속적으로 관리한다.

**다음 Rev. 4 예상 시점**: v0.6 Phase 3 완전 merge 후 (Stage 1+2+3 완료, historical audit 결과 통합).
