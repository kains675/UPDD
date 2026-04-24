# UPDD Upgrade Roadmap (Rev. 3.1)

_Rev. 3.1 (2026-04-24): Dual-file protocol added. Pathway E v0.7 roadmap extracted to parallel file `UPGRADE_v07_roadmap.md`. This file retains authority for v0.6 Phases, v0.6.8 Foundation, v0.7 Entry YAML, Phase 7-11, Completed history, Hardware/Coordination._
_Rev. 3 (2026-04-21): `CLAUDE.md` Rev. 3–4 + `UPDDnow.md` (v0.6.7, main `f8fe6f9`) + `UPDDfuture.md` rev 0.2 + Archi/SciVal/Coder 통합 consult._

**Dual-file protocol (Rev. 3.1)**:
- This file: v0.6 phases, v0.6.8 Foundation, v0.7 Entry YAML, Phase 7-11, Completed history, Hardware, Coordination, general roadmap
- `UPGRADE_v07_roadmap.md`: Pathway E (α/β/γ) phase-based entry/exit criteria + data flow + decision audit
- Conflict resolution: Pathway E items → `UPGRADE_v07_roadmap.md` wins; non-Pathway-E items → this file wins
- Both files updated together for Phase α/β/γ execution changes
- Archive backup: `archive/UPGRADE_rev3_pre_dual_file_20260424.md`

_Rev. 2 archive: `archive/CLAUDE_Rev2_20260418.md` 참조 (2026-04-19 기반)._

---

## v0.6 — Physical/Chemical Validity Recovery (대부분 완료, v0.6.7 shipped)

### Phase 1: Topology-based QM Region — ✅ 완료
- [x] `partition_qmmm()` topology mode (`run_qmmm.py`) — v0.6.1
- [x] `target_card.json` schema v0.6.2 — v0.6.2
- [x] `generate_target_card.py` 자동 생성기 — v0.6.2
- [x] Default mode 변경 — v0.6.2
- [x] **Rev. 3 추가**: `bicyclic` topology 4-branch 확장 (add_cyclic_bond + junction exclusion + cyclic_relaxation + Stage 1.5) — v0.6.6 (1SFI 검증, σ 15.18→8.28)
- [x] **Rev. 3 추가**: `numbering_convention` field (md_to_crystal_offset=3 for 6WGN) — Phase 3 Stage 2 merge
- [x] **Rev. 3 추가**: Budget 500 → 600 완화 — Rev. 3 decision

### Phase 2: 3-term Supermolecular Decomposition — ✅ 완료
- [x] `qm_int_kcal_frozen` 필드 (3-term 분해) — v0.6.5
- [x] `interaction_kcal` 보존 — v0.6.5
- [x] Senn-Thiel H-cap link atoms + parity tracking — v0.6.5
- [x] Baseline measurement (w4_1_s2/snap01 qm_int = −104.256 kcal/mol)

### Phase 3: Pre-SCF Charge Guard (R-15/R-16) — ✅ Stage 1+2 merged, ⏸ Stage 3 대기
- [x] **Stage 1** (commit `4034ee6`): `utils/charge_topology.py` + `utils/audit_charge_consistency.py` + run_qmmm guard, 23/23 tests
- [x] **Stage 2** (merge): Schema v0.6.4 (numbering_convention + whole-mode partitioning), auto-generated rationale
- [ ] **Stage 3 historical audit** (blocking v0.7 entry)
    - [ ] `utils/phase3_historical_audit.py` — audit_charge_consistency composer + git commit-era classifier
    - [ ] Output classifier: `{PASS, PARITY_ONLY, F2_AFFECTED, F1_AFFECTED, P1_AFFECTED, FAIL}` (UPDDfuture §D.1)
    - [ ] `outputs/_archive/v0.6_pre_1traj/` 신규 격리 (Fix 4 커밋 시)
    - [ ] 12-16h rerun for FAIL classified
    - Depends on: §A 1-traj Fix 4 완료 (protocol lock 후 audit)

### Phase 4: Pre-QM Filter — ✅ 완료
- [x] `utils/pre_qm_filter.py` 모듈 (SciVal approved)
- [x] `estimate_n_qm(target_card, binder_chain_len)` 함수
- [x] MD 진입 전 filtering 로직
- [x] `rejected_designs.json` 출력
- [x] 탈락 이유 분류: `n_qm_exceeded` / `chemistry_parity_incompatible` / `contact_residue_count`
- [x] Budget 500 → 600 완화 (Rev. 3)

### Phase 5: MM-GBSA + DFT Hybrid (🟡 v0.7 계획, Fix 4 의존)
- [ ] `utils/cluster_snapshots.py` — k-medoids representative selection
- [ ] `run_qmmm.py --representatives` 옵션
- [ ] `rank_results.py` hybrid integration (weighted ΔG)
- [ ] Uncertainty propagation (MM-GBSA stderr + DFT variance)
- [ ] **Precondition**: §A Fix 4 완료 + Phase 3 Stage 3 audit completed

### Phase 6: B10 Level Shift (⏸ 조건부 대기)
- [ ] Gate 조건: HOMO-LUMO near-degeneracy 재발 시. Path memory `pathology_homolumo_near_degeneracy_plateau.md` baseline 참조

### R-17 Cofactor Preservation (Stage 4, schema 0.6.5) — ✅ 완료 (Rev. 3 신규)
- [x] Schema upgrade path (`utils/utils_common.py::_upgrade_cofactor_entries_to_v065`)
- [x] Error classes (`utils/cofactor_errors.py`): CofactorMissing/ChargeSumMismatch/ParamMissing/ReinsertionError/CoordFrameError
- [x] **G1** preprocess cofactor presence (`utils/preprocess_target.py`)
- [x] **G2** AF2 reinsert Kabsch + RMSD ≤ 2.5 Å (`utils/reinsert_cofactors_post_af2.py`)
- [x] **G3** parameterize cofactor (`utils/parameterize_cofactor.py` — GAFF2 + Li-Merz 12-6-4)
- [x] **G4** MD cofactor bond inject (`run_restrained_md.py::inject_cofactor_bonds` + `register_cofactor_templates`)
- [x] **G5** snap strip cofactor-aware (`extract_snapshots.py::strip_solvent`)
- [x] **G6** QM/MM inject from snapshot (`run_qmmm.py::inject_cofactor_point_charges_from_snapshot`)
- [x] `target_card generator` v0.6.5 default + `_COFACTOR_FF_LIBRARY`
- [x] 15 tests (`tests/test_r17_cofactor_preservation.py`), 66/66 total
- [x] Keeper integrity memory: `integrity_cofactor_coord_frame.md` (declaration ≠ placement archetype)

### ncAA Universal Pipeline (v0.6.6–v0.6.7 신규, Rev. 3 신규 Phase)
- [x] **Generic amber14 parent walker** (`utils/parent_topology.py`, 510+ lines)
    - [x] VF2 graph isomorphism (parent template ↔ GAFF2 mol2)
    - [x] `_load_parent_template`, `parent_heavy_bonds`, `parent_neighbor_hint`, `parent_h_name_table`
    - [x] ACE/NME cap + free-acid OXT 자동 감지
    - [x] Multi-atom chained extension (P + 3 phosphate O's)
    - [x] GAFF type element override (`_atom_element`, `_GAFF_TYPE_ELEMENT`)
- [x] **Tier-1 ncAA 확장** (25/36 parameterized, commit `dc46690` + `05944fc`)
    - [x] N-methyl: NMA, NML, NMV, NMK, NMR, NMQ, MEA, SAR
    - [x] D-amino: DVA, DLE, DTR, DAR, DPR, DPN, DAL, DTY
    - [x] Side-chain: HYP, CHA, PFF, CLP, BRP, MTR
    - [x] Phospho PTMs: SEP, TPO, PTR
- [x] **CONECT integrity chain** (4 commits)
    - [x] `584d3e0` — chained-extension mutate geometry + unique H naming (HP1/H1P/H2P)
    - [x] `9c3f8a5` — external peptide-bond CONECT emission + DTY/SAR migration
    - [x] `a4c1047` — split_complex CONECT bucket-wise preservation
    - [x] `1333e17` — snap CONECT auto-injection (`_inject_ncaa_conect`)
- [x] **Hydrogen hygiene** (`f8fe6f9`)
    - [x] `_strip_spurious_ncaa_hydrogens` unconditional (H9/HX1 from free-acid NH2)
    - [x] snap atom-whitelist (pre-fix legacy H 제거)
- [ ] **11 deferred** (Silicon TMS/SIP/DPS, chain-length variants ORN/DAB/HAR/NLE, fused-ring TIC, Cα-quaternary AIB, alkyne PRA, azide AHA) — UPDDfuture §M

### MD Robustness — ✅ 완료 (v0.6.6–v0.6.7 신규)
- [x] DCD streaming to HDD (`run_restrained_md.py`)
- [x] GAFFTemplateGenerator cofactor registration (R-17 G4)
- [x] **Pass 1b seed-variant retry** (commit `7075416`): seed+13 auto-retry before dt=1fs fallback, 모든 topology 적용
- [x] **k=100 reset after Stage 1** (commit `09660d4`): cyclic_ss/linear ncAA NVT NaN fix
- [x] **Bicyclic topology 3-branch extension** (commit `09660d4`)
- [x] addSolvent compartment-aware (K⁺ intra / Na⁺ extra)
- [x] TIP3P-FB water, 0.15 M NaCl, padding 1.2 nm

### MM-GBSA Robustness — ✅ Rev. 3 Option α/β + P1 해결
- [x] **F1** Mg²⁺ strip 제거 (`_MMGBSA_STRIP_RESNAMES` 정리): 구조 cofactor 유지
- [x] **F2** Cyclic HTC disruption hack 철회 (`write_temp_pdb` +5Å+OXT 제거)
- [x] **P1** Under-minimization 해결: 500→5000 iter (env `UPDD_MMGBSA_MIN_ITER`, tol 10→1)
- [x] CustomGBForce parameter-ordering fix (GBn2: `charge/or/sr` 정확 처리, Mg/Zn/Ca/Mn Born radius)
- [x] amber14 12-6-4 ion FF (`amber/tip3p_HFE_multivalent.xml`)
- [x] Metal GB radius mbondi2 (Mg 1.18 / Zn 1.09 / Ca 1.37 / Mn 1.13 / Fe 1.08 / Cu 1.00 / Co 1.03 / Ni 1.00)
- [x] CUDA/CPU platform env var (`UPDD_MMGBSA_PLATFORM`)
- [x] split_complex + write_temp_pdb CONECT 전달 (commit `a4c1047`)

### Bundle B1-B9 (2026-04-19 기준, 변경 없음)
- [x] B2 `grids_level=2` — 2026-04-18
- [x] B3 `conv_tol=1e-7` — 2026-04-18
- [x] B5 `direct_scf_tol=1e-12` — 2026-04-18
- [x] B6 `init_guess=minao` (huckel opt-in) — 2026-04-19
- [x] B8 `min_ao_blksize=64` — 2026-04-18
- [x] B9 `cupy mempool 0.65` — 2026-04-18
- [x] DF auto-mode + threshold env var — 2026-04-18
- [-] **B1 deprecated** (rks_lowmem QMMMSCF 비호환)
- [ ] B7 conditional (T-Rank-1 smoke test 통과 후 opt-in)
- [ ] B10 pending (HOMO-LUMO near-degeneracy 재발 조건)

### v0.6 에서 마친 8 Target × ncAA End-to-End 검증 (Rev. 3 신규)
- [x] 2QKI (complement C3) cyclic_ss + compstatin W4→MTR (Cp4, MD+MM-GBSA, ⟨ΔG⟩=−13.63 w/ 5 snap 2 ns)
- [x] 1EBP (cytokine-R EPO) cyclic_ss + EMP1 W13→MTR (MD, 48.63s)
- [x] 7TL8 (surface) linear + Sa-D3 W6→MTR (MD, 97.53s)
- [x] 6WGN (oncogene KRAS-G12D) cyclic_htc + NMA + GNP+Mg²⁺ (MD, 66.45s; prior 4% 폭발 recovered)
- [x] 1YCR (PPI-oncology MDM2) linear + p53 L22→NML (MD, 20.46s)
- [x] 3IOL (Class-B GPCR GLP-1R) linear + GLP-1 L20→NML + Y19→PTR (MD both, 40.34s / 39.75s)
- [x] 2QKH (Class-B GPCR GIPR) linear + GIP W25→MTR (MD seed=7, 39.95s)
- [x] 1SFI (Ser protease β-trypsin) bicyclic + SFTI-1 WT (MD σ=8.28)

### v0.6.7 → v0.6.8 Calibration Queue (🔄 진행 중 2026-04-21 09:06 launch)
- [x] Queue script `/tmp/calib_queue_all.sh` (6 target × WT+ncAA × 3 rep × 5 ns × 10 snap = 36 runs)
- [x] Pass 1b seed retry 포함 (각 run 자동 복구)
- [x] Run 1 완료: 2QKI_WT_calib_s42 (재시도 seed=55) ⟨ΔG⟩ = −47.75 ± 17.44 (n=10)
- [ ] Run 2-36 진행 중 (ETA ~24-30h)
- [ ] 완료 후 SciVal + Keeper 검증 → v0.7 calibration entry 결정

### Deprecated (Rev. 2 에서 구현 안 함, 유지)
- [-] **DF-J/K on 16GB VRAM** (Wu 2024 arXiv:2404.09452 최대 168 atoms)
- [-] **rks_lowmem.RKS 전환** (QMMMSCF 구조적 비호환)
- [-] **Budget 500 → 300 축소** (Kulik 2016 위배)
- [-] **Huckel default init_guess** (w4_1_s2 발산)
- [-] **`mf.with_df.max_memory` 기반 J/K blksize** (gpu4pyscf 미참조)
- [-] **Rev. 3 신규**: **pure GAFF2 backbone for ncAA** (Khoury 2013 §2.3 backbone bias 3-5 kcal/mol; opt-in patch `UPDD_MTR_AMBER14_PATCH=1` 로만 amber14 type+charge overlay)

---

## v0.6.8 — MM-GBSA 1-Trajectory Migration + Audit (🔴 P0, Rev. 3 신규, UPDDfuture §A-D)

**Precondition**: UPDDfuture §A.4 Manual test 결과 Case A/B/C 결정 (Day 1). SciVal 🟡 CONDITIONAL APPROVE 에 따른 gate 수정.
**Status update (2026-04-22, v0.3.1)**: Phase 3 1-traj pilot COMPLETED (Runner 10/10 success); Case A both targets (2QKI_WT + 3IOL_WT); **Fix 4 GLOBAL deployment APPROVED** (SciVal `verdict_phase3_1traj_case_abc_20260422.md`). §A.2 Day 2-3 code implementation UNBLOCKED.

### v0.6.8 Foundation (2026-04-22 scaffold) — Post-Cp4-incident Route A' enablement
- [x] **Forward Plan Rev 0.2 integration to UPDDfuture v0.3.1** (2026-04-22, Archi in-place edit; `UPDATE.md` entry "UPDDfuture v0.3.1")
- [x] **4 user directive Fixes applied** (2026-04-22: Fix 1 "protocol error" retraction, Fix 2 Route A/A'/B/C 30-35% cons / 45-50% opt, Fix 3 §A.4 2×5 empirical, Fix 4 §G Gate 1-5 multi-layer)
- [x] **Phase 3 empirical 1-traj test** (2026-04-22, Runner 10/10 success; SciVal 🟢 verdict Case A both targets; `outputs/analysis/phase3_1traj_20260422/summary_20260422.json`)
- [x] **Agent Rev. 5 hardening complete** (2026-04-22: Keeper AGENT-01/02/03 + Path SRE active; `UPDATE.md` entry "Agent Rev. 5")
- [x] **Coder §A.2 Fix 4 1-traj implementation** (v0.6.8, 2026-04-22, commit `TBD`)
  - Opt-in via `UPDD_MMGBSA_PROTOCOL=1traj` (default `3traj` preserved for backward compat)
  - Regression 2QKI_WT_s42 snap01 ΔG_1traj = **−65.17 kcal/mol** (target −65.62 ± 0.5, \|Δ\|=0.459) ✅
  - σ reduction **7×** (n=10 production, σ_1traj=2.49 vs σ_3traj=17.44): Gate 2 noise_floor 🟢 resolved
  - Keeper 🟢 INTEGRITY PASS (8 gates, 0 at-risk consumers) + SciVal 🟡 CONDITIONAL MERGE APPROVED
  - Default flip deferred to post-Option Y (Gate 3 Student-t df=4 requires 5 replicates)
  - UPDATE.md entry: `[v0.6.8] 2026-04-22 — MM-GBSA 1-trajectory protocol (Fix 4, opt-in)`
- [x] **SciVal §A.4 post-Fix-4 regression verdict** (2026-04-22, `verdict_fix4_regression_v068_20260422.md`) — CONDITIONAL MERGE APPROVED, algorithmic parity verified, Case A signature preserved at n=10
- [x] **Coder §P0.2 MAD-any-of production merge** (2026-04-22, `mad_any_of_production_merge_20260422.md`) — `utils/mmgbsa_outlier_screen.py` 378 LOC + CLI wrapper + 79 sanitized JSON companions; schema field `sanitization_policy=mad_any_of_v1`; non-destructive (originals immutable); 3IOL_NML20_s101 raw +1.14 → sanitized −49.85 reproduction confirmed
- [x] **Runner §P0.1 2QKI_Cp4 Fix 4 5-rep rerun** (2026-04-22, `run_cp4_fix4_rerun_20260422.md`) — ⟨⟨ΔG⟩⟩ = −46.335, σ_between = 9.653 (4 existing DCDs + s19_reseed55 84% PARTIAL_SUCCESS recovery, total ~54 min GPU)
- [x] **Runner §P0.3 2QKI_WT Fix 4 5-rep rerun** (2026-04-22, `run_wt_fix4_rerun_p03_20260422.md`) — ⟨⟨ΔG⟩⟩ = −59.168, σ_between = 4.788 (5 seeds s7/s19/s42/s23/s101, all RC=0, 18 min GPU total)
- [x] **SciVal M3 sign flip formal verdict — Option B CONDITIONAL Case A** (2026-04-22, `verdict_m3_sign_flip_final_20260422.md`) — M3 ΔΔG = +12.833 ± 10.775, z_vs_lit = 1.558 < 2 literature-consistent; default flip CONDITIONAL APPROVE; n=10 required for unconditional; target-selective σ insight memory added
- [ ] **Default protocol flip** 3traj → 1traj (CONDITIONAL APPROVED per SciVal Option B; unconditional pending §A.8 n=10)
- [ ] 6WGN (GNP+Mg²⁺) 1-traj live test (post calibration-set completion; Keeper G6 structural-only verified)
- [x] Runner §A.3 sampling expansion (Option Y, 2026-04-22, 18/18 runs completed) — 5-rep achieved for 7 of 9 target_variant groups; Path Option Y verdict `pathology_option_y_5rep_20260422.md`; WT-absolute regime 9/9 PASS, ΔΔG 1/4 (2QKI_Cp4 [+3.49, +35.03])
- [ ] Phase 3 Stage 3 historical audit execution (blocks v0.7 entry; KRAS 6WGN w1_0_s1 re-evaluation deferred to v0.4)

### v0.7 Entry Prep (2026-04-22 Option Y results) — Post-Option-Y roadmap

> **Rev. 3.1 (2026-04-24)**: Track A/B outcomes complete. Pathway E adopted for v0.7 resolution strategy — see `UPGRADE_v07_roadmap.md` for phased α/β/γ detail. Track A/B/C/D below updated with final status + phase mapping.

- [x] **5 WT targets absolute ranking regime READY** (5-rep CI excludes 0): 1YCR_WT −70.03, 3IOL_WT −49.34, 2QKH_WT −45.14, 2QKI_WT −39.83, 1EBP_WT −37.49 (Path Option Y §7.1)
- [x] **2QKI_Cp4 ΔΔG Gate 3 PASS** (first ncAA pair in UPDD history): +19.26 [+3.49, +35.03] kcal/mol (Path Option Y §3)
- [x] UPDDfuture.md v0.3.1 → v0.4 synthesis (Archi, 2026-04-22 저녁) — §v0.4 Changelog + §A.7/A.8 + §G post-Option-Y + §K.6 Case X + §Q R10/R11
- [x] **2QKI_Cp4 Fix 4 rerun with same-protocol WT baseline (P0.1 + P0.3)** **P0** (2026-04-22, Runner) — both arms 5-rep complete; M3 ΔΔG = +12.833 ± 10.775 (z=1.558 literature-consistent); success criteria (§A.8): Genheden shift −25.77 kcal/mol for Cp4 ≥ 6 required ✅; σ 7× tighter at n=10 pattern (Fix 4 for Cp4) ✅; z_vs_lit 1.56 ≤ 2 ✅; Gate 3 CI at borderline (−0.54 excludes 0 marginally); SciVal M3 verdict: **Option B CONDITIONAL Case A**
- [ ] **2QKI_Cp4 + WT n=10 replicate expansion (§A.8, X.A' → X.A promotion)** **P1 DEPRIORITIZED** (v0.5 NEW, v0.6 deprioritized under X.B) — v0.6 Fix 4 symmetric n=10+n=8 aggregation shows X.B ESCALATE preserved against all Magotti 2009 anchors; sampling expansion alone cannot close z > 2 gap; retained for post-Track-A/B/C completion
- [ ] **Fix 4 5-rep on 4 remaining ncAA variants (§A.9, target-selective σ validation, R12)** **P1** (v0.5 NEW) — targets: 1EBP_MTR13, 1YCR_NML22, 3IOL_NML20, 2QKH_MTR25 (skip if calibration failed); reuse existing Option Y DCDs if available (~8h GPU); success: ≥2 of 4 show <1.5× σ reduction → R12 generalized OR all ≥3× → Cp4 was outlier
- [ ] **Option X trial 1YCR_NML22 (MD 5→20ns, 2 seeds, 72h GPU)** **P1** — requires explicit user approval (12h threshold exceeded); σ_between plateau empirically confirmed (n=3: 11.86 → n=5: 10.38, Δσ = −1.47 marginal); acceptance: σ_within drop factor ≈ √4 = 2×, σ_between stability
- [ ] **1EBP_MTR13 literature-unknown benchmark collaboration proposal** **P2** (external) — UPDD null prediction (ΔΔG = +3.63, z=0.70) → experimental Kd titration with ±5× bracket to test null; potential first-of-kind UPDD-first prediction for ncAA
- [ ] **Gate 5 Kendall τ re-computation post-Fix-4-ΔΔG** **P1** — 5 WT ranking stability + 2QKI_Cp4 ΔΔG pair (1-traj vs 3-traj) τ ≥ 0.9; design-level invariance for v0.7 entry
- [ ] **Gate 4 literature_distance anchor switched from −3.96 to Magotti 2009 SSOT pair** (v0.6 NEW) — primary: 2QKI_Cp4_primary −1.4 kcal/mol (ITC Table 3, 4W pure W4→1MeW, DOI: 10.1002/jmr.972); extended: 2QKI_Cp4_extended −3.0 kcal/mol (SPR Table 2, 4V9H→Cp4 evolution); retired registry for [−3.96, −7.8] preserved; `literature_ssot` references `verdict_cp4_literature_comparability_20260423.md`
- [x] **Track A: 4W reference calibration (direct Magotti 2009 comparison)** — ✅ **Completed 2026-04-23** — 4W MD 5-rep + Fix 4 MM-GBSA; ⟨⟨ΔG⟩⟩ = −62.85, σ_btwn = 5.00; 4W−4W9A delta = −0.785 kcal/mol → **H9A methodology gap negligible** (Rev. 6 Khoury expected ~0.3 kcal/mol; falsifies §K.6 Rank 6 hypothesis as primary X.B cause); aggregate: `outputs/analysis/track_a_4w_fix4_20260423.json`
  - → `UPGRADE_v07_roadmap.md` §"Reference Library Mapping / Phase α references" (4W baseline preserved as validation anchor)
- [x] **Track B: MM-PBSA Opt γ parallel track** — ✅ **FINAL 2026-04-24 00:02 KST** — Cp4 PBSA n=5 + WT PBSA n=5 complete; Cp4 ⟨⟨ΔG⟩⟩ = −9.04 (σ_btwn=7.12); WT ⟨⟨ΔG⟩⟩ = −12.96 (σ_btwn=5.77); **ΔΔG(PBSA) = +3.92 kcal/mol**, σ_ΔΔG = 9.16, CI95 = [−7.46, +15.30] (includes Magotti −1.4 AND −3.0); **z_vs_lit(−1.4) = 0.58** → **Tier X.A UNCONDITIONAL** 🟢; 78% reduction vs GB (+17.89); single-system n=5 caveat (σ_btwn=7.12 wide, s42 +1.68 outlier-ish)
  - → `UPGRADE_v07_roadmap.md` §"Phase α" (PBSA flip primary operational change) and §"Phase β" (multi-system validation required)
- [ ] **Track C: Option X MD 5→20ns for Cp4+WT (72h GPU exclusive window)** **reassess post Phase β** (status updated 2026-04-24) — originally addresses §K.6 root cause 2 (MD length insufficient); with Track B empirical X.A UNCONDITIONAL, MD length hypothesis priority downgraded; reassessment trigger: β exit scenario requires it (if multi-system fails)
  - → `UPGRADE_v07_roadmap.md` §"Phase β Decision Points" (β-D2 mid-phase decision may activate Track C) and §"Phase γ Failure Scenarios / Fallback path A"
- [ ] **Track D: Khoury RESP ff03 charges adoption for MTR/OMW** — **Pathway E γ primary deliverable** (status: P0 hotfix elevation per V+K CATR-CRITICAL 2026-04-23 evening; v0.6.2 or v0.7 context retired in favor of Pathway E γ scope)
  - → `UPGRADE_v07_roadmap.md` §γ-S1 (Khoury OMW XML production deploy)
  - → `UPGRADE_v07_roadmap.md` §γ-S2 (Full ncAA re-parameterize, Khoury-compatible subset: MTR/OMW, OMY, NMT, AIB, PAL, ALC, NMW, ALN, NMD, NMG, NMQ, NMR, NMH, OEY, NMA, NMY, NMC)
  - → `UPGRADE_v07_roadmap.md` §γ-S3 (Bonded parameter override validation)
  - → `UPGRADE_v07_roadmap.md` §γ-S4 (Cp4 Fix 4 n=5 under Khoury RESP + PBSA)
  - **Context**: UPDD MTR Σq = −0.1867 e confirmed (V 4-axis verification), 31/119 MM-GBSA runs affected (K audit 100% ncAA XML coverage). Khoury OMW (ff03 RESP, Σq=0.0, same 2QKI system) ready from `Reference/Khoury 2014/sb400168u_si_002/ffncaa.in` line 778-810.

### Pre-entry items — Phase Mapping (Rev. 3.1 added)

Existing Task tracker items mapped to Pathway E phases:

- **Pre-entry 1** (Task #13): Fix `parameterize_ncaa.py` MTR XML rename bug
  - → `UPGRADE_v07_roadmap.md` §"Existing UPGRADE.md Items — Phase Mapping / Phase α scope" (**included in α patch regenerate** — resolved as byproduct of α-S3 bulk regenerate)
- **Pre-entry 2** (Task #14): Compstatin ladder (4W9A/Cp20/4(1MeW)/Cp40) benchmark
  - → `UPGRADE_v07_roadmap.md` §"Phase β scope" (**β multi-system** via Phase 8 Compstatin Pilot in this file, now phase-gated by Pathway E)
- **Pre-entry 3** (Task #15): MD replicate statistical baseline
  - → `UPGRADE_v07_roadmap.md` §"Phase β scope" (**β n-scaling validation**); also this file §P1.4 "4 remaining WT Fix 4 5-rep" Gate 5 full-protocol τ=1.0 proof
- **Pre-entry 4** (Task #16): Cross-class benchmark (enzyme-inhibitor cyclic)
  - → `UPGRADE_v07_roadmap.md` §"Phase β scope" (**β multi-class** — 1SFI bicyclic β-trypsin variant included in β test matrix β-D1 Extended option)

### §A MM-GBSA 1-Trajectory Protocol
- [ ] **Day 1 Manual test (1h GPU, 큐 완료 후)**: 1-traj ΔG 재계산 on 기존 minimized complex
  - [ ] 2QKI w4_1_s2/snap01 (baseline sanity)
  - [ ] 6WGN w1_0_s1 (KRAS Case A/B/C decision)
  - [ ] **추가 (SciVal 권고)**: 1SFI bi-cyclic pilot (8-12h, Zn/Ca labile coordination)
  - [ ] **추가 (SciVal 권고)**: 2QKI Cp4 inter-replicate σ on 3 replicates (τ convergence 증거)
- [ ] **Fix 4 구현** (Day 2-3, Case A 가정):
    - [ ] `utils/run_mmgbsa.py::split_complex_1traj` (기존 `split_complex` + serial→bucket index sets)
    - [ ] `utils/run_mmgbsa.py::minimize_complex_and_extract` (기존 `calc_energy` 의 최소화 부분 분리)
    - [ ] `utils/run_mmgbsa.py::calc_energy_static` (no-minimize PE, `NoCutoff` + GBn2)
    - [ ] `utils/run_mmgbsa.py::calc_mmgbsa_1traj` — dispatcher
    - [ ] Env var `UPDD_MMGBSA_PROTOCOL={1traj|3traj}` (default `3traj` until Go verdict)
    - [ ] 예상 LOC: ~220 (180 new + 40 modified)
- [ ] **Schema v0.6.8 bump** (SciVal 권고 field 추가):
    - [ ] `mm_gbsa_protocol` (`1traj` / `3traj`)
    - [ ] `min_iter_complex`, `min_iter_receptor=0`, `min_iter_ligand=0`
    - [ ] `entropy_method` (`IE` / `NMA` / `none`)
    - [ ] `historical_compat` (bool, F2/F1/P1-era 구분)
    - [ ] `gb_model` (`igb=5` OBC2 or `igb=8` GBn2)
    - [ ] `salt_concentration_M` (default 0.15)
    - [ ] `sasa_model`, `surften_kcal_mol_A2` (nonpolar SASA term)
    - [ ] `interaction_entropy_window_ns` (IE averaging)
    - [ ] `receptor_ligand_extraction_mode` (`strict_static` / `zero_step_min`)
    - [ ] `switch_rmsd_gate` ({enabled, regions, threshold_median, threshold_max})
    - [ ] `replicates_id`, `replicate_aggregation` (`mean` / `median` / `block_bootstrap`)
    - [ ] `ff_cofactor_params_hash` (R-17 FF drift 방지)
    - [ ] `protocol_doi` / `protocol_hash`
- [ ] **12-snap regression test** (Day 4): 1-traj vs 3-traj 결과 비교 테이블
- [ ] **Unit tests**: OpenMM CPU platform 으로 queue 진행 중에도 테스트 가능

### §B Switch RMSD Monitoring Module
- [ ] 신규 `utils/trajectory_validation.py` (~160 LOC, MDTraj 의존)
    - [ ] `monitor_switch_rmsd(trajectory, regions)` per-frame Cα RMSD
    - [ ] Tier: safe (<2.0 Å) / moderate (2.0-2.5) / unstable (>2.5 max) / warn_max (>3.5)
    - [ ] SciVal 권고: **binding-pocket heavy-atom RMSD** (5 Å shell around ligand) 추가 — 비-GTPase target 에도 적용
- [ ] 호출 지점:
    - [ ] Post-MD `run_restrained_md.py::run_md()` tail → `<design>_switch_rmsd.json`
    - [ ] CLI standalone for post-hoc audit
    - [ ] `run_mmgbsa.py` 는 advisory only (queue 안전 위해 hard gate 금지, v0.6.8)
- [ ] Target-specific regions (target_card field `switch_region_monitoring`):
    - [ ] 6WGN KRAS: Switch I (32-38), Switch II (59-67), P-loop (10-17)
    - [ ] Generalization: 모든 GTPase target (RhoA, Rac1 v0.8+)

### §C Cofactor RMSD Threshold
- [ ] `utils/reinsert_cofactors_post_af2.py::DEFAULT_RMSD_THRESHOLD_A` 환경변수 opt-in
    - Current: 2.5 Å hardcoded (**UPDDfuture §C.1 "2.0 → 2.5" 서술 outdated — 이미 2.5 Å**)
    - Change: `float(os.environ.get("UPDD_G2_RMSD_THRESHOLD_A", "2.5"))`
- [ ] Phase 2 (v0.7): default 2.5 유지, tightening 은 empirical test 이후
- [ ] **SciVal 권고 추가**: §C.4 geometry sanity checks (bond-length/angle/clash) 를 G1 gate 에 추가

### §D Phase 3 Stage 3 Historical Audit
- [ ] 신규 `utils/phase3_historical_audit.py` (~240 LOC)
    - [ ] `audit_historical_results(output_root, archive_root)` 스캐너
    - [ ] `classify_result(json_path, git_blame)` — 6-bucket {PASS, PARITY_ONLY, F2_AFFECTED, F1_AFFECTED, P1_AFFECTED, FAIL}
    - [ ] Reuse `utils/audit_charge_consistency.py::audit_json_charges`
    - [ ] Git commit-era classifier (timestamp → F2/F1/P1 fix 시기)
    - [ ] `generate_rerun_plan(audit) -> dict`
    - [ ] CLI: `--output-root`, `--dry-run` (default), `--move-archive`
- [ ] `outputs/_archive/` 재구조화:
    - `v0.5_pre_R15_R16/`
    - `v0.6_wrong_charge/` (existing)
    - `v0.6_pre_1traj/` (Fix 4 커밋 시 신규)
- [ ] 12-16h rerun (PASS → no action; PARITY_ONLY → schema bump; FAIL → rerun)
- [ ] **v0.7 entry blocking**

---

## v0.7 — WT Counterfactual + Calibration (진입 조건 Rev. 3 강화)

### v0.7 진입 조건 (SciVal CONDITIONAL APPROVE 반영, Rev. 3)

YAML 형식 acceptance gates (UPDDfuture §G + SciVal Verdict 반영):

```yaml
v0.7_entry_criteria:
  # Energy magnitude (SciVal: 80→60 tighten, 80 warn)
  abs_dG_median_kcal_mol:
    block: 80  # REJECT
    warn:  60  # review
  abs_E_ligand_stdev_kcal_mol:
    max: 80
  design_stdev_ratio:
    max: 0.3   # σ/|μ|

  # Ranking validation (SciVal: tau 0.9 on N=4 임정의 수치적 문제 → Spearman ρ)
  kendall_tau_4pt_compstatin:
    min: 0.667  # ≤1 inversion out of 6 pairs
    requires: "ensemble-averaged ΔG ≥3 MD replicates"
  spearman_rho_5pt_ladder:
    min: 0.8
    note: "Cp20/30/40 extension 권고 — 4-point 에서 statistical 제약"
  correlation_r2:
    min: 0.5   # predicted vs log10(Kd)

  # Convergence
  converged_snaps_ratio:
    min: 0.8   # 10-snap 중 ≥8
  replicates: {min: 3}
  snapshots_per_replicate: {min: 10}

  # GTPase-specific (SciVal 🟢)
  switch_rmsd_median_A:
    max: 2.5
  switch_rmsd_max_A:
    max: 3.5   # SciVal 권고 (transient excursion catch)
  applies_to: ["6WGN", "RhoA targets", ...]

  # Phase 3 Stage 3
  phase3_audit_status:
    required: complete
  charge_consistency_audit:
    required: passed (모든 snap)

  # Protocol provenance (v0.6.8 schema)
  mm_gbsa_protocol: {required_value: "1traj"}
  ff_cofactor_params_hash: {required: frozen}

  # Failure action per gate
  failure_action:
    block_entry: [abs_dG_median, kendall_tau, phase3_audit, mm_gbsa_protocol]
    investigate: [correlation_r2, converged_snaps_ratio]
    rerun: [charge_consistency_audit]
    warn_and_fallback: [abs_dG_median_warn, switch_rmsd_median]
```

### Phase 7: WT Counterfactual
- [ ] `variants/{ncaa,wt}/` 디렉토리 구조
- [ ] `assess_wt_comparability()` 3-tier 판정
- [x] 1-Me-Trp (MTR) registry + generic parent walker — v0.6.7 (`dc46690`)
- [x] **Rev. 3 추가**: 24 추가 ncAA registry parent_residue (NML/NMV/... + phospho SEP/TPO/PTR + D-aa + 할로겐)

### Phase 8: Compstatin Pilot (Rev. 3 확장)
- [x] PDB 2QKI target 준비
- [x] **Rev. 3 2-point**: WT (4W9A) + Cp4 (W4→MTR4) MM-GBSA 검증 완료 (2ns × 5 snap; 5ns × 10 snap 재측정 진행 중)
- [ ] **Rev. 3 확장**: Cp20 (W4→MTR + 추가 변이, literature clarification 필요)
- [ ] **Rev. 3 확장**: Cp40 (DTY + SAR + NML + MTR on 7BEB, `prepare_cp40.py` 신규 스크립트 필요)
- [ ] Calibration fit (4 points)
- [ ] `calibration/v0.7.0-pilot-compstatin.json`

### Phase 9: Sa-D3 Pilot
- [x] PDB 7TL8 target 준비 + W6→MTR MD 검증
- [ ] 2 variants: N-Me, ΔN-Me (Rev. 2 원안 유지)
- [ ] Calibration 확장 (총 5 points)
- [ ] `calibration/v0.7.0-pilot-sad3.json`

### Phase 10 (Rev. 3 신규): Class B GPCR Pilot
- [x] PDB 3IOL (GLP-1R + GLP-1) + L20→NML + Y19→PTR MD 검증
- [x] PDB 2QKH (GIPR + GIP) + W25→MTR MD 검증
- [ ] Glucagon (5EE7) full-structure prep (MBP fusion + chain split 복잡성)
- [ ] CRF1R (3EHU) full-structure prep

### Phase 11 (Rev. 3 신규): Cross-class Benchmark (1SFI)
- [x] 1SFI bicyclic topology infrastructure (htc+ss 양 branch) — v0.6.6
- [x] WT 1SFI MD σ=8.28 검증
- [ ] NMK5 (Lys5→N-Me-Lys) variant MD
- [ ] NMR2 (Arg2→N-Me-Arg) variant MD

### Post-hoc Calibration Layer
- [ ] `utils/calibrate_results.py` 모듈
- [ ] `regime: ranking_only` → `regime: absolute_calibrated` 전환
- [ ] Config version 관리 (`calibration/v0.7.x-*.json`)

---

## v0.8 — Calibration Expansion + Publication Track

- [ ] **15-IgBP pilot** (Rev. 1 에서 연기)
- [ ] **MD replicates 확장** (3x → 9x short-trajectory, Hou 2011 권고)
- [ ] **Bootstrap CI** (Genheden 2015 권고)
- [ ] **Paper draft**
- [ ] **Rev. 3 추가**: **JoltQC Acceleration Stage 2** (opt-in hybrid mode)
    - [ ] Axis 3 DEFINITIVE PASS (Stage 7 re-run after correct-charge)
    - [ ] `UPDD_QMMM_BACKEND={pyscf|joltqc}` opt-in env var
    - [ ] Screening task default JoltQC, calibration default FP64
- [ ] **Rev. 3 신규 후보**: Option γ **Explicit MM-PBSA** (UPDDfuture §L)
    - [ ] Condition: 1-traj 로도 KRAS unresolved (§K Scenario B/C)
    - [ ] `utils/run_mmpbsa.py` (AMBER MMPBSA.py wrapper)
    - [ ] Cost: 10-20× MM-GBSA, implementation 1-2 week + validation 2 week

---

## v0.9+ — Future Tracks

### ncAA Registry 확장 (11 deferred, UPDDfuture §M)
- [ ] Silicon TMS/SIP/DPS (GAFF2 Si params + RESP, 2 week)
- [ ] Chain-length variants ORN/DAB/HAR/NLE (1 week)
- [ ] Cα-quaternary AIB (2 week, backbone dihedral 재검토)
- [ ] Click chemistry PRA/AHA (1 week)
- [ ] Fused-ring TIC (custom parent definition)

### Multi-MD Sampling (UPDDfuture §N)
- [ ] 9-replicate × 1.5 ns parallelism (Hou 2011 권고)
- [ ] Adaptive sampling / metadynamics (binding-mode multiplicity)
- [ ] Snap count 10 → 25+ per replicate

### Snap-level Runaway Filter (UPDDfuture §E)
- [ ] `utils/snap_quality.py::compute_interface_contact_count`
- [ ] `snap_level_runaway_check` (contact < 30% of initial → drop snap, not design)
- [ ] `rejected_snaps.json` / `rejected_designs.json` 분리

### JoltQC Full Integration (UPDDfuture §H, v1.0)
- [ ] Default acceleration after 2+ weeks stable hybrid mode

---

## Completed ✅ (Rev. 3 누적)

### v0.7 Pathway E adoption (2026-04-24 morning) — Rev. 3.1

- [x] **Pathway E decision** (User 2026-04-24 morning) — α/β/γ evolutionary phased approach; reframes X3-B single-verdict into 4-pathway option landscape with P-flavored α + S-flavored β + C-flavored γ
- [x] **UPGRADE_v07_roadmap.md creation** (545 lines, 2026-04-24) — Pathway E authoritative roadmap: Overview, Data Flow, Reference Library Mapping, Paper 1 Narrative Evolution, Phase α/β/γ entry/exit/decision criteria, Failure Scenarios & Rollback, Decision Audit Trail
- [x] **Track B Cp4 PBSA FINAL 5/5** (2026-04-24 00:02 KST) — Cp4 Fix 4 PBSA n=5 aggregation: ⟨⟨ΔG⟩⟩=−9.04, σ_btwn=7.12; ΔΔG(PBSA)=+3.92 kcal/mol, CI95=[−7.46,+15.30], z_vs_lit(−1.4)=0.58 → **X.A UNCONDITIONAL** (PBSA-regime, n=5 caveat); 78% reduction vs GB
- [x] **Archi Light Synthesis v3** (2026-04-24 00:55 KST) — 426 lines, 11 pages, 3-axis integration (X.B resolution / pipeline integrity / future strategy), 4-pathway comparison (P / B-revised / S / A-only), 12 [USER-DECIDE] questions, flags [X3-B-CHALLENGED-BY-TRACK-B] + [REFERENCE-TRACK-B-TENSION] + [MULTI-AXIS-PARADOX]; output: `outputs/synthesis/v07_light_synthesis_20260424.md`
- [x] **UPGRADE.md Rev. 3 → 3.1 dual-file protocol** (2026-04-24) — this entry; header dual-file note added; v0.7 Entry Prep cross-references to `UPGRADE_v07_roadmap.md` integrated; Track A/B/C/D status resolved + phase-mapped; Pre-entry 1-4 items phase-mapped; conflict resolution protocol documented; archive backup: `archive/UPGRADE_rev3_pre_dual_file_20260424.md`

### X.B investigation (2026-04-23 evening) — CATR-CRITICAL dual finding

- [x] **Σq V verification** (2026-04-23 evening) — 4-axis independent triangulation (Method 3 XML parse + Method 1 parmed prmtop + V2 10-seed determinism + V4 4-ncAA cross-check); **UPDD MTR Σq = −0.186700 e CONFIRMED** (10/10 seeds identical, zero variance); backbone N deviation −0.478 e vs amb14SB; 4/4 ncAAs non-neutral (MTR −0.187, MLE +0.165, NMA +0.160, PTR −1.80/−2.17 version drift); standard TRP control Σq=0.0000 confirms amb14SB patch applied to standard aa but NOT to ncAA path; verdict: `outputs/analysis/sigma_q_verification_20260423.md`
- [x] **amb14SB K audit** (2026-04-23 evening) — CATR-CRITICAL; `parameterize_ncaa.py` L1469 env-var gate `UPDD_MTR_AMBER14_PATCH=1` default OFF (no production orchestrator sets it); 24/24 production `params/MTR_gaff2.xml` unpatched; 25/27 `_md_input/MTR_gaff2.xml` unpatched (2 smoke-test patched); 31/119 MM-GBSA runs (26%, all MTR-variant) carry unpatched ncAA; X.B ΔΔG charge contribution estimated 5-15 kcal/mol (30-60% of +17.89 gap); Recovery Option A (default flip, 1-2h script + 15-20h re-runs) vs Option C (Khoury RESP, ~1.5 days); new gate PATCH-01 proposed; verdict: `outputs/audit/amb14sb_patch_audit_20260423.md`
- [x] **Reference Library construction** (2026-04-23 afternoon, user) — 7 papers downloaded to `/home/san/UPDD_proj/Reference/`: Khoury 2014 (ACS Synth Biol), Mulakala 2013 (JMGM), Weng-Hou 2019 (PCCP), Hasan 2025 (JPCB), Forouzesh 2021 (Molecules), Maffucci-Contini 2016 (JCIM), Kilburg-Gallicchio 2018 (Front Mol Biosci); supporting SI including Khoury 2014 SI sb400168u_si_002/ffncaa.in (147 frcmods + OMW block line 778-810)
- [x] **Khoury 2014 Coder Phase 1 comparison** (2026-04-23 evening) — UPDD MTR ≡ Khoury OMW (chemistry aligned, 27 atoms, Σq=0.0); charge comparison RMSD 0.163 e, max |Δq| 0.6129 e on backbone N; 10 atoms |Δq|>0.1e; 7 heavy-atom sign flips; Track D strong recommendation; report: `Reference/comparison_mtr_omw_20260423.md`
- [x] **Reference Deep Research 7 papers** (2026-04-23 evening, single dispatch ~4h) — DEPTH mode per-paper Section 1-7 analysis; comparative matrix; verdict **X3-B: A+C+W (Nwat=30 stack)** confidence MEDIUM-HIGH; per-paper: `outputs/analysis/reference_deep_research/paper_{1..7}_*.md`; integrated: `outputs/analysis/reference_deep_research_20260423.md`
  - **Status 2026-04-24**: X3-B reframed as one of 4 Pathway options per Archi Light v3 + Track B FINAL; X3-B preserved but not exclusive answer; see `UPGRADE_v07_roadmap.md` for Pathway E synthesis

### v0.6 Literature re-anchor (2026-04-23)
- [x] **SciVal literature re-verification 2026-04-23** (Magotti 2009 SSOT, retired −3.96) — `verdict_cp4_literature_comparability_20260423.md`; Janssen DOI mis-attribution exposed (authoritative 2QKI paper = Janssen 2007 JBC 10.1074/jbc.M704587200, structural only); 12 μM 4W9A Kd input unverified; new SSOT pair: −1.4 (pure, ITC Table 3, DOI: 10.1002/jmr.972) / −3.0 (extended, SPR Table 2)
- [x] **Partner independent Magotti 2009 extraction** (2026-04-23) — Table 3 ITC (4W Kd=140 nM, Cp4 Kd=15 nM) + Table 2 SPR (4V9H Kd=1641 nM, Cp4 Kd=11 nM) independently reproduced by user
- [x] **Fix 4 symmetric n=10 (WT) + n=8 (Cp4) aggregation** (`outputs/analysis/cp4_fix4_symmetric_n10_20260423.json`) — ΔΔG = +17.89, σ_ΔΔG = 9.77, 95% CI [+9.72, +26.06]; z_vs_lit 2.41 (−1.4) / 2.61 (−3.0) / 2.24 (sign gate) → X.B ESCALATE across all live anchors
- [x] **UPDDfuture.md v0.5 → v0.6 synthesis** (Archi, 2026-04-23) — §v0.6 Changelog + §G.1 Gate 4 dual Magotti anchors + retired registry + §K.6 root cause 4→6 + §Q R12 UPDATE + R13/R14 NEW + §F probability refinement + §9 v0.6 version history
- [x] **v0.6.1 refinement synthesis** (regime-split + hypothesis re-rank, 2026-04-23) — Archi Partner refinement integration: §F.1/F.2/F.3 regime-split (WT-absolute 70-80%/85-92%, ΔΔG 15-25%/25-35%, blended deprecated), §K.6 root cause H9A-adjusted re-rank (H1 FF TOP, H6 4W9A DOWN to 1.5%) + testability matrix, §G Gate 4 paradox explicit report_only reaffirmation + paper discussion draft; no new science, no code changes

### Pre-2026-04-21
- [x] v0.5 SSOT refactor (2026-04-17)
- [x] v0.5 Bug 1-3 fix (2026-04-17)
- [x] Path agent 도입 (2026-04-19 Rev. 1)
- [x] v0.6.0 QM/MM physical validity (2026-04-18 `14b0295`)
- [x] SCF bundle B2/B3/B5/B6/B8/B9 + DF auto-mode (2026-04-19)
- [x] R-11 regime field in result JSON (2026-04-19)
- [x] CLAUDE.md Rev. 2 반영 (2026-04-19)
- [x] UPGRADE.md Rev. 2 반영 (2026-04-19)

### v0.6.5 → v0.6.6 (2026-04-19 ~ 20)
- [x] **R-15/R-16 Stage 1** merge `4034ee6` (charge_topology + audit_charge_consistency + 23/23 tests)
- [x] **R-17 Stage 4** (cofactor preservation 6-gate + error classes + 15 tests)
- [x] CLAUDE.md Rev. 3 (JoltQC track + Rev. 3 Phase 3 Stage 1 반영)
- [x] CLAUDE.md Rev. 4 (Keeper agent 신설, 6-Agent taxonomy)
- [x] Schema v0.6.5 (cofactor_residues full fields)
- [x] F1 Mg²⁺ strip 복구 + F2 cyclic HTC proper FF + P1 5000-iter minimize
- [x] Ultrareview v22 P0 silent bugs 3-item fix (CustomGBForce param ordering 등)

### Agent Rev. 5 (2026-04-22) — Post-incident 6-agent hardening
- [x] **Backup**: `archive/claude_backups/pre_rev5_20260422/` (6 original agent files preserved)
- [x] **scival.md** (276 → 361): SIV + AHCD + PKI (Cp4 = −4.0 kcal/mol with DOI) + SRE thresholds + Cross-Agent Triangulation + Output Format 확장
- [x] **archi-architect.md** (221 → 272): SIV + AHCD (verified-before-recommend) + PKI + CATR (3-way code deploy, pipeline launch, strategic pivot)
- [x] **coder.md** (252 → 301): SIV + AHCD + F1/F2/P1 auto-BLOCK 금지 패턴 + Cross-Agent gates
- [x] **keeper.md** (255 → 303): SIV (ironic) + AHCD + **AGENT-01/02/03 self-surveillance gates** (agent file presence/size, memory dir, `.git/info/exclude`)
- [x] **path.md** (376 → 433): SIV + AHCD + 2026-04-22 statistical lesson (z=0.53 sampling inadequacy) + Sign confidence / Sampling threshold tables + Cross-Agent Escalation
- [x] **pipeline-runner.md** (269 → 329): SIV + AHCD + PARTIAL_SUCCESS rm-prohibition + forensic abort-state preservation + 6-step canonical launch sequence
- [x] Grep integrity check: 41 Rev. 5 signature matches across 6 files; `−7.8` appears only in FORBIDDEN / incident narrative contexts (no live claim)
- [x] UPDATE.md [Agent Rev. 5] entry 작성 (R-8)
- **Trigger**: 2026-04-14 `scival.md` silent deletion → Cp4 ΔΔG = −7.8 hallucination (2026-04-22 감지) → Forward Plan Rev 0.1 반려 → Rev 0.2 + Agent hardening
- **Source**: `AGENT_IMPROVEMENT_PLAN.md` (Rev. 5, 779 lines)
- **Verification (next session)**: Keeper AGENT-01/02/03 run; SciVal Cp4 ΔΔG self-test → −4.0 kcal/mol

### v0.6.7 (2026-04-21)
- [x] `dc46690` Generic amber14 parent walker + 14 ncAA parent_residue (LEU/VAL/PHE/ARG/TYR/SER/THR/LYS/GLN/GLY/ALA/PRO)
- [x] `05944fc` Phospho multi-atom chained extensions (SEP/TPO/PTR P+3O) + halogens (CLP/BRP) + D-aa extension (DPN/DAL)
- [x] `584d3e0` Chained-extension mutate geometry (tetrahedral) + unique phospho H naming (HP1/H1P/H2P collision 해결)
- [x] `7075416` Pass 1b seed-variant retry (auto recovery from probabilistic NaN)
- [x] `9c3f8a5` External peptide-bond CONECT emission for mid-chain ncAA + DTY D-Tyrosine + SAR→GLY parent migration
- [x] `a4c1047` MM-GBSA split_complex + write_temp_pdb CONECT bucket-wise preservation
- [x] `1333e17` Snapshot ncAA CONECT auto-injection (`_inject_ncaa_conect`: heavy bonds + parent H table + peptide bonds)
- [x] `f8fe6f9` Unconditional spurious-H strip + snap atom-whitelist filtering (MTR H9 artifact 해결 → Cp4 MM-GBSA 성공)
- [x] **8 target × ncAA pipeline** 검증: 2QKI/1EBP/7TL8/6WGN/1YCR/3IOL/2QKH/1SFI
- [x] **Class B GPCR 확장**: 3IOL (GLP-1R) + 2QKH (GIPR)
- [x] **MTR NVT NaN root cause** 해결: Stage 1 min 후 `k=100` reset universal
- [x] **Cp4 MM-GBSA 완료**: 2QKI W4→MTR ⟨ΔG⟩ = −13.63 kcal/mol (2ns × 5 snap, seed=19)
- [x] **Calibration queue launch** (6-target × 3-rep × 5ns × 10-snap, PID 1955081, ETA ~24-30h)
    - Run 1 완료: 2QKI_WT seed=55 ⟨ΔG⟩ = −47.75 ± 17.44 (n=10, σ/|μ|=0.37)
- [x] **UPDDnow.md 작성** (project root, 12 sections)
- [x] **UPDDfuture.md rev 0.2** 작성 (§A-R, 1297 lines)
- [x] **CLAUDE.md Rev. 4** Agent taxonomy 유지
- [x] **UPGRADE.md Rev. 3 반영** (본 문서)

---

## Hardware-Aware Design (Rev. 3 유지)

### 현재 (RTX 5070 Ti 16GB + 32GB RAM + PCIe 5.0 x16)
- Direct-SCF only (DF auto → off)
- n_qm ≤ 600 (Rev. 3 완화)
- 1 snap 당 30-60 min, 25 snap 전체 12-25 hours
- Pre-QM filter 로 탈락 design 조기 식별
- MD 5ns × 2.5 GB DCD 500ms/step (CUDA)

### RTX 5090 32GB 업그레이드 시 (자동 전환)
- DF auto → on (threshold 20 GB)
- n_qm ~ 700 까지 가능
- 1 snap 당 10-15 min, 25 snap 전체 4-6 hours
- 코드 수정 불필요 (env var 기반 자동 판정)

### H100 80GB (calibration 확장 단계)
- DF 안정 + replicate 확장 (9-rep)
- n_qm 1000+ 가능
- v0.8 publication track 의 하드웨어 기반
- JoltQC default 전환 이후 H100 활용도 2× 이상

---

## Coordination Call-outs (Rev. 3 신규, Archi consult 반영)

### SciVal Verdict Gap
- [ ] **1-trajectory protocol migration verdict** (§A, 🟡 CONDITIONAL APPROVE 이미 존재 but 공식 verdict 파일 생성 필요)
- [ ] **Switch RMSD threshold** verdict (§B, 2.0 safe / 2.5 moderate literature-derived)
- [ ] **4-point ladder τ gate 재규정** (SciVal: τ=0.9 on N=4 → impossible, τ=1.0 required, 권고 ρ≥0.8 on 5-point)
- [ ] **1SFI bi-cyclic 1-traj 타당성** (SciVal: 🔴 untested, 추가 pilot 필요)

### Coder Task Queue (GPU-free parallel during calibration queue)
- [ ] Day 1-2: §C env var (trivial) + §B `trajectory_validation.py` scaffold + §D audit composer
- [ ] Day 3-4: §A helpers (`split_complex_1traj`, `calc_energy_static`) on CPU platform
- [ ] Day 5-6 (큐 완료 후): §A.4 manual test GPU 1h + Go/No-Go + Fix 4 ship
- [ ] Day 7: §D classification run + SciVal/Keeper 통합 verdict

### Keeper Integrity Extension
- [ ] `integrity_1traj_coord_extraction.md` 신설 (receptor_from_min 추출 이 실제 minimized complex frame 과 일치하는지)
- [ ] Pre-flight gate: §A Fix 4 구현 후 **모든 snap** 에 대해 extraction-frame check

### Calibration Queue Disposition
- [ ] **현재 진행 중 queue** (3-traj protocol, 5000-iter) = pre-§A baseline
- [ ] §A Case A Go 시 → queue **finish-and-rerun** 결정 (schema v0.6.8 `historical_compat=false` 로 격리)
- [ ] §A Case B/C 시 → queue 결과 사용 + KRAS 별도 track

### Known Documentation Drift
- [ ] UPDDfuture §C "2.0 → 2.5" 설명 outdated (이미 2.5 hardcoded)
- [ ] UPDDfuture §F Cp4 "2ns" 기록은 Cp4_production dir 의 pilot MD; calibration queue 는 5ns × 10-snap
- [ ] UPDDfuture §M "11 defer" 목록에서 SAR 은 이미 parent_residue=GLY 추가되어 실제 defer 수는 10

---

**이 roadmap 은 진실의 원천이다. 각 Phase 완료 시 Archi 가 즉시 업데이트한다.**

**Rev. 3 다음 업데이트 트리거**:
- Calibration queue 완료 (6-target × 3-rep × 5ns × 10-snap aggregation 결과)
- §A.4 Manual 1-traj test Case A/B/C 결정
- Phase 3 Stage 3 audit 실행 결과 (PASS/PARITY_ONLY/FAIL 분포)
