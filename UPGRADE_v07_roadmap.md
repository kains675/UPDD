# v0.7 Roadmap — Pathway E (Evolutionary Phased Approach)

**Defined**: 2026-04-24 (post Track B PBSA final + V+K CATR-CRITICAL + Archi Light v3)
**Strategy**: Incremental phased improvement (P → β → γ)
**Source**: User decision 2026-04-24 morning + Archi Light v3 4-pathway analysis
**Supersedes**: Prior X3-B (A+C+W) single-path verdict (reframed)

---

## Overview

Pathway E addresses X.B resolution via incremental phased approach, leveraging Track B PBSA empirical breakthrough while preserving rigor path to root-cause fix:

- **Phase α** (~1 week): Immediate resolution via PBSA flip + charge fix (P-flavored)
- **Phase β** (~2-4 weeks): Multi-system validation + Paper 1 preprint (S-flavored)
- **Phase γ** (v0.8 or on-trigger): Khoury ff03 RESP root-cause fix (C-flavored)

Philosophy: **"Make it work, make it right, make it fast"**
- α: Make it work (empirical X.B resolution for Cp4)
- β: Make it right (multi-system validation, Paper 1 credibility)
- γ: Make it fast (literature-aligned, future-extensible)

---

## Data Flow Between Phases

```
Phase α output → Phase β input:
  - 31 patched ncAA XMLs (all Σq = 0.000 ± 1e-6)
  - PATCH-01 Keeper gate module
  - NVT smoke test results (4 ncAA classes validated)
  - Cp4 Fix 4 PBSA post-patch n=5 results
  - v0.6.x archive manifest

Phase β output → Phase γ input:
  - Multi-system PBSA evidence (3 systems)
  - Cross-system consistency metrics + decomposition analysis
  - Khoury OMW XML pipeline (validated, not yet applied)
  - Paper 1 preprint feedback (reviewer comments if available)
  - γ necessity assessment (LOW/MED/HIGH)

Phase γ output → v0.8 release:
  - Khoury charges applied (production)
  - FF mixing validated at scale
  - Multi-system re-verified under Khoury RESP
  - Paper 1 Journal-ready manuscript
  - KHOURY-01 Keeper gate integrated
```

Each phase's output is the next phase's gated input. Failures at any phase trigger rollback/fallback (see Phase-specific sections).

---

## Reference Library Mapping

Which Reference papers inform each phase:

**Phase α references** (`Reference/INDEX.md`):
- **Khoury 2014 Deep Dive** — patch table design, OMW atom topology
- V verification report (`outputs/analysis/sigma_q_verification_20260423.md`)
- K audit report (`outputs/audit/amb14sb_patch_audit_20260423.md`)

**Phase β references**:
- **Weng-Hou 2019** (C9CP01674K.pdf) — peptide-protein protocol standards
- **Hasan 2025** (MnM-W-MMGBSA.pdf) — Nwat considerations for cross-validation
- **Mulakala 2013** (Could MM-GBSA.pdf) — quantitative accuracy baseline
- **Forouzesh 2021** (Effective MM-GBSA Absolute.pdf) — absolute BFE reference

**Phase γ references**:
- **Khoury 2014 SI** (sb400168u_si_001.pdf, ffncaa.in, OMW.frcmod) — ff03 RESP charges
- **Maffucci-Contini 2016** (Improved Computation PPI.pdf) — Nwat=30 comparison
- **Kilburg-Gallicchio 2018** (Assessment Single Decoupling.pdf) — alternative if γ fails (SDM)

---

## Paper 1 Narrative Evolution

Paper 1 manuscript evolves across phases:

**After Phase α** (Preprint v0.1):
- Title: "Pipeline integrity audit reveals charge bug; PBSA protocol resolves quantitative gap for compstatin Cp4"
- Scope: Single-system (Cp4) resolution + pipeline bug discovery
- Suitable for: Internal documentation, early preprint

**After Phase β success** (Preprint v0.2 / Journal v1):
- Title: "Systematic MM-PBSA validation for cyclic peptide + ncAA systems"
- Scope: Multi-system evidence (3+1 systems), PBSA as default protocol
- Suitable for: bioRxiv, then Journal submission

**After Phase β partial** (Preprint v0.2 revised):
- Title: "System-specific MM-PBSA benefits; Khoury RESP for outlier systems"
- Scope: Mixed evidence, honest limitations discussion
- Suitable for: Preprint with γ plan disclosed

**After Phase γ** (Journal v1 / v2):
- Title: "Literature-aligned ff03 RESP + MM-PBSA for cyclic peptide + ncAA quantitative screening"
- Scope: Full validated methodology
- Suitable for: Journal submission (JCIM, JCTC, JPCB candidates)

---

## Existing UPGRADE.md Items — Phase Mapping

Existing items reorganized into phases:

**Phase α scope**:
- Pre-entry 1: MTR XML rename bug → **included in α patch regenerate**
- v0.6 + Rev. 6 CATR 3-layer ack (user pending) → **α entry prerequisite**

**Phase β scope**:
- Pre-entry 2: Compstatin ladder (4W9A/Cp20/4(1MeW)/Cp40) → **β multi-system**
- Pre-entry 3: MD replicate statistical baseline → **β n-scaling validation**
- Pre-entry 4: Cross-class benchmark (enzyme-inhibitor cyclic) → **β multi-class**
- P1.4 — 4 remaining WT Fix 4 5-rep → **β as needed for cross-validation**

**Phase γ scope**:
- Track D — Khoury RESP charges adopt for MTR/OMW → **γ primary deliverable**

**Phase-independent / hygiene**:
- (Add new items here as discovered during phase execution)

---

## Phase α — Immediate Resolution

**Goal**: X.B resolution + pipeline integrity + PBSA default operational
**Duration**: ~1 week wall clock
**Pathway alignment**: Pathway P (Protocol Flip) with mandatory A (charge fix)

### Phase α Decision Points

- **α-D1 (pre-dispatch)**: SciVal verdict on Option A (amb14SB patch default flip)
  - GO → proceed with α scope
  - CAUTIOUS-GO → smoke tests before full scope
  - NO-GO → escalate to γ-first (hybrid α+γ immediate)
  - Decided by: User + SciVal

- **α-D2 (mid-phase)**: NVT smoke test result across 4 ncAA classes
  - ALL PASS → proceed with bulk regenerate
  - PARTIAL FAIL → per-class investigation, potential γ elements pulled forward
  - ALL FAIL → SciVal emergency, likely γ acceleration
  - Decided by: User + SciVal

- **α-D3 (post-exit)**: β transition timing
  - Immediate β → exit criteria all met, proceed
  - Pause β → α exit met but Paper 1 α-complete preprint desired first
  - Rollback → exit criteria not met, debug
  - Decided by: User

### Phase α Entry Criteria

- [ ] **α-E1**: SciVal verdict on Option A received (α-D1 resolved)
- [ ] **α-E2**: v0.6.x runs archived (non-destructive, original preserved in `outputs/_archive/v0.6.x_charge_bias/` with manifest)
- [ ] **α-E3**: Keeper gate PATCH-01 designed (threshold: backbone N |Δq| < 0.02 e vs parent amb14SB reference)
- [ ] **α-E4**: 31 ncAA patch table extension script ready (extends MTR pattern to NML/MLE/NMA/NMV/NMK/NMR/NMQ/MEA/SAR/DVA/DLE/DTR/DAR/DPR/DPN/DAL/DTY/HYP/CHA/PFF/CLP/BRP/MTR/SEP/TPO/PTR)
- [ ] **α-E5**: UPDATE.md v0.7 entry draft prepared
- [ ] **α-E6**: v0.6.9 commit baseline (git branch `v0.7-alpha` forked from current main)
- [ ] **α-E7**: User ack of Pathway E strategy (this UPGRADE.md section)

### Phase α Scope

- [ ] **α-S1**: amb14SB patch default flip
  - Option A: `UPDD_MTR_AMBER14_PATCH=1` made default in env
  - Option B: remove env-var gate at `parameterize_ncaa.py:1469`, enforce always-on
  - (User decision at α-D1 based on SciVal)

- [ ] **α-S2**: Patch table extension to 31 ncAAs
  - Per-ncAA backbone patch tables
  - Coder dispatch: ~1-2h script + verification

- [ ] **α-S3**: Regenerate all 31 ncAA XMLs
  - Bulk script with Keeper PATCH-01 gate active
  - Σq verification for each
  - Fail-fast on any violation (CofactorParamMissingError pattern)

- [ ] **α-S4**: NVT smoke tests (4 ncAA classes × 100 ps)
  - N-methyl (NMA, MLE, NML)
  - D-amino acid (DVA, DLE, DAR)
  - Side-chain modified (MTR, CHA, HYP)
  - Phospho (SEP, TPO, PTR)
  - Each: 100 ps NVT warmup, energy/force monitoring
  - ~10 min GPU total

- [ ] **α-S5**: Keeper gate PATCH-01 integration
  - Pre-run check: backbone N charge validation
  - Follows R-17 COFACTOR-01 pattern
  - Per-class tolerance table (standard 0.02 e, phospho may need wider)

- [ ] **α-S6**: Cp4 Fix 4 re-run under patched charges + PBSA
  - Same seeds as Fix 4 baseline (s7, s42, s23, s101, s19)
  - PBSA protocol (not GB)
  - n=5 (matching Track B completeness)

- [ ] **α-S7**: Aggregation + tier classification
  - `outputs/analysis/cp4_alpha_pbsa_patched_20260424.json`
  - z_vs_lit computation + CI95
  - Tier assignment (X.A/X.B/X.C)

### Phase α Exit Criteria

- [ ] **α-X1**: All 31 ncAA XMLs: Σq = 0.000 ± 1e-6 (script-verified, automated)
- [ ] **α-X2**: PATCH-01 gate active in Keeper pre-run checks (integrated into production pipeline)
- [ ] **α-X3**: NVT 100 ps warmup stable across all 4 ncAA classes (no bonded-param violations, no mega-force frames)
- [ ] **α-X4**: Cp4 Fix 4 PBSA under patched charges satisfies:
  - Direction: ΔΔG positive reduced OR sign maintained vs Track B baseline (+3.92)
  - CI95 includes at least one literature anchor (-1.4 or -3.0 kcal/mol)
  - z_vs_lit < 2 (X.A regime) — **acceptable**
  - Note: If z_vs_lit slightly > 2 (X.B) but < 3 (not X.B ESCALATE), α still exits with flag for β priority
- [ ] **α-X5**: UPDATE.md v0.7.0-α entry committed to main branch
- [ ] **α-X6**: UPGRADE.md v0.7.0 Phase α checkboxes resolved
- [ ] **α-X7**: Paper 1 preprint abstract + α-complete draft outline ready
- [ ] **α-X8**: v0.7-alpha branch merged to main (or equivalent release tagging)

### Phase α Failure Modes & Rollback

**NVT warmup fails at production scale**:
- Rollback: Phase α halt at α-S4
- Fallback: SciVal emergency dispatch → Option C (Khoury hybrid) 가속
- Impact: γ elements pulled forward to α, timeline +1.5 days

**Post-patch Cp4 PBSA still X.B ESCALATE (z > 3)**:
- Rollback: Phase α completion flagged "partial success"
- Fallback: β immediate activation (multi-system 검증 시급)
- Impact: Paper 1 narrative shifts to "charge fix insufficient, multi-system required"

**Bulk regenerate fails (patch table incompatibility)**:
- Rollback: α-S2/S3 debug, per-ncAA investigation
- Fallback: Subset approach (prioritize MTR, OMW, common ncAAs first)
- Impact: Timeline +2-3 days, partial α exit

**v0.6.x archive corruption**:
- Rollback: Git reset from v0.6.9 commit baseline
- Fallback: Re-archive from git history
- Impact: Minor delay, data preserved via git

### Phase α → β Transition

**Automatic transition** if:
- All α-X1 through α-X8 met
- Cp4 Fix 4 PBSA post-patch shows z_vs_lit < 2

**Manual transition decision** (α-D3) if:
- Partial α success (some criteria fail)
- User prefers Paper 1 α-preprint first
- Rollback/fallback triggered

---

## Phase β — Systematic Validation

**Goal**: Multi-system evidence + Paper 1 publication readiness
**Duration**: ~2-4 weeks wall clock
**Pathway alignment**: S elements + P continuation

### Phase β Decision Points

- **β-D1 (entry)**: Test matrix 확정
  - Minimal (3 systems): 1EBP_MTR13, 7TL8_MTR6, 1YCR_NML22
  - Extended (4+ systems): Add 2QKH_MTR25 or cross-class (non-compstatin)
  - Decided by: User

- **β-D2 (mid-phase)**: If 1 of 3 multi-system fails
  - Continue β → collect more data
  - γ acceleration → root-cause investigation
  - Decided by: User + SciVal

- **β-D3 (exit)**: Paper 1 preprint submission timing
  - Submit now (β-complete) → bioRxiv with current evidence
  - Wait for γ → Journal-ready manuscript
  - Decided by: User

### Phase β Entry Criteria

- [ ] **β-E1**: Phase α exit criteria all met (α-X1 through α-X8)
- [ ] **β-E2**: Paper 1 preprint draft outline ready (from α-X7)
- [ ] **β-E3**: GPU budget ~55-75h reserved
- [ ] **β-E4**: Multi-system test matrix defined (β-D1 resolved):
  - 1EBP_MTR13 Fix 4 PBSA n=5
  - 7TL8_MTR6 Fix 4 PBSA n=5
  - 1YCR_NML22 Fix 4 PBSA n=5
  - (optional: additional systems per β-D1)
- [ ] **β-E5**: Khoury RESP preliminary pipeline designed (not applied yet)
  - OMW XML generation script (parmed prep)
  - Test seed identified
  - `Reference/Khoury 2014/sb400168u_si_002/ffncaa.in` line 778-810 OMW block parsed

### Phase β Scope

- [ ] **β-S1**: Multi-system PBSA Fix 4 (3 systems)
  - Same protocol as Track B Cp4 (Fix 4, PBSA, n=5)
  - Parallel execution where GPU allows
  - Per-system aggregation JSON

- [ ] **β-S2**: Cross-ncAA PBSA verification
  - MTR-containing systems cross-compared
  - NML (different ncAA class) as control
  - Consistency metrics: CI overlap, z_vs_lit distribution

- [ ] **β-S3**: Khoury OMW XML pipeline validated (not deployed)
  - Script: parmed conversion from ffncaa.in prepin to OpenMM XML
  - Test on 1 seed (s42 or s7 from Cp4)
  - Compare: Khoury charges vs current patched charges
  - Output: validation report, not production use

- [ ] **β-S4**: parmed conversion test
  - Full pipeline: Khoury prepin → AMBER prmtop → parmed → OpenMM XML
  - Success = NVT 100 ps stable
  - Output: validation for γ entry

- [ ] **β-S5**: Paper 1 preprint drafting
  - Based on Archi Light v3 §3.4 narrative
  - Incorporate α + β evidence
  - bioRxiv format ready

- [ ] **β-S6**: v0.7.1, v0.7.2 minor releases
  - Per-system results as milestones
  - UPGRADE.md checkbox updates

### Phase β Exit Criteria

- [ ] **β-X1**: Cross-system PBSA consistency documented with numerical criteria:
  - ≥2 of 3 systems show z_vs_lit < 2 (X.A regime), OR
  - CI95 overlap coefficient > 0.5 across systems, OR
  - Cross-system ΔΔG variance < σ_within-system × 1.5
- [ ] **β-X2**: Cross-system decomposition analysis complete
  - PBSA vs GB contribution per system
  - Charge-fix vs protocol-fix attribution
  - System-specific factors identified
- [ ] **β-X3**: Paper 1 preprint submitted to bioRxiv (if β-D3 = submit now)
- [ ] **β-X4**: Khoury RESP OpenMM XML pipeline validated (NVT stable, NaN-free)
- [ ] **β-X5**: UPGRADE.md v0.7.2+ items checkboxes resolved
- [ ] **β-X6**: Post-β γ necessity assessment documented:
  - LOW: all 3 systems X.A under P+A only → γ deferrable to v0.9
  - MEDIUM: mixed results → γ for outlier systems only
  - HIGH: multiple failures → γ immediate priority

### Phase β Scenarios

**Success scenario** (all 3 X.A):
- γ necessity: LOW
- Paper 1 narrative: "PBSA closes the gap for cyclic+ncAA" (strong claim)
- γ deferral possible (v0.9 or on-trigger)

**Partial scenario** (1-2 X.A):
- γ necessity: MEDIUM
- Paper 1 narrative: "System-specific PBSA benefits, Khoury RESP for outliers"
- γ for outlier systems specifically

**Failure scenario** (0-1 X.A):
- γ necessity: HIGH
- Paper 1 narrative: "Honest limitations, γ triggered immediately"
- γ entry accelerated

### Phase β → γ Transition

**Automatic γ trigger** if β failure scenario (0-1 systems X.A)
**User decision γ entry** if β partial/success (per β-D3 + γ trigger conditions)
**γ deferral to v0.9** if β success + Paper 1 reception positive

---

## Phase γ — Root-Cause Resolution

**Goal**: Khoury ff03 RESP full adoption + Journal publication
**Duration**: ~1-2 weeks work within v0.8 release cycle
**Pathway alignment**: Pathway C (Khoury hybrid)

### Phase γ Trigger Conditions

**Any of**:
- **γ-T1 (time-based)**: v0.8 release cycle planning begins → γ auto-scheduled
- **γ-T2 (evidence-based)**: β failure scenario (β-X6 = HIGH) → γ immediate
- **γ-T3 (user-driven)**: Explicit user decision for rigor upgrade (e.g., reviewer request post-preprint)

### Phase γ Decision Points

- **γ-D1 (entry)**: FF mixing risk acceptance
  - Accept: proceed with ff14SB backbone + ff03 NCAA
  - Escalate: consider alternative (SDM alchemical, ATM)
  - Decided by: User + SciVal (based on β-S3/S4 validation results)

- **γ-D2 (mid-phase)**: If FF mixing fails at production scale
  - Continue debug → investigate parameter compatibility
  - Alternative path → SDM alchemical, ATM, or full redesign
  - Decided by: User + SciVal

- **γ-D3 (exit)**: v0.8 release vs v0.9 defer
  - v0.8 immediate → γ deliverables complete, release now
  - v0.9 defer → γ partial success, polish in next cycle
  - Decided by: User

### Phase γ Entry Criteria

- [ ] **γ-E1**: β exit criteria met (at least partial)
- [ ] **γ-E2**: Khoury RESP pipeline validated (from β-X4)
- [ ] **γ-E3**: SciVal verdict on FF mixing (ff14SB backbone + ff03 NCAA)
- [ ] **γ-E4**: FF mixing pilot complete (1 seed × 10 ns, energy conservation verified)
- [ ] **γ-E5**: Paper 1 preprint feedback incorporated (if submitted)
- [ ] **γ-E6**: GPU budget 40-60h reserved

### Phase γ Scope

- [ ] **γ-S1**: Khoury OMW XML production deploy
  - From β-validated pipeline
  - All ncAAs with Khoury equivalent
  - `Reference/Khoury 2014/sb400168u_si_002/` → production `_md_input/`

- [ ] **γ-S2**: Full ncAA re-parameterize (Khoury-compatible subset)
  - MTR (OMW), OMY, NMT, AIB, PAL, ALC, NMW, ALN, NMD, NMG, NMQ, NMR, NMH, OEY, NMA, NMY, NMC
  - Non-Khoury ncAAs: keep α-patched (hybrid state)

- [ ] **γ-S3**: Bonded parameter override validation
  - GAFF2 bonded + Khoury bonded (OMW.frcmod)
  - Compatibility check with amb14SB backbone

- [ ] **γ-S4**: Cp4 Fix 4 n=5 under Khoury RESP + PBSA
  - Ultimate validation run
  - Compare: α (patched charges) vs γ (Khoury RESP)
  - ΔΔG distribution + CI

- [ ] **γ-S5**: Multi-system re-validation under Khoury RESP
  - Same β systems (1EBP_MTR13, 7TL8_MTR6, 1YCR_NML22)
  - Expected: γ results ≥ β results in rigor

- [ ] **γ-S6**: v0.8.0 release preparation
  - CHANGELOG
  - Migration guide from v0.7.x
  - Deprecation notices for α-patched charges

- [ ] **γ-S7**: Paper 1 Journal submission
  - Full manuscript with γ validation
  - Target: JCIM, JCTC, or JPCB

- [ ] **γ-S8**: Keeper gate KHOURY-01 integration
  - RESP consistency check (net charge, atom-type validity)
  - Follows PATCH-01 pattern

### Phase γ Exit Criteria

- [ ] **γ-X1**: Khoury charges applied to ≥MTR, OMW (critical ncAAs for compstatin)
- [ ] **γ-X2**: FF mixing validated at production scale (NVT stable, energy conservation, 10 ns per system)
- [ ] **γ-X3**: Cp4 X.B resolution confirmed under Khoury RESP + PBSA (z_vs_lit < 2)
- [ ] **γ-X4**: Multi-system consistency maintained (β X.A results preserved or improved)
- [ ] **γ-X5**: Paper 1 Journal-ready manuscript
- [ ] **γ-X6**: UPDD v0.8.0 released
- [ ] **γ-X7**: Keeper gate KHOURY-01 integrated (RESP consistency check)

### Phase γ Failure Scenarios

**FF mixing untenable at production**:
- Rollback: Stay at v0.7.x (Phase α state, Khoury opt-in only)
- Fallback path A: **SDM alchemical** (Kilburg-Gallicchio 2018)
  - OpenMM + AGBNP2 implicit solvent
  - Cyclic peptide convergence advantage
  - Paper 1 revision: "Alchemical path for high-rigor cyclic+ncAA ΔΔG"
- Fallback path B: **ATM (Alchemical Transfer Method)** (Gallicchio 2022)
  - Alternative alchemical approach
- Fallback path C: **v1.0+ complete redesign**
  - Different ncAA parameterization framework
  - Extended timeline, major refactor

**γ partial success** (some systems fail):
- Per-system Khoury adoption (not universal)
- Paper 1 honest report: "Khoury RESP for X subset, α-patched for Y subset"
- γ exit partial, backlog for v0.9+

### Phase γ Deferral Conditions

γ can be **deferred (not abandoned)** if:
- β fully succeeds (all 3 systems X.A under P+A only)
- Paper 1 preprint well-received (no reviewer requests for γ)
- Workload constraints (γ becomes v0.9 planning)

**Deferral tracking**:
- Remains in UPGRADE.md backlog (this section preserved)
- Reevaluate at every v0.8.x, v0.9 planning cycle
- **Not "forever deferred"** — explicit reassessment each cycle

---

## Cron Integration

Cron automation supports Phase E execution:

**Auto-monitoring** (every 3h cron):
- Phase α/β/γ entry/exit criteria status
- Active runs progress (GPU utilization, ETA)
- Pending user decisions queue

**Auto-flagging**:
- Phase transition ready (all exit criteria met)
- Phase failure scenario detected (rollback needed)
- User decision overdue (>24h since flag)

**Auto-reporting**:
- Daily summary to UPDATE.md pending entries
- Weekly phase status to `outputs/status/weekly_v07_status_YYYYMMDD.md`

---

## Success Metrics (Pathway E overall)

**End of Phase α** (~1 week):
- ✅ Pipeline integrity restored
- ✅ Cp4 X.B resolved (n=5 empirical)
- ✅ v0.7.0 released
- ✅ α-preprint abstract ready

**End of Phase β** (~1-2 months total):
- ✅ Multi-system evidence (3+ systems)
- ✅ Cross-system consistency documented
- ✅ Paper 1 preprint submitted (bioRxiv)
- ✅ v0.7.x stable

**End of Phase γ** (v0.8 release cycle):
- ✅ Literature-aligned ff03 RESP
- ✅ Full methodology validated
- ✅ Paper 1 Journal submission
- ✅ v0.8.0 released

**End of Pathway E** (v0.8-v0.9):
- ✅ Production-grade cyclic peptide + ncAA pipeline
- ✅ Publication-ready methodology
- ✅ Reference-aligned parameterization

---

## Decision Audit Trail

All Phase E decisions logged at:
- `outputs/decisions/phase_alpha_decisions.md` (α-D1, α-D2, α-D3)
- `outputs/decisions/phase_beta_decisions.md` (β-D1, β-D2, β-D3)
- `outputs/decisions/phase_gamma_decisions.md` (γ-D1, γ-D2, γ-D3)

Format per decision:
- Decision ID (α-D1 etc.)
- Date + timestamp
- Context + options
- User choice + rationale
- SciVal verdict (if applicable)
- Outcome assessment (post-hoc)

---

**End of v0.7 Roadmap — Pathway E**
