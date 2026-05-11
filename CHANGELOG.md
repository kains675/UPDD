# Changelog

All notable changes to UPDD are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project adheres to honest semantic-style versioning where v0.x process-oriented vs outcome-oriented releases are explicitly distinguished (see [v0.7] §"What this release does NOT claim").

## [v0.7.1] — 2026-05-11

UPDD v0.7.1 is the **Capability Level 1 module-suite + multi-layer verification closure release** on top of the v0.7 architectural integrity baseline. The numbers and Tier classifications established at v0.7 are preserved (post-CONECT-v55, post-L387-v54, post-PR-21 frozen baseline); v0.7.1 adds the four additional verification layers (CONECT wraparound, L387 four-layer PBC repair, intra-residue bond integrity, source/test integrity) plus the five-module Capability Level 1 user-facing suite, the simulated active-learning cycle on the six-target ncAA benchmark, and the reproducibility infrastructure required for preprint accompaniment. Today's small comment-hygiene fix (`fe4309b`) closes the v1 → v2 manuscript revision arc against the reviewer's six-point critique.

### Multi-layer PBC verification (Phase β extension of Phase α)

- **L387 v54 four-layer PBC repair** (`utils/extract_snapshots.py` save_snapshots layered fix, commit `34fe7c3`). Layer 0 `_ensure_intrachain_peptide_bonds` registers HETATM/ATOM-boundary peptide bonds into the mdtraj bond graph (mdtraj's PDB loader does not auto-register these); Layer 1 `_unwrap_frame_pbc` invokes `image_molecules` with an explicit anchor molecule chosen as the largest molecule from `topology.find_molecules()`, bypassing the default `guess_anchor_molecules` heuristic that silently fails on large protein–peptide complexes; Layer 2 `_unwrap_binder_manual` (legacy v53 retention) applies a minimum-image centroid shift for multi-chain complexes whose inter-chain centroid distance exceeds box/2. 70-system trajectory scan: 21 INTACT / 5 BROKEN / 2 AMBIGUOUS / 42 skipped pre-patch; 28-system re-extraction post-patch: 25 INTACT / 0 BROKEN / 3 AMBIGUOUS (the 3 AMBIGUOUS are authentic biophysical detachments handled by the §3.6 Light D Hybrid detachment classifier).
- **CONECT v55 atom-index wraparound fix** (`utils/extract_snapshots.py::_inject_ncaa_conect`, commit `ba85b39`). Pre-fix, atom-serial-based residue-map keys collided when wrapped water atoms (post-99 999 → 0, 1, 2, …) overwrote ncAA atom serials in the `by_serial` lookup, causing the residue-resname check to return `"HOH"` and skipping CONECT emission for MTR and other ncAAs. The v55 fix replaces atom-serial-keyed mapping with an atom-index-keyed map and adds hex-extended PDB serials (`int(s, 16) - 0xA0000 + 100000`) for systems with > 99 999 atoms. Cp4 was the sentinel system; the post-v55 ⟨⟨ΔG⟩⟩ anchor moved by Δ = −9.54 kcal/mol relative to the pre-v55 contaminated value (artifact correction, not methodology change). The post-CONECT Cp4 anchor at ⟨⟨ΔG⟩⟩ = −1.823 ± 6.656 kcal/mol (n = 8) lies inside the Magotti 2009 SSOT range (ITC −1.4 / SPR −3.0 kcal/mol).
- **PR-21 auto-diagnostic Level 1 intra-residue bond guard** (commit `c97b131`). Previously undocumented PBC defect class: a single ncAA backbone atom (e.g. MTR6_N) imaged across the periodic box while the rest of the residue stayed in place, leaving N ↔ CA = 113.15 Å (canonical 1.46 Å) while the inter-residue prev_C ↔ ncAA_N bond stayed at 1.40 Å. Layer 0–1–2 inter-residue check passed on this frame even though the residue was structurally broken. `intra_residue_bond_check` scans every covalent bond `(a, b)` with `a.residue.index == b.residue.index` and flags any frame whose bond distance exceeds `INTRA_RESIDUE_BOND_THRESHOLD_ANGSTROM = 5.0` Å. Critical implementation detail: mdtraj's PDB loader does not auto-create intra-residue bonds for HETATM ncAAs, so the scan must include explicit backbone-by-name iteration over the static `_BACKBONE_PAIRS = (("N", "CA"), ("CA", "C"), ("C", "O"))` map. The 7TL8_MTR6 s55 frame snap05_f63 is the empirical first-catch; removing the single defective frame restored σ_btwn from 10.78 → 3.72 (−65 %), CI from 0.84 → 0.39, and z_SE from 1.19 → 7.22.

### Capability Level 1 five-module suite (PR-NEW-A through E)

- **PR-NEW-A — Sequence Input Parser** (schema `sequence_parser/0.1`). Ingests parent sequence + ncAA modification specification, validates against the 25-entry ncAA registry, and emits a canonical design specification.
- **PR-NEW-B v0.2 — Branched ΔΔG Engine with Light D Hybrid detachment metric** (schema `branched_ddg/0.3`, `utils/branched_ddg.py`, commit `6ea98a1`). Paired (variant, WT) ensembles evaluated under identical Stage-4 protocols; advisory `detachment_metric` field (Light D Hybrid implementation) raises `high_detachment` warning for downstream consumers with per-frame distance cutoff at 5.0 Å (`utils/branched_ddg.py:1299`) and aggregate detachment-fraction warning level `_DETACH_WARN_HIGH = 0.5` (L149). Automated bit-identical WT control as a routine reproducibility benchmark.
- **PR-NEW-C — Iterative Cycle Manager** (schema `cycle_manager/0.2`, commit `5f90a05`). Cycle 0 / Cycle 1+ orchestration with multi-assay K_d inverse-variance aggregation and K_d plateau detection.
- **PR-NEW-D — User CLI Stage 1–2 wrappers** (commit `20a0473`). Subprocess wrappers around RFdiffusion + ProteinMPNN + AF2 for de novo discovery orchestration.
- **PR-NEW-E — HTML Report Generator** (schema `html_report/0.2`, commit `ffc7a2e`). Per-cycle ΔΔG ranking, σ-decomposition plots, WT control verification, cycle progression view, detachment banner, print CSS.

### Simulated active-learning cycle

- **`scripts/simulated_al_cycle.py`** (NEW, 645 LOC, schema `simulated_al_cycle/0.1`, commit `7e09e3d`). Active-learning cycle simulator over the six-target ncAA pair benchmark: 3 acquisition functions (UCB / EI / Maximum Entropy with phase-specific switching), calibration model with linear-residual → GPR transition at `cycle ≥ 3 AND n_cumulative ≥ 15` (RBF + WhiteKernel), coverage-probability metric per cycle, and 3 strategies in parallel (`al` / `random` / `greedy` baselines). SciVal C1–C7 pre-validation gates all implemented: Magotti single-anchor SSOT enforcement (Cp4 = −1.4 kcal/mol ITC), `interpretation_regime: "protocol_convergence_demonstration"`, mandatory `prospective_wet_lab_required: true` flag, coverage probability surfacing, σ-target-selectivity disclosure. 10/10 new pytest cases pass.

### Reproducibility infrastructure

- **`environment.yml`, `Dockerfile`, `.github/workflows/ci.yml`** added (commit `1440c28`). Conda env spec for the `qmmm` environment (Python 3.10 + openmm 8.4 + pytest 9.0.3), Docker image build, GitHub Actions CI workflow for the public-release preparation track.
- **Zenodo deposit metadata** for v0.7.1-paper1-v1 release (commits `a183628`, `659ecc8`). DOI: `10.5281/zenodo.20067323`.

### Code unification (Phase 1+2)

- **`utils/updd_cli.py` → `utils/updd_tui.py`** rename (commit `4b3b338`). Closes incident #99/#100 namespace-collision class permanently (the `scripts/updd_cli.py` ↔ `utils/updd_cli.py` import shadow that had caused 14 pytest failures in the regression suite). 18 `.bak` files + 23 closed-phase dispatch scripts archived under `archive/2026-05-05_pre_v1_freeze/`. 11 AmberTools MMPBSA scratch files moved to external drive (`/media/san/ExpDATA/UPDD_proj_Backup/mmpbsa_scratch_20260505/`). New `.gitignore` (3-line minimal: `_MMPBSA_*`, `mdout`, `*.mdcrd.0`).

### Comment hygiene (today)

- **`utils/parameterize_ncaa.py:1488–1492` stale comment refresh** (commit `fe4309b`, 2026-05-11). The comment block at L1488–1492 stated `OPT-IN via UPDD_MTR_AMBER14_PATCH=1` / `Default OFF`, which had been correct at v0.6.6 (pre-flip) but was contradicted by the actual code default at L1493 (`os.environ.get("UPDD_MTR_AMBER14_PATCH", "1") != "0"`, i.e., patch **ON** unless explicitly disabled) after the v0.7 default-flip described above. The contradiction was caught during the Paper 1 v2 adversarial verification of reviewer Critique 1 (production XMLs verified to carry q_N = −0.4157 e and |Σq| ≤ 1×10⁻⁶ across the 4-system sample audit, matching patch-ON state). Code at L1493 is unchanged; comment block now accurately describes the default-ON state, the q_N values for each gating state, and the audit-baseline location at `outputs/_archive/pre_amb14_patch_20260427/`.

### Sampling adequacy + ncAA charge handling (verification artifacts, not v0.7.1 production)

The Paper 1 v2 manuscript-revision arc (paper drafts retained as local-only per project memory policy; not committed in this changelog) produced a patch-on / patch-off MM-PBSA sensitivity sweep on two anchor systems as a separate verification artifact:

- **1EBP_MTR13** (Tier-1 funnel, n = 5/5 each state): ⟨⟨ΔG⟩⟩ envelope = −0.654 kcal/mol; σ_btwn 3.188 → 2.992 (−6 %); CI 0.156 → 0.142. Tier and sign preserved.
- **2QKI_Cp4** (Tier-3 multi-basin, n = 8/8 each state): ⟨⟨ΔG⟩⟩ envelope = −2.873 kcal/mol; σ_btwn 6.656 → 8.000 (+20 %); CI 3.651 → 1.704. Tier and sign preserved; both patch states retain the ⟨⟨ΔG⟩⟩ within or adjacent to the Magotti SSOT range.

The sweep is a methodology-disclosure artifact, not a v0.7.1 production deliverable; it does not amend the v0.7 / v0.7.1 frozen Table 1 numbers. The sensitivity envelope on Cp4 (−2.87 kcal/mol per family-level mean) is non-trivial relative to the Magotti SSOT magnitude (−3.0 to −1.4 kcal/mol) but does not invert sign or change Tier — the patch is rank-preserving on the validated subset. The sweep workspace is archived at `/media/san/ExpDATA/UPDD_proj_Backup/patchoff_sweep_20260511/` (3.5 GB, 18 items including aggregation script + envelope tables + 13 per-seed mmpbsa_summary.json files).

### Verification protocol — second-pass evidence

The v0.7.1 release builds on the v0.7 verification-protocol-as-contribution framing (catches its own bugs before publication) with second-pass empirical evidence:

1. **CONECT v55 caught the +7.72 → −1.823 Cp4 anchor correction** via the post-v54 PBC repair sweep, confirming that the pre-v55 +7.72 anchor was a serialisation artifact rather than a force-field methodology gap.
2. **PR-21 intra-residue bond guard caught snap05_f63** on 7TL8_MTR6, confirming that mdtraj's missing-bond default for HETATM ncAAs requires an explicit extraction-time contract, and demonstrating that the Convergence Index Framing B (QC-integrity auditor) is empirically triggerable.
3. **Today's comment-refresh closes a stale-documentation-vs-code gap** that the reviewer's Critique 1 surfaced. The audit infrastructure described in the v0.7 changelog continued to enforce the patch-ON state on all 24 production MTR XMLs (verified via 4-system sample audit during the Paper 1 v2 adversarial verification); only the comment text had drifted.

### Test coverage

294 pass / 9 skip / 1 fail (pre-existing `test_admet_filter::test_lipinski_constants_are_module_level` MW = 500 vs 1200 semantic mismatch — separate triage). Baseline preserved across the v0.7 → v0.7.1 progression. Today's comment-refresh commit is comment-only and the test baseline is unchanged.

### Known limitations (carried forward from v0.7 + additions)

- **Cp4 sign-significance to Tier-2 not yet achieved** (n = 8, z_SE = 0.78). Promotion requires N → 30 sampling expansion under realistic √N scaling, with the multi-basin caveat that the Convergence Index Framing A classification (CI = 3.65 → Tier-3) suggests σ_btwn structure may resist convergence under naïve √N scaling.
- **N-methyl ncAA-class backlog and MLE Σq = +0.165 caveat** for 3IOL_NML20 and 1YCR_NML22 unchanged from v0.7; bulk regeneration is the v0.8 charge-methodology work bundled with Khoury full RESP integration.
- **Khoury 2014 same-system R² = 0.388 reference value** is the conservative empirical reference for ncAA-extended PBSA-class methodology (best transparently reported same-system canonical+ncAA correlation in the published literature to date; canonical-only same-system MM-PB(GB)SA benchmarks reach R² ≈ 0.5–0.7 per Hou 2011 / Sun 2019 / Tuccinardi 2022).
- **Patch coverage TRP-derived only** (currently MTR) — N-methyl / PRO / Phospho overlays remain scheduled for v0.8.

### What v0.7.1 does NOT claim

- **NOT** a quantitative match across all six target families — only sign-significance demonstrated on 1EBP_MTR13 (X.B); 7TL8_MTR6 carries Tier-1 X.B (sign-significance only, no literature anchor cited — v2 relabel from v1 X.A); Cp4 carries magnitude consistency with Magotti SSOT but sub-threshold z_SE.
- **NOT** a closed-loop active-learning cycle with wet-lab K_d feedback — the simulated AL cycle is a synthetic-oracle methodology demonstration under explicit `interpretation_regime: "protocol_convergence_demonstration"` framing. The operational closed-loop is the v0.8 / v0.9 target.
- **NOT** the Khoury full RESP integration — the targeted backbone-N amb14SB overlay closes the Σq budget for the TRP-derived MTR class; full RESP refit of the side-chain GAFF2 charges remains deferred to v0.8.

## [v0.7] — 2026-04-28

UPDD v0.7 is an **architectural integrity release**. It establishes the methodologically defensible foundation for subsequent quantitative validation. The charge-correct baseline (post-Phase α) is the methodologically meaningful starting point for Phase γ improvement; v0.7 does not claim a quantitative match to experimental binding-affinity anchors.

### Major changes

- **`UPDD_MTR_AMBER14_PATCH` default flipped from OFF to ON** (Option A; `utils/parameterize_ncaa.py` L1469). The amb14SB Trp-derived charge overlay is now active by default for TRP-derived non-canonical amino acids (currently MTR; gate scope unchanged).
- **24 MTR production force-field XMLs regenerated** with backbone-charge sum Σq = 0 (machine-epsilon tolerance). Pre-patch state had Σq = −0.187 e on the MTR backbone N, biasing the Cp4 ΔΔG calculation by approximately +3.80 kcal/mol. Pre-patch XMLs preserved under `outputs/_archive/pre_amb14_patch_20260427/`.
- **Cache architecture cleanup (Option C v2)**. The MM-PBSA / MM-GBSA / restrained-MD readers (`scripts/run_mmpbsa.py`, `utils/run_mmgbsa.py`, `utils/run_restrained_md.py`) now prefer the per-output local `params/<xml_resname>_gaff2.xml` over the manifest's `/tmp/calib_params/...` indirection. Resolution is alias-aware (`manifest.get("xml_resname", ncaa_code)`), so `NML→MLE` resolves correctly. Legacy fallback to the manifest path is retained behind a `warnings.warn` (provenance-only safety net). The `/tmp/` cache is no longer authoritative for readers.
- **Audit infrastructure (PATCH-01 backbone-N gate + schema regression test)**. `utils/audit_charge_consistency.py` adds the PATCH-01 gate (strict mode raises `Patch01BackboneNError` on pre-patch state); `Patch-01_schema_regression` C5 module (7/7 pytest cases) blocks the Cron #34 `"none"` outlier pattern.

### Scientific validation

- **1EBP_MTR13 sign-significance promotion** (X.C → X.B): post-patch z_SE = 2.09 (vs pre-patch 1.88; threshold 2.0). The 1EBP_WT control is bit-identical pre/post (Δ = 0.00), demonstrating textbook intervention isolation — only MTR atoms changed, WT has no MTR, identical PBSA output.
- **σ_btwn vs σ_w decomposition framework** established (this release's primary scientific contribution beyond the architectural fix). The two ensemble-quality axes change in opposite directions under charge correction: between-replicate scatter σ_btwn tightens or stays stable (Cp4 7.12 → 5.50, −23 %; MTR13 2.82 → 3.00, +6 %) while within-replicate scatter σ_w loosens (Cp4 3.62 → 6.40, +77 %; MTR13 3.24 → 4.91, +52 %). Interpretation: charge fix transitions ensemble from artifact-narrow (replicates trapped in different bias-stabilised basins) to genuine-broad (correct electrostatics enabling proper sampling). This decomposition is generalisable beyond UPDD to any MM-PBSA pipeline with non-canonical amino-acid charge handling.
- **Multi-system σ_w within the Genheden 2015 strict regime** (< 10 kcal/mol): Cp4 σ_w = 6.40, 1EBP_MTR13 σ_w = 4.91.

### Strategic positioning

UPDD v0.7 establishes the **Capability Level 1** foundation: a decision support tool with non-canonical amino acid (ncAA) expansion as the primary differentiator. The release implements the **iterative cycle workflow** (Cycle 0 de novo discovery + Cycle 1+ wet-lab-feedback-driven refinement) with a **branched ΔΔG architecture** (same-protocol pairwise WT vs variant comparison). Capability Level 2 (calibrated quantitative, post-pharma access with wet lab data) and Level 3 (universal quantitative, long-term goal across targets) are explicitly deferred per the honest progression framework. All structures in the pipeline are computationally derived (AF2-predicted + MD); the same-protocol intervention isolation was verified at Phase α (1EBP_WT bit-identical control).

### Verification protocol

This release's verification design was empirically exercised twice during Phase α and caught both errors before production contamination:

1. **Auto-monitoring layer** (Cron #41, 2026-04-27 18:30 KST): caught a silent functional failure where the Step 6+7 v1 PBSA outputs were bit-identical to pre-patch baseline (impossible if charges differed). Root cause: manifest `xml_path` field pointing into the volatile `/tmp/calib_params/<run>/` cache that retained the pre-patch XML. Recovery: 4-location sync (cache + `_md_input/` + manifests + clean `_postpatch` directories) before re-running Step 6+7 v2.
2. **Pre-merge gate layer** (Option C v1 halt, 2026-04-28 ~00:50 KST): caught a structural inconsistency where the dispatch's audit probe used `f"{ncaa_code}_gaff2.xml"`, mismatching NML manifests whose `xml_resname` is `MLE` (10 BREAK candidates). Recovery (v2): 1-line fix `xml_resname = manifest.get("xml_resname", ncaa_code)` for alias support, then 12 min wall-clock to completion (well under the 30–45 min budget). Pre-merge audit re-ran with 37 OK / 0 BREAK.

Test coverage: pre-Option-C regression 109 pass / 9 skip / 0 fail; post-Option-C regression 117 pass / 9 skip / 0 fail (Δ = +8 matches the 8 new dedicated tests T1–T7 + legacy-compat). PATCH-01 gate: 24/24 MTR _test XMLs PASS. 4-method Σq verification: 24/24 PASS. Schema regression C5 module: 7/7 pytest PASS.

### Known limitations

- **Cp4 quantitative literature match unattained**. Post-patch ΔΔG = +7.72 kcal/mol (CI95 [−2.18, +17.62]) **excludes** the Magotti 2009 SSOT anchor (−1.4 kcal/mol ITC, −3.0 kcal/mol SPR; DOI 10.1002/jmr.972). Pre-patch ΔΔG = +3.92 (CI95 included both anchors) is reinterpreted as a methodologically conditioned result, not a robust quantitative match. The +3.80 kcal/mol shift is the magnitude of the charge-bias contribution to Cp4 ΔΔG within the PBSA + ff14SB methodology family.
- **Patch coverage is TRP-derived only**. The amb14SB charge overlay applies to TRP-parent ncAAs (currently MTR); N-methyl class (NMA, NML/MLE, NMV, NMK, NMR, NMQ, MEA, SAR), PRO class (DPR, HYP), and Phospho class (SEP, TPO, PTR) await dedicated overlays — PR-6 (N-methyl, 🔴 HIGH for v0.7.1), PR-7 (PRO), PR-8 (Phospho).
- **Cross-FF junction bonded parameters incomplete**. Bonds and angles are bridged at the CM–NE1 junction (sourced from DNA.OL15 5-methylcytosine analogues); torsions and improper dihedrals across the junction are not yet defined (verdict §3.10, SciVal Phase α V1 §2.2). Phase γ Khoury full RESP integration (PR-18) is the planned axis to address charge methodology comprehensively.
- **Sampling is far from the Genheden 2015 practical SE threshold**. Current n=5 gives SE ≈ 2.46 kcal/mol on Cp4 σ_btwn. SE < 1 requires n ≈ 30 for Cp4, ≈ 9 for 1EBP_MTR13. v0.7 explicitly admits sign-significance focus; quantitative refinement requires replicate expansion (Phase β).
- **Khoury 2014 same-system R² = 0.388 ceiling**. Even after Phase γ Khoury full RESP integration, ~3–4 kcal/mol Cp4 residual gap is expected (verdict §3.10); Phase γ' (explicit solvent + entropy / polarisable FF / QM-MM / alchemical) is the conditional axis if needed (v0.9).
- **Level 1 user-facing features pending** (verdict §5, PR-NEW-A through E). Sequence input parser, branched ΔΔG engine, iterative cycle manager, user CLI, and HTML report generator are scoped for v0.7.x progression. The branched ΔΔG architecture is implemented and verified at Phase α (manual workflow); user-facing automation is the next milestone.

### What v0.7 does NOT claim (cross-ref verdict §4.8)

- **NOT** a quantitative match to compstatin binding free energies — only sign-significance verified on one system (1EBP_MTR13).
- **NOT** a chemical-insight contribution beyond ranking — the R-11 ranking-only regime is intact; absolute ΔG values from this release are not calibrated against experiment and should not be quoted as predictive of Kd.
- **NOT** a resolution of all charge issues — PR-6 (N-methyl), PR-7 (PRO), PR-8 (Phospho) backlogs are explicit.
- **NOT** a guarantee that Phase γ Khoury RESP will close the literature gap — the Khoury R² = 0.388 ceiling is acknowledged.
- **NOT** a verification-protocol claim of bug-free pipelines — two errors were caught in Phase α; more may exist.
- **NOT** a Capability Level 2 or Level 3 release — quantitative calibration (Level 2) requires wet lab data partnership; universal cross-target prediction (Level 3) is a long-term goal currently unachievable by any field method.

v0.7 is honestly framed as a process-oriented foundation milestone (architectural integrity + sign-significance + σ decomposition methodology), not a final accuracy claim. Outcome-oriented quantitative validation is deferred to v0.8 (Phase γ Khoury full RESP) and conditionally v0.9 (Phase γ' methodology refinement). See `v07_user_verdict/v0.7 Section 4 finalized.md` §4.10 for the full significance statement.
