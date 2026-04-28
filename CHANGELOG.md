# Changelog

All notable changes to UPDD are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project adheres to honest semantic-style versioning where v0.x process-oriented vs outcome-oriented releases are explicitly distinguished (see [v0.7] §"What this release does NOT claim").

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
