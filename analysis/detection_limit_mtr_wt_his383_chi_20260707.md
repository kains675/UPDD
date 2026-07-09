# Detection-limit note — MTR−WT His383 CH-π effect (ncAA funnel)

**Date**: 2026-07-07 · **Regime**: R-11 ranking-only (`ddE_int ≠ ddG_bind`) · **Ties to**: power table #105, ADR-IDEA-01, ADR-IDEA-02
**One line**: The MTR/WT His383 CH-π methyl effect (~1 kcal) is **not recoverable as a quantitative funnel observable on the current FF+sampling**, for two coupled reasons (pose not sampled + effect below the variance floor). It is **not a dead end**: the ranking question is already answered by the endpoint QM anchor, and the quantitative barrier is effort/cost + FF-quality, not a hard physical wall.

---

## 1. Observable
ΔΔE_int(MTR−WT) for the crystallographic His383···residue-4 CH-π contact (2QKI). This is a gas/embedded 2-body interaction energy — **not** ΔΔG_bind (no solvation, entropy, reorganization, multi-residue pocket, environment induction). Every number below carries the R-11 ranking-only tag.

**Anchors (repo-verified 2026-07-07):**
| Method | ΔΔE_int(MTR−WT) kcal/mol | note |
|---|---|---|
| FF (ff14SB+GAFF2) | **−0.557** | = LJ −1.308 + Coul +0.751 (G0a); captures 54% of gold, "right for wrong reasons" |
| **gold FNO-DF-CCSD(T)/aVDZ** | **−1.029** | primary anchor |
| SAPT2+(3)/aVDZ | −1.138 | silver, gold±0.3 confirms |
| wB97X-D/def2-TZVP | −1.422 | over-binding proxy (+38% vs gold) — NOT a magnitude anchor |

SAPT decomposition: disp+exch −0.511, elst −0.390, ind −0.236 (dispersion-dominated).

## 2. Barrier 1 — GIGO (engaged pose not sampled) · EMPIRICAL, P1 census
Crystal ground truth: `2QKI_clean.pdb` Trp4↔His383 imidazole = **3.71 Å** (CH-π real). MD ensembles (mdtraj, PBC-on, backbone-only restraint ⇒ sidechain free ⇒ occupancy valid):

| Trajectory | His383 CH-π occ (≤4.5 Å) | residue-4 sidechain |
|---|---|---|
| Cp4 calib s83 (standard) | **0.0%** | in crystal pocket (Pro384/Thr382/Gly336) but His383 imid **7.4–14 Å** (CH-π relaxed) |
| Cp4 patchoff s443 | **0.0%** | relocated 20–50 Å (Arg275/Asn218/Tyr217) |
| Cp4 patchoff s457 | **0.0%** | relocated 17–53 Å (Ala403/Glu404/Arg394) |

All 2QKI DCDs (local + ExpDATA) are `_restrained`; no free ensemble exists (pipeline uses backbone-restrained MD only). Restraint anchors ncAA backbone N/CA/C only (`run_restrained_md.py:1218`) — sidechain free.
**Consequence**: per-snapshot QM rescoring of this ensemble scores a 7–50 Å geometry → confident null (≈0), not −1.0. Rescoring cannot recover the effect.

## 3. Barrier 2 — VAR (correction below detection floor) · #105
Correction magnitude = |FF − QM| = **0.47 (gold) to 0.87 (wB97X-D) kcal** — a fixed mean-shift, does not reduce σ.
#105 t-corrected E_min detection floor (n=6): **RBFE 2.70**, V3I 3.15, MM-PBSA-MTR13 6.72, Cp4 14.97. Correction is **3–6× below** the best floor. A <1 kcal mean-shift cannot cross a ≥2.7 kcal ranking floor.

## 4. The two barriers share one root — FF dispersion deficit
SAPT: FF under-binds the CH-π dispersion ~2× (captures 54%). P1: FF fails to hold the CH-π contact (relaxes/detaches). These are plausibly the **same cause** — an FF that under-binds aromatic dispersion neither deepens the energy nor holds the pose. So the sampling variance (huge for MTR: σ_btwn 2.99–6.36) is itself the pose-instability, and both barriers trace to the dispersion deficit.

## 5. Remaining methods (NOT a dead end) — ranked
1. **FF dispersion NBFIX (IDEA-02)** — the only lever attacking *both* barriers at once (re-stabilize pose ⇒ ↓GIGO, ↓σ ⇒ ↓VAR, deepen energy). P1 disengagement *re-motivates* it. Gated high-risk (G0a coupled-error: single-channel refit backfires; anisotropy). See [[ADR-IDEA-02]].
2. **Enhanced-sampling engaged ensemble** (umbrella/OPES/REST2 on the His383 CH-π CV) — generates the engaged ensemble, discriminates physics-vs-FF-artifact, measures engaged ΔΔG. v0.7+ project, still G-VAR-gated.
3. **Pose-restrained conditional ΔΔG** — restrain the contact in WT+MTR, compute conditional ΔΔG + restraint correction. Answers "given engagement," not unconditional binding.
4. **Cost-only brute force** — floor ∝ 1/√N. At σ_btwn≈1.4: resolve 0.9 kcal @ z≥2 ⇒ N≈10/endpoint ≈ **280 GPU-h**; 0.5 kcal ⇒ ~880 GPU-h. Physically possible, effort-bounded. **Fails for MTR's own σ (3–6)** until the pose is fixed (loops to #1/#2).

## 6. The honest physical possibility
If the His383 CH-π is genuinely **low-occupancy in solution** (MD suggests), the real ΔG_bind contribution is <1 kcal — **below meaningful detection by physics, not tooling**. The crystal 3.71 Å likely overstates the solution effect. In that branch "we cannot detect it" is the wrong framing; "it is not large enough to matter for binding" is the correct one.

## 7. The answer we DO have (ranking-only)
**gold −1.029** (ddE_int, R-11): MTR's N1-methyl deepens the His383 interaction vs WT by ~1 kcal, sign-confirmed across FF/SAPT0/SAPT2+(3)/wB97X-D/gold. This is the standing mechanistic ranking statement for MTR−WT. FF is a valid garbage-cut filter (correct sign, 54% captured).

## 8. Conclusion
Not "no methods left." (i) The **ranking question is already answered** (gold −1.029). (ii) Quantitative resolution is blocked by an **effort/FF-quality barrier, not a physical wall** — three named levers remain (FF dispersion fix = root, enhanced sampling, cost-scaling), the most promising being the dispersion NBFIX that P1 re-motivates. (iii) There is a **real physical possibility the effect is simply too small in solution to matter** — in which case the detection limit is the scientifically correct terminus, not a tooling gap. This note quantifies *why* the current ensemble cannot move the observable, converting "no pipeline built" into "we measured the limit and named what it would take to cross it."
