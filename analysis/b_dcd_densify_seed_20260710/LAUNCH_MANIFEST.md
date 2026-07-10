# B-DCD Launch Manifest

Date: 2026-07-10

This manifest records the preregistered commands for the W4A H18 densified
carved/uncarved expansion. Run dry-runs first; do not mix these outputs with the
older 400-cycle C1 roots.

## Shared Values

```bash
PY=/home/san/miniconda3/envs/atm/bin/python
LAUNCHER=scripts/trackb_inplace_rbfe_production.py
ROOT=outputs/_trackb/w4a_bdcd_h18_20260710
SEEDS=s7,s19,s23,s42,s83,s101,s127,s163,s199,s251
L1=0.025,0.05,0.1,0.2,0.3,0.4,0.5
L2=0,0.05,0.1,0.15,0.2,0.3,0.4,0.5
COMMON="--endpoints wt --seeds $SEEDS --directions dplus,dminus --twocopy --displacement-nm 4.0 --mutation w4a_trp_ala_res4 --lambda1-rampdown $L1 --lambda2-rampdown $L2 --n-windows-half 6 --softcore-band 2 --n-apex-bridge 0 --n-cycles 800 --md-steps-per-cycle 250 --platform CUDA --timestep-fs 1.0 --genuine-decouple-nm 1.2 --staged-min --reseed-perm-seed 20260618 --reseed-endpoint --reseed-endpoint-band-lambda2-max 0.25 --reseed-endpoint-equil-steps 2000 --mintimeid 600 --pool --device-index 0"
```

## Dry-Run Commands

```bash
$PY $LAUNCHER --dry-run --leg free  $COMMON --max-concurrent 2 --out-root $ROOT/carved_free   --carve-void-waters --carve-cutoff-nm 0.26
$PY $LAUNCHER --dry-run --leg bound $COMMON --max-concurrent 1 --out-root $ROOT/carved_bound  --carve-void-waters --carve-cutoff-nm 0.26
$PY $LAUNCHER --dry-run --leg free  $COMMON --max-concurrent 2 --out-root $ROOT/uncarved_free
$PY $LAUNCHER --dry-run --leg bound $COMMON --max-concurrent 1 --out-root $ROOT/uncarved_bound
```

## Production Commands

Use the same four commands without `--dry-run` after dry-run preregistration and
GPU checks pass.

Expected rough runtime from the 6-seed 400-cycle C1 run:

- free, carved or uncarved: at least several hours per 10-seed H18/800 arm;
- bound, carved or uncarved: roughly day-scale per 10-seed H18/800 arm on the
  16 GB 5070 Ti with bound concurrency fixed at 1.

## Post-Run Structural Gate

After each completed FE arm, run the local-water DCD sidecar on at least:

- W4A carved `s127 bound/dminus`
- W4A carved matched bound/dminus controls, including `s163`
- any new maximum-leverage seed from the densified scalar analysis
- corresponding uncarved bound/dminus controls if C1 remains INDETERMINATE

Then evaluate:

```bash
$PY analysis/b_dcd_densify_seed_20260710/bdcd_gate.py \
  --summary <post_densify_dcd_summary.json> \
  --phase post_densify \
  --json-out analysis/b_dcd_densify_seed_20260710/post_densify_gate.json \
  --report-out analysis/b_dcd_densify_seed_20260710/post_densify_gate.md
```
