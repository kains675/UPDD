# 1YCR WT Scaffold Seed Expansion

Status: `FROZEN BEFORE NEW-SCAFFOLD OUTPUT`

Date: 2026-07-16

Regime: scaffold generation only. No RBFE launch or free-energy interpretation
is authorized by this document.

## Objective

Add the five missing `1YCR_WT` scaffold seeds required for a future matched
W23A n=10 cohort:

`s83,s127,s163,s199,s251`

The existing canonical seeds are:

`s7,s19,s23,s42,s101`

## Provenance Audit

The five existing scaffolds were generated on 2026-04-21/22 after commit
`3d2a10c`. Their common input is byte-identical:

- source: `outputs/1YCR_WT_calib_s7/_md_input/1YCR_WT.pdb`
- SHA-256:
  `20f87791637e71a189e5e641f3ec5b807a845aca60bca6f1942ac879cda00cf7`
- raw structure: 818 atoms, chain A 85 residues, chain B 13 residues;
- no ncAA or required cofactor declaration.

The existing five `1YCR_WT_md.log` files each contain 100 numeric reporter rows
ending at simulation step 2,525,000. This is 25,000 equilibration steps plus a
2,500,000-step production request. The archived DCD checked for `s7` has 500
frames, matching the default `n_steps // 500` interval.

Historical generator blob at `3d2a10c`:

`af0c673571d906f0e616730646da5b0f8830e26c`

Current generator before launch:

- path: `utils/run_restrained_md.py`
- git blob: `fd3249e9c1feff9571097f140561d70a438ce202`
- SHA-256:
  `b27bc5794cf7bcd61682a67586077024c3d3e15b9a5acc4763718c4a22b5ce10`

The complete diff from the historical generator to the current generator has
four change classes:

1. scratch-directory routing;
2. GPU device selection externalized with default device 0;
3. ncAA manifest path handling, not consumed by `--ncaa_label none`;
4. an optional DCD interval environment override whose unset behavior is the
   historical default.

No WT force-field, system-construction, equilibration, integrator, barostat,
step-count, or default DCD cadence change is present in that diff. The new run
explicitly fixes device 0 and removes `UPDD_MD_DCD_INTERVAL` from the child
environment.

Residual limitation: a historical conda lockfile from 2026-04-21 is not
available. The current qmmm environment is recorded in the source inventory,
but exact package-level identity to the old run is unknown. This is reported as
provenance uncertainty, not silently treated as exact bitwise replay.

## Frozen Protocol

For each new seed, serially run:

- interpreter: `/home/san/miniconda3/envs/qmmm/bin/python`;
- generator: `utils/run_restrained_md.py`;
- input: a private copy of the frozen common input in
  `outputs/1YCR_WT_calib_s<seed>/_md_input/`;
- output: `outputs/1YCR_WT_calib_s<seed>/mdresult/`;
- `--steps 2500000`;
- `--topology linear`;
- `--ncaa_label none`;
- `--seed <seed>`;
- `--platform CUDA`;
- `--dt_fs 2.0`;
- default binder chain B, dispersion `auto`, and graph policy `strict`;
- `UPDD_MD_CUDA_DEVICE=0`;
- `UPDD_MD_DCD_INTERVAL` absent.

Despite the generator name, the WT path has no ncAA backbone restraint and is
an unrestrained WT MD path. That behavior matches the existing WT cohort and
must not be changed only for the added seeds.

The five cells run in the frozen seed order shown above. The queue starts only
after the W4A post-densify DCD cohort reaches `COMPLETE` and its frozen
structural gate has been evaluated. It never overlaps that GPU workload.

## Completion Gate Per Seed

Every seed must satisfy all conditions:

- child process completed and stdout contains one structured
  `md_complete` event with status `SUCCESS_WT`;
- batch summary records one successful WT MD and zero failure/partial result;
- input copy SHA-256 matches the frozen common input;
- no `EXPLODED`, partial-recovery, NaN, traceback, or CUDA error marker;
- `1YCR_WT_md.log` contains exactly 100 finite numeric rows;
- final reporter step is 2,525,000;
- temperature is finite and positive, and volume is finite and positive;
- `1YCR_WT_final.pdb` parses with finite coordinates and contains chains A/B;
- `1YCR_WT_restrained.dcd` contains exactly 500 frames and the same atom count
  as the final PDB;
- a seed manifest records command, code/input hashes, environment versions,
  output hashes, and validation fields.

`SUCCESS_WITH_WARNING`, `PARTIAL_SUCCESS`, or a final PDB recovered after an
explosion does not count as a valid scaffold for seed expansion.

## Restart And Stop Rules

- A valid completed seed is resume-skipped.
- An interrupted root without a final PDB or completed manifest is moved to
  `outputs/_archive/1ycr_wt_seed_expansion_20260716/` before a clean seed replay.
- A root containing a final PDB but lacking or failing its completed manifest
  is an ambiguous scientific result and stops the queue; it is not
  automatically replayed.
- Any generator/input hash drift, target-root collision, GPU overlap, validation
  failure, or environment setup failure stops the queue.
- Existing five scaffold roots and all prior W23A outputs remain untouched.

## Downstream Block

Creating these scaffolds does not authorize W23A RBFE production. W23A remains
blocked until the active W4A structural/transport branch is resolved and a
separate matched-cohort W23A preregistration passes its science and integrity
gates.
