# UPDD Local Control Center Preregistration

Status: `FROZEN BEFORE CONTROL-ENABLED SERVICE`

Date: 2026-07-16

Scope: local operational control and observability. This work must not change
any force field, Hamiltonian, lambda schedule, seed cohort, estimator, or
scientific acceptance criterion.

## ARCHI

The browser is a presentation client only. Compute ownership is:

```text
Streamlit UI (127.0.0.1)
  -> authenticated localhost HTTP API
    -> SQLite registry and append-only events
      -> allow-listed adapter
        -> systemd --user transient service/cgroup
          -> immutable argv job
```

Canonical state axes remain independent:

- execution: `DRAFT`, `PREFLIGHT`, `QUEUED`, `RUNNING`,
  `PAUSE_REQUESTED`, `PAUSED`, `COMPLETE`, `FAILED`, `STOPPED`, `ORPHANED`;
- artifact: `PRESENT`, `PARTIAL`, `MISSING`;
- science: `VALID`, `PARTIAL`, `INVALID`, `UNKNOWN`;
- ranking eligibility: `true`, `false`, or `null`.

The registry stores immutable canonical JSON specs and SHA-256 hashes. A
parameter change creates a new draft revision. Existing manual/tmux work is
imported read-only unless its identity and resume contract are proven.

Supported adapters, in rollout order:

1. W4A post-densify DCD sidecar;
2. 1YCR WT scaffold expansion;
3. Track B pool production;
4. Track A pipeline;
5. analysis commands;
6. legacy discovery, read-only only.

No adapter accepts shell text. Interpreters, script paths, cwd, output roots,
and resource names are absolute and allow-listed. The UI exposes no arbitrary
command, broad process kill, git commit, or git push operation.

## SCIVAL

Verdict: `CONDITIONAL APPROVE`.

Conditions:

- Control-token checks occur only at existing scientific unit boundaries:
  DCD cell, scaffold seed, Track B replicate/direction pool unit, or Track A
  pipeline step.
- A pause request drains the currently running unit and launches no new unit.
  It never uses `SIGSTOP`.
- With `UPDD_CONTROL_TOKEN` absent, runner behavior is unchanged.
- Resume requires the same spec hash and input digest. Only completed outputs
  passing the workflow's existing manifest/integrity gate may be skipped.
- Partial output is preserved under `_archive/`; it is never stitched or
  silently deleted.
- Track B resume must not archive or replay a valid completed replicate.
- A process exit code cannot set science validity or ranking eligibility.
- Current W4A H18 and DCD conclusions remain ranking/SIGN-only where
  applicable. The control center does not authorize new B transport tuning.

Stop on any control action that changes scientific parameters in place,
accepts a source/spec digest mismatch, starts a second host GPU lease, treats
legacy work as controllable without identity proof, or reports scientific
validity from process status alone.

## Control Semantics

Launch:

1. validate the adapter and canonical spec;
2. run adapter preflight;
3. freeze the spec hash;
4. acquire the named resource lease transactionally;
5. start the exact systemd user unit and record its cgroup/unit identity.

Pause and resume:

- `PAUSE_REQUESTED` is written atomically to a per-job token;
- the runner acknowledges `PAUSED` at a safe boundary and exits with code 75;
- systemd treats code 75 as a controlled success condition;
- resume revalidates spec and inputs before starting the same immutable argv;
- valid completed units are skipped; partial units use the existing archive
  contract.

Stop:

- requires a short-lived, job-bound confirmation nonce;
- targets only the recorded systemd unit/cgroup;
- records request and result as append-only events;
- a forced kill is a separate explicit action.

## Acceptance Criteria

1. Restarting UI or API does not stop a compute unit.
2. Reconciliation classifies a missing runtime as `ORPHANED`, `PAUSED`, or
   `COMPLETE` from persistent evidence.
3. Two jobs cannot hold the same GPU lease.
4. Pause drains the active unit and prevents the next unit from launching.
5. Resume requires the same spec/input digest and preserves partial output.
6. Declared, default, effective, and source parameters remain distinguishable.
7. Execution, artifact, science, and ranking states are independent.
8. Every action creates an append-only timestamped audit event containing job
   and spec identity, action, result, and reason.
9. No arbitrary shell, broad kill, automatic commit, or automatic push exists.
10. Read-only discovery is verified before control is enabled for an adapter.

## Verification Matrix

- unit tests: schema/hash/status transitions, event immutability, lease race,
  path allow-list, log offsets, stop nonce;
- integration tests: systemd-owned harmless multi-unit fixture, API restart
  survival, boundary pause/resume, exact stop, reconciliation;
- adapter tests: completed DCD and scaffold imports, Track B/Track A no-token
  parity, exact-spec completed-unit skip;
- browser tests: desktop and mobile screenshots, no overlap, nonblank data,
  controls disabled for read-only/completed jobs;
- live service check: localhost binding, bearer authentication, persistent DB,
  discovered current jobs, and host metrics.

No production GPU computation is authorized by this preregistration.
