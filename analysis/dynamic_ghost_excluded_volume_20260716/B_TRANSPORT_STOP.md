# B Transport Tuning Stop Record

Date: 2026-07-16 KST

Status: **ACTIVE STOP**

## Trigger

The frozen W4A post-densify structural gate completed 20/20 cells and returned
`PROMOTE_A_CARVE_REDESIGN`.

Primary trigger values:

```text
control median C = 0.334667
s127 H           = 0.497333
H-C              = 0.162667
H/C              = 1.486056
baseline trigger = true (C >= 0.15)
formal verdict   = PROMOTE_A_CARVE_REDESIGN
```

The trigger is cohort-wide: all nine non-s127 carved controls were above the
`0.15` floor. The amplification ratio being just below `1.5` does not change
the OR gate.

Frozen evidence hashes:

```text
080317080bf08a24de5ed0e44eb5e5ca6b5518572fc8cb97c3eadd70f9680513  analysis/w4a_postdensify_dcd_20260716/PREREGISTRATION.md
58a5d05bd001e15e7a76002a80c87321c8b56516e535a894500050e2b8acbcf2  analysis/w4a_postdensify_dcd_20260716/postdensify_gate.json
4d6bf722f83cc9baa600eef7e413cfacc4f3d98860a161d78a02516f70138aa0  analysis/w4a_postdensify_dcd_20260716/postdensify_dcd_summary.json
aeb5b696366454b0e54e4630d29cce0a650ca27765020668e320c571a56817d0  analysis/w4a_postdensify_dcd_20260716/source_inventory.json
```

Source inventory digest:
`d67c4fa420ccb1d6b170a4bef0fe480c4c0cc79395c006f198d48fa1293cbee9`.

## Stopped Actions

Until a new preregistration explicitly supersedes this record:

- do not add H18 windows, cycles, retries, or seed expansions to rescue B;
- do not launch W23A RBFE from the completed 1YCR scaffold cohort;
- do not promote H18 FE diagnostics, which remain non-converged;
- do not run forward-only as a decision estimator;
- do not apply `base=u1` capping;
- do not use a permanent ghost/cavity potential without a correction leg;
- do not overwrite, delete, or silently merge the completed B artifacts.

Allowed work is read-only reanalysis, artifact preservation, control-center
monitoring, and the separately preregistered A-branch P0/P1/P2 work after its
own six-pass authorization.

## Interpretation

This stop is not a claim that every B value is false. W23A remains a
bound-specific ranking/SIGN-only signal with prior caveats. The stop says the
current static carve and transport-tuning path cannot establish a
thermodynamically faithful quantitative result because dynamic water re-entry
persists across controls.

The next branch is
`analysis/dynamic_ghost_excluded_volume_20260716/PREREGISTRATION.md`.

No commit or push was performed.
