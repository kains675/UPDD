#!/usr/bin/env python
"""
T1 #1 N→30 expansion seed selector — deterministic numpy RNG audit trail.

Per verdict_t1_n30_expansion_20260522.md §C4: seed integers for the
Cp4 + WT expansion to N=30 must be drawn from a deterministic numpy RNG
(default_rng(SeedSequence(entropy=20260522))) so that the seed set is
reproducible from an audit trail and free from analyst cherry-picking.

Outputs JSON manifest:
  - seeds_cp4_new:  17 fresh Cp4 seed integers (avoids 15 existing Cp4 integers)
  - seeds_wt_new:   20 fresh WT seed integers (avoids 10 existing WT integers)
  - rng_entropy:    20260522 (the SeedSequence entropy used)
  - existing_cp4:   the 15 integers excluded from the Cp4 draw
  - existing_wt:    the 10 integers excluded from the WT draw
  - manifest_sha256: SHA-256 of the canonical JSON for cross-machine verification

Run:
  /home/san/miniconda3/envs/qmmm/bin/python scripts/t1_n30_seed_selector.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import List

import numpy as np

# ============================================================================
# Existing cohort integers (from filesystem 2026-05-22 enumeration)
# ============================================================================
# Cp4: 13 effective cohort + 2 LOAD_FAILED reextract entries (s113, s173) per
# t1_phase1_5_fresh.sh comment. All 15 integers excluded from the new draw.
# s19_reseed55 occupies integer 55 (the reseed value), not 19, since 19 is
# already used by the original s19 directory.
EXISTING_CP4 = sorted(
    [7, 19, 23, 42, 53, 55, 83, 89, 101, 113, 127, 163, 173, 199, 251]
)

# WT: 10 dirs, all effective
EXISTING_WT = sorted([7, 19, 23, 42, 83, 101, 127, 163, 199, 251])

# ============================================================================
# Selection parameters
# ============================================================================
RNG_ENTROPY = 20260522  # date of pre-launch verdict 2026-05-22 — audit anchor
SEED_RANGE_LOW = 1
SEED_RANGE_HIGH = 999  # avoid 4-digit; matches Phase 1.5 convention (s7-s251)
N_CP4_NEW = 17  # 30 - 13 effective
N_WT_NEW = 20   # 30 - 10


def draw_seeds(rng: np.random.Generator, n: int, exclude: List[int]) -> List[int]:
    """Draw n unique integers from [SEED_RANGE_LOW, SEED_RANGE_HIGH] excluding `exclude`."""
    excluded = set(exclude)
    picked: List[int] = []
    while len(picked) < n:
        candidate = int(rng.integers(SEED_RANGE_LOW, SEED_RANGE_HIGH + 1))
        if candidate in excluded:
            continue
        excluded.add(candidate)
        picked.append(candidate)
    return sorted(picked)


def main() -> int:
    ss = np.random.SeedSequence(entropy=RNG_ENTROPY)
    rng = np.random.default_rng(ss)

    # Cp4 draw first, then WT (sequential — RNG state carries over)
    seeds_cp4_new = draw_seeds(rng, N_CP4_NEW, EXISTING_CP4)
    seeds_wt_new = draw_seeds(rng, N_WT_NEW, EXISTING_WT)

    # Verify uniqueness internal
    assert len(set(seeds_cp4_new)) == N_CP4_NEW, "Cp4 duplicates"
    assert len(set(seeds_wt_new)) == N_WT_NEW, "WT duplicates"

    # Cp4 and WT integer sets are allowed to overlap (different output dirs),
    # consistent with Phase 1.5 baseline (e.g. s127 used in both Cp4 and WT).

    manifest = {
        "schema": "t1_n30_seed_selector/0.1",
        "generated_utc": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "verdict_anchor": "verdict_t1_n30_expansion_20260522.md §C4",
        "option": "β (Cp4 + WT both N→30 symmetric, per §C1)",
        "rng_entropy": RNG_ENTROPY,
        "rng_spec": "numpy default_rng(SeedSequence(entropy=20260522))",
        "seed_range": [SEED_RANGE_LOW, SEED_RANGE_HIGH],
        "existing_cp4_integers_excluded": EXISTING_CP4,
        "existing_wt_integers_excluded": EXISTING_WT,
        "n_cp4_new": N_CP4_NEW,
        "n_wt_new": N_WT_NEW,
        "seeds_cp4_new": seeds_cp4_new,
        "seeds_wt_new": seeds_wt_new,
        "cp4_target_total": len(EXISTING_CP4) - 2 + N_CP4_NEW,  # 13 effective + 17 = 30
        "wt_target_total": len(EXISTING_WT) + N_WT_NEW,         # 10 + 20 = 30
        "note_cp4_failed_dirs": [113, 173],
        "note_cp4_reseed55_alias": "s19_reseed55 directory, integer seed 55",
    }

    # Compute SHA-256 of canonical JSON for cross-machine verification
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    manifest["manifest_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()

    # Write manifest
    proj = Path(__file__).resolve().parent.parent
    out_dir = proj / "outputs" / "analysis" / "t1_n30_seed_selection"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "t1_n30_seed_manifest.json"
    with open(out_path, "w") as f:
        json.dump(manifest, f, indent=2)

    # Human-readable summary
    print("=" * 64)
    print("T1 #1 N→30 expansion seed selection")
    print("=" * 64)
    print(f"RNG entropy: {RNG_ENTROPY}")
    print(f"Cp4 new seeds (n={N_CP4_NEW}): {seeds_cp4_new}")
    print(f"WT  new seeds (n={N_WT_NEW}): {seeds_wt_new}")
    print(f"Cp4 target total: {manifest['cp4_target_total']} (13 effective + {N_CP4_NEW})")
    print(f"WT  target total: {manifest['wt_target_total']} (10 + {N_WT_NEW})")
    print(f"Manifest SHA-256: {manifest['manifest_sha256']}")
    print(f"Manifest path:    {out_path.relative_to(proj)}")
    print("=" * 64)
    print()
    print("Bash array form (paste into dispatcher):")
    print(f"  CP4_FRESH=({' '.join(str(s) for s in seeds_cp4_new)})")
    print(f"  WT_FRESH=({' '.join(str(s) for s in seeds_wt_new)})")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
