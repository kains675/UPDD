#!/usr/bin/env python
"""
backfill_md_status.py — Legacy run _md_status.json 역산 생성 (v0.5 Phase 7)

기존 mdresult/ 산출물 파일 패턴을 분석하여 _md_status.json 을 역산 생성한다.
QM/MM 재계산 없이 rank_results_qmmm 재실행만으로 올바른 ranking 을 확보하는 것이 목적.

사용법:
    python utils/backfill_md_status.py \\
        --project_dir outputs/6WGN_cyclic_htc_NMA_10_20-25 \\
        --topology cyclic_htc

    # Known-case 검증:
    python utils/backfill_md_status.py \\
        --project_dir outputs/6WGN_... \\
        --verify-known w4_0_s2=SUCCESS w1_0_s1=MARGINAL
"""

import os
import sys
import glob
import argparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils_common import (  # noqa: E402
    classify_md_outcome,
    load_md_status,
    write_md_status,
)


def discover_design_stems(md_dir):
    # type: (str) -> list
    """mdresult/ 디렉토리에서 고유 design stem 목록을 추출."""
    stems = set()
    for f in os.listdir(md_dir):
        if f.startswith("_") or f.startswith("."):
            continue
        if os.path.isdir(os.path.join(md_dir, f)):
            continue
        base = f
        for suffix in ["_final.pdb", "_restrained.dcd", "_md.log",
                        "_EXPLODED_dt2fs.log", "_EXPLODED_dt1fs.log",
                        "_EXPLODED_MD.log", "_EXPLODED_NVT.log",
                        "_EXPLODED_COLD_NVT.log", "_EXPLODED_MIN.log",
                        "_EXPLODED_CYCLIC_MIN.log", "_EXPLODED_CYCLIC_RELAX.log",
                        "_partial_md.pdb", "_partial_nvt.pdb"]:
            if base.endswith(suffix):
                stems.add(base[:-len(suffix)])
                break
    return sorted(stems)


def backfill(project_dir, topology="unknown"):
    # type: (str, str) -> dict
    """mdresult/ 파일 패턴 기반으로 _md_status.json 을 역산 생성."""
    md_dir = os.path.join(project_dir, "mdresult")
    if not os.path.isdir(md_dir):
        print(f"[!] mdresult/ not found: {md_dir}")
        sys.exit(1)

    stems = discover_design_stems(md_dir)
    if not stems:
        print(f"[!] No design stems found in {md_dir}")
        sys.exit(1)

    print(f"\n[Backfill] {len(stems)} designs found in {md_dir}")

    md_status = load_md_status(project_dir) or {
        "schema_version": "0.5.0",
        "project": os.path.basename(project_dir),
        "topology": topology,
        "created_by": "backfill_md_status.py",
        "designs": {},
    }

    tier_counts = {}  # type: dict
    for stem in stems:
        outcome = classify_md_outcome(md_dir, stem)
        md_status["designs"][stem] = outcome
        tier = outcome["tier"]
        tier_counts[tier] = tier_counts.get(tier, 0) + 1
        print(f"  {stem:40s} → {tier:10s} "
              f"(passes={len(outcome['passes'])}, "
              f"markers={outcome['stale_markers']})")

    write_md_status(project_dir, md_status)
    print(f"\n[Backfill] _md_status.json written to {project_dir}")
    print(f"  Tier summary: {tier_counts}")
    return md_status


def verify_known(md_status, known_cases):
    # type: (dict, list) -> bool
    """Known-case 검증. 'design=TIER' 형식의 리스트와 대조."""
    all_pass = True
    for case in known_cases:
        parts = case.split("=")
        if len(parts) != 2:
            print(f"  [!] Invalid format: {case} (expected design=TIER)")
            all_pass = False
            continue
        design, expected_tier = parts[0], parts[1].upper()

        found = False
        for key, entry in md_status.get("designs", {}).items():
            short_key = key.replace("design_", "")
            if short_key == design or key == design:
                actual_tier = entry.get("tier", "UNKNOWN")
                if actual_tier == expected_tier:
                    print(f"  [PASS] {design} → {actual_tier} (expected {expected_tier})")
                else:
                    print(f"  [FAIL] {design} → {actual_tier} (expected {expected_tier})")
                    all_pass = False
                found = True
                break

        if not found:
            print(f"  [FAIL] {design} not found in _md_status.json")
            all_pass = False

    return all_pass


def main():
    parser = argparse.ArgumentParser(
        description="Backfill _md_status.json from existing mdresult/ files"
    )
    parser.add_argument("--project_dir", required=True,
                        help="outputs/{project} directory")
    parser.add_argument("--topology", default="unknown",
                        help="Peptide topology (e.g. cyclic_htc)")
    parser.add_argument("--verify-known", nargs="*", default=None,
                        help="Known-case verification: design=TIER pairs")
    args = parser.parse_args()

    md_status = backfill(args.project_dir, args.topology)

    if args.verify_known:
        print(f"\n[Verify] Checking {len(args.verify_known)} known cases...")
        ok = verify_known(md_status, args.verify_known)
        if ok:
            print("  All known cases PASSED ✓")
        else:
            print("  Some known cases FAILED ✗")
            sys.exit(1)


if __name__ == "__main__":
    main()
