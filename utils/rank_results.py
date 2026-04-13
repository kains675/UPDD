#!/usr/bin/env python
"""
rank_results.py — AF2 결과 pLDDT / PAE / ipTM 기반 랭킹
pandas / tabulate 의존성 제거, stdlib csv 만 사용 (proteinmpnn 환경 호환)
"""
import os
import glob
import json
import csv
import argparse


def main():
    parser = argparse.ArgumentParser()
    # [FIX Bug2] UPDD.py 호출 규약에 맞게 arg 이름 수정 (--input_dir → --inputdir)
    parser.add_argument("--inputdir",  required=True,  help="AF2 results dir")
    parser.add_argument("--outputcsv", default=None,   help="순위 결과 저장 CSV 경로 (선택)")
    args = parser.parse_args()

    files = glob.glob(os.path.join(args.inputdir, "**", "*_scores_*.json"), recursive=True)
    results = []

    for f in files:
        try:
            with open(f, "r", encoding="utf-8") as fh:
                data = json.load(fh)

            plddt_list = data.get("plddt", [])
            if not plddt_list:
                plddt_list = data.get("mean_plddt", [])

            if isinstance(plddt_list, list) and plddt_list:
                plddt = sum(plddt_list) / len(plddt_list)
            elif isinstance(plddt_list, (int, float)):
                plddt = float(plddt_list)
            else:
                # [Diagnostics] pLDDT 필드가 list/scalar 어느 쪽으로도 인식되지 않으면
                # ColabFold 출력 포맷이 변경되었을 가능성을 의미한다. 타입 정보를
                # 함께 출력하여 운영자가 출처를 추적할 수 있게 한다.
                print(f"[WARNING] Missing or unrecognized plddt data in {f} "
                      f"(type={type(plddt_list).__name__}), defaulting to 0.0")
                plddt = 0.0

            pae_flat = [v for row in data.get("pae", []) for v in row]
            if pae_flat:
                pae = sum(pae_flat) / len(pae_flat)
            else:
                print(f"[WARNING] Missing pae data in {f}, defaulting to 0.0")
                pae = 0.0
            
            iptm     = data.get("iptm", None)

            name = os.path.basename(f).split("_scores")[0]
            results.append({
                "id":    name,
                "plddt": round(plddt, 3),
                "pae":   round(pae,   3),
                "iptm":  iptm,
            })
        except Exception as e:
            print(f"[!] 파일 파싱 실패: {f} → {e}")
            continue

    if not results:
        print("[!] No results found.")
        return

    # [FIX Bug1] pandas 없이 정렬 + 중복 제거
    results.sort(key=lambda x: x["plddt"], reverse=True)
    seen, deduped = set(), []
    for r in results:
        if r["id"] not in seen:
            seen.add(r["id"])
            deduped.append(r)

    # [FIX Bug3] tabulate 없이 직접 테이블 출력
    header = f"{'id':<40} {'plddt':>7} {'pae':>7} {'iptm':>6}"
    print("\nTOP 5 Candidates")
    print("-" * len(header))
    print(header)
    print("-" * len(header))
    for r in deduped[:5]:
        iptm_str = f"{r['iptm']:.3f}" if r["iptm"] is not None else "  N/A"
        print(f"{r['id']:<40} {r['plddt']:>7.3f} {r['pae']:>7.3f} {iptm_str:>6}")
    print("-" * len(header))

    if args.outputcsv:
        os.makedirs(os.path.dirname(args.outputcsv) or ".", exist_ok=True)
        with open(args.outputcsv, "w", newline="", encoding="utf-8") as csvf:
            writer = csv.DictWriter(csvf, fieldnames=["id", "plddt", "pae", "iptm"])
            writer.writeheader()
            writer.writerows(deduped)
        print(f"[+] 전체 결과 저장됨: {args.outputcsv}")


if __name__ == "__main__":
    main()
