#!/usr/bin/env python
"""
prepare_af2_input.py — MPNN FASTA → ColabFold input.csv 변환
pandas 의존성 제거, stdlib csv 만 사용 (proteinmpnn 환경 호환)
"""
import os
import glob
import csv
import argparse


def _parse_fasta(filepath):
    """[v3 11-1] 멀티라인 FASTA 파서 (제너레이터).

    표준 FASTA 규약에 따라 ">" 로 시작하는 헤더 다음에 여러 줄에 걸쳐 시퀀스가
    이어질 수 있도록 처리한다. 기존 i/i+1 페어 방식은 시퀀스가 단일 줄이라는
    가정에 의존하여, MPNN 출력이 줄바꿈을 포함할 경우 시퀀스가 절단되는 잠재
    버그가 있었다.

    Yields:
        Tuple[str, str]: (header, joined_sequence) — header 는 ">" 를 포함한
        원본 헤더, sequence 는 모든 시퀀스 줄을 공백 없이 이어붙인 결과.
    """
    header, seq_parts = None, []
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header, seq_parts = line, []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def main():
    parser = argparse.ArgumentParser()
    # [FIX Bug2] UPDD.py 호출 규약에 맞게 arg 이름 수정 (--input_dir → --inputdir)
    parser.add_argument("--inputdir",  required=True, help="MPNN seqs dir")
    parser.add_argument("--outputcsv", required=True, help="Output CSV path")
    args = parser.parse_args()

    fa_files    = glob.glob(os.path.join(args.inputdir, "**", "*.fa"),    recursive=True)
    fasta_files = glob.glob(os.path.join(args.inputdir, "**", "*.fasta"), recursive=True)
    files = fa_files + fasta_files

    if not files:
        print(f"[!] FASTA 파일을 찾을 수 없습니다: {args.inputdir}")
        return

    candidates = []

    for f in files:
        base_id = os.path.basename(f).rsplit(".", 1)[0]
        for sample_idx, (_header, seq_line) in enumerate(_parse_fasta(f)):
            design_id = f"{base_id}_s{sample_idx}"
            if "/" in seq_line:
                seqs       = seq_line.split("/")
                full_input = f"{seqs[0]}:{seqs[1]}"
            else:
                full_input = seq_line
            candidates.append({"id": design_id, "sequence": full_input})

    if not candidates:
        print("[!] 파싱 가능한 시퀀스가 없습니다.")
        return

    # [FIX Bug1] pandas → stdlib csv
    os.makedirs(os.path.dirname(args.outputcsv) or ".", exist_ok=True)
    with open(args.outputcsv, "w", newline="", encoding="utf-8") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=["id", "sequence"])
        writer.writeheader()
        writer.writerows(candidates)

    print(f"[+] ColabFold input 생성 완료: {args.outputcsv} ({len(candidates)} sequences)")


if __name__ == "__main__":
    main()
