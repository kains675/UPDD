# CLAUDE.md — execute_qmmm 전면 개선 (Adaptive + 태그 + 필터 + Resume)

## 범위
- `UPDD.py`의 `execute_qmmm` 함수 (line 1091~1163) 교체.
- `utils/run_qmmm.py`는 수정하지 않는다 (R-4).
- UPDATE.md 최상단 기록. R-1/R-2/R-3/R-4 준수.

## ⚠️ 안전 작업 지침
1. SMILES 검색 금지. 2. `view`로 줄 확인 후 수정. 3. pixi≠conda.

---

## 현재 문제 4가지

1. **OOM 연쇄 폭발**: 4 병렬 시 VRAM peak 초과 → 한 디자인 OOM → 나머지 연쇄 실패
2. **스냅샷 없는 디자인 순회**: 20개 전체를 순회하며 12개에 대해 amber_charges 로드 → "PDB 없음" 반복
3. **병렬 출력 뒤섞임**: 어느 디자인의 SCF인지 구별 불가
4. **Resume 없음**: 중단 후 재실행 시 이미 완료된 디자인도 재계산

## 수정: execute_qmmm 함수 전체 교체

```bash
# 현재 함수 위치 확인
grep -n "def execute_qmmm" ~/UPDD_proj/UPDD.py
# 기대: line 1091
```

**line 1091의 `def execute_qmmm(...)` 부터 line 1163의 `return qmmm_out` 까지 전체를 아래로 교체한다.**

⚠️ 함수 시그니처는 동일하게 유지한다 (호출부 수정 불필요):

```python
def execute_qmmm(snap_dir: str, qmmm_mode: str, ncaa_element: str, af2_out_dir: str,
                 parallel_workers: int, target_out_dir: str, binder_chain: str = "B") -> str:
    print_step("Step 11: QM/MM — wB97X-D / PySCF")
    qmmm_out = os.path.join(target_out_dir, "qmmm_results")
    os.makedirs(qmmm_out, exist_ok=True)
    elem_arg = ncaa_element if ncaa_element else "none"

    # ── Phase 1: 스냅샷 존재 디자인만 필터링 ──────────────────
    af2_ranking_path = os.path.join(af2_out_dir, "af2_ranking.csv")
    if not os.path.exists(af2_ranking_path):
        print("  [!] af2_ranking.csv 없음 — QM/MM 건너뜀")
        return qmmm_out

    ranking_data: List[Dict[str, Any]] = []
    with open(af2_ranking_path, 'r', encoding="utf-8") as f:
        for row in csv.DictReader(f):
            ranking_data.append({'id': row['id'], 'plddt': float(row['plddt'])})
    ranking_data.sort(key=lambda x: x['plddt'], reverse=True)

    valid_ids = []
    skipped_ids = []
    for item in ranking_data:
        snaps = glob.glob(os.path.join(snap_dir, f"{item['id']}*.pdb"))
        if snaps:
            valid_ids.append(item['id'])
        else:
            skipped_ids.append(item['id'])

    if skipped_ids:
        print(f"  [Filter] 스냅샷 없는 디자인 {len(skipped_ids)}개 제외")
    if not valid_ids:
        print("  [!] 유효한 스냅샷이 없습니다 — QM/MM 건너뜀")
        return qmmm_out

    # ── Phase 2: Resume — 이미 정상 완료된 디자인 스킵 ────────
    def _is_qmmm_done(design_id: str) -> bool:
        """해당 디자인의 모든 스냅샷 QM/MM이 정상 완료되었는지 확인."""
        snap_pdbs = sorted(glob.glob(os.path.join(snap_dir, f"{design_id}*.pdb")))
        if not snap_pdbs:
            return False
        for snap_pdb in snap_pdbs:
            stem = os.path.splitext(os.path.basename(snap_pdb))[0]
            json_path = os.path.join(qmmm_out, f"{stem}_qmmm_{qmmm_mode}.json")
            if not os.path.exists(json_path):
                return False
            try:
                with open(json_path, "r", encoding="utf-8") as jf:
                    data = json.load(jf)
                e_val = data.get("e_qmmm_hartree", data.get("e_qm_hartree", 0))
                if not e_val or e_val == 0:
                    return False  # OOM/에러로 생성된 빈 JSON
            except (json.JSONDecodeError, KeyError, TypeError):
                return False
        return True

    already_done = []
    need_calc = []
    for d_id in valid_ids:
        if _is_qmmm_done(d_id):
            already_done.append(d_id)
        else:
            need_calc.append(d_id)

    if already_done:
        print(f"  [Resume] 이미 완료된 디자인 {len(already_done)}개 스킵")
    if not need_calc:
        print(f"  [Resume] 모든 디자인 완료 — QM/MM 재계산 불필요")
    else:
        print(f"  [우선순위] {len(need_calc)}개 디자인 계산 예정 "
              f"(전체 {len(valid_ids)}개 중 {len(already_done)}개 완료)")

        # ── QM/MM 워커 정의 (태그 + OOM 감지) ────────────────
        base_cmd = [
            "conda", "run", "-n", ENV_NAMES.get("qmmm", "qmmm"),
            "--no-capture-output", "python", SCRIPT_PATHS["run_qmmm"],
        ]

        def qmmm_worker(design_id: str) -> dict:
            qmmm_args = [
                "--snapdir",      snap_dir,
                "--outputdir",    qmmm_out,
                "--ncaa_elem",    elem_arg,
                "--mode",         qmmm_mode,
                "--qm_xc",        "wb97xd",
                "--filter",       design_id,
                "--binder_chain", binder_chain,
            ]
            cmd = base_cmd + qmmm_args

            # 짧은 태그: "design_w4_1_s2" → "w4_1_s2"
            tag = design_id.replace("design_", "").split("_unrelaxed")[0]
            if len(tag) > 10:
                tag = tag[:10]

            try:
                proc = subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, bufsize=1
                )
                oom_detected = False
                for line in proc.stdout:
                    stripped = line.rstrip()
                    if stripped:
                        print(f"  [{tag}] {stripped}")
                        if "out of memory" in stripped.lower() or "cudaErrorMemoryAllocation" in stripped:
                            oom_detected = True
                proc.wait()

                if oom_detected:
                    return {"id": design_id, "status": "OOM",
                            "msg": f"CUDA OOM (returncode={proc.returncode})"}
                elif proc.returncode != 0:
                    return {"id": design_id, "status": "ERROR",
                            "msg": f"비정상 종료 (returncode={proc.returncode})"}
                else:
                    return {"id": design_id, "status": "SUCCESS", "msg": ""}
            except Exception as e:
                return {"id": design_id, "status": "ERROR", "msg": str(e)}

        # ── Phase 3: Adaptive 병렬 실행 (4→3→2→1) ────────────
        pending_ids = list(need_calc)
        current_workers = min(parallel_workers, len(pending_ids))
        all_success_ids = list(already_done)
        max_passes = 4

        for pass_num in range(1, max_passes + 1):
            if not pending_ids:
                break

            print(f"\n  🔥 [QM/MM Pass {pass_num}] {len(pending_ids)}개 디자인, "
                  f"{current_workers}개 GPU 워커")

            pass_results = []
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=current_workers) as executor:
                    future_map = {
                        executor.submit(qmmm_worker, d_id): d_id
                        for d_id in pending_ids
                    }
                    for future in concurrent.futures.as_completed(future_map):
                        result = future.result()
                        pass_results.append(result)
                        if result["status"] == "SUCCESS":
                            print(f"  ✅ {result['id']} 완료")
                        elif result["status"] == "OOM":
                            print(f"  ⚠️  {result['id']} OOM — 다음 pass에서 재시도")
                        else:
                            print(f"  ❌ {result['id']} 에러: {result['msg']}")
            except KeyboardInterrupt:
                print(f"\n  [!] 🛑 수동 개입! 병렬 작업을 강제 중단합니다.\n")
                break

            succeeded = [r for r in pass_results if r["status"] == "SUCCESS"]
            oom_failed = [r for r in pass_results if r["status"] == "OOM"]

            all_success_ids.extend([r["id"] for r in succeeded])

            if oom_failed:
                pending_ids = [r["id"] for r in oom_failed]
                # OOM 실패 디자인의 기존 결과 JSON 정리 (재계산 위해)
                for r in oom_failed:
                    for stale in glob.glob(os.path.join(qmmm_out, f"{r['id']}*_qmmm_*.json")):
                        os.remove(stale)
                new_workers = max(1, current_workers - 1)
                print(f"\n  [Adaptive] OOM {len(oom_failed)}개 감지. "
                      f"워커 축소: {current_workers} → {new_workers}")
                current_workers = new_workers
            else:
                pending_ids = []

        if pending_ids:
            print(f"  [!] {len(pending_ids)}개 디자인이 단일 워커에서도 OOM. "
                  f"VRAM 부족으로 건너뜁니다.")

        total = len(valid_ids)
        n_success = len(all_success_ids)
        n_oom = len(pending_ids)
        n_error = total - n_success - n_oom
        print(f"\n  [QM/MM 요약] 총 {total}개 | "
              f"✅ 성공 {n_success} | ⚠️ OOM {n_oom} | ❌ 에러 {n_error}")

    # ── [v55 유지] 개별 JSON → 통합 summary 병합 ─────────────
    all_results = []
    failed = []
    for jf in sorted(glob.glob(os.path.join(qmmm_out, f"*_qmmm_{qmmm_mode}.json"))):
        try:
            with open(jf, "r", encoding="utf-8") as f:
                result = json.load(f)
            all_results.append(result)
        except (json.JSONDecodeError, IOError) as _je:
            failed.append(os.path.basename(jf))
    summary_path = os.path.join(qmmm_out, "qmmm_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "method":    f"wb97xd/6-31G*",
            "n_calc":    len(all_results),
            "n_failed":  len(failed),
            "failed":    failed,
            "results":   all_results,
        }, f, indent=2)
    converged = [r for r in all_results if r.get("converged")]
    print(f"  [QM/MM 통합] {len(converged)}/{len(all_results)} 수렴, 요약: {summary_path}")

    return qmmm_out
```

---

## 검증

```bash
# 1. 문법 확인
python -c "import py_compile; py_compile.compile('UPDD.py', doraise=True)" && echo OK

# 2. 태그 출력 확인
grep "\[w[0-9]" ~/UPDD_proj/log/6WGN/updd_*.log | tail -5
# 기대: [w4_1_s2] cycle= ... 형태

# 3. 스냅샷 필터링 확인
grep "스냅샷 없는 디자인" ~/UPDD_proj/log/6WGN/updd_*.log
# 기대: "스냅샷 없는 디자인 12개 제외" (또는 해당 수)

# 4. Resume 확인
grep "Resume\|이미 완료" ~/UPDD_proj/log/6WGN/updd_*.log
# 기대: 이전 완료 디자인 스킵 메시지

# 5. Adaptive 확인
grep "Adaptive\|워커 축소\|Pass\|QM/MM 요약" ~/UPDD_proj/log/6WGN/updd_*.log

# 6. v55 통합 summary 확인
cat ~/UPDD_proj/outputs/6WGN_*/qmmm_results/qmmm_summary.json | python -m json.tool | head -10
```

---

## ⚠️ 절대 규칙
1. SciVal 거부 시 코드 배포 금지.
2. Runner 검증 없이 완료 선언 금지.
3. 단위 혼동 즉시 중단.
4. pixi≠conda.
5. UPDATE.md 기록 필수.
6. R-4: `run_qmmm.py`는 수정하지 않는다.
7. **캐시 무효화**: 수정 후 기존 qmmm_results/ 삭제 → 재실행.
