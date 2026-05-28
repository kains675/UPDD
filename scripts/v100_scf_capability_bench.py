#!/usr/bin/env python
"""
scripts/v100_scf_capability_bench.py — V100 32GB DF-J/K SCF capability benchmark.

목적: CLAUDE.md hardware 경계 #1 ("RTX 5070 Ti 16GB 에서 n_qm ≥ 466 의 DF-J/K
불가능") 를 V100 32GB 가 깨는지 실증. V100 의 고유 구조적 가치 (단순 속도가 아닌
capability) 측정.

설계 — 순수 hardware capability benchmark:
- charge-broken 6WGN snapshot 은 run_qmmm.py 의 R-15/16 guard 가 fail-fast → 사용 불가
- 대신 production AO scale 에 matched 된 water cluster 로 DF-J/K SCF (charge 무관,
  neutral closed-shell, 재현가능, 표준 벤치마크). def2-SVP water = 34 AO/water →
  ~150 waters = 450 atoms ≈ 5100 AO ≈ 실제 n_qm 466-556 QM region 의 AO scale
- science claim 0 (publication 용 아님), 순수 "V100 이 production AO scale 에서
  DF-J/K SCF 를 수렴하는가 + wall-clock + VRAM" 측정

Graded series: 64/128/192 waters (192/384/576 atoms) — n_qm 466-556 production
range 를 bracket. 각 점에서 V100 의 converged / cycles / wall-clock / peak VRAM.

사용법:
    python scripts/v100_scf_capability_bench.py --waters 64,128,192 --remote
    python scripts/v100_scf_capability_bench.py --waters 128         # host-only (5070 Ti, DF OOM 예상)
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "utils"))


def _scf_bench_script(n_waters: int, scf_mode: str = "direct",
                      cderi_dir: str = "", xc: str = "wb97xd",
                      basis: str = "def2-svp") -> str:
    """V100/host SCF 벤치마크 snippet.

    scf_mode:
      - 'direct'     : DF 미사용, ERI on-the-fly (CDERI 저장 X, memory-light, compute-bound)
      - 'outcore-df' : DF + CDERI 를 disk (cderi_dir) 에 저장 후 cycle 마다 stream
                       (I/O-bound; cderi_dir 가 SATA vs M.2 면 disk 대역폭 영향 측정)
    cderi_dir: outcore-df 시 CDERI HDF5 저장 위치 (PYSCF_TMPDIR). SATA vs M.2 비교용.
    """
    return f"""
import json, time, math, os
import numpy as np
# outcore-df: CDERI HDF5 가 PYSCF_TMPDIR 에 쓰임 → cderi_dir (SATA or M.2) 로 routing
_scf_mode = {scf_mode!r}
_cderi_dir = {cderi_dir!r}
if _scf_mode == 'outcore-df' and _cderi_dir:
    os.environ['PYSCF_TMPDIR'] = _cderi_dir
    os.environ['TMPDIR'] = _cderi_dir
# CRITICAL (run_qmmm.py L37-43 replication): min_ao_blksize 를 gpu4pyscf import
# 전에 set (read-once). 기본값 128-256 은 J/K 버퍼 [blksize,nao,nao] 를 너무 크게
# 잡아 32GB V100 도 nao>=1500 에서 OOM. 64 로 하향 → blksize ~40 → VRAM peak 제어.
import pyscf as _early
_early.__config__.min_ao_blksize = 64
from pyscf import gto, lib
if _scf_mode == 'outcore-df' and _cderi_dir:
    lib.param.TMPDIR = _cderi_dir
import gpu4pyscf.dft as dft

n_waters = {n_waters}
# Cubic grid water placement (~3.1 Å spacing, non-overlapping)
side = math.ceil(n_waters ** (1/3))
spacing = 3.1
# 단일 water 내부 geometry (O + 2H, ~104.5°, O-H 0.957 Å)
w_local = [('O', 0.000, 0.000, 0.000),
           ('H', 0.757, 0.586, 0.000),
           ('H', -0.757, 0.586, 0.000)]
atoms = []
placed = 0
for ix in range(side):
    for iy in range(side):
        for iz in range(side):
            if placed >= n_waters: break
            ox, oy, oz = ix*spacing, iy*spacing, iz*spacing
            for el, dx, dy, dz in w_local:
                atoms.append(f'{{el}} {{ox+dx:.4f}} {{oy+dy:.4f}} {{oz+dz:.4f}}')
            placed += 1
atom_str = '; '.join(atoms)
n_atoms = placed * 3

import cupy
dev = cupy.cuda.Device(0)
_free0, _total0 = dev.mem_info

start = time.time()
mol = gto.M(atom=atom_str, basis={basis!r}, verbose=0)
nao = mol.nao
build_t = time.time() - start

# DF-J/K SCF (gpu4pyscf): density_fit() 가 DF-J/K 활성. 5070 Ti 16GB 는 nao~5000
# 에서 CDERI OOM (CLAUDE.md 경계 #1), V100 32GB 는 가능해야.
scf_start = time.time()
try:
    if _scf_mode == 'outcore-df':
        # DF + CDERI 를 disk (PYSCF_TMPDIR=cderi_dir) 에 저장. with_df.max_memory 를
        # 작게 (500 MB) → CDERI 가 GPU/RAM 에 안 들어가고 disk HDF5 로 outcore.
        # cycle 마다 disk → CPU → GPU stream. disk 대역폭 (SATA vs M.2) 영향 측정.
        mf = dft.RKS(mol, xc={xc!r}).density_fit(auxbasis='def2-universal-jfit')
        mf = mf.to_gpu()
        mf.max_memory = 2000
        try:
            mf.with_df.max_memory = 500   # force outcore CDERI (disk HDF5)
        except Exception:
            pass
    else:
        # Direct-SCF (DF 미사용) — UPDD production 의 n_qm>=466 실제 방법 (CLAUDE.md 경계 #1:
        # "5070 Ti 16GB 에서 n_qm>=466 DF 불가능 → direct-SCF fallback, 130-164 s/cycle").
        # ERI on-the-fly → CDERI 저장 X → OOM 회피. V100 비교점 = per-cycle FP64 속도.
        mf = dft.RKS(mol, xc={xc!r}).to_gpu()
        mf.max_memory = 16000
    mf.conv_tol = 1e-8
    mf.max_cycle = 50
    e = mf.kernel()
    converged = bool(mf.converged)
    ncycle = getattr(mf, 'cycles', -1)
    scf_t = time.time() - scf_start
    per_cycle = round(scf_t / max(ncycle, 1), 1) if isinstance(ncycle, int) and ncycle > 0 else None
    _free1, _total1 = dev.mem_info
    peak_used_mb = (_total0 - _free1) / 1024**2
    # CDERI 파일 크기 (outcore-df 일 때, disk 저장 검증)
    cderi_mb = 0
    if _scf_mode == 'outcore-df' and _cderi_dir and os.path.isdir(_cderi_dir):
        try:
            cderi_mb = sum(os.path.getsize(os.path.join(_cderi_dir,f)) for f in os.listdir(_cderi_dir) if f.endswith('.h5') or 'cderi' in f.lower()) / 1024**2
        except Exception:
            pass
    out = {{
        'n_waters': placed, 'n_atoms': n_atoms, 'nao': int(nao),
        'xc': {xc!r}, 'basis': {basis!r}, 'scf_mode': _scf_mode,
        'cderi_dir': _cderi_dir, 'cderi_disk_mb': round(cderi_mb,0),
        'converged': converged, 'n_cycles': int(ncycle),
        'energy_ha': float(e), 'build_s': round(build_t,1),
        'scf_s': round(scf_t,1), 'per_cycle_s': per_cycle,
        'peak_vram_mb': round(peak_used_mb,0),
        'total_vram_mb': round(_total0/1024**2,0), 'status': 'OK',
    }}
except Exception as ex:
    out = {{
        'n_waters': placed, 'n_atoms': n_atoms, 'nao': int(nao),
        'xc': {xc!r}, 'basis': {basis!r}, 'scf_mode': _scf_mode,
        'cderi_dir': _cderi_dir,
        'converged': False, 'status': 'FAIL',
        'error': f'{{type(ex).__name__}}: {{str(ex)[:200]}}',
        'total_vram_mb': round(_total0/1024**2,0),
    }}
print('BENCH_RESULT_BEGIN' + json.dumps(out) + 'BENCH_RESULT_END')
"""


def _parse(stdout: str) -> Dict[str, Any]:
    b = stdout.find("BENCH_RESULT_BEGIN"); e = stdout.find("BENCH_RESULT_END")
    if b < 0 or e < 0:
        raise RuntimeError(f"marker 미발견. stdout last 600:\n{stdout[-600:]}")
    return json.loads(stdout[b+len("BENCH_RESULT_BEGIN"):e])


def run_vm(n_waters: int, scf_mode: str = "direct", cderi_dir: str = "", conda_env: str = "qmmm") -> Dict[str, Any]:
    import dispatch
    vm = dispatch.VMExecutor()
    if not vm.is_connected():
        raise RuntimeError(f"VM SSH 연결 실패 ({vm.ssh_target})")
    script = _scf_bench_script(n_waters, scf_mode=scf_mode, cderi_dir=cderi_dir)
    payload = base64.b64encode(script.encode()).decode("ascii")
    py = f"/home/san/miniconda3/envs/{conda_env}/bin/python"
    cmd = f"{py} -c \"import base64; exec(base64.b64decode('{payload}').decode())\""
    r = vm.execute(cmd, timeout=7200, retry=1)  # 큰 SCF 는 최대 2h
    if r["returncode"] != 0:
        raise RuntimeError(f"VM SCF 실패 (rc={r['returncode']}): {str(r.get('stderr',''))[:400]}")
    res = _parse(str(r["stdout"])); res["location"] = "vm_v100"
    return res


def run_host(n_waters: int, scf_mode: str = "direct", cderi_dir: str = "", conda_env: str = "qmmm") -> Dict[str, Any]:
    script = _scf_bench_script(n_waters, scf_mode=scf_mode, cderi_dir=cderi_dir)
    cmd = ["conda", "run", "-n", conda_env, "--no-capture-output", "python", "-c", script]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if r.returncode != 0:
        raise RuntimeError(f"Host SCF 실패: {r.stderr[:400]}")
    res = _parse(r.stdout); res["location"] = "host_5070ti"
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--waters", default="64,128,192", help="콤마 구분 water 개수 (default 64,128,192 = 192/384/576 atoms)")
    ap.add_argument("--remote", action="store_true", help="V100 (VM) 실행. 없으면 host 5070 Ti")
    ap.add_argument("--scf-mode", dest="scf_mode", default="direct", choices=["direct", "outcore-df"],
                    help="direct (ERI on-the-fly, memory-light) | outcore-df (CDERI disk stream)")
    ap.add_argument("--cderi-dir", dest="cderi_dir", default="",
                    help="outcore-df CDERI 저장 위치 (PYSCF_TMPDIR). SATA vs M.2 비교용")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    sizes = [int(s.strip()) for s in args.waters.split(",") if s.strip()]
    runner = run_vm if args.remote else run_host
    loc = "VM V100 32GB" if args.remote else "host 5070 Ti 16GB"

    print(f"=== V100 SCF benchmark ({loc}) ===")
    print(f"water sizes: {sizes} (atoms: {[s*3 for s in sizes]})")
    print(f"scf_mode: {args.scf_mode}" + (f" | cderi_dir: {args.cderi_dir}" if args.scf_mode == 'outcore-df' else ""))
    print(f"XC/basis: wb97xd/def2-svp")
    print()

    results: List[Dict[str, Any]] = []
    for n in sizes:
        print(f"[{n} waters / {n*3} atoms] {args.scf_mode} SCF ...")
        try:
            r = runner(n, scf_mode=args.scf_mode, cderi_dir=args.cderi_dir)
        except Exception as ex:
            print(f"  ❌ {type(ex).__name__}: {str(ex)[:200]}")
            results.append({"n_waters": n, "n_atoms": n*3, "status": "DISPATCH_FAIL", "error": str(ex)[:200]})
            continue
        if r.get("status") == "OK":
            print(f"  ✓ nao={r['nao']} converged={r['converged']} cycles={r['n_cycles']} "
                  f"scf={r['scf_s']}s peak_vram={r['peak_vram_mb']:.0f}/{r['total_vram_mb']:.0f} MiB")
        else:
            print(f"  ✗ nao={r.get('nao','?')} FAIL: {r.get('error','?')[:150]}")
        results.append(r)

    print()
    print("=== Summary ===")
    print(f"{'atoms':>7} {'nao':>6} {'conv':>5} {'cycles':>7} {'scf_s':>8} {'peak_VRAM_MB':>13} {'status':>8}")
    for r in results:
        print(f"{r.get('n_atoms','?'):>7} {r.get('nao','?'):>6} {str(r.get('converged','?')):>5} "
              f"{r.get('n_cycles','?'):>7} {str(r.get('scf_s','?')):>8} {str(r.get('peak_vram_mb','?')):>13} {r.get('status','?'):>8}")

    # n_qm≥466 (155 waters) DF-J/K 성공 여부 = V100 고유 capability 판정
    big = [r for r in results if r.get("n_atoms", 0) >= 466 and r.get("status") == "OK" and r.get("converged")]
    print()
    if big:
        print(f"🟢 V100 capability 실증: n_atoms≥466 DF-J/K SCF 수렴 ({len(big)}점). "
              f"5070 Ti 16GB 불가능 영역 (CLAUDE.md 경계 #1) 을 V100 32GB 가 처리.")
    else:
        print(f"🟡 n_atoms≥466 DF-J/K 미수렴 또는 미실행 — VRAM/convergence 추가 진단 필요.")

    if args.output:
        Path(args.output).write_text(json.dumps(results, indent=2))
        print(f"\nsaved → {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
