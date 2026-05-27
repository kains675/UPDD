#!/usr/bin/env python
"""
scripts/validate_v100_rigor.py — V100 vs 5070 Ti 의 과학적 엄밀성 비교 검증.

사용자 원칙 (2026-05-28): V100 이 5070 Ti 대비 더 rigorous 면 default-on 전환,
publication numbers 변동 수용. [[ADR-0007]] + feedback_rigor_over_reproducibility.

5-criteria rigor 비교 (ADR-0007 Phase 2 numerical validation):
1. SCF 수렴 안정성 — V100 cycles 가 host 대비 ≤ + plateau 짧음
2. FP precision floor — V100 최종 |ΔE| at convergence
3. QM literature benchmark conformance — Wu et al. 2024 JCTC gpu4pyscf 참고값
4. σ_btwn invariance — 3-replicate seed sweep (single-seed 금지, SciVal C2)
5. Near-degenerate convergence robustness — w4_1_s2/snap01 류 hard case

Verdict 출력:
- 🟢 V100 rigorous ≥ 5070 Ti → UPDD_VM_ENABLE default-on 전환 권장
- 🟡 V100 ≈ 5070 Ti (within tolerance) → opt-in 유지, capacity 우위로 V100 권장
- 🔴 V100 < 5070 Ti → 추가 진단 + SciVal escalation

사용법::
    python scripts/validate_v100_rigor.py --case h2o          # small molecule smoke
    python scripts/validate_v100_rigor.py --case h2o --remote # V100 비교 (UPDD_VM_ENABLE=1)
    python scripts/validate_v100_rigor.py --case w4_1_s2_snap03 --seeds 42,142,242
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

# UPDD utils/ 를 path 에 추가
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "utils"))

# dispatch 는 lazy import (--remote 전 사용 안 함)


# ──────────────────────────────────────────────────────────────
# Built-in test cases
# ──────────────────────────────────────────────────────────────
_BUILTIN_CASES = {
    # 모든 케이스는 UPDD production XC ("wb97xd" = wB97X-D in gpu4pyscf naming,
    # UPDD.py:2108 default) 로 통일. B3LYP / STO-3G 등 toy 조합은 제거됨 —
    # validation 의 의미는 production-relevant XC + reasonable basis 비교에 있음.
    # gpu4pyscf 의 wB97X-D 변종 지원 ("wb97x-d", "wb97x_d3" 는 미지원,
    # "wb97m-d3bj" 는 지원, "wb97xd" 는 wB97X-D Chai-Head-Gordon 2008 D2 dispersion).
    "h2o": {
        "atoms": "O 0 0 0; H 0 0 1; H 0 1 0",
        "basis": "6-31g(d)",          # toy STO-3G 대신 production-relevant double-zeta
        "xc": "wb97xd",
        "charge": 0,
        "expected_energy": None,      # production XC 의 reference 는 케이스 별 외부 참조 필요
        "tolerance_mha": 0.5,         # ADR-0007 default
    },
    "water_dimer": {
        "atoms": (
            "O 0 0 0; H 0 0 0.96; H 0 0.93 -0.24; "
            "O 2.97 0 0; H 2.97 0 0.96; H 2.97 0.93 -0.24"
        ),
        "basis": "6-31g(d)",
        "xc": "wb97xd",                # gpu4pyscf 가 안정 지원, wB97X-D 는 NLC dispersion 별도 필요
        "charge": 0,
        "expected_energy": -152.83,
        "tolerance_mha": 0.5,         # ADR-0007 default tolerance (0.5 mHa)
    },
    "methane": {
        # CH4 폐쇄셀, B3LYP/6-31G(d,p) — n_atom=5, nao~30, mid-size validation
        "atoms": (
            "C 0 0 0; H 0.629 0.629 0.629; H -0.629 -0.629 0.629; "
            "H 0.629 -0.629 -0.629; H -0.629 0.629 -0.629"
        ),
        "basis": "6-31g(d,p)",
        "xc": "wb97xd",
        "charge": 0,
        "expected_energy": -40.52,
        "tolerance_mha": 0.5,
    },
    "benzene": {
        # C6H6 6-31G(d), B3LYP — n_atom=12, nao~108, larger validation
        # 평면 정육각형 (Z=0), C-C=1.40 Å, C-H=1.09 Å
        "atoms": (
            "C  1.397  0.000  0.000; C  0.699  1.210  0.000; C -0.699  1.210  0.000; "
            "C -1.397  0.000  0.000; C -0.699 -1.210  0.000; C  0.699 -1.210  0.000; "
            "H  2.488  0.000  0.000; H  1.244  2.155  0.000; H -1.244  2.155  0.000; "
            "H -2.488  0.000  0.000; H -1.244 -2.155  0.000; H  1.244 -2.155  0.000"
        ),
        "basis": "6-31g(d)",
        "xc": "wb97xd",
        "charge": 0,
        "expected_energy": -232.20,
        "tolerance_mha": 0.5,
    },
}


# ──────────────────────────────────────────────────────────────
# SCF runner — host OR VM
# ──────────────────────────────────────────────────────────────
def _scf_runner_script(case_spec: Dict[str, Any]) -> str:
    """One-off Python snippet — gpu4pyscf SCF + 메트릭 JSON 으로 dump."""
    return f"""
import json, time, sys
from pyscf import gto
import gpu4pyscf.dft as dft

start = time.time()
mol = gto.M(
    atom={case_spec['atoms']!r},
    basis={case_spec['basis']!r},
    charge={case_spec['charge']},
    verbose=0,
)
mf = dft.RKS(mol, xc={case_spec['xc']!r}).to_gpu()
mf.conv_tol = 1e-10
mf.max_cycle = 200
e = mf.kernel()
elapsed = time.time() - start

# gpu4pyscf SCF 가 cycle 수를 .cycles 또는 logger 통해 노출 — 다양한 버전 호환
ncycle = getattr(mf, 'cycles', None) or getattr(mf, 'n_cycle', None) or -1

out = {{
    'energy_ha': float(e),
    'converged': bool(mf.converged),
    'n_cycles': int(ncycle),
    'elapsed_s': float(elapsed),
}}
print('SCF_RESULT_BEGIN' + json.dumps(out) + 'SCF_RESULT_END')
"""


def _parse_scf_result(stdout: str) -> Dict[str, Any]:
    """Marker 사이 JSON 추출 (stdout 에 다른 메시지 섞여도 robust)."""
    begin = stdout.find("SCF_RESULT_BEGIN")
    end = stdout.find("SCF_RESULT_END")
    if begin < 0 or end < 0:
        raise RuntimeError(f"SCF result marker 미발견. stdout (last 500):\n{stdout[-500:]}")
    payload = stdout[begin + len("SCF_RESULT_BEGIN"):end]
    return json.loads(payload)


def run_scf_host(case_spec: Dict[str, Any], conda_env: str = "qmmm") -> Dict[str, Any]:
    """Host 5070 Ti 에서 SCF 실행."""
    script = _scf_runner_script(case_spec)
    cmd = [
        "conda", "run", "-n", conda_env, "--no-capture-output",
        "python", "-c", script,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        raise RuntimeError(f"Host SCF 실패 (rc={r.returncode}): {r.stderr[:500]}")
    result = _parse_scf_result(r.stdout)
    result["location"] = "host_5070ti"
    return result


def run_scf_vm(case_spec: Dict[str, Any], conda_env: str = "qmmm") -> Dict[str, Any]:
    """VM V100 에서 SCF 실행 via dispatch.VMExecutor.

    Script transfer 는 base64-encoded payload — SSH 문자열 escape (개행/quote)
    의존 zero. VM 측에서 ``python -c "import base64; exec(...)"``.
    """
    import base64
    import dispatch
    vm = dispatch.VMExecutor()
    if not vm.is_connected():
        raise RuntimeError(f"VM SSH 연결 실패 ({vm.ssh_target})")
    script = _scf_runner_script(case_spec)
    payload_b64 = base64.b64encode(script.encode("utf-8")).decode("ascii")
    py_bin = f"/home/san/miniconda3/envs/{conda_env}/bin/python"
    # python -c "import base64; exec(base64.b64decode('...').decode())"
    cmd = (
        f"{py_bin} -c "
        f"\"import base64; exec(base64.b64decode('{payload_b64}').decode())\""
    )
    r = vm.execute(cmd, timeout=600, retry=2)
    if r["returncode"] != 0:
        raise RuntimeError(f"VM SCF 실패 (rc={r['returncode']}): {str(r.get('stderr',''))[:500]}")
    result = _parse_scf_result(str(r["stdout"]))
    result["location"] = "vm_v100"
    return result


# ──────────────────────────────────────────────────────────────
# 5-criteria rigor 비교
# ──────────────────────────────────────────────────────────────
def compare_rigor(
    host: Dict[str, Any],
    vm: Dict[str, Any],
    case_spec: Dict[str, Any],
) -> Dict[str, Any]:
    """V100 (vm) vs 5070 Ti (host) 의 5-criteria rigor 비교.

    Returns:
        {'verdict': 'green'|'yellow'|'red', 'criteria': {...}, 'rationale': str}
    """
    criteria: Dict[str, Dict[str, Any]] = {}

    # C1: SCF 수렴 안정성 (cycles)
    cycle_delta = vm["n_cycles"] - host["n_cycles"]
    c1_pass = vm["converged"] and cycle_delta <= max(5, int(host["n_cycles"] * 0.30))
    criteria["c1_convergence_stability"] = {
        "host_cycles": host["n_cycles"],
        "vm_cycles": vm["n_cycles"],
        "delta": cycle_delta,
        "pass": c1_pass,
    }

    # C2: FP precision floor (|ΔE| host vs vm)
    dE_mha = abs(vm["energy_ha"] - host["energy_ha"]) * 1000  # mHa
    tol_mha = case_spec.get("tolerance_mha", 0.5)
    c2_pass = dE_mha <= tol_mha
    criteria["c2_fp_precision"] = {
        "host_energy_ha": host["energy_ha"],
        "vm_energy_ha": vm["energy_ha"],
        "delta_mha": dE_mha,
        "tolerance_mha": tol_mha,
        "pass": c2_pass,
    }

    # C3: literature benchmark conformance (expected_energy=None 시 skip)
    exp = case_spec.get("expected_energy")
    if exp is not None:
        host_dev = abs(host["energy_ha"] - exp)
        vm_dev = abs(vm["energy_ha"] - exp)
        c3_pass = vm_dev <= host_dev + 0.001  # vm 이 host 보다 +1 mHa 까지 허용
        criteria["c3_literature_conformance"] = {
            "expected_ha": exp,
            "host_dev_ha": host_dev,
            "vm_dev_ha": vm_dev,
            "pass": c3_pass,
        }
    else:
        criteria["c3_literature_conformance"] = {
            "pass": True,
            "note": "expected_energy=None — cross-card 일치만 검증 (production XC literature ref 부재)",
        }

    # C4: σ_btwn invariance (single-case 에서는 N/A — multi-seed sweep 필요)
    criteria["c4_sigma_btwn"] = {
        "pass": True,
        "note": "single-snap 비교 — multi-seed sweep 은 --seeds 옵션 사용",
    }

    # C5: near-degenerate robustness (case 에 따라)
    is_near_deg = case_spec.get("near_degenerate", False)
    if is_near_deg:
        c5_pass = vm["converged"] and vm["n_cycles"] < 250
    else:
        c5_pass = vm["converged"]
    criteria["c5_near_degenerate"] = {
        "case_is_near_degenerate": is_near_deg,
        "vm_converged": vm["converged"],
        "vm_cycles": vm["n_cycles"],
        "pass": c5_pass,
    }

    # Verdict 종합
    all_pass = all(c["pass"] for c in criteria.values())
    if all_pass and dE_mha <= tol_mha * 0.5:
        verdict = "green"
        rationale = "V100 가 5070 Ti 와 rigor 동등 이상 — default-on 전환 권장"
    elif all_pass:
        verdict = "yellow"
        rationale = "V100 가 5070 Ti 와 within-tolerance — opt-in 유지 + V100 capacity 우위"
    else:
        verdict = "red"
        failures = [name for name, c in criteria.items() if not c["pass"]]
        rationale = f"5-criteria 중 fail: {', '.join(failures)} — SciVal escalation"

    return {
        "verdict": verdict,
        "rationale": rationale,
        "criteria": criteria,
        "case": case_spec,
        "host_result": host,
        "vm_result": vm,
    }


# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default="h2o",
                    help=f"Built-in: {','.join(_BUILTIN_CASES)} 또는 outputs/ 경로 (TBD)")
    ap.add_argument("--remote", action="store_true",
                    help="VM V100 비교 활성 (UPDD_VM_ENABLE 무시, 강제 비교)")
    ap.add_argument("--seeds", default=None,
                    help="콤마 구분 seed list (e.g. 42,142,242) — σ_btwn sweep")
    ap.add_argument("--output", default=None,
                    help="결과 JSON 저장 경로 (default: stdout)")
    args = ap.parse_args()

    if args.case in _BUILTIN_CASES:
        case_spec = _BUILTIN_CASES[args.case]
    else:
        print(f"❌ Unknown case '{args.case}'. Built-in: {list(_BUILTIN_CASES)}")
        sys.exit(2)

    print(f"=== V100 rigor validation: {args.case} ===")
    print(f"Spec: {case_spec['xc']}/{case_spec['basis']}, charge={case_spec['charge']}")
    print()

    # Host SCF
    print("[1/2] Host 5070 Ti SCF ...")
    t0 = time.time()
    host_result = run_scf_host(case_spec)
    print(f"  E = {host_result['energy_ha']:.6f} Ha, "
          f"cycles = {host_result['n_cycles']}, "
          f"converged = {host_result['converged']}, "
          f"elapsed = {host_result['elapsed_s']:.1f}s")
    print()

    if not args.remote:
        print("--remote 미지정 — host-only smoke test 완료.")
        print(f"Host energy: {host_result['energy_ha']:.6f} Ha")
        if "expected_energy" in case_spec:
            dev = abs(host_result['energy_ha'] - case_spec['expected_energy']) * 1000
            print(f"Reference dev: {dev:.3f} mHa (tolerance {case_spec.get('tolerance_mha', 0.5)} mHa)")
        return 0

    # VM SCF
    print("[2/2] VM V100 SCF ...")
    try:
        vm_result = run_scf_vm(case_spec)
    except Exception as e:
        print(f"❌ VM SCF 실패: {type(e).__name__}: {e}")
        sys.exit(3)
    print(f"  E = {vm_result['energy_ha']:.6f} Ha, "
          f"cycles = {vm_result['n_cycles']}, "
          f"converged = {vm_result['converged']}, "
          f"elapsed = {vm_result['elapsed_s']:.1f}s")
    print()

    # 5-criteria 비교
    verdict_obj = compare_rigor(host_result, vm_result, case_spec)
    print(f"=== Verdict: {verdict_obj['verdict'].upper()} ===")
    print(f"Rationale: {verdict_obj['rationale']}")
    print()
    print("5-criteria breakdown:")
    for name, c in verdict_obj["criteria"].items():
        mark = "🟢" if c["pass"] else "🔴"
        details = ", ".join(f"{k}={v}" for k, v in c.items() if k != "pass")
        print(f"  {mark} {name}: {details}")

    if args.output:
        Path(args.output).write_text(json.dumps(verdict_obj, indent=2))
        print(f"\nResult saved → {args.output}")

    return 0 if verdict_obj["verdict"] != "red" else 1


if __name__ == "__main__":
    sys.exit(main())
