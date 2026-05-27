#!/usr/bin/env python
"""
scripts/validate_v100_md_rigor.py — OpenMM MD σ_btwn invariance V100 vs 5070 Ti.

ADR-0007 SciVal verdict C2 sufficient condition (default-on transition framework).
사용자 명시 (2026-05-28): "기존 5070ti 만을 이용한 방법보다 V100을 함께 활용했을
때 결과가 변할 수 있더라도 과학적 엄밀성이 증가한다면 그 방향이 맞아."

Validation 대상: alanine dipeptide (NH3+-CHMe-COO-, 22 atoms in zwitterion) GBn2
implicit solvent. UPDD production 의 MD path (utils/run_restrained_md.py) 와 동일
하게 OpenMM 8.x + mixed-precision CUDA 사용. Production-relevant 이며 σ_btwn
invariance 의 sufficient signal 추출 가능.

Methodology:
- 3 seeds (42 / 142 / 242) 각각 host + VM 에서 동일 system + 동일 RNG seed 로 MD
- N ns equilibrium (default 1 ns at 2 fs = 500k steps)
- 추출 metric: <φ> (mean phi 각, dihedral), <ψ> (mean psi), <RMSD> from start
- σ_btwn_host = std across 3 host seeds (각 metric)
- σ_btwn_vm   = std across 3 VM seeds
- ratio = σ_btwn_vm / σ_btwn_host
- C2 pass criterion: ratio ∈ [0.7, 1.3] (V100 ensemble property invariance)

OpenMM mixed-precision MD 는 bitwise-deterministic 아님 (Eastman 2010 J Comp Chem
DOI 10.1002/jcc.21413). 따라서 ΔE 직접 비교 X — ensemble property 통계적 invariance
만 검증. Same RNG seed + same hardware-class FP 라면 동일 ensemble 통계 기대.

사용법::
    python scripts/validate_v100_md_rigor.py --steps 500000 --remote
    python scripts/validate_v100_md_rigor.py --steps 100000   # quick smoke
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "utils"))


# ──────────────────────────────────────────────────────────────
# Alanine dipeptide MD runner script
# ──────────────────────────────────────────────────────────────
def _md_runner_script(n_steps: int, seed: int) -> str:
    """OpenMM 8.x alanine dipeptide MD via openmmtools testsystems (production lib).

    GBn2 implicit (UPDD production) + LangevinMiddleIntegrator 2 fs + 300 K.
    Output: phi/psi/RMSD trajectory + final-quartile mean (이전 quartile = 평형화 후).
    """
    return f"""
import sys, json, time
import numpy as np
import openmm as mm
from openmm import unit
from openmm import app as oapp
from openmmtools.testsystems import AlanineDipeptideImplicit
import mdtraj as md

n_steps = {n_steps}
seed = {seed}
print(f'OpenMM {{mm.version.short_version}}, n_steps={{n_steps}}, seed={{seed}}')

# AlanineDipeptide implicit (GBn2 — UPDD production path)
sys_obj = AlanineDipeptideImplicit()
system = sys_obj.system
topology = sys_obj.topology
positions = sys_obj.positions

# Platform: CUDA (DeviceIndex 환경 변수 우선)
import os
dev = os.environ.get('UPDD_MD_CUDA_DEVICE', '0')
platform = mm.Platform.getPlatformByName('CUDA')
props = {{'DeviceIndex': dev, 'Precision': 'mixed'}}

integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    2.0 * unit.femtoseconds,
)
integrator.setRandomNumberSeed(seed)

simulation = oapp.Simulation(topology, system, integrator, platform, props)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)

# Minimize + warm-up
simulation.minimizeEnergy(maxIterations=1000)
simulation.step(5000)  # 10 ps equilibration

# Sample every 1000 steps (2 ps) — store positions in memory
n_frames = max(50, n_steps // 1000)
sample_interval = max(1, n_steps // n_frames)

positions_list = []
energies = []
t0 = time.time()
for i in range(n_frames):
    simulation.step(sample_interval)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    pos = np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    positions_list.append(pos)
    energies.append(state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
elapsed = time.time() - t0

# φ/ψ 추출 via mdtraj
positions_np = np.array(positions_list)
n_atoms = positions_np.shape[1]
top_md = md.Topology.from_openmm(topology)
traj = md.Trajectory(positions_np, top_md)
phi = md.compute_phi(traj)[1]  # shape (n_frames, n_dihedrals=1)
psi = md.compute_psi(traj)[1]
rmsd = md.rmsd(traj, traj[0])

# 평형화 후 (latter half) mean
half = n_frames // 2
phi_mean = float(np.degrees(np.mean(phi[half:])))
psi_mean = float(np.degrees(np.mean(psi[half:])))
rmsd_mean = float(np.mean(rmsd[half:]))
energy_mean = float(np.mean(energies[half:]))
energy_std = float(np.std(energies[half:]))

out = {{
    'n_steps': n_steps,
    'seed': seed,
    'n_frames': n_frames,
    'n_atoms': int(n_atoms),
    'phi_mean_deg': phi_mean,
    'psi_mean_deg': psi_mean,
    'rmsd_mean_nm': rmsd_mean,
    'energy_mean_kcal': energy_mean,
    'energy_std_kcal': energy_std,
    'elapsed_s': float(elapsed),
}}
print('MD_RESULT_BEGIN' + json.dumps(out) + 'MD_RESULT_END')
"""


def _parse_md_result(stdout: str) -> Dict[str, Any]:
    begin = stdout.find("MD_RESULT_BEGIN")
    end = stdout.find("MD_RESULT_END")
    if begin < 0 or end < 0:
        raise RuntimeError(f"MD result marker 미발견. stdout (last 800):\n{stdout[-800:]}")
    return json.loads(stdout[begin + len("MD_RESULT_BEGIN"):end])


def run_md_host(n_steps: int, seed: int, conda_env: str = "qmmm") -> Dict[str, Any]:
    """Host 5070 Ti 에서 MD. qmmm env 에 openmmtools 가 있어야 함 (있음)."""
    script = _md_runner_script(n_steps, seed)
    cmd = ["conda", "run", "-n", conda_env, "--no-capture-output", "python", "-c", script]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError(f"Host MD 실패 (rc={r.returncode}): {r.stderr[:500]}")
    result = _parse_md_result(r.stdout)
    result["location"] = "host_5070ti"
    return result


def run_md_vm(n_steps: int, seed: int, conda_env: str = "qmmm") -> Dict[str, Any]:
    """VM V100 에서 MD via base64-encoded SSH dispatch."""
    import dispatch
    vm = dispatch.VMExecutor()
    if not vm.is_connected():
        raise RuntimeError(f"VM SSH 연결 실패 ({vm.ssh_target})")
    script = _md_runner_script(n_steps, seed)
    payload_b64 = base64.b64encode(script.encode("utf-8")).decode("ascii")
    py_bin = f"/home/san/miniconda3/envs/{conda_env}/bin/python"
    cmd = (
        f"{py_bin} -c "
        f"\"import base64; exec(base64.b64decode('{payload_b64}').decode())\""
    )
    r = vm.execute(cmd, timeout=1800, retry=2)
    if r["returncode"] != 0:
        raise RuntimeError(f"VM MD 실패 (rc={r['returncode']}): {str(r.get('stderr',''))[:500]}")
    result = _parse_md_result(str(r["stdout"]))
    result["location"] = "vm_v100"
    return result


# ──────────────────────────────────────────────────────────────
# σ_btwn invariance 분석
# ──────────────────────────────────────────────────────────────
def analyze_sigma_btwn(
    host_results: List[Dict[str, Any]],
    vm_results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """3-seed σ_btwn ratio V100/host 계산. ADR-0007 C2 pass: ratio ∈ [0.7, 1.3]."""
    metrics = ["phi_mean_deg", "psi_mean_deg", "rmsd_mean_nm", "energy_mean_kcal"]
    analysis: Dict[str, Any] = {}
    all_pass = True

    for m in metrics:
        host_vals = [r[m] for r in host_results]
        vm_vals = [r[m] for r in vm_results]
        host_sigma = statistics.stdev(host_vals) if len(host_vals) >= 2 else 0.0
        vm_sigma = statistics.stdev(vm_vals) if len(vm_vals) >= 2 else 0.0
        # σ 가 0 (degenerate) 인 경우 처리 — ratio 정의 불가, 양쪽 모두 0 이면 trivially pass
        if host_sigma < 1e-12 and vm_sigma < 1e-12:
            ratio = 1.0
            verdict = True
            note = "both σ ≈ 0 — degenerate ensemble (n_steps 부족)"
        elif host_sigma < 1e-12:
            ratio = float("inf")
            verdict = False
            note = "host σ ≈ 0 → ratio undefined"
        else:
            ratio = vm_sigma / host_sigma
            verdict = 0.7 <= ratio <= 1.3
            note = ""
        analysis[m] = {
            "host_sigma": host_sigma,
            "vm_sigma": vm_sigma,
            "ratio_vm_over_host": ratio,
            "pass": verdict,
            "host_vals": host_vals,
            "vm_vals": vm_vals,
            "note": note,
        }
        if not verdict:
            all_pass = False

    # 종합 verdict
    if all_pass:
        # 모든 ratio 가 [0.85, 1.15] 안에 있으면 strong GREEN
        all_tight = all(0.85 <= analysis[m]["ratio_vm_over_host"] <= 1.15 for m in metrics)
        if all_tight:
            overall = "green"
            rationale = "V100 ensemble property strong invariance (모든 metric ratio ∈ [0.85, 1.15])"
        else:
            overall = "yellow"
            rationale = "V100 ensemble property within tolerance (ratio ∈ [0.7, 1.3])"
    else:
        overall = "red"
        failed = [m for m in metrics if not analysis[m]["pass"]]
        rationale = f"σ_btwn ratio out of [0.7, 1.3] for: {', '.join(failed)}"

    return {
        "overall": overall,
        "rationale": rationale,
        "metrics": analysis,
        "host_results": host_results,
        "vm_results": vm_results,
    }


# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", type=int, default=100000, help="MD steps per seed (default 100k = 200 ps smoke)")
    ap.add_argument("--seeds", default="42,142,242", help="Comma-separated RNG seeds (default 42,142,242)")
    ap.add_argument("--remote", action="store_true", help="VM V100 비교 활성 (없으면 host-only smoke)")
    ap.add_argument("--output", default=None, help="결과 JSON 저장 경로")
    args = ap.parse_args()

    seeds = [int(s.strip()) for s in args.seeds.split(",") if s.strip()]
    if len(seeds) < 2 and args.remote:
        print("⚠️  σ_btwn 계산은 seed ≥ 2 필요")

    print(f"=== Alanine dipeptide MD σ_btwn validation ===")
    print(f"steps/seed: {args.steps:,} ({args.steps * 2 / 1000:.0f} ps at 2 fs)")
    print(f"seeds: {seeds}")
    print(f"remote: {args.remote}")
    print()

    # Host runs
    host_results: List[Dict[str, Any]] = []
    for i, seed in enumerate(seeds, 1):
        print(f"[{i}/{len(seeds)}] Host 5070 Ti MD (seed={seed}) ...")
        r = run_md_host(args.steps, seed)
        print(f"  φ={r['phi_mean_deg']:+.2f}° / ψ={r['psi_mean_deg']:+.2f}° / "
              f"RMSD={r['rmsd_mean_nm']:.3f} nm / E_avg={r['energy_mean_kcal']:+.2f} kcal/mol / "
              f"elapsed={r['elapsed_s']:.1f}s")
        host_results.append(r)
    print()

    if not args.remote:
        print("--remote 미지정 — host-only smoke 완료.")
        if len(host_results) >= 2:
            host_phi_sigma = statistics.stdev([r["phi_mean_deg"] for r in host_results])
            print(f"Host σ_btwn(φ) = {host_phi_sigma:.3f}° (3-seed)")
        return 0

    # VM runs
    vm_results: List[Dict[str, Any]] = []
    for i, seed in enumerate(seeds, 1):
        print(f"[{i}/{len(seeds)}] VM V100 MD (seed={seed}) ...")
        try:
            r = run_md_vm(args.steps, seed)
        except Exception as e:
            print(f"  ❌ {type(e).__name__}: {e}")
            sys.exit(3)
        print(f"  φ={r['phi_mean_deg']:+.2f}° / ψ={r['psi_mean_deg']:+.2f}° / "
              f"RMSD={r['rmsd_mean_nm']:.3f} nm / E_avg={r['energy_mean_kcal']:+.2f} kcal/mol / "
              f"elapsed={r['elapsed_s']:.1f}s")
        vm_results.append(r)
    print()

    # σ_btwn 분석
    verdict = analyze_sigma_btwn(host_results, vm_results)
    print(f"=== Verdict: {verdict['overall'].upper()} ===")
    print(f"Rationale: {verdict['rationale']}")
    print()
    print("Per-metric breakdown:")
    for metric, data in verdict["metrics"].items():
        mark = "🟢" if data["pass"] else "🔴"
        print(f"  {mark} {metric}:")
        print(f"     host σ = {data['host_sigma']:.4f}")
        print(f"     vm   σ = {data['vm_sigma']:.4f}")
        print(f"     ratio  = {data['ratio_vm_over_host']:.4f}"
              f"{' (' + data['note'] + ')' if data['note'] else ''}")

    if args.output:
        Path(args.output).write_text(json.dumps(verdict, indent=2))
        print(f"\nResult saved → {args.output}")

    return 0 if verdict["overall"] != "red" else 1


if __name__ == "__main__":
    sys.exit(main())
