"""
v0.9 Dual-GPU Dispatch Module — host RTX 5070 Ti + VM Tesla V100 over SSH.

단일-파일 설계. 다음을 포함:
- Stage → GPU 라우팅 (``GPULocation`` enum + ``route_stage()``)
- OpenMM platform 선택 + DeviceIndex 명시 (``select_openmm_platform()``)
- SSH 기반 원격 실행기 + rsync 동기 (``VMExecutor``)
- Fallback chain (VM SSH 실패 → host, V100 OOM → 5070 Ti retry: ``with_fallback()``)

기본 동작: opt-in (``UPDD_VM_ENABLE`` env var).
- ``UPDD_VM_ENABLE=0`` (default) → 100% host 실행 (ChemRxiv v2/v3 재생산성 보장)
- ``UPDD_VM_ENABLE=1``           → Stages 9-12 가 VM V100 로 dispatch
  (Phase 2 numerical rigor validation 통과 후 default-on 검토,
   [[ADR-0007_dual_gpu_dispatch_architecture]] + rigor-priority 원칙)

Rigor-priority 원칙 (사용자 명시 2026-05-28):
V100 가 5070 Ti 대비 더 rigorous 한 것으로 증명되면 publication numbers 변동 수용.
SciVal C7 (publication-direct-use 금지) override.

환경변수:
- ``UPDD_CUDA_DEVICE`` (default "0"): cupy/gpu4pyscf device index
- ``UPDD_MMGBSA_CUDA_DEVICE`` (default "0"): OpenMM MM-GBSA DeviceIndex
- ``UPDD_MD_CUDA_DEVICE`` (default "0"): OpenMM MD DeviceIndex
- ``UPDD_VM_ENABLE`` (default "0"): "1" = dual-GPU dispatch 활성
- ``UPDD_VM_SSH_TARGET`` (default "san@192.168.122.155")
- ``UPDD_VM_PROJECT_ROOT`` (default "/home/san/UPDD_proj")
- ``UPDD_VM_SCRATCH_ROOT`` (default "/var/scratch")
- ``UPDD_VM_SYNC_METHOD`` (default "rsync")
"""
from __future__ import annotations

import enum
import os
import shlex
import subprocess
import time
from typing import Callable, Dict, List, Optional, Tuple


# ──────────────────────────────────────────────────────────────
# 환경 변수
# ──────────────────────────────────────────────────────────────
UPDD_CUDA_DEVICE = int(os.environ.get("UPDD_CUDA_DEVICE", "0"))
UPDD_MMGBSA_CUDA_DEVICE = int(os.environ.get("UPDD_MMGBSA_CUDA_DEVICE", "0"))
UPDD_MD_CUDA_DEVICE = int(os.environ.get("UPDD_MD_CUDA_DEVICE", "0"))

UPDD_VM_ENABLE = os.environ.get("UPDD_VM_ENABLE", "0") == "1"
UPDD_VM_SSH_TARGET = os.environ.get("UPDD_VM_SSH_TARGET", "san@192.168.122.155")
UPDD_VM_PROJECT_ROOT = os.environ.get("UPDD_VM_PROJECT_ROOT", "/home/san/UPDD_proj")
UPDD_VM_SCRATCH_ROOT = os.environ.get("UPDD_VM_SCRATCH_ROOT", "/home/san/scratch")
UPDD_VM_SYNC_METHOD = os.environ.get("UPDD_VM_SYNC_METHOD", "rsync")


# ──────────────────────────────────────────────────────────────
# GPU 라우팅
# ──────────────────────────────────────────────────────────────
class GPULocation(enum.Enum):
    HOST_5070TI = "host_5070ti"
    VM_V100 = "vm_v100"
    CPU_ONLY = "cpu_only"


# UPDD.py 의 13-stage 실제 stage 이름 + alias 모두 등록 (route_stage 의 robustness)
_STAGE_ROUTING: Dict[str, GPULocation] = {
    # AI inference / generation — host 5070 Ti (Blackwell FP16)
    "rfdiffusion":      GPULocation.HOST_5070TI,
    "rfdiff":           GPULocation.HOST_5070TI,
    "proteinmpnn":      GPULocation.HOST_5070TI,
    "mpnn":             GPULocation.HOST_5070TI,
    "alphafold2":       GPULocation.HOST_5070TI,
    "af2":              GPULocation.HOST_5070TI,
    "colabfold":        GPULocation.HOST_5070TI,
    # CPU 전용 stages
    "input":            GPULocation.CPU_ONLY,
    "grafting":         GPULocation.CPU_ONLY,
    "preprocess":       GPULocation.CPU_ONLY,
    "ncaa_mutation":    GPULocation.CPU_ONLY,
    "ncaa_mutate":      GPULocation.CPU_ONLY,
    "parameterize":     GPULocation.CPU_ONLY,
    "admet":            GPULocation.CPU_ONLY,
    "ranking":          GPULocation.CPU_ONLY,
    "pre_filter":       GPULocation.CPU_ONLY,
    "cofactor_reinsert":GPULocation.CPU_ONLY,
    # Chemistry compute — VM V100 (FP64 mature, Volta QM literature)
    "restrained_md":    GPULocation.VM_V100,
    "md":               GPULocation.VM_V100,
    "snapshot_extract": GPULocation.VM_V100,
    "snapshots":        GPULocation.VM_V100,
    "qmmm":             GPULocation.VM_V100,
    "qm_mm":            GPULocation.VM_V100,
    "mmgbsa":           GPULocation.VM_V100,
    "mm_gbsa":          GPULocation.VM_V100,
}


def route_stage(stage_name: str) -> GPULocation:
    """Stage 이름 → 실행 위치 결정.

    - ``UPDD_VM_ENABLE=0`` (default) 일 때 ``VM_V100`` 후보는 ``HOST_5070TI`` 로 강제.
      ChemRxiv v2/v3 reproducibility 보장 (opt-in 원칙).
    - 미등록 stage 는 ``HOST_5070TI`` default (보수적 — VM dispatch 는 명시적 routing 만).
    """
    base = _STAGE_ROUTING.get(stage_name.strip().lower(), GPULocation.HOST_5070TI)
    if not UPDD_VM_ENABLE and base == GPULocation.VM_V100:
        return GPULocation.HOST_5070TI
    return base


# ──────────────────────────────────────────────────────────────
# OpenMM platform 선택 helper
# ──────────────────────────────────────────────────────────────
def select_openmm_platform(
    preferred: str = "CUDA",
    device_id: Optional[int] = None,
    precision: str = "mixed",
) -> Tuple["openmm.Platform", Dict[str, str]]:  # type: ignore[name-defined]
    """OpenMM platform 선택 + DeviceIndex 명시 ``platformProperties`` 생성.

    Fallback chain: ``preferred`` → CUDA → OpenCL → CPU.

    Args:
        preferred: 선호 platform 이름 ("CUDA" / "OpenCL" / "CPU").
        device_id: CUDA device index. None 이면 ``UPDD_CUDA_DEVICE`` env var 사용.
        precision: CUDA precision ("mixed" default, "single" / "double" 도 가능).

    Returns:
        (platform, properties) — ``properties`` 는 ``Simulation()`` 의
        ``platformProperties`` 인자로 그대로 전달. Non-CUDA platform 은 빈 dict.
    """
    import openmm as mm  # lazy import (openmm 미설치 환경에서도 module load 가능)

    platform = None
    chosen_name: Optional[str] = None
    for name in (preferred, "CUDA", "OpenCL", "CPU"):
        if name is None:
            continue
        try:
            platform = mm.Platform.getPlatformByName(name)
            chosen_name = name
            break
        except Exception:
            continue
    if platform is None:
        platform = mm.Platform.getPlatformByName("CPU")
        chosen_name = "CPU"

    properties: Dict[str, str] = {}
    if chosen_name == "CUDA":
        dev = device_id if device_id is not None else UPDD_CUDA_DEVICE
        properties["DeviceIndex"] = str(dev)
        properties["Precision"] = precision
    return platform, properties


# ──────────────────────────────────────────────────────────────
# SSH 기반 VM Executor
# ──────────────────────────────────────────────────────────────
class VMExecutor:
    """V100 VM 원격 실행기. SSH 명령 + rsync 동기 + retry chain.

    호출 패턴::

        vm = VMExecutor()
        if not vm.is_connected():
            raise RuntimeError("VM 연결 실패")
        vm.sync_to_vm(local_case_dir)
        result = vm.execute(
            "conda run -n qmmm python utils/run_qmmm.py ...",
            cwd=os.path.join(vm.project_root, "outputs", basename),
            timeout=3600,
        )
        if result["returncode"] == 0:
            vm.sync_from_vm(remote_case_dir, local_case_dir,
                            patterns=["snapshots/", "qmmm_results/", "mmgbsa_results/"])
    """

    def __init__(
        self,
        ssh_target: str = UPDD_VM_SSH_TARGET,
        project_root: str = UPDD_VM_PROJECT_ROOT,
        sync_method: str = UPDD_VM_SYNC_METHOD,
    ) -> None:
        self.ssh_target = ssh_target
        self.project_root = project_root
        self.sync_method = sync_method

    # ── connectivity ───────────────────────────────────────────
    def is_connected(self, timeout: int = 5) -> bool:
        """SSH 연결 가능 여부 read-only check."""
        try:
            r = subprocess.run(
                [
                    "ssh",
                    "-o", f"ConnectTimeout={timeout}",
                    "-o", "StrictHostKeyChecking=no",
                    "-o", "LogLevel=ERROR",
                    self.ssh_target, "echo ok",
                ],
                capture_output=True, text=True, timeout=timeout + 5,
            )
            return r.returncode == 0 and "ok" in r.stdout
        except Exception:
            return False

    # ── execute ────────────────────────────────────────────────
    def execute(
        self,
        cmd: str,
        cwd: Optional[str] = None,
        timeout: int = 3600,
        env: Optional[Dict[str, str]] = None,
        retry: int = 3,
    ) -> Dict[str, object]:
        """SSH 명령 실행. Retry chain (5s/30s/120s) — 일시적 네트워크 단절 대응.

        Returns:
            ``{returncode, stdout, stderr, elapsed_s, retried}``.
            ``returncode == 0`` 만 성공으로 간주.
        """
        delays: List[int] = [5, 30, 120][:retry]
        env_prefix = " ".join(f"{k}={shlex.quote(v)}" for k, v in (env or {}).items())
        cd_prefix = f"cd {shlex.quote(cwd)} && " if cwd else ""
        full_cmd = f"{cd_prefix}{env_prefix} {cmd}".strip()

        start = time.time()
        retried = 0
        last_stderr = ""
        for delay in delays:
            try:
                r = subprocess.run(
                    [
                        "ssh",
                        "-o", "StrictHostKeyChecking=no",
                        "-o", "LogLevel=ERROR",
                        self.ssh_target, full_cmd,
                    ],
                    capture_output=True, text=True, timeout=timeout,
                )
                if r.returncode == 0:
                    return {
                        "returncode": 0,
                        "stdout": r.stdout,
                        "stderr": r.stderr,
                        "elapsed_s": time.time() - start,
                        "retried": retried,
                    }
                last_stderr = (r.stderr or "")[:500]
            except subprocess.TimeoutExpired:
                last_stderr = f"SSH timeout ({timeout}s)"
            except Exception as e:
                last_stderr = f"{type(e).__name__}: {str(e)[:300]}"
            print(
                f"  [VM] 실행 실패 (retry={retried + 1}/{len(delays)}) — "
                f"{delay}s 후 재시도. last_stderr={last_stderr[:120]}"
            )
            time.sleep(delay)
            retried += 1
        return {
            "returncode": -1,
            "stdout": "",
            "stderr": last_stderr or "Max retries exceeded",
            "elapsed_s": time.time() - start,
            "retried": retried,
        }

    # ── rsync sync ─────────────────────────────────────────────
    @staticmethod
    def _default_excludes() -> List[str]:
        # 대용량 파일 (DCD/NC/CHK) 은 boundary 별 selective sync 으로만 이동
        return [
            "*.dcd", "*.nc", "*.h5",
            "__pycache__/", ".git/", "_archive/", "_audit/",
            "outputs/", "archive/", "log/", "chkfile_archive/",
            "paper1*", "inputs/", "analysis/",
        ]

    def sync_to_vm(
        self,
        local_path: str,
        remote_path: Optional[str] = None,
        excludes: Optional[List[str]] = None,
        timeout: int = 600,
    ) -> bool:
        """rsync local → VM. local_path 끝 ``/`` 자동 부여 (디렉토리 내용 sync)."""
        if not local_path.endswith("/"):
            local_path += "/"
        if remote_path is None:
            remote_path = os.path.join(
                self.project_root, "outputs", os.path.basename(local_path.rstrip("/"))
            )
        if not remote_path.endswith("/"):
            remote_path += "/"
        # remote dir 자동 생성
        mk = self.execute(f"mkdir -p {shlex.quote(remote_path)}", retry=2, timeout=30)
        if mk["returncode"] != 0:
            print(f"  [VM rsync→] remote mkdir 실패: {mk['stderr'][:100]}")
            return False
        ex = excludes if excludes is not None else self._default_excludes()
        exclude_args: List[str] = []
        for e in ex:
            exclude_args.extend(["--exclude", e])
        cmd = [
            "rsync", "-az", "--delete-after",
            "-e", "ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR",
            *exclude_args,
            local_path,
            f"{self.ssh_target}:{remote_path}",
        ]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            return r.returncode == 0
        except Exception as e:
            print(f"  [VM rsync→] 실패: {type(e).__name__}: {str(e)[:100]}")
            return False

    def sync_from_vm(
        self,
        remote_path: str,
        local_path: str,
        patterns: Optional[List[str]] = None,
        timeout: int = 600,
    ) -> bool:
        """rsync VM → local. ``patterns`` 로 selective copy (e.g. ``["snapshots/", "qmmm_results/"]``).

        Patterns 미지정 시 전체 디렉토리 sync (대용량 파일 위험 — 권장 명시).
        """
        if not remote_path.endswith("/"):
            remote_path += "/"
        if not local_path.endswith("/"):
            local_path += "/"
        os.makedirs(local_path, exist_ok=True)

        cmd = [
            "rsync", "-az", "--prune-empty-dirs",
            "-e", "ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR",
        ]
        if patterns:
            for p in patterns:
                # include 디렉토리 + 그 안 모든 파일
                p_dir = p.rstrip("/") + "/"
                cmd.extend(["--include", p_dir, "--include", f"{p_dir}**"])
            cmd.extend(["--include", "*/", "--exclude", "*"])
        cmd.extend([
            f"{self.ssh_target}:{remote_path}",
            local_path,
        ])
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            return r.returncode == 0
        except Exception as e:
            print(f"  [VM rsync←] 실패: {type(e).__name__}: {str(e)[:100]}")
            return False


# ──────────────────────────────────────────────────────────────
# Fallback chain
# ──────────────────────────────────────────────────────────────
def with_fallback(
    primary_loc: GPULocation,
    primary_fn: Callable,
    host_fallback_fn: Callable,
    *args,
    **kwargs,
):
    """Primary location 시도 → 실패 시 host fallback.

    - ``primary_loc == VM_V100`` 일 때만 fallback chain 적용. 그 외는 그냥 primary 실행.
    - ``primary_fn`` 이 dict 를 반환하고 ``returncode != 0`` 이거나 예외 발생 시 fallback.
    - Fallback 도 실패하면 마지막 예외 그대로 raise.

    예시::

        result = with_fallback(
            GPULocation.VM_V100,
            lambda: vm.execute("python run_qmmm.py ...", timeout=3600),
            lambda: subprocess.run(["python", "run_qmmm.py", ...], check=True),
        )
    """
    if primary_loc != GPULocation.VM_V100:
        return primary_fn(*args, **kwargs)
    try:
        result = primary_fn(*args, **kwargs)
        if isinstance(result, dict) and result.get("returncode", 0) != 0:
            raise RuntimeError(
                f"VM execution returned rc={result.get('returncode')}: "
                f"{str(result.get('stderr', ''))[:300]}"
            )
        return result
    except Exception as e:
        print(
            f"  [VM→HOST fallback] VM 실패 — host 5070 Ti 로 재시도. "
            f"원인: {type(e).__name__}: {str(e)[:150]}"
        )
        return host_fallback_fn(*args, **kwargs)


# ──────────────────────────────────────────────────────────────
# 모듈 self-test (python -m utils.dispatch)
# ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys
    print(f"UPDD_VM_ENABLE        : {UPDD_VM_ENABLE}")
    print(f"UPDD_VM_SSH_TARGET    : {UPDD_VM_SSH_TARGET}")
    print(f"UPDD_VM_PROJECT_ROOT  : {UPDD_VM_PROJECT_ROOT}")
    print(f"UPDD_CUDA_DEVICE      : {UPDD_CUDA_DEVICE}")
    print()
    print("Stage routing:")
    for stage in ("rfdiffusion", "af2", "qmmm", "mmgbsa", "ranking"):
        print(f"  {stage:15s} → {route_stage(stage).value}")
    print()
    if UPDD_VM_ENABLE:
        vm = VMExecutor()
        ok = vm.is_connected()
        print(f"VM SSH connectivity   : {'OK' if ok else 'FAIL'}")
        sys.exit(0 if ok else 1)
