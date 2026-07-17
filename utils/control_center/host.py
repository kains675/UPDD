"""Read-only host, GPU, VM, and repository observations."""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import psutil

from . import settings


def _run(command: list[str], timeout: float = 5.0) -> tuple[int, str]:
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            timeout=timeout,
            check=False,
        )
        return result.returncode, result.stdout.strip()
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 127, f"{type(exc).__name__}: {exc}"


def gpu() -> dict[str, Any]:
    code, output = _run(
        [
            "nvidia-smi",
            "--query-gpu=name,utilization.gpu,memory.used,memory.total,power.draw,power.limit,temperature.gpu",
            "--format=csv,noheader,nounits",
        ]
    )
    if code or not output:
        return {"available": False, "error": output}
    values = [value.strip() for value in output.splitlines()[0].split(",")]
    if len(values) != 7:
        return {"available": False, "error": output}
    return {
        "available": True,
        "name": values[0],
        "utilization_pct": float(values[1]),
        "memory_used_mb": float(values[2]),
        "memory_total_mb": float(values[3]),
        "power_w": float(values[4]),
        "power_limit_w": float(values[5]),
        "temperature_c": float(values[6]),
    }


def vm() -> dict[str, Any]:
    code, output = _run(["virsh", "domstate", "v100-vm"])
    return {"name": "v100-vm", "state": output if code == 0 else "unknown", "error": None if code == 0 else output}


def repository() -> dict[str, Any]:
    _, branch = _run(["git", "-C", str(settings.REPO_ROOT), "branch", "--show-current"])
    _, head = _run(["git", "-C", str(settings.REPO_ROOT), "rev-parse", "--short", "HEAD"])
    code, status = _run(["git", "-C", str(settings.REPO_ROOT), "status", "--porcelain"], timeout=10)
    return {"branch": branch, "head": head, "dirty": bool(status) if code == 0 else None}


def snapshot() -> dict[str, Any]:
    memory = psutil.virtual_memory()
    swap = psutil.swap_memory()
    disk = shutil.disk_usage(settings.REPO_ROOT)
    return {
        "timestamp": time.time(),
        "gpu": gpu(),
        "memory": {
            "total_gb": memory.total / 2**30,
            "available_gb": memory.available / 2**30,
            "used_pct": memory.percent,
            "swap_used_gb": swap.used / 2**30,
            "swap_total_gb": swap.total / 2**30,
        },
        "disk": {
            "total_gb": disk.total / 2**30,
            "free_gb": disk.free / 2**30,
            "used_pct": 100.0 * disk.used / disk.total,
        },
        "vm": vm(),
        "git": repository(),
    }
