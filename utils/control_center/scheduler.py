"""Exact-unit systemd user scheduler; no shell and no broad process matching."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from . import settings
from .control_token import initialize
from .models import JobSpec


class SchedulerError(RuntimeError):
    pass


class SystemdScheduler:
    def __init__(self, state_dir: Path | None = None):
        self.state_dir = (state_dir or settings.state_root()).resolve()
        self.state_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def unit_name(job_id: str) -> str:
        safe = "".join(ch for ch in job_id.lower() if ch.isalnum())
        if not safe:
            raise ValueError("job id has no unit-safe characters")
        return f"updd-job-{safe[:48]}"

    def paths(self, job_id: str) -> dict[str, Path]:
        root = self.state_dir / "jobs" / job_id
        return {
            "root": root,
            "token": root / "control.json",
            "record": root / "wrapper.json",
            "log": root / "job.log",
        }

    def launch(self, job_id: str, spec: JobSpec) -> dict[str, str]:
        paths = self.paths(job_id)
        paths["root"].mkdir(parents=True, exist_ok=True)
        token = initialize(paths["token"], job_id, spec.spec_hash, spec.input_digest)
        unit = f"{self.unit_name(job_id)}-g{int(token['generation'])}"
        wrapper = (settings.REPO_ROOT / "utils/control_center/job_wrapper.py").resolve()
        command = [
            sys.executable,
            str(wrapper),
            "--record",
            str(paths["record"]),
            "--token",
            str(paths["token"]),
            "--job-id",
            job_id,
            "--spec-hash",
            spec.spec_hash,
            "--input-digest",
            spec.input_digest,
            "--",
            *spec.command,
        ]
        for key, value in spec.environment.items():
            if value is not None:
                command = ["/usr/bin/env", f"{key}={value}", *command]
        args = [
            "systemd-run",
            "--user",
            f"--unit={unit}",
            f"--description=UPDD control job {job_id}",
            "--property=Type=exec",
            "--property=KillMode=control-group",
            "--property=TimeoutStopSec=30s",
            "--property=SuccessExitStatus=75",
            f"--property=StandardOutput=append:{paths['log']}",
            f"--property=StandardError=append:{paths['log']}",
            f"--working-directory={spec.cwd}",
            *command,
        ]
        result = subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            timeout=20,
            check=False,
        )
        if result.returncode != 0:
            raise SchedulerError(result.stdout.strip() or f"systemd-run exit {result.returncode}")
        return {
            "unit_name": unit,
            "control_token": str(paths["token"]),
            "wrapper_record": str(paths["record"]),
            "primary_log": str(paths["log"]),
        }

    @staticmethod
    def show(unit: str) -> dict[str, str]:
        result = subprocess.run(
            [
                "systemctl",
                "--user",
                "show",
                unit,
                "--property=LoadState,ActiveState,SubState,Result,ExecMainStatus,MainPID",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            encoding="utf-8",
            timeout=10,
            check=False,
        )
        values: dict[str, str] = {}
        for line in result.stdout.splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                values[key] = value
        if result.returncode != 0 and not values:
            values["LoadState"] = "not-found"
        return values

    @classmethod
    def active(cls, unit: str) -> bool:
        return cls.show(unit).get("ActiveState") in {"activating", "active", "reloading"}

    @staticmethod
    def stop(unit: str, force: bool = False) -> None:
        if not unit.startswith("updd-job-"):
            raise SchedulerError(f"refusing non-UPDD unit: {unit}")
        if force:
            command = [
                "systemctl", "--user", "kill", "--signal=SIGKILL", "--kill-whom=all", unit
            ]
        else:
            command = ["systemctl", "--user", "stop", unit]
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            timeout=45,
            check=False,
        )
        if result.returncode != 0 and "not loaded" not in result.stdout.lower():
            raise SchedulerError(result.stdout.strip() or f"stop exit {result.returncode}")

    @staticmethod
    def wrapper_record(path: str | None) -> dict[str, Any]:
        if not path:
            return {}
        try:
            payload = json.loads(Path(path).read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError):
            return {}
        return payload if isinstance(payload, dict) else {}
