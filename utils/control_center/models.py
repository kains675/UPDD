"""Canonical immutable job model and orthogonal status axes."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Mapping

from . import settings


class ExecutionStatus(str, Enum):
    DRAFT = "DRAFT"
    PREFLIGHT = "PREFLIGHT"
    QUEUED = "QUEUED"
    RUNNING = "RUNNING"
    PAUSE_REQUESTED = "PAUSE_REQUESTED"
    PAUSED = "PAUSED"
    COMPLETE = "COMPLETE"
    FAILED = "FAILED"
    STOPPED = "STOPPED"
    ORPHANED = "ORPHANED"


class ArtifactStatus(str, Enum):
    PRESENT = "PRESENT"
    PARTIAL = "PARTIAL"
    MISSING = "MISSING"


class ScientificStatus(str, Enum):
    VALID = "VALID"
    PARTIAL = "PARTIAL"
    INVALID = "INVALID"
    UNKNOWN = "UNKNOWN"


TERMINAL_EXECUTION_STATUSES = {
    ExecutionStatus.COMPLETE,
    ExecutionStatus.FAILED,
    ExecutionStatus.STOPPED,
    ExecutionStatus.ORPHANED,
}


def _clean_json(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("job spec cannot contain non-finite numbers")
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _clean_json(item) for key, item in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [_clean_json(item) for item in value]
    raise TypeError(f"unsupported job-spec value: {type(value).__name__}")


def canonical_json(value: Any) -> str:
    return json.dumps(
        _clean_json(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True
    )


@dataclass(frozen=True)
class JobSpec:
    name: str
    campaign_id: str
    adapter: str
    interpreter: str
    script: str
    argv: tuple[str, ...] = ()
    cwd: str = str(settings.REPO_ROOT)
    output_root: str = ""
    resources: Mapping[str, Any] = field(default_factory=dict)
    dependencies: tuple[str, ...] = ()
    declared_parameters: Mapping[str, Any] = field(default_factory=dict)
    default_parameters: Mapping[str, Any] = field(default_factory=dict)
    effective_parameters: Mapping[str, Any] = field(default_factory=dict)
    parameter_sources: Mapping[str, Any] = field(default_factory=dict)
    input_digests: Mapping[str, str] = field(default_factory=dict)
    environment: Mapping[str, str | None] = field(default_factory=dict)
    log_paths: tuple[str, ...] = ()
    readonly: bool = False
    source_key: str | None = None
    schema: str = "updd_control_job_v1"

    def __post_init__(self) -> None:
        if self.schema != "updd_control_job_v1":
            raise ValueError(f"unsupported job schema: {self.schema}")
        if not self.name.strip() or not self.campaign_id.strip():
            raise ValueError("name and campaign_id are required")
        if self.adapter not in settings.ALLOWED_SCRIPTS:
            raise ValueError(f"adapter is not allow-listed: {self.adapter}")
        interpreter = Path(self.interpreter).expanduser()
        script = Path(self.script).expanduser()
        cwd = Path(self.cwd).expanduser()
        output = Path(self.output_root).expanduser() if self.output_root else None
        for label, path in (("interpreter", interpreter), ("script", script), ("cwd", cwd)):
            if not path.is_absolute():
                raise ValueError(f"{label} must be absolute: {path}")
        if output is not None and not output.is_absolute():
            raise ValueError(f"output_root must be absolute: {output}")
        if interpreter.resolve() not in settings.ALLOWED_INTERPRETERS:
            raise ValueError(f"interpreter is not allow-listed: {interpreter}")
        if script.resolve() not in settings.ALLOWED_SCRIPTS[self.adapter]:
            raise ValueError(f"script is not allow-listed for {self.adapter}: {script}")
        if any("\x00" in arg for arg in self.argv):
            raise ValueError("argv cannot contain NUL")
        unknown_env = set(self.environment) - settings.ALLOWED_ENVIRONMENT_KEYS
        if unknown_env:
            raise ValueError(f"environment keys are not allow-listed: {sorted(unknown_env)}")
        gpu = self.resources.get("gpu")
        if gpu is not None and gpu not in {"host_cuda0", "v100_vm", "selftest_gpu"}:
            raise ValueError(f"unknown GPU resource: {gpu}")
        _clean_json(self.to_dict())

    @property
    def command(self) -> list[str]:
        return [self.interpreter, self.script, *self.argv]

    @property
    def spec_hash(self) -> str:
        return hashlib.sha256(canonical_json(self.to_dict()).encode("ascii")).hexdigest()

    @property
    def input_digest(self) -> str:
        return hashlib.sha256(canonical_json(self.input_digests).encode("ascii")).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "name": self.name,
            "campaign_id": self.campaign_id,
            "adapter": self.adapter,
            "interpreter": str(Path(self.interpreter).resolve()),
            "script": str(Path(self.script).resolve()),
            "argv": list(self.argv),
            "cwd": str(Path(self.cwd).resolve()),
            "output_root": str(Path(self.output_root).resolve()) if self.output_root else "",
            "resources": _clean_json(self.resources),
            "dependencies": list(self.dependencies),
            "declared_parameters": _clean_json(self.declared_parameters),
            "default_parameters": _clean_json(self.default_parameters),
            "effective_parameters": _clean_json(self.effective_parameters),
            "parameter_sources": _clean_json(self.parameter_sources),
            "input_digests": _clean_json(self.input_digests),
            "environment": _clean_json(self.environment),
            "log_paths": list(self.log_paths),
            "readonly": self.readonly,
            "source_key": self.source_key,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "JobSpec":
        return cls(
            schema=str(payload.get("schema", "updd_control_job_v1")),
            name=str(payload["name"]),
            campaign_id=str(payload["campaign_id"]),
            adapter=str(payload["adapter"]),
            interpreter=str(payload["interpreter"]),
            script=str(payload["script"]),
            argv=tuple(str(value) for value in payload.get("argv", [])),
            cwd=str(payload.get("cwd", settings.REPO_ROOT)),
            output_root=str(payload.get("output_root", "")),
            resources=dict(payload.get("resources", {})),
            dependencies=tuple(str(value) for value in payload.get("dependencies", [])),
            declared_parameters=dict(payload.get("declared_parameters", {})),
            default_parameters=dict(payload.get("default_parameters", {})),
            effective_parameters=dict(payload.get("effective_parameters", {})),
            parameter_sources=dict(payload.get("parameter_sources", {})),
            input_digests={str(k): str(v) for k, v in payload.get("input_digests", {}).items()},
            environment={str(k): (None if v is None else str(v)) for k, v in payload.get("environment", {}).items()},
            log_paths=tuple(str(value) for value in payload.get("log_paths", [])),
            readonly=bool(payload.get("readonly", False)),
            source_key=None if payload.get("source_key") is None else str(payload["source_key"]),
        )
