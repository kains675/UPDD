"""Control-center orchestration independent from HTTP and Streamlit."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from . import host, settings
from .adapters import AdapterError, AdapterManager, legacy_processes
from .control_token import request_pause
from .models import ArtifactStatus, ExecutionStatus, JobSpec, ScientificStatus
from .registry import LeaseConflict, Registry
from .scheduler import SchedulerError, SystemdScheduler


class ServiceError(RuntimeError):
    pass


class ControlCenter:
    def __init__(
        self,
        registry: Registry | None = None,
        scheduler: SystemdScheduler | None = None,
        adapters: AdapterManager | None = None,
    ):
        self.registry = registry or Registry(settings.database_path())
        self.scheduler = scheduler or SystemdScheduler()
        self.adapters = adapters or AdapterManager()

    def discover(self) -> list[dict[str, Any]]:
        jobs = []
        for item in self.adapters.discover():
            job = self.registry.create_job(
                item.spec,
                execution_status=item.execution_status,
                artifact_status=item.artifact_status,
                scientific_status=item.scientific_status,
                ranking_eligible=item.ranking_eligible,
                progress=item.progress,
            )
            self.registry.update_runtime(job["job_id"], execution_status=item.execution_status, progress=item.progress)
            self.registry.update_evidence(
                job["job_id"],
                artifact_status=item.artifact_status,
                scientific_status=item.scientific_status,
                ranking_eligible=item.ranking_eligible,
                progress=item.progress,
            )
            jobs.append(self.registry.get_job(job["job_id"]))
        self.registry.append_event(None, "system", "DISCOVER", "OK", f"{len(jobs)} jobs observed")
        return jobs

    def register(self, payload: dict[str, Any], actor: str = "api") -> dict[str, Any]:
        spec = JobSpec.from_dict(payload)
        if not spec.readonly and spec.adapter != "self_test" and not (
            spec.source_key or ""
        ).startswith("preregistered:"):
            raise ServiceError("controllable jobs require a preregistered source_key")
        self.adapters.validate(spec, allow_readonly=True)
        job = self.registry.create_job(spec)
        self.registry.append_event(job["job_id"], actor, "REGISTER", "OK", "spec accepted")
        return job

    def preflight(self, job_id: str, actor: str = "ui") -> dict[str, Any]:
        job = self.registry.get_job(job_id)
        if job["readonly"]:
            raise ServiceError("read-only imported jobs cannot be controlled")
        if job["execution_status"] not in {
            ExecutionStatus.DRAFT.value,
            ExecutionStatus.PREFLIGHT.value,
        }:
            raise ServiceError(f"preflight requires DRAFT/PREFLIGHT, got {job['execution_status']}")
        spec = JobSpec.from_dict(job["spec"])
        if job["spec_hash"] != spec.spec_hash or job["input_digest"] != spec.input_digest:
            raise ServiceError("stored immutable spec hash mismatch")
        try:
            result = self.adapters.validate(spec)
        except Exception as exc:
            self.registry.append_event(job_id, actor, "PREFLIGHT", "FAIL", str(exc))
            raise
        self.registry.update_runtime(job_id, execution_status=ExecutionStatus.PREFLIGHT)
        self.registry.append_event(job_id, actor, "PREFLIGHT", "OK", "adapter validation passed", result)
        return result

    def _check_dependencies(self, spec: JobSpec) -> None:
        for dependency in spec.dependencies:
            if self.registry.get_job(dependency)["execution_status"] != ExecutionStatus.COMPLETE.value:
                raise ServiceError(f"dependency is not complete: {dependency}")

    def launch(self, job_id: str, actor: str = "ui", resume: bool = False) -> dict[str, Any]:
        job = self.registry.get_job(job_id)
        spec = JobSpec.from_dict(job["spec"])
        allowed = {ExecutionStatus.PREFLIGHT.value}
        if resume:
            allowed = {ExecutionStatus.PAUSED.value, ExecutionStatus.ORPHANED.value}
        if job["execution_status"] not in allowed:
            raise ServiceError(f"cannot {'resume' if resume else 'launch'} from {job['execution_status']}")
        self.adapters.validate(spec)
        if job["spec_hash"] != spec.spec_hash or job["input_digest"] != spec.input_digest:
            raise ServiceError("resume/launch identity mismatch")
        self._check_dependencies(spec)
        resource = spec.resources.get("gpu")
        try:
            if resource:
                self.registry.acquire_lease(str(resource), job_id)
            self.registry.update_runtime(job_id, execution_status=ExecutionStatus.QUEUED)
            self.registry.append_event(job_id, actor, "RESUME" if resume else "LAUNCH", "REQUESTED", "lease acquired")
            runtime = self.scheduler.launch(job_id, spec)
            job = self.registry.update_runtime(job_id, execution_status=ExecutionStatus.RUNNING, **runtime)
            self.registry.append_event(job_id, actor, "RESUME" if resume else "LAUNCH", "OK", runtime["unit_name"])
            return job
        except Exception as exc:
            if resource:
                self.registry.release_lease(str(resource), job_id)
            self.registry.update_runtime(job_id, execution_status=ExecutionStatus.FAILED, last_error=str(exc))
            self.registry.append_event(job_id, actor, "RESUME" if resume else "LAUNCH", "FAIL", str(exc))
            raise

    def pause(self, job_id: str, actor: str = "ui") -> dict[str, Any]:
        job = self.registry.get_job(job_id)
        if job["execution_status"] != ExecutionStatus.RUNNING.value:
            raise ServiceError(f"pause requires RUNNING, got {job['execution_status']}")
        if not job.get("control_token"):
            raise ServiceError("job has no cooperative control token")
        request_pause(Path(job["control_token"]), job_id, job["spec_hash"])
        updated = self.registry.update_runtime(job_id, execution_status=ExecutionStatus.PAUSE_REQUESTED)
        self.registry.append_event(job_id, actor, "PAUSE", "REQUESTED", "drain at next safe boundary")
        return updated

    def prepare_stop(self, job_id: str, actor: str = "ui", force: bool = False) -> dict[str, Any]:
        job = self.registry.get_job(job_id)
        if not job.get("unit_name"):
            raise ServiceError("job has no recorded systemd unit")
        action = "FORCE_STOP" if force else "STOP"
        nonce = self.registry.create_nonce(job_id, action)
        self.registry.append_event(job_id, actor, action, "PREPARED", "confirmation required")
        return {"job_id": job_id, "action": action, "nonce": nonce, "expires_in_s": 60}

    def confirm_stop(self, job_id: str, nonce: str, actor: str = "ui", force: bool = False) -> dict[str, Any]:
        action = "FORCE_STOP" if force else "STOP"
        if not self.registry.consume_nonce(nonce, job_id, action):
            raise ServiceError("invalid, expired, or consumed stop nonce")
        job = self.registry.get_job(job_id)
        try:
            self.scheduler.stop(str(job["unit_name"]), force=force)
            spec = JobSpec.from_dict(job["spec"])
            resource = spec.resources.get("gpu")
            if resource:
                self.registry.release_lease(str(resource), job_id)
            updated = self.registry.update_runtime(job_id, execution_status=ExecutionStatus.STOPPED)
            self.registry.append_event(job_id, actor, action, "OK", str(job["unit_name"]))
            return updated
        except Exception as exc:
            self.registry.append_event(job_id, actor, action, "FAIL", str(exc))
            raise

    def reconcile(self, actor: str = "system") -> list[dict[str, Any]]:
        reconciled = []
        for job in self.registry.list_jobs():
            spec = JobSpec.from_dict(job["spec"])
            observation = self.adapters.observe(spec)
            self.registry.update_evidence(
                job["job_id"],
                artifact_status=observation.get("artifact_status"),
                progress=observation.get("progress", {}),
            )
            current = job["execution_status"]
            if current not in {
                ExecutionStatus.QUEUED.value,
                ExecutionStatus.RUNNING.value,
                ExecutionStatus.PAUSE_REQUESTED.value,
            }:
                reconciled.append(self.registry.get_job(job["job_id"]))
                continue
            unit = job.get("unit_name")
            if unit and self.scheduler.active(unit):
                reconciled.append(self.registry.get_job(job["job_id"]))
                continue
            record = self.scheduler.wrapper_record(job.get("wrapper_record"))
            evidence = observation.get("execution_evidence")
            if record.get("status") == "PAUSED" or record.get("returncode") == 75:
                target = ExecutionStatus.PAUSED
            elif record.get("status") == "COMPLETE" or evidence == ExecutionStatus.COMPLETE:
                target = ExecutionStatus.COMPLETE
            elif record.get("status") == "FAILED":
                target = ExecutionStatus.FAILED
            else:
                target = ExecutionStatus.ORPHANED
            resource = spec.resources.get("gpu")
            if resource:
                self.registry.release_lease(str(resource), job["job_id"])
            self.registry.update_runtime(
                job["job_id"], execution_status=target, last_error=record.get("wrapper_error")
            )
            self.registry.append_event(
                job["job_id"], actor, "RECONCILE", "OK", f"{current} -> {target.value}", record
            )
            reconciled.append(self.registry.get_job(job["job_id"]))
        return reconciled

    def clone(self, job_id: str, overrides: dict[str, Any], actor: str = "ui") -> dict[str, Any]:
        source = self.registry.get_job(job_id)
        payload = dict(source["spec"])
        allowed = {"name", "campaign_id", "argv", "output_root", "declared_parameters", "effective_parameters", "parameter_sources", "input_digests", "environment", "source_key"}
        unknown = set(overrides) - allowed
        if unknown:
            raise ServiceError(f"clone overrides are not allowed: {sorted(unknown)}")
        payload.update(overrides)
        payload["readonly"] = False
        spec = JobSpec.from_dict(payload)
        self.adapters.validate(spec, allow_readonly=True)
        job = self.registry.create_job(spec, revision=int(source["spec_revision"]) + 1, parent_job_id=job_id)
        self.registry.append_event(job_id, actor, "CLONE", "OK", f"new draft {job['job_id']}")
        return job

    def read_log(self, job_id: str, offset: int = 0, limit: int = 131072, index: int = 0) -> dict[str, Any]:
        job = self.registry.get_job(job_id)
        spec = JobSpec.from_dict(job["spec"])
        candidates = ([job["primary_log"]] if job.get("primary_log") else []) + list(spec.log_paths)
        if not candidates or index < 0 or index >= len(candidates):
            return {"path": None, "offset": 0, "next_offset": 0, "text": ""}
        path = Path(candidates[index]).expanduser().resolve()
        allowed_roots = [settings.REPO_ROOT.resolve(), settings.state_root().resolve()]
        if not any(path == root or root in path.parents for root in allowed_roots):
            raise ServiceError("log path is outside allowed roots")
        if not path.is_file():
            return {"path": str(path), "offset": 0, "next_offset": 0, "text": ""}
        size = path.stat().st_size
        offset = max(0, min(int(offset), size))
        limit = max(1, min(int(limit), 1024 * 1024))
        with path.open("rb") as handle:
            handle.seek(offset)
            data = handle.read(limit)
        return {"path": str(path), "offset": offset, "next_offset": offset + len(data), "size": size, "text": data.decode("utf-8", errors="replace")}

    def snapshot(self, reconcile: bool = True) -> dict[str, Any]:
        if reconcile:
            self.reconcile()
        return {
            "schema": "updd_control_snapshot_v1",
            "generated_at": time.time(),
            "host": host.snapshot(),
            "jobs": self.registry.list_jobs(),
            "leases": self.registry.leases(),
            "events": self.registry.events(limit=200),
            "legacy_processes": legacy_processes(),
        }
