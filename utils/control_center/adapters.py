"""Allow-listed workflow adapters and read-only legacy discovery."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import psutil

from . import settings
from .models import ArtifactStatus, ExecutionStatus, JobSpec, ScientificStatus


def load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class Discovery:
    spec: JobSpec
    execution_status: ExecutionStatus
    artifact_status: ArtifactStatus
    scientific_status: ScientificStatus
    ranking_eligible: bool | None
    progress: dict[str, Any]


class AdapterError(RuntimeError):
    pass


class AdapterManager:
    def __init__(
        self,
        repo_root: Path = settings.REPO_ROOT,
        preregistered_root: Path | None = None,
    ):
        self.repo_root = repo_root.resolve()
        self.preregistered_root = (
            preregistered_root or self.repo_root / "analysis"
        ).resolve()

    def _scaffold_runner_revalidation(self) -> dict[str, str]:
        analysis = self.repo_root / "analysis/1ycr_wt_seed_expansion_20260716"
        inventory = load_json(analysis / "source_inventory.json")
        expected = str((inventory.get("runner") or {}).get("sha256", ""))
        runner = analysis / "expand_1ycr_wt_scaffolds.py"
        current = sha256(runner) if runner.is_file() else ""
        status = (
            "PASS"
            if expected and current == expected
            else "BLOCKED_CURRENT_RUNNER_DRIFT"
        )
        return {
            "current_runner_revalidation": status,
            "frozen_runner_sha256": expected,
            "current_runner_sha256": current,
        }

    def _validate_preregistered_spec(self, spec: JobSpec) -> None:
        prefix = "preregistered:"
        source_key = spec.source_key or ""
        if not source_key.startswith(prefix):
            raise AdapterError("controllable jobs require a preregistered source_key")
        relative = source_key[len(prefix):]
        if not relative or Path(relative).is_absolute():
            raise AdapterError("preregistered source_key must name a relative JSON file")
        declared_path = (self.preregistered_root / relative).resolve()
        if (
            declared_path.suffix != ".json"
            or not (
                declared_path == self.preregistered_root
                or self.preregistered_root in declared_path.parents
            )
        ):
            raise AdapterError("preregistered spec path escapes the catalog root")
        declared_payload = load_json(declared_path)
        if not declared_payload:
            raise AdapterError(f"preregistered spec is missing or invalid: {declared_path}")
        try:
            declared = JobSpec.from_dict(declared_payload)
        except (KeyError, TypeError, ValueError) as exc:
            raise AdapterError(f"invalid preregistered spec: {declared_path}: {exc}") from exc
        if declared.spec_hash != spec.spec_hash:
            raise AdapterError("registered job does not match the exact preregistered spec")

    def validate(self, spec: JobSpec, *, allow_readonly: bool = False) -> dict[str, Any]:
        if spec.readonly and not allow_readonly:
            raise AdapterError("read-only imported jobs cannot be controlled")
        for label, path in (
            ("interpreter", Path(spec.interpreter)),
            ("script", Path(spec.script)),
            ("cwd", Path(spec.cwd)),
        ):
            if not path.exists():
                raise AdapterError(f"{label} does not exist: {path}")
        if spec.adapter == "self_test":
            if os.environ.get("UPDD_CC_ENABLE_SELF_TEST") != "1":
                raise AdapterError("self-test adapter is disabled")
            return {"status": "PASS", "mode": "self_test"}
        if not spec.readonly:
            self._validate_preregistered_spec(spec)
        prereg = spec.input_digests.get("preregistration_sha256")
        if not spec.readonly and not prereg:
            raise AdapterError("controllable scientific jobs require a preregistration digest")
        for name, expected in spec.input_digests.items():
            if name.endswith("_path"):
                continue
            path_key = f"{name.removesuffix('_sha256')}_path"
            path_name = spec.input_digests.get(path_key)
            if name.endswith("_sha256") and path_name:
                path = Path(path_name)
                if not path.is_file() or sha256(path) != expected:
                    raise AdapterError(f"input digest mismatch: {name}")
        if spec.adapter in {"dcd_sidecar", "scaffold_md"}:
            inventory = load_json(
                self.repo_root
                / (
                    "analysis/w4a_postdensify_dcd_20260716/source_inventory.json"
                    if spec.adapter == "dcd_sidecar"
                    else "analysis/1ycr_wt_seed_expansion_20260716/source_inventory.json"
                )
            )
            expected = spec.input_digests.get("source_inventory")
            if expected and inventory.get("inventory_digest") != expected:
                raise AdapterError("frozen source inventory digest mismatch")
        return {
            "status": "PASS",
            "adapter": spec.adapter,
            "spec_hash": spec.spec_hash,
            "input_digest": spec.input_digest,
            "command": spec.command,
        }

    def observe(self, spec: JobSpec) -> dict[str, Any]:
        root = Path(spec.output_root) if spec.output_root else None
        if spec.adapter in {"dcd_sidecar", "scaffold_md", "self_test"} and root:
            state = load_json(root / "run_state.json")
            progress = load_json(root / "progress.json")
            if spec.adapter == "scaffold_md":
                progress = {**progress, **self._scaffold_runner_revalidation()}
            if spec.adapter == "self_test":
                progress = load_json(root / "progress.json")
                state = {"status": progress.get("status")}
            status = str(state.get("status", progress.get("status", ""))).upper()
            execution = {
                "COMPLETE": ExecutionStatus.COMPLETE,
                "PAUSED": ExecutionStatus.PAUSED,
                "FAILED": ExecutionStatus.FAILED,
            }.get(status)
            expected = int(progress.get("n_expected", progress.get("expected", 0)) or 0)
            completed = int(progress.get("n_completed", progress.get("completed", 0)) or 0)
            artifact = ArtifactStatus.MISSING
            if completed:
                artifact = ArtifactStatus.PRESENT if expected and completed == expected else ArtifactStatus.PARTIAL
            return {"execution_evidence": execution, "artifact_status": artifact, "progress": progress}
        if spec.adapter == "trackb_pool" and root:
            manifests = [
                path
                for path in root.glob("*/wt/*/rep*/run_manifest.json")
                if "_archive" not in path.parts
            ]
            expected = int(spec.effective_parameters.get("expected_replicates", 0) or 0)
            completed = len(manifests)
            return {
                "execution_evidence": ExecutionStatus.COMPLETE if expected and completed == expected else None,
                "artifact_status": ArtifactStatus.PRESENT if expected and completed == expected else (
                    ArtifactStatus.PARTIAL if completed else ArtifactStatus.MISSING
                ),
                "progress": {"completed": completed, "expected": expected, "eta_s": 0 if completed == expected else None},
            }
        if spec.adapter == "tracka" and root:
            status = load_json(root / "updd_status.json")
            completed = status.get("completed_steps", [])
            return {
                "execution_evidence": None,
                "artifact_status": ArtifactStatus.PARTIAL if status else ArtifactStatus.MISSING,
                "progress": {"completed_steps": completed, "raw": status},
            }
        return {"execution_evidence": None, "artifact_status": ArtifactStatus.MISSING, "progress": {}}

    def discover(self) -> list[Discovery]:
        found: list[Discovery] = []
        dcd = self._discover_dcd()
        if dcd:
            found.append(dcd)
        scaffold = self._discover_scaffold()
        if scaffold:
            found.append(scaffold)
        trackb = self._discover_h18()
        if trackb:
            found.append(trackb)
        found.extend(self._discover_tracka())
        return found

    def _discover_dcd(self) -> Discovery | None:
        analysis = self.repo_root / "analysis/w4a_postdensify_dcd_20260716"
        root = self.repo_root / "outputs/_trackb/w4a_postdensify_dcd_20260716"
        state = load_json(root / "run_state.json")
        if not state:
            return None
        progress = load_json(root / "progress.json")
        inventory = load_json(analysis / "source_inventory.json")
        gate = load_json(analysis / "postdensify_gate.json")
        status = ExecutionStatus.COMPLETE if state.get("status") == "COMPLETE" else ExecutionStatus.ORPHANED
        complete = int(state.get("n_completed", 0)) == int(state.get("n_expected", -1))
        spec = JobSpec(
            name="W4A post-densify DCD structural gate",
            campaign_id="w4a-postdensify-a-gate-20260716",
            adapter="dcd_sidecar",
            interpreter=str(settings.ATM_PYTHON),
            script=str(analysis / "postdensify_dcd.py"),
            cwd=str(self.repo_root),
            output_root=str(root),
            resources={"gpu": "host_cuda0", "ram_gb": 16, "cpu_slots": 4},
            declared_parameters={
                "cells": 20,
                "cycles": 50,
                "md_steps_per_cycle": 250,
                "site_distance_nm": 0.26,
                "control_floor": 0.15,
                "amplification_delta": 0.10,
                "amplification_ratio": 1.5,
            },
            effective_parameters={
                "cells": int(state.get("n_expected", 20)),
                "gate_verdict": gate.get("verdict", "UNKNOWN"),
                "branch": "STOP_B_PROMOTE_A" if gate.get("verdict") == "PROMOTE_A_CARVE_REDESIGN" else "UNKNOWN",
            },
            parameter_sources={
                "declared": str(analysis / "PREREGISTRATION.md"),
                "effective": str(analysis / "postdensify_gate.json"),
            },
            input_digests={"source_inventory": str(inventory.get("inventory_digest", ""))},
            log_paths=(str(root / "postdensify_dcd.log"), str(analysis / "postdensify_gate.md")),
            readonly=True,
            source_key="discovered:w4a-postdensify-dcd-20260716",
        )
        return Discovery(
            spec,
            status,
            ArtifactStatus.PRESENT if complete and gate else ArtifactStatus.PARTIAL,
            ScientificStatus.VALID if gate.get("verdict") else ScientificStatus.UNKNOWN,
            False,
            progress,
        )

    def _discover_scaffold(self) -> Discovery | None:
        analysis = self.repo_root / "analysis/1ycr_wt_seed_expansion_20260716"
        root = self.repo_root / "outputs/_trackb/1ycr_wt_seed_expansion_20260716"
        state = load_json(root / "run_state.json")
        if not state:
            return None
        progress = load_json(root / "progress.json")
        inventory = load_json(analysis / "source_inventory.json")
        runner_revalidation = self._scaffold_runner_revalidation()
        progress = {**progress, **runner_revalidation}
        complete = int(state.get("n_completed", 0)) == int(state.get("n_expected", -1))
        spec = JobSpec(
            name="1YCR WT scaffold seed expansion",
            campaign_id="1ycr-wt-scaffold-n10-20260716",
            adapter="scaffold_md",
            interpreter=str(settings.QMMM_PYTHON),
            script=str(analysis / "expand_1ycr_wt_scaffolds.py"),
            cwd=str(self.repo_root),
            output_root=str(root),
            resources={"gpu": "host_cuda0", "ram_gb": 8, "cpu_slots": 2},
            declared_parameters={"new_seeds": [83, 127, 163, 199, 251], "steps": 2500000, "dt_fs": 2.0},
            effective_parameters={
                "completed_seeds": progress.get("completed_seeds", []),
                "status": state.get("status"),
                "historical_artifact_validation": "PASS" if complete else "UNKNOWN",
                **runner_revalidation,
            },
            parameter_sources={"declared": str(analysis / "PREREGISTRATION.md"), "effective": str(root / "progress.json")},
            input_digests={"source_inventory": str(inventory.get("inventory_digest", ""))},
            log_paths=(str(root / "queue.log"), str(analysis / "scaffold_expansion_report.md")),
            readonly=True,
            source_key="discovered:1ycr-wt-seed-expansion-20260716",
        )
        return Discovery(
            spec,
            ExecutionStatus.COMPLETE if complete else ExecutionStatus.ORPHANED,
            ArtifactStatus.PRESENT if complete else ArtifactStatus.PARTIAL,
            ScientificStatus.VALID if complete else ScientificStatus.UNKNOWN,
            False,
            progress,
        )

    def _discover_h18(self) -> Discovery | None:
        root = self.repo_root / "outputs/_trackb/w4a_bdcd_h18_20260710_rerun1"
        preregs = sorted(root.glob("*/pre_registration.json"))
        if len(preregs) != 4:
            return None
        configs = {path.parent.name: load_json(path).get("config", {}) for path in preregs}
        manifests = [
            path for path in root.glob("*/wt/*/rep*/run_manifest.json") if "_archive" not in path.parts
        ]
        logs = sorted(root.glob("*.log"), key=lambda path: path.stat().st_mtime, reverse=True)
        spec = JobSpec(
            name="W4A H18 carved/uncarved matched cohort",
            campaign_id="w4a-h18-densify-n10-20260710",
            adapter="trackb_pool",
            interpreter=str(settings.ATM_PYTHON),
            script=str(self.repo_root / "scripts/trackb_inplace_rbfe_production.py"),
            cwd=str(self.repo_root),
            output_root=str(root),
            resources={"gpu": "host_cuda0", "ram_gb": 24, "cpu_slots": 4},
            declared_parameters={name: config for name, config in configs.items()},
            effective_parameters={
                "expected_replicates": 40,
                "completed_replicates": len(manifests),
                "hardened_transport_pass": "0/80 directions",
                "branch_decision": "STOP_B_TRANSPORT_PROMOTE_A",
            },
            parameter_sources={"declared": [str(path) for path in preregs], "science": "post-run hardened mixing audit"},
            input_digests={f"prereg_{path.parent.name}_sha256": sha256(path) for path in preregs},
            log_paths=tuple(str(path) for path in logs[:3]),
            readonly=True,
            source_key="discovered:w4a-h18-rerun1-20260710",
        )
        complete = len(manifests) == 40
        return Discovery(
            spec,
            ExecutionStatus.COMPLETE if complete else ExecutionStatus.ORPHANED,
            ArtifactStatus.PRESENT if complete else ArtifactStatus.PARTIAL,
            ScientificStatus.INVALID if complete else ScientificStatus.UNKNOWN,
            False,
            {
                "completed": len(manifests),
                "expected": 40,
                "eta_s": 0 if complete else None,
                "scientific_gate": "HARDENED_TRANSPORT_FAIL_0_OF_80" if complete else "UNKNOWN",
            },
        )

    def _discover_tracka(self) -> list[Discovery]:
        rows: list[Discovery] = []
        for path in sorted((self.repo_root / "outputs").glob("*/updd_status.json")):
            payload = load_json(path)
            if not payload:
                continue
            root = path.parent
            completed = payload.get("completed_steps", [])
            spec = JobSpec(
                name=f"Track A legacy: {root.name}",
                campaign_id=f"tracka-legacy-{root.name}",
                adapter="tracka",
                interpreter=str(settings.QMMM_PYTHON),
                script=str(self.repo_root / "UPDD.py"),
                cwd=str(self.repo_root),
                output_root=str(root),
                resources={"gpu": "host_cuda0"},
                effective_parameters={"completed_steps": completed},
                parameter_sources={"effective": str(path)},
                log_paths=(),
                readonly=True,
                source_key=f"discovered:tracka:{root.name}",
            )
            rows.append(
                Discovery(
                    spec,
                    ExecutionStatus.ORPHANED,
                    ArtifactStatus.PARTIAL,
                    ScientificStatus.UNKNOWN,
                    None,
                    {"completed_steps": completed},
                )
            )
        return rows


def legacy_processes() -> list[dict[str, Any]]:
    needles = ("UPDD.py", "trackb_inplace_rbfe_production.py", "postdensify_dcd.py", "run_restrained_md.py")
    rows = []
    for process in psutil.process_iter(["pid", "name", "cmdline", "create_time"]):
        try:
            command = " ".join(process.info.get("cmdline") or [])
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
        if command and any(needle in command for needle in needles):
            rows.append({"pid": process.pid, "name": process.info.get("name"), "command": command, "create_time": process.info.get("create_time"), "controllable": False})
    return rows
