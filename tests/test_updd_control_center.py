from __future__ import annotations

import json
import os
import sqlite3
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from utils.control_center import settings
from utils.control_center.adapters import AdapterError, AdapterManager
from utils.control_center.api import Server
from utils.control_center.control_token import (
    boundary_pause_requested,
    initialize,
    pause_pending,
    read,
    request_pause,
)
from utils.control_center.models import (
    ArtifactStatus,
    ExecutionStatus,
    JobSpec,
    ScientificStatus,
)
from utils.control_center.registry import LeaseConflict, Registry, RegistryError
from utils.control_center.scheduler import SystemdScheduler
from utils.control_center.service import ControlCenter, ServiceError


def self_test_spec(output_root: Path, *, name: str = "fixture", units: int = 8, sleep: float = 0.1) -> JobSpec:
    return JobSpec(
        name=name,
        campaign_id="control-center-tests",
        adapter="self_test",
        interpreter=str(settings.QMMM_PYTHON),
        script=str(settings.REPO_ROOT / "utils/control_center/self_test_job.py"),
        argv=("--output-root", str(output_root), "--units", str(units), "--sleep", str(sleep)),
        cwd=str(settings.REPO_ROOT),
        output_root=str(output_root),
        resources={"gpu": "selftest_gpu"},
        declared_parameters={"units": units, "sleep": sleep},
        effective_parameters={"units": units, "sleep": sleep},
        source_key=f"preregistered:selftest:{name}",
    )


def wait_until(predicate, timeout: float = 10.0, interval: float = 0.05):
    deadline = time.time() + timeout
    while time.time() < deadline:
        value = predicate()
        if value:
            return value
        time.sleep(interval)
    raise AssertionError("condition did not become true before timeout")


def test_job_spec_hash_is_canonical_and_rejects_unlisted_paths(tmp_path):
    first = self_test_spec(tmp_path / "a", name="canonical")
    second = JobSpec.from_dict(first.to_dict())
    assert first.spec_hash == second.spec_hash
    assert first.input_digest == second.input_digest
    payload = first.to_dict()
    payload["interpreter"] = "/bin/bash"
    with pytest.raises(ValueError, match="interpreter is not allow-listed"):
        JobSpec.from_dict(payload)


def test_job_spec_rejects_unlisted_environment(tmp_path):
    payload = self_test_spec(tmp_path / "a", name="env").to_dict()
    payload["environment"] = {"LD_PRELOAD": "/tmp/inject.so"}
    with pytest.raises(ValueError, match="environment keys"):
        JobSpec.from_dict(payload)


def test_registry_keeps_four_status_axes_independent(tmp_path):
    registry = Registry(tmp_path / "registry.sqlite3")
    job = registry.create_job(self_test_spec(tmp_path / "out", name="axes"))
    registry.update_runtime(job["job_id"], execution_status=ExecutionStatus.COMPLETE)
    result = registry.update_evidence(
        job["job_id"],
        artifact_status=ArtifactStatus.PRESENT,
        scientific_status=ScientificStatus.INVALID,
        ranking_eligible=False,
    )
    assert result["execution_status"] == "COMPLETE"
    assert result["artifact_status"] == "PRESENT"
    assert result["scientific_status"] == "INVALID"
    assert result["ranking_eligible"] is False


def test_events_are_append_only_at_database_level(tmp_path):
    registry = Registry(tmp_path / "registry.sqlite3")
    job = registry.create_job(self_test_spec(tmp_path / "out", name="events"))
    event_id = registry.append_event(job["job_id"], "test", "ACTION", "OK")
    con = sqlite3.connect(registry.path)
    with pytest.raises(sqlite3.IntegrityError, match="append-only"):
        con.execute("UPDATE events SET result='changed' WHERE event_id=?", (event_id,))
    con.close()


def test_job_spec_is_immutable_at_database_level(tmp_path):
    registry = Registry(tmp_path / "registry.sqlite3")
    job = registry.create_job(self_test_spec(tmp_path / "out", name="immutable"))
    con = sqlite3.connect(registry.path)
    with pytest.raises(sqlite3.IntegrityError, match="job spec is immutable"):
        con.execute("UPDATE jobs SET spec_json='{}' WHERE job_id=?", (job["job_id"],))
    con.close()


def test_controllable_source_key_cannot_silently_change_spec(tmp_path):
    registry = Registry(tmp_path / "registry.sqlite3")
    first = self_test_spec(tmp_path / "out", name="source-key", units=2)
    registry.create_job(first)
    changed = JobSpec.from_dict({**first.to_dict(), "argv": ["--output-root", str(tmp_path / "out"), "--units", "3"]})
    with pytest.raises(RegistryError, match="different immutable spec"):
        registry.create_job(changed)


def test_resource_lease_is_transactionally_exclusive(tmp_path):
    registry = Registry(tmp_path / "registry.sqlite3")
    one = registry.create_job(self_test_spec(tmp_path / "one", name="lease-one"))
    two = registry.create_job(self_test_spec(tmp_path / "two", name="lease-two"))
    registry.acquire_lease("host_cuda0", one["job_id"])
    with pytest.raises(LeaseConflict):
        registry.acquire_lease("host_cuda0", two["job_id"])
    assert registry.release_lease("host_cuda0", two["job_id"]) is False
    assert registry.release_lease("host_cuda0", one["job_id"]) is True


def test_stop_nonce_is_job_bound_expiring_and_single_use(tmp_path):
    registry = Registry(tmp_path / "registry.sqlite3")
    one = registry.create_job(self_test_spec(tmp_path / "one", name="nonce-one"))
    two = registry.create_job(self_test_spec(tmp_path / "two", name="nonce-two"))
    nonce = registry.create_nonce(one["job_id"], "STOP", ttl_s=60)
    assert registry.consume_nonce(nonce, two["job_id"], "STOP") is False
    assert registry.consume_nonce(nonce, one["job_id"], "FORCE_STOP") is False
    assert registry.consume_nonce(nonce, one["job_id"], "STOP") is True
    assert registry.consume_nonce(nonce, one["job_id"], "STOP") is False


def test_control_token_requires_identity_and_acknowledges_boundary(tmp_path, monkeypatch):
    path = tmp_path / "control.json"
    initialize(path, "job-1", "a" * 64, "b" * 64)
    request_pause(path, "job-1", "a" * 64)
    monkeypatch.setenv("UPDD_CONTROL_TOKEN", str(path))
    monkeypatch.setenv("UPDD_JOB_ID", "job-1")
    monkeypatch.setenv("UPDD_SPEC_HASH", "a" * 64)
    assert pause_pending() is True
    assert boundary_pause_requested({"completed": 3}) is True
    assert read(path)["acknowledged"] == "PAUSED"
    monkeypatch.setenv("UPDD_SPEC_HASH", "c" * 64)
    with pytest.raises(RuntimeError, match="identity/spec"):
        pause_pending()


def test_real_discovery_preserves_execution_science_separation(tmp_path):
    control = ControlCenter(
        Registry(tmp_path / "registry.sqlite3"),
        SystemdScheduler(tmp_path / "state"),
        AdapterManager(settings.REPO_ROOT),
    )
    jobs = {job["adapter"]: job for job in control.discover() if job["adapter"] != "tracka"}
    assert jobs["dcd_sidecar"]["execution_status"] == "COMPLETE"
    assert jobs["dcd_sidecar"]["scientific_status"] == "VALID"
    assert jobs["dcd_sidecar"]["ranking_eligible"] is False
    assert jobs["trackb_pool"]["execution_status"] == "COMPLETE"
    assert jobs["trackb_pool"]["artifact_status"] == "PRESENT"
    assert jobs["trackb_pool"]["scientific_status"] == "INVALID"
    assert jobs["trackb_pool"]["ranking_eligible"] is False
    assert (
        jobs["scaffold_md"]["progress"]["current_runner_revalidation"]
        == "BLOCKED_CURRENT_RUNNER_DRIFT"
    )


def test_readonly_discovered_job_cannot_preflight(tmp_path):
    control = ControlCenter(Registry(tmp_path / "db.sqlite3"), SystemdScheduler(tmp_path / "state"))
    job = next(row for row in control.discover() if row["adapter"] == "dcd_sidecar")
    with pytest.raises(Exception, match="read-only"):
        control.preflight(job["job_id"])


def test_launch_requires_successful_preflight(tmp_path, monkeypatch):
    monkeypatch.setenv("UPDD_CC_ENABLE_SELF_TEST", "1")
    control = ControlCenter(Registry(tmp_path / "db.sqlite3"), SystemdScheduler(tmp_path / "state"))
    job = control.register(self_test_spec(tmp_path / "out", name="preflight-required").to_dict())
    with pytest.raises(ServiceError, match="cannot launch from DRAFT"):
        control.launch(job["job_id"])


def test_scientific_registration_requires_exact_catalog_spec(tmp_path):
    catalog = tmp_path / "catalog"
    catalog.mkdir()
    manager = AdapterManager(settings.REPO_ROOT, preregistered_root=catalog)
    control = ControlCenter(
        Registry(tmp_path / "db.sqlite3"),
        SystemdScheduler(tmp_path / "state"),
        manager,
    )
    discovered = next(item for item in manager.discover() if item.spec.adapter == "dcd_sidecar")
    payload = discovered.spec.to_dict()
    payload["readonly"] = False
    payload["source_key"] = "preregistered:dcd_fixture.json"
    payload["input_digests"]["preregistration_sha256"] = "f" * 64
    declared = JobSpec.from_dict(payload)
    (catalog / "dcd_fixture.json").write_text(
        json.dumps(declared.to_dict(), sort_keys=True), encoding="utf-8"
    )

    accepted = control.register(declared.to_dict(), actor="test")
    assert accepted["spec_hash"] == declared.spec_hash

    changed = declared.to_dict()
    changed["argv"] = ["--analyze-only"]
    with pytest.raises(AdapterError, match="exact preregistered spec"):
        control.register(changed, actor="test")


def test_log_reads_are_offset_based_and_path_confined(tmp_path, monkeypatch):
    monkeypatch.setenv("UPDD_CC_ENABLE_SELF_TEST", "1")
    monkeypatch.setenv("UPDD_CC_STATE_ROOT", str(tmp_path / "state"))
    registry = Registry(tmp_path / "db.sqlite3")
    control = ControlCenter(registry, SystemdScheduler(tmp_path / "state"))
    spec = self_test_spec(tmp_path / "out", name="logs")
    job = registry.create_job(spec)
    log = tmp_path / "state/jobs" / job["job_id"] / "job.log"
    log.parent.mkdir(parents=True)
    log.write_text("first\nsecond\n", encoding="utf-8")
    registry.update_runtime(job["job_id"], primary_log=str(log))
    first = control.read_log(job["job_id"], offset=0, limit=6)
    second = control.read_log(job["job_id"], offset=first["next_offset"], limit=100)
    assert first["text"] == "first\n"
    assert second["text"] == "second\n"


def test_api_requires_bearer_token(tmp_path):
    control = ControlCenter(Registry(tmp_path / "db.sqlite3"), SystemdScheduler(tmp_path / "state"))
    server = Server(("127.0.0.1", 0), control, "secret")
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{server.server_port}"
    try:
        assert json.loads(urllib.request.urlopen(base + "/health").read())["status"] == "ok"
        with pytest.raises(urllib.error.HTTPError) as denied:
            urllib.request.urlopen(base + "/v1/jobs")
        assert denied.value.code == 401
        request = urllib.request.Request(base + "/v1/jobs", headers={"Authorization": "Bearer secret"})
        assert "jobs" in json.loads(urllib.request.urlopen(request).read())
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def test_systemd_job_survives_service_restart_pauses_and_resumes(tmp_path, monkeypatch):
    monkeypatch.setenv("UPDD_CC_ENABLE_SELF_TEST", "1")
    registry_path = tmp_path / "db.sqlite3"
    state = tmp_path / "state"
    scheduler = SystemdScheduler(state)
    first = ControlCenter(Registry(registry_path), scheduler)
    spec = self_test_spec(tmp_path / "out", name="restart-pause", units=12, sleep=0.12)
    job = first.register(spec.to_dict(), actor="test")
    first.preflight(job["job_id"], actor="test")
    running = first.launch(job["job_id"], actor="test")
    unit = running["unit_name"]
    try:
        wait_until(lambda: (json.loads((tmp_path / "out/progress.json").read_text())["completed"] >= 1) if (tmp_path / "out/progress.json").exists() else False)
        restarted = ControlCenter(Registry(registry_path), SystemdScheduler(state))
        assert restarted.scheduler.active(unit)
        assert restarted.reconcile(actor="restart")[0]["execution_status"] == "RUNNING"
        restarted.pause(job["job_id"], actor="test")
        wait_until(lambda: restarted.reconcile(actor="poll")[0]["execution_status"] == "PAUSED")
        paused = restarted.registry.get_job(job["job_id"])
        partial = paused["progress"]["completed"]
        assert 1 <= partial < 12
        assert restarted.registry.leases() == []
        restarted.launch(job["job_id"], actor="test", resume=True)
        wait_until(lambda: restarted.reconcile(actor="poll")[0]["execution_status"] == "COMPLETE", timeout=10)
        complete = restarted.registry.get_job(job["job_id"])
        assert complete["progress"]["completed"] == 12
        assert complete["artifact_status"] == "PRESENT"
        assert complete["scientific_status"] == "UNKNOWN"
    finally:
        current = first.registry.get_job(job["job_id"])
        if current.get("unit_name") and scheduler.active(current["unit_name"]):
            scheduler.stop(current["unit_name"], force=True)


def test_reconcile_after_runtime_loss_uses_persistent_evidence(tmp_path, monkeypatch):
    monkeypatch.setenv("UPDD_CC_ENABLE_SELF_TEST", "1")

    class InactiveScheduler(SystemdScheduler):
        @classmethod
        def active(cls, unit: str) -> bool:
            return False

    registry = Registry(tmp_path / "db.sqlite3")
    control = ControlCenter(registry, InactiveScheduler(tmp_path / "state"))

    def staged(name: str, status: str, completed: int, record: dict | None = None):
        root = tmp_path / name
        root.mkdir()
        (root / "progress.json").write_text(
            json.dumps({"status": status, "completed": completed, "expected": 4}),
            encoding="utf-8",
        )
        job = registry.create_job(self_test_spec(root, name=name, units=4))
        record_path = tmp_path / f"{name}.wrapper.json"
        if record is not None:
            record_path.write_text(json.dumps(record), encoding="utf-8")
        registry.update_runtime(
            job["job_id"],
            execution_status=ExecutionStatus.RUNNING,
            unit_name=f"updd-job-{name}",
            wrapper_record=str(record_path),
        )
        return job["job_id"]

    complete_id = staged("complete-evidence", "COMPLETE", 4)
    paused_id = staged("paused-record", "RUNNING", 2, {"status": "PAUSED", "returncode": 75})
    orphaned_id = staged("orphaned-no-record", "RUNNING", 1)

    rows = {row["job_id"]: row for row in control.reconcile(actor="startup-test")}
    assert rows[complete_id]["execution_status"] == "COMPLETE"
    assert rows[paused_id]["execution_status"] == "PAUSED"
    assert rows[orphaned_id]["execution_status"] == "ORPHANED"
    actions = [event["action"] for event in registry.events(limit=100)]
    assert actions.count("RECONCILE") == 3


def test_service_blocks_second_job_on_same_resource(tmp_path, monkeypatch):
    monkeypatch.setenv("UPDD_CC_ENABLE_SELF_TEST", "1")
    control = ControlCenter(Registry(tmp_path / "db.sqlite3"), SystemdScheduler(tmp_path / "state"))
    first = control.register(
        self_test_spec(tmp_path / "one", name="gpu-one", units=100, sleep=0.05).to_dict(),
        actor="test",
    )
    second = control.register(
        self_test_spec(tmp_path / "two", name="gpu-two", units=2, sleep=0.05).to_dict(),
        actor="test",
    )
    control.preflight(first["job_id"], actor="test")
    control.preflight(second["job_id"], actor="test")
    running = control.launch(first["job_id"], actor="test")
    try:
        wait_until(lambda: control.scheduler.active(running["unit_name"]))
        with pytest.raises(LeaseConflict, match="selftest_gpu is held"):
            control.launch(second["job_id"], actor="test")
        assert control.registry.leases()[0]["job_id"] == first["job_id"]
    finally:
        if control.scheduler.active(running["unit_name"]):
            prepared = control.prepare_stop(first["job_id"], actor="test", force=True)
            control.confirm_stop(
                first["job_id"], prepared["nonce"], actor="test", force=True
            )


def test_systemd_exact_stop_requires_confirmation(tmp_path, monkeypatch):
    monkeypatch.setenv("UPDD_CC_ENABLE_SELF_TEST", "1")
    control = ControlCenter(Registry(tmp_path / "db.sqlite3"), SystemdScheduler(tmp_path / "state"))
    spec = self_test_spec(tmp_path / "out", name="exact-stop", units=100, sleep=0.1)
    job = control.register(spec.to_dict(), actor="test")
    control.preflight(job["job_id"], actor="test")
    running = control.launch(job["job_id"], actor="test")
    try:
        wait_until(lambda: control.scheduler.active(running["unit_name"]))
        prepared = control.prepare_stop(job["job_id"], actor="test")
        with pytest.raises(ServiceError, match="invalid"):
            control.confirm_stop(job["job_id"], "wrong", actor="test")
        stopped = control.confirm_stop(job["job_id"], prepared["nonce"], actor="test")
        assert stopped["execution_status"] == "STOPPED"
        wait_until(lambda: not control.scheduler.active(running["unit_name"]))
        assert control.registry.leases() == []
    finally:
        if control.scheduler.active(running["unit_name"]):
            control.scheduler.stop(running["unit_name"], force=True)


def test_no_shell_or_broad_kill_surface_in_control_modules():
    roots = [settings.REPO_ROOT / "utils/control_center", settings.REPO_ROOT / "utils/updd_control_center.py"]
    text = "\n".join(
        path.read_text(encoding="utf-8")
        for root in roots
        for path in ([root] if root.is_file() else root.glob("*.py"))
    )
    assert "shell=True" not in text
    assert "pkill" not in text
    assert "git commit" not in text
    assert "git push" not in text
