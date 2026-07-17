"""Atomic cooperative pause tokens checked only at workflow boundaries."""

from __future__ import annotations

import fcntl
import json
import os
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator


PAUSE_EXIT_CODE = 75


def _now() -> float:
    return time.time()


@contextmanager
def _locked(path: Path) -> Iterator[None]:
    lock_path = path.with_suffix(path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _atomic_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def initialize(path: Path, job_id: str, spec_hash: str, input_digest: str) -> dict[str, Any]:
    payload = {
        "schema": "updd_control_token_v1",
        "job_id": job_id,
        "spec_hash": spec_hash,
        "input_digest": input_digest,
        "action": "RUN",
        "generation": 1,
        "updated_at": _now(),
        "acknowledged": None,
        "progress": None,
    }
    with _locked(path):
        if path.exists():
            current = read(path)
            for key in ("job_id", "spec_hash", "input_digest"):
                if current.get(key) != payload[key]:
                    raise ValueError(f"control token {key} mismatch")
            payload["generation"] = int(current.get("generation", 0)) + 1
        _atomic_write(path, payload)
    return payload


def read(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    if not isinstance(payload, dict):
        raise ValueError(f"invalid control token: {path}")
    return payload


def request_pause(path: Path, job_id: str, spec_hash: str) -> dict[str, Any]:
    with _locked(path):
        payload = read(path)
        if payload.get("job_id") != job_id or payload.get("spec_hash") != spec_hash:
            raise ValueError("control token identity mismatch")
        payload.update(
            action="PAUSE_REQUESTED",
            generation=int(payload.get("generation", 0)) + 1,
            updated_at=_now(),
            acknowledged=None,
        )
        _atomic_write(path, payload)
        return payload


def pause_pending() -> bool:
    """Read a control request without acknowledging the safe boundary yet."""
    token_name = os.environ.get("UPDD_CONTROL_TOKEN")
    if not token_name:
        return False
    path = Path(token_name).expanduser().resolve()
    payload = read(path)
    if (
        payload.get("job_id") != os.environ.get("UPDD_JOB_ID")
        or payload.get("spec_hash") != os.environ.get("UPDD_SPEC_HASH")
    ):
        raise RuntimeError("control token identity/spec mismatch")
    return payload.get("action") == "PAUSE_REQUESTED"


def boundary_pause_requested(progress: dict[str, Any] | None = None) -> bool:
    token_name = os.environ.get("UPDD_CONTROL_TOKEN")
    if not token_name:
        return False
    path = Path(token_name).expanduser().resolve()
    expected_job = os.environ.get("UPDD_JOB_ID")
    expected_hash = os.environ.get("UPDD_SPEC_HASH")
    with _locked(path):
        payload = read(path)
        if payload.get("job_id") != expected_job or payload.get("spec_hash") != expected_hash:
            raise RuntimeError("control token identity/spec mismatch at boundary")
        if payload.get("action") != "PAUSE_REQUESTED":
            return False
        payload.update(
            acknowledged="PAUSED",
            acknowledged_at=_now(),
            updated_at=_now(),
            progress=progress,
        )
        _atomic_write(path, payload)
    return True
