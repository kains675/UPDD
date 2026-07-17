"""Persistent systemd-owned wrapper that records child terminal evidence."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


PAUSE_EXIT_CODE = 75


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True)
    parser.add_argument("--token", required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--spec-hash", required=True)
    parser.add_argument("--input-digest", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = list(args.command)
    if command and command[0] == "--":
        command.pop(0)
    if not command or not Path(command[0]).is_absolute():
        parser.error("an absolute command is required after --")

    record_path = Path(args.record).expanduser().resolve()
    started = time.time()
    base = {
        "schema": "updd_control_wrapper_record_v1",
        "job_id": args.job_id,
        "spec_hash": args.spec_hash,
        "input_digest": args.input_digest,
        "command": command,
        "started_at": started,
        "wrapper_pid": os.getpid(),
        "status": "RUNNING",
    }
    atomic_json(record_path, base)
    env = os.environ.copy()
    env.update(
        UPDD_CONTROL_TOKEN=str(Path(args.token).expanduser().resolve()),
        UPDD_JOB_ID=args.job_id,
        UPDD_SPEC_HASH=args.spec_hash,
        UPDD_INPUT_DIGEST=args.input_digest,
    )
    try:
        result = subprocess.run(command, env=env, check=False)
        returncode = int(result.returncode)
        status = "PAUSED" if returncode == PAUSE_EXIT_CODE else (
            "COMPLETE" if returncode == 0 else "FAILED"
        )
        atomic_json(
            record_path,
            {
                **base,
                "status": status,
                "returncode": returncode,
                "ended_at": time.time(),
                "elapsed_s": time.time() - started,
            },
        )
        return returncode
    except BaseException as exc:
        atomic_json(
            record_path,
            {
                **base,
                "status": "FAILED",
                "returncode": None,
                "ended_at": time.time(),
                "elapsed_s": time.time() - started,
                "wrapper_error": f"{type(exc).__name__}: {exc}",
            },
        )
        raise


if __name__ == "__main__":
    raise SystemExit(main())
