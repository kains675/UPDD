"""Harmless boundary-aware fixture used by control-center integration tests."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
import time
from pathlib import Path

try:
    from .control_token import PAUSE_EXIT_CODE, boundary_pause_requested
except ImportError:
    from control_token import PAUSE_EXIT_CODE, boundary_pause_requested


def write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, indent=2)
            handle.write("\n")
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--units", type=int, default=3)
    parser.add_argument("--sleep", type=float, default=0.2)
    args = parser.parse_args()
    root = Path(args.output_root).resolve()
    progress_path = root / "progress.json"
    try:
        progress = json.loads(progress_path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        progress = {"completed": 0, "expected": args.units, "status": "RUNNING"}
    completed = int(progress.get("completed", 0))
    for unit in range(completed, args.units):
        time.sleep(args.sleep)
        completed = unit + 1
        progress = {"completed": completed, "expected": args.units, "status": "RUNNING"}
        write(progress_path, progress)
        print(f"self-test unit {completed}/{args.units}", flush=True)
        if boundary_pause_requested(progress):
            progress["status"] = "PAUSED"
            write(progress_path, progress)
            return PAUSE_EXIT_CODE
    progress["status"] = "COMPLETE"
    write(progress_path, progress)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
