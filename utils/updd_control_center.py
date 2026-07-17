#!/usr/bin/env python3
"""CLI entry point for the UPDD local control center."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utils.control_center import settings
from utils.control_center.api import main as api_main
from utils.control_center.registry import Registry
from utils.control_center.service import ControlCenter


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    serve = sub.add_parser("serve-api")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8765)
    serve.add_argument("--database", type=Path, default=settings.database_path())
    serve.add_argument("--token-file", type=Path, default=settings.api_token_path())
    sub.add_parser("discover")
    sub.add_parser("reconcile")
    sub.add_parser("snapshot")
    args = parser.parse_args(argv)
    if args.command == "serve-api":
        return api_main(
            [
                "--host", args.host,
                "--port", str(args.port),
                "--database", str(args.database),
                "--token-file", str(args.token_file),
            ]
        )
    control = ControlCenter(registry=Registry(settings.database_path()))
    if args.command == "discover":
        payload = {"jobs": control.discover()}
    elif args.command == "reconcile":
        payload = {"jobs": control.reconcile(actor="cli")}
    else:
        payload = control.snapshot()
    print(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
