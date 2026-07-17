"""Authenticated localhost JSON API implemented with the standard library."""

from __future__ import annotations

import argparse
import json
import os
import secrets
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

from . import settings
from .registry import Registry
from .service import ControlCenter


def ensure_token(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        token = path.read_text(encoding="utf-8").strip()
        if not token:
            raise RuntimeError(f"empty API token: {path}")
        os.chmod(path, 0o600)
        return token
    token = secrets.token_urlsafe(36)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(token + "\n")
    return token


class Handler(BaseHTTPRequestHandler):
    server_version = "UPDDControl/1"

    @property
    def control(self) -> ControlCenter:
        return self.server.control  # type: ignore[attr-defined]

    @property
    def token(self) -> str:
        return self.server.api_token  # type: ignore[attr-defined]

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def _json(self, status: int, payload: Any) -> None:
        body = json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def _authorized(self) -> bool:
        return secrets.compare_digest(self.headers.get("Authorization", ""), f"Bearer {self.token}")

    def _body(self) -> dict[str, Any]:
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError as exc:
            raise ValueError("invalid Content-Length") from exc
        if length > 65536:
            raise ValueError("request body is too large")
        if not length:
            return {}
        payload = json.loads(self.rfile.read(length).decode("utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("JSON body must be an object")
        return payload

    def _guard(self) -> bool:
        if not self._authorized():
            self._json(HTTPStatus.UNAUTHORIZED, {"error": "unauthorized"})
            return False
        return True

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/health":
            self._json(HTTPStatus.OK, {"status": "ok", "schema": "updd_control_api_v1"})
            return
        if not self._guard():
            return
        try:
            if parsed.path == "/v1/snapshot":
                self._json(HTTPStatus.OK, self.control.snapshot())
                return
            if parsed.path == "/v1/jobs":
                self._json(HTTPStatus.OK, {"jobs": self.control.registry.list_jobs()})
                return
            parts = parsed.path.strip("/").split("/")
            if len(parts) == 3 and parts[:2] == ["v1", "jobs"]:
                self._json(HTTPStatus.OK, self.control.registry.get_job(parts[2]))
                return
            if len(parts) == 4 and parts[:2] == ["v1", "jobs"] and parts[3] == "events":
                self._json(HTTPStatus.OK, {"events": self.control.registry.events(parts[2])})
                return
            if len(parts) == 4 and parts[:2] == ["v1", "jobs"] and parts[3] == "logs":
                query = parse_qs(parsed.query)
                payload = self.control.read_log(
                    parts[2],
                    offset=int(query.get("offset", [0])[0]),
                    limit=int(query.get("limit", [131072])[0]),
                    index=int(query.get("index", [0])[0]),
                )
                self._json(HTTPStatus.OK, payload)
                return
            self._json(HTTPStatus.NOT_FOUND, {"error": "not found"})
        except Exception as exc:
            self._json(HTTPStatus.BAD_REQUEST, {"error": f"{type(exc).__name__}: {exc}"})

    def do_POST(self) -> None:
        if not self._guard():
            return
        try:
            body = self._body()
            parsed = urlparse(self.path)
            if parsed.path == "/v1/discover":
                self._json(HTTPStatus.OK, {"jobs": self.control.discover()})
                return
            if parsed.path == "/v1/reconcile":
                self._json(HTTPStatus.OK, {"jobs": self.control.reconcile(actor="api")})
                return
            if parsed.path == "/v1/jobs":
                self._json(HTTPStatus.CREATED, self.control.register(body, actor="api"))
                return
            parts = parsed.path.strip("/").split("/")
            if len(parts) != 4 or parts[:2] != ["v1", "jobs"]:
                self._json(HTTPStatus.NOT_FOUND, {"error": "not found"})
                return
            job_id, action = parts[2], parts[3]
            if action == "preflight":
                payload = self.control.preflight(job_id)
            elif action == "launch":
                payload = self.control.launch(job_id)
            elif action == "pause":
                payload = self.control.pause(job_id)
            elif action == "resume":
                payload = self.control.launch(job_id, resume=True)
            elif action == "clone":
                payload = self.control.clone(job_id, body.get("overrides", {}))
            elif action == "stop-prepare":
                payload = self.control.prepare_stop(job_id, force=bool(body.get("force", False)))
            elif action == "stop-confirm":
                payload = self.control.confirm_stop(
                    job_id,
                    str(body.get("nonce", "")),
                    force=bool(body.get("force", False)),
                )
            else:
                self._json(HTTPStatus.NOT_FOUND, {"error": "unknown action"})
                return
            self._json(HTTPStatus.OK, payload)
        except Exception as exc:
            self._json(HTTPStatus.BAD_REQUEST, {"error": f"{type(exc).__name__}: {exc}"})


class Server(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, address: tuple[str, int], control: ControlCenter, token: str):
        super().__init__(address, Handler)
        self.control = control
        self.api_token = token


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--database", type=Path, default=settings.database_path())
    parser.add_argument("--token-file", type=Path, default=settings.api_token_path())
    args = parser.parse_args(argv)
    if args.host not in {"127.0.0.1", "::1", "localhost"}:
        parser.error("the control API is localhost-only")
    token = ensure_token(args.token_file.resolve())
    control = ControlCenter(registry=Registry(args.database))
    control.discover()
    control.reconcile(actor="startup")
    server = Server((args.host, args.port), control, token)
    print(f"UPDD control API listening on {args.host}:{args.port}", flush=True)
    try:
        server.serve_forever(poll_interval=0.5)
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
