"""SQLite registry with immutable specs, orthogonal state, events, and leases."""

from __future__ import annotations

import json
import secrets
import sqlite3
import time
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

from .models import ArtifactStatus, ExecutionStatus, JobSpec, ScientificStatus


SCHEMA_VERSION = 1


class RegistryError(RuntimeError):
    pass


class LeaseConflict(RegistryError):
    pass


class Registry:
    def __init__(self, path: Path):
        self.path = Path(path).expanduser().resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._initialize()

    def _connect(self) -> sqlite3.Connection:
        con = sqlite3.connect(self.path, timeout=30.0)
        con.row_factory = sqlite3.Row
        con.execute("PRAGMA foreign_keys=ON")
        con.execute("PRAGMA journal_mode=WAL")
        con.execute("PRAGMA synchronous=FULL")
        return con

    @contextmanager
    def transaction(self, immediate: bool = False) -> Iterator[sqlite3.Connection]:
        con = self._connect()
        try:
            con.execute("BEGIN IMMEDIATE" if immediate else "BEGIN")
            yield con
            con.commit()
        except Exception:
            con.rollback()
            raise
        finally:
            con.close()

    def _initialize(self) -> None:
        with self.transaction(immediate=True) as con:
            con.executescript(
                """
                CREATE TABLE IF NOT EXISTS metadata (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS jobs (
                    job_id TEXT PRIMARY KEY,
                    source_key TEXT UNIQUE,
                    parent_job_id TEXT REFERENCES jobs(job_id),
                    campaign_id TEXT NOT NULL,
                    name TEXT NOT NULL,
                    adapter TEXT NOT NULL,
                    spec_revision INTEGER NOT NULL,
                    spec_json TEXT NOT NULL,
                    spec_hash TEXT NOT NULL,
                    input_digest TEXT NOT NULL,
                    execution_status TEXT NOT NULL,
                    artifact_status TEXT NOT NULL,
                    scientific_status TEXT NOT NULL,
                    ranking_eligible INTEGER,
                    readonly INTEGER NOT NULL,
                    unit_name TEXT,
                    control_token TEXT,
                    wrapper_record TEXT,
                    primary_log TEXT,
                    progress_json TEXT NOT NULL DEFAULT '{}',
                    last_error TEXT,
                    created_at REAL NOT NULL,
                    updated_at REAL NOT NULL
                );
                CREATE TRIGGER IF NOT EXISTS jobs_spec_no_update
                BEFORE UPDATE OF
                    job_id,source_key,parent_job_id,campaign_id,name,adapter,
                    spec_revision,spec_json,spec_hash,input_digest,readonly,created_at
                ON jobs BEGIN SELECT RAISE(ABORT, 'job spec is immutable'); END;
                CREATE TABLE IF NOT EXISTS events (
                    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    job_id TEXT REFERENCES jobs(job_id),
                    timestamp REAL NOT NULL,
                    actor TEXT NOT NULL,
                    action TEXT NOT NULL,
                    result TEXT NOT NULL,
                    reason TEXT,
                    spec_hash TEXT,
                    details_json TEXT NOT NULL DEFAULT '{}'
                );
                CREATE TRIGGER IF NOT EXISTS events_no_update
                BEFORE UPDATE ON events BEGIN SELECT RAISE(ABORT, 'events are append-only'); END;
                CREATE TRIGGER IF NOT EXISTS events_no_delete
                BEFORE DELETE ON events BEGIN SELECT RAISE(ABORT, 'events are append-only'); END;
                CREATE TABLE IF NOT EXISTS artifacts (
                    artifact_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    job_id TEXT NOT NULL REFERENCES jobs(job_id),
                    kind TEXT NOT NULL,
                    path TEXT NOT NULL,
                    digest TEXT,
                    size INTEGER,
                    status TEXT NOT NULL,
                    observed_at REAL NOT NULL,
                    UNIQUE(job_id, kind, path)
                );
                CREATE TABLE IF NOT EXISTS dependencies (
                    job_id TEXT NOT NULL REFERENCES jobs(job_id),
                    dependency_job_id TEXT NOT NULL REFERENCES jobs(job_id),
                    PRIMARY KEY(job_id, dependency_job_id)
                );
                CREATE TABLE IF NOT EXISTS leases (
                    resource TEXT PRIMARY KEY,
                    job_id TEXT NOT NULL REFERENCES jobs(job_id),
                    acquired_at REAL NOT NULL
                );
                CREATE TABLE IF NOT EXISTS action_nonces (
                    nonce TEXT PRIMARY KEY,
                    job_id TEXT NOT NULL REFERENCES jobs(job_id),
                    action TEXT NOT NULL,
                    expires_at REAL NOT NULL,
                    consumed_at REAL
                );
                """
            )
            con.execute(
                "INSERT INTO metadata(key,value) VALUES('schema_version',?) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                (str(SCHEMA_VERSION),),
            )

    @staticmethod
    def _row(row: sqlite3.Row | None) -> dict[str, Any] | None:
        if row is None:
            return None
        payload = dict(row)
        for key in ("spec_json", "progress_json"):
            payload[key.removesuffix("_json")] = json.loads(payload.pop(key))
        payload["readonly"] = bool(payload["readonly"])
        ranking = payload["ranking_eligible"]
        payload["ranking_eligible"] = None if ranking is None else bool(ranking)
        return payload

    def create_job(
        self,
        spec: JobSpec,
        *,
        job_id: str | None = None,
        revision: int = 1,
        parent_job_id: str | None = None,
        execution_status: ExecutionStatus = ExecutionStatus.DRAFT,
        artifact_status: ArtifactStatus = ArtifactStatus.MISSING,
        scientific_status: ScientificStatus = ScientificStatus.UNKNOWN,
        ranking_eligible: bool | None = None,
        progress: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        job_id = job_id or str(uuid.uuid4())
        now = time.time()
        with self.transaction(immediate=True) as con:
            if spec.source_key:
                existing = con.execute(
                    "SELECT * FROM jobs WHERE source_key=?", (spec.source_key,)
                ).fetchone()
                if existing:
                    payload = self._row(existing)
                    if (
                        payload["spec_hash"] != spec.spec_hash
                        and not (payload["readonly"] and spec.readonly)
                    ):
                        raise RegistryError(
                            "source_key is already bound to a different immutable spec"
                        )
                    return payload  # type: ignore[return-value]
            con.execute(
                """INSERT INTO jobs(
                    job_id,source_key,parent_job_id,campaign_id,name,adapter,
                    spec_revision,spec_json,spec_hash,input_digest,
                    execution_status,artifact_status,scientific_status,
                    ranking_eligible,readonly,progress_json,created_at,updated_at
                ) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                (
                    job_id,
                    spec.source_key,
                    parent_job_id,
                    spec.campaign_id,
                    spec.name,
                    spec.adapter,
                    revision,
                    json.dumps(spec.to_dict(), sort_keys=True),
                    spec.spec_hash,
                    spec.input_digest,
                    execution_status.value,
                    artifact_status.value,
                    scientific_status.value,
                    None if ranking_eligible is None else int(ranking_eligible),
                    int(spec.readonly),
                    json.dumps(progress or {}, sort_keys=True),
                    now,
                    now,
                ),
            )
            for dependency in spec.dependencies:
                con.execute(
                    "INSERT OR IGNORE INTO dependencies(job_id,dependency_job_id) VALUES(?,?)",
                    (job_id, dependency),
                )
        self.append_event(job_id, "system", "CREATE", "OK", "immutable spec registered")
        return self.get_job(job_id)

    def get_job(self, job_id: str) -> dict[str, Any]:
        with self._connect() as con:
            row = con.execute("SELECT * FROM jobs WHERE job_id=?", (job_id,)).fetchone()
        payload = self._row(row)
        if payload is None:
            raise KeyError(job_id)
        return payload

    def list_jobs(self) -> list[dict[str, Any]]:
        with self._connect() as con:
            rows = con.execute("SELECT * FROM jobs ORDER BY created_at DESC").fetchall()
        return [self._row(row) for row in rows]  # type: ignore[list-item]

    def update_runtime(
        self,
        job_id: str,
        *,
        execution_status: ExecutionStatus | None = None,
        unit_name: str | None = None,
        control_token: str | None = None,
        wrapper_record: str | None = None,
        primary_log: str | None = None,
        progress: dict[str, Any] | None = None,
        last_error: str | None = None,
    ) -> dict[str, Any]:
        fields: list[str] = ["updated_at=?"]
        values: list[Any] = [time.time()]
        supplied = {
            "execution_status": None if execution_status is None else execution_status.value,
            "unit_name": unit_name,
            "control_token": control_token,
            "wrapper_record": wrapper_record,
            "primary_log": primary_log,
            "progress_json": None if progress is None else json.dumps(progress, sort_keys=True),
            "last_error": last_error,
        }
        for key, value in supplied.items():
            if value is not None:
                fields.append(f"{key}=?")
                values.append(value)
        values.append(job_id)
        with self.transaction(immediate=True) as con:
            result = con.execute(f"UPDATE jobs SET {','.join(fields)} WHERE job_id=?", values)
            if result.rowcount != 1:
                raise KeyError(job_id)
        return self.get_job(job_id)

    def update_evidence(
        self,
        job_id: str,
        *,
        artifact_status: ArtifactStatus | None = None,
        scientific_status: ScientificStatus | None = None,
        ranking_eligible: bool | None | object = ...,
        progress: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        fields = ["updated_at=?"]
        values: list[Any] = [time.time()]
        if artifact_status is not None:
            fields.append("artifact_status=?")
            values.append(artifact_status.value)
        if scientific_status is not None:
            fields.append("scientific_status=?")
            values.append(scientific_status.value)
        if ranking_eligible is not ...:
            fields.append("ranking_eligible=?")
            values.append(None if ranking_eligible is None else int(bool(ranking_eligible)))
        if progress is not None:
            fields.append("progress_json=?")
            values.append(json.dumps(progress, sort_keys=True))
        values.append(job_id)
        with self.transaction(immediate=True) as con:
            if con.execute(f"UPDATE jobs SET {','.join(fields)} WHERE job_id=?", values).rowcount != 1:
                raise KeyError(job_id)
        return self.get_job(job_id)

    def append_event(
        self,
        job_id: str | None,
        actor: str,
        action: str,
        result: str,
        reason: str | None = None,
        details: dict[str, Any] | None = None,
    ) -> int:
        spec_hash = None
        if job_id:
            try:
                spec_hash = self.get_job(job_id)["spec_hash"]
            except KeyError:
                spec_hash = None
        with self.transaction(immediate=True) as con:
            cursor = con.execute(
                """INSERT INTO events(
                    job_id,timestamp,actor,action,result,reason,spec_hash,details_json
                ) VALUES(?,?,?,?,?,?,?,?)""",
                (
                    job_id,
                    time.time(),
                    actor,
                    action,
                    result,
                    reason,
                    spec_hash,
                    json.dumps(details or {}, sort_keys=True),
                ),
            )
            return int(cursor.lastrowid)

    def events(self, job_id: str | None = None, limit: int = 200) -> list[dict[str, Any]]:
        limit = max(1, min(int(limit), 1000))
        with self._connect() as con:
            if job_id:
                rows = con.execute(
                    "SELECT * FROM events WHERE job_id=? ORDER BY event_id DESC LIMIT ?",
                    (job_id, limit),
                ).fetchall()
            else:
                rows = con.execute(
                    "SELECT * FROM events ORDER BY event_id DESC LIMIT ?", (limit,)
                ).fetchall()
        events = []
        for row in reversed(rows):
            item = dict(row)
            item["details"] = json.loads(item.pop("details_json"))
            events.append(item)
        return events

    def acquire_lease(self, resource: str, job_id: str) -> None:
        with self.transaction(immediate=True) as con:
            row = con.execute("SELECT job_id FROM leases WHERE resource=?", (resource,)).fetchone()
            if row and row["job_id"] != job_id:
                raise LeaseConflict(f"{resource} is held by {row['job_id']}")
            con.execute(
                "INSERT INTO leases(resource,job_id,acquired_at) VALUES(?,?,?) "
                "ON CONFLICT(resource) DO UPDATE SET job_id=excluded.job_id, acquired_at=excluded.acquired_at",
                (resource, job_id, time.time()),
            )

    def release_lease(self, resource: str, job_id: str) -> bool:
        with self.transaction(immediate=True) as con:
            result = con.execute(
                "DELETE FROM leases WHERE resource=? AND job_id=?", (resource, job_id)
            )
            return result.rowcount == 1

    def leases(self) -> list[dict[str, Any]]:
        with self._connect() as con:
            return [dict(row) for row in con.execute("SELECT * FROM leases ORDER BY resource")]

    def create_nonce(self, job_id: str, action: str, ttl_s: int = 60) -> str:
        nonce = secrets.token_urlsafe(24)
        with self.transaction(immediate=True) as con:
            con.execute(
                "INSERT INTO action_nonces(nonce,job_id,action,expires_at) VALUES(?,?,?,?)",
                (nonce, job_id, action, time.time() + ttl_s),
            )
        return nonce

    def consume_nonce(self, nonce: str, job_id: str, action: str) -> bool:
        now = time.time()
        with self.transaction(immediate=True) as con:
            row = con.execute(
                "SELECT * FROM action_nonces WHERE nonce=?", (nonce,)
            ).fetchone()
            if (
                row is None
                or row["job_id"] != job_id
                or row["action"] != action
                or row["consumed_at"] is not None
                or float(row["expires_at"]) < now
            ):
                return False
            con.execute("UPDATE action_nonces SET consumed_at=? WHERE nonce=?", (now, nonce))
            return True
