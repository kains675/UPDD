#!/usr/bin/env python3
"""Dense one-screen Streamlit UI for the UPDD local control center."""

from __future__ import annotations

import html
import json
import os
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests
import streamlit as st


API_URL = os.environ.get("UPDD_CC_API_URL", "http://127.0.0.1:8765").rstrip("/")
TOKEN_FILE = Path(
    os.environ.get(
        "UPDD_CC_API_TOKEN_FILE",
        str(Path.home() / ".local/state/updd-control-center/api.token"),
    )
).expanduser()


st.set_page_config(page_title="UPDD Control Center", page_icon=":material/science:", layout="wide")
st.markdown(
    """
    <style>
    :root { --cc-green:#16855b; --cc-amber:#b56a00; --cc-red:#c33b32; --cc-cyan:#087e8b; --cc-ink:#20252a; }
    .block-container { padding-top: 3.25rem; padding-bottom: 1rem; max-width: 100%; }
    h1 { font-size: 1.45rem !important; line-height: 1.7rem !important; letter-spacing:0 !important; margin:0 !important; }
    h2, h3 { letter-spacing:0 !important; }
    [data-testid="stMetric"] { border-top: 2px solid #d7dcdf; padding-top: .35rem; }
    [data-testid="stMetricLabel"] { font-size: .72rem; }
    [data-testid="stMetricValue"] { font-size: 1.05rem; }
    [data-testid="stSidebar"] { border-right: 1px solid #d7dcdf; }
    [data-testid="stAppDeployButton"] { display:none; }
    div[data-testid="stVerticalBlockBorderWrapper"] { border-radius: 6px !important; }
    .cc-title-row { display:flex; align-items:baseline; gap:.8rem; padding-bottom:.3rem; border-bottom:1px solid #d7dcdf; }
    .cc-subtle { color:#697177; font-size:.78rem; }
    .cc-metrics { display:grid; grid-template-columns:repeat(8,minmax(0,1fr)); gap:1rem; margin:.55rem 0 1rem; }
    .cc-metric { border-top:2px solid #d7dcdf; padding-top:.35rem; min-width:0; }
    .cc-metric b { display:block; color:#697177; font-size:.72rem; font-weight:500; overflow-wrap:anywhere; }
    .cc-metric span { display:block; color:#20252a; font-size:.9rem; margin-top:.18rem; overflow-wrap:anywhere; }
    .cc-axis { display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:.35rem; margin:.25rem 0 .55rem 0; }
    .cc-state { border-left:4px solid #7d878d; background:#f4f6f7; padding:.42rem .55rem; min-height:3.25rem; }
    .cc-state b { display:block; font-size:.72rem; color:#697177; font-weight:600; }
    .cc-state span { display:block; font-size:.88rem; overflow-wrap:anywhere; }
    .cc-good { border-left-color:var(--cc-green); }
    .cc-warn { border-left-color:var(--cc-amber); }
    .cc-bad { border-left-color:var(--cc-red); }
    .cc-info { border-left-color:var(--cc-cyan); }
    .stButton button { border-radius:5px; min-height:2.25rem; }
    code { letter-spacing:0 !important; }
    @media (max-width: 800px) {
      .cc-axis { grid-template-columns:repeat(2,minmax(0,1fr)); }
      .cc-metrics { grid-template-columns:repeat(2,minmax(0,1fr)); gap:.45rem .8rem; }
      .block-container { padding-top:3rem; padding-left:.55rem; padding-right:.55rem; }
    }
    </style>
    """,
    unsafe_allow_html=True,
)


class ApiError(RuntimeError):
    pass


def token() -> str:
    try:
        return TOKEN_FILE.read_text(encoding="utf-8").strip()
    except FileNotFoundError as exc:
        raise ApiError(f"API token unavailable: {TOKEN_FILE}") from exc


def api(method: str, path: str, payload: dict[str, Any] | None = None, params: dict[str, Any] | None = None) -> Any:
    try:
        response = requests.request(
            method,
            API_URL + path,
            headers={"Authorization": f"Bearer {token()}"},
            json=payload,
            params=params,
            timeout=12,
        )
    except requests.RequestException as exc:
        raise ApiError(str(exc)) from exc
    try:
        body = response.json()
    except ValueError as exc:
        raise ApiError(f"API returned HTTP {response.status_code}") from exc
    if response.status_code >= 400:
        raise ApiError(str(body.get("error", body)))
    return body


def fmt_gb(value: Any) -> str:
    return "-" if value is None else f"{float(value):.1f} GB"


def fmt_duration(seconds: Any) -> str:
    if seconds is None:
        return "-"
    seconds = max(0, int(float(seconds)))
    if seconds < 60:
        return f"{seconds}s"
    if seconds < 3600:
        return f"{seconds // 60}m"
    return f"{seconds // 3600}h {(seconds % 3600) // 60}m"


def state_class(value: Any) -> str:
    text = str(value)
    if text in {"COMPLETE", "PRESENT", "VALID", "True"}:
        return "cc-good"
    if text in {"FAILED", "STOPPED", "ORPHANED", "MISSING", "INVALID", "False"}:
        return "cc-bad"
    if text in {"RUNNING", "QUEUED", "PREFLIGHT"}:
        return "cc-info"
    return "cc-warn"


def state_axes(job: dict[str, Any]) -> None:
    values = [
        ("Execution", job["execution_status"]),
        ("Artifacts", job["artifact_status"]),
        ("Science", job["scientific_status"]),
        ("Ranking", "UNKNOWN" if job["ranking_eligible"] is None else str(job["ranking_eligible"])),
    ]
    blocks = "".join(
        f'<div class="cc-state {state_class(value)}"><b>{html.escape(label)}</b><span>{html.escape(str(value))}</span></div>'
        for label, value in values
    )
    st.markdown(f'<div class="cc-axis">{blocks}</div>', unsafe_allow_html=True)


def host_metrics(rows: list[tuple[str, str]]) -> None:
    blocks = "".join(
        f'<div class="cc-metric"><b>{html.escape(label)}</b><span>{html.escape(value)}</span></div>'
        for label, value in rows
    )
    st.markdown(f'<div class="cc-metrics">{blocks}</div>', unsafe_allow_html=True)


def progress_values(job: dict[str, Any]) -> tuple[int, int, float | None]:
    progress = job.get("progress") or {}
    completed = int(progress.get("n_completed", progress.get("completed", 0)) or 0)
    expected = int(progress.get("n_expected", progress.get("expected", 0)) or 0)
    eta = progress.get("eta_s")
    return completed, expected, eta


def flatten(prefix: str, value: Any, rows: list[dict[str, Any]]) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            flatten(f"{prefix}.{key}" if prefix else str(key), item, rows)
    else:
        rows.append({"parameter": prefix, "value": json.dumps(value, ensure_ascii=False) if isinstance(value, (list, tuple)) else value})


def parameter_table(spec: dict[str, Any]) -> pd.DataFrame:
    sources = {
        "declared": spec.get("declared_parameters", {}),
        "default": spec.get("default_parameters", {}),
        "effective": spec.get("effective_parameters", {}),
        "source": spec.get("parameter_sources", {}),
    }
    columns: dict[str, dict[str, Any]] = {}
    for column, values in sources.items():
        rows: list[dict[str, Any]] = []
        flatten("", values, rows)
        for row in rows:
            columns.setdefault(row["parameter"], {})[column] = row["value"]
    result = [{"parameter": key, **value} for key, value in sorted(columns.items())]
    return pd.DataFrame(
        result, columns=["parameter", "declared", "default", "effective", "source"]
    ).fillna("").astype(str)


def progress_matrix(job: dict[str, Any]) -> pd.DataFrame:
    progress = job.get("progress") or {}
    cells = progress.get("cells") or []
    rows = []
    for cell in cells:
        water = cell.get("water_min_nm") or {}
        validation = cell.get("validation") or {}
        rows.append(
            {
                "arm/seed": f"{cell.get('arm', '')} {cell.get('seed', '')}".strip(),
                "status": cell.get("status", "COMPLETE"),
                "frames": cell.get("n_frames", (validation.get("dcd") or {}).get("n_frames")),
                "contact<0.26": water.get("frac_lt_0p26"),
                "mixing": (
                    "PASS" if (cell.get("mixing") or {}).get("gate_passed") is True
                    else "FAIL" if (cell.get("mixing") or {}).get("gate_passed") is False
                    else "-"
                ),
                "elapsed": fmt_duration(cell.get("elapsed_s")),
            }
        )
    if rows:
        return pd.DataFrame(rows)
    completed, expected, eta = progress_values(job)
    return pd.DataFrame([{"completed": completed, "expected": expected, "ETA": fmt_duration(eta), "gate": progress.get("scientific_gate", "-")}])


def action(method: str, path: str, payload: dict[str, Any] | None = None) -> None:
    try:
        api(method, path, payload)
        st.session_state["action_message"] = "Action accepted"
        st.session_state["action_error"] = None
    except ApiError as exc:
        st.session_state["action_error"] = str(exc)
    st.rerun()


def choose_job() -> str | None:
    try:
        snapshot = api("GET", "/v1/snapshot")
    except ApiError as exc:
        st.error(f"Control API unavailable: {exc}")
        return None
    jobs = sorted(
        snapshot.get("jobs", []),
        key=lambda job: (
            0 if job["execution_status"] in {"RUNNING", "QUEUED", "PAUSE_REQUESTED", "PAUSED"} else 1,
            0 if job["adapter"] == "dcd_sidecar" else 1,
            0 if job["scientific_status"] == "INVALID" else 1,
            job["name"],
        ),
    )
    leases = snapshot.get("leases", [])
    with st.sidebar:
        st.subheader("Campaigns")
        filter_name = st.segmented_control(
            "Status",
            ["Active", "Issues", "Complete", "All"],
            default="All",
            label_visibility="collapsed",
        )
        filtered = jobs
        if filter_name == "Active":
            filtered = [job for job in jobs if job["execution_status"] in {"RUNNING", "QUEUED", "PAUSE_REQUESTED", "PAUSED"}]
        elif filter_name == "Issues":
            filtered = [job for job in jobs if job["execution_status"] in {"FAILED", "ORPHANED"} or job["scientific_status"] == "INVALID"]
        elif filter_name == "Complete":
            filtered = [job for job in jobs if job["execution_status"] == "COMPLETE"]
        labels = {f"{job['name']} · {job['execution_status']}": job["job_id"] for job in filtered}
        if not labels:
            st.caption("No matching jobs")
            return None
        current = st.session_state.get("selected_job")
        index = list(labels.values()).index(current) if current in labels.values() else 0
        selected_label = st.radio("Job", list(labels), index=index, label_visibility="collapsed")
        selected_id = labels[selected_label]
        st.session_state["selected_job"] = selected_id
        st.divider()
        st.caption(f"GPU lease: {leases[0]['job_id'][:8] if leases else 'free'}")
        st.caption(f"Legacy observed: {len(snapshot.get('legacy_processes', []))}")
        return selected_id


@st.fragment(run_every=3)
def dashboard(selected_id: str) -> None:
    try:
        snapshot = api("GET", "/v1/snapshot")
    except ApiError as exc:
        st.error(f"Control API unavailable: {exc}")
        return
    host = snapshot["host"]
    gpu = host.get("gpu", {})
    memory = host.get("memory", {})
    disk = host.get("disk", {})
    repo = host.get("git", {})
    jobs = snapshot.get("jobs", [])
    leases = snapshot.get("leases", [])
    queue_eta = sum(
        float((job.get("progress") or {}).get("eta_s") or 0)
        for job in jobs
        if job["execution_status"] in {"RUNNING", "QUEUED", "PAUSE_REQUESTED"}
    )

    st.markdown(
        f'<div class="cc-title-row"><h1>UPDD Control Center</h1><span class="cc-subtle">{html.escape(str(repo.get("branch", "-")))} · {html.escape(str(repo.get("head", "-")))} · {"dirty" if repo.get("dirty") else "clean"}</span></div>',
        unsafe_allow_html=True,
    )
    host_metrics(
        [
            ("GPU", f"{gpu.get('utilization_pct', 0):.0f}%" if gpu.get("available") else "-"),
            ("VRAM", f"{gpu.get('memory_used_mb', 0) / 1024:.1f}/{gpu.get('memory_total_mb', 0) / 1024:.1f} GB" if gpu.get("available") else "-"),
            ("GPU power/temp", f"{gpu.get('power_w', 0):.0f} W · {gpu.get('temperature_c', 0):.0f} C" if gpu.get("available") else "-"),
            ("RAM avail", fmt_gb(memory.get("available_gb"))),
            ("Swap", f"{memory.get('swap_used_gb', 0):.1f}/{memory.get('swap_total_gb', 0):.1f} GB"),
            ("Disk free", fmt_gb(disk.get("free_gb"))),
            ("VM", str(host.get("vm", {}).get("state", "unknown"))),
            ("Queue ETA", fmt_duration(queue_eta)),
        ]
    )

    if not jobs:
        st.info("No jobs registered")
        return
    job = next((row for row in jobs if row["job_id"] == selected_id), jobs[0])
    spec = job["spec"]
    center, controls = st.columns([3.4, 1.25], gap="large")
    with center:
        st.subheader(job["name"])
        st.caption(f"{job['adapter']} · unit {job.get('unit_name') or '-'} · job {job['job_id'][:8]}")
        revalidation = (job.get("progress") or {}).get("current_runner_revalidation")
        if revalidation == "BLOCKED_CURRENT_RUNNER_DRIFT":
            st.warning(
                "Historical artifacts remain valid, but current-code revalidation is blocked by runner hash drift."
            )
        state_axes(job)
        completed, expected, eta = progress_values(job)
        pcols = st.columns([3, 1, 1])
        if expected:
            pcols[0].progress(min(1.0, completed / expected), text=f"{completed} / {expected}")
        else:
            pcols[0].caption("No scalar progress")
        pcols[1].metric("ETA", fmt_duration(eta))
        pcols[2].metric("Revision", str(job["spec_revision"]))
        overview, parameters, audit = st.tabs(["Overview", "Parameters", "Audit"])
        with overview:
            st.dataframe(progress_matrix(job), width="stretch", hide_index=True, height=255)
        with parameters:
            st.dataframe(parameter_table(spec), width="stretch", hide_index=True, height=255)
        with audit:
            events = [event for event in snapshot.get("events", []) if event.get("job_id") == job["job_id"]]
            st.dataframe(pd.DataFrame(events), width="stretch", hide_index=True, height=255)

    with controls:
        st.subheader("Controls")
        readonly = bool(job["readonly"])
        status = job["execution_status"]
        st.code(job["spec_hash"][:16], language=None)
        if st.button("Preflight", icon=":material/fact_check:", width="stretch", disabled=readonly or status not in {"DRAFT", "PREFLIGHT"}):
            action("POST", f"/v1/jobs/{job['job_id']}/preflight")
        if st.button("Launch", icon=":material/play_arrow:", width="stretch", disabled=readonly or status != "PREFLIGHT"):
            action("POST", f"/v1/jobs/{job['job_id']}/launch")
        if st.button("Pause after unit", icon=":material/pause:", width="stretch", disabled=readonly or status != "RUNNING"):
            action("POST", f"/v1/jobs/{job['job_id']}/pause")
        if st.button("Resume", icon=":material/play_circle:", width="stretch", disabled=readonly or status not in {"PAUSED", "ORPHANED"}):
            action("POST", f"/v1/jobs/{job['job_id']}/resume")
        if st.button("Stop", icon=":material/stop_circle:", width="stretch", disabled=readonly or status not in {"RUNNING", "PAUSE_REQUESTED", "QUEUED"}):
            try:
                prepared = api("POST", f"/v1/jobs/{job['job_id']}/stop-prepare")
                st.session_state["stop_nonce"] = prepared
            except ApiError as exc:
                st.session_state["action_error"] = str(exc)
        prepared = st.session_state.get("stop_nonce")
        if prepared and prepared.get("job_id") == job["job_id"]:
            if st.button("Confirm stop", icon=":material/warning:", type="primary", width="stretch"):
                action("POST", f"/v1/jobs/{job['job_id']}/stop-confirm", {"nonce": prepared["nonce"]})
        if st.session_state.get("action_error"):
            st.error(st.session_state["action_error"])
        elif st.session_state.get("action_message"):
            st.success(st.session_state["action_message"])

    st.subheader("Live log")
    log_key = f"log:{job['job_id']}"
    if st.session_state.get("log_job") != job["job_id"]:
        st.session_state["log_job"] = job["job_id"]
        st.session_state[log_key] = {"offset": 0, "text": ""}
    log_state = st.session_state.setdefault(log_key, {"offset": 0, "text": ""})
    try:
        chunk = api("GET", f"/v1/jobs/{job['job_id']}/logs", params={"offset": log_state["offset"], "limit": 65536})
        log_state["offset"] = chunk["next_offset"]
        log_state["text"] = (log_state["text"] + chunk["text"])[-120000:]
        caption = f"{chunk.get('path') or '-'} · {chunk.get('next_offset', 0)} / {chunk.get('size', 0)} bytes"
        st.caption(caption)
        st.code(log_state["text"] or "(empty)", language="text", line_numbers=False)
    except ApiError as exc:
        st.warning(str(exc))


selected_job = choose_job()
if selected_job:
    dashboard(selected_job)
