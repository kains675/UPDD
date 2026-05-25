"""PySCF scratch directory configuration with SSD fallback.

PySCF reads ``PYSCF_TMPDIR`` (and ``TMPDIR`` as fallback) at import time into
``lib.param.TMPDIR``, so ``configure_pyscf_scratch()`` must be called BEFORE
any ``import pyscf`` / ``from pyscf import ...`` statement.

Resolution order:

1. ``UPDD_DISABLE_SCRATCH_AUTODETECT=1`` → no-op (PySCF default ``/tmp``).
2. ``PYSCF_TMPDIR`` already set → respect user override, no-op.
3. ``UPDD_SCRATCH_DIR`` set (non-empty) → use that as the preferred candidate.
4. Otherwise → preferred candidate ``DEFAULT_PREFERRED_SCRATCH``.

For the preferred candidate the function attempts ``mkdir -p`` then a
write-probe (``tempfile.NamedTemporaryFile`` create+delete). If both succeed,
``PYSCF_TMPDIR`` and ``TMPDIR`` are ``setdefault``-ed to that path. If either
step fails the function returns silently — PySCF then falls back to ``/tmp``
exactly as if this module had never been called.

The function returns a dict describing what happened so the caller (or tests)
can verify the decision path.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from typing import Dict, Optional, Set


DEFAULT_PREFERRED_SCRATCH = "/media/san/San/pyscf_scratch"
DEFAULT_PREFERRED_TMPDIR = "/media/san/San/tmp"
DEFAULT_PREFERRED_MMPBSA_ROOT = "/media/san/San/mmpbsa_scratch"
DEFAULT_CHKFILE_ARCHIVE_ROOT = "/media/san/ExpDATA/UPDD_proj_Backup/chkfile_archive"


def _probe_writable(path: str) -> bool:
    """Return True iff ``path`` exists and accepts a tempfile create+delete."""
    if not os.path.isdir(path):
        return False
    try:
        with tempfile.NamedTemporaryFile(dir=path, prefix=".updd_probe_", delete=True):
            pass
    except (OSError, PermissionError):
        return False
    return True


def configure_pyscf_scratch(
    preferred: Optional[str] = None,
    *,
    env: Optional[Dict[str, str]] = None,
    verbose: bool = True,
) -> Dict[str, object]:
    """Configure PYSCF_TMPDIR / TMPDIR with SSD-preferred + ``/tmp`` fallback.

    Args:
        preferred: candidate scratch path. ``None`` resolves via the
            ``UPDD_SCRATCH_DIR`` env var then ``DEFAULT_PREFERRED_SCRATCH``.
        env: mapping treated as the process environment. Defaults to
            ``os.environ``. Mutated in place to mirror real env behavior.
        verbose: emit a one-line stdout log describing the decision.

    Returns:
        ``{"action": <decision tag>, "scratch_dir": <effective path or None>,
        "reason": <human-readable explanation>}``.

        ``action`` is one of:
        ``"disabled"`` — ``UPDD_DISABLE_SCRATCH_AUTODETECT`` set.
        ``"respect_existing"`` — ``PYSCF_TMPDIR`` was already set.
        ``"configured"`` — env vars set to ``preferred``.
        ``"fallback"`` — preferred path unusable, leaving PySCF default.
    """
    if env is None:
        env = os.environ

    if env.get("UPDD_DISABLE_SCRATCH_AUTODETECT", "").strip() == "1":
        result = {
            "action": "disabled",
            "scratch_dir": None,
            "reason": "UPDD_DISABLE_SCRATCH_AUTODETECT=1",
        }
        if verbose:
            _log_decision(result)
        return result

    existing = env.get("PYSCF_TMPDIR", "").strip()
    if existing:
        result = {
            "action": "respect_existing",
            "scratch_dir": existing,
            "reason": f"PYSCF_TMPDIR already set to {existing!r}",
        }
        if verbose:
            _log_decision(result)
        return result

    if preferred is None:
        preferred = env.get("UPDD_SCRATCH_DIR", "").strip() or DEFAULT_PREFERRED_SCRATCH

    parent = os.path.dirname(preferred.rstrip("/")) or "/"
    if not os.path.isdir(parent):
        result = {
            "action": "fallback",
            "scratch_dir": None,
            "reason": f"parent {parent!r} not present (mount missing?)",
        }
        if verbose:
            _log_decision(result)
        return result

    try:
        os.makedirs(preferred, exist_ok=True)
    except (OSError, PermissionError) as exc:
        result = {
            "action": "fallback",
            "scratch_dir": None,
            "reason": f"mkdir {preferred!r} failed: {exc}",
        }
        if verbose:
            _log_decision(result)
        return result

    if not _probe_writable(preferred):
        result = {
            "action": "fallback",
            "scratch_dir": None,
            "reason": f"write probe failed at {preferred!r}",
        }
        if verbose:
            _log_decision(result)
        return result

    env.setdefault("PYSCF_TMPDIR", preferred)
    env.setdefault("TMPDIR", preferred)
    result = {
        "action": "configured",
        "scratch_dir": preferred,
        "reason": f"PYSCF_TMPDIR + TMPDIR set to {preferred!r}",
    }
    if verbose:
        _log_decision(result)
    return result


def _log_decision(result: Dict[str, object]) -> None:
    """Emit a one-line stdout note (mirrors run_qmmm.py [DIAG] style)."""
    print(
        f"[SCRATCH] action={result['action']} "
        f"scratch_dir={result['scratch_dir']!r} reason={result['reason']}",
        file=sys.stdout,
        flush=True,
    )


def configure_updd_tmpdir(
    preferred: Optional[str] = None,
    *,
    env: Optional[Dict[str, str]] = None,
    verbose: bool = True,
) -> Dict[str, object]:
    """Set ``TMPDIR`` to a dedicated SSD path (PYSCF-independent).

    Use this in entry points that DO NOT import pyscf (MD, MM-GBSA, MMPBSA,
    AF2/MPNN wrappers, ncAA parameterization, snapshot extraction). The hook
    must run before any library reads ``tempfile.gettempdir()`` at import.

    Resolution order:

    1. ``UPDD_DISABLE_SCRATCH_AUTODETECT=1`` → no-op (system default ``/tmp``).
    2. ``TMPDIR`` already set → respect user override, no-op.
    3. ``UPDD_TMPDIR`` set (non-empty) → use that as preferred candidate.
    4. Otherwise → preferred candidate ``DEFAULT_PREFERRED_TMPDIR``.

    Probe (mkdir + write-test) on the preferred path. Failure → silent fallback
    (system ``/tmp`` via PySCF/Python tempfile default).

    Args:
        preferred: candidate TMPDIR path. ``None`` resolves via env then default.
        env: process environment to mutate. Defaults to ``os.environ``.
        verbose: emit one-line stdout log.

    Returns:
        ``{"action": <tag>, "scratch_dir": <path or None>, "reason": <text>}``.
        ``action`` ∈ {``"disabled"``, ``"respect_existing"``, ``"configured"``,
        ``"fallback"``}.
    """
    if env is None:
        env = os.environ

    if env.get("UPDD_DISABLE_SCRATCH_AUTODETECT", "").strip() == "1":
        result = {
            "action": "disabled",
            "scratch_dir": None,
            "reason": "UPDD_DISABLE_SCRATCH_AUTODETECT=1",
        }
        if verbose:
            _log_decision(result)
        return result

    existing = env.get("TMPDIR", "").strip()
    if existing:
        result = {
            "action": "respect_existing",
            "scratch_dir": existing,
            "reason": f"TMPDIR already set to {existing!r}",
        }
        if verbose:
            _log_decision(result)
        return result

    if preferred is None:
        preferred = env.get("UPDD_TMPDIR", "").strip() or DEFAULT_PREFERRED_TMPDIR

    parent = os.path.dirname(preferred.rstrip("/")) or "/"
    if not os.path.isdir(parent):
        result = {
            "action": "fallback",
            "scratch_dir": None,
            "reason": f"parent {parent!r} not present (mount missing?)",
        }
        if verbose:
            _log_decision(result)
        return result

    try:
        os.makedirs(preferred, exist_ok=True)
    except (OSError, PermissionError) as exc:
        result = {
            "action": "fallback",
            "scratch_dir": None,
            "reason": f"mkdir {preferred!r} failed: {exc}",
        }
        if verbose:
            _log_decision(result)
        return result

    if not _probe_writable(preferred):
        result = {
            "action": "fallback",
            "scratch_dir": None,
            "reason": f"write probe failed at {preferred!r}",
        }
        if verbose:
            _log_decision(result)
        return result

    env.setdefault("TMPDIR", preferred)
    result = {
        "action": "configured",
        "scratch_dir": preferred,
        "reason": f"TMPDIR set to {preferred!r}",
    }
    if verbose:
        _log_decision(result)
    return result


def snapshot_dir_files(scratch_dir: str) -> Set[str]:
    """Return absolute paths of regular files in ``scratch_dir`` (snapshot
    used as the ``pre_existing`` baseline for ``archive_chkfiles_to_hdd``)."""
    if not os.path.isdir(scratch_dir):
        return set()
    out = set()
    try:
        for entry in os.listdir(scratch_dir):
            full = os.path.join(scratch_dir, entry)
            if os.path.isfile(full):
                out.add(full)
    except (OSError, PermissionError):
        pass
    return out


def archive_chkfiles_to_hdd(
    scratch_dir: str,
    snap_basename: str,
    pre_existing: Optional[Set[str]] = None,
    *,
    archive_root: Optional[str] = None,
    env: Optional[Dict[str, str]] = None,
    verbose: bool = False,
) -> Dict[str, object]:
    """Move chkfiles created since ``pre_existing`` baseline to HDD archive.

    Implements **R-7 preservation** (never delete raw data, only move) and the
    user-requested policy: chkfile retains analytical value (HOMO-LUMO
    diagnostics, Mulliken/RESP fits, JoltQC Axis 3 comparison, reviewer
    defense). Auto-archive to slow cold storage frees the SSD scratch for the
    next SCF without losing the binary checkpoint.

    Selection rules:
    - Operates on all regular files under ``scratch_dir`` that are NOT in
      ``pre_existing`` (so concurrent SCFs from other processes are isolated).
    - Skips files whose basename starts with ``.updd_probe_`` (write-test
      sentinels left by configure_*).

    Failure semantics:
    - Disabled by ``UPDD_DISABLE_SCRATCH_AUTODETECT=1`` (one switch turns off
      every SSD-routing helper in this module).
    - Disabled by ``UPDD_CHKFILE_ARCHIVE_TO_HDD=0`` (per-helper opt-out).
    - HDD parent missing / mkdir fails / individual ``shutil.move`` fails →
      that file stays on SSD (returned in ``left_behind``). Never raises.

    Args:
        scratch_dir: SSD scratch root (typically ``lib.param.TMPDIR``).
        snap_basename: snap identifier for archive subdir (e.g. ``"snap03"``).
        pre_existing: set returned by ``snapshot_dir_files`` BEFORE the SCF.
            None means "move every candidate" (unsafe under multi-process).
        archive_root: HDD archive root. None resolves via
            ``UPDD_CHKFILE_ARCHIVE_DIR`` env then
            ``DEFAULT_CHKFILE_ARCHIVE_ROOT``.
        env: process env. Defaults to ``os.environ``.
        verbose: emit one-line log per archive op.

    Returns:
        ``{"archived": [hdd_paths], "left_behind": [ssd_paths], "reason": str}``.
    """
    if env is None:
        env = os.environ

    if env.get("UPDD_DISABLE_SCRATCH_AUTODETECT", "").strip() == "1":
        return {"archived": [], "left_behind": [], "reason": "autodetect disabled"}
    if env.get("UPDD_CHKFILE_ARCHIVE_TO_HDD", "1").strip() == "0":
        return {"archived": [], "left_behind": [], "reason": "opt-out flag"}
    if not os.path.isdir(scratch_dir):
        return {
            "archived": [],
            "left_behind": [],
            "reason": f"scratch_dir {scratch_dir!r} not present",
        }

    if archive_root is None:
        archive_root = (
            env.get("UPDD_CHKFILE_ARCHIVE_DIR", "").strip()
            or DEFAULT_CHKFILE_ARCHIVE_ROOT
        )

    parent = os.path.dirname(archive_root.rstrip("/")) or "/"
    if not os.path.isdir(parent):
        return {
            "archived": [],
            "left_behind": [],
            "reason": f"HDD parent {parent!r} not present (mount missing?)",
        }

    target_dir = os.path.join(archive_root, snap_basename)
    try:
        os.makedirs(target_dir, exist_ok=True)
    except (OSError, PermissionError) as exc:
        return {
            "archived": [],
            "left_behind": [],
            "reason": f"mkdir {target_dir!r} failed: {exc}",
        }

    pre = pre_existing if pre_existing is not None else set()
    candidates = []
    try:
        for entry in os.listdir(scratch_dir):
            if entry.startswith(".updd_probe_"):
                continue
            full = os.path.join(scratch_dir, entry)
            if not os.path.isfile(full):
                continue
            if full in pre:
                continue
            candidates.append(full)
    except (OSError, PermissionError) as exc:
        return {
            "archived": [],
            "left_behind": [],
            "reason": f"listdir {scratch_dir!r} failed: {exc}",
        }

    archived = []
    left_behind = []
    for src in candidates:
        dst = os.path.join(target_dir, os.path.basename(src))
        try:
            shutil.move(src, dst)
            archived.append(dst)
            if verbose:
                print(f"[chkfile-archive] {src} -> {dst}", flush=True)
        except (OSError, PermissionError) as exc:
            left_behind.append(src)
            if verbose:
                print(f"[chkfile-archive] skip {src} ({exc})", flush=True)

    return {
        "archived": archived,
        "left_behind": left_behind,
        "reason": "ok",
    }


def resolve_mmpbsa_workdir(
    fallback_dir: str,
    subdir: str = "",
    *,
    preferred_root: Optional[str] = None,
    env: Optional[Dict[str, str]] = None,
    verbose: bool = False,
) -> str:
    """Return preferred MMPBSA work_dir path, SSD-rooted when available.

    MMPBSA.py uses its cwd as scratch directory for ``_MMPBSA_*.mdin/.mdcrd``
    intermediates. Per-snap concurrency can push burst write to GB-scale
    aggregate. Routing the work_dir onto a dedicated SSD (``DEFAULT_PREFERRED_MMPBSA_ROOT``)
    isolates that I/O from the root NVMe partition.

    Resolution order:

    1. ``UPDD_DISABLE_SCRATCH_AUTODETECT=1`` → return ``<fallback_dir>/<subdir>``.
    2. ``UPDD_MMPBSA_WORKDIR_ROOT`` set (non-empty) → use that as root.
    3. Otherwise → ``DEFAULT_PREFERRED_MMPBSA_ROOT``.

    The chosen root is created (mkdir + write probe). Failure → fallback path.

    Args:
        fallback_dir: original output_dir to use when SSD unavailable.
            The fallback target is ``<fallback_dir>/<subdir>`` (or
            ``<fallback_dir>/tmp/<subdir>`` if subdir is empty,
            preserving existing layout).
        subdir: per-snap label appended to the resolved root (e.g. ``"snap03"``).
            Empty string returns the root itself (with mkdir).
        preferred_root: caller override (None = env/default resolution).
        env: process environment. Defaults to ``os.environ``.
        verbose: emit a one-line log on the decision.

    Returns:
        Absolute path (mkdir-applied) of the chosen work_dir. Never raises.
    """
    if env is None:
        env = os.environ

    def _fallback() -> str:
        if subdir:
            target = os.path.join(fallback_dir, subdir)
        else:
            target = fallback_dir
        os.makedirs(target, exist_ok=True)
        return target

    if env.get("UPDD_DISABLE_SCRATCH_AUTODETECT", "").strip() == "1":
        chosen = _fallback()
        if verbose:
            print(f"[SCRATCH-MMPBSA] action=disabled work_dir={chosen!r}", flush=True)
        return chosen

    if preferred_root is None:
        preferred_root = (
            env.get("UPDD_MMPBSA_WORKDIR_ROOT", "").strip()
            or DEFAULT_PREFERRED_MMPBSA_ROOT
        )

    parent = os.path.dirname(preferred_root.rstrip("/")) or "/"
    if not os.path.isdir(parent):
        chosen = _fallback()
        if verbose:
            print(
                f"[SCRATCH-MMPBSA] action=fallback "
                f"reason=parent {parent!r} absent work_dir={chosen!r}",
                flush=True,
            )
        return chosen

    try:
        os.makedirs(preferred_root, exist_ok=True)
    except (OSError, PermissionError) as exc:
        chosen = _fallback()
        if verbose:
            print(
                f"[SCRATCH-MMPBSA] action=fallback "
                f"reason=mkdir failed: {exc} work_dir={chosen!r}",
                flush=True,
            )
        return chosen

    if not _probe_writable(preferred_root):
        chosen = _fallback()
        if verbose:
            print(
                f"[SCRATCH-MMPBSA] action=fallback "
                f"reason=write probe failed work_dir={chosen!r}",
                flush=True,
            )
        return chosen

    target = os.path.join(preferred_root, subdir) if subdir else preferred_root
    os.makedirs(target, exist_ok=True)
    if verbose:
        print(
            f"[SCRATCH-MMPBSA] action=configured work_dir={target!r}",
            flush=True,
        )
    return target
