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
import sys
import tempfile
from typing import Dict, Optional


DEFAULT_PREFERRED_SCRATCH = "/media/san/San/pyscf_scratch"


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
