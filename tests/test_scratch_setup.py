"""Tests for utils/scratch_setup.py — fallback semantics + env var resolution.

Covers the six decision branches of ``configure_pyscf_scratch``:
    1. UPDD_DISABLE_SCRATCH_AUTODETECT=1 → no-op
    2. PYSCF_TMPDIR pre-set → respect user override
    3. preferred dir writable → configured
    4. preferred dir not writable (read-only) → fallback
    5. preferred dir's parent missing (mount absent) → fallback
    6. UPDD_SCRATCH_DIR override picks a different candidate
"""

import os
import stat
import sys
import pytest


sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "utils"))
from scratch_setup import configure_pyscf_scratch, _probe_writable  # noqa: E402


@pytest.fixture
def clean_env():
    """Provide a fresh dict env unaffected by the real os.environ."""
    return {}


def test_disable_flag_short_circuits(clean_env):
    clean_env["UPDD_DISABLE_SCRATCH_AUTODETECT"] = "1"
    result = configure_pyscf_scratch(
        preferred="/nonexistent/should-not-be-touched",
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "disabled"
    assert result["scratch_dir"] is None
    assert "PYSCF_TMPDIR" not in clean_env
    assert "TMPDIR" not in clean_env


def test_respects_existing_pyscf_tmpdir(clean_env, tmp_path):
    clean_env["PYSCF_TMPDIR"] = "/user/already/picked"
    result = configure_pyscf_scratch(
        preferred=str(tmp_path / "scratch"),
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "respect_existing"
    assert result["scratch_dir"] == "/user/already/picked"
    assert clean_env["PYSCF_TMPDIR"] == "/user/already/picked"


def test_preferred_dir_writable_is_configured(clean_env, tmp_path):
    target = tmp_path / "pyscf_scratch"
    result = configure_pyscf_scratch(
        preferred=str(target),
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "configured"
    assert result["scratch_dir"] == str(target)
    assert clean_env["PYSCF_TMPDIR"] == str(target)
    assert clean_env["TMPDIR"] == str(target)
    assert target.is_dir()


def test_parent_missing_falls_back(clean_env):
    result = configure_pyscf_scratch(
        preferred="/definitely_not_a_mount_point_abc123/pyscf_scratch",
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "fallback"
    assert result["scratch_dir"] is None
    assert "PYSCF_TMPDIR" not in clean_env
    assert "TMPDIR" not in clean_env
    assert "not present" in result["reason"]


def test_readonly_dir_falls_back(clean_env, tmp_path):
    target = tmp_path / "readonly_scratch"
    target.mkdir()
    target.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        result = configure_pyscf_scratch(
            preferred=str(target),
            env=clean_env,
            verbose=False,
        )
        assert result["action"] == "fallback"
        assert result["scratch_dir"] is None
        assert "PYSCF_TMPDIR" not in clean_env
        assert "write probe failed" in result["reason"]
    finally:
        target.chmod(stat.S_IRWXU)


def test_env_var_override_picks_custom(clean_env, tmp_path):
    custom = tmp_path / "custom_scratch_via_env"
    clean_env["UPDD_SCRATCH_DIR"] = str(custom)
    result = configure_pyscf_scratch(env=clean_env, verbose=False)
    assert result["action"] == "configured"
    assert result["scratch_dir"] == str(custom)
    assert clean_env["PYSCF_TMPDIR"] == str(custom)


def test_empty_udpd_scratch_dir_uses_default(clean_env, tmp_path, monkeypatch):
    """Empty UPDD_SCRATCH_DIR must be treated as unset, defaulting to module default."""
    import scratch_setup
    fake_default = tmp_path / "fake_default_scratch"
    monkeypatch.setattr(scratch_setup, "DEFAULT_PREFERRED_SCRATCH", str(fake_default))
    clean_env["UPDD_SCRATCH_DIR"] = ""
    result = configure_pyscf_scratch(env=clean_env, verbose=False)
    assert result["action"] == "configured"
    assert result["scratch_dir"] == str(fake_default)
    assert clean_env["PYSCF_TMPDIR"] == str(fake_default)


def test_default_path_missing_mount_falls_back(clean_env, tmp_path, monkeypatch):
    """If the bundled default path's parent doesn't exist, fallback fires cleanly."""
    import scratch_setup
    monkeypatch.setattr(
        scratch_setup,
        "DEFAULT_PREFERRED_SCRATCH",
        "/_updd_definitely_no_mount_xyz/scratch",
    )
    result = configure_pyscf_scratch(env=clean_env, verbose=False)
    assert result["action"] == "fallback"
    assert result["scratch_dir"] is None
    assert "PYSCF_TMPDIR" not in clean_env


def test_setdefault_does_not_override_tmpdir(clean_env, tmp_path):
    """If TMPDIR is pre-set, configure_pyscf_scratch must not overwrite it."""
    preset_tmp = "/some/preset/tmpdir"
    clean_env["TMPDIR"] = preset_tmp
    target = tmp_path / "pyscf_only"
    result = configure_pyscf_scratch(
        preferred=str(target),
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "configured"
    assert clean_env["PYSCF_TMPDIR"] == str(target)
    assert clean_env["TMPDIR"] == preset_tmp


def test_probe_writable_helper_negative_cases(tmp_path):
    assert _probe_writable(str(tmp_path)) is True
    assert _probe_writable(str(tmp_path / "nonexistent_subdir")) is False
    ro = tmp_path / "ro"
    ro.mkdir()
    ro.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        assert _probe_writable(str(ro)) is False
    finally:
        ro.chmod(stat.S_IRWXU)


def test_log_format_when_verbose(clean_env, tmp_path, capsys):
    target = tmp_path / "verbose_log"
    configure_pyscf_scratch(
        preferred=str(target),
        env=clean_env,
        verbose=True,
    )
    captured = capsys.readouterr()
    assert "[SCRATCH]" in captured.out
    assert "action=configured" in captured.out
    assert str(target) in captured.out
