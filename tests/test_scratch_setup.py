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
from scratch_setup import (  # noqa: E402
    configure_pyscf_scratch,
    configure_updd_tmpdir,
    resolve_mmpbsa_workdir,
    _probe_writable,
)


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


# ---------------------------------------------------------------------------
# configure_updd_tmpdir — generic TMPDIR variant (PYSCF-independent)
# ---------------------------------------------------------------------------


def test_tmpdir_disable_flag_short_circuits(clean_env):
    clean_env["UPDD_DISABLE_SCRATCH_AUTODETECT"] = "1"
    result = configure_updd_tmpdir(
        preferred="/should-not-be-touched",
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "disabled"
    assert "TMPDIR" not in clean_env


def test_tmpdir_respects_existing(clean_env, tmp_path):
    clean_env["TMPDIR"] = "/user/already/picked"
    result = configure_updd_tmpdir(
        preferred=str(tmp_path / "ssd_tmp"),
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "respect_existing"
    assert clean_env["TMPDIR"] == "/user/already/picked"


def test_tmpdir_writable_is_configured(clean_env, tmp_path):
    target = tmp_path / "ssd_generic_tmp"
    result = configure_updd_tmpdir(
        preferred=str(target),
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "configured"
    assert result["scratch_dir"] == str(target)
    assert clean_env["TMPDIR"] == str(target)
    assert "PYSCF_TMPDIR" not in clean_env, "PYSCF_TMPDIR must NOT be touched by generic helper"
    assert target.is_dir()


def test_tmpdir_parent_missing_falls_back(clean_env):
    result = configure_updd_tmpdir(
        preferred="/no_such_mount_xyz/generic_tmp",
        env=clean_env,
        verbose=False,
    )
    assert result["action"] == "fallback"
    assert "TMPDIR" not in clean_env


def test_tmpdir_readonly_falls_back(clean_env, tmp_path):
    target = tmp_path / "readonly_tmp"
    target.mkdir()
    target.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        result = configure_updd_tmpdir(
            preferred=str(target),
            env=clean_env,
            verbose=False,
        )
        assert result["action"] == "fallback"
        assert "TMPDIR" not in clean_env
    finally:
        target.chmod(stat.S_IRWXU)


def test_tmpdir_env_var_override(clean_env, tmp_path):
    custom = tmp_path / "custom_via_UPDD_TMPDIR"
    clean_env["UPDD_TMPDIR"] = str(custom)
    result = configure_updd_tmpdir(env=clean_env, verbose=False)
    assert result["action"] == "configured"
    assert clean_env["TMPDIR"] == str(custom)


def test_tmpdir_empty_env_uses_default(clean_env, tmp_path, monkeypatch):
    import scratch_setup
    fake_default = tmp_path / "fake_default_tmp"
    monkeypatch.setattr(scratch_setup, "DEFAULT_PREFERRED_TMPDIR", str(fake_default))
    clean_env["UPDD_TMPDIR"] = ""
    result = configure_updd_tmpdir(env=clean_env, verbose=False)
    assert result["action"] == "configured"
    assert clean_env["TMPDIR"] == str(fake_default)


# ---------------------------------------------------------------------------
# resolve_mmpbsa_workdir — MMPBSA cwd SSD-preferred resolver
# ---------------------------------------------------------------------------


def test_mmpbsa_workdir_ssd_root_writable_with_subdir(clean_env, tmp_path):
    ssd_root = tmp_path / "ssd_mmpbsa"
    fallback = tmp_path / "fallback_out"
    fallback.mkdir()
    result = resolve_mmpbsa_workdir(
        fallback_dir=str(fallback),
        subdir="snap03",
        preferred_root=str(ssd_root),
        env=clean_env,
        verbose=False,
    )
    expected = ssd_root / "snap03"
    assert result == str(expected)
    assert expected.is_dir()


def test_mmpbsa_workdir_ssd_root_writable_no_subdir(clean_env, tmp_path):
    ssd_root = tmp_path / "ssd_mmpbsa_root_only"
    fallback = tmp_path / "fallback_out"
    fallback.mkdir()
    result = resolve_mmpbsa_workdir(
        fallback_dir=str(fallback),
        preferred_root=str(ssd_root),
        env=clean_env,
        verbose=False,
    )
    assert result == str(ssd_root)
    assert ssd_root.is_dir()


def test_mmpbsa_workdir_mount_absent_falls_back(clean_env, tmp_path):
    fallback = tmp_path / "fallback_out"
    fallback.mkdir()
    result = resolve_mmpbsa_workdir(
        fallback_dir=str(fallback),
        subdir="snap04",
        preferred_root="/no_such_mount_abc/scratch",
        env=clean_env,
        verbose=False,
    )
    expected = fallback / "snap04"
    assert result == str(expected)
    assert expected.is_dir()


def test_mmpbsa_workdir_readonly_falls_back(clean_env, tmp_path):
    ssd_root = tmp_path / "ssd_ro"
    ssd_root.mkdir()
    # readonly parent makes mkdir(ssd_root/subdir) impossible, but ssd_root
    # itself exists — we make ssd_root readonly to force write-probe failure
    ssd_root.chmod(stat.S_IRUSR | stat.S_IXUSR)
    fallback = tmp_path / "fallback_out"
    fallback.mkdir()
    try:
        result = resolve_mmpbsa_workdir(
            fallback_dir=str(fallback),
            subdir="snap05",
            preferred_root=str(ssd_root),
            env=clean_env,
            verbose=False,
        )
        # Expect fallback path used
        assert result == str(fallback / "snap05")
        assert (fallback / "snap05").is_dir()
    finally:
        ssd_root.chmod(stat.S_IRWXU)


def test_mmpbsa_workdir_env_var_override(clean_env, tmp_path):
    custom = tmp_path / "custom_mmpbsa_root"
    clean_env["UPDD_MMPBSA_WORKDIR_ROOT"] = str(custom)
    fallback = tmp_path / "fallback_out"
    fallback.mkdir()
    result = resolve_mmpbsa_workdir(
        fallback_dir=str(fallback),
        subdir="snap06",
        env=clean_env,
        verbose=False,
    )
    expected = custom / "snap06"
    assert result == str(expected)
    assert expected.is_dir()


def test_mmpbsa_workdir_disable_flag_forces_fallback(clean_env, tmp_path):
    ssd_root = tmp_path / "ssd_should_be_skipped"
    fallback = tmp_path / "fallback_out"
    fallback.mkdir()
    clean_env["UPDD_DISABLE_SCRATCH_AUTODETECT"] = "1"
    result = resolve_mmpbsa_workdir(
        fallback_dir=str(fallback),
        subdir="snap07",
        preferred_root=str(ssd_root),
        env=clean_env,
        verbose=False,
    )
    expected = fallback / "snap07"
    assert result == str(expected)
    assert expected.is_dir()
    assert not ssd_root.exists(), "SSD root must NOT be created when disabled"


def test_mmpbsa_workdir_empty_subdir_falls_back_to_root(clean_env, tmp_path):
    fallback = tmp_path / "fallback_only"
    result = resolve_mmpbsa_workdir(
        fallback_dir=str(fallback),
        subdir="",
        preferred_root="/no_such_mount_xyz/scratch",
        env=clean_env,
        verbose=False,
    )
    # subdir empty -> returns fallback_dir itself (mkdir'd)
    assert result == str(fallback)
    assert fallback.is_dir()


# ---------------------------------------------------------------------------
# Integration: configure_updd_tmpdir + configure_pyscf_scratch ordering
# ---------------------------------------------------------------------------


def test_orchestrator_order_separates_tmpdir_from_pyscf_tmpdir(clean_env, tmp_path):
    """Regression for UPDD.py call-order bug: TMPDIR (generic) must stay on
    /media/san/San/tmp even when configure_pyscf_scratch() also setdefault's
    TMPDIR to /media/san/San/pyscf_scratch. The orchestrator calls
    configure_updd_tmpdir() FIRST so that subsequent configure_pyscf_scratch()
    finds TMPDIR already set and only updates PYSCF_TMPDIR.
    """
    generic = tmp_path / "ssd_generic"
    pyscf = tmp_path / "ssd_pyscf"
    # Orchestrator order
    configure_updd_tmpdir(preferred=str(generic), env=clean_env, verbose=False)
    configure_pyscf_scratch(preferred=str(pyscf), env=clean_env, verbose=False)
    # TMPDIR must remain generic, PYSCF_TMPDIR points to pyscf subdir
    assert clean_env["TMPDIR"] == str(generic), \
        "TMPDIR shadowed by configure_pyscf_scratch — order bug regression"
    assert clean_env["PYSCF_TMPDIR"] == str(pyscf)


def test_reversed_order_shadows_tmpdir(clean_env, tmp_path):
    """Reversed order (pyscf first) DOES shadow TMPDIR — documents that the
    orchestrator order matters. If this test breaks, configure_pyscf_scratch
    semantics changed and the orchestrator hook needs review.
    """
    generic = tmp_path / "ssd_generic_b"
    pyscf = tmp_path / "ssd_pyscf_b"
    configure_pyscf_scratch(preferred=str(pyscf), env=clean_env, verbose=False)
    configure_updd_tmpdir(preferred=str(generic), env=clean_env, verbose=False)
    # configure_pyscf_scratch sets TMPDIR=pyscf_path too -> setdefault then
    # has no effect on TMPDIR, generic path NOT applied. This is the bug
    # the orchestrator guards against by ordering generic-first.
    assert clean_env["TMPDIR"] == str(pyscf)
    assert clean_env["PYSCF_TMPDIR"] == str(pyscf)
