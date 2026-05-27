"""
tests/test_dispatch.py — utils/dispatch.py 의 단일-파일 테스트.

Coverage:
- GPULocation enum + route_stage (opt-in 행동 + alias 등록)
- select_openmm_platform (CUDA / CPU fallback + DeviceIndex)
- VMExecutor mocked SSH (connectivity / execute / sync)
- with_fallback (VM 실패 → host 폴백)

V100 실 하드웨어 의존 X — mock 만 사용. 실 V100 통합 테스트는 별도 numerical validation.
"""
from __future__ import annotations

import importlib.util
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, List
from unittest import mock

import pytest


# ──────────────────────────────────────────────────────────────
# Module import (pytest-session pollution 방지, tests/test_updd_cli.py 패턴)
# ──────────────────────────────────────────────────────────────
_DISPATCH_PATH = Path(__file__).resolve().parent.parent / "utils" / "dispatch.py"


def _fresh_dispatch(env: Dict[str, str] | None = None):
    """Fresh dispatch module import w/ optional env override."""
    spec = importlib.util.spec_from_file_location("dispatch_test_isolated", _DISPATCH_PATH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    if env:
        with mock.patch.dict(os.environ, env, clear=False):
            spec.loader.exec_module(mod)
    else:
        spec.loader.exec_module(mod)
    return mod


# ──────────────────────────────────────────────────────────────
# GPULocation + route_stage
# ──────────────────────────────────────────────────────────────
class TestRouteStage:
    def test_ai_inference_routes_to_host(self):
        d = _fresh_dispatch({"UPDD_VM_ENABLE": "1"})
        assert d.route_stage("rfdiffusion") == d.GPULocation.HOST_5070TI
        assert d.route_stage("proteinmpnn") == d.GPULocation.HOST_5070TI
        assert d.route_stage("af2") == d.GPULocation.HOST_5070TI
        assert d.route_stage("alphafold2") == d.GPULocation.HOST_5070TI

    def test_chemistry_compute_routes_to_vm_when_enabled(self):
        d = _fresh_dispatch({"UPDD_VM_ENABLE": "1"})
        assert d.route_stage("restrained_md") == d.GPULocation.VM_V100
        assert d.route_stage("md") == d.GPULocation.VM_V100
        assert d.route_stage("qmmm") == d.GPULocation.VM_V100
        assert d.route_stage("mmgbsa") == d.GPULocation.VM_V100
        assert d.route_stage("snapshot_extract") == d.GPULocation.VM_V100

    def test_vm_disabled_forces_host_for_chemistry(self):
        d = _fresh_dispatch({"UPDD_VM_ENABLE": "0"})
        assert d.route_stage("qmmm") == d.GPULocation.HOST_5070TI
        assert d.route_stage("restrained_md") == d.GPULocation.HOST_5070TI
        assert d.route_stage("mmgbsa") == d.GPULocation.HOST_5070TI

    def test_cpu_only_stages(self):
        d = _fresh_dispatch({"UPDD_VM_ENABLE": "1"})
        assert d.route_stage("ranking") == d.GPULocation.CPU_ONLY
        assert d.route_stage("preprocess") == d.GPULocation.CPU_ONLY
        assert d.route_stage("ncaa_mutation") == d.GPULocation.CPU_ONLY

    def test_unknown_stage_defaults_to_host(self):
        """미등록 stage 는 보수적으로 HOST_5070TI."""
        d = _fresh_dispatch({"UPDD_VM_ENABLE": "1"})
        assert d.route_stage("unknown_xyz") == d.GPULocation.HOST_5070TI

    def test_case_insensitive(self):
        d = _fresh_dispatch({"UPDD_VM_ENABLE": "1"})
        assert d.route_stage("QMMM") == d.GPULocation.VM_V100
        assert d.route_stage("  RFdiffusion  ") == d.GPULocation.HOST_5070TI


# ──────────────────────────────────────────────────────────────
# select_openmm_platform
# ──────────────────────────────────────────────────────────────
class TestSelectOpenMMPlatform:
    def test_cuda_returns_device_index_property(self):
        d = _fresh_dispatch({"UPDD_CUDA_DEVICE": "0"})
        # openmm 미설치 환경 / CUDA 부재에서도 CPU fallback 까지는 작동
        platform, props = d.select_openmm_platform(preferred="CUDA", device_id=0)
        if platform.getName() == "CUDA":
            assert props["DeviceIndex"] == "0"
            assert props["Precision"] == "mixed"
        else:
            # CPU fallback — props 는 empty
            assert props == {}

    def test_device_id_explicit_override(self):
        d = _fresh_dispatch()
        platform, props = d.select_openmm_platform(preferred="CUDA", device_id=1)
        if platform.getName() == "CUDA":
            assert props["DeviceIndex"] == "1"

    def test_cpu_returns_empty_props(self):
        d = _fresh_dispatch()
        platform, props = d.select_openmm_platform(preferred="CPU")
        assert platform.getName() == "CPU"
        assert props == {}

    def test_env_var_default_device(self):
        d = _fresh_dispatch({"UPDD_CUDA_DEVICE": "0"})
        platform, props = d.select_openmm_platform(preferred="CUDA")  # device_id=None
        if platform.getName() == "CUDA":
            assert props["DeviceIndex"] == "0"


# ──────────────────────────────────────────────────────────────
# VMExecutor — mock subprocess
# ──────────────────────────────────────────────────────────────
class TestVMExecutor:
    def test_is_connected_success(self):
        d = _fresh_dispatch()
        vm = d.VMExecutor()
        fake_result = subprocess.CompletedProcess(args=[], returncode=0, stdout="ok\n", stderr="")
        with mock.patch("subprocess.run", return_value=fake_result):
            assert vm.is_connected() is True

    def test_is_connected_failure(self):
        d = _fresh_dispatch()
        vm = d.VMExecutor()
        fake_result = subprocess.CompletedProcess(args=[], returncode=255, stdout="", stderr="Connection refused")
        with mock.patch("subprocess.run", return_value=fake_result):
            assert vm.is_connected() is False

    def test_is_connected_exception(self):
        d = _fresh_dispatch()
        vm = d.VMExecutor()
        with mock.patch("subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="ssh", timeout=5)):
            assert vm.is_connected() is False

    def test_execute_returns_success_dict(self):
        d = _fresh_dispatch()
        vm = d.VMExecutor()
        fake_result = subprocess.CompletedProcess(args=[], returncode=0, stdout="output", stderr="")
        with mock.patch("subprocess.run", return_value=fake_result):
            r = vm.execute("echo test", timeout=5)
        assert r["returncode"] == 0
        assert r["stdout"] == "output"
        assert r["retried"] == 0

    def test_execute_retries_on_failure(self):
        d = _fresh_dispatch()
        vm = d.VMExecutor()
        # 3회 retry 모두 실패 시 returncode == -1
        fake_fail = subprocess.CompletedProcess(args=[], returncode=1, stdout="", stderr="fail")
        with mock.patch("subprocess.run", return_value=fake_fail), \
             mock.patch("time.sleep"):  # speed up
            r = vm.execute("false", timeout=5, retry=3)
        assert r["returncode"] == -1
        assert r["retried"] == 3

    def test_execute_succeeds_after_retry(self):
        d = _fresh_dispatch()
        vm = d.VMExecutor()
        # 첫 retry 실패, 두 번째 성공
        fake_fail = subprocess.CompletedProcess(args=[], returncode=1, stdout="", stderr="transient")
        fake_ok = subprocess.CompletedProcess(args=[], returncode=0, stdout="recovered", stderr="")
        with mock.patch("subprocess.run", side_effect=[fake_fail, fake_ok]), \
             mock.patch("time.sleep"):
            r = vm.execute("flaky", retry=3)
        assert r["returncode"] == 0
        assert r["retried"] == 1


# ──────────────────────────────────────────────────────────────
# with_fallback chain
# ──────────────────────────────────────────────────────────────
class TestWithFallback:
    def test_non_vm_loc_runs_primary_only(self):
        d = _fresh_dispatch()
        primary = mock.MagicMock(return_value={"returncode": 0, "stdout": "host"})
        fallback = mock.MagicMock(return_value={"returncode": 0, "stdout": "fallback"})
        result = d.with_fallback(d.GPULocation.HOST_5070TI, primary, fallback)
        assert result == {"returncode": 0, "stdout": "host"}
        primary.assert_called_once()
        fallback.assert_not_called()

    def test_vm_success_skips_fallback(self):
        d = _fresh_dispatch()
        primary = mock.MagicMock(return_value={"returncode": 0, "stdout": "vm ok"})
        fallback = mock.MagicMock()
        result = d.with_fallback(d.GPULocation.VM_V100, primary, fallback)
        assert result == {"returncode": 0, "stdout": "vm ok"}
        fallback.assert_not_called()

    def test_vm_fail_triggers_host_fallback(self):
        d = _fresh_dispatch()
        primary = mock.MagicMock(return_value={"returncode": -1, "stderr": "SSH timeout"})
        fallback = mock.MagicMock(return_value={"returncode": 0, "stdout": "host rescue"})
        result = d.with_fallback(d.GPULocation.VM_V100, primary, fallback)
        assert result == {"returncode": 0, "stdout": "host rescue"}
        primary.assert_called_once()
        fallback.assert_called_once()

    def test_vm_exception_triggers_host_fallback(self):
        d = _fresh_dispatch()
        primary = mock.MagicMock(side_effect=RuntimeError("VM crashed"))
        fallback = mock.MagicMock(return_value={"returncode": 0, "stdout": "host rescue"})
        result = d.with_fallback(d.GPULocation.VM_V100, primary, fallback)
        assert result == {"returncode": 0, "stdout": "host rescue"}
        fallback.assert_called_once()


# ──────────────────────────────────────────────────────────────
# 환경 변수 기본값
# ──────────────────────────────────────────────────────────────
class TestEnvDefaults:
    def test_vm_enable_default_off(self):
        d = _fresh_dispatch({"UPDD_VM_ENABLE": "0"})
        assert d.UPDD_VM_ENABLE is False

    def test_vm_enable_explicit_on(self):
        d = _fresh_dispatch({"UPDD_VM_ENABLE": "1"})
        assert d.UPDD_VM_ENABLE is True

    def test_cuda_device_default_zero(self):
        d = _fresh_dispatch({})
        assert d.UPDD_CUDA_DEVICE == 0
        assert d.UPDD_MMGBSA_CUDA_DEVICE == 0
        assert d.UPDD_MD_CUDA_DEVICE == 0
