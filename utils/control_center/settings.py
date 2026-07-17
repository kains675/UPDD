"""Filesystem and service defaults for the local control center."""

from __future__ import annotations

import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
QMMM_PYTHON = Path("/home/san/miniconda3/envs/qmmm/bin/python")
ATM_PYTHON = Path("/home/san/miniconda3/envs/atm/bin/python")
UI_PYTHON = Path("/home/san/miniconda3/envs/md-dashboard/bin/python")


def data_root() -> Path:
    return Path(
        os.environ.get(
            "UPDD_CC_DATA_ROOT",
            str(Path.home() / ".local/share/updd-control-center"),
        )
    ).expanduser().resolve()


def state_root() -> Path:
    return Path(
        os.environ.get(
            "UPDD_CC_STATE_ROOT",
            str(Path.home() / ".local/state/updd-control-center"),
        )
    ).expanduser().resolve()


def database_path() -> Path:
    return Path(
        os.environ.get("UPDD_CC_DATABASE", str(data_root() / "control_center.sqlite3"))
    ).expanduser().resolve()


def api_token_path() -> Path:
    return Path(
        os.environ.get("UPDD_CC_API_TOKEN_FILE", str(state_root() / "api.token"))
    ).expanduser().resolve()


def api_url() -> str:
    return os.environ.get("UPDD_CC_API_URL", "http://127.0.0.1:8765")


ALLOWED_INTERPRETERS = {QMMM_PYTHON.resolve(), ATM_PYTHON.resolve()}
ALLOWED_SCRIPTS = {
    "dcd_sidecar": {
        (REPO_ROOT / "analysis/w4a_postdensify_dcd_20260716/postdensify_dcd.py").resolve(),
    },
    "scaffold_md": {
        (REPO_ROOT / "analysis/1ycr_wt_seed_expansion_20260716/expand_1ycr_wt_scaffolds.py").resolve(),
    },
    "trackb_pool": {
        (REPO_ROOT / "scripts/trackb_inplace_rbfe_production.py").resolve(),
    },
    "tracka": {(REPO_ROOT / "UPDD.py").resolve()},
    "analysis": {
        (REPO_ROOT / "analysis/w4a_postdensify_dcd_20260716/postdensify_gate.py").resolve(),
    },
    "self_test": {(REPO_ROOT / "utils/control_center/self_test_job.py").resolve()},
}

ALLOWED_ENVIRONMENT_KEYS = {
    "CUDA_VISIBLE_DEVICES",
    "UPDD_MD_CUDA_DEVICE",
    "UPDD_MD_DCD_INTERVAL",
    "UPDD_MMGBSA_PLATFORM",
}
