"""[v3 13] tests/conftest.py — utils 디렉토리를 sys.path 에 등록하여 테스트가
프로젝트 루트(`pytest`) / tests/ 디렉토리(`pytest tests/`) 어디서 실행되어도
import 가 동일하게 동작하도록 한다.
"""
import os
import sys

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
_UTILS_DIR = os.path.join(_REPO_ROOT, "utils")
if _UTILS_DIR not in sys.path:
    sys.path.insert(0, _UTILS_DIR)
