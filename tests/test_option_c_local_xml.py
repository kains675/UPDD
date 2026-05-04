"""[Option C v2] Cache indirection deprecation — regression tests.

Background
----------
Prior to Option C, three sites resolved ncAA XML / hydrogens paths solely from
``manifest["xml_path"]`` / ``manifest["hydrogens_path"]`` written by
``utils/parameterize_ncaa.py``. Those manifest values can point to a volatile
cache (``/tmp/calib_params/<run>/<RES>_gaff2.xml``) which gets wiped on reboot
or during long-running campaigns. The on-disk SSOT is the per-run
``outputs/<run>/params/<resname>_gaff2.xml`` copy that the parameterizer also
emits.

Option C v2 changes the three readers (run_mmpbsa._load_forcefield,
run_mmgbsa __main__, run_restrained_md.load_ncaa_manifest) to **prefer the
local file** and only fall back to the manifest path (with ``warnings.warn``)
when the local file is missing.

The v2 patch is **alias-aware**: the on-disk filename uses
``manifest["xml_resname"]`` (e.g. ``MLE``) which can differ from
``manifest["ncaa_code"]`` (e.g. ``NML``). Falling back to ``ncaa_code`` when
``xml_resname`` is absent preserves legacy compat.

Test plan
---------
T1  run_mmpbsa source contains the alias-aware local-prefer pattern.
T2  run_mmgbsa source contains the alias-aware local-prefer pattern.
T3  run_restrained_md source contains the alias-aware local-prefer pattern.
T4  Algorithm: legacy fallback to manifest xml_path emits a warning.
T5  Algorithm: cache-invalidation — modifying the local XML changes what gets
    resolved on the next read (no stale cache binding).
T6  Algorithm: manifest xml_path pointing at /tmp/... is BYPASSED when local
    exists (Step 5.5 regression — the prior production behavior).
T7  Algorithm (NEW): NML/MLE alias case — ncaa_code='NML',
    xml_resname='MLE', local MLE_gaff2.xml correctly resolved.

Note on test design
-------------------
The patched ``run_mmpbsa._load_forcefield`` and the inline block in
``run_mmgbsa.__main__`` carry heavy OpenMM dependencies that are not
available in the test environment. Tests T4-T7 therefore implement the
**same alias-aware resolution algorithm** in a small reference helper
``_resolve()`` and exercise it on tmp_path fixtures. T1-T3 are textual
regression guards which verify that the production source still contains
the load-bearing pattern strings — these protect against accidental revert.
"""
from __future__ import annotations

import json
import os
import warnings

import pytest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
SITE1 = os.path.join(REPO_ROOT, "scripts", "run_mmpbsa.py")
SITE2 = os.path.join(REPO_ROOT, "utils", "run_mmgbsa.py")
SITE3 = os.path.join(REPO_ROOT, "utils", "run_restrained_md.py")


# ---------------------------------------------------------------------------
# Reference resolver — mirrors the alias-aware production logic exactly.
# Used by T4-T7 to keep the tests independent of OpenMM imports while still
# pinning the algorithmic contract.
# ---------------------------------------------------------------------------
def _resolve(params_dir: str, manifest: dict):
    """Return (xml_resolved, hydrogens_resolved, warnings_raised) for one
    manifest, using the alias-aware prefer-local-with-fallback pattern from
    Option C v2.
    """
    raised = []
    ncaa_code = manifest.get("ncaa_code", "")
    xml_resname = manifest.get("xml_resname", ncaa_code)  # NML→MLE alias
    if not xml_resname:
        return None, None, raised

    local_xml = os.path.join(params_dir, f"{xml_resname}_gaff2.xml")
    if os.path.exists(local_xml):
        xml_resolved = local_xml
    else:
        xml_path = manifest.get("xml_path", "")
        if xml_path and os.path.exists(xml_path):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                warnings.warn(
                    f"[Option C legacy fallback] {ncaa_code} (xml_resname={xml_resname}): "
                    f"local {local_xml} missing; using manifest xml_path={xml_path}.",
                    stacklevel=2,
                )
                raised.extend(w)
            xml_resolved = xml_path
        else:
            xml_resolved = None

    local_h = os.path.join(params_dir, f"{xml_resname}_hydrogens.xml")
    if os.path.exists(local_h):
        h_resolved = local_h
    else:
        h_path = manifest.get("hydrogens_path", "")
        if h_path and os.path.exists(h_path):
            h_resolved = h_path
        else:
            h_resolved = None

    return xml_resolved, h_resolved, raised


# ---------------------------------------------------------------------------
# T1-T3: Production source carries the alias-aware local-prefer pattern.
# ---------------------------------------------------------------------------
ALIAS_PATTERN = 'manifest.get("xml_resname", ncaa_code'
LOCAL_FNAME_PATTERN = '_gaff2.xml'
LEGACY_WARN_PATTERN = 'Option C legacy fallback'


@pytest.mark.parametrize(
    "site_path,label",
    [
        (SITE1, "scripts/run_mmpbsa.py"),
        (SITE2, "utils/run_mmgbsa.py"),
    ],
)
def test_site_contains_alias_aware_pattern(site_path, label):
    """T1, T2: the patched ForceField loaders use xml_resname alias and
    prefer the local params/<resname>_gaff2.xml path."""
    with open(site_path, "r", encoding="utf-8") as f:
        src = f.read()
    assert ALIAS_PATTERN in src, f"{label}: missing alias-aware xml_resname fallback"
    assert LOCAL_FNAME_PATTERN in src, f"{label}: missing local_xml resolution"
    assert LEGACY_WARN_PATTERN in src, f"{label}: missing legacy fallback warning"
    assert "local_xml = os.path.join(params_dir" in src, f"{label}: local_xml join missing"


def test_site3_contains_alias_aware_pattern():
    """T3: utils/run_restrained_md.py prefers the local xml/hydrogens but
    keeps xml_resname unchanged for caller contract."""
    with open(SITE3, "r", encoding="utf-8") as f:
        src = f.read()
    # xml_resname is still read from manifest (caller contract preserved).
    assert 'manifest.get("xml_resname")' in src, "Site 3: xml_resname caller contract changed"
    # Prefer-local for xml_path
    assert 'local_xml = os.path.join(params_dir_local' in src, \
        "Site 3: local_xml join missing"
    assert LEGACY_WARN_PATTERN in src, "Site 3: missing legacy fallback warning"
    # Prefer-local for hydrogens
    assert 'local_h = os.path.join(params_dir_local' in src, \
        "Site 3: local_h join missing"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def _write_xml(path: str, content: str = "<dummy/>") -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


def _make_manifest(params_dir: str, ncaa_code: str, xml_resname: str | None,
                   xml_path: str = "", hydrogens_path: str = "") -> str:
    manifest = {"ncaa_code": ncaa_code}
    if xml_resname is not None:
        manifest["xml_resname"] = xml_resname
    if xml_path:
        manifest["xml_path"] = xml_path
    if hydrogens_path:
        manifest["hydrogens_path"] = hydrogens_path
    mf_path = os.path.join(params_dir, f"{ncaa_code}_params_manifest.json")
    with open(mf_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    return mf_path


# ---------------------------------------------------------------------------
# T4: Legacy fallback path emits warning when local missing.
# ---------------------------------------------------------------------------
def test_legacy_fallback_warns(tmp_path):
    params = tmp_path / "params"
    params.mkdir()
    cache = tmp_path / "tmp_cache"
    cache.mkdir()
    legacy_xml = cache / "MTR_gaff2.xml"
    _write_xml(str(legacy_xml), "<legacy_v1/>")
    legacy_h = cache / "MTR_hydrogens.xml"
    _write_xml(str(legacy_h), "<legacy_h_v1/>")
    _make_manifest(
        str(params), ncaa_code="MTR", xml_resname="MTR",
        xml_path=str(legacy_xml), hydrogens_path=str(legacy_h),
    )
    manifest = json.load(open(os.path.join(str(params), "MTR_params_manifest.json")))
    xml_resolved, h_resolved, warns = _resolve(str(params), manifest)
    assert xml_resolved == str(legacy_xml)
    assert h_resolved == str(legacy_h)
    assert any("Option C legacy fallback" in str(w.message) for w in warns), \
        "Legacy fallback should emit a warning"


# ---------------------------------------------------------------------------
# T5: Cache invalidation — modifying local XML changes resolution.
# ---------------------------------------------------------------------------
def test_local_read_picks_up_modifications(tmp_path):
    params = tmp_path / "params"
    params.mkdir()
    local_xml = params / "MTR_gaff2.xml"
    _write_xml(str(local_xml), "<v1/>")
    _make_manifest(str(params), ncaa_code="MTR", xml_resname="MTR")
    manifest = json.load(open(os.path.join(str(params), "MTR_params_manifest.json")))

    xml_resolved, _, warns = _resolve(str(params), manifest)
    assert xml_resolved == str(local_xml)
    assert warns == [], "Local-prefer should not warn"
    with open(xml_resolved, "r", encoding="utf-8") as f:
        assert f.read() == "<v1/>"

    # Modify in place — next read sees new content (no stale cache).
    _write_xml(str(local_xml), "<v2/>")
    xml_resolved2, _, _ = _resolve(str(params), manifest)
    assert xml_resolved2 == str(local_xml)
    with open(xml_resolved2, "r", encoding="utf-8") as f:
        assert f.read() == "<v2/>"


# ---------------------------------------------------------------------------
# T6: /tmp manifest xml_path BYPASSED when local exists (Step 5.5 regression).
# ---------------------------------------------------------------------------
def test_local_bypasses_tmp_manifest_path(tmp_path):
    params = tmp_path / "params"
    params.mkdir()
    local_xml = params / "MTR_gaff2.xml"
    _write_xml(str(local_xml), "<local_canonical/>")
    # Pretend the manifest still references /tmp/calib_params/...
    tmp_cache = tmp_path / "tmp_calib_params"
    tmp_cache.mkdir()
    stale_xml = tmp_cache / "MTR_gaff2.xml"
    _write_xml(str(stale_xml), "<stale_cache/>")
    _make_manifest(
        str(params), ncaa_code="MTR", xml_resname="MTR",
        xml_path=str(stale_xml),
    )
    manifest = json.load(open(os.path.join(str(params), "MTR_params_manifest.json")))

    xml_resolved, _, warns = _resolve(str(params), manifest)
    assert xml_resolved == str(local_xml), \
        "Local file must take precedence over /tmp manifest xml_path"
    assert warns == [], "Bypass should not emit fallback warning"

    # Even if the /tmp file is later deleted, local resolution must be unaffected.
    os.remove(str(stale_xml))
    xml_resolved2, _, _ = _resolve(str(params), manifest)
    assert xml_resolved2 == str(local_xml)


# ---------------------------------------------------------------------------
# T7 (NEW): NML/MLE alias case.
# ---------------------------------------------------------------------------
def test_nml_mle_alias_resolves_to_mle_local(tmp_path):
    """ncaa_code='NML' but xml_resname='MLE' — the on-disk filename is
    MLE_gaff2.xml. Verify the alias-aware probe finds the MLE file."""
    params = tmp_path / "params"
    params.mkdir()
    mle_xml = params / "MLE_gaff2.xml"
    _write_xml(str(mle_xml), "<mle_local/>")
    mle_h = params / "MLE_hydrogens.xml"
    _write_xml(str(mle_h), "<mle_h/>")
    # Manifest uses ncaa_code=NML (legacy registry name) but ships
    # xml_resname=MLE (the actual XML residue label).
    _make_manifest(str(params), ncaa_code="NML", xml_resname="MLE")
    manifest = json.load(open(os.path.join(str(params), "NML_params_manifest.json")))

    xml_resolved, h_resolved, warns = _resolve(str(params), manifest)
    assert xml_resolved == str(mle_xml), \
        "Alias case must resolve to MLE_gaff2.xml, not NML_gaff2.xml"
    assert h_resolved == str(mle_h)
    assert warns == [], "Alias-aware local resolution must not warn"

    # Negative control: WITHOUT the alias-aware fallback (pre-v2 buggy path)
    # we would have probed NML_gaff2.xml and missed. Verify that file
    # genuinely doesn't exist — keeps the regression real.
    assert not (params / "NML_gaff2.xml").exists()


# ---------------------------------------------------------------------------
# Extra: legacy-compat — when xml_resname absent, fall back to ncaa_code.
# ---------------------------------------------------------------------------
def test_xml_resname_absent_falls_back_to_ncaa_code(tmp_path):
    params = tmp_path / "params"
    params.mkdir()
    local_xml = params / "MTR_gaff2.xml"
    _write_xml(str(local_xml), "<legacy_no_xml_resname/>")
    # Manifest without xml_resname (legacy/older parameterize_ncaa).
    _make_manifest(str(params), ncaa_code="MTR", xml_resname=None)
    manifest = json.load(open(os.path.join(str(params), "MTR_params_manifest.json")))
    assert "xml_resname" not in manifest, "Fixture sanity"

    xml_resolved, _, warns = _resolve(str(params), manifest)
    assert xml_resolved == str(local_xml)
    assert warns == []
