"""Guards for the manifest staleness check.

`results/data/SHA256SUMS` can go stale without anyone editing it. A merge
resolves the manifest independently of the files whose digests it records, so a
force-resolved conflict leaves a digest describing content that no longer
hashes to it. That happened in merge 1756d17: both sides had changed the
PROVENANCE.md line from the same base, the conflict was resolved by taking one
side, and PROVENANCE.md itself merged to content matching neither digest. CI
caught it only after the merge had landed.

Detection is therefore the thing to keep working, and `--check` is the cheap
form of it. These tests exercise it against a synthetic stale manifest rather
than trusting the happy path.
"""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "update_data_checksums.py"

if not MODULE_PATH.exists():  # sdist may omit the script
    pytest.skip("checksum script not present", allow_module_level=True)

SPEC = importlib.util.spec_from_file_location("update_data_checksums", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def _fixture(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    (data / "nested").mkdir(parents=True)
    (data / "alpha.json").write_text("{}")
    (data / "nested" / "beta.json").write_text("[]")
    (data / "PROVENANCE.md").write_text("first\n")
    return data


def _write_manifest(data: Path) -> None:
    (data / "SHA256SUMS").write_text("\n".join(mod.manifest_lines(data)) + "\n")


def test_check_passes_on_a_fresh_manifest(tmp_path):
    data = _fixture(tmp_path)
    _write_manifest(data)
    assert mod.check(data) == 0


def test_check_catches_a_changed_file(tmp_path):
    """The ordinary case: an artifact edited without regenerating."""
    data = _fixture(tmp_path)
    _write_manifest(data)
    (data / "alpha.json").write_text('{"changed": true}')
    assert mod.check(data) == 1


def test_check_catches_the_merge_signature(tmp_path, capsys):
    """A stale PROVENANCE.md digest alone, which is what a bad merge leaves."""
    data = _fixture(tmp_path)
    _write_manifest(data)
    manifest = data / "SHA256SUMS"
    stale = manifest.read_text().replace(
        hashlib.sha256((data / "PROVENANCE.md").read_bytes()).hexdigest(),
        "0" * 64,
    )
    manifest.write_text(stale)

    assert mod.check(data) == 1
    captured = capsys.readouterr()
    assert "PROVENANCE.md" in captured.err
    assert "make checksums" in captured.err
    assert "merge" in captured.err


def test_check_catches_an_unlisted_file(tmp_path):
    data = _fixture(tmp_path)
    _write_manifest(data)
    (data / "gamma.json").write_text("0")
    assert mod.check(data) == 1


def test_check_catches_a_deleted_file(tmp_path):
    data = _fixture(tmp_path)
    _write_manifest(data)
    (data / "nested" / "beta.json").unlink()
    assert mod.check(data) == 1


def test_check_reports_a_missing_manifest(tmp_path):
    data = _fixture(tmp_path)
    assert mod.check(data) == 1


def test_the_manifest_never_lists_itself(tmp_path):
    data = _fixture(tmp_path)
    _write_manifest(data)
    assert not any(line.endswith("SHA256SUMS") for line in mod.manifest_lines(data))


def test_the_published_manifest_is_current():
    """The real dataset, which is what CI would otherwise catch after merging."""
    data = ROOT / "results" / "data"
    if not data.exists():
        pytest.skip("published dataset not present (sdist build)")
    assert mod.check(data) == 0


def test_the_remedy_targets_exist():
    makefile = (ROOT / "Makefile").read_text()
    assert "\nchecksums:" in makefile
    assert "\nchecksums-check:" in makefile
    assert "--check" in makefile
