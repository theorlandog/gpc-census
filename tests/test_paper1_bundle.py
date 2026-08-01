"""The built Paper 1 supplement is minimal, closed, and cross-linked."""
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys
import zipfile

import pytest


REPO = pathlib.Path(__file__).resolve().parents[1]
BUILDER = REPO / "scripts/build_paper1_supplement.py"
BUNDLE_CHECKER = REPO / "scripts/verify_paper1_bundle_standalone.py"
EXPECTED_MEMBERS = {
    ".python-version",
    "LICENSE",
    "README.md",
    "SHA256SUMS",
    "SOURCE.json",
    "data/constraint_source_map.json",
    "data/constraints.json",
    "data/interference_certificates.json",
    "data/source_rows.json",
    "data/states.jsonl",
    "data/vertices/vertices_3_10.json",
    "data/vertices/vertices_3_6.json",
    "data/vertices/vertices_3_7.json",
    "data/vertices/vertices_3_8.json",
    "data/vertices/vertices_3_9.json",
    "data/vertices/vertices_4_10.json",
    "data/vertices/vertices_4_8.json",
    "data/vertices/vertices_4_9.json",
    "data/vertices/vertices_5_10.json",
    "pyproject.toml",
    "scripts/verify_constraint_source_standalone.py",
    "scripts/verify_interference_certificates_standalone.py",
    "scripts/verify_paper1_bundle_standalone.py",
    "scripts/verify_states_standalone.py",
    "scripts/verify_vertex_exhaustion_standalone.py",
    "uv.lock",
}
STATE_KEYS = {
    "classified",
    "closed_form",
    "denominator",
    "index",
    "integer_form",
    "system",
}


@pytest.fixture(scope="module")
def built_bundle(tmp_path_factory):
    directory = tmp_path_factory.mktemp("paper1-bundle")
    archive = directory / "supplement.zip"
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    environment["PYTHONNOUSERSITE"] = "1"
    proc = subprocess.run(
        [sys.executable, str(BUILDER), "--output", str(archive)],
        cwd=REPO,
        env=environment,
        capture_output=True,
        text=True,
        timeout=900,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    extracted = directory / "extracted"
    extracted.mkdir()
    with zipfile.ZipFile(archive) as bundle:
        bundle.extractall(extracted)
    return archive, extracted


def _copy_bundle(built_bundle, root: pathlib.Path) -> pathlib.Path:
    destination = root / "supplement"
    shutil.copytree(built_bundle[1], destination)
    return destination


def _refresh_checksums(root: pathlib.Path) -> None:
    lines = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "SHA256SUMS":
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(root)}")
    (root / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run(root: pathlib.Path) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    environment["PYTHONNOUSERSITE"] = "1"
    return subprocess.run(
        [sys.executable, str(BUNDLE_CHECKER), str(root)],
        cwd=REPO,
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )


def test_built_zip_has_the_exact_paper1_members(built_bundle):
    with zipfile.ZipFile(built_bundle[0]) as bundle:
        members = set(bundle.namelist())
    assert members == EXPECTED_MEMBERS


def test_bundle_crosslinks_accept_the_built_archive(built_bundle):
    proc = _run(built_bundle[1])
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "799 state-to-vertex links" in proc.stdout


def test_bundle_contains_only_exact_state_projection_keys(built_bundle):
    rows = [
        json.loads(line)
        for line in (built_bundle[1] / "data/states.jsonl").read_text().splitlines()
        if line.strip()
    ]
    assert len(rows) == 799
    assert all(set(record) == STATE_KEYS for record in rows)


def test_bundle_crosslinks_require_the_checksum_manifest(built_bundle, tmp_path):
    root = _copy_bundle(built_bundle, tmp_path)
    (root / "SHA256SUMS").unlink()
    proc = _run(root)
    assert proc.returncode != 0
    assert "mandatory SHA256SUMS" in proc.stderr


def test_bundle_crosslinks_reject_an_extra_file_even_when_checksummed(
    built_bundle, tmp_path
):
    root = _copy_bundle(built_bundle, tmp_path)
    extra = root / "data/vertices/vertices_3_11.json"
    extra.write_text("[]\n", encoding="utf-8")
    _refresh_checksums(root)
    proc = _run(root)
    assert proc.returncode != 0
    assert "archive path set mismatch" in proc.stderr


def test_bundle_crosslinks_reject_a_relabeled_certificate(built_bundle, tmp_path):
    root = _copy_bundle(built_bundle, tmp_path)
    path = root / "data/interference_certificates.json"
    artifact = json.loads(path.read_text())
    artifact["records"][0]["vertex_index"] += 1
    path.write_text(json.dumps(artifact), encoding="utf-8")
    _refresh_checksums(root)
    proc = _run(root)
    assert proc.returncode != 0
    assert "certificate" in proc.stderr.lower()


def test_bundle_crosslinks_reject_mutated_source_rows(built_bundle, tmp_path):
    root = _copy_bundle(built_bundle, tmp_path)
    path = root / "data/source_rows.json"
    artifact = json.loads(path.read_text())
    artifact["(4,9)"][0]["coeffs"][0] = 0
    path.write_text(json.dumps(artifact), encoding="utf-8")
    _refresh_checksums(root)
    proc = _run(root)
    assert proc.returncode != 0
    assert "source rows disagree" in proc.stderr


def test_bundle_crosslinks_reject_mutated_page_metadata(built_bundle, tmp_path):
    root = _copy_bundle(built_bundle, tmp_path)
    path = root / "data/constraint_source_map.json"
    artifact = json.loads(path.read_text())
    row = artifact["systems"]["(4,9)"]["rows"][0]
    row["pdf_page"] = 119
    row["thesis_page"] = 110
    path.write_text(json.dumps(artifact), encoding="utf-8")
    _refresh_checksums(root)
    proc = _run(root)
    assert proc.returncode != 0
    assert "source-map PDF page invalid" in proc.stderr


def test_bundle_crosslinks_reject_nonexact_state_fields(built_bundle, tmp_path):
    root = _copy_bundle(built_bundle, tmp_path)
    path = root / "data/states.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    rows[0]["secs"] = 0.3
    path.write_text(
        "".join(json.dumps(record, separators=(",", ":")) + "\n" for record in rows),
        encoding="utf-8",
    )
    _refresh_checksums(root)
    proc = _run(root)
    assert proc.returncode != 0
    assert "state projection contains nonexact" in proc.stderr


def test_bundle_crosslinks_validate_source_identity(built_bundle, tmp_path):
    root = _copy_bundle(built_bundle, tmp_path)
    path = root / "SOURCE.json"
    artifact = json.loads(path.read_text())
    artifact["source_commit"] = "not-a-git-object"
    path.write_text(json.dumps(artifact), encoding="utf-8")
    _refresh_checksums(root)
    proc = _run(root)
    assert proc.returncode != 0
    assert "source_commit" in proc.stderr
