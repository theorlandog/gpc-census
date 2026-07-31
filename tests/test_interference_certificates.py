import copy
import json
import subprocess
import sys
from pathlib import Path

from gpc_census.certify import verify_global_interference_certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "data" / "interference_certificates.json"
STANDALONE = ROOT / "scripts" / "verify_interference_certificates_standalone.py"


def load_artifact():
    return json.loads(ARTIFACT.read_text())


def test_exact_certificate_scope_and_summary():
    artifact = load_artifact()
    records = artifact["records"]
    assert artifact["schema_version"] == 1
    assert artifact["record_count"] == 12
    assert [(tuple(row["system"]), row["vertex_index"]) for row in records] == [
        ((3, 8), 10),
        ((3, 8), 14),
        ((3, 8), 15),
        ((3, 8), 19),
        ((3, 8), 20),
        ((3, 8), 24),
        ((3, 8), 28),
        ((3, 8), 29),
        ((3, 8), 30),
        ((3, 8), 31),
        ((3, 8), 32),
        ((4, 9), 65),
    ]
    assert artifact["summary"] == {
        "expanded_clauses": 5141,
        "farkas_vectors": 192,
        "maximum_proof_depth": 8,
        "proof_nodes": 226,
    }


def test_library_verifies_every_global_certificate():
    artifact = load_artifact()
    assert all(
        verify_global_interference_certificate(record)
        for record in artifact["records"]
    )


def test_standalone_verifier_accepts_artifact():
    result = subprocess.run(
        [sys.executable, str(STANDALONE), str(ARTIFACT)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "all 12 global interference certificates PASS" in result.stdout


def test_tampered_farkas_vector_is_rejected():
    record = copy.deepcopy(load_artifact()["records"][0])
    record["farkas_vectors"][0] = [0] * len(record["farkas_vectors"][0])
    assert not verify_global_interference_certificate(record)


def test_incomplete_branch_is_rejected():
    record = copy.deepcopy(load_artifact()["records"][-1])
    branch = next(node for node in record["proof"]["nodes"] if node["kind"] == "branch")
    branch["children"].pop()
    assert not verify_global_interference_certificate(record)


def test_standalone_verifier_rejects_tampered_summary(tmp_path):
    artifact = load_artifact()
    artifact["summary"]["proof_nodes"] -= 1
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(artifact))
    result = subprocess.run(
        [sys.executable, str(STANDALONE), str(tampered)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "FAIL: artifact summary" in result.stdout
