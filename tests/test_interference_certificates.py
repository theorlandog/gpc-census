import copy
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

from gpc_census.certify import (
    _system,
    farkas_interference,
    verify_global_interference_certificate,
)


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


def test_fixed_support_farkas_certificate_is_exact():
    """The single-support Farkas search must rationalize to an exact vector.

    v_B on the lone determinant (0,1,2,3): the only normalized weighting is
    w = 1, whose marginals are the indicator of that determinant, so the
    restricted system is infeasible and a certificate must exist.
    """
    spectrum = [Fraction(value, 23) for value in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
    support = [(0, 1, 2, 3)]
    y = farkas_interference(4, 9, spectrum, support=support)
    assert y is not None

    dets, a, b = _system(4, 9, spectrum)
    keep = [j for j, det in enumerate(dets) if det in set(support)]
    assert sum(value * target for value, target in zip(y, b)) > 0
    assert all(
        sum(y[i] * a[i][j] for i in range(len(a))) <= 0
        for j in keep
    )


def test_farkas_search_reports_no_certificate_on_a_feasible_support():
    """POWER TEST. A support carrying a design must not yield a certificate."""
    spectrum = [Fraction(1, 2)] * 4
    assert farkas_interference(2, 4, spectrum, support=[(0, 1), (2, 3)]) is None


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
