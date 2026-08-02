"""Exact acceptance tests for the pure-state RDM observability atlas."""
import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"

ARTIFACT = DATA / "rdm_observability_census.json"
STATES = DATA / "states.jsonl"


def load_artifact():
    return json.loads(ARTIFACT.read_text())


def test_exact_headline_counts():
    artifact = load_artifact()
    summary = artifact["summary"]
    assert artifact["conventions"]["pair"] == "P=(p_1<p_2) and A_P=a_{p_2}a_{p_1}"
    assert summary["records"] == 799
    assert summary["pair_incidence_full_column_rank"] == 799
    assert summary["gamma2_reconstructible"] == 771
    assert summary["gamma2_unresolved_by_this_certificate"] == 28
    assert summary["by_class"] == {
        "DESIGN-INT": {
            "records": 631,
            "gamma2_reconstructible": 603,
            "gamma2_unresolved": 28,
            "positive_singleton_kernel": 0,
        },
        "DESIGN-REAL": {
            "records": 12,
            "gamma2_reconstructible": 12,
            "gamma2_unresolved": 0,
            "positive_singleton_kernel": 0,
        },
        "INTERFERENCE": {
            "records": 156,
            "gamma2_reconstructible": 156,
            "gamma2_unresolved": 0,
            "positive_singleton_kernel": 63,
        },
    }


def test_every_reconstruction_certificate_has_a_tree():
    for row in load_artifact()["records"]:
        if row["pure_full_support_reconstruction_from_gamma2"]:
            assert len(row["gamma2_reconstruction_tree"]) == row["support_size"] - 1
            assert row["gamma2_reconstruction_tree"] == row[
                "first_connected_reconstruction_tree"
            ]
            assert row["gamma2_channels"]["unique_graph_connected"]
            assert row["pair_incidence"]["pivot_minor_determinant"] != 0
        else:
            assert row["gamma2_reconstruction_tree"] == []
        assert len(row["first_connected_reconstruction_tree"]) == row["support_size"] - 1


def test_all_interference_states_are_reconstructible():
    rows = [
        row
        for row in load_artifact()["records"]
        if row["classified"] == "INTERFERENCE"
    ]
    assert len(rows) == 156
    assert all(row["pure_full_support_reconstruction_from_gamma2"] for row in rows)


def test_positive_singleton_kernels_are_all_interference():
    artifact = load_artifact()
    summary = artifact["summary"]
    assert summary["positive_singleton_kernel"] == 63
    assert summary["positive_singleton_kernel_by_dimension"] == {
        "1": 47,
        "2": 15,
        "5": 1,
    }
    rows = [
        row
        for row in artifact["records"]
        if row["singleton_incidence"]["kernel_dimension"] > 0
    ]
    assert all(row["classified"] == "INTERFERENCE" for row in rows)


def test_unique_channel_hierarchy_reaches_every_support():
    assert load_artifact()["summary"]["first_connected_unique_channel_order"] == {
        "2": 771,
        "3": 22,
        "4": 5,
        "5": 1,
    }


def test_vb_has_an_eight_vertex_gamma2_tree():
    row = next(
        row
        for row in load_artifact()["records"]
        if row["system"] == "(4,9)" and row["index"] == 65
    )
    assert row["classified"] == "INTERFERENCE"
    assert row["support_size"] == 8
    assert row["singleton_incidence"]["kernel_dimension"] == 1
    assert len(row["gamma2_reconstruction_tree"]) == 7


def test_standalone_verifier():
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "verify_rdm_observability_standalone.py"),
            str(STATES),
            str(ARTIFACT),
        ],
        cwd=ROOT,
        check=True,
    )


def test_artifact_regeneration_is_byte_identical(tmp_path):
    regenerated = tmp_path / "rdm_observability_census.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "rdm_observability_census.py"),
            str(STATES),
            "--out",
            str(regenerated),
        ],
        cwd=ROOT,
        check=True,
    )
    assert regenerated.read_bytes() == ARTIFACT.read_bytes()
