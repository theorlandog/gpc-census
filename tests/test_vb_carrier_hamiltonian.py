"""Exact acceptance tests for the scheduled v_B carrier Hamiltonian."""
import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "vb_carrier_hamiltonian.json"
OBSERVABILITY = DATA / "rdm_observability_census.json"


def load_artifact():
    return json.loads(ARTIFACT.read_text())


def test_exact_carrier_driver_headline():
    artifact = load_artifact()
    checks = artifact["checks"]
    assert artifact["evidence"]["status"] == "EXACT"
    assert checks["canonical_basis_dimension"] == 126
    assert checks["carrier_dimension"] == 8
    assert checks["carrier_pair_kernel_symbolic_support_entries"] == 7
    assert checks["carrier_restricted_nonzero_entries"] == 5
    assert checks["carrier_leakage_entries"] == 0
    assert checks["carrier_schrodinger_residual_entries"] == 0
    assert checks["gamma2_pair_diagonal_direction_nonzero_entries"] == 8
    assert checks["gamma2_cubic_invariant_identity"]
    assert checks["gamma2_trace_identity"]


def test_common_full_one_rdm_driver_is_exact():
    artifact = load_artifact()
    checks = artifact["checks"]
    assert checks["common_rho_leakage_entries"] == 0
    assert checks["common_rho_schrodinger_residual_entries"] == 0
    assert checks["common_rho_matrix_residual_entries"] == 0
    assert artifact["common_full_one_rdm_lift"]["beta_rate"] == (
        "nu=(28-5*t)/(sqrt(85)*(4-t^2))"
    )


def test_clock_and_pair_kernel_scope_are_pinned():
    artifact = load_artifact()
    assert artifact["conic_clock"]["t_dot"] == "-s/sqrt(85)"
    assert artifact["conic_clock"]["s_dot"] == "sqrt(85)*(t-7/17)"
    assert artifact["checks"]["physical_radicals_uniformly_positive"]
    assert artifact["checks"]["wall_regular"]
    entries = artifact["carrier_hamiltonian"]["pair_kernel_symbolic_support"]
    assert len(entries) == 7
    assert {tuple(row["left"]) for row in entries} == {
        (0, 4),
        (0, 8),
        (1, 4),
        (1, 8),
        (3, 4),
        (3, 8),
    }
    assert any("not an autonomous" in claim for claim in artifact["nonclaims"])


def test_dynamic_certificate_matches_the_census_atlas():
    artifact = load_artifact()
    observability = json.loads(OBSERVABILITY.read_text())
    row = next(
        row
        for row in observability["records"]
        if row["system"] == "(4,9)" and row["index"] == 65
    )
    assert artifact["census_record"] == {
        "system": row["system"],
        "index": row["index"],
        "classified": row["classified"],
        "observability_artifact": "results/data/rdm_observability_census.json",
    }
    assert artifact["carrier"]["support"] == row["support_dets"]
    closure = artifact["exact_td2rdm_consequence"]
    assert closure["pair_incidence_pivots"] == row["pair_incidence"]["pivot_features"]
    assert closure["gamma2_reconstruction_tree"] == row[
        "gamma2_reconstruction_tree"
    ]


def test_standalone_verifier():
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "verify_vb_carrier_hamiltonian_standalone.py"),
            str(ARTIFACT),
        ],
        cwd=ROOT,
        check=True,
    )


def test_artifact_regeneration_is_byte_identical(tmp_path):
    regenerated = tmp_path / "vb_carrier_hamiltonian.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "vb_carrier_hamiltonian.py"),
            "--out",
            str(regenerated),
        ],
        cwd=ROOT,
        check=True,
    )
    assert regenerated.read_bytes() == ARTIFACT.read_bytes()
