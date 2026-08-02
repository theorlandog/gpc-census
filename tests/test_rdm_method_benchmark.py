import json
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "rdm_method_benchmark.json"


def summary():
    return {row["method"]: row for row in json.loads(ARTIFACT.read_text())["summary"]}


def test_benchmark_accuracy_and_compression():
    rows = summary()
    assert rows["sparse2rdm"]["max_reference_error"] < 1e-10
    assert rows["dense2rdm"]["max_reference_error"] < 1e-10
    assert rows["sparse2rdm"]["real_ode_variables"] < rows["fci"]["real_ode_variables"]
    assert rows["sparse2rdm"]["real_ode_variables"] * 20 < rows["dense2rdm"]["real_ode_variables"]


def test_claim_boundary_is_honest():
    rows = summary()
    assert rows["carrier"]["median_elapsed_seconds"] < rows["sparse2rdm"]["median_elapsed_seconds"]
    assert rows["fci"]["median_elapsed_seconds"] < rows["sparse2rdm"]["median_elapsed_seconds"]


def test_variable_counts_follow_from_the_representations():
    """These are the machine-independent part of the benchmark.

    Runtimes and RSS depend on the host; the variable counts, the shared step
    counts and the reference errors do not.
    """
    rows = summary()
    assert rows["fci"]["real_ode_variables"] == 2 * 210
    assert rows["carrier"]["real_ode_variables"] == 2 * 29
    assert rows["dense2rdm"]["real_ode_variables"] == 2 * 45 * 45
    assert rows["sparse2rdm"]["real_ode_variables"] == 29 + 2 * 69


def test_all_methods_integrate_the_same_problem():
    rows = summary()
    assert {row["median_nfev"] for row in rows.values()} == {314}
    assert {row["median_accepted_steps"] for row in rows.values()} == {26}


def test_reported_ratios_match_the_reported_summary():
    data = json.loads(ARTIFACT.read_text())
    rows = summary()
    ratios = data["ratios"]
    assert ratios["dense_vs_sparse_variable_count"] == (
        rows["dense2rdm"]["real_ode_variables"] / rows["sparse2rdm"]["real_ode_variables"]
    )
    assert ratios["fci_vs_sparse_variable_count"] == (
        rows["fci"]["real_ode_variables"] / rows["sparse2rdm"]["real_ode_variables"]
    )
    assert ratios["dense_vs_sparse_runtime"] > 1.0
