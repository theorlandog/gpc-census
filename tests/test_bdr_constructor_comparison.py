"""Guard the BDR constructor comparison artifact and its bridge conventions.

The artifact is the scored form of `docs/prereg_bdr_constructor_comparison.md`.
These tests re-derive the load-bearing numbers rather than restating them, and
pin the two claims that would silently rot: that no external verdict is taken
on citation, and that the benchmark is labeled INCOMPLETE while the external
backend cannot be run.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gpc_census.constraints import constraints
from gpc_census.generation.bdr import (
    constraint_from_tau,
    primitive_tau,
    tau_from_constraint,
)
from gpc_census.generation.model import IntegralConstraint

ARTIFACT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "data"
    / "bdr_constructor_comparison.json"
)


@pytest.fixture(scope="module")
def report() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("bdr_constructor_comparison.json is not present in this checkout")
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def systems(report: dict) -> dict:
    return {entry["system"]: entry for entry in report["systems"]}


def test_external_backend_is_not_claimed_to_have_run(report: dict) -> None:
    """The benchmark stays INCOMPLETE until the pinned backend actually runs."""
    completeness = report["benchmark_completeness"]
    assert completeness["external_backend_executed"] is False
    assert completeness["label"] == "INCOMPLETE"
    assert completeness["what_prevented_it"]
    assert report["external_pin"]["executed_here"] is False
    assert report["external_pin"]["pinned_revision"] == (
        "7ab347d9a7837d68d92bde9fac74606913f395ca"
    )


def test_no_external_runtime_or_candidate_counts_are_recorded(report: dict) -> None:
    """No BDR stage figure may appear while the backend is unrunnable.

    The prereg forbids estimating external numbers. A future session that runs
    the backend must add those fields deliberately, and this test is where it
    will notice.
    """
    for entry in report["systems"]:
        assert set(entry.keys()) >= {"full_path", "oracle_path"}
        assert "bdr_path" not in entry
        assert "external_path" not in entry


@pytest.mark.parametrize("system", ["(3,6)", "(3,7)", "(3,8)"])
def test_every_published_row_is_locally_certified(systems: dict, system: str) -> None:
    """Prediction P1: recertification, never citation."""
    ledger = systems[system]["oracle_path"]["recertification_ledger"]
    assert ledger
    for row in ledger:
        assert row["local_status"] == "certified", row
        assert row["local_certificate_replayed_and_bound"] is True
        assert row["local_certificate_sha256"]
        assert row["matched_local_orbit"]["matched"] is True


@pytest.mark.parametrize("system", ["(3,6)", "(3,7)", "(3,8)"])
def test_converted_tau_matches_the_published_row(systems: dict, system: str) -> None:
    """The stored conversion is re-derived from the published table, not trusted."""
    entry = systems[system]
    n, d = entry["particle_number"], entry["rank"]
    published = constraints(n, d)["inequalities"]
    ledger = entry["oracle_path"]["recertification_ledger"]
    assert len(ledger) == len(published)
    for row, table_row in zip(ledger, published, strict=True):
        constraint = IntegralConstraint(
            tuple(int(value) for value in table_row["coeffs"]), int(table_row["rhs"])
        )
        assert tuple(row["converted_tau"]) == tau_from_constraint(constraint, n)


@pytest.mark.parametrize("system", ["(3,7)", "(3,8)"])
def test_oracle_path_reproduces_the_published_system(systems: dict, system: str) -> None:
    """Predictions P2 and P3, plus the hybrid gate the brief demands."""
    entry = systems[system]
    regression = entry["oracle_path"]["known_system_regression"]
    assert regression["available"] is True
    assert regression["oriented_sets_equal"] is True
    assert regression["generated_only"] == []
    assert regression["published_only"] == []
    assert regression["passed"] is True
    assert entry["oracle_path"]["counts"]["redundant_rows"] == 0
    assert entry["oracle_path"]["counts"]["unresolved_rows"] == 0
    # The hybrid must reproduce the local retained system exactly.
    assert entry["cross_path_agreement"]["retained_tau_sets_equal"] is True
    assert entry["cross_path_agreement"]["full_only"] == []
    assert entry["cross_path_agreement"]["oracle_only"] == []


def test_rank6_is_excluded_from_the_ceiling_with_a_stated_reason(
    systems: dict,
) -> None:
    """A rank the composition layer rejects is excluded visibly, not dropped."""
    entry = systems["(3,6)"]
    assert entry["composition_supported"] is False
    assert entry["ceiling"] is None
    assert "ranks 7 through 13" in entry["composition_exclusion_reason"]


@pytest.mark.parametrize(
    ("system", "expected"),
    [
        (
            "(3,7)",
            {
                "unoriented_orbits": 19,
                "admissible_hyperplanes": 5341,
                "oriented_rows_total": 10682,
                "survivors_generated": 549,
                "structural_rows": 1,
                "certified_rows": 49,
                "symbolic_zero_rejected_rows": 499,
                "unresolved_rows": 0,
                "retained_rows": 4,
            },
        ),
        (
            "(3,8)",
            {
                "unoriented_orbits": 56,
                "admissible_hyperplanes": 166420,
                "oriented_rows_total": 332840,
                "survivors_generated": 6607,
                "structural_rows": 1,
                "certified_rows": 137,
                "symbolic_zero_rejected_rows": 6469,
                "unresolved_rows": 0,
                "retained_rows": 31,
            },
        ),
    ],
)
def test_full_path_reproduces_the_published_stage_counts(
    systems: dict, system: str, expected: dict
) -> None:
    """The comparison must sit on the same stage counts the constructor doc pins."""
    counts = systems[system]["full_path"]["counts"]
    for key, value in expected.items():
        assert counts[key] == value, (system, key, counts[key], value)


@pytest.mark.parametrize("system", ["(3,7)", "(3,8)"])
def test_count_conservation_holds(systems: dict, system: str) -> None:
    """A lost orbit or a dropped survivor cannot balance these counts."""
    counts = systems[system]["full_path"]["counts"]
    assert (
        counts["survivors_generated"] == counts["survivors_predicted_in_closed_form"]
    )
    assert (
        counts["trace_rejected_never_materialised"] + counts["survivors_generated"]
        == counts["oriented_rows_total"]
    )
    assert (
        counts["structural_rows"] + counts["screened_candidate_rows"]
        == counts["survivors_generated"]
    )
    assert (
        counts["certified_rows"]
        + counts["symbolic_zero_rejected_rows"]
        + counts["trace_rejected_rows"]
        + counts["inadmissible_rows"]
        + counts["unresolved_rows"]
        == counts["screened_candidate_rows"]
    )


def test_ceiling_is_recorded_with_its_bound_semantics(systems: dict) -> None:
    """The ceiling is a bound on any external generator, and says so."""
    for system in ("(3,7)", "(3,8)"):
        ceiling = systems[system]["ceiling"]
        assert ceiling["value"] > 1.0
        assert "upper bound" in ceiling["definition"]
        # An oracle that supplies the answer still pays for every retained row.
        assert ceiling["determinant_certifications_still_required"] == (
            systems[system]["oracle_path"]["counts"]["retained_rows"]
        )


def test_measured_ceiling_does_not_grow_from_rank_7_to_rank_8(systems: dict) -> None:
    """The finding that decides the sprint, pinned so a rerun cannot lose it.

    Two points are a trend, not a theorem, and the document says so. This test
    only fixes what was measured.
    """
    assert systems["(3,8)"]["ceiling"]["value"] < systems["(3,7)"]["ceiling"]["value"]


def test_prereg_predictions_are_all_scored(report: dict) -> None:
    """Every pre-registered prediction carries a verdict, including the FAIL."""
    verdicts = {item["id"]: item["verdict"] for item in report["predictions"]}
    assert set(verdicts) == {"P1", "P2", "P3", "P4", "P5", "P6"}
    assert all(
        value in {"PASS", "FAIL", "UNINFORMATIVE", "ALREADY-FALSIFIED"}
        for value in verdicts.values()
    )
    # P5 failed on the measuring machine. A rerun that flips it must be a
    # deliberate re-score, not a silent one.
    assert verdicts["P5"] == "FAIL"


@pytest.mark.parametrize("system", ["(3,6)", "(3,7)", "(3,8)"])
def test_tau_round_trip_is_orientation_preserving(system: str) -> None:
    """The convention map in the prereg, exercised on the published rows.

    A round trip is compared after primitive positive normalization only. A
    comparison that would need a sign flip to succeed is a mismatch.
    """
    n, d = (int(part) for part in system.strip("()").split(","))
    for table_row in constraints(n, d)["inequalities"]:
        constraint = IntegralConstraint(
            tuple(int(value) for value in table_row["coeffs"]), int(table_row["rhs"])
        )
        tau = tau_from_constraint(constraint, n)
        assert tau == primitive_tau(tau)
        recovered = constraint_from_tau(tau, n)
        assert tau_from_constraint(recovered, n) == tau
        assert len(tau) == d
