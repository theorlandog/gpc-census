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


def test_the_benchmark_is_complete_and_the_backend_ran(report: dict) -> None:
    """The label is retired on evidence: the pinned backend was measured.

    It was INCOMPLETE for a stated reason that turned out to be false, and then
    for a true one, that no BDR timing had been taken. This pins the fix.
    """
    completeness = report["benchmark_completeness"]
    assert completeness["external_backend_executed"] is True
    assert completeness["label"] == "COMPLETE"
    assert completeness["external_measurements"].endswith(
        "bdr_external_benchmark.json"
    )
    assert report["external_pin"]["executed_here"] is True
    assert report["external_pin"]["pinned_revision"] == (
        "7ab347d9a7837d68d92bde9fac74606913f395ca"
    )
    assert report["external_pin"]["license"].startswith("MIT")


def test_external_numbers_are_measured_not_estimated(report: dict) -> None:
    """Every BDR figure traces to the external artifact, never to an estimate.

    The pre-registration forbids estimating external numbers. The external half
    is post hoc and must say so, and it must never be scored against P1 to P6.
    """
    head = report["head_to_head"]
    assert head["rows"]
    assert "post hoc" in head["not_preregistered"]
    for entry in report["systems"]:
        assert set(entry.keys()) >= {"full_path", "oracle_path"}
        assert "bdr_path" not in entry
        assert "external_path" not in entry
    for prediction in report["predictions"]:
        assert "bdr" not in json.dumps(prediction).lower()


def test_the_measured_arm_is_named_and_its_limits_recorded(report: dict) -> None:
    """The ratio is meaningless without the arm it was measured on."""
    head = report["head_to_head"]
    assert "probabilistic" in head["measured_arm"]
    unavailable = head["exact_arm_unavailable"]
    assert unavailable["available"] is False
    assert unavailable["largest_variable_count_that_factors"] == 48
    assert unavailable["smallest_variable_count_that_fails"] == 56
    for system, needed in unavailable["variables_needed"].items():
        assert needed > unavailable["largest_variable_count_that_factors"], system
    assert "not a like-for-like" in head["caveat"]


def test_bdr_lands_under_the_non_circular_ceiling(report: dict) -> None:
    """The bound was measured before the backend ran; it held against it."""
    systems = {entry["system"]: entry for entry in report["systems"]}
    rank8 = systems["(3,8)"]
    ceiling = rank8["birationality_prefilter_path"][
        "speedup_if_external_enumeration_were_also_free"
    ]
    measured = report["head_to_head"]["rows"]["(3,8)"]["bdr_faster_by"]
    assert 1.0 < measured < ceiling


def test_bdr_reaches_the_filters_with_far_fewer_candidate_rows(report: dict) -> None:
    """The candidate-volume question that was OPEN, now a number."""
    row = report["head_to_head"]["rows"]["(3,8)"]
    assert row["local_candidate_rows_screened"] == 6606
    assert row["bdr_candidate_rows_reaching_the_filters"] == 173
    assert row["local_candidate_rows_screened"] > (
        30 * row["bdr_candidate_rows_reaching_the_filters"]
    )


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
