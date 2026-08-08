"""Guard the facet-index stability artifact.

Two claims here would rot silently if only restated in prose: that the
coefficient-one criterion is refuted in the converse direction, and that
trailing-zero padding is not validity preserving. Both are re-derived from the
published tables rather than read back from the artifact wherever that is
cheap.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

from gpc_census.constraints import constraints
from gpc_census.dataset import vertices
from gpc_census.generation.bdr import tau_from_constraint
from gpc_census.generation.cyclic_schubert import certify_cyclic_schubert
from gpc_census.generation.model import IntegralConstraint

ARTIFACT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "data"
    / "facet_index_stability.json"
)


@pytest.fixture(scope="module")
def report() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("facet_index_stability.json is not present in this checkout")
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def coefficients(report: dict) -> dict:
    return {entry["system"]: entry for entry in report["part_b_coefficient_census"]}


@pytest.fixture(scope="module")
def shapes(report: dict) -> dict:
    return {entry["system"]: entry for entry in report["part_a_shape_census"]}


def test_landing_agrees_with_the_scratch_tables(report: dict) -> None:
    """Part A is landed, so the landed values must equal what scratch found."""
    agreement = report["part_a_landing"]["scratch_agreement"]
    assert report["part_a_landing"]["verdict"] == "LANDED"
    for group, entries in agreement.items():
        assert entries, group
        for key, record in entries.items():
            assert record["agrees"] is True, (group, key, record)


@pytest.mark.parametrize("system", ["(3,7)", "(3,8)"])
def test_forward_direction_holds_every_facet_has_coefficient_one(
    coefficients: dict, system: str
) -> None:
    """The half that survives, and the half the prefilter's soundness rests on."""
    entry = coefficients[system]
    assert entry["criterion"]["forward_holds_every_facet_has_one"] is True
    assert entry["facets"]["coefficient_one"] == entry["facet_rows"]
    assert entry["facets"]["distinct_coefficients"] == [1]


@pytest.mark.parametrize(
    ("system", "expected_counterexamples"), [("(3,7)", 23), ("(3,8)", 36)]
)
def test_converse_is_refuted_with_counterexamples(
    coefficients: dict, system: str, expected_counterexamples: int
) -> None:
    """Coefficient one does NOT imply facet. This is the load-bearing finding."""
    entry = coefficients[system]
    assert entry["criterion"]["converse_holds_no_removed_row_has_one"] is False
    assert (
        entry["criterion"]["removed_rows_with_coefficient_one"]
        == expected_counterexamples
    )
    assert entry["criterion"]["separates_exactly"] is False
    assert entry["criterion"]["counterexamples"]


def test_a_stored_counterexample_replays_exactly(coefficients: dict) -> None:
    """Recompute one counterexample rather than trusting the stored value."""
    entry = coefficients["(3,7)"]
    example = entry["criterion"]["counterexamples"][0]
    constraint = example["constraint"]
    verdict = certify_cyclic_schubert(
        tuple(constraint["coeffs"]),
        int(constraint["rhs"]),
        particle_number=3,
        include_modular=False,
    )
    assert verdict.coefficient == 1
    assert verdict.coefficient == example["coefficient"]
    # And it is genuinely not a published facet.
    published = {
        tau_from_constraint(
            IntegralConstraint(tuple(row["coeffs"]), int(row["rhs"])), 3
        )
        for row in constraints(3, 7)["inequalities"]
    }
    assert tuple(example["tau"]) not in published


@pytest.mark.parametrize("system", ["(3,7)", "(3,8)"])
def test_no_row_is_unresolved_or_ill_formed(coefficients: dict, system: str) -> None:
    """A refutation built on unresolved rows would not be a refutation."""
    entry = coefficients[system]
    for group in ("facets", "removed"):
        assert entry[group]["unresolved"] == 0
        assert entry[group]["ill_formed"] == 0
        assert entry[group]["coefficient_zero"] == 0
    assert entry["facet_rows"] + entry["removed_rows"] == entry["certified_rows"]


@pytest.mark.parametrize("system", ["(3,7)", "(3,8)"])
def test_prefilter_is_sound_and_priced(coefficients: dict, system: str) -> None:
    """The prefilter must not change the answer, only the cost."""
    prefilter = coefficients[system]["coefficient_one_prefilter"]
    assert prefilter["added_after_scoring"] is True
    assert prefilter["retained_sets_identical"] is True
    assert prefilter["reduction_input_rows_after"] < (
        prefilter["reduction_input_rows_before"]
    )
    # It keeps every facet, so the filtered input is at least the facet count.
    assert prefilter["reduction_input_rows_after"] >= (
        coefficients[system]["facet_rows"]
    )
    assert prefilter["input_reduction_fraction"] > 0.4


def test_padding_is_not_validity_preserving(report: dict) -> None:
    """Seven of the eight decidable maps carry a violated padded row."""
    decidable = [
        entry
        for entry in report["part_a_padding_census"]
        if entry.get("validity_available")
    ]
    assert len(decidable) == 8
    failing = [entry for entry in decidable if entry["violated_at_higher_rank"] > 0]
    assert len(failing) == 7
    for entry in failing:
        assert entry["padding_is_a_valid_map"] is False
        assert entry["violation_examples"]


def test_a_padded_facet_is_never_merely_redundant(report: dict) -> None:
    """The sharper structural fact: violated or supporting, never interior."""
    for entry in report["part_a_padding_census"]:
        if entry.get("validity_available"):
            assert entry["valid_but_strictly_interior"] == 0, entry["map"]


def test_borland_dennis_padded_is_violated_at_rank_seven() -> None:
    """Re-derive the headline padding failure from the published data."""
    row = constraints(3, 6)["inequalities"][0]
    padded = tuple(int(v) for v in row["coeffs"]) + (0,)
    rhs = Fraction(int(row["rhs"]))
    worst = max(
        sum(c * Fraction(str(x)) for c, x in zip(padded, vertex["spectrum"], strict=True))
        for vertex in vertices(3, 7)
    )
    assert worst > rhs


def test_shape_compression_holds_across_every_published_system(
    report: dict,
) -> None:
    """Facet counts explode, shape counts do not."""
    for entry in report["part_a_shape_census"]:
        assert entry["distinct_sorted_tau_shapes"] <= entry["facet_count"]
        if entry["facet_count"] > 4:
            assert entry["shape_to_facet_ratio"] <= 0.30, entry["system"]


def test_entry_bound_does_not_transfer_across_particle_number(
    shapes: dict,
) -> None:
    """No bound uniform in N is claimed, and the data is why."""
    assert shapes["(3,10)"]["max_abs_tau_entry"] == 7
    assert shapes["(4,10)"]["max_abs_tau_entry"] > 7
    assert shapes["(5,10)"]["max_abs_tau_entry"] > shapes["(4,10)"]["max_abs_tau_entry"]


def test_shape_census_recomputes_from_the_published_tables(shapes: dict) -> None:
    """The census must not drift from the constraint tables it summarises."""
    for system, entry in shapes.items():
        n, d = entry["particle_number"], entry["rank"]
        rows = constraints(n, d)["inequalities"]
        assert entry["facet_count"] == len(rows)
        taus = [
            tau_from_constraint(
                IntegralConstraint(tuple(r["coeffs"]), int(r["rhs"])), n
            )
            for r in rows
        ]
        assert entry["max_abs_tau_entry"] == max(
            max(abs(v) for v in tau) for tau in taus
        ), system
        assert entry["distinct_sorted_tau_shapes"] == len(
            {tuple(sorted(tau)) for tau in taus}
        ), system


def test_all_predictions_scored(report: dict) -> None:
    verdicts = {item["id"]: item["verdict"] for item in report["predictions"]}
    assert set(verdicts) == {"P1", "P2", "P3", "P4", "P5", "P6"}
    assert all(value == "PASS" for value in verdicts.values())
