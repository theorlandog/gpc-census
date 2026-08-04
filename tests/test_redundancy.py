from __future__ import annotations

from dataclasses import replace

import pytest

from gpc_census.constraints import constraints
from gpc_census.generation.model import IntegralConstraint
from gpc_census.generation.redundancy import (
    certify_ordered_slice_redundancy,
    ordered_pauli_inequalities,
    verify_ordered_slice_reduction_certificate,
)


def _known_rows(d: int) -> tuple[IntegralConstraint, ...]:
    return tuple(
        IntegralConstraint(tuple(row["coeffs"]), row["rhs"])
        for row in constraints(3, d)["inequalities"]
    )


@pytest.mark.parametrize("d", range(7, 14))
def test_structural_rows_and_replay_across_target_ranks(d: int) -> None:
    structural_duplicate = IntegralConstraint((1,) + (0,) * (d - 1), 1)
    trace_duplicate = IntegralConstraint((1,) * d, 3)
    cutting_facet = IntegralConstraint((4,) + (0,) * (d - 1), 3)

    certificate = certify_ordered_slice_redundancy(
        (structural_duplicate, trace_duplicate, cutting_facet), d
    )

    assert len(ordered_pauli_inequalities(d)) == d + 1
    assert certificate.redundant_indices == (0, 1)
    assert certificate.retained_indices == (2,)
    assert certificate.retained_rows == (cutting_facet,)
    assert certificate.facets[0].violation_margin > 0
    assert verify_ordered_slice_reduction_certificate(certificate)


@pytest.mark.parametrize("d", range(7, 11))
def test_every_published_rank_row_is_an_exact_facet(d: int) -> None:
    rows = _known_rows(d)
    certificate = certify_ordered_slice_redundancy(rows, d)

    assert certificate.removals == ()
    assert certificate.retained_rows == rows
    assert len(certificate.facets) == len(rows)
    assert verify_ordered_slice_reduction_certificate(certificate)


def test_tampered_farkas_and_facet_payloads_fail_replay() -> None:
    d = 7
    duplicate = IntegralConstraint((1,) + (0,) * (d - 1), 1)
    facet = IntegralConstraint((4,) + (0,) * (d - 1), 3)
    certificate = certify_ordered_slice_redundancy((duplicate, facet), d)

    removal = certificate.removals[0]
    bad_removal = replace(removal, rhs_slack=removal.rhs_slack + 1)
    assert not verify_ordered_slice_reduction_certificate(
        replace(certificate, removals=(bad_removal,))
    )

    facet_certificate = certificate.facets[0]
    bad_facet = replace(
        facet_certificate,
        violation_margin=facet_certificate.violation_margin + 1,
    )
    assert not verify_ordered_slice_reduction_certificate(replace(certificate, facets=(bad_facet,)))
