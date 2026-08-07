"""Guards for the matching classification of the tangent determinant.

Two things are pinned here. First the classification itself, which has to agree
with the routine that actually decides vanishing, or it is describing a
different matrix. Second the consequence: the Hall test now runs BEFORE the
evaluation search rather than after it, and the verdict must be unchanged,
because the whole justification for the shortcut is that it reaches the same
claim by a shorter route.
"""
import collections

import pytest

from gpc_census.generation.orbit_native import enumerate_orbit_survivor_taus
from gpc_census.generation.ressayre import (
    admissible_candidate_from_tau,
    assess_candidate_by_evaluation,
    classify_tangent_determinant,
    tangent_determinant_is_identically_zero,
    verify_hall_violator,
)


def _classified(d):
    enumeration = enumerate_orbit_survivor_taus(3, d)
    tally = collections.Counter()
    records = []
    for tau in enumeration.survivor_taus:
        candidate = admissible_candidate_from_tau(3, tau)
        try:
            record = classify_tangent_determinant(3, d, candidate)
        except ValueError:          # no roots, nothing to classify
            continue
        tally[record["classification"]] += 1
        records.append((tau, candidate, record))
    return tally, records


def test_rank_7_classification_reproduces_the_published_partition():
    """474 structural and 25 algebraic zeros, 50 nonzero, as published."""
    tally, _ = _classified(7)
    assert tally["structural_zero"] == 474
    assert tally["exact_cancellation"] == 25
    assert tally["unique_matching"] == 16
    assert tally["unique_exposed_monomial"] == 24
    assert tally["nonzero_after_collision"] == 10
    assert sum(tally.values()) == 549
    zeros = tally["structural_zero"] + tally["exact_cancellation"]
    assert zeros == 499


def test_most_nonzero_determinants_need_no_summation():
    """The claim the audit exists to test, stated as a ratio not a count."""
    tally, _ = _classified(7)
    nonzero = (tally["unique_matching"] + tally["unique_exposed_monomial"]
               + tally["nonzero_after_collision"])
    cheap = tally["unique_matching"] + tally["unique_exposed_monomial"]
    assert nonzero == 50
    assert cheap / nonzero >= 0.75


def test_the_classification_agrees_with_the_routine_that_decides_vanishing():
    """A classifier that disagreed with the decider would be describing a
    different matrix, and every downstream number would inherit the error."""
    _, records = _classified(7)
    checked = 0
    for _tau, candidate, record in records:
        if record["classification"] == "structural_zero":
            assert verify_hall_violator(
                3, 7, candidate, tuple(record["hall_violator"]))
            continue
        # the expensive route, on the cases that actually need deciding
        assert (tangent_determinant_is_identically_zero(3, 7, candidate)
                == record["is_zero"])
        checked += 1
    assert checked == 75, checked


def test_matching_counts_are_consistent_with_the_monomials():
    _, records = _classified(7)
    for _tau, _candidate, record in records:
        if record["classification"] == "structural_zero":
            assert record["perfect_matchings"] == 0
            continue
        assert record["perfect_matchings"] >= record["distinct_monomials"] >= 1
        assert record["min_monomial_multiplicity"] >= 1
        if record["classification"] == "unique_matching":
            assert record["perfect_matchings"] == 1
        if record["classification"] == "nonzero_after_collision":
            assert record["min_monomial_multiplicity"] > 1


# --------------------------------------------------------------------------
# the shortcut in the evaluation search
# --------------------------------------------------------------------------

def test_proving_the_search_futile_returns_the_same_verdict_as_running_it():
    """The shortcut must be indistinguishable, not merely faster.

    Structurally zero rows previously consumed the whole trial-point budget and
    returned evaluation_unresolved. They now return the same verdict, with the
    same attempt count, from a Hall violator. If the two ever disagreed, stored
    manifests would silently change.
    """
    _, records = _classified(7)
    structural = 0
    for tau, candidate, record in records:
        if record["classification"] != "structural_zero":
            continue
        verdict = assess_candidate_by_evaluation(3, 7, candidate, attempts=32)
        assert verdict.status == "evaluation_unresolved"
        assert verdict.attempts == 32
        assert verdict.certificate is None
        # and the determinant really is zero, so no trial point could have won
        assert tangent_determinant_is_identically_zero(3, 7, candidate)
        structural += 1
        if structural >= 60:
            break
    assert structural == 60


def test_a_certified_row_still_finds_its_evaluation():
    """The direction that would break the pipeline: the shortcut must never
    intercept a row whose determinant is genuinely nonzero."""
    _, records = _classified(7)
    certified = 0
    for _tau, candidate, record in records:
        if record["is_zero"]:
            continue
        verdict = assess_candidate_by_evaluation(3, 7, candidate, attempts=32)
        assert verdict.status == "certified"
        assert verdict.certificate is not None
        assert verdict.certificate.determinant_value != 0
        certified += 1
    assert certified == 50


# --------------------------------------------------------------------------
# the polynomial-time certificate
# --------------------------------------------------------------------------

def test_a_separating_certificate_is_never_issued_for_a_zero_determinant():
    """The only direction that can break anything.

    A cheap test is allowed to say NO when the answer is yes; it is never
    allowed to say YES when the answer is no, because that would certify a
    facet that is not there. Every certificate issued at rank 7 is checked
    against the exact expansion.
    """
    from gpc_census.generation.ressayre import (
        separating_weight_certificate,
        verify_separating_weight_certificate,
    )

    _, records = _classified(7)
    issued = 0
    for _tau, candidate, record in records:
        if record["classification"] == "structural_zero":
            continue
        certificate = separating_weight_certificate(3, 7, candidate)
        if certificate is None:
            continue
        assert not record["is_zero"], "certificate issued on a vanishing determinant"
        assert verify_separating_weight_certificate(3, 7, candidate, certificate)
        issued += 1
    assert issued == 36


def test_the_certificate_reaches_most_of_what_the_expansion_can_reach():
    """36 of the 50 nonzero, against an upper bound of 40 from the audit.

    The bound is the count of determinants that HAVE a monomial with a single
    matching. The shortfall is the choice of weight vectors, not a limit of the
    method, and pinning both numbers keeps that distinction visible.
    """
    from gpc_census.generation.ressayre import separating_weight_certificate

    tally, records = _classified(7)
    exposable = tally["unique_matching"] + tally["unique_exposed_monomial"]
    issued = sum(
        1 for _t, candidate, record in records
        if record["classification"] != "structural_zero"
        and separating_weight_certificate(3, 7, candidate) is not None)
    assert exposable == 40
    assert issued == 36
    assert issued <= exposable


def test_a_tampered_certificate_is_rejected():
    """Replay must actually check, not merely parse."""
    from gpc_census.generation.ressayre import (
        separating_weight_certificate,
        verify_separating_weight_certificate,
    )

    _, records = _classified(7)
    tested = 0
    for _tau, candidate, record in records:
        if record["classification"] == "structural_zero":
            continue
        certificate = separating_weight_certificate(3, 7, candidate)
        if certificate is None:
            continue
        broken = dict(certificate)
        broken["matching_weight"] = int(certificate["matching_weight"]) + 1
        assert not verify_separating_weight_certificate(3, 7, candidate, broken)
        shuffled = dict(certificate)
        matching = list(certificate["matching"])
        if len(matching) > 1:
            matching[0], matching[1] = matching[1], matching[0]
            shuffled["matching"] = matching
            assert not verify_separating_weight_certificate(
                3, 7, candidate, shuffled)
        tested += 1
        if tested >= 8:
            break
    assert tested == 8


@pytest.mark.parametrize("d", [7])
def test_structural_share_is_the_reason_the_shortcut_pays(d):
    """If most rows were NOT structurally zero the reordering would be noise."""
    tally, _ = _classified(d)
    total = sum(tally.values())
    assert tally["structural_zero"] / total > 0.8
