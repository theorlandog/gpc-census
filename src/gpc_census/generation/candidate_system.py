"""Fail-closed composition of screened fixed-``N=3`` candidate systems.

The candidate screen and the ordered-slice reducer prove different facts.
The screen binds an exact Ressayre determinant witness to each accepted,
oriented candidate.  The reducer proves implications and relative facetness
for the polyhedron defined by the rows it receives.  This module composes the
two only in that order.  An unresolved or uncertified row is never allowed to
cut the polyhedron used for redundancy elimination.

A finalized result proves the following conditional statement: among the
supplied, determinant-certified candidates, the retained rows are an exact
irredundant facet description on the trace-three ordered Pauli slice.  It does
not prove that the upstream candidate search was exhaustive and therefore
does not, by itself, prove completeness of a generalized Pauli system.

Ranks 7 through 10 have an additional fail-closed regression gate.  Their
retained primitive oriented homogeneous rows must equal the published lookup
set exactly.  Ranks 11 through 13 have no lookup oracle and are labeled only
relative to the supplied candidate set.

For ``s`` screened rows, ``c`` certified rows, and ``f`` retained rows,
composition bookkeeping costs ``O(sd)`` plus replay of the supplied
determinant certificates.  Redundancy uses the ``O(c)`` linear programs and
``O(c^2 d)`` exact replay bound documented in :mod:`redundancy`.  The
potentially expensive cyclic-Schubert evaluator runs only ``f`` times rather
than ``c`` times.
"""

from __future__ import annotations

from dataclasses import dataclass

from gpc_census.constraints import constraints

from . import candidate_pipeline
from .bdr import (
    Tau,
    constraint_from_tau,
    is_pauli_lower_tau,
    primitive_tau,
    tau_from_constraint,
)
from .cyclic_schubert import (
    DEFAULT_MODULUS,
    CyclicSchubertVerdict,
    certify_cyclic_schubert,
    verify_exact_cyclic_schubert_certificate,
    verify_modular_cyclic_schubert_certificate,
)
from .full_dimension import verify_full_dimension_certificate
from .model import IntegralConstraint, RessayreCertificate
from .redundancy import (
    FacetWitnessCertificate,
    OrderedSliceReductionCertificate,
    certify_ordered_slice_redundancy,
    verify_ordered_slice_reduction_certificate,
)
from .ressayre import verify_tau_ressayre_certificate


FINALIZED_KNOWN_STATUS = "finalized_known_system_regression"
FINALIZED_RELATIVE_STATUS = "finalized_relative_to_supplied_candidates"
BLOCKED_UNRESOLVED_STATUS = "blocked_unresolved_candidates"
BLOCKED_KNOWN_MISMATCH_STATUS = "blocked_known_system_mismatch"


class CandidateSystemInvariantError(ValueError):
    """Raised when a screening object violates a composition invariant."""


@dataclass(frozen=True)
class KnownSystemRegression:
    """Exact comparison with the published rank-7 through rank-10 systems."""

    available: bool
    expected_count: int | None
    retained_count: int
    oriented_sets_equal: bool | None
    generated_only: tuple[Tau, ...]
    published_only: tuple[Tau, ...]
    lookup_source: str | None

    @property
    def passed(self) -> bool:
        """Return true when this optional gate does not block finalization."""
        return not self.available or self.oriented_sets_equal is True

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic exact regression payload."""
        return {
            "available": self.available,
            "expected_count": self.expected_count,
            "retained_count": self.retained_count,
            "oriented_sets_equal": self.oriented_sets_equal,
            "generated_only": [list(tau) for tau in self.generated_only],
            "published_only": [list(tau) for tau in self.published_only],
            "lookup_source": self.lookup_source,
            "passed": self.passed,
        }


@dataclass(frozen=True)
class FinalCandidateRow:
    """One retained row with validity, facet, and cross-validation bindings."""

    screening_index: int
    reduction_index: int
    tau: Tau
    constraint: IntegralConstraint
    ressayre_resolution_method: str
    ressayre_certificate: RessayreCertificate
    facet_certificate: FacetWitnessCertificate
    cyclic_schubert: CyclicSchubertVerdict

    def as_dict(self) -> dict[str, object]:
        """Return the replayable certificate bundle for this retained row."""
        return {
            "screening_index": self.screening_index,
            "reduction_index": self.reduction_index,
            "tau": list(self.tau),
            "constraint": self.constraint.as_dict(),
            "ressayre_resolution_method": self.ressayre_resolution_method,
            "ressayre_certificate": self.ressayre_certificate.as_dict(),
            "ordered_slice_facet_certificate": self.facet_certificate.as_dict(),
            "cyclic_schubert_cross_check": self.cyclic_schubert.as_dict(),
        }


@dataclass(frozen=True)
class CandidateSystemCounts:
    """Conservation counts across screening, reduction, and finalization."""

    normalized_tau_rows: int
    structural_rows: int
    screened_candidate_rows: int
    certified_rows: int
    exact_rejected_rows: int
    unresolved_rows: int
    redundancy_input_rows: int
    redundant_rows: int
    retained_rows: int
    retained_cyclic_cross_checks: int

    def as_dict(self) -> dict[str, int]:
        """Return stable count names for serialized artifacts."""
        return {
            "normalized_tau_rows": self.normalized_tau_rows,
            "structural_rows": self.structural_rows,
            "screened_candidate_rows": self.screened_candidate_rows,
            "certified_rows": self.certified_rows,
            "exact_rejected_rows": self.exact_rejected_rows,
            "unresolved_rows": self.unresolved_rows,
            "redundancy_input_rows": self.redundancy_input_rows,
            "redundant_rows": self.redundant_rows,
            "retained_rows": self.retained_rows,
            "retained_cyclic_cross_checks": self.retained_cyclic_cross_checks,
        }


@dataclass(frozen=True)
class CandidateSystemResult:
    """Diagnostic or finalized result of the fail-closed composition."""

    d: int
    status: str
    screening: candidate_pipeline.CandidateScreeningResult
    reduction: OrderedSliceReductionCertificate | None
    retained_rows: tuple[FinalCandidateRow, ...]
    known_system_regression: KnownSystemRegression | None
    counts: CandidateSystemCounts
    include_modular_cross_checks: bool

    @property
    def finalized(self) -> bool:
        """Whether a final relative facet system is present."""
        return self.status in {FINALIZED_KNOWN_STATUS, FINALIZED_RELATIVE_STATUS}

    @property
    def exit_status(self) -> int:
        """Return a stable process status for the composition outcome."""
        if self.finalized:
            return 0
        if self.status == BLOCKED_UNRESOLVED_STATUS:
            return 2
        if self.status == BLOCKED_KNOWN_MISMATCH_STATUS:
            return 3
        raise AssertionError(f"unknown candidate-system status: {self.status}")

    def as_dict(self) -> dict[str, object]:
        """Return a complete diagnostic artifact, finalized or blocked."""
        if self.finalized:
            final_system: dict[str, object] | None = {
                "particle_number": 3,
                "rank": self.d,
                "rows": [row.as_dict() for row in self.retained_rows],
            }
        else:
            final_system = None
        if self.status == BLOCKED_UNRESOLVED_STATUS:
            reduction_status = "not_run_due_to_unresolved_candidates"
            facet_status = "not_established"
        else:
            reduction_status = "proved_and_replayable_for_certified_rows"
            facet_status = "proved_relative_to_supplied_certified_candidates"
        if self.known_system_regression is None:
            known_status = "not_run_due_to_unresolved_candidates"
        elif not self.known_system_regression.available:
            known_status = "not_available_beyond_published_frontier"
        elif self.known_system_regression.passed:
            known_status = "exact_oriented_tau_set_equality"
        else:
            known_status = "exact_oriented_tau_set_mismatch"
        return {
            "particle_number": 3,
            "rank": self.d,
            "status": self.status,
            "finalized": self.finalized,
            "exit_status": self.exit_status,
            "include_modular_retained_cross_checks": (self.include_modular_cross_checks),
            "counts": self.counts.as_dict(),
            "screening": self.screening.as_dict(),
            "ordered_slice_reduction_certificate": (
                None if self.reduction is None else self.reduction.as_dict()
            ),
            "retained_candidate_rows": [row.as_dict() for row in self.retained_rows],
            "known_system_regression": (
                None
                if self.known_system_regression is None
                else self.known_system_regression.as_dict()
            ),
            "final_system": final_system,
            "proof_status": {
                "candidate_normalization_and_deduplication": "exact",
                "ressayre_certified_rows": "replayed_and_bound_to_oriented_tau",
                "unresolved_candidates": (
                    "none" if self.counts.unresolved_rows == 0 else "blocks_finalization"
                ),
                "exact_redundancy_elimination": reduction_status,
                "ordered_slice_relative_facetness": facet_status,
                "retained_cyclic_schubert": (
                    "cross_validation_performed_with_exact_status_preserved"
                    if self.finalized
                    else "not_run_before_finalization_gates"
                ),
                "known_system_regression": known_status,
                "candidate_exhaustiveness": "external_BDR_search_obligation_not_established",
                "system_completeness": "not_established_without_candidate_exhaustiveness",
                "polyhedral_certificate_scope": (
                    "relative_to_the_supplied_determinant_certified_candidates"
                ),
            },
        }


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CandidateSystemInvariantError(message)


def _certified_ressayre_certificate(
    row: candidate_pipeline.ScreenedCandidateRow,
) -> RessayreCertificate:
    """Return the exact certificate selected by the row's resolution method."""
    _require(row.status == "certified", "only certified rows may enter reduction")
    _require(
        row.certificate_replayed_and_bound is True,
        "certified row lacks its replay-and-binding flag",
    )
    assessment = row.ressayre_assessment
    _require(assessment is not None, "certified row lacks a Ressayre assessment")
    if row.resolution_method == "deterministic_evaluation":
        _require(
            assessment.status == "certified" and assessment.certificate is not None,
            "evaluation-certified row lacks its determinant certificate",
        )
        _require(
            row.symbolic_fallback is None,
            "evaluation-certified row unexpectedly has a symbolic fallback",
        )
        certificate = assessment.certificate
    elif row.resolution_method == "exact_symbolic_determinant":
        symbolic = row.symbolic_fallback
        _require(
            assessment.status == "evaluation_unresolved" and assessment.certificate is None,
            "symbolically certified row has inconsistent evaluation status",
        )
        _require(
            symbolic is not None
            and symbolic.status == "certified"
            and symbolic.certificate is not None
            and symbolic.certificate_replayed_and_bound is True,
            "symbolically certified row lacks its exact bound certificate",
        )
        certificate = symbolic.certificate
    else:
        raise CandidateSystemInvariantError(
            f"unknown certification method: {row.resolution_method}"
        )
    _require(
        verify_tau_ressayre_certificate(3, row.tau, certificate),
        "Ressayre certificate failed exact replay or oriented-tau binding",
    )
    return certificate


def _validate_noncertified_row(row: candidate_pipeline.ScreenedCandidateRow) -> None:
    """Check that rejected and unresolved rows carry no positive certificate."""
    _require(
        row.certificate_replayed_and_bound is not True,
        "noncertified row carries a positive replay-and-binding flag",
    )
    if row.status == "trace_rejected":
        _require(
            row.ressayre_assessment is not None
            and row.ressayre_assessment.status == "trace_rejected"
            and row.ressayre_assessment.certificate is None,
            "trace-rejected row has inconsistent assessment data",
        )
    elif row.status == "symbolic_zero_rejected":
        symbolic = row.symbolic_fallback
        _require(
            row.ressayre_assessment is not None
            and row.ressayre_assessment.status == "evaluation_unresolved"
            and row.ressayre_assessment.certificate is None,
            "symbolic-zero row has inconsistent evaluation data",
        )
        _require(
            symbolic is not None
            and symbolic.status == "symbolic_zero_rejected"
            and symbolic.certificate is None,
            "symbolic-zero row lacks its exact negative outcome",
        )
    elif row.status == "evaluation_unresolved":
        _require(
            row.ressayre_assessment is not None
            and row.ressayre_assessment.status == "evaluation_unresolved"
            and row.ressayre_assessment.certificate is None,
            "unresolved row has inconsistent evaluation data",
        )
        _require(
            row.symbolic_fallback is None
            or row.symbolic_fallback.status
            in {"resource_limit_unresolved", "exception_unresolved"},
            "unresolved row has an inconsistent symbolic outcome",
        )
    elif row.status == "inadmissible":
        _require(
            row.ressayre_assessment is None and row.inadmissibility_reason is not None,
            "inadmissible row lacks its exact reason",
        )
    else:
        raise CandidateSystemInvariantError(f"unknown screening status: {row.status}")


def _validate_screening(
    screening: candidate_pipeline.CandidateScreeningResult,
) -> tuple[tuple[int, candidate_pipeline.ScreenedCandidateRow], ...]:
    """Replay binding, conservation, uniqueness, and canonical-order invariants."""
    _require(
        isinstance(screening, candidate_pipeline.CandidateScreeningResult),
        "screening must be a CandidateScreeningResult",
    )
    d = screening.d
    _require(7 <= d <= 13, "candidate-system finalization supports ranks 7 through 13")
    _require(
        screening.full_dimension.d == d
        and screening.full_dimension.is_full_dimensional
        and verify_full_dimension_certificate(screening.full_dimension),
        "ambient full-dimension certificate failed exact replay",
    )
    _require(
        screening.cyclic_cross_checks_enabled is False,
        "candidate-wide cyclic checks must be deferred until after reduction",
    )
    _require(
        screening.include_modular_cross_checks is False,
        "candidate-wide modular checks must be deferred until after reduction",
    )

    normalized = screening.normalized_taus
    _require(
        normalized == tuple(sorted(set(normalized))), "normalized taus are not unique and sorted"
    )
    for tau in normalized:
        _require(len(tau) == d, "normalized tau has the wrong rank")
        _require(all(type(value) is int for value in tau), "tau entries must be integers")
        _require(primitive_tau(tau) == tau, "normalized tau is not primitive")

    expected_candidates = tuple(tau for tau in normalized if not is_pauli_lower_tau(tau))
    expected_structural = tuple(tau for tau in normalized if is_pauli_lower_tau(tau))
    _require(
        tuple(row.tau for row in screening.candidates) == expected_candidates,
        "candidate rows do not match the canonical normalized tau order",
    )
    _require(
        tuple(row.tau for row in screening.structural_rows) == expected_structural,
        "structural rows do not match the normalized tau partition",
    )

    for row in screening.structural_rows:
        _require(
            row.constraint == constraint_from_tau(row.tau, 3),
            "structural constraint is not bound to its tau",
        )
        _require(
            tau_from_constraint(row.constraint, 3) == row.tau,
            "structural affine row fails homogeneous round trip",
        )

    certified: list[tuple[int, candidate_pipeline.ScreenedCandidateRow]] = []
    status_counts = {
        "certified": 0,
        "trace_rejected": 0,
        "symbolic_zero_rejected": 0,
        "evaluation_unresolved": 0,
        "inadmissible": 0,
    }
    for index, row in enumerate(screening.candidates):
        _require(
            row.constraint == constraint_from_tau(row.tau, 3),
            "candidate constraint is not bound to its tau",
        )
        _require(
            tau_from_constraint(row.constraint, 3) == row.tau,
            "candidate affine row fails homogeneous round trip",
        )
        _require(row.status in status_counts, f"unknown screening status: {row.status}")
        status_counts[row.status] += 1
        _require(
            row.cyclic_schubert is None,
            "cyclic cross-check was not deferred on a screened candidate",
        )
        if row.status == "certified":
            _certified_ressayre_certificate(row)
            certified.append((index, row))
        else:
            _validate_noncertified_row(row)

    counts = screening.counts
    _require(
        counts.input_tau_rows >= 0 and counts.duplicate_or_rescaled_rows >= 0,
        "raw normalization counts must be nonnegative",
    )
    _require(counts.normalized_tau_rows == len(normalized), "normalized row count mismatch")
    _require(
        counts.structural_rows == len(screening.structural_rows),
        "structural row count mismatch",
    )
    _require(
        counts.nonstructural_candidate_rows == len(screening.candidates),
        "candidate row count mismatch",
    )
    _require(
        counts.normalized_tau_rows == counts.structural_rows + counts.nonstructural_candidate_rows,
        "normalized structural/candidate count conservation failed",
    )
    _require(
        counts.input_tau_rows == counts.normalized_tau_rows + counts.duplicate_or_rescaled_rows,
        "raw-to-normalized count conservation failed",
    )
    _require(counts.certified_rows == status_counts["certified"], "certified count mismatch")
    _require(
        counts.trace_rejected_rows == status_counts["trace_rejected"],
        "trace rejection count mismatch",
    )
    _require(
        counts.symbolic_zero_rejected_rows == status_counts["symbolic_zero_rejected"],
        "symbolic-zero rejection count mismatch",
    )
    _require(
        counts.evaluation_unresolved_rows == status_counts["evaluation_unresolved"],
        "unresolved count mismatch",
    )
    _require(
        counts.inadmissible_rows == status_counts["inadmissible"], "inadmissible count mismatch"
    )
    _require(
        counts.exact_rejected_rows
        == status_counts["trace_rejected"]
        + status_counts["symbolic_zero_rejected"]
        + status_counts["inadmissible"],
        "exact rejection count mismatch",
    )
    _require(
        counts.resolved_candidate_rows == status_counts["certified"] + counts.exact_rejected_rows,
        "resolved count mismatch",
    )
    _require(
        counts.resolved_candidate_rows + counts.evaluation_unresolved_rows
        == counts.nonstructural_candidate_rows,
        "screening outcome counts do not conserve candidates",
    )
    symbolic = tuple(
        row.symbolic_fallback for row in screening.candidates if row.symbolic_fallback is not None
    )
    _require(
        counts.symbolic_fallback_attempted_rows == len(symbolic),
        "symbolic fallback attempt count mismatch",
    )
    _require(
        counts.symbolic_fallback_certified_rows
        == sum(outcome.status == "certified" for outcome in symbolic),
        "symbolic certification count mismatch",
    )
    _require(
        counts.symbolic_fallback_resource_limited_rows
        == sum(outcome.status == "resource_limit_unresolved" for outcome in symbolic),
        "symbolic resource-limit count mismatch",
    )
    _require(
        counts.symbolic_fallback_exception_unresolved_rows
        == sum(outcome.status == "exception_unresolved" for outcome in symbolic),
        "symbolic exception count mismatch",
    )
    _require(
        counts.symbolic_fallback_attempted_rows
        == counts.symbolic_fallback_certified_rows
        + counts.symbolic_zero_rejected_rows
        + counts.symbolic_fallback_resource_limited_rows
        + counts.symbolic_fallback_exception_unresolved_rows,
        "symbolic fallback outcomes do not conserve attempts",
    )
    _require(
        counts.cyclic_schubert_cross_checked_rows == 0
        and counts.cyclic_schubert_deferred_rows == status_counts["certified"],
        "cyclic deferral counts do not conserve certified rows",
    )
    _require(
        counts.cyclic_schubert_exact_coefficient_rows == 0
        and counts.cyclic_schubert_nonzero_rows == 0
        and counts.cyclic_schubert_zero_rows == 0
        and counts.cyclic_schubert_coefficient_one_rows == 0
        and counts.cyclic_schubert_unavailable_rows == 0,
        "deferred cyclic screen contains nonzero result counts",
    )
    _require(
        screening.all_candidates_resolved == (status_counts["evaluation_unresolved"] == 0),
        "all-candidates-resolved flag is inconsistent",
    )
    return tuple(certified)


def _known_system_regression(
    d: int, retained_rows: tuple[IntegralConstraint, ...]
) -> KnownSystemRegression:
    retained_taus = {tau_from_constraint(row, 3) for row in retained_rows}
    if d > 10:
        return KnownSystemRegression(
            available=False,
            expected_count=None,
            retained_count=len(retained_rows),
            oriented_sets_equal=None,
            generated_only=(),
            published_only=(),
            lookup_source=None,
        )
    table = constraints(3, d)
    published_taus = {
        tau_from_constraint(
            IntegralConstraint(tuple(row["coeffs"]), int(row["rhs"])),
            3,
        )
        for row in table["inequalities"]
    }
    return KnownSystemRegression(
        available=True,
        expected_count=len(published_taus),
        retained_count=len(retained_rows),
        oriented_sets_equal=retained_taus == published_taus,
        generated_only=tuple(sorted(retained_taus - published_taus)),
        published_only=tuple(sorted(published_taus - retained_taus)),
        lookup_source=table["source"],
    )


def _cross_check_retained_row(
    constraint: IntegralConstraint,
    *,
    include_modular: bool,
    modulus: int,
) -> CyclicSchubertVerdict:
    verdict = certify_cyclic_schubert(
        constraint.coeffs,
        constraint.rhs,
        particle_number=3,
        include_modular=include_modular,
        modulus=modulus,
    )
    if verdict.exact_certificate is not None:
        _require(
            verify_exact_cyclic_schubert_certificate(verdict.exact_certificate),
            "retained exact cyclic-Schubert certificate failed replay",
        )
    if verdict.modular_certificate is not None:
        _require(
            verify_modular_cyclic_schubert_certificate(verdict.modular_certificate),
            "retained modular cyclic-Schubert certificate failed replay",
        )
    return verdict


def _counts(
    screening: candidate_pipeline.CandidateScreeningResult,
    reduction: OrderedSliceReductionCertificate | None,
    retained_cross_checks: int,
) -> CandidateSystemCounts:
    if reduction is None:
        redundancy_input = 0
        redundant = 0
        retained = 0
    else:
        redundancy_input = len(reduction.candidate_rows)
        redundant = len(reduction.removals)
        retained = len(reduction.retained_indices)
    return CandidateSystemCounts(
        normalized_tau_rows=screening.counts.normalized_tau_rows,
        structural_rows=screening.counts.structural_rows,
        screened_candidate_rows=screening.counts.nonstructural_candidate_rows,
        certified_rows=screening.counts.certified_rows,
        exact_rejected_rows=screening.counts.exact_rejected_rows,
        unresolved_rows=screening.counts.evaluation_unresolved_rows,
        redundancy_input_rows=redundancy_input,
        redundant_rows=redundant,
        retained_rows=retained,
        retained_cyclic_cross_checks=retained_cross_checks,
    )


def compose_candidate_system(
    screening: candidate_pipeline.CandidateScreeningResult,
    *,
    include_modular: bool = False,
    modulus: int = DEFAULT_MODULUS,
) -> CandidateSystemResult:
    """Screening-to-facet composition with exact replay and fail-closed gates.

    Unresolved candidates produce a serializable blocked result and no
    reduction.  A rank-7 through rank-10 lookup mismatch likewise preserves
    its exact reduction as diagnostics but exposes no ``final_system``.
    Operational certificate-search failures propagate without being converted
    into mathematical rejection claims.
    """
    if type(include_modular) is not bool:
        raise TypeError("include_modular must be bool")
    certified = _validate_screening(screening)
    if not screening.all_candidates_resolved:
        return CandidateSystemResult(
            d=screening.d,
            status=BLOCKED_UNRESOLVED_STATUS,
            screening=screening,
            reduction=None,
            retained_rows=(),
            known_system_regression=None,
            counts=_counts(screening, None, 0),
            include_modular_cross_checks=include_modular,
        )

    certified_rows = tuple(row.constraint for _index, row in certified)
    reduction = certify_ordered_slice_redundancy(certified_rows, screening.d)
    _require(
        verify_ordered_slice_reduction_certificate(reduction),
        "ordered-slice reduction certificate failed exact replay",
    )
    certified_indices = set(range(len(certified_rows)))
    removed_indices = set(reduction.redundant_indices)
    retained_indices = set(reduction.retained_indices)
    _require(
        removed_indices.isdisjoint(retained_indices)
        and removed_indices | retained_indices == certified_indices,
        "reduction indices do not partition certified candidates",
    )
    _require(
        len(reduction.removals) + len(reduction.retained_indices) == len(certified_rows),
        "reduction count conservation failed",
    )

    regression = _known_system_regression(screening.d, reduction.retained_rows)
    if not regression.passed:
        return CandidateSystemResult(
            d=screening.d,
            status=BLOCKED_KNOWN_MISMATCH_STATUS,
            screening=screening,
            reduction=reduction,
            retained_rows=(),
            known_system_regression=regression,
            counts=_counts(screening, reduction, 0),
            include_modular_cross_checks=include_modular,
        )

    facets = {facet.target_index: facet for facet in reduction.facets}
    final_rows: list[FinalCandidateRow] = []
    for reduction_index in reduction.retained_indices:
        screening_index, screened_row = certified[reduction_index]
        certificate = _certified_ressayre_certificate(screened_row)
        facet = facets.get(reduction_index)
        _require(facet is not None, "retained row lacks its exact facet witness")
        cyclic = _cross_check_retained_row(
            screened_row.constraint,
            include_modular=include_modular,
            modulus=modulus,
        )
        final_rows.append(
            FinalCandidateRow(
                screening_index=screening_index,
                reduction_index=reduction_index,
                tau=screened_row.tau,
                constraint=screened_row.constraint,
                ressayre_resolution_method=screened_row.resolution_method or "",
                ressayre_certificate=certificate,
                facet_certificate=facet,
                cyclic_schubert=cyclic,
            )
        )
    result = CandidateSystemResult(
        d=screening.d,
        status=(FINALIZED_KNOWN_STATUS if regression.available else FINALIZED_RELATIVE_STATUS),
        screening=screening,
        reduction=reduction,
        retained_rows=tuple(final_rows),
        known_system_regression=regression,
        counts=_counts(screening, reduction, len(final_rows)),
        include_modular_cross_checks=include_modular,
    )
    _require(
        len(result.retained_rows) == result.counts.retained_rows,
        "final retained-row binding count mismatch",
    )
    return result


__all__ = [
    "BLOCKED_KNOWN_MISMATCH_STATUS",
    "BLOCKED_UNRESOLVED_STATUS",
    "FINALIZED_KNOWN_STATUS",
    "FINALIZED_RELATIVE_STATUS",
    "CandidateSystemCounts",
    "CandidateSystemInvariantError",
    "CandidateSystemResult",
    "FinalCandidateRow",
    "KnownSystemRegression",
    "compose_candidate_system",
]
