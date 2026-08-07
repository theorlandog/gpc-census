"""Exact screening of supplied fixed-``N=3`` homogeneous candidates.

This module screens a candidate collection produced by the external
Bulois-Denis-Ressayre backend.  It does not generate that collection and
makes no claim that the supplied rows are exhaustive, facet defining, or
irredundant.

For each primitive, oriented, nonstructural row ``tau.lambda <= 0`` the
screen has five possible outcomes:

``certified``
    Admissibility and the Ressayre trace condition hold, and a replayed exact
    nonzero determinant evaluation is bound to the oriented row.
``trace_rejected``
    Admissibility holds but the exact trace counts differ.
``evaluation_unresolved``
    Admissibility and trace hold, but the finite deterministic evaluation set
    did not find a nonzero determinant and no exact symbolic fallback resolved
    the row.  This is inconclusive, never a vanishing proof.
``symbolic_zero_rejected``
    Admissibility and trace hold, and exact symbolic determinant construction
    proves that the determinant polynomial is identically zero.
``inadmissible``
    The row does not define a codimension-one exterior-weight hyperplane.

When requested, the cyclic-Schubert coefficient is evaluated only after a row
has a bound Ressayre determinant certificate.  It is an independent
cross-validation datum and never changes the Ressayre status.  The default
superset screen defers this potentially expensive test until exact redundancy
reduction.  In particular, a deferred test, coefficient zero, a positive
nonunit coefficient, and an unavailable cyclic test never become rejection
rules.

Correctness follows row by row.  Positive gcd normalization preserves the
oriented halfspace, and set deduplication removes only identical primitive
rows.  The structural pattern is checked explicitly.  Every positive
Ressayre decision carries a determinant evaluation that is replayed and bound
to the oriented primitive row.  Trace mismatch and inadmissibility are exact
negative decisions.  Finite evaluation exhaustion remains unresolved unless
the symbolic determinant routine either constructs a replayed nonzero witness
or returns exact zero after trace equality has been replayed.  Conservation
assertions then prove that every normalized row appears in exactly one outcome
class.

For ``r`` supplied rows of width ``d``, primitive normalization and sorting
cost ``O(r d log r)`` integer comparisons.  For one trace-matched candidate
with ``m`` weights on its hyperplane and a ``q`` by ``q`` tangent matrix,
``a`` deterministic evaluations use ``O(a(q m d + q^3))`` operations in the
current dense weight representation, apart from integer bit growth.  Symbolic
determinant expansion can have exponential expression growth, so it is
optional and has an explicit matrix order resource gate.  Cyclic-Schubert
work is also optional and normally deferred until redundancy reduction has
shrunk the candidate set.
"""

from __future__ import annotations

from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import partial
import multiprocessing as mp

from .bdr import (
    Tau,
    constraint_from_tau,
    is_pauli_lower_tau,
    normalize_tau_system,
    parse_python_tau_export,
)
from .cyclic_schubert import (
    DEFAULT_MODULUS,
    CyclicSchubertVerdict,
    _is_prime_below_2_64,
    certify_cyclic_schubert,
    verify_exact_cyclic_schubert_certificate,
    verify_modular_cyclic_schubert_certificate,
)
from .exterior import exterior_weights, negative_roots
from .full_dimension import (
    FullDimensionCertificate,
    full_dimension_certificate,
    verify_full_dimension_certificate,
)
from .model import IntegralConstraint, RessayreCertificate, RessayreEvaluationVerdict
from .ressayre import (
    admissible_candidate_from_tau,
    assess_tau_by_evaluation,
    certify_candidate,
    tangent_determinant_is_identically_zero,
    tangent_determinant_is_structurally_zero,
    verify_tau_ressayre_certificate,
)


@dataclass(frozen=True)
class StructuralPauliRow:
    """A separated dominant-chamber lower bound ``lambda_d >= 0``."""

    tau: Tau
    constraint: IntegralConstraint

    def as_dict(self) -> dict[str, object]:
        """Return the explicit structural proof payload."""
        return {
            "tau": list(self.tau),
            "constraint": self.constraint.as_dict(),
            "status": "structural_pauli_lower_bound",
            "proof": "tau is the primitive row (0,...,0,-1), hence lambda_d >= 0",
            "ressayre_assessment": None,
            "cyclic_schubert_cross_check": None,
        }


@dataclass(frozen=True)
class SymbolicFallbackOutcome:
    """Outcome of an attempted exact symbolic determinant fallback.

    ``resource_limit_unresolved`` and ``exception_unresolved`` are operational
    outcomes only.  Neither is mathematical evidence that the determinant
    vanishes.
    """

    status: str
    reason: str | None
    certificate: RessayreCertificate | None
    certificate_replayed_and_bound: bool | None

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-serializable fallback payload."""
        return {
            "status": self.status,
            "reason": self.reason,
            "certificate": (
                None if self.certificate is None else self.certificate.as_dict()
            ),
            "certificate_replayed_and_bound": self.certificate_replayed_and_bound,
        }


@dataclass(frozen=True)
class ScreenedCandidateRow:
    """One nonstructural row and its exact or explicitly unresolved outcome."""

    tau: Tau
    constraint: IntegralConstraint
    status: str
    ressayre_assessment: RessayreEvaluationVerdict | None
    inadmissibility_reason: str | None
    resolution_method: str | None
    certificate_replayed_and_bound: bool | None
    symbolic_fallback: SymbolicFallbackOutcome | None
    cyclic_schubert: CyclicSchubertVerdict | None

    @property
    def decision_class(self) -> str:
        """Coarsen statuses within the Ressayre screen.

        ``exact_rejected`` means exact failure of the implemented Ressayre
        criterion. It is not, by itself, a proof that the halfspace is false.
        """
        if self.status == "certified":
            return "certified"
        if self.status in {
            "trace_rejected",
            "inadmissible",
            "symbolic_zero_rejected",
        }:
            return "exact_rejected"
        if self.status == "evaluation_unresolved":
            return "unresolved"
        raise AssertionError(f"unknown candidate status: {self.status}")

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-serializable row payload."""
        if self.cyclic_schubert is not None:
            cyclic_status = "cross_checked"
        elif self.status == "certified":
            cyclic_status = "deferred"
        else:
            cyclic_status = "not_applicable"
        return {
            "tau": list(self.tau),
            "constraint": self.constraint.as_dict(),
            "status": self.status,
            "decision_class": self.decision_class,
            "inadmissibility_reason": self.inadmissibility_reason,
            "ressayre_assessment": (
                None
                if self.ressayre_assessment is None
                else self.ressayre_assessment.as_dict()
            ),
            "resolution_method": self.resolution_method,
            "certificate_replayed_and_bound": self.certificate_replayed_and_bound,
            "symbolic_fallback": (
                None if self.symbolic_fallback is None else self.symbolic_fallback.as_dict()
            ),
            "cyclic_schubert_cross_check": (
                None if self.cyclic_schubert is None else self.cyclic_schubert.as_dict()
            ),
            "cyclic_schubert_status": cyclic_status,
        }


@dataclass(frozen=True)
class CandidateScreeningCounts:
    """Conservation counts for normalization and every screening outcome."""

    input_tau_rows: int
    normalized_tau_rows: int
    duplicate_or_rescaled_rows: int
    structural_rows: int
    nonstructural_candidate_rows: int
    certified_rows: int
    trace_rejected_rows: int
    symbolic_zero_rejected_rows: int
    evaluation_unresolved_rows: int
    inadmissible_rows: int
    symbolic_fallback_attempted_rows: int
    symbolic_fallback_certified_rows: int
    symbolic_fallback_resource_limited_rows: int
    symbolic_fallback_exception_unresolved_rows: int
    cyclic_schubert_cross_checked_rows: int
    cyclic_schubert_deferred_rows: int
    cyclic_schubert_exact_coefficient_rows: int
    cyclic_schubert_nonzero_rows: int
    cyclic_schubert_zero_rows: int
    cyclic_schubert_coefficient_one_rows: int
    cyclic_schubert_unavailable_rows: int

    @property
    def exact_rejected_rows(self) -> int:
        """Return all exact negative decisions."""
        return (
            self.trace_rejected_rows
            + self.symbolic_zero_rejected_rows
            + self.inadmissible_rows
        )

    @property
    def resolved_candidate_rows(self) -> int:
        """Return candidates with a positive or exact negative decision."""
        return self.certified_rows + self.exact_rejected_rows

    def as_dict(self) -> dict[str, int]:
        """Return stable count names for artifact serialization."""
        return {
            "input_tau_rows": self.input_tau_rows,
            "normalized_tau_rows": self.normalized_tau_rows,
            "duplicate_or_rescaled_rows": self.duplicate_or_rescaled_rows,
            "structural_rows": self.structural_rows,
            "nonstructural_candidate_rows": self.nonstructural_candidate_rows,
            "certified_rows": self.certified_rows,
            "trace_rejected_rows": self.trace_rejected_rows,
            "symbolic_zero_rejected_rows": self.symbolic_zero_rejected_rows,
            "evaluation_unresolved_rows": self.evaluation_unresolved_rows,
            "inadmissible_rows": self.inadmissible_rows,
            "symbolic_fallback_attempted_rows": self.symbolic_fallback_attempted_rows,
            "symbolic_fallback_certified_rows": self.symbolic_fallback_certified_rows,
            "symbolic_fallback_resource_limited_rows": (
                self.symbolic_fallback_resource_limited_rows
            ),
            "symbolic_fallback_exception_unresolved_rows": (
                self.symbolic_fallback_exception_unresolved_rows
            ),
            "exact_rejected_rows": self.exact_rejected_rows,
            "resolved_candidate_rows": self.resolved_candidate_rows,
            "cyclic_schubert_cross_checked_rows": (
                self.cyclic_schubert_cross_checked_rows
            ),
            "cyclic_schubert_deferred_rows": self.cyclic_schubert_deferred_rows,
            "cyclic_schubert_exact_coefficient_rows": (
                self.cyclic_schubert_exact_coefficient_rows
            ),
            "cyclic_schubert_nonzero_rows": self.cyclic_schubert_nonzero_rows,
            "cyclic_schubert_zero_rows": self.cyclic_schubert_zero_rows,
            "cyclic_schubert_coefficient_one_rows": (
                self.cyclic_schubert_coefficient_one_rows
            ),
            "cyclic_schubert_unavailable_rows": self.cyclic_schubert_unavailable_rows,
        }


@dataclass(frozen=True)
class CandidateScreeningResult:
    """Deterministic screening payload for one supplied census system."""

    d: int
    normalized_taus: tuple[Tau, ...]
    candidates: tuple[ScreenedCandidateRow, ...]
    structural_rows: tuple[StructuralPauliRow, ...]
    counts: CandidateScreeningCounts
    full_dimension: FullDimensionCertificate
    ressayre_attempts: int
    cyclic_cross_checks_enabled: bool
    include_modular_cross_checks: bool
    symbolic_fallback_enabled: bool
    symbolic_max_determinant_order: int | None
    n: int = 3

    @property
    def all_candidates_resolved(self) -> bool:
        """Whether no finite determinant search remains inconclusive."""
        return self.counts.evaluation_unresolved_rows == 0

    def as_dict(self) -> dict[str, object]:
        """Return the complete screening artifact for one rank."""
        return {
            "particle_number": self.n,
            "rank": self.d,
            "counts": self.counts.as_dict(),
            "all_candidates_resolved": self.all_candidates_resolved,
            "ressayre_evaluation_attempt_limit": self.ressayre_attempts,
            "cyclic_cross_checks_enabled": self.cyclic_cross_checks_enabled,
            "include_modular_cross_checks": self.include_modular_cross_checks,
            "symbolic_fallback_enabled": self.symbolic_fallback_enabled,
            "symbolic_max_determinant_order": self.symbolic_max_determinant_order,
            "normalized_taus": [list(tau) for tau in self.normalized_taus],
            "ambient_moment_cone_full_dimension": {
                "established": self.full_dimension.is_full_dimensional,
                "certificate_replayed": True,
                "certificate": self.full_dimension.as_dict(),
            },
            "candidates": [row.as_dict() for row in self.candidates],
            "structural_pauli_rows": [row.as_dict() for row in self.structural_rows],
            "claim_scope": {
                "supplied_row_screening": "exact_or_explicitly_unresolved",
                "symbolic_zero_rejection": "exact_determinant_identity",
                "candidate_exhaustiveness": "not_established",
                "facetness": "not_established",
                "redundancy": "not_established",
                "cyclic_schubert_role": (
                    "optional_cross_validation_only_never_a_rejection_gate"
                ),
            },
        }


class UnresolvedCandidatesError(RuntimeError):
    """Raised by the optional fail gate while retaining the full result."""

    def __init__(self, result: CandidateScreeningResult):
        self.result = result
        count = result.counts.evaluation_unresolved_rows
        super().__init__(f"{count} candidate row(s) remain evaluation_unresolved")


def _coerce_raw_taus(raw_taus: Sequence[Sequence[int]]) -> tuple[Tau, ...]:
    rows = tuple(tuple(row) for row in raw_taus)
    if not rows:
        raise ValueError("raw_taus must not be empty")
    if any(not row for row in rows):
        raise ValueError("tau rows must not be empty")
    if any(type(value) is not int for row in rows for value in row):
        raise TypeError("tau entries must be integers, excluding bool")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("tau rows have inconsistent dimensions")
    if width < 3:
        raise ValueError("screening requires rank at least three")
    if any(not any(row) for row in rows):
        raise ValueError("tau rows must be nonzero")
    return rows


def _cross_check_cyclic(
    constraint: IntegralConstraint,
    *,
    include_modular: bool,
    modulus: int,
    particle_number: int = 3,
) -> CyclicSchubertVerdict:
    verdict = certify_cyclic_schubert(
        constraint.coeffs,
        constraint.rhs,
        particle_number=particle_number,
        include_modular=include_modular,
        modulus=modulus,
    )
    if verdict.exact_certificate is not None and not (
        verify_exact_cyclic_schubert_certificate(verdict.exact_certificate)
    ):
        raise AssertionError("generated exact cyclic-Schubert certificate failed replay")
    if verdict.modular_certificate is not None and not (
        verify_modular_cyclic_schubert_certificate(verdict.modular_certificate)
    ):
        raise AssertionError("generated modular cyclic-Schubert certificate failed replay")
    return verdict


def _attempt_symbolic_fallback(
    tau: Tau,
    assessment: RessayreEvaluationVerdict,
    *,
    max_determinant_order: int | None,
    particle_number: int = 3,
) -> SymbolicFallbackOutcome:
    """Try exact symbolic certification after finite evaluation exhaustion.

    The trace equality has already been checked by ``assessment``.  Therefore
    ``certify_candidate`` returning ``None`` proves that the exact determinant
    polynomial is zero.  Operational limits and implementation exceptions are
    recorded as unresolved outcomes instead of being converted to rejection.
    """
    if assessment.status != "evaluation_unresolved":
        raise ValueError("symbolic fallback requires evaluation_unresolved")
    try:
        candidate = admissible_candidate_from_tau(particle_number, tau)
        weights = exterior_weights(particle_number, len(tau))
        exact_below_count = sum(
            sum(left * right for left, right in zip(candidate.h, weight, strict=True))
            < candidate.z
            for weight in weights
        )
        exact_negative_root_count = sum(
            candidate.h[root[1]] - candidate.h[root[0]] < 0
            for root in negative_roots(len(tau))
        )
    except Exception as error:
        return SymbolicFallbackOutcome(
            status="exception_unresolved",
            reason=f"{type(error).__name__}: {error}",
            certificate=None,
            certificate_replayed_and_bound=None,
        )
    if (
        assessment.below_weight_count != exact_below_count
        or assessment.negative_root_count != exact_negative_root_count
    ):
        raise AssertionError("evaluation assessment trace counts failed exact replay")
    if exact_below_count != exact_negative_root_count:
        raise AssertionError("evaluation_unresolved assessment failed exact trace equality")
    # STRUCTURAL ZEROS FIRST. A support bipartite graph with no perfect
    # matching makes every permutation term of the determinant contain an empty
    # cell, so the determinant is identically zero and Hall's theorem hands back
    # a checkable violator. That decides 86.3 percent of the rank-7 rows that
    # reach here and 97.1 percent at rank 8, in O(E sqrt(V)) against a dynamic
    # program that is exponential in the matrix order: 0.08 seconds for all
    # 6,607 rank-8 survivors, against 224 seconds of DP.
    #
    # The outcome is deliberately INDISTINGUISHABLE from the one the dynamic
    # program returns, same status and same reason, so stored manifests stay
    # byte identical and this is a pure accelerator rather than a change of
    # verdict. The Hall violator is recorded by the orbit-native constructor,
    # which is where a certificate belongs; adding a field here would rewrite
    # every published row payload for no mathematical gain.
    #
    # It runs BEFORE the resource limit because it is never resource limited,
    # so a configured limit can now resolve rows it previously abandoned. With
    # the default limit of None nothing changes at all.
    try:
        structural, _violator = tangent_determinant_is_structurally_zero(
            particle_number, len(tau), candidate
        )
    except Exception as error:
        return SymbolicFallbackOutcome(
            status="exception_unresolved",
            reason=f"{type(error).__name__}: {error}",
            certificate=None,
            certificate_replayed_and_bound=None,
        )
    if structural:
        return SymbolicFallbackOutcome(
            status="symbolic_zero_rejected",
            reason="exact sparse Ressayre determinant is identically zero",
            certificate=None,
            certificate_replayed_and_bound=None,
        )

    determinant_order = exact_negative_root_count
    if (
        max_determinant_order is not None
        and determinant_order > max_determinant_order
    ):
        return SymbolicFallbackOutcome(
            status="resource_limit_unresolved",
            reason=(
                f"determinant order {determinant_order} exceeds configured limit "
                f"{max_determinant_order}"
            ),
            certificate=None,
            certificate_replayed_and_bound=None,
        )

    try:
        determinant_is_zero = tangent_determinant_is_identically_zero(
            particle_number, len(tau), candidate
        )
    except Exception as error:
        return SymbolicFallbackOutcome(
            status="exception_unresolved",
            reason=f"{type(error).__name__}: {error}",
            certificate=None,
            certificate_replayed_and_bound=None,
        )

    if determinant_is_zero:
        return SymbolicFallbackOutcome(
            status="symbolic_zero_rejected",
            reason="exact sparse Ressayre determinant is identically zero",
            certificate=None,
            certificate_replayed_and_bound=None,
        )
    try:
        certificate = certify_candidate(particle_number, len(tau), candidate)
    except Exception as error:
        return SymbolicFallbackOutcome(
            status="exception_unresolved",
            reason=f"{type(error).__name__}: {error}",
            certificate=None,
            certificate_replayed_and_bound=None,
        )
    if certificate is None:
        raise AssertionError(
            "sparse nonzero determinant disagrees with symbolic certification"
        )
    if not verify_tau_ressayre_certificate(particle_number, tau, certificate):
        raise AssertionError("symbolic Ressayre certificate failed replay or tau binding")
    return SymbolicFallbackOutcome(
        status="certified",
        reason=None,
        certificate=certificate,
        certificate_replayed_and_bound=True,
    )


def _validate_evaluation_assessment(
    assessment: RessayreEvaluationVerdict,
    *,
    attempt_limit: int,
) -> None:
    """Fail closed on an internally inconsistent tri-state assessment."""
    integer_fields = (
        assessment.attempts,
        assessment.below_weight_count,
        assessment.negative_root_count,
    )
    if any(type(value) is not int or value < 0 for value in integer_fields):
        raise AssertionError("Ressayre assessment counts must be nonnegative integers")
    trace_equal = assessment.below_weight_count == assessment.negative_root_count
    if assessment.status == "trace_rejected":
        if assessment.attempts != 0 or trace_equal or assessment.certificate is not None:
            raise AssertionError("trace_rejected assessment is internally inconsistent")
        return
    if assessment.status == "evaluation_unresolved":
        if (
            assessment.attempts != attempt_limit
            or not trace_equal
            or assessment.certificate is not None
        ):
            raise AssertionError("evaluation_unresolved assessment is internally inconsistent")
        return
    if assessment.status != "certified":
        raise AssertionError(f"unknown Ressayre assessment status: {assessment.status}")
    certificate = assessment.certificate
    if (
        not 1 <= assessment.attempts <= attempt_limit
        or not trace_equal
        or certificate is None
    ):
        raise AssertionError("certified assessment is internally inconsistent")
    if (
        len(certificate.below_indices) != assessment.below_weight_count
        or len(certificate.negative_roots) != assessment.negative_root_count
    ):
        raise AssertionError("certified assessment counts disagree with its certificate")


def _screen_candidate(
    tau: Tau,
    *,
    ressayre_attempts: int,
    cyclic_cross_check: bool,
    include_modular: bool,
    modulus: int,
    symbolic_fallback_enabled: bool,
    symbolic_max_determinant_order: int | None,
    particle_number: int = 3,
) -> ScreenedCandidateRow:
    constraint = constraint_from_tau(tau, particle_number)
    try:
        assessment = assess_tau_by_evaluation(
            particle_number,
            tau,
            attempts=ressayre_attempts,
        )
    except ValueError:
        try:
            admissible_candidate_from_tau(particle_number, tau)
        except ValueError as admissibility_error:
            return ScreenedCandidateRow(
                tau=tau,
                constraint=constraint,
                status="inadmissible",
                ressayre_assessment=None,
                inadmissibility_reason=str(admissibility_error),
                resolution_method="admissibility",
                certificate_replayed_and_bound=None,
                symbolic_fallback=None,
                cyclic_schubert=None,
            )
        raise

    _validate_evaluation_assessment(assessment, attempt_limit=ressayre_attempts)
    if assessment.status == "trace_rejected":
        return ScreenedCandidateRow(
            tau=tau,
            constraint=constraint,
            status=assessment.status,
            ressayre_assessment=assessment,
            inadmissibility_reason=None,
            resolution_method="trace_count",
            certificate_replayed_and_bound=None,
            symbolic_fallback=None,
            cyclic_schubert=None,
        )

    if assessment.status == "evaluation_unresolved":
        if not symbolic_fallback_enabled:
            return ScreenedCandidateRow(
                tau=tau,
                constraint=constraint,
                status="evaluation_unresolved",
                ressayre_assessment=assessment,
                inadmissibility_reason=None,
                resolution_method=None,
                certificate_replayed_and_bound=None,
                symbolic_fallback=None,
                cyclic_schubert=None,
            )
        symbolic = _attempt_symbolic_fallback(
            tau,
            assessment,
            max_determinant_order=symbolic_max_determinant_order,
            particle_number=particle_number,
        )
        if symbolic.status == "certified":
            cyclic = (
                _cross_check_cyclic(
                    constraint,
                    include_modular=include_modular,
                    modulus=modulus,
                    particle_number=particle_number,
                )
                if cyclic_cross_check
                else None
            )
            return ScreenedCandidateRow(
                tau=tau,
                constraint=constraint,
                status="certified",
                ressayre_assessment=assessment,
                inadmissibility_reason=None,
                resolution_method="exact_symbolic_determinant",
                certificate_replayed_and_bound=True,
                symbolic_fallback=symbolic,
                cyclic_schubert=cyclic,
            )
        if symbolic.status == "symbolic_zero_rejected":
            return ScreenedCandidateRow(
                tau=tau,
                constraint=constraint,
                status="symbolic_zero_rejected",
                ressayre_assessment=assessment,
                inadmissibility_reason=None,
                resolution_method="exact_symbolic_determinant",
                certificate_replayed_and_bound=None,
                symbolic_fallback=symbolic,
                cyclic_schubert=None,
            )
        if symbolic.status not in {
            "resource_limit_unresolved",
            "exception_unresolved",
        }:
            raise AssertionError(f"unknown symbolic fallback status: {symbolic.status}")
        return ScreenedCandidateRow(
            tau=tau,
            constraint=constraint,
            status="evaluation_unresolved",
            ressayre_assessment=assessment,
            inadmissibility_reason=None,
            resolution_method=None,
            certificate_replayed_and_bound=None,
            symbolic_fallback=symbolic,
            cyclic_schubert=None,
        )

    certificate = assessment.certificate
    if certificate is None:
        raise AssertionError("validated certified assessment lost its certificate")
    if not verify_tau_ressayre_certificate(particle_number, tau, certificate):
        raise AssertionError("Ressayre certificate failed replay or oriented tau binding")
    cyclic = (
        _cross_check_cyclic(
            constraint,
            include_modular=include_modular,
            modulus=modulus,
            particle_number=particle_number,
        )
        if cyclic_cross_check
        else None
    )
    return ScreenedCandidateRow(
        tau=tau,
        constraint=constraint,
        status="certified",
        ressayre_assessment=assessment,
        inadmissibility_reason=None,
        resolution_method="deterministic_evaluation",
        certificate_replayed_and_bound=True,
        symbolic_fallback=None,
        cyclic_schubert=cyclic,
    )


def _counts(
    *,
    input_count: int,
    normalized_count: int,
    candidates: tuple[ScreenedCandidateRow, ...],
    structural_count: int,
) -> CandidateScreeningCounts:
    cyclic = tuple(row.cyclic_schubert for row in candidates if row.cyclic_schubert)
    exact_cyclic = tuple(row for row in cyclic if row.coefficient is not None)
    symbolic = tuple(
        row.symbolic_fallback for row in candidates if row.symbolic_fallback is not None
    )
    return CandidateScreeningCounts(
        input_tau_rows=input_count,
        normalized_tau_rows=normalized_count,
        duplicate_or_rescaled_rows=input_count - normalized_count,
        structural_rows=structural_count,
        nonstructural_candidate_rows=len(candidates),
        certified_rows=sum(row.status == "certified" for row in candidates),
        trace_rejected_rows=sum(row.status == "trace_rejected" for row in candidates),
        symbolic_zero_rejected_rows=sum(
            row.status == "symbolic_zero_rejected" for row in candidates
        ),
        evaluation_unresolved_rows=sum(
            row.status == "evaluation_unresolved" for row in candidates
        ),
        inadmissible_rows=sum(row.status == "inadmissible" for row in candidates),
        symbolic_fallback_attempted_rows=len(symbolic),
        symbolic_fallback_certified_rows=sum(
            row.status == "certified" for row in symbolic
        ),
        symbolic_fallback_resource_limited_rows=sum(
            row.status == "resource_limit_unresolved" for row in symbolic
        ),
        symbolic_fallback_exception_unresolved_rows=sum(
            row.status == "exception_unresolved" for row in symbolic
        ),
        cyclic_schubert_cross_checked_rows=len(cyclic),
        cyclic_schubert_deferred_rows=(
            sum(row.status == "certified" for row in candidates) - len(cyclic)
        ),
        cyclic_schubert_exact_coefficient_rows=len(exact_cyclic),
        cyclic_schubert_nonzero_rows=sum(row.coefficient != 0 for row in exact_cyclic),
        cyclic_schubert_zero_rows=sum(row.coefficient == 0 for row in exact_cyclic),
        cyclic_schubert_coefficient_one_rows=sum(
            row.coefficient_one for row in exact_cyclic
        ),
        cyclic_schubert_unavailable_rows=len(cyclic) - len(exact_cyclic),
    )


def _verify_count_conservation(result: CandidateScreeningResult, input_count: int) -> None:
    counts = result.counts
    outcome_count = (
        counts.certified_rows
        + counts.trace_rejected_rows
        + counts.symbolic_zero_rejected_rows
        + counts.evaluation_unresolved_rows
        + counts.inadmissible_rows
    )
    if outcome_count != counts.nonstructural_candidate_rows:
        raise AssertionError("candidate outcome counts do not conserve rows")
    if counts.structural_rows + counts.nonstructural_candidate_rows != (
        counts.normalized_tau_rows
    ):
        raise AssertionError("normalized structural and candidate counts do not conserve rows")
    if counts.normalized_tau_rows + counts.duplicate_or_rescaled_rows != input_count:
        raise AssertionError("normalization counts do not conserve input rows")
    if (
        counts.cyclic_schubert_cross_checked_rows
        + counts.cyclic_schubert_deferred_rows
        != counts.certified_rows
    ):
        raise AssertionError("cyclic cross-check and deferral counts do not conserve rows")
    if counts.symbolic_fallback_certified_rows > counts.certified_rows:
        raise AssertionError("symbolic certifications exceed all certified rows")
    if (
        counts.symbolic_fallback_certified_rows
        + counts.symbolic_zero_rejected_rows
        + counts.symbolic_fallback_resource_limited_rows
        + counts.symbolic_fallback_exception_unresolved_rows
        != counts.symbolic_fallback_attempted_rows
    ):
        raise AssertionError("symbolic fallback outcome counts do not conserve attempts")
    if (
        counts.cyclic_schubert_exact_coefficient_rows
        + counts.cyclic_schubert_unavailable_rows
        != counts.cyclic_schubert_cross_checked_rows
    ):
        raise AssertionError("cyclic availability counts do not conserve cross-checks")
    if (
        counts.cyclic_schubert_nonzero_rows + counts.cyclic_schubert_zero_rows
        != counts.cyclic_schubert_exact_coefficient_rows
    ):
        raise AssertionError("cyclic coefficient counts do not conserve exact evaluations")


def screen_tau_candidates(
    raw_taus: Sequence[Sequence[int]],
    *,
    particle_number: int = 3,
    ressayre_attempts: int = 32,
    cyclic_cross_check: bool = False,
    include_modular: bool = False,
    modulus: int = DEFAULT_MODULUS,
    symbolic_fallback: bool | None = None,
    symbolic_max_determinant_order: int | None = None,
    workers: int = 1,
    require_resolved: bool = False,
) -> CandidateScreeningResult:
    """Normalize and screen supplied homogeneous rows for one census system.

    ``particle_number`` defaults to three, which is the fixed-``N=3`` generator
    program.  Nothing below this function is special to ``N = 3``: the weights,
    roots, tangent determinant, and full-dimension witness are all written for
    ``wedge^n C^d``, so the wider census systems screen by the same code path
    and produce the same kind of exact Ressayre certificate.

    ``require_resolved`` is an optional operational gate.  It is checked only
    after all rows have been screened, and its exception retains the complete
    result so unresolved rows can still be serialized and resumed.
    """
    if type(particle_number) is not int or particle_number <= 0:
        raise ValueError("particle_number must be a positive integer")
    if type(ressayre_attempts) is not int or ressayre_attempts <= 0:
        raise ValueError("ressayre_attempts must be a positive integer")
    if type(cyclic_cross_check) is not bool:
        raise TypeError("cyclic_cross_check must be bool")
    if type(include_modular) is not bool:
        raise TypeError("include_modular must be bool")
    if include_modular and not cyclic_cross_check:
        raise ValueError("include_modular requires cyclic_cross_check")
    if include_modular:
        if type(modulus) is not int:
            raise TypeError("modulus must be an integer, excluding bool")
        if modulus >= 2**64 or not _is_prime_below_2_64(modulus):
            raise ValueError("modulus must be a prime smaller than 2^64")
    if symbolic_fallback is not None and type(symbolic_fallback) is not bool:
        raise TypeError("symbolic_fallback must be bool or None")
    if (
        symbolic_max_determinant_order is not None
        and (
            type(symbolic_max_determinant_order) is not int
            or symbolic_max_determinant_order <= 0
        )
    ):
        raise ValueError("symbolic_max_determinant_order must be a positive integer or None")
    if type(require_resolved) is not bool:
        raise TypeError("require_resolved must be bool")
    if type(workers) is not int or workers <= 0:
        raise ValueError("workers must be a positive integer")
    raw = _coerce_raw_taus(raw_taus)
    normalized = normalize_tau_system(raw, particle_number)
    d = normalized.d
    if include_modular and modulus < d:
        raise ValueError("modulus must be at least the rank for distinct base residues")
    symbolic_fallback_enabled = d <= 13 if symbolic_fallback is None else symbolic_fallback

    ambient = full_dimension_certificate(d, particle_number)
    if not verify_full_dimension_certificate(ambient):
        raise AssertionError("ambient full-dimension certificate failed replay")

    structural_taus = tuple(
        tau for tau in normalized.raw_taus if is_pauli_lower_tau(tau)
    )
    candidate_taus = tuple(
        tau for tau in normalized.raw_taus if not is_pauli_lower_tau(tau)
    )
    structural = tuple(
        StructuralPauliRow(tau=tau, constraint=constraint_from_tau(tau, particle_number))
        for tau in structural_taus
    )
    screen_one = partial(
        _screen_candidate,
        particle_number=particle_number,
        ressayre_attempts=ressayre_attempts,
        cyclic_cross_check=cyclic_cross_check,
        include_modular=include_modular,
        modulus=modulus,
        symbolic_fallback_enabled=symbolic_fallback_enabled,
        symbolic_max_determinant_order=symbolic_max_determinant_order,
    )
    if workers == 1 or len(candidate_taus) < 2:
        candidates = tuple(screen_one(tau) for tau in candidate_taus)
    else:
        chunk_size = max(1, len(candidate_taus) // (workers * 8))
        start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=mp.get_context(start_method),
        ) as executor:
            candidates = tuple(
                executor.map(screen_one, candidate_taus, chunksize=chunk_size)
            )
    counts = _counts(
        input_count=len(raw),
        normalized_count=len(normalized.raw_taus),
        candidates=candidates,
        structural_count=len(structural),
    )
    result = CandidateScreeningResult(
        d=d,
        n=particle_number,
        normalized_taus=normalized.raw_taus,
        candidates=candidates,
        structural_rows=structural,
        counts=counts,
        full_dimension=ambient,
        ressayre_attempts=ressayre_attempts,
        cyclic_cross_checks_enabled=cyclic_cross_check,
        include_modular_cross_checks=include_modular,
        symbolic_fallback_enabled=symbolic_fallback_enabled,
        symbolic_max_determinant_order=symbolic_max_determinant_order,
    )
    _verify_count_conservation(result, len(raw))
    if require_resolved and not result.all_candidates_resolved:
        raise UnresolvedCandidatesError(result)
    return result


def screen_bdr_export(
    text: str,
    *,
    expected_rank: int | None = None,
    particle_number: int = 3,
    ressayre_attempts: int = 32,
    cyclic_cross_check: bool = False,
    include_modular: bool = False,
    modulus: int = DEFAULT_MODULUS,
    symbolic_fallback: bool | None = None,
    symbolic_max_determinant_order: int | None = None,
    workers: int = 1,
    require_resolved: bool = False,
) -> CandidateScreeningResult:
    """Parse one BDR Python export without executing it, then screen its rows."""
    raw_taus = parse_python_tau_export(text)
    if expected_rank is not None:
        if type(expected_rank) is not int or expected_rank < 3:
            raise ValueError("expected_rank must be an integer at least three")
        if any(len(tau) != expected_rank for tau in raw_taus):
            raise ValueError(f"BDR export does not contain rank-{expected_rank} rows")
    return screen_tau_candidates(
        raw_taus,
        particle_number=particle_number,
        ressayre_attempts=ressayre_attempts,
        cyclic_cross_check=cyclic_cross_check,
        include_modular=include_modular,
        modulus=modulus,
        symbolic_fallback=symbolic_fallback,
        symbolic_max_determinant_order=symbolic_max_determinant_order,
        workers=workers,
        require_resolved=require_resolved,
    )


__all__ = [
    "CandidateScreeningCounts",
    "CandidateScreeningResult",
    "ScreenedCandidateRow",
    "StructuralPauliRow",
    "SymbolicFallbackOutcome",
    "UnresolvedCandidatesError",
    "screen_bdr_export",
    "screen_tau_candidates",
]
