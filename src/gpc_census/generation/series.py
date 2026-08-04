"""Legacy coefficient-one regression harness for supplied final systems.

This module is intentionally stricter than the general candidate pipeline. It
accepts only rows whose selected cyclic-Schubert coefficient is one, which is
useful for reproducing the known final-system fixtures covered by its tests.
It is incomplete as a candidate screen and can reject a valid row with a
positive nonunit coefficient. New exhaustive-candidate and frontier work must
use :mod:`candidate_pipeline` followed by :mod:`candidate_system`, where the
cyclic result is cross-validation and never a validity gate.

The Bulois-Denis-Ressayre bridge performs only exact normalization.  This
module adds independent cyclic-Schubert and Ressayre determinant gates.  A
returned rank payload has therefore passed these deterministic steps:

1. primitive oriented tau normalization and deduplication;
2. conversion to the fixed-trace affine slice;
3. separation of structural Pauli lower bounds from nontrivial rows;
4. exact integer evaluation of every cyclic-Schubert coefficient;
5. rejection unless every coefficient is exactly one; and
6. exact admissibility, trace, and nonzero Ressayre determinant checks for
   every nonstructural row.

This layer does not claim that an upstream tau list is complete.  It certifies
specific arithmetic and representation-theoretic conditions for the rows that
it receives.  Exact redundancy and candidate exhaustiveness remain separate
proof obligations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .bdr import Tau, normalize_tau_system, tau_from_constraint
from .cyclic_schubert import (
    DEFAULT_MODULUS,
    CyclicSchubertVerdict,
    certify_cyclic_schubert,
)
from .model import IntegralConstraint, RessayreCertificate
from .ressayre import assess_tau_by_evaluation, verify_tau_ressayre_certificate


@dataclass(frozen=True)
class CertifiedAffineRow:
    """One normalized affine row with independently replayable checks."""

    tau: Tau
    constraint: IntegralConstraint
    schubert: CyclicSchubertVerdict
    row_kind: str
    ressayre: RessayreCertificate | None

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-serializable certificate payload."""
        return {
            "tau": list(self.tau),
            "constraint": self.constraint.as_dict(),
            "row_kind": self.row_kind,
            "cyclic_schubert": self.schubert.as_dict(),
            "ressayre": None if self.ressayre is None else self.ressayre.as_dict(),
        }


@dataclass(frozen=True)
class CertifiedRankCounts:
    """Normalization and certification counts for one rank."""

    input_tau_rows: int
    normalized_tau_rows: int
    duplicate_or_rescaled_rows: int
    inequality_rows: int
    structural_rows: int
    cyclic_schubert_coefficient_one_rows: int
    ressayre_certified_inequality_rows: int

    def as_dict(self) -> dict[str, int]:
        """Return counts with stable field names for artifact serialization."""
        return {
            "input_tau_rows": self.input_tau_rows,
            "normalized_tau_rows": self.normalized_tau_rows,
            "duplicate_or_rescaled_rows": self.duplicate_or_rescaled_rows,
            "inequality_rows": self.inequality_rows,
            "structural_rows": self.structural_rows,
            "cyclic_schubert_coefficient_one_rows": (
                self.cyclic_schubert_coefficient_one_rows
            ),
            "ressayre_certified_inequality_rows": (
                self.ressayre_certified_inequality_rows
            ),
        }


@dataclass(frozen=True)
class CertifiedRankSystem:
    """Deterministic coefficient-one payload for one ``(N, d)`` system."""

    particle_number: int
    d: int
    normalized_taus: tuple[Tau, ...]
    inequalities: tuple[CertifiedAffineRow, ...]
    structural_inequalities: tuple[CertifiedAffineRow, ...]
    counts: CertifiedRankCounts
    ressayre_attempts: int

    @property
    def all_coefficients_one(self) -> bool:
        """Return whether every normalized affine row has exact coefficient one."""
        rows = self.inequalities + self.structural_inequalities
        return len(rows) == self.counts.cyclic_schubert_coefficient_one_rows and all(
            row.schubert.coefficient_one for row in rows
        )

    @property
    def all_nonstructural_ressayre_certified(self) -> bool:
        """Return whether every nonstructural row has a bound exact witness."""
        return (
            len(self.inequalities)
            == self.counts.ressayre_certified_inequality_rows
            and all(row.ressayre is not None for row in self.inequalities)
            and all(row.ressayre is None for row in self.structural_inequalities)
        )

    def as_dict(self) -> dict[str, object]:
        """Return a stable, JSON-serializable per-rank payload."""
        return {
            "particle_number": self.particle_number,
            "rank": self.d,
            "counts": self.counts.as_dict(),
            "all_coefficients_one": self.all_coefficients_one,
            "all_nonstructural_ressayre_certified": (
                self.all_nonstructural_ressayre_certified
            ),
            "ressayre_evaluation_attempt_limit": self.ressayre_attempts,
            "normalized_taus": [list(tau) for tau in self.normalized_taus],
            "inequalities": [row.as_dict() for row in self.inequalities],
            "structural_inequalities": [
                row.as_dict() for row in self.structural_inequalities
            ],
        }


def _coerce_raw_taus(raw_taus: Sequence[Sequence[int]]) -> tuple[Tau, ...]:
    rows = tuple(tuple(row) for row in raw_taus)
    for row in rows:
        if any(type(value) is not int for value in row):
            raise TypeError("tau entries must be integers")
    return rows


def _certify_rows(
    rows: tuple[IntegralConstraint, ...],
    *,
    particle_number: int,
    include_modular: bool,
    modulus: int,
    require_ressayre: bool,
    ressayre_attempts: int,
) -> tuple[CertifiedAffineRow, ...]:
    certified: list[CertifiedAffineRow] = []
    for constraint in rows:
        tau = tau_from_constraint(constraint, particle_number)
        verdict = certify_cyclic_schubert(
            constraint.coeffs,
            constraint.rhs,
            particle_number=particle_number,
            include_modular=include_modular,
            modulus=modulus,
        )
        if not verdict.coefficient_one:
            raise ValueError(
                "cyclic-Schubert coefficient-one requirement failed for "
                f"{constraint.as_dict()}: status={verdict.status}, "
                f"coefficient={verdict.coefficient}"
            )
        ressayre: RessayreCertificate | None = None
        if require_ressayre:
            assessment = assess_tau_by_evaluation(
                particle_number,
                tau,
                attempts=ressayre_attempts,
            )
            ressayre = assessment.certificate
            if ressayre is None:
                raise ValueError(
                    "Ressayre certification failed for "
                    f"tau={tau}: status={assessment.status}, "
                    f"attempts={assessment.attempts}, "
                    f"below={assessment.below_weight_count}, "
                    f"roots={assessment.negative_root_count}"
                )
            if not verify_tau_ressayre_certificate(
                particle_number, tau, ressayre
            ):
                raise AssertionError("generated Ressayre certificate failed replay")
        certified.append(
            CertifiedAffineRow(
                tau=tau,
                constraint=constraint,
                schubert=verdict,
                row_kind=(
                    "nonstructural_candidate" if require_ressayre else "pauli_lower_bound"
                ),
                ressayre=ressayre,
            )
        )
    return tuple(certified)


def certify_tau_rank_system(
    raw_taus: Sequence[Sequence[int]],
    particle_number: int,
    *,
    include_modular: bool = False,
    modulus: int = DEFAULT_MODULUS,
    ressayre_attempts: int = 32,
) -> CertifiedRankSystem:
    """Normalize and exactly certify every row in one raw tau system.

    The function returns only after every normalized affine row has cyclic-
    Schubert coefficient exactly one and every nonstructural row has a bound,
    replayed Ressayre determinant witness.  Modular certificates can be
    included as deterministic cross-checks, but the coefficient-one decision
    always comes from the exact integer evaluator.

    Correctness is compositional.  ``normalize_tau_system`` preserves each
    oriented inequality on the trace slice and deterministically separates
    the Pauli lower row.  ``certify_cyclic_schubert`` then computes the exact
    characteristic-zero coefficient for each affine representative.  The
    Ressayre evaluator independently proves admissibility, trace balance, and
    determinant nonvanishing for each nonstructural row.  These facts prove
    row validity under the cited Ressayre theorem, but do not prove that the
    supplied list is exhaustive or irredundant.  Sorting and immutable tuple
    storage make output independent of input row order.

    For ``r`` raw rows, normalization costs ``O(r log r)`` comparisons plus
    linear gcd work.  Certification costs the sum of the memoized
    divided-difference evaluations documented in ``cyclic_schubert``, plus at
    most ``ressayre_attempts`` exact determinants per nonstructural row.
    """
    if ressayre_attempts <= 0:
        raise ValueError("ressayre_attempts must be positive")
    raw = _coerce_raw_taus(raw_taus)
    normalized = normalize_tau_system(raw, particle_number)
    inequalities = _certify_rows(
        normalized.inequalities,
        particle_number=particle_number,
        include_modular=include_modular,
        modulus=modulus,
        require_ressayre=True,
        ressayre_attempts=ressayre_attempts,
    )
    structural = _certify_rows(
        normalized.structural_inequalities,
        particle_number=particle_number,
        include_modular=include_modular,
        modulus=modulus,
        require_ressayre=False,
        ressayre_attempts=ressayre_attempts,
    )
    normalized_count = len(normalized.raw_taus)
    counts = CertifiedRankCounts(
        input_tau_rows=len(raw),
        normalized_tau_rows=normalized_count,
        duplicate_or_rescaled_rows=len(raw) - normalized_count,
        inequality_rows=len(inequalities),
        structural_rows=len(structural),
        cyclic_schubert_coefficient_one_rows=len(inequalities) + len(structural),
        ressayre_certified_inequality_rows=sum(
            row.ressayre is not None for row in inequalities
        ),
    )
    result = CertifiedRankSystem(
        particle_number=particle_number,
        d=normalized.d,
        normalized_taus=normalized.raw_taus,
        inequalities=inequalities,
        structural_inequalities=structural,
        counts=counts,
        ressayre_attempts=ressayre_attempts,
    )
    if not result.all_coefficients_one:
        raise AssertionError("successful rank payload contains an uncertified row")
    if not result.all_nonstructural_ressayre_certified:
        raise AssertionError("successful rank payload lacks a Ressayre certificate")
    return result


__all__ = [
    "CertifiedAffineRow",
    "CertifiedRankCounts",
    "CertifiedRankSystem",
    "certify_tau_rank_system",
]
