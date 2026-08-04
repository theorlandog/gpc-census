"""Immutable exact data models for constraint generation."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


Weight = tuple[int, ...]
Root = tuple[int, int]


@dataclass(frozen=True, order=True)
class IntegralConstraint:
    """An integral affine relation ``coeffs . lambda <= rhs`` or equality."""

    coeffs: tuple[int, ...]
    rhs: int

    def as_dict(self) -> dict[str, object]:
        return {"coeffs": list(self.coeffs), "rhs": self.rhs}


@dataclass(frozen=True, order=True)
class AdmissibleHyperplane:
    """A primitive pair ``(H, z)`` with ``sum(H)=0``.

    The hyperplane is ``H . x = z``.  ``witness`` gives affinely independent
    exterior-weight indices spanning it.
    """

    h: tuple[int, ...]
    z: int
    witness: tuple[int, ...]


@dataclass(frozen=True, order=True)
class RessayreCertificate:
    """Replayable nonvanishing certificate for ``H . lambda >= z``."""

    h: tuple[int, ...]
    z: int
    admissibility_witness: tuple[int, ...]
    on_indices: tuple[int, ...]
    below_indices: tuple[int, ...]
    negative_roots: tuple[Root, ...]
    determinant_evaluation: tuple[int, ...]
    determinant_value: int

    def as_dict(self) -> dict[str, object]:
        return {
            "h": list(self.h),
            "z": self.z,
            "inequality": "sum_i h_i * lambda_i >= z",
            "admissibility_witness": list(self.admissibility_witness),
            "on_indices": list(self.on_indices),
            "below_indices": list(self.below_indices),
            "negative_roots": [list(root) for root in self.negative_roots],
            "determinant_evaluation": list(self.determinant_evaluation),
            "determinant_value": self.determinant_value,
        }


@dataclass(frozen=True)
class RessayreEvaluationVerdict:
    """Tri-state result of deterministic determinant evaluation.

    ``trace_rejected`` is an exact rejection of the Ressayre trace condition.
    ``evaluation_unresolved`` means only that the finite deterministic point
    set found no witness and must never be interpreted as polynomial
    vanishing.  ``certified`` carries an exact nonzero determinant witness.
    """

    status: str
    attempts: int
    below_weight_count: int
    negative_root_count: int
    certificate: RessayreCertificate | None

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic JSON representation of the verdict."""
        return {
            "status": self.status,
            "attempts": self.attempts,
            "below_weight_count": self.below_weight_count,
            "negative_root_count": self.negative_root_count,
            "certificate": (
                None if self.certificate is None else self.certificate.as_dict()
            ),
        }


@dataclass(frozen=True)
class RessayreEnumeration:
    """Exact counts and accepted certificates from exhaustive enumeration."""

    admissible_hyperplane_count: int
    oriented_candidate_count: int
    trace_candidate_count: int
    certificates: tuple[RessayreCertificate, ...]


@dataclass(frozen=True)
class SchubertCertificate:
    """A finite divided-difference certificate for a Klyachko inequality."""

    test_spectrum: tuple[int, ...]
    source_permutation: tuple[int, ...]
    target_permutation: tuple[int, ...]
    induced_spectrum: tuple[tuple[int, int], ...]
    divided_difference_word: tuple[int, ...]
    schubert_coefficient: int
    constraint: IntegralConstraint

    def as_dict(self) -> dict[str, object]:
        return {
            "test_spectrum": list(self.test_spectrum),
            "source_permutation": list(self.source_permutation),
            "target_permutation": list(self.target_permutation),
            "induced_spectrum_value_multiplicity": [
                list(pair) for pair in self.induced_spectrum
            ],
            "divided_difference_word_zero_based": list(self.divided_difference_word),
            "schubert_coefficient": self.schubert_coefficient,
            "constraint": self.constraint.as_dict(),
        }


@dataclass(frozen=True)
class GeneratedSystem:
    """A generated GPC system plus its exact inner/outer closure evidence."""

    n: int
    d: int
    inequalities: tuple[IntegralConstraint, ...]
    equalities: tuple[IntegralConstraint, ...]
    structural_inequalities: tuple[IntegralConstraint, ...]
    outer_vertices: tuple[tuple[Fraction, ...], ...]
    inner_point_count: int
    outer_vertices_in_inner_cloud: bool
    enumeration: RessayreEnumeration
    schubert_cross_checks: tuple[SchubertCertificate, ...]
