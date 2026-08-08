"""Exact redundancy and facet certificates on the ordered pure-state slice.

This module treats a supplied list of affine rows as a polyhedral object.  It
does not prove that the rows are valid generalized Pauli constraints.  Given
that separate validity input, it proves exactly which supplied rows are
redundant and which define facets after adjoining the fixed structural system

``sum(lambda_i) = N``, ``lambda_0 <= 1``,
``lambda_0 >= ... >= lambda_(d-1)``, and ``lambda_(d-1) >= 0``.

For a target row ``a.x <= b``, redundancy is certified by rational Farkas
multipliers ``y_j >= 0`` and a free trace multiplier ``z`` satisfying

``a = sum_j y_j A_j + z * 1`` and
``b = sum_j y_j b_j + N*z + slack``, with ``slack >= 0``.

Every point satisfying the premise rows and the trace equation then satisfies
the target.  Removal certificates are replayed in order, so later certificates
may use rows that are removed later without creating a circular argument.

For every retained row, the certificate stores a rational point satisfying all
other retained rows but strictly violating the target.  One shared rational
point satisfies the final system strictly.  The segment between these two
points crosses the target hyperplane at a point strict for all other rows.
Thus the target face has relative codimension one in the ``d-1`` dimensional
trace slice and is a facet.  This also proves irredundancy directly.

Floating-point linear programming is used only to discover supports and
witnesses.  Every returned certificate is reconstructed over ``Q`` and checked
with integer and rational arithmetic.  Failure to reconstruct raises
``CertificateSearchError``; no numerical result is ever accepted as proof.

If ``r`` candidate rows are supplied, there are ``d+1`` structural rows.  The
implementation solves ``O(r)`` linear programs with ``O(r+d)`` constraints or
dual variables and ``O(d)`` equations.  Exact replay costs ``O(r^2 d)``
rational operations in the worst case, while a basic Farkas certificate has at
most ``d+1`` nonzero terms and the serialized payload has size ``O(rd)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

import numpy as np
import sympy as sp
from scipy.optimize import linprog

from .model import IntegralConstraint


Point = tuple[Fraction, ...]
MultiplierTerm = tuple[int, Fraction]


class CertificateSearchError(RuntimeError):
    """Raised when numerical discovery cannot be reconstructed exactly."""


def _fraction_payload(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _point_payload(point: Point) -> list[list[int]]:
    return [_fraction_payload(value) for value in point]


@dataclass(frozen=True)
class FarkasRemovalCertificate:
    """Exact implication proving one candidate redundant at its removal step."""

    target_index: int
    structural_terms: tuple[MultiplierTerm, ...]
    candidate_terms: tuple[MultiplierTerm, ...]
    trace_multiplier: Fraction
    rhs_slack: Fraction

    def as_dict(self) -> dict[str, object]:
        return {
            "target_index": self.target_index,
            "structural_terms": [
                {"index": index, "multiplier": _fraction_payload(multiplier)}
                for index, multiplier in self.structural_terms
            ],
            "candidate_terms": [
                {"index": index, "multiplier": _fraction_payload(multiplier)}
                for index, multiplier in self.candidate_terms
            ],
            "trace_multiplier": _fraction_payload(self.trace_multiplier),
            "rhs_slack": _fraction_payload(self.rhs_slack),
        }


@dataclass(frozen=True)
class FacetWitnessCertificate:
    """A rational witness that violates one retained row and no other one."""

    target_index: int
    violating_point: Point
    violation_margin: Fraction

    def as_dict(self) -> dict[str, object]:
        return {
            "target_index": self.target_index,
            "violating_point": _point_payload(self.violating_point),
            "violation_margin": _fraction_payload(self.violation_margin),
        }


@dataclass(frozen=True)
class OrderedSliceReductionCertificate:
    """Replayable reduction and facet proof for one ordered-slice system."""

    particle_number: int
    d: int
    candidate_rows: tuple[IntegralConstraint, ...]
    structural_rows: tuple[IntegralConstraint, ...]
    removals: tuple[FarkasRemovalCertificate, ...]
    retained_indices: tuple[int, ...]
    relative_interior_point: Point
    facets: tuple[FacetWitnessCertificate, ...]

    @property
    def retained_rows(self) -> tuple[IntegralConstraint, ...]:
        return tuple(self.candidate_rows[index] for index in self.retained_indices)

    @property
    def redundant_indices(self) -> tuple[int, ...]:
        return tuple(step.target_index for step in self.removals)

    def as_dict(self) -> dict[str, object]:
        return {
            "particle_number": self.particle_number,
            "rank": self.d,
            "candidate_rows": [row.as_dict() for row in self.candidate_rows],
            "structural_rows": [row.as_dict() for row in self.structural_rows],
            "removals": [step.as_dict() for step in self.removals],
            "retained_indices": list(self.retained_indices),
            "redundant_indices": list(self.redundant_indices),
            "relative_interior_point": _point_payload(self.relative_interior_point),
            "facets": [facet.as_dict() for facet in self.facets],
        }


def ordered_pauli_inequalities(d: int) -> tuple[IntegralConstraint, ...]:
    """Return the canonical Pauli and dominance rows for ``d`` orbitals.

    These are ``lambda_0 <= 1``, the dominance chain, and
    ``lambda_(d-1) >= 0``. None of them mentions the particle number: occupation
    numbers lie in ``[0, 1]`` and are ordered whatever ``N`` is, so the rank
    window here was an ``N=3`` habit rather than a mathematical limit. The
    particle number enters only through the trace equation, which is carried
    separately as the free multiplier.
    """
    if type(d) is not int or d < 2:
        raise ValueError("the ordered slice requires at least two orbitals")

    rows = [IntegralConstraint((1,) + (0,) * (d - 1), 1)]
    for index in range(d - 1):
        coefficients = [0] * d
        coefficients[index] = -1
        coefficients[index + 1] = 1
        rows.append(IntegralConstraint(tuple(coefficients), 0))
    rows.append(IntegralConstraint((0,) * (d - 1) + (-1,), 0))
    return tuple(rows)


def _validate_rows(rows: Sequence[IntegralConstraint], d: int) -> tuple[IntegralConstraint, ...]:
    result = tuple(rows)
    for row in result:
        if not isinstance(row, IntegralConstraint):
            raise TypeError("candidate rows must be IntegralConstraint instances")
        if len(row.coeffs) != d:
            raise ValueError("candidate row has the wrong rank")
        if any(type(value) is not int for value in row.coeffs) or type(row.rhs) is not int:
            raise TypeError("candidate rows must have integer coefficients")
        if not any(row.coeffs):
            raise ValueError("candidate row must have a nonzero normal")
    return result


def _lhs(row: IntegralConstraint, point: Point) -> Fraction:
    return sum(
        Fraction(coefficient) * coordinate
        for coefficient, coordinate in zip(row.coeffs, point, strict=True)
    )


def _satisfies(row: IntegralConstraint, point: Point) -> bool:
    return _lhs(row, point) <= row.rhs


def _strictly_satisfies(row: IntegralConstraint, point: Point) -> bool:
    return _lhs(row, point) < row.rhs


def _to_fraction(value: sp.Rational) -> Fraction:
    return Fraction(int(sp.numer(value)), int(sp.denom(value)))


def _rational_substitutions(
    expressions: tuple[sp.Expr, ...],
    symbols: tuple[sp.Symbol, ...],
    approximations: tuple[float, ...],
) -> tuple[tuple[Fraction, ...], ...]:
    free = set().union(*(expression.free_symbols for expression in expressions))
    results: list[tuple[Fraction, ...]] = []
    for denominator_bound in (1, 10**3, 10**6, 10**9, 10**12):
        substitutions: dict[sp.Symbol, sp.Rational] = {}
        for index, symbol in enumerate(symbols):
            if symbol not in free:
                continue
            approximation = Fraction(float(approximations[index])).limit_denominator(
                denominator_bound
            )
            substitutions[symbol] = sp.Rational(approximation.numerator, approximation.denominator)
        values = tuple(
            _to_fraction(sp.Rational(expression.subs(substitutions))) for expression in expressions
        )
        if values not in results:
            results.append(values)
    return tuple(results)


def _solve_supported_farkas(
    target: IntegralConstraint,
    premises: tuple[tuple[str, int, IntegralConstraint], ...],
    floating: np.ndarray,
    particle_number: int,
) -> tuple[tuple[Fraction, ...], Fraction, Fraction] | None:
    d = len(target.coeffs)
    premise_count = len(premises)
    configurations: list[tuple[tuple[int, ...], bool, bool]] = []
    for tolerance in (1e-7, 1e-9, 1e-11, 1e-13):
        support = tuple(
            index for index, value in enumerate(floating[:premise_count]) if value > tolerance
        )
        z_nonzero = abs(float(floating[premise_count])) > tolerance
        slack_nonzero = float(floating[premise_count + 1]) > tolerance
        for include_z in (True, z_nonzero):
            for include_slack in (True, slack_nonzero):
                configuration = (support, include_z, include_slack)
                if configuration not in configurations:
                    configurations.append(configuration)

    target_vector = sp.Matrix([*target.coeffs, target.rhs])
    for support, include_z, include_slack in configurations:
        columns: list[sp.Matrix] = []
        approximations: list[float] = []
        for premise_index in support:
            row = premises[premise_index][2]
            columns.append(sp.Matrix([*row.coeffs, row.rhs]))
            approximations.append(float(floating[premise_index]))
        if include_z:
            columns.append(sp.Matrix([*[1] * d, particle_number]))
            approximations.append(float(floating[premise_count]))
        if include_slack:
            columns.append(sp.Matrix([*[0] * d, 1]))
            approximations.append(float(floating[premise_count + 1]))
        if not columns:
            continue

        matrix = sp.Matrix.hstack(*columns)
        symbols = sp.symbols(f"q:{len(columns)}")
        solutions = sp.linsolve((matrix, target_vector), symbols)
        if solutions is sp.EmptySet:
            continue
        expressions = tuple(next(iter(solutions)))
        for values in _rational_substitutions(expressions, symbols, tuple(approximations)):
            offset = len(support)
            multipliers = [Fraction(0)] * premise_count
            for premise_index, value in zip(support, values[:offset], strict=True):
                multipliers[premise_index] = value
            z = values[offset] if include_z else Fraction(0)
            if include_z:
                offset += 1
            slack = values[offset] if include_slack else Fraction(0)
            if any(value < 0 for value in multipliers) or slack < 0:
                continue
            lhs_coefficients = tuple(
                sum(
                    multipliers[index] * premises[index][2].coeffs[coordinate]
                    for index in range(premise_count)
                )
                + z
                for coordinate in range(d)
            )
            lhs_rhs = (
                sum(multipliers[index] * premises[index][2].rhs for index in range(premise_count))
                + z * particle_number
                + slack
            )
            if lhs_coefficients == target.coeffs and lhs_rhs == target.rhs:
                return tuple(multipliers), z, slack
    return None


def _find_farkas_removal(
    target_index: int,
    target: IntegralConstraint,
    structural: tuple[IntegralConstraint, ...],
    candidates: tuple[IntegralConstraint, ...],
    active_other_indices: tuple[int, ...],
    particle_number: int,
) -> FarkasRemovalCertificate | None:
    premises = tuple(
        [("structural", index, row) for index, row in enumerate(structural)]
        + [("candidate", index, candidates[index]) for index in active_other_indices]
    )
    d = len(target.coeffs)
    columns = [[*row.coeffs, row.rhs] for _kind, _index, row in premises]
    columns.append([*[1] * d, particle_number])
    columns.append([*[0] * d, 1])
    equality_matrix = np.asarray(columns, dtype=float).T
    equality_rhs = np.asarray([*target.coeffs, target.rhs], dtype=float)
    premise_count = len(premises)
    result = linprog(
        np.zeros(premise_count + 2),
        A_eq=equality_matrix,
        b_eq=equality_rhs,
        bounds=[(0, None)] * premise_count + [(None, None), (0, None)],
        method="highs-ds",
    )
    if not result.success:
        return None

    exact = _solve_supported_farkas(target, premises, result.x, particle_number)
    if exact is None:
        return None
    multipliers, trace_multiplier, rhs_slack = exact
    structural_terms: list[MultiplierTerm] = []
    candidate_terms: list[MultiplierTerm] = []
    for premise, multiplier in zip(premises, multipliers, strict=True):
        if multiplier == 0:
            continue
        kind, index, _row = premise
        if kind == "structural":
            structural_terms.append((index, multiplier))
        else:
            candidate_terms.append((index, multiplier))
    return FarkasRemovalCertificate(
        target_index=target_index,
        structural_terms=tuple(structural_terms),
        candidate_terms=tuple(candidate_terms),
        trace_multiplier=trace_multiplier,
        rhs_slack=rhs_slack,
    )


def _exact_vertex_from_discovery(
    rows: tuple[IntegralConstraint, ...],
    approximation: np.ndarray,
    particle_number: int,
) -> Point | None:
    d = len(approximation)
    order = sorted(
        range(len(rows)),
        key=lambda index: abs(
            rows[index].rhs
            - sum(
                coefficient * approximation[coordinate]
                for coordinate, coefficient in enumerate(rows[index].coeffs)
            )
        ),
    )
    matrix = sp.Matrix([[1] * d])
    selected: list[int] = []
    current_rank = 1
    for index in order:
        extended = matrix.col_join(sp.Matrix([rows[index].coeffs]))
        extended_rank = int(extended.rank())
        if extended_rank == current_rank:
            continue
        matrix = extended
        selected.append(index)
        current_rank = extended_rank
        if current_rank == d:
            break
    if current_rank != d:
        return None
    rhs = sp.Matrix([particle_number, *(rows[index].rhs for index in selected)])
    solution = matrix.inv() * rhs
    return tuple(_to_fraction(sp.Rational(value)) for value in solution)


def _find_violating_point(
    target: IntegralConstraint,
    premises: tuple[IntegralConstraint, ...],
    particle_number: int,
) -> tuple[Point, Fraction] | None:
    d = len(target.coeffs)
    result = linprog(
        -np.asarray(target.coeffs, dtype=float),
        A_ub=np.asarray([row.coeffs for row in premises], dtype=float),
        b_ub=np.asarray([row.rhs for row in premises], dtype=float),
        A_eq=np.ones((1, d)),
        b_eq=np.asarray([particle_number], dtype=float),
        bounds=[(None, None)] * d,
        method="highs-ds",
    )
    if not result.success:
        return None
    point = _exact_vertex_from_discovery(premises, result.x, particle_number)
    if point is None or not all(_satisfies(row, point) for row in premises):
        return None
    margin = _lhs(target, point) - target.rhs
    if margin <= 0:
        return None
    return point, margin


def _find_relative_interior_point(
    rows: tuple[IntegralConstraint, ...], particle_number: int, d: int
) -> Point:
    augmented = np.asarray([[*row.coeffs, 1] for row in rows], dtype=float)
    result = linprog(
        np.asarray([*[0] * d, -1], dtype=float),
        A_ub=augmented,
        b_ub=np.asarray([row.rhs for row in rows], dtype=float),
        A_eq=np.asarray([[*[1] * d, 0]], dtype=float),
        b_eq=np.asarray([particle_number], dtype=float),
        bounds=[(None, None)] * d + [(0, 1)],
        method="highs-ds",
    )
    if not result.success or result.x[-1] <= 0:
        raise CertificateSearchError("the final system has no discovered strict point")

    for denominator_bound in (10**3, 10**6, 10**9, 10**12):
        coordinates = [
            Fraction(float(value)).limit_denominator(denominator_bound)
            for value in result.x[: d - 1]
        ]
        coordinates.append(Fraction(particle_number) - sum(coordinates))
        point = tuple(coordinates)
        if all(_strictly_satisfies(row, point) for row in rows):
            return point
    raise CertificateSearchError("could not reconstruct a rational strict point")


def _verify_farkas_step(
    step: FarkasRemovalCertificate,
    active: tuple[int, ...],
    certificate: OrderedSliceReductionCertificate,
) -> bool:
    if step.target_index not in active:
        return False
    target = certificate.candidate_rows[step.target_index]
    structural_indices = [index for index, _value in step.structural_terms]
    candidate_indices = [index for index, _value in step.candidate_terms]
    if len(structural_indices) != len(set(structural_indices)):
        return False
    if len(candidate_indices) != len(set(candidate_indices)):
        return False
    if any(not 0 <= index < len(certificate.structural_rows) for index in structural_indices):
        return False
    if any(index not in active or index == step.target_index for index in candidate_indices):
        return False
    if any(value <= 0 for _index, value in step.structural_terms):
        return False
    if any(value <= 0 for _index, value in step.candidate_terms):
        return False
    if step.rhs_slack < 0:
        return False

    d = certificate.d
    coefficients = [step.trace_multiplier for _index in range(d)]
    rhs = step.trace_multiplier * certificate.particle_number + step.rhs_slack
    for index, multiplier in step.structural_terms:
        row = certificate.structural_rows[index]
        coefficients = [
            value + multiplier * coefficient
            for value, coefficient in zip(coefficients, row.coeffs, strict=True)
        ]
        rhs += multiplier * row.rhs
    for index, multiplier in step.candidate_terms:
        row = certificate.candidate_rows[index]
        coefficients = [
            value + multiplier * coefficient
            for value, coefficient in zip(coefficients, row.coeffs, strict=True)
        ]
        rhs += multiplier * row.rhs
    return tuple(coefficients) == target.coeffs and rhs == target.rhs


def verify_ordered_slice_reduction_certificate(
    certificate: OrderedSliceReductionCertificate,
) -> bool:
    """Replay every implication, strict point, and facet witness over ``Q``."""
    try:
        if not 0 < certificate.particle_number < certificate.d:
            return False
        if certificate.structural_rows != ordered_pauli_inequalities(certificate.d):
            return False
        _validate_rows(certificate.candidate_rows, certificate.d)

        active = tuple(range(len(certificate.candidate_rows)))
        for step in certificate.removals:
            if not _verify_farkas_step(step, active, certificate):
                return False
            active = tuple(index for index in active if index != step.target_index)
        if active != certificate.retained_indices:
            return False

        final_rows = certificate.structural_rows + certificate.retained_rows
        interior = certificate.relative_interior_point
        if len(interior) != certificate.d or sum(interior) != certificate.particle_number:
            return False
        if not all(_strictly_satisfies(row, interior) for row in final_rows):
            return False

        facet_indices = tuple(facet.target_index for facet in certificate.facets)
        if facet_indices != certificate.retained_indices:
            return False
        if len(facet_indices) != len(set(facet_indices)):
            return False
        retained_set = set(certificate.retained_indices)
        for facet in certificate.facets:
            point = facet.violating_point
            if len(point) != certificate.d or sum(point) != certificate.particle_number:
                return False
            target = certificate.candidate_rows[facet.target_index]
            premises = certificate.structural_rows + tuple(
                certificate.candidate_rows[index]
                for index in certificate.retained_indices
                if index != facet.target_index
            )
            if facet.target_index not in retained_set:
                return False
            if not all(_satisfies(row, point) for row in premises):
                return False
            margin = _lhs(target, point) - target.rhs
            if margin <= 0 or margin != facet.violation_margin:
                return False
        return True
    except (IndexError, TypeError, ValueError, ZeroDivisionError):
        return False


def certify_ordered_slice_redundancy(
    candidate_rows: Sequence[IntegralConstraint],
    d: int,
    *,
    particle_number: int = 3,
) -> OrderedSliceReductionCertificate:
    """Eliminate redundant rows and prove every retained row is a facet.

    Rows are considered in input order.  A redundant row is removed
    immediately, and its Farkas certificate may use any row still active at
    that step.  A row with an exact violating witness remains nonredundant if
    later rows are removed because those removals only enlarge the comparison
    polyhedron.
    """
    if type(particle_number) is not int or not 0 < particle_number < d:
        raise ValueError("require an integer particle number with 0 < N < d")
    structural = ordered_pauli_inequalities(d)
    candidates = _validate_rows(candidate_rows, d)
    active = list(range(len(candidates)))
    removals: list[FarkasRemovalCertificate] = []
    facets: dict[int, FacetWitnessCertificate] = {}

    for target_index in range(len(candidates)):
        if target_index not in active:
            continue
        other_indices = tuple(index for index in active if index != target_index)
        removal = _find_farkas_removal(
            target_index,
            candidates[target_index],
            structural,
            candidates,
            other_indices,
            particle_number,
        )
        if removal is not None:
            removals.append(removal)
            active.remove(target_index)
            continue

        premises = structural + tuple(candidates[index] for index in other_indices)
        violating = _find_violating_point(candidates[target_index], premises, particle_number)
        if violating is None:
            raise CertificateSearchError(
                f"could not certify candidate row {target_index} as redundant or facet"
            )
        point, margin = violating
        facets[target_index] = FacetWitnessCertificate(
            target_index=target_index,
            violating_point=point,
            violation_margin=margin,
        )

    retained = tuple(active)
    final_rows = structural + tuple(candidates[index] for index in retained)
    interior = _find_relative_interior_point(final_rows, particle_number, d)
    result = OrderedSliceReductionCertificate(
        particle_number=particle_number,
        d=d,
        candidate_rows=candidates,
        structural_rows=structural,
        removals=tuple(removals),
        retained_indices=retained,
        relative_interior_point=interior,
        facets=tuple(facets[index] for index in retained),
    )
    if not verify_ordered_slice_reduction_certificate(result):
        raise AssertionError("constructed redundancy certificate failed exact replay")
    return result


__all__ = [
    "CertificateSearchError",
    "FacetWitnessCertificate",
    "FarkasRemovalCertificate",
    "OrderedSliceReductionCertificate",
    "certify_ordered_slice_redundancy",
    "ordered_pauli_inequalities",
    "verify_ordered_slice_reduction_certificate",
]
