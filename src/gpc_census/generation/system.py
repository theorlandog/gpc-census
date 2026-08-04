"""Hybrid exact assembly of certified GPC systems.

The first completed gate is ``(N,d)=(3,6)``.  Valid outer inequalities come
from Ressayre certificates.  Attainable inner points come independently from
plethysm.  Exact vertex enumeration closes the sandwich when every outer
vertex is one of the certified inner points.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Sequence
from fractions import Fraction

import sympy as sp

from gpc_census.plethysm import inner_spectra

from .model import GeneratedSystem, IntegralConstraint, RessayreCertificate
from .ressayre import enumerate_ressayre, verify_ressayre_certificate
from .schubert import (
    borland_dennis_schubert_certificates,
    verify_borland_dennis_schubert_certificates,
)


def _primitive_constraint(
    coeffs: Sequence[int], rhs: int, *, equality: bool
) -> IntegralConstraint:
    values = list(coeffs) + [rhs]
    divisor = 0
    for value in values:
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:
        raise ValueError("zero relation")
    values = [value // divisor for value in values]
    if equality and next(value for value in values if value) < 0:
        values = [-value for value in values]
    return IntegralConstraint(tuple(values[:-1]), values[-1])


def _canonical_equation_mod_trace(
    h: Sequence[int], z: int, n: int
) -> IntegralConstraint:
    """Choose a sparse primitive representative modulo ``sum(lambda)=n``."""
    candidates: list[IntegralConstraint] = []
    for shift in sorted({0, *(-value for value in h)}):
        candidates.append(
            _primitive_constraint(
                [value + shift for value in h], z + shift * n, equality=True
            )
        )

    def score(row: IntegralConstraint) -> tuple[object, ...]:
        values = row.coeffs + (row.rhs,)
        return (
            sum(value != 0 for value in row.coeffs),
            max(abs(value) for value in values),
            sum(abs(value) for value in values),
            row.coeffs,
            row.rhs,
        )

    return min(candidates, key=score)


def derive_equalities(
    certificates: Sequence[RessayreCertificate], n: int
) -> tuple[IntegralConstraint, ...]:
    """Derive equalities from oppositely oriented valid inequalities."""
    oriented = {(certificate.h, certificate.z) for certificate in certificates}
    seen: set[tuple[tuple[int, ...], int]] = set()
    equalities: set[IntegralConstraint] = set()
    for h, z in sorted(oriented):
        opposite = (tuple(-value for value in h), -z)
        if opposite not in oriented or (h, z) in seen:
            continue
        seen.add((h, z))
        seen.add(opposite)
        equalities.add(_canonical_equation_mod_trace(h, z, n))
    return tuple(
        sorted(
            equalities,
            key=lambda row: (
                next(index for index, value in enumerate(row.coeffs) if value),
                row.coeffs,
                row.rhs,
            ),
        )
    )


def _equality_reducer(
    equalities: Sequence[IntegralConstraint], d: int
) -> tuple[tuple[tuple[sp.Rational, ...], ...], tuple[int, ...]]:
    matrix = sp.Matrix(
        [list(row.coeffs) + [-row.rhs] for row in equalities]
    )
    reduced, pivots = matrix.rref()
    if any(pivot >= d for pivot in pivots):
        raise ValueError("inconsistent affine equalities")
    rows = tuple(
        tuple(sp.Rational(reduced[row, column]) for column in range(d + 1))
        for row in range(reduced.rows)
    )
    return rows, tuple(pivots)


def restrict_to_affine_hull(
    row: IntegralConstraint,
    equalities: Sequence[IntegralConstraint],
) -> IntegralConstraint | None:
    """Return a primitive ``<=`` representative modulo exact equalities."""
    d = len(row.coeffs)
    reduced_rows, pivots = _equality_reducer(equalities, d)
    vector = [sp.Rational(value) for value in row.coeffs] + [sp.Rational(-row.rhs)]
    for equation, pivot in zip(reduced_rows, pivots, strict=True):
        coefficient = vector[pivot]
        vector = [
            value - coefficient * equation_value
            for value, equation_value in zip(vector, equation, strict=True)
        ]
    if not any(vector):
        return None
    denominator = sp.ilcm(*(sp.denom(value) for value in vector))
    integers = [int(value * denominator) for value in vector]
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    integers = [value // divisor for value in integers]
    # The augmented convention is [coeffs, -rhs].  Positive scaling only is
    # allowed because this is an inequality.
    return IntegralConstraint(tuple(integers[:-1]), -integers[-1])


def structural_inequalities(
    d: int, equalities: Sequence[IntegralConstraint]
) -> tuple[IntegralConstraint, ...]:
    """Return Pauli, Weyl-chamber, and nonnegativity rows on the affine hull."""
    raw: list[IntegralConstraint] = [
        IntegralConstraint((1,) + (0,) * (d - 1), 1),
        IntegralConstraint((0,) * (d - 1) + (-1,), 0),
    ]
    for index in range(d - 1):
        coeffs = [0] * d
        coeffs[index] = -1
        coeffs[index + 1] = 1
        raw.append(IntegralConstraint(tuple(coeffs), 0))
    restricted = {
        result
        for row in raw
        if (result := restrict_to_affine_hull(row, equalities)) is not None
    }
    return tuple(sorted(restricted))


def _paired_hyperplanes(
    certificates: Sequence[RessayreCertificate],
) -> set[tuple[tuple[int, ...], int]]:
    oriented = {(certificate.h, certificate.z) for certificate in certificates}
    return {
        pair
        for pair in oriented
        if (tuple(-value for value in pair[0]), -pair[1]) in oriented
    }


def generated_outer_inequalities(
    certificates: Sequence[RessayreCertificate],
    equalities: Sequence[IntegralConstraint],
) -> tuple[IntegralConstraint, ...]:
    """Canonicalize all one-sided Ressayre inequalities on the affine hull."""
    paired = _paired_hyperplanes(certificates)
    rows: set[IntegralConstraint] = set()
    for certificate in certificates:
        if (certificate.h, certificate.z) in paired:
            continue
        # H.lambda >= z becomes (-H).lambda <= -z.
        raw = IntegralConstraint(tuple(-value for value in certificate.h), -certificate.z)
        restricted = restrict_to_affine_hull(raw, equalities)
        if restricted is not None:
            rows.add(restricted)
    return tuple(sorted(rows))


def _to_fraction(value: sp.Rational) -> Fraction:
    return Fraction(int(sp.numer(value)), int(sp.denom(value)))


def _affine_chart(
    equalities: Sequence[IntegralConstraint], d: int
) -> tuple[tuple[Fraction, ...], tuple[tuple[Fraction, ...], ...]]:
    """Return ``x = base + basis*y`` for an exact affine system."""
    augmented = sp.Matrix([list(row.coeffs) + [row.rhs] for row in equalities])
    reduced, pivots = augmented.rref()
    if any(pivot >= d for pivot in pivots):
        raise ValueError("inconsistent equalities")
    free = tuple(index for index in range(d) if index not in pivots)
    base = [Fraction(0) for _index in range(d)]
    basis = [[Fraction(0) for _index in free] for _coordinate in range(d)]
    for free_position, coordinate in enumerate(free):
        basis[coordinate][free_position] = Fraction(1)
    for row_index, pivot in enumerate(pivots):
        base[pivot] = _to_fraction(sp.Rational(reduced[row_index, d]))
        for free_position, coordinate in enumerate(free):
            basis[pivot][free_position] = _to_fraction(
                -sp.Rational(reduced[row_index, coordinate])
            )
    return tuple(base), tuple(tuple(row) for row in basis)


def enumerate_vertices(
    equalities: Sequence[IntegralConstraint],
    inequalities: Sequence[IntegralConstraint],
) -> tuple[tuple[Fraction, ...], ...]:
    """Enumerate bounded-polytope vertices by exact active-set intersections."""
    if not equalities and not inequalities:
        raise ValueError("empty system")
    d = len((equalities or inequalities)[0].coeffs)
    base, basis = _affine_chart(equalities, d)
    dimension = len(basis[0]) if basis else 0

    reduced_rows: list[tuple[tuple[Fraction, ...], Fraction]] = []
    for row in inequalities:
        normal = tuple(
            sum(Fraction(row.coeffs[index]) * basis[index][column] for index in range(d))
            for column in range(dimension)
        )
        rhs = Fraction(row.rhs) - sum(
            Fraction(row.coeffs[index]) * base[index] for index in range(d)
        )
        if any(normal):
            reduced_rows.append((normal, rhs))
        elif rhs < 0:
            raise ValueError("infeasible constant inequality")

    if dimension == 0:
        return (base,)
    vertices: set[tuple[Fraction, ...]] = set()
    for active in itertools.combinations(reduced_rows, dimension):
        matrix = sp.Matrix(
            [[sp.Rational(value.numerator, value.denominator) for value in row[0]] for row in active]
        )
        if matrix.det() == 0:
            continue
        rhs = sp.Matrix(
            [sp.Rational(row[1].numerator, row[1].denominator) for row in active]
        )
        solution = tuple(_to_fraction(value) for value in matrix.inv() * rhs)
        if any(
            sum(normal[index] * solution[index] for index in range(dimension)) > bound
            for normal, bound in reduced_rows
        ):
            continue
        vertex = tuple(
            base[index]
            + sum(basis[index][column] * solution[column] for column in range(dimension))
            for index in range(d)
        )
        vertices.add(vertex)
    return tuple(sorted(vertices, reverse=True))


def affine_dimension(points: Iterable[Sequence[Fraction]]) -> int:
    points = tuple(tuple(point) for point in points)
    if not points:
        return -1
    if len(points) == 1:
        return 0
    base = points[0]
    matrix = sp.Matrix(
        [
            [
                sp.Rational((value - base[index]).numerator, (value - base[index]).denominator)
                for index, value in enumerate(point)
            ]
            for point in points[1:]
        ]
    )
    return int(matrix.rank())


def _facet_rows(
    rows: Sequence[IntegralConstraint], vertices: Sequence[tuple[Fraction, ...]], dimension: int
) -> tuple[IntegralConstraint, ...]:
    facets: list[IntegralConstraint] = []
    for row in rows:
        tight = [
            vertex
            for vertex in vertices
            if sum(Fraction(c) * value for c, value in zip(row.coeffs, vertex, strict=True))
            == row.rhs
        ]
        if affine_dimension(tight) == dimension - 1:
            facets.append(row)
    return tuple(sorted(set(facets)))


def _satisfies(point: Sequence[Fraction], row: IntegralConstraint) -> bool:
    return sum(
        Fraction(coefficient) * value
        for coefficient, value in zip(row.coeffs, point, strict=True)
    ) <= row.rhs


def generate_3_6(*, inner_degree: int = 6) -> GeneratedSystem:
    """Generate and close the Borland-Dennis system from first principles.

    Stored constraints and stored vertices are not read.  Completeness is an
    exact sandwich: Ressayre gives a valid outer polytope, plethysm gives true
    inner points, and every exactly enumerated outer vertex occurs in that
    inner cloud.
    """
    n, d = 3, 6
    enumeration = enumerate_ressayre(n, d)
    if not all(
        verify_ressayre_certificate(n, d, certificate)
        for certificate in enumeration.certificates
    ):
        raise AssertionError("generated Ressayre certificate failed replay")

    equalities = derive_equalities(enumeration.certificates, n)
    if not equalities:
        raise AssertionError("Ressayre enumeration discovered no affine equalities")
    structural = structural_inequalities(d, equalities)
    resssayre_rows = generated_outer_inequalities(enumeration.certificates, equalities)
    outer_rows = tuple(sorted(set(structural) | set(resssayre_rows)))
    vertices = enumerate_vertices(equalities, outer_rows)

    inner = inner_spectra(n, d, inner_degree)
    inner_points = set(inner)
    if not all(_satisfies(point, row) for point in inner_points for row in outer_rows):
        raise AssertionError("certified inner point violates the generated outer system")
    closed = bool(vertices) and set(vertices) <= inner_points
    if not closed:
        raise AssertionError("inner/outer sandwich did not close")

    dimension = affine_dimension(vertices)
    structural_set = set(structural)
    nonstructural = tuple(row for row in resssayre_rows if row not in structural_set)
    inequalities = _facet_rows(nonstructural, vertices, dimension)
    schubert = borland_dennis_schubert_certificates()
    if not verify_borland_dennis_schubert_certificates(schubert):
        raise AssertionError("rank-6 Schubert certificates failed replay")
    if tuple(certificate.constraint for certificate in schubert[:3]) != equalities:
        raise AssertionError("Schubert pair rows and Ressayre equalities disagree")
    schubert_gpc = restrict_to_affine_hull(schubert[3].constraint, equalities)
    if schubert_gpc not in inequalities:
        raise AssertionError("Ressayre and Schubert rank-6 constraints disagree")
    return GeneratedSystem(
        n=n,
        d=d,
        inequalities=inequalities,
        equalities=equalities,
        structural_inequalities=structural,
        outer_vertices=vertices,
        inner_point_count=len(inner_points),
        outer_vertices_in_inner_cloud=closed,
        enumeration=enumeration,
        schubert_cross_checks=schubert,
    )
