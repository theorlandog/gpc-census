#!/usr/bin/env python3
"""Exact observable ranges for complex Hermitian triangles on a fixed 1-RDM fiber.

The fixed Borland--Dennis fiber has carrier amplitudes

    (sqrt(3/5), exp(i alpha)/2, sqrt(3/20) exp(i beta)).

For a rational Gaussian Hermitian 3x3 carrier matrix, the two-angle constrained
search reduces to a quartic plane curve with rational coefficients.  Its
critical-value discriminant is a degree-six polynomial over Q.  Under the
explicit nondegeneracy checks in :func:`solve_complex_triangle`, the smallest
and largest real roots of that sextic are the exact global range endpoints.

No floating-point optimizer is used to select an endpoint.  Decimal values are
rendered only after exact root isolation.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "data" / "complex_fixed_1rdm_observable_ranges.json"

WEIGHTS = (sp.Rational(3, 5), sp.Rational(1, 4), sp.Rational(3, 20))
SQRT15 = sp.sqrt(15)
EDGE_SCALES = (SQRT15 / 5, sp.Rational(3, 5), SQRT15 / 10)

D0 = (0, 1, 2)
D1 = (0, 3, 4)
D2 = (1, 3, 5)
CARRIER = (D0, D1, D2)

z = sp.symbols("z", real=True)
y = sp.symbols("y", real=True)
u = sp.symbols("u", real=True)


class ExactSolverError(RuntimeError):
    """The exact generic certificate conditions do not hold for an input."""


@dataclasses.dataclass(frozen=True)
class TriangleInvariants:
    """Gauge-invariant rational data for the branch-eliminated search."""

    offset: sp.Rational
    amplitude_norm: sp.Rational
    branch_direction_norm: sp.Rational
    mixed_product: sp.Rational
    branch_offset: sp.Rational
    signed_cycle_area: sp.Rational

    @property
    def gram_determinant(self) -> sp.Rational:
        return sp.factor(
            self.amplitude_norm * self.branch_direction_norm
            - self.mixed_product**2
        )


@dataclasses.dataclass(frozen=True)
class RootInterval:
    lower: sp.Rational
    upper: sp.Rational
    multiplicity: int

    def as_json(self) -> dict[str, object]:
        midpoint = (self.lower + self.upper) / 2
        return {
            "lower": rational_string(self.lower),
            "upper": rational_string(self.upper),
            "multiplicity": self.multiplicity,
            "midpoint_decimal": decimal_string(midpoint, 18),
        }


@dataclasses.dataclass(frozen=True)
class ExactRangeResult:
    matrix: sp.Matrix
    invariants: TriangleInvariants
    centered_critical_polynomial: sp.Poly
    energy_critical_polynomial: sp.Poly
    lift_denominator: sp.Poly
    lift_numerator: sp.Poly
    root_intervals: tuple[RootInterval, ...]

    @property
    def minimum(self) -> RootInterval:
        return self.root_intervals[0]

    @property
    def maximum(self) -> RootInterval:
        return self.root_intervals[-1]


def _as_exact_gaussian(value: sp.Expr) -> tuple[sp.Rational, sp.Rational]:
    expression = sp.sympify(value)
    real, imag = (sp.simplify(part) for part in expression.as_real_imag())
    if not real.is_Rational or not imag.is_Rational:
        raise ValueError(f"expected a rational Gaussian number, got {value!r}")
    return sp.Rational(real), sp.Rational(imag)


def validate_hermitian_matrix(matrix: Sequence[Sequence[sp.Expr]]) -> sp.Matrix:
    result = sp.Matrix(matrix)
    if result.shape != (3, 3):
        raise ValueError(f"expected a 3x3 matrix, got {result.shape}")
    for row in range(3):
        for column in range(3):
            _as_exact_gaussian(result[row, column])
            if sp.simplify(result[row, column] - sp.conjugate(result[column, row])):
                raise ValueError("matrix is not Hermitian")
        if _as_exact_gaussian(result[row, row])[1] != 0:
            raise ValueError("Hermitian diagonal entries must be real")
    return result


def gaussian_abs_squared(value: sp.Expr) -> sp.Rational:
    real, imag = _as_exact_gaussian(value)
    return real**2 + imag**2


def cycle_product(matrix: sp.Matrix) -> sp.Expr:
    """Return h01 h12 conjugate(h02), invariant under carrier rephasing."""
    return sp.expand_complex(matrix[0, 1] * matrix[1, 2] * sp.conjugate(matrix[0, 2]))


def triangle_invariants(matrix: Sequence[Sequence[sp.Expr]]) -> TriangleInvariants:
    h = validate_hermitian_matrix(matrix)
    h01, h02, h12 = h[0, 1], h[0, 2], h[1, 2]
    cycle_real, cycle_imag = _as_exact_gaussian(cycle_product(h))

    offset = sp.factor(sum(WEIGHTS[index] * h[index, index] for index in range(3)))
    if not offset.is_Rational:
        raise AssertionError("rational Hermitian diagonal produced a nonrational offset")

    # If a = 2 sqrt(w0 w1) h01, b = 2 sqrt(w0 w2) h02, and
    # c = 2 sqrt(w1 w2) h12, then after eliminating beta
    #
    #   E = E0 + a_vec . (cos alpha, sin alpha)
    #           +/- sqrt(R + d_vec . (cos alpha, sin alpha)).
    #
    # The following are |a_vec|^2, |d_vec|^2, a_vec.d_vec, and R.
    amplitude_norm = sp.Rational(3, 5) * gaussian_abs_squared(h01)
    branch_direction_norm = (
        sp.Rational(27, 125)
        * gaussian_abs_squared(h02)
        * gaussian_abs_squared(h12)
    )
    mixed_product = sp.Rational(9, 25) * cycle_real
    branch_offset = (
        sp.Rational(9, 25) * gaussian_abs_squared(h02)
        + sp.Rational(3, 20) * gaussian_abs_squared(h12)
    )
    signed_cycle_area = sp.Rational(9, 25) * cycle_imag

    result = TriangleInvariants(
        offset=sp.Rational(offset),
        amplitude_norm=sp.factor(amplitude_norm),
        branch_direction_norm=sp.factor(branch_direction_norm),
        mixed_product=sp.factor(mixed_product),
        branch_offset=sp.factor(branch_offset),
        signed_cycle_area=sp.factor(signed_cycle_area),
    )
    if sp.factor(result.gram_determinant - result.signed_cycle_area**2) != 0:
        raise AssertionError("cycle-area identity failed")
    return result


def branch_curve(invariants: TriangleInvariants) -> sp.Expr:
    """Return the exact centered-energy curve F(z,y)=0.

    Here y is the signed result of optimizing the beta phasor at fixed alpha,
    and z is E-E0.  If the cycle area is nonzero, this equation is equivalent
    to the original two-angle search rather than merely necessary.
    """
    a = invariants.amplitude_norm
    d = invariants.branch_direction_norm
    c = invariants.mixed_product
    r = invariants.branch_offset
    delta_squared = invariants.gram_determinant
    return sp.expand(
        d * (z - y) ** 2
        - 2 * c * (z - y) * (y**2 - r)
        + a * (y**2 - r) ** 2
        - delta_squared
    )


def _primitive_integer_poly(expression: sp.Expr, variable: sp.Symbol) -> sp.Poly:
    polynomial = sp.Poly(sp.expand(expression), variable, domain=sp.QQ)
    _, cleared = polynomial.clear_denoms(convert=True)
    coefficients = [int(coefficient) for coefficient in cleared.all_coeffs()]
    divisor = 0
    for coefficient in coefficients:
        divisor = math.gcd(divisor, abs(coefficient))
    divisor = max(divisor, 1)
    coefficients = [coefficient // divisor for coefficient in coefficients]
    if coefficients[0] < 0:
        coefficients = [-coefficient for coefficient in coefficients]
    degree = len(coefficients) - 1
    expression = sum(
        coefficient * variable ** (degree - index)
        for index, coefficient in enumerate(coefficients)
    )
    return sp.Poly(expression, variable, domain=sp.ZZ)


def critical_polynomial(invariants: TriangleInvariants) -> sp.Poly:
    """Compute the primitive sextic discriminant in centered energy."""
    curve = branch_curve(invariants)
    resultant = sp.resultant(curve, sp.diff(curve, y), y)
    polynomial = _primitive_integer_poly(resultant, z)
    if polynomial.degree() != 6:
        raise ExactSolverError(
            f"expected a degree-six critical polynomial, got degree {polynomial.degree()}"
        )
    return polynomial


def _normalize_polynomial_pair(
    first: sp.Expr, second: sp.Expr, variable: sp.Symbol
) -> tuple[sp.Poly, sp.Poly]:
    first_poly = sp.Poly(first, variable, domain=sp.QQ)
    second_poly = sp.Poly(second, variable, domain=sp.QQ)
    coefficients = first_poly.all_coeffs() + second_poly.all_coeffs()
    denominator = 1
    for coefficient in coefficients:
        denominator = sp.ilcm(denominator, int(sp.denom(coefficient)))
    integers = [int(coefficient * denominator) for coefficient in coefficients]
    divisor = 0
    for coefficient in integers:
        divisor = math.gcd(divisor, abs(coefficient))
    scale = sp.Rational(denominator, max(divisor, 1))
    normalized_first = sp.Poly(sp.expand(first * scale), variable, domain=sp.ZZ)
    normalized_second = sp.Poly(sp.expand(second * scale), variable, domain=sp.ZZ)
    leading = next(
        coefficient
        for coefficient in normalized_first.all_coeffs() + normalized_second.all_coeffs()
        if coefficient
    )
    if leading < 0:
        normalized_first = -normalized_first
        normalized_second = -normalized_second
    return normalized_first, normalized_second


def linear_lift_certificate(
    invariants: TriangleInvariants, polynomial: sp.Poly
) -> tuple[sp.Poly, sp.Poly]:
    """Return A(z), B(z) certifying y=-B(z)/A(z) at every sextic root."""
    curve = branch_curve(invariants)
    sequence = sp.subresultants(curve, sp.diff(curve, y), y)
    linear = next(
        (candidate for candidate in reversed(sequence) if sp.degree(candidate, y) == 1),
        None,
    )
    if linear is None:
        raise ExactSolverError("the subresultant sequence has no degree-one lift")
    as_poly = sp.Poly(linear, y)
    denominator = as_poly.coeff_monomial(y)
    numerator = as_poly.coeff_monomial(1)
    denominator_poly, numerator_poly = _normalize_polynomial_pair(
        denominator, numerator, z
    )

    if sp.gcd(denominator_poly, polynomial).degree() != 0:
        raise ExactSolverError(
            "the critical sextic and linear-lift denominator are not coprime"
        )

    rational_lift = -numerator_poly.as_expr() / denominator_poly.as_expr()
    curve_numerator = sp.together(curve.subs(y, rational_lift)).as_numer_denom()[0]
    derivative_numerator = sp.together(
        sp.diff(curve, y).subs(y, rational_lift)
    ).as_numer_denom()[0]
    if sp.rem(curve_numerator, polynomial.as_expr(), z) != 0:
        raise AssertionError("linear lift does not lie on the branch curve modulo sextic")
    if sp.rem(derivative_numerator, polynomial.as_expr(), z) != 0:
        raise AssertionError("linear lift is not critical modulo sextic")
    return denominator_poly, numerator_poly


def _refine_interval_to_decimal(
    polynomial: sp.Poly,
    lower: sp.Rational,
    upper: sp.Rational,
    digits: int,
) -> tuple[sp.Rational, sp.Rational]:
    scale = 10**digits
    lower_scaled = math.floor(sp.Rational(lower) * scale) - 1
    upper_scaled = math.ceil(sp.Rational(upper) * scale) + 1
    candidate_lower = sp.Rational(lower_scaled, scale)
    candidate_upper = sp.Rational(upper_scaled, scale)
    if polynomial.count_roots(candidate_lower, candidate_upper) != 1:
        raise AssertionError("decimal enclosure did not preserve one-root isolation")
    return candidate_lower, candidate_upper


def isolate_real_roots(polynomial: sp.Poly, digits: int = 14) -> tuple[RootInterval, ...]:
    exact_intervals = polynomial.intervals(eps=sp.Rational(1, 10 ** (digits + 4)))
    roots: list[RootInterval] = []
    for (lower, upper), multiplicity in exact_intervals:
        decimal_lower, decimal_upper = _refine_interval_to_decimal(
            polynomial, sp.Rational(lower), sp.Rational(upper), digits
        )
        roots.append(
            RootInterval(
                lower=decimal_lower,
                upper=decimal_upper,
                multiplicity=int(multiplicity),
            )
        )
    if not roots:
        raise ExactSolverError("critical polynomial has no real roots")
    return tuple(roots)


def solve_complex_triangle(
    matrix: Sequence[Sequence[sp.Expr]], *, isolation_digits: int = 14
) -> ExactRangeResult:
    """Solve the exact fixed-fiber range for a generic nonzero-flux triangle.

    The accepted class is exact and explicit: rational Gaussian Hermitian input,
    nonzero gauge-invariant cycle imaginary part, squarefree critical sextic, and
    a coprime linear subresultant lift.  These conditions are checked rather than
    assumed.  Degenerate zero/pi-flux cases are covered by the real-symmetric
    solver and by the unit-family bridge emitted in the artifact.
    """
    h = validate_hermitian_matrix(matrix)
    invariants = triangle_invariants(h)
    if invariants.signed_cycle_area == 0:
        raise ExactSolverError(
            "cycle flux is 0 or pi; use the real-symmetric branch solver"
        )
    if invariants.gram_determinant <= 0:
        raise AssertionError("nonzero cycle area must give a positive Gram determinant")

    centered = critical_polynomial(invariants)
    if sp.gcd(centered, centered.diff()).degree() != 0:
        raise ExactSolverError("critical sextic is not squarefree")
    denominator, numerator = linear_lift_certificate(invariants, centered)
    roots = isolate_real_roots(centered, digits=isolation_digits)

    # Translate from z=E-E0 to the physical energy variable E.
    energy_expression = sp.expand(centered.as_expr().subs(z, z - invariants.offset))
    energy_polynomial = _primitive_integer_poly(energy_expression, z)
    shifted_roots = tuple(
        RootInterval(
            lower=root.lower + invariants.offset,
            upper=root.upper + invariants.offset,
            multiplicity=root.multiplicity,
        )
        for root in roots
    )
    return ExactRangeResult(
        matrix=h,
        invariants=invariants,
        centered_critical_polynomial=centered,
        energy_critical_polynomial=energy_polynomial,
        lift_denominator=denominator,
        lift_numerator=numerator,
        root_intervals=shifted_roots,
    )


def unit_flux_critical_polynomial(parameter: sp.Expr = u) -> sp.Poly:
    """Critical sextic for unit edge magnitudes and cos(Phi)=parameter."""
    parameter = sp.sympify(parameter)
    expression = (
        1_000_000 * z**6
        + 2_000_000 * parameter * z**5
        + (1_110_000 * parameter**2 - 2_330_000) * z**4
        + (180_000 * parameter**3 - 6_060_000 * parameter) * z**3
        + (-5_164_200 * parameter**2 - 143_700) * z**2
        + (-1_798_200 * parameter**3 - 57_600 * parameter) * z
        - 218_700 * parameter**4
        + 6_831 * parameter**2
        + 3_969
    )
    free_symbols = sorted(expression.free_symbols - {z}, key=str)
    domain = sp.QQ if not free_symbols else sp.QQ.frac_field(*free_symbols)
    return sp.Poly(sp.expand(expression), z, domain=domain)


def unit_flux_matrix(cosine: sp.Rational, sine: sp.Rational) -> sp.Matrix:
    cosine = sp.Rational(cosine)
    sine = sp.Rational(sine)
    if cosine**2 + sine**2 != 1:
        raise ValueError("cosine and sine must lie on the rational unit circle")
    phase = cosine + sp.I * sine
    return sp.Matrix([[0, 1, 1], [1, 0, phase], [1, sp.conjugate(phase), 0]])


def rational_string(value: sp.Expr) -> str:
    value = sp.Rational(value)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


def decimal_string(value: sp.Expr, digits: int = 16) -> str:
    return str(sp.N(value, digits))


def gaussian_json(value: sp.Expr) -> dict[str, str]:
    real, imag = _as_exact_gaussian(value)
    return {"re": rational_string(real), "im": rational_string(imag)}


def matrix_json(matrix: sp.Matrix) -> list[list[dict[str, str]]]:
    return [
        [gaussian_json(matrix[row, column]) for column in range(3)]
        for row in range(3)
    ]


def polynomial_coefficients(polynomial: sp.Poly) -> list[str]:
    return [str(coefficient) for coefficient in polynomial.all_coeffs()]


def invariant_json(invariants: TriangleInvariants) -> dict[str, str]:
    return {
        "E0": rational_string(invariants.offset),
        "A": rational_string(invariants.amplitude_norm),
        "D": rational_string(invariants.branch_direction_norm),
        "C": rational_string(invariants.mixed_product),
        "R": rational_string(invariants.branch_offset),
        "delta": rational_string(invariants.signed_cycle_area),
        "A_D_minus_C_squared": rational_string(invariants.gram_determinant),
    }


def result_json(name: str, result: ExactRangeResult) -> dict[str, object]:
    cycle = cycle_product(result.matrix)
    roots = [root.as_json() for root in result.root_intervals]
    return {
        "name": name,
        "matrix": matrix_json(result.matrix),
        "cycle_product": gaussian_json(cycle),
        "invariants": invariant_json(result.invariants),
        "centered_critical_polynomial_coefficients_descending": polynomial_coefficients(
            result.centered_critical_polynomial
        ),
        "energy_critical_polynomial_coefficients_descending": polynomial_coefficients(
            result.energy_critical_polynomial
        ),
        "linear_lift": {
            "meaning": "At every sextic root z, y=-B(z)/A(z).",
            "A_coefficients_descending": polynomial_coefficients(result.lift_denominator),
            "B_coefficients_descending": polynomial_coefficients(result.lift_numerator),
            "coprime_to_critical_polynomial": True,
        },
        "real_roots_in_increasing_order": roots,
        "exact_range": {
            "minimum_root_index": 0,
            "minimum_isolating_interval": roots[0],
            "maximum_root_index": len(roots) - 1,
            "maximum_isolating_interval": roots[-1],
            "selection_method": (
                "smallest and largest real roots after exact Sturm isolation; "
                "the coprime linear lift proves every real root is physical"
            ),
        },
    }


def physical_operator_realization() -> dict[str, object]:
    return {
        "carrier_determinants": [list(determinant) for determinant in CARRIER],
        "pair_annihilator_convention": "A_ij=a_j*a_i for i<j",
        "fermionic_channel_signs": {
            "D0_to_D1_via_A34_dagger_A12": 1,
            "D0_to_D2_via_A35_dagger_A02": -1,
            "D1_to_D2_via_A15_dagger_A04": 1,
        },
        "formula": (
            "h00*n0*n1+h11*n0*n3+h22*n1*n3"
            "+conj(h01)*A34^dagger*A12+h01*A12^dagger*A34"
            "-conj(h02)*A35^dagger*A02-h02*A02^dagger*A35"
            "+conj(h12)*A15^dagger*A04+h12*A04^dagger*A15"
        ),
        "scope": (
            "This at-most-two-body operator preserves the three-determinant carrier "
            "and restricts to the requested Hermitian matrix."
        ),
    }


def build_artifact() -> dict[str, object]:
    quarter_flux = solve_complex_triangle(unit_flux_matrix(sp.Rational(0), sp.Rational(1)))
    pythagorean_flux = solve_complex_triangle(
        unit_flux_matrix(sp.Rational(3, 5), sp.Rational(4, 5))
    )

    plus_factorization = sp.factor(
        unit_flux_critical_polynomial(sp.Rational(1)).as_expr(), extension=SQRT15
    )
    pi_factorization = sp.factor(
        unit_flux_critical_polynomial(sp.Rational(-1)).as_expr(), extension=SQRT15
    )

    symbolic_coefficients = [
        sp.sstr(coefficient)
        for coefficient in unit_flux_critical_polynomial(u).all_coeffs()
    ]

    return {
        "schema_version": 1,
        "scope": (
            "Exact constrained-search ranges on the complete pinned Borland--Dennis "
            "fiber for rational Gaussian 3x3 Hermitian carrier Hamiltonians with "
            "nonzero gauge-invariant cycle flux and a certified generic linear lift."
        ),
        "target_1rdm_diagonal": ["17/20", "3/4", "3/5", "2/5", "1/4", "3/20"],
        "fiber_state": (
            "sqrt(3/5)|012>+(exp(i*alpha)/2)|034>+"
            "sqrt(3/20)*exp(i*beta)|135>"
        ),
        "theorem": {
            "phasors": [
                "a=(sqrt(15)/5)*h01",
                "b=(3/5)*h02",
                "c=(sqrt(15)/10)*h12",
            ],
            "fixed_alpha_elimination": (
                "E_sigma(alpha)=E0+Re(a*exp(i*alpha))"
                "+sigma*abs(b+c*exp(-i*alpha)), sigma in {-1,+1}"
            ),
            "cycle_invariant": "chi=h01*h12*conj(h02)",
            "rational_invariants": {
                "A": "(3/5)*abs(h01)^2",
                "D": "(27/125)*abs(h02)^2*abs(h12)^2",
                "C": "(9/25)*Re(chi)",
                "R": "(9/25)*abs(h02)^2+(3/20)*abs(h12)^2",
                "delta": "(9/25)*Im(chi)",
                "identity": "A*D-C^2=delta^2",
            },
            "branch_curve": (
                "F(z,y)=D(z-y)^2-2C(z-y)(y^2-R)+A(y^2-R)^2-(A*D-C^2)"
            ),
            "critical_equation": "Res_y(F,dF/dy)=constant*P_H(z)",
            "phase_recovery": [
                "U=z-y, V=y^2-R",
                "cos(alpha)=(d2*U-p2*V)/delta",
                "sin(alpha)=(-d1*U+p1*V)/delta",
                "exp(i*beta)=conj(b+c*exp(-i*alpha))/y",
            ],
            "global_selection": (
                "When delta!=0, P_H is squarefree, and the emitted linear lift is "
                "coprime, the exact interval is E0 plus the smallest and largest "
                "real roots of P_H."
            ),
        },
        "physical_operator_realization": physical_operator_realization(),
        "unit_edge_flux_family": {
            "matrix": "[[0,1,1],[1,0,exp(i*Phi)],[1,exp(-i*Phi),0]]",
            "parameter": "u=cos(Phi)",
            "critical_polynomial_variable": "z=E",
            "coefficients_descending_in_z": symbolic_coefficients,
            "symmetry": "P_{-u}(-z)=P_u(z)",
            "zero_flux_factorization": sp.sstr(plus_factorization),
            "pi_flux_factorization": sp.sstr(pi_factorization),
            "zero_flux_range": ["-1", "3*(2+sqrt(15))/10"],
            "pi_flux_range": ["-3*(2+sqrt(15))/10", "1"],
            "interpretation": (
                "The sextic continuously bridges the unfrustrated and pi-flux "
                "real-symmetric regressions.  The range is invariant under Phi -> -Phi."
            ),
        },
        "regressions": [
            result_json("quarter_flux", quarter_flux),
            result_json("pythagorean_flux_cos_3_5_sin_4_5", pythagorean_flux),
        ],
        "verification": {
            "generator": "scripts/complex_fixed_1rdm_observable_ranges.py",
            "standalone_verifier": (
                "scripts/verify_complex_fixed_1rdm_observable_ranges_standalone.py"
            ),
            "arithmetic": (
                "rational resultants, polynomial gcd/divisibility, and exact Sturm "
                "root counts; decimal midpoints never select an endpoint"
            ),
        },
        "limitations": [
            "The exact generic engine currently requires rational Gaussian matrix entries.",
            "Zero and pi cycle flux use the real-symmetric solver or the exact unit-family limits.",
            (
                "A nongeneric repeated sextic or failed coprime lift is reported, "
                "not guessed through numerics."
            ),
            "This benchmark solves one fixed three-determinant fiber, not an arbitrary GPC fiber.",
        ],
    }


def render_artifact() -> str:
    return json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"


def update_checksum_ledger(output: Path) -> None:
    ledger = output.parent / "SHA256SUMS"
    if not ledger.exists():
        return
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    target = f"{digest}  ./complex_fixed_1rdm_observable_ranges.json"
    lines = [
        line
        for line in ledger.read_text(encoding="utf-8").splitlines()
        if not line.endswith("  ./complex_fixed_1rdm_observable_ranges.json")
    ]
    lines.append(target)
    lines.sort(key=lambda line: line.split("  ", 1)[1] if "  " in line else line)
    ledger.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(list(argv) if argv is not None else None)

    rendered = render_artifact()
    if args.check:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            print(f"{args.output} is stale")
            return 1
        print(f"{args.output} is current")
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    update_checksum_ledger(args.output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
