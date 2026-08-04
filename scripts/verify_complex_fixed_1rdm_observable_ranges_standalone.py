#!/usr/bin/env python3
"""Standard-library verifier for complex fixed-1RDM observable ranges.

This script imports neither SymPy nor project code.  It reconstructs the exact
rational invariants, the quartic branch curve, its Sylvester resultant, the
linear critical-point lift, and every Sturm isolating interval in the generated
artifact.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT = ROOT / "results" / "data" / "complex_fixed_1rdm_observable_ranges.json"

Poly = tuple[Fraction, ...]  # ascending powers
Gaussian = tuple[Fraction, Fraction]


def frac(value: object) -> Fraction:
    return Fraction(str(value))


def trim(poly: Sequence[Fraction]) -> Poly:
    values = list(poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [Fraction(0)])


def zero() -> Poly:
    return (Fraction(0),)


def one() -> Poly:
    return (Fraction(1),)


def constant(value: Fraction | int) -> Poly:
    return (Fraction(value),)


def add(left: Poly, right: Poly) -> Poly:
    size = max(len(left), len(right))
    return trim(
        [
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(size)
        ]
    )


def neg(poly: Poly) -> Poly:
    return trim([-coefficient for coefficient in poly])


def sub(left: Poly, right: Poly) -> Poly:
    return add(left, neg(right))


def scale(poly: Poly, scalar: Fraction | int) -> Poly:
    scalar = Fraction(scalar)
    return trim([scalar * coefficient for coefficient in poly])


def mul(left: Poly, right: Poly) -> Poly:
    if left == zero() or right == zero():
        return zero()
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] += left_value * right_value
    return trim(result)


def power(poly: Poly, exponent: int) -> Poly:
    result = one()
    base = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = mul(result, base)
        base = mul(base, base)
        remaining >>= 1
    return result


def derivative(poly: Poly) -> Poly:
    if len(poly) <= 1:
        return zero()
    return trim([index * poly[index] for index in range(1, len(poly))])


def evaluate(poly: Poly, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def divmod_poly(dividend: Poly, divisor: Poly) -> tuple[Poly, Poly]:
    divisor = trim(divisor)
    if divisor == zero():
        raise ZeroDivisionError("polynomial division by zero")
    remainder = list(trim(dividend))
    quotient = [Fraction(0)] * max(1, len(remainder) - len(divisor) + 1)
    while len(remainder) >= len(divisor) and trim(remainder) != zero():
        degree = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[degree] += coefficient
        for index, value in enumerate(divisor):
            remainder[index + degree] -= coefficient * value
        remainder = list(trim(remainder))
    return trim(quotient), trim(remainder)


def monic(poly: Poly) -> Poly:
    poly = trim(poly)
    if poly == zero():
        return zero()
    return scale(poly, 1 / poly[-1])


def gcd_poly(left: Poly, right: Poly) -> Poly:
    left, right = trim(left), trim(right)
    while right != zero():
        _, remainder = divmod_poly(left, right)
        left, right = right, remainder
    return monic(left)


def primitive_integer(poly: Poly) -> tuple[int, ...]:
    poly = trim(poly)
    denominator_lcm = 1
    for coefficient in poly:
        denominator_lcm = math.lcm(denominator_lcm, coefficient.denominator)
    integers = [int(coefficient * denominator_lcm) for coefficient in poly]
    divisor = 0
    for coefficient in integers:
        divisor = math.gcd(divisor, abs(coefficient))
    divisor = max(divisor, 1)
    integers = [coefficient // divisor for coefficient in integers]
    if integers[-1] < 0:
        integers = [-coefficient for coefficient in integers]
    return tuple(integers)


def parse_poly_descending(values: Sequence[object]) -> Poly:
    return trim([frac(value) for value in reversed(values)])


def polynomial_sign_at_infinity(poly: Poly, positive: bool) -> int:
    leading = poly[-1]
    sign = 1 if leading > 0 else -1
    if not positive and (len(poly) - 1) % 2:
        sign = -sign
    return sign


def sign(value: Fraction) -> int:
    return (value > 0) - (value < 0)


def variations(signs: Iterable[int]) -> int:
    filtered = [value for value in signs if value]
    return sum(left != right for left, right in zip(filtered, filtered[1:]))


def sturm_sequence(poly: Poly) -> tuple[Poly, ...]:
    sequence = [trim(poly), derivative(poly)]
    while sequence[-1] != zero():
        _, remainder = divmod_poly(sequence[-2], sequence[-1])
        if remainder == zero():
            break
        sequence.append(neg(remainder))
    return tuple(sequence)


def sturm_variation_at(sequence: Sequence[Poly], value: Fraction) -> int:
    return variations(sign(evaluate(poly, value)) for poly in sequence)


def roots_between(sequence: Sequence[Poly], lower: Fraction, upper: Fraction) -> int:
    if lower >= upper:
        raise ValueError("invalid interval")
    if any(evaluate(poly, lower) == 0 for poly in sequence[:1]):
        raise ValueError("lower interval endpoint is a root")
    if any(evaluate(poly, upper) == 0 for poly in sequence[:1]):
        raise ValueError("upper interval endpoint is a root")
    return sturm_variation_at(sequence, lower) - sturm_variation_at(sequence, upper)


def total_real_roots(sequence: Sequence[Poly]) -> int:
    negative = variations(
        polynomial_sign_at_infinity(poly, positive=False) for poly in sequence
    )
    positive = variations(
        polynomial_sign_at_infinity(poly, positive=True) for poly in sequence
    )
    return negative - positive


def permutation_sign(permutation: Sequence[int]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant(matrix: Sequence[Sequence[Poly]]) -> Poly:
    size = len(matrix)
    result = zero()
    for permutation in itertools.permutations(range(size)):
        term = one()
        for row, column in enumerate(permutation):
            term = mul(term, matrix[row][column])
            if term == zero():
                break
        if term != zero():
            result = add(result, scale(term, permutation_sign(permutation)))
    return trim(result)


def sylvester_resultant(first_desc: Sequence[Poly], second_desc: Sequence[Poly]) -> Poly:
    first_degree = len(first_desc) - 1
    second_degree = len(second_desc) - 1
    size = first_degree + second_degree
    matrix = [[zero() for _ in range(size)] for _ in range(size)]
    for row in range(second_degree):
        for index, coefficient in enumerate(first_desc):
            matrix[row][row + index] = coefficient
    for local_row in range(first_degree):
        row = second_degree + local_row
        for index, coefficient in enumerate(second_desc):
            matrix[row][local_row + index] = coefficient
    return determinant(matrix)


def gaussian(value: dict[str, object]) -> Gaussian:
    return frac(value["re"]), frac(value["im"])


def g_add(left: Gaussian, right: Gaussian) -> Gaussian:
    return left[0] + right[0], left[1] + right[1]


def g_mul(left: Gaussian, right: Gaussian) -> Gaussian:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def g_conj(value: Gaussian) -> Gaussian:
    return value[0], -value[1]


def g_abs2(value: Gaussian) -> Fraction:
    return value[0] ** 2 + value[1] ** 2


def parse_matrix(record: dict[str, object]) -> tuple[tuple[Gaussian, ...], ...]:
    matrix = tuple(tuple(gaussian(entry) for entry in row) for row in record["matrix"])
    if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
        raise AssertionError("matrix is not 3x3")
    for row in range(3):
        if matrix[row][row][1] != 0:
            raise AssertionError("diagonal entry is not real")
        for column in range(3):
            if matrix[row][column] != g_conj(matrix[column][row]):
                raise AssertionError("matrix is not Hermitian")
    return matrix


def exact_invariants(matrix: tuple[tuple[Gaussian, ...], ...]) -> dict[str, Fraction]:
    cycle = g_mul(g_mul(matrix[0][1], matrix[1][2]), g_conj(matrix[0][2]))
    offset = (
        Fraction(3, 5) * matrix[0][0][0]
        + Fraction(1, 4) * matrix[1][1][0]
        + Fraction(3, 20) * matrix[2][2][0]
    )
    amplitude_norm = Fraction(3, 5) * g_abs2(matrix[0][1])
    direction_norm = (
        Fraction(27, 125) * g_abs2(matrix[0][2]) * g_abs2(matrix[1][2])
    )
    mixed = Fraction(9, 25) * cycle[0]
    branch_offset = (
        Fraction(9, 25) * g_abs2(matrix[0][2])
        + Fraction(3, 20) * g_abs2(matrix[1][2])
    )
    delta = Fraction(9, 25) * cycle[1]
    if amplitude_norm * direction_norm - mixed**2 != delta**2:
        raise AssertionError("A D - C^2 != delta^2")
    return {
        "E0": offset,
        "A": amplitude_norm,
        "D": direction_norm,
        "C": mixed,
        "R": branch_offset,
        "delta": delta,
        "A_D_minus_C_squared": delta**2,
    }


def branch_curve_coefficients(invariants: dict[str, Fraction]) -> list[Poly]:
    """Return F coefficients in ascending powers of y."""
    a = invariants["A"]
    d = invariants["D"]
    c = invariants["C"]
    r = invariants["R"]
    return [
        (a * r**2 - a * d + c**2, 2 * c * r, d),
        (-2 * c * r, -2 * d),
        (d - 2 * a * r, -2 * c),
        (2 * c,),
        (a,),
    ]


def y_derivative_coefficients(coefficients: Sequence[Poly]) -> list[Poly]:
    return [scale(coefficients[index], index) for index in range(1, len(coefficients))]


def evaluate_rational_lift_numerator(
    coefficients: Sequence[Poly], denominator: Poly, numerator: Poly
) -> Poly:
    """Numerator after substituting y=-B/A into sum f_i y^i."""
    degree = len(coefficients) - 1
    result = zero()
    negative_numerator = neg(numerator)
    for index, coefficient in enumerate(coefficients):
        term = mul(coefficient, power(negative_numerator, index))
        term = mul(term, power(denominator, degree - index))
        result = add(result, term)
    return trim(result)


def unit_flux_polynomial(cosine: Fraction) -> Poly:
    coefficients_descending = [
        Fraction(1_000_000),
        2_000_000 * cosine,
        1_110_000 * cosine**2 - 2_330_000,
        180_000 * cosine**3 - 6_060_000 * cosine,
        -5_164_200 * cosine**2 - 143_700,
        -1_798_200 * cosine**3 - 57_600 * cosine,
        -218_700 * cosine**4 + 6_831 * cosine**2 + 3_969,
    ]
    return tuple(Fraction(value) for value in reversed(coefficients_descending))


def apply_annihilation(state: tuple[int, ...], orbital: int) -> tuple[tuple[int, ...], int]:
    if orbital not in state:
        return state, 0
    position = state.index(orbital)
    return state[:position] + state[position + 1 :], -1 if position % 2 else 1


def apply_creation(state: tuple[int, ...], orbital: int) -> tuple[tuple[int, ...], int]:
    if orbital in state:
        return state, 0
    sign_value = -1 if sum(value < orbital for value in state) % 2 else 1
    return tuple(sorted(state + (orbital,))), sign_value


def apply_pair_hop(
    state: tuple[int, ...], created: tuple[int, int], annihilated: tuple[int, int]
) -> tuple[tuple[int, ...], int]:
    sign_value = 1
    current = state
    for orbital in annihilated:
        current, local = apply_annihilation(current, orbital)
        if not local:
            return state, 0
        sign_value *= local
    for orbital in reversed(created):
        current, local = apply_creation(current, orbital)
        if not local:
            return state, 0
        sign_value *= local
    return current, sign_value


def verify_physical_realization(artifact: dict[str, object]) -> None:
    physical = artifact["physical_operator_realization"]
    carrier = tuple(tuple(determinant) for determinant in physical["carrier_determinants"])
    expected = {
        (0, 1): ((3, 4), (1, 2), 1),
        (0, 2): ((3, 5), (0, 2), -1),
        (1, 2): ((1, 5), (0, 4), 1),
    }
    for (source, target), (created, annihilated, expected_sign) in expected.items():
        result, actual_sign = apply_pair_hop(carrier[source], created, annihilated)
        if result != carrier[target] or actual_sign != expected_sign:
            raise AssertionError("fermionic pair-channel sign mismatch")
        for other_index, determinant in enumerate(carrier):
            if other_index == source:
                continue
            result, actual_sign = apply_pair_hop(determinant, created, annihilated)
            if actual_sign and result not in carrier:
                raise AssertionError("pair channel leaks from the carrier")

    selectors = ((0, 1), (0, 3), (1, 3))
    signatures = [
        [int(set(pair).issubset(determinant)) for determinant in carrier]
        for pair in selectors
    ]
    if signatures != [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
        raise AssertionError("density-density diagonal selectors are not an identity minor")


def verify_regression(record: dict[str, object]) -> None:
    matrix = parse_matrix(record)
    invariants = exact_invariants(matrix)
    stored = {key: frac(value) for key, value in record["invariants"].items()}
    if invariants != stored:
        raise AssertionError(f"invariant mismatch for {record['name']}")
    if invariants["delta"] == 0:
        raise AssertionError("complex regression has zero cycle flux")

    curve_coefficients = branch_curve_coefficients(invariants)
    derivative_coefficients = y_derivative_coefficients(curve_coefficients)
    resultant = sylvester_resultant(
        list(reversed(curve_coefficients)), list(reversed(derivative_coefficients))
    )
    resultant_primitive = primitive_integer(resultant)

    critical = parse_poly_descending(
        record["centered_critical_polynomial_coefficients_descending"]
    )
    if primitive_integer(critical) != resultant_primitive:
        raise AssertionError(f"resultant mismatch for {record['name']}")
    if len(critical) - 1 != 6:
        raise AssertionError("critical polynomial is not sextic")
    if gcd_poly(critical, derivative(critical)) != one():
        raise AssertionError("critical polynomial is not squarefree")

    lift = record["linear_lift"]
    lift_denominator = parse_poly_descending(lift["A_coefficients_descending"])
    lift_numerator = parse_poly_descending(lift["B_coefficients_descending"])
    if gcd_poly(critical, lift_denominator) != one():
        raise AssertionError("lift denominator is not coprime to the sextic")
    curve_after_lift = evaluate_rational_lift_numerator(
        curve_coefficients, lift_denominator, lift_numerator
    )
    derivative_after_lift = evaluate_rational_lift_numerator(
        derivative_coefficients, lift_denominator, lift_numerator
    )
    if divmod_poly(curve_after_lift, critical)[1] != zero():
        raise AssertionError("lift does not satisfy F=0 modulo the sextic")
    if divmod_poly(derivative_after_lift, critical)[1] != zero():
        raise AssertionError("lift does not satisfy dF/dy=0 modulo the sextic")

    cycle = g_mul(g_mul(matrix[0][1], matrix[1][2]), g_conj(matrix[0][2]))
    stored_cycle = gaussian(record["cycle_product"])
    if cycle != stored_cycle:
        raise AssertionError("cycle product mismatch")
    if g_abs2(cycle) == 1:
        expected_family = unit_flux_polynomial(cycle[0])
        if primitive_integer(expected_family) != primitive_integer(critical):
            raise AssertionError("unit-flux family polynomial mismatch")
        reflected = tuple(
            coefficient * ((-1) ** index) for index, coefficient in enumerate(
                unit_flux_polynomial(-cycle[0])
            )
        )
        if reflected != unit_flux_polynomial(cycle[0]):
            raise AssertionError("P_{-u}(-z)=P_u(z) symmetry failed")

    intervals = [
        (frac(item["lower"]), frac(item["upper"]), int(item["multiplicity"]))
        for item in record["real_roots_in_increasing_order"]
    ]
    if any(left[1] >= right[0] for left, right in zip(intervals, intervals[1:])):
        raise AssertionError("root intervals overlap or are out of order")
    sequence = sturm_sequence(critical)
    if total_real_roots(sequence) != len(intervals):
        raise AssertionError("stored intervals do not exhaust the real roots")
    for lower, upper, multiplicity in intervals:
        if multiplicity != 1:
            raise AssertionError("regression root is not simple")
        if roots_between(sequence, lower, upper) != 1:
            raise AssertionError("stored interval is not a one-root Sturm enclosure")

    exact_range = record["exact_range"]
    if exact_range["minimum_root_index"] != 0:
        raise AssertionError("minimum is not the first real root")
    if exact_range["maximum_root_index"] != len(intervals) - 1:
        raise AssertionError("maximum is not the last real root")
    minimum = exact_range["minimum_isolating_interval"]
    maximum = exact_range["maximum_isolating_interval"]
    if (frac(minimum["lower"]), frac(minimum["upper"])) != intervals[0][:2]:
        raise AssertionError("minimum interval mismatch")
    if (frac(maximum["lower"]), frac(maximum["upper"])) != intervals[-1][:2]:
        raise AssertionError("maximum interval mismatch")


def verify_artifact(path: Path) -> None:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("schema_version") != 1:
        raise AssertionError("unsupported schema version")
    if artifact["target_1rdm_diagonal"] != [
        "17/20",
        "3/4",
        "3/5",
        "2/5",
        "1/4",
        "3/20",
    ]:
        raise AssertionError("target 1-RDM mismatch")
    verify_physical_realization(artifact)
    regressions = artifact["regressions"]
    if [record["name"] for record in regressions] != [
        "quarter_flux",
        "pythagorean_flux_cos_3_5_sin_4_5",
    ]:
        raise AssertionError("unexpected regression set")
    for record in regressions:
        verify_regression(record)

    # The unit-family limits reproduce the two real-symmetric checkpoints.
    plus = unit_flux_polynomial(Fraction(1))
    minus = unit_flux_polynomial(Fraction(-1))
    if evaluate(plus, Fraction(-1)) or evaluate(derivative(plus), Fraction(-1)):
        raise AssertionError("zero-flux interior minimum is not a double root")
    if evaluate(minus, Fraction(1)) or evaluate(derivative(minus), Fraction(1)):
        raise AssertionError("pi-flux interior maximum is not a double root")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", nargs="?", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    try:
        verify_artifact(args.artifact)
    except Exception as error:  # noqa: BLE001 - standalone verifier reports all failures
        print(f"FAIL: {error}")
        return 1
    print("all 2 complex fixed-1RDM observable-range certificates PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
