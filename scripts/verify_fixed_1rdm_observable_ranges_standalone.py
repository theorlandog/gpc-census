#!/usr/bin/env python3
"""Independent standard-library verifier for fixed_1rdm_observable_ranges.json."""
from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT = ROOT / "results" / "data" / "fixed_1rdm_observable_ranges.json"

SUPPORT = ((0, 1, 2), (0, 3, 4), (1, 3, 5))
WEIGHTS = (Fraction(3, 5), Fraction(1, 4), Fraction(3, 20))
PAIRING = ((0, 5), (1, 4), (2, 3))
FACET_ORBITALS = frozenset((0, 1, 3))
UNFRUSTRATED = ((0, 1, 1), (1, 0, 1), (1, 1, 0))
PI_FLUX = ((0, 1, 1), (1, 0, -1), (1, -1, 0))


@dataclass(frozen=True)
class Qsqrt15:
    """Exact element a + b*sqrt(15), with a,b rational."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __add__(self, other: object) -> "Qsqrt15":
        other = coerce(other)
        return Qsqrt15(self.a + other.a, self.b + other.b)

    def __radd__(self, other: object) -> "Qsqrt15":
        return self + other

    def __neg__(self) -> "Qsqrt15":
        return Qsqrt15(-self.a, -self.b)

    def __sub__(self, other: object) -> "Qsqrt15":
        return self + (-coerce(other))

    def __rsub__(self, other: object) -> "Qsqrt15":
        return coerce(other) - self

    def __mul__(self, other: object) -> "Qsqrt15":
        other = coerce(other)
        return Qsqrt15(
            self.a * other.a + 15 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def __rmul__(self, other: object) -> "Qsqrt15":
        return self * other

    def __truediv__(self, other: object) -> "Qsqrt15":
        other = coerce(other)
        norm = other.a * other.a - 15 * other.b * other.b
        if norm == 0:
            raise ZeroDivisionError
        return self * Qsqrt15(other.a / norm, -other.b / norm)

    def __pow__(self, exponent: int) -> "Qsqrt15":
        if exponent < 0:
            raise ValueError("negative powers are not needed")
        result = Qsqrt15(Fraction(1))
        base = self
        power = exponent
        while power:
            if power & 1:
                result *= base
            base *= base
            power >>= 1
        return result


def coerce(value: object) -> Qsqrt15:
    if isinstance(value, Qsqrt15):
        return value
    if isinstance(value, (int, Fraction)):
        return Qsqrt15(Fraction(value))
    raise TypeError(value)


def determinant_bareiss(matrix: list[list[int]]) -> int:
    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign *= -1
        pivot = work[column][column]
        for row in range(column + 1, size):
            for other in range(column + 1, size):
                numerator = (
                    work[row][other] * pivot
                    - work[row][column] * work[column][other]
                )
                if numerator % previous:
                    raise AssertionError("non-exact Bareiss division")
                work[row][other] = numerator // previous
        previous = pivot
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


def permutation_sign(sequence: tuple[int, ...]) -> int:
    inversions = sum(
        sequence[i] > sequence[j]
        for i in range(len(sequence))
        for j in range(i + 1, len(sequence))
    )
    return -1 if inversions % 2 else 1


def signed_channels(
    support: tuple[tuple[int, ...], ...], order: int
) -> dict[tuple[tuple[int, ...], tuple[int, ...]], dict[tuple[int, int], int]]:
    particles = len(support[0])
    terms = defaultdict(lambda: defaultdict(int))
    for left_index, left in enumerate(support):
        for right_index, right in enumerate(support):
            intersection = sorted(set(left).intersection(right))
            remainder_size = particles - order
            if len(intersection) < remainder_size:
                continue
            for remainder in itertools.combinations(intersection, remainder_size):
                remainder_set = set(remainder)
                left_modes = tuple(mode for mode in left if mode not in remainder_set)
                right_modes = tuple(mode for mode in right if mode not in remainder_set)
                sign = permutation_sign(left_modes + remainder)
                sign *= permutation_sign(right_modes + remainder)
                terms[(left_modes, right_modes)][(left_index, right_index)] += sign
    return {
        channel: {
            edge: coefficient
            for edge, coefficient in coefficients.items()
            if coefficient
        }
        for channel, coefficients in terms.items()
        if any(coefficients.values())
    }


def verify_selection_rule(payload: dict) -> None:
    determinants = list(itertools.combinations(range(6), 3))
    single_occupancy = [
        determinant
        for determinant in determinants
        if all(
            sum(orbital in determinant for orbital in pair) == 1
            for pair in PAIRING
        )
    ]
    facet_kernel = [
        determinant
        for determinant in determinants
        if 2 - sum(orbital in FACET_ORBITALS for orbital in determinant) == 0
    ]
    intersection = tuple(sorted(set(single_occupancy).intersection(facet_kernel)))
    assert len(determinants) == 20
    assert len(single_occupancy) == 8
    assert len(facet_kernel) == 9
    assert intersection == SUPPORT

    certificate = payload["selection_rule_certificate"]
    assert certificate["intersection_support"] == [list(row) for row in SUPPORT]
    assert certificate["equals_benchmark_support"] is True


def verify_one_body_fiber(payload: dict) -> None:
    incidence = [
        [int(orbital in determinant) for determinant in SUPPORT]
        for orbital in range(6)
    ]
    pivot_rows = [0, 1, 2]
    minor = [incidence[index] for index in pivot_rows]
    assert determinant_bareiss(minor) == 1

    gamma = tuple(
        sum(WEIGHTS[index] * incidence[orbital][index] for index in range(3))
        for orbital in range(6)
    )
    expected = (
        Fraction(17, 20),
        Fraction(3, 4),
        Fraction(3, 5),
        Fraction(2, 5),
        Fraction(1, 4),
        Fraction(3, 20),
    )
    assert gamma == expected
    assert all(gamma[index] > gamma[index + 1] for index in range(5))
    assert sum(gamma) == 3
    assert gamma[0] + gamma[5] == 1
    assert gamma[1] + gamma[4] == 1
    assert gamma[2] + gamma[3] == 1
    assert 2 - (gamma[0] + gamma[1] + gamma[3]) == 0

    # No support pair differs by one excitation, so all one-body off-diagonal
    # entries vanish identically and phases do not affect the complete 1-RDM.
    assert all(
        len(set(SUPPORT[left]).intersection(SUPPORT[right])) != 2
        for left, right in itertools.combinations(range(3), 2)
    )
    assert payload["target_1rdm"]["diagonal"] == [str(value) for value in gamma]
    assert payload["fiber"]["weights"] == [str(value) for value in WEIGHTS]


def verify_two_body_witnesses(payload: dict) -> None:
    channels = signed_channels(SUPPORT, 2)
    witnesses = payload["observable_realization"]["collision_free_witnesses"]
    assert len(witnesses) == 3
    seen = set()
    for witness in witnesses:
        channel = (
            tuple(witness["channel_left_modes"]),
            tuple(witness["channel_right_modes"]),
        )
        terms = channels[channel]
        assert len(terms) == 1
        ((oriented_edge, coefficient),) = terms.items()
        edge = tuple(sorted(oriented_edge))
        assert edge == tuple(witness["edge"])
        assert coefficient == witness["channel_coefficient"]
        assert Fraction(witness["operator_multiplier"]) * coefficient == 1
        seen.add(edge)
    assert seen == set(itertools.combinations(range(3), 2))


def verify_operator_decomposition(example: dict, expected: tuple[tuple[int, ...], ...]) -> None:
    decomposition = example["operator_decomposition"]
    reconstructed = [[Fraction(0) for _ in range(3)] for _ in range(3)]

    for term in decomposition["one_body_terms"]:
        coefficient = Fraction(term["coefficient"])
        orbital = int(term["orbital"])
        for index, determinant in enumerate(SUPPORT):
            reconstructed[index][index] += coefficient * int(orbital in determinant)

    channels = signed_channels(SUPPORT, 2)
    for term in decomposition["two_body_terms"]:
        channel = (
            tuple(term["channel_left_modes"]),
            tuple(term["channel_right_modes"]),
        )
        terms = channels[channel]
        assert len(terms) == 1
        ((oriented_edge, coefficient),) = terms.items()
        edge = tuple(sorted(oriented_edge))
        assert edge == tuple(term["edge"])
        matrix_element = Fraction(term["operator_multiplier"]) * coefficient
        left, right = edge
        reconstructed[left][right] += matrix_element
        reconstructed[right][left] += matrix_element

    expected_fractions = [[Fraction(value) for value in row] for row in expected]
    assert reconstructed == expected_fractions
    assert decomposition["reconstructed_carrier_matrix"] == [
        [str(value) for value in row] for row in expected
    ]


def cosine_gram(cos01: Qsqrt15, cos02: Qsqrt15, cos12: Qsqrt15) -> Qsqrt15:
    return (
        1
        + 2 * cos01 * cos02 * cos12
        - cos01 * cos01
        - cos02 * cos02
        - cos12 * cos12
    )


def verify_ranges(payload: dict) -> None:
    # sqrt(w_0)=sqrt(15)/5, sqrt(w_1)=1/2, sqrt(w_2)=sqrt(15)/10.
    radii = (
        Qsqrt15(Fraction(0), Fraction(1, 5)),
        Qsqrt15(Fraction(1, 2), Fraction(0)),
        Qsqrt15(Fraction(0), Fraction(1, 10)),
    )
    assert tuple(radius * radius for radius in radii) == tuple(
        Qsqrt15(weight) for weight in WEIGHTS
    )
    scales = (
        2 * radii[0] * radii[1],
        2 * radii[0] * radii[2],
        2 * radii[1] * radii[2],
    )
    assert scales == (
        Qsqrt15(Fraction(0), Fraction(1, 5)),
        Qsqrt15(Fraction(3, 5)),
        Qsqrt15(Fraction(0), Fraction(1, 10)),
    )

    examples = payload["examples"]
    unfrustrated = examples["unfrustrated_coherence_triangle"]
    pi_flux = examples["pi_flux_coherence_triangle"]
    verify_operator_decomposition(unfrustrated, UNFRUSTRATED)
    verify_operator_decomposition(pi_flux, PI_FLUX)

    unfrustrated_min_cosines = (
        Qsqrt15(Fraction(0), Fraction(-7, 30)),
        Qsqrt15(Fraction(-5, 6)),
        Qsqrt15(Fraction(0), Fraction(2, 15)),
    )
    assert cosine_gram(*unfrustrated_min_cosines) == Qsqrt15()
    unfrustrated_minimum = sum(
        (scale * cosine for scale, cosine in zip(scales, unfrustrated_min_cosines)),
        Qsqrt15(),
    )
    unfrustrated_maximum = sum(scales, Qsqrt15())
    assert unfrustrated_minimum == Qsqrt15(Fraction(-1))
    assert unfrustrated_maximum == Qsqrt15(Fraction(3, 5), Fraction(3, 10))

    first_range = unfrustrated["range"]
    assert first_range["minimum"] == "-1"
    assert first_range["maximum"] == "3*(2 + sqrt(15))/10"
    assert first_range["minimum_source"] == "interior stationary point"
    assert first_range["maximum_source"] == "alpha=0 endpoint"
    assert first_range["minimum_phase_cosines"] == {
        "cos_alpha": "-7*sqrt(15)/30",
        "cos_beta": "-5/6",
        "cos_beta_minus_alpha": "2*sqrt(15)/15",
    }
    assert first_range["candidate_reduction"]["candidate_count"] == 5

    pi_max_cosines = (
        Qsqrt15(Fraction(0), Fraction(7, 30)),
        Qsqrt15(Fraction(5, 6)),
        Qsqrt15(Fraction(0), Fraction(2, 15)),
    )
    assert cosine_gram(*pi_max_cosines) == Qsqrt15()
    pi_maximum = (
        scales[0] * pi_max_cosines[0]
        + scales[1] * pi_max_cosines[1]
        - scales[2] * pi_max_cosines[2]
    )
    pi_minimum = -sum(scales, Qsqrt15())
    assert pi_maximum == Qsqrt15(Fraction(1))
    assert pi_minimum == Qsqrt15(Fraction(-3, 5), Fraction(-3, 10))

    second_range = pi_flux["range"]
    assert second_range["minimum"] == "-3*(2 + sqrt(15))/10"
    assert second_range["maximum"] == "1"
    assert second_range["minimum_source"] == "alpha=pi endpoint"
    assert second_range["maximum_source"] == "interior stationary point"
    assert second_range["maximum_phase_cosines"] == {
        "cos_alpha": "7*sqrt(15)/30",
        "cos_beta": "5/6",
        "cos_beta_minus_alpha": "2*sqrt(15)/15",
    }
    assert second_range["candidate_reduction"]["candidate_count"] == 5

    polygon = unfrustrated["independent_polygon_certificate"]
    assert polygon["minimum"] == first_range["minimum"]
    assert polygon["maximum"] == first_range["maximum"]


def verify(path: Path) -> None:
    payload = json.loads(path.read_text())
    assert payload["schema_version"] == 2
    assert payload["evidence"]["status"] == "EXACT"
    assert payload["range_engine"]["maximum_candidate_count"] == 5
    verify_selection_rule(payload)
    verify_one_body_fiber(payload)
    verify_two_body_witnesses(payload)
    verify_ranges(payload)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", nargs="?", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    verify(args.artifact)
    print(
        "verified fixed 1-RDM range certificate: selection rule, fiber, "
        "operator realization, unfrustrated range, and pi-flux range"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
