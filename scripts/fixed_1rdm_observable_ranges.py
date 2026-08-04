#!/usr/bin/env python3
"""Exact two-body observable ranges on a pinned Borland--Dennis 1-RDM fiber.

The benchmark is the nondegenerate pinned point

    gamma = diag(17/20, 3/4, 3/5, 2/5, 1/4, 3/20).

The Borland--Dennis equalities and saturated generalized Pauli constraint select
exactly three determinants in the natural-orbital frame.  The diagonal 1-RDM
then fixes their squared amplitudes, leaving a two-torus of phases.  Every real
symmetric rational carrier Hamiltonian is realized by at-most-two-body terms
and its global range is obtained from four endpoint values plus at most one
interior critical value.

This is a global pure-state fixed-1RDM result for this target gamma, not merely
a range on an arbitrarily imposed support: the exact pinning selection rules
supply the support reduction and the occupations are strictly nondegenerate.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from pathlib import Path
from typing import Literal

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "data" / "fixed_1rdm_observable_ranges.json"
GENERATOR = "scripts/fixed_1rdm_observable_ranges.py"
SQRT15 = sp.sqrt(15)

# Natural-orbital labels are zero based.  In one-based notation this is the
# standard pinned Borland--Dennis carrier |123>, |145>, |246>.
BD_SUPPORT: tuple[tuple[int, ...], ...] = (
    (0, 1, 2),
    (0, 3, 4),
    (1, 3, 5),
)
BD_WEIGHTS: tuple[sp.Rational, ...] = (
    sp.Rational(3, 5),
    sp.Rational(1, 4),
    sp.Rational(3, 20),
)
BD_PAIRING: tuple[tuple[int, int], ...] = ((0, 5), (1, 4), (2, 3))
BD_FACET_ORBITALS = frozenset((0, 1, 3))

COHERENCE_TRIANGLE_MATRIX: tuple[tuple[sp.Rational, ...], ...] = (
    (sp.Rational(0), sp.Rational(1), sp.Rational(1)),
    (sp.Rational(1), sp.Rational(0), sp.Rational(1)),
    (sp.Rational(1), sp.Rational(1), sp.Rational(0)),
)
PI_FLUX_TRIANGLE_MATRIX: tuple[tuple[sp.Rational, ...], ...] = (
    (sp.Rational(0), sp.Rational(1), sp.Rational(1)),
    (sp.Rational(1), sp.Rational(0), sp.Rational(-1)),
    (sp.Rational(1), sp.Rational(-1), sp.Rational(0)),
)


def exact_string(value: sp.Expr) -> str:
    """Stable human-readable exact expression for JSON artifacts."""
    return sp.sstr(sp.factor(sp.radsimp(sp.simplify(value))))


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
    """Return oriented k-RDM channels in X_ab=conjugate(c_a)c_b variables."""
    particles = len(support[0])
    terms: dict[
        tuple[tuple[int, ...], tuple[int, ...]], dict[tuple[int, int], int]
    ] = defaultdict(lambda: defaultdict(int))
    for left_index, left in enumerate(support):
        left_set = set(left)
        for right_index, right in enumerate(support):
            intersection = sorted(left_set.intersection(right))
            remainder_size = particles - order
            if len(intersection) < remainder_size:
                continue
            for remainder in itertools.combinations(intersection, remainder_size):
                remainder_set = set(remainder)
                left_modes = tuple(mode for mode in left if mode not in remainder_set)
                right_modes = tuple(mode for mode in right if mode not in remainder_set)
                if len(left_modes) != order or len(right_modes) != order:
                    continue
                sign = permutation_sign(left_modes + remainder)
                sign *= permutation_sign(right_modes + remainder)
                terms[(left_modes, right_modes)][(left_index, right_index)] += sign

    reduced = {}
    for channel, coefficients in terms.items():
        nonzero = {edge: value for edge, value in coefficients.items() if value}
        if nonzero:
            reduced[channel] = nonzero
    return reduced


def one_body_incidence(
    support: tuple[tuple[int, ...], ...], orbitals: int
) -> sp.Matrix:
    return sp.Matrix(
        [
            [int(orbital in determinant) for determinant in support]
            for orbital in range(orbitals)
        ]
    )


def one_hop_pairs(support: tuple[tuple[int, ...], ...]) -> list[tuple[int, int]]:
    particles = len(support[0])
    return [
        (left, right)
        for left, right in itertools.combinations(range(len(support)), 2)
        if len(set(support[left]).intersection(support[right])) == particles - 1
    ]


def diagonal_1rdm(
    support: tuple[tuple[int, ...], ...],
    weights: tuple[sp.Rational, ...],
    orbitals: int,
) -> tuple[sp.Expr, ...]:
    incidence = one_body_incidence(support, orbitals)
    values = incidence * sp.Matrix(weights)
    return tuple(sp.simplify(value) for value in values)


def recover_weights(
    support: tuple[tuple[int, ...], ...], gamma_diagonal: tuple[sp.Expr, ...]
) -> tuple[tuple[sp.Expr, ...], dict]:
    incidence = one_body_incidence(support, len(gamma_diagonal))
    if incidence.rank() != len(support):
        raise ValueError("singleton incidence does not uniquely recover support weights")
    solution = sp.linsolve((incidence, sp.Matrix(gamma_diagonal)))
    rows = list(solution)
    if len(rows) != 1:
        raise ValueError("target diagonal has no unique support-weight solution")
    recovered = tuple(sp.simplify(value) for value in rows[0])

    # Pivot columns of A^T are independent rows of A.
    pivot_rows = tuple(int(index) for index in incidence.T.rref()[1])
    minor = incidence[list(pivot_rows), :]
    certificate = {
        "rank": int(incidence.rank()),
        "pivot_orbitals": list(pivot_rows),
        "pivot_minor_determinant": int(minor.det()),
    }
    return recovered, certificate


def borland_dennis_selection_certificate() -> dict:
    """Enumerate the exact support reduction from the BD equalities and facet."""
    determinants = list(itertools.combinations(range(6), 3))
    single_occupancy = [
        determinant
        for determinant in determinants
        if all(
            sum(orbital in determinant for orbital in pair) == 1
            for pair in BD_PAIRING
        )
    ]
    facet_kernel = [
        determinant
        for determinant in determinants
        if 2 - sum(orbital in BD_FACET_ORBITALS for orbital in determinant) == 0
    ]
    intersection = sorted(set(single_occupancy).intersection(facet_kernel))
    return {
        "ambient_determinants": len(determinants),
        "single_occupancy_pairs": [list(pair) for pair in BD_PAIRING],
        "single_occupancy_support": [list(row) for row in single_occupancy],
        "single_occupancy_support_size": len(single_occupancy),
        "facet_operator": "2 - (n_0 + n_1 + n_3)",
        "facet_kernel_support": [list(row) for row in facet_kernel],
        "facet_kernel_support_size": len(facet_kernel),
        "intersection_support": [list(row) for row in intersection],
        "intersection_support_size": len(intersection),
        "equals_benchmark_support": bool(tuple(intersection) == BD_SUPPORT),
    }


def complete_collision_free_witnesses(
    support: tuple[tuple[int, ...], ...], order: int = 2
) -> list[dict]:
    """Choose one collision-free Hermitian channel for every support pair."""
    channels = signed_channels(support, order)
    by_edge: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for (left_modes, right_modes), terms in sorted(channels.items()):
        # Keep one orientation of each Hermitian channel.
        if left_modes >= right_modes or len(terms) != 1:
            continue
        (left_index, right_index), coefficient = next(iter(terms.items()))
        if left_index == right_index:
            continue
        edge = tuple(sorted((left_index, right_index)))
        by_edge[edge].append(
            {
                "edge": list(edge),
                "channel_left_modes": list(left_modes),
                "channel_right_modes": list(right_modes),
                "channel_coefficient": coefficient,
                "operator_multiplier": exact_string(sp.Rational(1, coefficient)),
            }
        )

    witnesses = []
    for edge in itertools.combinations(range(len(support)), 2):
        candidates = by_edge.get(edge, [])
        if not candidates:
            raise ValueError(f"no collision-free order-{order} channel for edge {edge}")
        witnesses.append(candidates[0])
    return witnesses


def triangle_closes(weights: tuple[sp.Rational, sp.Rational, sp.Rational]) -> bool:
    """Decide the triangle inequality for square-root radii using rationals only."""
    ordered = sorted(weights, reverse=True)
    largest, middle, smallest = ordered
    # sqrt(largest) <= sqrt(middle) + sqrt(smallest).
    difference = largest - middle - smallest
    if difference <= 0:
        return True
    return bool(difference * difference <= 4 * middle * smallest)


def complete_triangle_range(
    weights: tuple[sp.Rational, sp.Rational, sp.Rational]
) -> dict:
    """Exact range of sum_{a<b}(X_ab+X_ba) at fixed positive weights."""
    if any(weight <= 0 for weight in weights):
        raise ValueError("the benchmark theorem assumes full support")
    if sum(weights) != 1:
        raise ValueError("weights must sum to one")

    radii = tuple(sp.sqrt(weight) for weight in weights)
    maximum = sp.factor(sum(radii) ** 2 - 1)
    closes = triangle_closes(weights)
    if closes:
        minimum = sp.Integer(-1)
    else:
        largest_index = max(range(3), key=lambda index: weights[index])
        deficit = radii[largest_index] - sum(
            radii[index] for index in range(3) if index != largest_index
        )
        minimum = sp.factor(deficit**2 - 1)

    # If the triangle closes, these are the exact pairwise cosines of a
    # minimizing phase configuration obtained from the law of cosines.
    minimum_cosines = None
    if closes:
        minimum_cosines = {}
        for left, right in itertools.combinations(range(3), 2):
            third = next(index for index in range(3) if index not in (left, right))
            cosine = (
                weights[third] - weights[left] - weights[right]
            ) / (2 * radii[left] * radii[right])
            minimum_cosines[f"{left}-{right}"] = exact_string(cosine)

    return {
        "polygon_closes": closes,
        "minimum": exact_string(minimum),
        "maximum": exact_string(maximum),
        "width": exact_string(maximum - minimum),
        "minimum_numeric": float(sp.N(minimum, 18)),
        "maximum_numeric": float(sp.N(maximum, 18)),
        "minimum_pair_phase_cosines": minimum_cosines,
        "maximum_phase_condition": "all three coefficient phases are equal",
    }


def _q15_parts(value: sp.Expr) -> tuple[sp.Rational, sp.Rational]:
    """Write a value in Q(sqrt(15)) as a+b*sqrt(15), exactly."""
    value = sp.radsimp(sp.simplify(value))
    conjugate = sp.radsimp(value.xreplace({SQRT15: -SQRT15}))
    rational = sp.simplify((value + conjugate) / 2)
    radical = sp.simplify((value - conjugate) / (2 * SQRT15))
    if rational.is_Rational is not True or radical.is_Rational is not True:
        raise ValueError(f"value is not in Q(sqrt(15)): {value}")
    return sp.Rational(rational), sp.Rational(radical)


def q15_sign(value: sp.Expr) -> int:
    """Return the exact sign of a Q(sqrt(15)) element."""
    rational, radical = _q15_parts(value)
    if rational == 0:
        return 0 if radical == 0 else (1 if radical > 0 else -1)
    if radical == 0:
        return 1 if rational > 0 else -1
    if rational > 0 and radical > 0:
        return 1
    if rational < 0 and radical < 0:
        return -1

    comparison = rational**2 - 15 * radical**2
    if comparison == 0:
        return 0
    if rational > 0:  # radical < 0
        return 1 if comparison > 0 else -1
    # rational < 0 and radical > 0
    return 1 if comparison < 0 else -1


def q15_abs(value: sp.Expr) -> sp.Expr:
    sign = q15_sign(value)
    return sp.Integer(0) if sign == 0 else sp.simplify(value if sign > 0 else -value)


def _validate_real_symmetric_matrix(
    matrix: tuple[tuple[sp.Expr, ...], ...] | list[list[sp.Expr]],
) -> sp.Matrix:
    hamiltonian = sp.Matrix(matrix)
    if hamiltonian.shape != (3, 3):
        raise ValueError("the pinned Borland-Dennis carrier matrix must be 3 by 3")
    if hamiltonian != hamiltonian.T:
        raise ValueError("the exact phase-elimination solver currently requires real symmetry")
    if any(sp.sympify(entry).is_Rational is not True for entry in hamiltonian):
        raise ValueError("carrier matrix entries must be rational")
    return hamiltonian.applyfunc(sp.Rational)


def carrier_operator_decomposition(
    matrix: tuple[tuple[sp.Expr, ...], ...] | list[list[sp.Expr]],
    witnesses: list[dict] | None = None,
) -> dict:
    """Realize a rational real-symmetric carrier matrix by <=2-body terms."""
    hamiltonian = _validate_real_symmetric_matrix(matrix)
    witnesses = (
        complete_collision_free_witnesses(BD_SUPPORT, order=2)
        if witnesses is None
        else witnesses
    )

    incidence = one_body_incidence(BD_SUPPORT, 6)
    pivot_orbitals = tuple(int(index) for index in incidence.T.rref()[1])
    pivot = incidence[list(pivot_orbitals), :]
    diagonal = sp.Matrix([hamiltonian[index, index] for index in range(3)])
    one_body_coefficients = pivot.T.inv() * diagonal

    one_body_terms = [
        {
            "orbital": orbital,
            "coefficient": exact_string(coefficient),
            "operator": f"n_{orbital}",
        }
        for orbital, coefficient in zip(pivot_orbitals, one_body_coefficients)
        if coefficient != 0
    ]

    by_edge = {tuple(witness["edge"]): witness for witness in witnesses}
    two_body_terms = []
    for left, right in itertools.combinations(range(3), 2):
        coefficient = hamiltonian[left, right]
        if coefficient == 0:
            continue
        witness = by_edge[(left, right)]
        channel_sign = sp.Rational(witness["channel_coefficient"])
        two_body_terms.append(
            {
                "edge": [left, right],
                "channel_left_modes": witness["channel_left_modes"],
                "channel_right_modes": witness["channel_right_modes"],
                "channel_coefficient": witness["channel_coefficient"],
                "operator_multiplier": exact_string(coefficient / channel_sign),
                "operator": (
                    "multiplier * (A_left^dagger A_right + "
                    "A_right^dagger A_left)"
                ),
            }
        )

    reconstructed = sp.zeros(3)
    for term in one_body_terms:
        orbital = int(term["orbital"])
        coefficient = sp.sympify(term["coefficient"])
        for determinant_index, determinant in enumerate(BD_SUPPORT):
            reconstructed[determinant_index, determinant_index] += (
                coefficient * int(orbital in determinant)
            )
    for term in two_body_terms:
        left, right = term["edge"]
        multiplier = sp.sympify(term["operator_multiplier"])
        channel_sign = sp.Rational(term["channel_coefficient"])
        reconstructed[left, right] += multiplier * channel_sign
        reconstructed[right, left] += multiplier * channel_sign
    if reconstructed != hamiltonian:
        raise AssertionError("second-quantized decomposition failed to reconstruct the matrix")

    return {
        "one_body_terms": one_body_terms,
        "two_body_terms": two_body_terms,
        "reconstructed_carrier_matrix": [
            [exact_string(reconstructed[row, column]) for column in range(3)]
            for row in range(3)
        ],
    }


def _phase_candidate(
    *,
    source: str,
    branch: Literal["upper", "lower"],
    value: sp.Expr,
    cos_alpha: sp.Expr,
    cos_beta: sp.Expr,
    cos_beta_minus_alpha: sp.Expr,
) -> dict:
    gram = sp.factor(
        1
        + 2 * cos_alpha * cos_beta * cos_beta_minus_alpha
        - cos_alpha**2
        - cos_beta**2
        - cos_beta_minus_alpha**2
    )
    if sp.simplify(gram) != 0:
        raise AssertionError(f"phase-cosine Gram determinant is not zero: {gram}")
    return {
        "source": source,
        "branch": branch,
        "value": sp.factor(sp.radsimp(value)),
        "cos_alpha": sp.factor(sp.radsimp(cos_alpha)),
        "cos_beta": sp.factor(sp.radsimp(cos_beta)),
        "cos_beta_minus_alpha": sp.factor(sp.radsimp(cos_beta_minus_alpha)),
    }


def _endpoint_candidate(
    *,
    offset: sp.Expr,
    coupling_alpha: sp.Expr,
    coupling_beta: sp.Expr,
    coupling_difference: sp.Expr,
    cos_alpha: int,
    branch: Literal["upper", "lower"],
) -> dict:
    combined = coupling_beta + cos_alpha * coupling_difference
    combined_sign = q15_sign(combined)
    # If the combined beta coefficient vanishes, beta=0 is a valid representative.
    alignment = 1 if combined_sign == 0 else combined_sign
    if branch == "lower":
        alignment *= -1
    cos_beta = sp.Integer(alignment)
    cos_difference = sp.Integer(cos_alpha * alignment)
    value = (
        offset
        + coupling_alpha * cos_alpha
        + coupling_beta * cos_beta
        + coupling_difference * cos_difference
    )
    return _phase_candidate(
        source="alpha=0 endpoint" if cos_alpha == 1 else "alpha=pi endpoint",
        branch=branch,
        value=value,
        cos_alpha=sp.Integer(cos_alpha),
        cos_beta=cos_beta,
        cos_beta_minus_alpha=cos_difference,
    )


def _select_candidate(candidates: list[dict], mode: Literal["minimum", "maximum"]) -> dict:
    selected = candidates[0]
    for candidate in candidates[1:]:
        difference = sp.simplify(candidate["value"] - selected["value"])
        sign = q15_sign(difference)
        if (mode == "minimum" and sign < 0) or (mode == "maximum" and sign > 0):
            selected = candidate
    return selected


def real_symmetric_triangle_range(
    matrix: tuple[tuple[sp.Expr, ...], ...] | list[list[sp.Expr]],
    weights: tuple[sp.Rational, sp.Rational, sp.Rational] = BD_WEIGHTS,
) -> dict:
    """Globally optimize a rational real-symmetric matrix on the phase torus.

    In the gauge c=(sqrt(w0),sqrt(w1)e^{i alpha},sqrt(w2)e^{i beta}),

        E = E0 + A cos(alpha) + B cos(beta) + C cos(beta-alpha).

    At fixed t=cos(alpha), eliminating beta gives E0+A*t +/-
    sqrt(B^2+C^2+2BC*t).  The upper branch is concave and the lower branch is
    convex, so two endpoints and at most one stationary point per relevant
    branch give the exact global interval.
    """
    hamiltonian = _validate_real_symmetric_matrix(matrix)
    if any(weight <= 0 for weight in weights) or sum(weights) != 1:
        raise ValueError("weights must be positive and normalized")

    radii = tuple(sp.sqrt(weight) for weight in weights)
    offset = sp.factor(
        sum(weights[index] * hamiltonian[index, index] for index in range(3))
    )
    coupling_alpha = sp.factor(2 * radii[0] * radii[1] * hamiltonian[0, 1])
    coupling_beta = sp.factor(2 * radii[0] * radii[2] * hamiltonian[0, 2])
    coupling_difference = sp.factor(
        2 * radii[1] * radii[2] * hamiltonian[1, 2]
    )

    # These exact comparisons are implemented in Q(sqrt(15)); rational input
    # and the chosen benchmark weights guarantee every candidate lies there.
    for value in (offset, coupling_alpha, coupling_beta, coupling_difference):
        _q15_parts(value)

    lower_candidates = [
        _endpoint_candidate(
            offset=offset,
            coupling_alpha=coupling_alpha,
            coupling_beta=coupling_beta,
            coupling_difference=coupling_difference,
            cos_alpha=endpoint,
            branch="lower",
        )
        for endpoint in (1, -1)
    ]
    upper_candidates = [
        _endpoint_candidate(
            offset=offset,
            coupling_alpha=coupling_alpha,
            coupling_beta=coupling_beta,
            coupling_difference=coupling_difference,
            cos_alpha=endpoint,
            branch="upper",
        )
        for endpoint in (1, -1)
    ]

    product = sp.factor(coupling_alpha * coupling_beta * coupling_difference)
    if product != 0:
        ratio = sp.factor(coupling_beta * coupling_difference / coupling_alpha)
        stationary_cos_alpha = sp.factor(
            (
                ratio**2
                - coupling_beta**2
                - coupling_difference**2
            )
            / (2 * coupling_beta * coupling_difference)
        )
        if (
            q15_sign(stationary_cos_alpha + 1) >= 0
            and q15_sign(1 - stationary_cos_alpha) >= 0
        ):
            product_sign = q15_sign(product)
            branch: Literal["upper", "lower"] = (
                "lower" if product_sign > 0 else "upper"
            )
            branch_sign = -1 if branch == "lower" else 1
            resultant = q15_abs(ratio)
            cos_beta = sp.factor(
                branch_sign
                * (coupling_beta + coupling_difference * stationary_cos_alpha)
                / resultant
            )
            cos_difference = sp.factor(
                branch_sign
                * (coupling_beta * stationary_cos_alpha + coupling_difference)
                / resultant
            )
            critical_value = sp.factor(
                offset
                + coupling_alpha * stationary_cos_alpha
                - coupling_beta * coupling_difference / coupling_alpha
            )
            candidate = _phase_candidate(
                source="interior stationary point",
                branch=branch,
                value=critical_value,
                cos_alpha=stationary_cos_alpha,
                cos_beta=cos_beta,
                cos_beta_minus_alpha=cos_difference,
            )
            if branch == "lower":
                lower_candidates.append(candidate)
            else:
                upper_candidates.append(candidate)

    minimum = _select_candidate(lower_candidates, "minimum")
    maximum = _select_candidate(upper_candidates, "maximum")
    if q15_sign(maximum["value"] - minimum["value"]) < 0:
        raise AssertionError("computed range is reversed")

    def serialize_candidate(candidate: dict) -> dict:
        return {
            "source": candidate["source"],
            "branch": candidate["branch"],
            "value": exact_string(candidate["value"]),
            "cos_alpha": exact_string(candidate["cos_alpha"]),
            "cos_beta": exact_string(candidate["cos_beta"]),
            "cos_beta_minus_alpha": exact_string(
                candidate["cos_beta_minus_alpha"]
            ),
        }

    return {
        "carrier_matrix": [
            [exact_string(hamiltonian[row, column]) for column in range(3)]
            for row in range(3)
        ],
        "fiber_expectation": (
            f"{exact_string(offset)} + {exact_string(coupling_alpha)} cos(alpha) "
            f"+ {exact_string(coupling_beta)} cos(beta) + "
            f"{exact_string(coupling_difference)} cos(beta-alpha)"
        ),
        "offset": exact_string(offset),
        "couplings": {
            "A_cos_alpha": exact_string(coupling_alpha),
            "B_cos_beta": exact_string(coupling_beta),
            "C_cos_beta_minus_alpha": exact_string(coupling_difference),
        },
        "candidate_reduction": {
            "minimum_candidates": [
                serialize_candidate(candidate) for candidate in lower_candidates
            ],
            "maximum_candidates": [
                serialize_candidate(candidate) for candidate in upper_candidates
            ],
            "candidate_count": len(lower_candidates) + len(upper_candidates),
        },
        "minimum": exact_string(minimum["value"]),
        "maximum": exact_string(maximum["value"]),
        "width": exact_string(maximum["value"] - minimum["value"]),
        "minimum_numeric": float(sp.N(minimum["value"], 18)),
        "maximum_numeric": float(sp.N(maximum["value"], 18)),
        "minimum_phase_cosines": {
            "cos_alpha": exact_string(minimum["cos_alpha"]),
            "cos_beta": exact_string(minimum["cos_beta"]),
            "cos_beta_minus_alpha": exact_string(
                minimum["cos_beta_minus_alpha"]
            ),
        },
        "maximum_phase_cosines": {
            "cos_alpha": exact_string(maximum["cos_alpha"]),
            "cos_beta": exact_string(maximum["cos_beta"]),
            "cos_beta_minus_alpha": exact_string(
                maximum["cos_beta_minus_alpha"]
            ),
        },
        "minimum_source": minimum["source"],
        "maximum_source": maximum["source"],
    }


def build() -> dict:
    gamma_diagonal = diagonal_1rdm(BD_SUPPORT, BD_WEIGHTS, 6)
    recovered, incidence_certificate = recover_weights(BD_SUPPORT, gamma_diagonal)
    selection = borland_dennis_selection_certificate()
    witnesses = complete_collision_free_witnesses(BD_SUPPORT, order=2)

    coherence_range = real_symmetric_triangle_range(COHERENCE_TRIANGLE_MATRIX)
    pi_flux_range = real_symmetric_triangle_range(PI_FLUX_TRIANGLE_MATRIX)
    polygon_range = complete_triangle_range(BD_WEIGHTS)
    coherence_decomposition = carrier_operator_decomposition(
        COHERENCE_TRIANGLE_MATRIX, witnesses
    )
    pi_flux_decomposition = carrier_operator_decomposition(
        PI_FLUX_TRIANGLE_MATRIX, witnesses
    )

    if coherence_range["minimum"] != polygon_range["minimum"]:
        raise AssertionError("critical-value and polygon minima disagree")
    if coherence_range["maximum"] != polygon_range["maximum"]:
        raise AssertionError("critical-value and polygon maxima disagree")

    # Exact generalized-Pauli checks in one-based lambda notation.
    pair_equalities = (
        gamma_diagonal[0] + gamma_diagonal[5],
        gamma_diagonal[1] + gamma_diagonal[4],
        gamma_diagonal[2] + gamma_diagonal[3],
    )
    gpc_slack = 2 - (gamma_diagonal[0] + gamma_diagonal[1] + gamma_diagonal[3])
    strictly_decreasing = all(
        gamma_diagonal[index] > gamma_diagonal[index + 1]
        for index in range(len(gamma_diagonal) - 1)
    )

    return {
        "schema_version": 2,
        "generated_by": GENERATOR,
        "project": (
            "certified two-body observable ranges on fermionic one-body "
            "moment-map fibers"
        ),
        "evidence": {
            "status": "EXACT",
            "arithmetic": "integer, rational, and quadratic-radical identities",
            "scope": (
                "Pure states in wedge^3 C^6 with the stated complete 1-RDM. "
                "Strictly nondegenerate occupations fix the natural-orbital frame; "
                "the Borland-Dennis equalities and saturated GPC reduce the global "
                "fiber to the certified three-determinant carrier. The general "
                "range engine covers rational real-symmetric carrier matrices."
            ),
        },
        "system": "wedge^3 C^6",
        "target_1rdm": {
            "matrix_form": "diagonal in the fixed natural-orbital frame",
            "diagonal": [exact_string(value) for value in gamma_diagonal],
            "trace": exact_string(sum(gamma_diagonal)),
            "strictly_decreasing": strictly_decreasing,
            "borland_dennis_pair_sums": [
                exact_string(value) for value in pair_equalities
            ],
            "saturated_gpc": "2 - (lambda_1 + lambda_2 + lambda_4) = 0",
            "gpc_slack": exact_string(gpc_slack),
        },
        "selection_rule_certificate": selection,
        "fiber": {
            "support_dets": [list(row) for row in BD_SUPPORT],
            "support_size": len(BD_SUPPORT),
            "one_hop_pairs": [list(edge) for edge in one_hop_pairs(BD_SUPPORT)],
            "one_hop_free": not one_hop_pairs(BD_SUPPORT),
            "singleton_incidence": incidence_certificate,
            "weights": [exact_string(weight) for weight in BD_WEIGHTS],
            "weights_recovered_from_target": [
                exact_string(value) for value in recovered
            ],
            "weights_unique": bool(recovered == BD_WEIGHTS),
            "parameterization": (
                "sqrt(3/5)|012> + exp(i alpha)/2 |034> "
                "+ sqrt(3/20) exp(i beta)|135>"
            ),
            "projective_real_dimension": 2,
        },
        "observable_realization": {
            "pair_annihilator_convention": "A_(p,q) = a_q a_p for p < q",
            "collision_free_witnesses": witnesses,
            "diagonal_realization": (
                "The invertible singleton-incidence minor on orbitals 0,1,2 "
                "realizes every carrier diagonal by a one-body operator."
            ),
            "off_diagonal_realization": (
                "The three witnesses realize every real carrier off-diagonal "
                "by signed Hermitian two-body pair hopping."
            ),
        },
        "range_engine": {
            "scope": "all rational real-symmetric 3 by 3 carrier matrices",
            "normal_form": "E0 + A cos(alpha) + B cos(beta) + C cos(beta-alpha)",
            "elimination": (
                "At t=cos(alpha), optimize beta exactly to obtain "
                "E0+A*t +/- sqrt(B^2+C^2+2BC*t)."
            ),
            "global_certificate": (
                "The plus branch is concave and the minus branch is convex on "
                "[-1,1]. Endpoints and the unique feasible stationary candidate "
                "therefore contain both global extrema."
            ),
            "maximum_candidate_count": 5,
            "coefficient_field_for_this_target": "Q(sqrt(15))",
        },
        "examples": {
            "unfrustrated_coherence_triangle": {
                "name": "complete collision-free two-body coherence triangle",
                "phasor_identity": (
                    "<W> = |sqrt(3/5) + exp(i alpha)/2 "
                    "+ sqrt(3/20) exp(i beta)|^2 - 1"
                ),
                "operator_decomposition": coherence_decomposition,
                "range": coherence_range,
                "independent_polygon_certificate": polygon_range,
            },
            "pi_flux_coherence_triangle": {
                "name": "signed coherence triangle with negative cycle product",
                "interpretation": (
                    "One negative carrier edge gives a gauge-invariant pi flux. "
                    "The interior critical point now supplies the maximum."
                ),
                "operator_decomposition": pi_flux_decomposition,
                "range": pi_flux_range,
            },
        },
        "theorems": [
            {
                "name": (
                    "complete-coherence polygon range on a fixed-support "
                    "one-body fiber"
                ),
                "hypotheses": [
                    "The determinant carrier is one-hop-free, so its 1-RDM is diagonal for every phase choice.",
                    "Singleton incidence has full column rank, so the fixed diagonal 1-RDM uniquely fixes positive weights w_a.",
                    "Every determinant pair has a collision-free two-body channel, with its fermionic sign absorbed into the observable.",
                ],
                "conclusion": (
                    "For W equal to the sum of all signed Hermitian coherence "
                    "channels, the exact range is "
                    "[(max(0, 2 r_max - sum_a r_a))^2 - 1, "
                    "(sum_a r_a)^2 - 1], where r_a=sqrt(w_a)."
                ),
                "benchmark_specialization": (
                    "The Borland-Dennis radii close a triangle, hence the minimum "
                    "is -1 and the maximum is 3(2+sqrt(15))/10."
                ),
            },
            {
                "name": "exact real-symmetric triangle constrained search",
                "hypotheses": [
                    "The fixed-1RDM fiber has the certified three-determinant phase-torus form.",
                    "The restricted carrier observable is real symmetric with rational entries.",
                ],
                "conclusion": (
                    "Both global extrema are selected exactly from the four "
                    "alpha=0,pi branch endpoints and at most one feasible interior "
                    "stationary value."
                ),
            },
        ],
        "nonclaims": [
            "Complex Hermitian carrier matrices with nontrivial cycle phase are not yet covered by the closed-form engine.",
            "No mixed-state constrained search is addressed.",
            "No claim is made that every fixed-1RDM fiber admits a collision-free complete graph.",
            "The general rank-(3,6) functional remains to be assembled from all strata and observables.",
        ],
        "next_steps": [
            "Extend exact elimination to rational complex Hermitian carrier matrices, whose cycle phase is gauge invariant.",
            "Enumerate all Borland-Dennis boundary strata and certify their fixed-1RDM fiber dimensions.",
            "Compare exact support-restricted extrema with unrestricted numerical constrained searches away from exact pinning.",
            "Add FCIDUMP/OpenFermion import and project physical two-body Hamiltonians onto certified fibers.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    payload = build()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.out.exists() or args.out.read_text() != text:
            raise SystemExit(f"artifact is stale: {args.out}")
        print(f"verified {args.out}")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    first = payload["examples"]["unfrustrated_coherence_triangle"]["range"]
    second = payload["examples"]["pi_flux_coherence_triangle"]["range"]
    print(
        f"wrote {args.out}: unfrustrated [{first['minimum']}, {first['maximum']}], "
        f"pi-flux [{second['minimum']}, {second['maximum']}]"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
