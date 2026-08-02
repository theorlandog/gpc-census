#!/usr/bin/env python3
"""Independently verify the exact v_B carrier-Hamiltonian artifact."""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import sympy as sp


SUPPORT = [
    (0, 1, 2, 4),
    (0, 1, 2, 8),
    (0, 1, 3, 5),
    (0, 2, 3, 6),
    (0, 3, 4, 7),
    (0, 3, 7, 8),
    (1, 3, 6, 8),
    (2, 3, 4, 8),
]
PAIRS = list(itertools.combinations(range(9), 2))
BASIS = list(itertools.combinations(range(9), 4))
PAIR_INCIDENCE_PIVOTS = [
    (0, 1),
    (0, 4),
    (0, 3),
    (0, 2),
    (0, 5),
    (1, 4),
    (1, 3),
    (2, 3),
]


def sign(sequence):
    inversions = sum(
        sequence[left] > sequence[right]
        for left in range(len(sequence))
        for right in range(left + 1, len(sequence))
    )
    return -1 if inversions % 2 else 1


def formal_star(expression):
    return sp.expand(expression).xreplace({sp.I: -sp.I})


def one_rdm(amplitudes):
    """Independently contract a pure carrier state to its 1-RDM."""
    support_pos = {determinant: index for index, determinant in enumerate(SUPPORT)}
    matrix = sp.zeros(9)
    for source_index, determinant in enumerate(SUPPORT):
        for removed_position, removed in enumerate(determinant):
            remainder = [mode for mode in determinant if mode != removed]
            for added in range(9):
                if added in remainder:
                    continue
                target = tuple(sorted(remainder + [added]))
                target_index = support_pos.get(target)
                if target_index is None:
                    continue
                matrix[added, removed] += (
                    formal_star(amplitudes[target_index])
                    * amplitudes[source_index]
                    * (-1) ** removed_position
                    * (-1) ** target.index(added)
                )
    return matrix


def two_rdm(amplitudes):
    """Independently contract a pure carrier state to its 2-RDM."""
    pair_pos = {pair: index for index, pair in enumerate(PAIRS)}
    vectors = {}
    for determinant, amplitude in zip(SUPPORT, amplitudes):
        for pair in itertools.combinations(determinant, 2):
            remainder = tuple(mode for mode in determinant if mode not in pair)
            vectors[(pair_pos[pair], pair_pos[remainder])] = sign(
                pair + remainder
            ) * amplitude
    matrix = sp.MutableSparseMatrix(len(PAIRS), len(PAIRS), {})
    for left in range(len(PAIRS)):
        for right in range(len(PAIRS)):
            value = sum(
                formal_star(vectors[(left, remainder)])
                * vectors[(right, remainder)]
                for remainder in range(len(PAIRS))
                if (left, remainder) in vectors and (right, remainder) in vectors
            )
            if value != 0:
                matrix[left, right] = value
    return matrix


def trace_cube(matrix):
    rows = {row: [] for row in range(matrix.rows)}
    for (row, column), value in matrix.todok().items():
        rows[row].append((column, value))
    result = sp.S.Zero
    for first in range(matrix.rows):
        for second, left in rows[first]:
            for third, middle in rows[second]:
                right = matrix[third, first]
                if right != 0:
                    result += left * middle * right
    return result


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        "results/data/vb_carrier_hamiltonian.json"
    )
    artifact = json.loads(path.read_text())
    assert artifact["generated_by"] == "scripts/vb_carrier_hamiltonian.py"
    assert artifact["evidence"]["status"] == "EXACT"
    assert artifact["carrier"]["support"] == [list(row) for row in SUPPORT]
    assert artifact["checks"] == {
        "canonical_basis_dimension": 126,
        "carrier_dimension": 8,
        "carrier_leakage_entries": 0,
        "carrier_pair_kernel_symbolic_support_entries": 7,
        "carrier_restricted_nonzero_entries": 5,
        "carrier_schrodinger_residual_entries": 0,
        "common_rho_leakage_entries": 0,
        "common_rho_matrix_residual_entries": 0,
        "common_rho_schrodinger_residual_entries": 0,
        "conic_reductions_exact": True,
        "gamma2_cubic_invariant_identity": True,
        "gamma2_pair_incidence_minor_determinant": 1,
        "gamma2_pair_diagonal_direction_nonzero_entries": 8,
        "gamma2_reconstruction_tree_connected": True,
        "gamma2_reconstruction_tree_edges": 7,
        "gamma2_trace_identity": True,
        "physical_radicals_uniformly_positive": True,
        "wall_regular": True,
    }
    assert artifact["conic_clock"]["equation"] == "s^2=119+70*t-85*t^2"
    assert artifact["census_record"] == {
        "system": "(4,9)",
        "index": 65,
        "classified": "INTERFERENCE",
        "observability_artifact": "results/data/rdm_observability_census.json",
    }
    assert artifact["conic_clock"]["t_dot"] == "-s/sqrt(85)"
    assert artifact["conic_clock"]["s_dot"] == "sqrt(85)*(t-7/17)"
    assert artifact["carrier_state"]["phase"] == (
        "e=(3+7*t-2*t^2+I*s)/(2*sqrt((1+t)*(8-t)*(2-t)*(2+t)))"
    )
    assert artifact["common_full_one_rdm_lift"]["one_rdm"] == (
        "(diag(20,14,14,14,5,4,4,4,13)+3*(|8><4|+|4><8|))/23"
    )

    t, s = sp.symbols("t s", real=True)
    wall = 119 + 70 * t - 85 * t**2
    sqrt85 = sp.sqrt(85)
    r0, r1 = sp.sqrt(1 + t), sp.sqrt(8 - t)
    r4, r5 = sp.sqrt(2 - t), sp.sqrt(2 + t)
    product = r0**2 * r1**2 * r4**2 * r5**2
    p = 3 + 7 * t - 2 * t**2
    td = -s / sqrt85
    sd = sqrt85 * (t - sp.Rational(7, 17))
    e = (p + sp.I * s) / (2 * r0 * r1 * r4 * r5)
    eb = (p - sp.I * s) / (2 * r0 * r1 * r4 * r5)
    omega = (85 * t**3 - 105 * t**2 + 12 * t + 364) / (
        2 * sqrt85 * product
    )
    a = td / (2 * r0 * r1)
    b = -td / (2 * r4 * r5)
    cap_a = -sp.I * a * eb
    cap_b = sp.I * b
    alpha = {
        0: (cap_a + cap_b) / 2,
        1: (cap_a - cap_b) / 2,
        3: -(cap_a - cap_b) / 2,
    }
    u = (7 * t - 5 + sp.I * s) / (6 * r4 * r5)
    ub = (7 * t - 5 - sp.I * s) / (6 * r4 * r5)
    nu = (28 - 5 * t) / (sqrt85 * (4 - t**2))

    def remainder(expression):
        numerator, _ = sp.fraction(sp.factor(sp.together(expression)))
        return sp.factor(
            sp.Poly(sp.expand(numerator), s, domain="EX")
            .rem(sp.Poly(s**2 - wall, s, domain="EX"))
            .as_expr()
        )

    def derivative(expression):
        return sp.diff(expression, t) * td + sp.diff(expression, s) * sd

    assert remainder(2 * s * sd - sp.diff(wall, t) * td) == 0
    assert remainder(e * eb - 1) == 0
    assert remainder(derivative(e) - sp.I * omega * e) == 0
    assert remainder(u * ub - 1) == 0
    assert remainder(derivative(u) - sp.I * nu * u) == 0

    kernel = {((1, 4), (1, 4)): -omega}
    for mode, coefficient in alpha.items():
        left, right = (mode, 8), (mode, 4)
        kernel[(left, right)] = coefficient
        kernel[(right, left)] = formal_star(coefficient)
    assert len(kernel) == 7
    assert artifact["carrier_hamiltonian"]["pair_kernel_symbolic_support"] == [
        {"left": [1, 4], "right": [1, 4], "coefficient": "-omega"},
        {"left": [0, 8], "right": [0, 4], "coefficient": "alpha_0"},
        {
            "left": [0, 4],
            "right": [0, 8],
            "coefficient": "conjugate(alpha_0)",
        },
        {"left": [1, 8], "right": [1, 4], "coefficient": "alpha_1"},
        {
            "left": [1, 4],
            "right": [1, 8],
            "coefficient": "conjugate(alpha_1)",
        },
        {"left": [3, 8], "right": [3, 4], "coefficient": "alpha_3"},
        {
            "left": [3, 4],
            "right": [3, 8],
            "coefficient": "conjugate(alpha_3)",
        },
    ]

    basis_pos = {determinant: index for index, determinant in enumerate(BASIS)}
    operator = {}
    for (left, right), coefficient in kernel.items():
        for rest in PAIRS:
            if set(left).intersection(rest) or set(right).intersection(rest):
                continue
            target = tuple(sorted(left + rest))
            source = tuple(sorted(right + rest))
            key = (basis_pos[target], basis_pos[source])
            operator[key] = operator.get(key, 0) + sign(left + rest) * sign(
                right + rest
            ) * coefficient

    carrier = {basis_pos[row] for row in SUPPORT}
    for (row, column), value in operator.items():
        if (row in carrier) != (column in carrier):
            assert remainder(value) == 0
    for (row, column), value in operator.items():
        assert remainder(value - formal_star(operator.get((column, row), 0))) == 0

    support_pos = [basis_pos[row] for row in SUPPORT]
    restricted = sp.Matrix(
        8,
        8,
        lambda row, column: operator.get((support_pos[row], support_pos[column]), 0),
    )
    amplitudes = sp.Matrix(
        [
            r0 * e / sp.sqrt(23),
            r1 / sp.sqrt(23),
            2 / sp.sqrt(23),
            sp.sqrt(3) / sp.sqrt(23),
            r4 / sp.sqrt(23),
            r5 / sp.sqrt(23),
            1 / sp.sqrt(23),
            sp.sqrt(2) / sp.sqrt(23),
        ]
    )
    for value in restricted * amplitudes - sp.I * amplitudes.applyfunc(derivative):
        assert remainder(value) == 0

    nonzero_restricted = sum(
        1
        for row in range(8)
        for column in range(8)
        if restricted[row, column] != 0
        and remainder(restricted[row, column]) != 0
    )
    assert nonzero_restricted == 5

    lifted_kernel = {
        key: (value if key[0] == key[1] else (u * value if key[0][1] == 8 else ub * value))
        for key, value in kernel.items()
    }
    lifted_operator = {}
    for (left, right), coefficient in lifted_kernel.items():
        for rest in PAIRS:
            if set(left).intersection(rest) or set(right).intersection(rest):
                continue
            target = tuple(sorted(left + rest))
            source = tuple(sorted(right + rest))
            key = (basis_pos[target], basis_pos[source])
            lifted_operator[key] = lifted_operator.get(key, 0) + sign(
                left + rest
            ) * sign(right + rest) * coefficient
    for index, determinant in enumerate(BASIS):
        if 8 in determinant:
            lifted_operator[(index, index)] = lifted_operator.get((index, index), 0) - nu
    for (row, column), value in lifted_operator.items():
        if (row in carrier) != (column in carrier):
            assert remainder(value) == 0

    lifted_restricted = sp.Matrix(
        8,
        8,
        lambda row, column: lifted_operator.get(
            (support_pos[row], support_pos[column]), 0
        ),
    )
    lifted_amplitudes = sp.Matrix(
        [
            value * (u if 8 in determinant else 1)
            for value, determinant in zip(amplitudes, SUPPORT, strict=True)
        ]
    )
    lifted_residual = (
        lifted_restricted * lifted_amplitudes
        - sp.I * lifted_amplitudes.applyfunc(derivative)
    )
    for value in lifted_residual:
        assert remainder(value) == 0

    rho = one_rdm(lifted_amplitudes)
    expected_rho = sp.diag(20, 14, 14, 14, 5, 4, 4, 4, 13) / 23
    expected_rho[8, 4] = sp.Rational(3, 23)
    expected_rho[4, 8] = sp.Rational(3, 23)
    for row in range(9):
        for column in range(9):
            assert remainder(rho[row, column] - expected_rho[row, column]) == 0

    gamma2 = two_rdm(amplitudes)
    gamma2_trace = sum(gamma2[index, index] for index in range(len(PAIRS)))
    assert remainder(gamma2_trace - 6) == 0
    gamma2_cubic = trace_cube(gamma2)
    gamma2_cubic_target = sp.Rational(6, 12167) * (2899 - 8 * t)
    assert remainder(gamma2_cubic - gamma2_cubic_target) == 0

    incidence = sp.Matrix(
        [
            [int(set(pair).issubset(determinant)) for determinant in SUPPORT]
            for pair in PAIR_INCIDENCE_PIVOTS
        ]
    )
    assert incidence.det() == 1
    weight_direction = sp.Matrix([1, -1, 0, 0, -1, 1, 0, 0])
    full_pair_incidence = sp.Matrix(
        [
            [int(set(pair).issubset(determinant)) for determinant in SUPPORT]
            for pair in PAIRS
        ]
    )
    assert sum(value != 0 for value in full_pair_incidence * weight_direction) == 8
    closure = artifact["exact_td2rdm_consequence"]
    assert closure["pair_incidence_pivots"] == [
        list(pair) for pair in PAIR_INCIDENCE_PIVOTS
    ]
    assert closure["pair_incidence_minor_determinant"] == 1
    assert closure["basis_free_correlation_witness"] == (
        "Tr(Gamma2^3)=6*(2899-8*t)/12167 and "
        "d Tr(Gamma2^3)/dt=-48/12167"
    )
    assert closure["two_rdm_trace"] == "Tr(Gamma2)=6"
    tree = closure["gamma2_reconstruction_tree"]
    assert len(tree) == 7
    adjacency = {index: set() for index in range(len(SUPPORT))}
    support_index = {determinant: index for index, determinant in enumerate(SUPPORT)}
    for edge in tree:
        left_pair = tuple(edge["left_modes"])
        right_pair = tuple(edge["right_modes"])
        contributions = []
        for remainder_pair in PAIRS:
            if set(left_pair).intersection(remainder_pair):
                continue
            if set(right_pair).intersection(remainder_pair):
                continue
            left_det = tuple(sorted(left_pair + remainder_pair))
            right_det = tuple(sorted(right_pair + remainder_pair))
            if left_det not in support_index or right_det not in support_index:
                continue
            contributions.append(
                (
                    support_index[left_det],
                    support_index[right_det],
                    sign(left_pair + remainder_pair)
                    * sign(right_pair + remainder_pair),
                )
            )
        assert contributions == [(edge["left"], edge["right"], edge["coefficient"])]
        assert {edge["parent"], edge["child"]} == {edge["left"], edge["right"]}
        adjacency[edge["parent"]].add(edge["child"])
        adjacency[edge["child"]].add(edge["parent"])
    reached = {0}
    frontier = [0]
    while frontier:
        current = frontier.pop()
        for neighbor in adjacency[current] - reached:
            reached.add(neighbor)
            frontier.append(neighbor)
    assert reached == set(range(len(SUPPORT)))

    print(
        "verified exact v_B carrier Hamiltonian: "
        "7 pair entries, zero leakage, both Schrodinger identities, fixed 1-RDM"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
