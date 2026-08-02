#!/usr/bin/env python3
"""Generate the exact scheduled two-body driver for the known v_B circle.

The certificate uses the signed conic coordinate

    s**2 = 119 + 70*t - 85*t**2

and an everywhere regular periodic clock.  A seven-entry Hermitian pair
kernel preserves the eight-determinant carrier and solves the Schrodinger
equation exactly.  A diagonal orbital lift gives the corresponding driver on
the common-full-1RDM loop.

Usage:
  python scripts/vb_carrier_hamiltonian.py
  python scripts/vb_carrier_hamiltonian.py --check
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "data" / "vb_carrier_hamiltonian.json"
GENERATOR = "scripts/vb_carrier_hamiltonian.py"

VB_SUPPORT = [
    (0, 1, 2, 4),
    (0, 1, 2, 8),
    (0, 1, 3, 5),
    (0, 2, 3, 6),
    (0, 3, 4, 7),
    (0, 3, 7, 8),
    (1, 3, 6, 8),
    (2, 3, 4, 8),
]
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
GAMMA2_RECONSTRUCTION_TREE = [
    {
        "parent": 0,
        "child": 1,
        "left": 0,
        "right": 1,
        "left_modes": [1, 4],
        "right_modes": [1, 8],
        "coefficient": 1,
    },
    {
        "parent": 0,
        "child": 2,
        "left": 0,
        "right": 2,
        "left_modes": [2, 4],
        "right_modes": [3, 5],
        "coefficient": 1,
    },
    {
        "parent": 0,
        "child": 3,
        "left": 0,
        "right": 3,
        "left_modes": [1, 4],
        "right_modes": [3, 6],
        "coefficient": -1,
    },
    {
        "parent": 0,
        "child": 7,
        "left": 0,
        "right": 7,
        "left_modes": [0, 1],
        "right_modes": [3, 8],
        "coefficient": -1,
    },
    {
        "parent": 1,
        "child": 6,
        "left": 1,
        "right": 6,
        "left_modes": [0, 2],
        "right_modes": [3, 6],
        "coefficient": -1,
    },
    {
        "parent": 2,
        "child": 4,
        "left": 2,
        "right": 4,
        "left_modes": [1, 5],
        "right_modes": [4, 7],
        "coefficient": -1,
    },
    {
        "parent": 2,
        "child": 5,
        "left": 2,
        "right": 5,
        "left_modes": [1, 5],
        "right_modes": [7, 8],
        "coefficient": -1,
    },
]
PAIRS = list(itertools.combinations(range(9), 2))
BASIS = list(itertools.combinations(range(9), 4))


def permutation_sign(sequence: tuple[int, ...]) -> int:
    inversions = sum(
        sequence[left] > sequence[right]
        for left in range(len(sequence))
        for right in range(left + 1, len(sequence))
    )
    return -1 if inversions % 2 else 1


def exact_objects() -> dict[str, object]:
    t, s = sp.symbols("t s", real=True)
    sqrt85 = sp.sqrt(85)
    wall = 119 + 70 * t - 85 * t**2
    phase_numerator = 3 + 7 * t - 2 * t**2

    # Keeping the four positive physical radicals factored makes all branch
    # conventions explicit and prevents unsafe automatic radical merging.
    r0 = sp.sqrt(1 + t)
    r1 = sp.sqrt(8 - t)
    r4 = sp.sqrt(2 - t)
    r5 = sp.sqrt(2 + t)
    amplitude_product = r0**2 * r1**2 * r4**2 * r5**2

    t_dot = -s / sqrt85
    s_dot = sqrt85 * (t - sp.Rational(7, 17))
    phase = (phase_numerator + sp.I * s) / (2 * r0 * r1 * r4 * r5)
    phase_bar = (phase_numerator - sp.I * s) / (2 * r0 * r1 * r4 * r5)
    omega = (
        85 * t**3 - 105 * t**2 + 12 * t + 364
    ) / (2 * sqrt85 * amplitude_product)

    transfer_01 = t_dot / (2 * r0 * r1)
    transfer_45 = -t_dot / (2 * r4 * r5)
    channel_01 = -sp.I * transfer_01 * phase_bar
    channel_45 = sp.I * transfer_45
    alpha0 = (channel_01 + channel_45) / 2
    alpha1 = (channel_01 - channel_45) / 2
    alpha3 = -alpha1

    lift_phase = (7 * t - 5 + sp.I * s) / (6 * r4 * r5)
    lift_phase_bar = (7 * t - 5 - sp.I * s) / (6 * r4 * r5)
    lift_rate = (28 - 5 * t) / (sqrt85 * (4 - t**2))

    amplitudes = [
        r0 * phase / sp.sqrt(23),
        r1 / sp.sqrt(23),
        sp.Integer(2) / sp.sqrt(23),
        sp.sqrt(3) / sp.sqrt(23),
        r4 / sp.sqrt(23),
        r5 / sp.sqrt(23),
        sp.Integer(1) / sp.sqrt(23),
        sp.sqrt(2) / sp.sqrt(23),
    ]
    lifted_amplitudes = [
        amplitude * (lift_phase if 8 in determinant else 1)
        for amplitude, determinant in zip(amplitudes, VB_SUPPORT, strict=True)
    ]

    return {
        "t": t,
        "s": s,
        "wall": wall,
        "phase_numerator": phase_numerator,
        "radicals": (r0, r1, r4, r5),
        "amplitude_product": amplitude_product,
        "t_dot": t_dot,
        "s_dot": s_dot,
        "phase": phase,
        "phase_bar": phase_bar,
        "omega": omega,
        "transfer_01": transfer_01,
        "transfer_45": transfer_45,
        "channel_01": channel_01,
        "channel_45": channel_45,
        "alpha0": alpha0,
        "alpha1": alpha1,
        "alpha3": alpha3,
        "lift_phase": lift_phase,
        "lift_phase_bar": lift_phase_bar,
        "lift_rate": lift_rate,
        "amplitudes": amplitudes,
        "lifted_amplitudes": lifted_amplitudes,
    }


def star(expression):
    """Formal conjugation on the real t,s chart with fixed positive radicals."""
    return sp.expand(expression).xreplace({sp.I: -sp.I})


def conic_remainder(expression, objects) -> sp.Expr:
    """Reduce an exact rational expression modulo s**2-wall(t)."""
    s = objects["s"]
    wall = objects["wall"]
    numerator, _ = sp.fraction(sp.factor(sp.together(expression)))
    polynomial = sp.Poly(sp.expand(numerator), s, domain="EX")
    divisor = sp.Poly(s**2 - wall, s, domain="EX")
    return sp.factor(polynomial.rem(divisor).as_expr())


def assert_conic_zero(expression, objects, label: str) -> None:
    remainder = conic_remainder(expression, objects)
    if remainder != 0:
        raise AssertionError(f"{label}: nonzero conic remainder {remainder}")


def clock_derivative(expression, objects):
    return sp.diff(expression, objects["t"]) * objects["t_dot"] + sp.diff(
        expression, objects["s"]
    ) * objects["s_dot"]


def pair_kernel(objects, lifted: bool = False) -> dict:
    phase_multiplier = objects["lift_phase"] if lifted else sp.Integer(1)
    kernel = {((1, 4), (1, 4)): -objects["omega"]}
    for mode, coefficient in (
        (0, objects["alpha0"]),
        (1, objects["alpha1"]),
        (3, objects["alpha3"]),
    ):
        left = (mode, 8)
        right = (mode, 4)
        value = phase_multiplier * coefficient
        kernel[(left, right)] = value
        kernel[(right, left)] = star(value)
    return kernel


def four_particle_operator(kernel: dict) -> dict[tuple[int, int], sp.Expr]:
    """Second-quantize a pair kernel on wedge^4 C^9."""
    basis_pos = {determinant: index for index, determinant in enumerate(BASIS)}
    entries: dict[tuple[int, int], sp.Expr] = {}
    for (left_pair, right_pair), coefficient in kernel.items():
        for remainder in PAIRS:
            if set(left_pair).intersection(remainder):
                continue
            if set(right_pair).intersection(remainder):
                continue
            target = tuple(sorted(left_pair + remainder))
            source = tuple(sorted(right_pair + remainder))
            sign = permutation_sign(left_pair + remainder)
            sign *= permutation_sign(right_pair + remainder)
            key = (basis_pos[target], basis_pos[source])
            entries[key] = entries.get(key, sp.S.Zero) + sign * coefficient
    return {key: sp.factor(value) for key, value in entries.items() if value != 0}


def add_one_body_number(operator: dict, mode: int, coefficient) -> dict:
    result = dict(operator)
    for index, determinant in enumerate(BASIS):
        if mode in determinant:
            key = (index, index)
            result[key] = sp.factor(result.get(key, sp.S.Zero) + coefficient)
    return result


def carrier_matrix(operator: dict) -> sp.Matrix:
    basis_pos = {determinant: index for index, determinant in enumerate(BASIS)}
    support_pos = [basis_pos[determinant] for determinant in VB_SUPPORT]
    return sp.Matrix(
        8,
        8,
        lambda row, column: operator.get(
            (support_pos[row], support_pos[column]), sp.S.Zero
        ),
    )


def verify_carrier_preservation(operator: dict, objects, label: str) -> int:
    basis_pos = {determinant: index for index, determinant in enumerate(BASIS)}
    carrier = {basis_pos[determinant] for determinant in VB_SUPPORT}
    violations = 0
    for (row, column), value in operator.items():
        if (column in carrier) != (row in carrier):
            if conic_remainder(value, objects) != 0:
                violations += 1
    if violations:
        raise AssertionError(f"{label}: {violations} carrier leakage entries")
    return violations


def one_rdm(amplitudes) -> sp.Matrix:
    support_pos = {determinant: index for index, determinant in enumerate(VB_SUPPORT)}
    matrix = sp.zeros(9)
    for source_index, determinant in enumerate(VB_SUPPORT):
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
                    star(amplitudes[target_index])
                    * amplitudes[source_index]
                    * (-1) ** removed_position
                    * (-1) ** target.index(added)
                )
    return matrix


def two_rdm(amplitudes) -> sp.MutableSparseMatrix:
    """Contract a pure carrier state to the pair-indexed 2-RDM."""
    pair_pos = {pair: index for index, pair in enumerate(PAIRS)}
    vectors: dict[tuple[int, int], sp.Expr] = {}
    for determinant, amplitude in zip(VB_SUPPORT, amplitudes, strict=True):
        for pair in itertools.combinations(determinant, 2):
            remainder = tuple(mode for mode in determinant if mode not in pair)
            key = (pair_pos[pair], pair_pos[remainder])
            vectors[key] = permutation_sign(pair + remainder) * amplitude
    matrix = sp.MutableSparseMatrix(len(PAIRS), len(PAIRS), {})
    for left in range(len(PAIRS)):
        for right in range(len(PAIRS)):
            value = sum(
                star(vectors[(left, remainder)]) * vectors[(right, remainder)]
                for remainder in range(len(PAIRS))
                if (left, remainder) in vectors and (right, remainder) in vectors
            )
            if value != 0:
                matrix[left, right] = value
    return matrix


def sparse_trace_cube(matrix: sp.MutableSparseMatrix) -> sp.Expr:
    rows: dict[int, list[tuple[int, sp.Expr]]] = {
        row: [] for row in range(matrix.rows)
    }
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


def verify_gamma2_reconstruction_certificate() -> dict:
    """Verify the carrier's rational pure-state reconstruction certificate."""
    incidence = sp.Matrix(
        [
            [int(set(pair).issubset(determinant)) for determinant in VB_SUPPORT]
            for pair in PAIR_INCIDENCE_PIVOTS
        ]
    )
    determinant = incidence.det()
    if determinant != 1:
        raise AssertionError(f"pair-incidence minor determinant is {determinant}")
    weight_direction = sp.Matrix([1, -1, 0, 0, -1, 1, 0, 0])
    full_pair_incidence = sp.Matrix(
        [
            [int(set(pair).issubset(determinant)) for determinant in VB_SUPPORT]
            for pair in PAIRS
        ]
    )
    pair_diagonal_direction = full_pair_incidence * weight_direction
    nonzero_pair_directions = sum(value != 0 for value in pair_diagonal_direction)
    if not nonzero_pair_directions:
        raise AssertionError("the v_B tangent does not change pair occupations")

    support_pos = {determinant: index for index, determinant in enumerate(VB_SUPPORT)}
    adjacency = {index: set() for index in range(len(VB_SUPPORT))}
    for edge in GAMMA2_RECONSTRUCTION_TREE:
        left_pair = tuple(edge["left_modes"])
        right_pair = tuple(edge["right_modes"])
        contributions = []
        for remainder in PAIRS:
            if set(left_pair).intersection(remainder):
                continue
            if set(right_pair).intersection(remainder):
                continue
            left_det = tuple(sorted(left_pair + remainder))
            right_det = tuple(sorted(right_pair + remainder))
            if left_det not in support_pos or right_det not in support_pos:
                continue
            contributions.append(
                (
                    support_pos[left_det],
                    support_pos[right_det],
                    permutation_sign(left_pair + remainder)
                    * permutation_sign(right_pair + remainder),
                )
            )
        expected = [(edge["left"], edge["right"], edge["coefficient"])]
        if contributions != expected:
            raise AssertionError(
                f"Gamma2 tree channel {left_pair},{right_pair}: {contributions}"
            )
        parent, child = edge["parent"], edge["child"]
        adjacency[parent].add(child)
        adjacency[child].add(parent)

    reached = {0}
    frontier = [0]
    while frontier:
        current = frontier.pop()
        for neighbor in adjacency[current] - reached:
            reached.add(neighbor)
            frontier.append(neighbor)
    if len(reached) != len(VB_SUPPORT):
        raise AssertionError("Gamma2 reconstruction tree is disconnected")
    return {
        "gamma2_pair_incidence_minor_determinant": int(determinant),
        "gamma2_pair_diagonal_direction_nonzero_entries": nonzero_pair_directions,
        "gamma2_reconstruction_tree_edges": len(GAMMA2_RECONSTRUCTION_TREE),
        "gamma2_reconstruction_tree_connected": True,
    }


def verify_exact_identities(objects) -> dict:
    t = objects["t"]
    s = objects["s"]
    wall = objects["wall"]
    phase = objects["phase"]
    phase_bar = objects["phase_bar"]
    lift_phase = objects["lift_phase"]
    lift_phase_bar = objects["lift_phase_bar"]

    assert_conic_zero(
        2 * s * objects["s_dot"] - sp.diff(wall, t) * objects["t_dot"],
        objects,
        "clock tangency",
    )
    assert_conic_zero(phase * phase_bar - 1, objects, "carrier phase norm")
    assert_conic_zero(
        clock_derivative(phase, objects) - sp.I * objects["omega"] * phase,
        objects,
        "carrier phase derivative",
    )
    assert_conic_zero(
        lift_phase * lift_phase_bar - 1, objects, "lift phase norm"
    )
    assert_conic_zero(
        clock_derivative(lift_phase, objects)
        - sp.I * objects["lift_rate"] * lift_phase,
        objects,
        "lift phase derivative",
    )

    state_norm = sum(
        star(amplitude) * amplitude for amplitude in objects["amplitudes"]
    )
    assert_conic_zero(state_norm - 1, objects, "carrier state norm")

    kernel = pair_kernel(objects)
    if len(kernel) != 7:
        raise AssertionError("the carrier pair kernel does not have seven entries")
    operator = four_particle_operator(kernel)
    leakage = verify_carrier_preservation(operator, objects, "carrier driver")
    restricted = carrier_matrix(operator)

    expected = sp.zeros(8)
    expected[0, 0] = -objects["omega"]
    expected[0, 1] = sp.I * objects["transfer_01"] * phase
    expected[1, 0] = objects["channel_01"]
    expected[4, 5] = sp.I * objects["transfer_45"]
    expected[5, 4] = -sp.I * objects["transfer_45"]
    for row in range(8):
        for column in range(8):
            assert_conic_zero(
                restricted[row, column] - expected[row, column],
                objects,
                f"carrier restriction ({row},{column})",
            )

    amplitudes = sp.Matrix(objects["amplitudes"])
    derivatives = sp.Matrix(
        [clock_derivative(amplitude, objects) for amplitude in amplitudes]
    )
    residual = restricted * amplitudes - sp.I * derivatives
    for row, value in enumerate(residual):
        assert_conic_zero(value, objects, f"carrier Schrodinger row {row}")

    lifted_kernel = pair_kernel(objects, lifted=True)
    lifted_operator = add_one_body_number(
        four_particle_operator(lifted_kernel), 8, -objects["lift_rate"]
    )
    lifted_leakage = verify_carrier_preservation(
        lifted_operator, objects, "common-rho driver"
    )
    lifted_restricted = carrier_matrix(lifted_operator)
    lifted_amplitudes = sp.Matrix(objects["lifted_amplitudes"])
    lifted_derivatives = sp.Matrix(
        [clock_derivative(amplitude, objects) for amplitude in lifted_amplitudes]
    )
    lifted_residual = lifted_restricted * lifted_amplitudes - sp.I * lifted_derivatives
    for row, value in enumerate(lifted_residual):
        assert_conic_zero(value, objects, f"common-rho Schrodinger row {row}")

    rho = one_rdm(objects["lifted_amplitudes"])
    expected_rho = sp.diag(20, 14, 14, 14, 5, 4, 4, 4, 13) / 23
    expected_rho[8, 4] = sp.Rational(3, 23)
    expected_rho[4, 8] = sp.Rational(3, 23)
    for row in range(9):
        for column in range(9):
            assert_conic_zero(
                rho[row, column] - expected_rho[row, column],
                objects,
                f"common-rho entry ({row},{column})",
            )

    gamma2 = two_rdm(objects["amplitudes"])
    assert_conic_zero(
        sum(gamma2[index, index] for index in range(len(PAIRS))) - 6,
        objects,
        "Gamma2 trace normalization",
    )
    gamma2_cubic = sparse_trace_cube(gamma2)
    gamma2_cubic_target = sp.Rational(6, 12167) * (2899 - 8 * t)
    assert_conic_zero(
        gamma2_cubic - gamma2_cubic_target,
        objects,
        "Gamma2 cubic invariant",
    )

    midpoint = sp.Rational(7, 17)
    radius = 18 * sp.sqrt(35) / 85
    signed_radius = 18 * sp.sqrt(119) / 17
    t_minus = midpoint - radius
    t_plus = midpoint + radius
    assert_conic_zero(
        wall - 85 * (t_plus - t) * (t - t_minus),
        objects,
        "wall factorization",
    )
    if sp.factor(signed_radius**2 - 85 * radius**2) != 0:
        raise AssertionError("periodic clock radii are inconsistent")
    positive_factors = (
        1 + t_minus,
        8 - t_plus,
        2 - t_plus,
        2 + t_minus,
    )
    if not all(sp.ask(sp.Q.positive(value)) is True for value in positive_factors):
        raise AssertionError("a physical radical is not uniformly positive")
    if not (
        sp.ask(sp.Q.negative(objects["phase_numerator"].subs(t, t_minus))) is True
        and sp.ask(sp.Q.positive(objects["phase_numerator"].subs(t, t_plus)))
        is True
    ):
        raise AssertionError("the wall phase signs are incorrect")

    checks = {
        "canonical_basis_dimension": len(BASIS),
        "carrier_dimension": len(VB_SUPPORT),
        "carrier_pair_kernel_symbolic_support_entries": len(kernel),
        "carrier_restricted_nonzero_entries": sum(
            conic_remainder(value, objects) != 0 for value in restricted
        ),
        "carrier_leakage_entries": leakage,
        "common_rho_leakage_entries": lifted_leakage,
        "carrier_schrodinger_residual_entries": 0,
        "common_rho_schrodinger_residual_entries": 0,
        "common_rho_matrix_residual_entries": 0,
        "gamma2_cubic_invariant_identity": True,
        "gamma2_trace_identity": True,
        "conic_reductions_exact": True,
        "physical_radicals_uniformly_positive": True,
        "wall_regular": True,
    }
    checks.update(verify_gamma2_reconstruction_certificate())
    return checks


def build() -> dict:
    objects = exact_objects()
    checks = verify_exact_identities(objects)
    return {
        "schema_version": 1,
        "generated_by": GENERATOR,
        "evidence": {
            "status": "EXACT",
            "method": (
                "SymPy exact arithmetic, signed determinant action, and polynomial "
                "reduction modulo s^2-(119+70*t-85*t^2)"
            ),
            "scope": (
                "The known v_B circle on its eight-determinant carrier in "
                "wedge^4 C^9, plus its specified common-full-1RDM orbital lift."
            ),
        },
        "conventions": {
            "orbital_labels": "zero-based",
            "pair_annihilator": "A_(i,j)=a_j*a_i for i<j",
            "pair_operator": "H_K=sum_(p,q) K[p,q]*A_p^dagger*A_q",
            "two_rdm_normalization": "Gamma2[p,q]=<A_p^dagger*A_q>, Tr(Gamma2)=6",
            "schrodinger_equation": "I*d psi/d phi=H(phi)*psi(phi)",
            "clock_orientation": "phi increases from 0 to 2*pi",
            "radical_branches": (
                "sqrt(1+t), sqrt(8-t), sqrt(2-t), and sqrt(2+t) are the "
                "positive roots on the physical interval"
            ),
        },
        "carrier": {
            "support": [list(determinant) for determinant in VB_SUPPORT],
            "weights": [
                "(1+t)/23",
                "(8-t)/23",
                "4/23",
                "3/23",
                "(2-t)/23",
                "(2+t)/23",
                "1/23",
                "2/23",
            ],
            "phase_carrier": 0,
        },
        "census_record": {
            "system": "(4,9)",
            "index": 65,
            "classified": "INTERFERENCE",
            "observability_artifact": "results/data/rdm_observability_census.json",
        },
        "conic_clock": {
            "equation": "s^2=119+70*t-85*t^2",
            "t_minus": "7/17-18*sqrt(35)/85",
            "t_plus": "7/17+18*sqrt(35)/85",
            "t_of_phi": "7/17+(18*sqrt(35)/85)*cos(phi)",
            "s_of_phi": "(18*sqrt(119)/17)*sin(phi)",
            "t_dot": "-s/sqrt(85)",
            "s_dot": "sqrt(85)*(t-7/17)",
            "wall_behavior": (
                "At both walls s=0, the hopping coefficients vanish while the "
                "diagonal phase rates remain finite."
            ),
        },
        "carrier_state": {
            "phase": (
                "e=(3+7*t-2*t^2+I*s)/(2*sqrt((1+t)*(8-t)*(2-t)*(2+t)))"
            ),
            "phase_rate": (
                "omega=(85*t^3-105*t^2+12*t+364)/"
                "(2*sqrt(85)*(1+t)*(8-t)*(2-t)*(2+t))"
            ),
            "phase_identity": "d e/d phi=I*omega*e",
        },
        "carrier_hamiltonian": {
            "operator": (
                "H_car=-omega*n_1*n_4+T+T^dagger, "
                "T=a_8^dagger*a_4*(alpha_0*n_0+alpha_1*n_1+alpha_3*n_3)"
            ),
            "definitions": {
                "a": "-s/(2*sqrt(85)*sqrt((1+t)*(8-t)))",
                "b": "s/(2*sqrt(85)*sqrt((2-t)*(2+t)))",
                "A": "-I*a*conjugate(e)",
                "B": "I*b",
                "alpha_0": "(A+B)/2",
                "alpha_1": "(A-B)/2",
                "alpha_3": "-alpha_1",
            },
            "pair_kernel_symbolic_support": [
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
            ],
            "carrier_action": {
                "D0_to_D1": "alpha_0+alpha_1=A",
                "D4_to_D5": "-(alpha_0+alpha_3)=-B",
                "sole_leakage_channel": (
                    "(1,3,4,6) to D6 has coefficient "
                    "-(alpha_1+alpha_3)=0"
                ),
            },
            "exact_statement": (
                "H_car preserves the complete carrier and satisfies "
                "I*d psi/d phi=H_car*psi exactly."
            ),
        },
        "common_full_one_rdm_lift": {
            "orbital_unitary": "U=exp(I*beta*n_8)",
            "u": "exp(I*beta)=(7*t-5+I*s)/(6*sqrt((2-t)*(2+t)))",
            "beta_rate": "nu=(28-5*t)/(sqrt(85)*(4-t^2))",
            "hamiltonian": "H_common=U*H_car*U^dagger-nu*n_8",
            "pair_kernel_rule": (
                "Multiply alpha_0, alpha_1, alpha_3 by u and keep their "
                "Hermitian conjugates."
            ),
            "canonical_pair_only_identity": "n_8=(1/3)*sum_(j=0)^7 n_j*n_8 on N=4",
            "one_rdm": (
                "(diag(20,14,14,14,5,4,4,4,13)+"
                "3*(|8><4|+|4><8|))/23"
            ),
            "exact_statement": (
                "The lifted state solves the scheduled Schrodinger equation "
                "and has the displayed complete 1-RDM for every phi."
            ),
        },
        "counterdiabatic_consequence": {
            "statement": (
                "If an at-most-two-body H_0(phi) obeys H_0(phi)*psi(phi)=0, "
                "then for phi_dot=v(tau), "
                "H_total=H_0(phi)+v(tau)*H_car(phi) exactly tracks the carrier "
                "state. Replace H_car by H_common for the common-full-1RDM lift."
            ),
            "body_order": "at most two-body under the stated hypothesis on H_0",
        },
        "exact_td2rdm_consequence": {
            "domain": (
                "Pure full-support states on the fixed eight-determinant v_B "
                "carrier and in the displayed orbital frame."
            ),
            "pair_incidence_pivots": [
                list(pair) for pair in PAIR_INCIDENCE_PIVOTS
            ],
            "pair_incidence_minor_determinant": 1,
            "gamma2_reconstruction_tree": GAMMA2_RECONSTRUCTION_TREE,
            "higher_rdm_map": "Gamma3=F_(S,3)(Gamma2), rational on this chart",
            "correlation_change": (
                "The weight tangent (1,-1,0,0,-1,1,0,0) has a nonzero "
                "pair-diagonal image, so Gamma2 changes while the lifted "
                "complete 1-RDM remains fixed."
            ),
            "basis_free_correlation_witness": (
                "Tr(Gamma2^3)=6*(2899-8*t)/12167 and "
                "d Tr(Gamma2^3)/dt=-48/12167"
            ),
            "two_rdm_trace": "Tr(Gamma2)=6",
            "closed_equation": (
                "d Gamma2/d phi=L_(H(phi))(Gamma2,F_(S,3)(Gamma2))"
            ),
            "exact_statement": (
                "Because the pair-incidence minor is unimodular, the seven "
                "listed collision-free Gamma2 channels form a connected tree, "
                "and H_car and H_common preserve the carrier, both schedules "
                "give exact closed nonautonomous TD-2RDM trajectories."
            ),
        },
        "checks": checks,
        "nonclaims": [
            "The schedule is not an autonomous state-independent Hamiltonian.",
            "Neighboring initial states are not asserted to keep the same 1-RDM.",
            "No spatial locality, experimental noise tolerance, or driver gap is proved.",
            "The pair-only expression for n_8 is asserted only on the N=4 sector.",
            "No statement is made about the full fixed-1RDM fiber or its dimension.",
            "The rational TD-2RDM closure is not asserted for mixed states or outside the certified carrier chart.",
            "The package certifies closure but does not yet implement an amplitude-free Gamma2 propagator.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    artifact = build()
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.out.exists() or args.out.read_text() != rendered:
            raise SystemExit(f"stale artifact: {args.out}")
        print(f"verified {args.out}")
        return 0
    args.out.write_text(rendered)
    print(f"wrote {args.out}: exact seven-entry carrier pair kernel")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
