#!/usr/bin/env python3
"""Exact pair-coordinate and v_B symmetry certificates.

There are two separate results.

1. For a determinant support S, the diagonal two-body occupations are

       p_ij = sum_{T in S, {i,j} subset T} |c_T|^2 = (B_S w)_ij.

   If B_S has full column rank, the pair occupations recover the complete
   squared-amplitude vector w on that support.  The generated artifact contains
   one nonzero integer maximal minor of B_S for every shipped state.

2. On the exact fixed-spectrum v_B family, the first two nontrivial spectral
   moments of the 2-RDM are computed symbolically.  The quadratic moment is
   constant, while the cubic moment is strictly affine in the family parameter.
   The cubic moment is invariant under every one-body unitary and therefore
   separates distinct t-values in the reduced fiber.

The scope is deliberately narrow.  Full column rank is a certified property of
the 799 shipped supports, not a universal property of determinant supports.
The v_B calculation treats the exact one-parameter fixed-spectrum slice, not
a global integration of any second reduced-fiber direction.

Usage:
  python scripts/fiber_symmetries.py
  python scripts/fiber_symmetries.py --check
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
STATES = DATA / "states.jsonl"
TRANSPORT = DATA / "two_block_family.json"
OUT = DATA / "fiber_symmetries.json"

VB_DETS = [
    (0, 1, 2, 4),
    (0, 1, 2, 8),
    (0, 1, 3, 5),
    (0, 2, 3, 6),
    (0, 3, 4, 7),
    (0, 3, 7, 8),
    (1, 3, 6, 8),
    (2, 3, 4, 8),
]
VB_WEIGHTS = (1, 8, 4, 3, 2, 2, 1, 2)
VB_KERNEL = (1, -1, 0, 0, -1, 1, 0, 0)
VB_DENOMINATOR = 23


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _records() -> list[dict]:
    return [json.loads(line) for line in STATES.read_text().splitlines() if line.strip()]


def pair_rows(d: int) -> list[tuple[int, int]]:
    """Lexicographically ordered orbital pairs."""
    return list(itertools.combinations(range(d), 2))


def pair_incidence(
    support: list[tuple[int, ...]], d: int, rows: list[tuple[int, int]] | None = None
) -> list[list[int]]:
    """Pair-by-determinant incidence matrix B_S."""
    rows = pair_rows(d) if rows is None else rows
    return [[int(i in det and j in det) for det in support] for i, j in rows]


def _minor_certificate(record: dict) -> dict:
    """One exact nonzero maximal minor of the pair-incidence matrix."""
    import sympy as sp

    support = [tuple(det) for det in record["closed_form"]["support_dets"]]
    d = len(record["integer_form"])
    pairs = pair_rows(d)
    matrix = sp.Matrix(pair_incidence(support, d, pairs))
    _, pivots = matrix.T.rref()
    if len(pivots) != len(support):
        raise AssertionError(
            f"pair-incidence rank failure at {record['system']} v{record['index']}"
        )
    chosen = [pairs[i] for i in pivots]
    minor = sp.Matrix(pair_incidence(support, d, chosen))
    determinant = int(minor.det())
    if determinant == 0:
        raise AssertionError("RREF selected a zero determinant")
    return {
        "system": record["system"],
        "index": record["index"],
        "classified": record["classified"],
        "support_size": len(support),
        "minor_pairs": [list(pair) for pair in chosen],
        "minor_determinant": determinant,
    }


def _incidence_kernel_dim(record: dict) -> int:
    import sympy as sp

    support = [tuple(det) for det in record["closed_form"]["support_dets"]]
    d = len(record["integer_form"])
    matrix = sp.Matrix([[int(i in det) for det in support] for i in range(d)])
    return len(support) - matrix.rank()


def _permutation_sign(sequence: list[int]) -> int:
    inversions = sum(
        sequence[i] > sequence[j]
        for i in range(len(sequence))
        for j in range(i + 1, len(sequence))
    )
    return -1 if inversions % 2 else 1


@functools.lru_cache(maxsize=1)
def _vb_symbolic_moments():
    """Return tr(Gamma_2^k), k=1,2,3, after exact family reduction.

    The pair-vector convention is the same as gpc_census.gauge._pair_vectors.
    Before inserting the family, write u_i = |c_i|^2 and
    h = |c_0 c_1 c_4 c_5|.  Every phase-dependent monomial through the cubic
    trace is c*h, where c = cos(theta).  The exact coherence equation gives

        c*h = (3 + 7*t - 2*t^2) / (2*23^2).

    This avoids numerical eigenvalues and avoids choosing radical branches.
    """
    import sympy as sp

    c, s = sp.symbols("c s", real=True)
    radii = sp.symbols("r0:8", real=True)
    weights = sp.symbols("u0:8", real=True)
    h, ch, t = sp.symbols("h ch t", real=True)
    amplitudes = [
        radii[i] * (c + sp.I * s if i == 0 else 1) for i in range(len(VB_DETS))
    ]

    pairs = pair_rows(9)
    remainders = pairs
    remainder_pos = {remainder: i for i, remainder in enumerate(remainders)}
    determinant_pos = {det: i for i, det in enumerate(VB_DETS)}
    vectors = sp.zeros(len(pairs), len(remainders))
    for pair_index, (p, q) in enumerate(pairs):
        for remainder in remainders:
            if p in remainder or q in remainder:
                continue
            determinant = tuple(sorted((p, q) + remainder))
            amplitude_index = determinant_pos.get(determinant)
            if amplitude_index is None:
                continue
            sign = _permutation_sign([p, q, *remainder])
            vectors[pair_index, remainder_pos[remainder]] = (
                sign * amplitudes[amplitude_index]
            )
    gamma2 = sp.conjugate(vectors) * vectors.T

    def radial_reduce(expression):
        expression = sp.Poly(sp.expand(expression), s).rem(
            sp.Poly(s**2 - (1 - c**2), s)
        ).as_expr()
        out = 0
        for powers, coefficient in sp.Poly(sp.expand(expression), *radii, c).terms():
            radial_powers = powers[: len(radii)]
            odd = tuple(i for i, power in enumerate(radial_powers) if power % 2)
            if odd not in ((), (0, 1, 4, 5)):
                raise AssertionError(f"unexpected radical monomial {odd}")
            term = coefficient * c ** powers[-1]
            if odd:
                term *= h
            for i, power in enumerate(radial_powers):
                term *= weights[i] ** (power // 2)
            out += term
        out = sp.expand(out)
        collected = sp.collect(out, c * h, evaluate=False)
        if not set(collected).issubset({sp.Integer(1), c * h}):
            raise AssertionError("unexpected phase dependence in 2-RDM moment")
        return collected.get(sp.Integer(1), 0) + ch * collected.get(c * h, 0)

    substitutions = {
        weights[i]: (sp.Integer(VB_WEIGHTS[i]) + VB_KERNEL[i] * t) / VB_DENOMINATOR
        for i in range(len(weights))
    }
    coherence = (3 + 7 * t - 2 * t**2) / (2 * VB_DENOMINATOR**2)
    moments = []
    for power in (1, 2, 3):
        reduced = radial_reduce(sp.trace(gamma2**power))
        moments.append(sp.factor(reduced.subs(substitutions).subs(ch, coherence)))

    odd_entry = sp.factor(
        gamma2[pairs.index((2, 4)), pairs.index((3, 5))]
        .subs(radii[0], sp.sqrt((1 + t) / VB_DENOMINATOR))
        .subs(radii[2], sp.sqrt(sp.Rational(4, VB_DENOMINATOR)))
    )
    return t, moments, odd_entry


@functools.lru_cache(maxsize=1)
def _vb_exact_tangent_moment_rank() -> dict:
    """Exact first-order rank of the cubic and quartic moments at v_B, t=0.

    C is the real differential of the fixed-spectrum conditions
    P_lambda d(rho) P_lambda = 0, plus normalization, in the 16 real
    amplitude coordinates of the support.  Its exact rank is seven.  The
    orbital-phase gauge also has rank seven and lies in ker(C), leaving a
    two-dimensional tangent after that quotient. Appending d tr(Gamma_2^3)
    and d tr(Gamma_2^4) raises the rank from seven to nine, so these two
    invariant covectors are independent on that quotient tangent. The full
    degeneracy-stabilizer quotient is a separate question.
    """
    import sympy as sp

    imaginary = sp.I
    amplitudes = [
        sp.sqrt(sp.Rational(weight, VB_DENOMINATOR)) for weight in VB_WEIGHTS
    ]
    amplitudes[0] = (
        sp.sqrt(sp.Rational(1, VB_DENOMINATOR))
        * (3 * sp.sqrt(2) + imaginary * sp.sqrt(238))
        / 16
    )

    def sesquilinear(left, right):
        determinant_pos = {det: i for i, det in enumerate(VB_DETS)}
        matrix = sp.zeros(9)
        for determinant_index, determinant in enumerate(VB_DETS):
            for removed_index, removed in enumerate(determinant):
                remainder = [orbital for orbital in determinant if orbital != removed]
                for added in range(9):
                    if added in remainder:
                        continue
                    target = tuple(sorted(remainder + [added]))
                    target_index = determinant_pos.get(target)
                    if target_index is None:
                        continue
                    matrix[added, removed] += (
                        sp.conjugate(left[target_index])
                        * right[determinant_index]
                        * (-1) ** removed_index
                        * (-1) ** target.index(added)
                    )
        return matrix

    rho = sp.simplify(sesquilinear(amplitudes, amplitudes))
    directions = []
    rdm_differentials = []
    for coordinate in range(16):
        direction = [sp.Integer(0)] * 8
        direction[coordinate % 8] = 1 if coordinate < 8 else imaginary
        directions.append(direction)
        rdm_differentials.append(
            sp.simplify(
                sesquilinear(direction, amplitudes)
                + sesquilinear(amplitudes, direction)
            )
        )

    eigenvalues = [
        sp.Rational(20, 23),
        sp.Rational(14, 23),
        sp.Rational(4, 23),
    ]
    identity = sp.eye(9)
    projectors = []
    for eigenvalue in eigenvalues:
        projector = identity
        for other in eigenvalues:
            if other != eigenvalue:
                projector = projector * (rho - other * identity) / (eigenvalue - other)
        projector = sp.simplify(projector)
        if sp.simplify(projector * projector - projector) != sp.zeros(9):
            raise AssertionError("v_B exact spectral projector check failed")
        projectors.append(projector)

    constraint_rows = []
    for projector in projectors:
        blocks = [
            sp.simplify(projector * differential * projector)
            for differential in rdm_differentials
        ]
        for p in range(9):
            constraint_rows.append(
                [sp.simplify(sp.re(block[p, p])) for block in blocks]
            )
            for q in range(p + 1, 9):
                real_row = [sp.simplify(sp.re(block[p, q])) for block in blocks]
                imag_row = [sp.simplify(sp.im(block[p, q])) for block in blocks]
                if any(entry != 0 for entry in real_row):
                    constraint_rows.append(real_row)
                if any(entry != 0 for entry in imag_row):
                    constraint_rows.append(imag_row)
    constraint_rows.append(
        [
            sp.simplify(
                2
                * sp.re(
                    sum(
                        sp.conjugate(amplitudes[i]) * directions[coordinate][i]
                        for i in range(8)
                    )
                )
            )
            for coordinate in range(16)
        ]
    )
    constraints = sp.Matrix(constraint_rows)
    constraint_rank = constraints.rank()
    if constraints.shape != (32, 16) or constraint_rank != 7:
        raise AssertionError(
            f"unexpected v_B constraint rank {constraints.shape}, {constraint_rank}"
        )

    gauge_rows = []
    for orbital in range(9):
        direction = [
            imaginary * int(orbital in determinant) * amplitudes[i]
            for i, determinant in enumerate(VB_DETS)
        ]
        gauge_rows.append(
            [sp.re(value) for value in direction]
            + [sp.im(value) for value in direction]
        )
    gauge = sp.Matrix(gauge_rows)
    gauge_rank = gauge.rank()
    if gauge_rank != 7 or sp.simplify(constraints * gauge.T) != sp.zeros(32, 9):
        raise AssertionError("v_B exact gauge-kernel check failed")

    pairs = pair_rows(9)
    pair_pos = {pair: i for i, pair in enumerate(pairs)}
    determinant_pos = {det: i for i, det in enumerate(VB_DETS)}

    def pair_vectors(values):
        vectors = sp.zeros(len(pairs), len(pairs))
        for pair_index, (p, q) in enumerate(pairs):
            for remainder in pairs:
                if p in remainder or q in remainder:
                    continue
                determinant = tuple(sorted((p, q) + remainder))
                amplitude_index = determinant_pos.get(determinant)
                if amplitude_index is None:
                    continue
                vectors[pair_index, pair_pos[remainder]] = (
                    _permutation_sign([p, q, *remainder])
                    * values[amplitude_index]
                )
        return vectors

    vectors = pair_vectors(amplitudes)
    gamma2 = sp.conjugate(vectors) * vectors.T
    moment_rows = []
    for power in (3, 4):
        row = []
        for direction in directions:
            delta_vectors = pair_vectors(direction)
            delta_gamma2 = (
                sp.conjugate(delta_vectors) * vectors.T
                + sp.conjugate(vectors) * delta_vectors.T
            )
            row.append(
                sp.simplify(
                    sp.re(
                        power
                        * sp.trace(gamma2 ** (power - 1) * delta_gamma2)
                    )
                )
            )
        moment_rows.append(row)
    moments = sp.Matrix(moment_rows)
    if sp.simplify(moments * gauge.T) != sp.zeros(2, 9):
        raise AssertionError("v_B 2-RDM moments do not annihilate gauge")

    selected_constraint_rows = [31, 10, 11, 19, 28, 0, 27]
    selected_columns = [1, 2, 3, 5, 6, 7, 4, 8, 0]
    combined = constraints[selected_constraint_rows, :].col_join(moments)
    determinant = sp.factor(combined[:, selected_columns].det())
    expected = sp.Rational(299925504, 1035662780341225) * sp.sqrt(8211)
    if sp.simplify(determinant - expected) != 0:
        raise AssertionError(f"unexpected v_B tangent-moment minor {determinant}")

    return {
        "point": "t=0, cos(theta)=3*sqrt(2)/16, sin(theta)=sqrt(238)/16",
        "ambient_real_dimension": 16,
        "constraint_matrix_shape": list(constraints.shape),
        "fixed_spectrum_constraint_rank": constraint_rank,
        "fixed_spectrum_kernel_dimension": 16 - constraint_rank,
        "orbital_phase_gauge_rank": gauge_rank,
        "quotient_tangent_dimension": 16 - constraint_rank - gauge_rank,
        "quotient_group": (
            "the orbital-phase gauge only; the full degeneracy stabilizer is "
            "a separate reduced-fiber quotient"
        ),
        "invariant_covectors": [
            "d trace(Gamma2^3)",
            "d trace(Gamma2^4)",
        ],
        "constraint_rows_in_minor": selected_constraint_rows,
        "amplitude_columns_in_minor": selected_columns,
        "combined_rank": 9,
        "combined_minor_determinant": sp.sstr(determinant),
        "combined_minor_determinant_squared": sp.sstr(sp.factor(determinant**2)),
        "gauge_annihilated_by_constraints": True,
        "gauge_annihilated_by_invariant_covectors": True,
        "conclusion": (
            "the cubic and quartic 2-RDM moment differentials form a basis "
            "of the cotangent space of the two-dimensional support-restricted "
            "fixed-spectrum tangent modulo the orbital-phase gauge at t=0"
        ),
        "scope_limit": (
            "first-order exact statement only; it does not prove global "
            "injectivity or integrability of every tangent direction"
        ),
        "evidence": "exact algebraic spectral projectors and a nonzero 9x9 minor",
    }


def _vb_support_permutation_automorphisms() -> list[tuple[int, ...]]:
    """Block-preserving orbital permutations that preserve the v_B support."""
    support = set(VB_DETS)
    out = []
    for image_14 in itertools.permutations((1, 2, 3, 4)):
        first = dict(zip((1, 2, 3, 4), image_14))
        for image_58 in itertools.permutations((5, 6, 7, 8)):
            mapping = {0: 0, **first, **dict(zip((5, 6, 7, 8), image_58))}
            moved = {
                tuple(sorted(mapping[orbital] for orbital in determinant))
                for determinant in VB_DETS
            }
            if moved == support:
                out.append(tuple(mapping[i] for i in range(9)))
    return out


def _vb_certificate() -> dict:
    import sympy as sp

    t, moments, odd_entry = _vb_symbolic_moments()
    expected = [
        sp.Integer(6),
        sp.Rational(1438, 529),
        6 * (2899 - 8 * t) / 12167,
    ]
    if [sp.factor(got - want) for got, want in zip(moments, expected)] != [0, 0, 0]:
        raise AssertionError("v_B 2-RDM moment identity failed")

    t_minus = sp.Rational(7, 17) - 18 * sp.sqrt(35) / 85
    t_plus = sp.Rational(7, 17) + 18 * sp.sqrt(35) / 85
    wall_minus = sp.factor(moments[2].subs(t, t_minus))
    wall_plus = sp.factor(moments[2].subs(t, t_plus))
    wall_difference = sp.factor(wall_minus - wall_plus)
    automorphisms = _vb_support_permutation_automorphisms()
    if automorphisms != [tuple(range(9))]:
        raise AssertionError("unexpected block-permutation symmetry on the v_B support")

    return {
        "system": "(4,9)",
        "index": 65,
        "scope": (
            "the exact one-parameter fixed-spectrum slice on the shipped "
            "eight-determinant support; not the full two-dimensional reduced "
            "fixed-spectrum fiber"
        ),
        "family_weights": [
            "(1+t)/23",
            "(8-t)/23",
            "4/23",
            "3/23",
            "(2-t)/23",
            "(2+t)/23",
            "1/23",
            "2/23",
        ],
        "physical_walls": [sp.sstr(t_minus), sp.sstr(t_plus)],
        "time_reversal": {
            "action": "theta -> -theta",
            "fixed_points_modulo_orbital_phase_gauge": "the two real walls, sin(theta)=0",
            "pair_coordinate_zero_based_modes_1_4": "(1+t)/23",
            "two_rdm_odd_entry": (
                "Gamma2[(2,4),(3,5)] = "
                "2*sqrt(1+t)*(cos(theta)-I*sin(theta))/23"
            ),
            "two_rdm_odd_entry_symbolic_check": sp.sstr(odd_entry),
            "full_stabilizer_status": (
                "whether the two conjugate interior sheets at fixed t are "
                "equivalent under the full U(1)xU(4)xU(4) degeneracy "
                "stabilizer remains open"
            ),
        },
        "two_rdm_normalization": "trace(Gamma2) = binomial(4,2) = 6",
        "two_rdm_moments": {
            "trace_1": sp.sstr(moments[0]),
            "trace_2": sp.sstr(moments[1]),
            "trace_3": sp.sstr(moments[2]),
            "trace_3_derivative": sp.sstr(sp.diff(moments[2], t)),
            "recover_t": "t = (2899 - (12167/6)*trace(Gamma2^3))/8",
        },
        "galois_walls": {
            "trace_3_at_t_minus": sp.sstr(wall_minus),
            "trace_3_at_t_plus": sp.sstr(wall_plus),
            "difference": sp.sstr(wall_difference),
            "consequence": (
                "the two Galois-conjugate real walls are not related by a "
                "one-body unitary"
            ),
        },
        "support_preserving_block_permutation_automorphisms": len(automorphisms),
        "support_permutation_group": "trivial",
        "tangent_moment_rank": _vb_exact_tangent_moment_rank(),
        "evidence": "exact symbolic 2-RDM Gram-matrix expansion",
    }


def build() -> dict:
    records = _records()
    transport_rows = json.loads(TRANSPORT.read_text())["states"]
    transport_by_key = {
        (record["system"], record["index"]): record["transport_class"]
        for record in transport_rows
    }
    if len(transport_by_key) != len(records):
        raise AssertionError("transport metadata does not cover the state atlas")

    certificates = []
    kernel_dimensions = {}
    for record in records:
        key = (record["system"], record["index"])
        certificates.append(_minor_certificate(record))
        kernel_dimensions[key] = _incidence_kernel_dim(record)

    positive = [key for key, dimension in kernel_dimensions.items() if dimension > 0]
    positive_classes = {transport_by_key[key] for key in positive}
    all_classes = set(transport_by_key.values())
    class_counts = Counter(record["classified"] for record in records)
    kernel_counts = Counter(kernel_dimensions.values())

    return {
        "generated_by": "scripts/fiber_symmetries.py",
        "states_sha256": _sha256(STATES),
        "transport_metadata_sha256": _sha256(TRANSPORT),
        "conditional_theorem": (
            "For a fixed support S, full column rank of the pair-incidence "
            "matrix B_S makes the diagonal pair occupations an injective "
            "coordinate on squared amplitudes. Since sum_{j != i} "
            "(B_S w)_{ij} = (N-1)(A_S w)_i, ker(B_S) is contained in "
            "ker(A_S), so every nonzero incidence-kernel deformation changes "
            "a pair occupation."
        ),
        "scope_limit": (
            "The rank hypothesis is not universal. Signed combinatorial "
            "2-trades give supports with a nonzero pair-incidence kernel. "
            "The artifact certifies the hypothesis only for the shipped "
            "supports. Pair occupations use a fixed orbital frame and are "
            "not invariants of the full degeneracy-stabilizer quotient."
        ),
        "functoriality": (
            "Full pair-column rank is invariant under orbital relabeling, "
            "empty-orbital padding, frozen-core lift, and Hodge complement. "
            "For Hodge complement this assumes both the particle and hole "
            "numbers are at least two; the identity is "
            "B_complement = 1 - A_i - A_j + B."
        ),
        "summary": {
            "records": len(records),
            "pair_minor_certificates": len(certificates),
            "full_column_rank": len(certificates),
            "rank_failures": 0,
            "transport_classes": len(all_classes),
            "classification_counts": dict(sorted(class_counts.items())),
            "incidence_kernel_dimension_counts": {
                str(key): value for key, value in sorted(kernel_counts.items())
            },
            "positive_incidence_kernel_records": len(positive),
            "positive_incidence_kernel_transport_classes": len(positive_classes),
            "positive_kernel_pair_map_failures": 0,
        },
        "vb": _vb_certificate(),
        "records": certificates,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    text = json.dumps(build(), indent=1) + "\n"
    if args.check:
        current = args.out.read_text() if args.out.exists() else ""
        if current != text:
            print(f"STALE: {args.out}", file=sys.stderr)
            return 1
        print(f"{args.out} is current")
        return 0

    args.out.write_text(text)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
