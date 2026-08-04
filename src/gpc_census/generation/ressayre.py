"""Exact Ressayre enumeration for multiplicity-free exterior representations.

The convention is that a certificate ``(H,z)`` proves

    H . lambda >= z.

Negative type-A roots are ``e_j-e_i`` with ``i < j``.  The determinant
matrix is the tangent map from root spaces with ``H_j-H_i < 0`` to weights
strictly below the affine hyperplane.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Sequence

import sympy as sp

from .exterior import exterior_root_action, exterior_weights, negative_roots
from .model import (
    AdmissibleHyperplane,
    RessayreCertificate,
    RessayreEvaluationVerdict,
    RessayreEnumeration,
    Root,
    Weight,
)


def _bareiss_det(matrix: Sequence[Sequence[int]]) -> int:
    """Return an exact integer determinant using fraction-free elimination."""
    data = [list(row) for row in matrix]
    size = len(data)
    if any(len(row) != size for row in data):
        raise ValueError("determinant requires a square matrix")
    if any(type(value) is not int for row in data for value in row):
        raise TypeError("determinant entries must be exact integers")
    if size == 0:
        return 1
    if size == 1:
        return data[0][0]

    sign = 1
    previous = 1
    for column in range(size - 1):
        if data[column][column] == 0:
            pivot_row = next(
                (row for row in range(column + 1, size) if data[row][column]), None
            )
            if pivot_row is None:
                return 0
            data[column], data[pivot_row] = data[pivot_row], data[column]
            sign = -sign
        pivot = data[column][column]
        for row in range(column + 1, size):
            for target_column in range(column + 1, size):
                numerator = (
                    data[row][target_column] * pivot
                    - data[row][column] * data[column][target_column]
                )
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("Bareiss division was not exact")
                data[row][target_column] = quotient
        previous = pivot
    return sign * data[-1][-1]


def _cofactor_nullvector(matrix: Sequence[Sequence[int]]) -> tuple[int, ...]:
    """Return the signed-maximal-minor nullvector of an ``(d-1) x d`` matrix."""
    rows = [list(row) for row in matrix]
    if not rows or len(rows[0]) != len(rows) + 1:
        raise ValueError("cofactor nullvector requires an (d-1) x d matrix")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("matrix rows have inconsistent widths")
    return tuple(
        (-1) ** column
        * _bareiss_det([row[:column] + row[column + 1 :] for row in rows])
        for column in range(width)
    )


def _primitive_pair(h: Sequence[int], z: int, *, orient: bool) -> tuple[tuple[int, ...], int]:
    values = list(h) + [z]
    divisor = 0
    for value in values:
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:
        raise ValueError("zero hyperplane")
    values = [value // divisor for value in values]
    if orient and next(value for value in values if value) < 0:
        values = [-value for value in values]
    return tuple(values[:-1]), values[-1]


def _modular_row_witness(
    rows: Sequence[Sequence[int]], target_rank: int
) -> tuple[int, ...] | None:
    """Find rows independent over a finite field, hence also over ``Q``.

    The exterior-weight incidence rows are small, and a nonzero minor modulo
    any prime is already an exact characteristic-zero independence proof.
    Trying several small primes avoids most symbolic rank computations while
    retaining a rational fallback for minors divisible by every tried prime.
    """
    if target_rank < 0:
        raise ValueError("target_rank must be nonnegative")
    if target_rank == 0:
        return ()
    integral_rows = tuple(tuple(row) for row in rows)
    if not integral_rows:
        return None
    width = len(integral_rows[0])
    if any(len(row) != width for row in integral_rows):
        raise ValueError("rank rows have inconsistent widths")
    if any(type(value) is not int for row in integral_rows for value in row):
        raise TypeError("rank rows must contain exact integers")

    for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        basis: dict[int, tuple[int, ...]] = {}
        selected: list[int] = []
        for row_index, row in enumerate(integral_rows):
            reduced = [value % prime for value in row]
            for pivot in sorted(basis):
                factor = reduced[pivot]
                if factor:
                    basis_row = basis[pivot]
                    reduced = [
                        (value - factor * basis_value) % prime
                        for value, basis_value in zip(reduced, basis_row, strict=True)
                    ]
            pivot = next((index for index, value in enumerate(reduced) if value), None)
            if pivot is None:
                continue
            inverse = pow(reduced[pivot], prime - 2, prime)
            normalized = tuple(value * inverse % prime for value in reduced)
            basis[pivot] = normalized
            selected.append(row_index)
            if len(selected) == target_rank:
                return tuple(selected)
    return None


def _rational_row_witness(
    rows: Sequence[Sequence[int]], target_rank: int
) -> tuple[int, ...] | None:
    """Find ``target_rank`` independent rows with an exact rational fallback."""
    integral_rows = tuple(tuple(row) for row in rows)
    modular = _modular_row_witness(integral_rows, target_rank)
    if modular is not None:
        return modular
    witness: list[int] = []
    rank = 0
    for index in range(len(integral_rows)):
        trial = witness + [index]
        trial_rank = sp.Matrix([integral_rows[row] for row in trial]).rank()
        if trial_rank > rank:
            witness.append(index)
            rank = trial_rank
        if rank == target_rank:
            return tuple(witness)
    return None


def enumerate_admissible_hyperplanes(n: int, d: int) -> tuple[AdmissibleHyperplane, ...]:
    """Exhaust all affine weight hyperplanes for ``wedge^n C^d``.

    Exterior weights span the trace hyperplane of dimension ``d-1``.  Every
    admissible hyperplane therefore contains ``d-1`` affinely independent
    weights.  Enumerating all such subsets is complete, although deliberately
    exponential.  Higher-rank callers will need a proven pruning backend.
    """
    weights = exterior_weights(n, d)
    hyperplanes: dict[tuple[tuple[int, ...], int], tuple[int, ...]] = {}
    for witness in itertools.combinations(range(len(weights)), d - 1):
        base = weights[witness[0]]
        difference_rows = [
            [weights[index][column] - base[column] for column in range(d)]
            for index in witness[1:]
        ]
        # Fix the U(1) ambiguity by requiring sum(H)=0.
        normal = _cofactor_nullvector(difference_rows + [[1] * d])
        if not any(normal):
            continue
        offset = sum(normal[index] * base[index] for index in range(d))
        h, z = _primitive_pair(normal, offset, orient=True)
        hyperplanes.setdefault((h, z), witness)
    return tuple(
        AdmissibleHyperplane(h=h, z=z, witness=hyperplanes[(h, z)])
        for h, z in sorted(hyperplanes)
    )


def _level_sets(
    weights: Sequence[Weight], h: Sequence[int], z: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    on: list[int] = []
    below: list[int] = []
    for index, weight in enumerate(weights):
        value = sum(left * right for left, right in zip(weight, h, strict=True))
        if value == z:
            on.append(index)
        elif value < z:
            below.append(index)
    return tuple(on), tuple(below)


def _negative_root_set(h: Sequence[int]) -> tuple[Root, ...]:
    return tuple(root for root in negative_roots(len(h)) if h[root[1]] - h[root[0]] < 0)


def _symbolic_determinant(
    weights: Sequence[Weight],
    on_indices: Sequence[int],
    below_indices: Sequence[int],
    roots: Sequence[Root],
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    variables = sp.symbols(f"x0:{len(on_indices)}")
    below_rows = {weights[index]: row for row, index in enumerate(below_indices)}
    matrix = sp.zeros(len(below_indices), len(roots))
    for column, root in enumerate(roots):
        for variable_index, weight_index in enumerate(on_indices):
            action = exterior_root_action(weights[weight_index], root)
            if action is None or action[0] not in below_rows:
                continue
            row = below_rows[action[0]]
            matrix[row, column] += action[1] * variables[variable_index]
    return sp.expand(matrix.det(method="domain-ge")), tuple(variables)


def _tangent_edge_matrix(
    weights: Sequence[Weight],
    on_indices: Sequence[int],
    below_indices: Sequence[int],
    roots: Sequence[Root],
) -> tuple[tuple[tuple[int, int] | None, ...], ...]:
    """Return tangent entries as ``(variable_index, sign)`` edge labels."""
    below_rows = {weights[index]: row for row, index in enumerate(below_indices)}
    matrix: list[list[tuple[int, int] | None]] = [
        [None for _root in roots] for _weight in below_indices
    ]
    for column, root in enumerate(roots):
        for variable_index, weight_index in enumerate(on_indices):
            action = exterior_root_action(weights[weight_index], root)
            if action is None or action[0] not in below_rows:
                continue
            row = below_rows[action[0]]
            if matrix[row][column] is not None:
                raise AssertionError("an exterior tangent cell has two preimages")
            matrix[row][column] = (variable_index, action[1])
    return tuple(tuple(row) for row in matrix)


def _sparse_determinant_coefficients(
    edge_matrix: Sequence[Sequence[tuple[int, int] | None]],
) -> dict[tuple[int, ...], int]:
    """Expand a sparse labeled determinant in exact integer arithmetic.

    A monomial is encoded by its sorted tuple of variable indices. Dynamic
    programming over used column masks combines equal monomials immediately,
    so exact cancellations are detected without constructing a large SymPy
    expression. An empty result is a replayable proof that the determinant
    polynomial is identically zero.
    """
    size = len(edge_matrix)
    if any(len(row) != size for row in edge_matrix):
        raise ValueError("sparse determinant requires a square matrix")
    states: dict[int, dict[tuple[int, ...], int]] = {0: {(): 1}}
    for row_index, row in enumerate(edge_matrix):
        next_states: dict[int, dict[tuple[int, ...], int]] = {}
        for used_columns, polynomial in states.items():
            if used_columns.bit_count() != row_index:
                raise AssertionError("determinant DP column mask has the wrong degree")
            for column, edge in enumerate(row):
                if edge is None or used_columns & (1 << column):
                    continue
                variable_index, edge_sign = edge
                inversion_sign = (
                    -1
                    if (used_columns >> (column + 1)).bit_count() % 2
                    else 1
                )
                target = next_states.setdefault(used_columns | (1 << column), {})
                factor = edge_sign * inversion_sign
                for monomial, coefficient in polynomial.items():
                    extended = tuple(sorted((*monomial, variable_index)))
                    updated = target.get(extended, 0) + factor * coefficient
                    if updated:
                        target[extended] = updated
                    else:
                        target.pop(extended, None)
        states = {
            mask: polynomial
            for mask, polynomial in next_states.items()
            if polynomial
        }
        if not states:
            return {}
    return states.get((1 << size) - 1, {})


def tangent_determinant_is_identically_zero(
    n: int,
    d: int,
    candidate: AdmissibleHyperplane,
) -> bool:
    """Decide exact vanishing of a trace-matched tangent determinant.

    This routine is intended as the fail-closed fallback after evaluated
    determinant discovery. It expands the sparse determinant by labeled
    matchings and therefore returns a mathematical zero decision, not a
    probabilistic generic-rank result.
    """
    weights, on_indices, below_indices, roots = _validate_admissible_candidate(
        n, d, candidate
    )
    if len(below_indices) != len(roots):
        raise ValueError("exact determinant fallback requires trace equality")
    edge_matrix = _tangent_edge_matrix(
        weights,
        on_indices,
        below_indices,
        roots,
    )
    return not _sparse_determinant_coefficients(edge_matrix)


def _nonzero_evaluation(
    polynomial: sp.Expr, variables: Sequence[sp.Symbol], degree_bound: int
) -> tuple[tuple[int, ...], int]:
    """Find a guaranteed exact evaluation of a known nonzero polynomial.

    At each variable, a nonzero polynomial of degree at most ``degree_bound``
    cannot vanish at every integer in ``0..degree_bound``.  Symbolic partial
    substitution therefore constructs a witness deterministically.
    """
    current = sp.expand(polynomial)
    evaluation: list[int] = []
    for variable in variables:
        selected: tuple[int, sp.Expr] | None = None
        for value in range(degree_bound + 1):
            candidate = sp.expand(current.subs(variable, value))
            if candidate != 0:
                selected = value, candidate
                break
        if selected is None:
            raise AssertionError("nonzero polynomial vanished on a sufficient grid")
        value, current = selected
        evaluation.append(value)
    if getattr(current, "free_symbols", set()):
        raise AssertionError("determinant evaluation left free variables")
    result = int(current)
    if result == 0:
        raise AssertionError("determinant witness is zero")
    return tuple(evaluation), result


def _evaluated_tangent_matrix(
    weights: Sequence[Weight],
    on_indices: Sequence[int],
    below_indices: Sequence[int],
    roots: Sequence[Root],
    evaluation: Sequence[int],
) -> list[list[int]]:
    below_rows = {weights[index]: row for row, index in enumerate(below_indices)}
    matrix = [[0 for _root in roots] for _weight in below_indices]
    for column, root in enumerate(roots):
        for variable_index, weight_index in enumerate(on_indices):
            action = exterior_root_action(weights[weight_index], root)
            if action is None or action[0] not in below_rows:
                continue
            matrix[below_rows[action[0]]][column] += action[1] * evaluation[variable_index]
    return matrix


def certify_candidate(
    n: int,
    d: int,
    candidate: AdmissibleHyperplane,
    *,
    orientation: int = 1,
) -> RessayreCertificate | None:
    """Return a certificate exactly when an oriented candidate is Ressayre."""
    if orientation not in (-1, 1):
        raise ValueError("orientation must be -1 or 1")
    weights = exterior_weights(n, d)
    h = tuple(orientation * value for value in candidate.h)
    z = orientation * candidate.z
    on_indices, below_indices = _level_sets(weights, h, z)
    roots = _negative_root_set(h)
    if len(below_indices) != len(roots):
        return None

    if not roots:
        evaluation: tuple[int, ...] = tuple(0 for _index in on_indices)
        determinant_value = 1
    else:
        determinant, variables = _symbolic_determinant(
            weights, on_indices, below_indices, roots
        )
        if determinant == 0:
            return None
        evaluation, determinant_value = _nonzero_evaluation(
            determinant, variables, len(roots)
        )
        replay = _bareiss_det(
            _evaluated_tangent_matrix(
                weights, on_indices, below_indices, roots, evaluation
            )
        )
        if replay != determinant_value:
            raise AssertionError("symbolic and replayed determinant values disagree")

    return RessayreCertificate(
        h=h,
        z=z,
        admissibility_witness=candidate.witness,
        on_indices=on_indices,
        below_indices=below_indices,
        negative_roots=roots,
        determinant_evaluation=evaluation,
        determinant_value=determinant_value,
    )


def admissible_candidate_from_tau(n: int, tau: Sequence[int]) -> AdmissibleHyperplane:
    """Convert ``tau.lambda <= 0`` into the centered Ressayre convention.

    On ``sum(lambda)=n``, set ``S=sum(tau)``, ``H_i=S-d*tau_i`` and
    ``z=n*S``.  Then ``H.lambda >= z`` is exactly the original inequality.
    The returned witness is selected greedily from the weights on the
    hyperplane and proves affine codimension one over the rationals.
    """
    if not tau or not any(tau):
        raise ValueError("tau must be a nonzero vector")
    if any(type(value) is not int for value in tau):
        raise TypeError("tau entries must be integers")
    d = len(tau)
    total = sum(tau)
    h, z = _primitive_pair(
        tuple(total - d * value for value in tau),
        n * total,
        orient=False,
    )
    weights = exterior_weights(n, d)
    on_indices, _below = _level_sets(weights, h, z)
    on_rows = tuple(tuple(weights[index]) + (1,) for index in on_indices)
    local_witness = _rational_row_witness(on_rows, d - 1)
    if local_witness is None:
        raise ValueError("tau does not define an admissible exterior-weight hyperplane")
    witness = tuple(on_indices[index] for index in local_witness)
    return AdmissibleHyperplane(h=h, z=z, witness=witness)


def _deterministic_evaluation(
    h: Sequence[int], z: int, size: int, attempt: int
) -> tuple[int, ...]:
    """Return a stable small-integer point for determinant evaluation."""
    modulus = 2**31
    state = attempt + 1
    for value in tuple(h) + (z,):
        state = (1_103_515_245 * state + value + 12_345) % modulus
    result: list[int] = []
    for _ in range(size):
        state = (1_103_515_245 * state + 12_345) % modulus
        result.append(1 + state % 97)
    return tuple(result)


def _validate_admissible_candidate(
    n: int,
    d: int,
    candidate: AdmissibleHyperplane,
) -> tuple[tuple[Weight, ...], tuple[int, ...], tuple[int, ...], tuple[Root, ...]]:
    """Validate the exact preconditions shared by evaluated certification.

    A determinant witness has mathematical force only for a primitive centered
    admissible hyperplane.  Checking those conditions before evaluation keeps
    the public certifier from returning superficially nonzero but unreplayable
    payloads for fabricated candidates.
    """
    if n <= 0 or d < n:
        raise ValueError("require 0 < n <= d")
    if len(candidate.h) != d or any(type(value) is not int for value in candidate.h):
        raise ValueError("candidate normal has the wrong dimension or type")
    if type(candidate.z) is not int or sum(candidate.h) != 0:
        raise ValueError("candidate must be an integral centered hyperplane")
    primitive_h, primitive_z = _primitive_pair(candidate.h, candidate.z, orient=False)
    if primitive_h != candidate.h or primitive_z != candidate.z:
        raise ValueError("candidate hyperplane must be primitive")

    weights = exterior_weights(n, d)
    on_indices, below_indices = _level_sets(weights, candidate.h, candidate.z)
    roots = _negative_root_set(candidate.h)
    if len(candidate.witness) != d - 1:
        raise ValueError("candidate needs d-1 affinely independent witness weights")
    if any(type(index) is not int for index in candidate.witness):
        raise ValueError("candidate witness indices must be integers")
    if any(index < 0 or index >= len(weights) for index in candidate.witness):
        raise ValueError("candidate witness index is out of range")
    if len(set(candidate.witness)) != d - 1:
        raise ValueError("candidate witness indices must be distinct")
    if not set(candidate.witness) <= set(on_indices):
        raise ValueError("candidate witness is not contained in the hyperplane")
    witness_rows = tuple(
        tuple(weights[index]) + (1,) for index in candidate.witness
    )
    if _rational_row_witness(witness_rows, d - 1) is None:
        raise ValueError("candidate witness is not affinely independent")
    return weights, on_indices, below_indices, roots


def assess_candidate_by_evaluation(
    n: int,
    d: int,
    candidate: AdmissibleHyperplane,
    *,
    attempts: int = 32,
) -> RessayreEvaluationVerdict:
    """Return a tri-state exact trace and determinant-evaluation verdict.

    A returned nonzero integer determinant is a proof that the Ressayre
    determinant polynomial is nonzero.  Exhausting the deterministic trial
    points is explicitly ``evaluation_unresolved`` and is never treated as
    proof that the polynomial vanishes.  A trace-count mismatch is the exact
    ``trace_rejected`` status.
    """
    if attempts <= 0:
        raise ValueError("attempts must be positive")
    weights, on_indices, below_indices, roots = _validate_admissible_candidate(
        n, d, candidate
    )
    if len(below_indices) != len(roots):
        return RessayreEvaluationVerdict(
            status="trace_rejected",
            attempts=0,
            below_weight_count=len(below_indices),
            negative_root_count=len(roots),
            certificate=None,
        )
    if not roots:
        evaluations = (tuple(0 for _index in on_indices),)
    else:
        evaluations = (
            _deterministic_evaluation(candidate.h, candidate.z, len(on_indices), attempt)
            for attempt in range(attempts)
        )
    for attempt_number, evaluation in enumerate(evaluations, start=1):
        determinant_value = _bareiss_det(
            _evaluated_tangent_matrix(
                weights,
                on_indices,
                below_indices,
                roots,
                evaluation,
            )
        )
        if determinant_value:
            certificate = RessayreCertificate(
                h=candidate.h,
                z=candidate.z,
                admissibility_witness=candidate.witness,
                on_indices=on_indices,
                below_indices=below_indices,
                negative_roots=roots,
                determinant_evaluation=evaluation,
                determinant_value=determinant_value,
            )
            return RessayreEvaluationVerdict(
                status="certified",
                attempts=attempt_number,
                below_weight_count=len(below_indices),
                negative_root_count=len(roots),
                certificate=certificate,
            )
    return RessayreEvaluationVerdict(
        status="evaluation_unresolved",
        attempts=attempts,
        below_weight_count=len(below_indices),
        negative_root_count=len(roots),
        certificate=None,
    )


def certify_candidate_by_evaluation(
    n: int,
    d: int,
    candidate: AdmissibleHyperplane,
    *,
    attempts: int = 32,
) -> RessayreCertificate | None:
    """Return an exact witness, or ``None`` for either noncertified status.

    Call :func:`assess_candidate_by_evaluation` when exact trace rejection
    must be distinguished from inconclusive finite evaluation.
    """
    return assess_candidate_by_evaluation(
        n,
        d,
        candidate,
        attempts=attempts,
    ).certificate


def certify_tau_by_evaluation(
    n: int, tau: Sequence[int], *, attempts: int = 32
) -> RessayreCertificate | None:
    """Certify the oriented homogeneous inequality ``tau.lambda <= 0``."""
    candidate = admissible_candidate_from_tau(n, tau)
    return certify_candidate_by_evaluation(
        n,
        len(tau),
        candidate,
        attempts=attempts,
    )


def assess_tau_by_evaluation(
    n: int, tau: Sequence[int], *, attempts: int = 32
) -> RessayreEvaluationVerdict:
    """Return the tri-state verdict for ``tau.lambda <= 0``."""
    candidate = admissible_candidate_from_tau(n, tau)
    return assess_candidate_by_evaluation(
        n,
        len(tau),
        candidate,
        attempts=attempts,
    )


def verify_tau_ressayre_certificate(
    n: int,
    tau: Sequence[int],
    certificate: RessayreCertificate,
) -> bool:
    """Replay a certificate and bind it to one oriented homogeneous row."""
    try:
        candidate = admissible_candidate_from_tau(n, tau)
    except (TypeError, ValueError):
        return False
    return (
        certificate.h == candidate.h
        and certificate.z == candidate.z
        and verify_ressayre_certificate(n, len(tau), certificate)
    )


def enumerate_ressayre(n: int, d: int) -> RessayreEnumeration:
    """Exhaust and certify all Ressayre elements from admissible hyperplanes."""
    candidates = enumerate_admissible_hyperplanes(n, d)
    trace_count = 0
    certificates: list[RessayreCertificate] = []
    weights = exterior_weights(n, d)
    for candidate in candidates:
        for orientation in (1, -1):
            h = tuple(orientation * value for value in candidate.h)
            z = orientation * candidate.z
            _on, below = _level_sets(weights, h, z)
            roots = _negative_root_set(h)
            if len(below) != len(roots):
                continue
            trace_count += 1
            certificate = certify_candidate(n, d, candidate, orientation=orientation)
            if certificate is not None:
                certificates.append(certificate)
    return RessayreEnumeration(
        admissible_hyperplane_count=len(candidates),
        oriented_candidate_count=2 * len(candidates),
        trace_candidate_count=trace_count,
        certificates=tuple(sorted(certificates)),
    )


def verify_ressayre_certificate(n: int, d: int, certificate: RessayreCertificate) -> bool:
    """Replay every exact condition in a Ressayre certificate."""
    try:
        if type(n) is not int or type(d) is not int or n <= 0 or d < n:
            return False
        integer_fields = (
            certificate.h
            + (certificate.z, certificate.determinant_value)
            + certificate.admissibility_witness
            + certificate.on_indices
            + certificate.below_indices
            + tuple(value for root in certificate.negative_roots for value in root)
            + certificate.determinant_evaluation
        )
        if any(type(value) is not int for value in integer_fields):
            return False
        weights = exterior_weights(n, d)
        if len(certificate.h) != d or sum(certificate.h) != 0:
            return False
        primitive_h, primitive_z = _primitive_pair(certificate.h, certificate.z, orient=False)
        if primitive_h != certificate.h or primitive_z != certificate.z:
            return False

        on_indices, below_indices = _level_sets(weights, certificate.h, certificate.z)
        roots = _negative_root_set(certificate.h)
        if on_indices != certificate.on_indices:
            return False
        if below_indices != certificate.below_indices:
            return False
        if roots != certificate.negative_roots:
            return False
        if len(below_indices) != len(roots):
            return False

        if len(certificate.admissibility_witness) != d - 1:
            return False
        if not set(certificate.admissibility_witness) <= set(on_indices):
            return False
        witness_matrix = sp.Matrix(
            [list(weights[index]) + [1] for index in certificate.admissibility_witness]
        )
        if witness_matrix.rank() != d - 1:
            return False
        on_matrix = sp.Matrix([list(weights[index]) + [1] for index in on_indices])
        if on_matrix.rank() != d - 1:
            return False

        if len(certificate.determinant_evaluation) != len(on_indices):
            return False
        determinant = _bareiss_det(
            _evaluated_tangent_matrix(
                weights,
                on_indices,
                below_indices,
                roots,
                certificate.determinant_evaluation,
            )
        )
        return determinant == certificate.determinant_value and determinant != 0
    except (IndexError, TypeError, ValueError):
        return False
