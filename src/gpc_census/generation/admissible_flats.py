"""Exact affine-flat enumeration for exterior-weight hyperplanes.

The brute-force reference enumerator visits every affine basis of every
hyperplane. This module instead builds the lattice of closed affine flats one
rank at a time. It is still an exhaustive reference algorithm, but it
deduplicates a flat as soon as it is discovered.

All closure calculations are performed over a finite field whose prime is
larger than the Hadamard bound for every relevant minor. Each augmented
``wedge^n`` weight has exactly ``n + 1`` unit entries, so every ``k`` by ``k``
minor has absolute value at most ``(n + 1) ** (k / 2)``. A prime larger than
the rank-``d`` bound therefore preserves every rational rank and closure
relation exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .exterior import exterior_weights
from .model import AdmissibleHyperplane
from .ressayre import _cofactor_nullvector, _primitive_pair


_KNOWN_PRIMES = (
    257,
    65_537,
    2_147_483_647,
    2_305_843_009_213_693_951,
)


@dataclass(frozen=True)
class FlatHyperplaneEnumeration:
    """An exhaustive affine-flat enumeration and its rank counts."""

    n: int
    d: int
    modulus: int
    hadamard_bound_squared: int
    flat_counts_by_rank: tuple[int, ...]
    hyperplanes: tuple[AdmissibleHyperplane, ...]

    def as_dict(self) -> dict[str, object]:
        """Return a compact deterministic proof manifest."""
        return {
            "particle_number": self.n,
            "rank": self.d,
            "rank_preserving_modulus": self.modulus,
            "hadamard_bound_squared": self.hadamard_bound_squared,
            "flat_counts_by_affine_matrix_rank": list(self.flat_counts_by_rank),
            "admissible_hyperplane_count": len(self.hyperplanes),
        }


@dataclass(frozen=True)
class _FlatState:
    mask: int
    basis_indices: tuple[int, ...]
    rref_rows: tuple[tuple[int, ...], ...]
    pivot_columns: tuple[int, ...]


def rank_preserving_modulus(n: int, d: int) -> tuple[int, int]:
    """Return a prime that preserves all augmented-weight ranks over ``Q``.

    The second return value is the squared Hadamard bound ``(n + 1) ** d``.
    Using the squared comparison avoids floating-point arithmetic.
    """
    if type(n) is not int or type(d) is not int or not (0 < n < d):
        raise ValueError("require integer parameters with 0 < n < d")
    bound_squared = (n + 1) ** d
    modulus = next(
        (prime for prime in _KNOWN_PRIMES if prime * prime > bound_squared),
        None,
    )
    if modulus is None:
        raise ValueError(
            "no configured prime exceeds the augmented-weight Hadamard bound"
        )
    return modulus, bound_squared


def _rref_mod(
    rows: tuple[tuple[int, ...], ...], modulus: int
) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]:
    """Return reduced row-echelon form over the configured prime field."""
    if not rows:
        return (), ()
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("RREF rows have inconsistent widths")
    matrix = [[value % modulus for value in row] for row in rows]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(width):
        selected = next(
            (
                row
                for row in range(pivot_row, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], modulus - 2, modulus)
        matrix[pivot_row] = [
            value * inverse % modulus for value in matrix[pivot_row]
        ]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (value - factor * basis_value) % modulus
                for value, basis_value in zip(
                    matrix[row], matrix[pivot_row], strict=True
                )
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return (
        tuple(tuple(row) for row in matrix[:pivot_row]),
        tuple(pivot_columns),
    )


def _residual_mod(
    row: tuple[int, ...],
    rref_rows: tuple[tuple[int, ...], ...],
    pivot_columns: tuple[int, ...],
    modulus: int,
) -> tuple[int, ...]:
    """Reduce one row modulo the span of an RREF basis."""
    residual = [value % modulus for value in row]
    for basis_row, pivot in zip(rref_rows, pivot_columns, strict=True):
        factor = residual[pivot]
        if factor:
            residual = [
                (value - factor * basis_value) % modulus
                for value, basis_value in zip(residual, basis_row, strict=True)
            ]
    return tuple(residual)


def _projective_key(residual: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    """Normalize a nonzero residual to its finite-field projective class."""
    pivot = next((index for index, value in enumerate(residual) if value), None)
    if pivot is None:
        raise AssertionError("an element outside a closed flat had zero residual")
    inverse = pow(residual[pivot], modulus - 2, modulus)
    return tuple(value * inverse % modulus for value in residual)


def _state_from_basis(
    mask: int,
    basis_indices: tuple[int, ...],
    points: tuple[tuple[int, ...], ...],
    modulus: int,
) -> _FlatState:
    rows = tuple(points[index] for index in basis_indices)
    rref_rows, pivots = _rref_mod(rows, modulus)
    if len(rref_rows) != len(basis_indices):
        raise AssertionError("flat basis lost rank modulo the rank-preserving prime")
    return _FlatState(mask, basis_indices, rref_rows, pivots)


def _extend_one_rank(
    states: dict[int, _FlatState],
    points: tuple[tuple[int, ...], ...],
    modulus: int,
) -> dict[int, _FlatState]:
    """Generate every one-rank extension of the supplied closed flats."""
    next_states: dict[int, _FlatState] = {}
    for state in states.values():
        residual_classes: dict[tuple[int, ...], tuple[int, int]] = {}
        for index, point in enumerate(points):
            if state.mask & (1 << index):
                continue
            key = _projective_key(
                _residual_mod(
                    point,
                    state.rref_rows,
                    state.pivot_columns,
                    modulus,
                ),
                modulus,
            )
            class_mask, representative = residual_classes.get(key, (0, index))
            residual_classes[key] = (class_mask | (1 << index), representative)

        for class_mask, representative in residual_classes.values():
            mask = state.mask | class_mask
            basis_indices = tuple(sorted((*state.basis_indices, representative)))
            previous = next_states.get(mask)
            if previous is not None and previous.basis_indices <= basis_indices:
                continue
            next_states[mask] = _state_from_basis(
                mask,
                basis_indices,
                points,
                modulus,
            )
    return next_states


@lru_cache(maxsize=None)
def enumerate_admissible_hyperplanes_by_flats(
    n: int, d: int
) -> FlatHyperplaneEnumeration:
    """Exhaust affine exterior-weight hyperplanes through closed flats.

    If ``F`` is a closed flat, two exterior weights outside ``F`` generate
    the same one-rank extension exactly when their residuals modulo
    ``span(F)`` are projectively equal. Grouping by that class enumerates
    every closure of rank ``rank(F) + 1`` once per parent. Induction from the
    empty flat therefore reaches every rank-``d-1`` flat, which is exactly an
    admissible affine hyperplane in the rank-``d`` augmented configuration.
    """
    modulus, bound_squared = rank_preserving_modulus(n, d)
    weights = exterior_weights(n, d)
    points = tuple(tuple(weight) + (1,) for weight in weights)

    states = {0: _FlatState(0, (), (), ())}
    flat_counts = [1]
    for _rank in range(1, d):
        states = _extend_one_rank(states, points, modulus)
        flat_counts.append(len(states))

    hyperplanes: dict[tuple[tuple[int, ...], int], AdmissibleHyperplane] = {}
    for state in states.values():
        witness = state.basis_indices
        base = weights[witness[0]]
        difference_rows = [
            [weights[index][column] - base[column] for column in range(d)]
            for index in witness[1:]
        ]
        normal = _cofactor_nullvector(difference_rows + [[1] * d])
        offset = sum(normal[index] * base[index] for index in range(d))
        h, z = _primitive_pair(normal, offset, orient=True)
        on_mask = sum(
            1 << index
            for index, weight in enumerate(weights)
            if sum(left * right for left, right in zip(h, weight, strict=True)) == z
        )
        if on_mask != state.mask:
            raise AssertionError("finite-field closure disagrees with exact hyperplane")
        candidate = AdmissibleHyperplane(h=h, z=z, witness=witness)
        previous = hyperplanes.setdefault((h, z), candidate)
        if previous != candidate and previous.witness != candidate.witness:
            raise AssertionError("distinct affine flats produced one hyperplane pair")

    return FlatHyperplaneEnumeration(
        n=n,
        d=d,
        modulus=modulus,
        hadamard_bound_squared=bound_squared,
        flat_counts_by_rank=tuple(flat_counts),
        hyperplanes=tuple(hyperplanes[key] for key in sorted(hyperplanes)),
    )
