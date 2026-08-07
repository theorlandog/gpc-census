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
    rows: tuple[tuple[int, ...], ...], modulus: int, inverse_table=None
) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]:
    """Return reduced row-echelon form over the configured prime field.

    ``inverse_table`` removes one modular exponentiation per pivot. The
    reduction runs once per parent flat, which is 1.2 million times at rank 8,
    so the exponentiations are worth removing even though each is cheap.
    """
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
        head = matrix[pivot_row][column]
        inverse = (int(inverse_table[head]) if inverse_table is not None
                   else pow(head, modulus - 2, modulus))
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


def _inverse_table(modulus: int):
    """Modular inverses ``1..modulus-1`` as an array, or None if too large.

    The inverse of a residual's leading entry is needed once per candidate
    point, which is the single hottest scalar operation in the induction
    (1.37 million ``pow`` calls at rank 7 alone). Every modulus this generator
    actually selects is small enough to tabulate: the Hadamard bound needs
    ``prime^2 > 4^d``, so 65537 serves every rank up to 15.
    """
    if modulus > 1 << 21:
        return None
    import numpy as np

    table = np.zeros(modulus, dtype=np.int64)
    table[1:] = [pow(value, modulus - 2, modulus) for value in range(1, modulus)]
    return table


def _incidence_masks(rows_out, weights, d: int) -> list[int]:
    """Which weights lie on each hyperplane, as bitmasks.

    This is the exact test ``<h, omega> == z`` for every (hyperplane, weight)
    pair, which is 9.3 million dot products at rank 8. As an integer matrix
    product it is one matmul, and it stays exact: the entries are bounded well
    inside int64 and the bound is asserted rather than assumed.
    """
    if not rows_out:
        return []
    try:
        import numpy as np
    except ModuleNotFoundError:                       # pragma: no cover
        return [
            sum(
                1 << index
                for index, weight in enumerate(weights)
                if sum(a * b for a, b in zip(h, weight, strict=True)) == z
            )
            for _, _, h, z in rows_out
        ]

    normals = np.array([h for _, _, h, _ in rows_out], dtype=np.int64)
    offsets = np.array([z for _, _, _, z in rows_out], dtype=np.int64)
    weight_matrix = np.asarray(weights, dtype=np.int64).T
    largest = int(np.abs(normals).max(initial=0))
    if largest * d >= 1 << 62:                        # never reached in practice
        raise AssertionError("hyperplane incidence product would overflow int64")
    on = (normals @ weight_matrix) == offsets[:, None]
    packed = np.packbits(on, axis=1, bitorder="little")
    return [int.from_bytes(row.tobytes(), "little") for row in packed]


def _extend_one_rank(
    states: dict[int, tuple[int, ...]],
    points: tuple[tuple[int, ...], ...],
    modulus: int,
    points_array=None,
    inverse_table=None,
) -> dict[int, tuple[int, ...]]:
    """Generate every one-rank extension of the supplied closed flats.

    A level is stored as ``mask -> basis_indices`` and NOT as a materialised
    row reduction. Carrying the RREF per flat costs about 870 of the 950 bytes
    a stored flat used to take, and it is redundant: the reduction is a
    function of the basis and is needed only while a flat is being extended.
    Recomputing it once per PARENT also replaces one reduction per CHILD, and
    children outnumber parents heavily at the wide middle ranks.

    The residual-and-key step is done for ALL candidate points of a parent at
    once. It is the dominant cost of the whole generator (59.8 million
    residuals at rank 8, about 50 per parent), and it is exactly a matrix
    operation: because the reduction is in reduced row-echelon form the
    residual is ``point - point[pivots] @ rref``, one matmul for the entire
    candidate block. The projective normalisation and the grouping key are
    vectorised with it.

    Rank validation is preserved: every flat is reduced when it is extended,
    and the final rank, which is never extended, is checked exactly by the
    hyperplane closure test in the caller.
    """
    if points_array is None or inverse_table is None:      # scalar fallback
        return _extend_one_rank_scalar(states, points, modulus)

    import numpy as np

    point_count, width = points_array.shape
    all_indices = np.arange(point_count)
    byte_width = (point_count + 7) // 8
    next_states: dict[int, tuple[int, ...]] = {}
    for mask, basis_indices in states.items():
        rows = tuple(points[index] for index in basis_indices)
        rref_rows, pivot_columns = _rref_mod(rows, modulus, inverse_table)
        if len(rref_rows) != len(basis_indices):
            raise AssertionError("flat basis lost rank modulo the rank-preserving prime")

        # The mask has one bit per exterior weight, so it exceeds 64 bits from
        # rank 9 on (84 weights). Shifting it by a numpy array would overflow a
        # C long; unpacking its bytes keeps the arbitrary-precision semantics.
        if mask:
            bits = np.unpackbits(
                np.frombuffer(mask.to_bytes(byte_width, "little"), dtype=np.uint8),
                bitorder="little",
            )[:point_count]
            candidates = all_indices[bits == 0]
        else:
            candidates = all_indices
        if candidates.size == 0:
            continue
        block = points_array[candidates]
        if pivot_columns:
            reduction = np.asarray(rref_rows, dtype=np.int64)
            residual = (block - block[:, pivot_columns] @ reduction) % modulus
        else:
            residual = block % modulus

        nonzero = residual != 0
        if not nonzero.any(axis=1).all():
            raise AssertionError("an element outside a closed flat had zero residual")
        leading = residual[np.arange(residual.shape[0]), nonzero.argmax(axis=1)]
        residual = (residual * inverse_table[leading][:, None]) % modulus

        raw = residual.astype("<u4").tobytes()
        stride = width * 4
        residual_classes: dict[bytes, tuple[int, int]] = {}
        for offset, index in enumerate(candidates.tolist()):
            key = raw[offset * stride:(offset + 1) * stride]
            class_mask, representative = residual_classes.get(key, (0, index))
            residual_classes[key] = (class_mask | (1 << index), representative)

        for class_mask, representative in residual_classes.values():
            child_mask = mask | class_mask
            child_basis = tuple(sorted((*basis_indices, representative)))
            previous = next_states.get(child_mask)
            if previous is not None and previous <= child_basis:
                continue
            next_states[child_mask] = child_basis
    return next_states


def iter_one_rank_children(
    states: dict[int, tuple[int, ...]],
    points: tuple[tuple[int, ...], ...],
    modulus: int,
    points_array=None,
    inverse_table=None,
):
    """Yield ``(child_mask, child_basis)`` one parent at a time.

    ``_extend_one_rank`` accumulates every child of every parent into one dict
    before anything downstream sees them, which is the right shape when the
    caller wants the whole level. It is the wrong shape for the ORBIT
    enumerator, which is about to collapse those children onto a few thousand
    canonical forms and therefore never needs the intermediate set: at (3,11)
    the level-9 orbits alone expand to millions of children, and holding them
    is what exhausted memory. Streaming caps the caller's peak at the number of
    ORBITS in the next level rather than the number of children.

    A child may be yielded more than once, by different parents or by the same
    parent through different points. Deduplication is the caller's business.
    """
    for parent, basis_indices in states.items():
        one = {parent: basis_indices}
        yield from _extend_one_rank(
            one, points, modulus, points_array, inverse_table).items()


def _extend_one_rank_scalar(
    states: dict[int, tuple[int, ...]],
    points: tuple[tuple[int, ...], ...],
    modulus: int,
) -> dict[int, tuple[int, ...]]:
    """Reference implementation of one rank extension, without numpy.

    Kept as the definition the vectorised path is checked against, and used
    when the modulus is too large to tabulate inverses.
    """
    next_states: dict[int, tuple[int, ...]] = {}
    for mask, basis_indices in states.items():
        rows = tuple(points[index] for index in basis_indices)
        rref_rows, pivot_columns = _rref_mod(rows, modulus)
        if len(rref_rows) != len(basis_indices):
            raise AssertionError("flat basis lost rank modulo the rank-preserving prime")

        residual_classes: dict[tuple[int, ...], tuple[int, int]] = {}
        for index, point in enumerate(points):
            if mask & (1 << index):
                continue
            key = _projective_key(
                _residual_mod(point, rref_rows, pivot_columns, modulus),
                modulus,
            )
            class_mask, representative = residual_classes.get(key, (0, index))
            residual_classes[key] = (class_mask | (1 << index), representative)

        for class_mask, representative in residual_classes.values():
            child_mask = mask | class_mask
            child_basis = tuple(sorted((*basis_indices, representative)))
            previous = next_states.get(child_mask)
            if previous is not None and previous <= child_basis:
                continue
            next_states[child_mask] = child_basis
    return next_states


def _level_path(checkpoint_dir, n: int, d: int, rank: int):
    from pathlib import Path
    return Path(checkpoint_dir) / f"flats_{n}_{d}_rank{rank}.txt"


def _save_level(checkpoint_dir, n: int, d: int, rank: int,
                states: dict[int, tuple[int, ...]]) -> None:
    """Write one closed-flat level as ``hexmask:i1,i2,...`` lines.

    Deliberately a plain text format rather than pickle: a checkpoint is data,
    and loading it must not be able to execute anything.
    """
    path = _level_path(checkpoint_dir, n, d, rank)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".partial")
    with tmp.open("w") as handle:
        for mask, basis in states.items():
            handle.write(f"{mask:x}:{','.join(str(i) for i in basis)}\n")
    tmp.replace(path)          # atomic, so a killed run leaves no half level


def _load_level(checkpoint_dir, n: int, d: int, rank: int):
    path = _level_path(checkpoint_dir, n, d, rank)
    if not path.exists():
        return None
    states: dict[int, tuple[int, ...]] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            mask_hex, _, basis = line.partition(":")
            states[int(mask_hex, 16)] = tuple(
                int(part) for part in basis.split(",") if part)
    return states


@lru_cache(maxsize=None)
def enumerate_admissible_hyperplanes_by_flats(
    n: int, d: int, checkpoint_dir: str | None = None
) -> FlatHyperplaneEnumeration:
    """Exhaust affine exterior-weight hyperplanes through closed flats.

    If ``F`` is a closed flat, two exterior weights outside ``F`` generate
    the same one-rank extension exactly when their residuals modulo
    ``span(F)`` are projectively equal. Grouping by that class enumerates
    every closure of rank ``rank(F) + 1`` once per parent. Induction from the
    empty flat therefore reaches every rank-``d-1`` flat, which is exactly an
    admissible affine hyperplane in the rank-``d`` augmented configuration.

    ``checkpoint_dir`` makes the level induction RESUMABLE. Each completed
    level is written out and reloaded on a later call, so a rank whose cost is
    hours rather than gigabytes can be finished across several runs. The
    enumeration is deterministic, so a resumed run and a fresh one produce the
    same lattice; ``test_checkpoint_round_trip_is_identical`` pins that.
    """
    modulus, bound_squared = rank_preserving_modulus(n, d)
    weights = exterior_weights(n, d)
    points = tuple(tuple(weight) + (1,) for weight in weights)

    inverse_table = _inverse_table(modulus)
    points_array = None
    if inverse_table is not None:
        import numpy as np

        points_array = np.asarray(points, dtype=np.int64)

    states: dict[int, tuple[int, ...]] = {0: ()}
    flat_counts = [1]
    for rank in range(1, d):
        restored = _load_level(checkpoint_dir, n, d, rank) if checkpoint_dir else None
        if restored is not None:
            states = restored
        else:
            states = _extend_one_rank(
                states, points, modulus, points_array, inverse_table)
            if checkpoint_dir:
                _save_level(checkpoint_dir, n, d, rank, states)
        flat_counts.append(len(states))

    rows_out: list[tuple[int, tuple[int, ...], tuple[int, ...], int]] = []
    for state_mask, witness in states.items():
        base = weights[witness[0]]
        difference_rows = [
            [weights[index][column] - base[column] for column in range(d)]
            for index in witness[1:]
        ]
        normal = _cofactor_nullvector(difference_rows + [[1] * d])
        offset = sum(normal[index] * base[index] for index in range(d))
        h, z = _primitive_pair(normal, offset, orient=True)
        rows_out.append((state_mask, witness, h, z))

    on_masks = _incidence_masks(rows_out, weights, d)

    hyperplanes: dict[tuple[tuple[int, ...], int], AdmissibleHyperplane] = {}
    for (state_mask, witness, h, z), on_mask in zip(rows_out, on_masks, strict=True):
        if on_mask != state_mask:
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
