"""Direct cocircuit enumeration of admissible exterior-weight hyperplanes.

``admissible_flats`` builds the whole lattice of closed affine flats and reads
the admissible hyperplanes off the top level; ``orbit_canonical`` runs the same
induction while holding only orbit representatives. Both reach codimension one
THROUGH every lower rank. This module reaches it directly.

WHY THAT IS LEGITIMATE, and not a heuristic shortcut. Every exterior weight
satisfies ``1 . chi_I = n``, so the affine hull of the slice misses the origin
and

    linear_rank(B_0) = affine_rank(B_0) + 1.

The Ressayre admissibility condition the generator enforces is that the zero
family has affine rank ``d-2``, hence linear rank ``d-1``, hence its span is a
linear hyperplane, which is exactly ``tau^perp``. No omitted weight can lie in
that span, since it would then also annihilate ``tau``. So ``B_0`` is a maximal
coplanar subset of linear rank ``d-1``: the zero set of a cocircuit of the
exterior-weight configuration. Conversely any maximal coplanar subset of linear
rank ``d-1`` has affine rank ``d-2``, so its primitive annihilator is
admissible. Admissible-hyperplane enumeration and cocircuit-zero-set
enumeration are the same problem.

The search is lexicographic subset reverse search with the two prunings that
make it finite in practice: a node whose omitted lexicographic predecessor
already lies in the current span cannot be the canonical representative of any
completion below it, and a node whose remaining suffix cannot lift the span to
rank ``d-1`` has no descendant at all. Ranks are computed over the same
rank-preserving prime ``admissible_flats`` uses, so modular rank equals
rational rank; the recovered normal is an exact integer cofactor vector and
every emitted zero set is re-derived from it in exact integer arithmetic before
the enumeration returns.

WHAT THIS IS FOR. It is an INDEPENDENT low-rank backend that shares no code
path with the flat lattice, so agreement at ``d = 5..8`` is a genuine second
derivation of those populations. It is unsymmetric, so it emits all 166,420
rank-8 hyperplanes rather than the 56 orbits, and from rank 9 on
``orbit_canonical`` is the enumerator to use. ``docs/cocircuit_chow_decoder.md``
carries the measured comparison.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .admissible_flats import rank_preserving_modulus
from .exterior import exterior_weights
from .ressayre import _cofactor_nullvector, _primitive_pair


class NodeLimitExceeded(RuntimeError):
    """The reverse search hit its node budget before finishing."""


@dataclass(frozen=True)
class CocircuitEnumeration:
    """Admissible hyperplanes found by direct cocircuit reverse search."""

    n: int
    d: int
    modulus: int
    weight_count: int
    node_count: int
    normals: tuple[tuple[int, ...], ...]
    zero_sets: tuple[tuple[int, ...], ...]

    @property
    def hyperplane_count(self) -> int:
        return len(self.normals)

    def as_dict(self) -> dict[str, object]:
        """Return a compact deterministic proof manifest."""
        return {
            "particle_number": self.n,
            "rank": self.d,
            "rank_preserving_modulus": self.modulus,
            "exterior_weights": self.weight_count,
            "reverse_search_nodes": self.node_count,
            "admissible_hyperplane_count": self.hyperplane_count,
        }

    def affine_pairs(self) -> tuple[tuple[tuple[int, ...], int], ...]:
        """Return each normal in the ``(h, z)`` convention of the flat lattice.

        On ``sum(lambda) = n`` the homogeneous form ``tau . omega = 0`` and the
        centered form ``h . omega = z`` with ``h_i = S - d tau_i`` and
        ``z = n S`` describe one hyperplane, so this is the exact bridge that
        lets the two enumerators be compared pair for pair and not only by
        count.
        """
        return tuple(
            _primitive_pair(
                tuple(sum(tau) - self.d * value for value in tau),
                self.n * sum(tau),
                orient=True,
            )
            for tau in self.normals
        )


def _primitive(vector) -> tuple[int, ...]:
    divisor = 0
    for value in vector:
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:                                      # pragma: no cover
        raise ValueError("zero cocircuit normal")
    values = tuple(value // divisor for value in vector)
    lead = next(value for value in values if value)
    return values if lead > 0 else tuple(-value for value in values)


def _suffix_rank(block, modulus: int, needed: int) -> int:
    """Rank of a residual block modulo the prime, capped at ``needed``.

    The only question asked of it is whether the remaining suffix can still
    lift the span to rank ``d-1``, so the elimination stops as soon as the
    answer is yes.
    """
    import numpy as np

    if needed <= 0:
        return 0
    block = block.copy()
    found = 0
    for column in range(block.shape[1]):
        nonzero = np.flatnonzero(block[:, column])
        if nonzero.size == 0:
            continue
        head = int(nonzero[0])
        pivot_row = (
            block[head] * pow(int(block[head, column]), modulus - 2, modulus)
        ) % modulus
        others = nonzero[1:]
        if others.size:
            block[others] = (
                block[others] - block[others, column][:, None] * pivot_row
            ) % modulus
        block[head] = 0
        found += 1
        if found >= needed:
            return found
    return found


def enumerate_cocircuit_hyperplanes(
    n: int, d: int, *, node_limit: int | None = None
) -> CocircuitEnumeration:
    """Enumerate admissible hyperplanes as cocircuit zero sets.

    Raises :class:`NodeLimitExceeded` rather than returning a truncated
    enumeration: a partial reverse search is not a population count and must
    never be mistaken for one.
    """
    import numpy as np

    modulus, _bound_squared = rank_preserving_modulus(n, d)
    weights = exterior_weights(n, d)
    columns = np.asarray(weights, dtype=np.int64)
    weight_count = len(weights)
    target_rank = d - 1

    node_count = 0
    normals: list[tuple[int, ...]] = []
    zero_sets: list[tuple[int, ...]] = []

    def recurse(selected, on_mask, maximum, residuals, in_span, rank, independent):
        nonlocal node_count

        # Prune 1: a smaller omitted weight already in the span is forced into
        # every maximal completion below this node, so the node cannot be the
        # lexicographic representative of any of them.
        if maximum > 0 and bool(np.any(in_span[:maximum] & ~on_mask[:maximum])):
            return 0
        # Prune 2: the suffix must still be able to lift the span to d-1.
        if rank < target_rank:
            suffix = residuals[maximum + 1 :]
            if suffix.shape[0] < target_rank - rank:
                return 0
            if _suffix_rank(suffix, modulus, target_rank - rank) < target_rank - rank:
                return 0

        node_count += 1
        if node_limit is not None and node_count > node_limit:
            raise NodeLimitExceeded(f"reverse search exceeded node_limit={node_limit}")

        descendants = 0
        for index in range(maximum + 1, weight_count):
            dependent = bool(in_span[index])
            if dependent:
                child_residuals, child_in_span = residuals, in_span
                child_rank, child_independent = rank, independent
            else:
                if rank == target_rank:
                    continue
                # One rank-one update reduces every residual modulo the new
                # span, so a child never re-reduces the whole configuration.
                row = residuals[index]
                pivot = int(np.flatnonzero(row)[0])
                row = row * pow(int(row[pivot]), modulus - 2, modulus) % modulus
                child_residuals = (
                    residuals - residuals[:, pivot][:, None] * row
                ) % modulus
                child_in_span = ~child_residuals.any(axis=1)
                child_rank = rank + 1
                child_independent = (*independent, index)

            on_mask[index] = True
            selected.append(index)
            descendants += recurse(
                selected,
                on_mask,
                index,
                child_residuals,
                child_in_span,
                child_rank,
                child_independent,
            )
            selected.pop()
            on_mask[index] = False

            # A dependent weight already lies in the span, so every maximal
            # completion of this branch contains it: one recursion suffices.
            if dependent:
                break

        if descendants:
            return descendants
        if rank == target_rank and not bool(np.any(in_span & ~on_mask)):
            normals.append(
                _primitive(
                    _cofactor_nullvector([list(weights[i]) for i in independent])
                )
            )
            zero_sets.append(tuple(selected))
            return 1
        return 0

    root_residuals = columns % modulus
    recurse(
        [],
        np.zeros(weight_count, dtype=bool),
        -1,
        root_residuals,
        ~root_residuals.any(axis=1),
        0,
        (),
    )

    _verify_zero_sets(n, d, normals, zero_sets)
    return CocircuitEnumeration(
        n=n,
        d=d,
        modulus=modulus,
        weight_count=weight_count,
        node_count=node_count,
        normals=tuple(normals),
        zero_sets=tuple(zero_sets),
    )


def _verify_zero_sets(n: int, d: int, normals, zero_sets) -> None:
    """Recheck every emitted zero set in exact integer arithmetic.

    The search runs modulo a rank-preserving prime, chosen above the Hadamard
    bound so the ranks are exact. The statement actually shipped, though, is
    that a specific integer normal annihilates a specific set of weights and no
    other. That is cheap to recheck outright, so it is rechecked outright
    rather than inferred from the choice of prime.
    """
    import numpy as np

    if not normals:
        return
    weights = np.asarray(exterior_weights(n, d), dtype=np.int64)
    grades = np.asarray(normals, dtype=np.int64) @ weights.T
    if int(np.abs(grades).max(initial=0)) >= 1 << 62:     # pragma: no cover
        raise AssertionError("cocircuit grade product would overflow int64")
    for row, zero_set in enumerate(zero_sets):
        recovered = tuple(int(value) for value in np.flatnonzero(grades[row] == 0))
        if recovered != zero_set:
            raise AssertionError(
                "modular cocircuit zero set disagrees with exact arithmetic"
            )
        if len(zero_set) < d - 1:                         # pragma: no cover
            raise AssertionError("an emitted cocircuit is below rank d-1")


def trace_survivors(enumeration: CocircuitEnumeration):
    """Yield ``(tau, positive_count)`` for oriented normals passing the trace test.

    The test is ``#{omega : tau . omega > 0} == #{i < j : tau_j > tau_i}``. In
    the centered convention ``h_i = S - d tau_i``, ``z = n S`` those two counts
    are exactly the repository's ``below_weight_count`` and
    ``negative_root_count``, so this is the same exact Ressayre trace test
    written in the homogeneous coordinates the cocircuit search produces.
    """
    import numpy as np

    if not enumeration.normals:
        return
    d = enumeration.d
    weights = np.asarray(exterior_weights(enumeration.n, d), dtype=np.int64)
    base = np.asarray(enumeration.normals, dtype=np.int64)
    oriented = np.concatenate([base, -base])
    positive = (oriented @ weights.T > 0).sum(axis=1)
    inversions = np.zeros(oriented.shape[0], dtype=np.int64)
    for left in range(d):
        for right in range(left + 1, d):
            inversions += oriented[:, right] > oriented[:, left]
    for row in np.flatnonzero(positive == inversions):
        yield tuple(int(value) for value in oriented[row]), int(positive[row])
