"""Core primitives for the GPC census pipeline."""

from __future__ import annotations

from itertools import combinations
from typing import Iterator


def slater_vertices(d: int, n: int) -> Iterator[tuple[int, ...]]:
    """Yield occupation-number vectors of Slater determinants for ``n`` fermions in ``d`` orbitals.

    These 0/1 vectors with exactly ``n`` ones are the vertices of the
    hypersimplex Delta(d, n), the classical Pauli polytope, and the
    starting point for a census of extremal states of the full moment
    polytope.
    """
    if d < 0 or n < 0:
        raise ValueError("d and n must be non-negative")
    if n > d:
        raise ValueError(f"cannot place {n} fermions in {d} orbitals")
    for occupied in combinations(range(d), n):
        vec = [0] * d
        for i in occupied:
            vec[i] = 1
        yield tuple(vec)
