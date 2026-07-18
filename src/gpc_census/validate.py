"""Structural validation invariants for vertex sets. Zero dependencies.

These checks caught two real errors during the original census; every
pipeline output should pass through them.
"""
from __future__ import annotations

from fractions import Fraction

Vertex = tuple[Fraction, ...]


def check_physical(verts: list[Vertex], n: int) -> list[str]:
    """Trace, ordering, and 0 <= lambda <= 1 for every vertex."""
    errs = []
    for v in verts:
        if sum(v) != n:
            errs.append(f"trace != {n}: {v}")
        if any(v[i] < v[i + 1] for i in range(len(v) - 1)):
            errs.append(f"not descending: {v}")
        if v and (v[0] > 1 or v[-1] < 0):
            errs.append(f"occupation outside [0,1]: {v}")
    return errs


def check_embedding(lower: list[Vertex], higher: list[Vertex]) -> list[str]:
    """Every higher-rank vertex with a trailing zero must be a lower-rank vertex."""
    lowset = set(lower)
    errs = []
    for v in higher:
        if v[-1] == 0 and v[:-1] not in lowset:
            errs.append(f"trailing-zero vertex not in lower rank: {v}")
    return errs


def check_selfdual(verts: list[Vertex], n: int, d: int) -> list[str]:
    """For half filling (2n == d) the vertex set must be complement-closed."""
    if 2 * n != d:
        return []
    vset = set(verts)
    errs = []
    for v in verts:
        dual = tuple(sorted((1 - x for x in v), reverse=True))
        if dual not in vset:
            errs.append(f"complement not a vertex: {v}")
    return errs
