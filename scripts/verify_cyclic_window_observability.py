#!/usr/bin/env python3
"""Exact checks for cyclic-window and glued-carrier observability claims.

This verifier uses only the Python standard library.  It checks the finite
linear-algebra identities by exact rational elimination and the selected
collision-free graph combinatorics for representative families.
"""

from __future__ import annotations

import argparse
import itertools
import math
from fractions import Fraction
from typing import Iterable

Det = frozenset[int]


def rank_q(matrix: list[list[int]]) -> int:
    """Rank over Q by exact Gaussian elimination."""
    if not matrix:
        return 0
    a = [[Fraction(x) for x in row] for row in matrix]
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        p = a[rank][col]
        a[rank] = [x / p for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                f = a[r][col]
                a[r] = [x - f * y for x, y in zip(a[r], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def cyclic_windows(d: int, n: int, offset: int = 0) -> list[Det]:
    return [
        frozenset(offset + ((k + j) % d) for j in range(n))
        for k in range(d)
    ]


def incidence(carrier: list[Det], order: int) -> list[list[int]]:
    orbitals = sorted(set().union(*carrier))
    rows = itertools.combinations(orbitals, order)
    return [[int(set(row).issubset(det)) for det in carrier] for row in rows]


def kronecker(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    return [
        [x * y for x in ar for y in br]
        for ar in a
        for br in b
    ]


def product_carrier(left: list[Det], right: list[Det]) -> list[Det]:
    return [frozenset(a | b) for a in left for b in right]


def cross_pair_incidence(left: list[Det], right: list[Det]) -> list[list[int]]:
    lo = sorted(set().union(*left))
    ro = sorted(set().union(*right))
    cols = [(a, b) for a in left for b in right]
    return [
        [int(i in a and j in b) for a, b in cols]
        for i in lo
        for j in ro
    ]


def cycle_edges(d: int) -> set[tuple[int, int]]:
    return {tuple(sorted((k, (k + 1) % d))) for k in range(d)}


def categorical_edges(d1: int, d2: int) -> set[tuple[int, int]]:
    e1, e2 = cycle_edges(d1), cycle_edges(d2)
    out: set[tuple[int, int]] = set()
    for a, ap in e1:
        for b, bp in e2:
            for u, v in [((a, b), (ap, bp)), ((a, bp), (ap, b))]:
                iu, iv = u[0] * d2 + u[1], v[0] * d2 + v[1]
                out.add(tuple(sorted((iu, iv))))
    return out


def connected(vertex_count: int, edges: Iterable[tuple[int, int]]) -> bool:
    adj = [set() for _ in range(vertex_count)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for w in adj[v] - seen:
            seen.add(w)
            stack.append(w)
    return len(seen) == vertex_count


def collision_free_edges(carrier: list[Det]) -> set[tuple[int, int]]:
    """The carrier's actual collision-free 2-RDM graph.

    A channel (P, Q) is collision-free when exactly one spectator R of size
    N-2, disjoint from P and Q, has both P|R and Q|R in the carrier: the
    channel then exposes a single coherence. This is the definition in
    docs/rdm_observability_census.md, recomputed here on the carrier itself
    rather than on the abstract cycle or categorical product, so Proposition 3
    and the graph half of Theorem 6 are checked and not merely illustrated.
    """
    index = {tuple(sorted(det)): i for i, det in enumerate(carrier)}
    members = {tuple(sorted(det)) for det in carrier}
    orbitals = sorted(set().union(*carrier))
    size = len(next(iter(carrier)))
    edges: set[tuple[int, int]] = set()
    for p, q in itertools.combinations(itertools.combinations(orbitals, 2), 2):
        touched = set(p) | set(q)
        spare = [o for o in orbitals if o not in touched]
        found = []
        for rest in itertools.combinations(spare, size - 2):
            left = tuple(sorted(set(p) | set(rest)))
            right = tuple(sorted(set(q) | set(rest)))
            if left in members and right in members and left != right:
                found.append((index[left], index[right]))
        if len(found) == 1:
            edges.add(tuple(sorted(found[0])))
    return edges


def verify_carrier_graph(carrier: list[Det], expected: bool) -> None:
    edges = collision_free_edges(carrier)
    assert connected(len(carrier), edges) == expected, (len(carrier), len(edges))


def verify_factor_range(max_d: int) -> None:
    for d in range(3, max_d + 1):
        for n in range(2, d):
            carrier = cyclic_windows(d, n)
            singleton_rank = rank_q(incidence(carrier, 1))
            expected = d - math.gcd(d, n) + 1
            assert singleton_rank == expected, (d, n, singleton_rank, expected)
            if math.gcd(d, n) == 1:
                assert rank_q(incidence(carrier, 2)) == d


def verify_product(d1: int, n1: int, d2: int, n2: int) -> None:
    left = cyclic_windows(d1, n1, 0)
    right = cyclic_windows(d2, n2, d1)
    a1, a2 = incidence(left, 1), incidence(right, 1)
    cross = cross_pair_incidence(left, right)
    assert cross == kronecker(a1, a2)
    product = product_carrier(left, right)
    if math.gcd(d1, n1) == math.gcd(d2, n2) == 1:
        assert rank_q(cross) == d1 * d2
        assert rank_q(incidence(product, 2)) == d1 * d2
    graph = categorical_edges(d1, d2)
    assert connected(d1 * d2, graph) == (d1 % 2 == 1 or d2 % 2 == 1)
    if (
        math.gcd(d1, n1) == math.gcd(d2, n2) == 1
        and (d1 % 2 == 1 or d2 % 2 == 1)
    ):
        # Theorem 6 concludes about the carrier's own collision-free graph,
        # so check that graph, not just the categorical product it contains.
        verify_carrier_graph(product, True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-d", type=int, default=14)
    args = parser.parse_args()

    verify_factor_range(args.max_d)
    # Proposition 3 and Theorem 4 on the carriers themselves.
    for d in range(4, 10):
        for n in range(2, d):
            if math.gcd(d, n) == 1:
                verify_carrier_graph(cyclic_windows(d, n), True)
    for case in [(5, 2, 5, 2), (5, 2, 7, 3), (6, 2, 5, 2), (4, 3, 5, 2)]:
        verify_product(*case)

    print(
        "verified cyclic singleton ranks through d=%d; "
        "coprime pair-rank implication; Kronecker cross incidence; "
        "categorical-product connectivity samples; and connected "
        "collision-free graphs on the coprime cyclic and glued carriers"
        % args.max_d
    )


if __name__ == "__main__":
    main()
