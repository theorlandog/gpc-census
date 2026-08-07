#!/usr/bin/env python3
"""Parabolic survivor traversal: incremental partition maintenance vs rebuild.

Status: REFUTED as a speedup. The reformulation is correct and the increment is
valid; it simply does not pay, because survivors do not share long prefixes.

The survivors of an orbit are the arrangements of the multiset h with a
prescribed inversion number, which are the minimal-length coset representatives
of S_d / W_h at Coxeter length B. What screening needs from each is the on/below
partition of the exterior weights, and the current path rebuilds all C(d,3)
triple sums for every survivor.

A fixed-length Gray code cannot supply the increment: an adjacent transposition
changes Coxeter length by exactly one, so two distinct length-B representatives
are never adjacent. The DFS that already generates them does supply it, because
consecutive survivors share a prefix: assigning position k determines exactly
the C(k,2) triples {i,j,k} with i<j<k, and backtracking retracts them.
"""
import collections
import itertools
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.generation import orbit_canonical as oc
from gpc_census.generation.orbit_native import enumerate_orbit_survivor_taus
from gpc_census.generation.ressayre import _level_sets, admissible_candidate_from_tau
from gpc_census.generation.exterior import exterior_weights

COUNT = collections.Counter()


def rebuild_partition(arrangement, z, n, d, triples):
    """What the current path does: every triple sum, for every survivor."""
    on, below = [], []
    for index, (i, j, k) in enumerate(triples):
        COUNT["rebuild_triple_sums"] += 1
        total = arrangement[i] + arrangement[j] + arrangement[k]
        if total == z:
            on.append(index)
        elif total < z:
            below.append(index)
    return tuple(on), tuple(below)


def incremental_survivors(h, weights_below, z, n, d, triple_index):
    """Generate survivors AND their partition, maintaining it along the DFS.

    Same recursion and same prune as generate_trace_survivors, with the on/below
    partition carried as a stack: placing position k evaluates only the C(k,2)
    triples it completes, and backtracking pops them.
    """
    counts = collections.Counter(h)
    values = sorted(counts)
    prefix = []
    on_stack, below_stack = [], []

    def max_inversions(total):
        return total * (total - 1) // 2 - sum(
            c * (c - 1) // 2 for c in counts.values())

    def descend(position, remaining, budget):
        if remaining == 0:
            if budget == 0:
                yield tuple(prefix), tuple(on_stack), tuple(below_stack)
            return
        smaller = 0
        for value in values:
            available = counts[value]
            if available:
                residual = budget - smaller
                counts[value] = available - 1
                if not counts[value]:
                    del counts[value]
                if 0 <= residual <= max_inversions(remaining - 1):
                    prefix.append(value)
                    # only the triples this position completes
                    added_on = added_below = 0
                    for i in range(position):
                        for j in range(i + 1, position):
                            COUNT["incremental_triple_sums"] += 1
                            total = prefix[i] + prefix[j] + value
                            if total == z:
                                on_stack.append(triple_index[(i, j, position)])
                                added_on += 1
                            elif total < z:
                                below_stack.append(triple_index[(i, j, position)])
                                added_below += 1
                    yield from descend(position + 1, remaining - 1, residual)
                    for _ in range(added_on):
                        on_stack.pop()
                    for _ in range(added_below):
                        below_stack.pop()
                    prefix.pop()
                counts[value] = available
            smaller += available

    yield from descend(0, len(h), weights_below)


def main():
    d = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    n = 3
    triples = list(itertools.combinations(range(d), n))
    triple_index = {t: i for i, t in enumerate(triples)}
    weights = exterior_weights(n, d)

    enumeration = enumerate_orbit_survivor_taus(n, d)
    # one oriented orbit per (h, z) representative
    oriented = []
    for record in enumeration.orbits:
        for entry in record.oriented:
            oriented.append((entry.h, entry.z, entry.weights_below))

    t0 = time.time()
    rebuilt = {}
    for h, z, below_count in oriented:
        for arrangement in oc.generate_trace_survivors(h, below_count):
            rebuilt[arrangement] = rebuild_partition(arrangement, z, n, d, triples)
    t_rebuild = time.time() - t0

    t0 = time.time()
    incremental = {}
    for h, z, below_count in oriented:
        for arrangement, on, below in incremental_survivors(
                h, below_count, z, n, d, triple_index):
            incremental[arrangement] = (tuple(sorted(on)), tuple(sorted(below)))
    t_incremental = time.time() - t0

    agree = (set(rebuilt) == set(incremental)
             and all(rebuilt[a] == incremental[a] for a in rebuilt))

    print(json.dumps({
        "d": d,
        "oriented_orbits": len(oriented),
        "survivors": len(rebuilt),
        "partitions_identical": agree,
        "rebuild_triple_sums": COUNT["rebuild_triple_sums"],
        "incremental_triple_sums": COUNT["incremental_triple_sums"],
        "triple_sum_ratio": round(
            COUNT["rebuild_triple_sums"] / max(COUNT["incremental_triple_sums"], 1), 2),
        "rebuild_seconds": round(t_rebuild, 2),
        "incremental_seconds": round(t_incremental, 2),
        "predicted_ratio_d_over_3": round(d / 3, 2),
    }, indent=1), flush=True)


main()
