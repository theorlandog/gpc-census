#!/usr/bin/env python3
"""Census-wide support geometry and RDM observability campaign.

Reads results/data/states.jsonl from gpc-census and independently computes:
  * singleton/pair incidence ranks;
  * support-polytope dimension and non-toric phase complexity kappa;
  * affine circuit dimension;
  * exact collision-free coherence graphs at every RDM order;
  * minimum certified reconstruction order;
  * overlap-distance and graph diagnostics.

No imports from the project generator are used.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter, defaultdict, deque
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATES = ROOT / "results/data/states.jsonl"
DEFAULT_OUT = ROOT / "results/data/census_support_geometry.json"


def parse_system(s: str) -> tuple[int, int]:
    n, d = s.strip("()").split(",")
    return int(n), int(d)


def incidence_rank(support: list[tuple[int, ...]], d: int, q: int) -> int:
    features = list(itertools.combinations(range(d), q))
    matrix = sp.Matrix([
        [int(set(feature).issubset(det)) for det in support]
        for feature in features
    ])
    return int(matrix.rank())


def inversion_sign(sequence: tuple[int, ...]) -> int:
    inv = sum(
        sequence[i] > sequence[j]
        for i in range(len(sequence))
        for j in range(i + 1, len(sequence))
    )
    return -1 if inv % 2 else 1


def channel_terms(
    support: list[tuple[int, ...]], d: int, n: int, q: int
) -> dict[tuple[tuple[int, ...], tuple[int, ...]], list[tuple[int, int, int]]]:
    """Enumerate exact unordered carrier coherences in ordered q-RDM channels."""
    index = {det: i for i, det in enumerate(support)}
    features = list(itertools.combinations(range(d), q))
    channels = {}
    remainder_size = n - q
    for P in features:
        setP = set(P)
        for Q in features:
            if P == Q:
                continue
            setQ = set(Q)
            terms: dict[tuple[int, int], int] = {}
            available = sorted(set(range(d)) - setP - setQ)
            for R in itertools.combinations(available, remainder_size):
                left = tuple(sorted(P + R))
                right = tuple(sorted(Q + R))
                if left not in index or right not in index or left == right:
                    continue
                a, b = index[left], index[right]
                sign = inversion_sign(P + R) * inversion_sign(Q + R)
                if a > b:
                    a, b = b, a
                terms[(a, b)] = terms.get((a, b), 0) + sign
            nonzero = [(a, b, s) for (a, b), s in terms.items() if s]
            if nonzero:
                channels[(P, Q)] = nonzero
    return channels


def graph_stats(vertex_count: int, edges: set[tuple[int, int]]) -> dict:
    adjacency = {i: set() for i in range(vertex_count)}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)

    components = []
    unseen = set(range(vertex_count))
    while unseen:
        root = min(unseen)
        queue = deque([root])
        component = set()
        while queue:
            u = queue.popleft()
            if u in component:
                continue
            component.add(u)
            queue.extend(adjacency[u] - component)
        unseen -= component
        components.append(component)

    diameter = None
    if len(components) == 1 and vertex_count:
        diameter = 0
        for source in range(vertex_count):
            dist = {source: 0}
            queue = deque([source])
            while queue:
                u = queue.popleft()
                for v in adjacency[u]:
                    if v not in dist:
                        dist[v] = dist[u] + 1
                        queue.append(v)
            diameter = max(diameter, max(dist.values()))

    return {
        "edge_count": len(edges),
        "component_count": len(components),
        "component_sizes": sorted((len(c) for c in components), reverse=True),
        "connected": len(components) == 1,
        "diameter": diameter,
        "cycle_rank": len(edges) - vertex_count + len(components),
        "min_degree": min((len(adjacency[i]) for i in adjacency), default=0),
        "max_degree": max((len(adjacency[i]) for i in adjacency), default=0),
    }


def pair_overlap_histogram(support: list[tuple[int, ...]], n: int) -> dict[str, int]:
    counts = Counter()
    for a, b in itertools.combinations(support, 2):
        overlap = len(set(a) & set(b))
        excitation_order = n - overlap
        counts[str(excitation_order)] += 1
    return dict(sorted(counts.items(), key=lambda kv: int(kv[0])))


def analyze_record(record: dict) -> dict:
    n, d = parse_system(record["system"])
    support = [tuple(det) for det in record["closed_form"]["support_dets"]]
    m = len(support)

    r1 = incidence_rank(support, d, 1)
    r2 = incidence_rank(support, d, 2)
    kappa = m - r1

    order_rows = {}
    minimum_order = None
    for q in range(2, n + 1):
        channels = channel_terms(support, d, n, q)
        unique_edges = {
            (terms[0][0], terms[0][1])
            for terms in channels.values()
            if len(terms) == 1
        }
        collided = sum(len(terms) > 1 for terms in channels.values())
        stats = graph_stats(m, unique_edges)
        order_rows[str(q)] = {
            **stats,
            "nonzero_oriented_channel_count": len(channels),
            "collided_oriented_channel_count": collided,
        }
        if minimum_order is None and stats["connected"]:
            minimum_order = q

    return {
        "system": record["system"],
        "index": record["index"],
        "classified": record["classified"],
        "carrier_size": m,
        "singleton_rank": r1,
        "support_polytope_dimension": r1 - 1,
        "non_toric_phase_complexity": kappa,
        "support_is_simplex": kappa == 0,
        "pair_incidence_rank": r2,
        "pair_incidence_full": r2 == m,
        "excitation_order_histogram": pair_overlap_histogram(support, n),
        "minimum_collision_free_reconstruction_order": minimum_order,
        "orders": order_rows,
    }


def summarize(rows: list[dict]) -> dict:
    by_class = defaultdict(list)
    for row in rows:
        by_class[row["classified"]].append(row)

    def aggregate(group):
        return {
            "records": len(group),
            "kappa_positive": sum(r["non_toric_phase_complexity"] > 0 for r in group),
            "kappa_zero": sum(r["non_toric_phase_complexity"] == 0 for r in group),
            "order2_certified": sum(
                r["minimum_collision_free_reconstruction_order"] == 2 for r in group
            ),
            "order2_unresolved": sum(
                r["minimum_collision_free_reconstruction_order"] != 2 for r in group
            ),
            "minimum_order_histogram": dict(Counter(
                str(r["minimum_collision_free_reconstruction_order"]) for r in group
            )),
            "kappa_histogram": dict(Counter(
                str(r["non_toric_phase_complexity"]) for r in group
            )),
        }

    unresolved = [
        r for r in rows if r["minimum_collision_free_reconstruction_order"] != 2
    ]
    positive_kappa = [r for r in rows if r["non_toric_phase_complexity"] > 0]

    return {
        "records": len(rows),
        "by_classification": {k: aggregate(v) for k, v in sorted(by_class.items())},
        "global": aggregate(rows),
        "order2_unresolved": unresolved,
        "positive_kappa_records": positive_kappa,
        "exact_implications": {
            "positive_kappa_implies_interference_in_this_census":
                all(r["classified"] == "INTERFERENCE" for r in positive_kappa),
            "positive_kappa_implies_order2_certified_in_this_census":
                all(r["minimum_collision_free_reconstruction_order"] == 2
                    for r in positive_kappa),
            "order2_unresolved_all_kappa_zero":
                all(r["non_toric_phase_complexity"] == 0 for r in unresolved),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--states",
        type=Path,
        default=DEFAULT_STATES,
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
    )
    args = parser.parse_args()

    records = [
        json.loads(line)
        for line in args.states.read_text().splitlines()
        if line.strip()
    ]
    rows = [analyze_record(record) for record in records]
    payload = {
        "schema_version": 1,
        "method": "independent exact support geometry and collision-free RDM scan",
        "records": rows,
        "summary": summarize(rows),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["summary"]["global"], indent=2))


if __name__ == "__main__":
    main()
