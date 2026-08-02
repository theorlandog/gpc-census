#!/usr/bin/env python3
"""Exact combinatorics for observable-manifold tree charts.

The mathematical theorem is support independent.  This script audits the
boundary stratification of selected coherence trees: on an active support A,
the tree chart leaves one relative phase per connected component, hence a
(c(A)-1)-torus ambiguity when the induced forest has c(A) components.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter, deque
from pathlib import Path


def path_tree(m: int) -> list[tuple[int, int]]:
    return [(i, i + 1) for i in range(m - 1)]


def star_tree(m: int) -> list[tuple[int, int]]:
    return [(0, i) for i in range(1, m)]


def balanced_tree(m: int) -> list[tuple[int, int]]:
    return [((i - 1) // 2, i) for i in range(1, m)]


def component_count(active: frozenset[int], edges: list[tuple[int, int]]) -> int:
    if not active:
        return 0
    adjacency = {v: set() for v in active}
    for u, v in edges:
        if u in active and v in active:
            adjacency[u].add(v)
            adjacency[v].add(u)
    unseen = set(active)
    count = 0
    while unseen:
        count += 1
        root = min(unseen)
        unseen.remove(root)
        queue = deque([root])
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if v in unseen:
                    unseen.remove(v)
                    queue.append(v)
    return count


def audit_tree(m: int, edges: list[tuple[int, int]], name: str) -> dict:
    by_support_size: dict[str, dict] = {}
    total = Counter()
    worst = {"ambiguity_dimension": -1, "active_support": []}
    for size in range(1, m + 1):
        counts = Counter()
        connected = 0
        for subset in itertools.combinations(range(m), size):
            active = frozenset(subset)
            components = component_count(active, edges)
            ambiguity = components - 1
            counts[ambiguity] += 1
            total[ambiguity] += 1
            connected += ambiguity == 0
            if ambiguity > worst["ambiguity_dimension"]:
                worst = {
                    "ambiguity_dimension": ambiguity,
                    "active_support": list(subset),
                }
        by_support_size[str(size)] = {
            "strata": sum(counts.values()),
            "connected_chart_strata": connected,
            "ambiguity_dimension_counts": {
                str(k): v for k, v in sorted(counts.items())
            },
        }
    return {
        "name": name,
        "vertices": m,
        "edges": [list(edge) for edge in edges],
        "full_support_dimension_real": 2 * m - 2,
        "full_support_chart_connected": component_count(frozenset(range(m)), edges) == 1,
        "boundary": by_support_size,
        "all_nonempty_strata_ambiguity_counts": {
            str(k): v for k, v in sorted(total.items())
        },
        "worst_boundary_stratum": worst,
    }


def build(max_vertices: int = 10) -> dict:
    audits = []
    for m in range(2, max_vertices + 1):
        audits.extend(
            [
                audit_tree(m, path_tree(m), "path"),
                audit_tree(m, star_tree(m), "star"),
                audit_tree(m, balanced_tree(m), "balanced"),
            ]
        )
    return {
        "schema_version": 1,
        "theorem": {
            "full_support_geometry": "interior simplex Delta^(m-1) times torus T^(m-1)",
            "real_dimension": "2m-2",
            "selected_tree_relations": "|z_uv|^2 = w_u w_v for each tree edge",
            "boundary_fiber": "T^(c(A)-1), where c(A) is the number of components of the induced tree on active support A",
            "rank_defect": "c(A)-1 in relative-phase directions",
        },
        "audits": audits,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-vertices", type=int, default=10)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/data/observable_manifold_geometry.json"),
    )
    args = parser.parse_args()
    artifact = build(args.max_vertices)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.out}: {len(artifact['audits'])} exact tree audits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
