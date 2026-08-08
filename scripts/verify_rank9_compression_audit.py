#!/usr/bin/env python3
"""Independent replay of the rank-7 to rank-9 full-population audit.

Two things are checked, and they are different in kind.

REPRODUCTION. The rank-7 and rank-8 blocks must agree, column by column, with
the previously published ``full_rank78_compression_audit.json``: survivor
totals, the Hall split, the hole distributions, the singleton-incidence counts,
the fusion-width distributions, and the identity of the seven one-hole
Hall-feasible rows. That artifact's population gates were re-read rather than
regenerated, so this agreement is what closes the gap: the numbers now come out
of the repository's own generator.

REPLAY. Every one-hole Hall-feasible row recorded at any rank is rebuilt from
its ``tau`` and re-decided here, by exact signed sparse expansion where the
matrix is small enough and by modular evaluation otherwise. A modular
certificate is re-evaluated at the stored point and must reproduce the stored
nonzero residue, which proves the symbolic determinant is nonzero over the
integers.

This file imports neither gpc_census nor any third-party package.

Run:

    uv run scripts/verify_rank9_compression_audit.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import pathlib
from collections import deque

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
AUDIT = DATA / "rank9_compression_audit.json"
PUBLISHED = DATA / "full_rank78_compression_audit.json"

# Columns that must agree between the regenerated audit and the published one.
REPRODUCED_SCALARS = (
    ("trace_survivors", "trace_survivors"),
    ("nonstructural_survivors", "nonstructural_survivors"),
    ("hall_feasible", "hall_feasible"),
    ("hall_zero", "hall_zero"),
    ("singleton_signature_count", "singleton_signature_count"),
    ("singleton_collision_groups", "singleton_collision_groups"),
    ("max_optimal_fusion_width", "max_optimal_fusion_width"),
)


def tangent_matrix(tau, n: int = 3):
    d = len(tau)
    determinants = list(itertools.combinations(range(d), n))
    on = [det for det in determinants if sum(tau[i] for i in det) == 0]
    below = [det for det in determinants if sum(tau[i] for i in det) > 0]
    roots = [(i, j) for i in range(d) for j in range(i + 1, d) if tau[i] < tau[j]]
    row_of = {det: index for index, det in enumerate(below)}
    matrix = [[None] * len(roots) for _ in below]
    for variable, det in enumerate(on):
        support = set(det)
        for column, (i, j) in enumerate(roots):
            if i not in support or j in support:
                continue
            occupied = list(det)
            annihilation = occupied.index(i)
            occupied.pop(annihilation)
            creation = sum(value < j for value in occupied)
            sign = -1 if (annihilation + creation) % 2 else 1
            target = tuple(sorted(occupied + [j]))
            row = row_of.get(target)
            if row is not None:
                matrix[row][column] = (variable, sign)
    return on, below, roots, matrix


def maximum_matching(matrix) -> int:
    adjacency = [[j for j, e in enumerate(row) if e is not None] for row in matrix]
    left = len(adjacency)
    right = len(matrix[0]) if matrix else 0
    left_match = [-1] * left
    right_match = [-1] * right
    distance = [0] * left

    def bfs() -> bool:
        queue = deque()
        found = False
        for u in range(left):
            if left_match[u] < 0:
                distance[u] = 0
                queue.append(u)
            else:
                distance[u] = 10**9
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                w = right_match[v]
                if w < 0:
                    found = True
                elif distance[w] == 10**9:
                    distance[w] = distance[u] + 1
                    queue.append(w)
        return found

    def dfs(u: int) -> bool:
        for v in adjacency[u]:
            w = right_match[v]
            if w < 0 or (distance[w] == distance[u] + 1 and dfs(w)):
                left_match[u] = v
                right_match[v] = u
                return True
        distance[u] = 10**9
        return False

    count = 0
    while bfs():
        for u in range(left):
            if left_match[u] < 0 and dfs(u):
                count += 1
    return count


def sparse_determinant(matrix) -> dict:
    size = len(matrix)
    states = {0: {(): 1}}
    for row in matrix:
        next_states: dict[int, dict] = {}
        for used, polynomial in states.items():
            for column, entry in enumerate(row):
                if entry is None or used & (1 << column):
                    continue
                variable, edge_sign = entry
                permutation_sign = -1 if (used >> (column + 1)).bit_count() % 2 else 1
                target = next_states.setdefault(used | (1 << column), {})
                for monomial, coefficient in polynomial.items():
                    extended = tuple(sorted((*monomial, variable)))
                    updated = (
                        target.get(extended, 0)
                        + edge_sign * permutation_sign * coefficient
                    )
                    if updated:
                        target[extended] = updated
                    else:
                        target.pop(extended, None)
        states = {mask: poly for mask, poly in next_states.items() if poly}
        if not states:
            return {}
    return states.get((1 << size) - 1, {})


def determinant_modulo(matrix, values, prime: int) -> int:
    size = len(matrix)
    grid = [
        [0 if e is None else (e[1] * values[e[0]]) % prime for e in row]
        for row in matrix
    ]
    result = 1
    for column in range(size):
        pivot = next((r for r in range(column, size) if grid[r][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            grid[column], grid[pivot] = grid[pivot], grid[column]
            result = -result % prime
        result = (result * grid[column][column]) % prime
        inverse = pow(grid[column][column], prime - 2, prime)
        for r in range(column + 1, size):
            factor = (grid[r][column] * inverse) % prime
            if factor:
                for c in range(column, size):
                    grid[r][c] = (grid[r][c] - factor * grid[column][c]) % prime
    return result % prime


def level_holes(tau) -> list[int]:
    levels = sorted(set(tau))
    step = 0
    for lower, upper in zip(levels, levels[1:]):
        step = math.gcd(step, upper - lower)
    step = step or 1
    normalized = [(value - levels[0]) // step for value in levels]
    return sorted(set(range(normalized[-1] + 1)) - set(normalized))


def check_reproduction(audit: dict, published: dict) -> dict:
    """Ranks 7 and 8 must match the published audit column by column."""
    checked = []
    for system in ("(3,7)", "(3,8)"):
        if system not in audit["systems"]:
            continue
        mine = audit["systems"][system]
        theirs = published["systems"][system]
        for ours, published_key in REPRODUCED_SCALARS:
            if mine[ours] != theirs[published_key]:
                raise AssertionError(
                    f"{system}.{ours}: regenerated {mine[ours]} against "
                    f"published {theirs[published_key]}"
                )
        if mine["hole_distribution_all"] != theirs["hole_distribution_all"]:
            raise AssertionError(f"{system}: hole distribution differs")
        if (
            mine["optimal_fusion_width_distribution"]
            != theirs["optimal_fusion_width_distribution_all"]
        ):
            raise AssertionError(f"{system}: optimal width distribution differs")
        if (
            mine["hall_zero_optimal_width_distribution"]
            != theirs["optimal_fusion_width_distribution_hall_zero"]
        ):
            raise AssertionError(f"{system}: Hall-zero width distribution differs")
        mine_taus = sorted(tuple(r["tau"]) for r in mine["one_hole_hall_feasible_audit"])
        their_taus = sorted(
            tuple(r["tau"]) for r in theirs["one_hole_hall_feasible_exact_audit"]
        )
        if mine_taus != their_taus:
            raise AssertionError(f"{system}: one-hole row identities differ")
        checked.append(system)
    return {"systems_reproduced": checked, "columns_per_system": len(REPRODUCED_SCALARS) + 4}


def check_replay(audit: dict) -> dict:
    """Rebuild and re-decide every recorded one-hole Hall-feasible row."""
    exact_zero = exact_nonzero = modular_nonzero = undecided = 0
    for system, record in sorted(audit["systems"].items()):
        n = record["particle_number"]
        for row in record["one_hole_hall_feasible_audit"]:
            tau = tuple(row["tau"])
            holes = level_holes(tau)
            if holes != row["missing_levels"] or len(holes) != 1:
                raise AssertionError(f"{system}: stored hole structure is wrong")
            on, below, roots, matrix = tangent_matrix(tau, n)
            if len(below) != row["matrix_size"] or len(on) != row["on_variables"]:
                raise AssertionError(f"{system}: stored matrix shape is wrong")
            if len(below) != len(roots):
                raise AssertionError(f"{system}: tangent matrix is not square")
            if maximum_matching(matrix) != len(below):
                raise AssertionError(f"{system}: row is not Hall-feasible on replay")

            verdict = row["determinant"]
            if verdict["method"] == "exact_sparse_expansion":
                coefficients = sparse_determinant(matrix)
                if bool(coefficients) != verdict["nonzero"]:
                    raise AssertionError(f"{system}: exact determinant verdict differs")
                if len(coefficients) != verdict["coefficient_count"]:
                    raise AssertionError(f"{system}: exact coefficient count differs")
                if coefficients:
                    exact_nonzero += 1
                else:
                    exact_zero += 1
            elif verdict["nonzero"]:
                residue = determinant_modulo(
                    matrix, verdict["variable_values"], verdict["prime"]
                )
                if residue != verdict["residue"] or residue == 0:
                    raise AssertionError(f"{system}: modular certificate does not replay")
                modular_nonzero += 1
            else:
                undecided += 1
    return {
        "exact_zero": exact_zero,
        "exact_nonzero": exact_nonzero,
        "modular_nonzero_certificates": modular_nonzero,
        "undecided_modular_samples": undecided,
        "note": (
            "Modular nonzero is a proof of nonvanishing over the integers. An "
            "undecided row is neither proved zero nor proved nonzero."
        ),
    }


def check_falsifiers(audit: dict) -> dict:
    """The falsifiers the frontier document named, evaluated on this population."""
    multi_hole_feasible = 0
    multi_hole_survivors = 0
    collisions = 0
    for record in audit["systems"].values():
        multi_hole_feasible += record["multi_hole_hall_feasible"]
        multi_hole_survivors += record["multi_hole_survivors"]
        collisions += record["singleton_collision_groups"]
        total = sum(int(v) for v in record["hole_distribution_all"].values())
        if total != record["trace_survivors"]:
            raise AssertionError("hole distribution does not cover the survivors")
        if record["hall_feasible"] + record["hall_zero"] != record["trace_survivors"]:
            raise AssertionError("Hall split does not cover the survivors")
    return {
        "multi_hole_hall_feasible": multi_hole_feasible,
        "multi_hole_survivors": multi_hole_survivors,
        "singleton_collision_groups": collisions,
        "at_most_one_hole_holds": multi_hole_feasible == 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=pathlib.Path, default=AUDIT)
    arguments = parser.parse_args()

    audit = json.loads(arguments.audit.read_text())
    published = json.loads(PUBLISHED.read_text())

    reproduction = check_reproduction(audit, published)
    replay = check_replay(audit)
    falsifiers = check_falsifiers(audit)

    print("rank-7 to rank-9 audit verifier: PASS")
    print(f"  reproduced published systems: {reproduction['systems_reproduced']}")
    print(
        "  one-hole rows replayed: "
        f"{replay['exact_zero']} exactly zero, "
        f"{replay['exact_nonzero']} exactly nonzero, "
        f"{replay['modular_nonzero_certificates']} modular nonzero, "
        f"{replay['undecided_modular_samples']} undecided"
    )
    print(f"  multi-hole Hall-feasible rows: {falsifiers['multi_hole_hall_feasible']}")
    print(f"  singleton collision groups: {falsifiers['singleton_collision_groups']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
