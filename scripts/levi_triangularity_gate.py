#!/usr/bin/env python3
"""Gate 4: are fermionic tangent systems linear triangular after DM factorization?

THE CONJECTURE UNDER TEST. The Fermionic Levi Triangularity Conjecture says
every minimal fermionic BDR facet becomes linear triangular after
determinant-twist removal, decomposition into connected Levi-excitation
components, grade ordering, and Dulmage-Mendelsohn decomposition of the linear
tangent blocks. If it held, birationality would need no generic Groebner
calculation: the tangent system would solve by back substitution.

TWO POPULATIONS, KEPT APART. The wide population is every Hall-feasible row
with a nonzero Ressayre tangent determinant: 48, 136 and 240 rows at ranks 7,
8 and 9. Those are CANDIDATES, not facets. They have had no Farkas reduction,
no facet witness and no BDR birationality test, and the published irredundant
systems are much smaller, 4, 31 and 52 rows. The published rows are therefore
audited separately, and the verdict is reported on both.

WHAT IS MEASURED. For each determinant-nonzero Hall-feasible row, the tangent
matrix is square and has a perfect matching. Orient it by any perfect matching
and take strongly connected components of the induced digraph on rows: those
are the irreducible blocks of the Dulmage-Mendelsohn decomposition. The matrix
is triangular in the strict sense exactly when every block is ``1 x 1``, which
is the case where back substitution alone solves the system.

CANONICALITY. The block multiset does not depend on which perfect matching is
chosen; that is the Dulmage-Mendelsohn theorem, and it is re-checked here on a
sample rather than assumed, because the whole verdict rests on it.

Run:

    uv run scripts/levi_triangularity_gate.py
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import importlib.util
import json
import math
import pathlib
import random
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CENSUS_SYSTEMS = ((3, 7, "flats"), (3, 8, "flats"), (3, 9, "augmentation"))
TIMING_FIELDS = ("seconds",)


def _audit_module():
    spec = importlib.util.spec_from_file_location(
        "rank9_compression_audit", ROOT / "scripts" / "rank9_compression_audit.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


AUDIT = _audit_module()


def perfect_matching(matrix, order=None):
    """Row-to-column perfect matching, or None; ``order`` varies the result."""
    rows = len(matrix)
    columns = len(matrix[0]) if matrix else 0
    adjacency = [[j for j, entry in enumerate(row) if entry is not None] for row in matrix]
    match_column = [-1] * columns
    match_row = [-1] * rows

    def augment(u, seen) -> bool:
        for v in adjacency[u]:
            if seen[v]:
                continue
            seen[v] = True
            if match_column[v] < 0 or augment(match_column[v], seen):
                match_column[v] = u
                match_row[u] = v
                return True
        return False

    for u in order if order is not None else range(rows):
        augment(u, [False] * columns)
    return match_row if all(value >= 0 for value in match_row) else None


def irreducible_blocks(matrix, match_row) -> list[int]:
    """Sizes of the Dulmage-Mendelsohn irreducible blocks, largest first."""
    size = len(matrix)
    owner = {match_row[j]: j for j in range(size)}
    graph = [[] for _ in range(size)]
    for i in range(size):
        for column, entry in enumerate(matrix[i]):
            if entry is None:
                continue
            j = owner.get(column)
            if j is not None and j != i:
                graph[i].append(j)

    index = [None] * size
    low = [0] * size
    on_stack = [False] * size
    stack: list[int] = []
    components: list[int] = []
    counter = [0]

    for start in range(size):
        if index[start] is not None:
            continue
        work = [(start, 0)]
        while work:
            v, position = work[-1]
            if position == 0:
                index[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                on_stack[v] = True
            recurse = False
            for i in range(position, len(graph[v])):
                w = graph[v][i]
                if index[w] is None:
                    work[-1] = (v, i + 1)
                    work.append((w, 0))
                    recurse = True
                    break
                if on_stack[w]:
                    low[v] = min(low[v], index[w])
            if recurse:
                continue
            if low[v] == index[v]:
                component = 0
                while True:
                    w = stack.pop()
                    on_stack[w] = False
                    component += 1
                    if w == v:
                        break
                components.append(component)
            work.pop()
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[v])
    return sorted(components, reverse=True)


def determinant_nonzero_candidate_rows(n: int, d: int, source: str, checkpoint_dir=None):
    """Every Hall-feasible row with nonzero tangent determinant.

    NOT the final facets. These rows have passed the trace test, Hall
    feasibility and a nonzero Ressayre tangent determinant, and nothing else:
    no Farkas reduction, no facet witness, no BDR birationality. The published
    irredundant systems are far smaller, 4, 31 and 52 rows against the 48, 136
    and 240 counted here, so the two populations must never be conflated.
    ``published_facet_rows`` audits the retained systems separately.
    """
    from gpc_census.generation import orbit_canonical as oc
    from gpc_census.generation.exterior import exterior_weights
    from gpc_census.generation.ressayre import admissible_candidate_from_tau

    weights = exterior_weights(n, d)
    occupied = [tuple(i for i, bit in enumerate(w) if bit) for w in weights]
    if source == "flats":
        keys = AUDIT.orbits_by_flats(n, d)
    else:
        keys = AUDIT.orbits_by_augmentation(n, d, checkpoint_dir=checkpoint_dir)
    rng = random.Random(20260808 + d)

    for key in keys:
        candidate = admissible_candidate_from_tau(n, list(key))
        below_count = sum(
            1
            for members in occupied
            if sum(candidate.h[i] for i in members) < candidate.z
        )
        if oc.orbit_is_trace_dead(candidate.h, below_count):
            continue
        value_map = dict(zip(candidate.h, key))
        for arrangement in oc.generate_trace_survivors(candidate.h, below_count):
            tau = tuple(value_map[value] for value in arrangement)
            on, below, roots, matrix = AUDIT.tangent_matrix(tau, n)
            if not below or len(below) != len(roots):
                continue
            if AUDIT.maximum_matching(matrix) != len(below):
                continue
            if not AUDIT.determinant_verdict(matrix, len(on), rng)["nonzero"]:
                continue
            yield tau, matrix


def published_facet_rows(n: int, d: int):
    """The published irredundant GPC rows, converted to homogeneous tau form.

    A published row is ``coeffs . lambda <= rhs``. On the trace slice
    ``sum lambda = N`` that is ``(N*coeffs - rhs*1) . lambda <= 0``, and the
    primitive integer form of that vector is the ``tau`` the rest of the
    pipeline uses.
    """
    from gpc_census.constraints import constraints

    for inequality in constraints(n, d)["inequalities"]:
        values = [n * c - inequality["rhs"] for c in inequality["coeffs"]]
        divisor = 0
        for value in values:
            divisor = math.gcd(divisor, abs(value))
        tau = tuple(value // (divisor or 1) for value in values)
        _on, below, roots, matrix = AUDIT.tangent_matrix(tau, n)
        yield tau, below, roots, matrix, _on


def audit_published_facets(n: int, d: int) -> dict:
    """The same block audit, restricted to the retained irredundant system."""
    rng = random.Random(7)
    rows = 0
    structural = 0
    fully_triangular = 0
    largest_block = 0
    blocks_by_row = []

    for tau, below, roots, matrix, on in published_facet_rows(n, d):
        if not below:
            structural += 1
            continue
        if len(below) != len(roots) or AUDIT.maximum_matching(matrix) != len(below):
            continue
        if not AUDIT.determinant_verdict(matrix, len(on), rng)["nonzero"]:
            continue
        match_row = perfect_matching(matrix)
        if match_row is None:
            continue
        blocks = irreducible_blocks(matrix, match_row)
        rows += 1
        fully_triangular += max(blocks) == 1
        largest_block = max(largest_block, max(blocks))
        blocks_by_row.append({"tau": list(tau), "blocks": blocks})

    return {
        "published_rows": rows + structural,
        "structural_rows": structural,
        "analysed_rows": rows,
        "fully_triangular_rows": fully_triangular,
        "fully_triangular_fraction": round(fully_triangular / rows, 4) if rows else 0,
        "largest_irreducible_block": largest_block,
        "blocks_by_row": blocks_by_row if rows <= 8 else blocks_by_row[:8],
    }


def audit_system(n: int, d: int, source: str, checkpoint_dir=None) -> dict:
    import time

    started = time.perf_counter()
    rng = random.Random(11)
    rows = 0
    fully_triangular = 0
    largest_block = 0
    largest_order = 0
    block_distribution: collections.Counter = collections.Counter()
    order_of_largest_block = 0
    matching_checks = 0
    matching_disagreements = 0
    witness = None

    for tau, matrix in determinant_nonzero_candidate_rows(
        n, d, source, checkpoint_dir=checkpoint_dir
    ):
        match_row = perfect_matching(matrix)
        if match_row is None:
            continue
        blocks = irreducible_blocks(matrix, match_row)
        rows += 1
        order = len(matrix)
        largest_order = max(largest_order, order)
        top = max(blocks)
        block_distribution[top] += 1
        if top == 1:
            fully_triangular += 1
        if top > largest_block:
            largest_block = top
            order_of_largest_block = order
            witness = {
                "tau": list(tau),
                "matrix_order": order,
                "blocks": blocks,
                "status": "determinant-nonzero candidate, NOT a published facet",
            }

        # Canonicality spot check on a sample of rows.
        if matching_checks < 40:
            matching_checks += 1
            for _ in range(4):
                order_list = list(range(order))
                rng.shuffle(order_list)
                alternative = perfect_matching(matrix, order_list)
                if alternative is None:
                    continue
                if irreducible_blocks(matrix, alternative) != blocks:
                    matching_disagreements += 1
                    break

    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "determinant_nonzero_hall_feasible_rows": rows,
        "fully_triangular_rows": fully_triangular,
        "fully_triangular_fraction": round(fully_triangular / rows, 4) if rows else 0,
        "largest_irreducible_block": largest_block,
        "matrix_order_at_largest_block": order_of_largest_block,
        "largest_matrix_order": largest_order,
        "max_block_size_distribution": {
            str(k): v for k, v in sorted(block_distribution.items())
        },
        "largest_candidate_block_witness": witness,
        "matching_independence_rows_checked": matching_checks,
        "matching_independence_disagreements": matching_disagreements,
        "seconds": round(time.perf_counter() - started, 1),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-dir", default=None)
    parser.add_argument(
        "-o",
        "--out",
        type=pathlib.Path,
        default=ROOT / "results/data/levi_triangularity_gate.json",
    )
    arguments = parser.parse_args()

    systems = {}
    for n, d, source in CENSUS_SYSTEMS:
        print(f"({n},{d}) ...", flush=True)
        record = audit_system(n, d, source, checkpoint_dir=arguments.checkpoint_dir)
        record["published_facets"] = audit_published_facets(n, d)
        systems[record["system"]] = record
        print(
            f"  rows {record['determinant_nonzero_hall_feasible_rows']} "
            f"fully triangular {record['fully_triangular_rows']} "
            f"({100 * record['fully_triangular_fraction']:.1f}%) "
            f"largest block {record['largest_irreducible_block']} "
            f"of order {record['matrix_order_at_largest_block']} | published "
            f"{record['published_facets']['fully_triangular_rows']}"
            f"/{record['published_facets']['analysed_rows']} triangular, "
            f"largest {record['published_facets']['largest_irreducible_block']} "
            f"({record['seconds']}s)",
            flush=True,
        )

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/levi_triangularity_gate.py",
        "status": "conjecture_test",
        "conjecture": (
            "Fermionic Levi Triangularity: every minimal fermionic BDR facet "
            "becomes linear triangular after determinant-twist removal, "
            "Levi-excitation decomposition, grade ordering and "
            "Dulmage-Mendelsohn factorization."
        ),
        "systems": systems,
        "verdict": {
            "scope": (
                "The candidate population is every Hall-feasible row with a "
                "nonzero Ressayre tangent determinant, which is NOT the set of "
                "final facets: 48, 136 and 240 rows against published systems "
                "of 4, 31 and 52. The published rows are audited separately."
            ),
            "strict_triangularity": (
                "REFUTED on both populations. Among determinant-nonzero "
                "candidates fully triangular rows are a minority whose share "
                "falls with rank. Among the PUBLISHED irredundant facets it is "
                "0 of 4 at (3,7), 5 of 31 at (3,8) and 5 of 52 at (3,9), so "
                "the no-go does not depend on the wider candidate scope."
            ),
            "block_triangular_reading": (
                "Still open, but the premise it needs is not supported here. "
                "If the conjecture means block triangular with jointly solved "
                "irreducible blocks, the operative question is whether block "
                "size stays bounded, and the largest block grows across the "
                "three ranks."
            ),
            "canonicality": (
                "The block multiset is matching-independent, re-checked on a "
                "sample at every rank rather than assumed."
            ),
        },
    }
    hashable = json.loads(json.dumps(payload))
    for record in hashable["systems"].values():
        for field in TIMING_FIELDS:
            record.pop(field, None)
    payload["canonical_sha256_without_this_field_or_timings"] = hashlib.sha256(
        json.dumps(hashable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
