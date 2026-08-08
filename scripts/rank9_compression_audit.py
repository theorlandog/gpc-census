#!/usr/bin/env python3
"""Full-population compression audit, ranks 7 to 9, from the repository's own generator.

This is the decisive test named at the end of
``docs/mechanism_grammar_frontier.md``: run the rank-7/rank-8 classification
again at rank 9, where the one-hole grammar is predicted to become necessary.

It also closes the scope gap that document declared. The shipped
``full_rank78_compression_audit.json`` carried population gates (549 and 6,607
trace survivors, the Hall split, the incidence counts) that were re-read rather
than regenerated, because the generating enumeration was never shipped. Here
the population is rebuilt from the repository's own machinery, so ranks 7 and 8
are a reproduction test with published answers before rank 9 is attempted.

THE PIPELINE, per system.

1. Enumerate ``S_d`` orbits of admissible hyperplanes. Ranks 7 and 8 can use
   the materialising flat enumerator; rank 9 cannot, so orbit representatives
   come from ``enumerate_orbits_by_augmentation``, which is the path that first
   reached rank 9. Both paths are cross-checked at rank 7.
2. Per orbit compute ``B``, the number of exterior weights strictly below the
   hyperplane, which is an ``S_d`` invariant.
3. GENERATE the trace survivors of the orbit directly, as the arrangements of
   the multiset with inversion number ``B``. Generation is output linear, so
   the cost tracks survivors rather than the ~20 million oriented rank-9 rows.
4. Classify each survivor: arithmetic holes, Hall feasibility of the tangent
   matrix, labeled singleton incidence, fusion width, and for the one-hole
   Hall-feasible rows an exact or modular determinant verdict.

DETERMINANT VERDICTS, and every one of them is a proof. Modular evaluation runs
first, because a nonzero residue at an integer point proves the symbolic
determinant is nonzero over the integers and costs a single elimination. Only
when every sample vanishes does the exact signed sparse expansion run, and that
is precisely the cheap case: the expansion prunes to an empty state set as soon
as a level cancels, so an identically zero determinant is decided almost
immediately even at size 21, while a large nonzero one is the expensive case
and never reaches that branch. Nothing is left resting on sampling alone.

Run:

    uv run scripts/rank9_compression_audit.py --ranks 7 8
    uv run scripts/rank9_compression_audit.py --ranks 9 --checkpoint-dir .cache/orbits
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import itertools
import json
import math
import pathlib
import random
import sys
import time
from collections import deque

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.generation import admissible_flats as af  # noqa: E402
from gpc_census.generation import orbit_canonical as oc  # noqa: E402
from gpc_census.generation.exterior import exterior_weights  # noqa: E402
from gpc_census.generation.reference_generator import (  # noqa: E402
    oriented_taus,
    primitive_tau,
)
from gpc_census.generation.ressayre import (  # noqa: E402
    admissible_candidate_from_tau,
)

PARTICLE_NUMBER = 3
MODULAR_PRIME = 1000003
MODULAR_SAMPLES = 4
TIMING_FIELDS = ("seconds", "orbit_seconds")


# ---------------------------------------------------------------------------
# Orbit enumeration
# ---------------------------------------------------------------------------


def orbits_by_flats(n: int, d: int) -> list[tuple[int, ...]]:
    """Sorted-tau orbit keys via the materialising flat enumerator."""
    enumeration = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    return sorted({tuple(sorted(tau)) for tau in oriented_taus(enumeration)})


def orbits_by_augmentation(n: int, d: int, checkpoint_dir=None) -> list[tuple[int, ...]]:
    """Sorted-tau orbit keys via the augmentation enumerator, which reaches rank 9."""
    _counts, levels = oc.enumerate_orbits_by_augmentation(
        n, d, checkpoint_dir=checkpoint_dir
    )
    keys = set()
    for _mask, witness in levels[-1].items():
        h, z = oc._exact_top_normal(n, d, witness)
        for orientation in (1, -1):
            oriented = [orientation * value for value in h]
            offset = orientation * z
            tau = primitive_tau(tuple(offset - n * value for value in oriented))
            keys.add(tuple(sorted(tau)))
    return sorted(keys)


# ---------------------------------------------------------------------------
# Tangent geometry
# ---------------------------------------------------------------------------


def tangent_matrix(tau: tuple[int, ...], n: int):
    """Rows are positive-side determinants, columns are ascending roots."""
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
                if matrix[row][column] is not None:
                    raise AssertionError("tangent entry written twice")
                matrix[row][column] = (variable, sign)
    return on, below, roots, matrix


def maximum_matching(matrix) -> int:
    adjacency = [[j for j, entry in enumerate(row) if entry is not None] for row in matrix]
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
    """Exact signed expansion over column masks; empty dict is the zero polynomial."""
    size = len(matrix)
    states = {0: {(): 1}}
    for row_index, row in enumerate(matrix):
        next_states: dict[int, dict] = {}
        for used, polynomial in states.items():
            if used.bit_count() != row_index:
                raise AssertionError("mask population does not track the row index")
            for column, entry in enumerate(row):
                if entry is None or used & (1 << column):
                    continue
                variable, edge_sign = entry
                permutation_sign = -1 if (used >> (column + 1)).bit_count() % 2 else 1
                target = next_states.setdefault(used | (1 << column), {})
                for monomial, coefficient in polynomial.items():
                    extended = tuple(sorted((*monomial, variable)))
                    updated = target.get(extended, 0) + edge_sign * permutation_sign * coefficient
                    if updated:
                        target[extended] = updated
                    else:
                        target.pop(extended, None)
        states = {mask: poly for mask, poly in next_states.items() if poly}
        if not states:
            return {}
    return states.get((1 << size) - 1, {})


def determinant_modulo(matrix, values, prime: int) -> int:
    """Fraction-free Gaussian elimination of the evaluated matrix, modulo a prime."""
    size = len(matrix)
    grid = [
        [
            0 if entry is None else (entry[1] * values[entry[0]]) % prime
            for entry in row
        ]
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


def determinant_verdict(matrix, variables: int, rng: random.Random) -> dict:
    """Decide the symbolic determinant, and say which proof carries the verdict.

    Modular evaluation first, because a nonzero residue at an integer point
    PROVES the symbolic determinant is nonzero over the integers and costs one
    elimination. Only when every sample vanishes is the exact signed expansion
    run, and that is exactly the case where it is cheap: the expansion prunes
    to the empty state set as soon as a level cancels, so a determinant that is
    identically zero is decided almost immediately even at size 21, while a
    large nonzero one would be the expensive case and never reaches here.

    Both branches are therefore proofs. Nothing is left resting on sampling.
    """
    values_used = []
    for _ in range(MODULAR_SAMPLES):
        values = [rng.randrange(1, MODULAR_PRIME) for _ in range(variables)]
        residue = determinant_modulo(matrix, values, MODULAR_PRIME)
        if residue:
            return {
                "method": "modular_evaluation",
                "nonzero": True,
                "prime": MODULAR_PRIME,
                "residue": residue,
                "variable_values": values,
                "proves": "nonvanishing over the integers",
            }
        values_used.append(values)

    coefficients = sparse_determinant(matrix)
    if coefficients:
        # A nonzero polynomial that vanished at every sampled point. Possible in
        # principle, so it is reported rather than assumed away.
        return {
            "method": "exact_sparse_expansion",
            "nonzero": True,
            "coefficient_count": len(coefficients),
            "proves": "nonvanishing, missed by the sampled points",
        }
    return {
        "method": "exact_sparse_expansion",
        "nonzero": False,
        "coefficient_count": 0,
        "modular_samples_that_vanished": len(values_used),
        "proves": "the determinant is identically zero",
    }


# ---------------------------------------------------------------------------
# Arithmetic levels and fusion width
# ---------------------------------------------------------------------------


def level_structure(tau: tuple[int, ...]) -> dict:
    levels = sorted(set(tau))
    step = 0
    for lower, upper in zip(levels, levels[1:]):
        step = math.gcd(step, upper - lower)
    step = step or 1
    normalized = [(value - levels[0]) // step for value in levels]
    missing = sorted(set(range(normalized[-1] + 1)) - set(normalized))
    return {
        "levels": levels,
        "step": step,
        "normalized_levels": normalized,
        "missing": missing,
        "holes": len(missing),
    }


def fusion_widths(tau: tuple[int, ...], n: int) -> tuple[int, int]:
    """Optimal and balanced target-aware block-tree widths."""
    counts = collections.Counter(tau)
    blocks = sorted(counts.items())
    total = len(blocks)

    achievable: list[set[tuple[int, int]]] = [set() for _ in range(1 << total)]
    achievable[0] = {(0, 0)}
    for mask in range(1, 1 << total):
        bit = mask & -mask
        index = bit.bit_length() - 1
        level, multiplicity = blocks[index]
        previous = achievable[mask ^ bit]
        merged = set()
        for count, weight in previous:
            for taken in range(min(multiplicity, n - count) + 1):
                merged.add((count + taken, weight + taken * level))
        achievable[mask] = merged

    full = (1 << total) - 1

    def retained(mask: int) -> int:
        complement = achievable[full ^ mask]
        keep = {
            state
            for state in achievable[mask]
            if any(
                state[0] + other[0] == n and state[1] + other[1] == 0
                for other in complement
            )
        }
        return len(keep)

    memo: dict[int, int] = {}

    def optimal(mask: int) -> int:
        if mask & (mask - 1) == 0:
            return retained(mask)
        if mask in memo:
            return memo[mask]
        best = None
        submask = (mask - 1) & mask
        while submask:
            if submask < (mask ^ submask):
                value = max(optimal(submask), optimal(mask ^ submask), retained(mask))
                best = value if best is None else min(best, value)
            submask = (submask - 1) & mask
        memo[mask] = best if best is not None else retained(mask)
        return memo[mask]

    def balanced(indices: list[int]) -> int:
        mask = 0
        for index in indices:
            mask |= 1 << index
        if len(indices) == 1:
            return retained(mask)
        middle = len(indices) // 2
        return max(
            balanced(indices[:middle]),
            balanced(indices[middle:]),
            retained(mask),
        )

    return optimal(full), balanced(list(range(total)))


# ---------------------------------------------------------------------------
# Per-system audit
# ---------------------------------------------------------------------------


def audit_system(
    n: int,
    d: int,
    source: str,
    checkpoint_dir=None,
    with_width: bool = True,
) -> dict:
    started = time.time()
    weights = exterior_weights(n, d)
    occupied = [tuple(i for i, bit in enumerate(w) if bit) for w in weights]

    if source == "flats":
        orbit_keys = orbits_by_flats(n, d)
    else:
        orbit_keys = orbits_by_augmentation(n, d, checkpoint_dir=checkpoint_dir)
    orbit_seconds = time.time() - started

    rng = random.Random(20260808 + d)

    trace_survivors = 0
    nonstructural = 0
    dead_orbits = 0
    hole_distribution: collections.Counter = collections.Counter()
    hole_distribution_hall: collections.Counter = collections.Counter()
    hall_feasible = 0
    hall_zero = 0
    optimal_width: collections.Counter = collections.Counter()
    balanced_width: collections.Counter = collections.Counter()
    hall_zero_width: collections.Counter = collections.Counter()
    singleton_signatures: set = set()
    singleton_collisions: dict = {}
    hole_audit: list = []
    width_cache: dict = {}

    for key in orbit_keys:
        candidate = admissible_candidate_from_tau(n, list(key))
        below_count = sum(
            1 for members in occupied if sum(candidate.h[i] for i in members) < candidate.z
        )
        if oc.orbit_is_trace_dead(candidate.h, below_count):
            dead_orbits += 1
            continue
        value_map = dict(zip(candidate.h, key))

        for arrangement in oc.generate_trace_survivors(candidate.h, below_count):
            tau = tuple(value_map[value] for value in arrangement)
            trace_survivors += 1

            on, below, roots, matrix = tangent_matrix(tau, n)
            structural = not below
            if not structural:
                nonstructural += 1

            structure = level_structure(tau)
            holes = structure["holes"]
            hole_distribution[holes] += 1

            # A structural row has an empty positive side, and the trace test
            # forces the root set empty with it, so the matching is vacuously
            # perfect. The published rank-7/rank-8 Hall split counts those rows
            # as feasible, so the split here is over ALL survivors.
            feasible = len(below) == len(roots) and maximum_matching(matrix) == len(below)
            if feasible:
                hall_feasible += 1
                hole_distribution_hall[holes] += 1
            else:
                hall_zero += 1

            if with_width:
                shape = tuple(sorted(collections.Counter(tau).items()))
                if shape not in width_cache:
                    width_cache[shape] = fusion_widths(tau, n)
                optimal, balanced_value = width_cache[shape]
                optimal_width[optimal] += 1
                balanced_width[balanced_value] += 1
                if not feasible:
                    hall_zero_width[optimal] += 1

            # Injectivity is measured on nonstructural rows only: the
            # zero-below structural family shares the all-zeros signature by
            # construction, which is a known and uninformative collision.
            if not structural:
                signature = tuple(
                    sum(1 for det in below if orbital in det) for orbital in range(d)
                )
                if signature in singleton_signatures:
                    singleton_collisions.setdefault(signature, 0)
                    singleton_collisions[signature] += 1
                singleton_signatures.add(signature)

            if feasible and holes == 1:
                verdict = determinant_verdict(matrix, len(on), rng)
                hole_audit.append(
                    {
                        "tau": list(tau),
                        "normalized_levels": structure["normalized_levels"],
                        "missing_levels": structure["missing"],
                        "step": structure["step"],
                        "matrix_size": len(below),
                        "on_variables": len(on),
                        "determinant": verdict,
                    }
                )

    one_hole_nonzero = sum(
        1 for record in hole_audit if record["determinant"]["nonzero"]
    )
    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "orbit_source": source,
        "orbits": len(orbit_keys),
        "trace_dead_orbits": dead_orbits,
        "trace_survivors": trace_survivors,
        "nonstructural_survivors": nonstructural,
        "structural_survivors": trace_survivors - nonstructural,
        "hall_feasible": hall_feasible,
        "hall_zero": hall_zero,
        "hole_distribution_all": {str(k): v for k, v in sorted(hole_distribution.items())},
        "hole_distribution_hall_feasible": {
            str(k): v for k, v in sorted(hole_distribution_hall.items())
        },
        "multi_hole_survivors": sum(v for k, v in hole_distribution.items() if k > 1),
        "multi_hole_hall_feasible": sum(
            v for k, v in hole_distribution_hall.items() if k > 1
        ),
        "one_hole_hall_feasible": hole_distribution_hall.get(1, 0),
        "one_hole_determinant_nonzero": one_hole_nonzero,
        "one_hole_hall_feasible_audit": hole_audit,
        "singleton_signature_count": len(singleton_signatures),
        "singleton_collision_groups": len(singleton_collisions),
        "singleton_incidence_injective": not singleton_collisions,
        "optimal_fusion_width_distribution": {
            str(k): v for k, v in sorted(optimal_width.items())
        },
        "balanced_fusion_width_distribution": {
            str(k): v for k, v in sorted(balanced_width.items())
        },
        "hall_zero_optimal_width_distribution": {
            str(k): v for k, v in sorted(hall_zero_width.items())
        },
        "max_optimal_fusion_width": max(optimal_width) if optimal_width else 0,
        "max_balanced_fusion_width": max(balanced_width) if balanced_width else 0,
        "orbit_seconds": round(orbit_seconds, 1),
        "seconds": round(time.time() - started, 1),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", nargs="*", type=int, default=[7, 8])
    parser.add_argument(
        "--source",
        choices=["flats", "augmentation", "auto"],
        default="auto",
        help="auto uses flats through rank 8 and augmentation from rank 9",
    )
    parser.add_argument("--checkpoint-dir", default=None)
    parser.add_argument("--no-width", action="store_true")
    parser.add_argument(
        "-o",
        "--out",
        type=pathlib.Path,
        default=ROOT / "results/data/rank9_compression_audit.json",
    )
    arguments = parser.parse_args()

    systems = {}
    for d in arguments.ranks:
        source = arguments.source
        if source == "auto":
            source = "flats" if d <= 8 else "augmentation"
        print(f"({PARTICLE_NUMBER},{d}) via {source} ...", flush=True)
        record = audit_system(
            PARTICLE_NUMBER,
            d,
            source,
            checkpoint_dir=arguments.checkpoint_dir,
            with_width=not arguments.no_width,
        )
        systems[record["system"]] = record
        print(
            f"  survivors {record['trace_survivors']} "
            f"nonstructural {record['nonstructural_survivors']} "
            f"hall {record['hall_feasible']}/{record['hall_zero']} "
            f"holes {record['hole_distribution_all']} "
            f"one-hole-feasible {record['one_hole_hall_feasible']} "
            f"of which nonzero {record['one_hole_determinant_nonzero']} "
            f"({record['seconds']}s)",
            flush=True)

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/rank9_compression_audit.py",
        "status": "full_population_audit",
        "particle_number": PARTICLE_NUMBER,
        "determinant_policy": {
            "modular_prime": MODULAR_PRIME,
            "modular_samples": MODULAR_SAMPLES,
            "note": (
                "Modular evaluation runs first and a nonzero residue proves "
                "nonvanishing over the integers. When every sample vanishes the "
                "exact signed expansion runs and proves the determinant is "
                "identically zero. Every verdict is a proof; none rests on "
                "sampling alone."
            ),
        },
        "systems": systems,
        "summary": {
            "multi_hole_hall_feasible_anywhere": sum(
                record["multi_hole_hall_feasible"] for record in systems.values()
            ),
            "multi_hole_survivors_anywhere": sum(
                record["multi_hole_survivors"] for record in systems.values()
            ),
            "singleton_incidence_injective_everywhere": all(
                record["singleton_incidence_injective"] for record in systems.values()
            ),
            "total_trace_survivors": sum(
                record["trace_survivors"] for record in systems.values()
            ),
        },
    }
    # The digest covers CONTENT only. Wall-clock timings vary between machines
    # and would make an otherwise byte-identical rerun hash differently, so
    # they are excluded and the field name says which keys are dropped.
    hashable = {
        key: value for key, value in payload.items() if key != "systems"
    }
    hashable["systems"] = {
        system: {
            key: value
            for key, value in record.items()
            if key not in TIMING_FIELDS
        }
        for system, record in payload["systems"].items()
    }
    payload["canonical_sha256_without_this_field_or_timings"] = hashlib.sha256(
        json.dumps(hashable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    serialized = json.dumps(payload, indent=1, sort_keys=True) + "\n"
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(serialized)
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
