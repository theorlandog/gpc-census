#!/usr/bin/env python3
"""Early go/no-go trials for compressed GPC construction.

The script performs three exact, standard-library-only experiments:

1. Projection/fusion-tree width over all 83 known system-specific Levi
   mechanisms in the supplied pilot and held-out artifacts.
2. A bounded exhaustive N=3 candidate-normal search for d=7,8,9 with primitive
   entries in [-5,5], followed by explicit trace-placement generation, Hall
   matching, and deterministic modular nonvanishing certificates.
3. Singleton-incidence tomography of known mechanisms and of every bounded
   trace placement.

The bounded candidate search is an early trial, not a completeness statement
about the repository's full candidate universe.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict, deque
from fractions import Fraction
from pathlib import Path


def canonical_hash(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def parse_system(system: str) -> tuple[int, int]:
    return tuple(int(part) for part in system.strip("()").split(","))


def load_mechanisms(pilot_path: Path, heldout_path: Path) -> list[dict]:
    pilot = json.loads(pilot_path.read_text())
    heldout = json.loads(heldout_path.read_text())
    mechanisms: list[dict] = []

    for system, payload in pilot["systems"].items():
        seen: dict[tuple[int, ...], dict] = {}
        for record in payload["records"]:
            seen.setdefault(tuple(sorted(record["tau"])), record)
        for shape, record in seen.items():
            mechanisms.append(
                {
                    "system": system,
                    "particle_number": int(payload["particle_number"]),
                    "rank": int(payload["rank"]),
                    "tau": list(shape),
                    "levels": [
                        [int(entry["value"]), int(entry["multiplicity"])]
                        for entry in record["levels"]
                    ],
                    "source": "pilot",
                }
            )

    for system, payload in heldout["systems"].items():
        for record in payload["signature_records"]:
            tau = record.get(
                "representative_tau",
                record.get("tau", record.get("shape")),
            )
            mechanisms.append(
                {
                    "system": system,
                    "particle_number": int(payload["particle_number"]),
                    "rank": int(payload["rank"]),
                    "tau": sorted(int(value) for value in tau),
                    "levels": [
                        [int(value), int(multiplicity)]
                        for value, multiplicity in record["levels"]
                    ],
                    "source": "heldout",
                }
            )

    if len(mechanisms) != 83:
        raise AssertionError(f"expected 83 mechanisms, found {len(mechanisms)}")
    return mechanisms


# ---------------------------------------------------------------------------
# Projection / fusion-tree width
# ---------------------------------------------------------------------------


def achievable_block_states(
    level: int,
    multiplicity: int,
    particle_number: int,
) -> set[tuple[int, int]]:
    return {
        (count, count * level)
        for count in range(min(multiplicity, particle_number) + 1)
    }


def subset_state_sets(
    levels: tuple[tuple[int, int], ...],
    particle_number: int,
) -> tuple[list[set[tuple[int, int]]], list[set[tuple[int, int]]]]:
    """All achievable and target-extendable states for every block subset."""

    block_count = len(levels)
    achievable: list[set[tuple[int, int]]] = [
        set() for _ in range(1 << block_count)
    ]
    achievable[0] = {(0, 0)}
    for mask in range(1, 1 << block_count):
        bit = mask & -mask
        index = bit.bit_length() - 1
        previous = mask ^ bit
        local = achievable_block_states(
            levels[index][0],
            levels[index][1],
            particle_number,
        )
        states: set[tuple[int, int]] = set()
        for left_count, left_sum in achievable[previous]:
            for right_count, right_sum in local:
                if left_count + right_count <= particle_number:
                    states.add(
                        (
                            left_count + right_count,
                            left_sum + right_sum,
                        )
                    )
        achievable[mask] = states

    full = (1 << block_count) - 1
    active: list[set[tuple[int, int]]] = []
    for mask in range(1 << block_count):
        complement = full ^ mask
        active.append(
            {
                (count, total)
                for count, total in achievable[mask]
                if (
                    particle_number - count,
                    -total,
                )
                in achievable[complement]
            }
        )
    return achievable, active


def optimal_tree_width(state_sets: list[set[tuple[int, int]]]) -> int:
    """Exact minimum width over all full binary trees on the block leaves."""

    full = len(state_sets) - 1
    dynamic = [0] * (full + 1)
    for mask in range(1, full + 1):
        local_width = len(state_sets[mask])
        if mask & (mask - 1) == 0:
            dynamic[mask] = local_width
            continue
        anchor = mask & -mask
        best = 10**9
        subset = (mask - 1) & mask
        while subset:
            if subset & anchor:
                other = mask ^ subset
                if other:
                    best = min(
                        best,
                        max(
                            local_width,
                            dynamic[subset],
                            dynamic[other],
                        ),
                    )
            subset = (subset - 1) & mask
        dynamic[mask] = best
    return dynamic[full]


def balanced_tree_width(state_sets: list[set[tuple[int, int]]]) -> int:
    """Exact minimum subject to recursively balanced block splits."""

    full = len(state_sets) - 1
    population = [0] * (full + 1)
    for mask in range(1, full + 1):
        population[mask] = population[mask >> 1] + (mask & 1)

    dynamic = [10**9] * (full + 1)
    for mask in range(1, full + 1):
        local_width = len(state_sets[mask])
        size = population[mask]
        if size == 1:
            dynamic[mask] = local_width
            continue
        left_size = size // 2
        right_size = size - left_size
        anchor = mask & -mask
        subset = (mask - 1) & mask
        while subset:
            other = mask ^ subset
            if (
                subset & anchor
                and population[subset] in {left_size, right_size}
                and population[other] in {left_size, right_size}
            ):
                dynamic[mask] = min(
                    dynamic[mask],
                    max(
                        local_width,
                        dynamic[subset],
                        dynamic[other],
                    ),
                )
            subset = (subset - 1) & mask
    return dynamic[full]


def optimal_chain_width(state_sets: list[set[tuple[int, int]]]) -> int:
    """Exact minimum over one-block-at-a-time projection chains."""

    full = len(state_sets) - 1
    dynamic = [0] * (full + 1)
    for mask in range(1, full + 1):
        local_width = len(state_sets[mask])
        if mask & (mask - 1) == 0:
            dynamic[mask] = local_width
            continue
        best = 10**9
        bits = mask
        while bits:
            bit = bits & -bits
            bits ^= bit
            best = min(
                best,
                max(local_width, dynamic[mask ^ bit]),
            )
        dynamic[mask] = best
    return dynamic[full]


def sorted_chain_width(state_sets: list[set[tuple[int, int]]], blocks: int) -> int:
    mask = 0
    width = 0
    for index in range(blocks):
        mask |= 1 << index
        width = max(width, len(state_sets[mask]))
    return width


def projection_width_audit(mechanisms: list[dict]) -> dict:
    records = []
    for mechanism in mechanisms:
        levels = tuple(
            (int(value), int(multiplicity))
            for value, multiplicity in mechanism["levels"]
        )
        achievable, active = subset_state_sets(
            levels,
            int(mechanism["particle_number"]),
        )
        record = {
            **mechanism,
            "zero_target_optimal_tree_width": optimal_tree_width(active),
            "zero_target_balanced_tree_width": balanced_tree_width(active),
            "zero_target_optimal_chain_width": optimal_chain_width(active),
            "zero_target_sorted_chain_width": sorted_chain_width(
                active,
                len(levels),
            ),
            "unpruned_full_state_count": len(achievable[-1]),
        }
        records.append(record)

    def distribution(field: str) -> dict[str, int]:
        return {
            str(key): value
            for key, value in sorted(
                Counter(record[field] for record in records).items()
            )
        }

    return {
        "definition": (
            "For a block subset S, retain exactly the pairs (k,s) achievable "
            "inside S that can be completed in the complement to the global "
            "target (N,0). Tree width is the maximum number of such states at "
            "a tree node."
        ),
        "mechanism_count": len(records),
        "optimal_tree_width_distribution": distribution(
            "zero_target_optimal_tree_width"
        ),
        "balanced_tree_width_distribution": distribution(
            "zero_target_balanced_tree_width"
        ),
        "optimal_chain_width_distribution": distribution(
            "zero_target_optimal_chain_width"
        ),
        "sorted_chain_width_distribution": distribution(
            "zero_target_sorted_chain_width"
        ),
        "max_optimal_tree_width": max(
            record["zero_target_optimal_tree_width"]
            for record in records
        ),
        "max_balanced_tree_width": max(
            record["zero_target_balanced_tree_width"]
            for record in records
        ),
        "max_optimal_chain_width": max(
            record["zero_target_optimal_chain_width"]
            for record in records
        ),
        "max_sorted_chain_width": max(
            record["zero_target_sorted_chain_width"]
            for record in records
        ),
        "max_unpruned_full_state_count": max(
            record["unpruned_full_state_count"]
            for record in records
        ),
        "records": records,
    }


def achievable_orbital_states(
    indices: frozenset[int],
    tau: tuple[int, ...],
    particle_number: int,
) -> set[tuple[int, int]]:
    states = {(0, 0)}
    for index in sorted(indices):
        value = tau[index]
        states |= {
            (count + 1, total + value)
            for count, total in list(states)
            if count < particle_number
        }
    return states


def balanced_contiguous_tree(indices: tuple[int, ...]):
    if len(indices) == 1:
        return indices[0]
    middle = len(indices) // 2
    return (
        balanced_contiguous_tree(indices[:middle]),
        balanced_contiguous_tree(indices[middle:]),
    )


def symmetric_extremes_tree(rank: int):
    groups: list[object] = [
        (index, rank - 1 - index)
        for index in range(rank // 2)
    ]
    if rank % 2:
        groups.append(rank // 2)

    def combine(values: list[object]):
        if len(values) == 1:
            return values[0]
        middle = len(values) // 2
        return (
            combine(values[:middle]),
            combine(values[middle:]),
        )

    return combine(groups)


def tree_nodes(tree) -> list[frozenset[int]]:
    nodes: list[frozenset[int]] = []

    def visit(node) -> frozenset[int]:
        if isinstance(node, int):
            result = frozenset({node})
            nodes.append(result)
            return result
        left = visit(node[0])
        right = visit(node[1])
        result = left | right
        nodes.append(result)
        return result

    visit(tree)
    return nodes


def explicit_tree_width(
    tau: tuple[int, ...],
    particle_number: int,
    tree,
) -> int:
    full = frozenset(range(len(tau)))
    cache: dict[frozenset[int], set[tuple[int, int]]] = {}

    def states(indices: frozenset[int]) -> set[tuple[int, int]]:
        if indices not in cache:
            cache[indices] = achievable_orbital_states(
                indices,
                tau,
                particle_number,
            )
        return cache[indices]

    width = 0
    for node in tree_nodes(tree):
        complement = full - node
        complement_states = states(complement)
        active = {
            (count, total)
            for count, total in states(node)
            if (
                particle_number - count,
                -total,
            )
            in complement_states
        }
        width = max(width, len(active))
    return width


def synthetic_projection_stress() -> list[dict]:
    systems = (
        (3, 10),
        (3, 20),
        (3, 40),
        (4, 10),
        (4, 20),
        (4, 40),
        (5, 10),
        (5, 20),
        (5, 40),
        (10, 20),
        (20, 40),
    )
    records = []
    for particle_number, rank in systems:
        target_index_sum = particle_number * (rank - 1) // 2
        tau = tuple(
            particle_number * index - target_index_sum
            for index in range(rank)
        )
        records.append(
            {
                "particle_number": particle_number,
                "rank": rank,
                "normal_family": (
                    "consecutive index levels tau_j = N*j - "
                    "floor(N*(d-1)/2)"
                ),
                "balanced_contiguous_width_upper_bound": explicit_tree_width(
                    tau,
                    particle_number,
                    balanced_contiguous_tree(tuple(range(rank))),
                ),
                "symmetric_extremes_width_upper_bound": explicit_tree_width(
                    tau,
                    particle_number,
                    symmetric_extremes_tree(rank),
                ),
            }
        )
    return records


# ---------------------------------------------------------------------------
# Bounded candidate search
# ---------------------------------------------------------------------------


def integer_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    matrix = [
        [Fraction(value) for value in row]
        for row in rows
    ]
    row_count = len(matrix)
    column_count = len(matrix[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (
                index
                for index in range(rank, row_count)
                if matrix[index][column] != 0
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        head = matrix[rank][column]
        matrix[rank] = [
            value / head for value in matrix[rank]
        ]
        for index in range(row_count):
            if index == rank or matrix[index][column] == 0:
                continue
            factor = matrix[index][column]
            matrix[index] = [
                left - factor * right
                for left, right in zip(matrix[index], matrix[rank])
            ]
        rank += 1
    return rank


def affine_rank(points: list[tuple[int, ...]]) -> int:
    if not points:
        return -1
    base = points[0]
    return integer_rank(
        [
            [left - right for left, right in zip(point, base)]
            for point in points[1:]
        ]
    )


def progression_data(tau: tuple[int, ...]) -> dict:
    values = sorted(set(tau))
    step = 0
    for left, right in zip(values, values[1:]):
        step = math.gcd(step, right - left)
    if step == 0:
        step = 1
    normalized = [
        (value - values[0]) // step
        for value in values
    ]
    missing = sorted(
        set(range(normalized[-1] + 1)) - set(normalized)
    )
    return {
        "step": step,
        "normalized_levels": normalized,
        "missing_levels": missing,
        "hole_count": len(missing),
    }


def multiset_maximum_inversions(tau: tuple[int, ...]) -> int:
    counts = Counter(tau)
    values = sorted(counts)
    return sum(
        counts[left] * counts[right]
        for index, left in enumerate(values)
        for right in values[index + 1 :]
    )


def candidate_shapes(
    particle_number: int,
    rank: int,
    lower: int,
    upper: int,
) -> list[dict]:
    determinants = list(
        itertools.combinations(range(rank), particle_number)
    )
    output = []
    for tau in itertools.combinations_with_replacement(
        range(lower, upper + 1),
        rank,
    ):
        if len(set(tau)) < 2:
            continue
        divisor = 0
        for value in tau:
            divisor = math.gcd(divisor, abs(value))
        if divisor != 1:
            continue

        on_points: list[tuple[int, ...]] = []
        positive: list[tuple[int, ...]] = []
        for determinant in determinants:
            total = sum(tau[index] for index in determinant)
            if total == 0:
                point = [0] * rank
                for index in determinant:
                    point[index] = 1
                on_points.append(tuple(point))
            elif total > 0:
                positive.append(determinant)

        if affine_rank(on_points) != rank - 2:
            continue

        maximum_inversions = multiset_maximum_inversions(tau)
        below_count = len(positive)
        if below_count > maximum_inversions:
            continue

        singleton = [0] * rank
        for determinant in positive:
            for index in determinant:
                singleton[index] += 1

        output.append(
            {
                "tau": list(tau),
                "on_weight_count": len(on_points),
                "below_weight_count": below_count,
                "maximum_inversions": maximum_inversions,
                "sorted_singleton_incidence": sorted(singleton),
                **progression_data(tau),
            }
        )
    return output


def inversion_number(sequence: tuple[int, ...]) -> int:
    return sum(
        sequence[left] > sequence[right]
        for left in range(len(sequence))
        for right in range(left + 1, len(sequence))
    )


def generate_multiset_permutations_at_inversion(
    tau: tuple[int, ...],
    target: int,
):
    counts_by_value = Counter(tau)
    values = tuple(sorted(counts_by_value))
    counts = [counts_by_value[value] for value in values]
    total = sum(counts)

    def recurse(
        prefix: list[int],
        remaining_target: int,
        remaining: int,
    ):
        if remaining == 0:
            if remaining_target == 0:
                yield tuple(prefix)
            return

        maximum = 0
        for left in range(len(values)):
            for right in range(left + 1, len(values)):
                maximum += counts[left] * counts[right]
        if remaining_target < 0 or remaining_target > maximum:
            return

        smaller = 0
        for index, value in enumerate(values):
            count = counts[index]
            if count:
                contribution = smaller
                counts[index] -= 1
                prefix.append(value)
                yield from recurse(
                    prefix,
                    remaining_target - contribution,
                    remaining - 1,
                )
                prefix.pop()
                counts[index] += 1
            smaller += count

    yield from recurse([], target, total)


def tangent_support_graph(
    tau: tuple[int, ...],
    particle_number: int,
):
    rank = len(tau)
    determinants = list(
        itertools.combinations(range(rank), particle_number)
    )
    on = [
        determinant
        for determinant in determinants
        if sum(tau[index] for index in determinant) == 0
    ]
    below = [
        determinant
        for determinant in determinants
        if sum(tau[index] for index in determinant) > 0
    ]
    roots = [
        (left, right)
        for left in range(rank)
        for right in range(left + 1, rank)
        if tau[left] < tau[right]
    ]

    on_set = set(on)
    adjacency = [[] for _ in below]
    for row_index, determinant in enumerate(below):
        support = set(determinant)
        for root_index, (left, right) in enumerate(roots):
            if right in support and left not in support:
                predecessor = tuple(
                    sorted((support - {right}) | {left})
                )
                if predecessor in on_set:
                    adjacency[row_index].append(root_index)
    return on, below, roots, adjacency


def maximum_matching(
    adjacency: list[list[int]],
    right_size: int,
) -> int:
    left_size = len(adjacency)
    left_match = [-1] * left_size
    right_match = [-1] * right_size
    distance = [0] * left_size
    infinity = 10**9

    def breadth_first() -> bool:
        queue: deque[int] = deque()
        found = False
        for left in range(left_size):
            if left_match[left] == -1:
                distance[left] = 0
                queue.append(left)
            else:
                distance[left] = infinity
        while queue:
            left = queue.popleft()
            for right in adjacency[left]:
                paired = right_match[right]
                if paired == -1:
                    found = True
                elif distance[paired] == infinity:
                    distance[paired] = distance[left] + 1
                    queue.append(paired)
        return found

    def depth_first(left: int) -> bool:
        for right in adjacency[left]:
            paired = right_match[right]
            if paired == -1 or (
                distance[paired] == distance[left] + 1
                and depth_first(paired)
            ):
                left_match[left] = right
                right_match[right] = left
                return True
        distance[left] = infinity
        return False

    matching = 0
    while breadth_first():
        for left in range(left_size):
            if left_match[left] == -1 and depth_first(left):
                matching += 1
    return matching


def tangent_symbolic_matrix(
    tau: tuple[int, ...],
    particle_number: int,
):
    rank = len(tau)
    determinants = list(
        itertools.combinations(range(rank), particle_number)
    )
    on = [
        determinant
        for determinant in determinants
        if sum(tau[index] for index in determinant) == 0
    ]
    below = [
        determinant
        for determinant in determinants
        if sum(tau[index] for index in determinant) > 0
    ]
    roots = [
        (left, right)
        for left in range(rank)
        for right in range(left + 1, rank)
        if tau[left] < tau[right]
    ]

    below_index = {
        determinant: index
        for index, determinant in enumerate(below)
    }
    entries = [
        [None] * len(roots)
        for _ in below
    ]
    for variable, determinant in enumerate(on):
        support = set(determinant)
        for root_index, (left, right) in enumerate(roots):
            if left not in support or right in support:
                continue
            reduced = list(determinant)
            annihilation_position = reduced.index(left)
            reduced.pop(annihilation_position)
            creation_position = sum(value < right for value in reduced)
            sign = (
                -1
                if (annihilation_position + creation_position) % 2
                else 1
            )
            image = tuple(sorted(reduced + [right]))
            row = below_index.get(image)
            if row is not None:
                entries[row][root_index] = (sign, variable)
    return on, below, roots, entries


def determinant_modulo(
    entries,
    values: list[int],
    prime: int,
) -> int:
    size = len(entries)
    matrix = [
        [
            0
            if entry is None
            else entry[0] * values[entry[1]] % prime
            for entry in row
        ]
        for row in entries
    ]
    determinant = 1
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if matrix[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = (
                matrix[pivot],
                matrix[column],
            )
            determinant = -determinant % prime
        head = matrix[column][column] % prime
        determinant = determinant * head % prime
        inverse = pow(head, prime - 2, prime)
        for row in range(column + 1, size):
            if matrix[row][column] == 0:
                continue
            factor = matrix[row][column] * inverse % prime
            for offset in range(column, size):
                matrix[row][offset] = (
                    matrix[row][offset]
                    - factor * matrix[column][offset]
                ) % prime
    return determinant


def modular_nonzero_certificate(
    tau: tuple[int, ...],
    particle_number: int,
    attempts: int = 32,
) -> dict | None:
    on, below, roots, entries = tangent_symbolic_matrix(
        tau,
        particle_number,
    )
    if len(below) != len(roots):
        return None
    adjacency = [
        [
            column
            for column, entry in enumerate(row)
            if entry is not None
        ]
        for row in entries
    ]
    if maximum_matching(adjacency, len(roots)) < len(roots):
        return None

    for prime in (1000003, 1000033, 1000037):
        for seed in range(1, attempts + 1):
            state = seed * 2654435761 + len(on) * 97
            values = []
            for variable in range(len(on)):
                state = (
                    1103515245 * state
                    + 12345
                    + variable * 7919
                ) & 0x7FFFFFFF
                values.append(1 + state % (prime - 1))
            residue = determinant_modulo(entries, values, prime)
            if residue:
                return {
                    "prime": prime,
                    "seed": seed,
                    "determinant_residue": residue,
                    "variable_values": values,
                    "matrix_size": len(roots),
                    "variables": len(on),
                }
    return None


def positive_singleton_incidence(
    tau: tuple[int, ...],
    particle_number: int,
) -> tuple[int, ...]:
    counts = [0] * len(tau)
    for determinant in itertools.combinations(
        range(len(tau)),
        particle_number,
    ):
        if sum(tau[index] for index in determinant) > 0:
            for index in determinant:
                counts[index] += 1
    return tuple(counts)


def bounded_candidate_audit() -> dict:
    systems = {}
    modular_certificates = []
    total_trace = 0
    total_hall = 0
    total_multihole_trace = 0
    total_multihole_hall = 0
    total_onehole_hall = 0
    total_nonstructural_placements = 0
    tomography_collisions = []

    for rank in (7, 8, 9):
        shapes = candidate_shapes(3, rank, -5, 5)
        trace_by_holes: Counter[int] = Counter()
        hall_by_holes: Counter[int] = Counter()
        hall_shapes_by_holes: Counter[int] = Counter()
        trace_placements = 0
        hall_placements = 0
        shape_records = []
        incidence_groups: dict[tuple[int, ...], list[dict]] = defaultdict(list)

        for shape_record in shapes:
            shape = tuple(shape_record["tau"])
            target_inversions = (
                shape_record["maximum_inversions"]
                - shape_record["below_weight_count"]
            )
            placements = list(
                generate_multiset_permutations_at_inversion(
                    shape,
                    target_inversions,
                )
            )
            feasible = 0
            maximum_match = 0
            for tau in placements:
                on, below, roots, adjacency = tangent_support_graph(tau, 3)
                if not (
                    len(below)
                    == len(roots)
                    == shape_record["below_weight_count"]
                ):
                    raise AssertionError("trace placement count mismatch")
                matching = maximum_matching(adjacency, len(roots))
                maximum_match = max(maximum_match, matching)
                if matching == len(roots):
                    feasible += 1
                    if shape_record["hole_count"] == 1:
                        total_onehole_hall += 1
                        certificate = modular_nonzero_certificate(tau, 3)
                        if certificate is not None:
                            modular_certificates.append(
                                {
                                    "rank": rank,
                                    "shape": list(shape),
                                    "tau": list(tau),
                                    **certificate,
                                }
                            )

                signature = positive_singleton_incidence(tau, 3)
                incidence_groups[signature].append(
                    {
                        "shape": list(shape),
                        "tau": list(tau),
                        "below_weight_count": (
                            shape_record["below_weight_count"]
                        ),
                    }
                )

            trace_placements += len(placements)
            hall_placements += feasible
            trace_by_holes[shape_record["hole_count"]] += len(placements)
            hall_by_holes[shape_record["hole_count"]] += feasible
            if feasible:
                hall_shapes_by_holes[shape_record["hole_count"]] += 1

            shape_records.append(
                {
                    **shape_record,
                    "trace_placements": len(placements),
                    "hall_feasible_placements": feasible,
                    "maximum_matching_size": maximum_match,
                }
            )

        cross_shape = []
        nonstructural_cross_shape = []
        for signature, group in incidence_groups.items():
            distinct_shapes = {
                tuple(record["shape"])
                for record in group
            }
            if len(distinct_shapes) > 1:
                cross_shape.append(
                    {
                        "signature": list(signature),
                        "records": group,
                    }
                )
                if any(
                    record["below_weight_count"] > 0
                    for record in group
                ):
                    nonstructural_cross_shape.append(
                        {
                            "signature": list(signature),
                            "records": group,
                        }
                    )

        total_trace += trace_placements
        total_hall += hall_placements
        total_multihole_trace += sum(
            count
            for holes, count in trace_by_holes.items()
            if holes >= 2
        )
        total_multihole_hall += sum(
            count
            for holes, count in hall_by_holes.items()
            if holes >= 2
        )
        total_nonstructural_placements += sum(
            len(group)
            for group in incidence_groups.values()
            if any(
                record["below_weight_count"] > 0
                for record in group
            )
        )
        tomography_collisions.extend(nonstructural_cross_shape)

        systems[f"(3,{rank})"] = {
            "entry_window": [-5, 5],
            "candidate_shapes": len(shapes),
            "trace_placements": trace_placements,
            "hall_feasible_placements": hall_placements,
            "trace_placements_by_hole_count": {
                str(key): value
                for key, value in sorted(trace_by_holes.items())
            },
            "hall_feasible_by_hole_count": {
                str(key): value
                for key, value in sorted(hall_by_holes.items())
            },
            "hall_feasible_shapes_by_hole_count": {
                str(key): value
                for key, value in sorted(hall_shapes_by_holes.items())
            },
            "singleton_incidence_signatures": len(incidence_groups),
            "singleton_cross_shape_collisions": len(cross_shape),
            "nonstructural_singleton_cross_shape_collisions": len(
                nonstructural_cross_shape
            ),
            "shape_records": shape_records,
        }

    return {
        "scope": (
            "Exhaustive among sorted primitive N=3 normals of ranks 7, 8, "
            "and 9 whose entries lie in [-5,5], whose zero-level exterior "
            "weights have affine rank d-2, and which possess explicit "
            "trace-compatible multiset placements."
        ),
        "systems": systems,
        "total_trace_placements": total_trace,
        "total_hall_feasible_placements": total_hall,
        "multi_hole_trace_placements": total_multihole_trace,
        "multi_hole_hall_feasible_placements": total_multihole_hall,
        "one_hole_hall_feasible_placements": total_onehole_hall,
        "one_hole_modular_nonzero_certificates": modular_certificates,
        "nonstructural_trace_placements": total_nonstructural_placements,
        "nonstructural_singleton_cross_shape_collisions": len(
            tomography_collisions
        ),
    }


def known_mechanism_tomography(mechanisms: list[dict]) -> dict:
    groups: dict[tuple[str, tuple[int, ...]], list[dict]] = defaultdict(list)
    for mechanism in mechanisms:
        tau = tuple(mechanism["tau"])
        signature = tuple(
            sorted(
                positive_singleton_incidence(
                    tau,
                    int(mechanism["particle_number"]),
                )
            )
        )
        groups[(mechanism["system"], signature)].append(mechanism)

    collisions = [
        {
            "system": system,
            "signature": list(signature),
            "shapes": [
                record["tau"]
                for record in group
            ],
        }
        for (system, signature), group in groups.items()
        if len({tuple(record["tau"]) for record in group}) > 1
    ]
    return {
        "mechanisms": len(mechanisms),
        "distinct_system_specific_signatures": len(groups),
        "cross_shape_collisions": collisions,
    }


def validate_tangent_convention(pilot_path: Path) -> dict:
    pilot = json.loads(pilot_path.read_text())
    checked = 0
    certified = 0
    failures = []
    for system, payload in pilot["systems"].items():
        particle_number = int(payload["particle_number"])
        for record in payload["records"]:
            if record["label"] == "STRUCTURAL":
                continue
            checked += 1
            tau = tuple(int(value) for value in record["tau"])
            on, below, roots, adjacency = tangent_support_graph(
                tau,
                particle_number,
            )
            trace_ok = len(below) == len(roots)
            hall_ok = (
                maximum_matching(adjacency, len(roots))
                == len(roots)
            )
            certificate = modular_nonzero_certificate(
                tau,
                particle_number,
                attempts=8,
            )
            if trace_ok and hall_ok and certificate is not None:
                certified += 1
            else:
                failures.append(
                    {
                        "system": system,
                        "index": record["index"],
                        "trace_ok": trace_ok,
                        "hall_ok": hall_ok,
                        "modular_nonzero": certificate is not None,
                    }
                )
    return {
        "nonstructural_pilot_facets_checked": checked,
        "passed_trace_hall_and_modular_nonzero": certified,
        "failures": failures,
    }


DATA = Path(__file__).resolve().parents[1] / "results" / "data"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pilot",
        type=Path,
        default=DATA / "levi_fusion_pilot.json",
    )
    parser.add_argument(
        "--heldout",
        type=Path,
        default=DATA / "levi_fusion_heldout_results.json",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DATA / "generator_compression_trial.json",
    )
    arguments = parser.parse_args()

    mechanisms = load_mechanisms(
        arguments.pilot,
        arguments.heldout,
    )
    payload = {
        "schema_version": 1,
        "status": "early_post_hoc_go_no_go_trial",
        "source_commit": (
            "058e0a34834b2ff3eacfab258d8030264db7d2b3"
        ),
        "input_sha256": {
            str(arguments.pilot): hashlib.sha256(
                arguments.pilot.read_bytes()
            ).hexdigest(),
            str(arguments.heldout): hashlib.sha256(
                arguments.heldout.read_bytes()
            ).hexdigest(),
        },
        "tangent_convention_validation": validate_tangent_convention(
            arguments.pilot
        ),
        "projection_width": projection_width_audit(mechanisms),
        "synthetic_projection_stress": synthetic_projection_stress(),
        "known_mechanism_singleton_tomography": (
            known_mechanism_tomography(mechanisms)
        ),
        "bounded_candidate_search": bounded_candidate_audit(),
    }
    payload["canonical_sha256_without_this_field"] = canonical_hash(
        payload
    )
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
