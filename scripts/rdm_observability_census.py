#!/usr/bin/env python3
"""Build exact support certificates for pure-state reconstruction from RDMs.

For a determinant carrier S, diagonal entries of the 2-RDM are the pair
incidence map applied to the coefficient weights.  Off-diagonal k-RDM entries
are signed sums of coefficient coherences.  This script uses integer and
rational arithmetic only to certify two sufficient conditions:

1. pair incidence has full column rank, so the 2-RDM recovers every weight;
2. collision-free 2-RDM coherences contain a connected spanning tree, so the
   2-RDM recovers every coefficient coherence on the pure full-support locus.

The emitted certificate is support-combinatorial.  It does not use floating
amplitudes, infer a mixed-state inverse, or assert that an unresolved support
is non-injective.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter, defaultdict, deque
from fractions import Fraction
from pathlib import Path


GENERATOR = "scripts/rdm_observability_census.py"


def parse_system(system: str) -> tuple[int, int]:
    particles, orbitals = system.strip("()").split(",")
    return int(particles), int(orbitals)


def permutation_sign(sequence: tuple[int, ...]) -> int:
    inversions = sum(
        sequence[i] > sequence[j]
        for i in range(len(sequence))
        for j in range(i + 1, len(sequence))
    )
    return -1 if inversions % 2 else 1


def independent_rows(matrix: list[list[int]]) -> list[int]:
    """Return original indices of an exact row basis."""
    if not matrix:
        return []
    work = [list(map(Fraction, row)) for row in matrix]
    labels = list(range(len(work)))
    pivot = 0
    for column in range(len(work[0])):
        chosen = next(
            (row for row in range(pivot, len(work)) if work[row][column]),
            None,
        )
        if chosen is None:
            continue
        work[pivot], work[chosen] = work[chosen], work[pivot]
        labels[pivot], labels[chosen] = labels[chosen], labels[pivot]
        value = work[pivot][column]
        work[pivot] = [entry / value for entry in work[pivot]]
        for row in range(pivot + 1, len(work)):
            if not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                entry - factor * basis
                for entry, basis in zip(work[row], work[pivot])
            ]
        pivot += 1
        if pivot == len(work):
            break
    return labels[:pivot]


def determinant_bareiss(matrix: list[list[int]]) -> int:
    """Compute an integer determinant without division roundoff."""
    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        chosen = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if chosen is None:
            return 0
        if chosen != column:
            work[column], work[chosen] = work[chosen], work[column]
            sign *= -1
        pivot = work[column][column]
        for row in range(column + 1, size):
            for other in range(column + 1, size):
                numerator = (
                    work[row][other] * pivot
                    - work[row][column] * work[column][other]
                )
                if numerator % previous:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][other] = numerator // previous
        previous = pivot
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


def incidence_certificate(
    support: list[tuple[int, ...]], orbitals: int, order: int
) -> dict:
    features = list(itertools.combinations(range(orbitals), order))
    matrix = [
        [int(set(feature).issubset(determinant)) for determinant in support]
        for feature in features
    ]
    pivots = independent_rows(matrix)
    rank = len(pivots)
    certificate = {
        "rank": rank,
        "column_count": len(support),
        "full_column_rank": rank == len(support),
        "pivot_features": [list(features[index]) for index in pivots],
    }
    if rank == len(support):
        minor = [matrix[index] for index in pivots]
        certificate["pivot_minor_determinant"] = determinant_bareiss(minor)
    return certificate


def signed_channels(
    support: list[tuple[int, ...]], order: int
) -> dict[tuple[tuple[int, ...], tuple[int, ...]], dict[tuple[int, int], int]]:
    """Return exact oriented k-RDM channels in coefficient-density variables."""
    particles = len(support[0])
    terms: dict[
        tuple[tuple[int, ...], tuple[int, ...]], dict[tuple[int, int], int]
    ] = defaultdict(lambda: defaultdict(int))
    for left_index, left in enumerate(support):
        left_set = set(left)
        for right_index, right in enumerate(support):
            intersection = sorted(left_set.intersection(right))
            remainder_size = particles - order
            if len(intersection) < remainder_size:
                continue
            for remainder in itertools.combinations(intersection, remainder_size):
                remainder_set = set(remainder)
                left_modes = tuple(mode for mode in left if mode not in remainder_set)
                right_modes = tuple(mode for mode in right if mode not in remainder_set)
                if len(left_modes) != order or len(right_modes) != order:
                    continue
                sign = permutation_sign(left_modes + remainder)
                sign *= permutation_sign(right_modes + remainder)
                terms[(left_modes, right_modes)][(left_index, right_index)] += sign
    reduced = {}
    for channel, coefficients in terms.items():
        nonzero = {edge: value for edge, value in coefficients.items() if value}
        if nonzero:
            reduced[channel] = nonzero
    return reduced


def unique_channel_graph(
    support: list[tuple[int, ...]], order: int
) -> tuple[dict, list[dict], list[list[int]]]:
    channels = signed_channels(support, order)
    witnesses = []
    adjacency: dict[int, set[int]] = defaultdict(set)
    multiplicities = Counter()
    for (left_modes, right_modes), terms in sorted(channels.items()):
        if left_modes >= right_modes or left_modes == right_modes:
            continue
        multiplicities[len(terms)] += 1
        if len(terms) != 1:
            continue
        (left_index, right_index), coefficient = next(iter(terms.items()))
        if left_index == right_index:
            continue
        witness = {
            "left": left_index,
            "right": right_index,
            "left_modes": list(left_modes),
            "right_modes": list(right_modes),
            "coefficient": coefficient,
        }
        witnesses.append(witness)
        adjacency[left_index].add(right_index)
        adjacency[right_index].add(left_index)

    unseen = set(range(len(support)))
    components = []
    while unseen:
        root = min(unseen)
        queue = deque([root])
        unseen.remove(root)
        component = []
        while queue:
            vertex = queue.popleft()
            component.append(vertex)
            for neighbor in sorted(adjacency[vertex]):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    queue.append(neighbor)
        components.append(sorted(component))

    tree = []
    if len(components) == 1:
        by_edge: dict[tuple[int, int], list[dict]] = defaultdict(list)
        for witness in witnesses:
            edge = tuple(sorted((witness["left"], witness["right"])))
            by_edge[edge].append(witness)
        root = 0
        visited = {root}
        queue = deque([root])
        while queue:
            parent = queue.popleft()
            candidates = []
            for edge, edge_witnesses in by_edge.items():
                if parent not in edge:
                    continue
                child = edge[1] if edge[0] == parent else edge[0]
                if child in visited:
                    continue
                candidates.append((child, edge_witnesses[0]))
            for child, witness in sorted(candidates, key=lambda item: item[0]):
                if child in visited:
                    continue
                visited.add(child)
                queue.append(child)
                tree.append({"parent": parent, "child": child, **witness})

    channel_summary = {
        "order": order,
        "off_diagonal_channels": sum(multiplicities.values()),
        "collision_free_channels": multiplicities[1],
        "collided_channels": sum(
            count for multiplicity, count in multiplicities.items() if multiplicity > 1
        ),
        "maximum_channel_multiplicity": max(multiplicities, default=0),
        "unique_graph_connected": len(components) == 1,
        "component_sizes": sorted((len(component) for component in components), reverse=True),
    }
    return channel_summary, tree, components


def record_certificate(record: dict) -> dict:
    particles, orbitals = parse_system(record["system"])
    support = [tuple(row) for row in record["closed_form"]["support_dets"]]
    singleton = incidence_certificate(support, orbitals, 1)
    pair = incidence_certificate(support, orbitals, 2)

    order_certificates = []
    tree = []
    first_connected_order = None
    gamma2_components = None
    for order in range(2, particles + 1):
        channel_summary, candidate_tree, components = unique_channel_graph(support, order)
        order_certificates.append(channel_summary)
        if order == 2:
            gamma2_components = components
        if channel_summary["unique_graph_connected"]:
            first_connected_order = order
            tree = candidate_tree
            break

    gamma2 = order_certificates[0]
    reconstructible = pair["full_column_rank"] and gamma2["unique_graph_connected"]
    one_hop_pairs = sum(
        len(set(left).intersection(right)) == particles - 1
        for left, right in itertools.combinations(support, 2)
    )
    return {
        "system": record["system"],
        "index": record["index"],
        "classified": record["classified"],
        "particles": particles,
        "orbitals": orbitals,
        "support_size": len(support),
        "support_dets": [list(determinant) for determinant in support],
        "one_hop_pairs": one_hop_pairs,
        "singleton_incidence": {
            **singleton,
            "kernel_dimension": len(support) - singleton["rank"],
        },
        "pair_incidence": pair,
        "gamma2_channels": gamma2,
        "gamma2_components": gamma2_components,
        "first_connected_unique_channel_order": first_connected_order,
        "orders_checked": order_certificates,
        "pure_full_support_reconstruction_from_gamma2": reconstructible,
        "gamma2_reconstruction_tree": tree if reconstructible else [],
        "first_connected_reconstruction_tree": tree,
    }


def summarize(records: list[dict]) -> dict:
    by_class = {}
    for class_name in sorted({record["classified"] for record in records}):
        rows = [record for record in records if record["classified"] == class_name]
        by_class[class_name] = {
            "records": len(rows),
            "gamma2_reconstructible": sum(
                record["pure_full_support_reconstruction_from_gamma2"] for record in rows
            ),
            "gamma2_unresolved": sum(
                not record["pure_full_support_reconstruction_from_gamma2"] for record in rows
            ),
            "positive_singleton_kernel": sum(
                record["singleton_incidence"]["kernel_dimension"] > 0 for record in rows
            ),
        }
    order_counts = Counter(
        record["first_connected_unique_channel_order"] for record in records
    )
    return {
        "records": len(records),
        "pair_incidence_full_column_rank": sum(
            record["pair_incidence"]["full_column_rank"] for record in records
        ),
        "gamma2_reconstructible": sum(
            record["pure_full_support_reconstruction_from_gamma2"] for record in records
        ),
        "gamma2_unresolved_by_this_certificate": sum(
            not record["pure_full_support_reconstruction_from_gamma2"]
            for record in records
        ),
        "positive_singleton_kernel": sum(
            record["singleton_incidence"]["kernel_dimension"] > 0 for record in records
        ),
        "positive_singleton_kernel_by_dimension": {
            str(dimension): count
            for dimension, count in sorted(
                Counter(
                    record["singleton_incidence"]["kernel_dimension"]
                    for record in records
                    if record["singleton_incidence"]["kernel_dimension"] > 0
                ).items()
            )
        },
        "one_hop_pairs": sum(record["one_hop_pairs"] for record in records),
        "first_connected_unique_channel_order": {
            str(order): count for order, count in sorted(order_counts.items())
        },
        "by_class": by_class,
    }


def build(states_path: Path) -> dict:
    raw = states_path.read_bytes()
    states = [json.loads(line) for line in raw.splitlines() if line.strip()]
    records = [record_certificate(record) for record in states]
    return {
        "schema_version": 1,
        "generated_by": GENERATOR,
        "input": {
            "path": "results/data/states.jsonl",
            "sha256": hashlib.sha256(raw).hexdigest(),
        },
        "evidence": {
            "status": "EXACT",
            "scope": (
                "Integer support combinatorics and exact rational row reduction. "
                "The inverse applies to pure states with every listed carrier amplitude nonzero."
            ),
        },
        "conventions": {
            "pair": "P=(p_1<p_2) and A_P=a_{p_2}a_{p_1}",
            "rdm": "Gamma_k[P,Q]=expectation(A_P^dagger A_Q)",
            "carrier_matrix": "X_ab=conjugate(c_a)c_b",
            "fermionic_sign": "epsilon(K,R)=(-1)^inversions(K concatenated with R)",
            "orbital_labels": "zero-based",
        },
        "theorem": {
            "name": "collision-free pure-state RDM reconstruction",
            "hypotheses": [
                "The pair-incidence matrix has full column rank.",
                "Collision-free 2-RDM coherence witnesses connect the support.",
                "The state is pure and has nonzero amplitude on every support determinant.",
                "The orbital frame and determinant carrier are fixed.",
            ],
            "conclusion": (
                "The 2-RDM rationally reconstructs the carrier rank-one matrix "
                "X_ab=conjugate(c_a)c_b and therefore every higher RDM."
            ),
            "nonclaims": [
                "No mixed-state inverse is asserted.",
                "An unresolved record is not asserted to be non-injective.",
                "No carrier-free or orbital-quotient inverse is asserted.",
                "No Hamiltonian is asserted to preserve a listed carrier.",
            ],
        },
        "summary": summarize(records),
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "states",
        nargs="?",
        type=Path,
        default=Path("results/data/states.jsonl"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/data/rdm_observability_census.json"),
    )
    args = parser.parse_args()
    artifact = build(args.states)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    summary = artifact["summary"]
    print(
        f"wrote {args.out}: {summary['gamma2_reconstructible']}/"
        f"{summary['records']} gamma2-reconstructible"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
