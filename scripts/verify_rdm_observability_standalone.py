#!/usr/bin/env python3
"""Independently verify the exact RDM observability census artifact.

This checker intentionally imports no project or generator code.  It verifies
the input hash, every incidence rank certificate, every collision-free tree
witness, every connected tree, and the published aggregate counts.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import sys
from collections import Counter, defaultdict, deque
from fractions import Fraction
from pathlib import Path


def sign(sequence):
    inversions = sum(
        sequence[i] > sequence[j]
        for i in range(len(sequence))
        for j in range(i + 1, len(sequence))
    )
    return -1 if inversions % 2 else 1


def rank(matrix):
    work = [list(map(Fraction, row)) for row in matrix]
    if not work:
        return 0
    pivot = 0
    for column in range(len(work[0])):
        chosen = next(
            (row for row in range(pivot, len(work)) if work[row][column]),
            None,
        )
        if chosen is None:
            continue
        work[pivot], work[chosen] = work[chosen], work[pivot]
        value = work[pivot][column]
        for row in range(pivot + 1, len(work)):
            if not work[row][column]:
                continue
            left = work[row][column]
            work[row] = [
                value * entry - left * basis
                for entry, basis in zip(work[row], work[pivot])
            ]
        pivot += 1
        if pivot == len(work):
            break
    return pivot


def determinant(matrix):
    work = [list(map(Fraction, row)) for row in matrix]
    value = Fraction(1)
    for column in range(len(work)):
        chosen = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if chosen is None:
            return 0
        if chosen != column:
            work[column], work[chosen] = work[chosen], work[column]
            value *= -1
        pivot = work[column][column]
        value *= pivot
        work[column] = [entry / pivot for entry in work[column]]
        for row in range(column + 1, len(work)):
            factor = work[row][column]
            work[row] = [
                entry - factor * basis
                for entry, basis in zip(work[row], work[column])
            ]
    if value.denominator != 1:
        raise AssertionError("integer minor acquired a fractional determinant")
    return int(value)


def channel_terms(support, order, left_modes, right_modes):
    particles = len(support[0])
    remainder_size = particles - order
    coefficients = defaultdict(int)
    for left_index, left in enumerate(support):
        for right_index, right in enumerate(support):
            intersection = sorted(set(left).intersection(right))
            for remainder in itertools.combinations(intersection, remainder_size):
                remainder_set = set(remainder)
                candidate_left = tuple(x for x in left if x not in remainder_set)
                candidate_right = tuple(x for x in right if x not in remainder_set)
                if candidate_left != left_modes or candidate_right != right_modes:
                    continue
                coefficients[(left_index, right_index)] += (
                    sign(candidate_left + remainder) * sign(candidate_right + remainder)
                )
    return {edge: value for edge, value in coefficients.items() if value}


def channel_graph(support, order):
    particles = len(support[0])
    rows = defaultdict(lambda: defaultdict(int))
    for left_index, left in enumerate(support):
        for right_index, right in enumerate(support):
            common = sorted(set(left).intersection(right))
            for remainder in itertools.combinations(common, particles - order):
                rest = set(remainder)
                left_modes = tuple(mode for mode in left if mode not in rest)
                right_modes = tuple(mode for mode in right if mode not in rest)
                if len(left_modes) != order or len(right_modes) != order:
                    continue
                rows[(left_modes, right_modes)][(left_index, right_index)] += (
                    sign(left_modes + remainder) * sign(right_modes + remainder)
                )

    adjacency = defaultdict(set)
    multiplicities = Counter()
    for (left_modes, right_modes), coefficients in rows.items():
        if left_modes >= right_modes:
            continue
        terms = {edge: value for edge, value in coefficients.items() if value}
        if not terms:
            continue
        multiplicities[len(terms)] += 1
        if len(terms) == 1:
            left_index, right_index = next(iter(terms))
            if left_index != right_index:
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

    summary = {
        "order": order,
        "off_diagonal_channels": sum(multiplicities.values()),
        "collision_free_channels": multiplicities[1],
        "collided_channels": sum(
            count for size, count in multiplicities.items() if size > 1
        ),
        "maximum_channel_multiplicity": max(multiplicities, default=0),
        "unique_graph_connected": len(components) == 1,
        "component_sizes": sorted((len(component) for component in components), reverse=True),
    }
    return summary, components


def verify_tree(support, order, tree):
    assert len(tree) == len(support) - 1
    edges = set()
    vertices = {0}
    tree_adjacency = defaultdict(set)
    for witness in tree:
        left_modes = tuple(witness["left_modes"])
        right_modes = tuple(witness["right_modes"])
        terms = channel_terms(support, order, left_modes, right_modes)
        expected = {(witness["left"], witness["right"]): witness["coefficient"]}
        assert terms == expected
        assert {witness["parent"], witness["child"]} == {
            witness["left"],
            witness["right"],
        }
        edge = frozenset((witness["parent"], witness["child"]))
        assert edge not in edges
        edges.add(edge)
        vertices.update(edge)
        tree_adjacency[witness["parent"]].add(witness["child"])
        tree_adjacency[witness["child"]].add(witness["parent"])
    assert vertices == set(range(len(support)))
    reached = {0}
    queue = deque([0])
    while queue:
        vertex = queue.popleft()
        for neighbor in tree_adjacency[vertex]:
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    assert reached == set(range(len(support)))


def verify_record(state, certificate):
    assert (state["system"], state["index"], state["classified"]) == (
        certificate["system"],
        certificate["index"],
        certificate["classified"],
    )
    support = [tuple(row) for row in state["closed_form"]["support_dets"]]
    assert support == [tuple(row) for row in certificate["support_dets"]]
    particles, orbitals = map(int, state["system"].strip("()").split(","))
    assert certificate["particles"] == particles
    assert certificate["orbitals"] == orbitals
    assert certificate["support_size"] == len(support)

    for order, key in ((1, "singleton_incidence"), (2, "pair_incidence")):
        features = list(itertools.combinations(range(orbitals), order))
        matrix = [
            [int(set(feature).issubset(det)) for det in support]
            for feature in features
        ]
        claimed = certificate[key]
        assert claimed["column_count"] == len(support)
        assert rank(matrix) == claimed["rank"]
        assert claimed["full_column_rank"] == (claimed["rank"] == len(support))
        if order == 1:
            assert claimed["kernel_dimension"] == len(support) - claimed["rank"]
        selected = [tuple(row) for row in claimed["pivot_features"]]
        selected_indices = [features.index(feature) for feature in selected]
        minor = [matrix[index] for index in selected_indices]
        assert rank(minor) == claimed["rank"]
        if claimed["full_column_rank"]:
            assert len(minor) == len(support)
            assert determinant(minor) == claimed["pivot_minor_determinant"]
            assert claimed["pivot_minor_determinant"] != 0

    one_hop_pairs = sum(
        len(set(left).intersection(right)) == particles - 1
        for left, right in itertools.combinations(support, 2)
    )
    assert one_hop_pairs == certificate["one_hop_pairs"]
    checked = []
    first_connected = None
    gamma2_components = None
    for order in range(2, particles + 1):
        channel_summary, components = channel_graph(support, order)
        checked.append(channel_summary)
        if order == 2:
            gamma2_components = components
        if channel_summary["unique_graph_connected"]:
            first_connected = order
            break
    assert checked == certificate["orders_checked"]
    assert checked[0] == certificate["gamma2_channels"]
    assert gamma2_components == certificate["gamma2_components"]
    assert first_connected == certificate["first_connected_unique_channel_order"]
    expected_reconstruction = (
        certificate["pair_incidence"]["full_column_rank"]
        and checked[0]["unique_graph_connected"]
    )
    assert expected_reconstruction == certificate[
        "pure_full_support_reconstruction_from_gamma2"
    ]
    first_tree = certificate["first_connected_reconstruction_tree"]
    verify_tree(support, first_connected, first_tree)
    gamma2_tree = certificate["gamma2_reconstruction_tree"]
    if certificate["pure_full_support_reconstruction_from_gamma2"]:
        assert gamma2_tree == first_tree
        verify_tree(support, 2, gamma2_tree)
    else:
        assert gamma2_tree == []
    return particles


def main():
    states_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/data/states.jsonl")
    artifact_path = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else Path("results/data/rdm_observability_census.json")
    )
    raw = states_path.read_bytes()
    states = [json.loads(line) for line in raw.splitlines() if line.strip()]
    artifact = json.loads(artifact_path.read_text())
    assert artifact["evidence"]["status"] == "EXACT"
    assert artifact["input"]["sha256"] == hashlib.sha256(raw).hexdigest()
    assert len(states) == len(artifact["records"]) == 799

    for state, certificate in zip(states, artifact["records"]):
        verify_record(state, certificate)

    rows = artifact["records"]
    summary = artifact["summary"]
    classes = sorted({row["classified"] for row in rows})
    expected_by_class = {}
    for class_name in classes:
        class_rows = [row for row in rows if row["classified"] == class_name]
        expected_by_class[class_name] = {
            "records": len(class_rows),
            "gamma2_reconstructible": sum(
                row["pure_full_support_reconstruction_from_gamma2"]
                for row in class_rows
            ),
            "gamma2_unresolved": sum(
                not row["pure_full_support_reconstruction_from_gamma2"]
                for row in class_rows
            ),
            "positive_singleton_kernel": sum(
                row["singleton_incidence"]["kernel_dimension"] > 0
                for row in class_rows
            ),
        }
    kernel_dimensions = Counter(
        row["singleton_incidence"]["kernel_dimension"]
        for row in rows
        if row["singleton_incidence"]["kernel_dimension"] > 0
    )
    first_orders = Counter(
        row["first_connected_unique_channel_order"] for row in rows
    )
    expected_summary = {
        "records": len(rows),
        "pair_incidence_full_column_rank": sum(
            row["pair_incidence"]["full_column_rank"] for row in rows
        ),
        "gamma2_reconstructible": sum(
            row["pure_full_support_reconstruction_from_gamma2"] for row in rows
        ),
        "gamma2_unresolved_by_this_certificate": sum(
            not row["pure_full_support_reconstruction_from_gamma2"] for row in rows
        ),
        "positive_singleton_kernel": sum(kernel_dimensions.values()),
        "positive_singleton_kernel_by_dimension": {
            str(dimension): count
            for dimension, count in sorted(kernel_dimensions.items())
        },
        "one_hop_pairs": sum(row["one_hop_pairs"] for row in rows),
        "first_connected_unique_channel_order": {
            str(order): count for order, count in sorted(first_orders.items())
        },
        "by_class": expected_by_class,
    }
    assert summary == expected_summary
    assert summary["pair_incidence_full_column_rank"] == 799
    assert summary["gamma2_reconstructible"] == 771
    assert summary["gamma2_unresolved_by_this_certificate"] == 28
    assert summary["positive_singleton_kernel"] == 63
    assert summary["one_hop_pairs"] == 305
    assert summary["positive_singleton_kernel_by_dimension"] == {
        "1": 47,
        "2": 15,
        "5": 1,
    }
    assert summary["first_connected_unique_channel_order"] == {
        "2": 771,
        "3": 22,
        "4": 5,
        "5": 1,
    }
    assert Counter(
        row["classified"]
        for row in rows
        if row["pure_full_support_reconstruction_from_gamma2"]
    ) == {"DESIGN-INT": 603, "DESIGN-REAL": 12, "INTERFERENCE": 156}
    assert all(
        row["classified"] == "INTERFERENCE"
        for row in rows
        if row["singleton_incidence"]["kernel_dimension"] > 0
    )
    print("all 799 RDM observability certificates PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
