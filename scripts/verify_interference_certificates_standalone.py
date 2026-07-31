#!/usr/bin/env python3
"""Standalone exact verifier for global no-design certificates.

This script uses only the Python standard library and shares no code with the
certificate generator. It verifies three layers:

1. Every integer Farkas vector has positive pairing with ``(spectrum, 1)``.
2. The vector and its equal-eigenvalue symmetry orbit define exact hitting
   clauses that every nonnegative incidence realization must meet.
3. A finite branch DAG proves that no one-hop-free determinant support can
   meet all of those clauses.

Usage:
    python3 scripts/verify_interference_certificates_standalone.py [artifact]
"""

import json
import math
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product


def incidence_system(n, d, spectrum):
    determinants = list(combinations(range(d), n))
    matrix = [
        [Fraction(int(mode in determinant)) for determinant in determinants]
        for mode in range(d)
    ]
    matrix.append([Fraction(1)] * len(determinants))
    return determinants, matrix, spectrum + [Fraction(1)]


def value_blocks(spectrum):
    grouped = defaultdict(list)
    order = []
    for position, value in enumerate(spectrum):
        if value not in grouped:
            order.append(value)
        grouped[value].append(position)
    return [grouped[value] for value in order]


def symmetry_permutations(spectrum):
    blocks = value_blocks(spectrum)
    choices = [tuple(permutations(block)) for block in blocks]
    output = []
    for images in product(*choices):
        permutation = list(range(len(spectrum)))
        for block, image in zip(blocks, images):
            for source, target in zip(block, image):
                permutation[source] = target
        output.append(tuple(permutation))
    return output


def determinant_permutation(determinants, mode_permutation):
    lookup = {determinant: i for i, determinant in enumerate(determinants)}
    return tuple(
        lookup[tuple(sorted(mode_permutation[mode] for mode in determinant))]
        for determinant in determinants
    )


def move_mask(mask, mapping):
    output = 0
    while mask:
        bit = mask & -mask
        output |= 1 << mapping[bit.bit_length() - 1]
        mask -= bit
    return output


def mask_members(mask):
    members = []
    while mask:
        bit = mask & -mask
        members.append(bit.bit_length() - 1)
        mask -= bit
    return members


def conflict_masks(determinants, n):
    output = [0] * len(determinants)
    for left_index, left in enumerate(determinants):
        left_set = set(left)
        for right_index in range(left_index + 1, len(determinants)):
            if len(left_set & set(determinants[right_index])) == n - 1:
                output[left_index] |= 1 << right_index
                output[right_index] |= 1 << left_index
    return output


def forbidden_by(selected, conflicts):
    forbidden = selected
    pending = selected
    while pending:
        bit = pending & -pending
        forbidden |= conflicts[bit.bit_length() - 1]
        pending -= bit
    return forbidden


def positive_columns(vector, matrix):
    output = 0
    for column in range(len(matrix[0])):
        pairing = sum(
            Fraction(vector[row]) * matrix[row][column]
            for row in range(len(matrix))
        )
        if pairing > 0:
            output |= 1 << column
    return output


def canonical(mask, determinant_maps):
    return min(move_mask(mask, mapping) for mapping in determinant_maps)


def verify_record(record):
    n, d = record["system"]
    spectrum = [Fraction(value) for value in record["spectrum"]]
    if len(spectrum) != d or sum(spectrum) != n:
        return False, "invalid spectrum"
    if record["symmetry_blocks"] != value_blocks(spectrum):
        return False, "incorrect symmetry blocks"
    if "integer_form" in record:
        denominator = record["denominator"]
        expected = [Fraction(value, denominator) for value in record["integer_form"]]
        if spectrum != expected:
            return False, "integer form disagrees with spectrum"

    determinants, matrix, target = incidence_system(n, d, spectrum)
    conflicts = conflict_masks(determinants, n)
    mode_permutations = symmetry_permutations(spectrum)
    determinant_maps = [
        determinant_permutation(determinants, permutation)
        for permutation in mode_permutations
    ]
    map_by_permutation = dict(zip(mode_permutations, determinant_maps))

    positive_sets = []
    all_clauses = set()
    for vector in record["farkas_vectors"]:
        if len(vector) != d + 1:
            return False, "wrong Farkas-vector length"
        divisor = 0
        for value in vector:
            divisor = math.gcd(divisor, abs(value))
        if divisor != 1:
            return False, "Farkas vector is not primitive"
        margin = sum(Fraction(value) * rhs for value, rhs in zip(vector, target))
        if margin <= 0:
            return False, "nonpositive Farkas margin"
        positive = positive_columns(vector, matrix)
        if positive == 0:
            return False, "empty Farkas clause"
        positive_sets.append(positive)
        for mapping in determinant_maps:
            all_clauses.add(move_mask(positive, mapping))
    if len(all_clauses) != record["expanded_clause_count"]:
        return False, "expanded clause count"

    proof = record["proof"]
    nodes = proof["nodes"]
    if proof["node_count"] != len(nodes):
        return False, "proof node count"
    root = proof["root"]
    if not 0 <= root < len(nodes) or nodes[root]["selected"]:
        return False, "invalid proof root"

    active = set()
    checked = set()
    maximum_depth = 0

    def check_node(node_id):
        nonlocal maximum_depth
        if node_id in checked:
            return True, "ok"
        if node_id in active or not 0 <= node_id < len(nodes):
            return False, "proof cycle or invalid node"
        active.add(node_id)
        node = nodes[node_id]
        chosen = node["selected"]
        if chosen != sorted(set(chosen)):
            return False, "selected determinants are not canonical lists"
        if any(not 0 <= value < len(determinants) for value in chosen):
            return False, "selected determinant out of range"
        selected = sum(1 << value for value in chosen)
        maximum_depth = max(maximum_depth, selected.bit_count())
        if any(conflicts[value] & selected for value in chosen):
            return False, "selected support has a one-hop conflict"
        forbidden = forbidden_by(selected, conflicts)

        source = node["clause"]
        vector_index = source["vector"]
        permutation = tuple(source["permutation"])
        if not 0 <= vector_index < len(positive_sets):
            return False, "invalid vector reference"
        mapping = map_by_permutation.get(permutation)
        if mapping is None:
            return False, "invalid symmetry reference"
        clause = move_mask(positive_sets[vector_index], mapping)
        if clause not in all_clauses or clause & selected:
            return False, "branch clause is invalid or already hit"
        available = clause & ~forbidden

        if node["kind"] == "dead":
            if available:
                return False, "dead leaf still has an available determinant"
        elif node["kind"] == "branch":
            children = node.get("children", [])
            candidates = [child["determinant"] for child in children]
            if candidates != mask_members(available):
                return False, "branch does not cover every available determinant"
            for child in children:
                child_id = child["node"]
                if not 0 <= child_id < len(nodes):
                    return False, "invalid child reference"
                target_selected = selected | (1 << child["determinant"])
                child_selected = sum(
                    1 << value for value in nodes[child_id]["selected"]
                )
                if child_selected.bit_count() != selected.bit_count() + 1:
                    return False, "child does not add exactly one determinant"
                if canonical(child_selected, determinant_maps) != canonical(
                    target_selected, determinant_maps
                ):
                    return False, "child is not symmetry-equivalent to its branch"
                ok, reason = check_node(child_id)
                if not ok:
                    return ok, reason
        else:
            return False, "unknown proof node kind"

        active.remove(node_id)
        checked.add(node_id)
        return True, "ok"

    ok, reason = check_node(root)
    if not ok:
        return ok, reason
    if len(checked) != len(nodes):
        return False, "unreachable proof nodes"
    if maximum_depth != proof["max_depth"]:
        return False, "proof depth"
    return True, "ok"


def main():
    path = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "results/data/interference_certificates.json"
    )
    artifact = json.load(open(path))
    if artifact.get("schema_version") != 1:
        print("FAIL: unsupported schema")
        return 1

    records = artifact.get("records", [])
    failures = []
    for record in records:
        ok, reason = verify_record(record)
        label = f"{tuple(record.get('system', []))} index {record.get('vertex_index')}"
        if ok:
            print(f"PASS {label}")
        else:
            print(f"FAIL {label}: {reason}")
            failures.append((label, reason))

    expected = artifact.get("record_count")
    if expected != len(records):
        print("FAIL: record count")
        failures.append(("artifact", "record count"))
    summary = {
        "expanded_clauses": sum(record["expanded_clause_count"] for record in records),
        "farkas_vectors": sum(len(record["farkas_vectors"]) for record in records),
        "maximum_proof_depth": max(
            (record["proof"]["max_depth"] for record in records),
            default=0,
        ),
        "proof_nodes": sum(record["proof"]["node_count"] for record in records),
    }
    if artifact.get("summary") != summary:
        print("FAIL: artifact summary")
        failures.append(("artifact", "summary"))
    if failures:
        return 1
    print(f"all {expected} global interference certificates PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
