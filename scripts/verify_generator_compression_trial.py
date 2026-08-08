#!/usr/bin/env python3
"""Independent lightweight verifier for the compression-trial artifact."""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path


def determinant_modulo(entries, values, prime):
    size = len(entries)
    matrix = [
        [
            0 if entry is None else entry[0] * values[entry[1]] % prime
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
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
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


def symbolic_matrix(tau, particle_number):
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
    entries = [[None] * len(roots) for _ in below]
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


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    artifact_path = root / "results" / "data" / "generator_compression_trial.json"
    payload = json.loads(artifact_path.read_text())

    digest = payload.pop("canonical_sha256_without_this_field")
    recomputed = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if digest != recomputed:
        raise AssertionError("canonical artifact digest mismatch")

    projection = payload["projection_width"]
    assert projection["mechanism_count"] == 83
    assert projection["max_optimal_tree_width"] == 5
    assert projection["max_balanced_tree_width"] == 6

    bounded = payload["bounded_candidate_search"]
    assert bounded["total_trace_placements"] == 65182
    assert bounded["total_hall_feasible_placements"] == 537
    assert bounded["multi_hole_trace_placements"] == 135
    assert bounded["multi_hole_hall_feasible_placements"] == 0
    assert bounded["nonstructural_singleton_cross_shape_collisions"] == 0

    certificates = bounded["one_hole_modular_nonzero_certificates"]
    assert len(certificates) == 3
    for certificate in certificates:
        tau = tuple(certificate["tau"])
        on, below, roots, entries = symbolic_matrix(tau, 3)
        assert len(roots) == certificate["matrix_size"]
        assert len(on) == certificate["variables"]
        assert len(below) == len(roots)
        residue = determinant_modulo(
            entries,
            certificate["variable_values"],
            certificate["prime"],
        )
        assert residue == certificate["determinant_residue"]
        assert residue != 0

    validation = payload["tangent_convention_validation"]
    assert validation["nonstructural_pilot_facets_checked"] == 49
    assert validation["passed_trace_hall_and_modular_nonzero"] == 49
    assert validation["failures"] == []

    print("compression trial verifier: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
