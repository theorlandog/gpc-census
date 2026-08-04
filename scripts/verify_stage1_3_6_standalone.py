#!/usr/bin/env python3
"""Dependency-free replay of the checked-in (3,6) Ressayre artifact."""

from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT = ROOT / "results" / "data" / "stage1_ressayre_3_6.json"


def weights() -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(int(index in occupied) for index in range(6))
        for occupied in itertools.combinations(range(6), 3)
    )


def matrix_rank(matrix: list[list[int]]) -> int:
    data = [[Fraction(value) for value in row] for row in matrix]
    if not data:
        return 0
    row = 0
    for column in range(len(data[0])):
        pivot = next(
            (index for index in range(row, len(data)) if data[index][column]), None
        )
        if pivot is None:
            continue
        data[row], data[pivot] = data[pivot], data[row]
        value = data[row][column]
        data[row] = [entry / value for entry in data[row]]
        for index in range(len(data)):
            if index == row or not data[index][column]:
                continue
            factor = data[index][column]
            data[index] = [
                entry - factor * pivot_entry
                for entry, pivot_entry in zip(data[index], data[row], strict=True)
            ]
        row += 1
        if row == len(data):
            break
    return row


def determinant(matrix: list[list[int]]) -> int:
    data = [list(row) for row in matrix]
    size = len(data)
    if size == 0:
        return 1
    if any(len(row) != size for row in data):
        raise AssertionError("nonsquare tangent matrix")
    if size == 1:
        return data[0][0]
    sign = 1
    previous = 1
    for column in range(size - 1):
        if data[column][column] == 0:
            pivot = next(
                (row for row in range(column + 1, size) if data[row][column]), None
            )
            if pivot is None:
                return 0
            data[column], data[pivot] = data[pivot], data[column]
            sign = -sign
        pivot_value = data[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                data[row][target] = (
                    data[row][target] * pivot_value
                    - data[row][column] * data[column][target]
                ) // previous
        previous = pivot_value
    return sign * data[-1][-1]


def root_action(
    weight: tuple[int, ...], root: tuple[int, int]
) -> tuple[tuple[int, ...], int] | None:
    i, j = root
    if not weight[i] or weight[j]:
        return None
    occupied = [index for index, value in enumerate(weight) if value]
    position = occupied.index(i)
    sign = -1 if position % 2 else 1
    remaining = occupied[:position] + occupied[position + 1 :]
    insertion = sum(index < j for index in remaining)
    if insertion % 2:
        sign = -sign
    target = remaining[:insertion] + [j] + remaining[insertion:]
    return tuple(int(index in target) for index in range(6)), sign


def verify_certificate(row: dict[str, object], exterior: tuple[tuple[int, ...], ...]) -> None:
    h = tuple(int(value) for value in row["h"])
    z = int(row["z"])
    assert len(h) == 6 and sum(h) == 0
    assert math.gcd(*(abs(value) for value in h + (z,))) == 1

    levels = [sum(a * b for a, b in zip(weight, h, strict=True)) for weight in exterior]
    on = tuple(index for index, value in enumerate(levels) if value == z)
    below = tuple(index for index, value in enumerate(levels) if value < z)
    roots = tuple(
        (i, j) for i in range(6) for j in range(i + 1, 6) if h[j] - h[i] < 0
    )
    assert on == tuple(row["on_indices"])
    assert below == tuple(row["below_indices"])
    assert roots == tuple(tuple(root) for root in row["negative_roots"])
    assert len(below) == len(roots)

    witness = tuple(row["admissibility_witness"])
    assert len(witness) == 5 and set(witness) <= set(on)
    assert matrix_rank([list(exterior[index]) + [1] for index in witness]) == 5
    assert matrix_rank([list(exterior[index]) + [1] for index in on]) == 5

    evaluation = tuple(int(value) for value in row["determinant_evaluation"])
    assert len(evaluation) == len(on)
    below_rows = {exterior[index]: target for target, index in enumerate(below)}
    tangent = [[0 for _root in roots] for _weight in below]
    for column, root in enumerate(roots):
        for variable, weight_index in enumerate(on):
            action = root_action(exterior[weight_index], root)
            if action is not None and action[0] in below_rows:
                tangent[below_rows[action[0]]][column] += action[1] * evaluation[variable]
    value = determinant(tangent)
    assert value == int(row["determinant_value"])
    assert value != 0


def verify(path: Path = DEFAULT_ARTIFACT) -> None:
    data = json.loads(path.read_text())
    assert data["system"] == "(3,6)"
    assert data["generator_inputs"]["stored_constraints_used"] is False
    assert data["generator_inputs"]["stored_vertices_used"] is False
    assert data["enumeration"] == {
        "admissible_hyperplanes": 362,
        "oriented_candidates": 724,
        "ressayre_elements": 21,
        "trace_candidates": 64,
    }

    exterior = weights()
    certificates = data["ressayre_certificates"]
    assert len(certificates) == 21
    for certificate in certificates:
        verify_certificate(certificate, exterior)

    assert data["generated_system"]["equalities"] == [
        {"coeffs": [1, 0, 0, 0, 0, 1], "rhs": 1},
        {"coeffs": [0, 1, 0, 0, 1, 0], "rhs": 1},
        {"coeffs": [0, 0, 1, 1, 0, 0], "rhs": 1},
    ]
    assert data["generated_system"]["inequalities"] == [
        {"coeffs": [0, 0, 0, 1, -1, -1], "rhs": 0}
    ]
    assert data["inner_outer_sandwich"]["outer_vertices"] == [
        ["1", "1", "1", "0", "0", "0"],
        ["1", "1/2", "1/2", "1/2", "1/2", "0"],
        ["3/4", "3/4", "1/2", "1/2", "1/4", "1/4"],
        ["1/2", "1/2", "1/2", "1/2", "1/2", "1/2"],
    ]
    assert data["inner_outer_sandwich"]["every_outer_vertex_is_plethysm_certified"]

    target = next(
        row
        for row in certificates
        if row["h"] == [-1, -1, 1, -1, 1, 1] and row["z"] == -1
    )
    assert target["below_indices"] == [1]
    assert target["negative_roots"] == [[2, 3]]
    assert target["determinant_value"] == 1


if __name__ == "__main__":
    verify()
    print("standalone Stage 1 (3,6) verification passed")
