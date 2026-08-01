#!/usr/bin/env python3
"""Independently verify exhaustion of the supplied polytope vertices.

This verifier has no project or third-party runtime dependency.  It rebuilds
each polytope from the supplied rational H-representation by an exact
incremental half-space algorithm and compares the resulting V-representation
with ``results/data/vertices``.

The completeness step is the elementary polytope identity

    vertices(P intersect H)
      = retained vertices(P)
        union {edge(P) intersect boundary(H): edge crosses boundary(H)}.

For a polytope given by all half-spaces processed so far, two distinct
vertices form an edge exactly when their common active rows, together with the
affine equalities, have coefficient rank d - 1.  Every solve, comparison, and
rank calculation below uses exact integer or ``fractions.Fraction`` arithmetic.

Run from any directory:

    python scripts/verify_vertex_exhaustion_standalone.py

The optional ``--system N,D`` flag restricts the check to one supplied system.
``--constraints`` and ``--vertices`` make the verifier portable outside a
source checkout, including from the archived verification bundle.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONSTRAINTS = ROOT / "src" / "gpc_census" / "data" / "constraints.json"
DEFAULT_VERTICES = ROOT / "results" / "data" / "vertices"
SYSTEMS = (
    (3, 6),
    (3, 7),
    (3, 8),
    (4, 8),
    (3, 9),
    (4, 9),
    (3, 10),
    (4, 10),
    (5, 10),
)


Vector = tuple[Fraction, ...]
Row = tuple[tuple[int, ...], int]


@dataclass(frozen=True)
class Vertex:
    point: Vector
    active: int


def _rank(rows: Iterable[Sequence[int | Fraction]], columns: int) -> int:
    """Exact row rank by fraction-free-input Gaussian elimination."""
    matrix = [[Fraction(value) for value in row] for row in rows]
    rank = 0
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [value / scale for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            scale = matrix[index][column]
            matrix[index] = [
                left - scale * right
                for left, right in zip(matrix[index], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def _solve_unique(rows: Sequence[Row], d: int) -> Vector | None:
    """Solve d affine equations in d variables, or return None if singular."""
    if len(rows) != d:
        raise ValueError(f"expected {d} equations, got {len(rows)}")
    matrix = [
        [Fraction(value) for value in coefficients] + [Fraction(rhs)]
        for coefficients, rhs in rows
    ]
    for column in range(d):
        pivot = next(
            (index for index in range(column, d) if matrix[index][column]),
            None,
        )
        if pivot is None:
            return None
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = matrix[column][column]
        matrix[column] = [value / scale for value in matrix[column]]
        for index in range(d):
            if index == column or not matrix[index][column]:
                continue
            scale = matrix[index][column]
            matrix[index] = [
                left - scale * right
                for left, right in zip(matrix[index], matrix[column])
            ]
    return tuple(matrix[index][-1] for index in range(d))


def _evaluate(row: Row, point: Vector) -> Fraction:
    coefficients, _rhs = row
    return sum(Fraction(coefficient) * value for coefficient, value in zip(coefficients, point))


def _scaled_point(
    point: Vector, cache: dict[Vector, tuple[tuple[int, ...], int]]
) -> tuple[tuple[int, ...], int]:
    """Return one common-denominator integer representation of a point."""
    cached = cache.get(point)
    if cached is not None:
        return cached
    denominator = math.lcm(*(value.denominator for value in point))
    numerators = tuple(
        value.numerator * (denominator // value.denominator) for value in point
    )
    result = (numerators, denominator)
    cache[point] = result
    return result


def _slack(
    row: Row,
    point: Vector,
    cache: dict[Vector, tuple[tuple[int, ...], int]],
) -> Fraction:
    """Exact rhs minus row value, using cached integer dot products."""
    numerators, denominator = _scaled_point(point, cache)
    coefficients, rhs = row
    numerator = rhs * denominator - sum(
        coefficient * value for coefficient, value in zip(coefficients, numerators)
    )
    return Fraction(numerator, denominator)


def _is_feasible(point: Vector, rows: Sequence[Row]) -> bool:
    return all(_evaluate(row, point) <= row[1] for row in rows)


def _independent_equalities(equalities: Sequence[Row], d: int) -> list[Row]:
    """Select a coefficient-row basis, checking discarded equations."""
    basis: list[Row] = []
    coefficient_rank = 0
    for row in equalities:
        candidate_rank = _rank(
            [coefficients for coefficients, _rhs in basis + [row]], d
        )
        if candidate_rank > coefficient_rank:
            basis.append(row)
            coefficient_rank = candidate_rank
            continue
        augmented_rank = _rank(
            [tuple(coefficients) + (rhs,) for coefficients, rhs in basis + [row]],
            d + 1,
        )
        if augmented_rank > coefficient_rank:
            raise AssertionError("inconsistent affine equalities")
    return basis


def _base_rows(d: int) -> list[Row]:
    rows: list[Row] = []
    first = [0] * d
    first[0] = 1
    rows.append((tuple(first), 1))
    for index in range(d - 1):
        order = [0] * d
        order[index] = -1
        order[index + 1] = 1
        rows.append((tuple(order), 0))
    last = [0] * d
    last[-1] = -1
    rows.append((tuple(last), 0))
    return rows


def _initial_vertices(equalities: Sequence[Row], rows: Sequence[Row], d: int) -> list[Vertex]:
    affine_dimension = d - len(equalities)
    points: dict[Vector, int] = {}
    for chosen in itertools.combinations(range(len(rows)), affine_dimension):
        point = _solve_unique(
            list(equalities) + [rows[index] for index in chosen], d
        )
        if point is None or not _is_feasible(point, rows):
            continue
        active = 0
        for index, row in enumerate(rows):
            if _evaluate(row, point) == row[1]:
                active |= 1 << index
        points[point] = active
    if not points:
        raise AssertionError("base polytope has no vertices")
    return [Vertex(point, active) for point, active in sorted(points.items(), reverse=True)]


def _iter_bits(bits: int):
    while bits:
        least = bits & -bits
        yield least.bit_length() - 1
        bits ^= least


def _add_halfspace(
    vertices: Sequence[Vertex],
    processed_rows: list[Row],
    new_row: Row,
    equality_normals: Sequence[Sequence[int]],
    d: int,
    rank_cache: dict[int, int],
    point_cache: dict[Vector, tuple[tuple[int, ...], int]],
) -> list[Vertex]:
    """Exact double-description update by crossing old edges."""
    row_index = len(processed_rows)
    boundary_bit = 1 << row_index
    signed = [_slack(new_row, vertex.point, point_cache) for vertex in vertices]
    kept_indices = [index for index, slack in enumerate(signed) if slack >= 0]
    inside_indices = [index for index, slack in enumerate(signed) if slack > 0]
    outside_indices = [index for index, slack in enumerate(signed) if slack < 0]
    if not kept_indices:
        raise AssertionError("half-space update produced an empty polytope")

    updated: dict[Vector, int] = {}
    for index in kept_indices:
        vertex = vertices[index]
        active = vertex.active | (boundary_bit if signed[index] == 0 else 0)
        updated[vertex.point] = active

    minimum_common = d - len(equality_normals) - 1

    def common_rank(bits: int) -> int:
        cached = rank_cache.get(bits)
        if cached is not None:
            return cached
        normals = list(equality_normals)
        normals.extend(processed_rows[index][0] for index in _iter_bits(bits))
        value = _rank(normals, d)
        rank_cache[bits] = value
        return value

    for inside_index in inside_indices:
        inside = vertices[inside_index]
        for outside_index in outside_indices:
            outside = vertices[outside_index]
            common = inside.active & outside.active
            if common.bit_count() < minimum_common or common_rank(common) != d - 1:
                continue
            inside_slack = signed[inside_index]
            outside_slack = signed[outside_index]
            parameter = inside_slack / (inside_slack - outside_slack)
            point = tuple(
                left + parameter * (right - left)
                for left, right in zip(inside.point, outside.point)
            )
            if not (0 < parameter < 1):
                raise AssertionError("invalid crossing parameter")
            active = common | boundary_bit
            previous = updated.get(point)
            if previous is not None and previous != active:
                # Coincident intersections from distinct old edges can expose
                # additional old rows.  Recompute the complete active set.
                active = 0
                for index, row in enumerate(processed_rows):
                    if _evaluate(row, point) == row[1]:
                        active |= 1 << index
                active |= boundary_bit
            updated[point] = active

    processed_rows.append(new_row)
    return [Vertex(point, active) for point, active in sorted(updated.items(), reverse=True)]


def _load_supplied(n: int, d: int, vertices_dir: Path) -> set[Vector]:
    path = vertices_dir / f"vertices_{n}_{d}.json"
    records = json.loads(path.read_text(encoding="utf-8"))
    points: set[Vector] = set()
    for record in records:
        point = tuple(Fraction(value) for value in record["spectrum"])
        denominator = int(record["denominator"])
        integer_form = tuple(int(value) for value in record["integer_form"])
        if point != tuple(Fraction(value, denominator) for value in integer_form):
            raise AssertionError(f"malformed integer form in {path.name}")
        if point in points:
            raise AssertionError(f"duplicate supplied vertex in {path.name}")
        points.add(point)
    return points


def verify_system(
    n: int,
    d: int,
    table: dict[str, dict],
    *,
    vertices_dir: Path = DEFAULT_VERTICES,
    verbose: bool = False,
) -> tuple[int, int]:
    system = table[f"{n},{d}"]
    equalities: list[Row] = [((1,) * d, n)]
    equalities.extend(
        (tuple(int(value) for value in row["coeffs"]), int(row["rhs"]))
        for row in system["equalities"]
    )
    equalities = _independent_equalities(equalities, d)
    equality_normals = [coefficients for coefficients, _rhs in equalities]

    processed_rows = _base_rows(d)
    vertices = _initial_vertices(equalities, processed_rows, d)
    initial_count = len(vertices)
    rank_cache: dict[int, int] = {}
    point_cache: dict[Vector, tuple[tuple[int, ...], int]] = {}
    gpc_rows: list[Row] = [
        (tuple(int(value) for value in row["coeffs"]), int(row["rhs"]))
        for row in system["inequalities"]
    ]
    pending_rows = list(enumerate(gpc_rows))
    for step in range(1, len(gpc_rows) + 1):
        # Half-space order does not affect the final intersection.  Cutting the
        # largest number of current vertices first keeps exact intermediate
        # V-representations small, especially for the (5,10) system.
        best_position = max(
            range(len(pending_rows)),
            key=lambda position: (
                sum(
                    _slack(pending_rows[position][1], vertex.point, point_cache) < 0
                    for vertex in vertices
                ),
                -pending_rows[position][0],
            ),
        )
        source_index, new_row = pending_rows.pop(best_position)
        vertices = _add_halfspace(
            vertices,
            processed_rows,
            new_row,
            equality_normals,
            d,
            rank_cache,
            point_cache,
        )
        if verbose:
            print(
                f"  ({n},{d}) row {step}/{len(gpc_rows)}: "
                f"source {source_index + 1}, {len(vertices)} vertices, "
                f"{len(rank_cache)} rank tests",
                flush=True,
            )

    # Final independent checks do not trust the incrementally maintained
    # active masks: recompute feasibility and vertex rank from all H rows.
    computed: set[Vector] = set()
    for vertex in vertices:
        if any(_slack(row, vertex.point, point_cache) < 0 for row in processed_rows):
            raise AssertionError("incremental enumeration emitted an infeasible point")
        if any(_evaluate(row, vertex.point) != row[1] for row in equalities):
            raise AssertionError("incremental enumeration violated an equality")
        active_normals = list(equality_normals)
        active_normals.extend(
            row[0]
            for row in processed_rows
            if _slack(row, vertex.point, point_cache) == 0
        )
        if _rank(active_normals, d) != d:
            raise AssertionError("incremental enumeration emitted a nonvertex")
        computed.add(vertex.point)
    supplied = _load_supplied(n, d, vertices_dir)
    missing = sorted(computed - supplied, reverse=True)
    extra = sorted(supplied - computed, reverse=True)
    if missing or extra:
        raise AssertionError(
            f"({n},{d}) supplied V-representation mismatch: "
            f"{len(missing)} missing, {len(extra)} extra"
        )
    return initial_count, len(computed)


def _parse_system(value: str) -> tuple[int, int]:
    try:
        n_text, d_text = value.split(",")
        result = (int(n_text), int(d_text))
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError("system must have form N,D") from error
    if result not in SYSTEMS:
        raise argparse.ArgumentTypeError(f"unknown supplied system {result}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", type=_parse_system, help="check only N,D")
    parser.add_argument(
        "--constraints",
        type=Path,
        default=DEFAULT_CONSTRAINTS,
        help="rational constraints.json input",
    )
    parser.add_argument(
        "--vertices",
        type=Path,
        default=DEFAULT_VERTICES,
        help="directory containing vertices_N_D.json files",
    )
    parser.add_argument("--verbose", action="store_true", help="report each H-row update")
    args = parser.parse_args()
    table = json.loads(args.constraints.read_text(encoding="utf-8"))
    selected = (args.system,) if args.system else SYSTEMS
    total = 0
    for n, d in selected:
        initial_count, final_count = verify_system(
            n,
            d,
            table,
            vertices_dir=args.vertices,
            verbose=args.verbose,
        )
        total += final_count
        print(
            f"PASS ({n},{d}): exact incremental H-to-V enumeration "
            f"{initial_count} -> {final_count} vertices"
        )
    print(f"PASS: {total} tabulated vertices exhaust {len(selected)} H-polytopes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
