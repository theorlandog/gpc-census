#!/usr/bin/env python3
"""Two-sided closure certificate: dimension, facetness, slack rank, duality.

The inner/outer sandwich is already proved in ``docs/census_inner_sandwich.md``:
every stored row replays as an exact Ressayre certificate, and every vertex of
the outer polytope is the spectrum of a certified state, so ``P = O`` at nine
systems without candidate exhaustiveness. That is the completeness half.

This script supplies the independent verification layer that sits on top of it,
in exact rational arithmetic and with nothing hardcoded:

1. AFFINE HULL, CERTIFIED. The dimension ``m`` is computed as the rank of the
   vertex difference matrix, never assumed. This matters: the trace slice at
   ``(3,6)`` has dimension 5 while the polytope has dimension 3, so two further
   affine relations hold and any check that assumed ``m = d-1`` would be wrong
   at the first system.

2. FACETNESS, NOT VALIDITY. A valid row need not be a facet. Each candidate row
   is kept only when its tight vertex set has affine dimension exactly
   ``m - 1``, which rules out redundancy, hidden equalities and lower faces at
   once. Candidate rows are the published GPC inequalities together with the
   chamber walls and the nonnegativity wall.

3. VERTEX RANK. Each listed vertex must have ``rank([E; A_active]) = d``, with
   ``E`` the affine-hull equations. That is what makes it a vertex rather than
   a relative interior point of a face.

4. SLACK MATRIX. ``S_{fv} = b_f - a_f . v`` must be nonnegative, must vanish
   exactly on the incidence relation, and must have rank ``m + 1``.

5. PARTICLE-HOLE PARITY, at half filling only. ``Theta(lambda)_i =
   1 - lambda_{d+1-i}`` maps ``(N,d)`` to ``(d-N,d)``, so it is an involution of
   a single system only when ``2N = d``. Among the census systems that is
   ``(3,6)``, ``(4,8)`` and ``(5,10)``; everywhere else it is a map between two
   different systems and no within-system parity relation is available.

Polar duality is not built explicitly. ``f_{m-1}(P) = f_0(P^o)`` and
``f_0(P) = f_{m-1}(P^o)`` say exactly that no facet is redundant and no listed
vertex fails to be a vertex, and checks 2 and 3 establish both directly, so
constructing the polar would re-derive them rather than test them.

Run:

    uv run scripts/closure_certificate.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys
import time
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SYSTEMS = ((3, 6), (3, 7), (3, 8), (4, 8), (3, 9), (4, 9), (3, 10), (4, 10), (5, 10))
TIMING_FIELDS = ("seconds",)


def row_reduce(rows: list[list[Fraction]]) -> list[list[Fraction]]:
    """Exact reduced row echelon form; returns the nonzero rows."""
    matrix = [list(row) for row in rows]
    pivot_row = 0
    columns = len(matrix[0]) if matrix else 0
    for column in range(columns):
        pivot = None
        for index in range(pivot_row, len(matrix)):
            if matrix[index][column]:
                pivot = index
                break
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for index in range(len(matrix)):
            if index != pivot_row and matrix[index][column]:
                factor = matrix[index][column]
                matrix[index] = [
                    a - factor * b for a, b in zip(matrix[index], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return [row for row in matrix if any(row)]


def rank(rows) -> int:
    rows = [list(row) for row in rows]
    if not rows:
        return 0
    return len(row_reduce(rows))


def affine_dimension(points) -> int:
    points = list(points)
    if len(points) <= 1:
        return 0
    base = points[0]
    return rank([[a - b for a, b in zip(point, base)] for point in points[1:]])


def affine_hull_equations(points) -> list[list[Fraction]]:
    """Rows ``(c, e)`` with ``c . lambda = e`` satisfied by every point."""
    points = list(points)
    dimension = len(points[0])
    base = points[0]
    directions = row_reduce(
        [[a - b for a, b in zip(point, base)] for point in points[1:]]
    )
    pivots = []
    for row in directions:
        for column, value in enumerate(row):
            if value:
                pivots.append(column)
                break
    free = [column for column in range(dimension) if column not in pivots]
    equations = []
    for index in free:
        normal = [Fraction(0)] * dimension
        normal[index] = Fraction(1)
        for row, pivot in zip(directions, pivots):
            normal[pivot] = -row[index]
        offset = sum(a * b for a, b in zip(normal, base))
        equations.append(normal + [offset])
    return equations


def load_vertices(n: int, d: int) -> list[list[Fraction]]:
    path = ROOT / "results" / "data" / "vertices" / f"vertices_{n}_{d}.json"
    records = json.loads(path.read_text())
    vertices = []
    for record in records:
        denominator = record["denominator"]
        vertices.append(
            [Fraction(value, denominator) for value in record["integer_form"]]
        )
    return vertices


def candidate_rows(n: int, d: int):
    """Every candidate facet: GPC rows, chamber walls, the nonnegativity wall."""
    from gpc_census.constraints import constraints

    rows = []
    for inequality in constraints(n, d)["inequalities"]:
        rows.append(
            {
                "kind": "gpc",
                "a": [Fraction(value) for value in inequality["coeffs"]],
                "b": Fraction(inequality["rhs"]),
            }
        )
    for i in range(d - 1):
        coefficients = [Fraction(0)] * d
        coefficients[i] = Fraction(-1)
        coefficients[i + 1] = Fraction(1)
        rows.append({"kind": "chamber", "a": coefficients, "b": Fraction(0)})
    coefficients = [Fraction(0)] * d
    coefficients[d - 1] = Fraction(-1)
    rows.append({"kind": "nonnegativity", "a": coefficients, "b": Fraction(0)})
    return rows


def particle_hole(vertex):
    return [1 - value for value in reversed(vertex)]


def audit_system(n: int, d: int) -> dict:
    started = time.perf_counter()
    vertices = load_vertices(n, d)
    dimension = affine_dimension(vertices)
    equations = affine_hull_equations(vertices)

    rows = candidate_rows(n, d)
    valid = 0
    # A facet is a maximal proper FACE, not a row. Two distinct rows can differ
    # by a relation of the affine hull and cut out the same face, so rows are
    # canonicalized by their tight vertex set before being counted. This is not
    # hypothetical: at (3,6) the hull carries two relations beyond the trace,
    # and six rows pass the codimension test while cutting only four faces.
    codimension_one_rows = 0
    by_face: dict[frozenset, dict] = {}
    for row in rows:
        slacks = [
            row["b"] - sum(a * v for a, v in zip(row["a"], vertex))
            for vertex in vertices
        ]
        if all(slack >= 0 for slack in slacks):
            valid += 1
        else:
            continue
        tight = frozenset(
            index for index, slack in enumerate(slacks) if slack == 0
        )
        if not tight:
            continue
        if affine_dimension([vertices[index] for index in tight]) == dimension - 1:
            codimension_one_rows += 1
            by_face.setdefault(tight, row)
    facets = list(by_face.values())

    # Vertex rank: rank([E; active normals]) must be the ambient d.
    true_vertices = 0
    facets_per_vertex = []
    for vertex in vertices:
        active = [
            facet["a"]
            for facet in facets
            if facet["b"] == sum(a * v for a, v in zip(facet["a"], vertex))
        ]
        facets_per_vertex.append(len(active))
        stacked = [equation[:-1] for equation in equations] + [
            list(normal) for normal in active
        ]
        if rank(stacked) == d:
            true_vertices += 1

    # Slack matrix.
    slack_rows = []
    incidence_consistent = True
    nonnegative = True
    for facet in facets:
        row = []
        for vertex in vertices:
            value = facet["b"] - sum(a * v for a, v in zip(facet["a"], vertex))
            if value < 0:
                nonnegative = False
            row.append(value)
        slack_rows.append(row)
    slack_rank = rank(slack_rows) if slack_rows else 0

    # The zero pattern must be exactly the incidence relation, by construction
    # of the slacks; this re-derives it independently as a guard against a
    # sign or indexing slip.
    for facet, row in zip(facets, slack_rows):
        for vertex, value in zip(vertices, row):
            tight = facet["b"] == sum(a * v for a, v in zip(facet["a"], vertex))
            if (value == 0) != tight:
                incidence_consistent = False

    record = {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "certified_affine_dimension": dimension,
        "affine_hull_equation_count": len(equations),
        "vertices": len(vertices),
        "candidate_rows": len(rows),
        "valid_rows": valid,
        "codimension_one_rows": codimension_one_rows,
        "facets": len(facets),
        "rows_collapsed_by_canonicalization": codimension_one_rows - len(facets),
        "facets_by_kind": {
            kind: sum(1 for facet in facets if facet["kind"] == kind)
            for kind in ("gpc", "chamber", "nonnegativity")
        },
        "vertices_with_full_active_rank": true_vertices,
        "min_facets_per_vertex": min(facets_per_vertex) if facets_per_vertex else 0,
        "max_facets_per_vertex": max(facets_per_vertex) if facets_per_vertex else 0,
        "slack_matrix_nonnegative": nonnegative,
        "slack_zero_pattern_is_incidence": incidence_consistent,
        "slack_rank": slack_rank,
        "slack_rank_equals_dimension_plus_one": slack_rank == dimension + 1,
    }

    if 2 * n == d:
        vertex_set = {tuple(vertex) for vertex in vertices}
        images = [particle_hole(vertex) for vertex in vertices]
        closed = all(tuple(image) in vertex_set for image in images)
        fixed = sum(
            1 for vertex, image in zip(vertices, images) if tuple(vertex) == tuple(image)
        )
        record["particle_hole"] = {
            "applies": True,
            "vertex_set_is_closed": closed,
            "fixed_vertices": fixed,
            "paired_vertex_orbits": (len(vertices) - fixed) // 2,
            "parity_holds": closed
            and fixed + 2 * ((len(vertices) - fixed) // 2) == len(vertices),
        }
    else:
        record["particle_hole"] = {
            "applies": False,
            "reason": (
                "Theta maps (N,d) to (d-N,d), so it is an involution of one "
                "system only at half filling, 2N = d."
            ),
        }

    record["seconds"] = round(time.perf_counter() - started, 1)
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--out",
        type=pathlib.Path,
        default=ROOT / "results/data/closure_certificate.json",
    )
    arguments = parser.parse_args()

    published = json.loads((ROOT / "docs" / "polytope_invariants.json").read_text())
    systems = {}
    mismatches = []
    for n, d in SYSTEMS:
        record = audit_system(n, d)
        reference = published["systems"].get(record["system"])
        if reference:
            record["published_facets"] = reference["facets"]
            record["published_vertices"] = reference["vertices"]
            record["published_dimension"] = reference["dim"]
            agree = (
                record["facets"] == reference["facets"]
                and record["vertices"] == reference["vertices"]
                and record["certified_affine_dimension"] == reference["dim"]
                and record["min_facets_per_vertex"] == reference["min_facets_per_vertex"]
                and record["max_facets_per_vertex"] == reference["max_facets_per_vertex"]
            )
            record["agrees_with_polytope_invariants"] = agree
            if not agree:
                mismatches.append(record["system"])
        systems[record["system"]] = record
        print(
            f"  {record['system']:8s} dim={record['certified_affine_dimension']} "
            f"facets={record['facets']} vertices={record['vertices']} "
            f"slack_rank={record['slack_rank']} "
            f"vertex_rank_ok={record['vertices_with_full_active_rank']}"
            f"/{record['vertices']} "
            f"agree={record.get('agrees_with_polytope_invariants')}",
            flush=True,
        )

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/closure_certificate.py",
        "status": "exact_verification_layer",
        "scope": (
            "Verification layer over the inner/outer sandwich of "
            "docs/census_inner_sandwich.md, which supplies completeness. This "
            "script does not re-prove P = O; it checks dimension, facetness, "
            "vertex rank, slack rank and half-filling parity exactly."
        ),
        "systems": systems,
        "summary": {
            "systems_audited": len(systems),
            "mismatches_against_polytope_invariants": mismatches,
            "all_slack_ranks_correct": all(
                record["slack_rank_equals_dimension_plus_one"]
                for record in systems.values()
            ),
            "all_vertices_have_full_active_rank": all(
                record["vertices_with_full_active_rank"] == record["vertices"]
                for record in systems.values()
            ),
            "half_filling_systems": [
                system
                for system, record in systems.items()
                if record["particle_hole"]["applies"]
            ],
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
