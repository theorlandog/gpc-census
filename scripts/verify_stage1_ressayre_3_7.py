#!/usr/bin/env python3
"""Independent verifier for the stage-1 fixed-N=3 rank-7 Ressayre manifest.

This checks the manifest against mathematics, not against the generator that
wrote it. Nothing here imports the generation package: the closed-flat lattice
is re-enumerated in exact rational arithmetic, the Farkas removals and facet
witnesses are replayed from the recorded multipliers, the Ressayre on/below
splits are recomputed from h and z against the exterior weights, and the final
retained rows are compared with the repository's own stored (3,7) table.

What this cannot check is the part that needs the generation modules: the
Ressayre tangent determinants, the trace-condition rejections, and the claim
that the enumeration was exhaustive in the first place. Those are recorded as
unverified here rather than silently treated as confirmed.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "results/data/stage1_ressayre_3_7.json"
sys.path.insert(0, str(ROOT / "src"))


def exterior_weights(n: int, d: int) -> list[tuple[int, ...]]:
    """Occupation vectors of the wedge^n C^d weight basis, lexicographic."""
    return [
        tuple(1 if i in subset else 0 for i in range(d))
        for subset in itertools.combinations(range(d), n)
    ]


def rref(rows: list[list[Fraction]]) -> tuple[list[tuple[Fraction, ...]], list[int]]:
    """Exact reduced row echelon form over Q."""
    matrix = [list(row) for row in rows]
    pivots: list[int] = []
    pivot_row = 0
    width = len(matrix[0]) if matrix else 0
    for column in range(width):
        selected = next((r for r in range(pivot_row, len(matrix)) if matrix[r][column]), None)
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        inverse = matrix[pivot_row][column]
        matrix[pivot_row] = [value / inverse for value in matrix[pivot_row]]
        for r in range(len(matrix)):
            if r != pivot_row and matrix[r][column]:
                factor = matrix[r][column]
                matrix[r] = [a - factor * b for a, b in zip(matrix[r], matrix[pivot_row])]
        pivots.append(column)
        pivot_row += 1
    return [tuple(row) for row in matrix[:pivot_row]], pivots


def residual(row, basis, pivots):
    out = list(row)
    for basis_row, pivot in zip(basis, pivots):
        factor = out[pivot]
        if factor:
            out = [a - factor * b for a, b in zip(out, basis_row)]
    return out


def closed_flat_counts(n: int, d: int) -> list[int]:
    """Count closed affine flats of the augmented weight configuration.

    Exact over Q. A flat is closed when it contains every weight lying in its
    affine span, and every rank-(k+1) closed flat is the closure of a rank-k
    closed flat plus one projective residual class, which is what makes the
    layer-by-layer construction exhaustive.
    """
    points = [tuple(Fraction(v) for v in w) + (Fraction(1),) for w in exterior_weights(n, d)]
    width = len(points)
    states: dict[int, tuple[int, ...]] = {0: ()}
    counts = [1]
    for _ in range(1, d):
        following: dict[int, tuple[int, ...]] = {}
        for mask, basis_indices in states.items():
            basis, pivots = rref([list(points[i]) for i in basis_indices]) if basis_indices else ([], [])
            classes: dict[tuple[Fraction, ...], int] = {}
            for i in range(width):
                if mask >> i & 1:
                    continue
                res = residual(list(points[i]), basis, pivots)
                lead = next(k for k, v in enumerate(res) if v)
                classes.setdefault(tuple(v / res[lead] for v in res), i)
            for representative in classes.values():
                extended = tuple(sorted(basis_indices + (representative,)))
                new_basis, new_pivots = rref([list(points[i]) for i in extended])
                new_mask = 0
                for i in range(width):
                    if not any(residual(list(points[i]), new_basis, new_pivots)):
                        new_mask |= 1 << i
                following.setdefault(new_mask, extended)
        states = following
        counts.append(len(states))
    return counts


def _fraction(pair) -> Fraction:
    return Fraction(pair[0], pair[1])


def check_hadamard(manifest, report):
    enumeration = manifest["candidate_enumeration"]
    n, d = enumeration["particle_number"], enumeration["rank"]
    bound = (n + 1) ** d
    prime = enumeration["rank_preserving_modulus"]
    report("hadamard_bound_squared", enumeration["hadamard_bound_squared"] == bound)
    report("prime_exceeds_hadamard_bound", prime * prime > bound)


def check_arithmetic(manifest, report):
    enumeration = manifest["candidate_enumeration"]
    partition = manifest["candidate_status_partition"]
    counts = manifest["candidate_system_counts"]
    report(
        "orientations_are_two_per_hyperplane",
        2 * enumeration["admissible_hyperplane_count"] == enumeration["orientation_count"],
    )
    report(
        "status_partition_is_complete",
        sum(v["count"] for v in partition.values()) == enumeration["orientation_count"],
    )
    report(
        "certified_splits_into_redundant_and_retained",
        counts["redundant_rows"] + counts["retained_rows"] == counts["certified_rows"],
    )
    report("no_unresolved_rows", counts["unresolved_rows"] == 0)


def check_full_dimension(manifest, report):
    block = manifest["ambient_full_dimension_certificate"]
    d = block["rank"]
    state = d
    coefficients = []
    for _ in range(len(block["coefficients"])):
        state = (1103515245 * state + 12345) % (2**31)
        coefficients.append(1 + state % 100)
    report("seed_rule_reproduces_coefficients", coefficients == block["coefficients"])
    report("coefficient_count_is_weight_count", len(coefficients) == math.comb(d, 3))
    report("expected_skew_rank_is_binom_d_2", block["expected_skew_rank"] == math.comb(d, 2))
    report("expected_symmetric_rank_is_d_times_d_plus_one_over_two",
           block["expected_symmetric_rank"] == d * (d + 1) // 2)
    report("skew_rank_attains_expectation",
           block["skew_rank_mod_prime"] == block["expected_skew_rank"])
    report("symmetric_rank_attains_expectation",
           block["symmetric_rank_mod_prime"] == block["expected_symmetric_rank"])


def check_farkas(manifest, report):
    reduction = manifest["ordered_slice_reduction_certificate"]
    candidates, structural = reduction["candidate_rows"], reduction["structural_rows"]
    d, n = reduction["rank"], reduction["particle_number"]
    failures = []
    for removal in reduction["removals"]:
        target = candidates[removal["target_index"]]
        coeffs = [Fraction(0)] * d
        rhs = Fraction(0)
        negative = False
        for term in removal["candidate_terms"]:
            mu = _fraction(term["multiplier"])
            negative |= mu < 0
            row = candidates[term["index"]]
            coeffs = [c + mu * a for c, a in zip(coeffs, row["coeffs"])]
            rhs += mu * row["rhs"]
        for term in removal["structural_terms"]:
            mu = _fraction(term["multiplier"])
            negative |= mu < 0
            row = structural[term["index"]]
            coeffs = [c + mu * a for c, a in zip(coeffs, row["coeffs"])]
            rhs += mu * row["rhs"]
        trace = _fraction(removal["trace_multiplier"])
        coeffs = [c + trace for c in coeffs]
        rhs += trace * n
        slack = _fraction(removal["rhs_slack"])
        if negative or coeffs != [Fraction(v) for v in target["coeffs"]] or rhs + slack != target["rhs"]:
            failures.append(removal["target_index"])
    report("every_removal_is_a_valid_farkas_combination", not failures, failures)
    removed = {r["target_index"] for r in reduction["removals"]}
    report("removed_set_matches_redundant_indices", sorted(removed) == reduction["redundant_indices"])
    report(
        "retained_are_exactly_the_unremoved",
        sorted(set(range(len(candidates))) - removed) == reduction["retained_indices"],
    )


def check_facets(manifest, report):
    reduction = manifest["ordered_slice_reduction_certificate"]
    candidates, structural = reduction["candidate_rows"], reduction["structural_rows"]
    n = reduction["particle_number"]
    ok = True
    for facet in reduction["facets"]:
        point = [_fraction(x) for x in facet["violating_point"]]
        target = candidates[facet["target_index"]]
        margin = _fraction(facet["violation_margin"])
        violates = sum(a * b for a, b in zip(target["coeffs"], point)) - target["rhs"] == margin
        others = all(
            sum(a * b for a, b in zip(candidates[i]["coeffs"], point)) <= candidates[i]["rhs"]
            for i in reduction["retained_indices"]
            if i != facet["target_index"]
        )
        slice_ok = all(
            sum(a * b for a, b in zip(row["coeffs"], point)) <= row["rhs"] for row in structural
        )
        ok &= violates and margin > 0 and others and slice_ok and sum(point) == n
    report("every_retained_row_has_a_genuine_facet_witness", ok)
    interior = [_fraction(x) for x in reduction["relative_interior_point"]]
    report("relative_interior_point_has_the_right_trace", sum(interior) == n)
    report(
        "relative_interior_point_is_strictly_feasible",
        all(
            sum(a * b for a, b in zip(candidates[i]["coeffs"], interior)) < candidates[i]["rhs"]
            for i in reduction["retained_indices"]
        ),
    )


def check_ressayre_rows(manifest, report):
    d = manifest["rank"]
    n = manifest["particle_number"]
    weights = exterior_weights(n, d)
    ok_on = ok_tau = True
    for row in manifest["retained_rows"]:
        certificate = row["ressayre_certificate"]
        h, z = certificate["h"], certificate["z"]
        values = [sum(h[i] * w[i] for i in range(d)) for w in weights]
        ok_on &= {k for k, v in enumerate(values) if v == z} == set(certificate["on_indices"])
        ok_on &= set(certificate["below_indices"]) <= {k for k, v in enumerate(values) if v < z}
        raw = [z - n * value for value in h]
        divisor = 0
        for value in raw:
            divisor = math.gcd(divisor, abs(value))
        ok_tau &= (([value // divisor for value in raw] if divisor else raw) == row["tau"])
    report("on_and_below_index_splits_recompute_from_h_and_z", ok_on)
    report("tau_is_the_primitive_form_of_z_minus_n_h", ok_tau)


def check_against_repository_table(manifest, report):
    try:
        from gpc_census.constraints import constraints
    except ImportError:
        report("retained_rows_match_repository_table", None)
        return
    d, n = manifest["rank"], manifest["particle_number"]
    stored = {(tuple(q["coeffs"]), q["rhs"]) for q in constraints(n, d)["inequalities"]}
    generated = {
        (tuple(r["constraint"]["coeffs"]), r["constraint"]["rhs"])
        for r in manifest["retained_rows"]
    }
    report("retained_rows_match_repository_table", stored == generated)


def check_flats(manifest, report):
    enumeration = manifest["candidate_enumeration"]
    counts = closed_flat_counts(enumeration["particle_number"], enumeration["rank"])
    report(
        "closed_flat_counts_reproduce_over_Q",
        counts == enumeration["flat_counts_by_affine_matrix_rank"],
        counts,
    )
    report(
        "top_flat_count_is_the_hyperplane_count",
        counts[-1] == enumeration["admissible_hyperplane_count"],
    )


UNVERIFIED_HERE = (
    "ressayre_tangent_determinants",
    "trace_condition_rejections",
    "exhaustiveness_of_the_candidate_run",
    "cyclic_schubert_coefficients",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--flats",
        action="store_true",
        help="also re-enumerate the closed-flat lattice exactly over Q (slow)",
    )
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())

    results: list[tuple[str, object]] = []

    def report(name, value, detail=None):
        results.append((name, value))
        mark = {True: "PASS", False: "FAIL", None: "SKIP"}[value if value in (True, False, None) else bool(value)]
        line = f"  {mark:4s} {name}"
        if detail is not None and value is not True:
            line += f"  {detail}"
        print(line)

    print("stage-1 (3,7) manifest verification")
    check_arithmetic(manifest, report)
    check_hadamard(manifest, report)
    check_full_dimension(manifest, report)
    check_farkas(manifest, report)
    check_facets(manifest, report)
    check_ressayre_rows(manifest, report)
    check_against_repository_table(manifest, report)
    if args.flats:
        check_flats(manifest, report)
    else:
        print("  SKIP closed_flat_counts_reproduce_over_Q  (pass --flats to run)")

    print("\nnot checkable without the generation package:")
    for item in UNVERIFIED_HERE:
        print(f"  {item}")

    failed = [name for name, value in results if value is False]
    print(f"\n{len(results) - len(failed)} of {len(results)} checks passed")
    if failed:
        print("failed: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
