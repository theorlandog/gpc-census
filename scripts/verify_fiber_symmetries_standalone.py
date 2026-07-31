#!/usr/bin/env python3
"""Standard-library verifier for results/data/fiber_symmetries.json."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"


def determinant_bareiss(matrix: list[list[int]]) -> int:
    """Exact integer determinant by fraction-free Bareiss elimination."""
    n = len(matrix)
    if n == 0:
        return 1
    if any(len(row) != n for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(n - 1):
        pivot_row = next((row for row in range(column, n) if work[row][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign *= -1
        pivot = work[column][column]
        for row in range(column + 1, n):
            for col in range(column + 1, n):
                numerator = pivot * work[row][col] - work[row][column] * work[column][col]
                if numerator % previous:
                    raise ArithmeticError("non-exact Bareiss division")
                work[row][col] = numerator // previous
        for row in range(column + 1, n):
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def rational_rank(matrix: list[list[int]]) -> int:
    """Exact rank over Q."""
    if not matrix:
        return 0
    work = [[Fraction(value) for value in row] for row in matrix]
    rows, columns = len(work), len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][column]
        work[rank] = [entry / value for entry in work[rank]]
        for row in range(rows):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                entry - factor * pivot_entry
                for entry, pivot_entry in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(
    artifact_path: Path = DATA / "fiber_symmetries.json",
    states_path: Path = DATA / "states.jsonl",
    transport_path: Path = DATA / "two_block_family.json",
) -> list[str]:
    errors = []
    artifact = json.loads(artifact_path.read_text())
    states = [
        json.loads(line) for line in states_path.read_text().splitlines() if line.strip()
    ]
    state_by_key = {(record["system"], record["index"]): record for record in states}
    transport_rows = json.loads(transport_path.read_text())["states"]
    transport_by_key = {
        (record["system"], record["index"]): record["transport_class"]
        for record in transport_rows
    }

    if artifact.get("states_sha256") != _sha256(states_path):
        errors.append("states SHA-256 mismatch")
    if artifact.get("transport_metadata_sha256") != _sha256(transport_path):
        errors.append("transport metadata SHA-256 mismatch")

    certificates = artifact.get("records", [])
    certificate_keys = set()
    kernel_dimensions = {}
    for certificate in certificates:
        key = (certificate.get("system"), certificate.get("index"))
        if key in certificate_keys:
            errors.append(f"duplicate certificate {key}")
            continue
        certificate_keys.add(key)
        state = state_by_key.get(key)
        if state is None:
            errors.append(f"certificate has no state {key}")
            continue
        support = [tuple(det) for det in state["closed_form"]["support_dets"]]
        d = len(state["integer_form"])
        pairs = [tuple(pair) for pair in certificate.get("minor_pairs", [])]
        if len(pairs) != len(support) or len(set(pairs)) != len(pairs):
            errors.append(f"{key}: minor row count or uniqueness failure")
            continue
        valid_pairs = set(itertools.combinations(range(d), 2))
        if not set(pairs).issubset(valid_pairs):
            errors.append(f"{key}: invalid orbital pair")
            continue
        minor = [
            [int(i in determinant and j in determinant) for determinant in support]
            for i, j in pairs
        ]
        determinant = determinant_bareiss(minor)
        if determinant == 0:
            errors.append(f"{key}: zero pair-incidence minor")
        if determinant != certificate.get("minor_determinant"):
            errors.append(f"{key}: recorded pair-incidence determinant mismatch")
        if certificate.get("support_size") != len(support):
            errors.append(f"{key}: support size mismatch")
        if certificate.get("classified") != state["classified"]:
            errors.append(f"{key}: classification mismatch")

        incidence = [
            [int(orbital in determinant) for determinant in support]
            for orbital in range(d)
        ]
        kernel_dimensions[key] = len(support) - rational_rank(incidence)

    if certificate_keys != set(state_by_key):
        errors.append("certificate target set does not equal the state atlas")

    summary = artifact.get("summary", {})
    positive = [key for key, dimension in kernel_dimensions.items() if dimension > 0]
    positive_classes = {transport_by_key[key] for key in positive}
    expected = {
        "records": len(states),
        "pair_minor_certificates": len(certificates),
        "full_column_rank": len(certificates),
        "rank_failures": 0,
        "transport_classes": len(set(transport_by_key.values())),
        "positive_incidence_kernel_records": len(positive),
        "positive_incidence_kernel_transport_classes": len(positive_classes),
        "positive_kernel_pair_map_failures": 0,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            errors.append(f"summary {key}: expected {value}, got {summary.get(key)}")

    class_counts = {}
    for state in states:
        class_counts[state["classified"]] = class_counts.get(state["classified"], 0) + 1
    if summary.get("classification_counts") != dict(sorted(class_counts.items())):
        errors.append("classification-count summary mismatch")

    kernel_counts = {}
    for dimension in kernel_dimensions.values():
        key = str(dimension)
        kernel_counts[key] = kernel_counts.get(key, 0) + 1
    if summary.get("incidence_kernel_dimension_counts") != dict(
        sorted(kernel_counts.items(), key=lambda item: int(item[0]))
    ):
        errors.append("incidence-kernel summary mismatch")

    vb = artifact.get("vb", {})
    if vb.get("support_preserving_block_permutation_automorphisms") != 1:
        errors.append("v_B support-permutation automorphism count mismatch")
    if vb.get("two_rdm_moments", {}).get("trace_3_derivative") != "-48/12167":
        errors.append("v_B cubic-moment derivative mismatch")
    if vb.get("galois_walls", {}).get("difference") != "1728*sqrt(35)/1034195":
        errors.append("v_B Galois-wall invariant difference mismatch")
    tangent = vb.get("tangent_moment_rank", {})
    if tangent.get("fixed_spectrum_constraint_rank") != 7:
        errors.append("v_B exact fixed-spectrum constraint rank mismatch")
    if tangent.get("orbital_phase_gauge_rank") != 7:
        errors.append("v_B exact gauge rank mismatch")
    if tangent.get("quotient_tangent_dimension") != 2:
        errors.append("v_B exact quotient tangent dimension mismatch")
    if tangent.get("combined_rank") != 9:
        errors.append("v_B cubic-quartic combined rank mismatch")
    if tangent.get("combined_minor_determinant") != (
        "299925504*sqrt(8211)/1035662780341225"
    ):
        errors.append("v_B cubic-quartic exact minor mismatch")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact", type=Path, default=DATA / "fiber_symmetries.json"
    )
    parser.add_argument("--states", type=Path, default=DATA / "states.jsonl")
    parser.add_argument(
        "--transport", type=Path, default=DATA / "two_block_family.json"
    )
    args = parser.parse_args()
    errors = verify(args.artifact, args.states, args.transport)
    if errors:
        for error in errors:
            print(f"FAIL: {error}")
        return 1
    print("fiber symmetry certificates verified exactly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
