#!/usr/bin/env python3
"""Replay the four (3,11) level-five Ressayre certificates independently.

Standard library only. This script imports nothing from ``gpc_census`` and
shares no code with the generation package that produced the certificates. It
rebuilds every object from the recorded ``(h, z)`` and checks it.

The repository's validation law is that a claim is not a certificate. The
rank-7 checkpoint already notes that the Ressayre tangent determinants sit
outside its independent verifier's trust boundary, so the rank-11 result would
inherit that gap without this file. What is replayed here:

- the 165 exterior weights of ``wedge^3 C^11``, generated from scratch;
- that the recorded centered pair ``(h, z)`` is the correct transform of the
  recorded constraint, checked both by the closed-form rule and by exact
  rational agreement of the two inequalities on random trace-three points;
- the on and below index sets, recomputed from ``(h, z)`` and the weights;
- the negative root set, recomputed from ``h``;
- the trace condition, that the below count equals the negative root count;
- affine independence of the recorded admissibility witness, by exact rational
  elimination;
- the tangent matrix, rebuilt with a fermionic sign derived by counting
  inversions rather than by the production package's position bookkeeping;
- the determinant, recomputed by exact fraction-free Gaussian elimination over
  ``Fraction``, which is a different algorithm from the Bareiss elimination the
  generator used, and compared to the recorded integer;
- that the certified inequality is satisfied by every spectrum the census has
  certified TRUE at (3,11), and violated by each candidate it refutes.

Exit status is 0 when every check passes and 1 otherwise.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT = ROOT / "results/data/level5_ressayre_3_11.json"

N = 3
D = 11


def weights_of_exterior_power(n: int, d: int) -> list[tuple[int, ...]]:
    """Determinant weights of wedge^n C^d as zero-one occupation vectors."""
    return [
        tuple(1 if index in occupied else 0 for index in range(d))
        for occupied in itertools.combinations(range(d), n)
    ]


def root_action(weight: tuple[int, ...], i: int, j: int):
    """Apply E_(j,i) = a_j^dagger a_i, returning (target, sign) or None.

    The sign is obtained by substituting ``j`` for ``i`` in the ordered
    occupation tuple and counting the inversions that resorting removes. That
    is a different derivation from tracking annihilation and creation
    positions separately, and it agrees only if both are right.
    """
    if not weight[i] or weight[j]:
        return None
    occupied = [index for index, value in enumerate(weight) if value]
    substituted = [j if index == i else index for index in occupied]
    inversions = sum(
        1
        for a in range(len(substituted))
        for b in range(a + 1, len(substituted))
        if substituted[a] > substituted[b]
    )
    sign = -1 if inversions % 2 else 1
    target = tuple(1 if index in substituted else 0 for index in range(len(weight)))
    return target, sign


def exact_determinant(matrix: list[list[int]]) -> int:
    """Exact determinant by Gaussian elimination over Fraction.

    Deliberately not the Bareiss routine the generator uses. Agreement of two
    different exact algorithms on the same matrix is the point.
    """
    rows = [[Fraction(value) for value in row] for row in matrix]
    size = len(rows)
    result = Fraction(1)
    for column in range(size):
        pivot = next((r for r in range(column, size) if rows[r][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            result = -result
        result *= rows[column][column]
        inverse = 1 / rows[column][column]
        rows[column] = [value * inverse for value in rows[column]]
        for r in range(column + 1, size):
            factor = rows[r][column]
            if factor:
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[column])]
    assert result.denominator == 1
    return int(result)


def affinely_independent(points: list[tuple[int, ...]]) -> bool:
    """Exact rational test that the points span an affine space of full rank."""
    if not points:
        return False
    base = points[0]
    rows = [
        [Fraction(a - b) for a, b in zip(point, base)] for point in points[1:]
    ]
    rank, columns = 0, len(base)
    for column in range(columns):
        pivot = next((r for r in range(rank, len(rows)) if rows[r][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        lead = rows[rank][column]
        rows[rank] = [value / lead for value in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][column]:
                factor = rows[r][column]
                rows[r] = [a - factor * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
    return rank == len(points) - 1


def primitive(values):
    divisor = 0
    for value in values:
        divisor = math.gcd(divisor, abs(value))
    return [value // divisor for value in values] if divisor else list(values)


class Checks:
    def __init__(self) -> None:
        self.passed = 0
        self.failed: list[str] = []

    def check(self, label: str, condition: bool) -> None:
        if condition:
            self.passed += 1
        else:
            self.failed.append(label)
        print(f"  [{'ok' if condition else 'FAIL'}] {label}")


def verify_certificate(row: dict, checks: Checks, rng: random.Random) -> None:
    family = row["family"]
    print(f"\n{family}: sum lambda_{row['subset']} <= {row['rhs']}")
    certificate = row["certificate"]
    h = certificate["h"]
    z = certificate["z"]
    coeffs = row["constraint"]["coeffs"]
    rhs = row["constraint"]["rhs"]

    # the constraint the artifact records is the row we claim to have proved
    expected = [1 if i + 1 in row["subset"] else 0 for i in range(D)]
    checks.check(f"{family}: recorded constraint is the claimed row",
                 coeffs == expected and rhs == row["rhs"])

    # tau_i = N * c_i - k, and (h, z) is its centered primitive transform
    tau = [N * c - rhs for c in coeffs]
    checks.check(f"{family}: recorded tau matches N*c - k", list(row["tau"]) == tau)
    total = sum(tau)
    centered = primitive([total - D * value for value in tau] + [N * total])
    checks.check(f"{family}: (h,z) is the centered transform of tau",
                 centered[:-1] == list(h) and centered[-1] == z)

    # the two inequalities must agree on the trace-N slice, checked exactly
    agree = True
    for _trial in range(500):
        parts = [Fraction(rng.randint(0, 40)) for _ in range(D)]
        weight = sum(parts)
        if not weight:
            continue
        lam = [N * part / weight for part in parts]
        left = sum(Fraction(c) * value for c, value in zip(coeffs, lam)) <= rhs
        right = sum(Fraction(c) * value for c, value in zip(h, lam)) >= z
        agree &= left == right
    checks.check(f"{family}: h.lambda >= z is the row on the trace-{N} slice", agree)

    weights = weights_of_exterior_power(N, D)
    checks.check(f"{family}: wedge^{N} C^{D} has {math.comb(D, N)} weights",
                 len(weights) == math.comb(D, N))

    on, below = [], []
    for index, weight in enumerate(weights):
        value = sum(a * b for a, b in zip(weight, h))
        if value == z:
            on.append(index)
        elif value < z:
            below.append(index)
    checks.check(f"{family}: on-set recomputed from (h,z)",
                 on == list(certificate["on_indices"]))
    checks.check(f"{family}: below-set recomputed from (h,z)",
                 below == list(certificate["below_indices"]))

    roots = [
        [i, j]
        for i in range(D)
        for j in range(i + 1, D)
        if h[j] - h[i] < 0
    ]
    checks.check(f"{family}: negative root set recomputed from h",
                 roots == [list(r) for r in certificate["negative_roots"]])
    checks.check(f"{family}: trace condition, {len(below)} below = {len(roots)} roots",
                 len(below) == len(roots))

    witness = list(certificate["admissibility_witness"])
    checks.check(f"{family}: witness has d-1 = {D - 1} weights on the hyperplane",
                 len(witness) == D - 1 and set(witness) <= set(on))
    checks.check(f"{family}: witness is affinely independent",
                 affinely_independent([weights[index] for index in witness]))

    evaluation = list(certificate["determinant_evaluation"])
    checks.check(f"{family}: one evaluation coordinate per on-weight",
                 len(evaluation) == len(on))
    below_row = {weights[index]: r for r, index in enumerate(below)}
    matrix = [[0] * len(roots) for _ in below]
    for column, (i, j) in enumerate(roots):
        for variable, weight_index in enumerate(on):
            action = root_action(weights[weight_index], i, j)
            if action is None or action[0] not in below_row:
                continue
            matrix[below_row[action[0]]][column] += action[1] * evaluation[variable]
    value = exact_determinant(matrix)
    checks.check(f"{family}: determinant replays to {certificate['determinant_value']}",
                 value == certificate["determinant_value"])
    checks.check(f"{family}: determinant is nonzero, so the row is Ressayre", value != 0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    payload = json.loads(args.artifact.read_text())
    eleven = payload["rank_eleven"]
    checks = Checks()
    rng = random.Random(20260804)

    print(f"replaying {eleven['certified_rows']} certificate(s) from {args.artifact}")
    certified = [row for row in eleven["rows"] if row["screen_status"] == "certified"]
    checks.check("all four claimed rows are certified", len(certified) == 4)
    for row in certified:
        verify_certificate(row, checks, rng)

    print("\nconsequences")
    subsets = [(row["family"], row["subset"]) for row in certified]
    attained = [
        [Fraction(x) for x in record["spectrum"]]
        for record in json.loads((ROOT / "docs/bracket_3_11.json").read_text())[
            "outer_vertices"
        ]
    ]
    settlement = json.loads((ROOT / "docs/bracket_3_11_settlement.json").read_text())
    true_indices = [
        c["index"] for c in settlement["candidates"] if c.get("verdict") == "TRUE-VERTEX"
    ]
    checks.check(
        "no certified row is violated by a TRUE (3,11) vertex",
        all(
            sum(attained[index][j - 1] for j in subset) <= 2
            for index in true_indices
            for _family, subset in subsets
        ),
    )
    checks.check(
        "the Slater point satisfies every certified row",
        all(
            sum(([Fraction(1)] * 3 + [Fraction(0)] * 8)[j - 1] for j in subset) <= 2
            for _family, subset in subsets
        ),
    )
    uniform = [Fraction(3, 11)] * 11
    checks.check(
        "the uniform (3/11)^11 point satisfies every certified row",
        all(
            sum(uniform[j - 1] for j in subset) <= 2 for _family, subset in subsets
        ),
    )
    for record in eleven["open_candidates"]:
        lam = [Fraction(x) for x in record["spectrum"]]
        cutting = [
            family for family, subset in subsets if sum(lam[j - 1] for j in subset) > 2
        ]
        checks.check(
            f"candidate {record['index']} is cut off by {cutting or 'nothing'}",
            cutting == [c["family"] for c in record["cut_by"]] and bool(cutting),
        )

    print(f"\n{checks.passed} checks passed, {len(checks.failed)} failed")
    for label in checks.failed:
        print(f"  FAILED: {label}")
    return 0 if not checks.failed else 1


if __name__ == "__main__":
    sys.exit(main())
