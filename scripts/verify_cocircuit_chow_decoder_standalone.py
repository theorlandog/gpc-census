#!/usr/bin/env python3
"""Standalone replay of the cocircuit and signed-Chow-decoder claims.

Imports NOTHING from ``gpc_census``. Everything below is rebuilt from
``itertools``, ``fractions`` and integer arithmetic, so agreement with
``results/data/cocircuit_chow_decoder.json`` is evidence about the mathematics
rather than about one implementation agreeing with itself.

Four replays, in increasing order of how much they would hurt if they failed.

1. C1, operationally. Admissible hyperplanes at `(3,5)` and `(3,6)` are
   enumerated by definition: every independent `(d-1)`-subset of the slice, its
   exact cofactor annihilator, and the zero set of that annihilator, which is
   by construction a maximal coplanar subset of rank `d-1`. No reverse search,
   no pruning, no modular arithmetic. The distinct annihilators are the
   hyperplanes, and their trace survivors are recounted the same way.

2. T1, EXHAUSTIVELY, at `(3,6)`. Every one of the 2^20 subfamilies of the slice
   is given its shadow, and every shadow that a threshold family realizes is
   checked to have exactly ONE preimage among all 1,048,576 subsets. This is
   the determination theorem evaluated rather than argued, and it is what makes
   any correct decoder's answer unique at that rank.

3. T4 on random normals: equal coordinates of `c_+` force the corresponding
   orbital transposition to fix `B_+` itself.

4. T5 on random normals: averaging the normal over the equal-degree blocks
   yields a separator that reproduces `B_+` exactly, whose block levels are
   strictly ordered the same way as `c_+`, and whose selected occupancy
   profiles form an upper ideal in prefix-sum dominance.

Then the same replays at arbitrary particle number, because none of the proofs
uses `N = 3`: brute force and exhaustive T1 at `(4,6)` with its `(2,6)` dual,
T4 and T5 at `(4,8)` and `(5,9)`, and the particle-hole shadow identity
`c(F^c) = m - c(F)` checked against the complement family itself. When
`results/data/cocircuit_chow_higher_n.json` is present, its dual-population
rows are cross-read against the `N=3` artifact and its half-filling wall is
checked to publish no partial count.

What this file does NOT replay: the full-population decoding at `(3,7)`,
`(3,8)` and `(5,8)`. Those runs need the decoder itself, and are covered by
`scripts/cocircuit_chow_decoder.py`, `scripts/cocircuit_chow_higher_n.py` and
their guard tests.

Run:

    python scripts/verify_cocircuit_chow_decoder_standalone.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import pathlib
import random
from fractions import Fraction

DEFAULT_ARTIFACT = (
    pathlib.Path(__file__).resolve().parents[1]
    / "results/data/cocircuit_chow_decoder.json"
)
DEFAULT_HIGHER_N_ARTIFACT = (
    pathlib.Path(__file__).resolve().parents[1]
    / "results/data/cocircuit_chow_higher_n.json"
)


def slice_weights(n: int, d: int):
    return list(itertools.combinations(range(d), n))


def primitive(vector):
    divisor = 0
    for value in vector:
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:
        return tuple(vector)
    values = tuple(value // divisor for value in vector)
    lead = next(value for value in values if value)
    return values if lead > 0 else tuple(-value for value in values)


def integer_determinant(matrix) -> int:
    """Exact determinant by fraction-free elimination."""
    size = len(matrix)
    if size == 0:
        return 1
    work = [[Fraction(value) for value in row] for row in matrix]
    determinant = Fraction(1)
    for column in range(size):
        pivot = next((r for r in range(column, size) if work[r][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        determinant *= work[column][column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(column + 1, size):
            if work[row][column]:
                factor = work[row][column]
                work[row] = [
                    left - factor * right
                    for left, right in zip(work[row], work[column])
                ]
    assert determinant.denominator == 1
    return int(determinant)


def cofactor_nullvector(rows):
    """Signed maximal minors of a `(d-1) x d` matrix: its exact annihilator."""
    width = len(rows[0])
    return tuple(
        (-1) ** column
        * integer_determinant([row[:column] + row[column + 1 :] for row in rows])
        for column in range(width)
    )


def brute_force_hyperplanes(n: int, d: int) -> dict:
    """Maximal rank-`d-1` coplanar subsets, straight from the definition.

    Every such subset contains an independent `(d-1)`-subset, whose exact
    annihilator is the normal, and the subset is that normal's zero set. So
    ranging over all `(d-1)`-subsets and keeping the distinct annihilators
    enumerates them all, once each.
    """
    weights = slice_weights(n, d)
    rows = [[1 if index in subset else 0 for index in range(d)] for subset in weights]
    found: dict[tuple[int, ...], tuple[int, ...]] = {}
    for choice in itertools.combinations(range(len(weights)), d - 1):
        tau = cofactor_nullvector([rows[index] for index in choice])
        if not any(tau):
            continue
        tau = primitive(tau)
        if tau in found:
            continue
        found[tau] = tuple(
            index
            for index, subset in enumerate(weights)
            if sum(tau[i] for i in subset) == 0
        )
    return found


def trace_counts(n: int, d: int, normals) -> tuple[int, int]:
    """Ressayre trace survivors and nonstructural survivors, both orientations."""
    weights = slice_weights(n, d)
    survivors = 0
    nonstructural = 0
    for base in normals:
        for sign in (1, -1):
            tau = tuple(sign * value for value in base)
            positive = sum(
                1 for subset in weights if sum(tau[i] for i in subset) > 0
            )
            inversions = sum(
                1
                for i in range(d)
                for j in range(i + 1, d)
                if tau[j] > tau[i]
            )
            if positive != inversions:
                continue
            survivors += 1
            if positive:
                nonstructural += 1
    return survivors, nonstructural


def orbit_key(tau):
    return min(tuple(sorted(tau)), tuple(sorted(-value for value in tau)))


def positive_family(tau, n: int, sign: int = 1) -> frozenset:
    return frozenset(
        subset
        for subset in slice_weights(n, len(tau))
        if sign * sum(tau[index] for index in subset) > 0
    )


def shadow(family, d: int) -> tuple[int, ...]:
    counts = [0] * d
    for subset in family:
        for index in subset:
            counts[index] += 1
    return tuple(counts)


def exhaustive_t1(n: int, d: int) -> int:
    """Every threshold shadow at `(n,d)` has exactly one preimage subfamily.

    Ranges over the entire power set of the slice, which is why it is confined
    to `(3,6)`: 2^20 subfamilies is a second of vectorised work and 2^35 is not
    a computation.
    """
    import numpy as np

    weights = slice_weights(n, d)
    count = len(weights)
    if count > 20:
        raise ValueError("exhaustive T1 is only affordable up to a 20-weight slice")
    incidence = np.array(
        [[1 if index in subset else 0 for index in range(d)] for subset in weights],
        dtype=np.int16,
    )
    universe = np.arange(1 << count, dtype=np.uint32)
    membership = np.unpackbits(
        universe.view(np.uint8).reshape(-1, 4), axis=1, bitorder="little"
    )[:, :count].astype(np.int16)
    shadows = membership @ incidence

    base = int(shadows.max()) + 1
    powers = base ** np.arange(d, dtype=np.int64)
    codes = shadows.astype(np.int64) @ powers
    values, counts = np.unique(codes, return_counts=True)
    multiplicity = dict(zip(values.tolist(), counts.tolist()))

    checked = 0
    for tau, _zero_set in brute_force_hyperplanes(n, d).items():
        for sign in (1, -1):
            oriented = tuple(sign * value for value in tau)
            family = positive_family(oriented, n)
            if not family:
                continue
            code = int(
                sum(
                    value * int(powers[index])
                    for index, value in enumerate(shadow(family, d))
                )
            )
            if multiplicity[code] != 1:
                raise AssertionError(
                    f"T1 failed: shadow of tau={oriented} has "
                    f"{multiplicity[code]} preimages"
                )
            checked += 1
    return checked


def blocks_of(counts):
    levels = sorted(set(counts), reverse=True)
    return [
        [index for index, value in enumerate(counts) if value == level]
        for level in levels
    ], levels


def profile_of(subset, blocks):
    return tuple(len(set(subset) & set(block)) for block in blocks)


def dominates(high, low) -> bool:
    running_high = running_low = 0
    for above, below in zip(high[:-1], low[:-1]):
        running_high += above
        running_low += below
        if running_high < running_low:
            return False
    return True


def check_t4(trials: int, n: int, d: int, seed: int) -> int:
    rng = random.Random(seed)
    checked = 0
    for _ in range(trials):
        tau = tuple(rng.randint(-6, 6) for _ in range(d))
        family = positive_family(tau, n)
        if not family:
            continue
        counts = shadow(family, d)
        for i in range(d):
            for j in range(i + 1, d):
                if counts[i] != counts[j]:
                    continue
                swapped = frozenset(
                    tuple(
                        sorted(
                            j if index == i else i if index == j else index
                            for index in subset
                        )
                    )
                    for subset in family
                )
                if swapped != family:
                    raise AssertionError(f"T4 failed at tau={tau} on ({i},{j})")
                checked += 1
    return checked


def check_t5(trials: int, n: int, d: int, seed: int) -> int:
    rng = random.Random(seed)
    checked = 0
    for _ in range(trials):
        tau = tuple(rng.randint(-6, 6) for _ in range(d))
        family = positive_family(tau, n)
        if not family:
            continue
        counts = shadow(family, d)
        blocks, _levels = blocks_of(counts)

        averaged = [Fraction(0)] * d
        for block in blocks:
            mean = Fraction(sum(tau[index] for index in block), len(block))
            for index in block:
                averaged[index] = mean
        rebuilt = frozenset(
            subset
            for subset in slice_weights(n, d)
            if sum(averaged[index] for index in subset) > 0
        )
        if rebuilt != family:
            raise AssertionError(f"T5 separator failed at tau={tau}")

        block_levels = [averaged[block[0]] for block in blocks]
        if any(
            block_levels[a] <= block_levels[a + 1]
            for a in range(len(block_levels) - 1)
        ):
            raise AssertionError(f"T5 block order failed at tau={tau}")

        selected = {profile_of(subset, blocks) for subset in family}
        rejected = {
            profile_of(subset, blocks) for subset in slice_weights(n, d)
        } - selected
        for chosen in selected:
            for other in rejected:
                if dominates(other, chosen):
                    raise AssertionError(f"T5 upper ideal failed at tau={tau}")
        checked += 1
    return checked


def check_particle_hole(trials: int, n: int, d: int, seed: int) -> int:
    """`c(F^c) = m - c(F)`, checked against the complement family itself."""
    rng = random.Random(seed)
    universe = frozenset(range(d))
    checked = 0
    for _ in range(trials):
        tau = tuple(rng.randint(-6, 6) for _ in range(d))
        family = positive_family(tau, n)
        if not family:
            continue
        counts = shadow(family, d)
        complement = frozenset(
            tuple(sorted(universe - frozenset(subset))) for subset in family
        )
        expected = tuple(len(family) - value for value in counts)
        if shadow(complement, d) != expected:
            raise AssertionError(f"particle-hole shadow failed at tau={tau}")
        checked += 1
    return checked


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=pathlib.Path, default=DEFAULT_ARTIFACT)
    parser.add_argument(
        "--higher-n-artifact",
        type=pathlib.Path,
        default=DEFAULT_HIGHER_N_ARTIFACT,
    )
    arguments = parser.parse_args()
    recorded = json.loads(arguments.artifact.read_text())["cocircuit"]

    for d in (5, 6):
        found = brute_force_hyperplanes(3, d)
        survivors, nonstructural = trace_counts(3, d, found)
        orbits = len({orbit_key(tau) for tau in found})
        entry = recorded[str(d)]
        assert len(found) == entry["admissible_hyperplanes"], (
            f"(3,{d}) brute force {len(found)} against artifact "
            f"{entry['admissible_hyperplanes']}"
        )
        assert survivors == entry["trace_survivors"]
        assert nonstructural == entry["nonstructural_trace_survivors"]
        assert orbits == entry["sd_unoriented_orbits"]
        print(
            f"(3,{d}) brute force: {len(found)} hyperplanes, {orbits} orbits, "
            f"{survivors} trace survivors, {nonstructural} nonstructural"
        )

    checked = exhaustive_t1(3, 6)
    print(f"T1 exhaustive at (3,6): {checked} threshold shadows, each with one "
          "preimage among all 1,048,576 subfamilies")

    for d in (7, 8):
        pairs = check_t4(60, 3, d, seed=40 + d)
        print(f"T4 verified on {pairs} equal-degree transpositions at (3,{d})")
        rows = check_t5(60, 3, d, seed=50 + d)
        print(f"T5 verified on {rows} random normals at (3,{d})")

    # ---------------------------------------------------------- arbitrary N
    higher = None
    if arguments.higher_n_artifact.exists():
        higher = json.loads(arguments.higher_n_artifact.read_text())

    found = brute_force_hyperplanes(4, 6)
    survivors, nonstructural = trace_counts(4, 6, found)
    orbits = len({orbit_key(tau) for tau in found})
    dual = brute_force_hyperplanes(2, 6)
    assert len(dual) == len(found), "particle-hole duals disagree at rank 6"
    print(
        f"(4,6) brute force: {len(found)} hyperplanes, {orbits} orbits, "
        f"{survivors} trace survivors, {nonstructural} nonstructural; "
        f"(2,6) dual agrees at {len(dual)}"
    )

    checked = exhaustive_t1(4, 6)
    print(f"T1 exhaustive at (4,6): {checked} threshold shadows, each with one "
          "preimage among all 32,768 subfamilies")

    for n, d in ((4, 8), (5, 9)):
        pairs = check_t4(40, n, d, seed=60 + d)
        print(f"T4 verified on {pairs} equal-degree transpositions at ({n},{d})")
        rows = check_t5(40, n, d, seed=70 + d)
        print(f"T5 verified on {rows} random normals at ({n},{d})")

    for n, d in ((5, 8), (7, 10)):
        rows = check_particle_hole(40, n, d, seed=80 + d)
        print(f"particle-hole shadow identity verified on {rows} families at ({n},{d})")

    if higher is not None:
        for system, row in sorted(higher.get("dual_populations", {}).items()):
            n, d = (int(part) for part in system.strip("()").split(","))
            partner = recorded.get(str(d))
            if partner is None:
                continue
            assert row["admissible_hyperplanes"] == partner["admissible_hyperplanes"]
            assert row["trace_survivors"] == partner["trace_survivors"]
            print(
                f"{system} matches its dual ({d - n},{d}) at "
                f"{row['admissible_hyperplanes']} hyperplanes and "
                f"{row['trace_survivors']} trace survivors"
            )
        wall = higher.get("half_filling_wall") or {}
        if wall and not wall.get("completed", True):
            assert wall.get("admissible_hyperplanes") is None, (
                "a non-completing reverse search must publish no population"
            )
            print(
                f"{wall['system']} recorded as not completing within "
                f"{wall['node_budget']} nodes, with no partial count"
            )

    print("STANDALONE VERIFICATION PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
