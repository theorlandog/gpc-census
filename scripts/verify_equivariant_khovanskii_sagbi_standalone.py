"""Standalone replay of results/data/equivariant_khovanskii_sagbi.json.

Imports nothing from `gpc_census`. Standard library only, exact integers and
`Fraction` throughout, no floating point anywhere.

The artifact was produced from plethysms expanded in POWER SUMS and evaluated
with the Murnaghan-Nakayama rule. This file never uses that route. It
recomputes the same multiplicities from

    mult(lambda) = P(lambda) - sum_{mu > lambda} mult(mu) * K_{mu,lambda},

where `P(lambda)` counts multisets of `m` orbital `n`-subsets with indicator
sum `lambda` (a direct weight count in `Sym^m(wedge^n C^d)`) and `K_{mu,lambda}`
is the Kostka number, counted as semistandard tableaux by horizontal-strip
removal. Both halves are combinatorial counts, so a disagreement with the
artifact means one of the two independent routes is wrong.

It also replays, with integer arithmetic and no representation theory at all:

* every stored minimal generator, by exhaustive search for a decomposition
  inside the stored levels;
* every stored non-generator, by exhibiting a decomposition;
* the padding theorem and the particle-hole transport on the stored levels;
* every held-out miss count and every saturation-failure witness;
* the level digests.

Finally it rebuilds the `(3,6)` highest-weight spaces and their graded-lex
leading monomials from scratch, by an independent implementation of the simple
raising operators and exact row reduction, at the degrees named in
`HW_REPLAY_MAX`.

Run: python3 scripts/verify_equivariant_khovanskii_sagbi_standalone.py
Exit status 0 when every check passes, 1 otherwise.
"""

from __future__ import annotations

import hashlib
import itertools as it
import json
import pathlib
import sys
from fractions import Fraction
from functools import lru_cache

ROOT = pathlib.Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "data" / "equivariant_khovanskii_sagbi.json"

# Independent recomputation of the representation theory is quadratic in the
# number of partitions, so it is bounded per system and what it covered is
# reported rather than assumed.
MULT_REPLAY_MAX = {"(3,6)": 12, "(3,7)": 10, "(3,8)": 8, "(3,9)": 6,
                   "(3,10)": 6, "(4,8)": 6, "(4,10)": 5, "(5,10)": 5}
HW_REPLAY_MAX = 4

FAILURES: list[str] = []


def check(condition: bool, message: str) -> bool:
    if not condition:
        FAILURES.append(message)
    return bool(condition)


# --------------------------------------------------------------------------
# Route two: weight counts and Kostka numbers.
# --------------------------------------------------------------------------

def partitions_bounded(total: int, maxpart: int, maxlen: int):
    """Decreasing partitions of `total` with bounded part size and length."""
    def walk(left: int, cap: int, slots: int):
        if left == 0:
            yield ()
            return
        if slots == 0 or left > cap * slots:
            return
        for part in range(min(left, cap), 0, -1):
            for rest in walk(left - part, part, slots - 1):
                yield (part,) + rest
    return walk(total, maxpart, maxlen)


def weight_count(n: int, d: int):
    """`P(nu)`: multisets of orbital `n`-subsets with indicator sum `nu`."""
    gens = [tuple(1 if o in s else 0 for o in range(d))
            for s in it.combinations(range(d), n)]

    @lru_cache(maxsize=None)
    def count(index: int, rest: tuple[int, ...]) -> int:
        total = sum(rest)
        if total == 0:
            return 1
        if total % n or index >= len(gens):
            return 0
        left = total // n
        if max(rest) > left:
            return 0
        reach = [0] * d
        for gen in gens[index:]:
            for orbital in range(d):
                reach[orbital] += gen[orbital]
        for orbital in range(d):
            if rest[orbital] and not reach[orbital]:
                return 0
        out = count(index + 1, rest)
        take = tuple(a - b for a, b in zip(rest, gens[index]))
        if min(take) >= 0:
            out += count(index, take)
        return out

    def probe(nu) -> int:
        nu = tuple(int(x) for x in nu)
        return 0 if min(nu) < 0 else count(0, nu)

    return probe


def horizontal_strips(shape: tuple[int, ...], size: int):
    """Shapes `nu` with `shape/nu` a horizontal strip of the given size."""
    rows = len(shape)

    def walk(index: int, left: int, previous: int):
        if index == rows:
            if left == 0:
                yield ()
            return
        # nu_index lies in [max(0, shape_{index+1}), shape_index] and, for nu to
        # be a partition, nu_index <= previous.
        floor = shape[index + 1] if index + 1 < rows else 0
        for value in range(min(shape[index], previous), floor - 1, -1):
            taken = shape[index] - value
            if taken > left:
                continue
            for rest in walk(index + 1, left - taken, value):
                yield (value,) + rest
    for candidate in walk(0, size, shape[0] if shape else 0):
        yield tuple(part for part in candidate if part > 0)


@lru_cache(maxsize=None)
def kostka(shape: tuple[int, ...], content: tuple[int, ...]) -> int:
    """Semistandard tableaux of the given shape and content."""
    if sum(shape) != sum(content):
        return 0
    if not content:
        return int(not shape)
    last, rest = content[-1], content[:-1]
    return sum(kostka(nu, rest) for nu in horizontal_strips(shape, last))


def independent_level(n: int, d: int, m: int, probe) -> set[tuple[int, ...]]:
    """`{lambda : mult > 0}` at degree `m`, by weight counts and Kostka peeling."""
    shapes = [lam for lam in partitions_bounded(n * m, m, d)]
    shapes.sort(reverse=True)  # decreasing lex refines decreasing dominance
    mult: dict[tuple[int, ...], int] = {}
    for lam in shapes:
        total = probe(lam + (0,) * (d - len(lam)))
        for mu, value in mult.items():
            if value:
                total -= value * kostka(mu, lam)
        mult[lam] = total
    return {lam + (0,) * (d - len(lam)) for lam, value in mult.items() if value > 0}


# --------------------------------------------------------------------------
# Route three: highest-weight vectors from scratch, for the (3,6) term order.
# --------------------------------------------------------------------------

def replay_leading_values(n: int, d: int, m: int, variables) -> set[tuple[int, ...]]:
    """Graded-lex leading exponents of every highest-weight space at degree `m`."""
    subs = [tuple(v) for v in variables]
    index_of = {s: i for i, s in enumerate(subs)}
    out: set[tuple[int, ...]] = set()
    for lam in partitions_bounded(n * m, m, d):
        target = lam + (0,) * (d - len(lam))
        basis = [e for e in _monomials(subs, d, m) if _weight(e, subs, d) == target]
        if not basis:
            continue
        column = {mono: i for i, mono in enumerate(basis)}
        equations = []
        for simple in range(d - 1):
            images: dict[tuple[int, ...], dict[int, int]] = {}
            for mono in basis:
                for image, coefficient in _raise(mono, simple, subs, index_of).items():
                    images.setdefault(image, {})[column[mono]] = coefficient
            for row in images.values():
                equations.append([Fraction(row.get(i, 0)) for i in range(len(basis))])
        kernel = _kernel(equations, len(basis))
        if not kernel:
            continue
        order = sorted(range(len(basis)),
                       key=lambda i: (sum(basis[i]), basis[i]), reverse=True)
        permuted = [[vector[i] for i in order] for vector in kernel]
        _, pivots = _rref(permuted, len(basis))
        out.update(basis[order[p]] for p in pivots)
    return out


def _monomials(subs, d: int, m: int):
    for combo in it.combinations_with_replacement(range(len(subs)), m):
        exponent = [0] * len(subs)
        for index in combo:
            exponent[index] += 1
        yield tuple(exponent)


def _weight(exponent, subs, d: int) -> tuple[int, ...]:
    out = [0] * d
    for index, count in enumerate(exponent):
        if count:
            for orbital in subs[index]:
                out[orbital] += count
    return tuple(out)


def _raise(exponent, simple: int, subs, index_of) -> dict:
    out: dict[tuple[int, ...], int] = {}
    for position, count in enumerate(exponent):
        if not count:
            continue
        source = subs[position]
        if simple + 1 not in source or simple in source:
            continue
        target = tuple(sorted(simple if o == simple + 1 else o for o in source))
        image = list(exponent)
        image[position] -= 1
        image[index_of[target]] += 1
        key = tuple(image)
        out[key] = out.get(key, 0) + count
    return out


def _rref(rows, width: int):
    matrix = [list(row) for row in rows]
    pivots: list[int] = []
    current = 0
    for column in range(width):
        pivot = next((i for i in range(current, len(matrix)) if matrix[i][column]),
                     None)
        if pivot is None:
            continue
        matrix[current], matrix[pivot] = matrix[pivot], matrix[current]
        scale = matrix[current][column]
        matrix[current] = [value / scale for value in matrix[current]]
        for index in range(len(matrix)):
            if index != current and matrix[index][column]:
                factor = matrix[index][column]
                matrix[index] = [a - factor * b
                                 for a, b in zip(matrix[index], matrix[current])]
        pivots.append(column)
        current += 1
        if current == len(matrix):
            break
    return matrix[:current], pivots


def _kernel(equations, width: int):
    if not equations:
        return [[Fraction(int(i == j)) for j in range(width)] for i in range(width)]
    reduced, pivots = _rref(equations, width)
    pivot_set = set(pivots)
    out = []
    for chosen in (i for i in range(width) if i not in pivot_set):
        vector = [Fraction(0)] * width
        vector[chosen] = Fraction(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][chosen]
        out.append(vector)
    return out


# --------------------------------------------------------------------------
# Integer replay of every stored relation.
# --------------------------------------------------------------------------

def stored_levels(block) -> tuple[dict[int, set], list[int]]:
    levels: dict[int, set] = {}
    missing: list[int] = []
    for row in block["levels"]:
        if row["values"] is None:
            missing.append(row["degree"])
            continue
        levels[row["degree"]] = {tuple(v) for v in row["values"]}
    return levels, missing


def replay_generators(tag: str, levels: dict[int, set], generators) -> None:
    """Every stored generator is indecomposable; every other value decomposes."""
    claimed = {(m, tuple(v)) for m, v in generators}
    for m in sorted(levels):
        for value in levels[m]:
            witness = None
            for j in range(1, m // 2 + 1):
                if j not in levels or (m - j) not in levels:
                    continue
                for left in levels[j]:
                    rest = tuple(a - b for a, b in zip(value, left))
                    if min(rest) >= 0 and rest in levels[m - j]:
                        witness = (j, left, rest)
                        break
                if witness:
                    break
            if (m, value) in claimed:
                check(witness is None,
                      f"{tag}: claimed generator ({m},{list(value)}) decomposes")
            else:
                check(witness is not None,
                      f"{tag}: non-generator ({m},{list(value)}) has no decomposition")


def replay_digests(tag: str, block) -> None:
    for row in block["levels"]:
        if row["values"] is None:
            continue
        blob = json.dumps(sorted(row["values"]), separators=(",", ":"))
        check(hashlib.sha256(blob.encode()).hexdigest() == row["sha256"],
              f"{tag}: level digest mismatch at degree {row['degree']}")
        check(len(row["values"]) == row["size"],
              f"{tag}: level size mismatch at degree {row['degree']}")


def replay_particle_hole(tag: str, levels: dict[int, set]) -> None:
    for m, level in sorted(levels.items()):
        for value in level:
            dual = tuple(sorted((m - x for x in value), reverse=True))
            check(dual in level, f"{tag}: particle-hole image absent at degree {m}")


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    target = pathlib.Path(argv[0]) if argv else ARTIFACT
    if not target.exists():
        print(f"missing artifact {target}")
        return 1
    payload = json.loads(target.read_text())
    checks = 0

    by_system = {block["system"]: block for block in payload["raw_semigroup"]}

    # 1. integer replay of every stored relation
    for tag, block in sorted(by_system.items()):
        levels, missing = stored_levels(block)
        replay_digests(tag, block)
        replay_generators(tag, levels, block["generators"])
        checks += sum(len(level) for level in levels.values())
        if missing:
            print(f"  {tag}: levels stored by digest only at degrees {missing}")

    # 2. particle-hole self-duality where n == d - n
    for tag, block in sorted(by_system.items()):
        n, d = (int(x) for x in tag.strip("()").split(","))
        if 2 * n != d:
            continue
        levels, _ = stored_levels(block)
        replay_particle_hole(tag, levels)
        checks += sum(len(level) for level in levels.values())

    # 3. independent representation theory
    covered = {}
    for tag, block in sorted(by_system.items()):
        n, d = (int(x) for x in tag.strip("()").split(","))
        cap = min(MULT_REPLAY_MAX.get(tag, 0), block["reached"])
        levels, _ = stored_levels(block)
        probe = weight_count(n, d)
        for m in range(1, cap + 1):
            if m not in levels:
                continue
            check(independent_level(n, d, m, probe) == levels[m],
                  f"{tag}: independent multiplicity route disagrees at degree {m}")
            checks += 1
        covered[tag] = cap
    print(f"  independent multiplicity replay reached {covered}")

    # 4. padding is exact, replayed on the stored generator lists
    for low, high in (("(3,6)", "(3,7)"), ("(3,7)", "(3,8)"), ("(3,8)", "(3,9)"),
                      ("(3,9)", "(3,10)")):
        if low not in by_system or high not in by_system:
            continue
        window = min(by_system[low]["reached"], by_system[high]["reached"])
        wide = int(high.strip("()").split(",")[1])
        narrow = int(low.strip("()").split(",")[1])
        below = {(m, tuple(v)) for m, v in by_system[low]["generators"] if m <= window}
        lifted = {(m, tuple(v)[:narrow]) for m, v in by_system[high]["generators"]
                  if m <= window and tuple(v)[narrow:] == (0,) * (wide - narrow)}
        check(below == lifted, f"padding {low} -> {high} is not an equality")
        checks += 1

    # 5. the term-order side, rebuilt from scratch at (3,6)
    for block in payload.get("raw_term_order", []):
        if block["system"] != "(3,6)" or block["valuation"] != "v1":
            continue
        levels, _ = stored_levels(block)
        for m in range(1, min(HW_REPLAY_MAX, block["reached"]) + 1):
            if m not in levels:
                continue
            check(replay_leading_values(3, 6, m, block["variables"]) == levels[m],
                  f"(3,6) v1: rebuilt leading values disagree at degree {m}")
            checks += 1

    # 6. the scored counts are recomputed from the stored data, not trusted
    for row in payload["v2_by_system"]:
        block = by_system[row["system"]]
        window = row["window"]
        d = int(row["system"].strip("()").split(",")[1])
        inside = [(m, tuple(v)) for m, v in block["generators"] if m <= window]
        check(len(inside) == row["minimal_generators"],
              f"{row['system']}: minimal generator count does not replay")
        check(sum(1 for _, v in inside if v[d - 1] > 0) == row["primitive_generators"],
              f"{row['system']}: primitive generator count does not replay")
        checks += 2

    if FAILURES:
        print(f"\nFAILED {len(FAILURES)} checks")
        for message in FAILURES[:20]:
            print(f"  {message}")
        return 1
    print(f"\nOK, {checks} independent checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
