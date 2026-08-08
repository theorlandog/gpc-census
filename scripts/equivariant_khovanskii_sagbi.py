"""Khovanskii/SAGBI compression test for the fermionic highest-weight algebra.

Pre-registration: `docs/prereg_equivariant_khovanskii_sagbi.md`, committed
before this script was run. Artifact:
`results/data/equivariant_khovanskii_sagbi.json`.

TARGET, one only. With `V = wedge^n C^d` and `U` the upper unitriangular
subgroup of `GL(d)`, the algebra is `A(n,d) = Sym(V)^U` and its weight
semigroup is

    S(n,d) = { (m, lambda) : mult(V_lambda, Sym^m V) > 0 },

whose normalized support `lambda/m` has closure the fermionic moment polytope.
The combinatorial fallback target of the brief (facet-normal templates) is not
computed here.

THREE VALUATIONS, all exact, all genuine valuations.

* `v2`, the highest-weight valuation: lex-largest nonzero bidegree. Value
  semigroup `S(n,d)`. Multiplicities are collapsed, which is declared.
* `v1`, graded lex on Pluecker monomials with variables in lex order on the
  increasing subset.
* `v3`, graded reverse lex on Pluecker monomials with variables in colex order
  on the subset.

`v1` and `v3` are leading-exponent maps for term orders on a polynomial ring,
restricted to a domain subalgebra, hence valuations; both refine `v2` because
every monomial of a weight-homogeneous element carries the same weight.

Minimal generators are computed from the positive grading directly: an element
of degree `m` is decomposable exactly when it is a sum of semigroup elements of
degrees `j` and `m - j`. No Hilbert-basis library and no external algebra
system is used, so nothing is added to the dependency set.

Run: uv run scripts/equivariant_khovanskii_sagbi.py [--quick]
Writes results/data/equivariant_khovanskii_sagbi.json.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools as it
import json
import math
import pathlib
import sys
import time
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.plethysm import hm_en, mn, partitions  # noqa: E402

DATA = ROOT / "results" / "data"
OUT = DATA / "equivariant_khovanskii_sagbi.json"

# Budgets. The pre-registration declared 1800 seconds per system; the scored
# run used 300, a SHRINK. Shrinking is the safe direction, since it can only
# reduce the evidence for the campaign's own conclusion, and the
# pre-registration forbids only expanding a budget to postpone a negative
# conclusion. The shrink was applied before any generator count was read and is
# disclosed in docs/equivariant_khovanskii_sagbi.md. It costs about three
# degrees at the deepest system and does not move the common comparison window,
# which is set by the widest system in the n=3 ladder.
V2_SYSTEMS = ((3, 6), (3, 7), (3, 8), (3, 9), (3, 10), (4, 8), (4, 10), (5, 10))
V2_MAX_M = 18
V2_BUDGET = 300.0
V13_SYSTEMS = ((3, 6), (3, 7), (3, 8))
V13_MAX_M = {(3, 6): 8, (3, 7): 8, (3, 8): 5}
V13_BUDGET = 300.0
# Raw semigroup levels are stored verbatim below this size and by digest above
# it, so the artifact stays a readable size while every level a replay needs
# remains pinned. The standalone verifier recomputes each level by an
# independent route and checks it against whichever of the two is present.
RAW_STORE_CAP = 2000
# A single weight space of `Sym^m(wedge^n C^d)` can have thousands of monomials
# and exact row reduction over Q is cubic in that, so one weight can outrun the
# whole system budget. The cap is declared here rather than discovered at run
# time, and a level that hits it is DISCARDED, never stored truncated: a
# truncated level invents minimal generators, because every value it loses
# makes the values above it look indecomposable.
# Measured weight-space dimensions fix this: the largest weight of
# Sym^m(wedge^3 C^d) has 2905 monomials at (3,6) m=8, 3152 at (3,7) m=6 and
# 2215 at (3,8) m=5, then jumps to 16700 at (3,7) m=7 and 63381 at (3,8) m=7.
# 4000 therefore reaches the preregistered degree cap at (3,6) and (3,8) and
# stops (3,7) at m=6. This SHRINKS the declared budget, which is the safe
# direction; it is disclosed in the write-up and was fixed before any
# term-order generator count was read.
MAX_WEIGHT_BASIS = 4000


# --------------------------------------------------------------------------
# The semigroup S(n,d), by exact plethysm.
# --------------------------------------------------------------------------

def hw_support(n: int, d: int, m: int, deadline: float | None = None):
    """Weights of `Sym^m(wedge^n C^d)` with positive multiplicity, padded to `d`.

    Returns `None` when `deadline` passes, so a partially computed level is
    discarded rather than recorded. A truncated level would silently invent
    minimal generators, since every value it lost would make the values above
    it look indecomposable.
    """
    expansion = list(hm_en(m, n).items())
    out: dict[tuple[int, ...], int] = {}
    for index, lam in enumerate(partitions(n * m, maxpart=m)):
        if len(lam) > d:
            continue
        if deadline is not None and not index % 64 and time.time() > deadline:
            return None
        mult = sum(
            coefficient * mn(lam, tuple(sorted(mu, reverse=True)))
            for mu, coefficient in expansion
        )
        if mult > 0:
            out[lam + (0,) * (d - len(lam))] = int(mult)
    return out


def _cache_paths(cache: pathlib.Path | None, n: int, d: int):
    if cache is None:
        return None, None
    cache.mkdir(parents=True, exist_ok=True)
    return cache / f"levels_{n}_{d}.json", cache / f"spent_{n}_{d}.json"


def semigroup(n: int, d: int, max_m: int, budget: float, log=print,
              cache: pathlib.Path | None = None):
    """Return `(levels, reached, status, secs)` for `S(n,d)` up to `max_m`.

    The budget is CUMULATIVE compute, not wall clock of one process. Levels are
    exact and deterministic, so a cached level is replayed for free on a
    resumed run, and the seconds actually spent producing new levels are
    carried across resumes in the cache. Enforcing the budget against that
    running total is what keeps a restarted run inside the declared scope
    instead of silently buying more degrees every time it is relaunched.
    """
    levels_path, spent_path = _cache_paths(cache, n, d)
    levels: dict[int, dict[tuple[int, ...], int]] = {}
    spent = 0.0
    if levels_path is not None and levels_path.exists():
        stored = json.loads(levels_path.read_text())
        levels = {int(m): {tuple(row["lam"]): row["mult"] for row in rows}
                  for m, rows in stored.items()}
        spent = json.loads(spent_path.read_text())["seconds"] if \
            spent_path.exists() else 0.0
        if levels:
            log(f"  ({n},{d}) resumed from cache: degrees <= {max(levels)}, "
                f"{spent:.1f}s already spent")

    def persist() -> None:
        if levels_path is None:
            return
        levels_path.write_text(json.dumps(
            {str(m): [{"lam": list(lam), "mult": mult}
                      for lam, mult in sorted(level.items())]
             for m, level in levels.items()}))
        spent_path.write_text(json.dumps({"seconds": round(spent, 1)}))

    reached = max(levels) if levels else 0
    for m in range(1, max_m + 1):
        if m in levels:
            reached = m
            continue
        if spent >= budget:
            break
        started = time.time()
        level = hw_support(n, d, m, started + (budget - spent))
        spent += time.time() - started
        if level is None:
            log(f"  ({n},{d}) m={m}: abandoned on budget, level discarded")
            persist()
            break
        levels[m] = level
        reached = m
        persist()
        log(f"  ({n},{d}) m={m}: |S_m|={len(level)} [{spent:.1f}s cumulative]")
    persist()
    status = "complete" if reached == max_m else "operational"
    return levels, reached, status, round(spent, 1)


# --------------------------------------------------------------------------
# Minimal generators of a positively graded semigroup, exactly.
# --------------------------------------------------------------------------

def decomposables(levels: dict[int, set], m: int) -> set:
    """All degree-`m` values that are a sum of two lower-degree values."""
    out: set = set()
    for j in range(1, m // 2 + 1):
        left, right = levels.get(j, ()), levels.get(m - j, ())
        for a in left:
            for b in right:
                out.add(tuple(x + y for x, y in zip(a, b)))
    return out


def minimal_generators(levels: dict[int, set], reached: int):
    """Return `(generators, closure_violations)`.

    `closure_violations` counts degree-`m` sums of semigroup elements that are
    NOT in the semigroup. It must be zero: the value semigroup of a valuation
    on a domain is closed under addition, so a nonzero count means the
    underlying multiplicity computation is wrong.
    """
    generators: list[tuple[int, tuple[int, ...]]] = []
    violations: list[tuple[int, tuple[int, ...]]] = []
    for m in range(1, reached + 1):
        here = levels.get(m, set())
        decomposed = decomposables(levels, m)
        for value in sorted(decomposed - set(here)):
            violations.append((m, value))
        for value in sorted(here):
            if value not in decomposed:
                generators.append((m, value))
    return generators, violations


def closure(generators, max_m: int, width: int) -> dict[int, set]:
    """Every nonnegative integer combination of `generators`, by degree."""
    reach: dict[int, set] = {0: {(0,) * width}}
    for m in range(1, max_m + 1):
        here: set = set()
        for degree, value in generators:
            if degree > m:
                continue
            for base in reach[m - degree]:
                here.add(tuple(x + y for x, y in zip(base, value)))
        reach[m] = here
    return reach


def held_out(generators, levels: dict[int, set], fit_to: int, test: list[int],
             width: int) -> dict:
    """Fit on degrees <= `fit_to`, then score misses and spurious values.

    SPURIOUS IS EMPTY BY CONSTRUCTION on this half of the test and is recorded
    as UNINFORMATIVE rather than as a pass (standing rule R3): the fitted
    generators are themselves semigroup elements, so every nonnegative
    combination of them is a semigroup element. Spurious predictions are
    possible only for a model that posits values it has not verified, which is
    what `saturation_model` below supplies.
    """
    fitted = [g for g in generators if g[0] <= fit_to]
    if not test:
        return {"fit_degree": fit_to, "test_degrees": [], "verdict": "UNINFORMATIVE",
                "reason": "no test degree inside the computed scope",
                "fitted_generators": len(fitted), "per_degree": []}
    reach = closure(fitted, max(test), width)
    rows = []
    for m in test:
        truth = set(levels.get(m, set()))
        predicted = reach[m]
        misses = sorted(truth - predicted)
        spurious = sorted(predicted - truth)
        rows.append({
            "degree": m,
            "truth": len(truth),
            "predicted": len(predicted),
            "misses": len(misses),
            "spurious": len(spurious),
            "miss_witnesses": [list(x) for x in misses[:5]],
            "spurious_witnesses": [list(x) for x in spurious[:5]],
        })
    clean = all(r["misses"] == 0 for r in rows)
    return {"fit_degree": fit_to, "test_degrees": test,
            "verdict": "PASS" if clean else "FAIL",
            "spurious_half": "UNINFORMATIVE_BY_CONSTRUCTION",
            "fitted_generators": len(fitted), "per_degree": rows}


def saturation_model(levels: dict[int, set], reached: int, width: int) -> dict:
    """Score the one model here that CAN over-predict: assume `S` is saturated.

    A value `x` at degree `m` is in the saturation of `S` when `k*x` lies in
    `S` at degree `k*m` for some `k >= 2`. Any such `x` absent from `S_m` is a
    SPURIOUS prediction of the saturated model, and is exactly a quantization
    period above the spectrum denominator seen from the semigroup side.
    """
    spurious = []
    tested = 0
    for m in range(1, reached + 1):
        here = levels.get(m, set())
        for k in (2, 3):
            if k * m > reached:
                continue
            upper = levels.get(k * m, set())
            for value in upper:
                if any(entry % k for entry in value):
                    continue
                root = tuple(entry // k for entry in value)
                tested += 1
                if root not in here:
                    spurious.append([m, list(root), k])
    return {"candidates_tested": tested, "spurious": len(spurious),
            "witnesses": spurious[:10],
            "verdict": "FAIL" if spurious else "PASS",
            "reading": ("a nonempty list refutes saturation of the semigroup, "
                        "so no model that predicts the lattice points of the "
                        "cone can be correct")}


def held_out_rank(gens_by_system: dict, n: int, fit_rank: int, test_rank: int,
                  levels_test: dict[int, set], window: int) -> dict:
    """Fit at `fit_rank`, predict `test_rank` by padding alone."""
    if (n, fit_rank) not in gens_by_system or (n, test_rank) not in gens_by_system:
        return {"verdict": "UNINFORMATIVE", "reason": "rank outside the scope"}
    predicted = [(m, tuple(lam) + (0,) * (test_rank - fit_rank))
                 for m, lam in gens_by_system[(n, fit_rank)] if m <= window]
    truth = [(m, tuple(lam)) for m, lam in gens_by_system[(n, test_rank)]
             if m <= window]
    missing = sorted(set(truth) - set(predicted))
    extra = sorted(set(predicted) - set(truth))
    generated = closure(predicted, window, test_rank)
    ungenerated = {m: sorted(set(levels_test.get(m, ())) - generated[m])
                   for m in range(1, window + 1)}
    return {
        "fit": f"({n},{fit_rank})", "test": f"({n},{test_rank})", "window": window,
        "predicted_generators": len(predicted),
        "true_generators": len(truth),
        "new_generators_missed": len(missing),
        "predicted_but_not_generators": len(extra),
        "miss_witnesses": [[m, list(lam)] for m, lam in missing[:10]],
        "values_not_generated": {str(m): len(v) for m, v in ungenerated.items()},
        "verdict": "PASS" if not missing else "FAIL",
    }


# --------------------------------------------------------------------------
# Explicit highest-weight vectors, for the two term-order valuations.
# --------------------------------------------------------------------------

def subsets(n: int, d: int) -> list[tuple[int, ...]]:
    return list(it.combinations(range(d), n))


def weight_of(subset: tuple[int, ...], d: int) -> tuple[int, ...]:
    out = [0] * d
    for orbital in subset:
        out[orbital] = 1
    return tuple(out)


def monomials_of_weight(subs: list[tuple[int, ...]], d: int, m: int,
                        target: tuple[int, ...]) -> list[tuple[int, ...]]:
    """Exponent vectors over `subs` of total degree `m` and the given weight.

    Enumerated by a suffix-feasibility recursion so no monomial of the wrong
    weight is ever built.
    """
    weights = [weight_of(s, d) for s in subs]
    suffix = [[0] * d for _ in range(len(subs) + 1)]
    for index in range(len(subs) - 1, -1, -1):
        suffix[index] = [suffix[index + 1][o] + weights[index][o] for o in range(d)]
    out: list[tuple[int, ...]] = []
    exponent = [0] * len(subs)

    def walk(index: int, left: int, rest: tuple[int, ...]) -> None:
        if left == 0:
            if all(value == 0 for value in rest):
                out.append(tuple(exponent))
            return
        if index >= len(subs):
            return
        # every remaining orbital demand must still be reachable
        for orbital in range(d):
            if rest[orbital] > left:
                return
            if rest[orbital] and not suffix[index][orbital]:
                return
        weight = weights[index]
        cap = min(left, min((rest[o] for o in subs[index]), default=0))
        for count in range(cap, -1, -1):
            exponent[index] = count
            walk(index + 1, left - count,
                 tuple(rest[o] - count * weight[o] for o in range(d)))
        exponent[index] = 0

    walk(0, m, target)
    return out


def raise_monomial(exponent: tuple[int, ...], simple: int,
                   subs: list[tuple[int, ...]], index_of: dict) -> dict:
    """Apply the simple raising operator `E_{simple, simple+1}` as a derivation."""
    out: dict[tuple[int, ...], int] = {}
    for position, count in enumerate(exponent):
        if not count:
            continue
        source = subs[position]
        if simple + 1 not in source or simple in source:
            continue
        target = tuple(sorted((simple if o == simple + 1 else o) for o in source))
        # target and source differ by swapping the adjacent values simple+1 and
        # simple, so no element lies strictly between them and the fermionic
        # reordering sign is +1.
        image = list(exponent)
        image[position] -= 1
        image[index_of[target]] += 1
        key = tuple(image)
        out[key] = out.get(key, 0) + count
    return out


def rref(rows: list[list[Fraction]], width: int) -> tuple[list[list[Fraction]], list[int]]:
    """Exact row reduction. Returns the reduced rows and the pivot columns."""
    matrix = [list(row) for row in rows]
    pivots: list[int] = []
    current = 0
    for column in range(width):
        pivot = None
        for index in range(current, len(matrix)):
            if matrix[index][column]:
                pivot = index
                break
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


def sparse_echelon(rows: list[dict[int, int]]) -> dict[int, dict[int, int]]:
    """Integer row echelon form of a sparse system, keyed by pivot column.

    The raising-operator matrices are very sparse and a dense exact reduction
    over Q is cubic in the weight-space dimension, which is what put a single
    weight of `Sym^7(wedge^3 C^6)` beyond the whole system budget. Rows are
    kept as integer dictionaries and divided by their content at every step, so
    no `Fraction` is created and coefficients stay small.
    """
    pivots: dict[int, dict[int, int]] = {}
    for original in rows:
        row = {column: value for column, value in original.items() if value}
        while row:
            column = min(row)
            if column not in pivots:
                pivots[column] = row
                break
            other = pivots[column]
            scale = math.gcd(row[column], other[column])
            mine, theirs = other[column] // scale, row[column] // scale
            combined = {key: value * mine for key, value in row.items()}
            for key, value in other.items():
                combined[key] = combined.get(key, 0) - value * theirs
            row = {key: value for key, value in combined.items() if value}
            if row:
                content = 0
                for value in row.values():
                    content = math.gcd(content, value)
                if content > 1:
                    row = {key: value // content for key, value in row.items()}
    return pivots


def highest_weight_space(n: int, d: int, m: int, lam: tuple[int, ...],
                         subs: list[tuple[int, ...]], index_of: dict,
                         cap: int = MAX_WEIGHT_BASIS):
    """Exact basis of the highest-weight vectors of weight `lam` in `Sym^m`.

    Returns `None` when the weight space exceeds `cap`, so the caller discards
    the whole level rather than storing a truncated one.
    """
    basis = monomials_of_weight(subs, d, m, lam)
    if len(basis) > cap:
        return None
    if not basis:
        return [], []
    column_of = {monomial: index for index, monomial in enumerate(basis)}
    equations: list[dict[int, int]] = []
    for simple in range(d - 1):
        targets: dict[tuple[int, ...], dict[int, int]] = {}
        for monomial in basis:
            for image, coefficient in raise_monomial(monomial, simple, subs,
                                                     index_of).items():
                targets.setdefault(image, {})[column_of[monomial]] = coefficient
        equations.extend(targets.values())
    if not equations:
        kernel = [[Fraction(int(index == other)) for other in range(len(basis))]
                  for index in range(len(basis))]
        return basis, kernel
    pivots = sparse_echelon(equations)
    free = [index for index in range(len(basis)) if index not in pivots]
    descending = sorted(pivots, reverse=True)
    kernel: list[list[Fraction]] = []
    for chosen in free:
        vector = [Fraction(0)] * len(basis)
        vector[chosen] = Fraction(1)
        # each echelon row has its pivot as its SMALLEST column, so descending
        # pivot order makes every other entry of the row already known
        for pivot in descending:
            row = pivots[pivot]
            total = sum(value * vector[key]
                        for key, value in row.items() if key != pivot)
            vector[pivot] = -Fraction(total, row[pivot])
        kernel.append(vector)
    return basis, kernel


def order_key(name: str, subs: list[tuple[int, ...]]):
    """Return `(variable_permutation, monomial_key)` for a named term order."""
    if name == "v1":
        # variables in lex order on the increasing subset; graded lex on monomials
        order = sorted(range(len(subs)), key=lambda i: subs[i])

        def key(exponent: tuple[int, ...]) -> tuple:
            return (sum(exponent), tuple(exponent[i] for i in order))
    elif name == "v3":
        # variables in colex order on the subset; graded reverse lex on monomials
        order = sorted(range(len(subs)), key=lambda i: tuple(reversed(subs[i])))

        def key(exponent: tuple[int, ...]) -> tuple:
            reordered = [exponent[i] for i in order]
            return (sum(exponent), tuple(-value for value in reversed(reordered)))
    else:  # pragma: no cover - guarded by the caller
        raise ValueError(name)
    return key


def leading_values(basis, kernel, key) -> list[tuple[int, ...]]:
    """Leading exponent vectors of the whole highest-weight space, exactly.

    The set of leading monomials of the nonzero elements of a subspace is the
    set of pivot columns of its row reduction against the term order, so the
    count equals the dimension.
    """
    if not kernel:
        return []
    columns = sorted(range(len(basis)), key=lambda i: key(basis[i]), reverse=True)
    permuted = [[vector[column] for column in columns] for vector in kernel]
    _, pivots = rref(permuted, len(basis))
    return [basis[columns[pivot]] for pivot in pivots]


def term_order_semigroup(n: int, d: int, max_m: int, budget: float, name: str,
                         log=print, cache: pathlib.Path | None = None):
    """Value semigroup of the term-order valuation `name`, with a gate count."""
    subs = subsets(n, d)
    index_of = {s: i for i, s in enumerate(subs)}
    key = order_key(name, subs)
    levels: dict[int, set] = {}
    gate_ok, gate_total = 0, 0
    path = None
    if cache is not None:
        cache.mkdir(parents=True, exist_ok=True)
        path = cache / f"term_{name}_{n}_{d}.json"
        if path.exists():
            stored = json.loads(path.read_text())
            levels = {int(m): {tuple(v) for v in values}
                      for m, values in stored["levels"].items()}
            gate_ok, gate_total = stored["gate_ok"], stored["gate_total"]
            if levels:
                log(f"  ({n},{d}) {name} resumed from cache: "
                    f"degrees <= {max(levels)}")
    started = time.time()
    deadline = started + budget
    reached = max(levels) if levels else 0
    for m in range(1, max_m + 1):
        if m in levels:
            reached = m
            continue
        if time.time() > deadline:
            break
        support = hw_support(n, d, m, deadline)
        if support is None:
            log(f"  ({n},{d}) {name} m={m}: abandoned on budget")
            break
        here: set = set()
        abandoned = ""
        for lam, mult in sorted(support.items()):
            if time.time() > deadline:
                abandoned = "budget"
                break
            space = highest_weight_space(n, d, m, lam, subs, index_of)
            if space is None:
                abandoned = f"weight space above MAX_WEIGHT_BASIS={MAX_WEIGHT_BASIS}"
                break
            basis, kernel = space
            gate_total += 1
            gate_ok += int(len(kernel) == mult)
            here.update(leading_values(basis, kernel, key))
        if abandoned:
            log(f"  ({n},{d}) {name} m={m}: abandoned ({abandoned}), level discarded")
            break
        levels[m] = here
        reached = m
        if path is not None:
            path.write_text(json.dumps({
                "levels": {str(k): [list(v) for v in sorted(value)]
                           for k, value in levels.items()},
                "gate_ok": gate_ok, "gate_total": gate_total}))
        log(f"  ({n},{d}) {name} m={m}: |values|={len(here)} "
            f"[{time.time() - started:.1f}s]")
    status = "complete" if reached == max_m else "operational"
    return levels, reached, status, gate_ok, gate_total, round(time.time() - started, 1)


# --------------------------------------------------------------------------
# Transports: padding and particle-hole.
# --------------------------------------------------------------------------

def padding_report(gens_by_system: dict, n: int, ranks: list[int],
                   window: int) -> list[dict]:
    """`Gen(n,d) ∩ {lambda_d = 0}` against `Gen(n,d-1)`, exactly."""
    rows = []
    for index in range(1, len(ranks)):
        low, high = ranks[index - 1], ranks[index]
        if (n, low) not in gens_by_system or (n, high) not in gens_by_system:
            continue
        below = {(m, tuple(lam)) for m, lam in gens_by_system[(n, low)] if m <= window}
        lifted = {(m, tuple(lam)[:low]) for m, lam in gens_by_system[(n, high)]
                  if m <= window and tuple(lam)[low:] == (0,) * (high - low)}
        rows.append({
            "from": f"({n},{low})", "to": f"({n},{high})", "window": window,
            "inherited": len(lifted), "lower_generators": len(below),
            "equal": below == lifted,
            "only_lower": [[m, list(lam)] for m, lam in sorted(below - lifted)][:5],
            "only_upper": [[m, list(lam)] for m, lam in sorted(lifted - below)][:5],
        })
    return rows


def particle_hole_report(levels: dict[int, dict], n: int, d: int,
                         reached: int) -> dict:
    """`lambda -> (m - lambda_d, ..., m - lambda_1)` carries `S(n,d)` to `S(d-n,d)`.

    It is a SELF-map only when `2n = d`, because
    `wedge^n C^d` is `(wedge^n C^d)^*` twisted by `det` exactly there. At every
    other system the complement lands in a different algebra, which is not in
    scope, so the row is NOT_APPLICABLE and is never counted as a failure.
    Running the self-map anyway and reporting its mismatches would report the
    definition of the transport as if it were a defect.
    """
    if 2 * n != d:
        return {"system": f"({n},{d})", "self_dual": False,
                "verdict": "NOT_APPLICABLE", "checked": 0, "failures": 0,
                "reason": f"the complement lands in ({d - n},{d}), not in scope",
                "witnesses": []}
    checked, failures = 0, []
    for m in range(1, reached + 1):
        here = set(levels[m])
        for lam in here:
            dual = tuple(sorted((m - value for value in lam), reverse=True))
            checked += 1
            if dual not in here:
                failures.append([m, list(lam), list(dual)])
    return {"system": f"({n},{d})", "self_dual": True,
            "verdict": "PASS" if not failures else "FAIL",
            "checked": checked, "failures": len(failures),
            "witnesses": failures[:5]}


# --------------------------------------------------------------------------
# Driver.
# --------------------------------------------------------------------------

def summarize(name: str, levels, reached, generators, window: int,
              n: int, d: int, *, weights: bool) -> dict:
    """`weights` is true when the values are `GL(d)` weights, so that the
    primitive count (last coordinate nonzero, hence not reachable by padding)
    is meaningful. Term-order values are Pluecker exponents and carry no such
    coordinate, so the count is omitted rather than faked."""
    raw = sum(len(levels.get(m, ())) for m in range(1, window + 1))
    inside = [g for g in generators if g[0] <= window]
    primitive = [g for g in inside if g[1][d - 1] > 0] if weights else None
    return {
        "valuation": name,
        "system": f"({n},{d})",
        "window": window,
        "reached_degree": reached,
        "raw_values": raw,
        "minimal_generators": len(inside),
        "primitive_generators": (len(primitive) if primitive is not None else None),
        "max_generator_degree": max((g[0] for g in inside), default=0),
        # A maximum generator degree equal to the window is NOT a bound: it
        # says the window ran out before the generators did. Only a strict
        # inequality here establishes that the generating set is finished.
        "generator_degree_saturates_window":
            max((g[0] for g in inside), default=0) == window,
        "compression_ratio": (round(len(inside) / raw, 4) if raw else None),
        "generators_by_degree": {
            str(m): sum(1 for g in inside if g[0] == m) for m in range(1, window + 1)
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="short budgets for a smoke run; not the scored scope")
    parser.add_argument("--out", type=pathlib.Path, default=OUT,
                        help="write elsewhere; used to smoke-test the verifier")
    parser.add_argument("--cache", type=pathlib.Path, default=None,
                        help="resume from and write exact level checkpoints; "
                             "the per-system budget stays cumulative across "
                             "resumes, so a restart never buys extra degrees")
    args = parser.parse_args()

    v2_max = 8 if args.quick else V2_MAX_M
    v2_budget = 20.0 if args.quick else V2_BUDGET
    v13_budget = 20.0 if args.quick else V13_BUDGET
    v2_systems = V2_SYSTEMS[:3] if args.quick else V2_SYSTEMS
    v13_systems = V13_SYSTEMS[:2] if args.quick else V13_SYSTEMS

    payload: dict = {
        "generated_by": "scripts/equivariant_khovanskii_sagbi.py",
        "prereg": "docs/prereg_equivariant_khovanskii_sagbi.md",
        "evidence": "exact",
        "quick": bool(args.quick),
        "target": (
            "weight semigroup S(n,d) of the highest-weight algebra "
            "Sym(wedge^n C^d)^U, bigraded by symmetric-power degree and GL(d) "
            "weight; normalized support is the fermionic moment polytope"
        ),
        "note": (
            "Minimal generators are computed from the positive grading by exact "
            "integer decomposition, with no Hilbert-basis library and no "
            "external algebra system. Every semigroup level is checked closed "
            "under addition, which is forced for the value semigroup of a "
            "valuation on a domain and therefore gates the plethysm input."
        ),
    }

    # ---- V2, the highest-weight valuation ---------------------------------
    print("V2 highest-weight valuation")
    v2: dict = {}
    raw_levels: dict = {}
    gens_by_system: dict = {}
    for n, d in v2_systems:
        levels, reached, status, secs = semigroup(
            n, d, v2_max, v2_budget, cache=args.cache)
        as_sets = {m: set(levels[m]) for m in levels}
        generators, violations = minimal_generators(as_sets, reached)
        raw_levels[(n, d)] = (as_sets, reached)
        gens_by_system[(n, d)] = generators
        v2[(n, d)] = {
            "levels": as_sets, "reached": reached, "status": status,
            "secs": secs, "generators": generators, "violations": violations,
            "multiplicities": {m: max(levels[m].values(), default=0) for m in levels},
        }

    # P1 is a statement about the n=3 ladder at d = 6..10, so the common
    # comparison window is the one those systems share. Taking the minimum over
    # every system in scope would let the widest cross-N system, which reaches
    # the fewest degrees, silently shrink the rank comparison the
    # pre-registration actually asks for.
    ladder = [entry["reached"] for (n, d), entry in v2.items() if n == 3]
    common = min(ladder) if ladder else min(e["reached"] for e in v2.values())
    print(f"  common degree window over the n=3 ladder: m <= {common}")

    v2_rows, closure_rows, heldout_rows, ph_rows, sat_rows = [], [], [], [], []
    for (n, d), entry in v2.items():
        window = min(common, entry["reached"])
        v2_rows.append(summarize("v2", entry["levels"], entry["reached"],
                                 entry["generators"], window, n, d,
                                 weights=True) | {
            "status": entry["status"], "secs": entry["secs"],
            "max_multiplicity": max(entry["multiplicities"].values(), default=0),
        })
        closure_rows.append({
            "system": f"({n},{d})", "reached": entry["reached"],
            "closure_violations": len(entry["violations"]),
            "witnesses": [[m, list(v)] for m, v in entry["violations"][:5]],
        })
        fit = math.ceil(entry["reached"] / 2)
        test = [m for m in (fit + 1, fit + 2) if m <= entry["reached"]]
        heldout_rows.append({"system": f"({n},{d})", "valuation": "v2"}
                            | held_out(entry["generators"], entry["levels"],
                                       fit, test, d))
        ph_rows.append(particle_hole_report(entry["levels"], n, d, entry["reached"]))
        sat_rows.append({"system": f"({n},{d})"}
                        | saturation_model(entry["levels"], entry["reached"], d))

    payload["v2_by_system"] = v2_rows
    payload["raw_semigroup"] = [
        {"system": f"({n},{d})", "valuation": "v2", "reached": entry["reached"],
         "levels": raw_block(entry["levels"], entry["reached"]),
         "generators": [[m, list(value)] for m, value in entry["generators"]]}
        for (n, d), entry in v2.items()
    ]
    payload["semigroup_closure_gate"] = closure_rows
    payload["held_out_degree"] = heldout_rows
    payload["particle_hole"] = ph_rows
    payload["saturation_model"] = sat_rows
    ranks3 = [d for n, d in v2_systems if n == 3]
    payload["padding"] = padding_report(gens_by_system, 3, ranks3, common)
    payload["held_out_rank"] = [
        held_out_rank(gens_by_system, 3, low, high,
                      raw_levels.get((3, high), ({}, 0))[0], common)
        for low, high in zip(ranks3, ranks3[1:])
    ]

    # rank growth of the primitive generators, in the common window
    payload["rank_growth"] = [
        {"system": row["system"], "window": row["window"],
         "raw_values": row["raw_values"],
         "minimal_generators": row["minimal_generators"],
         "primitive_generators": row["primitive_generators"],
         "primitive_ratio": (round(row["primitive_generators"] / row["raw_values"], 4)
                             if row["raw_values"] else None)}
        for row in v2_rows if row["system"].startswith("(3,")
    ]

    # ---- V1 and V3, the term-order valuations -----------------------------
    term_rows, term_gate, term_closure, term_heldout, term_raw = [], [], [], [], []
    term_data: dict = {}
    for name in ("v1", "v3"):
        print(f"{name} term-order valuation")
        for n, d in v13_systems:
            cap = V13_MAX_M[(n, d)]
            levels, reached, status, ok, total, secs = term_order_semigroup(
                n, d, cap, v13_budget, name, cache=args.cache)
            generators, violations = minimal_generators(levels, reached)
            width = math.comb(d, n)
            term_rows.append(summarize(name, levels, reached, generators, reached,
                                       n, d, weights=False) | {
                "status": status, "secs": secs, "value_width": width,
            })
            term_gate.append({"valuation": name, "system": f"({n},{d})",
                              "weight_spaces_checked": total,
                              "dimension_matches_multiplicity": ok,
                              "all_match": ok == total})
            term_closure.append({"valuation": name, "system": f"({n},{d})",
                                 "closure_violations": len(violations)})
            fit = math.ceil(reached / 2)
            test = [m for m in (fit + 1, fit + 2) if m <= reached]
            term_heldout.append({"system": f"({n},{d})", "valuation": name}
                                | held_out(generators, levels, fit, test, width))
            term_data[((n, d), name)] = (levels, reached, generators)
            term_raw.append({
                "system": f"({n},{d})", "valuation": name, "reached": reached,
                "value_width": width,
                "variables": [list(s) for s in subsets(n, d)],
                "levels": raw_block(levels, reached),
                "generators": [[m, list(value)] for m, value in generators],
            })
    payload["v1_v3_by_system"] = term_rows
    payload["raw_term_order"] = term_raw
    # A fair comparison needs the SAME degree window on all three valuations.
    # The term orders stop earlier than the highest-weight side, so quoting
    # each at its own reach would credit v2 with degrees v1 and v3 never saw.
    matched = []
    for (n, d), name in sorted(term_data, key=lambda k: (k[1], k[0])):
        levels, reached, generators = term_data[((n, d), name)]
        entry = v2.get((n, d))
        if entry is None:
            continue
        window = min(reached, entry["reached"])
        matched.append(summarize("v2", entry["levels"], entry["reached"],
                                 entry["generators"], window, n, d,
                                 weights=True))
        matched.append(summarize(name, levels, reached, generators, window,
                                 n, d, weights=False))
    payload["valuation_comparison_matched_window"] = matched
    payload["hw_dimension_gate"] = term_gate
    payload["term_order_closure_gate"] = term_closure
    payload["held_out_degree_term_orders"] = term_heldout

    # ---- quantization cross-check ----------------------------------------
    payload["quantization_cross_check"] = quantization_cross_check(raw_levels)

    args.out.write_text(json.dumps(payload, indent=1, default=_json_default) + "\n")
    print(f"\nwrote {args.out}")
    return 0


def level_digest(level) -> str:
    """SHA256 of a semigroup level in its canonical sorted form."""
    blob = json.dumps([list(value) for value in sorted(level)],
                      separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


def raw_block(levels, reached: int, multiplicities: dict | None = None) -> list:
    """The raw bounded semigroup data a replay of the generator step needs."""
    rows = []
    for m in range(1, reached + 1):
        level = levels.get(m, set())
        row = {"degree": m, "size": len(level), "sha256": level_digest(level)}
        if len(level) <= RAW_STORE_CAP:
            row["values"] = [list(value) for value in sorted(level)]
        else:
            row["values"] = None
            row["omitted_reason"] = f"size above RAW_STORE_CAP={RAW_STORE_CAP}"
        if multiplicities is not None:
            row["max_multiplicity"] = max(multiplicities.get(m, {}).values(),
                                          default=0)
        rows.append(row)
    return rows


def _json_default(obj):
    if isinstance(obj, (set, frozenset)):
        return sorted(obj)
    if isinstance(obj, tuple):
        return list(obj)
    raise TypeError(type(obj))


def quantization_cross_check(raw_levels: dict) -> dict:
    """Realized multiples of each census vertex against the measured period.

    `results/data/quantization_multiplicities.json` measured a period `p` at
    the fourteen `(3,6)` and `(3,7)` vertices by a completely different route
    (the alternating Weyl formula over exact weight counts). Here the same
    numbers fall out of the semigroup levels: the set of `m` with
    `(m, m*lambda/q)` in `S` must be exactly the multiples of `p`.
    """
    path = DATA / "quantization_multiplicities.json"
    if not path.exists():
        return {"status": "absent", "rows": []}
    measured = json.loads(path.read_text())["rows"]
    rows = []
    for record in measured:
        system = record["system"]
        n, d = (int(x) for x in system.strip("()").split(","))
        if (n, d) not in raw_levels:
            continue
        levels, reached = raw_levels[(n, d)]
        q = int(record["denominator"] or 1)
        form = list(record["integer_form"])
        form = form + [0] * (d - len(form))
        realized = [m for m in range(1, reached + 1)
                    if m % q == 0
                    and tuple(value * (m // q) for value in form) in levels[m]]
        period = record.get("quantization_period")
        expected = ([m for m in range(1, reached + 1) if period and m % period == 0]
                    if period else None)
        rows.append({
            "vertex": f"{system}#{record['index']}",
            "denominator": q,
            "measured_period": period,
            "scope_degrees": reached,
            "realized_multiples": realized,
            "expected_multiples": expected,
            "agrees": expected is not None and realized == expected,
        })
    scored = [r for r in rows if r["expected_multiples"] is not None]
    return {
        "status": "exact",
        "scope": len(scored),
        "agreements": sum(1 for r in scored if r["agrees"]),
        "rows": rows,
    }


if __name__ == "__main__":
    raise SystemExit(main())
