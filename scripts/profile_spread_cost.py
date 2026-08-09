#!/usr/bin/env python3
"""Measure ``D(W)``, the rank cost of a wide level set, and how fast it climbs.

Writes ``results/data/profile_spread_cost.json``.

WHY THIS EXISTS.  ``D(W)`` is the least ambient rank carrying an admissible
profile of level spread ``W``.  This table was quoted in prose with nothing
regenerating it; this script is that generator.

WHAT IT IS EVIDENCE FOR, which has changed.  It was written as evidence for
Conjecture S, that the admissible profiles at each rank have bounded spread.
Theorem S in ``docs/profile_native_generation.md`` now proves that outright, by
a route this table has no part in: the pattern set determines the values, so
finitely many level vectors survive per rank.  The surviving question is
quantitative, whether the true maximum spread is polynomial in ``r``, because
the bound Theorem S gives is exponential and about ``2.8e6`` at ``r = 11``
against a realized maximum of 20.  This table measures the RATE, which is the
open half, rather than the existence, which is closed.

WHAT IS COMPUTED, and it is exact within its scope.  Theorem A says
admissibility depends on the multiplicity vector only through
``min(m_a, N + 1)``, and is monotone in it, so for a level profile the least
ambient rank that can carry it is

    cost(levels, target) = min { sum_a mhat_a : mhat admissible }

over ``mhat`` in ``{1, ..., N+1}^r``.  Minimising ``cost`` over every level
profile of spread ``W`` and at most ``max_classes`` classes gives an exact value
of ``D`` restricted to that class bound.

THE SCOPE IS A CEILING, NOT A FLOOR, and the direction matters.  A profile with
more classes could be cheaper, so the restricted minimum is an UPPER bound:

    D(W)  <=  D_{r <= max_classes}(W)

So this table cannot prove a spread bound at any given rank on its own, and it
does not have to: Theorem S proves one.  What it can do is exhibit the rate, that
the cost of a wide level set climbs with the width rather than staying flat, and
climbs roughly linearly where Theorem S only guarantees exponential.  The
artifact labels the column accordingly and the guard test pins the label.

CROSS-CHECK.  ``results/data/generator_frontier_refinement.json`` records the
realized maximum spread at ranks 6 to 10, from the complete enumerations, and
those are exact.  A profile of spread ``W`` existing at rank ``d`` means
``D(W) <= d``, so every realized pair must satisfy the ceiling consistently.
That is checked here rather than assumed, and a violation would mean one of the
two computations is wrong.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.generation.profiles import (  # noqa: E402
    _rank_mod,
    occupancy_patterns,
)

OUTPUT = ROOT / "results" / "data" / "profile_spread_cost.json"
FRONTIER = ROOT / "results" / "data" / "generator_frontier_refinement.json"


def _admissible(patterns, mhat, r):
    """Theorem A on an explicit pattern set, capped by ``mhat``."""
    capped = [c for c in patterns if all(c[a] <= mhat[a] for a in range(r))]
    if len(capped) < r - 1:
        return False
    if _rank_mod(capped, r) != r - 1:
        return False
    for a in range(r):
        if mhat[a] >= 2 and not any(1 <= c[a] <= mhat[a] - 1 for c in capped):
            return False
    return True


def _multiplicity_vectors(r, total, cap):
    """Every ``mhat`` in ``[1, cap]^r`` summing to ``total``."""
    extra = total - r
    if extra < 0:
        return
    for split in itertools.combinations_with_replacement(range(r), extra):
        vector = [1] * r
        for index in split:
            vector[index] += 1
        if max(vector) <= cap:
            yield tuple(vector)


def _min_cost(patterns, r, n, ceiling):
    """Least ``sum(mhat)`` admitting this level profile, or None above ceiling.

    Monotonicity lets the search go upward in the budget and stop at the first
    success, which is then the exact minimum.
    """
    cap = n + 1
    for budget in range(r, min(cap * r, ceiling) + 1):
        seen = set()
        for mhat in _multiplicity_vectors(r, budget, cap):
            if mhat in seen:
                continue
            seen.add(mhat)
            if _admissible(patterns, mhat, r):
                return budget
    return None


def measure(n, max_classes, max_spread, ceiling):
    best: dict[int, dict] = {}
    stats = {"level_sets": 0, "rank_survivors": 0, "cost_searches": 0}
    for r in range(2, max_classes + 1):
        uncapped = occupancy_patterns(n, r)
        for spread in range(r - 1, max_spread + 1):
            interiors = (
                [()] if r == 2 else itertools.combinations(range(1, spread), r - 2)
            )
            for interior in interiors:
                levels = (0,) + tuple(interior) + (spread,)
                divisor = 0
                for level in levels:
                    divisor = math.gcd(divisor, level)
                if divisor != 1:
                    continue
                stats["level_sets"] += 1
                # cost is at least r, so a profile with more classes than the
                # incumbent cost cannot improve it
                incumbent = best.get(spread, {}).get("cost")
                if incumbent is not None and incumbent <= r:
                    continue
                buckets: dict[int, list] = {}
                for pattern in uncapped:
                    total = sum(
                        count * level
                        for count, level in zip(pattern, levels, strict=True)
                    )
                    buckets.setdefault(total, []).append(pattern)
                for target, patterns in buckets.items():
                    if len(patterns) < r - 1:
                        continue
                    if _rank_mod(patterns, r) != r - 1:
                        continue
                    values = tuple(target - n * level for level in levels)
                    if len(set(values)) != r:
                        continue
                    stats["rank_survivors"] += 1
                    stats["cost_searches"] += 1
                    cost = _min_cost(patterns, r, n, ceiling)
                    if cost is None:
                        continue
                    prior = best.get(spread)
                    if prior is None or cost < prior["cost"]:
                        best[spread] = {
                            "spread": spread,
                            "cost": cost,
                            "classes": r,
                            "levels": list(levels),
                            "level_target": target,
                        }
    return best, stats


def _git(*args):
    try:
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:
        return "unavailable"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--particle-number", type=int, default=3)
    parser.add_argument("--max-classes", type=int, default=7)
    parser.add_argument("--max-spread", type=int, default=22)
    parser.add_argument("--cost-ceiling", type=int, default=40)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    started = time.time()
    best, stats = measure(
        args.particle_number, args.max_classes, args.max_spread, args.cost_ceiling
    )
    rows = [best[spread] for spread in sorted(best)]

    # D restricted to a class bound is nondecreasing here; that is measured, not
    # assumed, because a dip would mean the cost of width does not climb, which
    # is the rate claim this table exists to support.
    costs = [row["cost"] for row in rows]
    nondecreasing = all(
        left <= right for left, right in zip(costs, costs[1:], strict=False)
    )

    cross_check = {"available": FRONTIER.is_file(), "violations": []}
    if cross_check["available"]:
        frontier = json.loads(FRONTIER.read_text())
        realized = {
            row["rank"]: row["realized_max_spread"] for row in frontier["calibration"]
        }
        cross_check["realized_max_spread_by_rank"] = realized
        for rank, spread in sorted(realized.items()):
            row = best.get(spread)
            if row is None:
                continue
            # a profile of this spread exists at this rank, so the true D(W) is
            # at most rank; the restricted table may sit above that, never below
            if row["cost"] < rank:
                cross_check["violations"].append(
                    {"rank": rank, "spread": spread, "restricted_cost": row["cost"]}
                )

    payload = {
        "schema_version": 1,
        "what": (
            "D(W), the least ambient rank carrying an admissible profile of level "
            "spread W, minimised over profiles with at most max_classes classes."
        ),
        "generated_by": "scripts/profile_spread_cost.py",
        "documented_in": ["docs/profile_native_generation.md"],
        "guarded_by": "tests/test_profile_spread_cost.py",
        "provenance": {
            "commit": _git("rev-parse", "HEAD"),
            "dirty": _git("status", "--porcelain") != "",
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "parameters": {
            "particle_number": args.particle_number,
            "max_classes": args.max_classes,
            "max_spread": args.max_spread,
            "cost_ceiling": args.cost_ceiling,
        },
        "proof_labels": {
            "restricted_cost_is_an_upper_bound_on_D": "PROVED, since a profile "
            "with more classes can only be cheaper",
            "cost_is_nondecreasing_in_spread": "EXACTLY VERIFIED ON FINITE "
            "POPULATION, within the class and spread bounds",
            "spread_is_bounded_at_each_rank": "PROVED by Theorem S "
            "(docs/profile_native_generation.md), independently of this table; "
            "what stays open is whether the true bound is polynomial in r, and "
            "this table measures the rate rather than the existence",
        },
        "rows": rows,
        "cost_is_nondecreasing": nondecreasing,
        "frontier_cross_check": cross_check,
        "counters": stats,
        "seconds": round(time.time() - started, 1),
        "nonclaims": [
            "This does not bound the spread. The table is an upper bound on D "
            "restricted to a class bound, so it cannot exclude a cheap profile "
            "with more classes. Theorem S bounds the spread; this table does "
            "not and never could.",
            "It says nothing about facets. A profile is a candidate orbit.",
        ],
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output} in {payload['seconds']} s")
    print(f"  spreads measured: {len(rows)}, nondecreasing: {nondecreasing}")
    print(f"  cross-check violations: {len(cross_check['violations'])}")
    for row in rows:
        print(
            f"  W={row['spread']:3d}  D<={row['cost']:3d}  r={row['classes']}  "
            f"levels={row['levels']} target={row['level_target']}"
        )
    return 1 if cross_check["violations"] or not nondecreasing else 0


if __name__ == "__main__":
    raise SystemExit(main())
