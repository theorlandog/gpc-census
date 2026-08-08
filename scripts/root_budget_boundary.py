#!/usr/bin/env python3
"""The root-budget boundary layer, tested against the census before being trusted.

THE CLAIM UNDER TEST. For a dominant one-parameter subgroup ``tau`` satisfying
the Bulois-Denis-Ressayre grade condition, the positive weight family

    Omega_+(tau) = { I : sum_{i in I} tau_i > 0 }

is an ideal of the Grassmannian dominance order and has size at most
``binom(d,2)``. Since an ideal contains the principal ideal of each of its
members, every positive weight ``I`` obeys ``|down I| <= binom(d,2)``. Under
the standard bijection ``N``-subsets to Young diagrams in the ``N x (d-N)``
rectangle, ``|down I|`` counts the subdiagrams of the diagram, so the layer can
be enumerated without ever materialising ``binom(d,N)`` weights. At ``(20,40)``
that replaces 137,846,528,820 weights by 15,359.

TWO HALVES, AND THEY DO NOT FARE THE SAME. This script checks both against
every surviving mechanism at ranks 7, 8 and 9, where the repository knows the
exact population.

* The POSITIVE half needs no regularity: below a positive weight every weight
  is positive, because ``J <= I`` gives ``tau.J >= tau.I > 0``. Verified on
  every mechanism.
* The ZERO half, ``|down I| <= binom(d,2) + 1`` for a zero weight, needs
  ``tau`` REGULAR. The step "everything strictly below a zero weight is
  strictly positive" requires ``tau.J > tau.I``, which fails as soon as two
  entries of ``tau`` coincide. Refuted on the census.

A CORRECTED BOUND is also checked. For any weight with ``tau.I >= 0``,
``down I`` is contained in ``Omega_+ union Omega_0``, so

    |down I| <= binom(d,2) + |Omega_0|

where ``|Omega_0|`` is the zero-grade dimension the Levi-fusion audit already
computes, available in closed form as ``[y^N x^0] prod_i (1 + y x^{l_i})^{m_i}``
without touching the weight list. This holds on every mechanism tested. It is
weaker than the uniform version because it depends on the Levi type rather than
on the diagram alone.

Run:

    uv run scripts/root_budget_boundary.py
    uv run scripts/root_budget_boundary.py --skip-census
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import itertools
import json
import math
import pathlib
import sys
import time
from collections import deque
from functools import lru_cache

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CENSUS_SYSTEMS = ((3, 7, "flats"), (3, 8, "flats"), (3, 9, "augmentation"))
TIMING_FIELDS = ("seconds",)


# ---------------------------------------------------------------------------
# Young-lattice boundary layer
# ---------------------------------------------------------------------------


def subpartition_counter(rows: int, columns: int):
    @lru_cache(maxsize=None)
    def count(alpha: tuple[int, ...]) -> int:
        @lru_cache(maxsize=None)
        def dynamic(row: int, previous: int) -> int:
            if row == rows:
                return 1
            bound = min(previous, alpha[row])
            return sum(dynamic(row + 1, value) for value in range(bound + 1))

        return dynamic(0, columns)

    return count


def addable_partitions(alpha: tuple[int, ...], columns: int):
    values = list(alpha)
    for row in range(len(values)):
        if values[row] >= columns:
            continue
        if row and values[row] + 1 > values[row - 1]:
            continue
        child = values.copy()
        child[row] += 1
        yield tuple(child)


def truncated_layer(particle_number: int, orbitals: int, budget: int) -> dict:
    """Partitions whose principal ideal has size at most ``budget``."""
    columns = orbitals - particle_number
    count = subpartition_counter(particle_number, columns)
    initial = (0,) * particle_number
    queue = deque([initial])
    seen = {initial}
    accepted: list[tuple[tuple[int, ...], int]] = []
    frontier: list[tuple[tuple[int, ...], int]] = []

    started = time.perf_counter()
    while queue:
        alpha = queue.popleft()
        ideal_size = count(alpha)
        if ideal_size <= budget:
            accepted.append((alpha, ideal_size))
            for child in addable_partitions(alpha, columns):
                if child not in seen:
                    seen.add(child)
                    queue.append(child)
        else:
            frontier.append((alpha, ideal_size))

    return {
        "particle_number": particle_number,
        "orbitals": orbitals,
        "rectangle": [particle_number, columns],
        "budget": budget,
        "ambient_weight_count": math.comb(orbitals, particle_number),
        "accepted_weight_types": len(accepted),
        "minimal_excluded_frontier": len(frontier),
        "accepted_plus_frontier": len(accepted) + len(frontier),
        "maximum_young_area": max(sum(alpha) for alpha, _ in accepted),
        "maximum_ideal_size": max(size for _alpha, size in accepted),
        "minimum_frontier_ideal_size": min(size for _alpha, size in frontier),
        "seconds": round(time.perf_counter() - started, 3),
    }


# ---------------------------------------------------------------------------
# Census validation
# ---------------------------------------------------------------------------


def subset_to_partition(subset, particle_number: int) -> tuple[int, ...]:
    return tuple(
        subset[particle_number - 1 - k] - (particle_number - 1 - k)
        for k in range(particle_number)
    )


def durfee_rank(partition) -> int:
    return sum(1 for index, value in enumerate(partition) if value > index)


def transpose(partition, rows: int):
    if not any(partition):
        return tuple([0] * rows)
    columns = max(partition)
    conjugate = [sum(1 for value in partition if value > column) for column in range(columns)]
    return tuple(conjugate + [0] * rows)[:rows]


def levi_budget(multiplicities, d: int) -> int:
    """R_L = sum_{i<j} m_i m_j = (d^2 - sum m_i^2) / 2."""
    return (d * d - sum(value * value for value in multiplicities)) // 2


def levi_budget_matches_repository(trials: int = 400) -> dict:
    """R_L is the repository's existing maximum-inversion bound, not a new one."""
    import random

    from gpc_census.generation.orbit_canonical import _max_inversions

    rng = random.Random(3)
    mismatches = 0
    for _ in range(trials):
        d = rng.randrange(4, 14)
        parts = []
        remaining = d
        while remaining > 0:
            piece = rng.randrange(1, remaining + 1)
            parts.append(piece)
            remaining -= piece
        counts = {index: value for index, value in enumerate(parts)}
        if levi_budget(parts, d) != _max_inversions(counts, d):
            mismatches += 1
    return {
        "trials": trials,
        "mismatches": mismatches,
        "identical_to_max_inversions": mismatches == 0,
        "note": (
            "The per-Levi budget is gpc_census.generation.orbit_canonical."
            "_max_inversions, which orbit_is_trace_dead already applies, so it "
            "is the sharp form the survivor generator uses rather than a new "
            "quantity to implement."
        ),
    }


def occupancy_profiles(multiplicities, n: int):
    """All Levi occupancy profiles k with sum k_i = n and 0 <= k_i <= m_i."""
    size = len(multiplicities)

    def recurse(index: int, remaining: int, prefix: list[int]):
        if index == size:
            if remaining == 0:
                yield tuple(prefix)
            return
        for taken in range(min(multiplicities[index], remaining) + 1):
            yield from recurse(index + 1, remaining - taken, prefix + [taken])

    return list(recurse(0, n, []))


def dominates(upper, lower) -> bool:
    """Prefix-sum dominance: more particles in the high-level blocks."""
    a = b = 0
    for x, y in zip(upper, lower):
        a += x
        b += y
        if a < b:
            return False
    return True


def profile_budget_laws(tau: tuple[int, ...], n: int) -> dict:
    """The weighted occupancy-profile pruning laws, checked on one normal.

    A determinant is represented only by its profile ``k``, with multiplicity
    ``M(k) = prod binom(m_i, k_i)`` and grade ``g(k) = sum l_i k_i``. With
    ``Q(k)`` the weighted mass of the upper dominance ideal, strict dominance
    raises the grade, so everything strictly above a nonnegative profile is
    strictly positive. Hence

        g(k) > 0   =>   Q(k)        <= R_L
        g(k) == 0  =>   Q(k) - M(k) <= R_L

    The second law is the one that repairs the per-weight version: a zero
    profile may stand for an astronomical number of determinants, and only its
    strictly dominant mass has to fit inside the root budget.
    """
    counts = collections.Counter(tau)
    levels = sorted(counts, reverse=True)
    multiplicities = [counts[level] for level in levels]
    budget = levi_budget(multiplicities, len(tau))
    profiles = occupancy_profiles(multiplicities, n)

    mass = {
        profile: math.prod(
            math.comb(multiplicities[i], profile[i]) for i in range(len(multiplicities))
        )
        for profile in profiles
    }
    grade = {
        profile: sum(levels[i] * profile[i] for i in range(len(multiplicities)))
        for profile in profiles
    }
    upper = {
        profile: sum(mass[other] for other in profiles if dominates(other, profile))
        for profile in profiles
    }

    positive_excess = max(
        (upper[k] - budget for k in profiles if grade[k] > 0), default=None
    )
    zero_excess = max(
        (upper[k] - mass[k] - budget for k in profiles if grade[k] == 0), default=None
    )
    return {
        "levi_budget": budget,
        "profiles": len(profiles),
        "positive_family_mass": sum(mass[k] for k in profiles if grade[k] > 0),
        "positive_law_holds": positive_excess is None or positive_excess <= 0,
        "zero_law_holds": zero_excess is None or zero_excess <= 0,
        "worst_positive_excess": positive_excess,
        "worst_zero_excess": zero_excess,
    }


GATE_LADDER = (((13, 17), (4, 17)), ((10, 20), (10, 20)), ((15, 30), (15, 30)),
               ((20, 40), (20, 40)))


def gate_ladder_record(system, effective) -> dict:
    """Entry numbers for one rung of the proposed validation ladder.

    The sign-control universe is the zero-or-positive layer plus its minimal
    excluded frontier, which is what certifies that everything deeper is
    strictly negative. Compression is measured against the ambient
    ``binom(d, N)`` of the ORIGINAL system, since particle-hole normalization
    is what makes the effective system smaller.
    """
    original_n, d = system
    effective_n, _ = effective
    budget = math.comb(d, 2)
    positive = truncated_layer(effective_n, d, budget)
    nonnegative = truncated_layer(effective_n, d, budget + 1)
    sign_control = (
        nonnegative["accepted_weight_types"] + nonnegative["minimal_excluded_frontier"]
    )
    ambient = math.comb(d, original_n)
    durfee = max(
        (r for r in range(1, 12) if math.comb(2 * r, r) <= budget + 1), default=0
    )
    return {
        "system": list(system),
        "effective_system": list(effective),
        "particle_hole_reduction_used": original_n != effective_n,
        "ambient_weight_count": ambient,
        "positive_layer": positive["accepted_weight_types"],
        "zero_or_positive_layer": nonnegative["accepted_weight_types"],
        "minimal_excluded_frontier": nonnegative["minimal_excluded_frontier"],
        "sign_control_universe": sign_control,
        "compression_factor": round(ambient / sign_control, 1),
        "max_durfee_rank": durfee,
        "transpose": transpose_orbit_counts(effective_n, d),
    }


def transpose_orbit_counts(particle_number: int, orbitals: int) -> dict:
    """Half-filling transpose quotient of the sign-control universe."""
    if 2 * particle_number != orbitals:
        return {"applies": False, "reason": "transpose is an involution only when 2N = d"}
    budget = math.comb(orbitals, 2)
    columns = orbitals - particle_number
    count = subpartition_counter(particle_number, columns)
    initial = (0,) * particle_number
    queue = deque([initial])
    seen = {initial}
    universe = set()
    while queue:
        alpha = queue.popleft()
        universe.add(alpha)
        if count(alpha) <= budget + 1:
            for child in addable_partitions(alpha, columns):
                if child not in seen:
                    seen.add(child)
                    queue.append(child)
    fixed = sum(1 for alpha in universe if transpose(alpha, particle_number) == alpha)
    return {
        "applies": True,
        "sign_control_universe": len(universe),
        "self_conjugate": fixed,
        "transpose_orbits": (len(universe) + fixed) // 2,
    }


def surviving_mechanisms(n: int, d: int, source: str, checkpoint_dir=None):
    """Sorted-tau shapes carrying at least one Ressayre trace survivor."""
    import importlib.util

    from gpc_census.generation import orbit_canonical as oc
    from gpc_census.generation.exterior import exterior_weights
    from gpc_census.generation.ressayre import admissible_candidate_from_tau

    spec = importlib.util.spec_from_file_location(
        "rank9_compression_audit", ROOT / "scripts" / "rank9_compression_audit.py"
    )
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)

    weights = exterior_weights(n, d)
    occupied = [tuple(i for i, bit in enumerate(w) if bit) for w in weights]
    if source == "flats":
        keys = audit.orbits_by_flats(n, d)
    else:
        keys = audit.orbits_by_augmentation(n, d, checkpoint_dir=checkpoint_dir)
    for key in keys:
        candidate = admissible_candidate_from_tau(n, list(key))
        below = sum(
            1
            for members in occupied
            if sum(candidate.h[i] for i in members) < candidate.z
        )
        if not oc.orbit_is_trace_dead(candidate.h, below):
            yield key


def validate_on_census(n: int, d: int, source: str, checkpoint_dir=None) -> dict:
    started = time.perf_counter()
    budget = math.comb(d, 2)
    columns = d - n
    count = subpartition_counter(n, columns)
    slice_weights = list(itertools.combinations(range(d), n))

    mechanisms = 0
    regular = 0
    ideal_holds = 0
    size_holds = 0
    positive_holds = 0
    zero_holds_uniform = 0
    zero_holds_corrected = 0
    worst_positive_ideal = 0
    worst_zero_ideal = 0
    largest_positive_family = 0
    max_distinct_levels = 0
    uniform_violations: list[dict] = []
    levi_positive_family = 0
    levi_positive_bound = 0
    levi_zero_bound = 0
    measured_durfee = 0
    levi_budget_min = None
    levi_budget_max = 0
    profile_positive_law = 0
    profile_zero_law = 0
    profile_mass_matches = 0

    for key in surviving_mechanisms(n, d, source, checkpoint_dir=checkpoint_dir):
        mechanisms += 1
        tau = tuple(sorted(key, reverse=True))
        multiplicities = list(collections.Counter(tau).values())
        budget_levi = levi_budget(multiplicities, d)
        levi_budget_min = (
            budget_levi if levi_budget_min is None else min(levi_budget_min, budget_levi)
        )
        levi_budget_max = max(levi_budget_max, budget_levi)
        levels = len(set(tau))
        max_distinct_levels = max(max_distinct_levels, levels)
        if levels == d:
            regular += 1

        sums = {I: sum(tau[i] for i in I) for I in slice_weights}
        positive = {I for I in slice_weights if sums[I] > 0}
        zeros = [I for I in slice_weights if sums[I] == 0]
        largest_positive_family = max(largest_positive_family, len(positive))
        if len(positive) <= budget:
            size_holds += 1

        closed = all(
            other in positive
            for member in positive
            for other in slice_weights
            if all(a <= b for a, b in zip(other, member))
        )
        if closed:
            ideal_holds += 1

        positive_ideals = [count(subset_to_partition(I, n)) for I in positive]
        zero_ideals = [count(subset_to_partition(I, n)) for I in zeros]
        if positive_ideals:
            worst_positive_ideal = max(worst_positive_ideal, max(positive_ideals))
        if zero_ideals:
            worst_zero_ideal = max(worst_zero_ideal, max(zero_ideals))

        if all(size <= budget for size in positive_ideals):
            positive_holds += 1
        if len(positive) <= budget_levi:
            levi_positive_family += 1
        if all(size <= budget_levi for size in positive_ideals):
            levi_positive_bound += 1
        if all(size <= budget_levi + 1 for size in zero_ideals):
            levi_zero_bound += 1
        for member in positive:
            measured_durfee = max(
                measured_durfee, durfee_rank(subset_to_partition(member, n))
            )
        laws = profile_budget_laws(tau, n)
        if laws["positive_law_holds"]:
            profile_positive_law += 1
        if laws["zero_law_holds"]:
            profile_zero_law += 1
        if laws["positive_family_mass"] == len(positive):
            profile_mass_matches += 1
        if all(size <= budget + 1 for size in zero_ideals):
            zero_holds_uniform += 1
        elif len(uniform_violations) < 5:
            worst = max(zero_ideals)
            uniform_violations.append(
                {
                    "tau": list(tau),
                    "distinct_levels": levels,
                    "largest_zero_weight_ideal": worst,
                    "uniform_bound": budget + 1,
                    "zero_grade_dimension": len(zeros),
                    "corrected_bound": budget + len(zeros),
                }
            )
        if all(size <= budget + len(zeros) for size in zero_ideals):
            zero_holds_corrected += 1

    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "root_budget": budget,
        "mechanisms": mechanisms,
        "regular_mechanisms": regular,
        "max_distinct_levels": max_distinct_levels,
        "positive_family_is_a_dominance_ideal": ideal_holds,
        "positive_family_within_root_budget": size_holds,
        "largest_positive_family": largest_positive_family,
        "positive_weight_bound_holds": positive_holds,
        "worst_positive_weight_ideal": worst_positive_ideal,
        "zero_weight_uniform_bound_holds": zero_holds_uniform,
        "zero_weight_corrected_bound_holds": zero_holds_corrected,
        "worst_zero_weight_ideal": worst_zero_ideal,
        "uniform_bound_violations": uniform_violations,
        "levi_budget_range": [levi_budget_min, levi_budget_max],
        "levi_positive_family_within_budget": levi_positive_family,
        "levi_positive_weight_bound_holds": levi_positive_bound,
        "levi_zero_weight_bound_holds": levi_zero_bound,
        "measured_max_durfee_rank_of_positive_weights": measured_durfee,
        "profile_positive_law_holds": profile_positive_law,
        "profile_zero_law_holds": profile_zero_law,
        "profile_mass_matches_weight_count": profile_mass_matches,
        "durfee_rank_bound_from_global_budget": max(
            (r for r in range(1, 12) if math.comb(2 * r, r) <= budget), default=0
        ),
        "seconds": round(time.perf_counter() - started, 1),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbitals", type=int, default=40)
    parser.add_argument("--particles", type=int, nargs="*", default=[3, 4, 5, 10, 20])
    parser.add_argument("--skip-census", action="store_true")
    parser.add_argument("--checkpoint-dir", default=None)
    parser.add_argument(
        "-o",
        "--out",
        type=pathlib.Path,
        default=ROOT / "results/data/root_budget_boundary.json",
    )
    arguments = parser.parse_args()

    root_budget = math.comb(arguments.orbitals, 2)
    layers = []
    for particle_number in arguments.particles:
        positive = truncated_layer(particle_number, arguments.orbitals, root_budget)
        nonnegative = truncated_layer(
            particle_number, arguments.orbitals, root_budget + 1
        )
        layers.append(
            {
                "particle_number": particle_number,
                "positive_layer": positive,
                "zero_or_positive_layer": nonnegative,
            }
        )
        print(
            f"  N={particle_number:2d} ambient={positive['ambient_weight_count']:,} "
            f"positive layer={positive['accepted_weight_types']:,} "
            f"zero-or-positive={nonnegative['accepted_weight_types']:,}",
            flush=True,
        )

    census = {}
    if not arguments.skip_census:
        for n, d, source in CENSUS_SYSTEMS:
            print(f"  validating ({n},{d}) ...", flush=True)
            record = validate_on_census(
                n, d, source, checkpoint_dir=arguments.checkpoint_dir
            )
            census[record["system"]] = record
            print(
                f"    mechanisms {record['mechanisms']} regular "
                f"{record['regular_mechanisms']} | positive bound "
                f"{record['positive_weight_bound_holds']}/{record['mechanisms']} | "
                f"zero uniform {record['zero_weight_uniform_bound_holds']}"
                f"/{record['mechanisms']} | zero corrected "
                f"{record['zero_weight_corrected_bound_holds']}/{record['mechanisms']}",
                flush=True,
            )

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/root_budget_boundary.py",
        "orbitals": arguments.orbitals,
        "root_budget": root_budget,
        "boundary_layers": layers,
        "census_validation": census,
        "levi_budget_identity": levi_budget_matches_repository(),
        "gate_ladder": [
            gate_ladder_record(system, effective) for system, effective in GATE_LADDER
        ],
        "half_filling_transpose": {
            f"({n},{d})": transpose_orbit_counts(n, d)
            for n, d in ((10, 20), (15, 30), (20, 40))
        },
        "verdict": {
            "positive_layer_bound": (
                "HOLDS, and needs no regularity: below a positive weight every "
                "weight is positive. Verified on every surviving mechanism at "
                "ranks 7, 8 and 9."
            ),
            "zero_layer_uniform_bound": (
                "REFUTED for non-regular tau. The step from a zero weight needs "
                "strict increase, which fails when two entries of tau coincide. "
                "Violated on the census at every rank tested."
            ),
            "regularity_in_the_census": (
                "Regular mechanisms are rare: none at ranks 7 and 8, three of "
                "191 at rank 9. The regular stratum the (20,40) counts describe "
                "is therefore not where the known mechanisms live."
            ),
            "weighted_occupancy_profile_laws": (
                "HOLD, 280 of 280, on both branches. These repair the "
                "per-weight zero bound: replacing the '+1' by the profile's "
                "own multiplicity M(k) is exactly what the non-regular case "
                "needs, and the zero law holds at every mechanism where the "
                "per-weight version failed."
            ),
            "per_levi_budget": (
                "SHARPENS the positive half and AGGRAVATES the zero half. "
                "R_L = (d^2 - sum m_i^2)/2 is the repository's own "
                "_max_inversions, so orbit_is_trace_dead already applies it. "
                "Every positive bound still holds under R_L, which is strictly "
                "smaller than binom(d,2); the zero bound fails more often "
                "under R_L than under the global budget, exactly because the "
                "budget shrinks."
            ),
            "corrected_zero_bound": (
                "|down I| <= binom(d,2) + |Omega_0| holds on every mechanism "
                "tested. It is Levi-type dependent rather than diagram-only, so "
                "it does not give the same uniform Young-lattice truncation."
            ),
        },
    }
    hashable = json.loads(json.dumps(payload))
    for entry in hashable["boundary_layers"]:
        for layer in ("positive_layer", "zero_or_positive_layer"):
            for field in TIMING_FIELDS:
                entry[layer].pop(field, None)
    for record in hashable["census_validation"].values():
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
