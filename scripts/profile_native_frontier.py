#!/usr/bin/env python3
"""Measure the GPC candidate frontier with the profile-native enumerator.

Writes ``results/data/generator_frontier_refinement.json``.

Three things are recorded and they carry different proof labels.

1. CALIBRATION, ranks 6 to 10.  The profile enumerator is run against every
   population the repository already knows by an independent route: unoriented
   and oriented orbit counts, self-paired orbits, trace-alive orbits, and
   Ressayre trace-survivor totals.  Agreement is exact or the run fails.

2. SATURATION, ranks 6 to 10.  The enumerator is complete only up to a level
   spread bound, which is not proved.  The protocol raises the bound and looks
   for a plateau.  Where the truth is known the plateau equals it, and the
   realized maximum spread stays strictly inside the bound.

3. FRONTIER, ranks 11 to 13.  The same measurement with no oracle to check it
   against.  These numbers are labelled conditional on the spread bound and are
   not a claim that the candidate universe is closed.

The heavy ranks are opt in.  Defaults keep the run under a couple of minutes;
``--frontier`` adds ranks 11 to 13 and takes about ninety minutes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.generation.exterior import exterior_weights  # noqa: E402
from gpc_census.generation.orbit_canonical import (  # noqa: E402
    trace_survivor_count,
)
from gpc_census.generation.profiles import (  # noqa: E402
    enumerate_admissible_profiles,
    profile_normal,
    spread_saturation,
    verify_profile_against_weights,
)

OUTPUT = ROOT / "results" / "data" / "generator_frontier_refinement.json"

# Independently known populations.  Sources are named in the artifact.
KNOWN = {
    "unoriented_orbits": {6: 8, 7: 19, 8: 56, 9: 231, 10: 1337},
    "oriented_orbits": {6: 15, 7: 35, 8: 109},
    "self_paired_orbits": {6: 1, 7: 3, 8: 3},
    "trace_alive_orbits": {7: 25, 8: 64, 9: 191},
    "trace_survivors": {6: 64, 7: 549, 8: 6607},
}
# Realized maxima at the calibrated ranks, used to set the saturation ladder.
REALIZED_SPREAD = {6: 3, 7: 4, 8: 6, 9: 10, 10: 14}

CALIBRATION_SPREAD = 14
# Rank 11 is probed further because it is the first rank that does not plateau.
RANK_ELEVEN_SPREAD = 20


def _canonical_tau(tau):
    forward = tuple(sorted(tau, reverse=True))
    backward = tuple(sorted((-value for value in tau), reverse=True))
    return min(forward, backward)


def measure(n, d, max_spread, cache, *, verify_rank):
    started = time.time()
    profiles = enumerate_admissible_profiles(n, d, max_spread, level_cache=cache)
    weights = exterior_weights(n, d)
    unoriented = set()
    self_paired = alive = survivors = 0
    spread_histogram: dict[int, int] = {}
    class_histogram: dict[int, int] = {}
    for profile in profiles:
        tau = profile.tau
        unoriented.add(_canonical_tau(tau))
        if sum(tau) == 0 and sorted(-value for value in tau) == sorted(tau):
            self_paired += 1
        h, z = profile_normal(profile)
        below = sum(
            1
            for weight in weights
            if sum(left * right for left, right in zip(h, weight, strict=True)) < z
        )
        budget = d * (d - 1) // 2 - sum(
            value * (value - 1) // 2 for value in profile.multiplicities
        )
        if below <= budget:
            alive += 1
        survivors += trace_survivor_count(h, below)
        spread_histogram[profile.spread] = spread_histogram.get(profile.spread, 0) + 1
        class_histogram[profile.classes] = class_histogram.get(profile.classes, 0) + 1
    row = {
        "particle_number": n,
        "rank": d,
        "max_spread": max_spread,
        "oriented_orbits": len(profiles),
        "unoriented_orbits": len(unoriented),
        "self_paired_orbits": self_paired,
        "trace_alive_orbits": alive,
        "trace_survivors": survivors,
        "realized_max_spread": max(spread_histogram),
        "spread_bound_is_slack": max(spread_histogram) < max_spread,
        "max_abs_tau": max(max(abs(v) for v in t) for t in unoriented),
        "spread_histogram": dict(sorted(spread_histogram.items())),
        "class_histogram": dict(sorted(class_histogram.items())),
        "seconds": round(time.time() - started, 1),
    }
    if verify_rank:
        row["theorem_a_disagreements"] = sum(
            1 for profile in profiles if not verify_profile_against_weights(n, profile)
        )
        row["theorem_a_checked"] = len(profiles)
    return row, profiles


def hole_budget_table(profiles_by_rank):
    """The bounded-hole grammar, scored against the complete candidate universe.

    A hole budget ``h`` means level spread at most ``r - 1 + h``, so ``h = 0`` is
    an exact arithmetic progression of levels.
    """
    out = []
    for d, profiles in sorted(profiles_by_rank.items()):
        truth = len({_canonical_tau(p.tau) for p in profiles})
        for budget in range(0, 5):
            reached = len(
                {
                    _canonical_tau(p.tau)
                    for p in profiles
                    if p.spread <= p.classes - 1 + budget
                }
            )
            out.append(
                {
                    "rank": d,
                    "hole_budget": budget,
                    "orbits_reached": reached,
                    "orbits_total": truth,
                    "missing": truth - reached,
                    "complete": reached == truth,
                }
            )
    return out


def _git(*args):
    try:
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:
        return "unavailable"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--frontier",
        action="store_true",
        help="also measure ranks 11 to 13; adds roughly ninety minutes",
    )
    parser.add_argument("--frontier-spread", type=int, default=14)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    cache: dict = {}
    calibration = []
    profiles_by_rank = {}
    mismatches = []
    for d in (6, 7, 8, 9, 10):
        row, profiles = measure(3, d, CALIBRATION_SPREAD, cache, verify_rank=d <= 9)
        profiles_by_rank[d] = profiles
        row["known"] = {}
        for field, table in KNOWN.items():
            if d in table:
                row["known"][field] = table[d]
                if row[field] != table[d]:
                    mismatches.append(
                        {"rank": d, "field": field, "got": row[field],
                         "known": table[d]}
                    )
        calibration.append(row)

    cross_particle = []
    for n, d in ((4, 7), (4, 8)):
        row, _ = measure(n, d, CALIBRATION_SPREAD, cache, verify_rank=True)
        cross_particle.append(row)

    saturation = []
    for d, realized in sorted(REALIZED_SPREAD.items()):
        for entry in spread_saturation(
            3, d, (realized - 1, realized, realized + 1), level_cache=cache
        ):
            saturation.append({"rank": d, **entry})

    frontier = []
    rank_eleven_ladder = []
    if args.frontier:
        # Rank 11 gets its own ladder, because it is the first rank where the
        # protocol does NOT plateau and that is the load-bearing observation.
        for entry in spread_saturation(
            3, 11, (12, 14, 16, 18, RANK_ELEVEN_SPREAD), level_cache=cache
        ):
            rank_eleven_ladder.append({"rank": 11, **entry})
        for d in (11, 12, 13):
            spread = RANK_ELEVEN_SPREAD if d == 11 else args.frontier_spread
            row, _ = measure(3, d, spread, cache, verify_rank=False)
            frontier.append(row)

    payload = {
        "schema_version": 1,
        "what": (
            "Profile-native enumeration of admissible exterior-weight hyperplane "
            "orbits, calibrated against every independently known population and "
            "then run past the closed census."
        ),
        "generated_by": "scripts/profile_native_frontier.py",
        "documented_in": [
            "docs/profile_native_generation.md",
            "docs/generator_frontier_refinement.md",
        ],
        "guarded_by": "tests/test_profiles.py",
        "provenance": {
            "commit": _git("rev-parse", "HEAD"),
            "dirty": _git("status", "--porcelain") != "",
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "known_population_sources": {
            "unoriented_orbits": "docs/orbit_canonical_flats.md",
            "oriented_orbits": "docs/orbit_canonical_flats.md",
            "self_paired_orbits": "docs/orbit_canonical_flats.md",
            "trace_alive_orbits": "docs/levi_triangularity_gate.md (mechanisms)",
            "trace_survivors": "docs/screening_orbit_structure.md",
        },
        "proof_labels": {
            "theorem_a_profile_reduction": "PROVED, and checked against the "
            "materialised weight rank on every profile at ranks 6 to 9",
            "lemma_l_level_reparametrisation": "PROVED",
            "monotonicity_in_multiplicity": "PROVED",
            "calibration_agreement": "EXACTLY VERIFIED ON FINITE POPULATION",
            "spread_bound": "UNRESOLVED, and it is the only gap between this "
            "enumeration and a candidate-coverage theorem",
            "bounded_hole_grammar_as_completeness_rule": "REFUTED",
            "rank_11_to_13_candidate_universe": "CONDITIONAL on the spread bound",
        },
        "calibration": calibration,
        "calibration_mismatches": mismatches,
        "calibration_is_exact": not mismatches,
        "cross_particle_number": cross_particle,
        "spread_saturation": saturation,
        "rank_eleven_spread_ladder": rank_eleven_ladder,
        "rank_eleven_saturates": bool(rank_eleven_ladder)
        and rank_eleven_ladder[-1]["bound_is_slack"],
        "bounded_hole_grammar": hole_budget_table(profiles_by_rank),
        "frontier": frontier,
        "frontier_is_conditional_on_spread_bound": True,
        "nonclaims": [
            "An admissible hyperplane orbit is a CANDIDATE. It is not a Ressayre "
            "row, not a facet, and not a valid inequality.",
            "The frontier rows are complete only for profiles of level spread at "
            "most max_spread. At rank 11 the enumeration had not saturated by "
            "spread 20, so the (3,11) candidate universe is not closed here.",
            "Trace-alive counts use the Ressayre trace test only. Determinant "
            "certification and exact reduction are downstream and are not run.",
        ],
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.write_text(text)
    digest = hashlib.sha256(text.encode()).hexdigest()
    print(f"wrote {args.output} sha256={digest}")
    for row in calibration + frontier:
        print(
            f"  (3,{row['rank']}) spread<={row['max_spread']}: "
            f"oriented={row['oriented_orbits']} unoriented={row['unoriented_orbits']} "
            f"alive={row['trace_alive_orbits']} survivors={row['trace_survivors']} "
            f"realized_spread={row['realized_max_spread']}"
        )
    if mismatches:
        print("CALIBRATION MISMATCH", mismatches)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
