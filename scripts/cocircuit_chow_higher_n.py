#!/usr/bin/env python3
"""The cocircuit route and the signed-Chow decoder at arbitrary particle number.

`docs/cocircuit_chow_decoder.md` proved C1, T4 and T5 and gated them at `N=3`.
None of the three proofs uses `N=3`: C1 needs only that every exterior weight
has coordinate sum `N != 0`, and T1 to T5 are fixed-weight statements. This
script tests that the CODE inherits the generality the mathematics has, which
is a different claim and is the one that can fail.

FOUR GATES, and the interesting ones are the two that push back.

1. DUAL POPULATIONS. `(4,7)` and `(5,8)` are the particle-hole duals of `(3,7)`
   and `(3,8)`, so their admissible-hyperplane populations must coincide with
   the `N=3` ones. Gated pair for pair against
   `enumerate_admissible_hyperplanes_by_flats`, not merely by count.

2. PUBLISHED HIGHER-N ROWS. Every representative normal the Levi-fusion audits
   stored at `(4,8)`, `(4,9)`, `(4,10)` and `(5,10)` is put through the whole
   chain: zero-family rank, the trace identity, both decoded sign families, and
   ORIENTED reconstruction of the primitive normal. These are census rows, so
   this is the gate that ties the general code to published data.

3. THE PARTICLE-HOLE REDUCTION, measured rather than assumed. Decoding the
   complement shadow in the `(d-N)` layer is exact. It does NOT shrink the
   profile problem: `k -> m - k` is a bijection between the two layers'
   profiles, so the variable count is invariant. What it does change is the
   arithmetic, and the branch count is measured here rather than asserted.

4. THE WALL. Near half filling the unsymmetric reverse search stops finishing.
   `(4,8)` is recorded as a NON-COMPLETION at a stated node budget. No partial
   hyperplane count is reported, because a truncated reverse search is not a
   population.

Plus a stress that is relevant rather than generic: exact trace survivors built
from degenerate level multisets through the repository's own
`generate_trace_survivors`, instead of uniform random normals most of which no
generator would ever see.

Run:

    uv run scripts/cocircuit_chow_higher_n.py
    uv run scripts/cocircuit_chow_higher_n.py --skip-dual --skip-wall
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import pathlib
import platform
import random
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.generation import orbit_canonical as oc  # noqa: E402
from gpc_census.generation.admissible_flats import (  # noqa: E402
    enumerate_admissible_hyperplanes_by_flats,
)
from gpc_census.generation.chow_decoder import (  # noqa: E402
    annihilator,
    decode_shadow,
    decode_signed_pair,
    shadow_of,
    threshold_family,
)
from gpc_census.generation.cocircuit import (  # noqa: E402
    NodeLimitExceeded,
    enumerate_cocircuit_hyperplanes,
    trace_survivors,
)
from gpc_census.generation.exterior import exterior_weights  # noqa: E402
from gpc_census.generation.ressayre import (  # noqa: E402
    _level_sets,
    admissible_candidate_from_tau,
)

DATA = ROOT / "results" / "data"
TIMING_FIELDS = ("seconds", "cocircuit_seconds", "flat_lattice_seconds")
DUAL_SYSTEMS = ((4, 7), (5, 8))
WALL_SYSTEM = (4, 8)
WALL_NODE_BUDGET = 2_000_000
PUBLISHED_SOURCES = ("levi_fusion_pilot", "levi_fusion_heldout_results")
STRESS_PLAN = ((5, 10, 24), (6, 12, 12))
STRESS_SEED = 20260809


def orbit_key(tau):
    """Unoriented `S_d` orbit key: coordinate multiset up to global sign."""
    return min(tuple(sorted(tau)), tuple(sorted(-value for value in tau)))


def inversions(tau) -> int:
    d = len(tau)
    return sum(
        1 for i in range(d) for j in range(i + 1, d) if tau[j] > tau[i]
    )


def content_free(tau) -> tuple[int, ...]:
    """Divide out the content, KEEPING the orientation.

    The decoder returns the primitive generator of the annihilator line,
    oriented to agree with the decoded positive family. A sampled normal that
    happens to be a positive multiple of that generator, `(-4,...,6,6,6)` say,
    defines the same families and is the same candidate; comparing it to the
    decoder's output without dividing the content out would score a correct
    reconstruction as a failure. That mistake was made once here and is the
    reason this helper exists rather than an inline gcd.
    """
    divisor = 0
    for value in tau:
        divisor = _gcd(divisor, abs(value))
    return tuple(value // divisor for value in tau) if divisor else tuple(tau)


def _gcd(left: int, right: int) -> int:
    from math import gcd

    return gcd(left, right)


def dual_gate(n: int, d: int) -> dict:
    """Enumerate at `(n,d)` and gate against the flat lattice pair for pair."""
    started = time.time()
    enumeration = enumerate_cocircuit_hyperplanes(n, d)
    cocircuit_seconds = round(time.time() - started, 1)

    survivors = list(trace_survivors(enumeration))
    nonstructural = [tau for tau, positive in survivors if positive > 0]

    started = time.time()
    reference = enumerate_admissible_hyperplanes_by_flats(n, d)
    flat_seconds = round(time.time() - started, 1)

    decode_failures = 0
    tau_failures = 0
    max_profile_variables = 0
    max_states = 0
    for tau in nonstructural:
        positive = threshold_family(tau, n, +1)
        negative = threshold_family(tau, n, -1)
        decoded = decode_signed_pair(
            shadow_of(positive, d), shadow_of(negative, d), n
        )
        if decoded is None or decoded.positive.family != positive:
            decode_failures += 1
            continue
        if decoded.negative.family != negative:
            decode_failures += 1
        if decoded.tau != tuple(tau):
            tau_failures += 1
        for side in (decoded.positive, decoded.negative):
            max_profile_variables = max(
                max_profile_variables, side.profile_variables
            )
            max_states = max(max_states, side.search_states)

    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "particle_hole_dual": f"({d - n},{d})",
        "exterior_weights": enumeration.weight_count,
        "admissible_hyperplanes": enumeration.hyperplane_count,
        "reverse_search_nodes": enumeration.node_count,
        "sd_unoriented_orbits": len(
            {orbit_key(tau) for tau in enumeration.normals}
        ),
        "trace_survivors": len(survivors),
        "nonstructural_trace_survivors": len(nonstructural),
        "flat_lattice_hyperplanes": len(reference.hyperplanes),
        "agrees_with_flat_lattice_pairwise": (
            set(enumeration.affine_pairs())
            == {(row.h, row.z) for row in reference.hyperplanes}
        ),
        "decode_failures": decode_failures,
        "tau_reconstruction_failures": tau_failures,
        "max_profile_variables": max_profile_variables,
        "max_search_states_per_side": max_states,
        "cocircuit_seconds": cocircuit_seconds,
        "flat_lattice_seconds": flat_seconds,
    }


def wall_probe(n: int, d: int, node_budget: int) -> dict:
    """Record whether the unsymmetric search finishes inside a node budget.

    Deliberately records only completion or non-completion. A truncated reverse
    search has visited some subtree and not others, so its hyperplane count is
    not an estimate of anything and must never be published as one.
    """
    started = time.time()
    try:
        enumeration = enumerate_cocircuit_hyperplanes(n, d, node_limit=node_budget)
    except NodeLimitExceeded:
        return {
            "system": f"({n},{d})",
            "particle_number": n,
            "rank": d,
            "exterior_weights": len(exterior_weights(n, d)),
            "node_budget": node_budget,
            "completed": False,
            "admissible_hyperplanes": None,
            "seconds": round(time.time() - started, 1),
            "note": (
                "no partial hyperplane count is recorded: a truncated reverse "
                "search is not a population"
            ),
        }
    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "exterior_weights": enumeration.weight_count,
        "node_budget": node_budget,
        "completed": True,
        "admissible_hyperplanes": enumeration.hyperplane_count,
        "reverse_search_nodes": enumeration.node_count,
        "seconds": round(time.time() - started, 1),
    }


def _published_records(payload):
    for system, entry in payload.get("systems", {}).items():
        rows = entry.get("records") or entry.get("signature_records") or []
        taus = [
            tuple(row["tau"]) if "tau" in row else tuple(row["representative_tau"])
            for row in rows
            if row.get("tau") or row.get("representative_tau")
        ]
        if taus:
            yield system, entry["particle_number"], entry["rank"], taus


def published_gate(minimum_particle_number: int = 4) -> dict:
    """Put every published higher-N representative through the whole chain."""
    results: dict[str, dict] = {}
    for name in PUBLISHED_SOURCES:
        path = DATA / f"{name}.json"
        payload = json.loads(path.read_text())
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        for system, n, d, taus in _published_records(payload):
            if n < minimum_particle_number:
                continue
            started = time.time()
            counts = {
                "span_is_hyperplane": 0,
                "trace_identity": 0,
                "decoded": 0,
                "tau_reconstructed_oriented": 0,
            }
            max_profile_variables = 0
            max_states = 0
            for tau in taus:
                positive = threshold_family(tau, n, +1)
                negative = threshold_family(tau, n, -1)
                if len(positive) == inversions(tau):
                    counts["trace_identity"] += 1
                zero = [
                    subset
                    for subset in itertools.combinations(range(d), n)
                    if sum(tau[index] for index in subset) == 0
                ]
                _recovered, rank = annihilator(zero, d)
                if rank == d - 1:
                    counts["span_is_hyperplane"] += 1
                decoded = decode_signed_pair(
                    shadow_of(positive, d), shadow_of(negative, d), n
                )
                if (
                    decoded is not None
                    and decoded.positive.family == positive
                    and decoded.negative.family == negative
                ):
                    counts["decoded"] += 1
                if decoded is not None and decoded.tau == tuple(tau):
                    counts["tau_reconstructed_oriented"] += 1
                if decoded is not None:
                    for side in (decoded.positive, decoded.negative):
                        max_profile_variables = max(
                            max_profile_variables, side.profile_variables
                        )
                        max_states = max(max_states, side.search_states)
            results[system] = {
                "system": system,
                "particle_number": n,
                "rank": d,
                "source_artifact": f"results/data/{name}.json",
                "source_sha256": digest,
                "representatives": len(taus),
                **counts,
                "max_profile_variables": max_profile_variables,
                "max_search_states_per_side": max_states,
                "seconds": round(time.time() - started, 1),
            }
    return results


def particle_hole_gate(systems, trials: int, seed: int) -> dict:
    """Measure what the particle-hole reduction actually buys.

    The claim under test is not correctness, which is a bijection, but cost.
    Both routes are run on the same shadows and compared field by field.
    """
    rng = random.Random(seed)
    records: dict[str, dict] = {}
    for n, d in systems:
        started = time.time()
        checked = 0
        answer_mismatches = 0
        variable_count_differs = 0
        literal_states = 0
        reduced_states = 0
        reduced_better = 0
        reduced_worse = 0
        for _ in range(trials):
            tau = tuple(rng.randrange(-5, 6) for _ in range(d))
            family = threshold_family(tau, n)
            if not family:
                continue
            shadow = shadow_of(family, d)
            literal = decode_shadow(shadow, n, particle_hole=False)
            reduced = decode_shadow(shadow, n, particle_hole=True)
            if literal is None or reduced is None:
                answer_mismatches += 1
                continue
            checked += 1
            if literal.family != family or reduced.family != family:
                answer_mismatches += 1
            if literal.profile_variables != reduced.profile_variables:
                variable_count_differs += 1
            literal_states += literal.search_states
            reduced_states += reduced.search_states
            if reduced.search_states < literal.search_states:
                reduced_better += 1
            elif reduced.search_states > literal.search_states:
                reduced_worse += 1
        records[f"({n},{d})"] = {
            "system": f"({n},{d})",
            "particle_number": n,
            "rank": d,
            "reduced_layer": min(n, d - n),
            "shadows_checked": checked,
            "answer_mismatches": answer_mismatches,
            "rows_where_variable_count_differs": variable_count_differs,
            "total_states_literal_layer": literal_states,
            "total_states_reduced_layer": reduced_states,
            "rows_reduced_cheaper": reduced_better,
            "rows_reduced_dearer": reduced_worse,
            "seconds": round(time.time() - started, 1),
        }
    return records


def exact_trace_survivor(n: int, d: int, rng) -> tuple | None:
    """Build one exact Ressayre trace survivor from a degenerate level multiset.

    Uniform random normals are the WRONG stress population: most of them fail
    the trace budget and no generator would ever see them. This instead samples
    a degenerate level pattern of the kind real candidates have, computes the
    orbit-invariant positive-weight count, and asks the repository's own
    survivor generator for an arrangement with exactly that inversion number.
    """
    level_count = rng.randint(2, min(6, d))
    levels = sorted(rng.sample(range(-6, 7), level_count))
    if not (levels[0] < 0 < levels[-1]):
        return None
    base = tuple(sorted((rng.choice(levels) for _ in range(d)), reverse=True))
    if len(set(base)) < 2 or not any(base):
        return None

    positive = sum(
        1
        for subset in itertools.combinations(range(d), n)
        if sum(base[index] for index in subset) > 0
    )
    try:
        candidate = admissible_candidate_from_tau(n, list(base))
    except (ValueError, TypeError):
        return None
    weights = exterior_weights(n, d)
    _on, below = _level_sets(weights, candidate.h, candidate.z)
    if len(below) != positive:                            # pragma: no cover
        raise AssertionError("the centered bridge disagreed on the below count")
    if oc.trace_survivor_count(candidate.h, positive) == 0:
        return None

    value_map = dict(zip(candidate.h, base, strict=True))
    for arrangement in oc.generate_trace_survivors(candidate.h, positive):
        tau = content_free(tuple(value_map[value] for value in arrangement))
        if inversions(tau) != positive:                   # pragma: no cover
            raise AssertionError("generated arrangement missed its inversion number")
        zero = [
            subset
            for subset in itertools.combinations(range(d), n)
            if sum(tau[index] for index in subset) == 0
        ]
        _recovered, rank = annihilator(zero, d)
        if rank == d - 1:
            return tau
    return None


def trace_survivor_stress(n: int, d: int, wanted: int, seed: int) -> dict:
    """Decode exact trace survivors at higher `N`, which is the relevant test."""
    rng = random.Random(seed)
    started = time.time()
    decoded_rows = 0
    failures = 0
    attempts = 0
    max_profile_variables = 0
    max_states = 0
    while decoded_rows + failures < wanted and attempts < 4000:
        attempts += 1
        tau = exact_trace_survivor(n, d, rng)
        if tau is None:
            continue
        positive = threshold_family(tau, n, +1)
        negative = threshold_family(tau, n, -1)
        decoded = decode_signed_pair(
            shadow_of(positive, d), shadow_of(negative, d), n
        )
        if (
            decoded is None
            or decoded.positive.family != positive
            or decoded.negative.family != negative
            or decoded.tau != tuple(tau)
        ):
            failures += 1
            continue
        decoded_rows += 1
        for side in (decoded.positive, decoded.negative):
            max_profile_variables = max(
                max_profile_variables, side.profile_variables
            )
            max_states = max(max_states, side.search_states)
    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "slice_size": len(exterior_weights(n, d)),
        "requested": wanted,
        "attempts": attempts,
        "exact_trace_survivors_decoded": decoded_rows,
        "failures": failures,
        "max_profile_variables": max_profile_variables,
        "max_search_states_per_side": max_states,
        "seed": seed,
        "seconds": round(time.time() - started, 1),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-dual", action="store_true")
    parser.add_argument("--skip-wall", action="store_true")
    parser.add_argument("--skip-stress", action="store_true")
    parser.add_argument("--particle-hole-trials", type=int, default=40)
    parser.add_argument(
        "-o",
        "--out",
        type=pathlib.Path,
        default=DATA / "cocircuit_chow_higher_n.json",
    )
    arguments = parser.parse_args()

    dual: dict[str, dict] = {}
    if not arguments.skip_dual:
        for n, d in DUAL_SYSTEMS:
            print(f"({n},{d}) dual gate ...", flush=True)
            record = dual_gate(n, d)
            dual[record["system"]] = record
            print(
                f"  hyperplanes {record['admissible_hyperplanes']} "
                f"orbits {record['sd_unoriented_orbits']} "
                f"survivors {record['trace_survivors']} "
                f"pairwise {record['agrees_with_flat_lattice_pairwise']} "
                f"decode failures {record['decode_failures']} "
                f"tau failures {record['tau_reconstruction_failures']}",
                flush=True,
            )

    print("published higher-N representatives ...", flush=True)
    published = published_gate()
    for system, record in sorted(published.items()):
        print(
            f"  {system} {record['representatives']} representatives, "
            f"trace {record['trace_identity']} "
            f"decoded {record['decoded']} "
            f"tau {record['tau_reconstructed_oriented']}",
            flush=True,
        )

    print("particle-hole reduction ...", flush=True)
    particle_hole = particle_hole_gate(
        ((5, 8), (6, 9), (7, 10)), arguments.particle_hole_trials, STRESS_SEED
    )
    for system, record in sorted(particle_hole.items()):
        print(
            f"  {system} variables differ on "
            f"{record['rows_where_variable_count_differs']} of "
            f"{record['shadows_checked']}, states "
            f"{record['total_states_literal_layer']} -> "
            f"{record['total_states_reduced_layer']}",
            flush=True,
        )

    stress: dict[str, dict] = {}
    if not arguments.skip_stress:
        for n, d, wanted in STRESS_PLAN:
            print(f"({n},{d}) exact trace-survivor stress ...", flush=True)
            record = trace_survivor_stress(n, d, wanted, STRESS_SEED)
            stress[record["system"]] = record
            print(
                f"  decoded {record['exact_trace_survivors_decoded']} "
                f"failures {record['failures']} "
                f"max states {record['max_search_states_per_side']}",
                flush=True,
            )

    wall = {}
    if not arguments.skip_wall:
        n, d = WALL_SYSTEM
        print(f"({n},{d}) wall probe at {WALL_NODE_BUDGET} nodes ...", flush=True)
        wall = wall_probe(n, d, WALL_NODE_BUDGET)
        print(f"  completed: {wall['completed']} ({wall['seconds']}s)", flush=True)

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/cocircuit_chow_higher_n.py",
        "status": "arbitrary_particle_number_gates",
        "machine": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "theorems": {
            "generality": (
                "C1 needs only that every exterior weight has coordinate sum "
                "N != 0, and T1 to T5 are fixed-weight statements, so none of "
                "them uses N = 3. Proved."
            ),
            "particle_hole_shadow": (
                "An n-uniform family F of size m with shadow c has complement "
                "family F^c in the (d-n) layer with shadow m - c. "
                "Complementation is a bijection between the layers, so "
                "decoding either side is exact. Proved."
            ),
        },
        "measured_not_proved": {
            "particle_hole_is_not_a_size_reduction": (
                "k -> m - k is a bijection between the two layers' occupancy "
                "profiles, so the profile VARIABLE COUNT is identical in both "
                "layers. The reduction changes only the coefficients and the "
                "target, and its effect on branch count is measured here. It "
                "is not a size reduction."
            ),
            "half_filling_wall": (
                "The unsymmetric reverse search does not finish at (4,8) "
                "inside the recorded node budget. This is an operational "
                "measurement of one implementation, not a theorem about the "
                "problem."
            ),
            "published_gate_scope": (
                "(4,8) is the complete 15-row pilot system. The (4,9), (4,10) "
                "and (5,10) rows are one representative per Levi signature "
                "class, not every facet row, so those three are "
                "representative-level gates."
            ),
        },
        "dual_populations": dual,
        "published_higher_n": published,
        "particle_hole": particle_hole,
        "trace_survivor_stress": stress,
        "half_filling_wall": wall,
    }

    def strip(record):
        if isinstance(record, dict):
            return {
                key: strip(value)
                for key, value in record.items()
                if key not in TIMING_FIELDS
            }
        return record

    payload["canonical_sha256_without_this_field_or_timings"] = hashlib.sha256(
        json.dumps(strip(payload), sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
