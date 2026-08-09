#!/usr/bin/env python3
"""Direct cocircuit discovery and constructive signed-Chow decoding.

TWO THINGS ARE MEASURED HERE, and they attack different halves of the fixed-N
generator.

DISCOVERY. `gpc_census.generation.cocircuit` enumerates admissible hyperplanes
as cocircuit zero sets, going straight to codimension one instead of building
the closed-flat lattice underneath it. That is legitimate because every
exterior weight has coordinate sum `n`, so admissibility (affine rank `d-2`) is
exactly linear rank `d-1`, which is exactly the zero set of a cocircuit. This
script gates that backend against the repository's sealed populations PAIR FOR
PAIR, not merely by count, and then reports what it costs against the two
enumerators the repository already has. The honest headline is in the cost
comparison, not in the agreement: the direct route is unsymmetric, so it is a
cross-check at low rank and not a replacement for `orbit_canonical`.

REPRESENTATION. `gpc_census.generation.chow_decoder` inverts a signed Chow
shadow. `docs/signed_chow_projection.md` proved that the pair `(c_+, c_-)`
determines the normal and recorded, as its one open operational gap, that no
decoder existed. T4 and T5 supply one: equal shadow coordinates force genuine
permutation symmetry of the positive family, and averaging a separator over
that symmetry gives a block-constant separator whose levels are `c_+`-ordered,
so the unknown is an upper ideal of occupancy profiles rather than a subset of
the slice. The decoder is run on every nonstructural trace survivor the
cocircuit backend finds, and separately on random threshold shadows at higher
rank, which is a stress test of decoding and NOT a count of GPC candidates.

Run:

    uv run scripts/cocircuit_chow_decoder.py
    uv run scripts/cocircuit_chow_decoder.py --ranks 5 6 7 --skip-cost
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
    decode_shadow,
    decode_signed_pair,
    shadow_of,
    threshold_family,
)
from gpc_census.generation.cocircuit import (  # noqa: E402
    enumerate_cocircuit_hyperplanes,
    trace_survivors,
)
from gpc_census.generation.exterior import exterior_weights  # noqa: E402
from gpc_census.generation.ressayre import (  # noqa: E402
    _level_sets,
    _negative_root_set,
    admissible_candidate_from_tau,
)

PARTICLE_NUMBER = 3
TIMING_FIELDS = ("seconds", "cocircuit_seconds", "flat_lattice_seconds",
                 "orbit_augmentation_seconds")
STRESS_PLAN = ((8, 500), (10, 300), (12, 200), (20, 50))
STRESS_SEED = 20260809


def orbit_key(tau) -> tuple[int, ...]:
    """Unoriented `S_d` orbit key: the coordinate multiset, up to global sign.

    A permutation acts on a primitive normal and on nothing else, and the flat
    is unoriented, so this needs no canonicalizer and is independent of the
    repository's own. `docs/orbit_canonical_flats.md` records the trap this
    avoids: forgetting the sign normalization answers the ORIENTED question and
    reports 15 orbits at rank 6 instead of 8.
    """
    forward = tuple(sorted(tau))
    backward = tuple(sorted(-value for value in tau))
    return min(forward, backward)


def signed_shadow(tau, n: int):
    """Return `(c_+, c_-)` for one normal, by one pass over the slice."""
    d = len(tau)
    positive = [0] * d
    negative = [0] * d
    for subset in itertools.combinations(range(d), n):
        total = sum(tau[index] for index in subset)
        if total > 0:
            for index in subset:
                positive[index] += 1
        elif total < 0:
            for index in subset:
                negative[index] += 1
    return tuple(positive), tuple(negative)


def repo_convention_survivors(n: int, d: int, normals):
    """Trace survivors decided by the repository's centered-convention test.

    The cocircuit backend applies the trace test in homogeneous coordinates:
    `#{omega : tau . omega > 0}` against `#{i < j : tau_j > tau_i}`. The
    repository applies it as `below_weight_count == negative_root_count` after
    the `h_i = S - d tau_i` bridge. They are the same test, and this recomputes
    the repository side so the claim is checked rather than asserted.
    """
    weights = exterior_weights(n, d)
    survivors = set()
    for tau in normals:
        for sign in (1, -1):
            oriented = tuple(sign * value for value in tau)
            candidate = admissible_candidate_from_tau(n, list(oriented))
            _on, below = _level_sets(weights, candidate.h, candidate.z)
            if len(below) == len(_negative_root_set(candidate.h)):
                survivors.add(oriented)
    return survivors


_ENUMERATIONS: dict[tuple[int, int], tuple[object, float]] = {}


def cached_cocircuit(n: int, d: int):
    """Enumerate once per rank, and remember what the first run cost.

    The gates and the cost comparison both need the same enumeration, and the
    flat-lattice enumerator it is timed against is itself memoised, so a second
    timing would compare a cold run with a warm one and report a speedup that
    is only a cache.
    """
    key = (n, d)
    if key not in _ENUMERATIONS:
        started = time.time()
        enumeration = enumerate_cocircuit_hyperplanes(n, d)
        _ENUMERATIONS[key] = (enumeration, round(time.time() - started, 1))
    return _ENUMERATIONS[key]


def cocircuit_gate(n: int, d: int, *, compare_convention: bool) -> dict:
    """Enumerate by cocircuit, then gate against the flat lattice."""
    enumeration, seconds = cached_cocircuit(n, d)

    survivors = list(trace_survivors(enumeration))
    nonstructural = [tau for tau, positive in survivors if positive > 0]

    reference = enumerate_admissible_hyperplanes_by_flats(n, d)
    reference_pairs = {(row.h, row.z) for row in reference.hyperplanes}
    cocircuit_pairs = set(enumeration.affine_pairs())

    record = {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "exterior_weights": enumeration.weight_count,
        "admissible_hyperplanes": enumeration.hyperplane_count,
        "reverse_search_nodes": enumeration.node_count,
        "rank_preserving_modulus": enumeration.modulus,
        "trace_survivors": len(survivors),
        "nonstructural_trace_survivors": len(nonstructural),
        "sd_unoriented_orbits": len({orbit_key(tau) for tau in enumeration.normals}),
        "flat_lattice_hyperplanes": len(reference.hyperplanes),
        "agrees_with_flat_lattice_pairwise": cocircuit_pairs == reference_pairs,
        "distinct_affine_pairs": len(cocircuit_pairs),
        "seconds": seconds,
    }
    if compare_convention:
        repo_side = repo_convention_survivors(n, d, enumeration.normals)
        record["agrees_with_centered_convention_test"] = repo_side == {
            tau for tau, _positive in survivors
        }
    else:
        record["agrees_with_centered_convention_test"] = None
    return record, nonstructural


def decoder_gate(n: int, d: int, normals) -> dict:
    """Decode both shadow sides of every nonstructural survivor."""
    started = time.time()
    decode_failures = 0
    tau_failures = 0
    max_profile_variables = 0
    max_states = 0
    total_states = 0
    sides = 0
    block_histogram: dict[int, int] = {}
    for tau in normals:
        positive_shadow, negative_shadow = signed_shadow(tau, n)
        decoded = decode_signed_pair(positive_shadow, negative_shadow, n)
        if decoded is None:
            decode_failures += 1
            continue
        if decoded.positive.family != threshold_family(tau, n, +1):
            decode_failures += 1
        if decoded.negative.family != threshold_family(tau, n, -1):
            decode_failures += 1
        # Oriented equality, not equality up to sign. The decoded positive
        # family fixes the orientation, so recovering `-tau` would be a
        # failure, and comparing sign-normalized representatives would hide it.
        if decoded.tau is None or decoded.tau != tuple(tau):
            tau_failures += 1
        for side in (decoded.positive, decoded.negative):
            max_profile_variables = max(max_profile_variables, side.profile_variables)
            max_states = max(max_states, side.search_states)
            total_states += side.search_states
            sides += 1
            blocks = len(side.blocks)
            block_histogram[blocks] = block_histogram.get(blocks, 0) + 1
    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "rows": len(normals),
        "decode_failures": decode_failures,
        "tau_reconstruction_failures": tau_failures,
        "max_profile_variables": max_profile_variables,
        "max_search_states_per_side": max_states,
        "mean_search_states_per_side": (
            round(total_states / sides, 4) if sides else 0
        ),
        "degree_block_histogram_both_sides": dict(sorted(block_histogram.items())),
        "seconds": round(time.time() - started, 1),
    }


def stress_gate(n: int, d: int, trials: int, seed: int) -> dict:
    """Decode random threshold shadows, to see how the search grows.

    THIS IS NOT A CANDIDATE COUNT. The normals are drawn uniformly and are not
    admissible, so nothing here says anything about `(3,20)` geometry. It
    measures only how the profile search behaves as the block structure widens.
    """
    rng = random.Random(seed)
    started = time.time()
    failures = 0
    max_states = 0
    max_profile_variables = 0
    total_states = 0
    scored = 0
    for _ in range(trials):
        tau = tuple(rng.randint(-8, 8) for _ in range(d))
        if not any(tau):
            continue
        family = threshold_family(tau, n, +1)
        decoded = decode_shadow(shadow_of(family, d), n)
        if decoded is None or decoded.family != family:
            failures += 1
            continue
        scored += 1
        max_states = max(max_states, decoded.search_states)
        max_profile_variables = max(
            max_profile_variables, decoded.profile_variables
        )
        total_states += decoded.search_states
    return {
        "system": f"({n},{d})",
        "rank": d,
        "trials": trials,
        "seed": seed,
        "failures": failures,
        "max_search_states": max_states,
        "mean_search_states": round(total_states / scored, 4) if scored else 0,
        "max_profile_variables": max_profile_variables,
        "seconds": round(time.time() - started, 1),
    }


def cost_comparison(n: int, d: int) -> dict:
    """Time the three enumerators that reach the same rank-`d-1` objects.

    Reported per machine and never averaged, per the sync invariant. The
    cocircuit and flat-lattice routes emit every hyperplane; the augmentation
    route emits orbit representatives and expands them, so the comparison is
    between what each one has to hold, not only between clocks.
    """
    enumeration, cocircuit_seconds = cached_cocircuit(n, d)

    started = time.time()
    flats = enumerate_admissible_hyperplanes_by_flats(n, d)
    flat_seconds = round(time.time() - started, 1)

    started = time.time()
    orbits, levels = oc.enumerate_orbits_by_augmentation(n, d)
    expanded = [
        sum(oc.orbit_size(mask, n, d) for mask in level) for level in levels
    ]
    orbit_seconds = round(time.time() - started, 1)

    return {
        "system": f"({n},{d})",
        "cocircuit_hyperplanes": enumeration.hyperplane_count,
        "cocircuit_reverse_search_nodes": enumeration.node_count,
        "cocircuit_seconds": cocircuit_seconds,
        "flat_lattice_flats_by_rank": list(flats.flat_counts_by_rank),
        "flat_lattice_stored_objects": sum(flats.flat_counts_by_rank),
        "flat_lattice_seconds": flat_seconds,
        "orbit_augmentation_orbits_by_rank": list(orbits),
        "orbit_augmentation_stored_objects": sum(orbits),
        "orbit_augmentation_hyperplanes": expanded[-1],
        "orbit_augmentation_seconds": orbit_seconds,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", nargs="*", type=int, default=[5, 6, 7, 8])
    parser.add_argument("--decoder-ranks", nargs="*", type=int, default=[6, 7, 8])
    parser.add_argument("--cost-rank", type=int, default=8)
    parser.add_argument("--skip-cost", action="store_true")
    parser.add_argument("--skip-stress", action="store_true")
    parser.add_argument(
        "-o",
        "--out",
        type=pathlib.Path,
        default=ROOT / "results/data/cocircuit_chow_decoder.json",
    )
    arguments = parser.parse_args()

    # BEFORE the gates, because the flat-lattice enumerator is memoised and a
    # warm call would time a dictionary lookup.
    cost = {}
    if not arguments.skip_cost:
        print(f"cost comparison at rank {arguments.cost_rank} ...", flush=True)
        cost = cost_comparison(PARTICLE_NUMBER, arguments.cost_rank)
        print(
            f"  cocircuit {cost['cocircuit_seconds']}s "
            f"flat lattice {cost['flat_lattice_seconds']}s "
            f"orbit augmentation {cost['orbit_augmentation_seconds']}s",
            flush=True,
        )

    cocircuit: dict[str, dict] = {}
    decoder: dict[str, dict] = {}
    for d in arguments.ranks:
        print(f"({PARTICLE_NUMBER},{d}) cocircuit ...", flush=True)
        record, nonstructural = cocircuit_gate(
            PARTICLE_NUMBER, d, compare_convention=True
        )
        cocircuit[str(d)] = record
        print(
            f"  hyperplanes {record['admissible_hyperplanes']} "
            f"nodes {record['reverse_search_nodes']} "
            f"orbits {record['sd_unoriented_orbits']} "
            f"survivors {record['trace_survivors']} "
            f"pairwise {record['agrees_with_flat_lattice_pairwise']} "
            f"({record['seconds']}s)",
            flush=True,
        )
        if d in arguments.decoder_ranks:
            print(f"({PARTICLE_NUMBER},{d}) decoder ...", flush=True)
            entry = decoder_gate(PARTICLE_NUMBER, d, nonstructural)
            decoder[str(d)] = entry
            print(
                f"  rows {entry['rows']} "
                f"decode failures {entry['decode_failures']} "
                f"tau failures {entry['tau_reconstruction_failures']} "
                f"max states {entry['max_search_states_per_side']} "
                f"({entry['seconds']}s)",
                flush=True,
            )

    stress: dict[str, dict] = {}
    if not arguments.skip_stress:
        for d, trials in STRESS_PLAN:
            print(f"stress d={d} ...", flush=True)
            stress[str(d)] = stress_gate(PARTICLE_NUMBER, d, trials, STRESS_SEED)
            print(
                f"  failures {stress[str(d)]['failures']} "
                f"max states {stress[str(d)]['max_search_states']}",
                flush=True,
            )

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/cocircuit_chow_decoder.py",
        "status": "independent_backend_and_constructive_decoder",
        "machine": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "theorems": {
            "C1_cocircuit_equivalence": (
                "Every exterior weight has coordinate sum n, so the affine hull "
                "of the slice misses the origin and affine rank d-2 is linear "
                "rank d-1. An admissible zero family is therefore a maximal "
                "coplanar subset of linear rank d-1, that is a cocircuit zero "
                "set, and conversely. Proved."
            ),
            "T4_block_symmetry": (
                "If c_+(tau)_i == c_+(tau)_j then transposing orbitals i and j "
                "fixes c_+, so by T1 it fixes B_+ itself. B_+ is a union of "
                "occupancy profiles over the equal-degree blocks. Proved."
            ),
            "T5_block_constant_separator": (
                "Averaging a separator over the product of symmetric groups on "
                "the equal-degree blocks preserves the strict-positive family "
                "and is block constant; T2 orders its levels by c_+; Abel "
                "summation then makes the selected profiles an upper ideal in "
                "prefix-sum dominance. Proved."
            ),
        },
        "measured_not_proved": {
            "decoder_complexity": (
                "No worst-case bound is claimed for the profile search. The "
                "stress block records how the maximum state count grows with "
                "rank and it grows; treat it as a measurement, never a bound."
            ),
            "symmetric_cocircuit_enumeration": (
                "The direct search here is unsymmetric, so it emits every "
                "hyperplane rather than every orbit. Cocircuit enumeration up "
                "to S_d is not implemented and rank 9 and above remain "
                "orbit_canonical's job."
            ),
            "downstream_gates_untouched": (
                "Cocircuit output is discovery data. Facetness, Hall, exact "
                "determinant, BDR birationality, redundancy and completeness "
                "are unchanged and still mandatory."
            ),
        },
        "cocircuit": cocircuit,
        "decoder": decoder,
        "decoder_stress": stress,
        "cost_comparison": cost,
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
