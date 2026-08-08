#!/usr/bin/env python3
"""Cross-artifact reconciliation of the mechanism-grammar frontier.

Four campaigns landed separately and their headline sentences look as though
they disagree:

* the Levi-fusion pilot and held-out audit refute a bounded Levi-block law;
* the arithmetic-level audit reports 78 of 83 mechanisms as exact arithmetic
  progressions;
* the bounded generator compression trial calls exact-progression-only a
  NO-GO;
* the full rank-7/rank-8 audit reports that every determinant-nonzero
  rank-7/rank-8 candidate has an exact arithmetic-progression level set.

This verifier reconciles them against each other rather than re-deriving each
one, and writes ``results/data/mechanism_grammar_synthesis.json``. It imports
neither gpc_census nor any third-party package.

Run:

    uv run scripts/verify_mechanism_grammar_synthesis.py
"""

from __future__ import annotations

import argparse
import json
import math
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
DEFAULT_OUT = DATA / "mechanism_grammar_synthesis.json"

ARITHMETIC = DATA / "arithmetic_level_audit.json"
TRIAL = DATA / "generator_compression_trial.json"
RANK78 = DATA / "full_rank78_compression_audit.json"
HELDOUT = DATA / "levi_fusion_heldout_results.json"
PILOT = DATA / "levi_fusion_pilot.json"

# The all-distinct (4,10) row that refutes the preregistered r <= 2N bound.
P1_WITNESS_TAU = (0, 1, 2, -3, 3, -6, -4, -5, -2, -1)
P1_WITNESS_N = 4


def zero_grade_descriptor(tau: tuple[int, ...], particle_number: int) -> dict:
    """Exact Levi blocks and zero-grade decomposition, recomputed from tau."""

    levels: dict[int, int] = {}
    for value in tau:
        levels[value] = levels.get(value, 0) + 1
    ordered = sorted(levels.items())
    heights = [height for height, _multiplicity in ordered]
    multiplicities = [multiplicity for _height, multiplicity in ordered]

    dimensions: list[int] = []
    for occupancies in product(*(range(m + 1) for m in multiplicities)):
        if sum(occupancies) != particle_number:
            continue
        if sum(h * n for h, n in zip(heights, occupancies)) != 0:
            continue
        dimension = 1
        for multiplicity, occupancy in zip(multiplicities, occupancies):
            dimension *= math.comb(multiplicity, occupancy)
        dimensions.append(dimension)

    return {
        "blocks": len(ordered),
        "fusion_rank": len(dimensions),
        "zero_grade_dimension": sum(dimensions),
        "ambient_dimension": math.comb(len(tau), particle_number),
    }


def rank_of(system: str) -> int:
    return int(system.strip("()").split(",")[1])


def check_population(arithmetic: dict, trial: dict) -> dict:
    """The four campaigns must be talking about the same 83 mechanisms."""

    records = arithmetic["records"]
    if len(records) != 83:
        raise AssertionError(f"expected 83 mechanisms, found {len(records)}")
    if trial["projection_width"]["mechanism_count"] != 83:
        raise AssertionError("trial mechanism count disagrees with the audit")
    if trial["known_mechanism_singleton_tomography"]["mechanisms"] != 83:
        raise AssertionError("tomography mechanism count disagrees with the audit")

    by_source: dict[str, int] = {}
    for record in records:
        by_source[record["source"]] = by_source.get(record["source"], 0) + 1
    if by_source != {"pilot": 13, "heldout": 70}:
        raise AssertionError(f"unexpected mechanism provenance split: {by_source}")

    return {
        "mechanisms": len(records),
        "pilot_mechanisms": by_source["pilot"],
        "heldout_mechanisms": by_source["heldout"],
        "structural_mechanisms": 1,
        "structural_note": (
            "The pilot contributes its (4,8) STRUCTURAL shape, so 82 of the 83 "
            "mechanisms are nonstructural. Every one-hole mechanism is "
            "nonstructural and held out."
        ),
        "agrees_with_trial_projection": True,
        "agrees_with_trial_tomography": True,
    }


def check_divisibility(arithmetic: dict) -> dict:
    """g | N is a theorem for primitive rows, so 83/83 is a consistency check."""

    records = arithmetic["records"]
    violations = [
        record for record in records if not record["step_divides_particle_number"]
    ]
    if violations:
        raise AssertionError(f"step does not divide N at {violations}")

    for record in records:
        levels = record["levels"]
        step = record["step"]
        level_gcd = 0
        for value in levels:
            level_gcd = math.gcd(level_gcd, abs(value))
        if level_gcd != 1:
            raise AssertionError(f"non-primitive level set: {levels}")
        if math.gcd(abs(levels[0]), step) != 1:
            raise AssertionError(f"offset not coprime to step: {levels}")
        normalized = record["normalized_levels"]
        rebuilt = [levels[0] + step * point for point in normalized]
        if rebuilt != levels:
            raise AssertionError(f"normalization does not rebuild levels: {levels}")

    return {
        "mechanisms_checked": len(records),
        "step_divides_particle_number": len(records),
        "primitive_level_sets": len(records),
        "offset_coprime_to_step": len(records),
        "status": "theorem_consequence_not_evidence",
        "statement": (
            "Levels in a+g*Z with gcd(a,g)=1 and a zero-grade allocation force "
            "g | N. Primitivity of a published row gives gcd(a,g)=1, so the "
            "83/83 observation is entailed rather than measured."
        ),
    }


def check_hole_budget(arithmetic: dict, trial: dict, rank78: dict) -> dict:
    """No population, in any of the three scopes, has a two-hole survivor."""

    summary = arithmetic["summary"]
    published = {
        "population": "published Levi mechanisms",
        "size": summary["mechanisms"],
        "exact_progressions": summary["complete_arithmetic_progressions"],
        "one_hole": summary["one_point_deleted_progressions"],
        "multi_hole": summary["more_than_one_point_missing"],
    }
    if published["multi_hole"] != 0:
        raise AssertionError("a published mechanism has two or more holes")

    bounded_search = trial["bounded_candidate_search"]
    bounded = {
        "population": "bounded N=3 Hall-feasible placements, ranks 7 to 9",
        "size": bounded_search["total_hall_feasible_placements"],
        "one_hole": bounded_search["one_hole_hall_feasible_placements"],
        "multi_hole": bounded_search["multi_hole_hall_feasible_placements"],
        "multi_hole_trace_placements_rejected_by_hall": bounded_search[
            "multi_hole_trace_placements"
        ],
    }
    if bounded["multi_hole"] != 0:
        raise AssertionError("a bounded Hall-feasible placement has two or more holes")
    if bounded["multi_hole_trace_placements_rejected_by_hall"] != 135:
        raise AssertionError("bounded multi-hole trace count changed")

    full_multi = 0
    full_one = 0
    full_total = 0
    for payload in rank78["systems"].values():
        for holes, count in payload["hole_distribution_all"].items():
            full_total += count
            if int(holes) == 1:
                full_one += count
            elif int(holes) > 1:
                full_multi += count
    full = {
        "population": "full rank-7/rank-8 trace survivors",
        "size": full_total,
        "one_hole": full_one,
        "multi_hole": full_multi,
    }
    if full_multi != 0:
        raise AssertionError("a full rank-7/rank-8 survivor has two or more holes")
    if full_total != rank78["summary"]["total_trace_survivors"]:
        raise AssertionError("hole distribution does not cover the survivor total")

    return {
        "populations": [published, bounded, full],
        "surviving_conjecture": (
            "Every Hall-feasible GPC candidate normal has arithmetic-progression "
            "levels with at most one missing point."
        ),
        "multi_hole_survivors_anywhere": 0,
    }


def check_hole_onset(arithmetic: dict, trial: dict, rank78: dict) -> dict:
    """The hole is unnecessary through rank 8 and necessary from rank 9."""

    exceptions = [
        record
        for record in arithmetic["records"]
        if record["missing_progression_points"]
    ]
    if len(exceptions) != 5:
        raise AssertionError(f"expected 5 one-hole mechanisms, found {len(exceptions)}")
    ranks = sorted({rank_of(record["system"]) for record in exceptions})
    if min(ranks) < 9:
        raise AssertionError(f"a published one-hole mechanism sits below rank 9: {ranks}")

    rank78_summary = rank78["summary"]
    if rank78_summary["one_hole_hall_feasible"] != 7:
        raise AssertionError("rank-7/rank-8 one-hole Hall-feasible count changed")
    if rank78_summary["one_hole_hall_feasible_exact_cancellation"] != 7:
        raise AssertionError("the seven rank-8 one-hole rows no longer all cancel")
    if not rank78_summary[
        "arithmetic_progression_suffices_after_exact_determinant_at_rank7_rank8"
    ]:
        raise AssertionError("rank-7/rank-8 exact-progression sufficiency flipped")

    certificates = trial["bounded_candidate_search"][
        "one_hole_modular_nonzero_certificates"
    ]
    if len(certificates) != 3:
        raise AssertionError("expected 3 one-hole nonvanishing certificates")
    for certificate in certificates:
        if len(certificate["tau"]) != 9:
            raise AssertionError("a one-hole certificate is not at rank 9")
        if certificate["determinant_residue"] == 0:
            raise AssertionError("a one-hole certificate has zero residue")

    return {
        "published_one_hole_mechanisms": len(exceptions),
        "published_one_hole_ranks": ranks,
        "published_one_hole_systems": sorted(
            record["system"] for record in exceptions
        ),
        "rank7_rank8_one_hole_hall_feasible": 7,
        "rank7_rank8_one_hole_determinant_nonzero": 0,
        "rank9_one_hole_modular_nonzero_certificates": len(certificates),
        "onset_rank": 9,
        "reconciliation": (
            "Exact arithmetic progressions are complete on the full rank-7 and "
            "rank-8 populations because all seven one-hole Hall-feasible rows "
            "have identically zero tangent determinant. From rank 9 the hole is "
            "load-bearing, established twice over: three modular-nonzero "
            "certificates in the bounded search, and five one-hole mechanisms "
            "among published facets at ranks 9 and 10. The trial NO-GO and the "
            "rank-7/rank-8 sufficiency are therefore consistent, separated by "
            "rank rather than in conflict."
        ),
    }


def check_bounded_full_agreement(trial: dict, rank78: dict) -> dict:
    """Two independently scoped rank-8 searches find the same one-hole shape."""

    full_shapes = {
        tuple(sorted(record["tau"]))
        for record in rank78["systems"]["(3,8)"]["one_hole_hall_feasible_exact_audit"]
    }
    bounded = [
        record
        for record in trial["bounded_candidate_search"]["systems"]["(3,8)"][
            "shape_records"
        ]
        if record["hole_count"] == 1 and record["hall_feasible_placements"] > 0
    ]
    bounded_shapes = {tuple(record["tau"]) for record in bounded}
    if full_shapes != bounded_shapes:
        raise AssertionError(
            f"rank-8 one-hole shapes disagree: {full_shapes} vs {bounded_shapes}"
        )
    if len(bounded) != 1:
        raise AssertionError("expected exactly one bounded one-hole shape at rank 8")
    placements = bounded[0]["hall_feasible_placements"]
    if placements != rank78["summary"]["one_hole_hall_feasible"]:
        raise AssertionError("one-hole placement counts disagree at rank 8")

    return {
        "shape": sorted(next(iter(full_shapes))),
        "normalized_levels": bounded[0]["normalized_levels"],
        "missing_levels": bounded[0]["missing_levels"],
        "hall_feasible_placements": placements,
        "note": (
            "The bounded window [-5,5] and the full repository trace-survivor "
            "population independently isolate the same single rank-8 one-hole "
            "shape with the same seven Hall-feasible placements."
        ),
    }


def check_width(trial: dict, rank78: dict) -> dict:
    projection = trial["projection_width"]
    rank78_summary = rank78["summary"]
    if projection["max_optimal_tree_width"] != 5:
        raise AssertionError("known-mechanism optimal width changed")
    if rank78_summary["max_optimal_fusion_width"] != 4:
        raise AssertionError("full rank-7/rank-8 optimal width changed")

    hall_zero_width_four = rank78["systems"]["(3,8)"][
        "optimal_fusion_width_distribution_hall_zero"
    ].get("4", 0)
    if hall_zero_width_four == 0:
        raise AssertionError("expected rejected rank-8 candidates at width four")

    return {
        "known_mechanism_max_optimal_width": projection["max_optimal_tree_width"],
        "known_mechanism_max_balanced_width": projection["max_balanced_tree_width"],
        "full_rank7_rank8_max_optimal_width": rank78_summary[
            "max_optimal_fusion_width"
        ],
        "rank8_hall_zero_candidates_at_width_four": hall_zero_width_four,
        "status": "evaluation_compression_only",
        "caveat": (
            "Width is an evaluation-cost result, not a facetness discriminator: "
            f"{hall_zero_width_four} rejected rank-8 Hall zeros also reach width "
            "four, so a small width does not select facets."
        ),
    }


def check_tomography(trial: dict, rank78: dict) -> dict:
    tomography = trial["known_mechanism_singleton_tomography"]
    if tomography["cross_shape_collisions"]:
        raise AssertionError("known-mechanism tomography has collisions")
    if tomography["distinct_system_specific_signatures"] != tomography["mechanisms"]:
        raise AssertionError("known-mechanism tomography is not injective")

    populations = []
    for system, payload in sorted(rank78["systems"].items()):
        survivors = payload["nonstructural_survivors"]
        signatures = payload["singleton_signature_count"]
        if survivors != signatures:
            raise AssertionError(f"{system} singleton incidence is not injective")
        if payload["singleton_collision_groups"] != 0:
            raise AssertionError(f"{system} reports singleton collision groups")
        populations.append(
            {
                "system": system,
                "nonstructural_survivors": survivors,
                "distinct_singleton_signatures": signatures,
            }
        )

    return {
        "known_mechanisms_separated": tomography["mechanisms"],
        "known_mechanism_collisions": 0,
        "full_population_injectivity": populations,
        "status": "no_theorem_yet",
        "caveat": (
            "Injectivity is measured on three finite populations. It licenses "
            "the incidence vector as a canonical hash and a fixed-weight "
            "Chow-parameter theorem attempt, not a proved separation."
        ),
    }


def check_block_count_refutation(heldout: dict) -> dict:
    """r <= 2N is refuted; f and dim W_0 are the quantities that stay small."""

    scorecard = heldout["scorecard_summary"]
    if scorecard["failed_predictions"] != ["P1_distinct_tau_levels_at_most_2N"]:
        raise AssertionError("the held-out scorecard no longer fails exactly P1")
    if len(scorecard["passed_predictions"]) != 5:
        raise AssertionError("the held-out scorecard no longer passes exactly P2-P6")

    witness = zero_grade_descriptor(P1_WITNESS_TAU, P1_WITNESS_N)
    if witness["blocks"] != 10:
        raise AssertionError("the (4,10) witness is not all-distinct")
    if witness["blocks"] <= 2 * P1_WITNESS_N:
        raise AssertionError("the (4,10) witness does not exceed 2N")
    if witness["fusion_rank"] != 9 or witness["zero_grade_dimension"] != 9:
        raise AssertionError(f"the (4,10) witness descriptor changed: {witness}")
    if witness["ambient_dimension"] != 210:
        raise AssertionError("the (4,10) ambient dimension changed")

    violating = sorted(
        system
        for system, payload in heldout["systems"].items()
        if payload["max_distinct_levels"] > 2 * payload["particle_number"]
    )
    if violating != ["(3,10)", "(4,10)"]:
        raise AssertionError(f"unexpected P1 violators: {violating}")

    return {
        "failed_prediction": "P1_distinct_tau_levels_at_most_2N",
        "violating_systems": violating,
        "witness_tau": list(P1_WITNESS_TAU),
        "witness_blocks": witness["blocks"],
        "witness_bound_2N": 2 * P1_WITNESS_N,
        "witness_fusion_rank": witness["fusion_rank"],
        "witness_zero_grade_dimension": witness["zero_grade_dimension"],
        "witness_ambient_dimension": witness["ambient_dimension"],
        "surviving_complexity_measure": (
            "Not the block count r, which reaches the orbital rank d. The small "
            "quantities are the fusion rank f and dim W_0, which stay at 9 out "
            "of an ambient 210 at the very row where r is largest."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--out", type=Path, default=DEFAULT_OUT)
    arguments = parser.parse_args()

    arithmetic = json.loads(ARITHMETIC.read_text())
    trial = json.loads(TRIAL.read_text())
    rank78 = json.loads(RANK78.read_text())
    heldout = json.loads(HELDOUT.read_text())

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/verify_mechanism_grammar_synthesis.py",
        "status": "cross_artifact_reconciliation",
        "inputs": [
            "results/data/arithmetic_level_audit.json",
            "results/data/full_rank78_compression_audit.json",
            "results/data/generator_compression_trial.json",
            "results/data/levi_fusion_heldout_results.json",
            "results/data/levi_fusion_pilot.json",
        ],
        "mechanism_population": check_population(arithmetic, trial),
        "divisibility": check_divisibility(arithmetic),
        "hole_budget": check_hole_budget(arithmetic, trial, rank78),
        "hole_onset": check_hole_onset(arithmetic, trial, rank78),
        "rank8_bounded_full_agreement": check_bounded_full_agreement(trial, rank78),
        "fusion_width": check_width(trial, rank78),
        "tomography": check_tomography(trial, rank78),
        "block_count_refutation": check_block_count_refutation(heldout),
        "evidence_boundary": (
            "Nothing here is a GPC proof. The reconciliation is exact on three "
            "finite populations: 83 published Levi mechanisms, the bounded N=3 "
            "normal window at ranks 7 to 9, and the full rank-7/rank-8 trace "
            "survivors. The one-hole grammar and the incidence separation both "
            "remain conjectures outside those windows, and fusion width is an "
            "evaluation cost rather than a facetness test."
        ),
    }

    arguments.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {arguments.out}")
    print("mechanism-grammar synthesis: PASS")
    print(
        "multi-hole survivors in any population: "
        f"{payload['hole_budget']['multi_hole_survivors_anywhere']}"
    )
    print(f"one-hole onset rank: {payload['hole_onset']['onset_rank']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
