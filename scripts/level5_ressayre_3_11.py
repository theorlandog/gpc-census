#!/usr/bin/env python3
"""Promote the four CLAIMED level-five (3,11) rows with exact Ressayre certificates.

The four rows that keep candidates 23, 26, 34 and 44 OPEN in the (3,11)
bracket are the "quadruple plus the last orbital" inequalities

    lambda_A + lambda_11 <= 2,    A in {2345}, {1346}, {1256}, {1247},

recorded in the census as AK-RMK-4.2.1 and KLY09-EXT-1/2/3. They come from
Klyachko 0904.2009 and CMP Remark 4.2.1 without proof, so the census tier is
CLAIMED and the four candidates that violate them stay OPEN.

This script runs the repository's own exact Ressayre screen on those four rows
and, separately, calibrates that screen at the two largest ranks where the
fixed-N=3 system is fully published.

Three independent things are computed and kept apart in the artifact.

1. CALIBRATION. At (3,9) and (3,10) every subset row of size 2..7 with right
   side 1 or 2 is screened, and every CERTIFIED row is then maximized by
   linear programming over the published system in
   ``src/gpc_census/data/constraints.json``. A left-side maximum at or below
   the right side proves the certified row is valid on the whole polytope, not
   merely at sampled points. This is the soundness test of the screen.

2. NEGATIVE CONTROL. Among the same screened rows, those that are
   demonstrably false (linear-programming maximum strictly above the right
   side) are split by screen status. Two of the rejection routes are cheap
   structural shortcuts that never form the tangent determinant: a row can be
   ``inadmissible``, meaning it is not a codimension-one exterior-weight
   hyperplane at all, or ``trace_rejected``, meaning the count of weights
   below it differs from the count of negative roots. Only the remaining
   statuses reach the determinant, and those are counted separately, because
   a control that every false row failed before the determinant would say
   nothing about the determinant. A false row that came back certified would
   refute the screen.

3. RANDOMIZED SWEEP. The calibration rows above all have zero-one
   coefficients, which is the shape of the rank-11 rows but a narrow slice of
   the row space. A seeded random sweep at rank 9 screens rows with
   coefficients in {-1,0,1,2} and right sides in {0,1,2,3}, and audits them
   the same way. It is the only part of the calibration that tests general
   coefficients.

4. THE RANK-11 RESULT. The four claimed rows are screened at (3,11), and every
   certified row is checked against all nineteen vertices the census has
   already certified TRUE at that rank.

What this does and does not establish. A Ressayre certificate is a validity
proof for the single inequality it certifies; it is not a claim of facetness,
irredundancy, or completeness of the (3,11) system, and the screen never makes
one. Validity is nevertheless enough to refute a candidate: a spectrum that
violates a valid inequality is outside the polytope.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import pulp

from gpc_census.generation.candidate_pipeline import screen_tau_candidates

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results/data/level5_ressayre_3_11.json"
BRACKET = ROOT / "docs/bracket_3_11.json"
SETTLEMENT = ROOT / "docs/bracket_3_11_settlement.json"
PUBLISHED = ROOT / "src/gpc_census/data/constraints.json"

N = 3
D11 = 11

# The four claimed level-five rows, as one-based orbital index sets with right
# side 2. The first four indices are the second-level quadruples of the (3,11)
# table; the fifth is the last orbital.
CLAIMED = (
    ("AK-RMK-4.2.1", (1, 2, 4, 7, 11)),
    ("KLY09-EXT-1", (2, 3, 4, 5, 11)),
    ("KLY09-EXT-2", (1, 3, 4, 6, 11)),
    ("KLY09-EXT-3", (1, 2, 5, 6, 11)),
)

CALIBRATION_RANKS = (9, 10)
CALIBRATION_SIZES = range(2, 8)
CALIBRATION_RHS = (1, 2)

# Statuses reached only after admissibility and the trace count both hold, so
# only these rows had the tangent determinant formed for them.
REACHED_DETERMINANT = frozenset(
    {"certified", "symbolic_zero_rejected", "evaluation_unresolved"}
)


def tau_of(subset, rhs: int, d: int) -> tuple[int, ...]:
    """Encode ``sum_{i in subset} lambda_i <= rhs`` in the pipeline's tau form.

    The convention, verified against the retained rank-7 rows, is
    ``tau_i = N * [i in A] - k``.
    """
    return tuple(N * (1 if i + 1 in subset else 0) - rhs for i in range(d))


def screen(rows, d: int, attempts: int = 32):
    """Return the exact screen status of each row, keyed by its tau."""
    taus = tuple(tau_of(subset, rhs, d) for subset, rhs in rows)
    result = screen_tau_candidates(
        taus,
        ressayre_attempts=attempts,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=1,
    ).as_dict()
    status = {tuple(row["tau"]): row.get("status") for row in result["candidates"]}
    return status, result["counts"], result["candidates"]


def published_system(d: int):
    entry = json.loads(PUBLISHED.read_text())[f"{N},{d}"]
    return [(list(map(F, q["coeffs"])), F(q["rhs"])) for q in entry["inequalities"]]


def lp_maximum(subset, d: int, system) -> float:
    """Maximize the row's left side over the published ordered polytope."""
    problem = pulp.LpProblem("row", pulp.LpMaximize)
    lam = [pulp.LpVariable(f"l{i}", lowBound=0, upBound=1) for i in range(d)]
    problem += pulp.lpSum(lam[i - 1] for i in subset)
    problem += pulp.lpSum(lam) == N
    for i in range(d - 1):
        problem += lam[i] >= lam[i + 1]
    for coeffs, rhs in system:
        problem += pulp.lpSum(float(c) * lam[i] for i, c in enumerate(coeffs)) <= float(rhs)
    problem.solve(pulp.PULP_CBC_CMD(msg=0))
    return pulp.value(problem.objective)


def calibrate(d: int, tolerance: float = 1e-7) -> dict:
    """Screen every small subset row at a published rank and audit the result."""
    system = published_system(d)
    table = {(tuple(coeffs), rhs) for coeffs, rhs in system}
    rows = [
        (subset, rhs)
        for size in CALIBRATION_SIZES
        for subset in itertools.combinations(range(1, d + 1), size)
        for rhs in CALIBRATION_RHS
    ]
    status, counts, _rows = screen(rows, d)
    by_tau = {}
    for subset, rhs in rows:
        by_tau.setdefault(tau_of(subset, rhs, d), (subset, rhs))

    certified, invalid, false_rows = [], [], []
    false_total = 0
    verbatim = 0
    for tau, outcome in status.items():
        subset, rhs = by_tau[tau]
        maximum = lp_maximum(subset, d, system)
        is_false = maximum > rhs + tolerance
        false_total += is_false
        if outcome == "certified":
            certified.append({
                "subset": list(subset),
                "rhs": rhs,
                "lp_maximum": maximum,
                "valid": not is_false,
            })
            row = tuple(1 if i + 1 in subset else 0 for i in range(d))
            verbatim += (row, F(rhs)) in table
            if is_false:
                invalid.append({"subset": list(subset), "rhs": rhs, "lp_maximum": maximum})
        if is_false:
            false_rows.append({
                "subset": list(subset),
                "rhs": rhs,
                "lp_maximum": maximum,
                "status": outcome,
            })

    reached = [r for r in false_rows if r["status"] in REACHED_DETERMINANT]
    return {
        "rank": d,
        "published_inequalities": len(system),
        "rows_screened": len(rows),
        "screen_counts": {k: v for k, v in counts.items() if v},
        "certified_rows": len(certified),
        "certified_verbatim_in_published_table": verbatim,
        "certified_rows_invalid_on_published_polytope": len(invalid),
        "invalid_rows": invalid,
        "negative_control": {
            "demonstrably_false_rows_screened": false_total,
            "false_row_statuses": dict(Counter(r["status"] for r in false_rows)),
            "false_rows_reaching_the_tangent_determinant": len(reached),
            "false_rows_certified": sum(
                1 for r in false_rows if r["status"] == "certified"
            ),
            "rows_reaching_the_determinant": reached,
        },
        "sound": not invalid,
        "certified": certified,
    }


SWEEP_RANK = 9
SWEEP_ROWS = 6000
SWEEP_SEED = 11
SWEEP_COEFFICIENTS = (-1, 0, 0, 1, 1, 2)
SWEEP_RHS = (0, 1, 2, 3)


def primitive(row):
    divisor = 0
    for value in row:
        divisor = math.gcd(divisor, abs(value))
    return tuple(value // divisor for value in row) if divisor else tuple(row)


def general_lp_maximum(coeffs, d: int, system) -> float:
    """Maximize an arbitrary integer row over the published ordered polytope."""
    problem = pulp.LpProblem("row", pulp.LpMaximize)
    lam = [pulp.LpVariable(f"l{i}", lowBound=0, upBound=1) for i in range(d)]
    problem += pulp.lpSum(float(c) * lam[i] for i, c in enumerate(coeffs))
    problem += pulp.lpSum(lam) == N
    for i in range(d - 1):
        problem += lam[i] >= lam[i + 1]
    for c, r in system:
        problem += pulp.lpSum(float(v) * lam[i] for i, v in enumerate(c)) <= float(r)
    problem.solve(pulp.PULP_CBC_CMD(msg=0))
    return pulp.value(problem.objective)


def randomized_sweep(tolerance: float = 1e-7) -> dict:
    """Screen seeded random general-coefficient rows and audit every verdict."""
    d = SWEEP_RANK
    system = published_system(d)
    rng = random.Random(SWEEP_SEED)
    by_tau = {}
    while len(by_tau) < SWEEP_ROWS:
        coeffs = [rng.choice(SWEEP_COEFFICIENTS) for _ in range(d)]
        rhs = rng.choice(SWEEP_RHS)
        if not any(coeffs):
            continue
        tau = primitive(tuple(N * c - rhs for c in coeffs))
        if not any(tau) or tau in by_tau:
            continue
        by_tau[tau] = (coeffs, rhs)

    result = screen_tau_candidates(
        tuple(by_tau),
        ressayre_attempts=32,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=1,
    ).as_dict()
    status = {primitive(tuple(r["tau"])): r.get("status") for r in result["candidates"]}

    certified, invalid, false_rows = [], [], []
    for tau, (coeffs, rhs) in by_tau.items():
        outcome = status.get(tau)
        if outcome is None:
            continue
        maximum = general_lp_maximum(coeffs, d, system)
        is_false = maximum > rhs + tolerance
        record = {
            "coefficients": coeffs,
            "rhs": rhs,
            "lp_maximum": maximum,
            "status": outcome,
        }
        if outcome == "certified":
            certified.append(record)
            if is_false:
                invalid.append(record)
        if is_false:
            false_rows.append(record)

    reached = [r for r in false_rows if r["status"] in REACHED_DETERMINANT]
    return {
        "rank": d,
        "seed": SWEEP_SEED,
        "coefficient_alphabet": sorted(set(SWEEP_COEFFICIENTS)),
        "right_sides": list(SWEEP_RHS),
        "rows_screened": len(by_tau),
        "screen_counts": {k: v for k, v in result["counts"].items() if v},
        "certified_rows": len(certified),
        "certified_rows_invalid_on_published_polytope": len(invalid),
        "invalid_rows": invalid,
        "negative_control": {
            "demonstrably_false_rows_screened": len(false_rows),
            "false_row_statuses": dict(Counter(r["status"] for r in false_rows)),
            "false_rows_reaching_the_tangent_determinant": len(reached),
            "false_rows_certified": sum(1 for r in false_rows if r["status"] == "certified"),
            "rows_reaching_the_determinant": reached,
        },
        "sound": not invalid,
        "certified": certified,
    }


def true_vertices() -> list[list[F]]:
    bracket = json.loads(BRACKET.read_text())["outer_vertices"]
    settlement = json.loads(SETTLEMENT.read_text())["candidates"]
    return [
        [F(x) for x in bracket[c["index"]]["spectrum"]]
        for c in settlement
        if c.get("verdict") == "TRUE-VERTEX"
    ]


def open_candidates() -> list[tuple[int, list[F]]]:
    """The candidates the pre-level-five settlement methods could not settle.

    Once this script's certificates land, those candidates carry a
    LEVEL5-RESSAYRE refutation in the settlement rather than an OPEN verdict,
    so keying off OPEN alone would find nothing and silently report an empty
    result. Accepting either keeps the two artifacts regenerable in any order.
    """
    bracket = json.loads(BRACKET.read_text())["outer_vertices"]
    settlement = json.loads(SETTLEMENT.read_text())["candidates"]
    return [
        (c["index"], [F(x) for x in bracket[c["index"]]["spectrum"]])
        for c in settlement
        if c.get("verdict") == "OPEN" or c.get("certificate") == "LEVEL5-RESSAYRE"
    ]


def rank_eleven() -> dict:
    rows = [(subset, 2) for _name, subset in CLAIMED]
    status, counts, screened = screen(rows, D11)
    attained = true_vertices()
    opens = open_candidates()

    # Export the full certificate payload so it can be replayed by a verifier
    # that shares no code with the generation package.
    by_tau_row = {tuple(row["tau"]): row for row in screened}

    records, certified_subsets = [], []
    for name, subset in CLAIMED:
        tau = tau_of(subset, 2, D11)
        outcome = status[tau]
        screened_row = by_tau_row[tau]
        record = {
            "family": name,
            "subset": list(subset),
            "rhs": 2,
            "tau": list(tau),
            "constraint": screened_row["constraint"],
            "census_tier_before": "CLAIMED",
            "screen_status": outcome,
        }
        assessment = screened_row.get("ressayre_assessment") or {}
        if outcome == "certified":
            certified_subsets.append((name, subset))
            record["census_tier_after"] = "PROVED"
            record["certificate"] = assessment["certificate"]
            record["below_weight_count"] = assessment["below_weight_count"]
            record["negative_root_count"] = assessment["negative_root_count"]
        records.append(record)

    violated_by_attained = [
        {"family": name, "subset": list(subset), "value": str(value)}
        for name, subset in certified_subsets
        for lam in attained
        if (value := sum(lam[j - 1] for j in subset)) > 2
    ]

    refuted = []
    for index, lam in opens:
        cutting = [
            {"family": name, "value": str(sum(lam[j - 1] for j in subset))}
            for name, subset in certified_subsets
            if sum(lam[j - 1] for j in subset) > 2
        ]
        refuted.append({
            "index": index,
            "spectrum": [str(x) for x in lam],
            "cut_by": cutting,
            "refuted": bool(cutting),
        })

    return {
        "rank": D11,
        "screen_counts": {k: v for k, v in counts.items() if v},
        "rows": records,
        "certified_rows": len(certified_subsets),
        "attained_vertices_checked": len(attained),
        "certified_rows_violated_by_an_attained_vertex": len(violated_by_attained),
        "violations": violated_by_attained,
        "open_candidates": refuted,
        "all_four_open_candidates_refuted": all(r["refuted"] for r in refuted),
    }


def build() -> dict:
    calibration = [calibrate(d) for d in CALIBRATION_RANKS]
    sweep = randomized_sweep()
    eleven = rank_eleven()
    audited = [*calibration, sweep]
    sound = all(c["sound"] for c in audited)
    controlled = all(
        c["negative_control"]["false_rows_certified"] == 0 for c in audited
    )
    settled = (
        sound
        and controlled
        and eleven["certified_rows"] == len(CLAIMED)
        and eleven["certified_rows_violated_by_an_attained_vertex"] == 0
        and eleven["all_four_open_candidates_refuted"]
    )
    return {
        "schema_version": 1,
        "system": "(3,11)",
        "generated_by": "scripts/level5_ressayre_3_11.py",
        "tau_convention": "tau_i = N * [i in A] - k for the row sum_{i in A} lambda_i <= k",
        "calibration": calibration,
        "randomized_sweep": sweep,
        "rank_eleven": eleven,
        "conclusion": {
            "screen_sound_at_published_ranks": sound,
            "negative_control_passed": controlled,
            "four_claimed_rows_certified": eleven["certified_rows"] == len(CLAIMED),
            "bracket_settles": settled,
            "tally_after": {"TRUE": 19, "REFUTED": 31, "OPEN": 0} if settled else None,
        },
        "nonclaims": [
            "A Ressayre certificate proves validity of the single inequality it "
            "certifies. It is not a claim of facetness or irredundancy, and the "
            "screen records both as not established.",
            "Nothing here claims the (3,11) system is complete. The bracket settles "
            "because valid inequalities refute the four remaining candidates, not "
            "because the system was enumerated.",
            "The randomized sweep is a sample, not an enumeration. It widens the "
            "coefficient shapes tested; it does not make the calibration complete "
            "over general rows.",
            "The calibration is complete at ranks 9 and 10 only. Soundness of the "
            "screen at rank 11 is inferred from the mathematics of the criterion "
            "and from the calibration, not proved by enumeration at rank 11.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    payload = build()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    for entry in (*payload["calibration"], payload["randomized_sweep"]):
        control = entry["negative_control"]
        print(
            f"rank {entry['rank']}: {entry['certified_rows']} certified, "
            f"{entry['certified_rows_invalid_on_published_polytope']} invalid; "
            f"{control['demonstrably_false_rows_screened']} false rows, "
            f"{control['false_rows_reaching_the_tangent_determinant']} reached the "
            f"determinant, {control['false_rows_certified']} certified"
        )
    eleven = payload["rank_eleven"]
    for row in eleven["rows"]:
        print(f"rank 11: {row['family']:>13s} {row['subset']} <= 2 -> {row['screen_status']}")
    for record in eleven["open_candidates"]:
        names = ", ".join(c["family"] for c in record["cut_by"]) or "nothing"
        print(f"  candidate {record['index']:>2}: cut by {names}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
