#!/usr/bin/env python3
"""Tier-A design search for the four OPEN (3,11) bracket candidates.

Two independent things are computed and kept apart in the artifact.

1. An exact evaluation of every candidate against the thirteen stable
   Grassmann partial families, recording which rows each violates and at which
   tier. This is exact rational arithmetic.

2. A search for a one-hop-free design realising each candidate. One-hop-free
   means no two support triples share exactly two orbitals, which is precisely
   a partial Steiner triple system; for such a support the 1-RDM is diagonal
   with lambda_p = sum over blocks containing p, so a feasible support and
   weight vector is an explicit attaining state.

   The search runs as a mixed-integer program over all 165 triples, and CBC
   works in floating point, so an infeasible verdict there is EVIDENCE, not a
   certificate. Candidate 44 additionally carries an exact counting proof of
   infeasibility, which is recorded separately and is a certificate.

Nothing here settles attainability. A candidate with no one-hop-free design
can still be attainable by an interference state; what the search rules out is
the design route.
"""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction as F
from pathlib import Path

import pulp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BRACKET = ROOT / "docs/bracket_3_11.json"
DEFAULT_OUT = ROOT / "results/data/rank11_open_candidate_designs.json"

OPEN_CANDIDATES = (23, 26, 34, 44)
D = 11

# The thirteen stable Grassmann partial families at (3,11), one-based orbital
# indices. Tiers follow docs/partial_families_3_11_3_12.json: PROVED is CMP
# 2008 Thm 4.2.1/4.3.1, CLAIMED is CMP Remark 4.2.1 or Klyachko 0904.2009 and
# is explicitly not a certificate.
FAMILIES = (
    ("1st-head l1+l11<=1+2/(r-1)", (1, 11), F(1) + F(2, D - 1), "PROVED"),
    ("1st-pair l2+l10<=1", (2, 10), F(1), "PROVED"),
    ("1st-pair l3+l9<=1", (3, 9), F(1), "PROVED"),
    ("1st-pair l4+l8<=1", (4, 8), F(1), "PROVED"),
    ("1st-pair l5+l7<=1", (5, 7), F(1), "PROVED"),
    ("2nd-quad 2345<=2", (2, 3, 4, 5), F(2), "PROVED"),
    ("2nd-quad 1346<=2", (1, 3, 4, 6), F(2), "PROVED"),
    ("2nd-quad 1256<=2", (1, 2, 5, 6), F(2), "PROVED"),
    ("2nd-quad 1247<=2", (1, 2, 4, 7), F(2), "PROVED"),
    ("AK-RMK-4.2.1 [1,2,4,7]+[11]<=2", (1, 2, 4, 7, 11), F(2), "CLAIMED"),
    ("KLY09-EXT-1 [2,3,4,5]+[11]<=2", (2, 3, 4, 5, 11), F(2), "CLAIMED"),
    ("KLY09-EXT-2 [1,3,4,6]+[11]<=2", (1, 3, 4, 6, 11), F(2), "CLAIMED"),
    ("KLY09-EXT-3 [1,2,5,6]+[11]<=2", (1, 2, 5, 6, 11), F(2), "CLAIMED"),
)


def spectra(bracket_path: Path) -> dict[int, list[F]]:
    vertices = json.loads(bracket_path.read_text())["outer_vertices"]
    return {i: [F(x) for x in vertices[i]["spectrum"]] for i in OPEN_CANDIDATES}


def evaluate_families(lam: list[F]) -> dict:
    violated = []
    for name, indices, rhs, tier in FAMILIES:
        lhs = sum(lam[j - 1] for j in indices)
        if lhs > rhs:
            violated.append({
                "family": name,
                "tier": tier,
                "lhs": str(lhs),
                "rhs": str(rhs),
                "excess": str(lhs - rhs),
            })
    return {
        "violated": violated,
        "violated_claimed_count": sum(v["tier"] == "CLAIMED" for v in violated),
        "satisfies_every_proved_family": not any(v["tier"] == "PROVED" for v in violated),
    }


def blocks_and_conflicts():
    blocks = list(itertools.combinations(range(D), 3))
    conflicts = [
        (i, j)
        for i, j in itertools.combinations(range(len(blocks)), 2)
        if len(set(blocks[i]) & set(blocks[j])) == 2
    ]
    return blocks, conflicts


def design_search(lam: list[F], blocks, conflicts, time_limit: int = 600) -> dict:
    """Is there a partial Steiner support and weights realising lambda?"""
    problem = pulp.LpProblem("one_hop_free_design", pulp.LpMinimize)
    weight = [pulp.LpVariable(f"w{i}", lowBound=0) for i in range(len(blocks))]
    chosen = [pulp.LpVariable(f"y{i}", cat="Binary") for i in range(len(blocks))]
    problem += 0
    for i in range(len(blocks)):
        problem += weight[i] <= chosen[i]
    for i, j in conflicts:
        problem += chosen[i] + chosen[j] <= 1
    for orbital in range(D):
        problem += pulp.lpSum(
            weight[i] for i, block in enumerate(blocks) if orbital in block
        ) == float(lam[orbital])
    problem += pulp.lpSum(weight) == 1
    status = pulp.LpStatus[problem.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=time_limit))]
    found = None
    if status == "Optimal":
        found = [
            {"block": list(blocks[i]), "weight": weight[i].value()}
            for i in range(len(blocks))
            if (weight[i].value() or 0) > 1e-9
        ]
    return {
        "solver": "CBC via pulp, floating point",
        "status": status,
        "one_hop_free_design_found": found is not None,
        "support": found,
        "evidence_only": True,
    }


def exact_no_design_proof_44(lam: list[F]) -> dict:
    """Exact counting proof that candidate 44 admits no one-hop-free design.

    Heavy set H is the five orbitals at 1/2, light set L the six at 1/12.
    Writing a, b, c, d for the total weight on blocks meeting H in 3, 2, 1, 0
    orbitals, the three linear conditions force a = 1/2 + c + 2d >= 1/2. At
    most two blocks fit inside a five-set under the Steiner condition, and each
    of the three cases then fails on light-orbital coverage.
    """
    heavy = [p for p in range(D) if lam[p] == F(1, 2)]
    light = [p for p in range(D) if lam[p] == F(1, 12)]
    if len(heavy) != 5 or len(light) != 6:
        return {"applies": False}

    # the maximum partial Steiner family inside a five-set, by enumeration
    inside = list(itertools.combinations(range(5), 3))
    largest = max(
        (
            k
            for k in range(len(inside) + 1)
            for family in itertools.combinations(inside, k)
            if all(len(set(a) & set(b)) <= 1 for a, b in itertools.combinations(family, 2))
        ),
        default=0,
    )

    # with two inside-H blocks meeting in one orbital, count the free H-pairs
    first, second = (0, 1, 2), (2, 3, 4)
    used = {frozenset(p) for block in (first, second) for p in itertools.combinations(block, 2)}
    free_pairs = [p for p in itertools.combinations(range(5), 2) if frozenset(p) not in used]

    return {
        "applies": True,
        "heavy_orbitals": heavy,
        "light_orbitals": light,
        "forced_inside_weight_at_least": str(F(1, 2)),
        "derivation": "3a+2b+c = 5/2 and a+b+c+d = 1 and b+2c+3d = 1/2 give a - c - 2d = 1/2",
        "max_blocks_inside_a_five_set": largest,
        "case_two_blocks": {
            "free_heavy_pairs": len(free_pairs),
            "light_orbitals_needing_positive_weight": len(light),
            "contradiction": len(free_pairs) < len(light),
            "reason": "a = 1/2 forces c = d = 0, so every |S cap H| = 2 block carries "
                      "exactly one light orbital; there are too few free heavy pairs",
        },
        "case_one_block": {
            "contradiction": True,
            "reason": "the single inside-H block has weight 1/2, saturating its three "
                      "orbitals; the other two heavy orbitals admit only the one block on "
                      "the remaining heavy pair, leaving five light orbitals at zero",
        },
        "case_zero_blocks": {
            "contradiction": True,
            "reason": "a = 0 contradicts a >= 1/2",
        },
        "conclusion": "no one-hop-free design realises candidate 44",
        "exact": True,
    }


def build(bracket_path: Path = DEFAULT_BRACKET, time_limit: int = 600) -> dict:
    lam_by_index = spectra(bracket_path)
    blocks, conflicts = blocks_and_conflicts()
    records = {}
    for index in OPEN_CANDIDATES:
        lam = lam_by_index[index]
        distinct = sorted(set(lam), reverse=True)
        record = {
            "index": index,
            "spectrum": [str(x) for x in lam],
            "trace": str(sum(lam)),
            "distinct_values": len(distinct),
            "multiplicities": {str(v): lam.count(v) for v in distinct},
            "families": evaluate_families(lam),
            "design_search": design_search(lam, blocks, conflicts, time_limit),
        }
        if index == 44:
            record["exact_no_design_proof"] = exact_no_design_proof_44(lam)
        records[str(index)] = record
    return {
        "schema_version": 1,
        "system": "(3,11)",
        "generated_by": "scripts/rank11_open_candidate_designs.py",
        "triples": len(blocks),
        "one_hop_conflict_pairs": len(conflicts),
        "candidates": records,
        "scope": (
            "A one-hop-free support is a partial Steiner triple system, for which the "
            "1-RDM is diagonal. Absence of such a design closes the design route only: "
            "a candidate with no design may still be attainable by an interference "
            "state, and none of this bears on the claimed level-five inequalities."
        ),
        "nonclaims": [
            "No candidate is settled TRUE or REFUTED by this artifact.",
            "The mixed-integer infeasibility verdicts come from a floating-point solver "
            "and are evidence; only candidate 44 carries an exact proof.",
            "Nothing here proves or disproves AK-RMK-4.2.1 or the KLY09 extensions.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bracket", type=Path, default=DEFAULT_BRACKET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--time-limit", type=int, default=600)
    args = parser.parse_args()
    payload = build(args.bracket, args.time_limit)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for index, record in sorted(payload["candidates"].items(), key=lambda kv: int(kv[0])):
        families = record["families"]
        print(
            f"idx {index:>2}: {record['distinct_values']}-valued, "
            f"violates {families['violated_claimed_count']} claimed row(s), "
            f"design {record['design_search']['status']}"
        )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
