"""Does a finite isotropy character explain the quantization period?

`docs/symplectic_invariants.md` measured a period `p` in the multiplicity
sequence of each census vertex, with `p = q` at ten of the fourteen `(3,6)` and
`(3,7)` vertices and `p = 2q` at the other four, and left the doubling
unexplained except at `(3,6)#3`.

This script tests the mechanism pre-registered in
`docs/prereg_quantization_character.md`: for a `lambda`-preserving permutation
`sigma` fixing the certified state's line, with `wedge^n(P_sigma) psi = c psi`,

    omega(sigma) = c / prod_i det(sigma|block_i)^{lambda_i}

must satisfy `omega^k = 1` whenever the multiplicity at `k` is nonzero, so the
admissible `k` lie in the multiples of the lcm of the least admissible `k` over
all such `sigma`. That lcm is the predicted period, and it is compared against
the measured one.

Everything is exact: `c` is an exact algebraic number and the character
condition is evaluated symbolically, never numerically.

Run: uv run scripts/quantization_character.py
Writes results/data/quantization_characters.json.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data"
OUT = DATA / "quantization_characters.json"

from gpc_census.symplectic import (  # noqa: E402
    levi_permutation_stabilizer,
    predicted_quantization_period,
    spectrum_blocks,
)


def census_predictions(states, measured_keys) -> list[dict]:
    """Character-predicted period for every census vertex.

    The character side needs only the attainer and costs seconds, where the
    multiplicity side capped the previous campaign at rank 7. These rows are
    therefore PREDICTED, not measured: `verified` is true only for the fourteen
    vertices whose period was independently measured in
    `quantization_multiplicities.json`, and those are the only rows that carry
    any evidential weight for the law itself.
    """
    rows = []
    for record in states.values():
        q = int(record["denominator"] or 1)
        t0 = time.time()
        elements = levi_permutation_stabilizer(record)
        predicted = predicted_quantization_period(elements, q)
        key = (record["system"], record["index"])
        rows.append({
            "system": record["system"],
            "index": record["index"],
            "integer_form": list(record["integer_form"]),
            "denominator": q,
            "evidence": "exact",
            "spectrum_blocks": [len(b) for b in spectrum_blocks(record["integer_form"])],
            "stabilizer_elements": len(elements),
            "predicted_period": predicted,
            "predicted_period_exceeds_denominator": bool(predicted and predicted != q),
            "verified": key in measured_keys,
            "secs": round(time.time() - t0, 2),
        })
    return sorted(rows, key=lambda r: (r["system"], r["index"]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", action="store_true",
                    help="also predict the period for all 799 census vertices")
    args = ap.parse_args()
    measured = json.loads((DATA / "quantization_multiplicities.json").read_text())
    states = {(r["system"], r["index"]): r
              for r in (json.loads(x)
                        for x in (DATA / "states.jsonl").read_text().splitlines())}

    rows = []
    started = time.time()
    for entry in measured["rows"]:
        key = (entry["system"], entry["index"])
        record = states[key]
        q = int(entry["denominator"] or 1)
        t0 = time.time()
        elements = levi_permutation_stabilizer(record)
        predicted = predicted_quantization_period(elements, q)
        nontrivial = [e for e in elements if e["nontrivial"]]
        doubling = [e for e in elements if e["least_admissible_k"] not in (None, q)]
        row = {
            "system": entry["system"],
            "index": entry["index"],
            "integer_form": entry["integer_form"],
            "denominator": q,
            "evidence": "exact",
            "spectrum_blocks": [len(b) for b in spectrum_blocks(entry["integer_form"])],
            "stabilizer_elements": len(elements),
            "nontrivial_stabilizer_elements": len(nontrivial),
            "predicted_period": predicted,
            "measured_period": entry["quantization_period"],
            "matches": predicted == entry["quantization_period"],
            "divides": bool(predicted and entry["quantization_period"]
                            and entry["quantization_period"] % predicted == 0),
            "period_doubling_elements": [
                {"sigma": e["sigma"], "c": e["c_pretty"],
                 "block_determinants": e["block_determinants"],
                 "least_admissible_k": e["least_admissible_k"]}
                for e in doubling[:4]
            ],
            "elements": [
                {"sigma": e["sigma"], "c": e["c_pretty"],
                 "block_determinants": e["block_determinants"],
                 "least_admissible_k": e["least_admissible_k"]}
                for e in elements
            ],
            "secs": round(time.time() - t0, 2),
        }
        rows.append(row)
        print(f"  {entry['system']}#{entry['index']} q={q} "
              f"stab={len(elements)} predicted={predicted} "
              f"measured={entry['quantization_period']} "
              f"{'MATCH' if row['matches'] else 'MISMATCH'} ({row['secs']}s)",
              flush=True)

    payload = {
        "generated_by": "scripts/quantization_character.py",
        "prereg": "docs/prereg_quantization_character.md",
        "explains": "results/data/quantization_multiplicities.json",
        "evidence": "exact",
        "note": (
            "Finite isotropy characters of the certified attainers, tested as the "
            "mechanism behind the quantization period. Only lambda-preserving "
            "PERMUTATION elements are searched, which is a subgroup of the full "
            "stabilizer, so the predicted period always divides the measured one "
            "and equality is the falsifiable content."
        ),
        "rows": rows,
        "census_predictions": (
            census_predictions(states, {(r["system"], r["index"]) for r in rows})
            if args.census else []),
        "scorecard": scorecard(rows),
        "wall_secs": round(time.time() - started, 1),
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    for name, verdict in payload["scorecard"].items():
        if isinstance(verdict, dict) and "verdict" in verdict:
            print(f"  {name}: {verdict['verdict']}")
    return 0


def scorecard(rows: list[dict]) -> dict:
    scope = {"scope": len(rows), "independent_scope": len(rows),
             "independent_distinct_states": len(rows)}

    p1_bad = [f"{r['system']}#{r['index']}" for r in rows if not r["divides"]]
    p2_bad = [f"{r['system']}#{r['index']}" for r in rows if not r["matches"]]

    doubled = {f"{r['system']}#{r['index']}" for r in rows
               if r["measured_period"] != r["denominator"]}
    carriers = {f"{r['system']}#{r['index']}" for r in rows
                if r["period_doubling_elements"]}
    identity_only = [f"{r['system']}#{r['index']}" for r in rows
                     if r["nontrivial_stabilizer_elements"] == 0]

    return {
        "P1_predicted_divides_measured": {
            "verdict": "PASS" if not p1_bad else "FAIL",
            "note": "correctness gate, forced by the derivation; not evidence",
            "violations": p1_bad, **scope,
        },
        "P2_predicted_equals_measured": {
            "verdict": "PASS" if not p2_bad else "FAIL",
            "violations": p2_bad,
            "matched": len(rows) - len(p2_bad), **scope,
        },
        "P3_doubling_iff_permutation_obstruction": {
            "verdict": "PASS" if doubled == carriers else "FAIL",
            "period_doubled_vertices": sorted(doubled),
            "vertices_carrying_an_obstruction": sorted(carriers),
            "symmetric_difference": sorted(doubled ^ carriers), **scope,
        },
        "P4_identity_always_present": {
            "verdict": "PASS" if all(r["stabilizer_elements"] >= 1 for r in rows) else "FAIL",
            "vertices_with_identity_only": identity_only, **scope,
        },
    }


if __name__ == "__main__":
    raise SystemExit(main())
