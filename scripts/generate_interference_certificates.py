#!/usr/bin/env python3
"""Generate exact global certificates for the paper's no-design claims.

The scope is the 11 negative cases in Proposition 1 and the v_B case in
Theorem 2. Discovery uses floating-point linear programming, but every vector
is rationalized and checked exactly before it enters the artifact. The branch
proof itself is purely combinatorial and exact.

Usage:
    uv run python scripts/generate_interference_certificates.py
    uv run python scripts/generate_interference_certificates.py --check
"""

import argparse
import json
import pathlib
import re
import sys
from fractions import Fraction

from gpc_census.certify import generate_global_interference_certificate


STATES = pathlib.Path("results/data/states.jsonl")
OUTPUT = pathlib.Path("results/data/interference_certificates.json")


def load_targets():
    targets = []
    for line in STATES.read_text().splitlines():
        record = json.loads(line)
        if record["classified"] != "INTERFERENCE":
            continue
        if record["system"] == "(3,8)" or (
            record["system"] == "(4,9)" and record["index"] == 65
        ):
            targets.append(record)
    targets.sort(key=lambda record: (record["system"], record["index"]))
    if len(targets) != 12:
        raise RuntimeError(f"expected 12 theorem/proposition targets, found {len(targets)}")
    return targets


def build_artifact():
    records = []
    for source in load_targets():
        n, d = map(int, re.findall(r"\d+", source["system"]))
        denominator = source["denominator"]
        integer_form = source["integer_form"]
        spectrum = [Fraction(value, denominator) for value in integer_form]
        record = generate_global_interference_certificate(
            n,
            d,
            spectrum,
            integer_form=integer_form,
            denominator=denominator,
            vertex_index=source["index"],
        )
        record["id"] = f"{n}-{d}-v{source['index']}"
        record["claim"] = "Proposition 1" if source["system"] == "(3,8)" else "Theorem 2"
        records.append(record)
        print(
            f"{record['id']}: {len(record['farkas_vectors'])} vectors, "
            f"{record['expanded_clause_count']} clauses, "
            f"{record['proof']['node_count']} proof nodes",
            flush=True,
        )

    return {
        "schema_version": 1,
        "record_count": len(records),
        "scope": {
            "claims": ["Proposition 1 negative cases", "Theorem 2"],
            "description": (
                "Exact disjunctive Farkas and branch certificates for every "
                "no-weighted-design claim carrying a numbered theorem or "
                "proposition in the paper."
            ),
        },
        "method": {
            "determinant_order": "lexicographic combinations",
            "symmetry": "all permutations inside equal-spectrum blocks",
            "arithmetic": "exact rational verification",
            "standalone_verifier": (
                "scripts/verify_interference_certificates_standalone.py"
            ),
        },
        "summary": {
            "farkas_vectors": sum(len(record["farkas_vectors"]) for record in records),
            "expanded_clauses": sum(record["expanded_clause_count"] for record in records),
            "proof_nodes": sum(record["proof"]["node_count"] for record in records),
            "maximum_proof_depth": max(
                record["proof"]["max_depth"] for record in records
            ),
        },
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=pathlib.Path, default=OUTPUT)
    args = parser.parse_args()

    artifact = build_artifact()
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            print(f"{args.output} is stale")
            return 1
        print(f"{args.output} is current")
        return 0

    args.output.write_text(rendered)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
