#!/usr/bin/env python3
"""Repair three stale numerical support fields in the exact state atlas.

The exact ``closed_form`` objects for these transported design records were
correct, but the convenience ``support`` arrays still described the numerical
discovery states that preceded exactification.  This one-time migration checks
the exact objects before replacing only the redundant numerical arrays.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


TARGETS = {
    ("(3,9)", 53): [
        [0, 1, 6],
        [0, 2, 7],
        [0, 3, 4],
        [1, 2, 4],
        [1, 3, 7],
        [2, 3, 5],
    ],
    ("(3,10)", 106): [
        [0, 1, 6],
        [0, 2, 7],
        [0, 3, 4],
        [1, 2, 4],
        [1, 3, 7],
        [2, 3, 5],
    ],
    ("(4,10)", 53): [
        [0, 1, 2, 7],
        [0, 1, 3, 8],
        [0, 1, 4, 5],
        [0, 2, 3, 5],
        [0, 2, 4, 8],
        [0, 3, 4, 6],
    ],
}
WEIGHTS = [2, 1, 3, 3, 1, 2]
DENOMINATOR = 12
PRETTY = ["sqrt(6)/6", "sqrt(3)/6", "1/2", "1/2", "sqrt(3)/6", "sqrt(6)/6"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("states", type=Path)
    args = parser.parse_args()

    rows = [json.loads(line) for line in args.states.read_text(encoding="utf-8").splitlines()]
    repaired: set[tuple[str, int]] = set()
    for record in rows:
        key = (record["system"], int(record["index"]))
        if key not in TARGETS:
            continue
        expected = {
            "den": DENOMINATOR,
            "weights": WEIGHTS,
            "pretty": PRETTY,
            "support_dets": TARGETS[key],
        }
        if record.get("closed_form") != expected:
            raise AssertionError(f"exact closed form drifted at {key}")
        record["support"] = [
            [determinant, math.sqrt(weight / DENOMINATOR), 0.0]
            for determinant, weight in zip(TARGETS[key], WEIGHTS)
        ]
        repaired.add(key)

    if repaired != set(TARGETS):
        raise AssertionError(f"missing targets: {sorted(set(TARGETS) - repaired)}")
    args.states.write_text(
        "".join(json.dumps(record, separators=(", ", ": ")) + "\n" for record in rows),
        encoding="utf-8",
    )
    print(f"repaired {len(repaired)} stale support arrays")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
