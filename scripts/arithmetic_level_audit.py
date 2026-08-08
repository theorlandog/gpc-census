#!/usr/bin/env python3
"""Arithmetic-level audit of the known GPC Levi mechanisms.

The audit takes the distinct ``tau`` levels of every mechanism recorded in
``results/data/levi_fusion_pilot.json`` and
``results/data/levi_fusion_heldout_results.json``, normalizes them by the gcd
of their successive differences, and records whether the normalized set is a
complete arithmetic progression, how many interior points are missing, and
whether the progression step divides the particle number.

The mechanism population is the set of distinct sorted-``tau`` shapes: the
pilot artifact stores one record per published row, so it is deduplicated by
sorted ``tau``, while the held-out artifact already stores one signature
record per shape.

Divisibility is a theorem rather than a measurement. If the distinct levels
lie in ``a + g*Z`` with ``gcd(a, g) = 1`` and some zero-grade allocation
exists, then ``N*a + g*S = 0`` for an integer ``S``, so ``g`` divides ``N*a``
and primitiveness gives ``g | N``. Because a published row is primitive, the
gcd of its distinct levels is 1, which forces ``gcd(a, g) = 1``. The audit
therefore reports ``step_divides_particle_number`` as a pipeline consistency
check, not as evidence for a pattern.

Run:

    uv run scripts/arithmetic_level_audit.py
    uv run scripts/arithmetic_level_audit.py --check
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
PILOT = DATA / "levi_fusion_pilot.json"
HELDOUT = DATA / "levi_fusion_heldout_results.json"
DEFAULT_OUT = DATA / "arithmetic_level_audit.json"

SCOPE = (
    "Distinct zero-grade Levi mechanisms from the rank-7/rank-8 pilot and the "
    "rank-9/rank-10 held-out artifact."
)


def load_mechanisms(pilot_path: Path, heldout_path: Path) -> list[dict]:
    """One entry per distinct sorted-tau shape, in artifact order."""

    pilot = json.loads(pilot_path.read_text())
    heldout = json.loads(heldout_path.read_text())
    mechanisms: list[dict] = []

    for system, payload in pilot["systems"].items():
        seen: dict[tuple[int, ...], dict] = {}
        for record in payload["records"]:
            seen.setdefault(tuple(sorted(record["tau"])), record)
        for record in seen.values():
            mechanisms.append(
                {
                    "system": system,
                    "particle_number": int(payload["particle_number"]),
                    "levels": [int(entry["value"]) for entry in record["levels"]],
                    "sectors": len(record["active_sectors"]),
                    "source": "pilot",
                }
            )

    for system, payload in heldout["systems"].items():
        for record in payload["signature_records"]:
            mechanisms.append(
                {
                    "system": system,
                    "particle_number": int(payload["particle_number"]),
                    "levels": [int(value) for value, _multiplicity in record["levels"]],
                    "sectors": len(record["allocations"]),
                    "source": "heldout",
                }
            )

    return mechanisms


def progression_step(levels: list[int]) -> int:
    step = 0
    for lower, upper in zip(levels, levels[1:]):
        step = math.gcd(step, upper - lower)
    return step or 1


def audit_mechanism(mechanism: dict) -> dict:
    levels = sorted(set(mechanism["levels"]))
    step = progression_step(levels)
    normalized = [(value - levels[0]) // step for value in levels]
    missing = sorted(set(range(normalized[-1] + 1)) - set(normalized))
    particle_number = mechanism["particle_number"]
    return {
        "is_complete_arithmetic_progression": not missing,
        "levels": levels,
        "missing_progression_points": missing,
        "normalized_levels": normalized,
        "particle_number": particle_number,
        "source": mechanism["source"],
        "step": step,
        "step_divides_particle_number": particle_number % step == 0,
        "system": mechanism["system"],
    }


def divisibility_witness(mechanism: dict) -> dict:
    """The hypotheses of the divisibility lemma, checked per mechanism."""

    levels = sorted(set(mechanism["levels"]))
    step = progression_step(levels)
    level_gcd = 0
    for value in levels:
        level_gcd = math.gcd(level_gcd, abs(value))
    return {
        "primitive_levels": level_gcd == 1,
        "offset_coprime_to_step": math.gcd(abs(levels[0]), step) == 1,
        "zero_grade_allocation_exists": mechanism["sectors"] >= 1,
        "step_divides_particle_number": mechanism["particle_number"] % step == 0,
    }


def summarize(records: list[dict]) -> dict:
    steps_by_particle_number: dict[str, dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    for record in records:
        bucket = steps_by_particle_number[str(record["particle_number"])]
        bucket[str(record["step"])] += 1
    return {
        "complete_arithmetic_progressions": sum(
            record["is_complete_arithmetic_progression"] for record in records
        ),
        "mechanisms": len(records),
        "more_than_one_point_missing": sum(
            len(record["missing_progression_points"]) > 1 for record in records
        ),
        "one_point_deleted_progressions": sum(
            len(record["missing_progression_points"]) == 1 for record in records
        ),
        "scope": SCOPE,
        "status": "post_hoc_exploratory_audit",
        "step_distribution_by_N": {
            key: dict(sorted(value.items(), key=lambda item: int(item[0])))
            for key, value in sorted(
                steps_by_particle_number.items(), key=lambda item: int(item[0])
            )
        },
        "step_divides_N": sum(
            record["step_divides_particle_number"] for record in records
        ),
    }


def build_payload(pilot_path: Path, heldout_path: Path) -> dict:
    mechanisms = load_mechanisms(pilot_path, heldout_path)
    records = [audit_mechanism(mechanism) for mechanism in mechanisms]

    for mechanism in mechanisms:
        witness = divisibility_witness(mechanism)
        if not all(witness.values()):
            raise AssertionError(
                f"divisibility lemma hypotheses failed for {mechanism}: {witness}"
            )

    return {
        "lemma": {
            "proof": (
                "For a zero-grade allocation, N*a + g*S = 0 for an integer S. "
                "Thus g divides N*a; gcd(a,g)=1 implies g divides N."
            ),
            "statement": (
                "If the distinct levels are a+gJ with gcd(a,g)=1 and at least "
                "one N-particle zero-grade allocation exists, then g divides N."
            ),
        },
        "records": records,
        "schema_version": 1,
        "summary": summarize(records),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, default=PILOT)
    parser.add_argument("--heldout", type=Path, default=HELDOUT)
    parser.add_argument("-o", "--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare with the stored artifact instead of writing it",
    )
    arguments = parser.parse_args()

    payload = build_payload(arguments.pilot, arguments.heldout)
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"

    if arguments.check:
        stored = json.loads(arguments.out.read_text())
        if stored != payload:
            raise AssertionError(f"{arguments.out} does not match a fresh regeneration")
        print(f"arithmetic level audit: {arguments.out.name} reproduces exactly")
        return 0

    arguments.out.write_text(serialized)
    summary = payload["summary"]
    print(f"wrote {arguments.out}")
    print(
        f"mechanisms={summary['mechanisms']} "
        f"exact_progressions={summary['complete_arithmetic_progressions']} "
        f"one_hole={summary['one_point_deleted_progressions']} "
        f"multi_hole={summary['more_than_one_point_missing']} "
        f"step_divides_N={summary['step_divides_N']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
