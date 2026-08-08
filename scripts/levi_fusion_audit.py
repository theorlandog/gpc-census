#!/usr/bin/env python3
"""Exact Levi block and zero-grade fusion audit for published GPC facets.

This script reads the exact row labels emitted by
``results/data/ordered_representation_stability.json``. For every primitive
homogeneous row ``tau . lambda <= 0``, it computes the eigenspace blocks of
``tau`` and the exact decomposition

    (wedge^N C^d)_0
      = direct_sum_{sum n_a=N, sum tau_a n_a=0}
          tensor_a wedge^{n_a} C^{m_a}.

The number of summands is called the fusion rank. The script also checks which
zero-grade sectors are occupied by shipped exact attaining states when their
1-RDM is diagonal in the shipped orbital frame.

The combinatorial row audit is exact and standard-library only. The state-frame
audit is explicitly conditional on the attainer gate and does not turn a
basis-dependent support statement into a basis-free fiber theorem.

Run:

    uv run scripts/levi_fusion_audit.py
    uv run scripts/levi_fusion_audit.py --systems 3,7 3,8 4,8
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
ORDERED = DATA / "ordered_representation_stability.json"
STATES = DATA / "states.jsonl"
ATTAINER_GATE = DATA / "attainer_gate_audit.json"
DEFAULT_OUT = DATA / "levi_fusion_audit.json"


def parse_system(text: str) -> tuple[int, int]:
    n, d = (int(part) for part in text.strip("()").split(","))
    return n, d


def canonical_json_hash(value: object) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def bounded_allocations(
    levels: tuple[int, ...],
    multiplicities: tuple[int, ...],
    particle_number: int,
) -> tuple[tuple[tuple[int, ...], int], ...]:
    """All zero-grade particle allocations and exact sector dimensions."""

    output: list[tuple[tuple[int, ...], int]] = []
    ranges = [range(multiplicity + 1) for multiplicity in multiplicities]
    for occupancy in product(*ranges):
        if sum(occupancy) != particle_number:
            continue
        if sum(level * count for level, count in zip(levels, occupancy)) != 0:
            continue
        dimension = math.prod(
            math.comb(multiplicity, count)
            for multiplicity, count in zip(multiplicities, occupancy)
        )
        output.append((tuple(occupancy), dimension))
    return tuple(output)


def module_factor(multiplicity: int, occupancy: int) -> str | None:
    if occupancy == 0:
        return None
    if occupancy == multiplicity:
        return f"det(C^{multiplicity})"
    if occupancy == 1:
        return f"C^{multiplicity}"
    return f"wedge^{occupancy}(C^{multiplicity})"


def module_name(
    multiplicities: tuple[int, ...],
    occupancy: tuple[int, ...],
) -> str:
    factors = [
        module_factor(multiplicity, count)
        for multiplicity, count in zip(multiplicities, occupancy)
    ]
    retained = [factor for factor in factors if factor is not None]
    return " tensor ".join(retained) if retained else "1"


def descriptor_from_tau(
    tau: tuple[int, ...],
    particle_number: int,
) -> dict[str, object]:
    counts = Counter(tau)
    levels = tuple(sorted(counts))
    multiplicities = tuple(counts[level] for level in levels)
    allocations = bounded_allocations(levels, multiplicities, particle_number)
    if not allocations:
        raise AssertionError("a facet row has an empty zero-grade exterior subspace")

    sectors = [
        {
            "occupancies": list(occupancy),
            "dimension": dimension,
            "module": module_name(multiplicities, occupancy),
        }
        for occupancy, dimension in allocations
    ]
    descriptor = {
        "levels": [
            {"value": level, "multiplicity": multiplicity}
            for level, multiplicity in zip(levels, multiplicities)
        ],
        "distinct_levels": len(levels),
        "fusion_rank": len(allocations),
        "zero_grade_dimension": sum(dimension for _occupancy, dimension in allocations),
        "one_dimensional_sectors": sum(
            dimension == 1 for _occupancy, dimension in allocations
        ),
        "active_sectors": sectors,
    }
    descriptor["oriented_signature_sha256"] = canonical_json_hash(
        {
            "particle_number": particle_number,
            "levels": descriptor["levels"],
            "active_sectors": sectors,
        }
    )

    negated_counts = Counter(-value for value in tau)
    negated_levels = tuple(sorted(negated_counts))
    negated_multiplicities = tuple(negated_counts[level] for level in negated_levels)
    negated_allocations = bounded_allocations(
        negated_levels,
        negated_multiplicities,
        particle_number,
    )
    positive_key = (
        tuple(zip(levels, multiplicities)),
        allocations,
    )
    negative_key = (
        tuple(zip(negated_levels, negated_multiplicities)),
        negated_allocations,
    )
    descriptor["unoriented_signature_sha256"] = canonical_json_hash(
        min(positive_key, negative_key)
    )
    return descriptor


def determinant_allocation(
    determinant: tuple[int, ...],
    tau: tuple[int, ...],
    levels: tuple[int, ...],
) -> tuple[int, ...]:
    index = {level: position for position, level in enumerate(levels)}
    counts = [0] * len(levels)
    for orbital in determinant:
        counts[index[tau[orbital]]] += 1
    return tuple(counts)


def load_state_data() -> tuple[dict[tuple[str, int], dict], dict[tuple[str, int], bool]]:
    if not STATES.exists() or not ATTAINER_GATE.exists():
        return {}, {}

    states = {
        (record["system"], int(record["index"])): record
        for record in (
            json.loads(line)
            for line in STATES.read_text().splitlines()
            if line.strip()
        )
    }
    gate_payload = json.loads(ATTAINER_GATE.read_text())
    gate = {
        (record["system"], int(record["index"])): bool(record["passes_gate"])
        for record in gate_payload["records"]
    }
    return states, gate


def tight_state_analysis(
    system: str,
    tau: tuple[int, ...],
    descriptor: dict[str, object],
    states: dict[tuple[str, int], dict],
    gate: dict[tuple[str, int], bool],
) -> dict[str, object]:
    candidates = [
        record
        for (record_system, _index), record in states.items()
        if record_system == system
    ]
    levels = tuple(entry["value"] for entry in descriptor["levels"])
    active = {
        tuple(entry["occupancies"])
        for entry in descriptor["active_sectors"]
    }

    records = []
    for state in candidates:
        integer_form = tuple(int(value) for value in state["integer_form"])
        if sum(left * right for left, right in zip(tau, integer_form)) != 0:
            continue
        key = (system, int(state["index"]))
        support = state.get("closed_form", {}).get("support_dets")
        if support is None:
            support = [entry[0] for entry in state.get("support", [])]
        determinants = [tuple(int(value) for value in det) for det in support]
        allocations = {
            determinant_allocation(determinant, tau, levels)
            for determinant in determinants
        }
        zero_grade = all(
            sum(tau[orbital] for orbital in determinant) == 0
            for determinant in determinants
        )
        records.append(
            {
                "index": int(state["index"]),
                "classified": state["classified"],
                "passes_attainer_gate": gate.get(key),
                "strict_occupation_nondegeneracy": all(
                    left != right
                    for left, right in zip(integer_form, integer_form[1:])
                ),
                "support_size": len(determinants),
                "support_is_zero_grade": zero_grade,
                "occupied_sector_count": len(allocations),
                "occupied_sectors": [list(value) for value in sorted(allocations)],
                "all_occupied_sectors_are_active": allocations <= active,
            }
        )

    gated = [record for record in records if record["passes_attainer_gate"] is True]
    strict = [
        record
        for record in gated
        if record["strict_occupation_nondegeneracy"]
    ]
    return {
        "tight_states": len(records),
        "tight_states_passing_attainer_gate": len(gated),
        "strict_nondegenerate_tight_states": len(strict),
        "gated_support_zero_grade": sum(
            record["support_is_zero_grade"] for record in gated
        ),
        "gated_occupied_sector_count_distribution": dict(
            sorted(Counter(record["occupied_sector_count"] for record in gated).items())
        ),
        "records": records,
        "scope_warning": (
            "Support-sector occupancy is a shipped-frame statement. It is used only "
            "when the attainer gate says the 1-RDM is diagonal in that frame, and it "
            "is not promoted to a basis-free fiber claim."
        ),
    }


def mean(values: list[int]) -> list[int]:
    if not values:
        return [0, 1]
    numerator = sum(values)
    denominator = len(values)
    divisor = math.gcd(numerator, denominator)
    return [numerator // divisor, denominator // divisor]


def summarize_system(records: list[dict]) -> dict[str, object]:
    nonstructural = [
        record for record in records if record["label"] != "STRUCTURAL"
    ]
    by_label: dict[str, dict] = {}
    for label in sorted({record["label"] for record in records}):
        selected = [record for record in records if record["label"] == label]
        by_label[label] = {
            "rows": len(selected),
            "fusion_rank_distribution": dict(
                sorted(Counter(record["fusion_rank"] for record in selected).items())
            ),
            "mean_fusion_rank": mean(
                [record["fusion_rank"] for record in selected]
            ),
            "single_sector_rows": sum(
                record["fusion_rank"] == 1 for record in selected
            ),
            "max_distinct_levels": max(
                record["distinct_levels"] for record in selected
            ),
        }

    oriented = {
        record["oriented_signature_sha256"] for record in nonstructural
    }
    unoriented = {
        record["unoriented_signature_sha256"] for record in nonstructural
    }
    primitive = {
        record["unoriented_signature_sha256"]
        for record in nonstructural
        if record["label"] == "PRIMITIVE"
    }
    inherited = {
        record["unoriented_signature_sha256"]
        for record in nonstructural
        if record["label"] in {"PADDING", "FROZEN-CORE"}
    }

    return {
        "rows": len(records),
        "nonstructural_rows": len(nonstructural),
        "fusion_rank_distribution": dict(
            sorted(Counter(record["fusion_rank"] for record in records).items())
        ),
        "max_fusion_rank": max(record["fusion_rank"] for record in records),
        "max_distinct_levels": max(
            record["distinct_levels"] for record in records
        ),
        "distinct_oriented_zero_grade_signatures": len(oriented),
        "distinct_unoriented_zero_grade_signatures": len(unoriented),
        "unoriented_signature_to_nonstructural_row_ratio": (
            [len(unoriented), len(nonstructural)]
            if nonstructural
            else [0, 1]
        ),
        "primitive_inherited_signature_overlap": len(primitive & inherited),
        "by_label": by_label,
    }


def score_heldout(systems: dict[str, dict]) -> dict[str, object]:
    heldout = {
        system: payload
        for system, payload in systems.items()
        if parse_system(system)[1] >= 9
        and payload.get("primary", False)
    }

    p1_violations = []
    p2_violations = []
    p3_violations = []
    p4_violations = []
    p5_violations = []
    p6_violations = []
    for system, payload in heldout.items():
        n, d = parse_system(system)
        summary = payload["summary"]
        if summary["max_distinct_levels"] > 2 * n:
            p1_violations.append(system)

        ratio_num, ratio_den = summary[
            "unoriented_signature_to_nonstructural_row_ratio"
        ]
        if 2 * ratio_num > ratio_den:
            p2_violations.append(system)

        by_label = summary["by_label"]
        primitive = by_label.get("PRIMITIVE")
        inherited_values = []
        for label in ("PADDING", "FROZEN-CORE"):
            entry = by_label.get(label)
            if entry:
                inherited_values.extend(
                    [record["fusion_rank"] for record in payload["rows"]]
                )
        if n == 3 and primitive:
            primitive_ranks = [
                record["fusion_rank"]
                for record in payload["rows"]
                if record["label"] == "PRIMITIVE"
            ]
            inherited_ranks = [
                record["fusion_rank"]
                for record in payload["rows"]
                if record["label"] in {"PADDING", "FROZEN-CORE"}
            ]
            if inherited_ranks and (
                Fraction(sum(primitive_ranks), len(primitive_ranks))
                <= Fraction(sum(inherited_ranks), len(inherited_ranks))
            ):
                p3_violations.append(system)

        if (
            any(
                record["label"] == "PRIMITIVE"
                for record in payload["rows"]
            )
            and any(
                record["label"] in {"PADDING", "FROZEN-CORE"}
                for record in payload["rows"]
            )
            and summary["primitive_inherited_signature_overlap"] == 0
        ):
            p4_violations.append(system)

        if n == 3 and summary["max_fusion_rank"] > d - 3:
            p5_violations.append(system)

        total_dimension = math.comb(d, n)
        if any(
            20 * record["zero_grade_dimension"] > 11 * total_dimension
            for record in payload["rows"]
            if record["label"] != "STRUCTURAL"
        ):
            p6_violations.append(system)

    def verdict(violations: list[str]) -> str:
        return "PASS" if not violations else "FAIL"

    return {
        "scope": sorted(heldout),
        "P1_distinct_tau_levels_at_most_2N": {
            "verdict": verdict(p1_violations),
            "violations": p1_violations,
        },
        "P2_zero_grade_mechanisms_compress_rows_by_at_least_2x": {
            "verdict": verdict(p2_violations),
            "evidence_role": "confirmatory_cross_realization_check",
            "violations": p2_violations,
        },
        "P3_N3_primitive_mean_fusion_exceeds_inherited": {
            "verdict": verdict(p3_violations),
            "violations": p3_violations,
        },
        "P4_primitive_and_inherited_mechanisms_overlap": {
            "verdict": verdict(p4_violations),
            "violations": p4_violations,
        },
        "P5_N3_max_fusion_rank_at_most_d_minus_3": {
            "verdict": verdict(p5_violations),
            "violations": p5_violations,
        },
        "P6_nonstructural_zero_grade_density_at_most_55_percent": {
            "verdict": verdict(p6_violations),
            "violations": p6_violations,
        },
        "status_note": (
            "P1 and P3-P6 were fixed after inspecting only d <= 8 pilot systems. "
            "P2 became confirmatory when an independent sorted-tau shape census "
            "landed on main before this audit was run at d = 9 and d = 10."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--systems",
        nargs="*",
        default=[],
        help="optional systems such as 3,7 3,8 4,8",
    )
    parser.add_argument("-o", "--out", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    ordered = json.loads(ORDERED.read_text())
    selected = {
        f"({value})" if not value.startswith("(") else value
        for value in args.systems
    }
    states, gate = load_state_data()
    core_map = ordered.get("cores", {})

    systems: dict[str, dict] = {}
    for system, row_entries in ordered["rows"].items():
        if selected and system not in selected:
            continue
        n, d = parse_system(system)
        records = []
        for row in row_entries:
            tau = tuple(int(value) for value in row["tau"])
            descriptor = descriptor_from_tau(tau, n)
            record = {
                "index": int(row["index"]),
                "label": row["label"],
                "tau": list(tau),
                "affine": row.get("affine"),
                **descriptor,
            }
            if system in core_map:
                record["core"] = core_map[system][int(row["index"])]
            if states:
                record["tight_state_analysis"] = tight_state_analysis(
                    system,
                    tau,
                    descriptor,
                    states,
                    gate,
                )
            records.append(record)

        systems[system] = {
            "particle_number": n,
            "rank": d,
            "primary": bool(ordered["systems"][system]["primary"]),
            "rows": records,
            "summary": summarize_system(records),
        }

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/levi_fusion_audit.py",
        "input": "results/data/ordered_representation_stability.json",
        "evidence": {
            "row_zero_grade_decomposition": "exact",
            "state_support_sector_occupancy": (
                "exact support arithmetic, conditional on the recorded attainer gate"
            ),
        },
        "definition": (
            "The fusion rank of tau is the number of bounded integer allocations "
            "n_a with sum n_a=N and sum h_a n_a=0 across the distinct tau "
            "eigenspaces. Each allocation contributes tensor_a wedge^{n_a} C^{m_a}."
        ),
        "systems": systems,
        "heldout_scorecard": score_heldout(systems),
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {output}")
    for system, data in systems.items():
        summary = data["summary"]
        print(
            f"{system}: rows={summary['rows']} "
            f"signatures={summary['distinct_unoriented_zero_grade_signatures']} "
            f"max_fusion={summary['max_fusion_rank']} "
            f"max_levels={summary['max_distinct_levels']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
