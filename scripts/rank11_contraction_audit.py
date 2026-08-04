#!/usr/bin/env python3
"""Recompute the four open (3,11) contraction floors on normalized states.

The historical attack omitted normalization from its objective. That omission
does not affect zero-residual attainability, but a positive optimum is then a
distance to the cone of unnormalized 1-RDMs. This audit uses the corrected
scale-invariant objective, includes an attainable rank-11 control, and records
the exact trace-hyperplane projection associated with candidate 44's claimed
level-5 inequality.

The output is numerical evidence, not a non-attainability proof. In particular,
the exact candidate-44 projection is conditional geometry for the claimed
hyperplane. Proving the inequality still requires a Schubert coefficient or an
exact SOS dual certificate.

Run:
  uv run python scripts/rank11_contraction_audit.py --tries 20

Writes results/data/rank11_contraction.json.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from fractions import Fraction

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from contraction_attack import attack_state  # noqa: E402
from gpc_census.classify import classify_full  # noqa: E402

SETTLEMENT = ROOT / "docs/bracket_3_11_settlement.json"
OUTPUT = ROOT / "results/data/rank11_contraction.json"


def candidate44_geometry():
    """Exact projection onto lambda_2+...+lambda_5+lambda_11=2 at fixed trace."""
    target = [Fraction(1, 2)] * 5 + [Fraction(1, 12)] * 6
    selected = {1, 2, 3, 4, 10}
    normal = [Fraction(int(i in selected)) for i in range(11)]
    mean = sum(normal) / 11
    trace_normal = [value - mean for value in normal]
    normal_squared = sum(value * value for value in trace_normal)
    violation = sum(target[i] for i in selected) - 2
    step = violation / normal_squared
    projection = [
        value - step * direction
        for value, direction in zip(target, trace_normal)
    ]
    distance_squared = violation * violation / normal_squared
    return {
        "indices_one_based": [2, 3, 4, 5, 11],
        "rhs": Fraction(2),
        "target_violation": violation,
        "trace_hyperplane_normal_squared": normal_squared,
        "distance_squared": distance_squared,
        "projection": projection,
    }


def _fraction_strings(values):
    return [str(value) for value in values]


def _run_case(record, tries, seed):
    out = {
        "index": record["index"],
        "verdict": record["verdict"],
        "integer_form": record["integer_form"],
        "denominator": record["denominator"],
        "runs": {},
    }
    for real in (True, False):
        result = attack_state(
            record["integer_form"],
            record["denominator"],
            n=3,
            real=real,
            tries=tries,
            seed=seed,
        )
        rho = result["one_rdm"]
        label = "real" if real else "complex"
        out["runs"][label] = {
            "normalized_residual": result["residual"],
            "eigenvalues": [float(value) for value in result["eigenvalues"]],
            "support_threshold_count": result["support_size"],
            "state_norm_error": float(
                abs(np.linalg.norm(result["amplitudes"]) - 1.0)
            ),
            "one_rdm_trace_error": float(abs(np.trace(rho).real - 3.0)),
        }
    return out



def _was_open_before_level5(record):
    """The four candidates the pre-level-five methods could not settle.

    Those four now carry a LEVEL5-RESSAYRE refutation in the settlement, so
    keying off the OPEN verdict alone would find nothing. Accepting either
    keeps this script correct against both states of that artifact.
    """
    return record.get("verdict") == "OPEN" or record.get("certificate") == "LEVEL5-RESSAYRE"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tries", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    settlement = json.loads(SETTLEMENT.read_text())
    opens = [
        record for record in settlement["candidates"]
        if _was_open_before_level5(record)
    ]
    control = next(
        record for record in settlement["candidates"] if record["index"] == 49
    )
    rows = [_run_case(record, args.tries, args.seed) for record in opens]
    control_row = _run_case(control, args.tries, args.seed)

    geometry = candidate44_geometry()
    projected_classification = classify_full(3, 11, geometry["projection"])
    candidate44 = next(row for row in rows if row["index"] == 44)
    projected = np.array([float(value) for value in geometry["projection"]])
    for run in candidate44["runs"].values():
        run["projection_spectrum_max_error"] = float(
            np.max(np.abs(np.array(run["eigenvalues"]) - projected))
        )
        run["projection_distance_squared_error"] = float(
            abs(
                run["normalized_residual"]
                - float(geometry["distance_squared"])
            )
        )

    report = {
        "note": "Normalized ansatz-free contraction audit for the four OPEN "
        "(3,11) bracket candidates. Positive floors are numerical evidence, "
        "not proof. The candidate-44 exact projection is conditional on the "
        "CLAIMED level-5 inequality.",
        "objective": "min over nonzero c of ||gamma(c)/||c||^2-diag(target)||_F^2",
        "tries_per_real_or_complex_run": args.tries,
        "seed": args.seed,
        "attainable_control": control_row,
        "open_candidates": rows,
        "candidate44_claimed_facet_geometry": {
            "inequality": "lambda_2+lambda_3+lambda_4+lambda_5+lambda_11<=2",
            "status": "CLAIMED-NOT-PROVED",
            "indices_one_based": geometry["indices_one_based"],
            "rhs": str(geometry["rhs"]),
            "target_violation": str(geometry["target_violation"]),
            "trace_hyperplane_normal_squared": str(
                geometry["trace_hyperplane_normal_squared"]
            ),
            "distance_squared": str(geometry["distance_squared"]),
            "projection_spectrum": _fraction_strings(geometry["projection"]),
            "projection_design_classification": projected_classification,
            "projection_attainment_status": "NUMERICAL-ONLY",
        },
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "control_floor": control_row["runs"]["complex"]["normalized_residual"],
        "open_floors": {
            str(row["index"]): row["runs"]["complex"]["normalized_residual"]
            for row in rows
        },
        "candidate44_exact_projection_distance_squared": str(
            geometry["distance_squared"]
        ),
    }, indent=2))


if __name__ == "__main__":
    main()
