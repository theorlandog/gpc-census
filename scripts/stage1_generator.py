#!/usr/bin/env python3
"""Generate replayable first-principles Stage 1 constraint certificates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from gpc_census.generation import generate_3_6


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "data" / "stage1_ressayre_3_6.json"


def payload() -> dict[str, object]:
    generated = generate_3_6()
    enumeration = generated.enumeration
    return {
        "schema_version": 1,
        "system": f"({generated.n},{generated.d})",
        "status": "proved first-principles reproduction",
        "generator_inputs": {
            "representation": "wedge^3 C^6",
            "stored_constraints_used": False,
            "stored_vertices_used": False,
            "plethysm_inner_degree": 6,
        },
        "convention": {
            "ressayre_inequality": "sum_i h_i * lambda_i >= z",
            "negative_root_pair_i_j": "e_j-e_i for zero-based i<j",
            "stored_constraint": "sum_i coeffs_i * lambda_i <= rhs",
        },
        "enumeration": {
            "admissible_hyperplanes": enumeration.admissible_hyperplane_count,
            "oriented_candidates": enumeration.oriented_candidate_count,
            "trace_candidates": enumeration.trace_candidate_count,
            "ressayre_elements": len(enumeration.certificates),
        },
        "generated_system": {
            "inequalities": [row.as_dict() for row in generated.inequalities],
            "equalities": [row.as_dict() for row in generated.equalities],
            "structural_inequalities_on_affine_hull": [
                row.as_dict() for row in generated.structural_inequalities
            ],
        },
        "inner_outer_sandwich": {
            "plethysm_inner_point_count": generated.inner_point_count,
            "outer_vertices": [
                [str(value) for value in vertex] for vertex in generated.outer_vertices
            ],
            "every_outer_vertex_is_plethysm_certified": (
                generated.outer_vertices_in_inner_cloud
            ),
        },
        "schubert_cross_checks": [
            certificate.as_dict() for certificate in generated.schubert_cross_checks
        ],
        "ressayre_certificates": [
            certificate.as_dict() for certificate in enumeration.certificates
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output}: "
        f"{result['enumeration']['ressayre_elements']} Ressayre certificates"
    )


if __name__ == "__main__":
    main()
