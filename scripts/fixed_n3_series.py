#!/usr/bin/env python3
"""Build deterministic exact checks from fixed-N=3 tau exports.

The heavy BDR search is intentionally external to this script.  Each input is
its Python export for one rank.  The file is parsed as syntax, never executed.
Every nonstructural row is checked by an exact Ressayre determinant witness
and the exact cyclic-Schubert evaluator.  Every rank also receives an exact
ambient moment-cone full-dimension witness.  Finally, rational ordered-slice
certificates prove that every supplied final-system row is irredundant and a
relative facet after adjoining the Pauli and dominance inequalities.

The ordered-slice certificate is a polyhedral statement about the supplied
rows.  It does not by itself establish GPC validity, candidate exhaustiveness,
or completeness of the resulting GPC system.  The remaining upstream search
obligation is recorded explicitly in the artifact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from gpc_census.constraints import constraints
from gpc_census.generation import certify_ordered_slice_redundancy
from gpc_census.generation.bdr import parse_python_tau_export, tau_from_constraint
from gpc_census.generation.full_dimension import require_full_dimension
from gpc_census.generation.model import IntegralConstraint
from gpc_census.generation.series import certify_tau_rank_system


UPSTREAM_URL = "https://github.com/ea-icj/moment_cone"
UPSTREAM_REVISION = "7ab347d9a7837d68d92bde9fac74606913f395ca"


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _lookup_regression(d: int, generated_rows: list[dict[str, object]]) -> dict[str, object]:
    if d > 10:
        return {"available": False}
    table = constraints(3, d)
    generated_keys = {
        tau_from_constraint(IntegralConstraint(tuple(row["coeffs"]), int(row["rhs"])), 3)
        for row in generated_rows
    }
    lookup_keys = {
        tau_from_constraint(IntegralConstraint(tuple(row["coeffs"]), int(row["rhs"])), 3)
        for row in table["inequalities"]
    }
    return {
        "available": True,
        "lookup_source": table["source"],
        "homogeneous_oriented_sets_equal": generated_keys == lookup_keys,
        "generated_only": [list(row) for row in sorted(generated_keys - lookup_keys)],
        "lookup_only": [list(row) for row in sorted(lookup_keys - generated_keys)],
    }


def build_rank_payload(
    d: int,
    export_path: Path,
    *,
    include_modular: bool,
    ressayre_attempts: int,
) -> dict[str, object]:
    data = export_path.read_bytes()
    raw_taus = parse_python_tau_export(data.decode())
    if any(len(row) != d for row in raw_taus):
        raise ValueError(f"{export_path} does not contain rank-{d} rows")
    system = certify_tau_rank_system(
        raw_taus,
        3,
        include_modular=include_modular,
        ressayre_attempts=ressayre_attempts,
    )
    if system.d != d:
        raise AssertionError("normalized system rank changed")
    system_payload = system.as_dict()
    rows = [entry["constraint"] for entry in system_payload["inequalities"]]
    reduction = certify_ordered_slice_redundancy(
        tuple(row.constraint for row in system.inequalities),
        d,
    )
    expected_retained = tuple(range(len(system.inequalities)))
    if reduction.retained_indices != expected_retained:
        raise ValueError(
            "exact ordered-slice reduction removed purported final-system "
            f"rows at indices {list(reduction.redundant_indices)}"
        )
    return {
        "rank": d,
        "source_export": {
            "filename": export_path.name,
            "sha256": _sha256(data),
            "bytes": len(data),
        },
        "ambient_moment_cone_full_dimension_certificate": (require_full_dimension(d).as_dict()),
        "system": system_payload,
        "ordered_slice_redundancy_certificate": reduction.as_dict(),
        "known_system_regression": _lookup_regression(d, rows),
        "proof_status": {
            "row_normalization": "proved_and_replayable",
            "full_dimension": "proved_and_replayable",
            "cyclic_schubert_coefficient_one": "proved_and_replayable",
            "ressayre_admissibility_trace_determinant": "proved_and_replayable",
            "exact_redundancy_elimination": "proved_and_replayable_for_supplied_rows",
            "ordered_slice_relative_facetness": "proved_and_replayable_for_supplied_rows",
            "candidate_exhaustiveness": "external_BDR_search_obligation",
            "system_completeness": "not_proved_here",
            "polyhedral_certificate_scope": (
                "does_not_by_itself_establish_GPC_validity_or_completeness"
            ),
        },
    }


def build_payload(
    inputs: list[tuple[int, Path]],
    *,
    include_modular: bool = False,
    ressayre_attempts: int = 32,
) -> dict[str, object]:
    ranks = [d for d, _path in inputs]
    if len(ranks) != len(set(ranks)):
        raise ValueError("each rank may be supplied only once")
    return {
        "schema_version": 3,
        "particle_number": 3,
        "upstream_generator": {
            "url": UPSTREAM_URL,
            "revision": UPSTREAM_REVISION,
            "license": "MIT",
        },
        "scope": {
            "ranks": sorted(ranks),
            "include_modular_cross_checks": include_modular,
            "ressayre_evaluation_attempt_limit": ressayre_attempts,
        },
        "rank_systems": [
            build_rank_payload(
                d,
                path,
                include_modular=include_modular,
                ressayre_attempts=ressayre_attempts,
            )
            for d, path in sorted(inputs)
        ],
    }


def _parse_input(values: list[list[str]]) -> list[tuple[int, Path]]:
    result: list[tuple[int, Path]] = []
    for rank_text, path_text in values:
        result.append((int(rank_text), Path(path_text)))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        nargs=2,
        action="append",
        metavar=("RANK", "EXPORT"),
        required=True,
        help="rank and BDR Python export path; repeat for every rank",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--include-modular", action="store_true")
    parser.add_argument("--ressayre-attempts", type=int, default=32)
    args = parser.parse_args()

    payload = build_payload(
        _parse_input(args.input),
        include_modular=args.include_modular,
        ressayre_attempts=args.ressayre_attempts,
    )
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
