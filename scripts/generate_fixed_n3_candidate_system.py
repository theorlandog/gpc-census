#!/usr/bin/env python3
"""Screen and exactly reduce one supplied fixed-N=3 BDR candidate export.

The command always writes a diagnostic artifact after successful screening.
Unresolved candidates block reduction and return exit status 2.  A retained
set mismatch against the published rank-7 through rank-10 system returns exit
status 3.  A successful rank-11 through rank-13 artifact remains conditional
on exhaustiveness of the supplied upstream candidate export.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from gpc_census.generation import candidate_pipeline
from gpc_census.generation.candidate_system import compose_candidate_system


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def build_payload(
    rank: int,
    export_path: Path,
    *,
    ressayre_attempts: int = 1,
    symbolic_fallback: bool | None = None,
    symbolic_max_determinant_order: int | None = None,
    workers: int = 1,
    include_modular: bool = False,
) -> dict[str, object]:
    """Screen, reduce, cross-check retained rows, and build one artifact."""
    data = export_path.read_bytes()
    screening = candidate_pipeline.screen_bdr_export(
        data.decode("utf-8"),
        expected_rank=rank,
        ressayre_attempts=ressayre_attempts,
        cyclic_cross_check=False,
        include_modular=False,
        symbolic_fallback=symbolic_fallback,
        symbolic_max_determinant_order=symbolic_max_determinant_order,
        workers=workers,
        require_resolved=False,
    )
    result = compose_candidate_system(
        screening,
        include_modular=include_modular,
    )
    fallback_request = (
        "auto_rank_at_most_13"
        if symbolic_fallback is None
        else "enabled"
        if symbolic_fallback
        else "disabled"
    )
    return {
        "schema_version": 1,
        "particle_number": 3,
        "rank": rank,
        "source_export": {
            "filename": export_path.name,
            "sha256": _sha256(data),
            "bytes": len(data),
        },
        "options": {
            "ressayre_evaluation_attempt_limit": ressayre_attempts,
            "symbolic_fallback_request": fallback_request,
            "symbolic_fallback_effective": screening.symbolic_fallback_enabled,
            "symbolic_max_determinant_order": symbolic_max_determinant_order,
            "worker_processes": workers,
            "candidate_wide_cyclic_cross_checks": False,
            "retained_row_cyclic_cross_checks": True,
            "include_modular_retained_cross_checks": include_modular,
        },
        "operational_gate": {
            "passed": result.finalized,
            "exit_status": result.exit_status,
            "status": result.status,
        },
        "candidate_system": result.as_dict(),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="BDR Python export")
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ressayre-attempts", type=int, default=1)
    parser.add_argument(
        "--symbolic-fallback",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "enable or disable exact symbolic fallback; by default it is enabled "
            "automatically through rank 13"
        ),
    )
    parser.add_argument(
        "--symbolic-max-determinant-order",
        type=int,
        default=None,
        help="leave larger symbolic determinants unresolved instead of expanding them",
    )
    parser.add_argument(
        "--include-modular",
        action="store_true",
        help="add modular cyclic-Schubert replay for retained rows only",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="parallel worker processes for independent candidate rows",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Write the diagnostic artifact and return its fail-closed status."""
    args = _parser().parse_args(argv)
    payload = build_payload(
        args.rank,
        args.input,
        ressayre_attempts=args.ressayre_attempts,
        symbolic_fallback=args.symbolic_fallback,
        symbolic_max_determinant_order=args.symbolic_max_determinant_order,
        workers=args.workers,
        include_modular=args.include_modular,
    )
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    exit_status = int(payload["operational_gate"]["exit_status"])
    if exit_status:
        print(
            f"candidate-system finalization blocked: "
            f"{payload['operational_gate']['status']}; "
            f"diagnostic artifact written to {args.output}",
            file=sys.stderr,
        )
    return exit_status


if __name__ == "__main__":
    raise SystemExit(main())
