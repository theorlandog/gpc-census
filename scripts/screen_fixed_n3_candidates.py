#!/usr/bin/env python3
"""Screen a supplied fixed-N=3 BDR candidate export exactly.

This command does not run the upstream candidate search.  It parses one
supplied ``brut_inequations`` export as syntax, primitive-normalizes and
deduplicates its rows, separates the structural Pauli lower bound, and records
an exact or explicitly unresolved outcome for every remaining row.

The resulting artifact does not claim candidate exhaustiveness, facetness,
irredundancy, or completeness of a GPC system.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from gpc_census.generation.candidate_pipeline import screen_bdr_export


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def build_payload(
    rank: int,
    export_path: Path,
    *,
    ressayre_attempts: int = 32,
    cyclic_cross_check: bool = False,
    include_modular: bool = False,
    symbolic_fallback: bool | None = None,
    symbolic_max_determinant_order: int | None = None,
    workers: int = 1,
    require_resolved: bool = False,
) -> dict[str, object]:
    """Build a deterministic candidate-screening artifact for one rank."""
    if type(require_resolved) is not bool:
        raise TypeError("require_resolved must be bool")
    data = export_path.read_bytes()
    result = screen_bdr_export(
        data.decode("utf-8"),
        expected_rank=rank,
        ressayre_attempts=ressayre_attempts,
        cyclic_cross_check=cyclic_cross_check,
        include_modular=include_modular,
        symbolic_fallback=symbolic_fallback,
        symbolic_max_determinant_order=symbolic_max_determinant_order,
        workers=workers,
        require_resolved=False,
    )
    resolved = result.all_candidates_resolved
    fallback_request = (
        "auto_rank_at_most_13"
        if symbolic_fallback is None
        else "enabled" if symbolic_fallback else "disabled"
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
            "cyclic_cross_check": cyclic_cross_check,
            "include_modular_cross_checks": include_modular,
            "symbolic_fallback_request": fallback_request,
            "symbolic_fallback_effective": result.symbolic_fallback_enabled,
            "symbolic_max_determinant_order": symbolic_max_determinant_order,
            "worker_processes": workers,
        },
        "operational_gate": {
            "require_resolved_requested": require_resolved,
            "passed": not require_resolved or resolved,
            "failure_exit_status": None if not require_resolved or resolved else 2,
        },
        "screening": result.as_dict(),
        "proof_status": {
            "normalization_and_deduplication": "exact",
            "full_dimension": (
                "proved_and_replayed"
                if result.full_dimension.is_full_dimensional
                else "not_established_by_certificate"
            ),
            "certified_rows": "replayed_exact_ressayre_certificates",
            "trace_and_symbolic_zero_rejections": "exact",
            "evaluation_unresolved_rows": "no_rejection_claim",
            "cyclic_schubert": "cross_validation_only_never_a_rejection_gate",
            "candidate_exhaustiveness": "not_established",
            "facetness": "not_established",
            "redundancy": "not_established",
            "system_completeness": "not_established",
        },
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="BDR Python export")
    parser.add_argument("--rank", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ressayre-attempts", type=int, default=32)
    parser.add_argument(
        "--cyclic-cross-check",
        action="store_true",
        help="run exact cyclic-Schubert cross-validation on every certified row",
    )
    parser.add_argument(
        "--include-modular",
        action="store_true",
        help="add modular cyclic replay; requires --cyclic-cross-check",
    )
    parser.add_argument(
        "--symbolic-fallback",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "enable or disable exact symbolic fallback after evaluation exhaustion; "
            "the default is automatic enablement for ranks at most 13"
        ),
    )
    parser.add_argument(
        "--symbolic-max-determinant-order",
        type=int,
        default=None,
        help="preserve larger symbolic determinants as unresolved without constructing them",
    )
    parser.add_argument(
        "--require-resolved",
        action="store_true",
        help="write the artifact but exit 2 if any candidate remains unresolved",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="parallel worker processes for independent candidate rows",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Write the artifact and return a process status suitable for automation."""
    parser = _parser()
    args = parser.parse_args(argv)
    if args.include_modular and not args.cyclic_cross_check:
        parser.error("--include-modular requires --cyclic-cross-check")
    payload = build_payload(
        args.rank,
        args.input,
        ressayre_attempts=args.ressayre_attempts,
        cyclic_cross_check=args.cyclic_cross_check,
        include_modular=args.include_modular,
        symbolic_fallback=args.symbolic_fallback,
        symbolic_max_determinant_order=args.symbolic_max_determinant_order,
        workers=args.workers,
        require_resolved=args.require_resolved,
    )
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if not payload["operational_gate"]["passed"]:
        unresolved = payload["screening"]["counts"]["evaluation_unresolved_rows"]
        print(
            f"{unresolved} candidate row(s) remain evaluation_unresolved; "
            f"artifact written to {args.output}",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
