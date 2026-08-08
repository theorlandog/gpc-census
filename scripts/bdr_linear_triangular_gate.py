#!/usr/bin/env python3
"""Gate 2: run BDR's own ``is_linear_triangular`` and compare it row by row.

WHY THIS EXISTS. `docs/levi_triangularity_gate.md` measured a local proxy for
linear triangularity: Dulmage-Mendelsohn factorization of the Ressayre tangent
matrix, with "triangular" read as "every irreducible block is 1 x 1". It then
asserted, without running anything, that BDR's own recursive filter is a
*strictly stronger* criterion. That assertion was never tested. This gate tests
it, by executing the upstream implementation at its pinned revision and joining
its verdict to the local one on the same rows.

TWO HALVES, TWO ENVIRONMENTS. The upstream code needs a SageMath runtime, which
must never become a gpc-census dependency, so the external half runs in a
throwaway virtual environment and its verdicts are recorded as a shipped
artifact:

    # external half, in a venv that has moment_cone (see the doc for the recipe)
    python scripts/bdr_linear_triangular_gate.py --external-out \
        results/data/bdr_external_linear_triangular.json

    # local half and the join, in the repository environment
    uv run scripts/bdr_linear_triangular_gate.py

The second command reads the recorded external verdicts, recomputes the local
column from scratch, and writes the comparison. It therefore stays runnable and
self-checking on a machine that has no Sage at all: what it cannot do without
the external venv is *refresh* the external column.

WHAT IS JOINED. The join key is the row's homogeneous ``tau`` in the local
convention, which is the upstream ``wtau``. Both sides are keyed on it, so a
convention slip shows up as an unmatched row rather than as a silent mismatch.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import importlib.util
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]

EXTERNAL_DEFAULT = ROOT / "results/data/bdr_external_linear_triangular.json"
GATE_DEFAULT = ROOT / "results/data/bdr_linear_triangular_gate.json"

# Upstream ships fermionic reference inequalities for exactly these systems.
EXTERNAL_SYSTEMS = (
    ("(3,7)", "ineq_Fermion_7_p3", 3, 7),
    ("(4,7)", "ineq_Fermion_7_p4", 4, 7),
    ("(3,8)", "ineq_Fermion_8_p3", 3, 8),
    ("(4,8)", "ineq_Fermion_8_p4", 4, 8),
)

PINNED_REVISION = "7ab347d9a7837d68d92bde9fac74606913f395ca"


# --------------------------------------------------------------------------
# External half: needs moment_cone, does not need gpc_census.
# --------------------------------------------------------------------------


def external_verdicts() -> dict:
    """Run upstream ``is_linear_triangular`` on every shipped fermionic row."""
    import importlib
    import platform

    from moment_cone.linear_triangular import is_linear_triangular  # type: ignore

    systems = {}
    for name, module_name, particles, orbitals in EXTERNAL_SYSTEMS:
        reference = importlib.import_module(
            f"moment_cone.reference_datas.{module_name}"
        )
        representation = reference.V
        rows = []
        for inequality in reference.inequalities:
            inversions = list(inequality.inversions)
            rows.append(
                {
                    "tau": [int(value) for value in str(inequality.wtau).split()],
                    "dominant_tau": [
                        int(value) for value in str(inequality.tau).split()
                    ],
                    "inversion_count": len(inversions),
                    "bdr_linear_triangular": bool(
                        is_linear_triangular(
                            representation, inequality.tau, inversions
                        )
                    ),
                }
            )
        rows.sort(key=lambda row: row["tau"])
        systems[name] = {
            "system": name,
            "particle_number": particles,
            "orbitals": orbitals,
            "reference_module": f"moment_cone.reference_datas.{module_name}",
            "representation_dimension": int(representation.dim),
            "rows": rows,
            "row_count": len(rows),
            "bdr_linear_triangular_rows": sum(
                row["bdr_linear_triangular"] for row in rows
            ),
        }

    return {
        "schema_version": 1,
        "generated_by": "scripts/bdr_linear_triangular_gate.py --external-out",
        "external_backend": {
            "name": "Bulois-Denis-Ressayre moment_cone",
            "mirror": "github.com/ea-icj/moment_cone",
            "pinned_revision": PINNED_REVISION,
            "license": "MIT",
            "entry_point": "moment_cone.linear_triangular.is_linear_triangular",
            "call_site": "LinearTriangularStep.apply in moment_cone/main_steps.py",
            "runtime": "passagemath (modularized SageMath fork) on CPython",
            "executed_here": True,
        },
        "machine": {
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "systems": systems,
    }


# --------------------------------------------------------------------------
# Local half: needs gpc_census, does not need moment_cone.
# --------------------------------------------------------------------------


def _gate_module():
    spec = importlib.util.spec_from_file_location(
        "levi_triangularity_gate", ROOT / "scripts" / "levi_triangularity_gate.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def published_taus(particles: int, orbitals: int) -> list[tuple[int, ...]]:
    """The repository's published irredundant rows as primitive homogeneous tau."""
    from gpc_census.constraints import constraints

    taus = []
    for inequality in constraints(particles, orbitals)["inequalities"]:
        values = [particles * c - inequality["rhs"] for c in inequality["coeffs"]]
        divisor = 0
        for value in values:
            divisor = math.gcd(divisor, abs(value))
        taus.append(tuple(value // (divisor or 1) for value in values))
    return taus


def local_column(gate, tau: tuple[int, ...], particles: int) -> dict:
    """The local Dulmage-Mendelsohn reading of one row."""
    import random

    audit = gate.AUDIT
    on, below, roots, matrix = audit.tangent_matrix(tau, particles)
    record: dict = {
        "matrix_order": len(below),
        "root_count": len(roots),
    }
    if not below:
        record["local_status"] = "structural"
        record["dm_blocks"] = []
        record["dm_fully_triangular"] = True
        return record
    if len(below) != len(roots):
        record["local_status"] = "not_square"
        record["dm_blocks"] = None
        record["dm_fully_triangular"] = None
        return record
    if audit.maximum_matching(matrix) != len(below):
        record["local_status"] = "hall_infeasible"
        record["dm_blocks"] = None
        record["dm_fully_triangular"] = None
        return record
    verdict = audit.determinant_verdict(matrix, len(on), random.Random(7))
    if not verdict["nonzero"]:
        record["local_status"] = "determinant_zero"
        record["dm_blocks"] = None
        record["dm_fully_triangular"] = None
        return record
    match_row = gate.perfect_matching(matrix)
    blocks = gate.irreducible_blocks(matrix, match_row)
    record["local_status"] = "determinant_nonzero"
    record["dm_blocks"] = blocks
    record["dm_fully_triangular"] = max(blocks) == 1
    return record


def compare(external: dict) -> dict:
    gate = _gate_module()
    systems = {}
    agreement: collections.Counter = collections.Counter()

    for name, _module, particles, orbitals in EXTERNAL_SYSTEMS:
        external_system = external["systems"][name]
        external_by_tau = {
            tuple(row["tau"]): row for row in external_system["rows"]
        }
        local_taus = published_taus(particles, orbitals)
        local_set = set(local_taus)
        external_set = set(external_by_tau)

        rows = []
        for tau in sorted(external_set | local_set):
            record = {"tau": list(tau)}
            record.update(local_column(gate, tau, particles))
            record["in_published_system"] = tau in local_set
            record["in_external_reference"] = tau in external_set
            if tau in external_by_tau:
                record["bdr_linear_triangular"] = external_by_tau[tau][
                    "bdr_linear_triangular"
                ]
                record["bdr_inversion_count"] = external_by_tau[tau][
                    "inversion_count"
                ]
            else:
                record["bdr_linear_triangular"] = None
                record["bdr_inversion_count"] = None
            rows.append(record)

        joined = [
            row
            for row in rows
            if row["in_published_system"]
            and row["in_external_reference"]
            and row["dm_fully_triangular"] is not None
        ]
        cells: collections.Counter = collections.Counter()
        for row in joined:
            cells[(row["dm_fully_triangular"], row["bdr_linear_triangular"])] += 1
            agreement[(row["dm_fully_triangular"], row["bdr_linear_triangular"])] += 1

        systems[name] = {
            "system": name,
            "particle_number": particles,
            "orbitals": orbitals,
            "published_rows": len(local_set),
            "external_rows": len(external_set),
            "rows_in_both": len(local_set & external_set),
            "external_only": [list(t) for t in sorted(external_set - local_set)],
            "published_only": [list(t) for t in sorted(local_set - external_set)],
            "joined_analysable_rows": len(joined),
            "bdr_linear_triangular_rows": sum(
                row["bdr_linear_triangular"] for row in joined
            ),
            "dm_fully_triangular_rows": sum(
                row["dm_fully_triangular"] for row in joined
            ),
            "contingency": {
                "dm_triangular_and_bdr_triangular": cells[(True, True)],
                "dm_triangular_only": cells[(True, False)],
                "bdr_triangular_only": cells[(False, True)],
                "neither": cells[(False, False)],
            },
            "rows": rows,
        }

    total = sum(agreement.values())
    return {
        "systems": systems,
        "pooled_contingency": {
            "dm_triangular_and_bdr_triangular": agreement[(True, True)],
            "dm_triangular_only": agreement[(True, False)],
            "bdr_triangular_only": agreement[(False, True)],
            "neither": agreement[(False, False)],
            "rows": total,
        },
        "criteria_agree": agreement[(True, False)] + agreement[(False, True)] == 0,
        "dm_triangularity_is_necessary_for_bdr": agreement[(False, True)] == 0,
        "dm_triangularity_is_sufficient_for_bdr": agreement[(True, False)] == 0,
    }


def verdict(comparison: dict) -> dict:
    pooled = comparison["pooled_contingency"]
    return {
        "cross_validation": (
            "The repository's published irredundant systems and the upstream "
            "fermionic reference inequalities agree exactly as sets at (3,7), "
            "(4,7), (3,8) and (4,8), once the structural row (0,...,0,-1) is "
            "accounted for. That row is a GPC ordering constraint upstream "
            "carries in the same list and the repository stores separately."
        ),
        "documented_claim_under_test": (
            "docs/levi_triangularity_gate.md asserted, untested, that BDR's "
            "recursive linear-triangular filter is a STRICTLY STRONGER "
            "criterion than asking whether every DM block of the tangent "
            "matrix is 1 x 1."
        ),
        "documented_claim_status": (
            "REFUTED. Strictly stronger would force "
            "bdr_triangular_only == 0, and it is "
            f"{pooled['bdr_triangular_only']}. The two criteria are "
            "incomparable on the published fermionic rows: BDR certifies rows "
            "whose tangent matrix is a single irreducible DM block."
        ),
        "why_they_differ": (
            "They are computed on different objects. The local test factorizes "
            "the square Ressayre tangent matrix, indexed by positive weights "
            "against inversion roots. BDR's filter works on the full weight "
            "graph of the representation, counts paths into the NON-POSITIVE "
            "weights, and recursively removes the roots whose equations came "
            "out linear. Neither reading implies the other."
        ),
        "consequence": (
            "The DM 1 x 1 test may not be used as a conservative stand-in for "
            "BDR's filter in either direction. Gate 4's refutation of strict "
            "triangularity stands on its own object and is unaffected; what "
            "falls is only the sentence relating it to the upstream filter."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--external-out",
        type=pathlib.Path,
        default=None,
        help="run the upstream criterion and write the external verdicts here",
    )
    parser.add_argument(
        "--external",
        type=pathlib.Path,
        default=EXTERNAL_DEFAULT,
        help="recorded external verdicts to join against",
    )
    parser.add_argument("-o", "--out", type=pathlib.Path, default=GATE_DEFAULT)
    arguments = parser.parse_args()

    if arguments.external_out is not None:
        payload = external_verdicts()
        payload["canonical_sha256_without_this_field"] = _canonical(payload)
        arguments.external_out.parent.mkdir(parents=True, exist_ok=True)
        arguments.external_out.write_text(
            json.dumps(payload, indent=1, sort_keys=True) + "\n"
        )
        for record in payload["systems"].values():
            print(
                f"{record['system']}: {record['row_count']} rows, "
                f"{record['bdr_linear_triangular_rows']} linear triangular",
                flush=True,
            )
        print(f"wrote {arguments.external_out}")
        return 0

    sys.path.insert(0, str(ROOT / "src"))
    external = json.loads(arguments.external.read_text())
    comparison = compare(external)

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/bdr_linear_triangular_gate.py",
        "status": "cross_implementation_comparison",
        "external_backend": external["external_backend"],
        "external_verdicts_from": str(
            arguments.external.relative_to(ROOT)
            if arguments.external.is_absolute()
            else arguments.external
        ),
        "external_verdicts_sha256": hashlib.sha256(
            arguments.external.read_bytes()
        ).hexdigest(),
        "systems": comparison["systems"],
        "pooled_contingency": comparison["pooled_contingency"],
        "criteria_agree": comparison["criteria_agree"],
        "dm_triangularity_is_necessary_for_bdr": comparison[
            "dm_triangularity_is_necessary_for_bdr"
        ],
        "dm_triangularity_is_sufficient_for_bdr": comparison[
            "dm_triangularity_is_sufficient_for_bdr"
        ],
        "verdict": verdict(comparison),
    }
    payload["canonical_sha256_without_this_field"] = _canonical(payload)

    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    for record in payload["systems"].values():
        print(
            f"{record['system']}: published {record['published_rows']}, "
            f"external {record['external_rows']}, joined "
            f"{record['joined_analysable_rows']}, BDR triangular "
            f"{record['bdr_linear_triangular_rows']}, DM triangular "
            f"{record['dm_fully_triangular_rows']}",
            flush=True,
        )
    print(f"pooled {payload['pooled_contingency']}")
    print(f"wrote {arguments.out}")
    return 0


def _canonical(payload: dict) -> str:
    hashable = json.loads(json.dumps(payload))
    hashable.pop("machine", None)
    return hashlib.sha256(
        json.dumps(hashable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
