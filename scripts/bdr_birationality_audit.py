#!/usr/bin/env python3
"""Phase 1: the three-way birationality audit, on candidates and on facets.

THE QUESTION. Gate 4 measured a local proxy for linear triangularity, and gate 2
showed that proxy is incomparable to BDR's own linear-triangular filter. Neither
is the thing that actually decides a facet. BDR's birationality test is, and it
had never been run here. This audit puts all three verdicts on the same rows and
asks which of them separates a facet from a mere candidate.

THE ROWS. The population is every Hall-feasible row with a nonzero Ressayre
tangent determinant at `(3,7)`, `(3,8)` and `(3,9)`: 48, 136 and 240 rows
against published irredundant systems of 4, 31 and 52. Those candidates have
passed the trace test, Hall feasibility and a nonzero determinant, and nothing
else. Being a published facet is the ground truth column, so each criterion can
be scored as a predictor rather than described.

THE COLUMNS.

- ``dm_fully_triangular``   every Dulmage-Mendelsohn block of the tangent
                            matrix is 1 x 1 (this repository, gate 4);
- ``bdr_linear_triangular`` the upstream recursive filter (gate 2);
- ``bdr_birational``        ``Is_Ram_contracted``, the upstream birationality
                            test, which is what its pipeline actually decides
                            on;
- ``is_published_facet``    ground truth.

THREE PASSES, TWO ENVIRONMENTS, because the upstream needs a SageMath runtime
that must never become a dependency of this package:

    # 1. the deterministic candidate list, in the repository environment
    uv run scripts/bdr_birationality_audit.py --candidates-out /tmp/cands.json

    # 2. the upstream verdicts, in a venv that has moment_cone
    #    (recipe in docs/bdr_linear_triangular_gate.md)
    python scripts/bdr_birationality_audit.py --candidates /tmp/cands.json \
        --external-out results/data/bdr_birationality_external.json

    # 3. the join, in the repository environment
    uv run scripts/bdr_birationality_audit.py

Pass 3 recomputes the local column from scratch and reads the external column
from the recorded artifact, so it runs on a machine with no Sage.

CAVEAT THAT TRAVELS WITH EVERY BDR VERDICT HERE. ``Is_Ram_contracted`` runs its
probabilistic arm, which is the upstream default; its exact arm does not run in
this environment, for the reason measured in
`results/data/bdr_external_benchmark.json`. A probabilistic birationality test
can in principle reject a valid row, so a `False` is weaker evidence than a
`True`.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]

CANDIDATES_DEFAULT = ROOT / "results/data/bdr_birationality_candidates.json"
EXTERNAL_DEFAULT = ROOT / "results/data/bdr_birationality_external.json"
AUDIT_DEFAULT = ROOT / "results/data/bdr_birationality_audit.json"

SYSTEMS = ((3, 7, "flats"), (3, 8, "flats"), (3, 9, "augmentation"))
PINNED_REVISION = "7ab347d9a7837d68d92bde9fac74606913f395ca"
SEED = 20260809


# --------------------------------------------------------------------------
# Pass 1 and 3: repository environment.
# --------------------------------------------------------------------------


def _gate_module():
    spec = importlib.util.spec_from_file_location(
        "levi_triangularity_gate", ROOT / "scripts" / "levi_triangularity_gate.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def published_taus(particles: int, orbitals: int) -> set[tuple[int, ...]]:
    from gpc_census.constraints import constraints

    taus = set()
    for inequality in constraints(particles, orbitals)["inequalities"]:
        values = [particles * c - inequality["rhs"] for c in inequality["coeffs"]]
        divisor = 0
        for value in values:
            divisor = math.gcd(divisor, abs(value))
        taus.add(tuple(value // (divisor or 1) for value in values))
    return taus


def candidate_rows(gate, particles: int, orbitals: int, source: str) -> list[dict]:
    """Determinant-nonzero Hall-feasible rows, with the local column attached."""
    published = published_taus(particles, orbitals)
    rows = []
    seen = set()
    for tau, matrix in gate.determinant_nonzero_candidate_rows(
        particles, orbitals, source
    ):
        if tau in seen:
            continue
        seen.add(tau)
        match_row = gate.perfect_matching(matrix)
        if match_row is None:
            continue
        blocks = gate.irreducible_blocks(matrix, match_row)
        rows.append(
            {
                "tau": list(tau),
                "matrix_order": len(matrix),
                "dm_blocks": blocks,
                "dm_fully_triangular": max(blocks) == 1,
                "is_published_facet": tau in published,
            }
        )

    # Published rows that never appear as determinant-nonzero candidates would
    # be invisible otherwise, and their absence is itself worth recording.
    covered = {tuple(row["tau"]) for row in rows}
    missing = sorted(published - covered)
    rows.sort(key=lambda row: row["tau"])
    return rows, [list(tau) for tau in missing]


def build_candidates() -> dict:
    gate = _gate_module()
    systems = {}
    for particles, orbitals, source in SYSTEMS:
        rows, missing = candidate_rows(gate, particles, orbitals, source)
        name = f"({particles},{orbitals})"
        systems[name] = {
            "system": name,
            "particle_number": particles,
            "orbitals": orbitals,
            "candidate_rows": len(rows),
            "published_facets_among_them": sum(row["is_published_facet"] for row in rows),
            "published_facets_not_candidates": missing,
            "rows": rows,
        }
        print(
            f"{name}: {len(rows)} candidates, "
            f"{systems[name]['published_facets_among_them']} of them published, "
            f"{len(missing)} published rows outside the candidate population",
            flush=True,
        )
    return {
        "schema_version": 1,
        "generated_by": "scripts/bdr_birationality_audit.py --candidates-out",
        "systems": systems,
    }


# --------------------------------------------------------------------------
# Pass 2: external environment.
# --------------------------------------------------------------------------


def external_verdicts(candidates: dict) -> dict:
    """Upstream linear-triangular and birationality verdicts, per row."""
    import platform
    import time

    from moment_cone import (  # type: ignore
        FermionRepresentation,
        Inequality,
        LinearGroup,
        Tau,
    )
    from moment_cone.linear_triangular import is_linear_triangular  # type: ignore
    from moment_cone.ramification import Is_Ram_contracted  # type: ignore

    systems = {}
    for name, record in candidates["systems"].items():
        group = LinearGroup([record["orbitals"]])
        representation = FermionRepresentation(
            group, particle_cnt=record["particle_number"], seed=SEED
        )
        representation.T_Pi_3D  # noqa: B018  (precomputed once, as the pipeline does)

        rows = []
        started = time.perf_counter()
        for row in record["rows"]:
            tau = tuple(row["tau"])
            entry = {"tau": list(tau)}
            try:
                inequality = Inequality.from_tau(Tau.from_flatten(list(tau), group))
                inversions = list(inequality.inversions)
                entry["inversion_count"] = len(inversions)
                entry["bdr_linear_triangular"] = bool(
                    is_linear_triangular(representation, inequality.tau, inversions)
                )
                entry["bdr_birational"] = bool(
                    Is_Ram_contracted(
                        inequality, representation, "probabilistic", "probabilistic"
                    )
                )
                entry["bdr_error"] = None
            except Exception as error:  # noqa: BLE001  (record, never swallow)
                entry["inversion_count"] = None
                entry["bdr_linear_triangular"] = None
                entry["bdr_birational"] = None
                entry["bdr_error"] = f"{type(error).__name__}: {error}"
            rows.append(entry)

        systems[name] = {
            "system": name,
            "rows": rows,
            "row_count": len(rows),
            "bdr_birational_rows": sum(row["bdr_birational"] is True for row in rows),
            "bdr_linear_triangular_rows": sum(
                row["bdr_linear_triangular"] is True for row in rows
            ),
            "errored_rows": sum(row["bdr_error"] is not None for row in rows),
            "seconds": round(time.perf_counter() - started, 1),
        }
        print(
            f"{name}: {len(rows)} rows, "
            f"{systems[name]['bdr_birational_rows']} birational, "
            f"{systems[name]['bdr_linear_triangular_rows']} linear triangular, "
            f"{systems[name]['errored_rows']} errored "
            f"({systems[name]['seconds']}s)",
            flush=True,
        )

    return {
        "schema_version": 1,
        "generated_by": "scripts/bdr_birationality_audit.py --external-out",
        "external_backend": {
            "name": "Bulois-Denis-Ressayre moment_cone",
            "mirror": "github.com/ea-icj/moment_cone",
            "pinned_revision": PINNED_REVISION,
            "license": "MIT",
            "entry_points": [
                "moment_cone.ramification.Is_Ram_contracted",
                "moment_cone.linear_triangular.is_linear_triangular",
            ],
            "birationality_arm": "probabilistic (the upstream default)",
            "executed_here": True,
        },
        "machine": {
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "systems": systems,
    }


# --------------------------------------------------------------------------
# Pass 3: the join and the scoring.
# --------------------------------------------------------------------------


def _score(rows: list[dict], column: str) -> dict:
    """Score one criterion as a predictor of facethood."""
    cells = {"tp": 0, "fp": 0, "fn": 0, "tn": 0, "unknown": 0}
    for row in rows:
        predicted = row.get(column)
        if predicted is None:
            cells["unknown"] += 1
            continue
        actual = row["is_published_facet"]
        if predicted and actual:
            cells["tp"] += 1
        elif predicted and not actual:
            cells["fp"] += 1
        elif not predicted and actual:
            cells["fn"] += 1
        else:
            cells["tn"] += 1
    decided = cells["tp"] + cells["fp"] + cells["fn"] + cells["tn"]
    facets = cells["tp"] + cells["fn"]
    predicted_true = cells["tp"] + cells["fp"]
    return {
        **cells,
        "decided_rows": decided,
        "recall_on_facets": (cells["tp"] / facets) if facets else None,
        "precision": (cells["tp"] / predicted_true) if predicted_true else None,
        "separates_facets_exactly": cells["fp"] == 0 and cells["fn"] == 0,
    }


def join(candidates: dict, external: dict) -> dict:
    gate = _gate_module()
    systems = {}
    pooled: list[dict] = []

    for name, record in candidates["systems"].items():
        external_rows = {
            tuple(row["tau"]): row for row in external["systems"][name]["rows"]
        }
        rows = []
        for row in record["rows"]:
            tau = tuple(row["tau"])
            # Recompute the local column rather than trust the candidate file.
            _on, below, roots, matrix = gate.AUDIT.tangent_matrix(
                tau, record["particle_number"]
            )
            match_row = gate.perfect_matching(matrix) if below else None
            blocks = (
                gate.irreducible_blocks(matrix, match_row)
                if match_row is not None
                else []
            )
            merged = {
                "tau": list(tau),
                "matrix_order": len(below),
                "dm_blocks": blocks,
                "dm_fully_triangular": bool(blocks) and max(blocks) == 1,
                "is_published_facet": row["is_published_facet"],
            }
            outside = external_rows.get(tau, {})
            merged["bdr_linear_triangular"] = outside.get("bdr_linear_triangular")
            merged["bdr_birational"] = outside.get("bdr_birational")
            merged["bdr_error"] = outside.get("bdr_error")
            rows.append(merged)

        pooled.extend(rows)
        systems[name] = {
            "system": name,
            "particle_number": record["particle_number"],
            "orbitals": record["orbitals"],
            "candidate_rows": len(rows),
            "published_facets_among_them": sum(row["is_published_facet"] for row in rows),
            "published_facets_not_candidates": record["published_facets_not_candidates"],
            "scores": {
                column: _score(rows, column)
                for column in (
                    "dm_fully_triangular",
                    "bdr_linear_triangular",
                    "bdr_birational",
                )
            },
            "rows": rows,
        }

    return {
        "systems": systems,
        "pooled_scores": {
            column: _score(pooled, column)
            for column in (
                "dm_fully_triangular",
                "bdr_linear_triangular",
                "bdr_birational",
            )
        },
        "pooled_rows": len(pooled),
    }


def verdict(joined: dict) -> dict:
    pooled = joined["pooled_scores"]
    birational = pooled["bdr_birational"]
    return {
        "question": (
            "Among determinant-nonzero Hall-feasible candidates, which of the "
            "three criteria separates a published facet from a non-facet?"
        ),
        "bdr_birationality": (
            f"recall {birational['recall_on_facets']}, precision "
            f"{birational['precision']}, false positives {birational['fp']}, "
            f"false negatives {birational['fn']}"
        ),
        "arm_caveat": (
            "Is_Ram_contracted ran its probabilistic arm, the upstream default; "
            "its exact arm does not run in this environment. A probabilistic "
            "birationality test can reject a valid row, so a False verdict is "
            "weaker evidence than a True one."
        ),
        "scope": (
            "Three systems, and the candidate population is not the full "
            "Ressayre-admissible population: it is already filtered by the "
            "trace test, Hall feasibility and a nonzero tangent determinant."
        ),
    }


def _canonical(payload: dict) -> str:
    hashable = json.loads(json.dumps(payload))
    hashable.pop("machine", None)
    for record in hashable.get("systems", {}).values():
        record.pop("seconds", None)
    return hashlib.sha256(
        json.dumps(hashable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates-out", type=pathlib.Path, default=None)
    parser.add_argument("--candidates", type=pathlib.Path, default=CANDIDATES_DEFAULT)
    parser.add_argument("--external-out", type=pathlib.Path, default=None)
    parser.add_argument("--external", type=pathlib.Path, default=EXTERNAL_DEFAULT)
    parser.add_argument("-o", "--out", type=pathlib.Path, default=AUDIT_DEFAULT)
    arguments = parser.parse_args()

    if arguments.external_out is not None:
        candidates = json.loads(arguments.candidates.read_text())
        payload = external_verdicts(candidates)
        payload["canonical_sha256_without_this_field_or_timings"] = _canonical(payload)
        arguments.external_out.parent.mkdir(parents=True, exist_ok=True)
        arguments.external_out.write_text(
            json.dumps(payload, indent=1, sort_keys=True) + "\n"
        )
        print(f"wrote {arguments.external_out}")
        return 0

    sys.path.insert(0, str(ROOT / "src"))

    if arguments.candidates_out is not None:
        payload = build_candidates()
        arguments.candidates_out.parent.mkdir(parents=True, exist_ok=True)
        arguments.candidates_out.write_text(
            json.dumps(payload, indent=1, sort_keys=True) + "\n"
        )
        print(f"wrote {arguments.candidates_out}")
        return 0

    candidates = json.loads(arguments.candidates.read_text())
    external = json.loads(arguments.external.read_text())
    joined = join(candidates, external)

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/bdr_birationality_audit.py",
        "status": "three_way_criterion_audit",
        "external_backend": external["external_backend"],
        "external_verdicts_sha256": hashlib.sha256(
            arguments.external.read_bytes()
        ).hexdigest(),
        "systems": joined["systems"],
        "pooled_scores": joined["pooled_scores"],
        "pooled_rows": joined["pooled_rows"],
        "verdict": verdict(joined),
    }
    payload["canonical_sha256_without_this_field_or_timings"] = _canonical(payload)

    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    for name, record in sorted(payload["systems"].items()):
        scores = record["scores"]["bdr_birational"]
        print(
            f"{name}: {record['candidate_rows']} candidates, "
            f"{record['published_facets_among_them']} facets, birational "
            f"tp {scores['tp']} fp {scores['fp']} fn {scores['fn']} "
            f"tn {scores['tn']}",
            flush=True,
        )
    print(f"pooled {payload['pooled_scores']['bdr_birational']}")
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
