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

WHICH DIRECTION THE RANDOMNESS FAVOURS, because getting this backwards inverts
every conclusion. ``Is_Ram_contracted`` SEARCHES FOR A NON-CONTRACTION WITNESS.
``is_not_contracted`` samples matrices over ``Q[i]`` and returns True the moment
one has full column rank, breaking early; ``Is_Ram_contracted`` returns
``False`` as soon as any divisor yields such a witness, and returns ``True``
only after the search is exhausted without finding one. Therefore:

    False  =  a witness was found        =>  NOT birational, and this is the
                                             direction carrying evidence
    True   =  no witness found in the
              random search              =>  Monte Carlo only; absence of a
                                             witness is not a proof

Upstream's ``random_deep`` defaults to 1, a single trial, so a True verdict at
the default depth is very weak indeed. This audit therefore fixes
``random_deep`` explicitly and repeats the whole sweep across several
independent seeds, and reports the two directions separately:

- ``non_birational_witnessed``  any trial returned False. Evidence bearing.
- ``no_obstruction_found``      every trial returned True. Monte Carlo.

The exact arm, which would settle it, does not run in this environment; see
`results/data/bdr_external_benchmark.json` for the measured reason.
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

# The probabilistic arm searches for a non-contraction witness. Upstream's
# default depth is 1, which is a single trial, so a True verdict at that depth
# is very weak. These are the sampling parameters this audit uses instead, and
# they are serialized into the artifact rather than left implicit.
SEEDS = (20260809, 20260810, 20260811)
RANDOM_DEEP = 8


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


def external_verdicts(candidates: dict, candidates_sha256: str) -> dict:
    """Upstream linear-triangular and birationality verdicts, per row.

    The birationality arm is a randomized SEARCH FOR A NON-CONTRACTION WITNESS,
    so False is the evidence-bearing direction and True only means the search
    came up empty. Every row is therefore evaluated once per seed at an explicit
    depth, and both directions are recorded separately.
    """
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
        representations = []
        for seed in SEEDS:
            representation = FermionRepresentation(
                group,
                particle_cnt=record["particle_number"],
                seed=seed,
                random_deep=RANDOM_DEEP,
            )
            assert representation.random_deep == RANDOM_DEEP
            representation.T_Pi_3D  # noqa: B018  (precomputed, as the pipeline does)
            representations.append((seed, representation))

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
                    is_linear_triangular(
                        representations[0][1], inequality.tau, inversions
                    )
                )
                trials = []
                for seed, representation in representations:
                    trials.append(
                        {
                            "seed": seed,
                            "random_deep": RANDOM_DEEP,
                            "birational": bool(
                                Is_Ram_contracted(
                                    inequality,
                                    representation,
                                    "probabilistic",
                                    "probabilistic",
                                )
                            ),
                        }
                    )
                entry["trials"] = trials
                # False means a witness was found, which is the strong verdict.
                entry["non_birational_witnessed"] = any(
                    not trial["birational"] for trial in trials
                )
                entry["no_obstruction_found"] = all(
                    trial["birational"] for trial in trials
                )
                entry["bdr_error"] = None
            except Exception as error:  # noqa: BLE001  (record, never swallow)
                entry["inversion_count"] = None
                entry["bdr_linear_triangular"] = None
                entry["trials"] = []
                entry["non_birational_witnessed"] = None
                entry["no_obstruction_found"] = None
                entry["bdr_error"] = f"{type(error).__name__}: {error}"
            rows.append(entry)

        witnessed = sum(row["non_birational_witnessed"] is True for row in rows)
        unstable = sum(
            1
            for row in rows
            if row["trials"]
            and len({trial["birational"] for trial in row["trials"]}) > 1
        )
        systems[name] = {
            "system": name,
            "rows": rows,
            "row_count": len(rows),
            "non_birational_witnessed_rows": witnessed,
            "no_obstruction_found_rows": sum(
                row["no_obstruction_found"] is True for row in rows
            ),
            "seed_disagreement_rows": unstable,
            "bdr_linear_triangular_rows": sum(
                row["bdr_linear_triangular"] is True for row in rows
            ),
            "errored_rows": sum(row["bdr_error"] is not None for row in rows),
            "seconds": round(time.perf_counter() - started, 1),
        }
        print(
            f"{name}: {len(rows)} rows, "
            f"{witnessed} witnessed NOT birational, "
            f"{systems[name]['no_obstruction_found_rows']} no obstruction found, "
            f"{unstable} disagreeing across seeds, "
            f"{systems[name]['errored_rows']} errored "
            f"({systems[name]['seconds']}s)",
            flush=True,
        )

    return {
        "schema_version": 2,
        "generated_by": "scripts/bdr_birationality_audit.py --external-out",
        "candidates_sha256": candidates_sha256,
        "external_backend": {
            "name": "Bulois-Denis-Ressayre moment_cone",
            "mirror": "github.com/ea-icj/moment_cone",
            "pinned_revision": PINNED_REVISION,
            "license": "MIT",
            "entry_points": [
                "moment_cone.ramification.Is_Ram_contracted",
                "moment_cone.linear_triangular.is_linear_triangular",
            ],
            "birationality_arm": "probabilistic (the upstream default arm)",
            "seeds": list(SEEDS),
            "random_deep": RANDOM_DEEP,
            "upstream_default_random_deep": 1,
            "which_direction_is_evidence": (
                "False. Is_Ram_contracted searches for a non-contraction "
                "witness and returns False as soon as one is found, so False "
                "certifies non-birationality up to exact arithmetic while True "
                "only records that the random search came up empty."
            ),
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


COLUMNS = (
    "dm_fully_triangular",
    "bdr_linear_triangular",
    "no_obstruction_found",
)


def check_binding(candidates_bytes: bytes, external: dict) -> None:
    """Refuse to score anything unless the two halves provably match.

    Fail closed. A missing external row would otherwise be scored as unknown
    and dropped from the denominators, so an incomplete external run would
    quietly look better than a complete one.
    """
    recorded = external.get("candidates_sha256")
    actual = hashlib.sha256(candidates_bytes).hexdigest()
    if recorded != actual:
        raise SystemExit(
            "external verdicts were produced against a different candidate "
            f"file: recorded {recorded}, supplied {actual}. Re-run "
            "--external-out against this candidates artifact."
        )

    candidates = json.loads(candidates_bytes)
    if set(candidates["systems"]) != set(external["systems"]):
        raise SystemExit(
            "system sets differ: "
            f"{sorted(candidates['systems'])} against "
            f"{sorted(external['systems'])}"
        )

    for name, record in candidates["systems"].items():
        wanted = {tuple(row["tau"]) for row in record["rows"]}
        got = {tuple(row["tau"]) for row in external["systems"][name]["rows"]}
        if wanted != got:
            raise SystemExit(
                f"{name}: external rows do not cover the candidates. "
                f"missing {sorted(wanted - got)[:3]}, "
                f"unexpected {sorted(got - wanted)[:3]}"
            )
        for row in external["systems"][name]["rows"]:
            if row.get("bdr_error") is not None:
                raise SystemExit(f"{name}: external row {row['tau']} errored")
            if row.get("no_obstruction_found") is None:
                raise SystemExit(f"{name}: external row {row['tau']} has no verdict")
            if not row.get("trials"):
                raise SystemExit(f"{name}: external row {row['tau']} has no trials")


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
            # check_binding has already proved every tau is present, so a
            # KeyError here is a bug rather than a tolerable gap.
            outside = external_rows[tau]
            merged["bdr_linear_triangular"] = outside["bdr_linear_triangular"]
            merged["non_birational_witnessed"] = outside["non_birational_witnessed"]
            merged["no_obstruction_found"] = outside["no_obstruction_found"]
            merged["trials"] = outside["trials"]
            rows.append(merged)

        pooled.extend(rows)
        systems[name] = {
            "system": name,
            "particle_number": record["particle_number"],
            "orbitals": record["orbitals"],
            "candidate_rows": len(rows),
            "published_facets_among_them": sum(row["is_published_facet"] for row in rows),
            "published_facets_not_candidates": record["published_facets_not_candidates"],
            "witnessed_non_birational_rows": sum(
                row["non_birational_witnessed"] for row in rows
            ),
            "witnessed_non_birational_facets": sum(
                row["non_birational_witnessed"] and row["is_published_facet"]
                for row in rows
            ),
            "seed_disagreement_rows": sum(
                len({trial["birational"] for trial in row["trials"]}) > 1
                for row in rows
            ),
            "scores": {column: _score(rows, column) for column in COLUMNS},
            "rows": rows,
        }

    return {
        "systems": systems,
        "pooled_scores": {column: _score(pooled, column) for column in COLUMNS},
        "pooled_rows": len(pooled),
        "strong_direction": {
            "meaning": (
                "Is_Ram_contracted returns False only when it has found a "
                "non-contraction witness, so False is the direction that "
                "carries evidence and True only records an empty search."
            ),
            "witnessed_non_birational_rows": sum(
                row["non_birational_witnessed"] for row in pooled
            ),
            "witnessed_non_birational_facets": sum(
                row["non_birational_witnessed"] and row["is_published_facet"]
                for row in pooled
            ),
            "seed_disagreement_rows": sum(
                len({trial["birational"] for trial in row["trials"]}) > 1
                for row in pooled
            ),
        },
    }


def verdict(joined: dict) -> dict:
    pooled = joined["pooled_scores"]
    weak = pooled["no_obstruction_found"]
    strong = joined["strong_direction"]
    return {
        "question": (
            "Among determinant-nonzero Hall-feasible candidates, which of the "
            "three criteria separates a published facet from a non-facet?"
        ),
        "which_direction_is_evidence": (
            "False. Is_Ram_contracted searches for a non-contraction witness "
            "and returns False the moment one is found, so False certifies "
            "non-birationality up to exact arithmetic. True records only that "
            "the random search came up empty and is Monte Carlo."
        ),
        "strong_result": (
            f"{strong['witnessed_non_birational_rows']} rows are witnessed NOT "
            f"birational, of which {strong['witnessed_non_birational_facets']} "
            "are published facets. Those verdicts carry evidence."
        ),
        "monte_carlo_result": (
            f"{weak['tp']} facets and {weak['fp']} non-facets survived the "
            "search with no obstruction found. Necessity is therefore "
            "EMPIRICALLY SUPPORTED, not proved, and insufficiency is a "
            "CANDIDATE claim: some of those non-facet survivors could be "
            "false positives of the random search rather than genuinely "
            "birational rows."
        ),
        "sampling": (
            f"{len(SEEDS)} independent seeds at random_deep={RANDOM_DEEP}, "
            f"against the upstream default of 1. "
            f"{strong['seed_disagreement_rows']} rows disagreed across seeds, "
            "which measures how much the verdict still moves with sampling."
        ),
        "what_would_settle_it": (
            "the exact arm of Is_Ram_contracted, which does not run in this "
            "environment; see results/data/bdr_external_benchmark.json."
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
        candidates_bytes = arguments.candidates.read_bytes()
        candidates = json.loads(candidates_bytes)
        payload = external_verdicts(
            candidates, hashlib.sha256(candidates_bytes).hexdigest()
        )
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

    candidates_bytes = arguments.candidates.read_bytes()
    external = json.loads(arguments.external.read_text())
    check_binding(candidates_bytes, external)
    candidates = json.loads(candidates_bytes)
    joined = join(candidates, external)

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/bdr_birationality_audit.py",
        "status": "three_way_criterion_audit",
        "external_backend": external["external_backend"],
        "external_verdicts_sha256": hashlib.sha256(
            arguments.external.read_bytes()
        ).hexdigest(),
        "candidates_sha256": hashlib.sha256(candidates_bytes).hexdigest(),
        "binding": (
            "scored only after proving the external verdicts were produced "
            "against this exact candidates artifact, that the tau sets are "
            "equal, and that no row errored or lacks a verdict"
        ),
        "systems": joined["systems"],
        "pooled_scores": joined["pooled_scores"],
        "pooled_rows": joined["pooled_rows"],
        "strong_direction": joined["strong_direction"],
        "sampling": {
            "seeds": list(SEEDS),
            "random_deep": RANDOM_DEEP,
            "upstream_default_random_deep": 1,
        },
        "verdict": verdict(joined),
    }
    payload["canonical_sha256_without_this_field_or_timings"] = _canonical(payload)

    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    for name, record in sorted(payload["systems"].items()):
        scores = record["scores"]["no_obstruction_found"]
        print(
            f"{name}: {record['candidate_rows']} candidates, "
            f"{record['published_facets_among_them']} facets | "
            f"witnessed NOT birational {record['witnessed_non_birational_rows']} "
            f"(facets {record['witnessed_non_birational_facets']}) | "
            f"no obstruction found tp {scores['tp']} fp {scores['fp']} "
            f"fn {scores['fn']} tn {scores['tn']}",
            flush=True,
        )
    print(f"pooled (Monte Carlo) {payload['pooled_scores']['no_obstruction_found']}")
    print(f"strong direction {payload['strong_direction']}")
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
