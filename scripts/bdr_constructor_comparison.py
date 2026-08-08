#!/usr/bin/env python3
"""Measure the oracle ceiling on any external candidate generator.

Pre-registration: ``docs/prereg_bdr_constructor_comparison.md``.
Narrative: ``docs/bdr_constructor_comparison.md``.
Artifact: ``results/data/bdr_constructor_comparison.json``.

The Bulois-Denis-Ressayre implementation cannot be executed here, so this
script does not benchmark it. It measures the quantity that bounds it.

The standing acceptance rule is that every row entering a published system
carries a locally replayed exact certificate. So the best case for ANY
external candidate generator is that it supplies exactly the minimal facet
rows for free. This script measures two paths per system:

``full``    the local orbit-native construction, stage attributed into
            enumeration, screening and composition;
``oracle``  local screening plus composition when the published minimal rows
            are supplied directly, with no local candidate generation at all.

Their ratio is an upper bound on the speedup available from replacing local
candidate generation by an external generator whose rows are still certified
locally. A perfect oracle is strictly stronger than BDR, which must emit at
least some rows this repository then rejects, so the ratio bounds BDR from
above without BDR being run.

The oracle path is also the Task 3 recertification ledger: every published row
is converted into the repository's primitive ``tau`` convention and put
through the local exact verifier, and no external verdict is accepted by
citation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from gpc_census.constraints import constraints  # noqa: E402
from gpc_census.generation.bdr import primitive_tau, tau_from_constraint  # noqa: E402
from gpc_census.generation.candidate_pipeline import screen_tau_candidates  # noqa: E402
from gpc_census.generation.candidate_system import compose_candidate_system  # noqa: E402
from gpc_census.generation.model import IntegralConstraint  # noqa: E402
from gpc_census.generation.orbit_native import (  # noqa: E402
    OrbitNativeGeneration,
    enumerate_orbit_survivor_taus,
)
from gpc_census.generation.ressayre import admissible_candidate_from_tau  # noqa: E402

DEFAULT_SYSTEMS = ((3, 6), (3, 7), (3, 8))
SCHEMA_VERSION = 1

# Pinned before the run; see the preregistration for the availability record.
EXTERNAL_PIN: dict[str, Any] = {
    "name": "Bulois-Denis-Ressayre moment-cone algorithm",
    "paper": "An Algorithm to compute the Kronecker cone and other moment cones",
    "authors": [
        "Michael Bulois (ICJ, AGL)",
        "Roland Denis (INSMI-CNRS, ICJ)",
        "Nicolas Ressayre (ICJ, AGL)",
    ],
    "preprint": "arXiv:2505.08812 (math.AG)",
    "hal_record": "hal-05044759, version v3",
    "code_repository": "ea-icj/moment_cone at plmlab.math.cnrs.fr",
    "landing_page": "https://ea-icj.github.io/",
    "pinned_revision": "7ab347d9a7837d68d92bde9fac74606913f395ca",
    "runtime_environment": "Python-Sage (SageMath)",
    "license": "not_verified_code_host_unreachable",
    "relevant_case": "fermionic cone, GL_d(C) acting on wedge^r C^d, here r = 3",
    "executed_here": False,
    "why_not_executed": [
        "plmlab.math.cnrs.fr denied by the session egress proxy (403 CONNECT)",
        "ea-icj.github.io, arxiv.org, cnrs.hal.science denied by the same policy",
        "moment-cone is not published on PyPI",
        "SageMath is not installed and must never be a package runtime dependency",
    ],
    "benchmark_label": "INCOMPLETE",
}


def _canonical_hash(payload: object) -> str:
    """Return a stable sha256 over a canonical JSON encoding."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _orbit_key(n: int, tau: tuple[int, ...]) -> tuple[tuple[int, ...], int]:
    """Return the unoriented orbit invariant of a row.

    The ``S_d`` orbit of an admissible hyperplane is determined by the multiset
    of its primitive normal entries together with its affine level, so the
    sorted pair ``(h, z)`` identifies the orbit a row belongs to.
    """
    candidate = admissible_candidate_from_tau(n, tau)
    return tuple(sorted(candidate.h)), candidate.z


def _orbit_keys_of_row(
    n: int, tau: tuple[int, ...]
) -> tuple[tuple[tuple[int, ...], int], ...]:
    """Return both oriented keys of the unoriented hyperplane behind ``tau``.

    Each unoriented hyperplane carries two taus, giving keys ``(h, z)`` and
    ``(-h, -z)``. The orbit enumerator stores one orientation per unoriented
    orbit and its choice is not a documented convention, so matching tests both
    keys for membership instead of imposing a canonical sign here. Picking one
    by an arbitrary rule silently fails to match whenever the enumerator kept
    the other, which is what a min-of-the-two rule did at rank 6.
    """
    forward = _orbit_key(n, tau)
    backward = _orbit_key(n, primitive_tau(tuple(-value for value in tau)))
    return (forward, backward)


def _published_rows(n: int, d: int) -> tuple[dict[str, Any], ...]:
    """Return the published nonstructural rows for one system."""
    table = constraints(n, d)
    return tuple(
        {
            "coeffs": tuple(int(value) for value in row["coeffs"]),
            "rhs": int(row["rhs"]),
            "source": table["source"],
        }
        for row in table["inequalities"]
    )


def _recertification_ledger(
    n: int,
    d: int,
    screening: Any,
    orbit_keys: dict[tuple[tuple[int, ...], int], int],
) -> list[dict[str, Any]]:
    """Bind each published row to its converted tau and local verdict.

    Task 3 of the brief: external row, converted tau, local status, local
    certificate hash, matched local orbit. The external verdict is never
    accepted by citation, so the status recorded here is the one the local
    exact screen produced.
    """
    by_tau = {row.tau: row for row in screening.candidates}
    ledger: list[dict[str, Any]] = []
    for published in _published_rows(n, d):
        constraint = IntegralConstraint(published["coeffs"], published["rhs"])
        tau = tau_from_constraint(constraint, n)
        row = by_tau.get(tau)
        keys = _orbit_keys_of_row(n, tau)
        matched_key = next((key for key in keys if key in orbit_keys), None)
        certificate = (
            None
            if row is None or row.ressayre_assessment is None
            else row.ressayre_assessment.certificate
        )
        ledger.append(
            {
                "external_row": {
                    "coeffs": list(published["coeffs"]),
                    "rhs": published["rhs"],
                    "published_source": published["source"],
                },
                "converted_tau": list(tau),
                "local_status": "row_absent_from_screen" if row is None else row.status,
                "local_resolution_method": None if row is None else row.resolution_method,
                "local_certificate_replayed_and_bound": (
                    None if row is None else row.certificate_replayed_and_bound
                ),
                "local_certificate_sha256": (
                    None if certificate is None else _canonical_hash(certificate.as_dict())
                ),
                "matched_local_orbit": {
                    "sorted_h": None if matched_key is None else list(matched_key[0]),
                    "z": None if matched_key is None else matched_key[1],
                    "orbit_index": (
                        None if matched_key is None else orbit_keys[matched_key]
                    ),
                    "matched": matched_key is not None,
                    "candidate_keys": [
                        {"sorted_h": list(key[0]), "z": key[1]} for key in keys
                    ],
                },
            }
        )
    return ledger


def _stage_counts(generation: OrbitNativeGeneration) -> dict[str, Any]:
    """Return the stage-by-stage counts of one full construction."""
    enumeration = generation.enumeration
    screening = generation.screening
    counts = screening.counts
    system = generation.candidate_system
    return {
        "unoriented_orbits": len(enumeration.orbits),
        "admissible_hyperplanes": enumeration.hyperplane_count,
        "oriented_rows_total": enumeration.oriented_row_count,
        "trace_rejected_never_materialised": enumeration.trace_rejected_count,
        "survivors_generated": len(enumeration.survivor_taus),
        "survivors_predicted_in_closed_form": enumeration.predicted_survivor_count,
        "structural_rows": counts.structural_rows,
        "screened_candidate_rows": counts.nonstructural_candidate_rows,
        "certified_rows": counts.certified_rows,
        "trace_rejected_rows": counts.trace_rejected_rows,
        "symbolic_zero_rejected_rows": counts.symbolic_zero_rejected_rows,
        "inadmissible_rows": counts.inadmissible_rows,
        "unresolved_rows": counts.evaluation_unresolved_rows,
        "symbolic_fallback_attempted_rows": counts.symbolic_fallback_attempted_rows,
        "redundancy_input_rows": system.counts.redundancy_input_rows,
        "redundant_rows": system.counts.redundant_rows,
        "retained_rows": system.counts.retained_rows,
    }


def _regression(system: Any) -> dict[str, Any]:
    """Return the published-system regression verdict of one composition."""
    regression = system.known_system_regression
    if regression is None:
        return {"available": False, "reason": "composition blocked before regression"}
    return {
        "available": regression.available,
        "expected_count": regression.expected_count,
        "retained_count": regression.retained_count,
        "oriented_sets_equal": regression.oriented_sets_equal,
        "generated_only": [list(tau) for tau in regression.generated_only],
        "published_only": [list(tau) for tau in regression.published_only],
        "lookup_source": regression.lookup_source,
        "passed": regression.passed,
    }


def _measure_screening_only(n: int, d: int, *, workers: int) -> dict[str, Any]:
    """Screen a system whose rank the composition layer does not accept.

    Rank 6 still answers the recertification question, which is whether the
    published rows certify locally. It cannot answer the ceiling question,
    because there is no composition stage to time.
    """
    enumeration_start = time.perf_counter()
    enumeration = enumerate_orbit_survivor_taus(n, d)
    enumeration_seconds = time.perf_counter() - enumeration_start
    if not enumeration.counts_are_conserved:
        raise AssertionError("orbit-native enumeration failed count conservation")

    screen_start = time.perf_counter()
    full_screening = screen_tau_candidates(
        enumeration.survivor_taus,
        particle_number=n,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=workers,
    )
    full_screening_seconds = time.perf_counter() - screen_start

    orbit_keys = {
        (tuple(sorted(record.h)), record.z): index
        for index, record in enumerate(enumeration.orbits)
    }
    oracle_taus = tuple(
        tau_from_constraint(IntegralConstraint(row["coeffs"], row["rhs"]), n)
        for row in _published_rows(n, d)
    )
    oracle_start = time.perf_counter()
    oracle_screening = screen_tau_candidates(
        oracle_taus,
        particle_number=n,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=workers,
    )
    oracle_screening_seconds = time.perf_counter() - oracle_start

    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "composition_supported": False,
        "composition_exclusion_reason": (
            "compose_candidate_system guards ranks 7 through 13; rank 6 is "
            "covered by the separate exhaustive reference gate in "
            "docs/stage1_ressayre_3_6.md, so this entry is screening only and "
            "is excluded from the ceiling"
        ),
        "full_path": {
            "method": "orbit_native_canonical_augmentation_screening_only",
            "counts": {
                "unoriented_orbits": len(enumeration.orbits),
                "admissible_hyperplanes": enumeration.hyperplane_count,
                "oriented_rows_total": enumeration.oriented_row_count,
                "trace_rejected_never_materialised": enumeration.trace_rejected_count,
                "survivors_generated": len(enumeration.survivor_taus),
                "survivors_predicted_in_closed_form": (
                    enumeration.predicted_survivor_count
                ),
                "structural_rows": full_screening.counts.structural_rows,
                "screened_candidate_rows": (
                    full_screening.counts.nonstructural_candidate_rows
                ),
                "certified_rows": full_screening.counts.certified_rows,
                "trace_rejected_rows": full_screening.counts.trace_rejected_rows,
                "symbolic_zero_rejected_rows": (
                    full_screening.counts.symbolic_zero_rejected_rows
                ),
                "inadmissible_rows": full_screening.counts.inadmissible_rows,
                "unresolved_rows": full_screening.counts.evaluation_unresolved_rows,
            },
            "seconds": {
                "enumeration": enumeration_seconds,
                "screening": full_screening_seconds,
                "composition": None,
                "total": enumeration_seconds + full_screening_seconds,
            },
        },
        "oracle_path": {
            "method": "published_minimal_rows_supplied_as_candidates_screening_only",
            "supplied_rows": len(oracle_taus),
            "counts": {
                "structural_rows": oracle_screening.counts.structural_rows,
                "screened_candidate_rows": (
                    oracle_screening.counts.nonstructural_candidate_rows
                ),
                "certified_rows": oracle_screening.counts.certified_rows,
                "trace_rejected_rows": oracle_screening.counts.trace_rejected_rows,
                "symbolic_zero_rejected_rows": (
                    oracle_screening.counts.symbolic_zero_rejected_rows
                ),
                "inadmissible_rows": oracle_screening.counts.inadmissible_rows,
                "unresolved_rows": oracle_screening.counts.evaluation_unresolved_rows,
            },
            "seconds": {
                "enumeration": 0.0,
                "screening": oracle_screening_seconds,
                "composition": None,
                "total": oracle_screening_seconds,
            },
            "known_system_regression": {
                "available": False,
                "reason": "composition layer does not accept this rank",
            },
            "recertification_ledger": _recertification_ledger(
                n, d, oracle_screening, orbit_keys
            ),
        },
        "ceiling": None,
    }


COMPOSITION_MIN_RANK = 7
COMPOSITION_MAX_RANK = 13


def _composition_supported(d: int) -> bool:
    """Whether ``compose_candidate_system`` accepts this rank.

    The composition layer guards ranks 7 through 13. Rank 6 is covered instead
    by the separate exhaustive reference gate in ``docs/stage1_ressayre_3_6.md``,
    so a rank-6 entry here is screening only and is excluded from the ceiling
    rather than quietly dropped.
    """
    return COMPOSITION_MIN_RANK <= d <= COMPOSITION_MAX_RANK


def measure_system(n: int, d: int, *, workers: int = 1) -> dict[str, Any]:
    """Measure the full and oracle paths for one system."""
    if not _composition_supported(d):
        return _measure_screening_only(n, d, workers=workers)
    full_start = time.perf_counter()
    enumeration = enumerate_orbit_survivor_taus(n, d)
    enumeration_seconds = time.perf_counter() - full_start
    if not enumeration.counts_are_conserved:
        raise AssertionError("orbit-native enumeration failed count conservation")

    screen_start = time.perf_counter()
    full_screening = screen_tau_candidates(
        enumeration.survivor_taus,
        particle_number=n,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=workers,
    )
    full_screening_seconds = time.perf_counter() - screen_start

    compose_start = time.perf_counter()
    full_system = compose_candidate_system(full_screening)
    full_compose_seconds = time.perf_counter() - compose_start
    full_seconds = time.perf_counter() - full_start

    full_generation = OrbitNativeGeneration(
        enumeration=enumeration,
        screening=full_screening,
        candidate_system=full_system,
    )

    # Orbits are keyed directly by their own primitive (h, z) pair, which is
    # the same invariant the row-side key computes from a tau.
    orbit_keys = {
        (tuple(sorted(record.h)), record.z): index
        for index, record in enumerate(enumeration.orbits)
    }

    oracle_taus = tuple(
        tau_from_constraint(IntegralConstraint(row["coeffs"], row["rhs"]), n)
        for row in _published_rows(n, d)
    )
    oracle_start = time.perf_counter()
    oracle_screening = screen_tau_candidates(
        oracle_taus,
        particle_number=n,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=workers,
    )
    oracle_screening_seconds = time.perf_counter() - oracle_start
    oracle_compose_start = time.perf_counter()
    oracle_system = compose_candidate_system(oracle_screening)
    oracle_compose_seconds = time.perf_counter() - oracle_compose_start
    oracle_seconds = time.perf_counter() - oracle_start

    # Third path, added AFTER the predictions were scored and labelled as such
    # in the artifact. It models the one step BDR has that this repository does
    # not: a birationality prefilter that removes determinant-zero rows before
    # they are screened. Candidate generation is left identical to the local
    # path, so this isolates the prefilter and nothing else. It is an upper
    # bound on that step, because a real prefilter cannot be free.
    prefilter_taus = tuple(
        row.tau for row in full_screening.candidates if row.status == "certified"
    )
    prefilter_start = time.perf_counter()
    prefilter_screening = screen_tau_candidates(
        prefilter_taus,
        particle_number=n,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=workers,
    )
    prefilter_screening_seconds = time.perf_counter() - prefilter_start
    prefilter_compose_start = time.perf_counter()
    prefilter_system = compose_candidate_system(prefilter_screening)
    prefilter_compose_seconds = time.perf_counter() - prefilter_compose_start
    prefilter_total = (
        enumeration_seconds + prefilter_screening_seconds + prefilter_compose_seconds
    )

    ledger = _recertification_ledger(n, d, oracle_screening, orbit_keys)

    full_taus = {
        tuple(row.tau) for row in full_system.retained_rows
    }
    oracle_row_taus = {tuple(row.tau) for row in oracle_system.retained_rows}

    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "composition_supported": True,
        "full_path": {
            "method": "orbit_native_canonical_augmentation",
            "counts": _stage_counts(full_generation),
            "seconds": {
                "enumeration": enumeration_seconds,
                "screening": full_screening_seconds,
                "composition": full_compose_seconds,
                "total": full_seconds,
            },
            "enumeration_share_of_total": enumeration_seconds / full_seconds,
            "status": full_system.status,
            "proved_complete": full_generation.proved_complete,
            "known_system_regression": _regression(full_system),
        },
        "oracle_path": {
            "method": "published_minimal_rows_supplied_as_candidates",
            "supplied_rows": len(oracle_taus),
            "counts": {
                "structural_rows": oracle_screening.counts.structural_rows,
                "screened_candidate_rows": (
                    oracle_screening.counts.nonstructural_candidate_rows
                ),
                "certified_rows": oracle_screening.counts.certified_rows,
                "trace_rejected_rows": oracle_screening.counts.trace_rejected_rows,
                "symbolic_zero_rejected_rows": (
                    oracle_screening.counts.symbolic_zero_rejected_rows
                ),
                "inadmissible_rows": oracle_screening.counts.inadmissible_rows,
                "unresolved_rows": (
                    oracle_screening.counts.evaluation_unresolved_rows
                ),
                "redundancy_input_rows": oracle_system.counts.redundancy_input_rows,
                "redundant_rows": oracle_system.counts.redundant_rows,
                "retained_rows": oracle_system.counts.retained_rows,
            },
            "seconds": {
                "enumeration": 0.0,
                "screening": oracle_screening_seconds,
                "composition": oracle_compose_seconds,
                "total": oracle_seconds,
            },
            "status": oracle_system.status,
            "known_system_regression": _regression(oracle_system),
            "recertification_ledger": ledger,
        },
        "birationality_prefilter_path": {
            "added_after_scoring": True,
            "method": "local_enumeration_plus_perfect_determinant_zero_prefilter",
            "models": (
                "the one BDR step with no local counterpart: a birationality "
                "filter that rejects determinant-zero rows before they are "
                "screened. Local candidate generation is unchanged, so this "
                "isolates the prefilter."
            ),
            "rows_screened": len(prefilter_taus),
            "rows_removed_by_prefilter": (
                full_screening.counts.nonstructural_candidate_rows
                - len(prefilter_taus)
            ),
            "seconds": {
                "enumeration": enumeration_seconds,
                "screening": prefilter_screening_seconds,
                "composition": prefilter_compose_seconds,
                "total": prefilter_total,
            },
            "speedup_over_full": full_seconds / prefilter_total,
            "speedup_if_external_enumeration_were_also_free": (
                full_seconds
                / (prefilter_screening_seconds + prefilter_compose_seconds)
            ),
            "why_that_is_the_non_circular_bound": (
                "an external generator cannot know which certified rows survive "
                "reduction without running the reduction, so it must emit at "
                "least the certified rows and this repository must screen and "
                "reduce them. The oracle ceiling is below this only because the "
                "oracle is handed the 31 final rows, which is circular at a new "
                "rank."
            ),
            "retained_rows": prefilter_system.counts.retained_rows,
            "known_system_regression": _regression(prefilter_system),
            "caveat": (
                "an upper bound on the prefilter step alone; it leaves "
                "enumeration and composition untouched, and a real prefilter "
                "costs something to run"
            ),
        },
        "ceiling": {
            "definition": (
                "full_total / oracle_total; an upper bound on the speedup from "
                "replacing local candidate generation by any external generator "
                "whose rows are still certified locally"
            ),
            "value": full_seconds / oracle_seconds,
            "candidate_rows_avoided": (
                _stage_counts(full_generation)["screened_candidate_rows"]
                - oracle_screening.counts.nonstructural_candidate_rows
            ),
            "determinant_certifications_still_required": (
                oracle_screening.counts.certified_rows
            ),
        },
        "cross_path_agreement": {
            "retained_tau_sets_equal": full_taus == oracle_row_taus,
            "full_only": [list(tau) for tau in sorted(full_taus - oracle_row_taus)],
            "oracle_only": [list(tau) for tau in sorted(oracle_row_taus - full_taus)],
        },
    }


def build_report(systems: tuple[tuple[int, int], ...], *, workers: int) -> dict[str, Any]:
    """Measure every requested system and score the pre-registered predictions."""
    measurements = [measure_system(n, d, workers=workers) for n, d in systems]
    by_system = {entry["system"]: entry for entry in measurements}

    ledger_rows = [
        row
        for entry in measurements
        for row in entry["oracle_path"]["recertification_ledger"]
    ]
    all_certified = all(row["local_status"] == "certified" for row in ledger_rows)
    composed = [entry for entry in measurements if entry["composition_supported"]]
    all_regressions_pass = bool(composed) and all(
        entry["oracle_path"]["known_system_regression"].get("passed") is True
        for entry in composed
    )
    no_removals = bool(composed) and all(
        entry["oracle_path"]["counts"]["redundant_rows"] == 0 for entry in composed
    )
    rank8 = by_system.get("(3,8)")
    if rank8 is not None and not rank8["composition_supported"]:
        rank8 = None

    predictions = [
        {
            "id": "P1",
            "statement": (
                "every published nonstructural row screens locally as certified"
            ),
            "verdict": "PASS" if all_certified else "FAIL",
            "scope": len(ledger_rows),
            "independent_scope": len(measurements),
            "measured": {
                "certified": sum(
                    1 for row in ledger_rows if row["local_status"] == "certified"
                ),
                "other_statuses": sorted(
                    {
                        row["local_status"]
                        for row in ledger_rows
                        if row["local_status"] != "certified"
                    }
                ),
            },
        },
        {
            "id": "P2",
            "statement": (
                "the oracle path reproduces the published system exactly at every rank"
            ),
            "verdict": "PASS" if all_regressions_pass else "FAIL",
            "scope": len(composed),
            "independent_scope": len(composed),
            "measured": {
                entry["system"]: entry["oracle_path"]["known_system_regression"]
                for entry in composed
            },
        },
        {
            "id": "P3",
            "statement": (
                "exact ordered-slice reduction removes zero rows from each "
                "published system supplied alone"
            ),
            "verdict": "PASS" if no_removals else "FAIL",
            "scope": len(composed),
            "independent_scope": len(composed),
            "measured": {
                entry["system"]: entry["oracle_path"]["counts"]["redundant_rows"]
                for entry in composed
            },
        },
        {
            "id": "P4",
            "statement": "oracle total for (3,8) is under five seconds",
            "verdict": (
                "UNINFORMATIVE"
                if rank8 is None
                else ("PASS" if rank8["oracle_path"]["seconds"]["total"] < 5.0 else "FAIL")
            ),
            "scope": 1,
            "independent_scope": 1,
            "measured": (
                None if rank8 is None else rank8["oracle_path"]["seconds"]["total"]
            ),
        },
        {
            "id": "P5",
            "statement": (
                "enumeration takes more than half the wall clock of the full "
                "(3,8) construction"
            ),
            "verdict": (
                "UNINFORMATIVE"
                if rank8 is None
                else (
                    "PASS"
                    if rank8["full_path"]["enumeration_share_of_total"] > 0.5
                    else "FAIL"
                )
            ),
            "scope": 1,
            "independent_scope": 1,
            "measured": (
                None if rank8 is None else rank8["full_path"]["enumeration_share_of_total"]
            ),
        },
        {
            "id": "P6",
            "statement": "the (3,8) oracle ceiling is at least five",
            "verdict": (
                "UNINFORMATIVE"
                if rank8 is None
                else ("PASS" if rank8["ceiling"]["value"] >= 5.0 else "FAIL")
            ),
            "scope": 1,
            "independent_scope": 1,
            "measured": None if rank8 is None else rank8["ceiling"]["value"],
        },
    ]

    return {
        "schema_version": SCHEMA_VERSION,
        "preregistration": "docs/prereg_bdr_constructor_comparison.md",
        "generator": "scripts/bdr_constructor_comparison.py",
        "external_pin": EXTERNAL_PIN,
        "machine": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "measurement_caveat": (
            "Wall clock on one machine, single worker, no equally optimized "
            "control for the external path because the external path could not "
            "be run. Timings are not averaged across machines."
        ),
        "systems": measurements,
        "predictions": predictions,
        "benchmark_completeness": {
            "external_backend_executed": False,
            "label": "INCOMPLETE",
            "what_prevented_it": EXTERNAL_PIN["why_not_executed"],
            "what_is_established_anyway": (
                "an upper bound on the speedup any external candidate generator "
                "can deliver under the local certification requirement, plus a "
                "local recertification of every published row in the overlapping "
                "systems"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "results" / "data" / "bdr_constructor_comparison.json",
        help="artifact path",
    )
    parser.add_argument(
        "--systems",
        default="3:6,3:7,3:8",
        help="comma separated n:d systems to measure",
    )
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    systems = tuple(
        (int(part.split(":")[0]), int(part.split(":")[1]))
        for part in args.systems.split(",")
        if part
    )
    report = build_report(systems, workers=args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    for prediction in report["predictions"]:
        print(f"{prediction['id']}: {prediction['verdict']}", file=sys.stderr)
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
