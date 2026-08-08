#!/usr/bin/env python3
"""What is stable about the facet system, and can an intrinsic criterion
replace exact reduction.

Pre-registration: `docs/prereg_facet_index_stability.md`.
Narrative: `docs/facet_index_stability.md`.
Artifact: `results/data/facet_index_stability.json`.

Part A, LANDING and not scored: the padding and shape census. Trailing-zero
padding is tested as a candidate map on facets, and the primitive `tau` normals
are measured for entry size and shape count across every published system. The
tables were computed in scratch before the pre-registration existed, so they
are landed with a generator and a guard test rather than reported as
discoveries. Only P6, the extension to `N != 3`, is predictive.

Part B, PREDICTIVE: the coefficient-one facetness converse. Exact
cyclic-Schubert evaluation is deferred past redundancy elimination in the
production pipeline, so the coefficients of the certified rows that reduction
REMOVES have never been computed. This computes them for every certified row at
`(3,7)` and `(3,8)` and cross-tabulates against the retained and redundant
partition.

Resource handling is fail closed. A cyclic evaluation exceeding its per-row
budget is recorded as `resource_limit_unresolved` and never becomes a
coefficient claim or evidence for either direction of the converse.
"""

from __future__ import annotations

import argparse
import json
import platform
import signal
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from gpc_census.constraints import constraints  # noqa: E402
from gpc_census.dataset import vertices  # noqa: E402
from gpc_census.generation.bdr import tau_from_constraint  # noqa: E402
from gpc_census.generation.candidate_pipeline import screen_tau_candidates  # noqa: E402
from gpc_census.generation.cyclic_schubert import certify_cyclic_schubert  # noqa: E402
from gpc_census.generation.model import IntegralConstraint  # noqa: E402
from gpc_census.generation.orbit_native import (  # noqa: E402
    enumerate_orbit_survivor_taus,
)
from gpc_census.generation.redundancy import (  # noqa: E402
    certify_ordered_slice_redundancy,
)

SCHEMA_VERSION = 1

# Every published system, so Part A is a census rather than a sample.
PUBLISHED_SYSTEMS = (
    (3, 6), (3, 7), (3, 8), (3, 9), (3, 10),
    (4, 7), (4, 8), (4, 9), (4, 10),
    (5, 8), (5, 9), (5, 10),
)
# Ranks where the composition layer accepts a full certified-row partition.
COEFFICIENT_SYSTEMS = ((3, 7), (3, 8))

# Recorded so the landing run can be checked against what scratch found rather
# than silently replacing it. See Part A of the pre-registration.
SCRATCH_PADDING_VIOLATIONS = {"(3,6)->(3,7)": 1, "(3,7)->(3,8)": 0,
                              "(3,8)->(3,9)": 21, "(3,9)->(3,10)": 8}
SCRATCH_MAX_TAU_ENTRY = {"(3,6)": 1, "(3,7)": 2, "(3,8)": 7,
                         "(3,9)": 7, "(3,10)": 7}
SCRATCH_SHAPE_COUNTS = {"(3,6)": 1, "(3,7)": 1, "(3,8)": 9,
                        "(3,9)": 8, "(3,10)": 16}


class _RowBudgetExceeded(Exception):
    """Raised when one cyclic evaluation exceeds its wall-clock budget."""


def _alarm(_signum: int, _frame: object) -> None:
    raise _RowBudgetExceeded


def _published_constraints(n: int, d: int) -> tuple[IntegralConstraint, ...]:
    return tuple(
        IntegralConstraint(tuple(int(v) for v in row["coeffs"]), int(row["rhs"]))
        for row in constraints(n, d)["inequalities"]
    )


def _exact_vertices(n: int, d: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(str(value)) for value in vertex["spectrum"])
        for vertex in vertices(n, d)
    )


def shape_census(n: int, d: int) -> dict[str, Any]:
    """Measure normal-entry size and shape count for one published system."""
    rows = _published_constraints(n, d)
    taus = [tau_from_constraint(row, n) for row in rows]
    shapes = {tuple(sorted(tau)) for tau in taus}
    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "facet_count": len(rows),
        "max_abs_tau_entry": max(max(abs(v) for v in tau) for tau in taus),
        "max_abs_affine_coefficient": max(max(abs(v) for v in r.coeffs) for r in rows),
        "max_rhs": max(r.rhs for r in rows),
        "distinct_sorted_tau_shapes": len(shapes),
        "shape_to_facet_ratio": len(shapes) / len(rows),
        "lookup_source": constraints(n, d)["source"],
    }


def _tight_padding():
    """Return the ordered-stability branch's tight padding extension.

    Imported rather than reimplemented so the two campaigns cannot drift.
    """
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    from ordered_representation_stability import pi_pad

    return pi_pad


def padding_census(n: int, d: int) -> dict[str, Any]:
    """Test trailing-zero padding of ``(n,d)`` facets against ``(n,d+1)``.

    Three outcomes per padded row, decided exactly over the rationals against
    the certified vertex list: violated by some vertex, valid and tight on the
    polytope, or valid and strictly interior. A padded row can also be valid
    and tight without being a facet of the higher system, so facet membership
    is reported separately from tightness.
    """
    lower = _published_constraints(n, d)
    higher_taus = {tau_from_constraint(row, n) for row in _published_constraints(n, d + 1)}
    try:
        higher_vertices = _exact_vertices(n, d + 1)
    except KeyError:
        # Some dual systems carry constraints but no precomputed vertex file.
        # Facet membership is still decidable; exact validity is not, and is
        # reported as unavailable rather than guessed.
        return {
            "map": f"({n},{d})->({n},{d + 1})",
            "lower_facets": len(lower),
            "validity_available": False,
            "unavailable_reason": f"no precomputed vertices for ({n},{d + 1})",
            "present_in_higher_facet_list": sum(
                1
                for row in lower
                if tau_from_constraint(
                    IntegralConstraint(row.coeffs + (0,), row.rhs), n
                )
                in higher_taus
            ),
        }

    # Reconciliation with the ordered-stability branch, added at merge time.
    # That campaign found a TIGHT padding extension that is total and exact
    # (docs/ordered_representation_stability.md), which reads as a
    # contradiction of the failure measured here until one notices the maps
    # differ: this one appends zero in affine coordinates, that one appends the
    # largest last entry keeping the row valid. Its function is imported rather
    # than reimplemented, so there stays one generator per number.
    tight_extension = _tight_padding()
    lower_vertices_available = True
    agrees_with_tight = 0

    violated = tight = interior = in_facet_list = 0
    examples: list[dict[str, Any]] = []
    for row in lower:
        padded = IntegralConstraint(row.coeffs + (0,), row.rhs)
        values = [
            sum(c * v for c, v in zip(padded.coeffs, vertex, strict=True))
            for vertex in higher_vertices
        ]
        is_violated = any(value > padded.rhs for value in values)
        is_tight = any(value == padded.rhs for value in values)
        padded_tau = tau_from_constraint(padded, n)
        is_facet = padded_tau in higher_taus
        if is_violated:
            violated += 1
        elif is_tight:
            tight += 1
        else:
            interior += 1
        if is_facet:
            in_facet_list += 1
        if lower_vertices_available:
            row_tau = tau_from_constraint(row, n)
            tight_tau = tight_extension(row_tau, [list(v) for v in higher_vertices])
            if tight_tau is not None and tuple(tight_tau) == tuple(padded_tau):
                agrees_with_tight += 1
        if is_violated and len(examples) < 3:
            worst = max(values)
            examples.append(
                {
                    "lower_row": {"coeffs": list(row.coeffs), "rhs": row.rhs},
                    "padded_tau": list(padded_tau),
                    "max_value_over_higher_vertices": str(worst),
                    "rhs": padded.rhs,
                    "excess": str(worst - padded.rhs),
                }
            )
    return {
        "map": f"({n},{d})->({n},{d + 1})",
        "lower_facets": len(lower),
        "validity_available": True,
        "violated_at_higher_rank": violated,
        "valid_and_tight": tight,
        "valid_but_strictly_interior": interior,
        "present_in_higher_facet_list": in_facet_list,
        "padding_is_a_valid_map": violated == 0,
        "agrees_with_tight_extension": agrees_with_tight,
        "facet_membership_matches_tight_agreement": (
            in_facet_list == agrees_with_tight
        ),
        "violation_examples": examples,
    }


def _certified_partition(n: int, d: int, workers: int) -> dict[str, Any]:
    """Build the certified rows of one system and split them by reduction."""
    enumeration = enumerate_orbit_survivor_taus(n, d)
    if not enumeration.counts_are_conserved:
        raise AssertionError("orbit-native enumeration failed count conservation")
    screening = screen_tau_candidates(
        enumeration.survivor_taus,
        particle_number=n,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=workers,
    )
    if not screening.all_candidates_resolved:
        raise AssertionError("screening left unresolved rows; refusing to partition")
    certified = tuple(
        row for row in screening.candidates if row.status == "certified"
    )
    reduction = certify_ordered_slice_redundancy(
        tuple(row.constraint for row in certified), d
    )
    retained = set(reduction.retained_indices)
    return {
        "certified": certified,
        "retained_indices": retained,
        "redundant_indices": set(reduction.redundant_indices),
        "screening_counts": screening.counts.as_dict(),
    }


def coefficient_census(n: int, d: int, *, workers: int, budget: float) -> dict[str, Any]:
    """Evaluate the cyclic coefficient of EVERY certified row, not only facets.

    The production pipeline defers this past reduction, so the removed rows
    have never been evaluated. Their coefficients are the converse direction of
    the AK2008 sharpness criterion.
    """
    partition = _certified_partition(n, d, workers)
    certified = partition["certified"]
    retained = partition["retained_indices"]

    rows: list[dict[str, Any]] = []
    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for index, row in enumerate(certified):
            is_facet = index in retained
            started = time.perf_counter()
            signal.setitimer(signal.ITIMER_REAL, budget)
            try:
                verdict = certify_cyclic_schubert(
                    row.constraint.coeffs,
                    row.constraint.rhs,
                    particle_number=n,
                    include_modular=False,
                )
                status = verdict.status
                coefficient = verdict.coefficient
            except _RowBudgetExceeded:
                status = "resource_limit_unresolved"
                coefficient = None
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0.0)
            rows.append(
                {
                    "certified_index": index,
                    "tau": list(row.tau),
                    "constraint": row.constraint.as_dict(),
                    "is_facet": is_facet,
                    "cyclic_status": status,
                    "coefficient": coefficient,
                    "coefficient_one": coefficient == 1,
                    "seconds": time.perf_counter() - started,
                }
            )
    finally:
        signal.signal(signal.SIGALRM, previous)

    facets = [row for row in rows if row["is_facet"]]
    removed = [row for row in rows if not row["is_facet"]]
    unresolved = [
        row for row in rows if row["cyclic_status"] == "resource_limit_unresolved"
    ]

    def _tally(subset: list[dict[str, Any]]) -> dict[str, Any]:
        return {
            "rows": len(subset),
            "coefficient_one": sum(1 for row in subset if row["coefficient_one"]),
            "coefficient_zero": sum(1 for row in subset if row["coefficient"] == 0),
            "coefficient_nonunit_nonzero": sum(
                1
                for row in subset
                if row["coefficient"] is not None
                and row["coefficient"] not in (0, 1)
            ),
            "unresolved": sum(
                1
                for row in subset
                if row["cyclic_status"] == "resource_limit_unresolved"
            ),
            "ill_formed": sum(
                1
                for row in subset
                if row["cyclic_status"]
                in {"degree_mismatch", "rhs_not_exterior_weight"}
            ),
            "distinct_coefficients": sorted(
                {row["coefficient"] for row in subset if row["coefficient"] is not None}
            ),
        }

    removed_with_one = [row for row in removed if row["coefficient_one"]]

    # Added AFTER the predictions were scored, and labelled as such. The
    # forward direction holding while the converse fails makes coefficient one
    # a SOUND PREFILTER rather than a criterion: it cannot discard a facet at a
    # rank where every facet has coefficient one, but it does not decide
    # facetness. This prices that prefilter on the reduction stage, which the
    # BDR comparison identified as the term that grows with facet count.
    all_rows = tuple(row.constraint for row in certified)
    kept = tuple(
        row.constraint
        for index, row in enumerate(certified)
        if rows[index]["coefficient_one"]
    )
    started = time.perf_counter()
    full_reduction = certify_ordered_slice_redundancy(all_rows, d)
    full_seconds = time.perf_counter() - started
    started = time.perf_counter()
    filtered_reduction = certify_ordered_slice_redundancy(kept, d)
    filtered_seconds = time.perf_counter() - started
    prefilter = {
        "added_after_scoring": True,
        "soundness_basis": (
            "every facet has coefficient one at this rank, measured not proved"
        ),
        "reduction_input_rows_before": len(all_rows),
        "reduction_input_rows_after": len(kept),
        "input_reduction_fraction": 1 - len(kept) / len(all_rows),
        "seconds_before": full_seconds,
        "seconds_after": filtered_seconds,
        "reduction_stage_speedup": full_seconds / filtered_seconds,
        "retained_sets_identical": (
            {tuple(row.coeffs) + (row.rhs,) for row in full_reduction.retained_rows}
            == {
                tuple(row.coeffs) + (row.rhs,)
                for row in filtered_reduction.retained_rows
            }
        ),
        "caveat": (
            "sound only where every facet has coefficient one, which is an "
            "observation at two ranks and not a theorem. Using it as a gate "
            "needs that forward direction proved."
        ),
    }

    return {
        "system": f"({n},{d})",
        "particle_number": n,
        "rank": d,
        "screening_counts": partition["screening_counts"],
        "certified_rows": len(rows),
        "facet_rows": len(facets),
        "removed_rows": len(removed),
        "facets": _tally(facets),
        "removed": _tally(removed),
        "criterion": {
            "statement": "coefficient one iff facet, on certified rows",
            "forward_holds_every_facet_has_one": (
                all(row["coefficient_one"] for row in facets) if facets else None
            ),
            "converse_holds_no_removed_row_has_one": len(removed_with_one) == 0,
            "removed_rows_with_coefficient_one": len(removed_with_one),
            "separates_exactly": (
                bool(facets)
                and all(row["coefficient_one"] for row in facets)
                and not removed_with_one
                and not unresolved
            ),
            "counterexamples": [
                {
                    "certified_index": row["certified_index"],
                    "tau": row["tau"],
                    "constraint": row["constraint"],
                    "coefficient": row["coefficient"],
                }
                for row in removed_with_one[:5]
            ],
        },
        "coefficient_one_prefilter": prefilter,
        "rows": rows,
    }


def build_report(*, workers: int, budget: float) -> dict[str, Any]:
    shapes = [shape_census(n, d) for n, d in PUBLISHED_SYSTEMS]
    by_system = {entry["system"]: entry for entry in shapes}

    paddings = []
    for n, d in PUBLISHED_SYSTEMS:
        if (n, d + 1) in PUBLISHED_SYSTEMS:
            paddings.append(padding_census(n, d))

    coefficients = [
        coefficient_census(n, d, workers=workers, budget=budget)
        for n, d in COEFFICIENT_SYSTEMS
    ]
    coeff_by_system = {entry["system"]: entry for entry in coefficients}

    rank8 = coeff_by_system.get("(3,8)")
    all_removed_with_one = sum(
        entry["criterion"]["removed_rows_with_coefficient_one"]
        for entry in coefficients
    )
    removed_tally_8 = None if rank8 is None else rank8["removed"]
    any_zero = any(
        entry["facets"]["coefficient_zero"] + entry["removed"]["coefficient_zero"] > 0
        for entry in coefficients
    )
    any_ill_formed = any(
        entry["facets"]["ill_formed"] + entry["removed"]["ill_formed"] > 0
        for entry in coefficients
    )
    non_n3 = [entry for entry in shapes if entry["particle_number"] != 3]
    p6_systems = [
        entry for entry in non_n3
        if entry["system"] in {"(4,8)", "(4,9)", "(4,10)", "(5,10)"}
    ]

    predictions = [
        {
            "id": "P1",
            "statement": "every retained row at (3,8) has exact coefficient 1",
            "verdict": (
                "UNINFORMATIVE"
                if rank8 is None
                else (
                    "PASS"
                    if rank8["criterion"]["forward_holds_every_facet_has_one"]
                    else "FAIL"
                )
            ),
            "scope": None if rank8 is None else rank8["facet_rows"],
            "independent_scope": 1,
            "measured": None if rank8 is None else rank8["facets"],
        },
        {
            "id": "P2",
            "statement": (
                "at least one redundant certified row at (3,7) or (3,8) has "
                "coefficient 1, refuting coefficient-one implies facet"
            ),
            "verdict": "PASS" if all_removed_with_one > 0 else "FAIL",
            "scope": sum(entry["removed_rows"] for entry in coefficients),
            "independent_scope": len(coefficients),
            "measured": {
                entry["system"]: entry["criterion"][
                    "removed_rows_with_coefficient_one"
                ]
                for entry in coefficients
            },
        },
        {
            "id": "P3",
            "statement": (
                "among redundant certified rows at (3,8), both coefficient-one "
                "and non-coefficient-one values occur"
            ),
            "verdict": (
                "UNINFORMATIVE"
                if removed_tally_8 is None
                else (
                    "PASS"
                    if removed_tally_8["coefficient_one"] > 0
                    and removed_tally_8["rows"] - removed_tally_8["coefficient_one"] > 0
                    else "FAIL"
                )
            ),
            "scope": None if removed_tally_8 is None else removed_tally_8["rows"],
            "independent_scope": 1,
            "measured": removed_tally_8,
        },
        {
            "id": "P4",
            "statement": "no certified row at either rank has exact coefficient 0",
            "verdict": "FAIL" if any_zero else "PASS",
            "scope": sum(entry["certified_rows"] for entry in coefficients),
            "independent_scope": len(coefficients),
            "measured": {
                entry["system"]: {
                    "facet_zeros": entry["facets"]["coefficient_zero"],
                    "removed_zeros": entry["removed"]["coefficient_zero"],
                }
                for entry in coefficients
            },
        },
        {
            "id": "P5",
            "statement": (
                "every certified row yields a well-formed cyclic problem, so no "
                "degree_mismatch and no rhs_not_exterior_weight"
            ),
            "verdict": "FAIL" if any_ill_formed else "PASS",
            "scope": sum(entry["certified_rows"] for entry in coefficients),
            "independent_scope": len(coefficients),
            "measured": {
                entry["system"]: {
                    "facet_ill_formed": entry["facets"]["ill_formed"],
                    "removed_ill_formed": entry["removed"]["ill_formed"],
                }
                for entry in coefficients
            },
        },
        {
            "id": "P6",
            "statement": (
                "for (4,8), (4,9), (4,10) and (5,10) the distinct sorted-tau "
                "shape count is at most one quarter of the facet count"
            ),
            "verdict": (
                "PASS"
                if p6_systems
                and all(entry["shape_to_facet_ratio"] <= 0.25 for entry in p6_systems)
                else "FAIL"
            ),
            "scope": len(p6_systems),
            "independent_scope": len(p6_systems),
            "measured": {
                entry["system"]: {
                    "facets": entry["facet_count"],
                    "shapes": entry["distinct_sorted_tau_shapes"],
                    "ratio": entry["shape_to_facet_ratio"],
                }
                for entry in p6_systems
            },
        },
    ]

    landing = {
        "verdict": "LANDED",
        "note": (
            "Part A tables were computed in scratch before the pre-registration "
            "existed, so under rule R4 they are landed with a generator and a "
            "guard test rather than scored as discoveries. Only P6 is "
            "predictive."
        ),
        "scratch_agreement": {
            "padding_violations": {
                entry["map"]: {
                    "scratch": SCRATCH_PADDING_VIOLATIONS.get(entry["map"]),
                    "landed": entry["violated_at_higher_rank"],
                    "agrees": SCRATCH_PADDING_VIOLATIONS.get(entry["map"])
                    == entry["violated_at_higher_rank"],
                }
                for entry in paddings
                if entry["map"] in SCRATCH_PADDING_VIOLATIONS
                and entry.get("validity_available")
            },
            "max_tau_entry": {
                system: {
                    "scratch": value,
                    "landed": by_system[system]["max_abs_tau_entry"],
                    "agrees": value == by_system[system]["max_abs_tau_entry"],
                }
                for system, value in SCRATCH_MAX_TAU_ENTRY.items()
            },
            "shape_counts": {
                system: {
                    "scratch": value,
                    "landed": by_system[system]["distinct_sorted_tau_shapes"],
                    "agrees": value == by_system[system]["distinct_sorted_tau_shapes"],
                }
                for system, value in SCRATCH_SHAPE_COUNTS.items()
            },
        },
    }

    return {
        "schema_version": SCHEMA_VERSION,
        "preregistration": "docs/prereg_facet_index_stability.md",
        "generator": "scripts/facet_index_stability.py",
        "machine": {
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "conditionality": (
            "The rank-9 and rank-10 constraint lists are conjecturally complete "
            "per AK2008 (results/data/PROVENANCE.md), so every statement resting "
            "on them is conditional on that completeness."
        ),
        "scope_caveat": (
            "The cyclic test is one selected coefficient of one recurrence, not "
            "the general Belkale-Kumar or Ressayre deformed product. A negative "
            "result refutes THAT criterion only."
        ),
        "part_a_shape_census": shapes,
        "part_a_padding_census": paddings,
        "part_a_landing": landing,
        "part_b_coefficient_census": coefficients,
        "predictions": predictions,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "results" / "data" / "facet_index_stability.json",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--per-row-budget",
        type=float,
        default=30.0,
        help="wall-clock seconds allowed for one cyclic evaluation",
    )
    args = parser.parse_args()

    report = build_report(workers=args.workers, budget=args.per_row_budget)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    for prediction in report["predictions"]:
        print(f"{prediction['id']}: {prediction['verdict']}", file=sys.stderr)
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
