#!/usr/bin/env python3
"""The external half of the constructor benchmark: BDR wall clock and stages.

WHY THIS EXISTS. `docs/bdr_constructor_comparison.md` was labelled INCOMPLETE
because no BDR runtime, candidate count or per-stage figure existed anywhere in
it. The stated reason was that the backend could not be fetched or run, which
was false and untested; see `docs/bdr_linear_triangular_gate.md`. This script
takes the measurements that were missing, so the label can be retired on
evidence rather than on a corrected excuse.

WHAT IS MEASURED. For each system, the upstream pipeline is run end to end and
every stage of `MomentConeStep.apply` is timed and counted: wall clock, CPU
time, and the size of the dataset entering and leaving each step.

ONLY THE PROBABILISTIC ARM RUNS HERE, and that is a measured fact rather than an
assumption. The upstream offers `symbolic` methods for dominancy and
birationality, which would be the like-for-like arm against this repository's
exact pipeline. They cannot run in this environment: `Is_Ram_contracted`
factorizes a univariate polynomial over the fraction field of `V.QV2`, which has
`2 * binom(d, N)` variables (70 at `(3,7)`, 112 at `(3,8)`), and libSingular in
this passagemath build raises `NotImplementedError` or crashes outright above
roughly 48 variables. `check_symbolic_arm_availability` below bisects that
threshold on a trivial polynomial `z^2 - 1`, so the limit is recorded as
evidence rather than asserted. `symbolic_int` is not a third option: both
`Is_Ram_contracted` and `is_not_contracted` reject it with a ValueError.

CONSEQUENCE, and it must travel with every number here. The measured arm is
BDR's default, which is what its authors ship and what a user would actually
run, so pricing it is legitimate. It is NOT a like-for-like control against
this repository, which is exact everywhere and ships a replayable certificate
per row. The comparison therefore favours BDR by exactly the cost of the
exactness it did not pay, and no conclusion may be drawn about which is faster
at equal rigour.

PROTOCOL, matching the local benchmark it will be compared with: one machine,
single worker (the upstream default executor is Sequential), repeated runs, no
averaging across machines. Reported per run so the spread is visible.

Run in a virtual environment that has `moment_cone` (see
`docs/bdr_linear_triangular_gate.md` for the recipe):

    python scripts/bdr_external_benchmark.py -o results/data/bdr_external_benchmark.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import statistics

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results/data/bdr_external_benchmark.json"

# The systems the local benchmark measures, plus (3,9) as a third point.
SYSTEMS = ((3, 7), (3, 8), (3, 9))
ARMS = ("probabilistic",)
REPEATS = 3
# Variable counts to bisect the libSingular factorization limit.
SYMBOLIC_PROBE_VARIABLES = (16, 32, 40, 48, 56, 70, 112)
SEED = 20260809

PINNED_REVISION = "7ab347d9a7837d68d92bde9fac74606913f395ca"


def check_symbolic_arm_availability() -> dict:
    """Bisect the libSingular limit that rules out the exact arm.

    The test polynomial is `z^2 - 1`, so a failure cannot be blamed on the size
    of the polynomials the pipeline actually builds. What varies is only the
    number of variables in the base ring.
    """
    import moment_cone  # noqa: F401  (initialises the sage stack)
    from sage.rings.polynomial.polynomial_ring_constructor import (  # type: ignore
        PolynomialRing,
    )
    from sage.rings.rational_field import QQ  # type: ignore

    probe = []
    for count in SYMBOLIC_PROBE_VARIABLES:
        base = PolynomialRing(QQ, [f"x{i}" for i in range(count)]).fraction_field()
        ring = PolynomialRing(base, "z")
        generator = ring.gen()
        try:
            (generator**2 - 1).factor()
        except Exception as error:  # noqa: BLE001  (any failure is a failure)
            probe.append(
                {"variables": count, "factors": False, "error": type(error).__name__}
            )
        else:
            probe.append({"variables": count, "factors": True, "error": None})

    working = [row["variables"] for row in probe if row["factors"]]
    failing = [row["variables"] for row in probe if not row["factors"]]
    return {
        "available": not failing,
        "test_polynomial": "z^2 - 1 over Frac(QQ[x_0..x_n])",
        "probe": probe,
        "largest_variable_count_that_factors": max(working) if working else None,
        "smallest_variable_count_that_fails": min(failing) if failing else None,
        "variables_needed": {
            "(3,7)": 70,
            "(3,8)": 112,
            "(3,9)": 168,
        },
        "why_it_matters": (
            "Is_Ram_contracted factorizes a univariate polynomial over the "
            "fraction field of V.QV2, which has 2*binom(d,N) variables. Every "
            "system measured here needs more than the limit, so the exact arm "
            "of the upstream birationality test cannot be run in this "
            "environment."
        ),
        "symbolic_int_is_not_an_alternative": (
            "Is_Ram_contracted and is_not_contracted both raise ValueError for "
            "symbolic_int; it is not accepted by these code paths."
        ),
    }


def _run_once(particles: int, orbitals: int, method: str) -> dict:
    """One end-to-end pipeline run, with every stage timed and counted."""
    from moment_cone import FermionRepresentation, LinearGroup  # type: ignore
    from moment_cone.main_steps import MomentConeStep  # type: ignore
    from moment_cone.task import Task  # type: ignore

    Task.all_tasks.clear()
    Task.quiet = True

    representation = FermionRepresentation(
        LinearGroup([orbitals]), particle_cnt=particles, seed=SEED
    )
    step = MomentConeStep(
        representation,
        quiet=True,
        store_steps=True,
        tpi_method=method,
        ram_schub_method=method,
        ram0_method=method,
    )
    dataset = step()

    stages = []
    for task in Task.all_tasks:
        if not Task.is_finished(task):
            continue
        wall_ns, cpu_ns = task.duration
        stages.append(
            {
                "stage": task.name,
                "wall_seconds": round(wall_ns / 1e9, 4),
                "cpu_seconds": round(cpu_ns / 1e9, 4),
            }
        )

    counts = []
    for inner in step.steps:
        output = getattr(inner, "output_dataset", None)
        if output is None:
            continue
        counts.append(
            {
                "stage": inner.name,
                "pending": len(list(output.pending())),
                "validated": len(list(output.validated())),
            }
        )

    total = next(
        (s["wall_seconds"] for s in stages if s["stage"] == "MomentConeStep"), None
    )
    return {
        "total_wall_seconds": total,
        "retained_rows": len(list(dataset.validated())),
        "stages": stages,
        "stage_counts": counts,
    }


def _measure(particles: int, orbitals: int, method: str) -> dict:
    runs = [_run_once(particles, orbitals, method) for _ in range(REPEATS)]
    totals = [run["total_wall_seconds"] for run in runs]
    retained = {run["retained_rows"] for run in runs}
    if len(retained) != 1:
        raise SystemExit(
            f"({particles},{orbitals}) {method}: unstable output across repeats, "
            f"retained {sorted(retained)}"
        )

    # The stage breakdown is reported from the median run, and the spread of
    # the totals is reported alongside so a single lucky run cannot be quoted.
    median_total = statistics.median(totals)
    representative = min(runs, key=lambda run: abs(run["total_wall_seconds"] - median_total))
    return {
        "system": f"({particles},{orbitals})",
        "particle_number": particles,
        "orbitals": orbitals,
        "method": method,
        "repeats": REPEATS,
        "retained_rows": retained.pop(),
        "total_wall_seconds_per_run": totals,
        "total_wall_seconds_median": round(median_total, 4),
        "total_wall_seconds_spread": round(max(totals) - min(totals), 4),
        "stages": representative["stages"],
        "stage_counts": representative["stage_counts"],
    }


def main() -> int:
    import platform

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--out", type=pathlib.Path, default=DEFAULT_OUT)
    arguments = parser.parse_args()

    symbolic_availability = check_symbolic_arm_availability()
    print(
        "exact arm available: "
        f"{symbolic_availability['available']}, factors up to "
        f"{symbolic_availability['largest_variable_count_that_factors']} "
        f"variables, fails from "
        f"{symbolic_availability['smallest_variable_count_that_fails']}",
        flush=True,
    )

    measurements: dict[str, dict] = {}
    for particles, orbitals in SYSTEMS:
        for method in ARMS:
            record = _measure(particles, orbitals, method)
            measurements.setdefault(record["system"], {})[method] = record
            print(
                f"{record['system']} {method}: "
                f"{record['total_wall_seconds_median']}s median, "
                f"{record['retained_rows']} retained, "
                f"spread {record['total_wall_seconds_spread']}s",
                flush=True,
            )

    payload = {
        "schema_version": 1,
        "generated_by": "scripts/bdr_external_benchmark.py",
        "status": "external_measurement",
        "not_preregistered": (
            "These measurements are post hoc. The pre-registration in "
            "docs/prereg_bdr_constructor_comparison.md scored P1 to P6 on the "
            "local half only, and nothing here is scored against it."
        ),
        "external_backend": {
            "name": "Bulois-Denis-Ressayre moment_cone",
            "mirror": "github.com/ea-icj/moment_cone",
            "pinned_revision": PINNED_REVISION,
            "license": "MIT",
            "executed_here": True,
        },
        "protocol": {
            "executor": "Sequential, one worker (the upstream default)",
            "repeats": REPEATS,
            "seed": SEED,
            "arms": list(ARMS),
            "exact_arm_unavailable": symbolic_availability,
            "caveat": (
                "Not an equally optimized control, and not a like-for-like one. "
                "The measured arm is BDR's probabilistic default while this "
                "repository is exact everywhere and ships a replayable "
                "certificate per row, so the comparison favours BDR by the cost "
                "of the exactness it did not pay. The exact arm could not be "
                "run here; see exact_arm_unavailable. Separately the stacks "
                "differ (CPython 3.11 under passagemath against CPython 3.14 "
                "with numpy and scipy) and the algorithms differ, so wall clock "
                "is an end-to-end cost on one machine, not a measure of "
                "algorithmic quality."
            ),
        },
        "machine": {
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "systems": measurements,
    }

    hashable = json.loads(json.dumps(payload))
    hashable.pop("machine", None)
    for system in hashable["systems"].values():
        for arm in system.values():
            for field in (
                "total_wall_seconds_per_run",
                "total_wall_seconds_median",
                "total_wall_seconds_spread",
                "stages",
            ):
                arm.pop(field, None)
    payload["canonical_sha256_without_this_field_or_timings"] = hashlib.sha256(
        json.dumps(hashable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
