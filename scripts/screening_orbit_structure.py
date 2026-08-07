#!/usr/bin/env python3
"""Is Ressayre screening constant on Weyl orbits? NO, and here is the structure.

THE HYPOTHESIS THAT FAILED. `docs/orbit_canonical_flats.md` reduced candidate
ENUMERATION from hyperplanes to `S_d` orbits, losslessly, because the moment
polytope is Weyl invariant. The natural next step was to screen one
representative per orbit instead of every row. That is FALSE. At rank 7, 25 of
the 35 oriented-tau orbits carry mixed verdicts, several containing both
`certified` and `trace_rejected`.

WHY IT FAILS, EXACTLY. The Ressayre trace test compares two counts:

    (a) negative roots (i<j) with h_j - h_i < 0
    (b) exterior weights strictly below the hyperplane

Count (b) is `S_d` INVARIANT: permuting modes permutes the weights, and "below"
is equivariant, so the count is a property of the orbit. Count (a) is exactly
the INVERSION NUMBER of the sequence `h`, which is not invariant at all: it
ranges over the whole Mahonian spread as `h` is permuted. The dominant chamber
`lambda_0 >= ... >= lambda_{d-1}` is a fundamental domain and therefore is NOT
Weyl stable, so a chamber-referenced condition cannot be an orbit invariant.
The Weyl invariance of the moment cone is a statement about the ambient cone,
not about the chamber-facing inequality list.

WHAT SURVIVES, AND IT IS USEFUL. The trace test decomposes into an orbit
invariant and a permutation statistic:

    trace test passes  <=>  inv(sigma h) == B,     B = #weights below, invariant.

So `B` is computed ONCE per orbit and every row costs only an inversion count,
never a sweep over the `C(d,3)` weights. Better, the surviving permutations are
those of a multiset with a prescribed inversion number, whose count is a
Gaussian multinomial and which can be generated directly rather than filtered
out of the full orbit.

Usage:
    python scripts/screening_orbit_structure.py --ranks 6 7 8
"""
from __future__ import annotations

import argparse
import collections
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.generation import admissible_flats as af  # noqa: E402
from gpc_census.generation.exterior import exterior_weights  # noqa: E402
from gpc_census.generation.reference_generator import oriented_taus  # noqa: E402
from gpc_census.generation.ressayre import (  # noqa: E402
    _negative_root_set,
    admissible_candidate_from_tau,
)


def inversions(h) -> int:
    return sum(1 for i in range(len(h)) for j in range(i + 1, len(h))
               if h[j] < h[i])


def weights_below(h, z, occupied) -> int:
    return sum(1 for members in occupied if sum(h[i] for i in members) < z)


def analyse(n: int, d: int) -> dict:
    started = time.time()
    weights = exterior_weights(n, d)
    occupied = [tuple(i for i, bit in enumerate(w) if bit) for w in weights]
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    taus = oriented_taus(enum)

    # group rows by orbit, which for a tau is just its sorted multiset
    by_orbit: dict[tuple, list] = collections.defaultdict(list)
    for tau in taus:
        by_orbit[tuple(sorted(tau))].append(tau)

    invariant_checks = 0
    invariant_failures = 0
    rows = survivors = 0
    orbit_records = []
    for key, members in by_orbit.items():
        candidates = [admissible_candidate_from_tau(n, tau) for tau in members]
        # B is computed ONCE, then asserted invariant across the orbit
        reference = weights_below(candidates[0].h, candidates[0].z, occupied)
        for candidate in candidates:
            invariant_checks += 1
            if weights_below(candidate.h, candidate.z, occupied) != reference:
                invariant_failures += 1
        counts = [len(_negative_root_set(c.h)) for c in candidates]
        passed = sum(1 for c in counts if c == reference)
        rows += len(members)
        survivors += passed
        orbit_records.append({
            "sorted_tau": list(key),
            "orbit_rows": len(members),
            "weights_below_invariant": reference,
            "inversion_counts": sorted(set(counts)),
            "trace_survivors": passed,
        })

    return {
        "n": n, "d": d,
        "oriented_rows": rows,
        "orbits": len(by_orbit),
        "trace_survivors": survivors,
        "survivor_fraction": round(survivors / rows, 6),
        "orbits_with_a_survivor": sum(1 for r in orbit_records
                                      if r["trace_survivors"]),
        "orbits_with_mixed_inversion_counts": sum(
            1 for r in orbit_records if len(r["inversion_counts"]) > 1),
        "weights_below_checks": invariant_checks,
        "weights_below_failures": invariant_failures,
        "weights_below_is_orbit_invariant": invariant_failures == 0,
        "seconds": round(time.time() - started, 1),
        "orbit_records": orbit_records,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranks", nargs="*", type=int, default=[6, 7])
    ap.add_argument("-o", "--out",
                    default=str(ROOT / "results/data/screening_orbit_structure.json"))
    args = ap.parse_args()

    blocks = []
    out = pathlib.Path(args.out)
    for d in args.ranks:
        print(f"(3,{d}): ...", flush=True)
        record = analyse(3, d)
        print(f"  rows {record['oriented_rows']}  orbits {record['orbits']}  "
              f"trace survivors {record['trace_survivors']} "
              f"({100 * record['survivor_fraction']:.2f}%)  "
              f"B invariant {record['weights_below_is_orbit_invariant']}  "
              f"({record['seconds']}s)", flush=True)
        blocks.append(record)
        payload = {
            "what": "Ressayre screening is NOT constant on Weyl orbits. The "
                    "trace test decomposes into an orbit invariant (weights "
                    "below) and a permutation statistic (the inversion number "
                    "of h), and only the first is orbit invariant.",
            "refuted": "screening one representative per orbit and transporting "
                       "the verdict; the dominant chamber is a fundamental "
                       "domain and therefore not Weyl stable, so a "
                       "chamber-referenced condition cannot be an orbit "
                       "invariant",
            "survives": "trace test passes iff inv(sigma h) == B with B "
                        "invariant, so B costs one evaluation per ORBIT and "
                        "each row costs only an inversion count",
            "ranks": {str(b["d"]): b for b in blocks},
            "summary": {
                "ranks": [b["d"] for b in blocks],
                "oriented_rows": {str(b["d"]): b["oriented_rows"] for b in blocks},
                "orbits": {str(b["d"]): b["orbits"] for b in blocks},
                "trace_survivors": {str(b["d"]): b["trace_survivors"]
                                    for b in blocks},
                "survivor_fraction": {str(b["d"]): b["survivor_fraction"]
                                      for b in blocks},
                "weights_below_invariant_everywhere": all(
                    b["weights_below_is_orbit_invariant"] for b in blocks),
                "orbits_with_mixed_inversion_counts": {
                    str(b["d"]): b["orbits_with_mixed_inversion_counts"]
                    for b in blocks},
            },
        }
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
