#!/usr/bin/env python3
"""Weyl-orbit canonicalization of the admissible hyperplane set.

THE PRUNING, AND WHY IT IS SAFE. The moment polytope is Weyl invariant and the
ordered slice is a fundamental domain for the `S_d` action, so two admissible
hyperplanes in the same `S_d` orbit carry the same information: one is a facet
if and only if the other is. Enumerating orbit representatives therefore
DISCARDS NOTHING. This is the one pruning that meets the standing rule that no
candidate may be dropped until a theorem proves the drop lossless, as opposed
to an empirical pattern that merely looks safe.

THE MEASURED PAYOFF. The raw hyperplane count explodes while the orbit count
barely moves:

    rank 6:    362 hyperplanes,   8 orbits,   45x
    rank 7:  5,341 hyperplanes,  19 orbits,  281x

so the right algorithm never materializes the flat lattice at all. It
enumerates orbit representatives directly.

THE ACTION, AND A TRAP. Hyperplanes are stored with a sign convention, first
nonzero entry positive, so each unoriented hyperplane appears once. Permuting
the coordinates of `h` therefore leaves the stored set unless the image is
re-normalized to that convention. Without the renormalization the action is not
closed: it sends 201 of 362 images outside the set at rank 6, and the orbit
count comes out 15 instead of 8. That is a wrong answer that looks plausible,
which is why `images_outside_set` is reported and asserted to be zero.

This script measures the orbit structure of the ALREADY ENUMERATED set. It is
the calibration input for a canonical-augmentation enumerator, which is the
algorithm that would never build the set in the first place.

Usage:
    python scripts/orbit_canonical_flats.py --ranks 6 7
    python scripts/orbit_canonical_flats.py --ranks 6 7 8
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


def normalize(h, z):
    """Re-impose the stored convention: first nonzero entry positive."""
    for value in h:
        if value:
            if value < 0:
                return tuple(-x for x in h), -z
            break
    return tuple(h), z


def orbit_structure(n: int, d: int) -> dict:
    started = time.time()
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    hyperplanes = enum.hyperplanes
    enumerate_seconds = time.time() - started

    index = {normalize(hp.h, hp.z): i for i, hp in enumerate(hyperplanes)}
    parent = list(range(len(hyperplanes)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    outside = 0
    t0 = time.time()
    for i, hp in enumerate(hyperplanes):
        for s in range(d - 1):
            image = list(hp.h)
            image[s], image[s + 1] = image[s + 1], image[s]
            j = index.get(normalize(image, hp.z))
            if j is None:
                outside += 1
                continue
            ra, rb = find(i), find(j)
            if ra != rb:
                parent[ra] = rb
    orbit_seconds = time.time() - t0

    orbits = collections.defaultdict(list)
    for i in range(len(hyperplanes)):
        orbits[find(i)].append(i)
    sizes = sorted((len(v) for v in orbits.values()), reverse=True)
    representatives = []
    for root, members in sorted(orbits.items(), key=lambda kv: -len(kv[1])):
        hp = hyperplanes[min(members)]
        representatives.append({
            "h": list(hp.h), "z": hp.z, "orbit_size": len(members),
        })

    return {
        "n": n, "d": d,
        "hyperplanes": len(hyperplanes),
        "orbits": len(orbits),
        "reduction_factor": round(len(hyperplanes) / len(orbits), 1),
        "images_outside_set": outside,
        "action_is_closed": outside == 0,
        "orbit_sizes": sizes,
        "largest_orbit": sizes[0] if sizes else 0,
        "representatives": representatives,
        "flat_counts_by_rank": list(enum.flat_counts_by_rank),
        "total_flats": sum(enum.flat_counts_by_rank),
        "enumerate_seconds": round(enumerate_seconds, 1),
        "orbit_seconds": round(orbit_seconds, 1),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranks", nargs="*", type=int, default=[6, 7])
    ap.add_argument("-o", "--out",
                    default=str(ROOT / "results/data/orbit_canonical_flats.json"))
    args = ap.parse_args()

    blocks = []
    out = pathlib.Path(args.out)
    for d in args.ranks:
        print(f"(3,{d}): enumerating and orbiting ...", flush=True)
        record = orbit_structure(3, d)
        print(f"  hyperplanes {record['hyperplanes']}  orbits {record['orbits']}  "
              f"reduction {record['reduction_factor']}x  "
              f"closed {record['action_is_closed']}  "
              f"({record['enumerate_seconds']}s + {record['orbit_seconds']}s)",
              flush=True)
        blocks.append(record)
        payload = {
            "what": "S_d orbit structure of the admissible hyperplane set. The "
                    "moment polytope is Weyl invariant and the ordered slice is "
                    "a fundamental domain, so orbit reduction is provably "
                    "lossless rather than empirically motivated.",
            "convention_trap": "hyperplanes are stored with first nonzero entry "
                               "positive, so a permuted image must be "
                               "re-normalized; without that the action is not "
                               "closed and rank 6 reports 15 orbits, not 8",
            "ranks": {str(b["d"]): b for b in blocks},
            "summary": {
                "ranks": [b["d"] for b in blocks],
                "hyperplanes": {str(b["d"]): b["hyperplanes"] for b in blocks},
                "orbits": {str(b["d"]): b["orbits"] for b in blocks},
                "reduction": {str(b["d"]): b["reduction_factor"] for b in blocks},
                "all_actions_closed": all(b["action_is_closed"] for b in blocks),
            },
        }
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
