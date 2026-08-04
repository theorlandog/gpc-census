#!/usr/bin/env python3
"""Orbit-closure (entanglement) moment polytopes for Lambda^3 C^d, d = 6 and 7.

Emits results/data/lambda3_c7_polytopes.json: one record per SLOCC class with
its certified inner hull, its certified outer H-polytope, the vertices of the
outer polytope, and an EXACT / BRACKET verdict saying whether the two met.

The d = 6 block is the validation gate: its four classes are the published
fermionic entanglement polytopes (Walter, Doran, Gross, Christandl, Science
340, 1205 (2013)), and the open orbit there must reproduce the ambient (3,6)
polytope exactly. The d = 7 block is the new computation.

Usage:
    python scripts/lambda3_c7_polytopes.py [-o results/data/lambda3_c7_polytopes.json]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from gpc_census import orbit  # noqa: E402


def build(d: int) -> dict:
    t0 = time.time()
    result = orbit.census(d)
    order = orbit.polytope_order(result)
    same = orbit.coincidences(result)
    exact = [k for k, r in result.items() if r["verdict"] == "EXACT"]
    return {
        "d": d,
        "n_classes": len(result),
        "n_exact": len(exact),
        "seconds": round(time.time() - t0, 1),
        "onehop_free_supports": len(orbit.onehop_free_supports(d)),
        "canonical_supports_by_class": {
            str(k): len(v) for k, v in orbit.canonical_supports_by_class(d).items()},
        "degeneration_order": {str(k): list(v)
                               for k, v in orbit.degeneration_order(d).items()},
        "polytope_containment": order,
        "polytope_coincidences": same,
        "classes": result,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", default="results/data/lambda3_c7_polytopes.json")
    ap.add_argument("--d", type=int, nargs="*", default=[6, 7])
    args = ap.parse_args()

    blocks = {}
    for d in args.d:
        print(f"computing Lambda^3 C^{d} ...", flush=True)
        blocks[str(d)] = build(d)
        b = blocks[str(d)]
        print(f"  {b['n_exact']}/{b['n_classes']} classes EXACT "
              f"({b['seconds']}s)", flush=True)
        for key, rec in b["classes"].items():
            flag = "" if rec["verdict"] == "EXACT" else (
                f"  ({len(rec['outer_vertices_not_reached'])} outer vertices "
                "not reached by the inner hull)")
            print(f"    {rec['label']:5s} orbit dim {key:>3s} rank {rec['rank']}  "
                  f"{len(rec['vertices']):2d} vertices  {rec['verdict']}{flag}")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "description": "Orbit-closure (entanglement) moment polytopes of "
                       "Lambda^3 C^d, one per SLOCC class",
        "method": "certified inner hull (torus orbits of one-hop-free supports, "
                  "plus the polytopes of certified degenerations) sandwiched "
                  "against a certified outer H-polytope (Schur-Horn/Newton "
                  "inequalities plus the ambient GPC system); verdict EXACT "
                  "when the two coincide, checked by exact rational LP",
        "blocks": blocks,
    }, indent=1) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
