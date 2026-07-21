#!/usr/bin/env python3
"""In-repo verification of the two facet laws, plus the depth triage score.

LAWS (see docs/RESEARCH.md, "Generator blueprint"):
  rhs LAW : every known GPC facet's rhs is an N-subset sum of its own sorted
            coefficient vector; its DEPTH (position among distinct subset
            sums, 1 = max = trivial) is small (<= 7 through rank 10).
  EDGE LAW: every facet's sorted coefficient vector is an EXTREMAL RAY of the
            arrangement cut on the ordered cone {a_1 >= ... >= a_d} by the
            subset-sum tie hyperplanes a_T = a_T', i.e. the tight conditions
            (ordering ties + subset-sum ties) have rank exactly d-2 -- a ray
            modulo the constant shift (1,...,1), under which the inequality
            system is invariant (Sum lambda = N).

TRIAGE SCORE (docs/RESEARCH.md, "Depth triage score"):
  total active depth of a vertex = sum of rhs depths over its saturated GPC
  facets. Calibration on certified (3,10) interference vertices:
  Spearman(depth, solve secs) = -0.625; median 8 s at depth > 12 vs 147 s at
  depth <= 12. Extremes predict the failure stratum of the open vertices
  (near-zero => search-bound, raise max_card/walltime; system-extreme high
  => ansatz-expressiveness-bound, route to the signed-cancellation
  extension). No universal cross-system threshold: compare within a system.

Usage:
  python scripts/facet_laws.py --verify            # both laws, all systems
  python scripts/facet_laws.py --triage 3 10       # depth score per vertex
"""
from __future__ import annotations
import argparse
import json
import sys
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from gpc_census.constraints import constraints, known_systems  # noqa: E402


def rank_q(rows):
    M = [[F(x) for x in r] for r in rows]
    m = len(M)
    n = len(M[0]) if m else 0
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, m) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
        if r == m:
            break
    return r


def facet_data(N, d):
    """Per inequality: (coeffs, rhs, depth, is_edge)."""
    subs = list(combinations(range(d), N))
    out = []
    for q in constraints(N, d)["inequalities"]:
        c = q["coeffs"]
        b = q["rhs"]
        a = sorted(c, reverse=True)
        sums = sorted({sum(a[i] for i in T) for T in subs}, reverse=True)
        in_sums = b in sums
        depth = (1 + sum(1 for v in sums if v > b)) if in_sums else None
        rows = []
        for i in range(d - 1):
            if a[i] == a[i + 1]:
                r = [0] * d
                r[i] = 1
                r[i + 1] = -1
                rows.append(r)
        bysum = {}
        for T in subs:
            bysum.setdefault(sum(a[i] for i in T), []).append(T)
        for g in bysum.values():
            T0 = g[0]
            for T in g[1:]:
                r = [0] * d
                for i in T0:
                    r[i] += 1
                for i in T:
                    r[i] -= 1
                if any(r):
                    rows.append(r)
        is_edge = (rank_q(rows) if rows else 0) == d - 2
        out.append((c, b, depth, is_edge))
    return out


def verify():
    tot = rhs_ok = edge_ok = 0
    print("system   #facets  rhs-law  edge-law  depth-distribution")
    for (N, d) in known_systems():
        data = facet_data(N, d)
        n = len(data)
        r = sum(1 for c, b, dep, e in data if dep is not None)
        e = sum(1 for c, b, dep, ed in data if ed)
        from collections import Counter
        dd = Counter(dep for c, b, dep, ed in data if dep is not None)
        print(f"({N},{d})   {n:4d}    {r}/{n}   {e}/{n}   "
              f"{dict(sorted(dd.items()))}")
        tot += n
        rhs_ok += r
        edge_ok += e
    print(f"\nTOTAL: rhs law {rhs_ok}/{tot}; edge law {edge_ok}/{tot}")
    return rhs_ok == tot and edge_ok == tot


def triage(N, d):
    data = facet_data(N, d)
    vfile = Path(__file__).resolve().parents[1] / "results" / "data" / \
        "vertices" / f"vertices_{N}_{d}.json"
    V = json.load(open(vfile))
    print(f"total active depth per vertex, ({N},{d})  "
          "(within-system comparison only)")
    scores = []
    for idx, v in enumerate(V):
        spec = [F(x) for x in v["spectrum"]]
        td = na = 0
        for c, b, dep, _e in data:
            if sum(F(ci) * x for ci, x in zip(c, spec)) == F(b):
                na += 1
                td += dep or 0
        scores.append((td, na, idx))
    for td, na, idx in sorted(scores, reverse=True):
        print(f"  v{idx:<4d} active {na:2d}  total-depth {td}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--triage", nargs=2, type=int, metavar=("N", "D"))
    args = ap.parse_args()
    if args.verify:
        raise SystemExit(0 if verify() else 1)
    if args.triage:
        triage(*args.triage)
        return
    ap.print_help()


if __name__ == "__main__":
    main()
