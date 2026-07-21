"""Exact combinatorial invariants of every determinate moment polytope, ranks
6-10, from the shipped V-representation (results/data/vertices) and the H-rep
(gpc_census.constraints). All rational, no convex-hull dependency.

Per system: affine dimension, vertex count f_0, facet count f_{dim-1}, the count
of simple vertices (on exactly dim facets), and the F = V + 1 relation the paper
flags at rank 8 (checked across all nine systems). Facets are identified from
the candidate H-rep (ordering walls lambda_i >= lambda_{i+1}, positivity
lambda_d >= 0, and the GPC inequalities), each accepted as a facet iff the
vertices saturating it span an affine subspace of dimension dim - 1.

Run: uv run scripts/polytope_invariants.py
Writes docs/polytope_invariants.json.
"""
from __future__ import annotations

import json
import pathlib
from fractions import Fraction as F

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SYSTEMS = [(3, 6), (3, 7), (3, 8), (4, 8), (3, 9), (4, 9), (3, 10), (4, 10), (5, 10)]


def affine_dim(points):
    """Exact affine dimension of a set of rational points (rows)."""
    if len(points) <= 1:
        return 0
    base = points[0]
    M = sp.Matrix([[a - b for a, b in zip(p, base)] for p in points[1:]])
    return M.rank()


def candidate_facets(n, d):
    """(coeffs, rhs) rows for a<=b half-spaces: ordering walls, positivity, GPCs."""
    from gpc_census.constraints import constraints
    rows = []
    for i in range(d - 1):                      # lambda_i - lambda_{i+1} >= 0
        c = [0] * d
        c[i], c[i + 1] = -1, 1
        rows.append((c, 0, f"order l{i+1}>=l{i+2}"))
    c = [0] * d                                 # lambda_d >= 0
    c[d - 1] = -1
    rows.append((c, 0, "pos l{}>=0".format(d)))
    sys_ = constraints(n, d)
    for k, q in enumerate(sys_["inequalities"]):
        rows.append(([F(x) for x in q["coeffs"]], F(q["rhs"]), f"gpc{k}"))
    return rows


def main() -> int:
    out = {"note": "Exact combinatorial invariants of the determinate moment "
                   "polytopes, ranks 6-10, from the shipped V- and H-reps.",
           "systems": {}}
    print(f"{'system':8}{'dim':>4}{'V':>5}{'F':>5}{'simple':>8}{'F=V+1':>7}")
    for n, d in SYSTEMS:
        verts = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
        pts = [[F(x) for x in v["spectrum"]] for v in verts]
        dim = affine_dim(pts)
        cands = candidate_facets(n, d)
        # incidence: which candidate rows each vertex saturates (equality holds)
        seen = {}                               # frozenset(vertices on facet) -> name
        for (c, rhs, name) in cands:
            vals = [sum(cc * pp for cc, pp in zip(c, p)) for p in pts]
            if any(val > rhs for val in vals):   # not a valid supporting half-space
                continue
            on = frozenset(i for i, val in enumerate(vals) if val == rhs)
            if not on:
                continue
            if affine_dim([pts[i] for i in on]) == dim - 1:   # facet-defining
                seen.setdefault(on, name)        # dedup coincident inequalities
        sat = list(seen.keys())                 # each facet once, by its vertex support
        nf = len(sat)
        per_vertex = [sum(1 for s in sat if i in s) for i in range(len(pts))]
        simple = sum(1 for k in per_vertex if k == dim)
        out["systems"][f"({n},{d})"] = {
            "dim": int(dim), "vertices": len(pts), "facets": nf,
            "simple_vertices": simple, "F_eq_V_plus_1": nf == len(pts) + 1,
            "min_facets_per_vertex": min(per_vertex),
            "max_facets_per_vertex": max(per_vertex)}
        print(f"({n},{d})".ljust(8) + f"{dim:>4}{len(pts):>5}{nf:>5}"
              f"{simple:>8}{str(nf==len(pts)+1):>7}")
    (ROOT / "docs" / "polytope_invariants.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print("wrote docs/polytope_invariants.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
