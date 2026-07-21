#!/usr/bin/env python3
"""Reproduction script for the two-lead arithmetic analysis of the census.

Lead A (lattice geometry): exact incidence / simpliciality and the edge-cone and
    normal-cone lattice invariants (Smith normal form, indices, facet-ray
    pairings) at every simple vertex, across all nine determinate systems. This
    is the source of the den = normal-cone-index law (Conjecture 1) and the
    isotypic-purity count.
Lead B (phase Galois structure): the gauge-invariant loop holonomies of every
    certified interference vertex, their minimal polynomials, and kernel-dim
    distribution. Galois groups for the degree-8 polynomials need PARI polgalois
    (pari-galdata); the polynomials themselves are printed.

Run from the repo root: uv run scripts/leads_analysis.py
Deps: sympy (a project runtime dependency). Optional: pari-gp for polgalois.

Note: full_hrep includes the equality constraints, so systems with equalities
(only the Borland-Dennis system (3,6), three equalities) are handled and the
simple-vertex census is complete (183, not the 179 obtained when (3,6) is
dropped). The affine hull is {sum = N} together with those equalities, so the
polytope dimension is d - rank(ones + equalities), not d - 1 in general.
"""
import json
import pathlib
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gpc_census.constraints import constraints  # noqa: E402

import sympy as sp  # noqa: E402
from sympy import I, Matrix, Symbol  # noqa: E402
from sympy.matrices.normalforms import smith_normal_form  # noqa: E402

SYSTEMS = [(3, 6), (3, 7), (3, 8), (4, 8), (3, 9), (4, 9), (3, 10), (4, 10), (5, 10)]


def full_hrep(N, d):
    """Inequalities (GPCs + ordering walls + nonnegativity) and equalities."""
    s = constraints(N, d)
    A = [[F(c) for c in q["coeffs"]] for q in s["inequalities"]]
    b = [F(q["rhs"]) for q in s["inequalities"]]
    labels = [("GPC", i) for i in range(len(A))]
    for i in range(d - 1):
        row = [F(0)] * d
        row[i], row[i + 1] = F(-1), F(1)
        A.append(row)
        b.append(F(0))
        labels.append(("ORD", i))
    row = [F(0)] * d
    row[d - 1] = F(-1)
    A.append(row)
    b.append(F(0))
    labels.append(("NN", d - 1))
    E = [[F(c) for c in q["coeffs"]] for q in s.get("equalities", [])]
    return A, b, labels, E


def rank_q(rows):
    M = [r[:] for r in rows]
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


def prim(iv):
    g = 0
    for x in iv:
        g = gcd(g, abs(x))
    return [x // g for x in iv] if g else iv


def _verdicts(N, d):
    out = {}
    for line in (ROOT / "results" / "data" / "states.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        if e["system"] == f"({N},{d})":
            out[tuple(e["integer_form"])] = e["classified"]
    return out


def lead1(N, d):
    """Return {dim, simple:[per-vertex lattice record]} for one system."""
    A, b, _labels, E = full_hrep(N, d)
    V = json.loads(
        (ROOT / "results" / "data" / "vertices" / f"vertices_{N}_{d}.json").read_text())
    verts = [[F(x) for x in v["spectrum"]] for v in V]
    iforms = [tuple(v["integer_form"]) for v in V]
    dens = [v["denominator"] for v in V]
    ones = [F(1)] * d
    aff = [ones] + E                       # affine hull: sum = N plus equalities
    dim = d - rank_q(aff)
    active = []
    for v in verts:
        act = [j for j in range(len(A)) if sum(a * x for a, x in zip(A[j], v)) == b[j]]
        active.append(set(act))
        assert rank_q(aff + [A[j] for j in act]) == d
    n = len(verts)
    edges = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            common = active[i] & active[j]
            if len(common) < dim - 1:
                continue
            if rank_q(aff + [A[k] for k in common]) == d - 1:
                edges[i].append(j)
                edges[j].append(i)
    verdict = _verdicts(N, d)
    simple = []
    for i in range(n):
        if len(edges[i]) != dim:
            continue
        rays = []
        for jj in edges[i]:
            vec = [verts[jj][k] - verts[i][k] for k in range(d)]
            den = 1
            for x in vec:
                den = den * x.denominator // gcd(den, x.denominator)
            rays.append(prim([int(x * den) for x in vec]))
        R = Matrix([r[:dim] for r in rays])      # first dim coords: unimodular chart
        Nrm = [prim([int(A[j][k]) - int(A[j][d - 1]) for k in range(dim)])
               for j in active[i]]
        W = Matrix(Nrm)
        pair = sorted({abs(sum(int(A[j][k]) * r[k] for k in range(d)))
                       for j in active[i] for r in rays
                       if sum(int(A[j][k]) * r[k] for k in range(d)) != 0})
        simple.append({
            "system": f"({N},{d})", "index": i, "iform": iforms[i], "den": dens[i],
            "verdict": verdict.get(iforms[i]), "edge_index": int(abs(R.det())),
            "snf": [int(smith_normal_form(R)[k, k]) for k in range(dim)],
            "normal_index": int(abs(W.det())) if W.shape[0] == W.shape[1] else None,
            "pairings": pair})
    return {"dim": dim, "n": n, "simple": simple}


def lead1_all():
    tot = 0
    cls = Counter()
    den_eq = Counter()
    den_seen = Counter()
    pure = 0
    for (N, d) in SYSTEMS:
        r = lead1(N, d)
        for s in r["simple"]:
            tot += 1
            k = "interference" if s["verdict"] == "INTERFERENCE" else "design"
            cls[k] += 1
            if len(s["pairings"]) == 1:
                pure += 1
            if s["normal_index"] is not None:
                den_seen[k] += 1
                den_eq[k] += (s["den"] == s["normal_index"])
        print(f"({N},{d}): dim {r['dim']}, {len(r['simple'])} simple", flush=True)
    print(f"\nTOTAL simple = {tot}  (interference {cls['interference']}, "
          f"design {cls['design']})")
    print(f"den == normal-cone index: interference "
          f"{den_eq['interference']}/{den_seen['interference']}, design "
          f"{den_eq['design']}/{den_seen['design']} "
          f"(violate {den_seen['design'] - den_eq['design']})")
    print(f"isotypically pure vertices: {pure}")


def lead2():
    x = Symbol("x")
    polys = defaultdict(list)
    kd = Counter()
    for line in (ROOT / "results" / "data" / "states.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        e = json.loads(line)
        cf = e.get("closed_form")
        if not cf or e["classified"] != "INTERFERENCE":
            continue
        dets = cf["support_dets"]
        thetas = [sp.arg(sp.sympify(p)) for p in cf["pretty"]]
        nm = max(max(t) for t in dets) + 1
        M = Matrix([[1 if o in t else 0 for o in range(nm)] for t in dets])
        K = M.T.nullspace()
        kd[len(K)] += 1
        for v in K:
            v = v * sp.lcm([sp.fraction(sp.nsimplify(vi))[1] for vi in v])
            inv = sp.simplify(sum(vi * th for vi, th in zip(list(v), thetas)))
            c, s = sp.simplify(sp.cos(inv)), sp.simplify(sp.sin(inv))
            if s == 0:
                continue
            mp = sp.minimal_polynomial(sp.simplify(c + I * s), x)
            polys[sp.sstr(mp)].append((e["system"], e["index"]))
    n_hol = sum(len(w) for w in polys.values())
    print(f"kernel dims: {dict(kd)}  "
          f"(loop-free {kd[0]}, with holonomy {sum(v for k, v in kd.items() if k)})")
    print(f"nontrivial holonomies: {n_hol}; distinct minimal polynomials: {len(polys)}")
    for mp, wh in sorted(polys.items()):
        print(f"  {mp}  [{len(wh)}] e.g. {wh[0]}")
    print("\nFor degree-8 polynomials run in PARI/GP: polgalois(<poly>)")


def main():
    lead1_all()
    print()
    lead2()


if __name__ == "__main__":
    main()
