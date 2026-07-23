#!/usr/bin/env python3
"""General reality solver for kdim-1 interference vertices, via the Cayley-Menger
(positivity-boundary) characterization of reality.

Reality is the boundary of a positivity region (docs/RESEARCH.md, "REALITY IS A
POSITIVITY BOUNDARY"): fixing the occupations pins each interference off-diagonal
to a Schur-Horn magnitude, the one-hop amplitudes form a closed polygon, and a REAL
extremal state is where every such polygon DEGENERATES -- each edge's Cayley-Menger
determinant vanishes. For a two-term edge the Cayley-Menger determinant is the wall
polynomial (tau*den^2 - p1 - p2)^2 - 4 p1 p2; a real state is a COMMON real zero of
all edges' determinants inside the positivity range. This generalizes
scripts/wall_test.py (single two-term block, one wall, one sign) to arbitrarily
many edges and both walls, closing the multi-class sign-compatibility gap by
letting verify_exact adjudicate the global sign assignment.

Candidate t values are the real roots of every edge's reality polynomial in the
positivity interval (the true simultaneous-degeneration point is a root of each);
verify_exact is the exact gate on the real-sign state, so no edge condition is
trusted implicitly. Diagnostic; scripts/ is ruff-excluded.

Run: .venv/bin/python scripts/wall_solver.py
"""
from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load():
    recs = {}
    for line in (ROOT / "results" / "data" / "states.jsonl").open():
        r = json.loads(line)
        recs[(r["system"], r["index"])] = r
    return recs


def solve(system, index, recs, max_signs=1 << 12):
    """Certify a real extremal state for one kdim-1 interference vertex, or not."""
    import sympy as sp
    from gpc_census.exactify import _active_offdiagonals, verify_exact

    r = recs[(system, index)]
    n, d = (int(x) for x in system.strip("()").split(","))
    cf = r.get("closed_form")
    if not cf or not cf.get("support_dets"):
        return "no state", None
    dets = [tuple(x) for x in cf["support_dets"]]
    w0 = [int(x) for x in cf["weights"]]
    den = int(cf["den"])
    spec = [Fraction(x, r["denominator"]) for x in r["integer_form"]]
    m = len(dets)
    M = sp.Matrix([[1 if o in T else 0 for o in range(d)] for T in dets])
    ns = M.T.nullspace()
    if len(ns) != 1:
        return f"kdim={len(ns)}", None
    u = ns[0]
    u = u * sp.lcm([c.q for c in u])
    u = [int(x) for x in u]
    t = sp.symbols("t")
    w = [sp.Integer(w0[i]) + t * u[i] for i in range(m)]
    occ, edges = _active_offdiagonals(d, dets, w0, den, spec)
    active = [e for e in edges if e[3] and e[2] is not None]

    # a 1-term active edge cannot re-phase, so its magnitude sqrt(w_i w_j)/den is
    # rigid; if that magnitude varies along the kernel it pins the deformation to a
    # point (no 1-parameter fiber), and the wall/reality method does not apply.
    for (_p, _q, _tgt, terms) in active:
        if len(terms) == 1:
            (i1, j1, _s1) = terms[0]
            if sp.diff(w[i1] * w[j1], t) != 0:
                return "rigid (1-term edge pins deformation)", None

    # candidate t: real roots of each edge's reality (Cayley-Menger) polynomial
    cands = set()
    for (_p, _q, tgt, terms) in active:
        tau = sp.Rational(tgt) * den ** 2
        if len(terms) == 2:
            (i1, j1, _s1), (i2, j2, _s2) = terms
            a, b = w[i1] * w[j1], w[i2] * w[j2]
            poly = sp.expand((tau - a - b) ** 2 - 4 * a * b)  # -16 Area^2
        elif len(terms) == 1:
            (i1, j1, _s1) = terms[0]
            poly = sp.expand(w[i1] * w[j1] - tau)             # single-term magnitude
        else:
            continue                                          # k>2: not prototyped
        for root in sp.solve(sp.Poly(poly, t), t):
            if root.is_real:
                cands.add(sp.nsimplify(root))

    # positivity interval
    los = [float(-Fraction(w0[i], u[i])) for i in range(m) if u[i] > 0]
    his = [float(-Fraction(w0[i], u[i])) for i in range(m) if u[i] < 0]
    tlo = max(los) if los else -1e9
    thi = min(his) if his else 1e9

    # gauge-reduced sign classes: fix the first amplitude +, enumerate the rest
    nfree = m - 1
    signs_iter = ([1] + list(s) for s in itertools.product((1, -1), repeat=nfree)) \
        if (1 << nfree) <= max_signs else None
    for tv in sorted(cands, key=float):
        if not (tlo - 1e-9 < float(tv) < thi + 1e-9):
            continue
        wv = [sp.nsimplify(w0[i] + tv * u[i]) for i in range(m)]
        if any(sp.re(x) <= 0 for x in wv):
            continue
        mags = [sp.sqrt(sp.Rational(1, den) * wv[i]) for i in range(m)]
        if signs_iter is None:
            continue
        for eps in ([1] + list(s) for s in itertools.product((1, -1), repeat=nfree)):
            amps = [eps[i] * mags[i] for i in range(m)]
            if verify_exact(n, d, spec, dets, amps):
                return "CERTIFIED", (sp.nsimplify(tv), tuple(eps))
    return "no real state", None


def main():
    recs = _load()
    jobs = [k for k, r in recs.items() if r.get("classified") == "INTERFERENCE"]
    cert, rigid, failed = [], [], []
    for key in jobs:
        status, res = solve(*key, recs=recs)
        if status.startswith("kdim") or status == "no state":
            continue
        if status == "CERTIFIED":
            cert.append((key, res[0]))
        elif status.startswith("rigid"):
            rigid.append(key)
        else:
            failed.append(key)
    genuine = len(cert) + len(failed)
    print(f"deforming (genuine 1-parameter) kdim-1 interference families: {genuine}")
    print(f"CERTIFIED real extremal state: {len(cert)}/{genuine}"
          + (f"  (uncertified: {sorted(failed)})" if failed else "  (ALL certified)"))
    print(f"rigid families (1-term edge pins the deformation; wall method N/A, "
          f"reality open): {len(rigid)} {sorted(rigid)}")
    for (sysx, tv) in sorted(cert):
        print(f"  {sysx[0]:8} v{sysx[1]:<4} t = {tv}")


if __name__ == "__main__":
    main()
