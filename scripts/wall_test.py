#!/usr/bin/env python3
"""Real extremal states at the fiber walls of interference vertices (exact).

Generalizes the v_B fiber result (scripts/vb_fiber_ideal.py) across the census.
A vertex with a one-dimensional support-incidence kernel carries a one-parameter
family of extremal states: displacing the shipped weights along the kernel cycle
keeps the diagonal 1-RDM fixed and re-solves the single interference off-diagonal.
Where the gauge-invariant holonomy reaches 0 or pi the state is REAL. For an
active two-term block the wall condition is a law of cosines,

    w_a w_b + w_c w_d + 2 sigma sqrt(w_a w_b w_c w_d) = target * den^2,

with sigma = +-1 the combined fermionic/holonomy sign; squaring gives a rational
polynomial in the kernel parameter t whose roots are the walls. At each real root
inside the positivity range the real-sign state is certified by the exact 1-RDM
characteristic-polynomial identity (gpc_census.exactify.verify_exact).

This confirms that real extremal states at quadratic-irrational weights are GENERIC
at single-block, kernel-dimension-1 interference vertices whose shipped states are
complex: the reality of an interference vertex is not exceptional. Diagnostic
research script; scripts/ is ruff-excluded.

Run: .venv/bin/python scripts/wall_test.py
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

# the 17 single-block, kernel-dim-1 interference vertices whose shipped states are
# complex and whose fiber walls certify real (v_B plus the corpus sweep)
DEFAULT = [("(4,9)", 65), ("(3,8)", 28), ("(3,8)", 31), ("(3,9)", 34),
           ("(3,9)", 37), ("(4,10)", 34), ("(4,10)", 37)] + \
          [("(3,10)", i) for i in (40, 49, 57, 73, 81, 92, 96, 99, 101, 108)]


def _load():
    recs = {}
    for line in (ROOT / "results" / "data" / "states.jsonl").open():
        r = json.loads(line)
        recs[(r["system"], r["index"])] = r
    return recs


def certify_wall(system, index, recs):
    """Return (status, t, field) for the real fiber wall of one vertex."""
    import sympy as sp
    from gpc_census.exactify import _active_offdiagonals, verify_exact

    r = recs[(system, index)]
    n, d = (int(x) for x in system.strip("()").split(","))
    cf = r["closed_form"]
    dets = [tuple(x) for x in cf["support_dets"]]
    w0 = [int(x) for x in cf["weights"]]
    den = int(cf["den"])
    spec = [Fraction(x, r["denominator"]) for x in r["integer_form"]]
    imag = max(abs(a[2]) for a in r["support"])  # shipped state complex?

    M = sp.Matrix([[1 if m in T else 0 for m in range(d)] for T in dets])
    ns = M.T.nullspace()
    if len(ns) != 1:
        return (f"kdim={len(ns)}", None, None, imag)
    u = ns[0]
    u = u * sp.lcm([c.q for c in u])
    u = [int(x) for x in u]
    t = sp.symbols("t")
    w = [sp.Integer(w0[i]) + t * u[i] for i in range(len(dets))]

    occ, edges = _active_offdiagonals(d, dets, w0, den, spec)
    active = [e for e in edges if e[3] and e[2] is not None and len(e[3]) == 2]
    if not active:
        return ("no active 2-term block", None, None, imag)

    # try both holonomy classes: all-plus (holonomy 0), or flip one odd-u det (pi)
    for holo in (+1, -1):
        eps = [1] * len(dets)
        if holo == -1:
            j = next((i for i in range(len(dets)) if abs(u[i]) % 2 == 1), None)
            if j is None:
                continue
            eps[j] = -1
        cand = set()
        for (p, q, tgt, terms) in active:
            (i1, j1, s1), (i2, j2, s2) = terms
            a, b = w[i1] * w[j1], w[i2] * w[j2]
            poly = sp.expand((sp.Rational(tgt) * den ** 2 - a - b) ** 2 - 4 * a * b)
            for root in sp.solve(sp.Poly(poly, t), t):
                if root.is_real:
                    cand.add(sp.nsimplify(root))
        for tv in sorted(cand, key=lambda z: float(z)):
            wv = [sp.nsimplify(w0[i] + tv * u[i]) for i in range(len(dets))]
            if any(sp.re(x) <= 0 for x in wv):
                continue
            amps = [eps[i] * sp.sqrt(sp.Rational(1, den) * wv[i]) for i in range(len(dets))]
            if not all(sp.im(sp.simplify(x)) == 0 for x in amps):
                continue
            if verify_exact(n, d, spec, dets, amps):
                return ("CERTIFIED", tv, None, imag)
    return ("no real wall found", None, None, imag)


def main():
    recs = _load()
    jobs = DEFAULT
    if len(sys.argv) == 3:
        jobs = [(sys.argv[1], int(sys.argv[2]))]
    ncert = 0
    for (system, index) in jobs:
        status, tv, _, imag = certify_wall(system, index, recs)
        if status == "CERTIFIED":
            ncert += 1
        cx = "complex" if imag > 1e-9 else "real"
        print(f"{system:8} v{index:<4} shipped={cx:7} {status:22}"
              f"{'  t = ' + str(tv) if tv is not None else ''}")
    print(f"\nCERTIFIED {ncert}/{len(jobs)} real fiber-wall states "
          "(each verify_exact True, all-real amplitudes)")


if __name__ == "__main__":
    main()
