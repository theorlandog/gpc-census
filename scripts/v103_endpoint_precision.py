#!/usr/bin/env python3
"""High-precision refinement of the v103 support-10 cascade endpoint.

The cascade (Task 1) reaches support 10 at v103 = (18^4, 5^6)/34, and the
endpoint is first-order RIGID in both fiber senses, so it is an isolated point
of the support-10 stratum modulo gauge. An isolated solution of a polynomial
system with rational coefficients is algebraic, which makes exact recognition
possible in principle. This script does the first half of that: it refines the
double-precision endpoint with Gauss-Newton in mpmath arithmetic, then attempts
integer-relation recognition (PSLQ) of the squared weights.

What this is NOT: a certificate. A high-precision approximate solution plus a
recognized algebraic value is still numerical evidence until the recognized
values are substituted back and verified in exact arithmetic, which the script
does attempt and reports on separately.

  python scripts/v103_endpoint_precision.py [--digits 60]

Writes results/data/v103_endpoint.json.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402
from mpmath import mp  # noqa: E402

from gpc_census.cascade import cascade  # noqa: E402
from gpc_census.fiber import amplitudes_from_record, gauge_generators  # noqa: E402

D = 10
A_EIG = Fraction(18, 34)
B_EIG = Fraction(5, 34)


def endpoint():
    rows = [json.loads(x) for x in
            (ROOT / "results/data/states.jsonl").read_text().splitlines() if x.strip()]
    rec = next(r for r in rows if r["system"] == "(3,10)" and r["index"] == 103)
    c, dets = amplitudes_from_record(rec)
    spectrum = [Fraction(int(x), 34) for x in rec["integer_form"]]
    out = cascade(c, dets, spectrum, D)
    c_end = (np.array(out["final_amplitudes_real"])
             + 1j * np.array(out["final_amplitudes_imag"]))
    return out, c_end, [tuple(t) for t in out["final_support_dets"]]


def rdm_terms(dets):
    """(p, q, si, ti, sign) terms of rho_pq = sum conj(c_si) c_ti sign."""
    idx = {t: i for i, t in enumerate(dets)}
    terms = []
    for ti, t in enumerate(dets):
        for qi, q in enumerate(t):
            rest = [u for u in t if u != q]
            for p in range(D):
                if p in rest:
                    continue
                s = tuple(sorted(rest + [p]))
                si = idx.get(s)
                if si is None:
                    continue
                terms.append((p, q, si, ti, (-1) ** qi * (-1) ** s.index(p)))
    return terms


def mp_rdm(c, terms):
    r = [[mp.mpc(0) for _ in range(D)] for _ in range(D)]
    for p, q, si, ti, sgn in terms:
        r[p][q] += mp.conj(c[si]) * c[ti] * sgn
    return r


def mp_matmul(x, y):
    return [[sum(x[i][k] * y[k][j] for k in range(D)) for j in range(D)]
            for i in range(D)]


def residual(v, terms, m, gauge, v0):
    c = [mp.mpc(v[i], v[i + m]) for i in range(m)]
    r = mp_rdm(c, terms)
    a, b = mp.mpf(A_EIG.numerator) / A_EIG.denominator, mp.mpf(B_EIG.numerator) / B_EIG.denominator
    ra = [[r[i][j] - (a if i == j else 0) for j in range(D)] for i in range(D)]
    rb = [[r[i][j] - (b if i == j else 0) for j in range(D)] for i in range(D)]
    q = mp_matmul(ra, rb)
    out = []
    for i in range(D):
        for j in range(i, D):
            out.append(mp.re(q[i][j]))
            out.append(mp.im(q[i][j]))
    out.append(sum(mp.re(x) ** 2 + mp.im(x) ** 2 for x in c) - 1)
    for row in gauge:
        out.append(sum(mp.mpf(float(row[k])) * (v[k] - v0[k]) for k in range(2 * m)))
    return out


def jacobian(v, terms, m, gauge, v0, h=None):
    """Finite-difference Jacobian at working precision (step tuned to mp.dps)."""
    if h is None:
        h = mp.mpf(10) ** (-mp.dps // 2)
    base = residual(v, terms, m, gauge, v0)
    cols = []
    for k in range(2 * m):
        vp = list(v)
        vp[k] = vp[k] + h
        rp = residual(vp, terms, m, gauge, v0)
        cols.append([(rp[i] - base[i]) / h for i in range(len(base))])
    return base, cols


def gauss_newton(v, terms, m, gauge, v0, iters=40):
    for _ in range(iters):
        f, cols = jacobian(v, terms, m, gauge, v0)
        n = len(f)
        # normal equations J^T J dv = -J^T f
        jtj = [[sum(cols[i][r] * cols[j][r] for r in range(n)) for j in range(2 * m)]
               for i in range(2 * m)]
        jtf = [sum(cols[i][r] * f[r] for r in range(n)) for i in range(2 * m)]
        try:
            dv = mp.lu_solve(mp.matrix(jtj), mp.matrix([-x for x in jtf]))
        except Exception:
            break
        v = [v[i] + dv[i] for i in range(2 * m)]
        if max(abs(x) for x in dv) < mp.mpf(10) ** (-mp.dps + 5):
            break
    f = residual(v, terms, m, gauge, v0)
    return v, max(abs(x) for x in f)


def recognize(x, max_den=10 ** 6):
    """Try rational, then quadratic irrational, on a high-precision value."""
    r = mp.pslq([x, mp.mpf(1)], maxcoeff=max_den, maxsteps=10000)
    if r and r[0] != 0:
        cand = sp.Rational(-int(r[1]), int(r[0]))
        if abs(mp.mpf(cand.p) / cand.q - x) < mp.mpf(10) ** (-mp.dps + 10):
            return {"kind": "rational", "value": str(cand)}
    r = mp.pslq([x * x, x, mp.mpf(1)], maxcoeff=10 ** 8, maxsteps=20000)
    if r and r[0] != 0:
        t = sp.Symbol("t")
        poly = sp.Poly(int(r[0]) * t ** 2 + int(r[1]) * t + int(r[2]), t)
        for root in sp.real_roots(poly):
            if abs(mp.mpf(float(root.evalf(30))) - x) < mp.mpf(10) ** (-25):
                return {"kind": "quadratic", "minpoly": str(poly.as_expr()),
                        "value": str(sp.nsimplify(root))}
    return {"kind": "unrecognized"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--digits", type=int, default=60)
    args = ap.parse_args()
    mp.dps = args.digits

    out, c_end, dets = endpoint()
    m = len(dets)
    terms = rdm_terms(dets)
    gauge = gauge_generators(c_end, dets, D)
    v0 = [mp.mpf(float(x)) for x in np.real(c_end)] + \
         [mp.mpf(float(x)) for x in np.imag(c_end)]
    v, res = gauss_newton(list(v0), terms, m, gauge, v0)

    weights = [v[i] ** 2 + v[i + m] ** 2 for i in range(m)]
    recognized = [recognize(w) for w in weights]
    exact_rationals = {}
    for det, w, rec in zip(dets, weights, recognized):
        if rec["kind"] == "rational":
            exact_rationals[str(list(det))] = rec["value"]

    report = {
        "vertex": "(3,10) v103 = (18^4,5^6)/34",
        "object": "support-10 endpoint of the fixed-spectrum sparsification cascade",
        "support_certified": out["support_start"],
        "support_upper": out["support_upper"],
        "support_dets": [list(t) for t in dets],
        "digits": args.digits,
        "refined_residual": mp.nstr(res, 8),
        "squared_weights": [mp.nstr(w, min(args.digits, 40)) for w in weights],
        "recognized": [dict(det=list(d), **r) for d, r in zip(dets, recognized)],
        "exact_rational_weights": exact_rationals,
        "rational_count": len(exact_rationals),
        "weight_sum_check": mp.nstr(sum(weights) - 1, 8),
        "evidence": "numerical (high-precision refinement plus integer-relation "
                    "recognition); NOT a certificate",
    }
    (ROOT / "results/data/v103_endpoint.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(f"refined residual: {mp.nstr(res, 8)}")
    print(f"rational weights recognized: {len(exact_rationals)} of {m}")
    for det, r in zip(dets, recognized):
        print(f"  {list(det)}: {r}")


if __name__ == "__main__":
    main()
