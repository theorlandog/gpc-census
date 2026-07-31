#!/usr/bin/env python3
"""Exactify the sparsification-cascade endpoints across the census.

The v103 support-10 endpoint turned out to be exactly recognizable: refine in
high precision, run integer-relation detection on the squared weights, gauge-fix
the phases, substitute back, and verify the characteristic-polynomial identity
in exact arithmetic. This script runs that recipe on every state the cascade
sparsified, so improved support bounds can be certified rather than left at [N].

A state is reported exact only when all three hold in exact arithmetic: the
squared weights are algebraic of degree at most two and sum to 1, the
gauge-fixed phases are multiples of pi (the state is real up to gauge), and
det(rho - x I) equals the target characteristic polynomial. Anything else is
reported with the reason it fell short, never as a bound.

  python scripts/exactify_cascade.py [--digits 50] [--jobs 4]

Writes results/data/cascade_exact.json.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp_pool
import pathlib
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402
from mpmath import mp  # noqa: E402

from gpc_census.fiber import gauge_generators  # noqa: E402

CASCADE = ROOT / "results/data/cascade.jsonl"
STATES = ROOT / "results/data/states.jsonl"


def rdm_terms(dets, d):
    idx = {t: i for i, t in enumerate(dets)}
    terms = []
    for ti, t in enumerate(dets):
        for qi, q in enumerate(t):
            rest = [u for u in t if u != q]
            for p in range(d):
                if p in rest:
                    continue
                s = tuple(sorted(rest + [p]))
                si = idx.get(s)
                if si is None:
                    continue
                terms.append((p, q, si, ti, (-1) ** qi * (-1) ** s.index(p)))
    return terms


def mp_residual(v, terms, m, d, lams, gauge, v0):
    c = [mp.mpc(v[i], v[i + m]) for i in range(m)]
    r = [[mp.mpc(0) for _ in range(d)] for _ in range(d)]
    for p, q, si, ti, sgn in terms:
        r[p][q] += mp.conj(c[si]) * c[ti] * sgn
    q_mat = [[mp.mpc(1) if i == j else mp.mpc(0) for j in range(d)] for i in range(d)]
    for lam in lams:
        shifted = [[r[i][j] - (lam if i == j else 0) for j in range(d)] for i in range(d)]
        q_mat = [[sum(q_mat[i][k] * shifted[k][j] for k in range(d)) for j in range(d)]
                 for i in range(d)]
    out = []
    for i in range(d):
        for j in range(i, d):
            out.append(mp.re(q_mat[i][j]))
            out.append(mp.im(q_mat[i][j]))
    out.append(sum(mp.re(x) ** 2 + mp.im(x) ** 2 for x in c) - 1)
    for row in gauge:
        out.append(sum(mp.mpf(float(row[k])) * (v[k] - v0[k]) for k in range(2 * m)))
    return out


def gauss_newton(v, terms, m, d, lams, gauge, v0, iters=30):
    h = mp.mpf(10) ** (-mp.dps // 2)
    for _ in range(iters):
        f = mp_residual(v, terms, m, d, lams, gauge, v0)
        cols = []
        for k in range(2 * m):
            vp = list(v)
            vp[k] = vp[k] + h
            fp = mp_residual(vp, terms, m, d, lams, gauge, v0)
            cols.append([(fp[i] - f[i]) / h for i in range(len(f))])
        n = len(f)
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
    return v, max(abs(x) for x in mp_residual(v, terms, m, d, lams, gauge, v0))


def independent_rows(a, rank):
    """Indices of `rank` linearly independent rows of a, by pivoted QR."""
    from scipy.linalg import qr

    _, _, piv = qr(a.T, pivoting=True, mode="economic")
    return sorted(int(i) for i in piv[:rank])


def as_rational(x, maxcoeff=10 ** 9):
    r = mp.pslq([x, mp.mpf(1)], maxcoeff=maxcoeff, maxsteps=20000)
    if not r or r[0] == 0:
        return None
    cand = sp.Rational(-int(r[1]), int(r[0]))
    if abs(mp.mpf(cand.p) / cand.q - x) < mp.mpf(10) ** (-mp.dps + 10):
        return cand
    return None


def as_quadratic(x, maxcoeff=10 ** 9):
    """Recognize a real degree-two algebraic number by PSLQ and select its root."""
    relation = mp.pslq(
        [x * x, x, mp.mpf(1)],
        maxcoeff=maxcoeff,
        maxsteps=20000,
    )
    if not relation or relation[0] == 0:
        return None
    a, b, c = (int(value) for value in relation)
    common = abs(int(sp.gcd(a, sp.gcd(b, c))))
    if common:
        a, b, c = a // common, b // common, c // common
    if a < 0:
        a, b, c = -a, -b, -c
    discriminant = sp.Integer(b * b - 4 * a * c)
    if discriminant < 0:
        return None
    roots = [
        sp.radsimp((-sp.Integer(b) + sign * sp.sqrt(discriminant)) / (2 * a))
        for sign in (1, -1)
    ]
    tolerance = mp.mpf(10) ** (-mp.dps + 10)
    matches = [
        root
        for root in roots
        if abs(mp.mpf(str(sp.N(root, mp.dps))) - x) < tolerance
    ]
    if len(matches) != 1:
        return None
    return matches[0]


def recognize_weight(x):
    """Return an exact rational or quadratic expression, or None."""
    rational = as_rational(x)
    return rational if rational is not None else as_quadratic(x)


def squarefree_part(value):
    """Signed square-free part of a nonzero integer."""
    value = int(value)
    sign = -1 if value < 0 else 1
    out = sign
    for prime, exponent in sp.factorint(abs(value)).items():
        if exponent % 2:
            out *= int(prime)
    return out


def exact_rho(amps, dets, d):
    idx = {t: i for i, t in enumerate(dets)}
    rho = sp.zeros(d, d)
    for ti, t in enumerate(dets):
        for qi, q in enumerate(t):
            rest = [u for u in t if u != q]
            for p in range(d):
                if p in rest:
                    continue
                s = tuple(sorted(rest + [p]))
                si = idx.get(s)
                if si is None:
                    continue
                rho[p, q] += amps[si] * amps[ti] * (-1) ** qi * (-1) ** s.index(p)
    return sp.simplify(rho)


def process(args):
    rec, digits, which = args
    mp.dps = digits
    system, index = rec["system"], rec["index"]
    d = int(system.strip("()").split(",")[1])
    pre = f"support_upper_{which}"
    if f"{pre}_dets" in rec:                      # two-locus schema
        dets = [tuple(t) for t in rec[f"{pre}_dets"]]
        c = (np.array(rec[f"{pre}_amplitudes_real"])
             + 1j * np.array(rec[f"{pre}_amplitudes_imag"]))
        upper = rec[pre]
    else:                                          # single-locus (spectrum) schema
        dets = [tuple(t) for t in rec["final_support_dets"]]
        c = (np.array(rec["final_amplitudes_real"])
             + 1j * np.array(rec["final_amplitudes_imag"]))
        upper = rec["support_upper"]
    m = len(dets)

    src = next(json.loads(x) for x in STATES.read_text().splitlines()
               if x.strip() and json.loads(x)["system"] == system
               and json.loads(x)["index"] == index)
    den = int(src["denominator"])
    spectrum = [Fraction(int(x), den) for x in src["integer_form"]]
    lams = sorted({mp.mpf(v.numerator) / v.denominator for v in spectrum}, reverse=True)

    out = {"system": system, "index": index,
           "locus": "natural_orbital" if which == "no" else "fixed_spectrum",
           "support_certified": rec["support_certified"],
           "support_upper": upper,
           "support_dets": [list(t) for t in dets]}
    try:
        terms = rdm_terms(dets, d)
        gauge = gauge_generators(c, dets, d)
        v0 = [mp.mpf(float(x)) for x in np.real(c)] + [mp.mpf(float(x)) for x in np.imag(c)]
        v, res = gauss_newton(list(v0), terms, m, d, lams, gauge, v0)
        out["refined_residual"] = mp.nstr(res, 6)
        if res > mp.mpf(10) ** (-digits + 20):
            out["status"] = "REFINEMENT-FAILED"
            return out

        weights = [v[i] ** 2 + v[i + m] ** 2 for i in range(m)]
        exact_weights = [recognize_weight(w) for w in weights]
        if any(weight is None for weight in exact_weights):
            out["status"] = "WEIGHTS-NOT-DEGREE-TWO"
            out["recognized"] = [
                str(weight) if weight is not None else None
                for weight in exact_weights
            ]
            return out
        out["squared_weights"] = [str(weight) for weight in exact_weights]
        out["squared_weight_minpolynomials"] = [
            str(sp.minpoly(weight)) for weight in exact_weights
        ]
        out["squared_weight_max_degree"] = max(
            int(sp.Poly(sp.minpoly(weight)).degree()) for weight in exact_weights
        )
        if sp.simplify(sum(exact_weights) - 1) != 0:
            out["status"] = "WEIGHTS-DO-NOT-SUM-TO-ONE"
            return out
        irrational = [
            weight for weight in exact_weights if not bool(weight.is_Rational)
        ]
        if irrational:
            primitive = sp.to_number_field(irrational)
            polynomial = primitive.minpoly
            discriminant = sp.discriminant(polynomial.as_expr())
            out["squared_weight_field"] = (
                f"Q(sqrt({squarefree_part(discriminant)}))"
            )
            out["squared_weight_field_generator"] = str(primitive.as_expr())
            out["squared_weight_field_minpoly"] = str(polynomial.as_expr())
        else:
            out["squared_weight_field"] = "Q"

        # Gauge-fix: is the state real up to the orbital-phase gauge? The gauge
        # acts by theta -> theta + A phi, so the question is whether theta lies
        # in image(A) + pi Z^m. A plain least-squares fit is the wrong test: it
        # spreads the residual over every coordinate instead of concentrating it
        # on the gauge-invariant directions. Fit exactly on a maximal
        # independent set of rows instead, so the residual lives only on the
        # remaining coordinates, where it must be a multiple of pi.
        theta = np.angle(np.array([complex(v[i], v[i + m]) for i in range(m)]))
        amat = np.array([[1.0 if o in t else 0.0 for o in range(d)] for t in dets])
        rank = int(np.linalg.matrix_rank(amat))
        basis = independent_rows(amat, rank)
        phi, *_ = np.linalg.lstsq(amat[basis], theta[basis], rcond=None)
        resid_theta = theta - amat @ phi
        if np.max(np.abs(resid_theta[basis])) > 1e-6:
            out["status"] = "GAUGE-FIT-FAILED"
            return out
        k = np.round(resid_theta / np.pi)
        if np.max(np.abs(resid_theta - k * np.pi)) > 1e-6:
            out["status"] = "STATE-NOT-REAL-UP-TO-GAUGE"
            out["gauge_fixed_phases_over_pi"] = [float(x) for x in resid_theta / np.pi]
            return out
        signs = [1 if int(x) % 2 == 0 else -1 for x in k]
        out["signs"] = signs

        amps = [
            sign * sp.sqrt(weight)
            for sign, weight in zip(signs, exact_weights)
        ]
        rho = exact_rho(amps, dets, d)
        x = sp.Symbol("x")
        charpoly = sp.expand(rho.charpoly(x).as_expr())
        target = sp.expand(sp.prod([(x - sp.Rational(v.numerator, v.denominator))
                                    for v in spectrum]))
        # The characteristic polynomial is CONJUGATION INVARIANT, so it is a
        # consistency check and never the acceptance criterion: it is satisfied
        # by every state in the U(d) orbit, including states whose 1-RDM is
        # rotated away from the diagonal. Acceptance requires attainment in the
        # census convention: rho exactly diagonal AND diag(rho) = lambda.
        out["charpoly_identity_exact"] = bool(sp.simplify(charpoly - target) == 0)
        offdiag = [sp.simplify(rho[i, j]) for i in range(d) for j in range(d)
                   if i != j]
        nonzero_off = [str(e) for e in offdiag if e != 0]
        out["rdm_offdiagonal_nonzero"] = nonzero_off[:8]
        out["rdm_is_exactly_diagonal"] = not nonzero_off
        diag = [sp.simplify(rho[i, i]) for i in range(d)]
        lam = [sp.Rational(v.numerator, v.denominator) for v in spectrum]
        out["rdm_diagonal"] = [str(e) for e in diag]
        out["diag_equals_lambda"] = bool(all(a == b for a, b in zip(diag, lam)))
        if not irrational:
            out["state_denominator"] = int(
                sp.ilcm(*[weight.q for weight in exact_weights])
            )

        if not out["charpoly_identity_exact"]:
            out["status"] = "CHARPOLY-MISMATCH"
            return out
        if not out["rdm_is_exactly_diagonal"]:
            # a genuine exact state on the spectrum locus, but NOT an attainer
            # in the diagonal sense: it bounds s_Q^free, not s_Q^NO
            out["status"] = "SPECTRUM-ONLY-NOT-DIAGONAL"
            out["bounds"] = "s_Q_free"
            out["evidence"] = "EXACT for the spectrum locus only"
            return out
        if not out["diag_equals_lambda"]:
            out["status"] = "DIAGONAL-BUT-WRONG-OCCUPATIONS"
            return out
        out["status"] = "CERTIFIED"
        out["bounds"] = "s_Q_NO"
        out["evidence"] = "EXACT"
    except Exception as exc:
        out["status"] = f"ERROR: {type(exc).__name__}: {exc}"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--digits", type=int, default=50)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    rows = [json.loads(x) for x in CASCADE.read_text().splitlines() if x.strip()]
    todo = []
    for r in rows:
        if "support_upper_free_dropped_by" in r:
            for which in ("no", "free"):
                if r.get(f"support_upper_{which}_dropped_by", 0) > 0:
                    todo.append((r, args.digits, which))
        elif r.get("support_dropped_by", 0) > 0:
            todo.append((r, args.digits, "free"))
    print(f"exactifying {len(todo)} cascade endpoints at {args.digits} digits "
          f"({sum(1 for x in todo if x[2] == 'no')} natural-orbital, "
          f"{sum(1 for x in todo if x[2] == 'free')} free-basis)", flush=True)

    if args.jobs == 1:
        results = [process(item) for item in todo]
    else:
        with mp_pool.Pool(args.jobs) as pool:
            results = pool.map(process, todo)

    results.sort(key=lambda r: (r["locus"], r["system"], r["index"]))
    certified = [r for r in results if r["status"] == "CERTIFIED"]
    spectrum_only = [r for r in results if r["status"] == "SPECTRUM-ONLY-NOT-DIAGONAL"]
    tally = {}
    for r in results:
        key = r["status"].split(":")[0]
        tally[key] = tally.get(key, 0) + 1
    report = {
        "note": "Exact certification of cascade endpoints with separate "
                "natural-orbital and fixed-spectrum statuses. CERTIFIED bounds "
                "s_Q_NO; SPECTRUM-ONLY-NOT-DIAGONAL bounds s_Q_free. Every "
                "other status is reported without an exact support claim.",
        "digits": args.digits,
        "endpoints_attempted": len(results),
        "certified_s_Q_NO": len(certified),
        "exact_spectrum_only_s_Q_free": len(spectrum_only),
        "acceptance_criterion": "rho exactly diagonal AND diag(rho) = lambda, "
                                "in exact arithmetic. The characteristic-"
                                "polynomial identity is conjugation invariant "
                                "and is recorded as a consistency check only.",
        "status_tally": tally,
        "certified_s_Q_NO_bounds": [
            {"system": r["system"], "index": r["index"],
             "support_certified_library": r["support_certified"],
             "support_certified_exact": r["support_upper"],
             "state_denominator": r.get("state_denominator"),
             "squared_weight_field": r["squared_weight_field"]}
            for r in certified],
        "exact_s_Q_free_bounds": [
            {"system": r["system"], "index": r["index"],
             "support_certified_library": r["support_certified"],
             "support_exact_free": r["support_upper"],
             "state_denominator": r.get("state_denominator"),
             "squared_weight_field": r["squared_weight_field"],
             "rdm_offdiagonal_nonzero": r["rdm_offdiagonal_nonzero"]}
            for r in spectrum_only],
        "results": results,
    }
    (ROOT / "results/data/cascade_exact.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in
                      ("endpoints_attempted", "certified_s_Q_NO",
                       "exact_spectrum_only_s_Q_free", "status_tally")}, indent=2))


if __name__ == "__main__":
    main()
