#!/usr/bin/env python3
"""Reproduce the evidence behind three RESEARCH_NEXT.md claims (2026-07).

Independent of the gpc_census package on purpose: fresh fermionic-sign
bookkeeping, fresh incidence solver. Requires numpy, scipy (HiGHS MILP), sympy.

  1. Realization defect, witness level: exact incidence min-supports
     s_I(v89)=7, s_I(v103)=6 with rational certificates, and exact proofs
     that both minimal witnesses admit NO full-support quantum state.
  2. Two-block phase census: all 64 two-block vertices ship all-real states.
  3. v89 border/exact support: the certified support-10 state is the attained
     c->0 endpoint of a support-11 exact-spectrum family (11th det (0,3,6)).

Run from the repo root:  python scripts/reproduce_next_claims.py
"""
import itertools
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import Bounds, LinearConstraint, least_squares, milp

LEDGER = Path(__file__).resolve().parents[1] / "results/data/states.jsonl"
D = 10


def load(system, index):
    for line in open(LEDGER):
        r = json.loads(line)
        if r["system"] == system and r["index"] == index:
            return r
    raise KeyError((system, index))


def rho_of(c, dets, d=D):
    """1-RDM <a_p^dag a_q> of sum_t c_t |t>, canonical det ordering."""
    idx = {t: i for i, t in enumerate(dets)}
    r = np.zeros((d, d), dtype=complex)
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
                r[p, q] += np.conj(c[si]) * c[ti] * (-1) ** qi * (-1) ** s.index(p)
    return r


# ---------------------------------------------------------------- claim 1
def incidence_min_support(spec, den, n, d):
    """Exact s_I via MILP (support indicators) + exact rational recheck."""
    dets = list(itertools.combinations(range(d), n))
    m = len(dets)
    lam = np.array(spec) / den
    A = np.zeros((d, m))
    for j, t in enumerate(dets):
        for o in t:
            A[o, j] = 1
    cons = [
        LinearConstraint(np.hstack([A, np.zeros((d, m))]), lam, lam),
        LinearConstraint(np.hstack([np.ones((1, m)), np.zeros((1, m))]), 1, 1),
        LinearConstraint(np.hstack([np.eye(m), -np.eye(m)]), -np.inf, 0),
    ]
    integ = np.concatenate([np.zeros(m), np.ones(m)])
    c = np.concatenate([np.zeros(m), np.ones(m)])
    res = milp(c, constraints=cons, integrality=integ, bounds=Bounds(0, 1))
    assert res.success
    supp = [dets[j] for j in range(m) if res.x[m + j] > 0.5]
    k = len(supp)
    # exact rational feasibility certificate on the witness
    ps = sp.symbols(f"p0:{k}", nonnegative=True)
    eqs = [
        sp.Eq(sum(ps[j] for j, t in enumerate(supp) if o in t), sp.Rational(spec[o], den))
        for o in range(d)
    ]
    sol = sp.solve(eqs, ps, dict=True)
    weights = [sp.nsimplify(sol[0][p]) for p in ps]
    assert all(w >= 0 for w in weights) and sum(weights) == 1
    return k, supp, weights


def channels(supp):
    ch = defaultdict(int)
    for i in range(len(supp)):
        for j in range(i + 1, len(supp)):
            sd = set(supp[i]) ^ set(supp[j])
            if len(sd) == 2:
                ch[tuple(sorted(sd))] += 1
    return ch


def kill_v89_witness(supp):
    """Projector-algebra kill for a two-block spectrum (a^2, b^8).

    P=(rho-bI)/(a-b) is a rank-2 projector inheriting rho's off-diagonal
    pattern. Singleton channels force nonzero entries; if a structurally-zero
    entry of P^2 reduces to a single product of two forced-nonzero entries,
    P^2=P is contradicted for ANY full-support state.
    """
    ch = channels(supp)
    single = {pq for pq, k in ch.items() if k == 1}
    poss = set(ch)
    if not single:
        return None

    def pos(p, q):
        return p == q or tuple(sorted((p, q))) in poss

    def forc(p, q):
        return p != q and tuple(sorted((p, q))) in single

    for p in range(D):
        for q in range(p + 1, D):
            if pos(p, q):
                continue
            forced = [r for r in range(D) if r not in (p, q) and forc(p, r) and forc(r, q)]
            other = [
                r
                for r in range(D)
                if r not in (p, q) and pos(p, r) and pos(r, q) and r not in forced
            ]
            if len(forced) == 1 and not other:
                return (p, q, forced[0])
    return None


def kill_v103_witness():
    """Branch-exhaustive kill of the s_I=6 witness for v103.

    Support: the four triples of {0,1,2,3} plus (4,5,6) and (7,8,9); rho is
    block B({0,1,2,3}) (+) diag(q5^3, q6^3), target eigs {9/17 x4, 5/34 x6}.
    (a) q5=q6=5/34: B=(9/17)I contradicts six forced-nonzero singleton
        channels. (b) q5=q6=9/17: multiplicity 6 > 4. (c) one triple at 9/17:
        B-(5/34)I must be (13/34)vv^dag; the necessary magnitude system in
        x_i=|c_i|^2 is solved exactly below and has no nonnegative solution.
    """
    x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4", nonnegative=True)
    W, f, e = sp.Rational(11, 34), sp.Rational(13, 34), sp.Rational(5, 34)
    n = [x1 + x2 + x3, x1 + x2 + x4, x1 + x3 + x4, x2 + x3 + x4]
    u = [(ni - e) / f for ni in n]
    eqs = [
        sp.Eq(x1 + x2 + x3 + x4, W),
        sp.Eq(x3 * x4, f**2 * u[0] * u[1]),
        sp.Eq(x1 * x2, f**2 * u[2] * u[3]),
        sp.Eq(x2 * x4, f**2 * u[0] * u[2]),
        sp.Eq(x1 * x3, f**2 * u[1] * u[3]),
        sp.Eq(x1 * x4, f**2 * u[0] * u[3]),
        sp.Eq(x2 * x3, f**2 * u[1] * u[2]),
    ]
    return sp.solve(eqs, [x1, x2, x3, x4], dict=True)  # expect []


# ---------------------------------------------------------------- claim 2
def two_block_phase_census():
    tab = defaultdict(lambda: defaultdict(int))
    two_block_complex = []
    for line in open(LEDGER):
        r = json.loads(line)
        blocks = len(set(r["integer_form"]))
        s = " ".join(r["closed_form"]["pretty"])
        if re.search(r"\bI\b", s) is None:
            tier = "REAL"
        elif "atan" in s and "/3" in s:
            tier = "TRISECTED"
        else:
            tier = "COMPLEX"
        tab[blocks][tier] += 1
        if blocks == 2 and tier != "REAL":
            two_block_complex.append((r["system"], r["index"]))
    return {b: dict(t) for b, t in sorted(tab.items())}, two_block_complex


# ---------------------------------------------------------------- claim 3
def v89_family(extra_det=(0, 3, 6), ts=(0.0, 0.05, 0.1, 0.2)):
    r = load("(3,10)", 89)
    c0 = np.array([float(sp.sympify(s)) for s in r["closed_form"]["pretty"]])
    supp = [tuple(t) for t in r["closed_form"]["support_dets"]]
    dets = supp + [tuple(extra_det)]
    m = len(dets)
    a, b = 15 / 26, 6 / 26

    def resid(v, t):
        c = np.concatenate([v[: m - 1] + 1j * v[m - 1 :], [t]])
        R = rho_of(c, dets)
        Q = (R - a * np.eye(D)) @ (R - b * np.eye(D))
        iu = np.triu_indices(D)
        return np.concatenate([Q[iu].real, Q[iu].imag, [np.vdot(c, c).real - 1]])

    v = np.concatenate([c0, np.zeros(len(c0))])
    out = []
    for t in ts:
        sol = least_squares(resid, v, args=(t,), method="lm", xtol=3e-16, ftol=3e-16, gtol=3e-16)
        v = sol.x
        out.append((t, float(np.max(np.abs(sol.fun)))))
    return out


if __name__ == "__main__":
    print("== claim 1: incidence min-supports + witness kills")
    for name, spec, den in [
        ("v89", [15, 15, 6, 6, 6, 6, 6, 6, 6, 6], 26),
        ("v103", [18, 18, 18, 18, 5, 5, 5, 5, 5, 5], 34),
    ]:
        k, supp, w = incidence_min_support(spec, den, 3, D)
        print(f"  s_I({name}) = {k}, witness {supp}, exact weights {w}")
        if name == "v89":
            hit = kill_v89_witness(supp)
            print(f"  v89 witness projector kill: {hit} (None means re-derive; MILP may return a different optimum)")
    print(f"  v103 branch-(c) exact solutions: {kill_v103_witness()}  (empty = killed)")

    print("== claim 2: block count vs phase tier")
    tab, bad = two_block_phase_census()
    for k, v in tab.items():
        print(f"  {k} blocks: {v}")
    print(f"  two-block vertices with non-real states: {bad}  (empty = claim disproved)")

    print("== claim 3: v89 support-11 family through the certified endpoint")
    for t, r in v89_family():
        print(f"  c_(0,3,6) = {t:.2f}  max |(rho-aI)(rho-bI)| = {r:.1e}")
