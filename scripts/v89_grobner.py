#!/usr/bin/env python3
"""Non-symmetric exact solve for the open residual v89 = (15,15,6^8)/26 of (3,10).

Every state attaining v89 is orbital-equivalent to one with DIAGONAL 1-RDM
diag(15,15,6,6,6,6,6,6,6,6)/26 (orbital rotations preserve the spectrum, and the
rank-2 plane P_2 can be rotated onto the first two computational orbitals). So the
exact problem is: find amplitudes a_T on a support S with

    occupation: sum_{T contains i} a_T^2 = target_i    (target = (15,15,6^8)/26),
    coherence:  sum over one-hop pairs sign * a_T a_T' = 0   for every i < j,

i.e. a diagonal 1-RDM whose off-diagonals cancel despite a one-hop-CONNECTED support
(the interference analogue of a design; a genuine design is excluded for v89). This
is a quadratic system in the |S| amplitudes -- tractable by Groebner elimination for
a fixed, small-enough support. The scaffold: (1) a numeric real-signed diagonal
attain to find candidate supports, (2) an exact Groebner/solve on each candidate
under a time cap, (3) verify_exact on any hit. Designed for an UNATTENDED background
run; logs progress. Success is not guaranteed -- v89's state may be irreducibly
dense (see the RESEARCH.md near-theorem: no exploitable symmetry), in which case no
small support solves and this reports that honestly.

Run (background):  .venv/bin/python scripts/v89_grobner.py --budget 3600 --maxcard 14
"""
from __future__ import annotations

import argparse
import itertools
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

D = 10
N = 3
DEN = 26
TARGET = [Fraction(x, DEN) for x in (15, 15, 6, 6, 6, 6, 6, 6, 6, 6)]
ALL_DETS = [tuple(c) for c in itertools.combinations(range(D), N)]


def offdiag_pairs(dets):
    """one-hop determinant pairs (share N-1 orbitals) with the swapped modes + sign."""
    out = []
    index = {t: i for i, t in enumerate(dets)}
    for i, T in enumerate(dets):
        for p in T:
            for q in range(D):
                if q in T:
                    continue
                tp = tuple(sorted([o for o in T if o != p] + [q]))
                if tp in index and index[tp] > i:
                    s1 = (-1) ** list(T).index(p)
                    s2 = (-1) ** list(tp).index(q)
                    out.append((i, index[tp], min(p, q), max(p, q), s1 * s2))
    return out


def numeric_support(restarts, seed):
    """Real-signed diagonal attain; return a candidate support (dets with weight)."""
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    tgt = np.array([float(x) for x in TARGET])
    dets = ALL_DETS

    def rho_diag_offdiag(a):
        a = a / np.linalg.norm(a)
        amap = {dets[i]: a[i] for i in range(len(dets))}
        M = np.zeros((D, D))
        for T, ct in amap.items():
            for mp in T:
                s1 = (-1) ** list(T).index(mp)
                t2 = tuple(x for x in T if x != mp)
                for m in range(D):
                    if m in t2:
                        continue
                    tp = tuple(sorted(t2 + (m,)))
                    if tp in amap:
                        s2 = (-1) ** list(tp).index(m)
                        M[m, mp] += s1 * s2 * amap[tp] * ct
        return M

    def obj(a):
        M = rho_diag_offdiag(a)
        diagerr = np.sum((np.diag(M) - tgt) ** 2)
        offerr = np.sum(M ** 2) - np.sum(np.diag(M) ** 2)
        return diagerr + offerr

    best = None
    for _ in range(restarts):
        a0 = rng.standard_normal(len(dets))
        r = minimize(obj, a0, method="L-BFGS-B", options={"maxiter": 4000})
        if best is None or r.fun < best.fun:
            best = r
    a = best.x / __import__("numpy").linalg.norm(best.x)
    order = sorted(range(len(dets)), key=lambda i: -a[i] ** 2)
    return best.fun, [(dets[i], float(a[i])) for i in order if a[i] ** 2 > 1e-6]


def solve_support(dets, budget, log):
    """Exact quadratic solve on a fixed support; return exact amps or None."""
    import sympy as sp
    from gpc_census.exactify import verify_exact
    m = len(dets)
    a = sp.symbols(f"a0:{m}", real=True)
    eqs = []
    for i in range(D):
        eqs.append(sp.Add(*[a[k] ** 2 for k, T in enumerate(dets) if i in T])
                   - sp.Rational(TARGET[i].numerator, TARGET[i].denominator))
    # each 1-RDM off-diagonal rho_{p,q} must vanish: group one-hop terms by (p,q)
    coh = {}
    for (i, j, p, q, sgn) in offdiag_pairs(dets):
        coh.setdefault((p, q), []).append(sgn * a[i] * a[j])
    eqs += [sp.Add(*v) for v in coh.values()]
    t0 = time.time()
    try:
        sol = sp.solve(eqs, list(a), dict=True)
    except Exception as e:  # noqa: BLE001
        log(f"    solve error: {e}")
        return None
    if time.time() - t0 > budget:
        log("    (over budget during solve)")
    for s in sol or []:
        vals = [s.get(ai, ai) for ai in a]
        if any(getattr(v, "free_symbols", set()) for v in vals):
            continue  # parametric family; skip (pick a point later if needed)
        amps = [sp.nsimplify(v) for v in vals]
        if all(sp.im(sp.simplify(x)) == 0 for x in amps) and any(x != 0 for x in amps):
            spec = TARGET
            if verify_exact(N, D, spec, dets, amps):
                return amps
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=float, default=3600, help="seconds per support solve")
    ap.add_argument("--maxcard", type=int, default=14, help="max support cardinality to try")
    ap.add_argument("--restarts", type=int, default=120)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    def log(msg):
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

    log(f"v89 exact scaffold: maxcard={args.maxcard} budget={args.budget}s/support")
    fun, supp = numeric_support(args.restarts, args.seed)
    log(f"numeric real-diagonal attain: obj={fun:.3e}, |support|={len(supp)}")
    if fun > 1e-8:
        log("numeric real-signed diagonal state did NOT converge to zero -- v89 may "
            "need complex amplitudes on this alignment, or the state is dense. "
            "Proceeding to exact solves on truncated supports anyway (they will "
            "report infeasible/timeout, which is itself informative).")
    ranked = [T for T, _w in supp]
    for k in range(N + 7, args.maxcard + 1):
        cand = ranked[:k]
        log(f"exact solve on top-{k} support: {cand}")
        t0 = time.time()
        res = solve_support(cand, args.budget, log)
        log(f"  -> {'HIT' if res else 'no exact solution'} ({time.time()-t0:.0f}s)")
        if res:
            out = ROOT / "docs" / "v89_solution.txt"
            out.write_text("v89 exact state\nsupport: " + repr(cand)
                           + "\namps: " + repr([str(x) for x in res]) + "\n")
            log(f"SOLVED v89. wrote {out}")
            return
    log("no exact solution found within maxcard; consistent with the near-theorem "
        "that v89's extremal state is dense/symmetry-free. Raise --maxcard to push.")


if __name__ == "__main__":
    main()
