#!/usr/bin/env python3
"""Signed/phased-design search, FAST version. Same semantics as
signed_design_generic.py (same skeleton counts, same hit definitions),
with three exact accelerations applied per skeleton, in order:

  P1 LONE-TERM PRUNE: a one-hop class with a single term can never cancel
     (one nonzero number), with signs OR phases. Reject immediately. O(n^2).
  P2 PAIR-MAGNITUDE PRUNE: a 2-term class cancels (signs or phases) only if
     the two magnitudes match: m1*sqrt(q1) == m2*sqrt(q2), i.e. q1==q2 and
     m1==m2. Mismatch => reject.
  P3 GF(2) SIGN PROPAGATION: each surviving 2-term class forces
     x_a x_b x_c x_d = -sign(co1*co2) -- a linear constraint over GF(2) on
     the sign bits. Gaussian elimination reduces the sign search from
     2^(n-1) to 2^(free); inconsistency => no sign solution (phases may
     still exist; falls through to the phase rung). Remaining >=3-term
     classes are checked exactly on the reduced assignment set.

Phase rung unchanged: polygon filter (exact necessary) + multi-start
numeric solve; PHASE-CANDIDATEs require exactification before any claim.

Usage identical to signed_design_generic.py:
  python signed_design_fast.py --ints 8,8,5,5,5,1,1,1,1,1 -N 3
  python signed_design_fast.py --ints 15,15,6,6,6,6,6,6,6,6 -N 3 --max-seconds 43200
Validation: reproduces v96 exactly (350,980 skeletons, 0 hits).
"""
from __future__ import annotations
import argparse
import time
from itertools import combinations


def sqfree(n):
    m, q, i = 1, n, 2
    while i * i <= q:
        while q % (i * i) == 0:
            q //= i * i
            m *= i
        i += 1
    return m, q


def hop(A, B):
    sA, sB = set(A), set(B)
    c = sA & sB
    if len(c) != len(A) - 1:
        return None, None
    (i,) = sA - c
    (j,) = sB - c
    return tuple(sorted((i, j))), (-1) ** A.index(i) * (-1) ** B.index(j)


def gf2_solve(constraints, n):
    """constraints: list of (bitmask over n vars, parity bit). Returns
    (free_var_indices, particular_assignments_generator) or None if
    inconsistent. Assignment = int bitmask of MINUS signs."""
    rows = [(m, p) for m, p in constraints]
    pivots = []
    for col in range(n):
        piv = None
        for ri in range(len(rows)):
            if ri in [p[0] for p in pivots]:
                continue
            if rows[ri][0] >> col & 1:
                piv = ri
                break
        if piv is None:
            continue
        pivots.append((piv, col))
        pm, pp = rows[piv]
        for ri in range(len(rows)):
            if ri != piv and rows[ri][0] >> col & 1:
                rows[ri] = (rows[ri][0] ^ pm, rows[ri][1] ^ pp)
    for m, p in rows:
        if m == 0 and p == 1:
            return None
    pivot_cols = {c for _r, c in pivots}
    free = [c for c in range(n) if c not in pivot_cols]

    def assignments():
        for bits in range(1 << len(free)):
            x = 0
            for i, c in enumerate(free):
                if bits >> i & 1:
                    x |= 1 << c
            for ri, c in pivots:
                m, p = rows[ri]
                v = p ^ bin(x & (m & ~(1 << c))).count("1") % 2
                if v:
                    x |= 1 << c
            yield x

    return free, assignments


def process_skeleton(supp, SF, do_phases, tol=1e-10):
    n = len(supp)
    pairterms = {}
    for a in range(n):
        for b in range(a + 1, n):
            pr, sg = hop(supp[a][0], supp[b][0])
            if pr is None:
                continue
            p = supp[a][1] * supp[b][1]
            if p not in SF:
                SF[p] = sqfree(p)
            m, q = SF[p]
            pairterms.setdefault(pr, []).append((a, b, sg * m, q))
    if not pairterms:
        return ("DESIGN", supp, None)
    two, big = [], []
    for terms in pairterms.values():
        if len(terms) == 1:
            return None                       # P1: lone term, dead for signs+phases
        if len(terms) == 2:
            (a, b, c1, q1), (c, d, c2, q2) = terms
            if q1 != q2 or abs(c1) != abs(c2):
                return None                   # P2: magnitude mismatch, dead
            two.append((a, b, c, d, c1, c2))
        else:
            big.append(terms)
    # P3: GF(2) system on sign bits (var 0 gauge-fixed to +)
    cons = [(1, 0)]  # x_0 = + (mask bit0, parity 0)
    for a, b, c, d, c1, c2 in two:
        mask = (1 << a) ^ (1 << b) ^ (1 << c) ^ (1 << d)
        parity = 0 if (c1 > 0) != (c2 > 0) else 1   # need product of 4 = -sgn(c1c2)
        cons.append((mask, parity))
    sol = gf2_solve(cons, n)
    if sol is not None:
        free, gen = sol
        if len(free) <= 22:
            for x in gen():
                e = [(-1) ** (x >> t & 1) for t in range(n)]
                ok = True
                for terms in big:
                    acc = {}
                    for a, b, co, q in terms:
                        acc[q] = acc.get(q, 0) + e[a] * e[b] * co
                    if any(acc.values()):
                        ok = False
                        break
                if ok:
                    return ("SIGNED-DESIGN", supp, tuple(e))
    if not do_phases:
        return None
    for terms in big:
        vals = [abs(co) * (q ** 0.5) for _a, _b, co, q in terms]
        if max(vals) > sum(vals) - max(vals) + 1e-12:
            return None
    import numpy as np
    from scipy.optimize import minimize
    classes = list(pairterms.values())

    def resid(phi):
        ph = np.concatenate(([0.0], phi))
        r = 0.0
        for terms in classes:
            z = 0j
            for a, b, co, q in terms:
                z += co * (q ** 0.5) * np.exp(1j * (ph[a] - ph[b]))
            r += (z * z.conjugate()).real
        return r

    rng = np.random.default_rng(0)
    best = None
    for _ in range(12):
        r = minimize(resid, rng.uniform(0, 6.283185307, n - 1),
                     method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
        if best.fun < tol:
            break
    if best.fun < tol:
        return ("PHASE-CANDIDATE", supp, tuple(best.x), best.fun)
    return None


def run(ints, N, do_phases=True, max_seconds=None):
    d = len(ints)
    D = sum(ints) // N
    dets = list(combinations(range(d), N))
    tight = sorted(range(d), key=lambda i: ints[i])
    rank = {m: r for r, m in enumerate(tight)}
    dets.sort(key=lambda T: min(rank[m] for m in T))
    SF = {}
    t0 = time.time()
    stats = {"skel": 0, "hits": [], "stopped": False}

    def rec(t, rem, left, acc):
        if max_seconds and time.time() - t0 > max_seconds:
            stats["stopped"] = True
            return True
        if left == 0:
            if all(r == 0 for r in rem):
                supp = [(dets[i], w) for i, w in acc]
                stats["skel"] += 1
                r = process_skeleton(supp, SF, do_phases)
                if r:
                    stats["hits"].append(r)
                    print("HIT:", r[0], flush=True)
                    for det, w in r[1]:
                        print("   det", det, "w", w)
                    print("  ", r[2] if len(r) > 2 else "")
                    if r[0] == "SIGNED-DESIGN":
                        return True
                if stats["skel"] % 200000 == 0:
                    print(f"  skeletons {stats['skel']}  "
                          f"{time.time()-t0:.0f}s", flush=True)
            return False
        if t == len(dets):
            return False
        T = dets[t]
        cap = min(left, *(rem[m] for m in T))
        for w in range(cap, -1, -1):
            nr = list(rem)
            for m in T:
                nr[m] -= w
            if all(nr[m] == 0 or any(m in dets[u] for u in range(t + 1, len(dets)))
                   for m in range(d)):
                if rec(t + 1, tuple(nr), left - w, acc + ([(t, w)] if w else [])):
                    return True
        return False

    rec(0, tuple(ints), D, [])
    tag = "STOPPED (walltime) -- PARTIAL" if stats["stopped"] else "COMPLETE"
    print(f"{tag}: skeletons {stats['skel']}, hits {len(stats['hits'])}, "
          f"{time.time()-t0:.0f}s", flush=True)
    return stats


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ints", required=True)
    ap.add_argument("-N", type=int, required=True)
    ap.add_argument("--no-phases", action="store_true")
    ap.add_argument("--max-seconds", type=float, default=None)
    a = ap.parse_args()
    ints = tuple(int(x) for x in a.ints.split(","))
    run(ints, a.N, do_phases=not a.no_phases, max_seconds=a.max_seconds)
