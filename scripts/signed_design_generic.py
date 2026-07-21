#!/usr/bin/env python3
"""Signed/phased-design search, generalized to any (N,d) vertex.

Same problem as signed_design_v96_full.py (see docs/RESEARCH.md, "v96
campaign"), for an arbitrary integer spectrum at m=1: integer weights with
mode sums EXACTLY the integer form, phases cancelling every one-hop class.

Enumeration: exhaustive DFS over determinants ordered so the tightest modes
bind first, with budget pruning. Reproduces the v96 tail-cover result
exactly (350,980 skeletons, 0 hits) as a cross-check. For roots with no
incidence-1 modes (v89, v103) the space is much larger: run with a
walltime and read the coverage line -- the DFS prints the frontier so a
partial run is a documented partial slice, not a silent cap.

Usage:
  python signed_design_generic.py --ints 5,5,5,5,2,1,1,1,1,1 -N 3
  python signed_design_generic.py --ints 9,9,5,5,5,2,1,1,1,1 -N 3      # v49
  python signed_design_generic.py --ints 19,19,10,10,6,6,6,6,1,1 -N 3  # v57
  python signed_design_generic.py --ints 15,15,6,6,6,6,6,6,6,6 -N 3    # v89
  python signed_design_generic.py --ints 18,18,18,18,5,5,5,5,5,5 -N 3  # v103
  python signed_design_generic.py --ints 14,14,14,14,10,10,10,2,1,1 -N 5  # (5,10) v261
Options: --no-phases (rung 1 only), --max-seconds W (walltime).
"""
from __future__ import annotations
import argparse
import time
from itertools import combinations, product


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
    for eps in product((1, -1), repeat=n - 1):
        e = (1,) + eps
        ok = True
        for terms in pairterms.values():
            acc = {}
            for a, b, co, q in terms:
                acc[q] = acc.get(q, 0) + e[a] * e[b] * co
            if any(acc.values()):
                ok = False
                break
        if ok:
            return ("SIGNED-DESIGN", supp, e)
    if not do_phases:
        return None
    for terms in pairterms.values():
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
    # order dets so those touching the tightest (smallest-budget) modes come first
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
                    if r[0] in ("SIGNED-DESIGN",):
                        return True
                if stats["skel"] % 100000 == 0:
                    print(f"  skeletons {stats['skel']}  "
                          f"{time.time()-t0:.0f}s", flush=True)
            return False
        if t == len(dets):
            return False
        T = dets[t]
        cap = min(left, *(rem[m] for m in T))
        # prune: remaining dets must be able to cover every remaining mode
        for w in range(cap, -1, -1):
            nr = list(rem)
            for m in T:
                nr[m] -= w
            # feasibility: any mode with positive residual must appear in a later det
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
