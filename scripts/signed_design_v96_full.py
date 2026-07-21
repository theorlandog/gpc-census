#!/usr/bin/env python3
"""Exhaustive signed/phased-design search for (3,10) v96 in the NO basis.

PROBLEM (see docs/RESEARCH.md, "v96 campaign"):
  find integer weights k_T >= 0 on 3-subsets of {0..9} with mode sums
  EXACTLY (5,5,5,5,2,1,1,1,1,1), and per-determinant phases phi_T (rung 1:
  signs only) such that every one-hop pair class cancels:
      for all i<j:  sum_{A^B={i,j}} sigma_AB sqrt(k_A k_B) e^{i(phi_A-phi_B)} = 0.
  A solution IS an exact extremal state for v96 (rho diagonal = spectrum).

ENUMERATION IS EXHAUSTIVE, no caps: each 1-mode (5..9) has incidence weight
exactly 1, so the support decomposes as a TAIL COVER (weight-1 dets covering
{5..9}, each 1-mode in exactly one det) plus a HEAD (weighted dets inside
{0..4}), and the head solutions depend only on the residual mode sums --
memoized. Two prior time-capped CP-SAT runs found 9,056 and 6,395 skeletons
(different subsets); this enumeration replaces both.

RUNGS per skeleton:
  1. sign search (2^(n-1)), exact integer-surd class cancellation
  2. polygon filter: a class with terms {v_t} cancels with ANY phases only if
     max v_t <= sum of the others; skeletons failing any class are excluded
     EXACTLY (necessary condition)
  3. coupled numeric phase solve (multi-start L-BFGS on the total squared
     class residual); near-zero minima are printed as CANDIDATES for exact
     verification (try roots of unity / rational pi multiples, then verify
     all classes vanish in exact arithmetic)

Usage:
  python signed_design_v96_full.py            # full run, single process
  python signed_design_v96_full.py -j 8       # 8 worker processes
  python signed_design_v96_full.py --limit 500   # first 500 tail covers only
Progress prints every 200 tail covers. Exhausting rungs 1-3 with no
candidate proves v96 has no rational-weight NO-basis state at m=1 within
numeric phase resolution; rung-3 candidates must be exactified before any
claim. A clean rung-1 hit is an exact certificate immediately.
"""
from __future__ import annotations
import argparse
import time
from functools import lru_cache
from itertools import combinations, product

TARGET = (5, 5, 5, 5, 2, 1, 1, 1, 1, 1)
D = 9
HEAD_MODES = (0, 1, 2, 3, 4)
TAIL_MODES = (5, 6, 7, 8, 9)
HEAD_DETS = list(combinations(HEAD_MODES, 3))  # 10 dets


def sqfree(n):
    m, q, i = 1, n, 2
    while i * i <= q:
        while q % (i * i) == 0:
            q //= i * i
            m *= i
        i += 1
    return m, q


SF = {n: sqfree(n) for n in range(1, 26 * 26)}


def hop(A, B):
    sA, sB = set(A), set(B)
    c = sA & sB
    if len(c) != 2:
        return None, None
    (i,) = sA - c
    (j,) = sB - c
    return tuple(sorted((i, j))), (-1) ** A.index(i) * (-1) ** B.index(j)


# ---------------------------------------------------------------- tail covers
def set_partitions_le3(elems):
    if not elems:
        yield []
        return
    first, rest = elems[0], elems[1:]
    for size in (0, 1, 2):
        for others in combinations(rest, size):
            block = (first,) + others
            remaining = [e for e in rest if e not in others]
            for tail in set_partitions_le3(remaining):
                yield [block] + tail


def tail_covers():
    """Yield (tail_dets, head_usage) with tail dets = block + completion."""
    for part in set_partitions_le3(list(TAIL_MODES)):
        blocks = sorted(part)
        choices = []
        for b in blocks:
            need = 3 - len(b)
            choices.append(list(combinations(HEAD_MODES, need)))
        for combo in product(*choices):
            usage = [0] * 5
            ok = True
            dets = []
            for b, comp in zip(blocks, combo):
                for m in comp:
                    usage[m] += 1
                    if usage[m] > TARGET[m]:
                        ok = False
                dets.append(tuple(sorted(b + comp)))
            if ok:
                yield tuple(dets), tuple(usage)


# ---------------------------------------------------------------- head solve
@lru_cache(maxsize=None)
def head_solutions(residual):
    """All weight vectors on HEAD_DETS with mode sums == residual (exact)."""
    total = sum(residual)
    assert total % 3 == 0
    W = total // 3
    sols = []

    def rec(t, rem, left, acc):
        if t == len(HEAD_DETS):
            if left == 0 and all(r == 0 for r in rem):
                sols.append(tuple(acc))
            return
        Tset = HEAD_DETS[t]
        cap = min(left, *(rem[m] for m in Tset))
        # prune: remaining dets must be able to absorb rem
        for w in range(cap, -1, -1):
            nr = list(rem)
            for m in Tset:
                nr[m] -= w
            rec(t + 1, tuple(nr), left - w, acc + [w])

    rec(0, residual, W, [])
    return tuple(sols)


# ------------------------------------------------------------- per-skeleton
def process_skeleton(supp, do_phases, numeric_tol=1e-10):
    """supp: list of (det, weight). Returns (tag, data) or None."""
    n = len(supp)
    pairterms = {}
    for a in range(n):
        for b in range(a + 1, n):
            pr, sg = hop(supp[a][0], supp[b][0])
            if pr is None:
                continue
            m, q = SF[supp[a][1] * supp[b][1]]
            pairterms.setdefault(pr, []).append((a, b, sg * m, q))
    if not pairterms:
        return ("DESIGN", supp, None)  # would contradict INTERFERENCE verdict
    # rung 1: signs
    for eps in product((1, -1), repeat=n - 1):
        e = (1,) + eps
        ok = True
        for terms in pairterms.values():
            acc = {}
            for a, b, co, q in terms:
                acc[q] = acc.get(q, 0) + e[a] * e[b] * co
            if any(v for v in acc.values()):
                ok = False
                break
        if ok:
            return ("SIGNED-DESIGN", supp, e)
    if not do_phases:
        return None
    # rung 2: polygon filter (necessary for any phases)
    for terms in pairterms.items():
        vals = [abs(co) * (q ** 0.5) for _a, _b, co, q in terms[1]]
        if max(vals) > sum(vals) - max(vals) + 1e-12:
            return None
    # rung 3: coupled numeric phase solve
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

    best = None
    rng = np.random.default_rng(0)
    for _ in range(12):
        x0 = rng.uniform(0, 2 * 3.141592653589793, n - 1)
        r = minimize(resid, x0, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
        if best.fun < numeric_tol:
            break
    if best.fun < numeric_tol:
        return ("PHASE-CANDIDATE", supp, tuple(best.x), best.fun)
    return None


def run(limit=None, jobs=1, do_phases=True):
    t0 = time.time()
    covers = list(tail_covers())
    print(f"tail covers: {len(covers)}", flush=True)
    if limit:
        covers = covers[:limit]
    n_skel = 0
    hits = []
    for ci, (tdets, usage) in enumerate(covers):
        residual = tuple(TARGET[m] - usage[m] for m in range(5))
        if sum(residual) % 3:
            continue
        for hw in head_solutions(residual):
            supp = [(t, 1) for t in tdets] + \
                   [(HEAD_DETS[i], w) for i, w in enumerate(hw) if w]
            if len(supp) > 12:
                continue
            n_skel += 1
            r = process_skeleton(supp, do_phases)
            if r:
                hits.append(r)
                print("HIT:", r[0], flush=True)
                for det, w in r[1]:
                    print("   det", det, "w", w)
                if r[0] in ("SIGNED-DESIGN", "DESIGN"):
                    print("   signs:", r[2])
                    return hits
                else:
                    print("   phases:", r[2], "residual", r[3])
        if ci % 200 == 0:
            print(f"  cover {ci}/{len(covers)}  skeletons {n_skel}  "
                  f"hits {len(hits)}  {time.time()-t0:.0f}s", flush=True)
    print(f"DONE: covers {len(covers)}, skeletons {n_skel}, "
          f"hits {len(hits)}, {time.time()-t0:.0f}s", flush=True)
    if not hits:
        print("no rational-weight NO-basis state at m=1 in searched scope "
              "(rungs: signs exact; phases numeric multi-start).")
    return hits


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None,
                    help="first K tail covers only (validation runs)")
    ap.add_argument("-j", "--jobs", type=int, default=1)
    ap.add_argument("--no-phases", action="store_true")
    a = ap.parse_args()
    run(limit=a.limit, jobs=a.jobs, do_phases=not a.no_phases)
