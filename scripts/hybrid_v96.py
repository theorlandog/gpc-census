#!/usr/bin/env python3
"""Hybrid-family search for (3,10) v96, driving the polygon-target solver.

See docs/RESEARCH.md, "v96 campaign" (remaining rung a). The pure
signed/phased-design family (all one-hop classes cancel, diagonal mode sums =
spectrum) is EXHAUSTIVELY EMPTY at m=1 (350,980 skeletons; signed_design_*).
The HYBRID family relaxes that: the state's 1-RDM is BLOCK diagonal, not fully
diagonal. Modes are grouped into degenerate 2-blocks; each block's two modes
carry an EQUAL occupation that splits, under an in-block rotation, into two
distinct spectrum eigenvalues, while every off-block one-hop class still cancels.

At den 9 the only degenerate split the v96 spectrum admits is {5, 1}: two modes
at occupation 3/9 whose 2x2 block [[3/9, x], [x*, 3/9]] has eigenvalues 5/9, 1/9,
forcing |x|^2 = (3/9)^2 - (5/9)(1/9) = 4/81 (the same target the census (3,10)
interference blocks carry). Pairing t of the four 5-modes with t of the five
1-modes (t = 1..4) gives block-diagonal mode sums that DIFFER from the spectrum
(block modes sum to 3, not 5 or 1), so this is skeleton space the pure-design
search never enumerated. t = 0 is that already-empty pure design.

For each block structure this enumerates the integer weight skeletons with the
block-diagonal mode sums (total weight 9, support <= 9 dets), pre-filters to
those where every block pair actually carries a one-hop term (else the block
magnitude cannot be realized), and hands each to
polygon_target.solve(targets={block pairs: 4/81}); every returned state is
verify_exact-certified, so a hit IS an exact extremal state for v96. Bounded by
--max-seconds and --limit; a capped run prints the frontier (a documented
partial slice, not a proof).

Usage:
  python scripts/hybrid_v96.py --t 1 --max-seconds 300
  python scripts/hybrid_v96.py --t 1,2,3,4 --limit 20000
"""
from __future__ import annotations

import argparse
import sys
import time
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import polygon_target as pt  # noqa: E402

N, D = 3, 10
DEN = 9
BLOCK_TSQ = None  # set lazily as sympy Rational(4, 81)


def block_structure(t):
    """Return (ints, blocks) for t degenerate {5,1} blocks.

    Canonical labels: modes 0..3 are eigenvalue-5, mode 4 is eigenvalue-2,
    modes 5..9 are eigenvalue-1. Block i pairs mode i (a 5) with mode 5+i (a 1);
    both then carry occupation 3. Remaining 5-modes and 1-modes are singletons.
    """
    ints = [0] * D
    for m in range(4):
        ints[m] = 3 if m < t else 5
    ints[4] = 2
    for m in range(5, 10):
        ints[m] = 3 if m < 5 + t else 1
    blocks = [(i, 5 + i) for i in range(t)]
    return ints, blocks


def skeletons(ints, max_support=9):
    """Yield weight supports [(det, w), ...] with mode sums == ints, total
    weight sum(ints)/N. Ordered DFS with feasibility pruning (same shape as
    signed_design_generic)."""
    dets = list(combinations(range(D), N))
    tight = sorted(range(D), key=lambda i: ints[i])
    rank = {m: r for r, m in enumerate(tight)}
    dets.sort(key=lambda T: min(rank[m] for m in T))
    W = sum(ints) // N

    def rec(t, rem, left, acc):
        if left == 0:
            if all(r == 0 for r in rem):
                yield [(dets[i], w) for i, w in acc]
            return
        if t == len(dets) or len(acc) > max_support:
            return
        T = dets[t]
        cap = min(left, *(rem[m] for m in T))
        for w in range(cap, -1, -1):
            nr = list(rem)
            for m in T:
                nr[m] -= w
            if all(nr[m] == 0
                   or any(m in dets[u] for u in range(t + 1, len(dets)))
                   for m in range(D)):
                yield from rec(t + 1, tuple(nr), left - w,
                               acc + ([(t, w)] if w else []))

    yield from rec(0, tuple(ints), W, [])


def block_has_terms(supp_dets, p, q):
    """True if some pair of support determinants differ exactly by swapping p,q
    (a one-hop term for class {p,q})."""
    sset = {tuple(t) for t in supp_dets}
    for t in sset:
        if p in t and q not in t:
            tp = tuple(sorted(tuple(x for x in t if x != p) + (q,)))
            if tp in sset:
                return True
    return False


def run(ts, max_seconds=None, limit=None, verbose=False, dump=None):
    import json

    import sympy as sp

    global BLOCK_TSQ
    BLOCK_TSQ = sp.Rational(4, 81)
    dumpf = open(dump, "w") if dump else None
    spec = [sp.Rational(v, DEN) for v in (5, 5, 5, 5, 2, 1, 1, 1, 1, 1)]
    t0 = time.time()
    total_skel = total_solved = 0
    hits = []
    stopped = False
    for t in ts:
        ints, blocks = block_structure(t)
        targets = {pq: BLOCK_TSQ for pq in blocks}
        n_skel = n_prefilt = n_solved = 0
        for supp in skeletons(ints):
            if max_seconds and time.time() - t0 > max_seconds:
                stopped = True
                break
            if limit and total_skel + n_skel >= limit:
                stopped = True
                break
            n_skel += 1
            dets = [d for d, _w in supp]
            ks = [w for _d, w in supp]
            # every block pair must carry a one-hop term or the magnitude is 0
            if any(not block_has_terms(dets, p, q) for (p, q) in blocks):
                continue
            n_prefilt += 1
            rec = pt.solve(N, D, spec, dets, ks, DEN, targets=targets)
            n_solved += 1
            if rec is not None:
                hits.append((t, supp, rec))
                print(f"HIT t={t}: {rec['pretty']}", flush=True)
                for d, w in supp:
                    print("   det", d, "w", w)
                if dumpf:
                    dumpf.write(json.dumps({
                        "t": t, "dets": [list(d) for d in dets], "ks": ks,
                        "den": DEN, "amplitudes": rec["amplitudes"],
                        "pretty": rec["pretty"]}) + "\n")
                    dumpf.flush()
        total_skel += n_skel
        total_solved += n_solved
        print(f"t={t}: skeletons {n_skel}, block-feasible {n_prefilt}, "
              f"solved {n_solved}, hits {len([h for h in hits if h[0]==t])}, "
              f"{time.time()-t0:.0f}s", flush=True)
        if stopped:
            break
    tag = "STOPPED (cap) -- PARTIAL" if stopped else "COMPLETE"
    print(f"\n{tag}: t={list(ts)}, skeletons {total_skel}, "
          f"solver calls {total_solved}, hits {len(hits)}, "
          f"{time.time()-t0:.0f}s", flush=True)
    return hits, stopped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t", default="1", help="comma list of block counts, e.g. 1,2,3,4")
    ap.add_argument("--max-seconds", type=float, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--dump", default=None, help="write each hit as JSONL")
    a = ap.parse_args()
    ts = [int(x) for x in a.t.split(",")]
    run(ts, max_seconds=a.max_seconds, limit=a.limit, dump=a.dump)


if __name__ == "__main__":
    main()
