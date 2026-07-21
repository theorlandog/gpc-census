#!/usr/bin/env python3
"""General hybrid-family search over the residual vertices.

Generalizes hybrid_v96.py to any vertex. Reuses the census's OWN block-ansatz
generator gpc_census.states.block_ansatze (2x2 blocks mixing two distinct
spectrum values with an integer Schur-Horn diagonal split) to enumerate block
structures, then for each enumerates the integer weight skeletons with that
block-diagonal degree vector and hands each to polygon_target.solve with the
block off-diagonal magnitudes as targets and every off-block one-hop class
cancelling.

This deliberately BYPASSES the min_block_count preflight, which is a false
negative on v96: it reports "no block ansatz feasible" (solve_vertex_exact_first
then bails in 0s) even though an off-block-hop-free support realizing the (5,1)
block demonstrably exists. Every returned state is verify_exact-certified.

Bounded per vertex by --max-seconds and --limit; prints the frontier per block
structure. Hits are dumped as JSONL for independent re-verification.

Usage:
  python scripts/hybrid_search.py --vertex 3,10,96 --max-seconds 300 --dump hits.jsonl
  python scripts/hybrid_search.py --all --max-seconds 120 --dump all_hits.jsonl
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import polygon_target as pt  # noqa: E402
from gpc_census.states import block_ansatze  # noqa: E402

# the 14 SOLVE-FAIL residual vertices (system, index, integer_form, den)
RESIDUAL = [
    ((3, 10), 40, (13, 9, 9, 9, 5, 5, 1, 1, 1, 1), 18),
    ((3, 10), 49, (9, 9, 5, 5, 5, 2, 1, 1, 1, 1), 13),
    ((3, 10), 57, (19, 19, 10, 10, 6, 6, 6, 6, 1, 1), 28),
    ((3, 10), 60, (8, 8, 5, 5, 5, 1, 1, 1, 1, 1), 12),
    ((3, 10), 73, (9, 9, 5, 5, 5, 5, 1, 1, 1, 1), 14),
    ((3, 10), 89, (15, 15, 6, 6, 6, 6, 6, 6, 6, 6), 26),
    ((3, 10), 96, (5, 5, 5, 5, 2, 1, 1, 1, 1, 1), 9),
    ((3, 10), 103, (18, 18, 18, 18, 5, 5, 5, 5, 5, 5), 34),
    ((4, 10), 60, (14, 11, 11, 6, 6, 6, 2, 2, 2, 0), 15),
    ((4, 10), 62, (14, 9, 9, 9, 9, 3, 3, 2, 2, 0), 15),
    ((4, 9), 40, (14, 11, 11, 6, 6, 6, 2, 2, 2), 15),
    ((4, 9), 42, (14, 9, 9, 9, 9, 3, 3, 2, 2), 15),
    ((5, 10), 113, (17, 17, 16, 8, 8, 8, 4, 4, 4, 4), 18),
    ((5, 10), 261, (14, 14, 14, 14, 10, 10, 10, 2, 1, 1), 18),
]


def skeletons(ints, N, d, max_support):
    dets = list(combinations(range(d), N))
    tight = sorted(range(d), key=lambda i: ints[i])
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
                   for m in range(d)):
                yield from rec(t + 1, tuple(nr), left - w,
                               acc + ([(t, w)] if w else []))

    yield from rec(0, tuple(ints), W, [])


def block_has_terms(dets, p, q):
    sset = {tuple(t) for t in dets}
    for t in sset:
        if p in t and q not in t:
            tp = tuple(sorted(tuple(x for x in t if x != p) + (q,)))
            if tp in sset:
                return True
    return False


def search_vertex(system, index, spec_int, den, max_blocks=1, max_seconds=None,
                  limit=None, dumpf=None, max_support=None):
    import sympy as sp

    N, d = system
    spectrum = [Fraction(x, den) for x in spec_int]
    W = sum(spec_int) // N
    if max_support is None:
        max_support = W  # weights are >= 1, so support <= total weight
    t0 = time.time()
    n_skel = n_solved = 0
    hits = []
    stopped = False
    # one lazy skeleton generator per block ansatz; round-robin a chunk from each
    # so no single ansatz starves the others (v96 cracks only on the (5,1) block).
    gens = []
    for nv, blocks in block_ansatze(N, d, spectrum, max_blocks=max_blocks):
        if not blocks:
            continue  # diagonal-only = the pure design, not the hybrid target
        block_pairs = [(u, v_) for (u, v_, _a, _b, _x2) in blocks]
        targets = {(u, v_): sp.Rational(x2, den * den)
                   for (u, v_, _a, _b, x2) in blocks}
        gens.append([skeletons(nv, N, d, max_support), nv, blocks,
                     block_pairs, targets, True])
    chunk = 400
    while gens and not hits:
        if max_seconds and time.time() - t0 > max_seconds:
            stopped = True
            break
        for g in gens:
            gen, nv, blocks, block_pairs, targets, alive = g
            done = 0
            for supp in gen:
                n_skel += 1
                done += 1
                dets = [tuple(dd) for dd, _w in supp]
                ks = [w for _dd, w in supp]
                if any(not block_has_terms(dets, p, q) for (p, q) in block_pairs):
                    if done >= chunk:
                        break
                    continue
                rec = pt.solve(N, d, spectrum, dets, ks, den, targets=targets)
                n_solved += 1
                if rec is not None:
                    hits.append((nv, blocks, dets, ks, rec))
                    if dumpf:
                        dumpf.write(json.dumps({
                            "system": list(system), "index": index,
                            "den": den, "spec": list(spec_int),
                            "nv": list(nv), "blocks": [list(b) for b in blocks],
                            "dets": [list(x) for x in dets], "ks": ks,
                            "amplitudes": rec["amplitudes"],
                            "pretty": rec["pretty"]}) + "\n")
                        dumpf.flush()
                    break
                if done >= chunk:
                    break
            else:
                g[5] = False  # generator exhausted
            if hits or (limit and n_skel >= limit):
                stopped = stopped or (limit and n_skel >= limit)
                break
        gens = [g for g in gens if g[5]]
    dt = time.time() - t0
    tag = "CRACKED" if hits else ("PARTIAL" if stopped else "no hit (exhausted slice)")
    print(f"({N},{d}) v{index}: {tag}  skeletons {n_skel}, solver calls "
          f"{n_solved}, {dt:.0f}s", flush=True)
    if hits:
        print(f"    -> {hits[0][4]['pretty']}", flush=True)
    return hits, stopped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vertex", help="N,d,index (one residual vertex)")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--max-blocks", type=int, default=1)
    ap.add_argument("--max-seconds", type=float, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--max-support", type=int, default=None,
                    help="cap the support cardinality (nonzero dets). Extremal "
                         "states are sparse, so a small cap (e.g. 10) shrinks "
                         "the skeleton search by orders of magnitude. Default: "
                         "the full weight W, which is intractable at high den.")
    ap.add_argument("--dump", default=None)
    a = ap.parse_args()
    dumpf = open(a.dump, "w") if a.dump else None
    if a.vertex:
        N, d, idx = (int(x) for x in a.vertex.split(","))
        row = next(r for r in RESIDUAL if r[0] == (N, d) and r[1] == idx)
        search_vertex(row[0], row[1], row[2], row[3],
                      max_blocks=a.max_blocks, max_seconds=a.max_seconds,
                      limit=a.limit, dumpf=dumpf, max_support=a.max_support)
    elif a.all:
        cracked = 0
        for (system, index, spec_int, den) in RESIDUAL:
            hits, _ = search_vertex(system, index, spec_int, den,
                                    max_blocks=a.max_blocks,
                                    max_seconds=a.max_seconds, limit=a.limit,
                                    dumpf=dumpf, max_support=a.max_support)
            cracked += bool(hits)
        print(f"\nCRACKED {cracked}/{len(RESIDUAL)} residual vertices "
              f"(this bounded slice)")
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
