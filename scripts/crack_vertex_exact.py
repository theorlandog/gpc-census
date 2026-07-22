#!/usr/bin/env python3
"""Crack one residual vertex with the filter-free production solver, full-block.

This is the configuration that cracked (3,10) v57 and (5,10) v113 where the
default census sweep (max_clique 3, max_cliques 1, low max_card) had left them
SOLVE-FAIL: the full block search (max_cliques=0), a wider clique
(max_clique=4), and a max_card above the natural denominator, on a long budget.
Every hit is gated by the exact characteristic-polynomial identity
(certify_tier_b) and dumped for independent re-verification and ledger
injection.

Usage:
  python scripts/crack_vertex_exact.py --n 3 --d 10 --idx 103 \
      --ints 18,18,18,18,5,5,5,5,5,5 --den 34 --maxcard 44 --budget 43200
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def crack(n, d, idx, ints, den, budget, maxcard, dump_dir):
    from gpc_census.states import solve_vertex_exact_first

    spec = [Fraction(x, den) for x in ints]
    t = time.time()
    rec = solve_vertex_exact_first(n, d, spec, max_card=maxcard, max_clique=4,
                                   max_cliques=0, clique_time_budget=budget,
                                   certify_tier_b=True)
    ok = bool(rec and rec.get("status") == "OK")
    ex = rec.get("exact") if rec else None
    tag = f"OK {ex.get('status') if ex else '?'}" if ok else "FAIL"
    print(f"v{idx} ({n},{d}) den{den}: {tag} ({time.time()-t:.0f}s)", flush=True)
    if ok and ex and ex.get("status") == "EXACT":
        out = Path(dump_dir) / f"v{idx}_production.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        json.dump({"system": [n, d], "index": idx, "den": den, "spec": list(ints),
                   "closed_form": ex, "support": rec.get("support")}, open(out, "w"))
        print(f"dumped {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--d", type=int, required=True)
    ap.add_argument("--idx", type=int, required=True)
    ap.add_argument("--ints", required=True, help="comma-separated integer spectrum")
    ap.add_argument("--den", type=int, required=True)
    ap.add_argument("--budget", type=float, default=43200, help="per clique-level seconds")
    ap.add_argument("--maxcard", type=int, default=44, help="max support cardinality")
    ap.add_argument("--dump-dir", default=str(ROOT / "docs" / "hybrid_cracks"))
    a = ap.parse_args()
    ints = tuple(int(x) for x in a.ints.split(","))
    crack(a.n, a.d, a.idx, ints, a.den, a.budget, a.maxcard, a.dump_dir)


if __name__ == "__main__":
    main()
