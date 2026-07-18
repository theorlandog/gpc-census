#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "scipy", "sympy"]
# ///
# Run inside the project env so gpc_census resolves from source:
#     uv run --project . scripts/solve_all.py --preflight
"""Tier A + Tier B campaign over every interference vertex in the census.

Run from the repository root:

    uv run scripts/solve_all.py             # all systems, small first
    uv run scripts/solve_all.py --preflight # v_B end to end, then stop
    uv run scripts/solve_all.py --systems 3_8,3_9

Appends one JSON record per vertex to results/data/states_interference.jsonl
with the numerical state (Tier A) and, when recognition and exact
verification succeed, the certified closed form (Tier B). Restartable:
already-recorded (system, index) pairs are skipped.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import pathlib
import sys
from fractions import Fraction

import numpy as np

from gpc_census.exactify import exactify
from gpc_census.states import _build, solve_vertex

DATA = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
OUT = DATA / "states_interference.jsonl"

# cheapest dimension first; (5,10) is the heavy tail
ORDER = ["3_8", "3_9", "4_8", "4_9", "3_10", "4_10", "5_10"]


def tasks(systems: list[str]):
    for tag in systems:
        f = DATA / "census" / f"census_{tag}_results.txt"
        if not f.exists():
            continue
        n, d = (int(x) for x in tag.split("_"))
        verts = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
        rows = [ln for ln in f.read_text().splitlines() if ln[:4].strip().isdigit()]
        for i, line in enumerate(rows):
            if "INTERFERENCE" in line:
                yield n, d, i, verts[i]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", default=",".join(ORDER))
    ap.add_argument("--preflight", action="store_true",
                    help="run v_B through both tiers and exit")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1,
                    help="parallel worker processes; 1 core per solve, BLAS pinned")
    args = ap.parse_args()
    np.random.seed(args.seed)

    if args.preflight:
        spec = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
        rec = solve_vertex(4, 9, spec)
        print("TierA", rec["status"], "support", rec.get("support_size"),
              "residual", rec.get("residual"))
        if rec["status"] != "OK":
            return 1
        ex = exactify(4, 9, spec, rec)
        print("TierB", ex["status"], ex.get("reason", ""))
        if ex["status"] == "EXACT":
            print("weights", ex["weights"], "/", ex["den"])
            for p in ex["pretty"]:
                print("  ", p)
        return 0 if ex["status"] == "EXACT" else 1

    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            try:
                r = json.loads(line)
                done.add((r["system"], r["index"]))
            except json.JSONDecodeError:
                pass

    todo = [(n, d, i, v) for n, d, i, v in tasks(args.systems.split(","))
            if (f"({n},{d})", i) not in done]
    with OUT.open("a") as out:
        if args.workers <= 1:
            cache = {}
            for t in todo:
                n, d = t[0], t[1]
                if (n, d) not in cache:
                    cache[(n, d)] = _build(d, n)
                _emit(out, _work(t, cache[(n, d)], args.seed))
        else:
            os.environ.setdefault("OMP_NUM_THREADS", "1")
            os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
            with mp.Pool(args.workers, initializer=_init_seed, initargs=(args.seed,)) as pool:
                for rec in pool.imap_unordered(_work_solo, todo):
                    _emit(out, rec)
    return 0


def _init_seed(seed: int) -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    np.random.seed(seed + os.getpid() % 10000)


def _work(t, built, seed):
    n, d, i, v = t
    np.random.seed(seed + 1000 * n + 10 * d + i)
    spec = [Fraction(s) for s in v["spectrum"]]
    rec = solve_vertex(n, d, spec, _built=built)
    rec["system"], rec["index"] = f"({n},{d})", i
    rec["integer_form"], rec["denominator"] = v["integer_form"], v["denominator"]
    if rec["status"] == "OK":
        rec["exact"] = exactify(n, d, spec, rec)
    return rec


_BUILT_CACHE: dict = {}


def _work_solo(t):
    n, d = t[0], t[1]
    if (n, d) not in _BUILT_CACHE:
        _BUILT_CACHE[(n, d)] = _build(d, n)
    return _work(t, _BUILT_CACHE[(n, d)], 0)


def _emit(out, rec):
    out.write(json.dumps(rec) + "\n")
    out.flush()
    tier_b = rec.get("exact", {}).get("status", "-")
    print(rec["system"], rec["index"], rec["status"], "support",
          rec.get("support_size"), "tierB", tier_b, flush=True)


if __name__ == "__main__":
    sys.exit(main())
