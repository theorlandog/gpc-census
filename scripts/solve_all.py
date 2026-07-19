"""Tier A + Tier B state-construction campaign over the classified census.

Run from the repository root:

    uv run scripts/solve_all.py             # interference vertices, small first
    uv run scripts/solve_all.py --all       # every vertex, routed by verdict
    uv run scripts/solve_all.py --preflight # v_B end to end, then stop
    uv run scripts/solve_all.py --systems 3_8,3_9

Each vertex is routed by the classification we already computed: DESIGN-INT
vertices are constructed directly from their design witness (no iterative
solve, exact by construction), while DESIGN-REAL and INTERFERENCE vertices go
through the state solver (cascade by default). Appends one JSON record per
vertex with the numerical state (Tier A) and, when recognition and exact
verification succeed, the certified closed form (Tier B). Default output is
results/data/states_interference.jsonl (interference only); --all writes the
full census to results/data/states.jsonl. Restartable: already-recorded
(system, index) pairs are skipped.
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

from checkpoint import Checkpointer, add_arguments
from gpc_census.exactify import exactify
from gpc_census.states import (_build, solve_design_vertex, solve_vertex,
                               solve_vertex_exact_first)

DATA = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
OUT = DATA / "states_interference.jsonl"
OUT_ALL = DATA / "states.jsonl"

# cheapest dimension first; (5,10) is the heavy tail
ORDER = ["3_8", "3_9", "4_8", "4_9", "3_10", "4_10", "5_10"]


def _verdict(line: str) -> str:
    for v in ("DESIGN-INT", "DESIGN-REAL", "INTERFERENCE"):
        if v in line:
            return v
    return "UNKNOWN"


def tasks(systems: list[str], all_vertices: bool = False):
    for tag in systems:
        f = DATA / "census" / f"census_{tag}_results.txt"
        if not f.exists():
            continue
        n, d = (int(x) for x in tag.split("_"))
        verts = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
        rows = [ln for ln in f.read_text().splitlines() if ln[:4].strip().isdigit()]
        for i, line in enumerate(rows):
            verdict = _verdict(line)
            if all_vertices or verdict == "INTERFERENCE":
                yield n, d, i, verts[i], verdict


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", default=",".join(ORDER))
    ap.add_argument("--preflight", action="store_true",
                    help="reconstruct and certify v_B end to end, then exit")
    ap.add_argument("--legacy-preflight", action="store_true",
                    help="Tier-A-only gate via the historical attain solver; "
                         "slow and does not certify interference vertices")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1,
                    help="parallel worker processes; 1 core per solve, BLAS pinned")
    ap.add_argument("--solver", choices=("cascade", "weights-first", "attain"),
                    default="cascade",
                    help="cascade (default): weights-first, then attain on its "
                         "failures, for the strongest available result per "
                         "vertex; weights-first: fast block-ansatz solve that "
                         "certifies closed forms but fails off its ansatz "
                         "family; attain: historical Tier-A alternating "
                         "projection (rarely certifies interference vertices)")
    ap.add_argument("--max-card", type=int, default=14,
                    help="max support size the weights-first solver enumerates")
    ap.add_argument("--max-blocks", type=int, default=2,
                    help="max 2x2 natural-basis blocks in a weights-first ansatz")
    ap.add_argument("--all", action="store_true",
                    help="process every vertex routed by verdict (DESIGN-INT "
                         "built directly from its witness), not just "
                         "interference; writes the full census to states.jsonl")
    add_arguments(ap, default_out=OUT)
    args = ap.parse_args()
    np.random.seed(args.seed)
    if args.all and args.out == str(OUT):
        args.out = str(OUT_ALL)  # full census goes to its own file by default
    global _SOLVER, _MAX_CARD, _MAX_BLOCKS
    _SOLVER, _MAX_CARD, _MAX_BLOCKS = args.solver, args.max_card, args.max_blocks

    if args.legacy_preflight:
        # historical calibration path: attain reaches the spectrum numerically
        # but lands on an arbitrary gauge, so exactify only certifies design
        # vertices (trivial phases), not interference vertices like v_B. Kept
        # for regression-checking attain itself; it is minutes-to-hours slow.
        spec = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
        rec = solve_vertex(4, 9, spec)
        print("TierA", rec["status"], "support", rec.get("support_size"),
              "residual", rec.get("residual"))
        if rec["status"] != "OK" or rec.get("residual", 1) > 1e-9:
            return 1
        ex = exactify(4, 9, spec, rec)
        print("TierB", ex["status"], ex.get("reason", ""))
        return 0

    if args.preflight:
        # the real gate: the weights-first solver must reconstruct v_B (Tier A,
        # exact) AND exactify to a certified closed form (Tier B). This is the
        # fast, reliable path and the one campaigns depend on; it completes in
        # minutes, unlike the legacy attain solver.
        spec = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
        rec = solve_vertex_exact_first(4, 9, spec, max_card=10,
                                       certify_tier_b=True)
        print("TierA", rec["status"], "support", rec.get("support_size"),
              "residual", rec.get("residual"))
        if rec["status"] != "OK" or rec.get("residual", 1) > 1e-12:
            print("preflight FAILED: v_B not reconstructed")
            return 1
        ex = rec.get("exact")
        print("TierB", ex["status"] if ex else "MISSING")
        if not ex or ex["status"] != "EXACT":
            print("preflight FAILED: v_B closed form not certified")
            return 1
        print("weights", ex["weights"], "/", ex["den"])
        for p in ex["pretty"]:
            print("  ", p)
        return 0

    out_path = pathlib.Path(args.out)
    ckpt = Checkpointer(out_path, interval=args.checkpoint_interval,
                        bucket=args.s3_bucket, key=args.s3_key, profile=args.s3_profile)
    ckpt.restore()  # a fresh spot instance resumes from the last S3 checkpoint

    done = set()
    if out_path.exists():
        for line in out_path.read_text().splitlines():
            try:
                r = json.loads(line)
                done.add((r["system"], r["index"]))
            except json.JSONDecodeError:
                pass

    todo = [t for t in tasks(args.systems.split(","), all_vertices=args.all)
            if (f"({t[0]},{t[1]})", t[2]) not in done]
    with out_path.open("a") as out:
        ckpt.attach(out)
        ckpt.install_signal_handlers()
        if args.workers <= 1:
            cache = {}
            for t in todo:
                n, d = t[0], t[1]
                if (n, d) not in cache:
                    cache[(n, d)] = _build(d, n)
                _emit(out, _work(t, cache[(n, d)], args.seed))
                ckpt.checkpoint()
        else:
            os.environ.setdefault("OMP_NUM_THREADS", "1")
            os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
            with mp.Pool(args.workers, initializer=_init_seed,
                         initargs=(args.seed, args.solver, args.max_card,
                                   args.max_blocks)) as pool:
                for rec in pool.imap_unordered(_work_solo, todo):
                    _emit(out, rec)
                    ckpt.checkpoint()
        ckpt.checkpoint(force=True)  # final flush + upload
    return 0


# solver selection, set in main() and inherited by workers (fork) or re-set in
# the pool initializer (spawn); defaults keep a bare import importable
_SOLVER = "cascade"
_MAX_CARD = 14
_MAX_BLOCKS = 2


def _init_seed(seed: int, solver="cascade", max_card=14, max_blocks=2) -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    np.random.seed(seed + os.getpid() % 10000)
    global _SOLVER, _MAX_CARD, _MAX_BLOCKS
    _SOLVER, _MAX_CARD, _MAX_BLOCKS = solver, max_card, max_blocks


def _weights_first(n, d, spec, built):
    return solve_vertex_exact_first(n, d, spec, max_card=_MAX_CARD,
                                    max_blocks=_MAX_BLOCKS,
                                    certify_tier_b=True, _built=built)


def _work(t, built, seed):
    n, d, i, v, verdict = t
    np.random.seed(seed + 1000 * n + 10 * d + i)
    spec = [Fraction(s) for s in v["spectrum"]]
    # route by the verdict we already computed: a DESIGN-INT vertex is built
    # directly from its design witness (one-hop-free support gives a diagonal
    # 1-RDM, exact by construction), no iterative solve. Everything else goes
    # to the state solver. cascade (default): try the fast block-ansatz solver,
    # which returns a certified closed form when the vertex is block-structured,
    # and fall back to attain on FAIL so we still record a Tier-A numeric state.
    # Neither solver covers every interference vertex; uncovered ones land in
    # Tier-A / TIER-C for the extended-ansatz or hand-analysis frontier.
    rec = None
    if verdict == "DESIGN-INT":
        rec = solve_design_vertex(n, d, spec)
    if rec is None:
        if _SOLVER == "attain":
            rec = solve_vertex(n, d, spec, _built=built)
        elif _SOLVER == "weights-first":
            rec = _weights_first(n, d, spec, built)
        else:  # cascade
            rec = _weights_first(n, d, spec, built)
            if rec.get("status") != "OK":
                rec = solve_vertex(n, d, spec, _built=built)
    rec["system"], rec["index"], rec["classified"] = f"({n},{d})", i, verdict
    rec["integer_form"], rec["denominator"] = v["integer_form"], v["denominator"]
    # ensure a Tier-B record even when the solver did not already certify one
    # (weights-first and the design builder attach rec["exact"] only lazily)
    if rec.get("status") == "OK" and "exact" not in rec:
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
