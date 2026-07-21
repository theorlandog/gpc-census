"""Build / resume the certified closed-form states dataset (states.jsonl).

Checkpointed and restartable: the --out file IS the checkpoint. Every solved
vertex is merged into it and the file is rewritten periodically, so you can stop
(Ctrl-C) and resume anytime by re-running the same command. Parallel across
--workers; each vertex is solved in a child process with a hard per-vertex
timeout so nothing can stall the run.

Resume semantics:
  default            : solve only vertices not yet in the checkpoint (adds the
                       unprocessed systems, keeps everything already done).
  --retry-uncertified: also re-attempt vertices recorded as SOLVE-FAIL /
                       NO-EXACT / TIMEOUT (e.g. with a bigger --clique-timeout on
                       stronger hardware); already-certified vertices are kept.

Records match results/data/states.jsonl and gpc_census.dataset.validate_states,
so validate before shipping:
    uv run scripts/validate_states.py <out>

Examples:
    # resume my dataset on your rig, re-attacking the hard tail with more time
    uv run scripts/build_states.py --out results/data/states.jsonl \\
        --retry-uncertified --clique-timeout 600 --workers 8
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import pathlib
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

DATA = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
SYSTEMS = ["3_9", "4_9", "3_10", "4_10", "5_10"]
CERTIFIED = {"EXACT", "EXACT-CONSTR"}


def _worker(n, d, spec_str, verdict, max_cliques, clique_budget, q, max_clique=3,
            max_card=16):
    import os
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    from gpc_census.states import (solve_design_real_vertex, solve_design_vertex,
                                   solve_vertex_exact_first)
    spec = [Fraction(s) for s in spec_str]
    if verdict == "DESIGN-INT":
        rec = solve_design_vertex(n, d, spec)
        if not rec or rec.get("status") != "OK":
            q.put(json.dumps({"status": "FAIL", "tierB": "DESIGN-FAIL"}))
            return
        sup = rec["support"]
        cf = {"den": rec["den"], "weights": rec["weights"],
              "pretty": [f"sqrt({w}/{rec['den']})" for w in rec["weights"]],
              "support_dets": [s[0] for s in sup]}
        q.put(json.dumps({"status": "OK", "tierB": "EXACT-CONSTR",
                          "support": sup, "closed_form": cf}))
        return
    if verdict == "DESIGN-REAL":
        rec = solve_design_real_vertex(n, d, spec)
        if rec and rec.get("status") == "OK":
            q.put(json.dumps({"status": "OK", "tierB": "EXACT",
                              "support": rec["support"],
                              "closed_form": rec["closed_form"]}))
            return
        # fall through to the general solver if no real design is found
    rec = solve_vertex_exact_first(n, d, spec, max_card=max_card,
                                   max_clique=max_clique,
                                   max_cliques=max_cliques,
                                   clique_time_budget=clique_budget,
                                   certify_tier_b=True)
    if not rec or rec.get("status") != "OK":
        q.put(json.dumps({"status": "FAIL", "tierB": "SOLVE-FAIL"}))
        return
    ex = rec.get("exact")
    tb = ex.get("status") if ex else "NO-EXACT"
    sup = rec.get("support")
    cf = None
    if ex and ex.get("status") == "EXACT":
        cf = {"den": ex["den"], "weights": ex["weights"], "pretty": ex["pretty"],
              "support_dets": [s[0] for s in sup]}
    q.put(json.dumps({"status": "OK", "tierB": tb, "support": sup, "closed_form": cf}))


def solve_one(task, max_cliques, clique_budget, cap, max_clique=3, max_card=16):
    n, d, i, verdict, spec_str = task
    q = mp.Queue()
    p = mp.Process(target=_worker,
                   args=(n, d, spec_str, verdict, max_cliques, clique_budget, q,
                         max_clique, max_card))
    t0 = time.time()
    p.start()
    p.join(cap)
    if p.is_alive():
        p.terminate()
        p.join(5)
        if p.is_alive():
            p.kill()
        res = {"status": "TIMEOUT", "tierB": "TIMEOUT"}
    else:
        try:
            res = json.loads(q.get_nowait())
        except Exception:
            res = {"status": "CRASH", "tierB": "CRASH"}
    return (n, d, i, verdict, spec_str, res, round(time.time() - t0, 1))


def _verdict(line):
    for v in ("DESIGN-INT", "DESIGN-REAL", "INTERFERENCE"):
        if v in line:
            return v
    return "UNKNOWN"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DATA / "states.jsonl"),
                    help="checkpoint file (read for resume, rewritten as it runs)")
    ap.add_argument("--systems", default=",".join(SYSTEMS))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max-cliques", type=int, default=1,
                    help="0 uses per-vertex capacity (full block search, slower)")
    ap.add_argument("--max-clique", type=int, default=3,
                    help="largest clique size (distinct eigenvalue classes mixed "
                         "in one block); 4 closes many k=3 SOLVE-FAILs, slower")
    ap.add_argument("--clique-timeout", type=float, default=300.0,
                    help="per clique-count-level budget; raise to certify slow vertices")
    ap.add_argument("--max-card", type=int, default=16,
                    help="max support cardinality (determinants); raise above the "
                         "vertex denominator to reach high-denominator supports "
                         "(e.g. (3,10) v89 den 26, v103 den 34)")
    ap.add_argument("--cap", type=float, default=0.0,
                    help="hard per-vertex kill seconds (default: 4x clique-timeout)")
    ap.add_argument("--retry-uncertified", action="store_true",
                    help="re-attempt SOLVE-FAIL / NO-EXACT / TIMEOUT vertices")
    args = ap.parse_args()
    out = pathlib.Path(args.out)
    cap = args.cap or (4 * args.clique_timeout)

    records = {}
    if out.exists():
        for ln in out.read_text().splitlines():
            if not ln.strip():
                continue
            try:
                r = json.loads(ln)
                records[(r["system"], r["index"])] = r
            except Exception:
                pass
    done = {k for k, r in records.items()
            if r.get("tierB") in CERTIFIED
            or not args.retry_uncertified}
    ncert = sum(1 for r in records.values() if r.get("tierB") in CERTIFIED)
    print(f"checkpoint {out}: {len(records)} records, {ncert} certified", flush=True)

    todo = []
    for tag in args.systems.split(","):
        n, d = (int(x) for x in tag.split("_"))
        verts = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
        rows = [ln for ln in (DATA / "census" / f"census_{tag}_results.txt")
                .read_text().splitlines() if ln[:4].strip().isdigit()]
        for i, line in enumerate(rows):
            if (f"({n},{d})", i) in done:
                continue
            todo.append((n, d, i, _verdict(line),
                         [str(s) for s in verts[i]["spectrum"]]))
    print(f"{len(todo)} vertices to solve on {args.workers} workers "
          f"(cap {cap:.0f}s/vertex)", flush=True)

    def flush():
        tmp = out.with_suffix(".tmp")
        tmp.write_text("\n".join(json.dumps(records[k]) for k in sorted(records)) + "\n")
        tmp.replace(out)  # atomic checkpoint

    n_new = 0
    last = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(solve_one, t, args.max_cliques, args.clique_timeout, cap,
                          args.max_clique, args.max_card)
                for t in todo]
        for fut in futs:
            n, d, i, verdict, spec_str, res, secs = fut.result()
            verts = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
            rec = {"system": f"({n},{d})", "index": i, "classified": verdict,
                   "integer_form": verts[i]["integer_form"],
                   "denominator": verts[i]["denominator"], "secs": secs}
            rec.update(res)
            records[(f"({n},{d})", i)] = rec
            n_new += 1
            tb = res.get("tierB")
            mark = "ok" if tb in CERTIFIED else f"GAP {tb}"
            print(f"  [{n_new}/{len(todo)}] ({n},{d}) idx {i:3d} {verdict:12s} "
                  f"{mark} [{secs}s]", flush=True)
            if time.time() - last > 30:   # checkpoint every 30s
                flush()
                last = time.time()
    flush()
    cert = sum(1 for r in records.values() if r.get("tierB") in CERTIFIED)
    print(f"done: {cert}/{len(records)} certified in {out}", flush=True)
    print("validate before shipping: "
          f"uv run scripts/validate_states.py {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
