#!/usr/bin/env python3
"""Minimize every certified state's support over its NATURAL-ORBITAL GAUGE.

The census library records the support of ONE point of each state's gauge
family U(m_1) x ... x U(m_k), never minimized over that family. This driver
does the minimization and reports, per vertex:

  - the gauge-invariant lower bound (profile classes, refined by per-class
    flattening ranks), which certifies gauge-minimality outright when it is
    attained,
  - the best support the search reached, closed EXACTLY on the orbit rather
    than by truncation, so the 1-RDM stays diag(lambda) identically,
  - the acceptance gate: 1-RDM exactly diagonal and equal to lambda, unit norm,
    and the 2-RDM spectrum unchanged (a gauge rotation cannot move it).

Everything here bounds s_Q^NO, the natural-orbital minimal support. It says
nothing about the free-basis minimum: see scripts/orbit_check.py.

  python scripts/gauge_census.py [--jobs N] [--limit N]

Writes results/data/gauge_min.jsonl and results/data/gauge_min_summary.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import pathlib
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from gpc_census.fiber import amplitudes_from_record, one_rdm
from gpc_census.gauge import (
    degenerate_blocks,
    flattening_ranks,
    greedy_gauge_drop,
    profile_classes,
    profile_lower_bound,
    two_rdm_spectrum,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
DIAG_TOL = 1e-10
SPEC_TOL = 1e-9


def run_one(record, tries=6):
    c, dets = amplitudes_from_record(record)
    lam = [v / record["denominator"] for v in record["integer_form"]]
    d = len(lam)
    blocks = degenerate_blocks(lam)
    ranks = flattening_ranks(c, dets, blocks)
    bound = int(sum(ranks.values()))
    classes = profile_classes(dets, blocks)
    t0 = time.time()

    out = {
        "system": record["system"],
        "index": record["index"],
        "classified": record["classified"],
        "support_library": len(dets),
        "profile_class_count": profile_lower_bound(dets, blocks),
        "gauge_lower_bound": bound,
        "block_sizes": [len(b) for b in blocks],
        "gauge_dim": int(sum(len(b) ** 2 for b in blocks)),
        "profile_classes": {",".join(map(str, k)): len(v) for k, v in classes.items()},
        "class_flattening_ranks": {",".join(map(str, k)): v for k, v in ranks.items()},
        "support_no_min": len(dets),
        "dropped_by": 0,
        "minimal_certified": len(dets) == bound,
        "rounds": [],
        "gate": {},
        "evidence": "EXACT lower bound; upper bound from an exactly closed "
                    "gauge rotation (residual reported), not a truncation",
    }
    if len(dets) == bound:
        out["seconds"] = round(time.time() - t0, 2)
        out["note"] = "library support attains the gauge lower bound: minimal, certified"
        return out

    cf, df, rounds = greedy_gauge_drop(c, dets, blocks, d, tries=tries)
    out["rounds"] = rounds
    out["support_no_min"] = len(df)
    out["dropped_by"] = len(dets) - len(df)
    out["minimal_certified"] = len(df) == bound
    if len(df) < len(dets):
        rho = one_rdm(cf, df, d)
        diag_err = float(np.max(np.abs(rho - np.diag(lam))))
        spec_dev = float(np.max(np.abs(two_rdm_spectrum(cf, df, d)
                                       - two_rdm_spectrum(c, dets, d))))
        out["gate"] = {
            "one_rdm_diagonal_error": diag_err,
            "norm_error": float(abs(np.linalg.norm(cf) - 1.0)),
            "two_rdm_spectrum_deviation": spec_dev,
            "passes": bool(diag_err < DIAG_TOL and spec_dev < SPEC_TOL),
        }
        out["final_support_dets"] = [list(t) for t in df]
        out["final_amplitudes_real"] = [float(x) for x in np.real(cf)]
        out["final_amplitudes_imag"] = [float(x) for x in np.imag(cf)]
    out["seconds"] = round(time.time() - t0, 2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    records = [json.loads(line) for line
               in (ROOT / "results/data/states.jsonl").read_text().splitlines()]
    if args.limit:
        records = records[: args.limit]

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(run_one, records, chunksize=4))

    out = ROOT / "results/data/gauge_min.jsonl"
    out.write_text("".join(json.dumps(r) + "\n" for r in rows))

    dropped = [r for r in rows if r["dropped_by"] > 0]
    tally = collections.Counter(r["classified"] for r in dropped)
    summary = {
        "note": "Support minimized over the natural-orbital gauge U(m_1) x ... x "
                "U(m_k). Bounds s_Q^NO, the census quantity. NOT the free-basis "
                "or orbit-minimal support.",
        "states": len(rows),
        "certified_gauge_minimal": sum(1 for r in rows if r["minimal_certified"]),
        "certified_without_search": sum(1 for r in rows
                                        if r["support_library"] == r["gauge_lower_bound"]),
        "states_that_gauge_sparsify": len(dropped),
        "amplitudes_removed": sum(r["dropped_by"] for r in dropped),
        "sparsifying_by_class": dict(tally),
        "gate_failures": [f"{r['system']} v{r['index']}" for r in dropped
                          if not r["gate"].get("passes")],
        "slack_histogram": dict(sorted(collections.Counter(
            r["support_library"] - r["gauge_lower_bound"] for r in rows).items())),
        "seconds": round(time.time() - t0, 1),
    }
    (ROOT / "results/data/gauge_min_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
