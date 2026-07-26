#!/usr/bin/env python3
"""Census-wide sparsification cascade (fiber program, Task 1).

Runs gpc_census.cascade over every certified state in results/data/states.jsonl
and records, per vertex:

  - support_upper: the smallest support reached on the fixed-spectrum fiber.
    The shipped support column is the CERTIFIED support, an upper bound on the
    minimal support s_Q; support_upper is a better upper bound wherever the
    cascade landed a drop. Neither column is ever a minimality claim.
  - both fiber tangent dimensions, each under a fiber-qualified key.
  - the cancellation-regime flag: exactly diagonal 1-RDM with at least one
    one-hop channel, i.e. every channel silent by cancellation.

All results are NUMERICAL (residuals reported). Writes results/data/cascade.jsonl
and results/data/cascade_summary.json.

  python scripts/census_cascade.py [--jobs N] [--limit N] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import pathlib
import sys
import time
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from gpc_census.cascade import cascade, diagonality_error  # noqa: E402
from gpc_census.fiber import (  # noqa: E402
    amplitudes_from_record, fixed_rho_tangent, fixed_spectrum_tangent,
    incidence_kernel_dim, one_rdm,
)

DIAG_TOL = 1e-9


def system_rank(system):
    n, d = system.strip("()").split(",")
    return int(n), int(d)


def one_hop_channels(dets):
    """Orbital pairs (p, q) joined by a one-hop determinant pair of the support."""
    out = set()
    for i in range(len(dets)):
        for j in range(i + 1, len(dets)):
            sd = set(dets[i]) ^ set(dets[j])
            if len(sd) == 2:
                out.add(tuple(sorted(sd)))
    return out


def analyse(record):
    n, d = system_rank(record["system"])
    c, dets = amplitudes_from_record(record)
    den = int(record["denominator"])
    spectrum = [Fraction(int(x), den) for x in record["integer_form"]]

    rho = one_rdm(c, dets, d)
    off = float(np.max(np.abs(rho - np.diag(np.diag(rho))))) if d else 0.0
    channels = one_hop_channels(dets)
    silent = off < DIAG_TOL and bool(channels)

    t0 = time.time()
    no = cascade(c, dets, spectrum, d, locus="natural_orbital")
    free = cascade(c, dets, spectrum, d, locus="fixed_spectrum")
    elapsed = time.time() - t0

    rec = {
        "system": record["system"],
        "index": record["index"],
        "classified": record["classified"],
        "denominator": den,
        "support_certified": len(dets),
        # the library state's own distance from the attainment convention:
        # 0 means its 1-RDM is diagonal with diagonal equal to lambda
        "start_diagonality_error": diagonality_error(c, dets, spectrum, d),
        "start_is_natural_orbital": diagonality_error(c, dets, spectrum, d) < DIAG_TOL,

        # s_Q^NO: support over states with rho = diag(lambda). Comparable with
        # the diagonal-incidence baseline s_I.
        "support_upper_no": no["support_upper"],
        "support_upper_no_dropped_by": len(dets) - no["support_upper"],
        "support_upper_no_provenance": no["method"],
        "support_upper_no_evidence": "numerical",
        "support_upper_no_blocked_reason": no["blocked_reason"],
        "support_upper_no_diagonality_error": no["final_diagonality_error"],
        "support_upper_no_dets": no["final_support_dets"],
        "support_upper_no_amplitudes_real": no["final_amplitudes_real"],
        "support_upper_no_amplitudes_imag": no["final_amplitudes_imag"],
        "support_upper_no_rounds": no["rounds"],

        # s_Q^free: support over the whole U(d)-saturated set. Lower bound is
        # majorization ONLY; never compare these against s_I.
        "support_upper_free": free["support_upper"],
        "support_upper_free_dropped_by": len(dets) - free["support_upper"],
        "support_upper_free_provenance": free["method"],
        "support_upper_free_evidence": "numerical",
        "support_upper_free_blocked_reason": free["blocked_reason"],
        "support_upper_free_diagonality_error": free["final_diagonality_error"],
        "support_upper_free_dets": free["final_support_dets"],
        "support_upper_free_amplitudes_real": free["final_amplitudes_real"],
        "support_upper_free_amplitudes_imag": free["final_amplitudes_imag"],
        "support_upper_free_rounds": free["rounds"],
        "free_endpoint_is_natural_orbital":
            free["final_diagonality_error"] < DIAG_TOL,

        "final_spectrum_error": free["final_spectrum_error"],
        "multiplicity_drift_seen": (no["multiplicity_drift_seen"]
                                    or free["multiplicity_drift_seen"]),
        "incidence_kernel_dim": incidence_kernel_dim(dets),
        "one_hop_channels": len(channels),
        "rdm_max_offdiagonal": off,
        "cancellation_regime": silent,
        "seconds": round(elapsed, 3),
    }
    rec.update({k: v for k, v in fixed_spectrum_tangent(c, dets, d).as_json().items()
                if k != "fiber"})
    rec.update({k: v for k, v in fixed_rho_tangent(c, dets, d).as_json().items()
                if k != "fiber"})
    return rec


def _worker(line):
    record = json.loads(line)
    try:
        return analyse(record)
    except Exception as exc:  # a failure must not lose the other 798
        return {"system": record["system"], "index": record["index"],
                "error": f"{type(exc).__name__}: {exc}"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (mp.cpu_count() or 2)))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=str(ROOT / "results" / "data"))
    args = ap.parse_args()

    lines = [x for x in (ROOT / "results/data/states.jsonl").read_text().splitlines()
             if x.strip()]
    if args.limit:
        lines = lines[: args.limit]

    started = time.time()
    with mp.Pool(args.jobs) as pool:
        results = []
        for i, rec in enumerate(pool.imap_unordered(_worker, lines, chunksize=4), 1):
            results.append(rec)
            if i % 50 == 0:
                print(f"  {i}/{len(lines)}  {time.time() - started:.0f}s", flush=True)

    results.sort(key=lambda r: (r["system"], r["index"]))
    out_dir = pathlib.Path(args.out)
    with open(out_dir / "cascade.jsonl", "w") as fh:
        for r in results:
            fh.write(json.dumps(r) + "\n")

    summary = summarize(results, time.time() - started)
    (out_dir / "cascade_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def summarize(results, elapsed):
    ok = [r for r in results if "error" not in r]
    no_drop = [r for r in ok if r["support_upper_no_dropped_by"] > 0]
    free_drop = [r for r in ok if r["support_upper_free_dropped_by"] > 0]
    silent = [r for r in ok if r["cancellation_regime"]]
    on_locus = [r for r in ok if r["start_is_natural_orbital"]]
    free_ends_diag = [r for r in free_drop if r["free_endpoint_is_natural_orbital"]]

    def reasons(key):
        out = {}
        for r in ok:
            out[r[key]] = out.get(r[key], 0) + 1
        return dict(sorted(out.items()))

    return {
        "states": len(results),
        "errors": [r for r in results if "error" in r],

        "note": "Two DIFFERENT minimal supports are reported. s_Q^NO is the "
                "minimum over states whose 1-RDM is diagonal with diagonal "
                "equal to lambda (the census attainment convention; s_I is a "
                "valid lower bound). s_Q^free is the minimum over the whole "
                "U(d)-saturated spectrum locus (majorization is its only valid "
                "lower bound). s_Q^free <= s_Q^NO, and free-basis endpoints "
                "must never be compared against s_I.",

        "library_states_on_natural_orbital_locus": len(on_locus),
        "library_states_off_natural_orbital_locus": len(ok) - len(on_locus),
        "library_off_locus_by_class": _by([r for r in ok
                                           if not r["start_is_natural_orbital"]],
                                          "classified"),

        "s_Q_NO": {
            "states_with_support_drop": len(no_drop),
            "total_amplitudes_removed": sum(r["support_upper_no_dropped_by"]
                                            for r in no_drop),
            "largest_drop": max((r["support_upper_no_dropped_by"] for r in ok),
                                default=0),
            "worst_diagonality_error_among_drops": max(
                (r["support_upper_no_diagonality_error"] for r in no_drop),
                default=0.0),
            "blocked_reasons": reasons("support_upper_no_blocked_reason"),
            "drops_by_class": _by(no_drop, "classified"),
            "drops_by_system": _by(no_drop, "system"),
        },

        "s_Q_free": {
            "states_with_support_drop": len(free_drop),
            "total_amplitudes_removed": sum(r["support_upper_free_dropped_by"]
                                            for r in free_drop),
            "largest_drop": max((r["support_upper_free_dropped_by"] for r in ok),
                                default=0),
            "endpoints_that_happen_to_be_diagonal": len(free_ends_diag),
            "worst_diagonality_error_among_drops": max(
                (r["support_upper_free_diagonality_error"] for r in free_drop),
                default=0.0),
            "blocked_reasons": reasons("support_upper_free_blocked_reason"),
            "drops_by_class": _by(free_drop, "classified"),
            "drops_by_system": _by(free_drop, "system"),
        },

        "cancellation_regime_members": [
            {"system": r["system"], "index": r["index"],
             "support_certified": r["support_certified"],
             "support_upper_no": r["support_upper_no"],
             "support_upper_free": r["support_upper_free"],
             "fixed_spectrum_tangent_dim": r["fixed_spectrum_tangent_dim"],
             "fixed_rho_tangent_dim": r["fixed_rho_tangent_dim"]}
            for r in silent],
        "cancellation_regime_count": len(silent),
        "evidence": "numerical: continuation residuals, not certified solutions",
        "seconds": round(elapsed, 1),
    }


def _by(rows, key):
    out = {}
    for r in rows:
        out[r[key]] = out.get(r[key], 0) + 1
    return dict(sorted(out.items()))


if __name__ == "__main__":
    main()
