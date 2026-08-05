#!/usr/bin/env python3
"""Blind first-principles reproduction of the fixed-N=3 systems: gates G2, G3.

The generator plan's validation ladder is

    G1  (3,6)             first-principles, sandwich-closed   [generate_3_6]
    G2  (3,7), (3,8)      blind, must match the A&K tables
    G3  (3,9), (3,10)     blind, must match the census H-data
    G4  (3,11)            the bracket
    G5  (3,12)+           frontier

and no new-rank output is worth anything until the gates below it close. This
script runs the coverage-first reference generator
(`generate_fixed_n3_reference_rank`) at the requested ranks and records, per
rank, the enumeration size, the screening partition, the retained system, and
the comparison against the stored table.

BLINDNESS. The generation path reads no stored constraints and no stored
vertices: it enumerates admissible affine flats, orients them into candidate
taus, screens each by the exact Ressayre criterion, and reduces on the ordered
Pauli slice. The stored table is consulted only afterwards, inside
`_known_system_regression`, to compare. A rank whose status is
``finalized_known_system_regression`` therefore reproduced the published
system without having been shown it.

COST. The flat enumeration dominates and grows steeply: 5341 admissible
hyperplanes at d = 7, 166420 at d = 8. Screening is cheap by comparison
(about 0.14 ms per tau). Resource exhaustion at a larger rank is an
operational result and carries no mathematical meaning, per the methodology
document's labelling rules.

Usage:
    python scripts/stage1_reference_series.py --ranks 7 8
    python scripts/stage1_reference_series.py --ranks 9 --timeout 0
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from gpc_census.generation import generate_fixed_n3_reference_rank  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results" / "data" / "stage1_reference_series.json"


def _counts(obj) -> dict:
    return {k: getattr(obj, k) for k in dir(obj)
            if not k.startswith("_") and isinstance(getattr(obj, k, None), int)}


def run_rank(d: int, workers: int, checkpoint_dir: str | None = None) -> dict:
    started = time.time()
    result = generate_fixed_n3_reference_rank(
        d, workers=workers, checkpoint_dir=checkpoint_dir)
    elapsed = time.time() - started
    cs = result.candidate_system
    reg = cs.known_system_regression
    rows = [
        {"coeffs": list(r.constraint.coeffs), "rhs": int(r.constraint.rhs)}
        for r in cs.retained_rows
    ]
    return {
        "d": d,
        "seconds": round(elapsed, 1),
        "admissible_hyperplanes": len(result.enumeration.hyperplanes),
        "oriented_taus": len(result.taus),
        "flat_counts_by_rank": list(result.enumeration.flat_counts_by_rank),
        "rank_preserving_modulus": result.enumeration.modulus,
        "hadamard_bound_squared": result.enumeration.hadamard_bound_squared,
        "screening_counts": _counts(result.screening.counts),
        "system_counts": _counts(cs.counts),
        "status": cs.status,
        "proved_complete": bool(result.proved_complete),
        "retained_rows": rows,
        "retained_count": len(rows),
        "known_system_regression": (reg.as_dict() if reg is not None else None),
        "gate_passed": bool(reg is not None and reg.passed),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranks", type=int, nargs="+", default=[7, 8])
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("-o", "--out", default=str(DEFAULT_OUT))
    ap.add_argument("--checkpoint-dir", default=None,
                    help="resume the flat enumeration across runs; a rank whose "
                         "cost is hours can then be finished incrementally")
    args = ap.parse_args()

    out_path = pathlib.Path(args.out)
    blocks: dict[str, dict] = {}
    if out_path.exists():          # accumulate across invocations
        blocks = json.loads(out_path.read_text()).get("ranks", {})

    for d in args.ranks:
        print(f"rank {d}: generating (workers={args.workers}) ...", flush=True)
        try:
            rec = run_rank(d, args.workers, args.checkpoint_dir)
        except Exception as exc:                     # operational, not a refutation
            rec = {"d": d, "status": "operational_failure",
                   "error": f"{type(exc).__name__}: {exc}", "gate_passed": False}
            print(f"  operational failure: {rec['error']}", flush=True)
        else:
            print(f"  {rec['seconds']}s  hyperplanes {rec['admissible_hyperplanes']}  "
                  f"retained {rec['retained_count']}  status {rec['status']}  "
                  f"gate_passed {rec['gate_passed']}", flush=True)
        blocks[str(d)] = rec

    payload = {
        "what": "Blind first-principles reproduction of the fixed-N=3 GPC "
                "systems, generator-plan gates G2 and G3.",
        "blindness": "the generation path reads no stored constraints or "
                     "vertices; the stored table is consulted only afterwards "
                     "for the regression comparison",
        "gate_ladder": {
            "G1": "(3,6) first-principles and sandwich-closed, generate_3_6",
            "G2": "(3,7), (3,8) blind against the Altunbulak-Klyachko tables",
            "G3": "(3,9), (3,10) blind against the census H-data",
        },
        "ranks": dict(sorted(blocks.items(), key=lambda kv: int(kv[0]))),
        "gates_passed": sorted(k for k, v in blocks.items() if v.get("gate_passed")),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
