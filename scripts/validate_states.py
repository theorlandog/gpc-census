"""Validate a states.jsonl atlas before it is committed.

Thin CLI over gpc_census.dataset.validate_states: every certified record's
closed form is re-checked independently of the solver that produced it, so the
shipped dataset is self-certifying. Exits nonzero if any certified record fails.

    uv run scripts/validate_states.py results/data/states.jsonl
"""
from __future__ import annotations

import json
import sys

from gpc_census.dataset import validate_states


def main(path: str) -> int:
    records = [json.loads(ln) for ln in open(path) if ln.strip()]
    n_cert = sum(1 for r in records
                 if r.get("tierB") in ("EXACT", "EXACT-CONSTR") and r.get("closed_form"))
    fails = validate_states(records)
    for sysname, idx, tb, why in fails:
        print(f"  FAIL {sysname} idx {idx} {tb}: {why}", flush=True)
    print(f"checked {n_cert} certified states; {len(fails)} failed")
    if fails:
        print("VALIDATION FAILED")
        return 1
    print("ALL CERTIFIED RECORDS VALID")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "results/data/states.jsonl"))
