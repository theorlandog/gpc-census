"""Close the 154 out-of-scope states by rewriting their support exactly.

`toral_component_obstruction.json` left 154 states out of scope: their 1-RDM is
not diagonal, so the standard-basis support incidence matrix describes the
wrong Levi. This rewrites each support in an exact eigenbasis of `rho` and
computes the obstruction there, over two eigenbasis variants whose lcm is a
tighter valid bound than either alone.

Run: uv run scripts/natural_orbital_toral.py
Writes results/data/natural_orbital_toral.json.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data"
OUT = DATA / "natural_orbital_toral.json"
CK = ROOT / "build" / "natural_orbital_toral.jsonl"

from gpc_census.symplectic import natural_orbital_toral_period  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip", nargs="*", default=[],
                    help="SYSTEM#INDEX to leave unattempted; named as unfinished")
    args = ap.parse_args()
    skip = set(args.skip)
    prev = json.loads((DATA / "toral_component_obstruction.json").read_text())
    targets = [r for r in prev["rows"] if not r["in_scope"]]
    states = {(x["system"], x["index"]): x
              for x in (json.loads(line)
                        for line in (DATA / "states.jsonl").read_text().splitlines())}
    CK.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if CK.exists():
        for line in CK.read_text().splitlines():
            row = json.loads(line)
            done[(row["system"], row["index"])] = row

    rows, started = [], time.time()
    with CK.open("a") as ck:
        for t in targets:
            key = (t["system"], t["index"])
            if f"{key[0]}#{key[1]}" in skip:
                continue
            if key in done:
                rows.append(done[key])
                continue
            t0 = time.time()
            res = natural_orbital_toral_period(states[key])
            row = {
                "system": key[0], "index": key[1], "denominator": t["denominator"],
                "integer_form": t["integer_form"], "classified": t["classified"],
                "evidence": "exact",
                "resolved": res is not None,
                "detail": res,
                "toral_period": (res or {}).get("combined"),
                "permutation_predicted_period": t["permutation_predicted_period"],
                "secs": round(time.time() - t0, 1),
            }
            from math import gcd
            a = row["permutation_predicted_period"]
            b = row["toral_period"]
            row["combined_predicted_period"] = a if not b else a * b // gcd(a, b)
            row["raises_prediction"] = row["combined_predicted_period"] != a
            rows.append(row)
            ck.write(json.dumps(row) + "\n")
            ck.flush()
            print(f"  {key[0]}#{key[1]} resolved={row['resolved']} "
                  f"toral={row['toral_period']} q={t['denominator']} "
                  f"raises={row['raises_prediction']} ({row['secs']}s)", flush=True)

    resolved = [r for r in rows if r["resolved"]]
    fires = [r for r in resolved if r["toral_period"] != r["denominator"]]
    raised = [r for r in rows if r["raises_prediction"]]
    payload = {
        "generated_by": "scripts/natural_orbital_toral.py",
        "extends": "results/data/toral_component_obstruction.json",
        "evidence": "exact",
        "note": (
            "The 154 states whose 1-RDM is not diagonal, with their supports "
            "rewritten in an exact eigenbasis of rho. Two eigenbasis variants are "
            "computed because the support is basis dependent inside a degenerate "
            "block; each gives a valid lower bound in the divisibility order and "
            "the lcm is reported."
        ),
        "rows": rows,
        "unfinished": sorted(skip),
        "summary": {
            "attempted": len(rows), "resolved": len(resolved),
            "unresolved": len(rows) - len(resolved),
            "obstruction_fires": len(fires),
            "fires": [f"{r['system']}#{r['index']}" for r in fires],
            "predictions_raised": len(raised),
            "all_interference": all(r["classified"] == "INTERFERENCE" for r in rows),
            "ratios_seen": sorted({r["toral_period"] // r["denominator"]
                                   for r in resolved}),
        },
        "wall_secs": round(time.time() - started, 1),
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(json.dumps(payload["summary"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
