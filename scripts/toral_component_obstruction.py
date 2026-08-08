"""The component group of the diagonal stabilizer, as a quantization obstruction.

`quantization_characters.json` searched only lambda-preserving PERMUTATIONS, so
it missed the obstruction carried by the DISCONNECTED part of the torus itself.
This script supplies that term by duality, without solving for any phase, and
corrects the census-wide predicted periods upward where it fires.

Scope is the 645 states whose 1-RDM is diagonal; the derivation needs
`lambda_o = sum_{S containing o} |psi_S|^2`. The other 154 are reported as
out of scope, not guessed at.

Run: uv run scripts/toral_component_obstruction.py
Writes results/data/toral_component_obstruction.json.
"""
from __future__ import annotations

import json
import pathlib
import sys
import time
from math import gcd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data"
OUT = DATA / "toral_component_obstruction.json"

from gpc_census.symplectic import toral_component_period  # noqa: E402


def main() -> int:
    states = [json.loads(x) for x in (DATA / "states.jsonl").read_text().splitlines()]
    chars = json.loads((DATA / "quantization_characters.json").read_text())
    perm = {(r["system"], r["index"]): r for r in chars["census_predictions"]}
    inv = json.loads((DATA / "symplectic_invariants.json").read_text())
    diagonal = {(r["system"], r["index"]): r["stabilizer"]["rho_is_diagonal"]
                for r in inv["rows"]}

    rows, started = [], time.time()
    for rec in states:
        key = (rec["system"], rec["index"])
        q = int(rec["denominator"] or 1)
        toral = toral_component_period(rec)
        before = perm[key]["predicted_period"]
        combined = before if toral is None else before * toral // gcd(before, toral)
        rows.append({
            "system": rec["system"], "index": rec["index"],
            "integer_form": list(rec["integer_form"]), "denominator": q,
            "classified": rec["classified"], "evidence": "exact",
            "rho_is_diagonal": diagonal[key],
            "in_scope": toral is not None,
            "toral_component_period": toral,
            "permutation_predicted_period": before,
            "combined_predicted_period": combined,
            "raises_prediction": bool(combined != before),
        })

    scope = [r for r in rows if r["in_scope"]]
    fires = [r for r in rows if r["in_scope"] and r["toral_component_period"] != r["denominator"]]
    raised = [r for r in rows if r["raises_prediction"]]
    payload = {
        "generated_by": "scripts/toral_component_obstruction.py",
        "corrects": "results/data/quantization_characters.json",
        "evidence": "exact",
        "note": (
            "Obstruction from the component group of the diagonal stabilizer, by "
            "Pontryagin duality: v_k = (-k lambda, k) must lie in the INTEGER "
            "rowspan of B, whose rows are [A_S | -1]. The gap between that lattice "
            "and its saturation IS the component group. Valid only where the 1-RDM "
            "is diagonal; the 154 other states are marked out of scope."
        ),
        "rows": rows,
        "summary": {
            "in_scope": len(scope),
            "out_of_scope": len(rows) - len(scope),
            "out_of_scope_all_non_diagonal": all(
                not r["rho_is_diagonal"] for r in rows if not r["in_scope"]),
            "obstruction_fires": len(fires),
            "fires_all_design_real": all(r["classified"] == "DESIGN-REAL" for r in fires),
            "predictions_raised": len(raised),
            "raised": [f"{r['system']}#{r['index']}" for r in raised],
            "ratios_seen": sorted({r["toral_component_period"] // r["denominator"]
                                   for r in scope}),
        },
        "wall_secs": round(time.time() - started, 1),
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(json.dumps(payload["summary"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
