#!/usr/bin/env python3
"""Reconciliation duty: re-derive the four claims this session corrects.

A parallel instance works this repository and may be citing the superseded
readings. Every item below is checked against the shipped ledger and the
current code, so the correction points at an artifact rather than at prose
(rule 2 of the fiber program). Writes results/data/reconciliation.json.

  python scripts/reconcile_claims.py
"""
from __future__ import annotations

import json
import pathlib
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import sympy as sp  # noqa: E402

from gpc_census.cascade import cascade  # noqa: E402
from gpc_census.fiber import (  # noqa: E402
    amplitudes_from_record, fixed_rho_tangent, fixed_spectrum_tangent,
)

LEDGER = ROOT / "results/data/states.jsonl"


def rows():
    return [json.loads(x) for x in LEDGER.read_text().splitlines() if x.strip()]


def get(all_rows, system, index):
    return next(r for r in all_rows if r["system"] == system and r["index"] == index)


def spectrum_of(rec):
    den = int(rec["denominator"])
    return [Fraction(int(x), den) for x in rec["integer_form"]]


def item_p4_rigidity(all_rows):
    """P4 scorecard: v103 'first order rigid' is a FIXED-RHO verdict only."""
    out = {}
    for idx in (89, 103):
        rec = get(all_rows, "(3,10)", idx)
        c, dets = amplitudes_from_record(rec)
        out[f"v{idx}"] = {
            "fixed_rho_tangent_dim": fixed_rho_tangent(c, dets, 10).tangent_dim,
            "fixed_spectrum_tangent_dim": fixed_spectrum_tangent(c, dets, 10).tangent_dim,
        }
    v103, v89 = out["v103"], out["v89"]
    return {
        "claim_superseded": "P4 scorecard reads v103 as first-order rigid without "
                            "naming the fiber",
        "correction": "v103 is rigid for the fixed-rho fiber and DEFORMABLE for "
                      "the fixed-spectrum fiber; v89 is rigid in both senses",
        "measured": out,
        "holds": (v103["fixed_rho_tangent_dim"] == 0
                  and v103["fixed_spectrum_tangent_dim"] > 0
                  and v89["fixed_rho_tangent_dim"] == 0
                  and v89["fixed_spectrum_tangent_dim"] == 0),
        "evidence": "numerical (first-order rank computation)",
    }


def item_isolated_support_14(all_rows):
    """v103's 'isolated support-14 point' is a landscape, not a fiber, statement."""
    rec = get(all_rows, "(3,10)", 103)
    c, dets = amplitudes_from_record(rec)
    out = cascade(c, dets, spectrum_of(rec), 10)
    return {
        "claim_superseded": "v103's certified state is an isolated support-14 "
                            "point, dropping any amplitude breaks feasibility",
        "correction": "true of the route-B search landscape and gauge slice, "
                      "false of the fixed-spectrum fiber",
        "measured": {
            "support_certified": out["support_start"],
            "support_upper": out["support_upper"],
            "dropped_dets": out["dropped_dets"],
            "final_certificate_residual": out["final_residual"],
            "final_spectrum_error": out["final_spectrum_error"],
            "blocked_reason": out["blocked_reason"],
            "floor_residual": out["floor_residual"],
        },
        "holds": out["support_upper"] < out["support_start"],
        "evidence": "numerical (minimal-polynomial continuation)",
    }


def item_support_column(all_rows):
    """The ledger support column is an UPPER bound on s_Q, strict at v103."""
    rec = get(all_rows, "(3,10)", 103)
    c, dets = amplitudes_from_record(rec)
    out = cascade(c, dets, spectrum_of(rec), 10)
    return {
        "claim_superseded": "the census support column is the minimal support",
        "correction": "it is the CERTIFIED support, an upper bound on the minimal "
                      "support s_Q, and it is known strict at v103",
        "measured": {"v103_certified_support": out["support_start"],
                     "v103_support_upper": out["support_upper"]},
        "holds": out["support_upper"] < out["support_start"],
        "evidence": "numerical",
    }


def item_design_real_ratio(all_rows):
    """Mined law '9 of 12 DESIGN-REAL at denominator ratio 2': ledger has 12/12."""
    ratios = {}
    for rec in all_rows:
        if rec["classified"] != "DESIGN-REAL":
            continue
        den = int(rec["denominator"])
        squares = [sp.nsimplify(sp.sympify(s)) ** 2 for s in rec["closed_form"]["pretty"]]
        state_den = sp.ilcm(*[sp.Rational(x).q for x in squares])
        ratios[f"{rec['system']} v{rec['index']}"] = str(sp.Rational(state_den, den))
    at_two = sum(1 for v in ratios.values() if v == "2")
    return {
        "claim_superseded": "9 of 12 DESIGN-REAL states sit at denominator ratio 2",
        "correction": f"the SHIPPED ledger has {at_two} of {len(ratios)} at ratio 2; "
                      "the mined count is representative dependent, exactly the "
                      "P2/P3/P4 caveat",
        "measured": {"ratios": ratios, "at_ratio_2": at_two, "total": len(ratios)},
        "holds": at_two == len(ratios) == 12,
        "evidence": "exact (rational arithmetic on the shipped closed forms)",
    }


def main():
    all_rows = rows()
    report = {
        "note": "Reconciliation of claims this session corrects. Each item names "
                "the superseded reading, the correction, and the measurement.",
        "items": {
            "p4_scorecard_v103_rigidity": item_p4_rigidity(all_rows),
            "v103_isolated_support_14": item_isolated_support_14(all_rows),
            "ledger_support_column": item_support_column(all_rows),
            "design_real_denominator_ratio": item_design_real_ratio(all_rows),
        },
    }
    report["all_hold"] = all(v["holds"] for v in report["items"].values())
    (ROOT / "results/data/reconciliation.json").write_text(
        json.dumps(report, indent=2) + "\n")
    for name, item in report["items"].items():
        print(f"{name}: {'CONFIRMED' if item['holds'] else 'NOT CONFIRMED'}")
        print(f"    {item['correction']}")
    print(f"all_hold: {report['all_hold']}")


if __name__ == "__main__":
    main()
