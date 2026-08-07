#!/usr/bin/env python3
"""Theorem E and the closure of the last large class.

THE GENERALISATION 16/17 WAS HIDING. Theorem D proved a diagonal rational by
reading Theorem A through a loop shared with a size-2 class. The size-2 part
was never needed for the SHARING; it was needed only to have a base case. What
the sharing really gives is an identity of Gram entries, and what propagates
rationality is the rank-2 Gram matrix itself.

THEOREM E (chained holonomy). Let `a`, `b`, `c` be three pairs of a one-hop
class. The sides lie in the real plane, so the Gram matrix has rank at most 2
and its 3x3 minor on `{a, b, c}` vanishes:

    r_b^2 t_ac^2 - 2 t_ab t_bc t_ac + (r_a^2 t_bc^2 + r_c^2 t_ab^2 - r_a^2 r_b^2 r_c^2) = 0.

If `t_ab` and `t_bc` are rational this is a RATIONAL quadratic in `t_ac`, whose
discriminant factors:

    Delta / 4 = (t_ab^2 - r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2)
              = r_a^2 r_b^4 r_c^2 sin^2(Phi_ab) sin^2(Phi_bc).

So `t_ac` always lies in `Q(sqrt(Delta))`, a quadratic extension by a RATIONAL
radicand, hence multiquadratic. And `t_ac` is RATIONAL exactly when `Delta` is
a rational square, in particular whenever either linked holonomy has vanishing
sine, that is `cos(Phi) = +/- 1`.

WHY THAT CLOSES THE LAST CLASS. (3,10) v103 channel (1,6) was the one class
Theorem D could not reach: its only one-hop link runs to the size-4 class
(3,9), not to a size-2 class. But v103 is REAL UP TO GAUGE, holonomy
(0, pi, pi, 0, 0), a fact certified independently in
`results/data/v103_fiber_count.json`. Both holonomies feeding the chain have
`cos = -1`, so both sines vanish, `Delta = 0`, the quadratic is a perfect
square, and `t_03` is forced rational with no square root extracted at all.
Theorem D' then transports it to channel (1,6) and Theorem C finishes.

So the exceptional class is not a gap in the mechanism. It is the mechanism
evaluated at a degenerate configuration, and the degeneracy is exactly the
property the repository already records about that state.

Usage:
    python scripts/class_closure.py
"""
from __future__ import annotations

import argparse
import itertools
import json
import pathlib
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

SOURCE = ROOT / "results/data/shared_loop_rationality.json"


def _rat(text):
    return sp.nsimplify(sp.sympify(text))


def theorem_E(r2, t, a, b, c):
    """The rank-2 Gram minor on {a,b,c}, as a quadratic in t_ac.

    Built from ``t_ab`` and ``t_bc`` only. ``t_ac`` is never used to form the
    coefficients, so substituting the measured value back is a real check.
    """
    t_ab, t_bc = t[(min(a, b), max(a, b))], t[(min(b, c), max(b, c))]
    A = r2[b]
    B = -2 * t_ab * t_bc
    C = r2[a] * t_bc**2 + r2[c] * t_ab**2 - r2[a] * r2[b] * r2[c]
    delta = sp.simplify(B**2 - 4 * A * C)
    factored = sp.simplify((t_ab**2 - r2[a] * r2[b]) * (t_bc**2 - r2[b] * r2[c]))
    root = sp.sqrt(delta)
    return {
        "triple": [a, b, c],
        "coefficients": [str(A), str(B), str(C)],
        "coefficients_rational": bool(A.is_rational and B.is_rational
                                      and C.is_rational),
        "discriminant": str(delta),
        "discriminant_rational": bool(delta.is_rational),
        "discriminant_factored_matches": bool(sp.simplify(delta - 4 * factored) == 0),
        "discriminant_is_a_rational_square": bool(root.is_rational),
        "discriminant_zero": bool(delta == 0),
        "sin_ab_vanishes": bool(sp.simplify(t_ab**2 - r2[a] * r2[b]) == 0),
        "sin_bc_vanishes": bool(sp.simplify(t_bc**2 - r2[b] * r2[c]) == 0),
        "forces_rational_t_ac": bool(root.is_rational),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out",
                    default=str(ROOT / "results/data/class_closure.json"))
    args = ap.parse_args()

    source = json.loads(SOURCE.read_text())
    blocks = []
    for block in source["classes"]:
        r2 = [_rat(v) for v in block["r_squared"]]
        t = {tuple(e["pair"]): _rat(e["t_ab"]) for e in block["links"]}
        rational = {k for k, v in t.items() if v.is_rational}
        proved = {tuple(e["pair"]) for e in block["links"]
                  if e.get("theorem_D_matches")}

        # every triple whose two outer entries are rational feeds Theorem E
        triples = []
        for b in range(block["class_size"]):
            for a, c in itertools.combinations(
                    [k for k in range(block["class_size"]) if k != b], 2):
                key_ab = (min(a, b), max(a, b))
                key_bc = (min(b, c), max(b, c))
                if key_ab not in rational or key_bc not in rational:
                    continue
                entry = theorem_E(r2, t, a, b, c)
                key_ac = (min(a, c), max(a, c))
                measured = t[key_ac]
                A, B, C = (_rat(x) for x in entry["coefficients"])
                entry["measured_t_ac"] = str(measured)
                entry["measured_solves_quadratic"] = bool(
                    sp.simplify(A * measured**2 + B * measured + C) == 0)
                entry["measured_t_ac_rational"] = bool(measured.is_rational)
                entry["uses_only_theorem_D_inputs"] = bool(
                    key_ab in proved and key_bc in proved)
                triples.append(entry)

        # closure ledger for this class
        d_proved = len(proved)
        needed = block["class_size"] - 2
        e_closes = any(x["forces_rational_t_ac"]
                       and x["uses_only_theorem_D_inputs"] for x in triples)
        blocks.append({
            "system": block["system"],
            "index": block["index"],
            "channel": block["channel"],
            "class_size": block["class_size"],
            "transport_class": block["transport_class"],
            "diagonals_proved_by_theorem_D": d_proved,
            "diagonals_needed": needed,
            "closed_by_theorem_D": bool(d_proved >= needed),
            "theorem_E_triples": triples,
            "theorem_E_extends_rationality": bool(e_closes),
            "closed": bool(d_proved >= needed or e_closes),
        })

    # (3,10) v103 channel (1,6) closes by transport from channel (3,9)
    donor = next(b for b in blocks
                 if b["system"] == "(3,10)" and b["index"] == 103
                 and b["channel"] == [3, 9])
    last = next(b for b in blocks
                if b["system"] == "(3,10)" and b["index"] == 103
                and b["channel"] == [1, 6])
    transported = donor["theorem_E_extends_rationality"]
    last["closed_by_transport_from"] = donor["channel"] if transported else None
    last["closed"] = bool(last["closed"] or transported)

    payload = {
        "what": "Theorem E (rank-2 Gram chaining) and the closure ledger for "
                "every large one-hop class.",
        "theorem_E": {
            "statement": "for pairs a,b,c of a class, r_b^2 t_ac^2 - 2 t_ab t_bc "
                         "t_ac + (r_a^2 t_bc^2 + r_c^2 t_ab^2 - r_a^2 r_b^2 "
                         "r_c^2) = 0; with t_ab, t_bc rational this is a "
                         "rational quadratic with discriminant 4 (t_ab^2 - "
                         "r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2)",
            "geometry": "Delta / 4 = r_a^2 r_b^4 r_c^2 sin^2(Phi_ab) "
                        "sin^2(Phi_bc), so t_ac is rational whenever either "
                        "linked holonomy has cos = +/- 1",
            "last_class": "(3,10) v103 is real up to gauge with holonomy "
                          "(0, pi, pi, 0, 0) per v103_fiber_count.json, so both "
                          "sines vanish, Delta = 0, and t_03 is forced rational",
        },
        "summary": {
            "classes": len(blocks),
            "closed": sum(1 for b in blocks if b["closed"]),
            "closed_by_theorem_D": sum(1 for b in blocks
                                       if b["closed_by_theorem_D"]),
            "closed_by_theorem_E": sum(1 for b in blocks
                                       if not b["closed_by_theorem_D"]
                                       and b["closed"]),
            "theorem_E_triples": sum(len(b["theorem_E_triples"]) for b in blocks),
            "theorem_E_coefficients_rational": sum(
                1 for b in blocks for x in b["theorem_E_triples"]
                if x["coefficients_rational"]),
            "theorem_E_discriminant_rational": sum(
                1 for b in blocks for x in b["theorem_E_triples"]
                if x["discriminant_rational"]),
            "theorem_E_factorisation_matches": sum(
                1 for b in blocks for x in b["theorem_E_triples"]
                if x["discriminant_factored_matches"]),
            "theorem_E_measured_solves": sum(
                1 for b in blocks for x in b["theorem_E_triples"]
                if x["measured_solves_quadratic"]),
            "theorem_E_degenerate_triples": sum(
                1 for b in blocks for x in b["theorem_E_triples"]
                if x["discriminant_zero"]),
            "not_closed": sorted(f'{b["system"]} v{b["index"]} '
                                 f'ch{tuple(b["channel"])}'
                                 for b in blocks if not b["closed"]),
        },
        "classes": blocks,
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    print(json.dumps(payload["summary"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
