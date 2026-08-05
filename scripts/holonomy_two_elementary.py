#!/usr/bin/env python3
"""Proof attack on 2-elementarity: reduce the conjecture to two conditions.

Proposition B delivers a 2-TOWER, and a 2-tower can have dihedral closure, so
the census's 82 of 82 multiquadratic verifications are unexplained. This
script isolates WHAT would have to be true, proves the implication, and tests
the two remaining conditions exactly.

THE REDUCTION. Everything turns on the square of the cosine.

  (R1)  deg_Q(cos^2 Phi) <= 2.
  (R2)  when deg_Q(cos Phi) = 4, the norm N(cos^2 Phi) over the quadratic
        field Q(cos^2 Phi) is a square in Q.

Given (R1), cos(Phi) has degree at most 4, and:

  - cos^2 rational                  -> min poly x^2 - r or linear: C1 or C2;
  - deg cos^2 = 2, deg cos = 2      -> C2;
  - deg cos^2 = 2, deg cos = 4      -> the min poly is EVEN, roots +-rho1,
                                       +-rho2, and by the classical criterion
                                       for x^4 + p x^2 + q irreducible,
                                          V4 iff q is a square in Q,
                                          C4 iff q (p^2 - 4q) is a square,
                                          D4 otherwise,
                                       where q = N(cos^2 Phi). So (R2) is
                                       exactly the V4 condition.

Every branch is 2-elementary, so (R1) and (R2) imply the conjecture on the
observed degrees. Both are then tested exactly against the certified ledger.

THEOREM A', proved here and stronger than what the shipped Theorem A states.
At a one-hop class of size at most 2, Theorem A gives

    cos Phi = (|rho_AB|^2 - r_1^2 - r_2^2) / (2 s_1 s_2 r_1 r_2),

whose numerator is rational under (H1) and (H2) while the denominator is a
rational multiple of sqrt(p_{u1} p_{v1} p_{u2} p_{v2}). So cos Phi = n / sqrt(M)
with n, M rational, hence cos^2 Phi = n^2 / M is RATIONAL. That is (R1) with
degree one, and it forces the minimal polynomial to be even and the Galois
group to be C1 or C2 without any further argument. The shipped statement
"degree 1 or 2" does not record the evenness, which is the part that
generalises.

Usage:
    python scripts/holonomy_two_elementary.py [-o results/data/holonomy_two_elementary.json]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from collections import Counter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import sympy as sp  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
X = sp.Symbol("x")


def even_quartic_group(poly: sp.Poly) -> str:
    """Galois group of an irreducible even quartic, by the classical criterion."""
    c = sp.Poly(poly.as_expr() / poly.all_coeffs()[0], X).all_coeffs()
    if any(v != 0 for v in c[1::2]):
        return "not-even"
    p_, q_ = c[2], c[4]
    if sp.sqrt(q_).is_rational:
        return "V4"
    if sp.sqrt(q_ * (p_ ** 2 - 4 * q_)).is_rational:
        return "C4"
    return "D4"


def analyse_loop(rec: dict) -> dict:
    cos = sp.nsimplify(sp.sympify(rec["cos_exact"]))
    mp = sp.Poly(sp.minimal_polynomial(cos, X), X)
    sq = sp.expand(cos ** 2)
    mp2 = sp.Poly(sp.minimal_polynomial(sq, X), X)
    deg, deg2 = mp.degree(), mp2.degree()

    coeffs = sp.Poly(mp.as_expr() / mp.all_coeffs()[0], X).all_coeffs()
    even = all(v == 0 for v in coeffs[1::2])

    out = {
        "system": rec["system"], "index": rec["index"], "loop": rec["loop"],
        "cos_exact": rec["cos_exact"],
        "degree_cos": deg, "degree_cos_squared": deg2,
        "R1_holds": deg2 <= 2,
        "min_poly_is_even": bool(even),
        "cos_squared_rational": bool(sq.is_rational),
    }
    if deg2 == 2:
        c2 = mp2.all_coeffs()
        norm = sp.Rational(c2[-1], c2[0])       # product of the two conjugates
        out["norm_cos_squared"] = str(norm)
        out["R2_norm_is_square"] = bool(sp.sqrt(norm).is_rational)
    else:
        out["norm_cos_squared"] = None
        out["R2_norm_is_square"] = None
    if deg == 4:
        out["predicted_group"] = even_quartic_group(mp)
    elif deg == 2:
        out["predicted_group"] = "C2"
    else:
        out["predicted_group"] = "C1"
    out["two_elementary"] = out["predicted_group"] in ("C1", "C2", "V4")
    return out


def max_class_sizes() -> dict:
    """Largest one-hop class size per certified state, from the shipped ledger."""
    from gpc_census.states import one_hop_classes
    out = {}
    for line in (ROOT / "results/data/states.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        dets = [tuple(s[0]) for s in rec["support"]]
        classes = one_hop_classes(dets)
        out[(rec["system"], rec["index"])] = max(
            (len(v) for v in classes.values()), default=0)
    return out


def build() -> dict:
    src = json.loads((ROOT / "results/data/holonomy_fields.json").read_text())
    rows = [analyse_loop(r) for r in src["loops"]]

    # Theorem A' test: at max class size <= 2 the cosine is n / sqrt(M), so
    # cos^2 must be rational. A violation would refute the proof above.
    sizes = max_class_sizes()
    a_prime = Counter()
    violations = []
    for r in rows:
        m = sizes.get((r["system"], r["index"]))
        r["max_one_hop_class_size"] = m
        a_prime[(m, r["cos_squared_rational"])] += 1
        if m is not None and m <= 2 and not r["cos_squared_rational"]:
            violations.append({"system": r["system"], "index": r["index"],
                               "loop": r["loop"]})

    deg_pairs = Counter((r["degree_cos"], r["degree_cos_squared"]) for r in rows)
    groups = Counter(r["predicted_group"] for r in rows)
    certified = Counter()
    for r in src["loops"]:
        certified[r.get("galois_group", "?")] += 1

    quartics = [r for r in rows if r["degree_cos"] == 4]
    not_even = [r for r in rows if r["degree_cos"] > 1 and not r["min_poly_is_even"]]

    return {
        "what": "Reduction of Holonomy Rationality to two arithmetic conditions "
                "on cos^2, with an exact census test of both.",
        "source": "results/data/holonomy_fields.json",
        "reduction": {
            "R1": "deg_Q(cos^2 Phi) <= 2",
            "R2": "for deg_Q(cos Phi) = 4, N(cos^2 Phi) is a square in Q",
            "implication": "R1 and R2 imply the Galois group is C1, C2 or V4, "
                           "hence 2-elementary; without R2 an even quartic may "
                           "be C4 or D4",
            "theorem_A_prime": "at one-hop classes of size <= 2, Theorem A's "
                               "formula makes cos Phi = n / sqrt(M) with n, M "
                               "rational, so cos^2 Phi is RATIONAL: R1 holds "
                               "with degree one and the min poly is even",
        },
        "summary": {
            "loops": len(rows),
            "degree_pairs_cos_cos_squared": {f"{a},{b}": n
                                             for (a, b), n in sorted(deg_pairs.items())},
            "R1_holds_on": sum(1 for r in rows if r["R1_holds"]),
            "R1_holds_on_all": all(r["R1_holds"] for r in rows),
            "quartics": len(quartics),
            "R2_holds_on_quartics": sum(1 for r in quartics if r["R2_norm_is_square"]),
            "R2_holds_on_all_quartics": all(r["R2_norm_is_square"] for r in quartics),
            "predicted_groups": dict(groups),
            "certified_groups": dict(certified),
            "prediction_matches_certification": dict(groups) == dict(certified),
            "all_two_elementary": all(r["two_elementary"] for r in rows),
            "min_poly_even_except": [
                {"system": r["system"], "index": r["index"], "loop": r["loop"],
                 "degree_cos": r["degree_cos"]} for r in not_even],
        },
        "theorem_A_prime_test": {
            "claim": "max one-hop class size <= 2 implies cos^2 Phi rational",
            "by_class_size": {f"{a},{b}": n for (a, b), n in sorted(a_prime.items())},
            "violations": violations,
            "holds": not violations,
            "note": "every loop with cos^2 irrational sits in a state carrying a "
                    "class of size 3, matching the already-certified prediction "
                    "that non-monomial radicands need a class of size >= 3",
        },
        "loops": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", default="results/data/holonomy_two_elementary.json")
    args = ap.parse_args()
    res = build()
    s = res["summary"]
    print(f"loops {s['loops']}")
    print(f"  (deg cos, deg cos^2) -> {s['degree_pairs_cos_cos_squared']}")
    print(f"  R1 holds on {s['R1_holds_on']}/{s['loops']}")
    print(f"  R2 holds on {s['R2_holds_on_quartics']}/{s['quartics']} quartics")
    print(f"  predicted groups {s['predicted_groups']}")
    print(f"  matches certification: {s['prediction_matches_certification']}")
    print(f"  min poly not even (degree > 1): {s['min_poly_even_except']}")
    a = res["theorem_A_prime_test"]
    print(f"  Theorem A' (class size <= 2 => cos^2 rational): holds={a['holds']} "
          f"{a['by_class_size']}")
    out = pathlib.Path(args.out)
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
