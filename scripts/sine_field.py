#!/usr/bin/env python3
"""Theorem F: the sine field, and why 17 of 18 triples degenerate.

`docs/class_closure.md` closed every large class, with 17 of 18 chained triples
closing because a sine vanished and exactly ONE closing because its
discriminant was a genuine rational square. The 1 carries the mechanism; the 17
are its degenerate limit.

THE OBJECT. Write `A_ab = Im(z_a conj(z_b))`, so

    z_a conj(z_b) = t_ab + i A_ab,     t_ab^2 + A_ab^2 = r_a^2 r_b^2,

and `A_ab = r_a r_b sin(Phi_ab)` is twice the signed area of the triangle
`(0, z_a, z_b)`. The Gram entry and the area are the real and imaginary parts of
one product, and Theorem E's discriminant is built from the AREAS:

    Delta / 4 = (t_ab^2 - r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2) = A_ab^2 A_bc^2.

So `sqrt(Delta)/2 = |A_ab A_bc|` and Theorem E forces `t_ac` rational exactly
when `A_ab A_bc` is rational, which happens exactly when `sin^2(Phi_ab)` and
`sin^2(Phi_bc)` share a squarefree part. That is a condition on `Q(sin Phi)`,
the SINE field, which is not the holonomy field `Q(cos Phi)` the arithmetic arc
has studied and which `holonomy_fields.json` already records separately as
`sin_radicands`.

PROPAGATION. The plane sine-addition identity reads

    r_b^2 A_ac = A_ab t_bc + t_ab A_bc,

so with `t_ab` and `t_bc` rational the sine field propagates: `A_ab` and `A_bc`
in `Q(sqrt d)` puts `A_ac` there too. If the sine field is a class invariant,
every Theorem E discriminant in that class is automatically a rational square
and no case-by-case squareness check is needed.

WHY THE 17 DEGENERATE. A state that is real up to gauge has every holonomy at
`0` or `pi`, so every sine is zero, every area vanishes, and every discriminant
is trivially a square. The 17 are that limit; the 1 shows the general shape.

Usage:
    python scripts/sine_field.py
"""
from __future__ import annotations

import argparse
import collections
import itertools
import json
import pathlib
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SOURCE = ROOT / "results/data/shared_loop_rationality.json"
HOLONOMY = ROOT / "results/data/holonomy_fields.json"
X = sp.Symbol("x")


def squarefree(value) -> int:
    """Squarefree part of a nonnegative rational, as an integer; 0 stays 0."""
    q = sp.Rational(value)
    if q == 0:
        return 0
    out = 1
    for part in (q.p, q.q):
        for prime, exponent in sp.factorint(abs(part)).items():
            if exponent % 2:
                out *= prime
    return out


def _degree(value) -> int | None:
    try:
        return sp.Poly(sp.minimal_polynomial(value, X), X).degree()
    except Exception:
        return None


def in_field_of(value, generator) -> bool | None:
    """Whether ``value`` lies in Q(generator), decided by degree comparison.

    ``Q(generator, value) == Q(generator)`` exactly when adjoining ``value``
    does not raise the degree, so the primitive element of the pair has the
    same degree as the generator alone.
    """
    try:
        base = _degree(generator)
        joint = _degree(sp.simplify(generator + sp.pi * 0 + value * 0 + value))
        if base is None or joint is None:
            return None
        together, _ = sp.primitive_element([sp.nsimplify(generator),
                                            sp.nsimplify(value)])
        return sp.Poly(together, X).degree() == base
    except Exception:
        return None


def analyse_class(block) -> dict:
    r2 = [sp.nsimplify(sp.sympify(v)) for v in block["r_squared"]]
    t = {tuple(e["pair"]): sp.nsimplify(sp.sympify(e["t_ab"]))
         for e in block["links"]}

    pairs = []
    for key, value in sorted(t.items()):
        a, b = key
        area2 = sp.radsimp(sp.simplify(r2[a] * r2[b] - value**2))
        entry = {
            "pair": list(key),
            "t_ab": str(value),
            "t_rational": bool(value.is_rational),
            "area_squared": str(area2),
            "area_squared_rational": bool(area2.is_rational),
            "area_rational": bool(sp.sqrt(area2).is_rational),
            "area_vanishes": bool(area2 == 0),
        }
        if area2.is_rational:
            entry["squarefree"] = squarefree(area2)
        pairs.append(entry)

    known = {p["squarefree"] for p in pairs
             if p.get("squarefree") not in (None, 0)}
    invariant = known.pop() if len(known) == 1 else (0 if not known else None)

    # P-F-4 on EVERY pair, including those with irrational t: is A/sqrt(d) in
    # Q(t_ab)? With d = 0 the class is degenerate and the test is vacuous.
    for entry in pairs:
        area2 = sp.nsimplify(sp.sympify(entry["area_squared"]))
        value = sp.nsimplify(sp.sympify(entry["t_ab"]))
        if invariant in (None, 0) or area2 == 0:
            entry["in_sine_field"] = None
            continue
        scaled = sp.radsimp(sp.simplify(sp.sqrt(area2 / invariant)))
        entry["area_over_sqrt_d"] = str(scaled)
        entry["in_sine_field"] = in_field_of(scaled, value)

    # Theorem F identities on every triple
    triples = []
    for b in range(block["class_size"]):
        for a, c in itertools.combinations(
                [k for k in range(block["class_size"]) if k != b], 2):
            key = lambda i, j: (min(i, j), max(i, j))  # noqa: E731
            t_ab, t_bc, t_ac = t[key(a, b)], t[key(b, c)], t[key(a, c)]
            A2_ab = sp.simplify(r2[a] * r2[b] - t_ab**2)
            A2_bc = sp.simplify(r2[b] * r2[c] - t_bc**2)
            A2_ac = sp.simplify(r2[a] * r2[c] - t_ac**2)
            delta = sp.simplify(4 * (t_ab**2 - r2[a] * r2[b])
                                * (t_bc**2 - r2[b] * r2[c]))
            # sine addition, up to the orientation sign of each signed area
            lhs = sp.simplify(r2[b] * sp.sqrt(A2_ac))
            candidates = [
                sp.simplify(sa * sp.sqrt(A2_ab) * t_bc + sb * t_ab * sp.sqrt(A2_bc))
                for sa in (1, -1) for sb in (1, -1)
            ]
            triples.append({
                "triple": [a, b, c],
                "delta_equals_area_product": bool(
                    sp.simplify(delta - 4 * A2_ab * A2_bc) == 0),
                "sine_addition_holds": bool(
                    any(sp.simplify(lhs - x) == 0 for x in candidates)),
                "both_t_rational": bool(t_ab.is_rational and t_bc.is_rational),
                "degenerate": bool(A2_ab == 0 or A2_bc == 0),
            })

    return {
        "system": block["system"],
        "index": block["index"],
        "channel": block["channel"],
        "class_size": block["class_size"],
        "transport_class": block["transport_class"],
        "sine_field_squarefree": invariant,
        "sine_field_is_uniform": invariant is not None,
        "pairs": pairs,
        "triples": triples,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out",
                    default=str(ROOT / "results/data/sine_field.json"))
    args = ap.parse_args()

    source = json.loads(SOURCE.read_text())
    holonomy = json.loads(HOLONOMY.read_text())
    recorded = collections.defaultdict(set)
    for loop in holonomy["loops"]:
        for radicand in loop["sin_radicands"]:
            for prime, exponent in sp.factorint(int(radicand)).items():
                if exponent % 2:
                    recorded[(loop["system"], loop["index"])].add(prime)

    blocks = []
    for block in source["classes"]:
        record = analyse_class(block)
        primes = set()
        d = record["sine_field_squarefree"]
        if d:
            primes = {p for p, e in sp.factorint(int(d)).items() if e % 2}
        seen = recorded.get((record["system"], record["index"]), set())
        record["sine_radicand_primes"] = sorted(primes)
        record["recorded_sin_radicand_primes"] = sorted(seen)
        record["matches_recorded_sin_radicands"] = bool(
            not primes or primes <= seen)
        blocks.append(record)
        print(f"{record['system']} v{record['index']} ch{record['channel']} "
              f"d={record['sine_field_squarefree']} "
              f"matches_recorded={record['matches_recorded_sin_radicands']}",
              flush=True)

    all_pairs = [p for b in blocks for p in b["pairs"]]
    tested = [p for p in all_pairs if p.get("in_sine_field") is not None]
    all_triples = [x for b in blocks for x in b["triples"]]

    payload = {
        "what": "Theorem F: Delta/4 = A_ab^2 A_bc^2 with A the signed area, so "
                "Theorem E's rationality condition is a statement about the "
                "SINE field Q(sin Phi), not the holonomy field Q(cos Phi).",
        "theorem_F": {
            "identity": "z_a conj(z_b) = t_ab + i A_ab with "
                        "t_ab^2 + A_ab^2 = r_a^2 r_b^2",
            "discriminant": "Delta/4 = A_ab^2 A_bc^2, so t_ac is rational iff "
                            "A_ab A_bc is rational iff sin^2(Phi_ab) and "
                            "sin^2(Phi_bc) share a squarefree part",
            "propagation": "r_b^2 A_ac = A_ab t_bc + t_ab A_bc, so with t_ab "
                           "and t_bc rational the sine field propagates",
            "degenerate_limit": "a state real up to gauge has every sine zero, "
                                "which is why 17 of 18 triples degenerate",
        },
        "summary": {
            "classes": len(blocks),
            "transport_classes": len({b["transport_class"] for b in blocks}),
            "sine_field_uniform": sum(1 for b in blocks
                                      if b["sine_field_is_uniform"]),
            "sine_field_values": sorted(
                {b["sine_field_squarefree"] for b in blocks
                 if b["sine_field_squarefree"] is not None}),
            "value_to_transport_classes": {
                str(v): len({b["transport_class"] for b in blocks
                             if b["sine_field_squarefree"] == v})
                for v in sorted({b["sine_field_squarefree"] for b in blocks
                                 if b["sine_field_squarefree"] is not None})},
            "degenerate_classes": sum(1 for b in blocks
                                      if b["sine_field_squarefree"] == 0),
            "matches_recorded_sin_radicands": sum(
                1 for b in blocks if b["matches_recorded_sin_radicands"]),
            "pairs": len(all_pairs),
            "pairs_with_rational_area": sum(1 for p in all_pairs
                                            if p["area_rational"]),
            "pairs_with_vanishing_area": sum(1 for p in all_pairs
                                             if p["area_vanishes"]),
            "in_sine_field_tested": len(tested),
            "in_sine_field_holds": sum(1 for p in tested if p["in_sine_field"]),
            "triples": len(all_triples),
            "delta_equals_area_product": sum(
                1 for x in all_triples if x["delta_equals_area_product"]),
            "sine_addition_holds": sum(1 for x in all_triples
                                       if x["sine_addition_holds"]),
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
