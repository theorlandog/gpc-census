#!/usr/bin/env python3
"""Separately implemented exact verification of the state atlas.

This checker is deliberately minimal and shares no code with the gpc-census
pipeline. For each record it:
  1. parses the symbolic amplitudes in closed_form.pretty with sympy,
  2. rebuilds the 1-RDM from scratch (independent fermionic-sign bookkeeping),
  3. verifies, in exact arithmetic, the characteristic-polynomial identity
         det(rho - lam I) == prod_m (n_m / D - lam)
     against the record's integer_form / denominator,
  4. checks normalization, support structure, and agreement among the three
     serialized amplitude representations,
  5. for design records, checks one-hop independence and the advertised
     integer-grid or rational off-grid arithmetic of the displayed witness.

Usage: python3 scripts/verify_states_standalone.py [states.jsonl]
By default the command also requires exactly 799 unique records.  Use
``--allow-subset`` only for focused debugging.
"""
import argparse
import bisect
import json
import re
import time
from itertools import combinations
from pathlib import Path

import sympy as sp


ALLOWED_EXPRESSION = re.compile(r"^[0-9A-Za-z_+*/().\- ]+$")
ALLOWED_NAMES = {
    "sqrt": sp.sqrt,
    "exp": sp.exp,
    "acos": sp.acos,
    "atan": sp.atan,
    "I": sp.I,
    "pi": sp.pi,
}


def parse_exact_expression(text):
    value = str(text)
    if not ALLOWED_EXPRESSION.fullmatch(value):
        raise ValueError(f"expression contains unsupported characters: {value!r}")
    names = set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", value))
    unknown = names - set(ALLOWED_NAMES)
    if unknown:
        raise ValueError(f"expression contains unsupported names: {sorted(unknown)}")
    expression = sp.sympify(value, locals=ALLOWED_NAMES)
    if expression.free_symbols:
        raise ValueError(f"expression contains free symbols: {expression.free_symbols}")
    return expression


def build_rdm(dets, amps, d):
    N = len(dets[0])
    rho = sp.zeros(d, d)
    for T, a in zip(dets, amps):
        aa = sp.simplify(sp.expand_complex(a * sp.conjugate(a)))
        for m in T:
            rho[m, m] += aa
    for (T1, a1), (T2, a2) in combinations(zip(dets, amps), 2):
        s1, s2 = set(T1), set(T2)
        if len(s1 & s2) == N - 1:
            i = (s1 - s2).pop()
            j = (s2 - s1).pop()
            T2s = sorted(T2)
            pj = T2s.index(j)
            rest = [x for x in T2s if x != j]
            pi = bisect.bisect_left(rest, i)
            sgn = (-1) ** (pj + pi)
            rho[i, j] += sgn * sp.conjugate(a1) * a2
            rho[j, i] += sgn * sp.conjugate(a2) * a1
    return rho


def one_hop_free(dets):
    N = len(dets[0])
    return all(len(set(a) & set(b)) < N - 1 for a, b in combinations(dets, 2))


def verify(rec):
    cf = rec["closed_form"]
    dets = [tuple(t) for t in cf["support_dets"]]
    try:
        amps = [
            sp.expand_complex(parse_exact_expression(p).rewrite(sp.cos))
            for p in cf["pretty"]
        ]
        weight_numerators = [parse_exact_expression(value) for value in cf["weights"]]
    except (TypeError, ValueError, sp.SympifyError) as error:
        return False, f"closed-form expression grammar: {error}"
    d = len(rec["integer_form"])
    D = rec["denominator"] or 1
    try:
        n = int(rec["system"].strip("()").split(",")[0])
    except (AttributeError, ValueError):
        return False, "malformed system key"

    if not dets or not (len(dets) == len(amps) == len(weight_numerators)):
        return False, "closed-form array lengths"
    if cf.get("den", 0) <= 0:
        return False, "closed-form denominator"
    if len(set(dets)) != len(dets):
        return False, "duplicate support determinant"
    if any(
        len(det) != n
        or tuple(sorted(det)) != det
        or len(set(det)) != n
        or any(index < 0 or index >= d for index in det)
        for det in dets
    ):
        return False, "malformed support determinant"
    top_support = rec.get("support")
    if top_support is not None and [tuple(item[0]) for item in top_support] != dets:
        return False, "top-level and closed-form supports disagree"

    weights = [sp.simplify(value / cf["den"]) for value in weight_numerators]
    for amp, weight in zip(amps, weights):
        amp_weight = sp.simplify(sp.expand_complex(amp * sp.conjugate(amp)))
        if sp.simplify(amp_weight - weight) != 0:
            return False, "pretty/weights/den disagreement"
        if weight.is_positive is not True:
            return False, "nonpositive serialized weight"

    norm = sp.simplify(sp.expand_complex(sum(a * sp.conjugate(a) for a in amps)))
    if sp.simplify(norm - 1) != 0:
        return False, "normalization"

    classification = rec["classified"]
    if classification.startswith("DESIGN"):
        if not one_hop_free(dets):
            return False, "design support not one-hop independent"
        scaled = [sp.simplify(weight * D) for weight in weights]
        if classification == "DESIGN-INT":
            if any(value.is_integer is not True or value.is_positive is not True for value in scaled):
                return False, "DESIGN-INT witness is not on the spectrum-denominator grid"
        elif classification == "DESIGN-REAL":
            if any(weight.is_rational is not True for weight in weights):
                return False, "DESIGN-REAL witness is not rational"
            if all(value.is_integer is True for value in scaled):
                return False, "DESIGN-REAL witness is not off the spectrum-denominator grid"
        else:
            return False, f"unknown design classification {classification!r}"

    rho = build_rdm(dets, amps, d)
    # lam is the eigenvalue variable of a Hermitian matrix, hence real; a
    # generic complex symbol makes expand_complex split it into re/im parts,
    # which blocks exact cancellation of nested-radical amplitude products.
    lam = sp.symbols("lam", real=True)
    cp = (rho - lam * sp.eye(d)).det()
    tgt = sp.prod([(sp.Rational(n, D) - lam) for n in rec["integer_form"]])
    diff = sp.expand(sp.expand_complex(cp - tgt))
    if diff != 0:
        diff = sp.simplify(diff)
    if diff != 0:
        # Nested-radical amplitudes (e.g. the real v_B library state, weights
        # in Q(sqrt(35))) defeat generic simplify. Decide exactly instead:
        # Poly normalization groups radical monomials the bare Add tree keeps
        # apart, and a surviving coefficient is zero iff its minimal
        # polynomial is x (an exact algebraic-number decision).
        x = sp.symbols("x")
        if not all(c == 0 or sp.minimal_polynomial(c, x) == x
                   for c in sp.Poly(diff, lam).all_coeffs()):
            return False, "characteristic-polynomial identity"
    return True, "ok"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("states", nargs="?", type=Path,
                        default=Path("results/data/states.jsonl"))
    parser.add_argument("--allow-subset", action="store_true",
                        help="do not require the complete 799-record atlas")
    args = parser.parse_args()
    records = [
        json.loads(line)
        for line in args.states.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not args.allow_subset and len(records) != 799:
        print(f"FAIL: expected 799 records, found {len(records)}")
        return 1
    keys = [(record.get("system"), record.get("index")) for record in records]
    if len(set(keys)) != len(keys):
        print("FAIL: duplicate system/index keys")
        return 1
    t0 = time.time()
    failures = []
    for k, rec in enumerate(records):
        ok, why = verify(rec)
        if not ok:
            failures.append((rec["system"], rec["index"], rec["integer_form"], why))
        if (k + 1) % 100 == 0:
            print(f"  {k + 1}/{len(records)} checked ({time.time() - t0:.0f}s)",
                  flush=True)
    print(f"checked {len(records)} records in {time.time() - t0:.0f}s")
    if failures:
        print(f"FAILURES ({len(failures)}):")
        for f in failures:
            print("  ", f)
        return 1
    print(f"all {len(records)} records PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
