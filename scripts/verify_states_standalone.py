#!/usr/bin/env python3
"""Independent exact verification of every closed-form state in states.jsonl.

This checker is deliberately minimal and shares no code with the gpc-census
pipeline. For each record it:
  1. parses the symbolic amplitudes in closed_form.pretty with sympy,
  2. rebuilds the 1-RDM from scratch (independent fermionic-sign bookkeeping),
  3. verifies, in exact arithmetic, the characteristic-polynomial identity
         det(rho - lam I) == prod_m (n_m / D - lam)
     against the record's integer_form / denominator,
  4. additionally checks normalization and, for DESIGN records, that the
     support is one-hop independent (so the 1-RDM is diagonal by inspection).

Usage: python3 scripts/verify_states_standalone.py [states.jsonl]
Exit code 0 iff all 799 records pass.
"""
import bisect
import json
import sys
import time
from itertools import combinations

import sympy as sp


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
    amps = [sp.expand_complex(sp.sympify(p).rewrite(sp.cos)) for p in cf["pretty"]]
    d = len(rec["integer_form"])
    D = rec["denominator"] or 1

    norm = sp.simplify(sp.expand_complex(sum(a * sp.conjugate(a) for a in amps)))
    if sp.simplify(norm - 1) != 0:
        return False, "normalization"

    if rec["classified"].startswith("DESIGN") and not one_hop_free(dets):
        return False, "design support not one-hop independent"

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
    path = sys.argv[1] if len(sys.argv) > 1 else "results/data/states.jsonl"
    records = [json.loads(line) for line in open(path)]
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
    print("all records PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
