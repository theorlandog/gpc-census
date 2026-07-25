#!/usr/bin/env python3
"""Ship the real extremal state for v_B as the library record.

The (4,9) library record for v_B = (20,14,14,14,14,4,4,4,4)/23 previously
shipped the rational-weight fiber point (squared weights (1,8,4,3,2,2,1,2)/23,
one phase-carrying channel). The same support carries a one-parameter family
k(t) = (1+t, 8-t, 4, 3, 2-t, 2+t, 1, 2)/23 along its incidence kernel, whose
reality wall 85 t^2 - 70 t - 119 = 0 has the interior root
t0 = 7/17 - (18/85) sqrt(35). At t0 the exchange channel closes with theta = 0
and the state is REAL with all squared weights in Q(sqrt(35)).

This script constructs that endpoint, proves the exact 1-RDM
characteristic-polynomial identity in exact arithmetic (the same gate every
census state passes), and rewrites the record in results/data/states.jsonl,
the feature_table.csv all_real flag, and the SHA256SUMS lines for the changed
files. It is idempotent; rerunning reproduces the same record.
"""
import csv
import hashlib
import json
import pathlib
import sys

import sympy as sp

DATA = pathlib.Path("results/data")
SYSTEM, INDEX = "(4,9)", 65
D, DEN = 9, 23
SUPPORT = [(0, 1, 2, 4), (0, 1, 2, 8), (0, 1, 3, 5), (0, 2, 3, 6),
           (0, 3, 4, 7), (0, 3, 7, 8), (1, 3, 6, 8), (2, 3, 4, 8)]
SPECTRUM = [20, 14, 14, 14, 14, 4, 4, 4, 4]


def build_rdm(dets, amps, d):
    import bisect
    from itertools import combinations
    N = len(dets[0])
    rho = sp.zeros(d, d)
    for T, a in zip(dets, amps):
        aa = sp.simplify(a * sp.conjugate(a))
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


def certify(amps):
    # real lam: the eigenvalue variable of a Hermitian matrix; see the same
    # declaration in scripts/verify_states_standalone.py
    lam = sp.symbols("lam", real=True)
    rho = build_rdm(SUPPORT, amps, D)
    cp = (rho - lam * sp.eye(D)).det()
    tgt = sp.prod([(sp.Rational(n, DEN) - lam) for n in SPECTRUM])
    diff = sp.expand(sp.expand_complex(cp - tgt))
    # Poly normalization groups radical monomials a bare Add tree keeps
    # apart; any surviving coefficient is zero iff its minimal polynomial
    # is x (an exact algebraic-number decision, no simplification heuristics)
    x = sp.symbols("x")
    return all(c == 0 or sp.minimal_polynomial(c, x) == x
               for c in sp.Poly(diff, lam).all_coeffs())


def main() -> int:
    t0 = sp.Rational(7, 17) - sp.Rational(18, 85) * sp.sqrt(35)
    ks = [1 + t0, 8 - t0, sp.Integer(4), sp.Integer(3),
          2 - t0, 2 + t0, sp.Integer(1), sp.Integer(2)]
    assert sp.simplify(sum(ks) - DEN) == 0
    assert all(k.evalf() > 0 for k in ks), "weights must stay interior"

    # theta = 0 endpoint: both exchange-channel contributions carry the same
    # orientation; search the channel amplitude signs for the exact identity.
    base = [sp.sqrt(sp.nsimplify(k) / DEN) for k in ks]
    for signs in ([1] * 8, [1, -1, 1, 1, 1, 1, 1, 1],
                  [1, 1, 1, 1, 1, -1, 1, 1], [1, -1, 1, 1, 1, -1, 1, 1]):
        amps = [s * a for s, a in zip(signs, base)]
        if certify(amps):
            break
    else:
        print("no sign pattern certifies; aborting with data untouched")
        return 1
    print("exact characteristic-polynomial identity PASSES for signs", signs)

    weights_exact = [sp.radsimp(sp.nsimplify(k * 1)) for k in ks]
    pretty = [str(sp.radsimp(a)) for a in amps]
    support_field = [[list(T), float(abs(a.evalf(20))),
                      0.0 if float(a.evalf(20)) >= 0 else float(sp.pi.evalf(20))]
                     for T, a in zip(SUPPORT, amps)]

    path = DATA / "states.jsonl"
    lines = path.read_text().splitlines()
    hit = 0
    for i, line in enumerate(lines):
        r = json.loads(line)
        if r["system"] == SYSTEM and r["index"] == INDEX:
            assert r["integer_form"] == SPECTRUM
            r["support"] = support_field
            r["closed_form"] = {
                "den": DEN,
                "weights": [str(w) for w in weights_exact],
                "pretty": pretty,
                "support_dets": [list(T) for T in SUPPORT],
                "family": "k(t)=(1+t,8-t,4,3,2-t,2+t,1,2)/23, real endpoint "
                          "t0=7/17-18*sqrt(35)/85 (wall 85*t**2-70*t-119)",
            }
            lines[i] = json.dumps(r)
            hit += 1
    assert hit == 1, f"expected exactly one v_B record, found {hit}"
    path.write_text("\n".join(lines) + "\n")

    ft = DATA / "feature_table.csv"
    rows = list(csv.reader(ft.open()))
    for row in rows[1:]:
        if row[0] == SYSTEM and row[1] == str(INDEX):
            row[9] = "1"  # all_real: the shipped state is now real
    with ft.open("w", newline="") as fh:
        csv.writer(fh).writerows(rows)

    sums = DATA / "SHA256SUMS"
    out = []
    for line in sums.read_text().splitlines():
        digest, name = line.split(None, 1)
        p = DATA / name.strip().lstrip("./")
        out.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {name.strip()}")
    sums.write_text("\n".join(out) + "\n")

    print("weights:", [str(w) for w in weights_exact])
    print("pretty :", pretty)
    print("states.jsonl, feature_table.csv, SHA256SUMS updated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
