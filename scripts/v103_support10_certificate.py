#!/usr/bin/env python3
"""EXACT certificate: v103 admits an extremal state of support 10.

Vertex v103 = (18^4, 5^6)/34 of Delta(3,10). Its certified library state has
support 14. The fixed-spectrum sparsification cascade (docs/cascade_census.md)
reaches support 10 numerically; high-precision refinement plus integer-relation
recognition (scripts/v103_endpoint_precision.py) identifies every squared weight
of that endpoint as a rational. This script closes the loop in EXACT arithmetic:
it rebuilds the 1-RDM from the exact rational weights and signs with independent
fermionic-sign bookkeeping and verifies the characteristic-polynomial identity

    det(rho - x I) = (x - 9/17)^4 (x - 5/34)^6

with no floating point anywhere. That makes s_Q(v103) <= 10 a CERTIFIED bound,
not a numerical one. It says nothing about whether 10 is minimal.

  python scripts/v103_support10_certificate.py

Writes results/data/v103_support10_certificate.json.
"""
from __future__ import annotations

import json
import pathlib

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = 10
N = 3

# The support-10 endpoint, gauge-fixed so every amplitude is real (its single
# loop holonomy is pi, so the state is real up to the orbital-phase gauge).
SUPPORT = [(1, 2, 3), (0, 2, 7), (0, 5, 6), (0, 1, 8), (0, 3, 4),
           (0, 4, 9), (3, 8, 9), (1, 6, 9), (3, 5, 9), (0, 1, 5)]
SQUARED_WEIGHTS = [sp.Rational(13, 34), sp.Rational(5, 34), sp.Rational(53, 442),
                   sp.Rational(195, 1802), sp.Rational(5, 68), sp.Rational(5, 68),
                   sp.Rational(35, 901), sp.Rational(1, 34), sp.Rational(18, 901),
                   sp.Rational(84, 11713)]
SIGNS = [1, 1, 1, -1, 1, 1, 1, 1, 1, 1]
TARGET = [sp.Rational(18, 34)] * 4 + [sp.Rational(5, 34)] * 6


def one_rdm(amplitudes, dets, d):
    """<a_p^dag a_q> in exact arithmetic, fermionic signs from scratch."""
    idx = {t: i for i, t in enumerate(dets)}
    rho = sp.zeros(d, d)
    for ti, t in enumerate(dets):
        for qi, q in enumerate(t):
            rest = [u for u in t if u != q]
            for p in range(d):
                if p in rest:
                    continue
                s = tuple(sorted(rest + [p]))
                si = idx.get(s)
                if si is None:
                    continue
                rho[p, q] += (amplitudes[si] * amplitudes[ti]
                              * (-1) ** qi * (-1) ** s.index(p))
    return sp.simplify(rho)


def main():
    checks = {}
    checks["support_size"] = len(SUPPORT)
    checks["normalization_exact"] = bool(sum(SQUARED_WEIGHTS) == 1)

    amps = [s * sp.sqrt(w) for s, w in zip(SIGNS, SQUARED_WEIGHTS)]
    rho = one_rdm(amps, SUPPORT, D)

    x = sp.Symbol("x")
    charpoly = sp.expand(rho.charpoly(x).as_expr())
    target = sp.expand(sp.prod([(x - lam) for lam in TARGET]))
    checks["charpoly_identity_exact"] = bool(sp.simplify(charpoly - target) == 0)
    checks["trace_equals_particle_number"] = bool(sp.simplify(sp.trace(rho) - N) == 0)
    checks["hermitian_exact"] = bool(sp.simplify(rho - rho.T) == sp.zeros(D, D))

    state_den = sp.ilcm(*[w.q for w in SQUARED_WEIGHTS])
    offdiag = sp.simplify(sp.Max(*[abs(rho[i, j]) for i in range(D)
                                   for j in range(D) if i != j]))

    report = {
        "vertex": "(3,10) index 103, spectrum (18,18,18,18,5,5,5,5,5,5)/34",
        "claim": "s_Q(v103) <= 10: an extremal state of Slater support 10 attains "
                 "this vertex exactly",
        "evidence": "EXACT (rational and radical arithmetic, no floating point)",
        "provenance": "fixed-spectrum sparsification cascade endpoint "
                      "(docs/cascade_census.md), weights recognized by "
                      "high-precision integer relation "
                      "(scripts/v103_endpoint_precision.py), verified here exactly",
        "certified_library_support": 14,
        "support_dets": [list(t) for t in SUPPORT],
        "squared_weights": [str(w) for w in SQUARED_WEIGHTS],
        "signs": SIGNS,
        "state_denominator": int(state_den),
        "state_denominator_factored": str(sp.factorint(state_den)),
        "library_state_denominator": 195364,
        "loop_holonomy": "pi (state is real up to the orbital-phase gauge)",
        "rdm_max_offdiagonal": str(offdiag),
        "checks": checks,
        "all_checks_pass": all(checks[k] for k in checks if k != "support_size"),
        "caveat": "an upper bound on the minimal support, not a minimality claim; "
                  "the endpoint is first-order rigid in both fiber senses but that "
                  "does not prove no sparser state exists",
    }
    (ROOT / "results/data/v103_support10_certificate.json").write_text(
        json.dumps(report, indent=2) + "\n")
    for k, v in checks.items():
        print(f"  {k}: {v}")
    print(f"all_checks_pass: {report['all_checks_pass']}")
    print(f"state denominator: {state_den} = {sp.factorint(state_den)}")


if __name__ == "__main__":
    main()
