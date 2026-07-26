#!/usr/bin/env python3
"""SUPERSEDED: the v103 support-10 state bounds s_Q^free, not s_Q^NO.

Vertex v103 = (18^4, 5^6)/34 of Delta(3,10), library support 14.

WHAT THIS SCRIPT ORIGINALLY CLAIMED, AND WHY IT WAS WRONG. It rebuilt the
support-10 cascade endpoint from recognized rational weights and verified
normalization, hermiticity, trace and the identity
det(rho - x I) = (x - 9/17)^4 (x - 5/34)^6 in exact arithmetic, then reported
s_Q(v103) <= 10 as CERTIFIED. Every one of those checks is correct, and none of
them establishes vertex attainment: the characteristic polynomial is invariant
under conjugation, so it is satisfied by every state in the U(d) orbit,
including states whose 1-RDM is rotated off the diagonal. This one is rotated.
Its 1-RDM has rho_{3,9} = -5/68 and rho_{1,6} = -sqrt(42)/221, and its diagonal
is strictly majorized by lambda rather than equal to it.

WHAT SURVIVES. The state is a genuine exact solution with the correct 1-RDM
spectrum and multiplicities, on support 10. So

    s_Q^free(v103) <= 10   EXACT,

where s_Q^free is the minimum support over the whole U(d)-saturated set
{psi : spec rho(psi) = lambda}. Its only valid lower bound is majorization.

WHAT DOES NOT. s_Q^NO(v103), the minimum over states whose 1-RDM is diagonal
with diagonal equal to lambda, is NOT bounded by this state. That is the census
attainment convention and the quantity the diagonal-incidence baseline s_I = 6
speaks about, so this state may not be placed in a chain with s_I. The
diagonality-constrained cascade does not sparsify v103 at all: it stays at
support 14.

This script now records the corrected status rather than the original claim.

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
    offdiag_max = sp.simplify(sp.Max(*[abs(rho[i, j]) for i in range(D)
                                       for j in range(D) if i != j]))

    diag = [sp.simplify(rho[i, i]) for i in range(D)]
    offdiag = [(i, j, sp.simplify(rho[i, j])) for i in range(D) for j in range(D)
               if i != j and sp.simplify(rho[i, j]) != 0]
    checks["rdm_is_exactly_diagonal"] = not offdiag
    checks["diag_equals_lambda"] = bool(diag == list(TARGET))

    report = {
        "status": "SUPERSEDED",
        "superseded_reason":
            "The acceptance criterion used here, the characteristic-polynomial "
            "identity, is invariant under conjugation and therefore cannot "
            "establish vertex attainment. This state's 1-RDM is NOT diagonal "
            "(rho_{3,9} = -5/68, rho_{1,6} = -sqrt(42)/221) and its diagonal is "
            "strictly majorized by lambda rather than equal to it. It is an "
            "exact state on the spectrum locus, so it bounds s_Q^free, not "
            "s_Q^NO. The diagonality-constrained cascade does not sparsify v103 "
            "at all (support stays 14).",
        "vertex": "(3,10) index 103, spectrum (18,18,18,18,5,5,5,5,5,5)/34",
        "claim_retained": "s_Q^free(v103) <= 10, EXACT: an exact state of "
                          "Slater support 10 has this 1-RDM spectrum with the "
                          "correct multiplicities",
        "claim_retracted": "s_Q^NO(v103) <= 10. Not supported by this state, "
                           "and not currently supported by anything.",
        "evidence": "EXACT for the spectrum-locus claim (rational and radical "
                    "arithmetic, no floating point)",
        "provenance": "fixed-spectrum sparsification cascade endpoint, weights "
                      "recognized by high-precision integer relation "
                      "(scripts/v103_endpoint_precision.py)",
        "certified_library_support": 14,
        "support_dets": [list(t) for t in SUPPORT],
        "squared_weights": [str(w) for w in SQUARED_WEIGHTS],
        "signs": SIGNS,
        "state_denominator": int(state_den),
        "state_denominator_factored": str(sp.factorint(state_den)),
        "library_state_denominator": 195364,
        "loop_holonomy": "pi (state is real up to the orbital-phase gauge)",
        "rdm_max_offdiagonal": str(offdiag_max),
        "rdm_nonzero_offdiagonals": [f"rho[{i},{j}] = {e}" for i, j, e in offdiag],
        "rdm_diagonal": [str(e) for e in diag],
        "target_lambda": [str(e) for e in TARGET],
        "checks": checks,
        "all_checks_pass": all(checks[k] for k in checks if k != "support_size"),
        "caveat": "an upper bound on s_Q^free, not a minimality claim, and not "
                  "an attainment claim in the census diagonal convention",
    }
    (ROOT / "results/data/v103_support10_certificate.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(f"status: {report['status']}")
    for k, v in checks.items():
        print(f"  {k}: {v}")
    print(f"  nonzero off-diagonals: {report['rdm_nonzero_offdiagonals']}")
    print(f"retained: {report['claim_retained']}")
    print(f"retracted: {report['claim_retracted']}")


if __name__ == "__main__":
    main()
