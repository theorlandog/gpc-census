#!/usr/bin/env python3
"""Certificate for s_Q^free(v103) <= 10, under the claim that actually holds.

The support-10 state over (3,10) v103 is correct and its arithmetic is exact.
What was wrong was the claim attached to it: it was first published as an
s_Q^NO bound, and s_Q^NO requires a 1-RDM that is exactly diag(lambda). This
state's 1-RDM has two nonzero off-diagonal entries, so it attains the
SPECTRUM in a rotated basis and bounds s_Q^free instead. The older artifact,
results/data/v103_support10_certificate.json, records that history and opens
with status SUPERSEDED.

This file re-certifies the SAME state under the CORRECT claim, so that the
manuscript's certificate table can cite an artifact that asserts what it
proves rather than one that announces its own supersession. Nothing about the
state changes; only the claim it is attached to.

What is certified, all in exact arithmetic:

1. the squared weights are rationals summing to 1;
2. the 1-RDM rebuilt from the closed form has characteristic polynomial
   (x - 9/17)^4 (x - 5/34)^6, so spec(rho) = lambda exactly;
3. the 1-RDM is NOT diagonal, and both nonzero off-diagonal entries are
   reported exactly, which is precisely why the bound is s_Q^free and not
   s_Q^NO;
4. the single loop holonomy is pi, so the state is real up to gauge.

Usage:
  python scripts/v103_sqfree_certificate.py
  python scripts/v103_sqfree_certificate.py --check
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SOURCE = DATA / "v103_support10_certificate.json"
OUT = DATA / "v103_sqfree_certificate.json"

D = 10
LAMBDA = [sp.Rational(9, 17)] * 4 + [sp.Rational(5, 34)] * 6


def one_rdm_exact(amps, dets, d):
    """1-RDM in exact arithmetic, with fermionic signs derived here."""
    index = {t: i for i, t in enumerate(dets)}
    rho = sp.zeros(d, d)
    for ti, t in enumerate(dets):
        for qi, q in enumerate(t):
            rest = [u for u in t if u != q]
            for p in range(d):
                if p in rest:
                    continue
                s = tuple(sorted(rest + [p]))
                si = index.get(s)
                if si is None:
                    continue
                rho[p, q] += (
                    sp.conjugate(amps[si]) * amps[ti] * (-1) ** qi * (-1) ** s.index(p)
                )
    return sp.simplify(rho)


def build():
    src = json.loads(SOURCE.read_text())
    dets = [tuple(t) for t in src["support_dets"]]
    weights = [sp.sympify(w) for w in src["squared_weights"]]
    signs = src["signs"]
    amps = [signs[i] * sp.sqrt(weights[i]) for i in range(len(weights))]

    rho = one_rdm_exact(amps, dets, D)
    x = sp.Symbol("x")
    charpoly = sp.factor(sp.expand((rho - x * sp.eye(D)).det()))
    target = sp.expand((sp.Rational(9, 17) - x) ** 4 * (sp.Rational(5, 34) - x) ** 6)
    spectrum_ok = sp.simplify(sp.expand(charpoly) - target) == 0

    offdiag = {}
    for i in range(D):
        for j in range(D):
            if i == j:
                continue
            value = sp.radsimp(sp.simplify(rho[i, j]))
            if value != 0:
                offdiag[f"rho[{i},{j}]"] = str(value)

    # the single loop holonomy
    incidence = sp.Matrix([[1 if o in t else 0 for o in range(D)] for t in dets])
    loops = incidence.T.nullspace()
    holonomies = []
    for vec in loops:
        vec = vec * sp.lcm([q.q for q in vec])
        vec = vec / sp.gcd([sp.Integer(q) for q in vec])
        k = [int(q) for q in vec]
        phase = sum(k[i] * sp.arg(amps[i]) for i in range(len(k)))
        holonomies.append(str(sp.simplify(phase)))

    return {
        "generated_by": "scripts/v103_sqfree_certificate.py",
        "status": "CURRENT",
        "vertex": "(3,10) v103, spectrum (18^4, 5^6)/34",
        "claim": "s_Q^free(v103) <= 10",
        "claim_is_not": (
            "s_Q^NO(v103) <= 10. That stronger claim is FALSE for this state "
            "and was the error the superseded artifact records: s_Q^NO needs a "
            "1-RDM exactly equal to diag(lambda), and this one is not diagonal."
        ),
        "supersedes_the_claim_in": "results/data/v103_support10_certificate.json",
        "note_on_the_superseded_file": (
            "That file is retained deliberately as documented history. The "
            "state in it is correct; only the claim attached to it was wrong."
        ),
        "evidence": "EXACT (rational and radical arithmetic throughout)",
        "support": len(dets),
        "support_dets": [list(t) for t in dets],
        "squared_weights": [str(w) for w in weights],
        "signs": signs,
        "weights_sum_to_one": bool(sp.simplify(sum(weights) - 1) == 0),
        "weights_all_rational": bool(all(w.is_Rational for w in weights)),
        "state_denominator": src["state_denominator"],
        "characteristic_polynomial": str(charpoly),
        "spectrum_equals_lambda": bool(spectrum_ok),
        "rdm_is_diagonal": not offdiag,
        "rdm_nonzero_offdiagonals": offdiag,
        "why_s_Q_free_and_not_s_Q_NO": (
            "The 1-RDM has two nonzero off-diagonal entries, "
            "rho[3,9] = -5/68 and rho[1,6] = -sqrt(42)/221. So the state "
            "attains the spectrum in a basis other than the natural-orbital "
            "one, and its support bounds s_Q^free. The s_Q^NO bound for this "
            "vertex remains the library value 14."
        ),
        "loop_holonomies": holonomies,
        "real_up_to_gauge": bool(
            all(sp.simplify(sp.cos(sp.sympify(h)) ** 2 - 1) == 0 for h in holonomies)
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    doc = build()
    text = json.dumps(doc, indent=1) + "\n"
    if a.check:
        cur = Path(a.out).read_text() if Path(a.out).exists() else ""
        if cur != text:
            print(f"STALE: {a.out}", file=sys.stderr)
            return 1
        print(f"{a.out} is current")
        return 0

    Path(a.out).write_text(text)
    print(f"wrote {a.out}")
    print(f"  claim: {doc['claim']}")
    print(f"  weights rational and sum to 1: {doc['weights_all_rational']} "
          f"{doc['weights_sum_to_one']}")
    print(f"  spectrum equals lambda (exact charpoly): {doc['spectrum_equals_lambda']}")
    print(f"  1-RDM diagonal: {doc['rdm_is_diagonal']}  "
          f"off-diagonals: {list(doc['rdm_nonzero_offdiagonals'])}")
    print(f"  real up to gauge: {doc['real_up_to_gauge']}  "
          f"holonomies {doc['loop_holonomies']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
