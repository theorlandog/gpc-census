"""The exact support-10 state over v103, verified independently of the script.

The census library ships a support-14 certified state for
v103 = (18^4,5^6)/34. The fixed-spectrum cascade reaches support 10, and the
endpoint's squared weights are exact rationals. This test rebuilds the 1-RDM in
exact arithmetic with its own fermionic-sign bookkeeping and re-checks the
characteristic-polynomial identity, so the certified bound s_Q(v103) <= 10 does
not rest on the script that produced it.
"""
import pathlib
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SUPPORT = [(1, 2, 3), (0, 2, 7), (0, 5, 6), (0, 1, 8), (0, 3, 4),
           (0, 4, 9), (3, 8, 9), (1, 6, 9), (3, 5, 9), (0, 1, 5)]
WEIGHTS = [sp.Rational(13, 34), sp.Rational(5, 34), sp.Rational(53, 442),
           sp.Rational(195, 1802), sp.Rational(5, 68), sp.Rational(5, 68),
           sp.Rational(35, 901), sp.Rational(1, 34), sp.Rational(18, 901),
           sp.Rational(84, 11713)]
SIGNS = [1, 1, 1, -1, 1, 1, 1, 1, 1, 1]
D = 10


def _rho_exact():
    """1-RDM by second quantization, signs derived here and not imported."""
    amps = [s * sp.sqrt(w) for s, w in zip(SIGNS, WEIGHTS)]
    pos = {t: i for i, t in enumerate(SUPPORT)}
    rho = sp.zeros(D, D)
    for col, det in enumerate(SUPPORT):
        for q in det:
            sign_out = (-1) ** det.index(q)
            reduced = tuple(u for u in det if u != q)
            for p in range(D):
                if p in reduced:
                    continue
                target = tuple(sorted(reduced + (p,)))
                row = pos.get(target)
                if row is None:
                    continue
                rho[p, q] += (amps[row] * amps[col] * sign_out
                              * (-1) ** target.index(p))
    return sp.simplify(rho)


def test_support10_state_is_normalized_exactly():
    assert sum(WEIGHTS) == 1


def test_support10_state_attains_v103_exactly():
    """det(rho - xI) = (x - 9/17)^4 (x - 5/34)^6 in exact arithmetic."""
    x = sp.Symbol("x")
    charpoly = sp.expand(_rho_exact().charpoly(x).as_expr())
    target = sp.expand((x - sp.Rational(9, 17)) ** 4 * (x - sp.Rational(5, 34)) ** 6)
    assert sp.simplify(charpoly - target) == 0


def test_support10_beats_the_certified_library_support():
    """The bound this state certifies is strictly better than the library's."""
    import json

    ledger = ROOT / "results" / "data" / "states.jsonl"
    if not ledger.exists():
        import pytest
        pytest.skip("no states atlas present")
    rec = next(json.loads(x) for x in ledger.read_text().splitlines()
               if x.strip() and json.loads(x)["system"] == "(3,10)"
               and json.loads(x)["index"] == 103)
    assert len(rec["closed_form"]["support_dets"]) == 14
    assert len(SUPPORT) < 14


def test_support10_state_is_a_new_denominator():
    """State denominator 46852 = 2^2 * 13 * 17 * 53, a new prime versus 195364."""
    den = sp.ilcm(*[w.q for w in WEIGHTS])
    assert den == 46852
    assert sp.factorint(den) == {2: 2, 13: 1, 17: 1, 53: 1}
