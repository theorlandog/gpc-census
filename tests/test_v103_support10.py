"""The v103 support-10 state: exact on the SPECTRUM locus, not an attainer.

CORRECTED. An earlier revision of this branch published this state as an exact
certificate for s_Q(v103) <= 10. That was wrong, and the error is instructive:
the acceptance criterion was the characteristic-polynomial identity, which is
invariant under conjugation and therefore blind to the basis. The state is a
genuine exact solution with the right 1-RDM SPECTRUM, so it bounds

    s_Q^free(v103) <= 10,

the minimum over the whole U(d)-saturated set. It does NOT bound

    s_Q^NO(v103),

the minimum over states whose 1-RDM is diagonal with diagonal equal to lambda,
which is the census attainment convention and the one the incidence baseline
s_I speaks about. Its 1-RDM has rho_{3,9} = -5/68 and a diagonal strictly
majorized by lambda.

The rejection itself is pinned in tests/test_attainment_criterion.py. What is
tested here is the positive content that survives.
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


def test_state_is_normalized_exactly():
    assert sum(WEIGHTS) == 1


def test_state_has_the_vertex_spectrum_exactly():
    """The positive content: spec(rho) = lambda exactly, so s_Q^free <= 10."""
    x = sp.Symbol("x")
    charpoly = sp.expand(_rho_exact().charpoly(x).as_expr())
    target = sp.expand((x - sp.Rational(9, 17)) ** 4 * (x - sp.Rational(5, 34)) ** 6)
    assert sp.simplify(charpoly - target) == 0


def test_state_bounds_free_support_not_natural_orbital_support():
    """The whole point of the correction, asserted rather than commented."""
    rho = _rho_exact()
    lam = [sp.Rational(18, 34)] * 4 + [sp.Rational(5, 34)] * 6
    offdiag = [sp.simplify(rho[i, j]) for i in range(D) for j in range(D) if i != j]
    assert any(e != 0 for e in offdiag), "not diagonal, so not an s_Q^NO witness"
    assert [sp.simplify(rho[i, i]) for i in range(D)] != lam
    assert len(SUPPORT) < 14  # still strictly better than the library support


def test_state_denominator_is_recorded_correctly():
    """46852 = 2^2 * 13 * 17 * 53. Arithmetic stands; its meaning is s_Q^free."""
    den = sp.ilcm(*[w.q for w in WEIGHTS])
    assert den == 46852
    assert sp.factorint(den) == {2: 2, 13: 1, 17: 1, 53: 1}


def test_shipped_certified_bounds_are_genuine_attainers():
    """Every CERTIFIED entry in cascade_exact.json must be diagonal at lambda.

    The strengthened checker accepts only states with rho exactly diagonal and
    diag(rho) = lambda. This re-derives that from scratch for a sample, so a
    regression to the conjugation-invariant criterion cannot ship.
    """
    import json
    from fractions import Fraction

    import pytest

    data = ROOT / "results" / "data" / "cascade_exact.json"
    states = ROOT / "results" / "data" / "states.jsonl"
    if not data.exists() or not states.exists():
        pytest.skip("no cascade dataset present")
    report = json.loads(data.read_text())
    certified = [r for r in report["results"] if r["status"] == "CERTIFIED"]
    if not certified:
        pytest.skip("no s_Q^NO bounds certified in the current dataset")

    ledger = {}
    for line in states.read_text().splitlines():
        if line.strip():
            rec = json.loads(line)
            ledger[(rec["system"], rec["index"])] = rec

    for r in sorted(certified, key=lambda z: z["support_upper"])[:4]:
        d = int(r["system"].strip("()").split(",")[1])
        dets = [tuple(t) for t in r["support_dets"]]
        weights = [sp.Rational(w) for w in r["squared_weights"]]
        assert sum(weights) == 1
        amps = [s * sp.sqrt(w) for s, w in zip(r["signs"], weights)]
        pos = {t: i for i, t in enumerate(dets)}
        rho = sp.zeros(d, d)
        for col, det in enumerate(dets):
            for q in det:
                sign_out = (-1) ** det.index(q)
                reduced = tuple(u for u in det if u != q)
                for p in range(d):
                    if p in reduced:
                        continue
                    target = tuple(sorted(reduced + (p,)))
                    row = pos.get(target)
                    if row is None:
                        continue
                    rho[p, q] += (amps[row] * amps[col] * sign_out
                                  * (-1) ** target.index(p))
        rho = sp.simplify(rho)
        src = ledger[(r["system"], r["index"])]
        den = int(src["denominator"])
        lam = [sp.Rational(int(v), den) for v in src["integer_form"]]
        offdiag = [sp.simplify(rho[i, j]) for i in range(d) for j in range(d)
                   if i != j]
        assert all(e == 0 for e in offdiag), (r["system"], r["index"], "not diagonal")
        assert [sp.simplify(rho[i, i]) for i in range(d)] == lam, (
            r["system"], r["index"], "diagonal != lambda")
        assert isinstance(Fraction(den), Fraction)
