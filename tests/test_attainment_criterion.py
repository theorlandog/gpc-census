"""Vertex attainment is a DIAGONAL-BASIS statement, and the checker must say so.

A spectrum-level check cannot establish attainment: det(rho - xI) is invariant
under conjugation, so it is satisfied by every state in the U(d) orbit,
including states whose 1-RDM is rotated away from the diagonal. These tests pin
the acceptance criterion actually used for the census convention:

    rho exactly diagonal, diag(rho) exactly equal to lambda, norm exactly 1.

The headline case is the support-10 state published for v103 in an earlier
revision of this branch. It satisfies normalization, hermiticity, trace and the
characteristic-polynomial identity exactly, and it is NOT a vertex attainer.
This test expects it to be REJECTED, so the weaker criterion cannot come back.
"""
import pathlib
import sys
from fractions import Fraction

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

D = 10
# the state published in results/data/v103_support10_certificate.json
SUPPORT = [(1, 2, 3), (0, 2, 7), (0, 5, 6), (0, 1, 8), (0, 3, 4),
           (0, 4, 9), (3, 8, 9), (1, 6, 9), (3, 5, 9), (0, 1, 5)]
WEIGHTS = [sp.Rational(13, 34), sp.Rational(5, 34), sp.Rational(53, 442),
           sp.Rational(195, 1802), sp.Rational(5, 68), sp.Rational(5, 68),
           sp.Rational(35, 901), sp.Rational(1, 34), sp.Rational(18, 901),
           sp.Rational(84, 11713)]
SIGNS = [1, 1, 1, -1, 1, 1, 1, 1, 1, 1]
LAMBDA = [sp.Rational(18, 34)] * 4 + [sp.Rational(5, 34)] * 6


def _rho_exact(weights, signs, support, d):
    amps = [s * sp.sqrt(w) for s, w in zip(signs, weights)]
    pos = {t: i for i, t in enumerate(support)}
    rho = sp.zeros(d, d)
    for col, det in enumerate(support):
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
                rho[p, q] += amps[row] * amps[col] * sign_out * (-1) ** target.index(p)
    return sp.simplify(rho)


def test_published_support10_state_passes_the_weak_checks():
    """It really does satisfy every check the old pipeline applied."""
    rho = _rho_exact(WEIGHTS, SIGNS, SUPPORT, D)
    x = sp.Symbol("x")
    assert sum(WEIGHTS) == 1
    assert sp.simplify(sp.trace(rho) - 3) == 0
    assert sp.simplify(rho - rho.T) == sp.zeros(D, D)
    charpoly = sp.expand(rho.charpoly(x).as_expr())
    target = sp.expand(sp.prod([(x - lam) for lam in LAMBDA]))
    assert sp.simplify(charpoly - target) == 0


def test_published_support10_state_is_rejected_as_an_attainer():
    """...and is still NOT a vertex attainer: its 1-RDM is not diagonal.

    Regression guard for the defect that shipped: the characteristic-polynomial
    identity is conjugation invariant, so passing it says nothing about the
    basis. rho_{3,9} = -5/68 here, and the diagonal is strictly majorized by
    lambda rather than equal to it.
    """
    rho = _rho_exact(WEIGHTS, SIGNS, SUPPORT, D)
    offdiag = [sp.simplify(rho[i, j]) for i in range(D) for j in range(D) if i != j]
    assert any(e != 0 for e in offdiag), "expected a NON-diagonal 1-RDM"
    assert sp.nsimplify(rho[3, 9]) == sp.Rational(-5, 68)
    diag = [sp.simplify(rho[i, i]) for i in range(D)]
    assert diag != LAMBDA
    assert sorted(diag) != sorted(LAMBDA)


def test_diagonality_error_flags_the_published_state():
    """The engine's attainment metric must be large on it, not tolerance noise."""
    from gpc_census.cascade import diagonality_error

    amps = [float(s) * float(sp.sqrt(w)) for s, w in zip(SIGNS, WEIGHTS)]
    err = diagonality_error(amps, SUPPORT, LAMBDA, D)
    assert err > 1e-3, err
    assert abs(err - float(sp.Rational(5, 68))) < 1e-9 or err > float(sp.Rational(5, 68))


def test_library_states_that_claim_attainment_are_diagonal():
    """Every library state recorded as natural-orbital really is one.

    The census ships both kinds: 645 records whose 1-RDM is exactly diagonal
    with diagonal equal to lambda, and 154 interference records whose 1-RDM is
    genuinely non-diagonal. This test checks the flag in cascade.jsonl agrees
    with a from-scratch recomputation, so the two populations cannot be
    silently merged.
    """
    import json

    import pytest

    from gpc_census.cascade import diagonality_error
    from gpc_census.fiber import amplitudes_from_record

    casc = ROOT / "results" / "data" / "cascade.jsonl"
    states = ROOT / "results" / "data" / "states.jsonl"
    if not casc.exists() or not states.exists():
        pytest.skip("no cascade dataset present")
    flags = {}
    for line in casc.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            if "start_is_natural_orbital" in r:
                flags[(r["system"], r["index"])] = r["start_is_natural_orbital"]
    if not flags:
        pytest.skip("cascade.jsonl predates the diagonality fields")

    checked = 0
    for line in states.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        key = (rec["system"], rec["index"])
        if key not in flags or checked >= 40:
            continue
        d = int(rec["system"].strip("()").split(",")[1])
        den = int(rec["denominator"])
        spectrum = [Fraction(int(v), den) for v in rec["integer_form"]]
        c, dets = amplitudes_from_record(rec)
        err = diagonality_error(c, dets, spectrum, d)
        assert (err < 1e-9) == flags[key], (key, err, flags[key])
        checked += 1
    assert checked > 0
