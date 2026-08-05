"""Guards for the 2-elementarity reduction.

The load-bearing claim is that two conditions on cos^2 DERIVE the certified
Galois groups. That is only worth anything if the derivation is exercised
against the certification rather than asserted alongside it, so the central
test compares the predicted group to the shipped `galois_group` field loop by
loop.

The classical even-quartic criterion is also tested on textbook polynomials
including a genuine D4 case, so a bug that always returns V4 cannot pass.
"""
import json
import pathlib

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "holonomy_two_elementary.json"
FIELDS = ROOT / "holonomy_fields.json"

import sys  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import holonomy_two_elementary as h2e  # noqa: E402

X = sp.Symbol("x")


def _art():
    return json.loads(ARTIFACT.read_text())


# --------------------------------------------------------------------------
# the classical criterion itself
# --------------------------------------------------------------------------

def test_even_quartic_criterion_separates_v4_c4_d4():
    # x^4 - 2x^2 - 1 is the standard D4 example, sqrt(1 + sqrt 2)
    assert h2e.even_quartic_group(sp.Poly(X**4 - 2*X**2 - 1, X)) == "D4"
    # x^4 - 10x^2 + 1 = min poly of sqrt2 + sqrt3, biquadratic, q = 1 a square
    assert h2e.even_quartic_group(sp.Poly(X**4 - 10*X**2 + 1, X)) == "V4"
    # x^4 + 5x^2 + 5 is cyclic (the 5th-root case): q(p^2-4q) = 5*5 = 25
    assert h2e.even_quartic_group(sp.Poly(X**4 + 5*X**2 + 5, X)) == "C4"
    # an odd polynomial is rejected rather than silently classified
    assert h2e.even_quartic_group(sp.Poly(X**4 + X**3 - 1, X)) == "not-even"


def test_criterion_is_scale_invariant():
    a = h2e.even_quartic_group(sp.Poly(X**4 - 10*X**2 + 1, X))
    b = h2e.even_quartic_group(sp.Poly(7*X**4 - 70*X**2 + 7, X))
    assert a == b == "V4"


# --------------------------------------------------------------------------
# the reduction, against the certified ledger
# --------------------------------------------------------------------------

def test_reduction_reproduces_every_certified_galois_group():
    art = _art()
    certified = {(r["system"], r["index"], r["loop"]): r["galois_group"]
                 for r in json.loads(FIELDS.read_text())["loops"]}
    assert len(art["loops"]) == len(certified) == 82
    for r in art["loops"]:
        key = (r["system"], r["index"], r["loop"])
        assert r["predicted_group"] == certified[key], key


def test_R1_holds_on_every_loop():
    art = _art()
    assert art["summary"]["R1_holds_on_all"]
    for r in art["loops"]:
        assert r["degree_cos_squared"] <= 2, (r["system"], r["index"])
        # R1 bounds the degree of the cosine itself at 4
        assert r["degree_cos"] <= 4


def test_R2_holds_on_every_quartic():
    art = _art()
    quartics = [r for r in art["loops"] if r["degree_cos"] == 4]
    assert len(quartics) == 13
    for r in quartics:
        assert r["R2_norm_is_square"], (r["system"], r["index"])
        assert r["min_poly_is_even"], (r["system"], r["index"])
        assert r["predicted_group"] == "V4"


def test_all_groups_are_two_elementary():
    art = _art()
    assert art["summary"]["all_two_elementary"]
    assert set(art["summary"]["predicted_groups"]) <= {"C1", "C2", "V4"}


# --------------------------------------------------------------------------
# Theorem A' and the boundary case
# --------------------------------------------------------------------------

def test_theorem_A_prime_has_no_violations():
    art = _art()
    t = art["theorem_A_prime_test"]
    assert t["holds"]
    assert t["violations"] == []
    # and the test has real content: some loops DO sit at class size <= 2
    small = [r for r in art["loops"] if r["max_one_hop_class_size"] == 2]
    assert small
    for r in small:
        assert r["cos_squared_rational"], (r["system"], r["index"])


def test_irrational_cos_squared_needs_a_large_class():
    """The contrapositive of Theorem A', which is what gives it bite."""
    art = _art()
    for r in art["loops"]:
        if not r["cos_squared_rational"]:
            assert r["max_one_hop_class_size"] >= 3, (r["system"], r["index"])


def test_the_two_non_even_loops_are_recorded_not_hidden():
    """'The min poly is always even' is FALSE; the exceptions must be listed."""
    art = _art()
    exc = art["summary"]["min_poly_even_except"]
    assert len(exc) == 2
    for e in exc:
        assert e["degree_cos"] == 2       # 2-elementarity survives trivially
    keys = {(e["system"], e["index"], e["loop"]) for e in exc}
    assert keys == {("(5,10)", 113, 1), ("(5,10)", 261, 1)}
    # they must still be classified, and correctly
    for r in art["loops"]:
        if (r["system"], r["index"], r["loop"]) in keys:
            assert r["predicted_group"] == "C2"


def test_theorem_A_prime_formula_on_a_worked_case():
    """cos = n / (2 sqrt M) with n, M rational forces cos^2 rational."""
    n, M = sp.Rational(3), sp.Rational(14)
    cos = n / (4 * sp.sqrt(M))                 # the v_B mechanism shape
    assert sp.expand(cos**2).is_rational
    mp = sp.Poly(sp.minimal_polynomial(cos, X), X)
    coeffs = sp.Poly(mp.as_expr() / mp.all_coeffs()[0], X).all_coeffs()
    assert all(v == 0 for v in coeffs[1::2])   # even
    assert mp.degree() == 2
