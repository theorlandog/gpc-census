"""Guards for Theorem E and the 17 of 17 closure.

The claim is that every large one-hop class has multiquadratic within-class
holonomy cosines, grounded in Theorem A alone. The tests that matter are the
ones that could catch the chain quietly consuming a measured rationality, or
the quadratic being circular:

- Theorem E's quadratic is rebuilt here from t_ab and t_bc and the measured
  t_ac substituted afterwards, so a certificate echoing its own input fails;
- the last class must close through a triple flagged as using ONLY
  Theorem-D-proved inputs, or the closure leans on something checked;
- the pre-registered failure must stay failed.
"""
import json
import pathlib

import pytest
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "class_closure.json"


def _art():
    return json.loads(ARTIFACT.read_text())


def _triples(art):
    return [(b, x) for b in art["classes"] for x in b["theorem_E_triples"]]


# --------------------------------------------------------------------------
# Theorem E, rebuilt rather than re-read
# --------------------------------------------------------------------------

def test_theorem_E_quadratic_is_rational_and_the_measured_root_solves_it():
    art = _art()
    pairs = _triples(art)
    assert len(pairs) == 18
    for block, x in pairs:
        A, B, C = (sp.nsimplify(sp.sympify(v)) for v in x["coefficients"])
        assert A.is_rational and B.is_rational and C.is_rational, block["system"]
        measured = sp.nsimplify(sp.sympify(x["measured_t_ac"]))
        assert sp.simplify(A * measured**2 + B * measured + C) == 0, (
            block["system"], block["index"], x["triple"])
        assert x["measured_solves_quadratic"]


def test_discriminant_factorisation_is_an_identity():
    """Delta/4 = (t_ab^2 - r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2), the geometric form."""
    art = _art()
    for block, x in _triples(art):
        assert x["discriminant_factored_matches"], (block["system"], x["triple"])
        assert x["discriminant_rational"]


def test_every_discriminant_is_a_rational_square():
    """Which is why every t_ac in these classes comes out rational."""
    art = _art()
    for block, x in _triples(art):
        delta = sp.nsimplify(sp.sympify(x["discriminant"]))
        assert sp.sqrt(delta).is_rational, (block["system"], x["triple"])
        assert x["discriminant_is_a_rational_square"]


# --------------------------------------------------------------------------
# the closure, and that it does not lean on anything checked
# --------------------------------------------------------------------------

def test_all_seventeen_classes_close():
    art = _art()
    summary = art["summary"]
    assert summary["classes"] == 17
    assert summary["closed"] == 17
    assert summary["not_closed"] == []
    assert summary["closed_by_theorem_D"] == 16
    assert summary["closed_by_theorem_E"] == 1


def test_the_last_class_closes_through_theorem_D_inputs_only():
    """The chain must be grounded in Theorem A, not in a measured rationality.

    (3,10) v103 ch(1,6) closes by transport from ch(3,9), and the triple that
    makes ch(3,9) supply it must be flagged as using only Theorem-D-proved
    entries. Without that flag the closure would be circular.
    """
    art = _art()
    donor = next(b for b in art["classes"]
                 if b["index"] == 103 and b["channel"] == [3, 9])
    last = next(b for b in art["classes"]
                if b["index"] == 103 and b["channel"] == [1, 6])

    grounded = [x for x in donor["theorem_E_triples"]
                if x["uses_only_theorem_D_inputs"] and x["forces_rational_t_ac"]]
    assert grounded, "no Theorem-E triple built from Theorem-D inputs alone"
    assert donor["theorem_E_extends_rationality"]
    assert last["closed"]
    assert last["closed_by_transport_from"] == [3, 9]
    assert not last["closed_by_theorem_D"], (
        "the point of this class is that Theorem D cannot reach it")


def test_the_last_class_closes_because_both_sines_vanish():
    """v103 is real up to gauge, so Delta is exactly 0, not merely a square."""
    art = _art()
    donor = next(b for b in art["classes"]
                 if b["index"] == 103 and b["channel"] == [3, 9])
    grounded = next(x for x in donor["theorem_E_triples"]
                    if x["uses_only_theorem_D_inputs"])
    assert grounded["discriminant_zero"]
    assert grounded["sin_ab_vanishes"] and grounded["sin_bc_vanishes"]
    assert sp.nsimplify(sp.sympify(grounded["discriminant"])) == 0


# --------------------------------------------------------------------------
# the pre-registered failure
# --------------------------------------------------------------------------

def test_degeneracy_is_the_rule_not_a_v103_peculiarity():
    """P-E-6 FAILED, as pre-registered.

    17 of 18 chained triples are degenerate, and they are not confined to
    v103. The write-up leans on this, so a change making degeneracy exclusive
    would silently convert a recorded refutation into a confirmation.
    """
    art = _art()
    assert art["summary"]["theorem_E_degenerate_triples"] == 17
    hosts = {(b["system"], b["index"]) for b, x in _triples(art)
             if x["discriminant_zero"]}
    assert ("(3,10)", 75) in hosts, "degeneracy is not exclusive to v103"
    assert ("(3,10)", 103) in hosts


def test_the_single_non_degenerate_triple_is_named():
    """Theorem E's other mode: a genuine rational square, not a vanishing sine."""
    art = _art()
    non_degenerate = [(b, x) for b, x in _triples(art)
                      if not x["discriminant_zero"]]
    assert len(non_degenerate) == 1
    block, x = non_degenerate[0]
    assert (block["system"], block["index"], tuple(block["channel"])) == (
        "(3,10)", 75, (2, 8))
    delta = sp.nsimplify(sp.sympify(x["discriminant"]))
    assert delta != 0
    assert sp.sqrt(delta).is_rational


@pytest.mark.parametrize("phrase", ["real up to gauge", "17 of 17"])
def test_the_write_up_states_the_mechanism(phrase):
    doc = (pathlib.Path(__file__).resolve().parents[1]
           / "docs" / "class_closure.md")
    if not doc.exists():
        pytest.skip("docs not present")
    assert phrase.lower() in doc.read_text().lower()
