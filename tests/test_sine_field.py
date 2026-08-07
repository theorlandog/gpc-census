"""Guards for Theorem F and the sine field.

The claim is that Theorem E's discriminant is a product of triangle areas, so
its rationality condition is membership in a common SINE field. The tests that
matter are the ones that could catch the invariant being vacuous or the field
test being broken:

- the field-membership helper must be able to FAIL, since the first run of the
  script reported 0 of 44 purely because a degree was read off the wrong symbol;
- the invariant must be tested on pairs with irrational t as well, or most
  classes contribute a single pair and the measurement says almost nothing;
- the cross-check against holonomy_fields.json sin_radicands is what ties the
  object to independently generated data.
"""
import json
import pathlib

import pytest
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "sine_field.json"


def _art():
    return json.loads(ARTIFACT.read_text())


def _pairs(art):
    return [(b, p) for b in art["classes"] for p in b["pairs"]]


# --------------------------------------------------------------------------
# Theorem F, the identities
# --------------------------------------------------------------------------

def test_discriminant_is_a_product_of_areas():
    art = _art()
    assert art["summary"]["triples"] == 60
    assert art["summary"]["delta_equals_area_product"] == 60
    for block in art["classes"]:
        for triple in block["triples"]:
            assert triple["delta_equals_area_product"], (
                block["system"], triple["triple"])


def test_sine_addition_identity_holds_everywhere():
    art = _art()
    assert art["summary"]["sine_addition_holds"] == 60
    for block in art["classes"]:
        for triple in block["triples"]:
            assert triple["sine_addition_holds"], (
                block["system"], triple["triple"])


def test_the_pythagorean_relation_is_exact_pair_by_pair():
    """t^2 + A^2 = r_a^2 r_b^2, recomputed rather than re-read."""
    source = json.loads((ROOT / "shared_loop_rationality.json").read_text())
    art = _art()
    by_key = {(b["system"], b["index"], tuple(b["channel"])): b
              for b in source["classes"]}
    checked = 0
    for block in art["classes"]:
        src = by_key[(block["system"], block["index"], tuple(block["channel"]))]
        r2 = [sp.nsimplify(sp.sympify(v)) for v in src["r_squared"]]
        for entry in block["pairs"]:
            a, b = entry["pair"]
            t = sp.nsimplify(sp.sympify(entry["t_ab"]))
            area2 = sp.nsimplify(sp.sympify(entry["area_squared"]))
            assert sp.simplify(t**2 + area2 - r2[a] * r2[b]) == 0, (
                block["system"], entry["pair"])
            checked += 1
    assert checked == 54


# --------------------------------------------------------------------------
# the invariant
# --------------------------------------------------------------------------

def test_the_sine_field_is_uniform_on_every_class():
    art = _art()
    assert art["summary"]["classes"] == 17
    assert art["summary"]["sine_field_uniform"] == 17
    assert art["summary"]["transport_classes"] == 11
    for block in art["classes"]:
        assert block["sine_field_is_uniform"], (block["system"], block["index"])


def test_membership_is_tested_on_irrational_t_pairs_too():
    """Otherwise most classes contribute one pair and the test is near vacuous.

    26 of the 54 pairs have rational t. The invariant is only interesting if
    the other 28 are tested as well.
    """
    art = _art()
    tested = [p for _, p in _pairs(art) if p.get("in_sine_field") is not None]
    assert art["summary"]["in_sine_field_tested"] == 44
    assert art["summary"]["in_sine_field_holds"] == 44
    irrational = [p for p in tested if not p["t_rational"]]
    assert len(irrational) >= 20, "the invariant must be probed off the rational pairs"
    for entry in tested:
        assert entry["in_sine_field"] is True


def test_the_field_membership_helper_can_fail():
    """The first run reported 0 of 44 because a degree was read off the wrong
    symbol. A test that only ever sees True would not have caught it, so the
    helper is exercised on a value that is genuinely NOT in the base field.
    """
    import importlib.util

    path = (pathlib.Path(__file__).resolve().parents[1]
            / "scripts" / "sine_field.py")
    if not path.exists():
        pytest.skip("script not shipped")
    spec = importlib.util.spec_from_file_location("sine_field", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module.in_field_of(sp.Rational(3, 4), sp.Rational(1, 2)) is True
    assert module.in_field_of(sp.sqrt(2), sp.Rational(1, 2)) is False
    assert module.in_field_of(sp.sqrt(2), sp.sqrt(2)) is True
    assert module.in_field_of(sp.sqrt(3), sp.sqrt(2)) is False


def test_no_nonzero_area_is_rational():
    """When the sine field is nondegenerate it is a genuine quadratic field."""
    art = _art()
    nonzero = [p for _, p in _pairs(art) if not p["area_vanishes"]]
    assert len(nonzero) == 44
    assert not [p for p in nonzero if p["area_rational"]]


# --------------------------------------------------------------------------
# the cross-check against independently generated data
# --------------------------------------------------------------------------

def test_squarefree_parts_are_recorded_sin_radicands():
    """holonomy_fields.json generates sin_radicands for a different purpose."""
    art = _art()
    assert art["summary"]["matches_recorded_sin_radicands"] == 17
    for block in art["classes"]:
        assert block["matches_recorded_sin_radicands"], block["system"]
        if block["sine_radicand_primes"]:
            assert set(block["sine_radicand_primes"]) <= set(
                block["recorded_sin_radicand_primes"])


def test_the_invisible_sine_field_is_v75_and_is_the_non_degenerate_class():
    """The finding the dissection was for.

    (3,10) v75 has every loop of degree 1, so the holonomy census records NO
    cosine radicands at all, yet its sines carry sqrt(15). It is also the class
    that produced the single non-degenerate Theorem E triple.
    """
    art = _art()
    differing = [b for b in art["classes"]
                 if b["sine_field_differs_from_cosine_field"]]
    assert len(differing) == 1
    block = differing[0]
    assert (block["system"], block["index"], tuple(block["channel"])) == (
        "(3,10)", 75, (2, 8))
    assert block["sine_field_squarefree"] == 15
    assert block["recorded_cos_radicand_primes"] == []

    closure = json.loads((ROOT / "class_closure.json").read_text())
    non_degenerate = [b for b in closure["classes"]
                      for x in b["theorem_E_triples"] if not x["discriminant_zero"]]
    assert len(non_degenerate) == 1
    assert non_degenerate[0]["index"] == 75


def test_degenerate_classes_are_the_real_up_to_gauge_state():
    """17 of 18 triples degenerate because v103 is real up to gauge."""
    art = _art()
    assert art["summary"]["degenerate_classes"] == 2
    degenerate = {(b["system"], b["index"]) for b in art["classes"]
                  if b["sine_field_squarefree"] == 0}
    assert degenerate == {("(3,10)", 103)}


@pytest.mark.parametrize("phrase", ["weakly", "sine field"])
def test_the_write_up_records_the_weak_pass(phrase):
    doc = (pathlib.Path(__file__).resolve().parents[1] / "docs" / "sine_field.md")
    if not doc.exists():
        pytest.skip("docs not present")
    assert phrase.lower() in doc.read_text().lower()
