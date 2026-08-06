"""Guards for Theorem C and the polygon diagonal scorecard.

The claim is that a one-hop class of size 3 with a rational diagonal has all
its within-class holonomy cosines in a MULTIQUADRATIC field, so (R1) holds and
C4 and D4 are excluded. The tests that matter are the ones that could catch it
being circular, vacuous, or an artifact of a failed symbolic simplification:

- the Gram identity must be re-derived live from rational data and checked
  against an independently measured t_ac, not re-read from the artifact;
- a degree that sympy failed to resolve must never be recorded as irrational,
  which is the bug the first run of this script actually had;
- the failed predictions must stay failed, since the write-up leans on them.
"""
import json
import pathlib

import pytest
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "polygon_diagonals.json"


def _art():
    return json.loads(ARTIFACT.read_text())


# --------------------------------------------------------------------------
# scope, re-derived rather than quoted
# --------------------------------------------------------------------------

def test_scope_is_the_seventeen_large_classes():
    art = _art()
    summary = art["summary"]
    assert summary["classes"] == 17
    assert summary["class_sizes"] == {"3": 16, "4": 1}
    # rule R2: the independent content is transport classes, not records
    assert summary["distinct_states"] == 16
    assert summary["transport_classes"] == 11
    assert summary["transport_classes"] < summary["classes"]


def test_no_degree_was_left_unresolved():
    """The first run recorded 29 unresolved degrees AS irrational.

    ``rational`` is tri-state for exactly this reason. If a value's degree
    cannot be computed the entry must say so rather than defaulting to a
    verdict, or the scorecard silently invents irrationality.
    """
    art = _art()
    assert art["summary"]["degree_unresolved"] == 0
    for block in art["classes"]:
        for subset in block["proper_subsets"]:
            assert subset["degree_resolved"], (block["system"], subset["subset"])
            assert subset["rational"] is not None
            assert subset["rational"] == (subset["degree_over_Q"] == 1)


# --------------------------------------------------------------------------
# Theorem C, re-derived live
# --------------------------------------------------------------------------

def test_gram_identity_is_exact_on_every_class():
    art = _art()
    summary = art["summary"]
    assert summary["theorem_C_applies"] == 17
    assert summary["theorem_C_gram_identity_holds"] == 17
    assert summary["theorem_C_all_coefficients_rational"] == 17
    assert summary["theorem_C_discriminant_rational"] == 17
    for block in art["classes"]:
        certificate = block["theorem_C"]
        assert certificate["applies"], block["system"]
        assert certificate["gram_identity_holds"], block["system"]


def test_gram_quadratic_replays_from_rational_data_alone():
    """Rebuild the quadratic from D^2, p, q and check the measured t_ac solves it.

    This is the non-circular version: the coefficients are recomputed here from
    the recorded rational scalars, and the root is the separately recorded
    measurement, so a certificate that merely echoed its own inputs would fail.
    """
    art = _art()
    checked = 0
    for block in art["classes"]:
        certificate = block["theorem_C"]
        if certificate.get("route") != "gram_triangulation":
            continue
        d2 = sp.Rational(certificate["D_squared"])
        p = sp.Rational(certificate["p"])
        q = sp.Rational(certificate["q"])
        a, b = certificate["diagonal_pair"]
        c = certificate["third"]
        r2 = [sp.Rational(v) for v in block["r_squared"]]
        t_ac = sp.nsimplify(sp.sympify(certificate["t_ac"]))
        value = (d2 * t_ac**2 - 2 * p * q * t_ac
                 - (r2[a] * r2[c] * d2 - r2[a] * p**2 - r2[c] * q**2))
        assert sp.simplify(value) == 0, (block["system"], block["index"], (a, b, c))
        # and the discriminant really is rational, which is the whole point
        discriminant = sp.nsimplify(sp.sympify(certificate["discriminant"]))
        assert discriminant.is_rational, block["system"]
        checked += 1
    assert checked == 14


def test_R1_holds_on_every_within_class_holonomy():
    """deg_Q(cos^2 Phi) <= 2, which is the open condition Theorem C targets."""
    art = _art()
    degrees = [c["cos_squared_degree"]
               for b in art["classes"] for c in b["within_class_cosines"]]
    assert len(degrees) == 54
    assert set(degrees) <= {1, 2}
    assert degrees.count(1) == 26
    assert degrees.count(2) == 28


def test_the_reduction_agrees_pair_by_pair():
    """t_ab rational iff cos^2 rational, which the write-up proves in advance.

    A disagreement would mean the exact arithmetic is wrong, not that the
    mathematics is interesting.
    """
    art = _art()
    for block in art["classes"]:
        for cosine in block["within_class_cosines"]:
            assert cosine["t_ab_rational"] == cosine["cos_squared_rational"], (
                block["system"], block["index"], cosine["pair"])


# --------------------------------------------------------------------------
# the failed predictions must stay failed
# --------------------------------------------------------------------------

def test_the_two_rdm_identification_holds_exactly():
    """P-B1-3 PASSED, and the failure of P-B1-2 is only meaningful if it did."""
    art = _art()
    assert art["summary"]["two_rdm_checked"] == 25
    assert art["summary"]["two_rdm_matches"] == 25


def test_H2_does_not_extend_to_the_two_rdm():
    """P-B1-2 FAILED: six spectator-mode 2-RDM entries are irrational.

    The write-up leans on this, so a later change that made them all rational
    would silently turn a recorded refutation into a confirmation.
    """
    art = _art()
    assert art["summary"]["mode_coherent_irrational"] == 6
    offenders = sorted({
        (b["system"], b["index"]) for b in art["classes"]
        for s in b["proper_subsets"]
        if s["mode_coherent"] and s["rational"] is False
    })
    assert offenders == [("(5,10)", 256), ("(5,10)", 261), ("(5,10)", 268)]


def test_mode_coherence_is_not_the_right_invariant():
    """P-B1-4 FAILED too: rationality does not track mode-coherence at all."""
    art = _art()
    summary = art["summary"]
    assert summary["mode_incoherent_rational"] == 11
    assert summary["mode_coherent_irrational"] == 6


# --------------------------------------------------------------------------
# the hypothesis is real, not automatic
# --------------------------------------------------------------------------

def test_the_rational_diagonal_hypothesis_is_not_vacuous():
    """If most pairs were rational the theorem would say nothing.

    28 of 54 within-class pairs are irrational and 14 of 17 classes carry
    exactly one rational diagonal, so the hypothesis is a real condition that
    happens to hold rather than a restatement of the setting.
    """
    art = _art()
    assert art["summary"]["classes_with_a_rational_diagonal"] == 17
    exactly_one = 0
    for block in art["classes"]:
        rational = [s for s in block["proper_subsets"]
                    if s["size"] == 2 and s["rational"] is True]
        assert rational, (block["system"], block["index"])
        exactly_one += len(rational) == 1
    assert exactly_one == 14


@pytest.mark.parametrize("key", ["POST-HOC", "out-of-sample"])
def test_the_write_up_labels_its_post_hoc_pattern(key):
    """17 of 17 was not pre-registered and must never be quoted as if it were."""
    doc = (pathlib.Path(__file__).resolve().parents[1]
           / "docs" / "polygon_diagonals.md")
    if not doc.exists():          # docs are not shipped in the sdist
        pytest.skip("docs not present")
    assert key.lower() in doc.read_text().lower()
