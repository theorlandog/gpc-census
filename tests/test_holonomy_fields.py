"""The holonomy field artifact, re-derived independently of the scan script.

results/data/holonomy_fields.json is the certifying artifact for the
Holonomy Rationality work (docs/holonomy_rationality.md). These tests pin its
acceptance criteria, and recompute the headline values at (3,10) v17 from the
ledger with their own gauge and kernel bookkeeping, so the 82/82 result does
not rest on the script that produced it.
"""
import json
import pathlib
import sys

import pytest
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"

ARTIFACT = DATA / "holonomy_fields.json"
LEDGER = DATA / "states.jsonl"


@pytest.fixture(scope="module")
def artifact():
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def ledger():
    return [json.loads(ln) for ln in LEDGER.read_text().splitlines() if ln.strip()]


# ---------------------------------------------------------------- acceptance


def test_scan_covers_every_loopy_interference_state(artifact, ledger):
    """The artifact's population is exactly the loopy INTERFERENCE states."""
    expected = set()
    total_loops = 0
    for rec in ledger:
        if rec["classified"] != "INTERFERENCE":
            continue
        dets = [tuple(t) for t in rec["closed_form"]["support_dets"]]
        d = max(m for T in dets for m in T) + 1
        m = sp.Matrix([[1 if o in T else 0 for o in range(d)] for T in dets])
        loops = len(m.T.nullspace())
        if loops:
            expected.add((rec["system"], rec["index"]))
            total_loops += loops

    seen = {(r["system"], r["index"]) for r in artifact["loops"]}
    assert seen == expected
    assert len(artifact["loops"]) == total_loops
    assert artifact["summary"]["loopy_states"] == 63
    assert artifact["summary"]["loops"] == 82


def test_every_degree_is_1_2_or_4(artifact):
    assert artifact["summary"]["degree_counts"] == {"1": 37, "2": 32, "4": 13}
    assert all(r["degree"] in (1, 2, 4) for r in artifact["loops"])


def test_every_galois_group_is_two_elementary(artifact):
    assert artifact["summary"]["all_two_elementary"]
    assert artifact["summary"]["galois_group_counts"] == {"C1": 37, "C2": 32, "V4": 13}
    for row in artifact["loops"]:
        assert row["two_elementary"], (row["system"], row["index"], row["loop"])


def test_every_relation_survived_double_precision_reverification(artifact):
    """The precision protocol: no loop ships without a re-verified relation."""
    s = artifact["summary"]
    assert s["all_pslq_reverified"]
    assert s["all_pslq_confirm_exact"]
    for row in artifact["loops"]:
        pslq = row["pslq"]
        assert pslq is not None
        assert pslq["verify_digits"] >= 2 * pslq["working_digits"]
        assert pslq["coeffs_ascending"] == row["min_poly_coeffs_ascending"]
        assert pslq["residual_log10"] < -(pslq["working_digits"] + pslq["verify_digits"]) / 2


def test_kernel_bases_are_saturated(artifact):
    """A finite-index sublattice would silently drop holonomy relations."""
    assert artifact["summary"]["all_kernels_saturated"]


def test_theorem_hypotheses_hold_across_the_census(artifact):
    """(H2) and (H3) of docs/holonomy_rationality.md."""
    h = artifact["channel_hypotheses"]
    assert h["H2_holds"], h["failures"]
    assert h["H3_holds"], h["failures"]
    assert h["failures"] == []
    assert h["channels"] == 84
    assert h["lemma_2_identity"] == "84 of 84"


# ---------------------------------------------------------------- the theorem


def test_size_two_classes_have_monomial_radicands(artifact):
    """The prediction Proposition B makes, scored in the artifact.

    A non-monomial radicand can only come from a channel class of size three
    or more; a state whose classes are all size two cannot have one.
    """
    from gpc_census.states import one_hop_classes

    ledger = {
        (r["system"], r["index"]): r
        for r in (
            json.loads(ln) for ln in LEDGER.read_text().splitlines() if ln.strip()
        )
    }
    violations = []
    for row in artifact["loops"]:
        rec = ledger[(row["system"], row["index"])]
        dets = [tuple(t) for t in rec["closed_form"]["support_dets"]]
        largest = max(len(p) for p in one_hop_classes(dets).values())
        outside = any(not m["in_weight_group"] for m in row["radicand_weight_group"])
        if outside and largest <= 2:
            violations.append((row["system"], row["index"], row["loop"]))
    assert violations == []


def test_the_monomial_subclause_is_false_as_stated(artifact):
    """Pinned so the falsification is not quietly lost in a later regeneration."""
    assert artifact["summary"]["all_radicands_in_weight_group"] is False
    assert artifact["summary"]["loops_with_radicand_outside_weight_group"] == 15


# ---------------------------------------------------------------- v17, by hand


def test_v17_holonomies_recomputed_from_the_ledger(ledger):
    """cos = 3/4 and cos = (sqrt 2 - sqrt 42)/8, derived here from scratch."""
    rec = next(r for r in ledger if r["system"] == "(3,10)" and r["index"] == 17)
    cf = rec["closed_form"]
    dets = [tuple(t) for t in cf["support_dets"]]
    amps = [sp.sympify(p) for p in cf["pretty"]]
    d = max(m for T in dets for m in T) + 1
    incidence = sp.Matrix([[1 if o in T else 0 for o in range(d)] for T in dets])
    phases = [sp.arg(a) for a in amps]

    found = {}
    for vec in incidence.T.nullspace():
        vec = vec * sp.lcm([q.q for q in vec])
        vec = vec / sp.gcd([sp.Integer(q) for q in vec])
        k = [int(q) for q in vec]
        assert (sp.Matrix([k]) * incidence).is_zero_matrix  # it is a loop
        holonomy = sum(k[i] * phases[i] for i in range(len(k)))
        cos = sp.simplify(sp.expand_trig(sp.cos(holonomy)))
        found[sp.simplify(sp.minimal_polynomial(cos, sp.Symbol("x")))] = cos

    x = sp.Symbol("x")
    assert set(found) == {4 * x - 3, 64 * x**4 - 88 * x**2 + 25}
    rational = found[4 * x - 3]
    assert sp.simplify(rational - sp.Rational(3, 4)) == 0
    quartic = found[64 * x**4 - 88 * x**2 + 25]
    assert abs(sp.N(sp.Abs(quartic) - (sp.sqrt(42) - sp.sqrt(2)) / 8, 40)) < 1e-35


def test_v17_quartic_is_klein_four(ledger):
    from sympy.polys.numberfields.galoisgroups import galois_group

    x = sp.Symbol("x")
    group, _ = galois_group(sp.Poly(64 * x**4 - 88 * x**2 + 25, x))
    assert group.order() == 4
    assert group.is_abelian
    assert all(g.order() in (1, 2) for g in group.elements)  # V4, not C4


# ---------------------------------------------------------------- trisection


def test_the_trisected_tier_contains_no_cubic_algebra(artifact):
    audit = artifact["trisection_audit"]
    assert audit["states"] == 13
    assert audit["no_cubic_anywhere"]
    assert audit["all_holonomies_degree_1_2_4"]
    for rec in audit["records"]:
        assert not rec["any_degree_divisible_by_three"]
        assert rec["all_degrees_are_powers_of_two"]


def test_the_printed_trisected_phase_reduces(ledger):
    """(pi + 3*atan(t))/3 is pi/3 + atan(t): the division by three cancels."""
    t = sp.sqrt(7)
    written = (sp.pi + 3 * sp.atan(t)) / 3
    assert sp.simplify(written - (sp.pi / 3 + sp.atan(t))) == 0
    # and the amplitude built from it is multiquadratic, not cubic
    amp = sp.exp(sp.I * written)
    value = sp.expand(sp.re(amp)) + sp.I * sp.expand(sp.im(amp))
    degree = sp.minimal_polynomial(sp.simplify(value), sp.Symbol("x")).as_poly().degree()
    assert degree % 3 != 0
    assert degree & (degree - 1) == 0


# ---------------------------------------------------------------- the outlier


def test_exactly_one_ledger_record_lacks_rational_weights(ledger):
    """Hypothesis (H1) is not vacuous: v_B is the single exception."""
    exceptions = [
        (r["system"], r["index"])
        for r in ledger
        if not all(isinstance(w, int) for w in r["closed_form"]["weights"])
    ]
    assert exceptions == [("(4,9)", 65)]

    rec = next(r for r in ledger if (r["system"], r["index"]) == ("(4,9)", 65))
    cf = rec["closed_form"]
    assert "family" in cf  # it carries the wall family it is an endpoint of

    # the weights are quadratic irrationals in Q(sqrt 35), and they are the
    # ones the pretty amplitudes square to
    for weight, pretty in zip(cf["weights"], cf["pretty"]):
        w = sp.sympify(weight)
        assert sp.simplify(sp.Abs(sp.sympify(pretty)) ** 2 - w / cf["den"]) == 0
    assert any(sp.sqrt(35) in sp.sympify(w).atoms(sp.Pow) for w in cf["weights"])

    # t0 is a root of the quadratic wall polynomial, so the Galois orbit has
    # size 2 and the weights cannot be rational
    t = sp.Symbol("t")
    t0 = sp.Rational(7, 17) - 18 * sp.sqrt(35) / 85
    assert sp.simplify((85 * t**2 - 70 * t - 119).subs(t, t0)) == 0
