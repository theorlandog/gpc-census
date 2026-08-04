"""Guards for the conditional 2-RDM compatibility body and its dual certificates.

The claim being pinned is that the observable-range solver is a support oracle
for a convex body, and that at carrier order three that body is exactly the
scaled complex elliptope. The gates are

1. the dual KKT elimination reproduces the shipped sextics, so the algebraic
   solver and the SDP are the same object seen twice,
2. dual PSD-ness selects the endpoints pointwise, with no appeal to having
   found every real root,
3. at order three the SDP equals the phase hull, and at order four it does not,
4. the exact rational primitives (Sturm, isolation, sign-at-root) behave.
"""
import json
import pathlib
from fractions import Fraction as F

import pytest

from gpc_census import elliptope as el

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "conditional_2rdm_body.json"
SEXTIC_ARTIFACT = ROOT / "complex_fixed_1rdm_observable_ranges.json"

W = [F(3, 5), F(1, 4), F(3, 20)]


def _body():
    return json.loads(ARTIFACT.read_text())


# --------------------------------------------------------------------------
# exact rational primitives
# --------------------------------------------------------------------------

def test_sturm_isolation_finds_every_real_root():
    # (t - 1)(t + 2)(t - 1/3) = t^3 + (2/3)t^2 - (7/3)t + 2/3
    p = [F(2, 3), F(-7, 3), F(2, 3), F(1)]
    got = el.isolate_real_roots(p)
    assert len(got) == 3
    # isolation only promises one root per interval; refine before comparing
    refined = [el.refine_root(p, lo, hi, F(1, 10 ** 12)) for lo, hi in got]
    approx = sorted(float((lo + hi) / 2) for lo, hi in refined)
    assert approx == pytest.approx([-2.0, 1 / 3, 1.0], abs=1e-9)
    # and the intervals are genuinely disjoint
    bounds = sorted(refined)
    assert all(bounds[i][1] <= bounds[i + 1][0] for i in range(len(bounds) - 1))


def test_sign_at_root_handles_a_minor_vanishing_at_the_root():
    p = [F(-1), F(0), F(1)]          # t^2 - 1, roots +-1
    lo, hi = F(0), F(2)
    # q = t - 1 vanishes exactly at the root t = 1
    assert el.sign_at_root([F(-1), F(1)], p, lo, hi)[0] == 0
    # q = t + 1 is strictly positive there
    assert el.sign_at_root([F(1), F(1)], p, lo, hi)[0] == 1


def test_refine_root_shrinks_the_interval():
    p = [F(-2), F(0), F(1)]           # t^2 - 2
    lo, hi = el.refine_root(p, F(1), F(2), width=F(1, 10 ** 12))
    assert hi - lo <= F(1, 10 ** 12)
    assert float(lo) < 2 ** 0.5 < float(hi)


# --------------------------------------------------------------------------
# the order-three theorem
# --------------------------------------------------------------------------

def test_extreme_rank_bound():
    assert el.extreme_rank_bound(3) == 1      # forces phase matrices
    assert el.extreme_rank_bound(4) == 2      # where the theorem stops
    assert el.extreme_rank_bound(9) == 3


def test_order_four_has_a_genuine_rank_two_extreme_point():
    np = pytest.importorskip("numpy")
    vs, R = el.rank_two_extreme_point_4x4()
    assert np.allclose(np.real(np.diag(R)), 1.0)
    assert np.linalg.eigvalsh(R).min() > -1e-9
    assert np.linalg.matrix_rank(R, tol=1e-9) == 2
    assert el.is_extreme_point(vs)


def test_sdp_is_exact_at_order_three_and_strict_above():
    body = _body()
    for name, rec in body["sdp_versus_phase_hull"].items():
        assert rec["exact_at_order_three"], name
        assert rec["gap"] < 1e-6, name
    b = body["order_four_boundary"]
    assert b["extreme_point_rank"] == 2
    assert b["is_extreme"]
    assert b["sdp_is_strict_outer_bound"]
    assert b["relaxation_gap"] > 1e-3


# --------------------------------------------------------------------------
# the dual certificates
# --------------------------------------------------------------------------

def test_kkt_elimination_reproduces_the_shipped_sextics():
    """The sextic is the KKT elimination ideal of the 3x3 elliptope SDP.

    This is what demotes the sextic from a bespoke trigonometric trick to the
    algebraic solution of a universal semidefinite program.
    """
    sp = pytest.importorskip("sympy")
    if not SEXTIC_ARTIFACT.exists():
        pytest.skip("sextic artifact not present")
    shipped = json.loads(SEXTIC_ARTIFACT.read_text())
    for reg in shipped["regressions"]:
        h = [[sp.sympify(x["re"]) + sp.I * sp.sympify(x["im"]) for x in row]
             for row in reg["matrix"]]
        cert = el.exact_certificates(h, W)
        want = [int(c) for c in
                reg["energy_critical_polynomial_coefficients_descending"]][::-1]
        assert cert["objective_polynomial_ascending"] == want, reg["name"]


def test_dual_certificates_agree_with_the_sextic_endpoints():
    if not SEXTIC_ARTIFACT.exists():
        pytest.skip("sextic artifact not present")
    sp = pytest.importorskip("sympy")
    shipped = json.loads(SEXTIC_ARTIFACT.read_text())
    for reg in shipped["regressions"]:
        h = [[sp.sympify(x["re"]) + sp.I * sp.sympify(x["im"]) for x in row]
             for row in reg["matrix"]]
        cert = el.exact_certificates(h, W)
        er = reg["exact_range"]
        for tag, key in (("maximum", "maximum_isolating_interval"),
                         ("minimum", "minimum_isolating_interval")):
            lo, hi = F(er[key]["lower"]), F(er[key]["upper"])
            got = F(cert[tag]["interval"]["lower"]), F(cert[tag]["interval"]["upper"])
            # the two isolating intervals must overlap: same real number
            assert got[1] >= lo and got[0] <= hi, (reg["name"], tag)


def test_every_certificate_has_a_valid_psd_witness():
    body = _body()
    for name, rec in body["certificates"].items():
        mx, mn = rec["maximum"], rec["minimum"]
        # maximum: Diag(y) - H PSD, so every principal minor is nonnegative
        assert all(s >= 0 for s in mx["diagonal_minor_signs"]), name
        assert all(s >= 0 for s in mx["two_by_two_minor_signs"]), name
        # minimum: H - Diag(y) PSD, which flips the odd-order minors only
        assert all(s <= 0 for s in mn["diagonal_minor_signs"]), name
        assert all(s >= 0 for s in mn["two_by_two_minor_signs"]), name
        assert mx["approx"] > mn["approx"], name
        # exactly one root certifies each endpoint
        tags = [r["certifies"] for r in rec["all_real_roots"]]
        assert tags.count("maximum") == 1, name
        assert tags.count("minimum") == 1, name


def test_degenerate_flux_is_covered_by_the_dual_route():
    """Zero and pi flux have Im(chi) = 0, which the sextic engine excludes.

    The dual certificate handles them uniformly, so this records a genuine
    capability gain rather than a restatement.
    """
    body = _body()
    degenerate = [n for n, r in body["certificates"].items() if r["degenerate_flux"]]
    assert set(degenerate) == {"zero_flux", "pi_flux"}
    for name in degenerate:
        rec = body["certificates"][name]
        assert rec["maximum"]["certifies"] == "maximum"
        assert rec["minimum"]["certifies"] == "minimum"


def test_real_symmetric_endpoints_match_the_known_closed_forms():
    sp = pytest.importorskip("sympy")
    body = _body()
    v = float(sp.Rational(3, 10) * (2 + sp.sqrt(15)))
    assert body["certificates"]["zero_flux"]["minimum"]["approx"] == pytest.approx(-1.0, abs=1e-9)
    assert body["certificates"]["zero_flux"]["maximum"]["approx"] == pytest.approx(v, abs=1e-9)
    assert body["certificates"]["pi_flux"]["minimum"]["approx"] == pytest.approx(-v, abs=1e-9)
    assert body["certificates"]["pi_flux"]["maximum"]["approx"] == pytest.approx(1.0, abs=1e-9)


def test_hidden_width_is_positive_and_symmetric_under_flux_reversal():
    body = _body()
    c = body["certificates"]
    assert all(r["hidden_width_approx"] > 0 for r in c.values())
    # conjugation symmetry: reversing the flux mirrors the range
    assert c["zero_flux"]["hidden_width_approx"] == pytest.approx(
        c["pi_flux"]["hidden_width_approx"], abs=1e-9)
    assert c["zero_flux"]["maximum"]["approx"] == pytest.approx(
        -c["pi_flux"]["minimum"]["approx"], abs=1e-9)
