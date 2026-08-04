import importlib.util
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "fixed_1rdm_observable_ranges.json"

MODULE_PATH = ROOT / "scripts" / "verify_fixed_1rdm_observable_ranges.py"
SPEC = importlib.util.spec_from_file_location("verify_fixed_1rdm", MODULE_PATH)
verify = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verify)


def artifact():
    return json.loads(ARTIFACT.read_text())


def diagonal():
    return [sp.Rational(x) for x in artifact()["target_1rdm"]["diagonal"]]


def test_target_is_a_pinned_borland_dennis_point():
    lam = diagonal()
    assert sum(lam) == 3
    assert lam[0] + lam[5] == lam[1] + lam[4] == lam[2] + lam[3] == 1
    assert 2 - (lam[0] + lam[1] + lam[3]) == 0
    assert all(lam[i] > lam[i + 1] for i in range(5))


def test_selection_rules_force_the_three_determinant_carrier():
    """The carrier is derived, not assumed, so derive it here too."""
    determinants, single, kernel, carrier = verify.selection_rule_carrier()
    assert len(determinants) == 20
    assert len(single) == 8
    assert len(kernel) == 9
    assert [list(D) for D in carrier] == [[0, 1, 2], [0, 3, 4], [1, 3, 5]]
    assert [list(D) for D in carrier] == artifact()["fiber"]["support_dets"]


def test_weights_are_uniquely_forced_by_the_diagonal():
    _determinants, _single, _kernel, carrier = verify.selection_rule_carrier()
    weights = verify.weights_from_diagonal(carrier, diagonal())
    assert weights is not None
    assert [str(w) for w in weights] == ["3/5", "1/4", "3/20"]
    assert sum(weights) == 1


def test_carrier_is_one_hop_free():
    """One-hop-freeness is what makes the 1-RDM diagonal for every phase."""
    _determinants, _single, _kernel, carrier = verify.selection_rule_carrier()
    assert verify.one_hop_pairs(carrier) == []


def test_phasor_identity_holds_because_the_radii_are_normalised():
    _determinants, _single, _kernel, carrier = verify.selection_rule_carrier()
    weights = verify.weights_from_diagonal(carrier, diagonal())
    A, B, C, radii = verify.couplings(weights)
    assert sp.simplify(sum(r**2 for r in radii)) == 1
    alpha, beta = sp.symbols("alpha beta", real=True)
    real_part = radii[0] + radii[1] * sp.cos(alpha) + radii[2] * sp.cos(beta)
    imag_part = radii[1] * sp.sin(alpha) + radii[2] * sp.sin(beta)
    phasor = sp.expand_trig(sp.expand(real_part**2 + imag_part**2 - 1))
    normal = sp.expand_trig(sp.expand(
        A * sp.cos(alpha) + B * sp.cos(beta) + C * sp.cos(beta - alpha)
    ))
    assert sp.simplify(phasor - normal) == 0


def test_both_ranges_recompute_by_the_branch_rule():
    """Compared symbolically: the branch rule returns a nested radical.

    sqrt(12*sqrt(15) + 51) denests to 6 + sqrt(15), so the recorded
    3(2 + sqrt(15))/10 is the same number in a tidier form.
    """
    _determinants, _single, _kernel, carrier = verify.selection_rule_carrier()
    weights = verify.weights_from_diagonal(carrier, diagonal())
    A, B, C, _radii = verify.couplings(weights)
    for name, block in artifact()["examples"].items():
        matrix = [[sp.sympify(v) for v in row] for row in block["range"]["carrier_matrix"]]
        low, high, _ = verify.exact_range(A * matrix[0][1], B * matrix[0][2], C * matrix[1][2])
        assert sp.simplify(low - sp.sympify(block["range"]["minimum"])) == 0, name
        assert sp.simplify(high - sp.sympify(block["range"]["maximum"])) == 0, name


def test_nested_radical_denests_as_claimed():
    assert sp.simplify(sp.sqrt(12 * sp.sqrt(15) + 51) - (6 + sp.sqrt(15))) == 0


def test_polygon_certificate_agrees_with_the_branch_rule():
    """Two independent routes to the same interval, not one route twice."""
    _determinants, _single, _kernel, carrier = verify.selection_rule_carrier()
    weights = verify.weights_from_diagonal(carrier, diagonal())
    _A, _B, _C, radii = verify.couplings(weights)
    low, high, closes = verify.polygon_range(radii)
    assert closes
    branch = artifact()["examples"]["unfrustrated_coherence_triangle"]["range"]
    assert sp.simplify(low - sp.sympify(branch["minimum"])) == 0
    assert sp.simplify(high - sp.sympify(branch["maximum"])) == 0


def test_every_recorded_phase_triple_is_realizable():
    """Cosines that do not come from real angles would be a fabricated witness."""
    for block in artifact()["examples"].values():
        for which in ("minimum_phase_cosines", "maximum_phase_cosines"):
            cosines = [
                sp.sympify(block["range"][which][k])
                for k in ("cos_alpha", "cos_beta", "cos_beta_minus_alpha")
            ]
            assert all(abs(sp.N(x)) <= 1 for x in cosines)
            assert verify.gram_defect(*cosines) == 0


def test_pi_flux_is_the_negation_of_the_unfrustrated_case():
    """Shifting both phases by pi flips the sign, so the ranges mirror.

    This is why the frustrated example is a real regression case rather than a
    rescaling: the interior stationary point moves from the minimum branch to
    the maximum branch.
    """
    examples = artifact()["examples"]
    plain = examples["unfrustrated_coherence_triangle"]["range"]
    flux = examples["pi_flux_coherence_triangle"]["range"]
    assert sp.simplify(sp.sympify(flux["maximum"]) + sp.sympify(plain["minimum"])) == 0
    assert sp.simplify(sp.sympify(flux["minimum"]) + sp.sympify(plain["maximum"])) == 0
    assert plain["minimum_source"] == "interior stationary point"
    assert flux["maximum_source"] == "interior stationary point"


def test_cycle_products_have_opposite_signs():
    examples = artifact()["examples"]
    def cycle(name):
        m = examples[name]["range"]["carrier_matrix"]
        return sp.sympify(m[0][1]) * sp.sympify(m[0][2]) * sp.sympify(m[1][2])
    assert cycle("unfrustrated_coherence_triangle") > 0
    assert cycle("pi_flux_coherence_triangle") < 0


def test_generator_reproduces_the_shipped_artifact():
    """The tripwire this replaces asserted the generator was absent.

    It has since been delivered, so the check becomes the stronger one: run it
    and require the shipped artifact back unchanged.
    """
    import subprocess

    before = ARTIFACT.read_bytes()
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "fixed_1rdm_observable_ranges.py")],
        check=True,
        capture_output=True,
    )
    assert ARTIFACT.read_bytes() == before
