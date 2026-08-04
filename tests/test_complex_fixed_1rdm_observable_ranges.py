from __future__ import annotations

import importlib.util
import itertools
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "complex_fixed_1rdm_observable_ranges.py"
SPEC = importlib.util.spec_from_file_location("complex_ranges", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = mod
SPEC.loader.exec_module(mod)

ARTIFACT = ROOT / "results" / "data" / "complex_fixed_1rdm_observable_ranges.json"
STANDALONE = (
    ROOT / "scripts" / "verify_complex_fixed_1rdm_observable_ranges_standalone.py"
)


def test_quarter_flux_exact_sextic_and_range():
    result = mod.solve_complex_triangle(mod.unit_flux_matrix(sp.Rational(0), sp.Rational(1)))
    assert result.centered_critical_polynomial.all_coeffs() == [
        1_000_000,
        0,
        -2_330_000,
        0,
        -143_700,
        0,
        3_969,
    ]
    assert result.minimum.upper < sp.Rational(-154578, 100000)
    assert result.minimum.lower > sp.Rational(-154579, 100000)
    assert result.maximum.lower > sp.Rational(154578, 100000)
    assert result.maximum.upper < sp.Rational(154579, 100000)


def test_pythagorean_flux_exact_sextic_and_asymmetric_range():
    matrix = mod.unit_flux_matrix(sp.Rational(3, 5), sp.Rational(4, 5))
    result = mod.solve_complex_triangle(matrix)
    assert result.centered_critical_polynomial.all_coeffs() == [
        6_250_000,
        7_500_000,
        -12_065_000,
        -22_482_000,
        -12_517_575,
        -2_643_570,
        -136_971,
    ]
    assert result.minimum.upper < sp.Rational(-134617, 100000)
    assert result.maximum.lower > sp.Rational(168507, 100000)
    assert abs(float(result.minimum.upper)) != pytest.approx(float(result.maximum.lower))


def test_cycle_invariant_controls_gram_determinant():
    matrix = mod.unit_flux_matrix(sp.Rational(3, 5), sp.Rational(4, 5))
    invariants = mod.triangle_invariants(matrix)
    assert mod.cycle_product(matrix) == sp.Rational(3, 5) + 4 * sp.I / 5
    assert invariants.amplitude_norm == sp.Rational(3, 5)
    assert invariants.branch_direction_norm == sp.Rational(27, 125)
    assert invariants.mixed_product == sp.Rational(27, 125)
    assert invariants.signed_cycle_area == sp.Rational(36, 125)
    assert invariants.gram_determinant == sp.Rational(1296, 15625)


def test_unit_flux_family_bridges_real_symmetric_limits():
    plus = mod.unit_flux_critical_polynomial(sp.Rational(1)).as_expr()
    minus = mod.unit_flux_critical_polynomial(sp.Rational(-1)).as_expr()
    assert sp.factor(plus.subs(mod.z, -1)) == 0
    assert sp.factor(sp.diff(plus, mod.z).subs(mod.z, -1)) == 0
    assert sp.factor(minus.subs(mod.z, 1)) == 0
    assert sp.factor(sp.diff(minus, mod.z).subs(mod.z, 1)) == 0
    assert sp.factor(
        mod.unit_flux_critical_polynomial(-mod.u).as_expr().subs(mod.z, -mod.z)
        - mod.unit_flux_critical_polynomial(mod.u).as_expr()
    ) == 0


def test_gauge_rephasing_leaves_cycle_and_range_polynomial_unchanged():
    matrix = mod.unit_flux_matrix(sp.Rational(3, 5), sp.Rational(4, 5))
    phases = [
        sp.Rational(3, 5) + 4 * sp.I / 5,
        sp.Rational(5, 13) + 12 * sp.I / 13,
        sp.Rational(7, 25) + 24 * sp.I / 25,
    ]
    unitary = sp.diag(*phases)
    transformed = sp.simplify(sp.conjugate(unitary.T) * matrix * unitary)
    original_result = mod.solve_complex_triangle(matrix)
    transformed_result = mod.solve_complex_triangle(transformed)
    assert mod.cycle_product(transformed) == mod.cycle_product(matrix)
    assert transformed_result.invariants == original_result.invariants
    assert (
        transformed_result.centered_critical_polynomial
        == original_result.centered_critical_polynomial
    )


def test_rational_diagonal_only_translates_the_exact_range():
    base = mod.unit_flux_matrix(sp.Rational(0), sp.Rational(1))
    diagonal = [sp.Rational(2), sp.Rational(-1), sp.Rational(3)]
    shifted = base + sp.diag(*diagonal)
    result = mod.solve_complex_triangle(shifted)
    offset = sum(weight * value for weight, value in zip(mod.WEIGHTS, diagonal))
    assert result.invariants.offset == offset
    assert result.minimum.lower - offset < 0 < result.maximum.upper - offset
    shifted_back = sp.Poly(
        sp.expand(
            result.energy_critical_polynomial.as_expr().subs(mod.z, mod.z + offset)
        ),
        mod.z,
    ).monic()
    assert shifted_back == result.centered_critical_polynomial.monic()


def test_linear_lift_is_exact_modulo_every_critical_root():
    result = mod.solve_complex_triangle(
        mod.unit_flux_matrix(sp.Rational(3, 5), sp.Rational(4, 5))
    )
    curve = mod.branch_curve(result.invariants)
    lift = -result.lift_numerator.as_expr() / result.lift_denominator.as_expr()
    curve_numerator = sp.together(curve.subs(mod.y, lift)).as_numer_denom()[0]
    critical_numerator = sp.together(
        sp.diff(curve, mod.y).subs(mod.y, lift)
    ).as_numer_denom()[0]
    assert sp.rem(curve_numerator, result.centered_critical_polynomial.as_expr(), mod.z) == 0
    assert (
        sp.rem(
            critical_numerator,
            result.centered_critical_polynomial.as_expr(),
            mod.z,
        )
        == 0
    )
    assert sp.gcd(result.lift_denominator, result.centered_critical_polynomial).degree() == 0


def test_branch_elimination_matches_direct_expectation_on_a_grid():
    matrix = np.array(
        [
            [0, 1, 1],
            [1, 0, 3 / 5 + 4j / 5],
            [1, 3 / 5 - 4j / 5, 0],
        ],
        dtype=complex,
    )
    result = mod.solve_complex_triangle(
        mod.unit_flux_matrix(sp.Rational(3, 5), sp.Rational(4, 5))
    )
    lower = float(result.minimum.lower)
    upper = float(result.maximum.upper)
    radii = np.sqrt(np.array([3 / 5, 1 / 4, 3 / 20]))
    for alpha, beta in itertools.product(
        np.linspace(0, 2 * math.pi, 41, endpoint=False),
        np.linspace(0, 2 * math.pi, 43, endpoint=False),
    ):
        state = radii * np.exp(1j * np.array([0.0, alpha, beta]))
        energy = float(np.real(np.vdot(state, matrix @ state)))
        assert lower - 1e-12 <= energy <= upper + 1e-12


def test_explicit_pair_channels_have_the_recorded_fermionic_signs():
    apply = mod.physical_operator_realization()["fermionic_channel_signs"]
    assert apply == {
        "D0_to_D1_via_A34_dagger_A12": 1,
        "D0_to_D2_via_A35_dagger_A02": -1,
        "D1_to_D2_via_A15_dagger_A04": 1,
    }


def test_zero_flux_and_invalid_inputs_are_rejected_exactly():
    with pytest.raises(mod.ExactSolverError, match="cycle flux"):
        mod.solve_complex_triangle(mod.unit_flux_matrix(sp.Rational(1), sp.Rational(0)))
    with pytest.raises(ValueError, match="Hermitian"):
        mod.solve_complex_triangle([[0, 1, 1], [0, 0, sp.I], [1, -sp.I, 0]])
    with pytest.raises(ValueError, match="rational Gaussian"):
        mod.solve_complex_triangle(
            [[0, sp.sqrt(2), 1], [sp.sqrt(2), 0, sp.I], [1, -sp.I, 0]]
        )


def test_generated_artifact_is_current_and_has_two_regressions():
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert artifact == mod.build_artifact()
    assert [record["name"] for record in artifact["regressions"]] == [
        "quarter_flux",
        "pythagorean_flux_cos_3_5_sin_4_5",
    ]
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--check"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_standard_library_verifier_passes_without_project_imports():
    completed = subprocess.run(
        [sys.executable, "-I", str(STANDALONE), str(ARTIFACT)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "all 2 complex fixed-1RDM observable-range certificates PASS" in completed.stdout


def test_flux_family_reproduces_the_real_symmetric_regressions():
    """Cross-artifact reconciliation with the real-symmetric benchmark.

    The unit-edge flux family claims to bridge the two real-symmetric
    coherence-triangle cases already shipped in
    results/data/fixed_1rdm_observable_ranges.json. Those were produced by a
    different engine (candidate reduction on the two branches), so agreement
    is a genuine independent check rather than a restatement. Pinned here
    because the two artifacts are otherwise free to drift apart.
    """
    sympy = pytest.importorskip("sympy")
    real_path = ARTIFACT.parent / "fixed_1rdm_observable_ranges.json"
    if not real_path.exists():
        pytest.skip("real-symmetric observable-range artifact not present")

    fam = json.loads(ARTIFACT.read_text())["unit_edge_flux_family"]
    u, z = sympy.symbols("u z")
    coeffs = [sympy.sympify(c) for c in fam["coefficients_descending_in_z"]]
    poly = sum(c * z ** (len(coeffs) - 1 - i) for i, c in enumerate(coeffs))

    # u = cos(Phi): +1 is the unfrustrated triangle, -1 the pi-flux one
    examples = json.loads(real_path.read_text())["examples"]
    for value, name in ((1, "unfrustrated_coherence_triangle"),
                        (-1, "pi_flux_coherence_triangle")):
        roots = sympy.real_roots(sympy.Poly(poly.subs(u, value), z))
        lo, hi = sympy.simplify(roots[0]), sympy.simplify(roots[-1])
        rng = examples[name]["range"]
        assert sympy.simplify(lo - sympy.sympify(rng["minimum"])) == 0, name
        assert sympy.simplify(hi - sympy.sympify(rng["maximum"])) == 0, name

    # the family must also specialise to the shipped quarter-flux sextic
    quarter = next(r for r in json.loads(ARTIFACT.read_text())["regressions"]
                   if r["name"] == "quarter_flux")
    qc = [sympy.Integer(c) for c in
          quarter["energy_critical_polynomial_coefficients_descending"]]
    qpoly = sum(c * z ** (len(qc) - 1 - i) for i, c in enumerate(qc))
    assert sympy.expand(poly.subs(u, 0) - qpoly) == 0
    # and satisfy its stated conjugation symmetry
    assert sympy.simplify(sympy.expand(poly.subs({u: -u, z: -z})) - sympy.expand(poly)) == 0
