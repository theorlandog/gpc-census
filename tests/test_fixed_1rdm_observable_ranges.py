import importlib.util
import json
import math
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "fixed_1rdm_observable_ranges.py"
VERIFIER = ROOT / "scripts" / "verify_fixed_1rdm_observable_ranges_standalone.py"
ARTIFACT = ROOT / "results" / "data" / "fixed_1rdm_observable_ranges.json"

spec = importlib.util.spec_from_file_location("fixed_ranges", SCRIPT)
fixed_ranges = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(fixed_ranges)


def test_generated_artifact_is_current():
    expected = json.dumps(fixed_ranges.build(), indent=2, sort_keys=True) + "\n"
    assert ARTIFACT.read_text() == expected


def test_borland_dennis_selection_rule_reduces_global_fiber_to_three_determinants():
    certificate = fixed_ranges.borland_dennis_selection_certificate()
    assert certificate["ambient_determinants"] == 20
    assert certificate["single_occupancy_support_size"] == 8
    assert certificate["facet_kernel_support_size"] == 9
    assert certificate["intersection_support_size"] == 3
    assert certificate["equals_benchmark_support"] is True


def test_target_1rdm_is_nondegenerate_pinned_and_fixes_the_weights():
    gamma = fixed_ranges.diagonal_1rdm(
        fixed_ranges.BD_SUPPORT, fixed_ranges.BD_WEIGHTS, 6
    )
    assert gamma == (
        sp.Rational(17, 20),
        sp.Rational(3, 4),
        sp.Rational(3, 5),
        sp.Rational(2, 5),
        sp.Rational(1, 4),
        sp.Rational(3, 20),
    )
    assert all(gamma[index] > gamma[index + 1] for index in range(5))
    assert gamma[0] + gamma[5] == 1
    assert gamma[1] + gamma[4] == 1
    assert gamma[2] + gamma[3] == 1
    assert 2 - (gamma[0] + gamma[1] + gamma[3]) == 0
    recovered, certificate = fixed_ranges.recover_weights(
        fixed_ranges.BD_SUPPORT, gamma
    )
    assert recovered == fixed_ranges.BD_WEIGHTS
    assert certificate == {
        "rank": 3,
        "pivot_orbitals": [0, 1, 2],
        "pivot_minor_determinant": 1,
    }
    assert fixed_ranges.one_hop_pairs(fixed_ranges.BD_SUPPORT) == []


def test_each_triangle_edge_has_a_collision_free_two_body_channel():
    witnesses = fixed_ranges.complete_collision_free_witnesses(
        fixed_ranges.BD_SUPPORT, order=2
    )
    assert witnesses == [
        {
            "edge": [0, 1],
            "channel_left_modes": [1, 2],
            "channel_right_modes": [3, 4],
            "channel_coefficient": 1,
            "operator_multiplier": "1",
        },
        {
            "edge": [0, 2],
            "channel_left_modes": [0, 2],
            "channel_right_modes": [3, 5],
            "channel_coefficient": -1,
            "operator_multiplier": "-1",
        },
        {
            "edge": [1, 2],
            "channel_left_modes": [0, 4],
            "channel_right_modes": [1, 5],
            "channel_coefficient": 1,
            "operator_multiplier": "1",
        },
    ]


def test_unfrustrated_exact_range_matches_independent_polygon_formula():
    general = fixed_ranges.real_symmetric_triangle_range(
        fixed_ranges.COHERENCE_TRIANGLE_MATRIX
    )
    polygon = fixed_ranges.complete_triangle_range(fixed_ranges.BD_WEIGHTS)
    assert general["minimum"] == polygon["minimum"] == "-1"
    assert general["maximum"] == polygon["maximum"] == "3*(2 + sqrt(15))/10"
    assert general["width"] == polygon["width"] == "(3*sqrt(15) + 16)/10"
    assert general["minimum_source"] == "interior stationary point"
    assert general["maximum_source"] == "alpha=0 endpoint"
    assert general["minimum_phase_cosines"] == {
        "cos_alpha": "-7*sqrt(15)/30",
        "cos_beta": "-5/6",
        "cos_beta_minus_alpha": "2*sqrt(15)/15",
    }
    assert general["candidate_reduction"]["candidate_count"] == 5


def test_pi_flux_reverses_which_side_has_the_interior_extremum():
    certificate = fixed_ranges.real_symmetric_triangle_range(
        fixed_ranges.PI_FLUX_TRIANGLE_MATRIX
    )
    assert certificate["minimum"] == "-3*(2 + sqrt(15))/10"
    assert certificate["maximum"] == "1"
    assert certificate["minimum_source"] == "alpha=pi endpoint"
    assert certificate["maximum_source"] == "interior stationary point"
    assert certificate["maximum_phase_cosines"] == {
        "cos_alpha": "7*sqrt(15)/30",
        "cos_beta": "5/6",
        "cos_beta_minus_alpha": "2*sqrt(15)/15",
    }
    assert certificate["candidate_reduction"]["candidate_count"] == 5


def test_diagonal_carrier_observable_is_constant_on_the_fiber():
    matrix = (
        (sp.Rational(1), sp.Rational(0), sp.Rational(0)),
        (sp.Rational(0), sp.Rational(2), sp.Rational(0)),
        (sp.Rational(0), sp.Rational(0), sp.Rational(3)),
    )
    certificate = fixed_ranges.real_symmetric_triangle_range(matrix)
    assert certificate["offset"] == "31/20"
    assert certificate["minimum"] == "31/20"
    assert certificate["maximum"] == "31/20"
    assert certificate["width"] == "0"
    assert certificate["candidate_reduction"]["candidate_count"] == 4


def test_every_rational_real_symmetric_carrier_matrix_has_a_two_body_realization():
    matrix = (
        (sp.Rational(1, 7), sp.Rational(2), sp.Rational(-1)),
        (sp.Rational(2), sp.Rational(-2, 5), sp.Rational(3, 2)),
        (sp.Rational(-1), sp.Rational(3, 2), sp.Rational(4, 9)),
    )
    decomposition = fixed_ranges.carrier_operator_decomposition(matrix)
    assert decomposition["reconstructed_carrier_matrix"] == [
        ["1/7", "2", "-1"],
        ["2", "-2/5", "3/2"],
        ["-1", "3/2", "4/9"],
    ]
    assert decomposition["one_body_terms"] == [
        {"orbital": 0, "coefficient": "-2/5", "operator": "n_0"},
        {"orbital": 1, "coefficient": "4/9", "operator": "n_1"},
        {"orbital": 2, "coefficient": "31/315", "operator": "n_2"},
    ]
    assert len(decomposition["two_body_terms"]) == 3


def _phase_grid_bounds(matrix, samples=160):
    weights = [float(weight) for weight in fixed_ranges.BD_WEIGHTS]
    radii = [math.sqrt(weight) for weight in weights]
    values = []
    for alpha_index in range(samples):
        alpha = 2 * math.pi * alpha_index / samples
        for beta_index in range(samples):
            beta = 2 * math.pi * beta_index / samples
            phases = (0.0, alpha, beta)
            value = sum(weights[index] * float(matrix[index][index]) for index in range(3))
            for left, right in ((0, 1), (0, 2), (1, 2)):
                value += (
                    2
                    * radii[left]
                    * radii[right]
                    * float(matrix[left][right])
                    * math.cos(phases[right] - phases[left])
                )
            values.append(value)
    return min(values), max(values)


def test_exact_solver_contains_dense_phase_grids_for_varied_rational_matrices():
    matrices = [
        fixed_ranges.COHERENCE_TRIANGLE_MATRIX,
        fixed_ranges.PI_FLUX_TRIANGLE_MATRIX,
        (
            (sp.Rational(1, 7), sp.Rational(2), sp.Rational(-1)),
            (sp.Rational(2), sp.Rational(-2, 5), sp.Rational(3, 2)),
            (sp.Rational(-1), sp.Rational(3, 2), sp.Rational(4, 9)),
        ),
        (
            (sp.Rational(-3, 8), sp.Rational(0), sp.Rational(5, 3)),
            (sp.Rational(0), sp.Rational(7, 11), sp.Rational(-2, 9)),
            (sp.Rational(5, 3), sp.Rational(-2, 9), sp.Rational(1, 6)),
        ),
    ]
    for matrix in matrices:
        certificate = fixed_ranges.real_symmetric_triangle_range(matrix)
        sampled_minimum, sampled_maximum = _phase_grid_bounds(matrix)
        assert certificate["minimum_numeric"] <= sampled_minimum + 1e-12
        assert sampled_maximum <= certificate["maximum_numeric"] + 1e-12


def test_qsqrt15_comparisons_are_exact_near_cancellation():
    assert fixed_ranges.q15_sign(4 - sp.sqrt(15)) == 1
    assert fixed_ranges.q15_sign(sp.sqrt(15) - 4) == -1
    assert fixed_ranges.q15_sign(sp.Rational(15, 4) - sp.sqrt(15)) == -1
    assert fixed_ranges.q15_sign(sp.sqrt(15) - sp.Rational(15, 4)) == 1
    assert fixed_ranges.q15_sign(sp.Integer(0)) == 0


def test_standard_library_verifier():
    subprocess.run(
        [sys.executable, str(VERIFIER), str(ARTIFACT)],
        check=True,
        capture_output=True,
        text=True,
    )


def test_artifact_matrix_entries_are_parseable_as_rationals():
    payload = json.loads(ARTIFACT.read_text())
    for example in payload["examples"].values():
        matrix = example["range"]["carrier_matrix"]
        assert all(Fraction(entry) == Fraction(entry) for row in matrix for entry in row)
