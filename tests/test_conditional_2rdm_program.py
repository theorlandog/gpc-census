"""Guards for steps 4 to 7 of the conditional 2-RDM program.

The load-bearing checks are the ones that would catch a silently wrong answer:

- the Slater-Condon carrier projection must reproduce the full CI submatrix
  exactly, since every downstream interval is computed from it;
- the fermionic channel signs must match the independently shipped artifact;
- the exact pair-current range must match the closed form, which is also the
  regression that caught the dual solver diverging on rank-deficient carriers;
- the preconditions (nondegenerate occupations, carrier is the dominant
  support) must be enforced, not assumed, because the quasipinning statement is
  false without them.
"""
import json
import pathlib
import subprocess
import sys

import numpy as np
import pytest

from gpc_census import elliptope as el
from gpc_census import pinned_physics as pp

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "conditional_2rdm_program.json"
SEXTIC = ROOT / "complex_fixed_1rdm_observable_ranges.json"

W = [0.6, 0.25, 0.15]


def _prog():
    return json.loads(ARTIFACT.read_text())


# --------------------------------------------------------------------------
# second-quantized machinery
# --------------------------------------------------------------------------

def test_operator_algebra_signs():
    # a_1^dag a_0 |0,2> = -|1,2>: one occupied orbital sits left of index 1
    got = pp.apply_ops((0, 2), [("c", 1), ("a", 0)])
    assert got == ((1, 2), 1)
    assert pp.apply_ops((0, 2), [("a", 1)]) is None       # not occupied
    assert pp.apply_ops((0, 2), [("c", 2)]) is None       # already occupied


def test_carrier_projection_matches_full_ci_exactly():
    h1, g2 = pp.tv_ring_integrals(6, 1.0, 2.0, gradient=0.35)
    H, dets = pp.ci_hamiltonian(h1, g2, 6, 3)
    idx = {D: i for i, D in enumerate(dets)}
    Hc = pp.carrier_projection(h1, g2)
    sub = np.array([[H[idx[A], idx[B]] for B in pp.DETS] for A in pp.DETS])
    assert np.abs(Hc - sub).max() == 0.0
    assert np.abs(Hc - Hc.conj().T).max() < 1e-14


def test_fermionic_channel_signs_match_the_shipped_artifact():
    if not SEXTIC.exists():
        pytest.skip("sextic artifact not present")
    want = json.loads(SEXTIC.read_text())["physical_operator_realization"][
        "fermionic_channel_signs"]
    pairs = {((0, 1), "D0_to_D1_via_A34_dagger_A12"),
             ((0, 2), "D0_to_D2_via_A35_dagger_A02"),
             ((1, 2), "D1_to_D2_via_A15_dagger_A04")}
    for (a, b), key in pairs:
        _, _, sign = pp.excitation_sign(pp.DETS[a], pp.DETS[b])
        assert sign == want[key], key


def test_fcidump_round_trip():
    import tempfile
    h1, g2 = pp.tv_ring_integrals(6, 1.0, 2.0, twist=0.4, gradient=0.35)
    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / "T.FCIDUMP"
        pp.write_fcidump(path, h1, g2, n_elec=3)
        h1b, g2b, norb = pp.read_fcidump(path)
    assert norb == 6
    assert np.allclose(h1, h1b)
    assert np.allclose(g2, g2b)


def test_one_rdm_and_natural_orbitals_are_consistent():
    h1, g2 = pp.tv_ring_integrals(6, 1.0, 2.0, gradient=0.35)
    H, dets = pp.ci_hamiltonian(h1, g2, 6, 3)
    psi = np.linalg.eigh(H)[1][:, 0]
    gamma = pp.one_rdm_ci(psi, dets, 6)
    assert np.abs(gamma - gamma.conj().T).max() < 1e-12
    lams, U = pp.natural_orbitals(gamma)
    assert abs(lams.sum() - 3.0) < 1e-10          # trace is the particle number
    assert all(-1e-12 <= x <= 1 + 1e-12 for x in lams)
    # the Borland-Dennis equalities hold identically for any three-fermion state
    for i in range(3):
        assert abs(lams[i] + lams[5 - i] - 1.0) < 1e-9
    # rotating to the natural-orbital basis preserves the norm and diagonalizes
    # the 1-RDM under this module's gamma[p,q] convention
    psi_no = pp.natural_orbital_coefficients(psi, dets, U)
    assert abs(np.linalg.norm(psi_no) - 1.0) < 1e-9
    gamma_no = pp.one_rdm_ci(psi_no, dets, 6)
    assert np.abs(gamma_no - np.diag(np.diag(gamma_no))).max() < 1e-12


def test_complex_natural_basis_uses_transpose_not_adjoint():
    h1, g2 = pp.tv_ring_integrals(
        6, 1.0, 2.0, gradient=7 / 20, boundary_phase=(4 + 3j) / 5,
    )
    H, dets = pp.ci_hamiltonian(h1, g2, 6, 3)
    psi = np.linalg.eigh(H)[1][:, 0]
    _, U = pp.natural_orbitals(pp.one_rdm_ci(psi, dets, 6))

    correct = pp.natural_orbital_coefficients(psi, dets, U)
    correct_gamma = pp.one_rdm_ci(correct, dets, 6)
    correct_offdiag = np.abs(
        correct_gamma - np.diag(np.diag(correct_gamma))
    ).max()

    old_wrong_convention = pp.rotate_ci(psi, dets, U.conj().T)
    wrong_gamma = pp.one_rdm_ci(old_wrong_convention, dets, 6)
    wrong_offdiag = np.abs(wrong_gamma - np.diag(np.diag(wrong_gamma))).max()
    assert correct_offdiag < 1e-12
    assert wrong_offdiag > 0.1


def test_complex_natural_basis_transforms_state_and_hamiltonian_together():
    h1, g2 = pp.tv_ring_integrals(
        6, 1.0, 2.0, gradient=7 / 20, boundary_phase=(4 + 3j) / 5,
    )
    H, dets = pp.ci_hamiltonian(h1, g2, 6, 3)
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    psi = eigenvectors[:, 0]
    _, U = pp.natural_orbitals(pp.one_rdm_ci(psi, dets, 6))
    psi_no = pp.natural_orbital_coefficients(psi, dets, U)
    h1_no, g2_no = pp.natural_basis_integrals(h1, g2, U)
    H_no, dets_no = pp.ci_hamiltonian(h1_no, g2_no, 6, 3)
    assert dets_no == dets
    assert np.max(np.abs(H_no @ psi_no - eigenvalues[0] * psi_no)) < 1e-12


def test_direct_boundary_phase_must_have_unit_modulus():
    with pytest.raises(ValueError, match="unit modulus"):
        pp.tv_ring_integrals(6, 1.0, 2.0, boundary_phase=1 + 1j)


# --------------------------------------------------------------------------
# block 4: flux and current
# --------------------------------------------------------------------------

def test_pair_current_range_matches_the_closed_form():
    """Regression: the dual solver used to diverge on this operator.

    J touches a single edge, so the optimal dual face is degenerate and an
    unrepaired SLSQP solve returned values of order 1e7 at some fluxes.
    """
    exact = 2 * np.sqrt(W[1] * W[2])
    for k in range(12):
        phi = 2 * np.pi * k / 12
        J = pp.current_carrier(phi)
        assert el.support_numeric(J, W, "max") == pytest.approx(exact, abs=1e-7)
        assert el.support_numeric(J, W, "min") == pytest.approx(-exact, abs=1e-7)


def test_current_operator_is_hermitian_and_is_minus_dH_dphi():
    for phi in (0.0, 0.7, 2.5):
        J = pp.current_carrier(phi)
        assert np.abs(J - J.conj().T).max() < 1e-14
        d = 1e-6
        num = -(pp.flux_carrier(phi + d) - pp.flux_carrier(phi - d)) / (2 * d)
        assert np.abs(J - num).max() < 1e-6


def test_hidden_current_is_nonzero_at_every_flux():
    """Same 1-RDM, different pair current: the physics claim of the block."""
    prog = _prog()["flux_and_current"]
    for row in prog["flux_sweep"]:
        lo, hi = row["current_range"]
        assert hi - lo > 0.7, row["phi"]
    lo, hi = prog["exact_current_range_approx"]
    assert hi == pytest.approx(2 * np.sqrt(W[1] * W[2]), abs=1e-12)
    assert lo == pytest.approx(-hi, abs=1e-12)


def test_branch_exchange_fluxes_include_the_degenerate_ends():
    prog = _prog()["flux_and_current"]
    us = sorted(round(b["u_approx"], 9) for b in prog["branch_exchange_fluxes"])
    assert -1.0 in us and 1.0 in us          # the real-symmetric limits
    assert len(us) >= 4                      # plus interior exchanges


# --------------------------------------------------------------------------
# blocks 5 and 6: projection and quasipinning
# --------------------------------------------------------------------------

def test_projection_block_is_exact_and_preconditions_are_enforced():
    b = _prog()["projection_and_quasipinning"]
    assert b["fcidump_round_trip"]
    assert b["carrier_projection_matches_full_ci_max_abs_error"] == 0.0
    assert b["preconditions"]["usable_points"] >= 4
    assert b["preconditions"]["usable_genuinely_complex_flux_points"] >= 4
    assert b["preconditions"]["usable_points"] < b["preconditions"]["total_points"]
    # points failing the preconditions must be reported, not silently dropped
    failing = [r for r in b["scan"] if not r["carrier_is_dominant_support"]]
    assert failing


def test_enlargement_bound_holds_where_its_preconditions_hold():
    b = _prog()["projection_and_quasipinning"]
    usable = [r for r in b["scan"] if r["pinned_preconditions_pass"]]
    assert usable
    for r in usable:
        assert r["true_energy_inside_enlarged"], r["V"]
        # and pinning genuinely narrows the prediction
        assert r["width_ratio_unconstrained_over_pinned"] > 1.0


def test_complex_flux_physical_benchmark_passes_every_frozen_gate():
    b = _prog()["projection_and_quasipinning"]["complex_flux_physical_benchmark"]
    assert b["status"] == "NUMERICAL_WITH_EXACT_MODEL_INPUT"
    assert all(b["gates"].values())
    assert b["exact_model_input"]["boundary_phase"] == "(4+3*i)/5"
    row = b["selected_scan_point"]
    assert row["ground_state_gap"] > 1.0
    assert row["min_occupation_gap"] > 1e-3
    assert row["natural_basis_offdiagonal_max"] < 1e-12
    assert row["off_carrier_weight"] < 0.01
    assert abs(row["carrier_cycle_product"][1]) > 1e-3
    assert row["true_to_carrier_energy_deviation"] <= row["enlargement_half_width"]


def test_complex_flux_standalone_verifier_passes():
    script = pathlib.Path(__file__).resolve().parents[1] / "scripts" / (
        "verify_complex_flux_physical_benchmark_standalone.py"
    )
    subprocess.run(
        [sys.executable, str(script), str(ARTIFACT)],
        check=True,
        capture_output=True,
        text=True,
    )


def test_enlargement_bound_is_monotone_in_its_inputs():
    base = pp.enlargement_bound(1.0, 0.1, 0.0, 1.0)
    assert pp.enlargement_bound(2.0, 0.1, 0.0, 1.0) > base      # larger ||H||
    assert pp.enlargement_bound(1.0, 0.2, 0.0, 1.0) > base      # larger eps
    assert pp.enlargement_bound(1.0, 0.1, 0.5, 1.0) > base      # weights moved
    assert pp.enlargement_bound(1.0, 0.0, 0.0, 1.0) == 0.0      # exactly pinned


def test_slack_controls_the_off_carrier_weight_on_usable_points():
    b = _prog()["projection_and_quasipinning"]
    lo, hi = b["slack_versus_off_carrier_weight"]["observed_ratio_range"]
    assert 0.5 < lo <= hi < 4.0


# --------------------------------------------------------------------------
# block 7: the atlas
# --------------------------------------------------------------------------

def test_atlas_covers_the_whole_census_and_tiers_it():
    a = _prog()["atlas"]
    assert a["n_records"] == 799
    assert a["tiers"]["exact_elliptope_sdp"] + a["tiers"]["elliptope_outer_bound"] == 799
    for rec in a["records"]:
        assert rec["exactly_spectrahedral"] == (rec["carrier_size"] <= 3)
        assert rec["method"] == ("exact_elliptope_sdp" if rec["carrier_size"] <= 3
                                 else "elliptope_outer_bound")
        assert (rec["single_excitation_pairs"] + rec["double_excitation_pairs"]
                + rec["inaccessible_pairs"]) == rec["coherence_pairs"]
        assert rec["coherence_graph_complete"] == (rec["inaccessible_pairs"] == 0)


def test_atlas_records_two_body_invisible_directions():
    """Most census carriers have coherences no two-body observable can see."""
    a = _prog()["atlas"]
    incomplete = [r for r in a["records"] if not r["coherence_graph_complete"]]
    assert incomplete, "expected supports with inaccessible coherence pairs"
    for r in incomplete:
        assert r["zero_width_two_body_dimension"] > r["carrier_size"]
