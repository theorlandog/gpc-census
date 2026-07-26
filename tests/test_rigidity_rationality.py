"""The rigidity-rationality census and the attainer-gate audit.

Two artifacts, one dependency: the rationality claim is measured on the
shipped representatives, and 154 of them are not natural-orbital
representatives, so the gate audit is what bounds how the claim may be
quoted. These tests pin both, including the scope limits.
"""
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"


@pytest.fixture(scope="module")
def gate():
    return json.loads((DATA / "attainer_gate_audit.json").read_text())


@pytest.fixture(scope="module")
def census():
    return json.loads((DATA / "rigidity_rationality.json").read_text())


# ------------------------------------------------------------------ the gate


def test_gate_failures_are_exactly_the_interference_records_bar_two(gate):
    s = gate["summary"]
    assert s["records"] == 799
    assert s["pass"] == 645
    assert s["fail"] == 154
    assert s["fail_by_class"] == {"INTERFERENCE": 154}
    assert s["pass_by_class"] == {
        "DESIGN-INT": 631,
        "DESIGN-REAL": 12,
        "INTERFERENCE": 2,
    }
    assert s["interference_passing_the_gate"] == ["(3,10) v103", "(3,10) v89"]


def test_every_record_still_attains_the_spectrum(gate):
    """The failures are a basis statement, not a wrong-state statement."""
    s = gate["summary"]
    assert s["all_records_attain_the_spectrum"]
    assert s["worst_spectrum_error"] < 1e-8
    assert s["worst_norm_error"] < 1e-8


def test_no_failure_is_merely_a_permuted_relabelling(gate):
    assert gate["summary"]["failures_that_are_merely_permuted"] == 0


def test_v_B_is_not_a_special_case(gate):
    """The handoff expected one violation; v_B is one of 154."""
    rows = {r["label"]: r for r in gate["records"]}
    vb = rows["(4,9) v65"]
    assert not vb["passes_gate"]
    assert vb["attains_spectrum"]
    assert abs(vb["max_offdiagonal"] - 3 / 23) < 1e-12


# --------------------------------------------------------------- the census


def test_rigid_implies_rational_with_no_counterexamples(census):
    t = census["all_records"]
    assert t["states"] == 799
    assert t["rigid"] == 740
    assert t["deformable"] == 59
    assert t["rigid_and_rational"] == 740
    assert t["rigid_and_irrational"] == 0
    assert t["counterexample_labels"] == []
    assert t["implication_holds"]


def test_the_only_irrational_record_is_deformable(census):
    t = census["all_records"]
    assert t["irrational_labels"] == ["(4,9) v65"]
    rows = {r["label"]: r for r in census["records"]}
    vb = rows["(4,9) v65"]
    assert not vb["rational_weights"]
    assert not vb["fixed_rho_rigid"]
    assert vb["weight_field_degree"] == 2
    assert vb["weight_field_generators"] == [35]


def test_deformability_does_not_force_irrationality(census):
    t = census["all_records"]
    assert t["deformable_and_rational"] == 58
    assert t["deformable_and_irrational"] == 1


def test_the_power_warning_is_recorded(census):
    """All discriminating power sits in the non-natural-orbital records."""
    rows = census["records"]
    gate_passing = [r for r in rows if r["passes_attainer_gate"]]
    assert len(gate_passing) == 645
    assert all(r["fixed_rho_rigid"] for r in gate_passing)
    deformable = [r for r in rows if not r["fixed_rho_rigid"]]
    assert len(deformable) == 59
    assert all(not r["passes_attainer_gate"] for r in deformable)
    assert "power_warning" in census
    assert census["gate_passing_records_only"]["deformable"] == 0


def test_weight_field_is_recomputed_not_read_off_the_encoding(census):
    """Rationality is decided from the values, not from int-versus-string."""
    for r in census["records"]:
        if r["rational_weights"]:
            assert r["weight_field_degree"] == 1
        else:
            assert r["weight_field_degree"] > 1
            assert r["weight_field_generators"]


# ------------------------------------------------- the v_B Galois orbit test


@pytest.fixture(scope="module")
def orbit():
    return json.loads((DATA / "vb_galois_orbit.json").read_text())


def test_vb_galois_orbit_has_size_two(orbit):
    """Both roots of the wall polynomial attain, exactly. Degree = orbit size."""
    assert orbit["orbit_size"] == 2
    assert orbit["verdict"] == "CONFIRMED"
    assert len(orbit["embeddings"]) == 2
    for e in orbit["embeddings"]:
        assert e["attains"]
        assert e["attains_exactly"]
        assert e["all_weights_positive"]
        assert e["attaining_sign_patterns"] == 64


def test_exactly_one_embedding_is_the_shipped_root(orbit):
    shipped = [e for e in orbit["embeddings"] if e["is_shipped_root"]]
    assert len(shipped) == 1
    assert shipped[0]["root"] == "7/17 - 18*sqrt(35)/85"


# ------------------------------------------------ the support-ratio scorecard


@pytest.fixture(scope="module")
def ratio():
    return json.loads((DATA / "no_support_ratio.json").read_text())


def test_ratio_above_one_implies_generic_interference(ratio):
    """The surviving one-directional diagnostic, zero exceptions."""
    above = [r for r in ratio["records"] if r["no_support_ratio"] > 1]
    assert len(above) == 141
    assert all(r["regime"] == "INTERFERENCE-GENERIC" for r in above)


def test_the_converse_fails_with_thirteen_flat_generics(ratio):
    assert ratio["P_NSR_1"]["verdict"] == "FAIL"
    assert ratio["P_NSR_1"]["exceptions"] == 13
    assert ratio["P_NSR_2"]["verdict"] == "PASS"
    flat = [
        r
        for r in ratio["records"]
        if r["regime"] == "INTERFERENCE-GENERIC" and r["no_support_ratio"] == 1
    ]
    assert len(flat) == 13
    # kernel freedom is necessary but not sufficient for staying flat
    assert all(r["incidence_kernel_dim"] >= 1 for r in flat)
    grew_with_kernel = [
        r
        for r in ratio["records"]
        if r["regime"] == "INTERFERENCE-GENERIC"
        and r["no_support_ratio"] > 1
        and r["incidence_kernel_dim"] >= 1
    ]
    assert grew_with_kernel  # so it is not sufficient


def test_both_cancellation_states_are_flat(ratio):
    cancel = [r for r in ratio["records"] if r["regime"] == "CANCELLATION"]
    assert len(cancel) == 2
    assert all(r["no_support_ratio"] == 1 for r in cancel)


def test_the_prereg_disclosure_is_recorded(ratio):
    """Not a blind prereg, and the artifact must say so."""
    assert "NOT a blind pre-registration" in ratio["disclosure"]


# --------------------------------------------------- gap (b), measured


@pytest.fixture(scope="module")
def uniqueness():
    return json.loads((DATA / "uniqueness_census.json").read_text())


def test_uniqueness_prereg_passes_after_the_loop_basis_fix(uniqueness):
    """Originally scored FAIL at v103. The exception was an artifact of a
    non-integer loop basis, which makes holonomy cosines depend on the phase
    branch; with an integer basis every state has one solution modulo gauge."""
    assert uniqueness["verdict"] == "PASS"
    assert uniqueness["exceptions"] == 0
    assert uniqueness["exception_labels"] == []
    assert uniqueness["sample_size"] == 12
    assert "ARTIFACT" in uniqueness["finding"]
    for r in uniqueness["records"]:
        assert r["distinct_modulo_gauge"] == 1
        assert r["distinct_weight_multisets"] == 1


def test_weight_multiset_is_unique_everywhere(uniqueness):
    """The failure is phase multiplicity, not weight multiplicity, so the
    Rigidity-Rationality target is untouched by it."""
    assert uniqueness["weights_verdict"] == "PASS"
    assert uniqueness["weight_exceptions"] == 0
    for r in uniqueness["records"]:
        assert r["distinct_weight_multisets"] == 1
        assert r["converged"] == r["starts"]


def test_the_solved_system_is_the_fixed_rho_fiber(uniqueness):
    """Targeting diag(lambda) would be infeasible for the 154 non-diagonal
    records, so the system must target the state's own rho."""
    assert "OWN 1-RDM" in uniqueness["system_solved"]
    assert uniqueness["evidence"].startswith("NUMERICAL")


# ------------------------------------------------- v103, the real number


@pytest.fixture(scope="module")
def fiber_count():
    return json.loads((DATA / "v103_fiber_count.json").read_text())


def test_v103_fiber_has_exactly_one_point(fiber_count):
    """Exhaustive search of the 5-torus of gauge-invariant holonomies, which is
    the right space once the weight vector is known to be common to every
    solution."""
    assert fiber_count["distinct_points"] == 1
    assert fiber_count["converged"] == fiber_count["starts"]
    assert fiber_count["starts"] >= 16807
    assert fiber_count["worst_reconstruction_error"] < 1e-9


def test_v103_single_holonomy_is_all_signs(fiber_count):
    """The one point has holonomy (0, pi, pi, 0, 0): the state is real up to
    gauge, as the library record says."""
    import math

    (holonomy,) = fiber_count["holonomies"]
    for value in holonomy:
        assert min(abs(value), abs(value - math.pi), abs(value - 2 * math.pi)) < 1e-6
