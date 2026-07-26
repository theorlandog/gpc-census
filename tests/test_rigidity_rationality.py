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
