"""Regression tests for the fixed-spectrum sparsification cascade.

The cascade lowers the UPPER BOUND on minimal support s_Q. These tests pin the
two behaviours that make its output trustworthy: a landed drop reproduces the
target spectrum with the right MULTIPLICITIES (not merely the right eigenvalue
set, which is all the minimal-polynomial certificate constrains), and a state
that cannot be sparsified on this path reports a floor rather than a drop.
"""
import json
import pathlib
import sys
from fractions import Fraction

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data" / "states.jsonl"


def _record(system, index):
    if not DATA.exists():
        pytest.skip("no states atlas present")
    for line in DATA.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            if r["system"] == system and r["index"] == index:
                return r
    raise AssertionError((system, index))


def _case(system, index, d):
    from gpc_census.fiber import amplitudes_from_record

    rec = _record(system, index)
    c, dets = amplitudes_from_record(rec)
    den = int(rec["denominator"])
    spectrum = [Fraction(int(x), den) for x in rec["integer_form"]]
    return c, dets, spectrum, d


def test_distinct_eigenvalues_dedups():
    """One factor per DISTINCT eigenvalue: the certificate stalls otherwise."""
    from gpc_census.cascade import distinct_eigenvalues

    spec = [Fraction(18, 34)] * 4 + [Fraction(5, 34)] * 6
    assert distinct_eigenvalues(spec) == [18 / 34, 5 / 34]


def test_v103_sparsifies_on_the_spectrum_locus_only():
    """v103 reaches support 10 on the SPECTRUM locus, and its endpoint is not diagonal.

    This bounds s_Q^free, not s_Q^NO. The endpoint's 1-RDM has off-diagonal
    5/68, so it does not attain the vertex in the census diagonal convention.
    """
    from gpc_census.cascade import cascade

    out = cascade(*_case("(3,10)", 103, 10), locus="fixed_spectrum")
    assert out["support_start"] == 14
    assert out["support_upper"] <= 12, out["support_upper"]
    assert out["final_residual"] < 1e-10
    assert out["final_spectrum_error"] < 1e-8
    assert out["locus"] == "fixed_spectrum"
    assert out["bounds"] == "s_Q_free"
    # the whole point of the correction: the endpoint left the diagonal locus
    assert out["final_diagonality_error"] > 1e-3
    for rnd in out["rounds"]:
        assert rnd["certificate_residual"] < 1e-10
        assert rnd["spectrum_error"] < 1e-8


def test_v103_does_not_sparsify_on_the_natural_orbital_locus():
    """Holding rho = diag(lambda), v103 does not sparsify at all.

    Guards the corrected claim: s_Q^NO(v103) <= 10 is NOT established.
    """
    from gpc_census.cascade import cascade

    out = cascade(*_case("(3,10)", 103, 10), locus="natural_orbital")
    assert out["support_upper"] == out["support_start"] == 14
    assert out["bounds"] == "s_Q_NO"
    assert out["final_diagonality_error"] < 1e-9


def test_off_locus_state_refuses_the_natural_orbital_cascade():
    """A library state with a non-diagonal 1-RDM cannot start on that locus."""
    from gpc_census.cascade import cascade

    out = cascade(*_case("(4,9)", 65, 9), locus="natural_orbital")
    assert out["blocked_reason"] == "start_off_locus"
    assert out["support_upper"] == out["support_start"]


def test_v89_does_not_sparsify_on_this_path():
    """v89 floors; the result must be a recorded floor, never a claimed drop."""
    from gpc_census.cascade import cascade

    out = cascade(*_case("(3,10)", 89, 10), locus="fixed_spectrum")
    assert out["support_upper"] == out["support_start"] == 10
    assert out["blocked_on_path"] is True
    assert out["blocked_reason"] in ("certificate_floor", "multiplicity_drift")
    assert out["floor_residual"] is not None


def test_spectrum_error_catches_multiplicity_drift():
    """The minimal polynomial pins the eigenvalue SET; multiplicities need this.

    Build a state on v103's spectrum values with the WRONG multiplicities and
    check the gate rejects it while the minimal-polynomial certificate does not.
    """
    from gpc_census.cascade import (certificate_residual, distinct_eigenvalues,
                                    spectrum_error)

    # a single Slater determinant has eigenvalues 1, 1, 1 and zeros: the
    # eigenvalue set of the spectrum (1,1,1,0,...) with multiplicities that
    # differ from the (2,1,1,1,...) target below
    dets = [(0, 1, 2)]
    c = np.array([1.0 + 0j])
    target = [Fraction(1)] * 3 + [Fraction(0)] * 7
    lams = distinct_eigenvalues(target)
    assert certificate_residual(c, dets, lams, 10) < 1e-12
    assert spectrum_error(c, dets, target, 10) < 1e-12
    wrong = [Fraction(1)] * 2 + [Fraction(0)] * 8  # same value set, wrong count
    assert spectrum_error(c, dets, wrong, 10) > 0.5


def test_cascade_reports_its_locus_and_evidence_level():
    """Every result names WHICH minimum it bounds and its evidence level."""
    from gpc_census.cascade import cascade

    for locus, bounds in (("natural_orbital", "s_Q_NO"),
                          ("fixed_spectrum", "s_Q_free")):
        out = cascade(*_case("(3,10)", 89, 10), locus=locus)
        assert out["locus"] == locus
        assert out["bounds"] == bounds
        assert out["evidence"] == "numerical"
        assert "method" in out and out["method"]
