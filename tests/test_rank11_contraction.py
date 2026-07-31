"""Exact geometry and numerical gates for the corrected rank-11 audit."""
from __future__ import annotations

import json
import pathlib
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from rank11_contraction_audit import candidate44_geometry  # noqa: E402


def test_candidate44_trace_hyperplane_projection_is_exact():
    geometry = candidate44_geometry()
    q = geometry["projection"]
    assert sum(q) == 3
    assert q == [
        Fraction(37, 72),
        Fraction(29, 60),
        Fraction(29, 60),
        Fraction(29, 60),
        Fraction(29, 60),
        Fraction(7, 72),
        Fraction(7, 72),
        Fraction(7, 72),
        Fraction(7, 72),
        Fraction(7, 72),
        Fraction(1, 15),
    ]
    assert sum(q[i - 1] for i in geometry["indices_one_based"]) == 2
    assert geometry["target_violation"] == Fraction(1, 12)
    assert geometry["trace_hyperplane_normal_squared"] == Fraction(30, 11)
    assert geometry["distance_squared"] == Fraction(11, 4320)


def test_shipped_rank11_audit_has_normalized_controls():
    report = json.loads(
        (ROOT / "results/data/rank11_contraction.json").read_text()
    )
    control = report["attainable_control"]
    assert control["index"] == 49
    for run in control["runs"].values():
        assert run["normalized_residual"] < 1e-18
        assert run["state_norm_error"] < 1e-12
        assert run["one_rdm_trace_error"] < 1e-12

    candidate44 = next(
        row for row in report["open_candidates"] if row["index"] == 44
    )
    assert report["candidate44_claimed_facet_geometry"][
        "projection_design_classification"
    ]["verdict"] == "INTERFERENCE"
    exact_distance = float(Fraction(11, 4320))
    for run in candidate44["runs"].values():
        assert abs(run["normalized_residual"] - exact_distance) < 1e-12
        assert run["projection_spectrum_max_error"] < 1e-7
