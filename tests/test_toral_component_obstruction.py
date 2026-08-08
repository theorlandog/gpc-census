"""Guards for the toral component-group obstruction."""
from __future__ import annotations

import json
import pathlib

import pytest

from gpc_census.symplectic import toral_component_period

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" if (ROOT / "tests" / "data").exists() else ROOT / "results" / "data"
ARTIFACT = DATA / "toral_component_obstruction.json"

pytestmark = pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact absent")


@pytest.fixture(scope="module")
def art():
    return json.loads(ARTIFACT.read_text())


def test_covers_the_census(art):
    assert len(art["rows"]) == 799
    assert len({(r["system"], r["index"]) for r in art["rows"]}) == 799


def test_scope_is_exactly_the_diagonal_states(art):
    """Out of scope must mean non-diagonal 1-RDM, with no exceptions."""
    assert art["summary"]["out_of_scope_all_non_diagonal"]
    for r in art["rows"]:
        if not r["in_scope"]:
            assert not r["rho_is_diagonal"], r
            assert r["toral_component_period"] is None, r


def test_obstruction_is_a_two_group(art):
    assert art["summary"]["ratios_seen"] == [1, 2]
    for r in art["rows"]:
        if r["in_scope"]:
            assert r["toral_component_period"] % r["denominator"] == 0, r
            assert r["toral_component_period"] // r["denominator"] in (1, 2), r


def test_it_fires_exactly_on_the_design_real_class(art):
    fires = [r for r in art["rows"]
             if r["in_scope"] and r["toral_component_period"] != r["denominator"]]
    assert len(fires) == 12
    assert all(r["classified"] == "DESIGN-REAL" for r in fires)
    assert art["summary"]["fires_all_design_real"]


def test_seven_predictions_are_raised(art):
    raised = [r for r in art["rows"] if r["raises_prediction"]]
    assert len(raised) == 7
    for r in raised:
        assert r["combined_predicted_period"] % r["permutation_predicted_period"] == 0, r
        assert r["combined_predicted_period"] > r["permutation_predicted_period"], r


def test_combined_period_is_the_lcm(art):
    from math import gcd
    for r in art["rows"]:
        a = r["permutation_predicted_period"]
        b = r["toral_component_period"]
        want = a if b is None else a * b // gcd(a, b)
        assert r["combined_predicted_period"] == want, r


def test_recompute_two_states(art):
    states = {(x["system"], x["index"]): x
              for x in (json.loads(line)
                        for line in (DATA / "states.jsonl").read_text().splitlines())}
    for key in [("(3,6)", 3), ("(4,9)", 35)]:
        row = next(r for r in art["rows"]
                   if (r["system"], r["index"]) == key)
        assert toral_component_period(states[key]) == row["toral_component_period"]
