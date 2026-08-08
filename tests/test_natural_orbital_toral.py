"""Guards for the natural-orbital extension of the toral obstruction."""
from __future__ import annotations

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" if (ROOT / "tests" / "data").exists() else ROOT / "results" / "data"
ARTIFACT = DATA / "natural_orbital_toral.json"

pytestmark = pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact absent")


@pytest.fixture(scope="module")
def art():
    return json.loads(ARTIFACT.read_text())


def test_targets_are_exactly_the_out_of_scope_states(art):
    prev = json.loads((DATA / "toral_component_obstruction.json").read_text())
    out = {f"{r['system']}#{r['index']}" for r in prev["rows"] if not r["in_scope"]}
    assert len(out) == 154
    covered = {f"{r['system']}#{r['index']}" for r in art["rows"]}
    assert covered | set(art["unfinished"]) == out
    assert not (covered & set(art["unfinished"]))


def test_every_target_is_interference(art):
    assert art["summary"]["all_interference"]
    assert all(r["classified"] == "INTERFERENCE" for r in art["rows"])


def test_all_attempted_rows_resolved(art):
    assert art["summary"]["unresolved"] == 0
    assert all(r["resolved"] for r in art["rows"])


def test_obstruction_never_fires_on_the_interference_class(art):
    assert art["summary"]["ratios_seen"] == [1]
    assert art["summary"]["obstruction_fires"] == 0
    for r in art["rows"]:
        assert r["toral_period"] == r["denominator"], r


def test_no_prediction_is_raised(art):
    assert art["summary"]["predictions_raised"] == 0
    for r in art["rows"]:
        assert r["combined_predicted_period"] == r["permutation_predicted_period"], r


def test_two_eigenbasis_variants_were_combined(art):
    for r in art["rows"]:
        d = r["detail"]
        assert "plain" in d and "mixed" in d, r
        assert d["combined"] == r["toral_period"], r
