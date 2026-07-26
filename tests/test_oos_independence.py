"""Provenance independence of the (3,11)/(3,12) bracket scorecard.

The holdout was scored at 38 rows, but 34 of them are transports or face
embeddings of rank-10 states. These tests pin the corrected accounting so a
regeneration cannot quietly restore the overstated scope, and check that every
verdict carries its independent scope beside its row scope (standing rule R2,
docs/prereg_template.md).
"""
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"

ARTIFACT = DATA / "oos_bracket.json"


@pytest.fixture(scope="module")
def artifact():
    return json.loads(ARTIFACT.read_text())


def test_independence_block_present_and_correct(artifact):
    i = artifact["provenance_independence"]
    assert i["total_rows"] == 38
    assert i["derived_rows"] == 34
    assert i["independent_rows"] == 4
    assert i["distinct_rank10_donors"] == 17
    assert i["distinct_independent_states"] == 2
    assert i["certificate_counts"] == {
        "EXPLICIT-STATE": 2,
        "FACE-EMBEDDING": 19,
        "STATE-TRANSPORT": 17,
    }


def test_every_measurement_traces_to_a_root_origin(artifact):
    per = artifact["provenance_independence"]["per_measurement"]
    assert len(per) == 38
    assert {p["label"] for p in per} == {m["label"] for m in artifact["measurements"]}
    for p in per:
        assert p["root_kind"] in ("census", "new")
        assert p["root_origin"]
        assert p["independent_of_rank10_corpus"] == (p["root_kind"] == "new")


def test_only_the_two_explicit_states_are_independent(artifact):
    per = artifact["provenance_independence"]["per_measurement"]
    indep = sorted(p["label"] for p in per if p["independent_of_rank10_corpus"])
    assert indep == [
        "(3,11) v22",
        "(3,11) v49",
        "(3,12) from (3,11) v22",
        "(3,12) from (3,11) v49",
    ]


def test_every_verdict_carries_independent_scope_and_a_correction(artifact):
    for name, v in artifact["scorecard"].items():
        assert "independent_scope" in v, name
        assert "independent_distinct_states" in v, name
        assert v["independent_scope"] <= v.get("scope", 38), name
        assert v["original_verdict"] in ("PASS", "FAIL"), name
        assert v["corrected_reading"], name
        assert v["correction_note"], name


def test_the_v43_anomaly_is_its_donor(artifact):
    """(3,11) v43 is (3,10) v108 padded, identical on every structural field."""
    per = {p["label"]: p for p in artifact["provenance_independence"]["per_measurement"]}
    assert per["(3,11) v43"]["root_origin"] == "(3,10) v108"
    assert not per["(3,11) v43"]["independent_of_rank10_corpus"]


def test_corrections_downgrade_the_three_weak_verdicts(artifact):
    sc = artifact["scorecard"]
    assert sc["B3_two_fibers_agree"]["corrected_reading"] == "ALREADY FALSIFIED AT RANK 10"
    assert sc["B4_no_cancellation_regime"]["corrected_reading"] == "UNINFORMATIVE"
    assert sc["B5_no_sparsification"]["corrected_reading"] == "INHERITED"
    assert sc["B1_unified_dimension_formula"]["corrected_reading"] == "PASS, SCOPE OVERSTATED"
    # the two that survive intact
    assert sc["B2_padding_invariance"]["corrected_reading"] == "PASS"
    assert sc["B6_numerics_at_new_rank"]["corrected_reading"] == "PASS"
