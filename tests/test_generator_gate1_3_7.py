"""Guard the (3,7) exhaustiveness gate and its scorecard.

The gate claims P = H for (3,7) by the inner/outer closure. The claim is a
composition of four exact facts, so the tests check each separately rather
than only the composed verdict, and they check the two controls, because a
closure that came from the ordered slice alone would prove nothing about the
Ressayre rows.
"""

import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "generator_gate1_3_7.json"

MODULE_PATH = ROOT / "scripts" / "generator_gate1_3_7.py"
if MODULE_PATH.exists():
    SPEC = importlib.util.spec_from_file_location("gate1", MODULE_PATH)
    mod = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(mod)
else:  # pragma: no cover
    mod = None


def artifact():
    if not ARTIFACT.exists():
        pytest.skip("gate 1 artifact not present")
    return json.loads(ARTIFACT.read_text())


def test_every_prediction_passed():
    scorecard = artifact()["scorecard"]
    failed = {k: v for k, v in scorecard.items() if v != "PASS"}
    assert not failed, f"gate 1 predictions not passing: {failed}"


def test_the_closure_is_a_composition_of_the_other_four():
    """G5 must not be assertable on its own."""
    s = artifact()["scorecard"]
    composed = all(
        s[k] == "PASS"
        for k in (
            "G1_ten_vertices",
            "G2_matches_reference",
            "G3_every_vertex_attained",
            "G4_every_row_certified_on_replay",
        )
    )
    assert (s["G5_closure_exact"] == "PASS") == composed


def test_the_vertex_count_matches_the_published_reference():
    data = artifact()
    assert data["vertex_count"] == data["reference_vertex_count"] == 10


def test_every_certificate_replayed_to_its_recorded_value():
    """Replayed by Fraction elimination, not the generator's Bareiss routine."""
    replays = artifact()["certificate_replays"]
    assert len(replays) == 4
    for r in replays:
        assert r["trace_condition"]
        assert r["determinant"] != 0
        assert r["determinant_matches_record"]
        assert r["on_matches_record"]
        assert r["below_matches_record"]
        assert r["roots_match_record"]


def test_every_vertex_has_an_exact_state_with_a_recomputed_spectrum():
    """G3 is non-circular: the spectrum comes from the support, not integer_form."""
    attain = artifact()["attainability"]
    assert len(attain) == 10
    for record in attain:
        assert record["ok"]
        assert record["tierB"] in {"EXACT", "EXACT-CONSTR"}
        assert record["recomputed_spectrum_matches"]


# ------------------------------------------------------------------ controls


def test_the_ressayre_rows_actually_cut_the_slice():
    """Without this the closure could be a property of the ordered slice."""
    data = artifact()
    assert data["ordered_slice_only_vertex_count"] == 13
    assert data["ordered_slice_only_vertex_count"] > data["vertex_count"]


def test_each_row_changes_the_polytope():
    drops = artifact()["single_row_drops"]
    assert len(drops) == 4
    for drop in drops:
        assert drop["load_bearing"]
        assert drop["vertex_set_changed"]


def test_the_load_bearing_test_is_a_set_comparison_not_a_count():
    """Row 3 is the case that makes the distinction matter.

    Dropping it leaves the vertex count at 10 while changing which ten they
    are. A count-based control would have scored it redundant, which is what a
    first version of the script did.
    """
    drops = {d["dropped_row"]: d for d in artifact()["single_row_drops"]}
    row3 = drops[3]
    assert row3["vertices"] == 10, "the count is unchanged for this row"
    assert row3["vertex_set_changed"], "but the vertex set is not"
    assert row3["load_bearing"]


# ----------------------------------------------------------------- nonclaims


def test_the_artifact_does_not_claim_candidate_exhaustiveness():
    nonclaims = " ".join(artifact()["nonclaims"])
    assert "exhaustiveness" in nonclaims
    assert "bypassed" in nonclaims


def test_the_manifest_completeness_field_is_untouched():
    """The closure route does not answer the enumeration question."""
    manifest = json.loads((DATA / "stage1_ressayre_3_7.json").read_text())
    assert "candidate_exhaustiveness" in manifest["proof_status"]


def test_the_artifact_records_that_the_states_are_not_independent():
    nonclaims = " ".join(artifact()["nonclaims"])
    assert "independent" in nonclaims


def test_the_gate_points_at_its_prereg():
    assert artifact()["prereg"] == "docs/prereg_generator_gate1_3_7.md"
    assert (ROOT / "docs" / "prereg_generator_gate1_3_7.md").exists()


# -------------------------------------------------- recomputed independently


@pytest.mark.skipif(mod is None, reason="gate script not present")
def test_the_ordered_slice_alone_really_has_thirteen_vertices():
    """Recompute rather than trust the artifact's own number."""
    assert len(mod.enumerate_vertices([])) == 13


@pytest.mark.skipif(mod is None, reason="gate script not present")
def test_the_enumerated_vertices_are_the_reference_vertices():
    manifest = json.loads((DATA / "stage1_ressayre_3_7.json").read_text())
    rows = mod.ressayre_rows(manifest)
    verts = set(mod.enumerate_vertices(rows))
    reference = {
        tuple(F(x) for x in v["spectrum"])
        for v in json.loads((DATA / "vertices" / "vertices_3_7.json").read_text())
    }
    assert verts == reference
