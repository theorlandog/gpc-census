"""The (3,11) closure: the four claimed level-five rows, and why to trust them.

The tests are split the way the artifact is. Calibration tests check that the
screen is sound where the answer is already published; negative-control tests
check that it is not vacuously sound; rank-eleven tests check the result and
its consequences. The last group asserts the nonclaims, because the value of
this artifact depends on it not overreaching.
"""

import importlib.util
import json
import subprocess
import sys
from fractions import Fraction as F
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "level5_ressayre_3_11.json"

MODULE_PATH = ROOT / "scripts" / "level5_ressayre_3_11.py"
SPEC = importlib.util.spec_from_file_location("level5", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)

FAMILIES = ("AK-RMK-4.2.1", "KLY09-EXT-1", "KLY09-EXT-2", "KLY09-EXT-3")


def artifact():
    return json.loads(ARTIFACT.read_text())


def calibration(rank):
    return next(c for c in artifact()["calibration"] if c["rank"] == rank)


# --------------------------------------------------------------- calibration


def test_the_screen_certifies_nothing_false_at_the_published_ranks():
    """The soundness test. Complete over the polytope, not sampled."""
    for rank in (9, 10):
        entry = calibration(rank)
        assert entry["certified_rows"] > 0
        assert entry["certified_rows_invalid_on_published_polytope"] == 0
        assert entry["invalid_rows"] == []
        assert entry["sound"]


def test_every_certified_row_is_valid_by_its_own_recorded_lp_maximum():
    """Recheck the per-row verdicts rather than trusting the summary count."""
    for rank in (9, 10):
        for row in calibration(rank)["certified"]:
            assert row["lp_maximum"] <= row["rhs"] + 1e-7
            assert row["valid"]


def test_the_calibration_actually_exercised_the_screen():
    """A screen that certified nothing would pass the soundness test vacuously."""
    assert calibration(9)["certified_rows"] == 19
    assert calibration(10)["certified_rows"] == 31
    assert calibration(9)["rows_screened"] == 984
    assert calibration(10)["rows_screened"] == 1914


def test_some_certified_rows_are_published_facets_and_most_are_not():
    """Validity is the claim; facetness is explicitly not.

    A screen returning only published facets would be suspicious in the other
    direction, since it would suggest it were reading the table.
    """
    for rank in (9, 10):
        entry = calibration(rank)
        assert 0 < entry["certified_verbatim_in_published_table"] < entry["certified_rows"]


# ----------------------------------------------------------- negative control


def test_the_calibration_screened_many_demonstrably_false_rows():
    for rank in (9, 10):
        control = calibration(rank)["negative_control"]
        assert control["demonstrably_false_rows_screened"] > 500


def test_no_false_row_was_ever_certified():
    for rank in (9, 10):
        assert calibration(rank)["negative_control"]["false_rows_certified"] == 0


def test_false_rows_reached_the_tangent_determinant_and_died_there():
    """Without this the control would only exercise the cheap shortcuts.

    ``inadmissible`` and ``trace_rejected`` are decided before the determinant
    is formed, so a control in which every false row took one of those routes
    would say nothing about the determinant step.
    """
    reached = sum(
        calibration(rank)["negative_control"]["false_rows_reaching_the_tangent_determinant"]
        for rank in (9, 10)
    )
    assert reached > 0
    for rank in (9, 10):
        for row in calibration(rank)["negative_control"]["rows_reaching_the_determinant"]:
            assert row["status"] in mod.REACHED_DETERMINANT
            assert row["status"] != "certified"
            assert row["lp_maximum"] > row["rhs"]


def test_the_determinant_route_is_a_minority_of_the_rejections():
    """Honest bookkeeping: the trace count does most of the rejecting.

    This is recorded so the artifact cannot be read as claiming the
    determinant was stress tested on hundreds of false rows.
    """
    for rank in (9, 10):
        control = calibration(rank)["negative_control"]
        assert (
            control["false_row_statuses"]["trace_rejected"]
            > control["false_rows_reaching_the_tangent_determinant"]
        )


# ------------------------------------------------------------ randomized sweep


def test_the_sweep_tests_general_coefficients_not_just_subset_rows():
    """The calibration rows are all zero-one; this is the only leg that is not."""
    sweep = artifact()["randomized_sweep"]
    assert sweep["coefficient_alphabet"] == [-1, 0, 1, 2]
    assert sweep["rows_screened"] == mod.SWEEP_ROWS
    assert sweep["seed"] == mod.SWEEP_SEED


def test_the_sweep_certifies_nothing_false():
    sweep = artifact()["randomized_sweep"]
    assert sweep["certified_rows"] > 0
    assert sweep["certified_rows_invalid_on_published_polytope"] == 0
    assert sweep["invalid_rows"] == []
    assert sweep["negative_control"]["false_rows_certified"] == 0
    assert sweep["sound"]
    for row in sweep["certified"]:
        assert row["lp_maximum"] <= row["rhs"] + 1e-7


def test_the_sweep_is_mostly_false_rows_so_the_control_has_teeth():
    control = artifact()["randomized_sweep"]["negative_control"]
    assert control["demonstrably_false_rows_screened"] > 1000


# ---------------------------------------------------------------- rank eleven


def test_all_four_claimed_rows_are_certified():
    rows = artifact()["rank_eleven"]["rows"]
    assert {r["family"] for r in rows} == set(FAMILIES)
    for row in rows:
        assert row["screen_status"] == "certified"
        assert row["census_tier_before"] == "CLAIMED"
        assert row["census_tier_after"] == "PROVED"
        assert row["rhs"] == 2
        assert len(row["subset"]) == 5
        assert row["subset"][-1] == 11


def test_each_row_is_a_second_level_quadruple_plus_the_last_orbital():
    quadruples = {(2, 3, 4, 5), (1, 3, 4, 6), (1, 2, 5, 6), (1, 2, 4, 7)}
    assert {tuple(r["subset"][:4]) for r in artifact()["rank_eleven"]["rows"]} == quadruples


def test_no_certified_row_is_violated_by_an_already_attained_spectrum():
    """The only soundness check available at a rank whose polytope is unknown."""
    eleven = artifact()["rank_eleven"]
    assert eleven["attained_vertices_checked"] == 19
    assert eleven["certified_rows_violated_by_an_attained_vertex"] == 0
    assert eleven["violations"] == []


def test_the_certified_rows_hold_at_the_slater_and_uniform_points():
    """Two attainable spectra the bracket does not have to supply."""
    points = ([F(1)] * 3 + [F(0)] * 8, [F(3, 11)] * 11)
    for _name, subset in mod.CLAIMED:
        for lam in points:
            assert sum(lam[j - 1] for j in subset) <= 2


def test_every_remaining_open_candidate_is_cut_off():
    records = artifact()["rank_eleven"]["open_candidates"]
    assert sorted(r["index"] for r in records) == [23, 26, 34, 44]
    for record in records:
        assert record["refuted"]
        assert record["cut_by"]
    assert artifact()["rank_eleven"]["all_four_open_candidates_refuted"]


def test_candidate_44_is_cut_by_exactly_the_one_row_it_stood_against():
    """44 was the cheapest falsification target in the census; it is closed.

    It violated a single claimed row, so attaining it would have refuted
    KLY09-EXT-1. That row is now certified, and the design route was already
    closed by exact counting, so nothing attains 44.
    """
    record = next(
        r for r in artifact()["rank_eleven"]["open_candidates"] if r["index"] == 44
    )
    assert [c["family"] for c in record["cut_by"]] == ["KLY09-EXT-1"]
    assert F(record["cut_by"][0]["value"]) == F(2) + F(1, 12)


def test_the_bracket_settles_at_nineteen_true_and_thirty_one_refuted():
    conclusion = artifact()["conclusion"]
    assert conclusion["screen_sound_at_published_ranks"]
    assert conclusion["negative_control_passed"]
    assert conclusion["four_claimed_rows_certified"]
    assert conclusion["bracket_settles"]
    assert conclusion["tally_after"] == {"TRUE": 19, "REFUTED": 31, "OPEN": 0}
    settlement = json.loads((ROOT / "docs" / "bracket_3_11_settlement.json").read_text())
    assert len(settlement["candidates"]) == 19 + 31


# ---------------------------------------------------- independent replay


VERIFIER = ROOT / "scripts" / "verify_level5_ressayre_3_11_standalone.py"


def run_verifier(artifact_path):
    return subprocess.run(
        [sys.executable, str(VERIFIER), "--artifact", str(artifact_path)],
        capture_output=True,
        text=True,
        check=False,
    )


def test_the_certificates_replay_under_a_verifier_that_shares_no_code():
    """The repository's validation law: a claim is not a certificate.

    The rank-7 checkpoint puts the Ressayre determinants outside its
    independent verifier's trust boundary, so without this the rank-11 result
    would rest on the generator checking itself. The standalone verifier uses
    the standard library only, derives the fermionic sign by counting
    inversions rather than by tracking creation and annihilation positions,
    and recomputes the determinant by Fraction elimination rather than by the
    generator's Bareiss routine.
    """
    if not VERIFIER.exists():
        pytest.skip("standalone verifier not present (sdist build)")
    result = run_verifier(ARTIFACT)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "0 failed" in result.stdout


def test_the_verifier_rejects_a_tampered_determinant(tmp_path):
    """A verifier that cannot fail proves nothing."""
    if not VERIFIER.exists():
        pytest.skip("standalone verifier not present (sdist build)")
    payload = artifact()
    row = next(r for r in payload["rank_eleven"]["rows"] if "certificate" in r)
    row["certificate"]["determinant_value"] += 1
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(payload))
    result = run_verifier(tampered)
    assert result.returncode == 1
    assert "determinant replays to" in result.stdout


def test_the_verifier_rejects_a_tampered_hyperplane(tmp_path):
    if not VERIFIER.exists():
        pytest.skip("standalone verifier not present (sdist build)")
    payload = artifact()
    row = next(r for r in payload["rank_eleven"]["rows"] if "certificate" in r)
    row["certificate"]["h"][0] += 1
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(payload))
    result = run_verifier(tampered)
    assert result.returncode == 1
    assert "centered transform of tau" in result.stdout


def test_every_certificate_carries_a_replayable_payload():
    for row in artifact()["rank_eleven"]["rows"]:
        assert row["screen_status"] == "certified"
        certificate = row["certificate"]
        assert len(certificate["h"]) == 11
        assert len(certificate["on_indices"]) == len(certificate["determinant_evaluation"])
        assert len(certificate["below_indices"]) == len(certificate["negative_roots"])
        assert certificate["determinant_value"] != 0
        assert row["below_weight_count"] == row["negative_root_count"]


# ------------------------------------------------------------------ nonclaims


def test_the_tau_convention_is_the_one_the_pipeline_decodes():
    """tau_i = N * [i in A] - k, checked against the module's own encoder."""
    assert mod.tau_of((1, 2, 4, 7, 11), 2, 11) == (1, 1, -2, 1, -2, -2, 1, -2, -2, -2, 1)
    assert mod.tau_of((2, 3, 4, 5), 2, 7) == (-2, 1, 1, 1, 1, -2, -2)


def test_the_artifact_claims_validity_only():
    nonclaims = artifact()["nonclaims"]
    assert any("facetness" in n for n in nonclaims)
    assert any("complete" in n for n in nonclaims)
    assert any("rank 11" in n for n in nonclaims)
    assert any("sample, not an enumeration" in n for n in nonclaims)


def test_the_settlement_artifact_is_left_as_the_record_of_its_own_method():
    """The closure is a new layer, not a rewrite of the contraction audit.

    scripts/settle_bracket_3_11.py still reproduces 19/27/4 by its own method,
    so that artifact must keep saying so.
    """
    settlement = json.loads((ROOT / "docs" / "bracket_3_11_settlement.json").read_text())
    still_open = sorted(
        c["index"] for c in settlement["candidates"] if c.get("verdict") == "OPEN"
    )
    assert still_open == [23, 26, 34, 44]
