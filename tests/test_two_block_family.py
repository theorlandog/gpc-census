"""The two-block family artifact and the cross-block classification result.

results/data/two_block_family.json certifies docs/two_block_family.md. These
tests pin the scored verdicts so a later regeneration cannot silently lose a
falsification, and recompute the load-bearing result (a cross-block channel is
necessary for the two fibers to disagree) from the ledger with independent
block and channel bookkeeping.
"""
import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"

ARTIFACT = DATA / "two_block_family.json"


@pytest.fixture(scope="module")
def artifact():
    return json.loads(ARTIFACT.read_text())


def test_the_family_has_exactly_four_members(artifact):
    fam = artifact["family"]
    assert len(fam) == 4
    assert {r["label"] for r in fam} == {
        "(3,10) v103",
        "(3,10) v89",
        "(3,10) v108",
        "(3,8) v32",
    }
    for r in fam:
        assert r["blocks"] == 2
        assert r["classified"] == "INTERFERENCE"


def test_family_channel_split_matches_the_mechanism(artifact):
    """v89 is rigid in both fibers because its only channel is within-block."""
    fam = {r["label"]: r for r in artifact["family"]}
    v89 = fam["(3,10) v89"]
    assert (v89["channels_within_block"], v89["channels_cross_block"]) == (1, 0)
    assert v89["gap"] == 0

    v103 = fam["(3,10) v103"]
    assert (v103["channels_within_block"], v103["channels_cross_block"]) == (2, 3)
    assert (v103["fixed_spectrum_tangent_dim"], v103["fixed_rho_tangent_dim"]) == (6, 0)

    v108 = fam["(3,10) v108"]
    assert (v108["channels_within_block"], v108["channels_cross_block"]) == (0, 1)
    assert v108["gap"] == 1

    v32 = fam["(3,8) v32"]
    assert v32["channels_cross_block"] == 1
    assert v32["incidence_kernel_dim"] == 0  # nothing for the condition to cut
    assert v32["gap"] == 0


def test_cross_block_channel_is_necessary_for_a_gap(artifact):
    """The load-bearing result, recomputed here from the per-state rows."""
    rows = artifact["states"]
    assert len(rows) == 799
    offenders = [r["label"] for r in rows if r["gap"] > 0 and r["X"] == 0]
    assert offenders == []
    assert sum(1 for r in rows if r["gap"] > 0) == 73
    assert artifact["scorecard"]["P_INV_3_forward"]["verdict"] == "PASS"
    assert artifact["scorecard"]["P_INV_3_forward"]["exceptions"] == 0


def test_cross_block_channel_is_not_sufficient(artifact):
    """The converse fails, and every exception is rigid for lack of freedom."""
    rows = artifact["states"]
    exceptions = [r for r in rows if r["X"] > 0 and r["gap"] == 0]
    assert len(exceptions) == 82
    for r in exceptions:
        assert r["incidence_kernel_dim"] == 0
        assert r["fixed_spectrum_tangent_dim"] == 0
        assert r["fixed_rho_tangent_dim"] == 0
    assert artifact["scorecard"]["P_INV_3_reverse"]["verdict"] == "FAIL"


def test_the_counting_rules_are_falsified(artifact):
    """Pinned so the falsification is not lost in a later regeneration."""
    sc = artifact["scorecard"]
    assert sc["P_INV_1"]["verdict"] == "FAIL"
    assert sc["P_INV_2"]["verdict"] == "FAIL"
    assert sc["P_INV_3"]["verdict"] == "FAIL"


def test_inversions_are_fourteen_not_one(artifact):
    """The handoff's premise said sample size one; it is 14 in 7 classes."""
    rows = artifact["states"]
    inv = [
        r
        for r in rows
        if r["fixed_rho_tangent_dim"] == 0 and r["fixed_spectrum_tangent_dim"] > 0
    ]
    assert len(inv) == 14
    assert len({r["transport_class"] for r in inv}) == 7
    assert all(r["X"] > 0 for r in inv)
    assert artifact["scorecard"]["P_INV_4"]["verdict"] == "PASS"
    assert artifact["scorecard"]["totals"]["inversion_states"] == 14


def test_empty_antecedent_is_uninformative_not_pass(artifact):
    """Standing rule R3: a prediction that could not fail does not PASS."""
    assert artifact["scorecard"]["P_INV_5"]["verdict"] == "UNINFORMATIVE"


def test_post_hoc_restriction_is_labeled(artifact):
    r = artifact["scorecard"]["P_INV_3_restricted"]
    assert r["verdict"] == "PASS"
    assert "POST-HOC" in r["registered_expectation"]
    assert "POST-HOC" in r["note"]


def test_v103_tangent_splits_three_weight_three_phase(artifact):
    """The lemma's formula (5-2) + (14-9-2) = 3 + 3, verified numerically."""
    split = {s["label"]: s for s in artifact["tangent_sector_split"]}
    v103 = split["(3,10) v103"]
    assert v103["tangent_dim"] == 6
    assert v103["weight_sector_directions"] == 3
    assert v103["phase_sector_beyond_gauge"] == 3
    assert v103["sectors_decouple"]

    # v108 is genuinely complex, so the real-state split does not apply
    v108 = split["(3,10) v108"]
    assert not v108["sectors_decouple"]
    assert v108["phase_sector_beyond_gauge"] is None


def test_two_block_is_not_the_operative_invariant(artifact):
    s = artifact["stratification"]
    assert s["two_block_designs_with_any_channel"] == 0
    assert s["two_block_with_cross_block_channel"] == 3
    assert s["two_block_vertices"] == 64
    # the gap population is spread across block counts, not concentrated at 2
    gap_dist = s["block_count_distribution"]["two_fiber_gap"]
    assert gap_dist["2"] == 2
    assert sum(gap_dist.values()) == 73
    assert max(gap_dist, key=lambda k: gap_dist[k]) == "4"


def test_independent_scope_is_reported_in_transport_classes(artifact):
    """Standing rule R2 applied to this scorecard, not just the bracket one."""
    sc = artifact["scorecard"]
    assert sc["totals"]["transport_classes"] == 403
    assert sc["totals"]["transport_classes_with_gap"] == 41
    for key in ("P_INV_1", "P_INV_3_forward", "P_INV_4"):
        assert sc[key]["independent_scope_transport_classes"] <= sc[key]["scope"]
