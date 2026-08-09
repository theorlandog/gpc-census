"""Guards for the spread-cost table.

The table measures how fast ``D(W)`` climbs. It never bounded the spread, and
since Theorem S proves a bound by an unrelated route it does not have to. The
thing most worth protecting is still its label: a ceiling read as a floor would
turn "the rate is visible" into "the bound is proved here", and the bound that
IS proved is exponential, so borrowing this table's linear-looking rate to stand
in for it would be the same error wearing a theorem's clothes.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "data" / "profile_spread_cost.json"
SCRIPT = ROOT / "scripts" / "profile_spread_cost.py"


def _artifact() -> dict:
    return json.loads(ARTIFACT.read_text())


def _module():
    if not SCRIPT.exists():  # sdist may omit scripts
        pytest.skip("generator script not present")
    spec = importlib.util.spec_from_file_location("profile_spread_cost", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_cost_climbs_with_spread_which_is_the_whole_point():
    artifact = _artifact()
    rows = artifact["rows"]
    assert len(rows) >= 20
    costs = [row["cost"] for row in rows]
    assert artifact["cost_is_nondecreasing"]
    assert costs == sorted(costs)
    # flat would leave the rate claim with no content at all
    assert costs[-1] > costs[0]


def test_the_ceiling_is_labelled_as_a_ceiling():
    labels = _artifact()["proof_labels"]
    assert labels["restricted_cost_is_an_upper_bound_on_D"].startswith("PROVED")
    assert labels["cost_is_nondecreasing_in_spread"].startswith(
        "EXACTLY VERIFIED ON FINITE POPULATION"
    )
    # the spread bound is proved, but NOT by this table, and the label has to
    # keep those two apart
    spread = labels["spread_is_bounded_at_each_rank"]
    assert spread.startswith("PROVED by Theorem S")
    assert "independently of this table" in spread
    assert any("does not bound the spread" in claim
               for claim in _artifact()["nonclaims"])


def test_agrees_with_the_frontier_enumeration():
    """A spread realized at rank d forces D(W) <= d, and the two runs agree."""
    check = _artifact()["frontier_cross_check"]
    assert check["available"]
    assert check["violations"] == []
    realized = check["realized_max_spread_by_rank"]
    assert {int(k): v for k, v in realized.items()} == {
        6: 3, 7: 4, 8: 6, 9: 10, 10: 14
    }


def test_every_row_is_a_replayable_witness():
    """Recompute the cost of each stored witness rather than reading it back."""
    module = _module()
    n = _artifact()["parameters"]["particle_number"]
    for row in _artifact()["rows"]:
        levels = tuple(row["levels"])
        target = row["level_target"]
        r = row["classes"]
        assert len(levels) == r
        assert levels[0] == 0 and levels[-1] == row["spread"]
        patterns = [
            pattern
            for pattern in module.occupancy_patterns(n, r)
            if sum(c * lev for c, lev in zip(pattern, levels, strict=True)) == target
        ]
        assert module._rank_mod(patterns, r) == r - 1
        assert module._min_cost(patterns, r, n, 4 * r) == row["cost"]


def test_a_cheaper_multiplicity_vector_really_does_not_exist():
    """Spot-check minimality directly: one below the reported cost must fail."""
    module = _module()
    n = _artifact()["parameters"]["particle_number"]
    checked = 0
    for row in _artifact()["rows"][:8]:
        levels = tuple(row["levels"])
        r = row["classes"]
        patterns = [
            pattern
            for pattern in module.occupancy_patterns(n, r)
            if sum(c * lev for c, lev in zip(pattern, levels, strict=True))
            == row["level_target"]
        ]
        below = row["cost"] - 1
        if below < r:
            continue
        for mhat in module._multiplicity_vectors(r, below, n + 1):
            assert not module._admissible(patterns, mhat, r)
        checked += 1
    assert checked >= 4
