"""Guards for the stabilizer-character explanation of the quantization period.

The expensive exact recomputation is re-run only for `(3,6)#3`, the vertex
whose period-4 mechanism is classical, which is enough to catch a broken
engine without making the suite slow.
"""
from __future__ import annotations

import json
import pathlib

import pytest

from gpc_census.symplectic import (
    levi_permutation_stabilizer,
    predicted_quantization_period,
    spectrum_blocks,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" if (ROOT / "tests" / "data").exists() else ROOT / "results" / "data"
ARTIFACT = DATA / "quantization_characters.json"
MEASURED = DATA / "quantization_multiplicities.json"


@pytest.fixture(scope="module")
def art():
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def states():
    return {(r["system"], r["index"]): r
            for r in (json.loads(x)
                      for x in (DATA / "states.jsonl").read_text().splitlines())}


def test_covers_every_measured_vertex(art):
    measured = json.loads(MEASURED.read_text())["rows"]
    assert len(art["rows"]) == len(measured) == 14
    keys = {(r["system"], r["index"]) for r in art["rows"]}
    assert keys == {(r["system"], r["index"]) for r in measured}
    assert all(r["evidence"] == "exact" for r in art["rows"])


def test_predicted_period_divides_measured(art):
    """Forced by the derivation: a subgroup can only under-predict."""
    bad = [f"{r['system']}#{r['index']}" for r in art["rows"] if not r["divides"]]
    assert not bad, bad


def test_predicted_period_equals_measured(art):
    bad = [f"{r['system']}#{r['index']}" for r in art["rows"] if not r["matches"]]
    assert not bad, bad


def test_period_doubling_is_exactly_the_obstructed_set(art):
    doubled = {f"{r['system']}#{r['index']}" for r in art["rows"]
               if r["measured_period"] != r["denominator"]}
    carriers = {f"{r['system']}#{r['index']}" for r in art["rows"]
                if r["period_doubling_elements"]}
    assert doubled == carriers
    assert doubled == {"(3,6)#3", "(3,7)#5", "(3,7)#7", "(3,7)#8"}


def test_block_exponents_use_the_block_occupation(art):
    """Regression: the exponent on block i is k*lambda_i, not k*form[i].

    Indexing `integer_form` by the block counter instead of by an orbital of
    that block manufactured obstructions at three vertices and broke the forced
    divisibility gate. Every block here must carry the occupation shared by its
    own orbitals.
    """
    for r in art["rows"]:
        blocks = spectrum_blocks(r["integer_form"])
        assert [len(b) for b in blocks] == r["spectrum_blocks"], r
        for b in blocks:
            assert len({r["integer_form"][o] for o in b}) == 1, r


def test_scorecard_verdicts_are_the_recorded_ones(art):
    card = art["scorecard"]
    for key in ("P1_predicted_divides_measured", "P2_predicted_equals_measured",
                "P3_doubling_iff_permutation_obstruction", "P4_identity_always_present"):
        assert card[key]["verdict"] == "PASS", key


def test_padding_carries_the_mechanism(art):
    """(3,7)#7 is (3,6)#3 with an empty orbital; its obstruction is the padded one."""
    src = next(r for r in art["rows"] if r["system"] == "(3,6)" and r["index"] == 3)
    dst = next(r for r in art["rows"] if r["system"] == "(3,7)" and r["index"] == 7)
    padded = {tuple(e["sigma"]) + (6,) for e in src["elements"]
              if e["least_admissible_k"] != src["denominator"]}
    got = {tuple(e["sigma"]) for e in dst["elements"]
           if e["least_admissible_k"] != dst["denominator"]}
    assert padded == got
    assert len(padded) == 18


def test_recompute_rank_six_uniform_vertex(art, states):
    row = next(r for r in art["rows"] if r["system"] == "(3,6)" and r["index"] == 3)
    elements = levi_permutation_stabilizer(states[("(3,6)", 3)])
    assert len(elements) == row["stabilizer_elements"] == 36
    assert predicted_quantization_period(elements, 2) == row["predicted_period"] == 4
    swap = next(e for e in elements if e["sigma"] == [3, 4, 5, 0, 1, 2])
    assert swap["c_pretty"] == "1"
    assert swap["block_determinants"] == [-1]
    assert swap["least_admissible_k"] == 4


def test_census_predictions_cover_every_vertex(art):
    rows = art["census_predictions"]
    assert len(rows) == 799
    assert len({(r["system"], r["index"]) for r in rows}) == 799
    assert all(r["evidence"] == "exact" for r in rows)
    assert sum(1 for r in rows if r["verified"]) == 14


def test_census_predictions_are_all_resolved(art):
    unresolved = [f"{r['system']}#{r['index']}" for r in art["census_predictions"]
                  if r["predicted_period"] is None]
    assert not unresolved, unresolved


def test_the_obstruction_is_exactly_a_two_group(art):
    """Across all 799 vertices the predicted period is q or 2q, never anything else."""
    ratios = set()
    for r in art["census_predictions"]:
        p, q = r["predicted_period"], r["denominator"]
        assert p % q == 0, r
        ratios.add(p // q)
    assert ratios == {1, 2}
    doubled = [r for r in art["census_predictions"] if r["predicted_period_exceeds_denominator"]]
    assert len(doubled) == 54


def test_verified_rows_agree_with_the_census_pass(art):
    census = {(r["system"], r["index"]): r for r in art["census_predictions"]}
    for row in art["rows"]:
        got = census[(row["system"], row["index"])]
        assert got["predicted_period"] == row["predicted_period"], row
        assert got["verified"] is True
