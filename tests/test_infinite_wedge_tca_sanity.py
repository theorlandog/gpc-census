"""Guards for the stable-semigroup bridge and the cross-rank facet measurement.

The load-bearing claims are re-derived here rather than read off the artifact:
the plethysm convention at a handful of weights, the four non-saturation
witnesses, one named zero-padding counterexample, and the Farkas implication at
the cheapest consecutive pair. Everything else is checked for internal
consistency, which is what catches a rewritten artifact.
"""
from __future__ import annotations

import json
import pathlib
from fractions import Fraction

import pytest

from gpc_census.constraints import constraints
from gpc_census.plethysm import hm_en, mn

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = (ROOT / "tests" / "data" if (ROOT / "tests" / "data").exists()
        else ROOT / "results" / "data")
SANITY = DATA / "infinite_wedge_tca_sanity.json"
DEPENDENCY = DATA / "infinite_wedge_tca_dependency.json"


@pytest.fixture(scope="module")
def art():
    return json.loads(SANITY.read_text())


@pytest.fixture(scope="module")
def graph():
    return json.loads(DEPENDENCY.read_text())


def _plethysm(lam, k, n=3):
    lam = tuple(x for x in lam if x > 0)
    if sum(lam) != n * k:
        return 0
    total = Fraction(0)
    for mu, coefficient in hm_en(k, n).items():
        character = mn(lam, tuple(sorted(mu, reverse=True)))
        if character:
            total += coefficient * character
    assert total.denominator == 1
    return int(total)


def _vertices(n, d):
    rows = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
    return [tuple(Fraction(int(v), int(rec["denominator"] or 1))
                  for v in rec["integer_form"]) for rec in rows]


# --------------------------------------------------------------------------
# convention
# --------------------------------------------------------------------------

def test_every_stored_multiplicity_row_replays_by_plethysm(art):
    """Spot check the convention on the cheap half of the rows."""
    checked = 0
    for row in art["convention"]["rows"]:
        if row["k"] > 6:
            continue
        n = int(row["system"].strip("()").split(",")[0])
        assert _plethysm(row["weight"], row["k"], n) == row["stored"], row
        checked += 1
    assert checked >= 30
    assert art["convention"]["all_agree"]
    assert not art["convention"]["mismatches"]


def test_the_convention_control_is_not_vacuous(art):
    """A gate every convention passes proves nothing; this one discriminates."""
    assert art["convention"]["control_is_load_bearing"]
    assert art["convention"]["control_pairs_where_a_wrong_convention_differs"] > 0


# --------------------------------------------------------------------------
# recovery
# --------------------------------------------------------------------------

def test_recovery_closes_both_inclusions(art):
    systems = {r["system"]: r for r in art["recovery"]["systems"]}
    assert set(systems) == {"(3,6)", "(3,7)"}
    for rec in systems.values():
        assert rec["points_violating_a_published_row"] == 0, rec["violations"]
        assert not rec["missing_vertices"], rec["missing_vertices"]
        assert rec["hull_equals_polytope"]
    # (3,7) genuinely needs degree 8; m <= 7 would miss a vertex.
    assert max(systems["(3,7)"]["first_degree_of_each_vertex"].values()) == 8


# --------------------------------------------------------------------------
# saturation
# --------------------------------------------------------------------------

def test_the_four_non_saturation_witnesses_replay(art):
    witnesses = art["saturation"]["non_saturation_witnesses"]
    assert [w["vertex"] for w in witnesses] == [
        "(3,6)#3", "(3,7)#5", "(3,7)#7", "(3,7)#8"]
    for witness in witnesses:
        n = int(witness["vertex"].strip("(").split(",")[0])
        q = witness["denominator"]
        # in the cone by construction, and outside the semigroup by measurement
        assert _plethysm(witness["integer_form"], q, n) == 0, witness
        assert witness["measured_period"] % q == 0
        assert witness["measured_period"] > q


def test_a_saturated_vertex_is_actually_saturated(art):
    """Control: the witnesses are not an artefact of every vertex vanishing."""
    assert art["saturation"]["saturated_at_k_equals_q"] == 10
    assert _plethysm([1, 1, 1, 0, 0, 0], 1, 3) == 1


# --------------------------------------------------------------------------
# padding, which is the campaign's refutation
# --------------------------------------------------------------------------

def test_named_zero_padding_counterexample():
    """`lambda_1 + lambda_8 <= 1` is a (3,8) row and fails at (3,9)."""
    row = ([1, 0, 0, 0, 0, 0, 0, 1], 1)
    assert row in [([int(c) for c in r["coeffs"]], int(r["rhs"]))
                   for r in constraints(3, 8)["inequalities"]]
    padded = row[0] + [0]
    worst = min(Fraction(row[1]) - sum(Fraction(c) * x for c, x in zip(padded, v))
                for v in _vertices(3, 9))
    assert worst == Fraction(-1, 4)


def test_padding_fails_at_the_pairs_the_note_reports(art):
    valid = {f"{p['from']} -> {p['to']}": (p["rows_still_valid"], p["rows_tested"])
             for p in art["padding"]["pairs"]}
    assert valid == {
        "(3,6) -> (3,7)": (0, 4),
        "(3,7) -> (3,8)": (4, 4),
        "(3,8) -> (3,9)": (10, 31),
        "(3,9) -> (3,10)": (44, 52),
        "(4,8) -> (4,9)": (8, 15),
        "(4,9) -> (4,10)": (19, 60),
    }
    assert art["scorecard"][
        "P4_zero_padding_preserves_validity_above_rank_6"]["verdict"] == "FAIL"


def test_no_order_preserving_embedding_rescues_it(art):
    """The refutation is not an artefact of choosing the trailing embedding."""
    unrescued = [p for p in art["insertion_post_hoc"]["pairs"]
                 if not p["every_row_has_some_valid_embedding"]]
    assert [p["from"] + " -> " + p["to"] for p in unrescued] == [
        "(3,6) -> (3,7)", "(3,8) -> (3,9)", "(3,9) -> (3,10)",
        "(4,8) -> (4,9)", "(4,9) -> (4,10)"]
    assert art["insertion_post_hoc"]["status"].startswith("POST-HOC")


# --------------------------------------------------------------------------
# restriction, which is the direction that works
# --------------------------------------------------------------------------

def test_restriction_is_complete_at_every_pair(art):
    pairs = art["restriction"]["pairs"]
    assert len(pairs) == 6
    for pair in pairs:
        assert pair["restriction_is_complete"], pair["lower_rows_not_implied"]


def test_restriction_recovers_the_borland_dennis_equalities():
    """The (3,7) system implies all three (3,6) equalities, from scratch.

    Trace three plus rows 1, 2 and 4 of (3,7) give `l_1 + l_6 >= 1`,
    `l_2 + l_5 >= 1` and `l_3 + l_4 >= 1`, whose sum is exactly three, so each
    is an equality. This is the whole restriction argument at its smallest.
    """
    rows = [([int(c) for c in r["coeffs"]], int(r["rhs"]))
            for r in constraints(3, 7)["inequalities"]]
    restricted = [(c[:6], b) for c, b in rows]
    pairs = {(0, 5), (1, 4), (2, 3)}
    found = set()
    for coeffs, rhs in restricted:
        zeros = tuple(i for i, c in enumerate(coeffs) if c == 0)
        if len(zeros) == 2 and all(c == 1 for c in coeffs if c) and rhs == 2:
            found.add(zeros)
    assert found == pairs
    for equality in constraints(3, 6)["equalities"]:
        coeffs = [int(c) for c in equality["coeffs"]]
        assert tuple(i for i, c in enumerate(coeffs) if c) in pairs
        assert int(equality["rhs"]) == 1


# --------------------------------------------------------------------------
# templates and the projection
# --------------------------------------------------------------------------

def test_template_stock_is_small_but_does_not_saturate(art):
    per_system = {r["system"]: r for r in art["templates"]["per_system"]}
    assert per_system["(3,7)"]["distinct_symmetric_templates"] == 1
    assert per_system["(3,10)"]["distinct_symmetric_templates"] == 16
    n3 = ["(3,6)", "(3,7)", "(3,8)", "(3,9)", "(3,10)"]
    assert [per_system[s]["new_symmetric_templates_vs_previous_rank"]
            for s in n3] == [1, 1, 8, 6, 9]
    assert [per_system[s]["max_symmetric_arity"] for s in n3] == [3, 4, 7, 7, 8]


def test_the_rank_40_projection_is_out_of_reach(art):
    totals = art["cost_projection"]["total_candidate_rows"]
    assert totals["20"] > 10 ** 9
    assert totals["40"] > 10 ** 12
    assert art["cost_projection"]["status"].startswith("PROJECTION")


# --------------------------------------------------------------------------
# the dependency graph
# --------------------------------------------------------------------------

def test_dependency_graph_is_well_formed(graph):
    allowed = set(graph["status_legend"])
    node_ids = {n["id"] for n in graph["nodes"]}
    source_ids = set(graph["sources"])
    arrow_ids = {a["id"] for a in graph["arrows"]}
    assert len(arrow_ids) == len(graph["arrows"])
    for arrow in graph["arrows"]:
        assert arrow["status"] in allowed, arrow
        for end in (arrow["from"], arrow["to"]):
            assert end in node_ids | source_ids | arrow_ids, (arrow["id"], end)
        for source in arrow.get("sources", []):
            assert source in source_ids, (arrow["id"], source)


def test_every_scope_audit_arrow_is_refuted(graph):
    """The audit's whole content is that none of the cited theorems applies."""
    audits = [a for a in graph["arrows"] if a.get("kind", "").startswith("scope audit")]
    assert len(audits) == 6
    for arrow in audits:
        assert arrow["status"] == "FALSE"
        assert arrow["obstruction"]


def test_the_missing_theorem_is_the_only_missing_arrow(graph):
    missing = [a for a in graph["arrows"] if a["status"] == "MISSING THEOREM"]
    assert [a["id"] for a in missing] == ["A9"]


def test_graph_verdicts_track_the_measurement(art, graph):
    by_id = {a["id"]: a for a in graph["arrows"]}
    assert by_id["A7"]["check"]["verdict"] == art["scorecard"][
        "P5_restriction_is_complete"]["verdict"]
    assert by_id["A8"]["counterexample"]["verdict"] == art["scorecard"][
        "P4_zero_padding_preserves_validity_above_rank_6"]["verdict"]
    assert by_id["A4rev"]["counterexample"]["witnesses"] == [
        w["vertex"] for w in art["saturation"]["non_saturation_witnesses"]]
