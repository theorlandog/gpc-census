import importlib.util
import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "hybrid_rank_generator_prototype.json"
STATES = DATA / "states.jsonl"

MODULE_PATH = ROOT / "scripts" / "hybrid_rank_generator_prototype.py"
SPEC = importlib.util.spec_from_file_location("hybrid_rank", MODULE_PATH)
hy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hy)


def test_bd_gpc_recovered():
    facets = hy.recover_facets(hy.BD_VERTICES)
    assert any(f["a"] == [1, -1, -1] and f["b"] == 0 for f in facets)


def test_prototype_vertices_are_the_certified_oracle():
    """A regression harness is only meaningful against the real oracle.

    Its whole purpose is reproducing the repository's certified vertices, so
    check the input, not only the recovered output.
    """
    rows = [json.loads(x) for x in STATES.read_text().splitlines() if x.strip()]
    oracle = sorted(
        tuple(Fraction(v, r["denominator"]) for v in r["integer_form"])
        for r in rows
        if r["system"] == "(3,6)"
    )
    lifted = sorted(tuple(Fraction(str(x)) for x in hy.lift_lambda(v)) for v in hy.BD_VERTICES)
    assert len(oracle) == 4
    assert lifted == oracle


def test_recovered_facets_cut_out_the_vertices():
    facets = hy.recover_facets(hy.BD_VERTICES)
    assert len(facets) == 4
    for vertex in hy.BD_VERTICES:
        tight = 0
        for f in facets:
            value = sum(f["a"][j] * vertex[j] for j in range(3))
            assert value <= f["b"]
            tight += value == f["b"]
        assert tight >= 3


def test_every_recovered_facet_is_labelled():
    """A fall-through label means the naming table drifted from the normals."""
    for f in hy.recover_facets(hy.BD_VERTICES):
        meaning = hy.facet_meaning(f["a"], f["b"])
        assert not meaning.startswith("("), (f["a"], f["b"], meaning)


def test_pinned_carrier_gamma2():
    S = hy.bd_selection_carrier()
    assert len(S) == 9
    assert all(len(set(D) & {0, 1, 3}) == 2 for D in S)
    assert hy.incidence_rank(S, 2) == 9
    assert hy.collision_free_graph(S, 2)["connected"]


def test_pinned_carrier_matches_the_artifact():
    carrier = json.loads(ARTIFACT.read_text())["pinned_carrier"]
    S = hy.bd_selection_carrier()
    assert carrier["size"] == len(S) == 9
    assert carrier["singleton_rank"] == hy.incidence_rank(S, 1)
    assert carrier["pair_incidence_rank"] == hy.incidence_rank(S, 2)
    assert carrier["non_toric_phase_complexity"] == len(S) - carrier["singleton_rank"]
    assert carrier["gamma2_graph"] == hy.collision_free_graph(S, 2)
    assert carrier["pure_full_support_gamma2_certified"] is True
