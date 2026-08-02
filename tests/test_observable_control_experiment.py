import importlib.util
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "observable_control_experiment.json"

MODULE_PATH = ROOT / "scripts" / "observable_control_experiment.py"
SPEC = importlib.util.spec_from_file_location("observable_control", MODULE_PATH)
oce = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(oce)

CENSUS_MODULE = ROOT / "scripts" / "census_support_geometry.py"
CSPEC = importlib.util.spec_from_file_location("census_support_geometry", CENSUS_MODULE)
csg = importlib.util.module_from_spec(CSPEC)
CSPEC.loader.exec_module(csg)


def test_disjoint_determinants_need_three_body_control_and_readout():
    support = [(0, 1, 2), (3, 4, 5)]
    assert oce.control_q(support, 3) == 3
    assert oce.observe_q(support, 3) == 3


def test_observability_never_precedes_control_examples():
    examples = [
        [(0, 1, 2), (0, 1, 3), (0, 4, 5)],
        [(0, 1, 2), (0, 3, 4), (1, 3, 5), (2, 4, 5)],
    ]
    for support in examples:
        assert oce.observe_q(support, 3) >= oce.control_q(support, 3)


def test_the_inequality_is_a_subgraph_theorem_not_a_measurement():
    """q_obs >= q_ctrl holds because one graph contains the other.

    A q-RDM channel for the pair (a,b) requires dist(D_a,D_b) <= q, so the
    collision-free observable graph at order q is a subgraph of the control
    graph at the same order. Check the containment itself, since that is the
    reason, rather than only its consequence.
    """
    universe = oce.dets(6, 3)
    for m in (2, 3, 4):
        for support in itertools.combinations(universe, m):
            support = list(support)
            for q in (1, 2, 3):
                observable, _ = oce.unique_edges(support, 3, q)
                control = {
                    (i, j)
                    for i, j in itertools.combinations(range(m), 2)
                    if oce.dist(support[i], support[j]) <= q
                }
                assert observable <= control


def test_sign_free_channel_model_matches_the_exact_signed_one():
    """The experiment drops fermionic signs; that has to be justified.

    Each carrier pair contributes exactly one plus-or-minus-one term to a given
    channel, so nothing can cancel and the sign-free edge set coincides with
    the exactly signed one. Verified against the census campaign's signed
    implementation.
    """
    universe = oce.dets(6, 3)
    for m in (2, 3, 4):
        for support in itertools.combinations(universe, m):
            support = list(support)
            row = {
                "system": "(3,6)",
                "index": 0,
                "classified": "X",
                "closed_form": {"support_dets": [list(x) for x in support]},
            }
            for q in (2, 3):
                channels = csg.channel_terms(support, 6, 3, q)
                signed = {(t[0][0], t[0][1]) for t in channels.values() if len(t) == 1}
                unsigned, _ = oce.unique_edges(support, 3, q)
                assert signed == unsigned, (support, q)
            exact = csg.analyze_record(row)["minimum_collision_free_reconstruction_order"]
            assert max(oce.observe_q(support, 3), 2) == exact


def test_artifact_summary_is_internally_consistent():
    summary = json.loads(ARTIFACT.read_text())["summary"]
    assert summary["q_observe_lt_control"] == 0
    assert summary["q_equal"] + summary["q_observe_gt_control"] == summary["supports"]
    assert sum(summary["pair_histogram"].values()) == summary["supports"]
    assert sum(summary["gap_histogram"].values()) == summary["supports"]
    for key in summary["pair_histogram"]:
        control, observe = (int(x) for x in key.split(","))
        assert observe >= control


def test_recorded_gaps_are_never_more_than_one_order():
    summary = json.loads(ARTIFACT.read_text())["summary"]
    assert set(summary["gap_histogram"]) <= {"0", "1"}
