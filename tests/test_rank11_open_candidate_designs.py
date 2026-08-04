import importlib.util
import itertools
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "rank11_open_candidate_designs.json"

MODULE_PATH = ROOT / "scripts" / "rank11_open_candidate_designs.py"
SPEC = importlib.util.spec_from_file_location("rank11_designs", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def artifact():
    return json.loads(ARTIFACT.read_text())


def test_the_four_open_candidates_are_the_ones_the_settlement_records():
    settlement = json.loads((ROOT / "docs" / "bracket_3_11_settlement.json").read_text())
    still_open = sorted(
        c["index"] for c in settlement["candidates"] if c.get("verdict") == "OPEN"
    )
    assert still_open == [23, 26, 34, 44]
    assert sorted(int(k) for k in artifact()["candidates"]) == still_open


def test_every_candidate_has_trace_three():
    for record in artifact()["candidates"].values():
        assert sum(F(x) for x in record["spectrum"]) == 3


def test_all_four_satisfy_every_proved_family():
    """The open verdicts rest on this: no PROVED row refutes any of them."""
    for record in artifact()["candidates"].values():
        assert record["families"]["satisfies_every_proved_family"]


def test_candidate_44_rests_on_a_single_claimed_row():
    """44 is the cheapest target because only one claimed row stands against it.

    The others carry two or three, so attaining 44 would falsify exactly one
    published claim rather than a cluster.
    """
    candidates = artifact()["candidates"]
    assert candidates["44"]["families"]["violated_claimed_count"] == 1
    only = candidates["44"]["families"]["violated"][0]
    assert only["family"].startswith("KLY09-EXT-1")
    assert only["tier"] == "CLAIMED"
    for index in ("23", "26", "34"):
        assert candidates[index]["families"]["violated_claimed_count"] >= 2


def test_each_candidate_misses_its_claimed_row_by_one_over_its_denominator():
    denominators = {"23": 7, "26": 27, "34": 17, "44": 12}
    for index, record in artifact()["candidates"].items():
        excesses = {F(v["excess"]) for v in record["families"]["violated"]}
        assert excesses == {F(1, denominators[index])}


def test_no_candidate_admits_a_one_hop_free_design():
    """Extends the settlement's note, which recorded this only for cand 23.

    The solver verdicts are floating point and are evidence; the exact
    certificate is candidate 44's counting proof, checked separately below.
    """
    for record in artifact()["candidates"].values():
        assert record["design_search"]["one_hop_free_design_found"] is False
        assert record["design_search"]["evidence_only"] is True


def test_one_hop_free_is_exactly_the_partial_steiner_condition():
    blocks, conflicts = mod.blocks_and_conflicts()
    assert len(blocks) == 165
    for i, j in conflicts:
        assert len(set(blocks[i]) & set(blocks[j])) == 2
    conflict_set = set(conflicts)
    for i, j in itertools.combinations(range(len(blocks)), 2):
        shares_two = len(set(blocks[i]) & set(blocks[j])) == 2
        assert shares_two == ((i, j) in conflict_set)


def test_candidate_44_exact_no_design_proof():
    """The one certificate here, independent of the floating-point solver."""
    proof = artifact()["candidates"]["44"]["exact_no_design_proof"]
    assert proof["applies"] and proof["exact"]
    assert len(proof["heavy_orbitals"]) == 5
    assert len(proof["light_orbitals"]) == 6
    assert F(proof["forced_inside_weight_at_least"]) == F(1, 2)
    assert proof["max_blocks_inside_a_five_set"] == 2
    assert proof["case_two_blocks"]["contradiction"]
    assert proof["case_two_blocks"]["free_heavy_pairs"] == 4
    assert proof["case_two_blocks"]["light_orbitals_needing_positive_weight"] == 6
    assert proof["case_one_block"]["contradiction"]
    assert proof["case_zero_blocks"]["contradiction"]


def test_the_forced_weight_identity_behind_the_proof():
    """a - c - 2d = 1/2 for candidate 44, recomputed from the spectrum alone."""
    lam = [F(x) for x in artifact()["candidates"]["44"]["spectrum"]]
    heavy = sum(x for x in lam if x == F(1, 2))
    light = sum(x for x in lam if x == F(1, 12))
    assert heavy == F(5, 2) and light == F(1, 2) and heavy + light == 3
    # 3a+2b+c = heavy, a+b+c+d = 1, so a - c - 2d = heavy - 2
    assert heavy - 2 == F(1, 2)


def test_at_most_two_blocks_fit_inside_a_five_set():
    inside = list(itertools.combinations(range(5), 3))
    best = 0
    for k in range(len(inside) + 1):
        for family in itertools.combinations(inside, k):
            if all(len(set(a) & set(b)) <= 1 for a, b in itertools.combinations(family, 2)):
                best = max(best, k)
    assert best == 2


def test_artifact_claims_no_verdict():
    """This layer narrows the four; it does not settle any of them."""
    data = artifact()
    assert any("settled" in n for n in data["nonclaims"])
    assert any("evidence" in n for n in data["nonclaims"])
    for record in data["candidates"].values():
        assert "verdict" not in record
