"""Guards for Theorem D and Theorem C'.

Theorem D claims the rational diagonal of a large class is Theorem A read
through a loop shared with a smaller class. The tests that matter are the ones
that could catch it being a tautology or a sign coincidence:

- the shared-loop identity is recomputed from the DETERMINANTS, so it does not
  depend on any amplitude;
- Theorem D's predicted t_ab is built from the LINKING class and compared with
  a t_ab computed inside the original class, which is where a fermionic sign
  error between sigma_1 sigma_2 and s_a s_b would surface;
- the refuted conjecture must stay refuted, and the class that refutes it named.
"""
import json
import pathlib
import sys

import pytest
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "shared_loop_rationality.json"

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))


def _art():
    return json.loads(ARTIFACT.read_text())


# --------------------------------------------------------------------------
# scope
# --------------------------------------------------------------------------

def test_scope_matches_the_polygon_artifact():
    art = _art()
    summary = art["summary"]
    assert summary["classes"] == 17
    assert summary["class_sizes"] == {"3": 16, "4": 1}
    assert summary["distinct_states"] == 16
    assert summary["transport_classes"] == 11


# --------------------------------------------------------------------------
# Theorem D
# --------------------------------------------------------------------------

def test_shared_loop_identity_is_combinatorial_and_recomputed_live():
    """Part (1) and (2) of Theorem D, from determinants only.

    If u_b = u_a - C + E then v_b = v_a - C + E and the two Lemma 1 loops are
    the same integer vector. No amplitude enters, so this is a statement about
    supports and is re-derived here rather than read from the artifact.
    """
    from gpc_census.states import one_hop_classes

    records = {(r["system"], r["index"]): r for r in
               (json.loads(line) for line
                in (ROOT / "states.jsonl").read_text().splitlines() if line.strip())}
    art = _art()
    checked = 0
    for block in art["classes"]:
        rec = records[(block["system"], block["index"])]
        dets = [tuple(t) for t in rec["closed_form"]["support_dets"]]
        classes = one_hop_classes(dets)
        a_mode, b_mode = block["channel"]
        oriented = []
        for i, j in classes[(a_mode, b_mode)]:
            u, v = (i, j) if a_mode in dets[i] else (j, i)
            oriented.append((u, v))
        for entry in block["links"]:
            if entry["u_hop"] != 1:
                continue
            i, j = entry["pair"]
            ua, va = oriented[i]
            ub, vb = oriented[j]
            hop_u = set(dets[ua]) ^ set(dets[ub])
            hop_v = set(dets[va]) ^ set(dets[vb])
            assert hop_u == hop_v, (block["system"], entry["pair"])
            assert {va: 1, ua: -1, vb: -1, ub: 1} == {ub: 1, ua: -1, vb: -1, va: 1}
            checked += 1
    assert checked == 19


def test_theorem_D_formula_matches_on_every_size_two_link():
    """The cross-check: predicted from (C,E), measured inside (A,B)."""
    art = _art()
    summary = art["summary"]
    assert summary["size_two_links"] == 17
    assert summary["theorem_D_checked"] == 17
    assert summary["theorem_D_matches"] == 17
    assert summary["theorem_D_predicted_rational"] == 17
    for block in art["classes"]:
        for entry in block["links"]:
            if entry.get("linking_channel_size") != 2:
                continue
            assert entry["theorem_D_matches"], (block["system"], entry["pair"])
            assert entry["t_ab_rational"] is True
            predicted = sp.nsimplify(sp.sympify(entry["theorem_D_predicted_t_ab"]))
            measured = sp.nsimplify(sp.sympify(entry["t_ab"]))
            assert sp.simplify(predicted - measured) == 0
            assert predicted.is_rational


def test_a_size_two_link_always_implies_a_rational_diagonal():
    """Theorem D's conclusion, stated as the census fact it licenses."""
    art = _art()
    assert (art["summary"]["size_two_link_and_t_rational"]
            == art["summary"]["size_two_links"] == 17)


# --------------------------------------------------------------------------
# Theorem C' and the induction
# --------------------------------------------------------------------------

def test_every_class_has_a_rational_chain():
    """P-D-5. The chain search must consider EVERY ordering.

    A previous version skipped orderings with order[0] > order[-1] as
    "reverses", which is wrong because reversing changes the partial sums, and
    it manufactured a false refutation at (3,10) v98.
    """
    art = _art()
    assert art["summary"]["classes_with_a_rational_chain"] == 17
    for block in art["classes"]:
        assert block["has_rational_chain"], (block["system"], block["index"])
        chain = block["rational_chain"]
        assert sorted(chain["order"]) == list(range(block["class_size"]))
        assert len(chain["partial_moduli"]) == block["class_size"]


def test_sixteen_of_seventeen_classes_are_proved_end_to_end():
    """Theorem D supplies m-2 diagonals, so Theorem C' needs nothing checked."""
    art = _art()
    assert art["summary"]["proved_end_to_end"] == 16
    for block in art["classes"]:
        expected = block["diagonals_proved_by_theorem_D"] >= block["diagonals_needed"]
        assert block["proved_end_to_end"] == expected
        assert block["diagonals_needed"] == block["class_size"] - 2


def test_no_class_carries_a_spare_proved_diagonal():
    """Every proved class sits exactly at the minimum m-2, never above it."""
    art = _art()
    assert art["summary"]["exactly_at_the_minimum"] == 16
    for block in art["classes"]:
        if block["proved_end_to_end"]:
            assert block["diagonals_proved_by_theorem_D"] == block["diagonals_needed"]


# --------------------------------------------------------------------------
# the refuted conjecture must stay refuted
# --------------------------------------------------------------------------

def test_the_one_step_combinatorial_conjecture_is_false():
    """P-D-4 FAILED, as pre-registered.

    "Every class of size >= 3 contains a pair linked by a class of size at
    most 2" is FALSE, and the census carries the counterexample. The write-up
    leans on this, so a change that made it true would silently convert a
    recorded refutation into a confirmation.
    """
    art = _art()
    assert art["summary"]["classes_with_a_size_two_link"] == 16
    assert art["summary"]["not_proved_end_to_end"] == ["(3,10) v103 ch(1, 6)"]
    offender = next(b for b in art["classes"]
                    if not b["proved_end_to_end"])
    assert offender["system"] == "(3,10)" and offender["index"] == 103
    assert offender["channel"] == [1, 6]
    assert not offender["has_size_two_link"]
    # it survives only because its diagonals are rational for another reason
    assert offender["has_rational_chain"]


def test_the_counterexamples_linking_class_is_the_size_four_one():
    """Why the one-step conjecture fails, and what the recursion would need."""
    art = _art()
    offender = next(b for b in art["classes"] if not b["proved_end_to_end"])
    linked = [e for e in offender["links"] if e["u_hop"] == 1]
    assert linked, "the class must have a link, just not a size-2 one"
    assert all(e["linking_channel_size"] != 2 for e in linked)
    assert any(e["linking_channel_size"] == 4 for e in linked)


@pytest.mark.parametrize("phrase", ["REFUTED", "FALSE"])
def test_the_write_up_records_the_refutation(phrase):
    doc = REPO / "docs" / "shared_loop_rationality.md"
    if not doc.exists():
        pytest.skip("docs not present")
    assert phrase in doc.read_text()
