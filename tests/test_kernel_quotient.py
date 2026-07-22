"""Feature 2: kernel quotienting of the weight search.

On a fixed support the integer weight vectors solving the degree system form a
coset of the incidence kernel; free-kernel translates are phase-equivalent, so
only one representative per (support, rigid-class products) is phase-solved.
Sound because rigidly-pinned distinctions keep separate keys and verify_exact
gates every hit. This test checks the kernel-dim primitive against the shipped
fiber-dimension census and that quotienting keeps a certified vertex certifying.
"""
import json
import pathlib
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data" / "states.jsonl"


def test_incidence_kernel_dim_matches_corpus():
    from gpc_census.states import _incidence_kernel_dim

    if not DATA.exists():
        import pytest
        pytest.skip("no states atlas present")
    rows = [json.loads(x) for x in DATA.read_text().splitlines() if x.strip()]
    # v_B has fiber (kernel) dimension 1; a Slater corner (single det) has 0
    vb = next(r for r in rows if r["system"] == "(4,9)" and r["index"] == 65)
    assert _incidence_kernel_dim([tuple(t) for t in vb["closed_form"]["support_dets"]]) == 1
    slater = next(r for r in rows if r["system"] == "(3,6)" and r["index"] == 0)
    assert _incidence_kernel_dim([tuple(t) for t in slater["closed_form"]["support_dets"]]) == 0


def test_quotient_keeps_certification_and_records_fiber_dim():
    import gpc_census.states as S

    # a small certified loopy interference vertex, fast to solve
    spec = [Fraction(x, 7) for x in (4, 4, 4, 4, 2, 1, 1, 1)]  # (3,8) v29
    assert S._KERNEL_QUOTIENT is True
    rec = S.solve_vertex_exact_first(3, 8, spec, max_card=20, max_clique=3,
                                     max_cliques=1, clique_time_budget=20,
                                     certify_tier_b=True)
    assert rec and rec.get("status") == "OK"
    assert rec.get("exact", {}).get("status") == "EXACT"
    assert "fiber_kernel_dim" in rec and rec["fiber_kernel_dim"] >= 0


def test_quotient_toggle_both_certify():
    import gpc_census.states as S

    spec = [Fraction(x, 7) for x in (4, 4, 4, 4, 2, 1, 1, 1)]
    try:
        for q in (True, False):
            S._KERNEL_QUOTIENT = q
            rec = S.solve_vertex_exact_first(3, 8, spec, max_card=20, max_clique=3,
                                             max_cliques=1, clique_time_budget=20,
                                             certify_tier_b=True)
            assert rec and rec.get("exact", {}).get("status") == "EXACT", f"quotient={q}"
    finally:
        S._KERNEL_QUOTIENT = True
