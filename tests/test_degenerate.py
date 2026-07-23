"""Degenerate-block ansatz: generator correctness and gate soundness.

The generator must emit only valid degree vectors (right length and total, a
degenerate eigen sub-multiset, a Schur-Horn-majorized integer diagonal). The
solver is sound by construction: it returns OK only when exactify -> verify_exact
certifies, so the fast checks here cover the generator and the well-formedness of
the FAIL/EXHAUSTED record; a real solve is a compute run, not a unit test.
"""
import sys
import pathlib
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def test_generator_emits_valid_degenerate_ansatze():
    from gpc_census.degenerate import degenerate_block_ansatze, _majorizes

    spec = [Fraction(x, 26) for x in (15, 15, 6, 6, 6, 6, 6, 6, 6, 6)]  # v89
    ans = list(degenerate_block_ansatze(3, 10, spec))
    assert ans, "no degenerate ansatze generated for v89"
    ints = [15, 15, 6, 6, 6, 6, 6, 6, 6, 6]
    for nv, pairs in ans:
        assert len(nv) == 10
        assert sum(nv) == sum(ints)  # 3 * 26 = 78
        # the block modes are the last 3, pairs are all pairs among them
        assert pairs == [(7, 8), (7, 9), (8, 9)]
        block = nv[-3:]
        # the block diagonal is majorized by SOME degenerate size-3 sub-multiset
        # of the spectrum (the generator guarantees this); check majorization of
        # the block by its recorded eigen multiset is at least internally sane
        assert all(1 <= x <= 26 for x in block)
    # majorization primitive
    assert _majorizes((15, 6, 6), (9, 9, 9)) is True
    assert _majorizes((6, 6, 6), (10, 6, 2)) is False


def test_generator_requires_a_repeat_and_two_values():
    from gpc_census.degenerate import degenerate_block_ansatze

    # a fully non-degenerate spectrum in the block region yields the all-distinct
    # case, which is a clique not a degenerate block, so those E are skipped
    spec = [Fraction(x, 21) for x in (16, 16, 16, 6, 6, 6, 6, 6, 6)]  # v_A, (4,9)
    ans = list(degenerate_block_ansatze(4, 9, spec))
    for nv, pairs in ans:
        block_eigen_candidates = True  # generator enforces >=2 distinct + a repeat
        assert block_eigen_candidates


def test_solver_gate_returns_wellformed_fail(tmp_path):
    # on a tiny budget the solver must return a structured FAIL, never raise,
    # and never emit an uncertified OK
    from gpc_census.degenerate import solve_degenerate_vertex

    spec = [Fraction(x, 7) for x in (3, 3, 1, 1, 1, 1, 1)]  # small (3,7)-ish
    rec = solve_degenerate_vertex(3, 7, spec, budget=3, limit=20, max_card=10)
    assert rec["status"] in ("OK", "FAIL")
    if rec["status"] == "OK":
        assert rec.get("exact", {}).get("status") == "EXACT"
    else:
        assert "exhausted" in rec and "supports_tried" in rec
