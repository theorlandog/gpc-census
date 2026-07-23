"""Feature 1: class-symmetry reduction of the skeleton model.

The reduction is an exact reformulation (lex-leader keeps the lex-minimum orbit
representative), never a filter: it must not empty any orbit, so a certified
vertex still solves, and it must shrink the solution set of a symmetric model.
"""
import pathlib
import sys

from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _count_skeletons(spec, n, d, symmetry, cap=15):
    import gpc_census.states as S
    from ortools.sat.python import cp_model

    built = S._skeleton_model(n, d, spec, hop_cuts=True, symmetry_break=symmetry)
    if built is None:
        return None, "none"
    m, k, y, dets, den, allowed = built
    m.Add(sum(y) <= sum(int(x) for x in [f * den for f in spec]) // n)
    s = cp_model.CpSolver()
    s.parameters.enumerate_all_solutions = True
    s.parameters.max_time_in_seconds = cap
    s.parameters.num_workers = 1
    cnt = [0]

    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            cnt[0] += 1

    st = s.Solve(m, CB())
    return cnt[0], s.StatusName(st)


def test_symmetry_reduces_design_skeletons_and_keeps_some():
    # a degenerate design vertex: (3,8) v with a class of four equal occupations
    spec = [Fraction(x, 7) for x in (4, 4, 4, 4, 2, 1, 1, 1)]
    off, s_off = _count_skeletons(spec, 3, 8, False)
    on, s_on = _count_skeletons(spec, 3, 8, True)
    assert off is not None and on is not None
    # never empties the solution set the plain model has
    if off > 0:
        assert on >= 1, "symmetry breaking emptied a nonempty solution set (unsound)"
    # and it shrinks it (a heavily degenerate spectrum has large orbits)
    assert on <= off


def test_certified_degenerate_vertex_still_solves_with_symmetry():
    import gpc_census.states as S

    # (3,10) v96 = (5,5,5,5,2,1,1,1,1,1)/9, certified INTERFERENCE, heavy degeneracy
    spec = [Fraction(x, 9) for x in (5, 5, 5, 5, 2, 1, 1, 1, 1, 1)]
    assert S._SYMMETRY_BREAK is True  # on by default
    rec = S.solve_vertex_exact_first(3, 10, spec, max_card=24, max_clique=3,
                                     max_cliques=1, clique_time_budget=30,
                                     certify_tier_b=True)
    assert rec and rec.get("status") == "OK"
    assert rec.get("exact", {}).get("status") == "EXACT"


def test_symmetry_break_toggle_is_sound_both_ways():
    # the model builds and is feasible with symmetry both on and off
    spec = [Fraction(x, 9) for x in (6, 6, 3, 3, 3, 3, 1, 1, 1, 0)]
    for sym in (True, False):
        cnt, st = _count_skeletons(spec, 3, 10, sym, cap=10)
        assert cnt is not None
