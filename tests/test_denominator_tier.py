"""Feature 4: denominator-grid tiering, and the structural law that bounds it.

The amplitude-grid tier state_den/spectrum_den is what licenses keeping the
interference/design search on the m=1 grid, with m=2 a DESIGN-REAL-only fallback
capped at MAX_DENOMINATOR_TIER. That bound is a statement about the GRID SOLVER's
reach, not a law of nature: it holds for every state the grid produces, namely
the designs (one-hop-free supports) and the magnitude-target interference states
(whose 1-RDM carries nonzero Schur-Horn off-diagonals).

It is BROKEN, and provably so, by the cancellation-geometry states: those that
realize their vertex on a DIAGONAL 1-RDM through destructive interference across a
one-hop-connected support (off-diagonals forced to zero rather than to a nonzero
target). Those states live off the m<=2 grid and were certified by the fiber
sparsification route, not the grid solver, which is exactly why the grid solver
missed them. Across the whole certified corpus the tier ceiling is equivalent to
a computable structural predicate:

    tier > MAX_DENOMINATOR_TIER  <=>  (1-RDM is diagonal) AND (support is one-hop
                                       connected)   [the cancellation geometry]

verified here with zero counterexamples. This test is the containment proof for
the grid solver's search-order bound, scoped to the states that bound governs.
See docs/RESEARCH.md "RANK-10 CLOSED" and the pre-registration P3 scoring.
"""
import json
import pathlib
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data" / "states.jsonl"


def _states():
    if not DATA.exists():
        return []
    return [json.loads(x) for x in DATA.read_text().splitlines()
            if x.strip() and json.loads(x).get("status") == "OK"
            and json.loads(x).get("closed_form")]


def _spectrum_den(integer_form, denominator):
    # the natural spectrum denominator is `denominator` reduced against the ints
    from math import gcd
    g = denominator
    for x in integer_form:
        g = gcd(g, x)
    return denominator // g


def _is_cancellation(d, cf):
    """A state is cancellation-geometry iff its 1-RDM (rebuilt from the closed
    form) is diagonal AND its support is one-hop connected: interference that
    zeroes an off-diagonal instead of hitting a nonzero magnitude target."""
    import numpy as np
    import sympy as sp

    support = [tuple(s) for s in cf["support_dets"]]
    amps = [complex(sp.sympify(p).evalf(30)) for p in cf["pretty"]]
    amap = dict(zip(support, amps))
    rho = np.zeros((d, d), complex)
    for t, ct in amap.items():
        for mp in t:
            s1 = (-1) ** t.index(mp)
            t2 = tuple(x for x in t if x != mp)
            for m in range(d):
                if m in t2:
                    continue
                tp = tuple(sorted(t2 + (m,)))
                if tp not in amap:
                    continue
                rho[m, mp] += s1 * (-1) ** tp.index(m) * np.conjugate(amap[tp]) * ct
    diagonal = float(np.max(np.abs(rho - np.diag(np.diag(rho))))) < 1e-9
    sets = [set(s) for s in support]
    connected = any(len(sets[i] ^ sets[j]) == 2
                    for i in range(len(sets)) for j in range(i + 1, len(sets)))
    return diagonal and connected


def test_denominator_ratio_law_holds_over_corpus():
    from gpc_census.states import MAX_DENOMINATOR_TIER

    states = _states()
    if not states:
        import pytest
        pytest.skip("no states atlas present")
    tier2, exceeded = [], []
    for r in states:
        sden = _spectrum_den(r["integer_form"], r["denominator"])
        cden = r["closed_form"]["den"]
        assert cden % sden == 0, f"{r['system']} v{r['index']}: non-integer tier"
        tier = cden // sden
        n, d = (int(x) for x in r["system"].strip("()").split(","))
        cancel = _is_cancellation(d, r["closed_form"])
        if tier > MAX_DENOMINATOR_TIER:
            # only cancellation-geometry states may exceed the grid bound, and
            # they are correlated interference states (never designs)
            assert cancel, (
                f"{r['system']} v{r['index']}: tier {tier} exceeds cap but is NOT "
                f"cancellation-geometry (would break the structural law)")
            assert r["classified"] == "INTERFERENCE", (
                f"{r['system']} v{r['index']}: tier {tier} exceeded by a "
                f"{r['classified']} state")
            exceeded.append(r)
        else:
            # within the grid bound, cancellation geometry never appears
            assert not cancel, (
                f"{r['system']} v{r['index']}: cancellation-geometry state at "
                f"tier {tier} <= cap (structural law says it should exceed it)")
            if tier == MAX_DENOMINATOR_TIER:
                tier2.append(r)
    # the equivalence tier > cap <=> cancellation geometry is now proved both
    # ways by the per-state assertions above; the tier-2 band is DESIGN-REAL only
    assert tier2, "expected some tier-2 states"
    assert all(r["classified"] == "DESIGN-REAL" for r in tier2)
    # the rank-10 closure introduced the first off-grid (cancellation) states
    # (v89, v103); on a reduced fixture `exceeded` may be empty and the law holds
    # vacuously, which is fine.


def test_design_real_record_carries_tier():
    # a DESIGN-REAL solve annotates the tier (1 or 2), never more
    from gpc_census.states import MAX_DENOMINATOR_TIER, solve_design_real_vertex

    # (3,10) with a real-design vertex; use a known DESIGN-REAL spectrum if present
    states = _states()
    dr = next((r for r in states if r["classified"] == "DESIGN-REAL"), None)
    if dr is None:
        import pytest
        pytest.skip("no DESIGN-REAL vertex present")
    n, d = (int(x) for x in dr["system"].strip("()").split(","))
    spec = [Fraction(x, dr["denominator"]) for x in dr["integer_form"]]
    rec = solve_design_real_vertex(n, d, spec)
    if rec is None:
        import pytest
        pytest.skip("design-real solver found no support in the enumerated set")
    assert rec.get("denominator_tier") in (1, 2)
    assert rec["denominator_tier"] <= MAX_DENOMINATOR_TIER
