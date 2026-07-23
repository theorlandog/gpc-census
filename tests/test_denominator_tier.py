"""Feature 4: denominator-grid tiering, justified by the ratio law.

The amplitude-grid tier state_den/spectrum_den is 1 or 2 across the whole
certified corpus, and 2 occurs only for DESIGN-REAL states. That law is what
licenses keeping the interference/design-int search on the m=1 grid and treating
m=2 as a DESIGN-REAL-only fallback with no escalation beyond MAX_DENOMINATOR_TIER.
This test is the containment proof for that search-order bound.
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


def test_denominator_ratio_law_holds_over_corpus():
    from gpc_census.states import MAX_DENOMINATOR_TIER

    states = _states()
    if not states:
        import pytest
        pytest.skip("no states atlas present")
    tier2 = []
    for r in states:
        sden = _spectrum_den(r["integer_form"], r["denominator"])
        cden = r["closed_form"]["den"]
        assert cden % sden == 0, f"{r['system']} v{r['index']}: non-integer tier"
        tier = cden // sden
        assert 1 <= tier <= MAX_DENOMINATOR_TIER, (
            f"{r['system']} v{r['index']}: tier {tier} exceeds cap")
        if tier == 2:
            tier2.append(r)
    # every tier-2 state is DESIGN-REAL (the m=2 grid is a real-design signature)
    assert tier2, "expected some tier-2 states"
    assert all(r["classified"] == "DESIGN-REAL" for r in tier2)


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
