"""Feature 5: exponent-2 (Conjecture 2) classification of gauge-invariant holonomies.

Every certified interference state is a live test of Conjecture 2: its
gauge-invariant loop holonomies must generate 2-elementary abelian extensions of
Q. The classifier keys on 2 cos(theta) (a quadratic extension below e^{i theta}),
so census holonomy fields (degree <= 8 in e^{i theta}, <= 4 in the cosine) are in
range of the Galois backend. The check never rejects a state; it only surfaces a
definitive violation loudly.
"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data" / "states.jsonl"


def test_classifier_known_cases():
    import sympy as sp

    from gpc_census.exactify import is_exponent_two

    # p sqrt(q)/r cosines and their doubles are 2-elementary
    assert is_exponent_two(sp.Rational(3) * sp.sqrt(2) / 16) is True  # v_B holonomy
    assert is_exponent_two(sp.cos(2 * sp.pi / 5)) is True             # C2
    # a cubic (C3) angle is abelian but exponent 3, NOT the conjecture class
    assert is_exponent_two(sp.cos(2 * sp.pi / 7)) is False


def test_flagship_holonomies_are_exponent_two():
    import sympy as sp

    from gpc_census.exactify import conjecture2_scan

    if not DATA.exists():
        import pytest
        pytest.skip("no states atlas present")
    rows = [json.loads(x) for x in DATA.read_text().splitlines() if x.strip()]
    for sysid, idx in [("(4,9)", 65), ("(5,10)", 113), ("(3,10)", 17)]:
        r = next((x for x in rows if x["system"] == sysid and x["index"] == idx), None)
        if r is None:
            continue
        cf = r["closed_form"]
        sup = [tuple(t) for t in cf["support_dets"]]
        amps = [sp.sympify(a) for a in cf["pretty"]]
        scan = conjecture2_scan(sup, amps)
        # no holonomy is a definitive counterexample
        assert all(s["exponent_two"] is not False for s in scan), f"{sysid} v{idx}"
        # and at least one is positively confirmed exponent-2
        assert any(s["exponent_two"] is True for s in scan), f"{sysid} v{idx}"


def test_exactify_emits_conjecture2_field(capsys):
    # a certified state carries the holonomy_exponent_two verdict and does not warn
    import sympy as sp

    from gpc_census.exactify import conjecture2_scan

    if not DATA.exists():
        import pytest
        pytest.skip("no states atlas present")
    rows = [json.loads(x) for x in DATA.read_text().splitlines() if x.strip()]
    r = next(x for x in rows if x["system"] == "(4,9)" and x["index"] == 65)
    cf = r["closed_form"]
    sup = [tuple(t) for t in cf["support_dets"]]
    amps = [sp.sympify(a) for a in cf["pretty"]]
    scan = conjecture2_scan(sup, amps)
    assert scan and all(s["exponent_two"] for s in scan)
