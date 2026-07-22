"""v96 is solvable: a census false negative, now resolved in the dataset.

(3,10) v96 = (5,5,5,5,2,1,1,1,1,1)/9 was recorded SOLVE-FAIL in the census, but
exact closed-form extremal states for it DO exist. They were found by the hybrid
family (scripts/hybrid_search.py driving scripts/polygon_target.py) on the
degenerate (5,1) block ansatz that the census's own block_ansatze generates --
which the earlier solve_vertex_exact_first missed because the unsound
support-filter preflight excluded the true-solution determinants.

With that filter removed the audit recovered v96 (and nine siblings), and the
shipped dataset now records it certified. This test pins both ends of the
resolution: the stored states pass the shipped exact gate verify_exact and an
independent from-scratch 1-RDM spectrum check, and states.jsonl now labels v96
OK with a closed form. See docs/RESEARCH.md, "v96 solved".
"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))


def test_stored_v96_states_are_exact():
    import verify_hybrid_state as vh

    from gpc_census.exactify import verify_exact  # noqa: F401  (import sanity)

    art = ROOT / "docs" / "hybrid_cracks" / "v96.jsonl"
    records = [json.loads(x) for x in art.read_text().splitlines() if x.strip()]
    assert records, "no stored v96 states"
    # verify the first stored state both ways (independent 1-RDM + shipped gate)
    import sympy as sp
    from fractions import Fraction

    h = records[0]
    N, d, den = h["system"][0], h["system"][1], h["den"]
    dets = [tuple(x) for x in h["dets"]]
    amps_sp = [sp.sympify(a) for a in h["amplitudes"]]
    amps_num = [complex(sp.N(a, 40)) for a in amps_sp]
    ev, herm, trace = vh.independent_spectrum(dets, amps_num, d, den)
    assert herm < 1e-9 and abs(trace - N) < 1e-9
    assert [round(x) for x in ev] == sorted(h["spec"], reverse=True)
    assert verify_exact(N, d, [Fraction(x, den) for x in h["spec"]], dets, amps_sp)


def test_v96_is_recorded_certified():
    # the audit merged the crack into the census: v96 is now OK with a closed form
    states = (ROOT / "results" / "data" / "states.jsonl")
    if not states.exists():
        return
    for line in states.read_text().splitlines():
        d = json.loads(line)
        if d.get("system") == "(3,10)" and d.get("index") == 96:
            assert d["status"] == "OK"
            assert d.get("closed_form"), "v96 certified but carries no closed form"
            return
    raise AssertionError("v96 not found in states.jsonl")
