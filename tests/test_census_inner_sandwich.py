"""Guards for the completeness-without-exhaustiveness argument.

The claim is P = O at every census system, proved by a sandwich whose two
halves are independent of candidate enumeration. The tests that matter are the
ones that could catch the argument being circular or vacuous:

- the outer half must be a REPLAYED certificate, not a citation, so the stored
  rows are re-screened live at (3,7);
- the inner half must not lean on the census vertex lists, since those were
  derived from the constraint tables; it leans on the attaining states, whose
  certificates reference no table, and one is re-verified live;
- a system with a resisting vertex must not be reported as closed.
"""
import json
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "census_inner_sandwich.json"

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))


def _art():
    return json.loads(ARTIFACT.read_text())


# --------------------------------------------------------------------------
# the sandwich closed, and on the whole N = 3 ladder
# --------------------------------------------------------------------------

def test_every_n3_system_closes():
    art = _art()
    for key in ("(3,6)", "(3,7)", "(3,8)", "(3,9)", "(3,10)"):
        rec = art["systems"].get(key)
        assert rec is not None, key
        assert rec["closure_proved"], key
        assert rec["resisting_vertices"] == [], key
        assert rec["attainment_failures"] == [], key
        assert rec["attained_vertices"] == rec["vertices"], key
    assert {"(3,6)", "(3,7)", "(3,8)", "(3,9)", "(3,10)"} <= set(art["closed_systems"])


def test_vertex_counts_match_the_published_census():
    """The enumeration is of O, and must land on the certified vertex counts."""
    art = _art()
    expected = {"(3,6)": 4, "(3,7)": 10, "(3,8)": 38, "(3,9)": 58, "(3,10)": 113}
    for key, count in expected.items():
        assert art["systems"][key]["vertices"] == count, key


def test_outer_half_is_certified_row_by_row():
    art = _art()
    for key in ("(3,7)", "(3,8)", "(3,9)", "(3,10)"):
        validity = art["systems"][key]["validity_P_subset_O"]
        assert validity["available"], key
        assert validity["all_valid"], key
        assert validity["unresolved"] == 0, key
        assert validity["exact_rejected"] == 0, key
        assert validity["certified"] + validity["structural"] == validity["rows"], key


# --------------------------------------------------------------------------
# the two halves, re-run live rather than re-read
# --------------------------------------------------------------------------

def test_stored_rows_screen_as_ressayre_certificates_live():
    """P subset O, replayed. A citation would not establish this."""
    from gpc_census.constraints import constraints
    from gpc_census.generation import screen_tau_candidates
    from gpc_census.generation.bdr import tau_from_constraint
    from gpc_census.generation.model import IntegralConstraint

    rows = constraints(3, 7)["inequalities"]
    taus = tuple(
        tau_from_constraint(IntegralConstraint(tuple(r["coeffs"]), int(r["rhs"])), 3)
        for r in rows
    )
    result = screen_tau_candidates(taus, ressayre_attempts=32,
                                   cyclic_cross_check=False,
                                   symbolic_fallback=True, workers=1)
    assert result.counts.certified_rows == len(rows) == 4
    assert result.counts.exact_rejected_rows == 0


def test_attaining_state_certificate_references_no_constraint_table():
    """O subset P, replayed on one state, which is where non-circularity lives."""
    import verify_states_standalone as vs

    records = [json.loads(line) for line
               in (ROOT / "states.jsonl").read_text().splitlines() if line.strip()]
    interference = [r for r in records
                    if r["system"] == "(3,10)" and r["classified"] == "INTERFERENCE"]
    assert interference
    ok, why = vs.verify(interference[0])
    assert ok, why


def test_vertex_enumeration_of_O_is_exact_and_reproducible():
    from fractions import Fraction

    from gpc_census import exact_polytope as ep
    from gpc_census import orbit

    amb_ineqs, amb_eqs = orbit.ambient_rows(3, 8)
    base_rows, eqs, base = orbit._chamber(8, 3, amb_eqs)
    vertices = ep.vertices_from_hrep(base, base_rows, eqs, amb_ineqs, 8)
    published = {tuple(Fraction(s) for s in r["spectrum"])
                 for r in json.loads((ROOT / "vertices/vertices_3_8.json").read_text())}
    assert set(vertices) == published
    assert len(vertices) == 38


# --------------------------------------------------------------------------
# the argument must be able to fail
# --------------------------------------------------------------------------

def test_a_resisting_vertex_blocks_closure():
    """A system is closed only when NOTHING resisted; check the guard exists."""
    art = _art()
    for key, rec in art["systems"].items():
        if rec.get("status") == "operational_failure":
            assert not rec["closure_proved"], key
            continue
        if rec["resisting_vertices"] or rec["attainment_failures"]:
            assert not rec["closure_proved"], key
        assert rec["closure_proved"] == (
            not rec["resisting_vertices"]
            and not rec["attainment_failures"]
            and (rec["validity_P_subset_O"].get("all_valid", False)
                 if rec["n"] == 3 else True)), key


def test_argument_block_states_its_own_non_circularity():
    art = _art()
    argument = art["argument"]
    assert "exhaustiveness" in argument["conclusion"]
    assert "circular" in argument["non_circularity"]


@pytest.mark.parametrize("key", ["(3,9)", "(3,10)"])
def test_the_expensive_ranks_are_the_new_content(key):
    """Ranks 9 and 10 had no in-repo completeness proof before this."""
    art = _art()
    rec = art["systems"][key]
    assert rec["closure_proved"]
    assert rec["outer_rows"] >= 52
    assert rec["states_reverified"]
