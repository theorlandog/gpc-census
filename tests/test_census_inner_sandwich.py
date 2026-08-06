"""Guards for the completeness-without-exhaustiveness argument.

The claim is P = O at every census system, proved by a sandwich whose two
halves are independent of candidate enumeration. The tests that matter are the
ones that could catch the argument being circular or vacuous:

- the outer half must be a REPLAYED certificate, not a citation, so the stored
  rows are re-screened live at (3,7);
- the inner half must not lean on the census vertex lists, since those were
  derived from the constraint tables; it leans on the attaining states, whose
  certificates reference no table, and one is re-verified live;
- a system with a resisting vertex must not be reported as closed;
- the outer half must certify at the system's OWN particle number, and must
  still reject perturbed rows there, or the general-N extension would be
  plumbing rather than proof.
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


def test_fast_adjacency_agrees_with_the_naive_definition():
    """The DD speedup must be an identity, not an approximation.

    ``cut`` memoises the exact rank on the common tight set and rejects pairs
    that cannot reach rank ``n-1`` on an integer count. Both are identities, so
    the cut it produces must equal the one produced by testing every pair with
    the unmemoised ``_adjacent``.
    """
    from gpc_census import exact_polytope as ep
    from gpc_census import orbit

    n, d = 3, 8
    amb_ineqs, amb_eqs = orbit.ambient_rows(n, d)
    base_rows, eqs, verts = orbit._chamber(d, n, amb_eqs)
    rows = list(base_rows)
    for row in amb_ineqs:
        a, c = row
        val = {v: sum(ai * xi for ai, xi in zip(a, v)) for v in verts}
        keep = [v for v in verts if val[v] >= c]
        drop = [v for v in verts if val[v] < c]
        naive = set(keep)
        for v in keep:
            if val[v] == c:
                continue
            for w in drop:
                if not ep._adjacent(rows, eqs, v, w, d):
                    continue
                theta = (c - val[w]) / (val[v] - val[w])
                naive.add(tuple(theta * vi + (1 - theta) * wi
                                for vi, wi in zip(v, w)))
        verts = ep.cut(verts, rows, eqs, row, d)
        assert set(verts) == naive
        rows.append(row)
        if not verts:
            break
    assert len(verts) == 38


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
            rec["inner_half_proved"] and rec["outer_half_proved"]), key
        assert rec["inner_half_proved"] == (
            not rec["resisting_vertices"] and not rec["attainment_failures"]), key


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


def test_the_inner_half_only_label_can_never_be_counted_as_closed():
    """The honesty invariant, asserted structurally rather than by system.

    An earlier version reported the N > 3 systems as closed while their outer
    half was only cited, because the screening pipeline was fixed-N=3. The
    pipeline now takes a particle number, so those systems close for real; the
    label that recorded the gap must still be impossible to confuse with a
    proof, whatever future system it applies to.
    """
    art = _art()
    for key, rec in art["systems"].items():
        if rec.get("status") == "operational_failure":
            continue
        assert rec["closure_proved"] == (
            rec["inner_half_proved"] and rec["outer_half_proved"]), key
        if rec["closure_status"] == "inner_half_only":
            assert not rec["outer_half_proved"], key
            assert not rec["closure_proved"], key
            assert key not in art["closed_systems"], key
            assert key in art["inner_half_only_systems"], key
    for key in art["inner_half_only_systems"]:
        assert key not in art["closed_systems"], key


def test_every_census_system_closes_including_the_wider_ones():
    """All nine census systems, not just the N = 3 ladder."""
    art = _art()
    expected = {"(3,6)", "(3,7)", "(3,8)", "(3,9)", "(3,10)",
                "(4,8)", "(4,9)", "(4,10)", "(5,10)"}
    assert set(art["closed_systems"]) == expected
    assert art["inner_half_only_systems"] == []
    for key in expected:
        rec = art["systems"][key]
        assert rec["outer_half_proved"] and rec["inner_half_proved"], key
        assert rec["validity_P_subset_O"]["particle_number"] == rec["n"], key


def test_wider_system_rows_screen_live_at_their_own_particle_number():
    """P subset O at N = 4, replayed.

    This is the half that used to be a citation. It is only worth anything if
    the general-N screening is discriminating, so a perturbation of every row
    must be rejected by the same call.
    """
    from gpc_census.constraints import constraints
    from gpc_census.generation import screen_tau_candidates
    from gpc_census.generation.bdr import tau_from_constraint
    from gpc_census.generation.model import IntegralConstraint

    rows = constraints(4, 8)["inequalities"]
    taus = tuple(
        tau_from_constraint(IntegralConstraint(tuple(r["coeffs"]), int(r["rhs"])), 4)
        for r in rows
    )
    good = screen_tau_candidates(taus, particle_number=4, ressayre_attempts=32,
                                 cyclic_cross_check=False, symbolic_fallback=True,
                                 workers=1)
    assert good.counts.certified_rows == len(rows) == 15
    assert good.counts.exact_rejected_rows == 0
    assert good.counts.evaluation_unresolved_rows == 0

    bumped = tuple((tau[0] + 1,) + tuple(tau[1:]) for tau in taus)
    bad = screen_tau_candidates(bumped, particle_number=4, ressayre_attempts=32,
                                cyclic_cross_check=False, symbolic_fallback=True,
                                workers=1)
    assert bad.counts.certified_rows == 0
    assert bad.counts.exact_rejected_rows == len(rows)


def test_screening_at_the_wrong_particle_number_does_not_certify():
    """The particle number is load bearing, not decorative.

    Screening the (4,8) rows as if they were N = 3 must not certify them all;
    otherwise the parameter would be doing nothing and the wider result would
    be an artifact of the plumbing.
    """
    from gpc_census.constraints import constraints
    from gpc_census.generation import screen_tau_candidates
    from gpc_census.generation.bdr import tau_from_constraint
    from gpc_census.generation.model import IntegralConstraint

    rows = constraints(4, 8)["inequalities"]
    taus = tuple(
        tau_from_constraint(IntegralConstraint(tuple(r["coeffs"]), int(r["rhs"])), 4)
        for r in rows
    )
    wrong = screen_tau_candidates(taus, particle_number=3, ressayre_attempts=32,
                                  cyclic_cross_check=False, symbolic_fallback=True,
                                  workers=1)
    assert wrong.counts.certified_rows < len(rows)


def test_full_dimension_witness_covers_the_wider_systems():
    """The BDR moment-cone theorem needs a nonempty interior at every system."""
    from gpc_census.generation.full_dimension import (
        require_full_dimension,
        verify_full_dimension_certificate,
    )

    for n, d in ((3, 7), (4, 8), (4, 9), (4, 10), (5, 10)):
        certificate = require_full_dimension(d, n)
        assert certificate.is_full_dimensional, (n, d)
        assert certificate.n == n
        assert verify_full_dimension_certificate(certificate), (n, d)
