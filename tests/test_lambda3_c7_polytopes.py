"""Guards for the orbit-closure (entanglement) moment polytope layer.

The validation law applies here as everywhere: no orbit polytope ships that
has not passed a structural gate. The gates are

1. the exact polytope tools reproduce the published ambient (3,7) vertex set
   from its H-representation,
2. the four Lambda^3 C^6 classes reproduce the published fermionic
   entanglement polytopes, and the open orbit reproduces the ambient (3,6)
   polytope exactly,
3. at d = 7 the open orbit reproduces the ambient (3,7) polytope exactly,
4. every certified inner point really is attained: its exact 1-RDM spectrum
   is recomputed from the state that is supposed to attain it,
5. inclusion of class polytopes is consistent with the certified degeneration
   order.

The d = 7 census is exercised through the shipped artifact rather than
recomputed, which keeps this file fast; ``scripts/lambda3_c7_polytopes.py``
regenerates it.
"""
import json
import pathlib
from fractions import Fraction as F

import pytest

from gpc_census import exact_polytope as ep
from gpc_census import orbit

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "lambda3_c7_polytopes.json"


def _ambient_vertices(n, d):
    return {tuple(F(s) for s in r["spectrum"])
            for r in json.loads((ROOT / "vertices" / f"vertices_{n}_{d}.json").read_text())}


# --------------------------------------------------------------------------
# gate 1: the exact polytope tools against the published ambient census
# --------------------------------------------------------------------------

def test_exact_tools_reproduce_ambient_3_7():
    amb_ineqs, amb_eqs = orbit.ambient_rows(3, 7)
    base_rows, eqs, base = orbit._chamber(7, 3, amb_eqs)
    verts = ep.vertices_from_hrep(base, base_rows, eqs, amb_ineqs, 7)
    assert set(verts) == _ambient_vertices(3, 7)


def test_simplex_feasibility_agrees_with_membership():
    pts = [(F(1), F(0)), (F(0), F(1)), (F(0), F(0))]
    assert ep.in_convex_hull((F(1, 3), F(1, 3)), pts)
    assert not ep.in_convex_hull((F(2, 3), F(2, 3)), pts)
    assert not ep.in_convex_hull((F(-1, 5), F(0)), pts)


def test_rank_handles_dependent_equality_block():
    # the (3,6) pairing equalities imply the trace equality; counting rows
    # instead of rank once produced an empty vertex set and a false EXACT
    _, amb_eqs = orbit.ambient_rows(3, 6)
    rows = [list(a) for a, _ in amb_eqs] + [[F(1)] * 6]
    assert ep.rank(rows) == 3


# --------------------------------------------------------------------------
# gate 2/3: the open orbit must reproduce the ambient polytope
# --------------------------------------------------------------------------

def test_d6_census_matches_published_entanglement_polytopes():
    res = orbit.census(6)
    assert all(r["verdict"] == "EXACT" for r in res.values())
    got = {r["label"]: {tuple(F(x) for x in v) for v in r["vertices"]}
           for r in res.values()}
    assert got["I"] == {(F(1), F(1), F(1), F(0), F(0), F(0))}
    assert got["II"] == {(F(1), F(1), F(1), F(0), F(0), F(0)),
                         (F(1), F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(0))}
    # the open orbit is the whole ambient (3,6) polytope
    assert got["IV"] == _ambient_vertices(3, 6)
    # the tangent class is strictly smaller, and strictly larger than rank 5
    assert got["III"] < got["IV"] | got["III"]
    assert got["II"] < got["III"]


def test_ressayre_manifest_agrees_with_ambient_rows():
    """Cross-check against the independently derived stage-1 Ressayre manifest.

    The manifest derives the (3,7) system from first principles; this layer
    reads it from the Altunbulak-Klyachko table. They must agree as an
    oriented set, and the manifest's rows alone must cut out the published
    vertices. This also pins the scope of that manifest: it is the AMBIENT
    cone, hence the OPEN-orbit polytope, and so bears on class IX only. It
    does not supply orbit-closure inequalities for the bracketed classes.
    """
    path = ROOT / "stage1_ressayre_3_7.json"
    if not path.exists():
        pytest.skip("stage-1 Ressayre manifest not present")
    rows = [r["constraint"] for r in json.loads(path.read_text())["retained_rows"]]
    amb_ineqs, amb_eqs = orbit.ambient_rows(3, 7)
    assert not amb_eqs
    mine = {(tuple(-x for x in a), -c) for a, c in amb_ineqs}
    theirs = {(tuple(F(x) for x in r["coeffs"]), F(r["rhs"])) for r in rows}
    assert mine == theirs

    extra = [([F(-x) for x in r["coeffs"]], F(-r["rhs"])) for r in rows]
    base_rows, eqs, base = orbit._chamber(7, 3, [])
    verts = set(ep.vertices_from_hrep(base, base_rows, eqs, extra, 7))
    assert verts == _ambient_vertices(3, 7)

    generic = json.loads(ARTIFACT.read_text())["blocks"]["7"]["classes"]["35"]
    assert {tuple(F(x) for x in v) for v in generic["vertices"]} == verts


def test_d7_open_orbit_is_the_ambient_polytope():
    block = json.loads(ARTIFACT.read_text())["blocks"]["7"]
    generic = block["classes"]["35"]
    assert generic["verdict"] == "EXACT"
    assert {tuple(F(x) for x in v) for v in generic["vertices"]} == _ambient_vertices(3, 7)


# --------------------------------------------------------------------------
# gate 4: certified inner points are really attained
# --------------------------------------------------------------------------

@pytest.mark.parametrize("dim", sorted(orbit.CLASSES_7))
def test_representative_spectrum_is_reproduced_two_ways(dim):
    psi = orbit.CLASSES_7[dim]["rep"]
    assert orbit.orbit_dimension(psi, 7) == dim
    assert orbit.support_rank(psi, 7) == orbit.CLASSES_7[dim]["rank"]
    # the diagonal shortcut and the characteristic polynomial must agree
    assert orbit.moment_point(psi, 7) == orbit.rational_spectrum(psi, 7)


def test_newton_cloud_points_are_attained_by_an_explicit_state():
    # A Newton weight w_S is the SQUARED amplitude, so the state realising it
    # is sum_S sqrt(w_S) e_S. Using coefficients whose squares are the weights
    # keeps the state rational too, and its spectrum must be the cloud point.
    psi = orbit.CLASSES_7[28]["rep"]
    sup = sorted(psi)
    cloud = orbit.newton_cloud(sup, 7, levels=(4, 5, 6, 7))
    for coeffs in ((1, 1, 1, 1), (2, 1, 1, 0), (3, 2, 1, 1)):
        weights = [c * c for c in coeffs]
        total = sum(weights)
        phi = {S: F(c) for S, c in zip(sup, coeffs) if c}
        expected = [F(0)] * 7
        for w, S in zip(weights, sup):
            for i in S:
                expected[i] += F(w, total)
        expected = tuple(sorted(expected, reverse=True))
        assert orbit.moment_point(phi, 7) == expected
        assert orbit.rational_spectrum(phi, 7) == expected
    # and the all-ones weighting really is one of the emitted cloud points
    even = [F(0)] * 7
    for S in sup:
        for i in S:
            even[i] += F(1, len(sup))
    assert tuple(sorted(even, reverse=True)) in cloud


def test_one_hop_free_supports_have_diagonal_rdm():
    for sup in orbit.canonical_supports_by_class(7)[35]:
        phi = {S: F(1) for S in sup}
        assert orbit.moment_point(phi, 7) is not None


def test_every_one_hop_free_support_has_trivial_coefficient_moduli():
    # this is what makes the support enumeration complete rather than sampled
    for sup in orbit.onehop_free_supports(6):
        assert orbit.coefficient_moduli(sup, 6) == 0


# --------------------------------------------------------------------------
# gate 5: structure of the shipped d = 7 result
# --------------------------------------------------------------------------

def test_artifact_structure_and_invariants():
    block = json.loads(ARTIFACT.read_text())["blocks"]["7"]
    assert block["n_classes"] == 9
    slater = (F(1), F(1), F(1), F(0), F(0), F(0), F(0))
    ambient = _ambient_vertices(3, 7)
    for key, rec in block["classes"].items():
        verts = [tuple(F(x) for x in v) for v in rec["vertices"]]
        assert verts, key
        for v in verts:
            assert sum(v) == 3, key
            assert all(v[i] >= v[i + 1] for i in range(6)), key
            assert v[0] <= 1 and v[-1] >= 0, key
            assert ep.in_convex_hull(v, sorted(ambient)), key
        # the closed orbit sits inside every orbit closure, so the Slater
        # spectrum belongs to every class polytope
        assert ep.in_convex_hull(slater, verts), key
        # a rank-r class cannot occupy more than r orbitals
        assert all(v[i] == 0 for v in verts for i in range(rec["rank"], 7)), key


def test_degeneration_order_implies_polytope_inclusion():
    block = json.loads(ARTIFACT.read_text())["blocks"]["7"]
    verts = {k: [tuple(F(x) for x in v) for v in r["vertices"]]
             for k, r in block["classes"].items()}
    for key, lower in block["degeneration_order"].items():
        for other in lower:
            # a certified degeneration forces containment of the polytopes
            for p in verts[str(other)]:
                assert ep.in_convex_hull(p, verts[key]), (key, other)


def test_bracketed_classes_declare_their_gap():
    block = json.loads(ARTIFACT.read_text())["blocks"]["7"]
    for key, rec in block["classes"].items():
        if rec["verdict"] == "EXACT":
            assert not rec["outer_vertices_not_reached"], key
        else:
            assert rec["verdict"] == "BRACKET"
            assert rec["outer_vertices_not_reached"], key
        # certificates must never contradict the outer system
        assert not rec["certificate_violations"], key
