"""Guards for the cross-realization invariant census.

Structural identities are checked on every row of the artifact; the expensive
exact recomputation is re-run only for the four `(3,6)` vertices, which is
enough to catch a broken engine without making the suite slow.
"""
from __future__ import annotations

import json
import pathlib

import pytest

from gpc_census.symplectic import (
    affine_dim,
    facet_supports,
    hull_lattice,
    stabilizer_invariants,
    tangent_cone_invariants,
    toral_stabilizer_dim_from_support,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" if (ROOT / "tests" / "data").exists() else ROOT / "results" / "data"
ARTIFACT = DATA / "symplectic_invariants.json"
QUANT = DATA / "quantization_multiplicities.json"


@pytest.fixture(scope="module")
def art():
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def states():
    return {(r["system"], r["index"]): r
            for r in (json.loads(x)
                      for x in (DATA / "states.jsonl").read_text().splitlines())}


def test_covers_the_whole_census(art):
    assert len(art["rows"]) == 799
    assert len({(r["system"], r["index"]) for r in art["rows"]}) == 799
    assert all(r["evidence"] == "exact" for r in art["rows"])


def test_orbit_and_stabilizer_dimensions_are_complementary(art):
    for r in art["rows"]:
        s = r["stabilizer"]
        assert s["stab_dim"] + s["orbit_dim"] == r["d"] ** 2, r


def test_stabilizer_sits_inside_the_commutant(art):
    """`U psi = e^{i t} psi` forces `U rho U^* = rho`, so the two agree exactly."""
    bad = [r for r in art["rows"]
           if r["stabilizer"]["stab_dim"] != r["stabilizer"]["commutant_stab_dim"]]
    assert not bad, [f"{r['system']}#{r['index']}" for r in bad]


def test_commutant_dimension_reproduces_the_census_spectrum(art):
    """Independent confirmation that the 1-RDM eigenvalues are the vertex."""
    bad = [r for r in art["rows"]
           if not r["stabilizer"]["commutant_dim_matches_multiplicities"]]
    assert not bad, [f"{r['system']}#{r['index']}" for r in bad]


def test_toral_stabilizer_is_the_combinatorial_prediction(art):
    """P4: the toral stabilizer is fixed by the support, with no amplitude input."""
    bad = [r for r in art["rows"]
           if r["stabilizer"]["toral_stab_dim"] != r["toral_stab_dim_predicted"]]
    assert not bad, [f"{r['system']}#{r['index']}" for r in bad]


def test_slater_vertices_are_grassmannian(art):
    seen = 0
    for r in art["rows"]:
        if set(r["integer_form"]) <= {0, r["denominator"] or 1}:
            seen += 1
            n, d = r["n"], r["d"]
            assert r["stabilizer"]["stab_dim"] == n ** 2 + (d - n) ** 2, r["system"]
            assert r["stabilizer"]["orbit_dim"] == 2 * n * (d - n), r["system"]
    assert seen == 9


def test_simple_vertices_match_the_published_polytope_invariants(art):
    published = json.loads((ROOT / "docs" / "polytope_invariants.json").read_text())
    for system, entry in published["systems"].items():
        rows = [r for r in art["rows"] if r["system"] == system]
        assert len(rows) == entry["vertices"], system
        assert rows[0]["dim"] == entry["dim"], system
        assert sum(1 for r in rows if r["cone"]["simple"]) == entry["simple_vertices"], system


def test_simple_iff_simplicial_tangent_cone(art):
    bad = [f"{r['system']}#{r['index']}" for r in art["rows"]
           if r["cone"]["simple"] != r["cone"]["simplicial_tangent_cone"]]
    assert not bad, bad


def test_cone_rays_are_primitive_and_span_correctly(art):
    from math import gcd
    for r in art["rows"]:
        cone = r["cone"]
        assert cone["edges"] == len(cone["rays"]) == len(cone["neighbours"])
        for ray in cone["rays"]:
            g = 0
            for x in ray:
                g = gcd(g, abs(x))
            assert g == 1, r
            assert sum(ray) == 0, r
        assert cone["ray_rank"] == r["dim"], r
        if cone["simplicial_tangent_cone"]:
            assert cone["lattice_index"] >= 1
            product = 1
            for f in cone["invariant_factors"]:
                product *= f
            assert product == cone["lattice_index"], r


def test_scorecard_verdicts_are_the_recorded_ones(art):
    card = art["scorecard"]
    assert card["P1_rho_diagonal_count"]["verdict"] == "PASS"
    assert card["P2_slater_vertices_are_grassmannian"]["verdict"] == "PASS"
    assert card["P3_stabilizer_never_trivial"]["verdict"] == "PASS"
    assert card["P4_toral_stabilizer_is_combinatorial"]["verdict"] == "PASS"
    assert card["P5_cone_types"]["verdict"] == "PASS"
    assert not card["internal_checks"]["stab_contained_in_commutant_violations"]
    assert not card["internal_checks"]["commutant_dim_vs_multiplicity_violations"]


@pytest.mark.parametrize("index", [0, 1, 2, 3])
def test_recompute_rank_six_exactly(art, states, index):
    row = next(r for r in art["rows"] if r["system"] == "(3,6)" and r["index"] == index)
    got = stabilizer_invariants(states[("(3,6)", index)])
    for key in ("stab_dim", "orbit_dim", "toral_stab_dim", "commutant_dim",
                "commutant_orbit_dim", "rho_is_diagonal"):
        assert got[key] == row["stabilizer"][key], (index, key)
    dets = states[("(3,6)", index)]["closed_form"]["support_dets"]
    assert toral_stabilizer_dim_from_support(dets, 6) == row["toral_stab_dim_predicted"]


def test_rank_six_cone_rays_regenerate(art):
    """Re-derive the `(3,6)` tangent cones from the V-rep through the library."""
    from fractions import Fraction as F

    verts = json.loads((DATA / "vertices" / "vertices_3_6.json").read_text())
    points = [[F(x) for x in v["spectrum"]] for v in verts]
    dim = affine_dim(points)
    facets = facet_supports(points, 3, 6, dim)
    lattice = hull_lattice(points, 6)
    for r in sorted((x for x in art["rows"] if x["system"] == "(3,6)"),
                    key=lambda x: x["index"]):
        cone = tangent_cone_invariants(points, r["index"], facets, lattice, dim)
        assert cone["lattice_index"] == r["cone"]["lattice_index"], r["index"]
        assert cone["rays"] == r["cone"]["rays"], r["index"]
        assert cone["invariant_factors"] == r["cone"]["invariant_factors"], r["index"]


@pytest.mark.skipif(not QUANT.exists(), reason="quantization artifact absent")
def test_quantization_multiplicities_are_zero_or_one(art):
    q = json.loads(QUANT.read_text())
    for row in q["rows"]:
        for entry in row["multiplicities"]:
            assert entry["multiplicity"] in (0, 1), row


@pytest.mark.skipif(not QUANT.exists(), reason="quantization artifact absent")
def test_rank_six_uniform_vertex_has_quantization_period_four():
    """The one census vertex whose quantization period beats its denominator.

    `lambda = (1/2, ..., 1/2)` in `(3,6)` has denominator 2, but the multiplicity
    is nonzero only at `k = 0 mod 4`, matching the classical fact that the
    `SL(6)` invariant ring of `wedge^3 C^6` is generated in degree 4.
    """
    q = json.loads(QUANT.read_text())
    row = next(r for r in q["rows"]
               if r["system"] == "(3,6)" and r["integer_form"] == [1, 1, 1, 1, 1, 1])
    assert row["denominator"] == 2
    for entry in row["multiplicities"]:
        assert entry["multiplicity"] == (1 if entry["k"] % 4 == 0 else 0), entry
    assert row.get("quantization_period") == 4
