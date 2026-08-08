"""Guards for the exact closure-certificate verification layer."""

from __future__ import annotations

import importlib.util
import json
import pathlib
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

CLOSURE_CERTIFICATE = DATA / "closure_certificate.json"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = _load("closure_certificate")


def _artifact() -> dict:
    return json.loads(CLOSURE_CERTIFICATE.read_text())


def test_every_system_agrees_with_the_published_invariants():
    artifact = _artifact()
    assert artifact["summary"]["mismatches_against_polytope_invariants"] == []
    assert artifact["summary"]["systems_audited"] == 9
    for record in artifact["systems"].values():
        assert record["agrees_with_polytope_invariants"]
        assert record["facets"] == record["published_facets"]
        assert record["vertices"] == record["published_vertices"]
        assert record["certified_affine_dimension"] == record["published_dimension"]


def test_all_799_census_vertices_are_true_vertices():
    artifact = _artifact()
    total = 0
    for record in artifact["systems"].values():
        assert record["vertices_with_full_active_rank"] == record["vertices"]
        total += record["vertices"]
    assert total == 799
    assert artifact["summary"]["all_vertices_have_full_active_rank"]


def test_slack_matrix_rank_is_dimension_plus_one():
    artifact = _artifact()
    for record in artifact["systems"].values():
        assert record["slack_rank"] == record["certified_affine_dimension"] + 1
        assert record["slack_matrix_nonnegative"]
        assert record["slack_zero_pattern_is_incidence"]
    assert artifact["summary"]["all_slack_ranks_correct"]


def test_borland_dennis_hull_is_degenerate_and_collapses_rows():
    """The one system where canonicalizing rows into faces is load bearing."""
    record = _artifact()["systems"]["(3,6)"]
    assert record["certified_affine_dimension"] == 3
    assert record["affine_hull_equation_count"] == 3
    assert record["codimension_one_rows"] == 6
    assert record["facets"] == 4
    assert record["rows_collapsed_by_canonicalization"] == 2


def test_no_other_system_collapses_rows():
    for system, record in _artifact()["systems"].items():
        if system == "(3,6)":
            continue
        assert record["rows_collapsed_by_canonicalization"] == 0
        assert record["affine_hull_equation_count"] == 1


def test_validity_is_weaker_than_facetness():
    """Valid rows that are not facets must exist and be dropped."""
    artifact = _artifact()["systems"]
    assert artifact["(3,7)"]["valid_rows"] > artifact["(3,7)"]["facets"]
    assert artifact["(3,6)"]["valid_rows"] > artifact["(3,6)"]["codimension_one_rows"]


def test_particle_hole_applies_only_at_half_filling():
    artifact = _artifact()
    assert artifact["summary"]["half_filling_systems"] == ["(3,6)", "(4,8)", "(5,10)"]
    for system, record in artifact["systems"].items():
        n, d = (int(part) for part in system.strip("()").split(","))
        parity = record["particle_hole"]
        assert parity["applies"] == (2 * n == d)
        if parity["applies"]:
            assert parity["vertex_set_is_closed"]
            assert parity["parity_holds"]
            assert (
                parity["fixed_vertices"] + 2 * parity["paired_vertex_orbits"]
                == record["vertices"]
            )


def test_particle_hole_orbit_counts():
    artifact = _artifact()["systems"]
    assert artifact["(3,6)"]["particle_hole"]["fixed_vertices"] == 4
    assert artifact["(4,8)"]["particle_hole"]["paired_vertex_orbits"] == 8
    assert artifact["(5,10)"]["particle_hole"]["paired_vertex_orbits"] == 142


def test_affine_dimension_is_computed_not_assumed():
    """(3,6) sits in a 5-dimensional trace slice but has dimension 3."""
    vertices = MODULE.load_vertices(3, 6)
    assert MODULE.affine_dimension(vertices) == 3
    equations = MODULE.affine_hull_equations(vertices)
    assert len(equations) == 3
    for equation in equations:
        normal, offset = equation[:-1], equation[-1]
        for vertex in vertices:
            assert sum(a * b for a, b in zip(normal, vertex)) == offset


def test_rank_helper_is_exact():
    rows = [[Fraction(1), Fraction(2)], [Fraction(2), Fraction(4)]]
    assert MODULE.rank(rows) == 1
    assert MODULE.rank([[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]]) == 2
