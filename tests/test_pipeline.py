"""The engine must reproduce the historical results in results/data."""
import json
import pathlib
import shutil
from fractions import Fraction

import pytest

from gpc_census.constraints import constraints
from gpc_census.validate import check_embedding, check_physical, check_selfdual

DATA = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
HAS_LRS = shutil.which("lrs") is not None

EXPECTED_COUNTS = {(3, 6): 4, (3, 7): 10, (3, 8): 38, (4, 8): 22,
                   (3, 9): 58, (4, 9): 103, (3, 10): 113, (4, 10): 159, (5, 10): 292}


def test_constraint_table_counts():
    expected = {(3, 9): 52, (4, 9): 60, (3, 10): 93, (4, 10): 125, (5, 10): 161,
                (3, 7): 4, (3, 8): 31, (4, 8): 15, (3, 6): 1}
    for (n, d), count in expected.items():
        assert len(constraints(n, d)["inequalities"]) == count


def test_duality_transport():
    dual = constraints(5, 8)
    assert len(dual["inequalities"]) == 31
    assert "dual" in dual["source"]
    with pytest.raises(KeyError):
        constraints(3, 11)


@pytest.mark.skipif(not HAS_LRS, reason="lrslib not installed")
@pytest.mark.parametrize("n,d", sorted(EXPECTED_COUNTS))
def test_vertex_enumeration_reproduces_census(n, d):
    from gpc_census.polytope import vertices
    verts = vertices(n, d)
    assert len(verts) == EXPECTED_COUNTS[(n, d)]
    assert check_physical(verts, n) == []
    assert check_selfdual(verts, n, d) == []


@pytest.mark.skipif(not HAS_LRS, reason="lrslib not installed")
def test_embedding_invariant():
    from gpc_census.polytope import vertices
    assert check_embedding(vertices(3, 9), vertices(3, 10)) == []


@pytest.mark.skipif(not HAS_LRS, reason="lrslib not installed")
def test_vertices_match_historical_files():
    from gpc_census.polytope import vertices
    for (n, d) in [(3, 9), (5, 10)]:
        hist = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
        hset = {tuple(Fraction(s) for s in o["spectrum"]) for o in hist}
        assert set(vertices(n, d)) == hset


def test_attain_borland_dennis():
    import numpy as np
    from gpc_census.states import attain
    np.random.seed(0)
    spec = [Fraction(x, 4) for x in (3, 3, 2, 2, 1, 1)]
    _, res, _ = attain(3, 6, spec)
    assert res < 1e-12


def test_classify_and_certify_design():
    from gpc_census.certify import verify_design
    from gpc_census.classify import classify
    spec = [Fraction(1), Fraction(1), Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
    r = classify(3, 6, spec)
    assert r["verdict"] == "DESIGN-INT" and r["solver"] in ("cpsat", "cbc")
    assert verify_design(3, 6, spec, r["witness"])


def test_classify_cbc_fallback(monkeypatch):
    import gpc_census.classify as cl
    monkeypatch.setattr(cl, "_HAVE_ORTOOLS", False)
    spec = [Fraction(1), Fraction(1), Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
    r = cl.classify(3, 6, spec)
    assert r["verdict"] == "DESIGN-INT" and r["solver"] == "cbc"


def test_farkas_certifies_vB():
    from gpc_census.certify import farkas_interference
    spec = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
    y = farkas_interference(4, 9, spec)
    assert y is not None


def test_classify_backend_reported():
    from gpc_census.classify import classify_full
    rec = classify_full(3, 6, [Fraction(x, 4) for x in (3, 3, 2, 2, 1, 1)])
    assert rec["verdict"] == "DESIGN-INT"
    assert rec["backend"] in ("cpsat", "cbc")
    assert "witness" in rec


def test_classify_cbc_arm(monkeypatch):
    import gpc_census.classify as cl
    monkeypatch.setattr(cl, "BACKEND", "cbc")
    rec = cl.classify_full(3, 6, [Fraction(x, 4) for x in (3, 3, 2, 2, 1, 1)])
    assert rec["verdict"] == "DESIGN-INT" and rec["backend"] == "cbc"
