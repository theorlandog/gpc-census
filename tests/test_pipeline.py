"""The engine must reproduce the historical results in results/data."""
import json
import pathlib
import shutil
from fractions import Fraction

import pytest

from gpc_census.constraints import constraints
from gpc_census.validate import check_embedding, check_physical, check_selfdual

DATA = pathlib.Path(__file__).resolve().parent / "data"
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
def test_classify_reproduces_vA_vB():
    from gpc_census.classify import classify_full
    va = classify_full(4, 9, [Fraction(x, 21) for x in (16, 16, 16, 6, 6, 6, 6, 6, 6)])
    vb = classify_full(4, 9, [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)])
    assert va["verdict"] == "DESIGN-INT"
    assert vb["verdict"] == "INTERFERENCE"


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


def test_exactify_recognition_layer():
    import math

    from gpc_census.exactify import recognize_cos, snap_moduli
    r = recognize_cos(3 / (4 * math.sqrt(14)))
    assert r is not None and str(r) in ("3*sqrt(14)/56", "sqrt(14)*3/56")
    assert snap_moduli([math.sqrt(4 / 23), math.sqrt(7 / 23)], 23) == [4, 7]
