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


def test_splits_exact_eigenvalues():
    import sympy as sp
    from gpc_census.states import _splits
    for e1, e2 in [(14, 4), (20, 14), (20, 4)]:
        splits = _splits(e1, e2, 23)
        assert splits, (e1, e2)
        for a, b, x2 in splits:
            m = sp.Matrix([[a, sp.sqrt(x2)], [sp.sqrt(x2), b]])
            assert set(m.eigenvals()) == {sp.Integer(e1), sp.Integer(e2)}


def test_cancellation_polygon_filter():
    from gpc_census.states import _cancellation_feasible
    # a lone one-hop pair puts a single uncancellable term on modes (3, 4)
    assert not _cancellation_feasible((((0, 1, 2, 3), 5), ((0, 1, 2, 4), 5)))
    # no one-hop pairs at all: nothing to cancel
    assert _cancellation_feasible((((0, 1, 2, 3), 5), ((0, 4, 5, 6), 5)))
    # the same lone term is fine when it is the target of a block ansatz
    assert _cancellation_feasible((((0, 1, 2, 3), 5), ((0, 1, 2, 4), 5)),
                                  {(3, 4): 5.0})
    # but not when the target magnitude is unreachable
    assert not _cancellation_feasible((((0, 1, 2, 3), 5), ((0, 1, 2, 4), 5)),
                                      {(3, 4): 50.0})
    # and a target pair with no hop terms is unreachable
    assert not _cancellation_feasible((((0, 1, 2, 3), 5), ((0, 4, 5, 6), 5)),
                                      {(3, 4): 5.0})


def test_orbit_dedup_maps():
    from gpc_census.states import _ansatz_maps, _orbit_key
    nv = [20, 14, 14, 14, 14, 4, 4, 4, 4]
    maps = _ansatz_maps(nv, [], 9)
    assert len(maps) == 576  # S4 x S4 on the degenerate tails
    i1 = (((0, 1, 2, 5), 3), ((0, 3, 4, 6), 2))
    i2 = (((0, 1, 3, 7), 3), ((0, 2, 4, 5), 2))  # block-permuted image of i1
    assert _orbit_key(i1, maps) == _orbit_key(i2, maps)
    i3 = (((0, 1, 2, 5), 3), ((0, 3, 4, 6), 1))  # new weight multiset: new orbit
    assert _orbit_key(i1, maps) != _orbit_key(i3, maps)


def test_phase_gradient_matches_finite_differences():
    import numpy as np
    from gpc_census.states import _build
    np.random.seed(3)
    spec = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
    dets, a = _build(9, 4)
    idx = {t: i for i, t in enumerate(dets)}
    supd = [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 3, 5), (0, 1, 4, 5),
            (2, 3, 6, 7), (2, 4, 6, 7)]
    sup = [idx[t] for t in supd]
    mod = np.array([(k / 23) ** 0.5 for k in (7, 4, 3, 2, 5, 2)])
    asub = a[:, :, sup][:, :, :, sup]
    target = np.diag([float(x) for x in spec]).astype(complex)

    def fg(th):
        c = mod * np.exp(1j * np.concatenate(([0.0], th)))
        rho = np.einsum("p,mnpq,q->mn", c.conj(), asub, c)
        m_ = rho - target
        en = float(np.real(np.sum(m_ * m_.conj())))
        g1 = np.einsum("mn,mnpq,q->p", m_.conj(), asub, c)
        return en, 4.0 * np.imag(c.conj() * g1)[1:]

    th = np.random.uniform(0, 2 * np.pi, len(sup) - 1)
    _, g = fg(th)
    h = 1e-6
    for j in range(len(th)):
        tp, tm = th.copy(), th.copy()
        tp[j] += h
        tm[j] -= h
        assert abs((fg(tp)[0] - fg(tm)[0]) / (2 * h) - g[j]) < 1e-6
    assert max(abs(g)) > 1e-3  # the check is nontrivial


# The exact v_B interference state, reconstructed by the weights-first block
# sweep: a single 2x2 block mixing a 14-mode and a 4-mode. Support, weights,
# and block are ground truth (phase residual 1e-32); the block split (6,12)
# has off-diagonal sqrt(16)/23 = 4/23 and eigenvalues 14/23, 4/23.
_VB_BLOCK = (4, 8, 6, 12, 16)
_VB_NV = [20, 14, 14, 14, 6, 4, 4, 4, 12]
_VB_SUPPORT = [(0, 1, 2, 5), (0, 1, 3, 4), (0, 1, 3, 8), (0, 2, 3, 7),
               (0, 2, 4, 6), (0, 2, 6, 8), (1, 2, 7, 8), (2, 3, 4, 8)]
_VB_WEIGHTS = (4, 1, 8, 3, 3, 1, 1, 2)


def test_exactify_gauge_fixes_and_certifies_vB():
    # a v_B interference realization has a pi/8 phase lattice; exactify must
    # gauge-fix the single-particle U(1)^d freedom (here scrambled by a random
    # orbital phase rotation) and still certify the exact closed form.
    import math

    import numpy as np

    from gpc_census.exactify import exactify
    spec = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
    supp = [(0, 1, 2, 7), (0, 1, 3, 5), (0, 1, 4, 6), (0, 1, 6, 8),
            (0, 2, 3, 4), (0, 2, 3, 8), (1, 2, 4, 8), (1, 3, 7, 8)]
    ks = [3, 4, 2, 2, 4, 5, 2, 1]
    eighths = [0, 0, 1, -1, -1, 1, 0, 0]  # interference phase in units of pi/8
    np.random.seed(1)
    phi = np.random.uniform(0, 2 * math.pi, 9)  # arbitrary orbital phase gauge
    record = {"support": [[list(t), (k / 23) ** 0.5,
                           float(m * math.pi / 8 + sum(phi[x] for x in t))]
                          for t, k, m in zip(supp, ks, eighths)]}
    ex = exactify(4, 9, spec, record)
    assert ex["status"] == "EXACT"
    assert ex["weights"] == ks and ex["den"] == 23


def test_degenerate_closure_is_sound_superset_for_vB():
    # the block-target support filter must never drop true support (validation
    # law): the degenerate-signature closure is a sound superset.
    from gpc_census.states import admissible_support
    spec = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
    # ground-truth support for one valid v_B block ansatz (a different split
    # from _VB_BLOCK, exercising the same closure)
    supp = [(0, 1, 2, 5), (0, 1, 3, 4), (0, 1, 3, 8), (0, 2, 3, 7),
            (0, 2, 4, 6), (0, 2, 6, 8), (1, 2, 7, 8), (2, 3, 4, 8)]
    adm = set(admissible_support(4, 9, spec, groups="degenerate"))
    assert all(t in adm for t in supp)


def test_operator_selection_rule_is_unsound_at_degeneracy():
    # guards the documented reason the operator eigenprojection is not used:
    # the exact v_B state is NOT in the b-eigenspace of dGamma(U diag(a) U^T).
    import numpy as np
    from gpc_census.states import (_active_facets, _build, _natural_rotation,
                                    phase_solve)
    spec = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
    built = _build(9, 4)
    dets, a_tensor = built
    didx = {t: i for i, t in enumerate(dets)}
    w = [0] * len(dets)
    for t, k in zip(_VB_SUPPORT, _VB_WEIGHTS):
        w[didx[t]] = k
    target = np.diag([v / 23 for v in _VB_NV]).astype(complex)
    u, v_, a_, b_, x2 = _VB_BLOCK
    target[u, v_] = target[v_, u] = (x2 ** 0.5) / 23
    np.random.seed(0)
    psi, res = phase_solve(4, 9, spec, dets, 23, w, _built=built, target=target,
                           tries=8)
    assert res < 1e-12  # the analytic phase solver hits the block target exactly
    U = _natural_rotation(9, [_VB_BLOCK])
    worst = 0.0
    for a, b in _active_facets(4, 9, spec):
        anat = U @ np.diag([float(c) for c in a]) @ U.T
        dgamma = np.einsum("mn,mnij->ij", anat, a_tensor).real
        worst = max(worst, float(np.linalg.norm(dgamma @ psi - b * psi)))
    assert worst > 0.1  # eigenprojection would require this to be ~0
