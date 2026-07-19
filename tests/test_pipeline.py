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


def test_design_vertex_built_from_witness():
    # a DESIGN-INT vertex needs no iterative solve: its design witness gives a
    # real, phase-free state that is exact by construction. v_A is DESIGN-INT.
    from gpc_census.exactify import exactify
    from gpc_census.states import solve_design_vertex
    vA = [Fraction(x, 21) for x in (16, 16, 16, 6, 6, 6, 6, 6, 6)]
    rec = solve_design_vertex(4, 9, vA)
    assert rec is not None and rec["status"] == "OK"
    assert all(ph == 0.0 for _, _, ph in rec["support"])  # phase-free
    assert exactify(4, 9, vA, rec)["status"] == "EXACT"
    # an interference vertex is not DESIGN-INT, so the builder declines it
    vB = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
    assert solve_design_vertex(4, 9, vB) is None


def test_min_block_count_predicts_budget():
    # the block-budget preflight: 0 for a design, the exact budget for a
    # single-block interference vertex, and None for a vertex outside the
    # block-ansatz family (fail-fast frontier flag).
    from gpc_census.states import min_block_count
    vA = [Fraction(x, 21) for x in (16, 16, 16, 6, 6, 6, 6, 6, 6)]
    vB = [Fraction(x, 23) for x in (20, 14, 14, 14, 14, 4, 4, 4, 4)]
    assert min_block_count(4, 9, vA) == 0          # design: one-hop-free support
    assert min_block_count(4, 9, vB, max_blocks=2) == 1   # v_B is single-block
    # (9:6:5:5:5:2:2:1:1) is a 4_9 interference vertex off the family <= 2 blocks
    frontier = [Fraction(x, 9) for x in (9, 6, 5, 5, 5, 2, 2, 1, 1)]
    assert min_block_count(4, 9, frontier, max_blocks=2) is None


def test_schur_horn_generalizes_splits():
    # the k=2 Schur-Horn diagonals must contain every _splits pair (both orders)
    from gpc_census.states import _schur_horn_diagonals, _splits
    sh = set(_schur_horn_diagonals([14, 4], 23))
    sp = {(a, b) for a, b, _ in _splits(14, 4, 23)}
    sp |= {(b, a) for a, b, _ in _splits(14, 4, 23)}
    assert sp <= sh


def test_esym_and_grad_matches_charpoly():
    # elementary symmetric polys equal the eigenvalue symmetric functions, and
    # the end-to-end use is gradient-checked in test_clique_phase_gradient
    import numpy as np

    from gpc_census.states import _esym_and_grad
    np.random.seed(1)
    m = np.random.randn(3, 3)
    b = (m + m.T) / 2
    ev, _ = _esym_and_grad(b)
    eig = np.linalg.eigvalsh(b)
    e1 = eig.sum()
    e2 = eig[0] * eig[1] + eig[0] * eig[2] + eig[1] * eig[2]
    e3 = float(np.prod(eig))
    assert np.allclose(ev, [e1, e2, e3])


def test_clique_phase_gradient_matches_finite_differences():
    # the analytic gradient L-BFGS uses in phase_solve_clique (off-clique
    # cancellation plus char-poly block matching) must match finite differences
    import numpy as np

    from gpc_census.states import _build, _esym_and_grad
    np.random.seed(2)
    dets, a = _build(9, 4)
    idx = {t: i for i, t in enumerate(dets)}
    supd = [(0, 1, 4, 6), (0, 2, 3, 4), (0, 4, 5, 7), (1, 2, 3, 4), (1, 4, 5, 8)]
    sup = [idx[t] for t in supd]
    mod = np.array([(w / 9) ** 0.5 for w in (2, 2, 1, 3, 1)])
    asub = a[:, :, sup][:, :, :, sup]
    clique = (0, 1, 4)
    tgt = _esym_and_grad(np.diag([e / 9 for e in (9, 6, 5)]))[0]
    offc = [(p, q) for p in range(9) for q in range(p + 1, 9)
            if not {p, q} <= set(clique)]

    def fg(th):
        c = mod * np.exp(1j * np.concatenate(([0.0], th)))
        rho = np.einsum("p,mnpq,q->mn", c.conj(), asub, c)
        en = 0.0
        g = np.zeros((9, 9), complex)
        for (p, q) in offc:
            en += 2 * abs(rho[p, q]) ** 2
            g[p, q] += 2 * rho[p, q]
            g[q, p] += 2 * rho[q, p]
        b = rho[np.ix_(clique, clique)].real
        ev, dev = _esym_and_grad(b)
        r = ev - tgt
        en += 10 * float(r @ r)
        dedb = 10 * sum(2 * r[j] * dev[j] for j in range(3))
        for ia, ma in enumerate(clique):
            for ib, mb in enumerate(clique):
                g[ma, mb] += dedb[ia, ib]
        g1 = np.einsum("mn,mnpq,q->p", g, asub, c)
        return en, 2 * np.imag(c.conj() * g1)[1:]

    th = np.random.uniform(0, 2 * np.pi, len(sup) - 1)
    _, grad = fg(th)
    h = 1e-6
    for j in range(len(th)):
        tp, tm = th.copy(), th.copy()
        tp[j] += h
        tm[j] -= h
        assert abs((fg(tp)[0] - fg(tm)[0]) / (2 * h) - grad[j]) < 1e-6
    assert max(abs(grad)) > 1e-3


def test_min_clique_count_flags_k3_vertex():
    # (9:6:5:5:5:2:2:1:1)/9 needs a 3-clique (no 2x2 configuration reaches it)
    from gpc_census.states import min_clique_count
    spec = [Fraction(x, 9) for x in (9, 6, 5, 5, 5, 2, 2, 1, 1)]
    assert min_clique_count(4, 9, spec) == 3


def test_exactify_certifies_real_k3_state():
    # a k=3 clique interference vertex can have a REAL closed form (large
    # Schur-Horn fiber), which the gauge-fixed exactify certifies with no
    # extension. idx 24 of (4,9), (9:6:5:5:5:2:2:1:1)/9, is one such vertex.
    import math

    from gpc_census.exactify import exactify
    spec = [Fraction(x, 9) for x in (9, 6, 5, 5, 5, 2, 2, 1, 1)]
    support = [(0, 1, 2, 5), (0, 1, 3, 4), (0, 2, 3, 7), (0, 2, 4, 5),
               (0, 2, 6, 8), (0, 3, 4, 8)]
    ks = [1, 1, 1, 1, 2, 3]  # |c|^2 * 9
    record = {"support": [[list(t), math.sqrt(k / 9), 0.0]
                          for t, k in zip(support, ks)]}
    ex = exactify(4, 9, spec, record)
    assert ex["status"] == "EXACT"
    assert ex["weights"] == ks and ex["den"] == 9


def test_recognize_algebraic_beyond_psqrtq():
    # PSLQ recognizes a degree-2 algebraic NOT of p*sqrt(q)/r form
    import math

    from gpc_census.exactify import recognize_algebraic
    x = (math.sqrt(5) - 1) / 4  # = cos(2pi/5), root of 4y^2 + 2y - 1
    r = recognize_algebraic(x)
    assert r is not None and abs(float(r) - x) < 1e-10
    y = (1 + math.sqrt(13)) / 6  # root of 3y^2 - y - 1
    r2 = recognize_algebraic(y)
    assert r2 is not None and abs(float(r2) - y) < 1e-10


def test_multi_clique_ansatze_disjoint():
    # two-clique ansatze must have disjoint modes and the right eigenvalue mixes
    from gpc_census.states import multi_clique_ansatze
    spec = [Fraction(x, 10) for x in (10, 7, 5, 5, 5, 2, 2, 2, 2)]
    seen = 0
    for nv, cliques in multi_clique_ansatze(4, 9, spec, sizes=(3,), n_cliques=2):
        assert len(cliques) == 2
        m0, m1 = set(cliques[0][0]), set(cliques[1][0])
        assert not (m0 & m1)  # disjoint modes
        for modes, evals in cliques:
            assert len(set(evals)) == len(evals)  # distinct eigenvalue classes
        seen += 1
        if seen >= 50:
            break
    assert seen > 0
