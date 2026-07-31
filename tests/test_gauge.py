"""Natural-orbital gauge minimization and the 2-RDM gate.

The invariance tests here are the load-bearing ones: they pin the claim that a
block-diagonal rotation is a gauge motion (1-RDM identically unchanged) and that
the 2-RDM spectrum is a basis-free fingerprint. Everything the module reports
about supports rests on those two facts.
"""
from __future__ import annotations

import itertools
import json
import pathlib

import numpy as np

from gpc_census.fiber import one_rdm
from gpc_census.gauge import (
    GaugeOrbit,
    compound_matrix,
    degenerate_blocks,
    descend,
    exterior_flattening_certificates,
    flattening_lower_bound,
    greedy_gauge_drop,
    profile_lower_bound,
    random_block_unitary,
    reachable_dets,
    same_state_up_to_orbital_rotation,
    solve_orbital_rotation,
    two_rdm_spectrum,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results/data"


def _v103():
    for line in (DATA / "states.jsonl").read_text().splitlines():
        r = json.loads(line)
        if r["system"] == "(3,10)" and r["index"] == 103:
            import sympy as sp

            cf = r["closed_form"]
            c = np.array([complex(sp.N(sp.sympify(s), 30)) for s in cf["pretty"]])
            dets = [tuple(t) for t in cf["support_dets"]]
            lam = [v / r["denominator"] for v in r["integer_form"]]
            return c, dets, lam
    raise AssertionError("v103 missing from the shipped ledger")


def test_compound_is_a_unitary_representation():
    """Gamma(U)Gamma(V) = Gamma(UV) and Gamma(U) unitary: the fermionic signs.

    Everything downstream is the compound action, so its sign bookkeeping is
    checked here structurally rather than against a second implementation of the
    same convention.
    """
    d, n = 7, 3
    dets = list(itertools.combinations(range(d), n))
    rng = np.random.default_rng(0)

    def haar(seed):
        z = (rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
        q, r = np.linalg.qr(z)
        return q * (np.diag(r) / np.abs(np.diag(r)))

    u, v = haar(0), haar(1)
    gu = compound_matrix(u, dets, dets)
    gv = compound_matrix(v, dets, dets)
    assert np.max(np.abs(gu.conj().T @ gu - np.eye(len(dets)))) < 1e-12
    assert np.max(np.abs(compound_matrix(u @ v, dets, dets) - gu @ gv)) < 1e-12


def test_block_rotation_is_a_gauge_motion():
    """1-RDM identically diag(lambda) along the whole orbit: no solver, no drift."""
    c, dets, lam = _v103()
    blocks = degenerate_blocks(lam)
    d = len(lam)
    target = reachable_dets(dets, blocks, d, 3)
    rng = np.random.default_rng(1)
    for _ in range(3):
        u = random_block_unitary(blocks, d, rng)
        rotated = compound_matrix(u, target, dets) @ c
        rho = one_rdm(rotated, target, d)
        assert np.max(np.abs(rho - np.diag(lam))) < 1e-12
        assert abs(np.linalg.norm(rotated) - 1.0) < 1e-12


def test_two_rdm_spectrum_is_basis_free():
    """The gate's premise: 2-RDM spectrum invariant under ANY one-body rotation."""
    c, dets, lam = _v103()
    d = len(lam)
    base = two_rdm_spectrum(c, dets, d)
    assert abs(base.sum() - 3.0) < 1e-10  # tr rho_2 = N(N-1)/2
    rng = np.random.default_rng(2)
    z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    q, r = np.linalg.qr(z)
    u = q * (np.diag(r) / np.abs(np.diag(r)))  # a FULL U(d) rotation, not a gauge one
    allbut = list(itertools.combinations(range(d), 3))
    rotated = compound_matrix(u, allbut, dets) @ c
    assert np.max(np.abs(two_rdm_spectrum(rotated, allbut, d) - base)) < 1e-10


def test_two_rdm_gate_separates_genuinely_different_states():
    """A different state must fail the gate, or the gate is vacuous."""
    c, dets, lam = _v103()
    d = len(lam)
    other = c.copy()
    other[0], other[1] = other[1], other[0]  # a different state, same moduli set
    same, dev = same_state_up_to_orbital_rotation(c, dets, other, dets, d)
    assert not same and dev > 1e-6


def test_profile_lemma_bounds_and_certifies():
    """Support >= occupied profile classes, with equality certifying minimality."""
    c, dets, lam = _v103()
    blocks = degenerate_blocks(lam)
    assert profile_lower_bound(dets, blocks) == 3
    single = [dets[0]]
    assert profile_lower_bound(single, blocks) == 1


def test_flattening_rank_is_a_gauge_invariant():
    """Per-class flattening ranks must not move along the orbit."""
    c, dets, lam = _v103()
    blocks = degenerate_blocks(lam)
    d = len(lam)
    before = flattening_lower_bound(c, dets, blocks)
    rng = np.random.default_rng(3)
    target = reachable_dets(dets, blocks, d, 3)
    u = random_block_unitary(blocks, d, rng)
    rotated = compound_matrix(u, target, dets) @ c
    live = [i for i in range(len(target)) if abs(rotated[i]) > 1e-12]
    after = flattening_lower_bound(rotated[live], [target[i] for i in live], blocks)
    assert before == after


def test_exterior_contraction_closes_v158_exactly():
    """The new interior flattening has exact rank 21, forcing support four."""
    import sympy as sp

    record = next(
        json.loads(line)
        for line in (DATA / "states.jsonl").read_text().splitlines()
        if (row := json.loads(line))["system"] == "(4,10)" and row["index"] == 158
    )
    amplitudes = [sp.sympify(value) for value in record["closed_form"]["pretty"]]
    dets = [tuple(t) for t in record["closed_form"]["support_dets"]]
    spectrum = [v / record["denominator"] for v in record["integer_form"]]
    certs = exterior_flattening_certificates(
        amplitudes,
        dets,
        degenerate_blocks(spectrum),
        exact=True,
    )
    certificate = certs[(4,)]
    assert certificate == {
        "lower_bound": 4,
        "rank": 21,
        "monomial_rank": 6,
        "q": [2],
        "shape": [45, 45],
        "arithmetic": "exact",
    }


def test_greedy_gauge_drop_recovers_a_known_reduction():
    """POWER TEST. A rank-1 class MUST collapse to a single determinant.

    Three determinants sharing a profile class whose flattening has rank 1 are a
    rank-1 tensor, so some element of U(4) on the second block collapses them to
    one. A search that cannot find this has no power and its null results mean
    nothing; this test is what licenses reading a no-drop result as evidence.
    """
    d, blocks = 6, [[0, 1], [2, 3, 4, 5]]
    dets = [(0, 2), (0, 3), (0, 4)]
    c = np.array([0.6, 0.6, np.sqrt(0.28)], dtype=complex)
    c /= np.linalg.norm(c)
    assert flattening_lower_bound(c, dets, blocks) == 1
    out, new_dets, rounds = greedy_gauge_drop(c, dets, blocks, d, tries=4)
    assert len(new_dets) == 1
    assert abs(abs(out[0]) - 1.0) < 1e-9
    assert all(r["closure_residual"] < 1e-12 for r in rounds)


def test_greedy_gauge_drop_respects_a_provable_obstruction():
    """CONTROL. A rank-2 flattening cannot drop below 2, and the search agrees."""
    d, blocks = 6, [[0, 1], [2, 3, 4, 5]]
    dets = [(0, 2), (1, 3)]
    c = np.array([0.6, 0.8], dtype=complex)
    assert flattening_lower_bound(c, dets, blocks) == 2
    _, new_dets, _ = greedy_gauge_drop(c, dets, blocks, d, tries=4)
    assert len(new_dets) == 2


def test_descent_does_not_walk_away_from_a_sparse_point():
    """Regression: a smoothed surrogate with a large epsilon rewards spreading.

    With the smoothing set comparable to the reweighting epsilon, growing a zero
    amplitude is free up to that epsilon and the descent leaves the sparsest
    point it was handed. It must not.
    """
    c, dets, lam = _v103()
    blocks = degenerate_blocks(lam)
    orbit = GaugeOrbit(c, dets, blocks, len(lam))
    _, out = descend(orbit, u0=None)
    assert int((np.abs(out) > 1e-7).sum()) <= len(dets)


def test_gradient_matches_finite_differences():
    """The analytic Riemannian gradient is the whole reason the search is cheap."""
    c, dets, lam = _v103()
    blocks = degenerate_blocks(lam)
    d = len(lam)
    orbit = GaugeOrbit(c, dets, blocks, d)
    rng = np.random.default_rng(4)
    u = random_block_unitary(blocks, d, rng)
    w = 1.0 / (np.abs(orbit.amplitudes(u)) + 1e-2)
    eps = 1e-2
    _, g, _ = orbit.value_grad(u, w, eps)
    h = 1e-6
    for a in (0, 5, len(orbit.gens) - 1):
        step = np.zeros(len(orbit.gens))
        step[a] = 1.0
        plus = orbit.value(orbit.step_unitary(u, step, h), w, eps)
        minus = orbit.value(orbit.step_unitary(u, step, -h), w, eps)
        assert abs((plus - minus) / (2 * h) - g[a]) < 1e-5


def test_v103_cascade_endpoint_is_the_library_state_rotated():
    """The endpoint is the SAME state in another basis, and the basis is not a gauge.

    Two claims at once, both load bearing for the corrected record: the shipped
    support-10 cascade endpoint lies on the U(d) orbit of the certified library
    state (so it is not a new attainer), and the connecting rotation is NOT
    block-diagonal, so the endpoint is not a natural-orbital representative and
    does not bound s_Q^NO.
    """
    c, dets, lam = _v103()
    d = len(lam)
    end = None
    for line in (DATA / "cascade.jsonl").read_text().splitlines():
        r = json.loads(line)
        if r["system"] == "(3,10)" and r["index"] == 103:
            end = r
            break
    ce = np.array(end["final_amplitudes_real"]) + 1j * np.array(end["final_amplitudes_imag"])
    de = [tuple(t) for t in end["final_support_dets"]]
    assert len(de) == 10
    same, dev = same_state_up_to_orbital_rotation(c, dets, ce, de, d)
    assert same and dev < 1e-12
    rho = one_rdm(ce, de, d)
    offdiag = np.max(np.abs(rho - np.diag(np.diag(rho))))
    assert offdiag > 1e-3  # NOT diagonal: 5/68
    u, resid = solve_orbital_rotation(c, dets, ce, de, d, blocks=None, tries=12)
    assert resid < 1e-10
    blocks = degenerate_blocks(lam)
    off_block = max(np.max(np.abs(u[np.ix_(blocks[0], blocks[1])])),
                    np.max(np.abs(u[np.ix_(blocks[1], blocks[0])])))
    assert off_block > 1e-3  # the rotation leaves the natural-orbital gauge


def test_acceptance_gate_accepts_the_library_state():
    """The gate must pass a genuine natural-orbital representative."""
    from gpc_census.validate import check_vertex_attainer

    c, dets, lam = _v103()
    assert check_vertex_attainer(c, dets, lam, len(lam)) == []


def test_acceptance_gate_rejects_a_non_diagonal_attainer():
    """The v103 support-10 state attains the spectrum but fails the gate.

    Its eigenvalues are lambda, so it is a legitimate attainer; its 1-RDM is not
    diagonal, so it is not a natural-orbital representative and the gate must
    say so rather than let its support be quoted against s_I.
    """
    from gpc_census.validate import check_vertex_attainer

    c, dets, lam = _v103()
    end = None
    for line in (DATA / "cascade.jsonl").read_text().splitlines():
        r = json.loads(line)
        if r["system"] == "(3,10)" and r["index"] == 103:
            end = r
            break
    ce = np.array(end["final_amplitudes_real"]) + 1j * np.array(end["final_amplitudes_imag"])
    de = [tuple(t) for t in end["final_support_dets"]]
    errs = check_vertex_attainer(ce, de, lam, len(lam))
    assert any("not diagonal" in e for e in errs)
    # and claiming it as a NEW attainer must fail on the 2-RDM spectrum too
    errs_new = check_vertex_attainer(ce, de, lam, len(lam), library=(c, dets),
                                     claimed_new=True)
    assert any("SAME state in different orbitals" in e for e in errs_new)
