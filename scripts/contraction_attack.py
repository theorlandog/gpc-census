#!/usr/bin/env python3
"""Ansatz-free numerical attainment of a GPC vertex via the contraction map.

Minimize ||gamma(Psi / ||Psi||) - gamma_0||^2 over the FULL Lambda^N C^d tensor
(every determinant amplitude free -- no support, sparsity, rationality, or block
restriction), with an analytic Wirtinger gradient. gamma(Psi / ||Psi||) is the 1-RDM,
assembled from the diagonal occupations and the one-hop (share N-1 orbitals)
off-diagonal contractions; gamma_0 = diag(spectrum). The analytic gradient is what
makes this converge to machine precision in seconds where support-restricted and
weak-gradient searches stall.

This attains the two open (3,10) holdouts v89 = (15,15,6^8)/26 and
v103 = (18^4,5^6)/34 to ~1e-16, complex AND real, from random starts -- see
docs/RESEARCH_ARCHIVE.md. IMPORTANT STATUS: numerical attainment is NOT certification. The
solutions are DENSE generic fiber points (support ~100-120 of 120), so they do not
yield a closed form and do NOT close the census (which requires an exact,
char-poly-gated state). This tool establishes attainability and the real/complex
distinction numerically; extracting a certifiable distinguished representative is a
separate, open step.

Usage: python scripts/contraction_attack.py 15,15,6,6,6,6,6,6,6,6 26 [--real] [--N 3]
"""
from __future__ import annotations

import itertools
import sys

import numpy as np
from scipy.optimize import minimize


def build(d, n):
    """Determinants and the sparse one-hop contraction index arrays for the 1-RDM."""
    dets = list(itertools.combinations(range(d), n))
    idx, aa, bb, sg = [], [], [], []
    for a, ad in enumerate(dets):                      # diagonal occupations
        for i in ad:
            idx.append(i * d + i)
            aa.append(a)
            bb.append(a)
            sg.append(1.0)
    for a, ad in enumerate(dets):                      # one-hop off-diagonals
        sa = set(ad)
        for b, bd in enumerate(dets):
            if a == b:
                continue
            common = sa & set(bd)
            if len(common) == n - 1:
                (i,) = sa - common
                (j,) = set(bd) - common
                idx.append(i * d + j)
                aa.append(a)
                bb.append(b)
                sg.append(float((-1) ** ad.index(i) * (-1) ** bd.index(j)))
    return dets, (np.array(idx), np.array(aa), np.array(bb), np.array(sg))


def normalized_objective_gradient(x, g0, contraction, real=False):
    """Scale-invariant distance and analytic gradient for a normalized state.

    ``contraction`` is the sparse tuple returned by :func:`build`. The input may
    have any nonzero norm. Normalizing inside the objective is load bearing:
    without it, a positive optimum is a distance to the cone of unnormalized
    1-RDMs rather than to the moment polytope.
    """
    idx, aa, bb, sg = contraction
    d2 = len(g0)
    nd = len(x) if real else len(x) // 2
    if real:
        c = np.asarray(x, dtype=float)
        norm2 = float(c @ c)
        if norm2 <= np.finfo(float).tiny:
            return float("inf"), np.zeros_like(c)
        raw = np.zeros(d2)
        np.add.at(raw, idx, sg * c[aa] * c[bb])
        rho = raw / norm2
        error = rho - g0
        grad = np.zeros(nd)
        np.add.at(grad, bb, 2 * error[idx] * sg * c[aa] / norm2)
        np.add.at(grad, aa, 2 * error[idx] * sg * c[bb] / norm2)
        overlap = float(error @ raw)
        grad -= 4 * overlap * c / norm2 ** 2
        return float(error @ error), grad

    c = np.asarray(x[:nd]) + 1j * np.asarray(x[nd:])
    norm2 = float(np.vdot(c, c).real)
    if norm2 <= np.finfo(float).tiny:
        return float("inf"), np.zeros_like(x)
    raw = np.zeros(d2, complex)
    np.add.at(raw, idx, sg * c[aa] * np.conj(c[bb]))
    rho = raw / norm2
    error = rho - g0
    wirtinger = np.zeros(nd, complex)
    np.add.at(
        wirtinger,
        bb,
        2 * np.conj(error[idx]) * sg * c[aa] / norm2,
    )
    overlap = float(np.vdot(error, raw).real)
    wirtinger -= 2 * overlap * c / norm2 ** 2
    grad = np.concatenate([2 * np.real(wirtinger), 2 * np.imag(wirtinger)])
    return float(np.vdot(error, error).real), grad


def attack_state(ints, den, n=3, real=False, tries=25, seed=0):
    """Return the best normalized contraction result, including its state."""
    d = len(ints)
    dets, contraction = build(d, n)
    nd = len(dets)
    g0 = np.diag(np.array(ints, dtype=float) / den).flatten()
    dim = nd if real else 2 * nd

    def objgrad(x):
        return normalized_objective_gradient(x, g0, contraction, real=real)

    best = (1e9, None)
    for k in range(tries):
        rng = np.random.default_rng(seed + k)
        x0 = rng.standard_normal(dim)
        x0 /= np.linalg.norm(x0)
        r = minimize(
            objgrad,
            x0,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 20000, "ftol": 1e-20, "gtol": 1e-16},
        )
        if r.fun < best[0]:
            best = (r.fun, r.x)
        if best[0] < 1e-19:
            break
    x = best[1]
    c = x if real else x[:nd] + 1j * x[nd:]
    c = c / np.linalg.norm(c)
    idx, aa, bb, sg = contraction
    e = np.zeros(d * d, complex)
    np.add.at(e, idx, sg * (c[aa] * c[bb] if real else c[aa] * np.conj(c[bb])))
    rho = e.reshape(d, d)
    ev = np.sort(np.linalg.eigvalsh((rho + rho.conj().T) / 2))[::-1]
    supp = int(np.sum(np.abs(c) ** 2 > 1e-6 * np.max(np.abs(c) ** 2)))
    return {
        "residual": float(best[0]),
        "eigenvalues": ev,
        "support_size": supp,
        "amplitudes": c,
        "determinants": dets,
        "one_rdm": rho,
    }


def attack(ints, den, n=3, real=False, tries=25, seed=0):
    """Return (best_residual, sorted_eigenvalues, support_size) for the attainment."""
    result = attack_state(ints, den, n=n, real=real, tries=tries, seed=seed)
    return result["residual"], result["eigenvalues"], result["support_size"]


if __name__ == "__main__":
    ints = [int(x) for x in sys.argv[1].split(",")]
    den = int(sys.argv[2])
    n = int(sys.argv[sys.argv.index("--N") + 1]) if "--N" in sys.argv else 3
    res, ev, supp = attack(ints, den, n=n, real="--real" in sys.argv)
    print(f"residual = {res:.3e}   support ~ {supp}   eig*den = {np.round(ev * den, 5)}")
    print("NOTE: numerically attained, NOT certified -- dense generic state, census unchanged.")
