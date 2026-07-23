#!/usr/bin/env python3
"""Ansatz-free numerical attainment of a GPC vertex via the contraction map.

Minimize ||gamma(Psi) - gamma_0||^2 over the FULL Lambda^N C^d tensor (every
determinant amplitude free -- no support, sparsity, rationality, or block
restriction), with an analytic Wirtinger gradient. gamma(Psi) is the 1-RDM,
assembled from the diagonal occupations and the one-hop (share N-1 orbitals)
off-diagonal contractions; gamma_0 = diag(spectrum). The analytic gradient is what
makes this converge to machine precision in seconds where support-restricted and
weak-gradient searches stall.

This attains the two open (3,10) holdouts v89 = (15,15,6^8)/26 and
v103 = (18^4,5^6)/34 to ~1e-16, complex AND real, from random starts -- see
docs/RESEARCH.md. IMPORTANT STATUS: numerical attainment is NOT certification. The
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


def attack(ints, den, n=3, real=False, tries=25, seed=0):
    """Return (best_residual, sorted_eigenvalues, support_size) for the attainment."""
    d = len(ints)
    dets, (idx, aa, bb, sg) = build(d, n)
    nd = len(dets)
    g0 = np.diag(np.array(ints) / den).flatten()

    if real:
        def objgrad(x):
            e = np.zeros(d * d)
            np.add.at(e, idx, sg * x[aa] * x[bb])
            e -= g0
            g = np.zeros(nd)
            np.add.at(g, bb, 2 * e[idx] * sg * x[aa])
            np.add.at(g, aa, 2 * e[idx] * sg * x[bb])
            return float(np.sum(e ** 2)), g
        dim = nd
    else:
        def objgrad(x):
            c = x[:nd] + 1j * x[nd:]
            e = np.zeros(d * d, complex)
            np.add.at(e, idx, sg * c[aa] * np.conj(c[bb]))
            e -= g0
            gz = np.zeros(nd, complex)
            np.add.at(gz, bb, 2 * np.conj(e[idx]) * sg * c[aa])
            return float(np.sum(np.abs(e) ** 2)), np.concatenate([2 * np.real(gz), 2 * np.imag(gz)])
        dim = 2 * nd

    best = (1e9, None)
    for k in range(tries):
        rng = np.random.default_rng(seed + k)
        r = minimize(objgrad, rng.standard_normal(dim) * 0.25, jac=True, method="L-BFGS-B",
                     options={"maxiter": 20000, "ftol": 1e-20, "gtol": 1e-16})
        if r.fun < best[0]:
            best = (r.fun, r.x)
        if best[0] < 1e-19:
            break
    x = best[1]
    c = x if real else x[:nd] + 1j * x[nd:]
    e = np.zeros(d * d, complex)
    np.add.at(e, idx, sg * (c[aa] * c[bb] if real else c[aa] * np.conj(c[bb])))
    rho = e.reshape(d, d) / np.sum(np.abs(c) ** 2)
    ev = np.sort(np.linalg.eigvalsh((rho + rho.conj().T) / 2))[::-1]
    supp = int(np.sum(np.abs(c) ** 2 > 1e-6 * np.max(np.abs(c) ** 2)))
    return best[0], ev, supp


if __name__ == "__main__":
    ints = [int(x) for x in sys.argv[1].split(",")]
    den = int(sys.argv[2])
    n = int(sys.argv[sys.argv.index("--N") + 1]) if "--N" in sys.argv else 3
    res, ev, supp = attack(ints, den, n=n, real="--real" in sys.argv)
    print(f"residual = {res:.3e}   support ~ {supp}   eig*den = {np.round(ev * den, 5)}")
    print("NOTE: numerically attained, NOT certified -- dense generic state, census unchanged.")
