"""Extremal state construction by alternating spectral projection.

Attains a target vertex spectrum with a complex pure state, then greedily
minimizes determinant support. Validated by reconstruction of the
(4,9) interference vertex v_B from a random start.
"""
from __future__ import annotations

from itertools import combinations


def _build(d: int, n: int):
    import numpy as np

    dets = list(combinations(range(d), n))
    idx = {t: i for i, t in enumerate(dets)}
    dim = len(dets)
    a = np.zeros((d, d, dim, dim), complex)
    for t in dets:
        for mp in t:
            s1 = (-1) ** t.index(mp)
            t2 = tuple(x for x in t if x != mp)
            for m in range(d):
                if m in t2:
                    continue
                tp = tuple(sorted(t2 + (m,)))
                s2 = (-1) ** tp.index(m)
                a[m, mp, idx[tp], idx[t]] = s1 * s2
    return dets, a


def attain(n: int, d: int, spectrum, mask=None, outer=250, tries=5, tol=1e-16, _built=None):
    """Find a complex pure state whose 1-RDM spectrum equals the target."""
    import numpy as np
    from scipy.optimize import minimize

    dets, a = _built if _built is not None else _build(d, n)
    dim = a.shape[2]
    lam = np.sort(np.array([float(x) for x in spectrum]))[::-1]
    best = (None, 1e9)
    for _ in range(tries):
        psi = np.random.randn(dim) + 1j * np.random.randn(dim)
        if mask is not None:
            psi *= mask
        psi /= np.linalg.norm(psi)
        spec_res = 1e9
        for _o in range(outer):
            rho = np.einsum("i,mnij,j->mn", psi.conj(), a, psi)
            e, u = np.linalg.eigh(rho)
            tgt = u[:, ::-1] @ np.diag(lam) @ u[:, ::-1].conj().T

            def fg(x):
                p = x[:dim] + 1j * x[dim:]
                if mask is not None:
                    p = p * mask
                nn = np.linalg.norm(p)
                p = p / max(nn, 1e-14)
                rho = np.einsum("i,mnij,j->mn", p.conj(), a, p)
                m_ = rho - tgt
                en = float(np.real(np.sum(m_ * m_.conj())))
                b = np.einsum("nm,mnij->ij", m_, a)
                g = b @ p
                g = g - np.vdot(p, g).real * p
                if mask is not None:
                    g = g * mask
                return en, np.concatenate([g.real, g.imag]) * 2 / max(nn, 1e-14)

            x0 = np.concatenate([psi.real, psi.imag])
            res = minimize(fg, x0, jac=True, method="L-BFGS-B",
                           options={"maxiter": 60, "ftol": 1e-18, "gtol": 1e-14})
            psi = res.x[:dim] + 1j * res.x[dim:]
            if mask is not None:
                psi = psi * mask
            psi /= np.linalg.norm(psi)
            rho = np.einsum("i,mnij,j->mn", psi.conj(), a, psi)
            e = np.linalg.eigvalsh(rho)[::-1].real
            spec_res = float(np.sum((e - lam) ** 2))
            if spec_res < tol:
                return psi, spec_res, dets
        if spec_res < best[1]:
            best = (psi, spec_res)
    return best[0], best[1], dets


def minimize_support(n: int, d: int, spectrum, psi, dets, res_tol=1e-12, _built=None):
    """Greedily remove determinants while the spectrum stays attainable."""
    import numpy as np

    built = _built if _built is not None else _build(d, n)
    mask = (np.abs(psi) > 1e-10).astype(float)
    improved = True
    while improved:
        improved = False
        order = sorted(np.where(mask > 0)[0], key=lambda i: abs(psi[i]))
        for i in order:
            m2 = mask.copy()
            m2[i] = 0
            if m2.sum() < 1:
                continue
            cand, res, _ = attain(n, d, spectrum, mask=m2, outer=120, tries=3, _built=built)
            if res < res_tol:
                psi, mask = cand, m2
                improved = True
                break
    sup = [i for i in range(len(psi)) if mask[i] and abs(psi[i]) > 1e-9]
    j0 = max(sup, key=lambda i: abs(psi[i]))
    psi = psi * np.exp(-1j * np.angle(psi[j0]))
    return psi, sup


def solve_vertex(n: int, d: int, spectrum, _built=None):
    """Attain a vertex spectrum and return the minimal-support state record."""
    import numpy as np

    built = _built if _built is not None else _build(d, n)
    dets = built[0]
    psi, res, _ = attain(n, d, spectrum, _built=built)
    if res > 1e-12:
        return {"status": "FAIL", "residual": res}
    psi, sup = minimize_support(n, d, spectrum, psi, dets, _built=built)
    return {"status": "OK", "residual": res, "support_size": len(sup),
            "support": [[list(dets[i]), float(abs(psi[i])), float(np.angle(psi[i]))]
                        for i in sup]}
