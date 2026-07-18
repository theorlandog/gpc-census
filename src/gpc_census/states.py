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


def attain(n: int, d: int, spectrum, mask=None, outer=250, tries=5, tol=1e-16, _built=None, psi0=None):
    """Find a complex pure state whose 1-RDM spectrum equals the target."""
    import numpy as np
    from scipy.optimize import minimize

    dets, a = _built if _built is not None else _build(d, n)
    dim = a.shape[2]
    lam = np.sort(np.array([float(x) for x in spectrum]))[::-1]
    best = (None, 1e9)
    for trial in range(tries):
        if psi0 is not None and trial == 0:
            psi = psi0.copy()
        else:
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
    """Sparsify by iterative hard thresholding, then greedy single removals."""
    import numpy as np

    built = _built if _built is not None else _build(d, n)
    mask = (np.abs(psi) > 1e-10).astype(float)
    # phase 1: threshold ladder with warm starts
    for frac in (0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
        thr = frac * np.max(np.abs(psi))
        m2 = ((np.abs(psi) > thr) & (mask > 0)).astype(float)
        if m2.sum() < 1 or m2.sum() == mask.sum():
            continue
        cand, res, _ = attain(n, d, spectrum, mask=m2, outer=150, tries=2,
                              _built=built, psi0=psi * m2)
        if res < res_tol:
            psi, mask = cand, m2
    # phase 2: greedy single removals (warm started)
    improved = True
    while improved:
        improved = False
        order = sorted(np.where(mask > 0)[0], key=lambda i: abs(psi[i]))
        for i in order:
            m2 = mask.copy()
            m2[i] = 0
            if m2.sum() < 1:
                continue
            cand, res, _ = attain(n, d, spectrum, mask=m2, outer=150, tries=2,
                                  _built=built, psi0=psi * m2)
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


def allowed_dets(n: int, d: int, spectrum):
    """Selection-rule support at a vertex: determinants saturating every
    active GPC facet, containing every fully occupied mode, avoiding every
    empty mode. Pinning of a facet sum(a_i lambda_i) = b forces the state
    onto configurations with sum(a_i, i in t) = b (Klyachko selection rule)."""
    from fractions import Fraction
    from itertools import combinations

    from .constraints import constraints

    spec = [Fraction(x) for x in spectrum]
    sys_ = constraints(n, d)
    active = []
    for iq in sys_["inequalities"] + sys_["equalities"]:
        if sum(c * s for c, s in zip(iq["coeffs"], spec)) == iq["rhs"]:
            active.append((iq["coeffs"], iq["rhs"]))
    full = [i for i, s in enumerate(spec) if s == 1]
    empty = [i for i, s in enumerate(spec) if s == 0]
    out = []
    for t in combinations(range(d), n):
        ts = set(t)
        if any(i not in ts for i in full) or any(i in ts for i in empty):
            continue
        if all(sum(a[i] for i in t) == b for a, b in active):
            out.append(t)
    return out


def attain_diag(n: int, d: int, spectrum, dets_sub, outer=2000, tries=8, tol=1e-24):
    """Attain rho = diag(spectrum) exactly, restricted to a determinant
    subset. Smooth objective, no eigendecomposition."""
    import numpy as np
    from scipy.optimize import minimize

    spec = np.array([float(x) for x in spectrum])
    dim = len(dets_sub)
    idx = {t: i for i, t in enumerate(dets_sub)}
    hops = []
    for t in dets_sub:
        for mp in t:
            s1 = (-1) ** t.index(mp)
            t2 = tuple(x for x in t if x != mp)
            for m in range(d):
                if m in t2:
                    continue
                tp = tuple(sorted(t2 + (m,)))
                if tp not in idx:
                    continue
                s2 = (-1) ** tp.index(m)
                hops.append((idx[t], idx[tp], m, mp, s1 * s2))
    occ = [[i for i, t in enumerate(dets_sub) if m in t] for m in range(d)]

    def fg(x):
        p = x[:dim] + 1j * x[dim:]
        nn = np.linalg.norm(p)
        p = p / max(nn, 1e-14)
        rho = np.zeros((d, d), complex)
        for m in range(d):
            rho[m, m] = np.sum(np.abs(p[occ[m]]) ** 2)
        for j, k, m, mp, s in hops:
            rho[m, mp] += s * np.conj(p[k]) * p[j]
        m_ = rho - np.diag(spec)
        en = float(np.real(np.sum(m_ * np.conj(m_))))
        g = np.zeros(dim, complex)
        for m in range(d):
            g[occ[m]] += m_[m, m].real * p[occ[m]]
        for j, k, m, mp, s in hops:
            g[j] += m_[m, mp] * s * p[k]
        g = g - np.vdot(p, g).real * p
        return en, np.concatenate([g.real, g.imag]) * 2 / max(nn, 1e-14)

    best = (None, 1e9)
    for _ in range(tries):
        p0 = np.random.randn(dim) + 1j * np.random.randn(dim)
        x0 = np.concatenate([p0.real, p0.imag])
        x0 /= np.linalg.norm(x0)
        r = minimize(fg, x0, jac=True, method="L-BFGS-B",
                     options={"maxiter": outer, "ftol": 1e-24, "gtol": 1e-16})
        if r.fun < best[1]:
            best = (r.x, r.fun)
        if best[1] < tol:
            break
    p = best[0][:dim] + 1j * best[0][dim:]
    p /= np.linalg.norm(p)
    j0 = int(np.argmax(np.abs(p)))
    p = p * np.exp(-1j * np.angle(p[j0]))
    return p, best[1]


def admissible_support(n: int, d: int, spectrum):
    """Selection rule at a vertex: pinning to each active facet (a, lam) = b
    forces support on determinants with sum(a_i, i in t) = b. Intersect over
    all active constraints. Exact rational arithmetic, no solver."""
    from fractions import Fraction as F

    from .constraints import constraints

    spectrum = [F(x) for x in spectrum]
    sys_ = constraints(n, d)
    dets = list(combinations(range(d), n))
    keep = list(range(len(dets)))
    for iq in sys_["inequalities"] + sys_["equalities"]:
        a, b = iq["coeffs"], iq["rhs"]
        if sum(F(c) * x for c, x in zip(a, spectrum)) != b:
            continue  # facet not active at this vertex
        keep = [j for j in keep if sum(a[i] for i in dets[j]) == b]
    return [dets[j] for j in keep]


def enumerate_weight_vectors(n: int, d: int, spectrum, cardinality: int, limit: int = 40,
                             support_filter=None):
    """Integer solutions of the degree system with exactly `cardinality` support
    determinants, one-hop pairs allowed (phases will handle cancellation)."""
    import math

    try:
        from ortools.sat.python import cp_model
    except ImportError as e:
        raise RuntimeError(
            "weights-first solve requires the cpsat extra: uv sync --extra cpsat"
        ) from e

    from fractions import Fraction as F
    spectrum = [F(x) for x in spectrum]
    den = 1
    for x in spectrum:
        den = den * x.denominator // math.gcd(den, x.denominator)
    nv = [int(x * den) for x in spectrum]
    dets = list(combinations(range(d), n))
    rows = [[j for j, t in enumerate(dets) if m in t] for m in range(d)]
    m = cp_model.CpModel()
    k = [m.NewIntVar(0, den, f"k{t}") for t in range(len(dets))]
    if support_filter is not None:
        allowed = set(support_filter)
        for j, t in enumerate(dets):
            if t not in allowed:
                m.Add(k[j] == 0)
    y = [m.NewBoolVar(f"y{t}") for t in range(len(dets))]
    for t in range(len(dets)):
        m.Add(k[t] <= den * y[t])
        m.Add(k[t] >= y[t])
    for mo, nm in enumerate(nv):
        m.Add(sum(k[j] for j in rows[mo]) == nm)
    m.Add(sum(y) == cardinality)
    out = []

    class _C(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            out.append([self.Value(x) for x in k])
            if len(out) >= limit:
                self.StopSearch()

    s = cp_model.CpSolver()
    s.parameters.enumerate_all_solutions = True
    s.parameters.max_time_in_seconds = 60
    s.Solve(m, _C())
    return dets, den, out


def phase_solve(n: int, d: int, spectrum, dets, den, weights, tries=6, _built=None):
    """Fix moduli sqrt(k/den) exactly; optimize phases only."""
    import numpy as np
    from scipy.optimize import minimize

    built = _built if _built is not None else _build(d, n)
    _, a = built
    dim = a.shape[2]
    sup = [i for i, k in enumerate(weights) if k > 0]
    mod = np.zeros(dim)
    for i in sup:
        mod[i] = (weights[i] / den) ** 0.5
    lam = np.sort(np.array([float(x) for x in spectrum]))[::-1]
    nfree = len(sup) - 1  # global phase gauge-fixed

    def psi_of(th):
        psi = np.zeros(dim, complex)
        psi[sup[0]] = mod[sup[0]]
        for j, i in enumerate(sup[1:]):
            psi[i] = mod[i] * np.exp(1j * th[j])
        return psi

    def f(th):
        psi = psi_of(th)
        rho = np.einsum("i,mnij,j->mn", psi.conj(), a, psi)
        e = np.linalg.eigvalsh(rho)[::-1].real
        return float(np.sum((e - lam) ** 2))

    best = (None, 1e9)
    for _ in range(tries):
        th0 = np.random.uniform(0, 2 * np.pi, nfree)
        r = minimize(f, th0, method="Nelder-Mead",
                     options={"maxiter": 4000, "xatol": 1e-13, "fatol": 1e-16})
        if r.fun < best[1]:
            best = (r.x, r.fun)
        if best[1] < 1e-16:
            break
    if best[0] is None:
        return None, 1e9
    return psi_of(best[0]), best[1]


def solve_vertex_exact_first(n: int, d: int, spectrum, max_card: int = 24,
                             _built=None):
    """Weights-first solve: enumerate integer weight skeletons by ascending
    support size, phase-solve each; moduli are on the natural grid by
    construction, so Tier B needs only phase recognition."""
    import numpy as np

    built = _built if _built is not None else _build(d, n)
    dets_all = built[0]
    adm = admissible_support(n, d, spectrum)
    for card in range(2, min(max_card, len(adm)) + 1):
        dets, den, sols = enumerate_weight_vectors(n, d, spectrum, card,
                                                   support_filter=adm)
        for w in sols:
            psi, res = phase_solve(n, d, spectrum, dets, den, w, _built=built)
            if res < 1e-15:
                sup = [i for i, k in enumerate(w) if k > 0]
                j0 = max(sup, key=lambda i: abs(psi[i]))
                psi = psi * np.exp(-1j * np.angle(psi[j0]))
                return {"status": "OK", "residual": res, "support_size": len(sup),
                        "weights": [w[i] for i in sup], "den": den,
                        "support": [[list(dets_all[i]), float(abs(psi[i])),
                                     float(np.angle(psi[i]))] for i in sup]}
    return {"status": "FAIL", "reason": f"no skeleton through cardinality {max_card}"}
