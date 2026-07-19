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


def admissible_support(n: int, d: int, spectrum, groups="degenerate"):
    """Selection rule at a vertex, in the state's natural basis: every
    active facet forces support on determinants achieving its bound.
    groups=None applies the strict rule (valid when the 1-RDM is diagonal
    in the canonical mode basis). Otherwise the strict set is closed under
    occupancy signatures over the given mode groups: "degenerate" uses the
    lambda-degenerate classes; an explicit list of mode collections covers
    rotated natural bases whose 2x2 blocks mix classes (see
    docs/RESEARCH.md, degeneracy lemma)."""
    from fractions import Fraction as F

    from .constraints import constraints

    spectrum = [F(x) for x in spectrum]
    sys_ = constraints(n, d)
    dets = list(combinations(range(d), n))
    active = []
    for iq in sys_["inequalities"] + sys_["equalities"]:
        a, b = iq["coeffs"], iq["rhs"]
        if sum(F(c) * x for c, x in zip(a, spectrum)) == b:
            active.append((a, b))
    strict = [t for t in dets
              if all(sum(a[i] for i in t) == b for a, b in active)]
    if groups is None:
        return strict
    if groups == "degenerate":
        groups = _degenerate_blocks(spectrum, d)
    gsets = [set(g) for g in groups]

    def sig(t):
        return tuple(sum(1 for m in t if m in g) for g in gsets)

    good = {sig(t) for t in strict}
    return [t for t in dets if sig(t) in good]


def _degenerate_blocks(spectrum, d: int):
    """Contiguous runs of equal spectrum values (vertex spectra are sorted)."""
    from fractions import Fraction as F

    spec = [F(x) for x in spectrum]
    blocks, start = [], 0
    for i in range(1, d + 1):
        if i == d or spec[i] != spec[start]:
            blocks.append(list(range(start, i)))
            start = i
    return blocks


def _ansatz_maps(nv, blocks, d: int, cap: int = 40320):
    """Mode permutations preserving the degree vector and the block
    structure: permutations within equal-nv classes of unblocked modes,
    composed with swaps of identical blocks (same diagonal split and
    off-diagonal). Returns index lists, or None when the group exceeds cap
    (fall back to no dedup)."""
    from itertools import permutations, product
    from math import factorial

    bmodes = set()
    for u, v, *_ in blocks:
        bmodes.add(u)
        bmodes.add(v)
    cls: dict = {}
    for m_ in range(d):
        if m_ not in bmodes:
            cls.setdefault(nv[m_], []).append(m_)
    bsig: dict = {}
    for i, (u, v, a_, b_, x2) in enumerate(blocks):
        bsig.setdefault((a_, b_, x2), []).append(i)
    size = 1
    for ms in cls.values():
        size *= factorial(len(ms))
    for grp in bsig.values():
        size *= factorial(len(grp))
    if size > cap:
        return None
    maps = []
    swap_groups = list(bsig.values())
    for combo in product(*(list(permutations(ms)) for ms in cls.values())):
        pm0 = list(range(d))
        for orig, img in zip(cls.values(), combo):
            for o, i2 in zip(orig, img):
                pm0[o] = i2
        for sw in product(*(list(permutations(g)) for g in swap_groups)):
            pm = list(pm0)
            for grp, img in zip(swap_groups, sw):
                for gi, gj in zip(grp, img):
                    pm[blocks[gi][0]] = blocks[gj][0]
                    pm[blocks[gi][1]] = blocks[gj][1]
            maps.append(pm)
    return maps


def _orbit_key(items, maps):
    """Canonical representative of a weighted support under block permutations.
    items: sorted tuple of (det, weight)."""
    if maps is None:
        return items
    return min(
        tuple(sorted((tuple(sorted(pm[x] for x in t)), kk) for t, kk in items))
        for pm in maps
    )


def _cancellation_feasible(items, xtargets=None) -> bool:
    """Necessary condition for phase attainment of a block target: each
    off-diagonal 1-RDM entry is a sum of hop terms with magnitudes
    sqrt(k_i k_j) and free phases, which can reach a prescribed resultant x
    (x = 0 off the blocks) only if the polygon inequality holds: the largest
    of the magnitudes and x is at most the sum of the rest. A lone term can
    never cancel, and a target pair with no hop terms is unreachable."""
    items = list(items)
    req = {tuple(sorted(p)): x for p, x in (xtargets or {}).items()}
    terms: dict = {}
    for ai in range(len(items)):
        t1, k1 = items[ai]
        s1 = set(t1)
        for bi in range(ai + 1, len(items)):
            t2, k2 = items[bi]
            diff = s1 ^ set(t2)
            if len(diff) == 2:
                terms.setdefault(tuple(sorted(diff)), []).append((k1 * k2) ** 0.5)
    for pair, mags in terms.items():
        sides = mags + ([req.pop(pair)] if pair in req else [])
        mx = max(sides)
        if mx - (sum(sides) - mx) > 1e-9:
            return False
    return all(x <= 1e-9 for x in req.values())


def _splits(e1: int, e2: int, den: int):
    """Integer diagonal splits (a, b, x2) of a 2x2 Hermitian block with
    exact eigenvalue numerators e1 > e2: a + b = e1 + e2, off-diagonal
    sqrt(x2)/den with x2 = a*b - e1*e2 > 0, occupations within [0, 1]."""
    out = []
    for a_ in range(1, (e1 + e2) // 2 + 1):
        b_ = e1 + e2 - a_
        if b_ > den:
            continue
        x2 = a_ * b_ - e1 * e2
        if x2 > 0:
            out.append((a_, b_, x2))
    return out


def block_ansatze(n: int, d: int, spectrum, max_blocks: int = 2):
    """Discrete family of exact 1-RDM targets at a vertex, canonical
    diagonal first, then block-diagonal targets whose 2x2 blocks mix two
    distinct spectrum values with an integer diagonal split (the general
    form of the historical v_B double-block ansatz; each block has exact
    eigenvalues e1/den, e2/den by construction). Yields (nv, blocks): the
    integer degree vector in canonical mode labeling, and a list of
    (u, v, a, b, x2) with block modes drawn from the tail of each value
    class. Assigning a to the higher class is a gauge choice; the swapped
    variant is related by the mode transposition (u v), which preserves
    the attained spectrum."""
    import math
    from fractions import Fraction as F
    from itertools import combinations_with_replacement, product

    spec = [F(x) for x in spectrum]
    den = 1
    for x in spec:
        den = den * x.denominator // math.gcd(den, x.denominator)
    nv0 = [int(x * den) for x in spec]
    yield nv0, []
    classes: dict = {}
    for m_, v in enumerate(nv0):
        classes.setdefault(v, []).append(m_)
    vals = sorted(classes, reverse=True)
    ptypes = [(v1, v2) for i, v1 in enumerate(vals) for v2 in vals[i + 1:]]
    for nb in range(1, max_blocks + 1):
        for chosen in combinations_with_replacement(range(len(ptypes)), nb):
            need: dict = {}
            for pi in chosen:
                for v in ptypes[pi]:
                    need[v] = need.get(v, 0) + 1
            if any(cnt > len(classes[v]) for v, cnt in need.items()):
                continue
            cursor = {v: len(classes[v]) - cnt for v, cnt in need.items()}
            assign = []
            for pi in chosen:
                e1, e2 = ptypes[pi]
                u = classes[e1][cursor[e1]]
                cursor[e1] += 1
                v_ = classes[e2][cursor[e2]]
                cursor[e2] += 1
                assign.append((u, v_, e1, e2))
            split_lists = [_splits(e1, e2, den) for (_, _, e1, e2) in assign]
            if any(not sl for sl in split_lists):
                continue
            for si in product(*(range(len(sl)) for sl in split_lists)):
                if any(chosen[i] == chosen[i - 1] and si[i] < si[i - 1]
                       for i in range(1, nb)):
                    continue  # mirror dedup for identical pair types
                nv = list(nv0)
                blocks = []
                for (u, v_, e1, e2), sl, j in zip(assign, split_lists, si):
                    a_, b_, x2 = sl[j]
                    nv[u], nv[v_] = a_, b_
                    blocks.append((u, v_, a_, b_, x2))
                yield nv, blocks


def _skeleton_model(n: int, d: int, spectrum, nv=None, support_filter=None,
                    require_hop_pairs=None):
    """CP-SAT model of the degree system: weights k on allowed determinants
    with prescribed mode sums, indicators y, and (for block ansatze) at least
    one support hop across each required mode pair. Returns
    (model, k, y, dets, den, allowed) or None when a required hop pair has no
    candidate hops inside the allowed support."""
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
    if nv is None:
        nv = [int(x * den) for x in spectrum]
    dets = list(combinations(range(d), n))
    allowed = list(support_filter) if support_filter is not None else dets
    rows = [[j for j, t in enumerate(allowed) if m in t] for m in range(d)]
    m = cp_model.CpModel()
    k = [m.NewIntVar(0, den, f"k{j}") for j in range(len(allowed))]
    y = [m.NewBoolVar(f"y{j}") for j in range(len(allowed))]
    for j in range(len(allowed)):
        m.Add(k[j] <= den * y[j])
        m.Add(k[j] >= y[j])
    for mo, nm in enumerate(nv):
        m.Add(sum(k[j] for j in rows[mo]) == nm)
    for u, v in (require_hop_pairs or []):
        pair = {u, v}
        hops = []
        for p in range(len(allowed)):
            sp_ = set(allowed[p])
            for q in range(p + 1, len(allowed)):
                if sp_ ^ set(allowed[q]) == pair:
                    z = m.NewBoolVar(f"z{u}_{v}_{p}_{q}")
                    # z == y[p] AND y[q]: fully determined, so solution
                    # enumeration does not multiply over free z assignments
                    m.Add(z <= y[p])
                    m.Add(z <= y[q])
                    m.Add(z >= y[p] + y[q] - 1)
                    hops.append(z)
        if not hops:
            return None
        m.Add(sum(hops) >= 1)
    return m, k, y, dets, den, allowed


def min_support_cardinality(n: int, d: int, spectrum, nv=None, support_filter=None,
                            require_hop_pairs=None, time_cap=30):
    """Exact minimum support size of the degree system, or None when the
    system is infeasible (the ansatz admits no skeleton at all). One solve
    replaces a ladder of per-cardinality infeasibility proofs."""
    from ortools.sat.python import cp_model

    built = _skeleton_model(n, d, spectrum, nv=nv, support_filter=support_filter,
                            require_hop_pairs=require_hop_pairs)
    if built is None:
        return None
    m, k, y, dets, den, allowed = built
    m.Minimize(sum(y))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_cap
    st = s.Solve(m)
    if st == cp_model.OPTIMAL:
        return int(sum(s.Value(x) for x in y))
    if st == cp_model.INFEASIBLE:
        return None
    return 2  # unknown within the cap: stay conservative


def enumerate_weight_vectors(n: int, d: int, spectrum, cardinality: int, limit: int = 40,
                             support_filter=None, prefilter=None, nv=None,
                             dedup_maps="auto", require_hop_pairs=None):
    """Integer solutions of the degree system with exactly `cardinality` support
    determinants, one-hop pairs allowed (phases will handle cancellation).
    Solutions are deduplicated to one representative per orbit of the mode
    permutations preserving the degree system (which preserve phase
    solvability), so `limit` counts genuinely distinct skeletons. An optional
    `prefilter(items)` drops skeletons before they count, where items is the
    sorted tuple of (det, weight) pairs. `nv` overrides the degree vector
    (block ansatze), `dedup_maps` the permutation group; `require_hop_pairs`
    forces at least one support hop across each given mode pair."""
    from ortools.sat.python import cp_model

    import math
    from fractions import Fraction as F

    den = 1
    for x in [F(x) for x in spectrum]:
        den = den * x.denominator // math.gcd(den, x.denominator)
    built = _skeleton_model(n, d, spectrum, nv=nv, support_filter=support_filter,
                            require_hop_pairs=require_hop_pairs)
    if built is None:
        return list(combinations(range(d), n)), den, []
    m, k, y, dets, den, allowed = built
    didx = {t: j for j, t in enumerate(dets)}
    if dedup_maps == "auto":
        nv_eff = nv if nv is not None else [int(F(x) * den) for x in spectrum]
        maps = _ansatz_maps(nv_eff, [], d)
    else:
        maps = dedup_maps
    m.Add(sum(y) == cardinality)
    out = []
    seen: set = set()

    class _C(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            items = tuple(sorted(
                (allowed[j], self.Value(kj)) for j, kj in enumerate(k) if self.Value(kj)))
            if prefilter is not None and not prefilter(items):
                return
            key = _orbit_key(items, maps)
            if key in seen:
                return
            seen.add(key)
            w = [0] * len(dets)
            for t, kk in items:
                w[didx[t]] = kk
            out.append(w)
            if len(out) >= limit:
                self.StopSearch()

    s = cp_model.CpSolver()
    s.parameters.enumerate_all_solutions = True
    s.parameters.max_time_in_seconds = 60
    s.Solve(m, _C())
    return dets, den, out


def phase_solve(n: int, d: int, spectrum, dets, den, weights, tries=6, _built=None,
                target=None):
    """Fix moduli sqrt(k/den) exactly; optimize phases only.

    The degree system pins diag(rho) to the target diagonal mode by mode,
    so attaining the spectrum is equivalent to hitting the target 1-RDM
    exactly: for the canonical diagonal target this is Ky Fan equality
    (block diagonality across degenerate blocks; a Hermitian block with
    all eigenvalues equal is scalar), and block targets carry the exact
    spectrum by construction. The objective |rho - target|_F^2 is a smooth
    quartic in the phases with analytic gradients, no eigendecomposition
    in the loop, immune to degeneracy flatness. The 1-RDM map is
    restricted to the support once, so each evaluation is O(d^2 s^2) for
    support size s. Returns (psi, spectral residual)."""
    import numpy as np
    from scipy.optimize import minimize

    built = _built if _built is not None else _build(d, n)
    _, a = built
    dim = a.shape[2]
    sup = [i for i, k in enumerate(weights) if k > 0]
    s_ = len(sup)
    mod = np.array([(weights[i] / den) ** 0.5 for i in sup])
    asub = a[:, :, sup][:, :, :, sup]  # (d, d, s, s)
    if target is None:
        target = np.diag([float(x) for x in spectrum]).astype(complex)
    lam = np.sort(np.array([float(x) for x in spectrum]))[::-1]

    def fg(th):
        c = mod * np.exp(1j * np.concatenate(([0.0], th)))  # gauge: th[sup0] = 0
        rho = np.einsum("p,mnpq,q->mn", c.conj(), asub, c)
        m_ = rho - target
        en = float(np.real(np.sum(m_ * m_.conj())))
        g1 = np.einsum("mn,mnpq,q->p", m_.conj(), asub, c)
        grad = 4.0 * np.imag(c.conj() * g1)
        return en, grad[1:]

    best = (None, 1e9)
    for _ in range(tries):
        th0 = np.random.uniform(0, 2 * np.pi, s_ - 1)
        r = minimize(fg, th0, jac=True, method="L-BFGS-B",
                     options={"maxiter": 1000, "ftol": 1e-22, "gtol": 1e-16})
        if r.fun < best[1]:
            best = (r.x, r.fun)
        if best[1] < 1e-22:
            break
    if best[0] is None:
        return None, 1e9
    c = mod * np.exp(1j * np.concatenate(([0.0], best[0])))
    rho = np.einsum("p,mnpq,q->mn", c.conj(), asub, c)
    e = np.linalg.eigvalsh(rho)[::-1].real
    psi = np.zeros(dim, complex)
    psi[sup] = c
    return psi, float(np.sum((e - lam) ** 2))


def solve_vertex_exact_first(n: int, d: int, spectrum, max_card: int = 24,
                             max_blocks: int = 2, _built=None):
    """Weights-first solve: enumerate integer weight skeletons by ascending
    support size and phase-solve each; moduli are on the natural grid by
    construction, so Tier B needs only phase recognition. Sweeps the block
    ansatz family (canonical diagonal target first, then 2x2 block targets
    that rotate the natural basis, the degeneracy lemma) with cardinality
    as the outer loop, so cheap sparse skeletons across all targets are
    tried before deep ones. Support filters: the strict selection rule for
    the diagonal target, signature closure over merged value classes for
    block targets. Returned support is in the ansatz mode labeling, which
    matches the canonical one up to a spectrum-preserving permutation."""
    import math
    from fractions import Fraction as F

    import numpy as np

    built = _built if _built is not None else _build(d, n)
    dets_all = built[0]
    spec = [F(x) for x in spectrum]
    den = 1
    for x in spec:
        den = den * x.denominator // math.gcd(den, x.denominator)
    nv0 = [int(x * den) for x in spec]
    classes: dict = {}
    for m_, v in enumerate(nv0):
        classes.setdefault(v, []).append(m_)

    ansatze = []
    for nv, blocks in block_ansatze(n, d, spectrum, max_blocks=max_blocks):
        if not blocks:
            adm = admissible_support(n, d, spectrum, groups=None)
            target = None
            xt = None
        else:
            merged = {v: frozenset([v]) for v in classes}
            for u, v_, a_, b_, x2 in blocks:
                key = merged[nv0[u]] | merged[nv0[v_]]
                for v in key:
                    merged[v] = key
            groups = [[m_ for m_ in range(d) if nv0[m_] in g]
                      for g in set(merged.values())]
            adm = admissible_support(n, d, spectrum, groups=groups)
            target = np.diag([v / den for v in nv]).astype(complex)
            for u, v_, a_, b_, x2 in blocks:
                target[u, v_] = target[v_, u] = (x2 ** 0.5) / den
            xt = {(u, v_): x2 ** 0.5 for u, v_, a_, b_, x2 in blocks}
        maps = _ansatz_maps(nv, blocks, d)
        hop_pairs = [(u, v_) for u, v_, a_, b_, x2 in blocks]

        def pref(items, _xt=xt):
            return _cancellation_feasible(items, _xt)

        mincard = min_support_cardinality(n, d, spectrum, nv=nv,
                                          support_filter=adm,
                                          require_hop_pairs=hop_pairs)
        if mincard is None:
            continue  # ansatz admits no skeleton at all
        ansatze.append((nv, blocks, adm, target, maps, pref, hop_pairs, mincard))

    for card in range(2, max_card + 1):
        for nv, blocks, adm, target, maps, pref, hop_pairs, mincard in ansatze:
            if card > len(adm) or card < mincard:
                continue
            dets, _den, sols = enumerate_weight_vectors(
                n, d, spectrum, card, support_filter=adm, prefilter=pref,
                nv=nv, dedup_maps=maps, require_hop_pairs=hop_pairs)
            for w in sols:
                psi, res = phase_solve(n, d, spectrum, dets, den, w,
                                       _built=built, target=target)
                if psi is None:
                    continue
                if 1e-12 <= res < 1e-9:
                    mask = np.array([1.0 if k > 0 else 0.0 for k in w])
                    psi2, res2, _ = attain(n, d, spectrum, mask=mask, outer=80,
                                           tries=1, _built=built, psi0=psi)
                    if res2 < res:
                        psi, res = psi2, res2
                if res < 1e-12:
                    sup = [i for i, k in enumerate(w) if k > 0]
                    j0 = max(sup, key=lambda i: abs(psi[i]))
                    psi = psi * np.exp(-1j * np.angle(psi[j0]))
                    return {"status": "OK", "residual": res,
                            "support_size": len(sup),
                            "weights": [w[i] for i in sup], "den": den,
                            "ansatz": {"nv": nv, "blocks": [list(b) for b in blocks]},
                            "support": [[list(dets_all[i]), float(abs(psi[i])),
                                         float(np.angle(psi[i]))] for i in sup]}
    return {"status": "FAIL", "reason": f"no skeleton through cardinality {max_card}"}
