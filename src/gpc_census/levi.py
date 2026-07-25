"""Levi-orbit diagnostics.  Corrected.

FIX 1 (critical): `levi_theorem_consistent` was `slater == slater`, a tautology.
    The gate could never fire.  Now the check is whether a one-dimensional
    branching sector's INDUCED occupation pattern actually reproduces the target
    spectrum -- a falsifiable statement.

FIX 2: `transverse_tangent_dim` is dominated by the ambient count 2H - d^2 and
    carries no structural signal.  The submersion RANK DEFICIENCY d^2 - rank(J)
    does: it is where the moment map psi -> rho degenerates.

FIX 3: fixed-rho fiber renamed; it is NOT the census's fixed-spectrum fiber.
FIX 4: support determinants sorted with fermionic sign tracking.
"""
from __future__ import annotations

import itertools as it
from itertools import combinations
from math import comb

import numpy as np
import sympy as sp


def spectrum_multiplicities(integer_form):
    counts, order = {}, []
    for v in integer_form:
        v = int(v)
        if v not in counts:
            counts[v] = 0
            order.append(v)
        counts[v] += 1
    return [counts[v] for v in order], order


def levi_dimension(m):
    return sum(int(x) ** 2 for x in m)


def block_sector_dimensions(n, mults):
    mults = tuple(int(m) for m in mults)
    out = []

    def rec(k, left, occ):
        if k == len(mults):
            if left == 0:
                dim = 1
                for m, nk in zip(mults, occ):
                    dim *= comb(m, nk)
                out.append({"block_occupations": list(occ), "dimension": dim,
                            "projective_line_possible": dim == 1})
            return
        for nk in range(max(0, left - sum(mults[k + 1:])), min(mults[k], left) + 1):
            rec(k + 1, left - nk, occ + [nk])

    rec(0, int(n), [])
    return out


def full_levi_line_possible(n, mults):
    return any(r["projective_line_possible"] for r in block_sector_dimensions(n, mults))


def is_slater_spectrum(integer_form, denominator):
    return all(int(x) in (0, int(denominator)) for x in integer_form)


def static_levi_record(record):
    """FIX 1: falsifiable theorem check."""
    n, d = (int(x) for x in record["system"].strip("()").split(","))
    lam = [int(x) for x in record["integer_form"]]
    den = int(record["denominator"])
    mults, values = spectrum_multiplicities(lam)
    slater = is_slater_spectrum(lam, den)

    sectors = block_sector_dimensions(n, mults)
    lines = [s for s in sectors if s["dimension"] == 1]
    # a 1-dim sector fills each block completely or leaves it empty
    # -> induced eigenvalue is den (full) or 0 (empty)
    induced = {tuple(den if nk == m else 0 for nk, m in zip(s["block_occupations"], mults))
               for s in lines}
    target = tuple(values)
    target_line_possible = target in induced          # NON-VACUOUS

    # hard obstruction (permutation level): no linear character of prod S_{m_k}
    hard = not any(sum(p) == n for p in
                   it.product(*[sorted({0, 1, m - 1, m}) for m in mults]))

    return {"system": record["system"], "index": int(record["index"]), "n": n, "d": d,
            "integer_form": lam, "denominator": den, "multiplicities": mults,
            "levi_dim": levi_dimension(mults), "slater_spectrum": slater,
            "representation_contains_full_levi_line": full_levi_line_possible(n, mults),
            "target_full_levi_line_possible": target_line_possible,
            "levi_theorem_consistent": target_line_possible == slater,
            "hard_obstruction": hard}


def determinant_basis(n, d):
    return list(combinations(range(d), n))


def _sort_det_with_sign(det):
    """FIX 4: sort a determinant index tuple, returning (sorted, parity sign)."""
    det = list(det)
    sign = 1
    for i in range(len(det)):
        for j in range(i + 1, len(det)):
            if det[i] > det[j]:
                det[i], det[j] = det[j], det[i]
                sign = -sign
    return tuple(det), sign


def state_vector(record):
    n, d = (int(x) for x in record["system"].strip("()").split(","))
    cf = record["closed_form"]
    basis = determinant_basis(n, d)
    pos = {b: i for i, b in enumerate(basis)}
    amps = [complex(sp.N(sp.sympify(x, locals={"I": sp.I}), 40)) for x in cf["pretty"]]
    psi = np.zeros(len(basis), dtype=complex)
    for det, a in zip(cf["support_dets"], amps):
        s, sgn = _sort_det_with_sign([int(x) for x in det])
        psi[pos[s]] += sgn * a
    psi /= np.linalg.norm(psi)
    return psi, basis, n, d


def annihilation_matrices(n, d, basis):
    lower = determinant_basis(n - 1, d)
    lp = {x: i for i, x in enumerate(lower)}
    mats = [np.zeros((len(lower), len(basis)), dtype=complex) for _ in range(d)]
    for col, det in enumerate(basis):
        for q in det:
            red = tuple(x for x in det if x != q)
            mats[q][lp[red], col] = (-1) ** det.index(q)
    return mats


def one_rdm(psi, ann):
    low = [a @ psi for a in ann]
    d = len(ann)
    return np.array([[np.vdot(low[p], low[q]) for q in range(d)] for p in range(d)])


def _herm_coords_matrix(d):
    """index map for d^2 real coords of a Hermitian dxd matrix"""
    idx = []
    for i in range(d):
        idx.append(("re", i, i))
    for i in range(d):
        for j in range(i + 1, d):
            idx.append(("re", i, j))
            idx.append(("im", i, j))
    return idx


def fiber_and_orbit(psi, basis, n, d, rank_rtol=1e-9, eig_tol=1e-8):
    """Levi-orbit diagnostics plus the FIXED-RHO fiber dimension.

    vectorized: build J from u_pq = A_p^dag A_q psi  (d^2 matvecs, not 2H*d)

    The dimension returned as fixed_rho_fiber_dim is of the fiber
    {psi : rho_1(psi) = rho}, over the FULL Fock space of (n, d), with only the
    global phase quotiented. It is NOT the census fixed-spectrum fiber and not
    comparable to gpc_census.fiber.fixed_spectrum_tangent, which is
    support-restricted and quotients the whole orbital-phase gauge. See the
    module docstring of gpc_census.fiber for the two-fiber contract.
    """
    ann = annihilation_matrices(n, d, basis)
    low = [a @ psi for a in ann]
    rho = np.array([[np.vdot(low[p], low[q]) for q in range(d)] for p in range(d)])
    H = len(basis)
    U = [[ann[p].conj().T @ low[q] for q in range(d)] for p in range(d)]
    coords = _herm_coords_matrix(d)
    J = np.zeros((len(coords), 2 * H))
    for r, (kind, p, q) in enumerate(coords):
        # drho[p,q] = x.(u_pq + conj(u_qp)) + y.(-i u_pq + i conj(u_qp))
        vx = U[p][q] + np.conj(U[q][p])
        vy = -1j * U[p][q] + 1j * np.conj(U[q][p])
        if kind == "re":
            J[r, :H] = vx.real
            J[r, H:] = vy.real
        else:
            J[r, :H] = vx.imag
            J[r, H:] = vy.imag
    s = np.linalg.svd(J, compute_uv=False)
    rank = int((s > rank_rtol * s[0]).sum()) if s.size else 0
    fixed_rho_fiber_dim = max(0, 2 * H - rank - 1)     # quotient global phase
    rank_deficiency = d * d - rank                      # FIX 2: the real signal

    # Levi orbit
    w, V = np.linalg.eigh(rho)
    o = np.argsort(w)[::-1]
    w = w[o]
    V = V[:, o]
    # one_rdm returns <a_p^dag a_q> = rho_1^T, whose eigenvectors are the
    # CONJUGATES of the natural orbitals.  Build Levi generators from the
    # natural orbitals themselves; without this, stabilizer/orbit ranks are
    # wrong for any state with complex amplitudes (real states unaffected).
    V = V.conj()
    blocks, cur = [], [0]
    for i in range(1, d):
        if abs(w[i] - w[cur[0]]) > eig_tol:
            blocks.append(cur)
            cur = [i]
        else:
            cur.append(i)
    blocks.append(cur)
    gens = []
    for blk in blocks:
        for ia, a in enumerate(blk):
            va = V[:, a]
            gens.append(np.outer(va, va.conj()))
            for b in blk[ia + 1:]:
                vb = V[:, b]
                E = np.outer(va, vb.conj())
                F = np.outer(vb, va.conj())
                gens.append(E + F)
                gens.append(-1j * (E - F))
    cols = []
    pos = {b: i for i, b in enumerate(basis)}
    for Hm in gens:
        out = np.zeros(H, dtype=complex)
        for col, det in enumerate(basis):
            amp = psi[col]
            if amp == 0:
                continue
            for q in det:
                rs = (-1) ** det.index(q)
                red = tuple(x for x in det if x != q)
                for p in range(d):
                    if p in red:
                        continue
                    tg = tuple(sorted(red + (p,)))
                    isg = (-1) ** tg.index(p)
                    out[pos[tg]] += Hm[p, q] * rs * isg * amp
        v = 1j * out
        cols.append(np.concatenate((v.real, v.imag)))
    Om = np.column_stack(cols)
    so = np.linalg.svd(Om, compute_uv=False)
    orbit_rank = int((so > rank_rtol * so[0]).sum()) if so.size else 0
    orbit_dim = max(0, orbit_rank - 1)
    mults = [len(b) for b in blocks]
    return dict(numerical_multiplicities=mults, levi_dim=levi_dimension(mults),
                levi_orbit_dim=orbit_dim,
                projective_stabilizer_dim=levi_dimension(mults) - orbit_dim,
                jacobian_rank=rank, rank_deficiency=rank_deficiency,
                fixed_rho_fiber_dim=fixed_rho_fiber_dim,
                transverse_dim=fixed_rho_fiber_dim - orbit_dim,
                hilbert_dim=H)
