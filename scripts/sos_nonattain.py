#!/usr/bin/env python3
"""Certified non-attainability via a charge-0 Hermitian-SOS (moment) relaxation.

This is the RANK-11 CERTIFICATE PATH recommended in docs/RESEARCH.md ("A
Lasserre/SOS moment-relaxation dual of min ||gamma(psi)-gamma_0||^2 yields a
CERTIFIED positive lower bound = rigorous non-attainability"), now runnable: an
SDP solver (Clarabel) is available in the venv. It computes a lower bound delta
on

    min_{||c||^2 = 1}  f(c),   f = sum_ij |gamma_ij(c) - g0_i delta_ij|^2,

the squared distance from the target spectrum diag(g0) to the nearest pure
N-representable 1-RDM. delta > 0 is a rigorous certificate that diag(g0) is
OUTSIDE the (N,d) polytope (non-attainable); the contraction attack only produces
an UNCONSTRAINED floor (evidence, not proof).

FORMULATION (gauge-reduced, exact cone). f is U(1)-gauge invariant
(c -> e^{i phi} c), so every monomial has charge |holo| - |anti| = 0, and f is a
sum of squared magnitudes of the SESQUILINEAR forms gamma_ij - g0 delta -- NOT a
holomorphic SOS. The correct cone squares elements of the charge-0 space
    B = {1} U {Re(c_a cbar_b): a<=b} U {Im(c_a cbar_b): a<b}   (size 1 + nd^2),
nd = C(d,N) the number of Slater determinants. The certificate is
    f - delta = B^T Q B + mu (||c||^2 - 1),   Q >> 0 (real symmetric),
with mu = sum_S mu_S B_S a free charge-0 degree-<=2 multiplier (this is exactly
the level-2 relaxation on the sphere). Coefficients are matched over complex
monomial keys (holo multiset, anti multiset) and handed to Clarabel natively
(min 1/2 x'Px + q'x s.t. Ax + s = b, s in cones): one PSDTriangleCone(1+nd^2) +
a ZeroCone of equalities. delta=0 is always feasible (Q = Gram of f, mu=0), so the
cone is correct; the multiplier lifts delta onto the sphere.

VALIDATED (self-test, run this file): (2,4), whose pure spectra are the degenerate
pairs (a,a,1-a,1-a):
  EXTER  (0.8,0.6,0.4,0.2) -> delta = 0.04  == exact squared distance to the
         nearest attainable (0.7,0.7,0.3,0.3); certified non-attainable.
  ATTAIN (0.7,0.7,0.3,0.3) -> delta ~ 0     (attainable, as it must be).

SCALING WALL (why symmetry reduction is REQUIRED, not optional). The single PSD
block has side 1 + nd^2: 401 for (3,6), 27226 for (3,11). Any dense-scaling
interior-point method forms an npsd x npsd Nesterov-Todd scaling with
npsd = side(side+1)/2, i.e. 80601^2 * 8 bytes ~ 52 GB already at (3,6). So the
UN-REDUCED relaxation is intractable beyond tiny systems. The target cand 44 of
(3,11), g0 = diag(6,6,6,6,6,1,1,1,1,1,1)/12 (numerical floor ~ 2.155e-3 ~ 1/464,
non-attainable -- see scripts/contraction_attack.py), has spectrum stabilizer
S5 x S6; block-diagonalizing Q under that group (Wedderburn / isotypic
decomposition of the fermionic determinant rep) collapses the 27226 side into a
handful of small irrep blocks. That symmetry-adapted solve is the next step; this
module supplies the (validated) un-reduced assembler and the exact target datum.

Run:  .venv/bin/python scripts/sos_nonattain.py            # (2,4) self-test
Diagnostic; scripts/ is ruff-excluded.
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "scripts")
from contraction_attack import build  # noqa: E402

R2 = np.sqrt(2.0)


def gamma_entries(d, N):
    """gamma_ij = sum_k s c_a conj(c_b): dict (i,j) -> {(a,b): signed int}."""
    dets, (IDX, A, B, SG) = build(d, N)
    nd = len(dets)
    G = defaultdict(lambda: defaultdict(int))
    for k in range(len(IDX)):
        i, j = divmod(int(IDX[k]), d)
        G[(i, j)][(int(A[k]), int(B[k]))] += int(round(SG[k]))
    return nd, G


def f_coeffs(d, N, g0):
    """f = sum_ij |gamma_ij - g0 delta|^2 as dict {(holo,anti): complex}."""
    nd, G = gamma_entries(d, N)
    F = defaultdict(complex)
    for i in range(d):
        for j in range(d):
            gij = dict(G.get((i, j), {}))
            L = [(a, b, float(s)) for (a, b), s in gij.items()]
            if i == j:
                L.append((-1, -1, -float(g0[i])))
            for (a1, b1, c1) in L:
                for (a2, b2, c2) in L:
                    holo = tuple(sorted(x for x in (a1, b2) if x >= 0))
                    anti = tuple(sorted(x for x in (b1, a2) if x >= 0))
                    F[(holo, anti)] += c1 * c2
    return nd, dict(F)


def charge0_basis(nd):
    """Real charge-0 basis: each element a real linear form in complex monomials
    c_h conj(c_a), returned as (label, [ (holo_tuple, anti_tuple, complex_coeff) ])."""
    B = [("const", [((), (), 1.0 + 0j)])]
    for a in range(nd):
        for b in range(a, nd):
            if a == b:
                B.append((("sym", a, a), [((a,), (a,), 1.0 + 0j)]))
            else:  # Re(c_a cbar_b)
                B.append((("sym", a, b), [((a,), (b,), 0.5 + 0j), ((b,), (a,), 0.5 + 0j)]))
    for a in range(nd):
        for b in range(a + 1, nd):  # Im(c_a cbar_b)
            B.append((("asym", a, b), [((a,), (b,), -0.5j), ((b,), (a,), 0.5j)]))
    return B


def _mul(t1, t2):
    (h1, a1, c1), (h2, a2, c2) = t1, t2
    return (tuple(sorted(h1 + h2)), tuple(sorted(a1 + a2)), c1 * c2)


def certify(d, N, g0, eps=1e-9, want_cert=False):
    """Maximize delta with f - delta = B^T Q B + mu(||c||^2-1), Q>>0. Native Clarabel."""
    import clarabel
    nd, F = f_coeffs(d, N, g0)
    B = charge0_basis(nd)
    nz = len(B)
    keys = {}

    def kidx(k):
        return keys.setdefault(k, len(keys))

    psd_pairs = [(i, j) for j in range(nz) for i in range(j + 1)]  # col-major svec
    npsd = len(psd_pairs)
    nmu = nz
    nx = npsd + nmu + 1
    idelta = nx - 1
    imu0 = npsd
    keymap = defaultdict(list)

    # sigma = B^T Q B ; svec entry (P<=R) scales the off-diagonal by sqrt2
    for col, (P, R) in enumerate(psd_pairs):
        scale = 1.0 if P == R else R2
        for t1 in B[P][1]:
            for t2 in B[R][1]:
                h, a, cc = _mul(t1, t2)
                keymap[(h, a)].append((col, scale * cc))
    # mu (||c||^2 - 1)
    for S in range(nz):
        col = imu0 + S
        for (h0, a0, c0) in B[S][1]:
            for k in range(nd):  # times +sum_k c_k cbar_k
                h, a, cc = _mul((h0, a0, c0), ((k,), (k,), 1.0 + 0j))
                keymap[(h, a)].append((col, cc))
            keymap[(h0, a0)].append((col, -c0))  # times -1
    keymap[((), ())].append((idelta, 1.0 + 0j))  # -delta on empty (LHS: f - delta)

    allkeys = sorted(set(keymap) | set(F))
    rows, cols, dat, b_vec = [], [], [], []
    rr = 0
    for k in allkeys:
        fv = F.get(k, 0)
        for part, target in ((np.real, float(np.real(fv))), (np.imag, float(np.imag(fv)))):
            for (col, cc) in keymap.get(k, []):
                v = part(cc)
                if abs(v) > 1e-14:
                    rows.append(rr); cols.append(col); dat.append(float(v))
            b_vec.append(target); rr += 1
    neq = rr
    for i in range(npsd):  # PSD slack: s = x_svec
        rows.append(neq + i); cols.append(i); dat.append(-1.0)

    A = sp.coo_matrix((dat, (rows, cols)), shape=(neq + npsd, nx)).tocsc()
    b = np.array(b_vec + [0.0] * npsd)
    P = sp.csc_matrix((nx, nx))
    q = np.zeros(nx); q[idelta] = -1.0
    cones = [clarabel.ZeroConeT(neq), clarabel.PSDTriangleConeT(nz)]
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = eps
    t0 = time.time()
    sol = clarabel.DefaultSolver(P, q, A, b, cones, settings).solve()
    out = dict(status=str(sol.status), delta=float(sol.x[idelta]), nz=nz,
               neq=neq, npsd=npsd, secs=time.time() - t0)
    if want_cert:
        out["residual"] = _identity_residual(sol.x, psd_pairs, B, nd, F, idelta, imu0)
    return out


def _identity_residual(x, psd_pairs, B, nd, F, idelta, imu0):
    """max_k | coeff_k( B^T Q B + mu(||c||^2-1) + delta*[empty] ) - F_k |, numeric check."""
    acc = defaultdict(complex)
    for col, (P, R) in enumerate(psd_pairs):
        scale = 1.0 if P == R else R2
        for t1 in B[P][1]:
            for t2 in B[R][1]:
                h, a, cc = _mul(t1, t2)
                acc[(h, a)] += x[col] * scale * cc
    for S in range(len(B)):
        v = x[imu0 + S]
        for (h0, a0, c0) in B[S][1]:
            for k in range(nd):
                h, a, cc = _mul((h0, a0, c0), ((k,), (k,), 1.0 + 0j))
                acc[(h, a)] += v * cc
            acc[(h0, a0)] += v * (-c0)
    acc[((), ())] += x[idelta]
    keysall = set(acc) | set(F)
    return max(abs(acc.get(k, 0) - F.get(k, 0)) for k in keysall)


def main():
    jobs = [("ATTAIN(2,4)", 4, 2, [0.7, 0.7, 0.3, 0.3]),
            ("EXTER (2,4)", 4, 2, [0.8, 0.6, 0.4, 0.2])]
    print("charge-0 Hermitian-SOS non-attainability certifier -- (2,4) self-test")
    print("(2,4) pure spectra are degenerate pairs (a,a,1-a,1-a); EXTER breaks it.\n")
    for tag, d, N, g0 in jobs:
        r = certify(d, N, g0, want_cert=True)
        verdict = ("NON-ATTAINABLE (delta>0)" if r["delta"] > 1e-6
                   else "attainable / boundary (delta~0)")
        print(f"{tag}: {r['status']}  delta = {r['delta']:.6g}  -> {verdict}")
        print(f"           identity residual {r['residual']:.2e}  "
              f"(nz {r['nz']}, {r['neq']} eq, {r['secs']:.2f}s)")
    print("\nExpected: EXTER delta = 0.04 (= squared distance to (0.7,0.7,0.3,0.3)),")
    print("          ATTAIN delta ~ 0. See module docstring for the (3,11) scaling wall.")


if __name__ == "__main__":
    main()
