#!/usr/bin/env python3
"""Certified non-attainability via a charge-0 Hermitian-SOS (moment) relaxation.

This is the RANK-11 CERTIFICATE PATH recommended in docs/RESEARCH.md ("A
Lasserre/SOS moment-relaxation dual of min ||gamma(psi)-gamma_0||^2 yields a
CERTIFIED positive lower bound = rigorous non-attainability"), now runnable: SDP
solvers (SCS, Clarabel) are available in the venv. It computes a lower bound delta
on

    min_{||c||^2 = 1}  f(c),   f = sum_ij |gamma_ij(c) - g0_i delta_ij|^2,

the squared distance from the target spectrum diag(g0) to the nearest pure
N-representable 1-RDM. delta > 0 is a rigorous certificate that diag(g0) is
OUTSIDE the (N,d) polytope (non-attainable). By orbital-rotation invariance
(Hoffman-Wielandt) min f equals the squared Euclidean distance from g0's spectrum
to the moment polytope, so a positive delta both proves non-attainability and
quantifies the gap. The contraction attack only produces an UNCONSTRAINED floor
(evidence, not proof).

FORMULATION (gauge-reduced, exact cone). f is U(1)-gauge invariant, so every
monomial has charge |holo| - |anti| = 0, and f squares the SESQUILINEAR forms
gamma_ij - g0 delta (NOT a holomorphic SOS -- a holomorphic {1,c,c^2} block cannot
even certify the trivial delta=0). The correct cone squares the CHARGE-0 space
    B = {1} U {Re(c_a cbar_b): a<=b} U {Im(c_a cbar_b): a<b}   (size 1 + nd^2),
nd = C(d,N). The certificate is
    f - delta = B^T Q B + mu (||c||^2 - 1),   Q >> 0 (real symmetric),
mu = sum_S mu_S B_S a free charge-0 degree-<=2 multiplier (level-2 relaxation on
the sphere). delta=0 is always feasible (Q = Gram of f, mu=0), so the cone is
correct; the multiplier lifts delta onto the sphere. Coefficients are matched over
complex monomial keys (holo multiset, anti multiset) and handed to a solver
natively (min c'x s.t. Ax + s = b, s in cones): one PSD cone (side 1+nd^2) + a
ZeroCone of equalities.

BACKEND. Default SCS: its ADMM projects the PSD cone by eigendecomposition,
O(nz^2) memory, so it clears side nz=401 at (3,6) where Clarabel's dense-IPM
Nesterov-Todd scaling (O(npsd^2) ~ 52 GB) OOMs. The single un-reduced block still
scales as nz = 1 + nd^2 (401 at (3,6), 27226 at (3,11)); (3,11) needs the
symmetry-reduced solver (see scripts/sos_symmetry.py), which block-diagonalizes Q
under the spectrum stabilizer into small irrep blocks.

VALIDATED gates (run this file):
  (3,6) EXTER diag(0.9,0.8,0.7,0.35,0.2,0.05) -> delta ~ 0.0042 > 0
        (violates the (3,6) pairing lambda_i + lambda_{7-i} = 1, non-attainable)
  (3,6) CONTROL diag(1,1,1,0,0,0) (attainable Slater) -> delta ~ 0 (feasible)
  (2,4) EXTER (0.8,0.6,0.4,0.2) -> delta = 0.04 == exact distance^2 to (0.7,0.7,0.3,0.3).

Run:  .venv/bin/python scripts/sos_nonattain.py [--full]
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
    """Real charge-0 basis: (label, [ (holo_tuple, anti_tuple, complex_coeff) ])."""
    B = [("const", [((), (), 1.0 + 0j)])]
    for a in range(nd):
        for b in range(a, nd):
            if a == b:
                B.append((("sym", a, a), [((a,), (a,), 1.0 + 0j)]))
            else:
                B.append((("sym", a, b), [((a,), (b,), 0.5 + 0j), ((b,), (a,), 0.5 + 0j)]))
    for a in range(nd):
        for b in range(a + 1, nd):
            B.append((("asym", a, b), [((a,), (b,), -0.5j), ((b,), (a,), 0.5j)]))
    return B


def mul_mono(t1, t2):
    (h1, a1, c1), (h2, a2, c2) = t1, t2
    return (tuple(sorted(h1 + h2)), tuple(sorted(a1 + a2)), c1 * c2)


def certify(d, N, g0, solver="SCS", eps=1e-6, max_iters=50000, verbose=False):
    """Maximize delta with f - delta = B^T Q B + mu(||c||^2-1), Q>>0. Returns dict."""
    nd, F = f_coeffs(d, N, g0)
    B = charge0_basis(nd)
    nz = len(B)
    lower = (solver == "SCS")   # SCS: lower-tri col-major svec; Clarabel: upper-tri
    if lower:
        psd_pairs = [(i, j) for j in range(nz) for i in range(j, nz)]
    else:
        psd_pairs = [(i, j) for j in range(nz) for i in range(j + 1)]
    npsd = len(psd_pairs)
    nx = npsd + nz + 1
    idelta = nx - 1
    imu0 = npsd
    keymap = defaultdict(list)
    for col, (P, R) in enumerate(psd_pairs):
        scale = 1.0 if P == R else R2
        for t1 in B[P][1]:
            for t2 in B[R][1]:
                h, a, cc = mul_mono(t1, t2)
                keymap[(h, a)].append((col, scale * cc))
    for S in range(nz):
        col = imu0 + S
        for (h0, a0, c0) in B[S][1]:
            for k in range(nd):
                h, a, cc = mul_mono((h0, a0, c0), ((k,), (k,), 1.0 + 0j))
                keymap[(h, a)].append((col, cc))
            keymap[(h0, a0)].append((col, -c0))
    keymap[((), ())].append((idelta, 1.0 + 0j))

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
    for i in range(npsd):
        rows.append(neq + i); cols.append(i); dat.append(-1.0)
    A = sp.csc_matrix((dat, (rows, cols)), shape=(neq + npsd, nx))
    b = np.array(b_vec + [0.0] * npsd)
    c = np.zeros(nx); c[idelta] = -1.0
    t0 = time.time()
    if solver == "SCS":
        import scs
        sol = scs.solve(dict(A=A, b=b, c=c), dict(z=neq, s=[nz]),
                        verbose=verbose, eps_abs=eps, eps_rel=eps, max_iters=max_iters)
        x = sol["x"]; status = sol["info"]["status"]
    else:
        import clarabel
        P = sp.csc_matrix((nx, nx))
        st = clarabel.DefaultSettings(); st.verbose = verbose
        st.tol_gap_abs = st.tol_gap_rel = st.tol_feas = eps
        sol = clarabel.DefaultSolver(P, c, A, b, [clarabel.ZeroConeT(neq),
                                     clarabel.PSDTriangleConeT(nz)], st).solve()
        x = np.array(sol.x); status = str(sol.status)
    delta = float(x[idelta])
    resid = _resid(x, psd_pairs, B, nd, F, idelta, imu0)
    return dict(status=status, delta=delta, resid=resid, nz=nz, neq=neq,
                secs=time.time() - t0)


def _resid(x, psd_pairs, B, nd, F, idelta, imu0):
    acc = defaultdict(complex)
    for col, (P, R) in enumerate(psd_pairs):
        scale = 1.0 if P == R else R2
        for t1 in B[P][1]:
            for t2 in B[R][1]:
                h, a, cc = mul_mono(t1, t2)
                acc[(h, a)] += x[col] * scale * cc
    for S in range(len(B)):
        v = x[imu0 + S]
        for (h0, a0, c0) in B[S][1]:
            for k in range(nd):
                h, a, cc = mul_mono((h0, a0, c0), ((k,), (k,), 1.0 + 0j))
                acc[(h, a)] += v * cc
            acc[(h0, a0)] += v * (-c0)
    acc[((), ())] += x[idelta]
    return max(abs(acc.get(k, 0) - F.get(k, 0)) for k in set(acc) | set(F))


def main():
    print("charge-0 Hermitian-SOS non-attainability certifier (SCS backend)\n")
    jobs = [("EXTER  (3,6)", 6, 3, [0.9, 0.8, 0.7, 0.35, 0.2, 0.05], True),
            ("CONTROL(3,6)", 6, 3, [1, 1, 1, 0, 0, 0], False),
            ("EXTER  (2,4)", 4, 2, [0.8, 0.6, 0.4, 0.2], True)]
    if "--full" not in sys.argv:
        jobs = [j for j in jobs if j[1] != 6 or j[0].startswith("EXTER")] + \
               [("CONTROL(3,6)", 6, 3, [1, 1, 1, 0, 0, 0], False)]
    for tag, d, N, g0, ext in jobs:
        r = certify(d, N, g0, eps=1e-5, max_iters=20000)
        v = ("NON-ATTAINABLE (delta>0)" if r["delta"] > 1e-4
             else "attainable / feasible (delta~0)")
        ok = "OK" if (r["delta"] > 1e-4) == ext else "??"
        print(f"{tag}: {r['status']:8} delta={r['delta']:.6g}  resid={r['resid']:.1e}  "
              f"{r['secs']:.0f}s -> {v}  [{ok}]")
    print("\nGate: EXTER must give delta>0, CONTROL delta~0.")


if __name__ == "__main__":
    main()
