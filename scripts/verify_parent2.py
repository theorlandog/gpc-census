#!/usr/bin/env python3
"""Stage F (rebuilt): explicit <=2-body parent Hamiltonian for the tracked loop.

Dictionary of hermitian number-conserving operators, real coefficients:
  - n_i (9), n_i n_j (36)                                  [diagonal]
  - correlated one-electron hops n_m (a8+ a4 + h.c.) and the i-quadrature,
    m in modes\{4,8\} (14), plus the bare hop quadratures (2)
  - two-electron pair hops for every 2-orbital difference channel that
    connects two support determinants, both quadratures.

For each phi on the winding grid:
  1. Theta*(phi) from the stage-D driver (the coupling circle minimizer).
  2. LSQ for the eigenstate condition Q H psi = 0 over the dictionary,
     with the driven quadrature pair (n_2-hop cos/sin) CONSTRAINED to
     (cos phi, sin phi) * g so the drive is explicit; all else free.
  3. Within the exact-eigenstate affine family, subgradient ascent on
     lambda_min of P_perp (H - E_psi) P_perp  ->  make psi the unique
     gapped ground state.
Report residual, fidelity of the true ground state to psi, gap, coefficient
continuity, and the Berry phase of the exact ground-state loop.
"""
from __future__ import annotations
import numpy as np
from itertools import combinations

DETS = [(0, 1, 2, 4), (0, 1, 2, 8), (0, 1, 3, 5), (0, 2, 3, 6),
        (0, 3, 4, 7), (0, 3, 7, 8), (1, 3, 6, 8), (2, 3, 4, 8)]
W0 = np.array([1, 8, 4, 3, 2, 2, 1, 2], float)
KER = np.array([1, -1, 0, 0, -1, 1, 0, 0], float)
DEN = 23

BASIS = [tuple(c) for c in combinations(range(9), 4)]
IDX = {T: i for i, T in enumerate(BASIS)}
DIM = len(BASIS)

def op_adag_a(p, q):
    M = np.zeros((DIM, DIM))
    for T in BASIS:
        if q not in T:
            continue
        s1 = (-1) ** T.index(q)
        rest = tuple(x for x in T if x != q)
        if p in rest:
            continue
        Tp = tuple(sorted(rest + (p,)))
        M[IDX[Tp], IDX[T]] += s1 * (-1) ** Tp.index(p)
    return M

def op_pairhop(p1, p2, q1, q2):
    """a+_{p1} a+_{p2} a_{q2} a_{q1} with p1<p2, q1<q2 (one orientation)."""
    M = np.zeros((DIM, DIM))
    for T in BASIS:
        if q1 not in T or q2 not in T:
            continue
        # apply a_{q1} then a_{q2}? use a+p1 a+p2 aq2 aq1 acting right-to-left
        s = 1.0
        cur = list(T)
        for q in (q1, q2):
            i = cur.index(q)
            s *= (-1) ** i
            cur.pop(i)
        okk = True
        for p in (p2, p1):
            if p in cur:
                okk = False
                break
            cur = sorted(cur + [p])
            s *= (-1) ** cur.index(p)
        if not okk:
            continue
        M[IDX[tuple(cur)], IDX[T]] += s
    return M

N1 = [op_adag_a(m, m) for m in range(9)]
HOP = op_adag_a(8, 4)

OPS = []          # hermitian matrices
NAMES = []
for i in range(9):
    OPS.append(np.diag(np.array([1.0 if i in T else 0 for T in BASIS])).astype(complex))
    NAMES.append(f'n{i}')
for i in range(9):
    for j in range(i + 1, 9):
        OPS.append(np.diag(np.array([1.0 if (i in T and j in T) else 0
                                     for T in BASIS])).astype(complex))
        NAMES.append(f'n{i}n{j}')
# hop quadratures
def herm_quads(M, tag):
    OPS.append((M + M.conj().T).astype(complex)); NAMES.append(tag + ':cos')
    OPS.append((1j * (M - M.conj().T)).astype(complex)); NAMES.append(tag + ':sin')
herm_quads(HOP, 'hop84')
for m in [0, 1, 2, 3, 5, 6, 7]:
    herm_quads(N1[m] @ HOP, f'n{m}*hop84')
# pair-hop channels between support dets
chan = set()
for a in range(8):
    for b in range(a + 1, 8):
        A, B = set(DETS[a]), set(DETS[b])
        if len(A - B) == 2:
            (q1, q2) = sorted(A - B)
            (p1, p2) = sorted(B - A)
            chan.add((p1, p2, q1, q2))
for (p1, p2, q1, q2) in sorted(chan):
    herm_quads(op_pairhop(p1, p2, q1, q2), f'ph({p1}{p2}<-{q1}{q2})')
DRIVE_C = NAMES.index('n2*hop84:cos')
DRIVE_S = NAMES.index('n2*hop84:sin')
NP_ = len(OPS)
print(f'dictionary: {NP_} hermitian operators '
      f'({len(chan)} pair-hop channels: {sorted(chan)})')

# fiber
def c_of_t(t):
    return (3 + 7*t - 2*t**2) / (2*np.sqrt((1+t)*(8-t)*(4-t*t)))
T_LO = 7/17 - 18*np.sqrt(35)/85
T_HI = 7/17 + 18*np.sqrt(35)/85
ts = np.linspace(T_LO, T_HI, 4001)
cs = c_of_t(ts)
def fiber_state(th, carrier=0):
    t = np.interp(np.cos(abs(th)), cs, ts)
    w = W0 + t * KER
    amp = np.sqrt(np.clip(w, 0, None) / DEN).astype(complex)
    amp[carrier] *= np.exp(1j * th)
    psi = np.zeros(DIM, complex)
    for i, T in enumerate(DETS):
        psi[IDX[T]] = amp[i]
    return psi

G = 1.0
def solve_eigcond(phi, psi):
    """LSQ over free coeffs with drive pinned to g(cos phi, sin phi)."""
    free = [a for a in range(NP_) if a not in (DRIVE_C, DRIVE_S)]
    Hdrive = -G * (np.cos(phi) * OPS[DRIVE_C] + np.sin(phi) * OPS[DRIVE_S])
    b0 = Hdrive @ psi
    Qm = lambda v: v - psi * np.vdot(psi, v)
    cols = np.array([Qm(OPS[a] @ psi) for a in free]).T
    bvec = -Qm(b0)
    Ar = np.vstack([cols.real, cols.imag])
    br = np.concatenate([bvec.real, bvec.imag])
    x, _, rank, sv = np.linalg.lstsq(Ar, br, rcond=None)
    resid = float(np.linalg.norm(Ar @ x - br))
    _, s, Vt = np.linalg.svd(Ar, full_matrices=True)
    tol = max(Ar.shape) * np.finfo(float).eps * (s[0] if s.size else 1)
    r = int(np.sum(s > tol))
    null = Vt[r:].T
    return free, Hdrive, x, resid, null

def build(free, Hdrive, x):
    H = Hdrive.copy()
    for a, xa in zip(free, x):
        H = H + xa * OPS[a]
    return H

def ascend_gap(free, Hdrive, x0, null, psi, y0=None, iters=250):
    """maximize lambda_min of P_perp(H - E_psi)P_perp over the eigen-family."""
    P = np.eye(DIM) - np.outer(psi, psi.conj())
    def spec_bottom(x):
        H = build(free, Hdrive, x)
        E = np.real(np.vdot(psi, H @ psi))
        Hp = P @ (H - E * np.eye(DIM)) @ P
        ev, V = np.linalg.eigh(Hp)
        # one zero mode is psi itself inside P... psi is annihilated by P: gives eigenvalue 0
        # exclude the psi direction: find eigvec ~ psi
        ov = np.abs(V.conj().T @ psi)
        keep = np.argsort(ov)  # smallest overlap with psi kept
        idxs = [i for i in np.argsort(ev) if ov[i] < 0.5]
        i0 = idxs[0]
        return ev[i0], V[:, i0], E
    y = np.zeros(null.shape[1]) if y0 is None else y0.copy()
    lr = 0.5
    lam, v, E = spec_bottom(x0 + null @ y)
    for it in range(iters):
        # subgradient: d lambda / d x_a = <v|O_a|v> - <psi|O_a|psi>  (E depends on x too)
        g = np.array([np.real(np.vdot(v, OPS[free[a]] @ v))
                      - np.real(np.vdot(psi, OPS[free[a]] @ psi))
                      for a in range(len(free))])
        gy = null.T @ g
        ng = np.linalg.norm(gy)
        if ng < 1e-12:
            break
        y2 = y + lr * gy / ng
        lam2, v2, E2 = spec_bottom(x0 + null @ y2)
        if lam2 > lam:
            y, lam, v, E = y2, lam2, v2, E2
            lr = min(lr * 1.3, 2.0)
        else:
            lr *= 0.5
            if lr < 1e-6:
                break
    return y, lam

def bargmann(loop):
    z = 1.0 + 0j
    for k in range(len(loop)):
        z *= np.vdot(loop[k], loop[(k + 1) % len(loop)])
    return -np.angle(z)

# tracked Theta*(phi)
NTH = 2048
thetas = np.linspace(-np.pi, np.pi, NTH, endpoint=False)
X2 = N1[2] @ HOP
Zc = np.array([np.vdot(fiber_state(th), X2 @ fiber_state(th)) for th in thetas])

phigrid = np.linspace(0, 2 * np.pi, 17)[:-1]
gs_loop = []
y_prev = None
print(f'{"phi/pi":>7} {"Th*/pi":>8} {"resid":>9} {"nulldim":>7} {"gap":>9} '
      f'{"GS fidelity":>11}')
ok_all = True
xs_full = []
for phi in phigrid:
    Ecurve = -2 * np.real(np.exp(1j * phi) * Zc)
    th_star = thetas[int(np.argmin(Ecurve))]
    psi = fiber_state(th_star)
    free, Hdrive, x0, resid, null = solve_eigcond(phi, psi)
    if y_prev is not None and null.shape[1] == y_prev.shape[0]:
        y0 = y_prev
    else:
        y0 = None
    y, lam = ascend_gap(free, Hdrive, x0, null, psi, y0)
    y_prev = y
    x = x0 + null @ y
    H = build(free, Hdrive, x)
    ev, V = np.linalg.eigh(H)
    fid = abs(np.vdot(V[:, 0], psi)) ** 2
    gap = ev[1] - ev[0]
    gs_loop.append(V[:, 0] * np.exp(-1j * np.angle(np.vdot(psi, V[:, 0]))))
    xf = np.zeros(NP_)
    for a, xa in zip(free, x):
        xf[a] = xa
    xf[DRIVE_C], xf[DRIVE_S] = -G * np.cos(phi), -G * np.sin(phi)
    xs_full.append(xf)
    ok = resid < 1e-9 and fid > 1 - 1e-9 and gap > 1e-3
    ok_all &= ok
    print(f'{phi/np.pi:7.3f} {th_star/np.pi:8.3f} {resid:9.2e} {null.shape[1]:7d} '
          f'{gap:9.4f} {fid:11.8f}{"" if ok else "  <-- FAIL"}')

xs_full = np.array(xs_full)
dx = np.linalg.norm(np.diff(np.vstack([xs_full, xs_full[:1]]), axis=0), axis=1)
print(f'coefficient continuity around the loop: max step {dx.max():.3f}, '
      f'mean {dx.mean():.3f} (coeff-vector norm ~ {np.linalg.norm(xs_full[0]):.2f})')
print(f'PARENT VERDICT (exact eigenstate + unique gapped ground state at all '
      f'{len(phigrid)} phis): {"PASS" if ok_all else "FAIL"}')
if ok_all:
    g = bargmann(gs_loop)
    print(f'Berry phase of the exact ground-state loop: {g:+.6f} rad = {g/np.pi:+.6f} pi')
    print('reference (fiber-loop, D0-carrier): -0.096502 pi ; D1-carrier: -0.686107 pi; '
          'tracked orientation is reversed (winding -1) so expect sign flip')
