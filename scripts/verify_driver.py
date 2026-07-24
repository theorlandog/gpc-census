#!/usr/bin/env python3
"""Verify the Hamiltonian driver for the v_B fiber anholonomy.

Stages:
A. Rebuild the fiber circle exactly; verify_exact certificates on both branches.
B. Monotonicity of cos Theta(t) on the physical interval (circle is a graph).
C. Independent recomputation of gamma_B: Berry connection integral and
   discrete Bargmann product, under BOTH lift conventions (phase on D0 vs D1),
   to pin down which convention the recorded -0.6861 pi belongs to.
D. Driver test: H(phi) = -(e^{i phi} n_2 a8+ a4 + h.c.). Constrained minimizer
   over the fiber circle vs phi: uniqueness, continuity, winding number.
E. Berry phase (Bargmann) of the tracked minimizer loop.
"""
from __future__ import annotations
import sys
from fractions import Fraction
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))

DETS = [(0, 1, 2, 4), (0, 1, 2, 8), (0, 1, 3, 5), (0, 2, 3, 6),
        (0, 3, 4, 7), (0, 3, 7, 8), (1, 3, 6, 8), (2, 3, 4, 8)]
W0 = np.array([1, 8, 4, 3, 2, 2, 1, 2], float)
KER = np.array([1, -1, 0, 0, -1, 1, 0, 0], float)
DEN = 23
SPEC_INT = (20, 14, 14, 14, 14, 4, 4, 4, 4)

def c_of_t(t):
    return (3 + 7*t - 2*t**2) / (2*np.sqrt((1+t)*(8-t)*(4-t*t)))

T_LO = 7/17 - 18*np.sqrt(35)/85
T_HI = 7/17 + 18*np.sqrt(35)/85

# ---------- Fock sector machinery (dim 126) ----------
from itertools import combinations
BASIS = [tuple(c) for c in combinations(range(9), 4)]
IDX = {T: i for i, T in enumerate(BASIS)}
DIM = len(BASIS)

def adag_a(p, q):
    """Matrix of a†_p a_q on the 4-particle sector, ordered-tuple convention."""
    M = np.zeros((DIM, DIM))
    for T in BASIS:
        if q not in T:
            continue
        s1 = (-1) ** T.index(q)
        rest = tuple(x for x in T if x != q)
        if p in rest:
            continue
        Tp = tuple(sorted(rest + (p,)))
        s2 = (-1) ** Tp.index(p)
        M[IDX[Tp], IDX[T]] += s1 * s2
    return M

N_OP = [adag_a(m, m) for m in range(9)]
HOP84 = adag_a(8, 4)                       # a8+ a4
X_OP = N_OP[2] @ HOP84                      # n_2 a8+ a4  (2 not in {4,8}: order-safe)

def one_rdm(psi):
    rho = np.zeros((9, 9), complex)
    for p in range(9):
        for q in range(9):
            rho[p, q] = np.vdot(psi, adag_a(p, q) @ psi)
    return rho

# ---------- fiber states ----------
def theta_of_t_grid(n=4001):
    ts = np.linspace(T_LO, T_HI, n)
    cs = c_of_t(ts)
    return ts, cs

def t_of_theta(theta_abs, ts, cs):
    """Invert cos Theta = c(t) given monotone c; theta_abs in [0, pi]."""
    c = np.cos(theta_abs)
    # c monotone in t (checked in stage B); interp requires ascending x
    if cs[0] < cs[-1]:
        return np.interp(c, cs, ts)
    return np.interp(c, cs[::-1], ts[::-1])

def fiber_state(theta_signed, ts, cs, carrier=0):
    """Closed lift: phase e^{i theta_signed} on DETS[carrier], others real."""
    t = t_of_theta(abs(theta_signed), ts, cs)
    w = W0 + t * KER
    amp = np.sqrt(np.clip(w, 0, None) / DEN).astype(complex)
    amp[carrier] *= np.exp(1j * theta_signed)
    psi = np.zeros(DIM, complex)
    for i, T in enumerate(DETS):
        psi[IDX[T]] = amp[i]
    return psi, t

# ================= Stage A: exact certificates =================
print('=== A. exact fiber certificates (verify_exact) ===')
import sympy as sp
from gpc_census.exactify import verify_exact

def cert(tv_rat, branch):
    t = sp.Rational(*tv_rat)
    w = [sp.Integer(int(W0[i])) + t * sp.Integer(int(KER[i])) for i in range(8)]
    c = (3 + 7*t - 2*t**2) / (2*sp.sqrt((1+t)*(8-t)*(4-t**2)))
    c = sp.simplify(c)
    ph = sp.exp(sp.I * branch * sp.acos(c))
    amps = [sp.sqrt(sp.Rational(int(w[i].p), int(w[i].q)) / DEN) if w[i].is_Rational
            else sp.sqrt(w[i] / DEN) for i in range(8)]
    amps = [sp.sqrt(w[i] / DEN) for i in range(8)]
    amps[0] = amps[0] * ph
    spec = [Fraction(x, DEN) for x in SPEC_INT]
    return verify_exact(4, 9, spec, DETS, amps)

for tv, br in [((0, 1), +1), ((0, 1), -1), ((1, 2), +1), ((-1, 2), -1)]:
    ok = cert(tv, br)
    print(f'  t={tv[0]}/{tv[1]} branch {br:+d}: verify_exact = {ok}')

# ================= Stage B: monotonicity / circle =================
print('=== B. circle structure ===')
ts, cs = theta_of_t_grid()
dc = np.diff(cs)
print(f'  walls: t_lo={T_LO:.10f}, t_hi={T_HI:.10f}, inside positivity (-1,2): '
      f'{T_LO > -1 and T_HI < 2}')
print(f'  c(t_lo)={cs[0]:+.6f}, c(t_hi)={cs[-1]:+.6f}')
print(f'  c(t) strictly monotone on [t_lo,t_hi]: {bool(np.all(dc > 0) or np.all(dc < 0))}'
      f'  (min|dc|={np.abs(dc).min():.2e})')
# spectrum check across the circle
worst = 0.0
spec_sorted = np.sort(np.array(SPEC_INT, float))[::-1]
for th in np.linspace(-np.pi, np.pi, 41):
    psi, _ = fiber_state(th, ts, cs)
    ev = np.sort(np.linalg.eigvalsh(one_rdm(psi)))[::-1] * DEN
    worst = max(worst, float(np.max(np.abs(ev - spec_sorted))))
print(f'  1-RDM spectrum match across circle: worst |dev| = {worst:.2e}')

# ================= Stage C: gamma_B, both carriers =================
print('=== C. Berry phase of the fiber loop (independent recomputation) ===')
def bargmann(loop_states):
    z = 1.0 + 0j
    n = len(loop_states)
    for k in range(n):
        z *= np.vdot(loop_states[k], loop_states[(k + 1) % n])
    return -np.angle(z)

for carrier, name in [(0, 'D0=(0,1,2,4), w0=1+t'), (1, 'D1=(0,1,2,8), w1=8-t')]:
    for N in (2000, 8000, 32000):
        thetas = np.linspace(-np.pi, np.pi, N, endpoint=False)
        loop = [fiber_state(th, ts, cs, carrier=carrier)[0] for th in thetas]
        g = bargmann(loop)
        conn = None
        if N == 32000:
            # connection integral: gamma = - oint w_carrier(t)/23 dTheta
            wc = np.array([W0[carrier] + KER[carrier] *
                           t_of_theta(abs(th), ts, cs) for th in thetas])
            conn = -np.trapezoid(np.r_[wc, wc[:1]] / DEN,
                                 np.r_[thetas, [np.pi]])
            print(f'  carrier {name}: Bargmann N={N}: gamma = {g:+.6f} rad '
                  f'= {g/np.pi:+.6f} pi ; connection integral = {conn:+.6f} rad '
                  f'= {conn/np.pi:+.6f} pi')
        else:
            print(f'  carrier {name}: Bargmann N={N}: gamma = {g:+.6f} rad '
                  f'= {g/np.pi:+.6f} pi')
print('  recorded value: -2.1555 rad = -0.6861 pi (i.e. +1.3139 pi mod 2pi)')

# ================= Stage D: driver tracking =================
print('=== D. driver tracking: H(phi) = -(e^{i phi} n_2 a8+ a4 + h.c.) ===')
NTH = 4096
thetas = np.linspace(-np.pi, np.pi, NTH, endpoint=False)
Z = np.empty(NTH, complex)   # z(theta) = <psi| X |psi>
STATES = []
for i, th in enumerate(thetas):
    psi, _ = fiber_state(th, ts, cs)
    STATES.append(psi)
    Z[i] = np.vdot(psi, X_OP @ psi)
print(f'  coherence |z| range over circle: [{np.abs(Z).min():.4f}, {np.abs(Z).max():.4f}]')

NPHI = 1441
phis = np.linspace(0, 2*np.pi, NPHI)
argmins = np.empty(NPHI, int)
gaps = np.empty(NPHI)        # global min vs best competing local min
for j, phi in enumerate(phis):
    E = -2*np.real(np.exp(1j*phi) * Z)
    k = int(np.argmin(E))
    argmins[j] = k
    # competing local minima away from the global one
    locmin = [i for i in range(NTH)
              if E[i] <= E[(i-1) % NTH] and E[i] <= E[(i+1) % NTH]]
    far = [i for i in locmin if min(abs(i-k), NTH-abs(i-k)) > NTH//16]
    gaps[j] = (min(E[i] for i in far) - E[k]) if far else np.inf
th_star = thetas[argmins]
# unwrap and winding
unw = np.unwrap(th_star)
winding = (unw[-1] - unw[0]) / (2*np.pi)
jumps = np.abs(np.diff(unw))
print(f'  max step in tracked Theta*(phi): {jumps.max():.4f} rad '
      f'(grid ~ {2*np.pi/NTH:.4f}; continuous tracking: {jumps.max() < 0.05})')
print(f'  winding number of Theta* over phi in [0,2pi): {winding:+.4f}')
finite = gaps[np.isfinite(gaps)]
print(f'  competing-basin margin: min over phi of (best competing local min - global min) '
      f'= {finite.min() if finite.size else float("inf"):.4f} '
      f'({np.sum(np.isfinite(gaps))}/{NPHI} phis have any competing local min)')

# ================= Stage E: Berry phase of the tracked loop =================
print('=== E. Berry phase of the driver-tracked loop ===')
loop = [STATES[k] for k in argmins[:-1]]   # phi=2pi duplicates phi=0
gd = bargmann(loop)
print(f'  Bargmann phase of tracked loop (D0-carrier lift): {gd:+.6f} rad = {gd/np.pi:+.6f} pi')
loopD1 = [fiber_state(thetas[k], ts, cs, carrier=1)[0] for k in argmins[:-1]]
gd1 = bargmann(loopD1)
print(f'  Bargmann phase of tracked loop (D1-carrier lift): {gd1:+.6f} rad = {gd1/np.pi:+.6f} pi')
