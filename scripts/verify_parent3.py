#!/usr/bin/env python3
"""Continuous parent family, both carrier lifts; Berry phase of driven loop."""
import numpy as np
from pathlib import Path
exec(open(Path(__file__).parent / 'verify_parent2.py').read().split('# tracked Theta')[0].replace(
    'print(f\'dictionary', '_p=(f\'dictionary'))  # reuse defs, mute print

NTH = 2048
thetas = np.linspace(-np.pi, np.pi, NTH, endpoint=False)
X2 = N1[2] @ HOP
Zc = np.array([np.vdot(fiber_state(th), X2 @ fiber_state(th)) for th in thetas])

def run(carrier, nphi=96):
    phig = np.linspace(0, 2*np.pi, nphi+1)[:-1]
    x_prev = None
    gs_loop, steps, gaps = [], [], []
    ok_all = True
    for phi in phig:
        Ec = -2*np.real(np.exp(1j*phi)*Zc)
        th_star = thetas[int(np.argmin(Ec))]
        psi = fiber_state(th_star, carrier=carrier)
        free, Hdrive, x0, resid, null = solve_eigcond(phi, psi)
        if x_prev is not None and null.shape[1] > 0:
            y = null.T @ (x_prev - x0)          # nearest point to previous coeffs
            x = x0 + null @ y
        else:
            x, y = x0, np.zeros(null.shape[1])
        H = build(free, Hdrive, x)
        ev, V = np.linalg.eigh(H)
        gap = ev[1]-ev[0]; fid = abs(np.vdot(V[:,0],psi))**2
        if gap < 0.05 or fid < 1-1e-9:          # ascend only if needed
            y, lam = ascend_gap(free, Hdrive, x0, null, psi, y0=y, iters=120)
            x = x0 + null @ y
            H = build(free, Hdrive, x)
            ev, V = np.linalg.eigh(H)
            gap = ev[1]-ev[0]; fid = abs(np.vdot(V[:,0],psi))**2
        ok_all &= (resid < 1e-9 and fid > 1-1e-9 and gap > 1e-3)
        if x_prev is not None:
            steps.append(np.linalg.norm(x - x_prev))
        x_prev = x
        gaps.append(gap)
        gs_loop.append(V[:,0]*np.exp(-1j*np.angle(np.vdot(psi, V[:,0]))))
    g = bargmann(gs_loop)
    print(f'carrier D{carrier}: all-phi PASS={ok_all}  min gap={min(gaps):.4f}  '
          f'max coeff step={max(steps):.4f}  '
          f'driven-loop Berry phase = {g:+.6f} rad = {g/np.pi:+.6f} pi')
    return g

print('=== continuous parent family, 96-point winding, both lifts ===')
g0 = run(0)
g1 = run(1)
print('fiber-loop references (stage C, forward orientation): '
      'D0-carrier -0.096502 pi, D1-carrier -0.686107 pi; driver winds -1 (reversed)')
