#!/usr/bin/env python3
"""Defining polynomial and genus of the v_B extremal fiber (derivation, not a fit).

The certified v_B state (4,9; spectrum (20,14,14,14,14,4,4,4,4)/23) sits in a
1-parameter family. Displace the eight support weights along the incidence-kernel
cycle v = (1,-1,0,0,-1,1,0,0): w(t) = w0 + t v. This keeps EVERY orbital
occupation fixed (the cycle leaves the diagonal 1-RDM invariant), so the 1-RDM
stays block-diagonal with a single nontrivial 2x2 block on modes (4,8), diagonal
(5,13)/23, which must split into the eigenvalues (14,4)/23. That pins the one
off-diagonal magnitude to 3/23 for all t. Exactly two determinant pairs drive it,
(D0,D1) and (D4,D5), with OPPOSITE fermionic signs, so the modulus condition is a
law of cosines with a relative minus, giving the exact defining polynomial

    F(t,c) = 4 c^2 (1+t)(8-t)(4-t^2) - (2 t^2 - 7 t - 3)^2 = 0,   c = cos(holonomy).

Substituting s = sqrt((1+t)(8-t)(4-t^2)) makes the curve birational to
s^2 = (1+t)(8-t)(2-t)(2+t), a quartic with four distinct roots -> a nonsingular
GENUS-1 elliptic curve. This script derives F, verifies it two ways (exact
verify_exact certificate at t=0 and t=1/2, numeric SPEC match across the domain),
and reports the genus and elliptic invariants. Diagnostic; scripts/ is
ruff-excluded.

Run: .venv/bin/python scripts/vb_fiber_ideal.py
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DETS = [(0, 1, 2, 4), (0, 1, 2, 8), (0, 1, 3, 5), (0, 2, 3, 6),
        (0, 3, 4, 7), (0, 3, 7, 8), (1, 3, 6, 8), (2, 3, 4, 8)]
W0 = [1, 8, 4, 3, 2, 2, 1, 2]
KER = [1, -1, 0, 0, -1, 1, 0, 0]
DEN = 23
SPEC_INT = (20, 14, 14, 14, 14, 4, 4, 4, 4)


def cos_branch(tv):
    """Physical holonomy-cosine branch c(t) as a sympy expression."""
    import sympy as sp
    return (3 + 7 * tv - 2 * tv ** 2) / (2 * sp.sqrt((1 + tv) * (8 - tv) * (4 - tv ** 2)))


def defining_polynomial():
    """F(t,c) with its bidegree and irreducibility over Q."""
    import sympy as sp
    t, c = sp.symbols("t c")
    F = sp.expand(4 * c ** 2 * (1 + t) * (8 - t) * (4 - t ** 2) - (2 * t ** 2 - 7 * t - 3) ** 2)
    return F, sp.degree(F, t), sp.degree(F, c), sp.Poly(F, t, c).is_irreducible


def verify_certificate(tv):
    """Exact verify_exact on the on-curve state at rational t (a proof)."""
    import sympy as sp
    from gpc_census.exactify import verify_exact
    w = [sp.Integer(W0[i]) + tv * KER[i] for i in range(8)]
    c = sp.simplify(cos_branch(tv))
    phase = [c + sp.I * sp.sqrt(1 - c ** 2)] + [sp.Integer(1)] * 7  # holonomy on D0
    amps = [sp.sqrt(sp.Rational(w[i], DEN)) * phase[i] for i in range(8)]
    spec = [Fraction(x, DEN) for x in SPEC_INT]
    return c, verify_exact(4, 9, spec, DETS, amps)


def numeric_domain_check():
    """Re-solved state matches SPEC eigenvalues across the real domain."""
    import numpy as np
    spec = np.array(sorted(SPEC_INT, reverse=True), float)
    worst = 0.0
    for k in range(-95, 200):
        tv = k / 100.0
        w = [W0[i] + tv * KER[i] for i in range(8)]
        rad = (1 + tv) * (8 - tv) * (4 - tv ** 2)
        if any(x < 0 for x in w) or rad <= 0:
            continue
        c = (3 + 7 * tv - 2 * tv ** 2) / (2 * np.sqrt(rad))
        if abs(c) > 1:
            continue
        th = np.arccos(c)
        amps = [np.sqrt(w[i] / DEN) * (np.exp(1j * th) if i == 0 else 1.0) for i in range(8)]
        amap = dict(zip(DETS, amps))
        rho = np.zeros((9, 9), complex)
        for T, ct in amap.items():
            for mp in T:
                s1 = (-1) ** T.index(mp)
                t2 = tuple(x for x in T if x != mp)
                for m in range(9):
                    if m in t2:
                        continue
                    tp = tuple(sorted(t2 + (m,)))
                    if tp not in amap:
                        continue
                    s2 = (-1) ** tp.index(m)
                    rho[m, mp] += s1 * s2 * np.conjugate(amap[tp]) * ct
        ev = np.sort(np.linalg.eigvalsh(rho))[::-1] * DEN
        worst = max(worst, float(np.max(np.abs(ev - spec))))
    return worst


def real_wall_states():
    """The arc endpoints |cos Theta| = 1: real extremal states over Q(sqrt35).

    Where the physical branch reaches cos Theta = +-1 the holonomy is 0 or pi, so
    every amplitude is real. Solving F(t, +-1) = 0 gives 85 t^2 - 70 t - 119 = 0,
    t = 7/17 -+ 18 sqrt35 / 85, both interior to the positivity range (-1, 2).
    Returns [(t, signs, verify_exact_ok, all_real)] for the two walls."""
    import sympy as sp
    from gpc_census.exactify import verify_exact
    t = sp.symbols("t")
    walls = sp.solve(sp.expand((2 * t ** 2 - 7 * t - 3) ** 2
                               - 4 * (1 + t) * (8 - t) * (4 - t ** 2)), t)
    spec = [Fraction(x, DEN) for x in SPEC_INT]
    out = []
    for tv in sorted(walls, key=float):
        w = [sp.nsimplify(W0[i] + tv * KER[i]) for i in range(8)]
        # cos Theta = -1 flips the D1 amplitude; cos Theta = +1 keeps all equal
        signs = [1, -1, 1, 1, 1, 1, 1, 1] if float(sp.re(cos_branch(tv))) < 0 \
            else [1, 1, 1, 1, 1, 1, 1, 1]
        amps = [signs[i] * sp.sqrt(sp.Rational(1, DEN) * w[i]) for i in range(8)]
        ok = verify_exact(4, 9, spec, DETS, amps)
        allreal = all(sp.im(sp.simplify(a)) == 0 for a in amps)
        out.append((sp.nsimplify(tv), signs, ok, allreal))
    return out


def paper_support_real_state():
    """Real extremal state for v_B on the PAPER's Theorem-3 support (Q(sqrt15)).

    The published psi_B lives on a different eight-determinant support (squared
    weights (4,2,7,3,2,2,1,2), block on orbitals (8,9) = (7,11)->(14,4),
    |rho_89|^2 = 21/23^2, cos gamma = 3/(4 sqrt14)). Its incidence cycle is
    (0,1,-1,0,-1,1,0,0); displacing the weights along it and demanding cos theta = 1
    (a real state) gives 109 t^2 - 110 t - 215 = 0, t = 55/109 -+ 42 sqrt15 / 109,
    both interior. This is the cleanest witness that v_B admits a real extremal
    state, on the paper's own support. Returns [(t, verify_exact_ok, all_real)]."""
    import sympy as sp
    from gpc_census.exactify import verify_exact
    dets = [(0, 1, 2, 5), (0, 1, 3, 7), (0, 1, 3, 8), (0, 2, 3, 4),
            (0, 2, 6, 7), (0, 2, 6, 8), (1, 2, 4, 7), (2, 3, 7, 8)]
    k0 = [4, 2, 7, 3, 2, 2, 1, 2]
    u = [0, 1, -1, 0, -1, 1, 0, 0]
    signs = [1, 1, 1, -1, 1, 1, -1, 1]  # psi_B signs; cos theta = 1 keeps them real
    spec = [Fraction(x, DEN) for x in SPEC_INT]
    walls = [sp.Rational(55, 109) - 42 * sp.sqrt(15) / 109,
             sp.Rational(55, 109) + 42 * sp.sqrt(15) / 109]
    out = []
    for tv in walls:
        k = [sp.nsimplify(k0[i] + tv * u[i]) for i in range(8)]
        amps = [signs[i] * sp.sqrt(sp.Rational(1, DEN) * k[i]) for i in range(8)]
        ok = verify_exact(4, 9, spec, dets, amps)
        allreal = all(sp.im(sp.simplify(a)) == 0 for a in amps)
        out.append((sp.nsimplify(tv), ok, allreal))
    return out


def genus_and_invariants():
    """Genus (via the birational quartic) and elliptic invariants I, J, j."""
    import sympy as sp
    t = sp.symbols("t")
    quartic = sp.expand((1 + t) * (8 - t) * (4 - t ** 2))
    roots = sp.roots(sp.Poly(quartic, t))
    squarefree = sp.gcd(quartic, sp.diff(quartic, t)) == 1
    a, b, c2, d, e = 1, -7, -12, 28, 32  # coeffs of quartic, high to low
    inv_i = 12 * a * e - 3 * b * d + c2 ** 2
    inv_j = 72 * a * c2 * e - 27 * a * d ** 2 - 27 * e * b ** 2 + 9 * b * c2 * d - 2 * c2 ** 3
    disc = 4 * inv_i ** 3 - inv_j ** 2
    j = sp.Rational(1728 * 4 * inv_i ** 3, disc)
    genus = 1 if (len(roots) == 4 and squarefree and disc != 0) else None
    return quartic, roots, squarefree, inv_i, inv_j, disc, j, genus


def main():
    import sympy as sp

    F, dt, dc, irr = defining_polynomial()
    print("defining polynomial of the v_B extremal fiber")
    print("  F(t,c) =", F)
    print(f"  bidegree (deg_t, deg_c) = ({dt}, {dc}); irreducible over Q: {irr}")
    print("  physical branch c(t) =", sp.simplify(cos_branch(sp.symbols('t'))))

    print("\nexact verify_exact certificates on the curve:")
    for tv in (sp.Integer(0), sp.Rational(1, 2)):
        c, ok = verify_certificate(tv)
        print(f"  t={tv}: cos(holonomy) = {c}   verify_exact = {ok}")

    worst = numeric_domain_check()
    print(f"\nnumeric SPEC match across the real domain t in (-1,2): max|eig-SPEC| = {worst:.2e}")

    print("\nreal extremal states at the arc endpoints (|cos Theta| = 1):")
    for tv, signs, ok, allreal in real_wall_states():
        print(f"  t={tv} ~ {float(tv):+.5f}: signs {tuple('+' if s>0 else '-' for s in signs)}"
              f"  verify_exact={ok}  all_real_amps={allreal}")
    print("  => v_B admits REAL extremal states over Q(sqrt35) (open question: YES)")

    print("\nreal extremal state on the PAPER's Theorem-3 support (Q(sqrt15)):")
    for tv, ok, allreal in paper_support_real_state():
        print(f"  t={tv} ~ {float(tv):+.5f}: verify_exact={ok}  all_real_amps={allreal}")

    quartic, roots, sf, inv_i, inv_j, disc, j, genus = genus_and_invariants()
    print("\ngenus via the birational quartic s^2 = (1+t)(8-t)(4-t^2):")
    print(f"  quartic = {quartic}; roots {dict(roots)}; squarefree {sf}")
    print(f"  quartic invariants I = {inv_i}, J = {inv_j}; 4I^3 - J^2 = {disc} (!= 0)")
    print(f"  j-invariant = {j} = {sp.factorint(j.p)} / {sp.factorint(j.q)}")
    print(f"  => nonsingular genus {genus} ELLIPTIC curve"
          " (the fiber is not rational)")


if __name__ == "__main__":
    main()
