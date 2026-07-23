#!/usr/bin/env python3
"""Symmetry reduction of the open residual vertex v89 = (15,15,6^8)/26 of (3,10).

v89 has symmetry S_2 x S_8: swap the two occupation-15 orbitals (head = {0,1}) and
permute the eight occupation-6 orbitals (octet = {2,...,9}). The 3-subset Slater
determinants fall into THREE orbits: {3 octet} (56), {1 head, 2 octet} (56),
{2 head, 1 octet} (8). A fully symmetric state assigns one amplitude per orbit
(x, y, z), collapsing the 120-determinant search to three parameters, and its 1-RDM
block-diagonalizes exactly:

    antisym-head eigenvalue  a - b,
    octet standard rep       c - d     (multiplicity 7),
    2x2 [[a+b, 4e],[4e, c+7d]]         (symmetric head and octet modes),

with a,b,c,d,e the five distinct RDM entries, each a quadratic form in (x,y,z).

RESULT (this script): the fully-symmetric all-positive ("design", trivial-rep)
ansatz PROVABLY CANNOT attain v89. The octet off-diagonal is d = 15x^2+12y^2+z^2,
always positive and large, so the multiplicity-7 octet eigenvalue is
c - d = 6x^2 + 2y^2, strictly LESS than (6/26)*(56x^2+56y^2+8z^2) for every real
(x,y,z) other than 0. But that mult-7 eigenvalue must equal 6/26 (it is part of the
eightfold 6). Hence no fully-symmetric design attains v89: the vertex requires
INTERFERENCE that reduces the octet off-diagonal, i.e. intra-orbit sign/phase
structure transforming under a nontrivial S_8 representation -- which is exactly why
rational-grid and design searches miss it. This is a constructive re-derivation of
v89's interference character and pinpoints the next ansatz. Diagnostic; scripts/ is
ruff-excluded.

Run: .venv/bin/python scripts/v89_symmetry.py
"""
from __future__ import annotations

import itertools

import sympy as sp

HEAD = [0, 1]
OCTET = list(range(2, 10))
D = 10


def orbits():
    a = [tuple(sorted(c)) for c in itertools.combinations(OCTET, 3)]
    b = [tuple(sorted((h,) + c)) for h in HEAD for c in itertools.combinations(OCTET, 2)]
    c = [tuple(sorted((0, 1, o))) for o in OCTET]
    return {"A": a, "B": b, "C": c}


def rdm_entry(m, mp, orb, cf):
    """Unnormalized rho[m, mp] as a polynomial in the orbit amplitudes."""
    amap = {T: cf[U] for U in orb for T in orb[U]}
    val = sp.Integer(0)
    for T, ct in amap.items():
        if mp not in T:
            continue
        s1 = (-1) ** list(T).index(mp)
        rest = tuple(a for a in T if a != mp)
        if m in rest:
            continue
        tp = tuple(sorted(rest + (m,)))
        if tp not in amap:
            continue
        s2 = (-1) ** list(tp).index(m)
        val += s1 * s2 * amap[tp] * ct
    return sp.expand(val)


def main():
    x, y, z = sp.symbols("x y z", real=True)
    orb = orbits()
    cf = {"A": x, "B": y, "C": z}
    a = rdm_entry(0, 0, orb, cf)
    b = rdm_entry(0, 1, orb, cf)
    c = rdm_entry(2, 2, orb, cf)
    dblk = rdm_entry(2, 3, orb, cf)
    e = rdm_entry(0, 2, orb, cf)
    n2 = 56 * x ** 2 + 56 * y ** 2 + 8 * z ** 2
    print(f"orbit sizes: |A|={len(orb['A'])} |B|={len(orb['B'])} |C|={len(orb['C'])}")
    print("RDM entries (unnormalized quadratic forms):")
    for nm, ex in [("a head-occ", a), ("b head-offdiag", b), ("c octet-occ", c),
                   ("d octet-offdiag", dblk), ("e head-octet", e)]:
        print(f"  {nm:16} = {ex}")

    octet_eig = sp.expand(c - dblk)                       # mult-7 octet eigenvalue * n2
    target = sp.expand(sp.Rational(6, 26) * n2)           # required value * n2
    gap = sp.expand(target - octet_eig)                   # > 0 for all real (x,y,z) != 0
    print(f"\nmult-7 octet eigenvalue * n2 = c - d = {octet_eig}")
    print(f"required (6/26) * n2          = {target}")
    print(f"required - actual             = {gap}  (positive-definite: "
          f"{all(v > 0 for v in sp.Poly(gap, x, y, z).coeffs())})")
    print("=> the octet eigenvalue is STRICTLY below 6/26 for every real amplitude;")
    print("   the trivial-rep (design) symmetric ansatz cannot attain v89.")
    print("   v89 needs interference (signed/phased octet, nontrivial S_8 rep).")


if __name__ == "__main__":
    main()
