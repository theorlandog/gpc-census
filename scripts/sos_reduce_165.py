#!/usr/bin/env python3
"""Cheap 165-dim symmetry reduction for the (3,11) cand 44 SOS certificate (G3).

The committed reduced solver (scripts/sos_symmetry_scaled.py) block-diagonalizes
the SOS Gram by an nz x nz Reynolds average + eigendecomposition on the charge-0
basis. That is O(nz^2 |G|) + O(nz^3): fine at (3,6) (nz = 381, 0.2 s), a wall at
cand 44 where the hop<=2 basis is nz = 17986 (measured: one Reynolds 199 s, one
eigh 41 s at nz = 4126 already, so ~3 h at nz = 17986). This module does the whole
block-structure computation in the 165-dim DETERMINANT rep instead (~2 s), which
is where the symmetry actually lives.

Two results, both matrix-free / cheap and both checked against the Burnside counts
recorded in docs/RESEARCH.md ("RANK-11 CERTIFICATE PATH"):

  1. decompose_V(): the 165 = C(11,3) signed determinant rep V under
     G = S5 (modes 0-4) x S6 (modes 5-10) into irrep copies, via a factorized
     Reynolds average on 165 x 165 (|S5| + |S6| = 840 conjugations, not |G| =
     86400). Reproduces V = 9 irreps, sum m_mu = 12, dim End_G(V) = sum m_mu^2 = 18.

  2. block_sizes(): the exact PSD block sizes of the reduced SOS Gram = the
     multiplicities of each S5 x S6 irrep nu in Herm(V) = V (x) Vbar, from exact
     Murnaghan-Nakayama characters (scripts/sn_characters.py):
         mult_nu = (1/|G|) sum_classes |c_a||c_b| chi_V(c_a,c_b)^2 chi_nu_a chi_nu_b.
     Reproduces 58 nonzero blocks, sum mult_nu^2 = 28884, LARGEST BLOCK = 58
     (the note's <= ~170 estimate; the true max is 58). So the un-reduced nz =
     27226 SDP reduces to 58 PSD blocks of side <= 58, trivially solvable.

This is the block structure the Clebsch-Gordan assembly of the reduced SDP targets;
the symmetry-adapted basis (combining V's irrep vectors into these nu-blocks) and
the exact rational Positivstellensatz endgame are the remaining pieces.

Run: .venv/bin/python scripts/sos_reduce_165.py
Diagnostic; scripts/ is ruff-excluded.
"""
from __future__ import annotations

import itertools
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from sn_characters import char_table, class_size  # noqa: E402
from sos_symmetry_scaled import signed_det_action  # noqa: E402


def det_basis(d, n):
    dets = [tuple(c) for c in itertools.combinations(range(d), n)]
    return dets, {t: i for i, t in enumerate(dets)}


def all_factor_actions(d, n, block, dets, idx):
    """signed (tgt, sgn) determinant actions for EVERY element of Sym(block)."""
    acts = []
    for pb in itertools.permutations(block):
        perm = list(range(d))
        for orig, img in zip(block, pb):
            perm[orig] = img
        acts.append(signed_det_action(perm, dets, idx))
    return acts


def reynolds(H, factor_acts):
    """Average H over G = prod factors by factorized signed-perm conjugation."""
    for acts in factor_acts:
        acc = np.zeros_like(H)
        for (tgt, sgn) in acts:
            acc[np.ix_(tgt, tgt)] += sgn[:, None] * H * sgn[None, :]
        H = acc / len(acts)
    return H


def decompose_V(d=11, n=3, factor_blocks=((0, 1, 2, 3, 4), (5, 6, 7, 8, 9, 10)), seed=0):
    """Return (irreps, factor_acts, dets). irreps = list of (dim, mult) in V."""
    dets, idx = det_basis(d, n)
    nd = len(dets)
    facts = [all_factor_actions(d, n, list(b), dets, idx) for b in factor_blocks]
    rng = np.random.default_rng(seed)
    H1 = rng.standard_normal((nd, nd)); H1 = H1 + H1.T
    A = reynolds(H1, facts)
    w, V = np.linalg.eigh(A)
    order = np.argsort(w); w = w[order]; V = V[:, order]
    comps, i = [], 0
    while i < nd:
        j = i + 1
        while j < nd and abs(w[j] - w[i]) < 1e-6:
            j += 1
        comps.append(V[:, i:j]); i = j
    H2 = rng.standard_normal((nd, nd)); H2 = H2 + H2.T
    B = reynolds(H2, facts)
    used = [False] * len(comps); irreps = []
    for k in range(len(comps)):
        if used[k]:
            continue
        dk = comps[k].shape[1]; mult = 1; used[k] = True
        for m in range(k + 1, len(comps)):
            if used[m] or comps[m].shape[1] != dk:
                continue
            M = comps[k].T @ B @ comps[m]; c2 = np.trace(M @ M.T) / dk
            if c2 > 1e-8 and np.allclose(M @ M.T, c2 * np.eye(dk), atol=1e-5 * max(c2, 1)):
                mult += 1; used[m] = True
        irreps.append((dk, mult))
    return sorted(irreps), facts, dets


def _class_rep(cycletype, base):
    perm, pos = {}, 0
    for L in cycletype:
        cyc = base[pos:pos + L]
        for i in range(L):
            perm[cyc[i]] = cyc[(i + 1) % L]
        pos += L
    return perm


def block_sizes(d=11, n=3, factor_blocks=((0, 1, 2, 3, 4), (5, 6, 7, 8, 9, 10))):
    """Exact PSD block sizes = mult of each S5 x S6 irrep in Herm(V). Returns
    (sorted [(label, mult)], sum_mult_sq, largest)."""
    dets = [tuple(c) for c in itertools.combinations(range(d), n)]
    na, nb = len(factor_blocks[0]), len(factor_blocks[1])
    base_a, base_b = list(factor_blocks[0]), list(factor_blocks[1])
    irr_a, cl_a, tab_a = char_table(na)
    irr_b, cl_b, tab_b = char_table(nb)
    G = 1
    for k in (na, nb):
        f = 1
        for x in range(2, k + 1):
            f *= x
        G *= f

    def chi_V(sig, tau):
        pm = {}; pm.update(sig); pm.update(tau); tr = 0
        for T in dets:
            if tuple(sorted(pm[x] for x in T)) == T:
                pT = [pm[x] for x in T]
                inv = sum(1 for a in range(n) for b in range(a + 1, n) if pT[a] > pT[b])
                tr += (-1) ** inv
        return tr

    chiV = {}
    for ca in cl_a:
        sig = _class_rep(ca, base_a)
        for cb in cl_b:
            chiV[(ca, cb)] = chi_V(sig, _class_rep(cb, base_b))

    out, tot2 = [], 0
    for va in irr_a:
        for vb in irr_b:
            m = Fraction(0)
            for ca in cl_a:
                sa = class_size(ca, na)
                for cb in cl_b:
                    m += sa * class_size(cb, nb) * chiV[(ca, cb)] ** 2 \
                        * tab_a[(va, ca)] * tab_b[(vb, cb)]
            m = m / G
            assert m.denominator == 1
            m = int(m)
            if m:
                out.append(((tuple(va), tuple(vb)), m)); tot2 += m * m
    out.sort(key=lambda x: -x[1])
    return out, tot2, out[0][1]


def main():
    t = time.time()
    irreps, _facts, _dets = decompose_V()
    print(f"V (3,11) 165-dim under S5 x S6: irreps (dim,mult) = {irreps}")
    print(f"  sum m_mu = {sum(m for _, m in irreps)} (expect 12); "
          f"dim End_G(V) = {sum(m * m for _, m in irreps)} (expect 18)")
    blocks, tot2, largest = block_sizes()
    print(f"reduced SOS Gram: {len(blocks)} PSD blocks, sum mult^2 = {tot2} "
          f"(expect 28884), largest block = {largest}")
    print(f"  top blocks: {[(a, b, m) for (a, b), m in blocks[:6]]}")
    print(f"done in {time.time() - t:.1f}s")


if __name__ == "__main__":
    main()
