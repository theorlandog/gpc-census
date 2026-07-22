"""Degenerate-block ansatz: fold kxk_degen's repeated-eigenvalue blocks into the
production family.

The 2x2 block_ansatze and the k-clique generators both mix DISTINCT eigenvalue
classes; neither produces a block that mixes a repeated eigenvalue (e.g. a 3x3
block with spectrum (15,6,6) for (3,10) v89 = (15,15,6,6,6,6,6,6,6,6)/26). This
module generates those degenerate block degree vectors (a sub-multiset with a
repeat, a Schur-Horn integer diagonal split across the block modes) and solves
each with the existing, certified machinery: enumerate one-hop-confined supports
with the shared skeleton model (so F1 symmetry reduction and F2 kernel
quotienting apply), attain a state by eigenvalue matching (attain is basis-free,
so a block-diagonal 1-RDM is fine, unlike the fixed-target phase_solve used for
2x2 blocks), and gate every hit with the exact characteristic-polynomial
identity (exactify -> verify_exact). Soundness rests entirely on that gate: a
wrong support or a spurious numeric solve simply fails verification.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd


def degenerate_block_ansatze(n: int, d: int, spectrum, block: int = 3):
    """Yield (nv, block_pairs) for every degenerate k-block on the last k modes.

    A block eigen sub-multiset E (size `block`, drawn from the spectrum, with at
    least one repeat and at least two distinct values, the degenerate case the
    2x2/clique families miss) plus an integer diagonal split (d_1,...,d_k)
    majorized by E gives a degree vector nv whose last k entries are the split
    and whose earlier entries are the remaining spectrum, sorted. block_pairs is
    every mode pair inside the block (each carries a nonzero off-diagonal).
    """
    spec = [Fraction(x) for x in spectrum]
    den = 1
    for x in spec:
        den = den * x.denominator // gcd(den, x.denominator)
    ints = [int(x * den) for x in spec]
    modes = list(range(d - block, d))
    pairs = [tuple(sorted(p)) for p in combinations(modes, 2)]
    seen_e = set()
    for E in combinations(sorted(ints, reverse=True), block):
        if len(set(E)) < 2 or len(set(E)) == block:
            continue  # need a repeat (degenerate) and not all-distinct (that is a clique)
        E = tuple(sorted(E, reverse=True))
        if E in seen_e:
            continue
        seen_e.add(E)
        rem = list(ints)
        for e in E:
            rem.remove(e)
        tot = sum(E)
        for split in _majorized_splits(E, den):
            if sum(split) != tot:
                continue
            nv = sorted(rem, reverse=True) + list(split)
            if len(nv) == d:
                yield nv, pairs


def _majorized_splits(E, den):
    """Integer diagonals (d_1..d_k), each in [1, den], with the same sum as E and
    majorized by E (necessary for a Hermitian block with spectrum E/den and that
    integer diagonal to exist, Schur-Horn). Canonical (non-increasing) only."""
    k = len(E)
    tot = sum(E)
    Esort = sorted(E, reverse=True)

    def rec(i, remaining, prev, acc):
        if i == k - 1:
            last = remaining
            if 1 <= last <= min(prev, den):
                cand = acc + [last]
                if _majorizes(Esort, cand):
                    yield tuple(cand)
            return
        hi = min(prev, den, remaining - (k - 1 - i))
        for v in range(hi, 0, -1):
            yield from rec(i + 1, remaining - v, v, acc + [v])

    yield from rec(0, tot, den, [])


def _majorizes(E, diag):
    """E majorizes diag: partial sums of sorted-desc E dominate those of diag."""
    ds = sorted(diag, reverse=True)
    se = sc = 0
    for i in range(len(E)):
        se += E[i]
        sc += ds[i]
        if sc > se:
            return False
    return True


def solve_degenerate_vertex(n: int, d: int, spectrum, budget: float = 600.0,
                            block: int = 3, max_card: int = 24, limit: int = 300):
    """Solve a vertex with the degenerate-block ansatz, certified by verify_exact.

    Returns an OK record with an exact closed form (record["exact"]) on success,
    else a FAIL/EXHAUSTED record. Reuses the shared skeleton enumeration (F1/F2
    apply) and attain + exactify; the exact identity is the only thing that
    certifies, so this can never emit an unsound state.
    """
    import time

    from . import states as S
    from .exactify import exactify

    built = S._build(d, n)
    dets_all = built[0]
    spec = [Fraction(x) for x in spectrum]
    den = 1
    for x in spec:
        den = den * x.denominator // gcd(den, x.denominator)

    t0 = time.time()
    n_supports = 0
    exhausted = True
    for nv, pairs in degenerate_block_ansatze(n, d, spectrum, block=block):
        if time.time() - t0 > budget:
            exhausted = False
            break
        mincard = S.min_support_cardinality(n, d, spectrum, nv=nv,
                                             require_hop_pairs=pairs, time_cap=2)
        if mincard is None:
            continue
        cap = min(max_card, len(dets_all))
        if cap < mincard:
            continue
        _, _den, sols = S.enumerate_weight_vectors(
            n, d, spectrum, mincard, nv=nv, require_hop_pairs=pairs,
            forbid_offtarget=True, max_cardinality=cap, limit=limit)
        tried: set = set()
        for w in sorted(sols, key=lambda w: sum(1 for k in w if k)):
            if time.time() - t0 > budget:
                exhausted = False
                break
            sup = tuple(i for i, k in enumerate(w) if k)
            key = (sup, tuple(sorted(S.rigid_class_products(
                [tuple(dets_all[i]) for i in sup], [w[i] for i in sup]).items())))
            if key in tried:
                continue
            tried.add(key)
            n_supports += 1
            mask = _mask(w, len(dets_all))
            psi, res, _ = S.attain(n, d, spectrum, mask=mask, outer=200, tries=3,
                                   _built=built)
            if psi is None or res >= 1e-12:
                continue
            rec = _record(n, d, spectrum, dets_all, w, psi, den, nv, pairs)
            ex = exactify(n, d, spectrum, rec)
            if ex.get("status") == "EXACT":
                rec["exact"] = ex
                rec["degenerate_block"] = {"nv": nv, "pairs": [list(p) for p in pairs]}
                return rec
    status = "EXHAUSTED" if exhausted else "TIMEOUT"
    return {"status": "FAIL", "reason": f"degenerate-block {status}",
            "exhausted": exhausted, "supports_tried": n_supports,
            "secs": round(time.time() - t0, 1)}


def _mask(w, m):
    import numpy as np

    return np.array([1.0 if (i < len(w) and w[i]) else 0.0 for i in range(m)])


def _record(n, d, spectrum, dets_all, w, psi, den, nv, pairs):
    import numpy as np

    sup = [i for i, k in enumerate(w) if k]
    j0 = max(sup, key=lambda i: abs(psi[i]))
    psi = psi * np.exp(-1j * np.angle(psi[j0]))
    return {"status": "OK", "residual": 0.0, "support_size": len(sup),
            "weights": [w[i] for i in sup], "den": den,
            "support": [[list(dets_all[i]), float(abs(psi[i])),
                         float(np.angle(psi[i]))] for i in sup]}
