#!/usr/bin/env python3
"""Stage-0 INNER hull generator for fermionic moment polytopes (GPC polytopes),
via Klyachko's practical algorithm, in exact arithmetic. No Sage, no lrcalc.

MATH
----
Klyachko (Altunbulak thesis, "Practical algorithm", thesis lines ~977-1076):
every irreducible GL(d) component S_lambda(C^d) appearing in S^m(wedge^N C^d)
certifies lambda/m as an ATTAINABLE occupation spectrum of an N-fermion pure
state (a point of the moment polytope Delta(N,d)). The multiplicity is

    mult(lambda) = < s_lambda , h_m[e_N] >

computed here by expanding the plethysm h_m[e_N] in the power-sum basis and
pairing with s_lambda via the Murnaghan-Nakayama rule:

    h_m[f]    = sum_{mu |- m} (1/z_mu) prod_i p_{mu_i}[f]
    p_k[e_N]  = e_N with the substitution p_j -> p_{jk}
    e_N (in p) by Newton's identity: e_n = (1/n) sum_{k=1}^{n} (-1)^{k-1} p_k e_{n-k}
    < f , s_lambda > = sum_mu c_mu * chi^lambda(mu)   for f = sum_mu c_mu p_mu

All coefficients are exact Fractions; characters are exact integers (MN rule
via beta-numbers). Every emitted point is a CERTIFIED attainable spectrum.

SCOPE AND CAVEATS (read before trusting output)
-----------------------------------------------
1. INNER only. This generates attainable points; it says nothing about
   validity of inequalities. The outer side (facet certification) is Stage 1
   (Theorem 3.2.1 / Schubert coefficient test) and is NOT here.
2. Convergence: a true vertex with denominator q generally requires m = q (or
   a multiple) before it appears. Validation record from the gpc-census repo:
     - (3,6), M<=6 : exact-verified extreme points == the 4 known vertices. PASS
     - (3,7), M<=7 : 9/10 known vertices present in the cloud;
       (3,7), M<=10: 10/10 present. (One vertex needs m in {8,9,10}.)
   So run M at least up to the largest expected vertex denominator, and treat
   anything labeled "cloud-extreme" as a CANDIDATE, not a certified vertex.
3. "Cloud-extreme" means: extreme point of conv(all generated points), with an
   exact separating-functional certificate (numeric LP proposes the functional,
   exact rational arithmetic verifies strict separation). A cloud-extreme point
   IS a certified member of the true polytope (it is attainable); whether it is
   a VERTEX of the true polytope is only guaranteed once the cloud has
   converged. Face-embedded vertices from lower ranks (pads: {lambda_d=0} face;
   frozen-core lifts: {lambda_1=1} face) ARE rigorously true vertices and
   should be merged in by the caller (see gpc-census transport tooling).
4. Cost: h_m[e_N] has p-support that grows with m; MN is memoized but character
   degree is N*m. At (3,7), M=10 took ~26 s and M=12 ~145 s on one core.
   The lambda loop is embarrassingly parallel (multiprocessing over lam).
5. The lambda_1 <= m filter encodes occupation <= 1 and is required.

USAGE
-----
    python plethysm_inner_hull.py --validate            # run the (3,6)/(3,7) gates
    python plethysm_inner_hull.py -N 3 -d 12 -M 6 -o out.json
    python plethysm_inner_hull.py -N 4 -d 11 -M 5 -o out.json

Output JSON: {"points": [{spectrum, m, mult}], "extreme": [...]} with all
spectra as exact fraction strings, sorted decreasing.
"""
from __future__ import annotations
import argparse
import json
from fractions import Fraction as F
from functools import lru_cache

# ----------------------------------------------------------------------------
# partitions and symmetric-function plumbing (p-basis dicts: key = partition
# as a decreasing tuple of part sizes, value = Fraction coefficient)
# ----------------------------------------------------------------------------

def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def zmu(mu):
    from collections import Counter
    z = 1
    for part, m in Counter(mu).items():
        z *= part ** m
        for i in range(1, m + 1):
            z *= i
    return z


def pmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(sorted(ka + kb, reverse=True))
            out[k] = out.get(k, F(0)) + va * vb
    return out


def padd(a, b, coeff=F(1)):
    for k, v in b.items():
        a[k] = a.get(k, F(0)) + coeff * v
    return a


@lru_cache(maxsize=None)
def e_in_p(n):
    """e_n expanded in the power-sum basis, exact. Newton:
    e_n = (1/n) * sum_{k=1}^{n} (-1)^{k-1} p_k e_{n-k}."""
    if n == 0:
        return (((), F(1)),)
    total = {}
    for k in range(1, n + 1):
        ek = dict(e_in_p(n - k))
        term = pmul({(k,): F(1)}, ek)
        padd(total, term, F((-1) ** (k - 1), n))
    return tuple(sorted(total.items()))


def pleth_pk(k, f_items):
    """p_k[f] for f given in p-basis: substitute p_j -> p_{jk}."""
    return {tuple(sorted((j * k for j in key), reverse=True)): v
            for key, v in f_items}


def hm_eN(m, N):
    """h_m[e_N] in the p-basis, exact."""
    eN = e_in_p(N)
    total = {}
    for mu in partitions(m):
        term = {(): F(1)}
        for part in mu:
            term = pmul(term, pleth_pk(part, eN))
        padd(total, term, F(1, zmu(mu)))
    return total


# ----------------------------------------------------------------------------
# Murnaghan-Nakayama (exact characters via beta-numbers)
# ----------------------------------------------------------------------------

@lru_cache(maxsize=None)
def mn(lam, mu):
    """chi^lam(mu). lam, mu decreasing tuples, sum(lam) == sum(mu)."""
    if not mu:
        return 1 if sum(lam) == 0 else 0
    k, rest = mu[0], mu[1:]
    lam = tuple(x for x in lam if x > 0)
    n = len(lam)
    if n == 0:
        return 0
    beta = [lam[i] + (n - 1 - i) for i in range(n)]
    bs = set(beta)
    total = 0
    for b in beta:
        if b - k >= 0 and (b - k) not in bs:
            ht = sum(1 for x in beta if b - k < x < b)
            nb = sorted((x if x != b else b - k for x in beta), reverse=True)
            newlam = tuple(x for x in (nb[i] - (n - 1 - i) for i in range(n))
                           if x > 0)
            total += (-1) ** ht * mn(newlam, rest)
    return total


# ----------------------------------------------------------------------------
# inner cloud and exact extreme-point certification
# ----------------------------------------------------------------------------

def inner_spectra(N, d, Mmax, verbose=False):
    """All certified attainable spectra lambda/m, m <= Mmax.
    Returns {spectrum_tuple_of_Fractions: (m, lam, multiplicity)}."""
    pts = {}
    for m in range(1, Mmax + 1):
        f = hm_eN(m, N)
        f_items = list(f.items())
        for lam in partitions(N * m, maxpart=m):      # lambda_1 <= m: occ <= 1
            if len(lam) > d:
                continue
            mult = 0
            for mu, c in f_items:
                ch = mn(lam, tuple(sorted(mu, reverse=True)))
                if ch:
                    mult += c * ch
            if mult > 0:
                spec = tuple([F(x, m) for x in lam] + [F(0)] * (d - len(lam)))
                if spec not in pts:
                    pts[spec] = (m, lam, int(mult))
        if verbose:
            print(f"  m={m}: cumulative points {len(pts)}")
    return pts


def extreme_points(pts):
    """Exact-certificate extreme points of conv(pts).
    Strategy: exact rref finds an affine coordinate chart; Qhull proposes the
    candidate list; a numeric LP proposes a separating functional per
    candidate; exact Fraction arithmetic verifies STRICT separation. Only
    exactly-verified candidates are returned (sound; possibly incomplete if a
    numeric step fails, which is reported)."""
    import numpy as np
    from scipy.spatial import ConvexHull
    from scipy.optimize import linprog

    keys = list(pts.keys())
    P = [list(s) for s in keys]
    base = P[0]
    diffs = [[a - b for a, b in zip(p, base)] for p in P[1:]]
    # exact rref for pivot columns (affine chart)
    M = [row[:] for row in diffs]
    rows, cols = len(M), len(M[0]) if M else 0
    piv, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f2 = M[i][c]
                M[i] = [x - f2 * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == rows:
            break
    dim = len(piv)
    Q = np.array([[float(p[c] - base[c]) for c in piv] for p in P])
    if len(P) > dim + 1:
        cand = sorted(set(ConvexHull(Q, qhull_options="Qx").vertices))
    else:
        cand = list(range(len(P)))

    def val(w, j):
        return sum(w[t] * (P[j][piv[t]] - base[piv[t]]) for t in range(dim))

    verified, unverified = [], []
    for i in cand:
        others = [j for j in cand if j != i]
        A = Q[others] - Q[i]
        res = linprog(np.zeros(dim), A_ub=A, b_ub=-np.ones(len(others)),
                      bounds=[(None, None)] * dim, method="highs")
        if not res.success:
            unverified.append(keys[i])
            continue
        w = [F(x).limit_denominator(10 ** 6) for x in res.x]
        fi = val(w, i)
        if all(val(w, j) < fi for j in others):
            verified.append(keys[i])
        else:
            unverified.append(keys[i])
    return dim, verified, unverified


# ----------------------------------------------------------------------------
# validation gates (run these after ANY edit; house rule: test-first)
# ----------------------------------------------------------------------------

def validate():
    ok = True
    # (3,6): exact vertex-set reproduction at M<=6
    known_36 = {
        (F(1), F(1), F(1), F(0), F(0), F(0)),
        (F(1), F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(0)),
        (F(3, 4), F(3, 4), F(1, 2), F(1, 2), F(1, 4), F(1, 4)),
        (F(1, 2),) * 6,
    }
    pts = inner_spectra(3, 6, 6)
    dim, vx, unv = extreme_points(pts)
    match = set(vx) == known_36 and not unv
    print(f"(3,6) M<=6: {len(pts)} pts, dim {dim}, extreme {len(vx)}, "
          f"exact match to known vertex set: {match}")
    ok &= match
    # (3,7): all ten known vertices present in the cloud by M<=10
    known_37 = {
        # true (3,7) vertex list (gpc-census vertices_3_7.json)
        (F(1), F(1), F(1), F(0), F(0), F(0), F(0)),
        (F(1), F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(0), F(0)),
        (F(1), F(1, 3), F(1, 3), F(1, 3), F(1, 3), F(1, 3), F(1, 3)),
        (F(3, 4), F(3, 4), F(1, 2), F(1, 2), F(1, 4), F(1, 4), F(0)),
        (F(5, 7), F(5, 7), F(3, 7), F(3, 7), F(3, 7), F(1, 7), F(1, 7)),
        (F(2, 3), F(2, 3), F(1, 3), F(1, 3), F(1, 3), F(1, 3), F(1, 3)),
        (F(3, 5), F(3, 5), F(3, 5), F(3, 5), F(1, 5), F(1, 5), F(1, 5)),
        (F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(0)),
        (F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(1, 4), F(1, 4)),
        (F(3, 7),) * 7,
    }
    pts7 = inner_spectra(3, 7, 10)
    present = known_37 & set(pts7)
    print(f"(3,7) M<=10: {len(pts7)} pts, known vertices present: "
          f"{len(present)}/10  ({'PASS' if len(present) == 10 else 'FAIL'})")
    ok &= (len(present) == 10)
    # generic e_N sanity: N=2 components must all be evenly paired
    pts2 = inner_spectra(2, 6, 5)
    from collections import Counter
    paired = all(all(m % 2 == 0 for v, m in Counter(
        x for x in s if x != 0).items()) for s in pts2)
    print(f"(2,6) M<=5: {len(pts2)} pts, all evenly paired "
          f"(N=2 theorem): {paired}")
    ok &= paired
    print("VALIDATION", "PASS" if ok else "FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-N", type=int, default=3)
    ap.add_argument("-d", type=int, default=12)
    ap.add_argument("-M", type=int, default=6)
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()
    if args.validate:
        raise SystemExit(0 if validate() else 1)
    pts = inner_spectra(args.N, args.d, args.M, verbose=args.verbose)
    dim, vx, unv = extreme_points(pts)
    print(f"({args.N},{args.d}) M<={args.M}: {len(pts)} certified attainable "
          f"spectra; dim {dim}; cloud-extreme (exact certificate): {len(vx)}"
          + (f"; UNVERIFIED candidates: {len(unv)}" if unv else ""))
    if args.out:
        json.dump({
            "system": f"({args.N},{args.d})",
            "stage": f"Stage-0 INNER, M<={args.M} (plethysm h_m[e_N])",
            "note": ("every point is a certified attainable spectrum; "
                     "cloud-extremality is NOT true-polytope vertexhood "
                     "until the cloud converges (M >= max vertex "
                     "denominator); merge face-embedded vertices from "
                     "lower ranks separately"),
            "points": [{"spectrum": [str(x) for x in s],
                        "m": pts[s][0], "mult": pts[s][2]} for s in pts],
            "extreme": [[str(x) for x in s] for s in sorted(vx, reverse=True)],
            "extreme_unverified": [[str(x) for x in s] for s in unv],
        }, open(args.out, "w"), indent=1)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
