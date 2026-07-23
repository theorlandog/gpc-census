#!/usr/bin/env python3
"""Fiber-curve census for kdim-1 single-class deforming interference states.

For each family the fiber curve is w^2 = Q(t), Q = product of the touched
class's weight-product factors. CORRECT reduction (an earlier radical-based
pass was WRONG -- see docs/RESEARCH.md retraction): divide Q by its square
factors, i.e. keep each linear factor to its multiplicity MOD 2. Genus is 1
iff the reduced polynomial has degree 3 or 4 (3 or 4 odd-multiplicity roots),
else 0. j-invariants are computed only for genuine genus-1 models.

Run: python scripts/fiber_curves.py  (from repo root; reads results/data/states.jsonl)
"""
import json
from collections import defaultdict

import sympy as sp

t = sp.Symbol('t')


def reduced_model(Q):
    """Return (genus, reduced_poly) with Q reduced modulo square factors."""
    P = sp.Poly(sp.expand(Q), t)
    const, factors = sp.factor_list(P.as_expr())
    red = sp.Integer(1) * const
    for f, m in factors:
        if m % 2 == 1:
            red *= f
    R = sp.Poly(sp.expand(red), t)
    dg = R.degree()
    return (1 if dg in (3, 4) else 0), R


def j_invariant(R):
    dg = R.degree()
    if dg == 4:
        a4, a3, a2, a1, a0 = (R.coeff_monomial(t**k) for k in (4, 3, 2, 1, 0))
        inv_i = 12 * a4 * a0 - 3 * a3 * a1 + a2**2
        inv_j = 72 * a4 * a2 * a0 - 27 * a4 * a1**2 - 27 * a3**2 * a0 + 9 * a3 * a2 * a1 - 2 * a2**3
        disc = 4 * inv_i**3 - inv_j**2
        if disc == 0:
            return None
        return sp.nsimplify(sp.Rational(6912) * inv_i**3 / disc)
    if dg == 3:
        a, b, c, d = (R.coeff_monomial(t**k) for k in (3, 2, 1, 0))
        big_b, big_c, big_d = b / a, c / a, d / a
        b2, b4, b6 = 4 * big_b, 2 * big_c, 4 * big_d
        c4 = b2**2 - 24 * b4
        c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
        disc = (c4**3 - c6**2) / 1728
        if disc == 0:
            return None
        return sp.nsimplify(sp.cancel(c4**3 / disc))
    return None


def families(states_path='results/data/states.jsonl'):
    for line in open(states_path):
        d = json.loads(line)
        cf = d.get('closed_form')
        if not cf or d['classified'] != 'INTERFERENCE':
            continue
        dets = [tuple(sorted(T)) for T in cf['support_dets']]
        n = len(dets)
        dd = max(max(T) for T in dets) + 1
        matrix = sp.Matrix([[1 if m in T else 0 for m in range(dd)] for T in dets])
        ker = matrix.T.nullspace()
        if len(ker) != 1:
            continue
        v = ker[0]
        v = v * sp.lcm([sp.fraction(sp.nsimplify(q))[1] for q in v])
        v = [int(q) for q in v]
        cls = defaultdict(list)
        for i in range(n):
            for j in range(i + 1, n):
                set_a, set_b = set(dets[i]), set(dets[j])
                if len(set_a & set_b) == len(dets[i]) - 1:
                    cls[tuple(sorted(set_a ^ set_b))].append((i, j))
        touched = [(p, tm) for p, tm in cls.items() if any(v[i] or v[j] for i, j in tm)]
        if len(touched) != 1 or len(touched[0][1]) != 2:
            continue
        (i1, j1), (i2, j2) = touched[0][1]
        k0 = cf['weights']
        q = (k0[i1] + t * v[i1]) * (k0[j1] + t * v[j1]) * \
            (k0[i2] + t * v[i2]) * (k0[j2] + t * v[j2])
        yield (d['system'], d['index']), q


if __name__ == '__main__':
    genus_count = defaultdict(int)
    js = defaultdict(list)
    for key, q in families():
        g, r = reduced_model(q)
        genus_count[g] += 1
        if g == 1:
            js[sp.sstr(j_invariant(r))].append(key)
    print('genus census (correct square-reduction):', dict(genus_count))
    print(f'{len(js)} distinct j-invariants among genuine genus-1 families:')
    for jv, ks in sorted(js.items(), key=lambda kv: -len(kv[1])):
        print(f'  j = {jv[:44]:44s} x{len(ks)}  {ks[:6]}')
