"""Exact certificates for census verdicts. Zero dependencies beyond scipy
for the Farkas search; every verification is exact rational arithmetic.

Design verdicts are certified by their witness (exact sum checks).
Interference verdicts are certified by a rational Farkas vector y with
y . b > 0 and y^T A <= 0 componentwise, proving no nonnegative real
weighting exists (which also rules out integer weightings).
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations


def _system(n: int, d: int, spectrum):
    dets = list(combinations(range(d), n))
    spectrum = [Fraction(x) for x in spectrum]
    # rows of A: per-mode sums then normalization; b: spectrum then 1
    a = [[Fraction(1) if m in t else Fraction(0) for t in dets] for m in range(d)]
    a.append([Fraction(1)] * len(dets))
    b = spectrum + [Fraction(1)]
    return dets, a, b


def verify_design(n: int, d: int, spectrum, witness) -> bool:
    """Exactly verify a weighted-design witness (integer or rational)."""
    dets, a, b = _system(n, d, spectrum)
    w = [Fraction(x).limit_denominator(10**9) if not isinstance(x, int) else Fraction(x)
         for x in witness]
    total = sum(w)
    if total == 0:
        return False
    w = [x / total for x in w]
    if any(x < 0 for x in w):
        return False
    for row, target in zip(a, b):
        if sum(r * x for r, x in zip(row, w)) != target:
            return False
    return True


def farkas_interference(n: int, d: int, spectrum, max_den: int = 10**6):
    """Search (float) and verify (exact) a Farkas certificate of infeasibility.

    Returns the rational certificate vector on success, None on failure.
    """
    import numpy as np
    from scipy.optimize import linprog

    dets, a, b = _system(n, d, spectrum)
    rows = len(a)
    af = np.array([[float(x) for x in row] for row in a])
    bf = np.array([float(x) for x in b])
    # find y maximizing y.b with y^T A <= 0, bounded
    res = linprog(-bf, A_ub=af.T, b_ub=np.zeros(af.shape[1]),
                  bounds=[(-1, 1)] * rows, method="highs")
    if not res.success or -res.fun <= 1e-9:
        return None
    y = [Fraction(v).limit_denominator(max_den) for v in res.x]
    if sum(yy * bb for yy, bb in zip(y, b)) <= 0:
        return None
    for j in range(len(dets)):
        if sum(y[i] * a[i][j] for i in range(rows)) > 0:
            return None
    return y
