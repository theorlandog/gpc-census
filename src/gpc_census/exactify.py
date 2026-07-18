"""Exactification of numerically attained states into closed forms.

Moduli: |c_t|^2 must equal k_t/den (integer design weights with phases),
so squared amplitudes snap to the natural-denominator grid exactly.
Phases: relative phases are snapped to a small algebraic lattice
(rational multiples of pi, then cosines of the form p*sqrt(q)/r) and the
reconstructed symbolic state is verified by an exact 1-RDM identity
check in sympy. A returned record is a certificate; recognition failures
fall through labeled for hand analysis.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations


def snap_moduli(psi_abs, den: int):
    """Snap |c|^2 to k/den exactly. Returns k list or None if off-grid."""
    ks = []
    for a in psi_abs:
        k = round(a * a * den)
        if abs(a * a * den - k) > 1e-6:
            return None
        ks.append(int(k))
    return ks


def recognize_cos(x: float, max_r: int = 60, max_q: int = 60):
    """Recognize x as p*sqrt(q)/r with small integers. Returns sympy expr or None."""
    import sympy as sp

    for q in range(1, max_q + 1):
        sq = q ** 0.5
        fr = Fraction(x / sq).limit_denominator(max_r)
        if fr.denominator <= max_r and abs(float(fr) * sq - x) < 1e-9:
            return sp.Rational(fr.numerator, fr.denominator) * sp.sqrt(q)
    return None


def recognize_phase(theta: float):
    """Recognize a phase angle: rational multiple of pi, else via its cosine."""
    import sympy as sp

    fr = Fraction(theta / 3.141592653589793).limit_denominator(24)
    if abs(float(fr) * 3.141592653589793 - theta) < 1e-9:
        return sp.pi * sp.Rational(fr.numerator, fr.denominator)
    c = recognize_cos(float(__import__("math").cos(theta)))
    if c is not None:
        s = recognize_cos(float(__import__("math").sin(theta)))
        if s is not None and sp.simplify(c**2 + s**2 - 1) == 0:
            return sp.atan2(s, c)
    return None


def verify_exact(n: int, d: int, spectrum, support, amps):
    """Exact 1-RDM identity check for a symbolic state. True is a proof."""
    import sympy as sp

    spectrum = [sp.Rational(Fraction(x).numerator, Fraction(x).denominator) for x in spectrum]
    dets = [tuple(s) for s, _ in zip((tuple(x) for x in support), amps)]
    amap = dict(zip(dets, amps))
    rho = sp.zeros(d, d)
    for t, ct in amap.items():
        for mp in t:
            s1 = (-1) ** t.index(mp)
            t2 = tuple(x for x in t if x != mp)
            for m in range(d):
                if m in t2:
                    continue
                tp = tuple(sorted(t2 + (m,)))
                if tp not in amap:
                    continue
                s2 = (-1) ** tp.index(m)
                rho[m, mp] += s1 * s2 * sp.conjugate(amap[tp]) * ct
    # spectrum check: char poly of rho equals prod (x - lambda_i)
    x = sp.symbols("x")
    p1 = sp.expand(rho.charpoly(x).as_expr())
    p2 = sp.expand(sp.prod([x - lv for lv in spectrum]))
    return sp.simplify(p1 - p2) == 0


def exactify(n: int, d: int, spectrum, record):
    """Tier B: numerical solve_vertex record -> exact certified state or labeled failure."""
    import math

    import sympy as sp

    den = 1
    for xx in spectrum:
        den = den * Fraction(xx).denominator // math.gcd(den, Fraction(xx).denominator)
    support = [tuple(s) for s, _, _ in record["support"]]
    mods = [a for _, a, _ in record["support"]]
    phases = [p for _, _, p in record["support"]]
    ks = snap_moduli(mods, den)
    if ks is None:
        return {"status": "TIER-C", "reason": "moduli off natural grid"}
    exact_amps = []
    for k, th in zip(ks, phases):
        ph = recognize_phase(th)
        if ph is None:
            return {"status": "TIER-C", "reason": f"unrecognized phase {th}"}
        exact_amps.append(sp.sqrt(sp.Rational(k, den)) * sp.exp(sp.I * ph))
    if not verify_exact(n, d, spectrum, support, exact_amps):
        return {"status": "TIER-C", "reason": "exact verification failed"}
    return {"status": "EXACT", "weights": ks, "den": den,
            "amplitudes": [sp.srepr(a) for a in exact_amps],
            "pretty": [str(sp.simplify(a)) for a in exact_amps]}
