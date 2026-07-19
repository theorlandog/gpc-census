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


def recognize_algebraic(x: float, maxdeg: int = 4, maxcoeff: int = 10 ** 5,
                        dps: int = 30):
    """Recognize a real number as a low-degree algebraic number via PSLQ on its
    powers: an integer relation among 1, x, ..., x^deg is a minimal polynomial,
    whose real root matching x is the exact value (a sympy radical or CRootOf).
    Loose by construction (float precision limits the reliable degree), so it is
    only ever a proposal; verify_exact is the gate that accepts or rejects it.
    Extends recognize_cos (p*sqrt(q)/r, the b = 0 degree-2 case) to general
    low-degree algebraics."""
    import mpmath as mp
    import sympy as sp

    mp.mp.dps = dps
    xm = mp.mpf(repr(x))
    for deg in range(2, maxdeg + 1):
        rel = mp.pslq([xm ** i for i in range(deg + 1)], maxcoeff=maxcoeff,
                      maxsteps=10 ** 5)
        if not rel or rel[-1] == 0:
            continue
        X = sp.Symbol("X")
        poly = sp.Poly([int(c) for c in reversed(rel)], X)
        try:
            roots = poly.all_roots()
        except (sp.PolynomialError, NotImplementedError):
            continue
        for r in roots:
            if r.is_real and abs(float(r.evalf(20)) - x) < 1e-8:
                return sp.nsimplify(r) if r.is_number else r
    return None


def recognize_phase(theta: float):
    """Recognize a phase angle: rational multiple of pi, then via its cosine on
    the p*sqrt(q)/r lattice, then via low-degree-algebraic cosine and sine
    (PSLQ). The algebraic branch is a proposal only; the caller's verify_exact
    accepts or rejects it, so a spurious PSLQ hit cannot certify a wrong state."""
    import math

    import sympy as sp

    fr = Fraction(theta / 3.141592653589793).limit_denominator(24)
    if abs(float(fr) * 3.141592653589793 - theta) < 1e-9:
        return sp.pi * sp.Rational(fr.numerator, fr.denominator)
    c = recognize_cos(float(math.cos(theta)))
    if c is not None:
        s = recognize_cos(float(math.sin(theta)))
        if s is not None and sp.simplify(c ** 2 + s ** 2 - 1) == 0:
            return sp.atan2(s, c)
    ca = recognize_algebraic(float(math.cos(theta)))
    if ca is not None:
        sa = recognize_algebraic(float(math.sin(theta)))
        if sa is not None and sp.simplify(ca ** 2 + sa ** 2 - 1) == 0:
            return sp.atan2(sa, ca)
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


def gauge_fix_phases(support, phases, d: int):
    """Remove the single-particle U(1)^d phase freedom before recognition.

    The state is defined only up to c_t -> c_t * prod_{m in t} exp(i phi_m)
    (an orbital phase rotation a_m -> exp(i phi_m) a_m), which conjugates the
    1-RDM by a diagonal unitary and so preserves its spectrum: verify_exact is
    invariant under it. A generic numerical solve lands on an arbitrary point
    of this orbit with scrambled phases; projecting them orthogonal to the
    gauge orbit (least squares with 2pi wrapping) leaves the gauge-invariant
    interference phases, which are the ones a closed form must express."""
    import math

    import numpy as np

    theta = np.array([float(p) for p in phases])
    a = np.array([[1.0 if m in t else 0.0 for m in range(d)] for t in support])
    phi = np.zeros(d)
    for _ in range(300):
        r = (theta + a @ phi + math.pi) % (2 * math.pi) - math.pi
        dphi, *_ = np.linalg.lstsq(a, -r, rcond=None)
        phi += dphi
        if np.linalg.norm(dphi) < 1e-15:
            break
    return list((theta + a @ phi + math.pi) % (2 * math.pi) - math.pi)


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
    # canonicalise the gauge before recognising phases; verify_exact below is
    # gauge-invariant, so the certificate stands in this frame
    phases = gauge_fix_phases(support, phases, d)
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
