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


def build_rho_symbolic(d: int, support, amps):
    """Exact symbolic 1-RDM of sum_t amps_t |t>, in rectangular (a + b i) form."""
    import sympy as sp

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
    # rectangular form up front: interference phases enter as exp(i*theta), and
    # sympy's charpoly block-factorization cannot order complex exp() factors
    # (raises on the symbolic <), so expand to a + b i before any determinant.
    return rho.applyfunc(lambda e: sp.expand_complex(sp.expand(e)))


def verify_exact(n: int, d: int, spectrum, support, amps):
    """Exact 1-RDM identity check for a symbolic state. True is a proof.

    Gauge-invariant: compares the characteristic polynomial of the 1-RDM to
    prod(x - lambda_i). Two robustness points matter for interference states
    whose amplitudes carry genuine phases exp(i*theta): the char poly is taken
    as the Berkowitz determinant of (x I - rho) (charpoly()'s block
    factorization tries to sort complex factors and raises), and each
    coefficient of the difference is reduced with expand_complex before the
    zero test (sp.simplify alone leaves true zeros like 1 + exp(2 i pi/3)
    unreduced, which silently rejected valid closed forms)."""
    import sympy as sp

    spectrum = [sp.Rational(Fraction(x).numerator, Fraction(x).denominator) for x in spectrum]
    rho = build_rho_symbolic(d, support, amps)
    x = sp.symbols("x")
    cp = (x * sp.eye(d) - rho).det(method="berkowitz")
    diff = sp.Poly(sp.expand(cp - sp.prod([x - lv for lv in spectrum])), x)
    return all(sp.simplify(sp.expand_complex(c)) == 0 for c in diff.all_coeffs())


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


def _active_offdiagonals(d: int, dets, ks, den, spectrum):
    """Off-diagonal 1-RDM entries forced nonzero by the support, with their exact
    Schur-Horn target magnitudes.

    Each modulus is fixed (|c_t|^2 = ks_t/den), so the diagonal occupations are
    known rationals. Wherever two orbitals p, q carry an equal occupation but the
    spectrum needs them split, the 2x2 block [[occ_p, x],[conj(x), occ_q]] must
    have two required eigenvalues lo, hi with lo + hi = occ_p + occ_q (trace) and
    lo*hi = occ_p*occ_q - |x|^2 (determinant), which forces
    |x| = sqrt(occ_p occ_q - lo*hi) exactly (an algebraic number, rational or a
    low-degree surd). This is the interference analogue of a design's diagonal
    1-RDM: the phases are no longer gauge, they must realize these off-diagonal
    magnitudes. Returns (occ, edges) with edges a list of (p, q, target_sq, terms)
    where terms are (i, j, sign): the determinant pair (i has q, j has p)
    contributing sign * sqrt(ks_i ks_j)/den * exp(i*(theta_i - theta_j))."""
    import sympy as sp

    dets = [tuple(t) for t in dets]
    index = {t: i for i, t in enumerate(dets)}
    occ = [sum(sp.Rational(ks[i], den) for i, t in enumerate(dets) if m in t)
           for m in range(d)]
    spec = sorted({sp.Rational(Fraction(v).numerator, Fraction(v).denominator)
                   for v in spectrum})
    edges = []
    for p in range(d):
        for q in range(p + 1, d):
            terms = []
            for j, t in enumerate(dets):
                if p in t and q not in t:
                    s1 = (-1) ** t.index(p)
                    tp = tuple(sorted(tuple(x for x in t if x != p) + (q,)))
                    if tp in index:
                        i = index[tp]
                        s2 = (-1) ** tp.index(q)
                        terms.append((i, j, s1 * s2))
            if not terms:
                continue
            tot = occ[p] + occ[q]
            # the block splits its (degenerate) diagonal into two spectrum values
            target_sq = None
            for lo in spec:
                hi = tot - lo
                if lo < hi and hi in spec:
                    cand = occ[p] * occ[q] - lo * hi
                    if cand > 0:
                        target_sq = cand
                        break
            edges.append((p, q, target_sq, terms))
    return occ, edges


def exactify_interference(n: int, d: int, spectrum, record):
    """Constructive Tier-B exactifier for interference vertices.

    The general exactify() recognizes each determinant's absolute (gauge-fixed)
    phase. That fails on interference corners whose phases are coupled polygon
    angles: the absolute phases are path sums of high algebraic degree, and a
    generic numeric solve lands on an arbitrary gauge, so per-phase recognition
    (and PSLQ) see only noise. This solver works in the pinned variables instead.
    The off-diagonal 1-RDM magnitudes are forced exactly by Schur-Horn on the
    degenerate occupation blocks (_active_offdiagonals); each active edge is a
    closed polygon whose sides are the fixed term moduli, so the relative phase
    is an exact arccos of a rational or a surd. Phases are propagated across edges
    (constraint propagation: an edge with a single unassigned determinant fixes
    it), and the closed form is accepted only if verify_exact certifies it, so a
    wrong branch cannot slip through. This is the same mechanism as v_B's
    cos(gamma) = 3/(4 sqrt(14)); v_B was the first visible instance, not a special
    case. Returns an EXACT record or None (caller falls through to TIER-C)."""
    import math

    import sympy as sp

    den = 1
    for xx in spectrum:
        den = den * Fraction(xx).denominator // math.gcd(den, Fraction(xx).denominator)
    dets = [tuple(s) for s, _, _ in record["support"]]
    mods = [a for _, a, _ in record["support"]]
    ks = snap_moduli(mods, den)
    if ks is None:
        return None
    occ, edges = _active_offdiagonals(d, dets, ks, den, spectrum)
    active = [e for e in edges if e[3]]
    if not active or any(e[2] is None for e in active):
        return None
    S = len(dets)

    def gmod(i, j):
        return sp.sqrt(sp.Rational(ks[i], den) * sp.Rational(ks[j], den))

    # Gauge: every determinant is real (phase 0) except one "solve" determinant
    # per active edge. Each edge is one magnitude equation |off-diagonal| = X, so
    # it pins exactly one phase; the U(1)^d orbital gauge lets the rest be zero.
    # Assign each edge a distinct solve determinant (one of its own determinants,
    # preferring one not shared with a still-unsolved edge), then propagate: an
    # edge is solvable once all but its solve determinant are fixed.
    theta = [sp.Integer(0)] * S
    solve_det = {}
    used = set()
    # how many active edges each determinant touches; a determinant private to one
    # edge can absorb that edge's constraint without perturbing any other, so the
    # data-flow order below always resolves. Prefer private determinants; edges
    # with none take a shared one and rely on iteration order.
    from collections import Counter
    edge_dets = [set(sum([[i, j] for i, j, _ in e[3]], [])) for e in active]
    glob = Counter(dd for es in edge_dets for dd in es)
    order = sorted(range(len(active)),
                   key=lambda e: min((glob[dd] for dd in edge_dets[e]), default=99))
    for e in order:
        (p, q, tsq, terms) = active[e]
        cand = [dd for (i, j, _) in terms for dd in (i, j)]
        pick = next((dd for dd in sorted(cand, key=lambda dd: glob[dd])
                     if dd not in used
                     and sum(dd in (i, j) for (i, j, _) in terms) == 1), None)
        if pick is None:
            return None
        solve_det[(p, q)] = pick
        used.add(pick)
        theta[pick] = None

    for _ in range(len(active) + 2):
        progressed = False
        for (p, q, tsq, terms) in active:
            free_det = solve_det[(p, q)]
            if theta[free_det] is not None:
                continue
            # terms not touching free_det must be fully fixed to propagate
            if any(theta[i] is None or theta[j] is None
                   for (i, j, s) in terms if free_det not in (i, j)):
                continue
            X = sp.sqrt(tsq)
            P = sum(s * gmod(i, j) * sp.exp(sp.I * (theta[i] - theta[j]))
                    for (i, j, s) in terms if free_det not in (i, j))
            P = sp.expand_complex(P)
            i, j, s = next((i, j, s) for (i, j, s) in terms if free_det in (i, j))
            g = gmod(i, j)
            Pabs2 = sp.expand_complex(sp.Abs(P) ** 2)
            rhs = sp.simplify((X ** 2 - Pabs2 - g ** 2) / (2 * s * g))
            if P == 0:
                phi = sp.Integer(0)  # single term: |x| forced, phase pure gauge
            else:
                phi = sp.arg(P) + sp.acos(sp.simplify(rhs / sp.Abs(P)))
            other = j if free_det == i else i
            theta[free_det] = sp.simplify((theta[other] + phi) if free_det == i
                                          else (theta[other] - phi))
            progressed = True
        if all(t is not None for t in theta) or not progressed:
            break
    theta = [sp.Integer(0) if t is None else t for t in theta]
    amps = [sp.sqrt(sp.Rational(ks[i], den)) * sp.exp(sp.I * theta[i]) for i in range(S)]
    if verify_exact(n, d, spectrum, dets, amps):
        return {"status": "EXACT", "weights": ks, "den": den,
                "amplitudes": [sp.srepr(a) for a in amps],
                "pretty": [str(sp.simplify(a)) for a in amps]}
    # Propagation cannot decouple edges that share every determinant (no edge has
    # a private phase to absorb its constraint). Those form a small coupled
    # system; solve it by a bounded joint search over the algebraic angles the
    # edge targets admit, on the union of the per-edge solve determinants. Each
    # candidate is still gated by verify_exact, so nothing uncertified is returned.
    import itertools

    seeds = sorted(set(solve_det.values()))
    if not seeds or len(seeds) > 3:
        return None
    cand = [sp.Integer(0), sp.pi, sp.pi / 3, -sp.pi / 3, 2 * sp.pi / 3, -2 * sp.pi / 3]
    for (p, q, tsq, terms) in active:
        gs = {sp.simplify(gmod(i, j)) for (i, j, s) in terms}
        X = sp.sqrt(tsq)
        for g in gs:
            for h in gs:
                r = sp.simplify((X ** 2 - g ** 2 - h ** 2) / (2 * g * h))
                if abs(sp.N(r)) <= 1:
                    cand += [sp.acos(r), -sp.acos(r)]
    cand = list(dict.fromkeys(cand))
    for combo in itertools.product(cand, repeat=len(seeds)):
        th = [sp.Integer(0)] * S
        for dd, val in zip(seeds, combo):
            th[dd] = val
        amps = [sp.sqrt(sp.Rational(ks[i], den)) * sp.exp(sp.I * th[i]) for i in range(S)]
        if verify_exact(n, d, spectrum, dets, amps):
            return {"status": "EXACT", "weights": ks, "den": den,
                    "amplitudes": [sp.srepr(a) for a in amps],
                    "pretty": [str(sp.simplify(a)) for a in amps]}
    return None


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
    gphases = gauge_fix_phases(support, phases, d)
    exact_amps = []
    recognized = True
    for k, th in zip(ks, gphases):
        ph = recognize_phase(th)
        if ph is None:
            recognized = False
            break
        exact_amps.append(sp.sqrt(sp.Rational(k, den)) * sp.exp(sp.I * ph))
    if recognized and verify_exact(n, d, spectrum, support, exact_amps):
        return {"status": "EXACT", "weights": ks, "den": den,
                "amplitudes": [sp.srepr(a) for a in exact_amps],
                "pretty": [str(sp.simplify(a)) for a in exact_amps]}
    # Per-phase recognition failed or did not certify. The absolute phases of an
    # interference corner are coupled polygon path-sums (high degree) even after
    # gauge fixing; solve instead in the pinned off-diagonal magnitudes.
    constructive = exactify_interference(n, d, spectrum, record)
    if constructive is not None:
        return constructive
    return {"status": "TIER-C", "reason": "interference phases not resolved"}
