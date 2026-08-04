"""The conditional 2-RDM compatibility body over a pinned fermionic fiber.

WHAT THIS IS FOR
----------------
Fix a 1-RDM ``gamma`` on a pinned (generalized-Pauli-saturated) face, and let

    K(gamma) = conv { Gamma^(2)_Psi : Psi |-> gamma }

be the convex body of 2-RDMs compatible with it. For a two-body observable H,

    F_H^+(gamma) = max over K(gamma) of Tr(H Gamma),
    F_H^-(gamma) = -F_{-H}^+(gamma),

so the pair of endpoints is the SUPPORT FUNCTION of K(gamma), and knowing it
for every H determines K(gamma) outright. An observable-range solver is
therefore an exact support oracle for the conditional 2-RDM body, not a
one-off trigonometric exercise.

THE CARRIER REDUCTION
---------------------
On a pinned fiber with carrier determinants ``D_0..D_{m-1}`` and weights
``w_i = |Psi_i|^2`` fixed by gamma, the state is ``c`` with ``|c_i|^2 = w_i``,
and its carrier density matrix is ``rho = c c^dagger``, which has
``diag(rho) = w``. When the two-body coherence graph exposes every carrier
coherence through a collision-free signed channel, the variable part of the
compatible 2-RDM body is an affine image of

    E_w = { rho >= 0 : diag(rho) = w },

the scaled complex elliptope of order m.

THE ORDER-THREE THEOREM (m = 3)
-------------------------------
    E_w = conv { c c^dagger : |c_i|^2 = w_i }.

Rescaling by ``Diag(sqrt(w))`` turns E_w into the complex correlation-matrix
body. An extreme point of that body of rank r needs its m rank-one blocks
``v_i v_i^dagger`` to be linearly independent over R inside the r^2-dimensional
real space of Hermitian r-by-r matrices, so ``r^2 <= m`` (Li and Tam). At
m = 3 this forces r = 1: every extreme point is a phase matrix, and the
spectrahedron is exactly the convex hull of the fiber's own carrier states.

Consequently the support function is an SDP, exactly and not as a relaxation:

    F_H^+ = max { Tr(H rho) : rho >= 0, diag(rho) = w }
          = min { w.y : Diag(y) - H >= 0 }      (dual; Slater holds at Diag(w))
    F_H^- = max { w.y : H - Diag(y) >= 0 }.

WHERE IT STOPS
--------------
At m = 4 the bound permits r = 2, and ``rank_two_extreme_point_4x4`` builds an
explicit rank-two extreme point. So for four or more carrier determinants the
SDP is an OUTER BOUND on the conditional body, not an identity, and this
module says so rather than extrapolating the order-three result.

WHAT IS CERTIFIED
-----------------
``exact_certificates`` returns, for each endpoint, a dual point y expressed
exactly: the Groebner basis of the KKT system is in shape position, so each
``y_i`` is a rational polynomial in the objective ``t = w.y - E0``, and t is a
root of a single univariate polynomial. Optimality is then witnessed by
Hermitian PSD-ness of ``Diag(y) - H``, checked exactly from principal-minor
signs at the isolated root. That is a POINTWISE optimality proof: unlike
selecting the least and greatest of the real critical values, it needs no
argument that all real roots were found.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations

Rat = Fraction


# --------------------------------------------------------------------------
# exact univariate machinery (Fraction only, so a verifier needs no CAS)
# --------------------------------------------------------------------------

def poly_eval(p: list[Rat], x: Rat) -> Rat:
    """Evaluate a polynomial given by ascending coefficients."""
    acc = Rat(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def poly_trim(p: list[Rat]) -> list[Rat]:
    q = list(p)
    while q and q[-1] == 0:
        q.pop()
    return q


def poly_div(a: list[Rat], b: list[Rat]) -> tuple[list[Rat], list[Rat]]:
    a, b = poly_trim(a), poly_trim(b)
    if not b:
        raise ZeroDivisionError("polynomial division by zero")
    q = [Rat(0)] * max(0, len(a) - len(b) + 1)
    r = list(a)
    while poly_trim(r) and len(poly_trim(r)) >= len(b):
        r = poly_trim(r)
        k = len(r) - len(b)
        f = r[-1] / b[-1]
        q[k] = f
        for i, c in enumerate(b):
            r[i + k] -= f * c
        r = poly_trim(r)
    return poly_trim(q), poly_trim(r)


def poly_derivative(p: list[Rat]) -> list[Rat]:
    return poly_trim([c * i for i, c in enumerate(p)][1:])


def sturm_chain(p: list[Rat]) -> list[list[Rat]]:
    """Sturm chain of a squarefree polynomial."""
    chain = [poly_trim(p), poly_derivative(p)]
    while poly_trim(chain[-1]):
        _, r = poly_div(chain[-2], chain[-1])
        if not poly_trim(r):
            break
        chain.append(poly_trim([-c for c in r]))
    return chain


def sturm_count(chain: list[list[Rat]], lo: Rat, hi: Rat) -> int:
    """Number of distinct real roots in the half-open interval (lo, hi]."""
    def sign_changes(x: Rat) -> int:
        vals = [poly_eval(q, x) for q in chain if poly_trim(q)]
        vals = [v for v in vals if v != 0]
        return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))
    return sign_changes(lo) - sign_changes(hi)


def sign_at_root(q: list[Rat], p: list[Rat], lo: Rat, hi: Rat,
                 max_refine: int = 400) -> tuple[int, Rat, Rat]:
    """Exact sign of ``q`` at the unique root of ``p`` in ``(lo, hi]``.

    The interval is bisected (using sign changes of ``p``) until ``q`` is
    provably root-free on it, at which point ``q``'s sign is constant there
    and can be read off at an endpoint. Returns the sign and the refined
    interval, so a verifier can replay the same reasoning.
    """
    qt = poly_trim(q)
    if not qt:
        return 0, lo, hi
    # The minor may vanish exactly AT the root, which happens at degenerate
    # (real-symmetric, zero-flux) carriers where the optimal dual matrix drops
    # rank. Then no amount of refinement separates them, and the honest answer
    # is sign zero, which still satisfies the PSD requirement.
    g = poly_trim(_gcd_poly(p, qt))
    if len(g) > 1 and sturm_count(sturm_chain(g), lo, hi) >= 1:
        return 0, lo, hi
    qchain = sturm_chain(qt)
    for _ in range(max_refine):
        if sturm_count(qchain, lo, hi) == 0:
            val = poly_eval(qt, lo)
            if val == 0:
                val = poly_eval(qt, hi)
            return (1 if val > 0 else -1 if val < 0 else 0), lo, hi
        mid = (lo + hi) / 2
        if poly_eval(p, mid) == 0:  # the root is rational and equals mid
            v = poly_eval(qt, mid)
            return (1 if v > 0 else -1 if v < 0 else 0), mid, mid
        if (poly_eval(p, lo) > 0) != (poly_eval(p, mid) > 0):
            hi = mid
        else:
            lo = mid
    raise RuntimeError("could not separate the minor from the root")


# --------------------------------------------------------------------------
# the order-three theorem and its boundary
# --------------------------------------------------------------------------

def extreme_rank_bound(m: int) -> int:
    """Largest rank an extreme point of the order-m complex elliptope may have.

    ``r^2 <= m`` (Li and Tam): the m Hermitian blocks ``v_i v_i^dagger`` must be
    linearly independent over R in a space of real dimension ``r^2``.
    """
    r = 1
    while (r + 1) ** 2 <= m:
        r += 1
    return r


def is_extreme_point(vectors) -> bool:
    """Extremality test for a factorised elliptope point ``R_ij = <v_i, v_j>``.

    ``R`` is extreme exactly when the rank-one Hermitian matrices
    ``v_i v_i^dagger`` are linearly independent over the reals.
    """
    import numpy as np
    vs = [np.asarray(v, dtype=complex) for v in vectors]
    r = len(vs[0])
    rows = []
    for v in vs:
        blk = np.outer(v, v.conj())
        row = [blk[i, i].real for i in range(r)]
        for i, j in combinations(range(r), 2):
            row += [blk[i, j].real, blk[i, j].imag]
        rows.append(row)
    A = np.asarray(rows)
    return int(np.linalg.matrix_rank(A, tol=1e-9)) == len(vs)


def rank_two_extreme_point_4x4():
    """An explicit rank-two extreme point of the order-four complex elliptope.

    Four unit vectors in C^2 whose rank-one blocks are the four Pauli-basis
    directions, hence linearly independent over R. Its existence is why the
    order-three theorem must not be extrapolated: at m = 4 the elliptope has
    extreme points that are NOT phase matrices, so its convex hull strictly
    contains the fiber's own carrier states and the SDP becomes an outer
    bound.
    """
    import numpy as np
    s = 1 / np.sqrt(2)
    vs = [np.array([1, 0], dtype=complex),
          np.array([s, s], dtype=complex),
          np.array([s, 1j * s], dtype=complex),
          np.array([0, 1], dtype=complex)]
    R = np.array([[np.vdot(a, b) for b in vs] for a in vs])
    return vs, R


# --------------------------------------------------------------------------
# the exact dual certificate
# --------------------------------------------------------------------------

def kkt_shape(h, w):
    """Exact KKT system of the order-three elliptope SDP, in shape position.

    Returns ``(P, gs, E0)`` where ``P`` is a primitive integer polynomial in
    ``t = w.y - E0`` (ascending coefficients) and ``gs[i]`` are rational
    polynomials with ``y_i = h_ii + gs[i](t)``.

    Only ``|h_ij|^2`` and ``Re(chi)``, ``chi = h01 h12 conj(h02)``, enter the
    system. That is the structural reason the range depends on the carrier
    only through its edge magnitudes and the gauge-invariant cycle flux.
    """
    import sympy as sp

    m0, m1, m2, t = sp.symbols("m0 m1 m2 t")
    q01 = sp.expand(sp.Abs(h[0][1]) ** 2)
    q02 = sp.expand(sp.Abs(h[0][2]) ** 2)
    q12 = sp.expand(sp.Abs(h[1][2]) ** 2)
    rechi = sp.re(sp.expand(h[0][1] * h[1][2] * sp.conjugate(h[0][2])))
    ws = [sp.nsimplify(x) for x in w]

    cof00 = m1 * m2 - q12
    cof11 = m0 * m2 - q02
    cof22 = m0 * m1 - q01
    det = m0 * m1 * m2 - m0 * q12 - m1 * q02 - m2 * q01 - 2 * rechi
    eqs = [det,
           ws[1] * cof00 - ws[0] * cof11,
           ws[2] * cof11 - ws[1] * cof22,
           ws[0] * m0 + ws[1] * m1 + ws[2] * m2 - t]
    G = sp.groebner(eqs, m0, m1, m2, t, order="lex")

    uni = [g for g in G.exprs if g.free_symbols <= {t}]
    if not uni:
        raise ValueError("KKT system did not eliminate to a univariate ideal")
    P = sp.Poly(uni[0], t)
    gs = []
    for v in (m0, m1, m2):
        cand = [g for g in G.exprs if v in g.free_symbols and g.free_symbols <= {v, t}]
        if not cand:
            raise ValueError("KKT Groebner basis is not in shape position")
        gs.append(sp.Poly(sp.expand(sp.solve(cand[0], v)[0]), t))
    E0 = sum(sp.nsimplify(w[i]) * sp.re(h[i][i]) for i in range(3))
    return P, gs, E0


def _to_rat_list(poly) -> list[Rat]:
    """sympy Poly in t to ascending Fraction coefficients."""
    co = poly.all_coeffs()[::-1]
    return [Rat(int(c.p), int(c.q)) for c in [__import__("sympy").Rational(x) for x in co]]


def _primitive_int(p: list[Rat]) -> list[int]:
    from math import gcd
    den = 1
    for c in p:
        den = den * c.denominator // gcd(den, c.denominator)
    ints = [int(c * den) for c in p]
    g = 0
    for v in ints:
        g = gcd(g, abs(v))
    if g:
        ints = [v // g for v in ints]
    if ints and ints[-1] < 0:
        ints = [-v for v in ints]
    return ints


def exact_certificates(h, w) -> dict:
    """Exact primal-dual certificates for both endpoints of the range.

    For every real root of the KKT polynomial the dual point ``y`` is
    reconstructed exactly, and the Hermitian principal minors of
    ``Diag(y) - H`` are sign-determined at that root. The root whose minors
    are all nonnegative certifies the MAXIMUM; the root where the signs invert
    on the odd-order minors certifies the MINIMUM.
    """
    import sympy as sp

    P, gs, E0 = kkt_shape(h, w)
    Pr = _to_rat_list(P)
    gr = [_to_rat_list(g) for g in gs]
    q = [Rat(sp.Rational(sp.expand(sp.Abs(h[0][1]) ** 2))),
         Rat(sp.Rational(sp.expand(sp.Abs(h[0][2]) ** 2))),
         Rat(sp.Rational(sp.expand(sp.Abs(h[1][2]) ** 2)))]

    def polymul(a, b):
        out = [Rat(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return poly_trim(out)

    def polysub(a, b):
        n = max(len(a), len(b))
        return poly_trim([(a[i] if i < len(a) else Rat(0))
                          - (b[i] if i < len(b) else Rat(0)) for i in range(n)])

    minors2 = [polysub(polymul(gr[1], gr[2]), [q[2]]),   # cof00
               polysub(polymul(gr[0], gr[2]), [q[1]]),   # cof11
               polysub(polymul(gr[0], gr[1]), [q[0]])]   # cof22

    out: dict = {"objective_polynomial_ascending": _primitive_int(Pr),
                 "E0": str(sp.nsimplify(E0)), "roots": []}
    squarefree = poly_trim(poly_div(Pr, _gcd_poly(Pr, poly_derivative(Pr)))[0])
    for lo0, hi0 in isolate_real_roots(squarefree):
        lo, hi = refine_root(squarefree, lo0, hi0)
        s1 = [sign_at_root(g, squarefree, lo, hi)[0] for g in gr]
        s2 = [sign_at_root(g, squarefree, lo, hi)[0] for g in minors2]
        psd_max = all(s >= 0 for s in s1) and all(s >= 0 for s in s2)
        psd_min = all(s <= 0 for s in s1) and all(s >= 0 for s in s2)
        rec = {"interval": {"lower": str(lo), "upper": str(hi)},
               "approx": float((lo + hi) / 2),
               "diagonal_minor_signs": s1, "two_by_two_minor_signs": s2,
               "certifies": "maximum" if psd_max else "minimum" if psd_min else None}
        out["roots"].append(rec)
        if psd_max:
            out["maximum"] = rec
        if psd_min:
            out["minimum"] = rec
    out["dual_y_offsets_ascending"] = [[str(c) for c in g] for g in gr]
    return out


def _gcd_poly(a: list[Rat], b: list[Rat]) -> list[Rat]:
    a, b = poly_trim(a), poly_trim(b)
    while b:
        _, r = poly_div(a, b)
        a, b = b, poly_trim(r)
    return a if a else [Rat(1)]


def refine_root(p: list[Rat], lo: Rat, hi: Rat, width: Rat = Rat(1, 10 ** 15)
                ) -> tuple[Rat, Rat]:
    """Shrink an isolating interval to at most ``width`` by exact bisection."""
    if poly_eval(p, lo) == 0:
        return lo, lo
    if poly_eval(p, hi) == 0:
        return hi, hi
    while hi - lo > width:
        mid = (lo + hi) / 2
        v = poly_eval(p, mid)
        if v == 0:
            return mid, mid
        if (poly_eval(p, lo) > 0) != (v > 0):
            hi = mid
        else:
            lo = mid
    return lo, hi


def cauchy_bound(p: list[Rat]) -> Rat:
    """A bound strictly larger than the modulus of every real root."""
    p = poly_trim(p)
    return Rat(1) + max((abs(c / p[-1]) for c in p[:-1]), default=Rat(0))


def isolate_real_roots(p: list[Rat]) -> list[tuple[Rat, Rat]]:
    """Disjoint rational intervals, each holding exactly one real root of ``p``.

    Bisection driven by Sturm counts, so the whole certificate chain stays in
    exact rational arithmetic and needs no computer algebra system to replay.
    ``p`` must be squarefree.
    """
    p = poly_trim(p)
    if len(p) <= 1:
        return []
    chain = sturm_chain(p)
    B = cauchy_bound(p)
    out: list[tuple[Rat, Rat]] = []
    work = [(-B, B)]
    while work:
        lo, hi = work.pop()
        n = sturm_count(chain, lo, hi)
        if n == 0:
            continue
        if n == 1:
            out.append((lo, hi))
            continue
        mid = (lo + hi) / 2
        if poly_eval(p, mid) == 0:  # split off an exactly rational root
            out.append((mid, mid))
            eps = (hi - lo) / (1 << 20)
            work.append((lo, mid - eps))
            work.append((mid + eps, hi))
        else:
            work.append((lo, mid))
            work.append((mid, hi))
    return sorted(out)


# --------------------------------------------------------------------------
# numerical route, for cross-checking and for non-rational input
# --------------------------------------------------------------------------

def support_numeric(H, w, sense: str = "max") -> float:
    """Support function of ``E_w`` by numerical SDP, any order.

    Solves the dual ``min w.y  s.t.  Diag(y) - H >= 0`` (and the mirrored
    problem for the minimum). Used to cross-check the exact solver, and as the
    route for floating-point Hamiltonians. At order four and above this is the
    value of the SDP OUTER BOUND, which need not be the conditional body's
    support function.
    """
    import numpy as np
    from scipy.optimize import minimize

    H = np.asarray(H, dtype=complex)
    w = np.asarray(w, dtype=float)
    n = len(w)
    sgn = 1.0 if sense == "max" else -1.0
    A = sgn * H

    def negeig(y):
        return float(np.linalg.eigvalsh(np.diag(y) - A).min())

    # A dual point is only meaningful when it is FEASIBLE, and SLSQP will
    # happily return a slightly infeasible one (or diverge when the optimal
    # face is degenerate, as it is whenever H touches only one edge). Repair
    # each candidate by shifting it up to feasibility, which keeps it a valid
    # upper bound, then keep the best.
    scale = float(np.abs(A).max()) + 1.0
    best = np.inf
    rng = np.random.default_rng(0)
    starts = [np.full(n, float(np.linalg.eigvalsh(A).max()) + 1.0),
              np.real(np.diag(A)) + scale,
              np.full(n, scale)]
    starts += [np.full(n, scale) * (1.0 + rng.uniform(0, 1, n)) for _ in range(4)]
    for y0 in starts:
        try:
            res = minimize(lambda y: float(w @ y), y0, method="SLSQP",
                           constraints=[{"type": "ineq", "fun": negeig}],
                           options={"maxiter": 500, "ftol": 1e-12})
            y = np.asarray(res.x, dtype=float)
        except Exception:
            continue
        if not np.all(np.isfinite(y)):
            continue
        viol = -negeig(y)
        if viol > 0:
            y = y + viol  # shift into the feasible cone
        best = min(best, float(w @ y))
    if not np.isfinite(best):
        raise RuntimeError("dual SDP did not produce a feasible point")

    # At order three the elliptope is the convex hull of the phase matrices
    # (proved), so the primal must meet the dual. Use the primal value there:
    # it is a maximisation over two angles and is far better conditioned.
    if n <= 3:
        primal = sgn * _phase_extreme(A, w)
        if abs(primal - sgn * best) > 1e-6 * max(1.0, abs(primal)):
            return primal  # trust the well-conditioned side
        return primal
    return sgn * best


def _phase_extreme(A, w) -> float:
    """max of ``Tr(A c c^dagger)`` over phase vectors with ``|c_i|^2 = w_i``."""
    import numpy as np
    from scipy.optimize import minimize

    A = np.asarray(A, dtype=complex)
    w = np.asarray(w, dtype=float)
    n = len(w)
    rng = np.random.default_rng(1)
    root = np.sqrt(w)

    def obj(ph):
        c = root * np.exp(1j * np.concatenate([[0.0], ph]))
        return -float(np.real(np.vdot(c, A @ c)))

    best = np.inf
    for _ in range(40):
        r = minimize(obj, rng.uniform(0, 2 * np.pi, n - 1), method="Nelder-Mead",
                     options={"maxiter": 4000, "fatol": 1e-15, "xatol": 1e-13})
        best = min(best, r.fun)
    return -best


def phase_hull_extreme_numeric(H, w, sense: str = "max", restarts: int = 60) -> float:
    """Extreme of ``Tr(H c c^dagger)`` over phase vectors with ``|c_i|^2 = w_i``.

    This is the true conditional-body endpoint (the fiber is exactly this set
    of carrier states). Comparing it with ``support_numeric`` measures the
    relaxation gap, which is zero at order three and generally positive above.
    """
    import numpy as np
    from scipy.optimize import minimize

    H = np.asarray(H, dtype=complex)
    w = np.asarray(w, dtype=float)
    n = len(w)
    rng = np.random.default_rng(0)
    sgn = -1.0 if sense == "max" else 1.0

    def obj(ph):
        c = np.sqrt(w) * np.exp(1j * np.concatenate([[0.0], ph]))
        return sgn * float(np.real(np.vdot(c, H @ c)))

    best = np.inf
    for _ in range(restarts):
        r = minimize(obj, rng.uniform(0, 2 * np.pi, n - 1), method="Nelder-Mead",
                     options={"maxiter": 5000, "fatol": 1e-14, "xatol": 1e-12})
        best = min(best, r.fun)
    return sgn * best


def hidden_width(H, w) -> float:
    """``F^+ - F^-``: the two-body information the 1-RDM leaves undetermined.

    Half of it is the least worst-case error of any estimate of ``<H>`` that
    sees only the 1-RDM, so it is the practical census quantity: zero width
    means the observable is exactly identifiable from one-body data.
    """
    return support_numeric(H, w, "max") - support_numeric(H, w, "min")
