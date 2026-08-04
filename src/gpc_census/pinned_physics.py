"""Physics on a pinned Borland-Dennis fiber: flux, currents, and real Hamiltonians.

``elliptope.py`` supplies the convex geometry (the conditional 2-RDM body and
its support function). This module supplies what that geometry is applied to:

- the unit-edge flux family and its pair-current operator, so the support
  oracle produces an exact energy/current joint body rather than one more
  energy curve;
- second-quantized machinery (occupation-string operator algebra, full CI in
  the three-fermion sector, natural orbitals) and the Slater-Condon projection
  of a Hamiltonian onto the three-determinant carrier, so the intervals attach
  to a recognizable Hamiltonian rather than to an abstract 3x3 matrix;
- the quasipinning enlargement bound, without which the exact solver applies
  only to exactly pinned input and is therefore useless on real data.

CONVENTIONS. Spin orbitals are indexed 0..d-1 and determinants are sorted
tuples of occupied indices. Two-electron integrals use the chemist ordering

    V = (1/2) sum_{pqrs} g[p,q,r,s] a_p^dag a_r^dag a_s a_q,

which is what an FCIDUMP file stores as ``g pq rs``. The Borland-Dennis pinned
carrier is ``DETS`` below, the standard three determinants selected by
saturation of ``lambda_4 <= lambda_5 + lambda_6``.
"""
from __future__ import annotations

from itertools import combinations

DETS = ((0, 1, 2), (0, 3, 4), (1, 3, 5))
"""The Borland-Dennis pinned carrier, in the natural-orbital basis."""


# --------------------------------------------------------------------------
# occupation-string operator algebra
# --------------------------------------------------------------------------

def apply_ops(det: tuple[int, ...], ops) -> tuple[tuple[int, ...], int] | None:
    """Apply a product of creation/annihilation operators to a determinant.

    ``ops`` is a sequence of ``(kind, index)`` with ``kind`` in ``{"c", "a"}``,
    applied RIGHT TO LEFT as written. Returns the resulting determinant and the
    fermionic sign, or None when the result is zero. Signs come from counting
    occupied orbitals to the left of the acted index, which is the Jordan-Wigner
    string, so no sign convention has to be guessed anywhere else.
    """
    occ = list(det)
    sign = 1
    for kind, idx in reversed(list(ops)):
        if kind == "a":
            if idx not in occ:
                return None
            pos = occ.index(idx)
            sign *= (-1) ** pos
            occ.pop(pos)
        else:
            if idx in occ:
                return None
            pos = sum(1 for o in occ if o < idx)
            sign *= (-1) ** pos
            occ.insert(pos, idx)
    return tuple(occ), sign


def excitation_sign(src: tuple[int, ...], dst: tuple[int, ...]):
    """For a double excitation ``src -> dst``, the removed/added pairs and sign.

    Returns ``(removed, added, sign)`` with ``sign`` defined by
    ``a_r^dag a_s^dag a_q a_p |src> = sign |dst>`` for
    ``removed = (p, q)`` and ``added = (r, s)``, or None if the determinants do
    not differ by exactly two orbitals.
    """
    rem = tuple(sorted(set(src) - set(dst)))
    add = tuple(sorted(set(dst) - set(src)))
    if len(rem) != 2 or len(add) != 2:
        return None
    p, q = rem
    r, s = add
    got = apply_ops(src, [("c", r), ("c", s), ("a", q), ("a", p)])
    if got is None:
        return None
    res, sign = got
    return rem, add, (sign if res == dst else -sign)


# --------------------------------------------------------------------------
# full CI in the N-fermion sector
# --------------------------------------------------------------------------

def determinants(d: int, n: int) -> list[tuple[int, ...]]:
    return list(combinations(range(d), n))


def ci_hamiltonian(h1, g2, d: int, n: int):
    """Full CI matrix of ``sum h1 a^dag a + (1/2) sum g2 a^dag a^dag a a``."""
    import numpy as np

    dets = determinants(d, n)
    index = {D: i for i, D in enumerate(dets)}
    H = np.zeros((len(dets), len(dets)), dtype=complex)
    for i, D in enumerate(dets):
        for p in range(d):
            for q in range(d):
                if h1[p][q] == 0:
                    continue
                got = apply_ops(D, [("c", p), ("a", q)])
                if got is None:
                    continue
                E, sg = got
                H[index[E], i] += sg * h1[p][q]
        for p in range(d):
            for q in range(d):
                for r in range(d):
                    for s in range(d):
                        v = g2[p][q][r][s]
                        if v == 0:
                            continue
                        got = apply_ops(D, [("c", p), ("c", r), ("a", s), ("a", q)])
                        if got is None:
                            continue
                        E, sg = got
                        H[index[E], i] += 0.5 * sg * v
    return H, dets


def one_rdm_ci(vec, dets, d: int):
    """1-RDM of a CI vector: ``gamma[p][q] = <Psi| a_p^dag a_q |Psi>``."""
    import numpy as np

    index = {D: i for i, D in enumerate(dets)}
    g = np.zeros((d, d), dtype=complex)
    for i, D in enumerate(dets):
        if vec[i] == 0:
            continue
        for p in range(d):
            for q in range(d):
                got = apply_ops(D, [("c", p), ("a", q)])
                if got is None:
                    continue
                E, sg = got
                g[p, q] += np.conj(vec[index[E]]) * sg * vec[i]
    return g


def natural_orbitals(gamma):
    """Descending occupations and the unitary to the natural-orbital basis."""
    import numpy as np

    vals, vecs = np.linalg.eigh(gamma)
    order = np.argsort(vals)[::-1]
    return np.real(vals[order]), vecs[:, order]


def rotate_ci(vec, dets, U):
    """Transport a CI vector through a one-body basis change.

    Acts by the third exterior power of ``U``: each determinant maps to the
    3x3 minors of ``U`` on the corresponding index sets.
    """
    import numpy as np

    index = {D: i for i, D in enumerate(dets)}
    out = np.zeros(len(dets), dtype=complex)
    for j, D in enumerate(dets):
        if vec[j] == 0:
            continue
        for E in dets:
            minor = U[np.ix_(list(E), list(D))]
            det = np.linalg.det(minor)
            if det != 0:
                out[index[E]] += det * vec[j]
    return out


# --------------------------------------------------------------------------
# Borland-Dennis pinning
# --------------------------------------------------------------------------

def borland_dennis_slack(lams) -> float:
    """``lambda_5 + lambda_6 - lambda_4`` (1-indexed), zero exactly when pinned.

    The (3,6) system's three equalities hold identically for any 1-RDM of a
    three-fermion state; the single inequality is what a state can saturate,
    and its slack is the distance to the pinned face.
    """
    return float(lams[4] + lams[5] - lams[3])


def carrier_weight(vec_no, dets) -> tuple[float, list[float]]:
    """Weight on the pinned carrier and the normalised carrier weights.

    Returns ``(off_carrier_norm_squared, w)`` where ``w`` sums to one.
    """
    index = {D: i for i, D in enumerate(dets)}
    amps = [vec_no[index[D]] for D in DETS]
    on = sum(abs(a) ** 2 for a in amps)
    eps2 = max(0.0, 1.0 - float(on))
    w = [float(abs(a) ** 2 / on) for a in amps] if on > 0 else [0.0, 0.0, 0.0]
    return eps2, w


def carrier_projection(h1, g2, dets_carrier=DETS):
    """Slater-Condon projection of a Hamiltonian onto the carrier.

    Diagonal entries are the ordinary Slater-Condon determinant energies; every
    off-diagonal entry of this carrier is a DOUBLE excitation, so it is a single
    signed antisymmetrized two-electron integral. Returns the Hermitian 3x3
    carrier matrix.
    """
    import numpy as np

    d = 1 + max(max(D) for D in dets_carrier)
    m = len(dets_carrier)
    H = np.zeros((m, m), dtype=complex)
    for i, D in enumerate(dets_carrier):
        acc = sum(h1[p][p] for p in D)
        for p in D:
            for q in D:
                acc += 0.5 * (g2[p][p][q][q] - g2[p][q][q][p])
        H[i, i] = acc
    for i, D in enumerate(dets_carrier):
        for j, E in enumerate(dets_carrier):
            if i >= j:
                continue
            info = excitation_sign(D, E)
            if info is None:
                raise ValueError("carrier determinants must differ by a double "
                                 "excitation")
            (p, q), (r, s), sg = info
            val = g2[r][p][s][q] - g2[r][q][s][p]
            H[j, i] = sg * val
            H[i, j] = np.conj(H[j, i])
    del d
    return H


# --------------------------------------------------------------------------
# the unit-edge flux family and its pair current
# --------------------------------------------------------------------------

def flux_carrier(phi):
    """``[[0,1,1],[1,0,e^{i phi}],[1,e^{-i phi},0]]``: unit edges, cycle flux phi."""
    import numpy as np
    e = np.exp(1j * phi)
    return np.array([[0, 1, 1], [1, 0, e], [1, np.conj(e), 0]], dtype=complex)


def current_carrier(phi):
    """``J = -dH/dphi`` for the unit-edge family: the pair-current operator."""
    import numpy as np
    e = np.exp(1j * phi)
    J = np.zeros((3, 3), dtype=complex)
    J[1, 2] = -1j * e
    J[2, 1] = 1j * np.conj(e)
    return J


def current_range_exact(w):
    """Exact range of the pair current: ``[-2 sqrt(w1 w2), 2 sqrt(w1 w2)]``.

    The current couples only to the fluxed edge, so its expectation is
    ``2 Re(conj(c1) c2 * phase)`` whose modulus is maximised at
    ``2 sqrt(w1 w2)`` by a free relative phase. Independent of the flux, which
    is why the hidden current is not a small perturbative effect.
    """
    import sympy as sp
    half = sp.sqrt(4 * sp.nsimplify(w[1]) * sp.nsimplify(w[2]))
    return -half, half


def joint_body_support(H, J, w, n_directions: int = 180):
    """Support function of the convex energy/current joint body.

    ``(Tr(H rho), Tr(J rho))`` is linear in rho, so the image of the elliptope
    is a convex planar body, and sweeping the direction recovers it. At carrier
    order three the elliptope is the conditional body, so this is exact.
    """
    import numpy as np

    from . import elliptope as el

    H = np.asarray(H, dtype=complex)
    J = np.asarray(J, dtype=complex)
    out = []
    for k in range(n_directions):
        th = 2 * np.pi * k / n_directions
        A = np.cos(th) * H + np.sin(th) * J
        out.append((float(th), float(el.support_numeric(A, w, "max"))))
    return out


def envelope_current(phi, w, sense: str = "max", delta: float = 1e-6) -> float:
    """``-dF^{+/-}/dphi``, the current carried by the extremizing state.

    By the envelope theorem the derivative of the optimal value equals the
    derivative of the objective at the optimizer, so this is exactly
    ``Tr(J rho)`` there wherever the optimizer is unique.
    """
    from . import elliptope as el
    hi = el.support_numeric(flux_carrier(phi + delta), w, sense)
    lo = el.support_numeric(flux_carrier(phi - delta), w, sense)
    return -float((hi - lo) / (2 * delta))


# --------------------------------------------------------------------------
# quasipinning enlargement
# --------------------------------------------------------------------------

def dual_optimal_y(H, w, sense: str = "max"):
    """The optimal dual point, whose sup-norm is the Lipschitz constant in w.

    ``F^+(w) = min { w.y : Diag(y) - H >= 0 }`` is a minimum of LINEAR
    functions of w over a feasible set that does not depend on w, so
    ``|F^+(w) - F^+(w')| <= ||w - w'||_1 * ||y||_inf``. The dual certificate is
    therefore also the sensitivity of the interval to the fiber weights.
    """
    import numpy as np
    from scipy.optimize import minimize

    H = np.asarray(H, dtype=complex)
    w = np.asarray(w, dtype=float)
    sgn = 1.0 if sense == "max" else -1.0
    A = sgn * H
    y0 = np.full(len(w), float(np.linalg.eigvalsh(A).max()) + 1.0)
    res = minimize(lambda y: float(w @ y), y0, method="SLSQP",
                   constraints=[{"type": "ineq",
                                 "fun": lambda y: float(
                                     np.linalg.eigvalsh(np.diag(y) - A).min())}],
                   options={"maxiter": 800, "ftol": 1e-14})
    return sgn * res.x


def enlargement_bound(H_norm: float, eps: float, dw_l1: float,
                      y_inf: float) -> float:
    """Half-width of the quasipinning enlargement of a pinned interval.

    Two independent effects, and conflating them understates the error:

    1. The state is not exactly on the carrier. With
       ``eps^2 = ||(I - P)Psi||^2`` and ``Phi = P Psi / ||P Psi||``,
       expanding ``<Psi|H|Psi>`` around ``<Phi|H|Phi>`` gives an error at most
       ``2 ||H|| (eps + eps^2)``.
    2. The projected state's own fiber weights differ from the pinned ones,
       moving the interval by at most ``||w - w0||_1 ||y||_inf`` because the
       dual optimum is the exact Lipschitz constant.

    Any quasipinning theorem supplying ``eps^2 <= C_gap D(lambda)`` plugs into
    the first term; none is assumed here, so the bound is stated in terms of
    the measured ``eps`` and is self-contained.
    """
    return 2.0 * H_norm * (eps + eps * eps) + dw_l1 * y_inf


# --------------------------------------------------------------------------
# a recognizable Hamiltonian family
# --------------------------------------------------------------------------

def tv_ring_integrals(n_orb: int, t: float, V: float, twist: float = 0.0,
                      gradient: float = 0.0):
    """Spinless fermions on a ring: hopping ``t``, nearest-neighbour ``V``.

    Preferred over the Hubbard model for this benchmark because a spinful
    three-electron ground state is at least a spin doublet, so its 1-RDM is not
    well defined. Spinless fermions avoid that, and the on-site ``gradient``
    breaks the ring's reflection symmetry, which otherwise forces PAIRS OF
    EQUAL natural occupations and leaves the pinned carrier ambiguous. The
    twist puts a Peierls phase on the wrap-around bond, the standard way a real
    model acquires genuinely complex hopping.
    """
    import numpy as np

    h1 = np.zeros((n_orb, n_orb), dtype=complex)
    g2 = np.zeros((n_orb,) * 4, dtype=complex)
    for i in range(n_orb):
        h1[i, i] += gradient * i
        j = (i + 1) % n_orb
        phase = np.exp(1j * twist) if j == 0 else 1.0
        h1[i, j] += -t * phase
        h1[j, i] += -t * np.conj(phase)
        g2[i, i, j, j] += V
        g2[j, j, i, i] += V
    return h1, g2


def hubbard_integrals(n_sites: int, t: float, U: float, twist: float = 0.0):
    """Spin-orbital integrals of a Hubbard ring with a twisted boundary.

    Spin orbital ``2*i + sigma`` is site ``i``, spin ``sigma``. CAUTION: with
    an odd electron number the ground state is at least a spin doublet, so its
    1-RDM is basis-dependent and the pinning diagnostics are meaningless
    without first fixing a sector. ``tv_ring_integrals`` is what the benchmark
    uses; this is kept because it is the recognizable model and is useful for
    exercising the integral pipeline.
    """
    import numpy as np

    d = 2 * n_sites
    h1 = np.zeros((d, d), dtype=complex)
    for i in range(n_sites):
        j = (i + 1) % n_sites
        phase = np.exp(1j * twist) if j == 0 and n_sites > 2 else 1.0
        for s in (0, 1):
            p, q = 2 * i + s, 2 * j + s
            h1[p, q] += -t * phase
            h1[q, p] += -t * np.conj(phase)
    g2 = np.zeros((d, d, d, d), dtype=complex)
    for i in range(n_sites):
        pu, pd = 2 * i, 2 * i + 1
        g2[pu, pu, pd, pd] = U
        g2[pd, pd, pu, pu] = U
    return h1, g2


def write_fcidump(path, h1, g2, n_elec: int, ms2: int = 1) -> None:
    """Write integrals in FCIDUMP format (spin-orbital, complex tolerated)."""
    import numpy as np

    h1 = np.asarray(h1)
    g2 = np.asarray(g2)
    d = h1.shape[0]
    lines = [f" &FCI NORB={d},NELEC={n_elec},MS2={ms2},", " ORBSYM=" + "1," * d,
             " ISYM=1,", " &END"]
    def fmt(v, p, q, r, s):
        # plain Python floats: numpy scalars repr as "np.float64(...)", which
        # no FCIDUMP reader will parse
        return (f"{float(v.real):.17g} {float(v.imag):.17g} "
                f"{p + 1} {q + 1} {r + 1} {s + 1}")

    for p in range(d):
        for q in range(d):
            for r in range(d):
                for s in range(d):
                    v = complex(g2[p, q, r, s])
                    if v != 0:
                        lines.append(fmt(v, p, q, r, s))
    for p in range(d):
        for q in range(d):
            v = complex(h1[p, q])
            if v != 0:
                lines.append(fmt(v, p, q, -1, -1))
    lines.append("0 0 0 0 0 0")
    open(path, "w").write("\n".join(lines) + "\n")


def read_fcidump(path):
    """Read the FCIDUMP written by ``write_fcidump``. Returns ``(h1, g2, norb)``."""
    import numpy as np

    text = open(path).read().splitlines()
    head = [ln for ln in text if "NORB" in ln.upper()]
    norb = int(head[0].upper().split("NORB=")[1].split(",")[0])
    h1 = np.zeros((norb, norb), dtype=complex)
    g2 = np.zeros((norb,) * 4, dtype=complex)
    started = False
    for ln in text:
        if not started:
            if ln.strip().upper().startswith("&END") or ln.strip() == "/":
                started = True
            continue
        parts = ln.split()
        if len(parts) != 6:
            continue
        re_, im_ = float(parts[0]), float(parts[1])
        p, q, r, s = (int(x) for x in parts[2:])
        if p == q == r == s == 0:
            continue
        if r == 0 and s == 0:
            h1[p - 1, q - 1] = complex(re_, im_)
        else:
            g2[p - 1, q - 1, r - 1, s - 1] = complex(re_, im_)
    return h1, g2, norb
