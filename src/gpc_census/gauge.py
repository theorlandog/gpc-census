"""Minimization over the NATURAL-ORBITAL GAUGE of a certified state.

Natural orbitals are unique only up to unitary rotation WITHIN degenerate
occupation eigenspaces. Every census vertex has a degenerate spectrum, so every
certified state carries a gauge freedom

    G(lambda) = U(m_1) x ... x U(m_k)

with m_k the eigenvalue multiplicities, acting on the state through the N-th
compound (exterior power) of the block-diagonal unitary. Every point of that
orbit has 1-RDM exactly diag(lambda): the rotation is block-diagonal in the
eigenspaces, so U diag(lambda) U^dag = diag(lambda) identically, with no
solver, no drift, and no gate needed. Every point is therefore an equally
legitimate natural-orbital representative of the SAME vertex, and the support
recorded in the census library is the support of one arbitrary point of this
family.

This module minimizes the support over that family. The quantity it bounds is

    s_Q^NO(lambda) = min support over states with rho exactly diag(lambda),

the census quantity, bounded below by s_I. It is NOT the orbit-minimal support
over the full U(d) orbit, which is a basis-free property of an individual state
and generally smaller; that one is bounded below by s_M and is computed by
`orbit_minimal` callers of the same machinery with `blocks=[[0..d-1]]`.

Method. The workhorse is the EXACT CLOSURE solve: for a candidate determinant,
solve directly for the block-diagonal U that annihilates its amplitude, by a
re-anchored Gauss-Newton whose Jacobian is the exact derivative of the compound
action. Success lands on the orbit itself rather than on a truncation of it, so
the resulting 1-RDM is diag(lambda) identically. `greedy_gauge_drop` applies
this determinant by determinant; `descend` is a reweighted-L1 Riemannian
descent used to propose candidates when the greedy order is unhelpful.

Finite rotations, not infinitesimal ones, are the point: at v103 all 42
off-diagonal U(4) x U(6) generators move the certified state off its support,
so a first-order method sees nothing at all (docs/prereg_gauge_minimization.md).
A search of this kind is only worth its null results if it can find a reduction
that provably exists, which is what the power test in tests/test_gauge.py
pins.

The profile lemma gives minimality certificates for free. A block-diagonal
unitary preserves each determinant's block profile (|T ∩ B_1|, ..., |T ∩ B_k|),
its compound preserves the norm carried by each profile class, so

    support of ANY gauge representative >= number of occupied profile classes.

A state whose support determinants lie in pairwise distinct profile classes is
exactly gauge-minimal with no computation at all.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.optimize import least_squares

# an amplitude below this (in a unit-norm state) is a candidate zero
ZERO_TOL = 1e-7
# the exact-closure solve must beat this to count as a gauge representative
CLOSURE_TOL = 1e-12
# gate on the reconstructed 1-RDM: max |rho - diag(lambda)|
DIAGONAL_TOL = 1e-10


def degenerate_blocks(spectrum, tol=1e-9):
    """Orbital index blocks of equal occupation, in the shipped orbital order.

    The census ships lambda descending and the certified state in the basis
    where rho = diag(lambda), so block k is a set of orbital indices.
    """
    vals = [float(v) for v in spectrum]
    blocks, cur = [], [0]
    for i in range(1, len(vals)):
        if abs(vals[i] - vals[cur[0]]) > tol:
            blocks.append(cur)
            cur = [i]
        else:
            cur.append(i)
    blocks.append(cur)
    return blocks


def block_profile(det, blocks):
    """(|T ∩ B_1|, ..., |T ∩ B_k|): the gauge-invariant class of a determinant."""
    sets = [set(b) for b in blocks]
    return tuple(sum(1 for o in det if o in s) for s in sets)


def profile_classes(support_dets, blocks):
    """Map profile -> list of support positions sharing it."""
    out: dict[tuple, list[int]] = {}
    for i, t in enumerate(support_dets):
        out.setdefault(block_profile(t, blocks), []).append(i)
    return out


def profile_lower_bound(support_dets, blocks):
    """Number of OCCUPIED profile classes: a gauge-invariant lower bound.

    Every gauge representative of this state has support at least this large,
    because the compound of a block-diagonal unitary preserves the norm carried
    by each profile class and a class carrying nonzero norm needs at least one
    support determinant. Equality certifies gauge-minimality exactly.
    """
    return len(profile_classes(support_dets, blocks))


def reachable_dets(support_dets, blocks, d, n):
    """Every determinant the gauge orbit can reach: same profile as some source."""
    live = {block_profile(t, blocks) for t in support_dets}
    return [t for t in combinations(range(d), n) if block_profile(t, blocks) in live]


def compound_matrix(u, target_dets, source_dets):
    """N-th compound: M[T, S] = det(u[T, S]), the action on determinant amplitudes.

    Gamma(u) psi has amplitudes M @ c. Gamma is a group homomorphism into the
    unitaries, so on the full determinant basis M is unitary.
    """
    u = np.asarray(u, dtype=complex)
    rows = np.array(target_dets, dtype=int)
    cols = np.array(source_dets, dtype=int)
    # (|T|, |S|, N, N) stack of submatrices, one batched determinant call
    sub = u[rows[:, None, :, None], cols[None, :, None, :]]
    return np.linalg.det(sub)


def block_diagonal(params, blocks, d):
    """exp(i H) for the block-diagonal Hermitian H described by `params`.

    Parameters are, per block, the strict upper triangle as (real, imag) pairs
    followed by the diagonal. The diagonal generates the U(1)^d orbital-phase
    gauge, which cannot change any |c_T|; it is carried anyway so the exact
    closure solve can use it to fix phases.
    """
    h = np.zeros((d, d), dtype=complex)
    k = 0
    for b in blocks:
        for a in range(len(b)):
            for c in range(a + 1, len(b)):
                h[b[a], b[c]] = params[k] + 1j * params[k + 1]
                h[b[c], b[a]] = params[k] - 1j * params[k + 1]
                k += 2
        for a in range(len(b)):
            h[b[a], b[a]] = params[k]
            k += 1
    w, v = np.linalg.eigh(h)
    return (v * np.exp(1j * w)) @ v.conj().T


def parameter_count(blocks):
    return sum(len(b) * len(b) for b in blocks)


def random_block_unitary(blocks, d, rng):
    """Haar-ish random element of U(m_1) x ... x U(m_k), embedded in U(d)."""
    u = np.eye(d, dtype=complex)
    for b in blocks:
        m = len(b)
        z = (rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))) / np.sqrt(2)
        q, r = np.linalg.qr(z)
        q = q * (np.diag(r) / np.abs(np.diag(r)))
        u[np.ix_(b, b)] = q
    return u


def _one_body_action(dets, d):
    """A[p, q] with (a_p^dag a_q psi)_T = sum_S A[p, q][T, S] c_S.

    Same fermionic convention as gpc_census.cascade.rdm_tensor, so the derived
    representation dGamma(H) = sum_pq H_pq a_p^dag a_q below is the exact
    generator of the compound action verified in tests.
    """
    idx = {t: i for i, t in enumerate(dets)}
    a = np.zeros((d, d, len(dets), len(dets)), dtype=float)
    for ti, t in enumerate(dets):
        for qi, q in enumerate(t):
            rest = [u for u in t if u != q]
            for p in range(d):
                if p in rest:
                    continue
                s = tuple(sorted(rest + [p]))
                si = idx.get(s)
                if si is None:
                    continue
                a[p, q, si, ti] += (-1) ** qi * (-1) ** s.index(p)
    return a


def block_generators(blocks):
    """Hermitian generators of the gauge algebra that can change |c_T|.

    Off-diagonal only: the diagonal generators are the U(1)^d orbital-phase
    gauge and act by phases alone, so they cannot move any modulus and are
    useless to a sparsity objective (they are carried separately by the exact
    closure solve, which does need them).
    """
    gens = []
    for b in blocks:
        for i in range(len(b)):
            for j in range(i + 1, len(b)):
                gens.append((b[i], b[j], "re"))
                gens.append((b[i], b[j], "im"))
    return gens


def _generator_matrix(gen, d):
    p, q, kind = gen
    h = np.zeros((d, d), dtype=complex)
    if kind == "re":
        h[p, q] = h[q, p] = 1.0
    else:
        h[p, q], h[q, p] = 1j, -1j
    return h


class GaugeOrbit:
    """The gauge orbit of one state, with exact gradients of sparsity surrogates.

    The gradient is taken at the current anchor U by the derived representation:
    d/dt Gamma(U exp(i t E)) psi_0 |_{t=0} = Gamma(U) dGamma(iE) psi_0, and
    dGamma(E) psi_0 is a fixed vector per generator because psi_0 is fixed. So a
    full gradient costs one compound evaluation, not one per parameter.
    """

    def __init__(self, amplitudes, support_dets, blocks, d, target=None):
        self.d = d
        self.blocks = [list(b) for b in blocks]
        self.dets = [tuple(int(o) for o in t) for t in support_dets]
        n = len(self.dets[0])
        self.target = (reachable_dets(self.dets, self.blocks, d, n)
                       if target is None else [tuple(int(o) for o in t) for t in target])
        pos = {t: i for i, t in enumerate(self.target)}
        self.c0 = np.zeros(len(self.target), dtype=complex)
        for t, amp in zip(self.dets, np.asarray(amplitudes, dtype=complex)):
            self.c0[pos[t]] = amp
        self.gens = block_generators(self.blocks)
        a = _one_body_action(self.target, d)
        # y_a = dGamma(i E_a) psi_0, one fixed vector per generator
        self.y = np.array([1j * np.einsum("pq,pqij,j->i",
                                          _generator_matrix(g, d), a, self.c0)
                           for g in self.gens])

    def amplitudes(self, u):
        return compound_matrix(u, self.target, self.target) @ self.c0

    def value_grad(self, u, weights, eps):
        """Reweighted smoothed-L1 value and its exact Riemannian gradient at u."""
        m = compound_matrix(u, self.target, self.target)
        c = m @ self.c0
        r = np.sqrt(np.abs(c) ** 2 + eps ** 2)
        val = float(np.sum(weights * r))
        dc = m @ self.y.T  # columns: d c / d theta_a
        g = (weights / r) @ np.real(np.conj(c)[:, None] * dc)
        return val, np.asarray(g, dtype=float), c

    def value(self, u, weights, eps):
        c = self.amplitudes(u)
        return float(np.sum(weights * np.sqrt(np.abs(c) ** 2 + eps ** 2)))

    def step_unitary(self, u, direction, t):
        h = np.zeros((self.d, self.d), dtype=complex)
        for coeff, gen in zip(direction, self.gens):
            if coeff:
                h += coeff * _generator_matrix(gen, self.d)
        w, v = np.linalg.eigh(h)
        return u @ ((v * np.exp(1j * t * w)) @ v.conj().T)


def descend(orbit, u0=None, reweight_rounds=8, iters=60, eps0=3e-2, eps_min=1e-7,
            smooth=1e-10):
    """Riemannian reweighted-L1 descent over the gauge group. Returns (u, c).

    Descent along geodesics U exp(i t H) with a line search, re-anchored every
    step so the gradient is exact rather than a small-angle approximation. This
    is what makes LARGE rotations reachable: the line search scans t over a full
    period, so a minimum a finite rotation away is found in one step where a
    first-order method sees nothing.

    The reweighting parameter (annealed from eps0) and the SMOOTHING parameter
    are separate on purpose. Smoothing the modulus as sqrt(|c|^2 + s^2) with s
    comparable to the reweighting epsilon makes growing a zero amplitude free up
    to s, so the surrogate rewards spreading and the descent walks away from the
    sparsest point it was handed. Keeping s at 1e-10 leaves an honest L1.
    """
    d = orbit.d
    u = np.eye(d, dtype=complex) if u0 is None else np.array(u0, dtype=complex)
    eps = eps0
    c = orbit.amplitudes(u)
    steps = (4.0, 2.0, 1.0, 0.5, 0.25, 0.12, 0.06, 0.03, 0.01, 0.003, 0.001)
    for _ in range(reweight_rounds):
        weights = 1.0 / (np.abs(c) + eps)
        prev, gprev, dprev = None, None, None
        for _ in range(iters):
            val, g, c = orbit.value_grad(u, weights, smooth)
            gn = float(np.linalg.norm(g))
            if gn < 1e-13:
                break
            # Polak-Ribiere conjugate gradient, restarted when it loses descent
            if gprev is None:
                direction = -g
            else:
                beta = max(0.0, float(g @ (g - gprev)) / float(gprev @ gprev))
                direction = -g + beta * dprev
                if float(direction @ g) >= 0:
                    direction = -g
            gprev, dprev = g, direction
            direction = direction / float(np.linalg.norm(direction))
            best = (val, 0.0)
            for t in steps:
                v = orbit.value(orbit.step_unitary(u, direction, t), weights,
                                smooth)
                if v < best[0]:
                    best = (v, t)
            if best[1] == 0.0:
                gprev = None  # line search failed: restart as steepest descent
                if prev is not None:
                    break
                prev = val
                continue
            u = orbit.step_unitary(u, direction, best[1])
            c = orbit.amplitudes(u)
            if prev is not None and abs(prev - best[0]) < 1e-14:
                break
            prev = best[0]
        eps = max(eps * 0.3, eps_min)
    c = orbit.amplitudes(u)
    return u, c


def close_exactly(c_source, source_dets, target_dets, blocks, d, keep,
                  u_start=None, tries=3, rng=None):
    """Solve for U in the gauge group annihilating every amplitude outside `keep`.

    This is what turns a sweep endpoint into a genuine gauge representative:
    rather than truncating small amplitudes (which perturbs the 1-RDM at first
    order in the leftover), it finds the block-diagonal U for which the rotated
    state is EXACTLY supported inside `keep`. Because U is block-diagonal the
    resulting 1-RDM is diag(lambda) identically.

    Returns (amplitudes over target_dets, residual) or (None, residual).
    """
    rng = np.random.default_rng(0) if rng is None else rng
    drop = [i for i, t in enumerate(target_dets) if t not in keep]
    if not drop:
        return None, 0.0
    orbit = GaugeOrbit(c_source, source_dets, blocks, d, target=target_dets)
    best = (None, float("inf"))
    for k in range(tries):
        u0 = (np.eye(d, dtype=complex) if u_start is None else np.asarray(u_start))
        if k:
            u0 = u0 @ random_block_unitary(blocks, d, rng)
        c, r = _gauss_newton_close(orbit, drop, u0)
        if r < best[1]:
            best = (c, r)
        if best[1] < CLOSURE_TOL:
            break
    return best


def _gauss_newton_close(orbit, drop, u, iters=80, lam0=1e-6):
    """Damped Gauss-Newton for "rotate these amplitudes to exactly zero".

    Re-anchored at every step, so the Jacobian is the EXACT derivative of the
    compound action (columns Gamma(U) dGamma(iE_a) psi_0) rather than a finite
    difference of it. That is what makes the closure solve cheap enough to run
    as the census-wide feasibility test: one compound evaluation per iteration
    instead of one per parameter.
    """
    lam = lam0
    c = orbit.amplitudes(u)
    best = float(np.max(np.abs(c[drop])))
    for _ in range(iters):
        m = compound_matrix(u, orbit.target, orbit.target)
        c = m @ orbit.c0
        r = c[drop]
        if float(np.max(np.abs(r))) < 1e-15:
            break
        j = (m @ orbit.y.T)[drop, :]
        jr = np.vstack([j.real, j.imag])
        rr = np.concatenate([r.real, r.imag])
        step = None
        for _ in range(12):
            a = jr.T @ jr + lam * np.eye(jr.shape[1])
            try:
                delta = -np.linalg.solve(a, jr.T @ rr)
            except np.linalg.LinAlgError:
                lam *= 10.0
                continue
            cand = orbit.step_unitary(u, delta, 1.0)
            val = float(np.max(np.abs(orbit.amplitudes(cand)[drop])))
            if val < best:
                step, best = cand, val
                lam = max(lam * 0.3, 1e-12)
                break
            lam *= 10.0
        if step is None:
            break
        u = step
    c = orbit.amplitudes(u)
    return c, float(np.max(np.abs(c[drop])))


def flattening_ranks(amplitudes, support_dets, blocks):
    """Per profile class, the rank of the class tensor's coarsest flattening.

    Inside a profile class the gauge acts as a TENSOR PRODUCT of exterior powers,
    one factor per block, so every flattening rank of the class tensor is a gauge
    invariant. The support of a tensor is at least the rank of any flattening, so

        min support over the gauge family >= sum over classes of (best rank),

    which refines the profile-class count (that is the special case of a rank-1
    flattening per occupied class). Returns {profile: rank}, the bound being the
    sum of the values.
    """
    blocks = [list(b) for b in blocks]
    dets = [tuple(int(o) for o in t) for t in support_dets]
    c = np.asarray(amplitudes, dtype=complex)
    out = {}
    for prof, idxs in profile_classes(dets, blocks).items():
        active = [k for k, b in enumerate(blocks) if 0 < prof[k] < len(b)]
        if not active:
            out[prof] = 1  # the class is one-dimensional
            continue
        if len(active) == 1:
            k = active[0]
            a, m = prof[k], len(blocks[k])
            if a in (1, m - 1):
                # Lambda^a(C^m) = C^m up to duality and U(m) is transitive on
                # its rays, so one determinant suffices for this class
                out[prof] = 1
            elif a == 2 or a == m - 2:
                # an antisymmetric form: Youla normal form gives exactly
                # rank/2 nonzero entries, an equality not just a bound
                out[prof] = _pfaffian_support(dets, c, idxs, blocks[k], a == 2)
            else:
                out[prof] = 1  # no cheap invariant for 3 <= a <= m-3
            continue
        # several active factors: the best single-factor-versus-rest flattening
        best = 1
        for k in active:
            rows, cols, ent = {}, {}, []
            for i in idxs:
                t = dets[i]
                key_a = tuple(o for o in t if o in blocks[k])
                key_b = tuple(o for o in t if o not in blocks[k])
                rows.setdefault(key_a, len(rows))
                cols.setdefault(key_b, len(cols))
                ent.append((rows[key_a], cols[key_b], c[i]))
            m = np.zeros((len(rows), len(cols)), dtype=complex)
            for r, s, v in ent:
                m[r, s] = v
            best = max(best, int(np.linalg.matrix_rank(m, tol=1e-9)))
        out[prof] = best
    return out


def _pfaffian_support(dets, c, idxs, block, direct):
    """Minimal support of an antisymmetric form under U(m): half its rank.

    A class whose single active factor sits in Lambda^2 (or, dually,
    Lambda^{m-2}) is an antisymmetric matrix; Youla's normal form makes it
    block-diagonal with rank/2 nonzero 2 x 2 blocks, which is exactly rank/2
    occupied determinants, and the rank is a unitary invariant. Equality, not
    just a bound.
    """
    pos = {o: i for i, o in enumerate(block)}
    m = len(block)
    a = np.zeros((m, m), dtype=complex)
    for i in idxs:
        occ = [o for o in dets[i] if o in pos]
        if direct:
            p, q = pos[occ[0]], pos[occ[1]]
        else:  # Lambda^{m-2}: index by the complementary pair
            missing = [o for o in block if o not in occ]
            p, q = pos[missing[0]], pos[missing[1]]
        a[p, q] += c[i]
        a[q, p] -= c[i]
    return max(1, int(np.linalg.matrix_rank(a, tol=1e-9)) // 2)


def flattening_lower_bound(amplitudes, support_dets, blocks):
    """Sum of per-class flattening ranks: a gauge-invariant support lower bound."""
    return int(sum(flattening_ranks(amplitudes, support_dets, blocks).values()))


def greedy_gauge_drop(amplitudes, support_dets, blocks, d, tries=4, seed=0,
                      max_rounds=None):
    """Greedily ask, for each support determinant: is there a gauge U killing it?

    For a candidate determinant T this solves for the block-diagonal U making
    every amplitude outside S \\ {T} vanish exactly. Success is an EXACT gauge
    representative of smaller support, not a truncation: the answer lands on the
    orbit, so its 1-RDM is diag(lambda) identically. Failure of the solve is not
    a proof that no such U exists, only that this multi-start did not find one.

    Returns (amplitudes, support_dets, rounds).
    """
    rng = np.random.default_rng(seed)
    dets = [tuple(int(o) for o in t) for t in support_dets]
    c = np.asarray(amplitudes, dtype=complex)
    n = len(dets[0])
    rounds = []
    limit = len(dets) if max_rounds is None else max_rounds
    for _ in range(limit):
        target = reachable_dets(dets, blocks, d, n)
        landed = False
        for j in np.argsort(np.abs(c)):
            keep = {t for i, t in enumerate(dets) if i != int(j)}
            got, resid = close_exactly(c, dets, target, blocks, d, keep,
                                       tries=tries, rng=rng)
            if got is None or resid > CLOSURE_TOL:
                continue
            live = [i for i, t in enumerate(target) if abs(got[i]) > ZERO_TOL]
            new_dets = [target[i] for i in live]
            if len(new_dets) >= len(dets):
                continue
            rounds.append({"dropped_det": list(dets[int(j)]),
                           "support_after": len(new_dets),
                           "closure_residual": resid})
            dets, c = new_dets, got[live]
            landed = True
            break
        if not landed:
            break
    return c, dets, rounds


# --- the 2-RDM gate -------------------------------------------------------

def _pair_vectors(c, dets, d):
    """v[(p<q)][R] = amplitude of the pair (p, q) over the (N-2)-remainder R."""
    n = len(dets[0])
    rest = list(combinations(range(d), n - 2))
    rpos = {r: i for i, r in enumerate(rest)}
    pairs = list(combinations(range(d), 2))
    v = np.zeros((len(pairs), len(rest)), dtype=complex)
    dpos = {t: i for i, t in enumerate(dets)}
    for pi, (p, q) in enumerate(pairs):
        for r in rest:
            if p in r or q in r:
                continue
            t = tuple(sorted((p, q) + r))
            j = dpos.get(t)
            if j is None:
                continue
            # sign of the permutation taking (p, q, r...) to sorted order
            seq = [p, q] + list(r)
            perm = sorted(range(len(seq)), key=lambda k: seq[k])
            sgn = 1
            seen = [False] * len(perm)
            for k in range(len(perm)):
                if seen[k]:
                    continue
                length, m = 0, k
                while not seen[m]:
                    seen[m] = True
                    m = perm[m]
                    length += 1
                if length % 2 == 0:
                    sgn = -sgn
            v[pi, rpos[r]] = sgn * c[j]
    return v


def two_rdm_spectrum(amplitudes, support_dets, d=None):
    """Eigenvalues of the 2-RDM, descending: a BASIS-FREE fingerprint of a state.

    The 2-RDM transforms by conjugation under any one-body unitary, so its
    spectrum is invariant along the whole U(d) orbit of a state. Two states with
    different 2-RDM spectra are genuinely different states; two with the same
    spectrum and a continuous path between them are the same state in different
    orbitals. This is the cheap gate that separates a NEW vertex attainer from a
    gauge restatement of the library one.
    """
    dets = [tuple(int(o) for o in t) for t in support_dets]
    c = np.asarray(amplitudes, dtype=complex)
    if d is None:
        d = max(o for t in dets for o in t) + 1
    if len(dets[0]) < 2:
        return np.zeros(0)
    v = _pair_vectors(c, dets, d)
    g = v.conj() @ v.T
    return np.sort(np.linalg.eigvalsh(g))[::-1]


def same_state_up_to_orbital_rotation(a_amps, a_dets, b_amps, b_dets, d=None,
                                      tol=1e-9):
    """2-RDM spectra agree: the two states are one state in two orbital bases.

    Returns (verdict, max deviation). A verdict of True is strong evidence, not
    proof: it says the states are indistinguishable by every basis-free
    two-body invariant of this kind. Proof requires exhibiting the one-body
    unitary (see `solve_orbital_rotation`).
    """
    sa = two_rdm_spectrum(a_amps, a_dets, d)
    sb = two_rdm_spectrum(b_amps, b_dets, d)
    n = max(len(sa), len(sb))
    sa = np.pad(sa, (0, n - len(sa)))
    sb = np.pad(sb, (0, n - len(sb)))
    dev = float(np.max(np.abs(sa - sb))) if n else 0.0
    return bool(dev <= tol), dev


def solve_orbital_rotation(a_amps, a_dets, b_amps, b_dets, d, blocks=None,
                           tries=24, seed=0):
    """Solve for the one-body U with Gamma(U) psi_a = psi_b (up to global phase).

    This is the exact version of the 2-RDM gate: it does not merely say the two
    states share every basis-free invariant tested, it exhibits the rotation.
    Pass `blocks` to restrict to the natural-orbital gauge subgroup; omit for
    the full U(d).
    """
    rng = np.random.default_rng(seed)
    blocks = [list(range(d))] if blocks is None else blocks
    a_dets = [tuple(int(o) for o in t) for t in a_dets]
    b_dets = [tuple(int(o) for o in t) for t in b_dets]
    n = len(a_dets[0])
    target = sorted(set(combinations(range(d), n)))
    tpos = {t: i for i, t in enumerate(target)}
    goal = np.zeros(len(target), dtype=complex)
    for t, amp in zip(b_dets, np.asarray(b_amps, dtype=complex)):
        goal[tpos[t]] = amp
    ca = np.asarray(a_amps, dtype=complex)
    n_par = parameter_count(blocks) + 1  # + global phase

    def resid(params):
        u = block_diagonal(params[:-1], blocks, d)
        v = np.exp(1j * params[-1]) * (compound_matrix(u, target, a_dets) @ ca)
        return np.concatenate([(v - goal).real, (v - goal).imag])

    best = (None, float("inf"))
    for k in range(tries):
        p0 = np.zeros(n_par) if k == 0 else rng.normal(scale=1.0, size=n_par)
        sol = least_squares(resid, p0, method="lm", xtol=3e-16, ftol=3e-16,
                            gtol=3e-16, max_nfev=6000)
        r = float(np.max(np.abs(resid(sol.x))))
        if r < best[1]:
            best = (block_diagonal(sol.x[:-1], blocks, d), r)
        if best[1] < 1e-12:
            break
    return best
