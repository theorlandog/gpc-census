"""First-order fiber geometry at a certified state: TWO DISTINCT FIBERS.

Over a vertex spectrum lambda there are two different sets of states, and
conflating them has caused repeated errors in this project (the levi.py audit,
the P4 scorecard reading of v103, a draft of the RESEARCH_NEXT diff). This
module makes the distinction structural: every entry point names its fiber, and
every dimension it reports is written under a fiber-qualified key.

FIXED-RHO fiber
    {psi : rho_1(psi) = rho} for the state's own 1-RDM rho. Conditions: every
    real coordinate of d(rho) vanishes. This is the fiber levi.py's
    fixed_rho_fiber_dim refers to.

FIXED-SPECTRUM fiber
    {psi : spec(rho_1(psi)) = lambda}, the census object. Conditions: only the
    block-diagonal parts P_i d(rho) P_i vanish, where P_i are the spectral
    projectors of rho. Off-block components rotate eigenvectors and move
    eigenvalues at second order only (see docs/silence_rank_lemma.md), so this
    fiber is at least as large as the fixed-rho one and is strictly larger
    whenever some off-block channel is deformable.

Both entry points report the non-gauge tangent dimension and the gauge
dimension separately. Gauge is the orbital-phase group U(1)^d acting by
c_t -> exp(i sum_{o in t} phi_o) c_t (the global phase is contained in it).
Gauge always preserves the spectrum, so it lies wholly inside the
fixed-spectrum kernel; it need not preserve rho itself, so for the fixed-rho
fiber only the part of the gauge span that lies inside the kernel is counted.

All dimensions here are FIRST ORDER (tangent-space) numbers. They are not
proofs about the fiber germ; a positive tangent dimension is upgraded to an
actual deformation only by continuation (see gpc_census.cascade) or by an
exact constant-rank argument.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

FIXED_RHO = "fixed_rho"
FIXED_SPECTRUM = "fixed_spectrum"


@dataclass(frozen=True)
class FiberTangent:
    """First-order tangent data for ONE named fiber.

    fiber: which fiber, "fixed_rho" or "fixed_spectrum". Never omit it, and
        never flatten these two into one unqualified dimension field.
    tangent_dim: dimension of the kernel modulo gauge (the reportable number).
    gauge_dim: dimension of the gauge span lying inside the kernel.
    kernel_dim: raw kernel dimension, gauge included.
    constraint_rank: rank of the linearized constraint map.
    ambient_dim: 2 * |support|, the real dimension of the amplitude space.
    """

    fiber: str
    tangent_dim: int
    gauge_dim: int
    kernel_dim: int
    constraint_rank: int
    ambient_dim: int

    def as_json(self) -> dict:
        """Fiber-qualified keys, safe to merge into a record.

        Emits e.g. fixed_spectrum_tangent_dim / fixed_spectrum_gauge_dim so a
        consumer cannot read a dimension without reading which fiber it is of.
        """
        out = {f"{self.fiber}_{k}": v for k, v in asdict(self).items()
               if k != "fiber"}
        out["fiber"] = self.fiber
        return out


def _rank(m, tol=1e-8):
    if m.size == 0:
        return 0
    s = np.linalg.svd(np.asarray(m, dtype=float), compute_uv=False)
    if s.size == 0 or s[0] == 0:
        return 0
    return int((s > max(tol * s[0], 1e-12)).sum())


def _nullspace(m, tol=1e-8):
    m = np.asarray(m, dtype=float)
    if m.size == 0:
        return np.eye(m.shape[1])
    u, s, vt = np.linalg.svd(m)
    cut = max(tol * s[0], 1e-12) if s.size and s[0] > 0 else 1e-12
    r = int((s > cut).sum())
    return vt[r:]


def _sesquilinear_rdm(x, y, dets, d):
    """<x| a_p^dag a_q |y> in the canonical determinant ordering."""
    idx = {t: i for i, t in enumerate(dets)}
    r = np.zeros((d, d), dtype=complex)
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
                r[p, q] += np.conj(x[si]) * y[ti] * (-1) ** qi * (-1) ** s.index(p)
    return r


def one_rdm(amplitudes, support_dets, d=None):
    """1-RDM <a_p^dag a_q> of sum_t c_t |t> restricted to the given support."""
    c, dets, d = _normalize(amplitudes, support_dets, d)
    return _sesquilinear_rdm(c, c, dets, d)


def _normalize(amplitudes, support_dets, d):
    dets = [tuple(int(o) for o in t) for t in support_dets]
    c = np.asarray(amplitudes, dtype=complex)
    if len(c) != len(dets):
        raise ValueError("amplitude count does not match support size")
    if d is None:
        d = max(o for t in dets for o in t) + 1
    return c, dets, int(d)


def _differentials(c, dets, d):
    """d(rho) for each of the 2m real deformation directions of the amplitudes.

    Returns a list of length 2m of complex dxd Hermitian matrices; direction k
    is the real part of c_k for k < m and the imaginary part for k >= m.
    """
    m = len(dets)
    out = []
    for k in range(2 * m):
        dc = np.zeros(m, dtype=complex)
        dc[k % m] = 1.0 if k < m else 1j
        out.append(_sesquilinear_rdm(dc, c, dets, d) + _sesquilinear_rdm(c, dc, dets, d))
    return out


def _herm_rows(mats, index_pairs):
    """Real coordinates of the listed Hermitian entries, one row per direction."""
    rows = []
    for dr in mats:
        row = []
        for p, q in index_pairs:
            if p == q:
                row.append(dr[p, p].real)
            else:
                row += [dr[p, q].real, dr[p, q].imag]
        rows.append(row)
    return np.array(rows, dtype=float)


def gauge_generators(amplitudes, support_dets, d=None):
    """Real 2m-vectors spanning the orbital-phase gauge action on the support.

    Generator o acts by c_t -> i c_t if orbital o is occupied in t. The global
    phase is the sum of all d generators divided by the particle number, so it
    is inside this span and must not be counted twice.
    """
    c, dets, d = _normalize(amplitudes, support_dets, d)
    rows = []
    for o in range(d):
        occ = np.array([1.0 if o in t else 0.0 for t in dets])
        v = 1j * occ * c
        rows.append(np.concatenate([v.real, v.imag]))
    return np.array(rows, dtype=float)


def spectral_blocks(rho, eig_tol=1e-8):
    """Spectral projector blocks of a Hermitian rho, as lists of eigen-indices.

    Returns (eigenvalues descending, eigenvectors as columns, block index
    lists). Blocks group eigenvalues equal to within eig_tol.
    """
    w, v = np.linalg.eigh(rho)
    order = np.argsort(w)[::-1]
    w, v = w[order], v[:, order]
    blocks, cur = [], [0]
    for i in range(1, len(w)):
        if abs(w[i] - w[cur[0]]) > eig_tol:
            blocks.append(cur)
            cur = [i]
        else:
            cur.append(i)
    blocks.append(cur)
    return w, v, blocks


def _tangent(constraint_rows, gauge, fiber, ambient, tol=1e-8):
    null = _nullspace(constraint_rows, tol=tol)
    kernel_dim = int(null.shape[0])
    # gauge directions that actually satisfy the constraints
    gauge_in_kernel = 0
    if gauge.size:
        gauge_dim_total = _rank(gauge, tol=tol)
        moved = _rank(constraint_rows @ gauge.T, tol=tol)
        gauge_in_kernel = max(0, gauge_dim_total - moved)
    return FiberTangent(
        fiber=fiber,
        tangent_dim=max(0, kernel_dim - gauge_in_kernel),
        gauge_dim=gauge_in_kernel,
        kernel_dim=kernel_dim,
        constraint_rank=int(constraint_rows.shape[0] and _rank(constraint_rows, tol=tol)),
        ambient_dim=ambient,
    )


def fixed_rho_tangent(amplitudes, support_dets, d=None, tol=1e-8):
    """First-order tangent to the FIXED-RHO fiber {psi : rho_1(psi) = rho}.

    Support-restricted: deformations are of the given amplitudes only, not of
    the ambient Fock space. Reports tangent_dim (modulo the gauge directions
    that preserve rho) and gauge_dim separately; see FiberTangent.

    This is NOT the census fixed-spectrum fiber. A zero here says the state is
    first-order rigid as a preimage of its own 1-RDM, which is compatible with
    a positive-dimensional fixed-spectrum fiber (v103 is exactly that case).
    """
    c, dets, d = _normalize(amplitudes, support_dets, d)
    diffs = _differentials(c, dets, d)
    pairs = [(p, q) for p in range(d) for q in range(p, d)]
    rows = _herm_rows(diffs, pairs)
    norm = np.array([[2 * np.real(np.vdot(c, dc)) for dc in _unit_dirs(len(dets))]]).T
    m = np.hstack([rows, norm]).T
    return _tangent(m, gauge_generators(c, dets, d), FIXED_RHO, 2 * len(dets), tol=tol)


def fixed_spectrum_tangent(amplitudes, support_dets, d=None, tol=1e-8,
                           eig_tol=1e-8, blocks=None):
    """First-order tangent to the FIXED-SPECTRUM fiber {psi : spec rho_1 = lambda}.

    Support-restricted, gauge quotiented. Conditions are the block-diagonal
    parts P_i d(rho) P_i in the eigenbasis of rho, one block per DISTINCT
    eigenvalue, plus the norm; off-block components are unconstrained at first
    order (docs/silence_rank_lemma.md). Pass explicit `blocks` (lists of
    orbital indices, valid only when rho is diagonal in the computational
    basis) to bypass the numerical eigendecomposition.

    This is the census fiber. It is the one the sparsification cascade moves
    along, and the one whose dimension the unified formula predicts.
    """
    c, dets, d = _normalize(amplitudes, support_dets, d)
    rho = _sesquilinear_rdm(c, c, dets, d)
    diffs = _differentials(c, dets, d)
    if blocks is None:
        _, vecs, blk = spectral_blocks(rho, eig_tol=eig_tol)
        # rotate every differential into the eigenbasis of rho
        diffs = [vecs.conj().T @ dr @ vecs for dr in diffs]
        groups = blk
    else:
        groups = [list(b) for b in blocks]
    pairs = [(p, q) for g in groups for i, p in enumerate(g) for q in g[i:]]
    rows = _herm_rows(diffs, pairs)
    norm = np.array([[2 * np.real(np.vdot(c, dc)) for dc in _unit_dirs(len(dets))]]).T
    m = np.hstack([rows, norm]).T
    return _tangent(m, gauge_generators(c, dets, d), FIXED_SPECTRUM,
                    2 * len(dets), tol=tol)


def _unit_dirs(m):
    for k in range(2 * m):
        dc = np.zeros(m, dtype=complex)
        dc[k % m] = 1.0 if k < m else 1j
        yield dc


def tangent_report(amplitudes, support_dets, d=None, **kw):
    """Both fibers at once, as a single JSON-ready dict with qualified keys."""
    a = fixed_rho_tangent(amplitudes, support_dets, d, tol=kw.get("tol", 1e-8))
    b = fixed_spectrum_tangent(amplitudes, support_dets, d, **kw)
    out = {k: v for k, v in a.as_json().items() if k != "fiber"}
    out.update({k: v for k, v in b.as_json().items() if k != "fiber"})
    out["fibers"] = [FIXED_RHO, FIXED_SPECTRUM]
    return out


def incidence_kernel_dim(support_dets):
    """Weight-sector incidence kernel dim: |S| minus rank of the incidence map.

    NOT a fiber dimension on its own. It is the ambient weight kernel K of
    docs/silence_rank_lemma.md, an UPPER bound on the fixed-spectrum weight
    sector before the within-block channel conditions are imposed.
    """
    dets = [tuple(t) for t in support_dets]
    if not dets:
        return 0
    d = max(o for t in dets for o in t) + 1
    a = np.array([[1.0 if o in t else 0.0 for o in range(d)] for t in dets])
    return len(dets) - int(np.linalg.matrix_rank(a))


def amplitudes_from_record(record):
    """(amplitudes, support_dets) from a shipped states.jsonl closed form."""
    import sympy as sp

    cf = record["closed_form"]
    c = np.array([complex(sp.N(sp.sympify(s), 30)) for s in cf["pretty"]])
    dets = [tuple(t) for t in cf["support_dets"]]
    return c, dets
