#!/usr/bin/env python3
"""Symmetry-reduced charge-0 Hermitian-SOS non-attainability certifier (G2/G3).

The un-reduced certifier (scripts/sos_nonattain.py) has a single PSD block of side
1 + nd^2 (401 at (3,6), 27226 at (3,11)); the latter is intractable. When the target
g0 has a spectrum stabilizer G (a product of symmetric groups permuting equal-
occupation modes), f, the sphere, and the SOS are all G-invariant, so the optimal
Gram Q can be taken G-invariant. This module block-diagonalizes it.

METHOD (numerical Wedderburn, character-free):
  * rho_B: the signed-permutation action of G on the charge-0 basis {1} U {Re/Im
    c_a cbar_b} (verified a genuine orthogonal rep).
  * invariant Q lives in the commutant; parametrize Q = sum_o y_o E_o over signed
    G-orbits of index pairs (E_o automatically invariant), with y_o = y_{o^T} so the
    Gram is symmetric.
  * block-diagonalize rho_B by eigendecomposing a random G-invariant operator
    (isotypic split) and linking equal-irrep copies via a second invariant
    (intertwiner alignment); Q >> 0 becomes small blocks Q_nu = T_nu^T Q T_nu >> 0.
  * match monomial coefficients (invariant => f is reproduced), maximize delta.

G2 (this file's __main__): on the S2xS4-symmetric EXTERIOR (3,6) point
diag(0.9,0.9,0.3,0.3,0.3,0.3) the reduced delta reproduces the un-reduced delta
(0.08823, agreement ~1e-6), validating the reduction. Blocks: 10, sizes up to 26.

SCALING NOTE: the dense orbit enumeration and dense invariant-operator averaging
here are O(nz^2 * |G|) and O(nz^2), fine at (3,6) (nz=401) but NOT at (3,11)
(nz=27226); cand 44 needs the matrix-free V-first variant (factorized S5xS6
projectors, see docs/RESEARCH.md). Diagnostic; scripts/ is ruff-excluded.
"""
import sys, itertools, time
sys.path.insert(0, "scripts")
import numpy as np
from collections import defaultdict
from sos_nonattain import charge0_basis, f_coeffs, mul_mono

R2 = np.sqrt(2.0)


def signed_det_action(perm, dets, idx):
    nd = len(dets); tgt = np.empty(nd, int); sgn = np.empty(nd, int)
    for i, T in enumerate(dets):
        img = [perm[x] for x in T]; s = sorted(img)
        inv = sum(1 for a in range(len(img)) for b in range(a+1, len(img)) if img[a] > img[b])
        tgt[i] = idx[tuple(s)]; sgn[i] = (-1)**inv
    return tgt, sgn


def basis_index(nd):
    B = charge0_basis(nd)
    mp = {lab: i for i, (lab, _t) in enumerate(B)}
    return B, mp


def rhoB(perm, dets, idx, nd, mp):
    dtgt, dsgn = signed_det_action(perm, dets, idx)
    nz = 1 + nd*nd; tgt = np.empty(nz, int); sgn = np.empty(nz, int)
    tgt[0] = 0; sgn[0] = 1
    for a in range(nd):
        for b in range(a, nd):
            i = mp[("sym", a, b)]; ta, tb = int(dtgt[a]), int(dtgt[b]); e = int(dsgn[a]*dsgn[b])
            lo, hi = (ta, tb) if ta <= tb else (tb, ta)
            tgt[i] = mp[("sym", lo, hi)]; sgn[i] = e
    for a in range(nd):
        for b in range(a+1, nd):
            i = mp[("asym", a, b)]; ta, tb = int(dtgt[a]), int(dtgt[b]); e = int(dsgn[a]*dsgn[b])
            if ta < tb:
                tgt[i] = mp[("asym", ta, tb)]; sgn[i] = e
            else:
                tgt[i] = mp[("asym", tb, ta)]; sgn[i] = -e
    return tgt, sgn


def group_elements(d, blocks):
    per = [list(itertools.permutations(b)) for b in blocks]
    for combo in itertools.product(*per):
        perm = list(range(d))
        for b, pb in zip(blocks, combo):
            for orig, img in zip(b, pb):
                perm[orig] = img
        yield perm


def build_rho(d, N, blocks):
    dets = [tuple(c) for c in itertools.combinations(range(d), N)]
    idx = {T: i for i, T in enumerate(dets)}; nd = len(dets)
    B, mp = basis_index(nd); nz = 1 + nd*nd
    gels = list(group_elements(d, blocks))
    perms = [rhoB(g, dets, idx, nd, mp) for g in gels]   # list of (tgt,sgn)
    return dets, nd, B, mp, nz, perms


def invariant_orbits(perms, nz):
    """Signed G-orbits of index-pairs -> list of orbit dicts {(P,R): sign_rel}."""
    seen = np.zeros((nz, nz), np.int8)
    orbits = []
    for P in range(nz):
        for R in range(nz):
            if seen[P, R]:
                continue
            # generate orbit of (P,R); track sign relative to this rep (rep sign +1)
            members = {}
            selfkill = False
            for (tgt, sgn) in perms:
                Pp, Rp = int(tgt[P]), int(tgt[R]); s = int(sgn[P]*sgn[R])
                if (Pp, Rp) in members:
                    if members[(Pp, Rp)] != s:
                        selfkill = True
                else:
                    members[(Pp, Rp)] = s
                seen[Pp, Rp] = 1
            if not selfkill:
                orbits.append(members)
    return orbits


def block_diagonalize(perms, nz, seed=0):
    """Return list of (dim, [copy_basis(nz x dim) aligned]) grouping irrep copies."""
    rng = np.random.default_rng(seed)
    def apply_conj(H):   # (1/|G|) sum_g P(g) H P(g)^T
        A = np.zeros((nz, nz))
        for (tgt, sgn) in perms:
            S = (sgn[:, None]*H*sgn[None, :])
            A[np.ix_(tgt, tgt)] += S
        return A/len(perms)
    H1 = rng.standard_normal((nz, nz)); H1 = H1+H1.T
    A = apply_conj(H1)
    w, V = np.linalg.eigh(A)
    # group eigenvectors by eigenvalue
    order = np.argsort(w); w = w[order]; V = V[:, order]
    groups = []; i = 0
    while i < nz:
        j = i+1
        while j < nz and abs(w[j]-w[i]) < 1e-6:
            j += 1
        groups.append(list(range(i, j))); i = j
    eigspaces = [V[:, g] for g in groups]  # each nz x dim, orthonormal, one irrep copy
    # second invariant to link copies of same irrep
    H2 = rng.standard_normal((nz, nz)); H2 = H2+H2.T
    B = apply_conj(H2)
    used = [False]*len(eigspaces); irreps = []
    for k in range(len(eigspaces)):
        if used[k]:
            continue
        Sk = eigspaces[k]; dk = Sk.shape[1]
        copies = [Sk]; used[k] = True
        for l in range(k+1, len(eigspaces)):
            if used[l]:
                continue
            Sl = eigspaces[l]
            if Sl.shape[1] != dk:
                continue
            M = Sk.T @ B @ Sl   # intertwiner dk x dk; nonzero & ~scalar*orthogonal iff same irrep
            # same irrep: M = c * Rotation (M M^T = c^2 I)
            MMt = M @ M.T
            c2 = np.trace(MMt)/dk
            if c2 > 1e-8 and np.allclose(MMt, c2*np.eye(dk), atol=1e-5*max(c2, 1)):
                # align Sl's internal basis to Sk via the isometry M/sqrt(c2)
                Q = M/np.sqrt(c2)                     # ~ orthogonal
                copies.append(Sl @ Q.T)               # now aligned to Sk internal indices
                used[l] = True
        irreps.append((dk, copies))
    return irreps


def reduced_certify(d, N, g0, blocks, solver="SCS", eps=1e-8):
    import cvxpy as cp
    dets, nd, B, mp, nz, perms = build_rho(d, N, blocks)
    _nd, F = f_coeffs(d, N, g0)
    irreps = block_diagonalize(perms, nz)
    orbits = invariant_orbits(perms, nz)
    norb = len(orbits)
    # transpose-orbit pairing (Gram must be symmetric: y_o == y_{o^T})
    pair2orb = {}
    for oi, members in enumerate(orbits):
        for (P, Rr) in members:
            pair2orb[(P, Rr)] = oi
    transp = np.arange(norb)
    for oi, members in enumerate(orbits):
        (P, Rr) = next(iter(members))
        transp[oi] = pair2orb[(Rr, P)]
    # orbit reps (one (P,R) per orbit, its sign is +1 by construction)
    # map: contribution of y_o to monomial key, and to each block Q_nu entry
    # (A) monomial coeff contributions of each orbit
    keymap = defaultdict(lambda: np.zeros(norb, complex))
    for oi, members in enumerate(orbits):
        for (P, Rr), s in members.items():
            # E_o entry at (P,R) = s ; sigma gets sum_{P,R} Q[P,R] B_P B_R
            for t1 in B[P][1]:
                for t2 in B[Rr][1]:
                    h, a, cc = mul_mono(t1, t2)
                    keymap[(h, a)][oi] += s*cc
    # (B) block contributions: Q_nu[a,b] = sum_o y_o * (Tnu^T E_o Tnu)[a,b]
    #     with E_o = sum_{(P,R)} s e_P e_R^T ; (Tnu^T E_o Tnu)[a,b] = sum s Tnu[P,a] Tnu[R,b]
    #     copies give Tnu columns; block index = copy index, internal index j shared.
    # Build per-irrep: list of copy-matrices (nz x dim). Reduced var Q_nu is (ncopy x ncopy).
    blockinfo = []
    for (dim, copies) in irreps:
        m = len(copies)
        blockinfo.append((dim, copies, m))
    # orbit->block linear maps: coeff_block[nu][a,b] as vector over orbits, per internal j (sum over j)
    # (Tnu^T E_o Tnu) summed appropriately: invariant Q acts equally on each internal index,
    # so Q_nu[a,b] (copy a,b) = (1/dim) sum_j sum_{(P,R),s} s copies[a][P,j] copies[b][R,j]
    block_orbit = []   # block_orbit[nu][a][b] = vector length norb
    for (dim, copies, m) in blockinfo:
        C = np.stack(copies, 0)          # m x nz x dim
        BO = np.zeros((m, m, norb))
        for oi, members in enumerate(orbits):
            # sum_{(P,R),s} s * C[a,P,:] . C[b,R,:]   (dot over internal j), /dim
            acc = np.zeros((m, m))
            for (P, Rr), s in members.items():
                acc += s*(C[:, P, :] @ C[:, Rr, :].T)   # m x m
            BO[:, :, oi] = acc/dim
        block_orbit.append(BO)
    y = cp.Variable(norb)
    # multiplier mu invariant: mu = sum over invariant basis vectors. Invariant basis of B =
    # G-orbit sums of single basis elements (trivial isotypic of the basis permutation action).
    # orbits of single indices:
    seen = np.zeros(nz, np.int8); bas_orbits = []
    for P in range(nz):
        if seen[P]:
            continue
        mem = {}
        kill = False
        for (tgt, sgn) in perms:
            Pp = int(tgt[P]); s = int(sgn[P])
            if Pp in mem:
                if mem[Pp] != s:
                    kill = True
            else:
                mem[Pp] = s
            seen[Pp] = 1
        if not kill:
            bas_orbits.append(mem)
    nmu = len(bas_orbits)
    mu = cp.Variable(nmu)
    # mu(||c||^2 - 1): for each basis-orbit (invariant combo v = sum s e_P),
    #   contribution to keys: v * (sum_k c_k cbar_k) - v
    mukey = defaultdict(lambda: np.zeros(nmu, complex))
    for mi, mem in enumerate(bas_orbits):
        for P, s in mem.items():
            for (h0, a0, c0) in B[P][1]:
                for k in range(nd):
                    h, a, cc = mul_mono((h0, a0, c0), ((k,), (k,), 1.0+0j))
                    mukey[(h, a)][mi] += s*cc
                mukey[(h0, a0)][mi] += s*(-c0)
    # assemble constraints (VECTORIZED). Q_nu = reshape(BOmat @ y) must be PSD.
    cons = []
    for i, BO in enumerate(block_orbit):
        m = BO.shape[0]
        BOmat = BO.reshape(m*m, norb)
        Qexpr = cp.reshape(BOmat @ y, (m, m), order="C")
        cons.append(Qexpr >> 0)
    delta = cp.Variable()
    allkeys = sorted(set(keymap) | set(mukey) | set(F))
    nk = len(allkeys)
    Kre = np.zeros((nk, norb)); Kim = np.zeros((nk, norb))
    Mre = np.zeros((nk, nmu)); Mim = np.zeros((nk, nmu))
    dvec = np.zeros(nk); rre = np.zeros(nk); rim = np.zeros(nk)
    for r, key in enumerate(allkeys):
        if key in keymap:
            Kre[r] = np.real(keymap[key]); Kim[r] = np.imag(keymap[key])
        if key in mukey:
            Mre[r] = np.real(mukey[key]); Mim[r] = np.imag(mukey[key])
        if key == ((), ()):
            dvec[r] = 1.0
        fv = F.get(key, 0); rre[r] = float(np.real(fv)); rim[r] = float(np.imag(fv))
    cons.append(Kre @ y + Mre @ mu + dvec*delta == rre)
    cons.append(Kim @ y + Mim @ mu == rim)
    # Gram symmetry: y_o == y_{o^T}
    for oi in range(norb):
        if transp[oi] > oi:
            cons.append(y[oi] == y[transp[oi]])
    prob = cp.Problem(cp.Maximize(delta), cons)
    t0 = time.time()
    prob.solve(solver=getattr(cp, solver), verbose=False)
    return dict(status=prob.status, delta=(None if delta.value is None else float(delta.value)),
                nz=nz, norb=norb, nblocks=len(blockinfo),
                blocksizes=[b[2] for b in blockinfo], secs=time.time()-t0)


def _unreduced(d, N, g0):
    from sos_nonattain import certify
    return certify(d, N, g0, eps=1e-6, max_iters=30000)["delta"]


if __name__ == "__main__":
    d, N = 6, 3
    blocks = [[0, 1], [2, 3, 4, 5]]
    g0 = [0.9, 0.9, 0.3, 0.3, 0.3, 0.3]   # S2xS4-symmetric EXTERIOR (violates pairing)
    print("G2 gate: symmetry-reduced SOS must reproduce the un-reduced delta.\n")
    r = reduced_certify(d, N, g0, blocks)
    du = _unreduced(d, N, g0)
    print(f"REDUCED   (3,6)/S2xS4: {r['status']} delta = {r['delta']:.8g}  "
          f"({r['nblocks']} blocks {r['blocksizes']}, {r['secs']:.0f}s)")
    print(f"UNREDUCED (3,6)      :          delta = {du:.8g}")
    print(f"|diff| = {abs(r['delta']-du):.2e}  -> "
          f"{'G2 PASS' if abs(r['delta']-du) < 1e-4 else 'G2 FAIL'}")
