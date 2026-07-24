#!/usr/bin/env python3
"""Scalable symmetry-reduced SOS for (3,11) cand 44 (G3).

The dense G2 machinery (scripts/sos_symmetry.py) is O(nz^2 |G|) in the block-
diagonalization and enumerates all nz^2 index pairs; both die at (3,11) where the
full charge-0 basis has nz = 1 + 165^2 = 27226 and |G| = |S5 x S6| = 86400. Three
scalable swaps make cand 44 reachable, each validated against G2 on (3,6)/S2xS4:

  1. SUPPORT RESTRICTION. Restrict the SOS Gram to the G-invariant charge-0
     sub-basis {1} U {diagonal} U {one-hop pairs} (nz ~ 4126). Restricting the
     Gram support only shrinks the SOS cone, so a PSD certificate here is a VALID
     global one -- delta > 0 remains rigorous (it is a term-sparsity pattern; if
     delta <= 0 the pattern is too tight and must be grown).
  2. FACTORIZED REYNOLDS. For G = A x B acting on DISJOINT modes (so A,B commute),
     Reynolds_G = Reynolds_A o Reynolds_B: |A| + |B| = 840 conjugations, not
     |A|*|B| = 86400. Verified to reproduce the dense average to 1e-15.
  3. GENERATOR-BFS ORBITS. Signed G-orbits of index pairs via BFS over the group
     generators (adjacent transpositions), O(nz^2 * gens) not O(nz^2 * |G|).

Everything downstream (invariant Gram = sum_o y_o E_o with y_o = y_{o^T},
coefficient matching, block PSD) is the validated G2 assembly.

Run: .venv/bin/python scripts/sos_symmetry_scaled.py [--cand44 | --test36]
Diagnostic; scripts/ is ruff-excluded.
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import defaultdict

import numpy as np

sys.path.insert(0, "scripts")
from sos_nonattain import f_coeffs, mul_mono  # noqa: E402

R2 = np.sqrt(2.0)


def hop_pairs(dets, max_hop):
    """unordered det-index pairs (a<b) with hop distance N-|a&b| <= max_hop."""
    nd = len(dets)
    N = len(dets[0])
    pairs = set()
    for i in range(nd):
        si = set(dets[i])
        for j in range(i + 1, nd):
            if N - len(si & set(dets[j])) <= max_hop:
                pairs.add((i, j))
    return pairs


def restricted_basis(dets, max_hop=1):
    """Charge-0 basis restricted to const + diagonal + pairs within max_hop.
    max_hop >= N gives the full charge-0 basis. Returns (labels, mp)."""
    nd = len(dets)
    oh = sorted(hop_pairs(dets, max_hop))
    labels = [("const", [((), (), 1.0 + 0j)])]
    for a in range(nd):                       # diagonal occupations |c_a|^2
        labels.append((("sym", a, a), [((a,), (a,), 1.0 + 0j)]))
    for (a, b) in oh:                          # Re(c_a cbar_b), one-hop
        labels.append((("sym", a, b), [((a,), (b,), 0.5 + 0j), ((b,), (a,), 0.5 + 0j)]))
    for (a, b) in oh:                          # Im(c_a cbar_b), one-hop
        labels.append((("asym", a, b), [((a,), (b,), -0.5j), ((b,), (a,), 0.5j)]))
    mp = {lab: i for i, (lab, _t) in enumerate(labels)}
    return labels, mp


def signed_det_action(perm, dets, idx):
    nd = len(dets); tgt = np.empty(nd, int); sgn = np.empty(nd, int)
    for i, T in enumerate(dets):
        img = [perm[x] for x in T]; s = sorted(img)
        inv = sum(1 for a in range(len(img)) for b in range(a + 1, len(img)) if img[a] > img[b])
        tgt[i] = idx[tuple(s)]; sgn[i] = (-1) ** inv
    return tgt, sgn


def rho_restricted(perm, dets, idx, labels, mp):
    """Signed action of a mode-permutation on the restricted basis: (tgt,sgn)."""
    dtgt, dsgn = signed_det_action(perm, dets, idx)
    nz = len(labels); tgt = np.empty(nz, int); sgn = np.empty(nz, int)
    for k, (lab, _t) in enumerate(labels):
        if lab == "const":
            tgt[k] = mp["const"]; sgn[k] = 1; continue
        kind, a, b = lab
        ta, tb = int(dtgt[a]), int(dtgt[b]); e = int(dsgn[a] * dsgn[b])
        if kind == "sym":
            lo, hi = (ta, tb) if ta <= tb else (tb, ta)
            tgt[k] = mp[("sym", lo, hi)]; sgn[k] = e
        else:
            if ta < tb:
                tgt[k] = mp[("asym", ta, tb)]; sgn[k] = e
            else:
                tgt[k] = mp[("asym", tb, ta)]; sgn[k] = -e
    return tgt, sgn


def factor_actions(d, N, factor_blocks, all_elements=True, max_hop=1):
    """Per-factor list of restricted-basis actions. all_elements: full factor group
    (for Reynolds); else adjacent-transposition generators (for BFS)."""
    dets = [tuple(c) for c in itertools.combinations(range(d), N)]
    idx = {T: i for i, T in enumerate(dets)}
    labels, mp = restricted_basis(dets, max_hop)
    out = []
    for fb in factor_blocks:
        acts = []
        if all_elements:
            gens = list(itertools.permutations(fb))
        else:
            gens = []
            for t in range(len(fb) - 1):     # adjacent transpositions
                pb = list(fb); pb[t], pb[t + 1] = pb[t + 1], pb[t]; gens.append(tuple(pb))
        for pb in gens:
            perm = list(range(d))
            for orig, img in zip(fb, pb):
                perm[orig] = img
            acts.append(rho_restricted(perm, dets, idx, labels, mp))
    # note: acts appended per element; regroup below
        out.append(acts)
    return dets, labels, mp, out


def reynolds_factorized(H, factor_perms):
    A = H
    for perms in factor_perms:
        A2 = np.zeros_like(A)
        for (tgt, sgn) in perms:
            A2[np.ix_(tgt, tgt)] += sgn[:, None] * A * sgn[None, :]
        A = A2 / len(perms)
    return A


def block_diagonalize(nz, factor_perms, seed=0):
    """Isotypic split + intertwiner alignment using factorized Reynolds."""
    rng = np.random.default_rng(seed)
    H1 = rng.standard_normal((nz, nz)); H1 = H1 + H1.T
    A = reynolds_factorized(H1, factor_perms)
    w, V = np.linalg.eigh(A)
    order = np.argsort(w); w = w[order]; V = V[:, order]
    groups = []; i = 0
    while i < nz:
        j = i + 1
        while j < nz and abs(w[j] - w[i]) < 1e-6:
            j += 1
        groups.append(list(range(i, j))); i = j
    eig = [V[:, g] for g in groups]
    H2 = rng.standard_normal((nz, nz)); H2 = H2 + H2.T
    B = reynolds_factorized(H2, factor_perms)
    used = [False] * len(eig); irreps = []
    for k in range(len(eig)):
        if used[k]:
            continue
        Sk = eig[k]; dk = Sk.shape[1]; copies = [Sk]; used[k] = True
        for l in range(k + 1, len(eig)):
            if used[l] or eig[l].shape[1] != dk:
                continue
            M = Sk.T @ B @ eig[l]; MMt = M @ M.T; c2 = np.trace(MMt) / dk
            if c2 > 1e-8 and np.allclose(MMt, c2 * np.eye(dk), atol=1e-5 * max(c2, 1)):
                copies.append(eig[l] @ (M / np.sqrt(c2)).T); used[l] = True
        irreps.append((dk, copies))
    return irreps


def bfs_orbits(nz, gen_perms):
    """Signed G-orbits of ordered index pairs via BFS over generators."""
    seen = np.zeros((nz, nz), np.int8)
    orbits = []
    for P0 in range(nz):
        for R0 in range(nz):
            if seen[P0, R0]:
                continue
            members = {(P0, R0): 1}; frontier = [(P0, R0, 1)]; seen[P0, R0] = 1
            kill = False
            while frontier:
                (P, R, s) = frontier.pop()
                for (tgt, sgn) in gen_perms:
                    Pp, Rp = int(tgt[P]), int(tgt[R]); s2 = s * int(sgn[P] * sgn[R])
                    if (Pp, Rp) in members:
                        if members[(Pp, Rp)] != s2:
                            kill = True
                    else:
                        members[(Pp, Rp)] = s2; seen[Pp, Rp] = 1; frontier.append((Pp, Rp, s2))
            if not kill:
                orbits.append(members)
    return orbits


def det_generators(d, N, factor_blocks):
    """Adjacent-transposition generators as (dtgt,dsgn) on the determinant basis."""
    dets = [tuple(c) for c in itertools.combinations(range(d), N)]
    idx = {T: i for i, T in enumerate(dets)}
    gens = []
    for fb in factor_blocks:
        for t in range(len(fb) - 1):
            perm = list(range(d)); perm[fb[t]], perm[fb[t + 1]] = fb[t + 1], fb[t]
            gens.append(signed_det_action(perm, dets, idx))
    return dets, gens


def key_image(key, dtgt):
    (h, a) = key
    return (tuple(sorted(int(dtgt[x]) for x in h)), tuple(sorted(int(dtgt[x]) for x in a)))


def canonize_keys(keys, det_gens):
    """Group keys into G-orbits (ignoring sign); return rep_of[key] and reps list."""
    rep_of = {}
    reps = []
    for k0 in keys:
        if k0 in rep_of:
            continue
        # BFS the orbit
        seen = {k0}; fr = [k0]
        while fr:
            k = fr.pop()
            for (dtgt, _s) in det_gens:
                k2 = key_image(k, dtgt)
                if k2 not in seen:
                    seen.add(k2); fr.append(k2)
        rep = min(seen)
        for k in seen:
            rep_of[k] = rep
        reps.append(rep)
    return rep_of, reps


def reduced_certify_direct(d, N, g0, factor_blocks, solver="SCS", eps=1e-5,
                           max_iters=40000, verbose=False, max_hop=1):
    """Direct block formulation: variables are the PSD blocks Q_nu (no orbit y),
    constraints matched per monomial-key ORBIT. sigma coeff of a key comes from the
    FEW basis pairs (P,Q) producing it -- no nz^2 member sum. sig built by nz^2 pass
    (fine for validation sizes)."""
    import cvxpy as cp
    t0 = time.time()
    dets, labels, mp, fperms = factor_actions(d, N, factor_blocks, all_elements=True, max_hop=max_hop)
    nz = len(labels); B = labels
    _nd, F = f_coeffs(d, N, g0)
    _dg, det_gens = det_generators(d, N, factor_blocks)
    tprep = time.time() - t0

    _d2, _l2, _m2, gens = factor_actions(d, N, factor_blocks, all_elements=False, max_hop=max_hop)
    genflat = [a for grp in gens for a in grp]
    t = time.time(); irreps = block_diagonalize(nz, fperms); tbd = time.time() - t

    # scalable sig: derive each key's (P,Q) via a det-pair lookup, no nz^2 pass.
    t = time.time()
    dpmap = defaultdict(list)                       # (h,a) content -> [(basis_idx, coeff)]
    for P in range(nz):
        for (h, a, c) in B[P][1]:
            dpmap[(h, a)].append((P, c))

    def _splits(seq):                              # distinct ordered (part_P, part_Q)
        n = len(seq); out = set()
        for m in range(1 << n):
            out.add((tuple(sorted(seq[i] for i in range(n) if not (m >> i) & 1)),
                     tuple(sorted(seq[i] for i in range(n) if (m >> i) & 1))))
        return out

    def key_pqlist(key):
        """all (P,Q,coeff) with B_P B_Q producing this key, via dpmap splits. Splits
        are DEDUPED so repeated indices (e.g. c_A^2) are not double-counted."""
        holo, anti = key
        acc = defaultdict(complex)
        for (hP, hQ) in _splits(holo):
            for (aP, aQ) in _splits(anti):
                for (P, cP) in dpmap.get((hP, aP), ()):
                    for (Q, cQ) in dpmap.get((hQ, aQ), ()):
                        acc[(P, Q)] += cP * cQ
        return [(P, Q, c) for (P, Q), c in acc.items()]

    # enumerate candidate keys: basis-orbit reps x all Q (covers every key orbit),
    # plus f and multiplier keys; then canonicalize to orbit reps.
    seenu = np.zeros(nz, np.int8); brep = []
    for P in range(nz):
        if seenu[P]:
            continue
        brep.append(P); fr = [P]; seenu[P] = 1
        while fr:
            Qi = fr.pop()
            for (tgt, _s) in genflat:
                Qp = int(tgt[Qi])
                if not seenu[Qp]:
                    seenu[Qp] = 1; fr.append(Qp)
    cand = set()
    for P0 in brep:
        for Q in range(nz):
            for t1 in B[P0][1]:
                for t2 in B[Q][1]:
                    h, a, _cc = mul_mono(t1, t2)
                    cand.add((h, a))
    cand |= set(F)
    tsig = time.time() - t

    # multiplier mu over single-index G-orbits (cheap); mukey[key] -> vec(nmu)
    seen = np.zeros(nz, np.int8); bas_orb = []
    for P in range(nz):
        if seen[P]:
            continue
        mem = {P: 1}; fr = [(P, 1)]; seen[P] = 1; kill = False
        while fr:
            (Qi, s) = fr.pop()
            for (tgt, sgn) in genflat:
                Qp = int(tgt[Qi]); s2 = s * int(sgn[Qi])
                if Qp in mem:
                    if mem[Qp] != s2:
                        kill = True
                else:
                    mem[Qp] = s2; seen[Qp] = 1; fr.append((Qp, s2))
        if not kill:
            bas_orb.append(mem)
    nmu = len(bas_orb)
    mukey = defaultdict(lambda: np.zeros(nmu, complex))
    for mi, mem in enumerate(bas_orb):
        for P, s in mem.items():
            for (h0, a0, c0) in B[P][1]:
                for k in range(len(dets)):
                    h, a, cc = mul_mono((h0, a0, c0), ((k,), (k,), 1.0 + 0j))
                    mukey[(h, a)][mi] += s * cc
                mukey[(h0, a0)][mi] += s * (-c0)

    # key orbits
    allkeys = cand | set(mukey)
    rep_of, reps = canonize_keys(allkeys, det_gens)

    # variables and per-block copy tensors
    import scipy.sparse as sp
    Cs = [np.stack(copies, 0) for (dim, copies) in irreps]   # each m x nz x d
    ms = [C.shape[0] for C in Cs]
    off = np.cumsum([0] + [m * m for m in ms])              # column offsets per block
    Qdim = int(off[-1])
    Qs = [cp.Variable((m, m), symmetric=True) for m in ms]
    mu = cp.Variable(nmu); delta = cp.Variable()
    cons = [Qs[i] >> 0 for i in range(len(Qs))]

    # VECTORIZED coefficient matching. Build sparse maps: qvec (Sum m^2) | mu | delta.
    # For each rep row: sum_i <Gr_i, Q_i> + mu.Re(mk) + delta[empty] = Re f ; imag row too.
    qr_r, qr_c, qr_v = [], [], []   # real rows vs qvec
    qi_r, qi_c, qi_v = [], [], []   # imag rows vs qvec
    Mr = np.zeros((len(reps), nmu)); Mi = np.zeros((len(reps), nmu))
    dvec = np.zeros(len(reps)); rre = np.zeros(len(reps)); rim = np.zeros(len(reps))
    for r, rep in enumerate(reps):
        pql = key_pqlist(rep)
        for i, C in enumerate(Cs):
            m = ms[i]; G = np.zeros((m, m), complex)
            for (P, Q, cc) in pql:
                G += cc * (C[:, P, :] @ C[:, Q, :].T)
            Gr = np.real(G); Gi = np.imag(G)
            fr = Gr.flatten(order="F"); fi = Gi.flatten(order="F")  # match cp.vec col-major
            nzr = np.nonzero(np.abs(fr) > 1e-12)[0]
            nzi = np.nonzero(np.abs(fi) > 1e-12)[0]
            for k in nzr:
                qr_r.append(r); qr_c.append(off[i] + k); qr_v.append(fr[k])
            for k in nzi:
                qi_r.append(r); qi_c.append(off[i] + k); qi_v.append(fi[k])
        mk = mukey.get(rep, None)
        if mk is not None:
            Mr[r] = np.real(mk); Mi[r] = np.imag(mk)
        if rep == ((), ()):
            dvec[r] = 1.0
        fv = F.get(rep, 0); rre[r] = float(np.real(fv)); rim[r] = float(np.imag(fv))
    nR = len(reps)
    Aqr = sp.csr_matrix((qr_v, (qr_r, qr_c)), shape=(nR, Qdim))
    Aqi = sp.csr_matrix((qi_v, (qi_r, qi_c)), shape=(nR, Qdim))
    qvec = cp.hstack([cp.vec(Q) for Q in Qs])
    cons.append(Aqr @ qvec + Mr @ mu + cp.multiply(dvec, delta) == rre)
    cons.append(Aqi @ qvec + Mi @ mu == rim)

    prob = cp.Problem(cp.Maximize(delta), cons)
    tasm = time.time()
    if solver == "SCS":
        prob.solve(solver=cp.SCS, verbose=verbose, eps=eps, max_iters=max_iters)
    else:
        prob.solve(solver=getattr(cp, solver), verbose=verbose)
    return dict(status=prob.status, delta=(None if delta.value is None else float(delta.value)),
                nz=nz, nreps=len(reps), nmu=nmu, blocksizes=[C.shape[0] for C in Cs],
                tprep=tprep, tbd=tbd, tsig=tsig, tsolve=time.time() - tasm)


def reduced_certify(d, N, g0, factor_blocks, solver="SCS", eps=1e-5, max_iters=40000,
                    verbose=False, max_hop=1):
    import cvxpy as cp
    t0 = time.time()
    dets, labels, mp, fperms = factor_actions(d, N, factor_blocks, all_elements=True, max_hop=max_hop)
    _dets2, _l2, _m2, gens = factor_actions(d, N, factor_blocks, all_elements=False, max_hop=max_hop)
    genflat = [a for grp in gens for a in grp]
    nz = len(labels)
    B = labels
    _nd, F = f_coeffs(d, N, g0)
    tprep = time.time() - t0

    t = time.time()
    irreps = block_diagonalize(nz, fperms)
    tbd = time.time() - t
    t = time.time()
    orbits = bfs_orbits(nz, genflat)
    torb = time.time() - t
    norb = len(orbits)

    pair2orb = {}
    for oi, mem in enumerate(orbits):
        for pr in mem:
            pair2orb[pr] = oi
    transp = np.arange(norb)
    for oi, mem in enumerate(orbits):
        (P, Rr) = next(iter(mem)); transp[oi] = pair2orb[(Rr, P)]

    keymap = defaultdict(lambda: np.zeros(norb, complex))
    for oi, mem in enumerate(orbits):
        for (P, Rr), s in mem.items():
            for t1 in B[P][1]:
                for t2 in B[Rr][1]:
                    h, a, cc = mul_mono(t1, t2)
                    keymap[(h, a)][oi] += s * cc

    blockinfo = [(dim, copies, len(copies)) for (dim, copies) in irreps]
    block_orbit = []
    for (dim, copies, m) in blockinfo:
        C = np.stack(copies, 0)
        BO = np.zeros((m, m, norb))
        for oi, mem in enumerate(orbits):
            acc = np.zeros((m, m))
            for (P, Rr), s in mem.items():
                acc += s * (C[:, P, :] @ C[:, Rr, :].T)
            BO[:, :, oi] = acc / dim
        block_orbit.append(BO)

    # invariant multiplier orbits (single-index signed G-orbits)
    seen = np.zeros(nz, np.int8); bas_orb = []
    for P in range(nz):
        if seen[P]:
            continue
        mem = {P: 1}; fr = [(P, 1)]; seen[P] = 1; kill = False
        while fr:
            (Q, s) = fr.pop()
            for (tgt, sgn) in genflat:
                Qp = int(tgt[Q]); s2 = s * int(sgn[Q])
                if Qp in mem:
                    if mem[Qp] != s2:
                        kill = True
                else:
                    mem[Qp] = s2; seen[Qp] = 1; fr.append((Qp, s2))
        if not kill:
            bas_orb.append(mem)
    nmu = len(bas_orb)
    mukey = defaultdict(lambda: np.zeros(nmu, complex))
    for mi, mem in enumerate(bas_orb):
        for P, s in mem.items():
            for (h0, a0, c0) in B[P][1]:
                for k in range(len(dets)):
                    kk = ("sym", k, k)
                    # ||c||^2 = sum_k |c_k|^2 ; multiply by that diagonal monomial
                    h, a, cc = mul_mono((h0, a0, c0), ((k,), (k,), 1.0 + 0j))
                    mukey[(h, a)][mi] += s * cc
                mukey[(h0, a0)][mi] += s * (-c0)

    y = cp.Variable(norb); mu = cp.Variable(nmu); delta = cp.Variable()
    cons = []
    for i, BO in enumerate(block_orbit):
        m = BO.shape[0]
        cons.append(cp.reshape(BO.reshape(m * m, norb) @ y, (m, m), order="C") >> 0)
    for oi in range(norb):
        if transp[oi] > oi:
            cons.append(y[oi] == y[transp[oi]])
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
    cons.append(Kre @ y + Mre @ mu + dvec * delta == rre)
    cons.append(Kim @ y + Mim @ mu == rim)
    prob = cp.Problem(cp.Maximize(delta), cons)
    tasm = time.time()
    if solver == "SCS":
        prob.solve(solver=cp.SCS, verbose=verbose, eps=eps, max_iters=max_iters)
    else:
        prob.solve(solver=getattr(cp, solver), verbose=verbose)
    return dict(status=prob.status, delta=(None if delta.value is None else float(delta.value)),
                nz=nz, norb=norb, nmu=nmu, nk=nk,
                blocksizes=[b[2] for b in blockinfo],
                tprep=tprep, tbd=tbd, torb=torb, tsolve=time.time() - tasm)


def main():
    if "--cand44" in sys.argv:
        d, N = 11, 3
        g0 = [6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1]; g0 = [x / 12 for x in g0]
        fb = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9, 10]]
        tag = "cand44 (3,11) S5xS6"
    else:
        d, N = 6, 3
        g0 = [0.9, 0.9, 0.3, 0.3, 0.3, 0.3]
        fb = [[0, 1], [2, 3, 4, 5]]
        tag = "test (3,6)/S2xS4 restricted"
    mh = 2 if "--cand44" in sys.argv else 1
    for a in sys.argv:
        if a.startswith("--hop="):
            mh = int(a.split("=")[1])
    fn = reduced_certify_direct if "--direct" in sys.argv else reduced_certify
    r = fn(d, N, g0, fb, verbose="--v" in sys.argv, max_hop=mh)
    tag += f" [hop<={mh}, {'direct' if '--direct' in sys.argv else 'orbit'}]"
    print(f"{tag}: {r['status']} delta = {r['delta']}", flush=True)
    print(f"  nz {r['nz']} blocks {r['blocksizes']} "
          f"reps/orbits {r.get('nreps', r.get('norb'))}", flush=True)
    print(f"  times: prep {r['tprep']:.1f}s blockdiag {r['tbd']:.1f}s "
          f"{'sig' if 'tsig' in r else 'orbits'} "
          f"{r.get('tsig', r.get('torb')):.1f}s solve {r['tsolve']:.1f}s", flush=True)
    return r
    print(f"  nz {r['nz']} orbits {r['norb']} mu {r['nmu']} keys {r['nk']} "
          f"blocks {r['blocksizes']}")
    print(f"  times: prep {r['tprep']:.1f}s blockdiag {r['tbd']:.1f}s orbits {r['torb']:.1f}s "
          f"solve {r['tsolve']:.1f}s")


if __name__ == "__main__":
    main()
