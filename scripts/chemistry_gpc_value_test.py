#!/usr/bin/env python3
"""Falsification test: is GPC geometry useful for quantum chemistry?

The question this script is built to answer is narrow and adversarial. GPC
geometry tells us which one-body spectra are possible. Quasipinning says a
physical spectrum often sits close to a polytope facet. The commercial claim
that follows would be: closeness to a facet certifies that the wavefunction
lives in the small determinant subspace that EXACT pinning would force, so a
chemist gets either a smaller CI space or a certificate on an active space.

Every piece of that chain is measured here rather than assumed:

  Experiment 0  the selection rule itself (pinned => determinant subspace),
                which turns out to need a nondegeneracy hypothesis that
                closed-shell chemistry violates exactly;
  Experiment 1  dimension, omitted weight and energy error of the GPC-selected
                subspace on real and model Hamiltonians;
  Experiment 2  GPC subspace against equal-dimension generic truncations;
  Experiment 3  the empirical relation between quasipinning distance D and
                omitted weight W_out, including a directed counterexample hunt;
  Experiment 4  what rigorous energy bound W_out actually buys;
  Experiment 5  how much of any reduction is already given by particle number,
                Sz, S^2, point group and saturated-occupation freezing.

CONVENTIONS. Spin orbitals are 0-indexed; determinants are sorted tuples of
occupied spin-orbital indices. Two-electron integrals are chemist ordered,

    V = (1/2) sum_{pqrs} g[p,q,r,s] a_p^dag a_r^dag a_s a_q,

matching gpc_census.pinned_physics. The one-body density matrix is
gamma[p,q] = <Psi| a_p^dag a_q |Psi>. GPC rows are stored as
``coeffs . lambda <= rhs`` with lambda sorted DESCENDING, so the slack is
``D = rhs - coeffs . lambda >= 0``.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time
from dataclasses import dataclass, field

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(REPO, "src") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "src"))

from gpc_census import constraints as gpc_constraints  # noqa: E402

# ---------------------------------------------------------------------------
# determinant algebra
# ---------------------------------------------------------------------------


def determinants(d: int, n: int) -> list[tuple[int, ...]]:
    return list(itertools.combinations(range(d), n))


class DetBasis:
    """The full C(d,n) determinant basis with cached one-body excitation maps."""

    _cache: dict[tuple[int, int], "DetBasis"] = {}

    def __new__(cls, d: int, n: int):
        key = (d, n)
        if key not in cls._cache:
            obj = super().__new__(cls)
            obj._init(d, n)
            cls._cache[key] = obj
        return cls._cache[key]

    def _init(self, d: int, n: int) -> None:
        self.d = d
        self.n = n
        self.dets = determinants(d, n)
        self.index = {D: i for i, D in enumerate(self.dets)}
        self.m = len(self.dets)
        self.occ = np.zeros((self.m, d), dtype=bool)
        for i, D in enumerate(self.dets):
            self.occ[i, list(D)] = True
        # a_p^dag a_q acting on each determinant: target index and sign
        tgt = np.full((d, d, self.m), -1, dtype=np.int32)
        sgn = np.zeros((d, d, self.m), dtype=np.int8)
        for i, D in enumerate(self.dets):
            occ = set(D)
            for q in D:
                rest = sorted(occ - {q})
                s_q = (-1) ** D.index(q)
                for p in range(d):
                    if p in rest:
                        continue
                    pos = sum(1 for o in rest if o < p)
                    new = tuple(sorted(rest + [p]))
                    tgt[p, q, i] = self.index[new]
                    sgn[p, q, i] = s_q * ((-1) ** pos)
        self.tgt = tgt
        self.sgn = sgn
        self._gen = None

    def generators(self) -> np.ndarray:
        """All a_p^dag a_q as one dense (d, d, m, m) array, built once."""
        if getattr(self, "_gen", None) is None:
            g = np.zeros((self.d, self.d, self.m, self.m))
            for p in range(self.d):
                for q in range(self.d):
                    g[p, q] = self.epq(p, q)
            self._gen = g
        return self._gen

    def epq(self, p: int, q: int) -> np.ndarray:
        """Dense matrix of a_p^dag a_q in the determinant basis."""
        mat = np.zeros((self.m, self.m))
        t = self.tgt[p, q]
        live = t >= 0
        mat[t[live], np.nonzero(live)[0]] = self.sgn[p, q][live]
        return mat

    # -- symmetry sectors -------------------------------------------------
    def sector(self, labels: np.ndarray, value) -> np.ndarray:
        """Indices of determinants whose summed orbital label equals ``value``."""
        tot = self.occ @ np.asarray(labels)
        return np.nonzero(np.isclose(tot, value))[0]


# ---------------------------------------------------------------------------
# Hamiltonian in a spin-orbital basis
# ---------------------------------------------------------------------------


@dataclass
class SpinOrbitalHamiltonian:
    h1: np.ndarray                      # (d, d)
    g2: np.ndarray                      # (d, d, d, d), chemist ordering
    e_core: float = 0.0
    label: str = ""

    @property
    def d(self) -> int:
        return self.h1.shape[0]

    def matrix(self, basis: DetBasis) -> np.ndarray:
        """Full CI matrix on ``basis``.

        Built from the one-body generators through
        ``a_p^dag a_r^dag a_s a_q = E_pq E_rs - delta_qr E_ps``, so every
        fermionic sign is fixed once inside E_pq and no Slater-Condon case
        analysis can get one wrong.
        """
        m = basis.m
        cplx = np.iscomplexobj(self.h1) or np.iscomplexobj(self.g2)
        dt = complex if cplx else float
        E = basis.generators().astype(dt)
        H = np.einsum("pq,pqij->ij", self.h1.astype(dt), E, optimize=True)
        G = np.einsum("pqrs,rsij->pqij", self.g2.astype(dt), E, optimize=True)
        H += 0.5 * np.einsum("pqik,pqkj->ij", E, G, optimize=True)
        contract = np.einsum("prrs->ps", self.g2).astype(dt)
        H -= 0.5 * np.einsum("ps,psij->ij", contract, E, optimize=True)
        H += self.e_core * np.eye(m, dtype=dt)
        return H

    def rotate(self, C: np.ndarray) -> "SpinOrbitalHamiltonian":
        """Integrals in the orbital basis whose creation ops are b_k^dag =
        sum_p C[p,k] a_p^dag."""
        h1 = C.conj().T @ self.h1 @ C
        g2 = np.einsum("pi,qj,rk,sl,pqrs->ijkl",
                       C.conj(), C, C.conj(), C, self.g2, optimize=True)
        return SpinOrbitalHamiltonian(h1, g2, self.e_core, self.label)


def one_rdm(vec: np.ndarray, basis: DetBasis) -> np.ndarray:
    d = basis.d
    g = np.zeros((d, d), dtype=complex)
    for p in range(d):
        for q in range(d):
            t = basis.tgt[p, q]
            live = t >= 0
            src = np.nonzero(live)[0]
            g[p, q] = np.sum(np.conj(vec[t[src]]) * basis.sgn[p, q][src] * vec[src])
    return g


def exterior_power(C: np.ndarray, basis: DetBasis) -> np.ndarray:
    """Matrix M with M[I,K] = det(C[I,K]); the n-th exterior power of C."""
    dets = np.array(basis.dets, dtype=int)
    rows = dets[:, None, :, None]          # (m, 1, n, 1)
    cols = dets[None, :, None, :]          # (1, m, 1, n)
    return np.linalg.det(C[rows, cols])    # broadcasts to (m, m, n, n)


def natural_orbitals(gamma: np.ndarray, spin: np.ndarray | None = None,
                     tol: float = 1e-9):
    """Descending occupations and the natural-orbital coefficient matrix C.

    When ``spin`` labels the spin orbitals and gamma is block diagonal in spin
    (the collinear case), each block is diagonalized separately so the natural
    orbitals stay spin pure. That choice matters: closed-shell spectra are
    exactly twofold degenerate, so eigh would otherwise return an arbitrary
    alpha/beta mixture inside every degenerate pair.
    """
    d = gamma.shape[0]
    C = np.zeros((d, d), dtype=complex)
    vals = np.zeros(d)
    if spin is None:
        w, V = np.linalg.eigh(gamma)
        vals, C = w, V
    else:
        off = 0.0
        for s in set(spin.tolist()):
            for t in set(spin.tolist()):
                if s != t:
                    blk = gamma[np.ix_(spin == s, spin == t)]
                    off = max(off, float(np.abs(blk).max(initial=0.0)))
        if off > 1e-8:
            w, V = np.linalg.eigh(gamma)
            vals, C = w, V
        else:
            for s in sorted(set(spin.tolist())):
                idx = np.nonzero(spin == s)[0]
                w, V = np.linalg.eigh(gamma[np.ix_(idx, idx)])
                vals[idx] = w
                C[np.ix_(idx, idx)] = V
    order = np.argsort(-vals, kind="stable")
    return np.real(vals[order]), C[:, order]


def to_natural_basis(vec: np.ndarray, basis: DetBasis, C: np.ndarray) -> np.ndarray:
    """CI coefficients in the basis created by b_k^dag = sum_p C[p,k] a_p^dag."""
    M = exterior_power(C, basis)
    return M.conj().T @ vec


# ---------------------------------------------------------------------------
# GPC layer
# ---------------------------------------------------------------------------

TRIVIAL_ROW = "trivial"
GPC_ROW = "gpc"


def classify_row(coeffs: list[int], rhs: int, d: int) -> str:
    """A row is trivial when it is a Pauli/ordering constraint.

    Pauli: lambda_i <= 1 (a single +1 with rhs 1) or -lambda_i <= 0.
    Ordering: lambda_{i+1} - lambda_i <= 0.
    Everything else is a genuine generalized Pauli constraint.
    """
    nz = [(i, c) for i, c in enumerate(coeffs) if c != 0]
    if len(nz) == 1:
        i, c = nz[0]
        if c == 1 and rhs == 1:
            return TRIVIAL_ROW
        if c == -1 and rhs == 0:
            return TRIVIAL_ROW
    if len(nz) == 2 and rhs == 0:
        (i, ci), (j, cj) = nz
        if j == i + 1 and ci == -1 and cj == 1:
            return TRIVIAL_ROW
    return GPC_ROW


@dataclass
class GpcSystem:
    n: int
    d: int
    rows: list[dict] = field(default_factory=list)   # coeffs, rhs, kind, name

    @classmethod
    def load(cls, n: int, d: int) -> "GpcSystem":
        sysd = gpc_constraints.constraints(n, d)
        rows = []
        for k, eq in enumerate(sysd["equalities"]):
            rows.append({"coeffs": list(eq["coeffs"]), "rhs": int(eq["rhs"]),
                         "kind": "equality", "name": f"eq{k}"})
        for k, iq in enumerate(sysd["inequalities"]):
            rows.append({"coeffs": list(iq["coeffs"]), "rhs": int(iq["rhs"]),
                         "kind": classify_row(iq["coeffs"], iq["rhs"], d),
                         "name": f"iq{k}"})
        # the trivial Pauli rows are not always tabulated; add them explicitly
        have = {(tuple(r["coeffs"]), r["rhs"]) for r in rows}
        for i in range(d):
            c = [0] * d
            c[i] = 1
            if (tuple(c), 1) not in have:
                rows.append({"coeffs": c, "rhs": 1, "kind": TRIVIAL_ROW,
                             "name": f"pauli_up{i}"})
            c = [0] * d
            c[i] = -1
            if (tuple(c), 0) not in have:
                rows.append({"coeffs": c, "rhs": 0, "kind": TRIVIAL_ROW,
                             "name": f"pauli_lo{i}"})
        return cls(n, d, rows)

    def slacks(self, lam: np.ndarray) -> np.ndarray:
        A = np.array([r["coeffs"] for r in self.rows], dtype=float)
        b = np.array([r["rhs"] for r in self.rows], dtype=float)
        return b - A @ lam

    def gpc_rows(self) -> list[int]:
        return [i for i, r in enumerate(self.rows) if r["kind"] == GPC_ROW]

    def trivial_rows(self) -> list[int]:
        return [i for i, r in enumerate(self.rows)
                if r["kind"] in (TRIVIAL_ROW, "equality")]


def row_values(row: dict, basis: DetBasis) -> np.ndarray:
    """D(I) = rhs - sum_{i in I} coeffs_i for every determinant I."""
    c = np.array(row["coeffs"], dtype=float)
    return row["rhs"] - basis.occ @ c


def selection_mask(rows: list[dict], basis: DetBasis) -> np.ndarray:
    """Determinants surviving the joint selection rule of ``rows``."""
    keep = np.ones(basis.m, dtype=bool)
    for r in rows:
        keep &= np.isclose(row_values(r, basis), 0.0)
    return keep


# ---------------------------------------------------------------------------
# variational subspace energetics and rigorous bounds
# ---------------------------------------------------------------------------


def subspace_energy(H: np.ndarray, mask: np.ndarray) -> float:
    idx = np.nonzero(mask)[0]
    if idx.size == 0:
        return float("nan")
    sub = H[np.ix_(idx, idx)]
    return float(np.linalg.eigvalsh(sub)[0])


def certified_energy_bound(w_out: float, e0: float, e_max: float) -> float:
    """Rigorous upper bound on E_S - E_0 given the omitted weight.

    With P the projector onto S, Q = 1 - P, Psi the exact ground state and
    w = <Psi|Q|Psi>, the trial vector PPsi/||PPsi|| gives

        E_S - E_0 <= <PPsi|(H - E_0)|PPsi> / (1 - w).

    Because (H - E_0)Psi = 0 the cross term vanishes exactly, so
    <PPsi|(H-E_0)|PPsi> = <QPsi|(H-E_0)|QPsi> <= w (E_max - E_0). This is a
    proof, not a fit: it needs only w, E_0 and the top of the spectrum.
    """
    if w_out >= 1.0:
        return float("inf")
    return w_out * (e_max - e0) / (1.0 - w_out)


# ---------------------------------------------------------------------------
# spin structure
# ---------------------------------------------------------------------------


def spin_labels(d: int) -> np.ndarray:
    """Spin orbital 2*i+s belongs to spatial orbital i with spin s."""
    return np.array([p % 2 for p in range(d)])


# ---------------------------------------------------------------------------
# baselines
# ---------------------------------------------------------------------------


def mask_topk(weights: np.ndarray, k: int, allowed: np.ndarray) -> np.ndarray:
    """Top-k determinants by ``weights``, restricted to ``allowed``."""
    m = weights.shape[0]
    out = np.zeros(m, dtype=bool)
    cand = np.nonzero(allowed)[0]
    if k >= cand.size:
        out[cand] = True
        return out
    order = cand[np.argsort(-weights[cand], kind="stable")]
    out[order[:k]] = True
    return out


def occupation_weights(lam: np.ndarray, basis: DetBasis) -> np.ndarray:
    """Independent-particle probability of each determinant, log scale.

    P(I) = prod_{i in I} lambda_i * prod_{i not in I} (1 - lambda_i). Ranking by
    this is the standard occupation-only guess at which determinants matter and
    is baseline D.
    """
    eps = 1e-300
    lo = np.log(np.clip(lam, eps, None))
    hi = np.log(np.clip(1.0 - lam, eps, None))
    return basis.occ @ lo + (~basis.occ) @ hi


def excitation_rank(basis: DetBasis, ref: int) -> np.ndarray:
    r = basis.dets[ref]
    return np.array([len(set(D) - set(r)) for D in basis.dets])


def random_mask(k: int, allowed: np.ndarray, must: int, rng) -> np.ndarray:
    out = np.zeros(allowed.shape[0], dtype=bool)
    cand = [i for i in np.nonzero(allowed)[0] if i != must]
    pick = rng.choice(len(cand), size=min(k - 1, len(cand)), replace=False)
    out[[cand[i] for i in pick]] = True
    out[must] = True
    return out


# ---------------------------------------------------------------------------
# the analysis of a single reference state
# ---------------------------------------------------------------------------


def degenerate_blocks(lam: np.ndarray, tol: float) -> list[list[int]]:
    blocks, cur = [], [0]
    for i in range(1, len(lam)):
        if abs(lam[i] - lam[i - 1]) <= tol:
            cur.append(i)
        else:
            blocks.append(cur)
            cur = [i]
    blocks.append(cur)
    return blocks


def tie_permutations(blocks: list[list[int]], cap: int = 4096):
    """Relabelings of the natural orbitals inside (near) degenerate blocks.

    The selection rule is a statement about orbital INDICES in the descending
    order, so ties leave it ambiguous. Closed-shell spectra are exactly twofold
    degenerate, which makes this ambiguity generic rather than exotic.
    """
    total = 1
    for b in blocks:
        total *= math.factorial(len(b))
    if total > cap:
        yield list(range(sum(len(b) for b in blocks)))
        return
    ranges = [list(itertools.permutations(b)) for b in blocks]
    for combo in itertools.product(*ranges):
        perm = [0] * sum(len(b) for b in blocks)
        for blk, pm in zip(blocks, combo):
            for old, new in zip(blk, pm):
                perm[old] = new
        yield perm


def permuted_vector(vec: np.ndarray, basis: DetBasis, perm: list[int]) -> np.ndarray:
    """CI vector after relabeling orbital ``i`` as ``perm[i]`` (sign tracked)."""
    out = np.zeros_like(vec)
    for i, D in enumerate(basis.dets):
        img = [perm[p] for p in D]
        order = np.argsort(img)
        sgn = permutation_sign(order)
        newdet = tuple(sorted(img))
        out[basis.index[newdet]] = sgn * vec[i]
    return out


def permutation_sign(order) -> int:
    order = list(order)
    seen = [False] * len(order)
    sign = 1
    for i in range(len(order)):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = order[j]
            ln += 1
        if ln % 2 == 0:
            sign = -sign
    return sign


@dataclass
class Reference:
    """A reference state and everything needed to analyse it."""
    name: str
    family: str                      # "chemistry" or "model"
    geometry: str
    n: int
    d: int
    ham: SpinOrbitalHamiltonian
    vec: np.ndarray                  # exact state, original spin-orbital basis
    e0: float
    e_max: float
    ms2: int | None = None           # 2*Ms, None for spinless models
    s_val: float | None = None       # total spin quantum number
    spinful: bool = True
    orbsym: np.ndarray | None = None  # irrep id per active spatial orbital
    state_irrep: int | None = None
    degeneracy: int = 1
    notes: str = ""


def stage_metrics(H, w, mask, e0, spin_ok, basis, no_spin, s_val):
    k = int(mask.sum())
    out = {"dim": k, "w_out": 1.0 - float(w[mask].sum())}
    out["dim_spin"] = (spin_adapted_dim_general(basis, no_spin, mask, s_val)
                       if spin_ok else k)
    if k:
        e = subspace_energy(H, mask)
        out["e"] = e
        out["de"] = e - e0
    else:
        out["e"] = float("nan")
        out["de"] = float("nan")
    return out


def analyse(ref: Reference, deg_tol: float = 1e-6, tau_list=(0.0,),
            rng=None, n_random: int = 40, sat_tol: float = 1e-8) -> dict:
    """Experiment 1/2/4/5 record for one reference state.

    The reductions are measured as a NESTED CHAIN, because the question is not
    how small a GPC subspace is but how much smaller it is than what particle
    number, Sz, S^2, point group and saturated occupations already give:

        A full C(d,N)  ->  B +Sz  ->  C +point group  ->  D +Pauli freezing
                       ->  E +GPC equalities  ->  F +GPC inequalities

    Only the D -> F step is attributable to nontrivial GPC geometry, and only
    the E -> F step needs quasipinning at all.
    """
    b = DetBasis(ref.d, ref.n)
    spin = spin_labels(ref.d) if ref.spinful else None
    gamma = one_rdm(ref.vec, b)
    lam, C = natural_orbitals(gamma, spin)
    vno = to_natural_basis(ref.vec, b, C)
    gno = one_rdm(vno, b)
    diag_err = float(np.abs(gno - np.diag(np.diag(gno))).max())
    lam_err = float(np.abs(np.real(np.diag(gno)) - lam).max())
    H = ref.ham.rotate(C).matrix(b)
    w = np.abs(vno) ** 2

    sysg = GpcSystem.load(ref.n, ref.d)
    slack = sysg.slacks(lam)
    gpc_idx = sysg.gpc_rows()
    eq_rows = [r for r in sysg.rows if r["kind"] == "equality"]
    pauli_idx = [i for i, r in enumerate(sysg.rows) if r["kind"] == TRIVIAL_ROW]

    # ---- stage B: Sz ------------------------------------------------------
    no_spin = None
    sec = np.ones(b.m, dtype=bool)
    spin_ok = False
    if ref.spinful:
        no_spin = np.array([int(np.argmax(np.abs(C[:, k]) ** 2)) % 2
                            for k in range(ref.d)])
        sz = np.where(no_spin == 0, 0.5, -0.5)
        sec = np.zeros(b.m, dtype=bool)
        sec[b.sector(sz, ref.ms2 / 2.0)] = True
        spin_ok = ref.s_val is not None and \
            int((no_spin == 0).sum()) == int((no_spin == 1).sum())
    idx = np.nonzero(sec)[0]
    ev = np.linalg.eigvalsh(H[np.ix_(idx, idx)])
    e0, e_max = float(ev[0]), float(ev[-1])

    # ---- stage C: point group --------------------------------------------
    pg_mask = sec.copy()
    pg_ok = False
    state_irrep = None
    if ref.orbsym is not None:
        irrep = _no_irreps(C, ref.orbsym, ref.spinful)
        if irrep is not None:
            prod = np.zeros(b.m, dtype=int)
            for k in range(ref.d):
                prod ^= np.where(b.occ[:, k], irrep[k], 0)
            # read the state's irrep off its own leading determinant rather than
            # assuming the totally symmetric one, then check the whole state
            # really sits in that sector before crediting the reduction
            state_irrep = int(prod[int(np.argmax(w))])
            cand = sec & (prod == state_irrep)
            if 1.0 - float(w[cand].sum()) < 1e-9:
                pg_ok = True
                pg_mask = cand

    blocks = degenerate_blocks(lam, deg_tol)
    perms = list(tie_permutations(blocks))
    pmaps = [permutation_map(b, p) for p in perms]

    def best_mask(rows):
        """Tie-optimal selection mask, intersected with the symmetry sector.

        Relabeling orbital i to index perm[i] sends determinant I to perm(I), so
        the rule read in the original labels is the mask pulled back along perm.
        Minimizing the omitted weight over ties is deliberately generous: a user
        without the exact state cannot know which tie-break to take.
        """
        if not rows:
            return pg_mask.copy(), list(range(ref.d)), 1.0 - float(w[pg_mask].sum())
        raw = selection_mask(rows, b)
        best = None
        for perm, pmap in zip(perms, pmaps):
            cand = pg_mask & raw[pmap]
            wout = 1.0 - float(w[cand].sum())
            if best is None or wout < best[2]:
                best = (cand, perm, wout)
        return best

    pauli_sat = [sysg.rows[i] for i in pauli_idx if slack[i] <= sat_tol]
    mask_d, _, _ = best_mask(pauli_sat)
    mask_e, _, _ = best_mask(pauli_sat + eq_rows)

    rec = {
        "name": ref.name, "family": ref.family, "geometry": ref.geometry,
        "n": ref.n, "d": ref.d, "notes": ref.notes,
        "e0_sector": e0, "e_max_sector": e_max, "degeneracy": ref.degeneracy,
        "lambda": [float(x) for x in lam],
        "diag_residual": diag_err, "lambda_residual": lam_err,
        "min_spectrum_gap": float(np.min(np.diff(-lam))),
        "min_spatial_occupation_gap": (
            float(np.min(np.diff(-2.0 * lam[0::2]))) if ref.spinful and ref.ms2 == 0
            and ref.d > 2 else None),
        "n_degenerate_blocks": len(blocks),
        "max_degenerate_block": max(len(x) for x in blocks),
        "pg_available": pg_ok, "state_irrep": state_irrep,
        "spin_adapted_available": spin_ok,
        "n_gpc_rows": len(gpc_idx), "n_gpc_equalities": len(eq_rows),
        "n_pauli_saturated": len(pauli_sat),
        "stage_A_full": {"dim": b.m},
        "stage_B_sz": stage_metrics(H, w, sec, e0, spin_ok, b, no_spin, ref.s_val),
        "stage_C_pg": stage_metrics(H, w, pg_mask, e0, spin_ok, b, no_spin, ref.s_val),
        "stage_D_pauli": stage_metrics(H, w, mask_d, e0, spin_ok, b, no_spin, ref.s_val),
        "stage_E_gpc_eq": stage_metrics(H, w, mask_e, e0, spin_ok, b, no_spin, ref.s_val),
        "stage_D2_freeze": {},
        "stage_F_gpc": [],
    }
    # baseline D in its hard form: freeze orbitals whose occupation is within
    # theta of saturation. This is what an active-space chemist actually does.
    for theta in (1e-3, 1e-2, 5e-2):
        keep = np.ones(b.m, dtype=bool)
        for i in range(ref.d):
            if lam[i] > 1.0 - theta:
                keep &= b.occ[:, i]
            elif lam[i] < theta:
                keep &= ~b.occ[:, i]
        mk = pg_mask & keep
        rec["stage_D2_freeze"][f"{theta:g}"] = stage_metrics(
            H, w, mk, e0, spin_ok, b, no_spin, ref.s_val)

    order = sorted(gpc_idx, key=lambda i: slack[i])
    rec["gpc_slacks_sorted"] = [float(slack[i]) for i in order[:8]]
    rec["nearest_gpc_slack"] = float(slack[order[0]]) if order else None
    rec["nearest_gpc_row"] = sysg.rows[order[0]]["coeffs"] if order else None
    rec["min_pauli_slack"] = float(min(slack[i] for i in pauli_idx)) if pauli_idx else None

    ref_det = int(np.argmax(w))
    exrank = excitation_rank(b, ref_det)
    occw = occupation_weights(lam, b)
    rng = rng or np.random.default_rng(0)

    # ORACLE FACET: pick the single row that minimizes the omitted weight using
    # the exact state. No user can do this, so it upper-bounds what any facet
    # selection rule could achieve and makes the negative result robust to the
    # choice of "nearest".
    if order:
        best = None
        for i in gpc_idx:
            mk, _, wo = best_mask(pauli_sat + eq_rows + [sysg.rows[i]])
            if int(mk.sum()) == 0:
                continue
            if best is None or wo < best[0]:
                best = (wo, i, mk)
        if best is not None:
            wo, i, mk = best
            oracle = stage_metrics(H, w, mk, e0, spin_ok, b, no_spin, ref.s_val)
            oracle["slack"] = float(slack[i])
            oracle["slack_rank"] = int(order.index(i))
            oracle["row"] = sysg.rows[i]["coeffs"]
            rec["oracle_facet"] = oracle

    for tau in tau_list:
        if not order:
            continue
        used = [sysg.rows[i] for i in gpc_idx if slack[i] <= slack[order[0]] + tau]
        mask, perm, _ = best_mask(pauli_sat + eq_rows + used)
        entry = stage_metrics(H, w, mask, e0, spin_ok, b, no_spin, ref.s_val)
        k = entry["dim"]
        entry.update({"tau": tau, "n_rows_used": len(used),
                      "max_slack_used": float(max(slack[i] for i in gpc_idx
                                                  if slack[i] <= slack[order[0]] + tau)),
                      "tie_perm": list(perm)})
        entry["certified_de_from_measured_w"] = certified_energy_bound(
            max(entry["w_out"], 0.0), e0, e_max)
        # can the omitted weight be bounded from the slack at all?
        dvals = np.array([row_values(r, b) for r in used])
        outside = pg_mask & ~mask
        entry["min_row_value_outside"] = (
            float(dvals[:, outside].min()) if outside.any() and len(used) else None)
        entry["one_sided_bound_available"] = bool(
            outside.any() and len(used) and dvals[:, outside].min() >= -1e-9)
        if k:
            m_top = mask_topk(w, k, pg_mask)
            entry["w_out_topk"] = 1.0 - float(w[m_top].sum())
            entry["de_topk"] = subspace_energy(H, m_top) - e0
            m_occ = mask_topk(occw, k, pg_mask)
            entry["w_out_occ"] = 1.0 - float(w[m_occ].sum())
            entry["de_occ"] = subspace_energy(H, m_occ) - e0
            m_ex = mask_topk(-exrank.astype(float), k, pg_mask)
            entry["w_out_exrank"] = 1.0 - float(w[m_ex].sum())
            entry["de_exrank"] = subspace_energy(H, m_ex) - e0
            rw, rd = [], []
            for _ in range(n_random):
                mr = random_mask(k, pg_mask, ref_det, rng)
                rw.append(1.0 - float(w[mr].sum()))
                rd.append(subspace_energy(H, mr) - e0)
            entry["w_out_random_median"] = float(np.median(rw))
            entry["de_random_median"] = float(np.median(rd))
            entry["w_out_random_best"] = float(np.min(rw))
            entry["de_random_best"] = float(np.min(rd))
        rec["stage_F_gpc"].append(entry)
    return rec


def permutation_map(basis: DetBasis, perm) -> np.ndarray:
    """Index array with ``pmap[i]`` the position of ``perm(dets[i])``.

    Pulling a selection mask back along the relabeling is then ``mask[pmap]``.
    """
    return np.array([basis.index[tuple(sorted(perm[p] for p in D))]
                     for D in basis.dets], dtype=np.int64)


_S2_CACHE: dict = {}


def _s2_eig(basis: DetBasis, no_spin_key: tuple):
    """Eigen-decomposition of S^2, cached per (basis, spin labelling).

    Alpha and beta natural orbitals are paired by position inside their spin
    block: for a spin eigenstate the two blocks of gamma share spatial
    eigenvectors, so equal-rank partners are the same spatial orbital.
    """
    key = (basis.d, basis.n, no_spin_key)
    if key in _S2_CACHE:
        return _S2_CACHE[key]
    no_spin = np.array(no_spin_key)
    a = np.nonzero(no_spin == 0)[0]
    bb = np.nonzero(no_spin == 1)[0]
    sz = np.where(no_spin == 0, 0.5, -0.5)
    szv = basis.occ @ sz
    S2 = np.diag(szv * szv + szv)
    for i in range(a.size):
        for j in range(a.size):
            S2 += basis.epq(bb[i], a[i]) @ basis.epq(a[j], bb[j])
    out = np.linalg.eigh(S2)
    _S2_CACHE[key] = out
    return out


def spin_adapted_dim_general(basis: DetBasis, no_spin: np.ndarray,
                             mask: np.ndarray, s_val: float,
                             tol: float = 1e-7) -> int:
    """Spin-adapted dimension inside ``mask`` for arbitrary spin labelling.

    The natural orbitals come in alpha/beta partners because the spin blocks of
    gamma are diagonalized separately and, for a spin eigenstate, carry the same
    spatial eigenvectors. Partners are matched by position within their spin
    block, which is the pairing S^2 needs.
    """
    idx = np.nonzero(mask)[0]
    if idx.size == 0:
        return 0
    a = np.nonzero(no_spin == 0)[0]
    bb = np.nonzero(no_spin == 1)[0]
    if a.size != bb.size:
        return int(mask.sum())
    wv, V = _s2_eig(basis, tuple(no_spin.tolist()))
    sel = np.abs(wv - s_val * (s_val + 1.0)) < 1e-6
    if not sel.any():
        return 0
    P = V[:, sel] @ V[:, sel].conj().T
    ev = np.linalg.eigvalsh(P[np.ix_(idx, idx)])
    return int(np.sum(ev > 1 - tol))


def _no_irreps(C: np.ndarray, orbsym: np.ndarray, spinful: bool):
    """Irrep of each natural orbital, or None if any of them mixes irreps."""
    d = C.shape[0]
    so_sym = np.repeat(orbsym, 2) if spinful else np.asarray(orbsym)
    out = np.zeros(d, dtype=int)
    for k in range(d):
        wgt = np.abs(C[:, k]) ** 2
        tot = {}
        for p in range(d):
            tot[so_sym[p]] = tot.get(so_sym[p], 0.0) + wgt[p]
        best = max(tot, key=tot.get)
        if tot[best] < 1 - 1e-6:
            return None
        out[k] = best
    return out


# ---------------------------------------------------------------------------
# reference-state construction
# ---------------------------------------------------------------------------


def ground_state(ham: SpinOrbitalHamiltonian, n: int, ms2: int | None,
                 spinful: bool = True, irrep_filter=None):
    """Exact ground state in the requested Sz sector, plus spectral extremes."""
    b = DetBasis(ham.d, n)
    H = ham.matrix(b)
    ev_all = np.linalg.eigvalsh(H)
    if spinful and ms2 is not None:
        spin = spin_labels(ham.d)
        sz = np.where(spin == 0, 0.5, -0.5)
        idx = b.sector(sz, ms2 / 2.0)
    else:
        idx = np.arange(b.m)
    if irrep_filter is not None:
        idx = np.array([i for i in idx if irrep_filter(b.dets[i])])
    sub = H[np.ix_(idx, idx)]
    ev, V = np.linalg.eigh(sub)
    vec = np.zeros(b.m, dtype=V.dtype)
    vec[idx] = V[:, 0]
    degen = int(np.sum(np.abs(ev - ev[0]) < 1e-9))
    return vec, float(ev[0]), float(ev_all[-1]), degen


def spinorbital_from_spatial(h1, eri, ecore=0.0, label=""):
    """Spin-orbital integrals from restricted spatial MO integrals."""
    ncas = h1.shape[0]
    d = 2 * ncas
    H1 = np.zeros((d, d))
    G2 = np.zeros((d,) * 4)
    for i in range(ncas):
        for j in range(ncas):
            H1[2 * i, 2 * j] = h1[i, j]
            H1[2 * i + 1, 2 * j + 1] = h1[i, j]
    for i in range(ncas):
        for j in range(ncas):
            for k in range(ncas):
                for ll in range(ncas):
                    v = eri[i, j, k, ll]
                    if v == 0:
                        continue
                    for s in (0, 1):
                        for t in (0, 1):
                            G2[2 * i + s, 2 * j + s, 2 * k + t, 2 * ll + t] = v
    return SpinOrbitalHamiltonian(H1, G2, ecore, label)


def cas_reference(name, geometry, atom, basis_set, ncas, nelecas, ms2,
                  spin_mult=None, charge=0, symmetry=True, unit="Angstrom",
                  notes=""):
    """CASCI reference state in an active space, ready for GPC analysis."""
    from pyscf import ao2mo, gto, mcscf, scf

    try:
        mol = gto.M(atom=atom, basis=basis_set, charge=charge,
                    spin=(spin_mult - 1) if spin_mult else 0,
                    symmetry=symmetry, unit=unit, verbose=0)
    except Exception:                                     # noqa: BLE001
        symmetry = False
        mol = gto.M(atom=atom, basis=basis_set, charge=charge,
                    spin=(spin_mult - 1) if spin_mult else 0,
                    symmetry=False, unit=unit, verbose=0)
    mf = scf.RHF(mol) if mol.spin == 0 else scf.ROHF(mol)
    mf.kernel()
    if not mf.converged:
        mf = mf.newton().run()
    cas = mcscf.CASCI(mf, ncas, nelecas)
    h1, ecore = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), ncas)
    ham = spinorbital_from_spatial(h1, eri, float(ecore), name)
    n = nelecas if isinstance(nelecas, int) else sum(nelecas)
    vec, e0, emax, degen = ground_state(ham, n, ms2)
    orbsym = None
    if symmetry:
        try:
            from pyscf import symm
            mo = mf.mo_coeff[:, cas.ncore:cas.ncore + ncas]
            orbsym = np.array(symm.label_orb_symm(mol, mol.irrep_id,
                                                  mol.symm_orb, mo))
        except Exception:
            orbsym = None
    s_val = ms2 / 2.0
    return Reference(name=name, family="chemistry", geometry=geometry,
                     n=n, d=2 * ncas, ham=ham, vec=vec, e0=e0, e_max=emax,
                     ms2=ms2, s_val=s_val, spinful=True, orbsym=orbsym,
                     degeneracy=degen,
                     notes=(notes + (f"; point group {mol.groupname}"
                                     if symmetry else "")).strip("; "))


def tv_ring_reference(name, geometry, n_orb, n_part, t, V, gradient=0.0,
                      notes=""):
    """Spinless fermions on a ring: the repo's own model family.

    Spinless removes the spin degeneracy entirely, so this is where a
    nondegenerate natural spectrum is available and the selection rule has its
    hypothesis satisfied. It reaches odd d, which collinear chemistry cannot.
    """
    from gpc_census.pinned_physics import tv_ring_integrals
    h1, g2 = tv_ring_integrals(n_orb, t, V, gradient=gradient)
    ham = SpinOrbitalHamiltonian(np.real(h1).copy(), np.real(g2).copy(), 0.0, name)
    vec, e0, emax, degen = ground_state(ham, n_part, None, spinful=False)
    return Reference(name=name, family="model", geometry=geometry,
                     n=n_part, d=n_orb, ham=ham, vec=vec, e0=e0, e_max=emax,
                     ms2=None, s_val=None, spinful=False, degeneracy=degen,
                     notes=notes)


def hubbard_reference(name, geometry, n_sites, n_elec, ms2, t, U,
                      ring=True, gradient=0.0, notes=""):
    d = 2 * n_sites
    h1 = np.zeros((d, d))
    g2 = np.zeros((d,) * 4)
    bonds = range(n_sites) if ring else range(n_sites - 1)
    for i in bonds:
        j = (i + 1) % n_sites
        for s in (0, 1):
            h1[2 * i + s, 2 * j + s] += -t
            h1[2 * j + s, 2 * i + s] += -t
    for i in range(n_sites):
        h1[2 * i, 2 * i] += gradient * i
        h1[2 * i + 1, 2 * i + 1] += gradient * i
        g2[2 * i, 2 * i, 2 * i + 1, 2 * i + 1] = U
        g2[2 * i + 1, 2 * i + 1, 2 * i, 2 * i] = U
    ham = SpinOrbitalHamiltonian(h1, g2, 0.0, name)
    vec, e0, emax, degen = ground_state(ham, n_elec, ms2)
    return Reference(name=name, family="model", geometry=geometry,
                     n=n_elec, d=d, ham=ham, vec=vec, e0=e0, e_max=emax,
                     ms2=ms2, s_val=ms2 / 2.0, spinful=True,
                     degeneracy=degen, notes=notes)


def pairing_reference(name, geometry, n_levels, n_pairs, spacing, G, notes=""):
    """Richardson pairing model, 2*n_pairs electrons in n_levels doublets."""
    d = 2 * n_levels
    h1 = np.zeros((d, d))
    g2 = np.zeros((d,) * 4)
    for i in range(n_levels):
        e = spacing * i
        h1[2 * i, 2 * i] = e
        h1[2 * i + 1, 2 * i + 1] = e
    for i in range(n_levels):
        for j in range(n_levels):
            g2[2 * i, 2 * j, 2 * i + 1, 2 * j + 1] += -G
            g2[2 * i + 1, 2 * j + 1, 2 * i, 2 * j] += -G
    ham = SpinOrbitalHamiltonian(h1, g2, 0.0, name)
    vec, e0, emax, degen = ground_state(ham, 2 * n_pairs, 0)
    return Reference(name=name, family="model", geometry=geometry,
                     n=2 * n_pairs, d=d, ham=ham, vec=vec, e0=e0, e_max=emax,
                     ms2=0, s_val=0.0, spinful=True, degeneracy=degen,
                     notes=notes)


def random_two_body_reference(name, n, d, seed, scale=1.0, notes=""):
    """Random number-conserving two-body Hamiltonian, spinless, real symmetric."""
    rng = np.random.default_rng(seed)
    h1 = rng.normal(size=(d, d))
    h1 = 0.5 * (h1 + h1.T)
    a = rng.normal(size=(d, d, d, d)) * scale
    g2 = 0.25 * (a + a.transpose(1, 0, 3, 2) + a.transpose(2, 3, 0, 1)
                 + a.transpose(3, 2, 1, 0))
    ham = SpinOrbitalHamiltonian(h1, g2, 0.0, name)
    vec, e0, emax, degen = ground_state(ham, n, None, spinful=False)
    return Reference(name=name, family="model", geometry=f"seed{seed}",
                     n=n, d=d, ham=ham, vec=vec, e0=e0, e_max=emax,
                     ms2=None, s_val=None, spinful=False, degeneracy=degen,
                     notes=notes)


# ---------------------------------------------------------------------------
# Experiment 0: does pinning imply the selection rule?
# ---------------------------------------------------------------------------


def census_selection_rule_audit() -> dict:
    """Replay the selection rule against every certified extremal state.

    The rule under test is the one the whole product idea rests on: if lambda
    saturates a valid row ``coeffs . lambda <= rhs`` then, in the natural
    orbital basis, every determinant carrying amplitude has
    ``D(I) = rhs - sum_{i in I} coeffs_i = 0``. The census ships exact attainers
    for all 799 vertices, so the rule is checkable rather than assumed.
    """
    path = os.path.join(REPO, "results", "data", "states.jsonl")
    from fractions import Fraction as Fr
    out = {"total": 0, "diagonal_gate": 0, "rule_holds": 0,
           "rule_holds_after_tie_relabel": 0, "degenerate": 0,
           "nondegenerate": 0, "per_system": {}}
    with open(path) as fh:
        rows = [json.loads(ln) for ln in fh]
    for r in rows:
        n, d = (int(x) for x in r["system"].strip("()").split(","))
        lam = [Fr(x, r["denominator"]) for x in r["integer_form"]]
        occ = [0.0] * d
        for det, re_, im_ in r["support"]:
            wgt = re_ * re_ + im_ * im_
            for i in det:
                occ[i] += wgt
        out["total"] += 1
        if not all(abs(occ[i] - float(lam[i])) < 1e-8 for i in range(d)):
            continue
        out["diagonal_gate"] += 1
        if len(set(lam)) == d:
            out["nondegenerate"] += 1
        else:
            out["degenerate"] += 1
        sysd = gpc_constraints.constraints(n, d)
        sat = [(e["coeffs"], e["rhs"]) for e in sysd["equalities"]]
        for iq in sysd["inequalities"]:
            if Fr(iq["rhs"]) - sum(Fr(c) * x for c, x in zip(iq["coeffs"], lam)) == 0:
                sat.append((iq["coeffs"], iq["rhs"]))
        supp = [tuple(x[0]) for x in r["support"]]
        ok = all(rhs - sum(c[i] for i in det) == 0
                 for c, rhs in sat for det in supp)
        key = r["system"]
        st = out["per_system"].setdefault(key, {"n": 0, "holds": 0, "holds_relabel": 0})
        st["n"] += 1
        if ok:
            out["rule_holds"] += 1
            st["holds"] += 1
            out["rule_holds_after_tie_relabel"] += 1
            st["holds_relabel"] += 1
            continue
        if _relabel_rescues(lam, sat, supp, d):
            out["rule_holds_after_tie_relabel"] += 1
            st["holds_relabel"] += 1
    return out


def _relabel_rescues(lam, sat, supp, d, cap: int = 200000) -> bool:
    """Is there a relabeling inside degenerate blocks that restores the rule?"""
    blocks, cur = [], [0]
    for i in range(1, d):
        if lam[i] == lam[i - 1]:
            cur.append(i)
        else:
            blocks.append(cur)
            cur = [i]
    blocks.append(cur)
    total = 1
    for blk in blocks:
        total *= math.factorial(len(blk))
    if total > cap:
        return False
    for combo in itertools.product(*[list(itertools.permutations(b)) for b in blocks]):
        perm = [0] * d
        for blk, pm in zip(blocks, combo):
            for old, new in zip(blk, pm):
                perm[old] = new
        if all(rhs - sum(c[perm[i]] for i in det) == 0
               for c, rhs in sat for det in supp):
            return True
    return False


def _gpc_functional(vec, basis, coeffs, spin=None):
    """(D, K) with D = -sum_i coeffs_i lambda_i and K the one-body operator.

    K = sum_i coeffs_i b_i^dag b_i in the natural-orbital basis of ``vec``.
    Because the natural orbitals are stationary to first order when lambda is
    nondegenerate, ``K`` doubles as the gradient generator: d/dpsi* of
    <psi|K|psi> is K psi with K held fixed.
    """
    gamma = one_rdm(vec, basis)
    lam, C = natural_orbitals(gamma, spin)
    Kmat = C @ np.diag(np.asarray(coeffs, dtype=float)) @ C.conj().T
    K = np.zeros((basis.m, basis.m), dtype=complex)
    for p in range(basis.d):
        for q in range(basis.d):
            if Kmat[p, q] != 0:
                K += Kmat[p, q] * basis.epq(p, q)
    return lam, np.real(K) if not np.iscomplexobj(vec) else K, C


def minimize_slack(n, d, row, seed=0, tries=6, maxiter=400):
    """Push a state onto a facet from random starts, then test the rule.

    Returns the best minimizer with its spectrum, its degeneracy pattern and
    the weight its natural-orbital expansion puts outside the selection-rule
    subspace. This is the nondegenerate test the census cannot supply, because
    every census vertex has a degenerate spectrum.
    """
    from scipy.optimize import minimize
    b = DetBasis(d, n)
    coeffs = np.array(row["coeffs"], dtype=float)
    rhs = float(row["rhs"])
    mask = selection_mask([row], b)
    best = None
    for t in range(tries):
        rng = np.random.default_rng(seed * 1000 + t)
        y0 = rng.normal(size=b.m)

        def fun(y):
            psi = y / np.linalg.norm(y)
            lam, K, _ = _gpc_functional(psi, b, coeffs)
            val = rhs - float(coeffs @ lam)
            g = -2.0 * (K @ psi - (psi @ K @ psi) * psi) / np.linalg.norm(y)
            return val, g

        res = minimize(fun, y0, jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "ftol": 1e-16, "gtol": 1e-12})
        psi = res.x / np.linalg.norm(res.x)
        lam, _, C = _gpc_functional(psi, b, coeffs)
        vno = to_natural_basis(psi, b, C)
        w = np.abs(vno) ** 2
        gap = float(np.min(np.diff(-lam)))
        entry = {"slack": rhs - float(coeffs @ lam),
                 "w_out": 1.0 - float(w[mask].sum()),
                 "min_gap": gap,
                 "lambda": [float(x) for x in lam]}
        if best is None or entry["slack"] < best["slack"]:
            best = entry
    return best


# ---------------------------------------------------------------------------
# Experiment 3: how tightly does quasipinning control the omitted weight?
# ---------------------------------------------------------------------------


def _wout_of(y, basis, row, mask, spin=None):
    psi = y / np.linalg.norm(y)
    gamma = one_rdm(psi, basis)
    lam, C = natural_orbitals(gamma, spin)
    vno = to_natural_basis(psi, basis, C)
    w = np.abs(vno) ** 2
    slack = float(row["rhs"] - np.array(row["coeffs"], dtype=float) @ lam)
    return slack, 1.0 - float(w[mask].sum()), lam


def slack_wout_frontier(n, d, row, w_targets, seed=0, maxiter=150, trials=3):
    """Smallest reachable quasipinning distance D at a prescribed omitted weight.

    A certification claim of the form ``W_out <= C * D`` is exactly the claim
    that ``D / W_out`` is bounded below. Minimizing D at fixed W_out measures
    that infimum directly, so the frontier either supplies the constant or
    refutes its existence.
    """
    from scipy.optimize import minimize
    b = DetBasis(d, n)
    coeffs = np.array(row["coeffs"], dtype=float)
    rhs = float(row["rhs"])
    mask = selection_mask([row], b)
    rng = np.random.default_rng(seed)
    out = []
    for w0 in w_targets:
        best = None
        for _trial in range(trials):
            base = rng.normal(size=b.m)
            base[~mask] = 0.0
            if np.linalg.norm(base) == 0:
                continue
            base /= np.linalg.norm(base)
            off = rng.normal(size=b.m)
            off[mask] = 0.0
            if np.linalg.norm(off) == 0:
                continue
            off /= np.linalg.norm(off)
            y0 = math.sqrt(1 - w0) * base + math.sqrt(w0) * off

            def fun(y):
                psi = y / np.linalg.norm(y)
                lam, K, _ = _gpc_functional(psi, b, coeffs)
                val = rhs - float(coeffs @ lam)
                g = -2.0 * (K @ psi - (psi @ K @ psi) * psi) / np.linalg.norm(y)
                return val, g

            def con(y):
                return _wout_of(y, b, row, mask)[1] - w0

            res = minimize(fun, y0, jac=True, method="SLSQP",
                           constraints=[{"type": "eq", "fun": con}],
                           options={"maxiter": maxiter, "ftol": 1e-14})
            slack, wout, lam = _wout_of(res.x, b, row, mask)
            if abs(wout - w0) > 0.15 * w0 + 1e-6:
                continue
            if best is None or slack < best[0]:
                best = (slack, wout, lam)
        if best is None:
            continue
        slack, wout, lam = best
        out.append({"w_target": w0, "w_out": wout, "slack": max(slack, 0.0),
                    "ratio_D_over_W": (max(slack, 0.0) / wout) if wout > 0 else None,
                    "min_gap": float(np.min(np.diff(-lam)))})
    return out


def perturbation_scan(n, d, row, seed=0, n_dirs=200):
    """(D, W_out) pairs from random perturbations of exactly pinned states."""
    b = DetBasis(d, n)
    mask = selection_mask([row], b)
    rng = np.random.default_rng(seed)
    pts = []
    for _ in range(n_dirs):
        base = rng.normal(size=b.m)
        base[~mask] = 0.0
        nb = np.linalg.norm(base)
        if nb == 0:
            continue
        base /= nb
        off = rng.normal(size=b.m)
        off[mask] = 0.0
        no = np.linalg.norm(off)
        if no == 0:
            continue
        off /= no
        for eps in (1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 2e-1):
            y = math.sqrt(1 - eps) * base + math.sqrt(eps) * off
            slack, wout, lam = _wout_of(y, b, row, mask)
            if wout <= 0:
                continue
            pts.append({"slack": max(slack, 0.0), "w_out": wout,
                        "ratio": max(slack, 0.0) / wout,
                        "min_gap": float(np.min(np.diff(-lam)))})
    return pts


# ---------------------------------------------------------------------------
# the system battery
# ---------------------------------------------------------------------------


def h_chain(nat, r, geom="linear"):
    if geom == "linear":
        return "; ".join(f"H 0 0 {i * r:.6f}" for i in range(nat))
    if geom == "ring":
        import math as _m
        rad = r / (2 * _m.sin(_m.pi / nat))
        return "; ".join(
            f"H {rad * _m.cos(2 * _m.pi * i / nat):.6f} "
            f"{rad * _m.sin(2 * _m.pi * i / nat):.6f} 0" for i in range(nat))
    raise ValueError(geom)


def chemistry_battery(quick: bool = False):
    """Reference states from real molecular Hamiltonians.

    Active spaces are chosen so that (N, d) has a trusted GPC system in the
    repository: d = 2 * n_active_orbitals, so d in {6, 8, 10}. Odd d is
    unreachable with collinear spin and is covered by the spinless models.
    """
    jobs = []
    rs = [0.8, 1.0, 1.4, 2.0, 3.0] if quick else [0.7, 0.9, 1.1, 1.4, 1.8, 2.4, 3.2]

    # (3,6): H3 in a minimal basis is a complete FCI, doublet
    for r in rs:
        jobs.append(dict(name="H3_linear", geometry=f"R={r}", atom=h_chain(3, r),
                         basis_set="sto-3g", ncas=3, nelecas=3, ms2=1,
                         spin_mult=2, notes="full FCI, minimal basis"))
    # (3,6) again: lithium atom, the textbook quasipinning system
    for b_ in (["sto-3g", "cc-pvdz"] if not quick else ["cc-pvdz"]):
        jobs.append(dict(name="Li_atom", geometry=b_, atom="Li 0 0 0",
                         basis_set=b_, ncas=3, nelecas=3, ms2=1, spin_mult=2,
                         notes="Borland-Dennis setting from the literature"))
    # (3,8) and (3,10)
    jobs.append(dict(name="Li_atom_cas34", geometry="cc-pvdz", atom="Li 0 0 0",
                     basis_set="cc-pvdz", ncas=4, nelecas=3, ms2=1, spin_mult=2))
    jobs.append(dict(name="Li_atom_cas35", geometry="cc-pvdz", atom="Li 0 0 0",
                     basis_set="cc-pvdz", ncas=5, nelecas=3, ms2=1, spin_mult=2))
    for r in rs[::2]:
        jobs.append(dict(name="H3_cas34", geometry=f"R={r}", atom=h_chain(3, r),
                         basis_set="6-31g", ncas=4, nelecas=3, ms2=1, spin_mult=2))
        jobs.append(dict(name="H3_cas35", geometry=f"R={r}", atom=h_chain(3, r),
                         basis_set="6-31g", ncas=5, nelecas=3, ms2=1, spin_mult=2))

    # (4,8): the static-correlation workhorses
    for r in rs:
        jobs.append(dict(name="H4_linear", geometry=f"R={r}", atom=h_chain(4, r),
                         basis_set="sto-3g", ncas=4, nelecas=4, ms2=0,
                         notes="full FCI, minimal basis"))
        jobs.append(dict(name="H4_ring", geometry=f"R={r}", atom=h_chain(4, r, "ring"),
                         basis_set="sto-3g", ncas=4, nelecas=4, ms2=0,
                         notes="square H4, degenerate frontier orbitals"))
    for r in ([1.6, 2.4, 3.4] if quick else [1.2, 1.6, 2.0, 2.4, 3.0, 4.0]):
        jobs.append(dict(name="LiH_cas44", geometry=f"R={r}",
                         atom=f"Li 0 0 0; H 0 0 {r}", basis_set="6-31g",
                         ncas=4, nelecas=4, ms2=0))
    for r in ([1.1, 1.6, 2.4] if quick else [1.0, 1.1, 1.3, 1.6, 2.0, 2.6]):
        jobs.append(dict(name="N2_cas44", geometry=f"R={r}",
                         atom=f"N 0 0 0; N 0 0 {r}", basis_set="6-31g",
                         ncas=4, nelecas=4, ms2=0, notes="N2 stretch"))
    for r in ([0.96, 1.6, 2.4] if quick else [0.96, 1.3, 1.6, 2.0, 2.6]):
        jobs.append(dict(name="H2O_cas44", geometry=f"ROH={r}",
                         atom=f"O 0 0 0; H 0 {0.7575 / 0.9584 * r:.4f} "
                              f"{0.5865 / 0.9584 * r:.4f}; "
                              f"H 0 {-0.7575 / 0.9584 * r:.4f} "
                              f"{0.5865 / 0.9584 * r:.4f}",
                         basis_set="6-31g", ncas=4, nelecas=4, ms2=0,
                         notes="symmetric OH stretch"))
    for r in ([1.4, 2.2] if quick else [1.0, 1.334, 1.7, 2.2, 3.0]):
        jobs.append(dict(name="BeH2_cas44", geometry=f"R={r}",
                         atom=f"Be 0 0 0; H 0 0 {r}; H 0 0 -{r}",
                         basis_set="6-31g", ncas=4, nelecas=4, ms2=0))
    # triplet controls: the spin degeneracy is lifted
    for r in rs[::2]:
        jobs.append(dict(name="H4_linear_triplet", geometry=f"R={r}",
                         atom=h_chain(4, r), basis_set="sto-3g", ncas=4,
                         nelecas=4, ms2=2, spin_mult=3,
                         notes="Ms=1 control, spectrum degeneracy lifted"))

    # (4,10) and (6,10)
    for r in ([1.4, 2.4] if quick else [0.9, 1.4, 1.8, 2.4, 3.2]):
        jobs.append(dict(name="H4_cas45", geometry=f"R={r}", atom=h_chain(4, r),
                         basis_set="6-31g", ncas=5, nelecas=4, ms2=0))
    for r in ([1.4, 2.4] if quick else [1.0, 1.4, 1.8, 2.4, 3.2]):
        jobs.append(dict(name="H6_cas65", geometry=f"R={r}", atom=h_chain(6, r),
                         basis_set="sto-3g", ncas=5, nelecas=6, ms2=0,
                         notes="(6,10), reached by particle-hole duality"))
    # (5,10)
    for r in ([1.2, 2.0] if quick else [0.9, 1.2, 1.6, 2.0, 2.8]):
        jobs.append(dict(name="H5_linear", geometry=f"R={r}", atom=h_chain(5, r),
                         basis_set="sto-3g", ncas=5, nelecas=5, ms2=1,
                         spin_mult=2, notes="full FCI, doublet, (5,10)"))
    return jobs


def model_battery(quick: bool = False):
    out = []
    vs = [0.5, 1.0, 2.0, 4.0, 8.0] if not quick else [1.0, 4.0]
    for (nn, dd) in [(3, 6), (3, 7), (3, 8), (3, 9), (3, 10),
                     (4, 8), (4, 9), (4, 10), (5, 10)]:
        for V in vs:
            out.append(("tv_ring", dict(
                name=f"tvring_{nn}_{dd}", geometry=f"V/t={V}", n_orb=dd,
                n_part=nn, t=1.0, V=V, gradient=0.17,
                notes="spinless fermions: nondegenerate spectrum by construction")))
    us = [1.0, 4.0, 8.0, 16.0] if not quick else [4.0]
    for U in us:
        out.append(("hubbard", dict(name="hubbard_4site", geometry=f"U/t={U}",
                                    n_sites=4, n_elec=4, ms2=0, t=1.0, U=U,
                                    ring=True, gradient=0.0)))
        out.append(("hubbard", dict(name="hubbard_3site", geometry=f"U/t={U}",
                                    n_sites=3, n_elec=3, ms2=1, t=1.0, U=U,
                                    ring=False, gradient=0.0)))
        out.append(("hubbard", dict(name="hubbard_5site", geometry=f"U/t={U}",
                                    n_sites=5, n_elec=5, ms2=1, t=1.0, U=U,
                                    ring=False, gradient=0.0)))
    for G in ([0.1, 0.3, 0.6, 1.0] if not quick else [0.3]):
        out.append(("pairing", dict(name="pairing_4x2", geometry=f"G={G}",
                                    n_levels=4, n_pairs=2, spacing=1.0, G=G)))
        out.append(("pairing", dict(name="pairing_5x2", geometry=f"G={G}",
                                    n_levels=5, n_pairs=2, spacing=1.0, G=G)))
    nseed = 4 if quick else 12
    for (nn, dd) in [(3, 6), (3, 7), (3, 8), (4, 8), (3, 9), (4, 9),
                     (3, 10), (4, 10), (5, 10)]:
        for sd in range(nseed):
            for sc in (0.25, 1.0):
                out.append(("random2b", dict(name=f"rand2b_{nn}_{dd}",
                                             n=nn, d=dd, seed=sd, scale=sc,
                                             notes="random number-conserving "
                                                   "two-body Hamiltonian")))
    return out


def build_reference(kind, spec):
    if kind == "tv_ring":
        return tv_ring_reference(**spec)
    if kind == "hubbard":
        return hubbard_reference(**spec)
    if kind == "pairing":
        return pairing_reference(**spec)
    if kind == "random2b":
        spec = dict(spec)
        spec["notes"] = spec.get("notes", "")
        return random_two_body_reference(**spec)
    raise ValueError(kind)


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

TAUS = (0.0, 0.005, 0.02, 0.05)


def run_batteries(quick: bool = False, chem: bool = True, model: bool = True):
    records, failures = [], []
    if chem:
        for spec in chemistry_battery(quick):
            t0 = time.time()
            try:
                ref = cas_reference(**spec)
                rec = analyse(ref, tau_list=TAUS)
                rec["seconds"] = round(time.time() - t0, 2)
                records.append(rec)
            except Exception as exc:                      # noqa: BLE001
                failures.append({"spec": spec, "error": repr(exc)})
    if model:
        for kind, spec in model_battery(quick):
            t0 = time.time()
            try:
                ref = build_reference(kind, spec)
                rec = analyse(ref, tau_list=TAUS)
                rec["seconds"] = round(time.time() - t0, 2)
                records.append(rec)
            except Exception as exc:                      # noqa: BLE001
                failures.append({"spec": {"kind": kind, **spec}, "error": repr(exc)})
    return records, failures


def experiment3(quick: bool = False):
    """(D, W_out) statistics: random perturbations, then a directed hunt."""
    out = {"perturbation": [], "frontier": [], "gap_binned": []}
    targets = [(3, 6), (3, 7), (3, 8), (4, 8), (3, 9), (4, 9)] if not quick \
        else [(3, 6), (3, 7)]
    for (n, d) in targets:
        sysg = GpcSystem.load(n, d)
        rows = sysg.gpc_rows()[: (2 if quick else 3)]
        for ri in rows:
            row = sysg.rows[ri]
            pts = perturbation_scan(n, d, row, seed=ri, n_dirs=40 if quick else 80)
            ratios = np.array([p["ratio"] for p in pts])
            gaps = np.array([p["min_gap"] for p in pts])
            wo = np.array([p["w_out"] for p in pts])
            out["perturbation"].append({
                "system": f"({n},{d})", "row": row["coeffs"], "rhs": row["rhs"],
                "n_points": len(pts),
                "ratio_min": float(ratios.min()), "ratio_median": float(np.median(ratios)),
                "ratio_max": float(ratios.max()),
                "w_out_range": [float(wo.min()), float(wo.max())],
            })
            for lo, hi in [(0, 1e-4), (1e-4, 1e-3), (1e-3, 1e-2), (1e-2, 1e-1), (1e-1, 10)]:
                sel = (gaps >= lo) & (gaps < hi)
                if sel.sum() >= 3:
                    out["gap_binned"].append({
                        "system": f"({n},{d})", "gap_lo": lo, "gap_hi": hi,
                        "n": int(sel.sum()),
                        "ratio_min": float(ratios[sel].min()),
                        "ratio_median": float(np.median(ratios[sel])),
                    })
            if (n, d) in ((3, 9), (4, 9)):
                # perturbation statistics only: the directed frontier is priced
                # at the four smallest systems, where it is affordable
                continue
            fr = slack_wout_frontier(n, d, row,
                                     [1e-3, 1e-2, 5e-2] if quick
                                     else [1e-4, 1e-3, 1e-2, 5e-2, 1.5e-1],
                                     seed=ri)
            for f in fr:
                f["system"] = f"({n},{d})"
                f["row"] = row["coeffs"]
            out["frontier"].extend(fr)
    return out


def experiment0(quick: bool = False):
    out = {"census_audit": census_selection_rule_audit(), "nondegenerate": []}
    targets = [(3, 6), (3, 7), (3, 8), (4, 8)] if not quick else [(3, 6), (3, 7)]
    for (n, d) in targets:
        sysg = GpcSystem.load(n, d)
        for ri in sysg.gpc_rows()[: (2 if quick else 6)]:
            row = sysg.rows[ri]
            best = minimize_slack(n, d, row, seed=ri, tries=3 if quick else 5)
            best["system"] = f"({n},{d})"
            best["row"] = row["coeffs"]
            best.pop("lambda", None)
            out["nondegenerate"].append(best)
    return out


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------


def make_plots(payload: dict, outdir: str) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:                                     # noqa: BLE001
        return []
    os.makedirs(outdir, exist_ok=True)
    recs = [r for r in payload.get("records", []) if r["degeneracy"] == 1]
    made = []

    def f_entry(r, tau=0.0):
        for e in r["stage_F_gpc"]:
            if e["tau"] == tau:
                return e
        return None

    # 1. dimension reduction, stage by stage
    fig, ax = plt.subplots(figsize=(8, 4.5))
    fams = sorted({(r["n"], r["d"]) for r in recs})
    xs, labels = [], []
    for k, nd in enumerate(fams):
        sub = [r for r in recs if (r["n"], r["d"]) == nd]
        vals = []
        for r in sub:
            e = f_entry(r)
            if e is None or e["dim"] == 0:
                continue
            vals.append((r["stage_A_full"]["dim"] / r["stage_B_sz"]["dim"],
                         r["stage_B_sz"]["dim"] / max(e["dim"], 1)))
        if not vals:
            continue
        v = np.array(vals)
        xs.append(k)
        labels.append(f"({nd[0]},{nd[1]})")
        ax.scatter([k - 0.12] * len(v), v[:, 0], s=14, c="tab:blue", alpha=.6)
        ax.scatter([k + 0.12] * len(v), v[:, 1], s=14, c="tab:red", alpha=.6)
    ax.axhline(10, ls="--", c="k", lw=.8)
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_ylabel("reduction factor")
    ax.set_title("blue: full -> Sz sector (ordinary symmetry).  "
                 "red: Sz sector -> GPC subspace.  dashed: 10x target")
    fig.tight_layout()
    path = os.path.join(outdir, "dimension_reduction.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    made.append(path)

    # 2. omitted norm against quasipinning distance
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for fam, col in (("chemistry", "tab:blue"), ("model", "tab:orange")):
        d_, w_ = [], []
        for r in recs:
            if r["family"] != fam:
                continue
            e = f_entry(r)
            if e is None or r["nearest_gpc_slack"] is None:
                continue
            d_.append(max(r["nearest_gpc_slack"], 1e-17))
            w_.append(max(e["w_out"], 1e-17))
        ax.scatter(d_, w_, s=16, alpha=.6, c=col, label=fam)
    fr = payload.get("experiment3_quasipinning", {}).get("frontier", [])
    if fr:
        ax.scatter([max(f["slack"], 1e-17) for f in fr],
                   [f["w_out"] for f in fr], marker="x", s=55, c="k",
                   label="adversarial frontier (min D at fixed W)")
    lim = [1e-16, 3]
    ax.plot(lim, lim, ls="--", c="grey", lw=.9, label="W = D")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-16, 3)
    ax.set_ylim(1e-16, 3)
    ax.set_xlabel("quasipinning distance D (nearest GPC facet)")
    ax.set_ylabel("omitted weight W_out")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(outdir, "omitted_norm_vs_quasipinning.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    made.append(path)

    # 3. energy error against subspace dimension
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for tau, mk in ((0.0, "o"), (0.02, "s")):
        x, y = [], []
        for r in recs:
            e = f_entry(r, tau)
            if e is None or e["dim"] == 0 or not np.isfinite(e["de"]):
                continue
            x.append(e["dim"] / r["stage_B_sz"]["dim"])
            y.append(max(abs(e["de"]) * 1000, 1e-9))
        ax.scatter(x, y, s=15, marker=mk, alpha=.6, label=f"GPC tau={tau}")
    ax.axhline(1.0, ls="--", c="k", lw=.9, label="1 mEh")
    ax.set_yscale("log")
    ax.set_xlabel("GPC dimension / symmetry-sector dimension")
    ax.set_ylabel("|E_S - E_0|  (mEh)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(outdir, "energy_error_vs_dimension.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    made.append(path)

    # 4. GPC against equal-dimension baselines
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    g, tk, oc, rd = [], [], [], []
    for r in recs:
        e = f_entry(r)
        if e is None or e["dim"] == 0 or "de_topk" not in e:
            continue
        g.append(max(abs(e["de"]) * 1000, 1e-9))
        tk.append(max(abs(e["de_topk"]) * 1000, 1e-9))
        oc.append(max(abs(e["de_occ"]) * 1000, 1e-9))
        rd.append(max(abs(e["de_random_median"]) * 1000, 1e-9))
    ax.scatter(g, tk, s=16, alpha=.6, label="top-k CI coefficients")
    ax.scatter(g, oc, s=16, alpha=.6, label="occupation ranking")
    ax.scatter(g, rd, s=16, alpha=.6, label="random subset (median)")
    lim = [1e-9, 1e4]
    ax.plot(lim, lim, ls="--", c="k", lw=.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("GPC subspace error (mEh)")
    ax.set_ylabel("equal-dimension baseline error (mEh)")
    ax.set_title("below the diagonal: the baseline beats GPC")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = os.path.join(outdir, "gpc_vs_baselines.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    made.append(path)
    return made


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="short battery")
    ap.add_argument("--out", default=os.path.join(REPO, "results", "data",
                                                  "chemistry_gpc_value_test.json"))
    ap.add_argument("--only", choices=["all", "e0", "e123", "e3"], default="all")
    ap.add_argument("--plots", default=os.path.join(REPO, "plots"))
    ap.add_argument("--merge", nargs="*", default=[],
                    help="merge previously written part files into the output")
    args = ap.parse_args(argv)

    payload = {"generated_by": "scripts/chemistry_gpc_value_test.py",
               "quick": args.quick}
    t0 = time.time()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    def checkpoint():
        """Write after every section: a late failure must not discard the work."""
        payload["seconds"] = round(time.time() - t0, 1)
        tmp = args.out + ".part"
        with open(tmp, "w") as fh:
            json.dump(payload, fh, indent=1, sort_keys=True)
        os.replace(tmp, args.out)

    for part in args.merge:
        with open(part) as fh:
            payload.update(json.load(fh))
    payload["selection_subspace_scaling"] = selection_subspace_census()
    checkpoint()
    if args.only in ("all", "e0"):
        payload["experiment0_selection_rule"] = experiment0(args.quick)
        checkpoint()
    if args.only in ("all", "e123"):
        recs, fails = run_batteries(args.quick)
        payload["records"] = recs
        payload["failures"] = fails
        for r in recs:
            r["simple_diagnostics"] = simple_diagnostics(np.array(r["lambda"]),
                                                         r["n"])
        checkpoint()
    payload["experiment5_spin_forced_facets"] = experiment5(args.quick)
    checkpoint()
    payload["experiment6_spin_polytope"] = experiment6(args.quick)
    checkpoint()
    if args.only in ("all", "e3"):
        payload["experiment3_quasipinning"] = experiment3(args.quick)
        checkpoint()
    made = make_plots(payload, args.plots)
    print(f"wrote {args.out} in {payload['seconds']}s; plots: {made}")
    return 0




# ---------------------------------------------------------------------------
# scaling: how large is a pinned subspace, a priori, as the rank grows?
# ---------------------------------------------------------------------------


def selection_subspace_census() -> list[dict]:
    """Size of the pinned selection-rule subspace of every known GPC facet.

    This is pure combinatorics on the facet coefficients and needs no state, so
    it settles the scaling question before any Hamiltonian is chosen: if the
    fraction |S| / C(d,N) does not fall as the active space grows, the method
    cannot become more valuable at larger rank.
    """
    out = []
    systems = [(3, 6), (3, 7), (3, 8), (4, 8), (3, 9), (4, 9),
               (3, 10), (4, 10), (5, 10), (5, 8), (6, 10), (7, 10)]
    for (n, d) in systems:
        sysg = GpcSystem.load(n, d)
        b = DetBasis(d, n)
        eq_rows = [r for r in sysg.rows if r["kind"] == "equality"]
        base = selection_mask(eq_rows, b) if eq_rows else np.ones(b.m, dtype=bool)
        fracs, negs = [], []
        for i in sysg.gpc_rows():
            row = sysg.rows[i]
            keep = base & np.isclose(row_values(row, b), 0.0)
            fracs.append(int(keep.sum()) / int(base.sum()))
            dv = row_values(row, b)[base]
            negs.append(float(dv.min()))
        out.append({
            "system": f"({n},{d})", "n": n, "d": d, "full_dim": b.m,
            "dim_after_equalities": int(base.sum()),
            "n_gpc_rows": len(fracs),
            "frac_kept_min": float(np.min(fracs)), "frac_kept_median": float(np.median(fracs)),
            "frac_kept_max": float(np.max(fracs)),
            "rows_with_nonnegative_values": int(sum(1 for x in negs if x >= -1e-9)),
            "most_negative_row_value": float(np.min(negs)),
        })
    return out


def simple_diagnostics(lam: np.ndarray, n: int) -> dict:
    """Occupation-only diagnostics that cost nothing beyond the 1-RDM.

    D_GPC is itself a linear functional of lambda, so it earns its keep only by
    beating these at predicting omitted weight or energy error.
    """
    lam = np.asarray(lam, dtype=float)
    p = np.clip(lam, 1e-15, 1 - 1e-15)
    entropy = float(-np.sum(p * np.log(p) + (1 - p) * np.log(1 - p)))
    return {
        "largest_virtual_occupation": float(lam[n]) if n < lam.size else 0.0,
        "smallest_occupied_occupation": float(lam[n - 1]),
        "occupation_entropy": entropy,
        "n_fractional_0p01": int(np.sum((lam > 0.01) & (lam < 0.99))),
        "one_minus_idempotency": float(np.sum(lam * (1 - lam))),
    }


# ---------------------------------------------------------------------------
# Experiment 5: which facets does the spin sector saturate for free?
# ---------------------------------------------------------------------------


def spin_forced_facets(n: int, d: int, ms2: int | None, n_samples: int = 200,
                       seed: int = 0, tol: float = 1e-9) -> dict:
    """Rows that RANDOM states of the sector already saturate.

    A facet saturated by every state of the fixed-(N, Sz) sector carries no
    information beyond particle number and Sz, so any reduction it licenses
    must be credited to spin, not to polytope geometry. The (3,6) case has a
    three-line proof: with n_alpha = 2 and n_beta = 1 in three spatial orbitals
    the state is a 3x3 bipartite pure state of one alpha HOLE and one beta
    particle, so lambda^alpha = 1 - lambda^beta and the Borland-Dennis slack is
    exactly max(1 - 2 p_1, 0) for the largest Schmidt coefficient p_1.
    """
    b = DetBasis(d, n)
    sysg = GpcSystem.load(n, d)
    spin = spin_labels(d) if ms2 is not None else None
    if ms2 is not None:
        sz = np.where(spin == 0, 0.5, -0.5)
        idx = b.sector(sz, ms2 / 2.0)
    else:
        idx = np.arange(b.m)
    if idx.size == 0:
        return {}
    rng = np.random.default_rng(seed)
    gpc_idx = sysg.gpc_rows()
    hits = np.zeros(len(gpc_idx), dtype=int)
    vals = np.array([row_values(sysg.rows[i], b) for i in gpc_idx])
    n_pinned, n_vacuous, nearest_slacks, n_sat_rows = 0, 0, [], []
    for _ in range(n_samples):
        v = np.zeros(b.m)
        v[idx] = rng.normal(size=idx.size)
        v /= np.linalg.norm(v)
        lam, C = natural_orbitals(one_rdm(v, b), spin)
        sl = sysg.slacks(lam)[gpc_idx]
        hits += (sl <= tol)
        nearest_slacks.append(float(sl.min()))
        # the sector in the state's OWN natural-orbital labelling
        if ms2 is None:
            secmask = np.ones(b.m, dtype=bool)
        else:
            nos = np.array([int(np.argmax(np.abs(C[:, k]) ** 2)) % 2
                            for k in range(d)])
            szn = np.where(nos == 0, 0.5, -0.5)
            secmask = np.zeros(b.m, dtype=bool)
            secmask[b.sector(szn, ms2 / 2.0)] = True
        sat = np.nonzero(sl <= tol)[0]
        n_sat_rows.append(int(sat.size))
        if sat.size:
            n_pinned += 1
            # a rule is VACUOUS when it removes no determinant of the sector
            if all(int((secmask & np.isclose(vals[i], 0.0)).sum())
                   == int(secmask.sum()) for i in sat):
                n_vacuous += 1
    return {
        "system": f"({n},{d})", "ms2": ms2, "sector_dim": int(idx.size),
        "n_gpc_rows": len(gpc_idx), "n_samples": n_samples,
        "always_saturated": int(np.sum(hits == n_samples)),
        "saturated_over_90pct": int(np.sum(hits >= 0.9 * n_samples)),
        "never_saturated": int(np.sum(hits == 0)),
        "always_saturated_rows": [sysg.rows[gpc_idx[i]]["coeffs"]
                                  for i in np.nonzero(hits == n_samples)[0][:6]],
        "frac_samples_exactly_pinned": n_pinned / n_samples,
        "frac_pinned_samples_with_only_vacuous_rules":
            (n_vacuous / n_pinned) if n_pinned else None,
        "median_saturated_rows_per_sample": float(np.median(n_sat_rows)),
        "median_nearest_slack": float(np.median(nearest_slacks)),
    }


def experiment5(quick: bool = False) -> list[dict]:
    out = []
    plan = [(3, 6, [1, 3]), (3, 8, [1, 3]), (4, 8, [0, 2, 4]),
            (3, 10, [1, 3]), (4, 10, [0, 2, 4]), (5, 10, [1, 3, 5]),
            (6, 10, [0, 2, 4]), (5, 8, [1, 3]),
            (3, 6, None), (3, 8, None), (4, 8, None), (5, 10, None)]
    ns = 60 if quick else 200
    for entry in plan:
        n, d, mss = entry
        for ms2 in (mss if isinstance(mss, list) else [mss]):
            if ms2 is not None and (ms2 > n or (n - ms2) % 2 or (n + ms2) // 2 > d // 2):
                continue
            rec = spin_forced_facets(n, d, ms2, n_samples=ns, seed=d * 31 + n)
            if rec:
                out.append(rec)
    return out




# ---------------------------------------------------------------------------
# Experiment 6: would a spin-adapted (GPC-spin) polytope change the answer?
# ---------------------------------------------------------------------------


def csf_basis(basis: DetBasis, s_val: float, ms2: int):
    """Orthonormal basis of the fixed-(S, Ms) subspace, as columns."""
    no_spin = spin_labels(basis.d)
    sz = np.where(no_spin == 0, 0.5, -0.5)
    idx = basis.sector(sz, ms2 / 2.0)
    wv, V = _s2_eig(basis, tuple(no_spin.tolist()))
    S2 = (V * wv) @ V.T
    ev, U = np.linalg.eigh(S2[np.ix_(idx, idx)])
    keep = np.abs(ev - s_val * (s_val + 1.0)) < 1e-8
    B = np.zeros((basis.m, int(keep.sum())))
    B[idx] = U[:, keep]
    return B


def spatial_occupations(vec, basis: DetBasis):
    """Descending spatial natural occupations in [0,2] and the alpha-block basis."""
    gamma = np.real(one_rdm(vec, basis))
    ga = gamma[0::2, 0::2]
    gb = gamma[1::2, 1::2]
    m, W = np.linalg.eigh(ga + gb)
    order = np.argsort(-m)
    return m[order], W[:, order]


def spin_polytope_support(basis: DetBasis, B, c, tries=4, seed=0, maxiter=300):
    """max c . n over pure states of the fixed spin sector, n the spatial NOs."""
    from scipy.optimize import minimize
    c = np.asarray(c, dtype=float)
    best = -np.inf
    for t in range(tries):
        rng = np.random.default_rng(seed * 97 + t)
        y0 = rng.normal(size=B.shape[1])

        def fun(y):
            psi = B @ (y / np.linalg.norm(y))
            occ, W = spatial_occupations(psi, basis)
            blk = W @ np.diag(c) @ W.T
            M = np.zeros((basis.d, basis.d))
            M[0::2, 0::2] = blk
            M[1::2, 1::2] = blk
            K = np.einsum("pq,pqij->ij", M, basis.generators())
            g = B.T @ (2.0 * (K @ psi - (psi @ K @ psi) * psi)) / np.linalg.norm(y)
            return -float(c @ occ), -g

        res = minimize(fun, y0, jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "ftol": 1e-15})
        best = max(best, -float(res.fun))
    return best


def spin_agnostic_support(n_elec, d, c):
    """max c . n over the SPIN-AGNOSTIC prediction, by exact LP.

    The prediction is the set of spatial occupation vectors whose doubled
    spectrum (n_i/2 twice, sorted) satisfies every published GPC row for
    (N, d). If the true fixed-S body is strictly smaller than this, a spin
    polytope carries information the existing census does not.
    """
    from scipy.optimize import linprog
    k = d // 2
    sysg = GpcSystem.load(n_elec, d)
    A, b = [], []
    for r in sysg.rows:
        # lambda_{2i-1} = lambda_{2i} = n_i / 2 under the sorted doubling
        coef = np.array([(r["coeffs"][2 * i] + r["coeffs"][2 * i + 1]) / 2.0
                         for i in range(k)])
        A.append(coef)
        b.append(float(r["rhs"]))
        if r["kind"] == "equality":
            A.append(-coef)
            b.append(-float(r["rhs"]))
    for i in range(k - 1):                       # sorted descending
        row = np.zeros(k)
        row[i + 1], row[i] = 1.0, -1.0
        A.append(row)
        b.append(0.0)
    res = linprog(-np.asarray(c, dtype=float), A_ub=np.array(A), b_ub=np.array(b),
                  A_eq=np.ones((1, k)), b_eq=[float(n_elec)],
                  bounds=[(0.0, 2.0)] * k, method="highs")
    return -float(res.fun) if res.success else None


def experiment6(quick: bool = False) -> dict:
    """Is the fixed-S occupation body strictly smaller than the census predicts?"""
    out = {"systems": []}
    plan = [(4, 8, 0.0, 0), (4, 10, 0.0, 0)] if not quick else [(4, 8, 0.0, 0)]
    n_dirs = 40 if quick else 120
    for (n, d, s_val, ms2) in plan:
        b = DetBasis(d, n)
        B = csf_basis(b, s_val, ms2)
        k = d // 2
        rng = np.random.default_rng(5)
        rows = []
        for j in range(n_dirs):
            c = rng.normal(size=k)
            c /= np.linalg.norm(c)
            h_spin = spin_polytope_support(b, B, c, tries=3, seed=j)
            h_gpc = spin_agnostic_support(n, d, c)
            if h_gpc is None:
                continue
            rows.append({"c": [float(x) for x in c], "h_spin": h_spin,
                         "h_gpc": h_gpc, "gap": h_gpc - h_spin})
        gaps = np.array([r["gap"] for r in rows])
        out["systems"].append({
            "system": f"({n},{d})", "s": s_val, "csf_dim": int(B.shape[1]),
            "n_directions": len(rows),
            "directions_with_strictly_smaller_spin_body":
                int(np.sum(gaps > 1e-6)),
            "gap_median": float(np.median(gaps)), "gap_max": float(gaps.max()),
            "gap_p90": float(np.percentile(gaps, 90)),
            "worst_direction": rows[int(np.argmax(gaps))]["c"] if rows else None,
        })
    return out


if __name__ == "__main__":
    raise SystemExit(main())
