#!/usr/bin/env python3
"""How many points does the v103 fixed-rho fiber have? Exhaustive count.

Random multi-start in the full 28-dimensional amplitude space does not
saturate: 800 starts reach 59 distinct points and keep climbing. That is the
wrong search space. Two facts collapse it:

1. Every solution carries the SAME weight vector (verified to 6.7e-14 across
   all converged solutions, results/data/uniqueness_census.json). So the
   weights are not unknowns; only the phases are.
2. Modulo the 9-dimensional orbital-phase gauge, the 14 phases leave exactly
   5 gauge-invariant degrees of freedom, one per loop.

So the fiber is the zero set of the channel conditions on a 5-TORUS, which is
compact and low-dimensional enough to search exhaustively rather than sample.

Coordinates. Choose vectors g_1..g_5 dual to an integer loop basis k_1..k_5,
so that theta = sum_j t_j g_j has loop holonomies Phi_i(theta) = t_i exactly.
Then t ranges over [0, 2*pi)^5 with period 2*pi in each coordinate, and the
channel conditions become, for each one-hop class (A,B) with oriented pairs
(u_a, v_a) and fermionic signs s_a,

    sum_a s_a sqrt(w_{u_a} w_{v_a}) exp(i(theta_{v_a} - theta_{u_a})) = 0,

which is 2 real equations per channel, 10 in total, in 5 real unknowns.

The polygon reading, which predicts the structure. A channel of size m is a
closed polygon with m fixed side lengths: m = 2 forces a single antiparallel
configuration, m = 3 is a rigid triangle admitting exactly its two
reflections, and m = 4 admits a one-parameter family of shapes. v103 has
sizes (2, 3, 2, 4, 2), so the discrete freedom is the triangle's reflection
and the quadrilateral's family is cut down by the lattice couplings: the five
channels contribute 8 within-channel holonomies constrained to a rank-5
lattice, leaving 3 relations.

Usage:
  python scripts/v103_fiber_count.py --grid 9      # coarse
  python scripts/v103_fiber_count.py --grid 13     # the reported count
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.fiber import one_rdm  # noqa: E402
from gpc_census.states import one_hop_classes  # noqa: E402

DATA = ROOT / "results" / "data"
OUT = DATA / "v103_fiber_count.json"
SYSTEM, INDEX = "(3,10)", 103


def setup():
    rec = next(
        json.loads(ln)
        for ln in (DATA / "states.jsonl").read_text().splitlines()
        if ln.strip() and json.loads(ln)["system"] == SYSTEM
        and json.loads(ln)["index"] == INDEX
    )
    dets = [tuple(int(o) for o in t) for t in rec["closed_form"]["support_dets"]]
    d = len(rec["integer_form"])
    n = len(dets)
    cf = rec["closed_form"]
    weights = np.array([w / cf["den"] for w in cf["weights"]], dtype=float)

    incidence = np.array([[1.0 if o in t else 0.0 for o in range(d)] for t in dets])

    # integer loop basis: primitive vectors spanning the left kernel
    import sympy as sp

    basis = []
    for vec in sp.Matrix(incidence.astype(int).tolist()).T.nullspace():
        vec = vec * sp.lcm([q.q for q in vec])
        vec = vec / sp.gcd([sp.Integer(q) for q in vec])
        basis.append([int(q) for q in vec])
    k = np.array(basis, dtype=float)  # (5, 14)

    # dual vectors: g with k @ g.T = I, taken in the row space of k
    g = np.linalg.pinv(k)  # (14, 5), satisfies k @ g = I
    assert np.allclose(k @ g, np.eye(len(basis)), atol=1e-9)

    # oriented channel pairs with fermionic signs
    channels = []
    for (a_mode, b_mode), pairs in sorted(one_hop_classes(dets).items()):
        oriented = []
        for i, j in pairs:
            u, v = (i, j) if a_mode in dets[i] else (j, i)
            sign = ((-1) ** dets[v].index(b_mode)) * ((-1) ** dets[u].index(a_mode))
            oriented.append((u, v, sign))
        channels.append(((a_mode, b_mode), oriented))
    return rec, dets, d, n, weights, k, g, channels


def make_residual(weights, g, channels):
    def residual(t):
        theta = g @ t
        out = []
        for _, oriented in channels:
            total = 0j
            for u, v, s in oriented:
                total += s * np.sqrt(weights[u] * weights[v]) * np.exp(
                    1j * (theta[v] - theta[u])
                )
            out.append(total.real)
            out.append(total.imag)
        return np.array(out)

    return residual


def canonical(t, k, g):
    """Loop holonomies wrapped into [0, 2*pi), the gauge-invariant signature."""
    return np.mod(k @ (g @ t), 2 * np.pi)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--grid", type=int, default=13)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args()

    rec, dets, d, n, weights, k, g, channels = setup()
    residual = make_residual(weights, g, channels)
    nloops = k.shape[0]

    grid = np.linspace(0, 2 * np.pi, a.grid, endpoint=False)
    found = []
    tried = converged = 0
    for start in product(grid, repeat=nloops):
        tried += 1
        sol = least_squares(residual, np.array(start), xtol=1e-14, ftol=1e-14,
                            gtol=1e-14)
        if np.max(np.abs(sol.fun)) > a.tol:
            continue
        converged += 1
        sig = canonical(sol.x, k, g)
        # cluster on the torus, and fold complex conjugation (t -> -t)
        dup = False
        for prev in found:
            for cand in (sig, np.mod(-sig, 2 * np.pi)):
                diff = np.abs(cand - prev)
                diff = np.minimum(diff, 2 * np.pi - diff)
                if np.max(diff) < 1e-5:
                    dup = True
                    break
            if dup:
                break
        if not dup:
            found.append(sig)

    # verify each survivor really reconstructs the 1-RDM
    lam = np.array([v / rec["denominator"] for v in rec["integer_form"]])
    worst = 0.0
    for sig in found:
        t = np.linalg.lstsq(k, sig, rcond=None)[0] if False else None
    for sig in found:
        theta = g @ sig  # k @ g = I, so this realizes holonomies = sig
        c = np.sqrt(weights) * np.exp(1j * theta)
        c = c / np.linalg.norm(c)
        rho = one_rdm(c, dets, d)
        worst = max(worst, float(np.max(np.abs(rho - np.diag(lam)))))

    doc = {
        "generated_by": "scripts/v103_fiber_count.py",
        "vertex": f"{SYSTEM} v{INDEX}",
        "method": (
            "Exhaustive grid over the 5-torus of gauge-invariant loop "
            "holonomies, with the weight vector held fixed at its verified "
            "common value, followed by least-squares refinement and clustering "
            "on the torus with complex conjugation folded in."
        ),
        "why_this_space": (
            "Every solution carries the same weight vector, so only phases "
            "are unknown, and modulo the 9-dimensional orbital-phase gauge "
            "the 14 phases leave exactly 5 degrees of freedom. The search "
            "space is a compact 5-torus, not the 28-dimensional amplitude "
            "space that random multi-start was failing to saturate."
        ),
        "channel_sizes": {str(list(kk)): len(v) for kk, v in channels},
        "grid_per_axis": a.grid,
        "starts": tried,
        "converged": converged,
        "distinct_points": len(found),
        "residual_tolerance": a.tol,
        "worst_reconstruction_error": worst,
        "holonomies": [[float(x) for x in sig] for sig in sorted(
            (list(s) for s in found), key=lambda s: tuple(np.round(s, 6))
        )],
        "evidence": "NUMERICAL",
    }
    print(f"grid {a.grid}^{nloops} = {tried} starts, {converged} converged")
    print(f"distinct points modulo gauge and conjugation: {len(found)}")
    print(f"worst 1-RDM reconstruction error: {worst:.2e}")
    if not a.no_write:
        Path(a.out).write_text(json.dumps(doc, indent=1) + "\n")
        print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
