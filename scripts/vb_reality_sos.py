#!/usr/bin/env python3
"""v_B reality over the degenerate-block rotations: the residual after holonomy.

scripts/vb_holonomy.py proves, exactly, that no DIAGONAL U(1)^d orbital-phase
gauge realifies psi_B. What it cannot see is the non-diagonal rotation of
degenerate natural orbitals: the 1-RDM of v_B has eigenvalues 20/23 (mult 1),
14/23 (mult 4), 4/23 (mult 4), so any basis in which the state stays a v_B
extremal state is reached by V in the stabilizer G = U(1) x U(4) x U(4). psi_B
has a real form iff

    mu = max_{V in G} |<psi_B | Lambda^4 V | conj psi_B>|   equals 1.

This script provides two layers, and is explicit about which is rigorous:

  LOWER BOUND (numerical, runnable here): a heavy multi-start maximizer over G
  (33 real parameters via matrix exponentials of the block Lie algebras). It
  returns the best overlap found -- a lower bound on mu. Its role is a stress
  test: if it climbs to 1, the complexity claim is REFUTED; if it plateaus well
  below 1 and is stable across seeds, that is strong (but not conclusive)
  evidence for complexity. This is what the paper's 0.9365/0.9917 figures are.

  UPPER BOUND (rigorous, needs an SDP solver): to PROVE mu < 1 one needs a
  certified upper bound. The objective |<psi, Lambda^4 V conj psi>|^2 is a
  polynomial of degree 8 in the real/imaginary entries of the two U(4) blocks,
  subject to the polynomial unitarity constraints V*V = I. A moment-SOS
  (Lasserre) relaxation at order 2-3 returns a certified upper bound; rounding
  the SDP dual to rationals and re-checking the SOS identity in exact arithmetic
  upgrades it to a proof. That requires cvxpy + an SDP backend (SCS/MOSEK),
  which are not vendored here; run `build_sos()` where they are installed. The
  relaxation is documented in build_sos() so it can be executed elsewhere.

Only the holonomy result (vb_holonomy.py) and the exhaustive signed-real
exclusion (interference8_1.py) are rigorous today; this file's numerical layer
is evidence, labeled as such.

Usage:
  python scripts/vb_reality_sos.py --seeds 200        # refine the lower bound
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "results" / "data" / "states.jsonl"
# v_B natural-orbital eigenspaces: mode 0 (20/23), modes 1-4 (14/23), 5-8 (4/23)
BLOCKS = [(0,), (1, 2, 3, 4), (5, 6, 7, 8)]


def _load_vb():
    import sympy as sp

    for ln in LEDGER.read_text().splitlines():
        if not ln.strip():
            continue
        r = json.loads(ln)
        if r["system"] == "(4,9)" and r["index"] == 65:
            cf = r["closed_form"]
            dets = [tuple(t) for t in cf["support_dets"]]
            amps = [complex(sp.N(sp.sympify(a), 40)) for a in cf["pretty"]]
            return dets, amps
    raise SystemExit("v_B (4,9) idx 65 not found")


def _block_gen(params, blocks, d):
    """Block-diagonal skew-Hermitian generator from real params (one u(k) per block)."""
    import numpy as np

    A = np.zeros((d, d), complex)
    off = 0
    for blk in blocks:
        k = len(blk)
        # k^2 real params -> a k x k skew-Hermitian: diagonal imag + off-diagonal
        H = np.zeros((k, k), complex)
        idx = off
        for i in range(k):
            H[i, i] = 1j * params[idx]
            idx += 1
        for i in range(k):
            for j in range(i + 1, k):
                H[i, j] = params[idx] + 1j * params[idx + 1]
                H[j, i] = -params[idx] + 1j * params[idx + 1]
                idx += 2
        off = idx
        for a, ma in enumerate(blk):
            for b, mb in enumerate(blk):
                A[ma, mb] = H[a, b]
    return A


def _nparams(blocks):
    return sum(len(b) * len(b) for b in blocks)


def overlap(dets, amps, blocks, d):
    """Return a function params -> |<psi, Lambda^4 V conj psi>| for V = expm(gen)."""
    import numpy as np
    from scipy.linalg import expm

    c = np.array(amps)
    S = [set(t) for t in dets]

    def f(params):
        V = expm(_block_gen(params, blocks, d))
        tot = 0j
        for a, Ta in enumerate(dets):
            for b, Tb in enumerate(dets):
                sub = V[np.ix_(list(Ta), list(Tb))]
                tot += np.conjugate(c[a]) * np.linalg.det(sub) * np.conjugate(c[b])
        return abs(tot)

    return f


def lower_bound(seeds=200, verbose=True):
    import numpy as np
    from scipy.optimize import minimize

    dets, amps = _load_vb()
    d = max(m for T in dets for m in T) + 1
    f = overlap(dets, amps, BLOCKS, d)
    n = _nparams(BLOCKS)
    best = 0.0
    rng = np.random.default_rng(0)
    for s in range(seeds):
        x0 = rng.uniform(-np.pi, np.pi, n)
        r = minimize(lambda x: -f(x), x0, method="Nelder-Mead",
                     options={"maxiter": 4000, "xatol": 1e-9, "fatol": 1e-12})
        if -r.fun > best:
            best = -r.fun
            if verbose:
                print(f"  seed {s}: new best overlap {best:.10f}", flush=True)
    print(f"\nreduced-group max overlap (lower bound, {seeds} seeds): {best:.10f}")
    print("  < 1 => consistent with complexity (evidence, not proof)."
          if best < 1 - 1e-6 else
          "  ~ 1 => REAL FORM MAY EXIST; complexity claim at risk. Investigate.")
    return best


def build_sos():
    """Certified upper bound on mu via a moment-SOS relaxation. Requires cvxpy+SCS.

    Formulation (execute where an SDP backend is installed):
      variables: real and imaginary parts of the two 4x4 blocks U1, U2 (32 reals),
        plus the U(1) phase (irrelevant to |.| and can be fixed to 1).
      constraints: U1^* U1 = I_4, U2^* U2 = I_4 (16 real quadratic eqs each).
      objective: maximize  |sum_{a,b} conj(c_a) det(V[Ta,Tb]) conj(c_b)|^2,
        a degree-8 polynomial in the 32 reals.
      relaxation: Lasserre order 2 or 3 over the constraint ideal; the optimal
        value is a certified UPPER bound on mu^2. If it is < 1, psi_B has no real
        form under any orbital rotation and v_B is a genuinely complex vertex.
      certificate: round the SDP dual multipliers to rationals and verify the SOS
        decomposition minus the objective is a nonnegative combination of the
        constraint polynomials in exact arithmetic (sympy), turning the numerical
        bound into a proof.
    """
    raise NotImplementedError(
        "SOS certificate requires cvxpy + an SDP backend (SCS/MOSEK), not vendored "
        "here. See the docstring for the exact relaxation to run where available.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=200)
    a = ap.parse_args()
    lower_bound(seeds=a.seeds)


if __name__ == "__main__":
    main()
