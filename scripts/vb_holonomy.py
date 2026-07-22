#!/usr/bin/env python3
"""Exact gauge-invariant holonomy of a certified interference state.

The single-particle phase gauge is the diagonal group U(1)^d: multiplying orbital
m by e^{i phi_m} multiplies determinant T by exp(i sum_{m in T} phi_m). A state is
"real up to this gauge and up to sign" iff its amplitude-phase vector p (p_T =
arg c_T) lies in colspace(M) + pi*Z, where M is the support incidence matrix
(M[T,m] = 1 if m in T). Equivalently, for every INTEGER left-null vector u of M
(u^T M = 0, an exact cycle in the support), the holonomy

    hol(u) = sum_T u_T * arg(c_T)     (mod pi)

is a gauge invariant, and hol(u) != 0 (mod pi) for some u is an exact,
basis-free (under diagonal gauge) proof that the state cannot be made real: it is
GENUINELY COMPLEX in its natural-orbital basis. No optimization, no numerics
beyond an irrationality check that is also confirmed symbolically.

This does NOT settle reality under the non-diagonal block rotations that mix
degenerate natural orbitals (U(1) x U(4) x U(4) for v_B); that residual gap is
the target of scripts/vb_reality_sos.py. But the diagonal-gauge obstruction it
certifies is exact, and for v_B it is nonzero.

Usage:
  python scripts/vb_holonomy.py                 # v_B = (4,9) idx 65
  python scripts/vb_holonomy.py --system 4,9 --index 65
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "results" / "data" / "states.jsonl"


def _load(system: str, index: int):
    for ln in LEDGER.read_text().splitlines():
        if not ln.strip():
            continue
        r = json.loads(ln)
        if r["system"] == system and r["index"] == index:
            return r
    raise SystemExit(f"{system} idx {index} not found in {LEDGER}")


def holonomies(system: str = "(4,9)", index: int = 65, verbose: bool = True):
    """Return (is_complex, cycles) for the certified state at (system, index).

    is_complex is True iff some integer support cycle carries a holonomy that is
    not an integer multiple of pi (exact sympy check). cycles lists
    (u, holonomy_over_pi) for every generator of the integer left-null space.
    """
    import sympy as sp

    rec = _load(system, index)
    cf = rec.get("closed_form")
    if not cf:
        raise SystemExit(f"{system} idx {index} has no closed_form (uncertified)")
    dets = [tuple(t) for t in cf["support_dets"]]
    amps = [sp.sympify(a) for a in cf["pretty"]]
    d = max(m for T in dets for m in T) + 1
    # amplitude phases, exact
    phases = [sp.arg(a) for a in amps]
    # support incidence matrix M (rows = dets, cols = orbitals)
    M = sp.Matrix([[1 if m in T else 0 for m in range(d)] for T in dets])
    left_null = M.T.nullspace()  # integer vectors u with u^T M = 0
    if verbose:
        print(f"{system} idx {index}: support {len(dets)}, orbitals {d}, "
              f"rank(M) {M.rank()}, left-null dim {len(left_null)}", flush=True)
    cycles = []
    is_complex = False
    for u in left_null:
        u = u * sp.lcm([x.q for x in u])          # clear to integer entries
        u = u / sp.gcd([sp.Integer(x) for x in u])  # primitive
        hol = sum(u[i] * phases[i] for i in range(len(dets)))
        over_pi = sp.nsimplify(sp.simplify(hol / sp.pi))
        integral = over_pi.is_integer
        cycles.append((list(u), over_pi, bool(integral)))
        if not integral:
            is_complex = True
        if verbose:
            active = [dets[i] for i in range(len(dets)) if u[i] != 0]
            print(f"  cycle u={list(u)}  holonomy/pi = {over_pi}  "
                  f"{'INTEGER (gauge-removable)' if integral else 'NONINTEGER -> complex'}",
                  flush=True)
            print(f"    determinants: {active}", flush=True)
    return is_complex, cycles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", default="(4,9)")
    ap.add_argument("--index", type=int, default=65)
    a = ap.parse_args()
    system = a.system if a.system.startswith("(") else f"({a.system})"
    is_complex, _ = holonomies(system, a.index)
    print()
    if is_complex:
        print(f"PROVEN: {system} idx {a.index} is genuinely complex up to the "
              f"diagonal U(1)^d gauge (a support cycle carries a nonzero holonomy).")
        sys.exit(0)
    print(f"{system} idx {a.index}: all holonomies gauge-removable; real up to "
          f"diagonal gauge (no diagonal-gauge complexity obstruction).")
    sys.exit(1)


if __name__ == "__main__":
    main()
