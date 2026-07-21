#!/usr/bin/env python3
"""Independently verify a stored hybrid extremal state realizes its vertex.

Two independent checks per state, both required:
  1. A from-scratch fermionic 1-RDM (explicit Jordan-Wigner second quantization,
     exact amplitudes evaluated at high precision) whose eigenvalue multiset must
     equal the vertex integer spectrum. Shares NO code with the solver or
     gpc_census.exactify.
  2. The shipped gauge-invariant exact gate gpc_census.exactify.verify_exact,
     the same certificate the 785 census states carry.

A state that passes both IS an exact closed-form extremal state for its vertex.
Reads the docs/hybrid_cracks/*.jsonl records (system, index, den, spec, dets,
amplitudes as sympy srepr).

Usage:
  python scripts/verify_hybrid_state.py docs/hybrid_cracks/v96.jsonl
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _occ(det, d):
    v = [0] * d
    for m in det:
        v[m] = 1
    return v


def _ann(o, p):
    if o[p] == 0:
        return 0, None
    return (-1) ** sum(o[:p]), o[:p] + [0] + o[p + 1:]


def _cre(o, p):
    if o[p] == 1:
        return 0, None
    return (-1) ** sum(o[:p]), o[:p] + [1] + o[p + 1:]


def independent_spectrum(dets, amps, d, den):
    """Numeric eigenvalues (times den) of the 1-RDM, built from first principles."""
    import numpy as np

    occs = [_occ(dd, d) for dd in dets]
    amap = {tuple(o): a for o, a in zip(occs, amps)}
    rho = np.zeros((d, d), dtype=complex)
    for o, c in zip(occs, amps):
        for p in range(d):
            s1, o1 = _ann(o, p)
            if s1 == 0:
                continue
            for q in range(d):
                s2, o2 = _cre(o1, q)
                if s2 == 0:
                    continue
                bra = amap.get(tuple(o2))
                if bra is not None:
                    rho[q, p] += np.conj(bra) * c * s1 * s2
    herm = float(np.max(np.abs(rho - rho.conj().T)))
    trace = float(np.real(np.trace(rho)))
    ev = sorted((np.linalg.eigvalsh(rho) * den).real, reverse=True)
    return ev, herm, trace


def verify_file(path):
    import numpy as np
    import sympy as sp
    from gpc_census.exactify import verify_exact

    n_ok = 0
    n_tot = 0
    for line in open(path):
        h = json.loads(line)
        n_tot += 1
        N, d = h["system"]
        den = h["den"]
        spec = sorted(h["spec"], reverse=True)
        dets = [tuple(x) for x in h["dets"]]
        amps_sp = [sp.sympify(a) for a in h["amplitudes"]]
        amps_num = [complex(sp.N(a, 40)) for a in amps_sp]
        ev, herm, trace = independent_spectrum(dets, amps_num, d, den)
        indep_ok = np.allclose(ev, spec, atol=1e-6)
        spec_frac = [Fraction(x, den) for x in h["spec"]]
        shipped_ok = verify_exact(N, d, spec_frac, dets, amps_sp)
        ok = indep_ok and shipped_ok
        n_ok += ok
        print(f"({N},{d})v{h['index']}: indep_1RDM={indep_ok} "
              f"shipped_verify_exact={shipped_ok} "
              f"herm={herm:.0e} trace={trace:.3f}")
    print(f"\n{n_ok}/{n_tot} states pass BOTH independent checks")
    return n_ok == n_tot


if __name__ == "__main__":
    raise SystemExit(0 if verify_file(sys.argv[1]) else 1)
