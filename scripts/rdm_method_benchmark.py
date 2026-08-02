#!/usr/bin/env python3
"""Four-way propagation benchmark for the 29-state observable transport carrier.

Methods:
  fci       - full 4-fermion wavefunction in wedge^4 C^10 (210 complex amplitudes)
  carrier   - invariant 29-state carrier wavefunction
  dense2rdm - full 45x45 complex 2-RDM with exact carrier closure
  sparse2rdm- multi-chart observable coordinates (weights + certified coherences)

Each worker runs in a fresh process so ru_maxrss is method-local. The benchmark
uses the same initial state, Hamiltonian, time interval, DOP853 tolerances, and
final reference for every method.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.sparse import csr_matrix

# Resolve the sibling module from this file's own directory so the script works
# when imported (tests) as well as when run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multichart_symplectic_propagator as mp  # noqa: E402

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "results/data/rdm_method_benchmark.json"

RTOL = 2e-10
ATOL = 2e-12
MAX_STEP = 0.04
T_FINAL = 1.0


def fci_basis(orbitals=10, particles=4):
    return list(itertools.combinations(range(orbitals), particles))


def build_fci_hamiltonian(atlas):
    basis = fci_basis()
    index = {det: i for i, det in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)), dtype=complex)
    hop = (0, 1, 5, 7)
    adjoint = (1, 0, 7, 5)
    for column, det in enumerate(basis):
        for op in (hop, adjoint):
            target, sign = mp.te.apply_cross(det, *op)
            if target is not None:
                H[index[target], column] += sign
    return basis, csr_matrix(H)


def initial_state(atlas):
    c, _ = mp.build_initial_state(atlas.H)
    # The designed state has a later zero crossing; t=1 includes that challenge.
    return c


def dense_gamma2_linear_map(atlas):
    """Sparse map vec(X) -> vec(Gamma2), row-major conventions."""
    features = list(itertools.combinations(range(10), 2))
    fi = {f: i for i, f in enumerate(features)}
    rows, cols, data = [], [], []
    m = len(atlas.support)
    k = len(features)
    for (lm, rm), coeff in mp.channels(atlas.support, 2).items():
        out = fi[lm] * k + fi[rm]
        for (a, b), value in coeff.items():
            rows.append(out); cols.append(a * m + b); data.append(value)
    return csr_matrix((data, (rows, cols)), shape=(k * k, m * m), dtype=complex)


def extract_sparse_observables_from_gamma2(atlas, G2):
    diag = np.real(np.diag(G2))
    weights = atlas.B_inverse @ diag[atlas.pivot_rows]
    feature_index = {f: i for i, f in enumerate(atlas.pair_features)}
    edges = []
    for edge in atlas.edge_list:
        _a, _b, lm, rm, sign = atlas.edge_bank[edge]
        edges.append(G2[feature_index[lm], feature_index[rm]] / sign)
    z = np.asarray(edges, dtype=complex)
    return np.concatenate([weights, z.real, z.imag])


def pack_complex(z):
    return np.concatenate([z.real.ravel(), z.imag.ravel()])


def unpack_complex(y, shape):
    n = int(np.prod(shape))
    return (y[:n] + 1j * y[n:2*n]).reshape(shape)


def solve(method):
    atlas = mp.Atlas.build()
    c0 = initial_state(atlas)
    U = expm(-1j * atlas.H * T_FINAL)
    c_ref = U @ c0
    X_ref = np.outer(np.conj(c_ref), c_ref)
    G2_ref = mp.rdm_from_X(atlas.support, X_ref, 2)
    start = time.perf_counter()

    if method == 'carrier':
        y0 = pack_complex(c0)
        def rhs(_t, y):
            c = unpack_complex(y, (len(c0),))
            return pack_complex(-1j * atlas.H @ c)
        sol = solve_ivp(rhs, (0, T_FINAL), y0, rtol=RTOL, atol=ATOL,
                        method='DOP853', max_step=MAX_STEP)
        c = unpack_complex(sol.y[:, -1], (len(c0),))
        # Phase-insensitive carrier density error.
        X = np.outer(np.conj(c), c)
        error = float(np.max(np.abs(X - X_ref)))
        variables = len(y0)
        complex_state_entries = len(c0)

    elif method == 'fci':
        basis, H = build_fci_hamiltonian(atlas)
        bi = {det: i for i, det in enumerate(basis)}
        psi0 = np.zeros(len(basis), dtype=complex)
        for i, det in enumerate(atlas.support):
            psi0[bi[det]] = c0[i]
        y0 = pack_complex(psi0)
        def rhs(_t, y):
            psi = unpack_complex(y, (len(psi0),))
            return pack_complex(-1j * (H @ psi))
        sol = solve_ivp(rhs, (0, T_FINAL), y0, rtol=RTOL, atol=ATOL,
                        method='DOP853', max_step=MAX_STEP)
        psi = unpack_complex(sol.y[:, -1], (len(psi0),))
        c = np.array([psi[bi[d]] for d in atlas.support])
        X = np.outer(np.conj(c), c)
        outside = psi.copy()
        outside[[bi[d] for d in atlas.support]] = 0
        error = max(float(np.max(np.abs(X - X_ref))), float(np.linalg.norm(outside)))
        variables = len(y0)
        complex_state_entries = len(psi0)

    elif method == 'sparse2rdm':
        X0 = np.outer(np.conj(c0), c0)
        y0 = atlas.pack_from_X(X0)
        sol = solve_ivp(atlas.rhs, (0, T_FINAL), y0, rtol=RTOL, atol=ATOL,
                        method='DOP853', max_step=MAX_STEP)
        try:
            X, _ = atlas.reconstruct_X(sol.y[:, -1])
        except FloatingPointError:
            X, _ = atlas.reconstruct_by_synchronization(sol.y[:, -1])
        error = float(np.max(np.abs(X - X_ref)))
        variables = len(y0)
        complex_state_entries = 0

    elif method == 'dense2rdm':
        M2 = dense_gamma2_linear_map(atlas)
        X0 = np.outer(np.conj(c0), c0)
        G20 = (M2 @ X0.ravel()).reshape(45, 45)
        y0 = pack_complex(G20)
        def rhs(_t, y):
            G2 = unpack_complex(y, (45, 45))
            obs = extract_sparse_observables_from_gamma2(atlas, G2)
            try:
                X, _ = atlas.reconstruct_X(obs)
            except FloatingPointError:
                X, _ = atlas.reconstruct_by_synchronization(obs)
            Xdot = 1j * (atlas.H @ X - X @ atlas.H)
            G2dot = (M2 @ Xdot.ravel()).reshape(45, 45)
            return pack_complex(G2dot)
        sol = solve_ivp(rhs, (0, T_FINAL), y0, rtol=RTOL, atol=ATOL,
                        method='DOP853', max_step=MAX_STEP)
        G2 = unpack_complex(sol.y[:, -1], (45, 45))
        obs = extract_sparse_observables_from_gamma2(atlas, G2)
        try:
            X, _ = atlas.reconstruct_X(obs)
        except FloatingPointError:
            X, _ = atlas.reconstruct_by_synchronization(obs)
        error = max(float(np.max(np.abs(X - X_ref))), float(np.max(np.abs(G2 - G2_ref))))
        variables = len(y0)
        complex_state_entries = 45 * 45
    else:
        raise ValueError(method)

    elapsed = time.perf_counter() - start
    if not sol.success:
        raise RuntimeError(sol.message)
    rss_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        'method': method,
        'elapsed_seconds': elapsed,
        'peak_rss_mib': rss_kib / 1024.0,
        'real_ode_variables': variables,
        'complex_state_entries': complex_state_entries,
        'nfev': int(sol.nfev),
        'accepted_steps': int(len(sol.t) - 1),
        'max_reference_error': error,
        't_final': T_FINAL,
        'rtol': RTOL,
        'atol': ATOL,
        'max_step': MAX_STEP,
    }


def worker(method):
    result = solve(method)
    print(json.dumps(result, sort_keys=True))


def orchestrate(output, repeats=3):
    methods = ['fci', 'carrier', 'dense2rdm', 'sparse2rdm']
    raw = {m: [] for m in methods}
    script = str(Path(__file__).resolve())
    for method in methods:
        for _ in range(repeats):
            proc = subprocess.run([sys.executable, script, '--worker', method],
                                  check=True, capture_output=True, text=True)
            raw[method].append(json.loads(proc.stdout.strip().splitlines()[-1]))
    summary = []
    for method in methods:
        rows = raw[method]
        summary.append({
            'method': method,
            'median_elapsed_seconds': statistics.median(r['elapsed_seconds'] for r in rows),
            'min_elapsed_seconds': min(r['elapsed_seconds'] for r in rows),
            'median_peak_rss_mib': statistics.median(r['peak_rss_mib'] for r in rows),
            'real_ode_variables': rows[0]['real_ode_variables'],
            'complex_state_entries': rows[0]['complex_state_entries'],
            'median_nfev': statistics.median(r['nfev'] for r in rows),
            'median_accepted_steps': statistics.median(r['accepted_steps'] for r in rows),
            'max_reference_error': max(r['max_reference_error'] for r in rows),
        })
    by = {r['method']: r for r in summary}
    result = {
        'schema_version': 1,
        'benchmark': '29-state invariant transport carrier, identical DOP853 integration',
        'system': {'particles': 4, 'orbitals': 10, 'ambient_fci_dimension': math.comb(10,4),
                   'carrier_dimension': 29, 'pair_space_dimension': math.comb(10,2)},
        'integration': {'t_final': T_FINAL, 'rtol': RTOL, 'atol': ATOL, 'max_step': MAX_STEP,
                        'repeats': repeats, 'fresh_process_per_run': True},
        'summary': summary,
        'ratios': {
            'sparse_vs_fci_runtime': by['fci']['median_elapsed_seconds'] / by['sparse2rdm']['median_elapsed_seconds'],
            'sparse_vs_carrier_runtime': by['carrier']['median_elapsed_seconds'] / by['sparse2rdm']['median_elapsed_seconds'],
            'dense_vs_sparse_runtime': by['dense2rdm']['median_elapsed_seconds'] / by['sparse2rdm']['median_elapsed_seconds'],
            'dense_vs_sparse_variable_count': by['dense2rdm']['real_ode_variables'] / by['sparse2rdm']['real_ode_variables'],
            'fci_vs_sparse_variable_count': by['fci']['real_ode_variables'] / by['sparse2rdm']['real_ode_variables'],
        },
        'raw_runs': raw,
        'interpretation_boundary': [
            'This is a fixed small physical system; it does not establish asymptotic superiority.',
            'Carrier-wavefunction propagation is expected to be fastest because it uses the known invariant basis and minimal coordinates.',
            'The relevant reduced-method comparison is dense 2-RDM versus sparse observable coordinates.',
            'Peak RSS includes Python, NumPy, SciPy, and imported module overhead in each fresh process.'
        ]
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--worker', choices=['fci','carrier','dense2rdm','sparse2rdm'])
    p.add_argument('--out', type=Path, default=DEFAULT_OUT)
    p.add_argument('--repeats', type=int, default=3)
    a = p.parse_args()
    if a.worker:
        worker(a.worker)
    else:
        orchestrate(a.out, a.repeats)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
