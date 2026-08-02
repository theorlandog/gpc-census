#!/usr/bin/env python3
"""Four-pronged exact reduced-dynamics campaign for cyclic-window carriers.

Provides:
1. exact entangling cross-block density interactions;
2. amplitude-free O(m) coordinate propagation and conditioning bounds;
3. scalable benchmarks against block-sector and full-FCI dimensions;
4. exact linear determinacy closure for unresolved census carriers.

The census scan imports support/channel routines from rdm_observability_census.py
and is intended to run from the repository root.
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import sys
import time
import tracemalloc
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable


def windows(d: int, n: int, offset: int = 0) -> list[tuple[int, ...]]:
    return [tuple(sorted(offset + ((k+j) % d) for j in range(n))) for k in range(d)]


def indicator_vector(d: int, n: int, orbital: int) -> list[int]:
    return [int(orbital in w) for w in windows(d, n)]


def interaction_energy(u: list[int], v: list[int], coupling: Fraction = Fraction(1)) -> list[list[Fraction]]:
    return [[coupling * a * b for b in v] for a in u]


def additive_energy_matrix(energy: list[list[Fraction]]) -> bool:
    """Return true iff E_ij = r_i + c_j, the non-entangling diagonal form."""
    rows, cols = len(energy), len(energy[0])
    return all(
        energy[i][j] + energy[0][0] == energy[i][0] + energy[0][j]
        for i in range(rows) for j in range(cols)
    )


def entanglement_witness(d1: int, n1: int, p: int, d2: int, n2: int, q: int) -> dict:
    """Exact 2x2 minor witness for exp(-it n_p n_q) entanglement.

    For full-support product coefficients, choose rows with u=0,1 and columns
    with v=0,1. The evolved coefficient matrix has a minor proportional to
    exp(-it)-1, hence Schmidt rank >=2 unless t is a multiple of 2*pi.
    """
    u, v = indicator_vector(d1, n1, p), indicator_vector(d2, n2, q)
    i0, i1 = u.index(0), u.index(1)
    j0, j1 = v.index(0), v.index(1)
    E = interaction_energy(u, v)
    return {
        "factors": [[d1, n1], [d2, n2]],
        "orbitals": [p, q],
        "selected_rows": [i0, i1],
        "selected_columns": [j0, j1],
        "energy_minor_cross_difference": str(E[i0][j0] + E[i1][j1] - E[i0][j1] - E[i1][j0]),
        "energy_is_additively_separable": additive_energy_matrix(E),
        "generic_schmidt_rank_lower_bound": 2,
        "exceptional_times": "coupling*t in 2*pi*Z",
    }


@dataclass
class ReducedCoordinates:
    """Pure-state coordinates without amplitudes: weights and rooted coherences."""
    weights: list[complex]
    root_coherences: list[complex]  # X_0j for j=1..m-1

    def reconstruct_x(self, a: int, b: int) -> complex:
        if a == b:
            return self.weights[a]
        x0a = self.weights[0] if a == 0 else self.root_coherences[a-1]
        x0b = self.weights[0] if b == 0 else self.root_coherences[b-1]
        # X_ab = conjugate(c_a)c_b = conjugate(X_0a) X_0b / w_0
        return x0a.conjugate() * x0b / self.weights[0]

    def propagate_diagonal(self, energies: list[float], t: float) -> "ReducedCoordinates":
        root = energies[0]
        return ReducedCoordinates(
            self.weights[:],
            [z * complex(math.cos(-(e-root)*t), math.sin(-(e-root)*t)) for z, e in zip(self.root_coherences, energies[1:])],
        )


def uniform_coordinates(m: int) -> ReducedCoordinates:
    w = 1.0 / m
    return ReducedCoordinates([complex(w)] * m, [complex(w)] * (m-1))


def local_path_condition_bound(path_edges: int, min_weight: float, relative_input_error: float) -> float:
    """First-order relative error bound for a path-product reconstruction.

    A path with L edges has L measured coherences and L-1 denominator weights.
    If each has relative error <= eps and every denominator is >= min_weight,
    the first-order relative condition coefficient is 2L-1. min_weight is
    included as an explicit chart-validity guard rather than hidden.
    """
    if min_weight <= 0:
        return math.inf
    return (2 * path_edges - 1) * relative_input_error


def benchmark(max_factors: int = 8, d: int = 5, n: int = 2) -> list[dict]:
    rows = []
    for q in range(1, max_factors + 1):
        m = d**q
        block_dim = math.comb(d, n)**q
        full_dim = math.comb(d*q, n*q)
        tracemalloc.start()
        t0 = time.perf_counter()
        coords = uniform_coordinates(m)
        # Cross-block n_0 n_0 interaction on the first two factors, represented
        # directly as diagonal carrier energies. Remaining factors are spectators.
        energies = [0.0] * m
        if q >= 2:
            u = indicator_vector(d, n, 0)
            for idx in range(m):
                a = (idx // (d ** (q-1))) % d
                b = (idx // (d ** (q-2))) % d
                energies[idx] = float(u[a] * u[b])
        moved = coords.propagate_diagonal(energies, 0.37)
        checksum = sum(z.real for z in moved.root_coherences[: min(1024, m-1)])
        elapsed = time.perf_counter() - t0
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        rows.append({
            "factors": q,
            "carrier_dimension": m,
            "block_sector_dimension": block_dim,
            "full_fci_dimension": full_dim,
            "full_fci_to_carrier_ratio": full_dim / m,
            "coordinate_count": 2*m-1,
            "elapsed_seconds": elapsed,
            "peak_bytes": peak,
            "checksum": checksum,
        })
    return rows


def rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    if not matrix:
        return [], []
    a = [row[:] for row in matrix]
    rows, cols = len(a), len(a[0])
    pivots, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        s = a[r][c]
        a[r] = [x/s for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                s = a[i][c]
                a[i] = [x-s*y for x,y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return a[:r], pivots


def linearly_identifiable_columns(matrix: list[list[int]]) -> set[int]:
    """Coordinates x_j uniquely determined by y=A x are e_j in row(A)."""
    if not matrix:
        return set()
    reduced, _ = rref([[Fraction(x) for x in row] for row in matrix])
    cols = len(matrix[0])
    identifiable = set()
    for j in range(cols):
        # e_j lies in rowspace iff adding it does not increase rank.
        _, piv0 = rref([row[:] for row in reduced])
        e = [Fraction(0)] * cols
        e[j] = 1
        _, piv1 = rref([row[:] for row in reduced] + [e])
        if len(piv1) == len(piv0):
            identifiable.add(j)
    return identifiable


def load_observability_module(repo: Path):
    path = repo / "scripts" / "rdm_observability_census.py"
    spec = importlib.util.spec_from_file_location("rdm_observability_census", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def determinacy_record(module, support: list[tuple[int, ...]]) -> dict:
    channels = module.signed_channels(support, 2)
    variables = [(a,b) for a in range(len(support)) for b in range(a+1, len(support))]
    index = {edge:i for i,edge in enumerate(variables)}
    rows = []
    for (left_modes, right_modes), terms in channels.items():
        if left_modes >= right_modes:
            continue
        row = [0] * len(variables)
        for (a,b), coeff in terms.items():
            if a == b:
                continue
            if a < b:
                row[index[(a,b)]] += coeff
            else:
                # A reversed coherence is the conjugate variable, which is not
                # a complex-linear function of the ones in `variables`, so the
                # row would be invalid rather than merely weaker. No channel in
                # the shipped census reaches this branch; assert rather than
                # discard the term silently if one ever does.
                raise AssertionError(
                    f"channel {(left_modes, right_modes)} mixes the conjugate "
                    f"coherence {(a, b)}; the linear model does not cover it"
                )
        if any(row):
            rows.append(row)
    identifiable = linearly_identifiable_columns(rows)
    adjacency = {i:set() for i in range(len(support))}
    for j in identifiable:
        a,b = variables[j]
        adjacency[a].add(b); adjacency[b].add(a)
    seen = set()
    components = []
    for root in range(len(support)):
        if root in seen: continue
        stack=[root]; seen.add(root); comp=[]
        while stack:
            v=stack.pop(); comp.append(v)
            for w in adjacency[v]:
                if w not in seen: seen.add(w); stack.append(w)
        components.append(sorted(comp))
    return {
        "support_size": len(support),
        "coherence_variables": len(variables),
        "channel_equations": len(rows),
        "linear_rank": len(rref([[Fraction(x) for x in row] for row in rows])[1]) if rows else 0,
        "linearly_identifiable_coherences": len(identifiable),
        "identifiable_graph_component_sizes": sorted((len(c) for c in components), reverse=True),
        "resolved_by_linear_plus_purity_path_closure": len(components) == 1,
    }


def scan_census(repo: Path) -> dict:
    module = load_observability_module(repo)
    states = [json.loads(line) for line in (repo / "results/data/states.jsonl").read_text().splitlines() if line.strip()]
    baseline = [module.record_certificate(record) for record in states]
    unresolved = [(record, cert) for record, cert in zip(states, baseline) if not cert["pure_full_support_reconstruction_from_gamma2"]]
    records=[]
    for record, cert in unresolved:
        support=[tuple(x) for x in record["closed_form"]["support_dets"]]
        result=determinacy_record(module, support)
        result.update({"system":record["system"], "index":record["index"], "classified":record["classified"]})
        records.append(result)
    return {
        "schema_version":1,
        "method":"exact full-channel linear identifiability followed by pure-state path closure",
        "baseline_unresolved":len(unresolved),
        "newly_resolved":sum(r["resolved_by_linear_plus_purity_path_closure"] for r in records),
        "remaining_unresolved":sum(not r["resolved_by_linear_plus_purity_path_closure"] for r in records),
        "records":records,
        "nonclaim":"Failure of this certificate is not non-injectivity; nonlinear rank-one equations may still resolve additional coherences.",
    }


def main() -> int:
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="cmd",required=True)
    e=sub.add_parser("entangle"); e.add_argument("--out",type=Path)
    b=sub.add_parser("benchmark"); b.add_argument("--max-factors",type=int,default=8); b.add_argument("--out",type=Path)
    c=sub.add_parser("closure"); c.add_argument("--repo",type=Path,default=Path(".")); c.add_argument("--out",type=Path,default=Path("results/data/rdm_determinacy_closure.json"))
    args=p.parse_args()
    if args.cmd=="entangle":
        result=entanglement_witness(5,2,0,5,2,0)
    elif args.cmd=="benchmark":
        result={"schema_version":1,"family":"q copies of S(5,2)","rows":benchmark(args.max_factors)}
    else:
        result=scan_census(args.repo)
    text=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if getattr(args,"out",None): args.out.parent.mkdir(parents=True,exist_ok=True); args.out.write_text(text)
    print(text,end="")
    return 0

if __name__=="__main__": raise SystemExit(main())
