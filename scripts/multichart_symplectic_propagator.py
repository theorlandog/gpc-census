#!/usr/bin/env python3
"""Multi-chart exact reduced propagator on the 29-determinant transport carrier.

The propagated variables are diagonal pair-incidence weights and a redundant
bank of collision-free 2-RDM coherences.  No wavefunction amplitudes are read
by the reduced RHS.  Several certified spanning trees are maintained; the
best-conditioned tree is selected independently at every RHS evaluation.
Cartesian coherences remain regular when an amplitude crosses zero.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

# Resolve the sibling module from this file's own directory so the script works
# when imported (tests) as well as when run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import transport_extensions as te  # noqa: E402

DEFAULT_OUT = (
    Path(__file__).resolve().parents[1] / "results/data/multichart_symplectic_propagator.json"
)


def support_and_hamiltonian():
    base = [tuple(sorted(a + b)) for a in te.windows(5, 2) for b in te.windows(5, 2, 5)]
    hop = (0, 1, 5, 7)
    adjoint = (1, 0, 7, 5)
    support = te.invariant_closure(base, [hop, adjoint])
    index = {det: i for i, det in enumerate(support)}
    H = np.zeros((len(support), len(support)), dtype=complex)
    for column, determinant in enumerate(support):
        for operator in (hop, adjoint):
            target, sign = te.apply_cross(determinant, *operator)
            if target is not None:
                if target not in index:
                    raise RuntimeError("selected carrier is not invariant")
                H[index[target], column] += sign
    if np.max(np.abs(H - H.conj().T)) > 1e-14:
        raise RuntimeError("Hamiltonian is not Hermitian")
    return support, H


def permutation_sign(sequence):
    inversions = sum(sequence[i] > sequence[j] for i in range(len(sequence)) for j in range(i + 1, len(sequence)))
    return -1 if inversions % 2 else 1


def channels(support, order):
    particles = len(support[0])
    terms = defaultdict(lambda: defaultdict(int))
    for a, left in enumerate(support):
        left_set = set(left)
        for b, right in enumerate(support):
            intersection = sorted(left_set.intersection(right))
            remainder_size = particles - order
            if len(intersection) < remainder_size:
                continue
            for remainder in itertools.combinations(intersection, remainder_size):
                remainder_set = set(remainder)
                lm = tuple(x for x in left if x not in remainder_set)
                rm = tuple(x for x in right if x not in remainder_set)
                terms[(lm, rm)][(a, b)] += permutation_sign(lm + remainder) * permutation_sign(rm + remainder)
    return {key: {edge: value for edge, value in coeff.items() if value} for key, coeff in terms.items() if any(coeff.values())}


def incidence_inverse(support, orbitals=10):
    features = list(itertools.combinations(range(orbitals), 2))
    B = np.array([[int(set(feature).issubset(det)) for det in support] for feature in features], dtype=float)
    rows = []
    rank = 0
    for row in range(len(features)):
        candidate = np.linalg.matrix_rank(B[rows + [row]], tol=1e-11)
        if candidate > rank:
            rows.append(row)
            rank = candidate
        if rank == len(support):
            break
    if rank != len(support):
        raise RuntimeError("pair incidence is not full column rank")
    return features, B, rows, np.linalg.inv(B[rows])


def unique_edge_bank(support):
    edge_to_witness = {}
    adjacency = [set() for _ in support]
    for (lm, rm), coeff in sorted(channels(support, 2).items()):
        if lm >= rm or len(coeff) != 1:
            continue
        (a, b), sign = next(iter(coeff.items()))
        if a == b:
            continue
        edge = tuple(sorted((a, b)))
        if edge not in edge_to_witness:
            if a <= b:
                edge_to_witness[edge] = (a, b, lm, rm, sign)
            else:
                edge_to_witness[edge] = (b, a, rm, lm, sign)
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen = {0}
    queue = deque([0])
    while queue:
        a = queue.popleft()
        for b in adjacency[a]:
            if b not in seen:
                seen.add(b)
                queue.append(b)
    if len(seen) != len(support):
        raise RuntimeError("collision-free graph is disconnected")
    return edge_to_witness, adjacency


def tree_bank(adjacency):
    """Deterministic bank of BFS/DFS spanning trees with different roots/orders."""
    n = len(adjacency)
    trees = []
    signatures = set()
    for root in range(n):
        for reverse in (False, True):
            for mode in ("bfs", "dfs"):
                parent = {root: -1}
                frontier = deque([root])
                while frontier:
                    vertex = frontier.popleft() if mode == "bfs" else frontier.pop()
                    neighbors = sorted(adjacency[vertex], reverse=reverse)
                    for neighbor in neighbors:
                        if neighbor in parent:
                            continue
                        parent[neighbor] = vertex
                        frontier.append(neighbor)
                if len(parent) != n:
                    continue
                edges = tuple(sorted(tuple(sorted((v, p))) for v, p in parent.items() if p >= 0))
                if edges not in signatures:
                    signatures.add(edges)
                    trees.append({"root": root, "parent": parent, "edges": edges})
    return trees


@dataclass
class Atlas:
    support: list
    H: np.ndarray
    pair_features: list
    B: np.ndarray
    pivot_rows: list
    B_inverse: np.ndarray
    edge_bank: dict
    adjacency: list
    trees: list
    edge_list: list
    edge_index: dict

    @classmethod
    def build(cls):
        support, H = support_and_hamiltonian()
        features, B, rows, inverse = incidence_inverse(support)
        bank, adjacency = unique_edge_bank(support)
        trees = tree_bank(adjacency)
        edge_list = sorted(bank)
        return cls(support, H, features, B, rows, inverse, bank, adjacency, trees, edge_list, {e: i for i, e in enumerate(edge_list)})

    def pack_from_X(self, X):
        weights = np.real(np.diag(X))
        z = np.array([X[a, b] for a, b in self.edge_list], dtype=complex)
        return np.concatenate([weights, z.real, z.imag])

    def unpack(self, y):
        m = len(self.support)
        e = len(self.edge_list)
        weights = np.array(y[:m], dtype=float)
        z = np.array(y[m:m + e]) + 1j * np.array(y[m + e:m + 2 * e])
        return weights, z

    def chart_score(self, tree, weights, floor=1e-15):
        # Every non-root parent is a denominator. Penalize small internal weights.
        children = defaultdict(list)
        for vertex, parent in tree["parent"].items():
            if parent >= 0:
                children[parent].append(vertex)
        denominators = [weights[v] for v in children if v != tree["root"] and children[v]]
        root_weight = weights[tree["root"]]
        minimum = min([root_weight] + denominators) if denominators else root_weight
        logarithmic = sum(-math.log10(max(abs(x), floor)) for x in [root_weight] + denominators)
        return (minimum, -logarithmic)

    def dynamic_tree(self, weights):
        """Maximum-conditioning spanning tree; small-weight vertices become leaves."""
        n = len(weights)
        root = int(np.argmax(weights))
        parent = {root: -1}
        visited = {root}
        while len(visited) < n:
            candidates = []
            for a in visited:
                for b in self.adjacency[a]:
                    if b in visited:
                        continue
                    # Prefer a large already-known parent and attach low-weight children late.
                    candidates.append((float(weights[a]), -float(weights[b]), -a, -b, a, b))
            if not candidates:
                raise RuntimeError("coherence graph disconnected during chart construction")
            *_, a, b = max(candidates)
            parent[b] = a
            visited.add(b)
        edges = tuple(sorted(tuple(sorted((v, p))) for v, p in parent.items() if p >= 0))
        return {"root": root, "parent": parent, "edges": edges, "dynamic": True}

    def choose_chart(self, weights, zero_tolerance=2e-13):
        valid = []
        for i, tree in enumerate(self.trees):
            children = defaultdict(list)
            for vertex, parent in tree["parent"].items():
                if parent >= 0:
                    children[parent].append(vertex)
            denominators = [tree["root"]] + [v for v in children if v != tree["root"] and children[v]]
            if all(weights[v] > zero_tolerance for v in denominators):
                valid.append(i)
        if valid:
            return max(valid, key=lambda i: self.chart_score(self.trees[i], weights))
        return None

    def reconstruct_X(self, y, chart_index=None, zero_tolerance=2e-13):
        weights, z = self.unpack(y)
        weights[np.abs(weights) < 5e-15] = 0.0
        if chart_index is None:
            chart_index = self.choose_chart(weights, zero_tolerance)
        tree = self.trees[chart_index] if chart_index is not None else self.dynamic_tree(weights)
        root = tree["root"]
        m = len(weights)
        X = np.zeros((m, m), dtype=complex)
        np.fill_diagonal(X, weights)
        edge_values = {}
        for edge, value in zip(self.edge_list, z):
            edge_values[edge] = value
            a, b = edge
            X[a, b] = value
            X[b, a] = np.conj(value)

        # Propagate root coherences through the chosen tree. Zero-weight leaves are
        # regular: their entire row is zero at the event and their selected edge
        # derivative carries the continuation information.
        children = defaultdict(list)
        for vertex, parent in tree["parent"].items():
            if parent >= 0:
                children[parent].append(vertex)
        root_coherence = {root: complex(weights[root])}
        queue = deque([root])
        while queue:
            parent = queue.popleft()
            for child in children[parent]:
                queue.append(child)
                edge = tuple(sorted((parent, child)))
                edge_value = edge_values[edge]
                if parent == root:
                    value = edge_value if edge[0] == root else np.conj(edge_value)
                elif abs(weights[parent]) > zero_tolerance:
                    oriented = edge_value if edge == (parent, child) else np.conj(edge_value)
                    value = root_coherence[parent] * oriented / weights[parent]
                elif abs(weights[child]) <= zero_tolerance:
                    value = 0.0j
                else:
                    raise FloatingPointError("chosen chart has a zero internal denominator")
                root_coherence[child] = value
                X[root, child] = value
                X[child, root] = np.conj(value)

        if abs(weights[root]) <= zero_tolerance:
            raise FloatingPointError("chosen chart root is singular")
        for a in range(m):
            for b in range(a + 1, m):
                if a == root or b == root:
                    continue
                if weights[a] <= zero_tolerance or weights[b] <= zero_tolerance:
                    value = 0.0j
                else:
                    value = np.conj(root_coherence[a]) * root_coherence[b] / weights[root]
                X[a, b] = value
                X[b, a] = np.conj(value)
        return X, chart_index if chart_index is not None else -1


    def reconstruct_by_synchronization(self, y, zero_tolerance=2e-13):
        """Boundary chart using all certified Cartesian coherences.

        Vertices with zero weight are omitted; their rows and columns are zero.
        The active induced coherence graph synchronizes relative phases without
        any wavefunction amplitudes or division by a vanishing weight.
        """
        weights, z = self.unpack(y)
        weights = np.maximum(weights, 0.0)
        active = {i for i, w in enumerate(weights) if w > zero_tolerance}
        X = np.zeros((len(weights), len(weights)), dtype=complex)
        np.fill_diagonal(X, weights)
        if not active:
            return X, -2
        edge_values = {edge: value for edge, value in zip(self.edge_list, z)}
        root = max(active, key=lambda i: weights[i])
        phase_factor = {root: 1.0 + 0.0j}
        queue = deque([root])
        while queue:
            a = queue.popleft()
            for b in self.adjacency[a]:
                if b not in active or b in phase_factor:
                    continue
                edge = tuple(sorted((a, b)))
                value = edge_values[edge] if edge == (a, b) else np.conj(edge_values[edge])
                scale = math.sqrt(weights[a] * weights[b])
                if scale <= zero_tolerance:
                    continue
                ratio = value / scale
                if abs(ratio) == 0:
                    continue
                phase_factor[b] = phase_factor[a] * ratio / abs(ratio)
                queue.append(b)
        if set(phase_factor) != active:
            missing = sorted(active - set(phase_factor))
            raise FloatingPointError(f"active coherence graph disconnected at boundary: {missing}")
        for a in active:
            for b in active:
                X[a, b] = math.sqrt(weights[a] * weights[b]) * np.conj(phase_factor[a]) * phase_factor[b]
        return X, -2

    def rhs(self, _time, y):
        try:
            X, _ = self.reconstruct_X(y)
        except FloatingPointError:
            X, _ = self.reconstruct_by_synchronization(y)
        # X_ab = conjugate(c_a)c_b; H is real symmetric here.
        derivative = 1j * (self.H @ X - X @ self.H)
        dw = np.real(np.diag(derivative))
        dz = np.array([derivative[a, b] for a, b in self.edge_list])
        return np.concatenate([dw, dz.real, dz.imag])


def rdm_from_X(support, X, order, orbitals=10):
    features = list(itertools.combinations(range(orbitals), order))
    fi = {f: i for i, f in enumerate(features)}
    result = np.zeros((len(features), len(features)), dtype=complex)
    for (lm, rm), coeff in channels(support, order).items():
        result[fi[lm], fi[rm]] = sum(value * X[a, b] for (a, b), value in coeff.items())
    return result


def projector_symplectic(X, dX, dY):
    return float(np.real(1j * np.trace(X @ (dX @ dY - dY @ dX))))


def reduced_symplectic_check(atlas, c0, t_final, trials=3, epsilon=2e-6):
    """Finite-difference symplectic test of the reduced propagator itself."""
    rng = np.random.default_rng(20260803)

    def propagate(c):
        X = np.outer(np.conj(c), c)
        y = atlas.pack_from_X(X)
        sol = solve_ivp(atlas.rhs, (0.0, t_final), y, rtol=4e-11, atol=4e-13,
                        method="DOP853", max_step=0.04)
        if not sol.success:
            raise RuntimeError(sol.message)
        try:
            return atlas.reconstruct_X(sol.y[:, -1])[0]
        except FloatingPointError:
            return atlas.reconstruct_by_synchronization(sol.y[:, -1])[0]

    X0 = np.outer(np.conj(c0), c0)
    X1 = propagate(c0)
    errors = []
    for _ in range(trials):
        u = rng.normal(size=len(c0)) + 1j * rng.normal(size=len(c0))
        v = rng.normal(size=len(c0)) + 1j * rng.normal(size=len(c0))
        u -= c0 * np.vdot(c0, u)
        v -= c0 * np.vdot(c0, v)
        u /= np.linalg.norm(u)
        v /= np.linalg.norm(v)
        cp_u = c0 + epsilon * u; cp_u /= np.linalg.norm(cp_u)
        cm_u = c0 - epsilon * u; cm_u /= np.linalg.norm(cm_u)
        cp_v = c0 + epsilon * v; cp_v /= np.linalg.norm(cp_v)
        cm_v = c0 - epsilon * v; cm_v /= np.linalg.norm(cm_v)
        dX0_u = (np.outer(np.conj(cp_u), cp_u) - np.outer(np.conj(cm_u), cm_u)) / (2 * epsilon)
        dX0_v = (np.outer(np.conj(cp_v), cp_v) - np.outer(np.conj(cm_v), cm_v)) / (2 * epsilon)
        dX1_u = (propagate(cp_u) - propagate(cm_u)) / (2 * epsilon)
        dX1_v = (propagate(cp_v) - propagate(cm_v)) / (2 * epsilon)
        omega0 = projector_symplectic(X0, dX0_u, dX0_v)
        omega1 = projector_symplectic(X1, dX1_u, dX1_v)
        errors.append(abs(omega1 - omega0))
    return max(errors)


def symplectic_random_pair_check(H, c0, t, trials=8):
    """Baseline: the exact unitary propagator preserves the projective form.

    This is a consistency floor rather than a test of the reduced flow. For a
    unitary U the two sides are the same inner product, so what it actually
    measures is the unitarity of expm at this step size. The reduced flow is
    tested by reduced_symplectic_check below.
    """
    U = expm(-1j * H * t)
    rng = np.random.default_rng(20260802)
    errors = []
    for _ in range(trials):
        u = rng.normal(size=len(c0)) + 1j * rng.normal(size=len(c0))
        v = rng.normal(size=len(c0)) + 1j * rng.normal(size=len(c0))
        # Horizontal projective tangents.
        u -= c0 * np.vdot(c0, u)
        v -= c0 * np.vdot(c0, v)
        u1, v1 = U @ u, U @ v
        omega0 = 2.0 * np.imag(np.vdot(u, v))
        omega1 = 2.0 * np.imag(np.vdot(u1, v1))
        errors.append(abs(omega1 - omega0))
    return max(errors)


def build_initial_state(H):
    # Choose a full-support state designed to force a true zero crossing in one
    # two-level hopping block: equal real amplitudes on a coupled pair evolve as
    # cos(t)-i sin(t), which does not vanish.  Instead use an eigen-combination
    # with a relative phase so one component crosses zero near pi/4.
    m = H.shape[0]
    c = (1.0 + 0.12 * np.cos(np.arange(m))) * np.exp(0.11j * np.arange(m))
    pairs = np.argwhere(np.triu(np.abs(H) > 0, 1))
    # Use a coupled pair whose vanishing endpoint is not an articulation point
    # of the collision-free coherence graph, so an overlapping boundary chart exists.
    a, b = map(int, pairs[1])
    sign = float(np.real(H[a, b]))
    c[a] = 1.0
    c[b] = -1.0j * sign
    c = c / np.linalg.norm(c)
    return c, (a, b)


def run(output):
    atlas = Atlas.build()
    c0, designed_pair = build_initial_state(atlas.H)
    X0 = np.outer(np.conj(c0), c0)
    y0 = atlas.pack_from_X(X0)
    times = np.array(sorted(set(np.linspace(0.0, 3.0, 121).tolist() + [math.pi / 4, 3 * math.pi / 4])))
    solution = solve_ivp(atlas.rhs, (times[0], times[-1]), y0, t_eval=times, rtol=2e-12, atol=2e-14, method="DOP853", max_step=0.025)
    if not solution.success:
        raise RuntimeError(solution.message)

    max_X = max_G2 = max_G3 = 0.0
    min_weight = 1.0
    switch_events = []
    zero_events = []
    samples = []
    previous_chart = None
    for k, t in enumerate(times):
        y = solution.y[:, k]
        weights, _ = atlas.unpack(y)
        chart = atlas.choose_chart(weights)
        try:
            X_reduced, used_chart = atlas.reconstruct_X(y, chart)
        except FloatingPointError:
            X_reduced, used_chart = atlas.reconstruct_by_synchronization(y)
        chart = used_chart
        U = expm(-1j * atlas.H * t)
        c_direct = U @ c0
        X_direct = np.outer(np.conj(c_direct), c_direct)
        G2_direct = rdm_from_X(atlas.support, X_direct, 2)
        G2_reduced = rdm_from_X(atlas.support, X_reduced, 2)
        G3_direct = rdm_from_X(atlas.support, X_direct, 3)
        G3_reduced = rdm_from_X(atlas.support, X_reduced, 3)
        x_error = float(np.max(np.abs(X_direct - X_reduced)))
        g2_error = float(np.max(np.abs(G2_direct - G2_reduced)))
        g3_error = float(np.max(np.abs(G3_direct - G3_reduced)))
        max_X = max(max_X, x_error)
        max_G2 = max(max_G2, g2_error)
        max_G3 = max(max_G3, g3_error)
        current_min = float(np.min(np.abs(c_direct) ** 2))
        min_weight = min(min_weight, current_min)
        if previous_chart is not None and chart != previous_chart:
            switch_events.append({"t": float(t), "from": int(previous_chart), "to": int(chart), "minimum_weight": current_min})
        previous_chart = chart
        near_zero = np.where(np.abs(c_direct) ** 2 < 2e-5)[0]
        if len(near_zero):
            zero_events.append({"t": float(t), "vertices": [int(i) for i in near_zero], "minimum_weight": current_min, "chart": int(chart)})
        if k % 5 == 0 or len(near_zero):
            samples.append({"t": float(t), "chart": int(chart), "minimum_weight": current_min, "max_X_error": x_error, "max_Gamma2_error": g2_error, "max_Gamma3_error": g3_error})

    # Check that the reduced initial data are genuinely only 2-RDM observables.
    G2_0 = rdm_from_X(atlas.support, X0, 2)
    pair_diag = np.real(np.diag(G2_0))
    weights_from_G2 = atlas.B_inverse @ pair_diag[atlas.pivot_rows]
    edge_from_G2 = []
    feature_index = {f: i for i, f in enumerate(atlas.pair_features)}
    for edge in atlas.edge_list:
        a, b, lm, rm, sign = atlas.edge_bank[edge]
        edge_from_G2.append(G2_0[feature_index[lm], feature_index[rm]] / sign)
    observable_initial_error = max(float(np.max(np.abs(weights_from_G2 - np.abs(c0) ** 2))), float(np.max(np.abs(np.array(edge_from_G2) - np.array([X0[a, b] for a, b in atlas.edge_list])))))

    symplectic_error = symplectic_random_pair_check(atlas.H, c0, times[-1])
    reduced_symplectic_error = reduced_symplectic_check(atlas, c0, 1.1)
    result = {
        "schema_version": 1,
        "method": "multi-chart Cartesian collision-free 2-RDM propagator",
        "carrier_dimension": len(atlas.support),
        "ambient_fci_dimension": math.comb(10, 4),
        "collision_free_edge_bank_size": len(atlas.edge_list),
        "certified_spanning_trees": len(atlas.trees),
        "reduced_real_variables": len(y0),
        "physical_projective_dimension": 2 * len(atlas.support) - 2,
        "integration": {"method": "DOP853", "rtol": 2e-12, "atol": 2e-14, "samples": len(times), "t_final": float(times[-1])},
        "designed_coupled_pair": list(designed_pair),
        "observable_initialization_max_error": observable_initial_error,
        "chart_switch_count": len(switch_events),
        "chart_switches": switch_events,
        "near_zero_samples": zero_events,
        "minimum_direct_weight": min_weight,
        "max_X_error": max_X,
        "max_Gamma2_error": max_G2,
        "max_Gamma3_error": max_G3,
        "direct_unitary_projective_symplectic_form_max_error": symplectic_error,
        "reduced_flow_projective_symplectic_form_max_error": reduced_symplectic_error,
        "wavefunction_use_boundary": "Wavefunction amplitudes are used only to initialize a test case and for independent direct-validation comparisons. The reduced RHS reads only 2-RDM-derived weights and collision-free coherences.",
        "samples": samples,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    result = run(args.out)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
