"""Exact certificates for census verdicts.

Design verdicts are certified by their witness. Interference is a disjunctive
system because every one-hop conflict requires at least one endpoint to have
zero weight. A global certificate therefore combines exact Farkas hitting
clauses with a finite branch proof over one-hop-free supports.

SciPy is used only to discover Farkas vectors. Every accepted vector and every
step of the exported branch proof are checked in exact rational arithmetic.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial, gcd, lcm


def _system(n: int, d: int, spectrum):
    dets = list(combinations(range(d), n))
    spectrum = [Fraction(x) for x in spectrum]
    # rows of A: per-mode sums then normalization; b: spectrum then 1
    a = [[Fraction(1) if m in t else Fraction(0) for t in dets] for m in range(d)]
    a.append([Fraction(1)] * len(dets))
    b = spectrum + [Fraction(1)]
    return dets, a, b


def _conflict_adjacency(dets, n: int) -> list[int]:
    adjacency = [0] * len(dets)
    for i, left in enumerate(dets):
        left_set = set(left)
        for j in range(i + 1, len(dets)):
            if len(left_set & set(dets[j])) == n - 1:
                adjacency[i] |= 1 << j
                adjacency[j] |= 1 << i
    return adjacency


def _equal_value_blocks(spectrum) -> list[list[int]]:
    grouped: dict[Fraction, list[int]] = defaultdict(list)
    order: list[Fraction] = []
    for i, value in enumerate(map(Fraction, spectrum)):
        if value not in grouped:
            order.append(value)
        grouped[value].append(i)
    return [grouped[value] for value in order]


def _mode_symmetries(spectrum, max_order: int = 100_000) -> list[tuple[int, ...]]:
    """Enumerate permutations inside equal-spectrum blocks."""
    blocks = _equal_value_blocks(spectrum)
    order = 1
    for block in blocks:
        order *= factorial(len(block))
    if order > max_order:
        raise ValueError(f"spectrum symmetry group has order {order}, above {max_order}")

    d = len(spectrum)
    out = []
    choices = [tuple(permutations(block)) for block in blocks]
    for images in product(*choices):
        perm = list(range(d))
        for block, image in zip(blocks, images):
            for source, target in zip(block, image):
                perm[source] = target
        out.append(tuple(perm))
    return out


def _determinant_maps(dets, mode_permutations):
    index = {det: i for i, det in enumerate(dets)}
    return [
        tuple(index[tuple(sorted(perm[m] for m in det))] for det in dets)
        for perm in mode_permutations
    ]


def _permute_mask(mask: int, mapping) -> int:
    out = 0
    while mask:
        bit = mask & -mask
        out |= 1 << mapping[bit.bit_length() - 1]
        mask -= bit
    return out


def _canonical_mask(mask: int, determinant_maps) -> int:
    return min(_permute_mask(mask, mapping) for mapping in determinant_maps)


def _banned_mask(selected: int, adjacency) -> int:
    banned = selected
    scan = selected
    while scan:
        bit = scan & -scan
        banned |= adjacency[bit.bit_length() - 1]
        scan -= bit
    return banned


def _mask_indices(mask: int) -> list[int]:
    out = []
    while mask:
        bit = mask & -mask
        out.append(bit.bit_length() - 1)
        mask -= bit
    return out


def _find_hitting_support(clauses, adjacency, determinant_maps) -> int | None:
    """Find an independent set hitting every clause, or prove none exists."""
    ordered = sorted(clauses, key=lambda mask: (mask.bit_count(), mask))
    failed: set[int] = set()

    def search(selected: int) -> int | None:
        selected = _canonical_mask(selected, determinant_maps)
        if selected in failed:
            return None
        banned = _banned_mask(selected, adjacency)
        best = None
        best_size = 0
        for clause in ordered:
            if clause & selected:
                continue
            available = clause & ~banned
            size = available.bit_count()
            if size == 0:
                failed.add(selected)
                return None
            if best is None or size < best_size:
                best = clause
                best_size = size
        if best is None:
            return selected

        available = best & ~banned
        while available:
            bit = available & -available
            available -= bit
            found = search(selected | bit)
            if found is not None:
                return found
        failed.add(selected)
        return None

    return search(0)


def _primitive_integer_vector(values) -> list[int]:
    values = [Fraction(value) for value in values]
    scale = 1
    for value in values:
        scale = lcm(scale, value.denominator)
    integers = [value.numerator * (scale // value.denominator) for value in values]
    divisor = 0
    for value in integers:
        divisor = gcd(divisor, abs(value))
    if divisor:
        integers = [value // divisor for value in integers]
    return integers


def _positive_mask(vector, a) -> int:
    mask = 0
    for j in range(len(a[0])):
        if sum(Fraction(vector[i]) * a[i][j] for i in range(len(a))) > 0:
            mask |= 1 << j
    return mask


def _search_farkas(a, b, keep, max_den: int = 10**6):
    import numpy as np
    from scipy.optimize import linprog

    rows = len(a)
    af = np.array([[float(row[j]) for j in keep] for row in a])
    bf = np.array([float(value) for value in b])
    constraints = af.T if keep else np.zeros((0, rows))
    result = linprog(
        -bf,
        A_ub=constraints,
        b_ub=np.zeros(len(keep)),
        bounds=[(-1, 1)] * rows,
        method="highs",
    )
    if not result.success or -result.fun <= 1e-9:
        return None

    denominators = (10**3, 10**4, 10**5, max_den)
    for denominator in denominators:
        vector = [
            Fraction(float(value)).limit_denominator(denominator)
            for value in result.x
        ]
        if sum(y * target for y, target in zip(vector, b)) <= 0:
            continue
        if all(
            sum(vector[i] * a[i][j] for i in range(rows)) <= 0
            for j in keep
        ):
            return vector
    return None


def verify_design(n: int, d: int, spectrum, witness) -> bool:
    """Exactly verify a weighted-design witness (integer or rational)."""
    dets, a, b = _system(n, d, spectrum)
    w = [Fraction(x).limit_denominator(10**9) if not isinstance(x, int) else Fraction(x)
         for x in witness]
    total = sum(w)
    if total == 0:
        return False
    w = [x / total for x in w]
    if any(x < 0 for x in w):
        return False
    for row, target in zip(a, b):
        if sum(r * x for r, x in zip(row, w)) != target:
            return False
    return True


def farkas_interference(n: int, d: int, spectrum, support=None, max_den: int = 10**6):
    """Search (float) and verify (exact) a Farkas certificate of infeasibility.

    Returns the rational certificate vector on success, None on failure.
    """
    dets, a, b = _system(n, d, spectrum)
    if support is not None:
        keep = [j for j, t in enumerate(dets) if t in set(map(tuple, support))]
    else:
        keep = list(range(len(dets)))
    return _search_farkas(a, b, keep, max_den=max_den)


def _build_branch_proof(
    clauses,
    clause_sources,
    adjacency,
    determinant_maps,
):
    ordered = sorted(clauses, key=lambda mask: (mask.bit_count(), mask))
    nodes: list[dict | None] = []
    memo: dict[int, int] = {}
    max_depth = 0

    def build(selected: int, depth: int) -> int:
        nonlocal max_depth
        selected = _canonical_mask(selected, determinant_maps)
        if selected in memo:
            return memo[selected]
        max_depth = max(max_depth, depth)
        node_id = len(nodes)
        memo[selected] = node_id
        nodes.append(None)

        banned = _banned_mask(selected, adjacency)
        best = None
        best_size = 0
        for clause in ordered:
            if clause & selected:
                continue
            available = clause & ~banned
            size = available.bit_count()
            if size == 0:
                vector_index, permutation = clause_sources[clause]
                nodes[node_id] = {
                    "kind": "dead",
                    "selected": _mask_indices(selected),
                    "clause": {
                        "vector": vector_index,
                        "permutation": list(permutation),
                    },
                }
                return node_id
            if best is None or size < best_size:
                best = clause
                best_size = size

        if best is None:
            raise RuntimeError("certificate clauses admit an independent hitting set")

        vector_index, permutation = clause_sources[best]
        children = []
        available = best & ~banned
        for determinant in _mask_indices(available):
            child_id = build(selected | (1 << determinant), depth + 1)
            children.append({"determinant": determinant, "node": child_id})
        nodes[node_id] = {
            "kind": "branch",
            "selected": _mask_indices(selected),
            "clause": {
                "vector": vector_index,
                "permutation": list(permutation),
            },
            "children": children,
        }
        return node_id

    root = build(0, 0)
    return {
        "root": root,
        "nodes": nodes,
        "node_count": len(nodes),
        "max_depth": max_depth,
    }


def generate_global_interference_certificate(
    n: int,
    d: int,
    spectrum,
    *,
    integer_form=None,
    denominator: int | None = None,
    vertex_index: int | None = None,
    max_iterations: int = 1000,
) -> dict:
    """Generate an exact global no-design certificate.

    Each exact vector ``y`` gives a hitting clause

        P_y = {T : y dot incidence(T) > 0}.

    If ``y dot (spectrum, 1) > 0``, every nonnegative incidence realization
    must use at least one determinant in ``P_y``. The exported branch DAG
    proves that no one-hop-free support can hit all clauses.
    """
    spectrum = [Fraction(value) for value in spectrum]
    dets, a, b = _system(n, d, spectrum)
    adjacency = _conflict_adjacency(dets, n)
    mode_permutations = _mode_symmetries(spectrum)
    determinant_maps = _determinant_maps(dets, mode_permutations)

    clauses: set[int] = set()
    clause_sources: dict[int, tuple[int, tuple[int, ...]]] = {}
    vectors: list[list[int]] = []

    for _ in range(max_iterations):
        support_mask = _find_hitting_support(clauses, adjacency, determinant_maps)
        if support_mask is None:
            break
        keep = _mask_indices(support_mask)
        candidate = _search_farkas(a, b, keep)
        if candidate is None:
            raise RuntimeError("independent hitting support is incidence-feasible")
        vector = _primitive_integer_vector(candidate)
        if sum(Fraction(y) * target for y, target in zip(vector, b)) <= 0:
            raise RuntimeError("discovered Farkas vector has nonpositive margin")
        positive = _positive_mask(vector, a)
        if positive & support_mask:
            raise RuntimeError("discovered Farkas clause does not cut its support")

        vector_index = len(vectors)
        new_clause = False
        for permutation, mapping in zip(mode_permutations, determinant_maps):
            clause = _permute_mask(positive, mapping)
            if clause not in clauses:
                clauses.add(clause)
                clause_sources[clause] = (vector_index, permutation)
                new_clause = True
        if not new_clause:
            raise RuntimeError("Farkas search repeated an existing clause orbit")
        vectors.append(vector)
    else:
        raise RuntimeError("global certificate generation exceeded its iteration limit")

    proof = _build_branch_proof(
        clauses,
        clause_sources,
        adjacency,
        determinant_maps,
    )
    record = {
        "system": [n, d],
        "spectrum": [str(value) for value in spectrum],
        "symmetry_blocks": _equal_value_blocks(spectrum),
        "farkas_vectors": vectors,
        "expanded_clause_count": len(clauses),
        "proof": proof,
    }
    if integer_form is not None:
        record["integer_form"] = list(integer_form)
    if denominator is not None:
        record["denominator"] = denominator
    if vertex_index is not None:
        record["vertex_index"] = vertex_index
    return record


def verify_global_interference_certificate(record: dict) -> bool:
    """Verify one exported global certificate using exact arithmetic."""
    try:
        n, d = record["system"]
        spectrum = [Fraction(value) for value in record["spectrum"]]
        if len(spectrum) != d:
            return False
        if record["symmetry_blocks"] != _equal_value_blocks(spectrum):
            return False
        if "integer_form" in record:
            denominator = record["denominator"]
            if spectrum != [Fraction(value, denominator) for value in record["integer_form"]]:
                return False

        dets, a, b = _system(n, d, spectrum)
        adjacency = _conflict_adjacency(dets, n)
        mode_permutations = _mode_symmetries(spectrum)
        determinant_maps = _determinant_maps(dets, mode_permutations)
        permutation_maps = dict(zip(mode_permutations, determinant_maps))

        vectors = record["farkas_vectors"]
        positives = []
        clauses = set()
        for vector in vectors:
            if len(vector) != d + 1:
                return False
            if sum(Fraction(y) * target for y, target in zip(vector, b)) <= 0:
                return False
            positive = _positive_mask(vector, a)
            if positive == 0:
                return False
            positives.append(positive)
            for mapping in determinant_maps:
                clauses.add(_permute_mask(positive, mapping))
        if len(clauses) != record["expanded_clause_count"]:
            return False

        proof = record["proof"]
        nodes = proof["nodes"]
        if proof["node_count"] != len(nodes):
            return False
        root = proof["root"]
        if not 0 <= root < len(nodes):
            return False
        if nodes[root]["selected"]:
            return False

        visiting: set[int] = set()
        visited: set[int] = set()

        def check(node_id: int) -> bool:
            if node_id in visited:
                return True
            if node_id in visiting or not 0 <= node_id < len(nodes):
                return False
            visiting.add(node_id)
            node = nodes[node_id]
            selected_indices = node["selected"]
            if selected_indices != sorted(set(selected_indices)):
                return False
            if any(not 0 <= value < len(dets) for value in selected_indices):
                return False
            selected = sum(1 << value for value in selected_indices)
            banned = _banned_mask(selected, adjacency)
            if any(adjacency[value] & selected for value in selected_indices):
                return False

            source = node["clause"]
            vector_index = source["vector"]
            permutation = tuple(source["permutation"])
            if not 0 <= vector_index < len(positives):
                return False
            mapping = permutation_maps.get(permutation)
            if mapping is None:
                return False
            clause = _permute_mask(positives[vector_index], mapping)
            if clause not in clauses or clause & selected:
                return False
            available = clause & ~banned

            if node["kind"] == "dead":
                if available:
                    return False
            elif node["kind"] == "branch":
                children = node.get("children", [])
                candidates = [child["determinant"] for child in children]
                if candidates != _mask_indices(available):
                    return False
                for child in children:
                    child_id = child["node"]
                    if not 0 <= child_id < len(nodes):
                        return False
                    target = selected | (1 << child["determinant"])
                    child_selected = sum(1 << value for value in nodes[child_id]["selected"])
                    if child_selected.bit_count() != selected.bit_count() + 1:
                        return False
                    if _canonical_mask(child_selected, determinant_maps) != _canonical_mask(
                        target, determinant_maps
                    ):
                        return False
                    if not check(child_id):
                        return False
            else:
                return False

            visiting.remove(node_id)
            visited.add(node_id)
            return True

        if not check(root):
            return False
        return len(visited) == len(nodes)
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False
