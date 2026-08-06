"""Dependency-free exact rational polytope tools.

Everything here is exact over the rationals: a Bland-rule simplex for
feasibility, a double-description vertex enumerator driven by an algebraic
adjacency test, and facet extraction from a vertex set. No numpy, no lrs.

``polytope.py`` enumerates the ambient census polytopes with lrs; this module
serves the orbit-closure layer (``orbit.py``), which has to enumerate many
small polytopes inside a single process and must not depend on an external
binary being installed.

Conventions. A polytope is carried as ``(ineqs, eqs)`` in ambient dimension
``n``:

- ``ineqs``: list of ``(a, c)`` meaning ``<a, x> >= c``,
- ``eqs``: list of ``(a, c)`` meaning ``<a, x> == c``.

The ``>=`` orientation is used throughout because every inequality this
project produces for an orbit closure is a lower bound on a partial sum of
occupations.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations

Vec = tuple[Fraction, ...]
Row = tuple[list[Fraction], Fraction]


# --------------------------------------------------------------------------
# exact linear algebra
# --------------------------------------------------------------------------

def rref(rows: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    """Reduced row echelon form and the list of pivot columns."""
    mat = [[Fraction(x) for x in r] for r in rows]
    if not mat:
        return [], []
    ncols = len(mat[0])
    pivots: list[int] = []
    r = 0
    for col in range(ncols):
        piv = next((i for i in range(r, len(mat)) if mat[i][col] != 0), None)
        if piv is None:
            continue
        mat[r], mat[piv] = mat[piv], mat[r]
        pv = mat[r][col]
        mat[r] = [x / pv for x in mat[r]]
        for i in range(len(mat)):
            if i != r and mat[i][col] != 0:
                f = mat[i][col]
                mat[i] = [x - f * y for x, y in zip(mat[i], mat[r])]
        pivots.append(col)
        r += 1
        if r == len(mat):
            break
    return mat[:r], pivots


def rank(rows: list[list[Fraction]]) -> int:
    """Exact rank of a rational matrix."""
    return len(rref(rows)[1])


def solve_square(rows: list[list[Fraction]], rhs: list[Fraction]) -> list[Fraction] | None:
    """Unique solution of a (possibly overdetermined) consistent system, else None."""
    aug = [list(r) + [b] for r, b in zip(rows, rhs)]
    red, piv = rref(aug)
    n = len(rows[0]) if rows else 0
    if n in piv:  # a pivot in the augmented column: inconsistent
        return None
    if len(piv) < n:  # underdetermined: not a unique point
        return None
    out = [Fraction(0)] * n
    for i, c in enumerate(piv):
        out[c] = red[i][n]
    return out


# --------------------------------------------------------------------------
# exact simplex (Phase I only: feasibility is all this project needs)
# --------------------------------------------------------------------------

def feasible_nonneg(mat: list[list[Fraction]], rhs: list[Fraction]) -> bool:
    """Exact: is there ``x >= 0`` with ``mat @ x == rhs``?

    Phase-I simplex with Bland's rule, so it terminates without cycling.
    Sizes here are tiny in rows and modest in columns, which is the regime
    Bland's rule is fine in.
    """
    if not mat:
        return all(v == 0 for v in rhs)
    m, n = len(mat), len(mat[0])
    a = [[Fraction(x) for x in row] for row in mat]
    b = [Fraction(v) for v in rhs]
    for i in range(m):  # artificials need a nonnegative right-hand side
        if b[i] < 0:
            a[i] = [-x for x in a[i]]
            b[i] = -b[i]
    # tableau columns: n originals, then m artificials
    tab = [a[i] + [Fraction(1) if j == i else Fraction(0) for j in range(m)] + [b[i]]
           for i in range(m)]
    basis = list(range(n, n + m))
    # cost row for min sum(artificials), expressed in nonbasic variables
    cost = [Fraction(0)] * (n + m + 1)
    for i in range(m):
        for j in range(n + m + 1):
            cost[j] -= tab[i][j]
    for j in range(n, n + m):
        cost[j] = Fraction(0)

    while True:
        enter = next((j for j in range(n + m) if cost[j] < 0), None)
        if enter is None:
            break
        leave, best = None, None
        for i in range(m):
            if tab[i][enter] > 0:
                ratio = tab[i][-1] / tab[i][enter]
                if best is None or ratio < best or (ratio == best and basis[i] < basis[leave]):
                    best, leave = ratio, i
        if leave is None:
            return True  # unbounded below on a nonnegative objective cannot happen
        pv = tab[leave][enter]
        tab[leave] = [x / pv for x in tab[leave]]
        for i in range(m):
            if i != leave and tab[i][enter] != 0:
                f = tab[i][enter]
                tab[i] = [x - f * y for x, y in zip(tab[i], tab[leave])]
        if cost[enter] != 0:
            f = cost[enter]
            cost = [x - f * y for x, y in zip(cost, tab[leave])]
        basis[leave] = enter
    return -cost[-1] == 0


def in_convex_hull(point: Vec, points: list[Vec]) -> bool:
    """Exact test: is ``point`` a convex combination of ``points``?"""
    if not points:
        return False
    n = len(point)
    mat = [[p[i] for p in points] for i in range(n)]
    mat.append([Fraction(1)] * len(points))
    rhs = list(point) + [Fraction(1)]
    return feasible_nonneg(mat, rhs)


# --------------------------------------------------------------------------
# double description
# --------------------------------------------------------------------------

def _tight(rows: list[Row], x: Vec) -> frozenset[int]:
    return frozenset(i for i, (a, c) in enumerate(rows)
                     if sum(ai * xi for ai, xi in zip(a, x)) == c)


def _adjacent(rows: list[Row], eqs: list[Row], v: Vec, w: Vec, n: int) -> bool:
    """Algebraic adjacency: the common tight face must be a single edge."""
    common = _tight(rows, v) & _tight(rows, w)
    mat = [list(eqs_a) for eqs_a, _ in eqs] + [list(rows[i][0]) for i in common]
    return rank(mat) == n - 1


def cut(verts: list[Vec], rows: list[Row], eqs: list[Row], new: Row, n: int) -> list[Vec]:
    """Intersect a V-represented polytope with the halfspace ``new``.

    ``rows`` are the inequalities already imposed (used for the adjacency
    test); ``new`` is appended by the caller after the call.

    The adjacency test dominates the cost, and the naive form recomputes a
    vertex's tight set once per pair and an exact rank once per pair. Three
    things fix that without changing a single answer, since all three are
    identities rather than heuristics:

    - each tight set is computed once per vertex, not once per pair;
    - a pair whose common tight rows cannot reach rank ``n-1`` at all, counted
      against the equalities, is rejected on an integer comparison;
    - the exact rank is memoised on the common tight set, which repeats
      heavily across pairs on a degenerate polytope.
    """
    a, c = new
    val = {v: sum(ai * xi for ai, xi in zip(a, v)) for v in verts}
    keep = [v for v in verts if val[v] >= c]
    drop = [v for v in verts if val[v] < c]
    if not drop:
        return verts

    eq_rows = [list(eq_a) for eq_a, _ in eqs]
    eq_rank = rank(eq_rows)
    tight = {v: _tight(rows, v) for v in verts}
    memo: dict[frozenset[int], bool] = {}

    def adjacent(v: Vec, w: Vec) -> bool:
        common = tight[v] & tight[w]
        if eq_rank + len(common) < n - 1:
            return False
        hit = memo.get(common)
        if hit is None:
            mat = eq_rows + [list(rows[i][0]) for i in common]
            hit = rank(mat) == n - 1
            memo[common] = hit
        return hit

    fresh: list[Vec] = []
    for v in keep:
        if val[v] == c:
            continue
        for w in drop:
            if not adjacent(v, w):
                continue
            # the point on segment [v, w] where <a, x> == c
            theta = (c - val[w]) / (val[v] - val[w])
            pt = tuple(theta * vi + (1 - theta) * wi for vi, wi in zip(v, w))
            fresh.append(pt)
    return sorted(set(keep) | set(fresh), reverse=True)


def vertices_from_hrep(base_verts: list[Vec], base_rows: list[Row], eqs: list[Row],
                       extra: list[Row], n: int) -> list[Vec]:
    """Vertices of ``base`` polytope further cut by the ``extra`` inequalities.

    ``base_verts`` must be the exact vertex set of the polytope described by
    ``base_rows`` together with ``eqs``. Each extra inequality is applied in
    turn by double description.
    """
    verts = list(base_verts)
    rows = list(base_rows)
    for row in extra:
        verts = cut(verts, rows, eqs, row, n)
        rows.append(row)
        if not verts:
            break
    return verts


def brute_force_vertices(rows: list[Row], eqs: list[Row], n: int) -> list[Vec]:
    """Vertex enumeration by tight-subset brute force.

    Only for small systems: it is used to seed the double description with the
    sorted-spectrum simplex, where the inequality count is single digit.
    """
    ineq_rows = [(list(a), c) for a, c in rows]
    eq_rows = [(list(a), c) for a, c in eqs]
    # the equality block can be rank deficient (in the (3,6) system the trace
    # equality is implied by the three pairing equalities), so count its rank
    # rather than its rows, or the tight-subset size comes out too small.
    need = n - rank([a for a, _ in eq_rows])
    out: set[Vec] = set()
    for pick in combinations(range(len(ineq_rows)), need):
        mat = [a for a, _ in eq_rows] + [ineq_rows[i][0] for i in pick]
        rhs = [c for _, c in eq_rows] + [ineq_rows[i][1] for i in pick]
        sol = solve_square(mat, rhs)
        if sol is None:
            continue
        if all(sum(ai * xi for ai, xi in zip(a, sol)) >= c for a, c in ineq_rows):
            out.add(tuple(sol))
    return sorted(out, reverse=True)


def null_space(rows: list[list[Fraction]], n: int) -> list[list[Fraction]]:
    """Basis of ``{a : rows @ a == 0}``, exactly."""
    red, piv = rref(rows) if rows else ([], [])
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        vec = [Fraction(0)] * n
        vec[f] = Fraction(1)
        for i, c in enumerate(piv):
            vec[c] = -red[i][f]
        basis.append(vec)
    return basis


def affine_hull(verts: list[Vec]) -> list[Row]:
    """Equations cutting out the affine hull of a vertex set.

    Reported alongside the facets so a low-dimensional polytope has a complete
    description: without it, an inequality list understates how thin the
    polytope is.
    """
    if not verts:
        return []
    n = len(verts[0])
    origin = verts[0]
    span = [[x - o for x, o in zip(v, origin)] for v in verts[1:]]
    out = []
    for a in null_space(span, n):
        out.append((a, sum(ai * oi for ai, oi in zip(a, origin))))
    return out


def facets(verts: list[Vec], rows: list[Row], eqs: list[Row], n: int) -> list[int]:
    """Indices of the inequalities in ``rows`` that are facets of the polytope.

    An inequality is a facet exactly when the vertices lying on it span an
    affine subspace of dimension one below the polytope's own.
    """
    if not verts:
        return []
    origin = verts[0]
    span = [[x - o for x, o in zip(v, origin)] for v in verts[1:]]
    dim = rank(span) if span else 0
    out = []
    for i, (a, c) in enumerate(rows):
        on = [v for v in verts if sum(ai * xi for ai, xi in zip(a, v)) == c]
        if not on:
            continue
        sub = [[x - o for x, o in zip(v, on[0])] for v in on[1:]]
        if (rank(sub) if sub else 0) == dim - 1:
            out.append(i)
    return out
