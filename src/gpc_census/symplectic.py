"""Cross-realization invariants of the census: lattice cones and stabilizers.

The census describes each extremal object through one presentation: the
occupation spectrum, its support combinatorics, and the facets it saturates.
This module computes two OTHER presentations of the same objects, exactly.

LATTICE presentation (vertices of the moment polytope).
    `tangent_cone_invariants` builds the tangent cone `T_v P` at a vertex from
    the edge directions at that vertex, reduces each ray to its primitive
    generator in the weight lattice `L_0 = {x in Z^d : sum x = 0}` restricted
    to the affine hull of `P`, and reports simpliciality, the Smith invariant
    factors of the ray matrix, the lattice index, and unimodularity.

SYMPLECTIC/LIE presentation (certified states).
    `stabilizer_invariants` computes `dim Stab_{U(d)}([psi])` for the exact
    certified state, the orbit dimension `d^2 - dim Stab`, the commutant
    `G_rho = {U : U rho U^* = rho}` of the state's own 1-RDM, and the dimension
    of the `G_rho` orbit through `[psi]`. Every number is an exact rank of a
    matrix over an algebraic number field, never a numerical rank.

FIBER QUALIFICATION, per the standing rule in `docs/RESEARCH.md`: `G_rho` is
the commutant of the state's OWN 1-RDM, so the stabilizer numbers here are
UNREDUCED and FIXED-RHO. Nothing in this module measures a fixed-spectrum
fiber, and nothing here reports a first-order tangent dimension of a fiber:
at these ranks `ker d(mu)` has dimension at least `2*binom(n,d) - d^2` for
trivial reasons, so such a number is a vacuous upper bound. See
`docs/prereg_symplectic_invariants.md`, "What is deliberately NOT claimed".
"""
from __future__ import annotations

import itertools as it
from fractions import Fraction
from math import gcd

import sympy as sp
from sympy.matrices.normalforms import invariant_factors

# --------------------------------------------------------------------------
# fermionic Fock-space plumbing
# --------------------------------------------------------------------------


def slater_basis(n: int, d: int) -> tuple[list[tuple[int, ...]], dict[tuple[int, ...], int]]:
    """Ordered basis of `wedge^n C^d`: the n-subsets of range(d), lexicographic."""
    b = list(it.combinations(range(d), n))
    return b, {s: i for i, s in enumerate(b)}


def excite(p: int, q: int, s: tuple[int, ...]) -> tuple[int, tuple[int, ...]] | None:
    """Act by `a_p^dagger a_q` on the sorted occupied set `s`.

    Returns `(sign, target)` or None when the result is zero. `s` is sorted, and
    the returned target is sorted, with the fermionic sign tracked through both
    the annihilation and the creation.
    """
    if q not in s:
        return None
    pos = s.index(q)
    sign = -1 if pos % 2 else 1
    rest = s[:pos] + s[pos + 1:]
    if p in rest:
        return None
    if sum(1 for o in rest if o < p) % 2:
        sign = -sign
    return sign, tuple(sorted(rest + (p,)))


def exact_amplitude(pretty: str) -> sp.Expr:
    """Parse one `closed_form.pretty` entry into an exact algebraic number.

    The census writes amplitudes as radicals and as `exp(I*...)` phases built
    from `acos` and `atan` of algebraic arguments. Both families are algebraic;
    `expand_complex` after `expand_func` puts them in `a + b*I` form with `a`
    and `b` radical expressions, which is what exact rank needs.
    """
    e = sp.sympify(pretty)
    return sp.expand(sp.simplify(sp.expand_complex(sp.expand_func(e))))


def state_vector(record: dict):
    """Exact amplitude vector of a census record, in the Slater basis."""
    n, d = (int(x) for x in record["system"].strip("()").split(","))
    basis, index = slater_basis(n, d)
    cf = record["closed_form"]
    psi: list[sp.Expr] = [sp.Integer(0)] * len(basis)
    for pretty, det in zip(cf["pretty"], cf["support_dets"]):
        psi[index[tuple(det)]] = exact_amplitude(pretty)
    return n, d, basis, index, psi


# --------------------------------------------------------------------------
# stabilizer / orbit invariants
# --------------------------------------------------------------------------


def _real_rank(vectors: list[list[sp.Expr]]) -> int:
    """Exact rank over R of complex vectors regarded as real vectors."""
    if not vectors:
        return 0
    rows = [[sp.re(x) for x in v] + [sp.im(x) for x in v] for v in vectors]
    return sp.Matrix(rows).rank()


def stabilizer_invariants(record: dict) -> dict:
    """Exact projective stabilizer and orbit dimensions for one census state.

    All dimensions are real dimensions. The projective orbit dimension is
    `rank_R span{d(pi)(X) psi : X in u(d)} - 1`: the span is real-orthogonal to
    `psi` because `<psi, d(pi)(X) psi>` is purely imaginary for anti-Hermitian
    `X`, and it always contains `i*psi` (take `X = i*I`), which is the direction
    projected out.
    """
    n, d, basis, index, psi = state_vector(record)
    h = len(basis)

    e: dict[tuple[int, int], list[sp.Expr]] = {}
    for p in range(d):
        for q in range(d):
            v = [sp.Integer(0)] * h
            for k, s in enumerate(basis):
                if psi[k] == 0:
                    continue
                out = excite(p, q, s)
                if out is None:
                    continue
                sign, target = out
                v[index[target]] += sign * psi[k]
            e[(p, q)] = [sp.expand(x) for x in v]

    rho = sp.Matrix(
        d, d,
        lambda p, q: sp.expand(sum(sp.conjugate(psi[k]) * e[(p, q)][k] for k in range(h))),
    )
    den = int(record["denominator"] or 1)
    target_rho = sp.diag(*[sp.Rational(int(x), den) for x in record["integer_form"]])
    rho_is_diagonal = sp.simplify(rho - target_rho) == sp.zeros(d, d)

    algebra, images = _u_d_basis(e, d)
    orbit_dim = _real_rank(images) - 1
    stab_dim = d * d - orbit_dim

    commutant = _commutant_coefficients(rho, algebra, d)
    commutant_dim = len(commutant)
    levi = [_combine(images, c) for c in commutant]
    levi_orbit_dim = (_real_rank(levi) - 1) if levi else 0

    mults = _multiplicities(record["integer_form"])
    predicted = sum(m * m for m in mults)

    toral = images[:d]
    toral_orbit_dim = _real_rank(toral) - 1
    toral_stab_dim = d - toral_orbit_dim

    return {
        "rho_is_diagonal": bool(rho_is_diagonal),
        "orbit_dim": int(orbit_dim),
        "stab_dim": int(stab_dim),
        "spectrum_multiplicities": mults,
        "commutant_dim": int(commutant_dim),
        "commutant_dim_matches_multiplicities": bool(commutant_dim == predicted),
        "commutant_orbit_dim": int(levi_orbit_dim),
        "commutant_stab_dim": int(commutant_dim - levi_orbit_dim),
        "toral_orbit_dim": int(toral_orbit_dim),
        "toral_stab_dim": int(toral_stab_dim),
    }


def _multiplicities(integer_form) -> list[int]:
    """Eigenvalue multiplicities of the vertex spectrum, as run lengths."""
    mults: list[int] = []
    for i, value in enumerate(integer_form):
        if i and value == integer_form[i - 1]:
            mults[-1] += 1
        else:
            mults.append(1)
    return mults


def _combine(images: list[list[sp.Expr]], coeffs) -> list[sp.Expr]:
    h = len(images[0])
    return [sp.expand(sum(c * images[j][k] for j, c in enumerate(coeffs) if c != 0))
            for k in range(h)]


def _u_d_basis(e, d):
    """A real basis of `u(d)` and the matching images `d(pi)(X) psi`.

    Order: the `d` diagonal generators `i*E_pp` first, so `images[:d]` is the
    image of the diagonal torus `u(1)^d`, then `E_pq - E_qp` and `i(E_pq+E_qp)`
    for `p < q`. The algebra elements are returned as explicit `d x d` matrices
    so the commutant can be solved for in the same basis.
    """
    algebra: list[sp.Matrix] = []
    images: list[list[sp.Expr]] = []
    for p in range(d):
        m = sp.zeros(d, d)
        m[p, p] = sp.I
        algebra.append(m)
        images.append([sp.I * x for x in e[(p, p)]])
    for p in range(d):
        for q in range(p + 1, d):
            m = sp.zeros(d, d)
            m[p, q], m[q, p] = 1, -1
            algebra.append(m)
            images.append([a - b for a, b in zip(e[(p, q)], e[(q, p)])])
            m = sp.zeros(d, d)
            m[p, q], m[q, p] = sp.I, sp.I
            algebra.append(m)
            images.append([sp.I * (a + b) for a, b in zip(e[(p, q)], e[(q, p)])])
    return algebra, images


def _commutant_coefficients(rho: sp.Matrix, algebra: list[sp.Matrix], d: int):
    """Real basis of `{X in u(d) : [X, rho^T] = 0}`, in coefficients of `algebra`.

    The transpose is not cosmetic. Differentiating `rho_pq = <psi|a_p^* a_q|psi>`
    along the flow of `X` gives

        d(rho) = [rho, X^T],

    so the subalgebra that preserves `rho`, and therefore contains the whole
    projective stabilizer, is the commutant of `rho^T`. For the 645 states whose
    1-RDM is real diagonal the two coincide identically; for the other 154 they
    are different subspaces of the same dimension, and using `rho` there gives a
    commutant that does not contain the stabilizer. The census caught this: the
    forced identity `stab_dim == commutant_stab_dim` failed at `(3,10)#73`.

    Solved as a linear system in the standard orbital basis, so no exact
    diagonalization of `rho` is ever needed.
    """
    rho_t = rho.T
    rows = []
    for x in algebra:
        c = sp.expand(x * rho_t - rho_t * x)
        rows.append([sp.re(c[i, j]) for i in range(d) for j in range(d)]
                    + [sp.im(c[i, j]) for i in range(d) for j in range(d)])
    return [list(v.T) for v in sp.Matrix(rows).T.nullspace()]


def toral_stabilizer_dim_from_support(support_dets: list[list[int]], d: int) -> int:
    """Combinatorial prediction for the toral stabilizer dimension.

    An orbital phase `phi` acts by `c_S -> exp(i sum_{o in S} phi_o) c_S`, so it
    fixes `[psi]` exactly when `A phi` is constant on the support, where `A` is
    the 0/1 support-incidence matrix. The solution space has dimension
    `d - rank([A_S - A_{S_0}])`. No amplitude enters.
    """
    if not support_dets:
        return d
    rows = []
    base = support_dets[0]
    for det in support_dets[1:]:
        row = [0] * d
        for o in det:
            row[o] += 1
        for o in base:
            row[o] -= 1
        rows.append(row)
    if not rows:
        return d
    return d - sp.Matrix(rows).rank()


def support_invariant(support_dets: list[list[int]], d: int) -> str:
    """Permutation invariant of the support hypergraph.

    Built from data that no relabelling of the orbitals can change: the orbital
    degree sequence, the multiset of pairwise determinant overlaps, and, per
    determinant and per orbital, the sorted profile of how it meets everything
    else.

    This is an INVARIANT, not a certified canonical form. Two non-isomorphic
    supports can still collide, which merges classes, so a count of distinct
    values is a LOWER bound on the number of support-isomorphism classes. That
    is the conservative direction for the independent-scope reporting the
    standing rules require: it understates independence rather than inflating
    it.

    An earlier version refined colours and then stored their RANKS. When
    refinement fully discretized, which is the common case here, every support
    of a given size collapsed to the same signature `(0, 1, ..., m-1)`, so the
    invariant reported far fewer classes than exist and the class-based
    scorecard rows were close to meaningless. The version below stores the
    structure itself, not its rank.
    """
    dets = [frozenset(s) for s in support_dets]
    degree = {o: sum(1 for s in dets if o in s) for o in range(d)}
    overlaps = {s: tuple(sorted(len(s & t) for t in dets if t is not s)) for s in dets}
    det_profile = sorted(
        (len(s), tuple(sorted(degree[o] for o in s)), overlaps[s]) for s in dets
    )
    orbital_profile = sorted(
        (degree[o], tuple(sorted(len(s) for s in dets if o in s)),
         tuple(sorted(overlaps[s] for s in dets if o in s)))
        for o in range(d)
    )
    pair_overlaps = sorted(
        len(dets[i] & dets[j])
        for i in range(len(dets)) for j in range(i + 1, len(dets))
    )
    return repr((len(dets), d, det_profile, orbital_profile, pair_overlaps))


# --------------------------------------------------------------------------
# lattice cone invariants
# --------------------------------------------------------------------------


def affine_dim(points) -> int:
    """Exact affine dimension of a set of rational points."""
    if len(points) <= 1:
        return 0
    base = points[0]
    return sp.Matrix([[a - b for a, b in zip(p, base)] for p in points[1:]]).rank()


def candidate_rows(n: int, d: int):
    """`(coeffs, rhs)` for every candidate supporting half-space `a.x <= b`.

    The ordering walls, positivity, and the tabulated GPC inequalities. Which of
    them are facets is decided by `facet_supports`.
    """
    from .constraints import constraints

    rows = []
    for i in range(d - 1):
        c = [Fraction(0)] * d
        c[i], c[i + 1] = Fraction(-1), Fraction(1)
        rows.append((c, Fraction(0)))
    c = [Fraction(0)] * d
    c[d - 1] = Fraction(-1)
    rows.append((c, Fraction(0)))
    for q in constraints(n, d)["inequalities"]:
        rows.append(([Fraction(x) for x in q["coeffs"]], Fraction(q["rhs"])))
    return rows


def facet_supports(points, n: int, d: int, dim: int) -> list[frozenset[int]]:
    """Vertex supports of the facets, deduplicated by support.

    A candidate half-space is a facet exactly when the vertices saturating it
    span an affine subspace of dimension `dim - 1`. Coincident inequalities
    collapse because the support is the key.
    """
    seen: dict[frozenset[int], None] = {}
    for coeffs, rhs in candidate_rows(n, d):
        vals = [sum(c * p for c, p in zip(coeffs, pt)) for pt in points]
        if any(v > rhs for v in vals):
            continue
        on = frozenset(i for i, v in enumerate(vals) if v == rhs)
        if on and affine_dim([points[i] for i in on]) == dim - 1:
            seen.setdefault(on, None)
    return list(seen)


def hull_lattice(points, d: int) -> list[list[int]]:
    """Saturated integer basis of the direction lattice of the affine hull."""
    base = points[0]
    diffs = sp.Matrix([[a - b for a, b in zip(p, base)] for p in points[1:]])
    rows = []
    for v in diffs.nullspace():
        lcm = 1
        for x in v:
            lcm = sp.ilcm(lcm, sp.Rational(x).q)
        rows.append([int(x * lcm) for x in v])
    if not rows:
        return [[1 if i == j else 0 for j in range(d)] for i in range(d)]
    return integer_kernel_basis(rows, d)


def integer_kernel_basis(rows: list[list[int]], d: int) -> list[list[int]]:
    """Saturated integer basis of `{y in Z^d : M y = 0}` for integer `M`.

    Column-style Hermite reduction: unimodular column operations are applied to
    `M` and mirrored on an identity matrix, so the columns of the mirror that
    end up annihilated by `M` are a basis of the kernel lattice. The result is
    saturated because the mirror stays unimodular throughout.
    """
    m = [list(r) for r in rows]
    v = [[1 if i == j else 0 for j in range(d)] for i in range(d)]

    def col_add(target: int, source: int, factor: int) -> None:
        for r in m:
            r[target] += factor * r[source]
        for r in v:
            r[target] += factor * r[source]

    def col_swap(a: int, b: int) -> None:
        for r in m:
            r[a], r[b] = r[b], r[a]
        for r in v:
            r[a], r[b] = r[b], r[a]

    pivot = 0
    for r in range(len(m)):
        nz = [c for c in range(pivot, d) if m[r][c] != 0]
        if not nz:
            continue
        col_swap(pivot, nz[0])
        while True:
            nz = [c for c in range(pivot + 1, d) if m[r][c] != 0]
            if not nz:
                break
            for c in nz:
                col_add(c, pivot, -(m[r][c] // m[r][pivot]))
            nz = [c for c in range(pivot + 1, d) if m[r][c] != 0]
            if nz:
                col_swap(pivot, nz[0])
        pivot += 1
        if pivot == d:
            break
    return [[v[i][c] for i in range(d)] for c in range(pivot, d)]


def primitive(vec: list[Fraction]) -> list[int]:
    """Scale a rational vector to the primitive integer vector on its ray."""
    lcm = 1
    for x in vec:
        lcm = lcm * x.denominator // gcd(lcm, x.denominator)
    ints = [int(x * lcm) for x in vec]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    return [x // g for x in ints] if g else ints


def tangent_cone_invariants(
    verts: list[list[Fraction]],
    which: int,
    facet_supports: list[frozenset[int]],
    lattice: list[list[int]],
    dim: int,
) -> dict:
    """Lattice invariants of the tangent cone at vertex `which`.

    Edges are found combinatorially: `[v, w]` is an edge of `P` exactly when the
    vertices saturating every facet through both `v` and `w` are exactly `v` and
    `w`. The tangent cone of a polytope at a vertex is generated by its edge
    directions, so those rays are the cone.

    `lattice` is a saturated integer basis of the direction lattice of the
    affine hull, so ray coordinates against it are the toric-geometry data:
    invariant factors, lattice index, and unimodularity.
    """
    at_v = {f for f, s in enumerate(facet_supports) if which in s}
    rays: list[list[int]] = []
    neighbours: list[int] = []
    for other in range(len(verts)):
        if other == which:
            continue
        common = at_v & {f for f, s in enumerate(facet_supports) if other in s}
        on_all = set(range(len(verts)))
        for f in common:
            on_all &= set(facet_supports[f])
        if on_all == {which, other}:
            neighbours.append(other)
            rays.append(primitive([a - b for a, b in zip(verts[other], verts[which])]))

    coords = sp.Matrix([_in_lattice(r, lattice) for r in rays]) if rays else sp.zeros(0, dim)
    simplicial = len(rays) == dim
    out = {
        "active_facets": len(at_v),
        "edges": len(rays),
        "simple": len(at_v) == dim,
        "simplicial_tangent_cone": simplicial,
        "neighbours": neighbours,
        "rays": rays,
        "ray_rank": int(coords.rank()) if rays else 0,
    }
    if simplicial:
        det = int(coords.det())
        out["lattice_index"] = abs(det)
        out["unimodular"] = abs(det) == 1
        out["invariant_factors"] = [int(x) for x in invariant_factors(coords, domain=sp.ZZ)]
    else:
        out["lattice_index"] = None
        out["unimodular"] = False
        out["invariant_factors"] = (
            [int(x) for x in invariant_factors(coords, domain=sp.ZZ)] if rays else [])
    return out


def _in_lattice(vec: list[int], lattice: list[list[int]]) -> list[int]:
    """Coordinates of a lattice vector against a saturated basis, exactly."""
    b = sp.Matrix(lattice).T
    sol = b.solve(sp.Matrix(vec))
    return [int(x) for x in sol]


# --------------------------------------------------------------------------
# finite isotropy characters and the quantization period
# --------------------------------------------------------------------------


def spectrum_blocks(integer_form) -> list[list[int]]:
    """Index blocks of equal occupation, the Levi factors of `L_lambda`."""
    blocks: list[list[int]] = []
    for i, value in enumerate(integer_form):
        if i and value == integer_form[i - 1]:
            blocks[-1].append(i)
        else:
            blocks.append([i])
    return blocks


def _sort_sign(values: list[int]) -> int:
    """Sign of the permutation that sorts `values` ascending."""
    sign = 1
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            if values[i] > values[j]:
                sign = -sign
    return sign


def _cycle_sign(perm: list[int]) -> int:
    seen = [False] * len(perm)
    sign = 1
    for i in range(len(perm)):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            j, length = perm[j], length + 1
        if length % 2 == 0:
            sign = -sign
    return sign


def levi_permutation_stabilizer(record: dict, max_period_multiple: int = 512) -> list[dict]:
    """`lambda`-preserving permutations that fix the certified state's LINE.

    Enumerates `sigma` in `S_d` preserving every block of equal occupation,
    which is exactly `S_d` intersected with the Levi `L_lambda`, and keeps those
    with `wedge^n(P_sigma) psi = c psi` for a scalar `c`. For each it records
    `c`, the block determinants (the sign of `sigma` restricted to each block),
    and the least `k` divisible by the spectrum denominator `q` with

        c^k * prod_i det_i^(-k * a_i / q) = 1,

    which is the character condition `omega(sigma)^k = 1` of
    `docs/prereg_quantization_character.md`. All arithmetic is exact.
    """
    n, d, basis, index, psi = state_vector(record)
    form = [int(x) for x in record["integer_form"]]
    q = int(record["denominator"] or 1)
    blocks = spectrum_blocks(form)

    support = {tuple(sorted(det)) for det in record["closed_form"]["support_dets"]}

    out: list[dict] = []
    for choice in it.product(*[list(it.permutations(b)) for b in blocks]):
        sigma = list(range(d))
        for block, image in zip(blocks, choice):
            for src, dst in zip(block, image):
                sigma[src] = dst

        # `P_sigma` is invertible, so `wedge^n(P_sigma) psi` is nonzero and can be
        # a multiple of `psi` only if `sigma` permutes the support. Testing that
        # first is a pure set operation and skips the exact amplitude work for
        # the overwhelming majority of candidates.
        if {tuple(sorted(sigma[o] for o in det)) for det in support} != support:
            continue

        moved = [sp.Integer(0)] * len(basis)
        for k, s in enumerate(basis):
            if psi[k] == 0:
                continue
            image = [sigma[o] for o in s]
            moved[index[tuple(sorted(image))]] += _sort_sign(image) * psi[k]

        scalar = None
        ok = True
        for k in range(len(basis)):
            if psi[k] == 0:
                if sp.simplify(moved[k]) != 0:
                    ok = False
                    break
                continue
            ratio = sp.simplify(moved[k] / psi[k])
            if scalar is None:
                scalar = ratio
            elif sp.simplify(ratio - scalar) != 0:
                ok = False
                break
        if not ok or scalar is None:
            continue

        dets = [_cycle_sign([block.index(sigma[o]) for o in block]) for block in blocks]
        # The exponent on block i is k*lambda_i, and lambda_i is the occupation
        # shared by the ORBITALS of that block, not `form[i]`: `form` is indexed
        # by orbital, `dets` by block. Indexing `form` by the block counter is
        # wrong as soon as the blocks are not all singletons, and it manufactures
        # obstructions that do not exist.
        block_numerators = [form[block[0]] for block in blocks]
        least = None
        for k in range(q, max_period_multiple * q + 1, q):
            value = scalar ** k
            for i, det in enumerate(dets):
                value *= sp.Integer(det) ** (-k * block_numerators[i] // q)
            if sp.simplify(value - 1) == 0:
                least = k
                break
        out.append({
            "sigma": sigma,
            "c": sp.srepr(sp.nsimplify(scalar)),
            "c_pretty": str(sp.nsimplify(scalar)),
            "block_determinants": dets,
            "least_admissible_k": least,
            "nontrivial": sigma != list(range(d)),
        })
    return out


def predicted_quantization_period(elements: list[dict], q: int) -> int | None:
    """lcm of the least admissible `k` over the isotropy elements found."""
    period = q
    for element in elements:
        least = element["least_admissible_k"]
        if least is None:
            return None
        period = period * least // gcd(period, least)
    return period


# --------------------------------------------------------------------------
# the symplectic slice (Marle-Guillemin-Sternberg local model)
# --------------------------------------------------------------------------


def symplectic_slice(record: dict) -> dict:
    """Exact symplectic slice at a certified state, with its validation gates.

    At `x = [psi]` in `P(wedge^n C^d)` with `K = U(d)` and moment map the 1-RDM,
    the symplectic slice is

        V_x = T_x(K.x)^omega / T_x(K_mu(x).x),

    where `T_x(K.x)^omega = ker d(mu)_x` and the denominator is the standard
    intersection `T_x(K.x)^omega and T_x(K.x)`. This is the object the previous
    campaign lacked: the RAW `ker d(mu)_x` is dominated by trivial ambient
    directions and is a vacuous upper bound on anything, while `dim V_x` is the
    honest local transverse dimension, and `dim V_x = 0` certifies outright that
    the reduced germ is a point.

    Two identities are checked on every state and returned, because both are
    forced and both caught bugs while this was written:

      `dim ker d(mu)_x == 2(H-1) - dim(K.x)`, since `im d(mu)_x = ann(k_x)`;
      `dim T_x(K_mu.x) == dim K_mu - dim K_x`, since `K_x` is contained in
      `K_mu`.

    NOT returned: a reduced-germ dimension. The obvious shortcut,
    `dim V_x - 2 dim(K_x . v)` at a generic `v` in `V_x`, is WRONG, because the
    formula needs a generic point of `mu_V^{-1}(0)` and a generic point of
    `V_x` has a strictly larger orbit. At the Slater vertex it returns -10.
    See `docs/symplectic_slice.md`.
    """
    n, d, basis, index, psi = state_vector(record)
    h = len(basis)

    e: dict[tuple[int, int], list[sp.Expr]] = {}
    for p in range(d):
        for q in range(d):
            v = [sp.Integer(0)] * h
            for k, s in enumerate(basis):
                if psi[k] == 0:
                    continue
                out = excite(p, q, s)
                if out:
                    v[index[out[1]]] += out[0] * psi[k]
            e[(p, q)] = [sp.expand(x) for x in v]

    algebra, images = _u_d_basis(e, d)
    rho = sp.Matrix(
        d, d,
        lambda p, q: sp.expand(sum(sp.conjugate(psi[k]) * e[(p, q)][k] for k in range(h))),
    )
    commutant = _commutant_coefficients(rho, algebra, d)

    real_psi = sp.Matrix([[sp.re(x) for x in psi] + [sp.im(x) for x in psi]])
    real_ipsi = sp.Matrix([[sp.re(sp.I * x) for x in psi] + [sp.im(sp.I * x) for x in psi]])
    n1, n2 = real_psi.dot(real_psi), real_ipsi.dot(real_ipsi)

    def to_tangent(vec):
        """Real coordinates of `vec`, projected onto the tangent space `psi^perp`."""
        r = sp.Matrix([[sp.re(x) for x in vec] + [sp.im(x) for x in vec]])
        r = r - (r.dot(real_psi) / n1) * real_psi - (r.dot(real_ipsi) / n2) * real_ipsi
        return [sp.simplify(x) for x in r]

    orbit = sp.Matrix([to_tangent(g) for g in images])
    orbit_dim = orbit.rank()
    stabilizer = [list(v.T) for v in orbit.T.nullspace()]

    ambient = []
    for k in range(h):
        for mult in (sp.Integer(1), sp.I):
            unit = [sp.Integer(0)] * h
            unit[k] = mult
            ambient.append(to_tangent(unit))
    span = sp.Matrix(ambient)
    # pivot COLUMNS of the transpose are the independent ROWS; using rref pivots
    # on `span` itself selects the wrong vectors and silently drops the rank of
    # d(mu) by one, which is how this was caught at the Slater vertex.
    tangent = [ambient[i] for i in span.T.rref()[1]]
    assert len(tangent) == 2 * h - 2, (len(tangent), 2 * h - 2)

    def d_mu(row):
        vec = [row[k] + sp.I * row[h + k] for k in range(h)]
        out = []
        for p in range(d):
            for q in range(p, d):
                a = sum(sp.conjugate(vec[j]) * e[(p, q)][j] for j in range(h))
                b = sum(sp.conjugate(vec[j]) * e[(q, p)][j] for j in range(h))
                val = sp.expand(a + sp.conjugate(b))
                out.append(sp.re(val))
                if p != q:
                    out.append(sp.im(val))
        return out

    jac = sp.Matrix([d_mu(t) for t in tangent])
    basis_matrix = sp.Matrix(tangent)
    kernel = [list(sp.Matrix([list(v.T)]) * basis_matrix) for v in jac.T.nullspace()]

    def combine_images(coeffs):
        out = [sp.Integer(0)] * h
        for j, c in enumerate(coeffs):
            if c == 0:
                continue
            for k in range(h):
                out[k] += c * images[j][k]
        return to_tangent(out)

    levi = sp.Matrix([combine_images(c) for c in commutant])
    levi_dim = levi.rank()

    reduced = sp.Matrix(kernel)
    if levi_dim:
        independent = sp.Matrix([levi.row(i) for i in levi.T.rref()[1]])
        projected = []
        for row in kernel:
            r = sp.Matrix([row])
            for i in range(independent.rows):
                q = independent.row(i)
                r = r - (r.dot(q) / q.dot(q)) * q
            projected.append([sp.simplify(x) for x in r])
        reduced = sp.Matrix(projected)
    slice_dim = len(reduced.T.rref()[1]) if kernel else 0

    return {
        "ambient_tangent_dim": 2 * h - 2,
        "orbit_dim": int(orbit_dim),
        "stab_dim": int(d * d - orbit_dim),
        "stabilizer_algebra_dim": len(stabilizer),
        "commutant_dim": len(commutant),
        "rank_d_mu": int(jac.rank()),
        "ker_d_mu_dim": len(kernel),
        "levi_orbit_dim": int(levi_dim),
        "slice_dim": int(slice_dim),
        "gate_rank_equals_orbit": bool(jac.rank() == orbit_dim),
        "gate_kernel_dim": bool(len(kernel) == 2 * h - 2 - orbit_dim),
        "gate_levi_dim": bool(levi_dim == len(commutant) - (d * d - orbit_dim)),
        "slice_certifies_point": bool(slice_dim == 0),
    }


def toral_component_period(record: dict, cap: int = 64) -> int | None:
    """Obstruction from the DISCONNECTED part of the diagonal stabilizer.

    `levi_permutation_stabilizer` searches permutations and so misses the
    component group of the torus itself. Duality supplies it without ever
    solving for a phase. Writing `u = (z, c)` in `U(1)^(d+1)`, a diagonal `h`
    fixes the line exactly when `u^B = 1`, where `B` has rows `[A_S | -1]` over
    the support, and the obstruction character is `u^(v_k)` with
    `v_k = (-k lambda, k)`. By Pontryagin duality that character is trivial on
    the whole stabilizer exactly when

        v_k lies in the INTEGER rowspan of B,

    and when it is not, it takes at least two values, so `k` is inadmissible.
    The integer rowspan, not its saturation, is the point: the gap between them
    is precisely the component group.

    SCOPE, and it is not optional. This is valid only where the certified
    state's 1-RDM is DIAGONAL, because the derivation needs `lambda` to be the
    diagonal of `rho`, i.e. `lambda_o = sum_{S containing o} |psi_S|^2`. That
    identity is what puts `v_k` in the RATIONAL rowspan in the first place:
    `[lambda | -1] = sum_S |psi_S|^2 [A_S | -1]`. Measured across the census,
    all 645 diagonal states satisfy it and all 49 states that fail it are among
    the 154 whose 1-RDM is not diagonal, with zero counterexamples. For those
    154 the support incidence matrix is written in the wrong basis and this
    function returns None rather than a wrong answer.

    Returns the least `k` divisible by the spectrum denominator that clears the
    obstruction, or None when out of scope or beyond `cap` multiples.
    """
    from sympy.matrices.normalforms import hermite_normal_form

    n, d = (int(x) for x in record["system"].strip("()").split(","))
    q = int(record["denominator"] or 1)
    form = [int(x) for x in record["integer_form"]]

    rows = []
    for det in record["closed_form"]["support_dets"]:
        row = [0] * (d + 1)
        for o in det:
            row[o] += 1
        row[d] = -1
        rows.append(row)
    m = sp.Matrix(rows).T

    probe = sp.Matrix([[-q * a // q for a in form] + [q]]).T
    if m.rank() != sp.Matrix.hstack(m, probe).rank():
        return None                       # out of scope: rho is not diagonal

    base = hermite_normal_form(m)
    for step in range(1, cap + 1):
        k = q * step
        v = sp.Matrix([[-k * a // q for a in form] + [k]]).T
        if hermite_normal_form(sp.Matrix.hstack(m, v)) == base:
            return k
    return None
