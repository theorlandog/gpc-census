"""Exact orbit-closure moment polytopes (entanglement polytopes) for trivectors.

The census layer (``constraints.py``, ``polytope.py``) computes the AMBIENT
moment polytope: every occupation spectrum reachable by some N-fermion pure
state. This module computes the finer object, one polytope per SLOCC class:

    Delta(psi) = { sorted spec(rho_phi) : phi in closure(GL_d . psi) }

which is a convex polytope (Kirwan, Brion, Ness), and is the fermionic
entanglement polytope of Walter, Doran, Gross and Christandl. Delta(psi)
depends only on the orbit, and shrinks as the orbit degenerates, so a
one-body spectrum outside Delta(psi) witnesses that the state is NOT in the
closure of psi's class.

Scope: ``Lambda^3 C^d`` for ``d <= 7``, where the number of GL_d orbits is
finite (Schouten's classification; ten orbits including zero at d = 7).

WHAT IS CERTIFIED, AND HOW
--------------------------
Both bounds are exact and rigorous; nothing here is sampled.

INNER (attainability). If phi lies in the orbit and its support is one-hop
free (no two support triples share two modes), then rho_phi is diagonal, and
so is rho for the whole positive torus orbit of phi. The moment image of the
closure of that torus orbit is exactly the Newton polytope conv{chi_S : S in
supp phi} (the standard toric fact). So every sorted point of that Newton
polytope is an attained spectrum. Each emitted point is therefore a CERTIFIED
member of Delta(psi), with the attaining state written down explicitly.

OUTER (validity). Let W be a subspace with dim W = m, and suppose every
monomial of phi in a basis adapted to W uses at least j factors from W.
Writing rho for the 1-RDM of any g.phi and using that the condition is
GL-covariant, tr(rho restricted to gW) >= j, and Ky Fan bounds that trace by
the sum of the m largest eigenvalues. Hence

    lambda_1 + ... + lambda_m >= j

is valid on Delta(psi), and stays valid on the closure by continuity. For a
coordinate subspace W = span(A) this reads j = min over support triples S of
|S and A|, which is pure support combinatorics. Every inequality emitted here
carries the (representative, A) pair that proves it.

SANDWICH. When the inner hull and the outer H-polytope coincide (checked by
exact rational LP: every vertex of the outer polytope is a convex combination
of certified attainable points), the polytope is EXACT, and the verdict is
recorded as such. Otherwise the result ships as a certified inner/outer
bracket, labelled honestly. No class result is upgraded past what its
certificate proves.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import combinations

from . import exact_polytope as ep

Triple = tuple[int, ...]
Trivector = dict[Triple, Fraction]
Vec = tuple[Fraction, ...]


# --------------------------------------------------------------------------
# trivector algebra over Q
# --------------------------------------------------------------------------

def triples(d: int) -> list[Triple]:
    """Sorted index triples, the standard basis of Lambda^3 C^d."""
    return list(combinations(range(d), 3))


def _det3(m: list[list[Fraction]]) -> Fraction:
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def act(g: list[list[Fraction]], psi: Trivector, d: int) -> Trivector:
    """Act by ``g`` in GL_d: ``(g.psi)_T = sum_S det(g[T, S]) psi_S``."""
    out: Trivector = {}
    items = [(S, c) for S, c in psi.items() if c != 0]
    for T in triples(d):
        acc = Fraction(0)
        for S, c in items:
            acc += _det3([[g[T[a]][S[b]] for b in range(3)] for a in range(3)]) * c
        if acc != 0:
            out[T] = acc
    return out


def transvection(d: int, a: int, b: int, c: Fraction) -> list[list[Fraction]]:
    """The elementary matrix ``I + c E_ab`` (``a != b``), an element of SL_d."""
    if a == b:
        raise ValueError("transvection needs a != b")
    g = [[Fraction(1) if i == j else Fraction(0) for j in range(d)] for i in range(d)]
    g[a][b] = Fraction(c)
    return g


def lie_derivation(a: int, b: int, psi: Trivector) -> Trivector:
    """``E_ab . psi``, the derivation action of gl_d on Lambda^3."""
    out: Trivector = {}
    for S, c in psi.items():
        if b not in S:
            continue
        seq = list(S)
        seq[S.index(b)] = a
        if len(set(seq)) < 3:
            continue
        sgn = 1
        arr = list(seq)
        for i in range(3):
            for j in range(i + 1, 3):
                if arr[i] > arr[j]:
                    arr[i], arr[j] = arr[j], arr[i]
                    sgn = -sgn
        T = tuple(arr)
        out[T] = out.get(T, Fraction(0)) + sgn * c
    return {k: v for k, v in out.items() if v != 0}


def orbit_dimension(psi: Trivector, d: int) -> int:
    """Dimension of the GL_d orbit of ``psi`` (the affine cone, so >= 1).

    The rank of ``X -> X.psi`` on gl_d. For Lambda^3 C^7 the nine nonzero
    orbits have nine distinct dimensions, so this single integer is a complete
    orbit invariant there; ``classify`` relies on that.
    """
    idx = {T: i for i, T in enumerate(triples(d))}
    rows: list[list[Fraction]] = []
    for a in range(d):
        for b in range(d):
            v = [Fraction(0)] * len(idx)
            if a == b:
                for S, c in psi.items():
                    if a in S:
                        v[idx[S]] += c
            else:
                for T, c in lie_derivation(a, b, psi).items():
                    v[idx[T]] += c
            rows.append(v)
    return ep.rank(rows)


def one_rdm(psi: Trivector, d: int) -> list[list[Fraction]]:
    """The unnormalised 1-RDM, ``rho_ij = <a_j psi | a_i psi>``; trace ``3|psi|^2``."""
    ann: list[dict[Triple, Fraction]] = []
    for i in range(d):
        col: dict[Triple, Fraction] = {}
        for S, c in psi.items():
            if i not in S:
                continue
            A = tuple(x for x in S if x != i)
            col[A] = col.get(A, Fraction(0)) + ((-1) ** S.index(i)) * c
        ann.append(col)
    rho = [[Fraction(0)] * d for _ in range(d)]
    for i in range(d):
        for j in range(d):
            rho[i][j] = sum(v * ann[j][A] for A, v in ann[i].items() if A in ann[j])
    return rho


def support_rank(psi: Trivector, d: int) -> int:
    """Least ``r`` with ``psi`` in ``Lambda^3 W``, ``dim W = r``: the rank of rho."""
    return ep.rank(one_rdm(psi, d))


def is_one_hop_free(support) -> bool:
    """No two support triples share exactly two modes (so rho is diagonal)."""
    sup = list(support)
    return all(len(set(a) & set(b)) != 2 for a, b in combinations(sup, 2))


def char_poly(mat: list[list[Fraction]]) -> list[Fraction]:
    """Characteristic polynomial coefficients, low degree first, exactly.

    Faddeev-LeVerrier, which stays inside the rationals and needs no
    factorisation machinery.
    """
    n = len(mat)
    coeffs = [Fraction(0)] * (n + 1)
    coeffs[n] = Fraction(1)
    M = [[Fraction(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        # M <- A M + c I
        AM = [[sum(mat[i][t] * M[t][j] for t in range(n)) for j in range(n)]
              for i in range(n)]
        for i in range(n):
            AM[i][i] += coeffs[n - k + 1]
        M = AM
        trace = sum(sum(mat[i][t] * M[t][i] for t in range(n)) for i in range(n))
        coeffs[n - k] = -trace / k
    return coeffs


def _rational_roots(coeffs: list[Fraction]) -> list[Fraction] | None:
    """All roots of a rational polynomial, or None unless it splits over Q.

    Rational root theorem plus deflation. Eigenvalues of a 1-RDM built from a
    rational state have small height, so divisor enumeration is cheap.
    """
    poly = list(coeffs)
    while poly and poly[-1] == 0:
        poly.pop()
    roots: list[Fraction] = []
    while len(poly) > 1:
        while len(poly) > 1 and poly[0] == 0:  # root at zero
            roots.append(Fraction(0))
            poly.pop(0)
        if len(poly) <= 1:
            break
        den = 1
        for c in poly:
            den = den * c.denominator // _gcd(den, c.denominator)
        ints = [int(c * den) for c in poly]
        g = 0
        for v in ints:
            g = _gcd(g, abs(v))
        if g:
            ints = [v // g for v in ints]
        a0, an = abs(ints[0]), abs(ints[-1])
        found = None
        for p in _divisors(a0):
            for q in _divisors(an):
                for cand in (Fraction(p, q), Fraction(-p, q)):
                    if sum(c * cand ** i for i, c in enumerate(poly)) == 0:
                        found = cand
                        break
                if found is not None:
                    break
            if found is not None:
                break
        if found is None:
            return None
        roots.append(found)
        # synthetic division by (x - found)
        new = [Fraction(0)] * (len(poly) - 1)
        carry = Fraction(0)
        for i in range(len(poly) - 1, 0, -1):
            carry = poly[i] + carry * found
            new[i - 1] = carry
        poly = new
    return roots


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def _divisors(n: int) -> list[int]:
    if n == 0:
        return [1]
    out, i = [], 1
    while i * i <= n:
        if n % i == 0:
            out.append(i)
            if i != n // i:
                out.append(n // i)
        i += 1
    return sorted(out)


def rational_spectrum(psi: Trivector, d: int) -> Vec | None:
    """Exact sorted occupation spectrum when it is rational, else None.

    Unlike ``moment_point`` this does not need a diagonal 1-RDM: it reads the
    spectrum off the characteristic polynomial. That is what lets states with
    a genuinely non-diagonal rho contribute certified points.
    """
    rho = one_rdm(psi, d)
    norm = sum(c * c for c in psi.values())
    if norm == 0:
        return None
    scaled = [[rho[i][j] / norm for j in range(d)] for i in range(d)]
    roots = _rational_roots(char_poly(scaled))
    if roots is None or len(roots) != d:
        return None
    return tuple(sorted(roots, reverse=True))


def moment_point(psi: Trivector, d: int) -> Vec | None:
    """Sorted occupation spectrum, exactly, or None when rho is not diagonal.

    A non-diagonal rho has the same spectrum but generally irrational
    eigenvalues; this module only emits points it can certify exactly, so it
    declines rather than approximating.
    """
    rho = one_rdm(psi, d)
    if any(rho[i][j] != 0 for i in range(d) for j in range(d) if i != j):
        return None
    norm = sum(c * c for c in psi.values())
    if norm == 0:
        return None
    return tuple(sorted((rho[i][i] / norm for i in range(d)), reverse=True))


# --------------------------------------------------------------------------
# the SLOCC class registry
# --------------------------------------------------------------------------

def _tri(*names: str) -> Trivector:
    """Build a trivector from 1-based labels, e.g. ``_tri("123", "456")``."""
    return {tuple(int(ch) - 1 for ch in nm): Fraction(1) for nm in names}


# Representatives are one-hop free, so each carries an exact torus cloud.
# Keyed by orbit dimension, which separates the classes completely.
CLASSES_7: dict[int, dict] = {
    13: {"label": "I", "rank": 3, "rep": _tri("123"),
         "name": "Slater (separable)"},
    20: {"label": "II", "rank": 5, "rep": _tri("123", "145"),
         "name": "rank-5 inherited"},
    21: {"label": "III", "rank": 7, "rep": _tri("123", "145", "167"),
         "name": "frozen mode times rank-6 form"},
    25: {"label": "IV", "rank": 6, "rep": _tri("123", "145", "246"),
         "name": "rank-6 tangent class"},
    26: {"label": "V", "rank": 6, "rep": _tri("123", "456"),
         "name": "rank-6 generic (GHZ-like)"},
    28: {"label": "VI", "rank": 7, "rep": _tri("123", "145", "167", "246"),
         "name": "full-rank, orbit dimension 28"},
    31: {"label": "VII", "rank": 7, "rep": _tri("123", "456", "147"),
         "name": "full-rank, orbit dimension 31"},
    34: {"label": "VIII", "rank": 7, "rep": _tri("123", "456", "147", "257"),
         "name": "full-rank, orbit dimension 34"},
    35: {"label": "IX", "rank": 7, "rep": _tri("123", "456", "147", "257", "367"),
         "name": "generic (G2 form, open orbit)"},
}

# The d = 6 classes, the validation gate against the published four.
CLASSES_6: dict[int, dict] = {
    10: {"label": "I", "rank": 3, "rep": _tri("123"), "name": "Slater (separable)"},
    15: {"label": "II", "rank": 5, "rep": _tri("123", "145"), "name": "rank-5"},
    19: {"label": "III", "rank": 6, "rep": _tri("123", "145", "246"), "name": "tangent class"},
    20: {"label": "IV", "rank": 6, "rep": _tri("123", "456"), "name": "generic (open orbit)"},
}


# The first certified non-toric inner point.  Squared amplitudes are stored as
# integers over a common denominator so the witness remains a small exact
# object.  The signs cancel every one-body coherence.  Its symbolic tangent
# rank is 34, which identifies class VIII because the nine nonzero d = 7
# orbits have distinct dimensions.
CLASS_VIII_INTERFERENCE_WITNESS = {
    "support": ((0, 1, 3), (0, 3, 4), (1, 4, 5),
                (2, 4, 6), (0, 2, 5), (1, 2, 6)),
    "integer_squared_weights": (2, 1, 1, 1, 2, 2),
    "signs": (1, 1, 1, 1, 1, -1),
    "denominator": 9,
    "spectrum": (Fraction(5, 9), Fraction(5, 9), Fraction(5, 9),
                 Fraction(1, 3), Fraction(1, 3), Fraction(1, 3),
                 Fraction(1, 3)),
    "orbit_dimension": 34,
}


def _symbolic_one_rdm(psi, d: int):
    """Exact 1-RDM over an algebraic SymPy coefficient field."""
    import sympy as sp

    ann = []
    for i in range(d):
        col = {}
        for S, c in psi.items():
            if i not in S:
                continue
            A = tuple(x for x in S if x != i)
            col[A] = col.get(A, sp.Integer(0)) + ((-1) ** S.index(i)) * c
        ann.append(col)
    rho = sp.zeros(d)
    for i in range(d):
        for j in range(d):
            rho[i, j] = sp.simplify(sum(
                sp.conjugate(v) * ann[j].get(A, sp.Integer(0))
                for A, v in ann[i].items()
            ))
    return rho


def _symbolic_orbit_dimension(psi, d: int) -> int:
    """Rank of ``gl_d -> Lambda^3 C^d`` over the witness's exact field."""
    import sympy as sp

    idx = {T: i for i, T in enumerate(triples(d))}
    cols = []
    for a in range(d):
        for b in range(d):
            v = [sp.Integer(0)] * len(idx)
            if a == b:
                for S, c in psi.items():
                    if a in S:
                        v[idx[S]] += c
            else:
                for T, c in lie_derivation(a, b, psi).items():
                    v[idx[T]] += c
            cols.append(sp.Matrix(v))
    return int(sp.Matrix.hstack(*cols).rank())


@lru_cache(maxsize=None)
def class_viii_interference_certificate() -> dict:
    """Verify and describe the exact class-VIII interference witness.

    This is deliberately recomputed from the six amplitudes rather than
    trusting the recorded target.  The two load-bearing checks are the exact
    diagonal 1-RDM and the exact orbit tangent rank.
    """
    import sympy as sp

    rec = CLASS_VIII_INTERFERENCE_WITNESS
    den = rec["denominator"]
    psi = {
        S: sign * sp.sqrt(sp.Rational(weight, den))
        for S, weight, sign in zip(
            rec["support"], rec["integer_squared_weights"], rec["signs"]
        )
    }
    norm = sp.simplify(sum(sp.conjugate(c) * c for c in psi.values()))
    rho = _symbolic_one_rdm(psi, 7)
    expected = sp.diag(*[sp.Rational(x.numerator, x.denominator)
                         for x in rec["spectrum"]])
    orbit_dim = _symbolic_orbit_dimension(psi, 7)
    if norm != 1 or rho != expected or orbit_dim != rec["orbit_dimension"]:
        raise AssertionError("class-VIII interference witness failed exact replay")
    return {
        "support": ["".join(str(i + 1) for i in S) for S in rec["support"]],
        "integer_squared_weights": list(rec["integer_squared_weights"]),
        "signs": list(rec["signs"]),
        "denominator": den,
        "spectrum": [str(x) for x in rec["spectrum"]],
        "norm": str(norm),
        "one_rdm": [[str(rho[i, j]) for j in range(7)] for i in range(7)],
        "characteristic_polynomial": str(sp.factor(rho.charpoly().as_expr())),
        "orbit_dimension": orbit_dim,
        "support_rank": int(rho.rank()),
        "one_hop_free": is_one_hop_free(rec["support"]),
    }


def classes(d: int) -> dict[int, dict]:
    if d == 6:
        return CLASSES_6
    if d == 7:
        return CLASSES_7
    raise KeyError(f"no orbit classification tabulated for Lambda^3 C^{d}")


def classify(psi: Trivector, d: int) -> int:
    """Orbit dimension of ``psi``, the key into ``classes(d)``."""
    dim = orbit_dimension(psi, d)
    if dim not in classes(d):
        raise ValueError(f"orbit dimension {dim} is not a tabulated class for d={d}")
    return dim


# --------------------------------------------------------------------------
# orbit exploration: the one-hop-free supports of each class
# --------------------------------------------------------------------------

@lru_cache(maxsize=None)
def onehop_free_supports(d: int) -> tuple[tuple[Triple, ...], ...]:
    """Every nonempty one-hop-free support in ``Lambda^3 C^d``.

    Grown by canonical extension, so each support is produced once. There are
    270 at d = 6 and 5595 at d = 7, the largest being the 7-triple Fano plane.
    """
    T = triples(d)
    out: list[tuple[Triple, ...]] = []
    cur: list[tuple[Triple, ...]] = [()]
    while cur:
        nxt = []
        for s in cur:
            start = T.index(s[-1]) + 1 if s else 0
            for i in range(start, len(T)):
                S = T[i]
                if all(len(set(S) & set(x)) != 2 for x in s):
                    nxt.append(s + (S,))
        out.extend(nxt)
        cur = nxt
    return tuple(out)


def coefficient_moduli(support, d: int) -> int:
    """``|support| - rank`` of the support's incidence matrix.

    Zero means the torus acts transitively on nonzero coefficient vectors over
    that support, so the all-ones trivector represents every state supported
    there, and the support alone pins down the GL_d orbit. This is zero for
    every one-hop-free support at d <= 7, which is what makes the enumeration
    below complete rather than sampled.
    """
    rows = [[Fraction(1) if i in S else Fraction(0) for i in range(d)] for S in support]
    return len(list(support)) - ep.rank(rows)


def _mode_degrees(support, d: int) -> tuple[int, ...]:
    deg = [0] * d
    for S in support:
        for i in S:
            deg[i] += 1
    return tuple(deg)


def canonical_support(support, d: int) -> tuple[Triple, ...]:
    """Canonical relabelling of a support under the mode-permutation group.

    Sorted occupation spectra do not see a relabelling of the modes, so for
    anything computed from sorted spectra (the Newton clouds) it is enough to
    keep one support per permutation class. Only permutations preserving the
    mode-degree sequence can give a smaller form, which cuts the search from
    ``d!`` to the product of the degree-multiplicity factorials.
    """
    deg = _mode_degrees(support, d)
    groups: dict[int, list[int]] = {}
    for i, g in enumerate(deg):
        groups.setdefault(g, []).append(i)
    order = sorted(groups)
    slots: list[int] = []
    for g in order:
        slots.extend(groups[g])
    best = None
    from itertools import permutations as _perms
    pools = [list(_perms(groups[g])) for g in order]

    def combine(idx: int, mapping: dict[int, int]):
        nonlocal best
        if idx == len(pools):
            form = tuple(sorted(tuple(sorted(mapping[i] for i in S)) for S in support))
            if best is None or form < best:
                best = form
            return
        base = sum(len(groups[order[k]]) for k in range(idx))
        for perm in pools[idx]:
            nxt = dict(mapping)
            for pos, src in enumerate(perm):
                nxt[src] = base + pos
            combine(idx + 1, nxt)

    combine(0, {})
    return best


@lru_cache(maxsize=None)
def canonical_supports_by_class(d: int) -> dict[int, tuple[tuple[Triple, ...], ...]]:
    """``supports_by_class`` deduplicated up to relabelling of the modes."""
    out: dict[int, tuple] = {}
    for dim, sups in supports_by_class(d).items():
        seen = {canonical_support(s, d) for s in sups}
        out[dim] = tuple(sorted(seen))
    return out


@lru_cache(maxsize=None)
def support_class_table(d: int) -> dict[frozenset, int]:
    """Lookup from a one-hop-free support to the orbit dimension it realises."""
    return {frozenset(s): dim for dim, sups in supports_by_class(d).items() for s in sups}


@lru_cache(maxsize=None)
def supports_by_class(d: int) -> dict[int, tuple[tuple[Triple, ...], ...]]:
    """One-hop-free supports bucketed by the orbit dimension they realise.

    Complete, not sampled: ``coefficient_moduli`` vanishes on every one-hop-free
    support, so the all-ones trivector is the only state supported there up to
    the torus, and its orbit dimension classifies the support outright. The
    moduli are verified here rather than assumed.
    """
    out: dict[int, list[tuple[Triple, ...]]] = {}
    for sup in onehop_free_supports(d):
        if coefficient_moduli(sup, d) != 0:
            raise AssertionError(f"support {sup} has nonzero coefficient moduli at d={d}")
        psi = {S: Fraction(1) for S in sup}
        out.setdefault(orbit_dimension(psi, d), []).append(sup)
    return {k: tuple(v) for k, v in sorted(out.items())}


# --------------------------------------------------------------------------
# inner bound: certified attainable spectra
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# certified degenerations
# --------------------------------------------------------------------------

def newton_faces(support) -> set[tuple[Triple, ...]]:
    """Sub-supports cut out by a torus limit: the faces of the Newton polytope.

    ``lim_{t->0} diag(t^a).phi`` keeps exactly the triples minimising
    ``<a, chi_S>``, so each face is an explicit one-parameter degeneration.
    """
    sup = sorted(support)
    d = 1 + max(i for S in sup for i in S)
    out: set[tuple[Triple, ...]] = set()
    for a in _weight_probe(d):
        vals = [sum(a[i] for i in S) for S in sup]
        lo = min(vals)
        out.add(tuple(S for S, v in zip(sup, vals) if v == lo))
    return out


def _weight_probe(d: int) -> list[tuple[int, ...]]:
    """Small integer weights, enough to expose every Newton face used here."""
    out = []
    for bits in range(3 ** d):
        a, x = [], bits
        for _ in range(d):
            a.append(x % 3)
            x //= 3
        out.append(tuple(a))
    return out


@lru_cache(maxsize=None)
def degeneration_order(d: int) -> dict[int, tuple[int, ...]]:
    """Classes certified to lie in each class's orbit closure.

    Two constructive sources, each an explicit one-parameter limit:

    - Newton faces of the class's one-hop-free supports (torus limits), and
    - ``pi(h.psi)``, the limit of ``diag(1,...,1,t).h.psi``, which drops every
      monomial using the last mode, scanned over a family of ``h``.

    Closed transitively. This is a certified LOWER bound on the closure order:
    a pair not listed is not thereby excluded, and ``census`` records which
    containments stay open.
    """
    by_class = canonical_supports_by_class(d)
    table = support_class_table(d)
    reach: dict[int, set[int]] = {}
    for dim, meta in sorted(classes(d).items()):
        got = {dim}
        for sup in by_class.get(dim, ()):
            for face in newton_faces(sup):
                if face:
                    got.add(table[frozenset(face)])
        for cs in _projection_probe(d):
            g = [[Fraction(1) if i == j else Fraction(0) for j in range(d)]
                 for i in range(d)]
            for j in range(d - 1):
                g[j][d - 1] = Fraction(cs[j])
            lim = {S: c for S, c in act(g, meta["rep"], d).items() if d - 1 not in S}
            if lim:
                got.add(orbit_dimension(lim, d))
        reach[dim] = got
    changed = True
    while changed:
        changed = False
        for dim in reach:
            grown = set(reach[dim])
            for other in list(reach[dim]):
                grown |= reach.get(other, set())
            if grown != reach[dim]:
                reach[dim] = grown
                changed = True
    return {k: tuple(sorted(v)) for k, v in reach.items()}


def _projection_probe(d: int) -> list[tuple[int, ...]]:
    """Coefficient vectors for the ``e_d -> e_d + v`` rank-drop degeneration."""
    out = [tuple(1 if k == j else 0 for k in range(d - 1)) for j in range(d - 1)]
    out += [tuple(1 for _ in range(d - 1)),
            tuple(j + 1 for j in range(d - 1)),
            tuple((-1) ** j * (j + 1) for j in range(d - 1)),
            tuple(p for p in (2, 3, 5, 7, 11, 13)[:d - 1])]
    return out


DEFAULT_LEVELS = (4, 5, 6, 7)
"""Lattice levels for the Newton clouds.

A level ``L`` grid contains every level dividing ``L``, so ``{4, 5, 6, 7}``
already covers ``L = 1..7``. Level 7 is what the Fano-plane support needs to
reach the maximally mixed spectrum ``(3/7)^7``.
"""


def newton_cloud(support, d: int, levels=DEFAULT_LEVELS) -> set[Vec]:
    """Sorted lattice points of the Newton polytope of a one-hop-free support.

    For weights ``w >= 0`` on the support with ``sum w = 1``, the occupations
    ``lambda_i = sum_{S ni i} w_S`` run over ``conv{chi_S : S in support}``,
    and because the support is one-hop free every one of those points is the
    exact, diagonal, occupation spectrum of the state ``sum_S sqrt(w_S) e_S``
    (the weight is the SQUARED amplitude, which is why the attaining state
    carries a square root while its spectrum stays rational). That state is
    the representative torus-scaled, so it lies in the orbit, and each point
    emitted here is a certified member of the orbit polytope.

    Integer ``w`` with ``sum w = L`` samples the Newton polytope on a lattice;
    the sandwich check, not this generator, is what establishes that enough of
    it was sampled.
    """
    sup = sorted(support)
    k = len(sup)
    out: set[Vec] = set()
    for L in levels:
        for cuts in combinations(range(L + k - 1), k - 1):
            w, prev = [], -1
            for c in cuts:
                w.append(c - prev - 1)
                prev = c
            w.append(L + k - 2 - prev)
            occ = [Fraction(0)] * d
            for wt, S in zip(w, sup):
                if wt:
                    for i in S:
                        occ[i] += Fraction(wt, L)
            out.add(tuple(sorted(occ, reverse=True)))
    return out


def inner_points(orbit_dim: int, d: int,
                 levels=DEFAULT_LEVELS) -> tuple[set[Vec], tuple]:
    """Certified attainable spectra for the orbit closure of a class.

    Every class certified to lie in the closure contributes its own points:
    a degeneration is an explicit limit inside the closure, so its whole
    polytope is included.
    """
    by_class = canonical_supports_by_class(d)
    supports: list = []
    for lower in degeneration_order(d).get(orbit_dim, (orbit_dim,)):
        supports.extend(by_class.get(lower, ()))
    pts: set[Vec] = set()
    for sup in supports:
        pts |= newton_cloud(sup, d, levels)
    if d == 7 and orbit_dim == CLASS_VIII_INTERFERENCE_WITNESS["orbit_dimension"]:
        # Exact non-toric point.  Replay before admitting it to the inner hull,
        # so a stale or edited witness fails closed rather than inflating the
        # certified polytope.
        class_viii_interference_certificate()
        pts.add(CLASS_VIII_INTERFERENCE_WITNESS["spectrum"])
    return pts, tuple(supports)


# --------------------------------------------------------------------------
# outer bound: certified valid inequalities
# --------------------------------------------------------------------------

def dominant_weights(d: int, cmax: int) -> list[tuple[int, ...]]:
    """Decreasing integer weights ``c_1 >= ... >= c_d = 0`` with ``c_1 <= cmax``.

    Normalising ``c_d = 0`` loses nothing: adding a constant to every entry
    shifts both sides of the inequality below by the same trace multiple.
    """
    out: list[tuple[int, ...]] = []

    def rec(prefix: list[int]):
        if len(prefix) == d - 1:
            out.append(tuple(prefix) + (0,))
            return
        hi = cmax if not prefix else prefix[-1]
        for v in range(hi, -1, -1):
            rec(prefix + [v])

    rec([])
    return out


def newton_inequalities(supports, d: int, cmax: int = 3) -> list[dict]:
    """Valid inequalities ``<c, lambda> >= min_S <c, chi_S>`` with their witnesses.

    Let ``c`` be decreasing and let ``phi`` be any state in the orbit, written
    in a basis where it has the given support. Then

        <c, diag(rho_phi)> = sum_S |phi_S|^2 <c, chi_S> >= min_S <c, chi_S>,

    while Schur-Horn gives ``diag(rho_phi)`` majorised by ``lambda``, and
    ``<c, .>`` is Schur-convex for decreasing ``c``, so
    ``<c, diag(rho_phi)> <= <c, lambda>``. Chaining the two proves the
    inequality on the orbit, and continuity carries it to the closure.

    Taking ``c = (1^m, 0^(d-m))`` recovers the Ky Fan / coordinate-subspace
    form ``lambda_1 + ... + lambda_m >= j``. General decreasing ``c`` is
    strictly stronger, because the right-hand side is concave in ``c`` and so
    is not implied by its values on the extreme rays. Relabelling the modes
    changes the right-hand side, so every support in the class is scanned and
    the strongest bound per ``c`` is kept.
    """
    best: dict[tuple[int, ...], tuple[int, tuple]] = {}
    weights = dominant_weights(d, cmax)
    sup_lists = [[[c for c in S] for S in sup] for sup in supports]
    for c in weights:
        if all(v == 0 for v in c):
            continue
        top = None
        arg = None
        for sup, raw in zip(supports, sup_lists):
            val = min(sum(c[i] for i in S) for S in raw)
            if top is None or val > top:
                top, arg = val, sup
        if top:
            best[c] = (top, arg)
    out = []
    for c, (rhs, sup) in sorted(best.items(), key=lambda kv: (-sum(kv[0]), kv[0])):
        out.append({
            "weight": list(c), "rhs": rhs,
            "witness_support": sorted("".join(str(i + 1) for i in S) for S in sup),
        })
    return out


# --------------------------------------------------------------------------
# assembling one class polytope
# --------------------------------------------------------------------------

def _chamber(d: int, n: int, extra_eqs: list[ep.Row]) -> tuple[list[ep.Row], list[ep.Row],
                                                               list[Vec]]:
    """Sorted-spectrum simplex ``1 >= l_1 >= ... >= l_d >= 0``, ``sum = n``.

    ``extra_eqs`` carries any equalities of the ambient constraint system (the
    (3,6) system has three), which have to sit in the equality block so the
    seed vertex enumeration and the adjacency test both see the right
    dimension.
    """
    rows: list[ep.Row] = [([Fraction(-1)] + [Fraction(0)] * (d - 1), Fraction(-1))]
    for i in range(d - 1):
        a = [Fraction(0)] * d
        a[i], a[i + 1] = Fraction(1), Fraction(-1)
        rows.append((a, Fraction(0)))
    a = [Fraction(0)] * d
    a[d - 1] = Fraction(1)
    rows.append((a, Fraction(0)))
    eqs: list[ep.Row] = [([Fraction(1)] * d, Fraction(n))] + list(extra_eqs)
    return rows, eqs, ep.brute_force_vertices(rows, eqs, d)


def ambient_rows(n: int, d: int) -> tuple[list[ep.Row], list[ep.Row]]:
    """The published ambient GPC system, in ``>=`` orientation.

    Returns ``(inequalities, equalities)``. Dropping the equalities silently
    inflates every orbit polytope, so both blocks are carried.
    """
    from .constraints import constraints
    sys_ = constraints(n, d)
    ineqs = [([Fraction(-c) for c in iq["coeffs"]], Fraction(-iq["rhs"]))
             for iq in sys_["inequalities"]]
    eqs = [([Fraction(c) for c in eq["coeffs"]], Fraction(eq["rhs"]))
           for eq in sys_["equalities"]]
    return ineqs, eqs


def orbit_polytope(psi: Trivector, d: int, levels=DEFAULT_LEVELS,
                   cmax: int = 3) -> dict:
    """Exact orbit-closure moment polytope of ``psi``, with its certificates."""
    n = 3
    dim = orbit_dimension(psi, d)
    amb_ineqs, amb_eqs = ambient_rows(n, d)
    base_rows, eqs, base_verts = _chamber(d, n, amb_eqs)
    pts, supports = inner_points(dim, d, levels=levels)

    # Inequalities may only be read off supports of THIS class: a support
    # belonging to a degeneration bounds the smaller polytope, not this one.
    ineqs = newton_inequalities(supports_by_class(d).get(dim, ()), d, cmax=cmax)
    extra: list[ep.Row] = []
    for rec in ineqs:
        extra.append(([Fraction(v) for v in rec["weight"]], Fraction(rec["rhs"])))
    extra += amb_ineqs

    verts = ep.vertices_from_hrep(base_verts, base_rows, eqs, extra, d)
    inner = sorted(pts, reverse=True)

    # sandwich: every outer vertex must be a convex combination of certified
    # attainable spectra, and every certified point must satisfy every
    # constraint (the latter is a consistency check on the certificates).
    unreached = [v for v in verts if not ep.in_convex_hull(v, inner)]
    all_rows = base_rows + extra
    violated = [p for p in inner
                if any(sum(ai * xi for ai, xi in zip(a, p)) < c for a, c in all_rows)
                or any(sum(ai * xi for ai, xi in zip(a, p)) != c for a, c in eqs)]

    exact = not unreached and not violated
    # One inequality per facet: a facet of a polytope is determined by the
    # vertices on it, so distinct inequalities cutting the same face are the
    # same facet and only one is reported.
    facet_idx, seen_faces = [], set()
    for i in ep.facets(verts, all_rows, eqs, d):
        a, c = all_rows[i]
        face = frozenset(j for j, v in enumerate(verts)
                         if sum(ai * xi for ai, xi in zip(a, v)) == c)
        if face not in seen_faces:
            seen_faces.add(face)
            facet_idx.append(i)
    return {
        "orbit_dimension": dim,
        "rank": support_rank(psi, d),
        "representative": sorted("".join(str(i + 1) for i in S) for S in psi),
        "degenerates_to": list(degeneration_order(d).get(dim, ())),
        "supports_explored": len(supports),
        "certified_points": len(inner),
        "interference_witnesses": (
            [class_viii_interference_certificate()]
            if d == 7 and dim == CLASS_VIII_INTERFERENCE_WITNESS["orbit_dimension"]
            else []
        ),
        "vertices": [[str(x) for x in v] for v in verts],
        "newton_inequalities": ineqs,
        "dimension": len(verts) and (d - len(ep.affine_hull(verts))),
        "affine_hull": [_describe(r, d, "==") for r in ep.affine_hull(verts)],
        "facets": [_describe(all_rows[i], d) for i in facet_idx],
        "verdict": "EXACT" if exact else "BRACKET",
        "outer_vertices_not_reached": [[str(x) for x in v] for v in unreached],
        "certificate_violations": [[str(x) for x in v] for v in violated],
    }


def _describe(row: ep.Row, d: int, rel: str = ">=") -> str:
    a, c = row
    terms = " + ".join(f"{ai}*l{i + 1}" for i, ai in enumerate(a) if ai != 0)
    return f"{terms} {rel} {c}"


def census(d: int = 7, levels=DEFAULT_LEVELS, cmax: int = 3) -> dict:
    """Compute the orbit-closure polytope of every SLOCC class of Lambda^3 C^d."""
    out = {}
    for dim, meta in sorted(classes(d).items()):
        rec = orbit_polytope(meta["rep"], d, levels=levels, cmax=cmax)
        rec["label"] = meta["label"]
        rec["name"] = meta["name"]
        if rec["rank"] != meta["rank"]:
            raise AssertionError(f"rank mismatch for class {meta['label']} at d={d}")
        if rec["orbit_dimension"] != dim:
            raise AssertionError(f"orbit dimension mismatch for class {meta['label']}")
        out[str(dim)] = rec
    return out


def polytope_order(result: dict) -> dict[str, list[str]]:
    """Which class polytopes contain which, from the computed vertex sets.

    Only meaningful between classes whose verdict is EXACT; a bracketed class
    contributes its outer polytope, so its containments are reported as
    upper bounds and flagged by the caller.
    """
    verts = {k: [tuple(Fraction(x) for x in v) for v in r["vertices"]]
             for k, r in result.items()}
    out: dict[str, list[str]] = {}
    for a, va in verts.items():
        inside = []
        for b, vb in verts.items():
            if a != b and all(ep.in_convex_hull(p, va) for p in vb):
                inside.append(b)
        out[a] = sorted(inside, key=int)
    return out


def coincidences(result: dict) -> list[list[str]]:
    """Groups of classes sharing one polytope: invisible to one-body spectra."""
    groups: dict[tuple, list[str]] = {}
    for k, r in result.items():
        key = tuple(tuple(v) for v in r["vertices"])
        groups.setdefault(key, []).append(k)
    return [sorted(v, key=int) for v in groups.values() if len(v) > 1]
