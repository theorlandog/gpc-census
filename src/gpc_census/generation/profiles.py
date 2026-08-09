"""Profile-native enumeration of admissible exterior-weight hyperplanes.

WHAT THIS REPLACES.  ``orbit_native.enumerate_orbit_survivor_taus`` is already
orbit native downstream of its first stage, but that first stage is
``orbit_canonical.enumerate_orbits_by_augmentation``, which walks the whole
lattice of closed flats one rank at a time.  The lattice is what costs: at
``(3,9)`` the top level holds 231 of 1,297 orbits and 78 percent of the time,
at ``(3,10)`` the peak level holds 2,705 orbits, and ``(3,11)`` has not been
finished by that route.  The orbits themselves are few.  Only the road to them
is long.

WHAT THIS BUILDS INSTEAD.  A codimension-one closed flat is a hyperplane, and
two hyperplanes lie in the same ``S_d`` orbit exactly when their primitive
normals are permutations of one another.  So an orbit IS a multiset of integers,
and the lattice never has to be entered at all.  This module enumerates those
multisets directly:

    value profile (multiplicities, values)
      -> exact admissibility test in R^r rather than R^d      (Theorem A)
      -> level reparametrization, which makes the search finite  (Lemma L)
      -> one exact (h, z) per top-level orbit

Downstream is untouched: ``orbit_native`` consumes orbit normals and emits the
taus that pass the Ressayre trace test.

THE TWO STATEMENTS THIS RESTS ON, both proved in
``docs/profile_native_generation.md`` and both checked here against the literal
definition rather than assumed.

Theorem A (profile reduction).  Let ``tau`` in Z^d have distinct values
``v_1 > ... > v_r`` with multiplicities ``m_1, ..., m_r``.  Let

    P(m) = { c in Z_{>=0}^r : sum_a c_a = N, c_a <= m_a }
    S    = { c in P(m) : sum_a c_a v_a = 0 }

Then the exterior weights lying on ``tau . x = 0`` have affine rank ``d - 1``,
which is exactly the Ressayre admissibility condition, if and only if

    (i)  rank_Q(S) = r - 1, and
    (ii) every class ``a`` with ``m_a >= 2`` is MOBILE, meaning that some
         ``c`` in ``S`` has ``1 <= c_a <= m_a - 1``.

Both sides are computed here, and ``tests/test_profiles.py`` pins their
agreement on every profile the enumerator visits at ranks 6 through 8 for
``N = 3`` and ``N = 4``.  Condition (ii) is not decoration: dropping it admits
``tau = (1,1,0,0,0,-2)`` at ``d = 6``, which has ``rank(S) = 2 = r - 1`` and
weight rank 2, not 5.

Lemma L (level reparametrization).  Suppose ``S`` is nonempty and ``tau`` is
primitive.  Put ``B = gcd_{a,b}(v_a - v_b)`` and ``l_a = (v_1 - v_a) / B``.  Then
the ``l_a`` are integers with ``l_1 = 0``, strictly increasing, of gcd one, every
``c`` in ``S`` has the SAME level sum ``t = N v_1 / B``, and

    v_a  is proportional to  t - N * l_a.

So a profile is named by a multiplicity composition, an increasing integer level
set, and one integer target.  A corollary drops out for free: ``gcd(v_1, B) = 1``
and ``B | N v_1`` give ``B | N``, which is the divisibility lemma of
``docs/arithmetic_level_audit.md``, so for ``N = 3`` the value gap is 1 or 3.

Theorem S (determination, and the spread bound).  Let ``V`` be the span of the
differences ``c - c'`` of elements of ``S``.  Condition (i) says exactly that
``dim V = r - 2``, and ``V`` is contained in both ``1^perp`` and ``v^perp``,
which already cut out an ``(r-2)``-dimensional space.  So

    V = 1^perp  n  v^perp,   hence   V^perp = span{1, v}.

The pattern set therefore DETERMINES the values, up to the reparametrization
``v -> a v + b 1`` that Lemma L already quotients out.  Since ``S`` is a subset
of the finite set ``P(m)``, which depends on ``m`` and not on the values, each
multiplicity vector carries finitely many admissible level vectors, and the
level SPREAD is bounded at every rank.  Cramer and Hadamard make that explicit:

    spread  <=  2 sqrt(r - 1) (n sqrt 2)^(r - 2),

which is ``spread_bound`` below.  ``level_vector_from_patterns`` runs the
determination backwards, recovering the levels from the pattern set alone.

WHAT THIS DOES AND DOES NOT LICENSE, stated once and carried into every artifact
this module writes.  The spread is bounded, so the search is finite and
candidate coverage is decidable at every rank.  The bound is far too large to
enumerate against: at ``n = 3, r = 11`` it is about ``2.8e6`` against a realized
maximum of 20.  The enumerator therefore still takes ``max_spread`` as an
argument, and its output is complete only for profiles within it, so ranks whose
``max_spread`` sits below ``spread_bound`` still report lower bounds.  See
``spread_saturation`` for the measurement protocol.  The realized maxima are
3, 4, 6, 10 and 14 at ranks 6 through 10, against known-complete orbit counts of
8, 19, 56, 231 and 1337.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from .exterior import exterior_weights


@dataclass(frozen=True, order=True)
class ValueProfile:
    """One ``S_d`` orbit of admissible hyperplanes, named by its normal multiset.

    ``multiplicities`` and ``values`` are aligned and ordered by strictly
    decreasing value.  ``levels`` and ``level_target`` are the Lemma L
    coordinates and are recorded so that a reader can replay the
    reparametrization without redoing the gcd.
    """

    particle_number: int
    multiplicities: tuple[int, ...]
    values: tuple[int, ...]
    levels: tuple[int, ...]
    level_target: int

    @property
    def rank(self) -> int:
        """Ambient rank ``d``, the total multiplicity."""
        return sum(self.multiplicities)

    @property
    def classes(self) -> int:
        """Number of distinct values ``r``."""
        return len(self.values)

    @property
    def spread(self) -> int:
        """Level spread ``l_r``, the quantity ``max_spread`` bounds."""
        return self.levels[-1]

    @property
    def tau(self) -> tuple[int, ...]:
        """The sorted primitive normal, one entry per mode."""
        out: list[int] = []
        for multiplicity, value in zip(self.multiplicities, self.values, strict=True):
            out.extend([value] * multiplicity)
        return tuple(out)

    def as_dict(self) -> dict[str, object]:
        return {
            "particle_number": self.particle_number,
            "multiplicities": list(self.multiplicities),
            "values": list(self.values),
            "levels": list(self.levels),
            "level_target": self.level_target,
            "rank": self.rank,
            "classes": self.classes,
            "spread": self.spread,
            "tau": list(self.tau),
        }


def _content(values) -> int:
    divisor = 0
    for value in values:
        divisor = math.gcd(divisor, abs(value))
    return divisor


def _primitive(values: tuple[int, ...]) -> tuple[int, ...]:
    divisor = _content(values)
    if divisor == 0:
        raise ValueError("the zero vector is not a hyperplane normal")
    return tuple(value // divisor for value in values)


@lru_cache(maxsize=None)
def occupancy_patterns(n: int, r: int) -> tuple[tuple[int, ...], ...]:
    """Every class occupancy ``c`` with ``sum(c) = n``, before multiplicity caps.

    A pattern records how many of the ``n`` occupied modes a weight takes from
    each value class.  Capping by ``m`` is done later, because the uncapped set
    depends only on ``r`` and is worth caching across compositions.
    """
    if type(n) is not int or type(r) is not int or n < 1 or r < 1:
        raise ValueError("require integer parameters with n >= 1 and r >= 1")
    out: list[tuple[int, ...]] = []
    for split in itertools.combinations_with_replacement(range(r), n):
        pattern = [0] * r
        for index in split:
            pattern[index] += 1
        out.append(tuple(pattern))
    return tuple(sorted(set(out)))


def capped_patterns(n: int, multiplicities: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Class occupancies realizable when class ``a`` has only ``m_a`` modes."""
    return tuple(
        pattern
        for pattern in occupancy_patterns(n, len(multiplicities))
        if all(count <= cap for count, cap in zip(pattern, multiplicities, strict=True))
    )


def zero_patterns(
    n: int, multiplicities: tuple[int, ...], values: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    """The set ``S`` of Theorem A: capped occupancies of zero value sum."""
    return tuple(
        pattern
        for pattern in capped_patterns(n, multiplicities)
        if sum(count * value for count, value in zip(pattern, values, strict=True)) == 0
    )


_RANK_PRIME = (1 << 61) - 1


def _rank_mod(rows, columns: int) -> int:
    """Rank over ``GF(_RANK_PRIME)``.

    The prime exceeds every relevant minor: a pattern row has entries summing to
    ``n``, so an ``r`` by ``r`` minor is at most ``n ** r`` by Hadamard, and the
    enumerator is only ever run with ``n ** r`` far below the prime.  A modular
    rank is then the rational rank.
    """
    matrix = [[value % _RANK_PRIME for value in row] for row in rows]
    rank = 0
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]), None
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], _RANK_PRIME - 2, _RANK_PRIME)
        head = [value * inverse % _RANK_PRIME for value in matrix[rank]]
        matrix[rank] = head
        for index in range(rank + 1, len(matrix)):
            factor = matrix[index][column]
            if factor:
                matrix[index] = [
                    (left - factor * right) % _RANK_PRIME
                    for left, right in zip(matrix[index], head, strict=True)
                ]
        rank += 1
        if rank == columns:
            break
    return rank


def rank_modulus_bound(n: int, r: int) -> tuple[int, int]:
    """Return the rank prime and the Hadamard bound it has to beat."""
    return _RANK_PRIME, n**r


def profile_is_admissible(
    n: int, multiplicities: tuple[int, ...], values: tuple[int, ...]
) -> bool:
    """Theorem A, evaluated in ``R^r`` rather than in ``R^d``."""
    r = len(multiplicities)
    if r != len(values):
        raise ValueError("multiplicities and values must be aligned")
    if len(set(values)) != r:
        raise ValueError("values must be distinct")
    if n**r >= _RANK_PRIME:
        raise ValueError("rank prime does not exceed the Hadamard bound")
    patterns = zero_patterns(n, multiplicities, values)
    if len(patterns) < r - 1:
        return False
    if _rank_mod(patterns, r) != r - 1:
        return False
    for index, cap in enumerate(multiplicities):
        if cap >= 2 and not any(1 <= pattern[index] <= cap - 1 for pattern in patterns):
            return False
    return True


def weight_affine_rank(n: int, tau: tuple[int, ...]) -> int:
    """Reference definition: affine rank of the weights on ``tau . x = 0``.

    This materializes ``binom(d, n)`` weights and exists only so that Theorem A
    can be checked against the thing it claims to compute.  The enumerator never
    calls it.
    """
    d = len(tau)
    rows = [
        tuple(1 if index in subset else 0 for index in range(d))
        for subset in itertools.combinations(range(d), n)
        if sum(tau[index] for index in subset) == 0
    ]
    if not rows:
        return 0
    return _rank_mod(rows, d)


# --- Theorem S: the pattern set determines the values ------------------------


def _nullspace(rows, columns: int) -> tuple[tuple[int, ...], ...]:
    """Exact rational nullspace of an integer matrix, as primitive integer vectors.

    Exact rather than modular: ``_rank_mod`` only has to certify a rank, but a
    recovered level vector is an answer, and a modular one would be a residue.
    """
    matrix = [[Fraction(value) for value in row] for row in rows]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (index for index in range(pivot_row, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        head = matrix[pivot_row]
        for index in range(len(matrix)):
            factor = matrix[index][column]
            if index != pivot_row and factor:
                matrix[index] = [
                    left - factor * right
                    for left, right in zip(matrix[index], head, strict=True)
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    basis: list[tuple[int, ...]] = []
    for free in (column for column in range(columns) if column not in pivots):
        vector = [Fraction(0)] * columns
        vector[free] = Fraction(1)
        for index, column in enumerate(pivots):
            vector[column] = -matrix[index][free]
        denominator = 1
        for value in vector:
            denominator = denominator * value.denominator // math.gcd(
                denominator, value.denominator
            )
        basis.append(_primitive(tuple(int(value * denominator) for value in vector)))
    return tuple(basis)


def level_vector_from_patterns(
    patterns: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], int] | None:
    """Theorem S run backwards: recover ``(levels, target)`` from ``S`` alone.

    The values never enter.  Returns ``None`` when the pattern set fails
    condition (i) of Theorem A, or when the recovered levels are not strictly
    increasing, which is the same failure seen from the other side: a pattern set
    that does not pin a hyperplane is not the zero set of one.
    """
    if not patterns:
        return None
    r = len(patterns[0])
    if any(len(pattern) != r for pattern in patterns):
        raise ValueError("patterns must all have one entry per value class")
    base = patterns[0]
    differences = [
        tuple(left - right for left, right in zip(pattern, base, strict=True))
        for pattern in patterns[1:]
    ]
    basis = _nullspace(differences, r)
    if len(basis) != 2:
        return None
    ones = (1,) * r
    witness = next(
        (vector for vector in basis if vector != ones and vector != (-1,) * r), None
    )
    if witness is None:
        return None
    raw = tuple(witness[0] - value for value in witness)
    divisor = _content(raw)
    if divisor == 0:
        return None
    levels = tuple(value // divisor for value in raw)
    if all(left > right for left, right in zip(levels, levels[1:])):
        levels = tuple(-value for value in levels)
    if any(left >= right for left, right in zip(levels, levels[1:])):
        return None
    targets = {
        sum(count * level for count, level in zip(pattern, levels, strict=True))
        for pattern in patterns
    }
    if len(targets) != 1:
        return None
    return levels, targets.pop()


def spread_bound(n: int, r: int) -> int:
    """Theorem S: no admissible profile with ``r`` classes exceeds this spread.

    Pick ``r - 2`` independent pattern differences as the rows of ``A`` and
    append the all-ones row.  The kernel of the result is the line through
    ``r v - (sum v) 1``, so Cramer writes a generator as the signed maximal
    minors of an ``(r-1)`` by ``r`` matrix whose pattern rows have Euclidean
    norm at most ``n sqrt 2`` (entries bounded by ``n``, absolute values summing
    to at most ``2n``) and whose last row has norm ``sqrt(r - 1)``.  Hadamard
    bounds each minor, and the spread is at most the difference of the first and
    last coordinates of that generator.

    True but enormous: 12 at ``(3, 3)``, about ``2.8e6`` at ``(3, 11)`` against a
    realized maximum of 20.  It settles finiteness, not tractability.
    """
    if type(n) is not int or type(r) is not int or n < 1 or r < 2:
        raise ValueError("require integer parameters with n >= 1 and r >= 2")
    squared = 4 * (r - 1) * 2 ** (r - 2) * n ** (2 * (r - 2))
    root = math.isqrt(squared)
    return root if root * root == squared else root + 1


def determination_holds(n: int, profile: ValueProfile) -> bool:
    """Check Theorem S on one profile: does ``S`` alone give back its levels?"""
    patterns = zero_patterns(n, profile.multiplicities, profile.values)
    return level_vector_from_patterns(patterns) == (
        profile.levels,
        profile.level_target,
    )


def profile_normal(profile: ValueProfile) -> tuple[tuple[int, ...], int]:
    """Return the primitive ``(h, z)`` with ``sum(h) = 0`` and ``h . omega = z``.

    ``h = d * tau - sum(tau) * 1`` and ``z = -n * sum(tau)``, so that
    ``h . omega - z = d * (tau . omega)`` for every weight, which makes the
    incidence sets of the two descriptions identical by construction.
    """
    tau = profile.tau
    d = len(tau)
    total = sum(tau)
    h = tuple(d * value - total for value in tau)
    z = -profile.particle_number * total
    divisor = _content(h + (z,))
    if divisor == 0:
        raise ValueError("degenerate profile has no hyperplane")
    return tuple(value // divisor for value in h), z // divisor


def _level_sets_of_spread(r: int, max_spread: int):
    """Increasing integer level sets ``0 = l_1 < ... < l_r <= max_spread``, gcd 1."""
    if r == 1:
        return
    for spread in range(r - 1, max_spread + 1):
        interiors = [()] if r == 2 else itertools.combinations(range(1, spread), r - 2)
        for interior in interiors:
            levels = (0,) + tuple(interior) + (spread,)
            if _content(levels) != 1:
                continue
            yield levels


def _compositions(total: int, parts: int):
    for cuts in itertools.combinations(range(1, total), parts - 1):
        bounds = (0,) + cuts + (total,)
        yield tuple(bounds[i + 1] - bounds[i] for i in range(parts))


@dataclass(frozen=True)
class _LevelProfile:
    levels: tuple[int, ...]
    target: int
    patterns: tuple[tuple[int, ...], ...]
    values: tuple[int, ...]


def _level_profiles(n: int, r: int, max_spread: int) -> tuple[_LevelProfile, ...]:
    """Stage one: level profiles whose UNCAPPED pattern set already has rank r-1.

    This stage does not depend on ``d``, so it is computed once per ``r`` and
    reused across ranks.  Capping can only remove patterns and rank is monotone,
    so a level profile failing here fails for every composition.
    """
    uncapped = occupancy_patterns(n, r)
    out: list[_LevelProfile] = []
    for levels in _level_sets_of_spread(r, max_spread):
        buckets: dict[int, list[tuple[int, ...]]] = {}
        for pattern in uncapped:
            total = sum(
                count * level for count, level in zip(pattern, levels, strict=True)
            )
            buckets.setdefault(total, []).append(pattern)
        for target, patterns in buckets.items():
            if len(patterns) < r - 1:
                continue
            if _rank_mod(patterns, r) != r - 1:
                continue
            values = _primitive(tuple(target - n * level for level in levels))
            if len(set(values)) != r:
                continue
            out.append(_LevelProfile(levels, target, tuple(patterns), values))
    return tuple(out)


def enumerate_admissible_profiles(
    n: int, d: int, max_spread: int, *, level_cache: dict | None = None
) -> tuple[ValueProfile, ...]:
    """Every admissible value profile at ``(n, d)`` of level spread at most ``max_spread``.

    Completeness is exact WITHIN the spread bound and unproved outside it.  The
    return value is sorted, so two runs are byte comparable.
    """
    if type(n) is not int or type(d) is not int or not (0 < n < d):
        raise ValueError("require integer parameters with 0 < n < d")
    if type(max_spread) is not int or max_spread < 1:
        raise ValueError("max_spread must be a positive integer")
    found: dict[tuple[tuple[int, ...], tuple[int, ...]], ValueProfile] = {}
    for r in range(2, d + 1):
        if n**r >= _RANK_PRIME:
            raise ValueError("rank prime does not exceed the Hadamard bound")
        key = (n, r, max_spread)
        if level_cache is not None and key in level_cache:
            stage_one = level_cache[key]
        else:
            stage_one = _level_profiles(n, r, max_spread)
            if level_cache is not None:
                level_cache[key] = stage_one
        if not stage_one:
            continue
        compositions = tuple(_compositions(d, r))
        for level_profile in stage_one:
            # Admissibility depends on m only through min(m_a, n + 1): patterns
            # cap at n, and mobility needs one spare mode beyond some pattern.
            memo: dict[tuple[int, ...], bool] = {}
            for multiplicities in compositions:
                capped = tuple(min(value, n + 1) for value in multiplicities)
                verdict = memo.get(capped)
                if verdict is None:
                    verdict = profile_is_admissible(n, capped, level_profile.values)
                    memo[capped] = verdict
                if not verdict:
                    continue
                found[(multiplicities, level_profile.values)] = ValueProfile(
                    particle_number=n,
                    multiplicities=multiplicities,
                    values=level_profile.values,
                    levels=level_profile.levels,
                    level_target=level_profile.target,
                )
    return tuple(sorted(found.values()))


def orbit_normals(
    n: int, d: int, max_spread: int, *, level_cache: dict | None = None
) -> tuple[tuple[tuple[int, ...], int], ...]:
    """One exact primitive ``(h, z)`` per admissible hyperplane orbit."""
    return tuple(
        profile_normal(profile)
        for profile in enumerate_admissible_profiles(
            n, d, max_spread, level_cache=level_cache
        )
    )


def verify_profile_against_weights(n: int, profile: ValueProfile) -> bool:
    """Check Theorem A against the materialized weight rank for one profile."""
    return weight_affine_rank(n, profile.tau) == profile.rank - 1


def incidence_agrees(n: int, profile: ValueProfile) -> bool:
    """Check that ``(h, z)`` cuts out exactly the weights that ``tau`` does."""
    h, z = profile_normal(profile)
    tau = profile.tau
    for weight in exterior_weights(n, len(tau)):
        on_h = sum(left * right for left, right in zip(h, weight, strict=True)) == z
        on_tau = sum(left * right for left, right in zip(tau, weight, strict=True)) == 0
        if on_h != on_tau:
            return False
    return True


def spread_saturation(
    n: int, d: int, spreads, *, level_cache: dict | None = None
) -> tuple[dict[str, object], ...]:
    """Run the enumerator at increasing spread bounds and report the plateau.

    The protocol is the only evidence this module offers that a spread bound is
    generous enough, and it is evidence rather than proof.  It is calibrated:
    at ranks 6 through 10 the plateau value equals the independently known orbit
    count, and the realized maximum spread stays strictly below the bound.
    """
    cache = {} if level_cache is None else level_cache
    out = []
    for max_spread in sorted(spreads):
        profiles = enumerate_admissible_profiles(n, d, max_spread, level_cache=cache)
        realized = max((profile.spread for profile in profiles), default=0)
        out.append(
            {
                "max_spread": max_spread,
                "profiles": len(profiles),
                "realized_max_spread": realized,
                "bound_is_slack": realized < max_spread,
            }
        )
    return tuple(out)


def enumerate_profile_survivor_taus(
    n: int, d: int, max_spread: int, *, level_cache: dict | None = None
):
    """Emit exactly the taus that pass the Ressayre trace test, orbit natively.

    This is ``orbit_native.enumerate_orbit_survivor_taus`` with its first stage
    replaced: orbit representatives come from the profile enumeration rather
    than from the closed-flat lattice.  Everything after that is the existing
    machinery, imported rather than reimplemented, so the two routes can be
    compared row for row.

    Returns ``(taus, counts)``.  ``counts`` carries the conservation data:
    survivors predicted by the Gaussian multinomial and survivors produced.  A
    disagreement raises, because a generator emitting the right number of wrong
    arrangements would silently replace the candidate list.

    One profile is one ORIENTED orbit, so orientations are not re-derived here.
    A self-paired orbit, meaning one whose normal multiset is closed under
    negation, appears once in the profile list and is therefore emitted once,
    which is the same accounting the lattice route does by hand.
    """
    from .bdr import primitive_tau
    from .orbit_canonical import generate_trace_survivors, trace_survivor_count
    from .ressayre import _level_sets

    profiles = enumerate_admissible_profiles(n, d, max_spread, level_cache=level_cache)
    weights = exterior_weights(n, d)
    taus: list[tuple[int, ...]] = []
    predicted = 0
    for profile in profiles:
        h, z = profile_normal(profile)
        below = len(_level_sets(weights, h, z)[1])
        expected = trace_survivor_count(h, below)
        predicted += expected
        produced = 0
        for arrangement in generate_trace_survivors(h, below):
            taus.append(primitive_tau(tuple(z - n * value for value in arrangement)))
            produced += 1
        if produced != expected:
            raise AssertionError(
                "generated survivor count disagrees with the Gaussian "
                "multinomial, which means one of the two routes is wrong"
            )
    if len(set(taus)) != len(taus):
        raise AssertionError(
            "profile-native generation produced duplicate taus, which means an "
            "orientation was counted twice"
        )
    counts = {
        "oriented_orbits": len(profiles),
        "predicted_survivors": predicted,
        "produced_survivors": len(taus),
        "max_spread": max_spread,
    }
    return tuple(taus), counts
