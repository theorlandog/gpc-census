"""Constructive inversion of a signed Chow shadow.

``docs/signed_chow_projection.md`` proves that the pair ``(c_+, c_-)``
determines an admissible normal (T1 to T3) and records, as the one thing it
does not supply, an ALGORITHM: "recovering ``B_+`` from a bare ``c_+`` is the
Chow-parameters inversion problem and is not solved here". This module supplies
it, using two further theorems proved in
``docs/cocircuit_chow_decoder.md``.

T4 (block symmetry). If ``c_i == c_j`` then the transposition of orbitals ``i``
and ``j`` fixes ``c_+``, so the permuted family has the same shadow, so T1
forces it to be ``B_+`` again. Hence ``B_+`` is invariant under the product of
symmetric groups on the equal-degree blocks of ``c_+``, and membership depends
only on a determinant's OCCUPANCY PROFILE ``k = (k_1, ..., k_r)`` across those
blocks. The unknown is a subset of the profiles, not a subset of the slice.

T5 (block-constant separator, and its order). Averaging any separator over that
product group leaves the strict-positive family unchanged and yields a
separator constant on each block, whose block levels are ordered by ``c_+``
through T2. Writing the blocks in decreasing ``c_+`` order, Abel summation
makes the separator value monotone in prefix-sum (Gale) dominance, so the
selected profiles form an UPPER IDEAL in that order.

The decoder solves the resulting exact 0/1 system

    sum_k q_a(k) x_k = c_a  for every block a,
    q_a(k) = C(m_a - 1, k_a - 1) * prod_{b != a} C(m_b, k_b),

with upper-ideal propagation and integer feasibility bounds. Every quantity is
an integer; nothing here is numerical. T1 makes the search an exact decoder:
any profile set satisfying the block equations induces a family whose shadow is
``c_+``, and T1 admits only one such family. Worst-case polynomial complexity is
NOT claimed, and the measured search sizes in
``results/data/cocircuit_chow_decoder.json`` are measurements rather than a
bound.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from fractions import Fraction

DEFAULT_MAX_STATES = 2_000_000


class DecoderBudgetError(RuntimeError):
    """The profile search exceeded its state budget."""


class ShadowUniquenessError(RuntimeError):
    """Two profile sets reproduced one shadow, which T1 forbids."""


@dataclass(frozen=True)
class ShadowDecoding:
    """One side of a signed shadow, decoded back into its family.

    ``particle_hole_layer`` is the particle number the profile search actually
    ran in. It differs from the caller's ``n`` exactly when the complementary
    layer was smaller, and it is recorded because ``profiles``,
    ``profile_variables`` and ``search_states`` all describe THAT layer, so a
    benchmark that ignores it would compare two different problems.
    """

    family: frozenset[tuple[int, ...]]
    profiles: tuple[tuple[int, ...], ...]
    blocks: tuple[tuple[int, ...], ...]
    block_degrees: tuple[int, ...]
    profile_variables: int
    search_states: int
    particle_hole_layer: int


@dataclass(frozen=True)
class SignedDecoding:
    """Both sides of a signed shadow, plus the normal they reconstruct."""

    positive: ShadowDecoding
    negative: ShadowDecoding
    zero_family: tuple[tuple[int, ...], ...]
    zero_span_rank: int
    tau: tuple[int, ...] | None


def primitive(vector) -> tuple[int, ...]:
    """Return the primitive integer representative with a positive lead."""
    values = [int(value) for value in vector]
    divisor = 0
    for value in values:
        divisor = math.gcd(divisor, abs(value))
    if divisor:
        values = [value // divisor for value in values]
    for value in values:
        if value:
            return tuple(values) if value > 0 else tuple(-v for v in values)
    return tuple(values)


def degree_blocks(shadow: tuple[int, ...]):
    """Split orbitals into equal-degree blocks, largest degree first.

    T5 fixes the order: the averaged separator is strictly decreasing across
    blocks in this order, which is what makes prefix-sum dominance the right
    order on profiles.
    """
    levels = sorted(set(shadow), reverse=True)
    blocks = tuple(
        tuple(index for index, value in enumerate(shadow) if value == level)
        for level in levels
    )
    return blocks, tuple(levels)


def profiles_with_caps(total: int, caps: tuple[int, ...]):
    """Every way to place ``total`` particles into capped blocks."""

    def walk(position: int, left: int, prefix: list[int]):
        if position == len(caps) - 1:
            if 0 <= left <= caps[position]:
                yield tuple([*prefix, left])
            return
        for taken in range(min(caps[position], left) + 1):
            yield from walk(position + 1, left - taken, [*prefix, taken])

    if not caps:
        return
    yield from walk(0, total, [])


def dominates(high: tuple[int, ...], low: tuple[int, ...]) -> bool:
    """Prefix-sum dominance: ``high`` moves particles to earlier blocks.

    Blocks are ordered by decreasing separator level, so by Abel summation
    ``sum_a k_a t_a`` is monotone in this order whenever ``t`` decreases. That
    is exactly why the selected profile set is an upper ideal.
    """
    running_high = running_low = 0
    for above, below in zip(high[:-1], low[:-1], strict=True):
        running_high += above
        running_low += below
        if running_high < running_low:
            return False
    return True


def profile_contribution(
    profile: tuple[int, ...], sizes: tuple[int, ...]
) -> tuple[int, ...]:
    """Degree that one profile orbit contributes to a vertex of each block."""
    contribution = []
    for position, (taken, size) in enumerate(zip(profile, sizes, strict=True)):
        if taken == 0:
            contribution.append(0)
            continue
        count = math.comb(size - 1, taken - 1)
        for other, (other_taken, other_size) in enumerate(
            zip(profile, sizes, strict=True)
        ):
            if other != position:
                count *= math.comb(other_size, other_taken)
        contribution.append(count)
    return tuple(contribution)


def _expand(blocks, profiles) -> frozenset[tuple[int, ...]]:
    """Turn a set of occupancy profiles into the family of determinants."""
    family: set[tuple[int, ...]] = set()
    for profile in profiles:
        choices = [
            tuple(itertools.combinations(block, taken))
            for block, taken in zip(blocks, profile, strict=True)
        ]
        for pieces in itertools.product(*choices):
            family.add(tuple(sorted(itertools.chain.from_iterable(pieces))))
    return frozenset(family)


def complement_shadow(shadow: tuple[int, ...], n: int):
    """Return the ``(d-n)``-layer shadow of the complement family, and its size.

    PARTICLE-HOLE REDUCTION. Let ``F`` be an ``n``-uniform family of size ``m``
    with shadow ``c``, and let ``F^c = {[d] \\ I : I in F}``, which is
    ``(d-n)``-uniform of the same size. Orbital ``i`` occurs in ``[d] \\ I``
    exactly when it does not occur in ``I``, so

        c'(F^c)_i = m - c(F)_i,        m = sum(c) / n.

    Complementation is a bijection between the two layers, so decoding ``c'``
    in the ``(d-n)`` layer and complementing the answer back is exact, not an
    approximation. It is worth doing when ``n > d - n``, because the profile
    count grows with the particle number.

    THE ORDER FLIPS, AND THAT IS WHY IT STILL WORKS. Since ``c' = m - c``, the
    equal-degree blocks are the same sets in the REVERSE order, which is what
    ``degree_blocks`` produces from ``c'`` on its own. Writing ``M_j`` and
    ``K_j`` for the block-size and profile prefix sums in the original order,
    the prefix sums of the complement profile in the reversed order are
    ``(M - M_j) - (n - K_j)``, which is increasing in ``K_j``. So dominance is
    preserved rather than reversed, and the complement family is an upper ideal
    in its own layer exactly when the original is one in its layer.
    """
    if sum(shadow) % n:
        return None, 0
    cardinality = sum(shadow) // n
    complement = tuple(cardinality - value for value in shadow)
    if any(value < 0 for value in complement):
        return None, cardinality
    return complement, cardinality


def decode_shadow(
    shadow: tuple[int, ...],
    n: int,
    *,
    max_states: int = DEFAULT_MAX_STATES,
    particle_hole: bool = True,
) -> ShadowDecoding | None:
    """Decode one positive Chow shadow into the unique family that induces it.

    Returns ``None`` when the exact profile equations have no solution, which
    is the honest answer for a vector that is not the shadow of any family.
    Raises ``DecoderBudgetError`` past ``max_states`` and
    ``ShadowUniquenessError`` if a second solution appears, which T1 forbids
    for a realizable shadow and which therefore acts as a self-check.

    ``particle_hole`` routes the search through the smaller layer when that is
    the complementary one; see :func:`complement_shadow`. Pass ``False`` to
    force the literal layer, which is what the guard test uses to check that
    the reduction changes the cost and not the answer.
    """
    d = len(shadow)
    if any(value < 0 for value in shadow):
        return None
    if sum(shadow) % n:
        return None

    if particle_hole and n > d - n:
        complement, _cardinality = complement_shadow(shadow, n)
        if complement is None:
            return None
        dual = decode_shadow(
            complement, d - n, max_states=max_states, particle_hole=False
        )
        if dual is None:
            return None
        universe = frozenset(range(d))
        return ShadowDecoding(
            family=frozenset(
                tuple(sorted(universe - frozenset(subset)))
                for subset in dual.family
            ),
            profiles=dual.profiles,
            blocks=dual.blocks,
            block_degrees=dual.block_degrees,
            profile_variables=dual.profile_variables,
            search_states=dual.search_states,
            particle_hole_layer=d - n,
        )

    if not any(shadow):
        return ShadowDecoding(
            family=frozenset(),
            profiles=(),
            blocks=(tuple(range(d)),),
            block_degrees=(0,),
            profile_variables=1,
            search_states=1,
            particle_hole_layer=n,
        )

    blocks, levels = degree_blocks(shadow)
    sizes = tuple(len(block) for block in blocks)
    profiles = tuple(profiles_with_caps(n, sizes))
    if not profiles:
        return None
    weight = tuple(profile_contribution(profile, sizes) for profile in profiles)

    # Upper-ideal closures, by T5: including a profile includes everything
    # dominating it, excluding one excludes everything it dominates.
    include_closure = tuple(
        tuple(j for j, other in enumerate(profiles) if dominates(other, profile))
        for profile in profiles
    )
    exclude_closure = tuple(
        tuple(j for j, other in enumerate(profiles) if dominates(profile, other))
        for profile in profiles
    )

    states = 0
    solutions: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()

    def propagate(status, index, value):
        status = list(status)
        stack = [(index, value)]
        while stack:
            position, decision = stack.pop()
            current = status[position]
            if current == decision:
                continue
            if current == -decision:
                return None
            status[position] = decision
            closure = (
                include_closure[position]
                if decision == 1
                else exclude_closure[position]
            )
            for other in closure:
                if status[other] == -decision:
                    return None
                if status[other] == 0:
                    stack.append((other, decision))
        return tuple(status)

    def bounds(status):
        attained = [0] * len(levels)
        reachable = [0] * len(levels)
        for index, decision in enumerate(status):
            if decision == 1:
                for block, value in enumerate(weight[index]):
                    attained[block] += value
                    reachable[block] += value
            elif decision == 0:
                for block, value in enumerate(weight[index]):
                    reachable[block] += value
        return attained, reachable

    def saturate(status):
        while True:
            attained, reachable = bounds(status)
            if any(
                attained[block] > levels[block] or reachable[block] < levels[block]
                for block in range(len(levels))
            ):
                return None
            forced = None
            for index, decision in enumerate(status):
                if decision != 0:
                    continue
                if any(
                    attained[block] + weight[index][block] > levels[block]
                    for block in range(len(levels))
                ):
                    forced = (index, -1)
                    break
                if any(
                    reachable[block] - weight[index][block] < levels[block]
                    for block in range(len(levels))
                ):
                    forced = (index, 1)
                    break
            if forced is None:
                return status
            status = propagate(status, *forced)
            if status is None:
                return None

    def branch_on(status):
        attained, _ = bounds(status)
        residual = [
            levels[block] - attained[block] for block in range(len(levels))
        ]
        return max(
            (index for index, decision in enumerate(status) if decision == 0),
            key=lambda index: (
                len(include_closure[index]) + len(exclude_closure[index]),
                sum(
                    min(weight[index][block], residual[block])
                    for block in range(len(levels))
                ),
                sum(weight[index]),
            ),
        )

    def visit(status):
        nonlocal states
        states += 1
        if states > max_states:
            raise DecoderBudgetError(
                f"profile search exceeded max_states={max_states}"
            )
        status = saturate(status)
        if status is None or status in seen:
            return
        seen.add(status)
        if all(decision != 0 for decision in status):
            attained, _ = bounds(status)
            if tuple(attained) == levels:
                solutions.append(status)
                if len(solutions) > 1:
                    raise ShadowUniquenessError(
                        "two profile sets induced one shadow, contradicting T1"
                    )
            return
        index = branch_on(status)
        for decision in (1, -1):
            child = propagate(status, index, decision)
            if child is not None:
                visit(child)

    visit(tuple(0 for _ in profiles))
    if not solutions:
        return None
    chosen = tuple(
        profile
        for profile, decision in zip(profiles, solutions[0], strict=True)
        if decision == 1
    )
    return ShadowDecoding(
        family=_expand(blocks, chosen),
        profiles=chosen,
        blocks=blocks,
        block_degrees=levels,
        profile_variables=len(profiles),
        search_states=states,
        particle_hole_layer=n,
    )


def threshold_family(tau, n: int, sign: int = 1) -> frozenset[tuple[int, ...]]:
    """The determinants on which ``sign * tau`` is strictly positive."""
    return frozenset(
        subset
        for subset in itertools.combinations(range(len(tau)), n)
        if sign * sum(tau[index] for index in subset) > 0
    )


def shadow_of(family, d: int) -> tuple[int, ...]:
    """The Chow shadow of a family: how often each orbital is occupied."""
    counts = [0] * d
    for subset in family:
        for index in subset:
            counts[index] += 1
    return tuple(counts)


def annihilator(zero_family, d: int):
    """Primitive generator of the annihilator of ``span(B_0)``, and its rank.

    Exact rational elimination, so a returned normal is the normal and not a
    fit. Returns ``(None, rank)`` when the span is not a hyperplane, which is
    the failure mode T3 is explicitly conditional on.
    """
    rows = [
        [Fraction(1 if index in set(subset) else 0) for index in range(d)]
        for subset in zero_family
    ]
    pivots: list[int] = []
    rank = 0
    for column in range(d):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column]), None
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [value / scale for value in rows[rank]]
        for row in range(len(rows)):
            if row != rank and rows[row][column]:
                factor = rows[row][column]
                rows[row] = [
                    value - factor * basis
                    for value, basis in zip(rows[row], rows[rank], strict=True)
                ]
        pivots.append(column)
        rank += 1
        if rank == len(rows):
            break
    free = [column for column in range(d) if column not in pivots]
    if len(free) != 1:
        return None, rank
    index = free[0]
    solution = [Fraction(0)] * d
    solution[index] = Fraction(1)
    for row, column in enumerate(pivots):
        solution[column] = -rows[row][index]
    denominator = 1
    for value in solution:
        denominator = math.lcm(denominator, value.denominator)
    return primitive(int(value * denominator) for value in solution), rank


def decode_signed_pair(
    positive_shadow: tuple[int, ...],
    negative_shadow: tuple[int, ...],
    n: int,
    *,
    max_states: int = DEFAULT_MAX_STATES,
    particle_hole: bool = True,
) -> SignedDecoding | None:
    """Decode both sides of a signed shadow and rebuild the primitive normal.

    The two sides are decoded independently, their complement in the slice is
    ``B_0`` by T3, and the exact annihilator of ``span(B_0)`` is the normal ray
    whenever that span is a hyperplane. The decoded positive family fixes the
    orientation, and both sides are re-derived from the recovered normal before
    it is returned, so nothing is trusted that was not recomputed.
    """
    d = len(positive_shadow)
    if len(negative_shadow) != d:
        raise ValueError("the two shadows must have the same length")
    positive = decode_shadow(
        positive_shadow, n, max_states=max_states, particle_hole=particle_hole
    )
    negative = decode_shadow(
        negative_shadow, n, max_states=max_states, particle_hole=particle_hole
    )
    if positive is None or negative is None:
        return None
    if positive.family & negative.family:
        return None
    slice_weights = frozenset(itertools.combinations(range(d), n))
    zero_family = tuple(
        sorted(slice_weights - positive.family - negative.family)
    )
    recovered, rank = annihilator(zero_family, d)
    if recovered is not None:
        if threshold_family(recovered, n, +1) != positive.family:
            recovered = tuple(-value for value in recovered)
        if (
            threshold_family(recovered, n, +1) != positive.family
            or threshold_family(recovered, n, -1) != negative.family
        ):
            recovered = None
    return SignedDecoding(
        positive=positive,
        negative=negative,
        zero_family=zero_family,
        zero_span_rank=rank,
        tau=recovered,
    )
