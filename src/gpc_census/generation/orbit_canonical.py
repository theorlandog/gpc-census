"""Weyl-orbit canonicalization of closed flats of the exterior weight set.

WHY THIS IS THE ONE PERMITTED PRUNING. The moment polytope is Weyl invariant
and the ordered slice is a fundamental domain for the ``S_d`` action, so two
admissible hyperplanes in the same ``S_d`` orbit are facets together or not at
all. Enumerating orbit representatives therefore discards nothing. Every other
pruning proposed for this search has been an empirical pattern, which the
standing rule forbids until a theorem makes it lossless; this one has the
theorem.

WHY IT MATTERS. The raw flat lattice explodes while the orbit count barely
moves: 362 hyperplanes in 8 orbits at rank 6, 5,341 in 19 at rank 7. An
enumerator that works orbit-wise never materialises the lattice.

THE REPRESENTATION. A flat is carried as a bitmask over the exterior weights,
one bit per ``n``-subset of ``d`` modes, exactly as
``generation.admissible_flats`` carries it. ``S_d`` permutes the modes, which
permutes the weights, which permutes the bits.

CANONICALIZATION. Brute force over ``d!`` is hopeless past rank 9 (``9! =
362,880``, ``11! = 4.0e7``), so modes are first split by an iterated incidence
signature and only permutations respecting that ordered partition are tried.
The signature is an invariant of the flat, so the refinement can never separate
two modes that some symmetry identifies; it only ever narrows the search.
Correctness therefore does not depend on the refinement being strong, only on
it being invariant, and the tests check the canonical form against orbit counts
computed independently by union-find.
"""

from __future__ import annotations

import itertools
from functools import lru_cache


@lru_cache(maxsize=None)
def weight_subsets(n: int, d: int) -> tuple[tuple[int, ...], ...]:
    """The exterior weights as ``n``-subsets, in the enumerator's bit order."""
    return tuple(itertools.combinations(range(d), n))


@lru_cache(maxsize=None)
def _subset_index(n: int, d: int) -> dict[tuple[int, ...], int]:
    return {subset: i for i, subset in enumerate(weight_subsets(n, d))}


def mask_to_subsets(mask: int, n: int, d: int) -> tuple[tuple[int, ...], ...]:
    """The weights present in a flat mask."""
    subsets = weight_subsets(n, d)
    return tuple(subsets[i] for i in range(len(subsets)) if mask >> i & 1)


def _signatures(present: tuple[tuple[int, ...], ...], d: int,
                rounds: int = 3) -> tuple[tuple, ...]:
    """An iterated, permutation-invariant signature per mode.

    Round zero is the degree of the mode in the flat. Each later round replaces
    a mode's signature by its own signature together with the sorted multiset
    of its co-members' signatures, which is the standard colour refinement and
    is invariant under relabelling by construction.
    """
    signature = [0] * d
    for mode in range(d):
        signature[mode] = sum(1 for subset in present if mode in subset)
    for _ in range(rounds):
        refined = []
        for mode in range(d):
            partners = []
            for subset in present:
                if mode in subset:
                    partners.append(tuple(sorted(signature[other]
                                                 for other in subset
                                                 if other != mode)))
            refined.append((signature[mode], tuple(sorted(partners))))
        # compress to small integers so the next round stays cheap
        order = {value: i for i, value in enumerate(sorted(set(refined)))}
        signature = [order[value] for value in refined]
    return tuple(signature)


def _cells(signature: tuple[int, ...], d: int) -> list[list[int]]:
    """Modes grouped by signature, cells ordered by signature value."""
    buckets: dict[int, list[int]] = {}
    for mode in range(d):
        buckets.setdefault(signature[mode], []).append(mode)
    return [buckets[key] for key in sorted(buckets)]


def _apply(present: tuple[tuple[int, ...], ...], mapping: list[int],
           n: int, d: int) -> int:
    index = _subset_index(n, d)
    mask = 0
    for subset in present:
        image = tuple(sorted(mapping[mode] for mode in subset))
        mask |= 1 << index[image]
    return mask


def canonical_form(mask: int, n: int, d: int,
                   with_automorphisms: bool = False):
    """Lexicographically minimal image of a flat mask under ``S_d``.

    Returns the canonical mask, or ``(mask, automorphism_count)`` when
    ``with_automorphisms`` is set. The automorphism count is the number of
    permutations fixing the flat, which gives the orbit size as
    ``d! / automorphisms`` without enumerating the orbit.
    """
    present = mask_to_subsets(mask, n, d)
    if not present:
        return (0, _factorial(d)) if with_automorphisms else 0

    cells = _cells(_signatures(present, d), d)
    best = None
    automorphisms = 0
    for choice in itertools.product(*(itertools.permutations(cell)
                                      for cell in cells)):
        mapping = [0] * d
        target = 0
        for cell, ordering in zip(cells, choice, strict=True):
            for mode in ordering:
                mapping[mode] = target
                target += 1
        image = _apply(present, mapping, n, d)
        if best is None or image < best:
            best = image
            automorphisms = 1
        elif image == best:
            automorphisms += 1
    if with_automorphisms:
        # permutations outside the refinement cannot fix the flat, because the
        # signature is an invariant, so counting inside the cells is complete
        return best, automorphisms
    return best


@lru_cache(maxsize=None)
def _factorial(k: int) -> int:
    result = 1
    for i in range(2, k + 1):
        result *= i
    return result


def orbit_size(mask: int, n: int, d: int) -> int:
    """Size of the ``S_d`` orbit of a flat, by orbit-stabilizer."""
    _, automorphisms = canonical_form(mask, n, d, with_automorphisms=True)
    return _factorial(d) // automorphisms


def canonical_representatives(masks, n: int, d: int) -> dict[int, int]:
    """Map canonical form -> one representative mask, for a set of flats."""
    out: dict[int, int] = {}
    for mask in masks:
        key = canonical_form(mask, n, d)
        if key not in out:
            out[key] = mask
    return out


def enumerate_orbits_by_augmentation(n: int, d: int, verbose: bool = False):
    """Enumerate ``S_d`` orbit representatives of closed flats, level by level.

    THE ALGORITHM. Hold only orbit representatives at each rank. Extend each
    representative by the repository's own one-rank closure step, canonicalize
    every child, and keep one per canonical form. The full lattice is never
    materialised.

    WHY IT LOSES NOTHING. Let ``O'`` be an orbit of rank ``k+1`` flats and pick
    ``F'`` in it. ``F'`` contains some closed rank-``k`` flat ``F``, and ``F``
    lies in some orbit whose stored representative is ``R = sigma F``. Then
    ``sigma F'`` contains ``R`` and still lies in ``O'``, so extending ``R``
    produces a member of ``O'``. Every orbit at level ``k+1`` is therefore
    reached from some representative at level ``k``.

    Returns ``(orbit_counts_by_rank, representatives_by_rank)`` where the
    counts are the number of orbits at each rank, starting with the empty flat.
    """
    from .admissible_flats import (
        _extend_one_rank,
        _inverse_table,
        rank_preserving_modulus,
    )
    from .exterior import exterior_weights

    modulus, _ = rank_preserving_modulus(n, d)
    weights = exterior_weights(n, d)
    points = tuple(tuple(weight) + (1,) for weight in weights)
    inverse_table = _inverse_table(modulus)
    points_array = None
    if inverse_table is not None:
        import numpy as np

        points_array = np.asarray(points, dtype=np.int64)

    states: dict[int, tuple[int, ...]] = {0: ()}
    counts = [1]
    levels = [dict(states)]
    for rank in range(1, d):
        children = _extend_one_rank(
            states, points, modulus, points_array, inverse_table)
        # keep one ACTUAL mask per canonical class, not the canonical mask,
        # because the extension step needs a genuine flat with its basis
        representatives: dict[int, tuple[int, ...]] = {}
        seen: set[int] = set()
        for mask, witness in children.items():
            key = canonical_form(mask, n, d)
            if key not in seen:
                seen.add(key)
                representatives[mask] = witness
        states = representatives
        counts.append(len(states))
        levels.append(dict(states))
        if verbose:
            print(f"  rank {rank}: {len(children)} children -> "
                  f"{len(states)} orbits", flush=True)
    return counts, levels


__all__ = [
    "canonical_form",
    "canonical_representatives",
    "enumerate_orbits_by_augmentation",
    "mask_to_subsets",
    "orbit_size",
    "weight_subsets",
]
