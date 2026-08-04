"""Weights, roots, and exact root actions for exterior representations."""

from __future__ import annotations

import itertools

from .model import Root, Weight


def exterior_weights(n: int, d: int) -> tuple[Weight, ...]:
    """Return the multiplicity-free determinant weights of ``wedge^n C^d``."""
    if not 0 < n < d:
        raise ValueError("exterior weight generation requires 0 < n < d")
    return tuple(
        tuple(int(index in occupied) for index in range(d))
        for occupied in itertools.combinations(range(d), n)
    )


def negative_roots(d: int) -> tuple[Root, ...]:
    """Return negative type-A roots as pairs ``(i,j)`` for ``e_j-e_i``.

    Here ``i < j`` and the root operator is ``E_(j,i)``.  This convention is
    paired with the dominant chamber ``lambda_0 >= ... >= lambda_(d-1)``.
    """
    return tuple((i, j) for i in range(d) for j in range(i + 1, d))


def exterior_root_action(weight: Weight, root: Root) -> tuple[Weight, int] | None:
    """Apply ``E_(j,i) = a_j^dagger a_i`` to a determinant weight.

    Return the target weight and its fermionic sign, or ``None`` when the
    action vanishes.
    """
    i, j = root
    if not 0 <= i < j < len(weight):
        raise ValueError("root must satisfy 0 <= i < j < d")
    if not weight[i] or weight[j]:
        return None

    occupied = [index for index, value in enumerate(weight) if value]
    annihilation_position = occupied.index(i)
    sign = -1 if annihilation_position % 2 else 1
    remaining = occupied[:annihilation_position] + occupied[annihilation_position + 1 :]
    creation_position = sum(index < j for index in remaining)
    if creation_position % 2:
        sign = -sign
    target_occupied = remaining[:creation_position] + [j] + remaining[creation_position:]
    target = tuple(int(index in target_occupied) for index in range(len(weight)))
    return target, sign
