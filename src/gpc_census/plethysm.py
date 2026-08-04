"""Exact plethysm primitives for fermionic moment-polytope inner points.

For ``m >= 1``, every irreducible ``GL(d)`` component indexed by ``lambda``
in ``Sym^m(wedge^n C^d)`` certifies the attainable spectrum ``lambda / m``.
The multiplicity is computed from ``h_m[e_n]`` in the power-sum basis and the
Murnaghan-Nakayama rule.  This module contains only exact arithmetic.
"""

from __future__ import annotations

from collections.abc import Iterator
from fractions import Fraction
from functools import lru_cache


PowerSum = dict[tuple[int, ...], Fraction]
Spectrum = tuple[Fraction, ...]


def partitions(n: int, maxpart: int | None = None) -> Iterator[tuple[int, ...]]:
    """Yield decreasing integer partitions of ``n``."""
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for part in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - part, part):
            yield (part,) + rest


def zmu(mu: tuple[int, ...]) -> int:
    """Return the standard symmetric-function factor ``z_mu``."""
    from collections import Counter

    out = 1
    for part, multiplicity in Counter(mu).items():
        out *= part**multiplicity
        for value in range(1, multiplicity + 1):
            out *= value
    return out


def pmul(left: PowerSum, right: PowerSum) -> PowerSum:
    """Multiply two power-sum expansions."""
    out: PowerSum = {}
    for left_key, left_value in left.items():
        for right_key, right_value in right.items():
            key = tuple(sorted(left_key + right_key, reverse=True))
            out[key] = out.get(key, Fraction(0)) + left_value * right_value
    return out


def padd(left: PowerSum, right: PowerSum, coefficient: Fraction = Fraction(1)) -> PowerSum:
    """Add ``coefficient * right`` to ``left`` in place."""
    for key, value in right.items():
        left[key] = left.get(key, Fraction(0)) + coefficient * value
    return left


@lru_cache(maxsize=None)
def e_in_p(n: int) -> tuple[tuple[tuple[int, ...], Fraction], ...]:
    """Expand the elementary symmetric function ``e_n`` in power sums."""
    if n == 0:
        return (((), Fraction(1)),)
    total: PowerSum = {}
    for k in range(1, n + 1):
        term = pmul({(k,): Fraction(1)}, dict(e_in_p(n - k)))
        padd(total, term, Fraction((-1) ** (k - 1), n))
    return tuple(sorted(total.items()))


def pleth_pk(k: int, expansion: tuple[tuple[tuple[int, ...], Fraction], ...]) -> PowerSum:
    """Apply the Adams operation ``p_j -> p_(j k)``."""
    return {
        tuple(sorted((part * k for part in key), reverse=True)): value
        for key, value in expansion
    }


def hm_en(m: int, n: int) -> PowerSum:
    """Expand the plethysm ``h_m[e_n]`` in the power-sum basis."""
    en = e_in_p(n)
    total: PowerSum = {}
    for mu in partitions(m):
        term: PowerSum = {(): Fraction(1)}
        for part in mu:
            term = pmul(term, pleth_pk(part, en))
        padd(total, term, Fraction(1, zmu(mu)))
    return total


@lru_cache(maxsize=None)
def mn(lam: tuple[int, ...], mu: tuple[int, ...]) -> int:
    """Return the irreducible character ``chi^lam(mu)`` exactly."""
    if not mu:
        return int(sum(lam) == 0)
    k, rest = mu[0], mu[1:]
    lam = tuple(value for value in lam if value > 0)
    nrows = len(lam)
    if nrows == 0:
        return 0
    beta = [lam[index] + nrows - 1 - index for index in range(nrows)]
    beta_set = set(beta)
    total = 0
    for value in beta:
        if value - k < 0 or value - k in beta_set:
            continue
        height = sum(1 for other in beta if value - k < other < value)
        next_beta = sorted(
            (other if other != value else value - k for other in beta), reverse=True
        )
        next_lam = tuple(
            part
            for part in (
                next_beta[index] - (nrows - 1 - index) for index in range(nrows)
            )
            if part > 0
        )
        total += (-1) ** height * mn(next_lam, rest)
    return total


def inner_spectra(
    n: int, d: int, max_degree: int, *, verbose: bool = False
) -> dict[Spectrum, tuple[int, tuple[int, ...], int]]:
    """Return all plethysm-certified attainable spectra through ``max_degree``.

    The value records the first symmetric power, highest weight, and positive
    multiplicity certifying each spectrum.  No convex-hull assertion is made.
    """
    points: dict[Spectrum, tuple[int, tuple[int, ...], int]] = {}
    for degree in range(1, max_degree + 1):
        expansion = list(hm_en(degree, n).items())
        for lam in partitions(n * degree, maxpart=degree):
            if len(lam) > d:
                continue
            multiplicity = sum(
                coefficient * mn(lam, tuple(sorted(mu, reverse=True)))
                for mu, coefficient in expansion
            )
            if multiplicity <= 0:
                continue
            spectrum = tuple(
                [Fraction(value, degree) for value in lam]
                + [Fraction(0)] * (d - len(lam))
            )
            points.setdefault(spectrum, (degree, lam, int(multiplicity)))
        if verbose:
            print(f"  m={degree}: cumulative points {len(points)}")
    return points


# Backward-compatible spelling used by the original Stage-0 script.
hm_eN = hm_en
