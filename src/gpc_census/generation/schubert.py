"""Independent flag-Schubert cross-checks for generated constraints.

This is intentionally not yet a general Schubert engine.  It implements the
closed-form divisor calculation for the Borland-Dennis GPC and exposes every
convention needed to turn it into a regression gate for the future engine.
"""

from __future__ import annotations

import itertools
from collections import Counter

import sympy as sp

from .model import IntegralConstraint, SchubertCertificate


def divided_difference(poly: sp.Expr, variables: tuple[sp.Symbol, ...], index: int) -> sp.Expr:
    """Apply the type-A simple divided difference ``partial_index``."""
    if not 0 <= index < len(variables) - 1:
        raise ValueError("divided-difference index out of range")
    left, right = variables[index], variables[index + 1]
    swapped = poly.xreplace({left: right, right: left})
    return sp.cancel((poly - swapped) / (left - right))


def _apply_word(
    poly: sp.Expr, variables: tuple[sp.Symbol, ...], word: tuple[int, ...]
) -> sp.Expr:
    for index in word:
        poly = sp.expand(divided_difference(poly, variables, index))
    return poly


def _cyclic_permutation(length: int, size: int) -> tuple[int, ...]:
    return tuple(range(2, length + 2)) + (1,) + tuple(range(length + 2, size + 1))


def borland_dennis_schubert_certificates() -> tuple[SchubertCertificate, ...]:
    """Return coefficient-one certificates for all four raw rank-6 rows."""
    variables = sp.symbols("x0:6")
    pair_test = (1, 1, 0, 0, 0, 0)
    pair_polynomial = sp.prod(
        variables[0] + variables[1] + variables[index] for index in range(2, 6)
    )
    pair_induced = tuple(
        sorted(
            Counter(
                sum(pair_test[index] for index in determinant)
                for determinant in itertools.combinations(range(6), 3)
            ).items(),
            reverse=True,
        )
    )
    if pair_induced != ((2, 4), (1, 12), (0, 4)):
        raise AssertionError("unexpected pair-test exterior spectrum")

    pair_data = (
        ((1, 6, 2, 3, 4, 5), (1, 2, 3, 4), IntegralConstraint((1, 0, 0, 0, 0, 1), 1)),
        ((2, 5, 1, 3, 4, 6), (1, 2, 3, 0), IntegralConstraint((0, 1, 0, 0, 1, 0), 1)),
        ((3, 4, 1, 2, 5, 6), (1, 2, 0, 1), IntegralConstraint((0, 0, 1, 1, 0, 0), 1)),
    )
    certificates: list[SchubertCertificate] = []
    for permutation, word, constraint in pair_data:
        coefficient = _apply_word(pair_polynomial, variables, word)
        if coefficient != 1:
            raise AssertionError("Borland-Dennis pair coefficient is not one")
        certificates.append(
            SchubertCertificate(
                test_spectrum=pair_test,
                source_permutation=permutation,
                target_permutation=_cyclic_permutation(4, 20),
                induced_spectrum=pair_induced,
                divided_difference_word=word,
                schubert_coefficient=1,
                constraint=constraint,
            )
        )

    gpc_test = (1, 1, 1, 0, 0, 0)
    gpc_polynomial = variables[0] + variables[1] + variables[2]
    gpc_coefficient = _apply_word(gpc_polynomial, variables, (2,))
    gpc_induced = tuple(
        sorted(
            Counter(
                sum(gpc_test[index] for index in determinant)
                for determinant in itertools.combinations(range(6), 3)
            ).items(),
            reverse=True,
        )
    )
    if gpc_coefficient != 1 or gpc_induced != ((3, 1), (2, 9), (1, 9), (0, 1)):
        raise AssertionError("unexpected nontrivial rank-6 Schubert calculation")
    certificates.append(
        SchubertCertificate(
            test_spectrum=gpc_test,
            source_permutation=(1, 2, 4, 3, 5, 6),
            target_permutation=_cyclic_permutation(1, 20),
            induced_spectrum=gpc_induced,
            divided_difference_word=(2,),
            schubert_coefficient=1,
            constraint=IntegralConstraint((1, 1, 0, 1, 0, 0), 2),
        )
    )
    return tuple(certificates)


def verify_borland_dennis_schubert_certificates(
    certificates: tuple[SchubertCertificate, ...],
) -> bool:
    """Replay all four independent rank-6 Schubert calculations."""
    try:
        return certificates == borland_dennis_schubert_certificates()
    except (AssertionError, TypeError, ValueError):
        return False
