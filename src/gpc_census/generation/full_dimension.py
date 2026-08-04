"""Exact full-dimension certificates for ``wedge^3 C^d``.

The Bulois-Denis-Ressayre moment-cone theorem assumes that the moment cone
has nonempty interior.  For a unitary representation this follows when one
can exhibit a vector with zero infinitesimal compact stabilizer.  This module
provides such an exhibit, with an exact finite-field certificate, for every
rank ``7 <= d <= 13`` used by the generator program.

Write an element of ``u(d)`` uniquely as ``A + i B``, where ``A`` is real
skew-symmetric and ``B`` is real symmetric.  The vector constructed here has
real integer coefficients.  Consequently

``(A + i B) psi = 0`` if and only if ``A psi = 0`` and ``B psi = 0``.

We form the two integer action matrices

``so(d) -> wedge^3 R^d, A |-> A psi`` and
``Sym(d) -> wedge^3 R^d, B |-> B psi``.

Full column rank modulo a prime exhibits a maximal integer minor that is
nonzero modulo that prime.  The same minor is therefore nonzero over the
rationals and the reals.  Full column rank of both matrices proves that the
compact Lie stabilizer of ``psi`` is zero.  The certificate payload records
the coefficient rule, prime, ranks, pivot-row minors, and their nonzero
determinants, so verification uses integer arithmetic only.

Rank six is deliberately rejected.  Independently of the chosen real vector
in this construction, ``dim Sym(6) = 21`` is larger than
``dim wedge^3 R^6 = 20``, so the symmetric action map cannot be injective.
The witness below has modular ranks 15 and 19, consistent with the separately
known lower-dimensional Borland-Dennis cone.

For ``m = binomial(d, 3)`` and ``n = dim u(d) = d^2``, dense matrix
construction takes ``O(m n) = O(d^5)`` integer operations and storage.
Gaussian elimination takes ``O(m n^2) = O(d^7)`` field operations; the two
minor determinants cost only ``O(d^6)`` more.  These bounds are deliberately
simple and deterministic because this certificate is computed once per rank.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import Final, Sequence


CERTIFICATE_PRIME: Final = 2_147_483_647
"""Prime used for all modular rank and minor computations."""

LCG_MODULUS: Final = 2**31
LCG_MULTIPLIER: Final = 1_103_515_245
LCG_INCREMENT: Final = 12_345
COEFFICIENT_MODULUS: Final = 100

SEED_RULE: Final = (
    "x_0=d; x_(k+1)=(1103515245*x_k+12345) mod 2^31; "
    "c_k=1+(x_(k+1) mod 100), with 3-subsets in lexicographic order"
)


@dataclass(frozen=True)
class FullDimensionCertificate:
    """Replayable zero-stabilizer certificate for ``wedge^3 C^d``.

    The two ``rank`` fields are ranks over ``F_prime``.  When they equal their
    respective column counts, the recorded nonzero minors prove the same full
    column ranks over the rationals and reals.
    """

    d: int
    prime: int
    seed_rule: str
    coefficients: tuple[int, ...]
    skew_rank: int
    symmetric_rank: int
    expected_skew_rank: int
    expected_symmetric_rank: int
    skew_pivot_rows: tuple[int, ...]
    symmetric_pivot_rows: tuple[int, ...]
    skew_minor_determinant_mod_prime: int | None
    symmetric_minor_determinant_mod_prime: int | None

    @property
    def is_full_dimensional(self) -> bool:
        """Whether the payload proves a zero compact Lie stabilizer."""
        return (
            self.skew_rank == self.expected_skew_rank
            and self.symmetric_rank == self.expected_symmetric_rank
            and self.skew_minor_determinant_mod_prime not in (None, 0)
            and self.symmetric_minor_determinant_mod_prime not in (None, 0)
        )

    def as_dict(self) -> dict[str, object]:
        """Return the complete deterministic replay payload."""
        return {
            "rank": self.d,
            "prime": self.prime,
            "seed_rule": self.seed_rule,
            "coefficients": list(self.coefficients),
            "skew_rank_mod_prime": self.skew_rank,
            "symmetric_rank_mod_prime": self.symmetric_rank,
            "expected_skew_rank": self.expected_skew_rank,
            "expected_symmetric_rank": self.expected_symmetric_rank,
            "skew_pivot_rows": list(self.skew_pivot_rows),
            "symmetric_pivot_rows": list(self.symmetric_pivot_rows),
            "skew_minor_determinant_mod_prime": self.skew_minor_determinant_mod_prime,
            "symmetric_minor_determinant_mod_prime": (
                self.symmetric_minor_determinant_mod_prime
            ),
            "is_full_dimensional": self.is_full_dimensional,
        }


def _coefficients(d: int) -> tuple[int, ...]:
    state = d
    result: list[int] = []
    for _ in range(d * (d - 1) * (d - 2) // 6):
        state = (LCG_MULTIPLIER * state + LCG_INCREMENT) % LCG_MODULUS
        result.append(1 + state % COEFFICIENT_MODULUS)
    return tuple(result)


def _matrix_unit_image(
    occupied: tuple[int, int, int], row: int, column: int
) -> tuple[tuple[int, int, int], int] | None:
    """Return ``E_(row,column) e_occupied`` in the ordered wedge basis."""
    if column not in occupied:
        return None
    if row == column:
        return occupied, 1
    if row in occupied:
        return None

    source_position = occupied.index(column)
    remaining = list(occupied)
    remaining.pop(source_position)
    target_position = sum(value < row for value in remaining)
    sign = -1 if (source_position + target_position) % 2 else 1
    remaining.insert(target_position, row)
    return (remaining[0], remaining[1], remaining[2]), sign


def _action_matrices(
    d: int, coefficients: Sequence[int]
) -> tuple[list[list[int]], list[list[int]]]:
    basis = tuple(combinations(range(d), 3))
    if len(coefficients) != len(basis):
        raise ValueError("coefficient count must equal binomial(d, 3)")
    basis_index = {occupied: index for index, occupied in enumerate(basis)}

    skew_column_count = d * (d - 1) // 2
    symmetric_column_count = d * (d + 1) // 2
    skew = [[0] * skew_column_count for _ in basis]
    symmetric = [[0] * symmetric_column_count for _ in basis]

    for diagonal in range(d):
        for source_index, occupied in enumerate(basis):
            if diagonal in occupied:
                symmetric[source_index][diagonal] += coefficients[source_index]

    skew_column = 0
    symmetric_column = d
    for a in range(d):
        for b in range(a + 1, d):
            for source_index, occupied in enumerate(basis):
                coefficient = coefficients[source_index]

                image = _matrix_unit_image(occupied, a, b)
                if image is not None:
                    target, sign = image
                    target_index = basis_index[target]
                    skew[target_index][skew_column] += sign * coefficient
                    symmetric[target_index][symmetric_column] += sign * coefficient

                image = _matrix_unit_image(occupied, b, a)
                if image is not None:
                    target, sign = image
                    target_index = basis_index[target]
                    skew[target_index][skew_column] -= sign * coefficient
                    symmetric[target_index][symmetric_column] += sign * coefficient

            skew_column += 1
            symmetric_column += 1

    return skew, symmetric


def _rank_and_pivot_rows_mod(
    matrix: Sequence[Sequence[int]], prime: int
) -> tuple[int, tuple[int, ...]]:
    if not matrix:
        return 0, ()
    column_count = len(matrix[0])
    if any(len(row) != column_count for row in matrix):
        raise ValueError("matrix rows must have equal length")

    work = [[entry % prime for entry in row] for row in matrix]
    row_ids = list(range(len(work)))
    rank = 0
    pivot_rows: list[int] = []
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        row_ids[rank], row_ids[pivot] = row_ids[pivot], row_ids[rank]
        pivot_rows.append(row_ids[rank])

        inverse = pow(work[rank][column], prime - 2, prime)
        work[rank] = [(entry * inverse) % prime for entry in work[rank]]
        for row in range(rank + 1, len(work)):
            factor = work[row][column]
            if factor:
                work[row] = [
                    (entry - factor * pivot_entry) % prime
                    for entry, pivot_entry in zip(work[row], work[rank], strict=True)
                ]
        rank += 1
        if rank == column_count:
            break
    return rank, tuple(pivot_rows)


def _determinant_mod(matrix: Sequence[Sequence[int]], prime: int) -> int:
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, prime - 2, prime)
        for row in range(column + 1, size):
            factor = work[row][column] * inverse % prime
            if factor:
                for entry_column in range(column, size):
                    work[row][entry_column] = (
                        work[row][entry_column] - factor * work[column][entry_column]
                    ) % prime
    return determinant % prime


def _full_column_minor(
    matrix: Sequence[Sequence[int]], rank: int, pivot_rows: Sequence[int], prime: int
) -> int | None:
    if not matrix:
        return 1
    column_count = len(matrix[0])
    if rank != column_count:
        return None
    minor = [matrix[row] for row in pivot_rows]
    return _determinant_mod(minor, prime)


@lru_cache(maxsize=None)
def full_dimension_certificate(d: int) -> FullDimensionCertificate:
    """Construct the deterministic exact certificate for rank ``d``.

    The returned payload can also represent failure.  In particular, rank
    six has symmetric rank 19 modulo the certificate prime against the
    required 21 and therefore has ``is_full_dimensional == False``.
    """
    if d < 3:
        raise ValueError("wedge^3 C^d requires d >= 3")

    coefficients = _coefficients(d)
    skew, symmetric = _action_matrices(d, coefficients)
    skew_rank, skew_pivots = _rank_and_pivot_rows_mod(skew, CERTIFICATE_PRIME)
    symmetric_rank, symmetric_pivots = _rank_and_pivot_rows_mod(symmetric, CERTIFICATE_PRIME)
    expected_skew_rank = d * (d - 1) // 2
    expected_symmetric_rank = d * (d + 1) // 2

    return FullDimensionCertificate(
        d=d,
        prime=CERTIFICATE_PRIME,
        seed_rule=SEED_RULE,
        coefficients=coefficients,
        skew_rank=skew_rank,
        symmetric_rank=symmetric_rank,
        expected_skew_rank=expected_skew_rank,
        expected_symmetric_rank=expected_symmetric_rank,
        skew_pivot_rows=skew_pivots,
        symmetric_pivot_rows=symmetric_pivots,
        skew_minor_determinant_mod_prime=_full_column_minor(
            skew, skew_rank, skew_pivots, CERTIFICATE_PRIME
        ),
        symmetric_minor_determinant_mod_prime=_full_column_minor(
            symmetric, symmetric_rank, symmetric_pivots, CERTIFICATE_PRIME
        ),
    )


def require_full_dimension(d: int) -> FullDimensionCertificate:
    """Return a certificate or reject a rank not certified by this witness."""
    certificate = full_dimension_certificate(d)
    if not certificate.is_full_dimensional:
        raise ValueError(
            f"wedge^3 C^{d} is not full-dimensional under the deterministic witness: "
            f"ranks are ({certificate.skew_rank}, {certificate.symmetric_rank}) "
            f"but ({certificate.expected_skew_rank}, "
            f"{certificate.expected_symmetric_rank}) are required"
        )
    return certificate


def verify_full_dimension_certificate(certificate: FullDimensionCertificate) -> bool:
    """Replay a certificate from its public payload using exact arithmetic."""
    if certificate.prime != CERTIFICATE_PRIME or certificate.seed_rule != SEED_RULE:
        return False
    if certificate.coefficients != _coefficients(certificate.d):
        return False
    expected = full_dimension_certificate(certificate.d)
    return certificate == expected


__all__ = [
    "CERTIFICATE_PRIME",
    "FullDimensionCertificate",
    "SEED_RULE",
    "full_dimension_certificate",
    "require_full_dimension",
    "verify_full_dimension_certificate",
]
