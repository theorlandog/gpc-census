from __future__ import annotations

from gpc_census.generation.admissible_flats import (
    enumerate_admissible_hyperplanes_by_flats,
    rank_preserving_modulus,
)
from gpc_census.generation.ressayre import (
    _sparse_determinant_coefficients,
    enumerate_admissible_hyperplanes,
)


def _pairs(enumeration):
    return {(candidate.h, candidate.z) for candidate in enumeration}


def test_hadamard_modulus_is_strictly_rank_preserving_through_rank_40() -> None:
    for d in range(4, 41):
        modulus, bound_squared = rank_preserving_modulus(3, d)
        assert modulus * modulus > bound_squared


def test_flat_enumerator_exactly_matches_rank6_basis_oracle() -> None:
    flat_result = enumerate_admissible_hyperplanes_by_flats(3, 6)
    basis_result = enumerate_admissible_hyperplanes(3, 6)

    assert flat_result.modulus == 257
    assert len(flat_result.hyperplanes) == 362
    assert _pairs(flat_result.hyperplanes) == _pairs(basis_result)


def test_rank7_flat_enumerator_reaches_the_complete_new_reference_universe() -> None:
    result = enumerate_admissible_hyperplanes_by_flats(3, 7)

    assert result.flat_counts_by_rank == (1, 35, 595, 4655, 15505, 19019, 5341)
    assert len(result.hyperplanes) == 5341
    assert len({(candidate.h, candidate.z) for candidate in result.hyperplanes}) == 5341


def test_sparse_labeled_determinant_combines_signs_and_exact_cancellations() -> None:
    nonzero = (
        ((0, 1), (1, 1)),
        ((1, 1), (0, 1)),
    )
    cancelling = (
        ((0, 1), (0, 1)),
        ((1, 1), (1, 1)),
    )

    assert _sparse_determinant_coefficients(nonzero) == {
        (0, 0): 1,
        (1, 1): -1,
    }
    assert _sparse_determinant_coefficients(cancelling) == {}
