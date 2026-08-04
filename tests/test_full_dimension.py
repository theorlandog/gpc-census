from __future__ import annotations

import pytest

from gpc_census.generation.full_dimension import (
    CERTIFICATE_PRIME,
    SEED_RULE,
    full_dimension_certificate,
    require_full_dimension,
    verify_full_dimension_certificate,
)


@pytest.mark.parametrize(
    ("d", "skew_rank", "symmetric_rank"),
    [
        (7, 21, 28),
        (8, 28, 36),
        (9, 36, 45),
        (10, 45, 55),
        (11, 55, 66),
        (12, 66, 78),
        (13, 78, 91),
    ],
)
def test_exact_full_dimension_certificates(d: int, skew_rank: int, symmetric_rank: int):
    certificate = require_full_dimension(d)

    assert certificate.prime == CERTIFICATE_PRIME
    assert certificate.seed_rule == SEED_RULE
    assert certificate.skew_rank == skew_rank
    assert certificate.symmetric_rank == symmetric_rank
    assert certificate.expected_skew_rank == skew_rank
    assert certificate.expected_symmetric_rank == symmetric_rank
    assert certificate.skew_minor_determinant_mod_prime not in (None, 0)
    assert certificate.symmetric_minor_determinant_mod_prime not in (None, 0)
    assert verify_full_dimension_certificate(certificate)


def test_rank_six_is_rejected_exactly():
    certificate = full_dimension_certificate(6)

    assert certificate.skew_rank == 15
    assert certificate.symmetric_rank == 19
    assert certificate.expected_skew_rank == 15
    assert certificate.expected_symmetric_rank == 21
    assert not certificate.is_full_dimensional
    assert certificate.symmetric_minor_determinant_mod_prime is None
    assert verify_full_dimension_certificate(certificate)
    with pytest.raises(ValueError, match="not full-dimensional"):
        require_full_dimension(6)
