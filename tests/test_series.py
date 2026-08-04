from __future__ import annotations

import json

import pytest

from gpc_census.generation.model import IntegralConstraint
from gpc_census.generation.series import certify_tau_rank_system


RANK7_TAUS = (
    (-2, 1, 1, 1, 1, -2, -2),
    (1, -2, 1, 1, -2, 1, -2),
    (1, 1, -2, -2, 1, 1, -2),
    (1, 1, -2, 1, -2, -2, 1),
    (0, 0, 0, 0, 0, 0, -1),
)


def test_rank7_composition_certifies_and_separates_every_row() -> None:
    raw = RANK7_TAUS + (tuple(2 * value for value in RANK7_TAUS[0]),)
    system = certify_tau_rank_system(raw, 3, include_modular=True)

    assert system.d == 7
    assert system.counts.as_dict() == {
        "input_tau_rows": 6,
        "normalized_tau_rows": 5,
        "duplicate_or_rescaled_rows": 1,
        "inequality_rows": 4,
        "structural_rows": 1,
        "cyclic_schubert_coefficient_one_rows": 5,
        "ressayre_certified_inequality_rows": 4,
    }
    assert system.all_coefficients_one
    assert system.all_nonstructural_ressayre_certified
    assert system.structural_inequalities[0].constraint == IntegralConstraint(
        (1, 1, 1, 1, 1, 1, 0), 3
    )
    for row in system.inequalities + system.structural_inequalities:
        assert row.schubert.coefficient == 1
        assert row.schubert.exact_certificate is not None
        assert row.schubert.modular_certificate is not None
        assert row.schubert.modular_certificate.residue == 1
    assert all(row.ressayre is not None for row in system.inequalities)
    assert system.structural_inequalities[0].ressayre is None
    assert system.structural_inequalities[0].row_kind == "pauli_lower_bound"

    reversed_system = certify_tau_rank_system(tuple(reversed(raw)), 3, include_modular=True)
    assert json.dumps(system.as_dict(), sort_keys=True, separators=(",", ":")) == json.dumps(
        reversed_system.as_dict(), sort_keys=True, separators=(",", ":")
    )


def test_rank7_composition_rejects_a_non_coefficient_one_row() -> None:
    with pytest.raises(ValueError, match="coefficient-one requirement failed"):
        certify_tau_rank_system(((1, 0, 0, 0, 0, 0, 0),), 3)
