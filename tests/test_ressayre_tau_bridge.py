from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest

from gpc_census.generation.bdr import tau_from_constraint
from gpc_census.generation.model import IntegralConstraint, RessayreCertificate
from gpc_census.generation.ressayre import (
    _modular_row_witness,
    admissible_candidate_from_tau,
    assess_candidate_by_evaluation,
    certify_tau_by_evaluation,
    verify_tau_ressayre_certificate,
    verify_ressayre_certificate,
)


ROOT = Path(__file__).parents[1]


@pytest.mark.parametrize("d", range(7, 11))
def test_every_stored_n3_gpc_has_an_exact_evaluated_ressayre_witness(d):
    table = json.loads(
        (ROOT / "src" / "gpc_census" / "data" / "constraints.json").read_text()
    )
    certificates = []
    for row in table[f"3,{d}"]["inequalities"]:
        tau = tau_from_constraint(
            IntegralConstraint(tuple(row["coeffs"]), row["rhs"]), 3
        )
        candidate = admissible_candidate_from_tau(3, tau)
        assert sum(candidate.h) == 0
        certificate = certify_tau_by_evaluation(3, tau)
        assert certificate is not None
        assert certificate.determinant_value != 0
        assert verify_ressayre_certificate(3, d, certificate)
        assert verify_tau_ressayre_certificate(3, tau, certificate)
        certificates.append(certificate)

    assert len(certificates) == len(table[f"3,{d}"]["inequalities"])


def test_tau_scaling_and_structural_pauli_certificate():
    tau = (-2, 1, 1, 1, 1, -2, -2)
    candidate = admissible_candidate_from_tau(3, tau)
    assert admissible_candidate_from_tau(3, tuple(5 * value for value in tau)) == candidate

    structural_tau = (0, 0, 0, 0, 0, 0, -1)
    structural = certify_tau_by_evaluation(3, structural_tau)
    assert structural is not None
    assert structural.negative_roots == ()
    assert structural.below_indices == ()
    assert structural.determinant_value == 1
    assert verify_tau_ressayre_certificate(3, structural_tau, structural)


def test_evaluation_search_never_turns_failure_into_a_negative_certificate():
    tau = (-2, 1, 1, 1, 1, -2, -2)
    assert certify_tau_by_evaluation(3, tau, attempts=1) is not None


def test_tau_binding_rejects_orientation_and_payload_tampering():
    tau = (-2, 1, 1, 1, 1, -2, -2)
    certificate = certify_tau_by_evaluation(3, tau)
    assert certificate is not None
    assert not verify_tau_ressayre_certificate(3, tuple(-x for x in tau), certificate)
    tampered = RessayreCertificate(
        h=certificate.h,
        z=certificate.z,
        admissibility_witness=certificate.admissibility_witness,
        on_indices=certificate.on_indices,
        below_indices=certificate.below_indices,
        negative_roots=certificate.negative_roots,
        determinant_evaluation=certificate.determinant_evaluation,
        determinant_value=certificate.determinant_value + 1,
    )
    assert not verify_tau_ressayre_certificate(3, tau, tampered)
    float_tampered = replace(
        certificate,
        determinant_evaluation=tuple(
            float(value) for value in certificate.determinant_evaluation
        ),
        determinant_value=float(certificate.determinant_value),
    )
    assert not verify_tau_ressayre_certificate(3, tau, float_tampered)


def test_candidate_validation_and_unresolved_status(monkeypatch):
    tau = (-2, 1, 1, 1, 1, -2, -2)
    candidate = admissible_candidate_from_tau(3, tau)
    scaled = replace(
        candidate,
        h=tuple(2 * value for value in candidate.h),
        z=2 * candidate.z,
    )
    with pytest.raises(ValueError, match="primitive"):
        assess_candidate_by_evaluation(3, 7, scaled)

    monkeypatch.setattr(
        "gpc_census.generation.ressayre._deterministic_evaluation",
        lambda _h, _z, size, _attempt: (0,) * size,
    )
    verdict = assess_candidate_by_evaluation(3, 7, candidate, attempts=1)
    assert verdict.status == "evaluation_unresolved"
    assert verdict.certificate is None


def test_modular_rank_witness_tries_multiple_characteristics():
    rows = ((1, 0), (0, 2))
    assert _modular_row_witness(rows, 2) == (0, 1)
    assert _modular_row_witness(((1, 2), (2, 4)), 2) is None
