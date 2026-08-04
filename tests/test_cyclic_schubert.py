from __future__ import annotations

import json

import pytest

from gpc_census.generation.cyclic_schubert import (
    DEFAULT_MODULUS,
    ExactCyclicSchubertCertificate,
    certify_cyclic_schubert,
    prepare_cyclic_schubert_problem,
    verify_exact_cyclic_schubert_certificate,
    verify_modular_cyclic_schubert_certificate,
)


def test_borland_dennis_gpc_has_exact_coefficient_one() -> None:
    verdict = certify_cyclic_schubert((1, 1, 0, 1, 0, 0), 2)

    assert verdict.status == "coefficient_one"
    assert verdict.coefficient == 1
    assert verdict.coefficient_one
    assert verdict.exact_certificate is not None
    assert verdict.modular_certificate is not None
    assert verdict.modular_certificate.modulus == DEFAULT_MODULUS
    assert verdict.modular_certificate.residue == 1
    assert verdict.modular_certificate.proves_nonzero
    assert verify_exact_cyclic_schubert_certificate(verdict.exact_certificate)
    assert verify_modular_cyclic_schubert_certificate(verdict.modular_certificate)


def test_nonunit_integer_coefficient_is_not_accepted() -> None:
    verdict = certify_cyclic_schubert((-2, -2, -1, 0, -1, -2), -4)

    assert verdict.status == "coefficient_nonunit"
    assert verdict.coefficient == 2
    assert not verdict.coefficient_one
    assert verdict.modular_certificate is not None
    assert verdict.modular_certificate.residue == 2
    assert verdict.modular_certificate.proves_nonzero


def test_trace_shift_and_positive_rescaling_preserve_coefficient() -> None:
    original = certify_cyclic_schubert((0, 1, 1, 1, 1, 0, 0), 2)
    shifted_scaled = certify_cyclic_schubert((6, 9, 9, 9, 9, 6, 6), 24)

    assert original.coefficient == shifted_scaled.coefficient == 1
    assert original.problem.source_permutation == shifted_scaled.problem.source_permutation
    assert original.problem.divided_difference_word == (
        shifted_scaled.problem.divided_difference_word
    )
    assert original.problem.upper_weight_subsets == shifted_scaled.problem.upper_weight_subsets


def test_degree_mismatch_is_a_deterministic_negative_verdict() -> None:
    verdict = certify_cyclic_schubert((1, 0, 0, 0), 0, particle_number=2)

    assert verdict.status == "degree_mismatch"
    assert verdict.coefficient is None
    assert not verdict.coefficient_one
    assert verdict.exact_certificate is None
    assert verdict.modular_certificate is None


def test_payload_is_deterministic_and_json_serializable() -> None:
    first = certify_cyclic_schubert((0, 1, 1, 1, 1, 0, 0), 2)
    second = certify_cyclic_schubert((0, 1, 1, 1, 1, 0, 0), 2)

    assert first == second
    assert json.dumps(first.as_dict(), sort_keys=True, separators=(",", ":")) == json.dumps(
        second.as_dict(), sort_keys=True, separators=(",", ":")
    )


def test_tampered_exact_certificate_fails_replay() -> None:
    verdict = certify_cyclic_schubert((1, 1, 0, 1, 0, 0), 2)
    assert verdict.exact_certificate is not None
    certificate = verdict.exact_certificate
    tampered = ExactCyclicSchubertCertificate(
        problem=certificate.problem,
        evaluation_base=certificate.evaluation_base,
        coefficient=certificate.coefficient + 1,
        evaluation_states=certificate.evaluation_states,
    )

    assert not verify_exact_cyclic_schubert_certificate(tampered)


def test_rejects_composite_modulus_and_inconsistent_problem() -> None:
    with pytest.raises(ValueError, match="prime"):
        certify_cyclic_schubert((1, 1, 0, 1, 0, 0), 2, modulus=15)

    problem = prepare_cyclic_schubert_problem((1, 1, 0, 1, 0, 0), 2)
    inconsistent = type(problem)(
        particle_number=problem.particle_number,
        coefficients=problem.coefficients,
        rhs=problem.rhs,
        dominant_coefficients=problem.dominant_coefficients,
        source_permutation=problem.source_permutation,
        divided_difference_word=(),
        upper_weight_subsets=problem.upper_weight_subsets,
        threshold_is_weight=problem.threshold_is_weight,
    )
    with pytest.raises(ValueError, match="inconsistent"):
        from gpc_census.generation.cyclic_schubert import (
            exact_cyclic_schubert_certificate,
        )

        exact_cyclic_schubert_certificate(inconsistent)
