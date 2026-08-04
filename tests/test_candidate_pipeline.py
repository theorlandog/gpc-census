from __future__ import annotations

from dataclasses import replace

import pytest

from gpc_census.generation import candidate_pipeline
from gpc_census.generation.candidate_pipeline import (
    UnresolvedCandidatesError,
    screen_bdr_export,
    screen_tau_candidates,
)
from gpc_census.generation.model import RessayreEvaluationVerdict
from gpc_census.generation.ressayre import assess_tau_by_evaluation


RANK7_CERTIFIED = (-2, 1, 1, 1, 1, -2, -2)
RANK7_COEFFICIENT_TWO = (-1, 0, 0, 1, -1, 1, 0)
RANK7_TRACE_REJECTED = (-2, -2, -2, -2, 1, 1, 1)
RANK7_INADMISSIBLE = (-1, -1, -1, -1, -1, -1, -1)
RANK7_STRUCTURAL = (0, 0, 0, 0, 0, 0, -1)
RANK5_SYMBOLIC_ZERO = (0, 0, 1, -1, 0)


def test_normalization_structural_split_and_default_cyclic_deferral():
    result = screen_tau_candidates(
        [RANK7_CERTIFIED, tuple(2 * value for value in RANK7_CERTIFIED), RANK7_STRUCTURAL]
    )

    assert result.counts.input_tau_rows == 3
    assert result.counts.normalized_tau_rows == 2
    assert result.counts.duplicate_or_rescaled_rows == 1
    assert result.counts.structural_rows == 1
    assert result.counts.certified_rows == 1
    assert result.counts.cyclic_schubert_cross_checked_rows == 0
    assert result.counts.cyclic_schubert_deferred_rows == 1
    assert result.candidates[0].cyclic_schubert is None
    assert result.candidates[0].as_dict()["cyclic_schubert_status"] == "deferred"
    assert result.structural_rows[0].tau == RANK7_STRUCTURAL
    assert result.full_dimension.is_full_dimensional
    assert result.as_dict()["ambient_moment_cone_full_dimension"]["certificate_replayed"]


def test_all_evaluation_outcomes_are_separate_and_unresolved_is_preserved(monkeypatch):
    real_assess = candidate_pipeline.assess_tau_by_evaluation

    def force_one_unresolved(n: int, tau: tuple[int, ...], *, attempts: int):
        if tau == RANK7_CERTIFIED:
            candidate = real_assess(n, tau, attempts=attempts)
            return RessayreEvaluationVerdict(
                status="evaluation_unresolved",
                attempts=attempts,
                below_weight_count=candidate.below_weight_count,
                negative_root_count=candidate.negative_root_count,
                certificate=None,
            )
        return real_assess(n, tau, attempts=attempts)

    monkeypatch.setattr(candidate_pipeline, "assess_tau_by_evaluation", force_one_unresolved)
    result = screen_tau_candidates(
        [
            RANK7_COEFFICIENT_TWO,
            RANK7_TRACE_REJECTED,
            RANK7_CERTIFIED,
            RANK7_INADMISSIBLE,
            RANK7_STRUCTURAL,
        ],
        cyclic_cross_check=True,
        symbolic_fallback=False,
    )
    by_tau = {row.tau: row for row in result.candidates}

    assert by_tau[RANK7_COEFFICIENT_TWO].status == "certified"
    assert by_tau[RANK7_TRACE_REJECTED].status == "trace_rejected"
    assert by_tau[RANK7_CERTIFIED].status == "evaluation_unresolved"
    assert by_tau[RANK7_INADMISSIBLE].status == "inadmissible"
    assert result.counts.certified_rows == 1
    assert result.counts.trace_rejected_rows == 1
    assert result.counts.evaluation_unresolved_rows == 1
    assert result.counts.inadmissible_rows == 1
    assert result.counts.exact_rejected_rows == 2


def test_positive_nonunit_cyclic_coefficient_never_rejects():
    result = screen_tau_candidates(
        [RANK7_COEFFICIENT_TWO],
        cyclic_cross_check=True,
        symbolic_fallback=False,
    )
    row = result.candidates[0]

    assert row.status == "certified"
    assert row.cyclic_schubert is not None
    assert row.cyclic_schubert.coefficient == 2
    assert result.counts.cyclic_schubert_nonzero_rows == 1
    assert result.counts.cyclic_schubert_coefficient_one_rows == 0


def test_exact_symbolic_zero_is_a_distinct_rejection():
    unresolved = screen_tau_candidates(
        [RANK5_SYMBOLIC_ZERO],
        ressayre_attempts=1,
        symbolic_fallback=False,
    )
    assert unresolved.candidates[0].status == "evaluation_unresolved"

    resolved = screen_tau_candidates(
        [RANK5_SYMBOLIC_ZERO],
        ressayre_attempts=1,
        symbolic_fallback=True,
    )
    row = resolved.candidates[0]
    assert row.status == "symbolic_zero_rejected"
    assert row.decision_class == "exact_rejected"
    assert row.symbolic_fallback is not None
    assert row.symbolic_fallback.status == "symbolic_zero_rejected"
    assert row.cyclic_schubert is None
    assert resolved.counts.symbolic_zero_rejected_rows == 1
    assert resolved.counts.evaluation_unresolved_rows == 0
    assert resolved.all_candidates_resolved


def test_sparse_exact_zero_fallback_does_not_require_sympy_determinant(monkeypatch):
    monkeypatch.setattr(
        candidate_pipeline,
        "certify_candidate",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("nonzero symbolic path must not run")
        ),
    )
    result = screen_tau_candidates(
        [RANK5_SYMBOLIC_ZERO],
        ressayre_attempts=1,
        symbolic_fallback=True,
    )
    assert result.candidates[0].status == "symbolic_zero_rejected"
    assert result.candidates[0].symbolic_fallback is not None
    assert "exact sparse" in result.candidates[0].symbolic_fallback.reason


def test_symbolic_fallback_can_certify_and_replay_after_evaluation_exhaustion(monkeypatch):
    monkeypatch.setattr(
        "gpc_census.generation.ressayre._deterministic_evaluation",
        lambda _h, _z, size, _attempt: (0,) * size,
    )
    result = screen_tau_candidates(
        [RANK7_CERTIFIED],
        ressayre_attempts=1,
        symbolic_fallback=True,
    )
    row = result.candidates[0]

    assert row.status == "certified"
    assert row.resolution_method == "exact_symbolic_determinant"
    assert row.ressayre_assessment is not None
    assert row.ressayre_assessment.status == "evaluation_unresolved"
    assert row.symbolic_fallback is not None
    assert row.symbolic_fallback.status == "certified"
    assert row.symbolic_fallback.certificate is not None
    assert row.symbolic_fallback.certificate_replayed_and_bound
    assert result.counts.symbolic_fallback_certified_rows == 1


def test_symbolic_resource_and_exception_outcomes_remain_unresolved(monkeypatch):
    monkeypatch.setattr(
        "gpc_census.generation.ressayre._deterministic_evaluation",
        lambda _h, _z, size, _attempt: (0,) * size,
    )
    limited = screen_tau_candidates(
        [RANK7_CERTIFIED],
        ressayre_attempts=1,
        symbolic_fallback=True,
        symbolic_max_determinant_order=1,
    )
    limited_row = limited.candidates[0]
    assert limited_row.status == "evaluation_unresolved"
    assert limited_row.symbolic_fallback is not None
    assert limited_row.symbolic_fallback.status == "resource_limit_unresolved"

    def exhaust_memory(*_args, **_kwargs):
        raise MemoryError("test resource exhaustion")

    monkeypatch.setattr(candidate_pipeline, "certify_candidate", exhaust_memory)
    failed = screen_tau_candidates(
        [RANK7_CERTIFIED],
        ressayre_attempts=1,
        symbolic_fallback=True,
    )
    failed_row = failed.candidates[0]
    assert failed_row.status == "evaluation_unresolved"
    assert failed_row.symbolic_fallback is not None
    assert failed_row.symbolic_fallback.status == "exception_unresolved"
    assert failed.counts.exact_rejected_rows == 0


def test_require_resolved_retains_complete_unresolved_result(monkeypatch):
    monkeypatch.setattr(
        "gpc_census.generation.ressayre._deterministic_evaluation",
        lambda _h, _z, size, _attempt: (0,) * size,
    )
    with pytest.raises(UnresolvedCandidatesError) as caught:
        screen_tau_candidates(
            [RANK7_CERTIFIED],
            ressayre_attempts=1,
            symbolic_fallback=False,
            require_resolved=True,
        )

    assert caught.value.result.counts.evaluation_unresolved_rows == 1
    assert caught.value.result.as_dict()["candidates"][0]["status"] == (
        "evaluation_unresolved"
    )


def test_tampered_evaluation_certificate_fails_closed(monkeypatch):
    assessment = assess_tau_by_evaluation(3, RANK7_CERTIFIED)
    assert assessment.certificate is not None
    tampered_certificate = replace(
        assessment.certificate,
        determinant_value=assessment.certificate.determinant_value + 1,
    )
    tampered = replace(assessment, certificate=tampered_certificate)
    monkeypatch.setattr(
        candidate_pipeline,
        "assess_tau_by_evaluation",
        lambda _n, _tau, *, attempts: replace(tampered, attempts=attempts),
    )

    with pytest.raises(AssertionError, match="replay or oriented tau binding"):
        screen_tau_candidates([RANK7_CERTIFIED], symbolic_fallback=False)


def test_internal_assessor_value_error_is_not_mislabeled_inadmissible(monkeypatch):
    def broken_assessor(*_args, **_kwargs):
        raise ValueError("internal evaluator failure")

    monkeypatch.setattr(candidate_pipeline, "assess_tau_by_evaluation", broken_assessor)
    with pytest.raises(ValueError, match="internal evaluator failure"):
        screen_tau_candidates([RANK7_CERTIFIED], symbolic_fallback=False)


def test_symbolic_zero_inference_replays_trace_counts(monkeypatch):
    real = assess_tau_by_evaluation(3, RANK5_SYMBOLIC_ZERO, attempts=1)
    inconsistent = replace(
        real,
        below_weight_count=real.below_weight_count + 1,
        negative_root_count=real.negative_root_count + 1,
    )
    monkeypatch.setattr(
        candidate_pipeline,
        "assess_tau_by_evaluation",
        lambda _n, _tau, *, attempts: replace(inconsistent, attempts=attempts),
    )

    with pytest.raises(AssertionError, match="trace counts failed exact replay"):
        screen_tau_candidates(
            [RANK5_SYMBOLIC_ZERO],
            ressayre_attempts=1,
            symbolic_fallback=True,
        )


@pytest.mark.parametrize("case", ["false_trace_rejection", "mismatched_certificate_counts"])
def test_tri_state_assessment_invariants_fail_closed(monkeypatch, case: str):
    real = assess_tau_by_evaluation(3, RANK7_CERTIFIED)
    assert real.certificate is not None
    if case == "false_trace_rejection":
        forged = RessayreEvaluationVerdict(
            status="trace_rejected",
            attempts=0,
            below_weight_count=real.below_weight_count,
            negative_root_count=real.below_weight_count,
            certificate=None,
        )
        match = "trace_rejected assessment is internally inconsistent"
    else:
        forged = replace(
            real,
            below_weight_count=real.below_weight_count + 1,
            negative_root_count=real.negative_root_count + 1,
        )
        match = "counts disagree with its certificate"
    monkeypatch.setattr(
        candidate_pipeline,
        "assess_tau_by_evaluation",
        lambda _n, _tau, *, attempts: replace(forged, attempts=forged.attempts),
    )

    with pytest.raises(AssertionError, match=match):
        screen_tau_candidates([RANK7_CERTIFIED], symbolic_fallback=False)


def test_rank6_full_dimension_failure_is_serialized_separately():
    result = screen_tau_candidates([(0, 0, 0, 0, 0, -1)])

    assert result.counts.structural_rows == 1
    assert result.counts.nonstructural_candidate_rows == 0
    assert not result.full_dimension.is_full_dimensional
    assert not result.as_dict()["ambient_moment_cone_full_dimension"]["established"]


def test_safe_bdr_parse_rank_binding_and_option_validation():
    text = "raise RuntimeError('must not execute')\nbrut_inequations = [(0, 0, 1, -1, 0)]"
    result = screen_bdr_export(text, expected_rank=5, symbolic_fallback=True)
    assert result.candidates[0].status == "symbolic_zero_rejected"

    with pytest.raises(ValueError, match="rank-6"):
        screen_bdr_export(text, expected_rank=6)
    with pytest.raises(ValueError, match="include_modular"):
        screen_tau_candidates([RANK7_CERTIFIED], include_modular=True)
    with pytest.raises(ValueError, match="prime smaller"):
        screen_tau_candidates(
            [RANK7_STRUCTURAL],
            cyclic_cross_check=True,
            include_modular=True,
            modulus=4,
        )
    with pytest.raises(ValueError, match="at least the rank"):
        screen_tau_candidates(
            [RANK7_STRUCTURAL],
            cyclic_cross_check=True,
            include_modular=True,
            modulus=5,
        )
    with pytest.raises(ValueError, match="positive integer"):
        screen_tau_candidates([RANK7_CERTIFIED], ressayre_attempts=0)
    with pytest.raises(ValueError, match="positive integer or None"):
        screen_tau_candidates(
            [RANK7_CERTIFIED], symbolic_max_determinant_order=0
        )
    with pytest.raises(TypeError, match="excluding bool"):
        screen_tau_candidates([(True, 0, 0)])
    with pytest.raises(ValueError, match="workers must be a positive integer"):
        screen_tau_candidates([RANK7_CERTIFIED], workers=0)


def test_parallel_screening_preserves_canonical_results():
    rows = (RANK7_CERTIFIED, RANK7_TRACE_REJECTED, RANK7_INADMISSIBLE)
    sequential = screen_tau_candidates(rows, ressayre_attempts=1, workers=1)
    parallel = screen_tau_candidates(rows, ressayre_attempts=1, workers=2)
    assert parallel == sequential
    with pytest.raises(TypeError, match="require_resolved"):
        screen_tau_candidates([RANK7_CERTIFIED], require_resolved=1)
