from __future__ import annotations

import importlib.util
import json
from dataclasses import replace
from pathlib import Path

import pytest

from gpc_census.generation import candidate_pipeline
from gpc_census.generation.candidate_system import (
    BLOCKED_KNOWN_MISMATCH_STATUS,
    BLOCKED_UNRESOLVED_STATUS,
    FINALIZED_KNOWN_STATUS,
    CandidateSystemInvariantError,
    compose_candidate_system,
)
from gpc_census.generation.redundancy import (
    verify_ordered_slice_reduction_certificate,
)


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "generate_fixed_n3_candidate_system.py"
RANK7_ROWS = (
    (-2, 1, 1, 1, 1, -2, -2),
    (1, -2, 1, 1, -2, 1, -2),
    (1, 1, -2, -2, 1, 1, -2),
    (1, 1, -2, 1, -2, -2, 1),
    (0, 0, 0, 0, 0, 0, -1),
)


def _screen(rows: tuple[tuple[int, ...], ...]):
    return candidate_pipeline.screen_tau_candidates(
        rows,
        ressayre_attempts=1,
        cyclic_cross_check=False,
        symbolic_fallback=None,
    )


def _write_export(path: Path, rows: tuple[tuple[int, ...], ...]) -> None:
    path.write_text(f"brut_inequations = {list(rows)!r}\n", encoding="utf-8")


def _script_module():
    spec = importlib.util.spec_from_file_location("generate_fixed_n3_candidate_system", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_rank7_composition_replays_and_matches_known_system() -> None:
    result = compose_candidate_system(_screen(RANK7_ROWS))

    assert result.status == FINALIZED_KNOWN_STATUS
    assert result.finalized
    assert result.exit_status == 0
    assert result.counts.certified_rows == 4
    assert result.counts.redundant_rows == 0
    assert result.counts.retained_rows == 4
    assert result.counts.retained_cyclic_cross_checks == 4
    assert result.known_system_regression is not None
    assert result.known_system_regression.passed
    assert result.reduction is not None
    assert verify_ordered_slice_reduction_certificate(result.reduction)
    assert all(row.cyclic_schubert.coefficient_one for row in result.retained_rows)
    payload = result.as_dict()
    assert payload["final_system"] is not None
    assert json.loads(json.dumps(payload, sort_keys=True))["finalized"]


def test_certified_structural_extra_is_removed_before_known_gate() -> None:
    pauli_upper = (2, -1, -1, -1, -1, -1, -1)
    result = compose_candidate_system(_screen(RANK7_ROWS + (pauli_upper,)))

    assert result.status == FINALIZED_KNOWN_STATUS
    assert result.counts.certified_rows == 5
    assert result.counts.redundant_rows == 1
    assert result.counts.retained_rows == 4
    assert result.known_system_regression is not None
    assert result.known_system_regression.generated_only == ()
    assert result.known_system_regression.published_only == ()


def test_known_rank_mismatch_blocks_final_system() -> None:
    screening = _screen((RANK7_ROWS[0], RANK7_ROWS[-1]))
    result = compose_candidate_system(screening)

    assert result.status == BLOCKED_KNOWN_MISMATCH_STATUS
    assert not result.finalized
    assert result.exit_status == 3
    assert result.reduction is not None
    assert result.known_system_regression is not None
    assert not result.known_system_regression.passed
    assert len(result.known_system_regression.published_only) == 3
    assert result.as_dict()["final_system"] is None


def test_unresolved_candidate_blocks_reduction_and_cli_writes_diagnostic(
    tmp_path: Path,
) -> None:
    unresolved_tau = (-1, -1, 2, 2, -1, -1, -1, -4, -4, -1)
    structural_tau = (0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
    screening = candidate_pipeline.screen_tau_candidates(
        (unresolved_tau, structural_tau),
        ressayre_attempts=1,
        cyclic_cross_check=False,
        symbolic_fallback=False,
    )
    result = compose_candidate_system(screening)

    assert result.status == BLOCKED_UNRESOLVED_STATUS
    assert result.exit_status == 2
    assert result.reduction is None
    assert result.retained_rows == ()
    assert result.as_dict()["final_system"] is None

    export = tmp_path / "rank10.py"
    output = tmp_path / "diagnostic.json"
    _write_export(export, (unresolved_tau, structural_tau))
    exit_status = _script_module().main(
        [
            "--input",
            str(export),
            "--rank",
            "10",
            "--output",
            str(output),
            "--no-symbolic-fallback",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert exit_status == 2
    assert payload["options"]["ressayre_evaluation_attempt_limit"] == 1
    assert payload["operational_gate"] == {
        "exit_status": 2,
        "passed": False,
        "status": BLOCKED_UNRESOLVED_STATUS,
    }
    assert payload["candidate_system"]["screening"]["all_candidates_resolved"] is False
    assert payload["candidate_system"]["ordered_slice_reduction_certificate"] is None
    assert payload["candidate_system"]["final_system"] is None


def test_binding_and_canonical_order_tampering_fail_closed() -> None:
    screening = _screen(RANK7_ROWS)
    first = screening.candidates[0]
    unbound = replace(first, certificate_replayed_and_bound=False)
    with pytest.raises(CandidateSystemInvariantError, match="replay-and-binding"):
        compose_candidate_system(
            replace(screening, candidates=(unbound, *screening.candidates[1:]))
        )

    with pytest.raises(CandidateSystemInvariantError, match="canonical normalized tau order"):
        compose_candidate_system(
            replace(
                screening,
                candidates=(
                    screening.candidates[1],
                    screening.candidates[0],
                    *screening.candidates[2:],
                ),
            )
        )
