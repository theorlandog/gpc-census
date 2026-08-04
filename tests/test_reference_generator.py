from __future__ import annotations

import json

from gpc_census.generation.reference_generator import (
    generate_fixed_n3_reference_rank,
    verify_reference_manifest,
)


EXPECTED_RANK7_TAUS = {
    (-2, 1, 1, 1, 1, -2, -2),
    (1, -2, 1, 1, -2, 1, -2),
    (1, 1, -2, -2, 1, 1, -2),
    (1, 1, -2, 1, -2, -2, 1),
}


def test_rank7_reference_generator_is_exhaustive_resolved_and_exact() -> None:
    result = generate_fixed_n3_reference_rank(7)

    assert result.proved_complete
    assert result.enumeration.flat_counts_by_rank == (
        1,
        35,
        595,
        4655,
        15505,
        19019,
        5341,
    )
    assert result.screening.counts.trace_rejected_rows == 10133
    assert result.screening.counts.symbolic_zero_rejected_rows == 499
    assert result.screening.counts.certified_rows == 49
    assert result.screening.counts.evaluation_unresolved_rows == 0
    assert result.candidate_system.counts.redundant_rows == 45
    assert {row.tau for row in result.candidate_system.retained_rows} == (
        EXPECTED_RANK7_TAUS
    )


def test_rank7_compact_manifest_round_trips_and_rejects_tampering() -> None:
    result = generate_fixed_n3_reference_rank(7)
    payload = json.loads(json.dumps(result.as_dict(), sort_keys=True))

    assert verify_reference_manifest(payload)
    payload["screening_counts"]["certified_rows"] += 1
    assert not verify_reference_manifest(payload)
