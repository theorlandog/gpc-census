from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "fixed_n3_series.py"


def _module():
    spec = importlib.util.spec_from_file_location("fixed_n3_series", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_rank7_payload_is_certified_and_matches_lookup(tmp_path: Path):
    export = tmp_path / "rank7.py"
    export.write_text(
        "brut_inequations = [\n"
        "(-2, 1, 1, 1, 1, -2, -2),\n"
        "(1, -2, 1, 1, -2, 1, -2),\n"
        "(1, 1, -2, -2, 1, 1, -2),\n"
        "(1, 1, -2, 1, -2, -2, 1),\n"
        "(0, 0, 0, 0, 0, 0, -1),\n"
        "]\n"
    )
    payload = _module().build_payload([(7, export)])
    row = payload["rank_systems"][0]

    assert payload["schema_version"] == 3
    assert row["system"]["counts"]["inequality_rows"] == 4
    assert row["system"]["counts"]["structural_rows"] == 1
    assert row["system"]["all_coefficients_one"]
    assert row["ambient_moment_cone_full_dimension_certificate"]["is_full_dimensional"]
    assert row["system"]["all_nonstructural_ressayre_certified"]
    assert all(entry["ressayre"] for entry in row["system"]["inequalities"])
    assert row["known_system_regression"]["homogeneous_oriented_sets_equal"]
    reduction = row["ordered_slice_redundancy_certificate"]
    assert reduction["candidate_rows"] == [
        entry["constraint"] for entry in row["system"]["inequalities"]
    ]
    assert reduction["retained_indices"] == [0, 1, 2, 3]
    assert reduction["redundant_indices"] == []
    assert reduction["removals"] == []
    assert len(reduction["facets"]) == 4
    assert (
        row["proof_status"]["exact_redundancy_elimination"]
        == "proved_and_replayable_for_supplied_rows"
    )
    assert (
        row["proof_status"]["ordered_slice_relative_facetness"]
        == "proved_and_replayable_for_supplied_rows"
    )
    assert row["proof_status"]["candidate_exhaustiveness"].endswith("obligation")
    assert row["proof_status"]["system_completeness"] == "not_proved_here"
    assert row["proof_status"]["polyhedral_certificate_scope"].startswith("does_not_by_itself")


def test_rank_payload_rejects_a_removed_final_row(tmp_path: Path) -> None:
    export = tmp_path / "rank7.py"
    export.write_text(
        "brut_inequations = [\n"
        "(-2, 1, 1, 1, 1, -2, -2),\n"
        "(1, -2, 1, 1, -2, 1, -2),\n"
        "(1, 1, -2, -2, 1, 1, -2),\n"
        "(1, 1, -2, 1, -2, -2, 1),\n"
        "(2, -1, -1, -1, -1, -1, -1),\n"
        "(0, 0, 0, 0, 0, 0, -1),\n"
        "]\n"
    )
    module = _module()

    with pytest.raises(ValueError, match=r"removed.*indices \[1\]"):
        module.build_rank_payload(
            7,
            export,
            include_modular=False,
            ressayre_attempts=32,
        )
