from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "screen_fixed_n3_candidates.py"
RANK7_CERTIFIED = (-2, 1, 1, 1, 1, -2, -2)
RANK7_STRUCTURAL = (0, 0, 0, 0, 0, 0, -1)
RANK5_SYMBOLIC_ZERO = (0, 0, 1, -1, 0)


def _module():
    spec = importlib.util.spec_from_file_location("screen_fixed_n3_candidates", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _write_export(path: Path, rows: list[tuple[int, ...]]) -> None:
    path.write_text(f"brut_inequations = {rows!r}\n", encoding="utf-8")


def test_payload_serializes_counts_certificates_and_conservative_scope(tmp_path: Path):
    export = tmp_path / "rank7.py"
    _write_export(
        export,
        [
            RANK7_CERTIFIED,
            tuple(2 * value for value in RANK7_CERTIFIED),
            RANK7_STRUCTURAL,
        ],
    )
    payload = _module().build_payload(7, export)
    screening = payload["screening"]

    assert payload["schema_version"] == 1
    assert screening["counts"]["input_tau_rows"] == 3
    assert screening["counts"]["duplicate_or_rescaled_rows"] == 1
    assert screening["counts"]["certified_rows"] == 1
    assert screening["counts"]["structural_rows"] == 1
    assert screening["counts"]["cyclic_schubert_cross_checked_rows"] == 0
    assert screening["counts"]["cyclic_schubert_deferred_rows"] == 1
    assert screening["ambient_moment_cone_full_dimension"]["established"]
    assert screening["ambient_moment_cone_full_dimension"]["certificate_replayed"]
    assert payload["options"]["symbolic_fallback_request"] == "auto_rank_at_most_13"
    assert payload["options"]["symbolic_fallback_effective"]
    assert payload["proof_status"]["candidate_exhaustiveness"] == "not_established"
    assert payload["proof_status"]["facetness"] == "not_established"
    assert payload["proof_status"]["system_completeness"] == "not_established"


def test_cli_require_resolved_writes_unresolved_artifact_then_exits_two(
    tmp_path: Path,
):
    export = tmp_path / "rank5.py"
    output = tmp_path / "screen.json"
    _write_export(export, [RANK5_SYMBOLIC_ZERO])

    status = _module().main(
        [
            "--input",
            str(export),
            "--rank",
            "5",
            "--output",
            str(output),
            "--ressayre-attempts",
            "1",
            "--no-symbolic-fallback",
            "--require-resolved",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert status == 2
    assert not payload["operational_gate"]["passed"]
    assert payload["operational_gate"]["failure_exit_status"] == 2
    assert payload["screening"]["counts"]["evaluation_unresolved_rows"] == 1
    assert payload["screening"]["candidates"][0]["status"] == "evaluation_unresolved"


def test_cli_default_symbolic_fallback_resolves_exact_zero(tmp_path: Path):
    export = tmp_path / "rank5.py"
    output = tmp_path / "screen.json"
    _write_export(export, [RANK5_SYMBOLIC_ZERO])

    status = _module().main(
        [
            "--input",
            str(export),
            "--rank",
            "5",
            "--output",
            str(output),
            "--ressayre-attempts",
            "1",
            "--require-resolved",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert status == 0
    assert payload["operational_gate"]["passed"]
    assert payload["screening"]["counts"]["symbolic_zero_rejected_rows"] == 1
    assert payload["screening"]["counts"]["evaluation_unresolved_rows"] == 0


def test_payload_builder_rejects_nonboolean_resolution_gate(tmp_path: Path):
    export = tmp_path / "rank5.py"
    _write_export(export, [RANK5_SYMBOLIC_ZERO])

    with pytest.raises(TypeError, match="require_resolved"):
        _module().build_payload(5, export, require_resolved=1)
