"""Guards for the external half of the constructor benchmark.

These pin the shape and the load-bearing caveats of the BDR measurements, not
the wall clock itself, which is machine dependent and excluded from the
artifact's canonical hash.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

BENCHMARK = DATA / "bdr_external_benchmark.json"
SYSTEMS = ("(3,7)", "(3,8)", "(3,9)")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = _load("bdr_external_benchmark")


def _artifact() -> dict:
    return json.loads(BENCHMARK.read_text())


def test_the_backend_was_executed_at_the_pinned_revision():
    backend = _artifact()["external_backend"]
    assert backend["executed_here"] is True
    assert backend["pinned_revision"] == MODULE.PINNED_REVISION
    assert backend["license"] == "MIT"


def test_every_system_was_measured_with_repeats():
    artifact = _artifact()
    assert set(artifact["systems"]) == set(SYSTEMS)
    assert artifact["protocol"]["repeats"] == MODULE.REPEATS >= 3
    for name in SYSTEMS:
        arm = artifact["systems"][name]["probabilistic"]
        assert len(arm["total_wall_seconds_per_run"]) == MODULE.REPEATS
        assert arm["total_wall_seconds_median"] > 0
        assert arm["stages"]
        assert arm["stage_counts"]


def test_the_pipeline_reproduces_the_published_row_counts():
    """Plus the structural row upstream carries and this repository does not."""
    expected = {"(3,7)": 4, "(3,8)": 32, "(3,9)": 53}
    systems = _artifact()["systems"]
    for name, rows in expected.items():
        assert systems[name]["probabilistic"]["retained_rows"] == rows


def test_only_the_probabilistic_arm_is_claimed():
    artifact = _artifact()
    assert artifact["protocol"]["arms"] == ["probabilistic"]
    for record in artifact["systems"].values():
        assert set(record) == {"probabilistic"}


def test_the_exact_arm_is_recorded_as_unavailable_with_evidence():
    """A limitation must be measured, not asserted."""
    unavailable = _artifact()["protocol"]["exact_arm_unavailable"]
    assert unavailable["available"] is False
    assert unavailable["largest_variable_count_that_factors"] == 48
    assert unavailable["smallest_variable_count_that_fails"] == 56
    probe = {row["variables"]: row["factors"] for row in unavailable["probe"]}
    assert probe[48] is True
    assert probe[56] is False
    for count, needed in unavailable["variables_needed"].items():
        assert needed >= 70, count
        assert needed > unavailable["largest_variable_count_that_factors"]


def test_the_caveat_names_the_asymmetry():
    caveat = _artifact()["protocol"]["caveat"]
    assert "not a like-for-like" in caveat
    assert "probabilistic" in caveat
    assert "certificate" in caveat


def test_the_measurements_are_marked_post_hoc():
    assert "post hoc" in _artifact()["not_preregistered"]


def test_the_birationality_stage_is_present_and_costly():
    """The stage with no local counterpart is measured, not inferred."""
    stages = {
        stage["stage"]: stage["wall_seconds"]
        for stage in _artifact()["systems"]["(3,8)"]["probabilistic"]["stages"]
    }
    assert "BirationalityStep" in stages
    assert stages["BirationalityStep"] > 0
    assert "PiDominancyStep" in stages
    assert "TauCandidatesStep" in stages


def test_candidate_volume_reaching_the_filters():
    """The number that closed the open candidate-generation question."""
    counts = {
        entry["stage"]: entry["pending"] + entry["validated"]
        for entry in _artifact()["systems"]["(3,8)"]["probabilistic"]["stage_counts"]
    }
    assert counts["TPiPreComputationStep"] == 173
    assert counts["TauCandidatesStep"] > counts["StabilizerConditionStep"]


def test_timings_are_excluded_from_the_canonical_hash():
    """A rerun that reproduces every count must hash identically."""
    import hashlib

    artifact = _artifact()
    stored = artifact.pop("canonical_sha256_without_this_field_or_timings")
    artifact.pop("machine", None)
    for system in artifact["systems"].values():
        for arm in system.values():
            for field in (
                "total_wall_seconds_per_run",
                "total_wall_seconds_median",
                "total_wall_seconds_spread",
                "stages",
            ):
                arm.pop(field, None)
    recomputed = hashlib.sha256(
        json.dumps(artifact, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert recomputed == stored
