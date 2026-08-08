"""Guards for the rank-7 to rank-9 full-population compression audit."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

RANK9_COMPRESSION_AUDIT = DATA / "rank9_compression_audit.json"
FULL_RANK78_COMPRESSION_AUDIT = DATA / "full_rank78_compression_audit.json"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _audit() -> dict:
    return json.loads(RANK9_COMPRESSION_AUDIT.read_text())


def _published() -> dict:
    return json.loads(FULL_RANK78_COMPRESSION_AUDIT.read_text())


def test_ranks_seven_and_eight_reproduce_the_published_population():
    """The gap the frontier document declared: these gates are now regenerated."""
    audit = _audit()
    published = _published()
    for system in ("(3,7)", "(3,8)"):
        mine = audit["systems"][system]
        theirs = published["systems"][system]
        assert mine["trace_survivors"] == theirs["trace_survivors"]
        assert mine["nonstructural_survivors"] == theirs["nonstructural_survivors"]
        assert mine["hall_feasible"] == theirs["hall_feasible"]
        assert mine["hall_zero"] == theirs["hall_zero"]
        assert mine["hole_distribution_all"] == theirs["hole_distribution_all"]
        assert mine["singleton_signature_count"] == theirs["singleton_signature_count"]
        assert (
            mine["optimal_fusion_width_distribution"]
            == theirs["optimal_fusion_width_distribution_all"]
        )


def test_published_survivor_totals():
    audit = _audit()
    assert audit["systems"]["(3,7)"]["trace_survivors"] == 549
    assert audit["systems"]["(3,8)"]["trace_survivors"] == 6607


def test_rank_nine_population():
    rank9 = _audit()["systems"]["(3,9)"]
    assert rank9["trace_survivors"] == 135343
    assert rank9["nonstructural_survivors"] == 135341
    assert rank9["hall_feasible"] + rank9["hall_zero"] == rank9["trace_survivors"]


def test_multi_hole_survivors_appear_at_rank_nine_and_hall_kills_them_all():
    """The conjecture's first real test: 5,310 chances to fail, taken zero times."""
    audit = _audit()
    assert audit["systems"]["(3,7)"]["multi_hole_survivors"] == 0
    assert audit["systems"]["(3,8)"]["multi_hole_survivors"] == 0
    rank9 = audit["systems"]["(3,9)"]
    assert rank9["multi_hole_survivors"] == 5310
    assert rank9["hole_distribution_all"]["2"] == 1227
    assert rank9["hole_distribution_all"]["3"] == 4083
    assert rank9["multi_hole_hall_feasible"] == 0
    assert audit["summary"]["multi_hole_hall_feasible_anywhere"] == 0


def test_one_hole_onset_is_rank_nine_on_the_full_population():
    audit = _audit()
    assert audit["systems"]["(3,7)"]["one_hole_hall_feasible"] == 0
    rank8 = audit["systems"]["(3,8)"]
    assert rank8["one_hole_hall_feasible"] == 7
    assert rank8["one_hole_determinant_nonzero"] == 0
    rank9 = audit["systems"]["(3,9)"]
    assert rank9["one_hole_hall_feasible"] == 32
    assert rank9["one_hole_determinant_nonzero"] == 3


def test_every_one_hole_verdict_is_a_proof():
    """No verdict may rest on sampling alone."""
    for record in _audit()["systems"].values():
        for row in record["one_hole_hall_feasible_audit"]:
            verdict = row["determinant"]
            assert "proves" in verdict, row
            if verdict["nonzero"]:
                assert verdict["proves"].startswith("nonvanishing")
            else:
                assert verdict["method"] == "exact_sparse_expansion"
                assert verdict["coefficient_count"] == 0


def test_singleton_incidence_injective_at_every_rank():
    audit = _audit()
    for record in audit["systems"].values():
        assert record["singleton_collision_groups"] == 0
        assert record["singleton_signature_count"] == record["nonstructural_survivors"]
    assert audit["summary"]["singleton_incidence_injective_everywhere"] is True


def test_optimal_fusion_width_grows_by_one_per_rank():
    audit = _audit()
    assert audit["systems"]["(3,7)"]["max_optimal_fusion_width"] == 3
    assert audit["systems"]["(3,8)"]["max_optimal_fusion_width"] == 4
    assert audit["systems"]["(3,9)"]["max_optimal_fusion_width"] == 5


def test_tangent_matrix_is_square_on_trace_survivors():
    """The trace test forces #positive-side determinants == #ascending roots."""
    module = _load("rank9_compression_audit")
    for record in _audit()["systems"].values():
        for row in record["one_hole_hall_feasible_audit"]:
            tau = tuple(row["tau"])
            on, below, roots, matrix = module.tangent_matrix(tau, 3)
            assert len(below) == len(roots) == row["matrix_size"]
            assert len(on) == row["on_variables"]
            assert module.maximum_matching(matrix) == len(below)


def test_verifier_replays():
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / "verify_rank9_compression_audit.py")],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "rank-7 to rank-9 audit verifier: PASS" in result.stdout
