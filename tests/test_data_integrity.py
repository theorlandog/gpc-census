"""Integrity checks for the census dataset. Zero runtime dependencies."""
import json
import pathlib

import pytest

# tests/data is a symlink to ../results/data (the single source of truth). In a
# built sdist the symlink is not shipped but results/data is, so fall back to it.
ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"

EXPECTED_VERTICES = {
    "vertices_3_9.json": 58,
    "vertices_4_9.json": 103,
    "vertices_3_10.json": 113,
    "vertices_4_10.json": 159,
    "vertices_5_10.json": 292,
}

EXPECTED_CENSUS = {
    "census_3_6_results.txt": (4, 0, 0),
    "census_3_7_results.txt": (10, 0, 0),
    "census_4_8_results.txt": (22, 0, 0),
    "census_3_9_results.txt": (58, 20, 1),
    "census_4_9_results.txt": (103, 16, 2),
    "census_3_10_results.txt": (113, 42, 2),
    "census_4_10_results.txt": (159, 25, 2),
    "census_5_10_results.txt": (292, 42, 4),
}


def test_vertex_counts():
    for name, count in EXPECTED_VERTICES.items():
        data = json.loads((ROOT / "vertices" / name).read_text())
        assert len(data) == count, name


def test_census_tallies():
    for name, (total, intf, real) in EXPECTED_CENSUS.items():
        text = (ROOT / "census" / name).read_text()
        rows = [ln for ln in text.splitlines() if ln[:4].strip().isdigit()]
        assert len(rows) == total, name
        assert sum("INTERFERENCE" in r for r in rows) == intf, name
        assert sum("DESIGN-REAL" in r for r in rows) == real, name


def test_constraint_counts():
    data = json.loads((ROOT / "constraints" / "ak_thesis_rank9_10_constraints.json").read_text())
    expected = {"(3,9)": 52, "(4,9)": 60, "(3,10)": 93, "(4,10)": 125, "(5,10)": 161}
    for system, count in expected.items():
        assert len(data[system]) == count, system


def test_fixture_matches_published_dataset():
    published = ROOT.parents[1] / "results" / "data"
    if not published.exists():
        pytest.skip("published dataset not present (sdist build)")
    fixture = {p.relative_to(ROOT): p.read_bytes() for p in ROOT.rglob("*") if p.is_file()}
    dataset = {
        p.relative_to(published): p.read_bytes() for p in published.rglob("*") if p.is_file()
    }
    assert fixture == dataset


def test_shipped_states_are_certified():
    # nothing ships unvalidated: every certified closed form in the atlas must
    # re-pass its check, independently of the solver that produced it
    import json

    from gpc_census.dataset import validate_states
    states = ROOT / "states.jsonl"
    if not states.exists():
        pytest.skip("no states atlas present")
    records = [json.loads(ln) for ln in states.read_text().splitlines() if ln.strip()]
    fails = validate_states(records)
    assert fails == [], f"invalid certified states: {fails}"
