"""Guards for gate 4: linear triangularity of fermionic tangent systems."""

from __future__ import annotations

import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

LEVI_TRIANGULARITY_GATE = DATA / "levi_triangularity_gate.json"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = _load("levi_triangularity_gate")


def _artifact() -> dict:
    return json.loads(LEVI_TRIANGULARITY_GATE.read_text())


def test_published_facets_are_mostly_not_triangular():
    """The load-bearing population: the actual irredundant facet system."""
    expected = {"(3,7)": (4, 0, 4), "(3,8)": (31, 5, 7), "(3,9)": (52, 5, 7)}
    systems = _artifact()["systems"]
    for system, (rows, triangular, largest) in expected.items():
        published = systems[system]["published_facets"]
        assert published["analysed_rows"] == rows
        assert published["fully_triangular_rows"] == triangular
        assert published["largest_irreducible_block"] == largest


def test_published_rank_seven_block_decompositions():
    """All four (3,7) facets fail strict triangularity."""
    published = _artifact()["systems"]["(3,7)"]["published_facets"]
    blocks = sorted(tuple(row["blocks"]) for row in published["blocks_by_row"])
    assert blocks == [(3, 1), (4,), (4,), (4,)]
    assert published["fully_triangular_rows"] == 0


def test_candidate_population_is_not_the_facet_system():
    """The two populations must never be conflated: 48/136/240 against 4/31/52."""
    systems = _artifact()["systems"]
    for system, candidates, published in (
        ("(3,7)", 48, 4),
        ("(3,8)", 136, 31),
        ("(3,9)", 240, 52),
    ):
        record = systems[system]
        assert record["determinant_nonzero_hall_feasible_rows"] == candidates
        assert record["published_facets"]["analysed_rows"] == published
        assert candidates > published


def test_the_largest_candidate_witness_is_not_called_a_facet():
    witness = _artifact()["systems"]["(3,9)"]["largest_candidate_block_witness"]
    assert witness["blocks"] == [15, 1]
    assert "NOT a published facet" in witness["status"]
    published = _artifact()["systems"]["(3,9)"]["published_facets"]
    assert published["largest_irreducible_block"] < max(witness["blocks"])


def test_strict_triangularity_is_a_minority_at_every_rank():
    expected = {"(3,7)": (48, 14), "(3,8)": (136, 20), "(3,9)": (240, 23)}
    systems = _artifact()["systems"]
    assert set(systems) == set(expected)
    for system, (rows, triangular) in expected.items():
        record = systems[system]
        assert record["determinant_nonzero_hall_feasible_rows"] == rows
        assert record["fully_triangular_rows"] == triangular
        assert record["fully_triangular_fraction"] < 0.5


def test_the_triangular_fraction_falls_with_rank():
    systems = _artifact()["systems"]
    fractions = [
        systems[system]["fully_triangular_fraction"]
        for system in ("(3,7)", "(3,8)", "(3,9)")
    ]
    assert fractions[0] > fractions[1] > fractions[2]


def test_largest_irreducible_blocks():
    systems = _artifact()["systems"]
    assert systems["(3,7)"]["largest_irreducible_block"] == 9
    assert systems["(3,8)"]["largest_irreducible_block"] == 14
    assert systems["(3,9)"]["largest_irreducible_block"] == 15


def test_rank_nine_witness_is_a_single_large_block():
    witness = _artifact()["systems"]["(3,9)"]["largest_candidate_block_witness"]
    assert witness["tau"] == [-1, -1, -1, -1, -1, 2, 2, -4, 2]
    assert witness["matrix_order"] == 16
    assert witness["blocks"] == [15, 1]
    assert sum(witness["blocks"]) == witness["matrix_order"]


def test_block_decomposition_is_matching_independent():
    """The verdict rests on canonicality, so it is checked, not assumed."""
    for record in _artifact()["systems"].values():
        assert record["matching_independence_rows_checked"] > 0
        assert record["matching_independence_disagreements"] == 0


def test_blocks_recompute_on_the_witness():
    """Recomputed from tau, not read from the artifact."""
    audit = _load("rank9_compression_audit")
    tau = (-1, -1, -1, -1, -1, 2, 2, -4, 2)
    _on, below, roots, matrix = audit.tangent_matrix(tau, 3)
    assert len(below) == len(roots) == 16
    match_row = MODULE.perfect_matching(matrix)
    assert match_row is not None
    assert MODULE.irreducible_blocks(matrix, match_row) == [15, 1]


def test_blocks_partition_the_matrix_order():
    """A sanity property of any decomposition into irreducible blocks."""
    audit = _load("rank9_compression_audit")
    for tau in [
        (-1, -1, -1, -1, -1, 2, 2, -4, 2),
        (2, 0, 0, 0, -1, -1, -1, -2),
    ]:
        _on, below, _roots, matrix = audit.tangent_matrix(tau, 3)
        match_row = MODULE.perfect_matching(matrix)
        if match_row is None:
            continue
        assert sum(MODULE.irreducible_blocks(matrix, match_row)) == len(below)


def test_verdict_records_the_refutation():
    verdict = _artifact()["verdict"]
    assert verdict["strict_triangularity"].startswith("REFUTED")
    assert "open" in verdict["block_triangular_reading"]
    assert "NOT the set of final facets" in verdict["scope"]
