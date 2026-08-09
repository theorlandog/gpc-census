"""Guards for the three-way birationality audit.

The load-bearing claim is the zero-false-negative one: every published facet
passes BDR's birationality test. The second claim is its converse failing, which
is what keeps redundancy elimination irreducible. Both are pinned here, and the
local column is recomputed rather than read back.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

AUDIT = DATA / "bdr_birationality_audit.json"
EXTERNAL = DATA / "bdr_birationality_external.json"
CANDIDATES = DATA / "bdr_birationality_candidates.json"

SYSTEMS = ("(3,7)", "(3,8)", "(3,9)")
CRITERIA = ("bdr_birational", "bdr_linear_triangular", "dm_fully_triangular")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = _load("bdr_birationality_audit")


def _artifact() -> dict:
    return json.loads(AUDIT.read_text())


def test_every_published_facet_is_birational():
    """The headline: zero false negatives, on every system and pooled."""
    artifact = _artifact()
    for name in SYSTEMS:
        score = artifact["systems"][name]["scores"]["bdr_birational"]
        assert score["fn"] == 0, name
        assert score["recall_on_facets"] == 1.0, name
    pooled = artifact["pooled_scores"]["bdr_birational"]
    assert pooled["fn"] == 0
    assert pooled["recall_on_facets"] == 1.0
    assert pooled["tp"] == 87


def test_birationality_is_not_sufficient():
    """If it were, redundancy elimination would not be the irreducible step."""
    pooled = _artifact()["pooled_scores"]["bdr_birational"]
    assert pooled["fp"] == 109
    assert 0.4 < pooled["precision"] < 0.5
    assert pooled["separates_facets_exactly"] is False


def test_precision_climbs_with_rank_while_recall_holds():
    systems = _artifact()["systems"]
    precisions = [systems[name]["scores"]["bdr_birational"]["precision"] for name in SYSTEMS]
    assert precisions[0] < precisions[1] < precisions[2]
    for name in SYSTEMS:
        assert systems[name]["scores"]["bdr_birational"]["recall_on_facets"] == 1.0


def test_both_triangularity_readings_are_poor_facet_tests():
    """Each misses roughly nine facets in ten."""
    pooled = _artifact()["pooled_scores"]
    for column in ("bdr_linear_triangular", "dm_fully_triangular"):
        score = pooled[column]
        assert score["fn"] > 70, column
        assert score["recall_on_facets"] < 0.2, column
        assert score["precision"] < 0.4, column


def test_birationality_beats_both_triangularity_readings():
    pooled = _artifact()["pooled_scores"]
    birational = pooled["bdr_birational"]
    for column in ("bdr_linear_triangular", "dm_fully_triangular"):
        assert birational["recall_on_facets"] > pooled[column]["recall_on_facets"]
        assert birational["precision"] > pooled[column]["precision"]


def test_every_published_facet_is_in_the_candidate_population():
    """Otherwise the scoring would be distorted by facets never enumerated."""
    artifact = _artifact()
    expected = {"(3,7)": (48, 4), "(3,8)": (136, 31), "(3,9)": (240, 52)}
    assert set(artifact["systems"]) == set(SYSTEMS)
    for name, (candidates, facets) in expected.items():
        record = artifact["systems"][name]
        assert record["candidate_rows"] == candidates
        assert record["published_facets_among_them"] == facets
        assert record["published_facets_not_candidates"] == []


def test_no_row_errored_in_the_external_pass():
    external = json.loads(EXTERNAL.read_text())
    for record in external["systems"].values():
        assert record["errored_rows"] == 0
        for row in record["rows"]:
            assert row["bdr_error"] is None
            assert row["bdr_birational"] in (True, False)


def test_scores_recompute_from_the_rows():
    """The stored contingency must match the rows it claims to summarise."""
    artifact = _artifact()
    for name in SYSTEMS:
        record = artifact["systems"][name]
        for column in CRITERIA:
            assert MODULE._score(record["rows"], column) == record["scores"][column]


def test_the_local_column_recomputes_on_a_named_row():
    """Recomputed from tau, not read from the artifact."""
    gate = _load("levi_triangularity_gate")
    rows = {tuple(r["tau"]): r for r in _artifact()["systems"]["(3,7)"]["rows"]}
    tau = (-2, 1, 1, 1, 1, -2, -2)
    stored = rows[tau]
    _on, below, _roots, matrix = gate.AUDIT.tangent_matrix(tau, 3)
    match_row = gate.perfect_matching(matrix)
    assert match_row is not None
    assert gate.irreducible_blocks(matrix, match_row) == stored["dm_blocks"]
    assert stored["is_published_facet"] is True
    assert stored["bdr_birational"] is True


def test_the_probabilistic_arm_is_named_in_the_verdict():
    """The caveat must ship with the numbers."""
    artifact = _artifact()
    assert "probabilistic" in artifact["external_backend"]["birationality_arm"]
    assert "probabilistic" in artifact["verdict"]["arm_caveat"]
    assert "weaker evidence" in artifact["verdict"]["arm_caveat"]


def test_external_artifact_is_bound_by_hash():
    artifact = _artifact()
    assert artifact["external_verdicts_sha256"] == hashlib.sha256(
        EXTERNAL.read_bytes()
    ).hexdigest()


def test_candidate_and_audit_populations_agree():
    candidates = json.loads(CANDIDATES.read_text())
    artifact = _artifact()
    for name in SYSTEMS:
        stored = {tuple(r["tau"]) for r in candidates["systems"][name]["rows"]}
        joined = {tuple(r["tau"]) for r in artifact["systems"][name]["rows"]}
        assert stored == joined
