"""Guards for the three-way birationality audit.

WHICH DIRECTION CARRIES EVIDENCE. `Is_Ram_contracted` searches for a
non-contraction witness and returns False the moment it finds one, so False is
the exact, certificate-bearing verdict and True only records that the random
search came up empty. An earlier version of this audit asserted the reverse,
which inverted both headlines. These guards pin the corrected reading: the
strong claim is that 228 rows are witnessed NOT birational and none of them is a
published facet; the Monte Carlo claim is separate and labelled as such.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

AUDIT = DATA / "bdr_birationality_audit.json"
EXTERNAL = DATA / "bdr_birationality_external.json"
CANDIDATES = DATA / "bdr_birationality_candidates.json"

SYSTEMS = ("(3,7)", "(3,8)", "(3,9)")
CRITERIA = ("no_obstruction_found", "bdr_linear_triangular", "dm_fully_triangular")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = _load("bdr_birationality_audit")


def _artifact() -> dict:
    return json.loads(AUDIT.read_text())


def test_no_published_facet_is_ever_witnessed_non_birational():
    """The strong direction: False is a witness, and no facet produced one."""
    artifact = _artifact()
    for name in SYSTEMS:
        record = artifact["systems"][name]
        assert record["witnessed_non_birational_facets"] == 0, name
    strong = artifact["strong_direction"]
    assert strong["witnessed_non_birational_facets"] == 0
    assert strong["witnessed_non_birational_rows"] == 228


def test_the_witnessed_rows_are_all_non_facets():
    """228 rows carry an actual obstruction, and every one is a non-facet."""
    artifact = _artifact()
    witnessed = [
        row
        for record in artifact["systems"].values()
        for row in record["rows"]
        if row["non_birational_witnessed"]
    ]
    assert len(witnessed) == 228
    assert not any(row["is_published_facet"] for row in witnessed)


def test_the_monte_carlo_direction_is_labelled_not_proved():
    """87 facets and 109 non-facets survived the search; that is not a proof."""
    artifact = _artifact()
    pooled = artifact["pooled_scores"]["no_obstruction_found"]
    assert pooled["tp"] == 87
    assert pooled["fp"] == 109
    assert pooled["fn"] == 0
    assert pooled["tn"] == 228
    assert pooled["unknown"] == 0
    assert pooled["recall_on_facets"] == 1.0
    verdict = artifact["verdict"]
    assert "EMPIRICALLY SUPPORTED, not proved" in verdict["monte_carlo_result"]
    assert "CANDIDATE claim" in verdict["monte_carlo_result"]


def test_the_evidence_direction_is_recorded_as_false():
    """The correction that motivated this revision must not silently regress."""
    artifact = _artifact()
    assert artifact["verdict"]["which_direction_is_evidence"].startswith("False")
    backend = artifact["external_backend"]
    assert backend["which_direction_is_evidence"].startswith("False")
    assert backend["upstream_default_random_deep"] == 1


def test_sampling_is_explicit_and_deeper_than_the_upstream_default():
    """random_deep must never be left to the upstream default of one trial."""
    artifact = _artifact()
    sampling = artifact["sampling"]
    assert sampling["random_deep"] == MODULE.RANDOM_DEEP >= 8
    assert sampling["seeds"] == list(MODULE.SEEDS)
    assert len(sampling["seeds"]) >= 3
    external = json.loads(EXTERNAL.read_text())
    assert external["external_backend"]["random_deep"] == sampling["random_deep"]
    assert external["external_backend"]["seeds"] == sampling["seeds"]


def test_every_row_carries_a_trial_per_seed():
    external = json.loads(EXTERNAL.read_text())
    for record in external["systems"].values():
        for row in record["rows"]:
            assert len(row["trials"]) == len(MODULE.SEEDS)
            assert {t["seed"] for t in row["trials"]} == set(MODULE.SEEDS)
            assert all(t["random_deep"] == MODULE.RANDOM_DEEP for t in row["trials"])


def test_the_verdicts_are_stable_across_seeds():
    """Zero disagreement is what makes the Monte Carlo half worth reporting."""
    artifact = _artifact()
    assert artifact["strong_direction"]["seed_disagreement_rows"] == 0
    for record in artifact["systems"].values():
        assert record["seed_disagreement_rows"] == 0


def test_the_join_is_bound_to_its_candidate_input():
    """Fail closed: the external half must name the exact candidates it saw."""
    artifact = _artifact()
    candidates_digest = hashlib.sha256(CANDIDATES.read_bytes()).hexdigest()
    assert artifact["candidates_sha256"] == candidates_digest
    external = json.loads(EXTERNAL.read_text())
    assert external["candidates_sha256"] == candidates_digest


def test_binding_rejects_a_mismatched_candidate_file():
    """The guard must actually refuse, not just record a hash."""
    external = json.loads(EXTERNAL.read_text())
    with pytest.raises(SystemExit) as raised:
        MODULE.check_binding(b'{"systems": {}}', external)
    assert "different candidate file" in str(raised.value)


def test_binding_rejects_a_missing_external_row():
    """A dropped row must fail, never be scored as unknown and excluded."""
    candidates_bytes = CANDIDATES.read_bytes()
    external = json.loads(EXTERNAL.read_text())
    external["systems"]["(3,7)"]["rows"] = external["systems"]["(3,7)"]["rows"][1:]
    with pytest.raises(SystemExit) as raised:
        MODULE.check_binding(candidates_bytes, external)
    assert "do not cover the candidates" in str(raised.value)


def test_binding_rejects_an_errored_row():
    candidates_bytes = CANDIDATES.read_bytes()
    external = json.loads(EXTERNAL.read_text())
    external["systems"]["(3,7)"]["rows"][0]["bdr_error"] = "boom"
    with pytest.raises(SystemExit) as raised:
        MODULE.check_binding(candidates_bytes, external)
    assert "errored" in str(raised.value)


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
            assert row["no_obstruction_found"] in (True, False)
            assert row["non_birational_witnessed"] in (True, False)


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
    assert stored["non_birational_witnessed"] is False


def test_the_probabilistic_arm_is_named_in_the_verdict():
    """The caveat must ship with the numbers."""
    artifact = _artifact()
    assert "probabilistic" in artifact["external_backend"]["birationality_arm"]
    assert "Monte Carlo" in artifact["verdict"]["which_direction_is_evidence"]
    assert "exact arm" in artifact["verdict"]["what_would_settle_it"]


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
