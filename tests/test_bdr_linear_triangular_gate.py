"""Guards for gate 2: BDR's own linear-triangular filter, run and compared.

Two things are pinned here. First, that the repository's published irredundant
systems agree as sets with the upstream inequalities at seven of the nine census
systems, which is the only external confirmation those systems have. Second,
that the upstream criterion and the local Dulmage-Mendelsohn reading are
incomparable, which refutes a sentence gate 4 shipped without testing.

The external column cannot be recomputed without a SageMath runtime, so it is
read from the recorded artifact. Everything on the local side is recomputed.
"""

from __future__ import annotations

import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

GATE = DATA / "bdr_linear_triangular_gate.json"
EXTERNAL = DATA / "bdr_external_linear_triangular.json"

SYSTEMS = ("(3,7)", "(4,7)", "(3,8)", "(4,8)", "(3,9)", "(4,9)", "(3,10)")

# Upstream carries the ordering row lambda_d >= 0 in the same list from d = 8 on.
STRUCTURAL_ROW = {d: [0] * (d - 1) + [-1] for d in (8, 9, 10)}


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = _load("bdr_linear_triangular_gate")


def _artifact() -> dict:
    return json.loads(GATE.read_text())


def _external() -> dict:
    return json.loads(EXTERNAL.read_text())


def test_the_computed_systems_came_from_the_pipeline():
    """(3,9), (4,9) and (3,10) are not shipped upstream; they were computed."""
    systems = _artifact()["systems"]
    for name in ("(3,9)", "(4,9)", "(3,10)"):
        assert "computed by the upstream pipeline" in systems[name]["external_source"]
    for name in ("(3,7)", "(4,7)", "(3,8)", "(4,8)"):
        assert "shipped reference data" in systems[name]["external_source"]


def test_borland_dennis_is_out_of_scope_upstream():
    """The one census system BDR declines, recorded rather than left silent."""
    out_of_scope = _external()["out_of_scope_upstream"]
    assert set(out_of_scope) == {"(3,6)"}
    assert "Condition C" in out_of_scope["(3,6)"]
    assert "(3,6)" not in _artifact()["systems"]


def test_bdr_triangular_share_collapses_as_systems_grow():
    """BDR validates a handful of rows however large the system gets."""
    systems = _artifact()["systems"]
    assert systems["(3,10)"]["published_rows"] == 93
    assert systems["(3,10)"]["bdr_linear_triangular_rows"] == 3
    for record in systems.values():
        assert record["bdr_linear_triangular_rows"] <= 5


def test_the_backend_was_actually_executed():
    """The correction that motivated this gate: it ran, it was not assumed."""
    backend = _external()["external_backend"]
    assert backend["executed_here"] is True
    assert backend["pinned_revision"] == MODULE.PINNED_REVISION
    assert backend["license"] == "MIT"
    assert backend["entry_point"].endswith("is_linear_triangular")


def test_published_systems_match_the_upstream_reference_sets():
    """Cross-validation: no published row is absent upstream, at seven systems."""
    systems = _artifact()["systems"]
    assert set(systems) == set(SYSTEMS)
    for name in SYSTEMS:
        assert systems[name]["published_only"] == []


def test_the_only_upstream_extra_row_is_the_structural_one():
    systems = _artifact()["systems"]
    assert systems["(3,7)"]["external_only"] == []
    assert systems["(4,7)"]["external_only"] == []
    for name in ("(3,8)", "(4,8)", "(3,9)", "(4,9)", "(3,10)"):
        record = systems[name]
        assert record["external_only"] == [STRUCTURAL_ROW[record["orbitals"]]]


def test_published_row_counts():
    expected = {
        "(3,7)": (4, 4),
        "(4,7)": (4, 4),
        "(3,8)": (31, 32),
        "(4,8)": (15, 16),
        "(3,9)": (52, 53),
        "(4,9)": (60, 61),
        "(3,10)": (93, 94),
    }
    systems = _artifact()["systems"]
    for name, (published, external) in expected.items():
        assert systems[name]["published_rows"] == published
        assert systems[name]["external_rows"] == external
        assert systems[name]["rows_in_both"] == published


def test_the_two_criteria_are_incomparable():
    """The refutation: neither direction of implication survives."""
    artifact = _artifact()
    pooled = artifact["pooled_contingency"]
    assert pooled["rows"] == 259
    assert pooled["dm_triangular_only"] == 16
    assert pooled["bdr_triangular_only"] == 12
    assert pooled["dm_triangular_and_bdr_triangular"] == 9
    assert pooled["neither"] == 222
    assert artifact["criteria_agree"] is False
    assert artifact["dm_triangularity_is_necessary_for_bdr"] is False
    assert artifact["dm_triangularity_is_sufficient_for_bdr"] is False


def test_the_strictly_stronger_claim_is_recorded_as_refuted():
    verdict = _artifact()["verdict"]
    assert verdict["documented_claim_status"].startswith("REFUTED")
    assert "STRICTLY STRONGER" in verdict["documented_claim_under_test"]


def test_per_system_bdr_triangular_counts():
    expected = {
        "(3,7)": 1,
        "(4,7)": 1,
        "(3,8)": 3,
        "(4,8)": 5,
        "(3,9)": 3,
        "(4,9)": 5,
        "(3,10)": 3,
    }
    systems = _artifact()["systems"]
    for name, count in expected.items():
        assert systems[name]["bdr_linear_triangular_rows"] == count


def test_the_rank_seven_witness_of_incomparability():
    """One named row: a single irreducible DM block that BDR still certifies."""
    rows = {tuple(row["tau"]): row for row in _artifact()["systems"]["(3,7)"]["rows"]}
    witness = rows[(-2, 1, 1, 1, 1, -2, -2)]
    assert witness["bdr_linear_triangular"] is True
    assert witness["dm_blocks"] == [4]
    assert witness["dm_fully_triangular"] is False


def test_local_column_recomputes_on_the_witness():
    """Recomputed from tau, not read from the artifact."""
    gate = _load("levi_triangularity_gate")
    record = MODULE.local_column(gate, (-2, 1, 1, 1, 1, -2, -2), 3)
    assert record["local_status"] == "determinant_nonzero"
    assert record["dm_blocks"] == [4]
    assert record["dm_fully_triangular"] is False


def test_published_taus_reproduce_the_artifact():
    """The join key is recomputed from the shipped constraints, not trusted."""
    systems = _artifact()["systems"]
    for name in SYSTEMS:
        record = systems[name]
        taus = set(MODULE.published_taus(record["particle_number"], record["orbitals"]))
        stored = {
            tuple(row["tau"]) for row in record["rows"] if row["in_published_system"]
        }
        assert taus == stored


def test_external_artifact_is_bound_by_hash():
    artifact = _artifact()
    import hashlib

    assert artifact["external_verdicts_sha256"] == hashlib.sha256(
        EXTERNAL.read_bytes()
    ).hexdigest()


def test_structural_row_is_triangular_on_both_readings():
    """A row with no positive weights is vacuously triangular either way."""
    rows = {tuple(row["tau"]): row for row in _artifact()["systems"]["(3,8)"]["rows"]}
    structural = rows[tuple(STRUCTURAL_ROW[8])]
    assert structural["local_status"] == "structural"
    assert structural["bdr_linear_triangular"] is True
    assert structural["in_published_system"] is False
