"""Keep the level-5 companion note in step with the artifacts it quotes.

The note states numbers in prose and in LaTeX tables, which the doc-table
emitter cannot own without mangling the manuscript. So it is guarded the other
way round: every load-bearing number in the note is re-derived here from the
artifact and asserted to appear in the text. A regenerated artifact that moves
one of them fails this test rather than silently contradicting the note.

Citation discipline is checked too. The build fails on an unresolved citation,
but it does not catch a bib entry that nothing cites, which the house rule
forbids.
"""

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
NOTE = ROOT / "results" / "report" / "level5" / "main.md"
BIB = ROOT / "results" / "report" / "level5" / "references.bib"
DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"

if not NOTE.exists():
    pytest.skip("level-5 note not present", allow_module_level=True)


@pytest.fixture(scope="module")
def text():
    return NOTE.read_text()


@pytest.fixture(scope="module")
def artifact():
    path = DATA / "level5_ressayre_3_11.json"
    if not path.exists():
        pytest.skip("level5 artifact not present")
    return json.loads(path.read_text())


def grouped(value: int) -> list[str]:
    """The note groups long integers, so accept both renderings."""
    plain = str(value)
    return [plain, f"{value:,}".replace(",", "\\,")]


def test_every_certified_determinant_is_quoted(text, artifact):
    rows = artifact["rank_eleven"]["rows"]
    assert len(rows) == 4
    for row in rows:
        value = row["certificate"]["determinant_value"]
        assert any(form in text for form in grouped(value)), (
            f"{row['family']}: determinant {value} is not quoted in the note"
        )


def test_the_trace_condition_numbers_are_quoted(text, artifact):
    for row in artifact["rank_eleven"]["rows"]:
        assert row["below_weight_count"] == row["negative_root_count"] == 10
        assert len(row["certificate"]["on_indices"]) == 60
    assert "60" in text and "165" in text


def test_the_calibration_table_matches_the_artifact(text, artifact):
    for entry in artifact["calibration"]:
        assert str(entry["certified_rows"]) in text
        assert entry["certified_rows_invalid_on_published_polytope"] == 0
        assert str(entry["certified_verbatim_in_published_table"]) in text
    assert "984" in text
    assert "1{,}914" in text or "1,914" in text


def test_the_negative_control_totals_match(text, artifact):
    false_total = sum(
        e["negative_control"]["demonstrably_false_rows_screened"]
        for e in artifact["calibration"]
    )
    assert false_total == 1857
    assert "1,857" in text
    reached = sum(
        e["negative_control"]["false_rows_reaching_the_tangent_determinant"]
        for e in artifact["calibration"]
    )
    assert reached == 6
    assert "six false rows" in text
    for entry in artifact["calibration"]:
        assert entry["negative_control"]["false_rows_certified"] == 0


def test_the_calibration_screened_row_count_is_right(text, artifact):
    """The note's 8,898 is the CALIBRATION total: 984 + 1,914 + 6,000.

    The four rank-11 rows are excluded on purpose. They are the certified
    result, not part of the control that certified nothing false, so folding
    them in would overstate what the calibration covered.
    """
    total = sum(e["rows_screened"] for e in artifact["calibration"])
    total += artifact["randomized_sweep"]["rows_screened"]
    assert total == 8898, total
    assert "8,898" in text
    assert len(artifact["rank_eleven"]["rows"]) == 4


def test_the_refutation_table_matches_the_artifact(text, artifact):
    for record in artifact["rank_eleven"]["open_candidates"]:
        assert record["refuted"]
        value = record["cut_by"][0]["value"]
        assert value in text, f"candidate {record['index']}: {value} not quoted"
        for cut in record["cut_by"]:
            assert cut["family"] in text


def test_the_note_states_what_is_not_proved(text):
    """The value of the note depends on it not overreaching."""
    lowered = text.lower()
    assert "facet" in lowered
    assert "irredundan" in lowered
    assert "complete" in lowered
    assert "validity is what is proved" in lowered or "validity is proved" in lowered


def test_the_note_records_the_weak_point_of_its_own_calibration(text):
    """The determinant step is the least stressed part; the note must say so."""
    assert "least" in text.lower()
    assert "trace condition does nearly all" in text.lower()


def test_no_attained_spectrum_violates_a_certified_row(artifact):
    eleven = artifact["rank_eleven"]
    assert eleven["attained_vertices_checked"] == 19
    assert eleven["certified_rows_violated_by_an_attained_vertex"] == 0


def test_every_bibliography_entry_is_cited(text):
    """House rule: no nocite, so an uncited entry is a defect."""
    keys = set(re.findall(r"^@\w+\{([^,]+),", BIB.read_text(), re.M))
    assert keys, "bibliography is empty"
    cited = set(re.findall(r"@([A-Za-z][A-Za-z0-9_:.-]*)", text))
    uncited = {k for k in keys if k not in cited}
    assert not uncited, f"bibliography entries never cited: {sorted(uncited)}"


def test_every_citation_resolves_to_a_bibliography_entry(text):
    keys = set(re.findall(r"^@\w+\{([^,]+),", BIB.read_text(), re.M))
    body = text.split("---", 2)[-1]
    cited = set(re.findall(r"\[@([A-Za-z][A-Za-z0-9_:.-]*)", body))
    cited |= set(re.findall(r";\s*@([A-Za-z][A-Za-z0-9_:.-]*)", body))
    missing = {c for c in cited if c not in keys}
    assert not missing, f"citations with no bibliography entry: {sorted(missing)}"


def test_the_note_uses_no_em_dashes():
    """House rule."""
    assert "—" not in NOTE.read_text()
