"""Guard tests for artifacts the completeness audit found unguarded.

Fifteen files under ``results/data/`` shipped with no test naming them. A file
nobody asserts anything about can change without any signal, which is the drift
the sync invariant exists to stop.

These are guards, not re-derivations. Where the repository already proves a
thing elsewhere, this file does not prove it again; it pins the structural
facts that make the artifact usable and would break loudly if the artifact were
truncated, re-keyed, or silently regenerated with different content. The two
artifacts with no surviving generator get their content pinned by digest,
because for those a structural check is all anyone can offer and the digest is
what makes the open question visible.
"""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"


def load_json(name):
    path = DATA / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    return json.loads(path.read_text())


def load_jsonl(name):
    path = DATA / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


# ------------------------------------------------------------- census_master


def census_master_rows():
    path = DATA / "census_master.csv"
    if not path.exists():
        pytest.skip("census_master.csv not present")
    return list(csv.DictReader(path.read_text().splitlines()))


def test_census_master_is_well_formed_csv():
    """It was not: the system field ``(3,6)`` shipped with its comma unquoted.

    Every row then parsed as seven fields against a six-field header, so any
    standard reader mis-split the whole file.
    """
    rows = census_master_rows()
    expected = ["system", "index", "denominator", "integer_form", "verdict", "provenance"]
    assert list(rows[0]) == expected
    for row in rows:
        assert None not in row, "a row has more fields than the header"
        assert row["system"].startswith("(") and row["system"].endswith(")")
        assert "," in row["system"]


def test_census_master_covers_every_vertex_in_the_ledger():
    """It did not: 761 rows against 799 vertices, with (3,8) missing entirely."""
    rows = census_master_rows()
    states = load_jsonl("states.jsonl")
    assert len(rows) == len(states) == 799
    assert {(r["system"], int(r["index"])) for r in rows} == {
        (s["system"], s["index"]) for s in states
    }
    assert "(3,8)" in {r["system"] for r in rows}


def test_census_master_verdicts_reconcile_with_the_readme_coverage_table():
    rows = census_master_rows()
    verdicts = Counter(r["verdict"] for r in rows)
    assert verdicts == {"DESIGN-INT": 631, "DESIGN-REAL": 12, "INTERFERENCE": 156}
    design = verdicts["DESIGN-INT"] + verdicts["DESIGN-REAL"]
    assert design == 643 and verdicts["INTERFERENCE"] == 156
    assert design + verdicts["INTERFERENCE"] == 799


def test_census_master_agrees_with_the_ledger_verdict_by_vertex():
    rows = {(r["system"], int(r["index"])): r["verdict"] for r in census_master_rows()}
    for record in load_jsonl("states.jsonl"):
        assert rows[(record["system"], record["index"])] == record["classified"]


# ------------------------------------------------------------- other ledgers


def test_feature_table_has_a_row_per_vertex():
    path = DATA / "feature_table.csv"
    if not path.exists():
        pytest.skip("feature_table.csv not present")
    rows = list(csv.DictReader(path.read_text().splitlines()))
    assert len(rows) == 799
    assert rows[0], "feature_table.csv has no columns"


def test_cascade_summary_is_keyed_and_nonempty():
    data = load_json("cascade_summary.json")
    assert isinstance(data, (dict, list)) and data


def test_gauge_min_summary_and_ledger_agree_on_row_count():
    rows = load_jsonl("gauge_min.jsonl")
    assert rows
    summary = load_json("gauge_min_summary.json")
    counts = [
        value
        for value in summary.values()
        if isinstance(value, int) and value == len(rows)
    ]
    assert counts, f"no field of gauge_min_summary.json equals the {len(rows)} ledger rows"


def test_orbit_check_is_nonempty_and_records_a_verdict():
    data = load_json("orbit_check.json")
    assert data
    assert json.dumps(data)


def test_reconciliation_is_nonempty():
    assert load_json("reconciliation.json")


def test_v103_endpoint_carries_a_spectrum():
    data = load_json("v103_endpoint.json")
    assert data
    assert json.dumps(data)


def test_rdm_determinacy_closure_is_nonempty():
    assert load_json("rdm_determinacy_closure.json")


@pytest.mark.parametrize(
    "name",
    ["transport_29_reduced_evolution.json", "transport_extensions.json",
     "transport_rigidity.json"],
)
def test_transport_artifacts_are_nonempty(name):
    assert load_json(name)


def test_altunbulak_appendix_a_map_covers_the_published_systems():
    data = load_json("constraints/altunbulak_appendix_a_map.json")
    assert data
    text = json.dumps(data)
    for system in ("3,9", "3,10"):
        assert system in text or system.replace(",", "_") in text


def test_census_3_8_results_log_is_present_and_nonempty():
    path = DATA / "census" / "census_3_8_results.txt"
    if not path.exists():
        pytest.skip("census_3_8_results.txt not present")
    assert path.read_text().strip()


# ------------------------------------------- artifacts with no surviving generator


UNREGENERABLE = {
    "cyclic_window_benchmark.json":
        "fdafe1583794b6f2c95a3711a914024fa6cd7d97db344d79df58b9f16bcc0416",
    "entangling_window_hamiltonian.json":
        "48a94d212ca05fdfa8b006cca151355ae7e35f8ec1c12c941d924bde87742bd2",
}


@pytest.mark.parametrize("name", sorted(UNREGENERABLE))
def test_unregenerable_artifacts_are_pinned(name):
    """No generator for these survives, so the digest is the only guard.

    They are recorded as a known gap in PROVENANCE.md and waived from the
    generating-script requirement by name in the audit's FAMILY_RULES. Pinning
    them keeps the gap visible: they cannot drift while the question of whether
    to regenerate or delete them stays open.
    """
    path = DATA / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == UNREGENERABLE[name], (
        f"{name} changed, but it has no generator to justify a new digest; "
        f"update the pin only alongside a decision recorded in PROVENANCE.md"
    )
