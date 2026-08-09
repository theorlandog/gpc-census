"""Enforce the completeness half of the sync invariant, in both directions.

Forward: every file under ``results/data/`` must carry a SHA256SUMS entry, a
PROVENANCE.md row, a named generating script, and at least one guard test.
An artifact missing any of those is an orphan: a number nobody can regenerate
or check, which is how published data quietly goes stale.

Reverse: every artifact a guard test names must be on disk. A test naming a
missing artifact SKIPS rather than fails, so without this direction a study can
ship its script, its documentation and its guards while the numbers those
guards would check were never committed. That happened: nine guards for the
chemistry study sat inert until its artifact landed, and the forward audit
reported clean throughout.

The audit found 52 orphans on first run, including two defects worth naming:
``census_master.csv`` had no generator, was malformed CSV, and was missing an
entire system block; and two fiber-era artifacts have no surviving generator at
all. This test is committed only after those were fixed or explicitly declared.
"""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "scripts" / "audit_data_completeness.py"

if not MODULE_PATH.exists():  # sdist may omit the audit script
    pytest.skip("completeness audit not present", allow_module_level=True)

SPEC = importlib.util.spec_from_file_location("data_completeness", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


@pytest.fixture(scope="module")
def report():
    if not (ROOT / "results" / "data").exists():
        pytest.skip("published dataset not present (sdist build)")
    return mod.audit()


def test_no_orphan_artifacts(report):
    orphans = report["orphans"]
    detail = "\n".join(
        f"  {o['artifact']}: missing {', '.join(o['missing'])}" for o in orphans
    )
    assert not orphans, f"{len(orphans)} orphan artifact(s):\n{detail}"


def test_the_audit_actually_looked_at_the_dataset(report):
    """A rule table that matched nothing would pass the orphan test vacuously."""
    assert report["artifacts_audited"] > 50
    assert report["checksum_entries"] >= report["artifacts_audited"]


def test_no_family_rule_has_gone_stale(report):
    """Every declared relaxation must still match a real file.

    A rule left behind after its artifacts were renamed or deleted is a hole
    that could silently absorb a future orphan.
    """
    assert report["stale_family_rules"] == []


def test_every_family_rule_declares_why_it_exists():
    """The relaxations are the weak point, so each one is documented in place."""
    for rule in mod.FAMILY_RULES:
        assert rule["why"].strip()
        assert rule["glob"].strip()


def test_the_generating_script_waiver_is_used_only_where_declared():
    """Exactly the two artifacts with no surviving generator may waive it."""
    waived = {
        rule["glob"]
        for rule in mod.FAMILY_RULES
        if rule.get("script") == mod.UNREGENERABLE_NO_GENERATOR
    }
    assert waived == {"cyclic_window_benchmark.json", "entangling_window_hamiltonian.json"}


def test_no_orphaned_guards(report):
    """Reverse direction: a guard naming a missing artifact only skips."""
    orphaned = report["orphaned_guards"]
    detail = "\n".join(
        f"  {o['artifact']}: named by {', '.join(o['referenced_by'])}"
        for o in orphaned
    )
    assert not orphaned, f"{len(orphaned)} orphaned guard(s):\n{detail}"


def test_the_reverse_check_actually_scanned_the_tests(report):
    """A pattern table that matched nothing would pass the check vacuously."""
    assert report["guard_references_checked"] > 40


def test_the_reverse_check_catches_a_missing_artifact(tmp_path, monkeypatch):
    """The failure it exists for, exercised rather than assumed.

    A test file that names an artifact absent from the data directory must be
    reported. Run against temporary directories so the real dataset is never
    touched.
    """
    present, absent = "present.json", "absent.json"
    data = tmp_path / "data"
    tests = tmp_path / "tests"
    data.mkdir()
    tests.mkdir()
    (data / present).write_text("{}")
    # Built by interpolation, never as a literal: this file is itself scanned
    # by the real audit, and a verbatim fixture would report as an orphan here.
    (tests / "test_fake.py").write_text(
        "\n".join('X = DATA / "%s"' % name for name in (present, absent))
    )
    monkeypatch.setattr(mod, "DATA", data)
    monkeypatch.setattr(mod, "TESTS", tests)

    findings = {finding["artifact"] for finding in mod.orphaned_guards()}
    assert findings == {absent}


def test_the_reverse_check_ignores_interpolated_family_names(tmp_path, monkeypatch):
    """f-string paths are the documented families; FAMILY_RULES covers those."""
    data = tmp_path / "data"
    tests = tmp_path / "tests"
    data.mkdir()
    tests.mkdir()
    interpolated = "vertices_{n}_{d}.json"
    (tests / "test_fake.py").write_text(
        'P = DATA / "vertices" / f"%s"' % interpolated
    )
    monkeypatch.setattr(mod, "DATA", data)
    monkeypatch.setattr(mod, "TESTS", tests)

    assert mod.orphaned_guards() == []


def test_guard_reference_exemptions_are_documented():
    """An exemption without a reason is a hole; the table must stay honest."""
    for name, reason in mod.GUARD_REFERENCE_EXEMPTIONS.items():
        assert name.strip()
        assert isinstance(reason, str) and reason.strip()


def test_wrong_directory_same_basename_is_orphaned(tmp_path, monkeypatch):
    """A same-named file elsewhere must not satisfy a directory-qualified ref.

    Resolving by basename would recreate the silent skip this direction exists
    to remove: the guard would look satisfied while checking nothing.
    """
    data = tmp_path / "data"
    tests = tmp_path / "tests"
    (data / "bar").mkdir(parents=True)
    tests.mkdir()
    (data / "bar" / "mine.json").write_text("{}")
    wanted = "foo/mine.json"
    (tests / "test_fake.py").write_text('X = DATA / "%s"' % wanted)
    monkeypatch.setattr(mod, "DATA", data)
    monkeypatch.setattr(mod, "TESTS", tests)

    findings = {f["artifact"]: f["why"] for f in mod.orphaned_guards()}
    assert set(findings) == {wanted}
    assert "same basename" in findings[wanted]


def test_ambiguous_basename_reference_is_orphaned(tmp_path, monkeypatch):
    """Two artifacts sharing a basename make a bare reference ambiguous."""
    data = tmp_path / "data"
    tests = tmp_path / "tests"
    (data / "one").mkdir(parents=True)
    (data / "two").mkdir(parents=True)
    tests.mkdir()
    (data / "one" / "shared.json").write_text("{}")
    (data / "two" / "shared.json").write_text("{}")
    (tests / "test_fake.py").write_text('X = DATA / "%s"' % "shared.json")
    monkeypatch.setattr(mod, "DATA", data)
    monkeypatch.setattr(mod, "TESTS", tests)

    findings = {f["artifact"]: f["why"] for f in mod.orphaned_guards()}
    assert set(findings) == {"shared.json"}
    assert "ambiguous" in findings["shared.json"]


def test_unique_basename_reference_still_resolves(tmp_path, monkeypatch):
    """The legitimate case: exactly one artifact carries that basename."""
    data = tmp_path / "data"
    tests = tmp_path / "tests"
    (data / "only").mkdir(parents=True)
    tests.mkdir()
    (data / "only" / "unique.json").write_text("{}")
    (tests / "test_fake.py").write_text('X = DATA / "%s"' % "unique.json")
    monkeypatch.setattr(mod, "DATA", data)
    monkeypatch.setattr(mod, "TESTS", tests)

    assert mod.orphaned_guards() == []


def test_single_quoted_references_are_scanned(tmp_path, monkeypatch):
    """A guard using single quotes must not slip past the scan."""
    data = tmp_path / "data"
    tests = tmp_path / "tests"
    data.mkdir()
    tests.mkdir()
    absent = "absent.json"
    (tests / "test_fake.py").write_text("X = DATA / %r" % absent)
    monkeypatch.setattr(mod, "DATA", data)
    monkeypatch.setattr(mod, "TESTS", tests)

    assert {f["artifact"] for f in mod.orphaned_guards()} == {absent}
