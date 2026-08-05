"""Enforce the completeness half of the sync invariant.

Every file under ``results/data/`` must carry a SHA256SUMS entry, a
PROVENANCE.md row, a named generating script, and at least one guard test.
An artifact missing any of those is an orphan: a number nobody can regenerate
or check, which is how published data quietly goes stale.

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
