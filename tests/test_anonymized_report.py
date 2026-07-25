"""The anonymized review copy cannot drift from the master manuscript.

results/report/anonymized/ is generated from results/report/main.md by
scripts/make_anonymized_report.py (see `make anonymize`). These tests fail
if the committed copy is stale relative to the master, or if any
author-identifying string survives anonymization.
"""
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "make_anonymized_report.py"
ANON = REPO / "results" / "report" / "anonymized"
SHIPPED = ("main.md", "references.bib", "institute-of-physics-numeric.csl")


def test_anonymized_copy_is_current(tmp_path):
    if not (SCRIPT.exists() and ANON.exists()):
        pytest.skip("source checkout artifacts not present (sdist build)")
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(tmp_path)],
        cwd=REPO, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    for name in SHIPPED:
        regen = (tmp_path / name).read_bytes()
        committed = (ANON / name).read_bytes()
        assert regen == committed, f"{name} is stale: run make anonymize"


def test_anonymized_copy_has_no_identifying_strings():
    if not ANON.exists():
        pytest.skip("no anonymized copy present (sdist build)")
    text = (ANON / "main.md").read_text().lower()
    for leak in ("orlando", "jamie", "orcid", "zenodo", "derry", "we thank"):
        assert leak not in text, leak
