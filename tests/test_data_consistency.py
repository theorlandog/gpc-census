"""The campaign-artifact gate the Paper 1 manuscript checker released."""
from __future__ import annotations

import json
import pathlib
import shutil
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from check_data_consistency import main  # noqa: E402

DATA = ROOT / "results/data"


def test_shipped_campaign_artifacts_are_consistent():
    assert main([]) == 0


@pytest.fixture
def artifact_copy(tmp_path):
    """A writable copy of the flat artifact files the checker reads."""
    target = tmp_path / "data"
    target.mkdir()
    for path in DATA.glob("*.json"):
        shutil.copy(path, target / path.name)
    shutil.copy(DATA / "states.jsonl", target / "states.jsonl")
    return target


def test_a_gauge_campaign_rerun_regression_is_caught(artifact_copy, capsys):
    """POWER TEST. This is the failure the checker exists for.

    Rerunning scripts/gauge_census.py against its own artifact once reported
    zero improvements and stripped every exact-arithmetic flag. Nothing in the
    Paper 1 checker reads these fields any more, so this gate has to.
    """
    path = artifact_copy / "gauge_min_summary.json"
    summary = json.loads(path.read_text())
    summary["strictly_improved_by_exterior_contractions"] = 0
    summary["newly_certified_by_exterior_contractions"] = 0
    summary["exterior_certificates_exact"] = 0
    path.write_text(json.dumps(summary, indent=2) + "\n")

    assert main(["--data", str(artifact_copy)]) == 1
    report = capsys.readouterr().out
    assert "strictly_improved_by_exterior_contractions" in report
    assert "exterior_certificates_exact" in report


def test_a_dropped_cascade_bound_is_caught(artifact_copy, capsys):
    path = artifact_copy / "cascade_exact.json"
    report = json.loads(path.read_text())
    report["exact_spectrum_only_s_Q_free"] = 69
    path.write_text(json.dumps(report, indent=2) + "\n")

    assert main(["--data", str(artifact_copy)]) == 1
    assert "cascade_exact" in capsys.readouterr().out


def test_a_v103_offdiagonal_regression_is_caught(artifact_copy, capsys):
    """Quoting only one of the two nonzero entries was a review finding."""
    path = artifact_copy / "v103_sqfree_certificate.json"
    certificate = json.loads(path.read_text())
    certificate["rdm_nonzero_offdiagonals"].pop("rho[1,6]", None)
    path.write_text(json.dumps(certificate, indent=2) + "\n")

    assert main(["--data", str(artifact_copy)]) == 1
    assert "rho[1,6]" in capsys.readouterr().out
