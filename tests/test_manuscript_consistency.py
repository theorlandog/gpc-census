"""The paper, the dataset, and the checksums cannot drift apart.

Wires scripts/check_manuscript_counts.py (the manuscript-dataset consistency
test) and the released SHA256SUMS into the default suite, so `make test` and
CI fail on any drift between results/report/main.md, results/data, and the
recorded checksums.
"""
import hashlib
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
DATA = REPO / "results" / "data"
SCRIPT = REPO / "scripts" / "check_manuscript_counts.py"


def test_manuscript_counts_match_dataset():
    if not (SCRIPT.exists() and (DATA / "states.jsonl").exists()):
        pytest.skip("source checkout artifacts not present (sdist build)")
    proc = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_sha256sums_match_data_files():
    sums = DATA / "SHA256SUMS"
    if not sums.exists():
        pytest.skip("no SHA256SUMS present (sdist build)")
    listed = {}
    for line in sums.read_text().splitlines():
        digest, name = line.split(None, 1)
        listed[name.strip().lstrip("./")] = digest
    on_disk = {
        str(p.relative_to(DATA)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in DATA.rglob("*")
        if p.is_file() and p.name != "SHA256SUMS"
    }
    assert on_disk == listed
