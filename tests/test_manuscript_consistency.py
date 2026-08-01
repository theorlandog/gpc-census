"""The paper, the dataset, and the checksums cannot drift apart.

Wires scripts/check_manuscript_counts.py (the manuscript-dataset consistency
test) and the released SHA256SUMS into the default suite, so `make test` and
CI fail on any drift between results/report/main.md, results/data, and the
recorded checksums.
"""
import hashlib
import pathlib
import runpy
import subprocess
import sys

import pytest
import sympy as sp

REPO = pathlib.Path(__file__).resolve().parents[1]
DATA = REPO / "results" / "data"
SCRIPT = REPO / "scripts" / "check_manuscript_counts.py"
STATE_VERIFIER = REPO / "scripts" / "verify_states_standalone.py"
MANUSCRIPT = REPO / "results" / "report" / "main.md"


def test_manuscript_counts_match_dataset():
    if not (SCRIPT.exists() and (DATA / "states.jsonl").exists()):
        pytest.skip("source checkout artifacts not present (sdist build)")
    proc = subprocess.run([sys.executable, str(SCRIPT)], cwd=REPO,
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_standalone_state_verifier_accepts_every_record():
    if not STATE_VERIFIER.exists():
        pytest.skip("standalone state verifier not present (sdist build)")
    proc = subprocess.run(
        [sys.executable, str(STATE_VERIFIER), str(DATA / "states.jsonl")],
        cwd=REPO,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "checked 799 records" in proc.stdout
    assert "all 799 records PASS" in proc.stdout


def test_standalone_state_verifier_rejects_an_incomplete_atlas(tmp_path):
    if not STATE_VERIFIER.exists():
        pytest.skip("standalone state verifier not present (sdist build)")
    first_record = (DATA / "states.jsonl").read_text().splitlines()[0] + "\n"
    subset = tmp_path / "states.jsonl"
    subset.write_text(first_record)
    proc = subprocess.run(
        [sys.executable, str(STATE_VERIFIER), str(subset)],
        cwd=REPO,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    assert proc.returncode == 1
    assert "expected 799 records" in proc.stdout


def test_vb_off_diagonal_pins_the_paper_rdm_convention():
    if not STATE_VERIFIER.exists():
        pytest.skip("standalone state verifier not present (sdist build)")
    build_rdm = runpy.run_path(str(STATE_VERIFIER))["build_rdm"]
    phase = (3 + sp.I * sp.sqrt(215)) / (4 * sp.sqrt(14))
    determinants = [
        (0, 1, 2, 5),
        (0, 1, 3, 7),
        (0, 1, 3, 8),
        (0, 2, 3, 4),
        (0, 2, 6, 7),
        (0, 2, 6, 8),
        (1, 2, 4, 7),
        (2, 3, 7, 8),
    ]
    amplitudes = [
        2,
        sp.sqrt(2),
        sp.sqrt(7) * phase,
        -sp.sqrt(3),
        sp.sqrt(2),
        sp.sqrt(2),
        -1,
        sp.sqrt(2),
    ]
    amplitudes = [value / sp.sqrt(23) for value in amplitudes]
    rho = build_rdm(determinants, amplitudes, 9)
    expected = (11 + sp.I * sp.sqrt(215)) / 92
    assert sp.simplify(rho[7, 8] - expected) == 0
    assert sp.simplify(rho[8, 7] - sp.conjugate(expected)) == 0

    manuscript = MANUSCRIPT.read_text()
    assert (
        r"\rho^{(1)}_{ij}(\psi)=\langle\psi|a_i^\dagger a_j|\psi\rangle"
        in manuscript
    )
    assert r"=\frac{11+i\sqrt{215}}{92}" in manuscript


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
