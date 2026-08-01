"""Regression gate for the dependency-free exact H-to-V verifier."""
from __future__ import annotations

import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
VERIFIER = ROOT / "scripts" / "verify_vertex_exhaustion_standalone.py"
CONSTRAINTS = ROOT / "src" / "gpc_census" / "data" / "constraints.json"
VERTICES = ROOT / "results" / "data" / "vertices"


def test_standalone_h_to_v_verifier_exhausts_all_nine_polytopes():
    result = subprocess.run(
        [
            sys.executable,
            str(VERIFIER),
            "--constraints",
            str(CONSTRAINTS),
            "--vertices",
            str(VERTICES),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "799 tabulated vertices exhaust 9 H-polytopes" in result.stdout
