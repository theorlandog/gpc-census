"""Replay the standalone held-out Levi audit."""

from __future__ import annotations

import pathlib
import subprocess
import sys


def test_standalone_heldout_replay():
    root = pathlib.Path(__file__).resolve().parents[1]
    script = root / "scripts" / "verify_levi_fusion_heldout_standalone.py"
    result = subprocess.run(
        [sys.executable, str(script)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "held-out Levi audit: PASS" in result.stdout
    assert "P1: FAIL at (3,10), (4,10)" in result.stdout
