"""Regression guard for the v96 signed-design enumerator.

The full exhaustive run (116,916 tail covers, 350,980 skeletons, 0 hits, ~40 s)
is the theorem-grade verification of "v96 admits no signed design at m=1"
(docs/RESEARCH.md, "v96 campaign"); reproduce it with
`python scripts/signed_design_v96_full.py --no-phases`. This test runs a fast
slice to guard the tail-cover / head-solve machinery against regressions.
"""
import pathlib
import subprocess
import sys


def test_v96_signed_design_slice_no_hit():
    root = pathlib.Path(__file__).resolve().parents[1]
    r = subprocess.run(
        [sys.executable, str(root / "scripts" / "signed_design_v96_full.py"),
         "--no-phases", "--limit", "5000"],
        capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "HIT:" not in r.stdout
    assert "hits 0" in r.stdout
