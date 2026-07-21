import pathlib
import subprocess
import sys


def test_facet_laws_verify():
    root = pathlib.Path(__file__).resolve().parents[1]
    r = subprocess.run(
        [sys.executable, str(root / "scripts" / "facet_laws.py"), "--verify"],
        capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "rhs law 542/542; edge law 542/542" in r.stdout
