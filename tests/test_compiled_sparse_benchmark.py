import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import compiled_sparse_benchmark as cb  # noqa: E402
import multichart_symplectic_propagator as mp  # noqa: E402

DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "compiled_sparse_benchmark.json"


def test_compiled_rhs_matches_reference_rhs():
    atlas = mp.Atlas.build()
    c, _ = mp.build_initial_state(atlas.H)
    y = atlas.pack_from_X(np.outer(np.conj(c), c))
    compiled = cb.CompiledSparseRHS(atlas)
    assert np.max(np.abs(compiled(0.0, y) - atlas.rhs(0.0, y))) < 1e-12


def test_compiled_rhs_agrees_along_a_trajectory():
    """One point could agree by accident; check the whole propagated path.

    Includes the amplitude zero at t = pi/4, where the compiled kernel's
    on-demand coherence formula and the reference reconstruction have to make
    the same boundary choice.
    """
    atlas = mp.Atlas.build()
    c0, _ = mp.build_initial_state(atlas.H)
    compiled = cb.CompiledSparseRHS(atlas)
    for t in (0.0, 0.3, np.pi / 4, 1.0, 2.0, 3 * np.pi / 4, 3.0):
        c = mp.expm(-1j * atlas.H * t) @ c0
        y = atlas.pack_from_X(np.outer(np.conj(c), c))
        assert np.max(np.abs(compiled(0.0, y) - atlas.rhs(0.0, y))) < 1e-11, t


def test_compiled_run_is_accurate():
    row = cb.run_one("compiled_sparse2rdm")
    assert row["error"] < 1e-9
    assert row["nfev"] > 0


def test_artifact_records_equal_work_and_equal_accuracy():
    """The speedup is only meaningful if both methods did the same work."""
    data = json.loads(ARTIFACT.read_text())
    rows = {row["method"]: row for row in data["summary"]}
    assert rows["sparse2rdm"]["nfev"] == rows["compiled_sparse2rdm"]["nfev"] == 314
    assert rows["sparse2rdm"]["steps"] == rows["compiled_sparse2rdm"]["steps"] == 26
    assert rows["sparse2rdm"]["variables"] == rows["compiled_sparse2rdm"]["variables"] == 167
    assert rows["sparse2rdm"]["max_error"] == rows["compiled_sparse2rdm"]["max_error"]
    assert rows["compiled_sparse2rdm"]["fallbacks"] == 0
    assert data["speedup_compiled_over_original"] > 1.0
