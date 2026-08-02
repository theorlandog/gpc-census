import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import exact_observable_benchmark as ex  # noqa: E402
import multichart_symplectic_propagator as mp  # noqa: E402

DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "exact_observable_benchmark.json"
CONTROL = DATA / "exact_observable_benchmark_with_control.json"


def test_exact_observable_matches_direct():
    assert ex.verify() < 1e-13


def test_analytic_rotation_is_exact_at_many_times():
    """One time could agree by accident; the rotation must be exact throughout.

    Includes both designed amplitude zeros and a time past the benchmark
    interval, so the agreement is not a property of t = 1.
    """
    atlas, c0, _X0, y0 = ex.build_case()
    prop = ex.ExactObservablePropagator(atlas)
    for t in (0.0, 0.4, np.pi / 4, 1.0, 2.0, 3 * np.pi / 4, 3.0, 7.5):
        X, _ = atlas.reconstruct_X(prop.propagate(y0, t))
        c = mp.expm(-1j * atlas.H * t) @ c0
        assert np.max(np.abs(X - np.outer(np.conj(c), c))) < 1e-13, t


def test_hamiltonian_really_is_four_disjoint_blocks_and_spectators():
    atlas, _c0, _X0, _y0 = ex.build_case()
    blocks = ex.ExactObservablePropagator(atlas).blocks
    assert len(blocks) == 4
    covered = {i for i, _j, _s in blocks} | {j for _i, j, _s in blocks}
    assert len(covered) == 8
    assert np.count_nonzero(atlas.H) == 8
    assert len(atlas.support) - len(covered) == 21
    for i, j, sign in blocks:
        assert atlas.H[i, j] == atlas.H[j, i] == sign


def test_control_update_is_exact_too():
    """The control has to be a real opponent, not a deliberately weak one."""
    atlas, c0, _X0, _y0 = ex.build_case()
    blocks = ex.ExactObservablePropagator(atlas).blocks
    updated = ex.specialized_carrier_update(blocks, c0, 1.0)
    reference = mp.expm(-1j * atlas.H * 1.0) @ c0
    assert np.max(np.abs(updated - reference)) < 1e-14


def test_control_artifact_reports_the_observable_method_losing():
    control = json.loads(CONTROL.read_text())
    assert control["specialized_carrier_speedup_over_observable"] > 1.0
    assert control["exact_observable_error"] < 1e-14
    assert "314" in control["work_accounting"]


def test_speedups_are_the_quotients_they_claim_to_be():
    art = json.loads(ARTIFACT.read_text())
    median = art["method"]["median_elapsed_seconds"]
    for key, method in (
        ("vs_carrier", "carrier"),
        ("vs_fci", "fci"),
        ("vs_dense2rdm", "dense2rdm"),
        ("vs_original_sparse2rdm", "sparse2rdm"),
    ):
        assert art["speedups"][key] == art["baseline_medians"][method] / median
    assert art["speedups"]["vs_compiled_sparse2rdm"] == art["compiled_sparse_median"] / median
    assert art["method"]["real_variables"] == 167
