import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "generalized_observable_benchmark.json"

MODULE_PATH = ROOT / "scripts" / "generalized_observable_benchmark.py"
SPEC = importlib.util.spec_from_file_location("gob", MODULE_PATH)
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


def test_exact_bloch_propagation():
    cells = m.initial_cells(9)
    b = m.cells_to_bloch(cells)
    w = 0.7 + 0.3 * np.arange(9) / 8
    got = m.propagate_observable(b, w, 1.3)
    ref = m.cells_to_bloch(m.propagate_factorized_amplitudes(cells, w, 1.3))
    assert np.max(np.abs(got - ref)) < 1e-14


def test_bloch_coordinates_are_the_pauli_expectations():
    cells = m.initial_cells(7)
    for coordinates, cell in zip(m.cells_to_bloch(cells), cells):
        rho = np.outer(cell, np.conj(cell))
        assert np.allclose(coordinates, [np.trace(rho @ S).real for S in (SX, SY, SZ)])
        assert abs(np.linalg.norm(coordinates) - 1.0) < 1e-12


def test_dense_reference_is_the_same_physics():
    """The structure-aware dense update must equal a real matrix exponential.

    Otherwise the speedup would be measured against a different calculation.
    """
    q = 6
    cells = m.initial_cells(q)
    omegas = 0.7 + 0.3 * np.arange(q) / (q - 1)
    psi = m.product_state(cells)
    dense = m.apply_local_rotations_dense(psi, omegas, 1.3)
    exact = m.expm(-1j * m.dense_hamiltonian(q, omegas) * 1.3) @ psi
    assert np.max(np.abs(dense - exact)) < 1e-12
    factorized = m.product_state(m.propagate_factorized_amplitudes(cells, omegas, 1.3))
    assert np.max(np.abs(factorized - dense)) < 1e-12


def test_the_trajectory_never_leaves_the_product_manifold():
    """This is what the exponential separation is really over.

    The Hamiltonian has no inter-cell coupling, so the dense carrier stays a
    product state and the 2**q amplitude vector is never needed.
    """
    for q in (4, 6, 8):
        cells = m.initial_cells(q)
        omegas = 0.7 + 0.3 * np.arange(q) / max(q - 1, 1)
        evolved = m.apply_local_rotations_dense(m.product_state(cells), omegas, 1.3)
        singular = np.linalg.svd(evolved.reshape(2, -1), compute_uv=False)
        assert int(np.sum(singular > 1e-12)) == 1


def test_scaling_artifact():
    rows = json.loads(ARTIFACT.read_text())["rows"]
    last = rows[-1]
    assert last["q"] == 18
    assert last["carrier_complex_amplitudes"] == 2**18
    assert last["observable_real_variables"] == 3 * 18
    assert max(row["max_bloch_error"] for row in rows) < 1e-14
    # The separation is structural; the exact factor is host dependent.
    assert last["observable_vs_dense_carrier_speedup"] > 100
    assert last["observable_vs_factorized_amplitude_speedup"] > 1


def test_artifact_matches_the_generator_schema():
    """The shipped artifact must be one this script could have written."""
    data = json.loads(ARTIFACT.read_text())
    required = {
        "q",
        "carrier_complex_amplitudes",
        "observable_real_variables",
        "observable_median_seconds",
        "factorized_amplitude_median_seconds",
        "dense_carrier_median_seconds",
        "observable_vs_dense_carrier_speedup",
        "observable_vs_factorized_amplitude_speedup",
        "max_bloch_error",
    }
    for r in data["rows"]:
        assert required <= set(r), required - set(r)
        assert r["observable_vs_dense_carrier_speedup"] == (
            r["dense_carrier_median_seconds"] / r["observable_median_seconds"]
        )
        assert r["observable_vs_factorized_amplitude_speedup"] == (
            r["factorized_amplitude_median_seconds"] / r["observable_median_seconds"]
        )
    assert "separation_scope" in data["claims"]


def test_dense_carrier_cost_actually_grows_exponentially():
    rows = {r["q"]: r for r in json.loads(ARTIFACT.read_text())["rows"]}
    assert rows[18]["dense_carrier_median_seconds"] > 50 * rows[10]["dense_carrier_median_seconds"]
    assert rows[18]["observable_median_seconds"] < 5 * rows[4]["observable_median_seconds"]
