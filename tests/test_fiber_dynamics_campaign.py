import importlib.util
import sys
from pathlib import Path

P = Path(__file__).parents[1] / "scripts" / "fiber_dynamics_campaign.py"
spec = importlib.util.spec_from_file_location("campaign", P)
mod = importlib.util.module_from_spec(spec)
# Registered before exec_module because @dataclass resolves its module through
# sys.modules; without this the import raises AttributeError.
sys.modules[spec.name] = mod
assert spec.loader is not None
spec.loader.exec_module(mod)


def test_cross_density_interaction_is_entangling():
    result = mod.entanglement_witness(5, 2, 0, 5, 2, 0)
    assert result["energy_minor_cross_difference"] == "1"
    assert not result["energy_is_additively_separable"]
    assert result["generic_schmidt_rank_lower_bound"] == 2


def test_amplitude_free_diagonal_propagation():
    coords = mod.uniform_coordinates(4)
    moved = coords.propagate_diagonal([0.0, 1.0, 2.0, 3.0], 0.2)
    assert moved.weights == coords.weights
    assert abs(abs(moved.root_coherences[2]) - 0.25) < 1e-14
    assert abs(moved.reconstruct_x(1, 3) - moved.root_coherences[0].conjugate() * moved.root_coherences[2] / moved.weights[0]) < 1e-14


def test_condition_bound():
    assert mod.local_path_condition_bound(1, 0.1, 1e-3) == 1e-3
    assert mod.local_path_condition_bound(4, 0.1, 1e-3) == 7e-3


def test_linear_identifiability_detects_collided_resolution():
    # y1=x0+x1, y2=x0-x1 determines both despite neither row being singleton.
    known = mod.linearly_identifiable_columns([[1, 1], [1, -1]])
    assert known == {0, 1}


def test_linear_identifiability_leaves_underdetermined_coordinates():
    assert mod.linearly_identifiable_columns([[1, 1]]) == set()
