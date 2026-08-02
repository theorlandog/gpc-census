import importlib.util
import json
import sys
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "multichart_symplectic_propagator.json"

MODULE = ROOT / "scripts" / "multichart_symplectic_propagator.py"
spec = importlib.util.spec_from_file_location("multichart", MODULE)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


@lru_cache(maxsize=1)
def atlas():
    """Building the atlas costs seconds, so share one across the module."""
    return mod.Atlas.build()


def test_atlas_has_multiple_certified_charts():
    a = atlas()
    assert len(a.support) == 29
    assert len(a.edge_list) >= 69
    assert len(a.trees) >= 100


def test_every_chart_is_a_spanning_tree_of_certified_edges():
    a = atlas()
    certified = set(a.edge_list)
    for tree in a.trees:
        assert len(tree["edges"]) == len(a.support) - 1
        assert set(tree["edges"]) <= certified
        assert len(set(tree["parent"])) == len(a.support)


def test_observable_initialization_and_reconstruction():
    a = atlas()
    c0, _ = mod.build_initial_state(a.H)
    X = np.outer(np.conj(c0), c0)
    y = a.pack_from_X(X)
    rebuilt, _ = a.reconstruct_X(y)
    assert np.max(np.abs(rebuilt - X)) < 1e-12


def test_boundary_chart_crosses_exact_zero():
    a = atlas()
    c0, _ = mod.build_initial_state(a.H)
    t = np.pi / 4
    c = mod.expm(-1j * a.H * t) @ c0
    X = np.outer(np.conj(c), c)
    y = a.pack_from_X(X)
    rebuilt, tag = a.reconstruct_by_synchronization(y)
    assert tag == -2
    assert np.min(np.abs(c) ** 2) < 1e-28
    assert np.max(np.abs(rebuilt - X)) < 1e-12


def test_tree_chart_also_survives_the_zero_because_it_is_a_leaf():
    """The selector never needs the fallback here, and the reason is structural.

    A vanishing weight is harmless when its determinant is a leaf of the
    selected tree: leaves are never reconstruction denominators, so the induced
    forest stays connected. This is why the artifact records only tree charts.
    """
    a = atlas()
    c0, _ = mod.build_initial_state(a.H)
    for t in (np.pi / 4, 3 * np.pi / 4):
        c = mod.expm(-1j * a.H * t) @ c0
        weights = np.abs(c) ** 2
        vanishing = [int(i) for i in np.where(weights < 1e-28)[0]]
        assert vanishing
        index = a.choose_chart(weights)
        assert index is not None
        tree = a.trees[index]
        children = defaultdict(list)
        for vertex, parent in tree["parent"].items():
            if parent >= 0:
                children[parent].append(vertex)
        for vertex in vanishing:
            assert vertex != tree["root"]
            assert not children[vertex]


def test_artifact_records_the_structural_counts():
    art = json.loads(ARTIFACT.read_text())
    a = atlas()
    assert art["carrier_dimension"] == len(a.support) == 29
    assert art["collision_free_edge_bank_size"] == len(a.edge_list) == 69
    assert art["certified_spanning_trees"] == len(a.trees) == 108
    assert art["reduced_real_variables"] == 29 + 2 * 69 == 167
    assert art["physical_projective_dimension"] == 2 * 29 - 2 == 56
    assert art["ambient_fci_dimension"] == 210
    assert art["max_X_error"] < 1e-12
    assert art["max_Gamma2_error"] < 1e-12
    assert art["max_Gamma3_error"] < 1e-12
    assert art["minimum_direct_weight"] < 1e-30
    # Only precomputed tree charts were needed; -1 and -2 are the fallbacks.
    assert {sample["chart"] for sample in art["samples"]}.isdisjoint({-1, -2})
