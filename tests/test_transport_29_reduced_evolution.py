"""Exact reduced evolution on the 29-determinant transport carrier."""
import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

P = ROOT / "scripts" / "transport_29_reduced_evolution.py"
spec = importlib.util.spec_from_file_location("transport_29_reduced_evolution", P)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
assert spec.loader is not None
spec.loader.exec_module(mod)


def test_carrier_and_hamiltonian():
    support, hop, adjoint = mod.build_support()
    hamiltonian = mod.hamiltonian(support, [hop, adjoint])
    assert len(support) == 29
    assert (hamiltonian - hamiltonian.conj().T).max() == 0


def test_observability_tree():
    support, _, _ = mod.build_support()
    _, _, rows = mod.pair_incidence(support)
    tree = mod.tree_witnesses(support)
    assert len(rows) == 29
    assert len(tree) == 28
