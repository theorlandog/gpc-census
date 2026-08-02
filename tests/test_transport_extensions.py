"""Exact transport-extension certificates for observable carriers."""
import importlib.util
import pathlib
import sys

P = pathlib.Path(__file__).parents[1] / "scripts" / "transport_extensions.py"
spec = importlib.util.spec_from_file_location("transport_extensions", P)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
assert spec.loader is not None
spec.loader.exec_module(mod)


def test_scheduled_tangent():
    result = mod.state_specific_tangent_certificate()
    assert result["exact_checks"]["leakage_zero"]
    assert result["exact_checks"]["tangent_nonprojective"]
    assert result["state_tangent_kernel_dimension"] > 0


def test_enlarged_closure_observable():
    result = mod.closure_certificate()
    assert result["enlarged_carrier_dimension"] == 29
    assert result["pair_incidence_full_column_rank"]
    assert result["collision_free_gamma2_connected"]
    assert result["invariant_under_selected_hermitian_hop"]


def test_active_families():
    for result in (
        mod.fixed_core_family(2, 5, 1),
        mod.fixed_core_family(2, 5, 2),
        mod.fixed_core_family(3, 6, 2),
    ):
        assert result["pair_incidence_full_column_rank"]
        assert result["collision_free_gamma2_connected"]
        assert result["guaranteed_active_offdiagonal_normalizer_dimension"] > 0
