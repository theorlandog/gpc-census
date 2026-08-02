import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.verify_cyclic_window_observability import (  # noqa: E402
    connected,
    categorical_edges,
    cyclic_windows,
    incidence,
    rank_q,
    verify_factor_range,
    verify_product,
)


def test_cyclic_rank_formula_small_range():
    verify_factor_range(12)


def test_sample_factor_ranks():
    assert rank_q(incidence(cyclic_windows(5, 2), 1)) == 5
    assert rank_q(incidence(cyclic_windows(7, 3), 1)) == 7
    assert rank_q(incidence(cyclic_windows(6, 3), 1)) == 4
    assert rank_q(incidence(cyclic_windows(6, 3), 2)) == 6


def test_product_incidence_and_connectivity():
    verify_product(5, 2, 5, 2)
    verify_product(5, 2, 7, 3)
    verify_product(6, 2, 5, 2)


def test_even_even_categorical_product_is_disconnected():
    assert not connected(24, categorical_edges(4, 6))


def test_carrier_collision_free_graphs_are_connected():
    """POWER TEST. Proposition 3 and Theorem 4 on the carriers themselves.

    The rank checks above and the categorical-product check are both abstract:
    neither builds a carrier's actual collision-free graph. This does, for the
    coprime cyclic windows the theorem covers.
    """
    from scripts.verify_cyclic_window_observability import (
        collision_free_edges,
        verify_carrier_graph,
    )

    for d, n in ((5, 2), (7, 3), (9, 4), (8, 3)):
        carrier = cyclic_windows(d, n)
        assert len(collision_free_edges(carrier)) >= d
        verify_carrier_graph(carrier, True)


def test_glued_carrier_graph_is_connected_not_just_its_categorical_subgraph():
    """Theorem 6 concludes about the product carrier's own graph."""
    from scripts.verify_cyclic_window_observability import (
        product_carrier,
        verify_carrier_graph,
    )

    for d1, n1, d2, n2 in ((3, 2, 5, 2), (5, 2, 5, 2), (4, 3, 5, 2)):
        product = product_carrier(
            cyclic_windows(d1, n1, 0), cyclic_windows(d2, n2, d1)
        )
        verify_carrier_graph(product, True)
