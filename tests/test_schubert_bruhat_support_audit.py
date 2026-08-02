import importlib.util
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "schubert_bruhat_support_audit.json"
CENSUS = DATA / "rdm_observability_census.json"

MODULE_PATH = ROOT / "scripts" / "schubert_bruhat_support_audit.py"
SPEC = importlib.util.spec_from_file_location("schubert_bruhat", MODULE_PATH)
sb = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sb)


def rows():
    return json.loads(ARTIFACT.read_text())["cyclic_windows"]


def test_kappa_is_gcd_minus_one_on_every_row():
    for r in rows():
        assert r["non_toric_phase_complexity"] == math.gcd(r["d"], r["N"]) - 1


def test_singleton_rank_matches_the_cyclic_formula():
    for d, n in ((6, 3), (8, 4), (9, 3), (10, 5), (7, 3)):
        window = sb.cyclic_windows(d, n)
        assert int(sb.singleton_matrix(window, d).rank()) == d - math.gcd(d, n) + 1


def test_support_polytope_dimension_is_rank_minus_one():
    """kappa = m - rank only equals (m-1) - dim P_S because dim P_S = rank - 1.

    Recomputed from differences of occupation vectors rather than assumed.
    """
    for d, n in ((6, 3), (8, 4), (9, 3)):
        window = sb.cyclic_windows(d, n)
        rank = int(sb.singleton_matrix(window, d).rank())
        base = [1 if p in window[0] else 0 for p in range(d)]
        deltas = [[(1 if p in D else 0) - base[p] for p in range(d)] for D in window[1:]]
        assert int(sb.sp.Matrix(deltas).rank()) == rank - 1


def test_only_nontrivial_cyclic_matroid_is_the_four_cycle():
    nontrivial = [
        (r["d"], r["N"]) for r in rows() if r["matroid_basis_set"] and 1 < r["N"] < r["d"] - 1
    ]
    assert nontrivial == [(4, 2)]
    assert sb.is_matroid_basis_set(sb.cyclic_windows(4, 2))[0]
    assert not sb.is_matroid_basis_set(sb.cyclic_windows(5, 2))[0]


def test_no_nontrivial_cyclic_gale_order_ideal():
    nontrivial = [
        (r["d"], r["N"]) for r in rows() if r["gale_order_ideal"] and 1 < r["N"] < r["d"] - 1
    ]
    assert nontrivial == []


def test_matroid_failure_witnesses_really_fail_exchange():
    """A recorded witness is worthless if it does not reproduce."""
    for r in rows():
        witness = r["matroid_failure_witness"]
        if witness is None or "A" not in witness:
            continue
        bases = {frozenset(D) for D in sb.cyclic_windows(r["d"], r["N"])}
        A, B = frozenset(witness["A"]), frozenset(witness["B"])
        x = witness["exchange_element"]
        assert A in bases and B in bases and x in A - B
        assert not any(frozenset((A - {x}) | {y}) in bases for y in B - A)


def test_gale_failure_witnesses_really_fail_the_ideal_property():
    for r in rows():
        witness = r["gale_failure_witness"]
        if witness is None:
            continue
        support = {tuple(D) for D in sb.cyclic_windows(r["d"], r["N"])}
        upper = tuple(witness["included"])
        lower = tuple(witness["missing_lower"])
        assert upper in support and lower not in support
        assert sb.gale_leq(lower, upper)


def test_glued_carriers_have_large_kappa_even_from_simplex_factors():
    """Section 9's warning, recomputed.

    Both factors are coprime cyclic windows with kappa = 0, yet the full glued
    carrier is far from a support simplex. The orbital-toric orbit is a proper
    subvariety of the glued carrier, so a product benchmark must say which one
    it propagates.
    """
    left = sb.cyclic_windows(5, 2)
    right = [tuple(x + 5 for x in D) for D in sb.cyclic_windows(5, 2)]
    assert int(sb.singleton_matrix(left, 5).rank()) == len(left)
    product = [tuple(sorted(set(a) | set(b))) for a in left for b in right]
    rank = int(sb.singleton_matrix(product, 10).rank())
    assert len(product) == 25
    assert len(product) - rank == 16


def test_the_28_census_exceptions_all_have_zero_kappa():
    """The doc's hypothesis, settled against the census.

    Section 7 suggests the order-two obstruction may sit in the internal phase
    sector. It does not: every unresolved record is a support simplex, and
    every positive-kappa record reconstructs at order two.
    """
    records = json.loads(CENSUS.read_text())["records"]

    def kappa(r):
        return r["support_size"] - r["singleton_incidence"]["rank"]

    unresolved = [r for r in records if r["first_connected_unique_channel_order"] != 2]
    positive = [r for r in records if kappa(r) > 0]
    assert len(unresolved) == 28
    assert len(positive) == 63
    assert all(kappa(r) == 0 for r in unresolved)
    assert all(r["first_connected_unique_channel_order"] == 2 for r in positive)
    # kappa > 0 first becomes possible at support size 6, and 26 of the 28
    # failures sit below that, so the association is confounded by size.
    assert min(r["support_size"] for r in positive) == 6
    assert sum(1 for r in unresolved if r["support_size"] < 6) == 26


def test_cyclic_windows_are_a_cyclic_shift_family():
    for d, n in ((7, 3), (9, 4)):
        window = sb.cyclic_windows(d, n)
        assert len(window) == d
        assert len(set(window)) == d
        for D in window:
            shifted = tuple(sorted((p + 1) % d for p in D))
            assert shifted in set(window)
        assert all(len(D) == n for D in window)


def test_gale_order_is_the_componentwise_order():
    assert sb.gale_leq((0, 1), (0, 2))
    assert not sb.gale_leq((0, 2), (0, 1))
    assert sb.gale_leq((0, 1, 2), (0, 1, 2))
    universe = list(itertools.combinations(range(5), 2))
    assert all(sb.gale_leq(subset, subset) for subset in universe)
