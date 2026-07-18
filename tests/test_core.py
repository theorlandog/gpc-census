from math import comb

import pytest

from gpc_census import slater_vertices


def test_vertex_count():
    assert len(list(slater_vertices(6, 3))) == comb(6, 3)


def test_vertices_are_binary_with_correct_sum():
    for vec in slater_vertices(5, 2):
        assert len(vec) == 5
        assert set(vec) <= {0, 1}
        assert sum(vec) == 2


def test_zero_orbitals():
    assert list(slater_vertices(0, 0)) == [()]


def test_too_many_fermions():
    with pytest.raises(ValueError):
        list(slater_vertices(2, 3))


def test_negative_arguments():
    with pytest.raises(ValueError):
        list(slater_vertices(-1, 0))
