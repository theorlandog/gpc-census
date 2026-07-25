"""Levi-orbit tests.

REGRESSION GUARD: the original `levi_theorem_consistent` reduced to
`slater == slater` and could never fail.  test_theorem_check_is_falsifiable
pins that it can now return False on malformed input.
"""
from gpc_census.levi import (
    block_sector_dimensions,
    full_levi_line_possible,
    levi_dimension,
    spectrum_multiplicities,
    static_levi_record,
)


def test_static_branching_empty_or_full_only():
    sectors = block_sector_dimensions(3, [2, 3])
    lines = [r for r in sectors if r["projective_line_possible"]]
    assert lines == [{"block_occupations": [0, 3], "dimension": 1,
                      "projective_line_possible": True}]
    assert full_levi_line_possible(3, [2, 3])
    assert not full_levi_line_possible(2, [4])


def test_static_record_slater_and_correlated():
    slater = {"system": "(2,4)", "index": 0, "integer_form": [1, 1, 0, 0],
              "denominator": 1}
    correlated = {"system": "(2,4)", "index": 1, "integer_form": [1, 1, 1, 1],
                  "denominator": 2}
    a = static_levi_record(slater)
    b = static_levi_record(correlated)
    assert a["multiplicities"] == [2, 2]
    assert a["levi_dim"] == 8
    assert a["representation_contains_full_levi_line"]
    assert a["target_full_levi_line_possible"]
    assert a["slater_spectrum"]
    assert a["levi_theorem_consistent"]

    assert b["multiplicities"] == [4]
    assert b["levi_dim"] == 16
    # FIXED: Lambda^2 C^4 has no one-dimensional sector, so this is False.
    # The previous suite asserted True here while another test asserted the
    # negation of the same quantity; both could not pass.
    assert not b["representation_contains_full_levi_line"]
    assert not b["target_full_levi_line_possible"]
    assert not b["slater_spectrum"]
    assert b["levi_theorem_consistent"]
    assert b["hard_obstruction"]


def test_theorem_check_is_falsifiable():
    """Malformed record: sum(integer_form) != n * denominator."""
    bad = {"system": "(1,4)", "index": 0, "integer_form": [3, 3, 0, 0],
           "denominator": 3}
    assert static_levi_record(bad)["levi_theorem_consistent"] is False


def test_levi_dimension():
    assert levi_dimension([2, 3]) == 13
    assert spectrum_multiplicities([2, 2, 1, 0, 0])[0] == [2, 1, 2]
