from __future__ import annotations

import json
from pathlib import Path

import pytest

from gpc_census.generation.bdr import (
    constraint_from_tau,
    normalize_tau_system,
    parse_python_tau_export,
    primitive_tau,
    tau_from_constraint,
)
from gpc_census.generation.model import IntegralConstraint


ROOT = Path(__file__).parents[1]


def test_rank7_tau_conversion_and_structural_split():
    raw = (
        (-2, 1, 1, 1, 1, -2, -2),
        (1, -2, 1, 1, -2, 1, -2),
        (1, 1, -2, -2, 1, 1, -2),
        (1, 1, -2, 1, -2, -2, 1),
        (0, 0, 0, 0, 0, 0, -1),
    )
    system = normalize_tau_system(raw, 3)
    assert system.d == 7
    assert len(system.inequalities) == 4
    assert system.structural_inequalities == (
        IntegralConstraint((1, 1, 1, 1, 1, 1, 0), 3),
    )

    table = json.loads(
        (ROOT / "src" / "gpc_census" / "data" / "constraints.json").read_text()
    )
    expected = {
        primitive_tau(
            tau_from_constraint(IntegralConstraint(tuple(row["coeffs"]), row["rhs"]), 3)
        )
        for row in table["3,7"]["inequalities"]
    }
    actual = {tau_from_constraint(row, 3) for row in system.inequalities}
    assert actual == expected


def test_all_full_dimensional_lookup_rows_round_trip_homogeneously():
    table = json.loads(
        (ROOT / "src" / "gpc_census" / "data" / "constraints.json").read_text()
    )
    for d in range(7, 11):
        for row in table[f"3,{d}"]["inequalities"]:
            original = IntegralConstraint(tuple(row["coeffs"]), row["rhs"])
            tau = tau_from_constraint(original, 3)
            assert tau_from_constraint(constraint_from_tau(tau, 3), 3) == tau


def test_parser_reads_literal_without_executing_other_statements():
    text = """
raise RuntimeError('must not run')
brut_inequations = [(1, -2, 1), (0, 0, -1)]
"""
    assert parse_python_tau_export(text) == ((1, -2, 1), (0, 0, -1))


@pytest.mark.parametrize(
    "text",
    [
        "x = []",
        "brut_inequations = []\nbrut_inequations = []",
        "brut_inequations = [(1, 2), (1,)]",
        "brut_inequations = [(1, 2.0)]",
    ],
)
def test_parser_rejects_malformed_exports(text: str):
    with pytest.raises(ValueError):
        parse_python_tau_export(text)
