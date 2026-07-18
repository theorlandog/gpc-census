"""Constraint systems for fermionic moment polytopes, ranks 6 through 10.

Version 1 ships these as a lookup table (see data/constraints.json) with
particle-hole duality transport for complement systems. Constraint
generation for d >= 11 is future work.
"""
from __future__ import annotations

import json
from importlib import resources
from typing import TypedDict


class Inequality(TypedDict):
    coeffs: list[int]
    rhs: int


class System(TypedDict):
    inequalities: list[Inequality]
    equalities: list[Inequality]
    source: str


def _table() -> dict[str, System]:
    with resources.files("gpc_census.data").joinpath("constraints.json").open() as f:
        return json.load(f)


def known_systems() -> list[tuple[int, int]]:
    """Directly tabulated (n, d) systems."""
    return sorted((int(k.split(",")[0]), int(k.split(",")[1])) for k in _table())


def constraints(n: int, d: int) -> System:
    """Return the constraint system for n fermions in d orbitals.

    Directly tabulated systems are returned as recorded. Complement systems
    (d - n, d) are derived by particle-hole duality: lambda_i -> 1 -
    lambda_{d+1-i} maps each valid inequality to a valid inequality.
    """
    key = f"{n},{d}"
    tab = _table()
    if key in tab:
        return tab[key]
    dual_key = f"{d - n},{d}"
    if dual_key in tab:
        sysd = tab[dual_key]

        def flip(iq: Inequality) -> Inequality:
            a = iq["coeffs"]
            return {"coeffs": [-a[d - 1 - i] for i in range(d)], "rhs": iq["rhs"] - sum(a)}

        return {
            "inequalities": [flip(iq) for iq in sysd["inequalities"]],
            "equalities": [flip(iq) for iq in sysd["equalities"]],
            "source": sysd["source"] + " (particle-hole dual)",
        }
    raise KeyError(
        f"No tabulated constraints for ({n},{d}); rank {d} exceeds the known frontier (d <= 10)"
    )
