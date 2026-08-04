"""Exact bridge to the Bulois-Denis-Ressayre moment-cone generator.

The external generator represents an inequality by an integral vector
``tau`` and the homogeneous relation ``tau . lambda <= 0``.  The census
uses affine relations ``a . lambda <= b`` on the trace slice
``sum(lambda) = particle_number``.  This module converts between the two
forms without floating point arithmetic and parses the external Python
export without executing it.

This is an interoperability layer, not a validity oracle.  A row parsed or
normalized here still needs an exact Ressayre or Schubert certificate.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from functools import reduce
from math import gcd

from .model import IntegralConstraint


Tau = tuple[int, ...]


def _content(values: tuple[int, ...]) -> int:
    """Return the positive gcd of ``values``, with zero content mapped to one."""

    return reduce(gcd, (abs(value) for value in values), 0) or 1


def primitive_tau(tau: tuple[int, ...]) -> Tau:
    """Remove a positive common factor while preserving inequality orientation."""

    if not tau or not any(tau):
        raise ValueError("tau must be a nonzero vector")
    divisor = _content(tau)
    return tuple(value // divisor for value in tau)


def tau_from_constraint(
    constraint: IntegralConstraint,
    particle_number: int,
) -> Tau:
    """Return the primitive homogeneous vector for an affine constraint.

    On the trace slice, ``a.lambda <= b`` is equivalent to
    ``(N*a - b*1).lambda <= 0``.  Positive rescaling does not change the
    oriented inequality.
    """

    if particle_number <= 0:
        raise ValueError("particle_number must be positive")
    tau = tuple(
        particle_number * value - constraint.rhs for value in constraint.coeffs
    )
    return primitive_tau(tau)


def constraint_from_tau(tau: Tau, particle_number: int) -> IntegralConstraint:
    """Choose a deterministic nonnegative affine representative of ``tau``.

    A positive multiplier is first chosen so that all entries have one
    residue modulo ``particle_number``.  Such a multiplier always exists,
    since multiplying by ``particle_number`` suffices.  The trace gauge is
    then fixed by requiring the smallest affine coefficient to be zero.
    """

    if particle_number <= 0:
        raise ValueError("particle_number must be positive")
    tau = primitive_tau(tau)
    scale = next(
        multiplier
        for multiplier in range(1, particle_number + 1)
        if len({(multiplier * value) % particle_number for value in tau}) == 1
    )
    scaled = tuple(scale * value for value in tau)
    rhs = (-scaled[0]) % particle_number
    coeffs = tuple((value + rhs) // particle_number for value in scaled)

    shift = -min(coeffs)
    coeffs = tuple(value + shift for value in coeffs)
    rhs += particle_number * shift

    divisor = _content(coeffs + (rhs,))
    return IntegralConstraint(
        tuple(value // divisor for value in coeffs),
        rhs // divisor,
    )


def is_pauli_lower_tau(tau: Tau) -> bool:
    """Return whether ``tau`` is the dominant-chamber row ``lambda_d >= 0``."""

    tau = primitive_tau(tau)
    return tau == (0,) * (len(tau) - 1) + (-1,)


def parse_python_tau_export(text: str) -> tuple[Tau, ...]:
    """Extract ``brut_inequations`` from an upstream Python export safely.

    Only a literal assignment is accepted.  Imports and all other statements
    in the export are ignored rather than executed.
    """

    tree = ast.parse(text)
    assignments = [
        node
        for node in tree.body
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        and (
            (
                isinstance(node, ast.Assign)
                and any(
                    isinstance(target, ast.Name) and target.id == "brut_inequations"
                    for target in node.targets
                )
            )
            or (
                isinstance(node, ast.AnnAssign)
                and isinstance(node.target, ast.Name)
                and node.target.id == "brut_inequations"
            )
        )
    ]
    if len(assignments) != 1:
        raise ValueError("expected exactly one brut_inequations assignment")
    assignment = assignments[0]
    value_node = assignment.value
    if value_node is None:
        raise ValueError("brut_inequations has no value")
    value = ast.literal_eval(value_node)
    if not isinstance(value, (list, tuple)):
        raise ValueError("brut_inequations must be a list or tuple")

    rows: list[Tau] = []
    width: int | None = None
    for raw_row in value:
        if not isinstance(raw_row, (list, tuple)) or not raw_row:
            raise ValueError("every tau row must be a nonempty list or tuple")
        if any(type(entry) is not int for entry in raw_row):
            raise ValueError("tau entries must be integers")
        row = tuple(raw_row)
        if width is None:
            width = len(row)
        elif len(row) != width:
            raise ValueError("tau rows have inconsistent dimensions")
        rows.append(row)
    return tuple(rows)


@dataclass(frozen=True)
class NormalizedTauSystem:
    """Canonical affine GPC rows and separated structural Pauli rows."""

    d: int
    particle_number: int
    raw_taus: tuple[Tau, ...]
    inequalities: tuple[IntegralConstraint, ...]
    structural_inequalities: tuple[IntegralConstraint, ...]


def normalize_tau_system(
    raw_taus: tuple[Tau, ...],
    particle_number: int,
) -> NormalizedTauSystem:
    """Deduplicate oriented tau rows and convert them to the trace slice."""

    if not raw_taus:
        raise ValueError("raw_taus must not be empty")
    d = len(raw_taus[0])
    if any(len(row) != d for row in raw_taus):
        raise ValueError("tau rows have inconsistent dimensions")

    unique = tuple(sorted({primitive_tau(row) for row in raw_taus}))
    structural_taus = tuple(row for row in unique if is_pauli_lower_tau(row))
    gpc_taus = tuple(row for row in unique if not is_pauli_lower_tau(row))
    return NormalizedTauSystem(
        d=d,
        particle_number=particle_number,
        raw_taus=unique,
        inequalities=tuple(sorted(constraint_from_tau(row, particle_number) for row in gpc_taus)),
        structural_inequalities=tuple(
            sorted(constraint_from_tau(row, particle_number) for row in structural_taus)
        ),
    )
