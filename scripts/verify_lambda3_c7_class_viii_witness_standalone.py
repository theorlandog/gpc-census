#!/usr/bin/env python3
"""Dependency-light independent replay of the class-VIII witness.

No project module is imported.  The verifier reconstructs exterior signs,
the exact 1-RDM, the characteristic polynomial, and the full 49-column orbit
tangent matrix directly from the six recorded amplitudes.
"""
from __future__ import annotations

import argparse
import itertools
import json
import pathlib

import sympy as sp


def triples(d: int):
    return list(itertools.combinations(range(d), 3))


def derivation(a: int, b: int, state):
    out = {}
    for det, coefficient in state.items():
        if b not in det:
            continue
        raw = list(det)
        raw[det.index(b)] = a
        if len(set(raw)) != 3:
            continue
        sign = 1
        for i in range(3):
            for j in range(i + 1, 3):
                if raw[i] > raw[j]:
                    raw[i], raw[j] = raw[j], raw[i]
                    sign *= -1
        key = tuple(raw)
        out[key] = out.get(key, sp.Integer(0)) + sign * coefficient
    return {key: value for key, value in out.items() if value != 0}


def one_rdm(state, d: int):
    annihilated = []
    for mode in range(d):
        column = {}
        for det, coefficient in state.items():
            if mode not in det:
                continue
            remainder = tuple(x for x in det if x != mode)
            sign = (-1) ** det.index(mode)
            column[remainder] = column.get(remainder, sp.Integer(0)) + sign * coefficient
        annihilated.append(column)
    rho = sp.zeros(d)
    for i in range(d):
        for j in range(d):
            rho[i, j] = sp.simplify(sum(
                sp.conjugate(value) * annihilated[j].get(key, sp.Integer(0))
                for key, value in annihilated[i].items()
            ))
    return rho


def orbit_dimension(state, d: int):
    basis = triples(d)
    index = {det: i for i, det in enumerate(basis)}
    columns = []
    for a in range(d):
        for b in range(d):
            column = [sp.Integer(0)] * len(basis)
            if a == b:
                for det, coefficient in state.items():
                    if a in det:
                        column[index[det]] += coefficient
            else:
                for det, coefficient in derivation(a, b, state).items():
                    column[index[det]] += coefficient
            columns.append(sp.Matrix(column))
    return int(sp.Matrix.hstack(*columns).rank())


def verify(path: pathlib.Path) -> int:
    payload = json.loads(path.read_text())
    certificate = payload["certificate"]
    support = [tuple(int(ch) - 1 for ch in text) for text in certificate["support"]]
    weights = certificate["integer_squared_weights"]
    signs = certificate["signs"]
    denominator = certificate["denominator"]
    state = {
        det: sign * sp.sqrt(sp.Rational(weight, denominator))
        for det, weight, sign in zip(support, weights, signs)
    }

    norm = sp.simplify(sum(sp.conjugate(value) * value for value in state.values()))
    assert norm == 1
    rho = one_rdm(state, 7)
    expected = sp.diag(
        sp.Rational(5, 9), sp.Rational(5, 9), sp.Rational(5, 9),
        sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(1, 3),
        sp.Rational(1, 3),
    )
    assert rho == expected
    characteristic = sp.factor(rho.charpoly().as_expr())
    symbol = rho.charpoly().gen
    expected_characteristic = sp.expand(
        (symbol - sp.Rational(5, 9)) ** 3
        * (symbol - sp.Rational(1, 3)) ** 4
    )
    assert sp.expand(characteristic - expected_characteristic) == 0
    assert rho.rank() == 7
    assert orbit_dimension(state, 7) == 34

    print("PASS: exact norm, 1-RDM, characteristic polynomial, rank, and orbit dimension")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "artifact", nargs="?",
        default="results/data/lambda3_c7_class_viii_witness.json",
    )
    args = parser.parse_args()
    return verify(pathlib.Path(args.artifact))


if __name__ == "__main__":
    raise SystemExit(main())
