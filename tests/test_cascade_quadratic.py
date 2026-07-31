"""Exact quadratic certificates for the last two cascade endpoints."""
from __future__ import annotations

import json
import pathlib

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data"
if not DATA.exists():
    DATA = ROOT / "results/data"


def _rho(amplitudes, determinants, d):
    positions = {det: i for i, det in enumerate(determinants)}
    rho = sp.zeros(d, d)
    for column, det in enumerate(determinants):
        for q_index, q in enumerate(det):
            rest = tuple(o for o in det if o != q)
            for p in range(d):
                if p in rest:
                    continue
                target = tuple(sorted(rest + (p,)))
                row = positions.get(target)
                if row is None:
                    continue
                rho[p, q] += (
                    amplitudes[row]
                    * amplitudes[column]
                    * (-1) ** q_index
                    * (-1) ** target.index(p)
                )
    return sp.simplify(rho)


def _records():
    report = json.loads((DATA / "cascade_exact.json").read_text())
    return {
        record["index"]: record
        for record in report["results"]
        if record["system"] == "(5,10)" and record["index"] in (144, 256)
    }


def test_last_two_endpoints_are_exact_over_one_quadratic_field():
    records = _records()
    assert set(records) == {144, 256}
    for record in records.values():
        assert record["status"] == "SPECTRUM-ONLY-NOT-DIAGONAL"
        assert record["squared_weight_field"] == "Q(sqrt(1065))"
        assert record["squared_weight_max_degree"] == 2
        weights = [sp.sympify(value) for value in record["squared_weights"]]
        assert sp.simplify(sum(weights) - 1) == 0
        assert all(weight.is_positive for weight in weights)


def test_quadratic_endpoints_pass_the_exact_spectrum_gate():
    spectra = {
        144: [sp.Rational(34, 37)] * 3
        + [sp.Rational(17, 37)] * 3
        + [sp.Rational(8, 37)] * 4,
        256: [sp.Rational(29, 37)] * 4
        + [sp.Rational(20, 37)] * 3
        + [sp.Rational(3, 37)] * 3,
    }
    x = sp.Symbol("x")
    for index, record in _records().items():
        weights = [sp.sympify(value) for value in record["squared_weights"]]
        amplitudes = [
            sign * sp.sqrt(weight)
            for sign, weight in zip(record["signs"], weights)
        ]
        determinants = [tuple(det) for det in record["support_dets"]]
        rho = _rho(amplitudes, determinants, 10)
        assert sp.simplify(rho.charpoly(x).as_expr() - sp.prod(
            x - value for value in spectra[index]
        )) == 0
        assert any(
            sp.simplify(rho[i, j]) != 0
            for i in range(10)
            for j in range(10)
            if i != j
        )
