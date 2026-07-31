"""Normalization and gradient gates for the ansatz-free contraction attack."""
from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from contraction_attack import (  # noqa: E402
    attack,
    build,
    normalized_objective_gradient,
)


def _finite_difference(fun, x, k, h=1e-6):
    plus = x.copy()
    minus = x.copy()
    plus[k] += h
    minus[k] -= h
    return (fun(plus)[0] - fun(minus)[0]) / (2 * h)


def test_normalized_objective_is_scale_invariant():
    """The objective is on projective state space, not the unnormalized cone."""
    _, contraction = build(5, 2)
    g0 = np.diag([0.8, 0.7, 0.2, 0.2, 0.1]).ravel()
    rng = np.random.default_rng(0)
    for real in (True, False):
        dim = 10 if real else 20
        x = rng.normal(size=dim)
        f, g = normalized_objective_gradient(x, g0, contraction, real=real)
        scaled, _ = normalized_objective_gradient(7.0 * x, g0, contraction, real=real)
        assert abs(f - scaled) < 1e-13
        assert abs(float(g @ x)) < 1e-11


def test_normalized_gradient_matches_finite_differences():
    """Pin the quotient-rule term that removes radial norm drift."""
    _, contraction = build(5, 2)
    g0 = np.diag([0.8, 0.7, 0.2, 0.2, 0.1]).ravel()
    rng = np.random.default_rng(1)
    for real in (True, False):
        dim = 10 if real else 20
        x = rng.normal(size=dim)

        def fun(v):
            return normalized_objective_gradient(v, g0, contraction, real=real)

        _, gradient = fun(x)
        for k in (0, dim // 3, dim - 1):
            assert abs(_finite_difference(fun, x, k) - gradient[k]) < 2e-7


def test_normalized_attack_retains_an_attainable_control():
    """A Slater spectrum remains a zero-distance control after normalization."""
    residual, eigenvalues, support = attack(
        [1, 1, 0, 0],
        1,
        n=2,
        real=True,
        tries=2,
        seed=2,
    )
    assert residual < 1e-18
    assert np.max(np.abs(eigenvalues - np.array([1.0, 1.0, 0.0, 0.0]))) < 1e-9
    assert support >= 1
