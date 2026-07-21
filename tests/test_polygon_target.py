"""Tests for the polygon-target solver (scripts/polygon_target.py).

The solver proposes exact phases and accepts a state only if
gpc_census.exactify.verify_exact certifies it, so these tests measure
COMPLETENESS (does it find a solution when one exists) and the two exact prunes;
soundness is guaranteed by the gate. See docs/RESEARCH.md, "v96 campaign".
"""
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import polygon_target as pt  # noqa: E402


def test_recertify_3_8_magnitude_path():
    # every certified (3,8) interference state must be re-solved from its
    # weights + spectrum alone (phases discarded): the magnitude/block path
    # plus the verify_exact gate, on the real census corpus
    ok, n = pt.recertify(3, 8)
    assert ok == n == 11


def test_cancellation_positive_four_cycle():
    # (|01> + |02> + |13> - |23>)/2 has spectrum (1/2,1/2,1/2,1/2); its two
    # 2-term one-hop classes must CANCEL (off-block target 0). The solver must
    # rediscover the sign that makes rho diagonal.
    dets = [(0, 1), (0, 2), (1, 3), (2, 3)]
    ks = [1, 1, 1, 1]
    den = 4
    spec = [Fraction(1, 2)] * 4
    targets = {pq: sp.Integer(0) for pq in pt.one_hop_classes(4, dets)}
    rec = pt.solve(2, 4, spec, dets, ks, den, targets=targets)
    assert rec is not None and rec["status"] == "EXACT"
    assert rec["weights"] == ks


def test_lone_pair_funnel_infeasible():
    # a single-term one-hop class forced to cancel is a lone nonzero vector:
    # exactly the P1 prune, must return None before any solve
    dets = [(0, 1), (0, 2)]
    ks = [1, 1]
    den = 2
    spec = [Fraction(1), Fraction(1, 2), Fraction(1, 2)]
    targets = {(1, 2): sp.Integer(0)}
    assert pt.solve(2, 3, spec, dets, ks, den, targets=targets) is None
