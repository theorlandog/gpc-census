"""Guards for the symplectic slice atlas.

The forced identities are checked on every completed row; the exact
recomputation is re-run only for the two rank-6 states, which is enough to
catch a broken engine without making the suite slow.
"""
from __future__ import annotations

import json
import pathlib

import pytest

from gpc_census.symplectic import symplectic_slice

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" if (ROOT / "tests" / "data").exists() else ROOT / "results" / "data"
ARTIFACT = DATA / "symplectic_slices.json"

pytestmark = pytest.mark.skipif(not ARTIFACT.exists(), reason="slice atlas absent")


@pytest.fixture(scope="module")
def art():
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def states():
    return {(r["system"], r["index"]): r
            for r in (json.loads(x)
                      for x in (DATA / "states.jsonl").read_text().splitlines())}


def test_all_forced_gates_hold(art):
    for key in ("rank_equals_orbit_violations", "kernel_dim_violations",
                "levi_dim_violations"):
        assert not art["gates"][key], (key, art["gates"][key])


def test_gate_arithmetic_is_internally_consistent(art):
    """Re-derive each gate from the recorded numbers, not from its own flag."""
    for r in art["rows"]:
        s = r["slice"]
        assert s["rank_d_mu"] == s["orbit_dim"], r
        assert s["ker_d_mu_dim"] == s["ambient_tangent_dim"] - s["orbit_dim"], r
        assert s["levi_orbit_dim"] == s["commutant_dim"] - s["stab_dim"], r
        assert s["slice_dim"] == s["ker_d_mu_dim"] - s["levi_orbit_dim"], r
        assert s["slice_dim"] >= 0, r
        assert s["stabilizer_algebra_dim"] == s["stab_dim"], r


def test_unfinished_states_are_named_not_dropped(art):
    attempted = 5
    assert len(art["rows"]) + len(art["unfinished"]) == attempted
    keys = {f"{r['system']}#{r['index']}" for r in art["rows"]}
    assert not (keys & set(art["unfinished"]))


def test_borland_dennis_vertex_is_certified_rigid(art):
    """slice_dim == 0 certifies the reduced germ is a point."""
    row = next((r for r in art["rows"]
                if r["system"] == "(3,6)" and r["index"] == 3), None)
    assert row is not None
    assert row["slice"]["slice_dim"] == 0
    assert row["slice"]["slice_certifies_point"] is True
    assert row["slice"]["ker_d_mu_dim"] == 19


def test_slater_shows_a_nonzero_slice_over_a_point_reduced_space(art):
    """The whole point: a nonzero slice does not mean a positive-dimensional germ."""
    row = next((r for r in art["rows"]
                if r["system"] == "(3,6)" and r["index"] == 0), None)
    assert row is not None
    assert row["slice"]["slice_dim"] == 20
    assert row["slice"]["ker_d_mu_dim"] == 20
    assert row["slice"]["levi_orbit_dim"] == 0
    assert row["slice"]["slice_certifies_point"] is False


@pytest.mark.parametrize("index,expected", [(0, 20), (3, 0)])
def test_recompute_rank_six_slices_exactly(art, states, index, expected):
    got = symplectic_slice(states[("(3,6)", index)])
    assert got["slice_dim"] == expected
    assert got["gate_rank_equals_orbit"]
    assert got["gate_kernel_dim"]
    assert got["gate_levi_dim"]
