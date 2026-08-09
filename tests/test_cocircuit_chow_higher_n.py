"""Guards for the cocircuit route and the Chow decoder at arbitrary `N`.

The live checks run at `(2,6)`, `(4,6)` and `(4,7)`, which recompute in a few
seconds. `(5,8)`, the published higher-N gates and the half-filling wall are
read from the artifact, which is what `scripts/cocircuit_chow_higher_n.py`
exists to produce.
"""

from __future__ import annotations

import itertools
import json
import pathlib
import random

import pytest

from gpc_census.generation.admissible_flats import (
    enumerate_admissible_hyperplanes_by_flats,
)
from gpc_census.generation.chow_decoder import (
    complement_shadow,
    decode_shadow,
    decode_signed_pair,
    shadow_of,
    threshold_family,
)
from gpc_census.generation.cocircuit import (
    enumerate_cocircuit_hyperplanes,
    trace_survivors,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
DOCS = ROOT / "docs"
ARTIFACT = DATA / "cocircuit_chow_higher_n.json"


def _artifact() -> dict:
    return json.loads(ARTIFACT.read_text())


def _orbit_key(tau):
    return min(tuple(sorted(tau)), tuple(sorted(-value for value in tau)))


def _profile(subset, blocks):
    return tuple(len(set(subset) & set(block)) for block in blocks)


# ------------------------------------------------- live, small higher-N gates


@pytest.mark.parametrize("n", [2, 4])
def test_particle_hole_duals_have_one_population_at_rank_six(n):
    """`(2,6)` and `(4,6)` are dual, so the populations must coincide."""
    enumeration = enumerate_cocircuit_hyperplanes(n, 6)
    assert enumeration.hyperplane_count == 112
    rows = list(trace_survivors(enumeration))
    assert len(rows) == 28
    assert sum(1 for _tau, positive in rows if positive > 0) == 26
    assert len({_orbit_key(tau) for tau in enumeration.normals}) == 6


@pytest.mark.parametrize("n", [2, 4])
def test_higher_n_agrees_with_the_flat_lattice_pair_for_pair(n):
    enumeration = enumerate_cocircuit_hyperplanes(n, 6)
    reference = enumerate_admissible_hyperplanes_by_flats(n, 6)
    assert set(enumeration.affine_pairs()) == {
        (row.h, row.z) for row in reference.hyperplanes
    }


def test_rank_seven_dual_reproduces_the_n3_population():
    """`(4,7)` must reproduce `(3,7)`: 5,341 hyperplanes, 19 orbits, 549 rows."""
    enumeration = enumerate_cocircuit_hyperplanes(4, 7)
    assert enumeration.hyperplane_count == 5341
    assert len({_orbit_key(tau) for tau in enumeration.normals}) == 19
    rows = list(trace_survivors(enumeration))
    assert len(rows) == 549
    assert sum(1 for _tau, positive in rows if positive > 0) == 547


def test_decoder_inverts_every_survivor_at_four_six():
    """The full `(4,6)` nonstructural population, decoded and reconstructed."""
    enumeration = enumerate_cocircuit_hyperplanes(4, 6)
    rows = 0
    for tau, positive in trace_survivors(enumeration):
        if positive == 0:
            continue
        rows += 1
        plus = threshold_family(tau, 4, +1)
        minus = threshold_family(tau, 4, -1)
        decoded = decode_signed_pair(shadow_of(plus, 6), shadow_of(minus, 6), 4)
        assert decoded is not None
        assert decoded.positive.family == plus
        assert decoded.negative.family == minus
        assert decoded.tau == tuple(tau)
    assert rows == 26


# ------------------------------------------------ the particle-hole reduction


def test_complement_shadow_is_the_complement_family_shadow():
    """`c(F^c) = m - c(F)`, checked against the complement family itself."""
    rng = random.Random(311)
    n, d = 5, 8
    for _ in range(40):
        tau = tuple(rng.randrange(-5, 6) for _ in range(d))
        family = threshold_family(tau, n)
        if not family:
            continue
        shadow = shadow_of(family, d)
        complement, cardinality = complement_shadow(shadow, n)
        assert cardinality == len(family)
        universe = frozenset(range(d))
        dual_family = frozenset(
            tuple(sorted(universe - frozenset(subset))) for subset in family
        )
        assert complement == shadow_of(dual_family, d)


def test_particle_hole_route_changes_the_cost_and_not_the_answer():
    """Both routes decode the same family, with the SAME variable count."""
    rng = random.Random(312)
    checked = 0
    for n, d in ((5, 8), (6, 9), (7, 10)):
        for _ in range(12):
            tau = tuple(rng.randrange(-5, 6) for _ in range(d))
            family = threshold_family(tau, n)
            if not family:
                continue
            shadow = shadow_of(family, d)
            literal = decode_shadow(shadow, n, particle_hole=False)
            reduced = decode_shadow(shadow, n, particle_hole=True)
            assert literal is not None and reduced is not None
            assert literal.family == family
            assert reduced.family == family
            # k -> m - k is a bijection between the layers' profiles, so the
            # reduction is never a size reduction. This pins that.
            assert literal.profile_variables == reduced.profile_variables
            assert literal.particle_hole_layer == n
            assert reduced.particle_hole_layer == min(n, d - n)
            checked += 1
    assert checked > 20


def test_particle_hole_is_a_no_op_below_half_filling():
    """Nothing to complement when the literal layer is already the smaller one."""
    tau = (4, 2, 1, -1, -2, -4, 0, 3)
    family = threshold_family(tau, 3)
    decoded = decode_shadow(shadow_of(family, 8), 3)
    assert decoded is not None
    assert decoded.particle_hole_layer == 3


# ---------------------------------------------------------------- T4/T5 at N>3


def test_t4_and_t5_hold_at_higher_particle_number():
    """Block symmetry and the upper ideal, on random `N = 4` and `N = 5` rows."""
    rng = random.Random(313)
    checked = 0
    for n, d in ((4, 8), (5, 9)):
        for _ in range(20):
            tau = tuple(rng.randrange(-5, 6) for _ in range(d))
            family = threshold_family(tau, n)
            if not family:
                continue
            counts = shadow_of(family, d)
            for i in range(d):
                for j in range(i + 1, d):
                    if counts[i] != counts[j]:
                        continue
                    swapped = frozenset(
                        tuple(
                            sorted(
                                j if index == i else i if index == j else index
                                for index in subset
                            )
                        )
                        for subset in family
                    )
                    assert swapped == family
            decoded = decode_shadow(counts, n, particle_hole=False)
            assert decoded is not None and decoded.family == family
            checked += 1
    assert checked > 20


# ---------------------------------------------------------------- the artifact


def test_artifact_dual_populations_match_their_n3_duals():
    dual = _artifact()["dual_populations"]
    expected = {"(4,7)": (5341, 19, 549, 547), "(5,8)": (166420, 56, 6607, 6605)}
    for system, (planes, orbits, survivors, nonstructural) in expected.items():
        row = dual[system]
        assert row["admissible_hyperplanes"] == planes
        assert row["flat_lattice_hyperplanes"] == planes
        assert row["sd_unoriented_orbits"] == orbits
        assert row["trace_survivors"] == survivors
        assert row["nonstructural_trace_survivors"] == nonstructural
        assert row["agrees_with_flat_lattice_pairwise"] is True
        assert row["decode_failures"] == 0
        assert row["tau_reconstruction_failures"] == 0


def test_artifact_dual_populations_equal_the_n3_artifact():
    """The dual gate and the `N=3` gate must be reading one population."""
    dual = _artifact()["dual_populations"]
    base = json.loads((DATA / "cocircuit_chow_decoder.json").read_text())["cocircuit"]
    for system, row in dual.items():
        n, d = (int(part) for part in system.strip("()").split(","))
        partner = base[str(d)]
        assert partner["particle_number"] == d - n
        for field in (
            "admissible_hyperplanes",
            "sd_unoriented_orbits",
            "trace_survivors",
            "nonstructural_trace_survivors",
        ):
            assert row[field] == partner[field]


def test_artifact_published_higher_n_is_clean():
    published = _artifact()["published_higher_n"]
    assert set(published) == {"(4,8)", "(4,9)", "(4,10)", "(5,10)"}
    total = 0
    for row in published.values():
        assert row["particle_number"] >= 4
        for field in (
            "span_is_hyperplane",
            "trace_identity",
            "decoded",
            "tau_reconstructed_oriented",
        ):
            assert row[field] == row["representatives"]
        assert row["source_sha256"]
        total += row["representatives"]
    assert total == 61


def test_artifact_published_digests_match_the_live_artifacts():
    """A stale digest would mean the gate ran against different data."""
    import hashlib

    for row in _artifact()["published_higher_n"].values():
        path = ROOT / row["source_artifact"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["source_sha256"]


def test_artifact_particle_hole_never_changes_the_variable_count():
    """The correction, pinned: the reduction is exact and is not a shrink."""
    records = _artifact()["particle_hole"]
    assert records
    for row in records.values():
        assert row["answer_mismatches"] == 0
        assert row["rows_where_variable_count_differs"] == 0
        assert row["shadows_checked"] > 0
        assert row["reduced_layer"] < row["particle_number"]


def test_artifact_trace_survivor_stress_has_no_failures():
    stress = _artifact()["trace_survivor_stress"]
    assert stress
    for row in stress.values():
        assert row["failures"] == 0
        assert row["exact_trace_survivors_decoded"] > 0
        assert row["max_profile_variables"] < row["slice_size"]


def test_artifact_records_the_wall_without_a_partial_count():
    """A truncated reverse search must never publish a hyperplane count."""
    wall = _artifact()["half_filling_wall"]
    assert wall["system"] == "(4,8)"
    assert wall["completed"] is False
    assert wall["admissible_hyperplanes"] is None
    assert wall["node_budget"] >= 1_000_000


def test_artifact_keeps_proved_and_measured_apart():
    artifact = _artifact()
    assert set(artifact["theorems"]) == {"generality", "particle_hole_shadow"}
    measured = artifact["measured_not_proved"]
    assert "particle_hole_is_not_a_size_reduction" in measured
    assert "half_filling_wall" in measured
    assert "published_gate_scope" in measured
    assert artifact["machine"]["platform"]


# ------------------------------------------------------------------ the docs


@pytest.mark.parametrize(
    "phrase",
    [
        "it is not a size reduction",
        "no partial hyperplane count",
        "representative-level gates",
    ],
)
def test_the_document_states_the_higher_n_limits(phrase):
    text = (DOCS / "cocircuit_chow_decoder.md").read_text().lower()
    assert phrase.lower() in text


def test_prereg_scores_round_two():
    text = (DOCS / "prereg_cocircuit_chow_decoder.md").read_text()
    assert "ROUND 2" in text
    for identifier in ("P9", "P10", "P11", "P12", "P13"):
        assert f"| {identifier} |" in text


def test_exhaustive_t1_holds_at_four_six():
    """T1 evaluated over the whole power set of the `(4,6)` slice.

    The slice has 15 weights, so all 32,768 subfamilies fit in a loop. Every
    shadow a threshold family realizes must have exactly one preimage, which is
    T1 checked rather than argued, at `N = 4`.
    """
    n, d = 4, 6
    weights = list(itertools.combinations(range(d), n))
    counts: dict[tuple[int, ...], int] = {}
    for size in range(len(weights) + 1):
        for chosen in itertools.combinations(weights, size):
            shadow = [0] * d
            for subset in chosen:
                for index in subset:
                    shadow[index] += 1
            key = tuple(shadow)
            counts[key] = counts.get(key, 0) + 1

    enumeration = enumerate_cocircuit_hyperplanes(n, d)
    checked = 0
    for base in enumeration.normals:
        for sign in (1, -1):
            tau = tuple(sign * value for value in base)
            family = threshold_family(tau, n)
            if not family:
                continue
            assert counts[shadow_of(family, d)] == 1
            checked += 1
    assert checked > 100
