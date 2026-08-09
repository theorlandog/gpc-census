"""Profile-native admissible-hyperplane enumeration.

Every number pinned here is one the repository already knew by an independent
route, so a passing run is a cross-algorithm regression rather than a
self-consistency check.  The sources are named beside each series.
"""

from __future__ import annotations

import itertools
import json
import pathlib

import pytest

from gpc_census.generation.admissible_flats import (
    enumerate_admissible_hyperplanes_by_flats,
)
from gpc_census.generation.exterior import exterior_weights
from gpc_census.generation.orbit_canonical import trace_survivor_count
from gpc_census.generation.profiles import (
    ValueProfile,
    capped_patterns,
    enumerate_admissible_profiles,
    incidence_agrees,
    occupancy_patterns,
    profile_is_admissible,
    profile_normal,
    spread_saturation,
    verify_profile_against_weights,
    weight_affine_rank,
    zero_patterns,
)

# docs/orbit_canonical_flats.md, unoriented S_d orbits of admissible hyperplanes
UNORIENTED_ORBITS = {6: 8, 7: 19, 8: 56, 9: 231, 10: 1337}
# docs/orbit_canonical_flats.md, oriented counts and self-paired orbits
ORIENTED_ORBITS = {6: 15, 7: 35, 8: 109}
SELF_PAIRED = {6: 1, 7: 3, 8: 3}
# docs/levi_triangularity_gate.md and docs/bdr_linear_triangular_gate.md,
# the trace-alive oriented orbits, there called mechanisms
ALIVE_ORBITS = {7: 25, 8: 64, 9: 191}
# docs/screening_orbit_structure.md, Ressayre trace survivors
TRACE_SURVIVORS = {6: 64, 7: 549, 8: 6607}

# Spread 14 covers every profile at ranks 6 through 10; see the saturation test.
SPREAD = 14


def _canonical_tau(tau):
    forward = tuple(sorted(tau, reverse=True))
    backward = tuple(sorted((-value for value in tau), reverse=True))
    return min(forward, backward)


def _tau_of_normal(h, z, n):
    tau = tuple(n * value - z for value in h)
    divisor = 0
    for value in tau:
        divisor = _gcd(divisor, abs(value))
    return _canonical_tau(tuple(value // divisor for value in tau))


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


@pytest.fixture(scope="module")
def cache():
    return {}


_PROFILES: dict[tuple[int, int, int], tuple] = {}


def profiles_for(n, d, spread, cache):
    """Enumerate once per (n, d, spread); the suite reuses the result."""
    key = (n, d, spread)
    if key not in _PROFILES:
        _PROFILES[key] = enumerate_admissible_profiles(
            n, d, spread, level_cache=cache
        )
    return _PROFILES[key]


@pytest.mark.parametrize("n,d", [(3, 6), (3, 7), (4, 7)])
def test_matches_reference_flat_enumeration_set_for_set(n, d, cache):
    """Exact set equality against the exhaustive closed-flat enumerator."""
    truth = {
        _tau_of_normal(plane.h, plane.z, n)
        for plane in enumerate_admissible_hyperplanes_by_flats(n, d).hyperplanes
    }
    mine = {
        _canonical_tau(profile.tau)
        for profile in profiles_for(n, d, SPREAD, cache)
    }
    assert mine == truth


@pytest.mark.parametrize("d", sorted(UNORIENTED_ORBITS))
def test_unoriented_orbit_counts(d, cache):
    profiles = profiles_for(3, d, SPREAD, cache)
    orbits = {_canonical_tau(profile.tau) for profile in profiles}
    assert len(orbits) == UNORIENTED_ORBITS[d]


@pytest.mark.parametrize("d", sorted(ORIENTED_ORBITS))
def test_oriented_orbit_and_self_paired_counts(d, cache):
    """One profile is one ORIENTED orbit, so the two counts must agree."""
    profiles = profiles_for(3, d, SPREAD, cache)
    assert len(profiles) == ORIENTED_ORBITS[d]
    self_paired = sum(
        1
        for profile in profiles
        if sum(profile.tau) == 0
        and sorted(-value for value in profile.tau) == sorted(profile.tau)
    )
    assert self_paired == SELF_PAIRED[d]
    assert len(profiles) == 2 * UNORIENTED_ORBITS[d] - self_paired


def _below_and_budget(profile: ValueProfile):
    h, z = profile_normal(profile)
    d = profile.rank
    weights = exterior_weights(profile.particle_number, d)
    below = sum(
        1
        for weight in weights
        if sum(left * right for left, right in zip(h, weight, strict=True)) < z
    )
    budget = d * (d - 1) // 2 - sum(
        value * (value - 1) // 2 for value in profile.multiplicities
    )
    return h, below, budget


@pytest.mark.parametrize("d", sorted(ALIVE_ORBITS))
def test_trace_alive_orbit_counts(d, cache):
    """The Levi-gate mechanism populations, reached without the flat lattice."""
    profiles = profiles_for(3, d, SPREAD, cache)
    alive = 0
    for profile in profiles:
        _, below, budget = _below_and_budget(profile)
        if below <= budget:
            alive += 1
    assert alive == ALIVE_ORBITS[d]


@pytest.mark.parametrize("d", sorted(TRACE_SURVIVORS))
def test_trace_survivor_totals(d, cache):
    """Closed-form survivor totals, which pin B as well as the orbit set."""
    total = 0
    for profile in profiles_for(3, d, SPREAD, cache):
        h, below, _ = _below_and_budget(profile)
        total += trace_survivor_count(h, below)
    assert total == TRACE_SURVIVORS[d]


@pytest.mark.parametrize("n,d", [(3, 6), (3, 7), (3, 8), (3, 9), (4, 7)])
def test_theorem_a_agrees_with_materialised_weight_rank(n, d, cache):
    """Theorem A is checked against the definition it replaces, not assumed."""
    profiles = profiles_for(n, d, SPREAD, cache)
    assert profiles
    for profile in profiles:
        assert verify_profile_against_weights(n, profile)


@pytest.mark.parametrize("n,d", [(3, 6), (3, 7), (4, 7)])
def test_profile_normal_cuts_out_the_same_weights(n, d, cache):
    for profile in profiles_for(n, d, SPREAD, cache):
        assert incidence_agrees(n, profile)


def test_mobility_is_necessary_not_decoration():
    """``tau = (1,1,0,0,0,-2)`` has pattern rank r-1 but weight rank 2, not 5.

    Dropping condition (ii) of Theorem A admits it.  The reference enumerator
    does not, so this row is a standing counterexample to the rank-only test.
    """
    multiplicities, values = (2, 3, 1), (1, 0, -2)
    patterns = zero_patterns(3, multiplicities, values)
    assert set(patterns) == {(0, 3, 0), (2, 0, 1)}
    # pattern rank is r - 1 = 2, so the rank half of Theorem A passes
    from gpc_census.generation.profiles import _rank_mod

    assert _rank_mod(patterns, 3) == 2
    # but class 0 has no pattern with 1 <= c_0 <= m_0 - 1 = 1
    assert not any(1 <= pattern[0] <= 1 for pattern in patterns)
    assert not profile_is_admissible(3, multiplicities, values)
    assert weight_affine_rank(3, (1, 1, 0, 0, 0, -2)) == 2


def test_lemma_l_value_gap_divides_the_particle_number():
    """``B | N`` falls out of Lemma L, so N=3 profiles have value gap 1 or 3."""
    for profile in enumerate_admissible_profiles(3, 8, SPREAD):
        values = profile.values
        gap = 0
        for left, right in itertools.combinations(values, 2):
            gap = _gcd(gap, abs(left - right))
        assert gap in (1, 3)
        # and every zero pattern shares one level sum, which is what names it
        patterns = zero_patterns(3, profile.multiplicities, values)
        sums = {
            sum(
                count * level
                for count, level in zip(pattern, profile.levels, strict=True)
            )
            for pattern in patterns
        }
        assert sums == {profile.level_target}


@pytest.mark.parametrize("d,realized", [(6, 3), (7, 4), (8, 6), (9, 10), (10, 14)])
def test_spread_saturation_plateaus_at_the_known_orbit_count(d, realized, cache):
    """The saturation protocol is calibrated where the truth is known.

    Raising the bound past the realized maximum must add nothing, and the
    plateau value must be the independently known orbit count.
    """
    rows = spread_saturation(3, d, (realized, realized + 1), level_cache=cache)
    assert rows[0]["profiles"] == rows[1]["profiles"]
    assert rows[0]["realized_max_spread"] == realized
    assert rows[1]["bound_is_slack"]
    orbits = {
        _canonical_tau(profile.tau) for profile in profiles_for(3, d, SPREAD, cache)
    }
    assert len(orbits) == UNORIENTED_ORBITS[d]


@pytest.mark.parametrize(
    "d,budget,found,truth",
    # REFUTATION: the bounded-hole arithmetic-level grammar is not a
    # candidate-completeness rule.  A hole budget h means spread <= r - 1 + h.
    [(8, 0, 47, 56), (9, 1, 210, 231), (9, 3, 231, 231), (10, 3, 1244, 1337)],
)
def test_bounded_hole_grammar_is_not_complete(d, budget, found, truth, cache):
    """docs/arithmetic_level_audit.md counts holes on FACETS; candidates differ.

    At most one hole holds on the mechanism population and fails badly on the
    candidate universe.  These four rows are the counterexample, kept so the
    grammar cannot be promoted to a pruning rule by accident.
    """
    orbits = {
        _canonical_tau(profile.tau)
        for profile in profiles_for(3, d, SPREAD, cache)
        if profile.spread <= profile.classes - 1 + budget
    }
    assert len(orbits) == found
    assert UNORIENTED_ORBITS[d] == truth
    if found != truth:
        assert found < truth


def test_pattern_helpers_are_capped_correctly():
    assert occupancy_patterns(3, 2) == ((0, 3), (1, 2), (2, 1), (3, 0))
    assert capped_patterns(3, (1, 4)) == ((0, 3), (1, 2))
    assert profile_is_admissible(3, (1, 3, 2), (1, 0, -1))
    with pytest.raises(ValueError):
        profile_is_admissible(3, (1, 2), (1, 1))


# --- guard for the shipped artifact -----------------------------------------

ARTIFACT = (
    pathlib.Path(__file__).resolve().parents[1]
    / "results"
    / "data"
    / "generator_frontier_refinement.json"
)


REALIZED_SPREAD_IN_ARTIFACT = {6: 3, 7: 4, 8: 6, 9: 10, 10: 14}


def _artifact() -> dict:
    if not ARTIFACT.is_file():
        pytest.skip("frontier artifact absent; run scripts/profile_native_frontier.py")
    return json.loads(ARTIFACT.read_text())


def test_generator_frontier_refinement_calibration_is_exact():
    """The shipped artifact must agree with every known population."""
    artifact = _artifact()
    assert artifact["calibration_is_exact"]
    assert artifact["calibration_mismatches"] == []
    by_rank = {row["rank"]: row for row in artifact["calibration"]}
    for d, expected in UNORIENTED_ORBITS.items():
        assert by_rank[d]["unoriented_orbits"] == expected
    for d, expected in ORIENTED_ORBITS.items():
        assert by_rank[d]["oriented_orbits"] == expected
    for d, expected in ALIVE_ORBITS.items():
        assert by_rank[d]["trace_alive_orbits"] == expected
    for d, expected in TRACE_SURVIVORS.items():
        assert by_rank[d]["trace_survivors"] == expected
    for d, realized in REALIZED_SPREAD_IN_ARTIFACT.items():
        assert by_rank[d]["realized_max_spread"] == realized
    # Theorem A is checked against the materialised definition, never assumed.
    for row in artifact["calibration"] + artifact["cross_particle_number"]:
        if "theorem_a_disagreements" in row:
            assert row["theorem_a_disagreements"] == 0
            assert row["theorem_a_checked"] > 0


def test_generator_frontier_refinement_keeps_its_labels_honest():
    """The spread bound stays UNRESOLVED and rank 11+ stays CONDITIONAL."""
    labels = _artifact()["proof_labels"]
    assert labels["spread_bound"].startswith("UNRESOLVED")
    assert labels["bounded_hole_grammar_as_completeness_rule"] == "REFUTED"
    assert labels["rank_11_to_13_candidate_universe"].startswith("CONDITIONAL")
    assert _artifact()["frontier_is_conditional_on_spread_bound"]


def test_generator_frontier_refinement_records_the_hole_refutation():
    rows = {
        (row["rank"], row["hole_budget"]): row
        for row in _artifact()["bounded_hole_grammar"]
    }
    assert rows[(10, 1)]["missing"] == 499
    assert rows[(10, 3)]["missing"] == 93
    assert rows[(9, 1)]["missing"] == 21
    assert rows[(8, 0)]["missing"] == 9
    assert not rows[(10, 4)]["complete"]


def test_generator_frontier_refinement_records_the_rank_eleven_non_saturation():
    """Rank 11 is the first rank where the plateau test fails; keep it visible."""
    artifact = _artifact()
    ladder = artifact.get("rank_eleven_spread_ladder", [])
    if not ladder:
        return  # artifact written without --frontier
    assert not artifact["rank_eleven_saturates"]
    counts = [row["profiles"] for row in ladder]
    assert counts == sorted(counts)
    assert counts[-1] > counts[0]
    # the realized maximum reaches the bound rather than sitting inside it
    assert not ladder[-1]["bound_is_slack"]
    frontier = {row["rank"]: row for row in artifact["frontier"]}
    assert frontier[11]["max_spread"] == ladder[-1]["max_spread"]
    assert not frontier[11]["spread_bound_is_slack"]


def test_profile_native_survivor_taus_match_the_lattice_route():
    """The two orbit sources must feed the pipeline identical rows.

    ``orbit_native`` walks the closed-flat lattice; ``profiles`` writes the
    orbits down.  Everything downstream is shared, so equality of the emitted
    tau SETS is the end-to-end check that the replacement is a replacement.
    """
    from gpc_census.generation.orbit_native import enumerate_orbit_survivor_taus
    from gpc_census.generation.profiles import enumerate_profile_survivor_taus

    for d, expected in ((6, 64), (7, 549), (8, 6607)):
        lattice = enumerate_orbit_survivor_taus(3, d)
        mine, counts = enumerate_profile_survivor_taus(3, d, SPREAD)
        assert counts["predicted_survivors"] == counts["produced_survivors"]
        assert len(mine) == expected
        assert set(mine) == set(lattice.survivor_taus)
        assert counts["oriented_orbits"] == ORIENTED_ORBITS[d]
