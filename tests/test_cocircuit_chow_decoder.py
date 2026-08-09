"""Guards for direct cocircuit enumeration and the signed-Chow decoder.

The live checks run at `(3,5)` and `(3,6)`, where a full recomputation costs
about a second. Ranks 7 and 8 are read from the artifact, which is what the
generating script and the standalone verifier exist to produce.
"""

from __future__ import annotations

import itertools
import json
import pathlib
import random
from fractions import Fraction

import pytest

from gpc_census.generation.admissible_flats import (
    enumerate_admissible_hyperplanes_by_flats,
)
from gpc_census.generation.chow_decoder import (
    DecoderBudgetError,
    decode_shadow,
    decode_signed_pair,
    degree_blocks,
    dominates,
    profile_contribution,
    shadow_of,
    threshold_family,
)
from gpc_census.generation.cocircuit import (
    NodeLimitExceeded,
    enumerate_cocircuit_hyperplanes,
    trace_survivors,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
DOCS = ROOT / "docs"
ARTIFACT = DATA / "cocircuit_chow_decoder.json"

# Sealed populations, restated here so a regression is caught by this file and
# not only by comparison with the artifact it also writes.
SEALED = {5: (30, 154, 10, 8), 6: (362, 2775, 64, 62)}


def _artifact() -> dict:
    return json.loads(ARTIFACT.read_text())


def _orbit_key(tau):
    return min(tuple(sorted(tau)), tuple(sorted(-value for value in tau)))


# --------------------------------------------------------------- enumeration


@pytest.mark.parametrize("d", sorted(SEALED))
def test_cocircuit_reproduces_the_sealed_populations(d):
    hyperplanes, nodes, survivors, nonstructural = SEALED[d]
    enumeration = enumerate_cocircuit_hyperplanes(3, d)
    assert enumeration.hyperplane_count == hyperplanes
    assert enumeration.node_count == nodes
    rows = list(trace_survivors(enumeration))
    assert len(rows) == survivors
    assert sum(1 for _tau, positive in rows if positive > 0) == nonstructural


@pytest.mark.parametrize("d", sorted(SEALED))
def test_cocircuit_agrees_with_the_flat_lattice_pair_for_pair(d):
    """The convention bridge is checked, not assumed.

    `docs/orbit_canonical_flats.md` records that a missing sign normalization
    answers a different and equally plausible question, so counts alone are not
    enough: the `(h, z)` pairs themselves must coincide.
    """
    enumeration = enumerate_cocircuit_hyperplanes(3, d)
    reference = enumerate_admissible_hyperplanes_by_flats(3, d)
    assert set(enumeration.affine_pairs()) == {
        (row.h, row.z) for row in reference.hyperplanes
    }


def test_cocircuit_orbits_match_the_orbit_canonical_artifact():
    """Orbit counts from bare coordinate multisets, no canonicalizer involved."""
    recorded = json.loads((DATA / "orbit_canonical_flats.json").read_text())
    checked = 0
    for d in sorted(SEALED):
        # (3,5) is below the orbit enumerator's recorded range; it is gated by
        # the standalone brute force instead.
        if str(d) not in recorded["augmentation"]:
            continue
        enumeration = enumerate_cocircuit_hyperplanes(3, d)
        orbits = {_orbit_key(tau) for tau in enumeration.normals}
        assert len(orbits) == recorded["augmentation"][str(d)]["top_orbits"]
        checked += 1
    assert checked


def test_every_emitted_zero_set_is_exactly_the_normals_zero_set():
    """The modular search is a device; the shipped statement is exact."""
    enumeration = enumerate_cocircuit_hyperplanes(3, 6)
    weights = list(itertools.combinations(range(6), 3))
    for tau, zero_set in zip(
        enumeration.normals, enumeration.zero_sets, strict=True
    ):
        exact = tuple(
            index
            for index, subset in enumerate(weights)
            if sum(tau[i] for i in subset) == 0
        )
        assert exact == zero_set
        assert len(zero_set) >= 5


def test_node_limit_raises_instead_of_truncating():
    """A partial reverse search must never be mistaken for a population."""
    with pytest.raises(NodeLimitExceeded):
        enumerate_cocircuit_hyperplanes(3, 6, node_limit=100)


# ------------------------------------------------------------------- T4 / T5


def test_t4_equal_degrees_force_symmetry():
    """Equal `c_+` coordinates make the transposition fix `B_+` itself."""
    rng = random.Random(404)
    n, d = 3, 7
    checked = 0
    for _ in range(60):
        tau = tuple(rng.randrange(-6, 7) for _ in range(d))
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
                checked += 1
    assert checked > 100


def test_t5_block_average_is_a_separator_with_ordered_levels():
    """Averaging over the blocks keeps `B_+` and strictly orders the levels."""
    rng = random.Random(505)
    n, d = 3, 7
    checked = 0
    for _ in range(60):
        tau = tuple(rng.randrange(-6, 7) for _ in range(d))
        family = threshold_family(tau, n)
        if not family:
            continue
        blocks, _levels = degree_blocks(shadow_of(family, d))
        averaged = [Fraction(0)] * d
        for block in blocks:
            mean = Fraction(sum(tau[index] for index in block), len(block))
            for index in block:
                averaged[index] = mean
        rebuilt = frozenset(
            subset
            for subset in itertools.combinations(range(d), n)
            if sum(averaged[index] for index in subset) > 0
        )
        assert rebuilt == family
        levels = [averaged[block[0]] for block in blocks]
        assert all(
            levels[a] > levels[a + 1] for a in range(len(levels) - 1)
        )
        checked += 1
    assert checked > 20


def test_t5_selected_profiles_form_an_upper_ideal():
    """Nothing dominating a selected profile may be rejected."""
    rng = random.Random(606)
    n, d = 3, 8
    for _ in range(40):
        tau = tuple(rng.randrange(-6, 7) for _ in range(d))
        family = threshold_family(tau, n)
        if not family:
            continue
        blocks, _levels = degree_blocks(shadow_of(family, d))

        def profile(subset, blocks=blocks):
            return tuple(len(set(subset) & set(block)) for block in blocks)

        selected = {profile(subset) for subset in family}
        rejected = {
            profile(subset) for subset in itertools.combinations(range(d), n)
        } - selected
        for chosen in selected:
            for other in rejected:
                assert not dominates(other, chosen)


def test_profile_contribution_counts_what_it_claims():
    """`q_a(k)` is the number of profile-`k` determinants through one orbital."""
    blocks = ((0, 1, 2), (3, 4), (5,))
    sizes = tuple(len(block) for block in blocks)
    for profile in itertools.product(range(4), repeat=3):
        if sum(profile) != 3 or any(
            taken > size for taken, size in zip(profile, sizes)
        ):
            continue
        expected = []
        for block in blocks:
            fixed = block[0]
            expected.append(
                sum(
                    1
                    for subset in itertools.combinations(range(6), 3)
                    if fixed in subset
                    and tuple(
                        len(set(subset) & set(other)) for other in blocks
                    )
                    == profile
                )
            )
        assert profile_contribution(profile, sizes) == tuple(expected)


# --------------------------------------------------------------- the decoder


def test_decoder_inverts_every_survivor_at_rank_six():
    """The full `(3,6)` nonstructural population, decoded and reconstructed."""
    enumeration = enumerate_cocircuit_hyperplanes(3, 6)
    rows = 0
    for tau, positive in trace_survivors(enumeration):
        if positive == 0:
            continue
        rows += 1
        plus = shadow_of(threshold_family(tau, 3, +1), 6)
        minus = shadow_of(threshold_family(tau, 3, -1), 6)
        decoded = decode_signed_pair(plus, minus, 3)
        assert decoded is not None
        assert decoded.positive.family == threshold_family(tau, 3, +1)
        assert decoded.negative.family == threshold_family(tau, 3, -1)
        assert decoded.zero_span_rank == 5
        assert decoded.tau == tuple(tau)
    assert rows == 62


def test_decoder_inverts_random_threshold_shadows():
    """Random normals at ranks 7 and 9, decoded from the bare shadow."""
    rng = random.Random(707)
    for d in (7, 9):
        checked = 0
        for _ in range(25):
            tau = tuple(rng.randrange(-6, 7) for _ in range(d))
            family = threshold_family(tau, 3)
            if not family:
                continue
            decoded = decode_shadow(shadow_of(family, d), 3)
            assert decoded is not None
            assert decoded.family == family
            checked += 1
        assert checked > 10


def test_decoder_returns_none_on_an_unrealizable_shadow():
    """A vector that no family induces gets `None`, not a guess.

    The last two are the interesting ones: they survive the cheap sum and
    sign checks and are rejected by the profile equations themselves.
    """
    assert decode_shadow((1, 0, 0, 0, 0, 0), 3) is None
    assert decode_shadow((-1, 2, 2, 2, 2, 1), 3) is None
    assert decode_shadow((3, 0, 0, 0, 0, 0), 3) is None
    assert decode_shadow((4, 4, 4, 4, 4, 4), 3) is None
    # A single block admits only the empty family and the whole slice.
    whole = decode_shadow((10,) * 6, 3)
    assert whole is not None and len(whole.family) == 20


def test_decoder_budget_is_enforced():
    tau = (5, 4, 3, 2, 1, 0, -1, -2, -3)
    shadow = shadow_of(threshold_family(tau, 3), 9)
    with pytest.raises(DecoderBudgetError):
        decode_shadow(shadow, 3, max_states=1)


def test_decoded_family_reproduces_the_shadow_it_was_given():
    """The decoder's own postcondition, checked independently of T1."""
    rng = random.Random(808)
    for _ in range(30):
        tau = tuple(rng.randrange(-5, 6) for _ in range(8))
        family = threshold_family(tau, 3)
        if not family:
            continue
        target = shadow_of(family, 8)
        decoded = decode_shadow(target, 3)
        assert decoded is not None
        assert shadow_of(decoded.family, 8) == target


# --------------------------------------------------------------- the artifact


def test_artifact_records_the_full_populations():
    artifact = _artifact()
    expected = {
        "5": (30, 154, 10, 8, None),
        "6": (362, 2775, 64, 62, 8),
        "7": (5341, 63136, 549, 547, 19),
        "8": (166420, 2042570, 6607, 6605, 56),
    }
    for key, (planes, nodes, survivors, nonstructural, orbits) in expected.items():
        row = artifact["cocircuit"][key]
        assert row["admissible_hyperplanes"] == planes
        assert row["reverse_search_nodes"] == nodes
        assert row["trace_survivors"] == survivors
        assert row["nonstructural_trace_survivors"] == nonstructural
        assert row["flat_lattice_hyperplanes"] == planes
        assert row["agrees_with_flat_lattice_pairwise"] is True
        assert row["agrees_with_centered_convention_test"] is True
        if orbits is not None:
            assert row["sd_unoriented_orbits"] == orbits


def test_artifact_populations_agree_with_the_sealed_series():
    """The same counts as the gate-2 reference series and the orbit ladder."""
    artifact = _artifact()["cocircuit"]
    series = json.loads((DATA / "stage1_reference_series.json").read_text())
    orbits = json.loads((DATA / "orbit_canonical_flats.json").read_text())
    for key, row in artifact.items():
        recorded = orbits["augmentation"].get(key)
        if recorded is not None:
            assert row["admissible_hyperplanes"] == recorded["hyperplanes"]
            assert row["sd_unoriented_orbits"] == recorded["top_orbits"]
    for key, entry in series["ranks"].items():
        if key in artifact:
            assert artifact[key]["admissible_hyperplanes"] == (
                entry["admissible_hyperplanes"]
            )
            assert artifact[key]["admissible_hyperplanes"] == (
                entry["flat_counts_by_rank"][-1]
            )


def test_artifact_decoder_has_zero_failures():
    artifact = _artifact()["decoder"]
    total = 0
    for row in artifact.values():
        assert row["decode_failures"] == 0
        assert row["tau_reconstruction_failures"] == 0
        total += row["rows"]
    assert total == 7214
    assert artifact["8"]["max_search_states_per_side"] < 100


def test_artifact_stress_is_labelled_and_grows():
    """The stress block must show growth, since no bound is claimed."""
    stress = _artifact()["decoder_stress"]
    assert all(row["failures"] == 0 for row in stress.values())
    ordered = [stress[key] for key in sorted(stress, key=int)]
    assert ordered[-1]["max_search_states"] > ordered[0]["max_search_states"]


def test_artifact_keeps_proved_and_measured_apart():
    artifact = _artifact()
    assert set(artifact["theorems"]) == {
        "C1_cocircuit_equivalence",
        "T4_block_symmetry",
        "T5_block_constant_separator",
    }
    measured = artifact["measured_not_proved"]
    assert "decoder_complexity" in measured
    assert "symmetric_cocircuit_enumeration" in measured
    assert "downstream_gates_untouched" in measured


def test_cost_comparison_records_a_machine_and_all_three_routes():
    artifact = _artifact()
    assert artifact["machine"]["platform"]
    cost = artifact["cost_comparison"]
    for field in (
        "cocircuit_seconds",
        "flat_lattice_seconds",
        "orbit_augmentation_seconds",
    ):
        assert cost[field] > 0
    assert cost["orbit_augmentation_hyperplanes"] == cost["cocircuit_hyperplanes"]
    assert (
        cost["orbit_augmentation_stored_objects"]
        < cost["flat_lattice_stored_objects"]
    )


# ------------------------------------------------------------------ the docs


@pytest.mark.parametrize(
    "phrase",
    [
        "not a replacement",
        "no polynomial worst-case bound",
        "regression targets, not new information",
        "Facetness, redundancy, birationality, completeness",
    ],
)
def test_the_document_states_what_is_not_claimed(phrase):
    assert phrase in (DOCS / "cocircuit_chow_decoder.md").read_text()


def test_the_signed_chow_document_records_that_decoding_is_now_solved():
    """The caveat that motivated this work must not still read as open."""
    text = (DOCS / "signed_chow_projection.md").read_text()
    assert "docs/cocircuit_chow_decoder.md" in text


def test_prereg_is_scored_against_the_artifact():
    text = (DOCS / "prereg_cocircuit_chow_decoder.md").read_text()
    assert "## SCORED" in text
    assert "INFORMED, not blind" in text
    for identifier in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"):
        assert f"| {identifier} |" in text
