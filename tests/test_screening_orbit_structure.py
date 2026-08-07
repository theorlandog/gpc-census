"""Guards for the refutation of screening orbit-invariance.

The write-up leans on a NEGATIVE result, so the tests have to keep it negative:
a later change that made screening look orbit-invariant would silently convert
a recorded refutation into a confirmation. They also pin the decomposition that
survives, since that is what the generator's cost model now rests on.
"""
import json
import pathlib

import pytest

from gpc_census.generation import admissible_flats as af
from gpc_census.generation.exterior import exterior_weights
from gpc_census.generation.ressayre import _negative_root_set

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "screening_orbit_structure.json"


def _art():
    return json.loads(ARTIFACT.read_text())


# --------------------------------------------------------------------------
# the refutation must stay refuted
# --------------------------------------------------------------------------

def test_screening_is_not_constant_on_orbits():
    """Orbits carry more than one inversion count, so no verdict transports."""
    art = _art()
    mixed = art["summary"]["orbits_with_mixed_inversion_counts"]
    assert mixed["7"] > 0
    assert mixed["7"] >= 20, "the refutation is supposed to be broad, not marginal"
    for record in art["ranks"].values():
        many = [r for r in record["orbit_records"]
                if len(r["inversion_counts"]) > 1]
        assert many, record["d"]


def test_inversion_counts_really_do_spread_across_an_orbit():
    """Recomputed live: permuting h moves the root count and fixes B."""
    import itertools

    n, d = 3, 6
    weights = exterior_weights(n, d)
    occupied = [tuple(i for i, bit in enumerate(w) if bit) for w in weights]
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    hp = enum.hyperplanes[100]

    below, roots = set(), set()
    for perm in itertools.permutations(range(d)):
        image = tuple(hp.h[perm[i]] for i in range(d))
        below.add(sum(1 for members in occupied
                      if sum(image[i] for i in members) < hp.z))
        roots.add(len(_negative_root_set(image)))
    assert len(below) == 1, "weights-below must be an orbit invariant"
    assert len(roots) > 1, "the root count must NOT be an orbit invariant"


# --------------------------------------------------------------------------
# the decomposition that survives
# --------------------------------------------------------------------------

def test_weights_below_is_an_orbit_invariant_on_every_row():
    art = _art()
    assert art["summary"]["weights_below_invariant_everywhere"]
    for record in art["ranks"].values():
        assert record["weights_below_failures"] == 0, record["d"]
        assert record["weights_below_checks"] == record["oriented_rows"]


def test_the_root_count_is_exactly_the_inversion_number():
    """Which is why it cannot be invariant, and why the row cost is O(d^2)."""
    def inversions(h):
        return sum(1 for i in range(len(h)) for j in range(i + 1, len(h))
                   if h[j] < h[i])

    for d in (6, 7):
        enum = af.enumerate_admissible_hyperplanes_by_flats(3, d)
        for hp in enum.hyperplanes[:300]:
            assert len(_negative_root_set(hp.h)) == inversions(hp.h)


def test_survivor_fractions_are_recorded_and_falling():
    art = _art()
    fractions = art["summary"]["survivor_fraction"]
    assert fractions["6"] > fractions["7"] > fractions["8"]
    assert fractions["8"] < 0.03


def test_rank_7_survivor_count_matches_the_published_partition():
    """10,682 rows minus 10,133 trace failures is 549, per the methodology doc."""
    art = _art()
    assert art["ranks"]["7"]["oriented_rows"] == 10682
    assert art["ranks"]["7"]["trace_survivors"] == 549


def test_orbit_counts_are_the_oriented_ones():
    """Each hyperplane gives two oriented tau rows, so tau orbits exceed the
    unoriented hyperplane orbits and are what screening actually consumes."""
    art = _art()
    orbits = art["summary"]["orbits"]
    assert orbits["6"] == 15
    assert orbits["7"] == 35
    assert orbits["8"] == 109
    # strictly between the unoriented count and twice it, because some orbits
    # are self-paired under orientation reversal
    for unoriented, key in ((8, "6"), (19, "7"), (56, "8")):
        assert unoriented < orbits[key] <= 2 * unoriented


@pytest.mark.parametrize("phrase", ["REFUTED", "fundamental domain"])
def test_the_write_up_states_the_refutation(phrase):
    doc = (pathlib.Path(__file__).resolve().parents[1]
           / "docs" / "screening_orbit_structure.md")
    if not doc.exists():
        pytest.skip("docs not present")
    assert phrase in doc.read_text()


def test_the_earlier_doc_carries_its_correction():
    """orbit_canonical_flats.md asserted orbit-invariance of validity."""
    doc = (pathlib.Path(__file__).resolve().parents[1]
           / "docs" / "orbit_canonical_flats.md")
    if not doc.exists():
        pytest.skip("docs not present")
    text = doc.read_text()
    assert "CORRECTION" in text
    assert "screening_orbit_structure" in text


# --------------------------------------------------------------------------
# the determinant stage: structural zeros decided by matching
# --------------------------------------------------------------------------

def test_structural_zeros_replay_and_a_sample_agrees_with_the_symbolic_route():
    """A Hall violator PROVES the determinant is identically zero.

    Every violator is replayed from the support alone, which is cheap, and a
    sample is cross-checked against the exponential symbolic route, which is
    not. The sample is what would catch a matching bug; replaying all of them
    is what would catch a certificate that merely echoed its own run.
    """
    from gpc_census.generation.reference_generator import oriented_taus
    from gpc_census.generation.ressayre import (
        _level_sets,
        admissible_candidate_from_tau,
        tangent_determinant_is_identically_zero,
        tangent_determinant_is_structurally_zero,
        verify_hall_violator,
    )

    n, d = 3, 7
    weights = exterior_weights(n, d)
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    structural = matching = cross_checked = 0
    for tau in oriented_taus(enum):
        candidate = admissible_candidate_from_tau(n, tau)
        _, below = _level_sets(weights, candidate.h, candidate.z)
        roots = _negative_root_set(candidate.h)
        if len(below) != len(roots) or not roots:
            continue
        is_structural, violator = tangent_determinant_is_structurally_zero(
            n, d, candidate)
        if not is_structural:
            matching += 1
            continue
        structural += 1
        assert verify_hall_violator(n, d, candidate, violator)
        if cross_checked < 40:
            assert tangent_determinant_is_identically_zero(n, d, candidate)
            cross_checked += 1
    assert structural == 474
    assert matching == 73, "the test must be able to say NO, not only YES"
    assert cross_checked == 40


def test_a_certified_row_is_never_called_structurally_zero():
    """The direction that would break the pipeline if it were wrong.

    A false structural verdict would silently discard a genuine facet, so
    every certified row in a prefix of the candidate list is checked to come
    back NOT structural.
    """
    from gpc_census.generation.reference_generator import oriented_taus
    from gpc_census.generation.ressayre import (
        _level_sets,
        admissible_candidate_from_tau,
        certify_candidate,
        tangent_determinant_is_structurally_zero,
    )

    n, d = 3, 7
    weights = exterior_weights(n, d)
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    certified = 0
    for tau in oriented_taus(enum)[:1500]:
        candidate = admissible_candidate_from_tau(n, tau)
        _, below = _level_sets(weights, candidate.h, candidate.z)
        roots = _negative_root_set(candidate.h)
        if len(below) != len(roots) or not roots:
            continue
        if certify_candidate(n, d, candidate) is None:
            continue
        certified += 1
        is_structural, violator = tangent_determinant_is_structurally_zero(
            n, d, candidate)
        assert not is_structural
        assert violator is None
    assert certified > 0
