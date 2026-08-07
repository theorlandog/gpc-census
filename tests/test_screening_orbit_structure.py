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
