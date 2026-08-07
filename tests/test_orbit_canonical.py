"""Guards for Weyl-orbit canonicalization and the augmentation enumerator.

Orbit reduction is the one pruning the standing rule permits, so the code that
performs it has to be verified rather than trusted. The tests that matter are
the ones that could catch a canonical form that is not actually invariant, or
an enumerator that silently loses orbits:

- the canonical form is checked against orbit counts computed independently by
  union-find over adjacent transpositions;
- orbit sizes from the orbit-stabilizer theorem must SUM to the materialised
  flat counts at every level, which is the check that would fail if the
  enumerator dropped an orbit;
- the sign-convention trap that gave 15 orbits instead of 8 is pinned.
"""
import json
import pathlib

import pytest

from gpc_census.generation import admissible_flats as af
from gpc_census.generation import orbit_canonical as oc

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "orbit_canonical_flats.json"


def _art():
    return json.loads(ARTIFACT.read_text())


def _hyperplane_masks(n, d):
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    subsets = oc.weight_subsets(n, d)
    masks = []
    for hp in enum.hyperplanes:
        mask = 0
        for i, subset in enumerate(subsets):
            if sum(hp.h[m] for m in subset) == hp.z:
                mask |= 1 << i
        masks.append(mask)
    return masks


# --------------------------------------------------------------------------
# the canonical form
# --------------------------------------------------------------------------

def test_canonical_form_is_permutation_invariant():
    import random

    n, d = 3, 6
    subsets = oc.weight_subsets(n, d)
    index = {s: i for i, s in enumerate(subsets)}
    random.seed(20260806)
    for _ in range(100):
        mask = random.getrandbits(len(subsets))
        perm = list(range(d))
        random.shuffle(perm)
        moved = 0
        for i, subset in enumerate(subsets):
            if mask >> i & 1:
                moved |= 1 << index[tuple(sorted(perm[m] for m in subset))]
        assert oc.canonical_form(mask, n, d) == oc.canonical_form(moved, n, d)


@pytest.mark.parametrize(("d", "orbits"), [(6, 8), (7, 19)])
def test_canonical_forms_reproduce_the_union_find_orbit_count(d, orbits):
    masks = _hyperplane_masks(3, d)
    assert len({oc.canonical_form(m, 3, d) for m in masks}) == orbits


@pytest.mark.parametrize(("d", "total"), [(6, 362), (7, 5341)])
def test_orbit_sizes_from_stabilizers_sum_to_the_total(d, total):
    """Orbit-stabilizer, checked against the materialised count.

    If the canonicalization merged two genuinely distinct orbits, or split one,
    this sum would not land on the hyperplane count.
    """
    masks = _hyperplane_masks(3, d)
    seen = {}
    for mask in masks:
        seen.setdefault(oc.canonical_form(mask, 3, d), mask)
    assert sum(oc.orbit_size(m, 3, d) for m in seen.values()) == total


def test_canonicalization_cost_is_bounded_by_the_automorphism_group():
    """The factorial tail must stay dead.

    Plain colour refinement followed by brute force inside the cells hit
    exactly ``d!`` on the most symmetric flats, which is the wrong shape and is
    what stalled rank 11. Individualization brings the leaf count down to the
    order of the automorphism group, which is a floor: every automorphism
    produces a distinct leaf, so no canonicalizer of this form can do better
    than 1x. This pins the ratio rather than a raw timing, so it is machine
    independent, and it would fail loudly if the individualization step were
    removed.
    """
    leaves = [0]
    real_apply = oc._apply

    def counting_apply(present, mapping, n, d):
        leaves[0] += 1
        return real_apply(present, mapping, n, d)

    d = 8
    _, levels = oc.enumerate_orbits_by_augmentation(3, d)
    oc._apply = counting_apply
    try:
        worst_ratio = 0.0
        worst_leaves = 0
        for level in levels[1:]:
            for mask in level:
                leaves[0] = 0
                _, automorphisms = oc.canonical_form(
                    mask, 3, d, with_automorphisms=True)
                worst_ratio = max(worst_ratio, leaves[0] / automorphisms)
                worst_leaves = max(worst_leaves, leaves[0])
    finally:
        oc._apply = real_apply
    assert worst_ratio <= 4.0, worst_ratio
    # the pre-individualization implementation hit 40320 here
    assert worst_leaves < oc._factorial(d) // 4, worst_leaves


# --------------------------------------------------------------------------
# the augmentation enumerator
# --------------------------------------------------------------------------

@pytest.mark.parametrize("d", [6, 7])
def test_augmentation_reproduces_the_full_lattice_at_every_level(d):
    """The decisive check: expanding orbits must rebuild the whole lattice.

    Not just the top level. An enumerator that lost an orbit at some
    intermediate rank would still land on the right hyperplane count if the
    loss happened to be reachable another way, so every level is compared.
    """
    counts, levels = oc.enumerate_orbits_by_augmentation(3, d)
    expanded = [sum(oc.orbit_size(m, 3, d) for m in level) for level in levels]
    materialised = list(
        af.enumerate_admissible_hyperplanes_by_flats(3, d).flat_counts_by_rank)
    assert expanded == materialised
    assert counts[-1] == len({oc.canonical_form(m, 3, d)
                              for m in _hyperplane_masks(3, d)})


def test_checkpointed_and_resumed_runs_agree_with_an_uninterrupted_one(tmp_path):
    """A restart must not change the answer, and must not silently skip a level.

    A (3,11) run was lost outright to a machine restart after reaching rank 8,
    so resumption is now part of the enumerator. The failure mode worth testing
    is not "does it crash" but "does a resumed run agree with a clean one",
    since a checkpoint written or read at the wrong rank would still produce a
    plausible ladder.
    """
    reference, _ = oc.enumerate_orbits_by_augmentation(3, 6)

    first, _ = oc.enumerate_orbits_by_augmentation(3, 6, checkpoint_dir=tmp_path)
    assert first == reference

    # every level is now on disk, so this run restores rather than computes
    resumed, levels = oc.enumerate_orbits_by_augmentation(
        3, 6, checkpoint_dir=tmp_path)
    assert resumed == reference
    expanded = [sum(oc.orbit_size(m, 3, 6) for m in level) for level in levels]
    assert expanded == list(
        af.enumerate_admissible_hyperplanes_by_flats(3, 6).flat_counts_by_rank)

    # a partial run: drop the tail and confirm it recomputes to the same place
    saved = sorted((tmp_path / "orbits").glob("*.txt"))
    assert len(saved) == 5, [p.name for p in saved]
    for path in saved[-2:]:
        path.unlink()
    partial, _ = oc.enumerate_orbits_by_augmentation(
        3, 6, checkpoint_dir=tmp_path)
    assert partial == reference


def test_checkpoints_do_not_collide_with_the_materialising_enumerator(tmp_path):
    """The level filenames are shared, so orbit levels live in a subdirectory.

    Reading a materialised level as if it were a set of representatives would
    put the whole lattice in memory, which is exactly what the orbit enumerator
    exists to avoid, and it would not look like an error.
    """
    oc.enumerate_orbits_by_augmentation(3, 6, checkpoint_dir=tmp_path)
    assert not list(tmp_path.glob("*.txt"))
    assert list((tmp_path / "orbits").glob("*.txt"))


@pytest.mark.parametrize(("d", "orbits", "total"), [
    (6, 8, 362), (7, 19, 5341), (8, 56, 166420)])
def test_the_codimension_one_key_agrees_with_the_hypergraph_canonicalizer(
        d, orbits, total):
    """The top level is decided by a sort, and it must decide the same way.

    Instrumenting (3,9) showed the top level costs 167 of 214 seconds, 78
    percent, because codimension-one flats carry most of the weight set and so
    have the largest automorphism groups. Replacing the search with the exact
    normal is only legitimate if it induces the SAME partition and the SAME
    orbit sizes, so both are checked against the expensive route rather than
    against a stored number.
    """
    _, levels = oc.enumerate_orbits_by_augmentation(3, d)
    top = levels[-1]
    assert len(top) == orbits
    cheap_total = 0
    for mask, witness in top.items():
        h, z = oc._exact_top_normal(3, d, witness)
        cheap = oc.top_level_orbit_size(h, z, d)
        assert cheap == oc.orbit_size(mask, 3, d), (h, z)
        cheap_total += cheap
    assert cheap_total == total

    # distinct orbits must get distinct keys, or the enumerator would merge them
    keys = {oc.top_level_canonical_key(*oc._exact_top_normal(3, d, w))
            for w in top.values()}
    assert len(keys) == orbits


def test_the_codimension_one_key_is_permutation_and_sign_invariant():
    """Both invariances matter: the flat is unoriented, so (h, z) and (-h, -z)
    describe the same flat and must not be given different keys."""
    import random

    random.seed(20260808)
    for _ in range(200):
        d = random.randint(3, 9)
        h = tuple(random.randint(-4, 4) for _ in range(d))
        z = random.randint(-5, 5)
        permuted = list(h)
        random.shuffle(permuted)
        assert (oc.top_level_canonical_key(h, z)
                == oc.top_level_canonical_key(tuple(permuted), z))
        assert (oc.top_level_canonical_key(h, z)
                == oc.top_level_canonical_key(tuple(-v for v in h), -z))


def test_orbit_counts_stay_tiny_against_the_materialised_lattice():
    """The whole point: the raw lattice explodes and the orbit lattice does not."""
    counts, _ = oc.enumerate_orbits_by_augmentation(3, 7)
    materialised = list(
        af.enumerate_admissible_hyperplanes_by_flats(3, 7).flat_counts_by_rank)
    assert counts == [1, 1, 3, 10, 25, 36, 19]
    assert max(counts) * 500 < max(materialised)


# --------------------------------------------------------------------------
# the artifact, and the trap
# --------------------------------------------------------------------------

def test_artifact_records_closed_actions_and_the_orbit_ladder():
    art = _art()
    summary = art["summary"]
    assert summary["all_actions_closed"]
    assert summary["orbits"]["6"] == 8
    assert summary["orbits"]["7"] == 19
    assert summary["orbits"]["8"] == 56
    assert summary["hyperplanes"]["8"] == 166420
    for record in art["ranks"].values():
        assert record["images_outside_set"] == 0, record["d"]


def test_rank_9_is_reached_by_augmentation_and_agrees_where_checkable():
    """Rank 9 is the rank the materialising enumerator could not finish.

    It died at level 6 having checkpointed levels 1 to 5, so the orbit
    expansion is compared against those five stored counts and then continues
    past the point of death.
    """
    art = _art()
    aug = art["augmentation"]["9"]
    assert aug["orbits_by_rank"] == [1, 1, 3, 12, 45, 146, 364, 494, 231]
    assert aug["expanded_by_rank"][:6] == [1, 84, 3486, 78274, 968121, 6383706]
    assert aug["hyperplanes"] == 10004154
    assert aug["top_orbits"] == 231
    assert aug["peak_orbits"] == 494


def test_rank_10_has_a_committed_artifact_and_not_only_a_write_up():
    """The write-up quoted rank 10 with nothing behind it.

    The measuring script computed augmentation and then wrote the payload only
    inside the --ranks loop, so requesting augmentation alone produced no file
    at all. That is how 889,205,792 reached the documentation unbacked. This
    pins the record, its provenance, and the checksum over it.
    """
    art = _art()
    aug = art["augmentation"]["10"]
    assert aug["orbits_by_rank"] == [1, 1, 3, 12, 49, 189, 691, 1840, 2705, 1337]
    assert aug["expanded_by_rank"] == [
        1, 120, 7140, 241150, 4823910, 56904687,
        379791510, 1315196575, 1994358765, 889205792]
    assert aug["hyperplanes"] == 889205792
    assert aug["top_orbits"] == 1337
    assert aug["peak_orbits"] == 2705
    # the top level is not the peak: the peak is codimension two
    assert aug["peak_orbits"] > aug["top_orbits"]


def test_recorded_timings_say_whether_they_timed_a_full_run():
    """A resumed run reports almost no wall time, which times nothing.

    Shipping that as a measurement would be worse than shipping no number, so
    every augmentation record carries how many levels it restored.
    """
    art = _art()
    for key, record in art["augmentation"].items():
        assert "timing_is_a_full_run" in record, key
        assert record["levels_total"] == int(key) - 1
        restored = record["levels_restored_from_checkpoint"]
        assert 0 <= restored <= record["levels_total"]
        assert record["timing_is_a_full_run"] == (restored == 0)


def test_the_artifact_carries_provenance_and_a_checksum():
    import hashlib
    import json as _json

    art = _art()
    assert set(art["machine"]) >= {"python", "platform", "processor", "cpu_count"}
    digest = hashlib.sha256(
        _json.dumps(art["augmentation"], sort_keys=True,
                    separators=(",", ":")).encode("utf-8")).hexdigest()
    assert digest == art["augmentation_sha256"]


def test_the_artifact_no_longer_justifies_the_reduction_by_the_refuted_claim():
    """The description said the reduction holds because the ordered slice is a
    fundamental domain. That is the reason screening is NOT invariant."""
    art = _art()
    what = art["what"]
    assert "STABLE" in what
    assert "does NOT make screening an orbit invariant" in what


def test_augmentation_agrees_with_the_materialised_lattice_where_both_ran():
    art = _art()
    agreement = art["summary"]["augmentation_agrees_with_materialised"]
    for d in ("6", "7", "8"):
        assert agreement[d] is True, d


def test_rank_8_flat_lattice_is_replicated():
    """Three independent derivations now agree on this vector."""
    art = _art()
    assert art["ranks"]["8"]["flat_counts_by_rank"] == [
        1, 56, 1540, 21420, 147630, 467082, 565208, 166420]


def test_the_sign_convention_trap_is_real():
    """Permuting h without re-normalizing breaks closure, giving 15 not 8.

    This is pinned because the wrong answer looks entirely plausible: it is a
    valid orbit count of a valid group action on a set that is simply not the
    hyperplane set.
    """
    enum = af.enumerate_admissible_hyperplanes_by_flats(3, 6)
    stored = {(hp.h, hp.z) for hp in enum.hyperplanes}
    outside = 0
    for hp in enum.hyperplanes:
        image = list(hp.h)
        image[0], image[1] = image[1], image[0]
        if (tuple(image), hp.z) not in stored:
            outside += 1
    assert outside > 0, "if this is 0 the convention changed and the trap is gone"
    # every stored hyperplane has its first nonzero entry positive
    for hp in enum.hyperplanes:
        first = next((v for v in hp.h if v), 0)
        assert first > 0


# --------------------------------------------------------------------------
# order-free invariants: what prunes, and what was refuted
# --------------------------------------------------------------------------

def test_closed_form_survivor_count_matches_brute_force():
    """The Gaussian multinomial predicts the trace survivors, per orbit.

    This is what makes the row scan avoidable in principle, so it is checked
    against an actual count on every rank-7 orbit rather than sampled.
    """
    import collections

    from gpc_census.generation.exterior import exterior_weights
    from gpc_census.generation.reference_generator import oriented_taus
    from gpc_census.generation.ressayre import (
        _level_sets,
        _negative_root_set,
        admissible_candidate_from_tau,
    )

    n, d = 3, 7
    weights = exterior_weights(n, d)
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    by_orbit = collections.defaultdict(list)
    for tau in oriented_taus(enum):
        by_orbit[tuple(sorted(tau))].append(tau)

    checked = 0
    for members in by_orbit.values():
        first = admissible_candidate_from_tau(n, members[0])
        below = len(_level_sets(weights, first.h, first.z)[1])
        actual = sum(
            1 for tau in members
            if len(_negative_root_set(admissible_candidate_from_tau(n, tau).h))
            == below)
        assert oc.trace_survivor_count(first.h, below) == actual
        checked += 1
    assert checked == 35


def test_the_trace_dead_prune_never_kills_a_survivor():
    """B above C(d,2) kills a whole orbit on one comparison.

    The direction that matters is safety: a prune that removed a live orbit
    would silently lose facets, so every killed orbit is confirmed to have no
    survivor at all.
    """
    import collections

    from gpc_census.generation.exterior import exterior_weights
    from gpc_census.generation.reference_generator import oriented_taus
    from gpc_census.generation.ressayre import (
        _level_sets,
        _negative_root_set,
        admissible_candidate_from_tau,
    )

    n, d = 3, 7
    weights = exterior_weights(n, d)
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    by_orbit = collections.defaultdict(list)
    for tau in oriented_taus(enum):
        by_orbit[tuple(sorted(tau))].append(tau)

    killed_orbits = killed_rows = 0
    for members in by_orbit.values():
        first = admissible_candidate_from_tau(n, members[0])
        below = len(_level_sets(weights, first.h, first.z)[1])
        if not oc.orbit_is_trace_dead(first.h, below):
            continue
        killed_orbits += 1
        killed_rows += len(members)
        for tau in members:
            candidate = admissible_candidate_from_tau(n, tau)
            assert len(_negative_root_set(candidate.h)) != below
    assert killed_orbits == 10
    assert killed_rows == 1694


def test_generated_survivors_are_exactly_the_filtered_survivors():
    """Generation must reproduce the SET, not merely the count.

    A generator that produced the right number of arrangements while getting
    the arrangements themselves wrong would pass a count check and silently
    replace the candidate list with a different one, so every rank-7 orbit is
    compared set against set with the rows that actually survive filtering.
    """
    import collections

    from gpc_census.generation.exterior import exterior_weights
    from gpc_census.generation.reference_generator import oriented_taus
    from gpc_census.generation.ressayre import (
        _level_sets,
        _negative_root_set,
        admissible_candidate_from_tau,
    )

    n, d = 3, 7
    weights = exterior_weights(n, d)
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    by_orbit = collections.defaultdict(list)
    for tau in oriented_taus(enum):
        by_orbit[tuple(sorted(tau))].append(tau)

    total_survivors = total_rows = 0
    for members in by_orbit.values():
        first = admissible_candidate_from_tau(n, members[0])
        below = len(_level_sets(weights, first.h, first.z)[1])
        filtered = set()
        for tau in members:
            candidate = admissible_candidate_from_tau(n, tau)
            if len(_negative_root_set(candidate.h)) == below:
                filtered.add(tuple(candidate.h))
        generated = set(oc.generate_trace_survivors(first.h, below))
        assert generated == filtered
        total_survivors += len(filtered)
        total_rows += len(members)
    # the published partition: 10,682 oriented rows, 549 trace survivors
    assert total_rows == 10682
    assert total_survivors == 549


def test_survivor_generation_wastes_no_search_nodes():
    """Output linearity is the claim, so it is measured, not asserted.

    The prune interval is the exact set of attainable inversion numbers, so
    every surviving branch must complete. If that reasoning were wrong the
    generator would still be CORRECT and merely slow, which no other test here
    would notice.
    """
    nodes = [0]
    real_max = oc._max_inversions

    def counting_max(counts, total):
        nodes[0] += 1
        return real_max(counts, total)

    oc._max_inversions = counting_max
    try:
        produced = len(list(oc.generate_trace_survivors((3, 1, 1, 2, 2, 0, 4), 9)))
    finally:
        oc._max_inversions = real_max
    assert produced > 0
    # a node per placement step of a produced arrangement, up to the pruned
    # children examined at each step; linear in the output, not in 7! = 5040
    assert nodes[0] <= 12 * produced


def test_survivor_generation_agrees_with_the_closed_form_off_the_lattice():
    """Random multisets, so the agreement is not an artifact of real data."""
    import random

    random.seed(20260807)
    for _ in range(60):
        size = random.randint(2, 7)
        h = tuple(random.randint(0, 3) for _ in range(size))
        for below in range(size * (size - 1) // 2 + 2):
            generated = list(oc.generate_trace_survivors(h, below))
            assert len(generated) == oc.trace_survivor_count(h, below)
            for arrangement in generated:
                assert sorted(arrangement) == sorted(h)
                assert sum(1 for i in range(size)
                           for j in range(i + 1, size)
                           if arrangement[j] < arrangement[i]) == below


def test_the_inversion_generating_function_is_a_real_distribution():
    """Total mass must be the multiset permutation count, and it must be able
    to be lopsided, or the prune above would be vacuous."""
    from math import factorial

    poly = oc.inversion_generating_function((3, 1, 2))
    assert sum(poly) == factorial(3)
    poly = oc.inversion_generating_function((1, 1, 2, 2))
    assert sum(poly) == 6
    assert oc.trace_survivor_count((1, 1, 1), 1) == 0
