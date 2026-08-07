"""Acceptance gates for the orbit-native constructor.

The fast components (orbit enumeration, direct survivor generation, the Hall
structural-zero filter) were each verified in isolation and none of them was on
the path that actually builds a system. This module pins the integration, and
the gate is deliberately not "does it run": it is that the orbit-native path
REPRODUCES the materialising reference exactly, because a constructor that
skips 98 percent of the rows is only trustworthy if it lands on the same system.

Three independent things are checked, because any one of them alone could pass
while the constructor silently lost candidates:

- the emitted survivor taus equal, as a SET, the rows the reference path keeps
  after its trace test;
- the retained constraints, the reduction certificate and the known-system
  regression are identical;
- count conservation closes, with the discarded rows accounted for in closed
  form by the Gaussian multinomial rather than by having been built.
"""
import json
import pathlib

import pytest

from gpc_census.generation import admissible_flats as af
from gpc_census.generation import orbit_native as on
from gpc_census.generation.exterior import exterior_weights
from gpc_census.generation.reference_generator import (
    generate_fixed_n3_reference_rank,
    oriented_taus,
)
from gpc_census.generation.ressayre import (
    _level_sets,
    _negative_root_set,
    admissible_candidate_from_tau,
)

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"


def _reference_survivor_taus(n, d):
    """The taus the materialising path keeps, computed the slow honest way."""
    weights = exterior_weights(n, d)
    enum = af.enumerate_admissible_hyperplanes_by_flats(n, d)
    kept = set()
    for tau in oriented_taus(enum):
        candidate = admissible_candidate_from_tau(n, tau)
        below = len(_level_sets(weights, candidate.h, candidate.z)[1])
        if len(_negative_root_set(candidate.h)) == below:
            kept.add(tau)
    return kept


def _strip_index(rows):
    """Drop the row's position in the screened list.

    It is bookkeeping, not mathematics, and it cannot agree: the reference path
    screens 10,682 rows at rank 7 where the orbit-native path screens 549.
    """
    return [{k: v for k, v in row.items() if k != "screening_index"}
            for row in rows]


# --------------------------------------------------------------------------
# the survivors themselves
# --------------------------------------------------------------------------

@pytest.mark.parametrize(("d", "hyperplanes", "rows", "survivors"), [
    (6, 362, 724, 64),
    (7, 5341, 10682, 549),
])
def test_generated_taus_are_exactly_the_reference_survivors(
        d, hyperplanes, rows, survivors):
    """Set equality, not count equality.

    A generator emitting the right NUMBER of wrong taus would produce a
    plausible system from the wrong candidates, and every count in the manifest
    would still look right.
    """
    enumeration = on.enumerate_orbit_survivor_taus(3, d)
    assert enumeration.hyperplane_count == hyperplanes
    assert enumeration.oriented_row_count == rows
    assert len(enumeration.survivor_taus) == survivors
    assert set(enumeration.survivor_taus) == _reference_survivor_taus(3, d)


@pytest.mark.parametrize(("d", "self_paired"), [(6, 1), (7, 3)])
def test_self_paired_orbits_are_detected_and_not_double_counted(d, self_paired):
    """The trap: an orbit that already contains its own reverse orientation.

    Processing both orientations of such an orbit emits every row twice. The
    counts are known independently: the oriented orbit totals are 15 and 35
    against unoriented 8 and 19, so exactly 1 and 3 orbits are self paired.
    """
    enumeration = on.enumerate_orbit_survivor_taus(3, d)
    assert sum(1 for r in enumeration.orbits if r.self_paired) == self_paired
    assert len(set(enumeration.survivor_taus)) == len(enumeration.survivor_taus)


def test_count_conservation_closes():
    """Coverage without materialising the rows that are thrown away.

    The left side is orbit-stabilizer plus the Gaussian multinomial, the right
    side is what the run did. A lost orbit or a dropped survivor cannot balance
    it.
    """
    enumeration = on.enumerate_orbit_survivor_taus(3, 7)
    assert enumeration.counts_are_conserved
    assert (enumeration.trace_rejected_count + len(enumeration.survivor_taus)
            == enumeration.oriented_row_count)
    assert enumeration.trace_rejected_count == 10133
    per_orbit = sum(2 * record.orbit_size for record in enumeration.orbits)
    assert per_orbit == enumeration.oriented_row_count


def test_orbit_representatives_carry_their_exact_hyperplane():
    """Each representative's exact (h, z) must reproduce its own flat.

    This is the join between the finite-field closure that found the flat and
    the exact rational normal the screening stage consumes; if they disagreed,
    the orbit would be describing a different hyperplane.
    """
    enumeration = on.enumerate_orbit_survivor_taus(3, 7)
    for record in enumeration.orbits:
        assert on._incidence_mask(3, 7, record.h, record.z) == record.canonical_mask
        assert record.automorphism_order * record.orbit_size == on._factorial(7)


# --------------------------------------------------------------------------
# the system, against the materialising reference
# --------------------------------------------------------------------------

def test_orbit_native_reproduces_the_reference_system_at_rank_7():
    """The gate. Same system, built from 549 rows instead of 10,682."""
    native = on.generate_orbit_native_rank(7)
    reference = generate_fixed_n3_reference_rank(7)
    assert native.proved_complete
    a, b = native.as_dict(), reference.as_dict()
    assert _strip_index(a["retained_rows"]) == _strip_index(b["retained_rows"])
    assert (a["ordered_slice_reduction_certificate"]
            == b["ordered_slice_reduction_certificate"])
    assert a["known_system_regression"] == b["known_system_regression"]
    certified = {tuple(r.tau) for r in native.screening.candidates
                 if r.status == "certified"}
    assert certified == {tuple(r.tau) for r in reference.screening.candidates
                         if r.status == "certified"}


def test_the_published_rank_7_partition_is_reproduced():
    """1 structural, 49 certified, 499 determinant zeros, 0 unresolved.

    Those are the numbers docs/fixed_n3_generator_methodology.md publishes for
    the materialising run, and they must survive skipping the 10,133 rejected
    rows entirely.
    """
    native = on.generate_orbit_native_rank(7)
    counts = native.screening.counts.as_dict()
    assert counts["structural_rows"] == 1
    assert counts["certified_rows"] == 49
    assert counts["symbolic_zero_rejected_rows"] == 499
    assert counts["evaluation_unresolved_rows"] == 0
    assert counts["trace_rejected_rows"] == 0, "rejected rows must not be built"


# --------------------------------------------------------------------------
# the Hall filter, now on the production path
# --------------------------------------------------------------------------

def test_hall_is_actually_reached_by_the_production_fallback():
    """A guard against the integration silently regressing.

    The filter is a pure accelerator with a deliberately identical outcome, so
    nothing in the manifest would change if it stopped being called. This
    asserts it is on the path by counting calls.
    """
    from gpc_census.generation import candidate_pipeline as cp

    calls = [0]
    real = cp.tangent_determinant_is_structurally_zero

    def counting(*args, **kwargs):
        calls[0] += 1
        return real(*args, **kwargs)

    cp.tangent_determinant_is_structurally_zero = counting
    try:
        native = on.generate_orbit_native_rank(7)
    finally:
        cp.tangent_determinant_is_structurally_zero = real
    # 499 rows reach the symbolic fallback at rank 7, and every one of them is
    # offered to the matching test before any determinant algebra runs
    assert calls[0] == 499
    assert native.screening.counts.as_dict()["symbolic_zero_rejected_rows"] == 499


def test_structural_and_algebraic_zeros_split_as_measured():
    """474 structural, 25 algebraic at rank 7, per screening_orbit_structure.

    The split is what justifies putting matching first, so it is checked
    against the production candidate set rather than a separate enumeration.
    """
    from gpc_census.generation.ressayre import (
        tangent_determinant_is_structurally_zero,
        verify_hall_violator,
    )

    native = on.generate_orbit_native_rank(7)
    structural = algebraic = 0
    for row in native.screening.candidates:
        if row.status != "symbolic_zero_rejected":
            continue
        candidate = admissible_candidate_from_tau(3, row.tau)
        is_structural, violator = tangent_determinant_is_structurally_zero(
            3, 7, candidate)
        if is_structural:
            assert verify_hall_violator(3, 7, candidate, violator)
            structural += 1
        else:
            algebraic += 1
    assert structural == 474
    assert algebraic == 25


# --------------------------------------------------------------------------
# the stale claims that were refuted
# --------------------------------------------------------------------------

@pytest.mark.parametrize("path", [
    "src/gpc_census/generation/orbit_canonical.py",
    "scripts/orbit_canonical_flats.py",
])
def test_no_module_claims_orbit_members_are_facets_together(path):
    """Screening is not an orbit invariant and the docstrings said it was."""
    source = pathlib.Path(__file__).resolve().parents[1] / path
    if not source.exists():
        pytest.skip("source not present")
    text = source.read_text()
    # the phrase may appear, but only inside the correction that denies it
    assert "so two\nadmissible hyperplanes in the same ``S_d`` orbit are facets" not in text
    assert "carry the same information: one is a facet" not in text
    if "facets together or not at all" in text:
        assert "does NOT say" in text or "It does NOT" in text
    assert "screening_orbit_structure" in text
    assert "STABLE" in text, "the reduction must stand on S_d stability instead"


def test_the_orbit_native_manifest_is_serializable_and_states_its_method():
    payload = on.enumerate_orbit_survivor_taus(3, 6).as_dict()
    json.dumps(payload)
    assert payload["method"] == "orbit_native_canonical_augmentation"
    assert payload["counts_are_conserved"]
    assert payload["materialised_rows_avoided"] == 724 - 64
