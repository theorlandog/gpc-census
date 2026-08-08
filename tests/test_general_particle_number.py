"""Guards for the general-``N`` finalizer, and for the gate that still fails.

The finalization layer was the fixed-``N=3`` choke point: certificate replay,
constraint conversion, the known-system lookup, the cyclic cross-check, the
serialized output and the Farkas trace multiplier all hardcoded three. This
pins what the parameterization achieved AND what it did not, because the
(4,7) particle-hole gate does not pass yet and a test suite that quietly
omitted that would misrepresent the state of the work.
"""
import json
import pathlib

import pytest

from gpc_census.constraints import constraints
from gpc_census.generation.bdr import tau_from_constraint
from gpc_census.generation.candidate_pipeline import screen_tau_candidates
from gpc_census.generation.model import IntegralConstraint
from gpc_census.generation.orbit_native import enumerate_orbit_survivor_taus
from gpc_census.generation.redundancy import (
    certify_ordered_slice_redundancy,
    ordered_pauli_inequalities,
    verify_ordered_slice_reduction_certificate,
)

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"


def _published(n, d):
    table = constraints(n, d)
    return tuple(IntegralConstraint(tuple(r["coeffs"]), int(r["rhs"]))
                 for r in table["inequalities"])


# --------------------------------------------------------------------------
# what the parameterization achieved
# --------------------------------------------------------------------------

def test_the_n3_manifest_is_still_byte_identical():
    """The acceptance condition: no stored N=3 payload may move."""
    from gpc_census.generation.reference_generator import (
        generate_fixed_n3_reference_rank,
    )

    stored = json.loads((ROOT / "stage1_ressayre_3_7.json").read_text())
    assert generate_fixed_n3_reference_rank(7).as_dict() == stored


@pytest.mark.parametrize("n", [3, 4])
def test_redundancy_reduces_the_published_table_at_its_own_particle_number(n):
    """The Farkas identity was already parameterized; only guards blocked it.

    Reducing a published system must return every one of its rows, since a
    published system is already irredundant. At N=4 this exercises the trace
    multiplier against sum(lambda) = 4 rather than 3.
    """
    rows = _published(n, 7)
    certificate = certify_ordered_slice_redundancy(rows, 7, particle_number=n)
    assert len(certificate.retained_indices) == len(rows)
    assert certificate.particle_number == n
    assert verify_ordered_slice_reduction_certificate(certificate)


def test_the_ordered_pauli_rows_do_not_mention_the_particle_number():
    """They are lambda_0 <= 1, dominance, lambda_last >= 0, for any N.

    The old rank window 7 <= d <= 13 was an N=3 habit, not a limit.
    """
    for d in (4, 6, 7, 15):
        rows = ordered_pauli_inequalities(d)
        assert len(rows) == d + 1
    with pytest.raises(ValueError):
        ordered_pauli_inequalities(1)


def test_the_particle_number_is_carried_not_assumed():
    from gpc_census.generation.candidate_system import compose_candidate_system

    enumeration = enumerate_orbit_survivor_taus(4, 7)
    screening = screen_tau_candidates(
        enumeration.survivor_taus, particle_number=4, symbolic_fallback=True)
    result = compose_candidate_system(screening)
    assert result.particle_number == 4
    assert result.as_dict()["particle_number"] == 4


def test_enumeration_and_screening_are_correct_at_n4():
    """Every published (4,7) row is produced and certified.

    This is the half of the gate that DOES pass, and it is what localises the
    failure below: candidate generation and Ressayre certification are right at
    N=4, so nothing upstream of reduction is losing the rows.
    """
    published = {tau_from_constraint(row, 4) for row in _published(4, 7)}
    enumeration = enumerate_orbit_survivor_taus(4, 7)
    assert published <= set(enumeration.survivor_taus)
    screening = screen_tau_candidates(
        enumeration.survivor_taus, particle_number=4, symbolic_fallback=True)
    certified = {tuple(row.tau) for row in screening.candidates
                 if row.status == "certified"}
    assert published <= certified
    assert screening.all_candidates_resolved


def test_particle_hole_duality_holds_at_the_enumeration_level():
    """(3,7) and (4,7) must enumerate identically, and they do.

    Complementation is a bijection between 3-subsets and 4-subsets of seven
    modes, so the flat lattices agree. Matching counts here are a confirmation
    rather than a symptom of N failing to propagate, which is worth pinning
    because they LOOK like a symptom.
    """
    three = enumerate_orbit_survivor_taus(3, 7)
    four = enumerate_orbit_survivor_taus(4, 7)
    assert three.hyperplane_count == four.hyperplane_count == 5341
    assert len(three.orbits) == len(four.orbits) == 19
    assert (len(three.survivor_taus) == len(four.survivor_taus)
            == three.predicted_survivor_count == 549)


# --------------------------------------------------------------------------
# what still fails: the GN-0 gate
# --------------------------------------------------------------------------

@pytest.mark.xfail(reason="GN-0 not met: reduction drops 3 of the 4 published "
                          "(4,7) rows although all four are certified",
                   strict=True)
def test_gn0_particle_hole_gate_at_4_7():
    """The gate, recorded as failing rather than omitted.

    Every published (4,7) row survives enumeration and screening, and the
    published table reduces to itself, yet composing the full 49-row certified
    set retains one row instead of four. The defect is therefore inside
    reduction against the larger candidate set, not in candidate generation.
    Marked strict so that fixing it fails this test and forces the update.
    """
    from gpc_census.generation.orbit_native import generate_orbit_native_rank

    result = generate_orbit_native_rank(7, particle_number=4)
    regression = result.as_dict()["known_system_regression"]
    assert regression["passed"], regression
