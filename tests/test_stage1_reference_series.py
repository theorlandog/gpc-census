"""Guards for the blind rank-generator gates G2 and G3.

The claim is that the generator reproduces the published fixed-N=3 systems
WITHOUT being shown them. A test that merely re-reads the artifact cannot
establish that, so the rank-7 gate is re-run live here (about 8 seconds) and
compared against `gpc_census.constraints`. Rank 8 takes about 7 minutes and is
checked from the shipped artifact instead.
"""
import json
import pathlib

import pytest

from gpc_census.constraints import constraints
from gpc_census.generation import generate_fixed_n3_reference_rank
from gpc_census.generation.bdr import tau_from_constraint
from gpc_census.generation.model import IntegralConstraint

ROOT = pathlib.Path(__file__).resolve().parent / "data"
if not ROOT.exists():
    ROOT = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
ARTIFACT = ROOT / "stage1_reference_series.json"


def _art():
    return json.loads(ARTIFACT.read_text())


def _oriented(rows) -> set:
    return {tau_from_constraint(IntegralConstraint(tuple(r["coeffs"]), int(r["rhs"])), 3)
            for r in rows}


# --------------------------------------------------------------------------
# the live gate
# --------------------------------------------------------------------------

def test_rank_7_generates_blind_and_matches_the_published_table():
    result = generate_fixed_n3_reference_rank(7)
    assert result.proved_complete
    cs = result.candidate_system
    assert cs.status == "finalized_known_system_regression"
    reg = cs.known_system_regression
    assert reg is not None and reg.passed
    assert reg.generated_only == () or list(reg.generated_only) == []
    assert reg.published_only == () or list(reg.published_only) == []

    generated = {tau_from_constraint(r.constraint, 3) for r in cs.retained_rows}
    published = _oriented(constraints(3, 7)["inequalities"])
    assert generated == published
    assert len(cs.retained_rows) == 4


def test_enumeration_is_deterministic_and_hadamard_bound_is_respected():
    result = generate_fixed_n3_reference_rank(7)
    enum = result.enumeration
    assert len(enum.hyperplanes) == 5341
    assert len(result.taus) == 2 * len(enum.hyperplanes)
    # the modulus must beat the Hadamard bound, or the finite-field closure is
    # not a proof of the characteristic-zero relations
    assert enum.modulus ** 2 > enum.hadamard_bound_squared


def test_flat_lattice_counts_are_pinned():
    """Regression on the level representation.

    The level dict was changed from carrying a materialised RREF per flat to
    carrying only ``mask -> basis_indices``, which cut peak memory by about
    4.4x at rank 8. The whole flat lattice must come out identical, so its
    per-rank counts are pinned here and in the artifact.
    """
    from gpc_census.generation import enumerate_admissible_hyperplanes_by_flats
    enum = enumerate_admissible_hyperplanes_by_flats(3, 7)
    assert list(enum.flat_counts_by_rank) == [1, 35, 595, 4655, 15505, 19019, 5341]
    assert sum(enum.flat_counts_by_rank) == 45151
    # rank 6 flats are exactly the admissible hyperplanes
    assert enum.flat_counts_by_rank[-1] == len(enum.hyperplanes)


def test_rank_8_flat_lattice_matches_the_artifact():
    art = _art()
    rec = art["ranks"]["8"]
    assert rec["flat_counts_by_rank"] == [1, 56, 1540, 21420, 147630, 467082,
                                          565208, 166420]
    assert rec["flat_counts_by_rank"][-1] == rec["admissible_hyperplanes"]


# --------------------------------------------------------------------------
# the shipped ladder
# --------------------------------------------------------------------------

def test_artifact_records_gate_G2_passing_at_both_ranks():
    art = _art()
    for d in ("7", "8"):
        rec = art["ranks"].get(d)
        assert rec is not None, d
        assert rec["gate_passed"], d
        assert rec["proved_complete"], d
        assert rec["status"] == "finalized_known_system_regression", d
    assert {"7", "8"} <= set(art["gates_passed"])


@pytest.mark.parametrize("d", [7, 8])
def test_retained_rows_agree_with_the_published_system(d):
    art = _art()
    rec = art["ranks"][str(d)]
    published = _oriented(constraints(3, d)["inequalities"])
    assert _oriented(rec["retained_rows"]) == published
    assert rec["retained_count"] == len(constraints(3, d)["inequalities"])


def test_screening_partition_is_exhaustive():
    """Every candidate is resolved: no row may be left in limbo."""
    art = _art()
    for d, rec in art["ranks"].items():
        if not rec.get("gate_passed"):
            continue
        c = rec["screening_counts"]
        assert c["evaluation_unresolved_rows"] == 0, d
        assert c["symbolic_fallback_exception_unresolved_rows"] == 0, d
        assert c["certified_rows"] + c["exact_rejected_rows"] + c["structural_rows"] \
            == c["input_tau_rows"], d
        assert c["normalized_tau_rows"] == rec["oriented_taus"], d


def test_cost_scaling_is_recorded():
    """The enumeration cost is the thing that gates the frontier."""
    art = _art()
    r7, r8 = art["ranks"]["7"], art["ranks"]["8"]
    assert r7["admissible_hyperplanes"] == 5341
    assert r8["admissible_hyperplanes"] == 166420
    # the growth is the reason G3 is not free
    assert r8["admissible_hyperplanes"] > 25 * r7["admissible_hyperplanes"]


def test_operational_failures_are_not_counted_as_gates():
    art = _art()
    for d, rec in art["ranks"].items():
        if rec.get("status") == "operational_failure":
            assert not rec["gate_passed"], d
            assert d not in art["gates_passed"]


def test_checkpoint_round_trip_is_identical():
    """A resumed enumeration must equal a fresh one, exactly.

    Checkpointing exists so a rank whose cost is hours can be finished across
    several runs. That is only safe if resuming cannot perturb the lattice, so
    the flat counts AND the hyperplane tuple are compared, not just the count.
    """
    import tempfile

    from gpc_census.generation import admissible_flats as af

    enum = af.enumerate_admissible_hyperplanes_by_flats
    enum.cache_clear()
    fresh = enum(3, 7)
    with tempfile.TemporaryDirectory() as tmp:
        enum.cache_clear()
        written = enum(3, 7, tmp)
        enum.cache_clear()
        resumed = enum(3, 7, tmp)
    enum.cache_clear()

    assert list(fresh.flat_counts_by_rank) == list(written.flat_counts_by_rank)
    assert list(fresh.flat_counts_by_rank) == list(resumed.flat_counts_by_rank)
    assert fresh.hyperplanes == written.hyperplanes == resumed.hyperplanes
    assert fresh.modulus == resumed.modulus


def test_checkpoint_files_are_plain_text_not_pickle():
    """A checkpoint is data; loading one must not be able to execute anything."""
    import tempfile

    from gpc_census.generation import admissible_flats as af

    enum = af.enumerate_admissible_hyperplanes_by_flats
    with tempfile.TemporaryDirectory() as tmp:
        enum.cache_clear()
        enum(3, 7, tmp)
        enum.cache_clear()
        files = sorted(pathlib.Path(tmp).glob("flats_3_7_rank*.txt"))
        assert len(files) == 6            # ranks 1..d-1
        head = files[0].read_text().splitlines()[0]
        mask_hex, sep, basis = head.partition(":")
        assert sep == ":"
        int(mask_hex, 16)                 # parses as hex
        assert all(part.isdigit() for part in basis.split(",") if part)
        # no leftover partial files from an interrupted write
        assert not list(pathlib.Path(tmp).glob("*.partial"))


def test_vectorised_extension_matches_the_scalar_definition_past_64_bits():
    """Regression: flat masks exceed 64 bits from rank 9 on.

    A mask carries one bit per exterior weight, and rank 9 has 84 of them, so
    ``mask >> numpy_array`` overflows a C long there while every rank the gates
    currently cover (56 weights at rank 8) stays inside it. This compares the
    vectorised extension against the scalar definition at rank 9, where the
    masks are genuinely 84 bits wide.
    """
    import numpy as np

    from gpc_census.generation import admissible_flats as af
    from gpc_census.generation.exterior import exterior_weights

    n, d = 3, 9
    modulus, _ = af.rank_preserving_modulus(n, d)
    weights = exterior_weights(n, d)
    points = tuple(tuple(w) + (1,) for w in weights)
    assert len(points) == 84                       # > 64, the whole point

    inverse_table = af._inverse_table(modulus)
    points_array = np.asarray(points, dtype=np.int64)
    states: dict[int, tuple[int, ...]] = {0: ()}
    for _ in range(2):
        vectorised = af._extend_one_rank(
            states, points, modulus, points_array, inverse_table)
        scalar = af._extend_one_rank_scalar(states, points, modulus)
        assert vectorised == scalar
        assert max(vectorised).bit_length() > 64
        states = vectorised
    assert len(states) == 3486


def test_inverse_table_covers_every_modulus_the_generator_selects():
    from gpc_census.generation import admissible_flats as af

    for d in range(6, 16):
        modulus, bound_squared = af.rank_preserving_modulus(3, d)
        assert modulus * modulus > bound_squared
        assert af._inverse_table(modulus) is not None, d
    table = af._inverse_table(257)
    for value in (1, 2, 3, 256):
        assert (int(table[value]) * value) % 257 == 1
