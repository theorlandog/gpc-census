"""Guards for the root-budget boundary layer and its census validation."""

from __future__ import annotations

import importlib.util
import json
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

ROOT_BUDGET_BOUNDARY = DATA / "root_budget_boundary.json"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = _load("root_budget_boundary")


def _artifact() -> dict:
    return json.loads(ROOT_BUDGET_BOUNDARY.read_text())


def _layer(particle_number: int, which: str) -> dict:
    for entry in _artifact()["boundary_layers"]:
        if entry["particle_number"] == particle_number:
            return entry[which]
    raise AssertionError(f"no layer for N={particle_number}")


def test_root_budget_is_binom_d_2():
    artifact = _artifact()
    assert artifact["root_budget"] == math.comb(artifact["orbitals"], 2) == 780


def test_twenty_forty_positive_layer_counts():
    """The headline reduction: 1.4e11 weights down to 15,359."""
    layer = _layer(20, "positive_layer")
    assert layer["ambient_weight_count"] == 137846528820
    assert layer["accepted_weight_types"] == 15359
    assert layer["maximum_young_area"] == 42
    assert _layer(20, "zero_or_positive_layer")["accepted_weight_types"] == 15371


def test_layer_counts_across_particle_numbers():
    expected = {3: 2587, 4: 4634, 5: 6525, 10: 12372, 20: 15359}
    for particle_number, count in expected.items():
        assert _layer(particle_number, "positive_layer")["accepted_weight_types"] == count


def test_truncated_layer_recomputes():
    """A small case is recomputed rather than read."""
    layer = MODULE.truncated_layer(3, 40, 780)
    assert layer["accepted_weight_types"] == 2587
    assert layer["maximum_ideal_size"] <= 780
    assert layer["minimum_frontier_ideal_size"] > 780


def test_positive_half_holds_on_every_mechanism():
    """No regularity needed, and verified on all 280 census mechanisms."""
    total = 0
    for record in _artifact()["census_validation"].values():
        mechanisms = record["mechanisms"]
        assert record["positive_family_is_a_dominance_ideal"] == mechanisms
        assert record["positive_family_within_root_budget"] == mechanisms
        assert record["positive_weight_bound_holds"] == mechanisms
        assert record["worst_positive_weight_ideal"] <= record["root_budget"]
        total += mechanisms
    assert total == 280


def test_positive_family_budget_is_attained():
    """The bound is tight: |Omega_+| reaches binom(9,2) at rank 9."""
    rank9 = _artifact()["census_validation"]["(3,9)"]
    assert rank9["largest_positive_family"] == rank9["root_budget"] == 36


def test_zero_half_is_refuted_off_the_regular_stratum():
    """The uniform zero bound fails, and the failures are recorded."""
    artifact = _artifact()["census_validation"]
    failures = {}
    for system, record in artifact.items():
        failures[system] = record["mechanisms"] - record["zero_weight_uniform_bound_holds"]
        assert record["worst_zero_weight_ideal"] > record["root_budget"] + 1
    assert failures == {"(3,7)": 4, "(3,8)": 10, "(3,9)": 8}
    assert sum(failures.values()) == 22


def test_regularity_is_nearly_absent_in_the_census():
    artifact = _artifact()["census_validation"]
    assert artifact["(3,7)"]["regular_mechanisms"] == 0
    assert artifact["(3,8)"]["regular_mechanisms"] == 0
    assert artifact["(3,9)"]["regular_mechanisms"] == 3
    assert sum(r["regular_mechanisms"] for r in artifact.values()) == 3


def test_corrected_bound_holds_everywhere():
    for record in _artifact()["census_validation"].values():
        assert record["zero_weight_corrected_bound_holds"] == record["mechanisms"]


def test_violations_carry_a_witness():
    """A refutation must name rows, not just a count."""
    for record in _artifact()["census_validation"].values():
        violations = record["uniform_bound_violations"]
        assert violations
        for violation in violations:
            assert violation["largest_zero_weight_ideal"] > violation["uniform_bound"]
            assert (
                violation["largest_zero_weight_ideal"] <= violation["corrected_bound"]
            )
            assert violation["distinct_levels"] < len(violation["tau"])


def test_per_levi_budget_is_the_repository_max_inversions():
    """The proposed refinement is a function the repository already applies."""
    identity = _artifact()["levi_budget_identity"]
    assert identity["identical_to_max_inversions"]
    assert identity["mismatches"] == 0
    assert identity["trials"] == 400
    assert MODULE.levi_budget([2, 2, 1], 5) == (25 - 9) // 2


def test_per_levi_budget_sharpens_positive_and_aggravates_zero():
    """A smaller budget is a stronger claim on both sides."""
    artifact = _artifact()["census_validation"]
    uniform_zero = 0
    levi_zero = 0
    for record in artifact.values():
        mechanisms = record["mechanisms"]
        assert record["levi_positive_family_within_budget"] == mechanisms
        assert record["levi_positive_weight_bound_holds"] == mechanisms
        assert record["levi_zero_weight_bound_holds"] < mechanisms
        assert (
            record["levi_zero_weight_bound_holds"]
            <= record["zero_weight_uniform_bound_holds"]
        )
        uniform_zero += mechanisms - record["zero_weight_uniform_bound_holds"]
        levi_zero += mechanisms - record["levi_zero_weight_bound_holds"]
    assert uniform_zero == 22
    assert levi_zero == 70


def test_durfee_bound_holds_with_slack():
    for record in _artifact()["census_validation"].values():
        measured = record["measured_max_durfee_rank_of_positive_weights"]
        assert measured == 2
        assert measured <= record["durfee_rank_bound_from_global_budget"] == 3


def test_half_filling_transpose_orbit_counts():
    transpose = _artifact()["half_filling_transpose"]
    assert transpose["(10,20)"]["transpose_orbits"] == 993
    assert transpose["(15,30)"]["transpose_orbits"] == 4249
    assert transpose["(20,40)"]["transpose_orbits"] == 11478
    for record in transpose.values():
        assert record["applies"]
        assert (
            record["sign_control_universe"] + record["self_conjugate"]
            == 2 * record["transpose_orbits"]
        )


def test_transpose_is_an_involution_of_the_rectangle():
    assert MODULE.transpose((3, 1, 0), 3) == (2, 1, 1)
    assert MODULE.transpose(MODULE.transpose((3, 1, 0), 3), 3) == (3, 1, 0)
    assert MODULE.transpose_orbit_counts(3, 10)["applies"] is False


def test_weighted_profile_laws_hold_everywhere():
    """The repair: what the per-weight zero bound could not do, profiles do."""
    total = 0
    for record in _artifact()["census_validation"].values():
        mechanisms = record["mechanisms"]
        assert record["profile_positive_law_holds"] == mechanisms
        assert record["profile_zero_law_holds"] == mechanisms
        assert record["profile_mass_matches_weight_count"] == mechanisms
        total += mechanisms
    assert total == 280


def test_profile_laws_hold_where_the_per_weight_bound_failed():
    """Every per-weight failure is covered by the profile formulation."""
    for record in _artifact()["census_validation"].values():
        assert record["levi_zero_weight_bound_holds"] < record["mechanisms"]
        assert record["profile_zero_law_holds"] == record["mechanisms"]


def test_profile_laws_recompute_on_a_degenerate_surviving_normal():
    """Recomputed, not read: a real (3,8) survivor with repeated levels."""
    laws = MODULE.profile_budget_laws((2, 0, 0, 0, -1, -1, -1, -2), 3)
    assert laws["positive_law_holds"]
    assert laws["zero_law_holds"]
    assert laws["positive_family_mass"] <= laws["levi_budget"]


def test_profile_laws_are_conditional_and_detect_trace_dead_normals():
    """The laws assume the root budget; violating them certifies trace-death.

    tau = (2,2,-1,-1,-1,-1,0,0) has |Omega_+| = 24 against R_L = 20, so no
    arrangement of it can pass the Ressayre trace test. The laws report failure,
    which is the correct verdict and doubles as a pruning test.
    """
    laws = MODULE.profile_budget_laws((2, 2, 0, 0, -1, -1, -1, -1), 3)
    assert laws["levi_budget"] == 20
    assert laws["positive_family_mass"] == 24 > laws["levi_budget"]
    assert not laws["positive_law_holds"]


def test_dominance_and_profiles():
    assert MODULE.dominates((2, 0, 0), (1, 1, 0))
    assert not MODULE.dominates((0, 1, 1), (1, 1, 0))
    profiles = MODULE.occupancy_profiles([2, 2], 2)
    assert set(profiles) == {(0, 2), (1, 1), (2, 0)}


def test_gate_ladder_table_reproduces():
    """Every entry number of the proposed validation ladder, recomputed."""
    expected = {
        (13, 17): (2380, 609, 3.9, 4),
        (10, 20): (184756, 1940, 95.2, 4),
        (15, 30): (155117520, 8409, 18446.6, 5),
        (20, 40): (137846528820, 22817, 6041395.8, 5),
    }
    ladder = {tuple(r["system"]): r for r in _artifact()["gate_ladder"]}
    assert set(ladder) == set(expected)
    for system, (ambient, sign, compression, durfee) in expected.items():
        record = ladder[system]
        assert record["ambient_weight_count"] == ambient
        assert record["sign_control_universe"] == sign
        assert record["compression_factor"] == compression
        assert record["max_durfee_rank"] == durfee
        assert (
            record["zero_or_positive_layer"] + record["minimal_excluded_frontier"]
            == sign
        )


def test_gate_ladder_particle_hole_normalization():
    """(13,17) is computed as (4,17); the half-filled rungs are self-dual."""
    ladder = {tuple(r["system"]): r for r in _artifact()["gate_ladder"]}
    assert ladder[(13, 17)]["effective_system"] == [4, 17]
    assert ladder[(13, 17)]["particle_hole_reduction_used"]
    assert not ladder[(13, 17)]["transpose"]["applies"]
    for system in ((10, 20), (15, 30), (20, 40)):
        record = ladder[system]
        assert record["effective_system"] == list(system)
        assert not record["particle_hole_reduction_used"]
        assert record["transpose"]["applies"]


def test_excitation_bound_holds_but_is_vacuous_on_the_census():
    """Proved and unviolated, yet it can bind at only 4 of 280 mechanisms."""
    artifact = _artifact()["census_validation"]
    can_bind = 0
    total = 0
    for record in artifact.values():
        assert record["excitation_bound_violations"] == 0
        assert record["max_positive_excitation_defect"] <= record[
            "structural_defect_ceiling"
        ]
        can_bind += record["excitation_bound_can_bind"]
        total += record["mechanisms"]
    assert total == 280
    assert can_bind == 4
    assert artifact["(3,9)"]["excitation_bound_can_bind"] == 0


def test_defect_is_bounded_by_particle_number_structurally():
    """delta <= N always, which is why the log2 bound cannot bind at small N."""
    assert MODULE.excitation_defect((1, 1, 1), [4, 4, 4]) == 3
    assert MODULE.excitation_defect((0, 3, 0), [2, 3, 2]) == 0
    assert MODULE.excitation_defect((2, 1), [2, 4]) == 1


def test_excitation_bound_first_bites_at_fourteen_orbitals():
    rows = _artifact()["excitation_vacuity_threshold"]
    binding = [r for r in rows if r["bound_can_bind"]]
    assert binding[0]["orbitals"] == 14
    assert binding[0]["half_filling_particles"] == 7
    for row in rows:
        if row["orbitals"] <= 12:
            assert not row["bound_can_bind"]


def test_verdict_records_both_halves():
    verdict = _artifact()["verdict"]
    assert verdict["positive_layer_bound"].startswith("HOLDS")
    assert verdict["zero_layer_uniform_bound"].startswith("REFUTED")
    assert verdict["weighted_occupancy_profile_laws"].startswith("HOLD")
    assert "VACUOUS" in verdict["levi_excitation_bound"]
