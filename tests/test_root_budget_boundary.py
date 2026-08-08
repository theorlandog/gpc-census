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


def test_verdict_records_both_halves():
    verdict = _artifact()["verdict"]
    assert verdict["positive_layer_bound"].startswith("HOLDS")
    assert verdict["zero_layer_uniform_bound"].startswith("REFUTED")
