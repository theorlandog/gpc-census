import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "census_support_geometry.json"
CENSUS = DATA / "rdm_observability_census.json"

MODULE_PATH = ROOT / "scripts" / "census_support_geometry.py"
SPEC = importlib.util.spec_from_file_location("census_support_geometry", MODULE_PATH)
csg = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(csg)


def record(dets, system="(3,6)"):
    return {
        "system": system,
        "index": 0,
        "classified": "DESIGN-INT",
        "closed_form": {"support_dets": dets},
    }


def test_disjoint_three_fermion_pair_is_not_order2_connected():
    row = csg.analyze_record(record([[0, 1, 2], [3, 4, 5]]))
    assert row["pair_incidence_full"]
    assert row["non_toric_phase_complexity"] == 0
    assert not row["orders"]["2"]["connected"]
    assert row["minimum_collision_free_reconstruction_order"] == 3


def test_single_support():
    row = csg.analyze_record(record([[0, 1, 2]]))
    assert row["non_toric_phase_complexity"] == 0
    assert row["minimum_collision_free_reconstruction_order"] == 2


def test_minimum_order_scan_starts_at_two():
    """The reported minimum is a minimum over q >= 2, by construction.

    Two determinants differing by one excitation are already connected by a
    one-body channel, but the campaign scans from q = 2 because order two is
    the question it exists to answer. Read the field with that bound in mind.
    """
    row = csg.analyze_record(record([[0, 1, 2], [0, 1, 3]]))
    assert row["minimum_collision_free_reconstruction_order"] == 2
    assert "1" not in row["orders"]


def test_kappa_is_singleton_incidence_nullity():
    row = csg.analyze_record(record([[0, 1, 2], [0, 3, 4], [1, 3, 5], [2, 4, 5]]))
    assert row["non_toric_phase_complexity"] == row["carrier_size"] - row["singleton_rank"]
    assert row["support_polytope_dimension"] == row["singleton_rank"] - 1


def test_excitation_histogram_counts_every_pair():
    row = csg.analyze_record(record([[0, 1, 2], [0, 1, 3], [3, 4, 5]]))
    assert sum(row["excitation_order_histogram"].values()) == 3


def test_campaign_reproduces_the_shipped_census_counts():
    """The campaign imports no project generator, so agreement is real.

    It must land on the same 771/28 split and the same per-class totals as
    results/data/rdm_observability_census.json.
    """
    summary = json.loads(ARTIFACT.read_text())["summary"]
    census = json.loads(CENSUS.read_text())["summary"]
    assert summary["records"] == census["records"] == 799
    assert summary["global"]["order2_certified"] == census["gamma2_reconstructible"] == 771
    assert summary["global"]["order2_unresolved"] == 28
    assert summary["global"]["kappa_positive"] == census["positive_singleton_kernel"] == 63
    for name, block in census["by_class"].items():
        got = summary["by_classification"][name]
        assert got["records"] == block["records"]
        assert got["order2_certified"] == block["gamma2_reconstructible"]


def test_exact_implications_hold_and_are_scoped_to_this_census():
    implications = json.loads(ARTIFACT.read_text())["summary"]["exact_implications"]
    assert all(implications.values())
    assert implications["positive_kappa_implies_interference_in_this_census"]
    assert implications["positive_kappa_implies_order2_certified_in_this_census"]
    assert implications["order2_unresolved_all_kappa_zero"]


def test_all_order2_exceptions_are_design_int_simplices():
    summary = json.loads(ARTIFACT.read_text())["summary"]
    unresolved = summary["order2_unresolved"]
    assert len(unresolved) == 28
    assert all(r["classified"] == "DESIGN-INT" for r in unresolved)
    assert all(r["non_toric_phase_complexity"] == 0 for r in unresolved)
    assert all(r["support_is_simplex"] for r in unresolved)


def test_positive_kappa_records_are_all_interference_and_order_two():
    positive = json.loads(ARTIFACT.read_text())["summary"]["positive_kappa_records"]
    assert len(positive) == 63
    assert all(r["classified"] == "INTERFERENCE" for r in positive)
    assert all(r["minimum_collision_free_reconstruction_order"] == 2 for r in positive)
