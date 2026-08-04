import importlib.util
import json
import math
import pathlib
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
MANIFEST = DATA / "stage1_ressayre_3_7.json"

MODULE_PATH = ROOT / "scripts" / "verify_stage1_ressayre_3_7.py"
SPEC = importlib.util.spec_from_file_location("verify_stage1", MODULE_PATH)
verify = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verify)


def manifest():
    return json.loads(MANIFEST.read_text())


def collect(check):
    """Run one check function and return its named results."""
    results = {}
    check(manifest(), lambda name, value, detail=None: results.__setitem__(name, value))
    return results


def test_internal_arithmetic():
    assert all(v for v in collect(verify.check_arithmetic).values())


def test_hadamard_bound_justifies_the_prime():
    assert all(v for v in collect(verify.check_hadamard).values())
    # the shorthand "prime > 2^7" and the coded "p^2 > (n+1)^d" agree
    assert math.isqrt((3 + 1) ** 7) == 128


def test_full_dimension_certificate_regenerates():
    assert all(v for v in collect(verify.check_full_dimension).values())


def test_every_farkas_removal_replays():
    assert all(v for v in collect(verify.check_farkas).values())


def test_every_facet_witness_is_genuine():
    assert all(v for v in collect(verify.check_facets).values())


def test_ressayre_rows_are_internally_consistent():
    assert all(v for v in collect(verify.check_ressayre_rows).values())


def test_retained_rows_match_the_repository_table():
    """The generated system must equal the stored Altunbulak-Klyachko table."""
    from gpc_census.constraints import constraints

    stored = {(tuple(q["coeffs"]), q["rhs"]) for q in constraints(3, 7)["inequalities"]}
    generated = {
        (tuple(r["constraint"]["coeffs"]), r["constraint"]["rhs"])
        for r in manifest()["retained_rows"]
    }
    assert len(generated) == 4
    assert stored == generated


def test_closed_flat_lattice_reproduces_on_a_small_rank():
    """Exercise the exact rational enumerator where it is cheap.

    The rank-7 lattice is verified by the script under --flats; here the same
    code path runs on wedge^3 C^6, whose top layer must again be the count of
    admissible hyperplanes.
    """
    counts = verify.closed_flat_counts(3, 6)
    assert counts[0] == 1
    assert counts[1] == math.comb(6, 3) == 20
    assert len(counts) == 6
    assert all(c > 0 for c in counts)


def test_farkas_multipliers_are_all_nonnegative():
    """A Farkas certificate with a negative multiplier proves nothing."""
    reduction = manifest()["ordered_slice_reduction_certificate"]
    for removal in reduction["removals"]:
        for key in ("candidate_terms", "structural_terms"):
            for term in removal[key]:
                assert Fraction(*term["multiplier"]) >= 0
        assert Fraction(*removal["rhs_slack"]) >= 0


def test_structural_rows_are_the_ordered_trace_slice():
    structural = manifest()["ordered_slice_reduction_certificate"]["structural_rows"]
    assert len(structural) == 8
    assert structural[0] == {"coeffs": [1, 0, 0, 0, 0, 0, 0], "rhs": 1}
    assert structural[-1] == {"coeffs": [0, 0, 0, 0, 0, 0, -1], "rhs": 0}
    for i in range(1, 7):
        row = structural[i]["coeffs"]
        assert row[i - 1] == -1 and row[i] == 1
        assert structural[i]["rhs"] == 0


def test_cyclic_schubert_cross_check_was_deferred_during_screening():
    """The screening counts and the retained rows tell different stories.

    Screening deferred all 49 rows and cross-checked none; the coefficient-one
    certificates appear only on the four retained rows afterwards. Pinning this
    so the manifest is not read as though 49 rows carried Schubert evidence.
    """
    counts = manifest()["screening_counts"]
    assert counts["cyclic_schubert_deferred_rows"] == 49
    assert counts["cyclic_schubert_cross_checked_rows"] == 0
    assert counts["cyclic_schubert_coefficient_one_rows"] == 0
    for row in manifest()["retained_rows"]:
        assert row["cyclic_schubert_cross_check"]["status"] == "coefficient_one"


def test_generated_stage1_manifest_is_current():
    """The tripwire this replaces asserted the generation package was absent.

    It has landed, so the check becomes the one it was standing in for: run
    the real generator and require the checked-in manifest back byte for byte.
    Byte identity, not semantic equality, because it also catches changed key
    ordering, dropped certificates and accidental float serialization.
    """
    import subprocess
    import sys
    import tempfile

    generator = ROOT / "scripts" / "generate_fixed_n3_reference_rank.py"
    assert generator.exists()
    assert (ROOT / "src" / "gpc_census" / "generation").is_dir()
    with tempfile.TemporaryDirectory() as directory:
        generated = pathlib.Path(directory) / "stage1_ressayre_3_7.json"
        subprocess.run(
            [sys.executable, str(generator), "--rank", "7", "--out", str(generated)],
            check=True,
            capture_output=True,
        )
        assert generated.read_bytes() == MANIFEST.read_bytes()


def test_completeness_is_now_reproducible_in_this_repository():
    """Every import of the generation package must resolve.

    Two of the fifteen modules arrived three separate times without their
    siblings; an import-time check is what distinguishes a package that is
    present from one that is merely partly present.
    """
    import importlib

    for name in (
        "admissible_flats", "bdr", "candidate_pipeline", "candidate_system",
        "cyclic_schubert", "exterior", "full_dimension", "model",
        "redundancy", "reference_generator", "ressayre", "schubert",
        "series", "system",
    ):
        importlib.import_module(f"gpc_census.generation.{name}")
    assert manifest()["proof_status"]["system_completeness"] == "proved"
