"""Guards for the mechanism-grammar synthesis and its five input artifacts."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

# Spelled literally so the results/data completeness audit can see the coverage.
ARITHMETIC_LEVEL_AUDIT = DATA / "arithmetic_level_audit.json"
FULL_RANK78_COMPRESSION_AUDIT = DATA / "full_rank78_compression_audit.json"
GENERATOR_COMPRESSION_TRIAL = DATA / "generator_compression_trial.json"
LEVI_FUSION_PILOT = DATA / "levi_fusion_pilot.json"
LEVI_FUSION_HELDOUT_RESULTS = DATA / "levi_fusion_heldout_results.json"
MECHANISM_GRAMMAR_SYNTHESIS = DATA / "mechanism_grammar_synthesis.json"

INPUT_ARTIFACTS = (
    ARITHMETIC_LEVEL_AUDIT,
    FULL_RANK78_COMPRESSION_AUDIT,
    GENERATOR_COMPRESSION_TRIAL,
    LEVI_FUSION_PILOT,
    LEVI_FUSION_HELDOUT_RESULTS,
)


def _load(name: str):
    script = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _artifact(path: pathlib.Path) -> dict:
    return json.loads(path.read_text())


def _run(script: str) -> str:
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / script)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_every_synthesis_input_is_present_and_declared():
    """The synthesis must name exactly the artifacts it actually reconciles."""
    declared = set(_artifact(MECHANISM_GRAMMAR_SYNTHESIS)["inputs"])
    for path in INPUT_ARTIFACTS:
        assert path.is_file(), path
        assert isinstance(_artifact(path), dict)
        assert f"results/data/{path.name}" in declared
    assert declared == {f"results/data/{path.name}" for path in INPUT_ARTIFACTS}


def test_arithmetic_audit_regenerates_exactly():
    """The stored audit must be reproducible from the two Levi artifacts."""
    module = _load("arithmetic_level_audit")
    payload = module.build_payload(module.PILOT, module.HELDOUT)
    assert payload == _artifact(ARITHMETIC_LEVEL_AUDIT)


def test_arithmetic_audit_headline_counts():
    summary = _artifact(ARITHMETIC_LEVEL_AUDIT)["summary"]
    assert summary["mechanisms"] == 83
    assert summary["complete_arithmetic_progressions"] == 78
    assert summary["one_point_deleted_progressions"] == 5
    assert summary["more_than_one_point_missing"] == 0
    assert summary["step_divides_N"] == 83


def test_divisibility_lemma_hypotheses_hold():
    """g | N is entailed, so its hypotheses must hold at every mechanism."""
    module = _load("arithmetic_level_audit")
    mechanisms = module.load_mechanisms(module.PILOT, module.HELDOUT)
    assert len(mechanisms) == 83
    for mechanism in mechanisms:
        witness = module.divisibility_witness(mechanism)
        assert all(witness.values()), (mechanism, witness)


def test_no_multi_hole_survivor_in_any_population():
    """The surviving conjecture's core count, across all three scopes."""
    synthesis = _artifact(MECHANISM_GRAMMAR_SYNTHESIS)
    budget = synthesis["hole_budget"]
    assert budget["multi_hole_survivors_anywhere"] == 0
    sizes = {}
    for population in budget["populations"]:
        assert population["multi_hole"] == 0
        sizes[population["population"]] = population["size"]
    assert sorted(sizes.values()) == [83, 537, 7156]


def test_one_hole_onset_is_rank_nine():
    onset = _artifact(MECHANISM_GRAMMAR_SYNTHESIS)["hole_onset"]
    assert onset["onset_rank"] == 9
    assert onset["published_one_hole_mechanisms"] == 5
    assert onset["published_one_hole_ranks"] == [9, 10]
    assert onset["rank7_rank8_one_hole_hall_feasible"] == 7
    assert onset["rank7_rank8_one_hole_determinant_nonzero"] == 0
    assert onset["rank9_one_hole_modular_nonzero_certificates"] == 3


def test_rank8_bounded_and_full_searches_agree():
    agreement = _artifact(MECHANISM_GRAMMAR_SYNTHESIS)[
        "rank8_bounded_full_agreement"
    ]
    assert agreement["shape"] == [-2, -1, -1, -1, 0, 0, 0, 2]
    assert agreement["normalized_levels"] == [0, 1, 2, 4]
    assert agreement["missing_levels"] == [3]
    assert agreement["hall_feasible_placements"] == 7


def test_block_count_refutation_witness():
    """r <= 2N fails while f and dim W_0 stay small at the same row."""
    module = _load("verify_mechanism_grammar_synthesis")
    descriptor = module.zero_grade_descriptor(module.P1_WITNESS_TAU, 4)
    assert descriptor["blocks"] == 10
    assert descriptor["blocks"] > 2 * 4
    assert descriptor["fusion_rank"] == 9
    assert descriptor["zero_grade_dimension"] == 9
    assert descriptor["ambient_dimension"] == 210


def test_width_is_not_a_facetness_discriminator():
    width = _artifact(MECHANISM_GRAMMAR_SYNTHESIS)["fusion_width"]
    assert width["known_mechanism_max_optimal_width"] == 5
    assert width["full_rank7_rank8_max_optimal_width"] == 4
    assert width["rank8_hall_zero_candidates_at_width_four"] == 2393
    assert width["status"] == "evaluation_compression_only"


def test_tomography_injective_on_full_populations():
    tomography = _artifact(MECHANISM_GRAMMAR_SYNTHESIS)["tomography"]
    assert tomography["known_mechanisms_separated"] == 83
    assert tomography["known_mechanism_collisions"] == 0
    counts = {
        entry["system"]: entry["nonstructural_survivors"]
        for entry in tomography["full_population_injectivity"]
    }
    assert counts == {"(3,7)": 547, "(3,8)": 6605}
    for entry in tomography["full_population_injectivity"]:
        assert entry["nonstructural_survivors"] == entry[
            "distinct_singleton_signatures"
        ]


def test_synthesis_verifier_replays():
    output = _run("verify_mechanism_grammar_synthesis.py")
    assert "mechanism-grammar synthesis: PASS" in output


def test_rank78_compact_verifier_replays():
    output = _run("verify_full_rank78_compression_audit.py")
    assert "full rank-7/rank-8 compact certificate verifier: PASS" in output


def test_generator_compression_trial_verifier_replays():
    output = _run("verify_generator_compression_trial.py")
    assert "compression trial verifier: PASS" in output
