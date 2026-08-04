from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from dataclasses import replace
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from gpc_census.constraints import constraints
from gpc_census.generation import generate_3_6
from gpc_census.generation.exterior import exterior_root_action, exterior_weights
from gpc_census.generation.model import IntegralConstraint
from gpc_census.generation.ressayre import verify_ressayre_certificate
from gpc_census.generation.schubert import verify_borland_dennis_schubert_certificates
from gpc_census.generation.system import restrict_to_affine_hull


ROOT = Path(__file__).parents[1]
ARTIFACT = ROOT / "results" / "data" / "stage1_ressayre_3_6.json"


@lru_cache(maxsize=1)
def generated():
    return generate_3_6()


def test_exterior_root_action_sign_convention():
    weights = exterior_weights(3, 6)
    source = (1, 1, 1, 0, 0, 0)
    target = (1, 1, 0, 1, 0, 0)
    assert len(weights) == 20
    assert exterior_root_action(source, (2, 3)) == (target, 1)
    assert exterior_root_action(target, (2, 3)) is None


def test_exhaustive_rank6_ressayre_counts_and_replay():
    result = generated().enumeration
    assert result.admissible_hyperplane_count == 362
    assert result.oriented_candidate_count == 724
    assert result.trace_candidate_count == 64
    assert len(result.certificates) == 21
    assert all(verify_ressayre_certificate(3, 6, row) for row in result.certificates)


def test_nontrivial_ressayre_sign_guard():
    result = generated().enumeration
    target = next(
        row
        for row in result.certificates
        if row.h == (-1, -1, 1, -1, 1, 1) and row.z == -1
    )
    weights = exterior_weights(3, 6)
    below = tuple(weights[index] for index in target.below_indices)
    assert below == ((1, 1, 0, 1, 0, 0),)
    assert target.negative_roots == ((2, 3),)
    assert target.determinant_value == 1
    assert not any(
        row.h == tuple(-value for value in target.h) and row.z == 1
        for row in result.certificates
    )


def test_rank6_system_is_derived_without_lookup():
    result = generated()
    assert result.equalities == (
        IntegralConstraint((1, 0, 0, 0, 0, 1), 1),
        IntegralConstraint((0, 1, 0, 0, 1, 0), 1),
        IntegralConstraint((0, 0, 1, 1, 0, 0), 1),
    )
    assert result.inequalities == (
        IntegralConstraint((0, 0, 0, 1, -1, -1), 0),
    )


def test_rank6_inner_outer_sandwich_closes_exactly():
    result = generated()
    expected = (
        (Fraction(1), Fraction(1), Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
        (
            Fraction(1),
            Fraction(1, 2),
            Fraction(1, 2),
            Fraction(1, 2),
            Fraction(1, 2),
            Fraction(0),
        ),
        (
            Fraction(3, 4),
            Fraction(3, 4),
            Fraction(1, 2),
            Fraction(1, 2),
            Fraction(1, 4),
            Fraction(1, 4),
        ),
        (Fraction(1, 2),) * 6,
    )
    assert result.outer_vertices == expected
    assert result.inner_point_count == 20
    assert result.outer_vertices_in_inner_cloud


def test_schubert_and_ressayre_routes_agree():
    result = generated()
    certificates = result.schubert_cross_checks
    assert verify_borland_dennis_schubert_certificates(certificates)
    assert len(certificates) == 4
    assert all(certificate.schubert_coefficient == 1 for certificate in certificates)
    assert tuple(certificate.constraint for certificate in certificates[:3]) == result.equalities
    assert restrict_to_affine_hull(certificates[3].constraint, result.equalities) in result.inequalities


def test_generated_system_matches_lookup_only_as_post_generation_regression():
    result = generated()
    legacy = constraints(3, 6)
    expected_inequalities = {
        IntegralConstraint(tuple(row["coeffs"]), row["rhs"])
        for row in legacy["inequalities"]
    }
    expected_equalities = {
        IntegralConstraint(tuple(row["coeffs"]), row["rhs"])
        for row in legacy["equalities"]
    }
    assert set(result.inequalities) == expected_inequalities
    assert set(result.equalities) == expected_equalities


def test_determinant_certificate_rejects_tampering():
    certificate = next(
        row for row in generated().enumeration.certificates if row.determinant_value != 1
    )
    tampered = replace(certificate, determinant_value=certificate.determinant_value + 1)
    assert not verify_ressayre_certificate(3, 6, tampered)


def test_checked_in_artifact_is_regenerated_byte_for_byte():
    module_path = ROOT / "scripts" / "stage1_generator.py"
    spec = importlib.util.spec_from_file_location("stage1_generator", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    expected = json.dumps(module.payload(), indent=2, sort_keys=True) + "\n"
    assert ARTIFACT.read_text() == expected


def test_standalone_artifact_verifier():
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "verify_stage1_3_6_standalone.py")],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "verification passed" in result.stdout
