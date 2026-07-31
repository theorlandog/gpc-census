"""Exact pair-coordinate and v_B reduced-fiber symmetry certificates."""
import copy
import json
import pathlib
import sys

import numpy as np
import pytest
import sympy as sp

from gpc_census.fiber import one_rdm
from gpc_census.gauge import _pair_vectors

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
DATA = ROOT / "results" / "data"

VB_DETS = [
    (0, 1, 2, 4),
    (0, 1, 2, 8),
    (0, 1, 3, 5),
    (0, 2, 3, 6),
    (0, 3, 4, 7),
    (0, 3, 7, 8),
    (1, 3, 6, 8),
    (2, 3, 4, 8),
]
VB_WEIGHTS = (1, 8, 4, 3, 2, 2, 1, 2)
VB_KERNEL = (1, -1, 0, 0, -1, 1, 0, 0)
VB_SPECTRUM = np.array([20, 14, 14, 14, 14, 4, 4, 4, 4]) / 23


def _family_amplitudes(t):
    """A family member, built without touching the generator's code path."""
    weights = [(VB_WEIGHTS[i] + VB_KERNEL[i] * t) / 23 for i in range(8)]
    q = (1 + t) * (8 - t) * (2 - t) * (2 + t)
    p = 3 + 7 * t - 2 * t * t
    amplitudes = np.array([np.sqrt(w) for w in weights], dtype=complex)
    amplitudes[0] *= np.exp(1j * np.arccos(p / (2 * np.sqrt(q))))
    return amplitudes


@pytest.fixture(scope="module")
def artifact():
    return json.loads((DATA / "fiber_symmetries.json").read_text())


def test_every_shipped_support_has_an_exact_pair_minor(artifact):
    from scripts.verify_fiber_symmetries_standalone import verify

    assert verify() == []
    summary = artifact["summary"]
    assert summary["records"] == 799
    assert summary["full_column_rank"] == 799
    assert summary["positive_incidence_kernel_records"] == 63
    assert summary["positive_incidence_kernel_transport_classes"] == 39


def test_pair_certificate_artifact_is_current(artifact):
    from scripts.fiber_symmetries import build

    assert build() == artifact


def test_tampered_pair_minor_is_rejected(tmp_path, artifact):
    from scripts.verify_fiber_symmetries_standalone import verify

    bad = copy.deepcopy(artifact)
    bad["records"][0]["minor_determinant"] = 0
    path = tmp_path / "tampered.json"
    path.write_text(json.dumps(bad))
    errors = verify(artifact_path=path)
    assert any("determinant mismatch" in error for error in errors)


def test_vb_cubic_two_rdm_moment_is_exact_and_injective():
    from scripts.fiber_symmetries import _vb_symbolic_moments

    t, moments, _ = _vb_symbolic_moments()
    expected = [
        sp.Integer(6),
        sp.Rational(1438, 529),
        6 * (2899 - 8 * t) / 12167,
    ]
    assert all(sp.simplify(got - want) == 0 for got, want in zip(moments, expected))
    assert sp.diff(moments[2], t) == sp.Rational(-48, 12167)


def test_vb_galois_walls_are_not_one_body_unitarily_equivalent(artifact):
    vb = artifact["vb"]
    assert vb["galois_walls"]["difference"] == "1728*sqrt(35)/1034195"
    assert vb["support_preserving_block_permutation_automorphisms"] == 1
    assert "remains open" in vb["time_reversal"]["full_stabilizer_status"]


def test_vb_moments_separate_the_exact_phase_gauge_quotient_tangent():
    from scripts.fiber_symmetries import _vb_exact_tangent_moment_rank

    certificate = _vb_exact_tangent_moment_rank()
    assert certificate["fixed_spectrum_constraint_rank"] == 7
    assert certificate["fixed_spectrum_kernel_dimension"] == 9
    assert certificate["orbital_phase_gauge_rank"] == 7
    assert certificate["quotient_tangent_dimension"] == 2
    assert certificate["combined_rank"] == 9
    assert certificate["combined_minor_determinant"] == (
        "299925504*sqrt(8211)/1035662780341225"
    )
    assert "first-order exact statement only" in certificate["scope_limit"]


def test_vb_moments_reproduce_independently_of_the_generator():
    """POWER TEST. A second derivation path for the shipped v_B identities.

    Everything else about v_B in this file routes through the generator's
    symbolic reduction. This rebuilds the family from the published weights and
    coherence equation, contracts it with the repository's own pair-vector code,
    and checks the artifact's rational moment formulas numerically.
    """
    artifact = json.loads((DATA / "fiber_symmetries.json").read_text())
    t = sp.Symbol("t")
    trace_2 = sp.sympify(artifact["vb"]["two_rdm_moments"]["trace_2"])
    trace_3 = sp.sympify(artifact["vb"]["two_rdm_moments"]["trace_3"])

    for value in (-0.5, 0.0, 0.5, 1.0, 1.5):
        c = _family_amplitudes(value)
        gram = _pair_vectors(c, VB_DETS, 9)
        gamma2 = gram.conj() @ gram.T
        powers = [
            float(np.trace(np.linalg.matrix_power(gamma2, k)).real)
            for k in (1, 2, 3)
        ]
        assert abs(powers[0] - 6) < 1e-12
        assert abs(powers[1] - float(trace_2)) < 1e-12
        assert abs(powers[2] - float(trace_3.subs(t, value))) < 1e-12


def test_vb_family_is_fixed_spectrum_and_not_fixed_diagonal():
    """The scope word is load bearing: rho keeps its spectrum, not its zeros."""
    off_diagonals = []
    for value in (0.0, 0.5, 1.0):
        rho = one_rdm(_family_amplitudes(value), VB_DETS, 9)
        eigenvalues = np.sort(np.linalg.eigvalsh((rho + rho.conj().T) / 2))[::-1]
        assert np.max(np.abs(eigenvalues - np.sort(VB_SPECTRUM)[::-1])) < 1e-12
        off_diagonals.append(rho[4, 8])
    assert abs(off_diagonals[0] - (-5 / 92 - 1j * np.sqrt(119) / 92)) < 1e-12
    assert abs(off_diagonals[0] - off_diagonals[1]) > 1e-6


def test_both_galois_walls_are_physical_attainers():
    """Pins the RESEARCH.md correction: both roots, not one, are physical."""
    root = np.sqrt(35)
    for value in (7 / 17 - 18 * root / 85, 7 / 17 + 18 * root / 85):
        weights = [(VB_WEIGHTS[i] + VB_KERNEL[i] * value) / 23 for i in range(8)]
        assert min(weights) > 0
        q = (1 + value) * (8 - value) * (2 - value) * (2 + value)
        p = 3 + 7 * value - 2 * value * value
        assert abs(abs(p / (2 * np.sqrt(q))) - 1) < 1e-12
        c = np.array([np.sqrt(w) for w in weights], dtype=complex)
        c[0] *= np.sign(p)
        rho = one_rdm(c, VB_DETS, 9)
        eigenvalues = np.sort(np.linalg.eigvalsh((rho + rho.conj().T) / 2))[::-1]
        assert np.max(np.abs(eigenvalues - np.sort(VB_SPECTRUM)[::-1])) < 1e-12
