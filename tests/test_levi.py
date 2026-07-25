"""Levi-orbit tests.

REGRESSION GUARD: the original `levi_theorem_consistent` reduced to
`slater == slater` and could never fail.  test_theorem_check_is_falsifiable
pins that it can now return False on malformed input.
"""
from gpc_census.levi import (
    block_sector_dimensions,
    full_levi_line_possible,
    levi_dimension,
    spectrum_multiplicities,
    static_levi_record,
)


def test_static_branching_empty_or_full_only():
    sectors = block_sector_dimensions(3, [2, 3])
    lines = [r for r in sectors if r["projective_line_possible"]]
    assert lines == [{"block_occupations": [0, 3], "dimension": 1,
                      "projective_line_possible": True}]
    assert full_levi_line_possible(3, [2, 3])
    assert not full_levi_line_possible(2, [4])


def test_static_record_slater_and_correlated():
    slater = {"system": "(2,4)", "index": 0, "integer_form": [1, 1, 0, 0],
              "denominator": 1}
    correlated = {"system": "(2,4)", "index": 1, "integer_form": [1, 1, 1, 1],
                  "denominator": 2}
    a = static_levi_record(slater)
    b = static_levi_record(correlated)
    assert a["multiplicities"] == [2, 2]
    assert a["levi_dim"] == 8
    assert a["representation_contains_full_levi_line"]
    assert a["target_full_levi_line_possible"]
    assert a["slater_spectrum"]
    assert a["levi_theorem_consistent"]

    assert b["multiplicities"] == [4]
    assert b["levi_dim"] == 16
    # FIXED: Lambda^2 C^4 has no one-dimensional sector, so this is False.
    # The previous suite asserted True here while another test asserted the
    # negation of the same quantity; both could not pass.
    assert not b["representation_contains_full_levi_line"]
    assert not b["target_full_levi_line_possible"]
    assert not b["slater_spectrum"]
    assert b["levi_theorem_consistent"]
    assert b["hard_obstruction"]


def test_theorem_check_is_falsifiable():
    """Malformed record: sum(integer_form) != n * denominator."""
    bad = {"system": "(1,4)", "index": 0, "integer_form": [3, 3, 0, 0],
           "denominator": 3}
    assert static_levi_record(bad)["levi_theorem_consistent"] is False


def test_levi_dimension():
    assert levi_dimension([2, 3]) == 13
    assert spectrum_multiplicities([2, 2, 1, 0, 0])[0] == [2, 1, 2]


def test_stabilizer_is_unitarily_invariant():
    """REGRESSION GUARD for the conjugate-natural-orbital bug.

    one_rdm returns <a_p^dag a_q> = rho_1^TRANSPOSE, so its eigenvectors are
    the conjugates of the natural orbitals.  Building Levi generators from
    those conjugates gives wrong orbit/stabilizer ranks for any state with
    complex amplitudes (real states are unaffected, which is why the census
    looked consistent at 798/799).  The stabilizer dimension is exactly
    invariant along the U(d)-orbit psi -> (Lambda^N U) psi; this test
    conjugates a Slater determinant by a COMPLEX unitary and pins stab == 8.
    """
    import numpy as np
    from scipy.linalg import qr

    from gpc_census.levi import determinant_basis, fiber_and_orbit

    basis = determinant_basis(2, 4)
    psi = np.zeros(len(basis), complex)
    psi[basis.index((0, 1))] = 1.0
    assert fiber_and_orbit(psi, basis, 2, 4)["projective_stabilizer_dim"] == 8

    rng = np.random.default_rng(3)
    z = (rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))) / np.sqrt(2)
    q, r = qr(z)
    U = q * (np.diag(r) / abs(np.diag(r)))
    phi = np.zeros_like(psi)
    for a, T in enumerate(basis):
        phi[a] = sum(np.linalg.det(U[np.ix_(list(T), list(S))]) * psi[i]
                     for i, S in enumerate(basis) if abs(psi[i]) > 1e-14)
    phi /= np.linalg.norm(phi)
    assert fiber_and_orbit(phi, basis, 2, 4)["projective_stabilizer_dim"] == 8


def test_one_rdm_is_rho1_transpose():
    """Pin the convention so the conj() in fiber_and_orbit stays justified:
    under psi -> (Lambda^N U) psi, one_rdm transforms as conj(U) rho U^T."""
    import numpy as np
    from scipy.linalg import qr

    from gpc_census.levi import annihilation_matrices, determinant_basis, one_rdm

    basis = determinant_basis(2, 4)
    psi = np.zeros(len(basis), complex)
    psi[basis.index((0, 1))] = 1.0
    rng = np.random.default_rng(5)
    z = (rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))) / np.sqrt(2)
    q, r = qr(z)
    U = q * (np.diag(r) / abs(np.diag(r)))
    phi = np.zeros_like(psi)
    for a, T in enumerate(basis):
        phi[a] = sum(np.linalg.det(U[np.ix_(list(T), list(S))]) * psi[i]
                     for i, S in enumerate(basis) if abs(psi[i]) > 1e-14)
    phi /= np.linalg.norm(phi)
    ann = annihilation_matrices(2, 4, basis)
    rho, rhop = one_rdm(psi, ann), one_rdm(phi, ann)
    assert np.linalg.norm(rhop - U.conj() @ rho @ U.T) < 1e-10
    assert np.linalg.norm(rhop - U @ rho @ U.conj().T) > 1e-2
