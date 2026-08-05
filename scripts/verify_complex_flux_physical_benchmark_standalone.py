#!/usr/bin/env python3
"""Independent replay of the complex-flux physical benchmark.

This verifier imports no ``gpc_census`` module.  It rebuilds the spinless
three-fermion Hamiltonian directly in occupation space, forms the exterior
power basis transformation from minors, and checks the selected benchmark in
``conditional_2rdm_program.json``.  The calculation is numerical, but every
model input is an exact integer or rational Gaussian number.
"""
from __future__ import annotations

import argparse
import json
import pathlib
from itertools import combinations

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[1]
CARRIER = ((0, 1, 2), (0, 3, 4), (1, 3, 5))


def apply_ops(det, ops):
    occ = list(det)
    sign = 1
    for kind, orbital in reversed(ops):
        if kind == "a":
            if orbital not in occ:
                return None
            position = occ.index(orbital)
            sign *= (-1) ** position
            occ.pop(position)
        else:
            if orbital in occ:
                return None
            position = sum(x < orbital for x in occ)
            sign *= (-1) ** position
            occ.insert(position, orbital)
    return tuple(occ), sign


def hamiltonian():
    """Direct occupation-space t-V ring Hamiltonian at the frozen input."""
    dets = list(combinations(range(6), 3))
    index = {det: i for i, det in enumerate(dets)}
    phase = (4.0 + 3.0j) / 5.0
    h1 = np.zeros((6, 6), dtype=complex)
    for p in range(6):
        h1[p, p] = (7.0 / 20.0) * p
        q = (p + 1) % 6
        edge_phase = phase if q == 0 else 1.0
        h1[p, q] += -edge_phase
        h1[q, p] += -np.conj(edge_phase)

    matrix = np.zeros((len(dets), len(dets)), dtype=complex)
    for column, det in enumerate(dets):
        for p in range(6):
            for q in range(6):
                if h1[p, q] == 0:
                    continue
                got = apply_ops(det, (("c", p), ("a", q)))
                if got is not None:
                    target, sign = got
                    matrix[index[target], column] += sign * h1[p, q]
        matrix[column, column] += 2.0 * sum(
            p in det and ((p + 1) % 6) in det for p in range(6)
        )
    return matrix, dets


def one_rdm(coefficients, dets):
    index = {det: i for i, det in enumerate(dets)}
    gamma = np.zeros((6, 6), dtype=complex)
    for column, det in enumerate(dets):
        for p in range(6):
            for q in range(6):
                got = apply_ops(det, (("c", p), ("a", q)))
                if got is None:
                    continue
                target, sign = got
                gamma[p, q] += (
                    np.conj(coefficients[index[target]])
                    * sign
                    * coefficients[column]
                )
    return gamma


def exterior_matrix(one_body, dets):
    """Matrix of the third exterior power in the determinant ordering."""
    matrix = np.zeros((len(dets), len(dets)), dtype=complex)
    for row, target in enumerate(dets):
        for column, source in enumerate(dets):
            matrix[row, column] = np.linalg.det(
                one_body[np.ix_(target, source)]
            )
    return matrix


def close(actual, expected, tolerance=5e-10):
    return bool(np.max(np.abs(np.asarray(actual) - np.asarray(expected))) < tolerance)


def verify(path: pathlib.Path) -> dict:
    artifact = json.loads(path.read_text())
    benchmark = artifact["projection_and_quasipinning"][
        "complex_flux_physical_benchmark"
    ]
    recorded = benchmark["selected_scan_point"]

    H, dets = hamiltonian()
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    psi = eigenvectors[:, 0]
    gamma = one_rdm(psi, dets)
    occupations, U = np.linalg.eigh(gamma)
    order = np.argsort(occupations)[::-1]
    occupations = np.real(occupations[order])
    U = U[:, order]

    # Since gamma[p,q] = <a_p^dag a_q>, the physical orbital columns are
    # conjugate(U).  W maps new determinant coefficients to the old basis.
    B = U.conj()
    W = exterior_matrix(B, dets)
    psi_no = W.conj().T @ psi
    H_no = W.conj().T @ H @ W
    gamma_no = one_rdm(psi_no, dets)
    offdiag = float(np.max(np.abs(gamma_no - np.diag(np.diag(gamma_no)))))

    index = {det: i for i, det in enumerate(dets)}
    carrier_indices = [index[det] for det in CARRIER]
    carrier_amplitudes = psi_no[carrier_indices]
    on_carrier = float(np.vdot(carrier_amplitudes, carrier_amplitudes).real)
    off_carrier = max(0.0, 1.0 - on_carrier)
    weights = np.abs(carrier_amplitudes) ** 2 / on_carrier
    phi = carrier_amplitudes / np.sqrt(on_carrier)
    Hc = H_no[np.ix_(carrier_indices, carrier_indices)]
    cycle = Hc[0, 1] * Hc[1, 2] * np.conj(Hc[0, 2])
    carrier_energy = float(np.vdot(phi, Hc @ phi).real)
    true_energy = float(eigenvalues[0])
    spectral_width = float(eigenvalues[-1] - eigenvalues[0])
    eps = np.sqrt(off_carrier)
    enlargement = spectral_width * (eps + eps * eps)
    top_three = {
        det for _, det in sorted(
            ((float(abs(a) ** 2), det) for a, det in zip(psi_no, dets)),
            reverse=True,
        )[:3]
    }

    # Negative control for the convention bug.  This was invisible at zero
    # flux because the natural orbitals could be chosen real.
    wrong = exterior_matrix(U.conj().T, dets) @ psi
    wrong_gamma = one_rdm(wrong, dets)
    wrong_offdiag = float(np.max(np.abs(
        wrong_gamma - np.diag(np.diag(wrong_gamma))
    )))

    checks = {
        "exact_input_recorded": benchmark["exact_model_input"] == {
            "particles": 3,
            "orbitals": 6,
            "hopping": "1",
            "nearest_neighbour_interaction": "2",
            "onsite_gradient": "7/20",
            "boundary_phase": "(4+3*i)/5",
            "boundary_phase_modulus_squared": "1",
            "time_reversal_is_broken": True,
        },
        "exterior_transform_unitary": close(W.conj().T @ W, np.eye(len(dets))),
        "transformed_eigenstate": close(H_no @ psi_no, true_energy * psi_no),
        "nondegenerate_ground_state": eigenvalues[1] - eigenvalues[0] > 1.0,
        "nondegenerate_occupations": min(np.diff(occupations) * -1) > 1e-3,
        "natural_basis_diagonal": offdiag < 1e-12,
        "dominant_carrier": top_three == set(CARRIER),
        "off_carrier_below_one_percent": off_carrier < 0.01,
        "genuine_carrier_flux": abs(cycle.imag) > 1e-3,
        "quasipinning_bound_direct": abs(true_energy - carrier_energy) <= enlargement,
        "wrong_conjugation_detected": wrong_offdiag > 0.1,
        "artifact_occupations_match": close(occupations, recorded["occupations"]),
        "artifact_carrier_weights_match": close(weights, recorded["carrier_weights"]),
        "artifact_cycle_match": close(
            [cycle.real, cycle.imag], recorded["carrier_cycle_product"]
        ),
        "artifact_scalar_metrics_match": close(
            [
                eigenvalues[1] - eigenvalues[0],
                min(np.diff(occupations) * -1),
                offdiag,
                off_carrier,
                carrier_energy,
                true_energy,
                enlargement,
            ],
            [
                recorded["ground_state_gap"],
                recorded["min_occupation_gap"],
                recorded["natural_basis_offdiagonal_max"],
                recorded["off_carrier_weight"],
                recorded["normalized_carrier_energy"],
                recorded["true_energy"],
                recorded["enlargement_half_width"],
            ],
        ),
    }
    checks = {name: bool(passed) for name, passed in checks.items()}
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise AssertionError(f"complex-flux replay failed: {failed}")
    return {
        "checks": checks,
        "ground_state_gap": float(eigenvalues[1] - eigenvalues[0]),
        "minimum_occupation_gap": float(min(np.diff(occupations) * -1)),
        "off_carrier_weight": off_carrier,
        "carrier_cycle_imaginary_part": float(cycle.imag),
        "natural_basis_offdiagonal_max": offdiag,
        "wrong_conjugation_offdiagonal_max": wrong_offdiag,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "artifact",
        nargs="?",
        type=pathlib.Path,
        default=ROOT / "results/data/conditional_2rdm_program.json",
    )
    args = parser.parse_args()
    result = verify(args.artifact)
    print(f"PASS {len(result['checks'])}/{len(result['checks'])} independent checks")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
