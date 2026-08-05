#!/usr/bin/env python3
"""Steps 4 to 7 of the conditional 2-RDM program.

Block 4  FLUX AND CURRENT. The support oracle applied to H(Phi) + s J(Phi)
         gives the exact convex energy/pair-current joint body, the envelope
         currents J_± = -dF^±/dPhi, and the discriminant fluxes where the
         optimizer branches exchange.

Block 5  PROJECTED HAMILTONIAN. Slater-Condon projection of a real
         second-quantized Hamiltonian onto the pinned carrier, over a
         parameter scan of a recognizable model (spinless t-V ring), with an
         FCIDUMP round trip.

Block 6  QUASIPINNING. The enlargement bound, with its two independent terms
         (off-carrier weight, and the shift of the fiber weights) and a
         numerical check that the true expectation lands inside.

Block 7  ATLAS. Every certified census support classified by carrier order,
         two-body accessibility of its coherences, and whether its
         compatibility body is exactly spectrahedral or only outer-bounded.

Usage:
    python scripts/conditional_2rdm_program.py [-o results/data/conditional_2rdm_program.json]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import tempfile
from itertools import combinations

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402

from gpc_census import elliptope as el  # noqa: E402
from gpc_census import pinned_physics as pp  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
W = [0.6, 0.25, 0.15]


# --------------------------------------------------------------------------
# block 4
# --------------------------------------------------------------------------

def block_flux_current(n_flux: int = 25, n_dirs: int = 120) -> dict:
    lo, hi = pp.current_range_exact([sp.Rational(3, 5), sp.Rational(1, 4),
                                     sp.Rational(3, 20)])
    sweep = []
    for k in range(n_flux):
        phi = 2 * np.pi * k / n_flux
        H = pp.flux_carrier(phi)
        J = pp.current_carrier(phi)
        Fhi = float(el.support_numeric(H, W, "max"))
        Flo = float(el.support_numeric(H, W, "min"))
        sweep.append({
            "phi": float(phi),
            "F_plus": Fhi, "F_minus": Flo, "hidden_width": Fhi - Flo,
            "current_at_max": float(pp.envelope_current(phi, W, "max")),
            "current_at_min": float(pp.envelope_current(phi, W, "min")),
            "current_range": [float(el.support_numeric(J, W, "min")),
                              float(el.support_numeric(J, W, "max"))],
        })

    # discriminant of the unit-edge family: where the root structure degenerates
    fam = json.loads((ROOT / "results/data/complex_fixed_1rdm_observable_ranges.json")
                     .read_text())["unit_edge_flux_family"]
    u, z = sp.symbols("u z")
    co = [sp.sympify(c) for c in fam["coefficients_descending_in_z"]]
    poly = sum(c * z ** (len(co) - 1 - i) for i, c in enumerate(co))
    disc = sp.factor(sp.discriminant(sp.Poly(poly, z)))
    branch = []
    for r in sp.solve(sp.Eq(disc, 0), u):
        if r.is_real and -1 <= float(r) <= 1:
            branch.append({"u": str(sp.nsimplify(r)), "u_approx": float(r),
                           "phi": float(sp.acos(float(r)))})

    phi0 = 1.0
    joint = pp.joint_body_support(pp.flux_carrier(phi0), pp.current_carrier(phi0),
                                  W, n_directions=n_dirs)
    return {
        "current_operator": "J = -dH/dPhi on the fluxed edge",
        "exact_current_range": [str(lo), str(hi)],
        "exact_current_range_approx": [float(lo), float(hi)],
        "current_range_is_flux_independent": True,
        "why": "the current couples only to the fluxed edge, so its expectation "
               "is 2 Re(conj(c1) c2 * phase) and a free relative phase attains "
               "2 sqrt(w1 w2) at every flux",
        "flux_sweep": sweep,
        "branch_exchange_fluxes": branch,
        "joint_body_support_at_phi_1": [{"theta": t, "value": v} for t, v in joint],
    }


# --------------------------------------------------------------------------
# blocks 5 and 6
# --------------------------------------------------------------------------

def _scan_point(V: float, boundary_phase: complex, gradient: float) -> dict:
    h1, g2 = pp.tv_ring_integrals(
        6, 1.0, V, gradient=gradient, boundary_phase=boundary_phase,
    )
    H, dets = pp.ci_hamiltonian(h1, g2, 6, 3)
    evals, evecs = np.linalg.eigh(H)
    psi = evecs[:, 0]
    gap = float(evals[1] - evals[0])
    gamma = pp.one_rdm_ci(psi, dets, 6)
    lams, U = pp.natural_orbitals(gamma)
    psi_no = pp.natural_orbital_coefficients(psi, dets, U)
    psi_no = psi_no / np.linalg.norm(psi_no)
    gamma_no = pp.one_rdm_ci(psi_no, dets, 6)
    natural_basis_offdiag = float(np.max(np.abs(
        gamma_no - np.diag(np.diag(gamma_no))
    )))
    eps2, w = pp.carrier_weight(psi_no, dets)
    slack = pp.borland_dennis_slack(lams)
    occ_gap = float(min(lams[i] - lams[i + 1] for i in range(5)))

    # the carrier must actually be the dominant support for the pinned
    # description to mean anything
    top = sorted(((float(abs(a) ** 2), D) for a, D in zip(psi_no, dets)),
                 reverse=True)[:3]
    carrier_matches = all(D in pp.DETS for _, D in top)

    # project the Hamiltonian in the natural-orbital basis
    h1_no, g2_no = pp.natural_basis_integrals(h1, g2, U)
    Hc = pp.carrier_projection(h1_no, g2_no)
    cycle = Hc[0, 1] * Hc[1, 2] * np.conj(Hc[0, 2])

    Fhi = float(el.support_numeric(Hc, w, "max"))
    Flo = float(el.support_numeric(Hc, w, "min"))
    true_E = float(np.real(np.vdot(psi, H @ psi)))

    carrier_index = {D: i for i, D in enumerate(dets)}
    carrier_amplitudes = np.array([psi_no[carrier_index[D]] for D in pp.DETS])
    carrier_amplitudes /= np.linalg.norm(carrier_amplitudes)
    carrier_E = float(np.real(np.vdot(carrier_amplitudes, Hc @ carrier_amplitudes)))

    # baseline: no 1-RDM information at all, only the Hamiltonian's spectrum
    unconstrained = float(evals[-1] - evals[0])

    y = pp.dual_optimal_y(Hc, w, "max")
    # The off-carrier error sees the full Hamiltonian, not merely P H P.
    # Centering by the middle of the spectral interval minimizes its norm
    # without changing expectation differences.
    Hnorm = unconstrained / 2
    eps = float(np.sqrt(eps2))
    eps_H = pp.enlargement_bound(Hnorm, eps, dw_l1=0.0, y_inf=float(np.abs(y).max()))
    occupation_nondegenerate = bool(occ_gap > 1e-8)
    pinned_preconditions = bool(
        gap > 1e-8
        and occupation_nondegenerate
        and carrier_matches
        and natural_basis_offdiag < 1e-10
    )
    return {
        "V": V,
        "boundary_phase": [float(boundary_phase.real), float(boundary_phase.imag)],
        "boundary_phase_label": (
            "zero_flux" if abs(boundary_phase.imag) < 1e-15
            else "rational_gaussian_flux_(4+3i)/5"
        ),
        "twist": float(np.angle(boundary_phase)),
        "gradient": gradient,
        "ground_state_gap": gap,
        "occupations": [float(x) for x in lams],
        "min_occupation_gap": occ_gap,
        "natural_occupations_nondegenerate": occupation_nondegenerate,
        "natural_basis_offdiagonal_max": natural_basis_offdiag,
        "borland_dennis_slack": float(slack),
        "off_carrier_weight": float(eps2),
        "slack_over_off_carrier_weight": float(slack / eps2) if eps2 > 1e-14 else None,
        "carrier_is_dominant_support": bool(carrier_matches),
        "pinned_preconditions_pass": pinned_preconditions,
        "carrier_weights": w,
        "carrier_cycle_product": [float(cycle.real), float(cycle.imag)],
        "genuinely_complex_carrier_flux": bool(abs(cycle.imag) > 1e-8),
        "pinned_interval": [Flo, Fhi],
        "pinned_width": Fhi - Flo,
        "normalized_carrier_energy": carrier_E,
        "carrier_energy_inside_pinned_interval": bool(Flo <= carrier_E <= Fhi),
        "true_to_carrier_energy_deviation": abs(true_E - carrier_E),
        "true_energy": true_E,
        "unconstrained_width": unconstrained,
        "width_ratio_unconstrained_over_pinned": (unconstrained / (Fhi - Flo)
                                                  if Fhi > Flo else None),
        "enlargement_half_width": float(eps_H),
        "true_energy_inside_enlarged": bool(Flo - eps_H <= true_E <= Fhi + eps_H),
    }


def block_projection_and_quasipinning() -> dict:
    scan = []
    phases = (1.0 + 0.0j, (4.0 + 3.0j) / 5.0)
    for V in (0.5, 1.0, 2.0, 3.0, 4.0, 6.0):
        for phase in phases:
            scan.append(_scan_point(V, phase, gradient=0.35))

    usable = [r for r in scan if r["pinned_preconditions_pass"]]
    complex_usable = [r for r in usable if r["genuinely_complex_carrier_flux"]]
    selected_complex = next(
        r for r in complex_usable
        if r["V"] == 2.0 and r["boundary_phase_label"].startswith("rational_gaussian")
    )
    # FCIDUMP round trip, on the integrals the scan actually used
    h1, g2 = pp.tv_ring_integrals(6, 1.0, 2.0, twist=0.0, gradient=0.35)
    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / "MODEL.FCIDUMP"
        pp.write_fcidump(path, h1, g2, n_elec=3)
        h1b, g2b, norb = pp.read_fcidump(path)
    roundtrip = bool(norb == 6 and np.allclose(h1, h1b) and np.allclose(g2, g2b))

    # the carrier projection must reproduce the full CI submatrix exactly
    Hfull, dets = pp.ci_hamiltonian(h1, g2, 6, 3)
    idx = {D: i for i, D in enumerate(dets)}
    Hc = pp.carrier_projection(h1, g2)
    sub = np.array([[Hfull[idx[A], idx[B]] for B in pp.DETS] for A in pp.DETS])
    proj_exact = float(np.abs(Hc - sub).max())

    return {
        "model": "spinless t-V ring, 6 orbitals, 3 fermions, on-site gradient "
                 "0.35 to break the reflection symmetry that would otherwise "
                 "force degenerate natural occupations",
        "fcidump_round_trip": roundtrip,
        "carrier_projection_matches_full_ci_max_abs_error": proj_exact,
        "scan": scan,
        "preconditions": {
            "note": "the pinned description is only meaningful where the "
                    "natural occupations are nondegenerate AND the carrier is "
                    "the dominant support; points failing either are reported "
                    "but not used to support the bound",
            "usable_points": len(usable),
            "usable_genuinely_complex_flux_points": len(complex_usable),
            "total_points": len(scan),
        },
        "complex_flux_physical_benchmark": {
            "status": "NUMERICAL_WITH_EXACT_MODEL_INPUT",
            "selected_scan_point": selected_complex,
            "exact_model_input": {
                "particles": 3,
                "orbitals": 6,
                "hopping": "1",
                "nearest_neighbour_interaction": "2",
                "onsite_gradient": "7/20",
                "boundary_phase": "(4+3*i)/5",
                "boundary_phase_modulus_squared": "1",
                "time_reversal_is_broken": True,
            },
            "gates": {
                "nondegenerate_ground_state": selected_complex["ground_state_gap"] > 1e-8,
                "nondegenerate_natural_occupations": selected_complex[
                    "natural_occupations_nondegenerate"
                ],
                "natural_basis_reconstruction": selected_complex[
                    "natural_basis_offdiagonal_max"
                ] < 1e-10,
                "dominant_pinned_carrier": selected_complex[
                    "carrier_is_dominant_support"
                ],
                "nonzero_gauge_invariant_carrier_flux": selected_complex[
                    "genuinely_complex_carrier_flux"
                ],
                "quasipinning_enlargement_contains_true_energy": selected_complex[
                    "true_energy_inside_enlarged"
                ],
            },
            "independent_verifier": (
                "scripts/verify_complex_flux_physical_benchmark_standalone.py"
            ),
        },
        "slack_versus_off_carrier_weight": {
            "observed_ratio_range": [
                min(r["slack_over_off_carrier_weight"] for r in usable),
                max(r["slack_over_off_carrier_weight"] for r in usable),
            ] if usable else None,
            "comment": "D(lambda) / eps^2 stays of order one across the usable "
                       "scan, so eps^2 is controlled linearly by the slack; "
                       "this is measured, not proved",
        },
        "enlargement_bound_holds_on_usable_points": all(
            r["true_energy_inside_enlarged"] for r in usable),
    }


# --------------------------------------------------------------------------
# block 7
# --------------------------------------------------------------------------

def block_atlas() -> dict:
    rows = [json.loads(line) for line
            in (ROOT / "results/data/states.jsonl").read_text().splitlines()
            if line.strip()]
    tiers: dict[str, int] = {}
    records = []
    for r in rows:
        dets = [tuple(s[0]) for s in r["support"]]
        m = len(dets)
        diffs = [len(set(a) - set(b)) for a, b in combinations(dets, 2)]
        singles = sum(1 for d in diffs if d == 1)
        doubles = sum(1 for d in diffs if d == 2)
        higher = sum(1 for d in diffs if d >= 3)
        accessible = singles + doubles
        exact = m <= 3
        tier = "exact_elliptope_sdp" if exact else "elliptope_outer_bound"
        tiers[tier] = tiers.get(tier, 0) + 1
        records.append({
            "system": r["system"], "index": r["index"], "class": r["classified"],
            "carrier_size": m,
            "phase_dimension": max(0, m - 1),
            "coherence_pairs": len(diffs),
            "single_excitation_pairs": singles,
            "double_excitation_pairs": doubles,
            "inaccessible_pairs": higher,
            "coherence_graph_complete": higher == 0,
            "two_body_accessible_fraction": (accessible / len(diffs)) if diffs else 1.0,
            "zero_width_two_body_dimension": m + 2 * higher,
            "exactly_spectrahedral": exact,
            "method": tier,
        })
    complete = sum(1 for x in records if x["coherence_graph_complete"])
    return {
        "n_records": len(records),
        "tiers": tiers,
        "carrier_size_histogram": {
            str(k): sum(1 for x in records if x["carrier_size"] == k)
            for k in sorted({x["carrier_size"] for x in records})},
        "coherence_graph_complete": complete,
        "hierarchy": {
            "m <= 3": "exact elliptope SDP (proved: extreme rank r^2 <= m forces "
                      "r = 1)",
            "m >= 4": "elliptope outer bound; exactness needs a separate proof "
                      "per support, and the order-four counterexample shows it "
                      "can genuinely fail",
        },
        "records": records,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", default="results/data/conditional_2rdm_program.json")
    args = ap.parse_args()

    print("block 4: flux and current ...", flush=True)
    b4 = block_flux_current()
    print(f"  exact current range {b4['exact_current_range']}, "
          f"{len(b4['branch_exchange_fluxes'])} branch-exchange fluxes")

    print("block 5/6: projection and quasipinning ...", flush=True)
    b56 = block_projection_and_quasipinning()
    print(f"  FCIDUMP round trip {b56['fcidump_round_trip']}, "
          f"projection error {b56['carrier_projection_matches_full_ci_max_abs_error']:.2e}")
    print(f"  usable scan points {b56['preconditions']['usable_points']}"
          f"/{b56['preconditions']['total_points']}, "
          f"enlargement bound holds: {b56['enlargement_bound_holds_on_usable_points']}")

    print("block 7: census atlas ...", flush=True)
    b7 = block_atlas()
    print(f"  {b7['n_records']} supports, tiers {b7['tiers']}")

    res = {
        "scope": "Steps 4-7 of the conditional 2-RDM program: flux/current joint "
                 "body, Slater-Condon projection, quasipinning enlargement, and "
                 "the census-wide carrier atlas.",
        "flux_and_current": b4,
        "projection_and_quasipinning": b56,
        "atlas": b7,
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
