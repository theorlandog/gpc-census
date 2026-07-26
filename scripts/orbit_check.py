#!/usr/bin/env python3
"""The 2-RDM gate: are the cascade endpoints NEW states, or the library states rotated?

The 2-RDM transforms by conjugation under a one-body unitary, so its spectrum is
invariant along the entire U(d) orbit of a state. That makes it a cheap
basis-free fingerprint and the decisive test for a claimed new vertex attainer:

    a state offered as a NEW attainer must differ from the library
    representative in 2-RDM spectrum; if it does not, it is the SAME state in
    different orbitals, which is a gauge statement and must be reported as one.

This script runs the gate on every sparsification-cascade endpoint, and for the
headline case (v103) goes further and solves for the connecting one-body
unitary explicitly, which upgrades "indistinguishable by this invariant" to
"here is the rotation".

It also records, per endpoint, whether that rotation could possibly be a
NATURAL-ORBITAL gauge rotation, by testing whether the endpoint's own 1-RDM is
diagonal with diagonal lambda. An endpoint that fails this test attains the
vertex spectrum in a NON-natural-orbital basis, so its support bounds the
orbit-minimal support of that state and NOT s_Q^NO.

  python scripts/orbit_check.py

Writes results/data/orbit_check.json.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

from gpc_census.fiber import amplitudes_from_record, one_rdm
from gpc_census.gauge import (
    degenerate_blocks,
    same_state_up_to_orbital_rotation,
    solve_orbital_rotation,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
DIAG_TOL = 1e-9


def main():
    states = {}
    for line in (ROOT / "results/data/states.jsonl").read_text().splitlines():
        r = json.loads(line)
        states[(r["system"], r["index"])] = r

    rows = []
    for line in (ROOT / "results/data/cascade.jsonl").read_text().splitlines():
        cas = json.loads(line)
        if cas["support_upper"] >= cas["support_certified"]:
            continue  # the cascade did not move this state
        key = (cas["system"], cas["index"])
        lib = states[key]
        c, dets = amplitudes_from_record(lib)
        lam = [v / lib["denominator"] for v in lib["integer_form"]]
        d = len(lam)
        ce = (np.array(cas["final_amplitudes_real"])
              + 1j * np.array(cas["final_amplitudes_imag"]))
        de = [tuple(t) for t in cas["final_support_dets"]]

        same, dev = same_state_up_to_orbital_rotation(c, dets, ce, de, d)
        rho = one_rdm(ce, de, d)
        diag_err = float(np.max(np.abs(rho - np.diag(lam))))
        offdiag = float(np.max(np.abs(rho - np.diag(np.diag(rho)))))
        rotation = None
        if same:
            # exhibit the rotation, not just the invariant agreeing: this is
            # what a reader needs to reproduce the sparser expansion
            u, resid = solve_orbital_rotation(c, dets, ce, de, d, blocks=None,
                                              tries=12)
            rotation = {
                "residual": resid,
                "verified": bool(resid < 1e-9),
                "unitarity_error": float(np.max(np.abs(
                    u.conj().T @ u - np.eye(d)))) if u is not None else None,
                "matrix_real": [[float(x) for x in row] for row in np.real(u)],
                "matrix_imag": [[float(x) for x in row] for row in np.imag(u)],
            }
        rows.append({
            "orbital_rotation": rotation,
            "system": cas["system"],
            "index": cas["index"],
            "support_library": cas["support_certified"],
            "support_endpoint": cas["support_upper"],
            "two_rdm_spectrum_deviation": dev,
            "same_state_as_library": bool(same),
            "endpoint_one_rdm_diagonal_error": diag_err,
            "endpoint_one_rdm_max_offdiagonal": offdiag,
            "endpoint_is_natural_orbital_basis": bool(diag_err < DIAG_TOL),
            "bounds": ("s_Q^NO (natural-orbital minimum)" if diag_err < DIAG_TOL
                       else "orbit-minimal support of the certified "
                            "representative, NOT s_Q^NO"),
        })

    v103 = next((r for r in rows if r["system"] == "(3,10)" and r["index"] == 103),
                None)
    explicit = None
    if v103 is not None:
        lib = states[("(3,10)", 103)]
        c, dets = amplitudes_from_record(lib)
        lam = [v / lib["denominator"] for v in lib["integer_form"]]
        cas = next(json.loads(line) for line
                   in (ROOT / "results/data/cascade.jsonl").read_text().splitlines()
                   if json.loads(line)["system"] == "(3,10)"
                   and json.loads(line)["index"] == 103)
        ce = (np.array(cas["final_amplitudes_real"])
              + 1j * np.array(cas["final_amplitudes_imag"]))
        de = [tuple(t) for t in cas["final_support_dets"]]
        u, resid = solve_orbital_rotation(c, dets, ce, de, 10, blocks=None, tries=40)
        blocks = degenerate_blocks(lam)
        off_block = float(max(np.max(np.abs(u[np.ix_(blocks[0], blocks[1])])),
                              np.max(np.abs(u[np.ix_(blocks[1], blocks[0])]))))
        explicit = {
            "claim": "the support-10 endpoint is the certified library state in "
                     "different orbitals, and the change of orbitals is NOT a "
                     "natural-orbital gauge rotation",
            "residual_of_Gamma_U_psi_lib_minus_psi_end": resid,
            "unitarity_error": float(np.max(np.abs(u.conj().T @ u - np.eye(10)))),
            "max_off_block_entry": off_block,
            "is_gauge_rotation": bool(off_block < 1e-9),
            "consequence": "bounds the orbit-minimal support of that state; does "
                           "NOT bound s_Q^NO, and must not be quoted next to s_I",
        }

    report = {
        "note": "2-RDM spectrum gate on every cascade endpoint. A matching "
                "spectrum plus a continuous path from the certified state is "
                "strong evidence of orbit membership, not proof; the v103 entry "
                "carries the explicit rotation, which is proof.",
        "endpoints_tested": len(rows),
        "same_state_as_library": sum(1 for r in rows if r["same_state_as_library"]),
        "new_states_found": sum(1 for r in rows if not r["same_state_as_library"]),
        "endpoints_in_natural_orbital_basis":
            sum(1 for r in rows if r["endpoint_is_natural_orbital_basis"]),
        "endpoints_not_in_natural_orbital_basis":
            sum(1 for r in rows if not r["endpoint_is_natural_orbital_basis"]),
        "worst_two_rdm_deviation": max((r["two_rdm_spectrum_deviation"]
                                        for r in rows), default=0.0),
        "v103_explicit_rotation": explicit,
        "endpoints": rows,
    }
    (ROOT / "results/data/orbit_check.json").write_text(
        json.dumps(report, indent=2) + "\n")
    for k, v in report.items():
        if k not in ("endpoints", "v103_explicit_rotation", "note"):
            print(f"{k}: {v}")
    if explicit:
        print("v103 explicit rotation:")
        for k, v in explicit.items():
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
