#!/usr/bin/env python3
"""The conditional 2-RDM compatibility body over the pinned Borland-Dennis fiber.

Promotes the interval solver to its actual object. The endpoints of a two-body
observable's range are the SUPPORT FUNCTION of

    K(gamma) = conv { Gamma^(2)_Psi : Psi |-> gamma },

so computing them for every H determines K(gamma). On this fiber the variable
part of K(gamma) is an affine image of the scaled complex elliptope
E_w = {rho >= 0, diag rho = w}, and at order three every extreme point of E_w
is a phase matrix, so the support function is an exact SDP rather than a
relaxation.

Emits results/data/conditional_2rdm_body.json:

- the order-three theorem with its extreme-rank bound,
- exact primal-dual KKT certificates for the zero-flux, pi-flux, quarter-flux
  and Pythagorean-flux carriers, each a dual point y with Diag(y) - H PSD
  proved from principal-minor signs at an isolated root,
- agreement of the SDP with the phase hull at order three,
- the order-four boundary: an explicit rank-two extreme point and a strictly
  positive relaxation gap, which is where the theorem stops.

Usage:
    python scripts/conditional_2rdm_body.py [-o results/data/conditional_2rdm_body.json]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import sympy as sp  # noqa: E402

from gpc_census import elliptope as el  # noqa: E402

W = [Fraction(3, 5), Fraction(1, 4), Fraction(3, 20)]
I = sp.I

CARRIERS = {
    "zero_flux": [[0, 1, 1], [1, 0, 1], [1, 1, 0]],
    "pi_flux": [[0, 1, 1], [1, 0, -1], [1, -1, 0]],
    "quarter_flux": [[0, 1, 1], [1, 0, I], [1, -I, 0]],
    "pythagorean_flux_cos_3_5_sin_4_5":
        [[0, 1, 1],
         [1, 0, sp.Rational(3, 5) + I * sp.Rational(4, 5)],
         [1, sp.Rational(3, 5) - I * sp.Rational(4, 5), 0]],
}


def _num(h):
    return [[complex(sp.N(sp.re(x), 30), sp.N(sp.im(x), 30)) for x in row] for row in h]


def build() -> dict:
    import numpy as np

    out: dict = {
        "scope": "Support function of the conditional 2-RDM body over the pinned "
                 "Borland-Dennis fiber, via the scaled complex elliptope.",
        "carrier_weights": [str(x) for x in W],
        "theorem": {
            "body": "K(gamma) variable part is an affine image of "
                    "E_w = {rho >= 0, diag(rho) = w}",
            "support_function": "F_H^+ = max{Tr(H rho) : rho in E_w} "
                                "= min{w.y : Diag(y) - H >= 0}",
            "extreme_rank_bound": "r^2 <= m (Li and Tam)",
            "order_three": "m = 3 forces r = 1, so E_w = conv{c c^dagger}: the "
                           "SDP is exact, not a relaxation",
            "rank_bound_by_order": {str(m): el.extreme_rank_bound(m)
                                    for m in (2, 3, 4, 5, 9)},
        },
        "certificates": {},
        "sdp_versus_phase_hull": {},
    }

    for name, h in CARRIERS.items():
        hs = [[sp.sympify(x) for x in row] for row in h]
        chi = sp.simplify(hs[0][1] * hs[1][2] * sp.conjugate(hs[0][2]))
        cert = el.exact_certificates(hs, W)
        mx, mn = cert.get("maximum"), cert.get("minimum")
        if mx is None or mn is None:
            raise AssertionError(f"no dual certificate found for {name}")
        rec = {
            "matrix": [[str(x) for x in row] for row in h],
            "cycle_product": {"re": str(sp.re(chi)), "im": str(sp.im(chi))},
            "degenerate_flux": bool(sp.im(chi) == 0),
            "objective_polynomial_ascending": cert["objective_polynomial_ascending"],
            "dual_y_offsets_ascending": cert["dual_y_offsets_ascending"],
            "maximum": mx,
            "minimum": mn,
            "all_real_roots": cert["roots"],
            "hidden_width_approx": float(mx["approx"] - mn["approx"]),
        }
        out["certificates"][name] = rec

        hn = _num(hs)
        wn = [float(x) for x in W]
        sdp_hi = float(el.support_numeric(hn, wn, "max"))
        sdp_lo = float(el.support_numeric(hn, wn, "min"))
        hull_hi = float(el.phase_hull_extreme_numeric(hn, wn, "max"))
        hull_lo = float(el.phase_hull_extreme_numeric(hn, wn, "min"))
        gap = float(max(abs(sdp_hi - hull_hi), abs(sdp_lo - hull_lo)))
        out["sdp_versus_phase_hull"][name] = {
            "sdp": [sdp_lo, sdp_hi],
            "phase_hull": [hull_lo, hull_hi],
            "gap": gap,
            "exact_at_order_three": bool(gap < 1e-6),
        }

    # order-four boundary
    vs, R = el.rank_two_extreme_point_4x4()
    V = np.array(vs).T
    P = V.conj().T @ np.linalg.pinv(V.conj().T)
    M = np.eye(4) - P
    H4 = -M
    w4 = [1.0] * 4
    sdp4 = float(el.support_numeric(H4, w4, "max"))
    hull4 = float(el.phase_hull_extreme_numeric(H4, w4, "max"))
    out["order_four_boundary"] = {
        "why": "r^2 <= m permits r = 2 at m = 4, so the elliptope acquires "
               "extreme points that are not phase matrices and the SDP becomes "
               "an outer bound on the conditional body",
        "extreme_point_rank": int(np.linalg.matrix_rank(R, tol=1e-9)),
        "is_extreme": bool(el.is_extreme_point(vs)),
        "diagonal": [float(np.real(R[i, i])) for i in range(4)],
        "witness_hamiltonian": "minus the projector onto the orthogonal "
                               "complement of the extreme point's range",
        "sdp_max": sdp4,
        "phase_hull_max": hull4,
        "relaxation_gap": float(sdp4 - hull4),
        "sdp_is_strict_outer_bound": bool(sdp4 - hull4 > 1e-6),
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", default="results/data/conditional_2rdm_body.json")
    args = ap.parse_args()

    res = build()
    for name, rec in res["certificates"].items():
        lo, hi = rec["minimum"]["approx"], rec["maximum"]["approx"]
        flag = " (degenerate flux; the sextic branch engine excludes this case)" \
            if rec["degenerate_flux"] else ""
        print(f"{name:34s} [{lo: .12f}, {hi: .12f}]  width {hi - lo:.12f}{flag}")
        agree = res["sdp_versus_phase_hull"][name]["exact_at_order_three"]
        print(f"{'':34s} SDP == phase hull: {agree}")
    b = res["order_four_boundary"]
    print(f"\norder-four boundary: rank-{b['extreme_point_rank']} extreme point, "
          f"relaxation gap {b['relaxation_gap']:.10f}")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
