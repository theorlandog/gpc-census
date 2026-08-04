#!/usr/bin/env python3
"""Independent verifier for the fixed-1RDM observable-range artifact.

Nothing is replayed from the artifact except the target 1-RDM and the two
carrier matrices. Everything else is recomputed from first principles in exact
symbolic arithmetic and compared against the recorded values by symbolic
equality, not string equality, so a different but equivalent radical form
still passes.

Recomputed here:

* the Borland-Dennis pair sums, the trace, and the saturated GPC;
* the selection-rule carrier, as the intersection of the single-occupancy
  determinants with the kernel of the facet operator;
* the carrier weights, solved from the six diagonal occupations;
* one-hop-freeness of the carrier;
* the couplings A, B, C and the phasor identity behind them;
* both exact ranges by the branch rule, and the polygon range independently;
* realizability of every recorded phase triple, via the cosine Gram identity.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT = ROOT / "results/data/fixed_1rdm_observable_ranges.json"

SINGLE_OCCUPANCY_PAIRS = ((0, 5), (1, 4), (2, 3))
FACET_ORBITALS = frozenset({0, 1, 3})
FACET_LEVEL = 2


def selection_rule_carrier():
    """The carrier forced by the two exact selection rules, from scratch."""
    determinants = list(itertools.combinations(range(6), 3))
    single = [
        D for D in determinants
        if all(len(set(D) & set(pair)) == 1 for pair in SINGLE_OCCUPANCY_PAIRS)
    ]
    kernel = [
        D for D in determinants
        if FACET_LEVEL - sum(o in FACET_ORBITALS for o in D) == 0
    ]
    return determinants, single, kernel, sorted(set(single) & set(kernel))


def weights_from_diagonal(carrier, diagonal):
    """Solve the six occupation equations for the carrier weights."""
    symbols = sp.symbols(f"w0:{len(carrier)}", positive=True)
    equations = [
        sp.Eq(sum(symbols[k] for k, D in enumerate(carrier) if p in D), diagonal[p])
        for p in range(len(diagonal))
    ]
    solutions = sp.solve(equations, symbols, dict=True)
    if len(solutions) != 1:
        return None
    return [solutions[0][s] for s in symbols]


def one_hop_pairs(carrier):
    return [
        (a, b)
        for a, b in itertools.combinations(range(len(carrier)), 2)
        if len(set(carrier[a]) & set(carrier[b])) == len(carrier[a]) - 1
    ]


def couplings(weights):
    r = [sp.sqrt(w) for w in weights]
    return (
        sp.simplify(2 * r[0] * r[1]),
        sp.simplify(2 * r[0] * r[2]),
        sp.simplify(2 * r[1] * r[2]),
        r,
    )


def exact_range(a, b, c):
    """Both global extrema of E0 + a cos(al) + b cos(be) + c cos(be-al).

    At fixed t = cos(al) the beta optimum is a t +/- sqrt(b^2+c^2+2bc t); the
    plus branch is concave and the minus branch convex on [-1,1], so the four
    endpoints plus at most one feasible interior stationary value contain both
    extrema.
    """
    candidates = []
    for t in (1, -1):
        radical = sp.sqrt(sp.simplify(b**2 + c**2 + 2 * b * c * t))
        candidates.append((sp.simplify(a * t + radical), "upper", t))
        candidates.append((sp.simplify(a * t - radical), "lower", t))
    interior = None
    if a != 0 and b != 0 and c != 0:
        t_star = sp.simplify(((b * c / a) ** 2 - b**2 - c**2) / (2 * b * c))
        if abs(sp.N(t_star)) <= 1:
            interior = sp.simplify(a * t_star - b * c / a)
            candidates.append((interior, "interior", t_star))
    values = [v for v, _, _ in candidates]
    return sp.simplify(sp.Min(*values)), sp.simplify(sp.Max(*values)), candidates


def polygon_range(radii):
    total = sp.simplify(sum(radii))
    largest = max(radii, key=lambda x: sp.N(x))
    lower = sp.simplify(sp.Max(0, 2 * largest - total) ** 2 - 1)
    return lower, sp.simplify(total**2 - 1), sp.N(2 * largest - total) <= 0


def gram_defect(cos_alpha, cos_beta, cos_difference):
    """Zero exactly when the three cosines come from real angles on a circle."""
    return sp.simplify(
        1 - cos_alpha**2 - cos_beta**2 - cos_difference**2
        + 2 * cos_alpha * cos_beta * cos_difference
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    artifact = json.loads(args.artifact.read_text())

    checks: list[tuple[str, bool]] = []

    def report(name, ok, detail=""):
        checks.append((name, bool(ok)))
        print(f"  {'PASS' if ok else 'FAIL':4s} {name}{('  ' + detail) if detail else ''}")

    diagonal = [sp.Rational(x) for x in artifact["target_1rdm"]["diagonal"]]
    print("target 1-RDM")
    report("trace_is_the_particle_number", sum(diagonal) == 3)
    report("occupations_strictly_decreasing",
           all(diagonal[i] > diagonal[i + 1] for i in range(5)))
    report("borland_dennis_pair_sums_are_one",
           diagonal[0] + diagonal[5] == 1
           and diagonal[1] + diagonal[4] == 1
           and diagonal[2] + diagonal[3] == 1)
    report("gpc_is_saturated", 2 - (diagonal[0] + diagonal[1] + diagonal[3]) == 0)

    print("selection rules")
    determinants, single, kernel, carrier = selection_rule_carrier()
    rule = artifact["selection_rule_certificate"]
    report("ambient_determinant_count", len(determinants) == rule["ambient_determinants"])
    report("single_occupancy_support", [list(D) for D in single] == rule["single_occupancy_support"])
    report("facet_kernel_support", [list(D) for D in kernel] == rule["facet_kernel_support"])
    report("intersection_is_the_certified_carrier",
           [list(D) for D in carrier] == rule["intersection_support"] == artifact["fiber"]["support_dets"])

    print("fiber")
    weights = weights_from_diagonal(carrier, diagonal)
    report("weights_uniquely_determined_by_the_diagonal", weights is not None)
    report("weights_match_the_artifact",
           [str(w) for w in weights] == artifact["fiber"]["weights"])
    report("weights_sum_to_one", sum(weights) == 1)
    report("carrier_is_one_hop_free", not one_hop_pairs(carrier))
    report("projective_fiber_is_a_two_torus",
           len(carrier) - 1 == artifact["fiber"]["projective_real_dimension"])

    print("couplings and phasor identity")
    A, B, C, radii = couplings(weights)
    alpha, beta = sp.symbols("alpha beta", real=True)
    real_part = radii[0] + radii[1] * sp.cos(alpha) + radii[2] * sp.cos(beta)
    imag_part = radii[1] * sp.sin(alpha) + radii[2] * sp.sin(beta)
    phasor = sp.expand_trig(sp.expand(real_part**2 + imag_part**2 - 1))
    normal = sp.expand_trig(sp.expand(
        A * sp.cos(alpha) + B * sp.cos(beta) + C * sp.cos(beta - alpha)
    ))
    report("squared_radii_sum_to_one", sp.simplify(sum(x**2 for x in radii)) == 1)
    report("phasor_form_equals_the_normal_form", sp.simplify(phasor - normal) == 0)

    print("exact ranges")
    for name in sorted(artifact["examples"]):
        block = artifact["examples"][name]["range"]
        matrix = [[sp.sympify(v) for v in row] for row in block["carrier_matrix"]]
        a, b, c = A * matrix[0][1], B * matrix[0][2], C * matrix[1][2]
        report(f"{name}.couplings_match",
               sp.simplify(a - sp.sympify(block["couplings"]["A_cos_alpha"])) == 0
               and sp.simplify(b - sp.sympify(block["couplings"]["B_cos_beta"])) == 0
               and sp.simplify(c - sp.sympify(block["couplings"]["C_cos_beta_minus_alpha"])) == 0)
        low, high, _ = exact_range(a, b, c)
        report(f"{name}.minimum",
               sp.simplify(low - sp.sympify(block["minimum"])) == 0, str(low))
        report(f"{name}.maximum",
               sp.simplify(high - sp.sympify(block["maximum"])) == 0, str(high))
        report(f"{name}.width",
               sp.simplify(sp.sympify(block["width"]) - (high - low)) == 0)
        report(f"{name}.numeric_endpoints_agree",
               abs(float(sp.N(low)) - block["minimum_numeric"]) < 1e-12
               and abs(float(sp.N(high)) - block["maximum_numeric"]) < 1e-12)
        for which in ("minimum_phase_cosines", "maximum_phase_cosines"):
            cosines = [sp.sympify(block[which][k])
                       for k in ("cos_alpha", "cos_beta", "cos_beta_minus_alpha")]
            report(f"{name}.{which}_are_realizable",
                   all(abs(sp.N(x)) <= 1 for x in cosines) and gram_defect(*cosines) == 0)

    print("independent polygon certificate")
    certificate = artifact["examples"]["unfrustrated_coherence_triangle"]["independent_polygon_certificate"]
    low_poly, high_poly, closes = polygon_range(radii)
    report("triangle_closes", closes and certificate["polygon_closes"] is True)
    report("polygon_minimum", sp.simplify(low_poly - sp.sympify(certificate["minimum"])) == 0)
    report("polygon_maximum", sp.simplify(high_poly - sp.sympify(certificate["maximum"])) == 0)
    unfrustrated = artifact["examples"]["unfrustrated_coherence_triangle"]["range"]
    report("polygon_agrees_with_the_branch_rule",
           sp.simplify(sp.sympify(unfrustrated["minimum"]) - low_poly) == 0
           and sp.simplify(sp.sympify(unfrustrated["maximum"]) - high_poly) == 0)

    print("pi flux")
    pi_flux = artifact["examples"]["pi_flux_coherence_triangle"]["range"]
    report("pi_flux_range_is_the_negated_unfrustrated_range",
           sp.simplify(sp.sympify(pi_flux["maximum"]) + sp.sympify(unfrustrated["minimum"])) == 0
           and sp.simplify(sp.sympify(pi_flux["minimum"]) + sp.sympify(unfrustrated["maximum"])) == 0)
    report("cycle_products_differ_in_sign",
           sp.sympify(unfrustrated["carrier_matrix"][0][1])
           * sp.sympify(unfrustrated["carrier_matrix"][0][2])
           * sp.sympify(unfrustrated["carrier_matrix"][1][2]) > 0
           and sp.sympify(pi_flux["carrier_matrix"][0][1])
           * sp.sympify(pi_flux["carrier_matrix"][0][2])
           * sp.sympify(pi_flux["carrier_matrix"][1][2]) < 0)

    failed = [name for name, ok in checks if not ok]
    print(f"\n{len(checks) - len(failed)} of {len(checks)} checks passed")
    if failed:
        print("failed: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
