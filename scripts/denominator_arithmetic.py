#!/usr/bin/env python3
"""Where a state denominator comes from: incidence solve, then pinning solve.

The state denominator of a certified extremal state is not the spectrum
denominator (that law was falsified; RESEARCH.md, P3). This script derives it
in two stages and reports which stage mints each prime.

Stage 1, INCIDENCE. Solve A_S p = lambda, sum p = 1 exactly over Q, where A_S
is the mode-incidence matrix of the support. The solution is an affine space
whose denominators are controlled by the invariant factors of A_S: when the
solve is unimodular (all nonzero invariant factors 1) the affine space is
defined over the spectrum denominator and introduces no new prime.

Stage 2, PINNING. The one-hop channel conditions cut the affine space down.
Each channel (A,B) of class size m contributes

    |rho_AB|^2 = sum_a p_{u_a} p_{v_a}
                 + 2 sum_{a<b} s_a s_b sqrt(p_{u_a} p_{v_a} p_{u_b} p_{v_b})
                               cos(Phi_ab)

(docs/holonomy_rationality.md, Lemma 2). A new prime enters the state
denominator through the pivot of this solve, not through the incidence step.

The worked example is v89, where the mechanism is visible in three lines: the
incidence solution is a line over the spectrum denominator 26, the silence
condition restricted to that line LOSES ITS QUADRATIC TERMS, and the surviving
linear equation 130 t - 1 = 0 mints the prime 5.

Usage:
  python scripts/denominator_arithmetic.py                # write the artifact
  python scripts/denominator_arithmetic.py --v89          # worked example only
  python scripts/denominator_arithmetic.py --v103
  python scripts/denominator_arithmetic.py --check
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import invariant_factors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.states import one_hop_classes  # noqa: E402

LEDGER = ROOT / "results" / "data" / "states.jsonl"
OUT = ROOT / "results" / "data" / "denominator_arithmetic.json"


def load_states(path: Path = LEDGER):
    return [json.loads(ln) for ln in path.read_text().splitlines() if ln.strip()]


def get(system: str, index: int):
    for r in load_states():
        if r["system"] == system and r["index"] == index:
            return r
    raise SystemExit(f"{system} idx {index} not in the ledger")


# --------------------------------------------------------------------------
# stage 1: the incidence solve


def incidence_data(rec):
    """Smith invariant factors, kernel dimension and solution space of A_S."""
    cf = rec["closed_form"]
    dets = [tuple(t) for t in cf["support_dets"]]
    n = len(dets)
    den_spec = rec["denominator"]
    lam = [sp.Rational(v, den_spec) for v in rec["integer_form"]]
    # d is the full orbital count from the spectrum, NOT max(mode in support)+1:
    # a mode the support never touches still carries the constraint lambda_m = 0,
    # and dropping those rows both loses constraints and misaligns the
    # normalization row against the right-hand side.
    d = len(lam)
    if any(m >= d for T in dets for m in T):
        raise SystemExit(f"support mode out of range for {rec['system']} v{rec['index']}")

    rows = [[1 if m in T else 0 for T in dets] for m in range(d)]
    rows.append([1] * n)
    A = sp.Matrix(rows)

    factors = [int(f) for f in invariant_factors(A) if f != 0]
    unimodular = all(f == 1 for f in factors)
    kernel_dim = n - A.rank()

    p = sp.symbols(f"p0:{n}")
    eqs = [sp.Eq(sum(A[r, t] * p[t] for t in range(n)), (lam + [sp.Integer(1)])[r])
           for r in range(A.rows)]
    sol = sp.solve(eqs, p, dict=True)
    sol = sol[0] if sol else {}
    free = [v for v in p if v not in sol]

    # denominators appearing in the particular solution
    dens = set()
    for v in sol.values():
        expr = v.subs({f: 0 for f in free})
        dens.add(int(sp.Rational(expr).q))
    solve_den = int(sp.ilcm(1, *dens)) if dens else 1

    return {
        "support": [list(t) for t in dets],
        "support_size": n,
        "spectrum_denominator": den_spec,
        "invariant_factors": factors,
        "unimodular": bool(unimodular),
        "incidence_kernel_dim": int(kernel_dim),
        "solution_denominator": solve_den,
        "incidence_stage_ratio": str(sp.Rational(solve_den, den_spec)),
        "incidence_stage_mints_nothing": bool(den_spec % solve_den == 0),
        "_symbols": p,
        "_solution": sol,
        "_free": free,
        "_dets": dets,
    }


# --------------------------------------------------------------------------
# stage 2: the pinning solve, worked exactly


def v89_worked_example():
    """The v89 minting derivation, exact and self-contained.

    The support carries exactly one one-hop class, of size 2, and the certified
    state is real, so the silence condition |rho_37|^2 = 0 reduces to
    p_a p_b = p_c p_d on the pair products. Restricted to the incidence line
    the two quadratic terms cancel identically and a linear equation survives.
    """
    rec = get("(3,10)", 89)
    inc = incidence_data(rec)
    p, sol, free, dets = inc["_symbols"], inc["_solution"], inc["_free"], inc["_dets"]
    if len(free) != 1:
        raise SystemExit(f"expected a one-parameter incidence line, got {len(free)}")
    t = free[0]

    classes = one_hop_classes(dets)
    if len(classes) != 1:
        raise SystemExit(f"expected one channel at v89, got {len(classes)}")
    (A_, B_), pairs = next(iter(classes.items()))
    if len(pairs) != 2:
        raise SystemExit(f"expected a size-2 class, got {len(pairs)}")

    (i1, j1), (i2, j2) = pairs
    on_line = {**sol, t: t}

    def pval(k):
        return sp.together(on_line.get(p[k], p[k]))

    # silence of a size-2 class at a real state: the pair products agree
    lhs = sp.expand(pval(i1) * pval(j1))
    rhs = sp.expand(pval(i2) * pval(j2))
    silence = sp.expand(lhs - rhs)
    poly = sp.Poly(silence, t)

    quadratic_cancels = poly.degree() <= 1
    root = sp.solve(sp.Eq(silence, 0), t)
    if len(root) != 1:
        raise SystemExit(f"expected a unique root, got {root}")
    t0 = sp.Rational(root[0])

    weights = {k: sp.Rational(sp.together(on_line.get(p[k], p[k])).subs(t, t0))
               for k in range(len(dets))}
    den_state = int(sp.ilcm(1, *[w.q for w in weights.values()]))

    cf = rec["closed_form"]
    ledger_weights = [sp.Rational(w, cf["den"]) for w in cf["weights"]]
    agrees = all(weights[k] == ledger_weights[k] for k in range(len(dets)))

    numer = sp.factor(sp.numer(sp.together(silence)))
    return {
        "vertex": "(3,10) v89",
        "spectrum": f"{rec['integer_form']} / {rec['denominator']}",
        "channel": [A_, B_],
        "channel_class_size": len(pairs),
        "incidence_solution_is_a_line": True,
        "free_parameter": str(t),
        "incidence_solution": {str(p[k]): str(sp.together(on_line.get(p[k], p[k])))
                               for k in range(len(dets))},
        "silence_condition_on_the_line": str(sp.factor(silence)),
        "quadratic_terms_cancel": bool(quadratic_cancels),
        "surviving_equation": f"{numer} = 0",
        "pinned_value": str(t0),
        "pinning_pivot": str(sp.nsimplify(1 / t0)),
        "state_denominator": den_state,
        "state_denominator_factored": {str(k): int(v)
                                       for k, v in sp.factorint(den_state).items()},
        "spectrum_denominator_factored": {str(k): int(v)
                                          for k, v in sp.factorint(rec["denominator"]).items()},
        "ratio": int(den_state // rec["denominator"]),
        "minted_primes": sorted(
            set(sp.factorint(den_state)) - set(sp.factorint(rec["denominator"]))
        ),
        "agrees_with_ledger": bool(agrees),
        "reading": (
            "The incidence solve is a line over the spectrum denominator 26 and "
            "mints nothing. The single silence condition, restricted to that "
            "line, loses its quadratic terms: both sides carry the same -t^2, "
            "so what survives is linear. Its pivot is 130, and the prime 5 in "
            "130 = 2 * 5 * 13 is minted there and nowhere else. No discriminant "
            "is involved at v89, which is why the state stays rational and real."
        ),
    }


def v103_worked_example():
    """The same derivation at v103, where the incidence solve leaves five free
    parameters and five channels do the pinning.

    v103's certified state is real up to gauge (every loop holonomy is an
    integer multiple of pi; results/data/holonomy_fields.json), so each channel
    condition |rho_AB|^2 = 0 reduces to a signed sum of square roots of pair
    products vanishing. The signs are read off the certified state, and the
    system is then solved exactly on the incidence affine space.
    """
    rec = get("(3,10)", 103)
    inc = incidence_data(rec)
    p, sol, free, dets = inc["_symbols"], inc["_solution"], inc["_free"], inc["_dets"]
    cf = rec["closed_form"]
    den = cf["den"]
    exact = [sp.Rational(w, den) for w in cf["weights"]]
    amps = [sp.sympify(s) for s in cf["pretty"]]

    on_space = {p[k]: sp.together(sol.get(p[k], p[k])) for k in range(len(dets))}

    equations = []
    for (A_, B_), pairs in sorted(one_hop_classes(dets).items()):
        oriented = []
        for i, j in pairs:
            u, v = (i, j) if A_ in dets[i] else (j, i)
            s = ((-1) ** dets[v].index(B_)) * ((-1) ** dets[u].index(A_))
            # sign of the real amplitude product, giving the term's sign
            term_sign = s * int(sp.sign(sp.re(sp.conjugate(amps[u]) * amps[v])))
            oriented.append((u, v, term_sign))
        # sum_a eta_a sqrt(p_u p_v) = 0, with eta from the certified state
        expr = sum(eta * sp.sqrt(on_space[p[u]] * on_space[p[v]])
                   for u, v, eta in oriented)
        at_certified = {f: exact[list(p).index(f)] for f in free}
        residual = sp.N(expr.subs(at_certified), 40)
        equations.append(
            {
                "channel": [A_, B_],
                "class_size": len(pairs),
                "oriented_pairs": [[u, v, eta] for u, v, eta in oriented],
                "vanishes_at_the_certified_state": bool(abs(residual) < sp.Float(10) ** -30),
                "residual": str(residual),
            }
        )

    # Solve the system exactly on the incidence space, by clearing the roots
    # pairwise. sympy handles the five-equation system directly in the squared
    # (product) form for the size-2 classes, which is what pins the primes.
    pinned = {}
    for eq in equations:
        if eq["class_size"] != 2:
            continue
        (u1, v1, e1), (u2, v2, e2) = eq["oriented_pairs"]
        # eta_1 sqrt(P1) + eta_2 sqrt(P2) = 0 with opposite signs means P1 = P2
        pinned[tuple(eq["channel"])] = sp.expand(
            on_space[p[u1]] * on_space[p[v1]] - on_space[p[u2]] * on_space[p[v2]]
        )

    solved = sp.solve([sp.Eq(v, 0) for v in pinned.values()], free, dict=True)

    # The same cancellation as at v89, one level up: two of the size-2 channel
    # equations carry identical quadratic parts, so their DIFFERENCE is linear
    # and it is there that the new prime appears as a literal coefficient.
    cancellations = []
    keys = sorted(pinned)
    for a in range(len(keys)):
        for b in range(a + 1, len(keys)):
            diff = sp.expand(sp.numer(sp.together(pinned[keys[a]] - pinned[keys[b]])))
            if diff == 0:
                continue
            poly = sp.Poly(diff, *free)
            if poly.total_degree() > 1:
                continue
            coeffs = [abs(int(c)) for c in poly.coeffs()]
            new_primes = sorted(
                {int(q) for c in coeffs for q in sp.factorint(c)}
                - set(sp.factorint(rec["denominator"]))
            )
            cancellations.append(
                {
                    "channels": [list(keys[a]), list(keys[b])],
                    "quadratic_parts_cancel": True,
                    "surviving_linear_equation": f"{sp.factor(diff)} = 0",
                    "coefficient_primes_beyond_the_spectrum_denominator": new_primes,
                }
            )

    ratio = sp.Rational(den, rec["denominator"])
    return {
        "vertex": "(3,10) v103",
        "spectrum": f"{rec['integer_form']} / {rec['denominator']}",
        "incidence_kernel_dim": inc["incidence_kernel_dim"],
        "incidence_unimodular": inc["unimodular"],
        "incidence_invariant_factors": inc["invariant_factors"],
        "incidence_solution_denominator": inc["solution_denominator"],
        "free_parameters": [str(f) for f in free],
        "channels": equations,
        "all_channels_silent_at_certified_state": all(
            e["vanishes_at_the_certified_state"] for e in equations
        ),
        "size_two_product_equations": {str(k): str(sp.factor(v)) for k, v in pinned.items()},
        "size_two_solution_branches": len(solved),
        "pairwise_quadratic_cancellations": cancellations,
        "state_denominator": den,
        "state_denominator_factored": {str(k): int(v) for k, v in sp.factorint(den).items()},
        "spectrum_denominator_factored": {
            str(k): int(v) for k, v in sp.factorint(rec["denominator"]).items()
        },
        "ratio": int(ratio),
        "ratio_factored": {str(k): int(v) for k, v in sp.factorint(ratio).items()},
        "minted_primes": sorted(
            set(sp.factorint(den)) - set(sp.factorint(rec["denominator"]))
        ),
        "reading": (
            "The incidence solve at v103 is unimodular and lands over the "
            "spectrum denominator 34, exactly as at v89, so it mints nothing. "
            "Everything beyond 34 comes from the pinning stage: the state "
            "denominator 195364 = 2^2 * 13^2 * 17^2 exceeds 34 = 2 * 17 by the "
            "ratio 5746 = 2 * 13^2 * 17, and the prime 13 appears in neither "
            "the spectrum nor the incidence solve. Unlike v89 the system is not "
            "a single linear equation: five channels act on a five-dimensional "
            "incidence kernel, three of the classes have size 2 and give "
            "product equations directly, while the size-3 and size-4 classes "
            "carry genuine square-root sums. 13 is minted by the pinning stage; "
            "attributing it to one individual pivot is not meaningful here "
            "because the five conditions are solved jointly."
        ),
        "prime_53_at_the_support_10_endpoint": (
            "NOT DERIVED. The support-10 endpoint of the v103 cascade has state "
            "denominator 46852 = 2^2 * 13 * 17 * 53 and its 1-RDM is not "
            "diagonal (off-diagonal 7.4e-2), so it is outside the cancellation "
            "regime and its channel conditions are not silence conditions. "
            "Deriving the 53 needs the diagonal support-10 recognition tracked "
            "as HANDOFF 2 item 4, which has not landed in this repository."
        ),
    }


# --------------------------------------------------------------------------
# the pre-registered census-wide test


def census_test():
    """Score P-DEN-2, P-DEN-3, P-DEN-4 of docs/prereg_denominator_arithmetic.md."""
    rows = []
    for rec in load_states():
        cf = rec["closed_form"]
        if not all(isinstance(w, int) for w in cf["weights"]):
            rows.append(
                {
                    "system": rec["system"],
                    "index": rec["index"],
                    "rational_weights": False,
                    "note": "excluded from the ratio predictions: no integer weights",
                }
            )
            continue
        inc = incidence_data(rec)
        den_state = cf["den"]
        den_spec = rec["denominator"]
        ratio = sp.Rational(den_state, den_spec)
        rows.append(
            {
                "system": rec["system"],
                "index": rec["index"],
                "classified": rec["classified"],
                "rational_weights": True,
                "spectrum_denominator": den_spec,
                "state_denominator": den_state,
                "ratio": str(ratio),
                "ratio_is_one": bool(ratio == 1),
                "unimodular": inc["unimodular"],
                "invariant_factors": inc["invariant_factors"],
                "incidence_kernel_dim": inc["incidence_kernel_dim"],
                "pinned": bool(inc["incidence_kernel_dim"] > 0),
            }
        )

    live = [r for r in rows if r.get("rational_weights")]
    p2_ex = [r for r in live if r["ratio_is_one"] and not r["unimodular"]]
    p3_ex = [r for r in live if r["unimodular"] and not r["ratio_is_one"]]
    unpinned = [r for r in live if not r["pinned"]]
    p4_ex = [r for r in unpinned if r["unimodular"] and not r["ratio_is_one"]]

    return {
        "prereg": "docs/prereg_denominator_arithmetic.md",
        "states_scored": len(live),
        "states_excluded_no_rational_weights": len(rows) - len(live),
        "excluded": [
            f"{r['system']} v{r['index']}" for r in rows if not r.get("rational_weights")
        ],
        "ratio_one_count": sum(1 for r in live if r["ratio_is_one"]),
        "unimodular_count": sum(1 for r in live if r["unimodular"]),
        "pinned_count": sum(1 for r in live if r["pinned"]),
        "P_DEN_2": {
            "claim": "ratio = 1 implies the incidence solve is unimodular",
            "verdict": "PASS" if not p2_ex else "FAIL",
            "exceptions": len(p2_ex),
            "exception_list": [f"{r['system']} v{r['index']}" for r in p2_ex[:20]],
        },
        "P_DEN_3": {
            "claim": "unimodular incidence solve implies ratio = 1",
            "verdict": "PASS" if not p3_ex else "FAIL",
            "exceptions": len(p3_ex),
            "exception_list": [f"{r['system']} v{r['index']}" for r in p3_ex[:20]],
            "all_exceptions_are_pinned": bool(all(r["pinned"] for r in p3_ex)),
        },
        "P_DEN_4": {
            "claim": (
                "among states with zero incidence kernel, unimodular implies ratio = 1"
            ),
            "verdict": "PASS" if not p4_ex else "FAIL",
            "population": len(unpinned),
            "exceptions": len(p4_ex),
            "exception_list": [f"{r['system']} v{r['index']}" for r in p4_ex[:20]],
        },
        "P_DEN_2_failure_mechanism": (
            "The rationale registered for P-DEN-2 was wrong. Unimodularity of "
            "A_S is SUFFICIENT for the solve to stay over the spectrum "
            "denominator but not NECESSARY: an invariant factor f > 1 obstructs "
            "only those right-hand sides outside the image sublattice, and the "
            "particular right-hand side (lambda, 1) can lie inside it anyway. "
            "The 18 exceptions are exactly states of that kind, with a "
            "non-unimodular A_S whose spectrum still happens to be reachable "
            "over den_spec."
        ),
        "P_DEN_3_failure_mechanism": (
            "Registered as an expected FAIL, and it failed as expected with the "
            "predicted shape: the single exception, v103, is pinned "
            "(incidence kernel 5), so the escape happened at the pinning stage, "
            "which unimodularity of the incidence step says nothing about."
        ),
        "ratio_not_one_split": {
            "status": "POST-HOC observation, not pre-registered, no out-of-sample test",
            "count": len(live) - sum(1 for r in live if r["ratio_is_one"]),
            "incidence_stage": {
                "description": (
                    "zero incidence kernel, so nothing is left to pin: the unique "
                    "incidence solution itself sits above the spectrum denominator"
                ),
                "states": [
                    f"{r['system']} v{r['index']}"
                    for r in live
                    if not r["ratio_is_one"] and not r["pinned"]
                ],
                "all_design_real": all(
                    r["classified"] == "DESIGN-REAL"
                    for r in live
                    if not r["ratio_is_one"] and not r["pinned"]
                ),
                "ratios": sorted(
                    {r["ratio"] for r in live if not r["ratio_is_one"] and not r["pinned"]}
                ),
            },
            "pinning_stage": {
                "description": (
                    "positive incidence kernel with the incidence stage minting "
                    "nothing: the whole ratio is made by the channel solve"
                ),
                "states": [
                    f"{r['system']} v{r['index']}"
                    for r in live
                    if not r["ratio_is_one"] and r["pinned"]
                ],
                "all_interference": all(
                    r["classified"] == "INTERFERENCE"
                    for r in live
                    if not r["ratio_is_one"] and r["pinned"]
                ),
                "ratios": sorted(
                    {r["ratio"] for r in live if not r["ratio_is_one"] and r["pinned"]}
                ),
            },
        },
        "records": live,
    }


def build():
    return {
        "generated_by": "scripts/denominator_arithmetic.py",
        "what": (
            "Two-stage derivation of the state denominator: the incidence solve "
            "and then the pinning solve, with the pre-registered census-wide "
            "test of which stage controls the ratio."
        ),
        "v89": v89_worked_example(),
        "v103": v103_worked_example(),
        "census": census_test(),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--v89", action="store_true")
    ap.add_argument("--v103", action="store_true")
    ap.add_argument("--census", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    if a.v89:
        print(json.dumps(v89_worked_example(), indent=2))
        return 0
    if a.v103:
        print(json.dumps(v103_worked_example(), indent=2))
        return 0
    if a.census:
        c = census_test()
        c.pop("records")
        print(json.dumps(c, indent=2))
        return 0

    doc = build()
    text = json.dumps(doc, indent=1) + "\n"
    if a.check:
        cur = Path(a.out).read_text() if Path(a.out).exists() else ""
        if cur != text:
            print(f"STALE: {a.out}", file=sys.stderr)
            return 1
        print(f"{a.out} is current")
        return 0

    Path(a.out).write_text(text)
    c = doc["census"]
    print(f"wrote {a.out}")
    print(f"  v89  ratio {doc['v89']['ratio']}  minted {doc['v89']['minted_primes']}"
          f"  quadratic cancels: {doc['v89']['quadratic_terms_cancel']}")
    print(f"  v103 ratio {doc['v103']['ratio']}  minted {doc['v103']['minted_primes']}")
    print(f"  census: {c['states_scored']} scored, "
          f"{c['states_excluded_no_rational_weights']} excluded {c['excluded']}")
    for key in ("P_DEN_2", "P_DEN_3", "P_DEN_4"):
        print(f"  {key}: {c[key]['verdict']}  ({c[key]['exceptions']} exceptions)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
