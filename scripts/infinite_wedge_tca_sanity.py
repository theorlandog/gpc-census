#!/usr/bin/env python3
"""Finite consequences of the stable highest-weight semigroup, in exact arithmetic.

Pre-registration: `docs/prereg_infinite_wedge_tca_sanity.md`, written and
committed before this script was run. Research note:
`docs/infinite_wedge_tca_bridge.md`. Artifact:
`results/data/infinite_wedge_tca_sanity.json`. Guard test:
`tests/test_infinite_wedge_tca_sanity.py`.

THE OBJECT. For a fixed particle number `N` write

    Sigma_{N,d} = { (lambda, k) : lambda a partition, ell(lambda) <= d, k >= 1,
                    mult(V_lambda^{GL_d}, Sym^k(wedge^N C^d)) > 0 }.

Because `Sym^k(wedge^N C^d)` decomposes with multiplicities equal to the
plethysm coefficients `<s_lambda, h_k[e_N]>`, which do not depend on `d`, and
because `S_lambda(C^d)` vanishes exactly when `ell(lambda) > d`, the whole
family is the length truncation of ONE `d`-independent monoid

    Sigma_N^stab = { (lambda, k) : <s_lambda, h_k[e_N]> > 0 },
    Sigma_{N,d}  = Sigma_N^stab intersect { ell(lambda) <= d }.

`{ell <= d}` is a face of the dominant cone, so the same holds after passing to
cones, and the normalized `k = 1` slice is the moment polytope. This script does
not prove that chain; it tests the six finite consequences that would break
first if any convention or functoriality assumption in it were wrong.

BLOCKS

A. convention   three routes to the same multiplicity at the 14 census vertices
                with a measured sequence, plus a control that the convention is
                load bearing.
B. recovery     conv of the normalized occurring highest weights equals the
                certified moment polytope at (3,6) and (3,7).
C. saturation   lattice points of cone(Sigma_{N,d}) that are NOT in Sigma_{N,d}.
D. padding      does a published row stay valid after zero padding up a rank?
E. restriction  does the rank d+1 system imply the rank d system after
                restricting to lambda_{d+1} = 0?
F. templates    padding templates and symmetric templates of the published
                rows, per system, and how many are new at each rank.

Everything is exact: Fractions, integer plethysm coefficients, and a rational
Phase-I simplex for the Farkas implications.

Run:
    uv run scripts/infinite_wedge_tca_sanity.py
    uv run scripts/infinite_wedge_tca_sanity.py --weyl-max-k 6 --quick
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gpc_census import exact_polytope as ep  # noqa: E402
from gpc_census.constraints import constraints  # noqa: E402
from gpc_census.plethysm import hm_en, inner_spectra, mn  # noqa: E402

DATA = ROOT / "results" / "data"
OUT = DATA / "infinite_wedge_tca_sanity.json"

LADDER_N3 = [(3, 6), (3, 7), (3, 8), (3, 9), (3, 10)]
LADDER_N4 = [(4, 8), (4, 9), (4, 10)]
ALL_SYSTEMS = LADDER_N3 + LADDER_N4 + [(5, 10)]

# Recovery degrees. (3,6) needs m up to the largest vertex denominator 4, and
# (3,7) needs m up to 7; the margins here are the ones the Stage-0 generator
# already records as sufficient.
RECOVERY = [(3, 6, 6), (3, 7, 10)]


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

def plethysm_multiplicity(lam: tuple[int, ...], k: int, n: int) -> int:
    """`<s_lam, h_k[e_n]>`, exact, and independent of the rank `d`."""
    lam = tuple(int(x) for x in lam if x > 0)
    if sum(lam) != n * k:
        return 0
    total = Fraction(0)
    for mu, coefficient in hm_en(k, n).items():
        character = mn(lam, tuple(sorted(mu, reverse=True)))
        if character:
            total += coefficient * character
    assert total.denominator == 1, "plethysm coefficient must be an integer"
    return int(total)


def census_vertices(n: int, d: int) -> list[tuple[Fraction, ...]]:
    rows = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
    out = []
    for rec in rows:
        den = int(rec["denominator"] or 1)
        out.append(tuple(Fraction(int(v), den) for v in rec["integer_form"]))
    return out


def published_rows(n: int, d: int) -> tuple[list[tuple[list[int], int]],
                                            list[tuple[list[int], int]]]:
    """`(inequalities, equalities)` as integer `c . lambda <= b` / `= b` rows."""
    system = constraints(n, d)
    ineqs = [([int(c) for c in r["coeffs"]], int(r["rhs"]))
             for r in system["inequalities"]]
    eqs = [([int(c) for c in r["coeffs"]], int(r["rhs"]))
           for r in system["equalities"]]
    return ineqs, eqs


def _normalize(coeffs: list[int], rhs: int) -> tuple[tuple[int, ...], int]:
    from math import gcd
    g = 0
    for c in coeffs:
        g = gcd(g, abs(int(c)))
    g = gcd(g, abs(int(rhs))) or 1
    return tuple(int(c) // g for c in coeffs), int(rhs) // g


def padding_template(coeffs: list[int], rhs: int) -> tuple[tuple[int, ...], int]:
    """Trailing zeros stripped, then divided by the positive gcd (OI notion)."""
    trimmed = list(coeffs)
    while trimmed and trimmed[-1] == 0:
        trimmed.pop()
    return _normalize(trimmed, rhs)


def symmetric_template(coeffs: list[int], rhs: int) -> tuple[tuple[int, ...], int]:
    """Multiset of nonzero coefficients, sorted, then gcd normalized (FI notion)."""
    nonzero = sorted(int(c) for c in coeffs if c)
    return _normalize(nonzero, rhs)


# ---------------------------------------------------------------------------
# exact Farkas implication
# ---------------------------------------------------------------------------

def implies(ineqs_ge: list[tuple[list[Fraction], Fraction]],
            eqs: list[tuple[list[Fraction], Fraction]],
            target: tuple[list[Fraction], Fraction]) -> bool:
    """Does ``{a.x >= c} and {e.x = f}`` imply ``g.x >= h``, over the rationals?

    Affine Farkas: the implication holds on a nonempty system exactly when

        g = sum y_i a_i + sum z_j e_j,     h = sum y_i c_i + sum z_j f_j - s

    is solvable with ``y >= 0``, ``s >= 0`` and ``z`` free. Free variables are
    split into a nonnegative difference, and the whole thing is handed to the
    project's exact Phase-I simplex.
    """
    g, h = target
    dim = len(g)
    columns: list[list[Fraction]] = []
    for a, c in ineqs_ge:                       # y_i >= 0
        columns.append([Fraction(v) for v in a] + [Fraction(c)])
    for e, f in eqs:                            # z_j = z+ - z-
        columns.append([Fraction(v) for v in e] + [Fraction(f)])
        columns.append([-Fraction(v) for v in e] + [-Fraction(f)])
    columns.append([Fraction(0)] * dim + [Fraction(-1)])   # slack s >= 0
    mat = [[col[r] for col in columns] for r in range(dim + 1)]
    rhs = [Fraction(v) for v in g] + [Fraction(h)]
    return ep.feasible_nonneg(mat, rhs)


def chamber_rows_ge(d: int) -> list[tuple[list[Fraction], Fraction]]:
    """`1 >= l_1 >= ... >= l_d >= 0` in `>=` orientation."""
    rows: list[tuple[list[Fraction], Fraction]] = [
        ([Fraction(-1)] + [Fraction(0)] * (d - 1), Fraction(-1))]
    for i in range(d - 1):
        a = [Fraction(0)] * d
        a[i], a[i + 1] = Fraction(1), Fraction(-1)
        rows.append((a, Fraction(0)))
    a = [Fraction(0)] * d
    a[d - 1] = Fraction(1)
    rows.append((a, Fraction(0)))
    return rows


# ---------------------------------------------------------------------------
# A. convention
# ---------------------------------------------------------------------------

def block_convention(weyl_max_k: int, weyl_budget: float) -> dict:
    stored = json.loads((DATA / "quantization_multiplicities.json").read_text())
    weyl_counters: dict[tuple[int, int], object] = {}
    rows, mismatches, control_hits = [], [], []
    weyl_done = 0
    weyl_started = time.time()

    for rec in stored["rows"]:
        n, d = (int(x) for x in rec["system"].strip("()").split(","))
        q = int(rec["denominator"] or 1)
        integer_form = [int(x) for x in rec["integer_form"]]
        for entry in rec["multiplicities"]:
            k = int(entry["k"])
            mu = tuple(x * k // q for x in integer_form)
            pleth = plethysm_multiplicity(mu, k, n)

            weyl = None
            if k <= weyl_max_k and time.time() - weyl_started < weyl_budget:
                if (n, d) not in weyl_counters:
                    from quantization_multiplicity import weight_counter
                    weyl_counters[(n, d)] = weight_counter(n, d)
                from quantization_multiplicity import multiplicity as weyl_multiplicity
                weyl = int(weyl_multiplicity(list(mu), n, d, weyl_counters[(n, d)]))
                weyl_done += 1

            # Controls: two standard wrong conventions on the same weight.
            conjugate = tuple(sum(1 for x in mu if x > j) for j in range(max(mu) or 1))
            reversed_complement = tuple(k - x for x in reversed(mu))
            control = {
                "conjugate": plethysm_multiplicity(conjugate, k, n),
                "reversed_complement": plethysm_multiplicity(reversed_complement, k, n),
            }
            if any(v != int(entry["multiplicity"]) for v in control.values()):
                control_hits.append(f"{rec['system']}#{rec['index']}@k={k}")

            row = {
                "system": rec["system"], "index": int(rec["index"]), "k": k,
                "weight": list(mu),
                "plethysm": pleth,
                "weyl_alternating_sum": weyl,
                "stored": int(entry["multiplicity"]),
                "control_wrong_conventions": control,
                "evidence": "exact",
                "provenance": "census vertex with a measured multiplicity sequence",
            }
            agree = pleth == int(entry["multiplicity"]) and (weyl is None or weyl == pleth)
            if not agree:
                mismatches.append(row)
            rows.append(row)

    return {
        "what": "three routes to mult(V_{k lambda}, Sym^k(wedge^N C^d))",
        "routes": ["Murnaghan-Nakayama plethysm pairing <s_lambda, h_k[e_N]>",
                   "Weyl alternating sum over exact weight counts",
                   "results/data/quantization_multiplicities.json"],
        "pairs": len(rows),
        "weyl_recomputed_pairs": weyl_done,
        "weyl_max_k": weyl_max_k,
        "mismatches": mismatches,
        "all_agree": not mismatches,
        "control_pairs_where_a_wrong_convention_differs": len(control_hits),
        "control_is_load_bearing": bool(control_hits),
        "rows": rows,
    }


# ---------------------------------------------------------------------------
# B. recovery
# ---------------------------------------------------------------------------

def block_recovery() -> dict:
    out = []
    for n, d, mmax in RECOVERY:
        started = time.time()
        points = inner_spectra(n, d, mmax)
        ineqs, eqs = published_rows(n, d)
        violations = []
        for spectrum in points:
            for coeffs, rhs in ineqs:
                if sum(Fraction(c) * x for c, x in zip(coeffs, spectrum)) > rhs:
                    violations.append({"spectrum": [str(x) for x in spectrum],
                                       "row": coeffs, "rhs": rhs})
                    break
            else:
                for coeffs, rhs in eqs:
                    if sum(Fraction(c) * x for c, x in zip(coeffs, spectrum)) != rhs:
                        violations.append({"spectrum": [str(x) for x in spectrum],
                                           "equality": coeffs, "rhs": rhs})
                        break
        verts = census_vertices(n, d)
        first_degree = {}
        missing = []
        for v in verts:
            hit = points.get(v)
            if hit is None:
                missing.append([str(x) for x in v])
            else:
                first_degree["/".join(str(x) for x in v)] = int(hit[0])
        out.append({
            "system": f"({n},{d})", "max_degree": mmax,
            "generated_points": len(points),
            "points_violating_a_published_row": len(violations),
            "violations": violations[:5],
            "census_vertices": len(verts),
            "census_vertices_occurring": len(verts) - len(missing),
            "missing_vertices": missing,
            "first_degree_of_each_vertex": first_degree,
            "hull_equals_polytope": not violations and not missing,
            "argument": ("every generated point is an attainable spectrum and "
                         "satisfies every certified row, so conv(points) is "
                         "inside Delta; every certified vertex occurs, so "
                         "Delta is inside conv(points)"),
            "seconds": round(time.time() - started, 1),
            "evidence": "exact",
        })
    return {"what": "conv of normalized occurring highest weights against the "
                    "certified moment polytope",
            "systems": out,
            "all_recovered": all(r["hull_equals_polytope"] for r in out)}


# ---------------------------------------------------------------------------
# C. saturation
# ---------------------------------------------------------------------------

def block_saturation(convention: dict) -> dict:
    by_vertex: dict[tuple[str, int], dict] = {}
    for row in convention["rows"]:
        by_vertex.setdefault((row["system"], row["index"]), {})[row["k"]] = row

    stored = json.loads((DATA / "quantization_multiplicities.json").read_text())
    witnesses, saturated = [], []
    for rec in stored["rows"]:
        q = int(rec["denominator"] or 1)
        key = (rec["system"], int(rec["index"]))
        at_q = by_vertex.get(key, {}).get(q)
        if at_q is None:
            continue
        record = {
            "vertex": f"{rec['system']}#{rec['index']}",
            "integer_form": [int(x) for x in rec["integer_form"]],
            "denominator": q,
            "lattice_point": {"lambda": [int(x) for x in rec["integer_form"]], "k": q},
            "multiplicity_at_k_equals_q": at_q["plethysm"],
            "measured_period": rec.get("quantization_period"),
            "evidence": "exact",
        }
        (witnesses if at_q["plethysm"] == 0 else saturated).append(record)

    characters = json.loads((DATA / "quantization_characters.json").read_text())
    predicted = [r for r in characters.get("census_predictions", [])
                 if r.get("predicted_period_exceeds_denominator")]
    return {
        "what": ("dominant lattice points of cone(Sigma_{N,d}) that are not in "
                 "Sigma_{N,d}; each one is a failure of saturation"),
        "why_the_point_is_in_the_cone": (
            "lambda_int / q is a certified vertex of Delta(N,d), so "
            "(lambda_int, q) = q * (lambda_int/q, 1) lies in cone(Sigma_{N,d}), "
            "and it is an integral dominant weight"),
        "measured_vertices": len(witnesses) + len(saturated),
        "non_saturation_witnesses": witnesses,
        "witness_count": len(witnesses),
        "saturated_at_k_equals_q": len(saturated),
        "predicted_census_wide_witnesses": len(predicted),
        "predicted_witness_vertices": [f"{r['system']}#{r['index']}" for r in predicted],
        "predicted_status": ("PREDICTED, not measured: the character side derives "
                             "divisibility at all 799 vertices and equality is "
                             "verified only at the 14 with a measured sequence"),
        "consequence": ("Sigma_{N,d} is not saturated, so the moment polytope "
                        "determines cone(Sigma_{N,d}) but not Sigma_{N,d}; a "
                        "finite generation statement about the semigroup is "
                        "strictly stronger than one about the cone"),
    }


# ---------------------------------------------------------------------------
# D. padding up the rank ladder
# ---------------------------------------------------------------------------

def consecutive_pairs(ladders: list[list[tuple[int, int]]]) -> list[tuple[int, int, int]]:
    """`(n, d, d+1)` for every consecutive pair inside a single-`N` ladder."""
    out = []
    for ladder in ladders:
        for (n, d), (n2, d2) in zip(ladder, ladder[1:]):
            if n == n2 and d2 == d + 1:
                out.append((n, d, d2))
    return out


def block_padding(ladders: list[list[tuple[int, int]]]) -> dict:
    pairs = []
    for n, d, d2 in consecutive_pairs(ladders):
        ineqs, eqs = published_rows(n, d)
        higher_ineqs, _ = published_rows(n, d2)
        higher_templates = {padding_template(c, b) for c, b in higher_ineqs}
        verts = census_vertices(n, d2)

        def score(rows, kind):
            out = []
            for coeffs, rhs in rows:
                padded = list(coeffs) + [0] * (d2 - d)
                worst, witness = None, None
                for v in verts:
                    slack = Fraction(rhs) - sum(Fraction(c) * x
                                                for c, x in zip(padded, v))
                    if worst is None or slack < worst:
                        worst, witness = slack, v
                out.append({
                    "kind": kind, "coeffs": coeffs, "rhs": rhs,
                    "padded_row_is_valid": bool(worst is not None and worst >= 0),
                    "min_slack": str(worst),
                    "worst_vertex": [str(x) for x in witness] if witness else None,
                    "padded_row_is_published_at_higher_rank":
                        padding_template(padded, rhs) in higher_templates,
                })
            return out

        rows = score(ineqs, "inequality") + score(eqs, "equality")
        pairs.append({
            "from": f"({n},{d})", "to": f"({n},{d2})",
            "rows_tested": len(rows),
            "rows_still_valid": sum(1 for r in rows if r["padded_row_is_valid"]),
            "rows_broken": [r for r in rows if not r["padded_row_is_valid"]],
            "rows_literally_published_higher":
                sum(1 for r in rows if r["padded_row_is_published_at_higher_rank"]),
            "padding_preserves_validity":
                all(r["padded_row_is_valid"] for r in rows),
            "evidence": "exact",
            "proof_strength": ("validity on the certified vertex set is validity "
                               "on the polytope, because P = O is proved at both "
                               "ranks in results/data/census_inner_sandwich.json"),
            "rows": rows,
        })
    return {"what": "does a published row stay valid after zero padding up one rank",
            "pairs": pairs}


# ---------------------------------------------------------------------------
# D2. every order-preserving embedding, not just the trailing one (POST-HOC)
# ---------------------------------------------------------------------------

def block_insertion(ladders: list[list[tuple[int, int]]]) -> dict:
    """Is SOME order-preserving embedding valid, when the trailing one is not?

    POST-HOC. This block was added after block D refuted zero padding, and is
    not pre-registered. Zero padding is the single embedding
    `[d] -> [d+1]`, `i -> i`, so a row that names the last orbital by position
    is being asked to name a different orbital after padding. The honest
    ordered-category question is whether ANY of the `d+1` order-preserving
    injections carries the row, which is exactly the OI action on rows.
    """
    pairs = []
    for n, d, d2 in consecutive_pairs(ladders):
        ineqs, _ = published_rows(n, d)
        verts = census_vertices(n, d2)
        rows = []
        for coeffs, rhs in ineqs:
            valid_positions = []
            for gap in range(d2):
                embedded = list(coeffs[:gap]) + [0] + list(coeffs[gap:])
                if all(sum(Fraction(c) * x for c, x in zip(embedded, v)) <= rhs
                       for v in verts):
                    valid_positions.append(gap)
            rows.append({
                "coeffs": coeffs, "rhs": rhs,
                "valid_insertion_positions": valid_positions,
                "some_embedding_is_valid": bool(valid_positions),
                "trailing_embedding_is_valid": (d2 - 1) in valid_positions,
            })
        lower_images = set()
        for coeffs, rhs in ineqs:
            for gap in range(d2):
                embedded = list(coeffs[:gap]) + [0] + list(coeffs[gap:])
                lower_images.add(_normalize(embedded, rhs))
        higher_ineqs, _ = published_rows(n, d2)
        oi_derived = sum(1 for c, b in higher_ineqs
                         if _normalize(c, b) in lower_images)
        pairs.append({
            "from": f"({n},{d})", "to": f"({n},{d2})",
            "rows_tested": len(rows),
            "rows_with_some_valid_embedding":
                sum(1 for r in rows if r["some_embedding_is_valid"]),
            "rows_with_no_valid_embedding":
                [r for r in rows if not r["some_embedding_is_valid"]],
            "every_row_has_some_valid_embedding":
                all(r["some_embedding_is_valid"] for r in rows),
            "higher_rows": len(higher_ineqs),
            "higher_rows_that_are_oi_images_of_a_lower_row": oi_derived,
            "higher_rows_genuinely_new_under_oi": len(higher_ineqs) - oi_derived,
            "evidence": "exact",
            "rows": rows,
        })
    return {
        "what": "does SOME order-preserving embedding carry a published row up "
                "one rank, when the trailing one does not, and how much of the "
                "higher system is an OI image of the lower one",
        "status": "POST-HOC, added after block D refuted zero padding; not "
                  "pre-registered and not scored",
        "pairs": pairs,
    }


# ---------------------------------------------------------------------------
# E. restriction down the rank ladder
# ---------------------------------------------------------------------------

def block_restriction(ladders: list[list[tuple[int, int]]]) -> dict:
    pairs = []
    for n, d, d2 in consecutive_pairs(ladders):
        higher_ineqs, higher_eqs = published_rows(n, d2)
        lower_ineqs, lower_eqs = published_rows(n, d)

        restricted_ge = [([-Fraction(c) for c in coeffs[:d]], -Fraction(rhs))
                         for coeffs, rhs in higher_ineqs]
        restricted_ge += chamber_rows_ge(d)
        restricted_eq = [([Fraction(c) for c in coeffs[:d]], Fraction(rhs))
                         for coeffs, rhs in higher_eqs]
        restricted_eq.append(([Fraction(1)] * d, Fraction(n)))

        unimplied = []
        for coeffs, rhs in lower_ineqs:
            target = ([-Fraction(c) for c in coeffs], -Fraction(rhs))
            if not implies(restricted_ge, restricted_eq, target):
                unimplied.append({"coeffs": coeffs, "rhs": rhs, "kind": "inequality"})
        for coeffs, rhs in lower_eqs:
            for sign in (1, -1):
                target = ([sign * Fraction(c) for c in coeffs], sign * Fraction(rhs))
                if not implies(restricted_ge, restricted_eq, target):
                    unimplied.append({"coeffs": coeffs, "rhs": rhs,
                                      "kind": "equality", "sign": sign})
        pairs.append({
            "from": f"({n},{d2})", "to": f"({n},{d})",
            "lower_rows_tested": len(lower_ineqs) + 2 * len(lower_eqs),
            "lower_rows_not_implied": unimplied,
            "restriction_is_complete": not unimplied,
            "evidence": "exact",
            "argument": ("Delta(N,d) = Delta(N,d+1) intersect {lambda_{d+1}=0} is a "
                         "theorem, so the restricted system contains Delta(N,d); "
                         "if it also implies every published rank-d row then the "
                         "two systems cut out the same set"),
        })
    return {"what": "does the rank d+1 system imply the rank d system by restriction",
            "pairs": pairs}


# ---------------------------------------------------------------------------
# F. templates
# ---------------------------------------------------------------------------

def block_templates(systems: list[tuple[int, int]]) -> dict:
    per_system, previous_pad, previous_sym = [], set(), set()
    last_n = None
    for n, d in systems:
        if n != last_n:
            previous_pad, previous_sym = set(), set()
            last_n = n
        ineqs, _ = published_rows(n, d)
        pads = {padding_template(c, b) for c, b in ineqs}
        syms = {symmetric_template(c, b) for c, b in ineqs}
        new_pad = pads - previous_pad
        new_sym = syms - previous_sym
        arities = [len([c for c in coeffs if c]) for coeffs, _ in ineqs]
        per_system.append({
            "system": f"({n},{d})",
            "published_inequalities": len(ineqs),
            "distinct_padding_templates": len(pads),
            "distinct_symmetric_templates": len(syms),
            "new_padding_templates_vs_previous_rank": len(new_pad),
            "new_symmetric_templates_vs_previous_rank": len(new_sym),
            "max_symmetric_arity": max(arities) if arities else 0,
            "symmetric_arity_histogram": {
                str(a): arities.count(a) for a in sorted(set(arities))},
        })
        previous_pad |= pads
        previous_sym |= syms
    return {
        "what": "template stock of the published rows, per rank",
        "padding_template": "coefficient vector with trailing zeros stripped, gcd "
                            "normalized; the ordered (OI) notion",
        "symmetric_template": "multiset of nonzero coefficients, gcd normalized; "
                              "the unordered (FI, ambient cone) notion",
        "why_arity_matters": (
            "if every facet of the ambient Weyl-invariant moment cone at every "
            "rank were an instance of a fixed finite template list of maximum "
            "symmetric arity a, the cone at rank d would have O(d^a) facets, "
            "because an FI orbit of an arity-a template has at most "
            "binom(d,a) * a! members. Unbounded arity is therefore the first "
            "observable symptom that no finite template list exists"),
        "per_system": per_system,
    }


# ---------------------------------------------------------------------------
# G. what a template list would cost at rank 20 and rank 40 (PROJECTION)
# ---------------------------------------------------------------------------

def block_cost_projection(n: int = 3, d: int = 10,
                          targets=(11, 12, 20, 40)) -> dict:
    """Instance counts if the rank-`d` templates were the whole template list.

    PROJECTION, not a prediction. It answers one question exactly: SUPPOSING the
    missing templating theorem held with exactly the symmetric templates
    observed at rank `d`, how many candidate rows would a constructor have to
    consider at rank 20 and rank 40? An arity-`a` template with repeated
    coefficients has `binom(D,a) * a! / prod(mult!)` ordered instances at rank
    `D`, because the ordered chamber distinguishes arrangements.
    """
    from math import comb, factorial

    ineqs, _ = published_rows(n, d)
    templates = sorted({symmetric_template(c, b) for c, b in ineqs})
    rows = []
    for coeffs, rhs in templates:
        arity = len(coeffs)
        repeats = 1
        for value in set(coeffs):
            repeats *= factorial(coeffs.count(value))
        arrangements = factorial(arity) // repeats
        rows.append({
            "template": list(coeffs), "rhs": rhs, "arity": arity,
            "distinct_arrangements": arrangements,
            "instances": {str(D): comb(D, arity) * arrangements for D in targets},
        })
    totals = {str(D): sum(r["instances"][str(D)] for r in rows) for D in targets}
    return {
        "what": "candidate-row count at higher rank IF the rank-%d symmetric "
                "templates were the complete template list" % d,
        "status": "PROJECTION, not a prediction, and an underestimate of the "
                  "template list by construction: new templates appear at every "
                  "rank measured",
        "source_system": f"({n},{d})",
        "templates": len(templates),
        "per_template": rows,
        "total_candidate_rows": totals,
        "comparison": "the orbit-canonical enumerator reaches 10,004,154 "
                      "admissible hyperplanes at (3,9) and did not finish (3,11) "
                      "before streaming was added, so a candidate count above "
                      "about 1e8 is not obviously an improvement",
    }


# ---------------------------------------------------------------------------
# scorecard
# ---------------------------------------------------------------------------

def scorecard(convention, recovery, saturation, padding, restriction, templates) -> dict:
    n3 = [r for r in templates["per_system"] if r["system"].startswith("(3,")]
    pad_pairs = padding["pairs"]
    above6 = [p for p in pad_pairs if p["from"] != "(3,6)"]
    pilot = [p for p in pad_pairs if p["from"] == "(3,6)"]

    def strictly_increasing(values):
        return all(b > a for a, b in zip(values, values[1:]))

    return {
        "P1_three_routes_agree": {
            "verdict": "PASS" if convention["all_agree"] else "FAIL",
            "kind": "forced correctness gate, uninformative as evidence",
            "scope": convention["pairs"],
            "independent_scope": 13,
            "independent_distinct_states": 13,
            "weyl_recomputed_pairs": convention["weyl_recomputed_pairs"],
            "mismatches": len(convention["mismatches"]),
        },
        "P2_convention_is_load_bearing": {
            "verdict": "PASS" if convention["control_is_load_bearing"] else "FAIL",
            "scope": convention["pairs"],
            "independent_scope": 13,
            "independent_distinct_states": 13,
            "pairs_where_a_wrong_convention_differs":
                convention["control_pairs_where_a_wrong_convention_differs"],
        },
        "P3_semigroup_is_not_saturated": {
            "verdict": "PASS" if saturation["witness_count"] == 4 else "FAIL",
            "kind": "ALREADY-KNOWN under R4; the numbers are in "
                    "results/data/quantization_multiplicities.json, only the "
                    "semigroup reading is new",
            "scope": saturation["measured_vertices"],
            "independent_scope": 13,
            "independent_distinct_states": 13,
            "witnesses": [w["vertex"] for w in saturation["non_saturation_witnesses"]],
        },
        "P4_zero_padding_preserves_validity_above_rank_6": {
            "verdict": "PASS" if all(p["padding_preserves_validity"] for p in above6)
                       else "FAIL",
            "scope": sum(p["rows_tested"] for p in above6),
            "independent_scope": len(above6),
            "independent_distinct_states": len(above6),
            "pairs": {p["from"] + " -> " + p["to"]: p["padding_preserves_validity"]
                      for p in pad_pairs},
            "disclosed_pilot_3_6_to_3_7": {
                "expected": "FAIL to propagate",
                "observed": [p["padding_preserves_validity"] for p in pilot],
            },
        },
        "P5_restriction_is_complete": {
            "verdict": "PASS" if all(p["restriction_is_complete"]
                                     for p in restriction["pairs"]) else "FAIL",
            "kind": "forced correctness gate on the published tables",
            "scope": sum(p["lower_rows_tested"] for p in restriction["pairs"]),
            "independent_scope": len(restriction["pairs"]),
            "independent_distinct_states": len(restriction["pairs"]),
            "pairs": {p["from"] + " -> " + p["to"]: p["restriction_is_complete"]
                      for p in restriction["pairs"]},
        },
        "P6_no_template_saturation_in_the_window": {
            "verdict": "PASS" if all(
                r["new_padding_templates_vs_previous_rank"] >= 1 for r in n3[1:]) and all(
                r["new_symmetric_templates_vs_previous_rank"] >= 1 for r in n3[1:])
                else "FAIL",
            "scope": sum(r["published_inequalities"] for r in n3),
            "independent_scope": len(n3) - 1,
            "independent_distinct_states": len(n3) - 1,
            "new_padding_templates": {r["system"]:
                                      r["new_padding_templates_vs_previous_rank"]
                                      for r in n3},
            "new_symmetric_templates": {r["system"]:
                                        r["new_symmetric_templates_vs_previous_rank"]
                                        for r in n3},
        },
        "P7_symmetric_arity_strictly_increases": {
            "verdict": "PASS" if strictly_increasing(
                [r["max_symmetric_arity"] for r in n3]) else "FAIL",
            "scope": sum(r["published_inequalities"] for r in n3),
            "independent_scope": len(n3),
            "independent_distinct_states": len(n3),
            "max_symmetric_arity": {r["system"]: r["max_symmetric_arity"] for r in n3},
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weyl-max-k", type=int, default=8,
                    help="recompute the Weyl alternating sum only up to this k")
    ap.add_argument("--weyl-budget", type=float, default=600.0,
                    help="seconds of Weyl recomputation before falling back to "
                         "the two-route check")
    ap.add_argument("--quick", action="store_true",
                    help="skip the recovery block, which dominates the runtime")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse the convention and recovery blocks of an "
                         "existing artifact instead of recomputing them")
    ap.add_argument("-o", "--out", default=str(OUT))
    args = ap.parse_args()

    started = time.time()
    cached = (json.loads(pathlib.Path(args.out).read_text())
              if args.reuse and pathlib.Path(args.out).exists() else None)

    print("A. convention ...", flush=True)
    convention = (cached["convention"] if cached
                  else block_convention(args.weyl_max_k, args.weyl_budget))
    print(f"   {convention['pairs']} pairs, all agree: {convention['all_agree']}",
          flush=True)

    print("B. recovery ...", flush=True)
    if cached:
        recovery = cached["recovery"]
    elif args.quick:
        recovery = {"what": "skipped (--quick)", "systems": [], "all_recovered": None}
    else:
        recovery = block_recovery()
    for rec in recovery["systems"]:
        print(f"   {rec['system']}: {rec['generated_points']} points, "
              f"hull equals polytope: {rec['hull_equals_polytope']}", flush=True)

    print("C. saturation ...", flush=True)
    saturation = block_saturation(convention)
    print(f"   {saturation['witness_count']} measured non-saturation witnesses",
          flush=True)

    print("D. padding ...", flush=True)
    padding = block_padding([LADDER_N3, LADDER_N4])
    for pair in padding["pairs"]:
        print(f"   {pair['from']} -> {pair['to']}: valid "
              f"{pair['rows_still_valid']}/{pair['rows_tested']}", flush=True)

    print("D2. order-preserving embeddings (post-hoc) ...", flush=True)
    insertion = block_insertion([LADDER_N3, LADDER_N4])
    for pair in insertion["pairs"]:
        print(f"   {pair['from']} -> {pair['to']}: some embedding valid for "
              f"{pair['rows_with_some_valid_embedding']}/{pair['rows_tested']}, "
              f"higher rows new under OI "
              f"{pair['higher_rows_genuinely_new_under_oi']}/{pair['higher_rows']}",
              flush=True)

    print("E. restriction ...", flush=True)
    restriction = block_restriction([LADDER_N3, LADDER_N4])
    for pair in restriction["pairs"]:
        print(f"   {pair['from']} -> {pair['to']}: complete "
              f"{pair['restriction_is_complete']}", flush=True)

    print("F. templates ...", flush=True)
    templates = block_templates(ALL_SYSTEMS)
    for rec in templates["per_system"]:
        print(f"   {rec['system']}: {rec['published_inequalities']} rows, "
              f"{rec['distinct_symmetric_templates']} symmetric templates, "
              f"max arity {rec['max_symmetric_arity']}", flush=True)

    print("G. cost projection ...", flush=True)
    projection = block_cost_projection()
    print(f"   candidate rows if the (3,10) templates were complete: "
          f"{projection['total_candidate_rows']}", flush=True)

    payload = {
        "generated_by": "scripts/infinite_wedge_tca_sanity.py",
        "prereg": "docs/prereg_infinite_wedge_tca_sanity.md",
        "note": "docs/infinite_wedge_tca_bridge.md",
        "dependency_artifact": "results/data/infinite_wedge_tca_dependency.json",
        "evidence": "exact",
        "scope": ("census systems only, ranks 6 to 10. Nothing here is a "
                  "statement about all d."),
        "convention": convention,
        "recovery": recovery,
        "saturation": saturation,
        "padding": padding,
        "insertion_post_hoc": insertion,
        "restriction": restriction,
        "templates": templates,
        "cost_projection": projection,
        "scorecard": scorecard(convention, recovery, saturation,
                               padding, restriction, templates),
        "wall_secs": round(time.time() - started, 1),
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")
    for name, block in payload["scorecard"].items():
        print(f"  {name}: {block['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
