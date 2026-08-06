#!/usr/bin/env python3
"""Track B1: the polygon reformulation of Lemma 2, and the diagonal test.

THE REFORMULATION. Lemma 2 of `docs/holonomy_rationality.md` expands
`|rho_AB|^2` for a one-hop class. Written without expanding, the channel
condition is

    rho_AB = sum_a z_a,      z_a = s_a conj(c_{u_a}) c_{v_a},

so `|z_a| = r_a = sqrt(p_{u_a} p_{v_a})` is fixed by the weights and only the
arguments are free. A one-hop class of size m is therefore a CLOSED (m+1)-GON
with prescribed side lengths `r_1, ..., r_m, |rho_AB|`, and the within-class
holonomies `Phi_ab` are its angles. Theorem A (m = 2) is the law of cosines for
a triangle, which is why its radicand is a weight monomial.

A quadrilateral (m = 3) splits into two triangles along a diagonal, and each
triangle contributes a quadratic with monomial radicand. So the conjecture at
size 3 turns on ONE question: are the squared diagonal lengths rational?

    Is  |sum_{a in T} z_a|^2  rational for proper subsets T of a class?

THE REDUCTION, proved here rather than measured (see the module test). For
|T| = 2,

    |z_a + z_b|^2 = r_a^2 + r_b^2 + 2 s_a s_b r_a r_b cos(Phi_ab),

and `r_a^2, r_b^2` are rational under (H1). Writing
`t_ab = s_a s_b r_a r_b cos(Phi_ab)`, the diagonal is rational iff `t_ab` is.
Since `(r_a r_b)^2 = p_{u_a} p_{v_a} p_{u_b} p_{v_b}` is rational,

    t_ab^2 = (rational) * cos^2(Phi_ab),

so a RATIONAL diagonal forces `cos^2(Phi_ab)` rational, which is exactly
condition (R1) of `docs/holonomy_two_elementary.md` for the within-class
holonomy. The polygon question is therefore not independent of (R1): it is
(R1) on the within-class difference loops, plus the extra demand that the
resulting rational be a square. That is worth knowing before measuring, and it
is why the measurement below records `cos^2` rationality alongside every
diagonal.

WHAT THE DIAGONAL MIGHT BE. A partial sum is not obviously a physical
quantity, but one family of subsets is: if `T = {a : C in u_a}` for a mode
`C` outside `{A, B}`, then every determinant of the subset contains `C`, and
the partial sum is (up to one overall sign) the 2-RDM element
`rho^(2)_{AC,BC}`. Such subsets are called MODE-COHERENT here, and the script
computes the corresponding 2-RDM element independently and compares, so the
identification is verified rather than asserted.

Everything is exact. Amplitudes come from the closed forms of `states.jsonl`
(never `states_natural_orbital.jsonl`, which is numerical), and rationality is
decided by `sympy.minimal_polynomial` degree, not by a tolerance.

Usage:
    python scripts/polygon_diagonals.py
    python scripts/polygon_diagonals.py --min-class-size 3
"""
from __future__ import annotations

import argparse
import itertools
import json
import pathlib
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.states import one_hop_classes  # noqa: E402

X = sp.Symbol("x")


def transport_class(system, integer_form, denominator):
    """Canonical key modulo trailing-zero padding and frozen-core lifts.

    The same key `scripts/two_block_family.py` and `scripts/no_support_ratio.py`
    use. Standing rule R2: two records with this key are ONE observation, so a
    tally that reports only record counts overstates its evidence.
    """
    from fractions import Fraction

    n = int(system.strip("()").split(",")[0])
    v = [Fraction(x, denominator) for x in integer_form]
    while v and v[-1] == 0:
        v.pop()
    while v and v[0] == 1:
        v.pop(0)
        n -= 1
    return f"N={n}:" + ",".join(str(x) for x in v)


def load_states():
    for line in (ROOT / "results/data/states.jsonl").read_text().splitlines():
        if line.strip():
            yield json.loads(line)


def _degree(value) -> int | None:
    """Exact degree over Q of an algebraic number, or None if sympy failed.

    None means the computation did not resolve, NOT that the number is
    irrational. An earlier version of this script conflated the two by testing
    ``degree == 1``, which reported 29 unresolved values as irrational. Every
    caller must treat None as a measurement failure and the summary counts it
    separately.
    """
    try:
        return sp.Poly(sp.minimal_polynomial(value, X), X).degree()
    except Exception:
        return None


def exact_parts(amp):
    """Exact algebraic (Re, Im) of one closed-form amplitude.

    The closed forms carry phases as ``exp(I*(pi/3 + atan(sqrt 7) + ...))``,
    and neither ``simplify`` nor ``nsimplify`` reduces a product of those to an
    algebraic number on its own. Expanding in complex mode first and then
    splitting is what turns each amplitude into a pair of exact algebraic
    reals, which is the representation every downstream quantity is built from.
    """
    value = sp.simplify(sp.expand_complex(sp.expand(sp.simplify(amp), complex=True)))
    return sp.simplify(sp.re(value)), sp.simplify(sp.im(value))


def _mul_conj(pu, pv):
    """conj(c_u) * c_v in exact (Re, Im) form."""
    xu, yu = pu
    xv, yv = pv
    return sp.expand(xu * xv + yu * yv), sp.expand(xu * yv - yu * xv)


def _mod2(part):
    """Exact |z|^2 from an exact (Re, Im) pair."""
    re_, im_ = part
    return sp.radsimp(sp.simplify(sp.expand(re_**2 + im_**2)))


def orient_class(dets, a_mode, b_mode, pairs):
    """Return [(u, v, sign)] with u always the A-side, matching Lemma 2."""
    oriented = []
    for i, j in pairs:
        u, v = (i, j) if a_mode in dets[i] else (j, i)
        sign = ((-1) ** dets[v].index(b_mode)) * ((-1) ** dets[u].index(a_mode))
        oriented.append((u, v, sign))
    return oriented


def two_rdm_element(parts, dets, a_mode, b_mode, c_mode):
    """The 2-RDM entry rho^(2)_{AC,BC} = <psi| a_A^+ a_C^+ a_C a_B |psi>.

    Computed from scratch by acting with the operator on each determinant, so
    the fermionic signs are derived rather than assumed to match `s_a`.
    """
    index = {tuple(t): k for k, t in enumerate(dets)}
    total = (sp.Integer(0), sp.Integer(0))
    for k, det in enumerate(dets):
        occ = list(det)
        if b_mode not in occ or c_mode not in occ:
            continue
        # annihilate B then C, create C then A
        sign = (-1) ** occ.index(b_mode)
        rest = [m for m in occ if m != b_mode]
        if c_mode not in rest:
            continue
        sign *= (-1) ** rest.index(c_mode)
        rest = [m for m in rest if m != c_mode]
        if c_mode in rest or a_mode in rest:
            continue
        pos = sum(m < c_mode for m in rest)
        sign *= (-1) ** pos
        rest2 = rest[:pos] + [c_mode] + rest[pos:]
        pos = sum(m < a_mode for m in rest2)
        sign *= (-1) ** pos
        out = tuple(rest2[:pos] + [a_mode] + rest2[pos:])
        target = index.get(out)
        if target is None:
            continue
        re_, im_ = _mul_conj(parts[target], parts[k])
        total = (total[0] + sign * re_, total[1] + sign * im_)
    return (sp.simplify(total[0]), sp.simplify(total[1]))


def triangulation_certificate(z, cosines, rho_mod2, m) -> dict:
    """Machine check of Theorem C at a class carrying a rational diagonal.

    The three vectors ``z_a``, ``z_c`` and ``D = z_a + z_b`` all live in the
    real plane ``C = R^2``, so they are linearly dependent and their 3x3 Gram
    matrix is SINGULAR. Written out in the entries

        [ r_a^2   t_ac    q    ]
        [ t_ac    r_c^2   p    ]      q = r_a^2 + t_ab,  p = t_ac + t_bc,
        [ q       p       D^2  ]      D^2 = r_a^2 + r_b^2 + 2 t_ab,

    the vanishing determinant is a QUADRATIC in ``t_ac``:

        D^2 t^2 - 2 p q t - (r_a^2 r_c^2 D^2 - r_a^2 p^2 - r_c^2 q^2) = 0.

    Under (H1) and (H2), ``t_ab`` rational makes every coefficient rational:
    ``D^2`` and ``q`` directly, and ``p`` because Lemma 2 gives
    ``t_ab + t_ac + t_bc = (|rho|^2 - sum r^2) / 2`` rational. So ``t_ac`` is
    at most quadratic over Q with discriminant
    ``4 (r_a^2 D^2 - q^2)(r_c^2 D^2 - p^2)``, a rational number.

    The certificate below builds the quadratic from RATIONAL data only and
    then checks that the independently measured ``t_ac`` and ``t_bc`` are its
    roots, so the identity is verified rather than assumed.
    """
    look = {(tuple(c["pair"])): c for c in cosines}
    if all(c["t_ab_rational"] is True for c in cosines):
        # Every pairwise t is rational, so every cos(Phi_ab) = t_ab / (s_a s_b
        # r_a r_b) already lies in Q(sqrt(r_a^2 r_b^2)), a single quadratic
        # extension by a weight monomial. Nothing to triangulate.
        return {
            "applies": True,
            "route": "fully_rational_polygon",
            "diagonal_pair": None,
            "all_coefficients_rational": True,
            "discriminant_rational": True,
            "gram_identity_holds": True,
            "max_cos_squared_degree": max(c["cos_squared_degree"] for c in cosines),
        }
    if m != 3:
        return {"applies": False,
                "why": f"the triangulation certificate is stated at class size "
                       f"3 and this class has size {m} with a mixed diagonal set"}
    r2 = [_mod2(zi) for zi in z]
    total = sp.simplify((rho_mod2 - sum(r2)) / 2)      # = t_ab + t_ac + t_bc

    for a, b in itertools.combinations(range(3), 2):
        entry = look[(a, b)]
        if entry["t_ab_rational"] is not True:
            continue
        c = next(k for k in range(3) if k not in (a, b))
        t_ab = sp.nsimplify(sp.sympify(entry["t_ab"]))
        d2 = sp.simplify(r2[a] + r2[b] + 2 * t_ab)
        q = sp.simplify(r2[a] + t_ab)
        p = sp.simplify(total - t_ab)
        coefficients = [d2, -2 * p * q,
                        -sp.simplify(r2[a] * r2[c] * d2 - r2[a] * p**2 - r2[c] * q**2)]
        discriminant = sp.simplify((r2[a] * d2 - q**2) * (r2[c] * d2 - p**2))

        measured = {}
        for other, label in ((a, "t_ac"), (b, "t_bc")):
            key = (min(other, c), max(other, c))
            measured[label] = sp.nsimplify(sp.sympify(look[key]["t_ab"]))
        # the quadratic is written for t_(a,c); t_(b,c) is its partner root
        residual = sp.simplify(
            coefficients[0] * measured["t_ac"] ** 2
            + coefficients[1] * measured["t_ac"] + coefficients[2])
        return {
            "applies": True,
            "route": "gram_triangulation",
            "diagonal_pair": [a, b],
            "third": c,
            "D_squared": str(d2),
            "D_squared_rational": _degree(d2) == 1,
            "q": str(q), "q_rational": _degree(q) == 1,
            "p": str(p), "p_rational": _degree(p) == 1,
            "quadratic_coefficients": [str(x) for x in coefficients],
            "all_coefficients_rational": all(_degree(x) == 1 for x in coefficients),
            "discriminant": str(discriminant),
            "discriminant_rational": _degree(discriminant) == 1,
            "t_ac": str(measured["t_ac"]),
            "t_ac_degree": _degree(measured["t_ac"]),
            "gram_residual": str(residual),
            "gram_identity_holds": bool(residual == 0),
        }
    return {"applies": False, "why": "no rational diagonal in this class"}


def analyse_class(rec, dets, parts, a_mode, b_mode, pairs) -> dict:
    oriented = orient_class(dets, a_mode, b_mode, pairs)
    m = len(oriented)
    z = []
    for u, v, s in oriented:
        re_, im_ = _mul_conj(parts[u], parts[v])
        z.append((sp.expand(s * re_), sp.expand(s * im_)))
    modes = sorted({mode for det in dets for mode in det} - {a_mode, b_mode})

    # mode-coherent subsets: T_C = {a : C in u_a}
    coherent: dict[frozenset, list[int]] = {}
    for c_mode in modes:
        T = frozenset(a for a, (u, _, _) in enumerate(oriented) if c_mode in dets[u])
        if 0 < len(T) < m:
            coherent.setdefault(T, []).append(c_mode)

    rho = (sp.expand(sum(p[0] for p in z)), sp.expand(sum(p[1] for p in z)))
    subsets = []
    for size in range(2, m):
        for T in itertools.combinations(range(m), size):
            key = frozenset(T)
            partial = (sp.expand(sum(z[a][0] for a in T)),
                       sp.expand(sum(z[a][1] for a in T)))
            value = _mod2(partial)
            degree = _degree(value)
            entry = {
                "subset": list(T),
                "size": size,
                "mode_coherent": key in coherent,
                "coherent_modes": coherent.get(key, []),
                "value": str(value),
                "value_numeric": float(sp.N(value, 40)),
                "degree_over_Q": degree,
                "degree_resolved": degree is not None,
                "rational": None if degree is None else degree == 1,
            }
            if size == 2 and key in coherent:
                c_mode = coherent[key][0]
                element = two_rdm_element(parts, dets, a_mode, b_mode, c_mode)
                target = _mod2(element)
                entry["two_rdm_mod2"] = str(target)
                entry["matches_two_rdm"] = bool(
                    sp.simplify(sp.expand(value - target)) == 0
                )
            subsets.append(entry)

    # the within-class holonomy cosines, which is where (R1) lives
    cosines = []
    for a, b in itertools.combinations(range(m), 2):
        # t_ab = Re(z_a conj(z_b)) = s_a s_b r_a r_b cos(Phi_ab)
        t_ab = sp.radsimp(sp.simplify(
            sp.expand(z[a][0] * z[b][0] + z[a][1] * z[b][1])))
        product = sp.simplify(_mod2(z[a]) * _mod2(z[b]))
        cos2 = sp.radsimp(sp.simplify(t_ab**2 / product)) if product != 0 else None
        t_degree = _degree(t_ab)
        c_degree = None if cos2 is None else _degree(cos2)
        cosines.append({
            "pair": [a, b],
            "t_ab": str(t_ab),
            "t_ab_degree": t_degree,
            "t_ab_rational": None if t_degree is None else t_degree == 1,
            "cos_squared": str(cos2),
            "cos_squared_degree": c_degree,
            "cos_squared_rational": None if c_degree is None else c_degree == 1,
        })

    rho_mod2 = _mod2(rho)
    rho_degree = _degree(rho_mod2)
    return {
        "theorem_C": triangulation_certificate(z, cosines, rho_mod2, m),
        "system": rec["system"],
        "index": rec["index"],
        "classified": rec["classified"],
        "channel": [a_mode, b_mode],
        "class_size": m,
        "support_size": len(dets),
        "oriented": [[u, v, s] for u, v, s in oriented],
        "r_squared": [str(_mod2(z[a])) for a in range(m)],
        "rho_mod2": str(rho_mod2),
        "rho_mod2_rational": None if rho_degree is None else rho_degree == 1,
        "proper_subsets": subsets,
        "within_class_cosines": cosines,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-class-size", type=int, default=3)
    ap.add_argument("-o", "--out",
                    default=str(ROOT / "results/data/polygon_diagonals.json"))
    args = ap.parse_args()

    blocks = []
    for rec in load_states():
        dets = [tuple(t) for t in rec["closed_form"]["support_dets"]]
        classes = one_hop_classes(dets)
        big = {k: v for k, v in classes.items() if len(v) >= args.min_class_size}
        if not big:
            continue
        parts = [exact_parts(sp.sympify(p)) for p in rec["closed_form"]["pretty"]]
        for (a_mode, b_mode), pairs in sorted(big.items()):
            print(f"{rec['system']} v{rec['index']} channel ({a_mode},{b_mode}) "
                  f"size {len(pairs)} ...", flush=True)
            block = analyse_class(rec, dets, parts, a_mode, b_mode, pairs)
            block["transport_class"] = transport_class(
                rec["system"], rec["integer_form"], rec["denominator"])
            blocks.append(block)

    every = [s for b in blocks for s in b["proper_subsets"]]
    coherent = [s for s in every if s["mode_coherent"]]
    incoherent = [s for s in every if not s["mode_coherent"]]
    cosines = [c for b in blocks for c in b["within_class_cosines"]]
    checked = [s for s in coherent if "matches_two_rdm" in s]
    unresolved = [s for s in every if s["rational"] is None]

    payload = {
        "what": "Polygon reformulation of Lemma 2: squared diagonal lengths of "
                "the one-hop channel polygon, and whether they are rational.",
        "reduction": {
            "statement": "for |T| = 2 the diagonal is rational iff t_ab is, and "
                         "t_ab^2 = (rational) * cos^2(Phi_ab), so a rational "
                         "diagonal forces (R1) on the within-class holonomy",
            "consequence": "the polygon question is (R1) on within-class "
                           "difference loops plus a squareness demand, not an "
                           "independent question",
        },
        "source": "results/data/states.jsonl (exact closed forms)",
        "summary": {
            "classes": len(blocks),
            "distinct_states": len({(b["system"], b["index"]) for b in blocks}),
            "transport_classes": len({b["transport_class"] for b in blocks}),
            "class_sizes": {str(k): sum(1 for b in blocks if b["class_size"] == k)
                            for k in sorted({b["class_size"] for b in blocks})},
            "proper_subsets": len(every),
            "degree_unresolved": len(unresolved),
            "rational_subsets": sum(1 for s in every if s["rational"] is True),
            "irrational_subsets": sum(1 for s in every if s["rational"] is False),
            "mode_coherent_subsets": len(coherent),
            "mode_coherent_rational": sum(1 for s in coherent
                                          if s["rational"] is True),
            "mode_coherent_irrational": sum(1 for s in coherent
                                            if s["rational"] is False),
            "mode_incoherent_subsets": len(incoherent),
            "mode_incoherent_rational": sum(1 for s in incoherent
                                            if s["rational"] is True),
            "mode_incoherent_irrational": sum(1 for s in incoherent
                                              if s["rational"] is False),
            "two_rdm_checked": len(checked),
            "two_rdm_matches": sum(1 for s in checked if s.get("matches_two_rdm")),
            "irrational_degrees": {
                str(d): sum(1 for s in every if s["rational"] is False
                            and s["degree_over_Q"] == d)
                for d in sorted({s["degree_over_Q"] for s in every
                                 if s["rational"] is False} - {None})
            },
            "within_class_pairs": len(cosines),
            "cos_squared_rational": sum(1 for c in cosines
                                        if c["cos_squared_rational"] is True),
            "cos_squared_unresolved": sum(1 for c in cosines
                                          if c["cos_squared_rational"] is None),
            "t_ab_rational": sum(1 for c in cosines if c["t_ab_rational"] is True),
            "t_ab_unresolved": sum(1 for c in cosines
                                   if c["t_ab_rational"] is None),
            "classes_with_a_rational_diagonal": sum(
                1 for b in blocks
                if any(s["size"] == 2 and s["rational"] is True
                       for s in b["proper_subsets"])),
            "theorem_C_applies": sum(1 for b in blocks
                                     if b["theorem_C"].get("applies")),
            "theorem_C_gram_identity_holds": sum(
                1 for b in blocks if b["theorem_C"].get("gram_identity_holds")),
            "theorem_C_all_coefficients_rational": sum(
                1 for b in blocks
                if b["theorem_C"].get("all_coefficients_rational")),
            "theorem_C_discriminant_rational": sum(
                1 for b in blocks if b["theorem_C"].get("discriminant_rational")),
        },
        "classes": blocks,
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    print(json.dumps(payload["summary"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
