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


def load_states():
    for line in (ROOT / "results/data/states.jsonl").read_text().splitlines():
        if line.strip():
            yield json.loads(line)


def _degree(value) -> int | None:
    """Exact degree over Q of an algebraic number, or None if it is not one."""
    try:
        return sp.Poly(sp.minimal_polynomial(sp.nsimplify(value), X), X).degree()
    except Exception:
        return None


def _exact_mod2(expr):
    """Exact |expr|^2 for a symbolic complex amplitude expression."""
    return sp.radsimp(sp.expand(sp.simplify(sp.expand(expr * sp.conjugate(expr)))))


def orient_class(dets, a_mode, b_mode, pairs):
    """Return [(u, v, sign)] with u always the A-side, matching Lemma 2."""
    oriented = []
    for i, j in pairs:
        u, v = (i, j) if a_mode in dets[i] else (j, i)
        sign = ((-1) ** dets[v].index(b_mode)) * ((-1) ** dets[u].index(a_mode))
        oriented.append((u, v, sign))
    return oriented


def two_rdm_element(amps, dets, a_mode, b_mode, c_mode):
    """The 2-RDM entry rho^(2)_{AC,BC} = <psi| a_A^+ a_C^+ a_C a_B |psi>.

    Computed from scratch by acting with the operator on each determinant, so
    the fermionic signs are derived rather than assumed to match `s_a`.
    """
    index = {tuple(t): k for k, t in enumerate(dets)}
    total = 0
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
        total += sign * sp.conjugate(amps[target]) * amps[k]
    return sp.simplify(total)


def analyse_class(rec, dets, amps, a_mode, b_mode, pairs) -> dict:
    oriented = orient_class(dets, a_mode, b_mode, pairs)
    m = len(oriented)
    z = [s * sp.conjugate(amps[u]) * amps[v] for u, v, s in oriented]
    modes = sorted({mode for det in dets for mode in det} - {a_mode, b_mode})

    # mode-coherent subsets: T_C = {a : C in u_a}
    coherent: dict[frozenset, list[int]] = {}
    for c_mode in modes:
        T = frozenset(a for a, (u, _, _) in enumerate(oriented) if c_mode in dets[u])
        if 0 < len(T) < m:
            coherent.setdefault(T, []).append(c_mode)

    rho = sp.simplify(sum(z))
    subsets = []
    for size in range(2, m):
        for T in itertools.combinations(range(m), size):
            key = frozenset(T)
            value = _exact_mod2(sum(z[a] for a in T))
            degree = _degree(value)
            entry = {
                "subset": list(T),
                "size": size,
                "mode_coherent": key in coherent,
                "coherent_modes": coherent.get(key, []),
                "value": str(value),
                "value_numeric": float(sp.re(sp.N(value, 40))),
                "degree_over_Q": degree,
                "rational": degree == 1,
            }
            if size == 2 and key in coherent:
                c_mode = coherent[key][0]
                element = two_rdm_element(amps, dets, a_mode, b_mode, c_mode)
                target = _exact_mod2(element)
                entry["two_rdm_mod2"] = str(target)
                entry["matches_two_rdm"] = bool(
                    sp.simplify(sp.nsimplify(value - target)) == 0
                )
            subsets.append(entry)

    # the within-class holonomy cosines, which is where (R1) lives
    cosines = []
    for a, b in itertools.combinations(range(m), 2):
        ua, va, sa = oriented[a]
        ub, vb, sb = oriented[b]
        t_ab = sp.simplify(sp.re(z[a] * sp.conjugate(z[b])))
        product = sp.simplify(_exact_mod2(z[a]) * _exact_mod2(z[b]))
        cos2 = sp.simplify(t_ab**2 / product) if product != 0 else None
        cosines.append({
            "pair": [a, b],
            "t_ab": str(t_ab),
            "t_ab_degree": _degree(t_ab),
            "t_ab_rational": _degree(t_ab) == 1,
            "cos_squared": str(cos2),
            "cos_squared_degree": _degree(cos2),
            "cos_squared_rational": _degree(cos2) == 1,
        })

    return {
        "system": rec["system"],
        "index": rec["index"],
        "classified": rec["classified"],
        "channel": [a_mode, b_mode],
        "class_size": m,
        "support_size": len(dets),
        "oriented": [[u, v, s] for u, v, s in oriented],
        "r_squared": [str(_exact_mod2(z[a])) for a in range(m)],
        "rho_mod2": str(_exact_mod2(rho)),
        "rho_mod2_rational": _degree(_exact_mod2(rho)) == 1,
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
        amps = [sp.sympify(p) for p in rec["closed_form"]["pretty"]]
        for (a_mode, b_mode), pairs in sorted(big.items()):
            print(f"{rec['system']} v{rec['index']} channel ({a_mode},{b_mode}) "
                  f"size {len(pairs)} ...", flush=True)
            blocks.append(analyse_class(rec, dets, amps, a_mode, b_mode, pairs))

    every = [s for b in blocks for s in b["proper_subsets"]]
    coherent = [s for s in every if s["mode_coherent"]]
    incoherent = [s for s in every if not s["mode_coherent"]]
    cosines = [c for b in blocks for c in b["within_class_cosines"]]
    checked = [s for s in coherent if "matches_two_rdm" in s]

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
            "class_sizes": {str(k): sum(1 for b in blocks if b["class_size"] == k)
                            for k in sorted({b["class_size"] for b in blocks})},
            "proper_subsets": len(every),
            "rational_subsets": sum(1 for s in every if s["rational"]),
            "irrational_subsets": sum(1 for s in every if not s["rational"]),
            "mode_coherent_subsets": len(coherent),
            "mode_coherent_rational": sum(1 for s in coherent if s["rational"]),
            "mode_incoherent_subsets": len(incoherent),
            "mode_incoherent_rational": sum(1 for s in incoherent if s["rational"]),
            "two_rdm_checked": len(checked),
            "two_rdm_matches": sum(1 for s in checked if s.get("matches_two_rdm")),
            "irrational_degrees": {
                str(d): sum(1 for s in every if not s["rational"]
                            and s["degree_over_Q"] == d)
                for d in sorted({s["degree_over_Q"] for s in every
                                 if not s["rational"]} - {None})
            },
            "within_class_pairs": len(cosines),
            "cos_squared_rational": sum(1 for c in cosines
                                        if c["cos_squared_rational"]),
            "t_ab_rational": sum(1 for c in cosines if c["t_ab_rational"]),
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
