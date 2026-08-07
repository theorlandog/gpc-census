#!/usr/bin/env python3
"""Theorem D and Theorem C': generalising the triangulation past size 3.

Two generalisations of `docs/polygon_diagonals.md`, in opposite directions.

THEOREM C' generalises the CONCLUSION from size 3 to arbitrary size. The
sides of a one-hop class all lie in the real plane, so the Gram matrix
`G = (t_ab)` with `t_aa = r_a^2` has RANK AT MOST 2. That single fact, not the
quadrilateral picture, is what makes triangulation work, and it does not care
about `m`. A class of size `m` has `m - 2` shape parameters after the closure
condition, so `m - 2` rational diagonals pin it down; the ordering that
supplies them is called a RATIONAL CHAIN here.

THEOREM D generalises the HYPOTHESIS, by proving it instead of checking it.
Take two pairs `a`, `b` of a class `(A, B)` whose A-side determinants are one
hop apart, `u_b = u_a - C + E`. Then `v_b = v_a - C + E` as well, because both
sides just swap A for B, so `(u_a, u_b)` and `(v_a, v_b)` are BOTH pairs of the
one-hop class `(C, E)`. Writing out Lemma 1 for each,

    loop of {a,b} in (A,B) = e_{v_a} - e_{u_a} - e_{v_b} + e_{u_b}
    loop of that pair in (C,E) = e_{u_b} - e_{u_a} - e_{v_b} + e_{v_a}

are the SAME integer vector. So the two classes share a loop, and its holonomy
can be read off in whichever class is more convenient. If `(C, E)` has size 2,
Theorem A applies there and gives

    cos(Phi) = (|rho_CE|^2 - p_{u_a} p_{u_b} - p_{v_a} p_{v_b}) / (2 sigma_1 sigma_2 sqrt(M)),
    M = p_{u_a} p_{u_b} p_{v_a} p_{v_b}.

The monomial `M` is exactly `(r_a r_b)^2` for the pair in `(A, B)`, so

    t_ab = s_a s_b r_a r_b cos(Phi) = s_a s_b (|rho_CE|^2 - p_{u_a} p_{u_b} - p_{v_a} p_{v_b}) / (2 sigma_1 sigma_2)

is RATIONAL under (H1) and (H2), with the square roots cancelling identically.
The rational diagonal is therefore not a coincidence: it is Theorem A, read
through a shared loop from a smaller class.

The consequence is an induction on class size. Theorem A settles size <= 2 and
supplies rational diagonals to the larger classes it links; Theorem C' converts
those into multiquadratic fields at any size. The residue is combinatorial and
carries no amplitudes at all: does every large class contain a pair linked by a
class of size at most 2?

Usage:
    python scripts/shared_loop_rationality.py
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
sys.path.insert(0, str(ROOT / "scripts"))

from gpc_census.states import one_hop_classes  # noqa: E402
from polygon_diagonals import (  # noqa: E402
    _degree,
    _mod2,
    _mul_conj,
    exact_parts,
    orient_class,
    transport_class,
)


def load_states():
    for line in (ROOT / "results/data/states.jsonl").read_text().splitlines():
        if line.strip():
            yield json.loads(line)


def _sign(dets, u, v, a_mode, b_mode):
    return ((-1) ** dets[v].index(b_mode)) * ((-1) ** dets[u].index(a_mode))


def linking_channel(dets, classes, ua, ub):
    """The one-hop class containing the unordered determinant pair (ua, ub)."""
    for (c_mode, e_mode), pairs in classes.items():
        for i, j in pairs:
            if {i, j} == {ua, ub}:
                return (c_mode, e_mode), len(pairs), (i, j)
    return None, None, None


def theorem_D_prediction(parts, dets, weights, channel, pairs, ua, ub, va, vb):
    """Theorem A on the linking class, converted into a prediction for t_ab.

    Returns the predicted RATIONAL value of t_ab up to the overall sign
    ``s_a s_b``, which the caller applies. Everything here is exact.
    """
    c_mode, e_mode = channel
    oriented = orient_class(dets, c_mode, e_mode, pairs)
    total = (sp.Integer(0), sp.Integer(0))
    for x, y, sigma in oriented:
        re_, im_ = _mul_conj(parts[x], parts[y])
        total = (total[0] + sigma * re_, total[1] + sigma * im_)
    rho_mod2 = _mod2(total)
    # the two pairs of the size-2 linking class, as determinant pairs
    (x1, y1, sig1), (x2, y2, sig2) = oriented
    bulk = weights[x1] * weights[y1] + weights[x2] * weights[y2]
    value = sp.radsimp(sp.simplify((rho_mod2 - bulk) / (2 * sig1 * sig2)))
    return value, rho_mod2


def analyse(rec, min_size: int) -> list[dict]:
    dets = [tuple(t) for t in rec["closed_form"]["support_dets"]]
    classes = one_hop_classes(dets)
    big = {k: v for k, v in classes.items() if len(v) >= min_size}
    if not big:
        return []
    parts = [exact_parts(sp.sympify(p)) for p in rec["closed_form"]["pretty"]]
    weights = [_mod2(p) for p in parts]

    out = []
    for (a_mode, b_mode), pairs in sorted(big.items()):
        oriented = orient_class(dets, a_mode, b_mode, pairs)
        z = []
        for u, v, s in oriented:
            re_, im_ = _mul_conj(parts[u], parts[v])
            z.append((sp.expand(s * re_), sp.expand(s * im_)))
        r2 = [_mod2(zi) for zi in z]

        links = []
        for i, j in itertools.combinations(range(len(oriented)), 2):
            ua, va, sa = oriented[i]
            ub, vb, sb = oriented[j]
            t_ab = sp.radsimp(sp.simplify(
                sp.expand(z[i][0] * z[j][0] + z[i][1] * z[j][1])))
            degree = _degree(t_ab)
            entry = {
                "pair": [i, j],
                "t_ab": str(t_ab),
                "t_ab_degree": degree,
                "t_ab_rational": None if degree is None else degree == 1,
                "u_hop": len(set(dets[ua]) ^ set(dets[ub])) // 2,
                "v_hop": len(set(dets[va]) ^ set(dets[vb])) // 2,
            }
            if entry["u_hop"] == 1:
                channel, size, _ = linking_channel(dets, classes, ua, ub)
                v_channel, v_size, _ = linking_channel(dets, classes, va, vb)
                entry["linking_channel"] = list(channel) if channel else None
                entry["linking_channel_size"] = size
                entry["v_side_same_channel"] = channel == v_channel
                # the shared-loop identity, checked as integer vectors
                loop_ab = {va: 1, ua: -1, vb: -1, ub: 1}
                loop_ce = {ub: 1, ua: -1, vb: -1, va: 1}
                entry["shared_loop"] = loop_ab == loop_ce
                if size == 2:
                    predicted, rho2 = theorem_D_prediction(
                        parts, dets, weights, channel, classes[channel],
                        ua, ub, va, vb)
                    predicted = sp.radsimp(sp.simplify(sa * sb * predicted))
                    entry["theorem_D_predicted_t_ab"] = str(predicted)
                    entry["theorem_D_predicted_rational"] = _degree(predicted) == 1
                    entry["theorem_D_matches"] = bool(
                        sp.simplify(sp.expand(predicted - t_ab)) == 0)
                    entry["linking_rho_mod2"] = str(rho2)
            links.append(entry)

        # Theorem C': a rational chain is an ordering whose partial sums all
        # have rational squared modulus. Search orderings for one.
        chain = None
        # Every ordering is searched. An earlier version skipped orderings with
        # order[0] > order[-1] to "dedupe reverses", which is wrong: reversing
        # an ordering changes which partial sums appear, so the reverse of a
        # chain is generally not a chain. That bug reported (3,10) v98 as having
        # no rational chain when it has one.
        for order in itertools.permutations(range(len(oriented))):
            partial = (sp.Integer(0), sp.Integer(0))
            moduli, good = [], True
            for k in order:
                partial = (sp.expand(partial[0] + z[k][0]),
                           sp.expand(partial[1] + z[k][1]))
                value = _mod2(partial)
                if _degree(value) != 1:
                    good = False
                    break
                moduli.append(str(value))
            if good:
                chain = {"order": list(order), "partial_moduli": moduli}
                break

        # Theorem C' needs m-2 rational diagonals; Theorem D proves one per
        # pair carrying a size-2 link. When D supplies enough, the class is
        # proved end to end from Theorem A with no checked arithmetic input.
        size = len(oriented)
        proved_diagonals = sum(1 for e in links
                               if e.get("linking_channel_size") == 2
                               and e.get("theorem_D_matches"))
        out.append({
            "system": rec["system"],
            "index": rec["index"],
            "channel": [a_mode, b_mode],
            "class_size": size,
            "diagonals_proved_by_theorem_D": proved_diagonals,
            "diagonals_needed": size - 2,
            "proved_end_to_end": bool(proved_diagonals >= size - 2),
            "exactly_at_the_minimum": bool(proved_diagonals == size - 2),
            "transport_class": transport_class(
                rec["system"], rec["integer_form"], rec["denominator"]),
            "r_squared": [str(v) for v in r2],
            "links": links,
            "rational_chain": chain,
            "has_rational_chain": chain is not None,
            "has_size_two_link": any(e.get("linking_channel_size") == 2
                                     for e in links),
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-class-size", type=int, default=3)
    ap.add_argument("-o", "--out",
                    default=str(ROOT / "results/data/shared_loop_rationality.json"))
    args = ap.parse_args()

    blocks = []
    for rec in load_states():
        found = analyse(rec, args.min_class_size)
        for block in found:
            print(f"{block['system']} v{block['index']} ch{block['channel']} "
                  f"size {block['class_size']} "
                  f"size2link={block['has_size_two_link']} "
                  f"chain={block['has_rational_chain']}", flush=True)
        blocks.extend(found)

    links = [e for b in blocks for e in b["links"]]
    hopped = [e for e in links if e["u_hop"] == 1]
    size_two = [e for e in hopped if e.get("linking_channel_size") == 2]
    checked = [e for e in size_two if "theorem_D_matches" in e]

    payload = {
        "what": "Theorem D (a shared loop makes the diagonal rational, by "
                "Theorem A on a smaller class) and Theorem C' (rank-2 Gram "
                "triangulation at arbitrary class size).",
        "theorems": {
            "D": "if u_b = u_a - C + E then the (A,B) loop of {a,b} equals the "
                 "(C,E) loop of (u_a,u_b); when (C,E) has size 2, Theorem A "
                 "there gives t_ab = s_a s_b (|rho_CE|^2 - p_ua p_ub - p_va "
                 "p_vb) / (2 sigma_1 sigma_2), which is RATIONAL because the "
                 "monomial M = (r_a r_b)^2 cancels identically",
            "C_prime": "the Gram matrix of the sides has rank at most 2 because "
                       "they lie in the real plane, so a rational chain of "
                       "partial sums forces every t_ab into a multiquadratic "
                       "field at any class size, not only size 3",
        },
        "summary": {
            "classes": len(blocks),
            "distinct_states": len({(b["system"], b["index"]) for b in blocks}),
            "transport_classes": len({b["transport_class"] for b in blocks}),
            "class_sizes": {str(k): sum(1 for b in blocks if b["class_size"] == k)
                            for k in sorted({b["class_size"] for b in blocks})},
            "pairs": len(links),
            "one_hop_linked_pairs": len(hopped),
            "shared_loop_holds": sum(1 for e in hopped if e["shared_loop"]),
            "v_side_same_channel": sum(1 for e in hopped
                                       if e["v_side_same_channel"]),
            "size_two_links": len(size_two),
            "theorem_D_checked": len(checked),
            "theorem_D_matches": sum(1 for e in checked if e["theorem_D_matches"]),
            "theorem_D_predicted_rational": sum(
                1 for e in checked if e["theorem_D_predicted_rational"]),
            "size_two_link_and_t_rational": sum(
                1 for e in size_two if e["t_ab_rational"]),
            "classes_with_a_size_two_link": sum(1 for b in blocks
                                                if b["has_size_two_link"]),
            "classes_with_a_rational_chain": sum(1 for b in blocks
                                                 if b["has_rational_chain"]),
            "proved_end_to_end": sum(1 for b in blocks if b["proved_end_to_end"]),
            "exactly_at_the_minimum": sum(1 for b in blocks
                                          if b["exactly_at_the_minimum"]),
            "not_proved_end_to_end": sorted(
                f'{b["system"]} v{b["index"]} ch{tuple(b["channel"])}'
                for b in blocks if not b["proved_end_to_end"]),
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
