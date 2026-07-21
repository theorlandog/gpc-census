#!/usr/bin/env python3
"""Polygon-target solver: exact phases for a fixed skeleton with mixed one-hop
targets (prescribed magnitudes in-block, zero off-block).

See docs/RESEARCH.md, "v96 campaign" (remaining rung (a), the hybrid family).
Given a support of determinants with fixed integer weights k_t (so moduli
|c_t| = sqrt(k_t/den)), find per-determinant phases theta_t so that every
one-hop pair class {p,q} meets its target:

    P_pq = sum_{A^B = {p,q}} sigma_AB sqrt(k_A k_B)/den exp(i(theta_A - theta_B))
    |P_pq|^2 = tsq_pq        (tsq_pq > 0: in-block Schur-Horn magnitude;
                              tsq_pq = 0: off-block, the class must CANCEL)

The magnitude and cancellation constraints are geometrically different, and the
solver treats them as such:

  MAGNITUDE (tsq > 0) is ONE real equation (the off-diagonal's own phase is free
    -- gauge inside the block). Reduced to a single unknown phase it is a law of
    cosines step: cos(theta_f - c) is pinned, theta_f = c +- acos(...), two
    reflection branches. This is the stage-3b off-diagonal exactifier
    (gpc_census.exactify.exactify_interference) as the single-class case.
  CANCELLATION (tsq = 0) is TWO real equations (Re P = Im P = 0). A 1-term class
    cannot cancel (one nonzero vector) -- the LONE-PAIR FUNNEL, an exact
    infeasibility prune (identical to the P1 prune in signed_design_fast.py). A
    2-term class cancels only if its two moduli match, then its phase difference
    is pinned exactly (no branch). A 3-term class is a RIGID TRIANGLE: the three
    fixed-length sides close, so the relative phases are pinned by the law of
    cosines up to reflection (two branches).

SAFETY: the solver only ever PROPOSES exact (symbolic) phase assignments; each
full assignment is accepted only if gpc_census.exactify.verify_exact certifies
it (gauge-invariant 1-RDM characteristic-polynomial identity). An incomplete
propagation can therefore only MISS a solution, never certify a wrong one.

SCOPE: classes of arity <= 3 are solved completely (this covers every class that
occurs in the 142 certified census interference states -- max arity 3). Classes
of arity >= 4 are polygons with internal freedom; they are handled by the same
bounded branch search but are not proven complete, and are flagged.

Usage:
  # auto targets (degenerate-pair Schur-Horn, reproduces the census exactifier):
  python scripts/polygon_target.py --recertify 3 10        # re-solve (3,10) inter.
  python scripts/polygon_target.py --recertify-all         # all 9 systems
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gpc_census.exactify import (  # noqa: E402
    _active_offdiagonals,
    verify_exact,
)


def one_hop_classes(d, dets):
    """All one-hop pair classes: {(p, q): [(i, j, sign), ...]}.

    Term (i, j, sign) means determinant j carries mode p, determinant i carries
    mode q (i is j with p replaced by q); sign = (-1)^pos(p in j) (-1)^pos(q in i),
    the Slater sign of that off-diagonal contribution. Same convention as
    gpc_census.exactify._active_offdiagonals.
    """
    dets = [tuple(t) for t in dets]
    index = {t: i for i, t in enumerate(dets)}
    classes = {}
    for p in range(d):
        for q in range(p + 1, d):
            terms = []
            for j, t in enumerate(dets):
                if p in t and q not in t:
                    s1 = (-1) ** t.index(p)
                    tp = tuple(sorted(tuple(x for x in t if x != p) + (q,)))
                    if tp in index:
                        i = index[tp]
                        s2 = (-1) ** tp.index(q)
                        terms.append((i, j, s1 * s2))
            if terms:
                classes[(p, q)] = terms
    return classes


def auto_targets(d, dets, ks, den, spectrum):
    """Schur-Horn targets for the degenerate-pair case (the census interference
    regime), via the vetted gpc_census.exactify._active_offdiagonals. Returns
    {(p, q): tsq} for every one-hop class; pairs the splitter leaves unassigned
    get target 0 (must cancel)."""
    import sympy as sp

    _occ, edges = _active_offdiagonals(d, [tuple(t) for t in dets], ks, den,
                                       spectrum)
    tgt = {}
    for (p, q, tsq, terms) in edges:
        tgt[(p, q)] = sp.Integer(0) if tsq is None else tsq
    return tgt


def _gmod(ks, den, i, j):
    import sympy as sp

    return sp.sqrt(sp.Rational(ks[i], den) * sp.Rational(ks[j], den))


def _class_value(theta, ks, den, terms):
    """Symbolic P for a class given a full (no None) theta assignment."""
    import sympy as sp

    P = sum(s * _gmod(ks, den, i, j) * sp.exp(sp.I * (theta[i] - theta[j]))
            for (i, j, s) in terms)
    return sp.expand_complex(P)


def _feasible_lone(ks, den, classes, targets):
    """Lone-pair funnel: the one rigorous 1-term prune.

    A cancellation class (tsq = 0) with a single term is a lone nonzero vector
    that cannot sum to zero, so the skeleton is infeasible -- identical to the
    P1 prune in signed_design_fast.py. A 1-term class with a NONZERO target is
    a lone in-block pair whose magnitude is fixed by its single term; that is
    not a phase constraint and not a hard reject here (verify_exact, which tests
    the spectrum rather than a chosen per-edge magnitude, is the real judge).
    Returns False only when a genuine cancellation lone pair is present.
    """
    for pq, terms in classes.items():
        if len(terms) == 1 and targets.get(pq, 0) == 0:
            return False
    return True


def _branch_candidates(ks, den, classes, targets):
    """Exact candidate phases for the bounded branch search, per class.

    For each class collect the relative-phase angles its geometry admits:
      - magnitude 2/3-term: cos(delta) pinned by law of cosines -> +-acos.
      - cancellation 2-term: antiparallel -> delta in {0, pi}.
      - cancellation 3-term: rigid triangle -> interior angles by law of cosines.
    These seed a joint search over the determinants that carry the coupling.
    """
    import sympy as sp

    base = [sp.Integer(0), sp.pi, sp.pi / 3, -sp.pi / 3,
            2 * sp.pi / 3, -2 * sp.pi / 3, sp.pi / 2, -sp.pi / 2]
    cand = list(base)
    for pq, terms in classes.items():
        tsq = targets.get(pq, sp.Integer(0))
        gs = [_gmod(ks, den, i, j) for (i, j, _s) in terms]
        for a in range(len(gs)):
            for b in range(len(gs)):
                if a == b:
                    continue
                ga, gb = gs[a], gs[b]
                # law of cosines against every other side length (incl. sqrt(tsq))
                thirds = [gs[c] for c in range(len(gs)) if c not in (a, b)]
                if tsq != 0:
                    thirds.append(sp.sqrt(tsq))
                for gc in thirds:
                    denom = 2 * ga * gb
                    if denom == 0:
                        continue
                    r = sp.simplify((gc ** 2 - ga ** 2 - gb ** 2) / denom)
                    if sp.Abs(sp.N(r)) <= 1:
                        cand += [sp.acos(r), -sp.acos(r),
                                 sp.pi - sp.acos(r), -(sp.pi - sp.acos(r))]
    # dedupe preserving order
    return list(dict.fromkeys(cand))


def solve(N, d, spectrum, dets, ks, den, targets=None, max_seeds=3,
          verbose=False):
    """Return an EXACT certified record for the skeleton, or None.

    spectrum: iterable of rationals/strings (the vertex occupation numbers).
    dets: support determinants (N-subsets). ks: integer weights. den: common
    denominator with |c_t|^2 = k_t/den. targets: {(p, q): tsq} magnitude^2 per
    one-hop class; None auto-derives the degenerate-pair Schur-Horn targets.
    """
    import sympy as sp

    dets = [tuple(t) for t in dets]
    S = len(dets)
    spec = [str(s) for s in spectrum]
    classes = one_hop_classes(d, dets)
    if targets is None:
        targets = auto_targets(d, dets, ks, den, spectrum)
    targets = {pq: sp.sympify(v) for pq, v in targets.items()}

    if not _feasible_lone(ks, den, classes, targets):
        if verbose:
            print("  infeasible: lone-pair funnel")
        return None

    # classes that still need a phase solve: arity >= 2 (1-term already gated)
    work = [(pq, terms) for pq, terms in classes.items() if len(terms) >= 2]

    # -- pass 1: single-free-determinant propagation (magnitude + cancellation) --
    solve_det = {}
    used = set()
    # order classes so those whose determinants are least shared resolve first
    from collections import Counter
    cdets = {pq: {x for (i, j, _s) in terms for x in (i, j)}
             for pq, terms in work}
    glob = Counter(x for s in cdets.values() for x in s)
    order = sorted(range(len(work)),
                   key=lambda e: min((glob[x] for x in cdets[work[e][0]]),
                                     default=99))
    ok_prop = True
    for e in order:
        pq, terms = work[e]
        pick = next((x for x in sorted({x for (i, j, _s) in terms
                                        for x in (i, j)}, key=lambda z: glob[z])
                     if x not in used
                     and sum(x in (i, j) for (i, j, _s) in terms) == 1), None)
        if pick is None:
            ok_prop = False
            break
        solve_det[pq] = pick
        used.add(pick)

    if ok_prop:
        rec = _propagate(N, d, spec, dets, ks, den, classes, targets,
                         solve_det, verbose)
        if rec is not None:
            return rec

    # -- pass 2: bounded exact branch search over the coupling determinants --
    seeds = sorted(used) if used else list(range(min(S, max_seeds)))
    if not seeds or len(seeds) > max_seeds:
        seeds = seeds[:max_seeds]
    cand = _branch_candidates(ks, den, classes, targets)
    if verbose:
        print(f"  branch search: {len(seeds)} seeds x {len(cand)} candidates")
    for combo in itertools.product(cand, repeat=len(seeds)):
        th = [sp.Integer(0)] * S
        for dd, val in zip(seeds, combo):
            th[dd] = val
        amps = [sp.sqrt(sp.Rational(ks[i], den)) * sp.exp(sp.I * th[i])
                for i in range(S)]
        if verify_exact(N, d, spectrum, dets, amps):
            return _record(ks, den, amps)
    return None


def _propagate(N, d, spec, dets, ks, den, classes, targets, solve_det, verbose):
    import sympy as sp

    S = len(dets)
    theta = [None if i in solve_det.values() else sp.Integer(0)
             for i in range(S)]
    work = [(pq, classes[pq]) for pq in solve_det]
    for _ in range(len(work) + 2):
        progressed = False
        for pq, terms in work:
            f = solve_det[pq]
            if theta[f] is not None:
                continue
            # every determinant in the class except f (including the free term's
            # other endpoint) must be fixed before f can be solved
            if any(theta[x] is None
                   for (i, j, _s) in terms for x in (i, j) if x != f):
                continue
            tsq = targets.get(pq, sp.Integer(0))
            P0 = sum(s * _gmod(ks, den, i, j) * sp.exp(sp.I * (theta[i] - theta[j]))
                     for (i, j, s) in terms if f not in (i, j))
            P0 = sp.expand_complex(P0)
            (i, j, s) = next((i, j, s) for (i, j, s) in terms if f in (i, j))
            g = _gmod(ks, den, i, j)
            other = j if f == i else i
            if tsq == 0:
                # cancellation: the free term must equal -P0, which needs the
                # fixed part's magnitude to match the free side (|P0| == g)
                if sp.simplify(sp.Abs(P0) ** 2 - g ** 2) != 0:
                    return None
                theta[f] = _pin_cancellation(theta, s, g, P0, i, j, f, other)
                if theta[f] is None:
                    return None
            else:
                X = sp.sqrt(tsq)
                Pabs2 = sp.expand_complex(sp.Abs(P0) ** 2)
                rhs = sp.simplify((X ** 2 - Pabs2 - g ** 2) / (2 * s * g))
                if P0 == 0:
                    phi = sp.Integer(0)
                else:
                    if sp.Abs(sp.N(rhs / sp.Abs(P0))) > 1:
                        return None
                    phi = sp.arg(P0) + sp.acos(sp.simplify(rhs / sp.Abs(P0)))
                theta[f] = sp.simplify((theta[other] + phi) if f == i
                                       else (theta[other] - phi))
            progressed = True
        if all(t is not None for t in theta) or not progressed:
            break
    theta = [sp.Integer(0) if t is None else t for t in theta]
    amps = [sp.sqrt(sp.Rational(ks[i], den)) * sp.exp(sp.I * theta[i])
            for i in range(S)]
    if verify_exact(N, d, [sp.Rational(Fraction(x).numerator, Fraction(x).denominator)
                           for x in spec], dets, amps):
        return _record(ks, den, amps)
    return None


def _pin_cancellation(theta, s, g, P0, i, j, f, other):
    """Set free phase so s*g*exp(i(theta_i - theta_j)) = -P0 exactly."""
    import sympy as sp

    if P0 == 0:
        return theta[other]
    want = sp.arg(-P0) - (sp.Integer(0) if s > 0 else sp.pi)
    # want = theta_i - theta_j ; theta[other] is the fixed one
    if f == i:
        return sp.simplify(theta[other] + want)
    return sp.simplify(theta[other] - want)


def _record(ks, den, amps):
    import sympy as sp

    return {"status": "EXACT", "weights": list(ks), "den": den,
            "amplitudes": [sp.srepr(a) for a in amps],
            "pretty": [str(sp.simplify(a)) for a in amps]}


# ------------------------------------------------------------------ CLI / recert
def _load_interference(system=None):
    out = []
    with open(ROOT / "results" / "data" / "states.jsonl") as f:
        for line in f:
            dd = json.loads(line)
            if dd["status"] != "OK" or dd["classified"] != "INTERFERENCE":
                continue
            if system is not None and dd["system"] != system:
                continue
            out.append(dd)
    return out


def recertify(N, d, limit=None):
    system = f"({N},{d})"
    states = _load_interference(system)
    if limit:
        states = states[:limit]
    ok = 0
    for st in states:
        cf = st["closed_form"]
        den = st["denominator"]
        spec = [Fraction(x, den) for x in st["integer_form"]]
        rec = solve(N, d, spec, cf["support_dets"], cf["weights"], den)
        good = rec is not None and rec["weights"] == cf["weights"]
        ok += good
        if not good:
            print(f"  MISS ({N},{d}) v{st['index']}")
    print(f"({N},{d}) interference re-solved: {ok}/{len(states)}")
    return ok, len(states)


def recertify_all():
    systems = sorted({s["system"] for s in _load_interference()},
                     key=lambda t: eval(t))
    tot_ok = tot = 0
    for sysstr in systems:
        N, d = eval(sysstr)
        ok, n = recertify(N, d)
        tot_ok += ok
        tot += n
    print(f"\nTOTAL interference re-solved by polygon-target solver: "
          f"{tot_ok}/{tot}")
    return tot_ok == tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recertify", nargs=2, type=int, metavar=("N", "D"))
    ap.add_argument("--recertify-all", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    if a.recertify_all:
        raise SystemExit(0 if recertify_all() else 1)
    if a.recertify:
        ok, n = recertify(*a.recertify, limit=a.limit)
        raise SystemExit(0 if ok == n else 1)
    ap.print_help()


if __name__ == "__main__":
    main()
