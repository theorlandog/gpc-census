#!/usr/bin/env python3
"""Exact support-geometry audit for fermionic determinant carriers.

No project repository imports are used.  The script tests:
  * singleton-incidence rank and non-toric phase complexity;
  * the matroid basis-exchange axiom;
  * Schubert/Gale order-ideal compatibility;
  * cyclic-window families.

Definitions use zero-based orbital labels.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path

import sympy as sp

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "results/data/schubert_bruhat_support_audit.json"

def all_determinants(d: int, n: int):
    return [tuple(x) for x in itertools.combinations(range(d), n)]

def cyclic_windows(d: int, n: int):
    return [tuple(sorted((k+j) % d for j in range(n))) for k in range(d)]

def singleton_matrix(support, d):
    return sp.Matrix([[int(p in D) for D in support] for p in range(d)])

def support_invariants(support, d):
    A = singleton_matrix(support, d)
    rank = int(A.rank())
    m = len(support)
    # All columns lie in sum x_i=N, hence affine dimension = rank(A)-1.
    affine_dim = rank - 1
    complexity = (m - 1) - affine_dim
    return {
        "carrier_size": m,
        "singleton_rank": rank,
        "support_polytope_affine_dimension": affine_dim,
        "non_toric_phase_complexity": complexity,
    }

def is_matroid_basis_set(support):
    bases = {frozenset(D) for D in support}
    if not bases:
        return False, None
    r = len(next(iter(bases)))
    if any(len(B) != r for B in bases):
        return False, {"reason": "nonuniform cardinality"}
    for A in sorted(bases, key=lambda x: tuple(sorted(x))):
        for B in sorted(bases, key=lambda x: tuple(sorted(x))):
            for x in sorted(A-B):
                witnesses = [
                    y for y in sorted(B-A)
                    if frozenset((A-{x})|{y}) in bases
                ]
                if not witnesses:
                    return False, {
                        "A": sorted(A), "B": sorted(B), "exchange_element": x
                    }
    return True, None

def gale_leq(I, J):
    return all(i <= j for i,j in zip(sorted(I), sorted(J)))

def is_gale_order_ideal(support, d, n):
    S = {tuple(D) for D in support}
    universe = all_determinants(d,n)
    for J in S:
        for I in universe:
            if gale_leq(I,J) and I not in S:
                return False, {"included": list(J), "missing_lower": list(I)}
    return True, None

def cyclic_audit(max_d=14):
    rows=[]
    for d in range(3,max_d+1):
        for n in range(1,d):
            S=cyclic_windows(d,n)
            inv=support_invariants(S,d)
            matroid,witness=is_matroid_basis_set(S)
            schubert,sw=is_gale_order_ideal(S,d,n)
            rows.append({
                "d":d,"N":n,
                **inv,
                "predicted_complexity_gcd_minus_one": math.gcd(d,n)-1,
                "matroid_basis_set":matroid,
                "matroid_failure_witness":witness,
                "gale_order_ideal":schubert,
                "gale_failure_witness":sw,
            })
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--max-d",type=int,default=14)
    ap.add_argument("--out",type=Path,default=DEFAULT_OUT)
    args=ap.parse_args()
    rows=cyclic_audit(args.max_d)
    payload={
        "schema_version":1,
        "cyclic_windows":rows,
        "summary":{
            "all_complexities_match_gcd_minus_one":
                all(r["non_toric_phase_complexity"] ==
                    r["predicted_complexity_gcd_minus_one"] for r in rows),
            "nontrivial_matroid_cases":[
                [r["d"],r["N"]] for r in rows
                if r["matroid_basis_set"] and 1 < r["N"] < r["d"]-1
            ],
            "nontrivial_gale_order_ideal_cases":[
                [r["d"],r["N"]] for r in rows
                if r["gale_order_ideal"] and 1 < r["N"] < r["d"]-1
            ],
        }
    }
    text=json.dumps(payload,indent=2)+"\n"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(f"wrote {args.out}: {len(rows)} cyclic-window rows")

if __name__=="__main__":
    main()
