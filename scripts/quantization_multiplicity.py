"""Quantization multiplicities at census vertices: the representation-theoretic
presentation of the same moment polytope.

For a vertex `lambda` of `P(n,d)` with denominator `q`, and `k` a positive
multiple of `q`, `k*lambda` is a dominant integral weight of `U(d)`. Its
multiplicity in `Sym^k(wedge^n C^d)` is the dimension of the space of sections
of the `k`-th power of the prequantum line bundle over the symplectic reduction
`M_lambda`. So the multiplicity sequence sees `M_lambda` from the
representation side, where the symplectic side is out of reach: a constant
sequence equal to 1 says `M_lambda` is a point, and growth like `k^r` says
`dim_C M_lambda = r`.

Multiplicities are computed by the alternating formula

    mult(mu) = sum_{w in S_d} sgn(w) * P_k(mu + rho - w rho),   rho = (d-1,...,0)

with `P_k(nu)` the weight multiplicity of `Sym^k(wedge^n C^d)`: the number of
multisets of `k` of the `binom(d,n)` occupation vectors summing to `nu`. Both
halves are exact integer counts, so the result is exact.

Run: uv run scripts/quantization_multiplicity.py [--systems 3,6 3,7] [--budget 900]
Writes results/data/quantization_multiplicities.json.
"""
from __future__ import annotations

import argparse
import itertools as it
import json
import pathlib
import time
from functools import lru_cache
from math import gcd

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
OUT = DATA / "quantization_multiplicities.json"


def weight_counter(n: int, d: int):
    """`P_k(nu)` for `Sym^k(wedge^n C^d)`, memoized over the generator suffix."""
    gens = [tuple(1 if o in s else 0 for o in range(d)) for s in it.combinations(range(d), n)]

    @lru_cache(maxsize=None)
    def count(i: int, rest: tuple[int, ...]) -> int:
        total = sum(rest)
        if total == 0:
            return 1
        if total % n:
            return 0
        k_left = total // n
        if i >= len(gens):
            return 0
        if max(rest) > k_left:
            return 0
        # every remaining generator from i on must cover the still-occupied orbitals
        reachable = [0] * d
        for g in gens[i:]:
            for o in range(d):
                reachable[o] += g[o]
        for o in range(d):
            if rest[o] and not reachable[o]:
                return 0
        g = gens[i]
        out = count(i + 1, rest)
        take = tuple(r - x for r, x in zip(rest, g))
        if min(take) >= 0:
            out += count(i, take)
        return out

    def p_k(nu) -> int:
        nu = tuple(int(x) for x in nu)
        if min(nu) < 0:
            return 0
        return count(0, nu)

    return p_k


def multiplicity(mu, n: int, d: int, p_k) -> int:
    """Multiplicity of the irreducible of highest weight `mu` in `Sym^k(...)`."""
    rho = list(range(d - 1, -1, -1))
    total = 0
    for perm in it.permutations(range(d)):
        sign = permutation_sign(perm)
        nu = [mu[i] + rho[i] - rho[perm[i]] for i in range(d)]
        if min(nu) < 0:
            continue
        total += sign * p_k(nu)
    return total


def permutation_sign(perm) -> int:
    seen = [False] * len(perm)
    sign = 1
    for i in range(len(perm)):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            length += 1
        if length % 2 == 0:
            sign = -sign
    return sign


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*", default=["3,6", "3,7"])
    ap.add_argument("--budget", type=float, default=900.0,
                    help="seconds per system before remaining k values are skipped")
    ap.add_argument("--max-k", type=int, default=24)
    args = ap.parse_args()

    results = []
    for spec in args.systems:
        n, d = (int(x) for x in spec.split(","))
        verts = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
        p_k = weight_counter(n, d)
        started = time.time()
        print(f"({n},{d}): {len(verts)} vertices", flush=True)
        rows = [{"system": f"({n},{d})", "index": i,
                 "integer_form": [int(x) for x in v["integer_form"]],
                 "denominator": int(v["denominator"] or 1),
                 "evidence": "exact", "multiplicities": []}
                for i, v in enumerate(verts)]
        # Round-robin over multiples, so a budget cutoff truncates every vertex at
        # the same multiple instead of starving the vertices late in the list.
        for step in range(1, args.max_k + 1):
            live = [r for r in rows if r["denominator"] * step <= args.max_k]
            if not live:
                break
            if time.time() - started > args.budget:
                print(f"  budget reached before multiple {step}", flush=True)
                break
            for r in live:
                k = r["denominator"] * step
                t0 = time.time()
                mu = [x * k // r["denominator"] for x in r["integer_form"]]
                m = multiplicity(mu, n, d, p_k)
                r["multiplicities"].append(
                    {"k": k, "multiplicity": int(m), "secs": round(time.time() - t0, 2)})
                print(f"  v{r['index']} k={k} mult={m} "
                      f"({r['multiplicities'][-1]['secs']}s)", flush=True)
        for r in rows:
            seq = r["multiplicities"]
            nonzero = [s["k"] for s in seq if s["multiplicity"]]
            period = 0
            for k in nonzero:
                period = gcd(period, k)
            r["all_one"] = bool(seq) and all(s["multiplicity"] == 1 for s in seq)
            r["nonzero_all_one"] = bool(nonzero) and all(
                s["multiplicity"] == 1 for s in seq if s["multiplicity"])
            r["quantization_period"] = period or None
            r["period_exceeds_denominator"] = bool(period and period != r["denominator"])
            r["period_consistent"] = bool(period) and all(
                (s["multiplicity"] > 0) == (s["k"] % period == 0) for s in seq)
            r["k_tested"] = [s["k"] for s in seq]
        results.extend(rows)

    tested = [r for r in results if r["k_tested"]]
    periodic = [f"{r['system']}#{r['index']}" for r in tested
                if r["period_exceeds_denominator"]]
    inconsistent = [f"{r['system']}#{r['index']}" for r in tested
                    if not r["period_consistent"]]
    untested = [f"{r['system']}#{r['index']}" for r in results if not r["k_tested"]]
    violations = [f"{r['system']}#{r['index']}" for r in tested if not r["all_one"]]
    zero = [f"{r['system']}#{r['index']}" for r in tested
            if any(s["multiplicity"] == 0 for s in r["multiplicities"])]
    payload = {
        "generated_by": "scripts/quantization_multiplicity.py",
        "prereg": "docs/prereg_symplectic_invariants.md",
        "evidence": "exact",
        "note": (
            "Multiplicity of highest weight k*lambda in Sym^k(wedge^n C^d) at census "
            "vertices, by the alternating Weyl formula over exact integer weight "
            "counts. A vertex whose multiplicity is 1 at every tested k has a "
            "zero-dimensional reduced space as far as the tested k can see; this is "
            "evidence, not a proof for all k."
        ),
        "rows": results,
        "scorecard": {
            "P6_vertices_are_multiplicity_free": {
                "verdict": "PASS" if (tested and not violations and not zero) else "FAIL",
                "scope": len(tested),
                "independent_scope": len(tested),
                "independent_distinct_states": len(tested),
                "vertices_untested_budget": untested,
                "violations": violations,
                "zero_multiplicity_vertices": zero,
            },
            "finding_quantization_period": {
                "note": (
                    "P6 was written assuming every admissible k realizes the vertex. "
                    "Where it fails, the failure is periodic, not erratic: the "
                    "multiplicity is 1 on multiples of a period p and 0 elsewhere, "
                    "with p a multiple of the spectrum denominator. p is a lattice "
                    "invariant of the vertex that the polytope does not record."
                ),
                "vertices_with_period_above_denominator": periodic,
                "period_model_inconsistent": inconsistent,
                "nonzero_multiplicities_all_one": all(
                    r["nonzero_all_one"] for r in tested),
            },
        },
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    print("P6:", payload["scorecard"]["P6_vertices_are_multiplicity_free"]["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
