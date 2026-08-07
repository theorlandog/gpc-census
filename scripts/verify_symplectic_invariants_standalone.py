"""Standalone replay of the cross-realization invariant census.

Imports NOTHING from `gpc_census`. Every quantity is recomputed from the
shipped data files by an independently written route, and compared against
`results/data/symplectic_invariants.json` and
`results/data/quantization_multiplicities.json`.

The routes are deliberately different from the generators:

  stabilizer   rebuilt from scratch from `closed_form.pretty`, exact rank of
               the real span of `d(pi)(X) psi` over `u(d)`;
  toral        solved as the nullspace of `[A | -1]` acting on `(phi, c)`,
               where the generator instead takes the rank of the row
               differences of `A`;
  cone index   edge rays derived from the H-representation by inverting the
               active facet normals at simple vertices, where the generator
               instead reads neighbours off the V-representation;
  multiplicity recomputed with a dictionary generating-function DP over
               weights instead of a memoized suffix DFS, and cross-checked
               against the Weyl dimension identity
               `sum_mu m(mu) dim V_mu = dim Sym^k(wedge^n C^d)`.

Run: uv run scripts/verify_symplectic_invariants_standalone.py [--sample 40]
Exit status 0 means every replayed quantity matched.
"""
from __future__ import annotations

import argparse
import itertools as it
import json
import pathlib
import random
from fractions import Fraction as F
from math import comb, gcd

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
CONSTRAINTS = ROOT / "src" / "gpc_census" / "data" / "constraints.json"

failures: list[str] = []
skipped: list[str] = []


def check(ok: bool, label: str) -> None:
    if not ok:
        failures.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label}")


def check_or_skip(got, want, label: str) -> None:
    """A replay route that declines to answer is a SKIP, never a pass or a fail."""
    if got is None:
        skipped.append(label)
        print(f"  skip {label} (H-rep route declined)")
        return
    check(got == want, label)


# ---------------------------------------------------------------- stabilizer


def replay_stabilizer(record):
    n, d = (int(x) for x in record["system"].strip("()").split(","))
    basis = list(it.combinations(range(d), n))
    at = {s: i for i, s in enumerate(basis)}
    psi = [sp.Integer(0)] * len(basis)
    for pretty, det in zip(record["closed_form"]["pretty"], record["closed_form"]["support_dets"]):
        e = sp.sympify(pretty)
        psi[at[tuple(det)]] = sp.expand(sp.simplify(sp.expand_complex(sp.expand_func(e))))

    def hop(p, q, s):
        if q not in s:
            return None
        i = s.index(q)
        sign = (-1) ** i
        rest = s[:i] + s[i + 1:]
        if p in rest:
            return None
        sign *= (-1) ** sum(1 for o in rest if o < p)
        return sign, tuple(sorted(rest + (p,)))

    col = {}
    for p in range(d):
        for q in range(d):
            v = [sp.Integer(0)] * len(basis)
            for k, s in enumerate(basis):
                if psi[k] == 0:
                    continue
                out = hop(p, q, s)
                if out:
                    sign, tgt = out
                    v[at[tgt]] += sign * psi[k]
            col[(p, q)] = v

    def real_rank(vecs):
        return sp.Matrix([[sp.re(x) for x in v] + [sp.im(x) for x in v] for v in vecs]).rank()

    torus = [[sp.I * x for x in col[(p, p)]] for p in range(d)]
    gens = list(torus)
    for p in range(d):
        for q in range(p + 1, d):
            gens.append([a - b for a, b in zip(col[(p, q)], col[(q, p)])])
            gens.append([sp.I * (a + b) for a, b in zip(col[(p, q)], col[(q, p)])])
    orbit = real_rank(gens) - 1
    return {"orbit_dim": int(orbit), "stab_dim": int(d * d - orbit),
            "toral_stab_dim": int(d - (real_rank(torus) - 1))}


def replay_toral_combinatorial(support_dets, d):
    """Nullspace of `(phi, c) -> A phi - c*1`, an independent route to `t_phi`."""
    rows = []
    for det in support_dets:
        row = [0] * d + [-1]
        for o in det:
            row[o] += 1
        rows.append(row)
    return len(sp.Matrix(rows).nullspace())


# ---------------------------------------------------------------- cone index


def hrep_rows(n, d):
    table = json.loads(CONSTRAINTS.read_text())
    key = f"{n},{d}"
    if key in table:
        sysd, flip = table[key], False
    else:
        sysd, flip = table[f"{d - n},{d}"], True

    def conv(q):
        a, b = [F(x) for x in q["coeffs"]], F(q["rhs"])
        if flip:
            a, b = [-a[d - 1 - i] for i in range(d)], b - sum(a)
        return a, b

    ineq = [conv(q) for q in sysd["inequalities"]]
    eq = [conv(q) for q in sysd["equalities"]]
    for i in range(d - 1):
        c = [F(0)] * d
        c[i], c[i + 1] = F(-1), F(1)
        ineq.append((c, F(0)))
    c = [F(0)] * d
    c[d - 1] = F(-1)
    ineq.append((c, F(0)))
    eq.append(([F(1)] * d, F(n)))
    return ineq, eq


def _primitive(vec):
    """Scale a rational vector to the primitive integer vector on its ray."""
    lcm = 1
    for x in vec:
        lcm = sp.ilcm(lcm, sp.Rational(x).q)
    ints = [int(x * lcm) for x in vec]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    return tuple(x // g for x in ints) if g else tuple(ints)


def replay_cone_index(row, points, lattice):
    """Edge rays at a simple vertex, from the H-rep: invert the active normals.

    The active normals must be the FACET normals themselves, deduplicated by
    direction. Replacing them with any basis of their span (an echelon
    `rowspace`, say) changes the dual basis and therefore the lattice index; an
    earlier version of this function did exactly that and disagreed with the
    generator on every simple vertex, with the telltale small factorial-looking
    determinants of an echelon basis.

    Deduplication is done on the H-rep side alone, by primitive direction in
    lattice coordinates, so this stays independent of the generator's V-rep
    facet test. If the surviving count is not `dim` the vertex is skipped rather
    than guessed at.
    """
    n, d, dim = row["n"], row["d"], row["dim"]
    v = points[row["index"]]
    ineq, _ = hrep_rows(n, d)
    lat = sp.Matrix(lattice).T                      # d x dim
    seen: dict[tuple, None] = {}
    for a, b in ineq:
        if sum(x * y for x, y in zip(a, v)) != b:
            continue
        restricted = [sum(a[i] * lat[i, j] for i in range(d)) for j in range(dim)]
        if all(x == 0 for x in restricted):
            continue
        seen.setdefault(_primitive(restricted), None)
    normals = list(seen)
    if len(normals) != dim:
        return None
    m = sp.Matrix(normals)
    if m.det() == 0:
        return None
    rays = []
    for j in range(dim):
        e = sp.zeros(dim, 1)
        e[j] = -1
        rays.append(list(_primitive(m.solve(e))))
    return abs(int(sp.Matrix(rays).det()))


# -------------------------------------------------------------- multiplicity


def weight_table(n, d, k):
    """All weight multiplicities of `Sym^k(wedge^n C^d)`, by generating-function DP."""
    gens = [tuple(1 if o in s else 0 for o in range(d)) for s in it.combinations(range(d), n)]
    cur = {(0,) * d: {0: 1}}
    for g in gens:
        nxt: dict[tuple[int, ...], dict[int, int]] = {}
        for w, counts in cur.items():
            for used, mult in counts.items():
                t, u, m = w, used, mult
                while u <= k:
                    nxt.setdefault(t, {})
                    nxt[t][u] = nxt[t].get(u, 0) + m
                    t = tuple(a + b for a, b in zip(t, g))
                    u += 1
                    if max(t) > k:
                        break
        cur = nxt
    return {w: c.get(k, 0) for w, c in cur.items() if c.get(k, 0)}


def weyl_dim(mu, d):
    num = den = 1
    for i in range(d):
        for j in range(i + 1, d):
            num *= mu[i] - mu[j] + j - i
            den *= j - i
    return num // den


def replay_multiplicity(mu, n, d, k, table):
    rho = list(range(d - 1, -1, -1))
    total = 0
    for perm in it.permutations(range(d)):
        seen, sign = [False] * d, 1
        for i in range(d):
            if seen[i]:
                continue
            j, ln = i, 0
            while not seen[j]:
                seen[j] = True
                j, ln = perm[j], ln + 1
            if ln % 2 == 0:
                sign = -sign
        nu = tuple(mu[i] + rho[i] - rho[perm[i]] for i in range(d))
        if min(nu) < 0:
            continue
        total += sign * table.get(nu, 0)
    return total


# --------------------------------------------------------------------- main


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20260807)
    args = ap.parse_args()

    art = json.loads((DATA / "symplectic_invariants.json").read_text())
    rows = art["rows"]
    states = {(r["system"], r["index"]): r
              for r in (json.loads(x) for x in (DATA / "states.jsonl").read_text().splitlines())}

    slater = [r for r in rows if set(r["integer_form"]) <= {0, r["denominator"] or 1}]
    rng = random.Random(args.seed)
    sample = slater + rng.sample([r for r in rows if r not in slater],
                                 min(args.sample, len(rows) - len(slater)))

    print(f"stabilizer replay on {len(sample)} rows "
          f"({len(slater)} Slater vertices, {len(sample) - len(slater)} sampled)")
    for row in sample:
        got = replay_stabilizer(states[(row["system"], row["index"])])
        want = row["stabilizer"]
        label = f"{row['system']}#{row['index']} stab={want['stab_dim']}"
        check(got["stab_dim"] == want["stab_dim"]
              and got["orbit_dim"] == want["orbit_dim"]
              and got["toral_stab_dim"] == want["toral_stab_dim"], label)

    print("\ntoral prediction, independent route")
    for row in sample:
        dets = states[(row["system"], row["index"])]["closed_form"]["support_dets"]
        got = replay_toral_combinatorial(dets, row["d"])
        check(got == row["toral_stab_dim_predicted"],
              f"{row['system']}#{row['index']} t_phi={got}")

    print("\ntangent cone lattice index at simple vertices, H-rep route")
    simple = [r for r in rows if r["cone"]["simplicial_tangent_cone"]]
    for row in simple[:args.sample]:
        verts = json.loads((DATA / "vertices" /
                            f"vertices_{row['n']}_{row['d']}.json").read_text())
        pts = [[F(x) for x in v["spectrum"]] for v in verts]
        lattice = lattice_basis(pts, row["d"])
        got = replay_cone_index(row, pts, lattice)
        check_or_skip(got, row["cone"]["lattice_index"],
                      f"{row['system']}#{row['index']} "
                      f"index={row['cone']['lattice_index']}")

    print("\nquantization multiplicities, independent DP")
    qpath = DATA / "quantization_multiplicities.json"
    if qpath.exists():
        q = json.loads(qpath.read_text())
        for r in q["rows"]:
            n, d = (int(x) for x in r["system"].strip("()").split(","))
            for entry in r["multiplicities"]:
                k = entry["k"]
                if k > 8:
                    continue
                table = weight_table(n, d, k)
                mu = [x * k // r["denominator"] for x in r["integer_form"]]
                got = replay_multiplicity(tuple(mu), n, d, k, table)
                check(got == entry["multiplicity"],
                      f"{r['system']}#{r['index']} k={k} mult={entry['multiplicity']}")
        for n, d, k in [(3, 6, 1), (3, 6, 2), (2, 4, 2), (2, 4, 3)]:
            table = weight_table(n, d, k)
            total = 0
            for mu in dominant_weights(d, n * k, k):
                m = replay_multiplicity(mu, n, d, k, table)
                if m:
                    total += m * weyl_dim(mu, d)
            check(total == comb(comb(d, n) + k - 1, k),
                  f"Weyl dimension identity ({n},{d}) k={k}: {total}")
    else:
        print("  (no quantization artifact present, skipped)")

    print(f"\n{len(failures)} failures, {len(skipped)} skipped")
    for f in failures:
        print("  FAIL", f)
    return 1 if failures else 0


def dominant_weights(d, total, cap):
    def rec(i, left, ceiling):
        if i == d - 1:
            if 0 <= left <= min(ceiling, cap):
                yield (left,)
            return
        for x in range(min(ceiling, cap, left), -1, -1):
            if left - x > x * (d - 1 - i):
                break
            yield from ((x,) + rest for rest in rec(i + 1, left - x, x))
    return rec(0, total, cap)


def lattice_basis(points, d):
    base = points[0]
    diffs = sp.Matrix([[a - b for a, b in zip(p, base)] for p in points[1:]])
    rows = []
    for v in diffs.nullspace():
        lcm = 1
        for x in v:
            lcm = sp.ilcm(lcm, sp.Rational(x).q)
        rows.append([int(x * lcm) for x in v])
    if not rows:
        return [[1 if i == j else 0 for j in range(d)] for i in range(d)]
    m = [list(r) for r in rows]
    v = [[1 if i == j else 0 for j in range(d)] for i in range(d)]

    def col_add(t, s, f):
        for r in m:
            r[t] += f * r[s]
        for r in v:
            r[t] += f * r[s]

    def col_swap(a, b):
        for r in m:
            r[a], r[b] = r[b], r[a]
        for r in v:
            r[a], r[b] = r[b], r[a]

    pivot = 0
    for r in range(len(m)):
        nz = [c for c in range(pivot, d) if m[r][c] != 0]
        if not nz:
            continue
        col_swap(pivot, nz[0])
        while True:
            nz = [c for c in range(pivot + 1, d) if m[r][c] != 0]
            if not nz:
                break
            for c in nz:
                col_add(c, pivot, -(m[r][c] // m[r][pivot]))
            nz = [c for c in range(pivot + 1, d) if m[r][c] != 0]
            if nz:
                col_swap(pivot, nz[0])
        pivot += 1
        if pivot == d:
            break
    return [[v[i][c] for i in range(d)] for c in range(pivot, d)]


if __name__ == "__main__":
    raise SystemExit(main())
