"""Standalone replay of the quantization-character result.

Imports NOTHING from `gpc_census`. The wedge action of a permutation is
recomputed through MINORS rather than through a sorting sign: for a linear map
`h` on `C^d`, the matrix of `wedge^n h` in the Slater basis is

    (wedge^n h)_{T,S} = det( h[T, S] ),

so building the `d x d` permutation matrix and taking `n x n` minors is an
independent route to the signs the generator obtains combinatorially.

The character condition is then re-derived and compared against
`results/data/quantization_characters.json`, and the predicted periods are
compared against the measured ones in
`results/data/quantization_multiplicities.json`.

Run: uv run scripts/verify_quantization_character_standalone.py
Exit status 0 means every replayed quantity matched.
"""
from __future__ import annotations

import itertools as it
import json
import pathlib
from math import gcd

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"

failures: list[str] = []


def check(ok: bool, label: str) -> None:
    if not ok:
        failures.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label}")


def exact_state(record):
    n, d = (int(x) for x in record["system"].strip("()").split(","))
    basis = list(it.combinations(range(d), n))
    at = {s: i for i, s in enumerate(basis)}
    psi = [sp.Integer(0)] * len(basis)
    for pretty, det in zip(record["closed_form"]["pretty"], record["closed_form"]["support_dets"]):
        e = sp.sympify(pretty)
        psi[at[tuple(det)]] = sp.expand(sp.simplify(sp.expand_complex(sp.expand_func(e))))
    return n, d, basis, at, psi


def wedge_by_minors(sigma, n, d, basis, at, psi):
    """`wedge^n(P_sigma) psi`, with every sign taken from an n x n minor."""
    p = sp.zeros(d, d)
    for j in range(d):
        p[sigma[j], j] = 1                       # P e_j = e_{sigma(j)}
    out = [sp.Integer(0)] * len(basis)
    for k, s in enumerate(basis):
        if psi[k] == 0:
            continue
        target = tuple(sorted(sigma[o] for o in s))
        minor = p[list(target), list(s)].det()
        out[at[target]] += minor * psi[k]
    return out


def blocks_of(form):
    out = []
    for i, v in enumerate(form):
        if i and v == form[i - 1]:
            out[-1].append(i)
        else:
            out.append([i])
    return out


def cycle_sign(perm):
    seen = [False] * len(perm)
    sign = 1
    for i in range(len(perm)):
        if seen[i]:
            continue
        j, length = i, 0
        while not seen[j]:
            seen[j] = True
            j, length = perm[j], length + 1
        if length % 2 == 0:
            sign = -sign
    return sign


def replay_vertex(record, q, limit=512):
    n, d, basis, at, psi = exact_state(record)
    form = [int(x) for x in record["integer_form"]]
    blocks = blocks_of(form)
    numerators = [form[b[0]] for b in blocks]

    period = q
    found = []
    for choice in it.product(*[list(it.permutations(b)) for b in blocks]):
        sigma = list(range(d))
        for block, image in zip(blocks, choice):
            for src, dst in zip(block, image):
                sigma[src] = dst
        moved = wedge_by_minors(sigma, n, d, basis, at, psi)
        scalar, ok = None, True
        for k in range(len(basis)):
            if psi[k] == 0:
                if sp.simplify(moved[k]) != 0:
                    ok = False
                    break
                continue
            ratio = sp.simplify(moved[k] / psi[k])
            if scalar is None:
                scalar = ratio
            elif sp.simplify(ratio - scalar) != 0:
                ok = False
                break
        if not ok or scalar is None:
            continue
        dets = [cycle_sign([b.index(sigma[o]) for o in b]) for b in blocks]
        least = None
        for k in range(q, limit * q + 1, q):
            value = scalar ** k
            for i, det in enumerate(dets):
                value *= sp.Integer(det) ** (-k * numerators[i] // q)
            if sp.simplify(value - 1) == 0:
                least = k
                break
        if least is None:
            return None, None
        found.append((tuple(sigma), least))
        period = period * least // gcd(period, least)
    return period, found


def main() -> int:
    art = json.loads((DATA / "quantization_characters.json").read_text())
    measured = {(r["system"], r["index"]): r for r in
                json.loads((DATA / "quantization_multiplicities.json").read_text())["rows"]}
    states = {(r["system"], r["index"]): r
              for r in (json.loads(x)
                        for x in (DATA / "states.jsonl").read_text().splitlines())}

    print(f"replaying {len(art['rows'])} vertices through the minor route")
    for row in art["rows"]:
        key = (row["system"], row["index"])
        q = row["denominator"]
        period, found = replay_vertex(states[key], q)
        label = f"{row['system']}#{row['index']} q={q} predicted={row['predicted_period']}"
        check(period == row["predicted_period"], label)
        check(len(found) == row["stabilizer_elements"],
              f"{row['system']}#{row['index']} elements={row['stabilizer_elements']}")
        check(period == measured[key]["quantization_period"],
              f"{row['system']}#{row['index']} predicted == measured period "
              f"({measured[key]['quantization_period']})")

    print(f"\n{len(failures)} failures")
    for f in failures:
        print("  FAIL", f)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
