"""Access to the precomputed census dataset shipped with the package.

The engine computes constraints, vertices, classifications, and closed-form
states; the results are stored under ``results/data`` (the single source of
truth) and embedded into the wheel at build time as package data
(``gpc_census/_data``). This module reads that data so the results are usable
as a library without recomputing, and falls back to the repository's
``results/data`` when running from a source checkout (where the embedded copy
does not exist).

Everything here is read-only precomputed output. To recompute or extend it, use
the engine directly (``gpc_census.constraints``, ``gpc_census.polytope``,
``gpc_census.classify``, ``gpc_census.states.certify_state``).
"""
from __future__ import annotations

import json
import pathlib
from functools import lru_cache


@lru_cache(maxsize=1)
def data_root() -> pathlib.Path:
    """Directory holding the precomputed dataset (embedded copy, else repo)."""
    try:
        from importlib.resources import files
        embedded = pathlib.Path(str(files("gpc_census") / "_data"))
        if embedded.is_dir():
            return embedded
    except (ModuleNotFoundError, FileNotFoundError, NotADirectoryError, TypeError):
        pass
    return pathlib.Path(__file__).resolve().parents[2] / "results" / "data"


def available_systems() -> list[tuple[int, int]]:
    """(n, d) pairs with a precomputed vertex enumeration on disk, sorted."""
    out = []
    vdir = data_root() / "vertices"
    if vdir.is_dir():
        for p in vdir.glob("vertices_*_*.json"):
            _, n, d = p.stem.split("_")
            out.append((int(n), int(d)))
    return sorted(out)


def vertices(n: int, d: int) -> list[dict]:
    """Precomputed vertices (extremal spectra) of the (n, d) moment polytope.

    Each record carries the exact spectrum, its integer form, and denominator.
    """
    p = data_root() / "vertices" / f"vertices_{n}_{d}.json"
    if not p.exists():
        raise KeyError(f"no precomputed vertices for ({n},{d})")
    return json.loads(p.read_text())


def classification(n: int, d: int) -> list[dict]:
    """Precomputed verdict per vertex: index, integer form, verdict label."""
    p = data_root() / "census" / f"census_{n}_{d}_results.txt"
    if not p.exists():
        raise KeyError(f"no precomputed classification for ({n},{d})")
    out = []
    for line in p.read_text().splitlines():
        if not line[:5].strip().isdigit():
            continue
        idx_str, rest = line.strip().split(" ", 1)
        lb = rest.index("[")
        rb = rest.index("]")
        vec = [int(x) for x in rest[lb + 1:rb].replace(",", " ").split()]
        verdict = rest[rb + 1:].strip()
        out.append({"index": int(idx_str), "integer_form": vec, "verdict": verdict})
    return out


def states(n: int | None = None, d: int | None = None,
           index: int | None = None) -> list[dict]:
    """Precomputed closed-form states, optionally filtered by system and index.

    Each record carries the classification, integer form, denominator, and,
    when the engine certified one, the closed form (natural denominator,
    integer weights, pretty symbolic amplitudes, support determinants).
    """
    p = data_root() / "states.jsonl"
    if not p.exists():
        return []
    want = f"({n},{d})" if n is not None and d is not None else None
    out = []
    for line in p.read_text().splitlines():
        if not line.strip():
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue  # tolerate a partial trailing line while the file is being rewritten
        if want is not None and rec.get("system") != want:
            continue
        if index is not None and rec.get("index") != index:
            continue
        out.append(rec)
    return out


def resolve_states(n: int, d: int, index: int | None = None,
                   mode: str = "hybrid", **knobs) -> list[dict]:
    """Resolve closed-form states in one of three modes (the mode, not the
    record, is the provenance: a record in the lookup is precomputed by
    definition):

    - ``mode="precompute"``: serve only the shipped lookup; vertices absent
      from it come back as ``{"index": i, "available": False}``.
    - ``mode="solve"``: ignore the lookup and recompute every vertex with the
      engine (independent reproduction, or to push bounds via ``knobs``).
    - ``mode="hybrid"`` (default): serve the lookup where present, solve the
      rest.

    ``knobs`` (max_card, max_blocks, max_clique, max_cliques) are passed to the
    engine for the solved vertices. solve and hybrid can be slow: they solve.
    """
    from fractions import Fraction
    if mode not in ("precompute", "solve", "hybrid"):
        raise ValueError(f"unknown mode {mode!r}")
    pre = {r["index"]: r for r in states(n, d)}
    verd = {r["index"]: r["verdict"] for r in classification(n, d)}
    out = []
    for i, v in enumerate(vertices(n, d)):
        if index is not None and i != index:
            continue
        if mode != "solve" and i in pre:
            out.append(pre[i])
            continue
        if mode == "precompute":
            out.append({"system": f"({n},{d})", "index": i, "available": False})
            continue
        from gpc_census.states import certify_state
        spec = [Fraction(s) for s in v["spectrum"]]
        rec = certify_state(n, d, spec, verdict=verd.get(i), **knobs) or {}
        rec.update({"system": f"({n},{d})", "index": i,
                    "integer_form": v["integer_form"], "denominator": v["denominator"]})
        out.append(rec)
    return out


def _design_ok(d, integer_form, cf):
    from itertools import combinations
    support = [tuple(s) for s in cf["support_dets"]]
    for a, b in combinations((set(s) for s in support), 2):
        if len(a ^ b) == 2:                       # one-hop connected: not a design
            return False, "support not one-hop free"
    inc = [0] * d
    for w, s in zip(cf["weights"], support):
        for m in s:
            inc[m] += w
    if sorted(inc, reverse=True) != sorted(integer_form, reverse=True):
        return False, "incidence sums != spectrum"
    return True, ""


def _exact_ok(d, spectrum, cf):
    import numpy as np
    import sympy as sp
    support = [tuple(s) for s in cf["support_dets"]]
    amps = [complex(sp.sympify(p).evalf(30)) for p in cf["pretty"]]
    amap = dict(zip(support, amps))
    rho = np.zeros((d, d), complex)
    for t, ct in amap.items():
        for mp in t:
            s1 = (-1) ** t.index(mp)
            t2 = tuple(x for x in t if x != mp)
            for m in range(d):
                if m in t2:
                    continue
                tp = tuple(sorted(t2 + (m,)))
                if tp not in amap:
                    continue
                rho[m, mp] += s1 * (-1) ** tp.index(m) * np.conjugate(amap[tp]) * ct
    eigs = sorted(np.linalg.eigvalsh(rho).real, reverse=True)
    want = sorted((float(x) for x in spectrum), reverse=True)
    err = max(abs(a - b) for a, b in zip(eigs, want))
    return (err < 1e-9), f"spectrum residual {err:.2e}"


def validate_states(records) -> list[tuple]:
    """Re-check every certified record's closed form, independently of the solver
    that produced it. Returns a list of ``(system, index, tierB, reason)`` for the
    records that fail; an empty list means the atlas is self-consistent.

    EXACT-CONSTR (integer designs) are proved combinatorially (one-hop-free
    support and matching incidence sums); EXACT (interference / real designs) are
    re-checked numerically (evaluate the closed form, rebuild the 1-RDM, require
    its eigenvalues to match the spectrum to 1e-9). The exact symbolic proof
    already ran at construction; this catches corruption before shipping.
    """
    from fractions import Fraction
    fails = []
    for r in records:
        tb = r.get("tierB")
        cf = r.get("closed_form")
        if tb not in ("EXACT", "EXACT-CONSTR") or not cf:
            continue
        n, d = (int(x) for x in r["system"].strip("()").split(","))
        if tb == "EXACT-CONSTR":
            ok, why = _design_ok(d, r["integer_form"], cf)
        else:
            spectrum = [Fraction(x, r["denominator"]) for x in r["integer_form"]]
            ok, why = _exact_ok(d, spectrum, cf)
        if not ok:
            fails.append((r["system"], r["index"], tb, why))
    return fails


def export(n: int, d: int) -> dict:
    """Everything precomputed for one system, as one machine-readable object:
    constraints, vertices, classification, and closed-form states."""
    from gpc_census.constraints import constraints as _constraints
    try:
        con = _constraints(n, d)
    except KeyError:
        con = None
    return {
        "system": {"n": n, "d": d},
        "constraints": con,
        "vertices": vertices(n, d),
        "classification": classification(n, d),
        "states": states(n, d),
    }
