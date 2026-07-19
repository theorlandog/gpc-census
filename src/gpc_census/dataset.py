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
