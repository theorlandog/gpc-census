"""Exact vertex enumeration of fermionic moment polytopes via lrslib."""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from fractions import Fraction
from pathlib import Path

from .constraints import constraints


def _hrep(n: int, d: int) -> str:
    sys_ = constraints(n, d)
    rows: list[str] = []
    for iq in sys_["inequalities"]:
        rows.append(" ".join([str(iq["rhs"])] + [str(-c) for c in iq["coeffs"]]))
    rows.append(" ".join(["1", "-1"] + ["0"] * (d - 1)))  # lambda_1 <= 1
    for i in range(d - 1):
        v = [0] * (d + 1)
        v[1 + i] = 1
        v[2 + i] = -1
        rows.append(" ".join(map(str, v)))
    v = [0] * (d + 1)
    v[d] = 1
    rows.append(" ".join(map(str, v)))  # lambda_d >= 0
    eqs = [" ".join([str(-iq["rhs"])] + [str(c) for c in iq["coeffs"]]) for iq in sys_["equalities"]]
    eqs.append(" ".join([f"-{n}"] + ["1"] * d))  # trace
    total = len(rows) + len(eqs)
    lin = " ".join(str(len(rows) + 1 + i) for i in range(len(eqs)))
    return (
        f"P\nH-representation\nlinearity {len(eqs)} {lin}\nbegin\n"
        f"{total} {d + 1} rational\n" + "\n".join(rows + eqs) + "\nend\n"
    )


def vertices(n: int, d: int) -> list[tuple[Fraction, ...]]:
    """Enumerate all vertices of the (n, d) polytope, exactly."""
    lrs = shutil.which("lrs")
    if lrs is None:
        raise RuntimeError("lrs binary not found; install lrslib (dnf install lrslib)")
    with tempfile.TemporaryDirectory() as tmp:
        ine = Path(tmp) / "p.ine"
        ine.write_text(_hrep(n, d))
        res = subprocess.run([lrs, str(ine)], capture_output=True, text=True, timeout=3600)
    verts: list[tuple[Fraction, ...]] = []
    buf: list[str] | None = None
    for line in res.stdout.splitlines():
        parts = line.split()
        if not parts or line.startswith("*"):
            continue
        if buf is None:
            if parts[0] == "1" and len(parts) >= 2:
                buf = parts[1:]
        else:
            buf += parts
        if buf is not None:
            if len(buf) == d:
                v = tuple(Fraction(x) for x in buf)
                if sum(v) == n:
                    verts.append(v)
                buf = None
            elif len(buf) > d:
                buf = None
    return sorted(set(verts), reverse=True)
