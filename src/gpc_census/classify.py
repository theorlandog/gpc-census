"""Design/interference classification of vertex spectra.

Backend selection happens once at import: the pinned reference solver
(ortools CP-SAT 9.15.6755) is used when importable at exactly that
version; otherwise the integer stage falls back to CBC through pulp.
Every verdict reports which backend produced it. Feasible verdicts are
verified against an exact rational witness check, so they are
certificates regardless of backend; infeasible verdicts from the CBC
path carry floating-point tolerances (see certify module docstring).
"""
from __future__ import annotations

import importlib.util
import math
from fractions import Fraction
from itertools import combinations

_PINNED = "9.15.6755"


def _detect() -> str:
    if importlib.util.find_spec("ortools") is None:
        return "cbc"
    try:
        import ortools

        ver = getattr(ortools, "__version__", None)
        if ver is None:
            from importlib.metadata import version

            ver = version("ortools")
        return "cpsat" if ver == _PINNED else "cbc"
    except Exception:
        return "cbc"


BACKEND = _detect()


def _verify_witness(weights: list[int], nv: list[int], den: int, rows) -> bool:
    if sum(weights) != den or any(w < 0 for w in weights):
        return False
    return all(sum(weights[j] for j in rows[m]) == nv[m] for m in range(len(nv)))


def classify(n: int, d: int, spectrum) -> str:
    return classify_full(n, d, spectrum)["verdict"]


def classify_full(n: int, d: int, spectrum) -> dict:
    import pulp

    spectrum = [Fraction(x) for x in spectrum]
    den = 1
    for x in spectrum:
        den = den * x.denominator // math.gcd(den, x.denominator)
    nv = [int(x * den) for x in spectrum]
    dets = list(combinations(range(d), n))
    rows = [[j for j, t in enumerate(dets) if m in t] for m in range(d)]

    int_feasible = None
    witness = None
    if BACKEND == "cpsat":
        from ortools.sat.python import cp_model

        m = cp_model.CpModel()
        k = [m.NewIntVar(0, den, f"k{j}") for j in range(len(dets))]
        m.Add(sum(k) == den)
        for mo in range(d):
            m.Add(sum(k[j] for j in rows[mo]) == nv[mo])
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = 120
        st = solver.Solve(m)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            int_feasible = True
            witness = [int(solver.Value(x)) for x in k]
        elif st == cp_model.INFEASIBLE:
            int_feasible = False
    if int_feasible is None:
        prob = pulp.LpProblem("int", pulp.LpMinimize)
        w = [pulp.LpVariable(f"k{j}", lowBound=0, upBound=den, cat="Integer")
             for j in range(len(dets))]
        prob += 0
        prob += pulp.lpSum(w) == den
        for mo in range(d):
            prob += pulp.lpSum(w[j] for j in rows[mo]) == nv[mo]
        prob.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=300))
        status = pulp.LpStatus[prob.status]
        if status == "Optimal":
            int_feasible = True
            witness = [int(round(v.value())) for v in w]
        elif status == "Infeasible":
            int_feasible = False
        else:
            return {"verdict": f"UNRESOLVED(int,{status})", "backend": BACKEND}

    if int_feasible:
        if witness is None or not _verify_witness(witness, nv, den, rows):
            return {"verdict": "UNRESOLVED(witness-failed)", "backend": BACKEND}
        return {"verdict": "DESIGN-INT", "backend": BACKEND, "witness": witness}

    prob = pulp.LpProblem("real", pulp.LpMinimize)
    w = [pulp.LpVariable(f"w{j}", lowBound=0) for j in range(len(dets))]
    prob += 0
    prob += pulp.lpSum(w) == 1
    for mo in range(d):
        prob += pulp.lpSum(w[j] for j in rows[mo]) == float(spectrum[mo])
    prob.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=300))
    status = pulp.LpStatus[prob.status]
    if status == "Optimal":
        return {"verdict": "DESIGN-REAL", "backend": BACKEND}
    if status == "Infeasible":
        return {"verdict": "INTERFERENCE", "backend": BACKEND}
    return {"verdict": f"UNRESOLVED(real,{status})", "backend": BACKEND}
