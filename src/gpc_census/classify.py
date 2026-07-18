"""Design/interference classification of vertex spectra.

Faithful port of scripts/census_engine.py: a weighted design is a
nonnegative weighting of determinants with prescribed mode sums whose
support is one-hop free (no two determinants sharing N-1 modes), the
condition that lets a phase-free superposition cancel all off-diagonal
one-body terms. Backend: pinned CP-SAT when importable at 9.15.6755,
otherwise CBC for both stages. Verdicts report their backend; feasible
verdicts are verified against an exact witness check including the
one-hop condition.
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
        from importlib.metadata import version

        return "cpsat" if version("ortools") == _PINNED else "cbc"
    except Exception:
        return "cbc"


BACKEND = _detect()


def _geometry(n: int, d: int):
    dets = list(combinations(range(d), n))
    rows = [[j for j, t in enumerate(dets) if m in t] for m in range(d)]
    conflicts = []
    for a in range(len(dets)):
        ta = set(dets[a])
        for b in range(a + 1, len(dets)):
            if len(ta & set(dets[b])) == n - 1:
                conflicts.append((a, b))
    return dets, rows, conflicts


def _verify_witness(weights, nv, den, rows, conflicts) -> bool:
    if any(w < 0 for w in weights) or sum(weights) != den:
        return False
    if any(sum(weights[j] for j in rows[m]) != nv[m] for m in range(len(nv))):
        return False
    return all(not (weights[a] and weights[b]) for a, b in conflicts)


def _int_stage(nv, den, dets, rows, conflicts):
    if BACKEND == "cpsat":
        from ortools.sat.python import cp_model

        m = cp_model.CpModel()
        k = [m.NewIntVar(0, den, f"k{t}") for t in range(len(dets))]
        y = [m.NewBoolVar(f"y{t}") for t in range(len(dets))]
        for t in range(len(dets)):
            m.Add(k[t] <= den * y[t])
            m.Add(k[t] >= y[t])
        for mo, nm in enumerate(nv):
            m.Add(sum(k[j] for j in rows[mo]) == nm)
        for a, b in conflicts:
            m.AddBoolOr([y[a].Not(), y[b].Not()])
        s = cp_model.CpSolver()
        s.parameters.max_time_in_seconds = 120
        s.parameters.num_workers = 2
        st = s.Solve(m)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            return "FEASIBLE", [int(s.Value(x)) for x in k]
        return ("INFEASIBLE", None) if st == cp_model.INFEASIBLE else ("UNKNOWN", None)
    import pulp

    prob = pulp.LpProblem("int", pulp.LpMinimize)
    k = [pulp.LpVariable(f"k{t}", 0, den, cat="Integer") for t in range(len(dets))]
    y = [pulp.LpVariable(f"y{t}", cat="Binary") for t in range(len(dets))]
    prob += 0
    for t in range(len(dets)):
        prob += k[t] <= den * y[t]
        prob += k[t] >= y[t]
    for mo, nm in enumerate(nv):
        prob += pulp.lpSum(k[j] for j in rows[mo]) == nm
    for a, b in conflicts:
        prob += y[a] + y[b] <= 1
    prob.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=300))
    status = pulp.LpStatus[prob.status]
    if status == "Optimal":
        return "FEASIBLE", [int(round(v.value())) for v in k]
    return ("INFEASIBLE", None) if status == "Infeasible" else (status, None)


def _real_stage(nv, dets, rows, conflicts, n):
    import pulp

    cap = float(sum(nv)) / n
    prob = pulp.LpProblem("real", pulp.LpMinimize)
    k = [pulp.LpVariable(f"k{t}", 0, cap) for t in range(len(dets))]
    y = [pulp.LpVariable(f"y{t}", cat="Binary") for t in range(len(dets))]
    prob += 0
    for t in range(len(dets)):
        prob += k[t] <= cap * y[t]
    for mo, nm in enumerate(nv):
        prob += pulp.lpSum(k[j] for j in rows[mo]) == nm
    for a, b in conflicts:
        prob += y[a] + y[b] <= 1
    prob.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=300))
    return pulp.LpStatus[prob.status]


def classify(n: int, d: int, spectrum) -> str:
    return classify_full(n, d, spectrum)["verdict"]


def classify_full(n: int, d: int, spectrum) -> dict:
    spectrum = [Fraction(x) for x in spectrum]
    den = 1
    for x in spectrum:
        den = den * x.denominator // math.gcd(den, x.denominator)
    nv = [int(x * den) for x in spectrum]
    dets, rows, conflicts = _geometry(n, d)

    di, witness = _int_stage(nv, den, dets, rows, conflicts)
    if di == "FEASIBLE":
        if not _verify_witness(witness, nv, den, rows, conflicts):
            return {"verdict": "UNRESOLVED(witness-failed)", "backend": BACKEND}
        return {"verdict": "DESIGN-INT", "backend": BACKEND, "witness": witness}
    dr = _real_stage(nv, dets, rows, conflicts, n)
    if dr == "Optimal":
        return {"verdict": "DESIGN-REAL", "backend": BACKEND}
    if dr == "Infeasible" and di == "INFEASIBLE":
        return {"verdict": "INTERFERENCE", "backend": BACKEND}
    return {"verdict": f"UNRESOLVED({di},{dr})", "backend": BACKEND}
