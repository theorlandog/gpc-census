"""Design/interference classification of vertex spectra.

Requires the classify extra: pip install gpc-census[classify]
Verdicts are solver certificates: DESIGN-INT means an integer weighted
design exists at the natural denominator, DESIGN-REAL means a real
nonnegative design exists but no integer one, INTERFERENCE means neither
exists and complex phase cancellation is forced.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations


def classify(n: int, d: int, spectrum: tuple[Fraction, ...]) -> str:
    try:
        import pulp
        from ortools.sat.python import cp_model
    except ImportError as e:
        raise RuntimeError("pip install gpc-census[classify] for solver support") from e
    den = 1
    for x in spectrum:
        den = den * x.denominator // __import__("math").gcd(den, x.denominator)
    nv = [int(x * den) for x in spectrum]
    dets = list(combinations(range(d), n))
    rows = [[j for j, t in enumerate(dets) if m in t] for m in range(d)]

    m = cp_model.CpModel()
    k = [m.NewIntVar(0, den, f"k{j}") for j in range(len(dets))]
    m.Add(sum(k) == den)
    for mo in range(d):
        m.Add(sum(k[j] for j in rows[mo]) == nv[mo])
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 120
    if solver.Solve(m) in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "DESIGN-INT"

    prob = pulp.LpProblem("real", pulp.LpMinimize)
    w = [pulp.LpVariable(f"w{j}", lowBound=0) for j in range(len(dets))]
    prob += 0
    prob += pulp.lpSum(w) == 1
    for mo in range(d):
        prob += pulp.lpSum(w[j] for j in rows[mo]) == float(spectrum[mo])
    prob.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=300))
    if pulp.LpStatus[prob.status] == "Optimal":
        return "DESIGN-REAL"
    if pulp.LpStatus[prob.status] == "Infeasible":
        return "INTERFERENCE"
    return f"UNRESOLVED({pulp.LpStatus[prob.status]})"
