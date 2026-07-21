"""Settle the rank-11 Stage-0 bracket for N=3 without generating the rank-11
constraint system, using four exact methods (all witnesses checked here):

  STATE-TRANSPORT  a certified census vertex with the same core (padding /
                   frozen-core lift) transports its state; verify_exact certifies.
  EXPLICIT-STATE   two hand constructions: the uniform (3/11)^11 as a Z11
                   difference design (base block {0,1,3}), and the frozen-core
                   lift of an N=2 paired state for (1, (1/5)^10).
  N2-PAIRING       lambda_1 = 1 pins a frozen core, so attainability reduces to
                   the residual N=2 spectrum; a 2-fermion 1-RDM has every nonzero
                   eigenvalue at even multiplicity (antisymmetric-matrix singular
                   values pair up), so an odd multiplicity refutes.
  ZERO-RESTRICTION trailing zeros force the state into wedge^3 H_{d'}; the
                   restricted point must satisfy the known (3, d') GPCs, else
                   refuted.

Writes docs/bracket_3_11_settlement.json. Reproduces 19 TRUE / 27 REFUTED / 4 OPEN.
"""
from __future__ import annotations

import json
import pathlib
import sys
from collections import Counter
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
sys.path.insert(0, str(ROOT / "scripts"))
from transport_states import core_of, core_state  # noqa: E402


def spectrum(intf, den):
    return [Fraction(x, den) for x in intf]


def try_transport(spec):
    """Certified census (N=3) vertex with the same core transports its state."""
    from gpc_census.exactify import verify_exact

    core = core_of([str(x) for x in spec])
    recs = [json.loads(ln) for ln in (DATA / "states.jsonl").read_text().splitlines() if ln.strip()]
    for r in recs:
        if not r["system"].startswith("(3,"):
            continue
        if r.get("tierB") not in {"EXACT", "EXACT-CONSTR"} or "closed_form" not in r:
            continue
        n, d = (int(x) for x in r["system"].strip("()").split(","))
        v = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())[r["index"]]
        if core_of(v["spectrum"]) != core:
            continue
        cdets, amps = core_state(r["closed_form"])
        ones = sum(1 for x in spec if x == 1)
        tdets = [list(range(ones)) + [o + ones for o in det] for det in cdets]
        if verify_exact(3, 11, spec, tdets, amps):
            return {"verdict": "TRUE-VERTEX", "certificate": "STATE-TRANSPORT",
                    "donor": {"system": r["system"], "index": r["index"]}}
    return None


def try_explicit(spec):
    import sympy as sp

    from gpc_census.exactify import verify_exact

    # uniform (3/11)^11 : Z11 difference design, base block {0,1,3}
    if spec == [Fraction(3, 11)] * 11:
        dets = [sorted(((0 + s) % 11, (1 + s) % 11, (3 + s) % 11)) for s in range(11)]
        amps = [sp.sqrt(sp.Rational(1, 11))] * 11
        if verify_exact(3, 11, spec, dets, amps):
            return {"verdict": "TRUE-VERTEX", "certificate": "EXPLICIT-STATE",
                    "note": "Z11 difference design, base block {0,1,3}; one-hop-free; DESIGN-INT",
                    "state": {"support_dets": dets, "weights": [1] * 11, "den": 11,
                              "pretty": ["sqrt(1/11)"] * 11}}
    # (1, (1/5)^10): frozen-core lift of an N=2 paired state
    if spec == [Fraction(1)] + [Fraction(1, 5)] * 10:
        dets = [[0, 1, 2], [0, 3, 4], [0, 5, 6], [0, 7, 8], [0, 9, 10]]
        amps = [sp.sqrt(sp.Rational(1, 5))] * 5
        if verify_exact(3, 11, spec, dets, amps):
            return {"verdict": "TRUE-VERTEX", "certificate": "EXPLICIT-STATE",
                    "note": "frozen-core lift of N=2 paired state",
                    "state": {"support_dets": dets, "weights": [1] * 5, "den": 5,
                              "pretty": ["sqrt(1/5)"] * 5}}
    return None


def try_n2_pairing(spec):
    if spec[0] != 1:
        return None
    resid = spec[1:]                       # residual N=2 spectrum after freezing the core
    mult = Counter(x for x in resid if x != 0)
    odd = sorted({x for x, m in mult.items() if m % 2 == 1})
    if odd:
        vals = ", ".join(str((x * x.denominator)) for x in odd)
        return {"verdict": "REFUTED", "certificate": "N2-PAIRING",
                "note": f"lambda_1=1 pins frozen core; residual N=2 spectrum has odd "
                        f"multiplicity of [{vals}] (violates the even-degeneracy theorem "
                        f"for 2-fermion 1-RDMs)"}
    return None


def try_zero_restriction(spec):
    from gpc_census import constraints as C

    zeros = sum(1 for x in spec if x == 0)
    if zeros == 0:
        return None
    dprime = 11 - zeros
    pt = [x for x in spec if x != 0]
    if len(pt) != dprime or (3, dprime) not in [(3, d) for d in range(6, 11)]:
        return None
    try:
        ineqs = C.constraints(3, dprime)["inequalities"]
    except Exception:
        return None
    for q in ineqs:
        if sum(Fraction(c) * x for c, x in zip(q["coeffs"], pt)) > Fraction(q["rhs"]):
            return {"verdict": "REFUTED", "certificate": "ZERO-RESTRICTION",
                    "note": f"trailing zeros force support into wedge^3 H_{dprime}; "
                            f"restricted point violates the known (3,{dprime}) system"}
    return None


def main() -> int:
    bracket = json.loads((ROOT / "docs" / "bracket_3_11.json").read_text())
    out = {"system": "(3,11)", "source_bracket": "docs/bracket_3_11.json",
           "method": "state transport + frozen-core/N=2 pairing + zero-restriction "
                     "membership + explicit constructions (exact arithmetic)",
           "candidates": []}
    tally: Counter = Counter()
    for i, cand in enumerate(bracket["outer_vertices"]):
        spec = spectrum(cand["integer_form"], cand["denominator"])
        rec = {"index": i, "integer_form": cand["integer_form"],
               "denominator": cand["denominator"]}
        settled = (try_transport(spec) or try_explicit(spec)
                   or try_n2_pairing(spec) or try_zero_restriction(spec))
        if settled is None:
            settled = {"verdict": "OPEN", "certificate": None,
                       "note": "genuinely rank-11 full-support candidate; needs Tier A "
                               "at d=11 or Stage 1"}
        rec.update(settled)
        tally[settled["verdict"]] += 1
        out["candidates"].append(rec)
    out["tally"] = {"TRUE-VERTEX": tally["TRUE-VERTEX"], "REFUTED": tally["REFUTED"],
                    "OPEN": tally["OPEN"]}
    (ROOT / "docs" / "bracket_3_11_settlement.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print(f"settled {tally['TRUE-VERTEX']} TRUE / {tally['REFUTED']} REFUTED / "
          f"{tally['OPEN']} OPEN of {len(out['candidates'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
