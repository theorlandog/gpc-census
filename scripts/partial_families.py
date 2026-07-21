"""Stable Grassmann partial-constraint families for N=3 at arbitrary rank, and
an exact-arithmetic validator against the certified census data.

These are the ONLY constraint data in print for ranks 11 and 12 (the complete
systems there were computed once, in Schilling et al PRA 97, 052503 (2018),
and never published; see docs/RESEARCH.md). They are NECESSARY conditions
only: an outer bound, not the polytope. Two provenance tiers, tagged in the
output:

  PROVED   first-kind pairs lambda_{k+1}+lambda_{r-k} <= 1 (CMP 2008 Thm 4.3.1)
           and the four second-kind quadruples (CMP 2008 Thm 4.2.1); valid GPCs
           at every rank.
  WEAK     the odd-rank head pair lambda_1+lambda_r <= 1 + 2/(r-1) (thesis,
           stated without printed proof); a genuine bound, weaker than proved.
  CLAIMED  the level-5 series extending each second-kind quadruple by the tail
           starting at lambda_11 (CMP 2008 Remark 4.2.1, proof deferred;
           Klyachko arXiv:0904.2009, unrefereed letter). NOT a certificate:
           do not settle a verdict on these alone.

Regenerates docs/partial_families_3_11_3_12.json and runs the consistency and
restriction checks. Run: uv run scripts/partial_families.py
"""
from __future__ import annotations

import json
import pathlib
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"

# second-kind quadruples, rank independent (1-indexed), CMP 2008 Thm 4.2.1
QUADS = [[2, 3, 4, 5], [1, 3, 4, 6], [1, 2, 5, 6], [1, 2, 4, 7]]
# the level-5 series index set 1,2,4,7,11,16,... (gaps 1,2,3,4,5,...)
SERIES = [1, 2, 4, 7, 11, 16, 22]


def families(r: int):
    """Return [(tag, name, indices, rhs)] for N=3 at rank r (1-indexed)."""
    out = []
    # first-kind pairs (a, r+1-a) <= 1, a = 1 .. r//2 (CMP 2008 Thm 4.3.1). For
    # even r the adjacent boundary pair (r/2, r/2+1) is degenerate and omitted;
    # for odd r the head pair (1, r) is not proved and is replaced by the
    # weakened thesis bound l1 + l_r <= 1 + 2/(r-1).
    for a in range(1, r // 2 + 1):
        b = r + 1 - a
        if a >= b:
            continue
        if r % 2 == 0 and a == r // 2:
            continue  # (r/2, r/2+1) adjacent boundary pair, omitted
        if r % 2 == 1 and a == 1:
            out.append(("WEAK", f"1st-head l1+l{r}<=1+2/(r-1)", [1, r], 1 + F(2, r - 1)))
        else:
            out.append(("PROVED", f"1st-pair l{a}+l{b}<=1", [a, b], F(1)))
    for q in QUADS:
        out.append(("PROVED", f"2nd-quad {''.join(map(str,q))}<=2", list(q), F(2)))
    # level-5 series: each quadruple gets the series tail beyond its own indices
    tail = [i for i in SERIES if i <= r and i >= 8]  # indices 11,16,... within rank
    if tail:
        # AK-RMK-4.2.1 is the (1,2,4,7)+tail case; the other three are KLY09-EXT
        for j, q in enumerate([[1, 2, 4, 7], [2, 3, 4, 5], [1, 3, 4, 6], [1, 2, 5, 6]]):
            name = "AK-RMK-4.2.1" if j == 0 else f"KLY09-EXT-{j}"
            out.append(("CLAIMED", f"{name} {q}+{tail}<=2", q + tail, F(2)))
    return out


def _val(spec, idxs):
    return sum(spec[i - 1] for i in idxs)


def _viol(spec, fam, tags):
    return [name for tag, name, idxs, rhs in fam if tag in tags and _val(spec, idxs) > rhs]


def main() -> int:
    fam11 = families(11)
    fam12 = families(12)
    doc = {
        "note": "Stable Grassmann partial families for N=3; NECESSARY conditions "
                "(outer bound), not the complete system. Tiers: PROVED (CMP 2008 "
                "Thms 4.2.1/4.3.1), WEAK (thesis, unprinted proof), CLAIMED (CMP "
                "Remark 4.2.1 / Klyachko 0904.2009, not proved, NOT a certificate).",
        "sources": ["Altunbulak-Klyachko, Comm. Math. Phys. 282, 287 (2008)",
                    "Klyachko, arXiv:0904.2009 (2009)"],
        "systems": {},
    }
    for r, fam in ((11, fam11), (12, fam12)):
        doc["systems"][f"(3,{r})"] = [
            {"tag": tag, "name": name, "indices": idxs, "rhs": str(rhs)}
            for tag, name, idxs, rhs in fam]
    (ROOT / "docs" / "partial_families_3_11_3_12.json").write_text(
        json.dumps(doc, indent=1) + "\n")

    # ---- validation ----
    ok = True
    s = json.loads((ROOT / "docs" / "bracket_3_11_settlement.json").read_text())
    trues = [c for c in s["candidates"] if c["verdict"] == "TRUE-VERTEX"]
    opens = [c for c in s["candidates"] if c["verdict"] == "OPEN"]

    nbad = 0
    for c in trues:
        spec = [F(x, c["denominator"]) for x in c["integer_form"]]
        v = _viol(spec, fam11, {"PROVED", "WEAK", "CLAIMED"})
        if v:
            nbad += 1
            print(f"  (3,11) TRUE idx {c['index']} VIOLATES {v}")
    print(f"[consistency] {len(trues)} certified (3,11) TRUE vertices: "
          f"{'all satisfy every family' if nbad == 0 else str(nbad)+' violations'}")
    ok &= nbad == 0

    print("[lead] OPEN candidates vs PROVED / CLAIMED:")
    for c in opens:
        spec = [F(x, c["denominator"]) for x in c["integer_form"]]
        vp = _viol(spec, fam11, {"PROVED", "WEAK"})
        vc = _viol(spec, fam11, {"CLAIMED"})
        print(f"  idx {c['index']}: proved={vp or 'none'}  claimed={vc or 'none'}")

    v10 = json.loads((DATA / "vertices" / "vertices_3_10.json").read_text())
    nbad = 0
    for v in v10:
        spec = [F(x) for x in v["spectrum"]] + [F(0)]  # {lambda_11=0} face point
        if _viol(spec, fam11, {"PROVED", "WEAK"}):
            nbad += 1
    print(f"[restriction] {{lambda_11=0}} on {len(v10)} certified (3,10) vertices: "
          f"{'all satisfy proved+weak' if nbad == 0 else str(nbad)+' violate'}")
    ok &= nbad == 0

    tv = ROOT / "docs" / "bracket_3_12_true_vertices.json"
    if tv.exists():
        vv = json.loads(tv.read_text())
        rows = vv["true_vertices"] if isinstance(vv, dict) else vv
        nbad = 0
        for row in rows:
            spec = [F(x) for x in (row["spectrum"] if isinstance(row, dict) else row)]
            if _viol(spec, fam12, {"PROVED", "WEAK", "CLAIMED"}):
                nbad += 1
        print(f"[consistency] {len(rows)} face-embedded (3,12) TRUE vertices: "
              f"{'all satisfy every family' if nbad == 0 else str(nbad)+' violations'}")
        ok &= nbad == 0

    print("VALIDATION PASS" if ok else "VALIDATION FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
