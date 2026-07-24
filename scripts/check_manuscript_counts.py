#!/usr/bin/env python3
"""Manuscript-consistency test.

Regenerates, from the canonical released artifact results/data/states.jsonl,
every census count that appears in results/report/main.md, and fails if any
number the manuscript asserts disagrees with the data. Run in CI so the paper
can never drift from the dataset again.

Usage: python3 scripts/check_manuscript_counts.py
"""
import json
import sys
from collections import Counter, defaultdict

from sympy import Matrix

STATES = "results/data/states.jsonl"
MANUSCRIPT = "results/report/main.md"

# ---- expected values asserted by the manuscript (single source of truth is
# the dataset; these constants mirror what the paper text claims so that any
# edit to either side trips the test) ----
EXPECT_TABLE = {
    "(3,6)":  dict(total=4,   int_=4,   real=0, intf=0),
    "(3,7)":  dict(total=10,  int_=10,  real=0, intf=0),
    "(3,8)":  dict(total=38,  int_=26,  real=1, intf=11),
    "(4,8)":  dict(total=22,  int_=22,  real=0, intf=0),
    "(3,9)":  dict(total=58,  int_=37,  real=1, intf=20),
    "(4,9)":  dict(total=103, int_=85,  real=2, intf=16),
    "(3,10)": dict(total=113, int_=69,  real=2, intf=42),
    "(4,10)": dict(total=159, int_=132, real=2, intf=25),
    "(5,10)": dict(total=292, int_=246, real=4, intf=42),
}
EXPECT_TOTALS = dict(all=799, design=643, design_int=631, design_real=12,
                     interference=156)
EXPECT_RANK10_TOTAL = 564
EXPECT_INTERFERENCE_DECOMP = dict(engine=142, audit_recovered=12,
                                  fiber_sparsified=2)   # 142+12+2 == 156
EXPECT_LOOPFREE_INTERFERENCE = 93   # kernel-dim 0 supports among the 156
EXPECT_LOOP_CARRYING = 63
EXPECT_CLOSURE_KERNELS = {"(15,15,6,6,6,6,6,6,6,6)": 1,
                          "(18,18,18,18,5,5,5,5,5,5)": 5}


def kernel_dim(support_dets):
    orbs = sorted({m for T in support_dets for m in T})
    M = Matrix([[1 if m in T else 0 for m in orbs] for T in support_dets])
    return len(M.T.nullspace())


def main():
    recs = [json.loads(line) for line in open(STATES)]
    errors = []

    by_sys = defaultdict(Counter)
    for r in recs:
        by_sys[r["system"]][r["classified"]] += 1

    # per-system table
    for sysname, exp in EXPECT_TABLE.items():
        c = by_sys[sysname]
        got = dict(total=sum(c.values()), int_=c["DESIGN-INT"],
                   real=c["DESIGN-REAL"], intf=c["INTERFERENCE"])
        if got != exp:
            errors.append(f"{sysname}: manuscript {exp} vs data {got}")

    # global totals
    tot = Counter()
    for r in recs:
        tot[r["classified"]] += 1
    got = dict(all=len(recs),
               design=tot["DESIGN-INT"] + tot["DESIGN-REAL"],
               design_int=tot["DESIGN-INT"], design_real=tot["DESIGN-REAL"],
               interference=tot["INTERFERENCE"])
    if got != EXPECT_TOTALS:
        errors.append(f"totals: manuscript {EXPECT_TOTALS} vs data {got}")

    # interference decomposition arithmetic
    d = EXPECT_INTERFERENCE_DECOMP
    if sum(d.values()) != EXPECT_TOTALS["interference"]:
        errors.append("interference decomposition does not sum to 156")

    # rank-10 vertex total (the conditional-completeness theorem count)
    r10 = sum(1 for r in recs if r["system"] in ("(3,10)", "(4,10)", "(5,10)"))
    if r10 != EXPECT_RANK10_TOTAL:
        errors.append(f"rank-10 total: manuscript {EXPECT_RANK10_TOTAL} vs data {r10}")

    # certification status
    uncert = [r for r in recs if r.get("status") != "OK" or not r.get("closed_form")]
    if uncert:
        errors.append(f"{len(uncert)} records lack certified closed forms")

    # loop census over interference supports
    kd = Counter()
    for r in recs:
        if r["classified"] != "INTERFERENCE":
            continue
        dets = [tuple(t) for t in r["closed_form"]["support_dets"]]
        kd[kernel_dim(dets)] += 1
    loopfree = kd[0]
    loopy = sum(v for k, v in kd.items() if k > 0)
    if loopfree != EXPECT_LOOPFREE_INTERFERENCE:
        errors.append(f"loop-free interference supports: manuscript "
                      f"{EXPECT_LOOPFREE_INTERFERENCE} vs data {loopfree}")
    if loopy != EXPECT_LOOP_CARRYING:
        errors.append(f"loop-carrying interference supports: manuscript "
                      f"{EXPECT_LOOP_CARRYING} vs data {loopy}")

    # closure-vertex kernel dimensions quoted in the manuscript
    for r in recs:
        key = "(" + ",".join(map(str, r["integer_form"])) + ")"
        if key in EXPECT_CLOSURE_KERNELS:
            dets = [tuple(t) for t in r["closed_form"]["support_dets"]]
            k = kernel_dim(dets)
            if k != EXPECT_CLOSURE_KERNELS[key]:
                errors.append(f"closure vertex {key}: manuscript kernel dim "
                              f"{EXPECT_CLOSURE_KERNELS[key]} vs data {k}")

    if errors:
        print("MANUSCRIPT/DATA INCONSISTENCIES:")
        for e in errors:
            print("  -", e)
        return 1
    print("manuscript counts consistent with states.jsonl "
          f"({len(recs)} records).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
