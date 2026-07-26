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

# ---- claims sourced from the companion data artifacts (galois tally,
# gauge minimization, sparsification cascade, v103 support-10 certificate);
# each block below re-reads the artifact the manuscript quotes ----
HOLONOMY = "results/data/holonomy_fields.json"
EXPECT_HOLONOMY = dict(loopy_states=63, loops=82,
                       degree_counts={"1": 37, "2": 32, "4": 13},
                       exp_degree_counts={"1": 8, "2": 29, "4": 32, "8": 13},
                       all_two_elementary=True, all_kernels_saturated=True,
                       loops_with_radicand_outside_weight_group=15)
GAUGE = "results/data/gauge_min_summary.json"
EXPECT_GAUGE = dict(states=799, certified_gauge_minimal=629,
                    states_that_gauge_sparsify=0)
CASCADE = "results/data/cascade_summary.json"
EXPECT_CASCADE = dict(states_with_support_drop=71,
                      total_amplitudes_removed=89,
                      multiplicity_drift_blocked=33,
                      drops_by_class={"INTERFERENCE": 71})
CASCADE_EXACT = "results/data/cascade_exact.json"
EXPECT_CASCADE_EXACT = dict(endpoints_attempted=71, certified=69)
ORBIT = "results/data/orbit_check.json"
EXPECT_ORBIT = dict(endpoints_tested=71, same_state_as_library=63,
                    new_states_found=8, endpoints_in_natural_orbital_basis=0)
V103_S10 = "results/data/v103_support10_certificate.json"
EXPECT_V103_S10 = dict(
    state_denominator=46852,
    squared_weights=["13/34", "5/34", "53/442", "195/1802", "5/68", "5/68",
                     "35/901", "1/34", "18/901", "84/11713"],
    support=10, certified_library_support=14)


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

    # galois tally quoted in the manuscript, from holonomy_fields.json
    h = json.load(open(HOLONOMY))["summary"]
    for k, want in EXPECT_HOLONOMY.items():
        if h.get(k) != want:
            errors.append(f"holonomy {k}: manuscript {want} vs data {h.get(k)}")

    # gauge-minimality numbers, from gauge_min_summary.json
    g = json.load(open(GAUGE))
    for k, want in EXPECT_GAUGE.items():
        if g.get(k) != want:
            errors.append(f"gauge {k}: manuscript {want} vs data {g.get(k)}")
    if g.get("states", 0) - g.get("certified_gauge_minimal", 0) != 170:
        errors.append("gauge: manuscript quotes 170 open, data disagrees")

    # sparsification-cascade numbers, from cascade_summary.json
    c = json.load(open(CASCADE))
    if c.get("states_with_support_drop") != EXPECT_CASCADE["states_with_support_drop"]:
        errors.append(f"cascade drops: manuscript "
                      f"{EXPECT_CASCADE['states_with_support_drop']} vs data "
                      f"{c.get('states_with_support_drop')}")
    if c.get("total_amplitudes_removed") != EXPECT_CASCADE["total_amplitudes_removed"]:
        errors.append(f"cascade amplitudes removed: manuscript "
                      f"{EXPECT_CASCADE['total_amplitudes_removed']} vs data "
                      f"{c.get('total_amplitudes_removed')}")
    if c.get("drops_by_class") != EXPECT_CASCADE["drops_by_class"]:
        errors.append(f"cascade drops by class: manuscript all-interference "
                      f"vs data {c.get('drops_by_class')}")
    if c.get("blocked_reasons", {}).get("multiplicity_drift") != \
            EXPECT_CASCADE["multiplicity_drift_blocked"]:
        errors.append(f"cascade multiplicity-drift count: manuscript "
                      f"{EXPECT_CASCADE['multiplicity_drift_blocked']} vs data "
                      f"{c.get('blocked_reasons', {}).get('multiplicity_drift')}")
    if not c.get("worst_residual_among_drops", 1) <= 2.3e-16:
        errors.append("cascade worst residual exceeds the manuscript bound")
    v103row = [m for m in c.get("cancellation_regime_members", [])
               if m.get("index") == 103]
    if not (v103row and v103row[0].get("support_certified") == 14
            and v103row[0].get("support_upper") == 10):
        errors.append("cascade v103 support 14 -> 10: data disagrees")

    # exact certification of the cascade bounds, from cascade_exact.json
    ce = json.load(open(CASCADE_EXACT))
    for k, want in EXPECT_CASCADE_EXACT.items():
        if ce.get(k) != want:
            errors.append(f"cascade_exact {k}: manuscript {want} vs data {ce.get(k)}")

    # 2-RDM orbit gate, from orbit_check.json
    o = json.load(open(ORBIT))
    for k, want in EXPECT_ORBIT.items():
        if o.get(k) != want:
            errors.append(f"orbit_check {k}: manuscript {want} vs data {o.get(k)}")

    # the v103 support-10 certificate the manuscript quotes verbatim
    v = json.load(open(V103_S10))
    if v.get("state_denominator") != EXPECT_V103_S10["state_denominator"]:
        errors.append(f"v103 support-10 denominator: manuscript "
                      f"{EXPECT_V103_S10['state_denominator']} vs data "
                      f"{v.get('state_denominator')}")
    if sorted(v.get("squared_weights", [])) != sorted(EXPECT_V103_S10["squared_weights"]):
        errors.append("v103 support-10 squared weights: manuscript list "
                      "disagrees with the certificate")
    if len(v.get("support_dets", [])) != EXPECT_V103_S10["support"]:
        errors.append("v103 support-10 support size: data disagrees")
    if v.get("certified_library_support") != EXPECT_V103_S10["certified_library_support"]:
        errors.append("v103 library support 14: data disagrees")
    if not v.get("all_checks_pass"):
        errors.append("v103 support-10 certificate reports a failing check")

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
