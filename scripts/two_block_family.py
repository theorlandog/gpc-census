#!/usr/bin/env python3
"""The two-block interference family, and what governs the two-fiber gap.

Two parts, both certified into results/data/two_block_family.json.

PART 1, the family. Exactly four census vertices are INTERFERENCE with a
two-valued spectrum. Two blocks is the minimal nontrivial degeneracy pattern,
both cancellation-regime states in the entire census live here, and the four
span every fiber behaviour currently known. This assembles one table: channel
census split into within-block and cross-block, both tangent dimensions,
holonomy fields, denominator arithmetic, gauge-minimality status, and cascade
behaviour.

PART 2, the scored predictions of docs/prereg_two_block_inversion.md. The
silence-rank lemma says first-order fixed-SPECTRUM weight deformations see
only WITHIN-block channel conditions while fixed-RHO deformations see all of
them, so cross-block channels are the entire first-order difference between
the two fibers. The predictions ask whether that mechanism is also a
classification: does the cross-block count X predict the gap

    GAP = fixed_spectrum_tangent_dim - fixed_rho_tangent_dim

state by state, across all 799 certified states?

Independent scope is reported over TRANSPORT CLASSES, not row counts, per
standing rule R2 of docs/prereg_template.md: padding and frozen-core lifts
carry a state's fiber data unchanged, so sibling rows are not independent
observations.

Usage:
  python scripts/two_block_family.py            # write the artifact
  python scripts/two_block_family.py --check    # fail if stale
  python scripts/two_block_family.py --family   # the four-vertex table only
  python scripts/two_block_family.py --score    # the scorecard only
"""
from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.states import one_hop_classes  # noqa: E402

DATA = ROOT / "results" / "data"
OUT = DATA / "two_block_family.json"


def _load(name):
    return [json.loads(ln) for ln in (DATA / name).read_text().splitlines() if ln.strip()]


def _keyed(name):
    return {(r["system"], r["index"]): r for r in _load(name)}


def particle_number(system: str) -> int:
    return int(system.strip("()").split(",")[0])


def transport_class(system: str, integer_form, denominator):
    """Canonical key modulo the two verdict-preserving functorial maps.

    Trailing-zero padding appends zeros; a frozen-core lift prepends an
    occupation 1 and raises N by one. Two states with the same key carry the
    same fiber data by construction, so they are one observation, not several.
    Particle-hole duality is NOT quotiented here (the census ships only one
    member of each dual pair among these six systems, except that (5,10) is
    self-dual); this is noted in the artifact rather than silently assumed.
    """
    n = particle_number(system)
    v = [Fraction(x, denominator) for x in integer_form]
    while v and v[-1] == 0:
        v.pop()
    while v and v[0] == 1:
        v.pop(0)
        n -= 1
    return f"N={n}:" + ",".join(str(x) for x in v)


def channel_split(record):
    """One-hop classes of the certified support, split by block relation.

    A class of mode pair (A, B) is WITHIN-BLOCK when lambda_A == lambda_B and
    CROSS-BLOCK otherwise. Cross-block classes are exactly the conditions the
    silence-rank lemma says fixed-spectrum deformations do not see.
    """
    lam = record["integer_form"]
    dets = [tuple(t) for t in record["closed_form"]["support_dets"]]
    within, cross = [], []
    for (a, b), pairs in sorted(one_hop_classes(dets).items()):
        entry = {"modes": [a, b], "class_size": len(pairs), "pairs": pairs}
        (within if lam[a] == lam[b] else cross).append(entry)
    return within, cross


def rows():
    """Per-state record joining spectrum, channels, and both fiber tangents."""
    states = _keyed("states.jsonl")
    cascade = _keyed("cascade.jsonl")
    out = []
    for key, rec in sorted(states.items()):
        casc = cascade[key]
        lam = rec["integer_form"]
        within, cross = channel_split(rec)
        fs = casc["fixed_spectrum_tangent_dim"]
        fr = casc["fixed_rho_tangent_dim"]
        out.append(
            {
                "system": key[0],
                "index": key[1],
                "label": f"{key[0]} v{key[1]}",
                "classified": rec["classified"],
                "blocks": len(set(lam)),
                "block_sizes": sorted(
                    (lam.count(v) for v in sorted(set(lam), reverse=True)), reverse=True
                ),
                "transport_class": transport_class(key[0], lam, rec["denominator"]),
                "support_certified": casc["support_certified"],
                "incidence_kernel_dim": casc["incidence_kernel_dim"],
                "channels_total": len(within) + len(cross),
                "channels_within_block": len(within),
                "channels_cross_block": len(cross),
                "X": len(cross),
                "fixed_spectrum_tangent_dim": fs,
                "fixed_rho_tangent_dim": fr,
                "gap": fs - fr,
                "cancellation_regime": casc["cancellation_regime"],
                "evidence": "numerical",
            }
        )
    return out


# --------------------------------------------------------------------------
# part 1: the four-vertex family


def family_table(all_rows):
    states = _keyed("states.jsonl")
    cascade = _keyed("cascade.jsonl")
    gauge = _keyed("gauge_min.jsonl")
    holo = json.loads((DATA / "holonomy_fields.json").read_text())
    by_state = {}
    for loop in holo["loops"]:
        by_state.setdefault((loop["system"], loop["index"]), []).append(loop)

    members = [
        r for r in all_rows if r["blocks"] == 2 and r["classified"] == "INTERFERENCE"
    ]
    out = []
    for r in members:
        key = (r["system"], r["index"])
        rec, casc, g = states[key], cascade[key], gauge[key]
        within, cross = channel_split(rec)
        loops = by_state.get(key, [])
        out.append(
            {
                **r,
                "spectrum": f"{rec['integer_form']} / {rec['denominator']}",
                "state_denominator": rec["closed_form"]["den"],
                "denominator_ratio": str(
                    Fraction(rec["closed_form"]["den"], rec["denominator"])
                ),
                "within_block_channels": within,
                "cross_block_channels": cross,
                "loops": len(loops),
                "holonomy_min_polys": [loop["min_poly"] for loop in loops],
                "holonomy_galois_groups": [loop["galois_group"] for loop in loops],
                "support_upper": casc["support_upper"],
                "cascade_drop": casc["support_certified"] - casc["support_upper"],
                "cascade_blocked_reason": casc["blocked_reason"],
                "gauge_lower_bound": g["gauge_lower_bound"],
                "gauge_support_no_min": g["support_no_min"],
                "gauge_minimal_certified": g["minimal_certified"],
                "profile_class_count": g["profile_class_count"],
            }
        )
    return sorted(out, key=lambda r: (-r["support_certified"], r["label"]))


# --------------------------------------------------------------------------
# part 2: the scored predictions


def _classes(sel):
    return sorted({r["transport_class"] for r in sel})


def score(all_rows):
    n = len(all_rows)
    gap_pos = [r for r in all_rows if r["gap"] > 0]
    x_pos = [r for r in all_rows if r["X"] > 0]
    inversions = [
        r
        for r in all_rows
        if r["fixed_rho_tangent_dim"] == 0 and r["fixed_spectrum_tangent_dim"] > 0
    ]
    fs_pos = [r for r in all_rows if r["fixed_spectrum_tangent_dim"] > 0]

    def verdict(ok, exceptions, claim, expectation, scope_rows, note=""):
        return {
            "claim": claim,
            "registered_expectation": expectation,
            "verdict": "PASS" if ok else "FAIL",
            "scope": len(scope_rows),
            "independent_scope_transport_classes": len(_classes(scope_rows)),
            "exceptions": len(exceptions),
            "exception_labels": [r["label"] for r in exceptions[:25]],
            "exception_transport_classes": len(_classes(exceptions)),
            "note": note,
        }

    p1_ex = [r for r in all_rows if r["gap"] != r["X"]]
    p2_ex = [r for r in all_rows if r["gap"] != 2 * r["X"]]
    p3_ex = [r for r in all_rows if (r["gap"] > 0) != (r["X"] > 0)]
    p4_ex = [r for r in inversions if r["X"] == 0]

    # the two halves of P-INV-3, scored separately because they came apart
    forward_ex = [r for r in all_rows if r["gap"] > 0 and r["X"] == 0]
    reverse_ex = [r for r in all_rows if r["X"] > 0 and r["gap"] == 0]

    two_block = [r for r in all_rows if r["blocks"] == 2]
    canc = [r for r in all_rows if r["cancellation_regime"]]
    p5_ex = [r for r in canc if r["blocks"] != 2]

    fam = [r for r in two_block if r["classified"] == "INTERFERENCE"]
    p6_ex = []
    for r in fam:
        rigid_both = r["fixed_spectrum_tangent_dim"] == 0 and r["fixed_rho_tangent_dim"] == 0
        predicted_rigid = r["incidence_kernel_dim"] == 0 or r["X"] == 0
        if rigid_both != predicted_rigid:
            p6_ex.append(r)

    return {
        "prereg": "docs/prereg_two_block_inversion.md",
        "scope_note": (
            "Scored on all 799 certified states. Independent scope is counted "
            "in TRANSPORT CLASSES (padding and frozen-core siblings collapse), "
            "per standing rule R2 of docs/prereg_template.md."
        ),
        "totals": {
            "states": n,
            "transport_classes": len(_classes(all_rows)),
            "states_with_gap": len(gap_pos),
            "transport_classes_with_gap": len(_classes(gap_pos)),
            "states_with_cross_block_channel": len(x_pos),
            "transport_classes_with_cross_block_channel": len(_classes(x_pos)),
            "inversion_states": len(inversions),
            "inversion_transport_classes": len(_classes(inversions)),
            "gap_ever_negative": sum(1 for r in all_rows if r["gap"] < 0),
        },
        "P_INV_1": verdict(
            not p1_ex, p1_ex, "GAP == X on every state", "FAIL", all_rows,
            "Failed as expected. A channel condition is complex, hence two real "
            "conditions, and the two fibers have different gauge dimensions, so "
            "a bare class count cannot equal a dimension.",
        ),
        "P_INV_2": verdict(
            not p2_ex, p2_ex, "GAP == 2 * X on every state", "uncertain", all_rows,
            "Failed. The factor-two bookkeeping would hold only if every "
            "cross-block condition were independent on the incidence kernel; "
            "with several channels on overlapping supports they are not.",
        ),
        "P_INV_3": verdict(
            not p3_ex, p3_ex, "GAP > 0 if and only if X > 0", "PASS", all_rows,
            "FAILED AS A BICONDITIONAL, but it failed in one direction only. "
            "See P_INV_3_forward and P_INV_3_reverse: the necessary direction "
            "holds with zero exceptions and is the load-bearing half.",
        ),
        "P_INV_3_forward": verdict(
            not forward_ex, forward_ex,
            "GAP > 0 implies X > 0 (a cross-block channel is NECESSARY for the "
            "two fibers to disagree)",
            "PASS (the half that survives)", all_rows,
            "Zero exceptions on all 799 states. This is the silence-rank "
            "lemma's cross-block clause promoted from a dimension formula "
            "verified at two states to a census-wide necessary condition.",
        ),
        "P_INV_3_reverse": verdict(
            not reverse_ex, reverse_ex,
            "X > 0 implies GAP > 0 (a cross-block channel is SUFFICIENT)",
            "PASS (registered inside P-INV-3)", all_rows,
            "Failed. Every exception has incidence kernel dimension 0 and is "
            "rigid in BOTH fibers: the cross-block condition exists but there "
            "is no first-order weight freedom for it to cut, so the gap is "
            "zero for want of anything to differ about.",
        ),
        "P_INV_3_restricted": verdict(
            not [r for r in fs_pos if (r["gap"] > 0) != (r["X"] > 0)],
            [r for r in fs_pos if (r["gap"] > 0) != (r["X"] > 0)],
            "restricted to states that deform at all in the fixed-spectrum "
            "fiber, GAP > 0 if and only if X > 0",
            "not registered in advance; POST-HOC",
            fs_pos,
            "POST-HOC, derived after seeing which direction failed. Reported "
            "because it is the exact biconditional the data supports, but it "
            "has had no out-of-sample test and must not be quoted as one.",
        ),
        "P_INV_4": verdict(
            not p4_ex, p4_ex, "every inversion state has X > 0", "PASS", inversions,
            "Holds on all 14 inversion states, 7 transport classes.",
        ),
        "P_INV_5": {
            "claim": "the cancellation regime occurs only at two-block spectra",
            "registered_expectation": "PASS in sample, UNINFORMATIVE as evidence",
            "verdict": "UNINFORMATIVE",
            "scope": len(canc),
            "independent_scope_transport_classes": len(_classes(canc)),
            "exceptions": len(p5_ex),
            "exception_labels": [r["label"] for r in p5_ex],
            "note": (
                "Both known cancellation-regime members are two-block and they "
                "are the only two known, so the prediction could not fail on "
                "this corpus. Scored UNINFORMATIVE per standing rule R3, not "
                "PASS. Stage 1 is what would test it."
            ),
        },
        "P_INV_6": verdict(
            not p6_ex, p6_ex,
            "within the four-vertex family, rigid in both fibers if and only if "
            "dim K == 0 or X == 0",
            "uncertain", fam,
            "Holds on all four: v89 has one channel and it is WITHIN-block "
            "(X = 0), v32 has incidence kernel 0, and both are rigid in both "
            "fibers; v103 and v108 have cross-block channels and positive "
            "kernels and neither is rigid. Four states, four transport "
            "classes, so this is a rule fitted to its own sample and is "
            "reported as a description of the family, not as a law.",
        ),
    }


def sector_split(system, index):
    """Decompose a fixed-spectrum tangent into weight and phase directions.

    The ambient real coordinates are ordered so that 0..m-1 perturb the real
    parts of the amplitudes (the WEIGHT sector) and m..2m-1 the imaginary
    parts (the PHASE sector); the gauge directions are pure phase. At a REAL
    state the silence-rank lemma's proof splits the two sectors orthogonally,
    so the tangent should decompose exactly; at a genuinely complex state it
    need not, and that is worth recording rather than assuming.
    """
    import numpy as np

    import gpc_census.fiber as fib

    rec = _keyed("states.jsonl")[(system, index)]
    c, dets = fib.amplitudes_from_record(rec)
    c, dets, d = fib._normalize(c, dets, None)
    m = len(dets)
    rho = fib._sesquilinear_rdm(c, c, dets, d)
    _, vecs, blocks = fib.spectral_blocks(rho)
    rotated = [vecs.conj().T @ D @ vecs for D in fib._differentials(c, dets, d)]
    pairs = [(p, q) for b in blocks for p in b for q in b if p <= q]
    constraint = np.hstack(
        [
            fib._herm_rows(rotated, pairs),
            np.array([[2 * np.real(np.vdot(c, dc)) for dc in fib._unit_dirs(m)]]).T,
        ]
    ).T
    null = fib._nullspace(constraint)
    gauge = fib.gauge_generators(c, dets, d)
    gauge_dim = int(np.linalg.matrix_rank(gauge, tol=1e-8)) if gauge.size else 0
    moved = int(np.linalg.matrix_rank(constraint @ gauge.T, tol=1e-8)) if gauge.size else 0
    gauge_in_kernel = max(0, gauge_dim - moved)

    weight_only = int(null.shape[0] - np.linalg.matrix_rank(null[:, m:], tol=1e-8))
    phase_only = int(null.shape[0] - np.linalg.matrix_rank(null[:, :m], tol=1e-8))
    decouples = weight_only + phase_only == int(null.shape[0])
    return {
        "label": f"{system} v{index}",
        "kernel_dim": int(null.shape[0]),
        "gauge_in_kernel": gauge_in_kernel,
        "tangent_dim": int(null.shape[0]) - gauge_in_kernel,
        "weight_sector_directions": weight_only,
        "phase_sector_directions": phase_only,
        # Only meaningful once the sectors decouple: the gauge directions are
        # pure phase, so subtracting them from the phase count assumes the
        # decomposition is a direct sum. At a complex state it is not, and the
        # subtraction would report a negative "dimension".
        "phase_sector_beyond_gauge": (phase_only - gauge_in_kernel) if decouples else None,
        "sectors_decouple": decouples,
        "note": (
            "Sectors decouple exactly at a REAL state, as the silence-rank "
            "lemma's proof requires: the weight part of a deformation is the "
            "real sector and the phase part the imaginary one. Where they do "
            "not decouple the state is genuinely complex, the lemma's "
            "real-state bookkeeping does not apply verbatim, and the "
            "weight/phase split is not defined."
            if not decouples
            else "Sectors decouple, as the lemma's real-state proof requires."
        ),
    }


def stratification(all_rows):
    """Is two-block the right family, or is the cross-block count?"""
    import collections

    def dist(sel):
        return dict(sorted(collections.Counter(r["blocks"] for r in sel).items()))

    gap = [r for r in all_rows if r["gap"] > 0]
    xpos = [r for r in all_rows if r["X"] > 0]
    inv = [
        r
        for r in all_rows
        if r["fixed_rho_tangent_dim"] == 0 and r["fixed_spectrum_tangent_dim"] > 0
    ]
    two = [r for r in all_rows if r["blocks"] == 2]
    two_design = [r for r in two if r["classified"] != "INTERFERENCE"]
    return {
        "question": (
            "Task 3 of the handoff: does the two-block family have a natural "
            "closure that grows the sample past four?"
        ),
        "block_count_distribution": {
            "all_states": dist(all_rows),
            "cross_block_channel_present": dist(xpos),
            "two_fiber_gap": dist(gap),
            "inversions": dist(inv),
        },
        "two_block_vertices": len(two),
        "two_block_by_class": dict(
            sorted(collections.Counter(r["classified"] for r in two).items())
        ),
        "two_block_designs_with_any_channel": sum(
            1 for r in two_design if r["channels_total"] > 0
        ),
        "two_block_with_cross_block_channel": sum(1 for r in two if r["X"] > 0),
        "rate_cross_block_two_block": f"{sum(1 for r in two if r['X'] > 0)}/{len(two)}",
        "rate_cross_block_all": f"{len(xpos)}/{len(all_rows)}",
        "answer": (
            "No. Two-block is not the operative invariant. The two-fiber gap "
            "population spans block counts 2 through 6 (2:2, 3:21, 4:31, 5:18, "
            "6:1) and two-block vertices are UNDER-represented in it: 3 of 64 "
            "two-block vertices carry a cross-block channel against 155 of 799 "
            "overall. The design end of the two-block family is trivial for a "
            "structural reason rather than an interesting one: all 60 two-block "
            "DESIGN vertices have zero one-hop channels, because design "
            "supports are one-hop free by definition, so X = 0 and no gap is "
            "possible. The four two-block interference vertices span the known "
            "fiber behaviours because they happen to span (X, dim K), not "
            "because two blocks causes anything. The natural closure of the "
            "family is the X-stratification, which already has 155 members in "
            "72 transport classes."
        ),
        "surviving_two_block_fact": (
            "Both known cancellation-regime states are two-block. That sample "
            "is size two and P_INV_5 is scored UNINFORMATIVE; it is the one "
            "reason two-block spectra remain worth targeting in Stage 1, and "
            "only for the cancellation-regime question."
        ),
    }


def build():
    all_rows = rows()
    return {
        "generated_by": "scripts/two_block_family.py",
        "what": (
            "The four two-block INTERFERENCE vertices in full, and the "
            "census-wide test of whether cross-block channel structure governs "
            "the fixed-rho / fixed-spectrum gap."
        ),
        "evidence": (
            "Channel splits and block structure are EXACT (combinatorial). "
            "Tangent dimensions are NUMERICAL, inherited from "
            "results/data/cascade.jsonl."
        ),
        "family": family_table(all_rows),
        "tangent_sector_split": [
            sector_split("(3,10)", 103),
            sector_split("(3,10)", 89),
            sector_split("(3,10)", 108),
            sector_split("(3,8)", 32),
        ],
        "stratification": stratification(all_rows),
        "scorecard": score(all_rows),
        "states": all_rows,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--family", action="store_true")
    ap.add_argument("--score", action="store_true")
    a = ap.parse_args()

    doc = build()
    if a.family:
        print(json.dumps(doc["family"], indent=2))
        return 0
    if a.score:
        print(json.dumps(doc["scorecard"], indent=2))
        return 0

    text = json.dumps(doc, indent=1) + "\n"
    if a.check:
        cur = Path(a.out).read_text() if Path(a.out).exists() else ""
        if cur != text:
            print(f"STALE: {a.out}", file=sys.stderr)
            return 1
        print(f"{a.out} is current")
        return 0

    Path(a.out).write_text(text)
    sc = doc["scorecard"]
    print(f"wrote {a.out}")
    print(f"  family members: {len(doc['family'])}")
    t = sc["totals"]
    print(f"  states {t['states']} in {t['transport_classes']} transport classes")
    print(f"  gap>0: {t['states_with_gap']} states / "
          f"{t['transport_classes_with_gap']} classes")
    print(f"  inversions: {t['inversion_states']} states / "
          f"{t['inversion_transport_classes']} classes")
    for sp in doc["tangent_sector_split"]:
        if sp["sectors_decouple"]:
            print(f"  {sp['label']:12} tangent {sp['tangent_dim']} = "
                  f"{sp['weight_sector_directions']} weight + "
                  f"{sp['phase_sector_beyond_gauge']} phase")
        else:
            print(f"  {sp['label']:12} tangent {sp['tangent_dim']}, sectors do "
                  f"NOT decouple (complex state); no weight/phase split")
    for key in ("P_INV_1", "P_INV_2", "P_INV_3", "P_INV_3_forward",
                "P_INV_3_reverse", "P_INV_3_restricted", "P_INV_4", "P_INV_5",
                "P_INV_6"):
        v = sc[key]
        print(f"  {key:20} {v['verdict']:14} exceptions {v['exceptions']:>3}"
              f"   (registered: {v['registered_expectation']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
