#!/usr/bin/env python3
"""The natural-orbital support ratio, and what it says about the trichotomy.

For each census record, compare the support of the sparse representative in
states.jsonl with the support of its natural-orbital representative in
states_natural_orbital.jsonl:

    no_support_ratio = s_NO / s_sparse

The motivating idea (docs/prereg_no_support_ratio.md) is that this one number
separates the three classes: DESIGN and CANCELLATION states carry sparsity and
diagonality in the SAME basis, so the ratio is 1, while generic INTERFERENCE
states must trade one for the other, so the ratio exceeds 1. The support
blowup would then be a MEASURE of interference rather than an artifact of how
the rotation was done.

DESIGN records already have a diagonal 1-RDM, so they do not appear in the
natural-orbital ledger and their ratio is 1 by definition rather than by
measurement. They are reported and excluded from the discriminating counts.

Usage:
  python scripts/no_support_ratio.py           # write the artifact
  python scripts/no_support_ratio.py --check
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.states import one_hop_classes  # noqa: E402

DATA = ROOT / "results" / "data"
OUT = DATA / "no_support_ratio.json"


def _lines(name):
    return [
        json.loads(ln) for ln in (DATA / name).read_text().splitlines() if ln.strip()
    ]


def transport_class(system, integer_form, denominator):
    """Same canonical key as scripts/two_block_family.py: modulo trailing-zero
    padding and frozen-core lifts."""
    n = int(system.strip("()").split(",")[0])
    v = [Fraction(x, denominator) for x in integer_form]
    while v and v[-1] == 0:
        v.pop()
    while v and v[0] == 1:
        v.pop(0)
        n -= 1
    return f"N={n}:" + ",".join(str(x) for x in v)


def build():
    sparse = {(r["system"], r["index"]): r for r in _lines("states.jsonl")}
    natural = {
        (r["system"], r["index"]): r for r in _lines("states_natural_orbital.jsonl")
    }
    cascade = {(r["system"], r["index"]): r for r in _lines("cascade.jsonl")}

    rows = []
    for key, rec in sorted(sparse.items()):
        casc = cascade[key]
        s_sparse = len(rec["closed_form"]["support_dets"])
        nat = natural.get(key)
        s_no = nat["support"] if nat else s_sparse
        dets = [tuple(t) for t in rec["closed_form"]["support_dets"]]
        classes = one_hop_classes(dets)
        cancellation = bool(casc["cancellation_regime"])
        rows.append(
            {
                "system": key[0],
                "index": key[1],
                "label": f"{key[0]} v{key[1]}",
                "classified": rec["classified"],
                "regime": (
                    "DESIGN"
                    if rec["classified"] != "INTERFERENCE"
                    else ("CANCELLATION" if cancellation else "INTERFERENCE-GENERIC")
                ),
                "transport_class": transport_class(
                    key[0], rec["integer_form"], rec["denominator"]
                ),
                "s_sparse": s_sparse,
                "s_natural_orbital": s_no,
                "no_support_ratio": s_no / s_sparse,
                "measured": bool(nat is not None),
                "one_hop_classes": len(classes),
                "max_class_size": max((len(p) for p in classes.values()), default=0),
                "incidence_kernel_dim": casc["incidence_kernel_dim"],
                "rdm_max_offdiagonal": casc["rdm_max_offdiagonal"],
                "loops": s_sparse - len({frozenset(t) for t in dets}) if False else None,
                "evidence": "NUMERICAL (natural-orbital support is a rotation result)",
            }
        )

    generic = [r for r in rows if r["regime"] == "INTERFERENCE-GENERIC"]
    cancel = [r for r in rows if r["regime"] == "CANCELLATION"]
    design = [r for r in rows if r["regime"] == "DESIGN"]

    flat_generic = [r for r in generic if r["no_support_ratio"] == 1]
    p1_exceptions = flat_generic + [r for r in cancel if r["no_support_ratio"] != 1]
    p2_exceptions = [r for r in rows if r["no_support_ratio"] < 1]

    def classes_of(sel):
        return len({r["transport_class"] for r in sel})

    # P-NSR-3: what distinguishes the flat generic records
    def profile(sel, field):
        return dict(sorted(Counter(r[field] for r in sel).items(), key=lambda kv: str(kv[0])))

    grew = [r for r in generic if r["no_support_ratio"] > 1]
    ratios = [r["no_support_ratio"] for r in generic]

    return {
        "generated_by": "scripts/no_support_ratio.py",
        "prereg": "docs/prereg_no_support_ratio.md",
        "evidence": "NUMERICAL",
        "disclosure": (
            "NOT a blind pre-registration. The aggregate counts (141 grew, 15 "
            "unchanged of 156) were already committed in "
            "results/data/natural_orbital_summary.json before the prediction "
            "was written, so P-NSR-1's failure was foreseeable. The identity "
            "and structure of the exceptions is the informative part."
        ),
        "populations": {
            "design": {"records": len(design), "transport_classes": classes_of(design)},
            "cancellation": {
                "records": len(cancel),
                "transport_classes": classes_of(cancel),
                "labels": [r["label"] for r in cancel],
            },
            "interference_generic": {
                "records": len(generic),
                "transport_classes": classes_of(generic),
            },
        },
        "P_NSR_1": {
            "claim": (
                "ratio == 1 exactly on DESIGN and CANCELLATION, and > 1 on "
                "every generic INTERFERENCE state"
            ),
            "registered_expectation": "FAIL, with roughly twelve exceptions",
            "verdict": "PASS" if not p1_exceptions else "FAIL",
            "exceptions": len(p1_exceptions),
            "exception_transport_classes": classes_of(p1_exceptions),
            "exception_labels": [r["label"] for r in p1_exceptions],
        },
        "P_NSR_2": {
            "claim": "ratio >= 1 on every record",
            "registered_expectation": "PASS",
            "verdict": "PASS" if not p2_exceptions else "FAIL",
            "exceptions": len(p2_exceptions),
            "exception_labels": [r["label"] for r in p2_exceptions],
        },
        "P_NSR_3": {
            "claim": (
                "the generic INTERFERENCE records with ratio 1 share a "
                "structural property"
            ),
            "verdict": "QUALITATIVE FINDING (post-hoc, no out-of-sample test)",
            "finding": (
                "A necessary condition is visible and a sufficient one is not. "
                "Every one of the 13 flat generic records has incidence kernel "
                "dimension at least 1 and a one-hop class of size at least 2, "
                "while 93 of the 141 that grow have kernel dimension 0 and "
                "maximum class size 1. But 48 records that grow also have "
                "kernel dimension at least 1, so kernel freedom is NECESSARY "
                "for a generic interference state to stay flat and is NOT "
                "SUFFICIENT. No law is claimed. The honest reading is that "
                "the ratio separates DESIGN and CANCELLATION from generic "
                "INTERFERENCE only in one direction: ratio > 1 implies "
                "generic interference, while ratio 1 does not imply "
                "design or cancellation."
            ),
            "flat_generic_records": len(flat_generic),
            "flat_generic_transport_classes": classes_of(flat_generic),
            "flat_generic_labels": [r["label"] for r in flat_generic],
            "comparison": {
                "flat": {
                    "one_hop_classes": profile(flat_generic, "one_hop_classes"),
                    "max_class_size": profile(flat_generic, "max_class_size"),
                    "incidence_kernel_dim": profile(flat_generic, "incidence_kernel_dim"),
                    "s_sparse": profile(flat_generic, "s_sparse"),
                },
                "grew": {
                    "one_hop_classes": profile(grew, "one_hop_classes"),
                    "max_class_size": profile(grew, "max_class_size"),
                    "incidence_kernel_dim": profile(grew, "incidence_kernel_dim"),
                },
            },
        },
        "ratio_distribution_generic": {
            "n": len(ratios),
            "min": min(ratios),
            "max": max(ratios),
            "mean": statistics.mean(ratios),
            "median": statistics.median(ratios),
            "equal_to_one": sum(1 for x in ratios if x == 1),
            "greater_than_one": sum(1 for x in ratios if x > 1),
        },
        "records": rows,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    doc = build()
    text = json.dumps(doc, indent=1) + "\n"
    if a.check:
        cur = Path(a.out).read_text() if Path(a.out).exists() else ""
        if cur != text:
            print(f"STALE: {a.out}", file=sys.stderr)
            return 1
        print(f"{a.out} is current")
        return 0

    Path(a.out).write_text(text)
    print(f"wrote {a.out}")
    for key in ("P_NSR_1", "P_NSR_2"):
        v = doc[key]
        print(f"  {key}: {v['verdict']}  exceptions {v['exceptions']}")
    d = doc["ratio_distribution_generic"]
    print(f"  generic interference ratio: min {d['min']:.2f} median "
          f"{d['median']:.2f} max {d['max']:.2f}; "
          f"{d['equal_to_one']} at exactly 1")
    p3 = doc["P_NSR_3"]
    print(f"  flat generic records: {p3['flat_generic_records']} in "
          f"{p3['flat_generic_transport_classes']} transport classes")
    print(f"    {p3['flat_generic_labels']}")
    print(f"    flat  max_class_size: {p3['comparison']['flat']['max_class_size']}")
    print(f"    grew  max_class_size: {p3['comparison']['grew']['max_class_size']}")
    print(f"    flat  kernel dim   : {p3['comparison']['flat']['incidence_kernel_dim']}")
    print(f"    grew  kernel dim   : {p3['comparison']['grew']['incidence_kernel_dim']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
