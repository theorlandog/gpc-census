#!/usr/bin/env python3
"""Run the attainer gate over every library record and report every violation.

The gate is `gpc_census.validate.check_vertex_attainer`: a state offered as
attaining a vertex must have a diagonal 1-RDM whose diagonal is lambda
entrywise, and unit norm. Conditions 1 and 2 are what make a state a
NATURAL-ORBITAL representative, so that its support bounds s_Q^NO rather than
the free-basis minimum s_Q^free (RESEARCH.md, "The support taxonomy").

WHAT THIS AUDIT FOUND (2026-07). 154 of the 799 shipped records fail the gate:
every INTERFERENCE record except v89 and v103. All 154 still attain the
SPECTRUM exactly and have unit norm; what they lack is diagonality. This is
not a newly introduced defect and not confined to one record:

- the library's own certificate is the characteristic-polynomial identity
  (scripts/verify_states_standalone.py), which certifies spec(rho) = lambda
  and never asserted diagonality except for design states, where it follows
  from one-hop freedom by inspection;
- the repository already records the same fact from the other side: the
  cancellation regime, defined as an exactly diagonal 1-RDM with silent
  channels, "has exactly TWO members across all 799 certified
  representatives, v89 and v103" (RESEARCH.md), and cascade.jsonl carries a
  per-state `rdm_max_offdiagonal` that is nonzero for the magnitude regime.

So the 154 are spectrum attainers in a non-natural-orbital basis BY
CONSTRUCTION. The consequence that does need fixing is a documentation one:
the support taxonomy says the library support column bounds s_Q^NO, and for
these 154 records it bounds only s_Q^free.

Usage:
  python scripts/attainer_gate_audit.py            # write the artifact
  python scripts/attainer_gate_audit.py --check    # fail if stale
  python scripts/attainer_gate_audit.py --summary  # counts only
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.fiber import amplitudes_from_record, one_rdm  # noqa: E402
from gpc_census.validate import check_vertex_attainer  # noqa: E402

DATA = ROOT / "results" / "data"
OUT = DATA / "attainer_gate_audit.json"
TOL = 1e-9


def audit_record(rec):
    amps, dets = amplitudes_from_record(rec)
    lam = [v / rec["denominator"] for v in rec["integer_form"]]
    d = len(lam)
    c = np.asarray(amps, dtype=complex)
    dets = [tuple(int(o) for o in t) for t in dets]

    rho = one_rdm(c, dets, d)
    off = float(np.max(np.abs(rho - np.diag(np.diag(rho)))))
    diag_err = float(np.max(np.abs(np.real(np.diag(rho)) - np.array(lam))))
    norm_err = abs(float(np.linalg.norm(c)) - 1.0)
    eig = np.sort(np.linalg.eigvalsh(rho))[::-1]
    spec_err = float(np.max(np.abs(eig - np.sort(np.array(lam))[::-1])))
    # is the diagonal merely permuted, i.e. a relabelling away from the gate?
    permuted = float(
        np.max(np.abs(np.sort(np.real(np.diag(rho)))[::-1] - np.sort(np.array(lam))[::-1]))
    )

    return {
        "system": rec["system"],
        "index": rec["index"],
        "label": f"{rec['system']} v{rec['index']}",
        "classified": rec["classified"],
        "support": len(dets),
        "max_offdiagonal": off,
        "diagonal_error": diag_err,
        "norm_error": norm_err,
        "spectrum_error": spec_err,
        "diagonal_is_permutation_of_lambda": bool(permuted < TOL),
        "passes_gate": bool(off < TOL and diag_err < TOL and norm_err < TOL),
        "attains_spectrum": bool(spec_err < 1e-8),
        "gate_errors": check_vertex_attainer(c, dets, lam, d=d),
        "evidence": "numerical",
    }


def build():
    records = [
        json.loads(ln)
        for ln in (DATA / "states.jsonl").read_text().splitlines()
        if ln.strip()
    ]
    rows = [audit_record(r) for r in records]
    failures = [r for r in rows if not r["passes_gate"]]
    passes = [r for r in rows if r["passes_gate"]]

    def by_class(sel):
        out = {}
        for r in sel:
            out[r["classified"]] = out.get(r["classified"], 0) + 1
        return dict(sorted(out.items()))

    interference = [r for r in rows if r["classified"] == "INTERFERENCE"]
    passing_interference = sorted(r["label"] for r in interference if r["passes_gate"])

    return {
        "generated_by": "scripts/attainer_gate_audit.py",
        "gate": "gpc_census.validate.check_vertex_attainer",
        "tolerance": TOL,
        "evidence": "numerical, cross-checked in exact arithmetic at (3,10) v17",
        "summary": {
            "records": len(rows),
            "pass": len(passes),
            "fail": len(failures),
            "pass_by_class": by_class(passes),
            "fail_by_class": by_class(failures),
            "all_records_attain_the_spectrum": all(r["attains_spectrum"] for r in rows),
            "worst_spectrum_error": max(r["spectrum_error"] for r in rows),
            "worst_norm_error": max(r["norm_error"] for r in rows),
            "failures_that_are_merely_permuted": sum(
                1 for r in failures if r["diagonal_is_permutation_of_lambda"]
            ),
            "interference_records": len(interference),
            "interference_passing_the_gate": passing_interference,
            "worst_offdiagonal_among_failures": max(
                (r["max_offdiagonal"] for r in failures), default=0.0
            ),
        },
        "reading": (
            "154 of 799 records fail, and they are exactly the INTERFERENCE "
            "records other than v89 and v103. Every one of them still attains "
            "the spectrum exactly and has unit norm; what they lack is a "
            "diagonal 1-RDM. None is merely a permuted relabelling. This is a "
            "property of how the interference library was built, not a "
            "regression: the shipped certificate is the "
            "characteristic-polynomial identity, which certifies "
            "spec(rho) = lambda, and the repository separately records that "
            "exactly two certified representatives have an exactly diagonal "
            "1-RDM. The consequence needing repair is a documentation one: the "
            "support column of these 154 records bounds s_Q^free, not s_Q^NO."
        ),
        "records": rows,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()

    doc = build()
    if a.summary:
        print(json.dumps(doc["summary"], indent=2))
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
    s = doc["summary"]
    print(f"wrote {a.out}")
    print(f"  records {s['records']}  pass {s['pass']}  fail {s['fail']}")
    print(f"  pass by class: {s['pass_by_class']}")
    print(f"  fail by class: {s['fail_by_class']}")
    print(f"  interference passing the gate: {s['interference_passing_the_gate']}")
    print(f"  all attain the spectrum: {s['all_records_attain_the_spectrum']} "
          f"(worst {s['worst_spectrum_error']:.2e})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
