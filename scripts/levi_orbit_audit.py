"""Levi-orbit census audit.

Runs the Levi-symmetry theorem over every certified state in
`results/data/states.jsonl` and enforces it as a hard gate.

Two layers:

STATIC (basis-independent, no state vector).  For each record compute
`static_levi_record`.  The theorem says: [psi] carries a one-dimensional
L_lambda-subrepresentation iff lambda is {0,1}-valued (a Slater spectrum).
`levi_theorem_consistent` is the falsifiable form of that equivalence -- it
compares whether a one-dimensional branching sector's INDUCED occupation
pattern reproduces the target spectrum against whether the spectrum is Slater.
Any record where these disagree fails the gate with SystemExit.

DYNAMIC (from the certified closed form).  Reconstruct psi, recompute the
1-RDM, and measure the Levi orbit dimension.  The theorem predicts
levi_orbit_dim == 0 exactly for the Slater states and > 0 for every correlated
state.  We also cross-check that the recorded spectrum's multiplicity pattern
matches the multiplicities of the recomputed 1-RDM (a data-integrity guard).

Run from the repo root with the package installed (make sync).  Needs numpy and
sympy.  ~15 s for the full census.
"""
from __future__ import annotations
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from gpc_census.levi import (
    fiber_and_orbit,
    spectrum_multiplicities,
    state_vector,
    static_levi_record,
)

ROOT = Path(__file__).resolve().parents[1]
STATES = ROOT / "results" / "data" / "states.jsonl"


def load_records(path=STATES):
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def main():
    records = load_records()
    n_records = len(records)

    theorem_failures = []
    multiplicity_mismatches = []
    stab_below_two = []
    slater_nonzero_orbit = []
    correlated_zero_orbit = []

    hard_by_class = Counter()
    total_by_class = Counter()
    n_slater = 0
    hard_records = []

    for rec in records:
        static = static_levi_record(rec)
        cls = rec.get("classified", "?")
        total_by_class[cls] += 1

        if not static["levi_theorem_consistent"]:
            theorem_failures.append((rec["system"], rec["index"]))
        if static["hard_obstruction"]:
            hard_by_class[cls] += 1
            hard_records.append((rec["system"], rec["index"], static["integer_form"],
                                 static["multiplicities"], cls))
        if static["slater_spectrum"]:
            n_slater += 1

        # dynamic
        psi, basis, n, d = state_vector(rec)
        info = fiber_and_orbit(psi, basis, n, d)

        recorded_mults = sorted(spectrum_multiplicities(static["integer_form"])[0])
        numeric_mults = sorted(info["numerical_multiplicities"])
        if recorded_mults != numeric_mults:
            multiplicity_mismatches.append((rec["system"], rec["index"],
                                            recorded_mults, numeric_mults))

        if info["projective_stabilizer_dim"] < 2:
            stab_below_two.append((rec["system"], rec["index"],
                                   info["projective_stabilizer_dim"]))

        if static["slater_spectrum"] and info["levi_orbit_dim"] != 0:
            slater_nonzero_orbit.append((rec["system"], rec["index"],
                                         info["levi_orbit_dim"]))
        if not static["slater_spectrum"] and info["levi_orbit_dim"] <= 0:
            correlated_zero_orbit.append((rec["system"], rec["index"]))

    n_hard = sum(hard_by_class.values())

    print(f"records                : {n_records}")
    print(f"slater spectra         : {n_slater}")
    print(f"non-slater spectra     : {n_records - n_slater}")
    print(f"hard obstruction       : {n_hard} / {n_records}")
    for cls in sorted(total_by_class):
        tot = total_by_class[cls]
        hard = hard_by_class[cls]
        print(f"  {cls:<13} {tot:>4}   {hard} hard ({100.0 * hard / tot:.1f}%)")
    if hard_records:
        print("hard records:")
        for system, index, form, mults, cls in hard_records:
            print(f"  {system} v{index} {form} blocks {tuple(mults)} [{cls}]")

    print("--- theorem gate ---")
    print(f"levi_theorem_consistent failures : {len(theorem_failures)}")
    print(f"multiplicity mismatches          : {len(multiplicity_mismatches)}")
    print(f"projective_stabilizer_dim < 2    : {len(stab_below_two)}")
    print(f"slater states with orbit != 0    : {len(slater_nonzero_orbit)}")
    print(f"correlated states with orbit <=0 : {len(correlated_zero_orbit)}")

    failed = False
    if theorem_failures:
        failed = True
        print(f"THEOREM VIOLATED on {theorem_failures[:5]} ...", file=sys.stderr)
    if slater_nonzero_orbit:
        failed = True
        print(f"SLATER ORBIT NONZERO on {slater_nonzero_orbit[:5]} ...", file=sys.stderr)
    if correlated_zero_orbit:
        failed = True
        print(f"CORRELATED ORBIT ZERO on {correlated_zero_orbit[:5]} ...", file=sys.stderr)
    if multiplicity_mismatches:
        failed = True
        print(f"MULTIPLICITY MISMATCH on {multiplicity_mismatches[:5]} ...", file=sys.stderr)

    if failed:
        raise SystemExit("Levi audit FAILED")
    print("Levi audit PASS")


if __name__ == "__main__":
    main()
