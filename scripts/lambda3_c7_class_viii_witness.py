#!/usr/bin/env python3
"""Emit the exact non-toric witness that closes class VIII of Lambda^3 C^7.

The target is the sole outer vertex previously missing from the certified
inner hull.  The witness is a signed six-determinant state with squared
weights (2,1,1,1,2,2)/9.  Exact one-body cancellation gives the target
spectrum, and exact tangent rank 34 puts the state in class VIII.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census import orbit  # noqa: E402


def build() -> dict:
    certificate = orbit.class_viii_interference_certificate()
    if certificate["one_hop_free"]:
        raise AssertionError("the class-VIII witness must exercise interference")
    return {
        "schema_version": 1,
        "system": "Lambda^3 C^7",
        "class": "VIII",
        "class_orbit_dimension": 34,
        "status": "EXACT",
        "claim": (
            "The spectrum (5/9,5/9,5/9,1/3,1/3,1/3,1/3) lies in the "
            "GL_7 orbit of class VIII. Adding it to the certified inner hull "
            "closes the existing inner/outer sandwich for that class."
        ),
        "discovery": {
            "method": (
                "integer incidence search at denominator 9 followed by exact "
                "signed cancellation of every one-body coherence channel"
            ),
            "skeletons_examined_before_hit": 590,
            "search_exhaustiveness_needed_for_claim": False,
        },
        "certificate": certificate,
        "independent_verifier": (
            "scripts/verify_lambda3_c7_class_viii_witness_standalone.py"
        ),
        "scope": (
            "This proves one exact orbit point and closes class VIII through "
            "the already-certified sandwich. It does not close classes VI or VII."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--out",
        default="results/data/lambda3_c7_class_viii_witness.json",
    )
    args = parser.parse_args()
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(build(), indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
