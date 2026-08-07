"""Exact symplectic slices at the census's landmark states.

The previous campaign reported no fiber tangent dimension because the raw
`ker d(mu)_x` is dominated by trivial ambient directions: at the `(3,6)` Slater
vertex, whose fixed-rho fiber is a single point, it has dimension 20.

The symplectic slice is the object that fixes this. With
`V_x = ker d(mu)_x / T_x(K_mu . x)`, a slice dimension of 0 certifies outright
that the reduced germ is a point, with no invariant theory needed.

Pilot scope, following the campaign brief: the Slater vertex, the uniform
Borland-Dennis vertex, v_B, v89 and v103. Results are checkpointed per state,
so a state that exceeds the compute budget leaves the others intact and is
reported as unfinished rather than silently dropped.

Run: uv run scripts/symplectic_slice_atlas.py
Writes results/data/symplectic_slices.json.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data"
OUT = DATA / "symplectic_slices.json"
CHECKPOINT = ROOT / "build" / "symplectic_slices_progress.jsonl"

from gpc_census.symplectic import symplectic_slice  # noqa: E402

PILOT = [
    ("(3,6)", 0, "Slater vertex; fixed-rho fiber is a single point"),
    ("(3,6)", 3, "uniform Borland-Dennis vertex (1/2)^6"),
    ("(4,9)", 65, "v_B, (20,14,14,14,14,4,4,4,4)/23"),
    ("(3,10)", 89, "v89, rank-10 cancellation-regime closure"),
    ("(3,10)", 103, "v103, rank-10 cancellation-regime closure"),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip", nargs="*", default=[],
                    help="states to leave unattempted, as SYSTEM#INDEX; they are "
                         "reported in the artifact's unfinished list, never dropped")
    args = ap.parse_args()
    skip = set(args.skip)
    states = {(r["system"], r["index"]): r
              for r in (json.loads(x)
                        for x in (DATA / "states.jsonl").read_text().splitlines())}
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if CHECKPOINT.exists():
        for line in CHECKPOINT.read_text().splitlines():
            row = json.loads(line)
            done[(row["system"], row["index"])] = row

    rows = []
    with CHECKPOINT.open("a") as ck:
        for system, index, label in PILOT:
            if f"{system}#{index}" in skip:
                continue
            if (system, index) in done:
                rows.append(done[(system, index)])
                continue
            record = states[(system, index)]
            t0 = time.time()
            slice_data = symplectic_slice(record)
            row = {
                "system": system, "index": index, "label": label,
                "integer_form": list(record["integer_form"]),
                "denominator": record["denominator"],
                "classified": record["classified"],
                "support_size": len(record["closed_form"]["support_dets"]),
                "evidence": "exact",
                "slice": slice_data,
                "secs": round(time.time() - t0, 1),
            }
            rows.append(row)
            ck.write(json.dumps(row) + "\n")
            ck.flush()
            print(f"  {system}#{index} slice_dim={slice_data['slice_dim']} "
                  f"ker={slice_data['ker_d_mu_dim']} gates="
                  f"{slice_data['gate_kernel_dim'] and slice_data['gate_levi_dim']} "
                  f"({row['secs']}s)", flush=True)

    attempted = {(s, i) for s, i, _ in PILOT}
    finished = {(r["system"], r["index"]) for r in rows}
    payload = {
        "generated_by": "scripts/symplectic_slice_atlas.py",
        "prereg": "docs/prereg_symplectic_slice.md",
        "evidence": "exact",
        "note": (
            "Symplectic slices at the census landmark states. A slice dimension of "
            "0 certifies that the reduced germ is a point. No reduced-germ dimension "
            "is reported where the slice is nonzero: the generic-orbit shortcut is "
            "wrong and returns -10 at the Slater vertex. See docs/symplectic_slice.md."
        ),
        "rows": rows,
        "unfinished": sorted(f"{s}#{i}" for s, i in attempted - finished),
        "gates": {
            "rank_equals_orbit_violations": [
                f"{r['system']}#{r['index']}" for r in rows
                if not r["slice"]["gate_rank_equals_orbit"]],
            "kernel_dim_violations": [
                f"{r['system']}#{r['index']}" for r in rows
                if not r["slice"]["gate_kernel_dim"]],
            "levi_dim_violations": [
                f"{r['system']}#{r['index']}" for r in rows
                if not r["slice"]["gate_levi_dim"]],
        },
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)} ({len(rows)} of {len(PILOT)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
