#!/usr/bin/env python3
"""Regenerate results/data/census_master.csv from the state ledger.

The file is the flat one-row-per-vertex summary of the census: system, vertex
index, spectrum, and classification verdict. It had no generator, which the
completeness audit caught, and the shipped copy carried two defects that a
generator makes impossible to reintroduce.

1. It was malformed CSV. The system field is written ``(3,6)``, whose comma was
   never quoted, so every row parsed as seven fields against a six-field
   header and any standard reader silently mis-split them. This writer uses
   ``csv``, which quotes the field.

2. It was missing the whole ``(3,8)`` block, 761 rows against 799 vertices.
   Deriving from ``states.jsonl``, which is complete by construction and
   guarded by its own tests, means the row count now cannot drift from the
   ledger.

The ``provenance`` column records how the verdict was obtained, and it is
deliberately coarse: ``direct`` for a classification carried in the ledger.
Anything finer belongs in the artifact that proves it, not in a summary table.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
STATES = DATA / "states.jsonl"
DEFAULT_OUT = DATA / "census_master.csv"

FIELDS = ("system", "index", "denominator", "integer_form", "verdict", "provenance")


def system_sort_key(system: str) -> tuple[int, int]:
    n, d = (int(part) for part in system.strip("()").split(","))
    return n, d


def rows() -> list[dict[str, str]]:
    records = [
        json.loads(line) for line in STATES.read_text().splitlines() if line.strip()
    ]
    records.sort(key=lambda r: (system_sort_key(r["system"]), r["index"]))
    out = []
    for record in records:
        denominator = record.get("denominator")
        out.append({
            "system": record["system"],
            "index": str(record["index"]),
            # the historical file left denominator 1 blank; keep that, it means
            # the integer form is already the spectrum
            "denominator": "" if denominator in (None, 1) else str(denominator),
            "integer_form": "[" + ",".join(str(x) for x in record["integer_form"]) + "]",
            "verdict": record["classified"],
            "provenance": "direct",
        })
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit 1 if the file on disk differs from a fresh build",
    )
    args = parser.parse_args()

    table = rows()
    lines = [",".join(FIELDS)]
    import io

    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, lineterminator="\n")
    writer.writerows(table)
    text = "\n".join(lines) + "\n" + buffer.getvalue()

    if args.check:
        current = args.out.read_text() if args.out.exists() else ""
        if current != text:
            print(f"{args.out} differs from a fresh build")
            return 1
        print(f"{args.out} is current: {len(table)} rows")
        return 0

    args.out.write_text(text)
    systems = sorted({row["system"] for row in table}, key=system_sort_key)
    print(f"wrote {len(table)} rows to {args.out} covering {len(systems)} systems")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
