#!/usr/bin/env python3
"""Generate a complete fixed-``N=3`` rank by exact affine-flat closure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from gpc_census.generation.reference_generator import (
    generate_fixed_n3_reference_rank,
)


ROOT = Path(__file__).resolve().parents[1]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rank", type=int, default=7)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--ressayre-attempts", type=int, default=32)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results/data/stage1_ressayre_3_7.json",
    )
    args = parser.parse_args(argv)
    result = generate_fixed_n3_reference_rank(
        args.rank,
        ressayre_attempts=args.ressayre_attempts,
        workers=args.workers,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result.as_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"rank {args.rank}: {len(result.enumeration.hyperplanes)} hyperplanes, "
        f"{result.screening.counts.certified_rows} certified rows, "
        f"{result.candidate_system.counts.retained_rows} facets"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
