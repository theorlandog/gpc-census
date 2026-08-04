#!/usr/bin/env python3
"""Regenerate and verify the compact ``(3,7)`` reference manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from gpc_census.generation.reference_generator import verify_reference_manifest


ROOT = Path(__file__).resolve().parents[1]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "artifact",
        nargs="?",
        type=Path,
        default=ROOT / "results/data/stage1_ressayre_3_7.json",
    )
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args(argv)
    payload = json.loads(args.artifact.read_text(encoding="utf-8"))
    if not verify_reference_manifest(payload, workers=args.workers):
        raise SystemExit("rank-7 reference manifest verification failed")
    print("rank-7 reference manifest verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
