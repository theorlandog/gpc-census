#!/usr/bin/env python3
"""Regenerate the complete SHA-256 manifest for results/data."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "data",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results/data",
    )
    args = parser.parse_args()
    manifest = args.data / "SHA256SUMS"
    paths = sorted(
        path for path in args.data.rglob("*") if path.is_file() and path != manifest
    )
    lines = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  ./{path.relative_to(args.data)}"
        for path in paths
    ]
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(lines)} checksums to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
