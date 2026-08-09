#!/usr/bin/env python3
"""Regenerate the complete SHA-256 manifest for results/data.

``--check`` verifies instead of writing, exiting 1 when the manifest is stale
and naming exactly which entries are wrong. That mode exists because the
manifest can go stale without anyone editing it: a merge resolves
``SHA256SUMS`` independently of the files whose digests it records, so a
force-resolved conflict leaves a recorded digest describing content that no
longer hashes to it. ``PROVENANCE.md`` is the usual casualty, since almost
every change touches it.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path


def manifest_lines(data: Path) -> list[str]:
    """The manifest the data directory should have, one line per file."""
    manifest = data / "SHA256SUMS"
    paths = sorted(
        path for path in data.rglob("*") if path.is_file() and path != manifest
    )
    return [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  ./{path.relative_to(data)}"
        for path in paths
    ]


def _entries(lines: list[str]) -> dict[str, str]:
    entries = {}
    for line in lines:
        if not line.strip():
            continue
        digest, name = line.split(None, 1)
        entries[name.strip().lstrip("*").removeprefix("./")] = digest
    return entries


def check(data: Path) -> int:
    manifest = data / "SHA256SUMS"
    if not manifest.is_file():
        print(f"{manifest} is missing; run `make checksums`", file=sys.stderr)
        return 1

    expected = _entries(manifest_lines(data))
    recorded = _entries(manifest.read_text(encoding="utf-8").splitlines())
    if expected == recorded:
        print(f"{len(expected)} checksums match {data}")
        return 0

    unlisted = sorted(set(expected) - set(recorded))
    absent = sorted(set(recorded) - set(expected))
    mismatched = sorted(
        name for name in set(expected) & set(recorded)
        if expected[name] != recorded[name]
    )
    print("results/data/SHA256SUMS is stale. Run `make checksums`.", file=sys.stderr)
    for label, names in (
        ("on disk but unlisted", unlisted),
        ("listed but absent", absent),
        ("digest mismatch", mismatched),
    ):
        if names:
            print(f"  {label}: {', '.join(names)}", file=sys.stderr)
    if mismatched == ["PROVENANCE.md"]:
        print(
            "  PROVENANCE.md alone is the signature of a merge that resolved "
            "the manifest textually while the file itself merged separately.",
            file=sys.stderr,
        )
    return 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "data",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results/data",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify the manifest and exit 1 if stale, without rewriting it",
    )
    args = parser.parse_args()

    if args.check:
        return check(args.data)

    manifest = args.data / "SHA256SUMS"
    lines = manifest_lines(args.data)
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(lines)} checksums to {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
