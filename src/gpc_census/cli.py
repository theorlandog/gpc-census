"""Command-line interface for gpc-census."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from gpc_census import __version__
from gpc_census.core import slater_vertices


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gpc-census",
        description=(
            "Construct exact extremal states for fermionic "
            "natural-occupation-number (moment) polytopes."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    vertices = subparsers.add_parser(
        "vertices",
        help="list Slater-determinant vertices of the Pauli polytope Delta(d, n)",
    )
    vertices.add_argument(
        "-d", "--orbitals", type=int, required=True, help="number of orbitals d"
    )
    vertices.add_argument(
        "-n", "--fermions", type=int, required=True, help="number of fermions n"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "vertices":
        try:
            for vec in slater_vertices(args.orbitals, args.fermions):
                print(" ".join(map(str, vec)))
        except ValueError as exc:
            parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
