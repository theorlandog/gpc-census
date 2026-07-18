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
    for name, helptext in [
        ("constraints", "print the constraint system for the (n, d) moment polytope"),
        ("polytope", "enumerate all vertices of the (n, d) moment polytope exactly (lrs)"),
        ("classify", "classify every vertex of (n, d) as design or interference"),
        ("solve", "construct a state attaining a vertex spectrum: --spectrum a,b,c,... over --den"),
    ]:
        sp = subparsers.add_parser(name, help=helptext)
        sp.add_argument("-d", "--orbitals", type=int, required=True)
        sp.add_argument("-n", "--fermions", type=int, required=True)
        if name == "solve":
            sp.add_argument("--spectrum", type=str, required=True,
                            help="integer occupations at the natural denominator, comma separated")
            sp.add_argument("--den", type=int, required=True, help="natural denominator")
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
    if args.command == "constraints":
        import json as _json
        from gpc_census.constraints import constraints
        print(_json.dumps(constraints(args.fermions, args.orbitals), indent=1))
    if args.command == "polytope":
        from gpc_census.polytope import vertices as pverts
        for v in pverts(args.fermions, args.orbitals):
            print(" ".join(str(x) for x in v))
    if args.command == "classify":
        from gpc_census.classify import classify
        from gpc_census.polytope import vertices as pverts
        for i, v in enumerate(pverts(args.fermions, args.orbitals)):
            print(i, [str(x) for x in v], classify(args.fermions, args.orbitals, v))
    if args.command == "solve":
        from fractions import Fraction
        from gpc_census.states import attain
        spec = [Fraction(int(x), args.den) for x in args.spectrum.split(",")]
        psi, res, dets = attain(args.fermions, args.orbitals, spec)
        print(f"residual {res}")
        for i, amp in enumerate(psi):
            if abs(amp) > 1e-9:
                print(dets[i], f"{abs(amp):.9f}", f"{__import__('numpy').angle(amp):+.9f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
