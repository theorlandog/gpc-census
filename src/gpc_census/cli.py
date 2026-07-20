"""Command-line interface for gpc-census.

Two modes, mirroring the project's dual nature as a dataset and an engine:

- Serve precomputed results (fast, no solve): ``export`` and ``states`` read the
  shipped census (constraints, vertices, classifications, closed-form states) in
  machine-readable JSON, so the results are usable as a library / data source.
- Recompute or extend with the engine: ``constraints``, ``polytope``,
  ``classify``, ``solve``, and ``states --recompute`` run the algorithms
  directly. The compute knobs (--max-card, --max-blocks, --max-clique,
  --max-cliques) bound the ansatz search, so more compute reaches further.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence

from gpc_census import __version__
from gpc_census.core import slater_vertices


def _add_nd(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("-d", "--orbitals", type=int, required=True,
                    help="number of orbitals d")
    sp.add_argument("-n", "--fermions", type=int, required=True,
                    help="number of fermions n")


def _add_json(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("--json", action="store_true", help="machine-readable JSON output")


def _add_knobs(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("--max-card", type=int, default=16,
                    help="max support cardinality the weights-first solver enumerates")
    sp.add_argument("--max-blocks", type=int, default=2,
                    help="max 2x2 natural-orbital blocks in an ansatz")
    sp.add_argument("--max-clique", type=int, default=3,
                    help="max clique (block) size k; >=3 enables the k>=3 solver")
    sp.add_argument("--max-cliques", type=int, default=1,
                    help="max number of disjoint cliques; 0 uses per-vertex capacity")
    sp.add_argument("--clique-timeout", type=float, default=60.0,
                    help="wall-clock budget in seconds for the k>=3 clique sweep "
                         "per vertex (raise it to let slow-but-solvable vertices finish)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gpc-census",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Construct exact extremal states for fermionic "
            "natural-occupation-number (moment) polytopes."
        ),
        epilog=(
            "exporting states (the precomputed closed forms):\n"
            "  gpc-census states -n 4 -d 9                 # every certified state, JSON\n"
            "  gpc-census states -n 4 -d 9 --index 65      # one vertex (v_B)\n"
            "  gpc-census export -n 4 -d 9 --kind states   # states only, export format\n"
            "  gpc-census export -n 4 -d 9                 # everything for (n, d)\n"
            "\n"
            "  add --source solve to recompute with the engine instead of the\n"
            "  shipped lookup, or --source hybrid to lookup-then-solve. -h on any\n"
            "  subcommand for its own options."
        ),
    )
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    v = sub.add_parser(
        "vertices",
        help="list Slater-determinant vertices of the Pauli polytope Delta(d, n)")
    _add_nd(v)

    con = sub.add_parser(
        "constraints", help="print the (n, d) moment-polytope constraint system")
    _add_nd(con)

    poly = sub.add_parser(
        "polytope", help="enumerate all vertices of the (n, d) moment polytope (lrs)")
    _add_nd(poly)
    _add_json(poly)

    cls = sub.add_parser(
        "classify", help="classify every vertex of (n, d) as design or interference")
    _add_nd(cls)
    _add_json(cls)

    sol = sub.add_parser(
        "solve", help="construct + certify a closed-form state for a spectrum")
    _add_nd(sol)
    sol.add_argument("--spectrum", type=str, required=True,
                     help="integer occupations at the natural denominator, comma separated")
    sol.add_argument("--den", type=int, required=True, help="natural denominator")
    _add_knobs(sol)
    _add_json(sol)

    exp = sub.add_parser(
        "export",
        help="dump precomputed results for (n, d) as JSON (constraints, "
             "vertices, classification, states)")
    _add_nd(exp)
    exp.add_argument("--kind", choices=("all", "constraints", "vertices",
                                        "classification", "states"),
                     default="all", help="which precomputed artifact to emit")

    st = sub.add_parser(
        "states",
        help="closed-form states for (n, d), with explicit provenance per record")
    _add_nd(st)
    st.add_argument("--index", type=int, default=None,
                    help="restrict to a single vertex index")
    st.add_argument("--source", choices=("precompute", "solve", "hybrid"),
                    default="precompute",
                    help="precompute: shipped lookup only (default); solve: run "
                         "the engine for every vertex; hybrid: lookup where "
                         "available, solve otherwise. Every record is tagged "
                         "with its source (precomputed / solved).")
    _add_knobs(st)
    _add_json(st)
    return parser


def _cmd_solve(args, parser) -> int:
    from fractions import Fraction
    from gpc_census.states import certify_state
    spec = [Fraction(int(x), args.den) for x in args.spectrum.split(",")]
    rec = certify_state(args.fermions, args.orbitals, spec,
                        max_card=args.max_card, max_blocks=args.max_blocks,
                        max_clique=args.max_clique, max_cliques=args.max_cliques,
                        clique_time_budget=args.clique_timeout)
    if args.json:
        print(json.dumps(rec, default=str, indent=1))
        return 0
    if rec is None or rec.get("status") != "OK":
        print("no state found", rec.get("reason", "") if rec else "")
        return 1
    ex = rec.get("exact") or {}
    if ex.get("status") == "EXACT":
        print(f"closed form (den {ex['den']}), weights {ex['weights']}")
        for det, amp in zip(ex.get("support_dets", [s[0] for s in rec["support"]]),
                            ex["pretty"]):
            print("  ", det, amp)
    else:
        print(f"numeric state, residual {rec.get('residual')}, "
              f"no closed form ({ex.get('status', 'n/a')})")
        for det, mod, ph in rec["support"]:
            print("  ", det, f"{mod:.9f}", f"{ph:+.9f}")
    return 0


def _cmd_states(args, parser) -> int:
    from gpc_census.dataset import resolve_states
    recs = resolve_states(args.fermions, args.orbitals, index=args.index,
                          mode=args.source, max_card=args.max_card,
                          max_blocks=args.max_blocks, max_clique=args.max_clique,
                          max_cliques=args.max_cliques,
                          clique_time_budget=args.clique_timeout)
    # provenance is the mode, not the record: precompute serves the lookup,
    # solve recomputes independently, hybrid serves the lookup and solves the rest
    envelope = {"mode": args.source, "states": recs}
    print(json.dumps(envelope, default=str, indent=1))
    return 0


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
    if args.command == "constraints":
        from gpc_census.constraints import constraints
        print(json.dumps(constraints(args.fermions, args.orbitals), indent=1))
        return 0
    if args.command == "polytope":
        from gpc_census.polytope import vertices as pverts
        vs = [[str(x) for x in v] for v in pverts(args.fermions, args.orbitals)]
        if args.json:
            print(json.dumps(vs, indent=1))
        else:
            for v in vs:
                print(" ".join(v))
        return 0
    if args.command == "classify":
        from gpc_census.classify import classify_full
        from gpc_census.polytope import vertices as pverts
        rows = []
        for i, v in enumerate(pverts(args.fermions, args.orbitals)):
            r = classify_full(args.fermions, args.orbitals, v)
            rows.append({"index": i, "spectrum": [str(x) for x in v],
                         "verdict": r["verdict"], "backend": r.get("backend")})
        if args.json:
            print(json.dumps(rows, indent=1))
        else:
            for row in rows:
                print(row["index"], row["spectrum"], row["verdict"], row["backend"])
        return 0
    if args.command == "solve":
        return _cmd_solve(args, parser)
    if args.command == "export":
        from gpc_census import dataset
        n, d = args.fermions, args.orbitals
        try:
            if args.kind == "all":
                payload = dataset.export(n, d)
            elif args.kind == "constraints":
                from gpc_census.constraints import constraints
                payload = constraints(n, d)
            elif args.kind == "vertices":
                payload = dataset.vertices(n, d)
            elif args.kind == "classification":
                payload = dataset.classification(n, d)
            else:
                payload = dataset.states(n, d)
        except KeyError as exc:
            parser.error(str(exc))
        print(json.dumps(payload, default=str, indent=1))
        return 0
    if args.command == "states":
        return _cmd_states(args, parser)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
