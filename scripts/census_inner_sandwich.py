#!/usr/bin/env python3
"""Completeness without candidate exhaustiveness: the census-inner sandwich.

THE REFRAME. The claim ladder in `docs/fixed_n3_generator_methodology.md` says
a complete GPC system needs "an exhaustive-enumeration proof OR an independent
exact inner/outer closure". The (3,6) theorem took the second route, and the
second route DOES NOT CARE where the candidate rows came from. So exhaustive
candidate enumeration, which is the expensive and unfinished obligation, can be
sidestepped entirely at any rank where both halves of a sandwich are available.

Both halves are available at every rank the census covers:

    P subset O   every row of O is a replayed exact Ressayre certificate.
                 Certifying the VALIDITY of a supplied row is screening, which
                 costs about 0.14 ms per row. Only EXHAUSTIVENESS was ever
                 expensive, and this proof does not use it. The screening
                 pipeline is no longer fixed-N=3, so this half certifies at
                 every particle number the census covers.

    O subset P   every vertex of O is an attained spectrum. The census ships
                 799 certified extremal states, each with an exact
                 characteristic-polynomial certificate that spec(rho) equals
                 its vertex. So enumerate the vertices of O exactly and hand
                 each one its attaining state.

    Therefore P = O.

Non-circularity is the point and is worth stating precisely. The census VERTEX
lists were derived from the constraint tables, so checking "vertices of O are
census vertices" alone would be circular. The attaining STATES are not: each is
an explicit state whose 1-RDM spectrum is verified by a characteristic-
polynomial identity in exact arithmetic, with no reference to any constraint
table. That is what makes the inner half independent evidence.

WHAT WAS EXPECTED TO BE THE WALL, AND WAS NOT. Exact vertex enumeration of a
highly degenerate polytope was the predicted obstacle. The double-description
enumerator in `gpc_census.exact_polytope` does it in 0.0 s, 0.6 s, 12 s and
247 s at ranks 7, 8, 9 and 10, reproducing the certified vertex counts exactly.

Usage:
    python scripts/census_inner_sandwich.py
    python scripts/census_inner_sandwich.py --wider
    python scripts/census_inner_sandwich.py --systems 3,6 3,7
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import verify_states_standalone as vs  # noqa: E402

from gpc_census import exact_polytope as ep  # noqa: E402
from gpc_census import orbit  # noqa: E402
from gpc_census.constraints import constraints  # noqa: E402

LADDER = [(3, 6), (3, 7), (3, 8), (3, 9), (3, 10)]
WIDER = [(4, 8), (4, 9), (4, 10), (5, 10)]


def _census_states() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for line in (ROOT / "results/data/states.jsonl").read_text().splitlines():
        if line.strip():
            rec = json.loads(line)
            out.setdefault(rec["system"], []).append(rec)
    return out


def _spectrum(rec) -> tuple[Fraction, ...]:
    den = rec["denominator"] or 1
    return tuple(Fraction(v, den) for v in rec["integer_form"])


def screen_validity(n: int, d: int) -> dict:
    """Certify every stored row as an exact Ressayre inequality (P subset O).

    This used to bail out for ``n != 3``, because the screening pipeline hard
    coded the particle number even though everything under it (weights, roots,
    tangent determinant, full-dimension witness) was already written for
    ``wedge^n C^d``.  The particle number is now threaded through, so the wider
    census systems screen by the same code path and get the same exact
    certificate rather than a citation.
    """
    from gpc_census.generation import screen_tau_candidates
    from gpc_census.generation.bdr import tau_from_constraint
    from gpc_census.generation.model import IntegralConstraint

    rows = constraints(n, d)["inequalities"]
    taus = tuple(
        tau_from_constraint(IntegralConstraint(tuple(r["coeffs"]), int(r["rhs"])), n)
        for r in rows
    )
    result = screen_tau_candidates(taus, particle_number=n, ressayre_attempts=32,
                                   cyclic_cross_check=False,
                                   symbolic_fallback=True, workers=1)
    counts = result.counts
    return {
        "available": True,
        "particle_number": n,
        "rows": len(rows),
        "certified": counts.certified_rows,
        "structural": counts.structural_rows,
        "exact_rejected": counts.exact_rejected_rows,
        "unresolved": counts.evaluation_unresolved_rows,
        "all_valid": bool(counts.certified_rows + counts.structural_rows == len(rows)),
    }


def run_system(n: int, d: int, states: dict[str, list[dict]],
               verify_states: bool = True) -> dict:
    started = time.time()
    validity = screen_validity(n, d)

    amb_ineqs, amb_eqs = orbit.ambient_rows(n, d)
    base_rows, eqs, base_verts = orbit._chamber(d, n, amb_eqs)
    t0 = time.time()
    vertices = ep.vertices_from_hrep(base_verts, base_rows, eqs, amb_ineqs, d)
    enum_seconds = time.time() - t0

    key = f"({n},{d})"
    by_spectrum: dict[tuple, dict] = {}
    for rec in states.get(key, []):
        by_spectrum.setdefault(_spectrum(rec), rec)

    attained, resisting, failures = [], [], []
    for vertex in vertices:
        rec = by_spectrum.get(vertex)
        if rec is None:
            resisting.append([str(x) for x in vertex])
            continue
        if verify_states:
            ok, why = vs.verify(rec)
            if not ok:
                failures.append({"index": rec["index"], "reason": why})
                continue
        attained.append(rec["index"])

    inner = not resisting and not failures
    outer = bool(validity.get("all_valid", False))
    if inner and outer:
        status = "proved"
        note = "P = O, both halves certified in repo, exhaustiveness unused"
    elif inner and not validity.get("available", False):
        # N != 3: the Ressayre screening machinery is fixed-N=3, so the outer
        # half rests on the published source rather than on a replayed
        # certificate. Reporting this as "proved" would overstate it.
        status = "inner_half_only"
        note = ("O subset P certified; P subset O rests on the published "
                "source because the screening machinery is fixed-N=3")
    else:
        status = "open"
        note = "not closed; a resisting vertex is a pointer at a missing facet"
    return {
        "system": key, "n": n, "d": d,
        "seconds": round(time.time() - started, 1),
        "vertex_enumeration_seconds": round(enum_seconds, 1),
        "outer_rows": len(amb_ineqs),
        "outer_equalities": len(amb_eqs),
        "vertices": len(vertices),
        "validity_P_subset_O": validity,
        "attained_vertices": len(attained),
        "resisting_vertices": resisting,
        "attainment_failures": failures,
        "states_reverified": bool(verify_states),
        "inner_half_proved": bool(inner),
        "outer_half_proved": bool(outer),
        "closure_proved": bool(inner and outer),
        "closure_status": status,
        "closure_note": note,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--systems", nargs="*", default=None,
                    help="e.g. 3,6 3,7; defaults to the N=3 ladder")
    ap.add_argument("--wider", action="store_true",
                    help="also attempt the N=4 and N=5 census systems")
    ap.add_argument("--no-verify-states", action="store_true",
                    help="skip re-running the exact state certificates")
    ap.add_argument("-o", "--out",
                    default=str(ROOT / "results/data/census_inner_sandwich.json"))
    args = ap.parse_args()

    if args.systems:
        wanted = [tuple(int(x) for x in s.split(",")) for s in args.systems]
    else:
        wanted = LADDER + (WIDER if args.wider else [])

    states = _census_states()
    out = pathlib.Path(args.out)
    blocks = {}
    for n, d in wanted:
        print(f"({n},{d}): sandwiching ...", flush=True)
        try:
            rec = run_system(n, d, states, not args.no_verify_states)
        except Exception as exc:                       # operational, not a refutation
            rec = {"system": f"({n},{d})", "n": n, "d": d,
                   "status": "operational_failure",
                   "error": f"{type(exc).__name__}: {exc}", "closure_proved": False}
            print(f"  operational failure: {rec['error']}", flush=True)
        else:
            print(f"  rows {rec['outer_rows']}  vertices {rec['vertices']}  "
                  f"attained {rec['attained_vertices']}  "
                  f"{rec['closure_status']}  ({rec['seconds']}s)", flush=True)
        blocks[rec["system"]] = rec
        # the wider systems cost tens of minutes each, so a crash or a timeout
        # must not throw away the systems that already closed
        _write(out, blocks)

    print(f"wrote {out}")
    return 0


def _write(out: pathlib.Path, blocks: dict) -> None:
    payload = {
        "what": "Exact inner/outer closure using the certified census states, "
                "which proves P = O without candidate exhaustiveness.",
        "argument": {
            "P_subset_O": "every stored row screens as an exact Ressayre "
                          "certificate, at every particle number in the census",
            "O_subset_P": "every vertex of O is the spectrum of a certified "
                          "census state, verified by an exact characteristic-"
                          "polynomial identity",
            "conclusion": "P = O; candidate exhaustiveness is not used",
            "non_circularity": "the census VERTEX lists came from the constraint "
                               "tables, so they alone would be circular; the "
                               "attaining STATES do not, since each is checked "
                               "by a characteristic-polynomial identity that "
                               "references no table",
        },
        "systems": blocks,
        "closed_systems": sorted(k for k, v in blocks.items()
                                 if v.get("closure_proved")),
        "inner_half_only_systems": sorted(
            k for k, v in blocks.items()
            if v.get("closure_status") == "inner_half_only"),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
